/*                    R T _ R E N D E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libbrlobol/rt_render.cpp
 *
 * A retained librt image renderer.  This deliberately does not use Obol or
 * Coin objects while rendering: those are sampled once by synchronize().
 */

#include "common.h"

#include "brlobol/database_source.h"
#include "brlobol/pick_detail.h"
#include "brlobol/rt_render.h"
#include "brlobol/view_controller.h"

#include "raytrace.h"

#include <Inventor/SbLine.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbPlane.h>
#include <Inventor/SbViewVolume.h>
#include <Inventor/nodes/SoCamera.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <thread>
#include <utility>

struct RtSourceKey {
    struct db_i *database = NULL;
    SbString instanceKey;
    SbString path;
    uint32_t sourceRevision = 0;

    bool operator==(const RtSourceKey &other) const
    {
	return database == other.database &&
	    instanceKey == other.instanceKey &&
	    sourceRevision == other.sourceRevision &&
	    path == other.path;
    }
};

struct RtSourceSnapshot {
    RtSourceKey key;
    SbBool visible = FALSE;
    SbBool selected = FALSE;
    SbBool highlighted = FALSE;
    SbBool colorOverride = FALSE;
    SbColor color = SbColor(0.8f, 0.8f, 0.8f);
    SbColor selectedColor = SbColor(1.0f, 1.0f, 1.0f);
    SbColor highlightedColor = SbColor(1.0f, 1.0f, 0.0f);
    SbBool materialColorValid = FALSE;
    SbColor materialColor = SbColor(0.8f, 0.8f, 0.8f);
    SbBool databaseMaterialColorValid = FALSE;
    SbColor databaseMaterialColor = SbColor(0.8f, 0.8f, 0.8f);
    float transparency = 0.0f;
    SbBool drawMatrixValid = FALSE;
    SbMatrix localToWorld = SbMatrix::identity();
    SbMatrix worldToLocal = SbMatrix::identity();
    SbPlane localClipPlanes[2];
    size_t localClipPlaneCount = 0;
};

struct RtPreparedCache {
    size_t snapshotIndex = 0;
    std::unique_ptr<BRLObolRtPickCache> cache;
};

struct RtTraceHit {
    BRLObolRtPickResult pick;
    size_t sourceIndex = 0;
};

static bool
samePresentation(const RtSourceSnapshot &a, const RtSourceSnapshot &b)
{
    return a.visible == b.visible &&
	a.selected == b.selected && a.highlighted == b.highlighted &&
	a.colorOverride == b.colorOverride && a.color == b.color &&
	a.selectedColor == b.selectedColor &&
	a.highlightedColor == b.highlightedColor &&
	a.materialColorValid == b.materialColorValid &&
	a.materialColor == b.materialColor &&
	a.databaseMaterialColorValid == b.databaseMaterialColorValid &&
	a.databaseMaterialColor == b.databaseMaterialColor &&
	std::fabs(a.transparency - b.transparency) <= 1.0e-6f &&
	a.drawMatrixValid == b.drawMatrixValid &&
	a.localToWorld.equals(b.localToWorld, 1.0e-6f);
}

static bool
sameClipPlanes(const SbPlane *a, size_t aCount, const SbPlane *b,
	size_t bCount)
{
    if (aCount != bCount)
	return false;
    for (size_t i = 0; i < aCount; ++i) {
	if (a[i] != b[i])
	    return false;
    }
    return true;
}

static SbColor
sourceColor(const RtSourceSnapshot &source, const BRLObolRtPickResult &hit)
{
    if (source.highlighted)
	return source.highlightedColor;
    if (source.selected)
	return source.selectedColor;
    if (source.colorOverride)
	return source.color;
    if (hit.detail.hasMaterialColor())
	return hit.detail.getMaterialColor();
    if (source.materialColorValid)
	return source.materialColor;
    if (source.databaseMaterialColorValid)
	return source.databaseMaterialColor;
    return source.color;
}

static unsigned char
colorByte(float value)
{
    const float clamped = std::max(0.0f, std::min(1.0f, value));
    return static_cast<unsigned char>(std::lround(clamped * 255.0f));
}

static unsigned int
clampWorkers(unsigned int workers)
{
    const unsigned int maximum = static_cast<unsigned int>(MAX_PSW);
    if (maximum == 0)
	return 1;
    return std::max(1u, std::min(workers, maximum));
}

static unsigned int
clampSamples(unsigned int samples)
{
    return std::max(1u, std::min(samples, 64u));
}

static bool
shouldCancel(const std::atomic_bool *cancelled)
{
    return cancelled && cancelled->load(std::memory_order_acquire);
}

struct BRLObolRtRenderer::Private {
    struct Source {
	/* Resources must outlive cache destruction: librt records each pointer in
	 * the rt_i and releases it while destroying the acceleration structure. */
	std::vector<std::unique_ptr<struct resource> > resources;
	std::unique_ptr<BRLObolRtPickCache> cache;
	RtSourceSnapshot snapshot;
    };

    std::vector<Source> sources;
    SbViewVolume viewVolume;
    SbPlane clipPlanes[2];
    size_t clipPlaneCount = 0;
    SbColor backgroundBottom = SbColor(0.0f, 0.0f, 0.0f);
    SbColor backgroundTop = SbColor(0.0f, 0.0f, 0.0f);
    uint64_t geometryRevision = 0;
    uint64_t presentationRevision = 0;
    SbMatrix viewAffine = SbMatrix::identity();
    SbMatrix viewProjection = SbMatrix::identity();
    SbBool lightingEnabled = TRUE;
    SbBool transparencyEnabled = TRUE;
    SbBool cameraReady = FALSE;
    SbBool presentationStateReady = FALSE;
};

BRLObolRtRenderSettings::BRLObolRtRenderSettings(void) :
    width(1),
    height(1),
    workers(1),
    samples(1),
    backgroundBottom(0.0f, 0.0f, 0.0f),
    backgroundTop(0.0f, 0.0f, 0.0f)
{
}

BRLObolRtRenderPlanes::BRLObolRtRenderPlanes(void)
{
    clear();
}

void
BRLObolRtRenderPlanes::clear(void)
{
    depth.clear();
    sourceIdentity.clear();
}

BRLObolRtSourceIdentity::BRLObolRtSourceIdentity(void)
{
    clear();
}

void
BRLObolRtSourceIdentity::clear(void)
{
    database = NULL;
    instanceKey = "";
    path = "";
    sourceRevision = 0;
}

BRLObolRtRenderStatus::BRLObolRtRenderStatus(void)
{
    clear();
}

void
BRLObolRtRenderStatus::clear(void)
{
    geometryRevision = 0;
    presentationRevision = 0;
    raysShot = 0;
    preparedSources = 0;
    width = 0;
    height = 0;
    complete = 0;
    cancelled = 0;
}

BRLObolRtRenderer::BRLObolRtRenderer(void) :
    p(new Private)
{
}

BRLObolRtRenderer::~BRLObolRtRenderer(void)
{
    delete p;
    p = NULL;
}

void
BRLObolRtRenderer::clear(void)
{
    if (!p)
	return;
    p->sources.clear();
    p->cameraReady = FALSE;
    p->presentationStateReady = FALSE;
    p->geometryRevision = 0;
    p->presentationRevision = 0;
}

SbBool
BRLObolRtRenderer::synchronize(BRLObolViewController *controller)
{
    if (!p || !controller || !controller->getCamera())
	return FALSE;

    const SbViewportRegion &region = controller->getViewportRegion();
    const SbVec2s size = region.getViewportSizePixels();
    if (size[0] <= 0 || size[1] <= 0)
	return FALSE;

    const float aspect = static_cast<float>(size[0]) /
	static_cast<float>(size[1]);
    const SbViewVolume nextViewVolume =
	controller->getCamera()->getViewVolume(aspect);
    SbMatrix nextViewAffine;
    SbMatrix nextViewProjection;
    nextViewVolume.getMatrices(nextViewAffine, nextViewProjection);
    const SbColor nextBackgroundBottom = controller->getBackgroundBottomColor();
    const SbColor nextBackgroundTop = controller->getBackgroundTopColor();
    SbPlane nextClipPlanes[2];
    const size_t nextClipPlaneCount =
	controller->getActiveClipPlanes(nextClipPlanes);
    const SbBool nextLightingEnabled = controller->isLightingEnabled();
    const SbBool nextTransparencyEnabled =
	controller->isTransparencyEnabled();
    const bool controllerPresentationChanged = !p->presentationStateReady ||
	!p->viewAffine.equals(nextViewAffine, 1.0e-6f) ||
	!p->viewProjection.equals(nextViewProjection, 1.0e-6f) ||
	p->backgroundBottom != nextBackgroundBottom ||
	p->backgroundTop != nextBackgroundTop ||
	!sameClipPlanes(p->clipPlanes, p->clipPlaneCount, nextClipPlanes,
	    nextClipPlaneCount) ||
	p->lightingEnabled != nextLightingEnabled ||
	p->transparencyEnabled != nextTransparencyEnabled;
    std::vector<RtSourceSnapshot> snapshots;
    const int sourceCount = controller->getDatabaseSourceCount();
    snapshots.reserve(sourceCount > 0 ? static_cast<size_t>(sourceCount) : 0);
    for (int i = 0; i < sourceCount; ++i) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (!source || !source->getDatabase() ||
	    source->path.getValue().getLength() == 0)
	    continue;

	RtSourceSnapshot snapshot;
	snapshot.key.database = source->getDatabase();
	snapshot.key.instanceKey = source->instanceKey.getValue();
	snapshot.key.path = source->path.getValue();
	snapshot.key.sourceRevision = static_cast<uint32_t>(
	    source->sourceRevision.getValue());
	snapshot.visible = source->visible.getValue();
	snapshot.selected = source->selected.getValue();
	snapshot.highlighted = source->highlighted.getValue();
	snapshot.colorOverride = source->colorOverride.getValue();
	snapshot.color = source->color.getValue();
	snapshot.selectedColor = source->selectedColor.getValue();
	snapshot.highlightedColor = source->highlightedColor.getValue();
	snapshot.materialColorValid = source->materialColorValid.getValue();
	snapshot.materialColor = source->materialColor.getValue();
	snapshot.databaseMaterialColorValid =
	    source->databaseMaterialColorValid.getValue();
	snapshot.databaseMaterialColor = source->databaseMaterialColor.getValue();
	snapshot.transparency = source->transparency.getValue();
	snapshot.drawMatrixValid = source->drawMatrixValid.getValue();
	if (snapshot.drawMatrixValid) {
	    snapshot.localToWorld = source->drawMatrix.getValue();
	    snapshot.worldToLocal = snapshot.localToWorld.inverse();
	}
	snapshot.localClipPlaneCount = nextClipPlaneCount;
	for (size_t planeIndex = 0; planeIndex < nextClipPlaneCount;
	    ++planeIndex) {
	    snapshot.localClipPlanes[planeIndex] = nextClipPlanes[planeIndex];
	    if (snapshot.drawMatrixValid)
		snapshot.localClipPlanes[planeIndex].transform(
		    snapshot.worldToLocal);
	}
	snapshots.push_back(snapshot);
    }

	/* Prepare every missing source before moving any retained entry.  A
	 * preparation failure must leave the last coherent retained scene intact. */
	std::vector<RtPreparedCache> prepared;
	prepared.reserve(snapshots.size());
	for (size_t i = 0; i < snapshots.size(); ++i) {
	    const RtSourceSnapshot &snapshot = snapshots[i];
	    const std::vector<Private::Source>::const_iterator found = std::find_if(
		p->sources.begin(), p->sources.end(),
		[&snapshot](const Private::Source &candidate) {
		    return candidate.snapshot.key == snapshot.key;
		});
	    if (found != p->sources.end())
		continue;

	    RtPreparedCache entry;
	    entry.snapshotIndex = i;
	    entry.cache.reset(new BRLObolRtPickCache);
	    std::vector<SbString> paths(1, snapshot.key.path);
	    if (!entry.cache->prepare(snapshot.key.database, paths))
		return FALSE;
	    prepared.push_back(std::move(entry));
	}

	std::vector<Private::Source> next;
	next.reserve(snapshots.size());
    bool geometryChanged = snapshots.size() != p->sources.size();
	bool presentationChanged = geometryChanged || controllerPresentationChanged;
	for (size_t snapshotIndex = 0; snapshotIndex < snapshots.size();
	    ++snapshotIndex) {
	const RtSourceSnapshot &snapshot = snapshots[snapshotIndex];
	std::vector<Private::Source>::iterator found = std::find_if(
	    p->sources.begin(), p->sources.end(),
	    [&snapshot](const Private::Source &candidate) {
		return candidate.snapshot.key == snapshot.key;
	    });
	if (found != p->sources.end()) {
	    presentationChanged = presentationChanged ||
		!samePresentation(found->snapshot, snapshot);
	    found->snapshot = snapshot;
	    next.push_back(std::move(*found));
	    continue;
	}

	Private::Source source;
	source.snapshot = snapshot;
	std::vector<RtPreparedCache>::iterator preparedCache =
	    std::find_if(prepared.begin(), prepared.end(), [snapshotIndex](
		const RtPreparedCache &candidate) {
		return candidate.snapshotIndex == snapshotIndex;
	    });
	if (preparedCache == prepared.end())
	    return FALSE;
	source.cache = std::move(preparedCache->cache);
	geometryChanged = TRUE;
	presentationChanged = TRUE;
	next.push_back(std::move(source));
    }

	p->sources.swap(next);
    p->viewVolume = nextViewVolume;
    p->viewAffine = nextViewAffine;
    p->viewProjection = nextViewProjection;
    p->backgroundBottom = nextBackgroundBottom;
    p->backgroundTop = nextBackgroundTop;
    p->clipPlaneCount = nextClipPlaneCount;
    for (size_t planeIndex = 0; planeIndex < p->clipPlaneCount;
	++planeIndex)
	p->clipPlanes[planeIndex] = nextClipPlanes[planeIndex];
    p->lightingEnabled = nextLightingEnabled;
    p->transparencyEnabled = nextTransparencyEnabled;
    p->cameraReady = TRUE;
    p->presentationStateReady = TRUE;
    if (geometryChanged)
	++p->geometryRevision;
    if (presentationChanged)
	++p->presentationRevision;
    return TRUE;
}

SbBool
BRLObolRtRenderer::render(const BRLObolRtRenderSettings &settings,
	std::vector<unsigned char> &rgb,
	BRLObolRtRenderStatus *status,
	const std::atomic_bool *cancelled)
{
    return renderRows(settings, rgb, 0, settings.height, status, cancelled);
}

SbBool
BRLObolRtRenderer::renderWithPlanes(const BRLObolRtRenderSettings &settings,
	std::vector<unsigned char> &rgb, BRLObolRtRenderPlanes &planes,
	BRLObolRtRenderStatus *status, const std::atomic_bool *cancelled)
{
    return renderRowsWithPlanes(settings, rgb, planes, 0, settings.height,
	status, cancelled);
}

SbBool
BRLObolRtRenderer::renderRows(const BRLObolRtRenderSettings &settings,
	std::vector<unsigned char> &rgb, unsigned int firstRow,
	unsigned int rowCount, BRLObolRtRenderStatus *status,
	const std::atomic_bool *cancelled)
{
    return renderRowsInternal(settings, rgb, NULL, firstRow, rowCount,
	status, cancelled);
}

SbBool
BRLObolRtRenderer::renderRowsWithPlanes(
	const BRLObolRtRenderSettings &settings,
	std::vector<unsigned char> &rgb, BRLObolRtRenderPlanes &planes,
	unsigned int firstRow, unsigned int rowCount,
	BRLObolRtRenderStatus *status, const std::atomic_bool *cancelled)
{
    return renderRowsInternal(settings, rgb, &planes, firstRow, rowCount,
	status, cancelled);
}

SbBool
BRLObolRtRenderer::renderRowsInternal(
	const BRLObolRtRenderSettings &settings,
	std::vector<unsigned char> &rgb, BRLObolRtRenderPlanes *planes,
	unsigned int firstRow, unsigned int rowCount,
	BRLObolRtRenderStatus *status, const std::atomic_bool *cancelled)
{
    BRLObolRtRenderStatus localStatus;
    if (status)
	*status = localStatus;
    if (!p || !p->cameraReady || settings.width == 0 || settings.height == 0) {
	return FALSE;
    }
    if (!rowCount || firstRow >= settings.height) {
	return FALSE;
    }
    const unsigned int lastRow = std::min(settings.height,
	firstRow + std::min(rowCount, settings.height - firstRow));

    const size_t pixelCount = static_cast<size_t>(settings.width) *
	static_cast<size_t>(settings.height);
    if (pixelCount > std::numeric_limits<size_t>::max() / 3u)
	return FALSE;

    const unsigned int workerCount = clampWorkers(settings.workers);
    const unsigned int sampleCount = clampSamples(settings.samples);
    const unsigned int activeWorkerCount = std::min(workerCount,
	lastRow - firstRow);
    const bool transparentCompositionRequired = p->transparencyEnabled &&
	std::any_of(p->sources.begin(), p->sources.end(),
	    [](const Private::Source &source) {
		return source.snapshot.transparency > 0.0f &&
		    source.snapshot.transparency < 1.0f;
	    });
    if (rgb.size() != pixelCount * 3u || firstRow == 0) {
	rgb.assign(pixelCount * 3u, 0u);
	for (unsigned int y = 0; y < settings.height; ++y) {
	    const float t = settings.height > 1 ?
		static_cast<float>(y) / static_cast<float>(settings.height - 1) :
		0.0f;
	    const SbColor background = settings.backgroundBottom * (1.0f - t) +
		settings.backgroundTop * t;
	    for (unsigned int x = 0; x < settings.width; ++x) {
		unsigned char *pixel = &rgb[(static_cast<size_t>(y) * settings.width + x) * 3u];
		pixel[0] = colorByte(background[0]);
		pixel[1] = colorByte(background[1]);
		pixel[2] = colorByte(background[2]);
	    }
	}
    }
    if (planes && (planes->depth.size() != pixelCount ||
	planes->sourceIdentity.size() != pixelCount || firstRow == 0)) {
	planes->depth.assign(pixelCount, 1.0f);
	planes->sourceIdentity.assign(pixelCount, 0u);
    }

    for (Private::Source &source : p->sources) {
	if (!source.cache || !source.cache->isReady())
	    continue;
	while (source.resources.size() < activeWorkerCount) {
	    std::unique_ptr<struct resource> resource(new struct resource);
	    *resource = RT_RESOURCE_INIT_ZERO;
	    if (!source.cache->initializeResource(resource.get(),
		static_cast<int>(source.resources.size())))
		return FALSE;
	    source.resources.push_back(std::move(resource));
	}
    }

    std::atomic<unsigned int> nextRow(firstRow);
    std::atomic<uint64_t> raysShot(0);
    std::atomic_bool stop(false);
    const auto traceLine = [&](const SbLine &line, unsigned int workerIndex,
	BRLObolRtPickResult &nearest, size_t &nearestSourceIndex,
	std::vector<RtTraceHit> *allHits = NULL) {
	nearest.clear();
	nearestSourceIndex = p->sources.size();
	if (allHits)
	    allHits->clear();
	for (size_t sourceIndex = 0; sourceIndex < p->sources.size();
	    ++sourceIndex) {
	    Private::Source &source = p->sources[sourceIndex];
	    if (!source.snapshot.visible ||
		(p->transparencyEnabled && source.snapshot.transparency >= 1.0f) ||
		!source.cache ||
		workerIndex >= source.resources.size())
		continue;
	    SbVec3f rayOrigin = line.getPosition();
	    SbVec3f rayDirection = line.getDirection();
	    if (source.snapshot.drawMatrixValid) {
		source.snapshot.worldToLocal.multVecMatrix(rayOrigin, rayOrigin);
		source.snapshot.worldToLocal.multDirMatrix(rayDirection, rayDirection);
	    }
	    if (rayDirection.length() <= 0.0f)
		continue;
	    rayDirection.normalize();
	    BRLObolRtPickResult candidate;
	    if (source.cache->pickRay(candidate, rayOrigin, rayDirection,
		source.snapshot.localClipPlanes,
		source.snapshot.localClipPlaneCount,
		source.resources[workerIndex].get()) && candidate.hit) {
		if (source.snapshot.drawMatrixValid) {
		    SbVec3f worldPoint;
		    SbVec3f worldNormal;
		    SbMatrix normalMatrix =
			source.snapshot.localToWorld.inverse().transpose();
		    source.snapshot.localToWorld.multVecMatrix(candidate.point,
			worldPoint);
		    normalMatrix.multDirMatrix(candidate.normal, worldNormal);
		    if (worldNormal.length() > 0.0f)
			worldNormal.normalize();
		    candidate.point = worldPoint;
		    candidate.normal = worldNormal;
		    candidate.distance =
			(worldPoint - line.getPosition()).length();
		    candidate.detail.setModelPoint(worldPoint);
		}
		if (allHits) {
		    RtTraceHit hit;
		    hit.pick = candidate;
		    hit.sourceIndex = sourceIndex;
		    allHits->push_back(hit);
		} else if (nearestSourceIndex == p->sources.size() ||
		    candidate.distance < nearest.distance) {
		    nearest = candidate;
		    nearestSourceIndex = sourceIndex;
		}
	    }
	    raysShot.fetch_add(1, std::memory_order_relaxed);
	}
	if (allHits && !allHits->empty()) {
	    std::stable_sort(allHits->begin(), allHits->end(),
		[](const RtTraceHit &a, const RtTraceHit &b) {
		    return a.pick.distance < b.pick.distance;
		});
	    nearest = allHits->front().pick;
	    nearestSourceIndex = allHits->front().sourceIndex;
	}
	return nearestSourceIndex < p->sources.size();
    };
    const auto shadeHit = [&](const SbLine &line,
	const BRLObolRtPickResult &hit, const RtSourceSnapshot &source) {
	SbVec3f normal = hit.normal;
	if (normal.length() > 0.0f) {
	    normal.normalize();
	    const SbVec3f toEye = -line.getDirection();
	    if (normal.dot(toEye) < 0.0f)
		normal = -normal;
	}
	const float diffuse = normal.length() > 0.0f ?
	    std::max(0.0f, normal.dot(-line.getDirection())) : 1.0f;
	const float illumination = p->lightingEnabled ?
	    0.20f + 0.80f * diffuse : 1.0f;
	return sourceColor(source, hit) * illumination;
    };
    const auto renderWorker = [&](unsigned int workerIndex) {
	std::vector<RtTraceHit> transparentHits;
	if (transparentCompositionRequired)
	    transparentHits.reserve(p->sources.size());
	while (!stop.load(std::memory_order_acquire)) {
	    if (shouldCancel(cancelled)) {
		stop.store(true, std::memory_order_release);
		break;
	    }
	    const unsigned int y = nextRow.fetch_add(1, std::memory_order_relaxed);
	    if (y >= lastRow)
		break;
	    for (unsigned int x = 0; x < settings.width; ++x) {
		if (shouldCancel(cancelled)) {
		    stop.store(true, std::memory_order_release);
		    break;
		}
		float accumulated[3] = {0.0f, 0.0f, 0.0f};
		for (unsigned int sample = 0; sample < sampleCount; ++sample) {
		    const unsigned int grid = static_cast<unsigned int>(std::ceil(
			std::sqrt(static_cast<double>(sampleCount))));
		    const float offsetX = (static_cast<float>(sample % grid) + 0.5f) /
			static_cast<float>(grid);
		    const float offsetY = (static_cast<float>(sample / grid) + 0.5f) /
			static_cast<float>(grid);
		    const float nx = (static_cast<float>(x) + offsetX) /
			static_cast<float>(settings.width);
		    const float ny = (static_cast<float>(y) + offsetY) /
			static_cast<float>(settings.height);
		    SbLine line;
		    p->viewVolume.projectPointToLine(SbVec2f(nx, ny), line);

		    BRLObolRtPickResult nearest;
		    size_t nearestSourceIndex = p->sources.size();
		    const SbBool hasNearest = traceLine(line, workerIndex, nearest,
			nearestSourceIndex, transparentCompositionRequired ?
			&transparentHits : NULL) ? TRUE : FALSE;
		    const RtSourceSnapshot *nearestSource = hasNearest ?
			&p->sources[nearestSourceIndex].snapshot : NULL;

		    const float backgroundT = settings.height > 1 ?
			static_cast<float>(y) / static_cast<float>(settings.height - 1) : 0.0f;
		    SbColor shaded = settings.backgroundBottom * (1.0f - backgroundT) +
			settings.backgroundTop * backgroundT;
		    if (transparentCompositionRequired) {
			for (std::vector<RtTraceHit>::reverse_iterator hit =
				transparentHits.rbegin(); hit != transparentHits.rend();
				++hit) {
			    const RtSourceSnapshot &source =
				p->sources[hit->sourceIndex].snapshot;
			    const float alpha = std::max(0.0f, std::min(1.0f,
				1.0f - source.transparency));
			    const SbColor lit = shadeHit(line, hit->pick, source);
			    shaded = lit * alpha + shaded * (1.0f - alpha);
			}
		    } else if (nearestSource) {
			shaded = shadeHit(line, nearest, *nearestSource);
		    }
		    accumulated[0] += shaded[0];
		    accumulated[1] += shaded[1];
		    accumulated[2] += shaded[2];
		}
		const float inverseSamples = 1.0f / static_cast<float>(sampleCount);
		unsigned char *pixel = &rgb[(static_cast<size_t>(y) * settings.width + x) * 3u];
		    pixel[0] = colorByte(accumulated[0] * inverseSamples);
		    pixel[1] = colorByte(accumulated[1] * inverseSamples);
		    pixel[2] = colorByte(accumulated[2] * inverseSamples);

		    if (planes) {
			const float centerX = (static_cast<float>(x) + 0.5f) /
			    static_cast<float>(settings.width);
			const float centerY = (static_cast<float>(y) + 0.5f) /
			    static_cast<float>(settings.height);
			SbLine centerLine;
			p->viewVolume.projectPointToLine(SbVec2f(centerX, centerY),
			    centerLine);
			BRLObolRtPickResult centerHit;
			size_t centerSourceIndex = p->sources.size();
			if (traceLine(centerLine, workerIndex, centerHit,
				centerSourceIndex)) {
			    SbVec3f screenPoint;
			    p->viewVolume.projectToScreen(centerHit.point, screenPoint);
			    const size_t planeIndex =
				static_cast<size_t>(y) * settings.width + x;
			    planes->depth[planeIndex] = std::max(0.0f, std::min(1.0f,
				screenPoint[2]));
			    planes->sourceIdentity[planeIndex] =
				static_cast<uint32_t>(centerSourceIndex + 1u);
			}
		    }
		}
	}
    };

    std::vector<std::thread> workers;
    workers.reserve(activeWorkerCount > 0 ? activeWorkerCount - 1u : 0u);
    for (unsigned int workerIndex = 1; workerIndex < activeWorkerCount;
	++workerIndex)
	workers.push_back(std::thread(renderWorker, workerIndex));
    renderWorker(0);
    for (std::thread &thread : workers)
	thread.join();

    localStatus.geometryRevision = p->geometryRevision;
    localStatus.presentationRevision = p->presentationRevision;
    localStatus.raysShot = raysShot.load(std::memory_order_relaxed);
    localStatus.preparedSources = static_cast<unsigned int>(p->sources.size());
    localStatus.width = settings.width;
    localStatus.height = settings.height;
    localStatus.cancelled = stop.load(std::memory_order_acquire) ? 1 : 0;
    localStatus.complete = localStatus.cancelled ? 0 : 1;
    if (status)
	*status = localStatus;
    return localStatus.complete ? TRUE : FALSE;
}

unsigned int
BRLObolRtRenderer::getPreparedSourceCount(void) const
{
    return p ? static_cast<unsigned int>(p->sources.size()) : 0u;
}

uint64_t
BRLObolRtRenderer::getGeometryRevision(void) const
{
    return p ? p->geometryRevision : 0u;
}

uint64_t
BRLObolRtRenderer::getPresentationRevision(void) const
{
    return p ? p->presentationRevision : 0u;
}

SbBool
BRLObolRtRenderer::getSourceIdentity(uint32_t identifier,
	BRLObolRtSourceIdentity &identity) const
{
    identity.clear();
    if (!p || identifier == 0u)
	return FALSE;
    const size_t sourceIndex = static_cast<size_t>(identifier - 1u);
    if (sourceIndex >= p->sources.size())
	return FALSE;
    const RtSourceKey &key = p->sources[sourceIndex].snapshot.key;
    identity.database = key.database;
    identity.instanceKey = key.instanceKey;
    identity.path = key.path;
    identity.sourceRevision = key.sourceRevision;
    return TRUE;
}
