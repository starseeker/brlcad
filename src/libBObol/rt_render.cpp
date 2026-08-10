/*                    R T _ R E N D E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libBObol/rt_render.cpp
 *
 * A retained librt image renderer.  This deliberately does not use Obol or
 * Coin objects while rendering: those are sampled once by synchronize().
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BPickDetail.h"
#include "BObol/BRtRender.h"
#include "BObol/BViewController.h"

#include "raytrace.h"

#include <Inventor/SbLine.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbPlane.h>
#include <Inventor/SbViewVolume.h>
#include <Inventor/actions/SoGetMatrixAction.h>
#include <Inventor/actions/SoSearchAction.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoLight.h>
#include <Inventor/nodes/SoPointLight.h>
#include <Inventor/nodes/SoSpotLight.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functional>
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
    SbString databaseMaterialShader;
    float transparency = 0.0f;
    SbBool drawMatrixValid = FALSE;
    SbMatrix localToWorld = SbMatrix::identity();
    SbMatrix worldToLocal = SbMatrix::identity();
    SbPlane localClipPlanes[2];
    size_t localClipPlaneCount = 0;
};

struct RtPreparedCache {
    size_t snapshotIndex = 0;
    std::unique_ptr<BObolRtPickCache> cache;
};

struct RtTraceHit {
    BObolRtPickResult pick;
    size_t sourceIndex = 0;
};

struct RtLightSnapshot {
    enum Kind {
	DIRECTIONAL = 0,
	POINT = 1,
	SPOT = 2
    };

    Kind kind = DIRECTIONAL;
    SbColor color = SbColor(1.0f, 1.0f, 1.0f);
    float intensity = 1.0f;
    SbVec3f position = SbVec3f(0.0f, 0.0f, 0.0f);
    SbVec3f direction = SbVec3f(0.0f, 0.0f, -1.0f);
    float dropOffRate = 0.0f;
    float cutOffAngle = static_cast<float>(M_PI) * 0.25f;
};

static bool
sameLight(const RtLightSnapshot &a, const RtLightSnapshot &b)
{
    return a.kind == b.kind && a.color == b.color &&
	std::fabs(a.intensity - b.intensity) <= 1.0e-6f &&
	(a.position - b.position).sqrLength() <= 1.0e-12f &&
	(a.direction - b.direction).sqrLength() <= 1.0e-12f &&
	std::fabs(a.dropOffRate - b.dropOffRate) <= 1.0e-6f &&
	std::fabs(a.cutOffAngle - b.cutOffAngle) <= 1.0e-6f;
}

static bool
sameLights(const std::vector<RtLightSnapshot> &a,
	const std::vector<RtLightSnapshot> &b)
{
    if (a.size() != b.size())
	return false;
    for (size_t i = 0; i < a.size(); ++i) {
	if (!sameLight(a[i], b[i]))
	    return false;
    }
    return true;
}

static std::vector<RtLightSnapshot>
controllerAuthoredLights(BObolViewController *controller,
	const SbViewportRegion &viewportRegion, size_t maximumAuthoredLights)
{
    std::vector<RtLightSnapshot> lights;

    if (!controller || !maximumAuthoredLights)
	return lights;

    /* Search both the geometry scene root and the sibling in-scene light group
     * (which holds database-derived lights, and is not part of the geometry
     * subtree).  Its light positions/directions are already world-space, so the
     * SoGetMatrixAction over that subtree yields identity and leaves them
     * unchanged. */
    SoNode *roots[2] = {
	controller->getRenderSceneRoot(),
	controller->getSceneLightsRoot()
    };
    for (int r = 0; r < 2 && lights.size() < maximumAuthoredLights; ++r) {
	SoNode *root = roots[r];
	if (!root)
	    continue;
	SoSearchAction search;
	search.setType(SoLight::getClassTypeId(), TRUE);
	search.setInterest(SoSearchAction::ALL);
	search.apply(root);
	SoPathList &paths = search.getPaths();
	for (int i = 0; i < paths.getLength() &&
	    lights.size() < maximumAuthoredLights; ++i) {
	    SoPath *path = paths[i];
	    SoNode *tail = path ? path->getTail() : NULL;
	    if (!tail || !tail->isOfType(SoLight::getClassTypeId()))
		continue;
	    SoLight *light = static_cast<SoLight *>(tail);
	    if (!light->on.getValue() || light->intensity.getValue() <= 0.0f)
		continue;

	    SoGetMatrixAction matrixAction(viewportRegion);
	    matrixAction.apply(path);
	    const SbMatrix &matrix = matrixAction.getMatrix();
	    RtLightSnapshot snapshot;
	    snapshot.color = light->color.getValue();
	    snapshot.intensity = std::max(0.0f,
		std::min(1.0f, light->intensity.getValue()));
	    if (tail->isOfType(SoSpotLight::getClassTypeId())) {
		SoSpotLight *spot = static_cast<SoSpotLight *>(tail);
		snapshot.kind = RtLightSnapshot::SPOT;
		matrix.multVecMatrix(spot->location.getValue(), snapshot.position);
		matrix.multDirMatrix(spot->direction.getValue(), snapshot.direction);
		snapshot.dropOffRate = std::max(0.0f,
		    std::min(1.0f, spot->dropOffRate.getValue()));
		snapshot.cutOffAngle = std::max(0.0f,
		    std::min(static_cast<float>(M_PI) * 0.5f,
			spot->cutOffAngle.getValue()));
	    } else if (tail->isOfType(SoPointLight::getClassTypeId())) {
		SoPointLight *point = static_cast<SoPointLight *>(tail);
		snapshot.kind = RtLightSnapshot::POINT;
		matrix.multVecMatrix(point->location.getValue(), snapshot.position);
	    } else if (tail->isOfType(SoDirectionalLight::getClassTypeId())) {
		SoDirectionalLight *directional =
		    static_cast<SoDirectionalLight *>(tail);
		snapshot.kind = RtLightSnapshot::DIRECTIONAL;
		matrix.multDirMatrix(directional->direction.getValue(),
		    snapshot.direction);
	    } else {
		continue;
	    }
	    if (snapshot.direction.length() > 0.0f)
		snapshot.direction.normalize();
	    lights.push_back(snapshot);
	}
    }
    return lights;
}

static std::vector<RtLightSnapshot>
controllerCameraLights(BObolViewController *controller)
{
    std::vector<RtLightSnapshot> lights;
    if (!controller)
	return lights;
    std::vector<BObolSceneLightRealization> realized;
    controller->getCameraLights(realized);
    lights.reserve(realized.size());
    for (size_t i = 0; i < realized.size(); i++) {
	RtLightSnapshot snapshot;
	snapshot.kind = RtLightSnapshot::DIRECTIONAL;
	snapshot.color = realized[i].color;
	snapshot.intensity = std::max(0.0f,
	    std::min(1.0f, realized[i].intensity));
	snapshot.direction = realized[i].direction;
	if (snapshot.direction.length() > 0.0f)
	    snapshot.direction.normalize();
	lights.push_back(snapshot);
    }
    return lights;
}

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
	a.databaseMaterialShader == b.databaseMaterialShader &&
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

struct RtMaterialModel {
    float ambient = 0.20f;
    float diffuse = 0.80f;
    float specular = 0.0f;
    float shine = 10.0f;
    float reflectance = 0.0f;
    float transmission = 0.0f;
    float refractiveIndex = 1.0f;
};

static bool
shaderTokenEqual(const char *token, size_t length, const char *expected)
{
    if (!token || !expected || std::strlen(expected) != length)
	return false;
    for (size_t i = 0; i < length; ++i) {
	if (std::tolower(static_cast<unsigned char>(token[i])) !=
	    std::tolower(static_cast<unsigned char>(expected[i])))
	    return false;
    }
    return true;
}

static const char *
shaderStage(const char *shader, const char *expected)
{
    if (!shader || !expected)
	return NULL;
    const char *cursor = shader;
    while (*cursor) {
	while (*cursor &&
	    !std::isalnum(static_cast<unsigned char>(*cursor)) &&
	    *cursor != '_')
	    ++cursor;
	const char *token = cursor;
	while (*cursor &&
	    (std::isalnum(static_cast<unsigned char>(*cursor)) ||
	     *cursor == '_'))
	    ++cursor;
	if (token == cursor)
	    break;
	if (shaderTokenEqual(token, static_cast<size_t>(cursor - token),
		expected))
	    return token;
    }
    return NULL;
}

static bool
shaderParameter(const char *shader, const char *name, float &value)
{
    if (!shader || !name)
	return false;

    const char *cursor = shader;
    while (*cursor) {
	while (*cursor &&
	    !std::isalnum(static_cast<unsigned char>(*cursor)) &&
	    *cursor != '_')
	    ++cursor;
	const char *token = cursor;
	while (*cursor &&
	    (std::isalnum(static_cast<unsigned char>(*cursor)) ||
	     *cursor == '_'))
	    ++cursor;
	if (token == cursor)
	    break;
	if (!shaderTokenEqual(token, static_cast<size_t>(cursor - token), name))
	    continue;

	while (*cursor && std::isspace(static_cast<unsigned char>(*cursor)))
	    ++cursor;
	if (*cursor == '=') {
	    ++cursor;
	    while (*cursor && std::isspace(static_cast<unsigned char>(*cursor)))
		++cursor;
	}
	char *end = NULL;
	const double parsed = std::strtod(cursor, &end);
	if (end != cursor && std::isfinite(parsed)) {
	    value = static_cast<float>(parsed);
	    return true;
	}
    }

    return false;
}

static const char *
effectiveShader(const RtSourceSnapshot &source,
	const BObolRtPickResult &hit)
{
    const SbString &primitiveShader = hit.detail.getMaterialShader();
    return primitiveShader.getLength() > 0 ? primitiveShader.getString() :
	source.databaseMaterialShader.getString();
}

static RtMaterialModel
materialModel(const RtSourceSnapshot &source, const BObolRtPickResult &hit)
{
    const char *shader = effectiveShader(source, hit);
    RtMaterialModel model;
    const char *plasticStage = shaderStage(shader, "plastic");
    const char *mirrorStage = shaderStage(shader, "mirror");
    const char *glassStage = shaderStage(shader, "glass");
    const bool plastic = plasticStage != NULL;
    const bool mirror = mirrorStage != NULL;
    const bool glass = glassStage != NULL;
    if (!plastic && !mirror && !glass)
	return model;
    const char *materialStage = glass ? glassStage :
	(mirror ? mirrorStage : plasticStage);

    /* Match the established BRL-CAD plastic-family vocabulary while keeping
     * this retained interactive model deliberately bounded. */
    model.ambient = 0.10f;
    if (mirror) {
	model.diffuse = 0.40f;
	model.specular = 0.60f;
	model.shine = 4.0f;
	model.reflectance = 0.75f;
	model.refractiveIndex = 1.65f;
    } else if (glass) {
	model.diffuse = 0.30f;
	model.specular = 0.70f;
	model.shine = 4.0f;
	model.reflectance = 0.10f;
	model.transmission = 0.80f;
	model.refractiveIndex = 1.65f;
    } else {
	model.diffuse = 0.30f;
	model.specular = 0.70f;
    }
    float parsed = 0.0f;
    if (shaderParameter(materialStage, "diffuse", parsed) ||
	shaderParameter(materialStage, "di", parsed))
	model.diffuse = std::max(0.0f, std::min(1.0f, parsed));
    if (shaderParameter(materialStage, "specular", parsed) ||
	shaderParameter(materialStage, "sp", parsed))
	model.specular = std::max(0.0f, std::min(1.0f, parsed));
    if (shaderParameter(materialStage, "shine", parsed) ||
	shaderParameter(materialStage, "sh", parsed))
	model.shine = std::max(1.0f, std::min(256.0f, parsed));
    if (shaderParameter(materialStage, "reflect", parsed) ||
	shaderParameter(materialStage, "re", parsed))
	model.reflectance = std::max(0.0f, std::min(1.0f, parsed));
    if (shaderParameter(materialStage, "transmit", parsed) ||
	shaderParameter(materialStage, "tr", parsed))
	model.transmission = std::max(0.0f, std::min(1.0f, parsed));
    if (shaderParameter(materialStage, "ri", parsed))
	model.refractiveIndex = std::max(1.0f, std::min(4.0f, parsed));

    const float secondaryEnergy = model.reflectance + model.transmission;
    if (secondaryEnergy > 1.0f) {
	model.reflectance /= secondaryEnergy;
	model.transmission /= secondaryEnergy;
    }
    return model;
}

static SbColor
sourceColor(const RtSourceSnapshot &source, const BObolRtPickResult &hit)
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

struct BObolRtRenderer::Private {
    struct Source {
	/* Resources must outlive cache destruction: librt records each pointer in
	 * the rt_i and releases it while destroying the acceleration structure. */
	std::vector<std::unique_ptr<struct resource> > resources;
	std::unique_ptr<BObolRtPickCache> cache;
	RtSourceSnapshot snapshot;
    };

    std::vector<Source> sources;
    std::vector<RtLightSnapshot> cameraLights;
    std::vector<RtLightSnapshot> authoredLights;
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
    float ambientIntensity = 0.18f;
    SbBool transparencyEnabled = TRUE;
    SbBool cameraReady = FALSE;
    SbBool presentationStateReady = FALSE;
};

BObolRtRenderSettings::BObolRtRenderSettings(void) :
    width(1),
    height(1),
    workers(1),
    samples(1),
    backgroundBottom(0.0f, 0.0f, 0.0f),
    backgroundTop(0.0f, 0.0f, 0.0f)
{
}

BObolRtRenderPlanes::BObolRtRenderPlanes(void)
{
    clear();
}

void
BObolRtRenderPlanes::clear(void)
{
    depth.clear();
    sourceIdentity.clear();
}

BObolRtSourceIdentity::BObolRtSourceIdentity(void)
{
    clear();
}

void
BObolRtSourceIdentity::clear(void)
{
    database = NULL;
    instanceKey = "";
    path = "";
    sourceRevision = 0;
}

BObolRtRenderStatus::BObolRtRenderStatus(void)
{
    clear();
}

void
BObolRtRenderStatus::clear(void)
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

BObolRtRenderer::BObolRtRenderer(void) :
    p(new Private)
{
}

BObolRtRenderer::~BObolRtRenderer(void)
{
    delete p;
    p = NULL;
}

void
BObolRtRenderer::clear(void)
{
    if (!p)
	return;
    p->sources.clear();
    p->cameraLights.clear();
    p->authoredLights.clear();
    p->cameraReady = FALSE;
    p->presentationStateReady = FALSE;
    p->geometryRevision = 0;
    p->presentationRevision = 0;
}

SbBool
BObolRtRenderer::synchronize(BObolViewController *controller)
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
    const float nextAmbientIntensity =
	controller->getLightingAmbientIntensity();
    const SbBool nextTransparencyEnabled =
	controller->isTransparencyEnabled();
    const std::vector<RtLightSnapshot> nextCameraLights =
	controllerCameraLights(controller);
    /* Match the interactive renderer's eight-light contract exactly.  Camera
     * lights are traversed first and therefore have priority; database lights
     * consume only the remaining slots. */
    const size_t maximumRtLights = 8u;
    const size_t maximumAuthoredLights = nextCameraLights.size() <
	maximumRtLights ? maximumRtLights - nextCameraLights.size() : 0u;
    const std::vector<RtLightSnapshot> nextAuthoredLights =
	controllerAuthoredLights(controller, region, maximumAuthoredLights);
    const bool controllerPresentationChanged = !p->presentationStateReady ||
	!p->viewAffine.equals(nextViewAffine, 1.0e-6f) ||
	!p->viewProjection.equals(nextViewProjection, 1.0e-6f) ||
	p->backgroundBottom != nextBackgroundBottom ||
	p->backgroundTop != nextBackgroundTop ||
	!sameClipPlanes(p->clipPlanes, p->clipPlaneCount, nextClipPlanes,
	    nextClipPlaneCount) ||
	p->lightingEnabled != nextLightingEnabled ||
	std::fabs(p->ambientIntensity - nextAmbientIntensity) > 1.0e-6f ||
	p->transparencyEnabled != nextTransparencyEnabled ||
	!sameLights(p->cameraLights, nextCameraLights) ||
	!sameLights(p->authoredLights, nextAuthoredLights);
    std::vector<RtSourceSnapshot> snapshots;
    const std::vector<SoBRLDatabaseSource *> renderSources =
	controller->getRenderDatabaseSources();
    snapshots.reserve(renderSources.size());
    for (size_t i = 0; i < renderSources.size(); ++i) {
	SoBRLDatabaseSource *source = renderSources[i];
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
	snapshot.databaseMaterialShader =
	    source->databaseMaterialShader.getValue();
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
	    entry.cache.reset(new BObolRtPickCache);
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
    p->ambientIntensity = nextAmbientIntensity;
    p->transparencyEnabled = nextTransparencyEnabled;
    p->cameraLights = nextCameraLights;
    p->authoredLights = nextAuthoredLights;
    p->cameraReady = TRUE;
    p->presentationStateReady = TRUE;
    if (geometryChanged)
	++p->geometryRevision;
    if (presentationChanged)
	++p->presentationRevision;
    return TRUE;
}

SbBool
BObolRtRenderer::render(const BObolRtRenderSettings &settings,
	std::vector<unsigned char> &rgb,
	BObolRtRenderStatus *status,
	const std::atomic_bool *cancelled)
{
    return renderRows(settings, rgb, 0, settings.height, status, cancelled);
}

SbBool
BObolRtRenderer::renderWithPlanes(const BObolRtRenderSettings &settings,
	std::vector<unsigned char> &rgb, BObolRtRenderPlanes &planes,
	BObolRtRenderStatus *status, const std::atomic_bool *cancelled)
{
    return renderRowsWithPlanes(settings, rgb, planes, 0, settings.height,
	status, cancelled);
}

SbBool
BObolRtRenderer::renderRows(const BObolRtRenderSettings &settings,
	std::vector<unsigned char> &rgb, unsigned int firstRow,
	unsigned int rowCount, BObolRtRenderStatus *status,
	const std::atomic_bool *cancelled)
{
    return renderRowsInternal(settings, rgb, NULL, firstRow, rowCount,
	status, cancelled);
}

SbBool
BObolRtRenderer::renderRowsWithPlanes(
	const BObolRtRenderSettings &settings,
	std::vector<unsigned char> &rgb, BObolRtRenderPlanes &planes,
	unsigned int firstRow, unsigned int rowCount,
	BObolRtRenderStatus *status, const std::atomic_bool *cancelled)
{
    return renderRowsInternal(settings, rgb, &planes, firstRow, rowCount,
	status, cancelled);
}

SbBool
BObolRtRenderer::renderRowsInternal(
	const BObolRtRenderSettings &settings,
	std::vector<unsigned char> &rgb, BObolRtRenderPlanes *planes,
	unsigned int firstRow, unsigned int rowCount,
	BObolRtRenderStatus *status, const std::atomic_bool *cancelled)
{
    BObolRtRenderStatus localStatus;
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
	BObolRtPickResult &nearest, size_t &nearestSourceIndex,
	std::vector<RtTraceHit> *allHits, float minimumDistance) {
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
	    float localMinimumDistance = minimumDistance;
	    if (source.snapshot.drawMatrixValid) {
		source.snapshot.worldToLocal.multVecMatrix(rayOrigin, rayOrigin);
		source.snapshot.worldToLocal.multDirMatrix(rayDirection, rayDirection);
		localMinimumDistance *= rayDirection.length();
	    }
	    if (rayDirection.length() <= 0.0f)
		continue;
	    rayDirection.normalize();
	    BObolRtPickResult candidate;
	    if (source.cache->pickRay(candidate, rayOrigin, rayDirection,
		source.snapshot.localClipPlanes,
		source.snapshot.localClipPlaneCount,
		source.resources[workerIndex].get(), localMinimumDistance) &&
		candidate.hit) {
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
    const auto shadeLocalHit = [&](const SbLine &line,
	const BObolRtPickResult &hit,
	const RtSourceSnapshot &source) -> SbColor {
	SbVec3f normal = hit.normal;
	SbVec3f toEye = -line.getDirection();
	if (toEye.length() > 0.0f)
	    toEye.normalize();
	if (normal.length() > 0.0f) {
	    normal.normalize();
	    if (normal.dot(toEye) < 0.0f)
		normal = -normal;
	}
	const SbColor base = sourceColor(source, hit);
	if (!p->lightingEnabled)
	    return base;
	const RtMaterialModel material = materialModel(source, hit);
	SbVec3f shaded = base * (material.ambient * p->ambientIntensity);
	const auto addLight = [&](SbVec3f toLight, const SbColor &lightColor,
		float intensity) {
	    if (intensity <= 0.0f || toLight.length() <= 0.0f)
		return;
	    toLight.normalize();
	    const float diffuse = normal.length() > 0.0f ?
		std::max(0.0f, normal.dot(toLight)) : 1.0f;
	    if (diffuse <= 0.0f)
		return;
	    SbVec3f halfVector = toLight + toEye;
	    float specularCosine = diffuse;
	    if (halfVector.length() > 0.0f && normal.length() > 0.0f) {
		halfVector.normalize();
		specularCosine = std::max(0.0f, normal.dot(halfVector));
	    }
	    const float diffuseIllumination = material.diffuse * diffuse *
		intensity;
	    const float specular = material.specular *
		std::pow(specularCosine, material.shine) * intensity;
	    const SbColor litBase(base[0] * lightColor[0],
		base[1] * lightColor[1], base[2] * lightColor[2]);
	    shaded += litBase * diffuseIllumination + lightColor * specular;
	};

	/* Camera lights are already normalized world-space snapshots of the same
	 * key/fill/rim nodes traversed by GL and OSMesa. */
	for (const RtLightSnapshot &light : p->cameraLights)
	    addLight(-light.direction, light.color, light.intensity);
	for (const RtLightSnapshot &light : p->authoredLights) {
	    SbVec3f toLight;
	    float intensity = light.intensity;
	    if (light.kind == RtLightSnapshot::DIRECTIONAL) {
		toLight = -light.direction;
	    } else {
		toLight = light.position - hit.point;
		if (light.kind == RtLightSnapshot::SPOT) {
		    SbVec3f fromLight = hit.point - light.position;
		    if (fromLight.length() <= 0.0f ||
			light.direction.length() <= 0.0f)
			continue;
		    fromLight.normalize();
		    const float coneCosine = light.direction.dot(fromLight);
		    const float cutoffCosine = std::cos(light.cutOffAngle);
		    if (coneCosine < cutoffCosine)
			continue;
		    if (light.dropOffRate > 0.0f)
			intensity *= std::pow(std::max(0.0f, coneCosine),
			    light.dropOffRate * 128.0f);
		}
	    }
	    addLight(toLight, light.color, intensity);
	}
	return SbColor(shaded[0], shaded[1], shaded[2]);
    };

    /* Secondary rays stay deliberately bounded for interactive use.  Three
     * levels are enough to enter and leave a transmissive solid and then
     * shade geometry behind it, while limiting a two-way branch to fifteen
     * total traces per primary sample. */
    static const unsigned int maximumSecondaryDepth = 3u;
    static const float secondaryRayBias = 1.0e-4f;
    const auto rayBackground = [&settings](const SbVec3f &direction) {
	const float t = std::max(0.0f, std::min(1.0f,
	    0.5f * (direction[1] + 1.0f)));
	return settings.backgroundBottom * (1.0f - t) +
	    settings.backgroundTop * t;
    };
    std::function<SbColor(const SbLine &, unsigned int, unsigned int,
	const SbColor &)> traceRadiance;
    std::function<SbColor(const SbLine &, const BObolRtPickResult &,
	size_t, unsigned int, unsigned int)> shadeMaterial;

    shadeMaterial = [&](const SbLine &line,
	const BObolRtPickResult &hit, size_t sourceIndex,
	unsigned int workerIndex, unsigned int depth) -> SbColor {
	if (sourceIndex >= p->sources.size())
	    return SbColor(0.0f, 0.0f, 0.0f);
	const RtSourceSnapshot &source = p->sources[sourceIndex].snapshot;
	const SbColor local = shadeLocalHit(line, hit, source);
	const RtMaterialModel material = materialModel(source, hit);
	if (depth >= maximumSecondaryDepth ||
	    (material.reflectance <= 0.0f && material.transmission <= 0.0f))
	    return local;

	SbVec3f incoming = line.getDirection();
	SbVec3f geometricNormal = hit.normal;
	if (incoming.length() <= 0.0f || geometricNormal.length() <= 0.0f)
	    return local;
	incoming.normalize();
	geometricNormal.normalize();
	const bool entering = incoming.dot(geometricNormal) < 0.0f;
	SbVec3f orientedNormal = entering ? geometricNormal : -geometricNormal;
	const float cosine = std::max(0.0f,
	    std::min(1.0f, -incoming.dot(orientedNormal)));

	float reflectionWeight = material.reflectance;
	float transmissionWeight = material.transmission;
	SbVec3f transmissionDirection(0.0f, 0.0f, 0.0f);
	if (transmissionWeight > 0.0f) {
	    const float eta = entering ? 1.0f / material.refractiveIndex :
		material.refractiveIndex;
	    const float discriminant = 1.0f - eta * eta *
		(1.0f - cosine * cosine);
	    if (discriminant < 0.0f) {
		/* Total internal reflection retains the transmitted energy. */
		reflectionWeight += transmissionWeight;
		transmissionWeight = 0.0f;
	    } else {
		transmissionDirection = incoming * eta + orientedNormal *
		    (eta * cosine - std::sqrt(discriminant));
		if (transmissionDirection.length() > 0.0f)
		    transmissionDirection.normalize();
		else
		    transmissionWeight = 0.0f;
	    }
	}
	reflectionWeight = std::max(0.0f, std::min(1.0f, reflectionWeight));
	transmissionWeight = std::max(0.0f, std::min(1.0f,
	    transmissionWeight));
	const float secondaryWeight = reflectionWeight + transmissionWeight;
	if (secondaryWeight > 1.0f) {
	    reflectionWeight /= secondaryWeight;
	    transmissionWeight /= secondaryWeight;
	}
	const float localWeight = std::max(0.0f,
	    1.0f - reflectionWeight - transmissionWeight);
	SbColor result = local * localWeight;

	if (reflectionWeight > 0.0f && !shouldCancel(cancelled)) {
	    SbVec3f reflectionDirection = incoming - orientedNormal *
		(2.0f * incoming.dot(orientedNormal));
	    if (reflectionDirection.length() > 0.0f) {
		reflectionDirection.normalize();
		SbLine reflectionRay;
		reflectionRay.setPosDir(hit.point +
		    reflectionDirection * secondaryRayBias,
		    reflectionDirection);
		result += traceRadiance(reflectionRay, workerIndex, depth + 1u,
		    rayBackground(reflectionDirection)) * reflectionWeight;
	    }
	}
	if (transmissionWeight > 0.0f && !shouldCancel(cancelled)) {
	    SbLine transmissionRay;
	    transmissionRay.setPosDir(hit.point +
		transmissionDirection * secondaryRayBias,
		transmissionDirection);
	    result += traceRadiance(transmissionRay, workerIndex, depth + 1u,
		rayBackground(transmissionDirection)) * transmissionWeight;
	}
	return result;
    };

    traceRadiance = [&](const SbLine &line, unsigned int workerIndex,
	unsigned int depth, const SbColor &background) -> SbColor {
	if (shouldCancel(cancelled)) {
	    stop.store(true, std::memory_order_release);
	    return background;
	}
	BObolRtPickResult hit;
	size_t sourceIndex = p->sources.size();
	if (!traceLine(line, workerIndex, hit, sourceIndex, NULL,
		depth > 0u ? secondaryRayBias : 0.0f))
	    return background;
	return shadeMaterial(line, hit, sourceIndex, workerIndex, depth);
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

		    BObolRtPickResult nearest;
		    size_t nearestSourceIndex = p->sources.size();
		    const SbBool hasNearest = traceLine(line, workerIndex, nearest,
			nearestSourceIndex, transparentCompositionRequired ?
			&transparentHits : NULL, 0.0f) ? TRUE : FALSE;
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
			    const SbColor lit = shadeMaterial(line, hit->pick,
				hit->sourceIndex, workerIndex, 0u);
			    shaded = lit * alpha + shaded * (1.0f - alpha);
			}
		    } else if (nearestSource) {
			shaded = shadeMaterial(line, nearest,
			    nearestSourceIndex, workerIndex, 0u);
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
			BObolRtPickResult centerHit;
			size_t centerSourceIndex = p->sources.size();
			if (traceLine(centerLine, workerIndex, centerHit,
				centerSourceIndex, NULL, 0.0f)) {
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
BObolRtRenderer::getPreparedSourceCount(void) const
{
    return p ? static_cast<unsigned int>(p->sources.size()) : 0u;
}

uint64_t
BObolRtRenderer::getGeometryRevision(void) const
{
    return p ? p->geometryRevision : 0u;
}

uint64_t
BObolRtRenderer::getPresentationRevision(void) const
{
    return p ? p->presentationRevision : 0u;
}

SbBool
BObolRtRenderer::getSourceIdentity(uint32_t identifier,
	BObolRtSourceIdentity &identity) const
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
