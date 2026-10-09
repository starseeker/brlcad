/* D A T A B A S E _ S O U R C E _ C O M P A C T _ A C C E S S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

/** @file database_source_compact_access.cpp
 *
 * Queries, sparse presentation mutations, edit extraction, and exact actions
 * over the retained compact occurrence registry.
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BDrawCache.h"
#include "BObol/BEvaluatedPoints.h"
#include "BObol/BExportAction.h"
#include "BObol/BLodMeshShape.h"
#include "BObol/BLodRealization.h"
#include "BObol/BLodService.h"
#include "BObol/BMaterialObject.h"
#include "BObol/BMeasureAction.h"
#include "BObol/BMeshLodCache.h"
#include "BObol/BMeshShape.h"
#include "BObol/BPickDetail.h"
#include "BObol/BSnapAction.h"
#include "BObol/BViewLod.h"
#include "BObol/BViewQuery.h"
#include "BObol/BVListShape.h"
#include "cad_assembly_private.h"
#include "cad_publication_private.h"
#include "compact_occurrence_registry_private.h"
#include "database_source_private.h"
#include "database_source_realization.h"
#include "performance_private.h"

#include "bg/line_layer.h"
#include "bg/pca.h"
#include "bg/trimesh.h"
#include "bg/vlist.h"
#include "bu/app.h"
#include "bu/color.h"
#include "bu/cv.h"
#include "bu/file.h"
#include "bu/hash.h"
#include "bu/list.h"
#include "bu/mapped_file.h"
#include "bu/parallel.h"
#include "bu/str.h"
#include "bu/datetime.h"
#include "bu/vls.h"
#include "nmg.h"
#include "raytrace.h"
#include "rt/func.h"
#include "rt/global.h"
#include "rt/db4.h"
#include "rt/nongeom.h"
#include "rt/db_fullpath.h"
#include "rt/eval_wireframe.h"
#include "rt/primitives/annot.h"
#include "rt/tree.h"
#include "rt/vlist.h"
#include "rt/view.h"
#include "wdb.h"

#include <Inventor/SbName.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/nodes/SoGroup.h>

#include <Inventor/nodes/SoMatrixTransform.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Obol/cad/CadProjectedProxy.h>
#include <Inventor/sensors/SoFieldSensor.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <inttypes.h>
#include <limits.h>
#include <limits>
#include <map>
#include <math.h>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

SbBool
SoBRLDatabaseSource::hasCompactInstanceIndex(void) const
{
    return (this->d->compactIndexActive && this->d->compactIndex &&
	    !this->d->compactIndex->entries.empty()) ? TRUE : FALSE;
}

SbBool
SoBRLDatabaseSource::isCompactOccurrenceRegistry(void) const
{
    return this->d->compactIndexActive && this->d->compactOccurrenceRegistry;
}

int
SoBRLDatabaseSource::getCompactInstanceCount(void) const
{
    if (!this->hasCompactInstanceIndex())
	return 0;
    return static_cast<int>(this->d->compactIndex->entries.size());
}

int
SoBRLDatabaseSource::getCompactSelectedInstanceCount(void) const
{
    if (!this->hasCompactInstanceIndex())
	return 0;
    return static_cast<int>(
	this->d->compactIndex->selectedInstances.size());
}

size_t
SoBRLDatabaseSource::getCompactExpectedInstanceCount(void) const
{
    const size_t current = this->d->compactIndex ?
	this->d->compactIndex->entries.size() : 0;
    if (this->d->compactExpectedInstanceCountCertified)
	return this->d->compactExpectedInstanceCount;
    return std::max(current, this->d->compactExpectedInstanceCount);
}

SbBool
SoBRLDatabaseSource::hasCompleteCompactInstancePopulation(void) const
{
    if (!this->d->compactExpectedInstanceCountCertified ||
	!this->d->compactIndex)
	return FALSE;
    const size_t current = this->d->compactIndex->entries.size();
    const size_t overviews = this->d->compactIndex->overviewCount;
    if (overviews > current)
	return FALSE;
    return current - overviews ==
	this->d->compactExpectedInstanceCount ? TRUE : FALSE;
}

SbBool
SoBRLDatabaseSource::getCompactSourceProfile(
    BObolCompactSourceProfile &profile) const
{
    profile = this->d->compactSourceProfile;
    return profile.valid;
}

int
SoBRLDatabaseSource::getCompactPartCount(void) const
{
    if (!this->hasCompactInstanceIndex())
	return 0;
    return static_cast<int>(this->d->compactIndex->parts.size());
}

SbBool
SoBRLDatabaseSource::getCompactInstanceHandle(
    int index, BObolCompactInstanceHandle &handle) const
{
    handle = BObolCompactInstanceHandle();
    if (!this->d->compactIndex || index < 0 ||
	static_cast<size_t>(index) >= this->d->compactIndex->entries.size())
	return FALSE;

    const BObolCompactInstanceEntry &entry =
	this->d->compactIndex->entries[static_cast<size_t>(index)];
    handle.sourceNodeId = this->d->compactHandleSourceId;
    handle.instanceWord0 = entry.instance.w0;
    handle.instanceWord1 = entry.instance.w1;
    return handle.isValid();
}

SbBool
SoBRLDatabaseSource::getCompactOccurrence(
    int index, BObolCompactOccurrence &occurrence) const
{
    occurrence = BObolCompactOccurrence();
    if (!this->d->compactIndex || index < 0 ||
	static_cast<size_t>(index) >= this->d->compactIndex->entries.size())
	return FALSE;

    const BObolCompactInstanceEntry &entry =
	this->d->compactIndex->entries[static_cast<size_t>(index)];
    occurrence.geometry = entry.geometry;
    occurrence.summary = entry.shapeSummary;
    occurrence.geometryTransform = entry.geometryTransform;
    occurrence.localTransform = entry.placementTransform;
    occurrence.viewDependentCsgGeometry = entry.viewDependentCsgGeometry;
    occurrence.lodBacked = entry.lodBacked;
    occurrence.sourceMeshRequestValid = entry.sourceMeshRequestValid;
    occurrence.sourceMeshRequest = entry.sourceMeshRequest;
    occurrence.occurrenceIndex = entry.occurrenceIndex;
    occurrence.booleanOperation = entry.booleanOperation;
    return occurrence.geometry ? TRUE : FALSE;
}

namespace {

struct compact_rectangle_clip_point {
    double v[4];
};

static compact_rectangle_clip_point
compact_rectangle_transform(const SbMatrix &matrix, const SbVec3f &point)
{
    const float *m = matrix[0];
    compact_rectangle_clip_point result;
    for (int column = 0; column < 4; ++column) {
	result.v[column] = static_cast<double>(point[0]) * m[column] +
	    static_cast<double>(point[1]) * m[4 + column] +
	    static_cast<double>(point[2]) * m[8 + column] + m[12 + column];
    }
    return result;
}

static double
compact_rectangle_plane_value(const compact_rectangle_clip_point &point,
	int plane)
{
    switch (plane) {
	case 0: return point.v[3] + point.v[0];
	case 1: return point.v[3] - point.v[0];
	case 2: return point.v[3] + point.v[1];
	case 3: return point.v[3] - point.v[1];
	case 4: return point.v[3] + point.v[2];
	default: return point.v[3] - point.v[2];
    }
}

static bool
compact_rectangle_overlaps(const SbBox3f &localBounds,
	const SbMatrix &localToWorld, const SbMatrix &viewProjection,
	float minimumX, float minimumY, float maximumX, float maximumY,
	SbVec3f &worldPoint, float &distance)
{
    if (localBounds.isEmpty())
	return false;

    SbMatrix localToClip = localToWorld;
    localToClip.multRight(viewProjection);
    const SbVec3f bmin = localBounds.getMin();
    const SbVec3f bmax = localBounds.getMax();
    bool allOutside[6] = {true, true, true, true, true, true};
    double projectedMinimumX = 0.0;
    double projectedMinimumY = 0.0;
    double projectedMaximumX = 0.0;
    double projectedMaximumY = 0.0;
    double nearestDepth = 0.0;
    bool projected = false;
    for (int z = 0; z < 2; ++z) {
	for (int y = 0; y < 2; ++y) {
	    for (int x = 0; x < 2; ++x) {
		const SbVec3f corner(x ? bmax[0] : bmin[0],
		    y ? bmax[1] : bmin[1], z ? bmax[2] : bmin[2]);
		const compact_rectangle_clip_point clip =
		    compact_rectangle_transform(localToClip, corner);
		for (int plane = 0; plane < 6; ++plane)
		    allOutside[plane] = allOutside[plane] &&
			compact_rectangle_plane_value(clip, plane) < 0.0;
		if (!std::isfinite(clip.v[3]) ||
		    std::fabs(clip.v[3]) < 1.0e-20)
		    continue;
		const double px = clip.v[0] / clip.v[3];
		const double py = clip.v[1] / clip.v[3];
		const double pz = clip.v[2] / clip.v[3];
		if (!std::isfinite(px) || !std::isfinite(py) ||
		    !std::isfinite(pz))
		    continue;
		if (!projected) {
		    projectedMinimumX = projectedMaximumX = px;
		    projectedMinimumY = projectedMaximumY = py;
		    nearestDepth = pz;
		    projected = true;
		} else {
		    projectedMinimumX = std::min(projectedMinimumX, px);
		    projectedMinimumY = std::min(projectedMinimumY, py);
		    projectedMaximumX = std::max(projectedMaximumX, px);
		    projectedMaximumY = std::max(projectedMaximumY, py);
		    nearestDepth = std::min(nearestDepth, pz);
		}
	    }
	}
    }
    /* Object selection spans the model depth represented by the view, as the
     * legacy bview selection prism did.  Near/far planes are renderer
     * precision controls and can lag a direct camera edit; letting them erase
     * otherwise valid screen-space candidates makes picking depend on an
     * unrelated clip synchronization detail. */
    for (int plane = 0; plane < 4; ++plane)
	if (allOutside[plane])
	    return false;
    if (!projected || projectedMaximumX < minimumX ||
	projectedMinimumX > maximumX || projectedMaximumY < minimumY ||
	projectedMinimumY > maximumY)
	return false;

    localToWorld.multVecMatrix(localBounds.getCenter(), worldPoint);
    distance = static_cast<float>(nearestDepth);
    return true;
}

}

int
SoBRLDatabaseSource::queryCompactRectangle(const SbMatrix &parentToWorld,
	const SbMatrix &viewProjection, const SbVec2s &viewportSize,
	float minimumX, float minimumY, float maximumX, float maximumY,
	std::vector<BObolViewPickRecord> &records) const
{
    if (!this->d->compactIndex || this->d->compactIndex->entries.empty())
	return -1;
    if (!this->visible.getValue())
	return 0;

    const size_t initialCount = records.size();
    for (const BObolCompactInstanceEntry &entry :
	 this->d->compactIndex->entries) {
	if (!entry.visible || !entry.selectable || !entry.geometry)
	    continue;
	SbBox3f localBounds = compact_part_geometry_bounds(entry.geometry);
	SbMatrix localToWorld = entry.localToSource;
	localToWorld.multRight(parentToWorld);
	if (entry.geometry->displayPlane) {
	    if (!Obol::cadDisplayPlaneTransform(*entry.geometry->displayPlane,
		    localToWorld, viewProjection, viewportSize, localToWorld))
		continue;
	    localBounds = Obol::cadPartGeometryBounds(*entry.geometry);
	}
	SbVec3f worldPoint;
	float distance = FLT_MAX;
	if (!compact_rectangle_overlaps(localBounds, localToWorld,
		viewProjection, minimumX, minimumY, maximumX, maximumY,
		worldPoint, distance))
	    continue;

	BObolViewPickRecord record;
	record.point = worldPoint;
	record.distance = distance;
	record.detail.setPath(entry.semantic.path);
	record.detail.setSourceInstanceKey(compact_instance_identity(entry));
	record.detail.setSourceName(entry.semantic.sourceName);
	record.detail.setSourceType(entry.semantic.sourceType);
	record.detail.setSourceId(entry.semantic.sourceId);
	record.detail.setRegionId(entry.semantic.regionId);
	record.detail.setAirCode(entry.semantic.airCode);
	record.detail.setMaterialId(entry.semantic.materialId);
	record.detail.setLos(entry.semantic.los);
	record.detail.setMaterialColor(entry.semantic.materialColorValid,
	    entry.semantic.materialColor);
	record.detail.setMaterialShader(entry.semantic.materialShader);
	record.detail.setEditIntent(entry.semantic.editIntentId,
	    entry.semantic.editIntentRole);
	record.detail.setModelPoint(worldPoint);
	record.detail.setPrimitive(SoBRLPickDetail::UNKNOWN, -1);
	records.push_back(record);
    }
    return static_cast<int>(records.size() - initialCount);
}

int
SoBRLDatabaseSource::querySourceRectangle(const SbMatrix &parentToWorld,
	const SbMatrix &viewProjection,
	float minimumX, float minimumY, float maximumX, float maximumY,
	std::vector<BObolViewPickRecord> &records) const
{
    /* An occurrence registry has more precise identities and bounds.  This
     * fallback must not add the source root alongside those occurrences. */
    if (this->d->compactIndex && !this->d->compactIndex->entries.empty())
	return -1;
    if (!this->visible.getValue())
	return 0;
    SbBox3f localBounds;
    if (!this->getSourceBounds(localBounds))
	return -1;

    SbMatrix localToWorld;
    localToWorld.makeIdentity();
    if (this->drawMatrixValid.getValue())
	localToWorld = this->drawMatrix.getValue();
    localToWorld.multRight(parentToWorld);
    SbVec3f worldPoint;
    float distance = FLT_MAX;
    if (!compact_rectangle_overlaps(localBounds, localToWorld,
	    viewProjection, minimumX, minimumY, maximumX, maximumY,
	    worldPoint, distance))
	return 0;

    BObolViewPickRecord record;
    record.point = worldPoint;
    record.distance = distance;
    record.detail.setPath(this->path.getValue());
    record.detail.setSourceInstanceKey(this->instanceKey.getValue());
    record.detail.setSourceName(this->displayName.getValue().getLength() > 0 ?
	this->displayName.getValue() : this->path.getValue());
    record.detail.setSourceId(this->sourceRevision.getValue());
    record.detail.setRegionId(this->databaseRegionId.getValue());
    record.detail.setAirCode(this->databaseAirCode.getValue());
    record.detail.setMaterialId(this->databaseMaterialId.getValue());
    record.detail.setLos(this->databaseLos.getValue());
    record.detail.setMaterialColor(
	this->databaseMaterialColorValid.getValue(),
	this->databaseMaterialColor.getValue());
    record.detail.setMaterialShader(this->databaseMaterialShader.getValue());
    record.detail.setModelPoint(worldPoint);
    record.detail.setPrimitive(SoBRLPickDetail::UNKNOWN, -1);
    records.push_back(record);
    return 1;
}

SbBool
SoBRLDatabaseSource::copyCompactWireGeometry(
    std::vector<SbVec3f> &points, std::vector<int32_t> &commands) const
{
    points.clear();
    commands.clear();
    if (!this->d->compactIndex)
	return FALSE;

    for (const BObolCompactInstanceEntry &entry :
	 this->d->compactIndex->entries) {
	if (!entry.visible || !entry.geometry || !entry.geometry->wire)
	    continue;

	const Obol::WireRep &wire = *entry.geometry->wire;
	const auto appendPoint = [&entry, &points, &commands](
	    const SbVec3f &point, int32_t command) {
	    SbVec3f transformed;
	    entry.localToSource.multVecMatrix(point, transformed);
	    points.push_back(transformed);
	    commands.push_back(command);
	};

	for (size_t i = 1; i < wire.segmentPoints.size(); i += 2) {
	    appendPoint(wire.segmentPoints[i - 1], 0);
	    appendPoint(wire.segmentPoints[i], 1);
	}
	for (const Obol::WirePolyline &polyline : wire.polylines) {
	    for (size_t i = 0; i < polyline.points.size(); i++)
		appendPoint(polyline.points[i], i == 0 ? 0 : 1);
	}
    }

    return points.empty() ? FALSE : TRUE;
}

SbBool
SoBRLDatabaseSource::getCompactStructuralPresentationCounts(
    int index, BObolLodCounts &counts) const
{
    counts.clear();
    if (!this->d->compactIndex || index < 0 ||
	static_cast<size_t>(index) >= this->d->compactIndex->entries.size())
	return FALSE;
    const BObolCompactInstanceEntry &entry =
	this->d->compactIndex->entries[static_cast<size_t>(index)];
    if (!entry.geometry || !entry.geometry->structuralProxy)
	return FALSE;
    counts = bobol_cad_geometry_counts(*entry.geometry);
    return TRUE;
}

SbBool
SoBRLDatabaseSource::copyCompactWireSegments(
    std::vector<BObolCompactWireSegment> &segments) const
{
    segments.clear();
    if (!this->d->compactIndex)
	return FALSE;

    for (const BObolCompactInstanceEntry &entry :
	 this->d->compactIndex->entries) {
	if (!entry.visible || !entry.geometry || !entry.geometry->wire)
	    continue;

	const Obol::WireRep &wire = *entry.geometry->wire;
	size_t segment_index = 0;
	const auto append_segment = [&entry, &wire, &segments, &segment_index](
	    const SbVec3f &local_start, const SbVec3f &local_end) {
	    BObolCompactWireSegment segment;
	    entry.localToSource.multVecMatrix(local_start, segment.start);
	    entry.localToSource.multVecMatrix(local_end, segment.end);

	    segment.color.setValue(entry.style.color[0], entry.style.color[1],
		entry.style.color[2]);
	    float alpha = entry.style.color[3];
	    segment.lineWidth = (std::max)(1.0f, entry.style.lineWidth);
	    segment.linePattern = entry.style.linePattern;
	    segment.linePatternFactor = entry.style.linePatternFactor;

	    if (!wire.styleRuns.empty()) {
		const Obol::WireStyle authored =
		    wire.styleAtSegment(segment_index);
		segment.lineWidth = Obol::cadWirePixelWidth(
		    segment.lineWidth, authored.widthScale);
		if (authored.patternValid) {
		    segment.linePattern = authored.linePattern;
		    segment.linePatternFactor = authored.linePatternFactor;
		}
		if (authored.colorValid && entry.style.useGeometryColor) {
		    segment.color.setValue(authored.color[0], authored.color[1],
			authored.color[2]);
		    alpha *= authored.color[3];
		}
	    }

	    alpha = (std::max)(0.0f, (std::min)(1.0f, alpha));
	    segment.transparency = 1.0f - alpha;
	    segments.push_back(segment);
	    segment_index++;
	};

	for (size_t i = 0; i + 1 < wire.segmentPoints.size(); i += 2)
	    append_segment(wire.segmentPoints[i], wire.segmentPoints[i + 1]);
	for (const Obol::WirePolyline &polyline : wire.polylines) {
	    for (size_t i = 1; i < polyline.points.size(); i++)
		append_segment(polyline.points[i - 1], polyline.points[i]);
	}
    }

    return segments.empty() ? FALSE : TRUE;
}

const BObolCompactInstanceEntry *
SoBRLDatabaseSource::findCompactInstanceEntry(
	const BObolCompactInstanceHandle &handle)
    const
{
    if (!this->d->compactIndex || !handle.isValid() ||
	handle.sourceNodeId != this->d->compactHandleSourceId)
	return NULL;
    Obol::InstanceId instance;
    instance.w0 = handle.instanceWord0;
    instance.w1 = handle.instanceWord1;
    const auto found = this->d->compactIndex->entryIndex.find(instance);
    if (found == this->d->compactIndex->entryIndex.end() ||
	found->second >= this->d->compactIndex->entries.size())
	return NULL;
    return &this->d->compactIndex->entries[found->second];
}

SbBool
SoBRLDatabaseSource::isCompactInstanceHandleValid(
    const BObolCompactInstanceHandle &handle) const
{
    return this->findCompactInstanceEntry(handle) ? TRUE : FALSE;
}

SbBool
SoBRLDatabaseSource::hasCompactInstanceKey(const char *occurrenceKey) const
{
    if (!this->d->compactIndex || !occurrenceKey || !occurrenceKey[0])
	return FALSE;

    const auto found =
	this->d->compactIndex->entryIndexByKey.find(occurrenceKey);
    if (found == this->d->compactIndex->entryIndexByKey.end() ||
	found->second >= this->d->compactIndex->entries.size())
	return FALSE;

    /* Keep the explicit string comparison as a defensive consistency check:
     * modern compact occurrence keys are labels derived from the instance ID,
     * not strings whose hash is the instance ID. */
    return bu_strcmp(compact_instance_identity(
	this->d->compactIndex->entries[found->second]).getString(),
	occurrenceKey) == 0 ? TRUE : FALSE;
}

SbBool
SoBRLDatabaseSource::getCompactInstanceIndex(
    const char *occurrenceKey, size_t &entryIndex) const
{
    entryIndex = 0;
    if (!this->d->compactIndex || !occurrenceKey || !occurrenceKey[0])
	return FALSE;
    const auto found =
	this->d->compactIndex->entryIndexByKey.find(occurrenceKey);
    if (found == this->d->compactIndex->entryIndexByKey.end() ||
	found->second >= this->d->compactIndex->entries.size())
	return FALSE;
    entryIndex = found->second;
    return TRUE;
}

SbBool
SoBRLDatabaseSource::getCompactInstanceIndex(
    const Obol::InstanceId &instance, size_t &entryIndex) const
{
    entryIndex = 0;
    if (!this->d->compactIndex || !instance.isValid())
	return FALSE;
    const auto found = this->d->compactIndex->entryIndex.find(instance);
    if (found == this->d->compactIndex->entryIndex.end() ||
	found->second >= this->d->compactIndex->entries.size())
	return FALSE;
    entryIndex = found->second;
    return TRUE;
}

uint64_t
SoBRLDatabaseSource::getCompactSourceRoutingId(void) const
{
    return this->d->routingId;
}

uint64_t
SoBRLDatabaseSource::getCompactPopulationEpoch(void) const
{
    return this->d->compactPopulationEpoch;
}

SbBool
SoBRLDatabaseSource::getCompactInstanceSummary(
    const BObolCompactInstanceHandle &handle,
    BObolCompactInstanceSummary &summary) const
{
    summary = BObolCompactInstanceSummary();
    const BObolCompactInstanceEntry *entry =
	this->findCompactInstanceEntry(handle);
    if (!entry)
	return FALSE;

    summary.valid = TRUE;
    summary.handle = handle;
    summary.path = entry->semantic.path;
    summary.sourceName = entry->semantic.sourceName;
    summary.sourceInstanceKey = compact_instance_identity(*entry);
    summary.geometryKind = entry->shapeSummary.geometryKind;
    if (entry->sourceMeshRequestValid) {
	summary.meshAssetPath = entry->sourceMeshRequest.meshAssetPath;
	summary.meshAssetName = entry->sourceMeshRequest.meshAssetName;
	summary.meshAssetBounds = entry->sourceMeshRequest.meshAssetBounds;
	summary.sourceContentHash =
	    entry->sourceMeshRequest.meshAssetContentHash;
	summary.sourceFaceCount = entry->sourceMeshRequest.faceCount;
	summary.sourcePointCount = entry->sourceMeshRequest.pointCount;
    }
    summary.localToSource = entry->localToSource;
    summary.geometryIdentity = entry->part.w0 ^
	(entry->part.w1 + 0x9e3779b97f4a7c15ULL +
	 (entry->part.w0 << 6) + (entry->part.w0 >> 2));
    summary.geometryRevision = entry->geometryRevision;
    summary.appearanceRevision = entry->appearanceRevision;
    summary.placementRevision = entry->placementRevision;
    summary.visibilityRevision = entry->visibilityRevision;
    summary.selectionRevision = entry->selectionRevision;
    summary.occurrenceIndex = entry->occurrenceIndex;
    summary.booleanOperation = entry->booleanOperation;
    summary.regionId = entry->semantic.regionId;
    summary.airCode = entry->semantic.airCode;
    summary.materialId = entry->semantic.materialId;
    summary.los = entry->semantic.los;
    summary.materialColorValid = entry->semantic.materialColorValid;
    summary.materialColor = entry->semantic.materialColor;
    summary.materialShader = entry->semantic.materialShader;
    summary.appearanceColorValid = entry->style.hasColorOverride ? TRUE : FALSE;
    summary.appearanceColor = SbColor(entry->style.color[0],
	entry->style.color[1], entry->style.color[2]);
    summary.lineStyle = entry->style.linePattern == 0xffffu ? 0 : 1;
    summary.lineWidth = entry->style.lineWidth > 0.0f ?
	static_cast<int>(entry->style.lineWidth + 0.5f) : 0;
    summary.transparency = 1.0f - entry->style.color[3];
    if (summary.transparency < 0.0f)
	summary.transparency = 0.0f;
    else if (summary.transparency > 1.0f)
	summary.transparency = 1.0f;
    summary.wireGeometry = entry->wireGeometry;
    summary.pointGeometry = entry->pointGeometry;
    summary.meshGeometry = entry->meshGeometry;
    summary.lodBacked = entry->lodBacked;
    summary.sourceMeshRequestValid = entry->sourceMeshRequestValid;
    summary.localBounds = compact_part_geometry_bounds(entry->geometry);
    summary.visible = entry->visible;
    summary.selectable = entry->selectable;
    summary.selected = entry->selected;
    summary.highlighted = entry->highlighted;
    return TRUE;
}

SbBool
SoBRLDatabaseSource::getCompactLodInstanceSummary(
    int index, BObolCompactLodInstanceSummary &summary) const
{
    summary = BObolCompactLodInstanceSummary();
    if (!this->d->compactIndex || index < 0 ||
	static_cast<size_t>(index) >= this->d->compactIndex->entries.size())
	return FALSE;

    const BObolCompactInstanceEntry &entry =
	this->d->compactIndex->entries[static_cast<size_t>(index)];
    summary.valid = TRUE;
    summary.path = entry.semantic.path;
    summary.sourceName = entry.semantic.sourceName;
    summary.sourceInstanceKey = compact_instance_identity(entry);
    if (entry.sourceMeshRequestValid) {
	summary.meshAssetPath = entry.sourceMeshRequest.meshAssetPath;
	summary.meshAssetName = entry.sourceMeshRequest.meshAssetName;
	summary.meshAssetBounds = entry.sourceMeshRequest.meshAssetBounds;
	summary.sourceContentHash =
	    entry.sourceMeshRequest.meshAssetContentHash;
	summary.sourceFaceCount = entry.sourceMeshRequest.faceCount;
	summary.sourcePointCount = entry.sourceMeshRequest.pointCount;
	summary.brepSource = BU_STR_EQUAL(
	    entry.sourceMeshRequest.sourceType.getString(), "brep") ?
	    TRUE : FALSE;
	summary.meshAssetTessellationAbsTol =
	    entry.sourceMeshRequest.meshAssetTessellationAbsTol;
	summary.meshAssetTessellationRelTol =
	    entry.sourceMeshRequest.meshAssetTessellationRelTol;
	summary.meshAssetTessellationNormTol =
	    entry.sourceMeshRequest.meshAssetTessellationNormTol;
    }
    summary.localToSource = entry.sourceMeshRequestValid ?
	compact_mesh_asset_matrix(this, entry) : entry.localToSource;
    summary.localBounds = entry.sourceMeshRequestValid &&
	!entry.sourceMeshRequest.meshAssetBounds.isEmpty() ?
	entry.sourceMeshRequest.meshAssetBounds :
	compact_part_geometry_bounds(entry.geometry);
    summary.presentationLocalToSource = entry.localToSource;
    summary.presentationLocalBounds =
	compact_part_geometry_bounds(entry.geometry);
    summary.presentationCornersValid =
	entry.geometry && Obol::cadPartGeometryProxyCorners(
	    *entry.geometry, summary.presentationCorners.data()) ? TRUE : FALSE;
    summary.meshGeometry = entry.meshGeometry;
    summary.lodBacked = entry.lodBacked;
    summary.sourceMeshRequestValid = entry.sourceMeshRequestValid;
    summary.visible = entry.visible;
    summary.selected = entry.selected;
    summary.highlighted = entry.highlighted;
    return TRUE;
}

SbBool
SoBRLDatabaseSource::getCompactResidentProgressiveSummary(
    int index, BObolCompactResidentProgressiveSummary &summary) const
{
    summary = BObolCompactResidentProgressiveSummary();
    if (!this->d->compactIndex || index < 0 ||
	static_cast<size_t>(index) >= this->d->compactIndex->entries.size())
	return FALSE;

    const BObolCompactInstanceEntry &entry =
	this->d->compactIndex->entries[static_cast<size_t>(index)];
    if (!entry.geometry || !entry.geometry->wire)
	return FALSE;
    const Obol::WireRep &wire = *entry.geometry->wire;
    if (!wire.isProgressive() || !wire.hasProgressiveErrorBounds() ||
	wire.progressiveCuts.size() > summary.primitiveCounts.size())
	return FALSE;

    summary.valid = TRUE;
    summary.wire = TRUE;
    summary.minimumCut = wire.progressiveMinimumCut;
    summary.residentCut = wire.progressiveResidentCut;
    for (size_t cut = wire.progressiveMinimumCut;
	 cut <= wire.progressiveResidentCut; ++cut) {
	summary.primitiveCounts[cut] = wire.segmentCountAtCut(
	    static_cast<uint8_t>(cut));
	summary.normalizedErrors[cut] = wire.normalizedErrorAtCut(
	    static_cast<uint8_t>(cut));
    }
    return TRUE;
}

SbBool
SoBRLDatabaseSource::getCompactLodProviderSummary(
    int index, BObolCompactLodProviderSummary &summary) const
{
    summary = BObolCompactLodProviderSummary();
    if (!this->d->compactIndex || index < 0 ||
	static_cast<size_t>(index) >= this->d->compactIndex->entries.size())
	return FALSE;
    const BObolCompactInstanceEntry &entry =
	this->d->compactIndex->entries[static_cast<size_t>(index)];
    if (!entry.sourceMeshRequestValid)
	return FALSE;
    if (this->d->compactStagedSourceStream) {
	summary.stagedSource =
	    this->d->compactStagedSourceStream->claimStagedSource(
		entry.sourceMeshRequest.stagedSource);
	if (!this->d->compactStagedSourceStream->stagedSourceByteCount())
	    this->d->compactStagedSourceStream.reset();
    }
    if (!summary.stagedSource)
	summary.stagedSource = entry.sourceMeshRequest.stagedSource.lock();
    summary.lodAvailable = entry.sourceMeshRequest.lodAvailable ?
	TRUE : FALSE;
    summary.lodActiveCut = entry.sourceMeshRequest.lodActiveCut;
    summary.lodFaceCount = entry.sourceMeshRequest.lodFaceCount;
    summary.lodPointCount = entry.sourceMeshRequest.lodPointCount;
    summary.lodHasNormals = entry.sourceMeshRequest.lodHasNormals ?
	TRUE : FALSE;
    return TRUE;
}

SbBool
SoBRLDatabaseSource::getCompactLodPlanningSummary(
    int index, BObolCompactLodPlanningSummary &summary) const
{
    summary = BObolCompactLodPlanningSummary();
    if (!this->d->compactIndex || index < 0 ||
	static_cast<size_t>(index) >= this->d->compactIndex->entries.size())
	return FALSE;

    const BObolCompactInstanceEntry &entry =
	this->d->compactIndex->entries[static_cast<size_t>(index)];
    summary.valid = TRUE;
    summary.sourceInstanceKey = compact_instance_identity(entry);
    summary.geometryRevision = entry.geometryRevision;
    summary.placementRevision = entry.placementRevision;
    if (entry.sourceMeshRequestValid) {
	summary.sourceContentHash =
	    entry.sourceMeshRequest.meshAssetContentHash;
	summary.sourceFaceCount = entry.sourceMeshRequest.faceCount;
	summary.sourcePointCount = entry.sourceMeshRequest.pointCount;
	summary.botSource = BU_STR_EQUAL(
	    entry.sourceMeshRequest.sourceType.getString(), "bot") ?
	    TRUE : FALSE;
	summary.brepSource = BU_STR_EQUAL(
	    entry.sourceMeshRequest.sourceType.getString(), "brep") ?
	    TRUE : FALSE;
	summary.meshAssetTessellationAbsTol =
	    entry.sourceMeshRequest.meshAssetTessellationAbsTol;
	summary.meshAssetTessellationRelTol =
	    entry.sourceMeshRequest.meshAssetTessellationRelTol;
	summary.meshAssetTessellationNormTol =
	    entry.sourceMeshRequest.meshAssetTessellationNormTol;
    }
    summary.localToSource = entry.sourceMeshRequestValid ?
	compact_mesh_asset_matrix(this, entry) : entry.localToSource;
    summary.localBounds = entry.sourceMeshRequestValid &&
	!entry.sourceMeshRequest.meshAssetBounds.isEmpty() ?
	entry.sourceMeshRequest.meshAssetBounds :
	compact_part_geometry_bounds(entry.geometry);
    summary.presentationLocalToSource = entry.localToSource;
    summary.presentationLocalBounds =
	compact_part_geometry_bounds(entry.geometry);
    summary.presentationCornersValid =
	entry.geometry && Obol::cadPartGeometryProxyCorners(
	    *entry.geometry, summary.presentationCorners.data()) ? TRUE : FALSE;
    summary.meshGeometry = entry.meshGeometry;
    summary.lodBacked = entry.lodBacked;
    summary.sourceMeshRequestValid = entry.sourceMeshRequestValid;
    summary.residentProgressiveGeometry =
	bobol_compact_geometry_is_resident_progressive(entry.geometry) ?
	    TRUE : FALSE;
    summary.visible = entry.visible;
    summary.selected = entry.selected;
    summary.highlighted = entry.highlighted;
    return TRUE;
}

SbBool
SoBRLDatabaseSource::getCompactLodPlanningSummaryForKey(
    const char *occurrenceKey, BObolCompactLodPlanningSummary &summary) const
{
    summary = BObolCompactLodPlanningSummary();
    if (!this->d->compactIndex || !occurrenceKey || !occurrenceKey[0])
	return FALSE;
    const auto found =
	this->d->compactIndex->entryIndexByKey.find(occurrenceKey);
    if (found == this->d->compactIndex->entryIndexByKey.end() ||
	found->second >= this->d->compactIndex->entries.size())
	return FALSE;
    return this->getCompactLodPlanningSummary(
	static_cast<int>(found->second), summary);
}


SbBool
SoBRLDatabaseSource::copyCompactInstanceEditGeometry(
    const BObolCompactInstanceHandle &handle,
    std::vector<SbVec3f> &points,
    std::vector<int32_t> &commands,
    BObolCompactInstanceSummary &summary) const
{
    points.clear();
    commands.clear();
    summary = BObolCompactInstanceSummary();

    const BObolCompactInstanceEntry *entry =
	this->findCompactInstanceEntry(handle);
    if (!entry || !entry->geometry ||
	!this->getCompactInstanceSummary(handle, summary))
	return FALSE;

    const auto appendPoint = [&entry, &points, &commands](
	const SbVec3f &point, int32_t command) {
	SbVec3f transformed;
	entry->localToSource.multVecMatrix(point, transformed);
	points.push_back(transformed);
	commands.push_back(command);
    };

    if (entry->geometry->wire) {
	const Obol::WireRep &wire = *entry->geometry->wire;
	for (size_t i = 1; i < wire.segmentPoints.size(); i += 2) {
	    appendPoint(wire.segmentPoints[i - 1], 0);
	    appendPoint(wire.segmentPoints[i], 1);
	}
	for (const Obol::WirePolyline &polyline : wire.polylines) {
	    for (size_t i = 0; i < polyline.points.size(); i++)
		appendPoint(polyline.points[i], i == 0 ? 0 : 1);
	}
    }

    if (entry->geometry->points) {
	const Obol::PointRep &pointRep = *entry->geometry->points;
	for (const SbVec3f &point : pointRep.positions)
	    appendPoint(point, 2);
    }

    /* Mesh-only compact occurrences still need an editable visual.  Build a
     * transient triangle-edge preview without adding a persistent mesh shape
     * to the compact index. */
    if (points.empty() && entry->geometry->shaded) {
	const Obol::TriMesh &mesh = *entry->geometry->shaded;
	for (size_t i = 0; i + 2 < mesh.indices.size(); i += 3) {
	    const uint32_t a = mesh.indices[i];
	    const uint32_t b = mesh.indices[i + 1];
	    const uint32_t c = mesh.indices[i + 2];
	    if (a >= mesh.positions.size() || b >= mesh.positions.size() ||
		c >= mesh.positions.size())
		continue;
	    appendPoint(mesh.positions[a], 0);
	    appendPoint(mesh.positions[b], 1);
	    appendPoint(mesh.positions[c], 1);
	    appendPoint(mesh.positions[a], 1);
	}
    }

    if (points.empty()) {
	summary = BObolCompactInstanceSummary();
	return FALSE;
    }
    return TRUE;
}

template <typename Visitor>
static void
compact_visit_entries_for_path(const BObolCompactInstanceIndex *index,
	const char *queryPath, SbBool includeDescendants, Visitor visitor)
{
    if (!index)
	return;

    const char *query = database_source_skip_leading_slash(
	queryPath ? queryPath : "");
    if (!query[0]) {
	for (size_t entryIndex = 0; entryIndex < index->entries.size();
		entryIndex++)
	    visitor(entryIndex);
	return;
    }

    const bool leafQuery = !strchr(query, '/') && !strchr(query, '@');
    if (leafQuery) {
	auto leafEntries = index->entryIndicesByLeaf.find(query);
	if (leafEntries != index->entryIndicesByLeaf.end()) {
	    for (size_t entryIndex : leafEntries->second)
		visitor(entryIndex);
	}
	if (!includeDescendants)
	    return;
    }

    const std::string pathKey(query);
    const size_t prefixLength = pathKey.size();
    auto entryIt = index->entryIndexByOrderedPath.lower_bound(pathKey);
    for (; entryIt != index->entryIndexByOrderedPath.end(); ++entryIt) {
	const char *candidate = entryIt->first.c_str();
	if (bu_strncmp(candidate, pathKey.c_str(), prefixLength) != 0)
	    break;
	const char suffix = candidate[prefixLength];
	if (!includeDescendants && suffix != '\0')
	    break;
	if (includeDescendants && suffix != '\0' && suffix != '/' &&
	    suffix != '@')
	    continue;
	const size_t entryIndex = entryIt->second;
	if (entryIndex >= index->entries.size())
	    continue;
	if (leafQuery && database_source_leaf_component(
		index->entries[entryIndex].semantic.path) == pathKey)
	    continue;
	visitor(entryIndex);
	if (!includeDescendants)
	    continue;
    }
}

template <typename Visitor>
static void
compact_visit_entries_for_path_match(const BObolCompactInstanceIndex *index,
	const char *queryPath, BObolCompactPathMatch match, Visitor visitor)
{
    if (!index)
	return;

    const char *query = database_source_skip_leading_slash(
	queryPath ? queryPath : "");
    if (!query[0]) {
	for (size_t entryIndex = 0; entryIndex < index->entries.size();
		entryIndex++)
	    visitor(entryIndex);
	return;
    }

    if (match == BOBOL_COMPACT_PATH_OBJECT) {
	const std::string object = database_source_leaf_component(
	    SbString(query));
	auto leafEntries = index->entryIndicesByLeaf.find(object);
	if (leafEntries == index->entryIndicesByLeaf.end())
	    return;
	for (size_t entryIndex : leafEntries->second)
	    visitor(entryIndex);
	return;
    }

    const std::string pathKey(query);
    if (match == BOBOL_COMPACT_PATH_EXACT) {
	auto entry = index->entryIndexByOrderedPath.find(pathKey);
	if (entry != index->entryIndexByOrderedPath.end() &&
	    entry->second < index->entries.size())
	    visitor(entry->second);
	return;
    }

    const size_t prefixLength = pathKey.size();
    auto entry = index->entryIndexByOrderedPath.lower_bound(pathKey);
    for (; entry != index->entryIndexByOrderedPath.end(); ++entry) {
	const char *candidate = entry->first.c_str();
	if (bu_strncmp(candidate, pathKey.c_str(), prefixLength))
	    break;
	const char suffix = candidate[prefixLength];
	if (suffix != '\0' && suffix != '/' && suffix != '@')
	    continue;
	if (entry->second < index->entries.size())
	    visitor(entry->second);
    }
}

int
SoBRLDatabaseSource::getCompactInstanceCountForPath(const char *queryPath,
    SbBool includeDescendants) const
{
    if (!this->d->compactIndex)
	return 0;
    int count = 0;
    compact_visit_entries_for_path(this->d->compactIndex, queryPath,
	includeDescendants, [&count](size_t UNUSED(entryIndex)) {
	    count++;
	});
    return count;
}

SbBool
SoBRLDatabaseSource::getCompactInstanceForPath(const char *queryPath,
    SbBool includeDescendants, SbBool visibleOnly,
    BObolCompactInstanceHandle &handle,
    BObolCompactInstanceSummary &summary) const
{
    handle = BObolCompactInstanceHandle();
    summary = BObolCompactInstanceSummary();
    if (!this->d->compactIndex)
	return FALSE;

    size_t matchIndex = this->d->compactIndex->entries.size();
    compact_visit_entries_for_path(this->d->compactIndex, queryPath,
	includeDescendants, [this, visibleOnly, &matchIndex](size_t entryIndex) {
	    if (visibleOnly &&
		!this->d->compactIndex->entries[entryIndex].visible)
		return;
	    if (entryIndex < matchIndex)
		matchIndex = entryIndex;
	});
    if (matchIndex >= this->d->compactIndex->entries.size() ||
	matchIndex > static_cast<size_t>(INT_MAX) ||
	!this->getCompactInstanceHandle(static_cast<int>(matchIndex), handle) ||
	!this->getCompactInstanceSummary(handle, summary)) {
	handle = BObolCompactInstanceHandle();
	summary = BObolCompactInstanceSummary();
	return FALSE;
    }
    return TRUE;
}

SbBool
SoBRLDatabaseSource::getCompactInstanceBoundsForPath(const char *queryPath,
    SbBool includeDescendants, SbBox3f &bounds) const
{
    bounds.makeEmpty();
    if (!this->d->compactIndex || !this->visible.getValue())
	return FALSE;

    compact_visit_entries_for_path(this->d->compactIndex, queryPath,
	includeDescendants, [this, &bounds](size_t entryIndex) {
	const BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[entryIndex];
	if (!entry.visible)
	    return;
	const SbBox3f localBounds = compact_part_geometry_bounds(entry.geometry);
	if (!localBounds.isEmpty())
	    bounds.extendBy(database_source_transform_bounds(localBounds,
		entry.localToSource));
    });
    return bounds.isEmpty() ? FALSE : TRUE;
}

BObolCompactMetadataEdit::BObolCompactMetadataEdit(const SoBRLCadAssembly::InstanceSemantic &current,
    int region, int air, int material, int lineOfSight, SbBool valid, const SbColor &rgb, const char *shader) :
    regionId(region), airCode(air), materialId(material), los(lineOfSight), colorValid(valid), color(rgb),
    shaderChanged(bu_strcmp(current.materialShader.getString(), shader) != 0)
{
    if (this->shaderChanged) {
	this->semanticShader = shader;
	this->summaryShader = shader;
    }
}

void
BObolCompactMetadataEdit::commit(BObolCompactInstanceEntry &entry)
{
    auto &semantic = entry.semantic;
    auto &summary = entry.shapeSummary;
    semantic.regionId = summary.regionId = this->regionId;
    semantic.airCode = summary.airCode = this->airCode;
    semantic.materialId = summary.materialId = this->materialId;
    semantic.los = summary.los = this->los;
    semantic.materialColorValid = summary.materialColorValid = this->colorValid;
    semantic.materialColor = summary.materialColor = this->color;
    if (this->shaderChanged) {
	semantic.materialShader = std::move(this->semanticShader);
	summary.materialShader = std::move(this->summaryShader);
    }
    compact_note_semantic_change(entry);
}

/* Prepare only changed scalar/style records and retained intent. Geometry,
 * request ownership and path indexes stay put. */
class BObolCompactEntryPublication {
public:
    struct Change {
	Change(BObolCompactInstanceEntry &target, size_t index) : entry(target), ordinal(index),
	    authoredVisible(target.authoredVisible), visible(target.visible), selectable(target.selectable),
	    selected(target.selected), authoredHighlighted(target.authoredHighlighted), highlighted(target.highlighted),
	    presentationVisibleValid(target.presentationVisibleValid), presentationVisible(target.presentationVisible),
	    presentationHighlightedValid(target.presentationHighlightedValid), presentationHighlighted(target.presentationHighlighted),
	    presentationTransparencyValid(target.presentationTransparencyValid), presentationTransparency(target.presentationTransparency),
	    normalStyle(target.normalStyle), selectedStyle(target.selectedStyle), highlightedStyle(target.highlightedStyle)
	{}
	int setMetadata(int region, int air, int materialId, int los, SbBool colorValid,
	    const SbColor &color, const SbString &shader)
	{
	    const auto &current = this->entry.semantic;
	    this->materialColorChanged = current.materialColorValid != colorValid ||
		(colorValid && !database_source_color_equal(current.materialColor, color));
	    if (current.regionId == region && current.airCode == air && current.materialId == materialId &&
		current.los == los && !this->materialColorChanged && current.materialShader == shader)
		return 0;
	    this->metadata.emplace(current, region, air, materialId, los, colorValid, color, shader.getString());
	    this->appearanceChanged = true;
	    return 1;
	}
	void commit(const SoBRLDatabaseSource &source)
	{
	    if (this->entry.visible != this->visible)
		this->entry.visibilityRevision = compact_next_revision(this->entry.visibilityRevision);
	    if (this->entry.selected != this->selected || this->entry.highlighted != this->highlighted ||
		this->entry.selectable != this->selectable)
		this->entry.selectionRevision = compact_next_revision(this->entry.selectionRevision);
	    if (this->appearanceChanged)
		this->entry.appearanceRevision = compact_next_revision(this->entry.appearanceRevision);
	    this->entry.authoredVisible = this->authoredVisible;
	    this->entry.visible = this->visible;
	    this->entry.selectable = this->selectable;
	    this->entry.selected = this->selected;
	    this->entry.authoredHighlighted = this->authoredHighlighted;
	    this->entry.highlighted = this->highlighted;
	    this->entry.presentationVisibleValid = this->presentationVisibleValid;
	    this->entry.presentationVisible = this->presentationVisible;
	    this->entry.presentationHighlightedValid = this->presentationHighlightedValid;
	    this->entry.presentationHighlighted = this->presentationHighlighted;
	    this->entry.presentationTransparencyValid = this->presentationTransparencyValid;
	    this->entry.presentationTransparency = this->presentationTransparency;
	    this->entry.normalStyle = this->normalStyle;
	    this->entry.selectedStyle = this->selectedStyle;
	    this->entry.highlightedStyle = this->highlightedStyle;
	    if (this->metadata)
		this->metadata->commit(this->entry);
	    if (this->materialColorChanged)
		compact_set_material_styles(this->entry, source, true);
	    this->entry.style = compact_effective_style(this->entry);
	    compact_sync_shape_display_summary(this->entry);
	}
	BObolCompactInstanceEntry &entry;
	size_t ordinal;
	SbBool authoredVisible, visible, selectable, selected, authoredHighlighted, highlighted;
	SbBool presentationVisibleValid, presentationVisible, presentationHighlightedValid, presentationHighlighted;
	SbBool presentationTransparencyValid;
	float presentationTransparency;
	Obol::InstanceStyle normalStyle, selectedStyle, highlightedStyle;
	std::optional<BObolCompactMetadataEdit> metadata;
	bool appearanceChanged = false, materialColorChanged = false;
	bool publishDisplay = true;
    };
    explicit BObolCompactEntryPublication(SoBRLDatabaseSource &target) :
	source(target), index(*target.d->compactIndex)
    {}
    template <typename Edit>
    void stage(size_t ordinal, Edit edit)
    {
	Change change(this->index.entries[ordinal], ordinal);
	const int changed = edit(change);
	if (!changed)
	    return;
	this->changes.push_back(std::move(change));
	this->count += changed;
    }
    SbBool allows(size_t ordinal)
    {
	if (!this->source.d->compactVisibilityFrontierActive)
	    return TRUE;
	if (this->allowed.empty())
	    this->allowed = compact_visibility_frontier_mask(this->index, *this->source.d, 0);
	return this->allowed[ordinal];
    }
    int publish(bool notify, SoBRLDatabaseSource::PublicationCommit committed = nullptr, void *context = nullptr)
    {
	if (this->changes.empty())
	    return 0;
	return this->publish(notify, [] {}, committed, context);
    }
    template <typename Commit>
    int publish(bool notify, Commit commit, SoBRLDatabaseSource::PublicationCommit committed = nullptr, void *context = nullptr)
    {
	// Traversal order follows semantic paths; membership lookup uses dense ordinals.
	std::vector<Change *> ordered;
	std::vector<size_t> changedEntries, visibilityEntries;
	ordered.reserve(this->changes.size());
	changedEntries.reserve(this->changes.size());
	for (auto &change : this->changes) {
	    ordered.push_back(&change);
	    if (change.publishDisplay)
		changedEntries.push_back(change.ordinal);
	    if (change.visible != change.entry.visible)
		visibilityEntries.push_back(change.ordinal);
	}
	std::sort(ordered.begin(), ordered.end(), [](const Change *a, const Change *b) { return a->ordinal < b->ordinal; });
	const auto hidden = [](const auto &entry) { return !entry.visible; };
	const auto selected = [](const auto &entry) { return entry.selected != FALSE; };
	const auto unpickable = [](const auto &entry) { return !entry.selectable; };
	const bool hiddenChanged = this->prepareMembership(this->index.hiddenInstances, hidden);
	const bool selectedChanged = this->prepareMembership(this->index.selectedInstances, selected);
	const bool unpickableChanged = this->prepareMembership(this->index.unpickableInstances, unpickable);
	const SbBool notifySelection = this->clearSourceSelection && this->source.selected.isNotifyEnabled();

	// Capacity growth and all string/path work precede this allocation-free commit.
	commit();
	if (this->clearSourceSelection) {
	    this->source.selected.enableNotify(FALSE);
	    this->source.selected = FALSE;
	    this->source.selected.enableNotify(notifySelection);
	}
	if (hiddenChanged)
	    this->commitMembership(this->index.hiddenInstances, ordered, hidden);
	if (selectedChanged)
	    this->commitMembership(this->index.selectedInstances, ordered, selected);
	if (unpickableChanged)
	    this->commitMembership(this->index.unpickableInstances, ordered, unpickable);
	for (auto &change : this->changes) {
	    change.commit(this->source);
	    compact_sync_instance_style(this->index, change.ordinal);
	}
	if (this->clearSourceSelection || !changedEntries.empty())
	    this->source.markCompiledAssemblyDirty();
	if (this->clearSourceSelection)
	    this->source.markCadBatchDirty();
	else if (!changedEntries.empty())
	    this->source.markCadBatchDirty(changedEntries);
	this->source.markDisplayMeshLodVisibilityDirty(visibilityEntries);
	if (committed)
	    committed(context);
	if (!visibilityEntries.empty() && getenv("BOBOL_LOD_TRACE_SOURCE_CONTRACT"))
	    bu_log("BObol LoD source contract visibility delta source=%p path=%s changed=%zu revision=%llu\n",
		static_cast<void *>(&this->source), this->source.path.getValue().getString(), visibilityEntries.size(),
		static_cast<unsigned long long>(this->source.getDisplayMeshLodVisibilityRevision()));
	std::exception_ptr failure;
	if (notifySelection) {
	    try { this->source.selected.touch(); }
	    catch (...) { failure = std::current_exception(); }
	}
	if (notify) {
	    try { this->source.touch(); }
	    catch (...) { if (!failure) failure = std::current_exception(); }
	}
	if (failure)
	    std::rethrow_exception(failure);
	return this->count;
    }
    static int setOverride(SoBRLDatabaseSource &source, const char *path, BObolCompactPathMatch match,
	BObolCompactOccurrenceRegistryState::PresentationOverride::Property property, SbBool state, float transparency);
    static int clearHighlightOverrides(SoBRLDatabaseSource &source);
    static int setFrontier(SoBRLDatabaseSource &source, SbBool active, SbBool defaultVisible,
	const std::vector<SbString> &paths, const std::vector<SbBool> *states,
	SoBRLDatabaseSource::PublicationCommit committed = nullptr,
	void *context = nullptr);
    static int setSelection(SoBRLDatabaseSource &source,
	const std::vector<SbString> &paths,
	SoBRLDatabaseSource::PublicationCommit committed = nullptr,
	void *context = nullptr);
    static int selectionDelta(SoBRLDatabaseSource &source,
	const std::vector<SbString> &added, const std::vector<SbString> &removed,
	SoBRLDatabaseSource::PublicationCommit committed = nullptr,
	void *context = nullptr);
private:
    void stageOverrides(size_t ordinal,
	const std::vector<BObolCompactOccurrenceRegistryState::PresentationOverride> &rules, bool visibility);
    void stageSelection(const std::unordered_map<size_t, SbBool> &targets);
    void retireAggregateSelection(const std::vector<SbString> &paths)
    {
	this->clearSourceSelection = !this->index.entries.empty() && !paths.empty() && this->source.selected.getValue();
    }
    template <typename Member>
    bool prepareMembership(std::vector<Obol::InstanceId> &instances, Member member)
    {
	bool changed = false;
	size_t additions = 0;
	for (const auto &change : this->changes) {
	    const bool next = member(change), old = member(change.entry);
	    changed = changed || next != old;
	    if (next && !old)
		++additions;
	}
	if (additions)
	    instances.reserve(instances.size() + additions);
	return changed;
    }
    template <typename Member>
    void commitMembership(std::vector<Obol::InstanceId> &instances,
	const std::vector<Change *> &ordered, Member member)
    {
	instances.erase(std::remove_if(instances.begin(), instances.end(), [&](const Obol::InstanceId &instance) {
	    const auto ordinal = this->index.entryIndex.find(instance);
	    if (ordinal == this->index.entryIndex.end())
		return false;
	    const auto change = std::lower_bound(ordered.begin(), ordered.end(), ordinal->second,
		[](const Change *candidate, size_t value) { return candidate->ordinal < value; });
	    return change != ordered.end() && (*change)->ordinal == ordinal->second && !member(**change);
	}), instances.end());
	for (const auto &change : this->changes)
	    if (member(change) && !member(change.entry))
		instances.push_back(change.entry.instance);
    }
    SoBRLDatabaseSource &source;
    BObolCompactInstanceIndex &index;
    std::vector<Change> changes;
    std::vector<SbBool> allowed;
    int count = 0;
    bool clearSourceSelection = false;
};

int
SoBRLDatabaseSource::setCompactInstanceDisplayStateForPath(const char *queryPath,
    SbBool includeDescendants,
    int visibleValid, SbBool nextVisible,
    int selectedValid, SbBool nextSelected,
    int highlightedValid, SbBool nextHighlighted)
{
    const char *query = database_source_skip_leading_slash(
	queryPath ? queryPath : "");
    const bool leafQuery = query[0] && !strchr(query, '/') &&
	!strchr(query, '@');
    const BObolCompactPathMatch match = leafQuery ?
	BOBOL_COMPACT_PATH_OBJECT :
	(includeDescendants ? BOBOL_COMPACT_PATH_SUBTREE :
	 BOBOL_COMPACT_PATH_EXACT);
    return this->setCompactInstanceDisplayStateForPathMatch(queryPath,
	match, visibleValid, nextVisible, selectedValid, nextSelected,
	highlightedValid, nextHighlighted);
}

int
SoBRLDatabaseSource::setCompactInstanceDisplayStateForPath(const char *queryPath,
    SbBool includeDescendants,
    int visibleValid, SbBool nextVisible,
    int selectedValid, SbBool nextSelected,
    int highlightedValid, SbBool nextHighlighted,
    PublicationCommit committed, void *context)
{
    const char *query = database_source_skip_leading_slash(
	queryPath ? queryPath : "");
    const bool leafQuery = query[0] && !strchr(query, '/') &&
	!strchr(query, '@');
    const BObolCompactPathMatch match = leafQuery ?
	BOBOL_COMPACT_PATH_OBJECT :
	(includeDescendants ? BOBOL_COMPACT_PATH_SUBTREE :
	 BOBOL_COMPACT_PATH_EXACT);
    return this->setCompactInstanceDisplayStateForPathMatch(queryPath,
	match, visibleValid, nextVisible, selectedValid, nextSelected,
	highlightedValid, nextHighlighted, committed, context);
}

int
SoBRLDatabaseSource::setCompactInstanceDisplayStateForPathMatch(
    const char *queryPath, BObolCompactPathMatch match,
    int visibleValid, SbBool nextVisible, int selectedValid, SbBool nextSelected,
    int highlightedValid, SbBool nextHighlighted)
{
    return this->setCompactInstanceDisplayStateForPathMatch(queryPath, match,
	visibleValid, nextVisible, selectedValid, nextSelected,
	highlightedValid, nextHighlighted, nullptr, nullptr);
}

int
SoBRLDatabaseSource::setCompactInstanceDisplayStateForPathMatch(
    const char *queryPath, BObolCompactPathMatch match,
    int visibleValid, SbBool nextVisible, int selectedValid, SbBool nextSelected,
    int highlightedValid, SbBool nextHighlighted,
    PublicationCommit committed, void *context)
{
    if (!this->d->compactIndex || (match != BOBOL_COMPACT_PATH_EXACT &&
	match != BOBOL_COMPACT_PATH_SUBTREE && match != BOBOL_COMPACT_PATH_OBJECT)) return 0;
    BObolCompactEntryPublication publication(*this);
    compact_visit_entries_for_path_match(this->d->compactIndex, queryPath, match, [&](size_t ordinal) {
	publication.stage(ordinal, [&](auto &change) {
	    int changed = 0;
	    const auto &entry = change.entry;
	    if (visibleValid && !(nextVisible && compact_retired_overview(entry)) && change.authoredVisible != nextVisible) {
		change.authoredVisible = nextVisible;
		change.visible = entry.presentationVisibleValid ? entry.presentationVisible : nextVisible;
		if (this->d->compactVisibilityFrontierActive) change.visible = change.visible && publication.allows(ordinal);
		++changed;
	    }
	    if (selectedValid && change.selected != nextSelected) { change.selected = nextSelected; ++changed; }
	    if (highlightedValid && change.authoredHighlighted != nextHighlighted) {
		change.authoredHighlighted = nextHighlighted;
		change.highlighted = entry.presentationHighlightedValid ? entry.presentationHighlighted : nextHighlighted;
		++changed;
	    }
	    return changed;
	});
    });
    return publication.publish(true, committed, context);
}

int
SoBRLDatabaseSource::setCompactInstanceTransparencyForPathMatch(
    const char *queryPath, BObolCompactPathMatch match, float nextTransparency)
{
    if (!this->d->compactIndex || (match != BOBOL_COMPACT_PATH_EXACT &&
	match != BOBOL_COMPACT_PATH_SUBTREE && match != BOBOL_COMPACT_PATH_OBJECT)) return 0;
    const float alpha = 1.0f - std::max(0.0f, std::min(1.0f, nextTransparency));
    BObolCompactEntryPublication publication(*this);
    compact_visit_entries_for_path_match(this->d->compactIndex, queryPath, match, [&](size_t ordinal) {
	publication.stage(ordinal, [&](auto &change) {
	    if (!database_source_float_different(change.normalStyle.color[3], alpha) &&
		!database_source_float_different(change.selectedStyle.color[3], alpha) &&
		!database_source_float_different(change.highlightedStyle.color[3], alpha)) return 0;
	    change.normalStyle.color[3] = change.selectedStyle.color[3] = change.highlightedStyle.color[3] = alpha;
	    change.appearanceChanged = true;
	    return 1;
	});
    });
    return publication.publish(true);
}

static bool
compact_presentation_path_matches(const BObolCompactInstanceEntry &entry,
	const SbString &queryPath, BObolCompactPathMatch match)
{
    const char *candidate = database_source_skip_leading_slash(
	entry.semantic.path.getString());
    const char *query = database_source_skip_leading_slash(
	queryPath.getString());
    if (!candidate || !query || !query[0])
	return !query || !query[0];
    if (match == BOBOL_COMPACT_PATH_OBJECT)
	return database_source_leaf_component(entry.semantic.path) ==
	    database_source_leaf_component(queryPath);

    const size_t queryLength = strlen(query);
    if (bu_strncmp(candidate, query, queryLength))
	return false;
    const char suffix = candidate[queryLength];
    if (match == BOBOL_COMPACT_PATH_EXACT)
	return suffix == '\0';
    return suffix == '\0' || suffix == '/' || suffix == '@';
}


bool
compact_presentation_override_same_key(
    const BObolCompactOccurrenceRegistryState::PresentationOverride &left,
    const BObolCompactOccurrenceRegistryState::PresentationOverride &right)
{
    return left.property == right.property && left.match == right.match &&
	database_source_string_equal(left.path, right.path.getString());
}


BObolCompactPresentationOverrideState
compact_presentation_override_state(const BObolCompactInstanceEntry &entry,
    const std::vector<BObolCompactOccurrenceRegistryState::PresentationOverride> &rules)
{
    BObolCompactPresentationOverrideState state;
    for (const auto &rule : rules) {
	if (!compact_presentation_path_matches(entry, rule.path, rule.match))
	    continue;
	switch (rule.property) {
	    case BObolCompactOccurrenceRegistryState::PresentationOverride::VISIBILITY:
		state.visibleValid = TRUE;
		state.visible = rule.state;
		break;
	    case BObolCompactOccurrenceRegistryState::PresentationOverride::HIGHLIGHT:
		state.highlightedValid = TRUE;
		state.highlighted = rule.state;
		break;
	    case BObolCompactOccurrenceRegistryState::PresentationOverride::TRANSPARENCY:
		state.transparencyValid = TRUE;
		state.transparency = rule.transparency;
		break;
	}
    }
    return state;
}


template <typename Presentation>
static void
compact_apply_presentation_rules(Presentation &presentation,
    const BObolCompactInstanceEntry &entry,
    const std::vector<BObolCompactOccurrenceRegistryState::PresentationOverride> &rules)
{
    const BObolCompactPresentationOverrideState state =
	compact_presentation_override_state(entry, rules);
    presentation.presentationVisibleValid = state.visibleValid;
    presentation.presentationVisible = state.visible;
    presentation.presentationHighlightedValid = state.highlightedValid;
    presentation.presentationHighlighted = state.highlighted;
    presentation.presentationTransparencyValid = state.transparencyValid;
    presentation.presentationTransparency = state.transparency;
}


/* Indexed frontier traversal serves existing populations. An unpublished leaf
 * needs the same path semantics before any live lookup can resolve it. */
static bool
compact_entry_matches_frontier(const BObolCompactInstanceEntry &entry,
    const SbString &path)
{
    return database_source_path_matches_frontier(entry.semantic.path, path.getString());
}

template <typename Visitor>
static void
compact_visit_selected_paths(const BObolCompactInstanceIndex &index,
    const std::vector<SbString> &paths, Visitor visit)
{
    for (const SbString &selectedPath : paths) {
	const char *query = database_source_skip_leading_slash(selectedPath.getString());
	if (!query[0])
	    continue;
	compact_visit_entries_for_path(&index, query, TRUE, visit);
	const std::string semantic = database_source_db_path_without_instance_suffixes(query);
	if (!semantic.empty() && semantic != query)
	    compact_visit_entries_for_path(&index, semantic.c_str(), TRUE, visit);
    }
}

static std::vector<SbBool>
compact_frontier_mask(const BObolCompactInstanceIndex &index, SbBool active, SbBool defaultVisible,
    const std::vector<SbString> &paths, const std::vector<SbBool> &states, size_t firstEntry)
{
    const size_t candidateCount = index.entries.size() - firstEntry;
    std::vector<SbBool> allowed(candidateCount, active ? defaultVisible : TRUE);
    if (active) {
	for (size_t i = 0; i < paths.size(); ++i) {
	    const SbBool visible = i < states.size() ? states[i] : TRUE;
	    compact_visit_entries_for_path(&index, paths[i].getString(), TRUE,
		[&allowed, visible, firstEntry](size_t ordinal) {
		    if (ordinal >= firstEntry && ordinal - firstEntry < allowed.size())
			allowed[ordinal - firstEntry] = visible;
		});
	}
    }
    return allowed;
}

std::vector<SbBool>
compact_visibility_frontier_mask(const BObolCompactInstanceIndex &index,
    const BObolCompactOccurrenceRegistryState &state, size_t firstEntry)
{
    return compact_frontier_mask(index, state.compactVisibilityFrontierActive,
	state.compactVisibilityFrontierDefault, state.compactVisibilityFrontier,
	state.compactVisibilityFrontierStates, firstEntry);
}

static void
compact_update_effective_presentation(BObolCompactInstanceEntry &entry,
    SbBool allowed, SbBool selected)
{
    const SbBool visible = compact_effective_authored_visibility(entry) && allowed;
    if (visible != entry.visible) {
	entry.visible = visible;
	entry.visibilityRevision = compact_next_revision(entry.visibilityRevision);
    }
    if (selected != entry.selected) {
	entry.selected = selected;
	entry.selectionRevision = compact_next_revision(entry.selectionRevision);
    }
    const SbBool highlighted = compact_effective_highlight(entry);
    if (highlighted != entry.highlighted) {
	entry.highlighted = highlighted;
	entry.selectionRevision = compact_next_revision(entry.selectionRevision);
    }
    const Obol::InstanceStyle style = compact_effective_style(entry);
    if (!compact_style_equal(style, entry.style)) {
	entry.style = style;
	entry.appearanceRevision = compact_next_revision(entry.appearanceRevision);
    }
}

void
compact_prepare_occurrence_presentation(BObolCompactInstanceEntry &entry,
    const BObolCompactOccurrenceRegistryState &state)
{
    compact_apply_presentation_rules(entry, entry, state.compactPresentationOverrides);
    SbBool allowed = state.compactVisibilityFrontierActive ?
	state.compactVisibilityFrontierDefault : TRUE;
    if (state.compactVisibilityFrontierActive) {
	for (size_t i = 0; i < state.compactVisibilityFrontier.size(); ++i) {
	    if (compact_entry_matches_frontier(entry, state.compactVisibilityFrontier[i]))
		allowed = i < state.compactVisibilityFrontierStates.size() ?
		    state.compactVisibilityFrontierStates[i] : TRUE;
	}
    }
    SbBool selected = entry.selected;
    for (const SbString &path : state.compactSelectedPaths) {
	const char *query = database_source_skip_leading_slash(path.getString());
	if (!query[0])
	    continue;
	bool matches = compact_entry_matches_frontier(entry, path);
	if (!matches) {
	    const std::string semantic = database_source_db_path_without_instance_suffixes(query);
	    matches = !semantic.empty() && semantic != query &&
		compact_entry_matches_frontier(entry, SbString(semantic.c_str()));
	}
	if (matches)
	    selected = TRUE;
    }
    compact_update_effective_presentation(entry, allowed, selected);
}

void
compact_prepare_registry_presentation(BObolCompactInstanceIndex &index,
    const BObolCompactOccurrenceRegistryState &state)
{
    /* Complete replacement uses the existing path indexes. Scanning every
     * selection path for every occurrence would make large selections
     * quadratic. The temporary masks replace the live reapplication masks. */
    const std::vector<SbBool> allowed = compact_visibility_frontier_mask(index, state, 0);
    std::vector<unsigned char> selected(index.entries.size(), 0);
    compact_visit_selected_paths(index, state.compactSelectedPaths,
	[&selected](size_t entryIndex) { selected[entryIndex] = 1; });
    index.hiddenInstances.clear();
    index.selectedInstances.clear();
    index.unpickableInstances.clear();
    for (size_t i = 0; i < index.entries.size(); ++i) {
	BObolCompactInstanceEntry &entry = index.entries[i];
	compact_apply_presentation_rules(entry, entry, state.compactPresentationOverrides);
	compact_update_effective_presentation(entry, allowed[i], selected[i] ? TRUE : FALSE);
	compact_sync_shape_summary_state(entry);
	index.instances[i].record.style = entry.style;
	if (!entry.visible)
	    index.hiddenInstances.push_back(entry.instance);
	if (entry.selected)
	    index.selectedInstances.push_back(entry.instance);
	if (!entry.selectable)
	    index.unpickableInstances.push_back(entry.instance);
    }
}

void
BObolCompactEntryPublication::stageOverrides(size_t ordinal,
    const std::vector<BObolCompactOccurrenceRegistryState::PresentationOverride> &rules, bool visibility)
{
    this->stage(ordinal, [&](Change &change) {
	const auto &entry = change.entry;
	compact_apply_presentation_rules(change, entry, rules);
	if (visibility)
	    change.visible = compact_presentation_visibility(change, compact_retired_overview(entry)) && this->allows(ordinal);
	change.highlighted = compact_presentation_highlight(change);
	change.appearanceChanged = !compact_style_equal(compact_presentation_style(change), entry.style);
	change.publishDisplay = change.visible != entry.visible || change.highlighted != entry.highlighted || change.appearanceChanged;
	return change.publishDisplay || change.presentationVisibleValid != entry.presentationVisibleValid ||
	    change.presentationVisible != entry.presentationVisible ||
	    change.presentationHighlightedValid != entry.presentationHighlightedValid ||
	    change.presentationHighlighted != entry.presentationHighlighted ||
	    change.presentationTransparencyValid != entry.presentationTransparencyValid ||
	    database_source_float_different(change.presentationTransparency, entry.presentationTransparency);
    });
}

int
BObolCompactEntryPublication::setOverride(SoBRLDatabaseSource &source, const char *path, BObolCompactPathMatch match,
    BObolCompactOccurrenceRegistryState::PresentationOverride::Property property, SbBool state, float transparency)
{
    using Rule = BObolCompactOccurrenceRegistryState::PresentationOverride;
    if (match != BOBOL_COMPACT_PATH_EXACT && match != BOBOL_COMPACT_PATH_SUBTREE && match != BOBOL_COMPACT_PATH_OBJECT)
	return 0;
    const char *query = database_source_skip_leading_slash(path ? path : "");
    const float normalizedTransparency = std::max(0.0f, std::min(1.0f, transparency));
    auto &current = source.d->compactPresentationOverrides;
    if (!current.empty()) {
	const auto &last = current.back();
	if (last.property == property && last.match == match && database_source_string_equal(last.path, query) &&
	    (property == Rule::TRANSPARENCY ? !database_source_float_different(last.transparency, normalizedTransparency) : last.state == state))
	    return 0;
    }
    Rule next;
    next.property = property;
    next.path = query;
    next.match = match;
    next.state = state;
    next.transparency = normalizedTransparency;
    std::vector<Rule> prepared;
    prepared.reserve(current.size() + 1);
    for (const auto &rule : current)
	if (!compact_presentation_override_same_key(rule, next)) prepared.push_back(rule);
    prepared.push_back(std::move(next));
    if (!source.d->compactIndex) {
	current.swap(prepared);
	source.touch();
	return 1;
    }
    BObolCompactEntryPublication publication(source);
    // Moving one rule to the end can only alter occurrences matching its key.
    compact_visit_entries_for_path_match(source.d->compactIndex, query, match, [&](size_t ordinal) {
	publication.stageOverrides(ordinal, prepared, property == Rule::VISIBILITY);
    });
    publication.publish(true, [&] { current.swap(prepared); });
    return 1;
}

int
BObolCompactEntryPublication::clearHighlightOverrides(SoBRLDatabaseSource &source)
{
    using Rule = BObolCompactOccurrenceRegistryState::PresentationOverride;
    auto &current = source.d->compactPresentationOverrides;
    if (std::none_of(current.begin(), current.end(), [](const Rule &rule) { return rule.property == Rule::HIGHLIGHT; }))
	return 0;
    std::vector<Rule> prepared;
    prepared.reserve(current.size());
    std::unordered_set<size_t> affected;
    for (const auto &rule : current) {
	if (rule.property != Rule::HIGHLIGHT) {
	    prepared.push_back(rule);
	    continue;
	}
	compact_visit_entries_for_path_match(source.d->compactIndex, rule.path.getString(), rule.match,
	    [&](size_t ordinal) { affected.insert(ordinal); });
    }
    if (!source.d->compactIndex) {
	current.swap(prepared);
	source.touch();
	return 1;
    }
    BObolCompactEntryPublication publication(source);
    for (size_t ordinal : affected) publication.stageOverrides(ordinal, prepared, false);
    publication.publish(true, [&] { current.swap(prepared); });
    return 1;
}

int
SoBRLDatabaseSource::setCompactInstanceVisibilityOverrideForPathMatch(
    const char *queryPath, BObolCompactPathMatch match, SbBool nextVisible)
{
    return BObolCompactEntryPublication::setOverride(*this, queryPath, match,
	BObolCompactOccurrenceRegistryState::PresentationOverride::VISIBILITY, nextVisible, 0.0f);
}

int
SoBRLDatabaseSource::setCompactInstanceHighlightOverrideForPathMatch(
    const char *queryPath, BObolCompactPathMatch match, SbBool nextHighlighted)
{
    return BObolCompactEntryPublication::setOverride(*this, queryPath, match,
	BObolCompactOccurrenceRegistryState::PresentationOverride::HIGHLIGHT, nextHighlighted, 0.0f);
}

int
SoBRLDatabaseSource::setCompactInstanceTransparencyOverrideForPathMatch(
    const char *queryPath, BObolCompactPathMatch match, float nextTransparency)
{
    return BObolCompactEntryPublication::setOverride(*this, queryPath, match,
	BObolCompactOccurrenceRegistryState::PresentationOverride::TRANSPARENCY, FALSE, nextTransparency);
}

int
SoBRLDatabaseSource::clearCompactInstanceHighlightOverrides(void)
{
    return BObolCompactEntryPublication::clearHighlightOverrides(*this);
}

int
BObolCompactEntryPublication::setFrontier(SoBRLDatabaseSource &source, SbBool active, SbBool defaultVisible,
    const std::vector<SbString> &paths, const std::vector<SbBool> *states,
    SoBRLDatabaseSource::PublicationCommit committed, void *context)
{
    auto &state = *source.d;
    if (states && states->size() != paths.size()) return 0;
    bool same = state.compactVisibilityFrontierActive == active && state.compactVisibilityFrontierDefault == defaultVisible &&
	state.compactVisibilityFrontier.size() == paths.size() && state.compactVisibilityFrontierStates.size() == paths.size();
    for (size_t i = 0; same && i < paths.size(); ++i)
	same = database_source_string_equal(state.compactVisibilityFrontier[i], paths[i].getString()) &&
	    state.compactVisibilityFrontierStates[i] == (states ? (*states)[i] : TRUE);
    if (same) return 0;
    std::vector<SbString> preparedPaths(paths);
    std::vector<SbBool> preparedStates = states ? *states : std::vector<SbBool>(paths.size(), TRUE);
    const auto commit = [&] {
	state.compactVisibilityFrontier.swap(preparedPaths);
	state.compactVisibilityFrontierStates.swap(preparedStates);
	state.compactVisibilityFrontierActive = active;
	state.compactVisibilityFrontierDefault = defaultVisible;
    };
    if (!state.compactIndex) {
	commit();
	if (committed) committed(context);
	source.touch();
	return 1;
    }
    const auto allowed = compact_frontier_mask(*state.compactIndex, active, defaultVisible, preparedPaths, preparedStates, 0);
    BObolCompactEntryPublication publication(source);
    for (size_t ordinal = 0; ordinal < allowed.size(); ++ordinal) {
	publication.stage(ordinal, [&](Change &change) {
	    const SbBool visible = compact_effective_authored_visibility(change.entry) && allowed[ordinal];
	    if (change.visible == visible) return 0;
	    change.visible = visible;
	    return 1;
	});
    }
    const int changed = publication.publish(true, commit, committed, context);
    return changed ? changed : 1;
}

int
SoBRLDatabaseSource::setCompactInstanceVisibilityFrontier(const std::vector<SbString> &paths)
{
    return BObolCompactEntryPublication::setFrontier(*this, TRUE, FALSE, paths, nullptr);
}

int
SoBRLDatabaseSource::setCompactInstanceVisibilityFrontier(
    const std::vector<SbString> &paths,
    PublicationCommit committed, void *context)
{
    return BObolCompactEntryPublication::setFrontier(*this, TRUE, FALSE,
	paths, nullptr, committed, context);
}

int
SoBRLDatabaseSource::setCompactInstanceVisibilityOverrides(
    const std::vector<SbString> &paths, const std::vector<SbBool> &states)
{
    return BObolCompactEntryPublication::setFrontier(*this, TRUE, TRUE, paths, &states);
}

int
SoBRLDatabaseSource::setCompactInstanceVisibilityOverrides(
    const std::vector<SbString> &paths, const std::vector<SbBool> &states,
    PublicationCommit committed, void *context)
{
    return BObolCompactEntryPublication::setFrontier(*this, TRUE, TRUE,
	paths, &states, committed, context);
}

int
SoBRLDatabaseSource::clearCompactInstanceVisibilityFrontier(void)
{
    return BObolCompactEntryPublication::setFrontier(*this, FALSE, FALSE, {}, nullptr);
}

int
SoBRLDatabaseSource::clearCompactInstanceVisibilityFrontier(
    PublicationCommit committed, void *context)
{
    return BObolCompactEntryPublication::setFrontier(*this, FALSE, FALSE,
	{}, nullptr, committed, context);
}

SbBool
SoBRLDatabaseSource::hasCompactInstanceVisibilityFrontier(void) const
{
    return this->d->compactVisibilityFrontierActive;
}

void
BObolCompactEntryPublication::stageSelection(const std::unordered_map<size_t, SbBool> &targets)
{
    for (const auto &target : targets) {
	this->stage(target.first, [&](Change &change) {
	    if (change.selected == target.second) return 0;
	    change.selected = target.second;
	    return 1;
	});
    }
}

int
BObolCompactEntryPublication::setSelection(SoBRLDatabaseSource &source,
    const std::vector<SbString> &paths,
    SoBRLDatabaseSource::PublicationCommit committed, void *context)
{
    auto &current = source.d->compactSelectedPaths;
    bool same = current.size() == paths.size();
    for (size_t i = 0; same && i < paths.size(); ++i)
	same = database_source_string_equal(current[i], paths[i].getString());
    if (same) return 0;
    std::vector<SbString> prepared(paths);
    if (!source.d->compactIndex) {
	current.swap(prepared);
	if (committed) committed(context);
	source.touch();
	return 1;
    }
    const auto &index = *source.d->compactIndex;
    std::unordered_map<size_t, SbBool> targets;
    for (const auto &instance : index.selectedInstances) {
	auto found = index.entryIndex.find(instance);
	if (found != index.entryIndex.end()) targets.emplace(found->second, FALSE);
    }
    compact_visit_selected_paths(index, prepared, [&](size_t ordinal) { targets[ordinal] = TRUE; });
    if (getenv("BOBOL_SELECTION_DEBUG")) {
	const size_t selectedCount = std::count_if(targets.begin(), targets.end(),
	    [](const auto &target) { return target.second != FALSE; });
	bu_log("[obol-selection] source=%s entries=%zu paths=%zu touched=%zu selected=%zu\n",
	    source.path.getValue().getString(), index.entries.size(), prepared.size(), targets.size(), selectedCount);
	for (const auto &path : prepared)
	    bu_log("[obol-selection]   path=%s\n", path.getString());
    }
    BObolCompactEntryPublication publication(source);
    publication.stageSelection(targets);
    publication.retireAggregateSelection(prepared);
    const int changed = publication.publish(true,
	[&] { current.swap(prepared); }, committed, context);
    return changed ? changed : 1;
}

int
SoBRLDatabaseSource::syncCompactInstanceSelectedPaths(const std::vector<SbString> &paths)
{
    return BObolCompactEntryPublication::setSelection(*this, paths);
}

int
SoBRLDatabaseSource::syncCompactInstanceSelectedPaths(
    const std::vector<SbString> &paths,
    PublicationCommit committed, void *context)
{
    return BObolCompactEntryPublication::setSelection(*this, paths,
	committed, context);
}

int
BObolCompactEntryPublication::selectionDelta(SoBRLDatabaseSource &source,
    const std::vector<SbString> &added, const std::vector<SbString> &removed,
    SoBRLDatabaseSource::PublicationCommit committed, void *context)
{
    if (added.empty() && removed.empty()) return 0;
    int frontierChanged = 0;
    std::unordered_set<std::string> removedPaths;
    for (const auto &path : removed) {
	const char *query = database_source_skip_leading_slash(path.getString());
	if (query[0]) removedPaths.insert(query);
    }
    auto &current = source.d->compactSelectedPaths;
    std::vector<SbString> prepared;
    prepared.reserve(current.size() + added.size());
    std::unordered_set<std::string> retained;
    for (const auto &path : current) {
	const char *query = database_source_skip_leading_slash(path.getString());
	if (query[0] && removedPaths.find(query) != removedPaths.end()) {
	    ++frontierChanged;
	    continue;
	}
	prepared.push_back(path);
	if (query[0]) retained.insert(query);
    }
    for (const auto &path : added) {
	const char *query = database_source_skip_leading_slash(path.getString());
	if (query[0] && retained.insert(query).second) {
	    prepared.push_back(path);
	    ++frontierChanged;
	}
    }
    if (!source.d->compactIndex) {
	if (frontierChanged) {
	    current.swap(prepared);
	    if (committed) committed(context);
	    source.touch();
	}
	return frontierChanged;
    }
    const auto &index = *source.d->compactIndex;
    // A removed child path may still be covered by a retained parent. Resolve
    // affected entries against final intent, using the streaming path semantics.
    std::unordered_map<size_t, SbBool> targets;
    const auto affected = [&](size_t ordinal) { targets.emplace(ordinal, FALSE); };
    compact_visit_selected_paths(index, added, affected);
    compact_visit_selected_paths(index, removed, affected);
    compact_visit_selected_paths(index, prepared, [&](size_t ordinal) {
	auto found = targets.find(ordinal);
	if (found != targets.end()) found->second = TRUE;
    });
    BObolCompactEntryPublication publication(source);
    publication.stageSelection(targets);
    publication.retireAggregateSelection(prepared);
    if (!frontierChanged && publication.changes.empty() && !publication.clearSourceSelection) return 0;
    const int changed = publication.publish(true,
	[&] { current.swap(prepared); }, committed, context);
    return changed ? changed : std::max(frontierChanged, int(publication.clearSourceSelection));
}

int
SoBRLDatabaseSource::applyCompactInstanceSelectionDelta(
    const std::vector<SbString> &added, const std::vector<SbString> &removed)
{
    return BObolCompactEntryPublication::selectionDelta(*this, added, removed);
}

int
SoBRLDatabaseSource::applyCompactInstanceSelectionDelta(
    const std::vector<SbString> &added,
    const std::vector<SbString> &removed,
    PublicationCommit committed, void *context)
{
    return BObolCompactEntryPublication::selectionDelta(*this, added,
	removed, committed, context);
}


int
SoBRLDatabaseSource::setCompactInstanceSelectableForPath(
    const char *queryPath, SbBool includeDescendants, SbBool nextSelectable)
{
    if (!this->d->compactIndex) return 0;
    BObolCompactEntryPublication publication(*this);
    compact_visit_entries_for_path(this->d->compactIndex, queryPath, includeDescendants, [&](size_t ordinal) {
	publication.stage(ordinal, [&](auto &change) {
	    if (change.selectable == nextSelectable) return 0;
	    change.selectable = nextSelectable;
	    return 1;
	});
    });
    return publication.publish(true);
}

int
SoBRLDatabaseSource::setCompactInstanceRegionIdForPath(
    const char *queryPath, SbBool includeDescendants, int regionId)
{
    return this->setCompactInstanceRegionIdForPath(queryPath, includeDescendants, regionId, nullptr, nullptr);
}

int
SoBRLDatabaseSource::setCompactInstanceRegionIdForPath(
    const char *queryPath, SbBool includeDescendants, int regionId, PublicationCommit committed, void *context)
{
    if (!this->d->compactIndex) return 0;
    BObolCompactEntryPublication publication(*this);
    compact_visit_entries_for_path(this->d->compactIndex, queryPath, includeDescendants, [&](size_t ordinal) {
	publication.stage(ordinal, [&](auto &change) {
	    const auto &current = change.entry.semantic;
	    return change.setMetadata(regionId, current.airCode, current.materialId, current.los,
		current.materialColorValid, current.materialColor, current.materialShader);
	});
    });
    return publication.publish(false, committed, context);
}

int
SoBRLDatabaseSource::setCompactInstanceRegionMetadataForPath(
    const char *queryPath, SbBool includeDescendants, int regionId, int airCode, int materialId, int los)
{
    if (!this->d->compactIndex) return 0;
    BObolCompactEntryPublication publication(*this);
    compact_visit_entries_for_path(this->d->compactIndex, queryPath, includeDescendants, [&](size_t ordinal) {
	publication.stage(ordinal, [&](auto &change) {
	    const auto &current = change.entry.semantic;
	    return change.setMetadata(regionId, airCode, materialId, los,
		current.materialColorValid, current.materialColor, current.materialShader);
	});
    });
    return publication.publish(false);
}

int
SoBRLDatabaseSource::setCompactInstanceMetadataForPath(const char *queryPath,
    SbBool includeDescendants, int regionId, int airCode, int materialId, int los,
    SbBool nextMaterialColorValid, const SbColor &nextMaterialColor, const SbString &materialShader)
{
    return this->setCompactInstanceMetadataForPath(queryPath, includeDescendants, regionId, airCode,
	materialId, los, nextMaterialColorValid, nextMaterialColor, materialShader, nullptr, nullptr);
}

int
SoBRLDatabaseSource::setCompactInstanceMetadataForPath(const char *queryPath,
    SbBool includeDescendants, int regionId, int airCode, int materialId, int los,
    SbBool nextMaterialColorValid, const SbColor &nextMaterialColor, const SbString &materialShader,
    PublicationCommit committed, void *context)
{
    if (!this->d->compactIndex) return 0;
    const SbColor normalizedColor = nextMaterialColorValid ? nextMaterialColor : SbColor(1.0f, 1.0f, 1.0f);
    BObolCompactEntryPublication publication(*this);
    compact_visit_entries_for_path(this->d->compactIndex, queryPath, includeDescendants, [&](size_t ordinal) {
	publication.stage(ordinal, [&](auto &change) {
	    return change.setMetadata(regionId, airCode, materialId, los, nextMaterialColorValid, normalizedColor, materialShader);
	});
    });
    return publication.publish(false, committed, context);
}

int
SoBRLDatabaseSource::setCompactSubtractLineStyle(int nextLineStyle)
{
    if (!this->d->compactIndex) return 0;
    constexpr uint16_t solidPattern = 0xffffu, subtractPattern = 0xcf33u;
    const uint16_t pattern = nextLineStyle ? subtractPattern : solidPattern;
    BObolCompactEntryPublication publication(*this);
    for (size_t ordinal = 0; ordinal < this->d->compactIndex->entries.size(); ++ordinal) {
	const auto &entry = this->d->compactIndex->entries[ordinal];
	if (entry.booleanOperation != BOOLEAN_SUBTRACT || !entry.wireGeometry) continue;
	publication.stage(ordinal, [&](auto &change) {
	    if (change.normalStyle.linePattern == pattern && change.selectedStyle.linePattern == pattern &&
		change.highlightedStyle.linePattern == pattern) return 0;
	    change.normalStyle.linePattern = change.selectedStyle.linePattern = change.highlightedStyle.linePattern = pattern;
	    change.appearanceChanged = true;
	    return 1;
	});
    }
    return publication.publish(false);
}
