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
    return std::max(current, this->d->compactExpectedInstanceCount);
}

SbBool
SoBRLDatabaseSource::getCompactSourceProfile(
    BObolCompactSourceProfile &profile) const
{
    profile = this->d->compactSourceProfile;
    return profile.valid;
}

void
SoBRLDatabaseSource::setCompactSourceProfile(
    const BObolCompactSourceProfile &profile)
{
    this->d->compactSourceProfile = profile;
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
	const SbMatrix &viewProjection,
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
	const SbBox3f localBounds = compact_part_geometry_bounds(entry.geometry);
	SbMatrix localToWorld = entry.localToSource;
	localToWorld.multRight(parentToWorld);
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
SoBRLDatabaseSource::setCompactInstanceDisplayStateForPathMatch(
    const char *queryPath, BObolCompactPathMatch match,
    int visibleValid, SbBool nextVisible,
    int selectedValid, SbBool nextSelected,
    int highlightedValid, SbBool nextHighlighted)
{
    if (!this->d->compactIndex)
	return 0;
    if (match != BOBOL_COMPACT_PATH_EXACT &&
	match != BOBOL_COMPACT_PATH_SUBTREE &&
	match != BOBOL_COMPACT_PATH_OBJECT)
	return 0;

    int changed = 0;
    bool anyVisibilityChanged = false;
    std::vector<size_t> visibilityChangedEntries;
    const bool frontierActive =
	this->d->compactVisibilityFrontierActive ? true : false;
    compact_visit_entries_for_path_match(this->d->compactIndex, queryPath,
	match, [this, visibleValid, nextVisible, selectedValid, nextSelected,
	highlightedValid, nextHighlighted, &changed, &anyVisibilityChanged,
	&visibilityChangedEntries, frontierActive](size_t entryIndex) {
	BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[entryIndex];
	bool visibilityChanged = false;
	bool selectionChanged = false;
	const bool retiredOverview =
	    BU_STR_EQUAL(entry.shapeSummary.recordRole.getString(),
		"lod-overview") && !entry.authoredVisible;
	if (visibleValid && !(nextVisible && retiredOverview) &&
	    entry.authoredVisible != nextVisible) {
	    entry.authoredVisible = nextVisible;
	    const SbBool effectiveVisible =
		compact_effective_authored_visibility(entry);
	    if (!frontierActive && entry.visible != effectiveVisible) {
		entry.visible = effectiveVisible;
		visibilityChanged = true;
		anyVisibilityChanged = true;
		visibilityChangedEntries.push_back(entryIndex);
	    }
	    changed++;
	}
	if (selectedValid && entry.selected != nextSelected) {
	    entry.selected = nextSelected;
	    selectionChanged = true;
	    changed++;
	}
	if (highlightedValid &&
	    entry.authoredHighlighted != nextHighlighted) {
	    entry.authoredHighlighted = nextHighlighted;
	    const SbBool effectiveHighlighted =
		compact_effective_highlight(entry);
	    if (entry.highlighted != effectiveHighlighted) {
		entry.highlighted = effectiveHighlighted;
		selectionChanged = true;
	    }
	    changed++;
	}
	if (visibilityChanged)
	    entry.visibilityRevision = compact_next_revision(
		entry.visibilityRevision);
	if (selectionChanged) {
	    entry.selectionRevision = compact_next_revision(
		entry.selectionRevision);
	    entry.style = compact_effective_style(entry);
	}
    });
    if (changed) {
	if (visibleValid && this->d->compactVisibilityFrontierActive)
	    (void)this->reapplyCompactInstanceVisibilityFrontier();
	else
	    this->rebuildCompactInstanceDisplayState(FALSE);
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
	if (anyVisibilityChanged &&
	    !this->d->compactVisibilityFrontierActive)
	    this->markDisplayMeshLodDirty(visibilityChangedEntries);
	this->touch();
    }
    return changed;
}


int
SoBRLDatabaseSource::setCompactInstanceTransparencyForPathMatch(
    const char *queryPath, BObolCompactPathMatch match,
    float nextTransparency)
{
    if (!this->d->compactIndex ||
	(match != BOBOL_COMPACT_PATH_EXACT &&
	 match != BOBOL_COMPACT_PATH_SUBTREE &&
	 match != BOBOL_COMPACT_PATH_OBJECT))
	return 0;
    nextTransparency = std::max(0.0f, std::min(1.0f, nextTransparency));
    const float alpha = 1.0f - nextTransparency;
    int changed = 0;
    compact_visit_entries_for_path_match(this->d->compactIndex, queryPath,
	match, [this, &changed, alpha](size_t entryIndex) {
	BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[entryIndex];
	if (!database_source_float_different(entry.normalStyle.color[3],
		alpha) &&
	    !database_source_float_different(entry.selectedStyle.color[3],
		alpha) &&
	    !database_source_float_different(entry.highlightedStyle.color[3],
		alpha))
	    return;
	entry.normalStyle.color[3] = alpha;
	entry.selectedStyle.color[3] = alpha;
	entry.highlightedStyle.color[3] = alpha;
	entry.style = compact_effective_style(entry);
	entry.appearanceRevision = compact_next_revision(
	    entry.appearanceRevision);
	changed++;
    });
    if (changed) {
	this->rebuildCompactInstanceDisplayState(FALSE);
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
	this->touch();
    }
    return changed;
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


static bool
compact_presentation_override_same_key(
    const BObolCompactOccurrenceRegistryState::PresentationOverride &left,
    const BObolCompactOccurrenceRegistryState::PresentationOverride &right)
{
    return left.property == right.property && left.match == right.match &&
	database_source_string_equal(left.path, right.path.getString());
}


static bool
compact_presentation_override_same_value(
    const BObolCompactOccurrenceRegistryState::PresentationOverride &left,
    const BObolCompactOccurrenceRegistryState::PresentationOverride &right)
{
    if (!compact_presentation_override_same_key(left, right))
	return false;
    if (left.property ==
	BObolCompactOccurrenceRegistryState::PresentationOverride::TRANSPARENCY)
	return !database_source_float_different(left.transparency,
	    right.transparency);
    return left.state == right.state;
}


int
SoBRLDatabaseSource::reapplyCompactInstancePresentationOverrides(
    size_t firstEntry)
{
    if (!this->d->compactIndex)
	return 0;
    BObolCompactInstanceIndex &index = *this->d->compactIndex;
    firstEntry = std::min(firstEntry, index.entries.size());
    int changed = 0;
    std::vector<size_t> styleChangedEntries;
    for (size_t entryIndex = firstEntry; entryIndex < index.entries.size();
	entryIndex++) {
	BObolCompactInstanceEntry &entry = index.entries[entryIndex];
	const SbBool oldPresentationVisibleValid =
	    entry.presentationVisibleValid;
	const SbBool oldPresentationVisible = entry.presentationVisible;
	const SbBool oldHighlighted = entry.highlighted;
	const SbBool oldPresentationTransparencyValid =
	    entry.presentationTransparencyValid;
	const float oldPresentationTransparency =
	    entry.presentationTransparency;
	const Obol::InstanceStyle oldStyle = entry.style;
	entry.presentationVisibleValid = FALSE;
	entry.presentationVisible = TRUE;
	entry.presentationHighlightedValid = FALSE;
	entry.presentationHighlighted = FALSE;
	entry.presentationTransparencyValid = FALSE;
	entry.presentationTransparency = 0.0f;
	for (const BObolCompactOccurrenceRegistryState::PresentationOverride &rule :
	     this->d->compactPresentationOverrides) {
	    if (!compact_presentation_path_matches(entry, rule.path,
		    rule.match))
		continue;
	    switch (rule.property) {
		case BObolCompactOccurrenceRegistryState::PresentationOverride::VISIBILITY:
		    entry.presentationVisibleValid = TRUE;
		    entry.presentationVisible = rule.state;
		    break;
		case BObolCompactOccurrenceRegistryState::PresentationOverride::HIGHLIGHT:
		    entry.presentationHighlightedValid = TRUE;
		    entry.presentationHighlighted = rule.state;
		    break;
		case BObolCompactOccurrenceRegistryState::PresentationOverride::TRANSPARENCY:
		    entry.presentationTransparencyValid = TRUE;
		    entry.presentationTransparency = rule.transparency;
		    break;
	    }
	}
	entry.highlighted = compact_effective_highlight(entry);
	entry.style = compact_effective_style(entry);
	if (entry.highlighted != oldHighlighted) {
	    entry.selectionRevision = compact_next_revision(
		entry.selectionRevision);
	    changed++;
	}
	if (!compact_style_equal(entry.style, oldStyle)) {
	    entry.appearanceRevision = compact_next_revision(
		entry.appearanceRevision);
	    styleChangedEntries.push_back(entryIndex);
	    changed++;
	}
	if (entry.presentationVisibleValid != oldPresentationVisibleValid ||
	    entry.presentationVisible != oldPresentationVisible ||
	    entry.presentationTransparencyValid !=
		oldPresentationTransparencyValid ||
	    database_source_float_different(entry.presentationTransparency,
		oldPresentationTransparency))
	    changed++;
	compact_sync_shape_summary_state(entry);
	if (entryIndex < index.instances.size())
	    index.instances[entryIndex].record.style = entry.style;
    }
    const int visibilityChanged =
	this->reapplyCompactInstanceVisibilityFrontier(firstEntry);
    if (!styleChangedEntries.empty()) {
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty(styleChangedEntries);
	this->touch();
    }
    return changed + visibilityChanged;
}


static int
compact_presentation_override_store(
    std::vector<BObolCompactOccurrenceRegistryState::PresentationOverride> &rules,
    const BObolCompactOccurrenceRegistryState::PresentationOverride &next)
{
    if (!rules.empty() &&
	compact_presentation_override_same_value(rules.back(), next))
	return 0;
    const std::vector<BObolCompactOccurrenceRegistryState::PresentationOverride>
	previous = rules;
    rules.erase(std::remove_if(rules.begin(), rules.end(),
	[&next](const BObolCompactOccurrenceRegistryState::PresentationOverride &rule) {
	    return compact_presentation_override_same_key(rule, next);
	}), rules.end());
    rules.push_back(next);
    if (previous.size() != rules.size())
	return 1;
    for (size_t i = 0; i < previous.size(); i++) {
	if (!compact_presentation_override_same_value(previous[i], rules[i]))
	    return 1;
    }
    return 0;
}


int
SoBRLDatabaseSource::setCompactInstanceVisibilityOverrideForPathMatch(
    const char *queryPath, BObolCompactPathMatch match, SbBool nextVisible)
{
    if (match != BOBOL_COMPACT_PATH_EXACT &&
	match != BOBOL_COMPACT_PATH_SUBTREE &&
	match != BOBOL_COMPACT_PATH_OBJECT)
	return 0;
    BObolCompactOccurrenceRegistryState::PresentationOverride rule;
    rule.property = BObolCompactOccurrenceRegistryState::PresentationOverride::VISIBILITY;
    rule.path = database_source_skip_leading_slash(
	queryPath ? queryPath : "");
    rule.match = match;
    rule.state = nextVisible;
    const int stored = compact_presentation_override_store(
	this->d->compactPresentationOverrides, rule);
    const int applied = stored ?
	this->reapplyCompactInstancePresentationOverrides(0) : 0;
    return stored || applied ? 1 : 0;
}


int
SoBRLDatabaseSource::setCompactInstanceHighlightOverrideForPathMatch(
    const char *queryPath, BObolCompactPathMatch match, SbBool nextHighlighted)
{
    if (match != BOBOL_COMPACT_PATH_EXACT &&
	match != BOBOL_COMPACT_PATH_SUBTREE &&
	match != BOBOL_COMPACT_PATH_OBJECT)
	return 0;
    BObolCompactOccurrenceRegistryState::PresentationOverride rule;
    rule.property = BObolCompactOccurrenceRegistryState::PresentationOverride::HIGHLIGHT;
    rule.path = database_source_skip_leading_slash(
	queryPath ? queryPath : "");
    rule.match = match;
    rule.state = nextHighlighted;
    const int stored = compact_presentation_override_store(
	this->d->compactPresentationOverrides, rule);
    const int applied = stored ?
	this->reapplyCompactInstancePresentationOverrides(0) : 0;
    return stored || applied ? 1 : 0;
}


int
SoBRLDatabaseSource::setCompactInstanceTransparencyOverrideForPathMatch(
    const char *queryPath, BObolCompactPathMatch match,
    float nextTransparency)
{
    if (match != BOBOL_COMPACT_PATH_EXACT &&
	match != BOBOL_COMPACT_PATH_SUBTREE &&
	match != BOBOL_COMPACT_PATH_OBJECT)
	return 0;
    BObolCompactOccurrenceRegistryState::PresentationOverride rule;
    rule.property = BObolCompactOccurrenceRegistryState::PresentationOverride::TRANSPARENCY;
    rule.path = database_source_skip_leading_slash(
	queryPath ? queryPath : "");
    rule.match = match;
    rule.transparency = std::max(0.0f,
	std::min(1.0f, nextTransparency));
    const int stored = compact_presentation_override_store(
	this->d->compactPresentationOverrides, rule);
    const int applied = stored ?
	this->reapplyCompactInstancePresentationOverrides(0) : 0;
    return stored || applied ? 1 : 0;
}


int
SoBRLDatabaseSource::clearCompactInstanceHighlightOverrides(void)
{
    const size_t previous = this->d->compactPresentationOverrides.size();
    this->d->compactPresentationOverrides.erase(std::remove_if(
	this->d->compactPresentationOverrides.begin(),
	this->d->compactPresentationOverrides.end(),
	[](const BObolCompactOccurrenceRegistryState::PresentationOverride &rule) {
	    return rule.property ==
		BObolCompactOccurrenceRegistryState::PresentationOverride::HIGHLIGHT;
	}), this->d->compactPresentationOverrides.end());
    if (previous == this->d->compactPresentationOverrides.size())
	return 0;
    (void)this->reapplyCompactInstancePresentationOverrides(0);
    return 1;
}


int
SoBRLDatabaseSource::reapplyCompactInstanceVisibilityFrontier(
    size_t firstEntry)
{
    if (!this->d->compactIndex)
	return 0;

    BObolCompactInstanceIndex &index = *this->d->compactIndex;
    firstEntry = std::min(firstEntry, index.entries.size());
    const size_t candidateCount = index.entries.size() - firstEntry;
    std::vector<SbBool> allowed(candidateCount,
	this->d->compactVisibilityFrontierActive ?
	    this->d->compactVisibilityFrontierDefault : TRUE);
    if (this->d->compactVisibilityFrontierActive) {
	for (size_t ruleIndex = 0;
	     ruleIndex < this->d->compactVisibilityFrontier.size();
	     ruleIndex++) {
	    const SbString &frontierPath =
		this->d->compactVisibilityFrontier[ruleIndex];
	    const SbBool ruleVisible =
		ruleIndex < this->d->compactVisibilityFrontierStates.size() ?
		this->d->compactVisibilityFrontierStates[ruleIndex] : TRUE;
	    compact_visit_entries_for_path(&index,
		frontierPath.getString(), TRUE,
		[&allowed, ruleVisible, firstEntry](size_t entryIndex) {
		    if (entryIndex >= firstEntry &&
			entryIndex - firstEntry < allowed.size())
			allowed[entryIndex - firstEntry] = ruleVisible;
		});
	}
    }

    int changed = 0;
    std::vector<size_t> changedEntries;
    for (size_t i = firstEntry; i < index.entries.size(); i++) {
	BObolCompactInstanceEntry &entry = index.entries[i];
	const SbBool nextVisible =
	    compact_effective_authored_visibility(entry) &&
	    allowed[i - firstEntry];
	if (entry.visible == nextVisible)
	    continue;
	entry.visible = nextVisible;
	entry.visibilityRevision = compact_next_revision(
	    entry.visibilityRevision);
	compact_sync_shape_summary_state(entry);
	changedEntries.push_back(i);
	changed++;
    }
    if (changed) {
	/* Visibility is an instance-set delta, not a geometry or style change.
	 * Rebuilding all display state here needlessly revisits every style and
	 * instance record and can make a one-occurrence erase look like a scene
	 * rebuild to the renderer. */
	if (changedEntries.size() > index.entries.size() / 4) {
	    index.hiddenInstances.clear();
	    for (const BObolCompactInstanceEntry &entry : index.entries) {
		if (!entry.visible)
		    index.hiddenInstances.push_back(entry.instance);
	    }
	} else {
	    std::unordered_set<Obol::InstanceId, std::hash<Obol::InstanceId>>
		changedInstances;
	    changedInstances.reserve(changedEntries.size());
	    for (const size_t entryIndex : changedEntries)
		changedInstances.insert(index.entries[entryIndex].instance);
	    index.hiddenInstances.erase(std::remove_if(
		index.hiddenInstances.begin(), index.hiddenInstances.end(),
		[&changedInstances](const Obol::InstanceId &instance) {
		    return changedInstances.find(instance) !=
			changedInstances.end();
		}), index.hiddenInstances.end());
	    for (const size_t entryIndex : changedEntries) {
		const BObolCompactInstanceEntry &entry =
		    index.entries[entryIndex];
		if (!entry.visible)
		    index.hiddenInstances.push_back(entry.instance);
	    }
	}
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty(changedEntries);
	/* A hierarchy visibility frontier changes which retained instances are
	 * presented; it does not change mesh inventory or immutable availability.
	 * Keep the view payload bound so an erase/redraw can hide and restore it in
	 * one presentation transaction.  Advancing the mesh-inventory revision
	 * here made the LoD delta pass remove hidden payloads and briefly commit an
	 * empty frame before an asynchronous redraw reattached the same resident
	 * asset.  Camera, policy, and resource-pressure revisions still perform the
	 * normal view allocation and may compact hidden data when warranted. */
	this->touch();
    }
    return changed;
}


int
SoBRLDatabaseSource::setCompactInstanceVisibilityFrontier(
    const std::vector<SbString> &paths)
{
    bool same = this->d->compactVisibilityFrontierActive &&
	this->d->compactVisibilityFrontierDefault == FALSE &&
	this->d->compactVisibilityFrontier.size() == paths.size() &&
	this->d->compactVisibilityFrontierStates.size() == paths.size();
    for (size_t i = 0; same && i < paths.size(); i++)
	same = database_source_string_equal(
	    this->d->compactVisibilityFrontier[i], paths[i].getString()) &&
	    this->d->compactVisibilityFrontierStates[i] == TRUE;
    if (same)
	return 0;

    this->d->compactVisibilityFrontier = paths;
    this->d->compactVisibilityFrontierStates.assign(paths.size(), TRUE);
    this->d->compactVisibilityFrontierDefault = FALSE;
    this->d->compactVisibilityFrontierActive = TRUE;
    const int changed = this->reapplyCompactInstanceVisibilityFrontier();
    return changed > 0 ? changed : 1;
}


int
SoBRLDatabaseSource::setCompactInstanceVisibilityOverrides(
    const std::vector<SbString> &paths,
    const std::vector<SbBool> &states)
{
    if (paths.size() != states.size())
	return 0;

    bool same = this->d->compactVisibilityFrontierActive &&
	this->d->compactVisibilityFrontierDefault == TRUE &&
	this->d->compactVisibilityFrontier.size() == paths.size() &&
	this->d->compactVisibilityFrontierStates.size() == states.size();
    for (size_t i = 0; same && i < paths.size(); i++) {
	same = database_source_string_equal(
	    this->d->compactVisibilityFrontier[i], paths[i].getString()) &&
	    this->d->compactVisibilityFrontierStates[i] == states[i];
    }
    if (same)
	return 0;

    this->d->compactVisibilityFrontier = paths;
    this->d->compactVisibilityFrontierStates = states;
    this->d->compactVisibilityFrontierDefault = TRUE;
    this->d->compactVisibilityFrontierActive = TRUE;
    const int changed = this->reapplyCompactInstanceVisibilityFrontier();
    return changed > 0 ? changed : 1;
}


int
SoBRLDatabaseSource::clearCompactInstanceVisibilityFrontier(void)
{
    if (!this->d->compactVisibilityFrontierActive)
	return 0;
    this->d->compactVisibilityFrontierActive = FALSE;
    this->d->compactVisibilityFrontierDefault = FALSE;
    this->d->compactVisibilityFrontier.clear();
    this->d->compactVisibilityFrontierStates.clear();
    const int changed = this->reapplyCompactInstanceVisibilityFrontier();
    return changed > 0 ? changed : 1;
}


SbBool
SoBRLDatabaseSource::hasCompactInstanceVisibilityFrontier(void) const
{
    return this->d->compactVisibilityFrontierActive;
}


int
SoBRLDatabaseSource::reapplyCompactInstanceSelectedPaths(
    size_t firstEntry)
{
    if (!this->d->compactIndex)
	return 0;

    BObolCompactInstanceIndex &index = *this->d->compactIndex;
    /* The compact records supersede the aggregate overview for selection
     * presentation.  Selection paths are reapplied below per occurrence; do
     * not leave the source proxy selected as a second, potentially lingering
     * white box. */
    if (!index.entries.empty() && !this->d->compactSelectedPaths.empty() &&
	this->selected.getValue())
	this->selected = FALSE;
    firstEntry = std::min(firstEntry, index.entries.size());
    const size_t candidateCount = index.entries.size() - firstEntry;
    std::vector<unsigned char> touched(candidateCount, 0);
    std::vector<unsigned char> selectedMask(candidateCount, 0);

    /* Existing selected entries are the only candidates that can become
     * deselected.  New selection paths contribute candidates through the
     * sorted path index, so selecting one row does not scan every occurrence
     * in a large assembly. */
    if (firstEntry == 0) {
	for (const Obol::InstanceId &instance : index.selectedInstances) {
	    auto found = index.entryIndex.find(instance);
	    if (found != index.entryIndex.end() &&
		found->second < touched.size())
		touched[found->second] = 1;
	}
    }

    const auto selectPath =
	[&index, &touched, &selectedMask, firstEntry](
	const char *queryPath) {
	if (!queryPath || !queryPath[0])
	    return;
	compact_visit_entries_for_path(&index, queryPath, TRUE,
	    [&touched, &selectedMask, firstEntry](size_t entryIndex) {
		if (entryIndex < firstEntry ||
		    entryIndex - firstEntry >= touched.size())
		    return;
		touched[entryIndex - firstEntry] = 1;
		selectedMask[entryIndex - firstEntry] = 1;
	    });
    };
    for (const SbString &selectedPath : this->d->compactSelectedPaths) {
	const char *selectedPathString = database_source_skip_leading_slash(
	    selectedPath.getString());
	if (!selectedPathString || !selectedPathString[0])
	    continue;
	selectPath(selectedPathString);
	const std::string semantic =
	    database_source_db_path_without_instance_suffixes(
		selectedPathString);
	if (!semantic.empty() && semantic != selectedPathString)
	    selectPath(semantic.c_str());
    }
    if (getenv("BOBOL_SELECTION_DEBUG")) {
	size_t touchedCount = 0;
	size_t selectedCount = 0;
	for (size_t i = 0; i < touched.size(); i++) {
	    touchedCount += touched[i] ? 1u : 0u;
	    selectedCount += selectedMask[i] ? 1u : 0u;
	}
	bu_log("[obol-selection] source=%s entries=%zu paths=%zu "
	    "touched=%zu selected=%zu\n",
	    this->path.getValue().getString(), index.entries.size(),
	    this->d->compactSelectedPaths.size(), touchedCount,
	    selectedCount);
	for (const SbString &selectedPath : this->d->compactSelectedPaths)
	    bu_log("[obol-selection]   path=%s\n",
		selectedPath.getString());
    }

    int changed = 0;
    std::vector<size_t> changedEntries;
    for (size_t i = firstEntry; i < index.entries.size(); i++) {
	if (!touched[i - firstEntry])
	    continue;
	BObolCompactInstanceEntry &entry = index.entries[i];
	const SbBool nextSelected =
	    selectedMask[i - firstEntry] ? TRUE : FALSE;
	if (entry.selected == nextSelected)
	    continue;
	entry.selected = nextSelected;
	entry.selectionRevision = compact_next_revision(
	    entry.selectionRevision);
	entry.style = compact_effective_style(entry);
	compact_sync_shape_summary_state(entry);
	if (i < index.instances.size() &&
	    index.instances[i].instance == entry.instance) {
	    index.instances[i].record.style = entry.style;
	} else {
	    auto found = std::find_if(index.instances.begin(),
		index.instances.end(),
		[&entry](const Obol::InstanceUpdate &update) {
		    return update.instance == entry.instance;
		});
	    if (found != index.instances.end())
		found->record.style = entry.style;
	}
	changedEntries.push_back(i);
	changed++;
    }
    if (changed) {
	std::unordered_set<Obol::InstanceId, std::hash<Obol::InstanceId>>
	    changedInstances;
	changedInstances.reserve(changedEntries.size());
	for (size_t entryIndex : changedEntries)
	    changedInstances.insert(index.entries[entryIndex].instance);

	std::vector<Obol::InstanceId> nextSelected;
	nextSelected.reserve(index.selectedInstances.size() +
	    changedEntries.size());
	for (const Obol::InstanceId &instance : index.selectedInstances) {
	    if (changedInstances.find(instance) == changedInstances.end())
		nextSelected.push_back(instance);
	}
	for (size_t entryIndex : changedEntries) {
	    const BObolCompactInstanceEntry &entry = index.entries[entryIndex];
	    if (entry.selected)
		nextSelected.push_back(entry.instance);
	}
	index.selectedInstances.swap(nextSelected);
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty(changedEntries);
	this->touch();
    }
    return changed;
}


int
SoBRLDatabaseSource::syncCompactInstanceSelectedPaths(
    const std::vector<SbString> &paths)
{
    bool same = this->d->compactSelectedPaths.size() == paths.size();
    for (size_t i = 0; same && i < paths.size(); i++)
	same = database_source_string_equal(
	    this->d->compactSelectedPaths[i], paths[i].getString());
    if (same)
	return 0;

    this->d->compactSelectedPaths = paths;
    const int changed = this->reapplyCompactInstanceSelectedPaths();
    return changed > 0 ? changed : 1;
}


int
SoBRLDatabaseSource::applyCompactInstanceSelectionDelta(
    const std::vector<SbString> &addedPaths,
    const std::vector<SbString> &removedPaths)
{
    int frontierChanged = 0;
    std::unordered_set<std::string> removedFrontier;
    removedFrontier.reserve(removedPaths.size());
    for (const SbString &removed : removedPaths) {
	const char *selectedPath = database_source_skip_leading_slash(
	    removed.getString());
	if (selectedPath && selectedPath[0])
	    removedFrontier.insert(selectedPath);
    }

    this->d->compactSelectedPaths.erase(std::remove_if(
	this->d->compactSelectedPaths.begin(),
	this->d->compactSelectedPaths.end(),
	[&removedFrontier, &frontierChanged](const SbString &existing) {
	    const char *selectedPath = database_source_skip_leading_slash(
		existing.getString());
	    if (!selectedPath || removedFrontier.find(selectedPath) ==
		removedFrontier.end())
		return false;
	    ++frontierChanged;
	    return true;
	}), this->d->compactSelectedPaths.end());

    std::unordered_set<std::string> existingFrontier;
    existingFrontier.reserve(this->d->compactSelectedPaths.size() +
	addedPaths.size());
    for (const SbString &existing : this->d->compactSelectedPaths) {
	const char *selectedPath = database_source_skip_leading_slash(
	    existing.getString());
	if (selectedPath && selectedPath[0])
	    existingFrontier.insert(selectedPath);
    }
    for (const SbString &added : addedPaths) {
	const char *selectedPath = database_source_skip_leading_slash(
	    added.getString());
	if (!selectedPath || !selectedPath[0])
	    continue;
	if (existingFrontier.insert(selectedPath).second) {
	    this->d->compactSelectedPaths.push_back(added);
	    frontierChanged++;
	}
    }

    if (!this->d->compactIndex) {
	if (frontierChanged)
	    this->touch();
	return frontierChanged;
    }
    BObolCompactInstanceIndex &index = *this->d->compactIndex;

    /* Resolve every path through the ordered occurrence index, but publish
     * the resulting style changes once.  Calling the single-path setter in a
     * loop rebuilt all N display records after each of K selected paths,
     * making a window selection O(N*K) on the GUI thread. */
    std::unordered_map<size_t, SbBool> targetStates;
    targetStates.reserve(addedPaths.size() + removedPaths.size());
    const auto setTargets = [&index, &targetStates](
	const std::vector<SbString> &paths, SbBool state) {
	for (const SbString &pathValue : paths) {
	    const char *selectedPath = pathValue.getString();
	    if (!selectedPath || !selectedPath[0])
		continue;
	    compact_visit_entries_for_path_match(&index, selectedPath,
		BOBOL_COMPACT_PATH_SUBTREE,
		[&targetStates, state](size_t entryIndex) {
		    targetStates[entryIndex] = state;
		});
	}
    };
    setTargets(removedPaths, FALSE);
    setTargets(addedPaths, TRUE);

    std::vector<size_t> changedEntries;
    changedEntries.reserve(targetStates.size());
    for (const auto &target : targetStates) {
	if (target.first >= index.entries.size())
	    continue;
	BObolCompactInstanceEntry &entry = index.entries[target.first];
	if (entry.selected == target.second)
	    continue;
	entry.selected = target.second;
	entry.selectionRevision = compact_next_revision(
	    entry.selectionRevision);
	entry.style = compact_effective_style(entry);
	compact_sync_shape_summary_state(entry);
	if (target.first < index.instances.size() &&
	    index.instances[target.first].instance == entry.instance) {
	    index.instances[target.first].record.style = entry.style;
	} else {
	    auto found = std::find_if(index.instances.begin(),
		index.instances.end(), [&entry](const Obol::InstanceUpdate &update) {
		    return update.instance == entry.instance;
		});
	    if (found != index.instances.end())
		found->record.style = entry.style;
	}
	changedEntries.push_back(target.first);
    }

    if (!changedEntries.empty()) {
	std::unordered_set<Obol::InstanceId, std::hash<Obol::InstanceId>>
	    changedInstances;
	changedInstances.reserve(changedEntries.size());
	for (size_t entryIndex : changedEntries)
	    changedInstances.insert(index.entries[entryIndex].instance);
	index.selectedInstances.erase(std::remove_if(
	    index.selectedInstances.begin(), index.selectedInstances.end(),
	    [&changedInstances](const Obol::InstanceId &instance) {
		return changedInstances.find(instance) != changedInstances.end();
	    }), index.selectedInstances.end());
	for (size_t entryIndex : changedEntries) {
	    const BObolCompactInstanceEntry &entry = index.entries[entryIndex];
	    if (entry.selected)
		index.selectedInstances.push_back(entry.instance);
	}
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty(changedEntries);
	this->touch();
    } else if (frontierChanged) {
	this->touch();
    }
    return !changedEntries.empty() ?
	static_cast<int>(changedEntries.size()) : frontierChanged;
}


int
SoBRLDatabaseSource::setCompactInstanceSelectableForPath(
    const char *queryPath, SbBool includeDescendants, SbBool nextSelectable)
{
    if (!this->d->compactIndex)
	return 0;
    int changed = 0;

    compact_visit_entries_for_path(this->d->compactIndex, queryPath,
	includeDescendants, [this, nextSelectable, &changed](size_t entryIndex) {
	BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[entryIndex];
	if (entry.selectable == nextSelectable)
	    return;
	entry.selectable = nextSelectable;
	entry.selectionRevision = compact_next_revision(
	    entry.selectionRevision);
	changed++;
    });
    if (changed) {
	this->rebuildCompactInstanceDisplayState(FALSE);
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
	this->touch();
    }
    return changed;
}

int
SoBRLDatabaseSource::setCompactInstanceRegionIdForPath(const char *queryPath,
    SbBool includeDescendants, int regionId)
{
    if (!this->d->compactIndex)
	return 0;
    int changed = 0;

    compact_visit_entries_for_path(this->d->compactIndex, queryPath,
	includeDescendants, [this, regionId, &changed](size_t entryIndex) {
	BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[entryIndex];
	if (entry.semantic.regionId == regionId)
	    return;
	entry.semantic.regionId = regionId;
	entry.appearanceRevision = compact_next_revision(
	    entry.appearanceRevision);
	compact_sync_shape_summary_state(entry);
	changed++;
    });
    if (changed) {
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
    }
    return changed;
}

int
SoBRLDatabaseSource::setCompactInstanceRegionMetadataForPath(
    const char *queryPath, SbBool includeDescendants, int regionId,
    int airCode, int materialId, int los)
{
    if (!this->d->compactIndex)
	return 0;
    int changed = 0;

    compact_visit_entries_for_path(this->d->compactIndex, queryPath,
	includeDescendants, [this, regionId, airCode, materialId, los,
	&changed](size_t entryIndex) {
	BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[entryIndex];
	if (entry.semantic.regionId == regionId &&
	     entry.semantic.airCode == airCode &&
	     entry.semantic.materialId == materialId &&
	     entry.semantic.los == los)
	    return;
	entry.semantic.regionId = regionId;
	entry.semantic.airCode = airCode;
	entry.semantic.materialId = materialId;
	entry.semantic.los = los;
	entry.appearanceRevision = compact_next_revision(
	    entry.appearanceRevision);
	compact_sync_shape_summary_state(entry);
	changed++;
    });
    if (changed) {
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
    }
    return changed;
}

int
SoBRLDatabaseSource::setCompactInstanceMetadataForPath(const char *queryPath,
    SbBool includeDescendants, int regionId, int airCode, int materialId,
    int los, SbBool nextMaterialColorValid, const SbColor &nextMaterialColor,
    const SbString &materialShader)
{
    if (!this->d->compactIndex)
	return 0;
    int changed = 0;

    compact_visit_entries_for_path(this->d->compactIndex, queryPath,
	includeDescendants, [this, regionId, airCode, materialId, los,
	nextMaterialColorValid, &nextMaterialColor, &materialShader,
	&changed](size_t entryIndex) {
	BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[entryIndex];
	const bool same = entry.semantic.regionId == regionId &&
	    entry.semantic.airCode == airCode &&
	    entry.semantic.materialId == materialId &&
	    entry.semantic.los == los &&
	    entry.semantic.materialColorValid == nextMaterialColorValid &&
	    (!nextMaterialColorValid || database_source_color_equal(
		entry.semantic.materialColor, nextMaterialColor)) &&
	    bu_strcmp(entry.semantic.materialShader.getString(),
		materialShader.getString()) == 0;
	if (same)
	    return;
	entry.semantic.regionId = regionId;
	entry.semantic.airCode = airCode;
	entry.semantic.materialId = materialId;
	entry.semantic.los = los;
	entry.semantic.materialColorValid = nextMaterialColorValid;
	entry.semantic.materialColor = nextMaterialColorValid ? nextMaterialColor :
	    SbColor(1.0f, 1.0f, 1.0f);
	entry.semantic.materialShader = materialShader;
	entry.normalStyle = compact_entry_style_from_source(this, entry,
	    FALSE, FALSE);
	entry.selectedStyle = compact_entry_style_from_source(this, entry,
	    TRUE, FALSE);
	entry.highlightedStyle = compact_entry_style_from_source(this, entry,
	    FALSE, TRUE);
	entry.style = compact_effective_style(entry);
	entry.appearanceRevision = compact_next_revision(
	    entry.appearanceRevision);
	compact_sync_shape_summary_state(entry);
	changed++;
    });
    if (changed) {
	this->rebuildCompactInstanceDisplayState(FALSE);
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
    }
    return changed;
}

int
SoBRLDatabaseSource::setCompactSubtractLineStyle(int nextLineStyle)
{
    if (!this->d->compactIndex)
	return 0;
    int changed = 0;
    const uint16_t pattern = nextLineStyle != 0 ? 0xcf33u : 0xffffu;
    for (BObolCompactInstanceEntry &entry : this->d->compactIndex->entries) {
	if (entry.booleanOperation != BOOLEAN_SUBTRACT || !entry.wireGeometry ||
	    entry.normalStyle.linePattern == pattern)
	    continue;
	entry.normalStyle.linePattern = pattern;
	entry.selectedStyle.linePattern = pattern;
	entry.highlightedStyle.linePattern = pattern;
	entry.style = compact_effective_style(entry);
	entry.appearanceRevision = compact_next_revision(
	    entry.appearanceRevision);
	changed++;
    }
    if (changed) {
	this->rebuildCompactInstanceDisplayState(FALSE);
	this->markCompiledAssemblyDirty();
	this->markCadBatchDirty();
    }
    return changed;
}
