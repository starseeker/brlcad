/*     D A T A B A S E _ S O U R C E _ I N T E R A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file database_source_interaction.cpp
 *
 * Compact database-source export, measurement, and snapping operations.
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BExportAction.h"
#include "BObol/BMeasureAction.h"
#include "BObol/BSnapAction.h"
#include "cad_assembly_private.h"
#include "database_source_private.h"

#include <algorithm>
#include <map>
#include <math.h>
#include <utility>
#include <vector>

static SbMatrix
compact_entry_local_to_world(const BObolCompactInstanceEntry &entry,
			     const SbMatrix &parentToWorld)
{
    SbMatrix matrix = entry.localToSource;
    matrix.multRight(parentToWorld);
    return matrix;
}

template <typename Callback>
static void
compact_for_each_point(const BObolCompactInstanceEntry &entry,
	Callback callback)
{
    if (!entry.geometry || !entry.geometry->points)
	return;
    const Obol::PointRep &points = *entry.geometry->points;
    for (size_t i = 0; i < points.positions.size(); i++) {
	const int primitiveIndex = i < points.pointIds.size() ?
	    static_cast<int>(points.pointIds[i]) : static_cast<int>(i);
	const SbBool colorValid = i < points.colorValid.size() &&
	    points.colorValid[i] && i < points.colors.size();
	const SbBool scaleValid = i < points.scaleValid.size() &&
	    points.scaleValid[i] && i < points.scales.size();
	const SbBool normalValid = i < points.normalValid.size() &&
	    points.normalValid[i] && i < points.normals.size();
	callback(primitiveIndex, points.positions[i], colorValid,
	    colorValid ? points.colors[i] : SbColor(1.0f, 1.0f, 1.0f),
	    scaleValid, scaleValid ? points.scales[i] : 0.0f,
	    normalValid, normalValid ? points.normals[i] :
	    SbVec3f(0.0f, 0.0f, 1.0f));
    }
}

template <typename Callback>
static void
compact_for_each_wire_segment(const BObolCompactInstanceEntry &entry,
	Callback callback)
{
    if (!entry.geometry || !entry.geometry->wire)
	return;
    const Obol::WireRep &wire = *entry.geometry->wire;
    int primitiveIndex = 0;
    for (size_t i = 0; i + 1 < wire.segmentPoints.size(); i += 2)
	callback(primitiveIndex++, wire.segmentPoints[i],
	    wire.segmentPoints[i + 1]);
    for (const Obol::WirePolyline &polyline : wire.polylines) {
	for (size_t i = 1; i < polyline.points.size(); i++)
	    callback(primitiveIndex++, polyline.points[i - 1],
		polyline.points[i]);
    }
}

static float
compact_transform_point_scale(const SbMatrix &matrix, float scale)
{
    SbVec3f x;
    SbVec3f y;
    SbVec3f z;
    matrix.multDirMatrix(SbVec3f(scale, 0.0f, 0.0f), x);
    matrix.multDirMatrix(SbVec3f(0.0f, scale, 0.0f), y);
    matrix.multDirMatrix(SbVec3f(0.0f, 0.0f, scale), z);
    return std::max(x.length(), std::max(y.length(), z.length()));
}

template <typename Callback>
static void
compact_for_each_triangle(const BObolCompactInstanceEntry &entry,
	Callback callback)
{
    if (!entry.geometry || !entry.geometry->shaded)
	return;
    const Obol::TriMesh &mesh = *entry.geometry->shaded;
    int primitiveIndex = 0;
    for (size_t i = 0; i + 2 < mesh.indices.size(); i += 3) {
	const uint32_t ia = mesh.indices[i];
	const uint32_t ib = mesh.indices[i + 1];
	const uint32_t ic = mesh.indices[i + 2];
	if (ia >= mesh.positions.size() || ib >= mesh.positions.size() ||
	    ic >= mesh.positions.size())
	    continue;
	callback(primitiveIndex++, mesh.positions[ia], mesh.positions[ib],
	    mesh.positions[ic], static_cast<int>(ia), static_cast<int>(ib),
	    static_cast<int>(ic));
    }
}

static SbVec3f
compact_closest_point_on_segment(const SbVec3f &query,
				 const SbVec3f &a,
				 const SbVec3f &b)
{
    SbVec3f ab = b - a;
    float len2 = ab.sqrLength();
    if (len2 <= 0.0f)
	return a;

    float t = (query - a).dot(ab) / len2;
    if (t < 0.0f)
	t = 0.0f;
    if (t > 1.0f)
	t = 1.0f;
    return a + ab * t;
}

static SbVec3f
compact_closest_point_on_triangle(const SbVec3f &query,
				  const SbVec3f &a,
				  const SbVec3f &b,
				  const SbVec3f &c)
{
    SbVec3f ab = b - a;
    SbVec3f ac = c - a;
    SbVec3f ap = query - a;
    float d1 = ab.dot(ap);
    float d2 = ac.dot(ap);
    if (d1 <= 0.0f && d2 <= 0.0f)
	return a;

    SbVec3f bp = query - b;
    float d3 = ab.dot(bp);
    float d4 = ac.dot(bp);
    if (d3 >= 0.0f && d4 <= d3)
	return b;

    float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
	float v = d1 / (d1 - d3);
	return a + ab * v;
    }

    SbVec3f cp = query - c;
    float d5 = ab.dot(cp);
    float d6 = ac.dot(cp);
    if (d6 >= 0.0f && d5 <= d6)
	return c;

    float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
	float w = d2 / (d2 - d6);
	return a + ac * w;
    }

    float va = d3 * d6 - d5 * d4;
    if (va <= 0.0f && (d4 - d3) >= 0.0f &&
	(d5 - d6) >= 0.0f) {
	float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
	return b + (c - b) * w;
    }

    float denom = 1.0f / (va + vb + vc);
    float v = vb * denom;
    float w = vc * denom;
    return a + ab * v + ac * w;
}

struct compact_measure_segment_record {
    SbString path;
    SbString editIntentId;
    SbString editIntentRole;
    int primitiveIndex;
    SbVec3f a;
    SbVec3f b;
};

struct compact_measure_angle_record {
    size_t segmentA;
    size_t segmentB;
    SbVec3f shared;
    float degrees;
};

static const float COMPACT_MEASURE_ANGLE_VERTEX_TOLERANCE = 1.0e-5f;

struct compact_measure_endpoint_cell {
    long long x;
    long long y;
    long long z;

    bool operator<(const compact_measure_endpoint_cell &other) const
    {
	if (this->x != other.x)
	    return this->x < other.x;
	if (this->y != other.y)
	    return this->y < other.y;
	return this->z < other.z;
    }
};

static compact_measure_endpoint_cell
compact_measure_make_endpoint_cell(long long x, long long y, long long z)
{
    compact_measure_endpoint_cell cell;
    cell.x = x;
    cell.y = y;
    cell.z = z;
    return cell;
}

static long long
compact_measure_endpoint_cell_coord(float value)
{
    return static_cast<long long>(
	floor(static_cast<double>(value) /
	      static_cast<double>(COMPACT_MEASURE_ANGLE_VERTEX_TOLERANCE)));
}

static compact_measure_endpoint_cell
compact_measure_endpoint_cell_for_point(const SbVec3f &point)
{
    return compact_measure_make_endpoint_cell(
	compact_measure_endpoint_cell_coord(point[0]),
	compact_measure_endpoint_cell_coord(point[1]),
	compact_measure_endpoint_cell_coord(point[2]));
}

static void
compact_measure_collect_angle_endpoint_candidates(
    const std::map<compact_measure_endpoint_cell,
		   std::vector<size_t> > &endpointMap,
    const SbVec3f &point,
    std::vector<size_t> &candidates)
{
    compact_measure_endpoint_cell center =
	compact_measure_endpoint_cell_for_point(point);
    for (int dx = -1; dx <= 1; dx++) {
	for (int dy = -1; dy <= 1; dy++) {
	    for (int dz = -1; dz <= 1; dz++) {
		compact_measure_endpoint_cell key =
		    compact_measure_make_endpoint_cell(center.x + dx,
						       center.y + dy,
						       center.z + dz);
		std::map<compact_measure_endpoint_cell,
		    std::vector<size_t> >::const_iterator it =
			endpointMap.find(key);
		if (it == endpointMap.end())
		    continue;
		candidates.insert(candidates.end(), it->second.begin(),
				  it->second.end());
	    }
	}
    }
}

static void
compact_measure_add_angle_endpoint(
    std::map<compact_measure_endpoint_cell,
	     std::vector<size_t> > &endpointMap,
    const SbVec3f &point,
    size_t segmentIndex)
{
    endpointMap[compact_measure_endpoint_cell_for_point(point)].push_back(
	segmentIndex);
}

static float
compact_measure_clamp_float(float value, float minValue, float maxValue)
{
    if (value < minValue)
	return minValue;
    if (value > maxValue)
	return maxValue;
    return value;
}

static SbBool
compact_measure_same_point(const SbVec3f &a, const SbVec3f &b)
{
    return (a - b).length() <= COMPACT_MEASURE_ANGLE_VERTEX_TOLERANCE;
}

static SbBool
compact_measure_shared_segment_vertex(
    const compact_measure_segment_record &sa,
    const compact_measure_segment_record &sb,
    SbVec3f &shared,
    SbVec3f &otherA,
    SbVec3f &otherB)
{
    if (compact_measure_same_point(sa.a, sb.a)) {
	shared = sa.a;
	otherA = sa.b;
	otherB = sb.b;
	return TRUE;
    }
    if (compact_measure_same_point(sa.a, sb.b)) {
	shared = sa.a;
	otherA = sa.b;
	otherB = sb.a;
	return TRUE;
    }
    if (compact_measure_same_point(sa.b, sb.a)) {
	shared = sa.b;
	otherA = sa.a;
	otherB = sb.b;
	return TRUE;
    }
    if (compact_measure_same_point(sa.b, sb.b)) {
	shared = sa.b;
	otherA = sa.a;
	otherB = sb.a;
	return TRUE;
    }
    return FALSE;
}

static SbBool
compact_measure_segment_angle_degrees(const SbVec3f &shared,
				      const SbVec3f &otherA,
				      const SbVec3f &otherB,
				      float &angleDegrees)
{
    SbVec3f va = otherA - shared;
    SbVec3f vb = otherB - shared;
    float lenA = va.length();
    float lenB = vb.length();
    if (lenA <= 0.0f || lenB <= 0.0f)
	return FALSE;

    float cosAngle = compact_measure_clamp_float(
	va.dot(vb) / (lenA * lenB), -1.0f, 1.0f);
    angleDegrees = acosf(cosAngle) * (180.0f / 3.14159265358979323846f);
    return TRUE;
}

static void
compact_measure_collect_connected_angles(
    const std::vector<compact_measure_segment_record> &measuredSegments,
    std::vector<compact_measure_angle_record> &angleRecords)
{
    if (measuredSegments.size() < 2)
	return;

    std::map<compact_measure_endpoint_cell, std::vector<size_t> > endpointMap;
    std::vector<std::pair<size_t, size_t> > connectedPairs;
    std::vector<size_t> candidates;

    for (size_t i = 0; i < measuredSegments.size(); i++) {
	candidates.clear();
	compact_measure_collect_angle_endpoint_candidates(endpointMap,
	    measuredSegments[i].a, candidates);
	compact_measure_collect_angle_endpoint_candidates(endpointMap,
	    measuredSegments[i].b, candidates);
	std::sort(candidates.begin(), candidates.end());
	candidates.erase(std::unique(candidates.begin(), candidates.end()),
			 candidates.end());
	for (size_t j = 0; j < candidates.size(); j++)
	    connectedPairs.push_back(std::make_pair(candidates[j], i));
	compact_measure_add_angle_endpoint(endpointMap, measuredSegments[i].a,
	    i);
	compact_measure_add_angle_endpoint(endpointMap, measuredSegments[i].b,
	    i);
    }

    std::sort(connectedPairs.begin(), connectedPairs.end());
    connectedPairs.erase(std::unique(connectedPairs.begin(),
				     connectedPairs.end()),
			 connectedPairs.end());
    for (size_t i = 0; i < connectedPairs.size(); i++) {
	size_t segmentA = connectedPairs[i].first;
	size_t segmentB = connectedPairs[i].second;
	SbVec3f shared;
	SbVec3f otherA;
	SbVec3f otherB;
	float degrees = 0.0f;
	if (!compact_measure_shared_segment_vertex(measuredSegments[segmentA],
		measuredSegments[segmentB], shared, otherA, otherB))
	    continue;
	if (!compact_measure_segment_angle_degrees(shared, otherA, otherB,
		degrees))
	    continue;
	compact_measure_angle_record record;
	record.segmentA = segmentA;
	record.segmentB = segmentB;
	record.shared = shared;
	record.degrees = degrees;
	angleRecords.push_back(record);
    }
}

int
SoBRLDatabaseSource::exportCompactInstances(SoBRLExportAction *action,
	const SbMatrix &parentToWorld)
{
    if (!action || !this->hasCompactInstanceIndex())
	return 0;
    if (!this->visible.getValue())
	return 1;

    for (size_t i = 0; i < this->d->compactIndex->entries.size(); i++) {
	const BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[i];
	if (!entry.visible)
	    continue;

	const SbMatrix localToWorld =
	    compact_entry_local_to_world(entry, parentToWorld);
	const BObolRealizedShapeSummary &summary = entry.shapeSummary;
	if (entry.lodBacked &&
	    action->getGeometryPolicy() == SoBRLExportAction::FULL_DETAIL &&
	    entry.sourceMeshRequestValid) {
	    BObolSourceMeshRequest request = entry.sourceMeshRequest;
	    request.localToWorld = localToWorld;
	    action->sourceBackedFullDetailRequests.push_back(request);
	    action->skippedLodDisplayMeshCount++;
	    continue;
	}
	if (entry.pointGeometry) {
	    compact_for_each_point(entry,
		[&](int primitiveIndex, const SbVec3f &localPoint,
		    SbBool pointColorValid, const SbColor &pointColor,
		    SbBool pointScaleValid, float pointScale,
		    SbBool pointNormalValid, const SbVec3f &pointNormal) {
		    SbVec3f worldPoint;
		    SbVec3f worldNormal = pointNormal;
		    localToWorld.multVecMatrix(localPoint, worldPoint);
		    if (pointNormalValid) {
			const SbMatrix normalMatrix =
			    localToWorld.inverse().transpose();
			normalMatrix.multDirMatrix(pointNormal, worldNormal);
			worldNormal.normalize();
		    }
		action->appendPoint(summary.path, summary.sourceName,
			summary.sourceType, summary.sourceId, summary.regionId,
			summary.airCode, summary.materialId, summary.los,
			summary.materialColorValid, summary.materialColor,
			summary.materialShader, primitiveIndex, entry.selected,
			entry.highlighted, summary.ghosted, summary.hiddenLine,
			summary.editEmphasis, summary.editIntentId,
			summary.editIntentRole, summary.lodPolicy,
			summary.colorOverride, summary.color, pointColorValid,
			pointColor, pointScaleValid,
			pointScaleValid ? compact_transform_point_scale(
			    localToWorld, pointScale) : 0.0f,
		    pointNormalValid, worldNormal, worldPoint);
		action->applyLastPointMetadata(summary);
		});
	}
	if (entry.wireGeometry) {
	    compact_for_each_wire_segment(entry,
		[&](int segmentIndex, const SbVec3f &localA,
		    const SbVec3f &localB) {
		SbVec3f worldA;
		SbVec3f worldB;
		localToWorld.multVecMatrix(localA, worldA);
		localToWorld.multVecMatrix(localB, worldB);
		action->appendLine(summary.path, summary.sourceName,
		    summary.sourceType, summary.sourceId, summary.regionId,
		    summary.airCode, summary.materialId, summary.los,
		    summary.materialColorValid, summary.materialColor,
		    summary.materialShader, segmentIndex, entry.selected,
		    entry.highlighted, summary.ghosted, summary.hiddenLine,
		    summary.editEmphasis, summary.lineStyle,
		    summary.lineWidth, summary.editIntentId,
		    summary.editIntentRole, summary.lodPolicy,
		    summary.colorOverride, summary.color,
		    worldA, worldB);
		action->applyLastLineMetadata(summary);
	    });
	    continue;
	}

	if (entry.meshGeometry) {
	    compact_for_each_triangle(entry,
		[&](int triangleIndex, const SbVec3f &a, const SbVec3f &b,
		    const SbVec3f &c, int ia, int ib, int ic) {
		SbVec3f worldA;
		SbVec3f worldB;
		SbVec3f worldC;
		localToWorld.multVecMatrix(a, worldA);
		localToWorld.multVecMatrix(b, worldB);
		localToWorld.multVecMatrix(c, worldC);
		action->appendTriangle(summary.path, summary.sourceName,
		    summary.sourceType, summary.sourceId, summary.regionId,
		    summary.airCode, summary.materialId, summary.los,
		    summary.materialColorValid, summary.materialColor,
		    summary.materialShader, triangleIndex, ia, ib, ic,
		    entry.selected, entry.highlighted, summary.ghosted,
		    summary.hiddenLine, summary.editEmphasis,
		    summary.editIntentId, summary.editIntentRole,
		    summary.lodPolicy, summary.lodAvailable,
		    summary.lodActiveCut, summary.lodFaceCount,
		    summary.lodPointCount, summary.lodOriginalPointCount,
		    summary.lodNormalCount, summary.lodHasSnappedPoints,
		    summary.lodHasNormals, summary.lodBoundsMin,
		    summary.lodBoundsMax, summary.colorOverride, summary.color,
		    worldA, worldB, worldC);
		action->applyLastTriangleMetadata(summary);
	    });
	}
    }
    return 1;
}

int
SoBRLDatabaseSource::measureCompactInstances(SoBRLMeasureAction *action,
	const SbMatrix &parentToWorld)
{
    if (!action || !this->hasCompactInstanceIndex())
	return 0;
    if (!this->visible.getValue())
	return 1;

    for (size_t i = 0; i < this->d->compactIndex->entries.size(); i++) {
	const BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[i];
	if (!entry.visible ||
	    !action->selectionAllows(entry.selected) ||
	    !action->highlightAllows(entry.highlighted))
	    continue;

	const SbMatrix localToWorld =
	    compact_entry_local_to_world(entry, parentToWorld);
	SbBool measuredShape = FALSE;
	if (entry.wireGeometry) {
	    std::vector<compact_measure_segment_record> measuredSegments;
	    compact_for_each_wire_segment(entry,
		[&](int segmentIndex, const SbVec3f &localA,
		    const SbVec3f &localB) {
		SbVec3f pointA =
		    action->pointForCoordinateSpace(localToWorld, localA);
		SbVec3f pointB =
		    action->pointForCoordinateSpace(localToWorld, localB);
		action->measureSegment(entry.semantic.path,
		    entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, segmentIndex,
		    pointA, pointB);
		if (action->angleComputationEnabled) {
		    compact_measure_segment_record record;
		    record.path = entry.semantic.path;
		    record.editIntentId = entry.semantic.editIntentId;
		    record.editIntentRole = entry.semantic.editIntentRole;
		    record.primitiveIndex = segmentIndex;
		    record.a = pointA;
		    record.b = pointB;
		    measuredSegments.push_back(record);
		}
		measuredShape = TRUE;
	    });
	    if (action->angleComputationEnabled) {
		std::vector<compact_measure_angle_record> angleRecords;
		compact_measure_collect_connected_angles(measuredSegments,
		    angleRecords);
		for (size_t j = 0; j < angleRecords.size(); j++) {
		    const compact_measure_angle_record &angle =
			angleRecords[j];
		    const compact_measure_segment_record &segmentA =
			measuredSegments[angle.segmentA];
		    const compact_measure_segment_record &segmentB =
			measuredSegments[angle.segmentB];
		    action->considerAngle(segmentA.path, segmentA.editIntentId,
			segmentA.editIntentRole, segmentA.primitiveIndex,
			segmentB.primitiveIndex, angle.shared, angle.degrees);
		}
	    }
	} else if (entry.meshGeometry) {
	    compact_for_each_triangle(entry,
		[&](int triangleIndex, const SbVec3f &a, const SbVec3f &b,
		    const SbVec3f &c, int ia, int ib, int ic) {
		SbVec3f pointA =
		    action->pointForCoordinateSpace(localToWorld, a);
		SbVec3f pointB =
		    action->pointForCoordinateSpace(localToWorld, b);
		SbVec3f pointC =
		    action->pointForCoordinateSpace(localToWorld, c);
		action->measureTriangle(entry.semantic.path,
		    entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, pointA, pointB,
		    pointC, ia, ib, ic);
		measuredShape = TRUE;
	    });
	}
	if (measuredShape)
	    action->shapeCount++;
    }
    return 1;
}

int
SoBRLDatabaseSource::snapCompactInstances(SoBRLSnapAction *action,
	const SbMatrix &parentToWorld)
{
    if (!action || !this->hasCompactInstanceIndex())
	return 0;
    if (!this->visible.getValue())
	return 1;

    const SbVec3f query = action->queryPoint;
    for (size_t i = 0; i < this->d->compactIndex->entries.size(); i++) {
	const BObolCompactInstanceEntry &entry =
	    this->d->compactIndex->entries[i];
	if (!entry.visible || !entry.selectable ||
	    !action->selectionAllows(entry.selected))
	    continue;

	const SbMatrix localToWorld =
	    compact_entry_local_to_world(entry, parentToWorld);
	SbBox3f centerBox;
	centerBox.makeEmpty();
	if (entry.pointGeometry) {
	    compact_for_each_point(entry,
		[&](int primitiveIndex, const SbVec3f &localPoint,
		    SbBool, const SbColor &, SbBool, float,
		    SbBool, const SbVec3f &) {
		    const SbVec3f point = action->pointForCoordinateSpace(
			localToWorld, localPoint);
		    centerBox.extendBy(point);
		    action->consider(SoBRLSnapAction::ENDPOINT,
			entry.semantic.path, entry.semantic.editIntentId,
			entry.semantic.editIntentRole, primitiveIndex, query,
			point);
		});
	}
	if (entry.wireGeometry) {
	    compact_for_each_wire_segment(entry,
		[&](int segmentIndex, const SbVec3f &localA,
		    const SbVec3f &localB) {
		SbVec3f pointA =
		    action->pointForCoordinateSpace(localToWorld, localA);
		SbVec3f pointB =
		    action->pointForCoordinateSpace(localToWorld, localB);
		centerBox.extendBy(pointA);
		centerBox.extendBy(pointB);
		action->consider(SoBRLSnapAction::ENDPOINT,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, segmentIndex, query,
		    pointA);
		action->consider(SoBRLSnapAction::ENDPOINT,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, segmentIndex, query,
		    pointB);
		action->consider(SoBRLSnapAction::LINE_NEAREST,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, segmentIndex, query,
		    compact_closest_point_on_segment(query, pointA, pointB));
		action->consider(SoBRLSnapAction::MIDPOINT,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, segmentIndex, query,
		    (pointA + pointB) * 0.5f);
	    });
	} else if (entry.meshGeometry) {
	    compact_for_each_triangle(entry,
		[&](int triangleIndex, const SbVec3f &a, const SbVec3f &b,
		    const SbVec3f &c, int ia, int ib, int ic) {
		SbVec3f pointA =
		    action->pointForCoordinateSpace(localToWorld, a);
		SbVec3f pointB =
		    action->pointForCoordinateSpace(localToWorld, b);
		SbVec3f pointC =
		    action->pointForCoordinateSpace(localToWorld, c);
		centerBox.extendBy(pointA);
		centerBox.extendBy(pointB);
		centerBox.extendBy(pointC);
		action->consider(SoBRLSnapAction::VERTEX,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, query, pointA, ia);
		action->consider(SoBRLSnapAction::VERTEX,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, query, pointB, ib);
		action->consider(SoBRLSnapAction::VERTEX,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, query, pointC, ic);
		action->consider(SoBRLSnapAction::FACE_NEAREST,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, query,
		    compact_closest_point_on_triangle(query, pointA, pointB,
					      pointC));
		action->consider(SoBRLSnapAction::EDGE_NEAREST,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, query,
		    compact_closest_point_on_segment(query, pointA, pointB),
		    -1, 0, ia, ib);
		action->consider(SoBRLSnapAction::EDGE_NEAREST,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, query,
		    compact_closest_point_on_segment(query, pointB, pointC),
		    -1, 1, ib, ic);
		action->consider(SoBRLSnapAction::EDGE_NEAREST,
		    entry.semantic.path, entry.semantic.editIntentId,
		    entry.semantic.editIntentRole, triangleIndex, query,
		    compact_closest_point_on_segment(query, pointC, pointA),
		    -1, 2, ic, ia);
	    });
	}
	if (!centerBox.isEmpty()) {
	    const SoBRLCadAssembly::InstanceSemantic &semantic =
		entry.semantic;
	    action->consider(SoBRLSnapAction::CENTER, semantic.path,
		semantic.editIntentId, semantic.editIntentRole, -1, query,
		centerBox.getCenter());
	}
    }
    return 1;
}

int
SoBRLDatabaseSource::prepareCompiledAssembly(void)
{
    return this->syncCompiledAssembly();
}

SbBool
SoBRLDatabaseSource::hasCompiledAssembly(void) const
{
    return (this->d->compiledAssembly && this->d->compiledAssemblyActive) ?
	   TRUE : FALSE;
}

int
SoBRLDatabaseSource::getCompiledAssemblyPartCount(void) const
{
    if (!this->d->compiledAssembly || !this->d->compiledAssemblyActive)
	return 0;
    return static_cast<int>(this->d->compiledAssembly->partCount());
}

int
SoBRLDatabaseSource::getCompiledAssemblyInstanceCount(void) const
{
    if (!this->d->compiledAssembly || !this->d->compiledAssemblyActive)
	return 0;
    return static_cast<int>(this->d->compiledAssembly->instanceCount());
}
