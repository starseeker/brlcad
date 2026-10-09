/*              L O D _ C O V E R A G E _ P R E V I E W . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "lod_coverage_preview_private.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_set>
#include <vector>

namespace {

static bool
coverage_source_domain(const struct BObolMeshLodData &data,
	SbVec3f &minimum, SbVec3f &maximum, SbVec3f &extent,
	size_t &positiveExtentAxes)
{
    minimum = SbVec3f(static_cast<float>(data.bmin[X]),
	static_cast<float>(data.bmin[Y]), static_cast<float>(data.bmin[Z]));
    maximum = SbVec3f(static_cast<float>(data.bmax[X]),
	static_cast<float>(data.bmax[Y]), static_cast<float>(data.bmax[Z]));
    extent = maximum - minimum;
    positiveExtentAxes = 0;
    for (size_t axis = 0; axis < 3; ++axis) {
	if (!std::isfinite(minimum[axis]) || !std::isfinite(maximum[axis]) ||
	    extent[axis] < 0.0f)
	    return false;
	if (extent[axis] > 0.0f)
	    ++positiveExtentAxes;
    }
    return true;
}

static size_t
coverage_cell_coordinate(float value, float minimum, float extent,
	size_t cellAxis)
{
    if (!(extent > 0.0f) || cellAxis == 0)
	return 0;
    const float normalized = (value - minimum) / extent;
    const float scaled = normalized * static_cast<float>(cellAxis);
    return static_cast<size_t>(std::max(0.0f, std::min(
	static_cast<float>(cellAxis - 1), scaled)));
}

static bool
coverage_occupancy(const struct BObolMeshLodData &data,
	const SbVec3f &minimum, const SbVec3f &extent, size_t cellAxis,
	std::vector<uint8_t> &occupied)
{
    if (!data.points || !data.point_count || !cellAxis ||
	cellAxis > std::numeric_limits<size_t>::max() / cellAxis)
	return false;
    const size_t cellPlane = cellAxis * cellAxis;
    if (cellPlane > std::numeric_limits<size_t>::max() / cellAxis)
	return false;
    occupied.assign(cellPlane * cellAxis, 0u);
    for (size_t index = 0; index < data.point_count; ++index) {
	const point_t &point = data.points[index];
	const float value[3] = {
	    static_cast<float>(point[X]), static_cast<float>(point[Y]),
	    static_cast<float>(point[Z])
	};
	if (!std::isfinite(value[X]) || !std::isfinite(value[Y]) ||
	    !std::isfinite(value[Z]))
	    return false;
	const size_t x = coverage_cell_coordinate(
	    value[X], minimum[X], extent[X], cellAxis);
	const size_t y = coverage_cell_coordinate(
	    value[Y], minimum[Y], extent[Y], cellAxis);
	const size_t z = coverage_cell_coordinate(
	    value[Z], minimum[Z], extent[Z], cellAxis);
	occupied[x + cellAxis * (y + cellAxis * z)] = 1u;
    }
    return std::find(occupied.begin(), occupied.end(), 1u) != occupied.end();
}

static bool
coverage_voxel_counts(const std::vector<uint8_t> &occupied, size_t cellAxis,
	bool wire, bool shaded, BObolLodCounts &counts)
{
    counts.clear();
    if (occupied.empty() || !cellAxis)
	return false;
    const int faceAxes[6] = {X, X, Y, Y, Z, Z};
    const int faceDirections[6] = {-1, 1, -1, 1, -1, 1};
    const uint8_t faceCorners[6][4] = {
	{0, 4, 6, 2}, {1, 3, 7, 5}, {0, 1, 5, 4},
	{2, 6, 7, 3}, {0, 2, 3, 1}, {4, 5, 7, 6}
    };
    size_t exposedFaceCount = 0;
    std::unordered_set<uint32_t> wireEdges;
    if (wire)
	wireEdges.reserve(occupied.size() * 3u);
    for (size_t z = 0; z < cellAxis; ++z) {
	for (size_t y = 0; y < cellAxis; ++y) {
	    for (size_t x = 0; x < cellAxis; ++x) {
		if (!occupied[x + cellAxis * (y + cellAxis * z)])
		    continue;
		for (size_t face = 0; face < 6; ++face) {
		    size_t neighbor[3] = {x, y, z};
		    const size_t axis = static_cast<size_t>(faceAxes[face]);
		    const int direction = faceDirections[face];
		    if ((direction < 0 && neighbor[axis] > 0) ||
			(direction > 0 && neighbor[axis] + 1 < cellAxis)) {
			neighbor[axis] = static_cast<size_t>(
			    static_cast<int>(neighbor[axis]) + direction);
			if (occupied[neighbor[X] + cellAxis *
			    (neighbor[Y] + cellAxis * neighbor[Z])])
			    continue;
		    }
		    ++exposedFaceCount;
		    if (!wire)
			continue;
		    const uint8_t *faceCorner = faceCorners[face];
		    for (size_t corner = 0; corner < 4; ++corner) {
			const uint8_t firstCorner = faceCorner[corner];
			const uint8_t secondCorner = faceCorner[(corner + 1) % 4];
			const uint8_t changedAxis = firstCorner ^ secondCorner;
			const size_t edgeAxis = changedAxis == 1 ? X :
			    changedAxis == 2 ? Y : Z;
			const size_t gridX = x + std::min<size_t>(
			    firstCorner & 1u, secondCorner & 1u);
			const size_t gridY = y + std::min<size_t>(
			    (firstCorner >> 1u) & 1u,
			    (secondCorner >> 1u) & 1u);
			const size_t gridZ = z + std::min<size_t>(
			    (firstCorner >> 2u) & 1u,
			    (secondCorner >> 2u) & 1u);
			const size_t gridAxis = cellAxis + 1u;
			wireEdges.insert(static_cast<uint32_t>(edgeAxis + 3u *
			    (gridX + gridAxis * (gridY + gridAxis * gridZ))));
		    }
		}
	    }
	}
    }
    if (!exposedFaceCount)
	return false;
    if (shaded) {
	counts.faceCount = exposedFaceCount * 2u;
	counts.pointCount = exposedFaceCount * 4u;
	counts.originalPointCount = counts.pointCount;
    }
    if (wire)
	counts.lineCount = wireEdges.size();
    return counts.faceCount || counts.lineCount;
}

static bool
coverage_voxel_builder(const struct BObolMeshLodData &data,
	const SbBox3f &bounds, int drawMode, size_t cellAxis,
	size_t renderCostAllowance, Obol::PartGeometryBuilder &geometry,
	BObolLodCounts &counts, size_t &renderCost)
{
    SbVec3f sourceMinimum;
    SbVec3f sourceMaximum;
    SbVec3f sourceExtent;
    size_t positiveExtentAxes = 0;
    if (!coverage_source_domain(data, sourceMinimum, sourceMaximum,
	    sourceExtent, positiveExtentAxes) || positiveExtentAxes != 3)
	return false;

    std::vector<uint8_t> occupied;
    if (!coverage_occupancy(data, sourceMinimum, sourceExtent, cellAxis,
	    occupied))
	return false;

    const bool wire = drawMode == BOBOL_LOD_DRAW_WIRE ||
	drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    const bool shaded = drawMode == BOBOL_LOD_DRAW_SHADED ||
	drawMode == BOBOL_LOD_DRAW_SHADED_BOTS ||
	drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    if (!wire && !shaded)
	return false;
    if (!coverage_voxel_counts(
	    occupied, cellAxis, wire, shaded, counts))
	return false;
    renderCost = bobol_lod_render_cost_units(counts, drawMode, 1);
    if (renderCost > renderCostAllowance)
	return false;

    Obol::WireRep wireRep;
    Obol::TriMesh mesh;
    wireRep.bounds.setBounds(sourceMinimum, sourceMaximum);
    mesh.bounds = wireRep.bounds;
    std::unordered_set<uint32_t> wireEdges;
    if (wire)
	wireEdges.reserve(occupied.size() * 3u);
    const int faceAxes[6] = {X, X, Y, Y, Z, Z};
    const int faceDirections[6] = {-1, 1, -1, 1, -1, 1};
    const uint8_t faceCorners[6][4] = {
	{0, 4, 6, 2}, {1, 3, 7, 5}, {0, 1, 5, 4},
	{2, 6, 7, 3}, {0, 2, 3, 1}, {4, 5, 7, 6}
    };
    uint32_t edgeId = 0;
    for (size_t z = 0; z < cellAxis; ++z) {
	for (size_t y = 0; y < cellAxis; ++y) {
	    for (size_t x = 0; x < cellAxis; ++x) {
		if (!occupied[x + cellAxis * (y + cellAxis * z)])
		    continue;
		const SbVec3f cellMinimum = sourceMinimum + SbVec3f(
		    sourceExtent[0] * static_cast<float>(x) / cellAxis,
		    sourceExtent[1] * static_cast<float>(y) / cellAxis,
		    sourceExtent[2] * static_cast<float>(z) / cellAxis);
		const SbVec3f cellMaximum = sourceMinimum + SbVec3f(
		    sourceExtent[0] * static_cast<float>(x + 1) / cellAxis,
		    sourceExtent[1] * static_cast<float>(y + 1) / cellAxis,
		    sourceExtent[2] * static_cast<float>(z + 1) / cellAxis);
		const SbVec3f corners[8] = {
		    SbVec3f(cellMinimum[0], cellMinimum[1], cellMinimum[2]),
		    SbVec3f(cellMaximum[0], cellMinimum[1], cellMinimum[2]),
		    SbVec3f(cellMinimum[0], cellMaximum[1], cellMinimum[2]),
		    SbVec3f(cellMaximum[0], cellMaximum[1], cellMinimum[2]),
		    SbVec3f(cellMinimum[0], cellMinimum[1], cellMaximum[2]),
		    SbVec3f(cellMaximum[0], cellMinimum[1], cellMaximum[2]),
		    SbVec3f(cellMinimum[0], cellMaximum[1], cellMaximum[2]),
		    SbVec3f(cellMaximum[0], cellMaximum[1], cellMaximum[2])
		};
		for (size_t face = 0; face < 6; ++face) {
		    size_t neighbor[3] = {x, y, z};
		    const size_t axis = static_cast<size_t>(faceAxes[face]);
		    const int direction = faceDirections[face];
		    if ((direction < 0 && neighbor[axis] > 0) ||
			(direction > 0 && neighbor[axis] + 1 < cellAxis)) {
			neighbor[axis] = static_cast<size_t>(
			    static_cast<int>(neighbor[axis]) + direction);
			if (occupied[neighbor[X] + cellAxis *
			    (neighbor[Y] + cellAxis * neighbor[Z])])
			    continue;
		    }
		    const uint8_t *faceCorner = faceCorners[face];
		    if (shaded) {
			const uint32_t first = static_cast<uint32_t>(
			    mesh.positions.size());
			for (size_t corner = 0; corner < 4; ++corner)
			    mesh.positions.push_back(corners[faceCorner[corner]]);
			mesh.indices.insert(mesh.indices.end(), {
			    first, first + 1, first + 2,
			    first, first + 2, first + 3
			});
		    }
		    if (wire) {
			for (size_t corner = 0; corner < 4; ++corner) {
			    const uint8_t firstCorner = faceCorner[corner];
			    const uint8_t secondCorner =
				faceCorner[(corner + 1) % 4];
			    const uint8_t changedAxis = firstCorner ^ secondCorner;
			    const size_t edgeAxis = changedAxis == 1 ? X :
				changedAxis == 2 ? Y : Z;
			    const size_t gridX = x + std::min<size_t>(
				firstCorner & 1u, secondCorner & 1u);
			    const size_t gridY = y + std::min<size_t>(
				(firstCorner >> 1u) & 1u,
				(secondCorner >> 1u) & 1u);
			    const size_t gridZ = z + std::min<size_t>(
				(firstCorner >> 2u) & 1u,
				(secondCorner >> 2u) & 1u);
			    const size_t gridAxis = cellAxis + 1u;
			    const uint32_t edgeKey = static_cast<uint32_t>(
				edgeAxis + 3u * (gridX + gridAxis *
				    (gridY + gridAxis * gridZ)));
			    if (!wireEdges.insert(edgeKey).second)
				continue;
			    wireRep.segmentPoints.push_back(corners[firstCorner]);
			    wireRep.segmentPoints.push_back(corners[secondCorner]);
			    wireRep.segmentIds.push_back(edgeId++);
			}
		    }
		}
	    }
	}
    }
    if (mesh.indices.empty() && wireRep.segmentPoints.empty())
	return false;
    if (!mesh.indices.empty())
	geometry.shaded = std::move(mesh);
    if (!wireRep.segmentPoints.empty())
	geometry.wire = std::move(wireRep);
    if (!bounds.isEmpty())
	geometry.conservativeBounds = bounds;
    geometry.subpixelProxyEligible = true;
    return true;
}

static bool
coverage_point_builder(const struct BObolMeshLodData &data,
	const SbBox3f &bounds, int drawMode, size_t renderCostAllowance,
	Obol::PartGeometryBuilder &geometry, BObolLodCounts &counts,
	size_t &renderCost)
{
    SbVec3f sourceMinimum;
    SbVec3f sourceMaximum;
    SbVec3f sourceExtent;
    size_t positiveExtentAxes = 0;
    if (!coverage_source_domain(data, sourceMinimum, sourceMaximum,
	    sourceExtent, positiveExtentAxes))
	return false;
    (void)positiveExtentAxes;

    constexpr size_t pointCellAxis =
	BOBOL_MESH_LOD_COVERAGE_PREVIEW_CELL_AXIS / 4u;
    static_assert(pointCellAxis > 0,
	"coverage point fallback needs a nonzero grid");
    constexpr size_t pointCellCount =
	pointCellAxis * pointCellAxis * pointCellAxis;
    std::vector<uint8_t> occupied(pointCellCount, 0u);
    std::vector<SbVec3f> representatives(occupied.size());
    size_t representativeCount = 0;
    for (size_t index = 0; index < data.point_count; ++index) {
	const point_t &point = data.points[index];
	const SbVec3f value(static_cast<float>(point[X]),
	    static_cast<float>(point[Y]), static_cast<float>(point[Z]));
	if (!std::isfinite(value[X]) || !std::isfinite(value[Y]) ||
	    !std::isfinite(value[Z]))
	    return false;
	const size_t x = coverage_cell_coordinate(
	    value[X], sourceMinimum[X], sourceExtent[X], pointCellAxis);
	const size_t y = coverage_cell_coordinate(
	    value[Y], sourceMinimum[Y], sourceExtent[Y], pointCellAxis);
	const size_t z = coverage_cell_coordinate(
	    value[Z], sourceMinimum[Z], sourceExtent[Z], pointCellAxis);
	const size_t cell = x + pointCellAxis * (y + pointCellAxis * z);
	if (occupied[cell])
	    continue;
	occupied[cell] = 1u;
	representatives[representativeCount++] = value;
    }
    representatives.resize(representativeCount);
    if (representatives.empty())
	return false;

    size_t admittedCount = representatives.size();
    BObolLodCounts candidateCounts;
    candidateCounts.pointCount = admittedCount;
    candidateCounts.originalPointCount = admittedCount;
    size_t candidateCost = bobol_lod_render_cost_units(
	candidateCounts, drawMode, 1);
    if (candidateCost > renderCostAllowance) {
	size_t low = 0;
	size_t high = admittedCount;
	while (low < high) {
	    const size_t middle = low + (high - low + 1u) / 2u;
	    candidateCounts.pointCount = middle;
	    candidateCounts.originalPointCount = middle;
	    if (bobol_lod_render_cost_units(candidateCounts, drawMode, 1) <=
		    renderCostAllowance)
		low = middle;
	    else
		high = middle - 1u;
	}
	admittedCount = low;
	if (!admittedCount)
	    return false;
	std::vector<SbVec3f> reduced;
	reduced.reserve(admittedCount);
	for (size_t index = 0; index < admittedCount; ++index) {
	    const size_t source = static_cast<size_t>(
		(static_cast<long double>(index) * representatives.size()) /
		admittedCount);
	    reduced.push_back(representatives[std::min(
		source, representatives.size() - 1u)]);
	}
	representatives = std::move(reduced);
	candidateCounts.pointCount = admittedCount;
	candidateCounts.originalPointCount = admittedCount;
	candidateCost = bobol_lod_render_cost_units(
	    candidateCounts, drawMode, 1);
    }

    Obol::PointRep points;
    points.positions = std::move(representatives);
    points.bounds.setBounds(sourceMinimum, sourceMaximum);
    geometry.points = std::move(points);
    if (!bounds.isEmpty())
	geometry.conservativeBounds = bounds;
    geometry.subpixelProxyEligible = true;
    counts = candidateCounts;
    renderCost = candidateCost;
    return true;
}

} /* namespace */

bool
bobol_lod_build_coverage_preview(
    const struct BObolMeshLodData &data,
    const SbBox3f &conservativeBounds,
    int drawMode,
    size_t renderCostAllowance,
    BObolLodCoveragePreviewBuild &result)
{
    result = BObolLodCoveragePreviewBuild();
    if (!data.points || !data.point_count)
	return false;

    SbVec3f sourceMinimum;
    SbVec3f sourceMaximum;
    SbVec3f sourceExtent;
    size_t positiveExtentAxes = 0;
    if (!coverage_source_domain(data, sourceMinimum, sourceMaximum,
	    sourceExtent, positiveExtentAxes))
	return false;

    constexpr std::array<size_t, 3> candidateAxes = {
	BOBOL_MESH_LOD_COVERAGE_PREVIEW_CELL_AXIS,
	BOBOL_MESH_LOD_COVERAGE_PREVIEW_CELL_AXIS / 2u,
	BOBOL_MESH_LOD_COVERAGE_PREVIEW_CELL_AXIS / 4u
    };
    static_assert(BOBOL_MESH_LOD_COVERAGE_PREVIEW_CELL_AXIS % 4u == 0,
	"coverage preview grids must nest exactly");
    if (positiveExtentAxes == 3 && drawMode != BOBOL_LOD_DRAW_POINTS) {
	for (const size_t cellAxis : candidateAxes) {
	    Obol::PartGeometryBuilder geometry;
	    BObolLodCounts counts;
	    size_t cost = 0;
	    if (!coverage_voxel_builder(data, conservativeBounds, drawMode,
		    cellAxis, renderCostAllowance, geometry, counts, cost))
		continue;
	    result.geometry = std::move(geometry);
	    result.counts = counts;
	    result.renderCost = cost;
	    result.cellAxis = cellAxis;
	    return true;
	}
    }

    if (!coverage_point_builder(data, conservativeBounds, drawMode,
	    renderCostAllowance, result.geometry, result.counts,
	    result.renderCost))
	return false;
    result.pointFallback = true;
    return true;
}
