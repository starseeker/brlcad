/*          T E S T _ L O D _ C O V E R A G E _ P R E V I E W . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "lod_coverage_preview_private.h"

#include "bu/app.h"

#include <Obol/cad/CadGeometry.h>
#include <Obol/cad/CadGeometryValidation.h>

#include <cstdint>
#include <cstdio>
#include <limits>
#include <vector>

static struct BObolMeshLodData
coverage_data(std::vector<fastf_t> &points, fastf_t xmax, fastf_t ymax,
	fastf_t zmax)
{
    struct BObolMeshLodData data = {};
    data.points = reinterpret_cast<const point_t *>(points.data());
    data.point_count = points.size() / 3u;
    data.points_orig = data.points;
    data.point_orig_count = data.point_count;
    VSET(data.bmin, 0.0, 0.0, 0.0);
    VSET(data.bmax, xmax, ymax, zmax);
    return data;
}

static int
test_adaptive_voxels(void)
{
    constexpr size_t axis = BOBOL_MESH_LOD_COVERAGE_PREVIEW_CELL_AXIS;
    std::vector<fastf_t> points;
    points.reserve(axis * axis * axis * 3u);
    for (size_t z = 0; z < axis; ++z) {
	for (size_t y = 0; y < axis; ++y) {
	    for (size_t x = 0; x < axis; ++x) {
		points.push_back(static_cast<fastf_t>(x) + 0.5);
		points.push_back(static_cast<fastf_t>(y) + 0.5);
		points.push_back(static_cast<fastf_t>(z) + 0.5);
	    }
	}
    }
    const struct BObolMeshLodData data = coverage_data(
	points, static_cast<fastf_t>(axis), static_cast<fastf_t>(axis),
	static_cast<fastf_t>(axis));
    const SbBox3f bounds(SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(static_cast<float>(axis), static_cast<float>(axis),
	    static_cast<float>(axis)));

    BObolLodCoveragePreviewBuild rich;
    if (!bobol_lod_build_coverage_preview(data, bounds,
	    BOBOL_LOD_DRAW_SHADED, SIZE_MAX, rich) || rich.pointFallback ||
	rich.cellAxis != axis || !rich.geometry.shaded ||
	rich.geometry.wire || !rich.counts.faceCount || !rich.renderCost) {
	std::fprintf(stderr, "FAIL: rich cold coverage grid was not selected\n");
	return 1;
    }
    const size_t richCost = rich.renderCost;
    if (!Obol::cadAdmitPartGeometry(std::move(rich.geometry))) {
	std::fprintf(stderr, "FAIL: rich cold coverage geometry was invalid\n");
	return 1;
    }

    BObolLodCoveragePreviewBuild medium;
    if (!bobol_lod_build_coverage_preview(data, bounds,
	    BOBOL_LOD_DRAW_SHADED, richCost - 1u, medium) ||
	medium.pointFallback || medium.cellAxis != axis / 2u ||
	medium.renderCost >= richCost) {
	std::fprintf(stderr, "FAIL: cold coverage did not adapt to 12^3\n");
	return 1;
    }
    const size_t mediumCost = medium.renderCost;

    BObolLodCoveragePreviewBuild coarse;
    if (!bobol_lod_build_coverage_preview(data, bounds,
	    BOBOL_LOD_DRAW_SHADED, mediumCost - 1u, coarse) ||
	coarse.pointFallback || coarse.cellAxis != axis / 4u ||
	coarse.renderCost >= mediumCost) {
	std::fprintf(stderr, "FAIL: cold coverage did not adapt to 6^3\n");
	return 1;
    }

    BObolLodCoveragePreviewBuild pointsOnly;
    if (!bobol_lod_build_coverage_preview(data, bounds,
	    BOBOL_LOD_DRAW_SHADED, coarse.renderCost - 1u, pointsOnly) ||
	!pointsOnly.pointFallback || pointsOnly.cellAxis != 0 ||
	!pointsOnly.geometry.points || pointsOnly.geometry.shaded ||
	pointsOnly.renderCost > coarse.renderCost - 1u) {
	std::fprintf(stderr, "FAIL: cold coverage did not fall back to points\n");
	return 1;
    }
    return 0;
}

static int
test_degenerate_extent_and_hard_budget(void)
{
    std::vector<fastf_t> planar;
    for (size_t y = 0; y < 24u; ++y) {
	for (size_t x = 0; x < 24u; ++x) {
	    planar.push_back(static_cast<fastf_t>(x) + 0.5);
	    planar.push_back(static_cast<fastf_t>(y) + 0.5);
	    planar.push_back(0.0);
	}
    }
    const struct BObolMeshLodData data = coverage_data(planar, 24.0, 24.0, 0.0);
    const SbBox3f bounds(SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(24.0f, 24.0f, 0.0f));

    BObolLodCoveragePreviewBuild preview;
    if (!bobol_lod_build_coverage_preview(data, bounds,
	    BOBOL_LOD_DRAW_SHADED, SIZE_MAX, preview) ||
	!preview.pointFallback || !preview.geometry.points ||
	preview.geometry.points->positions.empty() || preview.counts.faceCount ||
	preview.counts.pointCount != preview.geometry.points->positions.size()) {
	std::fprintf(stderr, "FAIL: planar cold source has no point coverage\n");
	return 1;
    }
    if (!Obol::cadAdmitPartGeometry(std::move(preview.geometry))) {
	std::fprintf(stderr, "FAIL: planar point coverage geometry was invalid\n");
	return 1;
    }

    BObolLodCounts onePoint;
    onePoint.pointCount = 1;
    onePoint.originalPointCount = 1;
    const size_t minimumCost = bobol_lod_render_cost_units(
	onePoint, BOBOL_LOD_DRAW_SHADED, 1);
    BObolLodCoveragePreviewBuild denied;
    if (!minimumCost || bobol_lod_build_coverage_preview(data, bounds,
	    BOBOL_LOD_DRAW_SHADED, minimumCost - 1u, denied)) {
	std::fprintf(stderr, "FAIL: cold coverage exceeded a hard render budget\n");
	return 1;
    }
    return 0;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argc > 0 ? argv[0] : "test_bobol_lod_coverage_preview");
    return test_adaptive_voxels() || test_degenerate_extent_and_hard_budget();
}
