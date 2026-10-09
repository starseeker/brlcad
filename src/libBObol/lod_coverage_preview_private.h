/*         L O D _ C O V E R A G E _ P R E V I E W _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_COVERAGE_PREVIEW_PRIVATE_H
#define LIBBOBOL_LOD_COVERAGE_PREVIEW_PRIVATE_H

#include "common.h"

#include "BObol/BLodRealization.h"
#include "BObol/BMeshLodCache.h"

#include <Obol/cad/CadGeometry.h>

#include <cstddef>

/* One budget-admitted, presentation-only whole-source summary.  A nonzero
 * cellAxis identifies voxel geometry; pointFallback identifies the bounded
 * point representation used for a degenerate source or when no voxel grid
 * fits. */
struct BObolLodCoveragePreviewBuild {
    Obol::PartGeometryBuilder geometry;
    BObolLodCounts counts;
    size_t renderCost = 0;
    size_t cellAxis = 0;
    bool pointFallback = false;
};

/* Select the richest 24/12/6 occupancy grid which fits renderCostAllowance.
 * If no grid fits, or the source extent is planar/linear, retain a spatially
 * stratified point summary when at least one point can fit.  The input arrays
 * are borrowed and output owns all renderer data. */
bool bobol_lod_build_coverage_preview(
    const struct BObolMeshLodData &data,
    const SbBox3f &conservativeBounds,
    int drawMode,
    size_t renderCostAllowance,
    BObolLodCoveragePreviewBuild &result);

#endif /* LIBBOBOL_LOD_COVERAGE_PREVIEW_PRIVATE_H */
