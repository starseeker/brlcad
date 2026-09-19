/*              D R A W _ C A C H E _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_DRAW_CACHE_PRIVATE_H
#define LIBBOBOL_DRAW_CACHE_PRIVATE_H

#include "common.h"

#include "vmath.h"

#include <stddef.h>
#include <stdint.h>

struct db_i;

void bobol_draw_cache_runtime_prepare(void);

/* Merge producer-certified canonical-asset OBB metadata into the compact
 * draw-asset record.  The ordinary AABB remains the coverage authority. */
int bobol_draw_lod_asset_oriented_bounds_publish(
    struct db_i *dbip, const char *name, uint64_t faceCount,
    uint64_t pointCount, const point_t boundsMin, const point_t boundsMax,
    const point_t orientedBounds[8]);

#endif /* LIBBOBOL_DRAW_CACHE_PRIVATE_H */
