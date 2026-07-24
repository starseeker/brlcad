/*         G E D _ D R A W _ S H A P E _ S T A T E _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file ged_draw_shape_state_private.h
 *
 * Registry/source-local draw-shape state layout.
 */

#ifndef LIBGED_GED_DRAW_SHAPE_STATE_PRIVATE_H
#define LIBGED_GED_DRAW_SHAPE_STATE_PRIVATE_H

#include "common.h"

#include <stdint.h>

#include "bg/defines.h"
#include "bn/tol.h"
#include "ged/draw.h"
#include "rt/db_fullpath.h"
#include "rt/view.h"

struct directory;
struct ged;

typedef struct ged_draw_shape_state {
    struct db_full_path s_fullpath;
    struct directory *leaf_dp;
    char *display_name;
    unsigned long long path_hash;
    int region_id;
    int aircode;
    int los;
    int material_id;
    uint64_t source_revision;
    uint64_t inputs_revision;
    uint64_t realized_source_revision;
    uint64_t realized_inputs_revision;
    ged_draw_stale_reason stale_reason;
    struct ged *gedp;
    void *u_data;
    void (*u_data_free)(void *);
    ged_draw_scene_handle source_ref;
    void *source_data;
    void (*source_data_free)(void *);
    struct bn_tol obol_snapshot_tol;
    struct bg_tess_tol obol_snapshot_ttol;
    int obol_snapshot_valid;
    size_t geometry_command_count;
    uint64_t geometry_revision;
} ged_draw_shape_state;


#endif /* LIBGED_GED_DRAW_SHAPE_STATE_PRIVATE_H */
