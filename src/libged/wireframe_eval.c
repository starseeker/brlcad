/*                 W I R E F R A M E _ E V A L . C
 * BRL-CAD
 *
 * Copyright (c) 1997-2026 United States Government as represented by
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this file; see the file named COPYING for more information.
 */
/** @addtogroup libged */
/** @{ */
/** @file libged/draw/wireframe_eval.c
 *
 * GED evaluated-wireframe provider adapter.
 *
 * The CSG/evaluated-wire implementation lives in libbrlobol so Obol-owned
 * database-source realization does not depend on retained BSG shape-ref
 * evaluation internals.  This file only resolves GED runtime/source state and
 * publishes the neutral line-command result.
 */
/** @} */

#include "common.h"

#include "brlobol/evaluated_wire.h"
#include "bg/defines.h"
#include "bn/tol.h"
#include "bu/vls.h"
#include "ged/draw.h"
#include "./bsg_ged_draw_private.h"


int
ged_draw_shape_ref_eval_wireframe(struct ged *gedp, ged_draw_shape_ref ref)
{
    if (!gedp || ged_draw_shape_ref_is_null(ref))
	return BRLCAD_OK;

    struct ged_draw_shape_source_snapshot source;
    if (!ged_draw_shape_ref_source_snapshot(gedp, ref, &source))
	return BRLCAD_OK;
    if (!source.dbip || !source.tol || !source.ttol)
	return BRLCAD_ERROR;

    const char *path = NULL;
    struct bu_vls ppath = BU_VLS_INIT_ZERO;
    if (source.fullpath) {
	db_path_to_vls(&ppath, source.fullpath);
	path = bu_vls_cstr(&ppath);
    } else {
	path = source.name;
    }

    point_t *line_points = NULL;
    int *line_commands = NULL;
    size_t line_count = 0;
    int ret = brlobol_evaluated_wire_evaluate_path_line_set(source.dbip,
	    path, source.tol, source.ttol, &line_points, &line_commands,
	    &line_count);
    bu_vls_free(&ppath);
    if (ret != BRLCAD_OK)
	return ret;

    if (!ged_draw_shape_ref_publish_line_set(gedp, ref,
	    (const point_t *)line_points, line_commands, line_count))
	ret = BRLCAD_ERROR;

    brlobol_evaluated_wire_line_set_free(line_points, line_commands);
    return ret;
}


int
ged_draw_obol_database_source_eval_wireframe_for_path(struct ged *gedp,
						      const char *path)
{
    if (!gedp || !path || !path[0])
	return BRLCAD_ERROR;

    struct ged_draw_obol_database_source_runtime runtime;
    memset(&runtime, 0, sizeof(runtime));
    if (!ged_draw_obol_database_source_runtime_for_path(gedp, path,
	    &runtime) || !runtime.valid || !runtime.dbip)
	return BRLCAD_ERROR;

    struct bn_tol tol;
    BN_TOL_INIT_SET_TOL(&tol);

    struct bg_tess_tol ttol;
    BG_TESS_TOL_INIT_SET_TOL(&ttol);
    if (runtime.tessellation_abs_tol >= 0.0)
	ttol.abs = runtime.tessellation_abs_tol;
    if (runtime.tessellation_rel_tol >= 0.0)
	ttol.rel = runtime.tessellation_rel_tol;
    if (runtime.tessellation_norm_tol >= 0.0)
	ttol.norm = runtime.tessellation_norm_tol;

    point_t *line_points = NULL;
    int *line_commands = NULL;
    size_t line_count = 0;
    int ret = brlobol_evaluated_wire_evaluate_path_line_set(runtime.dbip,
	    path, &tol, &ttol, &line_points, &line_commands, &line_count);
    if (ret == BRLCAD_OK &&
	    !ged_draw_obol_database_source_publish_line_set_for_path(gedp,
		path, (const point_t *)line_points, line_commands, line_count))
	ret = BRLCAD_ERROR;

    brlobol_evaluated_wire_line_set_free(line_points, line_commands);
    return ret;
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
