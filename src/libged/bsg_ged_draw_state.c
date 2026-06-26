/*                B S G _ G E D _ D R A W _ S T A T E . C
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
/** @file libged/b_s_g___g_e_d___d_r_a_w___s_t_a_t_e_._c.c
 *
 * GED draw-shape path and region state helpers.
 */

#include "common.h"

#include "ged.h"
#include "ged/draw.h"
#include "./ged_private.h"
#include "./bsg_ged_draw_private.h"

int
ged_draw_scene_ref_apply_path_state(struct ged *gedp,
				    bsg_scene_ref ref,
				    const struct db_full_path *path)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 0;

    if (!ged_draw_shape_state_ensure_scene_ref(gedp, ref))
	return 0;
    if (!ged_draw_scene_ref_set_indexed_fullpath(gedp, ref, path))
	return 0;

    if (path && path->fp_len > 0) {
	char *path_name = db_path_to_string(path);
	if (path_name) {
	    const char *semantic_path = path_name;
	    while (*semantic_path == '/')
		semantic_path++;
	    (void)ged_draw_scene_ref_set_draw_intent_path(ref,
		    semantic_path, GED_DRAW_MODE_WIRE);
	    bu_free(path_name, "ged draw shape intent path string");
	}
    }

    return 1;
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
