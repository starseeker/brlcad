/*              B S G _ G E D _ D R A W _ O V E R L A Y . C
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
/** @file bsg_ged_draw_overlay.c
 *
 * Typed GED overlay insertion helpers.
 */

#include "common.h"

#include "bu/log.h"

#include "ged.h"
#include "ged/draw.h"
#include "./ged_private.h"
#include "./bsg_ged_draw_private.h"


void
ged_draw_overlay_erase_name(struct ged *gedp, const char *name)
{
    void *view_ctx = ged_draw_active_view_ctx(gedp);
    bsg_scene_ref ref = ged_draw_view_context_overlay_scene_find(view_ctx, name);
    if (!ged_draw_scene_ref_is_null(ref))
	ged_draw_scene_ref_highlight_free_cb(ref);
    ged_draw_view_context_overlay_name_erase(view_ctx, name);
}


static int
_ged_overlay_publish_geometry(bsg_scene_ref ref,
			      const struct ged_draw_overlay_geometry *geometry)
{
    if (ged_draw_scene_ref_is_null(ref) || !geometry)
	return 0;

    switch (geometry->kind) {
	case GED_DRAW_OVERLAY_GEOMETRY_LINE_SET:
	    if (geometry->point_count != geometry->command_count)
		return 0;
	    return ged_draw_scene_ref_publish_line_set(ref,
		    geometry->points, geometry->commands, geometry->point_count);
	case GED_DRAW_OVERLAY_GEOMETRY_POINT_SET:
	    return ged_draw_scene_ref_publish_point_set(ref,
		    geometry->points, geometry->point_count);
	case GED_DRAW_OVERLAY_GEOMETRY_INDEXED_FACE_SET:
	    return ged_draw_scene_ref_publish_indexed_face_set(ref,
		    geometry->points, geometry->point_count,
		    geometry->normals, geometry->normal_count,
		    geometry->indices, geometry->index_count);
	default:
	    return 0;
    }
}


int
ged_draw_overlay_geometry_insert(struct ged *gedp, const char *name,
	const struct ged_draw_overlay_geometry *geometry,
	long int rgb, fastf_t transparency, int draw_mode, int csoltab,
	ged_draw_shape_ref *out)
{
    if (out)
	*out = GED_DRAW_SHAPE_REF_NULL;
    if (!gedp || !name || !geometry)
	return -1;
    void *view_ctx = ged_draw_active_view_ctx(gedp);
    if (!view_ctx)
	return 0;

    struct db_i *dbip = gedp->dbip;
    if (dbip == DBI_NULL)
	return 0;

    if (db_lookup(dbip, name, LOOKUP_QUIET) != RT_DIR_NULL) {
	bu_log("ged_draw_overlay_geometry_insert(%s) would clobber existing database entry, "
		"ignored\n", name);
	return -1;
    }

    ged_draw_overlay_erase_name(gedp, name);

    bsg_scene_ref overlay_scene = ged_draw_view_context_overlay_geometry_create(view_ctx, name,
	    geometry->kind);
    if (ged_draw_scene_ref_is_null(overlay_scene))
	return -1;
    (void)ged_draw_view_overlay_command_result_owner_set(overlay_scene,
	    gedp, name);
    (void)ged_draw_scene_ref_set_name(overlay_scene, name);

    if (!ged_draw_scene_ref_ensure_registry_entry(gedp, overlay_scene)) {
	ged_draw_scene_ref_release(overlay_scene);
	return -1;
    }

    if (!_ged_overlay_publish_geometry(overlay_scene, geometry)) {
	ged_draw_scene_ref_release(overlay_scene);
	return -1;
    }
    ged_draw_scene_ref_update_bounds_context(overlay_scene, view_ctx);

    if (!ged_draw_view_context_overlay_scene_append(view_ctx, overlay_scene)) {
	ged_draw_scene_ref_release(overlay_scene);
	return -1;
    }

    unsigned char basecolor[3] = {
	(unsigned char)((rgb >> 16) & 0xFF),
	(unsigned char)((rgb >>  8) & 0xFF),
	(unsigned char)(rgb & 0xFF)
    };
    (void)ged_draw_scene_ref_apply_overlay_geometry_attributes(overlay_scene,
	    basecolor, transparency, draw_mode);

    if (csoltab)
	color_soltab_scene_ref(gedp->dbip, overlay_scene);

    if (out)
	*out = ged_draw_shape_ref_from_scene_ref(gedp, overlay_scene);

    return 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
