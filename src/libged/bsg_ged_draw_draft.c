/*                B S G _ G E D _ D R A W _ D R A F T . C
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
/** @file libged/b_s_g___g_e_d___d_r_a_w___d_r_a_f_t_._c.c
 *
 * GED draw-shape draft construction and commit helpers.
 */

#include "common.h"

#include <stdlib.h>
#include <string.h>

#include "bu/ptbl.h"
#include "bu/str.h"
#include "bu/color.h"
#include "bu/hash.h"
#include "bu/parallel.h"
#include "bg/plot3.h"
#include "bg/clip.h"

#include "ged.h"
#include "ged/draw.h"
#include "rt/resource.h"
#include "rt/tree.h"
#include "rt/view.h"
#include "./ged_private.h"
#include "./bsg_ged_draw_private.h"

struct ged_draw_shape_draft {
    struct ged *gedp;
    void *view_ctx;
    bsg_scene_ref source_ref;
    bsg_scene_ref shape_ref;
    int committed;
};


static void
_ged_draw_shape_draft_destroy_refs(bsg_scene_ref source_ref, bsg_scene_ref shape_ref)
{
    if (!ged_draw_scene_ref_is_null(source_ref) &&
	    !ged_draw_scene_ref_is_null(shape_ref))
	(void)ged_draw_scene_ref_detach(shape_ref);
    if (!ged_draw_scene_ref_is_null(shape_ref))
	ged_draw_scene_ref_release(shape_ref);
    if (!ged_draw_scene_ref_is_null(source_ref))
	ged_draw_scene_ref_release(source_ref);
}


static void
_ged_draw_shape_draft_sync_aux_geometry(ged_draw_shape_draft *draft)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->source_ref) ||
	    ged_draw_scene_ref_is_null(draft->shape_ref))
	return;

    (void)ged_draw_scene_ref_sync_aux_display_state(draft->source_ref,
	    draft->shape_ref);
}


ged_draw_shape_draft *
ged_draw_shape_draft_create_context(struct ged *gedp, void *view_ctx, int registered)
{
    bsg_scene_ref source_ref = ged_draw_scene_ref_null();
    bsg_scene_ref shape_ref = ged_draw_scene_ref_null();
    if (!ged_draw_scene_ref_create_draft_pair(gedp, view_ctx, &source_ref,
	    &shape_ref))
	return NULL;

    ged_draw_shape_draft *draft;
    BU_GET(draft, ged_draw_shape_draft);
    draft->gedp = gedp;
    draft->view_ctx = view_ctx;
    draft->source_ref = source_ref;
    draft->shape_ref = shape_ref;
    draft->committed = 0;
    if (registered && !ged_draw_shape_draft_apply_db_object_marker(shape_ref)) {
	ged_draw_shape_draft_destroy(draft);
	return NULL;
    }
    return draft;
}


void
ged_draw_shape_draft_destroy(ged_draw_shape_draft *draft)
{
    if (!draft)
	return;
    if (!draft->committed)
	_ged_draw_shape_draft_destroy_refs(draft->source_ref, draft->shape_ref);
    draft->source_ref = ged_draw_scene_ref_null();
    draft->shape_ref = ged_draw_scene_ref_null();
    BU_PUT(draft, ged_draw_shape_draft);
}


static int
_ged_draw_shape_draft_apply_path_state(ged_draw_shape_draft *draft,
				       const struct db_full_path *path)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_apply_path_state(draft->gedp, draft->shape_ref,
	    path);
}


int
ged_draw_shape_draft_publish_line_set(ged_draw_shape_draft *draft,
				      const point_t *points,
				      const int *commands,
				      size_t point_count)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_publish_line_set(draft->shape_ref, points,
	    commands, point_count);
}


int
ged_draw_shape_draft_publish_bot_wireframe_line_set(ged_draw_shape_draft *draft,
						    const struct rt_bot_internal *bot)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !bot)
	return 0;
    return ged_draw_scene_ref_publish_bot_wireframe_line_set(draft->shape_ref,
	    bot);
}


int
ged_draw_shape_draft_publish_primitive_face_set(ged_draw_shape_draft *draft,
						struct rt_db_internal *ip,
						const struct bg_tess_tol *ttol,
						const struct bn_tol *tol,
						const struct rt_view_info *view_info)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !ip)
	return 0;
    return ged_draw_scene_ref_publish_primitive_face_set(draft->shape_ref, ip,
	    ttol, tol, view_info);
}


int
ged_draw_shape_draft_publish_brep_wireframe_line_set(ged_draw_shape_draft *draft,
						     const struct rt_brep_internal *bi,
						     const struct bn_tol *tol)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !bi)
	return 0;
    return ged_draw_scene_ref_publish_brep_wireframe_line_set(draft->shape_ref,
	    bi, tol);
}


int
ged_draw_shape_draft_publish_poly_wireframe_line_set(ged_draw_shape_draft *draft,
						     const struct rt_pg_internal *pg)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !pg)
	return 0;
    return ged_draw_scene_ref_publish_poly_wireframe_line_set(draft->shape_ref,
	    pg);
}


int
ged_draw_shape_draft_publish_nmg_region(ged_draw_shape_draft *draft,
					const struct nmgregion *r,
					int style)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !r)
	return 0;
    return ged_draw_scene_ref_geometry_publish_nmg_region(draft->shape_ref, r, style);
}


static int
_ged_draw_shape_draft_refresh_scene_bounds(ged_draw_shape_draft *draft)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_update_bounds_context(draft->shape_ref,
	    draft->view_ctx);
}


int
ged_draw_shape_draft_apply_known_bounds(ged_draw_shape_draft *draft,
					const point_t min,
					const point_t max)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    if (!min || !max)
	return 0;

    vect_t center;
    center[X] = (min[X] + max[X]) * 0.5;
    center[Y] = (min[Y] + max[Y]) * 0.5;
    center[Z] = (min[Z] + max[Z]) * 0.5;
    ged_draw_scene_ref_set_draw_center(draft->shape_ref, center);

    fastf_t size = max[X] - min[X];
    V_MAX(size, max[Y] - min[Y]);
    V_MAX(size, max[Z] - min[Z]);
    if (!ged_draw_scene_ref_set_draw_size(draft->shape_ref, size))
	return 0;

    return ged_draw_scene_ref_set_bounds(draft->shape_ref, min, max, 1);
}


int
ged_draw_shape_draft_publish_primitive_wireframe(ged_draw_shape_draft *draft,
						 struct rt_db_internal *ip,
						 const struct bg_tess_tol *ttol,
						 const struct bn_tol *tol,
						 void *view_ctx,
						 int adaptive)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return -1;
    struct rt_view_info view_info;
    const struct rt_view_info *view_infop = NULL;
    if (adaptive && view_ctx) {
	rt_view_context_info_get(&view_info, view_ctx);
	view_infop = &view_info;
    }
    int status = ged_draw_scene_ref_publish_primitive_wireframe(draft->shape_ref, ip,
	    ttol, tol, view_ctx, view_infop, adaptive);
    if (status >= 0 && !adaptive)
	(void)ged_draw_scene_ref_update_bounds_from_geometry(draft->shape_ref,
		NULL);
    return status;
}


int
ged_draw_shape_draft_apply_tree_result_state(ged_draw_shape_draft *draft,
					     int dashed,
					     int has_region,
					     int region_id,
					     int aircode,
					     int los,
					     int material_id,
					     const struct ged_draw_appearance_settings *settings)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !settings)
	return 0;
    (void)_ged_draw_shape_draft_refresh_scene_bounds(draft);
    if (!ged_draw_scene_ref_set_highlighted(draft->shape_ref, 0))
	return 0;
    if (!ged_draw_scene_ref_set_line_style(draft->shape_ref, dashed))
	return 0;
    if (!ged_draw_scene_ref_set_legacy_eval_flag(draft->shape_ref, 0))
	return 0;
    if (!has_region)
	return ged_draw_scene_ref_apply_display_settings(draft->shape_ref,
		settings);
    if (!ged_draw_scene_ref_set_legacy_region_id(draft->shape_ref,
	    region_id))
	return 0;
    if (!ged_draw_shape_draft_apply_registry_region(draft->gedp,
	    draft->shape_ref, region_id, aircode, los, material_id))
	return 0;
    return ged_draw_scene_ref_apply_display_settings(draft->shape_ref,
	    settings);
}


static int
_ged_draw_shape_draft_apply_name(ged_draw_shape_draft *draft, const char *name)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !name)
	return 0;
    if (!ged_draw_scene_ref_set_name(draft->shape_ref, name))
	return 0;
    if (!ged_draw_scene_ref_is_null(draft->source_ref)) {
	struct bu_vls source_name = BU_VLS_INIT_ZERO;
	bu_vls_sprintf(&source_name, "%s:source", name);
	(void)ged_draw_scene_ref_set_name(draft->source_ref,
		bu_vls_cstr(&source_name));
	bu_vls_free(&source_name);
    }
    return 1;
}


int
ged_draw_shape_draft_apply_tree_legacy_color(
	ged_draw_shape_draft *draft,
	const unsigned char wireframe_color_override[3],
	const struct db_tree_state *tsp)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;

    unsigned char bcolor[3] = {255, 0, 0};
    int user_color = 0;
    int default_color = 0;

    if (wireframe_color_override) {
	user_color = 1;
	bcolor[RED] = wireframe_color_override[RED];
	bcolor[GRN] = wireframe_color_override[GRN];
	bcolor[BLU] = wireframe_color_override[BLU];
    } else if (tsp) {
	if (tsp->ts_mater.ma_color_valid) {
	    bcolor[RED] = tsp->ts_mater.ma_color[RED] * 255.0;
	    bcolor[GRN] = tsp->ts_mater.ma_color[GRN] * 255.0;
	    bcolor[BLU] = tsp->ts_mater.ma_color[BLU] * 255.0;
	} else {
	    default_color = 1;
	}
    }

    if (!ged_draw_scene_ref_set_legacy_color_info(draft->shape_ref,
	    bcolor, user_color, default_color))
	return 0;
    color_soltab_scene_ref(tsp ? tsp->ts_dbip : NULL, draft->shape_ref);
    return 1;
}


static int
_ged_draw_shape_draft_attach_source_inputs(
	ged_draw_shape_draft *draft,
	struct db_i *dbip,
	const struct db_full_path *pathp,
	const struct bn_tol *tol,
	const struct bg_tess_tol *ttol)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;

    struct ged_draw_source_state_record record;
    memset(&record, 0, sizeof(record));
    record.dbip = dbip;
    record.fullpath = pathp;
    record.tol = tol;
    record.ttol = ttol;
    record.stale_reason = GED_DRAW_STALE_NONE;
    return ged_draw_scene_ref_attach_source_state_record(draft->shape_ref,
	    &record);
}


int
ged_draw_shape_draft_apply_path_source_state(
	ged_draw_shape_draft *draft,
	struct db_i *dbip,
	const struct db_full_path *pathp,
	const struct bn_tol *tol,
	const struct bg_tess_tol *ttol,
	int has_draw_mat,
	const mat_t draw_mat,
	const char *display_name)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !pathp)
	return 0;
    if (display_name && !_ged_draw_shape_draft_apply_name(draft, display_name))
	return 0;
    if (!_ged_draw_shape_draft_apply_path_state(draft, pathp))
	return 0;
    if (!_ged_draw_shape_draft_attach_source_inputs(draft, dbip, pathp, tol,
	    ttol))
	return 0;
    if (has_draw_mat &&
	    (!draw_mat || !ged_draw_scene_ref_set_draw_mat(draft->shape_ref,
		draw_mat)))
	return 0;
    return 1;
}


static int
_ged_draw_shape_draft_apply_material_rgb(ged_draw_shape_draft *draft,
					 const unsigned char rgb[3])
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !rgb)
	return 0;
    if (!ged_draw_scene_ref_set_color(draft->shape_ref, rgb))
	return 0;
    ged_draw_scene_ref_set_material_rgb(draft->shape_ref, rgb);
    return 1;
}


static int
_ged_draw_shape_draft_apply_material_rgb_revision(ged_draw_shape_draft *draft,
						  const unsigned char rgb[3])
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !rgb)
	return 0;
    if (!_ged_draw_shape_draft_apply_material_rgb(draft, rgb))
	return 0;
    return ged_draw_scene_ref_bump_appearance_revision(draft->shape_ref);
}


int
ged_draw_shape_draft_apply_late_display_state(
	ged_draw_shape_draft *draft,
	const struct db_full_path *path,
	int line_style,
	const struct ged_draw_appearance_settings *settings,
	const unsigned char material_rgb[3],
	int highlighted)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) ||
	    !path || !settings || !material_rgb)
	return 0;
    if (!_ged_draw_shape_draft_refresh_scene_bounds(draft))
	return 0;
    if (!_ged_draw_shape_draft_apply_path_state(draft, path))
	return 0;
    if (!ged_draw_scene_ref_set_line_style(draft->shape_ref, line_style))
	return 0;
    if (!ged_draw_scene_ref_apply_display_settings(draft->shape_ref,
	    settings))
	return 0;
    if (!_ged_draw_shape_draft_apply_material_rgb(draft, material_rgb))
	return 0;
    return ged_draw_scene_ref_set_highlighted(draft->shape_ref, highlighted);
}


int
ged_draw_shape_draft_apply_evaluated_path_display(
	ged_draw_shape_draft *draft,
	const struct ged_draw_appearance_settings *settings,
	const unsigned char material_rgb[3],
	int has_transform,
	const mat_t transform,
	int is_subtraction)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) ||
	    !settings || !material_rgb)
	return 0;

    if (!_ged_draw_shape_draft_apply_material_rgb(draft, material_rgb))
	return 0;
    if (has_transform && transform &&
	    !ged_draw_scene_ref_set_transform(draft->shape_ref, transform))
	return 0;
    if (!settings->draw_solid_lines_only &&
	    !ged_draw_scene_ref_set_line_style(draft->shape_ref,
		is_subtraction ? 1 : 0))
	return 0;
    if (settings->draw_non_subtract_only && is_subtraction &&
	    !ged_draw_scene_ref_set_visible(draft->shape_ref, 0))
	return 0;
    return ged_draw_scene_ref_apply_display_settings(draft->shape_ref,
	    settings);
}


int
ged_draw_shape_draft_apply_database_leaf_display(
	ged_draw_shape_draft *draft,
	const struct ged_draw_appearance_settings *settings,
	int bool_op,
	const unsigned char material_rgb[3],
	int has_draw_size,
	fastf_t draw_size)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) ||
	    !material_rgb)
	return 0;
    if (settings &&
	    !ged_draw_scene_ref_apply_display_settings(draft->shape_ref,
		settings))
	return 0;
    if ((!settings || !settings->draw_solid_lines_only) &&
	    !ged_draw_scene_ref_set_line_style(draft->shape_ref,
		(bool_op == 4) ? 1 : 0))
	return 0;
    if (!_ged_draw_shape_draft_apply_material_rgb_revision(draft,
	    material_rgb))
	return 0;
    if (has_draw_size &&
	    !ged_draw_scene_ref_set_draw_size(draft->shape_ref, draw_size))
	return 0;
    return 1;
}


ged_draw_shape_ref
ged_draw_shape_draft_commit_to_group(ged_draw_shape_draft *draft,
				     ged_draw_group_ref group_ref)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return GED_DRAW_SHAPE_REF_NULL;

    _ged_draw_shape_draft_sync_aux_geometry(draft);

    if (!ged_draw_group_ref_append_scene_ref(draft->gedp, group_ref, draft->shape_ref)) {
	ged_draw_shape_draft_destroy(draft);
	return GED_DRAW_SHAPE_REF_NULL;
    }
    ged_draw_shape_ref ref =
	ged_draw_shape_ref_from_scene_ref(draft->gedp, draft->shape_ref);
    if (ged_draw_shape_ref_is_null(ref)) {
	ged_draw_shape_draft_destroy(draft);
	return GED_DRAW_SHAPE_REF_NULL;
    }
    ref.scene_revision = 0;

    draft->committed = 1;
    ged_draw_shape_draft_destroy(draft);
    return ref;
}


int
ged_draw_shape_draft_commit_database_leaf_to_scene_ref(
	ged_draw_shape_draft *draft,
	bsg_scene_ref parent_ref)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) ||
	    ged_draw_scene_ref_is_null(parent_ref))
	return 0;

    _ged_draw_shape_draft_sync_aux_geometry(draft);

    if (!ged_draw_scene_ref_append_child(parent_ref, draft->source_ref)) {
	ged_draw_shape_draft_destroy(draft);
	return 0;
    }
    draft->committed = 1;
    ged_draw_shape_draft_destroy(draft);
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
