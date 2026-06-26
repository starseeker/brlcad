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
#include "bg/plot3.h"
#include "bg/clip.h"

#include "ged.h"
#include "ged/draw.h"
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

    size_t child_count = ged_draw_scene_ref_child_count(draft->source_ref);
    for (size_t i = 0; i < child_count; i++) {
	bsg_scene_ref child_ref = ged_draw_scene_ref_child_at(draft->source_ref, i);
	if (ged_draw_scene_ref_equal(child_ref, draft->shape_ref))
	    continue;
	ged_draw_scene_ref_copy_aux_display_state(child_ref,
		draft->shape_ref);
    }
}


ged_draw_shape_draft *
ged_draw_shape_draft_create_context(struct ged *gedp, void *view_ctx, int registered)
{
    (void)registered;

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


int
ged_draw_shape_draft_set_fullpath(ged_draw_shape_draft *draft,
				  const struct db_full_path *path)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_set_fullpath(draft->gedp, draft->shape_ref, path);
}


int
ged_draw_shape_draft_set_region(ged_draw_shape_draft *draft,
				int region_id,
				int aircode,
				int los,
				int material_id)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_set_region(draft->gedp, draft->shape_ref,
	    region_id, aircode, los, material_id);
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


int
ged_draw_shape_draft_clear_geometry(ged_draw_shape_draft *draft)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_geometry_clear(draft->shape_ref);
}


int
ged_draw_shape_draft_update_scene_bounds(ged_draw_shape_draft *draft)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_update_bounds_context(draft->shape_ref,
	    draft->view_ctx);
}


int
ged_draw_shape_draft_update_bounds_from_geometry(ged_draw_shape_draft *draft,
						 int *bad_cmd)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_update_bounds_from_geometry(draft->shape_ref, bad_cmd);
}


int
ged_draw_shape_draft_set_bounds_from_minmax(ged_draw_shape_draft *draft,
					    const point_t min,
					    const point_t max,
					    int set_node_bounds)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_set_bounds_from_minmax(draft->shape_ref, min, max,
	    set_node_bounds);
}


int
ged_draw_shape_draft_set_center_size(ged_draw_shape_draft *draft,
				     const point_t center,
				     fastf_t size)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !center)
	return 0;
    ged_draw_scene_ref_set_draw_center(draft->shape_ref, center);
    return ged_draw_scene_ref_set_draw_size(draft->shape_ref, size);
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
    return ged_draw_scene_ref_publish_primitive_wireframe(draft->shape_ref, ip,
	    ttol, tol, view_ctx, view_infop, adaptive);
}


int
ged_draw_shape_draft_set_highlighted(ged_draw_shape_draft *draft,
				     int highlighted)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_set_highlighted(draft->shape_ref, highlighted);
}


int
ged_draw_shape_draft_set_line_style(ged_draw_shape_draft *draft, int dashed)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_set_line_style(draft->shape_ref, dashed);
}


int
ged_draw_shape_draft_set_line_width(ged_draw_shape_draft *draft, int line_width)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_set_line_width(draft->shape_ref, line_width);
}


int
ged_draw_shape_draft_set_legacy_color_info(ged_draw_shape_draft *draft,
					   const unsigned char basecolor[3],
					   int user_color,
					   int default_color)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !basecolor)
	return 0;
    return ged_draw_scene_ref_set_legacy_color_info(draft->shape_ref,
	    basecolor, user_color, default_color);
}


int
ged_draw_shape_draft_set_legacy_uses_default_color(ged_draw_shape_draft *draft,
						   int default_color)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_set_legacy_uses_default_color(draft->shape_ref,
	    default_color);
}


int
ged_draw_shape_draft_set_legacy_eval_flag(ged_draw_shape_draft *draft,
					  int eflag)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_set_legacy_eval_flag(draft->shape_ref, eflag);
}


int
ged_draw_shape_draft_set_legacy_region_id(ged_draw_shape_draft *draft,
					  int region_id)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_set_legacy_region_id(draft->shape_ref,
	    region_id);
}


int
ged_draw_shape_draft_color_soltab(ged_draw_shape_draft *draft, struct db_i *dbip)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    color_soltab_scene_ref(dbip, draft->shape_ref);
    return 1;
}


int
ged_draw_shape_draft_set_name(ged_draw_shape_draft *draft, const char *name)
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
ged_draw_shape_draft_set_source_state(ged_draw_shape_draft *draft,
				      struct ged_draw_source_state *data)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    ged_draw_scene_ref_source_data_set(draft->shape_ref, data);
    return 1;
}


int
ged_draw_shape_draft_mark_db_object(ged_draw_shape_draft *draft)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_mark_db_object(draft->shape_ref);
}


int
ged_draw_shape_draft_apply_appearance_settings(ged_draw_shape_draft *draft,
					       const struct ged_draw_appearance_settings *settings)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !settings)
	return 0;

    return ged_draw_scene_ref_apply_appearance_settings(draft->shape_ref,
	    settings);
}


int
ged_draw_shape_draft_bump_appearance_revision(ged_draw_shape_draft *draft)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_bump_changed(draft->shape_ref);
}


int
ged_draw_shape_draft_set_material_rgb(ged_draw_shape_draft *draft,
				      unsigned char r,
				      unsigned char g,
				      unsigned char b)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    const unsigned char rgb[3] = {r, g, b};
    if (!ged_draw_scene_ref_set_color(draft->shape_ref, rgb))
	return 0;
    ged_draw_scene_ref_set_material_rgb(draft->shape_ref, rgb);
    return 1;
}


int
ged_draw_shape_draft_set_transform(ged_draw_shape_draft *draft, const mat_t mat)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !mat)
	return 0;
    return ged_draw_scene_ref_set_transform(draft->shape_ref, mat);
}


int
ged_draw_shape_draft_set_bounds(ged_draw_shape_draft *draft,
				const point_t min,
				const point_t max)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !min || !max)
	return 0;
    return ged_draw_scene_ref_set_bounds(draft->shape_ref, min, max, 1);
}


int
ged_draw_shape_draft_set_draw_size(ged_draw_shape_draft *draft, fastf_t size)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_set_draw_size(draft->shape_ref, size);
}


int
ged_draw_shape_draft_set_visible(ged_draw_shape_draft *draft, int visible)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_set_visible(draft->shape_ref, visible);
}


int
ged_draw_shape_draft_set_transparency(ged_draw_shape_draft *draft,
				      fastf_t transparency)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_set_transparency(draft->shape_ref,
	    transparency);
}


int
ged_draw_shape_draft_set_draw_mode(ged_draw_shape_draft *draft, int draw_mode)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    ged_draw_scene_ref_set_draw_mode(draft->shape_ref, draw_mode);
    return 1;
}


int
ged_draw_shape_draft_set_draw_mat(ged_draw_shape_draft *draft, const mat_t mat)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !mat)
	return 0;
    return ged_draw_scene_ref_set_draw_mat(draft->shape_ref, mat);
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


bsg_scene_ref
ged_draw_shape_draft_commit_to_scene_ref(ged_draw_shape_draft *draft,
					 bsg_scene_ref parent_ref)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) ||
	    ged_draw_scene_ref_is_null(parent_ref))
	return ged_draw_scene_ref_null();

    _ged_draw_shape_draft_sync_aux_geometry(draft);

    bsg_scene_ref shape_ref = draft->shape_ref;
    if (!ged_draw_scene_ref_append_child(parent_ref, draft->source_ref)) {
	ged_draw_shape_draft_destroy(draft);
	return ged_draw_scene_ref_null();
    }
    draft->committed = 1;
    ged_draw_shape_draft_destroy(draft);
    return shape_ref;
}


ged_draw_shape_ref
ged_draw_shape_draft_commit_to_last_group(ged_draw_shape_draft *draft)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return GED_DRAW_SHAPE_REF_NULL;

    _ged_draw_shape_draft_sync_aux_geometry(draft);

    ged_draw_append_scene_ref_to_last_group(draft->gedp, draft->shape_ref);
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


/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
