/*                B S G _ G E D _ D R A W _ R E F S . C
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
/** @file libged/b_s_g___g_e_d___d_r_a_w___r_e_f_s_._c.c
 *
 * Draw shape-ref mutation compatibility bridge.
 */

#include "common.h"

#include <stdlib.h>
#include <string.h>

#include "bu/malloc.h"
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

static rt_view_scene_ref
_ged_draw_scene_ref_to_rt(bsg_scene_ref ref)
{
    rt_view_scene_ref rt_ref = RT_VIEW_SCENE_REF_NULL_INIT;
    rt_ref.opaque = ged_draw_scene_ref_context(ref);
    return rt_ref;
}


rt_view_scene_ref
ged_draw_scene_ref_to_rt_view_ref(bsg_scene_ref ref)
{
    return _ged_draw_scene_ref_to_rt(ref);
}


bsg_scene_ref
ged_draw_scene_ref_from_rt_view_ref(rt_view_scene_ref ref)
{
    return ged_draw_scene_ref_from_context(ref.opaque);
}


bsg_scene_ref
ged_draw_scene_ref_active_scope(struct ged *gedp, void *view_ctx, int create,
				int fallback_root)
{
    if (!gedp)
	return ged_draw_scene_ref_null();

    if (view_ctx && rt_view_context_is_independent(view_ctx)) {
	bsg_scene_ref root_ref = ged_draw_scene_ref_null();
	if (create || fallback_root) {
	    root_ref = create ? ged_draw_ensure_root(gedp) :
		ged_draw_scene_ref_from_rt_view_ref(ged_scene_root_rt_ref(gedp));
	    if (create && ged_draw_scene_ref_is_null(root_ref))
		return ged_draw_scene_ref_null();
	}

	bsg_scene_ref scope_ref = ged_draw_scene_ref_from_rt_view_ref(
		rt_view_context_independent_scope_ref(view_ctx, create));
	if (!ged_draw_scene_ref_is_null(scope_ref))
	    return scope_ref;
	return fallback_root ? root_ref : ged_draw_scene_ref_null();
    }

    return create ? ged_draw_ensure_root(gedp) :
	ged_draw_scene_ref_from_rt_view_ref(ged_scene_root_rt_ref(gedp));
}


int
ged_draw_shape_ref_apply_qray_work_flag(struct ged *gedp, ged_draw_shape_ref ref, int wflag)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;
    ged_draw_scene_ref_apply_qray_work_flag(shape_ref, wflag);
    return 1;
}


void *
ged_draw_first_shape_context(struct ged *gedp)
{
    return ged_draw_scene_ref_context(ged_draw_first_shape_scene_ref(gedp));
}


void *
ged_draw_shape_ref_context(struct ged *gedp, ged_draw_shape_ref ref)
{
    return ged_draw_scene_ref_context(
	    ged_draw_registry_shape_scene_ref(gedp, ref));
}


void *
ged_draw_shape_ref_cache_context(struct ged *gedp, ged_draw_shape_ref ref)
{
    return ged_draw_scene_ref_context(
	    ged_draw_shape_scene_ref_from_cache_ref(gedp, ref));
}


void *
ged_draw_group_ref_context(struct ged *gedp, ged_draw_group_ref ref)
{
    return ged_draw_scene_ref_context(
	    ged_draw_registry_group_scene_ref(gedp, ref));
}


ged_draw_shape_ref
ged_draw_shape_ref_from_context(struct ged *gedp, void *shape_ctx)
{
    return ged_draw_shape_ref_from_scene_ref(gedp,
	    ged_draw_scene_ref_from_context(shape_ctx));
}


int
ged_draw_shape_ref_last_point(struct ged *gedp, ged_draw_shape_ref ref, point_t out)
{
    return ged_draw_scene_ref_last_point(
	    ged_draw_registry_shape_scene_ref(gedp, ref), out);
}


int
ged_draw_shape_ref_line_summary(struct ged *gedp,
				ged_draw_shape_ref ref,
				struct ged_draw_view_line_summary *out)
{
    return ged_draw_scene_ref_line_summary(
	    ged_draw_registry_shape_scene_ref(gedp, ref), out);
}


int
ged_draw_shape_ref_line_point_at(struct ged *gedp,
				 ged_draw_shape_ref ref,
				 size_t index,
				 point_t out)
{
    return ged_draw_scene_ref_line_point_at(
	    ged_draw_registry_shape_scene_ref(gedp, ref), index, out);
}


int
ged_draw_shape_ref_line_command_at(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   size_t index,
				   int *out)
{
    return ged_draw_scene_ref_line_command_at(
	    ged_draw_registry_shape_scene_ref(gedp, ref), index, out);
}


int
ged_draw_shape_context_line_summary(void *shape_ctx,
				    struct ged_draw_view_line_summary *out)
{
    return ged_draw_scene_ref_line_summary(
	    ged_draw_scene_ref_from_context(shape_ctx), out);
}


int
ged_draw_shape_context_line_point_at(void *shape_ctx,
				     size_t index,
				     point_t out)
{
    return ged_draw_scene_ref_line_point_at(
	    ged_draw_scene_ref_from_context(shape_ctx), index, out);
}


int
ged_draw_shape_context_line_command_at(void *shape_ctx,
				       size_t index,
				       int *out)
{
    return ged_draw_scene_ref_line_command_at(
	    ged_draw_scene_ref_from_context(shape_ctx), index, out);
}


int
ged_draw_shape_context_geometry_summary(void *shape_ctx,
					struct ged_draw_shape_geometry_summary *out)
{
    return ged_draw_scene_ref_geometry_summary(
	    ged_draw_scene_ref_from_context(shape_ctx), out);
}


int
ged_draw_shape_context_has_state(void *shape_ctx)
{
    if (!shape_ctx)
	return 0;

    struct ged_draw_shape_record_summary summary;
    return ged_draw_scene_ref_shape_record_summary(
	    ged_draw_scene_ref_from_context(shape_ctx), &summary);
}


void *
ged_draw_shape_context_source(void *shape_ctx)
{
    if (!shape_ctx)
	return NULL;

    return ged_draw_scene_ref_context(
	    ged_draw_scene_ref_source_owner(
		ged_draw_scene_ref_from_context(shape_ctx)));
}


const struct db_full_path *
ged_draw_scene_context_fullpath(void *scene_ctx)
{
    if (!scene_ctx)
	return NULL;

    struct ged_draw_scene_tree_summary summary;
    if (!ged_draw_scene_ref_tree_summary(
		ged_draw_scene_ref_from_context(scene_ctx), &summary))
	return NULL;
    return summary.fullpath;
}


int
ged_draw_group_context_dbpath(struct ged *gedp,
			      void *group_ctx,
			      struct db_full_path *out)
{
    (void)gedp;
    if (!group_ctx || !out)
	return -1;

    bsg_scene_ref group_ref = ged_draw_scene_ref_from_context(group_ctx);
    struct ged_draw_group_record_summary group_summary;
    if (ged_draw_scene_ref_group_record_summary(group_ref, &group_summary) &&
	    group_summary.is_overlay)
	return -1;

    const struct db_full_path *stored = ged_draw_scene_context_fullpath(group_ctx);
    if (!stored || stored->fp_len <= 0)
	return -1;
    db_dup_full_path(out, stored);
    return 0;
}


int
ged_draw_group_context_is_overlay(void *group_ctx)
{
    if (!group_ctx)
	return 0;

    struct ged_draw_group_record_summary group_summary;
    if (!ged_draw_scene_ref_group_record_summary(
	    ged_draw_scene_ref_from_context(group_ctx), &group_summary))
	return 0;
    return group_summary.is_overlay;
}


int
ged_draw_scene_context_display_summary(void *scene_ctx,
				       struct ged_draw_scene_display_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!scene_ctx)
	return 0;

    bsg_scene_ref scene_ref = ged_draw_scene_ref_from_context(scene_ctx);
    return ged_draw_scene_ref_display_summary(scene_ref, out);
}


int
ged_draw_scene_context_source_summary(void *scene_ctx,
				      struct ged_draw_database_source_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!scene_ctx)
	return 0;

    bsg_scene_ref scene_ref = ged_draw_scene_ref_from_context(scene_ctx);
    return ged_draw_scene_ref_source_summary(scene_ref, out);
}


int
ged_draw_scene_context_tree_summary(void *scene_ctx,
				    struct ged_draw_scene_tree_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!scene_ctx)
	return 0;

    bsg_scene_ref scene_ref = ged_draw_scene_ref_from_context(scene_ctx);
    return ged_draw_scene_ref_tree_summary(scene_ref, out);
}


struct ged_draw_context_child_at_ctx {
    size_t index;
    size_t cur;
    void *ctx;
};


static int
_ged_draw_context_child_at_cb(bsg_scene_ref child_ref, void *userdata)
{
    struct ged_draw_context_child_at_ctx *ctx =
	(struct ged_draw_context_child_at_ctx *)userdata;
    if (!ctx)
	return 0;
    if (ctx->cur == ctx->index) {
	ctx->ctx = ged_draw_scene_ref_context(child_ref);
	return 0;
    }
    ctx->cur++;
    return 1;
}


void *
ged_draw_scene_context_child_at(void *scene_ctx, size_t index)
{
    if (!scene_ctx)
	return NULL;

    struct ged_draw_context_child_at_ctx child_ctx;
    child_ctx.index = index;
    child_ctx.cur = 0;
    child_ctx.ctx = NULL;

    (void)ged_draw_scene_ref_foreach_child(
	    ged_draw_scene_ref_from_context(scene_ctx),
	    _ged_draw_context_child_at_cb,
	    &child_ctx);
    return child_ctx.ctx;
}


void *
ged_draw_scene_context_parent(void *scene_ctx)
{
    if (!scene_ctx)
	return NULL;

    return ged_draw_scene_ref_parent_context(
	    ged_draw_scene_ref_from_context(scene_ctx));
}


const char *
ged_draw_scene_context_name(void *scene_ctx)
{
    if (!scene_ctx)
	return NULL;

    struct ged_draw_scene_tree_summary summary;
    if (!ged_draw_scene_ref_tree_summary(
		ged_draw_scene_ref_from_context(scene_ctx), &summary))
	return NULL;

    return summary.name;
}


int
ged_draw_scene_context_subtree_bounds(void *scene_ctx,
				      vect_t *min,
				      vect_t *max,
				      int include_overlays)
{
    if (!scene_ctx || !min || !max)
	return 1;

    bsg_scene_ref scene_ref = ged_draw_scene_ref_from_context(scene_ctx);
    if (ged_draw_scene_ref_is_null(scene_ref))
	return 1;

    return ged_draw_scene_ref_subtree_bounds(scene_ref, min, max, include_overlays);
}


int
ged_draw_shape_ref_translate_geometry(struct ged *gedp, ged_draw_shape_ref ref,
				      const vect_t xlate)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    return ged_draw_scene_ref_translate_geometry(shape_ref, xlate);
}


int
ged_draw_shape_ref_set_center(struct ged *gedp, ged_draw_shape_ref ref,
			      const point_t center)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref) || !center)
	return 0;
    ged_draw_scene_ref_set_draw_center(shape_ref, center);
    return 1;
}


int
ged_draw_shape_ref_geometry_clear(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;
    return ged_draw_scene_ref_geometry_clear(shape_ref);
}


int
ged_draw_shape_ref_update_bounds_from_geometry(struct ged *gedp,
					       ged_draw_shape_ref ref,
					       int *bad_cmd)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    return ged_draw_scene_ref_update_bounds_from_geometry(shape_ref, bad_cmd);
}


int
ged_draw_shape_ref_source_snapshot(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   struct ged_draw_shape_source_snapshot *out)
{
    if (!out)
	return 0;
    memset(out, 0, sizeof(*out));

    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;

    return ged_draw_scene_ref_source_snapshot(shape_ref, out);
}


int
ged_draw_shape_ref_publish_line_set(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    const point_t *points,
				    const int *commands,
				    size_t point_count)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;
    return ged_draw_scene_ref_publish_line_set(shape_ref, points, commands,
	    point_count);
}


int
ged_draw_shape_ref_publish_indexed_face_set(struct ged *gedp,
					    ged_draw_shape_ref ref,
					    const point_t *points,
					    size_t point_count,
					    const vect_t *normals,
					    size_t normal_count,
					    const int *indices,
					    size_t index_count)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;
    return ged_draw_scene_ref_publish_indexed_face_set(shape_ref, points,
	    point_count, normals, normal_count, indices, index_count);
}


int
ged_draw_shape_ref_publish_primitive_wireframe(struct ged *gedp,
					       ged_draw_shape_ref ref,
					       struct rt_db_internal *ip,
					       const struct bg_tess_tol *ttol,
					       const struct bn_tol *tol,
					       void *view_ctx,
					       int adaptive)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return -1;
    struct rt_view_info view_info;
    const struct rt_view_info *view_infop = NULL;
    if (adaptive && view_ctx) {
	rt_view_context_info_get(&view_info, view_ctx);
	view_infop = &view_info;
    }
    return ged_draw_scene_ref_publish_primitive_wireframe(shape_ref, ip, ttol,
	    tol, view_ctx, view_infop, adaptive);
}


int
ged_draw_shape_ref_release(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;

    bsg_scene_ref release_ref = ged_draw_scene_ref_source_owner(shape_ref);
    if (ged_draw_scene_ref_is_null(release_ref))
	release_ref = shape_ref;
    (void)ged_draw_scene_ref_detach(release_ref);
    ged_draw_scene_ref_release(release_ref);
    return 1;
}


int
ged_draw_shape_ref_realize_context(struct ged *gedp, ged_draw_shape_ref ref,
				   void *view_ctx)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;
    ged_draw_scene_ref_realize(shape_ref, view_ctx);
    return 1;
}


int
ged_draw_shape_ref_set_visible(struct ged *gedp, ged_draw_shape_ref ref,
			       int visible)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    return ged_draw_scene_ref_set_visible(shape_ref, visible);
}


int
ged_draw_shape_ref_get_color(struct ged *gedp, ged_draw_shape_ref ref,
			     unsigned char rgb[3])
{
    if (!rgb)
	return 0;

    struct ged_draw_shape_material_summary material;
    if (!ged_draw_shape_ref_material_summary(gedp, ref, &material) ||
	    !material.valid)
	return 0;

    rgb[0] = material.material_color[0];
    rgb[1] = material.material_color[1];
    rgb[2] = material.material_color[2];
    return 1;
}


int
ged_draw_shape_ref_set_color(struct ged *gedp, ged_draw_shape_ref ref,
			     const unsigned char rgb[3])
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    return ged_draw_scene_ref_set_color(shape_ref, rgb);
}


int
ged_draw_shape_ref_material_summary(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    struct ged_draw_shape_material_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    return ged_draw_scene_ref_material_summary(shape_ref, out);
}


int
ged_draw_shape_ref_set_material_color(struct ged *gedp,
				      ged_draw_shape_ref ref,
				      const unsigned char rgb[3])
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref) || !rgb)
	return 0;
    ged_draw_scene_ref_set_material_rgb(shape_ref, rgb);
    return 1;
}


int
ged_draw_shape_ref_set_evaluated_region(struct ged *gedp,
					ged_draw_shape_ref ref,
					int evaluated_region)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    return ged_draw_scene_ref_set_legacy_eval_flag(shape_ref,
	    evaluated_region);
}


int
ged_draw_shape_ref_lod_ensure(struct ged *gedp, ged_draw_shape_ref ref,
			      void *first_view_ctx,
			      void **view_ctxs,
			      size_t view_ctx_count)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref) || !first_view_ctx)
	return 0;

    ged_draw_view_lod_policy policy;
    rt_view_context_lod_policy_get(&policy, first_view_ctx);
    struct ged_draw_scene_draw_state_summary draw_state = {0};
    if (!ged_draw_scene_ref_draw_state_summary(shape_ref, &draw_state))
	return 1;
    int candidate = draw_state.pipeline_candidate;
    struct ged_draw_database_source_summary source_summary;
    int source_backed =
	ged_draw_scene_ref_source_summary(
		ged_draw_scene_ref_source_owner(shape_ref), &source_summary) &&
	source_summary.has_state;
    if (!candidate && !source_backed)
	return 1;
    if (!candidate) {
	struct ged_draw_shape_source_snapshot source = {0};
	struct directory *dp = NULL;
	if (ged_draw_scene_ref_source_snapshot(shape_ref, &source))
	    dp = source.leaf_dp;
	int mode = draw_state.draw_mode;
	int mesh_ref = dp && policy.mesh_enabled &&
	    ((dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT &&
	      (mode == 0 || mode == 1)) ||
	     (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP &&
	      mode == 1));
	if (!mesh_ref && !(policy.csg_enabled && mode == 0))
	    return 1;
    }
    int realized = 0;
    for (size_t i = 0; i < view_ctx_count; i++) {
	if (view_ctxs && view_ctxs[i]) {
	    (void)ged_draw_scene_ref_realize_view_inputs_changed(shape_ref,
		    view_ctxs[i]);
	    realized = 1;
	}
    }
    if (!realized) {
	(void)ged_draw_scene_ref_realize_view_inputs_changed(shape_ref,
		first_view_ctx);
    }
    return 1;
}


void *
ged_draw_shape_ref_view_context(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    struct ged_draw_scene_draw_state_summary draw_state = {0};
    if (!ged_draw_scene_ref_draw_state_summary(shape_ref, &draw_state))
	return NULL;
    return draw_state.view_ctx;
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
