/*                B S G _ G E D _ D R A W _ V I E W . C
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this file; see the file named COPYING for more information.
 */
/** @file bsg_ged_draw_view.c
 *
 * BSG-to-RT view snapshot adapters for the transitional GED draw bridge.
 */

#include "common.h"

#include <string.h>

#include "bg/polygon.h"
#include "bg/line_layer.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "bsg/feature.h"
#include "bsg/geometry.h"
#include "bsg/hud.h"
#include "bsg/interaction.h"
#include "bsg/overlay.h"
#include "bsg/polygon.h"
#include "bsg/selection.h"
#include "bsg/snap_action.h"
#include "bsg/view_set.h"
#include "rt/primitives/sketch_legacy_bsg.h"
#include "rt/view_legacy_bsg.h"
#include "./bsg_ged_draw_private.h"

void
ged_draw_view_info_from_bsg(struct rt_view_info *view_info,
			    const struct bsg_view *view)
{
    rt_view_info_from_bsg(view_info, view);
}

fastf_t
ged_draw_view_perspective_from_bsg(const struct bsg_view *view)
{
    return rt_view_perspective_from_bsg(view);
}

fastf_t
ged_draw_view_scale_from_bsg(const struct bsg_view *view)
{
    return rt_view_scale_from_bsg(view);
}

int
ged_draw_view_lod_policy_from_bsg(ged_draw_view_lod_policy *policy,
				  const struct bsg_view *view)
{
    return rt_view_lod_policy_from_bsg(policy, view);
}

int
ged_draw_view_lod_policy_apply_bsg(struct bsg_view *view,
				   const ged_draw_view_lod_policy *policy)
{
    return rt_view_lod_policy_apply_bsg(view, policy);
}

int
ged_draw_view_lod_policy_apply_bsg_bot_threshold(
	struct bsg_view *view,
	const ged_draw_view_lod_policy *policy,
	size_t bot_threshold)
{
    if (!policy)
	return 0;

    ged_draw_view_lod_policy override_policy = *policy;
    override_policy.bot_threshold = bot_threshold;
    return ged_draw_view_lod_policy_apply_bsg(view, &override_policy);
}

int
ged_draw_view_autoview_default_bsg(struct bsg_view *view, int all_view_objs)
{
    return rt_view_autoview_bsg(view, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, all_view_objs);
}

int
ged_draw_view_hud_sync(struct bsg_view *view)
{
    return bsg_hud_sync(view);
}

static struct bsg_selection *
ged_draw_view_selection_get_bsg(struct bsg_view *view)
{
    return bsg_view_selection(view);
}

int
ged_draw_view_selection_available(struct bsg_view *view)
{
    return ged_draw_view_selection_get_bsg(view) ? 1 : 0;
}

size_t
ged_draw_view_selection_count(struct bsg_view *view)
{
    return rt_view_selection_count_bsg(view);
}

int
ged_draw_view_selection_path_foreach(struct bsg_view *view,
				     ged_draw_view_selection_path_cb cb,
				     void *data)
{
    if (!cb)
	return 0;

    struct bsg_selection *selection = ged_draw_view_selection_get_bsg(view);
    if (!selection)
	return 0;

    for (size_t i = 0; i < bsg_selection_count(selection); i++) {
	const struct bsg_interaction_record *record =
	    bsg_selection_record(selection, i);
	const char *path = bsg_interaction_record_path(record);
	if (path && path[0] && !cb(view, path, data))
	    return 0;
    }

    return 1;
}

int
ged_draw_view_selection_clear(struct bsg_view *view)
{
    return rt_view_selection_clear_bsg(view);
}

int
ged_draw_view_selection_contains_record(
	struct bsg_view *view,
	const struct bsg_interaction_record *record)
{
    struct bsg_selection *selection = ged_draw_view_selection_get_bsg(view);
    return selection ? bsg_selection_contains_record(selection, record) : 0;
}

int
ged_draw_view_selection_add_record(
	struct bsg_view *view,
	const struct bsg_interaction_record *record)
{
    struct bsg_selection *selection = ged_draw_view_selection_get_bsg(view);
    if (!selection || !record)
	return 0;

    if (!bsg_selection_contains_record(selection, record))
	bsg_selection_add_record(selection, record);
    return 1;
}

int
ged_draw_view_selection_set_record(
	struct bsg_view *view,
	const struct bsg_interaction_record *record)
{
    struct bsg_selection *selection = ged_draw_view_selection_get_bsg(view);
    if (!selection)
	return 0;

    bsg_selection_clear(selection);
    if (record)
	bsg_selection_add_record(selection, record);
    return 1;
}

static int
_ged_draw_view_snap_kind_bsg(enum ged_draw_view_snap_kind kind)
{
    switch (kind) {
	case GED_DRAW_VIEW_SNAP_GRID:
	    return RT_VIEW_SNAP_KIND_GRID_BSG;
	case GED_DRAW_VIEW_SNAP_ENDPOINT:
	    return RT_VIEW_SNAP_KIND_ENDPOINT_BSG;
	default:
	    return 0;
    }
}

int
ged_draw_view_snap_first_candidate(struct bsg_view *view,
				   const point_t sample,
				   enum ged_draw_view_snap_kind kind,
				   point_t candidate)
{
    if (!candidate)
	return 0;

    VSET(candidate, 0.0, 0.0, 0.0);

    int bsg_kind = _ged_draw_view_snap_kind_bsg(kind);
    if (!view || !sample || !bsg_kind)
	return 0;

    point_t bsg_sample = VINIT_ZERO;
    VMOVE(bsg_sample, sample);

    struct bsg_snap_result sres = {0};
    int ret = 0;
    if (rt_view_snap_candidates_bsg(view, bsg_sample, 0.0, bsg_kind, &sres) > 0) {
	VMOVE(candidate, sres.sr_candidates[0].sc_point);
	ret = 1;
    }
    bsg_snap_result_free(&sres);
    return ret;
}

struct bu_ptbl *
ged_draw_view_set_views_bsg(struct bsg_view_set *view_set)
{
    return rt_view_set_views_bsg(view_set);
}

void *
ged_draw_view_set_recycle_pool_bsg(struct bsg_view_set *view_set)
{
    return bsg_set_fsos(view_set);
}

int
ged_draw_view_is_independent_bsg(const struct bsg_view *view)
{
    return rt_view_is_independent_bsg(view);
}

bsg_scene_ref
ged_draw_view_independent_scope_ref_bsg(struct bsg_view *view, int create)
{
    return rt_view_independent_scope_ref_bsg(view, create);
}

void
ged_draw_view_set_lod_bounds_update(struct bsg_view *view)
{
    rt_view_lod_bounds_callback_set_bsg(view);
}

int
ged_draw_view_has_lod_bounds_update(const struct bsg_view *view)
{
    return rt_view_lod_bounds_callback_is_bsg(view);
}

void
ged_draw_scene_ref_realization_set_bsg_view_policy(bsg_scene_ref ref,
						   const struct bsg_view *view)
{
    if (!view)
	return;

    ged_draw_view_lod_policy policy;
    ged_draw_view_lod_policy_from_bsg(&policy, view);
    ged_draw_scene_ref_realization_set_view_policy(ref,
	    policy.csg_enabled,
	    rt_view_scale_from_bsg(view),
	    policy.bot_threshold,
	    policy.curve_scale,
	    policy.point_scale);
}

int
ged_draw_mesh_lod_load_view_scene_ref(struct rt_mesh_lod *lod,
				      bsg_scene_ref visibility_ref,
				      struct bsg_view *view,
				      int reset)
{
    return rt_mesh_lod_load_view_scene_ref_bsg(lod, visibility_ref, view, reset);
}

void
ged_draw_mesh_lod_free_scene_ref(bsg_scene_ref ref)
{
    rt_mesh_lod_free_scene_ref_bsg(ref);
}

int
ged_draw_view_feature_exists(struct bsg_view *view, const char *name)
{
    if (!view || !name)
	return 0;

    bsg_feature_ref ref = bsg_feature_find(view, name);
    return !bsg_feature_ref_is_null(ref);
}

int
ged_draw_view_feature_remove(struct bsg_view *view, const char *name)
{
    if (!view || !name)
	return 0;

    return bsg_feature_remove(view, name) ? 1 : 0;
}

struct ged_draw_view_feature_prefix_remove {
    const char *prefix;
    size_t prefix_len;
    char **names;
    size_t count;
    size_t capacity;
};

static int
_ged_draw_view_feature_prefix_remove_cb(bsg_feature_ref UNUSED(ref),
					const struct bsg_feature_record *record,
					void *data)
{
    struct ged_draw_view_feature_prefix_remove *ctx =
	(struct ged_draw_view_feature_prefix_remove *)data;
    const char *name = record ? record->name : NULL;
    if (!ctx || !name || strncmp(ctx->prefix, name, ctx->prefix_len))
	return 1;

    if (ctx->count == ctx->capacity) {
	size_t next_capacity = ctx->capacity ? ctx->capacity * 2 : 8;
	ctx->names = (char **)bu_realloc(ctx->names,
		next_capacity * sizeof(char *),
		"GED draw view feature prefix removal names");
	ctx->capacity = next_capacity;
    }

    ctx->names[ctx->count++] = bu_strdup(name);
    return 1;
}

int
ged_draw_view_features_remove_prefix(struct bsg_view *view, const char *prefix)
{
    if (!view || !prefix || !*prefix)
	return 0;

    struct ged_draw_view_feature_prefix_remove ctx;
    ctx.prefix = prefix;
    ctx.prefix_len = strlen(prefix);
    ctx.names = NULL;
    ctx.count = 0;
    ctx.capacity = 0;

    bsg_feature_visit(view, BSG_FEATURE_SCOPE_ALL,
	    _ged_draw_view_feature_prefix_remove_cb, &ctx);

    int removed = 0;
    for (size_t i = 0; i < ctx.count; i++) {
	if (bsg_feature_remove(view, ctx.names[i]))
	    removed++;
	bu_free(ctx.names[i], "GED draw view feature prefix removal name");
    }
    if (ctx.names)
	bu_free(ctx.names, "GED draw view feature prefix removal names");

    return removed;
}

int
ged_draw_view_feature_visible(struct bsg_view *view, const char *name)
{
    if (!view || !name)
	return 0;

    bsg_feature_ref ref = bsg_feature_find(view, name);
    struct bsg_feature_record rec;
    return (!bsg_feature_ref_is_null(ref) &&
	    bsg_feature_record_get(ref, &rec) && rec.visible) ? 1 : 0;
}

bsg_scene_ref
ged_draw_view_overlay_create(struct bsg_view *view, const char *name)
{
    if (!view || !name)
	return bsg_scene_ref_null();

    bsg_feature_ref ref = bsg_feature_create_overlay(view, name, 0);
    if (bsg_feature_ref_is_null(ref))
	return bsg_scene_ref_null();

    return bsg_feature_ref_as_scene(ref);
}

int
ged_draw_view_feature_depth(struct bsg_view *view,
			    const char *name,
			    int mode,
			    fastf_t *depth)
{
    if (depth)
	*depth = 0.0;
    if (!view || !name || !depth)
	return 0;

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    *depth = bsg_feature_view_depth(ref, view, mode);
    return 1;
}

struct ged_draw_view_feature_depth_visit {
    struct bsg_view *view;
    int mode;
    ged_draw_view_feature_depth_cb cb;
    void *data;
    int count;
};

static int
_ged_draw_view_feature_depth_visit_cb(bsg_feature_ref ref,
				      const struct bsg_feature_record *UNUSED(record),
				      void *data)
{
    struct ged_draw_view_feature_depth_visit *ctx =
	(struct ged_draw_view_feature_depth_visit *)data;
    if (!ctx || !ctx->cb)
	return 0;

    ctx->count++;
    return ctx->cb(bsg_feature_view_depth(ref, ctx->view, ctx->mode),
	    ctx->data);
}

int
ged_draw_view_feature_depth_foreach(struct bsg_view *view,
				    int mode,
				    ged_draw_view_feature_depth_cb cb,
				    void *data)
{
    if (!view || !cb)
	return 0;

    struct ged_draw_view_feature_depth_visit ctx;
    ctx.view = view;
    ctx.mode = mode;
    ctx.cb = cb;
    ctx.data = data;
    ctx.count = 0;
    bsg_feature_visit(view, BSG_FEATURE_SCOPE_ALL,
	    _ged_draw_view_feature_depth_visit_cb, &ctx);
    return ctx.count;
}

static void
_ged_draw_view_feature_style_to_bsg(struct bsg_feature_style *dst,
				    const struct ged_draw_view_feature_style *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    dst->visible = src->visible;
    dst->color_valid = src->color_valid;
    dst->color[0] = src->color[0];
    dst->color[1] = src->color[1];
    dst->color[2] = src->color[2];
    dst->line_width = src->line_width;
    dst->line_style = src->line_style;
    dst->arrow = src->arrow;
    dst->arrow_tip_length = src->arrow_tip_length;
    dst->arrow_tip_width = src->arrow_tip_width;
}

static void
_ged_draw_view_feature_style_from_bsg(struct ged_draw_view_feature_style *dst,
				      const struct bsg_feature_style *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    dst->visible = src->visible;
    dst->color_valid = src->color_valid;
    dst->color[0] = src->color[0];
    dst->color[1] = src->color[1];
    dst->color[2] = src->color[2];
    dst->line_width = src->line_width;
    dst->line_style = src->line_style;
    dst->arrow = src->arrow;
    dst->arrow_tip_length = src->arrow_tip_length;
    dst->arrow_tip_width = src->arrow_tip_width;
}

int
ged_draw_view_feature_style_get(struct bsg_view *view,
				const char *name,
				struct ged_draw_view_feature_style *style)
{
    if (!view || !name || !style)
	return 0;

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    struct bsg_feature_style bsg_style = BSG_FEATURE_STYLE_INIT;
    if (!bsg_feature_style_get(ref, &bsg_style))
	return 0;

    _ged_draw_view_feature_style_from_bsg(style, &bsg_style);
    return 1;
}

int
ged_draw_view_feature_style_apply(struct bsg_view *view,
				  const char *name,
				  const struct ged_draw_view_feature_style *style,
				  int recursive)
{
    if (!view || !name || !style)
	return 0;

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    struct bsg_feature_style bsg_style = BSG_FEATURE_STYLE_INIT;
    _ged_draw_view_feature_style_to_bsg(&bsg_style, style);
    return (recursive ?
	    bsg_feature_style_apply_recursive(ref, &bsg_style) :
	    bsg_feature_style_apply(ref, &bsg_style)) ? 1 : 0;
}

int
ged_draw_view_feature_realize(struct bsg_view *view, const char *name, int recursive)
{
    if (!view || !name)
	return 0;

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    return bsg_feature_realize(ref, view, recursive) ? 1 : 0;
}

int
ged_draw_view_indexed_face_set_replace(struct bsg_view *view,
				       const char *name,
				       int local,
				       const point_t *points,
				       size_t point_count,
				       const vect_t *normals,
				       size_t normal_count,
				       const int *indices,
				       size_t index_count,
				       const struct ged_draw_view_feature_style *style)
{
    if (!view || !name || !points || !point_count || !indices || !index_count)
	return 0;

    struct bsg_feature_style bsg_style = BSG_FEATURE_STYLE_INIT;
    const struct bsg_feature_style *style_ptr = NULL;
    if (style) {
	_ged_draw_view_feature_style_to_bsg(&bsg_style, style);
	style_ptr = &bsg_style;
    }

    bsg_feature_ref ref = bsg_feature_replace_indexed_face_set(view, name,
	    local, points, point_count, normals, normal_count, indices,
	    index_count, style_ptr);
    return bsg_feature_ref_is_null(ref) ? 0 : 1;
}

int
ged_draw_view_lines_replace(struct bsg_view *view,
			    const char *name,
			    int local,
			    const point_t *points,
			    const int *cmds,
			    size_t point_count,
			    const struct ged_draw_view_feature_style *style)
{
    if (!view || !name)
	return 0;

    bsg_feature_remove(view, name);
    if (!points || !point_count)
	return 1;

    bsg_feature_ref ref = bsg_feature_create_lines(view, name, local);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    if (!bsg_feature_points_replace(ref, BSG_FEATURE_LINES, points, cmds,
		point_count)) {
	bsg_feature_remove(view, name);
	return 0;
    }

    if (style) {
	struct bsg_feature_style bsg_style = BSG_FEATURE_STYLE_INIT;
	_ged_draw_view_feature_style_to_bsg(&bsg_style, style);
	if (!bsg_feature_style_apply(ref, &bsg_style)) {
	    bsg_feature_remove(view, name);
	    return 0;
	}
    }

    return 1;
}

int
ged_draw_view_line_layer_builder_replace(struct bsg_view *view,
					 const char *name,
					 int local,
					 const struct bg_line_layer_builder *builder)
{
    if (!view || !name)
	return 0;

    if (!builder || !bg_line_layer_builder_point_count(builder)) {
	bsg_feature_remove(view, name);
	return 1;
    }

    bsg_feature_ref ref = bsg_feature_replace_line_layer_builder(view, name,
	    local, (const struct bsg_line_layer_builder *)builder, NULL);
    return bsg_feature_ref_is_null(ref) ? 0 : 1;
}

static void
_ged_draw_view_line_layer_to_bsg(struct bsg_feature_line_layer *dst,
				 const struct ged_draw_view_line_layer_data *src)
{
    struct bsg_feature_line_layer init = BSG_FEATURE_LINE_LAYER_INIT;
    *dst = init;
    if (!src)
	return;

    dst->name = src->name;
    dst->points = src->points;
    dst->commands = src->commands;
    dst->point_count = src->point_count;
    _ged_draw_view_feature_style_to_bsg(&dst->style, &src->style);
}

int
ged_draw_view_line_layers_replace(struct bsg_view *view,
				  const char *name,
				  int local,
				  const struct ged_draw_view_line_layer_data *layers,
				  size_t layer_count,
				  const struct ged_draw_view_feature_style *style)
{
    if (!view || !name)
	return 0;

    size_t live_layers = 0;
    for (size_t i = 0; layers && i < layer_count; i++) {
	if (layers[i].points && layers[i].point_count)
	    live_layers++;
    }
    if (!layers || !live_layers) {
	bsg_feature_remove(view, name);
	return 1;
    }

    struct bsg_feature_line_layer *bsg_layers =
	(struct bsg_feature_line_layer *)bu_calloc(layer_count,
		sizeof(struct bsg_feature_line_layer), "GED draw view line layers");
    for (size_t i = 0; i < layer_count; i++)
	_ged_draw_view_line_layer_to_bsg(&bsg_layers[i], &layers[i]);

    struct bsg_feature_style bsg_style = BSG_FEATURE_STYLE_INIT;
    const struct bsg_feature_style *style_ptr = NULL;
    if (style) {
	_ged_draw_view_feature_style_to_bsg(&bsg_style, style);
	style_ptr = &bsg_style;
    }

    bsg_feature_ref ref = bsg_feature_replace_line_layers(view, name, local,
	    (const struct bsg_feature_line_layer *)bsg_layers, layer_count,
	    style_ptr);
    bu_free(bsg_layers, "GED draw view line layers");
    return bsg_feature_ref_is_null(ref) ? 0 : 1;
}

int
ged_draw_view_lines_create_model_annotation(struct bsg_view *view,
					    const char *name,
					    int local,
					    const point_t point)
{
    if (!view || !name || !point)
	return 0;

    bsg_feature_ref ref = bsg_feature_create_lines(view, name, local);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    bsg_feature_overlay_register_owner(ref, NULL,
	    BSG_OVERLAY_ROLE_MODEL,
	    BSG_OVERLAY_CLASS_USER_ANNOTATION,
	    BSG_OVERLAY_LC_PERSISTENT,
	    BSG_OVERLAY_ORDER_MODEL,
	    NULL, 0);

    int cmd = BSG_GEOMETRY_LINE_MOVE;
    if (!bsg_feature_points_replace(ref, BSG_FEATURE_LINES,
		(const point_t *)point, &cmd, 1)) {
	bsg_feature_remove(view, name);
	return 0;
    }

    return 1;
}

int
ged_draw_view_lines_append_point(struct bsg_view *view,
				 const char *name,
				 const point_t point)
{
    if (!view || !name || !point)
	return 0;

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    point_t *points = NULL;
    int *cmds = NULL;
    size_t point_count = 0;
    if (!bsg_feature_points_copy(ref, &points, &cmds, &point_count))
	return 0;

    point_t *npoints = (point_t *)bu_calloc(point_count + 1, sizeof(point_t),
	    "GED draw view line append points");
    int *ncmds = (int *)bu_calloc(point_count + 1, sizeof(int),
	    "GED draw view line append cmds");
    for (size_t i = 0; i < point_count; i++) {
	VMOVE(npoints[i], points[i]);
	ncmds[i] = cmds ? cmds[i] :
	    ((i == 0) ? BSG_GEOMETRY_LINE_MOVE : BSG_GEOMETRY_LINE_DRAW);
    }
    VMOVE(npoints[point_count], point);
    ncmds[point_count] = BSG_GEOMETRY_LINE_DRAW;

    int ret = bsg_feature_points_replace(ref, BSG_FEATURE_LINES,
	    (const point_t *)npoints, ncmds, point_count + 1);

    if (points)
	bu_free(points, "bsg feature points copy");
    if (cmds)
	bu_free(cmds, "bsg feature cmds copy");
    bu_free(npoints, "GED draw view line append points");
    bu_free(ncmds, "GED draw view line append cmds");

    return ret ? 1 : 0;
}

int
ged_draw_view_label_create(struct bsg_view *view,
			   const char *name,
			   int local,
			   const char *text,
			   const point_t point,
			   const point_t target,
			   int has_target)
{
    if (!view || !name || !text || !point || (has_target && !target))
	return 0;

    bsg_feature_ref ref = bsg_feature_create_label(view, name, local);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    bsg_feature_set_color(ref, 255, 255, 0);

    struct bsg_feature_label_data label = BSG_FEATURE_LABEL_DATA_INIT;
    label.text = text;
    VMOVE(label.point, point);
    if (has_target) {
	VMOVE(label.target, target);
	label.line_flag = 1;
    }

    return bsg_feature_labels_replace(ref, &label, 1) ? 1 : 0;
}

static void
_ged_draw_view_label_data_to_bsg(struct bsg_feature_label_data *dst,
				 const struct ged_draw_view_label_data *src)
{
    struct bsg_feature_label_data init = BSG_FEATURE_LABEL_DATA_INIT;
    *dst = init;
    if (!src)
	return;

    dst->text = src->text;
    VMOVE(dst->point, src->point);
    dst->color_valid = src->color_valid;
    dst->color[0] = src->color[0];
    dst->color[1] = src->color[1];
    dst->color[2] = src->color[2];
    dst->line_flag = src->line_flag;
    VMOVE(dst->target, src->target);
    dst->anchor = src->anchor;
    dst->arrow = src->arrow;
}

int
ged_draw_view_labels_replace(struct bsg_view *view,
			     const char *name,
			     int local,
			     const struct ged_draw_view_label_data *labels,
			     size_t label_count)
{
    if (!view || !name)
	return 0;

    bsg_feature_remove(view, name);
    if (!labels || !label_count)
	return 1;

    bsg_feature_ref ref = bsg_feature_create_label(view, name, local);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    struct bsg_feature_label_data *bsg_labels =
	(struct bsg_feature_label_data *)bu_calloc(label_count,
		sizeof(struct bsg_feature_label_data), "GED draw view label data");
    for (size_t i = 0; i < label_count; i++)
	_ged_draw_view_label_data_to_bsg(&bsg_labels[i], &labels[i]);

    int ret = bsg_feature_labels_replace(ref, bsg_labels, label_count);
    bu_free(bsg_labels, "GED draw view label data");
    if (!ret) {
	bsg_feature_remove(view, name);
	return 0;
    }

    return 1;
}

int
ged_draw_view_line_style_get(struct bsg_view *view,
			     const char *name,
			     struct ged_draw_view_line_style *style)
{
    if (!view || !name || !style)
	return 0;

    bsg_feature_ref ref = bsg_feature_find(view, name);
    struct bsg_feature_record rec;
    if (bsg_feature_ref_is_null(ref) || !bsg_feature_record_get(ref, &rec))
	return 0;

    style->color[0] = (int)rec.color[0];
    style->color[1] = (int)rec.color[1];
    style->color[2] = (int)rec.color[2];
    style->line_width = rec.line_width;
    return 1;
}

int
ged_draw_view_line_color_set(struct bsg_view *view,
			     const char *name,
			     int r,
			     int g,
			     int b)
{
    if (!view || !name)
	return 0;

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    bsg_feature_set_color(ref, r, g, b);
    return 1;
}

int
ged_draw_view_line_width_set(struct bsg_view *view,
			     const char *name,
			     int line_width)
{
    if (!view || !name)
	return 0;

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    bsg_feature_set_line_width(ref, line_width);
    return 1;
}

int
ged_draw_view_lines_points_copy(struct bsg_view *view,
				const char *name,
				point_t **points,
				size_t *point_count)
{
    if (points)
	*points = NULL;
    if (point_count)
	*point_count = 0;
    if (!view || !name || !points || !point_count)
	return 0;

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    return bsg_feature_points_copy(ref, points, NULL, point_count) ? 1 : 0;
}

int
ged_draw_view_tcl_lines_replace(struct bsg_view *view,
				const char *name,
				const point_t *points,
				size_t point_count,
				const struct ged_draw_view_line_style *style)
{
    if (!view || !name)
	return 0;

    if (point_count % 2)
	return 0;

    bsg_feature_remove(view, name);

    if (point_count < 2 || !points)
	return 1;

    bsg_feature_ref ref = bsg_feature_create_lines(view, name, 1);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    bsg_feature_overlay_register_owner(ref, NULL,
	    BSG_OVERLAY_ROLE_MODEL,
	    BSG_OVERLAY_CLASS_TCL_OVERLAY,
	    BSG_OVERLAY_LC_PER_COMMAND,
	    BSG_OVERLAY_ORDER_POST_TRANSPARENT,
	    NULL, 0);

    int *cmds = (int *)bu_calloc(point_count, sizeof(int),
	    "GED draw Tcl data line cmds");
    for (size_t i = 0; i + 1 < point_count; i += 2) {
	cmds[i] = BSG_GEOMETRY_LINE_MOVE;
	cmds[i + 1] = BSG_GEOMETRY_LINE_DRAW;
    }

    int ret = bsg_feature_points_replace(ref, BSG_FEATURE_LINES,
	    points, cmds, point_count);
    bu_free(cmds, "GED draw Tcl data line cmds");
    if (!ret)
	return 0;

    if (style)
	bsg_feature_set_color(ref, style->color[0], style->color[1],
		style->color[2]);
    if (style)
	bsg_feature_set_line_width(ref, style->line_width);
    bsg_feature_set_visible(ref, 1);

    return 1;
}

static void
_ged_draw_view_axes_to_bsg(struct bsg_axes *dst,
			   const struct ged_draw_view_axes_state *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    VMOVE(dst->axes_pos, src->position);
    dst->axes_size = src->size;
    dst->line_width = src->line_width;
    dst->axes_color[0] = src->color[0];
    dst->axes_color[1] = src->color[1];
    dst->axes_color[2] = src->color[2];
}

static void
_ged_draw_view_axes_from_bsg(struct ged_draw_view_axes_state *dst,
			     const struct bsg_axes *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    VMOVE(dst->position, src->axes_pos);
    dst->size = src->axes_size;
    dst->line_width = src->line_width;
    dst->color[0] = src->axes_color[0];
    dst->color[1] = src->axes_color[1];
    dst->color[2] = src->axes_color[2];
}

int
ged_draw_view_axes_create(struct bsg_view *view,
			  const char *name,
			  int local,
			  const struct ged_draw_view_axes_state *state)
{
    if (!view || !name || !state)
	return 0;

    bsg_feature_ref ref = bsg_feature_create_axes(view, name, local);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    struct bsg_axes axes;
    _ged_draw_view_axes_to_bsg(&axes, state);
    if (!bsg_feature_axes_state_replace(ref, &axes)) {
	bsg_feature_remove(view, name);
	return 0;
    }

    return 1;
}

int
ged_draw_view_axes_state_get(struct bsg_view *view,
			     const char *name,
			     struct ged_draw_view_axes_state *state)
{
    if (!view || !name || !state)
	return 0;

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    struct bsg_axes axes;
    if (!bsg_feature_axes_state_get(ref, &axes))
	return 0;

    _ged_draw_view_axes_from_bsg(state, &axes);
    return 1;
}

int
ged_draw_view_axes_state_replace(struct bsg_view *view,
				 const char *name,
				 const struct ged_draw_view_axes_state *state)
{
    if (!view || !name || !state)
	return 0;

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    struct bsg_axes axes;
    _ged_draw_view_axes_to_bsg(&axes, state);
    return bsg_feature_axes_state_replace(ref, &axes) ? 1 : 0;
}

static bsg_polygon_ref
_ged_draw_view_polygon_ref_to_bsg(ged_draw_view_polygon_ref ref)
{
    bsg_polygon_ref bsg_ref = { ref.token, ref.revision };
    return bsg_ref;
}

static ged_draw_view_polygon_ref
_ged_draw_view_polygon_ref_from_bsg(bsg_polygon_ref ref)
{
    ged_draw_view_polygon_ref ged_ref = { ref.token, ref.revision };
    return ged_ref;
}

int
ged_draw_view_polygon_ref_is_null(ged_draw_view_polygon_ref ref)
{
    return bsg_polygon_ref_is_null(_ged_draw_view_polygon_ref_to_bsg(ref));
}

ged_draw_view_polygon_ref
ged_draw_view_polygon_find(struct bsg_view *view, const char *name)
{
    return _ged_draw_view_polygon_ref_from_bsg(
	    bsg_view_polygon_find_ref(view, name));
}

ged_draw_view_polygon_ref
ged_draw_view_polygon_find_scoped(struct bsg_view *view,
				  const char *name,
				  int local_only)
{
    return _ged_draw_view_polygon_ref_from_bsg(
	    bsg_view_polygon_find_scoped_ref(view, name, local_only));
}

ged_draw_view_polygon_ref
ged_draw_view_polygon_create(struct bsg_view *view,
			     const char *name,
			     int local,
			     int type,
			     const point_t screen_point)
{
    point_t point;

    if (!screen_point)
	return GED_DRAW_VIEW_POLYGON_REF_NULL;

    VMOVE(point, screen_point);
    return _ged_draw_view_polygon_ref_from_bsg(
	    bsg_create_polygon_ref(view, name, local, type, &point));
}

ged_draw_view_polygon_ref
ged_draw_view_polygon_import_sketch(const char *name,
				    struct db_i *dbip,
				    struct directory *dp,
				    struct bsg_view *view,
				    int local)
{
    return _ged_draw_view_polygon_ref_from_bsg(
	    db_sketch_to_view_polygon_scoped_ref(name, dbip, dp, view, local));
}

int
ged_draw_view_polygon_export_sketch(struct db_i *dbip,
				    const char *name,
				    ged_draw_view_polygon_ref ref)
{
    return db_view_polygon_ref_to_sketch(dbip, name,
	    _ged_draw_view_polygon_ref_to_bsg(ref)) ? 1 : 0;
}

int
ged_draw_view_polygon_record_get(ged_draw_view_polygon_ref ref,
				 struct ged_draw_view_polygon_record *record)
{
    struct bsg_polygon_record bsg_record;

    if (!record)
	return 0;

    memset(record, 0, sizeof(*record));
    if (!bsg_polygon_record_get(_ged_draw_view_polygon_ref_to_bsg(ref),
	    &bsg_record))
	return 0;

    record->ref = _ged_draw_view_polygon_ref_from_bsg(bsg_record.ref);
    record->name = bsg_record.name;
    record->type = bsg_record.type;
    record->fill_flag = bsg_record.fill_flag;
    V2MOVE(record->fill_dir, bsg_record.fill_dir);
    record->fill_delta = bsg_record.fill_delta;
    record->fill_color = bsg_record.fill_color;
    record->edge_color[0] = bsg_record.edge_color[0];
    record->edge_color[1] = bsg_record.edge_color[1];
    record->edge_color[2] = bsg_record.edge_color[2];
    record->curr_contour_i = bsg_record.curr_contour_i;
    record->curr_point_i = bsg_record.curr_point_i;
    record->first_contour_open = bsg_record.first_contour_open;
    record->contour_count = bsg_record.contour_count;
    record->point_count = bsg_record.point_count;
    VMOVE(record->origin_point, bsg_record.origin_point);
    HMOVE(record->vp, bsg_record.vp);
    record->vZ = bsg_record.vZ;
    record->user_data = bsg_record.user_data;

    return 1;
}

int
ged_draw_view_polygon_has_data(ged_draw_view_polygon_ref ref)
{
    return bsg_polygon_data(_ged_draw_view_polygon_ref_to_bsg(ref)) ? 1 : 0;
}

int
ged_draw_view_polygon_update(ged_draw_view_polygon_ref ref,
			     struct bsg_view *view,
			     int op)
{
    return bsg_polygon_update(_ged_draw_view_polygon_ref_to_bsg(ref), view, op) ?
	1 : 0;
}

int
ged_draw_view_polygon_update_screen_pt(ged_draw_view_polygon_ref ref,
				       struct bsg_view *view,
				       int x,
				       int y,
				       int op)
{
    return bsg_polygon_update_screen_pt(_ged_draw_view_polygon_ref_to_bsg(ref),
	    view, x, y, op) ? 1 : 0;
}

int
ged_draw_view_polygon_set_current(ged_draw_view_polygon_ref ref,
				  long contour_i,
				  long point_i)
{
    return bsg_polygon_set_current(_ged_draw_view_polygon_ref_to_bsg(ref),
	    contour_i, point_i) ? 1 : 0;
}

int
ged_draw_view_polygon_set_contour_open(ged_draw_view_polygon_ref ref,
				       long contour_i,
				       int open)
{
    return bsg_polygon_set_contour_open(_ged_draw_view_polygon_ref_to_bsg(ref),
	    contour_i, open) ? 1 : 0;
}

int
ged_draw_view_polygon_set_all_contours_open(ged_draw_view_polygon_ref ref,
					    int open)
{
    return bsg_polygon_set_all_contours_open(
	    _ged_draw_view_polygon_ref_to_bsg(ref), open) ? 1 : 0;
}

int
ged_draw_view_polygon_area(ged_draw_view_polygon_ref ref,
			   struct bsg_view *view,
			   fastf_t *area)
{
    const struct bsg_polygon *poly =
	bsg_polygon_data(_ged_draw_view_polygon_ref_to_bsg(ref));

    if (!area)
	return 0;
    *area = 0.0;

    if (!poly || !view)
	return 0;

    *area = bg_find_polygon_area((struct bg_polygon *)&poly->polygon,
	    CLIPPER_MAX, (plane_t *)&poly->vp, rt_view_scale_from_bsg(view));
    return 1;
}

int
ged_draw_view_polygon_overlap(ged_draw_view_polygon_ref ref,
			      struct bsg_view *view,
			      const char *other_name,
			      const struct bn_tol *tol,
			      int *overlap)
{
    if (!overlap)
	return 0;
    *overlap = 0;

    if (!view || !other_name || !tol)
	return 0;

    const struct bsg_polygon *poly_a =
	bsg_polygon_data(_ged_draw_view_polygon_ref_to_bsg(ref));
    bsg_polygon_ref other_ref = bsg_view_polygon_find_ref(view, other_name);
    const struct bsg_polygon *poly_b = bsg_polygon_data(other_ref);
    if (!poly_a || !poly_b)
	return 0;

    *overlap = bg_polygons_overlap((struct bg_polygon *)&poly_a->polygon,
	    (struct bg_polygon *)&poly_b->polygon, (plane_t *)&poly_a->vp,
	    tol, rt_view_scale_from_bsg(view));
    return 1;
}

int
ged_draw_view_polygon_set_fill(ged_draw_view_polygon_ref ref,
			       int fill_flag,
			       fastf_t fill_slope_x,
			       fastf_t fill_slope_y,
			       fastf_t fill_density)
{
    return bsg_polygon_set_fill(_ged_draw_view_polygon_ref_to_bsg(ref),
	    fill_flag, fill_slope_x, fill_slope_y, fill_density) ? 1 : 0;
}

int
ged_draw_view_polygon_fill_color_get(ged_draw_view_polygon_ref ref,
				     struct bu_color *fill_color)
{
    return bsg_polygon_fill_color_get(_ged_draw_view_polygon_ref_to_bsg(ref),
	    fill_color) ? 1 : 0;
}

int
ged_draw_view_polygon_fill_color_set(ged_draw_view_polygon_ref ref,
				     const struct bu_color *fill_color)
{
    return bsg_polygon_fill_color_set(_ged_draw_view_polygon_ref_to_bsg(ref),
	    fill_color) ? 1 : 0;
}

int
ged_draw_view_polygon_csg(ged_draw_view_polygon_ref target,
			  struct bsg_view *view,
			  const char *other_name,
			  bg_clip_t op)
{
    if (!view || !other_name)
	return 0;

    bsg_polygon_ref other_ref = bsg_view_polygon_find_ref(view, other_name);
    return bsg_polygon_csg_ref(_ged_draw_view_polygon_ref_to_bsg(target),
	    other_ref, op) ? 1 : 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
