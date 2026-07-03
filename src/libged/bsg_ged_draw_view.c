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
#include "bu/color.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "bsg/feature.h"
#include "bsg/geometry.h"
#include "bsg/hud.h"
#include "bsg/interaction.h"
#include "bsg/overlay.h"
#include "bsg/polygon.h"
#include "bsg/selection.h"
#include "bsg/view_set.h"
#include "bsg/view_state.h"
#include "rt/view.h"
#include "rt/primitives/sketch_legacy_bsg.h"
#include "./bsg_ged_draw_private.h"
#include "./bsg_ged_draw_view_private.h"

int
ged_draw_view_context_hud_sync(void *view_ctx)
{
    return bsg_hud_sync((struct bsg_view *)view_ctx);
}

static struct bsg_selection *
ged_draw_view_selection_get_bsg(struct bsg_view *view)
{
    return bsg_view_selection(view);
}

static int
_ged_draw_view_selection_available(struct bsg_view *view)
{
    return ged_draw_view_selection_get_bsg(view) ? 1 : 0;
}

int
ged_draw_view_context_selection_available(void *view_ctx)
{
    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_selection_available(view_ctx);

    return _ged_draw_view_selection_available((struct bsg_view *)view_ctx);
}

static size_t
_ged_draw_view_selection_count(struct bsg_view *view)
{
    struct bsg_selection *selection = ged_draw_view_selection_get_bsg(view);
    return selection ? bsg_selection_count(selection) : 0;
}

size_t
ged_draw_view_context_selection_count(void *view_ctx)
{
    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_selection_count(view_ctx);

    return _ged_draw_view_selection_count((struct bsg_view *)view_ctx);
}

int
ged_draw_view_context_selection_path_foreach(
	void *view_ctx,
	ged_draw_view_context_selection_path_cb cb,
	void *data)
{
    if (!cb)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_selection_path_foreach(view_ctx,
		cb, data);

    struct bsg_selection *selection =
	ged_draw_view_selection_get_bsg((struct bsg_view *)view_ctx);
    if (!selection)
	return 0;

    for (size_t i = 0; i < bsg_selection_count(selection); i++) {
	const struct bsg_interaction_record *record =
	    bsg_selection_record(selection, i);
	const char *path = bsg_interaction_record_path(record);
	if (path && path[0] && !cb(view_ctx, path, data))
	    return 0;
    }

    return 1;
}

static int
_ged_draw_view_selection_clear(struct bsg_view *view)
{
    struct bsg_selection *selection = ged_draw_view_selection_get_bsg(view);
    if (!selection)
	return 0;

    bsg_selection_clear(selection);
    return 1;
}

int
ged_draw_view_context_selection_clear(void *view_ctx)
{
    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_selection_clear(view_ctx);

    return _ged_draw_view_selection_clear((struct bsg_view *)view_ctx);
}

static int
_ged_draw_view_selection_kind_bsg(enum ged_draw_view_selection_kind kind)
{
    switch (kind) {
	case GED_DRAW_VIEW_SELECTION_SELECTED_PATH:
	    return BSG_INTERACTION_SELECTED_PATH;
	case GED_DRAW_VIEW_SELECTION_HIGHLIGHTED_REF:
	    return BSG_INTERACTION_HIGHLIGHTED_REF;
	default:
	    return 0;
    }
}

static struct bsg_interaction_record *
_ged_draw_view_selection_record_create(struct bsg_view *view,
				       enum ged_draw_view_selection_kind kind,
				       const char *path)
{
    if (!path || !path[0])
	return NULL;

    int bsg_kind = _ged_draw_view_selection_kind_bsg(kind);
    if (!bsg_kind)
	return NULL;

    bsg_feature_ref feature = BSG_FEATURE_REF_NULL_INIT;
    return bsg_interaction_record_create_ref(view,
	    (bsg_interaction_kind)bsg_kind,
	    feature,
	    ged_draw_dbpath_skip_lead_slash(path),
	    NULL);
}

static int
_ged_draw_view_selection_contains_path(
	struct bsg_view *view,
	enum ged_draw_view_selection_kind kind,
	const char *path)
{
    struct bsg_selection *selection = ged_draw_view_selection_get_bsg(view);
    if (!selection)
	return 0;

    struct bsg_interaction_record *record =
	_ged_draw_view_selection_record_create(view, kind, path);
    if (!record)
	return 0;

    int ret = bsg_selection_contains_record(selection, record) ? 1 : 0;
    bsg_interaction_record_free(record);
    return ret;
}

int
ged_draw_view_context_selection_contains_path(
	void *view_ctx,
	enum ged_draw_view_selection_kind kind,
	const char *path)
{
    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_selection_contains_path(view_ctx,
		(int)kind, path);

    return _ged_draw_view_selection_contains_path((struct bsg_view *)view_ctx,
	    kind, path);
}

static int
_ged_draw_view_selection_add_path(
	struct bsg_view *view,
	enum ged_draw_view_selection_kind kind,
	const char *path)
{
    struct bsg_selection *selection = ged_draw_view_selection_get_bsg(view);
    if (!selection)
	return 0;

    struct bsg_interaction_record *record =
	_ged_draw_view_selection_record_create(view, kind, path);
    if (!record)
	return 0;

    if (!bsg_selection_contains_record(selection, record))
	bsg_selection_add_record(selection, record);
    bsg_interaction_record_free(record);
    return 1;
}

int
ged_draw_view_context_selection_add_path(
	void *view_ctx,
	enum ged_draw_view_selection_kind kind,
	const char *path)
{
    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_selection_add_path(view_ctx,
		(int)kind, path);

    return _ged_draw_view_selection_add_path((struct bsg_view *)view_ctx,
	    kind, path);
}

static int
_ged_draw_view_selection_set_path(
	struct bsg_view *view,
	enum ged_draw_view_selection_kind kind,
	const char *path)
{
    struct bsg_selection *selection = ged_draw_view_selection_get_bsg(view);
    if (!selection)
	return 0;

    bsg_selection_clear(selection);
    if (path && path[0]) {
	struct bsg_interaction_record *record =
	    _ged_draw_view_selection_record_create(view, kind, path);
	if (!record)
	    return 0;
	bsg_selection_add_record(selection, record);
	bsg_interaction_record_free(record);
    }
    return 1;
}

int
ged_draw_view_context_selection_set_path(
	void *view_ctx,
	enum ged_draw_view_selection_kind kind,
	const char *path)
{
    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_selection_set_path(view_ctx,
		(int)kind, path);

    return _ged_draw_view_selection_set_path((struct bsg_view *)view_ctx,
	    kind, path);
}

int
ged_draw_view_context_feature_exists(void *view_ctx, const char *name)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_feature_exists(view_ctx, name);

    bsg_feature_ref ref = bsg_feature_find(view, name);
    return !bsg_feature_ref_is_null(ref);
}

int
ged_draw_view_context_feature_remove(void *view_ctx, const char *name)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_feature_remove(view_ctx, name);

    return bsg_feature_remove(view, name) ? 1 : 0;
}

int
ged_draw_view_context_feature_summary(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_feature_summary *summary)
{
    if (!summary)
	return 0;

    memset(summary, 0, sizeof(*summary));
    if (!view_ctx || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_feature_summary(view_ctx, name,
		summary);

    bsg_feature_ref ref = bsg_feature_find((struct bsg_view *)view_ctx, name);
    if (bsg_feature_ref_is_null(ref))
	return 1;

    struct bsg_feature_record record;
    memset(&record, 0, sizeof(record));
    if (!bsg_feature_record_get(ref, &record))
	return 0;

    summary->exists = 1;
    summary->is_overlay = (record.family == BSG_FEATURE_OVERLAY);
    summary->is_label = (record.family == BSG_FEATURE_LABEL);
    summary->is_transient_preview =
	(record.family == BSG_FEATURE_TRANSIENT_PREVIEW);
    summary->visible = record.visible;
    summary->color[0] = record.color[0];
    summary->color[1] = record.color[1];
    summary->color[2] = record.color[2];
    summary->child_count = record.child_count;
    summary->geometry_command_count = record.geometry_command_count;

    return 1;
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
ged_draw_view_context_features_remove_prefix(void *view_ctx, const char *prefix)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !prefix || !*prefix)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return (int)ged_draw_obol_view_context_features_remove_prefix(
		view_ctx, prefix);

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
ged_draw_view_context_feature_visible(void *view_ctx, const char *name)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_feature_visible(view_ctx, name);

    bsg_feature_ref ref = bsg_feature_find(view, name);
    struct bsg_feature_record rec;
    return (!bsg_feature_ref_is_null(ref) &&
	    bsg_feature_record_get(ref, &rec) && rec.visible) ? 1 : 0;
}

int
ged_draw_view_context_feature_visible_set(void *view_ctx, const char *name, int visible)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_feature_visible_set(view_ctx, name,
		visible);

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    bsg_feature_set_visible(ref, visible ? 1 : 0);
    return 1;
}

int
ged_draw_view_context_feature_depth(void *view_ctx,
				    const char *name,
				    int mode,
				    fastf_t *depth)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (depth)
	*depth = 0.0;
    if (!view || !name || !depth)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_feature_depth(view_ctx, name, mode,
		depth);

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
ged_draw_view_context_feature_depth_foreach(
	void *view_ctx,
	int mode,
	ged_draw_view_feature_depth_cb cb,
	void *data)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !cb)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_feature_depth_foreach(view_ctx,
		mode, cb, data);

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
ged_draw_view_context_feature_style_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_feature_style *style)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name || !style)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_feature_style_get(view_ctx, name,
		style);

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
ged_draw_view_context_feature_style_apply(
	void *view_ctx,
	const char *name,
	const struct ged_draw_view_feature_style *style,
	int recursive)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name || !style)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_feature_style_apply(view_ctx, name,
		style, recursive);

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
ged_draw_view_context_feature_realize(void *view_ctx,
				      const char *name,
				      int recursive)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_feature_realize(view_ctx, name,
		recursive);

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    return bsg_feature_realize(ref, view, recursive) ? 1 : 0;
}

int
ged_draw_view_context_indexed_face_set_replace(
	void *view_ctx,
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
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name || !points || !point_count || !indices || !index_count)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_indexed_face_set_replace(view_ctx,
		name, local, points, point_count, normals, normal_count,
		indices, index_count, style);

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
ged_draw_view_context_lines_replace(void *view_ctx,
				    const char *name,
				    int local,
				    const point_t *points,
				    const int *cmds,
				    size_t point_count,
				    const struct ged_draw_view_feature_style *style)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_lines_replace(view_ctx, name,
		local, points, cmds, point_count, style);

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
ged_draw_view_context_tcl_polygons_replace(
	void *view_ctx,
	const char *name,
	const point_t *points,
	const int *cmds,
	size_t point_count,
	const struct ged_draw_view_feature_style *style)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_tcl_polygons_replace(view_ctx,
		name, points, cmds, point_count, style);

    bsg_feature_remove(view, name);
    if (!points || !cmds || !point_count)
	return 1;

    int *bsg_cmds = (int *)bu_calloc(point_count, sizeof(int),
	    "GED draw Tcl polygon feature cmds");
    for (size_t i = 0; i < point_count; i++) {
	bsg_cmds[i] = (cmds[i] == BG_GEOMETRY_LINE_MOVE) ?
	    BSG_GEOMETRY_LINE_MOVE : BSG_GEOMETRY_LINE_DRAW;
    }

    struct bsg_feature_style bsg_style = BSG_FEATURE_STYLE_INIT;
    const struct bsg_feature_style *bsg_stylep = NULL;
    if (style) {
	_ged_draw_view_feature_style_to_bsg(&bsg_style, style);
	bsg_stylep = &bsg_style;
    }

    bsg_feature_ref ref = bsg_feature_replace_lines(view, name, 1 /* local */,
	    points, point_count, bsg_stylep);
    int ret = 0;
    if (!bsg_feature_ref_is_null(ref)) {
	ret = bsg_feature_points_replace(ref, BSG_FEATURE_POLYGON, points,
		bsg_cmds, point_count) ? 1 : 0;
	if (ret) {
	    bsg_feature_overlay_register_owner(ref, NULL,
		    BSG_OVERLAY_ROLE_MODEL,
		    BSG_OVERLAY_CLASS_POLYGON_EDIT,
		    BSG_OVERLAY_LC_PER_TOOL,
		    BSG_OVERLAY_ORDER_POST_TRANSPARENT,
		    NULL, 0);
	}
    }
    bu_free(bsg_cmds, "GED draw Tcl polygon feature cmds");

    if (!ret)
	bsg_feature_remove(view, name);

    return ret;
}

int
ged_draw_view_context_line_layer_builder_replace(
	void *view_ctx,
	const char *name,
	int local,
	const struct bg_line_layer_builder *builder)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_line_layer_builder_replace(
		view_ctx, name, local, builder);

    if (!builder || !bg_line_layer_builder_point_count(builder)) {
	bsg_feature_remove(view, name);
	return 1;
    }

    bsg_feature_ref ref = bsg_feature_replace_line_layer_builder(view, name,
	    local, (const struct bsg_line_layer_builder *)builder, NULL);
    return bsg_feature_ref_is_null(ref) ? 0 : 1;
}

int
ged_draw_view_context_diagnostic_line_layer_builder_replace(
	void *view_ctx,
	const char *name,
	const struct bg_line_layer_builder *builder)
{
    if (!view_ctx || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_diagnostic_line_layer_builder_replace(
		view_ctx, name, builder);

    return ged_draw_view_context_line_layer_builder_replace(view_ctx, name,
	    0, builder);
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
ged_draw_view_context_line_layers_replace(void *view_ctx,
					  const char *name,
					  int local,
					  const struct ged_draw_view_line_layer_data *layers,
					  size_t layer_count,
					  const struct ged_draw_view_feature_style *style)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_line_layers_replace(view_ctx, name,
		local, layers, layer_count, style);

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
ged_draw_view_context_lines_create_model_annotation(
	void *view_ctx,
	const char *name,
	int local,
	const point_t point)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name || !point)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_lines_create_model_annotation(
		view_ctx, name, local, point);

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
ged_draw_view_context_lines_append_point(void *view_ctx,
					 const char *name,
					 const point_t point)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name || !point)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_lines_append_point(view_ctx, name,
		point);

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
ged_draw_view_context_label_create(void *view_ctx,
				   const char *name,
				   int local,
				   const char *text,
				   const point_t point,
				   const point_t target,
				   int has_target)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name || !text || !point || (has_target && !target))
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_label_create(view_ctx, name, local,
		text, point, target, has_target);

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
    dst->font_size = src->font_size;
}

int
ged_draw_view_context_labels_replace(
	void *view_ctx,
	const char *name,
	int local,
	const struct ged_draw_view_label_data *labels,
	size_t label_count)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_labels_replace(view_ctx, name,
		local, labels, label_count);

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
ged_draw_view_context_tcl_labels_replace(
	void *view_ctx,
	const char *name,
	int draw,
	const struct ged_draw_view_label_data *labels,
	size_t label_count)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_tcl_labels_replace(view_ctx, name,
		draw, labels, label_count);

    if (!draw || !labels || !label_count) {
	bsg_feature_remove(view, name);
	return 1;
    }

    return ged_draw_view_context_labels_replace(view_ctx, name, 1, labels,
	    label_count);
}

size_t
ged_draw_view_context_label_count(void *view_ctx, const char *name)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_label_count(view_ctx, name);

    bsg_feature_ref ref = bsg_feature_find(view, name);
    return bsg_feature_label_count(ref);
}

int
ged_draw_view_context_label_copy(void *view_ctx,
				 const char *name,
				 size_t index,
				 struct bu_vls *text,
				 point_t point,
				 unsigned char rgb[3])
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_label_copy(view_ctx, name, index,
		text, point, rgb);

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    return bsg_feature_label_copy(ref, index, text, point, rgb) ? 1 : 0;
}

int
ged_draw_view_context_label_point_set(void *view_ctx,
				      const char *name,
				      size_t index,
				      const point_t point)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name || !point)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_label_point_set(view_ctx, name,
		index, point);

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    return bsg_feature_label_point_set(ref, index, point) ? 1 : 0;
}

int
ged_draw_view_context_line_style_get(void *view_ctx,
				     const char *name,
				     struct ged_draw_view_line_style *style)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name || !style)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_line_style_get(view_ctx, name,
		style);

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
ged_draw_view_context_line_color_set(void *view_ctx,
				     const char *name,
				     int r,
				     int g,
				     int b)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_line_color_set(view_ctx, name, r,
		g, b);

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    bsg_feature_set_color(ref, r, g, b);
    return 1;
}

int
ged_draw_view_context_line_width_set(void *view_ctx,
				     const char *name,
				     int line_width)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_line_width_set(view_ctx, name,
		line_width);

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    bsg_feature_set_line_width(ref, line_width);
    return 1;
}

int
ged_draw_view_context_feature_points_copy(void *view_ctx,
					  const char *name,
					  point_t **points,
					  size_t *point_count)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (points)
	*points = NULL;
    if (point_count)
	*point_count = 0;
    if (!view || !name || !points || !point_count)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_feature_points_copy(view_ctx, name,
		points, point_count);

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    return bsg_feature_points_copy(ref, points, NULL, point_count) ? 1 : 0;
}

static int
ged_draw_view_line_command_from_bsg(int command)
{
    switch (command) {
	case BSG_GEOMETRY_LINE_MOVE:
	    return GED_DRAW_VIEW_LINE_MOVE;
	case BSG_GEOMETRY_LINE_DRAW:
	    return GED_DRAW_VIEW_LINE_DRAW;
	default:
	    return command;
    }
}

int
ged_draw_view_context_feature_line_command_at(
	void *view_ctx,
	const char *name,
	size_t index,
	int *out)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    point_t *points = NULL;
    int *cmds = NULL;
    size_t point_count = 0;
    int ret = 0;

    if (out)
	*out = 0;
    if (!view || !name || !out)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_feature_line_command_at(view_ctx,
		name, index, out);

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    if (!bsg_feature_points_copy(ref, &points, &cmds, &point_count))
	goto cleanup;
    if (index >= point_count || !cmds)
	goto cleanup;

    *out = ged_draw_view_line_command_from_bsg(cmds[index]);
    ret = 1;

cleanup:
    if (points)
	bu_free(points, "GED draw view feature command points");
    if (cmds)
	bu_free(cmds, "GED draw view feature commands");
    return ret;
}

int
ged_draw_view_context_lines_points_copy(void *view_ctx,
					const char *name,
					point_t **points,
					size_t *point_count)
{
    return ged_draw_view_context_feature_points_copy(view_ctx, name, points,
	    point_count);
}

int
ged_draw_view_context_tcl_lines_replace(
	void *view_ctx,
	const char *name,
	const point_t *points,
	size_t point_count,
	const struct ged_draw_view_line_style *style)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_tcl_lines_replace(view_ctx, name,
		points, point_count, style);

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

int
ged_draw_view_context_arrow_tip_get(void *view_ctx,
				    const char *name,
				    fastf_t *tip_length,
				    fastf_t *tip_width)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (tip_length)
	*tip_length = 0.0;
    if (tip_width)
	*tip_width = 0.0;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_arrow_tip_get(view_ctx, name,
		tip_length, tip_width);

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    return bsg_feature_arrow_tip_get(ref, tip_length, tip_width) ? 1 : 0;
}

int
ged_draw_view_context_arrow_tip_set(void *view_ctx,
				    const char *name,
				    fastf_t tip_length,
				    fastf_t tip_width)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_arrow_tip_set(view_ctx, name,
		tip_length, tip_width);

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    return bsg_feature_arrow_tip_set(ref, tip_length, tip_width) ? 1 : 0;
}

int
ged_draw_view_context_tcl_arrows_replace(
	void *view_ctx,
	const char *name,
	const point_t *points,
	size_t point_count,
	const struct ged_draw_view_feature_style *style)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_tcl_arrows_replace(view_ctx, name,
		points, point_count, style);

    if (point_count % 2)
	return 0;

    if (!points || !point_count) {
	if (!points && point_count)
	    return 0;
	bsg_feature_remove(view, name);
	return 1;
    }

    struct bsg_feature_style bsg_style = BSG_FEATURE_STYLE_INIT;
    const struct bsg_feature_style *bsg_stylep = NULL;
    if (style) {
	_ged_draw_view_feature_style_to_bsg(&bsg_style, style);
	bsg_stylep = &bsg_style;
    }

    bsg_feature_ref ref = bsg_feature_replace_arrow(view, name, 1 /* local */,
	    points, point_count, bsg_stylep);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    bsg_feature_overlay_register_owner(ref, NULL,
	    BSG_OVERLAY_ROLE_MODEL,
	    BSG_OVERLAY_CLASS_TCL_OVERLAY,
	    BSG_OVERLAY_LC_PER_COMMAND,
	    BSG_OVERLAY_ORDER_POST_TRANSPARENT,
	    NULL, 0);

    return 1;
}

int
ged_draw_view_context_feature_axes_centers_copy(
	void *view_ctx,
	const char *name,
	point_t **centers,
	size_t *center_count)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (centers)
	*centers = NULL;
    if (center_count)
	*center_count = 0;
    if (!view || !name || !centers || !center_count)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_feature_axes_centers_copy(view_ctx,
		name, centers, center_count);

    bsg_feature_ref ref = bsg_feature_find(view, name);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    return bsg_feature_axes_centers_copy(ref, centers, center_count) ? 1 : 0;
}

int
ged_draw_view_context_tcl_axes_replace(
	void *view_ctx,
	const char *name,
	const point_t *centers,
	size_t center_count,
	fastf_t half_axes_size,
	const struct ged_draw_view_feature_style *style)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name || !centers || !center_count)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_tcl_axes_replace(view_ctx, name,
		centers, center_count, half_axes_size, style);

    struct bsg_feature_style bsg_style = BSG_FEATURE_STYLE_INIT;
    const struct bsg_feature_style *bsg_stylep = NULL;
    if (style) {
	_ged_draw_view_feature_style_to_bsg(&bsg_style, style);
	bsg_stylep = &bsg_style;
    }

    bsg_feature_ref ref = bsg_feature_replace_axes(view, name, 1 /* local */,
	    centers, center_count, half_axes_size, bsg_stylep);
    if (bsg_feature_ref_is_null(ref))
	return 0;

    bsg_feature_overlay_register_owner(ref, NULL,
	    BSG_OVERLAY_ROLE_MODEL,
	    BSG_OVERLAY_CLASS_TCL_OVERLAY,
	    BSG_OVERLAY_LC_PER_COMMAND,
	    BSG_OVERLAY_ORDER_POST_TRANSPARENT,
	    NULL, 0);

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
ged_draw_view_context_axes_create(
	void *view_ctx,
	const char *name,
	int local,
	const struct ged_draw_view_axes_state *state)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name || !state)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_axes_create(view_ctx, name, local,
		state);

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
ged_draw_view_context_axes_state_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_axes_state *state)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name || !state)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_axes_state_get(view_ctx, name,
		state);

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
ged_draw_view_context_axes_state_replace(
	void *view_ctx,
	const char *name,
	const struct ged_draw_view_axes_state *state)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name || !state)
	return 0;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_axes_state_replace(view_ctx, name,
		state);

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

static rt_view_polygon_ref
_ged_draw_view_polygon_ref_to_rt(ged_draw_view_polygon_ref ref)
{
    rt_view_polygon_ref rt_ref = { ref.token, ref.revision };
    return rt_ref;
}

static ged_draw_view_polygon_ref
_ged_draw_view_polygon_ref_from_rt(rt_view_polygon_ref ref)
{
    ged_draw_view_polygon_ref ged_ref = { ref.token, ref.revision };
    return ged_ref;
}

static void
_ged_draw_view_polygon_record_from_rt(
	struct ged_draw_view_polygon_record *record,
	const struct rt_view_polygon_record *rt_record)
{
    if (!record || !rt_record)
	return;

    memset(record, 0, sizeof(*record));
    record->ref = _ged_draw_view_polygon_ref_from_rt(rt_record->ref);
    record->name = rt_record->name;
    record->type = rt_record->type;
    record->fill_flag = rt_record->fill_flag;
    V2MOVE(record->fill_dir, rt_record->fill_dir);
    record->fill_delta = rt_record->fill_delta;
    BU_COLOR_CPY(&record->fill_color, &rt_record->fill_color);
    record->edge_color[0] = rt_record->edge_color[0];
    record->edge_color[1] = rt_record->edge_color[1];
    record->edge_color[2] = rt_record->edge_color[2];
    record->curr_contour_i = rt_record->curr_contour_i;
    record->curr_point_i = rt_record->curr_point_i;
    record->first_contour_open = rt_record->first_contour_open;
    record->contour_count = rt_record->contour_count;
    record->point_count = rt_record->point_count;
    VMOVE(record->origin_point, rt_record->origin_point);
    HMOVE(record->vp, rt_record->vp);
    record->vZ = rt_record->vZ;
    record->user_data = rt_record->user_data;
}

static int
_ged_draw_view_polygon_edge_color(struct bu_color *edge_color,
	const struct ged_draw_view_polygon_record *record)
{
    if (!edge_color || !record)
	return 0;
    return bu_color_from_rgb_chars(edge_color, record->edge_color) ? 1 : 0;
}

int
ged_draw_view_polygon_ref_is_null(ged_draw_view_polygon_ref ref)
{
    return rt_view_polygon_ref_is_null(_ged_draw_view_polygon_ref_to_rt(ref));
}

ged_draw_view_polygon_ref
ged_draw_view_context_polygon_find(void *view_ctx, const char *name)
{
    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_polygon_find(view_ctx, name);

    return _ged_draw_view_polygon_ref_from_bsg(
	    bsg_view_polygon_find_ref((struct bsg_view *)view_ctx, name));
}

ged_draw_view_polygon_ref
ged_draw_view_context_polygon_find_scoped(void *view_ctx,
					  const char *name,
					  int local_only)
{
    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_polygon_find_scoped(view_ctx,
		name, local_only);

    return _ged_draw_view_polygon_ref_from_bsg(
	    bsg_view_polygon_find_scoped_ref((struct bsg_view *)view_ctx,
		name, local_only));
}

ged_draw_view_polygon_ref
ged_draw_view_context_polygon_create(void *view_ctx,
				     const char *name,
				     int local,
				     int type,
				     const point_t screen_point)
{
    point_t point;

    if (!screen_point)
	return GED_DRAW_VIEW_POLYGON_REF_NULL;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_polygon_create(view_ctx, name,
		local, type, screen_point);

    VMOVE(point, screen_point);
    return _ged_draw_view_polygon_ref_from_bsg(
	    bsg_create_polygon_ref((struct bsg_view *)view_ctx, name, local,
		type, &point));
}

ged_draw_view_polygon_ref
ged_draw_view_context_polygon_import_sketch(const char *name,
					    struct db_i *dbip,
					    struct directory *dp,
					    void *view_ctx,
					    int local)
{
    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	return ged_draw_obol_view_context_polygon_import_sketch(name,
		dbip, dp, view_ctx, local);

    return _ged_draw_view_polygon_ref_from_rt(
	    db_sketch_to_view_polygon_scoped_ref(name, dbip, dp, view_ctx,
		local));
}

int
ged_draw_view_polygon_export_sketch(struct db_i *dbip,
				    const char *name,
				    ged_draw_view_polygon_ref ref)
{
    return rt_view_polygon_export_sketch(dbip, name,
	    _ged_draw_view_polygon_ref_to_rt(ref)) ? 1 : 0;
}

int
ged_draw_view_polygon_record_get(ged_draw_view_polygon_ref ref,
				 struct ged_draw_view_polygon_record *record)
{
    struct rt_view_polygon_record rt_record;

    if (!record)
	return 0;

    if (!rt_view_polygon_record_get(_ged_draw_view_polygon_ref_to_rt(ref),
	    &rt_record))
	return 0;

    _ged_draw_view_polygon_record_from_rt(record, &rt_record);
    return 1;
}

int
ged_draw_view_polygon_has_data(ged_draw_view_polygon_ref ref)
{
    struct rt_view_polygon_record record;
    return rt_view_polygon_record_get(_ged_draw_view_polygon_ref_to_rt(ref),
	    &record) ? 1 : 0;
}

int
ged_draw_view_context_polygon_update(ged_draw_view_polygon_ref ref,
				     void *view_ctx,
				     int op)
{
    return rt_view_polygon_update_context(_ged_draw_view_polygon_ref_to_rt(ref),
	    view_ctx, op) ? 1 : 0;
}

int
ged_draw_view_context_polygon_update_screen_pt(ged_draw_view_polygon_ref ref,
					       void *view_ctx,
					       int x,
					       int y,
					       int op)
{
    return rt_view_polygon_update_screen_pt_context(
	    _ged_draw_view_polygon_ref_to_rt(ref), view_ctx, x, y, op) ?
	1 : 0;
}

int
ged_draw_view_polygon_set_current(ged_draw_view_polygon_ref ref,
				  long contour_i,
				  long point_i)
{
    if (ged_draw_obol_view_context_polygon_set_current(ref, contour_i,
		point_i))
	return 1;

    return bsg_polygon_set_current(_ged_draw_view_polygon_ref_to_bsg(ref),
	    contour_i, point_i) ? 1 : 0;
}

int
ged_draw_view_polygon_set_contour_open(ged_draw_view_polygon_ref ref,
				       long contour_i,
				       int open)
{
    if (ged_draw_obol_view_context_polygon_set_contour_open(ref,
		contour_i, open))
	return 1;

    return bsg_polygon_set_contour_open(_ged_draw_view_polygon_ref_to_bsg(ref),
	    contour_i, open) ? 1 : 0;
}

int
ged_draw_view_polygon_set_all_contours_open(ged_draw_view_polygon_ref ref,
					    int open)
{
    return rt_view_polygon_set_open(_ged_draw_view_polygon_ref_to_rt(ref),
	    open) ? 1 : 0;
}

int
ged_draw_view_context_polygon_area(ged_draw_view_polygon_ref ref,
				   void *view_ctx,
				   fastf_t *area)
{
    if (!area)
	return 0;
    *area = 0.0;

    if (!view_ctx)
	return 0;

    if (ged_draw_obol_view_context_polygon_area(ref, view_ctx, area))
	return 1;

    const struct bsg_polygon *poly =
	bsg_polygon_data(_ged_draw_view_polygon_ref_to_bsg(ref));
    if (!poly)
	return 0;

    *area = bg_find_polygon_area((struct bg_polygon *)&poly->polygon,
	    CLIPPER_MAX, (plane_t *)&poly->vp,
	    rt_view_context_scale_get(view_ctx));
    return 1;
}

int
ged_draw_view_context_polygon_overlap(ged_draw_view_polygon_ref ref,
				      void *view_ctx,
				      const char *other_name,
				      const struct bn_tol *tol,
				      int *overlap)
{
    if (!overlap)
	return 0;
    *overlap = 0;

    if (!view_ctx || !other_name || !tol)
	return 0;

    if (ged_draw_obol_view_context_polygon_overlap(ref, view_ctx,
		other_name, tol, overlap))
	return 1;

    const struct bsg_polygon *poly_a =
	bsg_polygon_data(_ged_draw_view_polygon_ref_to_bsg(ref));
    bsg_polygon_ref other_ref =
	bsg_view_polygon_find_ref((struct bsg_view *)view_ctx, other_name);
    const struct bsg_polygon *poly_b = bsg_polygon_data(other_ref);
    if (!poly_a || !poly_b)
	return 0;

    *overlap = bg_polygons_overlap((struct bg_polygon *)&poly_a->polygon,
	    (struct bg_polygon *)&poly_b->polygon, (plane_t *)&poly_a->vp,
	    tol, rt_view_context_scale_get(view_ctx));
    return 1;
}

int
ged_draw_view_polygon_set_fill(ged_draw_view_polygon_ref ref,
			       int fill_flag,
			       fastf_t fill_slope_x,
			       fastf_t fill_slope_y,
			       fastf_t fill_density)
{
    struct ged_draw_view_polygon_record record;
    struct bu_color edge_color;
    if (!ged_draw_view_polygon_record_get(ref, &record) ||
	    !_ged_draw_view_polygon_edge_color(&edge_color, &record))
	return 0;
    return rt_view_polygon_set_visual(_ged_draw_view_polygon_ref_to_rt(ref),
	    &edge_color, &record.fill_color, fill_slope_x, fill_slope_y,
	    fill_density, record.vZ, fill_flag) ? 1 : 0;
}

int
ged_draw_view_polygon_fill_color_get(ged_draw_view_polygon_ref ref,
				     struct bu_color *fill_color)
{
    struct ged_draw_view_polygon_record record;
    if (!fill_color || !ged_draw_view_polygon_record_get(ref, &record))
	return 0;
    BU_COLOR_CPY(fill_color, &record.fill_color);
    return 1;
}

int
ged_draw_view_polygon_fill_color_set(ged_draw_view_polygon_ref ref,
				     const struct bu_color *fill_color)
{
    struct ged_draw_view_polygon_record record;
    struct bu_color edge_color;
    if (!fill_color || !ged_draw_view_polygon_record_get(ref, &record) ||
	    !_ged_draw_view_polygon_edge_color(&edge_color, &record))
	return 0;
    return rt_view_polygon_set_visual(_ged_draw_view_polygon_ref_to_rt(ref),
	    &edge_color, fill_color, record.fill_dir[0],
	    record.fill_dir[1], record.fill_delta, record.vZ,
	    record.fill_flag) ? 1 : 0;
}

int
ged_draw_view_context_polygon_csg(ged_draw_view_polygon_ref target,
				  void *view_ctx,
				  const char *other_name,
				  bg_clip_t op)
{
    if (!view_ctx || !other_name)
	return 0;

    ged_draw_view_polygon_ref other_ref =
	ged_draw_view_context_polygon_find(view_ctx, other_name);
    if (ged_draw_view_polygon_ref_is_null(other_ref))
	return 0;

    return rt_view_polygon_csg(_ged_draw_view_polygon_ref_to_rt(target),
	    _ged_draw_view_polygon_ref_to_rt(other_ref), op) ? 1 : 0;
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
