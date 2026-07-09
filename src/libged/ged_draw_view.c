/*                    G E D _ D R A W _ V I E W . C
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
/** @file ged_draw_view.c
 *
 * GED view-feature wrappers backed exclusively by Obol/libbrlobol.
 *
 * Absent Obol view state is treated as a no-op/failure condition, not as
 * permission to synthesize alternate retained records.
 */

#include "common.h"

#include <string.h>

#include "bu/color.h"
#include "bu/malloc.h"
#include "bu/vls.h"
#include "./ged_draw_private.h"
#include "./ged_draw_view_private.h"

struct ged_draw_pick_record {
    struct bu_vls path;
    fastf_t hit_dist;
};

struct ged_draw_pick_result {
    struct ged_draw_pick_record *records;
    size_t count;
    size_t capacity;
};

static int
ged_draw_pick_result_reserve(struct ged_draw_pick_result *result,
			     size_t capacity)
{
    if (!result)
	return 0;
    if (capacity <= result->capacity)
	return 1;

    size_t new_capacity = result->capacity ? result->capacity : 4;
    while (new_capacity < capacity)
	new_capacity *= 2;

    result->records = (struct ged_draw_pick_record *)bu_realloc(
	    result->records, new_capacity * sizeof(struct ged_draw_pick_record),
	    "ged draw pick result records");
    for (size_t i = result->capacity; i < new_capacity; i++) {
	BU_VLS_INIT(&result->records[i].path);
	result->records[i].hit_dist = 0.0;
    }
    result->capacity = new_capacity;
    return 1;
}

struct ged_draw_pick_result *
ged_draw_pick_result_create(void)
{
    struct ged_draw_pick_result *result;
    BU_ALLOC(result, struct ged_draw_pick_result);
    result->records = NULL;
    result->count = 0;
    result->capacity = 0;
    return result;
}

void
ged_draw_pick_result_free(struct ged_draw_pick_result *result)
{
    if (!result)
	return;
    for (size_t i = 0; i < result->capacity; i++)
	bu_vls_free(&result->records[i].path);
    bu_free(result->records, "ged draw pick result records");
    bu_free(result, "ged draw pick result");
}

size_t
ged_draw_pick_result_count(const struct ged_draw_pick_result *result)
{
    return result ? result->count : 0;
}

int
ged_draw_pick_result_path(const struct ged_draw_pick_result *result,
			  size_t index,
			  struct bu_vls *path_out)
{
    if (!result || !path_out || index >= result->count)
	return 0;
    bu_vls_strcpy(path_out, bu_vls_cstr(&result->records[index].path));
    return 1;
}

fastf_t
ged_draw_pick_result_hit_dist(const struct ged_draw_pick_result *result,
			      size_t index)
{
    if (!result || index >= result->count)
	return 0.0;
    return result->records[index].hit_dist;
}

int
ged_draw_pick_result_append_path(struct ged_draw_pick_result *result,
				 const char *path,
				 fastf_t hit_dist)
{
    if (!result || !path || !path[0])
	return 0;
    if (!ged_draw_pick_result_reserve(result, result->count + 1))
	return 0;

    struct ged_draw_pick_record *record = &result->records[result->count++];
    bu_vls_strcpy(&record->path, path);
    record->hit_dist = hit_dist;
    return 1;
}

int
ged_draw_pick_result_append_copy(struct ged_draw_pick_result *dest,
				 const struct ged_draw_pick_result *src,
				 size_t index,
				 fastf_t hit_dist)
{
    if (!dest || !src || index >= src->count)
	return 0;
    return ged_draw_pick_result_append_path(dest,
	    bu_vls_cstr(&src->records[index].path), hit_dist);
}

struct ged_draw_pick_result *
ged_draw_pick_result_filter_first(const struct ged_draw_pick_result *src)
{
    struct ged_draw_pick_result *result = ged_draw_pick_result_create();
    if (!result)
	return NULL;
    if (src && src->count > 0)
	(void)ged_draw_pick_result_append_copy(result, src, 0,
		src->records[0].hit_dist);
    return result;
}

struct ged_draw_pick_result *
ged_draw_view_context_pick_point(void *view_ctx, int x, int y, int first_only)
{
    return ged_draw_obol_view_context_pick_point(view_ctx, x, y, first_only);
}

struct ged_draw_pick_result *
ged_draw_view_context_pick_nearest(void *view_ctx, int x, int y)
{
    return ged_draw_obol_view_context_pick_nearest(view_ctx, x, y);
}

struct ged_draw_pick_result *
ged_draw_view_context_pick_rect(void *view_ctx, int x0, int y0, int x1, int y1)
{
    return ged_draw_obol_view_context_pick_rect(view_ctx, x0, y0, x1, y1);
}

int
ged_draw_view_context_hud_sync(void *UNUSED(view_ctx))
{
    return 1;
}


int
ged_draw_view_context_selection_available(void *view_ctx)
{
    return ged_draw_obol_view_context_selection_available(view_ctx);
}


size_t
ged_draw_view_context_selection_count(void *view_ctx)
{
    return ged_draw_obol_view_context_selection_count(view_ctx);
}


int
ged_draw_view_context_selection_path_foreach(
	void *view_ctx,
	ged_draw_view_context_selection_path_cb cb,
	void *data)
{
    return ged_draw_obol_view_context_selection_path_foreach(view_ctx, cb,
	    data);
}


int
ged_draw_view_context_selection_clear(void *view_ctx)
{
    return ged_draw_obol_view_context_selection_clear(view_ctx);
}

int
ged_draw_view_context_selection_set_pick_result(
	void *view_ctx,
	const struct ged_draw_pick_result *result,
	ged_draw_view_selection_path_cb cb,
	void *data)
{
    if (!view_ctx)
	return 0;
    if (!ged_draw_view_context_selection_clear(view_ctx))
	return 0;
    if (!result)
	return 1;

    int ok = 1;
    for (size_t i = 0; i < result->count; i++) {
	const char *path = bu_vls_cstr(&result->records[i].path);
	if (!path || !path[0])
	    continue;
	if (!ged_draw_view_context_selection_add_path(view_ctx,
		GED_DRAW_VIEW_SELECTION_SELECTED_PATH, path))
	    ok = 0;
	if (cb)
	    cb(path, data);
    }
    return ok;
}


int
ged_draw_view_context_selection_contains_path(
	void *view_ctx,
	enum ged_draw_view_selection_kind kind,
	const char *path)
{
    return ged_draw_obol_view_context_selection_contains_path(view_ctx,
	    (int)kind, path);
}


int
ged_draw_view_context_selection_add_path(
	void *view_ctx,
	enum ged_draw_view_selection_kind kind,
	const char *path)
{
    return ged_draw_obol_view_context_selection_add_path(view_ctx,
	    (int)kind, path);
}


int
ged_draw_view_context_selection_set_path(
	void *view_ctx,
	enum ged_draw_view_selection_kind kind,
	const char *path)
{
    return ged_draw_obol_view_context_selection_set_path(view_ctx,
	    (int)kind, path);
}


int
ged_draw_view_context_feature_exists(void *view_ctx, const char *name)
{
    return ged_draw_obol_view_context_feature_exists(view_ctx, name);
}


int
ged_draw_view_context_feature_remove(void *view_ctx, const char *name)
{
    return ged_draw_obol_view_context_feature_remove(view_ctx, name);
}


int
ged_draw_view_context_feature_summary(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_feature_summary *summary)
{
    return ged_draw_obol_view_context_feature_summary(view_ctx, name,
	    summary);
}


size_t
ged_draw_view_context_feature_metadata_count(void *view_ctx, const char *name)
{
    return ged_draw_obol_view_context_feature_metadata_count(view_ctx, name);
}


int
ged_draw_view_context_feature_metadata_copy(
	void *view_ctx,
	const char *name,
	size_t index,
	struct bu_vls *key,
	struct bu_vls *value)
{
    return ged_draw_obol_view_context_feature_metadata_copy(view_ctx, name,
	    index, key, value);
}


size_t
ged_draw_view_context_feature_primitive_metadata_count(
	void *view_ctx,
	const char *name,
	int primitive)
{
    return ged_draw_obol_view_context_feature_primitive_metadata_count(
	    view_ctx, name, primitive);
}


int
ged_draw_view_context_feature_primitive_metadata_copy(
	void *view_ctx,
	const char *name,
	int primitive,
	size_t index,
	struct bu_vls *key,
	struct bu_vls *value)
{
    return ged_draw_obol_view_context_feature_primitive_metadata_copy(
	    view_ctx, name, primitive, index, key, value);
}


int
ged_draw_view_context_feature_pick_primitive_resolve(
	void *view_ctx,
	const char *picked_feature_name,
	int picked_primitive,
	int select,
	int highlight,
	struct bu_vls *feature_name,
	int *feature_primitive)
{
    return ged_draw_obol_view_context_feature_pick_primitive_resolve(
	    view_ctx, picked_feature_name, picked_primitive, select, highlight,
	    feature_name, feature_primitive);
}


int
ged_draw_view_context_feature_selected_primitives_replace(
	void *view_ctx,
	const char *name,
	const int *primitives,
	size_t primitive_count)
{
    return ged_draw_obol_view_context_feature_selected_primitives_replace(
	    view_ctx, name, primitives, primitive_count);
}


int
ged_draw_view_context_feature_highlighted_primitives_replace(
	void *view_ctx,
	const char *name,
	const int *primitives,
	size_t primitive_count)
{
    return ged_draw_obol_view_context_feature_highlighted_primitives_replace(
	    view_ctx, name, primitives, primitive_count);
}


size_t
ged_draw_view_context_feature_selected_primitive_count(
	void *view_ctx,
	const char *name)
{
    return ged_draw_obol_view_context_feature_selected_primitive_count(
	    view_ctx, name);
}


size_t
ged_draw_view_context_feature_highlighted_primitive_count(
	void *view_ctx,
	const char *name)
{
    return ged_draw_obol_view_context_feature_highlighted_primitive_count(
	    view_ctx, name);
}


int
ged_draw_view_context_feature_selected_primitive_at(
	void *view_ctx,
	const char *name,
	size_t index,
	int *primitive)
{
    return ged_draw_obol_view_context_feature_selected_primitive_at(
	    view_ctx, name, index, primitive);
}


int
ged_draw_view_context_feature_highlighted_primitive_at(
	void *view_ctx,
	const char *name,
	size_t index,
	int *primitive)
{
    return ged_draw_obol_view_context_feature_highlighted_primitive_at(
	    view_ctx, name, index, primitive);
}

int
ged_draw_view_feature_ref_is_null(ged_draw_view_feature_ref ref)
{
    return ged_draw_obol_view_feature_ref_is_null(ref);
}


int
ged_draw_view_context_edit_preview_publish_event(
	void *view_ctx,
	ged_draw_view_feature_ref feature,
	enum ged_draw_view_edit_preview_event event,
	const char *source_path)
{
    return ged_draw_obol_view_context_edit_preview_publish_event(view_ctx,
	    feature, event, source_path);
}


ged_draw_view_feature_ref
ged_draw_view_context_feature_overlay_ensure(
	void *view_ctx,
	const char *name,
	const void *owner,
	void *preview_ctx,
	const struct ged_draw_view_edit_preview_callbacks *callbacks,
	const char *source_path)
{
    return ged_draw_obol_view_context_feature_overlay_ensure(view_ctx, name,
	    owner, preview_ctx, callbacks, source_path);
}


ged_draw_view_feature_ref
ged_draw_view_context_feature_label_ensure(void *view_ctx,
					   const char *name,
					   const void *owner)
{
    return ged_draw_obol_view_context_feature_label_ensure(view_ctx, name,
	    owner);
}


void
ged_draw_view_feature_set_context(ged_draw_view_feature_ref ref, void *view_ctx)
{
    (void)ged_draw_obol_view_feature_set_context(ref, view_ctx);
}


void
ged_draw_view_feature_set_visible(ged_draw_view_feature_ref ref, int visible)
{
    (void)ged_draw_obol_view_feature_set_visible(ref, visible);
}


void
ged_draw_view_feature_set_color(ged_draw_view_feature_ref ref,
				int r,
				int g,
				int b)
{
    (void)ged_draw_obol_view_feature_set_color(ref, r, g, b);
}


int
ged_draw_view_feature_touch(ged_draw_view_feature_ref ref)
{
    return ged_draw_obol_view_feature_touch(ref);
}


int
ged_draw_view_feature_labels_replace(
	ged_draw_view_feature_ref ref,
	const struct ged_draw_view_feature_label *labels,
	size_t label_count)
{
    return ged_draw_obol_view_feature_labels_replace(ref, labels,
	    label_count);
}


int
ged_draw_view_feature_points_replace(
	ged_draw_view_feature_ref ref,
	enum ged_draw_view_feature_family family,
	const point_t *points,
	const int *cmds,
	size_t point_count)
{
    return ged_draw_obol_view_feature_points_replace(ref, family, points,
	    cmds, point_count);
}

int
ged_draw_view_feature_edit_preview_replace(
	ged_draw_view_feature_ref ref,
	const char *source_path,
	const char *edit_intent_id,
	const char *edit_intent_role,
	const point_t *points,
	const int *cmds,
	size_t point_count,
	uint32_t source_revision,
	uint32_t inputs_revision)
{
    return ged_draw_obol_view_feature_edit_preview_replace(ref, source_path,
	    edit_intent_id, edit_intent_role, points, cmds, point_count,
	    source_revision, inputs_revision);
}


int
ged_draw_view_feature_clear_geometry(ged_draw_view_feature_ref ref)
{
    return ged_draw_obol_view_feature_clear_geometry(ref);
}


int
ged_draw_view_context_features_remove_prefix(void *view_ctx, const char *prefix)
{
    return (int)ged_draw_obol_view_context_features_remove_prefix(view_ctx,
	    prefix);
}


int
ged_draw_view_context_feature_visible(void *view_ctx, const char *name)
{
    return ged_draw_obol_view_context_feature_visible(view_ctx, name);
}


int
ged_draw_view_context_feature_visible_set(void *view_ctx, const char *name, int visible)
{
    return ged_draw_obol_view_context_feature_visible_set(view_ctx, name,
	    visible);
}


int
ged_draw_view_context_feature_depth(void *view_ctx,
				    const char *name,
				    int mode,
				    fastf_t *depth)
{
    return ged_draw_obol_view_context_feature_depth(view_ctx, name, mode,
	    depth);
}


int
ged_draw_view_context_feature_depth_foreach(
	void *view_ctx,
	int mode,
	ged_draw_view_feature_depth_cb cb,
	void *data)
{
    return ged_draw_obol_view_context_feature_depth_foreach(view_ctx, mode,
	    cb, data);
}


int
ged_draw_view_context_feature_style_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_feature_style *style)
{
    return ged_draw_obol_view_context_feature_style_get(view_ctx, name,
	    style);
}


int
ged_draw_view_context_feature_style_apply(
	void *view_ctx,
	const char *name,
	const struct ged_draw_view_feature_style *style,
	int recursive)
{
    return ged_draw_obol_view_context_feature_style_apply(view_ctx, name,
	    style, recursive);
}


int
ged_draw_view_context_feature_realize(void *view_ctx,
				      const char *name,
				      int recursive)
{
    return ged_draw_obol_view_context_feature_realize(view_ctx, name,
	    recursive);
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
    return ged_draw_obol_view_context_indexed_face_set_replace(view_ctx,
	    name, local, points, point_count, normals, normal_count, indices,
	    index_count, style);
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
    return ged_draw_obol_view_context_lines_replace(view_ctx, name, local,
	    points, cmds, point_count, style);
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
    return ged_draw_obol_view_context_tcl_polygons_replace(view_ctx, name,
	    points, cmds, point_count, style);
}


int
ged_draw_view_context_data_polygons_replace(
	void *view_ctx,
	const char *name,
	int draw,
	const point_t *points,
	const int *cmds,
	size_t point_count,
	const struct ged_draw_view_feature_style *style)
{
    if (!view_ctx || !name)
	return 0;

    if (!draw || !points || !cmds || !point_count) {
	(void)ged_draw_view_context_feature_remove(view_ctx, name);
	return 1;
    }

    return ged_draw_view_context_tcl_polygons_replace(view_ctx, name, points,
	    cmds, point_count, style);
}


int
ged_draw_view_context_line_layer_builder_replace(
	void *view_ctx,
	const char *name,
	int local,
	const struct bg_line_layer_builder *builder)
{
    return ged_draw_obol_view_context_line_layer_builder_replace(view_ctx,
	    name, local, builder);
}


int
ged_draw_view_context_diagnostic_line_layer_builder_replace(
	void *view_ctx,
	const char *name,
	const struct bg_line_layer_builder *builder)
{
    return ged_draw_obol_view_context_diagnostic_line_layer_builder_replace(
	    view_ctx, name, builder);
}


int
ged_draw_view_context_line_layers_replace(void *view_ctx,
					  const char *name,
					  int local,
					  const struct ged_draw_view_line_layer_data *layers,
					  size_t layer_count,
					  const struct ged_draw_view_feature_style *style)
{
    return ged_draw_obol_view_context_line_layers_replace(view_ctx, name,
	    local, layers, layer_count, style);
}


int
ged_draw_view_context_lines_create_model_annotation(
	void *view_ctx,
	const char *name,
	int local,
	const point_t point)
{
    return ged_draw_obol_view_context_lines_create_model_annotation(
	    view_ctx, name, local, point);
}


int
ged_draw_view_context_lines_append_point(void *view_ctx,
					 const char *name,
					 const point_t point)
{
    return ged_draw_obol_view_context_lines_append_point(view_ctx, name,
	    point);
}


static int
ged_draw_view_context_annotation_line_replace(
	void *view_ctx,
	const char *name,
	int local,
	const point_t *points,
	size_t point_count,
	const struct ged_draw_view_feature_style *style)
{
    if (!view_ctx || !name)
	return 0;

    (void)ged_draw_view_context_feature_remove(view_ctx, name);
    if (!points || !point_count)
	return 1;

    if (!ged_draw_view_context_lines_create_model_annotation(view_ctx, name,
	    local, points[0]))
	return 0;

    for (size_t i = 1; i < point_count; i++) {
	if (!ged_draw_view_context_lines_append_point(view_ctx, name,
		points[i])) {
	    (void)ged_draw_view_context_feature_remove(view_ctx, name);
	    return 0;
	}
    }

    if (style && !ged_draw_view_context_feature_style_apply(view_ctx, name,
	    style, 0)) {
	(void)ged_draw_view_context_feature_remove(view_ctx, name);
	return 0;
    }

    return 1;
}


int
ged_draw_view_context_annotation_line_create(
	void *view_ctx,
	const char *name,
	int local,
	const point_t *points,
	size_t point_count,
	const struct ged_draw_view_feature_style *style,
	struct bu_vls *result)
{
    if (!view_ctx || !name || !points || !point_count)
	return 0;

    if (ged_draw_view_context_feature_exists(view_ctx, name)) {
	if (result)
	    bu_vls_printf(result, "View feature named %s already exists\n",
		    name);
	return 0;
    }

    if (!ged_draw_view_context_annotation_line_replace(view_ctx, name,
	    local, points, point_count, style)) {
	if (result)
	    bu_vls_printf(result, "Failed to create %s\n", name);
	return 0;
    }

    return 1;
}


int
ged_draw_view_context_annotation_line_append(
	void *view_ctx,
	const char *name,
	const point_t *points,
	size_t point_count,
	struct bu_vls *result)
{
    if (!view_ctx || !name || !points || !point_count)
	return 0;

    if (!ged_draw_view_context_feature_exists(view_ctx, name)) {
	if (result)
	    bu_vls_printf(result, "no view feature named %s\n", name);
	return 0;
    }

    for (size_t i = 0; i < point_count; i++) {
	if (!ged_draw_view_context_lines_append_point(view_ctx, name,
		points[i])) {
	    if (result)
		bu_vls_printf(result,
			"Failed to append point %zu to %s\n", i, name);
	    return 0;
	}
    }

    return 1;
}


int
ged_draw_view_context_annotation_line_remove_points(
	void *view_ctx,
	const char *name,
	size_t start,
	size_t count,
	struct bu_vls *result)
{
    point_t *points = NULL;
    size_t point_count = 0;

    if (!view_ctx || !name || !count)
	return 0;

    if (!ged_draw_view_context_feature_exists(view_ctx, name)) {
	if (result)
	    bu_vls_printf(result, "no view feature named %s\n", name);
	return 0;
    }

    if (!ged_draw_view_context_feature_points_copy(view_ctx, name, &points,
	    &point_count)) {
	if (result)
	    bu_vls_printf(result, "Failed to read points for %s\n", name);
	return 0;
    }

    if (start >= point_count || count > point_count - start) {
	if (result)
	    bu_vls_printf(result,
		    "Point range [%zu, %zu) is outside %s point count %zu\n",
		    start, start + count, name, point_count);
	if (points)
	    bu_free(points, "GED draw view annotation line copied points");
	return 0;
    }

    struct ged_draw_view_feature_style style =
	GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    int have_style = ged_draw_view_context_feature_style_get(view_ctx, name,
	    &style);
    size_t new_count = point_count - count;
    point_t *new_points = NULL;
    if (new_count) {
	new_points = (point_t *)bu_calloc(new_count, sizeof(point_t),
		"GED draw view annotation line rebuilt points");
	size_t out = 0;
	for (size_t i = 0; i < point_count; i++) {
	    if (i >= start && i < start + count)
		continue;
	    VMOVE(new_points[out], points[i]);
	    out++;
	}
    }

    int ret = ged_draw_view_context_annotation_line_replace(view_ctx, name,
	    1, (const point_t *)new_points, new_count,
	    have_style ? &style : NULL);

    if (new_points)
	bu_free(new_points, "GED draw view annotation line rebuilt points");
    if (points)
	bu_free(points, "GED draw view annotation line copied points");

    if (!ret && result)
	bu_vls_printf(result, "Failed to update %s\n", name);
    return ret ? 1 : 0;
}


int
ged_draw_view_context_annotation_line_clear(
	void *view_ctx,
	const char *name,
	struct bu_vls *result)
{
    if (!view_ctx || !name)
	return 0;

    if (!ged_draw_view_context_feature_exists(view_ctx, name)) {
	if (result)
	    bu_vls_printf(result, "no view feature named %s\n", name);
	return 0;
    }

    if (!ged_draw_view_context_feature_remove(view_ctx, name)) {
	if (result)
	    bu_vls_printf(result, "Failed to clear %s\n", name);
	return 0;
    }

    return 1;
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
    return ged_draw_obol_view_context_label_create(view_ctx, name, local,
	    text, point, target, has_target);
}


int
ged_draw_view_context_annotation_label_create(
	void *view_ctx,
	const char *name,
	int local,
	const char *text,
	const point_t point,
	const point_t target,
	int has_target,
	struct bu_vls *result)
{
    if (!view_ctx || !name || !text || !point || (has_target && !target))
	return 0;

    if (ged_draw_view_context_feature_exists(view_ctx, name)) {
	if (result)
	    bu_vls_printf(result, "View feature named %s already exists\n",
		    name);
	return 0;
    }

    if (!ged_draw_view_context_label_create(view_ctx, name, local, text,
	    point, target, has_target)) {
	if (result)
	    bu_vls_printf(result, "Failed to create %s\n", name);
	return 0;
    }

    return 1;
}


int
ged_draw_view_context_labels_replace(
	void *view_ctx,
	const char *name,
	int local,
	const struct ged_draw_view_label_data *labels,
	size_t label_count)
{
    return ged_draw_obol_view_context_labels_replace(view_ctx, name, local,
	    labels, label_count);
}


int
ged_draw_view_context_tcl_labels_replace(
	void *view_ctx,
	const char *name,
	int draw,
	const struct ged_draw_view_label_data *labels,
	size_t label_count)
{
    return ged_draw_obol_view_context_tcl_labels_replace(view_ctx, name,
	    draw, labels, label_count);
}


size_t
ged_draw_view_context_label_count(void *view_ctx, const char *name)
{
    return ged_draw_obol_view_context_label_count(view_ctx, name);
}


int
ged_draw_view_context_label_copy(void *view_ctx,
				 const char *name,
				 size_t index,
				 struct bu_vls *text,
				 point_t point,
				 unsigned char rgb[3])
{
    return ged_draw_obol_view_context_label_copy(view_ctx, name, index,
	    text, point, rgb);
}


int
ged_draw_view_context_label_point_set(void *view_ctx,
				      const char *name,
				      size_t index,
				      const point_t point)
{
    return ged_draw_obol_view_context_label_point_set(view_ctx, name, index,
	    point);
}


int
ged_draw_view_context_line_style_get(void *view_ctx,
				     const char *name,
				     struct ged_draw_view_line_style *style)
{
    return ged_draw_obol_view_context_line_style_get(view_ctx, name, style);
}


int
ged_draw_view_context_line_color_set(void *view_ctx,
				     const char *name,
				     int r,
				     int g,
				     int b)
{
    return ged_draw_obol_view_context_line_color_set(view_ctx, name, r, g,
	    b);
}


int
ged_draw_view_context_line_width_set(void *view_ctx,
				     const char *name,
				     int line_width)
{
    return ged_draw_obol_view_context_line_width_set(view_ctx, name,
	    line_width);
}


int
ged_draw_view_context_feature_points_copy(void *view_ctx,
					  const char *name,
					  point_t **points,
					  size_t *point_count)
{
    return ged_draw_obol_view_context_feature_points_copy(view_ctx, name,
	    points, point_count);
}


int
ged_draw_view_context_feature_line_command_at(
	void *view_ctx,
	const char *name,
	size_t index,
	int *out)
{
    return ged_draw_obol_view_context_feature_line_command_at(view_ctx, name,
	    index, out);
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
    return ged_draw_obol_view_context_tcl_lines_replace(view_ctx, name,
	    points, point_count, style);
}


int
ged_draw_view_context_data_lines_draw_get(
	void *view_ctx,
	const char *name)
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_feature_visible(view_ctx, name);
}


int
ged_draw_view_context_data_lines_draw_set(
	void *view_ctx,
	const char *name,
	int draw)
{
    if (!view_ctx || !name)
	return 0;

    if (!draw)
	(void)ged_draw_view_context_feature_remove(view_ctx, name);

    return 1;
}


int
ged_draw_view_context_data_lines_style_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_line_style *style)
{
    if (!view_ctx || !name || !style)
	return 0;
    return ged_draw_view_context_line_style_get(view_ctx, name, style);
}


int
ged_draw_view_context_data_lines_color_set(
	void *view_ctx,
	const char *name,
	int r,
	int g,
	int b)
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_line_color_set(view_ctx, name, r, g, b);
}


int
ged_draw_view_context_data_lines_line_width_set(
	void *view_ctx,
	const char *name,
	int line_width)
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_line_width_set(view_ctx, name, line_width);
}


int
ged_draw_view_context_data_lines_points_copy(
	void *view_ctx,
	const char *name,
	point_t **points,
	size_t *point_count)
{
    if (!view_ctx || !name || !points || !point_count)
	return 0;
    return ged_draw_view_context_lines_points_copy(view_ctx, name, points,
	    point_count);
}


int
ged_draw_view_context_data_lines_points_replace(
	void *view_ctx,
	const char *name,
	const point_t *points,
	size_t point_count,
	const struct ged_draw_view_line_style *style)
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_tcl_lines_replace(view_ctx, name, points,
	    point_count, style);
}


int
ged_draw_view_context_data_labels_draw_get(
	void *view_ctx,
	const char *name)
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_feature_exists(view_ctx, name);
}


int
ged_draw_view_context_data_labels_replace(
	void *view_ctx,
	const char *name,
	int draw,
	const struct ged_draw_view_label_data *labels,
	size_t label_count)
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_tcl_labels_replace(view_ctx, name, draw,
	    labels, label_count);
}


int
ged_draw_view_context_data_labels_color_get(
	void *view_ctx,
	const char *name,
	unsigned char rgb[3])
{
    if (rgb) {
	rgb[0] = 0;
	rgb[1] = 0;
	rgb[2] = 0;
    }
    if (!view_ctx || !name || !rgb)
	return 0;
    return ged_draw_view_context_label_copy(view_ctx, name, 0, NULL, NULL,
	    rgb);
}


size_t
ged_draw_view_context_data_labels_count(
	void *view_ctx,
	const char *name)
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_label_count(view_ctx, name);
}


int
ged_draw_view_context_data_labels_copy(
	void *view_ctx,
	const char *name,
	size_t index,
	struct bu_vls *text,
	point_t point,
	unsigned char rgb[3])
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_label_copy(view_ctx, name, index, text,
	    point, rgb);
}


int
ged_draw_view_context_data_labels_point_set(
	void *view_ctx,
	const char *name,
	size_t index,
	const point_t point)
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_label_point_set(view_ctx, name, index, point);
}


static void
ged_draw_data_arrow_style_default(struct ged_draw_view_feature_style *style)
{
    if (!style)
	return;
    struct ged_draw_view_feature_style init =
	GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    *style = init;
    style->visible = 1;
    style->color_valid = 1;
    style->color[0] = 255;
    style->color[1] = 255;
    style->color[2] = 0;
    style->line_width = 0;
    style->arrow = 1;
}


int
ged_draw_view_context_data_arrows_draw_get(
	void *view_ctx,
	const char *name)
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_feature_exists(view_ctx, name);
}


int
ged_draw_view_context_data_arrows_draw_set(
	void *view_ctx,
	const char *name,
	int draw)
{
    if (!view_ctx || !name)
	return 0;
    (void)ged_draw_view_context_feature_visible_set(view_ctx, name,
	    draw ? 1 : 0);
    return 1;
}


int
ged_draw_view_context_data_arrows_style_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_feature_style *style)
{
    if (!style)
	return 0;
    ged_draw_data_arrow_style_default(style);
    if (!view_ctx || !name)
	return 0;

    struct ged_draw_view_feature_style current =
	GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    if (!ged_draw_view_context_feature_style_get(view_ctx, name, &current))
	return 1;

    style->visible = current.visible;
    if (current.color_valid) {
	style->color_valid = 1;
	style->color[0] = current.color[0];
	style->color[1] = current.color[1];
	style->color[2] = current.color[2];
    }
    style->line_width = current.line_width;
    style->arrow_tip_length = current.arrow_tip_length;
    style->arrow_tip_width = current.arrow_tip_width;
    return 1;
}


int
ged_draw_view_context_data_arrows_color_set(
	void *view_ctx,
	const char *name,
	int r,
	int g,
	int b)
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_line_color_set(view_ctx, name, r, g, b);
}


int
ged_draw_view_context_data_arrows_line_width_set(
	void *view_ctx,
	const char *name,
	int line_width)
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_line_width_set(view_ctx, name, line_width);
}


int
ged_draw_view_context_data_arrows_points_copy(
	void *view_ctx,
	const char *name,
	point_t **points,
	size_t *point_count)
{
    if (!view_ctx || !name || !points || !point_count)
	return 0;
    return ged_draw_view_context_feature_points_copy(view_ctx, name, points,
	    point_count);
}


int
ged_draw_view_context_data_arrows_points_replace(
	void *view_ctx,
	const char *name,
	point_t *points,
	size_t point_count,
	const struct ged_draw_view_feature_style *style)
{
    if (!view_ctx || !name)
	return 0;
    if (point_count < 2) {
	(void)ged_draw_view_context_feature_remove(view_ctx, name);
	return 1;
    }
    return ged_draw_view_context_tcl_arrows_replace(view_ctx, name, points,
	    point_count, style);
}


int
ged_draw_view_context_data_arrows_tip_get(
	void *view_ctx,
	const char *name,
	fastf_t *tip_length,
	fastf_t *tip_width)
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_arrow_tip_get(view_ctx, name, tip_length,
	    tip_width);
}


int
ged_draw_view_context_data_arrows_tip_set(
	void *view_ctx,
	const char *name,
	fastf_t tip_length,
	fastf_t tip_width)
{
    if (!view_ctx || !name)
	return 0;
    if (!ged_draw_view_context_arrow_tip_get(view_ctx, name, NULL, NULL))
	return 1;
    return ged_draw_view_context_arrow_tip_set(view_ctx, name, tip_length,
	    tip_width);
}


int
ged_draw_view_context_arrow_tip_get(void *view_ctx,
				    const char *name,
				    fastf_t *tip_length,
				    fastf_t *tip_width)
{
    return ged_draw_obol_view_context_arrow_tip_get(view_ctx, name,
	    tip_length, tip_width);
}


int
ged_draw_view_context_arrow_tip_set(void *view_ctx,
				    const char *name,
				    fastf_t tip_length,
				    fastf_t tip_width)
{
    return ged_draw_obol_view_context_arrow_tip_set(view_ctx, name,
	    tip_length, tip_width);
}


int
ged_draw_view_context_tcl_arrows_replace(
	void *view_ctx,
	const char *name,
	point_t *points,
	size_t point_count,
	const struct ged_draw_view_feature_style *style)
{
    return ged_draw_obol_view_context_tcl_arrows_replace(view_ctx, name,
	    points, point_count, style);
}


int
ged_draw_view_context_feature_axes_centers_copy(
	void *view_ctx,
	const char *name,
	point_t **centers,
	size_t *center_count)
{
    return ged_draw_obol_view_context_feature_axes_centers_copy(view_ctx,
	    name, centers, center_count);
}


int
ged_draw_view_context_tcl_axes_replace(
	void *view_ctx,
	const char *name,
	point_t *centers,
	size_t center_count,
	fastf_t half_axes_size,
	const struct ged_draw_view_feature_style *style)
{
    return ged_draw_obol_view_context_tcl_axes_replace(view_ctx, name,
	    centers, center_count, half_axes_size, style);
}


static void
ged_draw_data_axes_style_default(struct ged_draw_view_feature_style *style)
{
    if (!style)
	return;
    struct ged_draw_view_feature_style init =
	GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    *style = init;
    style->visible = 1;
    style->color_valid = 1;
    style->color[0] = 255;
    style->color[1] = 255;
    style->color[2] = 0;
    style->line_width = 0;
}


int
ged_draw_view_context_data_axes_draw_get(
	void *view_ctx,
	const char *name)
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_feature_exists(view_ctx, name);
}


int
ged_draw_view_context_data_axes_draw_set(
	void *view_ctx,
	const char *name,
	int draw)
{
    if (!view_ctx || !name)
	return 0;
    (void)ged_draw_view_context_feature_visible_set(view_ctx, name,
	    draw ? 1 : 0);
    return 1;
}


int
ged_draw_view_context_data_axes_style_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_feature_style *style)
{
    if (!style)
	return 0;
    ged_draw_data_axes_style_default(style);
    if (!view_ctx || !name)
	return 0;

    struct ged_draw_view_feature_style current =
	GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    if (!ged_draw_view_context_feature_style_get(view_ctx, name, &current))
	return 1;

    style->visible = current.visible;
    if (current.color_valid) {
	style->color_valid = 1;
	style->color[0] = current.color[0];
	style->color[1] = current.color[1];
	style->color[2] = current.color[2];
    }
    style->line_width = current.line_width;
    return 1;
}


int
ged_draw_view_context_data_axes_color_set(
	void *view_ctx,
	const char *name,
	int r,
	int g,
	int b)
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_line_color_set(view_ctx, name, r, g, b);
}


int
ged_draw_view_context_data_axes_line_width_set(
	void *view_ctx,
	const char *name,
	int line_width)
{
    if (!view_ctx || !name)
	return 0;
    return ged_draw_view_context_line_width_set(view_ctx, name, line_width);
}


int
ged_draw_view_context_data_axes_half_size_get(
	void *view_ctx,
	const char *name,
	fastf_t *half_axes_size)
{
    if (half_axes_size)
	*half_axes_size = 1.0;
    if (!view_ctx || !name || !half_axes_size)
	return 0;

    point_t *points = NULL;
    size_t point_count = 0;
    int copied = ged_draw_view_context_feature_points_copy(view_ctx, name,
	    &points, &point_count);
    if (copied && point_count >= 2)
	*half_axes_size = (points[1][X] - points[0][X]) * 0.5;
    if (points)
	bu_free(points, "GED draw view data axes points copy");
    return 1;
}


int
ged_draw_view_context_data_axes_size_get(
	void *view_ctx,
	const char *name,
	fastf_t display_scale,
	fastf_t *size)
{
    if (size)
	*size = 0.0;
    if (!view_ctx || !name || !size)
	return 0;

    point_t *points = NULL;
    size_t point_count = 0;
    int copied = ged_draw_view_context_feature_points_copy(view_ctx, name,
	    &points, &point_count);
    if (copied && point_count >= 2) {
	fastf_t half = (points[1][X] - points[0][X]) * 0.5;
	*size = (display_scale > 0.0) ? (half * 2.0 / display_scale) : 0.0;
    }
    if (points)
	bu_free(points, "GED draw view data axes points copy");
    return 1;
}


int
ged_draw_view_context_data_axes_centers_copy(
	void *view_ctx,
	const char *name,
	point_t **centers,
	size_t *center_count)
{
    if (!view_ctx || !name || !centers || !center_count)
	return 0;
    return ged_draw_view_context_feature_axes_centers_copy(view_ctx, name,
	    centers, center_count);
}


int
ged_draw_view_context_data_axes_centers_replace(
	void *view_ctx,
	const char *name,
	point_t *centers,
	size_t center_count,
	fastf_t half_axes_size,
	const struct ged_draw_view_feature_style *style)
{
    if (!view_ctx || !name)
	return 0;
    if (!centers || !center_count) {
	(void)ged_draw_view_context_feature_remove(view_ctx, name);
	return 1;
    }
    return ged_draw_view_context_tcl_axes_replace(view_ctx, name, centers,
	    center_count, half_axes_size, style);
}


int
ged_draw_view_context_axes_create(
	void *view_ctx,
	const char *name,
	int local,
	const struct ged_draw_view_axes_state *state)
{
    return ged_draw_obol_view_context_axes_create(view_ctx, name, local,
	    state);
}


int
ged_draw_view_context_axes_state_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_axes_state *state)
{
    return ged_draw_obol_view_context_axes_state_get(view_ctx, name, state);
}


int
ged_draw_view_context_axes_state_replace(
	void *view_ctx,
	const char *name,
	const struct ged_draw_view_axes_state *state)
{
    return ged_draw_obol_view_context_axes_state_replace(view_ctx, name,
	    state);
}


int
ged_draw_view_context_annotation_axes_create(
	void *view_ctx,
	const char *name,
	int local,
	const struct ged_draw_view_axes_state *state,
	struct bu_vls *result)
{
    if (!view_ctx || !name || !state)
	return 0;

    if (ged_draw_view_context_feature_exists(view_ctx, name)) {
	if (result)
	    bu_vls_printf(result, "View feature named %s already exists\n",
		    name);
	return 0;
    }

    if (!ged_draw_view_context_axes_create(view_ctx, name, local, state)) {
	if (result)
	    bu_vls_printf(result, "Failed to set axes state for %s\n",
		    name);
	return 0;
    }

    return 1;
}


int
ged_draw_view_context_annotation_axes_state_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_axes_state *state,
	struct bu_vls *result)
{
    if (!view_ctx || !name || !state)
	return 0;

    if (!ged_draw_view_context_feature_exists(view_ctx, name)) {
	if (result)
	    bu_vls_printf(result, "View feature named %s does not exist\n",
		    name);
	return 0;
    }

    if (!ged_draw_view_context_axes_state_get(view_ctx, name, state)) {
	if (result)
	    bu_vls_printf(result, "View feature %s has no axes state\n",
		    name);
	return 0;
    }

    return 1;
}


int
ged_draw_view_context_annotation_axes_state_replace(
	void *view_ctx,
	const char *name,
	const struct ged_draw_view_axes_state *state,
	struct bu_vls *result)
{
    if (!view_ctx || !name || !state)
	return 0;

    if (!ged_draw_view_context_feature_exists(view_ctx, name)) {
	if (result)
	    bu_vls_printf(result, "View feature named %s does not exist\n",
		    name);
	return 0;
    }

    if (!ged_draw_view_context_axes_state_replace(view_ctx, name, state)) {
	if (result)
	    bu_vls_printf(result, "Failed to set axes state for %s\n",
		    name);
	return 0;
    }

    return 1;
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
    return ref.token ? 0 : 1;
}


ged_draw_view_polygon_ref
ged_draw_view_context_polygon_find(void *view_ctx, const char *name)
{
    return ged_draw_obol_view_context_polygon_find(view_ctx, name);
}


ged_draw_view_polygon_ref
ged_draw_view_context_polygon_find_scoped(void *view_ctx,
					  const char *name,
					  int local_only)
{
    return ged_draw_obol_view_context_polygon_find_scoped(view_ctx, name,
	    local_only);
}


ged_draw_view_polygon_ref
ged_draw_view_context_polygon_select(void *view_ctx,
				     const point_t model_point)
{
    return ged_draw_obol_view_context_polygon_select(view_ctx, model_point);
}


ged_draw_view_polygon_ref
ged_draw_view_context_polygon_create(void *view_ctx,
				     const char *name,
				     int local,
				     int type,
				     const point_t screen_point)
{
    return ged_draw_obol_view_context_polygon_create(view_ctx, name, local,
	    type, screen_point);
}


ged_draw_view_polygon_ref
ged_draw_view_context_polygon_dup(void *view_ctx,
				  const char *name,
				  const char *new_name)
{
    return ged_draw_obol_view_context_polygon_dup(view_ctx, name, new_name);
}


ged_draw_view_polygon_ref
ged_draw_view_context_polygon_import_sketch(const char *name,
					    struct db_i *dbip,
					    struct directory *dp,
					    void *view_ctx,
					    int local)
{
    return ged_draw_obol_view_context_polygon_import_sketch(name, dbip, dp,
	    view_ctx, local);
}


void
ged_draw_view_context_polygon_visit_records(
	void *view_ctx,
	ged_draw_view_polygon_record_cb callback,
	void *data)
{
    ged_draw_obol_view_context_polygon_visit_records(view_ctx, callback,
	    data);
}


size_t
ged_draw_view_context_polygon_snap_count(void *view_ctx,
	ged_draw_view_polygon_ref exclude)
{
    return ged_draw_obol_view_context_polygon_snap_count(view_ctx, exclude);
}


int
ged_draw_view_context_polygon_clear_point_selection(void *view_ctx)
{
    return ged_draw_obol_view_context_polygon_clear_point_selection(view_ctx);
}


int
ged_draw_view_context_polygon_snap_exclude_set(
	void *view_ctx,
	ged_draw_view_polygon_ref ref)
{
    return ged_draw_obol_view_context_polygon_snap_exclude_set(view_ctx, ref);
}


struct directory *
ged_draw_view_polygon_export_sketch(struct db_i *dbip,
				    const char *name,
				    ged_draw_view_polygon_ref ref)
{
    return ged_draw_obol_view_polygon_export_sketch(dbip, name, ref);
}


int
ged_draw_view_polygon_record_get(ged_draw_view_polygon_ref ref,
				 struct ged_draw_view_polygon_record *record)
{
    if (!record)
	return 0;

    return ged_draw_obol_view_polygon_record_get(ref, record);
}


int
ged_draw_view_polygon_has_data(ged_draw_view_polygon_ref ref)
{
    struct ged_draw_view_polygon_record ged_record;
    return ged_draw_view_polygon_record_get(ref, &ged_record) ? 1 : 0;
}


int
ged_draw_view_context_polygon_update(ged_draw_view_polygon_ref ref,
				     void *view_ctx,
				     int op)
{
    return ged_draw_obol_view_context_polygon_update(ref, view_ctx, op);
}


int
ged_draw_view_context_polygon_update_screen_pt(ged_draw_view_polygon_ref ref,
					       void *view_ctx,
					       int x,
					       int y,
					       int op)
{
    return ged_draw_obol_view_context_polygon_update_screen_pt(ref,
	    view_ctx, x, y, op);
}


int
ged_draw_view_polygon_move(ged_draw_view_polygon_ref ref,
			   point_t *current_point,
			   point_t *previous_point)
{
    return ged_draw_obol_view_polygon_move(ref, current_point,
	    previous_point);
}


int
ged_draw_view_polygon_set_name(ged_draw_view_polygon_ref ref,
			       const char *name)
{
    return ged_draw_obol_view_polygon_set_name(ref, name);
}


int
ged_draw_view_polygon_set_visual(ged_draw_view_polygon_ref ref,
				 const struct bu_color *edge_color,
				 const struct bu_color *fill_color,
				 fastf_t fill_slope_x,
				 fastf_t fill_slope_y,
				 fastf_t fill_density,
				 fastf_t vZ,
				 int fill_flag)
{
    return ged_draw_obol_view_polygon_set_visual(ref, edge_color, fill_color,
	    fill_slope_x, fill_slope_y, fill_density, vZ, fill_flag);
}


int
ged_draw_view_polygon_set_current(ged_draw_view_polygon_ref ref,
				  long contour_i,
				  long point_i)
{
    return ged_draw_obol_view_context_polygon_set_current(ref, contour_i,
	    point_i);
}


int
ged_draw_view_polygon_set_contour_open(ged_draw_view_polygon_ref ref,
				       long contour_i,
				       int open)
{
    return ged_draw_obol_view_context_polygon_set_contour_open(ref,
	    contour_i, open);
}


int
ged_draw_view_polygon_set_all_contours_open(ged_draw_view_polygon_ref ref,
					    int open)
{
    return ged_draw_obol_view_polygon_set_all_contours_open(ref, open);
}


int
ged_draw_view_polygon_close(ged_draw_view_polygon_ref ref)
{
    return ged_draw_obol_view_polygon_close(ref);
}


int
ged_draw_view_polygon_clear_selected_point(ged_draw_view_polygon_ref ref)
{
    return ged_draw_obol_view_polygon_clear_selected_point(ref);
}


int
ged_draw_view_polygon_remove(ged_draw_view_polygon_ref ref)
{
    return ged_draw_obol_view_polygon_remove(ref);
}


void *
ged_draw_view_polygon_user_data(ged_draw_view_polygon_ref ref)
{
    return ged_draw_obol_view_polygon_user_data(ref);
}


int
ged_draw_view_polygon_user_data_set(ged_draw_view_polygon_ref ref,
				    void *user_data)
{
    return ged_draw_obol_view_polygon_user_data_set(ref, user_data);
}


int
ged_draw_view_context_polygon_area(ged_draw_view_polygon_ref ref,
				   void *view_ctx,
				   fastf_t *area)
{
    return ged_draw_obol_view_context_polygon_area(ref, view_ctx, area);
}


int
ged_draw_view_context_polygon_overlap(ged_draw_view_polygon_ref ref,
				      void *view_ctx,
				      const char *other_name,
				      const struct bn_tol *tol,
				      int *overlap)
{
    return ged_draw_obol_view_context_polygon_overlap(ref, view_ctx,
	    other_name, tol, overlap);
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
    return ged_draw_view_polygon_set_visual(ref, &edge_color,
	    &record.fill_color, fill_slope_x, fill_slope_y, fill_density,
	    record.vZ, fill_flag);
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
    return ged_draw_view_polygon_set_visual(ref, &edge_color, fill_color,
	    record.fill_dir[0], record.fill_dir[1], record.fill_delta,
	    record.vZ, record.fill_flag);
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

    return ged_draw_view_polygon_csg(target, other_ref, op);
}


int
ged_draw_view_polygon_csg(ged_draw_view_polygon_ref target,
			  ged_draw_view_polygon_ref stencil,
			  bg_clip_t op)
{
    return ged_draw_obol_view_polygon_csg(target, stencil, op);
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
