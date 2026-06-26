/*            B S G _ G E D _ D R A W _ R E C O R D S . C
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
/** @file bsg_ged_draw_records.c
 *
 * Semantic GED draw-record construction and iteration.
 *
 * This file is intentionally record/ref-facing.  It may resolve the temporary
 * draw registry to retained nodes internally, but it must not grow new raw
 * node-storage dependencies.
 */

#include "common.h"

#include <string.h>

#include "bu/hash.h"
#include "bu/malloc.h"
#include "bu/sort.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "bsg/appearance.h"
#include "bsg/backend_scene.h"
#include "bsg/draw_intent.h"
#include "bsg/draw_source.h"
#include "bsg/export.h"
#include "bsg/render.h"
#include "bsg/scene_builder.h"
#include "ged/draw.h"
#include "rt/view.h"
#include "./ged_private.h"
#include "./bsg_ged_draw_private.h"


enum ged_draw_view_export_geometry_kind {
    GED_DRAW_VIEW_EXPORT_GEOMETRY_NONE = 0,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_MESH,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_LINE_SET,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_POINT_SET,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_FACE_SET,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_TEXT,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_IMAGE,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_OVERLAY,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_HUD,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_ANNOTATION,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_CSG_PROXY,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_BREP_PROXY,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_EDIT_PREVIEW,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_EXTERNAL_PROXY
};

#define GED_DRAW_VIEW_LINE_POINT_DRAW 12

struct ged_draw_view_export_detail {
    struct bu_vls path;
    struct ged_draw_view_db_object_record record;
    enum ged_draw_view_export_geometry_kind geometry_kind;
    struct {
	point_t *points;
	int *commands;
	int *indices;
	size_t point_count;
	size_t command_count;
	size_t index_count;
    } arrays;
    struct {
	point_t *points;
	int *indices;
	size_t point_count;
	size_t normal_count;
	size_t index_count;
	size_t face_count;
	int normals_per_index;
	int material_valid;
	int material_draw_mode;
	fastf_t material_transparency;
	int material_highlighted;
	unsigned char material_color[3];
    } surface;
    struct {
	int display_space;
	point_t *points;
	size_t point_count;
	size_t segment_count;
	int line_segment_valid;
	int line_start;
	int line_end;
	int text_segment_valid;
	int text_ref_point;
	char *text;
    } annotation;
};


static void
_ged_draw_view_export_detail_init(struct ged_draw_view_export_detail *detail)
{
    if (!detail)
	return;
    memset(detail, 0, sizeof(*detail));
    BU_VLS_INIT(&detail->path);
    MAT_IDN(detail->record.model_mat);
}


static void
_ged_draw_view_export_detail_free(struct ged_draw_view_export_detail *detail)
{
    if (!detail)
	return;
    bu_vls_free(&detail->path);
    if (detail->arrays.points)
	bu_free(detail->arrays.points, "ged draw view export array points");
    if (detail->arrays.commands)
	bu_free(detail->arrays.commands, "ged draw view export array commands");
    if (detail->arrays.indices)
	bu_free(detail->arrays.indices, "ged draw view export array indices");
    if (detail->surface.points)
	bu_free(detail->surface.points, "ged draw view export surface points");
    if (detail->surface.indices)
	bu_free(detail->surface.indices, "ged draw view export surface indices");
    if (detail->annotation.points)
	bu_free(detail->annotation.points,
		"ged draw view export annotation points");
    if (detail->annotation.text)
	bu_free(detail->annotation.text, "ged draw view export annotation text");
}


static const struct ged_draw_view_export_detail *
_ged_draw_view_export_detail_from_record(
	const struct ged_draw_view_db_object_record *rec)
{
    if (!rec || !rec->detail_token)
	return NULL;
    return (const struct ged_draw_view_export_detail *)rec->detail_token;
}


static int
_ged_draw_view_export_points_copy(point_t **dst,
				  const point_t *src,
				  size_t count,
				  const char *label)
{
    if (!dst)
	return 0;
    *dst = NULL;
    if (!count)
	return 1;
    if (!src)
	return 1;

    *dst = (point_t *)bu_calloc(count, sizeof(point_t), label);
    for (size_t i = 0; i < count; i++)
	VMOVE((*dst)[i], src[i]);
    return 1;
}


static int
_ged_draw_view_export_ints_copy(int **dst,
				const int *src,
				size_t count,
				const char *label)
{
    if (!dst)
	return 0;
    *dst = NULL;
    if (!count)
	return 1;
    if (!src)
	return 1;

    *dst = (int *)bu_calloc(count, sizeof(int), label);
    memcpy(*dst, src, count * sizeof(int));
    return 1;
}


static bsg_scene_ref
_draw_scene_for_shape_ref(struct ged *gedp, ged_draw_shape_ref ref)
{
    return ged_draw_registry_shape_scene_ref(gedp, ref);
}


static bsg_scene_ref
_draw_scene_for_group_ref(struct ged *gedp, ged_draw_group_ref ref)
{
    return ged_draw_registry_group_scene_ref(gedp, ref);
}


static int
_draw_foreach_shape_scene(bsg_scene_ref ref,
			  int (*cb)(bsg_scene_ref, void *),
			  void *userdata);

static const char *
_dbpath_skip_lead_slash(const char *s)
{
    if (s && *s == '/')
	return s + 1;
    return s;
}


static int
_ged_draw_shape_ref_selection_path(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   struct bu_vls *path_out)
{
    if (path_out)
	bu_vls_trunc(path_out, 0);
    if (!gedp || ged_draw_shape_ref_is_null(ref) || !path_out)
	return 0;

    bsg_scene_ref shape_ref = _draw_scene_for_shape_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;

    char *path = NULL;
    const char *source_path = NULL;
    ged_draw_shape_state *shape_data = ged_draw_scene_ref_shape_state(shape_ref);
    if (shape_data && shape_data->s_fullpath.fp_len > 0) {
	path = db_path_to_string(&shape_data->s_fullpath);
	source_path = path;
    } else {
	const struct db_full_path *fp = ged_draw_scene_ref_fullpath(shape_ref);
	if (fp) {
	    path = db_path_to_string(fp);
	    source_path = path;
	}
    }
    if ((!source_path || !source_path[0]) && shape_data &&
	    shape_data->display_name && shape_data->display_name[0])
	source_path = shape_data->display_name;
    const char *scene_name = ged_draw_scene_ref_name(shape_ref);
    if ((!source_path || !source_path[0]) && scene_name)
	source_path = scene_name;

    if (source_path && source_path[0])
	bu_vls_strcpy(path_out, _dbpath_skip_lead_slash(source_path));
    if (path)
	bu_free(path, "ged draw shape selection path");

    return bu_vls_strlen(path_out) > 0 ? 1 : 0;
}


static void
_draw_path_normalize(struct bu_vls *out, const char *path)
{
    if (!out)
	return;
    bu_vls_trunc(out, 0);
    if (!path)
	return;

    path = ged_draw_dbpath_skip_lead_slash(path);
    bu_vls_strcpy(out, path);
    while (bu_vls_strlen(out) > 0 &&
	    bu_vls_addr(out)[bu_vls_strlen(out) - 1] == '/')
	bu_vls_trunc(out, bu_vls_strlen(out) - 1);
}


static int
_draw_path_equal(const char *a, const char *b)
{
    if (!a || !b)
	return 0;
    a = ged_draw_dbpath_skip_lead_slash(a);
    b = ged_draw_dbpath_skip_lead_slash(b);
    return BU_STR_EQUAL(a, b);
}


static int
_draw_path_is_prefix(const char *prefix, const char *path)
{
    if (!prefix || !path)
	return 0;
    prefix = ged_draw_dbpath_skip_lead_slash(prefix);
    path = ged_draw_dbpath_skip_lead_slash(path);
    size_t plen = strlen(prefix);
    if (plen == 0)
	return 1;
    if (bu_strncmp(prefix, path, plen) != 0)
	return 0;
    return path[plen] == '\0' || path[plen] == '/';
}


struct _draw_path_leaf_ctx {
    bu_hash_tbl *paths;
    size_t count;
};


static union tree *
_draw_path_leaf_cb(struct db_tree_state *UNUSED(tsp),
		   const struct db_full_path *pathp,
		   struct rt_db_internal *UNUSED(ip),
		   void *client_data)
{
    struct _draw_path_leaf_ctx *ctx =
	(struct _draw_path_leaf_ctx *)client_data;
    if (!ctx || !ctx->paths || !pathp)
	return TREE_NULL;

    char *path = db_path_to_string(pathp);
    if (!path)
	return TREE_NULL;
    const char *key = ged_draw_dbpath_skip_lead_slash(path);
    if (key && key[0] && bu_hash_set(ctx->paths, (const uint8_t *)key,
	    strlen(key), (void *)1) == 1)
	ctx->count++;
    bu_free(path, "draw path leaf string");
    return TREE_NULL;
}


static size_t
_draw_path_expected_leaf_paths(struct ged *gedp, const char *path,
			       bu_hash_tbl *expected)
{
    if (!gedp || !gedp->dbip || !path || !expected)
	return 0;

    struct db_tree_state state;
    db_init_db_tree_state(&state, gedp->dbip);
    struct _draw_path_leaf_ctx ctx;
    ctx.paths = expected;
    ctx.count = 0;

    const char *av[1] = {path};
    int ret = db_walk_tree(gedp->dbip, 1, av, 1, &state,
	    NULL, NULL, _draw_path_leaf_cb, (void *)&ctx);
    db_free_db_tree_state(&state);

    if (ret < 0)
	return 0;
    return ctx.count;
}


static bsg_scene_ref
_draw_group_scene_ref_of_shape(bsg_scene_ref shape_ref)
{
    if (ged_draw_scene_ref_is_null(shape_ref))
	return bsg_scene_ref_null();

    bsg_scene_ref group_ref = ged_draw_scene_ref_parent(shape_ref);
    while (!ged_draw_scene_ref_is_null(group_ref)) {
	bsg_scene_ref parent_ref = ged_draw_scene_ref_parent(group_ref);
	if (ged_draw_scene_ref_is_null(parent_ref))
	    return group_ref;
	bsg_scene_ref grandparent_ref = ged_draw_scene_ref_parent(parent_ref);
	if (ged_draw_scene_ref_is_null(grandparent_ref))
	    return group_ref;
	group_ref = parent_ref;
    }
    return bsg_scene_ref_null();
}


static int
_draw_scene_ref_in_view_scope(bsg_scene_ref ref)
{
    bsg_scene_ref cur = ref;

    while (!ged_draw_scene_ref_is_null(cur)) {
	if (ged_draw_scene_ref_is_view_scope(cur))
	    return 1;
	cur = ged_draw_scene_ref_parent(cur);
    }

    return 0;
}


static int
_draw_group_scene_ref_in_view(bsg_scene_ref group_ref, void *view_ctx)
{
    if (ged_draw_scene_ref_is_null(group_ref))
	return 0;

    int in_view_scope = _draw_scene_ref_in_view_scope(group_ref);
    if (!view_ctx)
	return !in_view_scope;
    if (rt_view_context_is_independent(view_ctx))
	return in_view_scope && ged_draw_scene_ref_view_context(group_ref) == view_ctx;
    if (!in_view_scope)
	return 1;
    return ged_draw_scene_ref_view_context(group_ref) == view_ctx;
}


int
ged_draw_group_record_in_view(const struct ged_draw_group_record *rec,
			      void *view_ctx)
{
    if (!rec)
	return 0;
    if (!view_ctx)
	return !rec->in_view_scope;
    if (rt_view_context_is_independent(view_ctx))
	return rec->in_view_scope && rec->view == view_ctx;
    if (!rec->in_view_scope)
	return 1;
    return rec->view == view_ctx;
}


struct _draw_path_state_ctx {
    struct ged *gedp;
    void *view_ctx;
    const char *path;
    int mode;
    bu_hash_tbl *drawn_leaf_paths;
    int exact_shape;
    int descendant_shape;
};


static int
_draw_path_state_shape_scene_cb(bsg_scene_ref shape_ref, void *ud)
{
    struct _draw_path_state_ctx *ctx =
	(struct _draw_path_state_ctx *)ud;
    if (!ctx || ged_draw_scene_ref_is_null(shape_ref) ||
	    !ged_draw_scene_ref_visible(shape_ref))
	return 1;
    if (ctx->mode >= 0 && ged_draw_scene_ref_draw_mode(shape_ref) != ctx->mode)
	return 1;

    bsg_scene_ref group_ref = _draw_group_scene_ref_of_shape(shape_ref);
    if (ged_draw_group_scene_ref_is_overlay(group_ref) ||
	    !ged_draw_scene_ref_visible(group_ref) ||
	    !_draw_group_scene_ref_in_view(group_ref, ctx->view_ctx))
	return 1;

    const struct db_full_path *fullpath = ged_draw_scene_ref_fullpath(shape_ref);
    if (!fullpath)
	return 1;

    char *path = db_path_to_string(fullpath);
    if (!path)
	return 1;
    const char *key = ged_draw_dbpath_skip_lead_slash(path);

    if (_draw_path_equal(ctx->path, key))
	ctx->exact_shape = 1;
    if (_draw_path_is_prefix(ctx->path, key))
	ctx->descendant_shape = 1;
    if (key && key[0])
	(void)bu_hash_set(ctx->drawn_leaf_paths,
		(const uint8_t *)key, strlen(key), (void *)1);

    bu_free(path, "draw shape fullpath string");
    return 1;
}


int
ged_draw_path_state(struct ged *gedp,
		    void *view_ctx,
		    const char *path,
		    int mode)
{
    if (!gedp || !path)
	return 0;

    struct bu_vls norm = BU_VLS_INIT_ZERO;
    _draw_path_normalize(&norm, path);
    if (bu_vls_strlen(&norm) == 0) {
	bu_vls_free(&norm);
	return 0;
    }

    bu_hash_tbl *drawn = bu_hash_create(0);
    if (!drawn) {
	bu_vls_free(&norm);
	return 0;
    }

    struct _draw_path_state_ctx ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.gedp = gedp;
    ctx.view_ctx = view_ctx;
    ctx.path = bu_vls_cstr(&norm);
    ctx.mode = mode;
    ctx.drawn_leaf_paths = drawn;

    ged_draw_scene_root_foreach_shape(gedp, 0,
	    _draw_path_state_shape_scene_cb, &ctx);

    if (ctx.exact_shape) {
	bu_hash_destroy(drawn);
	bu_vls_free(&norm);
	return 1;
    }

    if (!ctx.descendant_shape) {
	bu_hash_destroy(drawn);
	bu_vls_free(&norm);
	return 0;
    }

    bu_hash_tbl *expected = bu_hash_create(0);
    if (!expected) {
	bu_hash_destroy(drawn);
	bu_vls_free(&norm);
	return 0;
    }

    size_t expected_count =
	_draw_path_expected_leaf_paths(gedp, bu_vls_cstr(&norm), expected);

    int state = 0;
    if (expected_count > 0) {
	size_t matched = 0;
	bu_hash_entry *e = NULL;
	while ((e = bu_hash_next(expected, e)) != NULL) {
	    uint8_t *key = NULL;
	    size_t key_len = 0;
	    if (bu_hash_key(e, &key, &key_len) == 0 &&
		    bu_hash_get(drawn, key, key_len))
		matched++;
	}

	if (ctx.exact_shape || matched == expected_count)
	    state = 1;
	else if (matched > 0 || ctx.descendant_shape)
	    state = 2;
    } else {
	if (ctx.exact_shape)
	    state = 1;
	else if (ctx.descendant_shape)
	    state = 2;
    }

    bu_hash_destroy(expected);
    bu_hash_destroy(drawn);
    bu_vls_free(&norm);
    return state;
}


struct _draw_path_list_ctx {
    struct ged *gedp;
    void *view_ctx;
    int mode;
    bu_hash_tbl *seen;
    bu_hash_tbl *state_cache;
    char **paths;
    size_t count;
    size_t cap;
};


static void
_draw_path_list_free(struct _draw_path_list_ctx *ctx)
{
    if (!ctx)
	return;
    for (size_t i = 0; i < ctx->count; i++)
	bu_free(ctx->paths[i], "draw path list string");
    if (ctx->paths)
	bu_free(ctx->paths, "draw path list");
    ctx->paths = NULL;
    ctx->count = 0;
    ctx->cap = 0;
    if (ctx->seen)
	bu_hash_destroy(ctx->seen);
    ctx->seen = NULL;
    if (ctx->state_cache)
	bu_hash_destroy(ctx->state_cache);
    ctx->state_cache = NULL;
}


static int
_draw_path_list_add(struct _draw_path_list_ctx *ctx, const char *path)
{
    if (!ctx || !ctx->seen || !path || !*path)
	return 1;

    path = ged_draw_dbpath_skip_lead_slash(path);
    if (!path || !*path)
	return 1;

    size_t len = strlen(path);
    if (bu_hash_set(ctx->seen, (const uint8_t *)path, len, (void *)1) != 1)
	return 1;

    if (ctx->count == ctx->cap) {
	size_t ncap = ctx->cap ? ctx->cap * 2 : 32;
	char **npaths = (char **)bu_realloc(ctx->paths,
		ncap * sizeof(char *), "draw path list");
	ctx->paths = npaths;
	ctx->cap = ncap;
    }
    ctx->paths[ctx->count++] = bu_strdup(path);
    return 1;
}


static int
_draw_path_strcmp(const void *a, const void *b, void *UNUSED(data))
{
    const char * const *sa = (const char * const *)a;
    const char * const *sb = (const char * const *)b;
    return bu_strcmp(*sa, *sb);
}


static int
_draw_path_list_add_collapsed(struct _draw_path_list_ctx *ctx, const char *path)
{
    if (!ctx || !path || !*path)
	return 1;
    path = ged_draw_dbpath_skip_lead_slash(path);

    for (size_t i = 0; i < ctx->count; i++) {
	if (_draw_path_equal(ctx->paths[i], path) ||
		_draw_path_is_prefix(ctx->paths[i], path))
	    return 1;
    }

    size_t write = 0;
    for (size_t i = 0; i < ctx->count; i++) {
	if (_draw_path_is_prefix(path, ctx->paths[i])) {
	    bu_hash_rm(ctx->seen, (const uint8_t *)ctx->paths[i],
		    strlen(ctx->paths[i]));
	    bu_free(ctx->paths[i], "draw path list collapsed descendant");
	    continue;
	}
	ctx->paths[write++] = ctx->paths[i];
    }
    ctx->count = write;

    return _draw_path_list_add(ctx, path);
}


static int
_draw_path_list_cached_state(struct _draw_path_list_ctx *ctx, const char *path)
{
    if (!ctx || !ctx->state_cache || !path)
	return 0;

    path = ged_draw_dbpath_skip_lead_slash(path);
    if (!path || !*path)
	return 0;

    size_t len = strlen(path);
    void *cached = bu_hash_get(ctx->state_cache, (const uint8_t *)path, len);
    if (cached)
	return (int)((uintptr_t)cached) - 1;

    int state = ged_draw_path_state(ctx->gedp, ctx->view_ctx, path, ctx->mode);
    (void)bu_hash_set(ctx->state_cache, (const uint8_t *)path, len,
	    (void *)(uintptr_t)(state + 1));
    return state;
}


static int
_draw_list_shape_scene_path_collapsed_cb(bsg_scene_ref shape_ref, void *ud)
{
    struct _draw_path_list_ctx *ctx = (struct _draw_path_list_ctx *)ud;
    if (!ctx || ged_draw_scene_ref_is_null(shape_ref) ||
	    !ged_draw_scene_ref_visible(shape_ref))
	return 1;
    if (ctx->mode >= 0 && ged_draw_scene_ref_draw_mode(shape_ref) != ctx->mode)
	return 1;

    bsg_scene_ref group_ref = _draw_group_scene_ref_of_shape(shape_ref);
    if (ged_draw_group_scene_ref_is_overlay(group_ref) ||
	    !ged_draw_scene_ref_visible(group_ref) ||
	    !_draw_group_scene_ref_in_view(group_ref, ctx->view_ctx))
	return 1;

    const struct db_full_path *fullpath = ged_draw_scene_ref_fullpath(shape_ref);
    if (!fullpath)
	return 1;

    char *path = db_path_to_string(fullpath);
    if (!path)
	return 1;
    const char *spath = ged_draw_dbpath_skip_lead_slash(path);
    if (!spath || !*spath) {
	bu_free(path, "draw list collapsed shape path");
	return 1;
    }

    struct bu_vls candidate = BU_VLS_INIT_ZERO;
    const char *best = spath;
    size_t len = strlen(spath);
    for (size_t i = 0; i <= len; i++) {
	if (spath[i] != '/' && spath[i] != '\0')
	    continue;
	bu_vls_trunc(&candidate, 0);
	bu_vls_strncpy(&candidate, spath, i);
	if (_draw_path_list_cached_state(ctx, bu_vls_cstr(&candidate)) == 1) {
	    best = bu_vls_cstr(&candidate);
	    break;
	}
    }

    (void)_draw_path_list_add_collapsed(ctx, best);
    bu_vls_free(&candidate);
    bu_free(path, "draw list collapsed shape path");
    return 1;
}


static int
_draw_list_shape_scene_path_cb(bsg_scene_ref shape_ref, void *ud)
{
    struct _draw_path_list_ctx *ctx = (struct _draw_path_list_ctx *)ud;
    if (!ctx || ged_draw_scene_ref_is_null(shape_ref) ||
	    !ged_draw_scene_ref_visible(shape_ref))
	return 1;
    if (ctx->mode >= 0 && ged_draw_scene_ref_draw_mode(shape_ref) != ctx->mode)
	return 1;

    bsg_scene_ref group_ref = _draw_group_scene_ref_of_shape(shape_ref);
    if (ged_draw_group_scene_ref_is_overlay(group_ref) ||
	    !ged_draw_scene_ref_visible(group_ref) ||
	    !_draw_group_scene_ref_in_view(group_ref, ctx->view_ctx))
	return 1;

    const struct db_full_path *fullpath = ged_draw_scene_ref_fullpath(shape_ref);
    if (!fullpath)
	return 1;

    char *path = db_path_to_string(fullpath);
    if (!path)
	return 1;
    (void)_draw_path_list_add(ctx, path);
    bu_free(path, "draw list shape path");
    return 1;
}


struct _draw_has_paths_ctx {
    void *view_ctx;
    int mode;
    int found;
};


static int
_draw_has_paths_shape_scene_cb(bsg_scene_ref shape_ref, void *ud)
{
    struct _draw_has_paths_ctx *ctx = (struct _draw_has_paths_ctx *)ud;
    if (!ctx || ged_draw_scene_ref_is_null(shape_ref) ||
	    !ged_draw_scene_ref_visible(shape_ref))
	return 1;
    if (ctx->mode >= 0 && ged_draw_scene_ref_draw_mode(shape_ref) != ctx->mode)
	return 1;

    bsg_scene_ref group_ref = _draw_group_scene_ref_of_shape(shape_ref);
    if (ged_draw_group_scene_ref_is_overlay(group_ref) ||
	    !ged_draw_scene_ref_visible(group_ref) ||
	    !_draw_group_scene_ref_in_view(group_ref, ctx->view_ctx))
	return 1;

    if (!ged_draw_scene_ref_fullpath(shape_ref))
	return 1;

    ctx->found = 1;
    return 0;
}


size_t
ged_draw_list_paths(struct ged *gedp,
		    void *view_ctx,
		    int mode,
		    int expanded,
		    struct bu_vls *result)
{
    if (!gedp || !result)
	return 0;

    struct _draw_path_list_ctx ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.gedp = gedp;
    ctx.view_ctx = view_ctx;
    ctx.mode = mode;
    ctx.seen = bu_hash_create(0);
    ctx.state_cache = bu_hash_create(0);
    if (!ctx.seen || !ctx.state_cache) {
	_draw_path_list_free(&ctx);
	return 0;
    }

    if (expanded)
	ged_draw_scene_root_foreach_shape(gedp, 0,
		_draw_list_shape_scene_path_cb, &ctx);
    else
	ged_draw_scene_root_foreach_shape(gedp, 0,
		_draw_list_shape_scene_path_collapsed_cb, &ctx);

    if (ctx.count > 1)
	bu_sort(ctx.paths, ctx.count, sizeof(char *), _draw_path_strcmp, NULL);

    for (size_t i = 0; i < ctx.count; i++)
	bu_vls_printf(result, "%s\n", ctx.paths[i]);

    size_t count = ctx.count;
    _draw_path_list_free(&ctx);
    return count;
}


int
ged_draw_has_paths(struct ged *gedp,
		   void *view_ctx,
		   int mode)
{
    if (!gedp)
	return 0;

    struct _draw_has_paths_ctx ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.view_ctx = view_ctx;
    ctx.mode = mode;

    ged_draw_scene_root_foreach_shape(gedp, 0,
	    _draw_has_paths_shape_scene_cb, &ctx);
    return ctx.found;
}


static int
_draw_foreach_shape_scene(bsg_scene_ref ref,
			  int (*cb)(bsg_scene_ref, void *),
			  void *userdata)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 1;

    if (ged_draw_scene_ref_is_shape(ref)) {
	if (ged_draw_shape_state_get_scene_ref(ref))
	    return (*cb)(ref, userdata);
	return 1;
    }

    for (size_t i = 0; i < ged_draw_scene_ref_child_count(ref); i++) {
	bsg_scene_ref child_ref = ged_draw_scene_ref_child_at(ref, i);
	if (!_draw_foreach_shape_scene(child_ref, cb, userdata))
	    return 0;
    }
    return 1;
}


static void
_ged_draw_fill_shape_record(struct ged *gedp,
			    bsg_scene_ref shape_ref,
			    struct ged_draw_shape_record *out)
{
    memset(out, 0, sizeof(*out));
    out->ref = ged_draw_shape_ref_from_scene_ref(gedp, shape_ref);
    out->group = ged_draw_group_ref_of_shape(gedp, out->ref);
    ged_draw_shape_state *shape_data = ged_draw_scene_ref_shape_state(shape_ref);
    if (shape_data) {
	out->fullpath = &shape_data->s_fullpath;
	out->display_name = shape_data->display_name;
	out->path_hash = shape_data->path_hash;
	out->source_revision = shape_data->source_revision;
	out->inputs_revision = shape_data->inputs_revision;
	out->realized_source_revision = shape_data->realized_source_revision;
	out->realized_inputs_revision = shape_data->realized_inputs_revision;
	out->stale = (shape_data->stale_reason != GED_DRAW_STALE_NONE ||
		shape_data->source_revision != shape_data->realized_source_revision ||
		shape_data->inputs_revision != shape_data->realized_inputs_revision);
	if (shape_data->stale_reason != GED_DRAW_STALE_NONE)
	    out->stale_reason =
		ged_draw_stale_reason_name(shape_data->stale_reason);
	else if (shape_data->source_revision !=
		shape_data->realized_source_revision)
	    out->stale_reason =
		ged_draw_stale_reason_name(GED_DRAW_STALE_SOURCE_CHANGED);
	else if (shape_data->inputs_revision !=
		shape_data->realized_inputs_revision)
	    out->stale_reason =
		ged_draw_stale_reason_name(GED_DRAW_STALE_VIEW_INPUT_CHANGED);
	else
	    out->stale_reason =
		ged_draw_stale_reason_name(GED_DRAW_STALE_NONE);
	if (shape_data->leaf_dp)
	    out->leaf_name = shape_data->leaf_dp->d_namep;
	else if (shape_data->s_fullpath.fp_len > 0 &&
		shape_data->s_fullpath.fp_names[shape_data->s_fullpath.fp_len - 1])
	    out->leaf_name =
		shape_data->s_fullpath.fp_names[shape_data->s_fullpath.fp_len - 1]->d_namep;
    }
    out->visible = ged_draw_scene_ref_visible(shape_ref);
    out->highlighted = ged_draw_scene_ref_highlighted(shape_ref);
    out->selected = 0;
    out->evaluated_region = ged_draw_scene_ref_legacy_eval_flag(shape_ref);
    void *active_view_ctx = ged_draw_active_view_ctx(gedp);
    if (active_view_ctx) {
	struct bu_vls selected_path = BU_VLS_INIT_ZERO;
	if (_ged_draw_shape_ref_selection_path(gedp, out->ref, &selected_path))
	    out->selected = ged_draw_view_context_selection_contains_path(
		    active_view_ctx,
		    GED_DRAW_VIEW_SELECTION_SELECTED_PATH,
		    bu_vls_cstr(&selected_path));
	bu_vls_free(&selected_path);
    }
    struct ged_draw_source_state *d = ged_draw_scene_ref_source_data(shape_ref);
    if (!shape_data && d) {
	out->source_revision = d->source_revision;
	out->inputs_revision = d->inputs_revision;
	out->realized_source_revision = d->realized_source_revision;
	out->realized_inputs_revision = d->realized_inputs_revision;
	out->stale = (d->stale_reason != GED_DRAW_STALE_NONE ||
		d->source_revision != d->realized_source_revision ||
		d->inputs_revision != d->realized_inputs_revision);
	if (d->stale_reason != GED_DRAW_STALE_NONE)
	    out->stale_reason = ged_draw_stale_reason_name(d->stale_reason);
	else if (d->source_revision != d->realized_source_revision)
	    out->stale_reason =
		ged_draw_stale_reason_name(GED_DRAW_STALE_SOURCE_CHANGED);
	else if (d->inputs_revision != d->realized_inputs_revision)
	    out->stale_reason =
		ged_draw_stale_reason_name(GED_DRAW_STALE_VIEW_INPUT_CHANGED);
	else
	    out->stale_reason = ged_draw_stale_reason_name(GED_DRAW_STALE_NONE);
    } else if (!shape_data) {
	out->stale_reason = ged_draw_stale_reason_name(GED_DRAW_STALE_NONE);
    }
    struct ged_draw_database_source_summary source_summary;
    if (ged_draw_scene_ref_source_summary(ged_draw_shape_source_ref(shape_ref),
	    &source_summary) && source_summary.has_state) {
	out->source_revision = source_summary.source_revision;
	out->inputs_revision = source_summary.inputs_revision;
	out->realized_source_revision =
	    source_summary.realized_source_revision;
	out->realized_inputs_revision =
	    source_summary.realized_inputs_revision;
	out->stale = source_summary.stale;
	out->stale_reason = source_summary.stale_reason;
    }
    out->drawn_revision = ged_draw_scene_ref_drawn_revision(shape_ref);
    out->transparency = ged_draw_scene_ref_transparency(shape_ref);
    int draw_mode = ged_draw_scene_ref_draw_mode(shape_ref);
    out->draw_mode = draw_mode;
    out->line_width = ged_draw_scene_ref_line_width(shape_ref);
    ged_draw_scene_ref_draw_center(shape_ref, out->center);
}


int
ged_draw_shape_record_get(struct ged *gedp,
			  ged_draw_shape_ref ref,
			  struct ged_draw_shape_record *out)
{
    if (!out)
	return 0;
    bsg_scene_ref shape_ref = _draw_scene_for_shape_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;
    _ged_draw_fill_shape_record(gedp, shape_ref, out);
    return 1;
}


int
ged_draw_view_selection_set_highlighted_shape_ref(struct ged *gedp,
						  void *view_ctx,
						  ged_draw_shape_ref ref)
{
    if (!gedp || !view_ctx)
	return 0;

    if (!ged_draw_view_context_selection_set_path(view_ctx,
	    GED_DRAW_VIEW_SELECTION_ALL, NULL))
	return 0;
    if (ged_draw_shape_ref_is_null(ref))
	return 1;

    struct bu_vls highlighted_path = BU_VLS_INIT_ZERO;
    if (!_ged_draw_shape_ref_selection_path(gedp, ref, &highlighted_path)) {
	bu_vls_free(&highlighted_path);
	return 0;
    }

    int ret = ged_draw_view_context_selection_add_path(view_ctx,
	    GED_DRAW_VIEW_SELECTION_HIGHLIGHTED_REF,
	    bu_vls_cstr(&highlighted_path));
    bu_vls_free(&highlighted_path);
    return ret;
}

int
ged_draw_view_selection_add_shape_ref_context(
	struct ged *gedp,
	void *view_ctx,
	ged_draw_shape_ref ref,
	void **selection_view_ctx,
	struct bu_vls *path)
{
    if (selection_view_ctx)
	*selection_view_ctx = NULL;
    if (path)
	bu_vls_trunc(path, 0);

    if (!gedp || ged_draw_shape_ref_is_null(ref))
	return 0;

    void *target = view_ctx;
    if (!ged_draw_view_context_selection_available(target))
	target = ged_draw_active_view_ctx(gedp);
    if (!ged_draw_view_context_selection_available(target))
	return 0;

    struct bu_vls selected_path = BU_VLS_INIT_ZERO;
    if (!_ged_draw_shape_ref_selection_path(gedp, ref, &selected_path)) {
	bu_vls_free(&selected_path);
	return 0;
    }

    if (path && bu_vls_strlen(&selected_path) > 0)
	bu_vls_strcpy(path, bu_vls_cstr(&selected_path));

    int ret = ged_draw_view_context_selection_add_path(target,
	    GED_DRAW_VIEW_SELECTION_SELECTED_PATH,
	    bu_vls_cstr(&selected_path));
    bu_vls_free(&selected_path);

    if (selection_view_ctx)
	*selection_view_ctx = target;
    return ret;
}


static int
_ged_draw_count_shape_cb(bsg_scene_ref UNUSED(ref), void *ud)
{
    int *count = (int *)ud;
    (*count)++;
    return 1;
}


struct _first_shape_scene_ctx {
    bsg_scene_ref result;
};


static int
_ged_draw_first_shape_scene_cb(bsg_scene_ref ref, void *ud)
{
    struct _first_shape_scene_ctx *ctx = (struct _first_shape_scene_ctx *)ud;
    ctx->result = ref;
    return 0;
}


int
ged_draw_group_record_get(struct ged *gedp,
			  ged_draw_group_ref ref,
			  struct ged_draw_group_record *out)
{
    if (!out)
	return 0;
    bsg_scene_ref group_ref = _draw_scene_for_group_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(group_ref))
	return 0;
    memset(out, 0, sizeof(*out));
    out->ref = ref;
    out->path = ged_draw_group_scene_ref_path(group_ref);
    out->view = ged_draw_scene_ref_view_context(group_ref);
    out->draw_mode = ged_draw_group_scene_ref_mode(group_ref);
    struct _first_shape_scene_ctx first_shape;
    first_shape.result = bsg_scene_ref_null();
    _draw_foreach_shape_scene(group_ref, _ged_draw_first_shape_scene_cb,
	    &first_shape);
    out->transparency = ged_draw_scene_ref_transparency(first_shape.result);
    out->visible = ged_draw_scene_ref_visible(group_ref);
    out->is_overlay = ged_draw_group_scene_ref_is_overlay(group_ref);
    out->is_view_scope = ged_draw_scene_ref_is_view_scope(group_ref);
    out->in_view_scope = _draw_scene_ref_in_view_scope(group_ref);
    out->is_view_source = ged_draw_scene_ref_is_view_source(group_ref);
    out->is_local_source = ged_draw_scene_ref_is_local_source(group_ref);
    _draw_foreach_shape_scene(group_ref, _ged_draw_count_shape_cb,
	    &out->shape_count);
    return 1;
}


struct _foreach_group_record_ctx {
    struct ged *gedp;
    int (*cb)(const struct ged_draw_group_record *, void *);
    void *userdata;
};


static int
_foreach_group_record_cb(bsg_scene_ref group_ref, void *ud)
{
    struct _foreach_group_record_ctx *ctx =
	(struct _foreach_group_record_ctx *)ud;
    struct ged_draw_group_record rec;
    if (!ged_draw_group_record_get(ctx->gedp,
	    ged_draw_group_ref_from_scene_ref(ctx->gedp, group_ref), &rec))
	return 1;
    return ctx->cb(&rec, ctx->userdata);
}


void
ged_draw_foreach_group_record(struct ged *gedp,
			      int (*cb)(const struct ged_draw_group_record *rec,
				  void *userdata),
			      void *userdata)
{
    if (!gedp || !cb)
	return;
    struct _foreach_group_record_ctx ctx;
    ctx.gedp = gedp;
    ctx.cb = cb;
    ctx.userdata = userdata;
    (void)ged_draw_scene_root_foreach_group(gedp, _foreach_group_record_cb,
	    &ctx);
}


struct _foreach_shape_record_ctx {
    struct ged *gedp;
    int (*cb)(const struct ged_draw_shape_record *, void *);
    void *userdata;
};


static int
_foreach_shape_record_cb(bsg_scene_ref ref, void *ud)
{
    struct _foreach_shape_record_ctx *ctx =
	(struct _foreach_shape_record_ctx *)ud;
    struct ged_draw_shape_record rec;
    _ged_draw_fill_shape_record(ctx->gedp, ref, &rec);
    return ctx->cb(&rec, ctx->userdata);
}


void
ged_draw_foreach_shape_record(struct ged *gedp,
			      int (*cb)(const struct ged_draw_shape_record *rec,
				  void *userdata),
			      void *userdata)
{
    if (!gedp || !cb)
	return;
    struct _foreach_shape_record_ctx ctx;
    ctx.gedp = gedp;
    ctx.cb = cb;
    ctx.userdata = userdata;
    (void)ged_draw_scene_root_foreach_shape(gedp, 0,
	    _foreach_shape_record_cb, &ctx);
}


static enum ged_draw_view_export_geometry_kind
_ged_draw_view_export_geometry_kind_from_bsg(bsg_render_geometry_kind kind)
{
    switch (kind) {
	case BSG_RENDER_GEOMETRY_MESH:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_MESH;
	case BSG_RENDER_GEOMETRY_LINE_SET:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET;
	case BSG_RENDER_GEOMETRY_INDEXED_LINE_SET:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_LINE_SET;
	case BSG_RENDER_GEOMETRY_POINT_SET:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_POINT_SET;
	case BSG_RENDER_GEOMETRY_INDEXED_FACE_SET:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_FACE_SET;
	case BSG_RENDER_GEOMETRY_TEXT:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_TEXT;
	case BSG_RENDER_GEOMETRY_IMAGE:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_IMAGE;
	case BSG_RENDER_GEOMETRY_OVERLAY:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_OVERLAY;
	case BSG_RENDER_GEOMETRY_HUD:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_HUD;
	case BSG_RENDER_GEOMETRY_ANNOTATION:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_ANNOTATION;
	case BSG_RENDER_GEOMETRY_CSG_PROXY:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_CSG_PROXY;
	case BSG_RENDER_GEOMETRY_BREP_PROXY:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_BREP_PROXY;
	case BSG_RENDER_GEOMETRY_EDIT_PREVIEW:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_EDIT_PREVIEW;
	case BSG_RENDER_GEOMETRY_EXTERNAL_PROXY:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_EXTERNAL_PROXY;
	case BSG_RENDER_GEOMETRY_NONE:
	default:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_NONE;
    }
}


static const char *
_ged_draw_view_export_type_name_from_bsg(
	bsg_render_geometry_role geometry_role)
{
    switch (geometry_role) {
	case BSG_RENDER_GEOMETRY_ROLE_AXES_WIDGET:
	    return "axes";
	case BSG_RENDER_GEOMETRY_ROLE_LINE_SET:
	    return "line";
	case BSG_RENDER_GEOMETRY_ROLE_TEXT_LABEL:
	    return "label";
	case BSG_RENDER_GEOMETRY_ROLE_POLYGON_REGION:
	    return "polygon";
	case BSG_RENDER_GEOMETRY_ROLE_DATABASE_OBJECT:
	    return "gobj";
	default:
	    return "object";
    }
}


static const char *
_ged_draw_view_export_geometry_name_from_bsg(
	bsg_render_geometry_kind geometry_kind)
{
    switch (geometry_kind) {
	case BSG_RENDER_GEOMETRY_MESH:
	    return "mesh";
	case BSG_RENDER_GEOMETRY_LINE_SET:
	    return "line-set";
	case BSG_RENDER_GEOMETRY_INDEXED_LINE_SET:
	    return "indexed-line-set";
	case BSG_RENDER_GEOMETRY_POINT_SET:
	    return "point-set";
	case BSG_RENDER_GEOMETRY_INDEXED_FACE_SET:
	    return "indexed-face-set";
	case BSG_RENDER_GEOMETRY_TEXT:
	    return "text";
	case BSG_RENDER_GEOMETRY_IMAGE:
	    return "image";
	case BSG_RENDER_GEOMETRY_OVERLAY:
	    return "overlay";
	case BSG_RENDER_GEOMETRY_HUD:
	    return "hud";
	case BSG_RENDER_GEOMETRY_ANNOTATION:
	    return "annotation";
	case BSG_RENDER_GEOMETRY_CSG_PROXY:
	    return "csg-proxy";
	case BSG_RENDER_GEOMETRY_BREP_PROXY:
	    return "brep-proxy";
	case BSG_RENDER_GEOMETRY_EDIT_PREVIEW:
	    return "edit-preview";
	case BSG_RENDER_GEOMETRY_EXTERNAL_PROXY:
	    return "external-proxy";
	case BSG_RENDER_GEOMETRY_NONE:
	default:
	    return "none";
    }
}


static int
_ged_draw_view_line_command_from_bsg(int command)
{
    switch (command) {
	case BSG_GEOMETRY_LINE_MOVE:
	    return GED_DRAW_VIEW_LINE_MOVE;
	case BSG_GEOMETRY_LINE_DRAW:
	    return GED_DRAW_VIEW_LINE_DRAW;
	case BSG_GEOMETRY_POINT_DRAW:
	    return GED_DRAW_VIEW_LINE_POINT_DRAW;
	default:
	    return command;
    }
}


static int
_ged_draw_view_line_command_from_detail(
	const struct ged_draw_view_export_detail *detail,
	size_t index)
{
    if (detail && index < detail->arrays.command_count &&
	    detail->arrays.commands)
	return detail->arrays.commands[index];
    return (index % 2) ? GED_DRAW_VIEW_LINE_DRAW : GED_DRAW_VIEW_LINE_MOVE;
}


static int
_ged_draw_view_export_detail_arrays_from_bsg(
	struct ged_draw_view_export_detail *detail,
	const struct bsg_render_geometry *geometry)
{
    if (!detail || !geometry)
	return 0;

    detail->arrays.point_count = geometry->arrays.point_count;
    detail->arrays.command_count = geometry->arrays.command_count;
    detail->arrays.index_count = geometry->arrays.index_count;
    if (!_ged_draw_view_export_points_copy(&detail->arrays.points,
	    (const point_t *)geometry->arrays.points,
	    geometry->arrays.point_count,
	    "ged draw view export array points"))
	return 0;
    if (!_ged_draw_view_export_ints_copy(&detail->arrays.commands,
	    geometry->arrays.commands, geometry->arrays.command_count,
	    "ged draw view export array commands"))
	return 0;
    for (size_t i = 0; i < detail->arrays.command_count; i++)
	detail->arrays.commands[i] =
	    _ged_draw_view_line_command_from_bsg(detail->arrays.commands[i]);
    if (!_ged_draw_view_export_ints_copy(&detail->arrays.indices,
	    geometry->arrays.indices, geometry->arrays.index_count,
	    "ged draw view export array indices"))
	return 0;

    return 1;
}


static int
_ged_draw_view_export_detail_surface_from_bsg(
	struct ged_draw_view_export_detail *detail,
	const struct bsg_render_geometry *geometry)
{
    if (!detail || !geometry)
	return 0;

    detail->surface.point_count = geometry->surface.point_count;
    detail->surface.normal_count = geometry->surface.normal_count;
    detail->surface.index_count = geometry->surface.index_count;
    detail->surface.face_count = geometry->surface.face_count;
    detail->surface.normals_per_index =
	(geometry->surface.normal_kind == BSG_RENDER_SURFACE_NORMALS_PER_INDEX);
    detail->surface.material_valid = geometry->surface.material_valid;
    detail->surface.material_draw_mode =
	geometry->surface.material.draw_mode;
    detail->surface.material_transparency =
	geometry->surface.material.transparency;
    detail->surface.material_highlighted =
	geometry->surface.material.highlighted;
    detail->surface.material_color[0] =
	geometry->surface.material.color[0];
    detail->surface.material_color[1] =
	geometry->surface.material.color[1];
    detail->surface.material_color[2] =
	geometry->surface.material.color[2];

    if (!_ged_draw_view_export_points_copy(&detail->surface.points,
	    (const point_t *)geometry->surface.points,
	    geometry->surface.point_count,
	    "ged draw view export surface points"))
	return 0;
    if (!_ged_draw_view_export_ints_copy(&detail->surface.indices,
	    geometry->surface.indices, geometry->surface.index_count,
	    "ged draw view export surface indices"))
	return 0;

    return 1;
}


static int
_ged_draw_view_export_detail_annotation_from_bsg(
	struct ged_draw_view_export_detail *detail,
	const struct bsg_render_geometry *geometry)
{
    if (!detail || !geometry)
	return 0;

    detail->annotation.display_space =
	(geometry->annotation.space == BSG_ANNOTATION_SPACE_DISPLAY);
    detail->annotation.point_count = geometry->annotation.point_count;
    detail->annotation.segment_count = geometry->annotation.segment_count;
    if (!_ged_draw_view_export_points_copy(&detail->annotation.points,
	    (const point_t *)geometry->annotation.points,
	    geometry->annotation.point_count,
	    "ged draw view export annotation points"))
	return 0;

    const struct bsg_annotation_segment *segments =
	geometry->annotation.segments;
    if (!segments)
	return 1;

    for (size_t i = 0; i < geometry->annotation.segment_count; i++) {
	const struct bsg_annotation_segment *seg = &segments[i];
	if (!detail->annotation.line_segment_valid &&
		seg->kind == BSG_ANNOTATION_SEGMENT_LINE) {
	    detail->annotation.line_segment_valid = 1;
	    detail->annotation.line_start = seg->data.line.start;
	    detail->annotation.line_end = seg->data.line.end;
	}
	if (!detail->annotation.text_segment_valid &&
		seg->kind == BSG_ANNOTATION_SEGMENT_TEXT) {
	    detail->annotation.text_segment_valid = 1;
	    detail->annotation.text_ref_point = seg->data.text.ref_pt;
	    if (seg->data.text.text)
		detail->annotation.text = bu_strdup(seg->data.text.text);
	}
    }

    return 1;
}


static int
_ged_draw_view_export_detail_geometry_from_bsg(
	struct ged_draw_view_export_detail *detail,
	const struct bsg_export_record *export_rec)
{
    if (!detail || !export_rec)
	return 0;

    detail->geometry_kind =
	_ged_draw_view_export_geometry_kind_from_bsg(
	    export_rec->geometry.kind);
    switch (detail->geometry_kind) {
	case GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET:
	case GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_LINE_SET:
	case GED_DRAW_VIEW_EXPORT_GEOMETRY_POINT_SET:
	    return _ged_draw_view_export_detail_arrays_from_bsg(detail,
		    &export_rec->geometry);
	case GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_FACE_SET:
	    return _ged_draw_view_export_detail_surface_from_bsg(detail,
		    &export_rec->geometry);
	case GED_DRAW_VIEW_EXPORT_GEOMETRY_ANNOTATION:
	    return _ged_draw_view_export_detail_annotation_from_bsg(detail,
		    &export_rec->geometry);
	default:
	    return 1;
    }
}


static int
_ged_draw_view_export_detail_from_bsg(
	struct ged_draw_view_export_detail *detail,
	const struct bsg_export_record *export_rec)
{
    if (!detail || !export_rec)
	return 0;

    _ged_draw_view_export_detail_init(detail);
    bu_vls_strcpy(&detail->path, bu_vls_cstr(&export_rec->path));

    struct ged_draw_view_db_object_record *rec = &detail->record;
    rec->path = bu_vls_cstr(&detail->path);
    rec->type_name =
	_ged_draw_view_export_type_name_from_bsg(
	    export_rec->source.geometry_role);
    rec->geometry_name =
	_ged_draw_view_export_geometry_name_from_bsg(
	    export_rec->geometry.kind);
    rec->draw_mode = export_rec->draw_mode;
    rec->transparency = export_rec->transparency;
    rec->evaluated_region = export_rec->evaluated_region;
    rec->is_database_source =
	(export_rec->source.scope == BSG_RENDER_SOURCE_SCOPE_DATABASE);
    rec->non_database_source = export_rec->non_database_source;
    rec->is_database_intent =
	(export_rec->source.draw_intent == BSG_RENDER_DRAW_INTENT_DATABASE);
    rec->is_local_source =
	(export_rec->source.scope == BSG_RENDER_SOURCE_SCOPE_VIEW_LOCAL);
    rec->is_view_source =
	(export_rec->source.scope == BSG_RENDER_SOURCE_SCOPE_VIEW_SHARED ||
	 export_rec->source.scope == BSG_RENDER_SOURCE_SCOPE_VIEW_LOCAL ||
	 export_rec->source.draw_intent == BSG_RENDER_DRAW_INTENT_OVERLAY ||
	 export_rec->source.draw_intent == BSG_RENDER_DRAW_INTENT_HUD);
    rec->highlighted = export_rec->highlighted;
    rec->selected = export_rec->selected;
    rec->visible = export_rec->visible;
    rec->line_style = export_rec->line_style;
    rec->color[0] = export_rec->color[0];
    rec->color[1] = export_rec->color[1];
    rec->color[2] = export_rec->color[2];
    MAT_COPY(rec->model_mat, export_rec->model_mat);
    VMOVE(rec->bounds_center, export_rec->bounds_center);
    rec->bounds_radius = export_rec->bounds_radius;
    rec->has_bounds = export_rec->has_bounds;
    rec->vlist_structure_count = export_rec->vlist_structure_count;
    rec->vlist_point_count = export_rec->vlist_point_count;
    rec->cache_identity = export_rec->cache_identity;
    rec->source_identity = export_rec->source.source_id;
    rec->detail_token = (uintptr_t)detail;

    if (!_ged_draw_view_export_detail_geometry_from_bsg(detail, export_rec)) {
	_ged_draw_view_export_detail_free(detail);
	return 0;
    }

    return 1;
}


static void
_ged_draw_foreach_view_export_record(void *view_ctx,
				     unsigned int query_flags,
				     unsigned int render_flags,
				     const char *glob,
				     int draw_mode,
				     ged_draw_view_db_object_record_cb cb,
				     void *userdata)
{
    struct bsg_view *v = (struct bsg_view *)view_ctx;

    if (!v || !cb)
	return;

    struct bsg_export_request request;
    bsg_export_request_init(&request, v);
    request.query_flags = query_flags;
    request.render_flags = render_flags;
    request.glob = glob;
    if (draw_mode >= 0)
	request.draw_mode = draw_mode;

    struct bsg_export_result *result = bsg_export_query(&request);
    if (!result)
	return;

    for (size_t i = 0; i < bsg_export_result_count(result); i++) {
	const struct bsg_export_record *export_rec =
	    bsg_export_result_get(result, i);
	if (!export_rec)
	    continue;

	struct ged_draw_view_export_detail detail;
	if (!_ged_draw_view_export_detail_from_bsg(&detail, export_rec))
	    continue;

	int keep_going = cb(&detail.record, userdata);
	_ged_draw_view_export_detail_free(&detail);
	if (!keep_going)
	    break;
    }

    bsg_export_result_free(result);
}


void
ged_draw_foreach_view_record_query(
	void *view_ctx,
	const struct ged_draw_view_record_query *query,
	ged_draw_view_db_object_record_cb cb,
	void *userdata)
{
    if (!view_ctx || !query)
	return;

    unsigned int query_flags = 0;
    if (query->flags & GED_DRAW_VIEW_RECORD_QUERY_VIEW_OBJECTS)
	query_flags |= BSG_EXPORT_QUERY_VIEW_OBJECTS;
    if (query->flags & GED_DRAW_VIEW_RECORD_QUERY_DB_OBJECTS)
	query_flags |= BSG_EXPORT_QUERY_DB_OBJECTS;
    if (query->flags & GED_DRAW_VIEW_RECORD_QUERY_LOCAL_ONLY)
	query_flags |= BSG_EXPORT_QUERY_LOCAL_ONLY;

    _ged_draw_foreach_view_export_record(view_ctx, query_flags,
	    BSG_RENDER_FLAG_PAYLOAD_PREPARE, query->glob, query->draw_mode,
	    cb, userdata);
}


void
ged_draw_foreach_view_db_object_record(void *view_ctx,
				       ged_draw_view_db_object_record_cb cb,
				       void *userdata)
{
    _ged_draw_foreach_view_export_record(view_ctx, BSG_EXPORT_QUERY_DB_OBJECTS,
	    BSG_RENDER_FLAG_PAYLOAD_PREPARE, NULL, BSG_EXPORT_DRAW_MODE_ANY,
	    cb, userdata);
}


void
ged_draw_foreach_visible_view_db_object_record(void *view_ctx,
					      ged_draw_view_db_object_record_cb cb,
					      void *userdata)
{
    _ged_draw_foreach_view_export_record(view_ctx,
	    BSG_EXPORT_QUERY_VISIBLE_ONLY | BSG_EXPORT_QUERY_DB_OBJECTS,
	    BSG_RENDER_FLAG_VISIBLE_ONLY | BSG_RENDER_FLAG_PAYLOAD_PREPARE,
	    NULL,
	    BSG_EXPORT_DRAW_MODE_ANY, cb, userdata);
}


void
ged_draw_foreach_visible_view_db_object_record_mode(
	void *view_ctx,
	int draw_mode,
	ged_draw_view_db_object_record_cb cb,
	void *userdata)
{
    _ged_draw_foreach_view_export_record(view_ctx,
	    BSG_EXPORT_QUERY_VISIBLE_ONLY | BSG_EXPORT_QUERY_DB_OBJECTS,
	    BSG_RENDER_FLAG_VISIBLE_ONLY | BSG_RENDER_FLAG_PAYLOAD_PREPARE,
	    NULL, draw_mode, cb, userdata);
}


void
ged_draw_foreach_visible_view_record(void *view_ctx,
				     ged_draw_view_db_object_record_cb cb,
				     void *userdata)
{
    _ged_draw_foreach_view_export_record(view_ctx, BSG_EXPORT_QUERY_VISIBLE_ONLY,
	    BSG_RENDER_FLAG_VISIBLE_ONLY | BSG_RENDER_FLAG_PAYLOAD_PREPARE,
	    NULL, BSG_EXPORT_DRAW_MODE_ANY, cb, userdata);
}


int
ged_draw_view_db_object_record_foreach_segment(
	const struct ged_draw_view_db_object_record *rec,
	ged_draw_view_segment_cb cb,
	void *userdata)
{
    const struct ged_draw_view_export_detail *detail =
	_ged_draw_view_export_detail_from_record(rec);
    size_t count = 0;

    if (!detail || !cb ||
	    detail->geometry_kind == GED_DRAW_VIEW_EXPORT_GEOMETRY_NONE)
	return 0;

    point_t last = VINIT_ZERO;
    point_t fin = VINIT_ZERO;

    if (detail->geometry_kind == GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET &&
	    detail->arrays.point_count && detail->arrays.points) {
	for (size_t i = 0; i < detail->arrays.point_count; i++) {
	    int cmd = _ged_draw_view_line_command_from_detail(detail, i);
	    if (cmd == GED_DRAW_VIEW_LINE_MOVE) {
		VMOVE(last, detail->arrays.points[i]);
		continue;
	    }
	    if (cmd == GED_DRAW_VIEW_LINE_DRAW) {
		VMOVE(fin, detail->arrays.points[i]);
	    } else {
		continue;
	    }
	    if (!cb(last, fin, userdata))
		return count;
	    count++;
	    VMOVE(last, fin);
	}
	return count;
    }

    if (detail->geometry_kind ==
	    GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_LINE_SET &&
	    detail->arrays.point_count && detail->arrays.points &&
	    detail->arrays.index_count && detail->arrays.indices) {
	int have_last = 0;
	for (size_t i = 0; i < detail->arrays.index_count; i++) {
	    int idx = detail->arrays.indices[i];
	    if (idx < 0) {
		have_last = 0;
		continue;
	    }
	    if ((size_t)idx >= detail->arrays.point_count)
		continue;
	    if (!have_last) {
		VMOVE(last, detail->arrays.points[idx]);
		have_last = 1;
		continue;
	    }
	    VMOVE(fin, detail->arrays.points[idx]);
	    if (!cb(last, fin, userdata))
		return count;
	    count++;
	    VMOVE(last, fin);
	}
	return count;
    }

    if (detail->geometry_kind ==
	    GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_FACE_SET &&
	    detail->surface.point_count && detail->surface.points &&
	    detail->surface.index_count && detail->surface.indices) {
	int first_idx = -1;
	int prev_idx = -1;
	size_t face_vertices = 0;

	for (size_t i = 0; i < detail->surface.index_count; i++) {
	    int idx = detail->surface.indices[i];
	    if (idx < 0) {
		if (face_vertices > 1 && first_idx >= 0 && prev_idx >= 0 &&
			prev_idx != first_idx) {
		    if (!cb(detail->surface.points[prev_idx],
			    detail->surface.points[first_idx], userdata))
			return count;
		    count++;
		}
		first_idx = -1;
		prev_idx = -1;
		face_vertices = 0;
		continue;
	    }
	    if ((size_t)idx >= detail->surface.point_count)
		continue;
	    if (!face_vertices) {
		first_idx = idx;
		prev_idx = idx;
		face_vertices = 1;
		continue;
	    }
	    if (prev_idx >= 0) {
		if (!cb(detail->surface.points[prev_idx],
			detail->surface.points[idx], userdata))
		    return count;
		count++;
	    }
	    prev_idx = idx;
	    face_vertices++;
	}

	if (face_vertices > 1 && first_idx >= 0 && prev_idx >= 0 &&
		prev_idx != first_idx) {
	    if (!cb(detail->surface.points[prev_idx],
		    detail->surface.points[first_idx], userdata))
		return count;
	    count++;
	}
    }

    return count;
}


int
ged_draw_view_db_object_record_foreach_point(
	const struct ged_draw_view_db_object_record *rec,
	ged_draw_view_point_cb cb,
	void *userdata)
{
    const struct ged_draw_view_export_detail *detail =
	_ged_draw_view_export_detail_from_record(rec);
    size_t count = 0;

    if (!detail || !cb ||
	    detail->geometry_kind == GED_DRAW_VIEW_EXPORT_GEOMETRY_NONE)
	return 0;

    if (detail->geometry_kind ==
	    GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_FACE_SET &&
	    detail->surface.point_count && detail->surface.points) {
	for (size_t i = 0; i < detail->surface.point_count; i++) {
	    if (!cb(detail->surface.points[i], userdata))
		return count;
	    count++;
	}
	return count;
    }

    if ((detail->geometry_kind == GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET ||
		detail->geometry_kind ==
		GED_DRAW_VIEW_EXPORT_GEOMETRY_POINT_SET ||
		detail->geometry_kind ==
		GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_LINE_SET) &&
	    detail->arrays.point_count && detail->arrays.points) {
	for (size_t i = 0; i < detail->arrays.point_count; i++) {
	    if (!cb(detail->arrays.points[i], userdata))
		return count;
	    count++;
	}
	return count;
    }

    return 0;
}


static int
_ged_draw_view_record_has_segment_cb(const point_t UNUSED(a),
				     const point_t UNUSED(b),
				     void *data)
{
    int *seen = (int *)data;
    if (seen)
	*seen = 1;
    return 0;
}


int
ged_draw_view_db_object_record_has_segments(
	const struct ged_draw_view_db_object_record *rec)
{
    int seen = 0;
    (void)ged_draw_view_db_object_record_foreach_segment(rec,
	    _ged_draw_view_record_has_segment_cb, &seen);
    return seen;
}


static const char *
_ged_draw_view_geometry_command_description(int cmd)
{
    switch (cmd) {
	case GED_DRAW_VIEW_LINE_MOVE:
	    return "line move";
	case GED_DRAW_VIEW_LINE_DRAW:
	    return "line draw";
	case GED_DRAW_VIEW_LINE_POINT_DRAW:
	    return "point draw";
	default:
	    return "unknown";
    }
}


void
ged_draw_view_db_object_record_geometry_report(
	const struct ged_draw_view_db_object_record *rec,
	struct bu_vls *out)
{
    const struct ged_draw_view_export_detail *detail =
	_ged_draw_view_export_detail_from_record(rec);
    if (!detail || !out)
	return;

    if ((detail->geometry_kind == GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET ||
		detail->geometry_kind ==
		GED_DRAW_VIEW_EXPORT_GEOMETRY_POINT_SET) &&
	    detail->arrays.point_count && detail->arrays.points) {
	for (size_t i = 0; i < detail->arrays.point_count; i++) {
	    int cmd = (detail->geometry_kind ==
		    GED_DRAW_VIEW_EXPORT_GEOMETRY_POINT_SET) ?
		GED_DRAW_VIEW_LINE_POINT_DRAW :
		_ged_draw_view_line_command_from_detail(detail, i);
	    const fastf_t *pt = detail->arrays.points[i];
	    bu_vls_printf(out, "  %s (%g, %g, %g)\n",
		    _ged_draw_view_geometry_command_description(cmd),
		    V3ARGS(pt));
	}
	return;
    }

    if (detail->geometry_kind ==
	    GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_LINE_SET &&
	    detail->arrays.point_count && detail->arrays.points &&
	    detail->arrays.index_count && detail->arrays.indices) {
	int have_current = 0;
	for (size_t i = 0; i < detail->arrays.index_count; i++) {
	    int idx = detail->arrays.indices[i];
	    if (idx < 0) {
		have_current = 0;
		continue;
	    }
	    if ((size_t)idx >= detail->arrays.point_count)
		continue;
	    int cmd = have_current ? GED_DRAW_VIEW_LINE_DRAW :
		GED_DRAW_VIEW_LINE_MOVE;
	    const fastf_t *pt = detail->arrays.points[idx];
	    bu_vls_printf(out, "  %s (%g, %g, %g)\n",
		    _ged_draw_view_geometry_command_description(cmd),
		    V3ARGS(pt));
	    have_current = 1;
	}
	return;
    }

    if (detail->geometry_kind ==
	    GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_FACE_SET &&
	    detail->surface.point_count && detail->surface.points &&
	    detail->surface.index_count && detail->surface.indices) {
	size_t face_vertex = 0;
	for (size_t i = 0; i < detail->surface.index_count; i++) {
	    int idx = detail->surface.indices[i];
	    if (idx < 0) {
		face_vertex = 0;
		continue;
	    }
	    if ((size_t)idx >= detail->surface.point_count)
		continue;
	    const fastf_t *pt = detail->surface.points[idx];
	    bu_vls_printf(out, "  indexed face %zu (%g, %g, %g)\n",
		    face_vertex, V3ARGS(pt));
	    face_vertex++;
	}
    }
}


int
ged_draw_view_db_object_record_annotation_summary(
	const struct ged_draw_view_db_object_record *rec,
	size_t point_index,
	struct ged_draw_view_annotation_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    const struct ged_draw_view_export_detail *detail =
	_ged_draw_view_export_detail_from_record(rec);
    if (!detail ||
	    detail->geometry_kind != GED_DRAW_VIEW_EXPORT_GEOMETRY_ANNOTATION)
	return 0;

    out->valid = 1;
    out->display_space = detail->annotation.display_space;
    out->point_count = detail->annotation.point_count;
    out->segment_count = detail->annotation.segment_count;
    out->cache_identity = rec->cache_identity;
    out->source_identity = rec->source_identity;
    if (point_index < detail->annotation.point_count &&
	    detail->annotation.points) {
	out->has_point = 1;
	VMOVE(out->point, detail->annotation.points[point_index]);
    }
    out->line_segment_valid = detail->annotation.line_segment_valid;
    out->line_start = detail->annotation.line_start;
    out->line_end = detail->annotation.line_end;
    out->text_segment_valid = detail->annotation.text_segment_valid;
    out->text_ref_point = detail->annotation.text_ref_point;
    out->text = detail->annotation.text;

    return 1;
}


int
ged_draw_view_db_object_record_line_summary(
	const struct ged_draw_view_db_object_record *rec,
	struct ged_draw_view_line_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    const struct ged_draw_view_export_detail *detail =
	_ged_draw_view_export_detail_from_record(rec);
    if (!detail ||
	    detail->geometry_kind != GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET)
	return 0;

    out->valid = 1;
    out->point_count = detail->arrays.point_count;
    out->cache_identity = rec->cache_identity;
    out->source_identity = rec->source_identity;

    return 1;
}


int
ged_draw_view_db_object_record_line_point_at(
	const struct ged_draw_view_db_object_record *rec,
	size_t index,
	point_t out)
{
    const struct ged_draw_view_export_detail *detail =
	_ged_draw_view_export_detail_from_record(rec);
    if (!out || !detail ||
	    detail->geometry_kind != GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET ||
	    index >= detail->arrays.point_count || !detail->arrays.points)
	return 0;

    VMOVE(out, detail->arrays.points[index]);
    return 1;
}


int
ged_draw_view_db_object_record_line_command_at(
	const struct ged_draw_view_db_object_record *rec,
	size_t index,
	int *out)
{
    const struct ged_draw_view_export_detail *detail =
	_ged_draw_view_export_detail_from_record(rec);
    if (!out || !detail ||
	    detail->geometry_kind != GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET ||
	    index >= detail->arrays.point_count)
	return 0;

    *out = _ged_draw_view_line_command_from_detail(detail, index);
    return 1;
}


int
ged_draw_view_db_object_record_surface_summary(
	const struct ged_draw_view_db_object_record *rec,
	struct ged_draw_view_surface_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    const struct ged_draw_view_export_detail *detail =
	_ged_draw_view_export_detail_from_record(rec);
    if (!detail ||
	    detail->geometry_kind !=
	    GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_FACE_SET)
	return 0;

    out->valid = 1;
    out->point_count = detail->surface.point_count;
    out->normal_count = detail->surface.normal_count;
    out->index_count = detail->surface.index_count;
    out->face_count = detail->surface.face_count;
    out->normals_per_index = detail->surface.normals_per_index;
    out->material_valid = detail->surface.material_valid;
    out->material_draw_mode = detail->surface.material_draw_mode;
    out->material_transparency = detail->surface.material_transparency;
    out->material_highlighted = detail->surface.material_highlighted;
    out->material_color[0] = detail->surface.material_color[0];
    out->material_color[1] = detail->surface.material_color[1];
    out->material_color[2] = detail->surface.material_color[2];
    out->cache_identity = rec->cache_identity;
    out->source_identity = rec->source_identity;

    return 1;
}


int
ged_draw_view_db_object_record_surface_index_at(
	const struct ged_draw_view_db_object_record *rec,
	size_t index,
	int *out)
{
    const struct ged_draw_view_export_detail *detail =
	_ged_draw_view_export_detail_from_record(rec);
    if (!out || !detail ||
	    detail->geometry_kind !=
	    GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_FACE_SET ||
	    index >= detail->surface.index_count ||
	    !detail->surface.indices)
	return 0;

    *out = detail->surface.indices[index];
    return 1;
}


static int
_ged_draw_view_rendered_object_summary_bsg(
	void *view_ctx,
	uint64_t cache_identity,
	struct ged_draw_view_rendered_object_summary *out)
{
    struct bsg_view *v = (struct bsg_view *)view_ctx;
    struct bsg_backend_scene *scene = NULL;
    int ret = 0;

    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!v || cache_identity == 0)
	return 0;

    scene = bsg_backend_scene_create();
    if (!scene)
	return 0;

    if (bsg_backend_scene_render_request(v, scene,
	    BSG_RENDER_FLAG_VISIBLE_ONLY) <= 0)
	goto cleanup;

    ret = 1;
    const struct bsg_backend_scene_node *node =
	bsg_backend_scene_find(scene, cache_identity);
    if (!node)
	goto cleanup;

    out->found = 1;
    out->source_identity = node->source_identity;
    out->material_draw_mode = node->material.draw_mode;
    out->material_transparency = node->material.transparency;
    out->selection_visible = node->selection.visible;
    out->selection_highlighted = node->selection.highlighted;

    if (node->geometry.kind == BSG_RENDER_GEOMETRY_INDEXED_FACE_SET) {
	out->is_indexed_face_set = 1;
	out->surface_point_count = node->geometry.surface.point_count;
	out->surface_normal_count = node->geometry.surface.normal_count;
	out->surface_index_count = node->geometry.surface.index_count;
	out->surface_face_count = node->geometry.surface.face_count;
	out->surface_normals_per_index =
	    (node->geometry.surface.normal_kind ==
	     BSG_RENDER_SURFACE_NORMALS_PER_INDEX);
	out->surface_material_valid =
	    node->geometry.surface.material_valid;
    } else if (node->geometry.kind == BSG_RENDER_GEOMETRY_ANNOTATION) {
	out->is_annotation = 1;
	out->annotation_segment_count =
	    node->geometry.annotation.segment_count;
    }

cleanup:
    bsg_backend_scene_destroy(scene);
    return ret;
}


int
ged_draw_view_rendered_object_summary(
	void *view_ctx,
	uint64_t cache_identity,
	struct ged_draw_view_rendered_object_summary *out)
{
    return _ged_draw_view_rendered_object_summary_bsg(view_ctx, cache_identity,
	    out);
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
