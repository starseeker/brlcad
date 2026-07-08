/*                G E D _ D R A W _ R E C O R D S . C
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
/** @file ged_draw_records.c
 *
 * Semantic GED draw-record construction and iteration.
 *
 * This file is intentionally record/ref-facing.  It may resolve the temporary
 * draw registry to retained nodes internally, but it must not grow new raw
 * node-storage dependencies.
 */

#include "common.h"

#include <stdio.h>
#include <string.h>

#include "bu/hash.h"
#include "bu/malloc.h"
#include "bu/sort.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "ged/db_index.h"
#include "ged/draw.h"
#include "rt/view.h"
#include "./ged_private.h"
#include "./ged_draw_private.h"
#include "./ged_draw_view_private.h"


static const struct ged_draw_view_export_detail *
_ged_draw_view_export_detail_from_record(
	const struct ged_draw_view_db_object_record *rec)
{
    if (!rec || !rec->detail_token)
	return NULL;
    return (const struct ged_draw_view_export_detail *)rec->detail_token;
}


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

    char *path = NULL;
    const char *source_path = NULL;
    struct ged_draw_shape_record_summary shape_summary = {0};
    memset(&shape_summary, 0, sizeof(shape_summary));
    if (ged_draw_shape_ref_record_summary(gedp, ref, &shape_summary) &&
	    shape_summary.fullpath && shape_summary.fullpath->fp_len > 0) {
	path = db_path_to_string(shape_summary.fullpath);
	source_path = path;
    }
    if ((!source_path || !source_path[0]) && shape_summary.display_name)
	source_path = shape_summary.display_name;

    if (source_path && source_path[0])
	bu_vls_strcpy(path_out, _dbpath_skip_lead_slash(source_path));
    if (path)
	bu_free(path, "ged draw shape selection path");

    return bu_vls_strlen(path_out) > 0 ? 1 : 0;
}


static void
_draw_path_without_instance_suffixes(struct bu_vls *out, const char *path)
{
    if (!out)
	return;
    bu_vls_trunc(out, 0);
    if (!path)
	return;

    path = ged_draw_dbpath_skip_lead_slash(path);
    for (const char *cp = path; *cp; cp++) {
	if (*cp == '@' && cp[1] >= '0' && cp[1] <= '9') {
	    while (cp[1] && cp[1] != '/')
		cp++;
	    continue;
	}
	bu_vls_putc(out, *cp);
    }
}


static void
_draw_path_normalize(struct bu_vls *out, const char *path)
{
    _draw_path_without_instance_suffixes(out, path);
    if (!out)
	return;
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


static const char *
_draw_path_leaf_name(const char *path)
{
    if (!path)
	return NULL;
    path = ged_draw_dbpath_skip_lead_slash(path);
    const char *slash = strrchr(path, '/');
    return slash ? slash + 1 : path;
}


static int
_draw_path_database_path_exists(struct ged *gedp, const char *path)
{
    if (!gedp || !path || !path[0])
	return 0;
    if (!ged_db_index_available(gedp))
	return 0;
    struct bu_vls db_path = BU_VLS_INIT_ZERO;
    _draw_path_normalize(&db_path, path);
    int ret = ged_db_index_path_exists(gedp, bu_vls_cstr(&db_path));
    bu_vls_free(&db_path);
    return ret;
}


struct _draw_path_leaf_ctx {
    bu_hash_tbl *paths;
    size_t count;
};


static union tree *
_draw_path_make_nop_tree(void)
{
    union tree *curtree;

    BU_GET(curtree, union tree);
    RT_TREE_INIT(curtree);
    curtree->tr_op = OP_NOP;
    return curtree;
}


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
    return _draw_path_make_nop_tree();
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

    struct bu_vls db_path = BU_VLS_INIT_ZERO;
    _draw_path_normalize(&db_path, path);
    const char *av[1] = {bu_vls_cstr(&db_path)};
    int ret = db_walk_tree(gedp->dbip, 1, av, 1, &state,
	    NULL, NULL, _draw_path_leaf_cb, (void *)&ctx);
    db_free_db_tree_state(&state);
    bu_vls_free(&db_path);

    if (ret < 0)
	return 0;
    return ctx.count;
}


static int
_draw_group_record_summary_in_view(
	const struct ged_draw_group_record_summary *group_summary,
	void *view_ctx)
{
    if (!group_summary)
	return 0;

    if (!view_ctx)
	return !group_summary->in_view_scope;
    if (rt_view_context_is_independent(view_ctx))
	return group_summary->in_view_scope &&
	    group_summary->view_ctx == view_ctx;
    if (!group_summary->in_view_scope)
	return 1;
    return group_summary->view_ctx == view_ctx;
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
    bu_hash_tbl *drawn_leaf_names;
    int exact_shape;
    int descendant_shape;
    int ancestor_group;
    int matching_leaf;
};


static int
_draw_obol_record_instance_in_view(
	const struct ged_draw_obol_database_source_record *record,
	void *view_ctx)
{
    if (!record || !record->database_path || !record->database_path[0])
	return 0;

    const char *instance_key = record->instance_key;
    if (!view_ctx || !rt_view_context_is_independent(view_ctx)) {
	if (!instance_key || !instance_key[0])
	    return 1;

	const char mode_marker[] = ":ged-draw-mode:";
	struct bu_vls base_key = BU_VLS_INIT_ZERO;
	const char *marker = NULL;
	const char *candidate = instance_key;
	while ((candidate = strstr(candidate, mode_marker)) != NULL) {
	    marker = candidate;
	    candidate++;
	}
	if (marker)
	    bu_vls_strncpy(&base_key, instance_key,
		    (size_t)(marker - instance_key));
	else
	    bu_vls_strcpy(&base_key, instance_key);
	int ret = _draw_path_equal(bu_vls_cstr(&base_key),
		record->database_path);
	bu_vls_free(&base_key);
	return ret;
    }

    const char *view_name = rt_view_context_name_get(view_ctx);
    char fallback[64] = {0};
    if (!view_name || !view_name[0]) {
	snprintf(fallback, sizeof(fallback), "%p", view_ctx);
	view_name = fallback;
    }

    struct bu_vls prefix = BU_VLS_INIT_ZERO;
    bu_vls_printf(&prefix, "ged-view:%s:", view_name);
    int ret = (instance_key && bu_strncmp(instance_key, bu_vls_cstr(&prefix),
	    bu_vls_strlen(&prefix)) == 0) ? 1 : 0;
    bu_vls_free(&prefix);
    return ret;
}


static int
_draw_obol_source_record_visible_in_view(
	struct ged *gedp,
	const struct ged_draw_obol_database_source_record *record,
	void *view_ctx,
	int mode,
	struct ged_draw_group_record_summary *group_out)
{
    if (group_out)
	memset(group_out, 0, sizeof(*group_out));
    if (!gedp || !record || !record->valid || !record->database_path ||
	    !record->database_path[0] || !record->visible)
	return 0;
    if (mode >= 0 && record->draw_mode != mode)
	return 0;
    if (!_draw_obol_record_instance_in_view(record, view_ctx))
	return 0;

    if (record->owner_group_path && record->owner_group_path[0]) {
	struct ged_draw_group_record_summary group_summary = {0};
	if (ged_draw_obol_group_record_summary_for_path(gedp,
		record->owner_group_path, &group_summary)) {
	    if (group_summary.is_overlay || !group_summary.visible)
		return 0;
	    if (group_out)
		*group_out = group_summary;
	}
    }

    return 1;
}


static const char *
_draw_obol_source_record_path(struct bu_vls *storage,
			      const char *source_path,
			      const struct ged_draw_group_record_summary *UNUSED(group_summary))
{
    if (storage)
	bu_vls_trunc(storage, 0);
    if (!source_path || !source_path[0])
	return NULL;

    if (!storage)
	return ged_draw_dbpath_skip_lead_slash(source_path);

    _draw_path_normalize(storage, source_path);
    return bu_vls_cstr(storage);
}


static int
_draw_path_state_obol_source_record_cb(struct ged *gedp,
	const struct ged_draw_obol_database_source_record *record,
	void *ud)
{
    struct _draw_path_state_ctx *ctx =
	(struct _draw_path_state_ctx *)ud;
    if (!ctx || !gedp || !record || !record->database_path ||
	    !record->database_path[0])
	return 1;

    struct ged_draw_group_record_summary group_summary = {0};
    if (!_draw_obol_source_record_visible_in_view(gedp, record,
	    ctx->view_ctx, ctx->mode, &group_summary))
	return 1;
    if (group_summary.path && !_draw_path_equal(ctx->path, group_summary.path) &&
	    _draw_path_is_prefix(group_summary.path, ctx->path))
	ctx->ancestor_group = 1;

    struct bu_vls record_path_storage = BU_VLS_INIT_ZERO;
    const char *key = _draw_obol_source_record_path(&record_path_storage,
	    record->database_path, &group_summary);
    const char *query_leaf = _draw_path_leaf_name(ctx->path);
    const char *source_leaf = _draw_path_leaf_name(key);
    if (query_leaf && source_leaf && BU_STR_EQUAL(query_leaf, source_leaf))
	ctx->matching_leaf = 1;
    if (_draw_path_equal(ctx->path, key))
	ctx->exact_shape = 1;
    if (_draw_path_is_prefix(ctx->path, key))
	ctx->descendant_shape = 1;
    if (key && key[0]) {
	(void)bu_hash_set(ctx->drawn_leaf_paths,
		(const uint8_t *)key, strlen(key), (void *)1);
	if (source_leaf && source_leaf[0] && ctx->drawn_leaf_names)
	    (void)bu_hash_set(ctx->drawn_leaf_names,
		    (const uint8_t *)source_leaf, strlen(source_leaf),
		    (void *)1);
    }
    bu_vls_free(&record_path_storage);

    return 1;
}


static int
_draw_path_state_shape_ref_cb(ged_draw_shape_ref ref, void *ud)
{
    struct _draw_path_state_ctx *ctx =
	(struct _draw_path_state_ctx *)ud;
    if (!ctx || ged_draw_shape_ref_is_null(ref))
	return 1;

    struct ged_draw_shape_record_summary shape_summary = {0};
    if (!ged_draw_shape_ref_record_summary(ctx->gedp, ref, &shape_summary) ||
	    !shape_summary.fullpath || !shape_summary.visible)
	return 1;
    if (ctx->mode >= 0 && shape_summary.draw_mode != ctx->mode)
	return 1;

    struct ged_draw_group_record_summary group_summary;
    if (!ged_draw_shape_ref_owning_group_record_summary(
		ctx->gedp, ref, &group_summary) ||
	    group_summary.is_overlay || !group_summary.visible ||
	    !_draw_group_record_summary_in_view(&group_summary, ctx->view_ctx))
	return 1;

    char *path = db_path_to_string(shape_summary.fullpath);
    if (!path)
	return 1;
    const char *key = ged_draw_dbpath_skip_lead_slash(path);

    if (_draw_path_equal(ctx->path, key))
	ctx->exact_shape = 1;
    if (_draw_path_is_prefix(ctx->path, key))
	ctx->descendant_shape = 1;
    if (key && key[0]) {
	(void)bu_hash_set(ctx->drawn_leaf_paths,
		(const uint8_t *)key, strlen(key), (void *)1);
	const char *leaf = _draw_path_leaf_name(key);
	if (leaf && leaf[0] && ctx->drawn_leaf_names)
	    (void)bu_hash_set(ctx->drawn_leaf_names,
		    (const uint8_t *)leaf, strlen(leaf), (void *)1);
    }

    bu_free(path, "draw shape fullpath string");
    return 1;
}


static int
_draw_path_state_eval(struct ged *gedp,
		      const char *path,
		      struct _draw_path_state_ctx *ctx)
{
    if (!gedp || !path || !ctx || !ctx->drawn_leaf_paths)
	return 0;

    if (ctx->exact_shape)
	return 1;
    if (ctx->ancestor_group && ctx->matching_leaf &&
	    _draw_path_database_path_exists(gedp, path))
	return 1;

    if (!ctx->descendant_shape)
	return 0;

    bu_hash_tbl *expected = bu_hash_create(0);
    if (!expected)
	return 0;

    size_t expected_count =
	_draw_path_expected_leaf_paths(gedp, path, expected);

    int state = 0;
    if (expected_count > 0) {
	size_t matched = 0;
	bu_hash_entry *e = NULL;
	while ((e = bu_hash_next(expected, e)) != NULL) {
	    uint8_t *key = NULL;
	    size_t key_len = 0;
	    if (bu_hash_key(e, &key, &key_len) != 0)
		continue;
	    if (bu_hash_get(ctx->drawn_leaf_paths, key, key_len)) {
		matched++;
		continue;
	    }
	    struct bu_vls expected_path = BU_VLS_INIT_ZERO;
	    bu_vls_strncpy(&expected_path, (const char *)key, key_len);
	    const char *leaf = _draw_path_leaf_name(bu_vls_cstr(&expected_path));
	    if (leaf && leaf[0] && ctx->drawn_leaf_names &&
		    bu_hash_get(ctx->drawn_leaf_names,
			(const uint8_t *)leaf, strlen(leaf)))
		matched++;
	    bu_vls_free(&expected_path);
	}

	if (ctx->exact_shape || matched == expected_count)
	    state = 1;
	else if (matched > 0 || ctx->descendant_shape)
	    state = 2;
    } else {
	if (ctx->exact_shape)
	    state = 1;
	else if (ctx->descendant_shape)
	    state = 2;
    }

    bu_hash_destroy(expected);
    return state;
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
    bu_hash_tbl *drawn_names = bu_hash_create(0);
    if (!drawn || !drawn_names) {
	if (drawn)
	    bu_hash_destroy(drawn);
	if (drawn_names)
	    bu_hash_destroy(drawn_names);
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
    ctx.drawn_leaf_names = drawn_names;

    if (ged_draw_obol_scene_controller_full_synced(gedp)) {
	int obol_status = ged_draw_obol_database_source_records_foreach(gedp, 1,
		_draw_path_state_obol_source_record_cb, &ctx);
	if (obol_status >= 0) {
	    int state = _draw_path_state_eval(gedp, bu_vls_cstr(&norm),
		    &ctx);
	    bu_hash_destroy(drawn);
	    bu_hash_destroy(drawn_names);
	    bu_vls_free(&norm);
	    return state;
	}
    }

    ged_draw_source_root_foreach_shape_ref(gedp, 0,
	    _draw_path_state_shape_ref_cb, &ctx);

    int state = _draw_path_state_eval(gedp, bu_vls_cstr(&norm), &ctx);
    bu_hash_destroy(drawn);
    bu_hash_destroy(drawn_names);
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
_draw_list_shape_ref_path_collapsed_cb(ged_draw_shape_ref ref, void *ud)
{
    struct _draw_path_list_ctx *ctx = (struct _draw_path_list_ctx *)ud;
    if (!ctx || ged_draw_shape_ref_is_null(ref))
	return 1;

    struct ged_draw_shape_record_summary shape_summary;
    if (!ged_draw_shape_ref_record_summary(ctx->gedp, ref, &shape_summary) ||
	    !shape_summary.fullpath || !shape_summary.visible)
	return 1;
    if (ctx->mode >= 0 && shape_summary.draw_mode != ctx->mode)
	return 1;

    struct ged_draw_group_record_summary group_summary;
    if (!ged_draw_shape_ref_owning_group_record_summary(
		ctx->gedp, ref, &group_summary) ||
	    group_summary.is_overlay || !group_summary.visible ||
	    !_draw_group_record_summary_in_view(&group_summary, ctx->view_ctx))
	return 1;

    char *path = db_path_to_string(shape_summary.fullpath);
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
_draw_list_shape_ref_path_cb(ged_draw_shape_ref ref, void *ud)
{
    struct _draw_path_list_ctx *ctx = (struct _draw_path_list_ctx *)ud;
    if (!ctx || ged_draw_shape_ref_is_null(ref))
	return 1;

    struct ged_draw_shape_record_summary shape_summary;
    if (!ged_draw_shape_ref_record_summary(ctx->gedp, ref, &shape_summary) ||
	    !shape_summary.fullpath || !shape_summary.visible)
	return 1;
    if (ctx->mode >= 0 && shape_summary.draw_mode != ctx->mode)
	return 1;

    struct ged_draw_group_record_summary group_summary;
    if (!ged_draw_shape_ref_owning_group_record_summary(
		ctx->gedp, ref, &group_summary) ||
	    group_summary.is_overlay || !group_summary.visible ||
	    !_draw_group_record_summary_in_view(&group_summary, ctx->view_ctx))
	return 1;

    char *path = db_path_to_string(shape_summary.fullpath);
    if (!path)
	return 1;
    (void)_draw_path_list_add(ctx, path);
    bu_free(path, "draw list shape path");
    return 1;
}


static int
_draw_list_obol_source_record_cb(struct ged *gedp,
	const struct ged_draw_obol_database_source_record *record,
	void *ud)
{
    struct _draw_path_list_ctx *ctx = (struct _draw_path_list_ctx *)ud;
    if (!ctx || !record || !record->database_path ||
	    !record->database_path[0])
	return 1;

    struct ged_draw_group_record_summary group_summary = {0};
    if (!_draw_obol_source_record_visible_in_view(gedp, record,
	    ctx->view_ctx, ctx->mode, &group_summary))
	return 1;

    struct bu_vls record_path_storage = BU_VLS_INIT_ZERO;
    const char *record_path = _draw_obol_source_record_path(
	    &record_path_storage, record->database_path, &group_summary);
    (void)_draw_path_list_add(ctx, record_path);
    bu_vls_free(&record_path_storage);
    return 1;
}


static int
_draw_list_obol_source_record_collapsed_cb(struct ged *gedp,
	const struct ged_draw_obol_database_source_record *record,
	void *ud)
{
    struct _draw_path_list_ctx *ctx = (struct _draw_path_list_ctx *)ud;
    if (!ctx || !record || !record->database_path ||
	    !record->database_path[0])
	return 1;

    struct ged_draw_group_record_summary group_summary = {0};
    if (!_draw_obol_source_record_visible_in_view(gedp, record,
	    ctx->view_ctx, ctx->mode, &group_summary))
	return 1;

    struct bu_vls record_path_storage = BU_VLS_INIT_ZERO;
    struct bu_vls candidate = BU_VLS_INIT_ZERO;
    const char *record_path = _draw_obol_source_record_path(
	    &record_path_storage, record->database_path, &group_summary);
    const char *best = record_path;
    if (record_path && record_path[0]) {
	size_t len = strlen(record_path);
	for (size_t i = 0; i <= len; i++) {
	    if (record_path[i] != '/' && record_path[i] != '\0')
		continue;
	    bu_vls_trunc(&candidate, 0);
	    bu_vls_strncpy(&candidate, record_path, i);
	    if (_draw_path_list_cached_state(ctx, bu_vls_cstr(&candidate)) == 1) {
		best = bu_vls_cstr(&candidate);
		break;
	    }
	}
    }

    (void)_draw_path_list_add_collapsed(ctx, best);
    bu_vls_free(&candidate);
    bu_vls_free(&record_path_storage);
    return 1;
}


struct _draw_has_paths_ctx {
    struct ged *gedp;
    void *view_ctx;
    int mode;
    int found;
};


static int
_draw_has_paths_obol_source_record_cb(struct ged *gedp,
	const struct ged_draw_obol_database_source_record *record,
	void *ud)
{
    struct _draw_has_paths_ctx *ctx = (struct _draw_has_paths_ctx *)ud;
    if (!ctx || !gedp || !record || !record->database_path ||
	    !record->database_path[0])
	return 1;

    if (!_draw_obol_source_record_visible_in_view(gedp, record,
	    ctx->view_ctx, ctx->mode, NULL))
	return 1;

    ctx->found = 1;
    return 0;
}


static int
_draw_has_paths_shape_ref_cb(ged_draw_shape_ref ref, void *ud)
{
    struct _draw_has_paths_ctx *ctx = (struct _draw_has_paths_ctx *)ud;
    if (!ctx || ged_draw_shape_ref_is_null(ref))
	return 1;

    struct ged_draw_shape_record_summary shape_summary;
    if (!ged_draw_shape_ref_record_summary(ctx->gedp, ref, &shape_summary) ||
	    !shape_summary.fullpath || !shape_summary.visible)
	return 1;
    if (ctx->mode >= 0 && shape_summary.draw_mode != ctx->mode)
	return 1;

    struct ged_draw_group_record_summary group_summary;
    if (!ged_draw_shape_ref_owning_group_record_summary(
		ctx->gedp, ref, &group_summary) ||
	    group_summary.is_overlay || !group_summary.visible ||
	    !_draw_group_record_summary_in_view(&group_summary, ctx->view_ctx))
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

    int obol_status = -1;
    if (ged_draw_obol_scene_controller_full_synced(gedp)) {
	obol_status = ged_draw_obol_database_source_records_foreach(gedp, 1,
		expanded ? _draw_list_obol_source_record_cb :
		    _draw_list_obol_source_record_collapsed_cb, &ctx);
    }

    if (obol_status < 0) {
	if (expanded)
	    ged_draw_source_root_foreach_shape_ref(gedp, 0,
		    _draw_list_shape_ref_path_cb, &ctx);
	else
	    ged_draw_source_root_foreach_shape_ref(gedp, 0,
		    _draw_list_shape_ref_path_collapsed_cb, &ctx);
    }

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
    ctx.gedp = gedp;
    ctx.view_ctx = view_ctx;
    ctx.mode = mode;

    if (ged_draw_obol_scene_controller_full_synced(gedp)) {
	int obol_status = ged_draw_obol_database_source_records_foreach(gedp, 1,
		_draw_has_paths_obol_source_record_cb, &ctx);
	if (obol_status >= 0)
	    return ctx.found;
    }

    ged_draw_source_root_foreach_shape_ref(gedp, 0,
	    _draw_has_paths_shape_ref_cb, &ctx);
    return ctx.found;
}


static void
_ged_draw_fill_shape_record(struct ged *gedp,
			    ged_draw_shape_ref ref,
			    struct ged_draw_shape_record *out)
{
    memset(out, 0, sizeof(*out));
    out->ref = ref;
    out->group = ged_draw_group_ref_of_shape(gedp, out->ref);
    struct ged_draw_shape_record_summary shape_summary = {0};
    if (ged_draw_shape_ref_record_summary(gedp, ref, &shape_summary)) {
	out->fullpath = shape_summary.fullpath;
	out->display_name = shape_summary.display_name;
	out->leaf_name = shape_summary.leaf_name;
	out->path_hash = shape_summary.path_hash;
	out->source_revision = shape_summary.source_revision;
	out->inputs_revision = shape_summary.inputs_revision;
	out->realized_source_revision =
	    shape_summary.realized_source_revision;
	out->realized_inputs_revision =
	    shape_summary.realized_inputs_revision;
	out->stale = shape_summary.stale;
	out->stale_reason = shape_summary.stale_reason;
    }
    out->visible = shape_summary.visible;
    out->highlighted = shape_summary.highlighted;
    out->selected = 0;
    out->evaluated_region = shape_summary.evaluated_region;
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
    if (!out->stale_reason)
	out->stale_reason = ged_draw_stale_reason_name(GED_DRAW_STALE_NONE);

    struct ged_draw_database_source_summary source_summary;
    if (ged_draw_shape_ref_database_source_summary(gedp, ref,
	    &source_summary) &&
	    source_summary.has_state) {
	out->source_revision = source_summary.source_revision;
	out->inputs_revision = source_summary.inputs_revision;
	out->realized_source_revision =
	    source_summary.realized_source_revision;
	out->realized_inputs_revision =
	    source_summary.realized_inputs_revision;
	out->stale = source_summary.stale;
	out->stale_reason = source_summary.stale_reason;
    }
    out->drawn_revision = shape_summary.drawn_revision;
    out->transparency = shape_summary.transparency;
    out->draw_mode = shape_summary.draw_mode;
    out->line_width = shape_summary.line_width;
    struct ged_draw_scene_display_summary display_summary;
    if (ged_draw_shape_ref_display_summary(gedp, ref, &display_summary) &&
	    display_summary.valid) {
	out->visible = display_summary.visible;
	out->highlighted = display_summary.highlighted;
	out->transparency = display_summary.transparency;
	out->draw_mode = display_summary.draw_mode;
	out->line_width = display_summary.line_width;
    }
    VMOVE(out->center, shape_summary.center);
}


int
ged_draw_shape_record_get(struct ged *gedp,
			  ged_draw_shape_ref ref,
			  struct ged_draw_shape_record *out)
{
    if (!out)
	return 0;
    if (ged_draw_shape_ref_is_null(ref))
	return 0;
    struct ged_draw_shape_record_summary summary;
    if (!ged_draw_shape_ref_record_summary(gedp, ref, &summary))
	return 0;
    _ged_draw_fill_shape_record(gedp, ref, out);
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


int
ged_draw_group_record_get(struct ged *gedp,
			  ged_draw_group_ref ref,
			  struct ged_draw_group_record *out)
{
    if (!out)
	return 0;
    if (ged_draw_group_ref_is_null(ref))
	return 0;
    memset(out, 0, sizeof(*out));
    struct ged_draw_group_record_summary group_summary;
    if (!ged_draw_group_ref_record_summary(gedp, ref, &group_summary))
	return 0;
    out->ref = ref;
    out->path = group_summary.path;
    struct ged_draw_scene_tree_summary tree_summary;
    if (ged_draw_group_ref_tree_summary(gedp, ref, &tree_summary))
	out->fullpath = tree_summary.fullpath;
    out->view = group_summary.view_ctx;
    if (!ged_draw_group_ref_appearance_settings(gedp, ref, &out->appearance))
	out->appearance = (struct ged_draw_appearance_settings)GED_DRAW_APPEARANCE_SETTINGS_INIT;
    out->draw_mode = group_summary.draw_mode;
    out->transparency = group_summary.transparency;
    out->visible = group_summary.visible;
    out->is_overlay = group_summary.is_overlay;
    out->is_view_scope = group_summary.is_view_scope;
    out->in_view_scope = group_summary.in_view_scope;
    out->is_view_source = group_summary.is_view_source;
    out->is_local_source = group_summary.is_local_source;
    (void)ged_draw_group_ref_shape_count(gedp, ref, &out->shape_count);
    return 1;
}


struct _foreach_group_record_ctx {
    struct ged *gedp;
    int (*cb)(const struct ged_draw_group_record *, void *);
    void *userdata;
};


static int
_foreach_group_record_cb(ged_draw_group_ref group_ref, void *ud)
{
    struct _foreach_group_record_ctx *ctx =
	(struct _foreach_group_record_ctx *)ud;
    struct ged_draw_group_record rec;
    if (!ged_draw_group_record_get(ctx->gedp, group_ref, &rec))
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
    (void)ged_draw_source_root_foreach_group_ref(gedp,
	    _foreach_group_record_cb,
	    &ctx);
}


struct _foreach_shape_record_ctx {
    struct ged *gedp;
    int (*cb)(const struct ged_draw_shape_record *, void *);
    void *userdata;
};


static int
_foreach_shape_record_cb(ged_draw_shape_ref ref, void *ud)
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
    (void)ged_draw_source_root_foreach_shape_ref(gedp, 0,
	    _foreach_shape_record_cb, &ctx);
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
	query_flags |= GED_DRAW_VIEW_EXPORT_QUERY_VIEW_OBJECTS;
    if (query->flags & GED_DRAW_VIEW_RECORD_QUERY_DB_OBJECTS)
	query_flags |= GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS;
    if (query->flags & GED_DRAW_VIEW_RECORD_QUERY_LOCAL_ONLY)
	query_flags |= GED_DRAW_VIEW_EXPORT_QUERY_LOCAL_ONLY;

    ged_draw_view_context_foreach_export_record(view_ctx, query_flags,
	    GED_DRAW_VIEW_EXPORT_RENDER_PAYLOAD_PREPARE, query->glob,
	    query->draw_mode, cb, userdata);
}


void
ged_draw_foreach_view_db_object_record(void *view_ctx,
				       ged_draw_view_db_object_record_cb cb,
				       void *userdata)
{
    ged_draw_view_context_foreach_export_record(view_ctx,
	    GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS,
	    GED_DRAW_VIEW_EXPORT_RENDER_PAYLOAD_PREPARE, NULL,
	    GED_DRAW_VIEW_EXPORT_DRAW_MODE_ANY, cb, userdata);
}


void
ged_draw_foreach_visible_view_db_object_record(void *view_ctx,
					      ged_draw_view_db_object_record_cb cb,
					      void *userdata)
{
    ged_draw_view_context_foreach_export_record(view_ctx,
	    GED_DRAW_VIEW_EXPORT_QUERY_VISIBLE_ONLY |
	    GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS,
	    GED_DRAW_VIEW_EXPORT_RENDER_VISIBLE_ONLY |
	    GED_DRAW_VIEW_EXPORT_RENDER_PAYLOAD_PREPARE,
	    NULL, GED_DRAW_VIEW_EXPORT_DRAW_MODE_ANY, cb, userdata);
}


void
ged_draw_foreach_visible_view_db_object_record_mode(
	void *view_ctx,
	int draw_mode,
	ged_draw_view_db_object_record_cb cb,
	void *userdata)
{
    ged_draw_view_context_foreach_export_record(view_ctx,
	    GED_DRAW_VIEW_EXPORT_QUERY_VISIBLE_ONLY |
	    GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS,
	    GED_DRAW_VIEW_EXPORT_RENDER_VISIBLE_ONLY |
	    GED_DRAW_VIEW_EXPORT_RENDER_PAYLOAD_PREPARE,
	    NULL, draw_mode, cb, userdata);
}


void
ged_draw_foreach_visible_view_record(void *view_ctx,
				     ged_draw_view_db_object_record_cb cb,
				     void *userdata)
{
    ged_draw_view_context_foreach_export_record(view_ctx,
	    GED_DRAW_VIEW_EXPORT_QUERY_VISIBLE_ONLY,
	    GED_DRAW_VIEW_EXPORT_RENDER_VISIBLE_ONLY |
	    GED_DRAW_VIEW_EXPORT_RENDER_PAYLOAD_PREPARE,
	    NULL, GED_DRAW_VIEW_EXPORT_DRAW_MODE_ANY, cb, userdata);
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


struct _ged_draw_rendered_object_summary_ctx {
    uint64_t cache_identity;
    struct ged_draw_view_rendered_object_summary *out;
    int visited;
};


static int
_ged_draw_rendered_object_summary_record_cb(
	const struct ged_draw_view_db_object_record *rec,
	void *userdata)
{
    struct _ged_draw_rendered_object_summary_ctx *ctx =
	(struct _ged_draw_rendered_object_summary_ctx *)userdata;
    if (!ctx || !rec)
	return 1;

    ctx->visited = 1;
    if (rec->cache_identity != ctx->cache_identity)
	return 1;

    struct ged_draw_view_rendered_object_summary *out = ctx->out;
    out->found = 1;
    out->source_identity = rec->source_identity;
    out->material_draw_mode = rec->draw_mode;
    out->material_transparency = rec->transparency;
    out->selection_visible = rec->visible;
    out->selection_highlighted = rec->highlighted;

    struct ged_draw_view_surface_summary surface;
    if (ged_draw_view_db_object_record_surface_summary(rec, &surface)) {
	out->is_indexed_face_set = 1;
	out->surface_point_count = surface.point_count;
	out->surface_normal_count = surface.normal_count;
	out->surface_index_count = surface.index_count;
	out->surface_face_count = surface.face_count;
	out->surface_normals_per_index = surface.normals_per_index;
	out->surface_material_valid = surface.material_valid;
	out->material_draw_mode = surface.material_draw_mode;
	out->material_transparency = surface.material_transparency;
    } else {
	struct ged_draw_view_annotation_summary annotation;
	if (ged_draw_view_db_object_record_annotation_summary(rec, 0,
		&annotation)) {
	    out->is_annotation = 1;
	    out->annotation_segment_count = annotation.segment_count;
	}
    }

    return 0;
}


int
ged_draw_view_rendered_object_summary(
	void *view_ctx,
	uint64_t cache_identity,
	struct ged_draw_view_rendered_object_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!view_ctx || cache_identity == 0)
	return 0;

    struct _ged_draw_rendered_object_summary_ctx ctx;
    ctx.cache_identity = cache_identity;
    ctx.out = out;
    ctx.visited = 0;
    ged_draw_foreach_visible_view_record(view_ctx,
	    _ged_draw_rendered_object_summary_record_cb, &ctx);

    return ctx.visited;
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
