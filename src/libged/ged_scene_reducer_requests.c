/*             G E D _ S C E N E _ R E D U C E R _ R E Q U E S T S . C
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
/** @file libged/ged_scene_reducer_requests.c
 *
 * Semantic scene reduction, database invalidation, and redraw synchronization.
 */

#include "common.h"

#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include "bu/ptbl.h"
#include "bu/str.h"
#include "bu/datetime.h"
#include "bu/color.h"
#include "bu/hash.h"
#include "bu/malloc.h"
#include "bg/clip.h"
#include "bv.h"

#include "ged.h"
#include "ged/draw.h"
#include "ged/display.h"
#include "rt/calc.h"
#include "rt/view.h"
#include "ged/selection.h"
#include "./ged_private.h"
#include "./ged_draw_private.h"
#include "./ged_draw_view_private.h"
#include "./ged_scene_backend_private.h"

struct ged_draw_db_update_ctx {
    struct ged *gedp;
    const char *path;
    ged_draw_stale_reason reason;
    int changed;
};


struct ged_draw_observer_entry {
    ged_scene_observer_token token;
    ged_scene_observer_func_t scene_func;
    void *client_data;
    int active;
};


struct ged_scene_delta {
    enum ged_scene_delta_kind kind;
    struct ged_view_context *view_ctx;
    const char *const *paths;
    size_t path_count;
    int affects_all;
    size_t group_count;
    size_t shape_count;
    int presentation_only;
    uint64_t revision_before;
    uint64_t revision_after;
    const char *diagnostic;
};


static struct ged_drawable *
_ged_draw_gdp(struct ged *gedp)
{
    return (gedp && gedp->i) ? gedp->i->ged_gdp : NULL;
}


uint64_t
ged_scene_revision_advance(struct ged *gedp)
{
    struct ged_drawable *gdp = _ged_draw_gdp(gedp);
    if (!gdp)
	return 0;
    gdp->gd_draw_rev++;
    return gdp->gd_draw_rev;
}


static struct ged_view_context *
_ged_draw_shared_fallback_view_ctx(struct ged *gedp)
{
    if (!gedp)
	return NULL;

    struct ged_view_context *active = ged_draw_active_view_ctx(gedp);
    if (active && !ged_view_context_is_independent(active))
	return active;

    struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
    if (views) {
	for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	    struct ged_view_context *view_ctx =
		(struct ged_view_context *)BU_PTBL_GET(views, i);
	    if (view_ctx && !ged_view_context_is_independent(view_ctx))
		return view_ctx;
	}
    }

    return active;
}


static void
_ged_draw_observers_init_if_needed(struct ged_drawable *gdp)
{
    if (!gdp || gdp->gd_draw_observers_init)
	return;
    BU_PTBL_INIT(&gdp->gd_draw_observers);
    gdp->gd_draw_observers_init = 1;
    gdp->gd_draw_next_observer_token = 1;
    gdp->gd_draw_observer_dispatch_depth = 0;
}


static void
_ged_draw_observer_entry_free(struct ged_draw_observer_entry *entry)
{
    if (!entry)
	return;
    bu_free(entry, "ged draw observer entry");
}


static struct ged_draw_observer_entry *
_ged_draw_observer_find(struct ged_drawable *gdp,
			ged_scene_observer_token token)
{
    if (!gdp || !gdp->gd_draw_observers_init || !token)
	return NULL;
    for (size_t i = 0; i < BU_PTBL_LEN(&gdp->gd_draw_observers); i++) {
	struct ged_draw_observer_entry *entry =
	    (struct ged_draw_observer_entry *)BU_PTBL_GET(
		&gdp->gd_draw_observers, i);
	if (entry && entry->token == token)
	    return entry;
    }
    return NULL;
}


static int
_ged_draw_observers_have_active(struct ged *gedp)
{
    struct ged_drawable *gdp = _ged_draw_gdp(gedp);
    if (!gdp || !gdp->gd_draw_observers_init)
	return 0;
    for (size_t i = 0; i < BU_PTBL_LEN(&gdp->gd_draw_observers); i++) {
	struct ged_draw_observer_entry *entry =
	    (struct ged_draw_observer_entry *)BU_PTBL_GET(
		&gdp->gd_draw_observers, i);
	if (entry && entry->active && entry->scene_func)
	    return 1;
    }
    return 0;
}


static void
_ged_draw_observers_prune(struct ged_drawable *gdp)
{
    if (!gdp || !gdp->gd_draw_observers_init ||
	gdp->gd_draw_observer_dispatch_depth > 0)
	return;

    size_t i = 0;
    while (i < BU_PTBL_LEN(&gdp->gd_draw_observers)) {
	struct ged_draw_observer_entry *entry =
	    (struct ged_draw_observer_entry *)BU_PTBL_GET(
		&gdp->gd_draw_observers, i);
	if (entry && (!entry->active || !entry->scene_func)) {
	    bu_ptbl_rm(&gdp->gd_draw_observers, (long *)entry);
	    _ged_draw_observer_entry_free(entry);
	    continue;
	}
	i++;
    }
}


static void
_ged_scene_observers_dispatch(struct ged *gedp,
			      const struct ged_scene_delta *delta)
{
    struct ged_drawable *gdp = _ged_draw_gdp(gedp);
    if (!gdp || !gdp->gd_draw_observers_init || !delta)
	return;

    size_t len = BU_PTBL_LEN(&gdp->gd_draw_observers);
    gdp->gd_draw_observer_dispatch_depth++;
    for (size_t i = 0; i < len; i++) {
	struct ged_draw_observer_entry *entry =
	    (struct ged_draw_observer_entry *)BU_PTBL_GET(
		&gdp->gd_draw_observers, i);
	if (!entry || !entry->active || !entry->scene_func)
	    continue;
	entry->scene_func(gedp, delta, entry->client_data);
    }
    gdp->gd_draw_observer_dispatch_depth--;
    _ged_draw_observers_prune(gdp);
}


void
ged_scene_delta_dispatch_internal(
    struct ged *gedp,
    enum ged_scene_delta_kind kind,
    struct ged_view_context *view_ctx,
    const char *const *paths,
    size_t path_count,
    int affects_all,
    size_t group_count,
    size_t shape_count,
    int presentation_only,
    uint64_t revision_before,
    uint64_t revision_after,
    const char *diagnostic)
{
    struct ged_scene_delta delta;
    delta.kind = kind;
    delta.view_ctx = view_ctx;
    delta.paths = paths;
    delta.path_count = path_count;
    delta.affects_all = affects_all ? 1 : 0;
    delta.group_count = group_count;
    delta.shape_count = shape_count;
    delta.presentation_only = presentation_only ? 1 : 0;
    delta.revision_before = revision_before;
    delta.revision_after = revision_after;
    delta.diagnostic = diagnostic ? diagnostic : "";
    _ged_scene_observers_dispatch(gedp, &delta);
}


static enum ged_scene_delta_kind
_ged_scene_delta_kind_from_request(
    const struct ged_scene_reducer_request *transaction)
{
    if (!transaction)
	return GED_SCENE_DELTA_UNKNOWN;
    switch (transaction->kind) {
	case GED_SCENE_REDUCER_DRAW:
	    return GED_SCENE_DELTA_DRAW;
	case GED_SCENE_REDUCER_ERASE:
	case GED_SCENE_REDUCER_ERASE_PREFIX:
	    return GED_SCENE_DELTA_ERASE;
	case GED_SCENE_REDUCER_REDRAW:
	    return GED_SCENE_DELTA_REDRAW;
	case GED_SCENE_REDUCER_VISIBILITY:
	    return GED_SCENE_DELTA_VISIBILITY;
	case GED_SCENE_REDUCER_HIGHLIGHT:
	case GED_SCENE_REDUCER_HIGHLIGHTS_CLEAR:
	case GED_SCENE_REDUCER_HIGHLIGHT_OCCURRENCE:
	    return GED_SCENE_DELTA_HIGHLIGHT;
	case GED_SCENE_REDUCER_MATERIAL_CHANGED:
	case GED_SCENE_REDUCER_TRANSPARENCY:
	case GED_SCENE_REDUCER_DEFAULT_DRAW_MODE:
	    return GED_SCENE_DELTA_STYLE;
	case GED_SCENE_REDUCER_STALE_SOURCE:
	case GED_SCENE_REDUCER_SOURCE_UPDATED:
	case GED_SCENE_REDUCER_SOURCE_RENAMED:
	case GED_SCENE_REDUCER_SOURCE_REFERENCES_REMOVED:
	    return GED_SCENE_DELTA_SOURCE;
	case GED_SCENE_REDUCER_CLEAR:
	case GED_SCENE_REDUCER_CLEAR_SCOPE:
	    return GED_SCENE_DELTA_CLEAR;
	case GED_SCENE_REDUCER_TEARDOWN:
	    return GED_SCENE_DELTA_TEARDOWN;
	case GED_SCENE_REDUCER_NONE:
	default:
	    return GED_SCENE_DELTA_UNKNOWN;
    }
}


static void
_ged_draw_observers_dispatch(struct ged *gedp,
			     const struct ged_scene_reducer_request *txn,
			     const struct ged_scene_reducer_result *result,
			     const char *resolved_path)
{
    if (!txn || !result || result->status <= 0)
	return;

    const char *local_paths[2] = {NULL, NULL};
    const char *const *paths = NULL;
    size_t path_count = 0;
    if (txn->paths && txn->path_count > 0) {
	paths = (const char *const *)txn->paths;
	path_count = (size_t)txn->path_count;
    } else {
	if (resolved_path && resolved_path[0])
	    local_paths[path_count++] = resolved_path;
	if (txn->new_path && txn->new_path[0] &&
	    (!path_count || !BU_STR_EQUAL(local_paths[0], txn->new_path)))
	    local_paths[path_count++] = txn->new_path;
	paths = local_paths;
    }

    struct ged_scene_delta delta;
    delta.kind = _ged_scene_delta_kind_from_request(txn);
    delta.view_ctx = txn->view;
    delta.paths = paths;
    delta.path_count = path_count;
    delta.affects_all = path_count ? 0 : 1;
    delta.group_count = result->affected_groups > 0 ?
	(size_t)result->affected_groups : 0;
    delta.shape_count = result->affected_shapes > 0 ?
	(size_t)result->affected_shapes : 0;
    delta.presentation_only = result->presentation_only ? 1 : 0;
    delta.revision_before = result->scene_revision_before;
    delta.revision_after = result->scene_revision_after;
    delta.diagnostic = BU_VLS_IS_INITIALIZED(&result->errors) ?
	bu_vls_cstr(&result->errors) : "";
    _ged_scene_observers_dispatch(gedp, &delta);
}


ged_scene_observer_token
ged_scene_observer_add(struct ged *gedp,
		       ged_scene_observer_func_t callback,
		       void *client_data)
{
    struct ged_drawable *gdp = _ged_draw_gdp(gedp);
    if (!gdp || !callback)
	return 0;

    _ged_draw_observers_init_if_needed(gdp);

    struct ged_draw_observer_entry *entry =
	(struct ged_draw_observer_entry *)bu_calloc(1, sizeof(*entry),
	    "ged scene observer entry");
    entry->token = gdp->gd_draw_next_observer_token++;
    if (!entry->token)
	entry->token = gdp->gd_draw_next_observer_token++;
    entry->scene_func = callback;
    entry->client_data = client_data;
    entry->active = 1;
    bu_ptbl_ins(&gdp->gd_draw_observers, (long *)entry);
    return entry->token;
}


int
ged_scene_observer_remove(struct ged *gedp, ged_scene_observer_token token)
{
    struct ged_drawable *gdp = _ged_draw_gdp(gedp);
    struct ged_draw_observer_entry *entry =
	_ged_draw_observer_find(gdp, token);
    if (!entry || !entry->active)
	return 0;

    entry->active = 0;
    entry->scene_func = NULL;
    if (gdp->gd_draw_observer_dispatch_depth > 0)
	return 1;

    bu_ptbl_rm(&gdp->gd_draw_observers, (long *)entry);
    _ged_draw_observer_entry_free(entry);
    return 1;
}


enum ged_scene_delta_kind
ged_scene_delta_kind(const struct ged_scene_delta *delta)
{
    return delta ? delta->kind : GED_SCENE_DELTA_UNKNOWN;
}


struct ged_view_context *
ged_scene_delta_view(const struct ged_scene_delta *delta)
{
    return delta ? delta->view_ctx : NULL;
}


size_t
ged_scene_delta_path_count(const struct ged_scene_delta *delta)
{
    return delta ? delta->path_count : 0;
}


const char *
ged_scene_delta_path_at(const struct ged_scene_delta *delta, size_t index)
{
    return delta && delta->paths && index < delta->path_count ?
	delta->paths[index] : NULL;
}


int
ged_scene_delta_affects_all(const struct ged_scene_delta *delta)
{
    return delta && delta->affects_all ? 1 : 0;
}


size_t
ged_scene_delta_group_count(const struct ged_scene_delta *delta)
{
    return delta ? delta->group_count : 0;
}


size_t
ged_scene_delta_shape_count(const struct ged_scene_delta *delta)
{
    return delta ? delta->shape_count : 0;
}


int
ged_scene_delta_is_presentation_only(const struct ged_scene_delta *delta)
{
    return delta && delta->presentation_only ? 1 : 0;
}


uint64_t
ged_scene_delta_revision_before(const struct ged_scene_delta *delta)
{
    return delta ? delta->revision_before : 0;
}


uint64_t
ged_scene_delta_revision_after(const struct ged_scene_delta *delta)
{
    return delta ? delta->revision_after : 0;
}


const char *
ged_scene_delta_diagnostic(const struct ged_scene_delta *delta)
{
    return delta && delta->diagnostic ? delta->diagnostic : "";
}


void
ged_scene_observers_free(struct ged *gedp)
{
    struct ged_drawable *gdp = _ged_draw_gdp(gedp);
    if (!gdp || !gdp->gd_draw_observers_init)
	return;

    for (size_t i = 0; i < BU_PTBL_LEN(&gdp->gd_draw_observers); i++) {
	struct ged_draw_observer_entry *entry =
	    (struct ged_draw_observer_entry *)BU_PTBL_GET(
		&gdp->gd_draw_observers, i);
	_ged_draw_observer_entry_free(entry);
    }
    bu_ptbl_free(&gdp->gd_draw_observers);
    gdp->gd_draw_observers_init = 0;
    gdp->gd_draw_next_observer_token = 1;
    gdp->gd_draw_observer_dispatch_depth = 0;
}


static const char *
_ged_draw_component_name(const char *name)
{
    if (!name)
	return NULL;

    name = ged_draw_dbpath_skip_lead_slash(name);
    const char *basename = strrchr(name, '/');
    return basename ? basename + 1 : name;
}


static int
_dbfullpath_has_component(const struct db_full_path *fp, const char *name)
{
    if (!fp || !name)
	return 0;
    const char *basename = _ged_draw_component_name(name);
    if (!basename || !*basename)
	return 0;
    for (size_t i = 0; i < fp->fp_len; i++) {
	struct directory *dp = fp->fp_names[i];
	if (dp && BU_STR_EQUAL(dp->d_namep, basename))
	    return 1;
    }
    return 0;
}


static int
_ged_draw_path_has_component(const char *path, const char *name)
{
    if (!path || !name)
	return 0;

    path = ged_draw_dbpath_skip_lead_slash(path);
    const char *basename = _ged_draw_component_name(name);
    if (!*path || !basename || !*basename)
	return 0;

    size_t namelen = strlen(basename);
    const char *p = path;
    while (*p) {
	while (*p == '/')
	    p++;
	if (!*p)
	    break;
	const char *slash = strchr(p, '/');
	size_t len = slash ? (size_t)(slash - p) : strlen(p);
	if (len == namelen && bu_strncmp(p, basename, len) == 0)
	    return 1;
	if (!slash)
	    break;
	p = slash + 1;
    }

    return 0;
}


static int
_ged_draw_shape_record_has_component(const struct ged_draw_shape_record *rec,
				     const char *name)
{
    if (!rec || !name)
	return 0;

    if (rec->display_name &&
	_ged_draw_path_has_component(rec->display_name, name))
	return 1;
    if (rec->fullpath && _dbfullpath_has_component(rec->fullpath, name))
	return 1;

    const char *basename = _ged_draw_component_name(name);
    if (rec->leaf_name && basename && BU_STR_EQUAL(rec->leaf_name, basename))
	return 1;

    return 0;
}


static int
_ged_draw_mark_db_change_shape_ref(struct ged_draw_db_update_ctx *ctx,
				   ged_draw_shape_ref ref)
{
    if (!ctx || ged_draw_shape_ref_is_null(ref))
	return 1;

    if (ged_draw_shape_ref_mark_database_source_changed(ctx->gedp, ref,
	    ctx->reason))
	ctx->changed++;
    return 1;
}


static int
_ged_draw_mark_db_change_cb(const struct ged_draw_shape_record *rec,
			    void *userdata)
{
    struct ged_draw_db_update_ctx *ctx =
	(struct ged_draw_db_update_ctx *)userdata;
    if (!ctx || !rec)
	return 1;
    if (ctx->path && !_ged_draw_shape_record_has_component(rec, ctx->path))
	return 1;
    return _ged_draw_mark_db_change_shape_ref(ctx, rec->ref);
}


static int
_ged_draw_mark_db_change_index_cb(ged_draw_shape_ref ref, void *userdata)
{
    struct ged_draw_db_update_ctx *ctx =
	(struct ged_draw_db_update_ctx *)userdata;
    return _ged_draw_mark_db_change_shape_ref(ctx, ref);
}


int
ged_draw_source_mark_changed(struct ged *gedp,
			      const char *path,
			      ged_draw_stale_reason reason)
{
    if (!gedp)
	return 0;
    struct ged_draw_db_update_ctx ctx;
    ctx.gedp = gedp;
    ctx.path = path ? ged_draw_dbpath_skip_lead_slash(path) : NULL;
    ctx.reason = reason ? reason : GED_DRAW_STALE_SOURCE_CHANGED;
    ctx.changed = 0;
    int indexed = ctx.path ?
		  ged_draw_shape_ref_index_for_component(gedp, ctx.path,
		      _ged_draw_mark_db_change_index_cb, &ctx) : -1;
    if (indexed < 0) {
	struct ged_drawable *gdp = _ged_draw_gdp(gedp);
	if (gdp)
	    gdp->gd_draw_index_slow_path_shape_scans++;
	ged_draw_foreach_shape_record(gedp, _ged_draw_mark_db_change_cb,
				      &ctx);
    }
    return ctx.changed;
}


static int
_ged_draw_set_default_mode(struct ged *gedp, ged_draw_mode mode)
{
    if (!gedp)
	return 0;
    if (gedp->i->ged_gdp->gd_shaded_mode == (int)mode)
	return 0;
    gedp->i->ged_gdp->gd_shaded_mode = (int)mode;
    return 1;
}


static int
ged_draw_txn_view_array(struct ged *gedp,
			struct ged_view_context *view_ctx,
			struct ged_view_context ***view_ctxs_out)
{
    if (view_ctxs_out)
	*view_ctxs_out = NULL;
    if (!gedp || !view_ctx || !view_ctxs_out)
	return 0;

    size_t count = 0;
    if (ged_view_context_is_independent(view_ctx)) {
	count = 1;
    } else {
	struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
	for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	    struct ged_view_context *local_view_ctx =
		(struct ged_view_context *)BU_PTBL_GET(views, i);
	    if (local_view_ctx && !ged_view_context_is_independent(local_view_ctx))
		count++;
	}
	if (!count)
	    count = 1;
    }

    struct ged_view_context **out = (struct ged_view_context **)bu_calloc(
	    count, sizeof(struct ged_view_context *),
	    "scene reducer view context array");
    if (ged_view_context_is_independent(view_ctx)) {
	out[0] = view_ctx;
    } else {
	size_t idx = 0;
	struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
	for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	    struct ged_view_context *local_view_ctx =
		(struct ged_view_context *)BU_PTBL_GET(views, i);
	    if (local_view_ctx && !ged_view_context_is_independent(local_view_ctx))
		out[idx++] = local_view_ctx;
	}
	if (!idx)
	    out[idx++] = view_ctx;
	count = idx;
    }

    *view_ctxs_out = out;
    return (int)count;
}


int
ged_draw_autoview_for_transaction(struct ged *gedp,
				  struct ged_view_context *view_ctx,
				  const char **draw_paths,
				  int draw_count,
				  int allow_database_fallback)
{
    struct ged_view_context **view_ctxs = NULL;
    size_t view_ctx_count = ged_draw_txn_view_array(gedp, view_ctx, &view_ctxs);
    if (!view_ctx_count || !view_ctxs)
	return 0;

    vect_t database_min, database_max;
    int have_database_bounds = 0;
    if (allow_database_fallback && gedp->dbip &&
	draw_paths && draw_count > 0) {
	struct bu_vls messages = BU_VLS_INIT_ZERO;
	have_database_bounds = rt_obj_bounds(&messages, gedp->dbip,
	    draw_count, draw_paths, 0, database_min, database_max) ==
	    BRLCAD_OK;
	if (have_database_bounds &&
	    (!isfinite(database_min[X]) || !isfinite(database_min[Y]) ||
	     !isfinite(database_min[Z]) || !isfinite(database_max[X]) ||
	     !isfinite(database_max[Y]) || !isfinite(database_max[Z])))
	    have_database_bounds = 0;
	bu_vls_free(&messages);
    }

    int adjusted = 0;
    for (size_t i = 0; i < view_ctx_count; i++) {
	if (!view_ctxs[i])
	    continue;
	if (have_database_bounds &&
	    bv_autoview_bounds(bv_context_view((struct bv_context *)view_ctxs[i]),
		BV_AUTOVIEW_SCALE_DEFAULT, database_min, database_max)) {
	    adjusted++;
	}
    }

    bu_free(view_ctxs, "scene reducer view context array");
    return adjusted;
}


int
ged_draw_apply_draw_inner(struct ged *gedp,
			  const struct ged_scene_reducer_request *txn,
			  const char *path,
			  struct ged_scene_reducer_result *result)
{
    if (!gedp || !txn)
	return -1;

    struct ged_view_context *view_ctx = txn->view ? txn->view :
		     _ged_draw_shared_fallback_view_ctx(gedp);
    if (!view_ctx)
	return -1;

    const char *single_path[1] = {NULL};
    const char **paths = txn->paths;
    int path_count = txn->path_count;
    if ((!paths || path_count <= 0) && path) {
	single_path[0] = path;
	paths = single_path;
	path_count = 1;
    }
    if (!paths || path_count <= 0)
	return 0;

    const char **draw_paths = (const char **)bu_calloc(path_count + 1,
			      sizeof(const char *), "scene reducer paths");
    int draw_count = 0;
    for (int i = 0; i < path_count; i++) {
	if (paths[i] && paths[i][0])
	    draw_paths[draw_count++] = paths[i];
    }
    if (draw_count <= 0) {
	bu_free((void *)draw_paths, "scene reducer paths");
	return 0;
    }
    /* The semantic reducer records the requested roots here.  Heavy geometry
     * realization is performed once by the private backend adapter after the
     * semantic revision is committed and before application observers run. */
    if (result) {
	if (result->affected_groups <= 0)
	    result->affected_groups = draw_count;
	if (result->affected_shapes <= 0)
	    result->affected_shapes = draw_count;
    }
    bu_free((void *)draw_paths, "scene reducer paths");
    return draw_count;
}


struct ged_scene_reducer_request
ged_scene_reducer_request_make(ged_scene_reducer_operation kind, const char *path)
{
    struct ged_scene_reducer_request txn;
    memset(&txn, 0, sizeof(txn));
    txn.kind = kind;
    txn.path = path;
    txn.mode = -1;
    txn.stale_reason = GED_DRAW_STALE_SOURCE_CHANGED;
    return txn;
}


struct ged_scene_reducer_request
ged_scene_reducer_request_make_value(ged_scene_reducer_operation kind,
				const char *path,
				double value)
{
    struct ged_scene_reducer_request txn = ged_scene_reducer_request_make(kind, path);
    txn.value = value;
    return txn;
}


void
ged_scene_reducer_result_init(struct ged_scene_reducer_result *result)
{
    if (!result)
	return;
    memset(result, 0, sizeof(*result));
    BU_VLS_INIT(&result->errors);
}


static void
_ged_scene_reducer_paths_clear(struct ged_scene_reducer_result *result)
{
    size_t i;

    if (!result)
	return;
    for (i = 0; i < result->path_count; i++) {
	bu_free(result->paths[i], "scene reducer result path");
	result->paths[i] = NULL;
    }
    result->path_count = 0;
}


void
ged_scene_reducer_result_free(struct ged_scene_reducer_result *result)
{
    if (!result)
	return;
    _ged_scene_reducer_paths_clear(result);
    if (result->paths)
	bu_free(result->paths, "scene reducer result paths");
    if (BU_VLS_IS_INITIALIZED(&result->errors))
	bu_vls_free(&result->errors);
    memset(result, 0, sizeof(*result));
}


static void
_ged_scene_reducer_result_prepare(struct ged_scene_reducer_result *result,
			     struct ged *gedp)
{
    if (!result)
	return;
    _ged_scene_reducer_paths_clear(result);
    if (!BU_VLS_IS_INITIALIZED(&result->errors)) {
	BU_VLS_INIT(&result->errors);
    } else {
	bu_vls_trunc(&result->errors, 0);
    }
    result->status = 0;
    result->affected_groups = 0;
    result->affected_shapes = 0;
    result->stale_count = 0;
    result->redrawn_count = 0;
    result->error_count = 0;
    result->presentation_only = 0;
    result->progressive_data_complete = 0;
    result->scene_revision_before = ged_draw_scene_revision(gedp);
    result->scene_revision_after = result->scene_revision_before;
}


static const char *
_ged_scene_reducer_path(struct ged *gedp,
		   const struct ged_scene_reducer_request *txn,
		   struct bu_vls *tmp)
{
    if (!txn)
	return NULL;
    if (txn->path)
	return txn->path;
    if (!txn->dfp)
	return NULL;
    char *s = db_path_to_string(txn->dfp);
    if (!s)
	return NULL;
    bu_vls_sprintf(tmp, "%s", ged_draw_dbpath_skip_lead_slash(s));
    bu_free(s, "ged_scene reducer path");
    (void)gedp;
    return bu_vls_cstr(tmp);
}


static void
_ged_scene_reducer_note_name(struct ged_scene_reducer_result *result,
			const char *name)
{
    size_t i;

    if (!result || !name || !*name)
	return;
    for (i = 0; i < result->path_count; i++) {
	if (BU_STR_EQUAL(result->paths[i], name))
	    return;
    }
    if (result->path_count == result->path_capacity) {
	size_t new_capacity = result->path_capacity ?
	    result->path_capacity * 2 : 8;
	result->paths = (char **)bu_realloc(result->paths,
	    new_capacity * sizeof(char *),
	    "scene reducer result paths");
	memset(result->paths + result->path_capacity, 0,
	    (new_capacity - result->path_capacity) * sizeof(char *));
	result->path_capacity = new_capacity;
    }
    result->paths[result->path_count++] = bu_strdup(name);
}


static int
_ged_scene_reducer_kind_changes_scene(enum ged_scene_reducer_operation kind)
{
    switch (kind) {
	case GED_SCENE_REDUCER_DRAW:
	case GED_SCENE_REDUCER_ERASE:
	case GED_SCENE_REDUCER_ERASE_PREFIX:
	case GED_SCENE_REDUCER_TEARDOWN:
	case GED_SCENE_REDUCER_CLEAR:
	case GED_SCENE_REDUCER_CLEAR_SCOPE:
	case GED_SCENE_REDUCER_TRANSPARENCY:
	case GED_SCENE_REDUCER_VISIBILITY:
	case GED_SCENE_REDUCER_HIGHLIGHT:
	case GED_SCENE_REDUCER_HIGHLIGHTS_CLEAR:
	case GED_SCENE_REDUCER_HIGHLIGHT_OCCURRENCE:
	case GED_SCENE_REDUCER_MATERIAL_CHANGED:
	case GED_SCENE_REDUCER_REDRAW:
	case GED_SCENE_REDUCER_STALE_SOURCE:
	case GED_SCENE_REDUCER_SOURCE_UPDATED:
	case GED_SCENE_REDUCER_SOURCE_RENAMED:
	case GED_SCENE_REDUCER_SOURCE_REFERENCES_REMOVED:
	    return 1;
	case GED_SCENE_REDUCER_NONE:
	case GED_SCENE_REDUCER_DEFAULT_DRAW_MODE:
	default:
	    return 0;
    }
}


static void
_ged_scene_reducer_bump_revision_if_needed(struct ged *gedp,
				      const struct ged_scene_reducer_request *txn,
				      int ret,
				      uint64_t rev0)
{
    if (!gedp || !txn || ret <= 0 || !_ged_scene_reducer_kind_changes_scene(txn->kind))
	return;
    if (ged_draw_scene_revision(gedp) != rev0)
	return;
    if (gedp->i && gedp->i->ged_gdp)
	gedp->i->ged_gdp->gd_draw_rev++;
}


static int
_ged_scene_reduce_inner(struct ged *gedp,
				  const struct ged_scene_reducer_request *txn,
				  struct ged_scene_reducer_result *result,
				  const char *path)
{
    if (!gedp || !txn)
	return 0;

    uint64_t rev0 = gedp->i->ged_gdp->gd_draw_rev;
    int ret = 0;

    switch (txn->kind) {
	case GED_SCENE_REDUCER_ERASE:
	    if (!path)
		return 0;
	    ret = ged_draw_frontier_erase_path(gedp, path, txn->view,
		    txn->mode, 0, result);
	    if (result) {
		result->affected_groups = ret;
		result->affected_shapes = ret;
	    }
	    break;
	case GED_SCENE_REDUCER_ERASE_PREFIX:
	    if (!path)
		return 0;
	    ret = ged_draw_frontier_erase_path(gedp, path, txn->view,
		    txn->mode, 1, result);
	    if (result) {
		result->affected_groups = ret;
		result->affected_shapes = ret;
	    }
	    break;
	case GED_SCENE_REDUCER_TEARDOWN:
	    ret = ged_draw_frontier_clear(gedp, NULL, -1);
	    if (!ret)
		ret = 1;
	    ged_draw_highlighted_shape_ref_invalidate(gedp);
	    if (result)
		result->affected_groups = ret;
	    break;
	case GED_SCENE_REDUCER_CLEAR: {
	    ret = ged_draw_frontier_clear(gedp, NULL, -1);
	    ged_draw_highlighted_shape_ref_invalidate(gedp);
	    if (result)
		result->affected_groups = ret;
	    break;
	}
	case GED_SCENE_REDUCER_CLEAR_SCOPE:
	    ret = ged_draw_frontier_clear(gedp, txn->view, txn->mode);
	    if (result)
		result->affected_groups = ret;
	    break;
	case GED_SCENE_REDUCER_TRANSPARENCY:
	    if (!path)
		return 0;
	    ret = ged_draw_frontier_presentation_set(gedp, txn, path);
	    if (result)
		result->affected_shapes = ret;
	    break;
	case GED_SCENE_REDUCER_DEFAULT_DRAW_MODE:
	    ret = _ged_draw_set_default_mode(gedp, (ged_draw_mode)((int)txn->value));
	    break;
	case GED_SCENE_REDUCER_VISIBILITY:
	    if (!path)
		return 0;
	    ret = ged_draw_frontier_presentation_set(gedp, txn, path);
	    if (result) {
		result->affected_groups = ret;
		result->affected_shapes = ret;
	    }
	    break;
	case GED_SCENE_REDUCER_HIGHLIGHT:
	    if (!path)
		return 0;
	    ret = ged_draw_frontier_presentation_set(gedp, txn, path);
	    if (ret > 0)
		(void)ged_draw_highlight_revision_advance(gedp);
	    if (result)
		result->affected_shapes = ret;
	    break;
	case GED_SCENE_REDUCER_HIGHLIGHTS_CLEAR: {
	    const uint64_t highlight_rev = ged_draw_highlight_revision(gedp);
	    ged_draw_highlighted_shape_ref_invalidate(gedp);
	    ret = ged_draw_highlight_revision(gedp) != highlight_rev ? 1 : 0;
	    if (ged_draw_frontier_highlights_clear(gedp) > 0) {
		if (ged_draw_highlight_revision(gedp) == highlight_rev)
		    (void)ged_draw_highlight_revision_advance(gedp);
		ret = 1;
	    }
	    if (result)
		result->affected_shapes = ret;
	    break;
	}
	case GED_SCENE_REDUCER_HIGHLIGHT_OCCURRENCE: {
	    const uint64_t highlight_rev = ged_draw_highlight_revision(gedp);
	    if (ged_draw_shape_ref_is_null(txn->shape_ref)) {
		ged_draw_set_highlighted_shape_ref(gedp,
		    GED_DRAW_SHAPE_REF_NULL);
		ret = ged_draw_highlight_revision(gedp) != highlight_rev ? 1 : 0;
	    } else {
		ret = ged_draw_shape_set_highlighted(gedp, txn->shape_ref,
		    !ZERO(txn->value));
	    }
	    if (result)
		result->affected_shapes = ret > 0 ? 1 : 0;
	    break;
	}
	case GED_SCENE_REDUCER_MATERIAL_CHANGED:
	    ged_draw_bump_material_revision(gedp);
	    ret = 1;
	    if (result)
		result->affected_shapes = ged_draw_frontier_root_count(gedp);
	    break;
	case GED_SCENE_REDUCER_REDRAW:
	    ret = ged_draw_frontier_source_affected(gedp, path,
		(const char *const *)txn->paths,
		txn->path_count > 0 ? (size_t)txn->path_count : 0,
		txn->view, txn->mode);
	    if (result)
		result->redrawn_count = (ret > 0) ? ret : 0;
	    break;
	case GED_SCENE_REDUCER_STALE_SOURCE:
	    ret = ged_draw_frontier_source_affected(gedp, path,
		(const char *const *)txn->paths,
		txn->path_count > 0 ? (size_t)txn->path_count : 0,
		txn->view, txn->mode);
	    if (result)
		result->stale_count = ret;
	    break;
	case GED_SCENE_REDUCER_SOURCE_UPDATED:
	    ret = ged_draw_frontier_source_affected(gedp, path,
		(const char *const *)txn->paths,
		txn->path_count > 0 ? (size_t)txn->path_count : 0,
		txn->view, txn->mode);
	    if (result) {
		if (txn->removed)
		    result->affected_groups = (ret > 0) ? ret : 0;
		else if (txn->redraw)
		    result->redrawn_count = (ret > 0) ? ret : 0;
		else
		    result->stale_count = (ret > 0) ? ret : 0;
	    }
	    break;
	case GED_SCENE_REDUCER_SOURCE_RENAMED:
	    if (!path || !txn->new_path)
		return 0;
	    ret = ged_draw_frontier_source_affected(gedp, path, NULL, 0,
		txn->view, txn->mode);
	    {
		const int renamed = ged_draw_frontier_source_rename(gedp,
		    path, txn->new_path);
		if (renamed > ret)
		    ret = renamed;
	    }
	    if (result)
		result->affected_groups = (ret > 0) ? ret : 0;
	    break;
	case GED_SCENE_REDUCER_SOURCE_REFERENCES_REMOVED:
	    if (!path)
		return 0;
	    ret = ged_draw_frontier_source_affected(gedp, path,
		(const char *const *)txn->paths,
		txn->path_count > 0 ? (size_t)txn->path_count : 0,
		txn->view, txn->mode);
	    if (result) {
		result->affected_groups = ret;
		result->affected_shapes = ret;
	    }
	    break;
	case GED_SCENE_REDUCER_DRAW:
	    ret = ged_draw_frontier_absorb_draw(gedp, txn, path, result);
	    if (ret == 0)
		ret = ged_draw_apply_draw_inner(gedp, txn, path, result);
	    break;
	case GED_SCENE_REDUCER_NONE:
	default:
	    ret = 0;
	    break;
    }

    _ged_scene_reducer_bump_revision_if_needed(gedp, txn, ret, rev0);

    if (result && ret) {
	if (path) {
	    _ged_scene_reducer_note_name(result, path);
	} else if (txn->paths && txn->path_count > 0) {
	    for (int i = 0; i < txn->path_count; i++)
		_ged_scene_reducer_note_name(result, txn->paths[i]);
	}
    }
    return ret;
}


int
ged_scene_reduce(struct ged *gedp,
			   const struct ged_scene_reducer_request *txn,
			   struct ged_scene_reducer_result *result)
{
    if (!gedp || !txn) {
	if (result) {
	    _ged_scene_reducer_result_prepare(result, gedp);
	    result->status = 0;
	}
	return 0;
    }

    struct ged_scene_reducer_result local_result;
    int use_local_result = 0;
    if (!result && (_ged_draw_observers_have_active(gedp) ||
	    _ged_scene_reducer_kind_changes_scene(txn->kind))) {
	ged_scene_reducer_result_init(&local_result);
	result = &local_result;
	use_local_result = 1;
    }

    struct bu_vls path_tmp = BU_VLS_INIT_ZERO;
    const char *path = _ged_scene_reducer_path(gedp, txn, &path_tmp);
    if (result)
	_ged_scene_reducer_result_prepare(result, gedp);

    int reducer_ret = _ged_scene_reduce_inner(gedp, txn, result,
	path);
    int ret = reducer_ret;
    if (reducer_ret > 0)
	ged_draw_frontier_note_transaction(gedp, txn, path);
    if (result) {
	/* Give the active scene adapter a chance to handle backend-owned scene
	 * state even when the renderer-neutral frontier has no matching draw
	 * intent.  Custom plugin objects and retained Obol-only sources are valid
	 * scene citizens; requiring the semantic reducer to recognize them first
	 * made their update/erase transactions unreachable.
	 *
	 * A negative reducer result is a malformed/failed semantic transaction,
	 * not an extension point.  Otherwise the reducer and backend contribute
	 * independently to one committed result. */
	result->status = reducer_ret;
	result->scene_revision_after = ged_draw_scene_revision(gedp);
	if (reducer_ret < 0) {
	    result->error_count = 1;
	    bu_vls_printf(&result->errors, "scene reducer failed");
	    if (path)
		bu_vls_printf(&result->errors, ": %s", path);
	}
	if (reducer_ret >= 0 && _ged_scene_reducer_kind_changes_scene(txn->kind)) {
	    const int backend_ret =
		ged_scene_backend_apply_private(gedp, txn, result);
	    if (backend_ret > 0) {
		if (ret == 0)
		    ret = 1;
		/* Incremental frontier records are an optimization for the active
		 * adapter, not semantic state.  A late/replacement adapter receives
		 * the authoritative snapshot, so never replay stale records on its
		 * first later transaction. */
		ged_draw_frontier_visibility_changes_discard(gedp);
	    }
	}
	if (ret > 0) {
	    _ged_scene_reducer_bump_revision_if_needed(gedp, txn, ret,
		result->scene_revision_before);
	    result->status = ret;
	    result->scene_revision_after = ged_draw_scene_revision(gedp);
	}
	if (ret > 0 && _ged_scene_reducer_kind_changes_scene(txn->kind))
	    (void)ged_selection_present_private(gedp);
	if (ret > 0 && txn->kind == GED_SCENE_REDUCER_DRAW && txn->autoview) {
	    const struct ged_draw_appearance_settings default_appearance =
		GED_DRAW_APPEARANCE_SETTINGS_INIT;
	    const struct ged_draw_appearance_settings *appearance =
		txn->appearance ?
		(const struct ged_draw_appearance_settings *)txn->appearance :
		&default_appearance;
	    if (!appearance->defer_leaf_expansion) {
		const char *single_path[1] = {path};
		const char **draw_paths = txn->paths;
		int draw_count = txn->path_count;
		if ((!draw_paths || draw_count <= 0) && path) {
		    draw_paths = single_path;
		    draw_count = 1;
		}
		struct ged_view_context *autoview_ctx = txn->view ? txn->view :
		    _ged_draw_shared_fallback_view_ctx(gedp);
		(void)ged_draw_autoview_for_transaction(gedp, autoview_ctx,
		    draw_paths, draw_count, 1);
	    }
	}
	_ged_draw_observers_dispatch(gedp, txn, result, path);
    }

    bu_vls_free(&path_tmp);
    if (use_local_result)
	ged_scene_reducer_result_free(&local_result);
    return ret;
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
