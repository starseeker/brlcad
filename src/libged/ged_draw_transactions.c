/*                    G E D _ D R A W _ T R A N S A C T I O N S . C
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
/** @file libged/ged_draw_transactions.c
 *
 * Draw transactions, database invalidation, and redraw synchronization.
 */

#include "common.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "bu/ptbl.h"
#include "bu/str.h"
#include "bu/color.h"
#include "bu/hash.h"
#include "bu/malloc.h"
#include "bg/clip.h"
#include "brlobol/draw_cache.h"
#include "bv.h"

#include "ged.h"
#include "ged/draw.h"
#include "ged/draw_obol.h"
#include "rt/calc.h"
#include "rt/view.h"
#include "ged/selection_state.h"
#include "../librt/librt_private.h"
#include "./ged_private.h"
#include "./ged_draw_private.h"
#include "./ged_draw_view_private.h"

struct ged_draw_db_update_ctx {
    struct ged *gedp;
    const char *path;
    ged_draw_stale_reason reason;
    int changed;
};


struct ged_draw_transparency_ctx {
    struct ged *gedp;
    struct directory **dpp;
    fastf_t transparency;
    int changed;
};


struct ged_draw_visibility_ctx {
    struct ged *gedp;
    const char *path;
    const char *basename;
    int visible;
    int changed;
};


struct ged_draw_highlight_ctx {
    struct ged *gedp;
    const char *basename;
    int highlighted;
    int matched;
};


struct ged_draw_redraw_path_ctx {
    struct ged *gedp;
    struct db_full_path *obj_path;
    void *view_ctx;
    int found;
    int ret;
};


struct ged_draw_redraw_shape_entry {
    ged_draw_shape_ref ref;
    void *view_ctx;
};


struct ged_draw_redraw_source_ctx {
    struct ged *gedp;
    const char *path;
    void *view_ctx;
    struct bu_ptbl shape_refs;
};


struct ged_draw_reexpand_group_entry {
    ged_draw_group_ref ref;
    char *path;
    void *view_ctx;
    struct ged_draw_appearance_settings appearance;
};


struct ged_draw_reexpand_source_ctx {
    struct ged *gedp;
    const char *path;
    void *view_ctx;
    struct bu_ptbl groups;
};


struct ged_draw_lod_finalize_ctx {
    struct ged *gedp;
    void *first_view_ctx;
    void **view_ctxs;
    size_t view_ctx_count;
    int ensured;
};


struct ged_draw_rename_ctx {
    struct ged *gedp;
    const char *old_path;
    const char *new_path;
    int changed;
};


struct ged_draw_rename_source_entry {
    char *old_path;
    char *new_path;
};


struct ged_draw_rename_source_ctx {
    const char *old_path;
    const char *new_path;
    struct bu_ptbl entries;
};


struct ged_draw_observer_entry {
    ged_draw_observer_token token;
    ged_draw_transaction_observer_func_t func;
    void *client_data;
    int active;
};


static struct ged_drawable *
_ged_draw_gdp(struct ged *gedp)
{
    return (gedp && gedp->i) ? gedp->i->ged_gdp : NULL;
}


static void *
_ged_draw_shared_fallback_view_ctx(struct ged *gedp)
{
    if (!gedp)
	return NULL;

    void *active = ged_draw_active_view_ctx(gedp);
    if (active && !ged_view_context_is_independent(active))
	return active;

    struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
    if (views) {
	for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	    void *view_ctx = BU_PTBL_GET(views, i);
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
			ged_draw_observer_token token)
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
	if (entry && entry->active && entry->func)
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
	if (entry && (!entry->active || !entry->func)) {
	    bu_ptbl_rm(&gdp->gd_draw_observers, (long *)entry);
	    _ged_draw_observer_entry_free(entry);
	    continue;
	}
	i++;
    }
}


static void
_ged_draw_observers_dispatch(struct ged *gedp,
			     const struct ged_draw_transaction *txn,
			     const struct ged_draw_transaction_result *result)
{
    struct ged_drawable *gdp = _ged_draw_gdp(gedp);
    if (!gdp || !gdp->gd_draw_observers_init || !txn || !result ||
	result->status <= 0)
	return;

    size_t len = BU_PTBL_LEN(&gdp->gd_draw_observers);
    gdp->gd_draw_observer_dispatch_depth++;
    for (size_t i = 0; i < len; i++) {
	struct ged_draw_observer_entry *entry =
	    (struct ged_draw_observer_entry *)BU_PTBL_GET(
		&gdp->gd_draw_observers, i);
	if (!entry || !entry->active || !entry->func)
	    continue;
	entry->func(gedp, txn, result, entry->client_data);
    }
    gdp->gd_draw_observer_dispatch_depth--;
    _ged_draw_observers_prune(gdp);
}


ged_draw_observer_token
ged_draw_observer_add(struct ged *gedp,
		      ged_draw_transaction_observer_func_t func,
		      void *client_data)
{
    struct ged_drawable *gdp = _ged_draw_gdp(gedp);
    if (!gdp || !func)
	return 0;

    _ged_draw_observers_init_if_needed(gdp);

    struct ged_draw_observer_entry *entry =
	(struct ged_draw_observer_entry *)bu_calloc(1, sizeof(*entry),
	    "ged draw observer entry");
    entry->token = gdp->gd_draw_next_observer_token++;
    if (!entry->token)
	entry->token = gdp->gd_draw_next_observer_token++;
    entry->func = func;
    entry->client_data = client_data;
    entry->active = 1;
    bu_ptbl_ins(&gdp->gd_draw_observers, (long *)entry);
    return entry->token;
}


int
ged_draw_observer_remove(struct ged *gedp, ged_draw_observer_token token)
{
    struct ged_drawable *gdp = _ged_draw_gdp(gedp);
    struct ged_draw_observer_entry *entry =
	_ged_draw_observer_find(gdp, token);
    if (!entry || !entry->active)
	return 0;

    entry->active = 0;
    entry->func = NULL;
    if (gdp->gd_draw_observer_dispatch_depth > 0)
	return 1;

    bu_ptbl_rm(&gdp->gd_draw_observers, (long *)entry);
    _ged_draw_observer_entry_free(entry);
    return 1;
}


void
ged_draw_observers_free(struct ged *gedp)
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
_ged_draw_path_replace_component(struct bu_vls *out,
				 const char *path,
				 const char *old_name,
				 const char *new_name)
{
    if (!out || !path || !old_name || !new_name)
	return 0;

    path = ged_draw_dbpath_skip_lead_slash(path);
    old_name = ged_draw_dbpath_skip_lead_slash(old_name);
    new_name = ged_draw_dbpath_skip_lead_slash(new_name);
    size_t old_len = strlen(old_name);
    if (!*path || !old_len || !*new_name)
	return 0;

    int changed = 0;
    int first = 1;
    bu_vls_trunc(out, 0);
    const char *p = path;
    while (*p) {
	const char *slash = strchr(p, '/');
	size_t len = slash ? (size_t)(slash - p) : strlen(p);
	if (!first)
	    bu_vls_putc(out, '/');
	if (len == old_len && bu_strncmp(p, old_name, old_len) == 0) {
	    bu_vls_strcat(out, new_name);
	    changed = 1;
	} else {
	    bu_vls_strncat(out, p, len);
	}
	first = 0;
	if (!slash)
	    break;
	p = slash + 1;
    }
    return changed;
}


static int
_ged_draw_rename_group_cb(const struct ged_draw_group_record *rec, void *userdata)
{
    struct ged_draw_rename_ctx *ctx =
	(struct ged_draw_rename_ctx *)userdata;
    if (!ctx || !ctx->gedp || !ctx->old_path || !ctx->new_path ||
	!rec || !rec->path)
	return 1;

    const char *path = ged_draw_dbpath_skip_lead_slash(rec->path);
    const char *old_path = ged_draw_dbpath_skip_lead_slash(ctx->old_path);
    const char *new_path = ged_draw_dbpath_skip_lead_slash(ctx->new_path);
    struct bu_vls updated = BU_VLS_INIT_ZERO;
    if (!_ged_draw_path_replace_component(&updated, path, old_path, new_path)) {
	bu_vls_free(&updated);
	return 1;
    }

    struct db_full_path dfp;
    db_full_path_init(&dfp);
    if (db_string_to_path(&dfp, ctx->gedp->dbip, bu_vls_cstr(&updated)) == 0) {
	if (ged_draw_group_ref_set_dbpath(ctx->gedp, rec->ref, &dfp))
	    ctx->changed++;
    }
    db_free_full_path(&dfp);
    bu_vls_free(&updated);
    return 1;
}


static void
_ged_draw_rename_source_entry_free(struct ged_draw_rename_source_entry *entry)
{
    if (!entry)
	return;
    if (entry->old_path)
	bu_free(entry->old_path, "draw rename old source path");
    if (entry->new_path)
	bu_free(entry->new_path, "draw rename new source path");
    bu_free(entry, "draw rename source entry");
}


static int
_ged_draw_collect_obol_source_rename_cb(struct ged *UNUSED(gedp),
					const char *source_path,
					void *userdata)
{
    struct ged_draw_rename_source_ctx *ctx =
	(struct ged_draw_rename_source_ctx *)userdata;
    if (!ctx || !source_path || !source_path[0] ||
	!ctx->old_path || !ctx->new_path)
	return 1;

    struct bu_vls updated = BU_VLS_INIT_ZERO;
    if (!_ged_draw_path_replace_component(&updated, source_path,
					  ctx->old_path, ctx->new_path)) {
	bu_vls_free(&updated);
	return 1;
    }

    const char *old_path = ged_draw_dbpath_skip_lead_slash(source_path);
    const char *new_path = bu_vls_cstr(&updated);
    if (!new_path || !new_path[0] || BU_STR_EQUAL(old_path, new_path)) {
	bu_vls_free(&updated);
	return 1;
    }

    struct ged_draw_rename_source_entry *entry =
	(struct ged_draw_rename_source_entry *)bu_calloc(1, sizeof(*entry),
	    "draw rename source entry");
    entry->old_path = bu_strdup(old_path);
    entry->new_path = bu_strdup(new_path);
    bu_ptbl_ins(&ctx->entries, (long *)entry);

    bu_vls_free(&updated);
    return 1;
}


static int
_ged_draw_apply_obol_component_source_renames(struct ged *gedp,
	const char *old_path,
	const char *new_path)
{
    if (!gedp || !old_path || !new_path)
	return 0;

    struct ged_draw_rename_source_ctx ctx;
    ctx.old_path = old_path;
    ctx.new_path = new_path;
    bu_ptbl_init(&ctx.entries, 8, "draw rename source entries");

    int status = ged_draw_obol_database_source_paths_foreach(gedp, 1,
		 _ged_draw_collect_obol_source_rename_cb, &ctx);
    int changed = 0;
    if (status >= 0) {
	unsigned long long revision = ged_draw_scene_revision(gedp) + 1;
	for (size_t i = 0; i < BU_PTBL_LEN(&ctx.entries); i++) {
	    struct ged_draw_rename_source_entry *entry =
		(struct ged_draw_rename_source_entry *)BU_PTBL_GET(&ctx.entries,
		    i);
	    if (entry && ged_draw_obol_database_source_rename_for_path(gedp,
		    entry->old_path, entry->new_path, revision))
		changed++;
	}
    }

    for (size_t i = 0; i < BU_PTBL_LEN(&ctx.entries); i++)
	_ged_draw_rename_source_entry_free(
	    (struct ged_draw_rename_source_entry *)BU_PTBL_GET(&ctx.entries,
		i));
    bu_ptbl_free(&ctx.entries);

    return changed;
}


static int
_ged_draw_apply_database_rename(struct ged *gedp,
				const char *old_path,
				const char *new_path)
{
    if (!gedp || !old_path || !new_path)
	return 0;

    struct ged_draw_rename_ctx ctx;
    ctx.gedp = gedp;
    ctx.old_path = old_path;
    ctx.new_path = new_path;
    ctx.changed = ged_draw_obol_database_source_rename_for_path(gedp,
		  old_path, new_path, ged_draw_scene_revision(gedp) + 1) ? 1 : 0;
    ctx.changed += _ged_draw_apply_obol_component_source_renames(gedp,
		   old_path, new_path);
    ged_draw_foreach_group_record(gedp, _ged_draw_rename_group_cb, &ctx);
    return ctx.changed;
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
ged_draw_mark_database_change(struct ged *gedp,
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
_ged_draw_set_transparency_cb(const struct ged_draw_shape_record *rec,
			      void *userdata)
{
    struct ged_draw_transparency_ctx *ctx =
	(struct ged_draw_transparency_ctx *)userdata;
    if (!ctx || !ctx->dpp)
	return 1;

    if (!rec || !rec->fullpath)
	return 1;

    struct directory **tmp_dpp = ctx->dpp;
    size_t i = 0;
    for (i = 0;
	 i < rec->fullpath->fp_len && *tmp_dpp != RT_DIR_NULL;
	 ++i, ++tmp_dpp) {
	if (rec->fullpath->fp_names[i] != *tmp_dpp)
	    break;
    }

    if (*tmp_dpp != RT_DIR_NULL)
	return 1;

    if (!ZERO(rec->transparency - ctx->transparency) &&
	ged_draw_shape_ref_set_transparency(ctx->gedp, rec->ref,
					    ctx->transparency)) {
	ctx->changed++;
    }

    return 1;
}


static int
_ged_draw_set_transparency(struct ged *gedp, const char *path,
			   fastf_t transparency)
{
    if (!gedp || !path)
	return 0;

    struct directory **dpp = _ged_build_dpp(gedp, path);
    if (!dpp)
	return 0;

    struct ged_draw_transparency_ctx ctx;
    ctx.gedp = gedp;
    ctx.dpp = dpp;
    ctx.transparency = transparency;
    ctx.changed = 0;

    ged_draw_foreach_shape_record(gedp, _ged_draw_set_transparency_cb, &ctx);
    bu_free((void *)dpp, "_ged_draw_set_transparency: directory pointers");
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
_ged_draw_visibility_group_cb(const struct ged_draw_group_record *rec,
			      void *userdata)
{
    struct ged_draw_visibility_ctx *ctx =
	(struct ged_draw_visibility_ctx *)userdata;
    if (!ctx || !rec)
	return 1;

    const char *tail = rec->path;
    if (tail) {
	const char *p = tail;
	while (*p) {
	    if (*p == '/')
		tail = p + 1;
	    p++;
	}
    }
    if (!((rec->path && BU_STR_EQUAL(rec->path, ctx->path)) ||
	  (tail && BU_STR_EQUAL(tail, ctx->basename))))
	return 1;

    if (rec->visible != ctx->visible &&
	ged_draw_group_ref_set_visible(ctx->gedp, rec->ref,
				       ctx->visible)) {
	ctx->changed++;
    }

    return 1;
}


static int
_ged_draw_visibility_shape_cb(const struct ged_draw_shape_record *rec,
			      void *userdata)
{
    struct ged_draw_visibility_ctx *ctx =
	(struct ged_draw_visibility_ctx *)userdata;
    if (!ctx || !rec || !rec->fullpath)
	return 1;

    int match = 0;
    for (size_t i = 0; i < rec->fullpath->fp_len; i++) {
	struct directory *dp = rec->fullpath->fp_names[i];
	if (dp && BU_STR_EQUAL(dp->d_namep, ctx->basename)) {
	    match = 1;
	    break;
	}
    }
    if (!match)
	return 1;

    if (rec->visible != ctx->visible &&
	ged_draw_shape_ref_set_visible(ctx->gedp, rec->ref,
				       ctx->visible)) {
	ctx->changed++;
    }

    return 1;
}


static int
_ged_draw_set_visibility(struct ged *gedp, const char *path, int visible)
{
    if (!gedp || !path)
	return 0;

    const char *basename = strrchr(path, '/');
    basename = basename ? basename + 1 : path;

    struct ged_draw_visibility_ctx ctx;
    ctx.gedp = gedp;
    ctx.path = ged_draw_dbpath_skip_lead_slash(path);
    ctx.basename = basename;
    ctx.visible = visible ? 1 : 0;
    ctx.changed = 0;

    ged_draw_foreach_group_record(gedp, _ged_draw_visibility_group_cb, &ctx);
    ged_draw_foreach_shape_record(gedp, _ged_draw_visibility_shape_cb, &ctx);
    return ctx.changed;
}


static int
_ged_draw_highlight_shape_cb(const struct ged_draw_shape_record *rec,
			     void *userdata)
{
    struct ged_draw_highlight_ctx *ctx =
	(struct ged_draw_highlight_ctx *)userdata;
    if (!ctx || !rec || !rec->fullpath)
	return 1;

    for (size_t i = 0; i < rec->fullpath->fp_len; i++) {
	struct directory *dp = rec->fullpath->fp_names[i];
	if (!dp || !BU_STR_EQUAL(dp->d_namep, ctx->basename))
	    continue;
	ctx->matched++;
	(void)ged_draw_shape_set_highlighted(ctx->gedp, rec->ref,
					     ctx->highlighted);
	break;
    }

    return 1;
}


static int
_ged_draw_set_highlight(struct ged *gedp, const char *path, int highlighted)
{
    if (!gedp || !path)
	return 0;

    const char *basename = strrchr(path, '/');
    basename = basename ? basename + 1 : path;

    struct ged_draw_highlight_ctx ctx;
    ctx.gedp = gedp;
    ctx.basename = basename;
    ctx.highlighted = highlighted ? 1 : 0;
    ctx.matched = 0;

    ged_draw_foreach_shape_record(gedp, _ged_draw_highlight_shape_cb, &ctx);
    ged_draw_highlighted_shape_ref_invalidate(gedp);
    return ctx.matched;
}


static int
_ged_draw_redraw_group_ref(struct ged *gedp, ged_draw_group_ref ref,
			   int skip_subtractions, void *view_ctx)
{
    if (!gedp || ged_draw_group_ref_is_null(ref))
	return -1;

    struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
    if (!wdbp)
	return -1;

    return ged_draw_group_ref_redraw_wireframe(gedp, ref, gedp->dbip,
	    wdbp->wdb_initial_tree_state.ts_tol,
	    wdbp->wdb_initial_tree_state.ts_ttol,
	    view_ctx ? view_ctx : ged_draw_active_view_ctx(gedp),
	    skip_subtractions);
}


int
ged_draw_redraw_group_ref(struct ged *gedp, ged_draw_group_ref ref,
			  int skip_subtractions)
{
    return _ged_draw_redraw_group_ref(gedp, ref, skip_subtractions, NULL);
}


static int
_ged_draw_redraw_all_cb(const struct ged_draw_group_record *rec, void *ud)
{
    struct ged_draw_redraw_path_ctx *ctx =
	(struct ged_draw_redraw_path_ctx *)ud;
    if (!rec)
	return 1;
    if (ctx->view_ctx && !ged_draw_group_record_in_view(rec, ctx->view_ctx))
	return 1;
    void *redraw_view_ctx = ctx->view_ctx ? ctx->view_ctx :
			    ((rec->in_view_scope && rec->view) ? rec->view :
			     ged_draw_active_view_ctx(ctx->gedp));
    int ret = _ged_draw_redraw_group_ref(ctx->gedp, rec->ref, 0,
					 redraw_view_ctx);
    if (ret < 0) {
	ctx->ret = -1;
	return 0;
    }
    ctx->found++;
    return 1;
}


static int
_ged_draw_redraw_path_cb(const struct ged_draw_group_record *rec, void *ud)
{
    struct ged_draw_redraw_path_ctx *ctx =
	(struct ged_draw_redraw_path_ctx *)ud;
    if (!rec)
	return 1;
    if (ctx->view_ctx && !ged_draw_group_record_in_view(rec, ctx->view_ctx))
	return 1;
    if (!rec->fullpath || rec->fullpath->fp_len <= 0) {
	ctx->ret = -1;
	return 0;
    }

    if (db_full_path_match_top(rec->fullpath, ctx->obj_path)) {
	ctx->found = 1;
	void *redraw_view_ctx = ctx->view_ctx ? ctx->view_ctx :
				((rec->in_view_scope && rec->view) ? rec->view :
				 ged_draw_active_view_ctx(ctx->gedp));
	ctx->ret = _ged_draw_redraw_group_ref(ctx->gedp, rec->ref, 0,
					      redraw_view_ctx);
	return 0;
    }

    return 1;
}


static int
_ged_draw_redraw(struct ged *gedp, const char *path, void *view_ctx)
{
    if (!gedp)
	return -1;

    if (!path) {
	struct ged_draw_redraw_path_ctx ctx;
	ctx.gedp = gedp;
	ctx.obj_path = NULL;
	ctx.view_ctx = view_ctx;
	ctx.found = 0;
	ctx.ret = 0;
	ged_draw_foreach_group_record(gedp, _ged_draw_redraw_all_cb, &ctx);
	if (ctx.ret < 0)
	    return -1;
	return ctx.found;
    }

    struct db_full_path obj_path;
    db_full_path_init(&obj_path);
    int ret = db_string_to_path(&obj_path, gedp->dbip,
				ged_draw_dbpath_skip_lead_slash(path));
    if (ret < 0)
	return -1;

    struct ged_draw_redraw_path_ctx ctx;
    ctx.gedp = gedp;
    ctx.obj_path = &obj_path;
    ctx.view_ctx = view_ctx;
    ctx.found = 0;
    ctx.ret = 0;
    ged_draw_foreach_group_record(gedp, _ged_draw_redraw_path_cb, &ctx);
    db_free_full_path(&obj_path);

    if (ctx.ret < 0)
	return -1;
    return ctx.found ? 1 : 0;
}


static int
_ged_draw_redraw_source_add_shape(struct ged_draw_redraw_source_ctx *ctx,
				  ged_draw_shape_ref ref)
{
    if (!ctx || ged_draw_shape_ref_is_null(ref))
	return 1;

    struct ged_draw_group_record grec;
    if (!ged_draw_group_record_get(ctx->gedp,
				   ged_draw_group_ref_of_shape(ctx->gedp, ref), &grec))
	return 1;
    if (grec.is_overlay)
	return 1;
    if (ctx->view_ctx && !ged_draw_group_record_in_view(&grec, ctx->view_ctx))
	return 1;

    struct ged_draw_redraw_shape_entry *entry =
	(struct ged_draw_redraw_shape_entry *)bu_calloc(1, sizeof(*entry),
	    "redraw source shape entry");
    entry->ref = ref;
    entry->view_ctx = ctx->view_ctx ? ctx->view_ctx :
		      ((grec.in_view_scope && grec.view) ? grec.view :
		       ged_draw_active_view_ctx(ctx->gedp));
    bu_ptbl_ins(&ctx->shape_refs, (long *)entry);
    return 1;
}


static int
_ged_draw_redraw_source_shape_cb(const struct ged_draw_shape_record *rec,
				 void *userdata)
{
    struct ged_draw_redraw_source_ctx *ctx =
	(struct ged_draw_redraw_source_ctx *)userdata;
    if (!ctx || !rec)
	return 1;
    if (ctx->path && !_ged_draw_shape_record_has_component(rec, ctx->path))
	return 1;

    return _ged_draw_redraw_source_add_shape(ctx, rec->ref);
}


static int
_ged_draw_redraw_source_index_cb(ged_draw_shape_ref ref, void *userdata)
{
    struct ged_draw_redraw_source_ctx *ctx =
	(struct ged_draw_redraw_source_ctx *)userdata;
    if (!ctx || ged_draw_shape_ref_is_null(ref))
	return 1;

    return _ged_draw_redraw_source_add_shape(ctx, ref);
}


static int
_ged_draw_redraw_source(struct ged *gedp, const char *path,
			void *view_ctx)
{
    if (!gedp || !gedp->dbip)
	return -1;

    struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
    if (!wdbp)
	return -1;

    struct ged_draw_redraw_source_ctx ctx;
    ctx.gedp = gedp;
    ctx.path = path ? ged_draw_dbpath_skip_lead_slash(path) : NULL;
    ctx.view_ctx = view_ctx;
    bu_ptbl_init(&ctx.shape_refs, 8, "redraw source shape refs");

    int indexed = ctx.path ?
		  ged_draw_shape_ref_index_for_component(gedp, ctx.path,
		      _ged_draw_redraw_source_index_cb, &ctx) : -1;
    if (indexed < 0) {
	struct ged_drawable *gdp = _ged_draw_gdp(gedp);
	if (gdp)
	    gdp->gd_draw_index_slow_path_shape_scans++;
	ged_draw_foreach_shape_record(gedp, _ged_draw_redraw_source_shape_cb,
				      &ctx);
    }

    int redrawn = 0;
    int failed = 0;
    for (size_t i = 0; i < BU_PTBL_LEN(&ctx.shape_refs); i++) {
	struct ged_draw_redraw_shape_entry *entry =
	    (struct ged_draw_redraw_shape_entry *)BU_PTBL_GET(&ctx.shape_refs, i);
	if (!entry)
	    continue;
	entry->ref.scene_revision = 0;
	int ret = ged_draw_shape_ref_redraw_wireframe(gedp, entry->ref,
		  gedp->dbip, wdbp->wdb_initial_tree_state.ts_tol,
		  wdbp->wdb_initial_tree_state.ts_ttol,
		  entry->view_ctx ? entry->view_ctx :
		  ged_draw_active_view_ctx(gedp), 0);
	if (ret < 0)
	    failed = 1;
	else
	    redrawn++;
	bu_free(entry, "redraw source shape entry");
    }
    bu_ptbl_free(&ctx.shape_refs);

    return failed ? -1 : redrawn;
}


static int
ged_draw_txn_view_array(struct ged *gedp,
			void *view_ctx,
			void ***view_ctxs_out)
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
	    void *local_view_ctx = BU_PTBL_GET(views, i);
	    if (local_view_ctx && !ged_view_context_is_independent(local_view_ctx))
		count++;
	}
	if (!count)
	    count = 1;
    }

    void **out = (void **)bu_calloc(count, sizeof(void *),
				    "draw transaction view context array");
    if (ged_view_context_is_independent(view_ctx)) {
	out[0] = view_ctx;
    } else {
	size_t idx = 0;
	struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
	for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	    void *local_view_ctx = BU_PTBL_GET(views, i);
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


static int
ged_draw_prepare_views_for_transaction(struct ged *gedp,
				       void *view_ctx)
{
    void **view_ctxs = NULL;
    size_t view_ctx_count = ged_draw_txn_view_array(gedp, view_ctx, &view_ctxs);
    if (!view_ctx_count || !view_ctxs)
	return 0;

    for (size_t i = 0; i < view_ctx_count; i++) {
	if (view_ctxs[i])
	    ged_draw_view_context_lod_bounds_callback_set(view_ctxs[i]);
    }

    bu_free(view_ctxs, "draw transaction view context array");
    return (int)view_ctx_count;
}


static int
ged_draw_autoview_for_transaction(struct ged *gedp,
				  void *view_ctx,
				  const char **draw_paths,
				  int draw_count,
				  int allow_database_fallback)
{
    (void)draw_paths;
    (void)draw_count;

    void **view_ctxs = NULL;
    size_t view_ctx_count = ged_draw_txn_view_array(gedp, view_ctx, &view_ctxs);
    if (!view_ctx_count || !view_ctxs)
	return 0;

    vect_t obol_min, obol_max;
    int obol_empty = 1;
    int have_obol_bounds =
	ged_draw_obol_scene_database_autoview_bounds(gedp, &obol_min,
	    &obol_max, &obol_empty) && !obol_empty;

    int adjusted = 0;
    for (size_t i = 0; i < view_ctx_count; i++) {
	if (!view_ctxs[i])
	    continue;
	if (have_obol_bounds &&
	    bv_autoview_bounds(bv_context_view((struct bv_context *)view_ctxs[i]),
			       BV_AUTOVIEW_SCALE_DEFAULT, obol_min, obol_max)) {
	    adjusted++;
	} else if (allow_database_fallback) {
	    /* The historical RT autoview fallback was BSG-only.  Obol/libbv
	     * callers require explicit bounds; leave the view unchanged when
	     * none are available. */
	}
    }

    bu_free(view_ctxs, "draw transaction view context array");
    return adjusted;
}


static int
ged_draw_finalize_lod_shape_cb(const struct ged_draw_shape_record *rec,
			       void *userdata)
{
    struct ged_draw_lod_finalize_ctx *ctx =
	(struct ged_draw_lod_finalize_ctx *)userdata;
    if (!ctx || !rec)
	return 1;
    if (ged_draw_shape_ref_lod_ensure(ctx->gedp, rec->ref,
				      ctx->first_view_ctx, ctx->view_ctxs, ctx->view_ctx_count))
	ctx->ensured++;
    return 1;
}


static int
ged_draw_finalize_lod_for_transaction(struct ged *gedp,
				      void *view_ctx)
{
    if (!gedp || !view_ctx)
	return 0;

    void **view_ctxs = NULL;
    size_t view_ctx_count = ged_draw_txn_view_array(gedp, view_ctx, &view_ctxs);
    if (!view_ctx_count || !view_ctxs)
	return 0;

    struct ged_draw_lod_finalize_ctx ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.gedp = gedp;
    ctx.first_view_ctx = view_ctx;
    ctx.view_ctxs = view_ctxs;
    ctx.view_ctx_count = view_ctx_count;

    ged_draw_foreach_shape_record(gedp, ged_draw_finalize_lod_shape_cb, &ctx);

    bu_free(view_ctxs, "draw transaction view context array");
    return ctx.ensured;
}


static int
_ged_draw_apply_draw(struct ged *gedp,
		     const struct ged_draw_transaction *txn,
		     const char *path,
		     struct ged_draw_transaction_result *result)
{
    if (!gedp || !txn)
	return -1;

    void *view_ctx = txn->view ? txn->view :
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
			      sizeof(const char *), "draw transaction paths");
    int draw_count = 0;
    for (int i = 0; i < path_count; i++) {
	if (paths[i] && paths[i][0])
	    draw_paths[draw_count++] = paths[i];
    }
    if (draw_count <= 0) {
	bu_free((void *)draw_paths, "draw transaction paths");
	return 0;
    }
    struct ged_draw_appearance_settings neutral_settings =
	    GED_DRAW_APPEARANCE_SETTINGS_INIT;
    if (txn->appearance)
	neutral_settings =
	    *(const struct ged_draw_appearance_settings *)txn->appearance;

    void *saved_view = ged_draw_active_view_ctx(gedp);
    ged_draw_active_view_ctx_set(gedp, view_ctx);
    (void)ged_draw_prepare_views_for_transaction(gedp, view_ctx);

    int ret = ged_draw_obol_apply_draw_paths(gedp, view_ctx, draw_paths,
		  draw_count, &neutral_settings, result);
    if (ret < 0) {
	bu_vls_printf(gedp->ged_result_str,
		      "Obol draw provider failed for draw mode %d\n",
		      neutral_settings.draw_mode);
    }
    ged_draw_active_view_ctx_set(gedp, saved_view);

    if (ret < 0) {
	bu_free((void *)draw_paths, "draw transaction paths");
	return -1;
    }

    if (txn->autoview)
	(void)ged_draw_autoview_for_transaction(gedp, view_ctx, draw_paths,
						draw_count, 0);
    (void)ged_selection_draw_sync(gedp, NULL);

    if (result) {
	if (result->affected_groups <= 0)
	    result->affected_groups = draw_count;
	if (result->affected_shapes <= 0)
	    result->affected_shapes = draw_count;
    }
    bu_free((void *)draw_paths, "draw transaction paths");
    return draw_count;
}


static int
_ged_draw_source_is_comb(struct ged *gedp, const char *path)
{
    if (!gedp || !gedp->dbip || !path)
	return 0;

    const char *name = _ged_draw_component_name(path);
    if (!name || !*name)
	return 0;

    struct directory *dp = db_lookup(gedp->dbip, name, LOOKUP_QUIET);
    return (dp && (dp->d_flags & RT_DIR_COMB)) ? 1 : 0;
}


static int
_ged_draw_reexpand_group_seen(const struct ged_draw_reexpand_source_ctx *ctx,
			      ged_draw_group_ref ref)
{
    if (!ctx || ged_draw_group_ref_is_null(ref))
	return 1;

    for (size_t i = 0; i < BU_PTBL_LEN(&ctx->groups); i++) {
	const struct ged_draw_reexpand_group_entry *entry =
	    (const struct ged_draw_reexpand_group_entry *)BU_PTBL_GET(
		&ctx->groups, i);
	if (entry && entry->ref.token == ref.token)
	    return 1;
    }

    return 0;
}


static int
_ged_draw_reexpand_group_add(struct ged_draw_reexpand_source_ctx *ctx,
			     const struct ged_draw_group_record *rec)
{
    if (!ctx || !rec || !rec->path || rec->is_overlay ||
	ged_draw_group_ref_is_null(rec->ref))
	return 0;
    if (ctx->view_ctx && !ged_draw_group_record_in_view(rec, ctx->view_ctx))
	return 0;
    if (_ged_draw_reexpand_group_seen(ctx, rec->ref))
	return 0;

    struct ged_draw_reexpand_group_entry *entry =
	(struct ged_draw_reexpand_group_entry *)bu_calloc(1, sizeof(*entry),
	    "reexpand group entry");
    entry->ref = rec->ref;
    entry->path = bu_strdup(rec->path);
    entry->view_ctx = ctx->view_ctx ? ctx->view_ctx :
		      ((rec->in_view_scope && rec->view) ? rec->view :
		       ged_draw_active_view_ctx(ctx->gedp));
    struct ged_draw_appearance_settings appearance =
	    GED_DRAW_APPEARANCE_SETTINGS_INIT;
    if (!ged_draw_group_ref_appearance_settings(ctx->gedp, rec->ref,
	    &appearance)) {
	appearance.draw_mode = rec->draw_mode;
	appearance.transparency = rec->transparency;
    }
    entry->appearance = appearance;
    bu_ptbl_ins(&ctx->groups, (long *)entry);
    return 1;
}


static int
_ged_draw_reexpand_group_add_ref(struct ged_draw_reexpand_source_ctx *ctx,
				 ged_draw_group_ref ref)
{
    if (!ctx || ged_draw_group_ref_is_null(ref))
	return 0;

    struct ged_draw_group_record rec;
    if (!ged_draw_group_record_get(ctx->gedp, ref, &rec))
	return 0;
    return _ged_draw_reexpand_group_add(ctx, &rec);
}


static int
_ged_draw_reexpand_source_group_cb(const struct ged_draw_group_record *rec,
				   void *userdata)
{
    struct ged_draw_reexpand_source_ctx *ctx =
	(struct ged_draw_reexpand_source_ctx *)userdata;
    if (!ctx || !rec || !rec->path)
	return 1;
    if (ctx->path && !_ged_draw_path_has_component(rec->path, ctx->path))
	return 1;
    (void)_ged_draw_reexpand_group_add(ctx, rec);
    return 1;
}


static int
_ged_draw_reexpand_source_shape_cb(const struct ged_draw_shape_record *rec,
				   void *userdata)
{
    struct ged_draw_reexpand_source_ctx *ctx =
	(struct ged_draw_reexpand_source_ctx *)userdata;
    if (!ctx || !rec)
	return 1;
    if (ctx->path && !_ged_draw_shape_record_has_component(rec, ctx->path))
	return 1;

    struct ged_draw_group_record grec;
    if (!ged_draw_group_record_get(ctx->gedp, rec->group, &grec))
	return 1;
    (void)_ged_draw_reexpand_group_add(ctx, &grec);
    return 1;
}


static int
_ged_draw_reexpand_source_group_index_cb(ged_draw_group_ref ref,
	void *userdata)
{
    struct ged_draw_reexpand_source_ctx *ctx =
	(struct ged_draw_reexpand_source_ctx *)userdata;
    if (!ctx || ged_draw_group_ref_is_null(ref))
	return 1;

    (void)_ged_draw_reexpand_group_add_ref(ctx, ref);
    return 1;
}


static int
_ged_draw_reexpand_source_shape_index_cb(ged_draw_shape_ref ref,
	void *userdata)
{
    struct ged_draw_reexpand_source_ctx *ctx =
	(struct ged_draw_reexpand_source_ctx *)userdata;
    if (!ctx || ged_draw_shape_ref_is_null(ref))
	return 1;

    (void)_ged_draw_reexpand_group_add_ref(ctx,
					   ged_draw_group_ref_of_shape(ctx->gedp, ref));
    return 1;
}


static int
_ged_draw_reexpand_source_groups(struct ged *gedp, const char *path,
				 void *view_ctx)
{
    if (!gedp)
	return -1;

    struct ged_draw_reexpand_source_ctx ctx;
    ctx.gedp = gedp;
    ctx.path = path ? ged_draw_dbpath_skip_lead_slash(path) : NULL;
    ctx.view_ctx = view_ctx;
    bu_ptbl_init(&ctx.groups, 8, "reexpand source groups");

    int groups_indexed = ctx.path ?
			 ged_draw_group_ref_index_for_component(gedp, ctx.path,
			     _ged_draw_reexpand_source_group_index_cb, &ctx) : -1;
    int shapes_indexed = ctx.path ?
			 ged_draw_shape_ref_index_for_component(gedp, ctx.path,
			     _ged_draw_reexpand_source_shape_index_cb, &ctx) : -1;

    if (groups_indexed < 0) {
	struct ged_drawable *gdp = _ged_draw_gdp(gedp);
	if (gdp)
	    gdp->gd_draw_index_slow_path_group_scans++;
	ged_draw_foreach_group_record(gedp,
				      _ged_draw_reexpand_source_group_cb, &ctx);
    }
    if (shapes_indexed < 0) {
	struct ged_drawable *gdp = _ged_draw_gdp(gedp);
	if (gdp)
	    gdp->gd_draw_index_slow_path_shape_scans++;
	ged_draw_foreach_shape_record(gedp,
				      _ged_draw_reexpand_source_shape_cb, &ctx);
    }

    int reexpanded = 0;
    int failed = 0;
    for (size_t i = 0; i < BU_PTBL_LEN(&ctx.groups); i++) {
	struct ged_draw_reexpand_group_entry *entry =
	    (struct ged_draw_reexpand_group_entry *)BU_PTBL_GET(&ctx.groups, i);
	if (!entry)
	    continue;

	struct ged_draw_appearance_settings settings = entry->appearance;

	struct ged_draw_transaction txn =
	    ged_draw_transaction_make(GED_DRAW_TXN_DRAW, entry->path);
	txn.view = entry->view_ctx;
	txn.appearance = &settings;
	txn.mode = settings.draw_mode;
	txn.autoview = 0;
	if (_ged_draw_apply_draw(gedp, &txn, entry->path, NULL) < 0)
	    failed = 1;
	else
	    reexpanded++;

	bu_free(entry->path, "reexpand group path");
	bu_free(entry, "reexpand group entry");
    }
    bu_ptbl_free(&ctx.groups);

    return failed ? -1 : reexpanded;
}


struct ged_draw_transaction
ged_draw_transaction_make(ged_draw_transaction_kind kind, const char *path)
{
    struct ged_draw_transaction txn;
    memset(&txn, 0, sizeof(txn));
    txn.kind = kind;
    txn.path = path;
    txn.mode = -1;
    txn.stale_reason = GED_DRAW_STALE_SOURCE_CHANGED;
    return txn;
}


struct ged_draw_transaction
ged_draw_transaction_make_value(ged_draw_transaction_kind kind,
				const char *path,
				double value)
{
    struct ged_draw_transaction txn = ged_draw_transaction_make(kind, path);
    txn.value = value;
    return txn;
}


void
ged_draw_transaction_result_init(struct ged_draw_transaction_result *result)
{
    if (!result)
	return;
    memset(result, 0, sizeof(*result));
    BU_VLS_INIT(&result->names);
    BU_VLS_INIT(&result->errors);
}


void
ged_draw_transaction_result_free(struct ged_draw_transaction_result *result)
{
    if (!result)
	return;
    if (BU_VLS_IS_INITIALIZED(&result->names))
	bu_vls_free(&result->names);
    if (BU_VLS_IS_INITIALIZED(&result->errors))
	bu_vls_free(&result->errors);
    memset(result, 0, sizeof(*result));
}


static void
_ged_draw_txn_result_prepare(struct ged_draw_transaction_result *result,
			     struct ged *gedp)
{
    if (!result)
	return;
    if (!BU_VLS_IS_INITIALIZED(&result->names)) {
	BU_VLS_INIT(&result->names);
    } else {
	bu_vls_trunc(&result->names, 0);
    }
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
    result->scene_revision_before = ged_draw_scene_revision(gedp);
    result->scene_revision_after = result->scene_revision_before;
}


static const char *
_ged_draw_txn_path(struct ged *gedp,
		   const struct ged_draw_transaction *txn,
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
    bu_free(s, "ged_draw transaction path");
    (void)gedp;
    return bu_vls_cstr(tmp);
}


static void
_ged_draw_txn_note_name(struct ged_draw_transaction_result *result,
			const char *name)
{
    if (!result || !name || !*name)
	return;
    if (bu_vls_strlen(&result->names)) {
	bu_vls_putc(&result->names, ' ');
    }
    bu_vls_printf(&result->names, "%s", name);
}


static int
_ged_draw_txn_kind_changes_scene(enum ged_draw_transaction_kind kind)
{
    switch (kind) {
	case GED_DRAW_TXN_DRAW:
	case GED_DRAW_TXN_ERASE:
	case GED_DRAW_TXN_ERASE_PREFIX:
	case GED_DRAW_TXN_TEARDOWN:
	case GED_DRAW_TXN_CLEAR:
	case GED_DRAW_TXN_CLEAR_SCOPE:
	case GED_DRAW_TXN_TRANSPARENCY:
	case GED_DRAW_TXN_VISIBILITY:
	case GED_DRAW_TXN_HIGHLIGHT:
	case GED_DRAW_TXN_REFRESH_MATERIAL_COLORS:
	case GED_DRAW_TXN_REDRAW:
	case GED_DRAW_TXN_STALE_SOURCE:
	case GED_DRAW_TXN_SOURCE_UPDATED:
	case GED_DRAW_TXN_SOURCE_RENAMED:
	case GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED:
	    return 1;
	case GED_DRAW_TXN_NONE:
	case GED_DRAW_TXN_DEFAULT_DRAW_MODE:
	case GED_DRAW_TXN_MATERIAL_CHANGED:
	default:
	    return 0;
    }
}


static void
_ged_draw_txn_bump_revision_if_needed(struct ged *gedp,
				      const struct ged_draw_transaction *txn,
				      int ret,
				      uint64_t rev0)
{
    if (!gedp || !txn || ret <= 0 || !_ged_draw_txn_kind_changes_scene(txn->kind))
	return;
    if (ged_draw_scene_revision(gedp) != rev0)
	return;
    if (gedp->i && gedp->i->ged_gdp)
	gedp->i->ged_gdp->gd_draw_rev++;
}


static int
_ged_draw_apply_transaction_inner(struct ged *gedp,
				  const struct ged_draw_transaction *txn,
				  struct ged_draw_transaction_result *result,
				  const char *path)
{
    if (!gedp || !txn)
	return 0;

    uint64_t rev0 = gedp->i->ged_gdp->gd_draw_rev;
    int had_groups0 = ged_draw_has_groups(gedp);
    int ret = 0;

    switch (txn->kind) {
	case GED_DRAW_TXN_ERASE:
	    if (!path)
		return 0;
	    if (txn->view || txn->mode >= 0)
		ret = ged_draw_erase_path_string_scoped(gedp, path,
							txn->view, txn->mode);
	    else
		ret = ged_draw_erase_path_string(gedp, path);
	    if (result) {
		result->affected_groups = ret;
		result->affected_shapes = ret;
	    }
	    break;
	case GED_DRAW_TXN_ERASE_PREFIX:
	    if (!path)
		return 0;
	    if (txn->view || txn->mode >= 0)
		ret = ged_draw_erase_path_prefix_string_scoped(gedp, path,
		      txn->view, txn->mode);
	    else
		ret = ged_draw_erase_path_prefix_string(gedp, path);
	    if (result) {
		result->affected_groups = ret;
		result->affected_shapes = ret;
	    }
	    break;
	case GED_DRAW_TXN_TEARDOWN:
	case GED_DRAW_TXN_CLEAR: {
	    ged_draw_clear(gedp);
	    ret = (had_groups0 || rev0 != 0) ? 1 : 0;
	    if (result)
		result->affected_groups = had_groups0 ? 1 : 0;
	    break;
	}
	case GED_DRAW_TXN_CLEAR_SCOPE:
	    ret = ged_draw_clear_view(gedp, txn->view);
	    if (result)
		result->affected_groups = ret;
	    break;
	case GED_DRAW_TXN_TRANSPARENCY:
	    if (!path)
		return 0;
	    ret = _ged_draw_set_transparency(gedp, path, (fastf_t)txn->value);
	    if (result)
		result->affected_shapes = ret;
	    break;
	case GED_DRAW_TXN_DEFAULT_DRAW_MODE:
	    ret = _ged_draw_set_default_mode(gedp, (ged_draw_mode)((int)txn->value));
	    break;
	case GED_DRAW_TXN_VISIBILITY:
	    if (!path)
		return 0;
	    ret = _ged_draw_set_visibility(gedp, path, !ZERO(txn->value));
	    if (result) {
		result->affected_groups = ret;
		result->affected_shapes = ret;
	    }
	    break;
	case GED_DRAW_TXN_HIGHLIGHT:
	    if (!path)
		return 0;
	    ret = _ged_draw_set_highlight(gedp, path, !ZERO(txn->value));
	    if (result)
		result->affected_shapes = ret;
	    break;
	case GED_DRAW_TXN_MATERIAL_CHANGED:
	    ged_draw_bump_material_revision(gedp);
	    ret = 1;
	    if (result)
		result->affected_shapes = ged_draw_shape_count(gedp);
	    break;
	case GED_DRAW_TXN_REFRESH_MATERIAL_COLORS:
	    ged_draw_refresh_material_colors(gedp);
	    ret = 1;
	    if (result)
		result->affected_shapes = ged_draw_shape_count(gedp);
	    break;
	case GED_DRAW_TXN_REDRAW:
	    ret = _ged_draw_redraw(gedp, path, txn->view);
	    if (ret >= 0)
		(void)ged_draw_finalize_lod_for_transaction(gedp,
			txn->view ? txn->view : ged_draw_active_view_ctx(gedp));
	    if (result)
		result->redrawn_count = (ret > 0) ? ret : 0;
	    break;
	case GED_DRAW_TXN_STALE_SOURCE:
	    ret = ged_draw_mark_database_change(gedp, path, txn->stale_reason);
	    if (result)
		result->stale_count = ret;
	    break;
	case GED_DRAW_TXN_SOURCE_UPDATED:
	    ret = ged_draw_apply_database_update(gedp, path, txn->removed,
						 txn->redraw);
	    if (result) {
		if (txn->removed)
		    result->affected_groups = (ret > 0) ? ret : 0;
		else if (txn->redraw)
		    result->redrawn_count = (ret > 0) ? ret : 0;
		else
		    result->stale_count = (ret > 0) ? ret : 0;
	    }
	    break;
	case GED_DRAW_TXN_SOURCE_RENAMED:
	    if (!path || !txn->new_path)
		return 0;
	    ret = _ged_draw_apply_database_rename(gedp, path, txn->new_path);
	    if (result)
		result->affected_groups = (ret > 0) ? ret : 0;
	    break;
	case GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED:
	    if (!path)
		return 0;
	    ret = ged_draw_erase_nonroot_component_string_scoped(gedp, path,
		  txn->view, txn->mode);
	    if (result) {
		result->affected_groups = ret;
		result->affected_shapes = ret;
	    }
	    break;
	case GED_DRAW_TXN_DRAW:
	    ret = _ged_draw_apply_draw(gedp, txn, path, result);
	    break;
	case GED_DRAW_TXN_NONE:
	default:
	    ret = 0;
	    break;
    }

    _ged_draw_txn_bump_revision_if_needed(gedp, txn, ret, rev0);

    if (result && ret) {
	if (path) {
	    _ged_draw_txn_note_name(result, path);
	} else if (txn->paths && txn->path_count > 0) {
	    for (int i = 0; i < txn->path_count; i++)
		_ged_draw_txn_note_name(result, txn->paths[i]);
	}
    }
    return ret;
}


int
ged_draw_apply_transaction(struct ged *gedp,
			   const struct ged_draw_transaction *txn,
			   struct ged_draw_transaction_result *result)
{
    if (!gedp || !txn) {
	if (result) {
	    _ged_draw_txn_result_prepare(result, gedp);
	    result->status = 0;
	}
	return 0;
    }

    (void)ged_draw_obol_scene_controller_ensure_owned(gedp, 1);
    if (ged_draw_obol_scene_controller_full_synced(gedp))
	(void)ged_draw_source_root_attach_view_contexts(gedp,
		ged_draw_active_view_ctx(gedp), ged_view_set_views_ctx(gedp));

    struct ged_draw_transaction_result local_result;
    int use_local_result = 0;
    if (!result && (_ged_draw_observers_have_active(gedp) ||
	    _ged_draw_txn_kind_changes_scene(txn->kind))) {
	ged_draw_transaction_result_init(&local_result);
	result = &local_result;
	use_local_result = 1;
    }

    struct bu_vls path_tmp = BU_VLS_INIT_ZERO;
    const char *path = _ged_draw_txn_path(gedp, txn, &path_tmp);
    if (result)
	_ged_draw_txn_result_prepare(result, gedp);

    int ret = _ged_draw_apply_transaction_inner(gedp, txn, result, path);
    if (result) {
	result->status = ret;
	result->scene_revision_after = ged_draw_scene_revision(gedp);
	if (ret < 0) {
	    result->error_count = 1;
	    bu_vls_printf(&result->errors, "draw transaction failed");
	    if (path)
		bu_vls_printf(&result->errors, ": %s", path);
	} else if (ret == 0 && _ged_draw_txn_kind_changes_scene(txn->kind)) {
	    int obol_ret =
		ged_draw_obol_scene_sync_attached_transaction(gedp, txn,
		    result);
	    if (obol_ret > 0) {
		ret = obol_ret;
		result->status = ret;
		_ged_draw_txn_bump_revision_if_needed(gedp, txn, ret,
						      result->scene_revision_before);
		result->scene_revision_after = ged_draw_scene_revision(gedp);
	    }
	}
	if (ret >= 0 && _ged_draw_txn_kind_changes_scene(txn->kind))
	    (void)ged_selection_draw_sync(gedp, NULL);
	_ged_draw_observers_dispatch(gedp, txn, result);
    }

    bu_vls_free(&path_tmp);
    if (use_local_result)
	ged_draw_transaction_result_free(&local_result);
    return ret;
}


int
ged_draw_apply_database_update(struct ged *gedp,
			       const char *path,
			       int removed,
			       int redraw)
{
    if (!gedp)
	return 0;
    if (removed) {
	if (!path)
	    return 0;
	return ged_draw_erase_component_string_scoped(gedp,
		ged_draw_dbpath_skip_lead_slash(path), NULL, -1);
    }

    int marked = ged_draw_mark_database_change(gedp, path,
		 GED_DRAW_STALE_SOURCE_CHANGED);
    if (!redraw)
	return marked;

    if (path && _ged_draw_source_is_comb(gedp, path)) {
	int reexpanded = _ged_draw_reexpand_source_groups(gedp, path, NULL);
	if (reexpanded < 0)
	    return reexpanded;
	if (reexpanded > 0)
	    return marked ? marked : reexpanded;
    }

    int redrawn = _ged_draw_redraw_source(gedp, path, NULL);
    if (redrawn < 0)
	return redrawn;
    return marked ? marked : redrawn;
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
