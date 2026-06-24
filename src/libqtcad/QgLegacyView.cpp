/*                    Q G L E G A C Y V I E W . C P P
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
/** @file QgLegacyView.cpp
 *
 * Neutral qtcad helpers for the staged legacy view handle.
 */

#include "common.h"

#include <new>

#include "bu/malloc.h"
#include "bu/vls.h"
#include "dm.h"
#include "ged.h"
#include "ged/bsg_ged_draw.h"
#include "rt/edit_legacy_bsg.h"
#include "QgLegacyViewContext.h"
#include "QgLegacyViewDm.h"
#include "QgLegacyViewSketch.h"
#include "qtcad/QgLegacyView.h"

extern "C" {
#include "rt/view_legacy_bsg.h"
}

static_assert(QG_LEGACY_VIEW_DRAW_MODE_WIRE == GED_DRAW_MODE_WIRE,
	"qtcad draw mode constants must track GED draw modes");
static_assert(QG_LEGACY_VIEW_DRAW_MODE_SHADED_BOTS ==
	GED_DRAW_MODE_SHADED_BOTS,
	"qtcad draw mode constants must track GED draw modes");

static struct dm *
qg_legacy_dm_to_dm(qg_legacy_dm *dmp)
{
    return reinterpret_cast<struct dm *>(dmp);
}

static qg_legacy_dm *
qg_legacy_dm_from_dm(struct dm *dmp)
{
    return reinterpret_cast<qg_legacy_dm *>(dmp);
}

static struct fb *
qg_legacy_fb_to_fb(qg_legacy_fb *ifp)
{
    return reinterpret_cast<struct fb *>(ifp);
}

static qg_legacy_fb *
qg_legacy_fb_from_fb(struct fb *ifp)
{
    return reinterpret_cast<qg_legacy_fb *>(ifp);
}
static_assert(QG_LEGACY_VIEW_DRAW_MODE_SHADED == GED_DRAW_MODE_SHADED,
	"qtcad draw mode constants must track GED draw modes");
static_assert(QG_LEGACY_VIEW_DRAW_MODE_EVAL_WIRE == GED_DRAW_MODE_EVAL_WIRE,
	"qtcad draw mode constants must track GED draw modes");
static_assert(QG_LEGACY_VIEW_DRAW_MODE_HIDDEN_LINE ==
	GED_DRAW_MODE_HIDDEN_LINE,
	"qtcad draw mode constants must track GED draw modes");
static_assert(QG_LEGACY_VIEW_DRAW_MODE_EVAL_POINTS ==
	GED_DRAW_MODE_EVAL_POINTS,
	"qtcad draw mode constants must track GED draw modes");
static_assert((int)QG_LEGACY_VIEW_FEATURE_UNKNOWN ==
	(int)RT_VIEW_FEATURE_UNKNOWN,
	"qtcad feature constants must track neutral RT feature constants");
static_assert((int)QG_LEGACY_VIEW_FEATURE_TRANSIENT_PREVIEW ==
	(int)RT_VIEW_FEATURE_TRANSIENT_PREVIEW,
	"qtcad feature constants must track neutral RT feature constants");
static_assert((int)QG_LEGACY_VIEW_EDIT_PREVIEW_BEGIN ==
	(int)RT_VIEW_EDIT_PREVIEW_BEGIN,
	"qtcad edit-preview constants must track neutral RT constants");
static_assert((int)QG_LEGACY_VIEW_EDIT_PREVIEW_DISCARD ==
	(int)RT_VIEW_EDIT_PREVIEW_DISCARD,
	"qtcad edit-preview constants must track neutral RT constants");

struct qg_legacy_view_bounds_update_state {
    rt_view_bounds_update_callback_bsg_t callback = nullptr;
};

struct qg_legacy_view_draw_observer_bridge {
    qg_legacy_view_draw_observer_callback_t callback = nullptr;
    void *client_data = nullptr;
    uintptr_t ged_token = 0;
};

struct qg_legacy_view_draw_transaction_storage {
    qg_legacy_view_draw_transaction_kind kind = QG_LEGACY_VIEW_DRAW_TXN_NONE;
    const char *path = nullptr;
    const char *new_path = nullptr;
    const char **paths = nullptr;
    int path_count = 0;
    const struct db_full_path *dfp = nullptr;
    qg_legacy_view *view = nullptr;
    int mode = -1;
    const qg_legacy_view_draw_appearance *appearance = nullptr;
    int appearance_draw_mode = -1;
    int autoview = 0;
    double value = 0.0;
    int stale_reason = GED_DRAW_STALE_SOURCE_CHANGED;
    int removed = 0;
    int redraw = 0;
};

struct qg_legacy_view_draw_result_storage {
    int status = 0;
    int affected_groups = 0;
    int affected_shapes = 0;
    int stale_count = 0;
    int redrawn_count = 0;
    int error_count = 0;
    uint64_t scene_revision_before = 0;
    uint64_t scene_revision_after = 0;
    struct bu_vls names = BU_VLS_INIT_ZERO;
    struct bu_vls errors = BU_VLS_INIT_ZERO;
};

struct qg_legacy_view_sketch_lines {
    void *ctx = nullptr;
};

struct qg_legacy_view_draw_appearance {
    struct ged_draw_appearance_settings settings;
};

static_assert(sizeof(qg_legacy_view_draw_transaction) >=
	sizeof(qg_legacy_view_draw_transaction_storage),
	"qtcad draw transaction opaque storage is too small");
static_assert(alignof(qg_legacy_view_draw_transaction) >=
	alignof(qg_legacy_view_draw_transaction_storage),
	"qtcad draw transaction opaque storage is under-aligned");

static qg_legacy_view_pick_result *
qg_pick_result_from_rt(void *result)
{
    return reinterpret_cast<qg_legacy_view_pick_result *>(result);
}

static qg_legacy_view_draw_transaction_storage *
qg_draw_transaction_storage(qg_legacy_view_draw_transaction *txn)
{
    return txn ?
	reinterpret_cast<qg_legacy_view_draw_transaction_storage *>(
		txn->storage) :
	nullptr;
}

static const qg_legacy_view_draw_transaction_storage *
qg_draw_transaction_storage_const(const qg_legacy_view_draw_transaction *txn)
{
    return txn ?
	reinterpret_cast<const qg_legacy_view_draw_transaction_storage *>(
		txn->storage) :
	nullptr;
}

static qg_legacy_view_draw_result_storage *
qg_draw_result_storage(qg_legacy_view_draw_transaction_result *result)
{
    return result ?
	reinterpret_cast<qg_legacy_view_draw_result_storage *>(result->impl) :
	nullptr;
}

static const qg_legacy_view_draw_result_storage *
qg_draw_result_storage_const(
	const qg_legacy_view_draw_transaction_result *result)
{
    return result ?
	reinterpret_cast<const qg_legacy_view_draw_result_storage *>(
		result->impl) :
	nullptr;
}

static rt_view_feature_ref
qg_feature_ref_to_rt(qg_legacy_view_feature_ref ref)
{
    rt_view_feature_ref rt_ref = { ref.token, ref.revision };
    return rt_ref;
}

static qg_legacy_view_feature_ref
qg_feature_ref_from_rt(rt_view_feature_ref ref)
{
    qg_legacy_view_feature_ref qg_ref = { ref.token, ref.revision };
    return qg_ref;
}

static enum rt_view_feature_family
qg_feature_family_to_rt(qg_legacy_view_feature_family family)
{
    switch (family) {
	case QG_LEGACY_VIEW_FEATURE_TRANSIENT_PREVIEW:
	    return RT_VIEW_FEATURE_TRANSIENT_PREVIEW;
	case QG_LEGACY_VIEW_FEATURE_UNKNOWN:
	default:
	    return RT_VIEW_FEATURE_UNKNOWN;
    }
}

static enum rt_view_edit_preview_event
qg_edit_preview_event_to_rt(qg_legacy_view_edit_preview_event event)
{
    switch (event) {
	case QG_LEGACY_VIEW_EDIT_PREVIEW_BEGIN:
	    return RT_VIEW_EDIT_PREVIEW_BEGIN;
	case QG_LEGACY_VIEW_EDIT_PREVIEW_UPDATE:
	    return RT_VIEW_EDIT_PREVIEW_UPDATE;
	case QG_LEGACY_VIEW_EDIT_PREVIEW_COMMIT:
	    return RT_VIEW_EDIT_PREVIEW_COMMIT;
	case QG_LEGACY_VIEW_EDIT_PREVIEW_CANCEL:
	    return RT_VIEW_EDIT_PREVIEW_CANCEL;
	case QG_LEGACY_VIEW_EDIT_PREVIEW_REPLACE_SOURCE:
	    return RT_VIEW_EDIT_PREVIEW_REPLACE_SOURCE;
	case QG_LEGACY_VIEW_EDIT_PREVIEW_DISCARD:
	    return RT_VIEW_EDIT_PREVIEW_DISCARD;
	default:
	    return RT_VIEW_EDIT_PREVIEW_UPDATE;
    }
}

static void
qg_edit_preview_callbacks_to_rt(
	struct rt_view_edit_preview_callbacks *rt_callbacks,
	const qg_legacy_view_edit_preview_callbacks *callbacks)
{
    if (!rt_callbacks || !callbacks)
	return;

    rt_callbacks->revision_cb = callbacks->revision_cb;
    rt_callbacks->update_cb = callbacks->update_cb;
    rt_callbacks->pick_cb = callbacks->pick_cb;
}

static ged_draw_transaction_kind
qg_draw_transaction_kind_to_ged(qg_legacy_view_draw_transaction_kind kind)
{
    switch (kind) {
	case QG_LEGACY_VIEW_DRAW_TXN_DRAW:
	    return GED_DRAW_TXN_DRAW;
	case QG_LEGACY_VIEW_DRAW_TXN_ERASE:
	    return GED_DRAW_TXN_ERASE;
	case QG_LEGACY_VIEW_DRAW_TXN_ERASE_PREFIX:
	    return GED_DRAW_TXN_ERASE_PREFIX;
	case QG_LEGACY_VIEW_DRAW_TXN_REDRAW:
	    return GED_DRAW_TXN_REDRAW;
	case QG_LEGACY_VIEW_DRAW_TXN_VISIBILITY:
	    return GED_DRAW_TXN_VISIBILITY;
	case QG_LEGACY_VIEW_DRAW_TXN_HIGHLIGHT:
	    return GED_DRAW_TXN_HIGHLIGHT;
	case QG_LEGACY_VIEW_DRAW_TXN_MATERIAL_CHANGED:
	    return GED_DRAW_TXN_MATERIAL_CHANGED;
	case QG_LEGACY_VIEW_DRAW_TXN_REFRESH_MATERIAL_COLORS:
	    return GED_DRAW_TXN_REFRESH_MATERIAL_COLORS;
	case QG_LEGACY_VIEW_DRAW_TXN_TRANSPARENCY:
	    return GED_DRAW_TXN_TRANSPARENCY;
	case QG_LEGACY_VIEW_DRAW_TXN_DEFAULT_DRAW_MODE:
	    return GED_DRAW_TXN_DEFAULT_DRAW_MODE;
	case QG_LEGACY_VIEW_DRAW_TXN_STALE_SOURCE:
	    return GED_DRAW_TXN_STALE_SOURCE;
	case QG_LEGACY_VIEW_DRAW_TXN_SOURCE_UPDATED:
	    return GED_DRAW_TXN_SOURCE_UPDATED;
	case QG_LEGACY_VIEW_DRAW_TXN_SOURCE_RENAMED:
	    return GED_DRAW_TXN_SOURCE_RENAMED;
	case QG_LEGACY_VIEW_DRAW_TXN_SOURCE_REFERENCES_REMOVED:
	    return GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED;
	case QG_LEGACY_VIEW_DRAW_TXN_CLEAR:
	    return GED_DRAW_TXN_CLEAR;
	case QG_LEGACY_VIEW_DRAW_TXN_CLEAR_SCOPE:
	    return GED_DRAW_TXN_CLEAR_SCOPE;
	case QG_LEGACY_VIEW_DRAW_TXN_TEARDOWN:
	    return GED_DRAW_TXN_TEARDOWN;
	case QG_LEGACY_VIEW_DRAW_TXN_NONE:
	default:
	    return GED_DRAW_TXN_NONE;
    }
}

static qg_legacy_view_draw_transaction_kind
qg_draw_transaction_kind_from_ged(ged_draw_transaction_kind kind)
{
    switch (kind) {
	case GED_DRAW_TXN_DRAW:
	    return QG_LEGACY_VIEW_DRAW_TXN_DRAW;
	case GED_DRAW_TXN_ERASE:
	    return QG_LEGACY_VIEW_DRAW_TXN_ERASE;
	case GED_DRAW_TXN_ERASE_PREFIX:
	    return QG_LEGACY_VIEW_DRAW_TXN_ERASE_PREFIX;
	case GED_DRAW_TXN_REDRAW:
	    return QG_LEGACY_VIEW_DRAW_TXN_REDRAW;
	case GED_DRAW_TXN_VISIBILITY:
	    return QG_LEGACY_VIEW_DRAW_TXN_VISIBILITY;
	case GED_DRAW_TXN_HIGHLIGHT:
	    return QG_LEGACY_VIEW_DRAW_TXN_HIGHLIGHT;
	case GED_DRAW_TXN_MATERIAL_CHANGED:
	    return QG_LEGACY_VIEW_DRAW_TXN_MATERIAL_CHANGED;
	case GED_DRAW_TXN_REFRESH_MATERIAL_COLORS:
	    return QG_LEGACY_VIEW_DRAW_TXN_REFRESH_MATERIAL_COLORS;
	case GED_DRAW_TXN_TRANSPARENCY:
	    return QG_LEGACY_VIEW_DRAW_TXN_TRANSPARENCY;
	case GED_DRAW_TXN_DEFAULT_DRAW_MODE:
	    return QG_LEGACY_VIEW_DRAW_TXN_DEFAULT_DRAW_MODE;
	case GED_DRAW_TXN_STALE_SOURCE:
	    return QG_LEGACY_VIEW_DRAW_TXN_STALE_SOURCE;
	case GED_DRAW_TXN_SOURCE_UPDATED:
	    return QG_LEGACY_VIEW_DRAW_TXN_SOURCE_UPDATED;
	case GED_DRAW_TXN_SOURCE_RENAMED:
	    return QG_LEGACY_VIEW_DRAW_TXN_SOURCE_RENAMED;
	case GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED:
	    return QG_LEGACY_VIEW_DRAW_TXN_SOURCE_REFERENCES_REMOVED;
	case GED_DRAW_TXN_CLEAR:
	    return QG_LEGACY_VIEW_DRAW_TXN_CLEAR;
	case GED_DRAW_TXN_CLEAR_SCOPE:
	    return QG_LEGACY_VIEW_DRAW_TXN_CLEAR_SCOPE;
	case GED_DRAW_TXN_TEARDOWN:
	    return QG_LEGACY_VIEW_DRAW_TXN_TEARDOWN;
	case GED_DRAW_TXN_NONE:
	default:
	    return QG_LEGACY_VIEW_DRAW_TXN_NONE;
    }
}

static void
qg_draw_transaction_from_ged(qg_legacy_view_draw_transaction *qtxn,
	const struct ged_draw_transaction *gtxn)
{
    if (!qtxn)
	return;

    qg_legacy_view_draw_transaction_init(qtxn,
	    QG_LEGACY_VIEW_DRAW_TXN_NONE, nullptr);
    qg_legacy_view_draw_transaction_storage *storage =
	qg_draw_transaction_storage(qtxn);
    if (!storage)
	return;
    if (!gtxn)
	return;

    storage->kind = qg_draw_transaction_kind_from_ged(gtxn->kind);
    storage->path = gtxn->path;
    storage->new_path = gtxn->new_path;
    storage->paths = gtxn->paths;
    storage->path_count = gtxn->path_count;
    storage->dfp = gtxn->dfp;
    storage->view = qg_legacy_view_from_context(gtxn->view);
    storage->mode = gtxn->mode;
    if (gtxn->appearance) {
	const struct ged_draw_appearance_settings *appearance =
	    (const struct ged_draw_appearance_settings *)gtxn->appearance;
	storage->appearance_draw_mode = appearance->draw_mode;
    }
    storage->autoview = gtxn->autoview;
    storage->value = gtxn->value;
    storage->stale_reason = gtxn->stale_reason;
    storage->removed = gtxn->removed;
    storage->redraw = gtxn->redraw;
}

static struct ged_draw_transaction
qg_draw_transaction_to_ged(const qg_legacy_view_draw_transaction *qtxn)
{
    const qg_legacy_view_draw_transaction_storage *storage =
	qg_draw_transaction_storage_const(qtxn);
    struct ged_draw_transaction gtxn = ged_draw_transaction_make(
	    storage ? qg_draw_transaction_kind_to_ged(storage->kind) :
	    GED_DRAW_TXN_NONE,
	    storage ? storage->path : nullptr);
    if (!storage)
	return gtxn;

    gtxn.new_path = storage->new_path;
    gtxn.paths = storage->paths;
    gtxn.path_count = storage->path_count;
    gtxn.dfp = storage->dfp;
    gtxn.view = qg_legacy_view_to_context(storage->view);
    gtxn.mode = storage->mode;
    gtxn.appearance = storage->appearance ? &storage->appearance->settings :
	nullptr;
    gtxn.autoview = storage->autoview;
    gtxn.value = storage->value;
    gtxn.stale_reason = (ged_draw_stale_reason)storage->stale_reason;
    gtxn.removed = storage->removed;
    gtxn.redraw = storage->redraw;
    return gtxn;
}

static void
qg_draw_result_copy_from_ged(qg_legacy_view_draw_transaction_result *qresult,
	const struct ged_draw_transaction_result *gresult)
{
    qg_legacy_view_draw_result_storage *storage =
	qg_draw_result_storage(qresult);
    if (!storage || !gresult)
	return;

    storage->status = gresult->status;
    storage->affected_groups = gresult->affected_groups;
    storage->affected_shapes = gresult->affected_shapes;
    storage->stale_count = gresult->stale_count;
    storage->redrawn_count = gresult->redrawn_count;
    storage->error_count = gresult->error_count;
    storage->scene_revision_before = gresult->scene_revision_before;
    storage->scene_revision_after = gresult->scene_revision_after;

    bu_vls_trunc(&storage->names, 0);
    if (BU_VLS_IS_INITIALIZED(&gresult->names))
	bu_vls_strcat(&storage->names, bu_vls_cstr(&gresult->names));

    bu_vls_trunc(&storage->errors, 0);
    if (BU_VLS_IS_INITIALIZED(&gresult->errors))
	bu_vls_strcat(&storage->errors, bu_vls_cstr(&gresult->errors));
}

static void
qg_draw_observer_bridge_callback(struct ged *gedp,
	const struct ged_draw_transaction *gtxn,
	const struct ged_draw_transaction_result *gresult,
	void *client_data)
{
    qg_legacy_view_draw_observer_bridge *bridge =
	static_cast<qg_legacy_view_draw_observer_bridge *>(client_data);
    if (!bridge || !bridge->callback)
	return;

    qg_legacy_view_draw_transaction qtxn;
    qg_draw_transaction_from_ged(&qtxn, gtxn);

    qg_legacy_view_draw_transaction_result qresult;
    qg_legacy_view_draw_result_init(&qresult);
    qg_draw_result_copy_from_ged(&qresult, gresult);

    bridge->callback(gedp, &qtxn, &qresult, bridge->client_data);

    qg_legacy_view_draw_result_free(&qresult);
}

static void *
qg_pick_result_to_rt(qg_legacy_view_pick_result *result)
{
    return reinterpret_cast<void *>(result);
}

static const void *
qg_pick_result_const_to_rt(const qg_legacy_view_pick_result *result)
{
    return reinterpret_cast<const void *>(result);
}

qg_legacy_view *
qg_legacy_view_local_create(const char *name)
{
    void *view_ctx = rt_view_context_create_bsg();
    if (name)
	rt_view_context_name_set_bsg(view_ctx, name);
    return qg_legacy_view_from_context(view_ctx);
}

void
qg_legacy_view_local_free(qg_legacy_view *view)
{
    /* Canvas-local views may outlive GED scene cleanup in tests and app
     * shutdown; keep the historical qtcad release behavior until the lower
     * retained view source is retired. */
    rt_view_context_release_storage_bsg(qg_legacy_view_to_context(view));
}

void
qg_legacy_view_local_destroy(qg_legacy_view *view)
{
    void *view_ctx = qg_legacy_view_to_context(view);
    if (!view_ctx)
	return;

    rt_view_context_free_bsg(view_ctx);
}

struct rt_edit *
qg_legacy_view_sketch_edit_create(qg_legacy_view *view,
	struct db_full_path *dfp,
	struct db_i *dbip,
	struct bn_tol *tol)
{
    void *view_ctx = qg_legacy_view_to_context(view);
    if (!view_ctx || !dfp || !dbip || !tol)
	return nullptr;
    return rt_edit_create_bsg(dfp, dbip, tol, view_ctx);
}

struct qg_legacy_view_sketch_lines *
qg_legacy_view_sketch_lines_create(qg_legacy_view *view,
	const char *name,
	unsigned char r,
	unsigned char g,
	unsigned char b)
{
    void *view_ctx = qg_legacy_view_to_context(view);
    if (!view_ctx || !name)
	return nullptr;

    struct qg_legacy_view_sketch_lines *lines = new qg_legacy_view_sketch_lines;
    lines->ctx = rt_view_context_line_set_create_bsg(view_ctx, name, r, g, b);

    return lines;
}

int
qg_legacy_view_sketch_lines_is_null(
	const struct qg_legacy_view_sketch_lines *lines)
{
    return (!lines || rt_view_line_set_context_is_null_bsg(lines->ctx));
}

int
qg_legacy_view_sketch_lines_set_points(
	struct qg_legacy_view_sketch_lines *lines,
	const point_t *points,
	const int *commands,
	size_t point_count)
{
    if (qg_legacy_view_sketch_lines_is_null(lines))
	return 0;
    return rt_view_line_set_context_set_points_bsg(lines->ctx, points,
	    commands, point_count);
}

void
qg_legacy_view_sketch_lines_destroy(
	struct qg_legacy_view_sketch_lines *lines)
{
    if (!lines)
	return;
    rt_view_line_set_context_destroy_bsg(lines->ctx);
    delete lines;
}

int
qg_legacy_view_dimensions_set(qg_legacy_view *view, int width, int height)
{
    return rt_view_context_dimensions_set_bsg(qg_legacy_view_to_context(view),
	    width, height);
}

int
qg_legacy_view_width_get(const qg_legacy_view *view)
{
    return rt_view_context_width_from_bsg(qg_legacy_view_to_context(view));
}

int
qg_legacy_view_height_get(const qg_legacy_view *view)
{
    return rt_view_context_height_from_bsg(qg_legacy_view_to_context(view));
}

int
qg_legacy_view_unit_conversion_set(qg_legacy_view *view,
	double local2base,
	double base2local)
{
    return rt_view_context_unit_conversion_set_bsg(qg_legacy_view_to_context(view),
	    local2base, base2local);
}

int
qg_legacy_view_name_set(qg_legacy_view *view, const char *name)
{
    return rt_view_context_name_set_bsg(qg_legacy_view_to_context(view), name);
}

qg_legacy_view *
qg_legacy_view_ged_active_get(struct ged *gedp)
{
    return qg_legacy_view_from_context(ged_view_active_ctx(gedp));
}

void
qg_legacy_view_ged_active_set(struct ged *gedp, qg_legacy_view *view)
{
    ged_view_active_ctx_set(gedp, qg_legacy_view_to_context(view));
}

struct db_i *
qg_legacy_view_ged_database(struct ged *gedp)
{
    return gedp ? gedp->dbip : nullptr;
}

int
qg_legacy_view_ged_view_set_add(struct ged *gedp, qg_legacy_view *view)
{
    if (!gedp || !view)
	return 0;

    return rt_view_set_context_add_bsg(ged_view_set_ctx(gedp),
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_ged_view_set_remove(struct ged *gedp, qg_legacy_view *view)
{
    if (!gedp || !view)
	return 0;

    return rt_view_set_context_remove_bsg(ged_view_set_ctx(gedp),
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_ged_view_set_attach(struct ged *gedp, qg_legacy_view *view)
{
    if (!gedp || !view)
	return 0;

    return rt_view_context_view_set_attach_bsg(qg_legacy_view_to_context(view),
	    ged_view_set_ctx(gedp));
}

qg_legacy_view *
qg_legacy_view_draw_transaction_view_get(
	const qg_legacy_view_draw_transaction *txn)
{
    const qg_legacy_view_draw_transaction_storage *storage =
	qg_draw_transaction_storage_const(txn);
    return storage ? storage->view : nullptr;
}

void
qg_legacy_view_draw_transaction_view_set(qg_legacy_view_draw_transaction *txn,
	qg_legacy_view *view)
{
    qg_legacy_view_draw_transaction_storage *storage =
	qg_draw_transaction_storage(txn);
    if (storage)
	storage->view = view;
}

void
qg_legacy_view_draw_transaction_init(qg_legacy_view_draw_transaction *txn,
	qg_legacy_view_draw_transaction_kind kind,
	const char *path)
{
    if (!txn)
	return;

    *txn = qg_legacy_view_draw_transaction{};
    qg_legacy_view_draw_transaction_storage *storage =
	new (static_cast<void *>(txn->storage))
	qg_legacy_view_draw_transaction_storage;
    if (!storage)
	return;

    storage->kind = kind;
    storage->path = path;
}

void
qg_legacy_view_draw_transaction_paths_set(
	qg_legacy_view_draw_transaction *txn,
	const char **paths,
	int path_count)
{
    qg_legacy_view_draw_transaction_storage *storage =
	qg_draw_transaction_storage(txn);
    if (!storage)
	return;

    storage->paths = paths;
    storage->path_count = (path_count > 0) ? path_count : 0;
}

void
qg_legacy_view_draw_transaction_autoview_set(
	qg_legacy_view_draw_transaction *txn,
	int autoview)
{
    qg_legacy_view_draw_transaction_storage *storage =
	qg_draw_transaction_storage(txn);
    if (storage)
	storage->autoview = autoview ? 1 : 0;
}

qg_legacy_view_draw_appearance *
qg_legacy_view_draw_appearance_create(int draw_mode)
{
    qg_legacy_view_draw_appearance *appearance =
	new qg_legacy_view_draw_appearance;
    appearance->settings = GED_DRAW_APPEARANCE_SETTINGS_INIT;
    if (draw_mode >= 0)
	appearance->settings.draw_mode = draw_mode;
    return appearance;
}

void
qg_legacy_view_draw_appearance_destroy(
	qg_legacy_view_draw_appearance *appearance)
{
    delete appearance;
}

void
qg_legacy_view_draw_transaction_appearance_set(
	qg_legacy_view_draw_transaction *txn,
	const qg_legacy_view_draw_appearance *appearance)
{
    qg_legacy_view_draw_transaction_storage *storage =
	qg_draw_transaction_storage(txn);
    if (!storage)
	return;

    storage->appearance = appearance;
    storage->appearance_draw_mode = appearance ? appearance->settings.draw_mode :
	-1;
}

int
qg_legacy_view_draw_path_state(struct ged *gedp, qg_legacy_view *view,
	const char *path, int mode)
{
    return ged_draw_path_state(gedp, qg_legacy_view_to_context(view), path, mode);
}

int
qg_legacy_view_draw_has_paths(struct ged *gedp, qg_legacy_view *view, int mode)
{
    return ged_draw_has_paths(gedp, qg_legacy_view_to_context(view), mode);
}

size_t
qg_legacy_view_draw_list_paths(struct ged *gedp, qg_legacy_view *view,
	int mode, int expanded, struct bu_vls *paths)
{
    return ged_draw_list_paths(gedp, qg_legacy_view_to_context(view), mode,
	    expanded, paths);
}

uint64_t
qg_legacy_view_draw_scene_revision(struct ged *gedp)
{
    return ged_draw_scene_revision(gedp);
}

int
qg_legacy_view_scene_attached(const qg_legacy_view *view)
{
    return rt_view_context_scene_attached_bsg(qg_legacy_view_to_context(view));
}

uintptr_t
qg_legacy_view_draw_observer_add(struct ged *gedp,
	qg_legacy_view_draw_observer_callback_t callback,
	void *client_data)
{
    if (!gedp || !callback)
	return 0;

    qg_legacy_view_draw_observer_bridge *bridge =
	new qg_legacy_view_draw_observer_bridge;
    bridge->callback = callback;
    bridge->client_data = client_data;
    bridge->ged_token = ged_draw_observer_add(gedp,
	    qg_draw_observer_bridge_callback, bridge);
    if (!bridge->ged_token) {
	delete bridge;
	return 0;
    }

    return reinterpret_cast<uintptr_t>(bridge);
}

int
qg_legacy_view_draw_observer_remove(struct ged *gedp, uintptr_t token)
{
    qg_legacy_view_draw_observer_bridge *bridge =
	reinterpret_cast<qg_legacy_view_draw_observer_bridge *>(token);
    if (!bridge)
	return 0;

    int ret = ged_draw_observer_remove(gedp, bridge->ged_token);
    if (ret) {
	delete bridge;
	return ret;
    }

    return 0;
}

int
qg_legacy_view_draw_redraw(struct ged *gedp)
{
    qg_legacy_view_draw_transaction txn;
    qg_legacy_view_draw_transaction_init(&txn,
	    QG_LEGACY_VIEW_DRAW_TXN_REDRAW, NULL);
    return qg_legacy_view_draw_transaction_apply(gedp, &txn, NULL);
}

int
qg_legacy_view_draw_transaction_has_view(
	const qg_legacy_view_draw_transaction *txn)
{
    const qg_legacy_view_draw_transaction_storage *storage =
	qg_draw_transaction_storage_const(txn);
    return (storage && storage->view) ? 1 : 0;
}

int
qg_legacy_view_draw_transaction_is_view_only(
	const qg_legacy_view_draw_transaction *txn)
{
    const qg_legacy_view_draw_transaction_storage *storage =
	qg_draw_transaction_storage_const(txn);
    return (storage &&
	    (storage->kind == QG_LEGACY_VIEW_DRAW_TXN_MATERIAL_CHANGED ||
	     storage->kind == QG_LEGACY_VIEW_DRAW_TXN_REFRESH_MATERIAL_COLORS)) ?
	1 : 0;
}

qg_legacy_view_draw_transaction_kind
qg_legacy_view_draw_transaction_kind_get(
	const qg_legacy_view_draw_transaction *txn)
{
    const qg_legacy_view_draw_transaction_storage *storage =
	qg_draw_transaction_storage_const(txn);
    if (!storage)
	return QG_LEGACY_VIEW_DRAW_TXN_NONE;

    return storage->kind;
}

const char *
qg_legacy_view_draw_transaction_path(const qg_legacy_view_draw_transaction *txn)
{
    const qg_legacy_view_draw_transaction_storage *storage =
	qg_draw_transaction_storage_const(txn);
    return storage ? storage->path : nullptr;
}

int
qg_legacy_view_draw_transaction_path_count(
	const qg_legacy_view_draw_transaction *txn)
{
    const qg_legacy_view_draw_transaction_storage *storage =
	qg_draw_transaction_storage_const(txn);
    return (storage && storage->paths && storage->path_count > 0) ?
	storage->path_count : 0;
}

const char *
qg_legacy_view_draw_transaction_path_at(
	const qg_legacy_view_draw_transaction *txn,
	int index)
{
    const qg_legacy_view_draw_transaction_storage *storage =
	qg_draw_transaction_storage_const(txn);
    if (!storage || !storage->paths || index < 0 ||
	    index >= storage->path_count)
	return nullptr;

    return storage->paths[index];
}

int
qg_legacy_view_draw_transaction_mode(struct ged *gedp,
	const qg_legacy_view_draw_transaction *txn)
{
    const qg_legacy_view_draw_transaction_storage *storage =
	qg_draw_transaction_storage_const(txn);
    int mode = QG_LEGACY_VIEW_DRAW_MODE_ANY;
    if (storage && storage->appearance_draw_mode >= 0)
	mode = storage->appearance_draw_mode;
    else if (storage && storage->kind == QG_LEGACY_VIEW_DRAW_TXN_DRAW &&
	    storage->mode >= 0)
	mode = storage->mode;
    if (mode < 0 && gedp)
	mode = ged_draw_default_mode(gedp);
    return mode;
}

void
qg_legacy_view_draw_result_init(qg_legacy_view_draw_transaction_result *result)
{
    if (!result)
	return;

    *result = qg_legacy_view_draw_transaction_result{};
    qg_legacy_view_draw_result_storage *storage =
	new qg_legacy_view_draw_result_storage;
    result->impl = reinterpret_cast<uintptr_t>(storage);
}

void
qg_legacy_view_draw_result_free(qg_legacy_view_draw_transaction_result *result)
{
    if (!result)
	return;

    qg_legacy_view_draw_result_storage *storage =
	qg_draw_result_storage(result);
    if (storage) {
	if (BU_VLS_IS_INITIALIZED(&storage->names))
	    bu_vls_free(&storage->names);
	if (BU_VLS_IS_INITIALIZED(&storage->errors))
	    bu_vls_free(&storage->errors);
	delete storage;
    }
    *result = qg_legacy_view_draw_transaction_result{};
}

int
qg_legacy_view_draw_transaction_apply(struct ged *gedp,
	const qg_legacy_view_draw_transaction *txn,
	qg_legacy_view_draw_transaction_result *result)
{
    struct ged_draw_transaction gtxn = qg_draw_transaction_to_ged(txn);
    if (!result)
	return ged_draw_apply_transaction(gedp, &gtxn, NULL);

    struct ged_draw_transaction_result gresult;
    ged_draw_transaction_result_init(&gresult);
    int ret = ged_draw_apply_transaction(gedp, &gtxn, &gresult);
    qg_draw_result_copy_from_ged(result, &gresult);
    ged_draw_transaction_result_free(&gresult);
    return ret;
}

int
qg_legacy_view_draw_result_status(
	const qg_legacy_view_draw_transaction_result *result)
{
    const qg_legacy_view_draw_result_storage *storage =
	qg_draw_result_storage_const(result);
    return storage ? storage->status : 0;
}

uint64_t
qg_legacy_view_draw_result_scene_revision_after(
	const qg_legacy_view_draw_transaction_result *result)
{
    const qg_legacy_view_draw_result_storage *storage =
	qg_draw_result_storage_const(result);
    return storage ? storage->scene_revision_after : 0;
}

const char *
qg_legacy_view_draw_result_names(
	const qg_legacy_view_draw_transaction_result *result)
{
    const qg_legacy_view_draw_result_storage *storage =
	qg_draw_result_storage_const(result);
    if (!storage || !BU_VLS_IS_INITIALIZED(&storage->names) ||
	    !bu_vls_strlen(&storage->names))
	return nullptr;

    return bu_vls_cstr(&storage->names);
}

const char *
qg_legacy_view_draw_result_errors(
	const qg_legacy_view_draw_transaction_result *result)
{
    const qg_legacy_view_draw_result_storage *storage =
	qg_draw_result_storage_const(result);
    if (!storage || !BU_VLS_IS_INITIALIZED(&storage->errors) ||
	    !bu_vls_strlen(&storage->errors))
	return nullptr;

    return bu_vls_cstr(&storage->errors);
}

int
qg_legacy_view_draw_paths_apply(struct ged *gedp,
	qg_legacy_view *view,
	const char **paths,
	int path_count,
	int autoview)
{
    qg_legacy_view_draw_transaction txn;
    qg_legacy_view_draw_transaction_init(&txn,
	    QG_LEGACY_VIEW_DRAW_TXN_DRAW, NULL);
    qg_legacy_view_draw_transaction_view_set(&txn, view);
    qg_legacy_view_draw_transaction_autoview_set(&txn, autoview);
    qg_legacy_view_draw_transaction_paths_set(&txn, paths, path_count);

    qg_legacy_view_draw_transaction_result result;
    qg_legacy_view_draw_result_init(&result);
    int ret = qg_legacy_view_draw_transaction_apply(gedp, &txn, &result);
    qg_legacy_view_draw_result_free(&result);
    return ret;
}

int
qg_legacy_view_draw_erase_path_apply(struct ged *gedp,
	qg_legacy_view *view,
	const char *path)
{
    qg_legacy_view_draw_transaction txn;
    qg_legacy_view_draw_transaction_init(&txn,
	    QG_LEGACY_VIEW_DRAW_TXN_ERASE, path);
    qg_legacy_view_draw_transaction_view_set(&txn, view);

    qg_legacy_view_draw_transaction_result result;
    qg_legacy_view_draw_result_init(&result);
    int ret = qg_legacy_view_draw_transaction_apply(gedp, &txn, &result);
    qg_legacy_view_draw_result_free(&result);
    return ret;
}

int
qg_legacy_view_draw_highlight_shape_by_name(struct ged *gedp, const char *name)
{
    return name ? ged_draw_highlight_shape_ref_by_name(gedp, name).token != 0 :
	0;
}

int
qg_legacy_view_draw_highlight_clear(struct ged *gedp)
{
    if (!gedp)
	return 0;
    ged_draw_set_highlighted_shape_ref(gedp, GED_DRAW_SHAPE_REF_NULL);
    return 1;
}

int
qg_legacy_view_refresh_request(qg_legacy_view *view, unsigned flags)
{
    return rt_view_context_refresh_request_bsg(qg_legacy_view_to_context(view),
	    flags);
}

uint32_t
qg_legacy_view_refresh_consume(qg_legacy_view *view)
{
    return rt_view_context_refresh_consume_bsg(qg_legacy_view_to_context(view));
}

int
qg_legacy_view_refresh_complete(qg_legacy_view *view)
{
    return rt_view_context_refresh_complete_bsg(qg_legacy_view_to_context(view));
}

struct qg_legacy_view_bounds_update_state *
qg_legacy_view_bounds_update_suspend(qg_legacy_view *view)
{
    void *view_ctx = qg_legacy_view_to_context(view);
    rt_view_bounds_update_callback_bsg_t callback =
	rt_view_context_bounds_update_callback_from_bsg(view_ctx);
    if (!callback)
	return nullptr;

    auto *state = new qg_legacy_view_bounds_update_state;
    state->callback = callback;
    rt_view_context_bounds_update_callback_set_bsg(view_ctx, nullptr);
    return state;
}

int
qg_legacy_view_bounds_update_restore(qg_legacy_view *view,
	struct qg_legacy_view_bounds_update_state *state,
	int refresh_bounds)
{
    if (!state)
	return 0;

    int restored = 0;
    void *view_ctx = qg_legacy_view_to_context(view);
    if (view_ctx) {
	rt_view_context_bounds_update_callback_set_bsg(view_ctx, state->callback);
	if (refresh_bounds)
	    rt_view_context_bounds_update_callback_call_bsg(view_ctx);
	restored = 1;
    }

    delete state;
    return restored;
}

int
qg_legacy_view_adjust(qg_legacy_view *view,
	int dx,
	int dy,
	point_t keypoint,
	int mode,
	unsigned long long flags)
{
    return rt_view_context_adjust_bsg(qg_legacy_view_to_context(view), dx, dy,
	    keypoint, mode, flags);
}

double
qg_legacy_view_local2base_get(const qg_legacy_view *view)
{
    return rt_view_context_local2base_from_bsg(qg_legacy_view_to_context(view));
}

double
qg_legacy_view_base2local_get(const qg_legacy_view *view)
{
    return rt_view_context_base2local_from_bsg(qg_legacy_view_to_context(view));
}

void
qg_legacy_view_info_get(const qg_legacy_view *view, struct rt_view_info *info)
{
    rt_view_context_info_from_bsg(info, qg_legacy_view_to_context(view));
}

const char *
qg_legacy_view_name_get(const qg_legacy_view *view)
{
    return rt_view_context_name_from_bsg(qg_legacy_view_to_context(view));
}

int
qg_legacy_view_aet_get(const qg_legacy_view *view, vect_t aet)
{
    return rt_view_context_aet_from_bsg(aet, qg_legacy_view_to_context(view));
}

int
qg_legacy_view_center_get(const qg_legacy_view *view, mat_t center)
{
    return rt_view_context_center_from_bsg(center,
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_rotation_get(const qg_legacy_view *view, mat_t rotation)
{
    return rt_view_context_rotation_from_bsg(rotation,
	    qg_legacy_view_to_context(view));
}

fastf_t
qg_legacy_view_scale_get(const qg_legacy_view *view)
{
    return rt_view_context_scale_from_bsg(qg_legacy_view_to_context(view));
}

int
qg_legacy_view_scale_set(qg_legacy_view *view, fastf_t scale)
{
    return rt_view_context_scale_set_bsg(qg_legacy_view_to_context(view), scale);
}

fastf_t
qg_legacy_view_perspective_get(const qg_legacy_view *view)
{
    return rt_view_context_perspective_from_bsg(qg_legacy_view_to_context(view));
}

fastf_t *
qg_legacy_view_scale_storage_get(qg_legacy_view *view)
{
    return rt_view_context_scale_storage_from_bsg(qg_legacy_view_to_context(view));
}

int
qg_legacy_view_display_manager_set(qg_legacy_view *view, qg_legacy_dm *dmp)
{
    return rt_view_context_display_manager_set_bsg(qg_legacy_view_to_context(view),
	    qg_legacy_dm_to_dm(dmp));
}

qg_legacy_dm *
qg_legacy_view_display_manager_get(qg_legacy_view *view)
{
    return qg_legacy_dm_from_dm(reinterpret_cast<struct dm *>(
		rt_view_context_display_manager_from_bsg(
		    qg_legacy_view_to_context(view))));
}

int
qg_legacy_view_dm_bind(qg_legacy_view *view, qg_legacy_dm *dmp)
{
    if (!view || !dmp)
	return 0;

    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    dm_set_vp(dm, qg_legacy_view_scale_storage_get(view));
    return qg_legacy_view_display_manager_set(view, dmp);
}

int
qg_legacy_view_dm_sync_dimensions(qg_legacy_view *view, qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!view || !dm)
	return 0;

    return qg_legacy_view_dimensions_set(view, dm_get_width(dm),
	    dm_get_height(dm));
}

unsigned long long
qg_legacy_view_dm_hash(qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_hash(dm) : 0ULL;
}

int
qg_legacy_view_dm_native_repaint_pending_get(qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_get_native_repaint_pending(dm) : 0;
}

void
qg_legacy_view_dm_native_repaint_pending_set(qg_legacy_dm *dmp, int pending)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (dm)
	dm_set_native_repaint_pending(dm, pending);
}

int
qg_legacy_view_dm_configure_window(qg_legacy_dm *dmp, int force)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_configure_win(dm, force) : 0;
}

int
qg_legacy_view_dm_width_get(qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_get_width(dm) : 0;
}

int
qg_legacy_view_dm_height_get(qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_get_height(dm) : 0;
}

void
qg_legacy_view_dm_dimensions_set(qg_legacy_dm *dmp, int width, int height)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!dm)
	return;

    dm_set_width(dm, width);
    dm_set_height(dm, height);
}

int
qg_legacy_view_dm_framebuffer_setup_existing(qg_legacy_fb *ifp,
	qg_legacy_dm *dmp)
{
    struct fb *fb = qg_legacy_fb_to_fb(ifp);
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!fb || !dm)
	return 0;

    struct fb_platform_specific *fbps = fb_get_platform_specific(FB_QTGL_MAGIC);
    if (!fbps)
	return 0;

    fbps->data = (void *)dm;
    fb_setup_existing(fb, dm_get_width(dm), dm_get_height(dm), fbps);
    fb_put_platform_specific(fbps);
    return 1;
}

qg_legacy_dm *
qg_legacy_view_dm_open_qtgl(void *context)
{
    const char *acmd = "attach";
    return qg_legacy_dm_from_dm(dm_open(context, nullptr, "qtgl", 1, &acmd));
}

qg_legacy_dm *
qg_legacy_view_dm_open_swrast(qg_legacy_view *view, void *context)
{
    const char *acmd = "attach";
    struct dm *dmp = dm_open(qg_legacy_view_dm_open_context(view), nullptr,
	    "swrast", 1, &acmd);
    if (dmp)
	dm_set_udata(dmp, context);
    return qg_legacy_dm_from_dm(dmp);
}

int
qg_legacy_view_dm_close(qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_close(dm) : 0;
}

int
qg_legacy_view_dm_setup_qtgl(qg_legacy_dm *dmp, fastf_t zmin, fastf_t zmax)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!dm)
	return 0;

    dm_set_pathname(dm, "QTDM");
    fastf_t windowbounds[6] = { -1, 1, -1, 1, zmin, zmax };
    dm_set_win_bounds(dm, windowbounds);
    return 1;
}

int
qg_legacy_view_dm_setup_swrast(qg_legacy_dm *dmp, fastf_t zmin, fastf_t zmax)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!dm)
	return 0;

    dm_configure_win(dm, 0);
    dm_set_pathname(dm, "SWDM");
    dm_set_zbuffer(dm, 1);
    fastf_t windowbounds[6] = { -1, 1, -1, 1, zmin, zmax };
    dm_set_win_bounds(dm, windowbounds);
    return 1;
}

int
qg_legacy_view_dm_background_get(qg_legacy_dm *dmp, unsigned char bg1[3],
	unsigned char bg2[3])
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!dm || !bg1 || !bg2)
	return 0;

    unsigned char *dm_bg1 = nullptr;
    unsigned char *dm_bg2 = nullptr;
    dm_get_bg(&dm_bg1, &dm_bg2, dm);
    if (!dm_bg1 || !dm_bg2)
	return 0;

    bg1[0] = dm_bg1[0];
    bg1[1] = dm_bg1[1];
    bg1[2] = dm_bg1[2];
    bg2[0] = dm_bg2[0];
    bg2[1] = dm_bg2[1];
    bg2[2] = dm_bg2[2];
    return 1;
}

int
qg_legacy_view_dm_background_set(qg_legacy_dm *dmp, const unsigned char bg1[3],
	const unsigned char bg2[3])
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!dm || !bg1 || !bg2)
	return 0;

    dm_set_bg(dm, bg1[0], bg1[1], bg1[2], bg2[0], bg2[1], bg2[2]);
    return 1;
}

int
qg_legacy_view_dm_background_restore(qg_legacy_dm *dmp)
{
    unsigned char bg1[3] = {0, 0, 0};
    unsigned char bg2[3] = {0, 0, 0};
    if (!qg_legacy_view_dm_background_get(dmp, bg1, bg2))
	return 0;

    return qg_legacy_view_dm_background_set(dmp, bg1, bg2);
}

int
qg_legacy_view_dm_draw_begin(qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_draw_begin(dm) : 0;
}

int
qg_legacy_view_dm_draw_end(qg_legacy_dm *dmp)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return dm ? dm_draw_end(dm) : 0;
}

int
qg_legacy_view_dm_display_image_get(qg_legacy_dm *dmp, unsigned char **image,
	int copy, int release)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    return (dm && image) ? dm_get_display_image(dm, image, copy, release) : 1;
}

void
qg_legacy_view_dm_display_image_release(unsigned char *image)
{
    if (image)
	bu_free(image, "copy of backend image");
}

int
qg_legacy_view_dm_load_current_model2view(qg_legacy_dm *dmp,
	const qg_legacy_view *view, int which_eye_or_mode)
{
    struct dm *dm = qg_legacy_dm_to_dm(dmp);
    if (!dm || !view)
	return 0;

    mat_t model2view;
    if (!qg_legacy_view_model2view_get(view, model2view))
	return 0;

    dm_loadmatrix(dm, model2view, which_eye_or_mode);
    return 1;
}

qg_legacy_fb *
qg_legacy_view_framebuffer_raw_create(const char *type)
{
    struct fb *ifp = fb_raw(type);
    if (ifp)
	fb_set_standalone(ifp, 0);
    return qg_legacy_fb_from_fb(ifp);
}

qg_legacy_fb *
qg_legacy_view_framebuffer_handle_from_raw(void *ifp)
{
    return reinterpret_cast<qg_legacy_fb *>(ifp);
}

QgGL *
qg_legacy_view_framebuffer_qtgl_canvas_get(qg_legacy_fb *ifp)
{
    struct fb *fb = qg_legacy_fb_to_fb(ifp);
    struct dm *dm = fb ? fb_get_dm(fb) : nullptr;
    return dm ? reinterpret_cast<QgGL *>(dm_get_ctx(dm)) : nullptr;
}

QgSW *
qg_legacy_view_framebuffer_swrast_canvas_get(qg_legacy_fb *ifp)
{
    struct fb *fb = qg_legacy_fb_to_fb(ifp);
    struct dm *dm = fb ? fb_get_dm(fb) : nullptr;
    return dm ? reinterpret_cast<QgSW *>(dm_get_udata(dm)) : nullptr;
}

int
qg_legacy_view_framebuffer_release(qg_legacy_fb *ifp, int initialized)
{
    struct fb *fb = qg_legacy_fb_to_fb(ifp);
    if (!fb || fb_get_standalone(fb))
	return 0;

    if (initialized)
	return fb_close_existing(fb);

    fb_put(fb);
    return 1;
}

int
qg_legacy_view_framebuffer_standalone_get(qg_legacy_fb *ifp)
{
    struct fb *fb = qg_legacy_fb_to_fb(ifp);
    return fb ? fb_get_standalone(fb) : 0;
}

int
qg_legacy_view_framebuffer_configure(qg_legacy_fb *ifp, int width, int height)
{
    struct fb *fb = qg_legacy_fb_to_fb(ifp);
    return fb ? fb_configure_window(fb, width, height) : 0;
}

void *
qg_legacy_view_dm_open_context(qg_legacy_view *view)
{
    return qg_legacy_view_to_context(view);
}

int
qg_legacy_view_dm_draw(qg_legacy_view *view)
{
    void *view_ctx = qg_legacy_view_to_context(view);
    if (!view_ctx)
	return 0;

    dm_draw_objs(view_ctx);
    return 1;
}

const char *
qg_legacy_view_dm_init_messages(void)
{
    return dm_init_msgs();
}

int
qg_legacy_view_feature_ref_is_null(qg_legacy_view_feature_ref ref)
{
    return rt_view_feature_ref_is_null_bsg(qg_feature_ref_to_rt(ref));
}

int
qg_legacy_view_edit_preview_publish_event(qg_legacy_view *view,
	qg_legacy_view_feature_ref feature,
	qg_legacy_view_edit_preview_event event,
	const char *source_path)
{
    return rt_view_context_edit_preview_publish_event_bsg(
	    qg_legacy_view_to_context(view),
	    qg_feature_ref_to_rt(feature),
	    qg_edit_preview_event_to_rt(event),
	    source_path);
}

qg_legacy_view_feature_ref
qg_legacy_view_feature_overlay_ensure(qg_legacy_view *view,
	const char *name,
	const void *owner,
	void *preview_ctx,
	const qg_legacy_view_edit_preview_callbacks *callbacks,
	const char *source_path)
{
    struct rt_view_edit_preview_callbacks rt_callbacks =
	RT_VIEW_EDIT_PREVIEW_CALLBACKS_INIT;
    const struct rt_view_edit_preview_callbacks *rt_callbacksp = nullptr;
    if (callbacks) {
	qg_edit_preview_callbacks_to_rt(&rt_callbacks, callbacks);
	rt_callbacksp = &rt_callbacks;
    }

    return qg_feature_ref_from_rt(rt_view_context_feature_overlay_ensure_bsg(
		qg_legacy_view_to_context(view),
		name,
		owner,
		preview_ctx,
		rt_callbacksp,
		source_path));
}

qg_legacy_view_feature_ref
qg_legacy_view_feature_label_ensure(qg_legacy_view *view,
	const char *name,
	const void *owner)
{
    return qg_feature_ref_from_rt(rt_view_context_feature_label_ensure_bsg(
			qg_legacy_view_to_context(view), name, owner));
}

int
qg_legacy_view_feature_remove(qg_legacy_view *view, const char *name)
{
    return rt_view_context_feature_remove_bsg(qg_legacy_view_to_context(view),
	    name);
}

void
qg_legacy_view_feature_set_view(qg_legacy_view_feature_ref ref,
	qg_legacy_view *view)
{
    rt_view_feature_set_context_bsg(qg_feature_ref_to_rt(ref),
	    qg_legacy_view_to_context(view));
}

void
qg_legacy_view_feature_set_visible(qg_legacy_view_feature_ref ref,
	int visible)
{
    rt_view_feature_set_visible_bsg(qg_feature_ref_to_rt(ref), visible);
}

void
qg_legacy_view_feature_set_color(qg_legacy_view_feature_ref ref,
	int r,
	int g,
	int b)
{
    rt_view_feature_set_color_bsg(qg_feature_ref_to_rt(ref), r, g, b);
}

int
qg_legacy_view_feature_touch(qg_legacy_view_feature_ref ref)
{
    return rt_view_feature_touch_bsg(qg_feature_ref_to_rt(ref));
}

int
qg_legacy_view_feature_labels_replace(qg_legacy_view_feature_ref ref,
	const qg_legacy_view_feature_label *labels,
	size_t label_count)
{
    if (qg_legacy_view_feature_ref_is_null(ref))
	return 0;
    if (label_count && !labels)
	return 0;

    struct rt_view_feature_label *rt_labels = nullptr;
    if (label_count > 0) {
	rt_labels = (struct rt_view_feature_label *)bu_calloc(label_count,
		sizeof(struct rt_view_feature_label),
		"qtcad edit labels");
	for (size_t i = 0; i < label_count; i++) {
	    rt_labels[i].text = labels[i].text;
	    VMOVE(rt_labels[i].point, labels[i].point);
	    rt_labels[i].color_valid = labels[i].color_valid;
	    rt_labels[i].color[0] = labels[i].color[0];
	    rt_labels[i].color[1] = labels[i].color[1];
	    rt_labels[i].color[2] = labels[i].color[2];
	}
    }

    int ret = rt_view_feature_labels_replace_bsg(qg_feature_ref_to_rt(ref),
	    rt_labels, label_count);
    if (rt_labels)
	bu_free(rt_labels, "qtcad edit labels");
    return ret;
}

int
qg_legacy_view_feature_points_replace(qg_legacy_view_feature_ref ref,
	qg_legacy_view_feature_family family,
	const point_t *points,
	const int *cmds,
	size_t point_count)
{
    return rt_view_feature_points_replace_bsg(qg_feature_ref_to_rt(ref),
	    qg_feature_family_to_rt(family),
	    point_count ? points : NULL,
	    point_count ? cmds : NULL,
	    point_count);
}

int
qg_legacy_view_feature_clear_geometry(qg_legacy_view_feature_ref ref)
{
    return rt_view_feature_clear_geometry_bsg(qg_feature_ref_to_rt(ref));
}

unsigned long long
qg_legacy_view_hash(const qg_legacy_view *view)
{
    return rt_view_context_hash_bsg(qg_legacy_view_to_context(view));
}

int
qg_legacy_view_unique_object_name(struct bu_vls *name, const char *seed,
	qg_legacy_view *view)
{
    return rt_view_context_unique_object_name_bsg(name, seed,
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_framebuffer_mode_get(const qg_legacy_view *view)
{
    return rt_view_context_framebuffer_mode_from_bsg(
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_framebuffer_mode_set(qg_legacy_view *view, int mode)
{
    return rt_view_context_framebuffer_mode_set_bsg(
	    qg_legacy_view_to_context(view), mode);
}

int
qg_legacy_view_lod_policy_get(const qg_legacy_view *view,
	struct rt_view_lod_policy *policy)
{
    return rt_view_context_lod_policy_from_bsg(policy,
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_lod_policy_apply(qg_legacy_view *view,
	const struct rt_view_lod_policy *policy)
{
    return rt_view_context_lod_policy_apply_bsg(qg_legacy_view_to_context(view),
	    policy);
}

int
qg_legacy_view_lod_policy_copy(qg_legacy_view *dst, const qg_legacy_view *src)
{
    return rt_view_context_lod_policy_copy_bsg(qg_legacy_view_to_context(dst),
	    qg_legacy_view_to_context(src));
}

int
qg_legacy_view_autoview_default(qg_legacy_view *view, int all_view_objs)
{
    return rt_view_context_autoview_bsg(qg_legacy_view_to_context(view),
	    RT_VIEW_AUTOVIEW_SCALE_DEFAULT, all_view_objs);
}

int
qg_legacy_view_lod_bounds_update(qg_legacy_view *view)
{
    return rt_view_context_lod_bounds_update_bsg(qg_legacy_view_to_context(view));
}

int
qg_legacy_view_adc_state_get(const qg_legacy_view *view,
	struct rt_view_adc_state *state)
{
    return rt_view_context_adc_state_from_bsg(state,
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_adc_state_set(qg_legacy_view *view,
	const struct rt_view_adc_state *state)
{
    return rt_view_context_adc_state_set_bsg(qg_legacy_view_to_context(view),
	    state);
}

int
qg_legacy_view_center_dot_state_get(const qg_legacy_view *view,
	struct rt_view_other_state *state)
{
    return rt_view_context_center_dot_state_from_bsg(state,
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_center_dot_state_set(qg_legacy_view *view,
	const struct rt_view_other_state *state)
{
    return rt_view_context_center_dot_state_set_bsg(
	    qg_legacy_view_to_context(view), state);
}

int
qg_legacy_view_grid_state_get(const qg_legacy_view *view,
	struct rt_view_grid_state *state)
{
    return rt_view_context_grid_state_from_bsg(state,
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_grid_state_set(qg_legacy_view *view,
	const struct rt_view_grid_state *state)
{
    return rt_view_context_grid_state_set_bsg(qg_legacy_view_to_context(view),
	    state);
}

int
qg_legacy_view_model_axes_state_get(const qg_legacy_view *view,
	struct rt_view_axes_state *state)
{
    return rt_view_context_model_axes_state_from_bsg(state,
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_model_axes_state_set(qg_legacy_view *view,
	const struct rt_view_axes_state *state)
{
    return rt_view_context_model_axes_state_set_bsg(
	    qg_legacy_view_to_context(view), state);
}

int
qg_legacy_view_scale_overlay_state_get(const qg_legacy_view *view,
	struct rt_view_other_state *state)
{
    return rt_view_context_scale_overlay_state_from_bsg(state,
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_scale_overlay_state_set(qg_legacy_view *view,
	const struct rt_view_other_state *state)
{
    return rt_view_context_scale_overlay_state_set_bsg(
	    qg_legacy_view_to_context(view), state);
}

int
qg_legacy_view_view_axes_state_get(const qg_legacy_view *view,
	struct rt_view_axes_state *state)
{
    return rt_view_context_view_axes_state_from_bsg(state,
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_view_axes_state_set(qg_legacy_view *view,
	const struct rt_view_axes_state *state)
{
    return rt_view_context_view_axes_state_set_bsg(
	    qg_legacy_view_to_context(view), state);
}

int
qg_legacy_view_params_state_get(const qg_legacy_view *view,
	struct rt_view_params_state *state)
{
    return rt_view_context_params_state_from_bsg(state,
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_params_state_set(qg_legacy_view *view,
	const struct rt_view_params_state *state)
{
    return rt_view_context_params_state_set_bsg(qg_legacy_view_to_context(view),
	    state);
}

int
qg_legacy_view_snap_source_view_only_set(qg_legacy_view *view)
{
    return rt_view_context_snap_source_flags_set_bsg(qg_legacy_view_to_context(view),
	    RT_VIEW_SNAP_VIEW_BSG);
}

int
qg_legacy_view_snap_source_db_set(qg_legacy_view *view)
{
    return rt_view_context_snap_source_flags_set_bsg(qg_legacy_view_to_context(view),
	    RT_VIEW_SNAP_DB_BSG);
}

int
qg_legacy_view_snap_lines_get(const qg_legacy_view *view)
{
    return rt_view_context_snap_lines_from_bsg(qg_legacy_view_to_context(view));
}

int
qg_legacy_view_snap_lines_set(qg_legacy_view *view, int enabled)
{
    return rt_view_context_snap_lines_set_bsg(qg_legacy_view_to_context(view),
	    enabled);
}

int
qg_legacy_view_db_snap_enabled(const qg_legacy_view *view)
{
    const void *view_ctx = qg_legacy_view_to_context(view);
    if (!view_ctx || !rt_view_context_snap_lines_from_bsg(view_ctx))
	return 0;

    int flags = rt_view_context_snap_source_flags_from_bsg(view_ctx);
    if (flags == RT_VIEW_SNAP_TCL_BSG)
	return 0;

    return !flags || (flags & RT_VIEW_SNAP_DB_BSG);
}

unsigned
qg_legacy_view_snap_kind_mask_get(const qg_legacy_view *view)
{
    unsigned kinds = 0;
    unsigned long long bsg_kinds =
	rt_view_context_snap_kind_mask_from_bsg(qg_legacy_view_to_context(view));

    if (bsg_kinds & RT_VIEW_SNAP_KIND_GRID_BSG)
	kinds |= QG_LEGACY_VIEW_SNAP_KIND_GRID;
    if (bsg_kinds & RT_VIEW_SNAP_KIND_ENDPOINT_BSG)
	kinds |= QG_LEGACY_VIEW_SNAP_KIND_ENDPOINT;
    if (bsg_kinds & RT_VIEW_SNAP_KIND_MIDPOINT_BSG)
	kinds |= QG_LEGACY_VIEW_SNAP_KIND_MIDPOINT;
    if (bsg_kinds & RT_VIEW_SNAP_KIND_INTERSECTION_BSG)
	kinds |= QG_LEGACY_VIEW_SNAP_KIND_INTERSECTION;
    if (bsg_kinds & RT_VIEW_SNAP_KIND_PERPENDICULAR_BSG)
	kinds |= QG_LEGACY_VIEW_SNAP_KIND_PERPENDICULAR;
    if (bsg_kinds & RT_VIEW_SNAP_KIND_TANGENT_BSG)
	kinds |= QG_LEGACY_VIEW_SNAP_KIND_TANGENT;
    if (bsg_kinds & RT_VIEW_SNAP_KIND_OVERLAY_HANDLE_BSG)
	kinds |= QG_LEGACY_VIEW_SNAP_KIND_OVERLAY_HANDLE;

    return kinds;
}

fastf_t
qg_legacy_view_size_get(const qg_legacy_view *view)
{
    return rt_view_context_size_from_bsg(qg_legacy_view_to_context(view));
}

int
qg_legacy_view_size_set(qg_legacy_view *view, fastf_t size)
{
    return rt_view_context_size_set_bsg(qg_legacy_view_to_context(view), size);
}

double
qg_legacy_view_snap_tolerance_factor_get(const qg_legacy_view *view)
{
    return rt_view_context_snap_tolerance_factor_from_bsg(
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_snap_exclude_clear(qg_legacy_view *view)
{
    return rt_view_context_snap_exclude_feature_clear_bsg(
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_center_vec_set(qg_legacy_view *view, const point_t center)
{
    return rt_view_context_center_vec_set_bsg(qg_legacy_view_to_context(view),
	    center);
}

int
qg_legacy_view_aet_set(qg_legacy_view *view, const vect_t aet)
{
    return rt_view_context_aet_set_bsg(qg_legacy_view_to_context(view), aet);
}

int
qg_legacy_view_update(qg_legacy_view *view)
{
    return rt_view_context_update_bsg(qg_legacy_view_to_context(view));
}

int
qg_legacy_view_screen_to_view(qg_legacy_view *view, fastf_t *vx,
	fastf_t *vy, fastf_t sx, fastf_t sy)
{
    return rt_view_context_screen_to_view_from_bsg(vx, vy,
	    qg_legacy_view_to_context(view), sx, sy);
}

int
qg_legacy_view_screen_point_get(qg_legacy_view *view, point_t point,
	fastf_t sx, fastf_t sy)
{
    return rt_view_context_screen_point_from_bsg(point,
	    qg_legacy_view_to_context(view), sx, sy);
}

int
qg_legacy_view_current_point_get(const qg_legacy_view *view, point_t point)
{
    return rt_view_context_current_point_from_bsg(point,
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_current_point_set(qg_legacy_view *view, const point_t point)
{
    return rt_view_context_current_point_set_bsg(qg_legacy_view_to_context(view),
	    point);
}

int
qg_legacy_view_mouse_state_set(qg_legacy_view *view, int x, int y)
{
    return rt_view_context_mouse_state_set_bsg(qg_legacy_view_to_context(view),
	    x, y);
}

int
qg_legacy_view_view2model_get(const qg_legacy_view *view, mat_t view2model)
{
    rt_view_context_view2model_from_bsg(view2model,
	    qg_legacy_view_to_context(view));
    return view && view2model;
}

int
qg_legacy_view_view2model_set(qg_legacy_view *view, const mat_t view2model)
{
    return rt_view_context_view2model_set_bsg(qg_legacy_view_to_context(view),
	    view2model);
}

int
qg_legacy_view_model2view_get(const qg_legacy_view *view, mat_t model2view)
{
    rt_view_context_model2view_from_bsg(model2view,
	    qg_legacy_view_to_context(view));
    return view && model2view;
}

int
qg_legacy_view_ray_from_screen(qg_legacy_view *view, int sx, int sy,
	point_t origin, vect_t direction)
{
    void *view_ctx = qg_legacy_view_to_context(view);
    if (!view_ctx || !origin || !direction)
	return 0;

    fastf_t vx = -FLT_MAX;
    fastf_t vy = -FLT_MAX;
    if (!rt_view_context_screen_to_view_from_bsg(&vx, &vy, view_ctx, sx, sy))
	return 0;

    point_t vpnt;
    point_t mpnt;
    VSET(vpnt, vx, vy, 0);

    mat_t view2model;
    rt_view_context_view2model_from_bsg(view2model, view_ctx);
    MAT4X3PNT(mpnt, view2model, vpnt);

    mat_t view_rotation;
    rt_view_context_rotation_from_bsg(view_rotation, view_ctx);
    VMOVEN(direction, view_rotation + 8, 3);
    VUNITIZE(direction);
    VSCALE(direction, direction, rt_view_context_radius_from_bsg(view_ctx));
    VADD2(origin, mpnt, direction);
    VUNITIZE(direction);
    VSCALE(direction, direction, -1);
    return 1;
}

int
qg_legacy_view_interactive_rect_state_get(qg_legacy_view *view,
	struct rt_view_interactive_rect_state *state)
{
    return rt_view_context_interactive_rect_state_from_bsg(state,
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_interactive_rect_state_set(qg_legacy_view *view,
	const struct rt_view_interactive_rect_state *state)
{
    return rt_view_context_interactive_rect_state_set_bsg(
	    qg_legacy_view_to_context(view), state);
}

qg_legacy_view_pick_result *
qg_legacy_view_pick_point(qg_legacy_view *view, int x, int y, int first_only)
{
    return qg_pick_result_from_rt(rt_view_context_pick_point_bsg(
	    qg_legacy_view_to_context(view), x, y, first_only));
}

qg_legacy_view_pick_result *
qg_legacy_view_pick_nearest(qg_legacy_view *view, int x, int y)
{
    return qg_pick_result_from_rt(rt_view_context_pick_nearest_bsg(
	    qg_legacy_view_to_context(view), x, y));
}

qg_legacy_view_pick_result *
qg_legacy_view_pick_rect(qg_legacy_view *view, int x0, int y0, int x1, int y1)
{
    return qg_pick_result_from_rt(rt_view_context_pick_rect_bsg(
	    qg_legacy_view_to_context(view), x0, y0, x1, y1));
}

qg_legacy_view_pick_result *
qg_legacy_view_pick_result_create(void)
{
    return qg_pick_result_from_rt(rt_view_pick_result_context_create_bsg());
}

void
qg_legacy_view_pick_result_free(qg_legacy_view_pick_result *result)
{
    rt_view_pick_result_context_free_bsg(qg_pick_result_to_rt(result));
}

size_t
qg_legacy_view_pick_result_count(const qg_legacy_view_pick_result *result)
{
    return rt_view_pick_result_context_count_bsg(qg_pick_result_const_to_rt(
		result));
}

int
qg_legacy_view_pick_result_path(const qg_legacy_view_pick_result *result,
	size_t index, struct bu_vls *path_out)
{
    return rt_view_pick_result_context_path_bsg(qg_pick_result_const_to_rt(result),
	    index, path_out);
}

fastf_t
qg_legacy_view_pick_result_hit_dist(const qg_legacy_view_pick_result *result,
	size_t index)
{
    return rt_view_pick_result_context_hit_dist_bsg(
	    qg_pick_result_const_to_rt(result), index);
}

int
qg_legacy_view_pick_result_append_copy(qg_legacy_view_pick_result *dest,
	const qg_legacy_view_pick_result *src, size_t index, fastf_t hit_dist)
{
    return rt_view_pick_result_context_append_copy_bsg(
	    qg_pick_result_to_rt(dest), qg_pick_result_const_to_rt(src), index,
	    hit_dist);
}

qg_legacy_view_pick_result *
qg_legacy_view_pick_result_filter_first(
	const qg_legacy_view_pick_result *src)
{
    return qg_pick_result_from_rt(rt_view_pick_result_context_filter_first_bsg(
	    qg_pick_result_const_to_rt(src)));
}

int
qg_legacy_view_selection_set_pick_result_ref(qg_legacy_view *view,
	const qg_legacy_view_pick_result *result,
	qg_legacy_view_selection_path_callback_t callback,
	void *data)
{
    return rt_view_context_selection_set_pick_result_context_bsg(
	    qg_legacy_view_to_context(view), qg_pick_result_const_to_rt(result),
	    callback, data);
}

int
qg_legacy_view_selection_clear(qg_legacy_view *view)
{
    return rt_view_context_selection_clear_bsg(qg_legacy_view_to_context(view));
}

int
qg_legacy_view_polygon_ref_is_null(rt_view_polygon_ref ref)
{
    return rt_view_polygon_ref_is_null_bsg(ref);
}

int
qg_legacy_view_polygon_record_get(rt_view_polygon_ref ref,
	struct rt_view_polygon_record *record)
{
    return rt_view_polygon_record_get_bsg(ref, record);
}

rt_view_polygon_ref
qg_legacy_view_polygon_create(qg_legacy_view *view, int type,
	point_t *first_point)
{
    return rt_view_context_polygon_create_bsg(qg_legacy_view_to_context(view), type,
	    first_point);
}

rt_view_polygon_ref
qg_legacy_view_polygon_select(qg_legacy_view *view, point_t *current_point)
{
    return rt_view_context_polygon_select_bsg(qg_legacy_view_to_context(view),
	    current_point);
}

rt_view_polygon_ref
qg_legacy_view_polygon_find(qg_legacy_view *view, const char *name)
{
    return rt_view_context_polygon_find_bsg(qg_legacy_view_to_context(view), name);
}

rt_view_polygon_ref
qg_legacy_view_polygon_dup(qg_legacy_view *view, const char *name,
	const char *new_name)
{
    return rt_view_context_polygon_dup_bsg(qg_legacy_view_to_context(view), name,
	    new_name);
}

void
qg_legacy_view_polygon_visit_records(qg_legacy_view *view,
	qg_legacy_view_polygon_record_callback_t callback,
	void *data)
{
    rt_view_context_polygon_visit_records_bsg(qg_legacy_view_to_context(view),
	    callback, data);
}

size_t
qg_legacy_view_polygon_snap_count(qg_legacy_view *view,
	rt_view_polygon_ref exclude)
{
    return rt_view_context_polygon_snap_count_bsg(qg_legacy_view_to_context(view),
	    exclude);
}

int
qg_legacy_view_polygon_clear_point_selection(qg_legacy_view *view)
{
    return rt_view_context_polygon_clear_point_selection_bsg(
	    qg_legacy_view_to_context(view));
}

int
qg_legacy_view_polygon_update(rt_view_polygon_ref ref, qg_legacy_view *view,
	int utype)
{
    return rt_view_polygon_update_context_bsg(ref, qg_legacy_view_to_context(view),
	    utype);
}

int
qg_legacy_view_polygon_update_screen_point(rt_view_polygon_ref ref,
	qg_legacy_view *view,
	int x,
	int y,
	int utype)
{
    return rt_view_polygon_update_screen_pt_context_bsg(ref,
	    qg_legacy_view_to_context(view), x, y, utype);
}

int
qg_legacy_view_polygon_move(rt_view_polygon_ref ref, point_t *current_point,
	point_t *previous_point)
{
    return rt_view_polygon_move_bsg(ref, current_point, previous_point);
}

int
qg_legacy_view_polygon_set_name(rt_view_polygon_ref ref, const char *name)
{
    return rt_view_polygon_set_name_bsg(ref, name);
}

int
qg_legacy_view_polygon_set_view(rt_view_polygon_ref ref, qg_legacy_view *view)
{
    return rt_view_polygon_set_context_bsg(ref, qg_legacy_view_to_context(view));
}

int
qg_legacy_view_polygon_set_visual(rt_view_polygon_ref ref,
	const struct bu_color *edge_color,
	const struct bu_color *fill_color,
	fastf_t fill_slope_x,
	fastf_t fill_slope_y,
	fastf_t fill_density,
	fastf_t vZ,
	int fill_flag)
{
    return rt_view_polygon_set_visual_bsg(ref, edge_color, fill_color,
	    fill_slope_x, fill_slope_y, fill_density, vZ, fill_flag);
}

int
qg_legacy_view_polygon_set_open(rt_view_polygon_ref ref, int open)
{
    return rt_view_polygon_set_open_bsg(ref, open);
}

int
qg_legacy_view_polygon_close(rt_view_polygon_ref ref)
{
    return rt_view_polygon_close_bsg(ref);
}

int
qg_legacy_view_polygon_clear_selected_point(rt_view_polygon_ref ref)
{
    return rt_view_polygon_clear_selected_point_bsg(ref);
}

int
qg_legacy_view_polygon_remove(rt_view_polygon_ref ref)
{
    return rt_view_polygon_remove_bsg(ref);
}

void *
qg_legacy_view_polygon_user_data(rt_view_polygon_ref ref)
{
    return rt_view_polygon_user_data_bsg(ref);
}

int
qg_legacy_view_polygon_user_data_set(rt_view_polygon_ref ref, void *user_data)
{
    return rt_view_polygon_user_data_set_bsg(ref, user_data);
}

int
qg_legacy_view_polygon_csg(rt_view_polygon_ref target,
	rt_view_polygon_ref stencil,
	bg_clip_t op)
{
    return rt_view_polygon_csg_bsg(target, stencil, op);
}

rt_view_polygon_ref
qg_legacy_view_polygon_import_sketch(const char *name, struct db_i *dbip,
	struct directory *dp,
	qg_legacy_view *view)
{
    return rt_view_polygon_import_sketch_context_bsg(name, dbip, dp,
	    qg_legacy_view_to_context(view));
}

struct directory *
qg_legacy_view_polygon_export_sketch(struct db_i *dbip, const char *name,
	rt_view_polygon_ref ref)
{
    return rt_view_polygon_export_sketch_bsg(dbip, name, ref);
}

int
qg_legacy_view_polygon_snap_exclude_set(qg_legacy_view *view,
	rt_view_polygon_ref ref)
{
    return rt_view_context_polygon_snap_exclude_set_bsg(
	    qg_legacy_view_to_context(view), ref);
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
