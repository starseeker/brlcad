/*             Q G E D _ E D I T _ P R E V I E W _ U T I L . C P P
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

#include "common.h"

#include <algorithm>

#include "bu/malloc.h"
#include "bu/ptbl.h"
#include "ged/draw.h"
#include "ged/defines.h"
#include "ged/view_edit.h"
#include "qtcad/QgPluginContext.h"

#include "qged_edit_preview_util.h"

static ged_view_edit_ref
qged_edit_feature_ref_to_ged(struct qged_edit_feature_ref ref)
{
    ged_view_edit_ref ged_ref = {
	ref.owner, ref.id, ref.generation
    };
    return ged_ref;
}


static struct qged_edit_feature_ref
qged_edit_feature_ref_from_ged(struct ged_view_context *view_ctx,
	ged_view_edit_ref ref)
{
    struct qged_edit_feature_ref qged_ref = {
	view_ctx, ref.owner, ref.id, ref.generation
    };
    return qged_ref;
}


static enum ged_view_edit_geometry_family
qged_edit_feature_family_to_ged(enum qged_edit_feature_family family)
{
    switch (family) {
	case QGED_EDIT_FEATURE_TRANSIENT_PREVIEW:
	    return GED_VIEW_EDIT_GEOMETRY_TRANSIENT_PREVIEW;
	case QGED_EDIT_FEATURE_UNKNOWN:
	default:
	    return GED_VIEW_EDIT_GEOMETRY_UNKNOWN;
    }
}


static enum ged_view_edit_preview_event
qged_edit_preview_event_to_ged(enum qged_edit_preview_event event)
{
    switch (event) {
	case QGED_EDIT_PREVIEW_BEGIN:
	    return GED_VIEW_EDIT_PREVIEW_BEGIN;
	case QGED_EDIT_PREVIEW_UPDATE:
	    return GED_VIEW_EDIT_PREVIEW_UPDATE;
	case QGED_EDIT_PREVIEW_COMMIT:
	    return GED_VIEW_EDIT_PREVIEW_COMMIT;
	case QGED_EDIT_PREVIEW_CANCEL:
	    return GED_VIEW_EDIT_PREVIEW_CANCEL;
	case QGED_EDIT_PREVIEW_REPLACE_SOURCE:
	    return GED_VIEW_EDIT_PREVIEW_REPLACE_SOURCE;
	case QGED_EDIT_PREVIEW_DISCARD:
	    return GED_VIEW_EDIT_PREVIEW_DISCARD;
	default:
	    return GED_VIEW_EDIT_PREVIEW_UPDATE;
    }
}


struct bv_context *
qged_edit_view_context(const QgPluginContext *ctx)
{
    return ctx ? static_cast<struct bv_context *>(ctx->activeViewContext()) : nullptr;
}


struct ged_view_context *
qged_edit_ged_view_context(const QgPluginContext *ctx)
{
    return ged_view_context_from_bv(qged_edit_view_context(ctx));
}


std::vector<struct ged_view_context *>
qged_edit_ged_view_contexts(const QgPluginContext *ctx)
{
    std::vector<struct ged_view_context *> result;
    struct ged *gedp = ctx ? ctx->getGed() : nullptr;
    if (!gedp)
	return result;
    struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
    if (views) {
	for (size_t i = 0; i < BU_PTBL_LEN(views); ++i) {
	    struct ged_view_context *view =
		reinterpret_cast<struct ged_view_context *>(BU_PTBL_GET(views, i));
	    if (view && std::find(result.begin(), result.end(), view) ==
		    result.end())
		result.push_back(view);
	}
    }
    struct ged_view_context *active = ged_view_active_ctx(gedp);
    if (active && std::find(result.begin(), result.end(), active) ==
	    result.end())
	result.push_back(active);
    return result;
}


int
qged_edit_feature_ref_is_null(struct qged_edit_feature_ref ref)
{
    return ged_view_edit_ref_is_null(
	    qged_edit_feature_ref_to_ged(ref));
}


int
qged_edit_preview_publish_event(const QgPluginContext *ctx,
				struct qged_edit_feature_ref feature,
				const char *feature_name,
				enum qged_edit_preview_event event,
				const char *source_path)
{
    struct ged_view_edit_transaction transaction =
	ged_view_edit_transaction_default();
    transaction.event = qged_edit_preview_event_to_ged(event);
    transaction.feature = qged_edit_feature_ref_to_ged(feature);
    transaction.feature_name = feature_name;
    transaction.source_path = source_path;
    return ged_view_edit_transaction_apply(
	qged_edit_ged_view_context(ctx), &transaction, NULL);
}


int
qged_edit_primitive_preview_apply_all(const QgPluginContext *ctx,
	const char *feature_name, enum qged_edit_preview_event event,
	const char *source_path, struct db_i *dbip,
	struct rt_db_internal *internal, const fastf_t *matrix,
	const struct bg_tess_tol *ttol, const struct bn_tol *tol,
	int r, int g, int b)
{
    struct ged *gedp = ctx ? ctx->getGed() : nullptr;
    if (!gedp || !feature_name || !feature_name[0])
	return 0;
    struct ged_view_edit_transaction transaction =
	ged_view_edit_transaction_default();
    transaction.event = qged_edit_preview_event_to_ged(event);
    transaction.feature_name = feature_name;
    transaction.source_path = source_path;
    transaction.dbip = dbip;
    transaction.internal = internal;
    transaction.matrix = matrix;
    transaction.ttol = ttol;
    transaction.tol = tol;
    transaction.color_valid = 1;
    transaction.color[0] = static_cast<unsigned char>(std::clamp(r, 0, 255));
    transaction.color[1] = static_cast<unsigned char>(std::clamp(g, 0, 255));
    transaction.color[2] = static_cast<unsigned char>(std::clamp(b, 0, 255));
    return ged_view_edit_transaction_apply_all(gedp, &transaction);
}


int
qged_edit_feature_remove_all(const QgPluginContext *ctx, const char *name)
{
    if (!name || !name[0])
	return 0;
    int removed = 0;
    for (struct ged_view_context *view : qged_edit_ged_view_contexts(ctx))
	removed += ged_view_feature_remove(view, name);
    return removed;
}


int
qged_edit_feature_labels_replace_all(const QgPluginContext *ctx,
	const char *name, const struct rt_point_labels *point_labels,
	int label_count, const unsigned char color[3])
{
    if (!name || !name[0] || label_count < 0 ||
	(label_count && (!point_labels || !color)))
	return 0;
    struct ged_view_feature_label *labels = nullptr;
    if (label_count) {
	labels = static_cast<struct ged_view_feature_label *>(bu_calloc(
	    static_cast<size_t>(label_count),
	    sizeof(struct ged_view_feature_label), "qged edit labels all views"));
	for (int i = 0; i < label_count; ++i) {
	    labels[i].text = point_labels[i].str;
	    VMOVE(labels[i].point, point_labels[i].pt);
	    labels[i].color_valid = 1;
	    labels[i].color[0] = color[0];
	    labels[i].color[1] = color[1];
	    labels[i].color[2] = color[2];
	}
    }
    int updated = 0;
    for (struct ged_view_context *view : qged_edit_ged_view_contexts(ctx)) {
	ged_view_edit_ref ref = ged_view_edit_label_ensure(view, name);
	if (!ged_view_edit_ref_is_null(ref) &&
	    ged_view_edit_labels_replace(view, ref, labels,
		static_cast<size_t>(label_count)))
	    updated++;
    }
    if (labels)
	bu_free(labels, "qged edit labels all views");
    return updated;
}


int
qged_edit_extrude_preview_apply_all(const QgPluginContext *ctx,
	const char *feature_name, enum qged_edit_preview_event event,
	const char *source_path, const struct rt_extrude_internal *extrude,
	const struct bg_tess_tol *ttol)
{
    if (!feature_name || !feature_name[0] || !extrude)
	return 0;
    int updated = 0;
    for (struct ged_view_context *view : qged_edit_ged_view_contexts(ctx)) {
	ged_view_edit_ref gedRef = ged_view_edit_overlay_ensure(view,
	    feature_name, source_path);
	qged_edit_feature_ref ref = qged_edit_feature_ref_from_ged(view, gedRef);
	if (qged_edit_feature_ref_is_null(ref) ||
	    !qged_edit_feature_replace_extrude_wireframe(ref,
		QGED_EDIT_FEATURE_TRANSIENT_PREVIEW, extrude, ttol))
	    continue;
	(void)ged_view_edit_color_set(view, gedRef, 255, 255, 255);
	(void)ged_view_edit_visible_set(view, gedRef, 1);
	(void)ged_view_edit_preview_publish_event(view, gedRef,
	    qged_edit_preview_event_to_ged(event), source_path);
	updated++;
    }
    return updated;
}


int
qged_edit_revolve_preview_apply_all(const QgPluginContext *ctx,
	const char *feature_name, enum qged_edit_preview_event event,
	const char *source_path, const struct rt_revolve_internal *revolve,
	const struct bg_tess_tol *ttol)
{
    if (!feature_name || !feature_name[0] || !revolve)
	return 0;
    int updated = 0;
    for (struct ged_view_context *view : qged_edit_ged_view_contexts(ctx)) {
	ged_view_edit_ref gedRef = ged_view_edit_overlay_ensure(view,
	    feature_name, source_path);
	qged_edit_feature_ref ref = qged_edit_feature_ref_from_ged(view, gedRef);
	if (qged_edit_feature_ref_is_null(ref) ||
	    !qged_edit_feature_replace_revolve_wireframe(ref,
		QGED_EDIT_FEATURE_TRANSIENT_PREVIEW, revolve, ttol))
	    continue;
	(void)ged_view_edit_color_set(view, gedRef, 255, 255, 255);
	(void)ged_view_edit_visible_set(view, gedRef, 1);
	(void)ged_view_edit_preview_publish_event(view, gedRef,
	    qged_edit_preview_event_to_ged(event), source_path);
	updated++;
    }
    return updated;
}


struct qged_edit_feature_ref
qged_edit_feature_overlay_ensure(const QgPluginContext *ctx,
				 const char *name,
				 const char *source_path)
{
	struct ged_view_context *view_ctx = qged_edit_ged_view_context(ctx);
    return qged_edit_feature_ref_from_ged(view_ctx,
	    ged_view_edit_overlay_ensure(
		view_ctx,
		name, source_path));
}


struct qged_edit_feature_ref
qged_edit_feature_label_ensure(const QgPluginContext *ctx,
			       const char *name)
{
	struct ged_view_context *view_ctx = qged_edit_ged_view_context(ctx);
    return qged_edit_feature_ref_from_ged(view_ctx,
	    ged_view_edit_label_ensure(view_ctx, name));
}


int
qged_edit_feature_remove(const QgPluginContext *ctx, const char *name)
{
    return ged_view_feature_remove(qged_edit_ged_view_context(ctx),
	    name);
}


int
qged_edit_feature_remove_ref(struct qged_edit_feature_ref ref)
{
    return ged_view_edit_remove_ref(ref.view_ctx,
	qged_edit_feature_ref_to_ged(ref));
}


void
qged_edit_feature_set_visible(struct qged_edit_feature_ref ref, int visible)
{
	ged_view_edit_visible_set(ref.view_ctx,
	    qged_edit_feature_ref_to_ged(ref),
	    visible);
}


void
qged_edit_feature_set_color(struct qged_edit_feature_ref ref, int r, int g, int b)
{
	ged_view_edit_color_set(ref.view_ctx,
	    qged_edit_feature_ref_to_ged(ref),
	    r, g, b);
}


int
qged_edit_feature_touch(struct qged_edit_feature_ref ref)
{
    return ged_view_edit_touch(ref.view_ctx,
	qged_edit_feature_ref_to_ged(ref));
}


int
qged_edit_feature_labels_replace(struct qged_edit_feature_ref ref,
				 const struct rt_point_labels *point_labels,
				 int label_count,
				 const unsigned char color[3])
{
    if (qged_edit_feature_ref_is_null(ref))
	return 0;

    struct ged_view_feature_label *labels = NULL;
    if (label_count > 0) {
	labels = (struct ged_view_feature_label *)bu_calloc(
		(size_t)label_count, sizeof(struct ged_view_feature_label),
		"qged edit labels");
	for (int i = 0; i < label_count; i++) {
	    labels[i].text = point_labels[i].str;
	    VMOVE(labels[i].point, point_labels[i].pt);
	    labels[i].color_valid = 1;
	    labels[i].color[0] = color[0];
	    labels[i].color[1] = color[1];
	    labels[i].color[2] = color[2];
	}
    }

    int ret = ged_view_edit_labels_replace(
	    ref.view_ctx,
	    qged_edit_feature_ref_to_ged(ref),
	    labels, label_count > 0 ? (size_t)label_count : 0);
    if (labels)
	bu_free(labels, "qged edit labels");
    return ret;
}


int
qged_edit_preview_lines_replace(struct qged_edit_feature_ref ref,
				enum qged_edit_feature_family family,
				const struct qged_edit_preview_lines *lines)
{
    size_t count = lines ? lines->count : 0;
    return ged_view_edit_points_replace(
	    ref.view_ctx,
	    qged_edit_feature_ref_to_ged(ref),
	    qged_edit_feature_family_to_ged(family),
	    count ? (const point_t *)lines->points : NULL,
	    count ? lines->cmds : NULL,
	    count);
}


int
qged_edit_feature_clear_geometry(struct qged_edit_feature_ref ref)
{
    return ged_view_edit_geometry_clear(
	    ref.view_ctx,
	    qged_edit_feature_ref_to_ged(ref));
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
