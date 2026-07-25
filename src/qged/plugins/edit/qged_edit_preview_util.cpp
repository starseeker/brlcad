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

#include "bu/malloc.h"
#include "ged/draw.h"
#include "qtcad/QgPluginContext.h"

#include "qged_edit_preview_util.h"

static ged_view_feature_ref
qged_edit_feature_ref_to_ged(struct qged_edit_feature_ref ref)
{
    ged_view_feature_ref ged_ref = {
	ref.owner, ref.id, ref.generation
    };
    return ged_ref;
}


static struct qged_edit_feature_ref
qged_edit_feature_ref_from_ged(ged_view_feature_ref ref)
{
    struct qged_edit_feature_ref qged_ref = {
	ref.owner, ref.id, ref.generation
    };
    return qged_ref;
}


static enum ged_view_feature_family
qged_edit_feature_family_to_ged(enum qged_edit_feature_family family)
{
    switch (family) {
	case QGED_EDIT_FEATURE_TRANSIENT_PREVIEW:
	    return GED_VIEW_FEATURE_TRANSIENT_PREVIEW;
	case QGED_EDIT_FEATURE_UNKNOWN:
	default:
	    return GED_VIEW_FEATURE_UNKNOWN;
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


int
qged_edit_feature_ref_is_null(struct qged_edit_feature_ref ref)
{
    return ged_view_feature_ref_is_null(
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
	GED_VIEW_EDIT_TRANSACTION_INIT;
    transaction.event = qged_edit_preview_event_to_ged(event);
    transaction.feature = qged_edit_feature_ref_to_ged(feature);
    transaction.feature_name = feature_name;
    transaction.source_path = source_path;
    return ged_view_feature_edit_transaction_apply(
	qged_edit_ged_view_context(ctx), &transaction, NULL);
}


struct qged_edit_feature_ref
qged_edit_feature_overlay_ensure(const QgPluginContext *ctx,
				 const char *name,
				 const void *owner,
				 const char *source_path)
{
    return qged_edit_feature_ref_from_ged(
	    ged_view_feature_overlay_ensure(
		qged_edit_ged_view_context(ctx),
		name, owner, source_path));
}


struct qged_edit_feature_ref
qged_edit_feature_label_ensure(const QgPluginContext *ctx,
			       const char *name,
			       const void *owner)
{
    return qged_edit_feature_ref_from_ged(
	    ged_view_feature_label_ensure(qged_edit_ged_view_context(ctx),
		name, owner));
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
    return ged_view_feature_remove_ref(qged_edit_feature_ref_to_ged(ref));
}


void
qged_edit_feature_set_view(struct qged_edit_feature_ref ref, const QgPluginContext *ctx)
{
    ged_view_feature_set_context(qged_edit_feature_ref_to_ged(ref),
	    qged_edit_ged_view_context(ctx));
}


void
qged_edit_feature_set_visible(struct qged_edit_feature_ref ref, int visible)
{
    ged_view_feature_set_visible(qged_edit_feature_ref_to_ged(ref),
	    visible);
}


void
qged_edit_feature_set_color(struct qged_edit_feature_ref ref, int r, int g, int b)
{
    ged_view_feature_set_color(qged_edit_feature_ref_to_ged(ref),
	    r, g, b);
}


int
qged_edit_feature_touch(struct qged_edit_feature_ref ref)
{
    return ged_view_feature_touch(qged_edit_feature_ref_to_ged(ref));
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

    int ret = ged_view_feature_labels_replace(
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
    return ged_view_feature_points_replace(
	    qged_edit_feature_ref_to_ged(ref),
	    qged_edit_feature_family_to_ged(family),
	    count ? (const point_t *)lines->points : NULL,
	    count ? lines->cmds : NULL,
	    count);
}


int
qged_edit_feature_clear_geometry(struct qged_edit_feature_ref ref)
{
    return ged_view_feature_clear_geometry(
	    qged_edit_feature_ref_to_ged(ref));
}

int
qged_edit_feature_replace_primitive_wireframe(
	struct qged_edit_feature_ref ref,
	struct db_i *dbip,
	struct rt_db_internal *ip,
	const mat_t mat,
	const struct bg_tess_tol *ttol,
	const struct bn_tol *tol)
{
    return ged_view_feature_primitive_wireframe_replace(
	    qged_edit_feature_ref_to_ged(ref), dbip, ip, mat, ttol, tol);
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
