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
#include "qtcad/QgLegacyView.h"
#include "qtcad/QgPluginContext.h"

#include "qged_edit_preview_util.h"

static qg_legacy_view_feature_ref
qged_edit_feature_ref_to_legacy(struct qged_edit_feature_ref ref)
{
    qg_legacy_view_feature_ref legacy_ref = { ref.token, ref.revision };
    return legacy_ref;
}


static struct qged_edit_feature_ref
qged_edit_feature_ref_from_legacy(qg_legacy_view_feature_ref ref)
{
    struct qged_edit_feature_ref qged_ref = { ref.token, ref.revision };
    return qged_ref;
}


static qg_legacy_view_feature_family
qged_edit_feature_family_to_legacy(enum qged_edit_feature_family family)
{
    switch (family) {
	case QGED_EDIT_FEATURE_TRANSIENT_PREVIEW:
	    return QG_LEGACY_VIEW_FEATURE_TRANSIENT_PREVIEW;
	case QGED_EDIT_FEATURE_UNKNOWN:
	default:
	    return QG_LEGACY_VIEW_FEATURE_UNKNOWN;
    }
}


static qg_legacy_view_edit_preview_event
qged_edit_preview_event_to_legacy(enum qged_edit_preview_event event)
{
    switch (event) {
	case QGED_EDIT_PREVIEW_BEGIN:
	    return QG_LEGACY_VIEW_EDIT_PREVIEW_BEGIN;
	case QGED_EDIT_PREVIEW_UPDATE:
	    return QG_LEGACY_VIEW_EDIT_PREVIEW_UPDATE;
	case QGED_EDIT_PREVIEW_COMMIT:
	    return QG_LEGACY_VIEW_EDIT_PREVIEW_COMMIT;
	case QGED_EDIT_PREVIEW_CANCEL:
	    return QG_LEGACY_VIEW_EDIT_PREVIEW_CANCEL;
	case QGED_EDIT_PREVIEW_REPLACE_SOURCE:
	    return QG_LEGACY_VIEW_EDIT_PREVIEW_REPLACE_SOURCE;
	case QGED_EDIT_PREVIEW_DISCARD:
	    return QG_LEGACY_VIEW_EDIT_PREVIEW_DISCARD;
	default:
	    return QG_LEGACY_VIEW_EDIT_PREVIEW_UPDATE;
    }
}


static qg_legacy_view *
qged_edit_view(const QgPluginContext *ctx)
{
    return ctx ? ctx->activeView() : nullptr;
}


int
qged_edit_feature_ref_is_null(struct qged_edit_feature_ref ref)
{
    return ref.token ? 0 : 1;
}


static void
qged_edit_preview_callbacks_to_legacy(
	qg_legacy_view_edit_preview_callbacks *dst,
	const struct qged_edit_preview_callbacks *callbacks)
{
    if (!dst || !callbacks)
	return;

    dst->revision_cb = callbacks->revision_cb;
    dst->update_cb = callbacks->update_cb;
    dst->pick_cb = callbacks->pick_cb;
}


int
qged_edit_preview_publish_event(const QgPluginContext *ctx,
				struct qged_edit_feature_ref feature,
				enum qged_edit_preview_event event,
				const char *source_path)
{
    return qg_legacy_view_edit_preview_publish_event(qged_edit_view(ctx),
	    qged_edit_feature_ref_to_legacy(feature),
	    qged_edit_preview_event_to_legacy(event), source_path);
}


struct qged_edit_feature_ref
qged_edit_feature_overlay_ensure(const QgPluginContext *ctx,
				 const char *name,
				 const void *owner,
				 void *preview_ctx,
				 const struct qged_edit_preview_callbacks *callbacks,
				 const char *source_path)
{
    qg_legacy_view_edit_preview_callbacks legacy_callbacks =
	QG_LEGACY_VIEW_EDIT_PREVIEW_CALLBACKS_INIT;
    const qg_legacy_view_edit_preview_callbacks *legacy_callbacks_ptr = NULL;

    if (callbacks) {
	qged_edit_preview_callbacks_to_legacy(&legacy_callbacks, callbacks);
	legacy_callbacks_ptr = &legacy_callbacks;
    }

    return qged_edit_feature_ref_from_legacy(
	    qg_legacy_view_feature_overlay_ensure(qged_edit_view(ctx),
		name, owner, preview_ctx, legacy_callbacks_ptr, source_path));
}


struct qged_edit_feature_ref
qged_edit_feature_label_ensure(const QgPluginContext *ctx,
			       const char *name,
			       const void *owner)
{
    return qged_edit_feature_ref_from_legacy(
	    qg_legacy_view_feature_label_ensure(qged_edit_view(ctx),
		name, owner));
}


int
qged_edit_feature_remove(const QgPluginContext *ctx, const char *name)
{
    return qg_legacy_view_feature_remove(qged_edit_view(ctx), name);
}


void
qged_edit_feature_set_view(struct qged_edit_feature_ref ref, const QgPluginContext *ctx)
{
    qg_legacy_view_feature_set_view(qged_edit_feature_ref_to_legacy(ref),
	    qged_edit_view(ctx));
}


void
qged_edit_feature_set_visible(struct qged_edit_feature_ref ref, int visible)
{
    qg_legacy_view_feature_set_visible(qged_edit_feature_ref_to_legacy(ref),
	    visible);
}


void
qged_edit_feature_set_color(struct qged_edit_feature_ref ref, int r, int g, int b)
{
    qg_legacy_view_feature_set_color(qged_edit_feature_ref_to_legacy(ref),
	    r, g, b);
}


int
qged_edit_feature_touch(struct qged_edit_feature_ref ref)
{
    return qg_legacy_view_feature_touch(qged_edit_feature_ref_to_legacy(ref));
}


int
qged_edit_feature_labels_replace(struct qged_edit_feature_ref ref,
				 const struct rt_point_labels *point_labels,
				 int label_count,
				 const unsigned char color[3])
{
    if (qged_edit_feature_ref_is_null(ref))
	return 0;

    qg_legacy_view_feature_label *labels = NULL;
    if (label_count > 0) {
	labels = (qg_legacy_view_feature_label *)bu_calloc((size_t)label_count,
		sizeof(qg_legacy_view_feature_label), "qged edit labels");
	for (int i = 0; i < label_count; i++) {
	    labels[i].text = point_labels[i].str;
	    VMOVE(labels[i].point, point_labels[i].pt);
	    labels[i].color_valid = 1;
	    labels[i].color[0] = color[0];
	    labels[i].color[1] = color[1];
	    labels[i].color[2] = color[2];
	}
    }

    int ret = qg_legacy_view_feature_labels_replace(
	    qged_edit_feature_ref_to_legacy(ref),
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
    return qg_legacy_view_feature_points_replace(
	    qged_edit_feature_ref_to_legacy(ref),
	    qged_edit_feature_family_to_legacy(family),
	    count ? (const point_t *)lines->points : NULL,
	    count ? lines->cmds : NULL,
	    count);
}


int
qged_edit_feature_clear_geometry(struct qged_edit_feature_ref ref)
{
    return qg_legacy_view_feature_clear_geometry(
	    qged_edit_feature_ref_to_legacy(ref));
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
