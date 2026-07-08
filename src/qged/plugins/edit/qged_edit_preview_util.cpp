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
#include "QgLegacyViewContext.h"
#include "qtcad/QgPluginContext.h"

#include "qged_edit_preview_util.h"

static ged_draw_view_feature_ref
qged_edit_feature_ref_to_ged(struct qged_edit_feature_ref ref)
{
    ged_draw_view_feature_ref ged_ref = { ref.token, ref.revision };
    return ged_ref;
}


static struct qged_edit_feature_ref
qged_edit_feature_ref_from_ged(ged_draw_view_feature_ref ref)
{
    struct qged_edit_feature_ref qged_ref = { ref.token, ref.revision };
    return qged_ref;
}


static enum ged_draw_view_feature_family
qged_edit_feature_family_to_ged(enum qged_edit_feature_family family)
{
    switch (family) {
	case QGED_EDIT_FEATURE_TRANSIENT_PREVIEW:
	    return GED_DRAW_VIEW_FEATURE_TRANSIENT_PREVIEW;
	case QGED_EDIT_FEATURE_UNKNOWN:
	default:
	    return GED_DRAW_VIEW_FEATURE_UNKNOWN;
    }
}


static enum ged_draw_view_edit_preview_event
qged_edit_preview_event_to_ged(enum qged_edit_preview_event event)
{
    switch (event) {
	case QGED_EDIT_PREVIEW_BEGIN:
	    return GED_DRAW_VIEW_EDIT_PREVIEW_BEGIN;
	case QGED_EDIT_PREVIEW_UPDATE:
	    return GED_DRAW_VIEW_EDIT_PREVIEW_UPDATE;
	case QGED_EDIT_PREVIEW_COMMIT:
	    return GED_DRAW_VIEW_EDIT_PREVIEW_COMMIT;
	case QGED_EDIT_PREVIEW_CANCEL:
	    return GED_DRAW_VIEW_EDIT_PREVIEW_CANCEL;
	case QGED_EDIT_PREVIEW_REPLACE_SOURCE:
	    return GED_DRAW_VIEW_EDIT_PREVIEW_REPLACE_SOURCE;
	case QGED_EDIT_PREVIEW_DISCARD:
	    return GED_DRAW_VIEW_EDIT_PREVIEW_DISCARD;
	default:
	    return GED_DRAW_VIEW_EDIT_PREVIEW_UPDATE;
    }
}


static qg_legacy_view *
qged_edit_view(const QgPluginContext *ctx)
{
    return ctx ? ctx->activeView() : nullptr;
}


static void *
qged_edit_view_context(const QgPluginContext *ctx)
{
    return qg_legacy_view_to_context(qged_edit_view(ctx));
}


int
qged_edit_feature_ref_is_null(struct qged_edit_feature_ref ref)
{
    return ged_draw_view_feature_ref_is_null(
	    qged_edit_feature_ref_to_ged(ref));
}


static void
qged_edit_preview_callbacks_to_ged(
	struct ged_draw_view_edit_preview_callbacks *dst,
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
    return ged_draw_view_context_edit_preview_publish_event(
	    qged_edit_view_context(ctx),
	    qged_edit_feature_ref_to_ged(feature),
	    qged_edit_preview_event_to_ged(event), source_path);
}


struct qged_edit_feature_ref
qged_edit_feature_overlay_ensure(const QgPluginContext *ctx,
				 const char *name,
				 const void *owner,
				 void *preview_ctx,
				 const struct qged_edit_preview_callbacks *callbacks,
				 const char *source_path)
{
    struct ged_draw_view_edit_preview_callbacks ged_callbacks =
	GED_DRAW_VIEW_EDIT_PREVIEW_CALLBACKS_INIT;
    const struct ged_draw_view_edit_preview_callbacks *ged_callbacks_ptr = NULL;

    if (callbacks) {
	qged_edit_preview_callbacks_to_ged(&ged_callbacks, callbacks);
	ged_callbacks_ptr = &ged_callbacks;
    }

    return qged_edit_feature_ref_from_ged(
	    ged_draw_view_context_feature_overlay_ensure(
		qged_edit_view_context(ctx),
		name, owner, preview_ctx, ged_callbacks_ptr, source_path));
}


struct qged_edit_feature_ref
qged_edit_feature_label_ensure(const QgPluginContext *ctx,
			       const char *name,
			       const void *owner)
{
    return qged_edit_feature_ref_from_ged(
	    ged_draw_view_context_feature_label_ensure(qged_edit_view_context(ctx),
		name, owner));
}


int
qged_edit_feature_remove(const QgPluginContext *ctx, const char *name)
{
    return ged_draw_view_context_feature_remove(qged_edit_view_context(ctx),
	    name);
}


void
qged_edit_feature_set_view(struct qged_edit_feature_ref ref, const QgPluginContext *ctx)
{
    ged_draw_view_feature_set_context(qged_edit_feature_ref_to_ged(ref),
	    qged_edit_view_context(ctx));
}


void
qged_edit_feature_set_visible(struct qged_edit_feature_ref ref, int visible)
{
    ged_draw_view_feature_set_visible(qged_edit_feature_ref_to_ged(ref),
	    visible);
}


void
qged_edit_feature_set_color(struct qged_edit_feature_ref ref, int r, int g, int b)
{
    ged_draw_view_feature_set_color(qged_edit_feature_ref_to_ged(ref),
	    r, g, b);
}


int
qged_edit_feature_touch(struct qged_edit_feature_ref ref)
{
    return ged_draw_view_feature_touch(qged_edit_feature_ref_to_ged(ref));
}


int
qged_edit_feature_labels_replace(struct qged_edit_feature_ref ref,
				 const struct rt_point_labels *point_labels,
				 int label_count,
				 const unsigned char color[3])
{
    if (qged_edit_feature_ref_is_null(ref))
	return 0;

    struct ged_draw_view_feature_label *labels = NULL;
    if (label_count > 0) {
	labels = (struct ged_draw_view_feature_label *)bu_calloc(
		(size_t)label_count, sizeof(struct ged_draw_view_feature_label),
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

    int ret = ged_draw_view_feature_labels_replace(
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
    return ged_draw_view_feature_points_replace(
	    qged_edit_feature_ref_to_ged(ref),
	    qged_edit_feature_family_to_ged(family),
	    count ? (const point_t *)lines->points : NULL,
	    count ? lines->cmds : NULL,
	    count);
}


int
qged_edit_feature_clear_geometry(struct qged_edit_feature_ref ref)
{
    return ged_draw_view_feature_clear_geometry(
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
