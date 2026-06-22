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
#include "bsg/feature.h"
#include "bsg/interaction.h"
#include "qtcad/QgLegacyViewBsg.h"
#include "qtcad/QgPluginContext.h"

#include "qged_edit_preview_util.h"


struct qged_edit_preview_callback_bridge {
    void *preview_ctx;
    struct qged_edit_preview_callbacks callbacks;
};


static bsg_feature_ref
qged_edit_feature_ref_to_bsg(struct qged_edit_feature_ref ref)
{
    bsg_feature_ref bsg_ref = { ref.token, ref.revision };
    return bsg_ref;
}


static struct qged_edit_feature_ref
qged_edit_feature_ref_from_bsg(bsg_feature_ref ref)
{
    struct qged_edit_feature_ref qged_ref = { ref.token, ref.revision };
    return qged_ref;
}


static enum bsg_feature_family
qged_edit_feature_family_to_bsg(enum qged_edit_feature_family family)
{
    switch (family) {
	case QGED_EDIT_FEATURE_TRANSIENT_PREVIEW:
	    return BSG_FEATURE_TRANSIENT_PREVIEW;
	case QGED_EDIT_FEATURE_UNKNOWN:
	default:
	    return BSG_FEATURE_UNKNOWN;
    }
}


static struct bsg_view *
qged_edit_view(const QgPluginContext *ctx)
{
    return ctx ? qg_legacy_view_to_bsg(ctx->activeView()) : nullptr;
}


int
qged_edit_feature_ref_is_null(struct qged_edit_feature_ref ref)
{
    return ref.token ? 0 : 1;
}


static void
qged_edit_preview_bridge_free(void *bridge_ctx)
{
    if (bridge_ctx)
	bu_free(bridge_ctx, "qged edit preview callback bridge");
}


static uint64_t
qged_edit_preview_bridge_revision(void *bridge_ctx)
{
    struct qged_edit_preview_callback_bridge *bridge =
	(struct qged_edit_preview_callback_bridge *)bridge_ctx;
    if (!bridge || !bridge->callbacks.revision_cb)
	return 0;
    return bridge->callbacks.revision_cb(bridge->preview_ctx);
}


static int
qged_edit_preview_bridge_update(void *bridge_ctx, struct bsg_view *UNUSED(v))
{
    struct qged_edit_preview_callback_bridge *bridge =
	(struct qged_edit_preview_callback_bridge *)bridge_ctx;
    if (!bridge || !bridge->callbacks.update_cb)
	return 0;
    return bridge->callbacks.update_cb(bridge->preview_ctx);
}


static int
qged_edit_preview_bridge_pick(void *bridge_ctx, struct bsg_view *UNUSED(v),
			      int x, int y, void *pick_out)
{
    struct qged_edit_preview_callback_bridge *bridge =
	(struct qged_edit_preview_callback_bridge *)bridge_ctx;
    if (!bridge || !bridge->callbacks.pick_cb)
	return 0;
    return bridge->callbacks.pick_cb(bridge->preview_ctx, x, y, pick_out);
}


static void
qged_edit_preview_callbacks_to_bsg(struct bsg_edit_preview_ops *ops,
				   const struct qged_edit_preview_callbacks *callbacks,
				   void *preview_ctx)
{
    if (!ops || !callbacks)
	return;
    if (!callbacks->revision_cb && !callbacks->update_cb && !callbacks->pick_cb)
	return;

    struct qged_edit_preview_callback_bridge *bridge =
	(struct qged_edit_preview_callback_bridge *)bu_calloc(1,
		sizeof(struct qged_edit_preview_callback_bridge),
		"qged edit preview callback bridge");
    bridge->preview_ctx = preview_ctx;
    bridge->callbacks = *callbacks;

    ops->preview_ctx = bridge;
    ops->owns_preview_ctx = 1;
    ops->revision_cb = callbacks->revision_cb ?
	qged_edit_preview_bridge_revision : NULL;
    ops->update_cb = callbacks->update_cb ?
	qged_edit_preview_bridge_update : NULL;
    ops->pick_cb = callbacks->pick_cb ?
	qged_edit_preview_bridge_pick : NULL;
    ops->free_cb = qged_edit_preview_bridge_free;
}


static bsg_edit_preview_op
qged_edit_preview_event_to_bsg(enum qged_edit_preview_event event)
{
    switch (event) {
	case QGED_EDIT_PREVIEW_BEGIN:
	    return BSG_EDIT_PREVIEW_BEGIN;
	case QGED_EDIT_PREVIEW_UPDATE:
	    return BSG_EDIT_PREVIEW_UPDATE;
	case QGED_EDIT_PREVIEW_COMMIT:
	    return BSG_EDIT_PREVIEW_COMMIT;
	case QGED_EDIT_PREVIEW_CANCEL:
	    return BSG_EDIT_PREVIEW_CANCEL;
	case QGED_EDIT_PREVIEW_REPLACE_SOURCE:
	    return BSG_EDIT_PREVIEW_REPLACE_SOURCE;
	case QGED_EDIT_PREVIEW_DISCARD:
	    return BSG_EDIT_PREVIEW_DISCARD;
	default:
	    return BSG_EDIT_PREVIEW_UPDATE;
    }
}


int
qged_edit_preview_publish_event(const QgPluginContext *ctx,
				struct qged_edit_feature_ref feature,
				enum qged_edit_preview_event event,
				const char *source_path)
{
    struct bsg_view *v = qged_edit_view(ctx);
    struct bsg_interaction_record *record =
	bsg_interaction_edit_preview_record(v, qged_edit_feature_ref_to_bsg(feature),
		qged_edit_preview_event_to_bsg(event), source_path);
    if (!record)
	return 0;
    bsg_interaction_record_free(record);
    return 1;
}


struct qged_edit_feature_ref
qged_edit_feature_overlay_ensure(const QgPluginContext *ctx,
				 const char *name,
				 const void *owner,
				 void *preview_ctx,
				 const struct qged_edit_preview_callbacks *callbacks,
				 const char *source_path)
{
    struct bsg_view *v = qged_edit_view(ctx);
    if (!v || !name)
	return QGED_EDIT_FEATURE_REF_NULL;

    bsg_feature_ref ref = bsg_feature_find(v, name);
    if (bsg_feature_ref_is_null(ref)) {
	ref = bsg_feature_create_overlay(v, name, 1 /* local */);
	if (!bsg_feature_ref_is_null(ref)) {
	    bsg_feature_overlay_register_owner(ref, owner,
		    BSG_OVERLAY_ROLE_MODEL,
		    BSG_OVERLAY_CLASS_EDIT_HANDLE,
		    BSG_OVERLAY_LC_PER_TOOL,
		    BSG_OVERLAY_ORDER_POST_TRANSPARENT,
		    NULL, 0);

	    if (callbacks) {
		struct bsg_edit_preview_ops ops = BSG_EDIT_PREVIEW_OPS_INIT;
		qged_edit_preview_callbacks_to_bsg(&ops, callbacks, preview_ctx);
		int attached = bsg_feature_edit_preview_attach(ref, preview_ctx, NULL, &ops);
		if (!attached && ops.owns_preview_ctx && ops.preview_ctx)
		    qged_edit_preview_bridge_free(ops.preview_ctx);
		qged_edit_preview_publish_event(ctx,
			qged_edit_feature_ref_from_bsg(ref),
			QGED_EDIT_PREVIEW_BEGIN, source_path);
	    }
	}
    }

    return qged_edit_feature_ref_from_bsg(ref);
}


struct qged_edit_feature_ref
qged_edit_feature_label_ensure(const QgPluginContext *ctx,
			       const char *name,
			       const void *owner)
{
    struct bsg_view *v = qged_edit_view(ctx);
    if (!v || !name)
	return QGED_EDIT_FEATURE_REF_NULL;

    bsg_feature_ref ref = bsg_feature_find(v, name);
    if (bsg_feature_ref_is_null(ref))
	ref = bsg_feature_create_label(v, name, 1 /* local */);
    if (!bsg_feature_ref_is_null(ref)) {
	bsg_feature_overlay_register_owner(ref, owner,
		BSG_OVERLAY_ROLE_MODEL,
		BSG_OVERLAY_CLASS_EDIT_HANDLE,
		BSG_OVERLAY_LC_PER_TOOL,
		BSG_OVERLAY_ORDER_POST_TRANSPARENT,
		NULL, 1);
    }

    return qged_edit_feature_ref_from_bsg(ref);
}


int
qged_edit_feature_remove(const QgPluginContext *ctx, const char *name)
{
    struct bsg_view *v = qged_edit_view(ctx);
    return (v && name) ? bsg_feature_remove(v, name) : 0;
}


void
qged_edit_feature_set_view(struct qged_edit_feature_ref ref, const QgPluginContext *ctx)
{
    struct bsg_view *v = qged_edit_view(ctx);
    if (!qged_edit_feature_ref_is_null(ref))
	bsg_feature_set_view(qged_edit_feature_ref_to_bsg(ref), v);
}


void
qged_edit_feature_set_visible(struct qged_edit_feature_ref ref, int visible)
{
    if (!qged_edit_feature_ref_is_null(ref))
	bsg_feature_set_visible(qged_edit_feature_ref_to_bsg(ref), visible);
}


void
qged_edit_feature_set_color(struct qged_edit_feature_ref ref, int r, int g, int b)
{
    if (!qged_edit_feature_ref_is_null(ref))
	bsg_feature_set_color(qged_edit_feature_ref_to_bsg(ref), r, g, b);
}


int
qged_edit_feature_touch(struct qged_edit_feature_ref ref)
{
    return qged_edit_feature_ref_is_null(ref) ? 0 :
	bsg_feature_edit_preview_touch(qged_edit_feature_ref_to_bsg(ref));
}


int
qged_edit_feature_labels_replace(struct qged_edit_feature_ref ref,
				 const struct rt_point_labels *point_labels,
				 int label_count,
				 const unsigned char color[3])
{
    if (qged_edit_feature_ref_is_null(ref))
	return 0;

    struct bsg_feature_label_data *labels = NULL;
    if (label_count > 0) {
	labels = (struct bsg_feature_label_data *)bu_calloc((size_t)label_count,
		sizeof(struct bsg_feature_label_data), "qged edit labels");
	for (int i = 0; i < label_count; i++) {
	    labels[i].text = point_labels[i].str;
	    VMOVE(labels[i].point, point_labels[i].pt);
	    labels[i].color_valid = 1;
	    labels[i].color[0] = color[0];
	    labels[i].color[1] = color[1];
	    labels[i].color[2] = color[2];
	    labels[i].anchor = BSG_ANCHOR_AUTO;
	}
    }

    int ret = bsg_feature_labels_replace(qged_edit_feature_ref_to_bsg(ref),
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
    return bsg_feature_points_replace(qged_edit_feature_ref_to_bsg(ref),
	    qged_edit_feature_family_to_bsg(family),
	    count ? (const point_t *)lines->points : NULL,
	    count ? lines->cmds : NULL,
	    count);
}


int
qged_edit_feature_clear_geometry(struct qged_edit_feature_ref ref)
{
    return bsg_feature_points_replace(qged_edit_feature_ref_to_bsg(ref),
	    BSG_FEATURE_UNKNOWN, NULL, NULL, 0);
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
