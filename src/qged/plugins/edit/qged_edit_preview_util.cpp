/*             Q G E D _ E D I T _ P R E V I E W _ U T I L . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */

#include "common.h"

#include <algorithm>

#include "bu/malloc.h"
#include "bu/ptbl.h"
#include "ged/draw.h"
#include "ged/view_edit.h"
#include "qtcad/QgPluginContext.h"
#include "rt/misc.h"

#include "qged_edit_preview_util.h"


struct ged_view_context *
qged_edit_ged_view_context(const QgPluginContext *ctx)
{
    struct bv_context *view = ctx ?
	static_cast<struct bv_context *>(ctx->activeViewContext()) : nullptr;
    return ged_view_context_from_bv(view);
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
		reinterpret_cast<struct ged_view_context *>(
		    BU_PTBL_GET(views, i));
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
	const char *name, const struct rt_point_labels *pointLabels,
	int labelCount, const unsigned char color[3])
{
    if (!name || !name[0] || labelCount < 0 ||
	(labelCount && (!pointLabels || !color)))
	return 0;

    struct ged_view_feature_label *labels = nullptr;
    if (labelCount) {
	labels = static_cast<struct ged_view_feature_label *>(bu_calloc(
	    static_cast<size_t>(labelCount),
	    sizeof(struct ged_view_feature_label),
	    "qged edit labels all views"));
	for (int i = 0; i < labelCount; ++i) {
	    labels[i].text = pointLabels[i].str;
	    VMOVE(labels[i].point, pointLabels[i].pt);
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
		static_cast<size_t>(labelCount)))
	    updated++;
    }
    if (labels)
	bu_free(labels, "qged edit labels all views");
    return updated;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
