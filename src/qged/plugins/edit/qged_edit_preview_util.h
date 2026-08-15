/*             Q G E D _ E D I T _ P R E V I E W _ U T I L . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */

#ifndef QGED_EDIT_PREVIEW_UTIL_H
#define QGED_EDIT_PREVIEW_UTIL_H

#include <vector>

class QgPluginContext;

struct ged_view_context;
struct rt_point_labels;

/** Return the active view, or null when the plugin has no attached view. */
struct ged_view_context *
qged_edit_ged_view_context(const QgPluginContext *ctx);

/** Return each live GED view exactly once, including the active view. */
std::vector<struct ged_view_context *>
qged_edit_ged_view_contexts(const QgPluginContext *ctx);

/** Remove a named plugin-owned presentation feature from every live view. */
int
qged_edit_feature_remove_all(const QgPluginContext *ctx, const char *name);

/** Replace a named set of edit labels in every live view. */
int
qged_edit_feature_labels_replace_all(const QgPluginContext *ctx,
	const char *name, const struct rt_point_labels *labels,
	int label_count, const unsigned char color[3]);

#endif /* QGED_EDIT_PREVIEW_UTIL_H */

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
