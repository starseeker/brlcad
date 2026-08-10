/* V I E W _ F E A T U R E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged/view_feature.h */

#ifndef GED_VIEW_FEATURE_H
#define GED_VIEW_FEATURE_H

#include "ged/view_feature_types.h"

__BEGIN_DECLS

/** Remove a named retained feature from @p view_ctx. */
GED_EXPORT extern int
ged_view_feature_remove(struct ged_view_context *view_ctx, const char *name);

/** Remove all retained features whose names start with @p prefix. */
GED_EXPORT extern int
ged_view_feature_remove_prefix(struct ged_view_context *view_ctx,
	const char *prefix);

/** Return non-zero when @p view_ctx contains a feature named @p name. */
GED_EXPORT extern int
ged_view_feature_exists(struct ged_view_context *view_ctx, const char *name);

/** Return non-zero when the named feature is currently visible. */
GED_EXPORT extern int
ged_view_feature_visible(struct ged_view_context *view_ctx, const char *name);

/** Set the visibility of a named retained feature. */
GED_EXPORT extern int
ged_view_feature_visible_set(struct ged_view_context *view_ctx,
	const char *name, int visible);

/** Copy the effective style of a named retained feature. */
GED_EXPORT extern int
ged_view_feature_style_get(struct ged_view_context *view_ctx,
	const char *name, struct ged_view_feature_style *style);

/** Apply @p style to a named feature and optionally to its children. */
GED_EXPORT extern int
ged_view_feature_style_apply(struct ged_view_context *view_ctx,
	const char *name, const struct ged_view_feature_style *style,
	int recursive);

/** Ensure the named feature, and optionally its children, is realized. */
GED_EXPORT extern int
ged_view_feature_realize(struct ged_view_context *view_ctx,
	const char *name, int recursive);

/** Return the number of metadata entries associated with a feature. */
GED_EXPORT extern size_t
ged_view_feature_metadata_count(struct ged_view_context *view_ctx,
	const char *name);

/** Copy one feature metadata key/value pair into caller-owned strings. */
GED_EXPORT extern int
ged_view_feature_metadata_copy(struct ged_view_context *view_ctx,
	const char *name, size_t index, struct bu_vls *key,
	struct bu_vls *value);

/** Return the metadata count for one feature primitive. */
GED_EXPORT extern size_t
ged_view_feature_primitive_metadata_count(struct ged_view_context *view_ctx,
	const char *name, int primitive);

/** Copy one primitive metadata key/value pair into caller-owned strings. */
GED_EXPORT extern int
ged_view_feature_primitive_metadata_copy(struct ged_view_context *view_ctx,
	const char *name, int primitive, size_t index, struct bu_vls *key,
	struct bu_vls *value);

/** Resolve a picked child primitive to its owning retained feature. */
GED_EXPORT extern int
ged_view_feature_pick_primitive_resolve(struct ged_view_context *view_ctx,
	const char *picked_feature_name, int picked_primitive, int select,
	int highlight, struct bu_vls *feature_name, int *feature_primitive);

/** Replace the selected primitive set of a retained feature. */
GED_EXPORT extern int
ged_view_feature_set_selection(struct ged_view_context *view_ctx,
	const char *name, const int *primitives, size_t primitive_count);

/** Replace the highlighted primitive set of a retained feature. */
GED_EXPORT extern int
ged_view_feature_set_highlights(struct ged_view_context *view_ctx,
	const char *name, const int *primitives, size_t primitive_count);

/** Return the number of selected primitives in a retained feature. */
GED_EXPORT extern size_t
ged_view_feature_selection_count(struct ged_view_context *view_ctx,
	const char *name);

/** Return the number of highlighted primitives in a retained feature. */
GED_EXPORT extern size_t
ged_view_feature_highlight_count(struct ged_view_context *view_ctx,
	const char *name);

/** Copy the selected primitive index at @p index. */
GED_EXPORT extern int
ged_view_feature_selection_at(struct ged_view_context *view_ctx,
	const char *name, size_t index, int *primitive);

/** Copy the highlighted primitive index at @p index. */
GED_EXPORT extern int
ged_view_feature_highlight_at(struct ged_view_context *view_ctx,
	const char *name, size_t index, int *primitive);

/** Copy a renderer-neutral summary of a named retained feature. */
GED_EXPORT extern int
ged_view_feature_get_summary(struct ged_view_context *view_ctx,
	const char *name, struct ged_view_feature_summary *summary);

__END_DECLS

#endif /* GED_VIEW_FEATURE_H */
