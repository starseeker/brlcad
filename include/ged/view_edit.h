/*                      V I E W _ E D I T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged/view_edit.h
 *
 * View-owned retained edit previews.  Every operation that resolves a
 * reference also receives its owning context, making stale and cross-view
 * rejection explicit and eliminating process-global owner lookup.
 */

#ifndef GED_VIEW_EDIT_H
#define GED_VIEW_EDIT_H

#include "ged/view_edit_types.h"
#include "ged/view_feature_types.h"

struct ged;
struct ged_view_context;

__BEGIN_DECLS

/** Return non-zero when @p ref is the null edit reference. */
GED_EXPORT extern int ged_view_edit_ref_is_null(ged_view_edit_ref ref);

/** Ensure a local edit-preview line feature and return its value reference. */
GED_EXPORT extern ged_view_edit_ref ged_view_edit_overlay_ensure(
    struct ged_view_context *view_ctx, const char *name, const void *owner,
    const char *source_path);

/** Ensure a local edit label feature and return its value reference. */
GED_EXPORT extern ged_view_edit_ref ged_view_edit_label_ensure(
    struct ged_view_context *view_ctx, const char *name, const void *owner);

/** Remove exactly the edit feature identified by @p ref. */
GED_EXPORT extern int ged_view_edit_remove_ref(
    struct ged_view_context *view_ctx, ged_view_edit_ref ref);

/** Set visibility after validating @p ref against @p view_ctx. */
GED_EXPORT extern int ged_view_edit_visible_set(
    struct ged_view_context *view_ctx, ged_view_edit_ref ref, int visible);

/** Set an edit feature's RGB color after owner validation. */
GED_EXPORT extern int ged_view_edit_color_set(
    struct ged_view_context *view_ctx, ged_view_edit_ref ref,
    int r, int g, int b);

/** Advance the retained revision of an edit feature. */
GED_EXPORT extern int ged_view_edit_touch(
    struct ged_view_context *view_ctx, ged_view_edit_ref ref);

/** Replace edit label geometry with a deep copy of @p labels. */
GED_EXPORT extern int ged_view_edit_labels_replace(
    struct ged_view_context *view_ctx, ged_view_edit_ref ref,
    const struct ged_view_feature_label *labels, size_t label_count);

/** Replace edit line/point geometry with copied arrays. */
GED_EXPORT extern int ged_view_edit_points_replace(
    struct ged_view_context *view_ctx, ged_view_edit_ref ref,
    enum ged_view_edit_geometry_family family, const point_t *points,
    const int *commands, size_t point_count);

/** Replace an edit preview and its semantic identity/revision metadata. */
GED_EXPORT extern int ged_view_edit_preview_replace(
    struct ged_view_context *view_ctx, ged_view_edit_ref ref,
    const char *source_path, const char *edit_intent_id,
    const char *edit_intent_role, const point_t *points,
    const int *commands, size_t point_count, uint32_t source_revision,
    uint32_t inputs_revision);

/** Remove all geometry while preserving the validated edit feature. */
GED_EXPORT extern int ged_view_edit_geometry_clear(
    struct ged_view_context *view_ctx, ged_view_edit_ref ref);

/** Plot one primitive into the edit preview using optional transform/tolerances. */
GED_EXPORT extern int ged_view_edit_primitive_wireframe_replace(
    struct ged_view_context *view_ctx, ged_view_edit_ref ref,
    struct db_i *dbip, struct rt_db_internal *internal, const mat_t matrix,
    const struct bg_tess_tol *ttol, const struct bn_tol *tol);

/** Publish a live or terminal preview event; terminal events retire it. */
GED_EXPORT extern int ged_view_edit_preview_publish_event(
    struct ged_view_context *view_ctx, ged_view_edit_ref ref,
    enum ged_view_edit_preview_event event, const char *source_path);

/** Apply one typed edit transaction to @p view_ctx. */
GED_EXPORT extern int ged_view_edit_transaction_apply(
    struct ged_view_context *view_ctx,
    const struct ged_view_edit_transaction *transaction,
    ged_view_edit_ref *feature_out);

/** Apply one typed edit transaction to every view owned by @p gedp. */
GED_EXPORT extern int ged_view_edit_transaction_apply_all(
    struct ged *gedp, const struct ged_view_edit_transaction *transaction);

__END_DECLS

#endif /* GED_VIEW_EDIT_H */
