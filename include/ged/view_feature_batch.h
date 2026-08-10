/*        V I E W _ F E A T U R E _ B A T C H . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged/view_feature_batch.h
 *
 * Staged retained view-feature publication.  A batch owns copies of related
 * feature geometry and metadata until one generation-checked commit.  The
 * commit applies the operations in staged order before returning; callbacks
 * report each applied feature.
 */

#ifndef GED_VIEW_FEATURE_BATCH_H
#define GED_VIEW_FEATURE_BATCH_H

#include "ged/view_feature_types.h"

__BEGIN_DECLS

/** Begin an owner-scoped staged feature batch for @p view_ctx. */
GED_EXPORT extern struct ged_view_feature_batch *
ged_view_feature_batch_begin(
    struct ged_view_context *view_ctx,
    const struct ged_view_feature_batch_desc *desc);

/** Stage removal of batch-owned features whose names start with @p prefix. */
GED_EXPORT extern size_t
ged_view_feature_batch_remove_prefix(
    struct ged_view_feature_batch *batch,
    const char *prefix);

/** Stage line layers accumulated by a libbg line-layer builder. */
GED_EXPORT extern int
ged_view_feature_batch_line_layer_builder_replace(
    struct ged_view_feature_batch *batch,
    const char *name,
    const struct bg_line_layer_builder *builder,
    const struct ged_view_feature_style *style);

/** Stage a named collection of independently styled line layers. */
GED_EXPORT extern int
ged_view_feature_batch_line_layers_replace(
    struct ged_view_feature_batch *batch,
    const char *name,
    const struct ged_annotation_line_layer *layers,
    size_t layer_count,
    const struct ged_view_feature_style *style);

/** Stage a named line set.  Point and command arrays are copied. */
GED_EXPORT extern int
ged_view_feature_batch_line_set_replace(
    struct ged_view_feature_batch *batch,
    const char *name,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    const struct ged_view_feature_style *style);

/** Stage a named point set.  The point array is copied. */
GED_EXPORT extern int
ged_view_feature_batch_point_set_replace(
    struct ged_view_feature_batch *batch,
    const char *name,
    const point_t *points,
    size_t point_count,
    const struct ged_view_feature_style *style);

/** Stage a named indexed triangle surface; all arrays are copied. */
GED_EXPORT extern int
ged_view_feature_batch_indexed_face_set_replace(
    struct ged_view_feature_batch *batch,
    const char *name,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count,
    const struct ged_view_feature_style *style);

/** Stage one named HUD label. */
GED_EXPORT extern int
ged_view_feature_batch_hud_label_replace(
    struct ged_view_feature_batch *batch,
    const char *name,
    const struct ged_diagnostic_hud_label *label);

/** Replace metadata associated with a staged named feature. */
GED_EXPORT extern int
ged_view_feature_batch_metadata_replace(
    struct ged_view_feature_batch *batch,
    const char *name,
    const struct ged_view_feature_metadata *metadata,
    size_t metadata_count);

/** Replace metadata associated with one primitive of a staged feature. */
GED_EXPORT extern int
ged_view_feature_batch_primitive_metadata_replace(
    struct ged_view_feature_batch *batch,
    const char *name,
    int primitive,
    const struct ged_view_feature_metadata *metadata,
    size_t metadata_count);

/** Publish the staged changes in order and destroy @p batch. */
GED_EXPORT extern int
ged_view_feature_batch_commit(struct ged_view_feature_batch *batch);

/** Discard staged changes and destroy @p batch. */
GED_EXPORT extern void
ged_view_feature_batch_abort(struct ged_view_feature_batch *batch);

__END_DECLS

#endif /* GED_VIEW_FEATURE_BATCH_H */
