/* R E S U L T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged/result.h */

#ifndef GED_RESULT_H
#define GED_RESULT_H

#include "common.h"
#include "ged/draw_types.h"

__BEGIN_DECLS

GED_EXPORT extern struct ged_result_scene *
ged_result_begin(
    struct ged_view_context *view_ctx,
    const struct ged_result_desc *desc);

GED_EXPORT extern size_t
ged_result_features_remove_prefix(
    struct ged_result_scene *scene,
    const char *prefix);

GED_EXPORT extern int
ged_result_line_layer_builder_replace(
    struct ged_result_scene *scene,
    const char *name,
    const struct bg_line_layer_builder *builder,
    const struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_result_line_layers_replace(
    struct ged_result_scene *scene,
    const char *name,
    const struct ged_annotation_line_layer *layers,
    size_t layer_count,
    const struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_result_line_set_replace(
    struct ged_result_scene *scene,
    const char *name,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    const struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_result_point_set_replace(
    struct ged_result_scene *scene,
    const char *name,
    const point_t *points,
    size_t point_count,
    const struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_result_indexed_face_set_replace(
    struct ged_result_scene *scene,
    const char *name,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count,
    const struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_result_hud_label_replace(
    struct ged_result_scene *scene,
    const char *name,
    const struct ged_diagnostic_hud_label *label);

GED_EXPORT extern int
ged_result_custom_node_replace(
    struct ged_result_scene *scene,
    const char *name,
    ged_result_custom_node_cb node_cb,
    void *node_cb_data,
    const struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_result_feature_metadata_replace(
    struct ged_result_scene *scene,
    const char *name,
    const struct ged_result_metadata *metadata,
    size_t metadata_count);

GED_EXPORT extern int
ged_result_feature_primitive_metadata_replace(
    struct ged_result_scene *scene,
    const char *name,
    int primitive,
    const struct ged_result_metadata *metadata,
    size_t metadata_count);

GED_EXPORT extern int
ged_result_commit(struct ged_result_scene *scene);

GED_EXPORT extern void
ged_result_abort(struct ged_result_scene *scene);

__END_DECLS

#endif /* GED_RESULT_H */

