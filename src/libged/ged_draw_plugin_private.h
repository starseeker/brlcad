/*           G E D _ D R A W _ P L U G I N _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged_draw_plugin_private.h
 *
 * Audited, non-installed drawing ABI for in-tree dynamic GED command modules.
 * These declarations are not public libged API.  Additions must also be
 * reviewed in ged_draw_plugin_api.symbols.
 */

#ifndef LIBGED_GED_DRAW_PLUGIN_PRIVATE_H
#define LIBGED_GED_DRAW_PLUGIN_PRIVATE_H

#include "common.h"

#include <stdio.h>

#include "ged.h"
#include "ged/view_feature_types.h"

__BEGIN_DECLS

struct bg_line_layer_builder;
struct db_full_path;
struct ged_uplot_stream;
struct rt_edit;

GED_EXPORT extern int ged_view_feature_gobject_create(
    struct ged *gedp, struct ged_view_context *view_ctx,
    const char *db_path, const char *gobject_name, struct bu_vls *result);

GED_EXPORT extern int _ged_view_feature_batch_publish_uplot(struct ged *gedp,
    FILE *fp, const char *name, double char_size, int mode,
    const char *owner_id, const char *owner_role, const char *remove_prefix,
    const char *result_kind, uint64_t generation);
GED_EXPORT extern int _ged_view_feature_batch_publish_uplot_files(
    struct ged *gedp, const char * const *files, size_t file_count,
    const char *name, double char_size, int mode, const char *owner_id,
    const char *owner_role, const char *remove_prefix,
    const char *result_kind, uint64_t generation);
GED_EXPORT extern struct ged_uplot_stream *_ged_uplot_stream_create(
    double char_size, int mode);
GED_EXPORT extern int _ged_uplot_stream_process(
    struct ged_uplot_stream *stream, FILE *fp, int command);
GED_EXPORT extern int _ged_view_feature_batch_publish_uplot_stream(
    struct ged *gedp, struct ged_uplot_stream *stream, const char *name,
    const char *owner_id, const char *owner_role, const char *remove_prefix,
    const char *result_kind, uint64_t generation);
GED_EXPORT extern void _ged_uplot_stream_free(
    struct ged_uplot_stream *stream);
GED_EXPORT extern int _ged_view_feature_batch_publish_line_layer_builder(
    struct ged *gedp, const char *name,
    const struct bg_line_layer_builder *builder, const char *owner_id,
    const char *owner_role, const char *remove_prefix,
    const char *result_kind, uint64_t generation);
GED_EXPORT extern int _ged_view_feature_batch_publish_line_set(
    struct ged *gedp, const char *name, const point_t *points,
    const int *cmds, size_t point_count,
    const struct ged_view_feature_style *style, const char *owner_id,
    const char *owner_role, const char *remove_prefix,
    const char *result_kind, uint64_t generation);
GED_EXPORT extern int _ged_view_feature_batch_publish_indexed_face_set(
    struct ged *gedp, const char *name, const point_t *points,
    size_t point_count, const vect_t *normals, size_t normal_count,
    const int *indices, size_t index_count,
    const struct ged_view_feature_style *style, const char *owner_id,
    const char *owner_role, const char *remove_prefix,
    const char *result_kind, uint64_t generation);
GED_EXPORT extern int _ged_view_feature_batch_remove_prefix(
    struct ged *gedp, const char *prefix, const char *owner_id,
    const char *owner_role, uint64_t generation);

GED_EXPORT extern int ged_view_data_lines(
    struct ged *gedp, int argc, const char *argv[]);

GED_EXPORT extern struct rt_edit *ged_edit_buf_get(
    struct ged *gedp, const struct db_full_path *dfp);
GED_EXPORT extern void ged_edit_buf_set(
    struct ged *gedp, const struct db_full_path *dfp, struct rt_edit *state);
GED_EXPORT extern int ged_edit_buf_promote(
    struct ged *gedp, const struct db_full_path *dfp);
GED_EXPORT extern void ged_edit_buf_abandon(
    struct ged *gedp, const struct db_full_path *dfp);
GED_EXPORT extern void ged_edit_buf_flush(struct ged *gedp);

__END_DECLS

#endif /* LIBGED_GED_DRAW_PLUGIN_PRIVATE_H */
