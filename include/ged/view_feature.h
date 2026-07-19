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

#include "common.h"
#include "ged/draw_types.h"

__BEGIN_DECLS

/** Replace a retained screen-space axes overlay for one view. */
GED_EXPORT extern int
ged_annotation_hud_axes_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const struct bv_axes_state *axes,
    const mat_t rotation);

/** Replace retained screen-space line geometry expressed in view coordinates. */
GED_EXPORT extern int
ged_annotation_hud_lines_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    const struct ged_view_feature_style *style);

/** Replace a retained screen-space label batch expressed in view coordinates. */
GED_EXPORT extern int
ged_annotation_hud_labels_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const struct ged_annotation_label *labels,
    size_t label_count,
    const struct ged_view_feature_style *style);

/** Replace retained styled screen-space line layers in view coordinates. */
GED_EXPORT extern int
ged_annotation_hud_line_layers_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const struct ged_annotation_line_layer *layers,
    size_t layer_count,
    const struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_view_feature_remove(struct ged_view_context *view_ctx,
				     const char *name);

/** Remove exactly the feature identified by a retained feature reference. */
GED_EXPORT extern int
ged_view_feature_remove_ref(ged_view_feature_ref ref);

GED_EXPORT extern int
ged_view_feature_remove_prefix(struct ged_view_context *view_ctx,
	const char *prefix);

GED_EXPORT extern int
ged_view_feature_exists(struct ged_view_context *view_ctx,
				     const char *name);

GED_EXPORT extern int
ged_view_feature_visible(struct ged_view_context *view_ctx,
				      const char *name);

GED_EXPORT extern int
ged_view_feature_visible_set(struct ged_view_context *view_ctx,
	const char *name,
	int visible);

GED_EXPORT extern int
ged_view_feature_style_get(
    struct ged_view_context *view_ctx,
    const char *name,
    struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_view_feature_style_apply(
    struct ged_view_context *view_ctx,
    const char *name,
    const struct ged_view_feature_style *style,
    int recursive);

GED_EXPORT extern int
ged_view_feature_realize(struct ged_view_context *view_ctx,
				      const char *name,
				      int recursive);

/**
 * Command-result metadata is an ordered key/value collection.  Diagnostic
 * producers use `result.schema` to identify a versioned schema and
 * `result.severity` values `info`, `warning`, `error`, or `mixed`.  The
 * currently defined command schemas are:
 *
 * * `brlcad.rtcheck.overlap.v1`
 * * `brlcad.nirt.query-ray.v1`
 * * `brlcad.check.overlap.v1`
 * * `brlcad.gqa.{overlap,gap,adjacent-air,exposed-air,volume-sample}.v1`
 *
 * Per-primitive diagnostic records use `result.primitive.kind`,
 * `segment.start_mm`, and `segment.end_mm`.  Ray results additionally use
 * `hit.entry_mm` and `hit.exit_mm`; check overlap records add
 * `overlap.region_a`, `overlap.region_b`, and `overlap.depth_mm`; NIRT
 * partition records may add `nirt.partition.parity`.  Consumers must use the
 * neutral copy APIs below rather than libBObol feature-store records.
 */
GED_EXPORT extern size_t
ged_view_feature_metadata_count(struct ged_view_context *view_ctx,
	const char *name);

GED_EXPORT extern int
ged_view_feature_metadata_copy(struct ged_view_context *view_ctx,
	const char *name,
	size_t index,
	struct bu_vls *key,
	struct bu_vls *value);

GED_EXPORT extern size_t
ged_view_feature_primitive_metadata_count(struct ged_view_context *view_ctx,
	const char *name,
	int primitive);

GED_EXPORT extern int
ged_view_feature_primitive_metadata_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    int primitive,
    size_t index,
    struct bu_vls *key,
    struct bu_vls *value);

GED_EXPORT extern int
ged_view_feature_pick_primitive_resolve(
    struct ged_view_context *view_ctx,
    const char *picked_feature_name,
    int picked_primitive,
    int select,
    int highlight,
    struct bu_vls *feature_name,
    int *feature_primitive);

GED_EXPORT extern int
ged_view_feature_set_selection(
    struct ged_view_context *view_ctx,
    const char *name,
    const int *primitives,
    size_t primitive_count);

GED_EXPORT extern int
ged_view_feature_set_highlights(
    struct ged_view_context *view_ctx,
    const char *name,
    const int *primitives,
    size_t primitive_count);

GED_EXPORT extern size_t
ged_view_feature_selection_count(struct ged_view_context *view_ctx,
	const char *name);

GED_EXPORT extern size_t
ged_view_feature_highlight_count(struct ged_view_context *view_ctx,
	const char *name);

GED_EXPORT extern int
ged_view_feature_selection_at(struct ged_view_context *view_ctx,
	const char *name,
	size_t index,
	int *primitive);

GED_EXPORT extern int
ged_view_feature_highlight_at(struct ged_view_context *view_ctx,
	const char *name,
	size_t index,
	int *primitive);

GED_EXPORT extern int
ged_view_feature_ref_is_null(ged_view_feature_ref ref);

/* Live-edit events refresh preview callbacks/metadata.  Terminal events
 * (`COMMIT`, `CANCEL`, `DISCARD`) retire the transient preview feature so the
 * canonical scene can take over. */
GED_EXPORT extern int
ged_view_feature_edit_preview_publish_event(
    struct ged_view_context *view_ctx,
    ged_view_feature_ref feature,
    enum ged_view_edit_preview_event event,
    const char *source_path);

GED_EXPORT extern int
ged_view_feature_edit_transaction_apply(
    struct ged_view_context *view_ctx,
    const struct ged_view_edit_transaction *transaction,
    ged_view_feature_ref *feature_out);

GED_EXPORT extern int
ged_draw_edit_transaction_apply(
    struct ged *gedp,
    const struct ged_view_edit_transaction *transaction);

GED_EXPORT extern ged_view_feature_ref
ged_view_feature_overlay_ensure(
    struct ged_view_context *view_ctx,
    const char *name,
    const void *owner,
    const char *source_path);

GED_EXPORT extern ged_view_feature_ref
ged_view_feature_label_ensure(struct ged_view_context *view_ctx,
					   const char *name,
					   const void *owner);

GED_EXPORT extern void
ged_view_feature_set_context(ged_view_feature_ref ref,
				  struct ged_view_context *view_ctx);

GED_EXPORT extern void
ged_view_feature_set_visible(ged_view_feature_ref ref,
				  int visible);

GED_EXPORT extern void
ged_view_feature_set_color(ged_view_feature_ref ref,
				int r,
				int g,
				int b);

GED_EXPORT extern int
ged_view_feature_touch(ged_view_feature_ref ref);

GED_EXPORT extern int
ged_view_feature_labels_replace(
    ged_view_feature_ref ref,
    const struct ged_view_feature_label *labels,
    size_t label_count);

GED_EXPORT extern int
ged_view_feature_points_replace(
    ged_view_feature_ref ref,
    enum ged_view_feature_family family,
    const point_t *points,
    const int *cmds,
    size_t point_count);

GED_EXPORT extern int
ged_view_feature_edit_preview_replace(
    ged_view_feature_ref ref,
    const char *source_path,
    const char *edit_intent_id,
    const char *edit_intent_role,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    uint32_t source_revision,
    uint32_t inputs_revision);

GED_EXPORT extern int
ged_view_feature_clear_geometry(ged_view_feature_ref ref);

GED_EXPORT extern int
ged_view_feature_primitive_wireframe_replace(
    ged_view_feature_ref ref,
    struct db_i *dbip,
    struct rt_db_internal *ip,
    const mat_t mat,
    const struct bg_tess_tol *ttol,
    const struct bn_tol *tol);

GED_EXPORT extern int
ged_view_feature_gobject_create(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char *db_path,
    const char *gobject_name,
    struct bu_vls *result);

GED_EXPORT extern int
ged_view_feature_depth(struct ged_view_context *view_ctx,
				    const char *name,
				    int mode,
				    fastf_t *depth);

GED_EXPORT extern int
ged_view_feature_depth_foreach(
    struct ged_view_context *view_ctx,
    int mode,
    ged_view_feature_depth_cb cb,
    void *data);

GED_EXPORT extern void
ged_view_feature_info_get(struct bv_view_info *view_info,
			       const struct ged_view_context *view_ctx);

GED_EXPORT extern int
ged_view_feature_indexed_face_set_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count,
    const struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_annotation_lines_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    const struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_annotation_lines_create_model_annotation(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const point_t point);

GED_EXPORT extern int
ged_annotation_lines_append_point(struct ged_view_context *view_ctx,
	const char *name,
	const point_t point);

GED_EXPORT extern int
ged_annotation_line_create(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const point_t *points,
    size_t point_count,
    const struct ged_view_feature_style *style,
    struct bu_vls *result);

GED_EXPORT extern int
ged_annotation_line_append(
    struct ged_view_context *view_ctx,
    const char *name,
    const point_t *points,
    size_t point_count,
    struct bu_vls *result);

GED_EXPORT extern int
ged_annotation_line_erase(
    struct ged_view_context *view_ctx,
    const char *name,
    size_t start,
    size_t count,
    struct bu_vls *result);

GED_EXPORT extern int
ged_annotation_line_clear(
    struct ged_view_context *view_ctx,
    const char *name,
    struct bu_vls *result);

GED_EXPORT extern int
ged_annotation_line_layer_builder_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const struct bg_line_layer_builder *builder);

GED_EXPORT extern int
ged_annotation_diagnostic_line_layer_builder_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const struct bg_line_layer_builder *builder);

GED_EXPORT extern int
ged_annotation_line_layers_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const struct ged_annotation_line_layer *layers,
    size_t layer_count,
    const struct ged_view_feature_style *style);


GED_EXPORT extern int
ged_annotation_tcl_polygons_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    const struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_annotation_data_polygons_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    int draw,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    const struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_annotation_tcl_labels_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    int draw,
    const struct ged_annotation_label *labels,
    size_t label_count);

GED_EXPORT extern int
ged_annotation_labels_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const struct ged_annotation_label *labels,
    size_t label_count);

GED_EXPORT extern int
ged_annotation_label_set(struct ged_view_context *view_ctx,
				   const char *name,
				   int local,
				   const char *text,
				   const point_t point,
				   const point_t target,
				   int has_target);

GED_EXPORT extern int
ged_annotation_label_create(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const char *text,
    const point_t point,
    const point_t target,
    int has_target,
    struct bu_vls *result);

GED_EXPORT extern size_t
ged_annotation_label_count(struct ged_view_context *view_ctx,
				  const char *name);

GED_EXPORT extern int
ged_annotation_label_copy(struct ged_view_context *view_ctx,
				 const char *name,
				 size_t index,
				 struct bu_vls *text,
				 point_t point,
				 unsigned char rgb[3]);

GED_EXPORT extern int
ged_annotation_label_point_set(struct ged_view_context *view_ctx,
				      const char *name,
				      size_t index,
				      const point_t point);

GED_EXPORT extern int
ged_annotation_line_color_set(struct ged_view_context *view_ctx,
				     const char *name,
				     int r,
				     int g,
				     int b);

GED_EXPORT extern int
ged_annotation_line_width_set(struct ged_view_context *view_ctx,
				     const char *name,
				     int line_width);

GED_EXPORT extern int
ged_annotation_line_style_get(
    struct ged_view_context *view_ctx,
    const char *name,
    struct ged_draw_view_line_style *style);

GED_EXPORT extern int
ged_view_feature_points_copy(struct ged_view_context *view_ctx,
	const char *name,
	point_t **points,
	size_t *point_count);

GED_EXPORT extern int
ged_view_feature_line_command_at(
    struct ged_view_context *view_ctx,
    const char *name,
    size_t index,
    int *out);

GED_EXPORT extern int
ged_annotation_lines_points_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    point_t **points,
    size_t *point_count);

GED_EXPORT extern int
ged_annotation_tcl_lines_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const point_t *points,
    size_t point_count,
    const struct ged_draw_view_line_style *style);

GED_EXPORT extern int
ged_annotation_data_lines_draw_get(
    struct ged_view_context *view_ctx,
    const char *name);

GED_EXPORT extern int
ged_annotation_data_lines_draw_set(
    struct ged_view_context *view_ctx,
    const char *name,
    int draw);

GED_EXPORT extern int
ged_annotation_data_lines_style_get(
    struct ged_view_context *view_ctx,
    const char *name,
    struct ged_draw_view_line_style *style);

GED_EXPORT extern int
ged_annotation_data_lines_color_set(
    struct ged_view_context *view_ctx,
    const char *name,
    int r,
    int g,
    int b);

GED_EXPORT extern int
ged_annotation_data_lines_line_width_set(
    struct ged_view_context *view_ctx,
    const char *name,
    int line_width);

GED_EXPORT extern int
ged_annotation_data_lines_points_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    point_t **points,
    size_t *point_count);

GED_EXPORT extern int
ged_annotation_data_lines_points_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const point_t *points,
    size_t point_count,
    const struct ged_draw_view_line_style *style);

GED_EXPORT extern int
ged_annotation_data_labels_draw_get(
    struct ged_view_context *view_ctx,
    const char *name);

GED_EXPORT extern int
ged_annotation_data_labels_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    int draw,
    const struct ged_annotation_label *labels,
    size_t label_count);

GED_EXPORT extern int
ged_annotation_data_labels_color_get(
    struct ged_view_context *view_ctx,
    const char *name,
    unsigned char rgb[3]);

GED_EXPORT extern size_t
ged_annotation_data_labels_count(
    struct ged_view_context *view_ctx,
    const char *name);

GED_EXPORT extern int
ged_annotation_data_labels_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    size_t index,
    struct bu_vls *text,
    point_t point,
    unsigned char rgb[3]);

GED_EXPORT extern int
ged_annotation_data_labels_point_set(
    struct ged_view_context *view_ctx,
    const char *name,
    size_t index,
    const point_t point);

GED_EXPORT extern int
ged_annotation_data_arrows_draw_get(
    struct ged_view_context *view_ctx,
    const char *name);

GED_EXPORT extern int
ged_annotation_data_arrows_draw_set(
    struct ged_view_context *view_ctx,
    const char *name,
    int draw);

GED_EXPORT extern int
ged_annotation_data_arrows_style_get(
    struct ged_view_context *view_ctx,
    const char *name,
    struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_annotation_data_arrows_color_set(
    struct ged_view_context *view_ctx,
    const char *name,
    int r,
    int g,
    int b);

GED_EXPORT extern int
ged_annotation_data_arrows_line_width_set(
    struct ged_view_context *view_ctx,
    const char *name,
    int line_width);

GED_EXPORT extern int
ged_annotation_data_arrows_points_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    point_t **points,
    size_t *point_count);

GED_EXPORT extern int
ged_annotation_data_arrows_points_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    point_t *points,
    size_t point_count,
    const struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_annotation_data_arrows_tip_get(
    struct ged_view_context *view_ctx,
    const char *name,
    fastf_t *tip_length,
    fastf_t *tip_width);

GED_EXPORT extern int
ged_annotation_data_arrows_tip_set(
    struct ged_view_context *view_ctx,
    const char *name,
    fastf_t tip_length,
    fastf_t tip_width);

GED_EXPORT extern int
ged_annotation_arrow_tip_get(struct ged_view_context *view_ctx,
				    const char *name,
				    fastf_t *tip_length,
				    fastf_t *tip_width);

GED_EXPORT extern int
ged_annotation_arrow_tip_set(struct ged_view_context *view_ctx,
				    const char *name,
				    fastf_t tip_length,
				    fastf_t tip_width);

GED_EXPORT extern int
ged_annotation_tcl_arrows_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    point_t *points,
    size_t point_count,
    const struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_view_feature_axes_centers_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    point_t **centers,
    size_t *center_count);

GED_EXPORT extern int
ged_annotation_tcl_axes_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    point_t *centers,
    size_t center_count,
    fastf_t half_axes_size,
    const struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_annotation_data_axes_draw_get(
    struct ged_view_context *view_ctx,
    const char *name);

GED_EXPORT extern int
ged_annotation_data_axes_draw_set(
    struct ged_view_context *view_ctx,
    const char *name,
    int draw);

GED_EXPORT extern int
ged_annotation_data_axes_style_get(
    struct ged_view_context *view_ctx,
    const char *name,
    struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_annotation_data_axes_color_set(
    struct ged_view_context *view_ctx,
    const char *name,
    int r,
    int g,
    int b);

GED_EXPORT extern int
ged_annotation_data_axes_line_width_set(
    struct ged_view_context *view_ctx,
    const char *name,
    int line_width);

GED_EXPORT extern int
ged_annotation_data_axes_half_size_get(
    struct ged_view_context *view_ctx,
    const char *name,
    fastf_t *half_axes_size);

GED_EXPORT extern int
ged_annotation_data_axes_size_get(
    struct ged_view_context *view_ctx,
    const char *name,
    fastf_t display_scale,
    fastf_t *size);

GED_EXPORT extern int
ged_annotation_data_axes_centers_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    point_t **centers,
    size_t *center_count);

GED_EXPORT extern int
ged_annotation_data_axes_centers_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    point_t *centers,
    size_t center_count,
    fastf_t half_axes_size,
    const struct ged_view_feature_style *style);

GED_EXPORT extern int
ged_annotation_axes_set(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const struct ged_annotation_axes *state);

GED_EXPORT extern int
ged_annotation_axes_get(
    struct ged_view_context *view_ctx,
    const char *name,
    struct ged_annotation_axes *state);

GED_EXPORT extern int
ged_annotation_axes_update(
    struct ged_view_context *view_ctx,
    const char *name,
    const struct ged_annotation_axes *state);

GED_EXPORT extern int
ged_annotation_axes_create(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const struct ged_annotation_axes *state,
    struct bu_vls *result);

GED_EXPORT extern int
ged_annotation_axes_state_get(
    struct ged_view_context *view_ctx,
    const char *name,
    struct ged_annotation_axes *state,
    struct bu_vls *result);

GED_EXPORT extern int
ged_annotation_axes_state_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const struct ged_annotation_axes *state,
    struct bu_vls *result);


GED_EXPORT extern int
ged_view_feature_get_summary(
    struct ged_view_context *view_ctx,
    const char *name,
    struct ged_view_feature_summary *summary);

__END_DECLS

#endif /* GED_VIEW_FEATURE_H */
