/*          B S G _ G E D _ D R A W _ V I E W _ P R I V A T E . H
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
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file bsg_ged_draw_view_private.h
 *
 * Narrow BSG-backed view adapter surface for libged command code.
 *
 * The full bsg_ged_draw_private.h bridge exposes retained draw-scene internals.
 * This header is limited to view-scoped adapters still needed while GED commands
 * migrate away from direct BSG view state.
 */

#ifndef LIBGED_BSG_GED_DRAW_VIEW_PRIVATE_H
#define LIBGED_BSG_GED_DRAW_VIEW_PRIVATE_H

#include "common.h"

#include "bg/polygon_types.h"
#include "bu/color.h"
#include "bsg/scene_builder.h"
#include "bsg/scene_object.h"
#include "ged/bsg_ged_draw.h"
#include "ged/defines.h"
#include "rt/view.h"
#include "vmath.h"

__BEGIN_DECLS

struct bsg_view;
struct bsg_view_set;
struct bg_line_layer_builder;
struct bn_tol;
struct bu_ptbl;
struct bu_vls;
struct db_i;
struct directory;
struct rt_mesh_lod;

typedef int (*ged_draw_view_selection_path_cb)(struct bsg_view *view,
					       const char *path,
					       void *data);

enum ged_draw_view_selection_kind {
    GED_DRAW_VIEW_SELECTION_ALL = -1,
    GED_DRAW_VIEW_SELECTION_SELECTED_PATH = 0,
    GED_DRAW_VIEW_SELECTION_HIGHLIGHTED_REF = 1
};

GED_EXPORT extern void ged_draw_view_info_from_bsg(struct rt_view_info *view_info,
						   const struct bsg_view *view);
GED_EXPORT extern void ged_draw_view_context_info_from_bsg(struct rt_view_info *view_info,
							   const void *view_ctx);
GED_EXPORT extern const char *ged_draw_view_context_name_from_bsg(const void *view_ctx);
GED_EXPORT extern fastf_t ged_draw_view_perspective_from_bsg(const struct bsg_view *view);
GED_EXPORT extern fastf_t ged_draw_view_context_perspective_from_bsg(
	const void *view_ctx);
GED_EXPORT extern fastf_t ged_draw_view_scale_from_bsg(const struct bsg_view *view);
GED_EXPORT extern fastf_t ged_draw_view_context_scale_from_bsg(
	const void *view_ctx);
GED_EXPORT extern int ged_draw_view_context_obb_from_bsg(
	const void *view_ctx,
	point_t center,
	vect_t extent1,
	vect_t extent2,
	vect_t extent3);
GED_EXPORT extern int ged_draw_view_lod_policy_from_bsg(ged_draw_view_lod_policy *policy,
							const struct bsg_view *view);
GED_EXPORT extern int ged_draw_view_context_lod_policy_from_bsg(
	ged_draw_view_lod_policy *policy,
	const void *view_ctx);
GED_EXPORT extern int ged_draw_view_lod_policy_apply_bsg(struct bsg_view *view,
							 const ged_draw_view_lod_policy *policy);
GED_EXPORT extern int ged_draw_view_context_lod_policy_apply_bsg(
	void *view_ctx,
	const ged_draw_view_lod_policy *policy);
GED_EXPORT extern int ged_draw_view_lod_policy_apply_bsg_bot_threshold(
	struct bsg_view *view,
	const ged_draw_view_lod_policy *policy,
	size_t bot_threshold);
GED_EXPORT extern int ged_draw_view_context_lod_policy_apply_bsg_bot_threshold(
	void *view_ctx,
	const ged_draw_view_lod_policy *policy,
	size_t bot_threshold);
GED_EXPORT extern int ged_draw_view_autoview_default_bsg(struct bsg_view *view,
							 int all_view_objs);
GED_EXPORT extern int ged_draw_view_hud_sync(struct bsg_view *view);
GED_EXPORT extern int ged_draw_view_context_hud_sync(void *view_ctx);
GED_EXPORT extern int ged_draw_view_selection_available(struct bsg_view *view);
GED_EXPORT extern size_t ged_draw_view_selection_count(struct bsg_view *view);
GED_EXPORT extern size_t ged_draw_view_context_selection_count(void *view_ctx);
GED_EXPORT extern int ged_draw_view_selection_path_foreach(struct bsg_view *view,
							   ged_draw_view_selection_path_cb cb,
							   void *data);
GED_EXPORT extern int ged_draw_view_context_selection_path_foreach(
	void *view_ctx,
	ged_draw_view_context_selection_path_cb cb,
	void *data);
GED_EXPORT extern int ged_draw_view_selection_clear(struct bsg_view *view);
GED_EXPORT extern int ged_draw_view_context_selection_clear(void *view_ctx);
GED_EXPORT extern int ged_draw_view_selection_contains_path(
	struct bsg_view *view,
	enum ged_draw_view_selection_kind kind,
	const char *path);
GED_EXPORT extern int ged_draw_view_selection_add_path(
	struct bsg_view *view,
	enum ged_draw_view_selection_kind kind,
	const char *path);
GED_EXPORT extern int ged_draw_view_selection_set_path(
	struct bsg_view *view,
	enum ged_draw_view_selection_kind kind,
	const char *path);
GED_EXPORT extern int ged_draw_view_snap_first_candidate(struct bsg_view *view,
							 const point_t sample,
							 enum ged_draw_view_snap_kind kind,
							 point_t candidate);
GED_EXPORT extern int ged_draw_view_context_snap_first_candidate(void *view_ctx,
								 const point_t sample,
								 enum ged_draw_view_snap_kind kind,
								 point_t candidate);
GED_EXPORT extern struct bu_ptbl *ged_draw_view_set_views_bsg(struct bsg_view_set *view_set);
GED_EXPORT extern void *ged_draw_view_set_recycle_pool_bsg(struct bsg_view_set *view_set);
GED_EXPORT extern int ged_draw_view_is_independent_bsg(const struct bsg_view *view);
GED_EXPORT extern bsg_scene_ref ged_draw_view_independent_scope_ref_bsg(struct bsg_view *view,
									int create);
GED_EXPORT extern void ged_draw_view_set_lod_bounds_update(struct bsg_view *view);
GED_EXPORT extern int ged_draw_view_has_lod_bounds_update(const struct bsg_view *view);
GED_EXPORT extern int ged_draw_mesh_lod_load_view_scene_ref(struct rt_mesh_lod *lod,
							    bsg_scene_ref visibility_ref,
							    struct bsg_view *view,
							    int reset);
GED_EXPORT extern int ged_draw_mesh_lod_load_view_scene_ref_context(struct rt_mesh_lod *lod,
								    bsg_scene_ref visibility_ref,
								    void *view_ctx,
								    int reset);
GED_EXPORT extern void ged_draw_mesh_lod_free_scene_ref(bsg_scene_ref ref);
GED_EXPORT extern bsg_scene_ref ged_draw_view_scene_root_ref(const struct bsg_view *v);
GED_EXPORT extern int ged_draw_view_has_scene_root(const struct bsg_view *v);
GED_EXPORT extern int ged_draw_view_feature_exists(struct bsg_view *view,
						   const char *name);
GED_EXPORT extern int ged_draw_view_context_feature_exists(void *view_ctx,
							   const char *name);
GED_EXPORT extern int ged_draw_view_feature_remove(struct bsg_view *view,
						   const char *name);
GED_EXPORT extern int ged_draw_view_context_feature_remove(void *view_ctx,
							   const char *name);
GED_EXPORT extern int ged_draw_view_features_remove_prefix(struct bsg_view *view,
							   const char *prefix);
GED_EXPORT extern int ged_draw_view_context_features_remove_prefix(void *view_ctx,
								   const char *prefix);
GED_EXPORT extern int ged_draw_view_feature_visible(struct bsg_view *view,
						    const char *name);
GED_EXPORT extern int ged_draw_view_context_feature_visible(void *view_ctx,
							    const char *name);
GED_EXPORT extern int ged_draw_view_feature_visible_set(struct bsg_view *view,
							const char *name,
							int visible);
GED_EXPORT extern int ged_draw_view_context_feature_visible_set(void *view_ctx,
								const char *name,
								int visible);
GED_EXPORT extern bsg_scene_ref ged_draw_view_overlay_scene_find(struct bsg_view *view,
								 const char *name);
GED_EXPORT extern void ged_draw_view_overlay_name_erase(struct bsg_view *view,
							const char *name);
GED_EXPORT extern int ged_draw_view_overlay_scene_append(struct bsg_view *view,
							 bsg_scene_ref scene);
GED_EXPORT extern int ged_draw_view_overlay_command_result_owner_set(bsg_scene_ref scene,
								     const void *owner,
								     const char *source_path);
GED_EXPORT extern bsg_scene_ref ged_draw_view_overlay_create(struct bsg_view *view,
							     const char *name);
GED_EXPORT extern bsg_scene_ref ged_draw_view_context_overlay_create(void *view_ctx,
								     const char *name);
GED_EXPORT extern int ged_draw_view_feature_depth(struct bsg_view *view,
						  const char *name,
						  int mode,
						  fastf_t *depth);
GED_EXPORT extern int ged_draw_view_context_feature_depth(void *view_ctx,
							  const char *name,
							  int mode,
							  fastf_t *depth);
GED_EXPORT extern int ged_draw_view_feature_depth_foreach(struct bsg_view *view,
							  int mode,
							  ged_draw_view_feature_depth_cb cb,
							  void *data);
GED_EXPORT extern int ged_draw_view_context_feature_depth_foreach(
	void *view_ctx,
	int mode,
	ged_draw_view_feature_depth_cb cb,
	void *data);
GED_EXPORT extern int ged_draw_view_feature_style_get(struct bsg_view *view,
						      const char *name,
						      struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_context_feature_style_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_feature_style_apply(struct bsg_view *view,
							const char *name,
							const struct ged_draw_view_feature_style *style,
							int recursive);
GED_EXPORT extern int ged_draw_view_context_feature_style_apply(
	void *view_ctx,
	const char *name,
	const struct ged_draw_view_feature_style *style,
	int recursive);
GED_EXPORT extern int ged_draw_view_feature_realize(struct bsg_view *view,
						    const char *name,
						    int recursive);
GED_EXPORT extern int ged_draw_view_context_feature_realize(void *view_ctx,
							    const char *name,
							    int recursive);
GED_EXPORT extern int ged_draw_view_indexed_face_set_replace(struct bsg_view *view,
							     const char *name,
							     int local,
							     const point_t *points,
							     size_t point_count,
							     const vect_t *normals,
							     size_t normal_count,
							     const int *indices,
							     size_t index_count,
							     const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_context_indexed_face_set_replace(
	void *view_ctx,
	const char *name,
	int local,
	const point_t *points,
	size_t point_count,
	const vect_t *normals,
	size_t normal_count,
	const int *indices,
	size_t index_count,
	const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_lines_replace(struct bsg_view *view,
						  const char *name,
						  int local,
						  const point_t *points,
						  const int *cmds,
						  size_t point_count,
						  const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_context_lines_replace(
	void *view_ctx,
	const char *name,
	int local,
	const point_t *points,
	const int *cmds,
	size_t point_count,
	const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_tcl_polygons_replace(struct bsg_view *view,
							 const char *name,
							 const point_t *points,
							 const int *cmds,
							 size_t point_count,
							 const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_context_tcl_polygons_replace(
	void *view_ctx,
	const char *name,
	const point_t *points,
	const int *cmds,
	size_t point_count,
	const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_line_layer_builder_replace(struct bsg_view *view,
							       const char *name,
							       int local,
							       const struct bg_line_layer_builder *builder);
GED_EXPORT extern int ged_draw_view_context_line_layer_builder_replace(
	void *view_ctx,
	const char *name,
	int local,
	const struct bg_line_layer_builder *builder);
GED_EXPORT extern int ged_draw_view_line_layers_replace(struct bsg_view *view,
							const char *name,
							int local,
							const struct ged_draw_view_line_layer_data *layers,
							size_t layer_count,
							const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_context_line_layers_replace(
	void *view_ctx,
	const char *name,
	int local,
	const struct ged_draw_view_line_layer_data *layers,
	size_t layer_count,
	const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_lines_create_model_annotation(struct bsg_view *view,
								  const char *name,
								  int local,
								  const point_t point);
GED_EXPORT extern int ged_draw_view_context_lines_create_model_annotation(
	void *view_ctx,
	const char *name,
	int local,
	const point_t point);
GED_EXPORT extern int ged_draw_view_lines_append_point(struct bsg_view *view,
						       const char *name,
						       const point_t point);
GED_EXPORT extern int ged_draw_view_context_lines_append_point(void *view_ctx,
							       const char *name,
							       const point_t point);
GED_EXPORT extern int ged_draw_view_label_create(struct bsg_view *view,
						 const char *name,
						 int local,
						 const char *text,
						 const point_t point,
						 const point_t target,
						 int has_target);
GED_EXPORT extern int ged_draw_view_context_label_create(void *view_ctx,
							 const char *name,
							 int local,
							 const char *text,
							 const point_t point,
							 const point_t target,
							 int has_target);
GED_EXPORT extern int ged_draw_view_labels_replace(struct bsg_view *view,
						   const char *name,
						   int local,
						   const struct ged_draw_view_label_data *labels,
						   size_t label_count);
GED_EXPORT extern int ged_draw_view_context_labels_replace(
	void *view_ctx,
	const char *name,
	int local,
	const struct ged_draw_view_label_data *labels,
	size_t label_count);
GED_EXPORT extern int ged_draw_view_tcl_labels_replace(struct bsg_view *view,
						       const char *name,
						       int draw,
						       const struct ged_draw_view_label_data *labels,
						       size_t label_count);
GED_EXPORT extern int ged_draw_view_context_tcl_labels_replace(
	void *view_ctx,
	const char *name,
	int draw,
	const struct ged_draw_view_label_data *labels,
	size_t label_count);
GED_EXPORT extern size_t ged_draw_view_label_count(struct bsg_view *view,
						   const char *name);
GED_EXPORT extern size_t ged_draw_view_context_label_count(void *view_ctx,
							   const char *name);
GED_EXPORT extern int ged_draw_view_label_copy(struct bsg_view *view,
					       const char *name,
					       size_t index,
					       struct bu_vls *text,
					       point_t point,
					       unsigned char rgb[3]);
GED_EXPORT extern int ged_draw_view_context_label_copy(void *view_ctx,
						       const char *name,
						       size_t index,
						       struct bu_vls *text,
						       point_t point,
						       unsigned char rgb[3]);
GED_EXPORT extern int ged_draw_view_label_point_set(struct bsg_view *view,
						    const char *name,
						    size_t index,
						    const point_t point);
GED_EXPORT extern int ged_draw_view_context_label_point_set(void *view_ctx,
							    const char *name,
							    size_t index,
							    const point_t point);
GED_EXPORT extern int ged_draw_view_line_style_get(struct bsg_view *view,
						   const char *name,
						   struct ged_draw_view_line_style *style);
GED_EXPORT extern int ged_draw_view_context_line_style_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_line_style *style);
GED_EXPORT extern int ged_draw_view_line_color_set(struct bsg_view *view,
						   const char *name,
						   int r,
						   int g,
						   int b);
GED_EXPORT extern int ged_draw_view_context_line_color_set(void *view_ctx,
							   const char *name,
							   int r,
							   int g,
							   int b);
GED_EXPORT extern int ged_draw_view_line_width_set(struct bsg_view *view,
						   const char *name,
						   int line_width);
GED_EXPORT extern int ged_draw_view_context_line_width_set(void *view_ctx,
							   const char *name,
							   int line_width);
GED_EXPORT extern int ged_draw_view_feature_points_copy(struct bsg_view *view,
							const char *name,
							point_t **points,
							size_t *point_count);
GED_EXPORT extern int ged_draw_view_context_feature_points_copy(void *view_ctx,
								const char *name,
								point_t **points,
								size_t *point_count);
GED_EXPORT extern int ged_draw_view_context_feature_line_command_at(
	void *view_ctx,
	const char *name,
	size_t index,
	int *out);
GED_EXPORT extern int ged_draw_view_lines_points_copy(struct bsg_view *view,
						      const char *name,
						      point_t **points,
						      size_t *point_count);
GED_EXPORT extern int ged_draw_view_context_lines_points_copy(
	void *view_ctx,
	const char *name,
	point_t **points,
	size_t *point_count);
GED_EXPORT extern int ged_draw_view_tcl_lines_replace(struct bsg_view *view,
						      const char *name,
						      const point_t *points,
						      size_t point_count,
						      const struct ged_draw_view_line_style *style);
GED_EXPORT extern int ged_draw_view_context_tcl_lines_replace(
	void *view_ctx,
	const char *name,
	const point_t *points,
	size_t point_count,
	const struct ged_draw_view_line_style *style);
GED_EXPORT extern int ged_draw_view_arrow_tip_get(struct bsg_view *view,
						  const char *name,
						  fastf_t *tip_length,
						  fastf_t *tip_width);
GED_EXPORT extern int ged_draw_view_context_arrow_tip_get(void *view_ctx,
							  const char *name,
							  fastf_t *tip_length,
							  fastf_t *tip_width);
GED_EXPORT extern int ged_draw_view_arrow_tip_set(struct bsg_view *view,
						  const char *name,
						  fastf_t tip_length,
						  fastf_t tip_width);
GED_EXPORT extern int ged_draw_view_context_arrow_tip_set(void *view_ctx,
							  const char *name,
							  fastf_t tip_length,
							  fastf_t tip_width);
GED_EXPORT extern int ged_draw_view_tcl_arrows_replace(struct bsg_view *view,
						       const char *name,
						       const point_t *points,
						       size_t point_count,
						       const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_context_tcl_arrows_replace(
	void *view_ctx,
	const char *name,
	const point_t *points,
	size_t point_count,
	const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_feature_axes_centers_copy(struct bsg_view *view,
							      const char *name,
							      point_t **centers,
							      size_t *center_count);
GED_EXPORT extern int ged_draw_view_context_feature_axes_centers_copy(
	void *view_ctx,
	const char *name,
	point_t **centers,
	size_t *center_count);
GED_EXPORT extern int ged_draw_view_tcl_axes_replace(struct bsg_view *view,
						     const char *name,
						     const point_t *centers,
						     size_t center_count,
						     fastf_t half_axes_size,
						     const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_context_tcl_axes_replace(
	void *view_ctx,
	const char *name,
	const point_t *centers,
	size_t center_count,
	fastf_t half_axes_size,
	const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_axes_create(struct bsg_view *view,
						const char *name,
						int local,
						const struct ged_draw_view_axes_state *state);
GED_EXPORT extern int ged_draw_view_context_axes_create(
	void *view_ctx,
	const char *name,
	int local,
	const struct ged_draw_view_axes_state *state);
GED_EXPORT extern int ged_draw_view_axes_state_get(struct bsg_view *view,
						   const char *name,
						   struct ged_draw_view_axes_state *state);
GED_EXPORT extern int ged_draw_view_context_axes_state_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_axes_state *state);
GED_EXPORT extern int ged_draw_view_axes_state_replace(struct bsg_view *view,
						       const char *name,
						       const struct ged_draw_view_axes_state *state);
GED_EXPORT extern int ged_draw_view_context_axes_state_replace(
	void *view_ctx,
	const char *name,
	const struct ged_draw_view_axes_state *state);
GED_EXPORT extern int ged_draw_view_polygon_ref_is_null(ged_draw_view_polygon_ref ref);
GED_EXPORT extern ged_draw_view_polygon_ref ged_draw_view_polygon_find(struct bsg_view *view,
								       const char *name);
GED_EXPORT extern ged_draw_view_polygon_ref ged_draw_view_context_polygon_find(void *view_ctx,
									       const char *name);
GED_EXPORT extern ged_draw_view_polygon_ref ged_draw_view_polygon_find_scoped(struct bsg_view *view,
									      const char *name,
									      int local_only);
GED_EXPORT extern ged_draw_view_polygon_ref ged_draw_view_context_polygon_find_scoped(void *view_ctx,
										      const char *name,
										      int local_only);
GED_EXPORT extern ged_draw_view_polygon_ref ged_draw_view_polygon_create(struct bsg_view *view,
									 const char *name,
									 int local,
									 int type,
									 const point_t screen_point);
GED_EXPORT extern ged_draw_view_polygon_ref ged_draw_view_context_polygon_create(void *view_ctx,
										 const char *name,
										 int local,
										 int type,
										 const point_t screen_point);
GED_EXPORT extern ged_draw_view_polygon_ref
ged_draw_view_polygon_import_sketch(const char *name,
				    struct db_i *dbip,
				    struct directory *dp,
				    struct bsg_view *view,
				    int local);
GED_EXPORT extern ged_draw_view_polygon_ref
ged_draw_view_context_polygon_import_sketch(const char *name,
					    struct db_i *dbip,
					    struct directory *dp,
					    void *view_ctx,
					    int local);
GED_EXPORT extern int ged_draw_view_polygon_export_sketch(struct db_i *dbip,
							 const char *name,
							 ged_draw_view_polygon_ref ref);
GED_EXPORT extern int ged_draw_view_polygon_record_get(ged_draw_view_polygon_ref ref,
						       struct ged_draw_view_polygon_record *record);
GED_EXPORT extern int ged_draw_view_polygon_has_data(ged_draw_view_polygon_ref ref);
GED_EXPORT extern int ged_draw_view_polygon_update(ged_draw_view_polygon_ref ref,
						   struct bsg_view *view,
						   int op);
GED_EXPORT extern int ged_draw_view_context_polygon_update(ged_draw_view_polygon_ref ref,
							   void *view_ctx,
							   int op);
GED_EXPORT extern int ged_draw_view_polygon_update_screen_pt(ged_draw_view_polygon_ref ref,
							     struct bsg_view *view,
							     int x,
							     int y,
							     int op);
GED_EXPORT extern int ged_draw_view_context_polygon_update_screen_pt(ged_draw_view_polygon_ref ref,
								     void *view_ctx,
								     int x,
								     int y,
								     int op);
GED_EXPORT extern int ged_draw_view_polygon_set_current(ged_draw_view_polygon_ref ref,
							long contour_i,
							long point_i);
GED_EXPORT extern int ged_draw_view_polygon_set_contour_open(ged_draw_view_polygon_ref ref,
							     long contour_i,
							     int open);
GED_EXPORT extern int ged_draw_view_polygon_set_all_contours_open(ged_draw_view_polygon_ref ref,
								  int open);
GED_EXPORT extern int ged_draw_view_polygon_area(ged_draw_view_polygon_ref ref,
						 struct bsg_view *view,
						 fastf_t *area);
GED_EXPORT extern int ged_draw_view_context_polygon_area(ged_draw_view_polygon_ref ref,
							 void *view_ctx,
							 fastf_t *area);
GED_EXPORT extern int ged_draw_view_polygon_overlap(ged_draw_view_polygon_ref ref,
						    struct bsg_view *view,
						    const char *other_name,
						    const struct bn_tol *tol,
						    int *overlap);
GED_EXPORT extern int ged_draw_view_context_polygon_overlap(ged_draw_view_polygon_ref ref,
							    void *view_ctx,
							    const char *other_name,
							    const struct bn_tol *tol,
							    int *overlap);
GED_EXPORT extern int ged_draw_view_polygon_set_fill(ged_draw_view_polygon_ref ref,
						     int fill_flag,
						     fastf_t fill_slope_x,
						     fastf_t fill_slope_y,
						     fastf_t fill_density);
GED_EXPORT extern int ged_draw_view_polygon_fill_color_get(ged_draw_view_polygon_ref ref,
							   struct bu_color *fill_color);
GED_EXPORT extern int ged_draw_view_polygon_fill_color_set(ged_draw_view_polygon_ref ref,
							   const struct bu_color *fill_color);
GED_EXPORT extern int ged_draw_view_polygon_csg(ged_draw_view_polygon_ref target,
						struct bsg_view *view,
						const char *other_name,
						bg_clip_t op);
GED_EXPORT extern int ged_draw_view_context_polygon_csg(ged_draw_view_polygon_ref target,
							void *view_ctx,
							const char *other_name,
							bg_clip_t op);

__END_DECLS

#endif /* LIBGED_BSG_GED_DRAW_VIEW_PRIVATE_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
