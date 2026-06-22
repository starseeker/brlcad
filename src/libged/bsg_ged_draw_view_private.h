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
#include "ged/defines.h"
#include "rt/view.h"
#include "vmath.h"

__BEGIN_DECLS

struct bsg_interaction_record;
struct bsg_view;
struct bsg_view_set;
struct bg_line_layer_builder;
struct bn_tol;
struct bu_ptbl;
struct db_i;
struct directory;
struct rt_mesh_lod;

typedef struct rt_view_lod_policy ged_draw_view_lod_policy;
typedef int (*ged_draw_view_feature_depth_cb)(fastf_t depth, void *data);
typedef int (*ged_draw_view_selection_path_cb)(struct bsg_view *view,
					       const char *path,
					       void *data);

enum ged_draw_view_snap_kind {
    GED_DRAW_VIEW_SNAP_GRID = 1,
    GED_DRAW_VIEW_SNAP_ENDPOINT = 2
};

struct ged_draw_view_line_style {
    int color[3];
    int line_width;
};

struct ged_draw_view_axes_state {
    point_t position;
    fastf_t size;
    int line_width;
    int color[3];
};

struct ged_draw_view_feature_style {
    int visible;
    int color_valid;
    unsigned char color[3];
    int line_width;
    int line_style;
    int arrow;
    fastf_t arrow_tip_length;
    fastf_t arrow_tip_width;
};

#define GED_DRAW_VIEW_FEATURE_STYLE_INIT { -1, 0, {0, 0, 0}, -1, -1, -1, -1.0, -1.0 }

struct ged_draw_view_label_data {
    const char *text;
    point_t point;
    int color_valid;
    unsigned char color[3];
    int line_flag;
    point_t target;
    int anchor;
    int arrow;
};

#define GED_DRAW_VIEW_LABEL_DATA_INIT { NULL, VINIT_ZERO, 0, {0, 0, 0}, 0, VINIT_ZERO, 0, 0 }

struct ged_draw_view_line_layer_data {
    const char *name;
    const point_t *points;
    const int *commands;
    size_t point_count;
    struct ged_draw_view_feature_style style;
};

#define GED_DRAW_VIEW_LINE_LAYER_DATA_INIT { NULL, NULL, NULL, 0, GED_DRAW_VIEW_FEATURE_STYLE_INIT }

enum ged_draw_view_polygon_type {
    GED_DRAW_VIEW_POLYGON_GENERAL = 0,
    GED_DRAW_VIEW_POLYGON_CIRCLE = 1,
    GED_DRAW_VIEW_POLYGON_ELLIPSE = 2,
    GED_DRAW_VIEW_POLYGON_RECTANGLE = 3,
    GED_DRAW_VIEW_POLYGON_SQUARE = 4
};

enum ged_draw_view_polygon_update {
    GED_DRAW_VIEW_POLYGON_UPDATE_DEFAULT = 0,
    GED_DRAW_VIEW_POLYGON_UPDATE_PROPS_ONLY = 1,
    GED_DRAW_VIEW_POLYGON_UPDATE_PT_SELECT = 2,
    GED_DRAW_VIEW_POLYGON_UPDATE_PT_SELECT_CLEAR = 3,
    GED_DRAW_VIEW_POLYGON_UPDATE_PT_MOVE = 4,
    GED_DRAW_VIEW_POLYGON_UPDATE_PT_APPEND = 5
};

typedef struct ged_draw_view_polygon_ref {
    uintptr_t token;
    uint64_t revision;
} ged_draw_view_polygon_ref;

#define GED_DRAW_VIEW_POLYGON_REF_NULL_INIT {0, 0}

#ifdef __cplusplus
#  define GED_DRAW_VIEW_POLYGON_REF_NULL ged_draw_view_polygon_ref{0, 0}
#else
#  define GED_DRAW_VIEW_POLYGON_REF_NULL ((ged_draw_view_polygon_ref){0, 0})
#endif

struct ged_draw_view_polygon_record {
    ged_draw_view_polygon_ref ref;
    const char *name;
    int type;
    int fill_flag;
    vect2d_t fill_dir;
    fastf_t fill_delta;
    struct bu_color fill_color;
    unsigned char edge_color[3];
    long curr_contour_i;
    long curr_point_i;
    int first_contour_open;
    size_t contour_count;
    size_t point_count;
    point_t origin_point;
    plane_t vp;
    fastf_t vZ;
    void *user_data;
};

GED_EXPORT extern void ged_draw_view_info_from_bsg(struct rt_view_info *view_info,
						   const struct bsg_view *view);
GED_EXPORT extern fastf_t ged_draw_view_perspective_from_bsg(const struct bsg_view *view);
GED_EXPORT extern fastf_t ged_draw_view_scale_from_bsg(const struct bsg_view *view);
GED_EXPORT extern int ged_draw_view_lod_policy_from_bsg(ged_draw_view_lod_policy *policy,
							const struct bsg_view *view);
GED_EXPORT extern int ged_draw_view_lod_policy_apply_bsg(struct bsg_view *view,
							 const ged_draw_view_lod_policy *policy);
GED_EXPORT extern int ged_draw_view_lod_policy_apply_bsg_bot_threshold(
	struct bsg_view *view,
	const ged_draw_view_lod_policy *policy,
	size_t bot_threshold);
GED_EXPORT extern int ged_draw_view_autoview_default_bsg(struct bsg_view *view,
							 int all_view_objs);
GED_EXPORT extern int ged_draw_view_hud_sync(struct bsg_view *view);
GED_EXPORT extern int ged_draw_view_selection_available(struct bsg_view *view);
GED_EXPORT extern size_t ged_draw_view_selection_count(struct bsg_view *view);
GED_EXPORT extern int ged_draw_view_selection_path_foreach(struct bsg_view *view,
							   ged_draw_view_selection_path_cb cb,
							   void *data);
GED_EXPORT extern int ged_draw_view_selection_clear(struct bsg_view *view);
GED_EXPORT extern int ged_draw_view_selection_contains_record(
	struct bsg_view *view,
	const struct bsg_interaction_record *record);
GED_EXPORT extern int ged_draw_view_selection_add_record(
	struct bsg_view *view,
	const struct bsg_interaction_record *record);
GED_EXPORT extern int ged_draw_view_selection_set_record(
	struct bsg_view *view,
	const struct bsg_interaction_record *record);
GED_EXPORT extern int ged_draw_view_snap_first_candidate(struct bsg_view *view,
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
GED_EXPORT extern void ged_draw_mesh_lod_free_scene_ref(bsg_scene_ref ref);
GED_EXPORT extern bsg_scene_ref ged_draw_view_scene_root_ref(const struct bsg_view *v);
GED_EXPORT extern int ged_draw_view_has_scene_root(const struct bsg_view *v);
GED_EXPORT extern int ged_draw_view_feature_exists(struct bsg_view *view,
						   const char *name);
GED_EXPORT extern int ged_draw_view_feature_remove(struct bsg_view *view,
						   const char *name);
GED_EXPORT extern int ged_draw_view_features_remove_prefix(struct bsg_view *view,
							   const char *prefix);
GED_EXPORT extern int ged_draw_view_feature_visible(struct bsg_view *view,
						    const char *name);
GED_EXPORT extern bsg_scene_ref ged_draw_view_overlay_create(struct bsg_view *view,
							     const char *name);
GED_EXPORT extern int ged_draw_view_feature_depth(struct bsg_view *view,
						  const char *name,
						  int mode,
						  fastf_t *depth);
GED_EXPORT extern int ged_draw_view_feature_depth_foreach(struct bsg_view *view,
							  int mode,
							  ged_draw_view_feature_depth_cb cb,
							  void *data);
GED_EXPORT extern int ged_draw_view_feature_style_get(struct bsg_view *view,
						      const char *name,
						      struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_feature_style_apply(struct bsg_view *view,
							const char *name,
							const struct ged_draw_view_feature_style *style,
							int recursive);
GED_EXPORT extern int ged_draw_view_feature_realize(struct bsg_view *view,
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
GED_EXPORT extern int ged_draw_view_lines_replace(struct bsg_view *view,
						  const char *name,
						  int local,
						  const point_t *points,
						  const int *cmds,
						  size_t point_count,
						  const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_line_layer_builder_replace(struct bsg_view *view,
							       const char *name,
							       int local,
							       const struct bg_line_layer_builder *builder);
GED_EXPORT extern int ged_draw_view_line_layers_replace(struct bsg_view *view,
							const char *name,
							int local,
							const struct ged_draw_view_line_layer_data *layers,
							size_t layer_count,
							const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_view_lines_create_model_annotation(struct bsg_view *view,
								  const char *name,
								  int local,
								  const point_t point);
GED_EXPORT extern int ged_draw_view_lines_append_point(struct bsg_view *view,
						       const char *name,
						       const point_t point);
GED_EXPORT extern int ged_draw_view_label_create(struct bsg_view *view,
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
GED_EXPORT extern int ged_draw_view_line_style_get(struct bsg_view *view,
						   const char *name,
						   struct ged_draw_view_line_style *style);
GED_EXPORT extern int ged_draw_view_line_color_set(struct bsg_view *view,
						   const char *name,
						   int r,
						   int g,
						   int b);
GED_EXPORT extern int ged_draw_view_line_width_set(struct bsg_view *view,
						   const char *name,
						   int line_width);
GED_EXPORT extern int ged_draw_view_lines_points_copy(struct bsg_view *view,
						      const char *name,
						      point_t **points,
						      size_t *point_count);
GED_EXPORT extern int ged_draw_view_tcl_lines_replace(struct bsg_view *view,
						      const char *name,
						      const point_t *points,
						      size_t point_count,
						      const struct ged_draw_view_line_style *style);
GED_EXPORT extern int ged_draw_view_axes_create(struct bsg_view *view,
						const char *name,
						int local,
						const struct ged_draw_view_axes_state *state);
GED_EXPORT extern int ged_draw_view_axes_state_get(struct bsg_view *view,
						   const char *name,
						   struct ged_draw_view_axes_state *state);
GED_EXPORT extern int ged_draw_view_axes_state_replace(struct bsg_view *view,
						       const char *name,
						       const struct ged_draw_view_axes_state *state);
GED_EXPORT extern int ged_draw_view_polygon_ref_is_null(ged_draw_view_polygon_ref ref);
GED_EXPORT extern ged_draw_view_polygon_ref ged_draw_view_polygon_find(struct bsg_view *view,
								       const char *name);
GED_EXPORT extern ged_draw_view_polygon_ref ged_draw_view_polygon_find_scoped(struct bsg_view *view,
									      const char *name,
									      int local_only);
GED_EXPORT extern ged_draw_view_polygon_ref ged_draw_view_polygon_create(struct bsg_view *view,
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
GED_EXPORT extern int ged_draw_view_polygon_export_sketch(struct db_i *dbip,
							 const char *name,
							 ged_draw_view_polygon_ref ref);
GED_EXPORT extern int ged_draw_view_polygon_record_get(ged_draw_view_polygon_ref ref,
						       struct ged_draw_view_polygon_record *record);
GED_EXPORT extern int ged_draw_view_polygon_has_data(ged_draw_view_polygon_ref ref);
GED_EXPORT extern int ged_draw_view_polygon_update(ged_draw_view_polygon_ref ref,
						   struct bsg_view *view,
						   int op);
GED_EXPORT extern int ged_draw_view_polygon_update_screen_pt(ged_draw_view_polygon_ref ref,
							     struct bsg_view *view,
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
GED_EXPORT extern int ged_draw_view_polygon_overlap(ged_draw_view_polygon_ref ref,
						    struct bsg_view *view,
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
