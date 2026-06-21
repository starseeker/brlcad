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

#include "bsg/scene_builder.h"
#include "bsg/scene_object.h"
#include "ged/defines.h"
#include "rt/view.h"
#include "vmath.h"

__BEGIN_DECLS

struct bsg_selection;
struct bsg_view;
struct bsg_view_set;
struct bu_ptbl;
struct rt_mesh_lod;

typedef struct rt_view_lod_policy ged_draw_view_lod_policy;

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
GED_EXPORT extern struct bsg_selection *ged_draw_view_selection_bsg(struct bsg_view *view);
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
GED_EXPORT extern int ged_draw_view_feature_visible(struct bsg_view *view,
						    const char *name);
GED_EXPORT extern bsg_scene_ref ged_draw_view_overlay_create(struct bsg_view *view,
							     const char *name);
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
