/*                    V I E W _ L E G A C Y _ B S G . H
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
/** @file rt/view_legacy_bsg.h
 *
 * Transitional adapters for callers that still own BSG view state while RT
 * primitive callbacks move to RT-owned view snapshots.
 */

#ifndef RT_VIEW_LEGACY_BSG_H
#define RT_VIEW_LEGACY_BSG_H

#include "common.h"
#include "bsg/scene_builder.h"
#include "bsg/view_state.h"
#include "rt/defines.h"
#include "rt/view.h"

__BEGIN_DECLS

struct bsg_view;
struct bsg_adc_state;
struct bsg_axes;
struct bsg_grid_state;
struct bsg_interactive_rect_state;
struct bsg_other_state;
struct bsg_params_state;
struct bsg_selection;
struct bsg_snap_result;
struct rt_mesh_lod;
struct bsg_view_knobs;
struct bsg_view_set;
struct bu_data_hash_state;
struct bu_ptbl;
struct bu_vls;

typedef void (*rt_view_bounds_update_callback_bsg_t)(struct bsg_view *);

enum rt_view_knobs_category_bsg {
    RT_VIEW_KNOBS_ALL_BSG = 0,
    RT_VIEW_KNOBS_RATE_BSG = 1,
    RT_VIEW_KNOBS_ABS_BSG = 2
};

#define RT_VIEW_REFRESH_VIEW_BSG BSG_VIEW_REFRESH_VIEW
#define RT_VIEW_REFRESH_DRAW_BSG BSG_VIEW_REFRESH_DRAW
#define RT_VIEW_REFRESH_EDIT_BSG BSG_VIEW_REFRESH_EDIT
#define RT_VIEW_REFRESH_FRAMEBUFFER_BSG BSG_VIEW_REFRESH_FRAMEBUFFER
#define RT_VIEW_REFRESH_OVERLAY_BSG BSG_VIEW_REFRESH_OVERLAY
#define RT_VIEW_REFRESH_FORCE_BSG BSG_VIEW_REFRESH_FORCE
#define RT_VIEW_REFRESH_ALL_BSG BSG_VIEW_REFRESH_ALL

#define RT_VIEW_SNAP_SHARED_BSG BSG_SNAP_SHARED
#define RT_VIEW_SNAP_LOCAL_BSG BSG_SNAP_LOCAL
#define RT_VIEW_SNAP_DB_BSG BSG_SNAP_DB
#define RT_VIEW_SNAP_VIEW_BSG BSG_SNAP_VIEW
#define RT_VIEW_SNAP_TCL_BSG BSG_SNAP_TCL

#define RT_VIEW_CLEAR_DB_BSG BSG_SOURCE_DB
#define RT_VIEW_CLEAR_VIEW_BSG BSG_SOURCE_VIEW
#define RT_VIEW_CLEAR_LOCAL_BSG BSG_SOURCE_LOCAL

#define RT_VIEW_SNAP_KIND_GRID_BSG BSG_SNAP_KIND_GRID
#define RT_VIEW_SNAP_KIND_ENDPOINT_BSG BSG_SNAP_KIND_ENDPOINT
#define RT_VIEW_SNAP_KIND_MIDPOINT_BSG BSG_SNAP_KIND_MIDPOINT
#define RT_VIEW_SNAP_KIND_INTERSECTION_BSG BSG_SNAP_KIND_INTERSECTION
#define RT_VIEW_SNAP_KIND_PERPENDICULAR_BSG BSG_SNAP_KIND_PERPENDICULAR
#define RT_VIEW_SNAP_KIND_TANGENT_BSG BSG_SNAP_KIND_TANGENT
#define RT_VIEW_SNAP_KIND_OVERLAY_HANDLE_BSG BSG_SNAP_KIND_OVERLAY_HANDLE

RT_EXPORT extern void
rt_view_info_from_bsg(struct rt_view_info *info, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_width_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_height_from_bsg(const struct bsg_view *v);

RT_EXPORT extern fastf_t
rt_view_radius_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_dimensions_set_bsg(struct bsg_view *v, int width, int height);

RT_EXPORT extern int
rt_view_screen_to_view_from_bsg(fastf_t *fx,
				fastf_t *fy,
				struct bsg_view *v,
				fastf_t x,
				fastf_t y);

RT_EXPORT extern int
rt_view_screen_point_from_bsg(point_t point,
			      struct bsg_view *v,
			      fastf_t x,
			      fastf_t y);

RT_EXPORT extern int
rt_view_current_point_from_bsg(point_t point, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_current_point_set_bsg(struct bsg_view *v, const point_t point);

RT_EXPORT extern int
rt_view_previous_mouse_from_bsg(fastf_t *x, fastf_t *y, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_previous_mouse_set_bsg(struct bsg_view *v, fastf_t x, fastf_t y);

struct rt_view_mouse_delta_settings {
    fastf_t min_delta;
    fastf_t max_delta;
    fastf_t rotate_scale;
    fastf_t scale_scale;
};

#define RT_VIEW_MOUSE_DELTA_SETTINGS_INIT { 0.0, 0.0, 0.0, 0.0 }

RT_EXPORT extern int
rt_view_mouse_delta_settings_from_bsg(struct rt_view_mouse_delta_settings *settings,
				      const struct bsg_view *v);

RT_EXPORT extern int
rt_view_mouse_state_set_bsg(struct bsg_view *v, int x, int y);

RT_EXPORT extern int
rt_view_unique_object_name_bsg(struct bu_vls *oname,
			       const char *seed,
			       struct bsg_view *v);

RT_EXPORT extern int
rt_view_interactive_rect_from_bsg(struct bsg_interactive_rect_state *record,
				  const struct bsg_view *v);

RT_EXPORT extern int
rt_view_interactive_rect_set_bsg(struct bsg_view *v,
				 const struct bsg_interactive_rect_state *record);

RT_EXPORT extern int
rt_view_adc_from_bsg(struct bsg_adc_state *record,
		     const struct bsg_view *v);

RT_EXPORT extern int
rt_view_adc_set_bsg(struct bsg_view *v,
		    const struct bsg_adc_state *record);

RT_EXPORT extern int
rt_view_grid_from_bsg(struct bsg_grid_state *record,
		      const struct bsg_view *v);

RT_EXPORT extern int
rt_view_grid_set_bsg(struct bsg_view *v,
		     const struct bsg_grid_state *record);

RT_EXPORT extern int
rt_view_model_axes_from_bsg(struct bsg_axes *record,
			    const struct bsg_view *v);

RT_EXPORT extern int
rt_view_model_axes_set_bsg(struct bsg_view *v,
			   const struct bsg_axes *record);

RT_EXPORT extern int
rt_view_view_axes_from_bsg(struct bsg_axes *record,
			   const struct bsg_view *v);

RT_EXPORT extern int
rt_view_view_axes_set_bsg(struct bsg_view *v,
			  const struct bsg_axes *record);

RT_EXPORT extern int
rt_view_center_dot_from_bsg(struct bsg_other_state *record,
			    const struct bsg_view *v);

RT_EXPORT extern int
rt_view_center_dot_set_bsg(struct bsg_view *v,
			   const struct bsg_other_state *record);

RT_EXPORT extern int
rt_view_scale_overlay_from_bsg(struct bsg_other_state *record,
			       const struct bsg_view *v);

RT_EXPORT extern int
rt_view_scale_overlay_set_bsg(struct bsg_view *v,
			      const struct bsg_other_state *record);

RT_EXPORT extern int
rt_view_params_from_bsg(struct bsg_params_state *record,
			const struct bsg_view *v);

RT_EXPORT extern int
rt_view_params_set_bsg(struct bsg_view *v,
		       const struct bsg_params_state *record);

RT_EXPORT extern int
rt_view_refresh_request_bsg(struct bsg_view *v, uint32_t flags);

RT_EXPORT extern int
rt_view_refresh_dirty_from_bsg(const struct bsg_view *v);

RT_EXPORT extern uint32_t
rt_view_refresh_consume_bsg(struct bsg_view *v);

RT_EXPORT extern int
rt_view_refresh_complete_bsg(struct bsg_view *v);

RT_EXPORT extern int
rt_view_refresh_enabled_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_refresh_enabled_set_bsg(struct bsg_view *v, int enabled);

RT_EXPORT extern int
rt_view_refresh_suppressed_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_refresh_suppress_begin_bsg(struct bsg_view *v);

RT_EXPORT extern int
rt_view_refresh_suppress_end_bsg(struct bsg_view *v);

RT_EXPORT extern int
rt_view_refresh_drawn_count_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_refresh_drawn_count_set_bsg(struct bsg_view *v, int count);

RT_EXPORT extern int
rt_view_orientation_quat_from_bsg(quat_t orientation,
				  const struct bsg_view *v);

RT_EXPORT extern int
rt_view_aet_from_bsg(vect_t aet, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_aet_set_bsg(struct bsg_view *v, const vect_t aet);

RT_EXPORT extern int
rt_view_aet_state_set_bsg(struct bsg_view *v, const vect_t aet);

RT_EXPORT extern int
rt_view_update_bsg(struct bsg_view *v);

#define RT_VIEW_AUTOVIEW_SCALE_DEFAULT -1

RT_EXPORT extern int
rt_view_autoview_bsg(struct bsg_view *v, fastf_t scale, int all_view_objs);

RT_EXPORT extern int
rt_view_autoview_bounds_bsg(struct bsg_view *v,
			    fastf_t scale,
			    const point_t min,
			    const point_t max);

#define RT_VIEW_ADJUST_IDLE   0x000ULL
#define RT_VIEW_ADJUST_ROT    0x001ULL
#define RT_VIEW_ADJUST_TRANS  0x002ULL
#define RT_VIEW_ADJUST_SCALE  0x004ULL
#define RT_VIEW_ADJUST_CENTER 0x008ULL
#define RT_VIEW_ADJUST_CON_X  0x010ULL
#define RT_VIEW_ADJUST_CON_Y  0x020ULL
#define RT_VIEW_ADJUST_CON_Z  0x040ULL
#define RT_VIEW_ADJUST_CON_GRID 0x080ULL
#define RT_VIEW_ADJUST_CON_LINES 0x100ULL

RT_EXPORT extern int
rt_view_adjust_bsg(struct bsg_view *v,
		   int dx,
		   int dy,
		   point_t keypoint,
		   int mode,
		   unsigned long long flags);

RT_EXPORT extern unsigned long long
rt_view_hash_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_snap_candidates_bsg(struct bsg_view *v,
			    point_t sample,
			    double tol,
			    unsigned long long kinds,
			    struct bsg_snap_result *out);

RT_EXPORT extern int
rt_view_snap_point_2d_bsg(struct bsg_view *v,
			  fastf_t *vx,
			  fastf_t *vy,
			  unsigned long long kinds);

RT_EXPORT extern int
rt_view_snap_grid_2d_bsg(struct bsg_view *v, fastf_t *vx, fastf_t *vy);

RT_EXPORT extern struct bsg_selection *
rt_view_selection_bsg(struct bsg_view *v);

RT_EXPORT extern const struct bsg_selection *
rt_view_selection_const_bsg(const struct bsg_view *v);

RT_EXPORT extern struct bu_ptbl *
rt_view_set_views_bsg(struct bsg_view_set *s);

RT_EXPORT extern struct bsg_view *
rt_view_set_find_view_bsg(struct bsg_view_set *s, const char *name);

RT_EXPORT extern void
rt_view_set_init_bsg(struct bsg_view_set *s);

RT_EXPORT extern void
rt_view_set_free_bsg(struct bsg_view_set *s);

RT_EXPORT extern void
rt_view_init_bsg(struct bsg_view *v, struct bsg_view_set *s);

RT_EXPORT extern void
rt_view_free_bsg(struct bsg_view *v);

RT_EXPORT extern void
rt_view_set_add_view_bsg(struct bsg_view_set *s, struct bsg_view *v);

RT_EXPORT extern void
rt_view_set_remove_view_bsg(struct bsg_view_set *s, struct bsg_view *v);

RT_EXPORT extern int
rt_view_knobs_reset_bsg(struct bsg_view *v, int category);

RT_EXPORT extern int
rt_view_knob_state_reset_bsg(struct bsg_view_knobs *knobs, int category);

RT_EXPORT extern unsigned long long
rt_view_knob_state_hash_bsg(struct bsg_view_knobs *knobs,
			    struct bu_data_hash_state *state);

RT_EXPORT extern unsigned long long
rt_view_knobs_hash_bsg(struct bsg_view *v,
		       struct bu_data_hash_state *state);

RT_EXPORT extern int
rt_view_knobs_cmd_process_bsg(vect_t *rvec,
			      int *do_rot,
			      vect_t *tvec,
			      int *do_tran,
			      struct bsg_view *v,
			      const char *cmd,
			      fastf_t factor,
			      char origin,
			      int model_flag,
			      int incr_flag);

RT_EXPORT extern int
rt_view_knobs_translate_bsg(struct bsg_view *v,
			    const vect_t tvec,
			    int model_flag);

RT_EXPORT extern int
rt_view_knobs_rotate_bsg(struct bsg_view *v,
			 const vect_t rvec,
			 char origin,
			 char coords,
			 const matp_t obj_rot,
			 const pointp_t pvt_pt);

RT_EXPORT extern int
rt_view_knobs_update_rate_flags_bsg(struct bsg_view *v);

RT_EXPORT extern int
rt_view_is_independent_bsg(const struct bsg_view *v);

RT_EXPORT extern bsg_scene_ref
rt_view_independent_scope_ref_bsg(struct bsg_view *v, int create);

RT_EXPORT extern void
rt_view_independent_scope_destroy_bsg(struct bsg_view *v);

RT_EXPORT extern void
rt_view_init_copy_bsg(struct bsg_view *dest,
		      const struct bsg_view *src,
		      struct bsg_view_set *s);

RT_EXPORT extern void
rt_view_tclcad_data_init_bsg(struct bsg_data_tclcad *d);

RT_EXPORT extern size_t
rt_view_clear_bsg(struct bsg_view *v, int flags);

RT_EXPORT extern fastf_t
rt_view_perspective_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_perspective_set_bsg(struct bsg_view *v, fastf_t perspective);

RT_EXPORT extern int
rt_view_model2view_from_bsg(mat_t model2view, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_model2view_set_bsg(struct bsg_view *v, const mat_t model2view);

RT_EXPORT extern int
rt_view_view2model_from_bsg(mat_t view2model, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_view2model_set_bsg(struct bsg_view *v, const mat_t view2model);

RT_EXPORT extern int
rt_view_pmodel2view_from_bsg(mat_t pmodel2view, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_pmodel2view_set_bsg(struct bsg_view *v, const mat_t pmodel2view);

RT_EXPORT extern int
rt_view_pmat_from_bsg(mat_t pmat, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_pmat_set_bsg(struct bsg_view *v, const mat_t pmat);

RT_EXPORT extern int
rt_view_rotation_from_bsg(mat_t rotation, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_rotation_set_bsg(struct bsg_view *v, const mat_t rotation);

RT_EXPORT extern int
rt_view_center_from_bsg(mat_t center, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_center_vec_set_bsg(struct bsg_view *v, const point_t center);

RT_EXPORT extern int
rt_view_plane_from_bsg(plane_t *plane, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_lod_policy_from_bsg(struct rt_view_lod_policy *policy,
			    const struct bsg_view *v);

RT_EXPORT extern int
rt_view_lod_policy_apply_bsg(struct bsg_view *v,
			     const struct rt_view_lod_policy *policy);

RT_EXPORT extern int
rt_view_lod_policy_copy_bsg(struct bsg_view *dst, const struct bsg_view *src);

RT_EXPORT extern void
rt_view_lod_bounds_update_bsg(struct bsg_view *v);

RT_EXPORT extern void
rt_view_lod_bounds_callback_set_bsg(struct bsg_view *v);

RT_EXPORT extern int
rt_view_lod_bounds_callback_is_bsg(const struct bsg_view *v);

RT_EXPORT extern rt_view_bounds_update_callback_bsg_t
rt_view_bounds_update_callback_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_bounds_update_callback_set_bsg(struct bsg_view *v,
				       rt_view_bounds_update_callback_bsg_t callback);

RT_EXPORT extern int
rt_view_bounds_update_callback_call_bsg(struct bsg_view *v);

RT_EXPORT extern fastf_t
rt_view_scale_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_scale_set_bsg(struct bsg_view *v, fastf_t scale);

RT_EXPORT extern fastf_t *
rt_view_scale_storage_from_bsg(struct bsg_view *v);

RT_EXPORT extern int
rt_view_scale_state_set_bsg(struct bsg_view *v,
			    fastf_t scale,
			    fastf_t initial_scale,
			    fastf_t absolute_scale,
			    fastf_t size,
			    fastf_t inverse_size);

RT_EXPORT extern fastf_t
rt_view_initial_scale_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_initial_scale_set_bsg(struct bsg_view *v, fastf_t scale);

RT_EXPORT extern fastf_t
rt_view_absolute_scale_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_absolute_scale_set_bsg(struct bsg_view *v, fastf_t scale);

RT_EXPORT extern fastf_t
rt_view_local2base_from_bsg(const struct bsg_view *v);

RT_EXPORT extern fastf_t
rt_view_base2local_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_unit_conversion_set_bsg(struct bsg_view *v,
				fastf_t local2base,
				fastf_t base2local);

RT_EXPORT extern fastf_t
rt_view_size_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_size_set_bsg(struct bsg_view *v, fastf_t size);

RT_EXPORT extern fastf_t
rt_view_inverse_size_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_eye_pos_from_bsg(point_t eye_pos, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_eye_pos_set_bsg(struct bsg_view *v, const point_t eye_pos);

RT_EXPORT extern int
rt_view_keypoint_from_bsg(point_t keypoint, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_keypoint_set_bsg(struct bsg_view *v, const point_t keypoint);

RT_EXPORT extern char
rt_view_rotate_about_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_rotate_about_set_bsg(struct bsg_view *v, char rotate_about);

RT_EXPORT extern char
rt_view_coord_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_coord_set_bsg(struct bsg_view *v, char coord);

RT_EXPORT extern int
rt_view_snap_lines_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_snap_lines_set_bsg(struct bsg_view *v, int enabled);

RT_EXPORT extern int
rt_view_snap_source_flags_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_snap_source_flags_set_bsg(struct bsg_view *v, int flags);

RT_EXPORT extern unsigned long long
rt_view_snap_kind_mask_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_snap_exclude_feature_from_bsg(bsg_feature_ref *ref,
				      const struct bsg_view *v);

RT_EXPORT extern int
rt_view_snap_exclude_feature_set_bsg(struct bsg_view *v, bsg_feature_ref ref);

RT_EXPORT extern int
rt_view_snap_exclude_feature_clear_bsg(struct bsg_view *v);

RT_EXPORT extern unsigned long long
rt_view_prepare_tcl_snap_bsg(struct bsg_view *v);

RT_EXPORT extern int
rt_view_center_linesnap_bsg(struct bsg_view *v);

RT_EXPORT extern int
rt_view_zclip_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_zclip_set_bsg(struct bsg_view *v, int zclip);

RT_EXPORT extern int
rt_view_framebuffer_mode_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_framebuffer_mode_set_bsg(struct bsg_view *v, int mode);

RT_EXPORT extern int
rt_view_cleared_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_cleared_set_bsg(struct bsg_view *v, int cleared);

RT_EXPORT extern int
rt_view_settings_shared_bsg(const struct bsg_view *a, const struct bsg_view *b);

RT_EXPORT extern double
rt_view_snap_tolerance_factor_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_snap_tolerance_factor_set_bsg(struct bsg_view *v, double factor);

RT_EXPORT extern int
rt_mesh_lod_load_view_scene_ref_bsg(struct rt_mesh_lod *lod,
				    bsg_scene_ref visibility_ref,
				    struct bsg_view *v,
				    int reset);

RT_EXPORT extern void
rt_mesh_lod_free_scene_ref_bsg(bsg_scene_ref ref);

__END_DECLS

#endif /* RT_VIEW_LEGACY_BSG_H */
