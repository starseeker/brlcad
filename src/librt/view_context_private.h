/*              V I E W _ C O N T E X T _ P R I V A T E . H
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
/** @file view_context_private.h
 *
 * Internal native/libbv-backed RT view context helpers.  Public RT view
 * context entry points remain in view_legacy_bsg.c while retained BSG
 * adapters are being split away.
 */

#ifndef LIBRT_VIEW_CONTEXT_PRIVATE_H
#define LIBRT_VIEW_CONTEXT_PRIVATE_H

#include "common.h"

#include "rt/view.h"

__BEGIN_DECLS

int _rt_view_context_native_is(const void *ctx);
int _rt_view_set_context_native_is(const void *view_set_ctx);

void *_rt_view_context_native_create(void);
void *_rt_view_context_native_create_with_set(void *view_set_ctx);
void *_rt_view_context_native_create_copy_with_set(const void *src_ctx,
	void *view_set_ctx);
void _rt_view_context_native_free(void *ctx);
int _rt_view_context_native_release_storage(void *ctx);

void *_rt_view_set_context_native_create(void);
void _rt_view_set_context_native_destroy(void *view_set_ctx);
void _rt_view_set_context_native_init(void *view_set_ctx);
void _rt_view_set_context_native_free(void *view_set_ctx);
struct bu_ptbl *_rt_view_set_context_native_views(void *view_set_ctx);
void *_rt_view_set_context_native_recycle_pool(void *view_set_ctx);
void *_rt_view_set_context_native_find_view(void *view_set_ctx,
	const char *name);
int _rt_view_set_context_native_add(void *view_set_ctx, void *view_ctx);
int _rt_view_set_context_native_remove(void *view_set_ctx, void *view_ctx);
int _rt_view_context_native_view_set_attach(void *ctx, void *view_set_ctx);

void *_rt_view_context_native_user_data_get(const void *ctx);
int _rt_view_context_native_user_data_set(void *ctx, void *user_data);
int _rt_view_context_native_tclcad_data_set(void *ctx, void *tcl_data);
int _rt_view_context_native_callbacks_set(void *ctx, struct bu_ptbl *callbacks);
int _rt_view_context_native_update_callback_set(void *ctx,
	rt_view_context_update_callback_t callback, void *data);
int _rt_view_context_native_name_set(void *ctx, const char *name);
const char *_rt_view_context_native_name_get(const void *ctx);
int _rt_view_context_native_width_get(const void *ctx);
int _rt_view_context_native_height_get(const void *ctx);
int _rt_view_context_native_dimensions_set(void *ctx, int width, int height);
void *_rt_view_context_native_display_manager_get(const void *ctx);
int _rt_view_context_native_display_manager_set(void *ctx, void *dmp);
int _rt_view_context_native_edit_matrix_set(void *ctx, matp_t edit_mat);
int _rt_view_context_native_edit_matrix_clear(void *ctx);
uint64_t _rt_view_context_native_frame_revision_get(const void *ctx);
uint64_t _rt_view_context_native_frame_revision_bump(void *ctx);
int _rt_view_context_native_model_matrices_identity(void *ctx);
int _rt_view_context_native_frametime_set(void *ctx, uint64_t frametime);

int _rt_view_context_native_refresh_request(void *ctx, uint32_t flags);
int _rt_view_context_native_refresh_dirty_get(const void *ctx);
uint32_t _rt_view_context_native_refresh_consume(void *ctx);
int _rt_view_context_native_refresh_complete(void *ctx);
int _rt_view_context_native_refresh_enabled_get(const void *ctx);
int _rt_view_context_native_refresh_enabled_set(void *ctx, int enabled);
int _rt_view_context_native_refresh_suppressed_get(const void *ctx);
int _rt_view_context_native_refresh_suppress_begin(void *ctx);
int _rt_view_context_native_refresh_suppress_end(void *ctx);
int _rt_view_context_native_refresh_drawn_count_get(const void *ctx);
int _rt_view_context_native_refresh_drawn_count_set(void *ctx, int count);

fastf_t _rt_view_context_native_scale_get(const void *ctx);
int _rt_view_context_native_scale_set(void *ctx, fastf_t scale);
fastf_t *_rt_view_context_native_scale_storage_get(void *ctx);
int _rt_view_context_native_scale_state_set(void *ctx, fastf_t scale,
	fastf_t initial_scale, fastf_t absolute_scale, fastf_t size,
	fastf_t inverse_size);
fastf_t _rt_view_context_native_initial_scale_get(const void *ctx);
int _rt_view_context_native_initial_scale_set(void *ctx, fastf_t scale);
fastf_t _rt_view_context_native_absolute_scale_get(const void *ctx);
int _rt_view_context_native_absolute_scale_set(void *ctx, fastf_t scale);
int _rt_view_context_native_unit_conversion_set(void *ctx,
	fastf_t local2base, fastf_t base2local);
fastf_t _rt_view_context_native_local2base_get(const void *ctx);
fastf_t _rt_view_context_native_base2local_get(const void *ctx);
fastf_t _rt_view_context_native_size_get(const void *ctx);
int _rt_view_context_native_size_set(void *ctx, fastf_t size);
fastf_t _rt_view_context_native_inverse_size_get(const void *ctx);
fastf_t _rt_view_context_native_perspective_get(const void *ctx);
int _rt_view_context_native_perspective_set(void *ctx, fastf_t perspective);
fastf_t _rt_view_context_native_radius_get(const void *ctx);

void _rt_view_context_native_info_get(struct rt_view_info *info,
	const void *ctx);
int _rt_view_context_native_obb_get(const void *ctx, point_t center,
	vect_t extent1, vect_t extent2, vect_t extent3);
int _rt_view_context_native_model2view_get(mat_t model2view, const void *ctx);
int _rt_view_context_native_model2view_set(void *ctx, const mat_t model2view);
int _rt_view_context_native_view2model_get(mat_t view2model, const void *ctx);
int _rt_view_context_native_view2model_set(void *ctx, const mat_t view2model);
int _rt_view_context_native_pmodel2view_get(mat_t pmodel2view,
	const void *ctx);
int _rt_view_context_native_pmodel2view_set(void *ctx,
	const mat_t pmodel2view);
int _rt_view_context_native_pmat_get(mat_t pmat, const void *ctx);
int _rt_view_context_native_pmat_set(void *ctx, const mat_t pmat);
int _rt_view_context_native_rotation_get(mat_t rotation, const void *ctx);
int _rt_view_context_native_rotation_set(void *ctx, const mat_t rotation);
int _rt_view_context_native_center_get(mat_t center, const void *ctx);
int _rt_view_context_native_center_mat_set(void *ctx, const mat_t center);
int _rt_view_context_native_center_set(void *ctx, const point_t center);
int _rt_view_context_native_plane_get(plane_t *plane, const void *ctx);
int _rt_view_context_native_orientation_quat_get(quat_t orientation,
	const void *ctx);
int _rt_view_context_native_aet_get(vect_t aet, const void *ctx);
int _rt_view_context_native_aet_set(void *ctx, const vect_t aet);
int _rt_view_context_native_aet_state_set(void *ctx, const vect_t aet);
int _rt_view_context_native_update(void *ctx);
int _rt_view_context_native_autoview_bounds(void *ctx, fastf_t scale,
	const point_t min, const point_t max);
int _rt_view_context_native_adjust(void *ctx, int dx, int dy,
	point_t keypoint, int mode, unsigned long long flags);
unsigned long long _rt_view_context_native_hash(const void *ctx);
size_t _rt_view_context_native_clear(void *ctx, int flags);

rt_view_scene_ref _rt_view_context_native_independent_scope_ref(void *ctx,
	int create);
int _rt_view_context_native_is_independent(const void *ctx);
int _rt_view_context_native_independent_scope_is_null(void *ctx,
	int create);
void _rt_view_context_native_independent_scope_destroy(void *ctx);
rt_view_scene_ref _rt_view_context_native_scene_root_ref(const void *ctx);
int _rt_view_context_native_scene_root_ref_attach(void *ctx,
	rt_view_scene_ref root_ref);
int _rt_view_context_native_scene_attached(const void *ctx);
int _rt_view_context_native_measure_candidates(void *ctx, point_t a,
	point_t b, struct rt_view_measure_result *out);

int _rt_view_context_native_screen_to_view(fastf_t *fx, fastf_t *fy,
	void *ctx, fastf_t x, fastf_t y);
int _rt_view_context_native_screen_point_get(point_t point, void *ctx,
	fastf_t x, fastf_t y);
int _rt_view_context_native_current_point_get(point_t point,
	const void *ctx);
int _rt_view_context_native_current_point_set(void *ctx,
	const point_t point);
int _rt_view_context_native_previous_mouse_get(fastf_t *x, fastf_t *y,
	const void *ctx);
int _rt_view_context_native_previous_mouse_set(void *ctx, fastf_t x,
	fastf_t y);
int _rt_view_context_native_mouse_delta_settings_get(
	struct rt_view_mouse_delta_settings *settings, const void *ctx);
int _rt_view_context_native_mouse_state_set(void *ctx, int x, int y);
int _rt_view_context_native_interactive_rect_state_get(
	struct rt_view_interactive_rect_state *record, const void *ctx);
int _rt_view_context_native_interactive_rect_state_set(void *ctx,
	const struct rt_view_interactive_rect_state *record);
int _rt_view_context_native_adc_state_get(struct rt_view_adc_state *record,
	const void *ctx);
int _rt_view_context_native_adc_state_set(void *ctx,
	const struct rt_view_adc_state *record);
int _rt_view_context_native_grid_state_get(struct rt_view_grid_state *record,
	const void *ctx);
int _rt_view_context_native_grid_state_set(void *ctx,
	const struct rt_view_grid_state *record);
int _rt_view_context_native_model_axes_state_get(
	struct rt_view_axes_state *record, const void *ctx);
int _rt_view_context_native_model_axes_state_set(void *ctx,
	const struct rt_view_axes_state *record);
int _rt_view_context_native_view_axes_state_get(
	struct rt_view_axes_state *record, const void *ctx);
int _rt_view_context_native_view_axes_state_set(void *ctx,
	const struct rt_view_axes_state *record);
int _rt_view_context_native_center_dot_state_get(
	struct rt_view_other_state *record, const void *ctx);
int _rt_view_context_native_center_dot_state_set(void *ctx,
	const struct rt_view_other_state *record);
int _rt_view_context_native_scale_overlay_state_get(
	struct rt_view_other_state *record, const void *ctx);
int _rt_view_context_native_scale_overlay_state_set(void *ctx,
	const struct rt_view_other_state *record);
int _rt_view_context_native_params_state_get(
	struct rt_view_params_state *record, const void *ctx);
int _rt_view_context_native_params_state_set(void *ctx,
	const struct rt_view_params_state *record);

int _rt_view_context_native_eye_pos_get(point_t eye_pos, const void *ctx);
int _rt_view_context_native_eye_pos_set(void *ctx, const point_t eye_pos);
int _rt_view_context_native_keypoint_get(point_t keypoint, const void *ctx);
int _rt_view_context_native_keypoint_set(void *ctx, const point_t keypoint);
char _rt_view_context_native_rotate_about_get(const void *ctx);
int _rt_view_context_native_rotate_about_set(void *ctx, char rotate_about);
char _rt_view_context_native_coord_get(const void *ctx);
int _rt_view_context_native_coord_set(void *ctx, char coord);
int _rt_view_context_native_zclip_get(const void *ctx);
int _rt_view_context_native_zclip_set(void *ctx, int zclip);
int _rt_view_context_native_framebuffer_mode_get(const void *ctx);
int _rt_view_context_native_framebuffer_mode_set(void *ctx, int mode);
int _rt_view_context_native_cleared_get(const void *ctx);
int _rt_view_context_native_cleared_set(void *ctx, int cleared);
int _rt_view_context_native_settings_shared(const void *a, const void *b);
int _rt_view_context_native_snap_lines_get(const void *ctx);
int _rt_view_context_native_snap_lines_set(void *ctx, int enabled);
int _rt_view_context_native_snap_source_flags_get(const void *ctx);
int _rt_view_context_native_snap_source_flags_set(void *ctx, int flags);
unsigned long long _rt_view_context_native_snap_kind_mask_get(const void *ctx);
int _rt_view_context_native_snap_exclude_feature_clear(void *ctx);
unsigned long long _rt_view_context_native_prepare_tcl_snap(void *ctx);
double _rt_view_context_native_snap_tolerance_factor_get(const void *ctx);
int _rt_view_context_native_snap_tolerance_factor_set(void *ctx, double factor);

int _rt_view_context_native_knobs_reset(void *ctx, int category);
int _rt_view_context_native_knobs_state_get(struct rt_view_knobs *knobs,
	const void *ctx);
int _rt_view_context_native_knobs_state_set(void *ctx,
	const struct rt_view_knobs *knobs);
unsigned long long _rt_view_context_native_knobs_hash(void *ctx,
	struct bu_data_hash_state *state);
int _rt_view_context_native_knobs_cmd_process(vect_t *rvec, int *do_rot,
	vect_t *tvec, int *do_tran, void *ctx, const char *cmd,
	fastf_t factor, char origin, int model_flag, int incr_flag);
int _rt_view_context_native_knobs_translate(void *ctx, const vect_t tvec,
	int model_flag);
int _rt_view_context_native_knobs_rotate(void *ctx, const vect_t rvec,
	char origin, char coords, const matp_t obj_rot, const pointp_t pvt_pt);
int _rt_view_context_native_knobs_update_rate_flags(void *ctx);
int _rt_view_context_native_knob_values_get(struct rt_view_knob_values *values,
	const void *ctx);
int _rt_view_context_native_knobs_calibrate(void *ctx);

int _rt_view_context_native_lod_policy_get(struct rt_view_lod_policy *policy,
	const void *ctx);
int _rt_view_context_native_lod_policy_apply(void *ctx,
	const struct rt_view_lod_policy *policy);
int _rt_view_context_native_lod_policy_copy(void *dst_ctx,
	const void *src_ctx);
int _rt_view_context_native_lod_bounds_update(void *ctx);
int _rt_view_context_native_lod_bounds_callback_set(void *ctx);
int _rt_view_context_native_lod_bounds_callback_is(const void *ctx);
void *_rt_view_context_native_bounds_update_suspend(void *ctx);
int _rt_view_context_native_bounds_update_restore(void *ctx,
	void *state_ctx, int refresh_bounds);

__END_DECLS

#endif /* LIBRT_VIEW_CONTEXT_PRIVATE_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
