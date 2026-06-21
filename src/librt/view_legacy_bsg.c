/*             V I E W _ L E G A C Y _ B S G . C
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
/** @file rt/view_legacy_bsg.c
 *
 * Transitional BSG view adapters for callers that still own bsg_view state.
 */

#include "common.h"

#include <string.h>

#include "bg/plane.h"
#include "bn/qmath.h"
#include "bsg/faceplate.h"
#include "bsg/lod.h"
#include "bsg/snap.h"
#include "bsg/util.h"
#include "bsg/view_state.h"
#include "rt/view_legacy_bsg.h"

struct bsg_mesh_lod *
_rt_mesh_lod_bsg(struct rt_mesh_lod *lod);

void
rt_view_info_from_bsg(struct rt_view_info *info, const struct bsg_view *v)
{
    if (!info)
	return;

    rt_view_info_init(info);
    if (!v)
	return;

    info->width = v->gv_width;
    info->height = v->gv_height;
    info->size = v->gv_size;

    struct rt_view_lod_policy policy = RT_VIEW_LOD_POLICY_INIT;
    if (rt_view_lod_policy_from_bsg(&policy, v)) {
	info->lod.scale = policy.scale;
	info->lod.curve_scale = policy.curve_scale;
	info->lod.point_scale = policy.point_scale;
	info->lod.bot_threshold = policy.bot_threshold;
    }
    rt_view_info_sanitize(info);
}

int
rt_view_width_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_width : 0;
}

int
rt_view_height_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_height : 0;
}

fastf_t
rt_view_radius_from_bsg(const struct bsg_view *v)
{
    return v ? v->radius : 1.0;
}

int
rt_view_dimensions_set_bsg(struct bsg_view *v, int width, int height)
{
    if (!v)
	return 0;

    v->gv_width = width;
    v->gv_height = height;
    return 1;
}

int
rt_view_screen_to_view_from_bsg(fastf_t *fx,
				fastf_t *fy,
				struct bsg_view *v,
				fastf_t x,
				fastf_t y)
{
    if (fx)
	*fx = 0.0;
    if (fy)
	*fy = 0.0;
    if (!fx || !fy || !v)
	return 0;

    return bsg_screen_to_view(v, fx, fy, x, y) == 0;
}

int
rt_view_screen_point_from_bsg(point_t point,
			      struct bsg_view *v,
			      fastf_t x,
			      fastf_t y)
{
    if (point)
	VSETALL(point, 0.0);
    if (!point || !v)
	return 0;

    point_t model_point = VINIT_ZERO;
    if (bsg_screen_pt(&model_point, x, y, v))
	return 0;

    VMOVE(point, model_point);
    return 1;
}

int
rt_view_interactive_rect_from_bsg(struct bsg_interactive_rect_state *record,
				  const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(struct bsg_interactive_rect_state));
    if (!record || !v)
	return 0;

    return bsg_view_interactive_rect_get(v, record);
}

int
rt_view_interactive_rect_set_bsg(struct bsg_view *v,
				 const struct bsg_interactive_rect_state *record)
{
    if (!v || !record)
	return 0;

    bsg_view_interactive_rect_set(v, record);
    return 1;
}

int
rt_view_adc_from_bsg(struct bsg_adc_state *record,
		     const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(struct bsg_adc_state));
    if (!record || !v)
	return 0;

    return bsg_view_adc_get(v, record);
}

int
rt_view_adc_set_bsg(struct bsg_view *v,
		    const struct bsg_adc_state *record)
{
    if (!v || !record)
	return 0;

    bsg_view_adc_set(v, record);
    return 1;
}

int
rt_view_grid_from_bsg(struct bsg_grid_state *record,
		      const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(struct bsg_grid_state));
    if (!record || !v)
	return 0;

    return bsg_view_grid_get(v, record);
}

int
rt_view_grid_set_bsg(struct bsg_view *v,
		     const struct bsg_grid_state *record)
{
    if (!v || !record)
	return 0;

    bsg_view_grid_set(v, record);
    return 1;
}

int
rt_view_model_axes_from_bsg(struct bsg_axes *record,
			    const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(struct bsg_axes));
    if (!record || !v)
	return 0;

    return bsg_view_model_axes_get(v, record);
}

int
rt_view_model_axes_set_bsg(struct bsg_view *v,
			   const struct bsg_axes *record)
{
    if (!v || !record)
	return 0;

    bsg_view_model_axes_set(v, record);
    return 1;
}

int
rt_view_view_axes_from_bsg(struct bsg_axes *record,
			   const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(struct bsg_axes));
    if (!record || !v)
	return 0;

    return bsg_view_view_axes_get(v, record);
}

int
rt_view_view_axes_set_bsg(struct bsg_view *v,
			  const struct bsg_axes *record)
{
    if (!v || !record)
	return 0;

    bsg_view_view_axes_set(v, record);
    return 1;
}

int
rt_view_center_dot_from_bsg(struct bsg_other_state *record,
			    const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(struct bsg_other_state));
    if (!record || !v)
	return 0;

    return bsg_view_center_dot_get(v, record);
}

int
rt_view_center_dot_set_bsg(struct bsg_view *v,
			   const struct bsg_other_state *record)
{
    if (!v || !record)
	return 0;

    bsg_view_center_dot_set(v, record);
    return 1;
}

int
rt_view_scale_overlay_from_bsg(struct bsg_other_state *record,
			       const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(struct bsg_other_state));
    if (!record || !v)
	return 0;

    return bsg_view_scale_state_get(v, record);
}

int
rt_view_scale_overlay_set_bsg(struct bsg_view *v,
			      const struct bsg_other_state *record)
{
    if (!v || !record)
	return 0;

    bsg_view_scale_state_set(v, record);
    return 1;
}

int
rt_view_params_from_bsg(struct bsg_params_state *record,
			const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(struct bsg_params_state));
    if (!record || !v)
	return 0;

    return bsg_view_params_get(v, record);
}

int
rt_view_params_set_bsg(struct bsg_view *v,
		       const struct bsg_params_state *record)
{
    if (!v || !record)
	return 0;

    bsg_view_params_set(v, record);
    return 1;
}

int
rt_view_refresh_request_bsg(struct bsg_view *v, uint32_t flags)
{
    if (!v)
	return 0;

    bsg_view_refresh_request(v, flags);
    return 1;
}

int
rt_view_refresh_dirty_from_bsg(const struct bsg_view *v)
{
    return bsg_view_refresh_dirty(v);
}

uint32_t
rt_view_refresh_consume_bsg(struct bsg_view *v)
{
    return bsg_view_refresh_consume(v);
}

int
rt_view_refresh_complete_bsg(struct bsg_view *v)
{
    if (!v)
	return 0;

    bsg_view_refresh_complete(v);
    return 1;
}

int
rt_view_refresh_enabled_from_bsg(const struct bsg_view *v)
{
    return bsg_view_refresh_enabled(v);
}

int
rt_view_refresh_enabled_set_bsg(struct bsg_view *v, int enabled)
{
    if (!v)
	return 0;

    bsg_view_refresh_set_enabled(v, enabled);
    return 1;
}

int
rt_view_refresh_suppressed_from_bsg(const struct bsg_view *v)
{
    return bsg_view_refresh_suppressed(v);
}

int
rt_view_refresh_suppress_begin_bsg(struct bsg_view *v)
{
    if (!v)
	return 0;

    bsg_view_refresh_suppress_begin(v);
    return 1;
}

int
rt_view_refresh_suppress_end_bsg(struct bsg_view *v)
{
    if (!v)
	return 0;

    bsg_view_refresh_suppress_end(v);
    return 1;
}

int
rt_view_refresh_drawn_count_from_bsg(const struct bsg_view *v)
{
    return bsg_view_refresh_drawn_count(v);
}

int
rt_view_refresh_drawn_count_set_bsg(struct bsg_view *v, int count)
{
    if (!v)
	return 0;

    bsg_view_refresh_set_drawn_count(v, count);
    return 1;
}

int
rt_view_orientation_quat_from_bsg(quat_t orientation, const struct bsg_view *v)
{
    mat_t identity;

    if (!orientation)
	return 0;

    if (!v) {
	MAT_IDN(identity);
	quat_mat2quat(orientation, identity);
	return 0;
    }

    quat_mat2quat(orientation, v->gv_rotation);
    return 1;
}

int
rt_view_aet_from_bsg(vect_t aet, const struct bsg_view *v)
{
    if (!aet)
	return 0;

    if (!v) {
	VSETALL(aet, 0.0);
	return 0;
    }

    VMOVE(aet, v->gv_aet);
    return 1;
}

int
rt_view_aet_set_bsg(struct bsg_view *v, const vect_t aet)
{
    if (!v || !aet)
	return 0;

    bsg_view_set_aet(v, aet);
    return 1;
}

int
rt_view_aet_state_set_bsg(struct bsg_view *v, const vect_t aet)
{
    if (!v || !aet)
	return 0;

    VMOVE(v->gv_aet, aet);
    return 1;
}

fastf_t
rt_view_perspective_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_perspective : 0.0;
}

int
rt_view_perspective_set_bsg(struct bsg_view *v, fastf_t perspective)
{
    if (!v)
	return 0;

    bsg_view_set_perspective(v, perspective);
    return 1;
}

int
rt_view_model2view_from_bsg(mat_t model2view, const struct bsg_view *v)
{
    if (!model2view)
	return 0;

    if (!v) {
	MAT_IDN(model2view);
	return 0;
    }

    MAT_COPY(model2view, v->gv_model2view);
    return 1;
}

int
rt_view_model2view_set_bsg(struct bsg_view *v, const mat_t model2view)
{
    if (!v || !model2view)
	return 0;

    MAT_COPY(v->gv_model2view, model2view);
    return 1;
}

int
rt_view_view2model_from_bsg(mat_t view2model, const struct bsg_view *v)
{
    if (!view2model)
	return 0;

    if (!v) {
	MAT_IDN(view2model);
	return 0;
    }

    MAT_COPY(view2model, v->gv_view2model);
    return 1;
}

int
rt_view_view2model_set_bsg(struct bsg_view *v, const mat_t view2model)
{
    if (!v || !view2model)
	return 0;

    MAT_COPY(v->gv_view2model, view2model);
    return 1;
}

int
rt_view_pmodel2view_from_bsg(mat_t pmodel2view, const struct bsg_view *v)
{
    if (!pmodel2view)
	return 0;

    if (!v) {
	MAT_IDN(pmodel2view);
	return 0;
    }

    MAT_COPY(pmodel2view, v->gv_pmodel2view);
    return 1;
}

int
rt_view_pmodel2view_set_bsg(struct bsg_view *v, const mat_t pmodel2view)
{
    if (!v || !pmodel2view)
	return 0;

    MAT_COPY(v->gv_pmodel2view, pmodel2view);
    return 1;
}

int
rt_view_pmat_from_bsg(mat_t pmat, const struct bsg_view *v)
{
    if (!pmat)
	return 0;

    if (!v) {
	MAT_IDN(pmat);
	return 0;
    }

    MAT_COPY(pmat, v->gv_pmat);
    return 1;
}

int
rt_view_pmat_set_bsg(struct bsg_view *v, const mat_t pmat)
{
    if (!v || !pmat)
	return 0;

    MAT_COPY(v->gv_pmat, pmat);
    return 1;
}

int
rt_view_rotation_from_bsg(mat_t rotation, const struct bsg_view *v)
{
    if (!rotation)
	return 0;

    if (!v) {
	MAT_IDN(rotation);
	return 0;
    }

    MAT_COPY(rotation, v->gv_rotation);
    return 1;
}

int
rt_view_rotation_set_bsg(struct bsg_view *v, const mat_t rotation)
{
    if (!v || !rotation)
	return 0;

    bsg_view_set_rotation(v, rotation);
    return 1;
}

int
rt_view_center_from_bsg(mat_t center, const struct bsg_view *v)
{
    if (!center)
	return 0;

    if (!v) {
	MAT_IDN(center);
	return 0;
    }

    MAT_COPY(center, v->gv_center);
    return 1;
}

int
rt_view_center_vec_set_bsg(struct bsg_view *v, const point_t center)
{
    if (!v || !center)
	return 0;

    bsg_view_set_center_vec(v, center);
    return 1;
}

int
rt_view_plane_from_bsg(plane_t *plane, const struct bsg_view *v)
{
    point_t center;
    vect_t normal;

    if (!plane || !v)
	return -1;

    MAT_DELTAS_GET_NEG(center, v->gv_center);
    VMOVEN(normal, v->gv_rotation + 8, 3);
    VUNITIZE(normal);
    VSCALE(normal, normal, -1.0);

    return bg_plane_pt_nrml(plane, center, normal);
}

int
rt_view_lod_policy_from_bsg(struct rt_view_lod_policy *policy,
			    const struct bsg_view *v)
{
    if (!policy)
	return 0;

    rt_view_lod_policy_init(policy);
    if (!v)
	return 0;

    struct bsg_lod_source_policy_settings bsg_policy;
    memset(&bsg_policy, 0, sizeof(bsg_policy));
    bsg_policy.policy = BSG_LOD_AUTO;
    bsg_policy.scale = 1.0;
    bsg_policy.curve_scale = 1.0;
    bsg_policy.point_scale = 1.0;
    if (!bsg_view_lod_source_policy_get(v, &bsg_policy))
	return 0;

    policy->policy = (int)bsg_policy.policy;
    policy->forced_level = bsg_policy.forced_level;
    policy->mesh_enabled = bsg_policy.mesh_enabled ? 1 : 0;
    policy->csg_enabled = bsg_policy.csg_enabled ? 1 : 0;
    policy->zoom_refresh = bsg_policy.zoom_refresh ? 1 : 0;
    policy->bot_threshold = bsg_policy.bot_threshold;
    policy->scale = bsg_policy.scale;
    policy->curve_scale = bsg_policy.curve_scale;
    policy->point_scale = bsg_policy.point_scale;
    rt_view_lod_policy_sanitize(policy);
    return 1;
}

int
rt_view_lod_policy_apply_bsg(struct bsg_view *v,
			     const struct rt_view_lod_policy *policy)
{
    if (!v || !policy)
	return 0;

    struct rt_view_lod_policy sanitized = *policy;
    rt_view_lod_policy_sanitize(&sanitized);

    struct bsg_lod_source_policy_settings bsg_policy;
    memset(&bsg_policy, 0, sizeof(bsg_policy));
    bsg_policy.policy = (bsg_lod_policy)sanitized.policy;
    bsg_policy.forced_level = sanitized.forced_level;
    bsg_policy.mesh_enabled = sanitized.mesh_enabled ? 1 : 0;
    bsg_policy.csg_enabled = sanitized.csg_enabled ? 1 : 0;
    bsg_policy.zoom_refresh = sanitized.zoom_refresh ? 1 : 0;
    bsg_policy.scale = sanitized.scale;
    bsg_policy.bot_threshold = sanitized.bot_threshold;
    bsg_policy.curve_scale = sanitized.curve_scale;
    bsg_policy.point_scale = sanitized.point_scale;
    bsg_view_lod_source_policy_set(v, &bsg_policy);
    return 1;
}

int
rt_view_lod_policy_copy_bsg(struct bsg_view *dst, const struct bsg_view *src)
{
    struct rt_view_lod_policy policy = RT_VIEW_LOD_POLICY_INIT;
    if (!dst || !rt_view_lod_policy_from_bsg(&policy, src))
	return 0;
    return rt_view_lod_policy_apply_bsg(dst, &policy);
}

void
rt_view_lod_bounds_update_bsg(struct bsg_view *v)
{
    if (v)
	bsg_view_bounds(v);
}

void
rt_view_lod_bounds_callback_set_bsg(struct bsg_view *v)
{
    if (v)
	v->gv_bounds_update = &bsg_view_bounds;
}

int
rt_view_lod_bounds_callback_is_bsg(const struct bsg_view *v)
{
    return (v && v->gv_bounds_update == &bsg_view_bounds) ? 1 : 0;
}

rt_view_bounds_update_callback_bsg_t
rt_view_bounds_update_callback_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_bounds_update : NULL;
}

int
rt_view_bounds_update_callback_set_bsg(struct bsg_view *v,
				       rt_view_bounds_update_callback_bsg_t callback)
{
    if (!v)
	return 0;

    v->gv_bounds_update = callback;
    return 1;
}

int
rt_view_bounds_update_callback_call_bsg(struct bsg_view *v)
{
    rt_view_bounds_update_callback_bsg_t callback =
	rt_view_bounds_update_callback_from_bsg(v);

    if (!callback)
	return 0;

    (*callback)(v);
    return 1;
}

fastf_t
rt_view_scale_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_scale : 1.0;
}

int
rt_view_scale_set_bsg(struct bsg_view *v, fastf_t scale)
{
    if (!v)
	return 0;

    bsg_view_set_scale(v, scale);
    return 1;
}

fastf_t *
rt_view_scale_storage_from_bsg(struct bsg_view *v)
{
    return v ? &v->gv_scale : NULL;
}

int
rt_view_scale_state_set_bsg(struct bsg_view *v,
			    fastf_t scale,
			    fastf_t initial_scale,
			    fastf_t absolute_scale,
			    fastf_t size,
			    fastf_t inverse_size)
{
    if (!v)
	return 0;

    v->gv_scale = scale;
    v->gv_i_scale = initial_scale;
    v->gv_a_scale = absolute_scale;
    v->gv_size = size;
    v->gv_isize = inverse_size;
    return 1;
}

fastf_t
rt_view_initial_scale_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_i_scale : 1.0;
}

int
rt_view_initial_scale_set_bsg(struct bsg_view *v, fastf_t scale)
{
    if (!v)
	return 0;

    v->gv_i_scale = scale;
    return 1;
}

fastf_t
rt_view_absolute_scale_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_a_scale : 0.0;
}

int
rt_view_absolute_scale_set_bsg(struct bsg_view *v, fastf_t scale)
{
    if (!v)
	return 0;

    v->gv_a_scale = scale;
    return 1;
}

fastf_t
rt_view_local2base_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_local2base : 1.0;
}

fastf_t
rt_view_base2local_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_base2local : 1.0;
}

int
rt_view_unit_conversion_set_bsg(struct bsg_view *v,
				fastf_t local2base,
				fastf_t base2local)
{
    if (!v)
	return 0;

    v->gv_local2base = local2base;
    v->gv_base2local = base2local;
    return 1;
}

fastf_t
rt_view_size_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_size : 1.0;
}

int
rt_view_size_set_bsg(struct bsg_view *v, fastf_t size)
{
    if (!v)
	return 0;

    bsg_view_set_size(v, size);
    return 1;
}

fastf_t
rt_view_inverse_size_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_isize : 1.0;
}

int
rt_view_eye_pos_from_bsg(point_t eye_pos, const struct bsg_view *v)
{
    if (!eye_pos)
	return 0;

    if (!v) {
	VSET(eye_pos, 0.0, 0.0, 1.0);
	return 0;
    }

    VMOVE(eye_pos, v->gv_eye_pos);
    return 1;
}

int
rt_view_eye_pos_set_bsg(struct bsg_view *v, const point_t eye_pos)
{
    if (!v || !eye_pos)
	return 0;

    VMOVE(v->gv_eye_pos, eye_pos);
    return 1;
}

int
rt_view_keypoint_from_bsg(point_t keypoint, const struct bsg_view *v)
{
    if (!keypoint)
	return 0;

    if (!v) {
	VSETALL(keypoint, 0.0);
	return 0;
    }

    VMOVE(keypoint, v->gv_keypoint);
    return 1;
}

int
rt_view_keypoint_set_bsg(struct bsg_view *v, const point_t keypoint)
{
    if (!v || !keypoint)
	return 0;

    VMOVE(v->gv_keypoint, keypoint);
    return 1;
}

char
rt_view_rotate_about_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_rotate_about : 'v';
}

int
rt_view_rotate_about_set_bsg(struct bsg_view *v, char rotate_about)
{
    if (!v)
	return 0;

    v->gv_rotate_about = rotate_about;
    return 1;
}

char
rt_view_coord_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_coord : 'v';
}

int
rt_view_coord_set_bsg(struct bsg_view *v, char coord)
{
    if (!v)
	return 0;

    v->gv_coord = coord;
    return 1;
}

int
rt_view_snap_lines_from_bsg(const struct bsg_view *v)
{
    return bsg_view_snap_lines(v);
}

int
rt_view_snap_lines_set_bsg(struct bsg_view *v, int enabled)
{
    if (!v)
	return 0;

    bsg_view_set_snap_lines(v, enabled);
    return 1;
}

int
rt_view_snap_source_flags_from_bsg(const struct bsg_view *v)
{
    return bsg_view_snap_source_flags(v);
}

int
rt_view_snap_source_flags_set_bsg(struct bsg_view *v, int flags)
{
    if (!v)
	return 0;

    bsg_view_set_snap_source_flags(v, flags);
    return 1;
}

unsigned long long
rt_view_snap_kind_mask_from_bsg(const struct bsg_view *v)
{
    return (unsigned long long)bsg_view_snap_kind_mask(v);
}

unsigned long long
rt_view_prepare_tcl_snap_bsg(struct bsg_view *v)
{
    return (unsigned long long)bsg_view_prepare_tcl_snap(v);
}

int
rt_view_center_linesnap_bsg(struct bsg_view *v)
{
    if (!v)
	return 0;

    bsg_view_center_linesnap(v);
    return 1;
}

int
rt_view_zclip_from_bsg(const struct bsg_view *v)
{
    return bsg_view_zclip(v);
}

int
rt_view_zclip_set_bsg(struct bsg_view *v, int zclip)
{
    if (!v)
	return 0;

    bsg_view_set_zclip(v, zclip);
    return 1;
}

int
rt_view_framebuffer_mode_from_bsg(const struct bsg_view *v)
{
    return bsg_view_framebuffer_mode(v);
}

int
rt_view_framebuffer_mode_set_bsg(struct bsg_view *v, int mode)
{
    if (!v)
	return 0;

    bsg_view_set_framebuffer_mode(v, mode);
    return 1;
}

int
rt_view_cleared_from_bsg(const struct bsg_view *v)
{
    return bsg_view_cleared(v);
}

int
rt_view_cleared_set_bsg(struct bsg_view *v, int cleared)
{
    if (!v)
	return 0;

    bsg_view_set_cleared(v, cleared);
    return 1;
}

int
rt_view_settings_shared_bsg(const struct bsg_view *a, const struct bsg_view *b)
{
    return bsg_view_settings_shared(a, b);
}

double
rt_view_snap_tolerance_factor_from_bsg(const struct bsg_view *v)
{
    return bsg_view_snap_tolerance_factor(v);
}

int
rt_view_snap_tolerance_factor_set_bsg(struct bsg_view *v, double factor)
{
    if (!v)
	return 0;

    bsg_view_set_snap_tolerance_factor(v, factor);
    return 1;
}

int
rt_mesh_lod_load_view_scene_ref_bsg(struct rt_mesh_lod *lod,
				    bsg_scene_ref visibility_ref,
				    struct bsg_view *v,
				    int reset)
{
    struct bsg_mesh_lod *bsg_lod = _rt_mesh_lod_bsg(lod);
    if (!bsg_lod)
	return -1;

    return bsg_mesh_lod_load_view_scene_ref(bsg_lod, visibility_ref, v, reset);
}

void
rt_mesh_lod_free_scene_ref_bsg(bsg_scene_ref ref)
{
    bsg_mesh_lod_free_scene_ref(ref);
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
