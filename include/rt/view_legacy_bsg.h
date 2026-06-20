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
#include "rt/defines.h"
#include "rt/view.h"

__BEGIN_DECLS

struct bsg_view;
struct rt_mesh_lod;

typedef void (*rt_view_bounds_update_callback_bsg_t)(struct bsg_view *);

RT_EXPORT extern void
rt_view_info_from_bsg(struct rt_view_info *info, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_width_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_height_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_dimensions_set_bsg(struct bsg_view *v, int width, int height);

RT_EXPORT extern int
rt_view_orientation_quat_from_bsg(quat_t orientation,
				  const struct bsg_view *v);

RT_EXPORT extern int
rt_view_aet_from_bsg(vect_t aet, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_aet_set_bsg(struct bsg_view *v, const vect_t aet);

RT_EXPORT extern int
rt_view_aet_state_set_bsg(struct bsg_view *v, const vect_t aet);

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
