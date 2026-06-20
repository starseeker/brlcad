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

RT_EXPORT extern void
rt_view_info_from_bsg(struct rt_view_info *info, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_orientation_quat_from_bsg(quat_t orientation,
				  const struct bsg_view *v);

RT_EXPORT extern int
rt_view_aet_from_bsg(vect_t aet, const struct bsg_view *v);

RT_EXPORT extern fastf_t
rt_view_perspective_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_model2view_from_bsg(mat_t model2view, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_view2model_from_bsg(mat_t view2model, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_rotation_from_bsg(mat_t rotation, const struct bsg_view *v);

RT_EXPORT extern int
rt_view_center_from_bsg(mat_t center, const struct bsg_view *v);

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

RT_EXPORT extern fastf_t
rt_view_scale_from_bsg(const struct bsg_view *v);

RT_EXPORT extern int
rt_view_eye_pos_from_bsg(point_t eye_pos, const struct bsg_view *v);

RT_EXPORT extern int
rt_mesh_lod_load_view_scene_ref_bsg(struct rt_mesh_lod *lod,
				    bsg_scene_ref visibility_ref,
				    struct bsg_view *v,
				    int reset);

RT_EXPORT extern void
rt_mesh_lod_free_scene_ref_bsg(bsg_scene_ref ref);

__END_DECLS

#endif /* RT_VIEW_LEGACY_BSG_H */
