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

#include "bn/qmath.h"
#include "bsg/lod.h"
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

fastf_t
rt_view_perspective_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_perspective : 0.0;
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

fastf_t
rt_view_scale_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_scale : 1.0;
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
