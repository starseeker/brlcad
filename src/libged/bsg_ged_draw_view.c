/*                B S G _ G E D _ D R A W _ V I E W . C
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this file; see the file named COPYING for more information.
 */
/** @file bsg_ged_draw_view.c
 *
 * BSG-to-RT view snapshot adapters for the transitional GED draw bridge.
 */

#include "common.h"

#include "rt/view_legacy_bsg.h"
#include "./bsg_ged_draw_private.h"

void
ged_draw_view_info_from_bsg(struct rt_view_info *view_info,
			    const struct bsg_view *view)
{
    rt_view_info_from_bsg(view_info, view);
}

fastf_t
ged_draw_view_perspective_from_bsg(const struct bsg_view *view)
{
    return rt_view_perspective_from_bsg(view);
}

fastf_t
ged_draw_view_scale_from_bsg(const struct bsg_view *view)
{
    return rt_view_scale_from_bsg(view);
}

int
ged_draw_view_lod_policy_from_bsg(ged_draw_view_lod_policy *policy,
				  const struct bsg_view *view)
{
    return rt_view_lod_policy_from_bsg(policy, view);
}

int
ged_draw_view_lod_policy_apply_bsg(struct bsg_view *view,
				   const ged_draw_view_lod_policy *policy)
{
    return rt_view_lod_policy_apply_bsg(view, policy);
}

int
ged_draw_view_lod_policy_apply_bsg_bot_threshold(
	struct bsg_view *view,
	const ged_draw_view_lod_policy *policy,
	size_t bot_threshold)
{
    if (!policy)
	return 0;

    ged_draw_view_lod_policy override_policy = *policy;
    override_policy.bot_threshold = bot_threshold;
    return ged_draw_view_lod_policy_apply_bsg(view, &override_policy);
}

void
ged_draw_view_set_lod_bounds_update(struct bsg_view *view)
{
    rt_view_lod_bounds_callback_set_bsg(view);
}

int
ged_draw_view_has_lod_bounds_update(const struct bsg_view *view)
{
    return rt_view_lod_bounds_callback_is_bsg(view);
}

void
ged_draw_scene_ref_realization_set_bsg_view_policy(bsg_scene_ref ref,
						   const struct bsg_view *view)
{
    if (!view)
	return;

    ged_draw_view_lod_policy policy;
    ged_draw_view_lod_policy_from_bsg(&policy, view);
    ged_draw_scene_ref_realization_set_view_policy(ref,
	    policy.csg_enabled,
	    rt_view_scale_from_bsg(view),
	    policy.bot_threshold,
	    policy.curve_scale,
	    policy.point_scale);
}

int
ged_draw_mesh_lod_load_view_scene_ref(struct rt_mesh_lod *lod,
				      bsg_scene_ref visibility_ref,
				      struct bsg_view *view,
				      int reset)
{
    return rt_mesh_lod_load_view_scene_ref_bsg(lod, visibility_ref, view, reset);
}

void
ged_draw_mesh_lod_free_scene_ref(bsg_scene_ref ref)
{
    rt_mesh_lod_free_scene_ref_bsg(ref);
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
