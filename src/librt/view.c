/*                  V I E W . C
 * BRL-CAD
 *
 * Copyright (c) 2008-2026 United States Government as represented by
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
/** @file rt/view.c
 *
 * Information routines to support LoD primitive realization.
 *
 */

#include "common.h"

#include <math.h>
#include <stdlib.h>

#include "librt_private.h"
#include "rt/view.h"

void
rt_view_info_init(struct rt_view_info *info)
{
    bv_view_info_init(info);
}

void
rt_view_info_sanitize(struct rt_view_info *info)
{
    bv_view_info_sanitize(info);
}

void
rt_view_lod_policy_init(struct rt_view_lod_policy *policy)
{
    bv_lod_policy_init(policy);
}

void
rt_view_lod_policy_sanitize(struct rt_view_lod_policy *policy)
{
    bv_lod_policy_sanitize(policy);
}

fastf_t
rt_view_lod_curve_scale(const struct rt_view_info *v)
{
    return bv_view_lod_curve_scale(v);
}

size_t
rt_view_lod_bot_threshold(const struct rt_view_info *v)
{
    return bv_view_lod_bot_threshold(v);
}

fastf_t
rt_view_avg_sample_spacing(const struct rt_view_info *info)
{
    return bv_view_avg_sample_spacing(info);
}

fastf_t
rt_view_solid_point_spacing(const struct rt_view_info *info, fastf_t solid_width)
{
    return bv_view_solid_point_spacing(info, solid_width);
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
