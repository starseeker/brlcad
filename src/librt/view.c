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

#include <stdlib.h>

#include "librt_private.h"
#include "rt/view.h"

void
rt_view_info_init(struct rt_view_info *info)
{
    struct rt_view_info defaults = RT_VIEW_INFO_INIT;
    if (info)
	*info = defaults;
}

void
rt_view_info_sanitize(struct rt_view_info *info)
{
    if (!info)
	return;

    if (info->width <= 0)
	info->width = 1;
    if (info->height <= 0)
	info->height = 1;
    if (info->size <= SMALL_FASTF)
	info->size = 1.0;
    if (info->lod.scale <= SMALL_FASTF)
	info->lod.scale = 1.0;
    if (info->lod.curve_scale <= SMALL_FASTF)
	info->lod.curve_scale = 1.0;
    if (info->lod.point_scale <= SMALL_FASTF)
	info->lod.point_scale = 1.0;
}

void
rt_view_lod_policy_init(struct rt_view_lod_policy *policy)
{
    struct rt_view_lod_policy defaults = RT_VIEW_LOD_POLICY_INIT;
    if (policy)
	*policy = defaults;
}

void
rt_view_lod_policy_sanitize(struct rt_view_lod_policy *policy)
{
    if (!policy)
	return;

    if (policy->scale <= SMALL_FASTF)
	policy->scale = 1.0;
    if (policy->curve_scale <= SMALL_FASTF)
	policy->curve_scale = 1.0;
    if (policy->point_scale <= SMALL_FASTF)
	policy->point_scale = 1.0;
}

static struct rt_view_lod_settings
view_lod_policy(const struct rt_view_info *info)
{
    struct rt_view_lod_settings policy = RT_VIEW_LOD_SETTINGS_INIT;
    if (info)
	policy = info->lod;
    if (policy.curve_scale <= SMALL_FASTF)
	policy.curve_scale = 1.0;
    if (policy.point_scale <= SMALL_FASTF)
	policy.point_scale = 1.0;
    if (policy.scale <= SMALL_FASTF)
	policy.scale = 1.0;
    return policy;
}

fastf_t
rt_view_lod_curve_scale(const struct rt_view_info *v)
{
    return view_lod_policy(v).curve_scale;
}

size_t
rt_view_lod_bot_threshold(const struct rt_view_info *v)
{
    return view_lod_policy(v).bot_threshold;
}

static fastf_t
view_avg_size(const struct rt_view_info *info)
{
    fastf_t view_aspect, x_size, y_size;

    if (!info || info->width <= 0 || info->height <= 0 || info->size <= SMALL_FASTF)
	return 1.0;

    view_aspect = (fastf_t)info->width / info->height;
    x_size = info->size;
    y_size = x_size / view_aspect;

    return (x_size + y_size) / 2.0;
}

fastf_t
rt_view_avg_sample_spacing(const struct rt_view_info *info)
{
    fastf_t avg_view_size, avg_view_samples;

    if (!info || info->width <= 0 || info->height <= 0)
	return 1.0;

    avg_view_size = view_avg_size(info);
    avg_view_samples = (info->width + info->height) / 2.0;

    return avg_view_size / avg_view_samples;
}

fastf_t
rt_view_solid_point_spacing(const struct rt_view_info *info, fastf_t solid_width)
{
    fastf_t radius, avg_view_size, avg_sample_spacing;
    point2d_t p1, p2;

    if (solid_width < SQRT_SMALL_FASTF)
	solid_width = SQRT_SMALL_FASTF;

    avg_view_size = view_avg_size(info);

    /* Now, for the sake of simplicity we're going to make
     * several assumptions:
     *  - our samples represent a grid of square pixels
     *  - a circle with a diameter half the width of the solid is a
     *    good proxy for the kind of curve that will be plotted
     */
    radius = solid_width / 4.0;
    if (avg_view_size < solid_width) {
	/* If the solid is larger than the view, it is
	 * probably only partly visible and likely isn't the
	 * primary focus of the user. We'll cap the point
	 * spacing and avoid wasting effort.
	 */
	radius = avg_view_size / 4.0;
    }

    /* We imagine our representative circular curve lying in
     * the XY plane centered at the origin.
     *
     * Suppose we're viewing the circle head on, and that the
     * apex of the curve (0, radius) lies just inside the
     * top edge of a pixel. Here we place a plotted point p1.
     *
     * As we continue clockwise around the circle we pass
     * through neighboring pixels in the same row, until we
     * vertically drop a distance equal to the pixel spacing,
     * in which case we just barely enter a pixel in the next
     * row. Here we place a plotted point p2 (y = radius -
     * avg_sample_spacing).
     *
     * In theory, the line segment between p1 and p2 passes
     * through all the same pixels that the actual curve does,
     * and thus produces the exact same rasterization as if
     * the curve between p1 and p2 was approximated with an
     * infinite number of line segments.
     *
     * We assume that the distance between p1 and p2 is the
     * maximum point sampling distance we can use for the
     * curve which will give a perfect rasterization, i.e.
     * the same rasterization as if we chose a point distance
     * of 0.
     */
    p1[Y] = radius;
    p1[X] = 0.0;

    avg_sample_spacing = rt_view_avg_sample_spacing(info);
    if (avg_sample_spacing < radius) {
	p2[Y] = radius - (avg_sample_spacing);
    } else {
	/* no particular reason other than symmetry, just need
	 * to prevent sqrt(negative).
	 */
	p2[Y] = radius;
    }
    p2[X] = sqrt((radius * radius) - (p2[Y] * p2[Y]));

    return DIST_PNT2_PNT2(p1, p2) / view_lod_policy(info).point_scale;
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
