/*                       V I E W _ I N F O . C
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
/** @file librt/tests/view_info.c
 *
 * Unit tests for RT-owned view snapshots.
 */

#include "common.h"

#include <math.h>
#include <stdio.h>

#include "bu/app.h"
#include "rt/view.h"

static int
fastf_equal(fastf_t a, fastf_t b)
{
    return fabs(a - b) <= SMALL_FASTF;
}

static int
fastf_near(fastf_t a, fastf_t b, fastf_t tol)
{
    return fabs(a - b) <= tol;
}

static int
check_view_info(const char *label,
		const struct rt_view_info *info,
		int width,
		int height,
		fastf_t size,
		fastf_t scale,
		fastf_t curve_scale,
		fastf_t point_scale,
		size_t bot_threshold)
{
    if (!info) {
	printf("FAIL: %s info is NULL\n", label);
	return 1;
    }

    if (info->width != width || info->height != height ||
	    !fastf_equal(info->size, size)) {
	printf("FAIL: %s dimensions/size: got %dx%d size %g\n",
		label, info->width, info->height, info->size);
	return 1;
    }

    if (!fastf_equal(info->lod.scale, scale) ||
	    !fastf_equal(info->lod.curve_scale, curve_scale) ||
	    !fastf_equal(info->lod.point_scale, point_scale) ||
	    info->lod.bot_threshold != bot_threshold) {
	printf("FAIL: %s lod: got scale=%g curve=%g point=%g threshold=%zu\n",
		label, info->lod.scale, info->lod.curve_scale,
		info->lod.point_scale, info->lod.bot_threshold);
	return 1;
    }

    return 0;
}

static int
test_defaults(void)
{
    struct rt_view_info info;

    info.width = 23;
    info.height = 42;
    info.size = 99.0;
    info.lod.scale = 7.0;
    info.lod.curve_scale = 8.0;
    info.lod.point_scale = 9.0;
    info.lod.bot_threshold = 100;

    rt_view_info_init(NULL);
    rt_view_info_init(&info);

    return check_view_info("defaults", &info, 1, 1, 1.0, 1.0, 1.0, 1.0, 0);
}

static int
test_sanitize(void)
{
    struct rt_view_info info;

    info.width = 0;
    info.height = -10;
    info.size = 0.0;
    info.lod.scale = 0.0;
    info.lod.curve_scale = -2.0;
    info.lod.point_scale = 0.0;
    info.lod.bot_threshold = 33;

    rt_view_info_sanitize(NULL);
    rt_view_info_sanitize(&info);

    return check_view_info("sanitize", &info, 1, 1, 1.0, 1.0, 1.0, 1.0, 33);
}

static int
test_lod_policy_defaults(void)
{
    struct rt_view_lod_policy policy;

    policy.policy = RT_VIEW_LOD_OFF;
    policy.forced_level = 9;
    policy.mesh_enabled = 1;
    policy.csg_enabled = 1;
    policy.zoom_refresh = 1;
    policy.bot_threshold = 111;
    policy.scale = 7.0;
    policy.curve_scale = 8.0;
    policy.point_scale = 9.0;

    rt_view_lod_policy_init(NULL);
    rt_view_lod_policy_init(&policy);
    if (policy.policy != RT_VIEW_LOD_AUTO ||
	    policy.forced_level != 0 ||
	    policy.mesh_enabled != 0 ||
	    policy.csg_enabled != 0 ||
	    policy.zoom_refresh != 0 ||
	    policy.bot_threshold != 0 ||
	    !fastf_equal(policy.scale, 1.0) ||
	    !fastf_equal(policy.curve_scale, 1.0) ||
	    !fastf_equal(policy.point_scale, 1.0)) {
	printf("FAIL: lod policy defaults\n");
	return 1;
    }

    policy.scale = 0.0;
    policy.curve_scale = -2.0;
    policy.point_scale = 0.0;
    policy.bot_threshold = 222;
    rt_view_lod_policy_sanitize(NULL);
    rt_view_lod_policy_sanitize(&policy);
    if (!fastf_equal(policy.scale, 1.0) ||
	    !fastf_equal(policy.curve_scale, 1.0) ||
	    !fastf_equal(policy.point_scale, 1.0) ||
	    policy.bot_threshold != 222) {
	printf("FAIL: lod policy sanitize\n");
	return 1;
    }

    return 0;
}

static int
test_view_metrics(void)
{
    struct rt_view_info info = RT_VIEW_INFO_INIT;
    struct rt_view_info bad_policy = RT_VIEW_INFO_INIT;
    fastf_t point_spacing;

    info.width = 200;
    info.height = 100;
    info.size = 20.0;
    info.lod.curve_scale = 3.0;
    info.lod.point_scale = 2.0;
    info.lod.bot_threshold = 42;

    if (!fastf_equal(rt_view_lod_curve_scale(&info), 3.0) ||
	    rt_view_lod_bot_threshold(&info) != 42) {
	printf("FAIL: view metric lod policy\n");
	return 1;
    }

    if (!fastf_equal(rt_view_avg_sample_spacing(&info), 0.1)) {
	printf("FAIL: view average sample spacing: got %g\n",
		rt_view_avg_sample_spacing(&info));
	return 1;
    }

    point_spacing = rt_view_solid_point_spacing(&info, 4.0);
    if (!fastf_near(point_spacing, sqrt(0.2) / 2.0, 1.0e-6)) {
	printf("FAIL: view solid point spacing: got %g\n", point_spacing);
	return 1;
    }

    bad_policy.lod.curve_scale = 0.0;
    bad_policy.lod.point_scale = -2.0;
    bad_policy.lod.scale = 0.0;
    bad_policy.lod.bot_threshold = 55;
    if (!fastf_equal(rt_view_lod_curve_scale(&bad_policy), 1.0) ||
	    rt_view_lod_bot_threshold(&bad_policy) != 55) {
	printf("FAIL: view metric bad lod policy\n");
	return 1;
    }

    if (!fastf_equal(rt_view_lod_curve_scale(NULL), 1.0) ||
	    rt_view_lod_bot_threshold(NULL) != 0 ||
	    !fastf_equal(rt_view_avg_sample_spacing(NULL), 1.0)) {
	printf("FAIL: view metric null defaults\n");
	return 1;
    }

    return 0;
}

int
main(int argc, char *argv[])
{
    bu_setprogname(argv[0]);

    if (argc != 1) {
	printf("Usage: %s\n", argv[0]);
	return 1;
    }

    if (test_defaults())
	return 1;
    if (test_sanitize())
	return 1;
    if (test_lod_policy_defaults())
	return 1;
    if (test_view_metrics())
	return 1;

    return 0;
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
