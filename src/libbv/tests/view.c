/*                   L I B B V / V I E W . C
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

#include "common.h"

#include <stdio.h>

#include "bv.h"

static int
near_fastf(fastf_t a, fastf_t b)
{
    return NEAR_EQUAL(a, b, 1.0e-6);
}

static int
near_point(const point_t a, fastf_t x, fastf_t y, fastf_t z)
{
    return near_fastf(a[X], x) && near_fastf(a[Y], y) && near_fastf(a[Z], z);
}

static int
fail(const char *msg)
{
    fprintf(stderr, "%s\n", msg);
    return 1;
}

int
main(int UNUSED(argc), const char **UNUSED(argv))
{
    struct bv *v = bv_create();
    struct bv *v2 = bv_create();
    struct bv_set *set = bv_set_create();
    point_t min;
    point_t max;
    point_t center;
    point_t model;
    vect_t rvec = VINIT_ZERO;
    vect_t tvec = VINIT_ZERO;
    fastf_t sx = 1.0;
    fastf_t sy = 1.0;
    void *sentinel = (void *)v;
    int do_rot = 0;
    int do_tran = 0;

    if (!bv_is_valid(v))
	return fail("new view is not valid");
    if (!near_fastf(v->scale, BV_DEFAULT_SCALE) || !near_fastf(v->size, 2.0 * BV_DEFAULT_SCALE))
	return fail("new view scale/size defaults are incorrect");
    if (!bv_refresh_enabled_get(v) || bv_refresh_dirty_get(v))
	return fail("new view refresh defaults are incorrect");

    if (!bv_dimensions_set(v, 800, 400))
	return fail("failed to set view dimensions");
    if (!bv_user_data_set(v, sentinel) || bv_user_data_get(v) != sentinel)
	return fail("failed to set view user data");
    if (!bv_screen_to_view(&sx, &sy, v, 400.0, 200.0))
	return fail("failed to convert screen center to view");
    if (!near_fastf(sx, 0.0) || !near_fastf(sy, 0.0))
	return fail("screen center did not map to view origin");

    VSET(min, -1.0, -2.0, -3.0);
    VSET(max, 3.0, 2.0, 1.0);
    if (!bv_autoview_bounds(v, BV_AUTOVIEW_SCALE_DEFAULT, min, max))
	return fail("autoview bounds failed");
    if (!bv_center_get(center, v) || !near_point(center, 1.0, 0.0, -1.0))
	return fail("autoview center is incorrect");
    if (!near_fastf(v->scale, 2.0) || !near_fastf(v->size, 4.0))
	return fail("autoview scale/size are incorrect");
    if (!bv_screen_to_model(model, v, 400.0, 200.0) || !near_point(model, 1.0, 0.0, -1.0))
	return fail("screen center did not map to model center");

    if (bv_knobs_cmd_process(&rvec, &do_rot, &tvec, &do_tran, v, "aX", 10.0, 'v', 0, 0) != BRLCAD_OK)
	return fail("absolute X knob command failed");
    if (!do_tran || do_rot || !near_fastf(tvec[X], 10.0))
	return fail("absolute X knob command did not report the expected translation");

    if (!bv_refresh_request(v, 0x3) || !bv_refresh_dirty_get(v))
	return fail("refresh request did not mark view dirty");
    if (bv_refresh_consume(v) != 0x3 || bv_refresh_dirty_get(v))
	return fail("refresh consume did not return and clear dirty flags");
    if (!bv_refresh_suppress_begin(v) || !bv_refresh_request(v, 0x4) ||
	    bv_refresh_dirty_get(v) || !bv_refresh_suppressed_get(v))
	return fail("refresh suppression did not hold dirty flags");
    if (!bv_refresh_suppress_end(v) || bv_refresh_suppressed_get(v))
	return fail("refresh suppression did not clear");
    if (!bv_refresh_enabled_set(v, 0) || bv_refresh_enabled_get(v) ||
	    !bv_refresh_request(v, 0x8) || bv_refresh_dirty_get(v))
	return fail("disabled refresh did not hold dirty flags");
    if (!bv_refresh_enabled_set(v, 1) || !bv_refresh_drawn_count_set(v, 7) ||
	    bv_refresh_drawn_count_get(v) != 7)
	return fail("refresh drawn count state failed");
    if (!bv_frametime_set(v, 1234) || bv_frametime_get(v) != 1234 ||
	    !bv_zclip_set(v, 1) || !bv_zclip_get(v) ||
	    !bv_framebuffer_mode_set(v, 2) || bv_framebuffer_mode_get(v) != 2 ||
	    !bv_cleared_set(v, 1) || !bv_cleared_get(v))
	return fail("view bookkeeping state failed");

    if (!bv_name_set(v, "primary") || !bv_set_add(set, v))
	return fail("failed to add named view to set");
    if (bv_set_find(set, "primary") != v)
	return fail("failed to find named view in set");
    if (!bv_set_remove(set, v) || bv_set_find(set, "primary") != NULL)
	return fail("failed to remove view from set");

    if (!bv_copy(v2, v) || !bv_center_get(center, v2) || !near_point(center, 1.0, 0.0, -1.0))
	return fail("view copy did not preserve model center");
    if (bv_user_data_get(v2) != sentinel ||
	    bv_refresh_drawn_count_get(v2) != 7 ||
	    bv_frametime_get(v2) != 1234 || !bv_zclip_get(v2) ||
	    bv_framebuffer_mode_get(v2) != 2 ||
	    !bv_cleared_get(v2))
	return fail("view copy did not preserve bookkeeping state");

    bv_set_destroy(set);
    bv_destroy(v2);
    bv_destroy(v);

    return 0;
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
