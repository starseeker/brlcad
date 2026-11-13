/*                     V L I S T . C
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

#include <math.h>

#include "vmath.h"
#include "bg/vlist.h"
#include "bu/app.h"
#include "bu/log.h"


static int
expect_point(const char *name, const point_t actual, fastf_t x, fastf_t y, fastf_t z)
{
    if (!NEAR_EQUAL(actual[X], x, SMALL_FASTF) ||
	!NEAR_EQUAL(actual[Y], y, SMALL_FASTF) ||
	!NEAR_EQUAL(actual[Z], z, SMALL_FASTF)) {
	bu_log("FAIL: %s is (%g, %g, %g), expected (%g, %g, %g)\n",
	       name, V3ARGS(actual), x, y, z);
	return 1;
    }

    return 0;
}


int
main(int UNUSED(argc), char *argv[])
{
    struct bu_list head;
    struct bu_list copied;
    struct bu_list vlfree;
    point_t bmin;
    point_t bmax;
    point_t pts[8] = {
	{0.0, 0.0, 0.0},
	{1.0, 0.0, 0.0},
	{1.0, 1.0, 0.0},
	{0.0, 1.0, 0.0},
	{0.0, 0.0, 1.0},
	{1.0, 0.0, 1.0},
	{1.0, 1.0, 1.0},
	{0.0, 1.0, 1.0}
    };
    size_t length = 0;
    int dispmode = 0;
    int ret = 0;

    bu_setprogname(argv[0]);

    BU_LIST_INIT(&head);
    BU_LIST_INIT(&copied);
    BU_LIST_INIT(&vlfree);

    bg_vlist_arb8(&head, &vlfree, pts);
    if (bg_vlist_cmd_cnt((bg_vlist *)&head) != 18) {
	bu_log("FAIL: bg_vlist_arb8 did not emit 18 commands\n");
	ret = 1;
    }

    VSET(bmin, INFINITY, INFINITY, INFINITY);
    VSET(bmax, -INFINITY, -INFINITY, -INFINITY);
    if (bg_vlist_bbox(&head, &bmin, &bmax, &length, &dispmode) != 0) {
	bu_log("FAIL: bg_vlist_bbox rejected bg_vlist_arb8 output\n");
	ret = 1;
    }
    if (length != 18) {
	bu_log("FAIL: bg_vlist_bbox length is %zu, expected 18\n", length);
	ret = 1;
    }
    if (dispmode != 0) {
	bu_log("FAIL: bg_vlist_bbox reported display-mode data for simple wire output\n");
	ret = 1;
    }
    ret |= expect_point("bmin", bmin, 0.0, 0.0, 0.0);
    ret |= expect_point("bmax", bmax, 1.0, 1.0, 1.0);

    bg_vlist_copy(&vlfree, &copied, &head);
    if (bg_vlist_cmd_cnt((bg_vlist *)&copied) != 18) {
	bu_log("FAIL: bg_vlist_copy did not preserve command count\n");
	ret = 1;
    }

    BG_FREE_VLIST(&vlfree, &head);
    BG_FREE_VLIST(&vlfree, &copied);

    struct bg_vlblock *vbp = bg_vlblock_init(&vlfree, 2);
    struct bu_list *red = bg_vlblock_find(vbp, 255, 0, 0);
    if (!red) {
	bu_log("FAIL: bg_vlblock_find returned NULL\n");
	ret = 1;
    } else {
	BG_ADD_VLIST(&vlfree, red, pts[0], BG_VLIST_LINE_MOVE);
    }
    bg_vlblock_free(vbp);
    bg_vlist_cleanup(&vlfree);

    return ret;
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
