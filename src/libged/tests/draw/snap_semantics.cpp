/*              S N A P _ P A R I T Y . C P P
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
/** @file snap_semantics.cpp
 *
 * Verify libbv grid snapping for a headless passive view.  Scene/object snap
 * behavior is covered by Obol/qtcad and GED draw tests; this test keeps the
 * display-independent grid math covered without exercising the retired RT
 * view-context snap facade.
 */

#include "common.h"

#include <cstdio>
#include <cmath>

#include <bu.h>
#include "view_test_util.h"

#define ASSERT(cond) do { \
    nchecks++; \
    if (!(cond)) { \
	bu_log("FAIL [%s:%d] %s\n", __FILE__, __LINE__, #cond); \
	nfails++; \
    } \
} while (0)

static int nchecks = 0;
static int nfails  = 0;


int
main(int UNUSED(ac), char *av[])
{
    bu_setprogname(av[0]);

    struct bv_context *ctx = bv_context_create();
    ASSERT(ctx != NULL);
    struct bv *view = bv_context_view(ctx);
    ASSERT(view != NULL);
    ASSERT(bv_dimensions_set(view, 512, 512));

    struct bv_grid_state grid = BV_GRID_STATE_INIT;
    ASSERT(bv_grid_state_get(&grid, view));
    grid.res_h = 1.0;
    grid.res_v = 1.0;
    ASSERT(bv_grid_state_set(view, &grid));

    /* ------------------------------------------------------------------
     * Test 1: snap view-space coordinates to the projected model grid.
     * ------------------------------------------------------------------ */
    {
	fastf_t vx = 0.401;
	fastf_t vy = -0.317;
	const fastf_t start_x = vx;
	const fastf_t start_y = vy;
	ASSERT(bv_snap_grid_2d(view, &vx, &vy));
	ASSERT(std::isfinite((double)vx));
	ASSERT(std::isfinite((double)vy));
	ASSERT(fabs(vx - start_x) <= 0.5 / bv_scale_get(view) + VUNITIZE_TOL);
	ASSERT(fabs(vy - start_y) <= 0.5 / bv_scale_get(view) + VUNITIZE_TOL);

	fastf_t qx = vx * bv_scale_get(view) * bv_base2local_get(view);
	fastf_t qy = vy * bv_scale_get(view) * bv_base2local_get(view);
	ASSERT(fabs(qx - round(qx)) < 1.0e-9);
	ASSERT(fabs(qy - round(qy)) < 1.0e-9);
    }

    /* ------------------------------------------------------------------
     * Test 2: origin snap stays finite and stable.
     * ------------------------------------------------------------------ */
    {
	fastf_t vx = 0.0, vy = 0.0;
	ASSERT(bv_snap_grid_2d(view, &vx, &vy));
	ASSERT(std::isfinite((double)vx));
	ASSERT(std::isfinite((double)vy));
	ASSERT(fabs(vx) < 1.0e-12);
	ASSERT(fabs(vy) < 1.0e-12);
    }

    /* ------------------------------------------------------------------
     * Test 3: NULL guards — must not crash.
     * ------------------------------------------------------------------ */
    {
	fastf_t vx = 0.0, vy = 0.0;
	ASSERT(bv_snap_grid_2d(NULL, &vx, &vy) == 0);
	ASSERT(bv_snap_grid_2d(view, NULL, &vy) == 0);
	ASSERT(bv_snap_grid_2d(view, &vx, NULL) == 0);
    }

    bv_context_destroy(ctx);

    bu_log("snap semantic checks: %d checks, %d failures\n", nchecks, nfails);
    return nfails ? 1 : 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
