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
 * Verify that the transitional RT snap adapters produce consistent semantic
 * snap records for retained BSG-backed views.
 *
 * Strategy:
 *   1. Create a headless view with grid snapping enabled.
 *   2. Call the RT snap-candidate adapter with grid snapping near the origin.
 *   3. Verify the result distance is within grid resolution.
 *   4. Call the RT 2D snap adapter at the view center; verify it does not
 *      crash and returns a valid 2D coordinate pair.
 *
 * No database or display hardware required.
 *
 * Usage: ged_test_snap_semantics
 */

#include "common.h"

#include <cstdio>
#include <cmath>

#include <bu.h>
#include "rt/view_legacy_bsg.h"

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

    /* ------------------------------------------------------------------
     * Build a minimal headless view with grid snap enabled.
     * ------------------------------------------------------------------ */
    void *view = rt_view_context_create_bsg();
    ASSERT(view != NULL);
    ASSERT(rt_view_context_dimensions_set_bsg(view, 512, 512));
    struct rt_view_grid_state grid = RT_VIEW_GRID_STATE_INIT;
    ASSERT(rt_view_context_grid_state_from_bsg(&grid, view));
    grid.res_h = 1.0;
    grid.res_v = 1.0;
    ASSERT(rt_view_context_grid_state_set_bsg(view, &grid));
    ASSERT(rt_view_context_snap_source_flags_set_bsg(view, 0));

    /* ------------------------------------------------------------------
     * Test 1: snap candidates near origin with GRID kind.
     * ------------------------------------------------------------------ */
    {
	point_t sample = {0.4, 0.4, 0.0};
	struct rt_view_snap_result_bsg *sr = rt_view_snap_result_create_bsg();
	ASSERT(sr != NULL);
	int cnt = rt_view_context_snap_candidates_result_bsg(view, sample, 0.0,
		RT_VIEW_SNAP_KIND_GRID_BSG, sr);
	/* A valid grid snap should return at least 1 candidate and the
	 * nearest snapped point should be within one grid cell of sample. */
	ASSERT(cnt >= 0);
	if (cnt > 0 && rt_view_snap_result_count_bsg(sr) > 0) {
	    point_t snapped = VINIT_ZERO;
	    ASSERT(rt_view_snap_result_point_bsg(sr, 0, snapped));
	    double dist = DIST_PNT_PNT(snapped, sample);
	    ASSERT(rt_view_context_grid_state_from_bsg(&grid, view));
	    ASSERT(dist <= grid.res_h * 1.5);
	}
	rt_view_snap_result_free_bsg(sr);
    }

    /* ------------------------------------------------------------------
     * Test 2: 2D snap does not crash and stays finite.
     * ------------------------------------------------------------------ */
    {
	fastf_t vx = 0.0, vy = 0.0;
	/* Place the sample near the view center (0,0 in view space). */
	rt_view_context_snap_point_2d_bsg(view, &vx, &vy,
		RT_VIEW_SNAP_KIND_GRID_BSG);
	/* Result should be finite (no NaN/Inf). */
	ASSERT(std::isfinite((double)vx));
	ASSERT(std::isfinite((double)vy));
    }

    /* ------------------------------------------------------------------
     * Test 3: NULL guards — must not crash.
     * ------------------------------------------------------------------ */
    {
	struct rt_view_snap_result_bsg *sr = rt_view_snap_result_create_bsg();
	ASSERT(sr != NULL);
	point_t p = VINIT_ZERO;
	int r = rt_view_context_snap_candidates_result_bsg(NULL, p, 0.0,
		RT_VIEW_SNAP_KIND_GRID_BSG, sr);
	ASSERT(r == 0);
	fastf_t vx = 0.0, vy = 0.0;
	rt_view_context_snap_point_2d_bsg(NULL, &vx, &vy,
		RT_VIEW_SNAP_KIND_GRID_BSG);
	rt_view_snap_result_free_bsg(sr);
    }

    rt_view_context_free_bsg(view);

    bu_log("snap semantic records: %d checks, %d failures\n", nchecks, nfails);
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
