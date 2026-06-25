/*           M E A S U R E _ P A R I T Y . C P P
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
/** @file measure_semantics.cpp
 *
 * Verify that the transitional RT measure adapter returns semantic distance
 * records consistent with DIST_PNT_PNT for known point pairs.
 *
 * Strategy:
 *   1. Create a headless view.
 *   2. Call the RT measure adapter for two known points.
 *   3. If the result is valid, verify distance matches DIST_PNT_PNT(a, b) within
 *      a small epsilon (measure candidates are allowed to refine the distance,
 *      but must not be wildly inconsistent with the straight-line distance).
 *   4. Verify NULL guards do not crash.
 *
 * No database or display hardware required.
 *
 * Usage: ged_test_measure_semantics
 */

#include "common.h"

#include <cstdio>
#include <cmath>

#include <bu.h>
#include "rt/view.h"

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
     * Build a minimal headless view.
     * ------------------------------------------------------------------ */
    void *view = rt_view_context_create();
    ASSERT(view != NULL);
    ASSERT(rt_view_context_dimensions_set(view, 512, 512));

    /* ------------------------------------------------------------------
     * Test 1: known distance.
     * ------------------------------------------------------------------ */
    {
	point_t a = {0.0, 0.0, 0.0};
	point_t b = {3.0, 4.0, 0.0};   /* classic 3-4-5 triangle: dist = 5.0 */
	double expected = DIST_PNT_PNT(a, b);
	ASSERT(fabs(expected - 5.0) < 1e-9);

	struct rt_view_measure_result mr = RT_VIEW_MEASURE_RESULT_INIT;
	int rc = rt_view_context_measure_candidates(view, a, b, &mr);
	ASSERT(rc >= 0);
	if (mr.valid) {
	    /* distance must not deviate from straight-line by more than
	     * 1 % (candidates may snap, but not teleport). */
	    double err = fabs(mr.distance - expected);
	    ASSERT(err / expected < 0.01);
	}
    }

    /* ------------------------------------------------------------------
     * Test 2: zero-length segment - must not crash, distance >= 0.
     * ------------------------------------------------------------------ */
    {
	point_t a = {1.0, 2.0, 3.0};
	point_t b = {1.0, 2.0, 3.0};
	struct rt_view_measure_result mr = RT_VIEW_MEASURE_RESULT_INIT;
	int rc = rt_view_context_measure_candidates(view, a, b, &mr);
	ASSERT(rc >= 0);
	if (mr.valid)
	    ASSERT(mr.distance >= 0.0);
    }

    /* ------------------------------------------------------------------
     * Test 3: NULL guards - must not crash.
     * ------------------------------------------------------------------ */
    {
	point_t a = VINIT_ZERO;
	point_t b = VINIT_ZERO;
	/* Passing NULL for out must return 0 without crashing. */
	int rc = rt_view_context_measure_candidates(view, a, b, NULL);
	ASSERT(rc == 0);
	/* Passing NULL view is also valid (falls back to no projection). */
	struct rt_view_measure_result mr = RT_VIEW_MEASURE_RESULT_INIT;
	rc = rt_view_context_measure_candidates(NULL, a, b, &mr);
	ASSERT(rc == 0);   /* zero-length segment returns 0 regardless of view */
    }

    rt_view_context_free(view);

    bu_log("measure semantic records: %d checks, %d failures\n", nchecks, nfails);
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
