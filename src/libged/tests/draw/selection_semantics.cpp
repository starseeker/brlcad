/*          S E L E C T I O N _ P A R I T Y . C P P
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
/** @file selection_semantics.cpp
 *
 * Verify that view selection is consistent with opaque GED pick results
 * after a GED draw+pick cycle.
 *
 * Strategy:
 *   1. Open moss.g, draw the scene.
 *   2. Run a screen-center point pick; collect an opaque GED pick result.
 *   3. Apply the pick result through the GED selection API.
 *   4. Verify selected path callback data and selection count.
 *   5. Clear selection through the GED view-context API; verify count
 *      drops to 0.
 *
 * No display hardware required (headless GED).
 *
 * Usage: ged_test_selection_semantics <directory-containing-moss.g>
 */

#include "common.h"

#include <fstream>
#include <cstring>
#include <string>

#include <bu.h>
#include <ged.h>
#include <ged/draw.h>
#include "rt/view.h"
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

struct selected_path_state {
    int count;
    std::string last_path;
};

static void
selection_path_callback(const char *path, void *data)
{
    struct selected_path_state *state = (struct selected_path_state *)data;
    if (!state || !path || !path[0])
	return;
    state->count++;
    state->last_path = path;
}


int
main(int ac, char *av[])
{
    bu_setprogname(av[0]);

    if (ac != 2) {
	bu_log("Usage: %s <directory-containing-moss.g>\n", av[0]);
	return 1;
    }
    if (!bu_file_directory(av[1])) {
	bu_log("ERROR: [%s] is not a directory.\n", av[1]);
	return 2;
    }

    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    char lcache[MAXPATHLEN] = {0};
    bu_dir(lcache, MAXPATHLEN, BU_DIR_CURR,
	   "ged_selection_semantics_cache", NULL);
    bu_mkdir(lcache);
    bu_setenv("BU_DIR_CACHE", lcache, 1);

    /* Copy moss.g into a working file. */
    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", av[1]);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    if (!orig.is_open()) {
	bu_log("ERROR: cannot open %s\n", bu_vls_cstr(&fname));
	bu_vls_free(&fname);
	return 3;
    }

    struct bu_vls wname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&wname, "%s/moss_selection_semantics_tmp.g",
		   bu_dir(NULL, 0, BU_DIR_CURR, NULL));
    {
	std::ofstream dst(bu_vls_cstr(&wname), std::ios::binary);
	dst << orig.rdbuf();
    }
    bu_vls_free(&fname);

    struct ged *gedp = ged_open("db", bu_vls_cstr(&wname), 1);
    bu_vls_free(&wname);
    if (!gedp) {
	bu_log("ERROR: ged_open failed\n");
	return 4;
    }

    /* Draw the scene. */
    const char *draw_av[] = {"draw", "all.g", NULL};
    ged_exec(gedp, 2, draw_av);

    /* Obtain the default view. */
    struct ged_view_context *v = ged_view_active_ctx(gedp);
    ASSERT(v != NULL);

    if (!v) {
	ged_close(gedp);
	return nfails ? 1 : 0;
    }

    /* ------------------------------------------------------------------
     * Test 1: selection lifecycle round-trip
     * ------------------------------------------------------------------ */
    ASSERT(ged_view_selection_count(v) == 0);

    /* ------------------------------------------------------------------
     * Test 2: pick then add to selection
     * ------------------------------------------------------------------ */
    int cx = bv_width_get(DRAW_TEST_BV_CONST(v)) / 2;
    int cy = bv_height_get(DRAW_TEST_BV_CONST(v)) / 2;
    struct ged_pick_result *pr =
	ged_pick_point(v, cx, cy, 1);
    if (pr && ged_pick_result_count(pr) > 0) {
	struct selected_path_state path_state = {0, std::string()};
	ASSERT(ged_view_selection_set_pick(v, pr,
		    selection_path_callback, &path_state));
	ASSERT(path_state.count > 0);
	ASSERT(!path_state.last_path.empty());
	ASSERT(ged_view_selection_count(v) == 1);
    }

    /* Clear and verify */
    ASSERT(ged_view_selection_clear(v));
    ASSERT(ged_view_selection_count(v) == 0);

    ged_pick_result_free(pr);

    ged_close(gedp);

    bu_log("selection semantic records: %d checks, %d failures\n", nchecks, nfails);
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
