/*              M G E D _ S H A D E D _ M O D E . C P P
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
/** @file mged_shaded_mode.cpp
 *
 * MGED shaded-mode, display-state, and draw-record regression test.
 *
 * Verifies the lighting / zbuffer / zclip plumbing used by
 * mged_shaded_mode_helper and the GED draw-mode pipeline:
 *
 *   Test 1 — endpoint renderer-policy and view-state round-trip
 *             Exercises what mged_shaded_mode_helper does:
 *               dm set zbuffer $val; dm set zclip $val; dm set lighting $val
 *             Verifies the GED command updates Obol renderer and libbv state.
 *
 *   Test 2 — shaded draw-mode pipeline (draw -m1, draw -m2) produces expected
 *             Obol-rendered output relative to wireframe (draw -m0).
 *             This test uses the headless Obol/Coin off-screen renderer so it
 *             does not require display hardware or the legacy retained-DM
 *             rendering path.
 *
 *   Test 3 — Obol-rendered output stability across zap/redraw
 *             Draw once, capture the non-black pixel count.  Zap and re-draw.
 *             The second count must match the first within a small tolerance.
 *
 * Uses the headless Obol/Coin off-screen renderer for image output; no display
 * hardware required.
 *
 * Usage: ged_test_mged_shaded_mode <directory-containing-moss.g>
 */

#include "common.h"

#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include <bu.h>
#include "bu/opt.h"
#include "rt/view.h"
#include "view_test_util.h"
#include <ged.h>
#include "ged/draw.h"
#include "brlobol/display_endpoint.h"
#include "brlobol/view_controller.h"

extern "C" void ged_changed_callback(struct db_i *UNUSED(dbip), struct directory *dp, int mode, void *u_data);
extern "C" long draw_test_count_nonblack_pixels(const char *filename);
extern "C" int draw_test_images_differ(const char *a, const char *b, int offmany_threshold);

/* ---- shared helpers ------------------------------------------------------ */

static void
do_refresh(struct ged *gedp)
{
    void *v = ged_view_active_ctx(gedp);
    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
    txn.view = v;
    ged_draw_apply_transaction(gedp, &txn, NULL);
}

/* Screen grabs are a point-in-time assertion, so settle the small test model
 * rather than accidentally comparing its shared deferred root each time. */
static long
capture_obol_pixels(struct ged *gedp, const char *filename,
	const char *label, int require_nonempty)
{
    do_refresh(gedp);
    const char *sg_av[3] = {"screengrab", filename, NULL};
    if (ged_exec_screengrab(gedp, 2, sg_av) & BRLCAD_ERROR) {
	bu_log("FAIL: screengrab failed for %s: %s\n", label,
		bu_vls_cstr(gedp->ged_result_str));
	return -1;
    }

    long pixels = draw_test_count_nonblack_pixels(filename);
    if (pixels <= 0 && require_nonempty) {
	bu_log("FAIL: Obol capture for %s is empty (%ld non-black pixels)\n",
		label, pixels);
	return -1;
    }

    bu_log("Obol capture for %s: %ld non-black pixels\n", label, pixels);
    return pixels;
}

static struct ged *
open_gedp(const char *gfile, int width, int height)
{
    struct ged *gedp = ged_open("db", gfile, 1);
    if (!gedp) return NULL;

    db_add_changed_clbk(gedp->dbip, &ged_changed_callback, (void *)gedp);

    void *v = ged_view_active_ctx(gedp);
    bv_dimensions_set(DRAW_TEST_BV(v), width, height);
    const char *s_av[16] = {
	"dm", "open", "--host", "headless", "--renderer", "sw", NULL
    };
    if (ged_exec_dm(gedp, 6, s_av) != BRLCAD_OK) {
	ged_close(gedp);
	return NULL;
    }
    bv_unit_conversion_set(DRAW_TEST_BV(v), gedp->dbip->dbi_local2base, gedp->dbip->dbi_base2local);

    s_av[0] = "ae"; s_av[1] = "35"; s_av[2] = "25"; s_av[3] = NULL;
    ged_exec_ae(gedp, 3, s_av);

    return gedp;
}

/* ========================================================================== */
/* Test 1: dm_set/get_light, dm_set/get_zbuffer, dm_set/get_zclip round-trip  */
/* Simulates the state flips that mged_shaded_mode_helper issues:              */
/*   dm set zbuffer $val; dm set zclip $val; dm set lighting $val             */
/* ========================================================================== */
static int
test_dm_lighting_flags(const char *datadir)
{
    bu_log("\n--- Test 1: dm lighting/zbuffer/zclip flag round-trip ---\n");
    bu_log("(simulates mged_shaded_mode_helper code path)\n");

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", datadir);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmp("smb_t1.g", std::ios::binary);
    tmp << orig.rdbuf();
    orig.close(); tmp.close();
    bu_vls_free(&fname);

    struct ged *gedp = open_gedp("smb_t1.g", 256, 256);
    if (!gedp) { bu_log("FAIL: ged_open failed\n"); bu_file_delete("smb_t1.g"); return 1; }

    void *v = ged_view_active_ctx(gedp);
    int fail = 0;

    /* --- endpoint-native ged_exec_dm path used by applications ------------ */
    brlobol_display_endpoint_t *endpoint =
	ged_view_context_display_endpoint_get(v);
    BRLObolViewController *controller = endpoint ?
	static_cast<BRLObolViewController *>(
	    brlobol_display_endpoint_controller(endpoint)) : NULL;
    if (!controller) {
	bu_log("FAIL: Obol DM did not publish a view endpoint\n");
	fail++;
    } else {
	const char *s_av[5] = {"dm", "set", "lighting", "1", NULL};
	if (ged_exec_dm(gedp, 4, s_av) != BRLCAD_OK ||
	    !controller->isLightingEnabled()) {
	    bu_log("FAIL: dm set lighting 1 did not enable Obol lighting\n");
	    fail++;
	}
	s_av[2] = "zbuffer";
	if (ged_exec_dm(gedp, 4, s_av) != BRLCAD_OK ||
	    !controller->isDepthTestEnabled()) {
	    bu_log("FAIL: dm set zbuffer 1 did not enable Obol depth testing\n");
	    fail++;
	}
	s_av[2] = "zclip";
	if (ged_exec_dm(gedp, 4, s_av) != BRLCAD_OK ||
	    !bv_zclip_get(DRAW_TEST_BV(v))) {
	    bu_log("FAIL: dm set zclip 1 did not enable libbv clipping\n");
	    fail++;
	}

	s_av[3] = "0";
	s_av[2] = "lighting";
	if (ged_exec_dm(gedp, 4, s_av) != BRLCAD_OK ||
	    controller->isLightingEnabled()) {
	    bu_log("FAIL: dm set lighting 0 did not disable Obol lighting\n");
	    fail++;
	}
	s_av[2] = "zbuffer";
	if (ged_exec_dm(gedp, 4, s_av) != BRLCAD_OK ||
	    controller->isDepthTestEnabled()) {
	    bu_log("FAIL: dm set zbuffer 0 did not disable Obol depth testing\n");
	    fail++;
	}
	s_av[2] = "zclip";
	if (ged_exec_dm(gedp, 4, s_av) != BRLCAD_OK ||
	    bv_zclip_get(DRAW_TEST_BV(v))) {
	    bu_log("FAIL: dm set zclip 0 did not disable libbv clipping\n");
	    fail++;
	}
	if (!fail)
	    bu_log("PASS: GED dm policy aliases update typed owners\n");
    }

    ged_close(gedp);
    bu_file_delete("smb_t1.g");
    return fail;
}

/* ========================================================================== */
/* Test 2: shaded draw modes produce different output from wireframe           */
/* ========================================================================== */
static int
test_shaded_mode_draw(const char *datadir)
{
    bu_log("\n--- Test 2: shaded draw modes produce non-empty and distinct output ---\n");

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", datadir);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmp("smb_t2.g", std::ios::binary);
    tmp << orig.rdbuf();
    orig.close(); tmp.close();
    bu_vls_free(&fname);

    struct ged *gedp = open_gedp("smb_t2.g", 512, 512);
    if (!gedp) { bu_file_delete("smb_t2.g"); return 1; }

    const char *s_av[8] = {NULL};
    int fail = 0;

    /* --- wireframe (mode 0) ----------------------------------------------- */
    s_av[0] = "draw"; s_av[1] = "-m0"; s_av[2] = "all.g"; s_av[3] = NULL;
    ged_exec_draw(gedp, 3, s_av);
    s_av[0] = "autoview"; s_av[1] = NULL;
    ged_exec_autoview(gedp, 1, s_av);
    if (!draw_test_obol_progressive_drain(gedp, ged_view_active_ctx(gedp),
	    500, 1)) {
	bu_log("FAIL: wireframe progressive realization did not settle\n");
	fail++;
    }

    long pix_m0 = capture_obol_pixels(gedp, "smb_t2_m0.png",
	    "wireframe mode 0", 1);
    bu_log("Wireframe (mode 0): %ld non-black pixels\n", pix_m0);
    if (pix_m0 <= 0) {
	bu_log("FAIL: wireframe rendered nothing\n");
	fail++;
    }

    /* --- shaded bots only (mode 1) ---------------------------------------- */
    s_av[0] = "Z"; s_av[1] = NULL;
    ged_exec_Z(gedp, 1, s_av);

    s_av[0] = "draw"; s_av[1] = "-m1"; s_av[2] = "all.g"; s_av[3] = NULL;
    ged_exec_draw(gedp, 3, s_av);
    s_av[0] = "autoview"; s_av[1] = NULL;
    ged_exec_autoview(gedp, 1, s_av);
    if (!draw_test_obol_progressive_drain(gedp, ged_view_active_ctx(gedp),
	    500, 1)) {
	bu_log("FAIL: shaded-bots progressive realization did not settle\n");
	fail++;
    }

    long pix_m1 = capture_obol_pixels(gedp, "smb_t2_m1.png",
	    "shaded bots mode 1", 1);
    bu_log("Shaded bots only (mode 1): %ld non-black pixels\n", pix_m1);
    if (pix_m1 <= 0) {
	bu_log("FAIL: shaded-bots mode rendered nothing\n");
	fail++;
    }

    /* --- shaded all (mode 2) ---------------------------------------------- */
    s_av[0] = "Z"; s_av[1] = NULL;
    ged_exec_Z(gedp, 1, s_av);

    s_av[0] = "draw"; s_av[1] = "-m2"; s_av[2] = "all.g"; s_av[3] = NULL;
    ged_exec_draw(gedp, 3, s_av);
    s_av[0] = "autoview"; s_av[1] = NULL;
    ged_exec_autoview(gedp, 1, s_av);
    if (!draw_test_obol_progressive_drain(gedp, ged_view_active_ctx(gedp),
	    500, 1)) {
	bu_log("FAIL: shaded-all progressive realization did not settle\n");
	fail++;
    }

    long pix_m2 = capture_obol_pixels(gedp, "smb_t2_m2.png",
	    "shaded all mode 2", 1);
    bu_log("Shaded all (mode 2): %ld non-black pixels\n", pix_m2);
    if (pix_m2 <= 0) {
	bu_log("FAIL: shaded-all mode rendered nothing\n");
	fail++;
    }

    /* shaded modes should produce visually different output from wireframe */
    if (pix_m1 > 0 && pix_m0 > 0) {
	if (draw_test_images_differ("smb_t2_m0.png", "smb_t2_m1.png", 10) > 0)
	    bu_log("PASS: mode 1 (shaded bots) differs from mode 0 (wireframe)\n");
	else
	    bu_log("INFO: mode 1 and mode 0 produced the same pixels "
		   "(may be OK if moss.g has no BoT primitives)\n");
    }
    if (pix_m2 > 0 && pix_m0 > 0) {
	if (draw_test_images_differ("smb_t2_m0.png", "smb_t2_m2.png", 10) > 0)
	    bu_log("PASS: mode 2 (shaded all) differs from mode 0 (wireframe)\n");
	else {
	    bu_log("FAIL: mode 2 (shaded all) and wireframe are identical — "
		   "shaded path is broken\n");
	    fail++;
	}
    }

    bu_file_delete("smb_t2_m0.png");
    bu_file_delete("smb_t2_m1.png");
    bu_file_delete("smb_t2_m2.png");
    bu_file_delete("smb_t2.g");
    ged_close(gedp);
    return fail;
}

/* ========================================================================== */
/* Test 3: Obol-rendered output is stable across zap/redraw                    */
/* ========================================================================== */
static int
test_obol_render_stability(const char *datadir)
{
    bu_log("\n--- Test 3: Obol render stability across zap/redraw ---\n");

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", datadir);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmp("smb_t4.g", std::ios::binary);
    tmp << orig.rdbuf();
    orig.close(); tmp.close();
    bu_vls_free(&fname);

    struct ged *gedp = open_gedp("smb_t4.g", 512, 512);
    if (!gedp) { bu_file_delete("smb_t4.g"); return 1; }

    int fail = 0;
    const char *s_av[4] = {NULL};

    /* Draw once */
    s_av[0] = "draw"; s_av[1] = "all.g"; s_av[2] = NULL;
    ged_exec_draw(gedp, 2, s_av);
    s_av[0] = "autoview"; s_av[1] = NULL;
    ged_exec_autoview(gedp, 1, s_av);

    long pix_first = capture_obol_pixels(gedp, "smb_t4_first.png",
	    "first draw", 1);
    bu_log("First render: %ld non-black pixels\n", pix_first);

    /* Zap and re-draw identically */
    s_av[0] = "Z"; s_av[1] = NULL;
    ged_exec_Z(gedp, 1, s_av);

    s_av[0] = "draw"; s_av[1] = "all.g"; s_av[2] = NULL;
    ged_exec_draw(gedp, 2, s_av);
    s_av[0] = "autoview"; s_av[1] = NULL;
    ged_exec_autoview(gedp, 1, s_av);

    long pix_second = capture_obol_pixels(gedp, "smb_t4_second.png",
	    "second draw after zap/redraw", 1);
    bu_log("Second render (after zap+redraw): %ld non-black pixels\n", pix_second);

    /* The rendered pixel counts should be the same */
    long diff = pix_first > pix_second ? pix_first - pix_second
	                                : pix_second - pix_first;
    if (diff > 500) {
	bu_log("FAIL: pixel count after zap/redraw differs by %ld from first render "
	       "(first=%ld second=%ld)\n", diff, pix_first, pix_second);
	fail++;
    } else {
	bu_log("PASS: pixel count consistent across zap/redraw (%ld vs %ld)\n",
	       pix_first, pix_second);
    }

    bu_file_delete("smb_t4_first.png");
    bu_file_delete("smb_t4_second.png");
    bu_file_delete("smb_t4.g");
    ged_close(gedp);
    return fail;
}

/* ========================================================================== */
/* main                                                                        */
/* ========================================================================== */

int
main(int argc, char *argv[])
{
    bu_setprogname(argv[0]);

    if (argc != 2) {
	bu_log("Usage: %s <directory-containing-moss.g>\n", argv[0]);
	return 1;
    }
    const char *datadir = argv[1];

    int failures = 0;
    failures += test_dm_lighting_flags(datadir);
    failures += test_shaded_mode_draw(datadir);
    failures += test_obol_render_stability(datadir);

    if (failures == 0) {
	bu_log("\nAll MGED shaded-mode tests PASSED (3/3)\n");
    } else {
	bu_log("\n%d MGED shaded-mode test(s) FAILED\n", failures);
    }
    return failures;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
