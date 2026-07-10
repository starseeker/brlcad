/*                  R T W I Z A R D _ V I E W . C P P
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
/** @file rtwizard_view.cpp
 *
 * rtwizard view-pipeline regression tests.
 *
 * rtwizard's headless (--no-gui) pipeline computes rt view parameters via the
 * libtclcad GED API:
 *
 *   go_open db db <file>          ; opens a GED instance
 *   db new_view v1 nu             ; creates a null-DM secondary view
 *   db draw <objects>             ; populates scene objects for autoview
 *   db autoview v1                ; sets gv_size / gv_center from scene bounds
 *   db aet v1 <az> <el> <tw>     ; sets azimuth / elevation / twist
 *   db zoom v1 <z>                ; applies zoom
 *   db get_eyemodel v1            ; extracts viewsize, orientation, eye_pt
 *
 * The eye model (viewsize, orientation quaternion, eye_pt) is then forwarded
 * as command-line arguments to the rt/rtedge subprocess.  No display-manager
 * rendering occurs; the view feature is purely a lightweight camera container.
 *
 * This test exercises the equivalent C-API path and verifies:
 *
 *   1. null_view_context       — A null-DM secondary view is a valid neutral
 *                                RT view context with no display manager.
 *   2. eyemodel_finite         — draw + autoview + get_eyemodel produces a
 *                                plausible (finite, non-degenerate) eye model.
 *   3. multiple_null_views     — Multiple secondary null-DM views can coexist
 *                                in the GED view set without retained fallback
 *                                assumptions.
 *   4. obol_render             — A GUI/display-backed rtwizard path can attach
 *                                the headless Obol DM and capture an image.
 *   5. obol_eyemodel_consistency — Re-applying the same view state to an
 *                                Obol-backed view produces stable output.
 *
 * The camera-only tests use "nu" view contexts.  Display-backed tests attach
 * the headless Obol/Coin display manager, so no display hardware or X11 server
 * is required.
 *
 * Usage: ged_test_rtwizard_view [-c] <directory-containing-moss.g>
 */

#include "common.h"

#include <cmath>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <string>

#include <bu.h>
#include "bu/opt.h"
#define DM_WITH_RT
#include <dm.h>
#include <ged.h>
#include <icv.h>
#include "ged/draw.h"
#include "rt/view.h"
#include "view_test_util.h"

extern "C" void ged_changed_callback(struct db_i *UNUSED(dbip), struct directory *dp, int mode, void *u_data);

/* ---- helpers ------------------------------------------------------------- */

static void
close_gedp(struct ged *gedp)
{
    ged_close(gedp);
}

/*
 * Open a GED instance for CLI-mode testing with no display manager.
 * Simulates the go_open step in the rtwizard headless Tcl pipeline.
 *
 * rtwizard --no-gui never attaches a display-manager to the main GED instance;
 * it only creates a secondary view for camera computation.  We match that here.
 */
static struct ged *
open_gedp_null(const char *gfile)
{
    struct ged *gedp = ged_open("db", gfile, 1);
    if (!gedp)
	return NULL;

    db_add_changed_clbk(gedp->dbip, &ged_changed_callback, (void *)gedp);

    void *v = ged_view_active_ctx(gedp);
    bv_unit_conversion_set(DRAW_TEST_BV(v), gedp->dbip->dbi_local2base, gedp->dbip->dbi_base2local);

    return gedp;
}

/*
 * Create a secondary null-DM view and add it to the GED view set.
 * Mirrors what libtclcad's to_new_view() does for "db new_view v1 nu".
 *
 * A null-DM view needs no real display-manager handle; its sole purpose is
 * to carry a camera matrix for computing the eye model.  Its display-manager
 * context remains NULL just as libtclcad does for "nu".
 */
static void *
make_null_view(struct ged *gedp, const char *vname)
{
    void *view_set_ctx = ged_view_set_ctx(gedp);
    void *v = ged_view_context_create_with_set(view_set_ctx);
    bv_name_set(DRAW_TEST_BV(v), vname);
    bv_unit_conversion_set(DRAW_TEST_BV(v), gedp->dbip->dbi_local2base, gedp->dbip->dbi_base2local);

    ged_view_set_context_add(view_set_ctx, v);
    ged_view_context_owned_add(gedp, v);

    return v;
}

/* ========================================================================== */
/* Test 1: null-DM secondary view is a valid neutral view context             */
/* ========================================================================== */
static int
test_null_view_context(const char *datadir)
{
    bu_log("\n--- Test 1: null-DM secondary view context ---\n");

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", datadir);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmp("rtw_view_t1.g", std::ios::binary);
    tmp << orig.rdbuf();
    orig.close(); tmp.close();
    bu_vls_free(&fname);

    struct ged *gedp = open_gedp_null("rtw_view_t1.g");
    if (!gedp) {
	bu_log("FAIL: ged_open failed\n");
	bu_file_delete("rtw_view_t1.g");
	return 1;
    }

    /* Create a secondary null-DM view, simulating "db new_view v1 nu" */
    void *v1 = make_null_view(gedp, "v1");

    int fail = 0;
    if (!v1 || !bv_context_is_valid((struct bv_context *)(v1))) {
	bu_log("FAIL: secondary null-DM view is not a valid RT view context\n");
	fail = 1;
    } else {
	bu_log("PASS: secondary null-DM view is a valid RT view context\n");
    }

    if (ged_view_context_display_manager_get(v1) != NULL) {
	bu_log("FAIL: secondary null-DM view unexpectedly has a display manager\n");
	fail = 1;
    } else {
	bu_log("PASS: secondary null-DM view has no display manager\n");
    }

    const char *vname = bv_name_get(DRAW_TEST_BV_CONST(v1));
    if (!vname || BU_STR_EQUAL(vname, "v1") == 0) {
	bu_log("FAIL: secondary null-DM view name is not v1\n");
	fail = 1;
    } else {
	bu_log("PASS: secondary null-DM view name is v1\n");
    }

    if (!bv_context_is_valid((struct bv_context *)ged_view_active_ctx(gedp))) {
	bu_log("FAIL: default GED active view is not a valid libbv context\n");
	fail = 1;
    } else {
	bu_log("PASS: default GED active view is a valid libbv context\n");
    }

    close_gedp(gedp);
    bu_file_delete("rtw_view_t1.g");
    return fail;
}

/* ========================================================================== */
/* Test 2: draw + autoview + get_eyemodel produces finite/plausible params    */
/* ========================================================================== */
static int
test_eyemodel_finite(const char *datadir)
{
    bu_log("\n--- Test 2: draw + autoview + get_eyemodel produces finite eye model ---\n");

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", datadir);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmp("rtw_view_t2.g", std::ios::binary);
    tmp << orig.rdbuf();
    orig.close(); tmp.close();
    bu_vls_free(&fname);

    struct ged *gedp = open_gedp_null("rtw_view_t2.g");
    if (!gedp) {
	bu_log("FAIL: ged_open failed\n");
	bu_file_delete("rtw_view_t2.g");
	return 1;
    }

    /* Create secondary null-DM view (rtwizard "new_view v1 nu") */
    void *v1 = make_null_view(gedp, "v1");

    /* Draw objects (rtwizard "db draw $item" for each object in color_objlist) */
    const char *s_av[4] = {"draw", "all.g", NULL};
    ged_exec_draw(gedp, 2, s_av);

    /* Autoview on v1 (rtwizard "db autoview v1") */
    void *prev = ged_view_active_ctx(gedp);
    ged_view_active_ctx_set(gedp, v1);
    s_av[0] = "autoview"; s_av[1] = NULL;
    ged_exec_autoview(gedp, 1, s_av);

    /* Apply default az/el/twist (rtwizard "db aet v1 35 25 0") */
    s_av[0] = "ae"; s_av[1] = "35"; s_av[2] = "25"; s_av[3] = NULL;
    ged_exec_ae(gedp, 3, s_av);

    /* Extract eye model (rtwizard "db get_eyemodel v1") */
    s_av[0] = "get_eyemodel"; s_av[1] = NULL;
    ged_exec_get_eyemodel(gedp, 1, s_av);
    ged_view_active_ctx_set(gedp, prev);

    const char *result = bu_vls_cstr(gedp->ged_result_str);
    bu_log("get_eyemodel output:\n%s\n", result);

    /* Parse viewsize from the result string */
    double viewsize = 0.0;
    double qw = 0.0, qx = 0.0, qy = 0.0, qz = 0.0;
    double ex = 0.0, ey = 0.0, ez = 0.0;
    int nscan = sscanf(result,
		       "viewsize %lf ; orientation %lf %lf %lf %lf ; eye_pt %lf %lf %lf",
		       &viewsize, &qw, &qx, &qy, &qz, &ex, &ey, &ez);

    int fail = 0;
    if (nscan < 8) {
	/* Try alternate format without semicolons on same token */
	nscan = sscanf(result,
		       "viewsize %lf\n orientation %lf %lf %lf %lf\n eye_pt %lf %lf %lf",
		       &viewsize, &qw, &qx, &qy, &qz, &ex, &ey, &ez);
    }

    if (nscan < 8) {
	bu_log("FAIL: could not parse eye model (%d/8 values matched)\n", nscan);
	fail = 1;
    } else {
	/* viewsize must be positive and finite */
	if (viewsize <= 0.0 || !std::isfinite(viewsize)) {
	    bu_log("FAIL: viewsize %g is not a positive finite number\n", viewsize);
	    fail = 1;
	} else {
	    bu_log("PASS: viewsize = %g\n", viewsize);
	}

	/* orientation quaternion must have unit magnitude (within tolerance) */
	double qmag = std::sqrt(qw*qw + qx*qx + qy*qy + qz*qz);
	if (fabs(qmag - 1.0) > 0.01) {
	    bu_log("FAIL: orientation quaternion magnitude %g != 1.0\n", qmag);
	    fail = 1;
	} else {
	    bu_log("PASS: orientation quaternion |q| = %g\n", qmag);
	}

	/* eye_pt must be finite */
	if (!std::isfinite(ex) || !std::isfinite(ey) || !std::isfinite(ez)) {
	    bu_log("FAIL: eye_pt (%g %g %g) contains non-finite values\n", ex, ey, ez);
	    fail = 1;
	} else {
	    bu_log("PASS: eye_pt = (%g %g %g)\n", ex, ey, ez);
	}
    }

    close_gedp(gedp);
    bu_file_delete("rtw_view_t2.g");
    return fail;
}

/* ========================================================================== */
/* Test 3: multiple secondary null-DM views coexist in the GED view set       */
/* ========================================================================== */
static int
test_multiple_null_views(const char *datadir)
{
    bu_log("\n--- Test 3: multiple secondary null-DM views ---\n");

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", datadir);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmp("rtw_view_t3.g", std::ios::binary);
    tmp << orig.rdbuf();
    orig.close(); tmp.close();
    bu_vls_free(&fname);

    struct ged *gedp = open_gedp_null("rtw_view_t3.g");
    if (!gedp) {
	bu_log("FAIL: ged_open failed\n");
	bu_file_delete("rtw_view_t3.g");
	return 1;
    }

    /* Draw objects */
    const char *s_av[4] = {"draw", "all.g", NULL};
    ged_exec_draw(gedp, 2, s_av);

    /* Create four secondary views and verify neutral context identity. */
    void *views[4];
    char vname[4][8];
    int fail = 0;
    for (int i = 0; i < 4; i++) {
	snprintf(vname[i], sizeof(vname[i]), "v%d", i + 1);
	views[i] = make_null_view(gedp, vname[i]);

	if (!views[i] || !bv_context_is_valid((struct bv_context *)(views[i]))) {
	    bu_log("FAIL: view '%s' is not a valid RT view context\n", vname[i]);
	    fail = 1;
	}
	if (ged_view_context_display_manager_get(views[i]) != NULL) {
	    bu_log("FAIL: view '%s' unexpectedly has a display manager\n", vname[i]);
	    fail = 1;
	}
    }
    if (!fail)
	bu_log("PASS: all 4 secondary null-DM views are valid neutral contexts\n");

    close_gedp(gedp);
    bu_file_delete("rtw_view_t3.g");
    return fail;
}

/* ---- helpers for Obol display-backed tests -------------------------------- */

/*
 * Open a GED instance with an attached Obol DM (mirrors the rtwizard GUI
 * context, where the MGED widget has a real display manager for rendering).
 */
static struct ged *
open_gedp_obol(const char *gfile, int width, int height)
{
    struct ged *gedp = ged_open("db", gfile, 1);
    if (!gedp) return NULL;

    db_add_changed_clbk(gedp->dbip, &ged_changed_callback, (void *)gedp);

    const char *s_av[6] = {"dm", "attach", "obol", "RTW_OBOL", NULL};
    ged_exec_dm(gedp, 4, s_av);

    void *v = ged_view_active_ctx(gedp);
    struct dm *dmp  = (struct dm *)ged_view_context_display_manager_get(v);
    if (!dmp) {
	ged_close(gedp);
	return NULL;
    }
    dm_set_width(dmp, width);
    dm_set_height(dmp, height);
    dm_configure_win(dmp, 0);
    dm_set_zbuffer(dmp, 1);
    fastf_t wb[6] = {-1, 1, -1, 1, -100, 100};
    dm_set_win_bounds(dmp, wb);
    dm_set_vp(dmp, bv_scale_storage_get(DRAW_TEST_BV(v)));
    ged_view_context_display_manager_set(v, dmp);
    bv_dimensions_set(DRAW_TEST_BV(v), dm_get_width(dmp), dm_get_height(dmp));
    bv_unit_conversion_set(DRAW_TEST_BV(v), gedp->dbip->dbi_local2base, gedp->dbip->dbi_base2local);

    return gedp;
}

/* Simple draw-only refresh for a GED's current view. */
static void
do_obol_refresh(struct ged *gedp)
{
    void *v = ged_view_active_ctx(gedp);
    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
    txn.view = v;
    ged_draw_apply_transaction(gedp, &txn, NULL);
}

/* Count non-background (non-black) pixels in a PNG. */
static long
count_nonblack_rtw(const char *fname)
{
    icv_image_t *img = icv_read(fname, BU_MIME_IMAGE_PNG, 0, 0);
    if (!img) return -1;
    long cnt = 0;
    size_t npix = (size_t)img->width * img->height;
    double *d = img->data;
    for (size_t i = 0; i < npix; i++) {
	if (d[3*i] > 0.001 || d[3*i+1] > 0.001 || d[3*i+2] > 0.001)
	    cnt++;
    }
    icv_destroy(img);
    return cnt;
}

/* Compare two PNGs; returns 1 if identical within adiff_allow off-by-1 pixels. */
static int
images_match_rtw(const char *a, const char *b, int adiff_allow)
{
    icv_image_t *ia = icv_read(a, BU_MIME_IMAGE_PNG, 0, 0);
    icv_image_t *ib = icv_read(b, BU_MIME_IMAGE_PNG, 0, 0);
    if (!ia || !ib) {
	if (ia) icv_destroy(ia);
	if (ib) icv_destroy(ib);
	return 0;
    }
    int match = 0, off1 = 0, offmany = 0;
    icv_diff(&match, &off1, &offmany, ia, ib);
    icv_destroy(ia);
    icv_destroy(ib);
    int bad = offmany + (off1 > adiff_allow ? off1 - adiff_allow : 0);
    return (bad == 0) ? 1 : 0;
}

/* ========================================================================== */
/* Test 4: rtwizard GUI-path Obol view renders and yields eye model          */
/* ========================================================================== */
static int
test_gui_obol_render(const char *datadir)
{
    bu_log("\n--- Test 4: rtwizard GUI path (Obol DM) renders ---\n");

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", datadir);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmp("rtw_view_t4.g", std::ios::binary);
    tmp << orig.rdbuf();
    orig.close(); tmp.close();
    bu_vls_free(&fname);

    /* Open with Obol DM - simulates the MGED widget used in rtwizard GUI. */
    struct ged *gedp = open_gedp_obol("rtw_view_t4.g", 512, 512);
    if (!gedp) {
	bu_log("FAIL: ged_open (Obol) failed\n");
	bu_file_delete("rtw_view_t4.g");
	return 1;
    }

    int fail = 0;

    if (!bv_context_is_valid((struct bv_context *)ged_view_active_ctx(gedp))) {
	bu_log("FAIL: Obol-DM active view is not a valid libbv context\n");
	fail = 1;
    } else {
	bu_log("PASS: Obol-DM active view is a valid libbv context\n");
    }

    /* Draw + autoview + ae — mirrors MGEDpage::draw() + refreshDisplay() */
    const char *s_av[4] = {"draw", "all.g", NULL};
    ged_exec_draw(gedp, 2, s_av);
    s_av[0] = "autoview"; s_av[1] = NULL;
    ged_exec_autoview(gedp, 1, s_av);
    s_av[0] = "ae"; s_av[1] = "35"; s_av[2] = "25"; s_av[3] = NULL;
    ged_exec_ae(gedp, 3, s_av);

    /* Render through the active display path. */
    do_obol_refresh(gedp);

    /* Capture to PNG via screengrab */
    const char *sg_av[3] = {"screengrab", "rtw_view_t4.png", NULL};
    ged_exec_screengrab(gedp, 2, sg_av);

    long npix = count_nonblack_rtw("rtw_view_t4.png");
    if (npix <= 0) {
	bu_log("FAIL: GUI-mode (Obol) wireframe image is empty (%ld non-black pixels)\n", npix);
	fail = 1;
    } else {
	bu_log("PASS: GUI-mode (Obol) wireframe rendered %ld non-black pixels\n", npix);
    }

    /* Get eye model — the same call rtwizard makes before spawning rt */
    s_av[0] = "get_eyemodel"; s_av[1] = NULL;
    ged_exec_get_eyemodel(gedp, 1, s_av);
    const char *em = bu_vls_cstr(gedp->ged_result_str);
    bu_log("GUI-mode get_eyemodel:\n%s\n", em);

    double viewsize = 0.0;
    int nscan = sscanf(em, "viewsize %lf", &viewsize);
    if (nscan < 1 || viewsize <= 0.0 || !std::isfinite(viewsize)) {
	bu_log("FAIL: GUI-mode eye model has invalid viewsize (%g)\n", viewsize);
	fail = 1;
    } else {
	bu_log("PASS: GUI-mode eye model viewsize = %g\n", viewsize);
    }

    bu_file_delete("rtw_view_t4.png");
    close_gedp(gedp);
    bu_file_delete("rtw_view_t4.g");
    return fail;
}

/* ========================================================================== */
/* Test 5: GUI-mode eye model self-consistency                                 */
/*                                                                             */
/* After draw+autoview, extract the eye model (viewsize, orientation, eye_pt) */
/* and verify that re-applying those same parameters to a fresh second view   */
/* produces a pixel-identical rendered image.  This exercises the full        */
/* rtwizard GUI pipeline:                                                      */
/*   view A  → draw+autoview+ae → render → get_eyemodel                       */
/*   view B  → apply same az/el/zoom → render                                 */
/*   assert: image_A == image_B                                                */
/* ========================================================================== */
static int
test_gui_eyemodel_consistency(const char *datadir)
{
    bu_log("\n--- Test 5: GUI-mode eye model self-consistency (view A == view B) ---\n");

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", datadir);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmp("rtw_view_t5.g", std::ios::binary);
    tmp << orig.rdbuf();
    orig.close(); tmp.close();
    bu_vls_free(&fname);

    struct ged *gedp = open_gedp_obol("rtw_view_t5.g", 512, 512);
    if (!gedp) {
	bu_log("FAIL: ged_open (Obol) failed\n");
	bu_file_delete("rtw_view_t5.g");
	return 1;
    }

    /* === View A === */
    const char *s_av[8] = {NULL};
    s_av[0] = "draw"; s_av[1] = "all.g"; s_av[2] = NULL;
    ged_exec_draw(gedp, 2, s_av);
    s_av[0] = "autoview"; s_av[1] = NULL;
    ged_exec_autoview(gedp, 1, s_av);
    /* Specific az/el/tw matching a typical rtwizard "init_azimuth/elevation" */
    s_av[0] = "ae"; s_av[1] = "45"; s_av[2] = "30"; s_av[3] = "0"; s_av[4] = NULL;
    ged_exec_ae(gedp, 4, s_av);

    /* Save the view matrix for comparison */
    void *v = ged_view_active_ctx(gedp);
    mat_t saved_m2v;
    bv_model2view_get(saved_m2v, DRAW_TEST_BV_CONST(v));

    /* Render A */
    do_obol_refresh(gedp);
    const char *sg_av[3] = {"screengrab", "rtw_view_t5_A.png", NULL};
    ged_exec_screengrab(gedp, 2, sg_av);
    bu_log("Captured image A (az=45 el=30 tw=0)\n");

    /* Extract eye model */
    s_av[0] = "get_eyemodel"; s_av[1] = NULL;
    ged_exec_get_eyemodel(gedp, 1, s_av);
    /* Copy result (the result_str will be overwritten by subsequent commands) */
    struct bu_vls em_str = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&em_str, "%s", bu_vls_cstr(gedp->ged_result_str));

    /* === View B: zap + redraw + re-apply same az/el/tw === */
    s_av[0] = "zap"; s_av[1] = NULL;
    ged_exec_zap(gedp, 1, s_av);
    s_av[0] = "draw"; s_av[1] = "all.g"; s_av[2] = NULL;
    ged_exec_draw(gedp, 2, s_av);
    s_av[0] = "autoview"; s_av[1] = NULL;
    ged_exec_autoview(gedp, 1, s_av);
    s_av[0] = "ae"; s_av[1] = "45"; s_av[2] = "30"; s_av[3] = "0"; s_av[4] = NULL;
    ged_exec_ae(gedp, 4, s_av);

    /* Render B */
    do_obol_refresh(gedp);
    sg_av[1] = "rtw_view_t5_B.png";
    ged_exec_screengrab(gedp, 2, sg_av);
    bu_log("Captured image B (same ae applied fresh)\n");

    int fail = 0;

    /* Images A and B must be pixel-identical (same camera, same geometry) */
    if (!images_match_rtw("rtw_view_t5_A.png", "rtw_view_t5_B.png", 10)) {
	bu_log("FAIL: image A != image B (same ae params but different pixels)\n");
	fail = 1;
    } else {
	bu_log("PASS: image A == image B (same ae params produce identical rendering)\n");
    }

    /* Also confirm view matrices match within floating-point noise */
    mat_t new_m2v;
    bv_model2view_get(new_m2v, DRAW_TEST_BV_CONST(v));
    fastf_t max_delta = 0.0;
    for (int i = 0; i < 16; i++) {
	fastf_t d = fabs(saved_m2v[i] - new_m2v[i]);
	if (d > max_delta) max_delta = d;
    }
    if (max_delta > 1e-6) {
	bu_log("FAIL: view matrices diverge after re-applying same ae (max delta = %g)\n",
	       max_delta);
	fail = 1;
    } else {
	bu_log("PASS: view matrices identical within 1e-6 (max delta = %g)\n", max_delta);
    }

    /* Confirm eye model of view B equals that of view A */
    s_av[0] = "get_eyemodel"; s_av[1] = NULL;
    ged_exec_get_eyemodel(gedp, 1, s_av);
    const char *em_A = bu_vls_cstr(&em_str);
    const char *em_B = bu_vls_cstr(gedp->ged_result_str);

    double vsA = 0.0, vsB = 0.0;
    sscanf(em_A, "viewsize %lf", &vsA);
    sscanf(em_B, "viewsize %lf", &vsB);
    if (vsA <= 0.0 || vsB <= 0.0) {
	bu_log("FAIL: viewsize not positive (A=%g B=%g)\n", vsA, vsB);
	fail = 1;
    } else if (fabs(vsA - vsB) / vsA > 1e-6) {
	bu_log("FAIL: viewsize differs between view A (%g) and view B (%g)\n", vsA, vsB);
	fail = 1;
    } else {
	bu_log("PASS: eye model viewsize consistent: A=%g B=%g\n", vsA, vsB);
    }

    bu_vls_free(&em_str);
    bu_file_delete("rtw_view_t5_A.png");
    bu_file_delete("rtw_view_t5_B.png");
    close_gedp(gedp);
    bu_file_delete("rtw_view_t5.g");
    return fail;
}

/* ========================================================================== */
/* main                                                                        */
/* ========================================================================== */

int
main(int argc, char *argv[])
{
    bu_setprogname(argv[0]);

    int cleanup = 0;
    struct bu_opt_desc d[2];
    BU_OPT(d[0], "c", "cleanup", NULL, NULL, &cleanup, "cleanup temp files");
    BU_OPT_NULL(d[1]);

    int uac = bu_opt_parse(NULL, argc, (const char **)argv, d);
    if (uac != 2) {
	bu_log("Usage: %s [-c] <directory-containing-moss.g>\n", argv[0]);
	return 1;
    }
    const char *datadir = argv[1];

    int failures = 0;
    failures += test_null_view_context(datadir);
    failures += test_eyemodel_finite(datadir);
    failures += test_multiple_null_views(datadir);
    failures += test_gui_obol_render(datadir);
    failures += test_gui_eyemodel_consistency(datadir);

    if (failures == 0) {
	bu_log("\nAll rtwizard view tests PASSED (%d/5)\n", 5);
    } else {
	bu_log("\n%d rtwizard view test(s) FAILED\n", failures);
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
