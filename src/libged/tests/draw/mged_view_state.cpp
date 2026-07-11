/*              M G E D _ V I E W _ S T A T E . C P P
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
/** @file mged_view_state.cpp
 *
 * MGED view-state regression tests not covered by the image baselines.
 *
 * Verifies behavior that historically lived in retained draw/view state but
 * is now reached through neutral GED/RT view facades:
 *
 *   1. Highlight state updates drawn shape records.
 *   2. Edit matrix state is accepted by the neutral view context.
 *   3. Obol DM attachments use per-view controller attachments.
 *   4. Obol DMs do not inherit stale standalone rt framebuffer devices.
 *   5. GL display hosts can own a per-view Obol render endpoint.
 *
 * Uses the headless Obol/Coin off-screen renderer; no display hardware
 * required.
 *
 * Usage: ged_test_mged_view_state [-c] <directory-containing-moss.g>
 */

#include "common.h"

#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include <bu.h>
#include "bu/opt.h"
#include "rt/edit.h"
#include "rt/view.h"
#include "view_test_util.h"
#define DM_WITH_RT
#include <dm.h>
#include <ged.h>
#include "ged/draw.h"
#include "ged/draw_obol.h"
#include "brlobol/view_attachment.h"
#include "brlobol/scene_controller.h"
#include "brlobol/view_controller.h"
#include "../../ged_private.h"
#include "../../ged_draw_view_private.h"

extern "C" void ged_changed_callback(struct db_i *UNUSED(dbip), struct directory *dp, int mode, void *u_data);
extern "C" long draw_test_count_nonblack_pixels(const char *filename);
extern "C" int draw_test_images_differ(const char *a, const char *b, int offmany_threshold);

/* Full redraw through the neutral GED draw transaction path. */
static void
do_full_refresh(struct ged *gedp)
{
    void *v = ged_view_active_ctx(gedp);
    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
    txn.view = v;
    ged_draw_apply_transaction(gedp, &txn, NULL);
}

static int
capture_screengrab_nonempty(struct ged *gedp, const char *filename,
	const char *label)
{
    const char *sg_av[3] = {"screengrab", filename, NULL};
    if (ged_exec_screengrab(gedp, 2, sg_av) & BRLCAD_ERROR) {
	bu_log("FAIL: screengrab failed for %s: %s\n", label,
		bu_vls_cstr(gedp->ged_result_str));
	return 1;
    }

    long pixels = draw_test_count_nonblack_pixels(filename);
    if (pixels <= 0) {
	bu_log("FAIL: Obol capture for %s is empty (%ld non-black pixels)\n",
		label, pixels);
	return 1;
    }

    bu_log("PASS: Obol capture for %s produced %ld non-black pixels\n",
	    label, pixels);
    return 0;
}

static BRLObolViewController *
obol_controller_for_view(void *view_ctx)
{
    return (BRLObolViewController *)ged_draw_obol_controller_opaque_for_view(
	    view_ctx);
}

static int
configure_obol_view(struct ged *gedp, void *view_ctx, int width, int height)
{
    struct dm *dmp = (struct dm *)ged_view_context_display_manager_get(view_ctx);
    if (!dmp)
	return 0;

    dm_set_width(dmp, width);
    dm_set_height(dmp, height);
    dm_configure_win(dmp, 0);
    dm_set_zbuffer(dmp, 1);
    fastf_t wb[6] = {-1, 1, -1, 1, -100, 100};
    dm_set_win_bounds(dmp, wb);
    dm_set_vp(dmp, bv_scale_storage_get(DRAW_TEST_BV(view_ctx)));
    ged_view_context_display_manager_set(view_ctx, dmp);
    bv_dimensions_set(DRAW_TEST_BV(view_ctx), dm_get_width(dmp), dm_get_height(dmp));
    bv_unit_conversion_set(DRAW_TEST_BV(view_ctx), gedp->dbip->dbi_local2base, gedp->dbip->dbi_base2local);
    return 1;
}

/* ---- common GED + Obol setup --------------------------------------------- */

static struct ged *
open_gedp(const char *gfile, int width, int height)
{
    struct ged *gedp = ged_open("db", gfile, 1);
    if (!gedp) return NULL;

    db_add_changed_clbk(gedp->dbip, &ged_changed_callback, (void *)gedp);

    /* Attach the headless Obol display manager. */
    const char *s_av[16] = {NULL};
    s_av[0] = "dm"; s_av[1] = "attach"; s_av[2] = "obol"; s_av[3] = "OBOL"; s_av[4] = NULL;
    ged_exec_dm(gedp, 4, s_av);

    void *v = ged_view_active_ctx(gedp);
    struct dm *dmp  = (struct dm *)ged_view_context_display_manager_get(v);
    if (!dmp) {
	ged_close(gedp);
	return NULL;
    }
    if (!ged_draw_obol_controller_opaque_for_view(v)) {
	bu_log("FAIL: Obol DM attach did not associate a GED view controller\n");
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

    s_av[0] = "ae"; s_av[1] = "35"; s_av[2] = "25"; s_av[3] = NULL;
    ged_exec_ae(gedp, 3, s_av);

    return gedp;
}

/* Test 1: Highlight state updates drawn shape records                        */
/* ========================================================================== */
static int
test_highlight_state(const char *datadir)
{
    bu_log("\n--- Test 1: highlight state updates drawn shape records ---\n");

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", datadir);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmp("mged_view_state_t1.g", std::ios::binary);
    tmp << orig.rdbuf();
    orig.close(); tmp.close();
    bu_vls_free(&fname);

    struct ged *gedp = open_gedp("mged_view_state_t1.g", 512, 512);
    if (!gedp) { bu_file_delete("mged_view_state_t1.g"); return 1; }

    const char *s_av[8] = {NULL};
    s_av[0] = "draw"; s_av[1] = "all.g"; s_av[2] = NULL;
    ged_exec_draw(gedp, 2, s_av);
    s_av[0] = "autoview"; s_av[1] = NULL;
    ged_exec_autoview(gedp, 1, s_av);

    do_full_refresh(gedp);

    int fail = 0;
    uint64_t highlight_rev0 = ged_draw_highlight_revision(gedp);
    fail += capture_screengrab_nonempty(gedp,
	    "mged_view_state_t1_base.png", "initial highlight state");

    ged_draw_set_highlight_state(gedp, 1);
    do_full_refresh(gedp);
    fail += capture_screengrab_nonempty(gedp,
	    "mged_view_state_t1_highlight.png", "highlighted state");

    uint64_t highlight_rev1 = ged_draw_highlight_revision(gedp);
    if (highlight_rev1 <= highlight_rev0) {
	bu_log("FAIL: highlight revision did not advance after highlight set\n");
	fail = 1;
    } else {
	bu_log("PASS: highlight revision advanced after highlight set\n");
    }

    int highlight_diff = draw_test_images_differ(
	    "mged_view_state_t1_base.png",
	    "mged_view_state_t1_highlight.png", 10);
    if (highlight_diff < 0) {
	bu_log("FAIL: could not compare initial and highlighted Obol captures\n");
	fail = 1;
    } else if (!highlight_diff) {
	bu_log("FAIL: highlighted Obol capture did not visibly differ from initial capture\n");
	fail = 1;
    } else {
	bu_log("PASS: highlighted Obol capture visibly differs from initial capture\n");
    }

    ged_draw_set_highlight_state(gedp, 0);
    do_full_refresh(gedp);
    fail += capture_screengrab_nonempty(gedp,
	    "mged_view_state_t1_clear.png", "cleared highlight state");

    uint64_t highlight_rev2 = ged_draw_highlight_revision(gedp);
    if (highlight_rev2 <= highlight_rev1) {
	bu_log("FAIL: highlight revision did not advance after highlight clear\n");
	fail = 1;
    } else {
	bu_log("PASS: highlight revision advanced after highlight clear\n");
    }

    int clear_diff = draw_test_images_differ(
	    "mged_view_state_t1_base.png",
	    "mged_view_state_t1_clear.png", 10);
    if (clear_diff < 0) {
	bu_log("FAIL: could not compare initial and cleared Obol captures\n");
	fail = 1;
    } else if (clear_diff) {
	bu_log("FAIL: cleared Obol capture differs from initial capture\n");
	fail = 1;
    } else {
	bu_log("PASS: cleared Obol capture matches initial capture\n");
    }

    bu_file_delete("mged_view_state_t1_base.png");
    bu_file_delete("mged_view_state_t1_highlight.png");
    bu_file_delete("mged_view_state_t1_clear.png");
    bu_file_delete("mged_view_state_t1.g");
    ged_close(gedp);
    return fail;
}

/* ========================================================================== */
/* Test 2: librt edit view snapshots consume libbv view contexts              */
/* ========================================================================== */
static int
test_edit_context_snapshot(const char *datadir)
{
    bu_log("\n--- Test 2: edit context view snapshot ---\n");

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", datadir);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmp("mged_view_state_t2.g", std::ios::binary);
    tmp << orig.rdbuf();
    orig.close(); tmp.close();
    bu_vls_free(&fname);

    struct ged *gedp = open_gedp("mged_view_state_t2.g", 512, 512);
    if (!gedp) { bu_file_delete("mged_view_state_t2.g"); return 1; }

    const char *s_av[8] = {NULL};
    s_av[0] = "draw"; s_av[1] = "all.g"; s_av[2] = NULL;
    ged_exec_draw(gedp, 2, s_av);
    s_av[0] = "autoview"; s_av[1] = NULL;
    ged_exec_autoview(gedp, 1, s_av);

    void *v = ged_view_active_ctx(gedp);
    bv_scale_set(DRAW_TEST_BV(v), 321.0);
    bv_context_update((struct bv_context *)v, BV_CONTEXT_CHANGED_VIEW);

    int fail = 0;
    struct rt_edit_view ev;
    rt_edit_view_from_context(&ev, v);
    if (!NEAR_ZERO(ev.gv_scale - 321.0, VUNITIZE_TOL)) {
	bu_log("FAIL: rt_edit_view_from_context did not snapshot bv scale "
		"(got %g)\n", ev.gv_scale);
	fail = 1;
    } else {
	bu_log("PASS: rt_edit_view_from_context snapshots bv scale\n");
    }

    struct bn_tol tol;
    BN_TOL_INIT(&tol);
    struct rt_edit *edit = rt_edit_create_context(NULL, NULL, &tol, v);
    if (!edit || !edit->vp ||
	    !NEAR_ZERO(edit->vp->gv_scale - 321.0, VUNITIZE_TOL)) {
	bu_log("FAIL: rt_edit_create_context did not initialize edit view "
		"state from bv context\n");
	fail = 1;
    } else {
	bu_log("PASS: rt_edit_create_context initializes edit view state "
		"from bv context\n");
    }
    rt_edit_destroy(edit);

    bu_file_delete("mged_view_state_t2.g");
    ged_close(gedp);
    return fail;
}

/* ========================================================================== */
/* Test 3: multiple Obol DMs stay attached to their own view controllers       */
/* ========================================================================== */
static int
test_multi_obol_dm_attachment(const char *datadir)
{
    bu_log("\n--- Test 3: multi-view Obol DM attachment ---\n");

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", datadir);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmp("mged_view_state_t3.g", std::ios::binary);
    tmp << orig.rdbuf();
    orig.close(); tmp.close();
    bu_vls_free(&fname);

    struct ged *gedp = open_gedp("mged_view_state_t3.g", 384, 384);
    if (!gedp) {
	bu_file_delete("mged_view_state_t3.g");
	return 1;
    }

    int fail = 0;
    void *v0 = ged_view_active_ctx(gedp);
    bv_name_set(DRAW_TEST_BV(v0), "V0");

    void *view_set_ctx = ged_view_set_ctx(gedp);
    void *v1 = ged_view_context_create_with_set(view_set_ctx);
    if (!v1) {
	bu_log("FAIL: secondary view creation failed\n");
	fail = 1;
    } else {
	bv_name_set(DRAW_TEST_BV(v1), "V1");
	bv_unit_conversion_set(DRAW_TEST_BV(v1), gedp->dbip->dbi_local2base, gedp->dbip->dbi_base2local);
	ged_view_set_context_add(view_set_ctx, v1);
	ged_view_context_owned_add(gedp, v1);

	const char *dm_av[7] = {"dm", "attach", "-V", "V1", "obol",
	    "OBOL1", NULL};
	if (ged_exec_dm(gedp, 6, dm_av) != BRLCAD_OK ||
		!configure_obol_view(gedp, v1, 384, 384)) {
	    bu_log("FAIL: secondary Obol DM attachment failed: %s\n",
		    bu_vls_cstr(gedp->ged_result_str));
	    fail = 1;
	}
    }

    BRLObolViewController *v0_controller = obol_controller_for_view(v0);
    BRLObolViewController *v1_controller = obol_controller_for_view(v1);
    if (!v0_controller || !v1_controller ||
	v0_controller == v1_controller) {
	bu_log("FAIL: views should have distinct Obol controllers\n");
	fail = 1;
    }
    if (!fail) {
	const char *fast_av[5] = {"dm", "set", "software_wire", "fast", NULL};
	const int fast_ret = ged_exec_dm(gedp, 4, fast_av);
	if (fast_ret != BRLCAD_OK ||
		v0_controller->getSoftwareWireMode() !=
		    BRLObolViewController::SOFTWARE_WIRE_FAST ||
		v1_controller->getSoftwareWireMode() !=
		    BRLObolViewController::SOFTWARE_WIRE_AUTO) {
	    bu_log("FAIL: dm set software_wire fast did not update only the active "
		    "view (ret=%d v0=%d v1=%d result=%s)\n", fast_ret,
		    (int)v0_controller->getSoftwareWireMode(),
		    (int)v1_controller->getSoftwareWireMode(),
		    bu_vls_cstr(gedp->ged_result_str));
	    fail = 1;
	}
	const char *query_av[4] = {"dm", "set", "software_wire", NULL};
	if (!fail && (ged_exec_dm(gedp, 3, query_av) != BRLCAD_OK ||
		!BU_STR_EQUAL(bu_vls_cstr(gedp->ged_result_str), "fast"))) {
	    bu_log("FAIL: dm set software_wire did not report fast: %s\n",
		    bu_vls_cstr(gedp->ged_result_str));
	    fail = 1;
	}
	const char *named_av[7] = {
	    "dm", "set", "-d", "OBOL1", "software_wire", "quality", NULL
	};
	if (!fail && (ged_exec_dm(gedp, 6, named_av) != BRLCAD_OK ||
		v0_controller->getSoftwareWireMode() !=
		    BRLObolViewController::SOFTWARE_WIRE_FAST ||
		v1_controller->getSoftwareWireMode() !=
		    BRLObolViewController::SOFTWARE_WIRE_QUALITY)) {
	    bu_log("FAIL: named dm software_wire setting was not view-local\n");
	    fail = 1;
	}
	const char *bad_av[5] = {"dm", "set", "software_wire", "turbo", NULL};
	if (!fail && ged_exec_dm(gedp, 4, bad_av) != BRLCAD_ERROR) {
	    bu_log("FAIL: dm set software_wire accepted an invalid mode\n");
	    fail = 1;
	}
	if (!fail)
	    bu_log("PASS: dm software_wire modes are queryable and view-local\n");
    }
    if (!fail) {
	BRLObolViewAttachment *v0_attachment =
	    ged_view_context_obol_attachment(v0);
	BRLObolViewAttachment *v1_attachment =
	    ged_view_context_obol_attachment(v1);
	if (!v0_attachment || !v1_attachment ||
		v0_attachment != v0_controller->getViewAttachment() ||
		v1_attachment != v1_controller->getViewAttachment() ||
		v0_attachment == v1_attachment) {
	    bu_log("FAIL: GED view contexts should share attachment state "
		    "with their attached Obol controllers\n");
	    fail = 1;
	} else {
	    ged_draw_view_lod_policy lod_policy = BV_LOD_POLICY_INIT;
	    lod_policy.csg_enabled = 1;
	    lod_policy.bot_threshold = 42;
	    lod_policy.scale = 2.5;
	    if (!ged_draw_view_context_lod_policy_apply(v1, &lod_policy)) {
		bu_log("FAIL: applying per-view LoD policy should succeed\n");
		fail = 1;
	    } else {
		struct bv_lod_policy controller_policy = BV_LOD_POLICY_INIT;
		v1_controller->getViewAttachment()->getLodPolicy(
		    &controller_policy);
		if (!controller_policy.csg_enabled ||
			controller_policy.bot_threshold != 42 ||
			!NEAR_ZERO(controller_policy.scale - 2.5,
			    VUNITIZE_TOL)) {
		    bu_log("FAIL: GED LoD policy updates should be visible "
			    "through the attached controller attachment "
			    "(csg=%d mesh=%d threshold=%zu scale=%g)\n",
			    controller_policy.csg_enabled,
			    controller_policy.mesh_enabled,
			    controller_policy.bot_threshold,
			    controller_policy.scale);
		    fail = 1;
		}
	    }
	}
    }

    const char *draw_av[3] = {"draw", "all.g", NULL};
    if (!fail && ged_exec_draw(gedp, 2, draw_av) != BRLCAD_OK) {
	bu_log("FAIL: shared draw failed: %s\n",
		bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }

    do_full_refresh(gedp);

    SoBRLSceneController *shared_scene = ged_draw_obol_scene_controller(gedp);
    int shared_sources = shared_scene ?
	shared_scene->getDatabaseSourceCount() : 0;
    if (shared_sources <= 0) {
	bu_log("FAIL: GED Obol scene should own shared draw state "
		"(sources=%d)\n", shared_sources);
	fail = 1;
    } else if (!v0_controller->getRenderSceneRoot() ||
	    !v1_controller->getRenderSceneRoot()) {
	bu_log("FAIL: Obol DMs should have per-view render scene roots\n");
	fail = 1;
    } else {
	bu_log("PASS: GED Obol scene owns shared draw state and both DMs "
		"render through per-view roots (sources=%d)\n",
		shared_sources);
    }

    unsigned char *bridge_image = NULL;
    int bridge_ret = ged_draw_obol_view_display_image(gedp, v0,
	    &bridge_image, 1, 0);
    if (bridge_ret != 1 || !bridge_image) {
	bu_log("FAIL: GED Obol display-image bridge should render an "
		"attached Obol view (ret=%d)\n", bridge_ret);
	fail = 1;
    } else {
	struct dm *v0_dmp =
	    (struct dm *)ged_view_context_display_manager_get(v0);
	size_t pixel_count = (size_t)dm_get_width(v0_dmp) *
	    (size_t)dm_get_height(v0_dmp);
	size_t lit_pixels = 0;
	for (size_t i = 0; i < pixel_count; i++) {
	    const unsigned char *p = bridge_image + i * 3;
	    if (p[0] || p[1] || p[2])
		lit_pixels++;
	}
	if (!lit_pixels) {
	    bu_log("FAIL: GED Obol display-image bridge returned an "
		    "empty image\n");
	    fail = 1;
	} else {
	    bu_log("PASS: GED Obol display-image bridge rendered %zu "
		    "non-black pixels\n", lit_pixels);
	}
    }
    if (bridge_image)
	bu_free(bridge_image, "GED Obol display-image bridge test image");

    ged_close(gedp);
    bu_file_delete("mged_view_state_t3.g");
    return fail;
}

/* Test 4: Obol DM attach clears stale standalone rt framebuffer cache        */
/* ========================================================================== */
static int
test_obol_rt_framebuffer_cache(const char *datadir)
{
    bu_log("\n--- Test 4: Obol DM attach clears stale rt framebuffer cache ---\n");

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", datadir);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmp("mged_view_state_t4.g", std::ios::binary);
    tmp << orig.rdbuf();
    orig.close(); tmp.close();
    bu_vls_free(&fname);

    struct ged *gedp = ged_open("db", "mged_view_state_t4.g", 1);
    if (!gedp) {
	bu_file_delete("mged_view_state_t4.g");
	return 1;
    }

    int fail = 0;
    ged_rt_fb_set(gedp, "/dev/swrast");
    const char *stale_fb = ged_rt_fb_get(gedp);
    if (!stale_fb || !BU_STR_EQUAL(stale_fb, "/dev/swrast")) {
	bu_log("FAIL: test setup could not seed stale rt framebuffer device\n");
	fail = 1;
    }

    const char *dm_av[5] = {"dm", "attach", "obol", "OBOLRTFB", NULL};
    if (!fail && ged_exec_dm(gedp, 4, dm_av) != BRLCAD_OK) {
	bu_log("FAIL: could not attach Obol DM for rt framebuffer cache test: %s\n",
		bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }

    const char *fb_after_attach = ged_rt_fb_get(gedp);
    if (!fail && fb_after_attach && fb_after_attach[0]) {
	bu_log("FAIL: Obol DM attach should clear stale rt framebuffer "
		"device, got %s\n", fb_after_attach);
	fail = 1;
    } else if (!fail) {
	bu_log("PASS: Obol DM attach cleared stale rt framebuffer device\n");
    }

    ged_rt_fb_set(gedp, "/dev/qtgl");
    ged_rt_fb_refresh(gedp);
    const char *fb_after_refresh = ged_rt_fb_get(gedp);
    if (!fail && fb_after_refresh && fb_after_refresh[0]) {
	bu_log("FAIL: rt framebuffer refresh on Obol DM should clear "
		"standalone device, got %s\n", fb_after_refresh);
	fail = 1;
    } else if (!fail) {
	bu_log("PASS: rt framebuffer refresh on Obol DM keeps standalone "
		"device clear\n");
    }

    ged_close(gedp);
    bu_file_delete("mged_view_state_t4.g");
    return fail;
}

/* ========================================================================== */
/* Test 5: non-Obol GL hosts can render through an owned Obol view endpoint   */
/* ========================================================================== */
static int
test_owned_render_endpoint(const char *datadir)
{
    bu_log("\n--- Test 5: owned Obol render endpoint ---\n");

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", datadir);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmp("mged_view_state_t5.g", std::ios::binary);
    tmp << orig.rdbuf();
    orig.close(); tmp.close();
    bu_vls_free(&fname);

    struct ged *gedp = ged_open("db", "mged_view_state_t5.g", 1);
    if (!gedp) {
	bu_file_delete("mged_view_state_t5.g");
	return 1;
    }

    int fail = 0;
    void *view_ctx = ged_view_active_ctx(gedp);
    bv_dimensions_set(DRAW_TEST_BV(view_ctx), 384, 384);
    bv_unit_conversion_set(DRAW_TEST_BV(view_ctx),
	gedp->dbip->dbi_local2base, gedp->dbip->dbi_base2local);

    if (!ged_draw_obol_render_endpoint_ensure_for_view(gedp, view_ctx, 1)) {
	bu_log("FAIL: could not create an owned Obol render endpoint\n");
	fail = 1;
    }

    BRLObolViewController *controller = obol_controller_for_view(view_ctx);
    if (!fail && (!controller || !controller->getRenderSceneRoot())) {
	bu_log("FAIL: owned endpoint did not bind a per-view render root\n");
	fail = 1;
    }

    const char *av[3] = {"draw", "all.g", NULL};
    if (!fail && ged_exec_draw(gedp, 2, av) != BRLCAD_OK) {
	bu_log("FAIL: owned-endpoint draw failed: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    av[0] = "autoview";
    if (!fail && ged_exec_autoview(gedp, 1, av) != BRLCAD_OK) {
	bu_log("FAIL: owned-endpoint autoview failed: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    if (!fail)
	do_full_refresh(gedp);

    unsigned char *image = NULL;
    if (!fail) {
	controller->setViewportSize(384, 384);
	controller->syncCameraFromViewContext(view_ctx);
	int ret = controller->renderToImage(&image, 1, 0);
	if (ret != BRLCAD_OK || !image) {
	    bu_log("FAIL: owned endpoint did not render image data\n");
	    fail = 1;
	}
    }

    if (!fail) {
	size_t lit_pixels = 0;
	for (size_t i = 0; i < 384u * 384u; i++) {
	    const unsigned char *p = image + i * 3;
	    if (p[0] || p[1] || p[2])
		lit_pixels++;
	}
	if (!lit_pixels) {
	    bu_log("FAIL: owned endpoint rendered an empty image\n");
	    fail = 1;
	} else {
	    bu_log("PASS: owned endpoint rendered %zu non-black pixels\n",
		lit_pixels);
	}
    }

    unsigned char *deferred_root_image = NULL;
    unsigned char *deferred_full_image = NULL;
    if (!fail) {
	const char *erase_av[2] = {"erase", "all.g"};
	const char *deferred_av[3] = {
	    "draw", "--defer-leaf-expansion", "all.g"
	};
	if (ged_exec_erase(gedp, 2, erase_av) != BRLCAD_OK ||
	    ged_exec_draw(gedp, 3, deferred_av) != BRLCAD_OK) {
	    bu_log("FAIL: owned-endpoint deferred draw failed: %s\n",
		bu_vls_cstr(gedp->ged_result_str));
	    fail = 1;
	}
    }

    BRLObolProgressiveStatus root_status;
    if (!fail &&
	(controller->renderToImage(&deferred_root_image, 1, 0, NULL, NULL,
	    &root_status) != BRLCAD_OK || !deferred_root_image ||
	 !root_status.hasMore || root_status.changed ||
	 memcmp(image, deferred_root_image, 384u * 384u * 3u) == 0)) {
	bu_log("FAIL: deferred root frame should be distinct and request refinement\n");
	fail = 1;
    }

    BRLObolProgressiveStatus full_status;
    if (!fail &&
	(controller->renderToImage(&deferred_full_image, 1, 0, NULL, NULL,
	    &full_status) != BRLCAD_OK || !deferred_full_image ||
	 full_status.hasMore || !full_status.changed ||
	 memcmp(image, deferred_full_image, 384u * 384u * 3u) != 0)) {
	bu_log("FAIL: deferred refinement should reproduce the direct image exactly\n");
	fail = 1;
    } else if (!fail) {
	bu_log("PASS: deferred root refines to the exact direct image\n");
    }

    if (image)
	bu_free(image, "owned Obol render endpoint test image");
    if (deferred_root_image)
	bu_free(deferred_root_image, "owned Obol deferred root test image");
    if (deferred_full_image)
	bu_free(deferred_full_image, "owned Obol deferred full test image");
    ged_close(gedp);
    bu_file_delete("mged_view_state_t5.g");
    return fail;
}

/* main                                                                        */
/* ========================================================================== */

int
main(int argc, char *argv[])
{
    bu_setprogname(argv[0]);

    int cleanup = 0;
    struct bu_opt_desc d[3];
    BU_OPT(d[0], "c", "cleanup", NULL, NULL, &cleanup, "cleanup temp files");
    BU_OPT_NULL(d[1]);

    int uac = bu_opt_parse(NULL, argc, (const char **)argv, d);
    if (uac != 2) {
	bu_log("Usage: %s [-c] <directory-containing-moss.g>\n", argv[0]);
	return 1;
    }
    const char *datadir = argv[1];

    int failures = 0;
    failures += test_highlight_state(datadir);
    failures += test_edit_context_snapshot(datadir);
    failures += test_multi_obol_dm_attachment(datadir);
    failures += test_obol_rt_framebuffer_cache(datadir);
    failures += test_owned_render_endpoint(datadir);
    failures += test_owned_render_endpoint(datadir);

    if (failures == 0) {
	bu_log("\nAll MGED view-state tests PASSED (%d/6)\n", 6);
    } else {
	bu_log("\n%d MGED view-state test(s) FAILED\n", failures);
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
