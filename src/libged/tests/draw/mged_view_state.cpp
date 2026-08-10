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
 *   6. The session framebuffer capture provider moves between endpoints.
 *
 * Uses the headless Obol/Coin off-screen renderer; no display hardware
 * required.
 *
 * Usage: ged_test_mged_view_state [-c] <directory-containing-moss.g>
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include <cmath>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include <bu.h>
#include <icv.h>
#include "bu/opt.h"
#include "rt/edit.h"
#include "rt/view.h"
#include "view_test_util.h"
#include <ged.h>
#include "ged/draw.h"
#include "ged/display.h"
#include "BObol/BDisplayEndpoint.h"
#include "BObol/BInit.h"
#include "BObol/BViewAttachment.h"
#include "BObol/BSceneController.h"
#include "BObol/BViewController.h"
#include <Inventor/nodes/SoGroup.h>
#include "../../ged_private.h"
#include "../../ged_draw_view_private.h"

extern "C" void ged_changed_callback(struct db_i *UNUSED(dbip), struct directory *dp, int mode, void *u_data);
extern "C" long draw_test_count_nonblack_pixels(const char *filename);
extern "C" int draw_test_images_differ(const char *a, const char *b, int offmany_threshold);

/* Full redraw through the neutral GED draw transaction path. */
static void
do_full_refresh(struct ged *gedp)
{
    struct ged_view_context *v = ged_view_active_ctx(gedp);
    struct ged_scene_redraw_request request;
    ged_scene_redraw_request_init(&request);
    request.view = v;
    (void)ged_scene_redraw(gedp, &request, NULL);
}

static double
rgb_image_ssim(const unsigned char *a, const unsigned char *b,
	size_t width, size_t height)
{
    if (!a || !b || !width || !height)
	return -1.0;

    icv_image_t *ia = icv_create(width, height, ICV_COLOR_SPACE_RGB);
    icv_image_t *ib = icv_create(width, height, ICV_COLOR_SPACE_RGB);
    if (!ia || !ib) {
	if (ia)
	    icv_destroy(ia);
	if (ib)
	    icv_destroy(ib);
	return -1.0;
    }

    const size_t channel_count = width * height * 3u;
    for (size_t i = 0; i < channel_count; i++) {
	ia->data[i] = static_cast<double>(a[i]) / 255.0;
	ib->data[i] = static_cast<double>(b[i]) / 255.0;
    }
    const double score = icv_adiff(ia, ib, ICV_DIFF_SSIM);
    icv_destroy(ia);
    icv_destroy(ib);
    return score;
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

static BObolViewController *
obol_controller_for_view(struct ged_view_context *view_ctx)
{
    bobol_display_endpoint_t *endpoint =
	ged_view_context_obol_endpoint_get(view_ctx);
    return endpoint ? static_cast<BObolViewController *>(
	bobol_display_endpoint_controller(endpoint)) : NULL;
}

static int
obol_interlay_is_between_scene_and_local_root(
	BObolViewController *controller, SoNode *shared_root)
{
    if (!controller || !shared_root)
	return 0;
    SoNode *lod_root = controller->getRenderSceneRoot();
    if (!lod_root || !lod_root->isOfType(SoGroup::getClassTypeId()))
	return 0;
    SoGroup *lod_group = static_cast<SoGroup *>(lod_root);
    if (lod_group->getNumChildren() != 1)
	return 0;
    SoNode *composition_node = lod_group->getChild(0);
    if (!composition_node ||
	!composition_node->isOfType(SoGroup::getClassTypeId()))
	return 0;
    SoGroup *composition = static_cast<SoGroup *>(composition_node);
    return composition->getNumChildren() == 3 &&
	composition->getChild(0) == shared_root &&
	composition->getChild(1) == controller->getFramebufferInterlayRoot() &&
	composition->getChild(2) == controller->getSceneRoot();
}

static int
configure_obol_view(struct ged *gedp, struct ged_view_context *view_ctx, int width, int height)
{
    bobol_display_endpoint_t *endpoint =
	ged_view_context_obol_endpoint_get(view_ctx);
    if (!endpoint)
	return 0;

    bv_dimensions_set(DRAW_TEST_BV(view_ctx), width, height);
    bv_unit_conversion_set(DRAW_TEST_BV(view_ctx), gedp->dbip->dbi_local2base, gedp->dbip->dbi_base2local);
    return bobol_display_endpoint_resize(endpoint, (unsigned int)width,
	    (unsigned int)height, 1.0);
}

/* ---- common GED + Obol setup --------------------------------------------- */

static struct ged *
open_gedp(const char *gfile, int width, int height)
{
    struct ged *gedp = ged_open("db", gfile, 1);
    if (!gedp) return NULL;

    db_add_changed_clbk(gedp->dbip, &ged_changed_callback, (void *)gedp);

    /* Open a headless Obol endpoint without a compatibility DM. */
    const char *s_av[16] = {NULL};
    struct ged_view_context *v = ged_view_active_ctx(gedp);
    bv_dimensions_set(DRAW_TEST_BV(v), width, height);
    s_av[0] = "dm"; s_av[1] = "open"; s_av[2] = "--host";
    s_av[3] = "headless"; s_av[4] = "--renderer"; s_av[5] = "sw";
    s_av[6] = NULL;
    if (ged_exec_dm(gedp, 6, s_av) != BRLCAD_OK) {
	ged_close(gedp);
	return NULL;
    }
    if (!obol_controller_for_view(v)) {
	bu_log("FAIL: Obol endpoint open did not associate a GED view controller\n");
	ged_close(gedp);
	return NULL;
    }
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
    /* Establish the baseline only after the asynchronous source realization
     * settles.  Highlighting itself requests realization, so comparing its
     * final frame to an initial proxy frame would test scheduling, not style. */
    if (!draw_test_obol_progressive_drain(gedp, ged_view_active_ctx(gedp),
	    500, 1)) {
	bu_log("FAIL: initial Obol realization did not settle for highlight test\n");
	fail = 1;
    }
    uint64_t highlight_rev0 = ged_scene_highlight_revision(gedp);
    fail += capture_screengrab_nonempty(gedp,
	    "mged_view_state_t1_base.png", "initial highlight state");

    struct ged_scene_path_request highlight_request;
    ged_scene_path_request_init(&highlight_request);
    highlight_request.path = "box.s";
    highlight_request.match = GED_SCENE_PATH_MATCH_OBJECT;
    (void)ged_scene_highlight_set(gedp, &highlight_request, 1, NULL);
    do_full_refresh(gedp);
    fail += capture_screengrab_nonempty(gedp,
	    "mged_view_state_t1_highlight.png", "highlighted state");

    uint64_t highlight_rev1 = ged_scene_highlight_revision(gedp);
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

    (void)ged_scene_highlights_clear(gedp, NULL);
    do_full_refresh(gedp);
    fail += capture_screengrab_nonempty(gedp,
	    "mged_view_state_t1_clear.png", "cleared highlight state");

    uint64_t highlight_rev2 = ged_scene_highlight_revision(gedp);
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

    struct ged_view_context *v = ged_view_active_ctx(gedp);
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
    struct ged_view_context *v0 = ged_view_active_ctx(gedp);
    bv_name_set(DRAW_TEST_BV(v0), "V0");

    struct ged_view_set *view_set_ctx = ged_view_set_ctx(gedp);
    struct ged_view_context *v1 =
	ged_view_context_create_with_set(view_set_ctx);
    if (!v1) {
	bu_log("FAIL: secondary view creation failed\n");
	fail = 1;
    } else {
	bv_name_set(DRAW_TEST_BV(v1), "V1");
	bv_unit_conversion_set(DRAW_TEST_BV(v1), gedp->dbip->dbi_local2base, gedp->dbip->dbi_base2local);
	ged_view_set_context_add(view_set_ctx, v1);
	ged_view_context_owned_add(gedp, v1);

	bv_dimensions_set(DRAW_TEST_BV(v1), 384, 384);
	const char *dm_av[9] = {"dm", "open", "-V", "V1", "--host",
	    "headless", "--renderer", "sw", NULL};
	if (ged_exec_dm(gedp, 8, dm_av) != BRLCAD_OK ||
		!configure_obol_view(gedp, v1, 384, 384)) {
	    bu_log("FAIL: secondary Obol DM attachment failed: %s\n",
		    bu_vls_cstr(gedp->ged_result_str));
	    fail = 1;
	}
    }

    BObolViewController *v0_controller = obol_controller_for_view(v0);
    BObolViewController *v1_controller = obol_controller_for_view(v1);
    if (!v0_controller || !v1_controller ||
	v0_controller == v1_controller) {
	bu_log("FAIL: views should have distinct Obol controllers\n");
	fail = 1;
    }
    if (!fail) {
	const char *fast_av[5] = {
	    "dm", "set", "controller.software_wire", "fast", NULL
	};
	const int fast_ret = ged_exec_dm(gedp, 4, fast_av);
	if (fast_ret != BRLCAD_OK ||
		v0_controller->getSoftwareWireMode() !=
		    BObolViewController::SOFTWARE_WIRE_FAST ||
		v1_controller->getSoftwareWireMode() !=
		    BObolViewController::SOFTWARE_WIRE_AUTO) {
	    bu_log("FAIL: dm set controller.software_wire fast did not update only the active "
		    "view (ret=%d v0=%d v1=%d result=%s)\n", fast_ret,
		    (int)v0_controller->getSoftwareWireMode(),
		    (int)v1_controller->getSoftwareWireMode(),
		    bu_vls_cstr(gedp->ged_result_str));
	    fail = 1;
	}
	const char *query_av[4] = {
	    "dm", "get", "controller.software_wire", NULL
	};
	if (!fail && (ged_exec_dm(gedp, 3, query_av) != BRLCAD_OK ||
		!BU_STR_EQUAL(bu_vls_cstr(gedp->ged_result_str), "fast"))) {
	    bu_log("FAIL: dm get controller.software_wire did not report fast: %s\n",
		    bu_vls_cstr(gedp->ged_result_str));
	    fail = 1;
	}
	const char *named_av[7] = {
	    "dm", "set", "-V", "V1", "controller.software_wire", "quality", NULL
	};
	if (!fail && (ged_exec_dm(gedp, 6, named_av) != BRLCAD_OK ||
		v0_controller->getSoftwareWireMode() !=
		    BObolViewController::SOFTWARE_WIRE_FAST ||
		v1_controller->getSoftwareWireMode() !=
		    BObolViewController::SOFTWARE_WIRE_QUALITY)) {
	    bu_log("FAIL: named dm controller.software_wire setting was not view-local\n");
	    fail = 1;
	}
	const char *bad_av[5] = {
	    "dm", "set", "controller.software_wire", "turbo", NULL
	};
	if (!fail && ged_exec_dm(gedp, 4, bad_av) != BRLCAD_ERROR) {
	    bu_log("FAIL: dm set controller.software_wire accepted an invalid mode\n");
	    fail = 1;
	}
	if (!fail)
	    bu_log("PASS: typed software-wire modes are queryable and view-local\n");
    }
    if (!fail) {
	BObolViewAttachment *v0_attachment =
	    ged_view_context_obol_attachment(v0);
	BObolViewAttachment *v1_attachment =
	    ged_view_context_obol_attachment(v1);
	if (!v0_attachment || !v1_attachment ||
		v0_attachment != v0_controller->getViewAttachment() ||
		v1_attachment != v1_controller->getViewAttachment() ||
		v0_attachment == v1_attachment) {
	    bu_log("FAIL: GED view contexts should share attachment state "
		    "with their attached Obol controllers\n");
	    fail = 1;
	} else {
	    ged_view_lod_policy lod_policy = BV_LOD_POLICY_INIT;
	    lod_policy.csg_enabled = 1;
	    lod_policy.bot_threshold = 42;
	    lod_policy.scale = 2.5;
	    if (!ged_view_lod_policy_apply(v1, &lod_policy)) {
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
    if (!fail && !draw_test_obol_progressive_drain(gedp, v0, 500, 1)) {
	bu_log("FAIL: shared progressive draw did not settle before image capture\n");
	fail = 1;
    }

    do_full_refresh(gedp);

    SoNode *shared_root = NULL;
    SoNode *v0_render_root = v0_controller->getRenderSceneRoot();
    if (v0_render_root && v0_render_root->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *render_group = static_cast<SoGroup *>(v0_render_root);
	SoNode *composition_node = render_group->getNumChildren() == 1 ?
	    render_group->getChild(0) : NULL;
	if (composition_node &&
	    composition_node->isOfType(SoGroup::getClassTypeId())) {
	    SoGroup *composition = static_cast<SoGroup *>(composition_node);
	    if (composition->getNumChildren() > 0)
		shared_root = composition->getChild(0);
	}
    }
    BObolSceneController shared_scene(shared_root);
    int shared_sources = shared_scene.getDatabaseSourceCount();
    if (shared_sources <= 0) {
	bu_log("FAIL: GED Obol scene should own shared draw state "
		"(sources=%d)\n", shared_sources);
	fail = 1;
    } else if (!v0_controller->getRenderSceneRoot() ||
	    !v1_controller->getRenderSceneRoot()) {
	bu_log("FAIL: Obol DMs should have per-view render scene roots\n");
	fail = 1;
	} else if (!obol_interlay_is_between_scene_and_local_root(v0_controller,
		shared_root) ||
	    !obol_interlay_is_between_scene_and_local_root(v1_controller,
		shared_root)) {
	bu_log("FAIL: Obol framebuffer interlay should separate shared geometry "
		"from view-local features\n");
	fail = 1;
    } else {
	bu_log("PASS: GED Obol scene owns shared draw state and both DMs "
		"render through per-view roots (sources=%d)\n",
		shared_sources);
    }

    unsigned char *bridge_image = NULL;
    int bridge_ret = ged_view_display_image_capture(gedp, v0,
	    &bridge_image, 1, 0);
    if (bridge_ret != 1 || !bridge_image) {
	bu_log("FAIL: GED Obol display-image bridge should render an "
		"attached Obol view (ret=%d)\n", bridge_ret);
	fail = 1;
    } else {
	size_t pixel_count = (size_t)bv_width_get(DRAW_TEST_BV(v0)) *
	    (size_t)bv_height_get(DRAW_TEST_BV(v0));
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

/* Test 4: non-Obol GL hosts can render through an owned Obol view endpoint   */
/* ========================================================================== */
static int
test_owned_render_endpoint(const char *datadir)
{
    bu_log("\n--- Test 4: owned Obol render endpoint ---\n");

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
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    bv_dimensions_set(DRAW_TEST_BV(view_ctx), 384, 384);
    bv_unit_conversion_set(DRAW_TEST_BV(view_ctx),
	gedp->dbip->dbi_local2base, gedp->dbip->dbi_base2local);

    if (!ged_view_context_display_endpoint_ensure(view_ctx)) {
	bu_log("FAIL: could not create an owned Obol render endpoint\n");
	fail = 1;
    }

    BObolViewController *controller = obol_controller_for_view(view_ctx);
    if (!fail && (!controller || !controller->getRenderSceneRoot())) {
	bu_log("FAIL: owned endpoint did not bind a per-view render root\n");
	fail = 1;
    }
    if (!fail)
	controller->setRenderContextManager(bobol_headless_context_manager());

    const char *mged_lighting_av[5] = {
	"view", "lighting", "profile", "mged", NULL
    };
    struct bv_lighting_state lighting_state;
    if (!fail && (ged_exec(gedp, 4, mged_lighting_av) != BRLCAD_OK ||
	!bv_lighting_state_get(&lighting_state, DRAW_TEST_BV(view_ctx)) ||
	lighting_state.profile != BV_LIGHTING_MGED ||
	controller->getLightingProfile() != BObolViewController::LIGHTING_MGED)) {
	bu_log("FAIL: view lighting profile did not select shared MGED policy: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    const char *studio_lighting_av[5] = {
	"view", "lighting", "profile", "studio", NULL
    };
    if (!fail && (ged_exec(gedp, 4, studio_lighting_av) != BRLCAD_OK ||
	!bv_lighting_state_get(&lighting_state, DRAW_TEST_BV(view_ctx)) ||
	lighting_state.profile != BV_LIGHTING_STUDIO ||
	controller->getLightingProfile() != BObolViewController::LIGHTING_STUDIO)) {
	bu_log("FAIL: view lighting profile did not restore shared studio policy: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }

    const char *status_av[3] = {"dm", "status", NULL};
    if (!fail && (ged_exec_dm(gedp, 2, status_av) != BRLCAD_OK ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str), "renderer=auto") ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str), "resolved=auto") ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str), "width="))) {
	bu_log("FAIL: dm status could not report a DM-less endpoint: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }

    const char *get_renderer_av[4] = {
	"dm", "get", "endpoint.renderer", NULL
    };
    if (!fail && (ged_exec_dm(gedp, 3, get_renderer_av) != BRLCAD_OK ||
	    !BU_STR_EQUAL(bu_vls_cstr(gedp->ged_result_str), "auto"))) {
	bu_log("FAIL: dm get did not read a typed endpoint property: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }

    const char *set_bg_av[5] = {
	"dm", "set", "controller.background.bottom", "0.1/0.2/0.3", NULL
    };
    if (!fail && ged_exec_dm(gedp, 4, set_bg_av) != BRLCAD_OK) {
	bu_log("FAIL: dm set did not write a typed color property: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    const char *get_bg_av[4] = {
	"dm", "get", "controller.background.bottom", NULL
    };
    if (!fail && (ged_exec_dm(gedp, 3, get_bg_av) != BRLCAD_OK ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str), "0.1"))) {
	bu_log("FAIL: typed background property did not round trip: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }

    const char *set_headlight_color_av[5] = {
	"dm", "set", "renderer.headlight.color", "0.25/0.5/0.75", NULL
    };
    if (!fail && (ged_exec_dm(gedp, 4, set_headlight_color_av) !=
	BRLCAD_OK || controller->getHeadlightColor() !=
	SbColor(0.25f, 0.5f, 0.75f))) {
	bu_log("FAIL: dm set did not update typed headlight color: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    const char *set_headlight_intensity_av[5] = {
	"dm", "set", "renderer.headlight.intensity", "0.4", NULL
    };
    if (!fail && (ged_exec_dm(gedp, 4, set_headlight_intensity_av) !=
	BRLCAD_OK || std::fabs(controller->getHeadlightIntensity() - 0.4f) >
	1.0e-6f)) {
	bu_log("FAIL: dm set did not update typed headlight intensity: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    const char *get_headlight_color_av[4] = {
	"dm", "get", "renderer.headlight.color", NULL
    };
    if (!fail && (ged_exec_dm(gedp, 3, get_headlight_color_av) !=
	BRLCAD_OK || !strstr(bu_vls_cstr(gedp->ged_result_str), "0.25"))) {
	bu_log("FAIL: typed headlight color did not round trip: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }

    if (!fail && ged_view_framebuffer_backend_ensure(gedp,
	    view_ctx) != BRLCAD_OK) {
	bu_log("FAIL: could not create an endpoint-backed framebuffer stream\n");
	fail = 1;
    }
    const char *set_underlay_av[5] = {
	"dm", "set", "composition.framebuffer.mode", "underlay", NULL
    };
    if (!fail && (ged_exec_dm(gedp, 4, set_underlay_av) != BRLCAD_OK ||
	    bv_framebuffer_mode_get(DRAW_TEST_BV(view_ctx)) !=
	    BV_FRAMEBUFFER_MODE_UNDERLAY)) {
	bu_log("FAIL: typed framebuffer composition did not select underlay: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    const char *get_framebuffer_mode_av[4] = {
	"dm", "get", "composition.framebuffer.mode", NULL
    };
    if (!fail && (ged_exec_dm(gedp, 3, get_framebuffer_mode_av) != BRLCAD_OK ||
	    !BU_STR_EQUAL(bu_vls_cstr(gedp->ged_result_str), "underlay"))) {
	bu_log("FAIL: typed framebuffer composition did not round trip: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    set_underlay_av[3] = "interlay";
    if (!fail && (ged_exec_dm(gedp, 4, set_underlay_av) != BRLCAD_OK ||
	    bv_framebuffer_mode_get(DRAW_TEST_BV(view_ctx)) !=
	    BV_FRAMEBUFFER_MODE_INTERLAY)) {
	bu_log("FAIL: typed framebuffer composition did not select interlay: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    if (!fail && (ged_exec_dm(gedp, 3, get_framebuffer_mode_av) != BRLCAD_OK ||
	    !BU_STR_EQUAL(bu_vls_cstr(gedp->ged_result_str), "interlay"))) {
	bu_log("FAIL: typed framebuffer composition did not round trip interlay: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    set_underlay_av[3] = "invalid";
    if (!fail && ged_exec_dm(gedp, 4, set_underlay_av) != BRLCAD_ERROR) {
	bu_log("FAIL: typed framebuffer composition accepted an invalid mode\n");
	fail = 1;
    }
    set_underlay_av[3] = "off";
    if (!fail && (ged_exec_dm(gedp, 4, set_underlay_av) != BRLCAD_OK ||
	    bv_framebuffer_mode_get(DRAW_TEST_BV(view_ctx)) !=
	    BV_FRAMEBUFFER_MODE_OFF)) {
	bu_log("FAIL: typed framebuffer composition did not disable display\n");
	fail = 1;
    }

    const char *renderer_none_av[4] = {"dm", "renderer", "none", NULL};
    if (!fail && ged_exec_dm(gedp, 3, renderer_none_av) != BRLCAD_OK) {
	bu_log("FAIL: dm renderer could not select none: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    const char *renderer_query_av[3] = {"dm", "renderer", NULL};
    if (!fail && (ged_exec_dm(gedp, 2, renderer_query_av) != BRLCAD_OK ||
	    !BU_STR_EQUAL(bu_vls_cstr(gedp->ged_result_str), "none"))) {
	bu_log("FAIL: dm renderer did not report its selected engine: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    renderer_none_av[2] = "rt";
	if (!fail && (ged_exec_dm(gedp, 3, renderer_none_av) != BRLCAD_OK ||
	    ged_exec_dm(gedp, 2, renderer_query_av) != BRLCAD_OK ||
	    !BU_STR_EQUAL(bu_vls_cstr(gedp->ged_result_str), "rt"))) {
	bu_log("FAIL: dm renderer could not select retained rt: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    renderer_none_av[2] = "auto";
    if (!fail && ged_exec_dm(gedp, 3, renderer_none_av) != BRLCAD_OK)
	fail = 1;
    if (!fail)
	bu_log("PASS: dm endpoint status, typed properties, and renderer selection\n");

    const char *wire_av[5] = {
	"dm", "set", "controller.software_wire", "fast", NULL
    };
    if (!fail && (ged_exec_dm(gedp, 4, wire_av) != BRLCAD_OK ||
	    controller->getSoftwareWireMode() !=
		BObolViewController::SOFTWARE_WIRE_FAST)) {
	bu_log("FAIL: typed software-wire property could not configure an Obol view: %s\n",
		bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    } else if (!fail) {
	bu_log("PASS: typed software-wire property configures an Obol view\n");
    }
    wire_av[3] = "auto";
    if (!fail && ged_exec_dm(gedp, 4, wire_av) != BRLCAD_OK)
	fail = 1;

    const char *av[3] = {"draw", "--eager-leaf-expansion", "all.g"};
    if (!fail && ged_exec_draw(gedp, 3, av) != BRLCAD_OK) {
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
	const char *deferred_av[2] = {"draw", "all.g"};
	if (ged_exec_erase(gedp, 2, erase_av) != BRLCAD_OK ||
	    ged_exec_draw(gedp, 2, deferred_av) != BRLCAD_OK) {
	    bu_log("FAIL: owned-endpoint automatic deferred draw failed: %s\n",
		bu_vls_cstr(gedp->ged_result_str));
	    fail = 1;
	}
    }

    BObolProgressiveStatus root_status;
    if (!fail &&
	(controller->renderToImage(&deferred_root_image, 1, 0, NULL, NULL,
	    &root_status) != BRLCAD_OK || !deferred_root_image)) {
	bu_log("FAIL: automatic deferred draw did not render a root frame\n");
	fail = 1;
    }

    BObolProgressiveStatus full_status;
    if (!fail && root_status.hasMore) {
	const int drained =
	    draw_test_obol_progressive_drain(gedp, view_ctx, 500, 1);
	const int rendered =
	    controller->renderToImage(&deferred_full_image, 1, 0, NULL, NULL,
		&full_status);
	const double score = deferred_full_image ?
	    rgb_image_ssim(image, deferred_full_image, 384u, 384u) : -1.0;
	if (!drained || rendered != BRLCAD_OK || !deferred_full_image ||
	    full_status.hasMore || score < 0.98) {
	    bu_log("FAIL: deferred refinement did not reproduce a visually equivalent terminal image (SSIM %.9g)\n",
		score);
	    fail = 1;
	} else if (!fail) {
	    bu_log("PASS: interactive draw published %s then refined to a visually equivalent terminal image (SSIM %.9g)\n",
		root_status.changed ? "a progressive root" : "pending final data",
		score);
	}
    } else if (!fail) {
	/* A small job may complete before the first frame.  Prefer its final
	 * geometry instead of forcing a visible intermediate proxy frame. */
	const double root_score =
	    rgb_image_ssim(image, deferred_root_image, 384u, 384u);
	const int rendered =
	    controller->renderToImage(&deferred_full_image, 1, 0, NULL, NULL,
		&full_status);
	const double full_score = deferred_full_image ?
	    rgb_image_ssim(image, deferred_full_image, 384u, 384u) : -1.0;
	if (root_score < 0.98 || rendered != BRLCAD_OK ||
	    !deferred_full_image || full_status.hasMore || full_score < 0.98) {
	    bu_log("FAIL: completed deferred work did not use a visually equivalent terminal image (root SSIM %.9g, final SSIM %.9g)\n",
		root_score, full_score);
	    fail = 1;
	} else {
	    bu_log("PASS: interactive draw used visually equivalent completed geometry immediately (SSIM %.9g)\n",
		full_score);
	}
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

/* ========================================================================== */
/* Test 5: endpoint-native dm host lifecycle                                  */
/* ========================================================================== */
static int
test_endpoint_dm_lifecycle(const char *datadir)
{
    bu_log("\n--- Test 5: endpoint-native dm host lifecycle ---\n");

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", datadir);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmp("mged_view_state_t6.g", std::ios::binary);
    tmp << orig.rdbuf();
    orig.close(); tmp.close();
    bu_vls_free(&fname);

    struct ged *gedp = ged_open("db", "mged_view_state_t6.g", 1);
    if (!gedp) {
	bu_file_delete("mged_view_state_t6.g");
	return 1;
    }

    int fail = 0;
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    bv_dimensions_set(DRAW_TEST_BV(view_ctx), 0, 0);
    if (ged_view_context_obol_endpoint_get(view_ctx)) {
	bu_log("FAIL: lifecycle test unexpectedly started with an endpoint\n");
	fail = 1;
    }

    const char *open_av[7] = {
	"dm", "open", "--host", "headless", "--renderer", "sw", NULL
    };
    if (!fail && ged_exec_dm(gedp, 6, open_av) != BRLCAD_OK) {
	bu_log("FAIL: dm open could not create a headless software endpoint: %s\n",
		bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }

    bobol_display_endpoint_t *endpoint =
	ged_view_context_obol_endpoint_get(view_ctx);
    void *controller = endpoint ?
	bobol_display_endpoint_controller(endpoint) : NULL;
    if (!fail && (!endpoint || !controller ||
	    bobol_display_endpoint_render_engine_get(endpoint) !=
		BOBOL_RENDER_ENGINE_SW ||
	    !bobol_display_endpoint_host_factory_name(endpoint) ||
	    !BU_STR_EQUAL(bobol_display_endpoint_host_factory_name(endpoint),
		"headless"))) {
	bu_log("FAIL: dm open did not establish the requested endpoint state\n");
	fail = 1;
    }
    struct bv_display_property_value endpoint_width =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    struct bv_display_property_value endpoint_height =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    if (!fail && (bv_width_get(DRAW_TEST_BV(view_ctx)) != 512 ||
	bv_height_get(DRAW_TEST_BV(view_ctx)) != 512 ||
	bobol_display_endpoint_property_get(endpoint, "endpoint.width",
	    &endpoint_width) != BV_DISPLAY_PROPERTY_OK ||
	bobol_display_endpoint_property_get(endpoint, "endpoint.height",
	    &endpoint_height) != BV_DISPLAY_PROPERTY_OK ||
	endpoint_width.uint_value != 512 || endpoint_height.uint_value != 512)) {
	bu_log("FAIL: dm open did not seed a usable canonical viewport\n");
	fail = 1;
    }

    const char *host_av[3] = {"dm", "host", NULL};
    if (!fail && (ged_exec_dm(gedp, 2, host_av) != BRLCAD_OK ||
	    !BU_STR_EQUAL(bu_vls_cstr(gedp->ged_result_str), "headless"))) {
	bu_log("FAIL: dm host did not report headless: %s\n",
		bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }

    const char *diagnostics_av[3] = {"dm", "diagnostics", NULL};
    if (!fail && (ged_exec_dm(gedp, 2, diagnostics_av) != BRLCAD_OK ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str), "host=headless") ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str), "renderer=sw") ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str),
		"renderer.resolved=sw") ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str),
		"property.endpoint.renderer type=enum access=rw") ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str),
		"property.endpoint.renderer.resolved type=enum access=r"))) {
	bu_log("FAIL: dm diagnostics omitted endpoint state: %s\n",
		bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }

    const char *bg_set_av[9] = {
	"dm", "bg", "0", "0", "32", "64", "0", "0", NULL
    };
    if (!fail && ged_exec_dm(gedp, 8, bg_set_av) != BRLCAD_OK) {
	bu_log("FAIL: dm bg could not set an endpoint gradient: %s\n",
		bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    const char *bg_get_av[3] = {"dm", "bg", NULL};
    if (!fail && (ged_exec_dm(gedp, 2, bg_get_av) != BRLCAD_OK ||
	    !BU_STR_EQUAL(bu_vls_cstr(gedp->ged_result_str),
		"0/0/32->64/0/0\n"))) {
	bu_log("FAIL: dm bg endpoint gradient did not round trip: %s\n",
		bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }

    const char *policy_av[5] = {"dm", "set", "view.zclip", "1", NULL};
    if (!fail && (ged_exec_dm(gedp, 4, policy_av) != BRLCAD_OK ||
	    !bv_zclip_get(DRAW_TEST_BV(view_ctx)))) {
	bu_log("FAIL: dm set view.zclip did not update the owning bv view\n");
	fail = 1;
    }
    BObolViewController *view_controller =
	static_cast<BObolViewController *>(controller);
    policy_av[2] = "renderer.depth_test";
    policy_av[3] = "0";
    if (!fail && (ged_exec_dm(gedp, 4, policy_av) != BRLCAD_OK ||
	    view_controller->isDepthTestEnabled())) {
	bu_log("FAIL: dm set renderer.depth_test did not update Obol depth state\n");
	fail = 1;
    }
    policy_av[2] = "renderer.lighting";
    if (!fail && (ged_exec_dm(gedp, 4, policy_av) != BRLCAD_OK ||
	    view_controller->isLightingEnabled())) {
	bu_log("FAIL: dm set renderer.lighting did not update Obol lighting state\n");
	fail = 1;
    }
    policy_av[3] = "!";
    if (!fail && (ged_exec_dm(gedp, 4, policy_av) != BRLCAD_OK ||
	    !view_controller->isLightingEnabled())) {
	bu_log("FAIL: dm set renderer.lighting ! did not toggle typed lighting state\n");
	fail = 1;
    }
    policy_av[2] = "renderer.depth_cue";
    policy_av[3] = "1";
    if (!fail && (ged_exec_dm(gedp, 4, policy_av) != BRLCAD_OK ||
	    !view_controller->isDepthCueEnabled())) {
	bu_log("FAIL: dm set renderer.depth_cue did not update Obol fog state\n");
	fail = 1;
    }
    const char *retired_aliases[] = {
	"software_wire", "zclip", "zbuffer", "lighting", "depthcue"
    };
    for (size_t i = 0; !fail && i < sizeof(retired_aliases) /
	    sizeof(retired_aliases[0]); i++) {
	const char *retired_av[5] = {
	    "dm", "set", retired_aliases[i], "1", NULL
	};
	if (ged_exec_dm(gedp, 4, retired_av) != BRLCAD_ERROR) {
	    bu_log("FAIL: dm retained legacy property alias %s\n",
		retired_aliases[i]);
	    fail = 1;
	}
    }

    const char *renderer_none_av[4] = {"dm", "renderer", "none", NULL};
    if (!fail && (ged_exec_dm(gedp, 3, renderer_none_av) != BRLCAD_OK ||
	    view_controller->getRenderRoot())) {
	bu_log("FAIL: dm renderer none did not retain and suppress the graphical scene\n");
	fail = 1;
    }
    const char *renderer_diagnostic_av[4] = {
	"dm", "renderer", "diagnostic", NULL
    };
    if (!fail && (ged_exec_dm(gedp, 3, renderer_diagnostic_av) != BRLCAD_OK ||
	    view_controller->getRenderRoot() ||
	    ged_exec_dm(gedp, 2, diagnostics_av) != BRLCAD_OK ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str), "renderer=diagnostic") ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str),
		"property.endpoint.renderer.supported type=string access=r") ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str),
		"property.renderer.diagnostic.summary type=string access=r") ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str), "scene=available"))) {
	bu_log("FAIL: diagnostic renderer did not traverse and report typed state: %s\n",
		bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    const char *renderer_sw_av[4] = {"dm", "renderer", "sw", NULL};
    if (!fail && (ged_exec_dm(gedp, 3, renderer_sw_av) != BRLCAD_OK ||
	    !view_controller->getRenderRoot())) {
	bu_log("FAIL: leaving diagnostic renderer did not restore graphical presentation\n");
	fail = 1;
    }

    const char *close_av[3] = {"dm", "close", NULL};
    if (!fail && ged_exec_dm(gedp, 2, close_av) != BRLCAD_OK) {
	bu_log("FAIL: dm close rejected an open endpoint: %s\n",
		bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    if (!fail && (ged_view_context_obol_endpoint_get(view_ctx) != endpoint ||
	    bobol_display_endpoint_controller(endpoint) != controller ||
	    bobol_display_endpoint_host_factory_name(endpoint))) {
	bu_log("FAIL: dm close did not retain endpoint scene identity\n");
	fail = 1;
    }
    if (!fail && (ged_exec_dm(gedp, 2, host_av) != BRLCAD_OK ||
	    !BU_STR_EQUAL(bu_vls_cstr(gedp->ged_result_str), "unbound"))) {
	bu_log("FAIL: dm host did not report the closed endpoint as unbound\n");
	fail = 1;
    }

    const char *attach_bad_av[4] = {"dm", "attach", "obol", NULL};
    if (!fail && ged_exec_dm(gedp, 3, attach_bad_av) != BRLCAD_ERROR) {
	bu_log("FAIL: GED dm still exposes the compatibility attach command\n");
	fail = 1;
    }
    if (!fail && bobol_display_endpoint_host_factory_name(endpoint)) {
	bu_log("FAIL: rejected dm attach changed the endpoint host\n");
	fail = 1;
    }

    const char *list_av[3] = {"dm", "list", NULL};
    if (!fail && (ged_exec_dm(gedp, 2, list_av) != BRLCAD_OK ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str),
		bv_name_get(DRAW_TEST_BV(view_ctx))))) {
	bu_log("FAIL: bare dm list did not report retained endpoints\n");
	fail = 1;
    }
    open_av[3] = "no-such-host";
    if (!fail && ged_exec_dm(gedp, 6, open_av) != BRLCAD_ERROR) {
	bu_log("FAIL: dm open accepted an unknown host factory\n");
	fail = 1;
    }
    if (!fail && (ged_view_context_obol_endpoint_get(view_ctx) != endpoint ||
	    bobol_display_endpoint_controller(endpoint) != controller ||
	    bobol_display_endpoint_host_factory_name(endpoint))) {
	bu_log("FAIL: failed dm open changed retained endpoint identity\n");
	fail = 1;
    }

    open_av[3] = "headless";
    open_av[5] = "hw";
    if (!fail && ged_exec_dm(gedp, 6, open_av) != BRLCAD_ERROR) {
	bu_log("FAIL: dm open accepted hardware rendering on headless host\n");
	fail = 1;
    }
    open_av[5] = "sw";
    if (!fail && ged_exec_dm(gedp, 6, open_av) != BRLCAD_OK)
	fail = 1;
    if (!fail && ged_exec_dm(gedp, 2, close_av) != BRLCAD_OK)
	fail = 1;
    open_av[5] = "auto";
    if (!fail && (ged_exec_dm(gedp, 6, open_av) != BRLCAD_OK ||
	    bobol_display_endpoint_render_engine_get(endpoint) !=
		BOBOL_RENDER_ENGINE_AUTO ||
	    bobol_display_endpoint_render_engine_resolved_get(endpoint) !=
		BOBOL_RENDER_ENGINE_SW ||
	    ged_exec_dm(gedp, 2, diagnostics_av) != BRLCAD_OK ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str),
		"renderer=auto\n") ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str),
		"renderer.resolved=sw"))) {
	bu_log("FAIL: dm auto did not resolve headless presentation to software: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    if (!fail && ged_exec_dm(gedp, 2, close_av) != BRLCAD_OK)
	fail = 1;
    open_av[5] = "rt";
    if (!fail && (ged_exec_dm(gedp, 6, open_av) != BRLCAD_OK ||
	    ged_exec_dm(gedp, 2, diagnostics_av) != BRLCAD_OK ||
	    !strstr(bu_vls_cstr(gedp->ged_result_str), "renderer=rt"))) {
	bu_log("FAIL: dm open did not select retained rt on headless host: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
	fail = 1;
    }
    if (!fail && ged_exec_dm(gedp, 2, close_av) != BRLCAD_OK)
	fail = 1;

    if (!fail)
	bu_log("PASS: dm open/close/host/diagnostics use retained endpoints\n");
    ged_close(gedp);
    bu_file_delete("mged_view_state_t6.g");
    return fail;
}

/* ========================================================================== */
/* Test 6: session framebuffer capture provider belongs to one endpoint       */
/* ========================================================================== */
static int
test_framebuffer_capture_provider_rebind(const char *datadir)
{
    bu_log("\n--- Test 6: framebuffer capture-provider rebinding ---\n");

    struct bu_vls fname = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fname, "%s/moss.g", datadir);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmp("mged_view_state_t7.g", std::ios::binary);
    tmp << orig.rdbuf();
    orig.close(); tmp.close();
    bu_vls_free(&fname);

    struct ged *gedp = open_gedp("mged_view_state_t7.g", 320, 240);
    if (!gedp) {
	bu_file_delete("mged_view_state_t7.g");
	return 1;
    }

    int fail = 0;
    struct ged_view_context *first_view = ged_view_active_ctx(gedp);
    bobol_display_endpoint_t *first_endpoint =
	ged_view_context_obol_endpoint_get(first_view);
    if (!first_endpoint ||
	ged_view_framebuffer_backend_ensure(gedp,
	    first_view) != BRLCAD_OK) {
	bu_log("FAIL: could not bind the first framebuffer endpoint\n");
	fail = 1;
    }

    unsigned char *pixels = NULL;
    size_t size = 0;
    unsigned int width = 0;
    unsigned int height = 0;
    unsigned int components = 0;
    if (!fail && !bobol_display_endpoint_capture_plane(first_endpoint,
	    BOBOL_CAPTURE_FRAMEBUFFER, &pixels, &size, &width, &height,
	    &components)) {
	bu_log("FAIL: first endpoint did not expose the framebuffer provider\n");
	fail = 1;
    }
    if (pixels)
	bu_free(pixels, "first framebuffer capture");

    struct ged_view_set *view_set_ctx = ged_view_set_ctx(gedp);
    struct ged_view_context *second_view =
	ged_view_context_create_with_set(view_set_ctx);
    if (!second_view) {
	bu_log("FAIL: secondary framebuffer view creation failed\n");
	fail = 1;
    } else {
	bv_name_set(DRAW_TEST_BV(second_view), "framebuffer-secondary");
	bv_dimensions_set(DRAW_TEST_BV(second_view), 320, 240);
	bv_unit_conversion_set(DRAW_TEST_BV(second_view),
	    gedp->dbip->dbi_local2base, gedp->dbip->dbi_base2local);
	ged_view_set_context_add(view_set_ctx, second_view);
	ged_view_context_owned_add(gedp, second_view);
	const char *open_av[] = {"dm", "open", "-V",
	    "framebuffer-secondary", "--host", "headless", "--renderer",
	    "sw", NULL};
	const int endpoint_opened = ged_exec_dm(gedp, 8, open_av);
	const int bridge_bound = endpoint_opened == BRLCAD_OK ?
	    ged_view_framebuffer_backend_ensure(gedp,
		second_view) : BRLCAD_ERROR;
	if (endpoint_opened != BRLCAD_OK || bridge_bound != BRLCAD_OK) {
	    bu_log("FAIL: could not bind the second framebuffer endpoint "
		"(open=%d bridge=%d result=%s)\n", endpoint_opened,
		bridge_bound, bu_vls_cstr(gedp->ged_result_str));
	    fail = 1;
	}
    }

    bobol_display_endpoint_t *second_endpoint = second_view ?
	ged_view_context_obol_endpoint_get(second_view) : NULL;
    pixels = NULL;
    size = 0;
    width = 0;
    height = 0;
    components = 0;
    if (!fail && bobol_display_endpoint_capture_plane(first_endpoint,
	BOBOL_CAPTURE_FRAMEBUFFER, &pixels, &size, &width, &height,
	&components)) {
	bu_log("FAIL: previous endpoint retained the session framebuffer provider\n");
	if (pixels) {
	    bu_free(pixels, "stale framebuffer capture");
	    pixels = NULL;
	}
	fail = 1;
    }
    if (!fail && !bobol_display_endpoint_capture_plane(second_endpoint,
	BOBOL_CAPTURE_FRAMEBUFFER, &pixels, &size, &width, &height,
	&components)) {
	bu_log("FAIL: second endpoint did not receive the framebuffer provider\n");
	fail = 1;
    }
    if (pixels)
	bu_free(pixels, "second framebuffer capture");

    /* Endpoint replacement can also happen for one existing view.  Keep the
     * original endpoint alive for this assertion: the replacement must clear
     * its provider before the record stops owning it. */
    if (!fail && ged_view_framebuffer_backend_ensure(gedp,
	first_view) != BRLCAD_OK) {
	bu_log("FAIL: could not return the framebuffer bridge to the first view\n");
	fail = 1;
    }
    if (!fail && !ged_view_context_obol_endpoint_set(first_view,
	first_endpoint, 0)) {
	bu_log("FAIL: could not preserve the first endpoint for same-view handoff\n");
	fail = 1;
    }
    bobol_display_endpoint_t *replacement_endpoint = NULL;
    if (!fail) {
	replacement_endpoint = bobol_display_endpoint_create(NULL, 0);
	if (!replacement_endpoint || !ged_view_context_obol_endpoint_set(
		first_view, replacement_endpoint, 1) ||
	    ged_view_framebuffer_backend_ensure(gedp,
		first_view) != BRLCAD_OK) {
	    bu_log("FAIL: could not rebind framebuffer after same-view endpoint replacement\n");
	    if (replacement_endpoint &&
		ged_view_context_obol_endpoint_get(first_view) !=
		replacement_endpoint)
		bobol_display_endpoint_destroy(replacement_endpoint);
	    replacement_endpoint = NULL;
	    fail = 1;
	}
    }
    pixels = NULL;
    if (!fail && bobol_display_endpoint_capture_plane(first_endpoint,
	BOBOL_CAPTURE_FRAMEBUFFER, &pixels, &size, &width, &height,
	&components)) {
	bu_log("FAIL: same-view endpoint replacement retained the old capture provider\n");
	if (pixels)
	    bu_free(pixels, "same-view stale framebuffer capture");
	fail = 1;
    }
    pixels = NULL;
    if (!fail && (!replacement_endpoint ||
	!bobol_display_endpoint_capture_plane(replacement_endpoint,
	    BOBOL_CAPTURE_FRAMEBUFFER, &pixels, &size, &width, &height,
	    &components))) {
	bu_log("FAIL: same-view endpoint replacement did not bind its capture provider\n");
	fail = 1;
    }
    if (pixels)
	bu_free(pixels, "same-view framebuffer capture");

    ged_view_framebuffer_release(gedp);
    pixels = NULL;
	if (!fail && replacement_endpoint &&
	bobol_display_endpoint_capture_plane(replacement_endpoint,
	BOBOL_CAPTURE_FRAMEBUFFER, &pixels, &size, &width, &height,
	&components)) {
	bu_log("FAIL: framebuffer release left a capture provider behind\n");
	if (pixels)
	    bu_free(pixels, "released framebuffer capture");
	fail = 1;
    }

    if (!fail)
	bu_log("PASS: framebuffer capture provider moves between endpoints\n");
    ged_close(gedp);
    bobol_display_endpoint_destroy(first_endpoint);
    bu_file_delete("mged_view_state_t7.g");
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
    failures += test_owned_render_endpoint(datadir);
    failures += test_endpoint_dm_lifecycle(datadir);
    failures += test_framebuffer_capture_provider_rebind(datadir);

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
