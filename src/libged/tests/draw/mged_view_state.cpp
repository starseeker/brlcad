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
#include "rt/view.h"
#define DM_WITH_RT
#include <dm.h>
#include "dm/obol.h"
#include <ged.h>
#include "ged/draw.h"
#include "ged/draw_obol.h"
#include "brlobol/scene_controller.h"
#include "brlobol/view_controller.h"

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
    struct dm *dmp = (struct dm *)rt_view_context_display_manager_get(view_ctx);
    return dmp ? (BRLObolViewController *)dm_obol_controller(dmp) : NULL;
}

static int
configure_obol_view(struct ged *gedp, void *view_ctx, int width, int height)
{
    struct dm *dmp = (struct dm *)rt_view_context_display_manager_get(view_ctx);
    if (!dmp)
	return 0;

    dm_set_width(dmp, width);
    dm_set_height(dmp, height);
    dm_configure_win(dmp, 0);
    dm_set_zbuffer(dmp, 1);
    fastf_t wb[6] = {-1, 1, -1, 1, -100, 100};
    dm_set_win_bounds(dmp, wb);
    dm_set_vp(dmp, rt_view_context_scale_storage_get(view_ctx));
    rt_view_context_display_manager_set(view_ctx, dmp);
    rt_view_context_dimensions_set(view_ctx, dm_get_width(dmp),
	    dm_get_height(dmp));
    rt_view_context_unit_conversion_set(view_ctx,
	    gedp->dbip->dbi_local2base,
	    gedp->dbip->dbi_base2local);
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
    struct dm *dmp  = (struct dm *)rt_view_context_display_manager_get(v);
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
    dm_set_vp(dmp, rt_view_context_scale_storage_get(v));
    rt_view_context_display_manager_set(v, dmp);
    rt_view_context_dimensions_set(v, dm_get_width(dmp), dm_get_height(dmp));
    rt_view_context_unit_conversion_set(v,
	gedp->dbip->dbi_local2base,
	gedp->dbip->dbi_base2local);

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
/* Test 2: Edit matrix state is accepted by the neutral view context           */
/* ========================================================================== */
static int
test_edit_matrix(const char *datadir)
{
    bu_log("\n--- Test 2: edit matrix state set/clear ---\n");

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
    ged_draw_set_highlight_state(gedp, 1);

    rt_view_context_edit_matrix_clear(v);

    /* Apply a clear translation through the neutral edit-matrix state. */
    mat_t edit_mat;
    MAT_IDN(edit_mat);
    edit_mat[3] = 10.0;

    int edit_matrix_ready = rt_view_context_edit_matrix_set(v, edit_mat);
    int fail = 0;
    if (!edit_matrix_ready) {
	bu_log("FAIL: edit-matrix set failed on the active view context\n");
	fail = 1;
    } else {
	bu_log("PASS: edit-matrix set accepted by the active view context\n");
    }

    if (!rt_view_context_edit_matrix_clear(v)) {
	bu_log("FAIL: edit-matrix clear failed on the active view context\n");
	fail = 1;
    } else {
	bu_log("PASS: edit-matrix clear accepted by the active view context\n");
    }

    ged_draw_set_highlight_state(gedp, 0);
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
    rt_view_context_name_set(v0, "V0");

    void *view_set_ctx = ged_view_set_ctx(gedp);
    void *v1 = rt_view_context_create_with_set(view_set_ctx);
    if (!v1) {
	bu_log("FAIL: secondary view creation failed\n");
	fail = 1;
    } else {
	rt_view_context_name_set(v1, "V1");
	rt_view_context_unit_conversion_set(v1,
		gedp->dbip->dbi_local2base,
		gedp->dbip->dbi_base2local);
	rt_view_set_context_add(view_set_ctx, v1);
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

    ged_close(gedp);
    bu_file_delete("mged_view_state_t3.g");
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
    failures += test_edit_matrix(datadir);
    failures += test_multi_obol_dm_attachment(datadir);

    if (failures == 0) {
	bu_log("\nAll MGED view-state tests PASSED (%d/3)\n", 3);
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
