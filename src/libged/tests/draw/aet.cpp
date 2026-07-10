/*                         A E T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2018-2026 United States Government as represented by
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
/** @file aet.cpp
 *
 * Attempts to test the stability of views across various operations
 * such as resizing.
 *
 */

#include "common.h"

#include <stdio.h>
#include <fstream>
#include <bu.h>
#include <icv.h>
#define DM_WITH_RT
#include <dm.h>
#include <ged.h>
#include <ged/draw.h>
#include <rt/view.h>
#include "view_test_util.h"

extern "C" int unpack_apng(const char *src_dir, const char *apng_name, const char *out_dir, const char *prefix);
void
dm_refresh(struct ged *gedp, int vnum)
{
    struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
    void *v = views ? BU_PTBL_GET(views, vnum) : NULL;
    if (!v)
	return;
    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
    txn.view = v;
    ged_draw_apply_transaction(gedp, &txn, NULL);

    struct dm *dmp = (struct dm *)ged_view_context_display_manager_get(v);
    /* Ensure rendering goes to this view's DM context, not the last-active one.
     * With multiple DMs each view has its own OSMesa context; without making the
     * correct context current here dm_set_bg and dm_draw_objs will operate on
     * whichever context was last activated, leaving this buffer empty. */
    dm_make_current(dmp);
    unsigned char *dm_bg1;
    unsigned char *dm_bg2;
    dm_get_bg(&dm_bg1, &dm_bg2, dmp);
    dm_set_bg(dmp, dm_bg1[0], dm_bg1[1], dm_bg1[2], dm_bg2[0], dm_bg2[1], dm_bg2[2]);
    dm_set_native_repaint_pending(dmp, 0);
    dm_draw_objs(v);
    dm_draw_end(dmp);
}

int
img_cmp(int vnum, int id, struct ged *gedp, const char *cdir, int soft_fail)
{
    icv_image_t *ctrl, *timg;
    struct bu_vls tname = BU_VLS_INIT_ZERO;
    struct bu_vls cname = BU_VLS_INIT_ZERO;
    if (id <= 0) {
	bu_vls_sprintf(&tname, "aet_clear.png");
	bu_vls_sprintf(&cname, "%s/empty.png", cdir);
    } else {
	bu_vls_sprintf(&tname, "aet_%02d_%03d.png", vnum, id);
	bu_vls_sprintf(&cname, "%s/aet_%02d_%03d_ctrl.png", cdir, vnum, id);
    }

    dm_refresh(gedp, vnum);

    struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
    void *v = views ? BU_PTBL_GET(views, vnum) : NULL;
    if (!v)
	bu_exit(EXIT_FAILURE, "Invalid view specifier: %d\n", vnum);
    struct dm *dmp = (struct dm *)ged_view_context_display_manager_get(v);

    const char *s_av[4] = {NULL};
    s_av[0] = "screengrab";
    s_av[1] = "-D";
    s_av[2] = bu_vls_cstr(dm_get_pathname(dmp));
    s_av[3] = bu_vls_cstr(&tname);
    if (ged_exec_screengrab(gedp, 4, s_av) & BRLCAD_ERROR) {
	bu_log("Failed to grab screen for DM %s\n", bu_vls_cstr(dm_get_pathname(dmp)));
	bu_vls_free(&tname);
	return BRLCAD_ERROR;
    }

    timg = icv_read(bu_vls_cstr(&tname), BU_MIME_IMAGE_PNG, 0, 0);
    if (!timg) {
	if (soft_fail) {
	    bu_log("Failed to read %s\n", bu_vls_cstr(&tname));
	    bu_vls_free(&tname);
	    return BRLCAD_ERROR;
	}
	bu_exit(EXIT_FAILURE, "failed to read %s\n", bu_vls_cstr(&tname));
    }
    ctrl = icv_read(bu_vls_cstr(&cname), BU_MIME_IMAGE_PNG, 0, 0);
    if (!ctrl) {
	if (soft_fail) {
	    bu_log("Failed to read %s\n", bu_vls_cstr(&cname));
	    bu_vls_free(&tname);
	    bu_vls_free(&cname);
	    return BRLCAD_ERROR;
	}
	bu_exit(EXIT_FAILURE, "failed to read %s\n", bu_vls_cstr(&cname));
    }
    bu_vls_free(&cname);
    int matching_cnt = 0;
    int off_by_1_cnt = 0;
    int off_by_many_cnt = 0;
    int iret = icv_diff(&matching_cnt, &off_by_1_cnt, &off_by_many_cnt, ctrl,timg);
    if (iret) {
	if (soft_fail) {
	    bu_log("%d wireframe diff failed.  %d matching, %d off by 1, %d off by many\n", id, matching_cnt, off_by_1_cnt, off_by_many_cnt);
	    icv_destroy(ctrl);
	    icv_destroy(timg);
	    return BRLCAD_ERROR;
	}
	bu_exit(EXIT_FAILURE, "%d wireframe diff failed.  %d matching, %d off by 1, %d off by many\n", id, matching_cnt, off_by_1_cnt, off_by_many_cnt);
    }

    icv_destroy(ctrl);
    icv_destroy(timg);

    // Image comparison done and successful
    bu_file_delete(bu_vls_cstr(&tname));
    bu_vls_free(&tname);

    return BRLCAD_OK;
}

int
main(int ac, char *av[]) {
    struct ged *gedp;
    struct bu_vls fname = BU_VLS_INIT_ZERO;
    int need_help = 0;
    int run_unstable_tests = 0;
    int soft_fail = 0;
    int ret = BRLCAD_OK;

    bu_setprogname(av[0]);

    struct bu_opt_desc d[4];
    BU_OPT(d[0], "h", "help",            "", NULL,          &need_help, "Print help and exit");
    BU_OPT(d[1], "U", "enable-unstable", "", NULL, &run_unstable_tests, "Test drawing routines known to differ between build configs/platforms.");
    BU_OPT(d[2], "c", "continue",        "", NULL,          &soft_fail, "Continue testing if a failure is encountered.");
    BU_OPT_NULL(d[3]);

    /* Done with program name */
    (void)bu_opt_parse(NULL, ac, (const char **)av, d);

    if (!bu_file_directory(av[1])) {
	printf("ERROR: [%s] is not a directory.  Expecting control image directory\n", av[1]);
	return 2;
    }

    /* Enable all the experimental logic */
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    /* Use a local working-directory cache so we do not pollute the user's
     * real BRL-CAD cache and so the test is fully self-contained. */
    char lcache[MAXPATHLEN] = {0};
    bu_dir(lcache, MAXPATHLEN, BU_DIR_CURR, "ged_aet_test_cache", NULL);
    bu_mkdir(lcache);
    bu_setenv("BU_DIR_CACHE", lcache, 1);

    unpack_apng(av[1], "aet_00.apng", lcache, "aet_00_");
    unpack_apng(av[1], "aet_01.apng", lcache, "aet_01_");
    unpack_apng(av[1], "aet_02.apng", lcache, "aet_02_");
    unpack_apng(av[1], "aet_03.apng", lcache, "aet_03_");

    if (!bu_file_exists(av[1], NULL)) {
	printf("ERROR: [%s] does not exist, expecting .g file\n", av[1]);
	return 2;
    }

    /* Use a local working-directory cache so we do not pollute the user's
     * real BRL-CAD cache and so the test is fully self-contained. */
    bu_dir(lcache, MAXPATHLEN, BU_DIR_CURR, "ged_aet_test_cache", NULL);
    bu_mkdir(lcache);
    bu_setenv("BU_DIR_CACHE", lcache, 1);

    /* make a temporary copy of moss */
    bu_vls_sprintf(&fname, "%s/moss.g", av[1]);
    std::ifstream orig(bu_vls_cstr(&fname), std::ios::binary);
    std::ofstream tmpg("moss_aet_tmp.g", std::ios::binary);
    tmpg << orig.rdbuf();
    orig.close();
    tmpg.close();

    /* Open the temp file */
    const char *s_av[15] = {NULL};
    gedp = ged_open("db", "moss_aet_tmp.g", 1);

    bu_setenv("DM_SWRAST", "1", 1);

    // We don't want the default GED views for this test
    void *view_set_ctx = ged_view_set_ctx(gedp);
    ged_view_set_context_remove(view_set_ctx, NULL);

    // Set up the views.  Unlike the other drawing tests, we are explicitly
    // out to test the behavior of multiple views and dms, so we need to
    // set up multiples.  We'll start out with four non-independent views,
    // to mimic the most common multi-dm/view display - a Quad view widget.
    // Each view will get its own attached swrast DM.
    void *views[4];
    for (size_t i = 0; i < 4; i++) {
	char view_name[16];
	snprintf(view_name, sizeof(view_name), "V%zd", i);
	views[i] = ged_view_context_create_with_set(view_set_ctx);
	if (!i)
	    ged_view_active_ctx_set(gedp, views[i]);
	bv_name_set(DRAW_TEST_BV(views[i]), view_name);
	ged_view_set_context_add(view_set_ctx, views[i]);
	ged_view_context_owned_add(gedp, views[i]);

	/* To generate images that will allow us to check if the drawing
	 * is proceeding as expected, we use the swrast off-screen dm. */
	struct bu_vls dm_name = BU_VLS_INIT_ZERO;
	s_av[0] = "dm";
	s_av[1] = "attach";
	s_av[2] = "-V";
	s_av[3] = view_name;
	s_av[4] = "swrast";
	bu_vls_sprintf(&dm_name, "SW%zd", i);
	s_av[5] = bu_vls_cstr(&dm_name);
	s_av[6] = NULL;
	ged_exec_dm(gedp, 6, s_av);
	bu_vls_free(&dm_name);

	struct dm *dmp = (struct dm *)ged_view_context_display_manager_get(views[i]);
	dm_set_width(dmp, 512);
	dm_set_height(dmp, 512);

	dm_configure_win(dmp, 0);
	dm_set_zbuffer(dmp, 1);

	// See QtSW.cpp...
	fastf_t windowbounds[6] = { -1, 1, -1, 1, -100, 100 };
	dm_set_win_bounds(dmp, windowbounds);

	dm_set_vp(dmp, bv_scale_storage_get(DRAW_TEST_BV(views[i])));
	ged_view_context_display_manager_set(views[i], dmp);
	bv_dimensions_set(DRAW_TEST_BV(views[i]), dm_get_width(dmp), dm_get_height(dmp));
	bv_unit_conversion_set(DRAW_TEST_BV(views[i]), gedp->dbip->dbi_local2base, gedp->dbip->dbi_base2local);
    }

    /* Set distinct view az/el for each of the four quad views.  For
     * this test we are deliberately testing view settings that have
     * the potential to be challenging in "gimbal lock" positions in
     * multiples of 90 degrees and using non-zero twist components. */
    vect_t aet = VINIT_ZERO;
    VSET(aet, 0, 0, 90);
    bv_aet_set(DRAW_TEST_BV(views[0]), aet);
    ged_view_context_update(views[0]);

    VSET(aet, 90, 90, 180);
    bv_aet_set(DRAW_TEST_BV(views[1]), aet);
    ged_view_context_update(views[1]);

    VSET(aet, -90, 270, -90);
    bv_aet_set(DRAW_TEST_BV(views[2]), aet);
    ged_view_context_update(views[2]);

    VSET(aet, 270, -180, 90);
    bv_aet_set(DRAW_TEST_BV(views[3]), aet);
    ged_view_context_update(views[3]);


    /************************************************************************/
    /* Draw a wireframe.  Because we're doing tests of view and dm
     * manipulation, this is the only draw command needed for this particular
     * setup - we just need some non-empty geometry displayed. */
    bu_log("Basic drawing test...\n");
    s_av[0] = "draw";
    s_av[1] = "-m0";
    s_av[2] = "all.g";
    s_av[3] = NULL;
    ged_exec_draw(gedp, 4, s_av);

    // Sanity
    ret += img_cmp(0, 1, gedp, lcache, soft_fail);
    ret += img_cmp(1, 1, gedp, lcache, soft_fail);
    ret += img_cmp(2, 1, gedp, lcache, soft_fail);
    ret += img_cmp(3, 1, gedp, lcache, soft_fail);

    // Resize dm to larger dimensions
    bu_log("Resize to 600x600...\n");
    int len = 600;
    for (size_t i = 0; i < 4; i++) {
	struct dm *dmp = (struct dm *)ged_view_context_display_manager_get(views[i]);
	dm_set_width(dmp, len);
	dm_set_height(dmp, len);
	bv_dimensions_set(DRAW_TEST_BV(views[i]), len, len);
	dm_configure_win(dmp, 0);
	// NOTE:  deliberately not resetting aet here - we want to see if it is
	// stable without adjustment.
	ged_view_context_update(views[i]);
    }
    ret += img_cmp(0, 2, gedp, lcache, soft_fail);
    ret += img_cmp(1, 2, gedp, lcache, soft_fail);
    ret += img_cmp(2, 2, gedp, lcache, soft_fail);
    ret += img_cmp(3, 2, gedp, lcache, soft_fail);

    // Shrink back to default dimensions
    bu_log("Shrink to 512x512...\n");
    len = 512;
    for (size_t i = 0; i < 4; i++) {
	struct dm *dmp = (struct dm *)ged_view_context_display_manager_get(views[i]);
	dm_set_width(dmp, len);
	dm_set_height(dmp, len);
	bv_dimensions_set(DRAW_TEST_BV(views[i]), len, len);
	dm_configure_win(dmp, 0);
	// NOTE:  deliberately not resetting aet here - we want to see if it is
	// stable without adjustment.
	ged_view_context_update(views[i]);
    }
    ret += img_cmp(0, 1, gedp, lcache, soft_fail);
    ret += img_cmp(1, 1, gedp, lcache, soft_fail);
    ret += img_cmp(2, 1, gedp, lcache, soft_fail);
    ret += img_cmp(3, 1, gedp, lcache, soft_fail);

    // Cycle through a bunch of resizes
    bu_log("Cycle through multiple resizes...\n");
    for (int i = 513; i < 600; i++) {
	for (size_t j = 0; j < 4; j++) {
	    struct dm *dmp = (struct dm *)ged_view_context_display_manager_get(views[j]);
	    dm_set_width(dmp, i);
	    dm_set_height(dmp, i);
	    bv_dimensions_set(DRAW_TEST_BV(views[j]), i, i);
	    dm_configure_win(dmp, 0);
	    // NOTE:  deliberately not resetting aet here - we want to see if it is
	    // stable without adjustment.
	    ged_view_context_update(views[j]);
	}
    }
    len = 512;
    for (size_t i = 0; i < 4; i++) {
	struct dm *dmp = (struct dm *)ged_view_context_display_manager_get(views[i]);
	dm_set_width(dmp, len);
	dm_set_height(dmp, len);
	bv_dimensions_set(DRAW_TEST_BV(views[i]), len, len);
	dm_configure_win(dmp, 0);
	// NOTE:  deliberately not resetting aet here - we want to see if it is
	// stable without adjustment.
	ged_view_context_update(views[i]);
    }
    ret += img_cmp(0, 1, gedp, lcache, soft_fail);
    ret += img_cmp(1, 1, gedp, lcache, soft_fail);
    ret += img_cmp(2, 1, gedp, lcache, soft_fail);
    ret += img_cmp(3, 1, gedp, lcache, soft_fail);

    ged_close(gedp);

    return ret;
}


// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
