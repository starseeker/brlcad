/*            V I E W _ C O M M A N D . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */

#include "common.h"

#include <cmath>
#include <cstdio>
#include <string>

#include <bu.h>
#include <ged.h>
#include <ged/draw.h>
#include <rt/view_legacy_bsg.h>

#define ASSERT(cond) do { \
    nchecks++; \
    if (!(cond)) { \
	bu_log("FAIL [%s:%d] %s\n", __FILE__, __LINE__, #cond); \
	nfails++; \
    } \
} while (0)

static int nchecks = 0;
static int nfails = 0;

static int
run_view(struct ged *gedp, int argc, const char **argv)
{
    int ret = ged_exec_view(gedp, argc, argv);
    return ret;
}

static std::string
result_str(struct ged *gedp)
{
    const char *r = bu_vls_cstr(gedp->ged_result_str);
    return (r) ? std::string(r) : std::string();
}

static void
assert_view_ok(struct ged *gedp, int argc, const char **argv, int line)
{
    nchecks++;
    int ret = run_view(gedp, argc, argv);
    if (ret != BRLCAD_OK) {
	bu_log("FAIL [%s:%d] view command returned %d: ", __FILE__, line, ret);
	for (int i = 0; i < argc; i++)
	    bu_log("%s%s", i ? " " : "", argv[i]);
	bu_log("\n  result: %s\n", bu_vls_cstr(gedp->ged_result_str));
	nfails++;
    }
}

#define ASSERT_VIEW_OK(gedp, argc, argv) assert_view_ok((gedp), (argc), (argv), __LINE__)

static void
refresh_scene_records(struct ged *gedp, void *v)
{
    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
    txn.view = v;
    ASSERT(ged_draw_apply_transaction(gedp, &txn, NULL) >= 0);
}

static void
test_command_report_record_consistency(struct ged *gedp, void *v)
{
    const char *draw_av[] = {"draw", "all.g", NULL};
    ASSERT(ged_exec_draw(gedp, 2, draw_av) == BRLCAD_OK);
    refresh_scene_records(gedp, v);

    const char *who_av[] = {"who", NULL};
    ASSERT(ged_exec_who(gedp, 1, who_av) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("all.g") != std::string::npos);

    const char *who_solids_av[] = {"who", "solids", "1", NULL};
    ASSERT(ged_exec_who(gedp, 3, who_solids_av) == BRLCAD_OK);
    std::string solids = result_str(gedp);
    ASSERT(solids.find("all.g") != std::string::npos);
    ASSERT(solids.find("cent=") != std::string::npos);

    struct rt_view_render_export_consistency_bsg consistency =
	RT_VIEW_RENDER_EXPORT_CONSISTENCY_BSG_INIT;
    ASSERT(rt_view_context_render_export_consistency_bsg(v, "all.g",
	    &consistency));
    ASSERT(consistency.export_record_found);
    ASSERT(consistency.render_item_found);
    ASSERT(consistency.export_render_consistent);
    ASSERT(consistency.backend_node_found);
    ASSERT(consistency.export_backend_consistent);

    struct rt_view_pick_result_bsg *pick =
	rt_view_context_pick_semantic_path_bsg(v, "all.g");
    struct bu_vls pick_path = BU_VLS_INIT_ZERO;
    ASSERT(pick != NULL);
    ASSERT(rt_view_pick_result_count_bsg(pick) > 0);
    if (pick && rt_view_pick_result_count_bsg(pick) > 0) {
	ASSERT(rt_view_pick_result_path_bsg(pick, 0, &pick_path));
	ASSERT(BU_STR_EQUAL(bu_vls_cstr(&pick_path), "all.g"));
    }

    point_t sample = VINIT_ZERO;
    struct rt_view_snap_result_bsg *snap = rt_view_snap_result_create_bsg();
    ASSERT(snap != NULL);
    int snap_count = snap ? rt_view_context_snap_candidates_result_bsg(v,
	    sample, 1.0, RT_VIEW_SNAP_KIND_ENDPOINT_BSG, snap) : 0;
    ASSERT(snap_count >= 0);
    if (snap_count > 0 && rt_view_snap_result_count_bsg(snap) > 0) {
	struct bu_vls snap_path = BU_VLS_INIT_ZERO;
	ASSERT(rt_view_snap_result_source_path_bsg(snap, 0, &snap_path));
	ASSERT(bu_vls_strlen(&snap_path) > 0);
	bu_vls_free(&snap_path);
    }

    point_t a = VINIT_ZERO;
    point_t b = {1.0, 0.0, 0.0};
    struct rt_view_measure_result measure = RT_VIEW_MEASURE_RESULT_INIT;
    ASSERT(rt_view_context_measure_candidates_bsg(v, a, b, &measure) == 1);
    ASSERT(measure.valid);
    ASSERT(fabs(measure.distance - 1.0) < 1.0e-9);

    rt_view_snap_result_free_bsg(snap);
    bu_vls_free(&pick_path);
    rt_view_pick_result_free_bsg(pick);
}

int
main(int argc, const char **argv)
{
    bu_setprogname(argv[0]);

    if (argc != 2)
	bu_exit(EXIT_FAILURE, "Usage: ged_test_view_command <directory-containing-moss.g>\n");

    struct bu_vls gpath = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&gpath, "%s/moss.g", argv[1]);
    struct ged *gedp = ged_open("db", bu_vls_cstr(&gpath), 1);
    ASSERT(gedp != NULL);
    if (!gedp)
	return EXIT_FAILURE;

    void *view_set_ctx = ged_view_set_ctx(gedp);
    ASSERT(rt_view_set_context_remove_bsg(view_set_ctx, NULL));
    void *views[2] = {NULL, NULL};
    for (int i = 0; i < 2; i++) {
	char view_name[16];
	snprintf(view_name, sizeof(view_name), "V%d", i);
	views[i] = rt_view_context_create_bsg();
	ASSERT(views[i] != NULL);
	ASSERT(rt_view_context_name_set_bsg(views[i], view_name));
	ASSERT(rt_view_context_dimensions_set_bsg(views[i], 640, 480));
	ASSERT(rt_view_set_context_add_bsg(view_set_ctx, views[i]));
	ged_view_context_owned_add(gedp, views[i]);
	if (!i)
	    ged_view_active_ctx_set(gedp, views[i]);
    }

    test_command_report_record_consistency(gedp, views[0]);

    const char *c0[] = {"view", "obj", "create", "u_line", "line", "create", "0", "0", "0", NULL};
    ASSERT(run_view(gedp, 9, c0) == BRLCAD_OK);

    const char *c1[] = {"view", "obj", "info", "u_line", "type", NULL};
    ASSERT(run_view(gedp, 5, c1) == BRLCAD_OK);

    const char *c2[] = {"view", "obj", "set", "u_line", "draw", "0", NULL};
    ASSERT(run_view(gedp, 6, c2) == BRLCAD_OK);
    const char *c3[] = {"view", "obj", "info", "u_line", "draw", NULL};
    ASSERT(run_view(gedp, 5, c3) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("DOWN") != std::string::npos);
    const char *c3a[] = {"view", "obj", "set", "u_line", "draw", "1", NULL};
    ASSERT(run_view(gedp, 6, c3a) == BRLCAD_OK);
    const char *c3b[] = {"view", "obj", "info", "u_line", "draw", NULL};
    ASSERT(run_view(gedp, 5, c3b) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("UP") != std::string::npos);

    const char *c3c[] = {"view", "obj", "set", "u_line", "color", "10/20/30", NULL};
    ASSERT(run_view(gedp, 6, c3c) == BRLCAD_OK);
    const char *c3d[] = {"view", "obj", "info", "u_line", "color", NULL};
    ASSERT(run_view(gedp, 5, c3d) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("10/20/30") != std::string::npos);

    const char *c4[] = {"view", "obj", "list", "u_*", NULL};
    ASSERT(run_view(gedp, 4, c4) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("u_line") != std::string::npos);

    const char *c5[] = {"view", "obj", "set", "u_line", "arrow", "1", NULL};
    ASSERT(run_view(gedp, 6, c5) == BRLCAD_OK);
    const char *c5a[] = {"view", "obj", "set", "u_line", "update", "12", "34", NULL};
    ASSERT(run_view(gedp, 7, c5a) == BRLCAD_OK);
    const char *c5b[] = {"view", "vZ", "-N", "u_line", NULL};
    ASSERT(run_view(gedp, 4, c5b) == BRLCAD_OK);
    ASSERT(!result_str(gedp).empty());

    const char *a0[] = {"view", "obj", "create", "u_axes", "axes", "create", "1", "2", "3", NULL};
    ASSERT_VIEW_OK(gedp, 9, a0);
    const char *a1[] = {"view", "obj", "create", "u_axes", "axes", "pos", NULL};
    ASSERT_VIEW_OK(gedp, 6, a1);
    ASSERT(result_str(gedp).find("1.000000 2.000000 3.000000") != std::string::npos);
    const char *a2[] = {"view", "obj", "create", "u_axes", "axes", "pos", "4", "5", "6", NULL};
    ASSERT_VIEW_OK(gedp, 9, a2);
    const char *a3[] = {"view", "obj", "create", "u_axes", "axes", "size", "12.5", NULL};
    ASSERT_VIEW_OK(gedp, 7, a3);
    const char *a4[] = {"view", "obj", "create", "u_axes", "axes", "size", NULL};
    ASSERT_VIEW_OK(gedp, 6, a4);
    ASSERT(result_str(gedp).find("12.500000") != std::string::npos);
    const char *a5[] = {"view", "obj", "create", "u_axes", "axes", "line_width", "3", NULL};
    ASSERT_VIEW_OK(gedp, 7, a5);
    const char *a6[] = {"view", "obj", "create", "u_axes", "axes", "line_width", NULL};
    ASSERT_VIEW_OK(gedp, 6, a6);
    ASSERT(result_str(gedp).find("3") != std::string::npos);
    const char *a7[] = {"view", "obj", "create", "u_axes", "axes", "axes_color", "10", "20", "30", NULL};
    ASSERT_VIEW_OK(gedp, 9, a7);
    const char *a8[] = {"view", "obj", "create", "u_axes", "axes", "axes_color", NULL};
    ASSERT_VIEW_OK(gedp, 6, a8);
    ASSERT(result_str(gedp).find("10 20 30") != std::string::npos);

    const char *p0[] = {"view", "obj", "create", "u_poly", "polygon", "create", "10", "10", NULL};
    ASSERT_VIEW_OK(gedp, 8, p0);
    const char *p1[] = {"view", "obj", "create", "u_poly", "polygon", "append", "30", "10", NULL};
    ASSERT_VIEW_OK(gedp, 8, p1);
    const char *p2[] = {"view", "obj", "create", "u_poly", "polygon", "append", "30", "30", NULL};
    ASSERT_VIEW_OK(gedp, 8, p2);
    const char *p3[] = {"view", "obj", "create", "u_poly", "polygon", "append", "10", "30", NULL};
    ASSERT_VIEW_OK(gedp, 8, p3);
    const char *p4[] = {"view", "obj", "create", "u_poly", "polygon", "close", NULL};
    ASSERT_VIEW_OK(gedp, 6, p4);
    const char *p5[] = {"view", "obj", "create", "u_poly", "polygon", "fill", "1", "0", "0.1", NULL};
    ASSERT_VIEW_OK(gedp, 9, p5);
    const char *p6[] = {"view", "obj", "create", "u_poly", "polygon", "fill_color", "10/20/30", NULL};
    ASSERT_VIEW_OK(gedp, 7, p6);
    const char *p7[] = {"view", "obj", "create", "u_poly", "polygon", "area", NULL};
    ASSERT_VIEW_OK(gedp, 6, p7);
    ASSERT(!result_str(gedp).empty());
    rt_view_polygon_ref poly_ref =
	rt_view_context_polygon_find_bsg(views[0], "u_poly");
    ASSERT(!rt_view_polygon_ref_is_null_bsg(poly_ref));
    struct rt_view_polygon_record poly_rec = {};
    ASSERT(rt_view_polygon_record_get_bsg(poly_ref, &poly_rec));
    ASSERT(poly_rec.type == RT_VIEW_POLYGON_GENERAL);
    ASSERT(poly_rec.contour_count == 1);
    ASSERT(poly_rec.point_count == 4);
    ASSERT(poly_rec.first_contour_open == 0);
    ASSERT(poly_rec.fill_flag == 1);

    const char *c6[] = {"view", "-V", "V0", "obj", "-L", "create", "l_line", "line", "create", "0", "0", "0", NULL};
    ASSERT(run_view(gedp, 12, c6) == BRLCAD_OK);
    const char *c7[] = {"view", "-V", "V0", "obj", "list", NULL};
    ASSERT(run_view(gedp, 5, c7) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("l_line") != std::string::npos);
    const char *c8[] = {"view", "-V", "V1", "obj", "list", NULL};
    ASSERT(run_view(gedp, 5, c8) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("l_line") == std::string::npos);

    const char *c11[] = {"view", "obj", "-g", "all.g", "create", "g2", NULL};
    ASSERT(run_view(gedp, 6, c11) == BRLCAD_OK);
    const char *c12[] = {"view", "obj", "remove", "g2", NULL};
    ASSERT(run_view(gedp, 4, c12) == BRLCAD_OK);

    bu_vls_free(&gpath);
    ged_close(gedp);

    bu_log("view_command: %d checks, %d failures\n", nchecks, nfails);
    return nfails ? EXIT_FAILURE : EXIT_SUCCESS;
}
