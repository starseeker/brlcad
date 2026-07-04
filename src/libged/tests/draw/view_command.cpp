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
#include <rt/view.h>

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
assert_feature_point_count(void *view, const char *name, size_t expected, int line)
{
    nchecks++;
    point_t *points = NULL;
    size_t point_count = 0;
    int copied = ged_draw_view_context_feature_points_copy(view, name,
	    &points, &point_count);
    if (!copied || point_count != expected) {
	bu_log("FAIL [%s:%d] %s point count expected %zu, got %zu\n",
		__FILE__, line, name, expected, point_count);
	nfails++;
    }
    if (points)
	bu_free(points, "view command test copied feature points");
}

#define ASSERT_FEATURE_POINT_COUNT(view, name, expected) assert_feature_point_count((view), (name), (expected), __LINE__)

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

    struct rt_view_render_export_consistency consistency =
	RT_VIEW_RENDER_EXPORT_CONSISTENCY_INIT;
    ASSERT(rt_view_context_render_export_consistency(v, "all.g",
	    &consistency));
    ASSERT(consistency.export_record_found);
    ASSERT(consistency.render_item_found);
    ASSERT(consistency.export_render_consistent);
    ASSERT(consistency.backend_node_found);
    ASSERT(consistency.export_backend_consistent);

    void *pick = rt_view_context_pick_semantic_path(v, "all.g");
    struct bu_vls pick_path = BU_VLS_INIT_ZERO;
    ASSERT(pick != NULL);
    ASSERT(rt_view_pick_result_context_count(pick) > 0);
    if (pick && rt_view_pick_result_context_count(pick) > 0) {
	ASSERT(rt_view_pick_result_context_path(pick, 0, &pick_path));
	ASSERT(BU_STR_EQUAL(bu_vls_cstr(&pick_path), "all.g"));
    }

    point_t sample = VINIT_ZERO;
    void *snap = rt_view_snap_result_context_create();
    ASSERT(snap != NULL);
    int snap_count = snap ? rt_view_context_snap_candidates_result(v,
	    sample, 1.0, RT_VIEW_SNAP_KIND_ENDPOINT, snap) : 0;
    ASSERT(snap_count >= 0);
    if (snap_count > 0 && rt_view_snap_result_context_count(snap) > 0) {
	struct bu_vls snap_path = BU_VLS_INIT_ZERO;
	ASSERT(rt_view_snap_result_context_source_path(snap, 0, &snap_path));
	ASSERT(bu_vls_strlen(&snap_path) > 0);
	bu_vls_free(&snap_path);
    }

    point_t a = VINIT_ZERO;
    point_t b = {1.0, 0.0, 0.0};
    struct rt_view_measure_result measure = RT_VIEW_MEASURE_RESULT_INIT;
    ASSERT(rt_view_context_measure_candidates(v, a, b, &measure) == 1);
    ASSERT(measure.valid);
    ASSERT(fabs(measure.distance - 1.0) < 1.0e-9);

    rt_view_snap_result_context_free(snap);
    bu_vls_free(&pick_path);
    rt_view_pick_result_context_free(pick);
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
    ASSERT(rt_view_set_context_remove(view_set_ctx, NULL));
    void *views[2] = {NULL, NULL};
    for (int i = 0; i < 2; i++) {
	char view_name[16];
	snprintf(view_name, sizeof(view_name), "V%d", i);
	views[i] = rt_view_context_create();
	ASSERT(views[i] != NULL);
	ASSERT(rt_view_context_name_set(views[i], view_name));
	ASSERT(rt_view_context_dimensions_set(views[i], 640, 480));
	ASSERT(rt_view_set_context_add(view_set_ctx, views[i]));
	ged_view_context_owned_add(gedp, views[i]);
	if (!i)
	    ged_view_active_ctx_set(gedp, views[i]);
    }

    test_command_report_record_consistency(gedp, views[0]);

    const char *c0[] = {"view", "annotation", "line", "create", "u_line", "0", "0", "0", "1", "0", "0", NULL};
    ASSERT(run_view(gedp, 11, c0) == BRLCAD_OK);
    ASSERT_FEATURE_POINT_COUNT(views[0], "u_line", 2);

    const char *lc0[] = {"view", "annotation", "line", "create", "u_line_edit", "0", "0", "0", "1", "0", "0", NULL};
    ASSERT_VIEW_OK(gedp, 11, lc0);
    ASSERT_FEATURE_POINT_COUNT(views[0], "u_line_edit", 2);
    const char *lc1[] = {"view", "annotation", "line", "append", "u_line_edit", "2", "0", "0", "3", "0", "0", NULL};
    ASSERT_VIEW_OK(gedp, 11, lc1);
    ASSERT_FEATURE_POINT_COUNT(views[0], "u_line_edit", 4);
    const char *lc2[] = {"view", "annotation", "line", "remove", "u_line_edit", "1", "2", NULL};
    ASSERT_VIEW_OK(gedp, 7, lc2);
    ASSERT_FEATURE_POINT_COUNT(views[0], "u_line_edit", 2);
    const char *lc3[] = {"view", "annotation", "line", "clear", "u_line_edit", NULL};
    ASSERT_VIEW_OK(gedp, 5, lc3);
    ASSERT(ged_draw_view_context_feature_exists(views[0], "u_line_edit") == 0);

    const char *c1[] = {"view", "object", "info", "u_line", "type", NULL};
    ASSERT(run_view(gedp, 5, c1) == BRLCAD_OK);

    const char *c2[] = {"view", "object", "hide", "u_line", NULL};
    ASSERT(run_view(gedp, 4, c2) == BRLCAD_OK);
    const char *c3[] = {"view", "object", "info", "u_line", "visible", NULL};
    ASSERT(run_view(gedp, 5, c3) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("DOWN") != std::string::npos);
    const char *c3a[] = {"view", "object", "show", "u_line", NULL};
    ASSERT(run_view(gedp, 4, c3a) == BRLCAD_OK);
    const char *c3b[] = {"view", "object", "info", "u_line", "visible", NULL};
    ASSERT(run_view(gedp, 5, c3b) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("UP") != std::string::npos);

    const char *c3c[] = {"view", "object", "style", "set", "u_line", "color", "10/20/30", NULL};
    ASSERT(run_view(gedp, 7, c3c) == BRLCAD_OK);
    const char *c3d[] = {"view", "object", "info", "u_line", "color", NULL};
    ASSERT(run_view(gedp, 5, c3d) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("10/20/30") != std::string::npos);

    const char *c4[] = {"view", "object", "list", "u_*", NULL};
    ASSERT(run_view(gedp, 4, c4) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("u_line") != std::string::npos);

    const char *c5[] = {"view", "object", "style", "set", "u_line", "arrow", "1", NULL};
    ASSERT(run_view(gedp, 7, c5) == BRLCAD_OK);
    const char *c5a[] = {"view", "object", "realize", "u_line", "12", "34", NULL};
    ASSERT(run_view(gedp, 6, c5a) == BRLCAD_OK);
    const char *c5b[] = {"view", "vZ", "-N", "u_line", NULL};
    ASSERT(run_view(gedp, 4, c5b) == BRLCAD_OK);
    ASSERT(!result_str(gedp).empty());

    const char *a0[] = {"view", "annotation", "axes", "create", "u_axes", "1", "2", "3", NULL};
    ASSERT_VIEW_OK(gedp, 8, a0);
    const char *a1[] = {"view", "annotation", "axes", "pos", "u_axes", NULL};
    ASSERT_VIEW_OK(gedp, 5, a1);
    ASSERT(result_str(gedp).find("1.000000 2.000000 3.000000") != std::string::npos);
    const char *a2[] = {"view", "annotation", "axes", "pos", "u_axes", "4", "5", "6", NULL};
    ASSERT_VIEW_OK(gedp, 8, a2);
    const char *a3[] = {"view", "annotation", "axes", "size", "u_axes", "12.5", NULL};
    ASSERT_VIEW_OK(gedp, 6, a3);
    const char *a4[] = {"view", "annotation", "axes", "size", "u_axes", NULL};
    ASSERT_VIEW_OK(gedp, 5, a4);
    ASSERT(result_str(gedp).find("12.500000") != std::string::npos);
    const char *a5[] = {"view", "annotation", "axes", "line_width", "u_axes", "3", NULL};
    ASSERT_VIEW_OK(gedp, 6, a5);
    const char *a6[] = {"view", "annotation", "axes", "line_width", "u_axes", NULL};
    ASSERT_VIEW_OK(gedp, 5, a6);
    ASSERT(result_str(gedp).find("3") != std::string::npos);
    const char *a7[] = {"view", "annotation", "axes", "axes_color", "u_axes", "10", "20", "30", NULL};
    ASSERT_VIEW_OK(gedp, 8, a7);
    const char *a8[] = {"view", "annotation", "axes", "axes_color", "u_axes", NULL};
    ASSERT_VIEW_OK(gedp, 5, a8);
    ASSERT(result_str(gedp).find("10 20 30") != std::string::npos);

    const char *p0[] = {"view", "polygon", "create", "u_poly", "10", "10", NULL};
    ASSERT_VIEW_OK(gedp, 6, p0);
    const char *p1[] = {"view", "polygon", "append", "u_poly", "30", "10", NULL};
    ASSERT_VIEW_OK(gedp, 6, p1);
    const char *p2[] = {"view", "polygon", "append", "u_poly", "30", "30", NULL};
    ASSERT_VIEW_OK(gedp, 6, p2);
    const char *p3[] = {"view", "polygon", "append", "u_poly", "10", "30", NULL};
    ASSERT_VIEW_OK(gedp, 6, p3);
    const char *p4[] = {"view", "polygon", "close", "u_poly", NULL};
    ASSERT_VIEW_OK(gedp, 4, p4);
    const char *p5[] = {"view", "polygon", "fill", "u_poly", "1", "0", "0.1", NULL};
    ASSERT_VIEW_OK(gedp, 7, p5);
    const char *p6[] = {"view", "polygon", "fill_color", "u_poly", "10/20/30", NULL};
    ASSERT_VIEW_OK(gedp, 5, p6);
    const char *p7[] = {"view", "polygon", "area", "u_poly", NULL};
    ASSERT_VIEW_OK(gedp, 4, p7);
    ASSERT(!result_str(gedp).empty());
    rt_view_polygon_ref poly_ref =
	rt_view_context_polygon_find(views[0], "u_poly");
    ASSERT(!rt_view_polygon_ref_is_null(poly_ref));
    struct rt_view_polygon_record poly_rec = {};
    ASSERT(rt_view_polygon_record_get(poly_ref, &poly_rec));
    ASSERT(poly_rec.type == RT_VIEW_POLYGON_GENERAL);
    ASSERT(poly_rec.contour_count == 1);
    ASSERT(poly_rec.point_count == 4);
    ASSERT(poly_rec.first_contour_open == 0);
    ASSERT(poly_rec.fill_flag == 1);

    const char *c6[] = {"view", "-V", "V0", "annotation", "-L", "line", "create", "l_line", "0", "0", "0", NULL};
    ASSERT(run_view(gedp, 11, c6) == BRLCAD_OK);
    const char *c7[] = {"view", "-V", "V0", "object", "list", NULL};
    ASSERT(run_view(gedp, 5, c7) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("l_line") != std::string::npos);
    const char *c8[] = {"view", "-V", "V1", "object", "list", NULL};
    ASSERT(run_view(gedp, 5, c8) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("l_line") == std::string::npos);

    const char *c11[] = {"view", "db", "add", "all.g", "--as", "g2", NULL};
    ASSERT(run_view(gedp, 6, c11) == BRLCAD_OK);
    const char *c12[] = {"view", "db", "delete", "g2", NULL};
    ASSERT(run_view(gedp, 4, c12) == BRLCAD_OK);

    bu_vls_free(&gpath);
    ged_close(gedp);

    bu_log("view_command: %d checks, %d failures\n", nchecks, nfails);
    return nfails ? EXIT_FAILURE : EXIT_SUCCESS;
}
