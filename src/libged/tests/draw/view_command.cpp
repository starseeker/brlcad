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

#include <Inventor/SoViewport.h>
#include "brlobol/measure_action.h"
#include "brlobol/view_controller.h"
#include <bu.h>
#include <ged.h>
#include <ged/draw.h>
#include <ged/draw_obol.h>
#include "view_test_util.h"

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
assert_feature_overlay(void *view, const char *name, int expected, int line)
{
    nchecks++;
    struct ged_draw_view_feature_summary summary =
	GED_DRAW_VIEW_FEATURE_SUMMARY_INIT;
    int copied = ged_draw_view_context_feature_summary(view, name, &summary);
    if (!copied || !summary.exists || summary.is_overlay != expected) {
	bu_log("FAIL [%s:%d] %s overlay expected %d, got exists=%d overlay=%d\n",
		__FILE__, line, name, expected, summary.exists,
		summary.is_overlay);
	nfails++;
    }
}

#define ASSERT_FEATURE_OVERLAY(view, name, expected) assert_feature_overlay((view), (name), (expected), __LINE__)

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

    struct ged_draw_render_export_consistency consistency =
	GED_DRAW_RENDER_EXPORT_CONSISTENCY_INIT;
    ASSERT(ged_draw_view_context_render_export_consistency(gedp, v, "all.g",
	    &consistency));
    ASSERT(consistency.export_record_found);
    ASSERT(consistency.render_item_found);
    ASSERT(consistency.export_render_consistent);
    ASSERT(consistency.backend_node_found);
    ASSERT(consistency.export_backend_consistent);

    struct ged_draw_pick_result *pick =
	ged_draw_view_context_pick_semantic_path(gedp, v, "all.g");
    struct bu_vls pick_path = BU_VLS_INIT_ZERO;
    ASSERT(pick != NULL);
    ASSERT(ged_draw_pick_result_count(pick) > 0);
    if (pick && ged_draw_pick_result_count(pick) > 0) {
	ASSERT(ged_draw_pick_result_path(pick, 0, &pick_path));
	ASSERT(BU_STR_EQUAL(bu_vls_cstr(&pick_path), "all.g"));
    }

    point_t sample = VINIT_ZERO;
    point_t snap_candidate = VINIT_ZERO;
    int snap_count = ged_draw_view_context_snap_first_candidate(v, sample,
	    GED_DRAW_VIEW_SNAP_ENDPOINT, snap_candidate);
    ASSERT(snap_count >= 0);
    if (snap_count > 0) {
	ASSERT(std::isfinite((double)snap_candidate[X]));
	ASSERT(std::isfinite((double)snap_candidate[Y]));
	ASSERT(std::isfinite((double)snap_candidate[Z]));
    }

    BRLObolViewController *controller = ged_draw_obol_controller(gedp);
    ASSERT(controller != NULL);
    if (controller && controller->getViewport() &&
	    controller->getViewport()->getRoot()) {
	SoBRLMeasureAction measure;
	measure.setGeometryPolicy(SoBRLMeasureAction::DISPLAY_LEVEL);
	measure.apply(controller->getViewport()->getRoot());
	ASSERT(measure.hasSegments());
	ASSERT(measure.getSegmentCount() > 0);
	ASSERT(measure.getTotalLength() > 0.0f);
	ASSERT(!measure.getBounds().isEmpty());
    }

    bu_vls_free(&pick_path);
    ged_draw_pick_result_free(pick);
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
    ASSERT(ged_view_set_context_remove(view_set_ctx, NULL));
    void *views[2] = {NULL, NULL};
    for (int i = 0; i < 2; i++) {
	char view_name[16];
	snprintf(view_name, sizeof(view_name), "V%d", i);
	views[i] = ged_view_context_create();
	ASSERT(views[i] != NULL);
	ASSERT(bv_name_set(DRAW_TEST_BV(views[i]), view_name));
	ASSERT(bv_dimensions_set(DRAW_TEST_BV(views[i]), 640, 480));
	ASSERT(ged_view_set_context_add(view_set_ctx, views[i]));
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
    ASSERT_FEATURE_OVERLAY(views[0], "u_line_edit", 1);
    const char *lc1[] = {"view", "annotation", "line", "append", "u_line_edit", "2", "0", "0", "3", "0", "0", NULL};
    ASSERT_VIEW_OK(gedp, 11, lc1);
    ASSERT_FEATURE_POINT_COUNT(views[0], "u_line_edit", 4);
    const char *lc2[] = {"view", "annotation", "line", "remove", "u_line_edit", "1", "2", NULL};
    ASSERT_VIEW_OK(gedp, 7, lc2);
    ASSERT_FEATURE_POINT_COUNT(views[0], "u_line_edit", 2);
    ASSERT_FEATURE_OVERLAY(views[0], "u_line_edit", 1);
    const char *lc3[] = {"view", "annotation", "line", "clear", "u_line_edit", NULL};
    ASSERT_VIEW_OK(gedp, 5, lc3);
    ASSERT(ged_draw_view_context_feature_exists(views[0], "u_line_edit") == 0);

    const char *c1[] = {"view", "feature", "info", "u_line", "type", NULL};
    ASSERT(run_view(gedp, 5, c1) == BRLCAD_OK);
    const char *c1a[] = {"view", "feature", "info", "u_line", "kind", NULL};
    ASSERT(run_view(gedp, 5, c1a) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("lines") != std::string::npos);
    const char *c1b[] = {"view", "feature", "info", "u_line", "scope", NULL};
    ASSERT(run_view(gedp, 5, c1b) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("shared") != std::string::npos);
    const char *c1c[] = {"view", "feature", "info", "u_line", "overlay_class", NULL};
    ASSERT(run_view(gedp, 5, c1c) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("user-annotation") != std::string::npos);
    const char *c1d[] = {"view", "feature", "info", "u_line", "lifecycle", NULL};
    ASSERT(run_view(gedp, 5, c1d) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("persistent") != std::string::npos);
    const char *c1e[] = {"view", "feature", "info", "u_line", "command_result", NULL};
    ASSERT(run_view(gedp, 5, c1e) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("0") != std::string::npos);

    const char *c2[] = {"view", "feature", "hide", "u_line", NULL};
    ASSERT(run_view(gedp, 4, c2) == BRLCAD_OK);
    const char *c3[] = {"view", "feature", "info", "u_line", "visible", NULL};
    ASSERT(run_view(gedp, 5, c3) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("DOWN") != std::string::npos);
    const char *c3a[] = {"view", "feature", "show", "u_line", NULL};
    ASSERT(run_view(gedp, 4, c3a) == BRLCAD_OK);
    const char *c3b[] = {"view", "feature", "info", "u_line", "visible", NULL};
    ASSERT(run_view(gedp, 5, c3b) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("UP") != std::string::npos);

    const char *c3c[] = {"view", "feature", "style", "set", "u_line", "color", "10/20/30", NULL};
    ASSERT(run_view(gedp, 7, c3c) == BRLCAD_OK);
    const char *c3d[] = {"view", "feature", "info", "u_line", "color", NULL};
    ASSERT(run_view(gedp, 5, c3d) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("10/20/30") != std::string::npos);

    const char *c4[] = {"view", "feature", "list", "u_*", NULL};
    ASSERT(run_view(gedp, 4, c4) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("u_line") != std::string::npos);

    const char *c5[] = {"view", "feature", "style", "set", "u_line", "arrow", "1", NULL};
    ASSERT(run_view(gedp, 7, c5) == BRLCAD_OK);
    const char *c5a[] = {"view", "feature", "realize", "u_line", "12", "34", NULL};
    ASSERT(run_view(gedp, 6, c5a) == BRLCAD_OK);
    const char *c5b[] = {"view", "vZ", "-N", "u_line", NULL};
    ASSERT(run_view(gedp, 4, c5b) == BRLCAD_OK);
    ASSERT(!result_str(gedp).empty());

    const char *l0[] = {"view", "annotation", "label", "create", "u_label", "note", "1", "2", "3", NULL};
    ASSERT_VIEW_OK(gedp, 9, l0);
    ASSERT(ged_draw_view_context_label_count(views[0], "u_label") == 1);
    const char *l1[] = {"view", "feature", "info", "u_label", "type", NULL};
    ASSERT_VIEW_OK(gedp, 5, l1);
    ASSERT(!result_str(gedp).empty());

    const char *dl0[] = {"data_lines", "points", "{0 0 0} {1 0 0}", NULL};
    ASSERT_VIEW_OK(gedp, 3, dl0);
    ASSERT_FEATURE_POINT_COUNT(views[0], "_tcl_data_lines", 2);
    const char *dl1[] = {"data_lines", "draw", NULL};
    ASSERT_VIEW_OK(gedp, 2, dl1);
    ASSERT(result_str(gedp).find("1") != std::string::npos);
    const char *dl2[] = {"data_lines", "color", "11", "22", "33", NULL};
    ASSERT_VIEW_OK(gedp, 5, dl2);
    const char *dl3[] = {"data_lines", "color", NULL};
    ASSERT_VIEW_OK(gedp, 2, dl3);
    ASSERT(result_str(gedp).find("11 22 33") != std::string::npos);
    const char *dl4[] = {"data_lines", "line_width", "4", NULL};
    ASSERT_VIEW_OK(gedp, 3, dl4);
    const char *dl5[] = {"data_lines", "line_width", NULL};
    ASSERT_VIEW_OK(gedp, 2, dl5);
    ASSERT(result_str(gedp).find("4") != std::string::npos);
    const char *dl6[] = {"data_lines", "draw", "0", NULL};
    ASSERT_VIEW_OK(gedp, 3, dl6);
    ASSERT(ged_draw_view_context_feature_exists(views[0], "_tcl_data_lines") == 0);

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
    ged_draw_view_polygon_ref poly_ref =
	ged_draw_view_context_polygon_find(views[0], "u_poly");
    ASSERT(!ged_draw_view_polygon_ref_is_null(poly_ref));
    struct ged_draw_view_polygon_record poly_rec = {};
    ASSERT(ged_draw_view_polygon_record_get(poly_ref, &poly_rec));
    ASSERT(poly_rec.type == GED_DRAW_VIEW_POLYGON_GENERAL);
    ASSERT(poly_rec.contour_count == 1);
    ASSERT(poly_rec.point_count == 4);
    ASSERT(poly_rec.first_contour_open == 0);
    ASSERT(poly_rec.fill_flag == 1);

    const char *c6[] = {"view", "-V", "V0", "annotation", "-L", "line", "create", "l_line", "0", "0", "0", NULL};
    ASSERT(run_view(gedp, 11, c6) == BRLCAD_OK);
    const char *c7[] = {"view", "-V", "V0", "feature", "list", NULL};
    ASSERT(run_view(gedp, 5, c7) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("l_line") != std::string::npos);
    const char *c7a[] = {"view", "-V", "V0", "feature", "info", "l_line", "scope", NULL};
    ASSERT(run_view(gedp, 7, c7a) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("local") != std::string::npos);
    const char *c8[] = {"view", "-V", "V1", "feature", "list", NULL};
    ASSERT(run_view(gedp, 5, c8) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("l_line") == std::string::npos);

    const char *c11[] = {"view", "db", "add", "all.g", "--as", "g2", NULL};
    ASSERT(run_view(gedp, 6, c11) == BRLCAD_OK);
    const char *c11a[] = {"view", "feature", "style", "set", "g2", "color", "20/30/40", NULL};
    ASSERT(run_view(gedp, 7, c11a) == BRLCAD_OK);
    const char *c11b[] = {"view", "feature", "style", "get", "g2", "color", NULL};
    ASSERT(run_view(gedp, 6, c11b) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("20/30/40") != std::string::npos);
    const char *c11c[] = {"view", "feature", "style", "set", "g2", "arrow", "1", NULL};
    ASSERT(run_view(gedp, 7, c11c) == BRLCAD_ERROR);
    const char *c12[] = {"view", "db", "delete", "g2", NULL};
    ASSERT(run_view(gedp, 4, c12) == BRLCAD_OK);

    bu_vls_free(&gpath);
    ged_close(gedp);

    bu_log("view_command: %d checks, %d failures\n", nchecks, nfails);
    return nfails ? EXIT_FAILURE : EXIT_SUCCESS;
}
