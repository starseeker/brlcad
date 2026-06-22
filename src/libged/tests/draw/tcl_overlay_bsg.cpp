/*              T C L _ O V E R L A Y _ B S G . C P P
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
/** @file tcl_overlay_bsg.cpp
 *
 * Phase T1/T4 (drawing_stack_modernization) — structural regression test.
 *
 * Verifies that the BSG view-scope feature lifecycle used by the T1
 * sync helpers (arrows, lines, labels, axes, polygons) behaves
 * correctly:
 *
 *   • bsg_feature_create_lines creates a named local-scope feature
 *   • bsg_feature_find         locates it by name
 *   • bsg_feature_visit        reports stable feature records
 *   • bsg_feature_remove       deletes it; subsequent find returns NULL
 *   • dm_draw_objs              can be called headlessly (NULL dmp) without
 *                               crashing — the T2-final call site in
 *                               go_refresh_draw must be no-op safe.
 *
 * This test does NOT attach a display manager.  It is intentionally headless
 * so that the BSG feature lifecycle contract is enforced independently of
 * any rendering backend.
 *
 * Usage: ged_test_tcl_overlay_bsg <directory-with-moss.g>
 */

#include "common.h"

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fstream>

#include <bu.h>
#include "bg/line_layer.h"
#include "bsg/feature.h"
#include "bsg/geometry.h"
#include "bsg/tcl_data.h"
#include "bsg/util.h"
#include "bsg/scene_object.h"
#include "bsg/draw_source.h"
#include "bsg/view_state.h"
#include <dm.h>
#include <ged.h>

#include "../../bsg_ged_draw_view_private.h"

#define ASSERT(cond) do { \
    nchecks++; \
    if (!(cond)) { \
	bu_log("FAIL [%s:%d] %s\n", __FILE__, __LINE__, #cond); \
	nfails++; \
    } \
} while (0)

static int nchecks = 0;
static int nfails  = 0;

/* Visitor that counts features reached via bsg_feature_visit. */
static int
_count_visit_cb(bsg_feature_ref /*ref*/, const struct bsg_feature_record *record, void *ud)
{
    int *cnt = (int *)ud;
    ASSERT(record != NULL);
    (*cnt)++;
    return 1;
}

int
main(int ac, char *av[])
{
    bu_setprogname(av[0]);

    if (ac != 2) {
	bu_exit(EXIT_FAILURE, "Usage: %s <directory-containing-moss.g>\n", av[0]);
    }

    /* ------------------------------------------------------------------ *
     * Open a headless GED session.                                        *
     * ------------------------------------------------------------------ */
    struct bu_vls fname = BU_VLS_INIT_ZERO;
    struct bu_vls moss = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&moss, "%s/moss.g", av[1]);
    char tmpname[MAXPATHLEN];
    FILE *fp = bu_temp_file(tmpname, MAXPATHLEN);
    if (!fp) {
	bu_log("failed to create temp db path: %s\n", std::strerror(errno));
	bu_vls_free(&moss);
	bu_vls_free(&fname);
	return 1;
    }
    fclose(fp);
    bu_vls_sprintf(&fname, "%s", tmpname);
    {
	/* This test is headless, but ged_open still requires a valid .g file. */
	std::ifstream orig(bu_vls_cstr(&moss), std::ios::binary);
	std::ofstream tmpg(bu_vls_cstr(&fname), std::ios::binary);
	if (!orig.good() || !tmpg.good()) {
	    bu_log("failed to prepare tmp db: %s\n", bu_vls_cstr(&fname));
	    bu_vls_free(&moss);
	    bu_vls_free(&fname);
	    return 1;
	}
	tmpg << orig.rdbuf();
	orig.close();
	tmpg.close();
    }
    struct ged *gedp = ged_open("db", bu_vls_cstr(&fname), 1);
    bu_vls_free(&moss);
    if (!gedp) {
	bu_log("ged_open failed\n");
	bu_file_delete(bu_vls_cstr(&fname));
	bu_vls_free(&fname);
	return 1;
    }

    bu_log("=== TCL feature BSG lifecycle ===\n");

    struct bsg_view *v = gedp->ged_gvp;
    ASSERT(v != NULL);
    ASSERT(bsg_view_scene_attached(v));

    /* ------------------------------------------------------------------ *
 * [1] create: bsg_feature_create_lines must return a non-NULL scene  *
     *     ref and register it in the local scope.                        *
     * ------------------------------------------------------------------ */
    bu_log("[1] bsg_feature_create_lines...\n");
    const char *tname = "_tcl_test_feature";
    bsg_feature_ref obj_ref = bsg_feature_create_lines(v, tname, 1 /*local*/);
    ASSERT(!bsg_feature_ref_is_null(obj_ref));
    if (bsg_feature_ref_is_null(obj_ref)) goto done;

    /* ------------------------------------------------------------------ *
     * [2] find: bsg_feature_find must locate the object by name.         *
     * ------------------------------------------------------------------ */
    bu_log("[2] bsg_feature_find...\n");
    {
	bsg_feature_ref found_ref = bsg_feature_find(v, tname);
	ASSERT(!bsg_feature_ref_is_null(found_ref));
	ASSERT(found_ref.token == obj_ref.token);
	struct bsg_feature_record record;
	ASSERT(bsg_feature_record_get(found_ref, &record) == 1);
	ASSERT(record.family == BSG_FEATURE_LINES);
	ASSERT(record.scope == BSG_FEATURE_SCOPE_LOCAL);
    }

    /* ------------------------------------------------------------------ *
     * [3] feature geometry: add typed line data and verify non-empty.     *
     * ------------------------------------------------------------------ */
    bu_log("[3] bsg_feature_points_replace...\n");
    {
	point_t points[2] = {{0, 0, 0}, {1, 0, 0}};
	int cmds[2] = {BSG_GEOMETRY_LINE_MOVE, BSG_GEOMETRY_LINE_DRAW};
	ASSERT(bsg_feature_points_replace(obj_ref, BSG_FEATURE_LINES, points, cmds, 2) == 1);
	point_t *copied_points = NULL;
	int *copied_cmds = NULL;
	size_t point_count = 0;
	ASSERT(bsg_feature_points_copy(obj_ref, &copied_points, &copied_cmds, &point_count) == 1);
	ASSERT(point_count == 2);
	if (point_count == 2 && copied_points && copied_cmds) {
	    ASSERT(VNEAR_EQUAL(copied_points[0], points[0], VUNITIZE_TOL));
	    ASSERT(VNEAR_EQUAL(copied_points[1], points[1], VUNITIZE_TOL));
	    ASSERT(copied_cmds[0] == BSG_GEOMETRY_LINE_MOVE);
	    ASSERT(copied_cmds[1] == BSG_GEOMETRY_LINE_DRAW);
	}
	if (copied_points)
	    bu_free(copied_points, "feature test copied points");
	if (copied_cmds)
	    bu_free(copied_cmds, "feature test copied commands");
    }

    /* ------------------------------------------------------------------ *
     * [4] set_color / set_line_width / set_visible typed setters.        *
     * ------------------------------------------------------------------ */
    bu_log("[4] typed setters...\n");
    {
	bsg_feature_set_color(obj_ref, 255, 128, 0);
	bsg_feature_set_line_width(obj_ref, 2);
	bsg_feature_set_visible(obj_ref, 1);
	struct bsg_feature_style style = BSG_FEATURE_STYLE_INIT;
	ASSERT(bsg_feature_style_get(obj_ref, &style) == 1);
	ASSERT(style.color_valid == 1);
	ASSERT(style.color[0] == 255 && style.color[1] == 128 && style.color[2] == 0);
	ASSERT(style.line_width == 2);
	ASSERT(style.visible == 1);
    }

    /* ------------------------------------------------------------------ *
     * [5] GED data-arrow facade: replace, style, tip, and removal.       *
     * ------------------------------------------------------------------ */
    bu_log("[5] GED data-arrow facade...\n");
    {
	const char *aname = "_tcl_test_arrows_adapter";
	point_t arrow_points[2] = {{0, 0, 0}, {0, 1, 0}};
	struct ged_draw_view_feature_style arrow_style = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	arrow_style.visible = 1;
	arrow_style.color_valid = 1;
	arrow_style.color[0] = 10;
	arrow_style.color[1] = 20;
	arrow_style.color[2] = 30;
	arrow_style.line_width = 3;
	arrow_style.arrow = 1;
	arrow_style.arrow_tip_length = 12.0;
	arrow_style.arrow_tip_width = 5.0;

	ASSERT(ged_draw_view_tcl_arrows_replace(v, aname, arrow_points, 2, &arrow_style) == 1);
	ASSERT(ged_draw_view_feature_exists(v, aname) == 1);

	point_t *copied_points = NULL;
	size_t point_count = 0;
	ASSERT(ged_draw_view_feature_points_copy(v, aname, &copied_points, &point_count) == 1);
	ASSERT(point_count == 2);
	if (point_count == 2 && copied_points) {
	    ASSERT(VNEAR_EQUAL(copied_points[0], arrow_points[0], VUNITIZE_TOL));
	    ASSERT(VNEAR_EQUAL(copied_points[1], arrow_points[1], VUNITIZE_TOL));
	}
	if (copied_points)
	    bu_free(copied_points, "GED data-arrow facade copied points");

	struct ged_draw_view_feature_style current = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	ASSERT(ged_draw_view_feature_style_get(v, aname, &current) == 1);
	ASSERT(current.visible == 1);
	ASSERT(current.color_valid == 1);
	ASSERT(current.color[0] == 10 && current.color[1] == 20 && current.color[2] == 30);
	ASSERT(current.line_width == 3);

	fastf_t tip_length = 0.0;
	fastf_t tip_width = 0.0;
	ASSERT(ged_draw_view_arrow_tip_get(v, aname, &tip_length, &tip_width) == 1);
	ASSERT(NEAR_EQUAL(tip_length, 12.0, SMALL_FASTF));
	ASSERT(NEAR_EQUAL(tip_width, 5.0, SMALL_FASTF));

	ASSERT(ged_draw_view_arrow_tip_set(v, aname, 14.0, 6.0) == 1);
	ASSERT(ged_draw_view_arrow_tip_get(v, aname, &tip_length, &tip_width) == 1);
	ASSERT(NEAR_EQUAL(tip_length, 14.0, SMALL_FASTF));
	ASSERT(NEAR_EQUAL(tip_width, 6.0, SMALL_FASTF));

	ASSERT(ged_draw_view_feature_visible_set(v, aname, 0) == 1);
	ASSERT(ged_draw_view_feature_visible(v, aname) == 0);
	ASSERT(ged_draw_view_line_color_set(v, aname, 1, 2, 3) == 1);
	ASSERT(ged_draw_view_line_width_set(v, aname, 4) == 1);
	ASSERT(ged_draw_view_feature_style_get(v, aname, &current) == 1);
	ASSERT(current.color[0] == 1 && current.color[1] == 2 && current.color[2] == 3);
	ASSERT(current.line_width == 4);

	ASSERT(ged_draw_view_feature_remove(v, aname) == 1);
	ASSERT(ged_draw_view_feature_exists(v, aname) == 0);
    }

    /* ------------------------------------------------------------------ *
     * [6] GED data-axes facade: replace, centers, style, and removal.    *
     * ------------------------------------------------------------------ */
    bu_log("[6] GED data-axes facade...\n");
    {
	const char *xname = "_tcl_test_axes_adapter";
	point_t axes_centers[2] = {{1, 2, 3}, {4, 5, 6}};
	struct ged_draw_view_feature_style axes_style = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	axes_style.visible = 1;
	axes_style.color_valid = 1;
	axes_style.color[0] = 40;
	axes_style.color[1] = 50;
	axes_style.color[2] = 60;
	axes_style.line_width = 7;

	ASSERT(ged_draw_view_tcl_axes_replace(v, xname, axes_centers, 2, 2.5, &axes_style) == 1);
	ASSERT(ged_draw_view_feature_exists(v, xname) == 1);

	point_t *copied_centers = NULL;
	size_t center_count = 0;
	ASSERT(ged_draw_view_feature_axes_centers_copy(v, xname, &copied_centers, &center_count) == 1);
	ASSERT(center_count == 2);
	if (center_count == 2 && copied_centers) {
	    ASSERT(VNEAR_EQUAL(copied_centers[0], axes_centers[0], VUNITIZE_TOL));
	    ASSERT(VNEAR_EQUAL(copied_centers[1], axes_centers[1], VUNITIZE_TOL));
	}
	if (copied_centers)
	    bu_free(copied_centers, "GED data-axes facade copied centers");

	point_t *copied_points = NULL;
	size_t point_count = 0;
	ASSERT(ged_draw_view_feature_points_copy(v, xname, &copied_points, &point_count) == 1);
	ASSERT(point_count == 12);
	if (point_count >= 2 && copied_points)
	    ASSERT(NEAR_EQUAL((copied_points[1][X] - copied_points[0][X]) * 0.5, 2.5, SMALL_FASTF));
	if (copied_points)
	    bu_free(copied_points, "GED data-axes facade copied points");

	struct ged_draw_view_feature_style current = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	ASSERT(ged_draw_view_feature_style_get(v, xname, &current) == 1);
	ASSERT(current.visible == 1);
	ASSERT(current.color_valid == 1);
	ASSERT(current.color[0] == 40 && current.color[1] == 50 && current.color[2] == 60);
	ASSERT(current.line_width == 7);

	ASSERT(ged_draw_view_feature_visible_set(v, xname, 0) == 1);
	ASSERT(ged_draw_view_feature_visible(v, xname) == 0);
	ASSERT(ged_draw_view_line_color_set(v, xname, 4, 5, 6) == 1);
	ASSERT(ged_draw_view_line_width_set(v, xname, 8) == 1);
	ASSERT(ged_draw_view_feature_style_get(v, xname, &current) == 1);
	ASSERT(current.color[0] == 4 && current.color[1] == 5 && current.color[2] == 6);
	ASSERT(current.line_width == 8);

	ASSERT(ged_draw_view_feature_remove(v, xname) == 1);
	ASSERT(ged_draw_view_feature_exists(v, xname) == 0);
    }

    /* ------------------------------------------------------------------ *
     * [7] GED Tcl polygon facade: typed outline replacement.             *
     * ------------------------------------------------------------------ */
    bu_log("[7] GED Tcl polygon facade...\n");
    {
	const char *pname = "_tcl_test_polygon_adapter";
	point_t polygon_points[4] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 0, 0}};
	int polygon_cmds[4] = {
	    BG_GEOMETRY_LINE_MOVE,
	    BG_GEOMETRY_LINE_DRAW,
	    BG_GEOMETRY_LINE_DRAW,
	    BG_GEOMETRY_LINE_DRAW
	};
	struct ged_draw_view_feature_style polygon_style = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	polygon_style.visible = 1;
	polygon_style.color_valid = 1;
	polygon_style.color[0] = 70;
	polygon_style.color[1] = 80;
	polygon_style.color[2] = 90;
	polygon_style.line_width = 2;
	polygon_style.line_style = 1;

	ASSERT(ged_draw_view_tcl_polygons_replace(v, pname, polygon_points,
		polygon_cmds, 4, &polygon_style) == 1);
	ASSERT(ged_draw_view_feature_exists(v, pname) == 1);

	bsg_feature_ref pref = bsg_feature_find(v, pname);
	point_t *copied_points = NULL;
	int *copied_cmds = NULL;
	size_t point_count = 0;
	ASSERT(bsg_feature_points_copy(pref, &copied_points, &copied_cmds, &point_count) == 1);
	ASSERT(point_count == 4);
	if (point_count == 4 && copied_points && copied_cmds) {
	    ASSERT(VNEAR_EQUAL(copied_points[0], polygon_points[0], VUNITIZE_TOL));
	    ASSERT(VNEAR_EQUAL(copied_points[2], polygon_points[2], VUNITIZE_TOL));
	    ASSERT(copied_cmds[0] == BSG_GEOMETRY_LINE_MOVE);
	    ASSERT(copied_cmds[1] == BSG_GEOMETRY_LINE_DRAW);
	    ASSERT(copied_cmds[3] == BSG_GEOMETRY_LINE_DRAW);
	}
	if (copied_points)
	    bu_free(copied_points, "GED Tcl polygon facade copied points");
	if (copied_cmds)
	    bu_free(copied_cmds, "GED Tcl polygon facade copied commands");

	struct ged_draw_view_feature_style current = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	ASSERT(ged_draw_view_feature_style_get(v, pname, &current) == 1);
	ASSERT(current.visible == 1);
	ASSERT(current.color_valid == 1);
	ASSERT(current.color[0] == 70 && current.color[1] == 80 && current.color[2] == 90);
	ASSERT(current.line_width == 2);
	ASSERT(current.line_style == 1);

	ASSERT(ged_draw_view_feature_remove(v, pname) == 1);
	ASSERT(ged_draw_view_feature_exists(v, pname) == 0);
    }

    /* ------------------------------------------------------------------ *
     * [8] GED Tcl label facade: replace, copy, point update, removal.    *
     * ------------------------------------------------------------------ */
    bu_log("[8] GED Tcl label facade...\n");
    {
	const char *lname = "_tcl_test_labels_adapter";
	struct ged_draw_view_label_data labels[2] = {
	    GED_DRAW_VIEW_LABEL_DATA_INIT,
	    GED_DRAW_VIEW_LABEL_DATA_INIT
	};
	labels[0].text = "label one";
	VSET(labels[0].point, 1.0, 2.0, 3.0);
	labels[0].color_valid = 1;
	labels[0].color[0] = 11;
	labels[0].color[1] = 22;
	labels[0].color[2] = 33;
	labels[1].text = "label two";
	VSET(labels[1].point, 4.0, 5.0, 6.0);
	labels[1].color_valid = 1;
	labels[1].color[0] = 44;
	labels[1].color[1] = 55;
	labels[1].color[2] = 66;

	ASSERT(ged_draw_view_tcl_labels_replace(v, lname, 1, labels, 2) == 1);
	ASSERT(ged_draw_view_feature_exists(v, lname) == 1);
	ASSERT(ged_draw_view_label_count(v, lname) == 2);

	struct bu_vls text = BU_VLS_INIT_ZERO;
	point_t copied_point = VINIT_ZERO;
	unsigned char rgb[3] = {0, 0, 0};
	ASSERT(ged_draw_view_label_copy(v, lname, 0, &text, copied_point, rgb) == 1);
	ASSERT(BU_STR_EQUAL(bu_vls_cstr(&text), "label one"));
	ASSERT(VNEAR_EQUAL(copied_point, labels[0].point, VUNITIZE_TOL));
	ASSERT(rgb[0] == 11 && rgb[1] == 22 && rgb[2] == 33);
	bu_vls_free(&text);

	point_t moved = {7.0, 8.0, 9.0};
	ASSERT(ged_draw_view_label_point_set(v, lname, 1, moved) == 1);
	ASSERT(ged_draw_view_label_copy(v, lname, 1, NULL, copied_point, rgb) == 1);
	ASSERT(VNEAR_EQUAL(copied_point, moved, VUNITIZE_TOL));
	ASSERT(rgb[0] == 44 && rgb[1] == 55 && rgb[2] == 66);

	ASSERT(ged_draw_view_tcl_labels_replace(v, lname, 0, labels, 2) == 1);
	ASSERT(ged_draw_view_feature_exists(v, lname) == 0);
	ASSERT(ged_draw_view_tcl_labels_replace(v, lname, 1, NULL, 0) == 1);
	ASSERT(ged_draw_view_feature_exists(v, lname) == 0);
    }

    /* ------------------------------------------------------------------ *
     * [9] visit: bsg_feature_visit with BSG_FEATURE_SCOPE_LOCAL must      *
     *     reach at least the one object we created.                      *
     * ------------------------------------------------------------------ */
    bu_log("[9] bsg_feature_visit (local scope)...\n");
    {
	int cnt = 0;
	bsg_feature_visit(v, BSG_FEATURE_SCOPE_LOCAL, _count_visit_cb, &cnt);
	ASSERT(cnt >= 1);
    }

    /* ------------------------------------------------------------------ *
     * [10] dm_draw_objs headless: must not crash when dmp is NULL.       *
     *     This mirrors the T2-final call in go_refresh_draw.             *
     * ------------------------------------------------------------------ */
    bu_log("[10] dm_draw_objs headless (NULL dmp)...\n");
    {
	struct dm *saved_dmp = (struct dm *)v->dmp;
	v->dmp = NULL;
	dm_draw_objs(v);   /* must be a no-op, not a crash */
	v->dmp = saved_dmp;
    }

    /* ------------------------------------------------------------------ *
     * [11] remove: bsg_feature_remove must delete the named object;      *
     *     a subsequent find must return NULL.                            *
     * ------------------------------------------------------------------ */
    bu_log("[11] bsg_feature_remove...\n");
    {
	int r = bsg_feature_remove(v, tname);
	ASSERT(r == 1);
	ASSERT(bsg_feature_ref_is_null(bsg_feature_find(v, tname)));
    }

    /* ------------------------------------------------------------------ *
     * [12] remove idempotency: removing a non-existent name is safe.     *
     * ------------------------------------------------------------------ */
    bu_log("[12] remove idempotency...\n");
    {
	int r = bsg_feature_remove(v, "_tcl_test_nonexistent");
	(void)r; /* return value may differ by impl; the call must not crash */
    }

    /* ------------------------------------------------------------------ *
     * [13] multi-object: create data/sdata objects for all four overlay *
     *     slots (arrows, lines, labels, polygons) and verify they are    *
     *     independently addressable and removable.                       *
     * ------------------------------------------------------------------ */
    bu_log("[13] multi-object feature slots...\n");
    {
	const char *slots[] = {
	    "_tcl_data_arrows",    "_tcl_sdata_arrows",
	    "_tcl_data_lines",     "_tcl_sdata_lines",
	    "_tcl_data_labels",    "_tcl_sdata_labels",
	    "_tcl_data_axes",      "_tcl_sdata_axes",
	    "_tcl_data_polygons",  "_tcl_sdata_polygons",
	    NULL
	};
	/* Create all slots */
	for (int k = 0; slots[k]; k++) {
	    bsg_feature_ref ref = bsg_feature_create_lines(v, slots[k], 1);
	    ASSERT(!bsg_feature_ref_is_null(ref));
	}
	/* Verify all are findable */
	for (int k = 0; slots[k]; k++) {
	    ASSERT(!bsg_feature_ref_is_null(bsg_feature_find(v, slots[k])));
	}
	/* Remove all slots */
	for (int k = 0; slots[k]; k++) {
	    bsg_feature_remove(v, slots[k]);
	    ASSERT(bsg_feature_ref_is_null(bsg_feature_find(v, slots[k])));
	}
    }

done:
    ged_close(gedp);
    bu_file_delete(bu_vls_cstr(&fname));
    bu_vls_free(&fname);

    bu_log("Result: %d checks, %d failures\n", nchecks, nfails);
    return (nfails > 0) ? 1 : 0;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
