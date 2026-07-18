/*                     T E S T _ G Q A . C
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
/** @file test_gqa.c
 *
 * Extract plotting data from a gqa run
 *
 */

#include "common.h"

#include <stdio.h>
#include <string.h>
#include <bu.h>
#include "bg/plot3.h"
#include <ged.h>
#include <ged/draw.h>
#include <ged/draw_obol.h>
#include "wdb.h"

struct gqa_segment_writer {
    FILE *fp;
    size_t count;
};

struct gqa_record_capture {
    struct gqa_segment_writer *writer;
    const char *plot_fname;
    const char *feature_path;
    int rgb[3];
    int write_plot;
};

static int
write_gqa_segment(const point_t a, const point_t b, void *data)
{
    struct gqa_segment_writer *writer = (struct gqa_segment_writer *)data;
    if (writer->fp)
	pdv_3line(writer->fp, a, b);
    writer->count++;
    return 1;
}

static int
capture_gqa_record(const struct ged_draw_view_db_object_record *rec, void *data)
{
    struct gqa_record_capture *capture = (struct gqa_record_capture *)data;
    struct gqa_segment_writer *writer = capture ? capture->writer : NULL;

    if (!writer || !capture || !capture->feature_path || !rec ||
	    !rec->path || !strstr(rec->path, capture->feature_path))
	return 1;
    if (!rec->vlist_point_count)
	return 1;

    if (capture->write_plot && !writer->fp) {
	writer->fp = fopen(capture->plot_fname, "wb");
	if (!writer->fp)
	    bu_exit(EXIT_FAILURE, "Could not open %s for writing\n",
		    capture->plot_fname);
	pl_color(writer->fp, capture->rgb[0], capture->rgb[1],
		capture->rgb[2]);
    }
    printf("found %s;\n", rec->path);
    (void)ged_draw_view_db_object_record_foreach_segment(rec,
	    write_gqa_segment, writer);

    return 1;
}

static void
require_command_result_metadata(struct ged_view_context *view_ctx, const char *feature_name,
	const char *key, const char *value)
{
    size_t count = ged_view_feature_metadata_count(view_ctx,
	    feature_name);
    size_t i;
    for (i = 0; i < count; i++) {
	struct bu_vls found_key = BU_VLS_INIT_ZERO;
	struct bu_vls found_value = BU_VLS_INIT_ZERO;
	int copied = ged_view_feature_metadata_copy(view_ctx,
		feature_name, i, &found_key, &found_value);
	int matched = copied &&
	    BU_STR_EQUAL(bu_vls_cstr(&found_key), key) &&
	    BU_STR_EQUAL(bu_vls_cstr(&found_value), value);
	bu_vls_free(&found_key);
	bu_vls_free(&found_value);
	if (matched)
	    return;
    }
    bu_exit(EXIT_FAILURE, "%s metadata %s=%s was not published.\n",
	    feature_name, key, value);
}

static void
require_command_result_primitive_metadata(struct ged_view_context *view_ctx,
	const char *feature_name,
	int primitive,
	const char *key,
	const char *value)
{
    size_t count = ged_view_feature_primitive_metadata_count(
	view_ctx, feature_name, primitive);
    size_t i;
    for (i = 0; i < count; i++) {
	struct bu_vls found_key = BU_VLS_INIT_ZERO;
	struct bu_vls found_value = BU_VLS_INIT_ZERO;
	int copied = ged_view_feature_primitive_metadata_copy(
	    view_ctx, feature_name, primitive, i, &found_key, &found_value);
	int matched = copied &&
	    BU_STR_EQUAL(bu_vls_cstr(&found_key), key) &&
	    BU_STR_EQUAL(bu_vls_cstr(&found_value), value);
	bu_vls_free(&found_key);
	bu_vls_free(&found_value);
	if (matched)
	    return;
    }
    bu_exit(EXIT_FAILURE,
	"%s primitive %d metadata %s=%s was not published.\n",
	feature_name, primitive, key, value);
}

static void
require_command_result_label(struct ged_view_context *view_ctx, const char *feature_name,
	const char *text_fragment)
{
    struct ged_view_feature_summary feature_summary =
	GED_VIEW_FEATURE_SUMMARY_INIT;
    if (!ged_view_feature_get_summary(view_ctx, feature_name,
	    &feature_summary) ||
	    !feature_summary.exists ||
	    !feature_summary.is_overlay ||
	    !feature_summary.is_label ||
	    !feature_summary.is_command_result ||
	    !feature_summary.metadata_count)
	bu_exit(EXIT_FAILURE,
		"%s was not published as a command-result label.\n",
		feature_name);

    if (ged_annotation_label_count(view_ctx, feature_name) != 1)
	bu_exit(EXIT_FAILURE,
		"%s did not publish exactly one label.\n",
		feature_name);

    struct bu_vls text = BU_VLS_INIT_ZERO;
    point_t point = VINIT_ZERO;
    unsigned char rgb[3] = {0, 0, 0};
    if (!ged_annotation_label_copy(view_ctx, feature_name, 0,
	    &text, point, rgb)) {
	bu_vls_free(&text);
	bu_exit(EXIT_FAILURE,
		"%s label could not be read back.\n",
		feature_name);
    }
    if (!strstr(bu_vls_cstr(&text), text_fragment) ||
	    rgb[0] != 255 || rgb[1] != 255 || rgb[2] != 0) {
	bu_vls_free(&text);
	bu_exit(EXIT_FAILURE,
		"%s label content was not preserved.\n",
		feature_name);
    }
    bu_vls_free(&text);
}

static void
require_command_result_feature(struct ged_view_context *view_ctx, const char *feature_name)
{
    struct ged_view_feature_summary feature_summary =
	GED_VIEW_FEATURE_SUMMARY_INIT;
    if (!ged_view_feature_get_summary(view_ctx, feature_name,
	    &feature_summary) ||
	    !feature_summary.exists ||
	    !feature_summary.is_overlay ||
	    !feature_summary.is_command_result ||
	    !feature_summary.metadata_count)
	bu_exit(EXIT_FAILURE,
		"%s was not published as command-result scene content.\n",
		feature_name);
}

static size_t
count_gqa_feature_segments(struct ged_view_context *view_ctx, const char *feature_name,
	const char *plot_fname, const int rgb[3], int write_plot)
{
    struct gqa_segment_writer writer = {NULL, 0};
    struct ged_draw_view_record_query query;
    memset(&query, 0, sizeof(query));
    query.flags = GED_DRAW_VIEW_RECORD_QUERY_VIEW_OBJECTS;
    query.draw_mode = -1;
    struct gqa_record_capture record_capture = {&writer, plot_fname,
	feature_name, {255, 255, 255}, write_plot};
    if (rgb) {
	record_capture.rgb[0] = rgb[0];
	record_capture.rgb[1] = rgb[1];
	record_capture.rgb[2] = rgb[2];
    }
    ged_draw_foreach_view_record_query(view_ctx, &query, capture_gqa_record,
	    &record_capture);
    if (writer.fp)
	fclose(writer.fp);
    return writer.count;
}

static int
make_gqa_comb(struct rt_wdb *wdbp, const char *name, const char **members,
	const int *ops, size_t member_count, int region, int air, int material)
{
    struct wmember wm;
    size_t i;

    BU_LIST_INIT(&wm.l);
    for (i = 0; i < member_count; i++) {
	if (!mk_addmember(members[i], &wm.l, NULL, ops[i]))
	    return 0;
    }
    return mk_comb(wdbp, name, &wm.l, region, NULL, NULL, NULL,
	    0, air, material, 0, 0, 0, 0) == 0;
}

static struct ged *
open_gqa_result_fixture(void)
{
    char tmppath[MAXPATHLEN] = {0};
    FILE *fp = bu_temp_file(tmppath, MAXPATHLEN);
    if (!fp)
	return NULL;
    fclose(fp);

    struct rt_wdb *wdbp = wdb_fopen(tmppath);
    if (!wdbp)
	return NULL;

    point_t box1_min = {0.0, 0.0, 0.0};
    point_t box1_max = {1000.0, 1000.0, 1000.0};
    point_t box2_min = {100.0, 100.0, 100.0};
    point_t box2_max = {900.0, 900.0, 900.0};
    point_t box3_min = {100.0, 100.0, 100.0};
    point_t box3_max = {1000.0, 900.0, 900.0};
    point_t box4_min = {50.0, 50.0, 50.0};
    point_t box4_max = {950.0, 500.0, 950.0};

    int ok = mk_rpp(wdbp, "box1.s", box1_min, box1_max) == 0 &&
	mk_rpp(wdbp, "box2.s", box2_min, box2_max) == 0 &&
	mk_rpp(wdbp, "box3.s", box3_min, box3_max) == 0 &&
	mk_rpp(wdbp, "box4.s", box4_min, box4_max) == 0;

    if (ok) {
	const char *members[] = {"box1.s", "box2.s"};
	const int ops[] = {WMOP_UNION, WMOP_SUBTRACT};
	ok = make_gqa_comb(wdbp, "closed_box.r", members, ops, 2, 1, 0,
		5);
    }
    if (ok) {
	const char *members[] = {"box1.s", "box3.s"};
	const int ops[] = {WMOP_UNION, WMOP_SUBTRACT};
	ok = make_gqa_comb(wdbp, "open_box.r", members, ops, 2, 1, 0,
		5);
    }
    if (ok) {
	const char *members[] = {"box3.s"};
	const int ops[] = {WMOP_UNION};
	ok = make_gqa_comb(wdbp, "exposed_air.r", members, ops, 1, 1, 2,
		2);
    }
    if (ok) {
	const char *members[] = {"exposed_air.r", "open_box.r"};
	const int ops[] = {WMOP_UNION, WMOP_UNION};
	ok = make_gqa_comb(wdbp, "exposed_air.g", members, ops, 2, 0, 0,
		0);
    }
    if (ok) {
	const char *members[] = {"box2.s", "box4.s"};
	const int ops[] = {WMOP_UNION, WMOP_UNION};
	ok = make_gqa_comb(wdbp, "adj_air1.r", members, ops, 2, 1, 3, 2);
    }
    if (ok) {
	const char *members[] = {"box2.s", "box4.s"};
	const int ops[] = {WMOP_UNION, WMOP_SUBTRACT};
	ok = make_gqa_comb(wdbp, "adj_air2.r", members, ops, 2, 1, 4, 2);
    }
    if (ok) {
	const char *members[] = {"closed_box.r", "adj_air1.r",
	    "adj_air2.r"};
	const int ops[] = {WMOP_UNION, WMOP_UNION, WMOP_UNION};
	ok = make_gqa_comb(wdbp, "adj_air.g", members, ops, 3, 0, 0, 0);
    }
    if (ok) {
	const char *members[] = {"closed_box.r", "adj_air2.r"};
	const int ops[] = {WMOP_UNION, WMOP_UNION};
	ok = make_gqa_comb(wdbp, "gap.g", members, ops, 2, 0, 0, 0);
    }

    wdb_close(wdbp);
    if (!ok) {
	bu_file_delete(tmppath);
	return NULL;
    }

    return ged_open("db", tmppath, 1);
}

static void
run_gqa_result_case(struct ged *gedp, const char *analysis_flag,
	const char *object_name, const char *feature_name,
	const char *result_kind, const char *schema, const char *severity)
{
    char grid_arg[] = "250mm-50mm";
    const char *gqa[] = {
	"gqa", "-P", "1", analysis_flag, "-g", grid_arg,
	object_name, NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    if (ged_exec_gqa(gedp, 7, gqa) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "GQA %s failed: %s\n", feature_name,
		bu_vls_cstr(gedp->ged_result_str));

    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx)
	bu_exit(EXIT_FAILURE, "No active GED view context available.\n");
    require_command_result_feature(view_ctx, feature_name);
    require_command_result_metadata(view_ctx, feature_name, "result.owner",
	    "gqa");
    require_command_result_metadata(view_ctx, feature_name, "result.kind",
	    result_kind);
    require_command_result_metadata(view_ctx, feature_name, "result.schema",
	    schema);
    require_command_result_metadata(view_ctx, feature_name,
	    "result.severity", severity);
    if (!count_gqa_feature_segments(view_ctx, feature_name, NULL, NULL, 0))
	bu_exit(EXIT_FAILURE, "No GQA segments exported for %s.\n",
		feature_name);
    require_command_result_primitive_metadata(view_ctx, feature_name, 0,
	"result.schema", schema);
    require_command_result_primitive_metadata(view_ctx, feature_name, 0,
	"result.severity", severity);
}

static void
run_gqa_volume_result_case(struct ged *gedp)
{
    char grid_arg[] = "250mm-50mm";
    char prefix[MAXPATHLEN] = {0};
    FILE *fp = bu_temp_file(prefix, MAXPATHLEN);
    if (!fp)
	bu_exit(EXIT_FAILURE, "Could not create GQA volume plot prefix.\n");
    fclose(fp);
    bu_file_delete(prefix);

    const char *gqa[] = {
	"gqa", "-P", "1", "-Av", "-g", grid_arg, "-p", prefix,
	"closed_box.r", NULL
    };
    bu_vls_trunc(gedp->ged_result_str, 0);
    if (ged_exec_gqa(gedp, 9, gqa) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "GQA volume samples failed: %s\n",
		bu_vls_cstr(gedp->ged_result_str));

    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx)
	bu_exit(EXIT_FAILURE, "No active GED view context available.\n");

    require_command_result_feature(view_ctx, "gqa::volume-samples");
    require_command_result_metadata(view_ctx, "gqa::volume-samples",
	    "result.owner", "gqa");
    require_command_result_metadata(view_ctx, "gqa::volume-samples",
	    "result.kind", "volume-sample");
    require_command_result_metadata(view_ctx, "gqa::volume-samples",
	    "result.schema", "brlcad.gqa.volume-sample.v1");
    require_command_result_metadata(view_ctx, "gqa::volume-samples",
	    "result.severity", "info");
    if (!count_gqa_feature_segments(view_ctx, "gqa::volume-samples", NULL,
	    NULL, 0))
	bu_exit(EXIT_FAILURE,
		"No GQA volume sample segments exported.\n");
    require_command_result_primitive_metadata(view_ctx,
	"gqa::volume-samples", 0, "result.severity", "info");

    struct bu_vls plot_file = BU_VLS_INIT_ZERO;
    bu_vls_printf(&plot_file, "%svolume.plot3", prefix);
    bu_file_delete(bu_vls_cstr(&plot_file));
    bu_vls_free(&plot_file);
}

int
main(int ac, char *av[]) {
    struct ged *gedp;
    const char *gqa_plot_fname = "gqa_ovlps.plot3";
    const char *gqa[4] = {"gqa", "-Aop", "ovlp", NULL};
    char check_grid_arg[] = "250mm-50mm";
    const char *check[8] = {
	"check", "overlaps", "-P", "1", "-o", "-g", check_grid_arg, "ovlp"
    };
    const int overlap_rgb[3] = {255, 255, 0};

    bu_setprogname(av[0]);

    if (ac != 2) {
	printf("Usage: %s file.g\n", av[0]);
	return 1;
    }
    if (!bu_file_exists(av[1], NULL)) {
	printf("ERROR: [%s] does not exist, expecting .g file\n", av[1]);
	return 2;
    }

    gedp = ged_open("db", av[1], 1);
    struct ged_draw_transaction init_txn =
	ged_draw_transaction_make(GED_DRAW_TXN_NONE, NULL);
    if (ged_draw_apply_transaction(gedp, &init_txn, NULL) < 0)
	bu_exit(EXIT_FAILURE, "Could not initialize draw state for GQA test.\n");
    ged_exec_gqa(gedp, 3, gqa);
    printf("%s\n", bu_vls_cstr(gedp->ged_result_str));

    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx)
	bu_exit(EXIT_FAILURE, "No active GED view context available.\n");

    require_command_result_feature(view_ctx, "gqa::overlaps");
    require_command_result_metadata(view_ctx, "gqa::overlaps",
	    "result.owner", "gqa");
    require_command_result_metadata(view_ctx, "gqa::overlaps",
	    "result.kind", "overlap");
    require_command_result_metadata(view_ctx, "gqa::overlaps",
	    "result.schema", "brlcad.gqa.overlap.v1");
    require_command_result_metadata(view_ctx, "gqa::overlaps",
	    "result.severity", "error");
    require_command_result_label(view_ctx, "gqa::overlaps::summary",
	    "gqa overlaps:");
    require_command_result_metadata(view_ctx, "gqa::overlaps::summary",
	    "result.kind", "overlap-summary");
    size_t overlap_count = count_gqa_feature_segments(view_ctx,
	    "gqa::overlaps", gqa_plot_fname, overlap_rgb, 1);

    if (overlap_count) {
	printf("Writing plot data to %s for inspection with overlay command\n", gqa_plot_fname);
    } else {
	bu_exit(EXIT_FAILURE, "No GQA plotting data found.\n");
    }
    require_command_result_primitive_metadata(view_ctx, "gqa::overlaps", 0,
	"result.severity", "error");

    bu_vls_trunc(gedp->ged_result_str, 0);
    if (ged_exec(gedp, 8, check) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "check overlaps failed: %s\n",
		bu_vls_cstr(gedp->ged_result_str));
    require_command_result_feature(view_ctx, "check::overlaps");
    require_command_result_metadata(view_ctx, "check::overlaps",
	    "result.owner", "check");
    require_command_result_metadata(view_ctx, "check::overlaps",
	    "result.kind", "overlap");
    require_command_result_metadata(view_ctx, "check::overlaps",
	    "result.schema", "brlcad.check.overlap.v1");
    require_command_result_metadata(view_ctx, "check::overlaps",
	    "result.severity", "error");
    require_command_result_label(view_ctx, "check::overlaps::summary",
	    "check overlaps:");
    require_command_result_metadata(view_ctx, "check::overlaps::summary",
	    "result.kind", "overlap-summary");
    if (!count_gqa_feature_segments(view_ctx, "check::overlaps", NULL,
	    NULL, 0))
	bu_exit(EXIT_FAILURE, "No check overlap segments exported.\n");
    require_command_result_primitive_metadata(view_ctx, "check::overlaps", 0,
	"result.schema", "brlcad.check.overlap.v1");
    require_command_result_primitive_metadata(view_ctx, "check::overlaps", 0,
	"result.severity", "error");

    ged_close(gedp);

    gedp = open_gqa_result_fixture();
    if (!gedp)
	bu_exit(EXIT_FAILURE, "Could not create GQA result fixture database.\n");
    init_txn = ged_draw_transaction_make(GED_DRAW_TXN_NONE, NULL);
    if (ged_draw_apply_transaction(gedp, &init_txn, NULL) < 0)
	bu_exit(EXIT_FAILURE, "Could not initialize draw state for the GQA result fixture.\n");

    run_gqa_result_case(gedp, "-Ag", "gap.g", "gqa::gaps", "gap",
	"brlcad.gqa.gap.v1", "warning");
    run_gqa_result_case(gedp, "-Aa", "adj_air.g", "gqa::adjacent-air",
	    "adjacent-air", "brlcad.gqa.adjacent-air.v1", "warning");
    run_gqa_result_case(gedp, "-Ae", "exposed_air.g", "gqa::exposed-air",
	    "exposed-air", "brlcad.gqa.exposed-air.v1", "warning");
    run_gqa_volume_result_case(gedp);

    ged_close(gedp);

    return 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
