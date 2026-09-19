/*                    L I N E _ L A Y E R . C
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

#include "common.h"

#include "bu/str.h"

#include <stdio.h>
#include <string.h>

#include "vmath.h"
#include "bg/line_layer.h"
#include "bg/plot3.h"
#include "bu/app.h"
#include "bu/log.h"


static int
expect_point(const char *name, const point_t actual, fastf_t x, fastf_t y, fastf_t z)
{
    if (!NEAR_EQUAL(actual[X], x, SMALL_FASTF) ||
	!NEAR_EQUAL(actual[Y], y, SMALL_FASTF) ||
	!NEAR_EQUAL(actual[Z], z, SMALL_FASTF)) {
	bu_log("FAIL: %s is (%g, %g, %g), expected (%g, %g, %g)\n",
	       name, V3ARGS(actual), x, y, z);
	return 1;
    }

    return 0;
}


static int
expect_plot3_text(FILE *fp, const char *label, const char *expected)
{
    char got[256] = {0};
    size_t nread = 0;

    if (fflush(fp) != 0) {
	bu_log("FAIL: fflush failed for %s plot3 stream\n", label);
	return 1;
    }
    rewind(fp);

    nread = fread(got, 1, sizeof(got) - 1, fp);
    if (ferror(fp)) {
	bu_log("FAIL: fread failed for %s plot3 stream\n", label);
	return 1;
    }
    got[nread] = '\0';

    if (bu_strcmp(got, expected) != 0) {
	bu_log("FAIL: %s plot3 stream is:\n%s\nexpected:\n%s\n",
	       label, got, expected);
	return 1;
    }

    return 0;
}


static int
expect_plot3_export(const struct bg_line_layer_builder *builder,
		    const struct bg_line_layer *layer)
{
    const char *layer_expected = "C 255 0 0\nO 0 1 2\nQ 3 4 5\n";
    const char *builder_expected =
	"C 255 0 0\n"
	"O 0 1 2\n"
	"Q 3 4 5\n"
	"C 0 0 255\n"
	"X 6 7 8\n";
    FILE *fp = NULL;
    int old_mode = pl_getOutputMode();
    int ret = 0;

    fp = tmpfile();
    if (!fp) {
	bu_log("FAIL: tmpfile failed for line-layer plot3 export test\n");
	return 1;
    }

    pl_setOutputMode(PL_OUTPUT_MODE_TEXT);
    if (!bg_line_layer_write_plot3(fp, layer)) {
	bu_log("FAIL: bg_line_layer_write_plot3 rejected a valid layer\n");
	ret = 1;
    } else {
	ret |= expect_plot3_text(fp, "single-layer", layer_expected);
    }
    fclose(fp);

    fp = tmpfile();
    if (!fp) {
	pl_setOutputMode(old_mode);
	bu_log("FAIL: tmpfile failed for line-layer builder plot3 export test\n");
	return 1;
    }

    if (!bg_line_layer_builder_write_plot3(fp, builder)) {
	bu_log("FAIL: bg_line_layer_builder_write_plot3 rejected a valid builder\n");
	ret = 1;
    } else {
	ret |= expect_plot3_text(fp, "builder", builder_expected);
    }
    fclose(fp);
    pl_setOutputMode(old_mode);

    if (bg_line_layer_write_plot3(NULL, layer) ||
	bg_line_layer_write_plot3(stdout, NULL) ||
	bg_line_layer_builder_write_plot3(NULL, builder) ||
	bg_line_layer_builder_write_plot3(stdout, NULL)) {
	bu_log("FAIL: line-layer plot3 export accepted invalid inputs\n");
	ret = 1;
    }

    return ret;
}


int
main(int UNUSED(argc), char *argv[])
{
    struct bg_line_layer_builder *builder = NULL;
    struct bg_line_layer *red = NULL;
    const struct bg_line_layer *first = NULL;
    point_t p0 = {0.0, 1.0, 2.0};
    point_t p1 = {3.0, 4.0, 5.0};
    point_t p2 = {6.0, 7.0, 8.0};
    point_t got = VINIT_ZERO;
    unsigned char r = 0;
    unsigned char g = 0;
    unsigned char b = 0;
    int cmd = -1;
    int ret = 0;

    bu_setprogname(argv[0]);

    builder = bg_line_layer_builder_create();
    if (!builder) {
	bu_log("FAIL: bg_line_layer_builder_create returned NULL\n");
	return 1;
    }

    if (bg_line_layer_builder_layer_count(builder) != 0 ||
	bg_line_layer_builder_point_count(builder) != 0) {
	bu_log("FAIL: new line-layer builder is not empty\n");
	ret = 1;
    }

    if (!bg_line_layer_builder_add(builder, 255, 0, 0, p0, BG_GEOMETRY_LINE_MOVE) ||
	!bg_line_layer_builder_add(builder, 255, 0, 0, p1, BG_GEOMETRY_LINE_DRAW) ||
	!bg_line_layer_builder_add(builder, 0, 0, 255, p2, BG_GEOMETRY_POINT_DRAW)) {
	bu_log("FAIL: bg_line_layer_builder_add rejected valid commands\n");
	ret = 1;
    }

    if (bg_line_layer_builder_layer_count(builder) != 2 ||
	bg_line_layer_builder_point_count(builder) != 3) {
	bu_log("FAIL: unexpected line-layer builder counts\n");
	ret = 1;
    }

    red = bg_line_layer_builder_find(builder, 255, 0, 0);
    if (!red) {
	bu_log("FAIL: bg_line_layer_builder_find returned NULL for an existing layer\n");
	ret = 1;
    }

    if (red && (!bg_line_layer_color(red, &r, &g, &b) ||
	    r != 255 || g != 0 || b != 0)) {
	bu_log("FAIL: red layer color is (%u, %u, %u)\n", r, g, b);
	ret = 1;
    }

    if (red && (bg_line_layer_point_count(red) != 2 ||
	    !bg_line_layer_point_at(red, 1, got) ||
	    bg_line_layer_command_at(red, 1, &cmd) == 0 ||
	    cmd != BG_GEOMETRY_LINE_DRAW)) {
	bu_log("FAIL: red layer did not preserve second point/command\n");
	ret = 1;
    }
    ret |= expect_point("red[1]", got, 3.0, 4.0, 5.0);

    if (red && bg_line_layer_add(red, p0, -1)) {
	bu_log("FAIL: bg_line_layer_add accepted an invalid command\n");
	ret = 1;
    }

    first = bg_line_layer_builder_layer_at(builder, 0);
    if (!first || !bg_line_layer_name(first) ||
	!bg_line_layer_points(first) ||
	!bg_line_layer_commands(first)) {
	bu_log("FAIL: line-layer direct accessors did not return populated data\n");
	ret = 1;
    }

    if (bg_line_layer_builder_layer_at(builder, 2)) {
	bu_log("FAIL: bg_line_layer_builder_layer_at accepted an out-of-range index\n");
	ret = 1;
    }

    if (red)
	ret |= expect_plot3_export(builder, red);

    bg_line_layer_builder_clear(builder);
    if (bg_line_layer_builder_layer_count(builder) != 0 ||
	bg_line_layer_builder_point_count(builder) != 0) {
	bu_log("FAIL: bg_line_layer_builder_clear did not reset counts\n");
	ret = 1;
    }

    bg_line_layer_builder_free(builder);
    return ret;
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
