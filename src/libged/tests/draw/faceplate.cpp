/*                    F A C E P L A T E . C P P
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
/** @file faceplate.cpp
 *
 * Testing routines for drawing "built-in" faceplate view elements
 *
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include <stdio.h>
#include <fstream>

#include <bu.h>
#include <BObol/BDisplayEndpoint.h>
#include <icv.h>
#include <imgstream/fbserv.h>
#include "rt/view.h"
#include "view_test_util.h"
#include <ged.h>
#include <ged/draw.h>
#include <ged/display.h>

#define ADIFF_THRES 0.99
/* The retained wire renderer preserves projected topology but does not emit
 * the legacy display manager's identical primitive tessellation. */
#define WIREFRAME_ADIFF_THRES 0.98

extern "C" int img_cmp(int id, struct ged *gedp, const char *cdir, bool clear_scene, bool clear_image, int soft_fail, fastf_t approximate_check, const char *clear_root, const char *img_root);
extern "C" int img_not_empty(int id, struct ged *gedp, const char *cdir, bool clear_scene, bool clear_image, int soft_fail, const char *clear_root, const char *img_root);
extern "C" int unpack_apng(const char *src_dir, const char *apng_name, const char *out_dir, const char *prefix);

static bool
framebuffer_matches(bobol_display_endpoint_t *endpoint,
	const unsigned char *expected, size_t expected_size,
	unsigned int expected_width, unsigned int expected_height)
{
    unsigned char *pixels = NULL;
    size_t size = 0;
    unsigned int width = 0;
    unsigned int height = 0;
    unsigned int components = 0;
    bool matched = endpoint &&
	bobol_display_endpoint_capture_plane(endpoint,
	    BOBOL_CAPTURE_FRAMEBUFFER, &pixels, &size, &width, &height,
	    &components) && pixels && size == expected_size &&
	width == expected_width && height == expected_height &&
	components == 3 && memcmp(pixels, expected, expected_size) == 0;
    if (!matched)
	bu_log("framebuffer comparison: captured=%ux%u components=%u bytes=%zu, expected=%ux%u components=3 bytes=%zu, pixels=%p\n",
	    width, height, components, size, expected_width, expected_height,
	    expected_size, (void *)pixels);
    if (!matched && pixels && size == expected_size &&
	width == expected_width && height == expected_height && components == 3) {
	size_t mismatched = 0;
	size_t flipped_mismatched = 0;
	const size_t row_size = (size_t)width * 3;
	for (unsigned int y = 0; y < height; y++) {
	    const unsigned char *actual_row = pixels + (size_t)y * row_size;
	    const unsigned char *expected_row = expected + (size_t)y * row_size;
	    const unsigned char *flipped_row =
		expected + (size_t)(height - 1 - y) * row_size;
	    for (size_t x = 0; x < row_size; x++) {
		mismatched += actual_row[x] != expected_row[x];
		flipped_mismatched += actual_row[x] != flipped_row[x];
	    }
	}
	bu_log("framebuffer comparison: mismatched=%zu, vertically-flipped mismatched=%zu, bytes=%zu\n",
	    mismatched, flipped_mismatched, expected_size);
    }
    if (pixels)
	bu_free(pixels, "faceplate framebuffer comparison");
    return matched;
}

static void
wait_for_progressive_draw(struct ged *gedp, struct ged_view_context *view_ctx)
{
    if (!draw_test_obol_progressive_drain(gedp, view_ctx, 2000, 1))
	bu_exit(EXIT_FAILURE,
	    "Obol progressive realization did not settle before baseline capture\n");
}

static struct ged_view_feature_batch *
faceplate_feature_batch(struct ged_view_context *view_ctx)
{
    struct ged_view_feature_batch_desc desc = ged_view_feature_batch_desc_default();
    desc.owner_id = "faceplate-test";
    desc.owner_role = "test";
    desc.generation = 1;
    return ged_view_feature_batch_begin(view_ctx, &desc);
}

static int
faceplate_hud_axes_replace(struct ged_view_context *view_ctx, const char *name,
	const struct bv_axes_state *axes, const mat_t rotation)
{
    struct ged_view_feature_batch *batch = faceplate_feature_batch(view_ctx);
    if (!batch)
	return 0;
    if (!ged_view_feature_batch_hud_axes_replace(batch, name, axes, rotation)) {
	ged_view_feature_batch_abort(batch);
	return 0;
    }
    return ged_view_feature_batch_commit(batch);
}

static int
faceplate_hud_lines_replace(struct ged_view_context *view_ctx, const char *name,
	const point_t *points, const int *cmds, size_t point_count,
	const struct ged_view_feature_style *style)
{
    struct ged_view_feature_batch *batch = faceplate_feature_batch(view_ctx);
    if (!batch)
	return 0;
    if (!ged_view_feature_batch_hud_line_set_replace(batch, name, points,
	    cmds, point_count, style)) {
	ged_view_feature_batch_abort(batch);
	return 0;
    }
    return ged_view_feature_batch_commit(batch);
}

static int
faceplate_hud_labels_replace(struct ged_view_context *view_ctx, const char *name,
	const struct ged_view_feature_label *labels, size_t label_count,
	const struct ged_view_feature_style *style)
{
    struct ged_view_feature_batch *batch = faceplate_feature_batch(view_ctx);
    if (!batch)
	return 0;
    if (!ged_view_feature_batch_hud_labels_replace(batch, name, labels,
	    label_count, style)) {
	ged_view_feature_batch_abort(batch);
	return 0;
    }
    return ged_view_feature_batch_commit(batch);
}

static int
faceplate_hud_line_layers_replace(struct ged_view_context *view_ctx,
	const char *name, const struct ged_view_feature_line_layer *layers,
	size_t layer_count, const struct ged_view_feature_style *style)
{
    struct ged_view_feature_batch *batch = faceplate_feature_batch(view_ctx);
    if (!batch)
	return 0;
    if (!ged_view_feature_batch_hud_line_layers_replace(batch, name, layers,
	    layer_count, style)) {
	ged_view_feature_batch_abort(batch);
	return 0;
    }
    return ged_view_feature_batch_commit(batch);
}

int
main(int ac, char *av[]) {
    struct ged *gedp;
    struct bu_vls fname = BU_VLS_INIT_ZERO;
    int need_help = 0;
    int run_unstable_tests = 0;
    int soft_fail = 0;
    int keep_images = 0;
    int ret = BRLCAD_OK;

    bu_setprogname(av[0]);

    struct bu_opt_desc d[5];
    BU_OPT(d[0], "h", "help",            "", NULL,          &need_help, "Print help and exit");
    BU_OPT(d[1], "U", "enable-unstable", "", NULL, &run_unstable_tests, "Test drawing routines known to differ between build configs/platforms.");
    BU_OPT(d[2], "c", "continue",        "", NULL,          &soft_fail, "Continue testing if a failure is encountered.");
    BU_OPT(d[3], "k", "keep",            "", NULL,        &keep_images, "Keep images generated by the run.");
    BU_OPT_NULL(d[4]);

    /* Done with program name */
    int uac = bu_opt_parse(NULL, ac, (const char **)av, d);
    if (uac == -1 || need_help)

    if (ac != 2)
	bu_exit(EXIT_FAILURE, "%s [-h] [-U] <directory>", av[0]);

    if (!bu_file_directory(av[1])) {
	printf("ERROR: [%s] is not a directory.  Expecting control image directory\n", av[1]);
	return 2;
    }

    bool clear_images = !keep_images;

    /* Enable all the experimental logic */
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    /* Use a local working-directory cache so we do not pollute the user's
     * real BRL-CAD cache and so the test is fully self-contained. */
    char lcache[MAXPATHLEN] = {0};
    char runtime_cache[MAXPATHLEN] = {0};
    bu_dir(lcache, MAXPATHLEN, BU_DIR_CURR, "ged_draw_test_fp_cache", NULL);
    bu_mkdir(lcache);
    bu_dir(runtime_cache, MAXPATHLEN, BU_DIR_CURR, "ged_draw_test_fp_cache",
	   "cache", NULL);
    bu_mkdir(runtime_cache);
    /* Cache maintenance must not erase extracted image controls. */
    bu_setenv("BU_DIR_CACHE", runtime_cache, 1);

    if (!bu_file_exists(av[1], NULL)) {
	printf("ERROR: [%s] does not exist, expecting .g file\n", av[1]);
	return 2;
    }

    unpack_apng(av[1], "faceplate.apng", lcache, "fp");

    /* Open the temp file, then dbconcat argv[1] into it */
    bu_vls_sprintf(&fname, "%s/moss.g", av[1]);
    gedp = ged_open("db", bu_vls_cstr(&fname), 1);

    /* Image baselines use the GED-owned headless Obol render endpoint. */
    const char *s_av[15] = {NULL};
    struct ged_view_context *v = ged_view_active_ctx(gedp);
    if (draw_test_obol_view_init(gedp, v, 512, 512) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "failed to initialize headless Obol render endpoint\n");

    /***** Sanity - basic wireframe draw *****/
    bu_log("Testing basic db wireframe draw...\n");
    s_av[0] = "draw";
    s_av[1] = "all.g";
    s_av[2] = NULL;
    ged_exec_draw(gedp, 2, s_av);

    wait_for_progressive_draw(gedp, v);

    s_av[0] = "autoview";
    s_av[1] = NULL;
    ged_exec_autoview(gedp, 1, s_av);

    s_av[0] = "ae";
    s_av[1] = "35";
    s_av[2] = "25";
    s_av[3] = NULL;
    ged_exec_ae(gedp, 3, s_av);
    ret += img_cmp(1, gedp, lcache, true, clear_images, soft_fail,
	WIREFRAME_ADIFF_THRES, "faceplate_clear", "fp");

    // Check that everything is in fact cleared
    ret += img_cmp(0, gedp, lcache, false, clear_images, soft_fail, 0, "faceplate_clear", "fp");
    bu_log("Done.\n");

    /***** Center Dot *****/
    bu_log("Testing center dot...\n");
    s_av[0] = "view";
    s_av[1] = "faceplate";
    s_av[2] = "center_dot";
    s_av[3] = "1";
    s_av[4] = NULL;
    ged_exec_view(gedp, 4, s_av);
    ret += img_cmp(2, gedp, lcache, false, clear_images, soft_fail, ADIFF_THRES, "faceplate_clear", "fp");

    // Check that turning off works
    s_av[3] = "0";
    ged_exec_view(gedp, 4, s_av);
    ret += img_cmp(0, gedp, lcache, false, clear_images, soft_fail, 0, "faceplate_clear", "fp");
    bu_log("Done.\n");

    /***** Grid *****/
    bu_log("Testing grid...\n");

    s_av[0] = "view";
    s_av[1] = "faceplate";
    s_av[2] = "grid";
    s_av[3] = "1";
    s_av[4] = NULL;
    ged_exec_view(gedp, 4, s_av);
    ret += img_not_empty(3, gedp, lcache, false, clear_images, soft_fail, "faceplate_clear", "fp");

    // Exercise the draw subcommand too: it stages grid state before the
    // endpoint provider applies the visibility property.
    s_av[3] = "draw";
    s_av[4] = "0";
    s_av[5] = NULL;
    ged_exec_view(gedp, 5, s_av);
    ret += img_cmp(0, gedp, lcache, false, clear_images, soft_fail, 0, "faceplate_clear", "fp");
    bu_log("Done.\n");

    /***** Params *****/
    bu_log("Testing parameters reporting...\n");

    // First make the font smaller so we can see all params
    s_av[0] = "view";
    s_av[1] = "faceplate";
    s_av[2] = "params";
    s_av[3] = "font_size";
    s_av[4] = "10";
    s_av[5] = NULL;
    ged_exec_view(gedp, 5, s_av);

    s_av[0] = "view";
    s_av[1] = "faceplate";
    s_av[2] = "params";
    s_av[3] = "1";
    s_av[4] = NULL;
    ged_exec_view(gedp, 4, s_av);
    ret += img_not_empty(4, gedp, lcache, false, clear_images, soft_fail, "faceplate_clear", "fp");

    bu_log("Testing turning on frames per second reporting...\n");

    // So we don't get random values here, override the timing variable values
    bv_frametime_set(DRAW_TEST_BV(v), 1000000000);

    s_av[0] = "view";
    s_av[1] = "faceplate";
    s_av[2] = "params";
    s_av[3] = "fps";
    s_av[4] = "1";
    s_av[5] = NULL;
    ged_exec_view(gedp, 5, s_av);
    ret += img_not_empty(5, gedp, lcache, false, clear_images, soft_fail, "faceplate_clear", "fp");

    // Check that turning off works
    s_av[3] = "0";
    ged_exec_view(gedp, 4, s_av);
    ret += img_cmp(0, gedp, lcache, false, clear_images, soft_fail, 0, "faceplate_clear", "fp");
    bu_log("Done.\n");

    // Restore default font size
    s_av[0] = "view";
    s_av[1] = "faceplate";
    s_av[2] = "params";
    s_av[3] = "font_size";
    s_av[4] = "20";
    s_av[5] = NULL;
    ged_exec_view(gedp, 5, s_av);



    /***** Scale *****/
    bu_log("Testing scale reporting...\n");
    s_av[0] = "view";
    s_av[1] = "faceplate";
    s_av[2] = "scale";
    s_av[3] = "1";
    s_av[4] = NULL;
    ged_exec_view(gedp, 4, s_av);
    ret += img_not_empty(6, gedp, lcache, false, clear_images, soft_fail, "faceplate_clear", "fp");

    // Check that turning off works
    s_av[3] = "0";
    ged_exec_view(gedp, 4, s_av);
    ret += img_cmp(0, gedp, lcache, false, clear_images, soft_fail, 0, "faceplate_clear", "fp");
    bu_log("Done.\n");


    /***** View axes *****/
    bu_log("Testing view axes drawing...\n");
    s_av[0] = "view";
    s_av[1] = "faceplate";
    s_av[2] = "view_axes";
    s_av[3] = "1";
    s_av[4] = NULL;
    ged_exec_view(gedp, 4, s_av);
    ret += img_not_empty(7, gedp, lcache, false, clear_images, soft_fail, "faceplate_clear", "fp");

    // Check that turning off works
    s_av[3] = "0";
    ged_exec_view(gedp, 4, s_av);
    ret += img_cmp(0, gedp, lcache, false, clear_images, soft_fail, 0, "faceplate_clear", "fp");
    bu_log("Done.\n");

    /***** Model axes *****/
    bu_log("Testing model axes drawing...\n");
    s_av[0] = "view";
    s_av[1] = "faceplate";
    s_av[2] = "model_axes";
    s_av[3] = "1";
    s_av[4] = NULL;
    ged_exec_view(gedp, 4, s_av);
    ret += img_not_empty(8, gedp, lcache, false, clear_images, soft_fail, "faceplate_clear", "fp");

    // Check that turning off works
    s_av[3] = "0";
    ged_exec_view(gedp, 4, s_av);
    ret += img_cmp(0, gedp, lcache, false, clear_images, soft_fail, 0, "faceplate_clear", "fp");
    bu_log("Done.\n");

    /***** Application-owned retained HUD axes *****/
    struct bv_axes_state edit_axes = BV_AXES_STATE_INIT;
    mat_t edit_rotation;
    MAT_IDN(edit_rotation);
    edit_axes.draw = 1;
    edit_axes.axes_size = 0.25;
    edit_axes.line_width = 2;
    edit_axes.label_flag = 1;
    VSET(edit_axes.axes_color, 255, 128, 0);
    VSET(edit_axes.label_color, 255, 255, 255);
    if (!faceplate_hud_axes_replace(v, "_test/edit_axes",
	    &edit_axes, edit_rotation) ||
	!ged_view_feature_exists(v, "_test/edit_axes/lines") ||
	!ged_view_feature_exists(v, "_test/edit_axes/labels"))
	bu_exit(EXIT_FAILURE, "retained HUD axes publication failed\n");
    edit_axes.draw = 0;
    if (!faceplate_hud_axes_replace(v, "_test/edit_axes",
	    &edit_axes, edit_rotation) ||
	ged_view_feature_exists(v, "_test/edit_axes/lines") ||
	ged_view_feature_exists(v, "_test/edit_axes/labels"))
	bu_exit(EXIT_FAILURE, "retained HUD axes removal failed\n");

    point_t hud_points[2] = {{-0.25, -0.25, 0.0}, {0.25, 0.25, 0.0}};
    int hud_cmds[2] = {GED_DRAW_VIEW_LINE_MOVE, GED_DRAW_VIEW_LINE_DRAW};
    struct ged_view_feature_style hud_style =
	ged_view_feature_style_default();
    hud_style.visible = 1;
    hud_style.color_valid = 1;
    VSET(hud_style.color, 32, 192, 255);
    if (!faceplate_hud_lines_replace(v, "_test/hud_lines",
	    (const point_t *)hud_points, hud_cmds, 2, &hud_style) ||
	!ged_view_feature_exists(v, "_test/hud_lines"))
	bu_exit(EXIT_FAILURE, "retained HUD line publication failed\n");
    if (!faceplate_hud_lines_replace(v, "_test/hud_lines",
	    NULL, NULL, 0, NULL) ||
	ged_view_feature_exists(v, "_test/hud_lines"))
	bu_exit(EXIT_FAILURE, "retained HUD line removal failed\n");

    struct ged_view_feature_label hud_label =
	{};
    hud_label.text = "retained HUD";
    VSET(hud_label.point, -0.5, 0.5, 0.0);
    hud_label.color_valid = 1;
    VSET(hud_label.color, 255, 196, 32);
    hud_label.font_size = 14.0;
    if (!faceplate_hud_labels_replace(v, "_test/hud_labels",
	    &hud_label, 1, &hud_style) ||
	!ged_view_feature_exists(v, "_test/hud_labels"))
	bu_exit(EXIT_FAILURE, "retained HUD label publication failed\n");
    if (!faceplate_hud_labels_replace(v, "_test/hud_labels",
	    NULL, 0, NULL) ||
	ged_view_feature_exists(v, "_test/hud_labels"))
	bu_exit(EXIT_FAILURE, "retained HUD label removal failed\n");

    struct ged_view_feature_line_layer hud_layers[2] = {
	ged_view_feature_line_layer_default(),
	ged_view_feature_line_layer_default()
    };
    hud_layers[0].name = "first";
    hud_layers[0].points = (const point_t *)hud_points;
    hud_layers[0].commands = hud_cmds;
    hud_layers[0].point_count = 2;
    hud_layers[0].style = hud_style;
    hud_layers[1] = hud_layers[0];
    hud_layers[1].name = "second";
    hud_layers[1].style.color[0] = 255;
    if (!faceplate_hud_line_layers_replace(v,
	    "_test/hud_line_layers", hud_layers, 2, &hud_style) ||
	!ged_view_feature_exists(v, "_test/hud_line_layers"))
	bu_exit(EXIT_FAILURE, "retained HUD line-layer publication failed\n");
    if (!faceplate_hud_line_layers_replace(v,
	    "_test/hud_line_layers", NULL, 0, NULL) ||
	ged_view_feature_exists(v, "_test/hud_line_layers"))
	bu_exit(EXIT_FAILURE, "retained HUD line-layer removal failed\n");

    /***** Interactive rectangle *****/
    bu_log("Testing interactive rectangle drawing...\n");
    struct bv_interactive_rect_state interactive_rect =
	BV_INTERACTIVE_RECT_STATE_INIT;
    interactive_rect.x = -0.50;
    interactive_rect.y = -0.25;
    interactive_rect.width = 1.00;
    interactive_rect.height = 0.50;
    interactive_rect.line_width = 1;
    interactive_rect.line_style = 0;
    VSET(interactive_rect.color, 32, 192, 255);
    if (!bv_interactive_rect_state_set(DRAW_TEST_BV(v), &interactive_rect))
	bu_exit(EXIT_FAILURE, "failed to configure interactive rectangle\n");

    struct bv_display_property_value interactive_value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    interactive_value.type = BV_DISPLAY_PROPERTY_UINT;
    interactive_value.uint_value = 3;
    if (ged_view_context_display_property_set(v,
	    "view.interactive.rectangle.line_width", &interactive_value) !=
	BV_DISPLAY_PROPERTY_OK)
	bu_exit(EXIT_FAILURE, "failed to set interactive rectangle line width\n");
    interactive_value.uint_value = 1;
    if (ged_view_context_display_property_set(v,
	    "view.interactive.rectangle.line_style", &interactive_value) !=
	BV_DISPLAY_PROPERTY_OK)
	bu_exit(EXIT_FAILURE, "failed to set interactive rectangle line style\n");
    struct bv_display_property_value interactive_readback =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    if (ged_view_context_display_property_get(v,
	    "view.interactive.rectangle.line_style", &interactive_readback) !=
	BV_DISPLAY_PROPERTY_OK || interactive_readback.uint_value != 1)
	bu_exit(EXIT_FAILURE, "interactive rectangle line-style readback failed\n");
    interactive_value.uint_value = 2;
    if (ged_view_context_display_property_set(v,
	    "view.interactive.rectangle.line_style", &interactive_value) !=
	BV_DISPLAY_PROPERTY_INVALID)
	bu_exit(EXIT_FAILURE, "invalid interactive rectangle line style was accepted\n");
    interactive_value.type = BV_DISPLAY_PROPERTY_COLOR3;
    interactive_value.color3[0] = 1.0;
    interactive_value.color3[1] = 0.25;
    interactive_value.color3[2] = 0.0;
    if (ged_view_context_display_property_set(v,
	    "view.interactive.rectangle.color", &interactive_value) !=
	BV_DISPLAY_PROPERTY_OK)
	bu_exit(EXIT_FAILURE, "failed to set interactive rectangle color\n");
    interactive_value.type = BV_DISPLAY_PROPERTY_BOOL;
    interactive_value.bool_value = 1;
    if (ged_view_context_display_property_set(v,
	    "view.interactive.rectangle.visible", &interactive_value) !=
	BV_DISPLAY_PROPERTY_OK)
	bu_exit(EXIT_FAILURE, "failed to show interactive rectangle\n");
    if (ged_view_faceplate_sync(gedp, v) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "failed to synchronize interactive rectangle\n");
    if (!ged_view_feature_exists(v,
	    "_faceplate/interactive_rect"))
	bu_exit(EXIT_FAILURE, "interactive rectangle was not retained\n");
    ret += img_not_empty(10, gedp, lcache, false, clear_images, soft_fail,
	"faceplate_clear", "fp");

    interactive_value.bool_value = 0;
    if (ged_view_context_display_property_set(v,
	    "view.interactive.rectangle.visible", &interactive_value) !=
	BV_DISPLAY_PROPERTY_OK)
	bu_exit(EXIT_FAILURE, "failed to hide interactive rectangle\n");
    ret += img_cmp(0, gedp, lcache, false, clear_images, soft_fail, 0,
	"faceplate_clear", "fp");
    bu_log("Done.\n");

    /***** Framebuffer *****/
    bu_log("Testing framebuffer...\n");
    struct bu_vls fb_img = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fb_img, "%s/moss.png", av[1]);
    if (ged_view_framebuffer_backend_ensure(gedp, v) !=
	BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "failed to initialize Obol framebuffer backend\n");
    icv_image_t *fb_source = icv_read(bu_vls_cstr(&fb_img),
	BU_MIME_IMAGE_PNG, 0, 0);
    if (!fb_source)
	bu_exit(EXIT_FAILURE, "failed to read framebuffer baseline image\n");
    /* Match png-fb's default conversion of untagged input into a linear
     * framebuffer. */
    fb_source->gamma_corr = 0.5f;
    unsigned char *fb_pixels = icv_data2uchar(fb_source);
    const char *png_fb_av[3] = {
	"png2fb", bu_vls_cstr(&fb_img), NULL
    };
    if (!fb_pixels || ged_exec(gedp, 2, png_fb_av) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "png2fb failed to publish framebuffer image: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));

    struct bu_vls fb_capture_path = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fb_capture_path, "%s/faceplate_fb_capture.png", lcache);
    const char *fb_sg_av[4] = {
	"screengrab", "-F", bu_vls_cstr(&fb_capture_path), NULL
    };
    if (ged_exec_screengrab(gedp, 3, fb_sg_av) != BRLCAD_OK ||
	!bu_file_exists(bu_vls_cstr(&fb_capture_path), NULL))
	bu_exit(EXIT_FAILURE, "endpoint framebuffer screengrab failed: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));

    bobol_display_endpoint_t *endpoint =
	ged_view_context_obol_endpoint_get(v);
    const size_t expected_size = fb_source->width * fb_source->height * 3;
    if (!framebuffer_matches(endpoint, fb_pixels, expected_size,
	    (unsigned int)fb_source->width,
	    (unsigned int)fb_source->height))
	bu_exit(EXIT_FAILURE,
	    "png2fb endpoint plane does not match decoded image\n");

    struct bu_vls fb_pix_path = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&fb_pix_path, "%s/faceplate_fb.pix", lcache);
    bu_file_delete(bu_vls_cstr(&fb_pix_path));
    char fb_width[32] = {0};
    char fb_height[32] = {0};
    snprintf(fb_width, sizeof(fb_width), "%zu", fb_source->width);
    snprintf(fb_height, sizeof(fb_height), "%zu", fb_source->height);
    const char *fb2pix_av[7] = {
	"fb2pix", "-w", fb_width, "-n", fb_height,
	bu_vls_cstr(&fb_pix_path), NULL
    };
    if (ged_exec(gedp, 6, fb2pix_av) != BRLCAD_OK ||
	!bu_file_exists(bu_vls_cstr(&fb_pix_path), NULL))
	bu_exit(EXIT_FAILURE, "fb2pix endpoint export failed: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));

    const char *fbclear_av[2] = {"fbclear", NULL};
    if (ged_exec(gedp, 1, fbclear_av) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "fbclear failed: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));
    const char *pix2fb_av[7] = {
	"pix2fb", "-w", fb_width, "-n", fb_height,
	bu_vls_cstr(&fb_pix_path), NULL
    };
    if (ged_exec(gedp, 6, pix2fb_av) != BRLCAD_OK ||
	!framebuffer_matches(endpoint, fb_pixels, expected_size,
	    (unsigned int)fb_source->width,
	    (unsigned int)fb_source->height))
	bu_exit(EXIT_FAILURE, "pix2fb endpoint round trip failed: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));

    if (ged_exec(gedp, 1, fbclear_av) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "pre-overlay fbclear failed\n");
    const char *overlay_fb_av[4] = {
	"overlay", "-F", bu_vls_cstr(&fb_img), NULL
    };
    if (ged_exec(gedp, 3, overlay_fb_av) != BRLCAD_OK ||
	!framebuffer_matches(endpoint, fb_pixels, expected_size,
	    (unsigned int)fb_source->width,
	    (unsigned int)fb_source->height))
	bu_exit(EXIT_FAILURE, "overlay -F endpoint import failed: %s\n",
	    bu_vls_cstr(gedp->ged_result_str));

    bu_file_delete(bu_vls_cstr(&fb_pix_path));
    bu_vls_free(&fb_pix_path);
    bu_file_delete(bu_vls_cstr(&fb_capture_path));
    bu_vls_free(&fb_capture_path);
    bu_free(fb_pixels, "faceplate framebuffer pixels");
    icv_destroy(fb_source);
    bu_vls_free(&fb_img);

    s_av[0] = "view";
    s_av[1] = "faceplate";
    s_av[2] = "fb";
    s_av[3] = "1";
    s_av[4] = NULL;
    ged_exec_view(gedp, 4, s_av);
    ret += img_cmp(9, gedp, lcache, false, clear_images, soft_fail,
	ADIFF_THRES, "faceplate_clear", "fp");

    // Check that turning off works
    s_av[3] = "0";
    ged_exec_view(gedp, 4, s_av);
    ret += img_cmp(0, gedp, lcache, false, clear_images, soft_fail, 0, "faceplate_clear", "fp");

    // Re-enable and make sure clear works
    s_av[3] = "1";
    ged_exec_view(gedp, 4, s_av);
    ret += img_cmp(9, gedp, lcache, false, clear_images, soft_fail,
	ADIFF_THRES, "faceplate_clear", "fp");

    if (ged_exec(gedp, 1, fbclear_av) != BRLCAD_OK)
	bu_exit(EXIT_FAILURE, "failed to clear Obol framebuffer backend\n");
    ret += img_cmp(0, gedp, lcache, false, clear_images, soft_fail, 0, "faceplate_clear", "fp");

    s_av[0] = "view";
    s_av[1] = "faceplate";
    s_av[2] = "fb";
    s_av[3] = "0";
    s_av[4] = NULL;
    ged_exec_view(gedp, 4, s_av);

    bu_log("Done.\n");


    ged_view_framebuffer_release(gedp);
    unsigned char *captured_pixels = NULL;
    size_t captured_size = 0;
    unsigned int captured_width = 0;
    unsigned int captured_height = 0;
    unsigned int captured_components = 0;
    if (bobol_display_endpoint_capture_plane(endpoint,
	    BOBOL_CAPTURE_FRAMEBUFFER, &captured_pixels, &captured_size,
	    &captured_width, &captured_height, &captured_components)) {
	if (captured_pixels)
	    bu_free(captured_pixels, "unexpected framebuffer capture");
	bu_exit(EXIT_FAILURE,
	    "released framebuffer bridge left an endpoint capture callback\n");
    }
    ged_close(gedp);

    if (!keep_images)
	bu_dirclear(lcache);

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
