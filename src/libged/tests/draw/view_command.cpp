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

#include "ged/display_obol_private.h"

#include <cmath>
#include <cstdio>
#include <string>

#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include "bv.h"
#include "BObol/BDisplayEndpoint.h"
#include "BObol/BViewController.h"
#include <bu.h>
#include <ged.h>
#include <ged/draw.h>
#include <ged/display.h>
#include <ged/view_feature_batch.h>
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
assert_feature_point_count(struct ged_view_context *view, const char *name,
    size_t expected, int line)
{
    nchecks++;
    point_t *points = NULL;
    size_t point_count = 0;
    int copied = ged_view_feature_points_copy(view, name,
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
assert_feature_overlay(struct ged_view_context *view, const char *name,
    int expected, int line)
{
    nchecks++;
    struct ged_view_feature_summary summary =
	GED_VIEW_FEATURE_SUMMARY_INIT;
    int copied = ged_view_feature_get_summary(view, name, &summary);
    if (!copied || !summary.exists || summary.is_overlay != expected) {
	bu_log("FAIL [%s:%d] %s overlay expected %d, got exists=%d overlay=%d\n",
		__FILE__, line, name, expected, summary.exists,
		summary.is_overlay);
	nfails++;
    }
}

#define ASSERT_FEATURE_OVERLAY(view, name, expected) assert_feature_overlay((view), (name), (expected), __LINE__)

static void
refresh_scene_records(struct ged *gedp, struct ged_view_context *v)
{
    struct ged_scene_redraw_request request;
    ged_scene_redraw_request_init(&request);
    request.view = v;
    ASSERT(ged_scene_redraw(gedp, &request, NULL) == GED_SCENE_OK);
}

static void
assert_endpoint_background_color(bobol_display_endpoint_t *endpoint,
	const char *property_name, double r, double g, double b)
{
    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    ASSERT(bobol_display_endpoint_property_get(endpoint, property_name,
	&value) == BV_DISPLAY_PROPERTY_OK);
    ASSERT(value.type == BV_DISPLAY_PROPERTY_COLOR3);
    ASSERT(std::fabs(value.color3[0] - r) < 0.0001);
    ASSERT(std::fabs(value.color3[1] - g) < 0.0001);
    ASSERT(std::fabs(value.color3[2] - b) < 0.0001);
}

static void
test_display_endpoint_slot(struct ged *gedp, struct ged_view_context *view)
{
    struct bv_background_state background = BV_BACKGROUND_STATE_INIT;
    VSET(background.bottom, 12, 34, 56);
    VSET(background.top, 78, 90, 123);
    ASSERT(bv_background_state_set(DRAW_TEST_BV(view), &background));

    /* A policy selected while the view is headless must initialize a later
     * endpoint rather than being replaced by the controller default. */
    struct bv_shading_state shading = BV_SHADING_STATE_INIT;
    shading.normal_style = BV_NORMAL_SMOOTH;
    shading.normal_crease_angle = 37.5;
    ASSERT(bv_shading_state_set(DRAW_TEST_BV(view), &shading));

    bobol_display_endpoint_t *owned =
	bobol_display_endpoint_create(NULL, 0);
    ASSERT(owned != NULL);
    void *owned_controller = bobol_display_endpoint_controller(owned);
    ASSERT(owned_controller != NULL);
    ASSERT(ged_view_context_obol_endpoint_set(view, owned, 1));
    ASSERT(ged_view_context_obol_endpoint_get(view) == owned);
    ASSERT(bobol_display_endpoint_controller(
	ged_view_context_obol_endpoint_get(view)) == owned_controller);
    BObolViewController *view_controller =
	static_cast<BObolViewController *>(owned_controller);
    ASSERT(view_controller->getNormalStyle() ==
	BObolViewLodState::NORMAL_SMOOTH);
    ASSERT(std::fabs(view_controller->getNormalCreaseAngle() - 37.5f) <
	0.0001f);
    assert_endpoint_background_color(owned, "controller.background.bottom",
	12.0 / 255.0, 34.0 / 255.0, 56.0 / 255.0);
    assert_endpoint_background_color(owned, "controller.background.top",
	78.0 / 255.0, 90.0 / 255.0, 123.0 / 255.0);

    struct bv_display_property_value background_property =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    background_property.type = BV_DISPLAY_PROPERTY_COLOR3;
    VSET(background_property.color3, 0.2, 0.4, 0.6);
    ASSERT(bobol_display_endpoint_property_set(owned,
	"controller.background.bottom", &background_property) ==
	BV_DISPLAY_PROPERTY_OK);
    VSET(background_property.color3, 0.1, 0.3, 0.5);
    ASSERT(bobol_display_endpoint_property_set(owned,
	"controller.background.top", &background_property) ==
	BV_DISPLAY_PROPERTY_OK);

    /* Re-registering an endpoint must not import stale passive view state. */
    ASSERT(ged_view_context_obol_endpoint_set(view, owned, 1));
    assert_endpoint_background_color(owned, "controller.background.bottom",
	0.2, 0.4, 0.6);
    assert_endpoint_background_color(owned, "controller.background.top",
	0.1, 0.3, 0.5);

    struct bv_display_property_value perspective =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    perspective.type = BV_DISPLAY_PROPERTY_DOUBLE;
    perspective.double_value = 45.0;
    ASSERT(ged_view_context_display_property_set(view, "view.perspective",
	&perspective) == BV_DISPLAY_PROPERTY_OK);
    ASSERT(std::fabs(bv_perspective_get(DRAW_TEST_BV(view)) - 45.0) <
	0.0001);
    perspective = BV_DISPLAY_PROPERTY_VALUE_INIT;
    ASSERT(ged_view_context_display_property_get(view, "view.perspective",
	&perspective) == BV_DISPLAY_PROPERTY_OK);
    ASSERT(perspective.type == BV_DISPLAY_PROPERTY_DOUBLE);
    ASSERT(std::fabs(perspective.double_value - 45.0) < 0.0001);
    ASSERT(bobol_display_endpoint_view_sync(owned, view));
    ASSERT(static_cast<BObolViewController *>(owned_controller)->getCamera()
	->isOfType(SoPerspectiveCamera::getClassTypeId()));
    perspective.type = BV_DISPLAY_PROPERTY_DOUBLE;
    perspective.double_value = -1.0;
    ASSERT(ged_view_context_display_property_set(view, "view.perspective",
	&perspective) == BV_DISPLAY_PROPERTY_INVALID);
    perspective.double_value = 0.0;
    ASSERT(ged_view_context_display_property_set(view, "view.perspective",
	&perspective) == BV_DISPLAY_PROPERTY_OK);
    mat_t pmat;
    ASSERT(bv_pmat_get(pmat, DRAW_TEST_BV(view)));
    ASSERT(std::fabs(pmat[0] - 1.0) < 0.0001 &&
	std::fabs(pmat[5] - 1.0) < 0.0001 &&
	std::fabs(pmat[14]) < 0.0001);

    struct bv_display_property_value navigation =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    navigation.type = BV_DISPLAY_PROPERTY_DOUBLE;
    navigation.double_value = -12.0;
    ASSERT(ged_view_context_display_property_set(view,
	"view.navigation.min_delta", &navigation) ==
	BV_DISPLAY_PROPERTY_OK);
    navigation.double_value = 12.0;
    ASSERT(ged_view_context_display_property_set(view,
	"view.navigation.max_delta", &navigation) ==
	BV_DISPLAY_PROPERTY_OK);
    navigation.double_value = 0.6;
    ASSERT(ged_view_context_display_property_set(view,
	"view.navigation.rotate_scale", &navigation) ==
	BV_DISPLAY_PROPERTY_OK);
    navigation.double_value = 3.0;
    ASSERT(ged_view_context_display_property_set(view,
	"view.navigation.scale_scale", &navigation) ==
	BV_DISPLAY_PROPERTY_OK);
    navigation = BV_DISPLAY_PROPERTY_VALUE_INIT;
    ASSERT(ged_view_context_display_property_get(view,
	"view.navigation.rotate_scale", &navigation) ==
	BV_DISPLAY_PROPERTY_OK &&
	navigation.type == BV_DISPLAY_PROPERTY_DOUBLE &&
	std::fabs(navigation.double_value - 0.6) < 0.0001);
    navigation.type = BV_DISPLAY_PROPERTY_DOUBLE;
    navigation.double_value = 0.0;
    ASSERT(ged_view_context_display_property_set(view,
	"view.navigation.scale_scale", &navigation) ==
	BV_DISPLAY_PROPERTY_INVALID);

    const char *font_properties[] = {
	"view.faceplate.params.font_size",
	"view.faceplate.center_dot.font_size",
	"view.faceplate.scale.font_size"
    };
    struct bv_display_property_value font_size =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    font_size.type = BV_DISPLAY_PROPERTY_UINT;
    font_size.uint_value = 18;
    for (size_t i = 0; i < sizeof(font_properties) / sizeof(font_properties[0]); i++) {
	ASSERT(ged_view_context_display_property_set(view, font_properties[i],
	    &font_size) == BV_DISPLAY_PROPERTY_OK);
	font_size = BV_DISPLAY_PROPERTY_VALUE_INIT;
	ASSERT(ged_view_context_display_property_get(view, font_properties[i],
	    &font_size) == BV_DISPLAY_PROPERTY_OK);
	ASSERT(font_size.type == BV_DISPLAY_PROPERTY_UINT);
	ASSERT(font_size.uint_value == 18);
	font_size.type = BV_DISPLAY_PROPERTY_UINT;
	font_size.uint_value = 97;
	ASSERT(ged_view_context_display_property_set(view, font_properties[i],
	    &font_size) == BV_DISPLAY_PROPERTY_INVALID);
	font_size.uint_value = 18;
    }

    const char *color_properties[] = {
	"view.faceplate.params.color",
	"view.faceplate.center_dot.color",
	"view.faceplate.scale.color"
    };
    struct bv_display_property_value color =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    color.type = BV_DISPLAY_PROPERTY_COLOR3;
    VSET(color.color3, 17.0 / 255.0, 34.0 / 255.0, 51.0 / 255.0);
    for (size_t i = 0; i < sizeof(color_properties) / sizeof(color_properties[0]); i++) {
	ASSERT(ged_view_context_display_property_set(view, color_properties[i],
	    &color) == BV_DISPLAY_PROPERTY_OK);
	color = BV_DISPLAY_PROPERTY_VALUE_INIT;
	ASSERT(ged_view_context_display_property_get(view, color_properties[i],
	    &color) == BV_DISPLAY_PROPERTY_OK);
	ASSERT(color.type == BV_DISPLAY_PROPERTY_COLOR3);
	ASSERT(std::fabs(color.color3[0] - 17.0 / 255.0) < 0.0001);
	ASSERT(std::fabs(color.color3[1] - 34.0 / 255.0) < 0.0001);
	ASSERT(std::fabs(color.color3[2] - 51.0 / 255.0) < 0.0001);
	color.type = BV_DISPLAY_PROPERTY_COLOR3;
	VSET(color.color3, 1.1, 0.0, 0.0);
	ASSERT(ged_view_context_display_property_set(view, color_properties[i],
	    &color) == BV_DISPLAY_PROPERTY_INVALID);
	VSET(color.color3, 17.0 / 255.0, 34.0 / 255.0, 51.0 / 255.0);
    }

    const char *adc_color_properties[] = {
	"view.faceplate.adc.line_color",
	"view.faceplate.adc.tick_color"
    };
    for (size_t i = 0; i < sizeof(adc_color_properties) /
	 sizeof(adc_color_properties[0]); i++) {
	ASSERT(ged_view_context_display_property_set(view, adc_color_properties[i],
	    &color) == BV_DISPLAY_PROPERTY_OK);
	color = BV_DISPLAY_PROPERTY_VALUE_INIT;
	ASSERT(ged_view_context_display_property_get(view, adc_color_properties[i],
	    &color) == BV_DISPLAY_PROPERTY_OK);
	ASSERT(color.type == BV_DISPLAY_PROPERTY_COLOR3);
	ASSERT(std::fabs(color.color3[0] - 17.0 / 255.0) < 0.0001);
	ASSERT(std::fabs(color.color3[1] - 34.0 / 255.0) < 0.0001);
	ASSERT(std::fabs(color.color3[2] - 51.0 / 255.0) < 0.0001);
	color.type = BV_DISPLAY_PROPERTY_COLOR3;
	VSET(color.color3, 17.0 / 255.0, 34.0 / 255.0, 51.0 / 255.0);
    }
    struct bv_display_property_value adc_line_width =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    adc_line_width.type = BV_DISPLAY_PROPERTY_UINT;
    adc_line_width.uint_value = 3;
    ASSERT(ged_view_context_display_property_set(view,
	"view.faceplate.adc.line_width", &adc_line_width) ==
	BV_DISPLAY_PROPERTY_OK);
    adc_line_width = BV_DISPLAY_PROPERTY_VALUE_INIT;
    ASSERT(ged_view_context_display_property_get(view,
	"view.faceplate.adc.line_width", &adc_line_width) ==
	BV_DISPLAY_PROPERTY_OK);
    ASSERT(adc_line_width.type == BV_DISPLAY_PROPERTY_UINT);
    ASSERT(adc_line_width.uint_value == 3);

    const char *faceplate_visibility_properties[] = {
	"view.faceplate.adc.visible",
	"view.faceplate.center_dot.visible",
	"view.faceplate.grid.visible",
	"view.faceplate.model_axes.visible",
	"view.faceplate.scale.visible",
	"view.faceplate.view_axes.visible",
	"view.faceplate.params.visible",
	"view.faceplate.params.size",
	"view.faceplate.params.center",
	"view.faceplate.params.azimuth",
	"view.faceplate.params.elevation",
	"view.faceplate.params.twist",
	"view.faceplate.params.fps"
    };
    struct bv_display_property_value visibility =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    visibility.type = BV_DISPLAY_PROPERTY_BOOL;
    visibility.bool_value = 0;
    for (size_t i = 0;
	i < sizeof(faceplate_visibility_properties) /
	    sizeof(faceplate_visibility_properties[0]); i++) {
	ASSERT(ged_view_context_display_property_set(view,
	    faceplate_visibility_properties[i], &visibility) ==
	    BV_DISPLAY_PROPERTY_OK);
	visibility = BV_DISPLAY_PROPERTY_VALUE_INIT;
	ASSERT(ged_view_context_display_property_get(view,
	    faceplate_visibility_properties[i], &visibility) ==
	    BV_DISPLAY_PROPERTY_OK);
	ASSERT(visibility.type == BV_DISPLAY_PROPERTY_BOOL);
	ASSERT(visibility.bool_value == 0);
	visibility.type = BV_DISPLAY_PROPERTY_BOOL;
	visibility.bool_value = 2;
	ASSERT(ged_view_context_display_property_set(view,
	    faceplate_visibility_properties[i], &visibility) ==
	    BV_DISPLAY_PROPERTY_INVALID);
	visibility.bool_value = 0;
    }

    const char *adc_enable[] = {"adc", "draw", "1", NULL};
    ASSERT(ged_exec_adc(gedp, 3, adc_enable) == BRLCAD_OK);
    visibility = BV_DISPLAY_PROPERTY_VALUE_INIT;
    ASSERT(ged_view_context_display_property_get(view,
	"view.faceplate.adc.visible", &visibility) ==
	BV_DISPLAY_PROPERTY_OK);
    ASSERT(visibility.type == BV_DISPLAY_PROPERTY_BOOL);
    ASSERT(visibility.bool_value == 1);
    const char *adc_disable[] = {"adc", "draw", "0", NULL};
    ASSERT(ged_exec_adc(gedp, 3, adc_disable) == BRLCAD_OK);
    visibility = BV_DISPLAY_PROPERTY_VALUE_INIT;
    ASSERT(ged_view_context_display_property_get(view,
	"view.faceplate.adc.visible", &visibility) ==
	BV_DISPLAY_PROPERTY_OK);
    ASSERT(visibility.bool_value == 0);

    const char *grid_bool_properties[] = {
	"view.faceplate.grid.adaptive",
	"view.faceplate.grid.snap"
    };
    for (size_t i = 0; i < sizeof(grid_bool_properties) /
	 sizeof(grid_bool_properties[0]); i++) {
	visibility = BV_DISPLAY_PROPERTY_VALUE_INIT;
	visibility.type = BV_DISPLAY_PROPERTY_BOOL;
	visibility.bool_value = 1;
	ASSERT(ged_view_context_display_property_set(view,
	    grid_bool_properties[i], &visibility) ==
	    BV_DISPLAY_PROPERTY_OK);
	visibility = BV_DISPLAY_PROPERTY_VALUE_INIT;
	ASSERT(ged_view_context_display_property_get(view,
	    grid_bool_properties[i], &visibility) ==
	    BV_DISPLAY_PROPERTY_OK);
	ASSERT(visibility.type == BV_DISPLAY_PROPERTY_BOOL);
	ASSERT(visibility.bool_value == 1);
    }
    const char *grid_double_properties[] = {
	"view.faceplate.grid.anchor.x",
	"view.faceplate.grid.anchor.y",
	"view.faceplate.grid.anchor.z",
	"view.faceplate.grid.resolution.horizontal",
	"view.faceplate.grid.resolution.vertical"
    };
    for (size_t i = 0; i < sizeof(grid_double_properties) /
	 sizeof(grid_double_properties[0]); i++) {
	struct bv_display_property_value grid_value =
	    BV_DISPLAY_PROPERTY_VALUE_INIT;
	grid_value.type = BV_DISPLAY_PROPERTY_DOUBLE;
	grid_value.double_value = 2.0 + static_cast<double>(i);
	ASSERT(ged_view_context_display_property_set(view,
	    grid_double_properties[i], &grid_value) ==
	    BV_DISPLAY_PROPERTY_OK);
	grid_value = BV_DISPLAY_PROPERTY_VALUE_INIT;
	ASSERT(ged_view_context_display_property_get(view,
	    grid_double_properties[i], &grid_value) ==
	    BV_DISPLAY_PROPERTY_OK);
	ASSERT(grid_value.type == BV_DISPLAY_PROPERTY_DOUBLE);
	ASSERT(std::fabs(grid_value.double_value -
	    (2.0 + static_cast<double>(i))) < 0.0001);
    }
    const char *grid_uint_properties[] = {
	"view.faceplate.grid.major.horizontal",
	"view.faceplate.grid.major.vertical"
    };
    for (size_t i = 0; i < sizeof(grid_uint_properties) /
	 sizeof(grid_uint_properties[0]); i++) {
	struct bv_display_property_value grid_value =
	    BV_DISPLAY_PROPERTY_VALUE_INIT;
	grid_value.type = BV_DISPLAY_PROPERTY_UINT;
	grid_value.uint_value = 3 + i;
	ASSERT(ged_view_context_display_property_set(view,
	    grid_uint_properties[i], &grid_value) ==
	    BV_DISPLAY_PROPERTY_OK);
	grid_value = BV_DISPLAY_PROPERTY_VALUE_INIT;
	ASSERT(ged_view_context_display_property_get(view,
	    grid_uint_properties[i], &grid_value) ==
	    BV_DISPLAY_PROPERTY_OK);
	ASSERT(grid_value.type == BV_DISPLAY_PROPERTY_UINT);
	ASSERT(grid_value.uint_value == 3 + i);
    }
    struct bv_display_property_value grid_color =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    grid_color.type = BV_DISPLAY_PROPERTY_COLOR3;
    VSET(grid_color.color3, 0.1, 0.2, 0.3);
    ASSERT(ged_view_context_display_property_set(view,
	"view.faceplate.grid.color", &grid_color) ==
	BV_DISPLAY_PROPERTY_OK);
    grid_color = BV_DISPLAY_PROPERTY_VALUE_INIT;
    ASSERT(ged_view_context_display_property_get(view,
	"view.faceplate.grid.color", &grid_color) ==
	BV_DISPLAY_PROPERTY_OK);
    ASSERT(grid_color.type == BV_DISPLAY_PROPERTY_COLOR3);
    ASSERT(std::fabs(grid_color.color3[0] - 26.0 / 255.0) < 0.0001);
    ASSERT(std::fabs(grid_color.color3[1] - 51.0 / 255.0) < 0.0001);
    ASSERT(std::fabs(grid_color.color3[2] - 77.0 / 255.0) < 0.0001);

    const char *grid_snap[] = {"grid", "snap", "0", NULL};
    ASSERT(ged_exec_grid(gedp, 3, grid_snap) == BRLCAD_OK);
    visibility = BV_DISPLAY_PROPERTY_VALUE_INIT;
    ASSERT(ged_view_context_display_property_get(view,
	"view.faceplate.grid.snap", &visibility) ==
	BV_DISPLAY_PROPERTY_OK);
    ASSERT(visibility.bool_value == 0);
    const char *grid_anchor[] = {"grid", "anchor", "7", "8", "9", NULL};
    ASSERT(ged_exec_grid(gedp, 5, grid_anchor) == BRLCAD_OK);
    struct bv_display_property_value anchor =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    ASSERT(ged_view_context_display_property_get(view,
	"view.faceplate.grid.anchor.z", &anchor) ==
	BV_DISPLAY_PROPERTY_OK);
    ASSERT(std::fabs(anchor.double_value - 9.0) < 0.0001);
    const char *grid_color_command[] = {
	"grid", "color", "10", "20", "30", NULL
    };
    ASSERT(ged_exec_grid(gedp, 5, grid_color_command) == BRLCAD_OK);
    grid_color = BV_DISPLAY_PROPERTY_VALUE_INIT;
    ASSERT(ged_view_context_display_property_get(view,
	"view.faceplate.grid.color", &grid_color) ==
	BV_DISPLAY_PROPERTY_OK);
    ASSERT(std::fabs(grid_color.color3[0] - 10.0 / 255.0) < 0.0001);

    const char *faceplate_grid_color[] = {
	"view", "faceplate", "grid", "color", "11", "22", "33", NULL
    };
    ASSERT_VIEW_OK(gedp, 7, faceplate_grid_color);
    grid_color = BV_DISPLAY_PROPERTY_VALUE_INIT;
    ASSERT(ged_view_context_display_property_get(view,
	"view.faceplate.grid.color", &grid_color) ==
	BV_DISPLAY_PROPERTY_OK);
    ASSERT(std::fabs(grid_color.color3[2] - 33.0 / 255.0) < 0.0001);

    const char *axes_names[] = {"model_axes", "view_axes"};
    const char *axes_bool_suffixes[] = {
	"position_only", "labels.visible", "triple_color", "ticks.visible"
    };
    const char *axes_double_suffixes[] = {
	"position.x", "position.y", "position.z", "size", "ticks.interval"
    };
    const char *axes_uint_suffixes[] = {
	"line_width", "ticks.length", "ticks.major_length",
	"ticks.per_major", "ticks.threshold"
    };
    const char *axes_color_suffixes[] = {
	"color", "labels.color", "ticks.color", "ticks.major_color"
    };
    for (const char *axes_name : axes_names) {
	for (const char *suffix : axes_bool_suffixes) {
	    const std::string property = "view.faceplate." +
		std::string(axes_name) + "." + suffix;
	    struct bv_display_property_value axes_value =
		BV_DISPLAY_PROPERTY_VALUE_INIT;
	    axes_value.type = BV_DISPLAY_PROPERTY_BOOL;
	    axes_value.bool_value = 1;
	    ASSERT(ged_view_context_display_property_set(view, property.c_str(),
		&axes_value) == BV_DISPLAY_PROPERTY_OK);
	    axes_value = BV_DISPLAY_PROPERTY_VALUE_INIT;
	    ASSERT(ged_view_context_display_property_get(view, property.c_str(),
		&axes_value) == BV_DISPLAY_PROPERTY_OK);
	    ASSERT(axes_value.type == BV_DISPLAY_PROPERTY_BOOL);
	    ASSERT(axes_value.bool_value == 1);
	}
	for (size_t i = 0; i < sizeof(axes_double_suffixes) /
	     sizeof(axes_double_suffixes[0]); i++) {
	    const std::string property = "view.faceplate." +
		std::string(axes_name) + "." + axes_double_suffixes[i];
	    struct bv_display_property_value axes_value =
		BV_DISPLAY_PROPERTY_VALUE_INIT;
	    axes_value.type = BV_DISPLAY_PROPERTY_DOUBLE;
	    axes_value.double_value = 2.0 + static_cast<double>(i);
	    ASSERT(ged_view_context_display_property_set(view, property.c_str(),
		&axes_value) == BV_DISPLAY_PROPERTY_OK);
	    axes_value = BV_DISPLAY_PROPERTY_VALUE_INIT;
	    ASSERT(ged_view_context_display_property_get(view, property.c_str(),
		&axes_value) == BV_DISPLAY_PROPERTY_OK);
	    ASSERT(axes_value.type == BV_DISPLAY_PROPERTY_DOUBLE);
	    ASSERT(std::fabs(axes_value.double_value -
		(2.0 + static_cast<double>(i))) < 0.0001);
	}
	for (size_t i = 0; i < sizeof(axes_uint_suffixes) /
	     sizeof(axes_uint_suffixes[0]); i++) {
	    const std::string property = "view.faceplate." +
		std::string(axes_name) + "." + axes_uint_suffixes[i];
	    struct bv_display_property_value axes_value =
		BV_DISPLAY_PROPERTY_VALUE_INIT;
	    axes_value.type = BV_DISPLAY_PROPERTY_UINT;
	    axes_value.uint_value = 2 + i;
	    ASSERT(ged_view_context_display_property_set(view, property.c_str(),
		&axes_value) == BV_DISPLAY_PROPERTY_OK);
	    axes_value = BV_DISPLAY_PROPERTY_VALUE_INIT;
	    ASSERT(ged_view_context_display_property_get(view, property.c_str(),
		&axes_value) == BV_DISPLAY_PROPERTY_OK);
	    ASSERT(axes_value.type == BV_DISPLAY_PROPERTY_UINT);
	    ASSERT(axes_value.uint_value == 2 + i);
	}
	for (const char *suffix : axes_color_suffixes) {
	    const std::string property = "view.faceplate." +
		std::string(axes_name) + "." + suffix;
	    struct bv_display_property_value axes_value =
		BV_DISPLAY_PROPERTY_VALUE_INIT;
	    axes_value.type = BV_DISPLAY_PROPERTY_COLOR3;
	    VSET(axes_value.color3, 0.2, 0.4, 0.6);
	    ASSERT(ged_view_context_display_property_set(view, property.c_str(),
		&axes_value) == BV_DISPLAY_PROPERTY_OK);
	    axes_value = BV_DISPLAY_PROPERTY_VALUE_INIT;
	    ASSERT(ged_view_context_display_property_get(view, property.c_str(),
		&axes_value) == BV_DISPLAY_PROPERTY_OK);
	    ASSERT(axes_value.type == BV_DISPLAY_PROPERTY_COLOR3);
	    ASSERT(std::fabs(axes_value.color3[0] - 51.0 / 255.0) < 0.0001);
	    ASSERT(std::fabs(axes_value.color3[2] - 153.0 / 255.0) < 0.0001);
	}
    }
    const char *faceplate_model_axes_color[] = {
	"view", "faceplate", "model_axes", "axes_color", "11", "22", "33", NULL
    };
    ASSERT_VIEW_OK(gedp, 7, faceplate_model_axes_color);
    struct bv_display_property_value axes_color =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    ASSERT(ged_view_context_display_property_get(view,
	"view.faceplate.model_axes.color", &axes_color) ==
	BV_DISPLAY_PROPERTY_OK);
    ASSERT(std::fabs(axes_color.color3[1] - 22.0 / 255.0) < 0.0001);
    const char *faceplate_view_axes_tick[] = {
	"view", "faceplate", "view_axes", "tick_length", "9", NULL
    };
    ASSERT_VIEW_OK(gedp, 5, faceplate_view_axes_tick);
    struct bv_display_property_value axes_tick =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    ASSERT(ged_view_context_display_property_get(view,
	"view.faceplate.view_axes.ticks.length", &axes_tick) ==
	BV_DISPLAY_PROPERTY_OK);
    ASSERT(axes_tick.uint_value == 9);

    bobol_display_endpoint_t *borrowed =
	bobol_display_endpoint_create(NULL, 0);
    ASSERT(borrowed != NULL);
    ASSERT(ged_view_context_obol_endpoint_set(view, borrowed, 0));
    ASSERT(ged_view_context_obol_endpoint_get(view) == borrowed);
    ASSERT(bobol_display_endpoint_controller(
	ged_view_context_obol_endpoint_get(view)) ==
	bobol_display_endpoint_controller(borrowed));
    ASSERT(ged_view_context_obol_endpoint_set(view, NULL, 0));
    bobol_display_endpoint_destroy(borrowed);

    /* Restore this view's production default for the following command
     * coverage. */
    shading = BV_SHADING_STATE_INIT;
    ASSERT(bv_shading_state_set(DRAW_TEST_BV(view), &shading));
}

static void
test_command_report_record_consistency(struct ged *gedp,
    struct ged_view_context *v)
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

    struct ged_pick_result *pick =
	ged_pick_semantic_path(gedp, v, "all.g");
    struct bu_vls pick_path = BU_VLS_INIT_ZERO;
    ASSERT(pick != NULL);
    ASSERT(ged_pick_result_count(pick) > 0);
    if (pick && ged_pick_result_count(pick) > 0) {
	ASSERT(ged_pick_result_path(pick, 0, &pick_path));
	ASSERT(BU_STR_EQUAL(bu_vls_cstr(&pick_path), "all.g"));
    }

    point_t sample = VINIT_ZERO;
    point_t snap_candidate = VINIT_ZERO;
    int snap_count = ged_view_selection_snap(v, sample,
	    GED_SELECTION_SNAP_ENDPOINT, snap_candidate);
    ASSERT(snap_count >= 0);
    if (snap_count > 0) {
	ASSERT(std::isfinite((double)snap_candidate[X]));
	ASSERT(std::isfinite((double)snap_candidate[Y]));
	ASSERT(std::isfinite((double)snap_candidate[Z]));
    }

    bu_vls_free(&pick_path);
    ged_pick_result_free(pick);
}

static void
test_shading_policy_command(struct ged *gedp,
    struct ged_view_context *primary, struct ged_view_context *secondary)
{
    struct bv_shading_state shading = BV_SHADING_STATE_INIT;
    ASSERT(bv_shading_state_get(&shading, DRAW_TEST_BV(primary)));
    ASSERT(shading.normal_style == BV_NORMAL_AUTHORED);
    ASSERT(std::fabs(shading.normal_crease_angle - 60.0) < 0.0001);

    const char *status[] = {"view", "shading", NULL};
    ASSERT(run_view(gedp, 2, status) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("normals authored") != std::string::npos);
    ASSERT(result_str(gedp).find("crease 60") != std::string::npos);

    const char *smooth[] = {"view", "shading", "normals", "smooth", NULL};
    ASSERT(run_view(gedp, 4, smooth) == BRLCAD_OK);
    ASSERT(bv_shading_state_get(&shading, DRAW_TEST_BV(primary)));
    ASSERT(shading.normal_style == BV_NORMAL_SMOOTH);

    const char *crease[] = {"view", "shading", "crease", "37.5", NULL};
    ASSERT(run_view(gedp, 4, crease) == BRLCAD_OK);
    ASSERT(bv_shading_state_get(&shading, DRAW_TEST_BV(primary)));
    ASSERT(std::fabs(shading.normal_crease_angle - 37.5) < 0.0001);

    const char *flat[] = {"view", "shading", "normals", "flat", NULL};
    ASSERT(run_view(gedp, 4, flat) == BRLCAD_OK);
    ASSERT(bv_shading_state_get(&shading, DRAW_TEST_BV(primary)));
    ASSERT(shading.normal_style == BV_NORMAL_FLAT);

    /* Shading is view-local.  Selecting V1 must not alter V0. */
    const char *secondary_smooth[] = {
	"view", "-V", "V1", "shading", "normals", "smooth", NULL
    };
    ASSERT(run_view(gedp, 6, secondary_smooth) == BRLCAD_OK);
    ASSERT(bv_shading_state_get(&shading, DRAW_TEST_BV(secondary)));
    ASSERT(shading.normal_style == BV_NORMAL_SMOOTH);
    ASSERT(bv_shading_state_get(&shading, DRAW_TEST_BV(primary)));
    ASSERT(shading.normal_style == BV_NORMAL_FLAT);

    const char *bad_style[] = {
	"view", "shading", "normals", "rounded", NULL
    };
    ASSERT(run_view(gedp, 4, bad_style) == BRLCAD_ERROR);
    ASSERT(bv_shading_state_get(&shading, DRAW_TEST_BV(primary)));
    ASSERT(shading.normal_style == BV_NORMAL_FLAT);

    const char *bad_crease[] = {
	"view", "shading", "crease", "181", NULL
    };
    ASSERT(run_view(gedp, 4, bad_crease) == BRLCAD_ERROR);
    ASSERT(bv_shading_state_get(&shading, DRAW_TEST_BV(primary)));
    ASSERT(std::fabs(shading.normal_crease_angle - 37.5) < 0.0001);

    /* Leave the test view at the production default. */
    const char *authored[] = {
	"view", "shading", "normals", "authored", NULL
    };
    const char *default_crease[] = {
	"view", "shading", "crease", "60", NULL
    };
    ASSERT(run_view(gedp, 4, authored) == BRLCAD_OK);
    ASSERT(run_view(gedp, 4, default_crease) == BRLCAD_OK);
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

    struct ged_view_set *view_set_ctx = ged_view_set_ctx(gedp);
    ASSERT(ged_view_set_context_remove(view_set_ctx, NULL));
    struct ged_view_context *views[2] = {NULL, NULL};
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

    test_display_endpoint_slot(gedp, views[0]);
    test_shading_policy_command(gedp, views[0], views[1]);
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
    ASSERT(ged_view_feature_exists(views[0], "u_line_edit") == 0);

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

    const char *c1json[] = {"view", "feature", "--json", "info", "u_line", NULL};
    ASSERT(run_view(gedp, 5, c1json) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("{\"name\":\"u_line\"") == 0);
    ASSERT(result_str(gedp).find("\"kind\":\"lines\"") != std::string::npos);
    ASSERT(result_str(gedp).find("\"scope\":\"shared\"") != std::string::npos);

    const char *c1missing[] = {"view", "feature", "info", "no_such_feature", NULL};
    ASSERT(run_view(gedp, 4, c1missing) == BRLCAD_ERROR);
    ASSERT(result_str(gedp).find("No view feature named no_such_feature") !=
	std::string::npos);

    struct ged_view_feature_batch_desc result_desc = GED_VIEW_FEATURE_BATCH_DESC_INIT;
    result_desc.owner_id = "view-command-test";
    result_desc.owner_role = "command";
    result_desc.run_id = "run-1";
    result_desc.generation = 1;
    struct ged_view_feature_batch *aborted_batch =
	ged_view_feature_batch_begin(views[0], &result_desc);
    point_t aborted_points[2] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    int aborted_commands[2] = {0, 1};
    ASSERT(aborted_batch != NULL);
    ASSERT(aborted_batch && ged_view_feature_batch_line_set_replace(
	aborted_batch, "aborted_line", aborted_points, aborted_commands, 2,
	NULL));
    ASSERT(!ged_view_feature_exists(views[0], "aborted_line"));
    ged_view_feature_batch_abort(aborted_batch);
    ASSERT(!ged_view_feature_exists(views[0], "aborted_line"));

    struct ged_view_feature_batch *result_scene = ged_view_feature_batch_begin(views[0],
	&result_desc);
    point_t result_points[2] = {{0.0, 0.0, 0.0}, {1.0, 1.0, 0.0}};
    int result_commands[2] = {0, 1};
    ASSERT(result_scene != NULL);
    ASSERT(result_scene && ged_view_feature_batch_line_set_replace(result_scene,
	"result_line", result_points, result_commands, 2, NULL));
    ASSERT(!ged_view_feature_exists(views[0], "result_line"));
    ASSERT(result_scene && ged_view_feature_batch_commit(result_scene));
    const char *c1result[] = {"view", "feature", "--json", "info", "result_line", NULL};
    ASSERT(run_view(gedp, 5, c1result) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("\"command_result\":true") !=
	std::string::npos);

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
    const char *c4json[] = {"view", "feature", "--json", "list", "u_*", NULL};
    ASSERT(run_view(gedp, 5, c4json) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("[\"") == 0);
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
    ASSERT(ged_view_feature_label_count(views[0], "u_label") == 1);
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
    ASSERT(ged_view_feature_exists(views[0], "_tcl_data_lines") == 0);

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
    ged_view_polygon_ref poly_ref =
	ged_view_polygon_find(views[0], "u_poly");
    ASSERT(!ged_view_polygon_ref_is_null(poly_ref));
    struct ged_view_polygon_record poly_rec = {};
    ASSERT(ged_view_polygon_record_get(views[0], poly_ref, &poly_rec));
    ASSERT(poly_rec.type == GED_VIEW_POLYGON_GENERAL);
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

    const char *amb_shared[] = {"view", "annotation", "line", "create",
	"ambiguous_name", "0", "0", "0", "1", "0", "0", NULL};
    ASSERT_VIEW_OK(gedp, 11, amb_shared);
    const char *amb_local[] = {"view", "-V", "V0", "polygon", "-L",
	"create", "ambiguous_name", "10", "10", NULL};
    ASSERT_VIEW_OK(gedp, 9, amb_local);
    const char *amb_info[] = {"view", "feature", "info", "ambiguous_name", NULL};
    ASSERT(run_view(gedp, 4, amb_info) == BRLCAD_ERROR);
    ASSERT(result_str(gedp).find("ambiguous (2 matches:") != std::string::npos);
    const char *amb_local_info[] = {"view", "feature", "-L", "info",
	"ambiguous_name", "scope", NULL};
    ASSERT(run_view(gedp, 6, amb_local_info) == BRLCAD_OK);
    ASSERT(result_str(gedp).find("local") != std::string::npos);

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
