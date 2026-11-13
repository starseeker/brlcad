/*                           A X E S . C
 * BRL-CAD
 *
 * Copyright (c) 2000-2026 United States Government as represented by
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
/** @addtogroup libtclcad */
/** @{ */
/** @file libtclcad/view/axes.c
 *
 */
/** @} */

#include "common.h"

#include "ged/display_obol_private.h"
#include "bv.h"
#include "BObol/BDisplayEndpoint.h"
#include "bu/units.h"
#include "ged.h"
#include "ged/view.h"
#include "ged/view_feature_batch.h"
#include "rt/view.h"
#include "tclcad.h"

/* Private headers */
#include "ged/draw.h"
#include "../tclcad_private.h"
#include "../view/view.h"

/* TclCAD owns data-axes command state.  GED receives immutable retained
 * publication batches and is never queried as a command-state database. */
#define BVDAS_DEFAULT_VIEW_WIDTH 512

static fastf_t
_tclcad_data_axes_display_scale(struct ged_view_context *view_ctx)
{
    const struct bv *view =
	bv_context_view_const((const struct bv_context *)view_ctx);
    const int width = bv_context_width_get(
	(const struct bv_context *)view_ctx);
    return bv_size_get(view) /
	(fastf_t)(width > 0 ? width : BVDAS_DEFAULT_VIEW_WIDTH);
}

int
tclcad_data_axes_publish(struct ged_view_context *view_ctx, const char *name,
	const struct tclcad_data_axes_state *state)
{
    if (!view_ctx || !name || !state)
	return 0;

    struct ged_view_feature_batch_desc desc = ged_view_feature_batch_desc_default();
    desc.owner_id = "tclcad-axes";
    desc.owner_role = "tcl-overlay";
    desc.overlay_class = GED_VIEW_FEATURE_OVERLAY_CLASS_TCL_OVERLAY;
    desc.local = 1;
    struct ged_view_feature_batch *batch =
	ged_view_feature_batch_begin(view_ctx, &desc);
    if (!batch)
	return 0;

    struct ged_view_feature_style style = ged_view_feature_style_default();
    style.visible = state->draw ? 1 : 0;
    style.selectable = 1;
    style.color_valid = 1;
    style.color[0] = (unsigned char)state->color[0];
    style.color[1] = (unsigned char)state->color[1];
    style.color[2] = (unsigned char)state->color[2];
    style.line_width = state->line_width;
    const point_t *centers = state->draw && state->num_points > 0 ?
	(const point_t *)state->points : NULL;
    const size_t count = centers ? (size_t)state->num_points : 0;
    const fastf_t half_size = state->size > 0.0 ?
	state->size * 0.5 * _tclcad_data_axes_display_scale(view_ctx) : 0.0;
    if (!ged_view_feature_batch_axes_replace(batch, name, centers, count,
	    half_size, &style)) {
	ged_view_feature_batch_abort(batch);
	return 0;
    }
    return ged_view_feature_batch_commit(batch);
}

static int
tclcad_axes_visibility_endpoint_set(struct ged_view_context *view_ctx,
	const char *property_name, int enabled)
{
    if (!view_ctx || !property_name || enabled < 0 || enabled > 1)
	return BV_DISPLAY_PROPERTY_INVALID;

    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    value.type = BV_DISPLAY_PROPERTY_BOOL;
    value.bool_value = enabled;
    return ged_view_context_display_property_set(view_ctx, property_name,
	&value);
}

static int
tclcad_axes_endpoint_property_set(struct ged_view_context *view_ctx, const char *name,
	const struct bv_display_property_value *value)
{
    return ged_view_context_display_property_set(view_ctx, name, value) ==
	BV_DISPLAY_PROPERTY_OK;
}

static int
tclcad_axes_endpoint_style_set(struct ged_view_context *view_ctx,
	const char *visibility_property, const char *field,
	const struct bv_axes_state *axes)
{
    if (!view_ctx || !visibility_property || !field || !axes)
	return 0;
    const char *prefix = BU_STR_EQUAL(visibility_property,
	"view.faceplate.model_axes.visible") ?
	"view.faceplate.model_axes." :
	"view.faceplate.view_axes.";
    char property[128] = {0};
    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;

    if (BU_STR_EQUAL(field, "axes_size")) {
	value.type = BV_DISPLAY_PROPERTY_DOUBLE;
	value.double_value = axes->axes_size;
	snprintf(property, sizeof(property), "%ssize", prefix);
	return tclcad_axes_endpoint_property_set(view_ctx, property, &value);
    }
    if (BU_STR_EQUAL(field, "axes_pos")) {
	const char *coordinates[] = {"position.x", "position.y", "position.z"};
	value.type = BV_DISPLAY_PROPERTY_DOUBLE;
	for (int axis = 0; axis < 3; axis++) {
	    value.double_value = axes->axes_pos[axis];
	    snprintf(property, sizeof(property), "%s%s", prefix,
		coordinates[axis]);
	    if (!tclcad_axes_endpoint_property_set(view_ctx, property, &value))
		return 0;
	}
	return 1;
    }
    if (BU_STR_EQUAL(field, "axes_color") ||
	BU_STR_EQUAL(field, "label_color") ||
	BU_STR_EQUAL(field, "tick_color") ||
	BU_STR_EQUAL(field, "tick_major_color")) {
	const int *color = BU_STR_EQUAL(field, "axes_color") ?
	    axes->axes_color : BU_STR_EQUAL(field, "label_color") ?
	    axes->label_color : BU_STR_EQUAL(field, "tick_color") ?
	    axes->tick_color : axes->tick_major_color;
	const char *suffix = BU_STR_EQUAL(field, "axes_color") ? "color" :
	    BU_STR_EQUAL(field, "label_color") ? "labels.color" :
	    BU_STR_EQUAL(field, "tick_color") ? "ticks.color" :
	    "ticks.major_color";
	value.type = BV_DISPLAY_PROPERTY_COLOR3;
	for (int axis = 0; axis < 3; axis++)
	    value.color3[axis] = color[axis] / 255.0;
	snprintf(property, sizeof(property), "%s%s", prefix, suffix);
	return tclcad_axes_endpoint_property_set(view_ctx, property, &value);
    }

    const int is_bool = BU_STR_EQUAL(field, "pos_only") ||
	BU_STR_EQUAL(field, "tick_enable") ||
	BU_STR_EQUAL(field, "triple_color");
    if (is_bool) {
	const int enabled = BU_STR_EQUAL(field, "pos_only") ? axes->pos_only :
	    BU_STR_EQUAL(field, "tick_enable") ? axes->tick_enabled :
	    axes->triple_color;
	const char *suffix = BU_STR_EQUAL(field, "pos_only") ?
	    "position_only" : BU_STR_EQUAL(field, "tick_enable") ?
	    "ticks.visible" : "triple_color";
	value.type = BV_DISPLAY_PROPERTY_BOOL;
	value.bool_value = enabled ? 1 : 0;
	snprintf(property, sizeof(property), "%s%s", prefix, suffix);
	return tclcad_axes_endpoint_property_set(view_ctx, property, &value);
    }

    const int number = BU_STR_EQUAL(field, "line_width") ?
	axes->line_width : BU_STR_EQUAL(field, "tick_length") ?
	axes->tick_length : BU_STR_EQUAL(field, "tick_major_length") ?
	axes->tick_major_length : BU_STR_EQUAL(field, "ticks_per_major") ?
	axes->ticks_per_major : BU_STR_EQUAL(field, "tick_threshold") ?
	axes->tick_threshold : -1;
    if (number >= 0) {
	const char *suffix = BU_STR_EQUAL(field, "line_width") ?
	    "line_width" : BU_STR_EQUAL(field, "tick_length") ?
	    "ticks.length" : BU_STR_EQUAL(field, "tick_major_length") ?
	    "ticks.major_length" : BU_STR_EQUAL(field, "ticks_per_major") ?
	    "ticks.per_major" : "ticks.threshold";
	value.type = BV_DISPLAY_PROPERTY_UINT;
	value.uint_value = (uint64_t)number;
	snprintf(property, sizeof(property), "%s%s", prefix, suffix);
	return tclcad_axes_endpoint_property_set(view_ctx, property, &value);
    }
    if (BU_STR_EQUAL(field, "tick_interval")) {
	value.type = BV_DISPLAY_PROPERTY_DOUBLE;
	value.double_value = axes->tick_interval;
	snprintf(property, sizeof(property), "%sticks.interval", prefix);
	return tclcad_axes_endpoint_property_set(view_ctx, property, &value);
    }
    return 0;
}

int
to_axes(struct ged *gedp,
	struct ged_view_context *view_ctx,
	struct bv_axes_state *gasp,
	const char *visibility_property,
	int argc,
	const char *argv[],
	const char *usage)
{

    if (BU_STR_EQUAL(argv[2], "draw")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->draw);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int i;

	    if (bu_sscanf(argv[3], "%d", &i) != 1)
		goto bad;

	    if (ged_view_context_obol_endpoint_get(view_ctx)) {
		if (tclcad_axes_visibility_endpoint_set(view_ctx,
			visibility_property, i) != BV_DISPLAY_PROPERTY_OK)
		    goto bad;
	    } else {
		gasp->draw = i ? 1 : 0;
	    }

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "axes_size")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%lf", gasp->axes_size);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    double size; /* must be double for scanf */

	    if (bu_sscanf(argv[3], "%lf", &size) != 1)
		goto bad;

	    gasp->axes_size = size;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "axes_pos")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%lf %lf %lf",
			  V3ARGS(gasp->axes_pos));
	    return BRLCAD_OK;
	}

	if (argc == 6) {
	    double x, y, z; /* must be double for scanf */

	    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
		bu_sscanf(argv[4], "%lf", &y) != 1 ||
		bu_sscanf(argv[5], "%lf", &z) != 1)
		goto bad;

	    VSET(gasp->axes_pos, x, y, z);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "axes_color")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d %d %d",
			  V3ARGS(gasp->axes_color));
	    return BRLCAD_OK;
	}

	if (argc == 6) {
	    int r, g, b;

	    /* set background color */
	    if (bu_sscanf(argv[3], "%d", &r) != 1 ||
		bu_sscanf(argv[4], "%d", &g) != 1 ||
		bu_sscanf(argv[5], "%d", &b) != 1)
		goto bad;

	    /* validate color */
	    if (r < 0 || 255 < r ||
		g < 0 || 255 < g ||
		b < 0 || 255 < b)
		goto bad;

	    VSET(gasp->axes_color, r, g, b);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "label_color")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d %d %d",
			  V3ARGS(gasp->label_color));
	    return BRLCAD_OK;
	}

	if (argc == 6) {
	    int r, g, b;

	    /* set background color */
	    if (bu_sscanf(argv[3], "%d", &r) != 1 ||
		bu_sscanf(argv[4], "%d", &g) != 1 ||
		bu_sscanf(argv[5], "%d", &b) != 1)
		goto bad;

	    /* validate color */
	    if (r < 0 || 255 < r ||
		g < 0 || 255 < g ||
		b < 0 || 255 < b)
		goto bad;

	    VSET(gasp->label_color, r, g, b);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "line_width")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->line_width);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int line_width;

	    if (bu_sscanf(argv[3], "%d", &line_width) != 1)
		goto bad;

	    gasp->line_width = line_width;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "pos_only")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->pos_only);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int i;

	    if (bu_sscanf(argv[3], "%d", &i) != 1)
		goto bad;

	    if (i)
		gasp->pos_only = 1;
	    else
		gasp->pos_only = 0;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "tick_color")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d %d %d",
			  V3ARGS(gasp->tick_color));
	    return BRLCAD_OK;
	}

	if (argc == 6) {
	    int r, g, b;

	    /* set background color */
	    if (bu_sscanf(argv[3], "%d", &r) != 1 ||
		bu_sscanf(argv[4], "%d", &g) != 1 ||
		bu_sscanf(argv[5], "%d", &b) != 1)
		goto bad;

	    /* validate color */
	    if (r < 0 || 255 < r ||
		g < 0 || 255 < g ||
		b < 0 || 255 < b)
		goto bad;

	    VSET(gasp->tick_color, r, g, b);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "tick_enable")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->tick_enabled);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int i;

	    if (bu_sscanf(argv[3], "%d", &i) != 1)
		goto bad;

	    if (i)
		gasp->tick_enabled = 1;
	    else
		gasp->tick_enabled = 0;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "tick_interval")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%f", gasp->tick_interval);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int tick_interval;

	    if (bu_sscanf(argv[3], "%d", &tick_interval) != 1)
		goto bad;

	    gasp->tick_interval = tick_interval;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "tick_length")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->tick_length);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int tick_length;

	    if (bu_sscanf(argv[3], "%d", &tick_length) != 1)
		goto bad;

	    gasp->tick_length = tick_length;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "tick_major_color")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d %d %d",
			  V3ARGS(gasp->tick_major_color));
	    return BRLCAD_OK;
	}

	if (argc == 6) {
	    int r, g, b;

	    /* set background color */
	    if (bu_sscanf(argv[3], "%d", &r) != 1 ||
		bu_sscanf(argv[4], "%d", &g) != 1 ||
		bu_sscanf(argv[5], "%d", &b) != 1)
		goto bad;

	    /* validate color */
	    if (r < 0 || 255 < r ||
		g < 0 || 255 < g ||
		b < 0 || 255 < b)
		goto bad;

	    VSET(gasp->tick_major_color, r, g, b);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "tick_major_length")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->tick_major_length);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int tick_major_length;

	    if (bu_sscanf(argv[3], "%d", &tick_major_length) != 1)
		goto bad;

	    gasp->tick_major_length = tick_major_length;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "ticks_per_major")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->ticks_per_major);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int ticks_per_major;

	    if (bu_sscanf(argv[3], "%d", &ticks_per_major) != 1)
		goto bad;

	    gasp->ticks_per_major = ticks_per_major;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "tick_threshold")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->tick_threshold);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int tick_threshold;

	    if (bu_sscanf(argv[3], "%d", &tick_threshold) != 1)
		goto bad;

	    if (tick_threshold < 1)
		tick_threshold = 1;

	    gasp->tick_threshold = tick_threshold;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "triple_color")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->triple_color);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int i;

	    if (bu_sscanf(argv[3], "%d", &i) != 1)
		goto bad;

	    if (i)
		gasp->triple_color = 1;
	    else
		gasp->triple_color = 0;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

bad:
    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
    return BRLCAD_ERROR;
}

int
go_data_axes(Tcl_Interp *interp,
	     struct ged *gedp,
	     struct ged_view_context *draw_view_ctx,
	     int argc,
	     const char *argv[],
	     const char *usage)
{
    int ret;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 2 || 5 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }
    to_refresh_suppress_all_begin(current_top);

    ret = to_data_axes_func(interp, gedp, draw_view_ctx, argc, argv);
    to_refresh_suppress_all_end(current_top);
    if (ret & BRLCAD_ERROR)
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);

    return ret;
}

int
to_data_axes(struct ged *gedp,
	     int argc,
	     const char *argv[],
	     ged_func_ptr UNUSED(func),
	     const char *usage,
	     int UNUSED(maxargs))
{
    struct ged_view_context *view_ctx;
    int ret;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 3 || 6 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    view_ctx = ged_view_find_ctx(gedp, argv[1]);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* shift the command name to argv[1] before calling to_data_axes_func */
    argv[1] = argv[0];
    ret = to_data_axes_func(current_top->to_interp, gedp, view_ctx, argc-1, argv+1);
    if (ret == BRLCAD_ERROR)
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);

    return ret;
}

int
to_data_axes_func(Tcl_Interp *interp,
		  struct ged *gedp,
		  struct ged_view_context *view_ctx,
		  int argc,
		  const char *argv[])
{
    const int staged = argv[0][0] == 's';
    const char *feature_name = staged ? "_tcl_sdata_axes" : "_tcl_data_axes";
    tclcad_view_state *view_state =
	tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!view_state)
	return BRLCAD_ERROR;
    struct tclcad_data_axes_state *state = staged ?
	&view_state->gv_sdata_axes : &view_state->gv_data_axes;

    if (BU_STR_EQUAL(argv[1], "draw")) {
	if (argc == 2) {
	    bu_vls_printf(gedp->ged_result_str, "%d", state->draw);
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int i;

	    if (bu_sscanf(argv[2], "%d", &i) != 1)
		goto bad;

	    state->draw = i ? 1 : 0;
	    (void)tclcad_data_axes_publish(view_ctx, feature_name, state);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "color")) {
	if (argc == 2) {
	    bu_vls_printf(gedp->ged_result_str, "%d %d %d",
		state->color[0], state->color[1], state->color[2]);
	    return BRLCAD_OK;
	}

	if (argc == 5) {
	    int r, g, b;

	    /* set background color */
	    if (bu_sscanf(argv[2], "%d", &r) != 1 ||
		bu_sscanf(argv[3], "%d", &g) != 1 ||
		bu_sscanf(argv[4], "%d", &b) != 1)
		goto bad;

	    /* validate color */
	    if (r < 0 || 255 < r ||
		g < 0 || 255 < g ||
		b < 0 || 255 < b)
		goto bad;

	    VSET(state->color, r, g, b);
	    (void)tclcad_data_axes_publish(view_ctx, feature_name, state);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "line_width")) {
	if (argc == 2) {
	    bu_vls_printf(gedp->ged_result_str, "%d", state->line_width);
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int line_width;

	    if (bu_sscanf(argv[2], "%d", &line_width) != 1)
		goto bad;

	    state->line_width = line_width;
	    (void)tclcad_data_axes_publish(view_ctx, feature_name, state);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "size")) {
	if (argc == 2) {
	    bu_vls_printf(gedp->ged_result_str, "%lf", state->size);
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    double size; /* must be double for scanf */

	    if (bu_sscanf(argv[2], "%lf", &size) != 1)
		goto bad;

	    state->size = (fastf_t)size;
	    (void)tclcad_data_axes_publish(view_ctx, feature_name, state);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "points")) {
	register int i;

	if (argc == 2) {
	    for (int j = 0; j < state->num_points; j++)
		bu_vls_printf(gedp->ged_result_str, " {%lf %lf %lf} ",
		    V3ARGS(state->points[j]));
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int ac;
	    const char **av;

	    if (Tcl_SplitList(interp, argv[2], &ac, &av) != TCL_OK) {
		bu_vls_printf(gedp->ged_result_str, "%s", Tcl_GetStringResult(interp));
		return BRLCAD_ERROR;
	    }

	    if (state->points) {
		bu_free(state->points, "TclCAD axes points");
		state->points = NULL;
	    }
	    state->num_points = 0;

	    if (ac < 1) {
		(void)tclcad_data_axes_publish(view_ctx, feature_name, state);
		to_refresh_view(view_ctx);
		Tcl_Free((char *)av);
		return BRLCAD_OK;
	    }

	    /* Parse new center points into temporary local array. */
	    point_t *pts = (point_t *)bu_calloc(ac, sizeof(point_t), "axes points");
	    for (i = 0; i < ac; ++i) {
		double scan[3];

		if (bu_sscanf(av[i], "%lf %lf %lf", &scan[X], &scan[Y], &scan[Z]) != 3) {
		    bu_vls_printf(gedp->ged_result_str, "bad data point - %s\n", av[i]);
		    bu_free(pts, "axes points");
		    Tcl_Free((char *)av);
		    return BRLCAD_ERROR;
		}
		VMOVE(pts[i], scan);
	    }

	    state->points = pts;
	    state->num_points = ac;
	    (void)tclcad_data_axes_publish(view_ctx, feature_name, state);
	    Tcl_Free((char *)av);
	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}
    }

bad:
    return BRLCAD_ERROR;
}

int
to_model_axes(struct ged *gedp,
	      int argc,
	      const char *argv[],
	      ged_func_ptr UNUSED(func),
	      const char *usage,
	      int UNUSED(maxargs))
{
    struct ged_view_context *view_ctx;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 3 || 6 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    view_ctx = ged_view_find_ctx(gedp, argv[1]);
    if (!view_ctx) {
        bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
        return BRLCAD_ERROR;
    }

    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    if (!view)
	return BRLCAD_ERROR;

    struct bv_axes_state axes;
    if (!bv_model_axes_state_get(&axes, view))
	return BRLCAD_ERROR;
    int ret = to_axes(gedp, view_ctx, &axes,
	"view.faceplate.model_axes.visible", argc, argv, usage);
    if (ret == BRLCAD_OK) {
	if (ged_view_context_obol_endpoint_get(view_ctx)) {
	    if (!(argc == 4 && BU_STR_EQUAL(argv[2], "draw")) && argc > 3 &&
		!tclcad_axes_endpoint_style_set(view_ctx,
		    "view.faceplate.model_axes.visible", argv[2], &axes))
		return BRLCAD_ERROR;
	} else {
	    bv_model_axes_state_set(view, &axes);
	}
    }
    return ret;
}

int
go_view_axes(struct ged *gedp,
	     struct ged_view_context *draw_view_ctx,
	     int argc,
	     const char *argv[],
	     const char *usage)
{
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 3 || 6 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    struct bv *view = bv_context_view((struct bv_context *)draw_view_ctx);
    if (!view)
	return BRLCAD_ERROR;

    struct bv_axes_state axes;
    if (!bv_view_axes_state_get(&axes, view))
	return BRLCAD_ERROR;
    int ret = to_axes(gedp, draw_view_ctx, &axes,
	"view.faceplate.view_axes.visible", argc, argv, usage);
    if (ret == BRLCAD_OK) {
	if (ged_view_context_obol_endpoint_get(draw_view_ctx)) {
	    if (!(argc == 4 && BU_STR_EQUAL(argv[2], "draw")) && argc > 3 &&
		!tclcad_axes_endpoint_style_set(draw_view_ctx,
		    "view.faceplate.view_axes.visible", argv[2], &axes))
		return BRLCAD_ERROR;
	} else {
	    bv_view_axes_state_set(view, &axes);
	}
    }
    return ret;
}


int
to_view_axes(struct ged *gedp,
	     int argc,
	     const char *argv[],
	     ged_func_ptr UNUSED(func),
	     const char *usage,
	     int UNUSED(maxargs))
{
    struct ged_view_context *view_ctx;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 3 || 6 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    view_ctx = ged_view_find_ctx(gedp, argv[1]);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    if (!view)
	return BRLCAD_ERROR;

    struct bv_axes_state axes;
    if (!bv_view_axes_state_get(&axes, view))
	return BRLCAD_ERROR;
    int ret = to_axes(gedp, view_ctx, &axes,
	"view.faceplate.view_axes.visible", argc, argv, usage);
    if (ret == BRLCAD_OK) {
	if (ged_view_context_obol_endpoint_get(view_ctx)) {
	    if (!(argc == 4 && BU_STR_EQUAL(argv[2], "draw")) && argc > 3 &&
		!tclcad_axes_endpoint_style_set(view_ctx,
		    "view.faceplate.view_axes.visible", argv[2], &axes))
		return BRLCAD_ERROR;
	} else {
	    bv_view_axes_state_set(view, &axes);
	}
    }
    return ret;
}



/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
