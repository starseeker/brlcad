/*                        A X E S . C
 * BRL-CAD
 *
 * Copyright (c) 2008-2026 United States Government as represented by
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
/** @file libged/view/faceplate/faceplate_axes.cpp
 *
 * Commands for HUD axes
 *
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "bu/cmd.h"
#include "bu/color.h"
#include "bu/opt.h"
#include "bu/vls.h"
#include "bv.h"
#include "BObol/BDisplayEndpoint.h"

#include "../../ged_private.h"
#include "../ged_view.h"
#include "./faceplate.h"

struct _ged_fp_axes_info {
    struct _ged_view_info *gd;
    struct bv_axes_state *a;
};

struct _fp_axes_property_names {
    const char *position[3];
    const char *size;
    const char *line_width;
    const char *color;
    const char *position_only;
    const char *labels_visible;
    const char *labels_color;
    const char *triple_color;
    const char *ticks_visible;
    const char *ticks_length;
    const char *ticks_major_length;
    const char *ticks_interval;
    const char *ticks_per_major;
    const char *ticks_threshold;
    const char *ticks_color;
    const char *ticks_major_color;
};

static const struct _fp_axes_property_names model_axes_properties = {
    {"view.faceplate.model_axes.position.x",
	"view.faceplate.model_axes.position.y",
	"view.faceplate.model_axes.position.z"},
    "view.faceplate.model_axes.size",
    "view.faceplate.model_axes.line_width",
    "view.faceplate.model_axes.color",
    "view.faceplate.model_axes.position_only",
    "view.faceplate.model_axes.labels.visible",
    "view.faceplate.model_axes.labels.color",
    "view.faceplate.model_axes.triple_color",
    "view.faceplate.model_axes.ticks.visible",
    "view.faceplate.model_axes.ticks.length",
    "view.faceplate.model_axes.ticks.major_length",
    "view.faceplate.model_axes.ticks.interval",
    "view.faceplate.model_axes.ticks.per_major",
    "view.faceplate.model_axes.ticks.threshold",
    "view.faceplate.model_axes.ticks.color",
    "view.faceplate.model_axes.ticks.major_color"
};

static const struct _fp_axes_property_names view_axes_properties = {
    {"view.faceplate.view_axes.position.x",
	"view.faceplate.view_axes.position.y",
	"view.faceplate.view_axes.position.z"},
    "view.faceplate.view_axes.size",
    "view.faceplate.view_axes.line_width",
    "view.faceplate.view_axes.color",
    "view.faceplate.view_axes.position_only",
    "view.faceplate.view_axes.labels.visible",
    "view.faceplate.view_axes.labels.color",
    "view.faceplate.view_axes.triple_color",
    "view.faceplate.view_axes.ticks.visible",
    "view.faceplate.view_axes.ticks.length",
    "view.faceplate.view_axes.ticks.major_length",
    "view.faceplate.view_axes.ticks.interval",
    "view.faceplate.view_axes.ticks.per_major",
    "view.faceplate.view_axes.ticks.threshold",
    "view.faceplate.view_axes.ticks.color",
    "view.faceplate.view_axes.ticks.major_color"
};

static int
_fp_axes_endpoint_property_set(struct ged_view_context *view_ctx, const char *name,
	const struct bv_display_property_value *value)
{
    bobol_display_endpoint_t *endpoint =
	ged_view_context_obol_endpoint_get(view_ctx);
    return endpoint && bobol_display_endpoint_property_set(endpoint, name,
	value) == BV_DISPLAY_PROPERTY_OK;
}

static int
_fp_axes_endpoint_bool_set(struct ged_view_context *view_ctx, const char *name, int enabled)
{
    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    value.type = BV_DISPLAY_PROPERTY_BOOL;
    value.bool_value = enabled ? 1 : 0;
    return _fp_axes_endpoint_property_set(view_ctx, name, &value);
}

static int
_fp_axes_endpoint_double_set(struct ged_view_context *view_ctx, const char *name, double number)
{
    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    value.type = BV_DISPLAY_PROPERTY_DOUBLE;
    value.double_value = number;
    return _fp_axes_endpoint_property_set(view_ctx, name, &value);
}

static int
_fp_axes_endpoint_uint_set(struct ged_view_context *view_ctx, const char *name, int number)
{
    if (number < 0)
	return 0;
    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    value.type = BV_DISPLAY_PROPERTY_UINT;
    value.uint_value = (uint64_t)number;
    return _fp_axes_endpoint_property_set(view_ctx, name, &value);
}

static int
_fp_axes_endpoint_color_set(struct ged_view_context *view_ctx, const char *name,
	const int color[3])
{
    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    value.type = BV_DISPLAY_PROPERTY_COLOR3;
    for (int i = 0; i < 3; i++)
	value.color3[i] = color[i] / 255.0;
    return _fp_axes_endpoint_property_set(view_ctx, name, &value);
}

static int
_fp_axes_endpoint_state_apply(struct ged_view_context *view_ctx,
	const struct _fp_axes_property_names *properties,
	const struct bv_axes_state *before, const struct bv_axes_state *after)
{
    if (!properties || !before || !after)
	return 0;
    for (int axis = 0; axis < 3; axis++) {
	if (!NEAR_EQUAL(before->axes_pos[axis], after->axes_pos[axis],
		SMALL_FASTF) &&
	    !_fp_axes_endpoint_double_set(view_ctx, properties->position[axis],
		after->axes_pos[axis]))
	    return 0;
    }
    if (!NEAR_EQUAL(before->axes_size, after->axes_size, SMALL_FASTF) &&
	!_fp_axes_endpoint_double_set(view_ctx, properties->size,
	    after->axes_size))
	return 0;
    if (before->line_width != after->line_width &&
	!_fp_axes_endpoint_uint_set(view_ctx, properties->line_width,
	    after->line_width))
	return 0;
    if (memcmp(before->axes_color, after->axes_color,
	sizeof(after->axes_color)) != 0 &&
	!_fp_axes_endpoint_color_set(view_ctx, properties->color,
	    after->axes_color))
	return 0;
    if (before->pos_only != after->pos_only &&
	!_fp_axes_endpoint_bool_set(view_ctx, properties->position_only,
	    after->pos_only))
	return 0;
    if (before->label_flag != after->label_flag &&
	!_fp_axes_endpoint_bool_set(view_ctx, properties->labels_visible,
	    after->label_flag))
	return 0;
    if (memcmp(before->label_color, after->label_color,
	sizeof(after->label_color)) != 0 &&
	!_fp_axes_endpoint_color_set(view_ctx, properties->labels_color,
	    after->label_color))
	return 0;
    if (before->triple_color != after->triple_color &&
	!_fp_axes_endpoint_bool_set(view_ctx, properties->triple_color,
	    after->triple_color))
	return 0;
    if (before->tick_enabled != after->tick_enabled &&
	!_fp_axes_endpoint_bool_set(view_ctx, properties->ticks_visible,
	    after->tick_enabled))
	return 0;
    if (before->tick_length != after->tick_length &&
	!_fp_axes_endpoint_uint_set(view_ctx, properties->ticks_length,
	    after->tick_length))
	return 0;
    if (before->tick_major_length != after->tick_major_length &&
	!_fp_axes_endpoint_uint_set(view_ctx, properties->ticks_major_length,
	    after->tick_major_length))
	return 0;
    if (!NEAR_EQUAL(before->tick_interval, after->tick_interval,
	SMALL_FASTF) && !_fp_axes_endpoint_double_set(view_ctx,
	properties->ticks_interval, after->tick_interval))
	return 0;
    if (before->ticks_per_major != after->ticks_per_major &&
	!_fp_axes_endpoint_uint_set(view_ctx, properties->ticks_per_major,
	    after->ticks_per_major))
	return 0;
    if (before->tick_threshold != after->tick_threshold &&
	!_fp_axes_endpoint_uint_set(view_ctx, properties->ticks_threshold,
	    after->tick_threshold))
	return 0;
    if (memcmp(before->tick_color, after->tick_color,
	sizeof(after->tick_color)) != 0 &&
	!_fp_axes_endpoint_color_set(view_ctx, properties->ticks_color,
	    after->tick_color))
	return 0;
    if (memcmp(before->tick_major_color, after->tick_major_color,
	sizeof(after->tick_major_color)) != 0 &&
	!_fp_axes_endpoint_color_set(view_ctx, properties->ticks_major_color,
	    after->tick_major_color))
	return 0;
    return 1;
}

int
_fp_axes_cmd_size(void *bs, int argc, const char **argv)
{
    struct _ged_fp_axes_info *ainfo = (struct _ged_fp_axes_info *)bs;
    struct _ged_view_info *gd = ainfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_axes size [#]";
    const char *purpose_string = "adjust axes size";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_axes_state *a = ainfo->a;
     if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%f\n", a->axes_size);
	return BRLCAD_OK;
    }

    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    fastf_t val;
    if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[0], (void *)&val) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    a->axes_size = val;

    return BRLCAD_OK;
}

int
_fp_axes_cmd_linewidth(void *bs, int argc, const char **argv)
{

    struct _ged_fp_axes_info *ainfo = (struct _ged_fp_axes_info *)bs;
    struct _ged_view_info *gd = ainfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_axes linewidth [#]";
    const char *purpose_string = "adjust axes line width";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_axes_state *a = ainfo->a;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", a->line_width);
	return BRLCAD_OK;
    }

    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int val;
    if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&val) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    if (val < 1) {
	bu_vls_printf(gedp->ged_result_str, "Smallest supported value is 1\n");
	return BRLCAD_ERROR;
    }

    a->line_width = val;

    return BRLCAD_OK;
}

int
_fp_axes_cmd_pos_only(void *bs, int argc, const char **argv)
{

    struct _ged_fp_axes_info *ainfo = (struct _ged_fp_axes_info *)bs;
    struct _ged_view_info *gd = ainfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_axes pos_only [0|1]";
    const char *purpose_string = "enable/disable axes decorations";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_axes_state *a = ainfo->a;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", a->pos_only);
	return BRLCAD_OK;
    }

     if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    int val;
    if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&val) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    val = (val) ? 1 : 0;

    a->pos_only = val;

    return BRLCAD_OK;
}

int
_fp_axes_cmd_fp_axes_color(void *bs, int argc, const char **argv)
{

    struct _ged_fp_axes_info *ainfo = (struct _ged_fp_axes_info *)bs;
    struct _ged_view_info *gd = ainfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_axes color [r/g/b]";
    const char *purpose_string = "get/set color of axes";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_axes_state *a = ainfo->a;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d %d %d\n", a->axes_color[0], a->axes_color[1], a->axes_color[2]);
	return BRLCAD_OK;
    }

    // For color need either 1 or 3 non-subcommand args
    if (argc != 1 && argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    struct bu_color c;
    int opt_ret = bu_opt_color(NULL, argc, (const char **)argv, (void *)&c);
    if (opt_ret != 1 && opt_ret != 3) {
	bu_vls_printf(gedp->ged_result_str, "Invalid color specifier\n");
	return BRLCAD_ERROR;
    }

    bu_color_to_rgb_ints(&c, &a->axes_color[0], &a->axes_color[1], &a->axes_color[2]);

    return BRLCAD_OK;
}

int
_fp_axes_cmd_label(void *bs, int argc, const char **argv)
{

    struct _ged_fp_axes_info *ainfo = (struct _ged_fp_axes_info *)bs;
    struct _ged_view_info *gd = ainfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_axes label [0|1]";
    const char *purpose_string = "enable/disable text labels for axes";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_axes_state *a = ainfo->a;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", a->label_flag);
	return BRLCAD_OK;
    }

    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int val;
    if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&val) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    val = (val) ? 1 : 0;

    a->label_flag = val;

    return BRLCAD_OK;
}

int
_fp_axes_cmd_label_color(void *bs, int argc, const char **argv)
{

    struct _ged_fp_axes_info *ainfo = (struct _ged_fp_axes_info *)bs;
    struct _ged_view_info *gd = ainfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_axes label_color [r/g/b]";
    const char *purpose_string = "get/set color of text labels for axes";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_axes_state *a = ainfo->a;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d %d %d\n", a->label_color[0], a->label_color[1], a->label_color[2]);
	return BRLCAD_OK;
    }

    // For color need either 1 or 3 non-subcommand args
    if (argc != 1 && argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    struct bu_color c;
    int opt_ret = bu_opt_color(NULL, argc, (const char **)argv, (void *)&c);
    if (opt_ret != 1 && opt_ret != 3) {
	bu_vls_printf(gedp->ged_result_str, "Invalid color specifier\n");
	return BRLCAD_ERROR;
    }

    bu_color_to_rgb_ints(&c, &a->label_color[0], &a->label_color[1], &a->label_color[2]);


    return BRLCAD_OK;
}

int
_fp_axes_cmd_triple_color(void *bs, int argc, const char **argv)
{

    struct _ged_fp_axes_info *ainfo = (struct _ged_fp_axes_info *)bs;
    struct _ged_view_info *gd = ainfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_axes triple_color [0|1]";
    const char *purpose_string = "enable/disable tri-color mode for axes coloring";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_axes_state *a = ainfo->a;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", a->triple_color);
	return BRLCAD_OK;
    }

    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int val;
    if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&val) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    val = (val) ? 1 : 0;

    a->triple_color = val;

    return BRLCAD_OK;
}

int
_fp_axes_cmd_tick(void *bs, int argc, const char **argv)
{

    struct _ged_fp_axes_info *ainfo = (struct _ged_fp_axes_info *)bs;
    struct _ged_view_info *gd = ainfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_axes tick [0|1]";
    const char *purpose_string = "enable/disable axes tick drawing";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_axes_state *a = ainfo->a;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", a->tick_enabled);
	return BRLCAD_OK;
    }

    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int val;
    if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&val) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    val = (val) ? 1 : 0;

    a->tick_enabled = val;

    return BRLCAD_OK;
}

int
_fp_axes_cmd_tick_length(void *bs, int argc, const char **argv)
{

    struct _ged_fp_axes_info *ainfo = (struct _ged_fp_axes_info *)bs;
    struct _ged_view_info *gd = ainfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_axes tick_length [#]";
    const char *purpose_string = "get/set tick length";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_axes_state *a = ainfo->a;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", a->tick_length);
	return BRLCAD_OK;
    }

    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int val;
    if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&val) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    if (val < 1) {
	bu_vls_printf(gedp->ged_result_str, "Smallest supported value is 1\n");
	return BRLCAD_ERROR;
    }

    a->tick_length = val;

    return BRLCAD_OK;
}

int
_fp_axes_cmd_tick_major_length(void *bs, int argc, const char **argv)
{

    struct _ged_fp_axes_info *ainfo = (struct _ged_fp_axes_info *)bs;
    struct _ged_view_info *gd = ainfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_axes tick_major_length [#]";
    const char *purpose_string = "get/set tick major length";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_axes_state *a = ainfo->a;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", a->tick_major_length);
	return BRLCAD_OK;
    }

    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int val;
    if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&val) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    if (val < 1) {
	bu_vls_printf(gedp->ged_result_str, "Smallest supported value is 1\n");
	return BRLCAD_ERROR;
    }

    a->tick_major_length = val;

    return BRLCAD_OK;
}

int
_fp_axes_cmd_tick_interval(void *bs, int argc, const char **argv)
{

    struct _ged_fp_axes_info *ainfo = (struct _ged_fp_axes_info *)bs;
    struct _ged_view_info *gd = ainfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_axes tick_interval [#]";
    const char *purpose_string = "get/set tick interval";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_axes_state *a = ainfo->a;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%f\n", a->tick_interval);
	return BRLCAD_OK;
    }

    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    fastf_t val;
    if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[0], (void *)&val) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    a->tick_interval = val;

    return BRLCAD_OK;
}

int
_fp_axes_cmd_ticks_per_major(void *bs, int argc, const char **argv)
{

    struct _ged_fp_axes_info *ainfo = (struct _ged_fp_axes_info *)bs;
    struct _ged_view_info *gd = ainfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_axes ticks_per_major [#]";
    const char *purpose_string = "get/set ticks per major";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_axes_state *a = ainfo->a;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", a->ticks_per_major);
	return BRLCAD_OK;
    }

    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int val;
    if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&val) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    if (val < 0) {
	bu_vls_printf(gedp->ged_result_str, "Smallest supported value is 0\n");
	return BRLCAD_ERROR;
    }

    a->ticks_per_major = val;

    return BRLCAD_OK;
}

int
_fp_axes_cmd_tick_threshold(void *bs, int argc, const char **argv)
{

    struct _ged_fp_axes_info *ainfo = (struct _ged_fp_axes_info *)bs;
    struct _ged_view_info *gd = ainfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_axes tick_threshold [#]";
    const char *purpose_string = "get/set tick threshold";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_axes_state *a = ainfo->a;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", a->tick_threshold);
	return BRLCAD_OK;
    }

    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int val;
    if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&val) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    if (val < 0) {
	bu_vls_printf(gedp->ged_result_str, "Smallest supported value is 0\n");
	return BRLCAD_ERROR;
    }

    a->tick_threshold = val;

    return BRLCAD_OK;
}

int
_fp_axes_cmd_tick_color(void *bs, int argc, const char **argv)
{

    struct _ged_fp_axes_info *ainfo = (struct _ged_fp_axes_info *)bs;
    struct _ged_view_info *gd = ainfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_axes tick_color [r/g/b]";
    const char *purpose_string = "get/set color of ticks";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_axes_state *a = ainfo->a;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d %d %d\n", a->tick_color[0], a->tick_color[1], a->tick_color[2]);
	return BRLCAD_OK;
    }

    // For color need either 1 or 3 non-subcommand args
    if (argc != 1 && argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    struct bu_color c;
    int opt_ret = bu_opt_color(NULL, argc, (const char **)argv, (void *)&c);
    if (opt_ret != 1 && opt_ret != 3) {
	bu_vls_printf(gedp->ged_result_str, "Invalid color specifier\n");
	return BRLCAD_ERROR;
    }

    bu_color_to_rgb_ints(&c, &a->tick_color[0], &a->tick_color[1], &a->tick_color[2]);


    return BRLCAD_OK;
}

int
_fp_axes_cmd_tick_major_color(void *bs, int argc, const char **argv)
{

    struct _ged_fp_axes_info *ainfo = (struct _ged_fp_axes_info *)bs;
    struct _ged_view_info *gd = ainfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_axes tick_major_color [r/g/b]";
    const char *purpose_string = "get/set tick_major_color";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_axes_state *a = ainfo->a;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d %d %d\n", a->tick_major_color[0], a->tick_major_color[1], a->tick_major_color[2]);
	return BRLCAD_OK;
    }

    // For color need either 1 or 3 non-subcommand args
    if (argc != 1 && argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    struct bu_color c;
    int opt_ret = bu_opt_color(NULL, argc, (const char **)argv, (void *)&c);
    if (opt_ret != 1 && opt_ret != 3) {
	bu_vls_printf(gedp->ged_result_str, "Invalid color specifier\n");
	return BRLCAD_ERROR;
    }

    bu_color_to_rgb_ints(&c, &a->tick_major_color[0], &a->tick_major_color[1], &a->tick_major_color[2]);

    return BRLCAD_OK;
}

const struct bu_cmdtab _fp_axes_cmds[] = {
    { "size",              _fp_axes_cmd_size},
    { "line_width",        _fp_axes_cmd_linewidth},
    { "pos_only",          _fp_axes_cmd_pos_only},
    { "axes_color",        _fp_axes_cmd_fp_axes_color},
    { "label",             _fp_axes_cmd_label},
    { "label_color",       _fp_axes_cmd_label_color},
    { "triple_color",      _fp_axes_cmd_triple_color},
    { "tick",              _fp_axes_cmd_tick},
    { "tick_length",       _fp_axes_cmd_tick_length},
    { "tick_major_length", _fp_axes_cmd_tick_major_length},
    { "tick_interval",     _fp_axes_cmd_tick_interval},
    { "ticks_per_major",   _fp_axes_cmd_ticks_per_major},
    { "tick_threshold",    _fp_axes_cmd_tick_threshold},
    { "tick_color",        _fp_axes_cmd_tick_color},
    { "tick_major_color",  _fp_axes_cmd_tick_major_color},
    { (char *)NULL,      NULL}
};

int
_fp_cmd_model_axes(void *bs, int argc, const char **argv)
{
    int help = 0;
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);

    const char *usage_string = "view faceplate model_axes subcmd [args]";
    const char *purpose_string = "manipulate view axes";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, ": no view current in GED");
	return BRLCAD_ERROR;
    }
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);

    // We know we're the axes command - start processing args
    argc--; argv++;

    if (argc == 1) {
	struct bv_axes_state axes;
	if (!bv_model_axes_state_get(&axes, view))
	    return BRLCAD_ERROR;
	if (BU_STR_EQUAL("1", argv[0])) {
	    if (ged_view_context_obol_endpoint_get(view_ctx))
		return _fp_bool_property_set(gedp, view_ctx,
		    "view.faceplate.model_axes.visible", 1);
	    axes.draw = 1;
	    bv_model_axes_state_set(view, &axes);
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL("0", argv[0])) {
	    if (ged_view_context_obol_endpoint_get(view_ctx))
		return _fp_bool_property_set(gedp, view_ctx,
		    "view.faceplate.model_axes.visible", 0);
	    axes.draw = 0;
	    bv_model_axes_state_set(view, &axes);
	    return BRLCAD_OK;
	}
    }

    // See if we have any high level options set
    struct bu_opt_desc d[2];
    BU_OPT(d[0], "h", "help",  "",  NULL,  &help,      "Print help");
    BU_OPT_NULL(d[1]);

    gd->gopts = d;

    // High level options are only defined prior to the subcommand
    int cmd_pos = -1;
    for (int i = 0; i < argc; i++) {
	if (bu_cmd_valid(_fp_axes_cmds, argv[i]) == BRLCAD_OK) {
	    cmd_pos = i;
	    break;
	}
    }

    int acnt = (cmd_pos >= 0) ? cmd_pos : argc;
    (void)bu_opt_parse(NULL, acnt, argv, d);

    struct bv_axes_state axes;
    if (!bv_model_axes_state_get(&axes, view))
	return BRLCAD_ERROR;
    const struct bv_axes_state initial_axes = axes;

    struct _ged_fp_axes_info ainfo;
    ainfo.gd = gd;
    ainfo.a = &axes;

    int ret = _ged_subcmd_exec(gedp, d, _fp_axes_cmds, "view faceplate model_axes", "[options] subcommand [args]", (void *)&ainfo, argc, argv, help, cmd_pos);
    if (ret == BRLCAD_OK) {
	if (ged_view_context_obol_endpoint_get(view_ctx)) {
	    if (!_fp_axes_endpoint_state_apply(view_ctx, &model_axes_properties,
		&initial_axes, &axes))
		return BRLCAD_ERROR;
	} else {
	    bv_model_axes_state_set(view, &axes);
	}
    }
    return ret;
}

int
_fp_cmd_view_axes(void *bs, int argc, const char **argv)
{
    int help = 0;
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);

    const char *usage_string = "view faceplate view_axes subcmd [args]";
    const char *purpose_string = "manipulate view axes";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, ": no view current in GED");
	return BRLCAD_ERROR;
    }
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);


    // We know we're the axes command - start processing args
    argc--; argv++;

    if (argc == 1) {
	struct bv_axes_state axes;
	if (!bv_view_axes_state_get(&axes, view))
	    return BRLCAD_ERROR;
	if (BU_STR_EQUAL("1", argv[0])) {
	    if (ged_view_context_obol_endpoint_get(view_ctx))
		return _fp_bool_property_set(gedp, view_ctx,
		    "view.faceplate.view_axes.visible", 1);
	    axes.draw = 1;
	    bv_view_axes_state_set(view, &axes);
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL("0", argv[0])) {
	    if (ged_view_context_obol_endpoint_get(view_ctx))
		return _fp_bool_property_set(gedp, view_ctx,
		    "view.faceplate.view_axes.visible", 0);
	    axes.draw = 0;
	    bv_view_axes_state_set(view, &axes);
	    return BRLCAD_OK;
	}
    }

    // See if we have any high level options set
    struct bu_opt_desc d[2];
    BU_OPT(d[0], "h", "help",  "",  NULL,  &help,      "Print help");
    BU_OPT_NULL(d[1]);

    gd->gopts = d;

    // High level options are only defined prior to the subcommand
    int cmd_pos = -1;
    for (int i = 0; i < argc; i++) {
	if (bu_cmd_valid(_fp_axes_cmds, argv[i]) == BRLCAD_OK) {
	    cmd_pos = i;
	    break;
	}
    }

    int acnt = (cmd_pos >= 0) ? cmd_pos : argc;
    (void)bu_opt_parse(NULL, acnt, argv, d);

    struct bv_axes_state axes;
    if (!bv_view_axes_state_get(&axes, view))
	return BRLCAD_ERROR;
    const struct bv_axes_state initial_axes = axes;

    struct _ged_fp_axes_info ainfo;
    ainfo.gd = gd;
    ainfo.a = &axes;

    int ret = _ged_subcmd_exec(gedp, d, _fp_axes_cmds, "view faceplate view_axes", "[options] subcommand [args]", (void *)&ainfo, argc, argv, help, cmd_pos);
    if (ret == BRLCAD_OK) {
	if (ged_view_context_obol_endpoint_get(view_ctx)) {
	    if (!_fp_axes_endpoint_state_apply(view_ctx, &view_axes_properties,
		&initial_axes, &axes))
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
