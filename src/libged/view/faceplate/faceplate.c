/*                      F A C E P L A T E . C
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
/** @file libged/faceplate/faceplate.c
 *
 * Controls for view elements (center dot, model axes, view axes,
 * etc.) that are built into BRL-CAD views.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "bu/cmd.h"
#include "bu/color.h"
#include "bu/opt.h"
#include "bu/vls.h"
#include "bv.h"
#include "ged/draw.h"
#include "ged/draw_obol.h"

#include "../../ged_private.h"
#include "../ged_view.h"
#include "./faceplate.h"

#ifndef COMMA
#  define COMMA ','
#endif

#define HELPFLAG "--print-help"
#define PURPOSEFLAG "--print-purpose"

struct _ged_fp_info {
    struct ged *gedp;
    long verbosity;
    const struct bu_cmdtab *cmds;
    struct bu_opt_desc *gopts;
};

static int
_fp_cmd_msgs(void *bs, int argc, const char **argv, const char *us, const char *ps)
{
    struct _ged_fp_info *gd = (struct _ged_fp_info *)bs;
    if (argc == 2 && BU_STR_EQUAL(argv[1], HELPFLAG)) {
	bu_vls_printf(gd->gedp->ged_result_str, "%s\n%s\n", us, ps);
	return 1;
    }
    if (argc == 2 && BU_STR_EQUAL(argv[1], PURPOSEFLAG)) {
	bu_vls_printf(gd->gedp->ged_result_str, "%s\n", ps);
	return 1;
    }
    return 0;
}

static int
_fp_color_property_set(struct ged *gedp, void *view_ctx,
	const char *property_name, int argc, const char **argv)
{
    struct bu_color color;
    int rgb[3] = {0, 0, 0};
    const int opt_ret = bu_opt_color(NULL, argc, argv, &color);
    if (opt_ret != 1 && opt_ret != 3) {
	bu_vls_printf(gedp->ged_result_str, "invalid color specification\n");
	return BRLCAD_ERROR;
    }
    bu_color_to_rgb_ints(&color, &rgb[0], &rgb[1], &rgb[2]);

    struct brlobol_endpoint_property_value value =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    value.type = BRLOBOL_ENDPOINT_PROPERTY_COLOR3;
    for (int i = 0; i < 3; i++)
	value.color3[i] = rgb[i] / 255.0;
    if (ged_view_context_display_property_set(view_ctx, property_name,
	&value) != BRLOBOL_ENDPOINT_PROPERTY_OK) {
	bu_vls_printf(gedp->ged_result_str,
	    "active view has no Obol faceplate color policy\n");
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

int
_fp_bool_property_set(struct ged *gedp, void *view_ctx,
	const char *property_name, int enabled)
{
    struct brlobol_endpoint_property_value value =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    value.type = BRLOBOL_ENDPOINT_PROPERTY_BOOL;
    value.bool_value = enabled;
    if (ged_view_context_display_property_set(view_ctx, property_name,
	&value) != BRLOBOL_ENDPOINT_PROPERTY_OK) {
	bu_vls_printf(gedp->ged_result_str,
	    "active view has no Obol faceplate visibility policy\n");
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

static int
_fp_bool_argument(struct ged *gedp, const char *argument, int *enabled)
{
    if (BU_STR_EQUAL("1", argument)) {
	*enabled = 1;
	return 1;
    }
    if (BU_STR_EQUAL("0", argument)) {
	*enabled = 0;
	return 1;
    }
    bu_vls_printf(gedp->ged_result_str,
	"value %s is invalid - valid values are 0 or 1\n", argument);
    return 0;
}


int
_fp_cmd_center_dot(void *ds, int argc, const char **argv)
{
    const char *usage_string = "faceplate [options] center_dot [0|1] [color r/g/b]";
    const char *purpose_string = "Enable/disable center dot and set its color.";
    if (_fp_cmd_msgs(ds, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    argc--; argv++;

    struct _ged_fp_info *gd = (struct _ged_fp_info *)ds;
    struct ged *gedp = gd->gedp;
    void *view_ctx = ged_view_active_ctx(gedp);
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    struct bv_other_state center_dot;
    if (!bv_center_dot_state_get(&center_dot, view))
	return BRLCAD_ERROR;

    if (!argc) {
	if (gd->verbosity) {
	    bu_vls_printf(gedp->ged_result_str, "%d (%d/%d/%d)", center_dot.gos_draw,
		    center_dot.gos_line_color[0], center_dot.gos_line_color[1],
		    center_dot.gos_line_color[2]);
	} else {
	    bu_vls_printf(gedp->ged_result_str, "%d", center_dot.gos_draw);
	}
	return BRLCAD_OK;
    }

    if (argc == 1) {
	int enabled = 0;
	if (!_fp_bool_argument(gedp, argv[0], &enabled))
	    return BRLCAD_ERROR;
	if (ged_view_context_display_endpoint_get(view_ctx))
	    return _fp_bool_property_set(gedp, view_ctx,
		"view.faceplate.center_dot.visible", enabled);
	center_dot.gos_draw = enabled;
	bv_center_dot_state_set(view, &center_dot);
	return BRLCAD_OK;
    }

    if (argc > 1) {
	if (!BU_STR_EQUAL("color", argv[0])) {
	    bu_vls_printf(gedp->ged_result_str, "unknown subcommand %s\n", argv[0]);
	    return BRLCAD_ERROR;
	}
	argc--; argv++;
	return _fp_color_property_set(gedp, view_ctx,
	    "view.faceplate.center_dot.color", argc, argv);
    }

    bu_vls_printf(gedp->ged_result_str, "invalid command\n");

    return BRLCAD_OK;
}

int
_fp_cmd_fb(void *ds, int argc, const char **argv)
{
    const char *usage_string = "faceplate [options] fb [0|1|2|3]";
    const char *purpose_string = "Report/set framebuffer mode";
    if (_fp_cmd_msgs(ds, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    argc--; argv++;

    struct _ged_fp_info *gd = (struct _ged_fp_info *)ds;
    struct ged *gedp = gd->gedp;
    void *view_ctx = ged_view_active_ctx(gedp);
    if (!argc) {
	struct brlobol_endpoint_property_value value =
	    BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	if (ged_view_context_display_property_get(view_ctx,
		"composition.framebuffer.mode", &value) !=
	    BRLOBOL_ENDPOINT_PROPERTY_OK) {
	    bu_vls_printf(gedp->ged_result_str,
		"active view has no Obol framebuffer composition policy\n");
	    return BRLCAD_ERROR;
	}
	if (BU_STR_EQUAL(value.string_value, "underlay"))
	    bu_vls_printf(gedp->ged_result_str, "2");
	else if (BU_STR_EQUAL(value.string_value, "interlay"))
	    bu_vls_printf(gedp->ged_result_str, "3");
	else if (BU_STR_EQUAL(value.string_value, "overlay"))
	    bu_vls_printf(gedp->ged_result_str, "1");
	else
	    bu_vls_printf(gedp->ged_result_str, "0");
	return BRLCAD_OK;
    }

    if (argc == 1) {
	const char *mode = NULL;
	if (BU_STR_EQUAL("2", argv[0]))
	    mode = "underlay";
	else if (BU_STR_EQUAL("3", argv[0]))
	    mode = "interlay";
	else if (BU_STR_EQUAL("1", argv[0]))
	    mode = "overlay";
	else if (BU_STR_EQUAL("0", argv[0]))
	    mode = "off";
	if (mode) {
	    struct brlobol_endpoint_property_value value =
		BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	    value.type = BRLOBOL_ENDPOINT_PROPERTY_ENUM;
	    value.string_value = mode;
	    if (ged_view_context_display_property_set(view_ctx,
		"composition.framebuffer.mode", &value) ==
		BRLOBOL_ENDPOINT_PROPERTY_OK)
		return BRLCAD_OK;
	    bu_vls_printf(gedp->ged_result_str,
		"active view has no Obol framebuffer composition policy\n");
	    return BRLCAD_ERROR;
	}
	bu_vls_printf(gedp->ged_result_str, "value %s is invalid - valid values are 0, 1, 2 and 3\n", argv[0]);
	return BRLCAD_ERROR;
    }

    bu_vls_printf(gedp->ged_result_str, "invalid command\n");

    return BRLCAD_OK;
}

int
_fp_cmd_scale(void *ds, int argc, const char **argv)
{
    const char *usage_string = "faceplate [options] scale [0|1] [color r/g/b]";
    const char *purpose_string = "Enable/disable scale and set its color.";
    if (_fp_cmd_msgs(ds, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    argc--; argv++;

    struct _ged_fp_info *gd = (struct _ged_fp_info *)ds;
    struct ged *gedp = gd->gedp;
    void *view_ctx = ged_view_active_ctx(gedp);
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    struct bv_other_state scale_state;
    if (!bv_scale_overlay_state_get(&scale_state, view))
	return BRLCAD_ERROR;

    if (!argc) {
	if (gd->verbosity) {
	    bu_vls_printf(gedp->ged_result_str, "%d (%d/%d/%d)", scale_state.gos_draw,
		    scale_state.gos_line_color[0], scale_state.gos_line_color[1],
		    scale_state.gos_line_color[2]);
	} else {
	    bu_vls_printf(gedp->ged_result_str, "%d", scale_state.gos_draw);
	}
	return BRLCAD_OK;
    }

    if (argc == 1) {
	int enabled = 0;
	if (!_fp_bool_argument(gedp, argv[0], &enabled))
	    return BRLCAD_ERROR;
	if (ged_view_context_display_endpoint_get(view_ctx))
	    return _fp_bool_property_set(gedp, view_ctx,
		"view.faceplate.scale.visible", enabled);
	scale_state.gos_draw = enabled;
	bv_scale_overlay_state_set(view, &scale_state);
	return BRLCAD_OK;
    }

    if (argc > 1) {
	if (!BU_STR_EQUAL("color", argv[0])) {
	    bu_vls_printf(gedp->ged_result_str, "unknown subcommand %s\n", argv[0]);
	    return BRLCAD_ERROR;
	}
	argc--; argv++;
	return _fp_color_property_set(gedp, view_ctx,
	    "view.faceplate.scale.color", argc, argv);
    }

    bu_vls_printf(gedp->ged_result_str, "invalid command\n");

    return BRLCAD_OK;
}

int
_fp_cmd_params(void *ds, int argc, const char **argv)
{
    const char *usage_string = "faceplate [options] params [0|1] [size 0/1] [center 0/1] [az 0/1] [el 0/1] [tw 0/1] [fps 0/1] [color r/g/b] [font_size 5..96]";
    const char *purpose_string = "Enable/disable params and set color/font size.";
    if (_fp_cmd_msgs(ds, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    argc--; argv++;

    struct _ged_fp_info *gd = (struct _ged_fp_info *)ds;
    struct ged *gedp = gd->gedp;
    void *view_ctx = ged_view_active_ctx(gedp);
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    struct bv_params_state params;
    if (!bv_params_state_get(&params, view))
	return BRLCAD_ERROR;

    if (!argc) {
	if (gd->verbosity) {
	    bu_vls_printf(gedp->ged_result_str, "%d[%d] (%d/%d/%d) [%d %d %d %d %d %d]",
		    params.draw,
		    params.font_size,
		    V3ARGS(params.color),
		    params.draw_size,
		    params.draw_center,
		    params.draw_az,
		    params.draw_el,
		    params.draw_tw,
		    params.draw_fps
		    );
	} else {
	    bu_vls_printf(gedp->ged_result_str, "%d", params.draw);
	}
	return BRLCAD_OK;
    }

    if (argc == 1) {
	int enabled = 0;
	if (BU_STR_EQUAL("1", argv[0]) || BU_STR_EQUAL("0", argv[0])) {
	    if (!_fp_bool_argument(gedp, argv[0], &enabled))
		return BRLCAD_ERROR;
	    if (ged_view_context_display_endpoint_get(view_ctx))
		return _fp_bool_property_set(gedp, view_ctx,
		    "view.faceplate.params.visible", enabled);
	    params.draw = enabled;
	    bv_params_state_set(view, &params);
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL("size", argv[0])) {
	    bu_vls_printf(gedp->ged_result_str, "%d", params.draw_size);
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL("center", argv[0])) {
	    bu_vls_printf(gedp->ged_result_str, "%d", params.draw_center);
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL("az", argv[0])) {
	    bu_vls_printf(gedp->ged_result_str, "%d", params.draw_az);
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL("el", argv[0])) {
	    bu_vls_printf(gedp->ged_result_str, "%d", params.draw_el);
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL("tw", argv[0])) {
	    bu_vls_printf(gedp->ged_result_str, "%d", params.draw_tw);
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL("fps", argv[0])) {
	    bu_vls_printf(gedp->ged_result_str, "%d", params.draw_fps);
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL("font_size", argv[0])) {
	    bu_vls_printf(gedp->ged_result_str, "%d", params.font_size);
	    return BRLCAD_OK;
	}
	bu_vls_printf(gedp->ged_result_str, "input %s is invalid\n", argv[0]);
	return BRLCAD_ERROR;
    }

    if (argc > 1) {
	if (BU_STR_EQUAL("color", argv[0])) {
	    argc--; argv++;
	    return _fp_color_property_set(gedp, view_ctx,
		"view.faceplate.params.color", argc, argv);
	}
	const char *property_name = NULL;
	int *draw_flag = NULL;
	if (BU_STR_EQUAL("size", argv[0])) {
	    property_name = "view.faceplate.params.size";
	    draw_flag = &params.draw_size;
	} else if (BU_STR_EQUAL("center", argv[0])) {
	    property_name = "view.faceplate.params.center";
	    draw_flag = &params.draw_center;
	} else if (BU_STR_EQUAL("az", argv[0])) {
	    property_name = "view.faceplate.params.azimuth";
	    draw_flag = &params.draw_az;
	} else if (BU_STR_EQUAL("el", argv[0])) {
	    property_name = "view.faceplate.params.elevation";
	    draw_flag = &params.draw_el;
	} else if (BU_STR_EQUAL("tw", argv[0])) {
	    property_name = "view.faceplate.params.twist";
	    draw_flag = &params.draw_tw;
	} else if (BU_STR_EQUAL("fps", argv[0])) {
	    property_name = "view.faceplate.params.fps";
	    draw_flag = &params.draw_fps;
	}
	if (property_name) {
	    int enabled = 0;
	    if (!_fp_bool_argument(gedp, argv[1], &enabled))
		return BRLCAD_ERROR;
	    if (ged_view_context_display_endpoint_get(view_ctx))
		return _fp_bool_property_set(gedp, view_ctx, property_name, enabled);
	    *draw_flag = enabled;
	    bv_params_state_set(view, &params);
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL("font_size", argv[0])) {
	    argc--; argv++;
	    int fsize = 0;
	    struct bu_vls msg = BU_VLS_INIT_ZERO;
	    if (bu_opt_int(&msg, argc, argv, &fsize) != 1 ||
		fsize < 5 || fsize > 96) {
		bu_vls_printf(gedp->ged_result_str, "invalid font size specification\n");
		bu_vls_free(&msg);
		return BRLCAD_ERROR;
	    }
	    bu_vls_free(&msg);
	    struct brlobol_endpoint_property_value value =
		BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	    value.type = BRLOBOL_ENDPOINT_PROPERTY_UINT;
	    value.uint_value = (uint64_t)fsize;
	    if (ged_view_context_display_property_set(view_ctx,
		"view.faceplate.params.font_size", &value) !=
		BRLOBOL_ENDPOINT_PROPERTY_OK) {
		bu_vls_printf(gedp->ged_result_str,
		    "active view has no Obol faceplate font policy\n");
		return BRLCAD_ERROR;
	    }
	    return BRLCAD_OK;
	}
	bu_vls_printf(gedp->ged_result_str, "unknown subcommand %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    bu_vls_printf(gedp->ged_result_str, "invalid command\n");

    return BRLCAD_OK;
}

const struct bu_cmdtab _fp_cmds[] = {
    { "center_dot",      _fp_cmd_center_dot},
    { "fb",              _fp_cmd_fb},
    { "grid",            _fp_cmd_grid},
    { "irect",           _fp_cmd_irect},
    { "model_axes",      _fp_cmd_model_axes},
    { "params",          _fp_cmd_params},
    { "scale",           _fp_cmd_scale},
    { "view_axes",       _fp_cmd_view_axes},
    { (char *)NULL,      NULL}
};

int
ged_faceplate_core(struct ged *gedp, int argc, const char *argv[])
{
    int help = 0;
    struct _ged_fp_info gd;
    gd.gedp = gedp;
    gd.cmds = _fp_cmds;
    gd.verbosity = 0;

    // Sanity
    if (UNLIKELY(!gedp || !argc || !argv)) {
	return BRLCAD_ERROR;
    }

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    // We know we're the faceplate command - start processing args
    argc--; argv++;

    // See if we have any high level options set
    struct bu_opt_desc d[3];
    BU_OPT(d[0], "h", "help",    "",  NULL,               &help,         "Print help");
    BU_OPT(d[1], "v", "verbose", "",  &bu_opt_incr_long,  &gd.verbosity, "Verbose output");
    BU_OPT_NULL(d[2]);

    gd.gopts = d;

    int ac = bu_opt_parse(NULL, argc, argv, d);

    if (!ac || help) {
	_ged_subcmd_help(gedp, (struct bu_opt_desc *)d, (const struct bu_cmdtab *)_fp_cmds,
	       	"view faceplate", "[options] subcommand [args]", &gd, 0, NULL);
	return BRLCAD_OK;
    }

    void *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, ": no current view set");
	return BRLCAD_ERROR;
    }

    int ret;
    if (bu_cmd(_fp_cmds, ac, argv, 0, (void *)&gd, &ret) == BRLCAD_OK) {
	if (ret == BRLCAD_OK)
	    (void)ged_draw_obol_view_context_faceplate_sync(gedp, view_ctx);
	return ret;
    } else {
	bu_vls_printf(gedp->ged_result_str, "subcommand %s not defined", argv[0]);
    }

    return BRLCAD_ERROR;
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
