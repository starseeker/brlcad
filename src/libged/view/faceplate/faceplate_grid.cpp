/*                        G R I D . C
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
/** @file libged/view/faceplate/faceplate_grid.cpp
 *
 * Commands for HUD grid overlay
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

struct _ged_fp_grid_info {
    struct _ged_view_info *gd;
    struct bv_grid_state *g;
    int draw_property_set;
};


static int
_fp_grid_endpoint_property_set(struct ged_view_context *view_ctx, const char *name,
	const struct bv_display_property_value *value)
{
    bobol_display_endpoint_t *endpoint =
	ged_view_context_obol_endpoint_get(view_ctx);
    return endpoint && bobol_display_endpoint_property_set(endpoint, name,
	value) == BV_DISPLAY_PROPERTY_OK;
}

static int
_fp_grid_endpoint_bool_set(struct ged_view_context *view_ctx, const char *name, int enabled)
{
    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    value.type = BV_DISPLAY_PROPERTY_BOOL;
    value.bool_value = enabled ? 1 : 0;
    return _fp_grid_endpoint_property_set(view_ctx, name, &value);
}

static int
_fp_grid_endpoint_double_set(struct ged_view_context *view_ctx, const char *name, double number)
{
    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    value.type = BV_DISPLAY_PROPERTY_DOUBLE;
    value.double_value = number;
    return _fp_grid_endpoint_property_set(view_ctx, name, &value);
}

static int
_fp_grid_endpoint_uint_set(struct ged_view_context *view_ctx, const char *name, int number)
{
    if (number < 0)
	return 0;
    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    value.type = BV_DISPLAY_PROPERTY_UINT;
    value.uint_value = (uint64_t)number;
    return _fp_grid_endpoint_property_set(view_ctx, name, &value);
}

static int
_fp_grid_endpoint_state_apply(struct ged_view_context *view_ctx,
	const struct bv_grid_state *before, const struct bv_grid_state *after)
{
    if (!before || !after)
	return 0;
    if (before->adaptive != after->adaptive &&
	!_fp_grid_endpoint_bool_set(view_ctx, "view.faceplate.grid.adaptive",
	    after->adaptive))
	return 0;
    if (before->snap != after->snap &&
	!_fp_grid_endpoint_bool_set(view_ctx, "view.faceplate.grid.snap",
	    after->snap))
	return 0;

    const char *anchor_properties[] = {
	"view.faceplate.grid.anchor.x",
	"view.faceplate.grid.anchor.y",
	"view.faceplate.grid.anchor.z"
    };

    for (int axis = 0; axis < 3; axis++) {
	if (!NEAR_EQUAL(before->anchor[axis], after->anchor[axis],
		SMALL_FASTF) &&
	    !_fp_grid_endpoint_double_set(view_ctx, anchor_properties[axis],
		after->anchor[axis]))
	    return 0;
    }

    if (!NEAR_EQUAL(before->res_h, after->res_h, SMALL_FASTF) &&
	!_fp_grid_endpoint_double_set(view_ctx,
	    "view.faceplate.grid.resolution.horizontal", after->res_h))
	return 0;

    if (!NEAR_EQUAL(before->res_v, after->res_v, SMALL_FASTF) &&
	!_fp_grid_endpoint_double_set(view_ctx,
	    "view.faceplate.grid.resolution.vertical", after->res_v))
	return 0;
    if (before->res_major_h != after->res_major_h &&
	!_fp_grid_endpoint_uint_set(view_ctx,
	    "view.faceplate.grid.major.horizontal", after->res_major_h))
	return 0;
    if (before->res_major_v != after->res_major_v &&
	!_fp_grid_endpoint_uint_set(view_ctx,
	    "view.faceplate.grid.major.vertical", after->res_major_v))
	return 0;
    if (memcmp(before->color, after->color, sizeof(after->color)) != 0) {
	struct bv_display_property_value value =
	    BV_DISPLAY_PROPERTY_VALUE_INIT;
	value.type = BV_DISPLAY_PROPERTY_COLOR3;
	for (int axis = 0; axis < 3; axis++)
	    value.color3[axis] = after->color[axis] / 255.0;
	if (!_fp_grid_endpoint_property_set(view_ctx,
		"view.faceplate.grid.color", &value))
	    return 0;
    }
    return 1;
}

int
_fp_grid_cmd_draw(void *bs, int argc, const char **argv)
{

    struct _ged_fp_grid_info *ginfo = (struct _ged_fp_grid_info *)bs;
    struct _ged_view_info *gd = ginfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate grid draw [0|1]";
    const char *purpose_string = "enable/disable grid drawing";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_grid_state *g = ginfo->g;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", g->draw);
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
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    if (ged_view_context_obol_endpoint_get(view_ctx)) {
	if (_fp_bool_property_set(gedp, view_ctx,
		"view.faceplate.grid.visible", val) != BRLCAD_OK)
	    return BRLCAD_ERROR;
	ginfo->draw_property_set = 1;
	return BRLCAD_OK;
    }
    g->draw = val;

    return BRLCAD_OK;
}

int
_fp_grid_cmd_snap(void *bs, int argc, const char **argv)
{

    struct _ged_fp_grid_info *ginfo = (struct _ged_fp_grid_info *)bs;
    struct _ged_view_info *gd = ginfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate grid snap [0|1]";
    const char *purpose_string = "enable/disable grid snapping";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_grid_state *g = ginfo->g;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", g->snap);
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

    g->snap = val;

    return BRLCAD_OK;
}

int
_fp_grid_cmd_anchor(void *bs, int argc, const char **argv)
{
    struct _ged_fp_grid_info *ginfo = (struct _ged_fp_grid_info *)bs;
    struct _ged_view_info *gd = ginfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate grid anchor [# # #]";
    const char *purpose_string = "adjust grid size";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_grid_state *g = ginfo->g;
     if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%g %g %g\n", V3ARGS(g->anchor));
	return BRLCAD_OK;
    }

    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    vect_t val;
    int ret = bu_opt_vect_t(NULL, argc, (const char **)&argv[0], (void *)&val);
    if (ret != 1 && ret != 3) {
	bu_vls_printf(gedp->ged_result_str, "Invalid specification\n");
	return BRLCAD_ERROR;
    }

    VMOVE(g->anchor, val);

    return BRLCAD_OK;
}

int
_fp_grid_cmd_res_h(void *bs, int argc, const char **argv)
{

    struct _ged_fp_grid_info *ginfo = (struct _ged_fp_grid_info *)bs;
    struct _ged_view_info *gd = ginfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate grid res_ [#]";
    const char *purpose_string = "adjust grid line width";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_grid_state *g = ginfo->g;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%g\n", g->res_h);
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

    if (val < BN_TOL_DIST) {
	bu_vls_printf(gedp->ged_result_str, "Smallest supported value is %f\n", BN_TOL_DIST);
	return BRLCAD_ERROR;
    }

    g->res_h = val;

    return BRLCAD_OK;
}

int
_fp_grid_cmd_res_v(void *bs, int argc, const char **argv)
{

    struct _ged_fp_grid_info *ginfo = (struct _ged_fp_grid_info *)bs;
    struct _ged_view_info *gd = ginfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate grid res_ [#]";
    const char *purpose_string = "adjust grid line width";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_grid_state *g = ginfo->g;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%g\n", g->res_v);
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

    if (val < BN_TOL_DIST) {
	bu_vls_printf(gedp->ged_result_str, "Smallest supported value is %f\n", BN_TOL_DIST);
	return BRLCAD_ERROR;
    }

    g->res_v = val;

    return BRLCAD_OK;
}

int
_fp_grid_cmd_res_major_h(void *bs, int argc, const char **argv)
{

    struct _ged_fp_grid_info *ginfo = (struct _ged_fp_grid_info *)bs;
    struct _ged_view_info *gd = ginfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate grid res_ [#]";
    const char *purpose_string = "adjust grid line width";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_grid_state *g = ginfo->g;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", g->res_major_h);
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

    g->res_major_h = val;

    return BRLCAD_OK;
}

int
_fp_grid_cmd_res_major_v(void *bs, int argc, const char **argv)
{

    struct _ged_fp_grid_info *ginfo = (struct _ged_fp_grid_info *)bs;
    struct _ged_view_info *gd = ginfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate grid res_major_v [#]";
    const char *purpose_string = "adjust grid line width";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_grid_state *g = ginfo->g;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", g->res_major_v);
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

    g->res_major_v = val;

    return BRLCAD_OK;
}

int
_fp_grid_cmd_color(void *bs, int argc, const char **argv)
{

    struct _ged_fp_grid_info *ginfo = (struct _ged_fp_grid_info *)bs;
    struct _ged_view_info *gd = ginfo->gd;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view faceplate [model|view]_grid color [r/g/b]";
    const char *purpose_string = "get/set color of grid";
    if (_view_cmd_msgs((void *)gd, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bv_grid_state *g = ginfo->g;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d %d %d\n", g->color[0], g->color[1], g->color[2]);
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

    bu_color_to_rgb_ints(&c, &g->color[0], &g->color[1], &g->color[2]);

    return BRLCAD_OK;
}

const struct bu_cmdtab _fp_grid_cmds[] = {
    { "draw",        _fp_grid_cmd_draw},
    { "snap",        _fp_grid_cmd_snap},
    { "anchor",      _fp_grid_cmd_anchor},
    { "res_h",       _fp_grid_cmd_res_h},
    { "res_v",       _fp_grid_cmd_res_v},
    { "res_major_h", _fp_grid_cmd_res_major_h},
    { "res_major_v", _fp_grid_cmd_res_major_v},
    { "color",       _fp_grid_cmd_color},
    { (char *)NULL,      NULL}
};

int
_fp_cmd_grid(void *bs, int argc, const char **argv)
{
    int help = 0;
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);

    const char *usage_string = "view faceplate grid subcmd [args]";
    const char *purpose_string = "manipulate faceplate grid overlay";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, ": no view current in GED");
	return BRLCAD_ERROR;
    }
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);

    // We know we're the grid command - start processing args
    argc--; argv++;

    if (argc == 1) {
	struct bv_grid_state grid;
	if (!bv_grid_state_get(&grid, view))
	    return BRLCAD_ERROR;
	if (BU_STR_EQUAL("1", argv[0])) {
	    if (ged_view_context_obol_endpoint_get(view_ctx))
		return _fp_bool_property_set(gedp, view_ctx,
		    "view.faceplate.grid.visible", 1);
	    grid.draw = 1;
	    bv_grid_state_set(view, &grid);
	    return BRLCAD_OK;
	}
	if (BU_STR_EQUAL("0", argv[0])) {
	    if (ged_view_context_obol_endpoint_get(view_ctx))
		return _fp_bool_property_set(gedp, view_ctx,
		    "view.faceplate.grid.visible", 0);
	    grid.draw = 0;
	    bv_grid_state_set(view, &grid);
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
	if (bu_cmd_valid(_fp_grid_cmds, argv[i]) == BRLCAD_OK) {
	    cmd_pos = i;
	    break;
	}
    }

    int acnt = (cmd_pos >= 0) ? cmd_pos : argc;
    (void)bu_opt_parse(NULL, acnt, argv, d);

    struct bv_grid_state grid;
    if (!bv_grid_state_get(&grid, view))
	return BRLCAD_ERROR;
    const struct bv_grid_state initial_grid = grid;

    struct _ged_fp_grid_info ginfo;
    ginfo.gd = gd;
    ginfo.g = &grid;
    ginfo.draw_property_set = 0;

    int ret = _ged_subcmd_exec(gedp, d, _fp_grid_cmds, "view faceplate grid", "[options] subcommand [args]", (void *)&ginfo, argc, argv, help, cmd_pos);
    if (ret == BRLCAD_OK && !ginfo.draw_property_set) {
	if (ged_view_context_obol_endpoint_get(view_ctx)) {
	    if (!_fp_grid_endpoint_state_apply(view_ctx, &initial_grid, &grid))
		return BRLCAD_ERROR;
	} else {
	    bv_grid_state_set(view, &grid);
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
