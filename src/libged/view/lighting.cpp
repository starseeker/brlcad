/*                      L I G H T I N G . C P P
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
/** @file libged/view/lighting.cpp
 *
 * The "view lighting" subcommand: enable/disable the camera-driven headlight
 * and in-scene (database) lights, and control whether the headlight tracks the
 * camera.  Settings are stored per-view in struct bv_lighting_state and pushed
 * to the live Obol renderer via ged_draw_obol_lighting_sync().
 *
 */

#include "common.h"

#include <stdlib.h>
#include <string.h>

#include "bu/cmd.h"
#include "bu/opt.h"
#include "bu/vls.h"
#include "bv.h"
#include "ged/draw_obol.h"

#include "../ged_private.h"
#include "./ged_view.h"

#define HELPFLAG "--print-help"
#define PURPOSEFLAG "--print-purpose"

struct _ged_lgt_info {
    struct ged *gedp;
    long verbosity;
    const struct bu_cmdtab *cmds;
    struct bu_opt_desc *gopts;
};

static int
_lgt_cmd_msgs(void *bs, int argc, const char **argv, const char *us, const char *ps)
{
    struct _ged_lgt_info *gd = (struct _ged_lgt_info *)bs;
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
_lgt_bool_argument(struct ged *gedp, const char *argument, int *enabled)
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

/* Shared handler for the boolean lighting toggles.  field_offset selects which
 * int in struct bv_lighting_state is read/written. */
static int
_lgt_bool_field(void *ds, int argc, const char **argv,
	const char *usage_string, const char *purpose_string,
	size_t field_offset)
{
    if (_lgt_cmd_msgs(ds, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    struct _ged_lgt_info *gd = (struct _ged_lgt_info *)ds;
    struct ged *gedp = gd->gedp;
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    struct bv_lighting_state lighting;
    if (!view || !bv_lighting_state_get(&lighting, view))
	return BRLCAD_ERROR;
    int *field = (int *)((char *)&lighting + field_offset);

    if (!argc) {
	bu_vls_printf(gedp->ged_result_str, "%d", *field);
	return BRLCAD_OK;
    }

    if (argc == 1) {
	int enabled = 0;
	if (!_lgt_bool_argument(gedp, argv[0], &enabled))
	    return BRLCAD_ERROR;
	*field = enabled;
	if (!bv_lighting_state_set(view, &lighting))
	    return BRLCAD_ERROR;
	return BRLCAD_OK;
    }

    bu_vls_printf(gedp->ged_result_str, "invalid command\n");
    return BRLCAD_ERROR;
}

static int
_lgt_cmd_headlight(void *ds, int argc, const char **argv)
{
    return _lgt_bool_field(ds, argc, argv,
	"lighting headlight [0|1]",
	"Enable/disable the camera-driven headlight.",
	offsetof(struct bv_lighting_state, headlight_enabled));
}

static int
_lgt_cmd_tracking(void *ds, int argc, const char **argv)
{
    return _lgt_bool_field(ds, argc, argv,
	"lighting tracking [0|1]",
	"Headlight follows the camera (1) vs stays scene-fixed (0).",
	offsetof(struct bv_lighting_state, headlight_tracks_camera));
}

static int
_lgt_cmd_scene(void *ds, int argc, const char **argv)
{
    return _lgt_bool_field(ds, argc, argv,
	"lighting scene [0|1]",
	"Enable/disable in-scene (database) light sources.",
	offsetof(struct bv_lighting_state, scene_lights_enabled));
}

static int
_lgt_cmd_offset(void *ds, int argc, const char **argv)
{
    const char *usage_string = "lighting offset [x y z]";
    const char *purpose_string =
	"Get/set the eye-space headlight offset direction.";
    if (_lgt_cmd_msgs(ds, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    struct _ged_lgt_info *gd = (struct _ged_lgt_info *)ds;
    struct ged *gedp = gd->gedp;
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    struct bv_lighting_state lighting;
    if (!view || !bv_lighting_state_get(&lighting, view))
	return BRLCAD_ERROR;

    if (!argc) {
	bu_vls_printf(gedp->ged_result_str, "%g %g %g",
	    lighting.headlight_offset[0], lighting.headlight_offset[1],
	    lighting.headlight_offset[2]);
	return BRLCAD_OK;
    }
    if (argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    lighting.headlight_offset[0] = atof(argv[0]);
    lighting.headlight_offset[1] = atof(argv[1]);
    lighting.headlight_offset[2] = atof(argv[2]);
    if (!bv_lighting_state_set(view, &lighting))
	return BRLCAD_ERROR;
    return BRLCAD_OK;
}

static const struct bu_cmdtab _lgt_cmds[] = {
    { "headlight",  _lgt_cmd_headlight},
    { "tracking",   _lgt_cmd_tracking},
    { "scene",      _lgt_cmd_scene},
    { "offset",     _lgt_cmd_offset},
    { (char *)NULL, NULL}
};

int
ged_lighting_core(struct ged *gedp, int argc, const char *argv[])
{
    int help = 0;
    struct _ged_lgt_info gd;
    gd.gedp = gedp;
    gd.cmds = _lgt_cmds;
    gd.verbosity = 0;

    if (UNLIKELY(!gedp || !argc || !argv))
	return BRLCAD_ERROR;

    bu_vls_trunc(gedp->ged_result_str, 0);

    /* consume the "lighting" command word */
    argc--; argv++;

    struct bu_opt_desc d[3];
    BU_OPT(d[0], "h", "help",    "",  NULL,               &help,         "Print help");
    BU_OPT(d[1], "v", "verbose", "",  &bu_opt_incr_long,  &gd.verbosity, "Verbose output");
    BU_OPT_NULL(d[2]);
    gd.gopts = d;

    int ac = bu_opt_parse(NULL, argc, argv, d);

    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, ": no current view set");
	return BRLCAD_ERROR;
    }

    /* No subcommand: report the full lighting state (or help with -h). */
    if (!ac) {
	if (help) {
	    _ged_subcmd_help(gedp, (struct bu_opt_desc *)d,
		(const struct bu_cmdtab *)_lgt_cmds,
		"view lighting", "[options] subcommand [args]", &gd, 0, NULL);
	    return BRLCAD_OK;
	}
	struct bv *view = bv_context_view((struct bv_context *)view_ctx);
	struct bv_lighting_state lighting;
	if (!view || !bv_lighting_state_get(&lighting, view))
	    return BRLCAD_ERROR;
	bu_vls_printf(gedp->ged_result_str,
	    "headlight %d\ntracking %d\nscene %d\noffset %g %g %g",
	    lighting.headlight_enabled, lighting.headlight_tracks_camera,
	    lighting.scene_lights_enabled,
	    lighting.headlight_offset[0], lighting.headlight_offset[1],
	    lighting.headlight_offset[2]);
	return BRLCAD_OK;
    }

    int ret;
    if (bu_cmd(_lgt_cmds, ac, argv, 0, (void *)&gd, &ret) == BRLCAD_OK) {
	if (ret == BRLCAD_OK) {
	    (void)ged_draw_obol_lighting_sync(gedp, view_ctx);
	    /* A lighting change alters the rendered appearance without changing
	     * geometry or the view matrix, so nothing else marks the view dirty.
	     * Request a view refresh and invoke the application's refresh handler
	     * so mged/qged repaint immediately instead of on the next nudge. */
	    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
	    if (view)
		(void)bv_refresh_request(view, GED_VIEW_REFRESH_VIEW);
	    ged_refresh_cb(gedp);
	}
	return ret;
    }

    bu_vls_printf(gedp->ged_result_str, "subcommand %s not defined", argv[0]);
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
