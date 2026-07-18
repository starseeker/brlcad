/*                         K N O B . C
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
/** @file libged/view/knob.c
 *
 * The view knob subcommand.
 *
 */

#include "common.h"

#include <string.h>
#include "bv.h"
#include "../ged_private.h"
#include "./ged_view.h"

static void
print_knob_vals(struct bu_vls *o, struct ged_view_context *view_ctx)
{
    struct bv_knob_values values;

    if (!o || !bv_knob_values_get(&values,
	    bv_context_view_const((const struct bv_context *)view_ctx)))
	return;

    bu_vls_printf(o, "rate - rotation:\n");
    bu_vls_printf(o, " x = %f\n", values.rate_rotation[X]);
    bu_vls_printf(o, " y = %f\n", values.rate_rotation[Y]);
    bu_vls_printf(o, " z = %f\n", values.rate_rotation[Z]);
    bu_vls_printf(o, "rate - translation (skew):\n");
    bu_vls_printf(o, " X = %f\n", values.rate_translation[X]);
    bu_vls_printf(o, " Y = %f\n", values.rate_translation[Y]);
    bu_vls_printf(o, " Z = %f\n", values.rate_translation[Z]);
    bu_vls_printf(o, "rate - scale:\n");
    bu_vls_printf(o, " S = %f\n", values.rate_scale);

    bu_vls_printf(o, "absolute - rotation:\n");
    bu_vls_printf(o, " ax = %f\n", values.absolute_rotation[X]);
    bu_vls_printf(o, " ay = %f\n", values.absolute_rotation[Y]);
    bu_vls_printf(o, " az = %f\n", values.absolute_rotation[Z]);
    bu_vls_printf(o, "absolute - translation (skew):\n");
    bu_vls_printf(o, " aX = %f\n", values.absolute_translation[X]);
    bu_vls_printf(o, " aY = %f\n", values.absolute_translation[Y]);
    bu_vls_printf(o, " aZ = %f\n", values.absolute_translation[Z]);
    bu_vls_printf(o, "absolute - scale:\n");
    bu_vls_printf(o, " aS = %f\n", values.absolute_scale);
}

int
ged_knob_core(struct ged *gedp, int argc, const char *argv[])
{
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);
    GED_CHECK_VIEW(gedp, BRLCAD_ERROR);
    /* Make sure the view coordinate conversion values match the database */
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    bv_unit_conversion_set(view,
	(gedp->dbip) ? gedp->dbip->dbi_local2base : 1.0,
	(gedp->dbip) ? gedp->dbip->dbi_base2local : 1.0);

    static const char *usage =
	"knob [options] x # y # z # for rotation in degrees (may specify one or more of x, y, z)\n"
	"knob [options] S #for scale\n"
       	"knob [options] X # Y # Z # for translation (a.k.a slew) (rates, range -1..+1, may specify one or more of X, Y, Z))\n"
	"knob [options] ax # ay # az # for absolute rotation in degrees (may specify one or more of ax, ay, az)\n"
	"knob [options] aS # for absolute scale.\n"
	"knob [options] aX, aY, aZ for absolute translate (a.k.a slew) (may specify one or more of ax, ay, az).\n"
	"knob calibrate to set current absolute translate (a.k.a slew) to 0\n"
	"knob zero|clear|reset to cancel motion\n";

    int print_help = 0;
    int incr_flag = 0;  /* interpret values as increments */
    int view_flag = 0;  /* manipulate view using view coords */
    int model_flag = 0; /* manipulate view using model coords */
    char origin = 'v';

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct bu_opt_desc d[7];
    BU_OPT(d[0], "h", "help",         "",        NULL,          &print_help,   "Print help and exit");
    BU_OPT(d[1], "?", "",             "",        NULL,          &print_help,   "");
    BU_OPT(d[2], "i", "incr",         "",        NULL,          &incr_flag,    "Interpret values as increments");
    BU_OPT(d[3], "v", "view-coords",  "",        NULL,          &view_flag,    "Manipulate view using view coordinates.  Conflicts with 'm' If neither v or m options are specified, current gv_coord setting from parent view is used.");
    BU_OPT(d[4], "m", "model-coords", "",        NULL,          &model_flag,   "Manipulate view using model coordinates.  Conflicts with 'v'");
    BU_OPT(d[5], "o", "origin",       "<char>",  &bu_opt_char,  &origin,       "Specify origin mode - may be 'e' (eye_pt), 'm' (model origin), 'v' (view origin - default), or 'k' (keypoint - set by keypoint cmd).");
    BU_OPT_NULL(d[6]);

    /* parse standard options */
    int opt_ret = bu_opt_parse(NULL, argc, argv, d);

    if (print_help) {
        _ged_cmd_help(gedp, usage, d);
        return BRLCAD_OK;
    }

    /* adjust argc to match the leftovers of the options parsing */
    argc = opt_ret;
    if (argc == 1) {
	// print current values
	print_knob_vals(gedp->ged_result_str, view_ctx);
	return BRLCAD_OK;
    }

    // Make sure model flag is set, if that's what either the user
    // options or view settings indicate we should use
    if (model_flag && view_flag) {
        _ged_cmd_help(gedp, usage, d);
	return BRLCAD_ERROR;
    }
    if (!model_flag && !view_flag)
	model_flag = (bv_coord_get(view) == 'm') ? 1 : 0;

    int do_tran = 0;
    int do_rot = 0;
    vect_t tvec;
    vect_t rvec;
    VSETALL(tvec, 0.0);
    VSETALL(rvec, 0.0);

    for (--argc, ++argv; argc; --argc, ++argv) {

	const char *cmd = *argv;

	if (BU_STR_EQUAL(cmd, "zap") || BU_STR_EQUAL(cmd, "zero") ||
		BU_STR_EQUAL(cmd, "clear") || BU_STR_EQUAL(cmd, "stop")) {
	    // Per MGED, this command seems to reset the rate entries
	    // but not the absolute entries.
	    bv_knobs_reset(&view->knobs, BV_KNOBS_RATE);
	    continue;
	}
	if (BU_STR_EQUAL(cmd, "calibrate")) {
	    bv_knobs_calibrate(view);
	    continue;
	}

	if (argc < 2) {
	    _ged_cmd_help(gedp, usage, d);
	    return BRLCAD_ERROR;
	}

	// Get our number
	struct bu_vls fmsg = BU_VLS_INIT_ZERO;
	fastf_t f;
	if (bu_opt_fastf_t(&fmsg, 1, (const char **)&argv[1], (void *)&f) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "Unable to parse numerical argument: '%s':%s\n", argv[1], bu_vls_cstr(&fmsg));
	    bu_vls_free(&fmsg);
	    return BRLCAD_ERROR;
	}
	bu_vls_free(&fmsg);

	// Read in the numerical argument, advance the array
	--argc;	++argv;

	// Process the actual command
	int kp_ret = bv_knobs_cmd_process(
		&rvec, &do_rot, &tvec, &do_tran,
		view, cmd, f,
		origin, model_flag, incr_flag);

	if (kp_ret != BRLCAD_OK) {
	    _ged_cmd_help(gedp, usage, d);
	    return BRLCAD_ERROR;
	}
    }

    if (do_tran) {
	bv_knobs_translate(view, tvec, model_flag);
    }

    if (do_rot) {
	point_t keypoint;
	const pointp_t pvt_pt = (origin == 'k') ? keypoint : NULL;

	if (pvt_pt)
	    bv_keypoint_get(keypoint, view);

	// Note - we don't (currently) support 'o' coords here, so the obj_rot matrix is always NULL.
	bv_knobs_rotate(view, rvec, origin,
	    (model_flag ? 'm' : 'v'), NULL, pvt_pt);
    }

    bv_knobs_update_rate_flags(view);

    return BRLCAD_OK;
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
