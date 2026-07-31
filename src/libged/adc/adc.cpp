/*                         A D C . C P P
 * BRL-CAD
 *
 * Copyright (c) 1985-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file libged/adc.cpp
 *
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>


#include "vmath.h"
#include "bv.h"
#include "ged.h"
#include "ged/view.h"
#include "BObol/BDisplayEndpoint.h"

void ged_calc_adc_pos(const struct ged_view_context *view_ctx, struct bv_adc_state *adc);
void ged_calc_adc_a1(const struct ged_view_context *view_ctx, struct bv_adc_state *adc);
void ged_calc_adc_a2(const struct ged_view_context *view_ctx, struct bv_adc_state *adc);
void ged_calc_adc_dst(const struct ged_view_context *view_ctx, struct bv_adc_state *adc);


/* An attached view publishes ADC visibility through its endpoint.  The ADC
 * geometry remains passive per-view state, so an unattached view can still be
 * configured before a host is created. */
static int
adc_draw_set(struct ged *gedp, struct ged_view_context *view_ctx, struct bv_adc_state *adc,
	int enabled)
{
    if (!gedp || !view_ctx || !adc)
	return BRLCAD_ERROR;

    if (ged_view_context_obol_endpoint_get(view_ctx)) {
	struct bv_display_property_value value =
	    BV_DISPLAY_PROPERTY_VALUE_INIT;
	value.type = BV_DISPLAY_PROPERTY_BOOL;
	value.bool_value = enabled ? 1 : 0;
	if (ged_view_context_display_property_set(view_ctx,
		"view.faceplate.adc.visible", &value) !=
	    BV_DISPLAY_PROPERTY_OK) {
	    bu_vls_printf(gedp->ged_result_str,
		"active view has no Obol ADC visibility policy\n");
	    return BRLCAD_ERROR;
	}
    }

    adc->draw = enabled ? 1 : 0;
    return BRLCAD_OK;
}

static void
adc_vls_print(const struct ged_view_context *view_ctx, const struct bv_adc_state *adc, fastf_t base2local, struct bu_vls *out_vp)
{
    fastf_t view_scale = bv_scale_get(bv_context_view_const((const struct bv_context *)view_ctx));

    bu_vls_printf(out_vp, "draw = %d\n", adc->draw);
    bu_vls_printf(out_vp, "a1 = %.15e\n", adc->a1);
    bu_vls_printf(out_vp, "a2 = %.15e\n", adc->a2);
    bu_vls_printf(out_vp, "dst = %.15e\n", adc->dst * view_scale * base2local);
    bu_vls_printf(out_vp, "odst = %d\n", adc->dv_dist);
    bu_vls_printf(out_vp, "hv = %.15e %.15e\n",
		  adc->pos_grid[X] * view_scale * base2local,
		  adc->pos_grid[Y] * view_scale * base2local);
    bu_vls_printf(out_vp, "xyz = %.15e %.15e %.15e\n",
		  adc->pos_model[X] * base2local,
		  adc->pos_model[Y] * base2local,
		  adc->pos_model[Z] * base2local);
    bu_vls_printf(out_vp, "x = %d\n", adc->dv_x);
    bu_vls_printf(out_vp, "y = %d\n", adc->dv_y);
    bu_vls_printf(out_vp, "anchor_pos = %d\n", adc->anchor_pos);
    bu_vls_printf(out_vp, "anchor_a1 = %d\n", adc->anchor_a1);
    bu_vls_printf(out_vp, "anchor_a2 = %d\n", adc->anchor_a2);
    bu_vls_printf(out_vp, "anchor_dst = %d\n", adc->anchor_dst);
    bu_vls_printf(out_vp, "anchorpoint_a1 = %.15e %.15e %.15e\n",
		  adc->anchor_pt_a1[X] * base2local,
		  adc->anchor_pt_a1[Y] * base2local,
		  adc->anchor_pt_a1[Z] * base2local);
    bu_vls_printf(out_vp, "anchorpoint_a2 = %.15e %.15e %.15e\n",
		  adc->anchor_pt_a2[X] * base2local,
		  adc->anchor_pt_a2[Y] * base2local,
		  adc->anchor_pt_a2[Z] * base2local);
    bu_vls_printf(out_vp, "anchorpoint_dst = %.15e %.15e %.15e\n",
		  adc->anchor_pt_dst[X] * base2local,
		  adc->anchor_pt_dst[Y] * base2local,
		  adc->anchor_pt_dst[Z] * base2local);
}


static void
adc_usage(struct bu_vls *vp, const char *name)
{
    bu_vls_printf(vp, "Usage: %s \n", name);
    bu_vls_printf(vp, "%s", " adc vars                      print a list of all variables (i.e. var = val)\n");
    bu_vls_printf(vp, "%s", " adc draw [0|1]                set or get the draw parameter\n");
    bu_vls_printf(vp, "%s", " adc a1   [#]                  set or get angle1\n");
    bu_vls_printf(vp, "%s", " adc a2   [#]                  set or get angle2\n");
    bu_vls_printf(vp, "%s", " adc dst  [#]                  set or get radius (distance) of tick\n");
    bu_vls_printf(vp, "%s", " adc odst [#]                  set or get radius (distance) of tick (+-2047)\n");
    bu_vls_printf(vp, "%s", " adc hv   [# #]                set or get position (grid coordinates)\n");
    bu_vls_printf(vp, "%s", " adc xyz  [# # #]              set or get position (model coordinates)\n");
    bu_vls_printf(vp, "%s", " adc x [#]                     set or get horizontal position (+-2047)\n");
    bu_vls_printf(vp, "%s", " adc y [#]                     set or get vertical position (+-2047)\n");
    bu_vls_printf(vp, "%s", " adc dh #                      add to horizontal position (grid coordinates)\n");
    bu_vls_printf(vp, "%s", " adc dv #                      add to vertical position (grid coordinates)\n");
    bu_vls_printf(vp, "%s", " adc dx #                      add to X position (model coordinates)\n");
    bu_vls_printf(vp, "%s", " adc dy #                      add to Y position (model coordinates)\n");
    bu_vls_printf(vp, "%s", " adc dz #                      add to Z position (model coordinates)\n");
    bu_vls_printf(vp, "%s", " adc anchor_pos [0|1]          anchor ADC to current position in model coordinates\n");
    bu_vls_printf(vp, "%s", " adc anchor_a1  [0|1]          anchor angle1 to go through anchorpoint_a1\n");
    bu_vls_printf(vp, "%s", " adc anchor_a2  [0|1]          anchor angle2 to go through anchorpoint_a2\n");
    bu_vls_printf(vp, "%s", " adc anchor_dst [0|1]          anchor tick distance to go through anchorpoint_dst\n");
    bu_vls_printf(vp, "%s", " adc anchorpoint_a1  [# # #]   set or get anchor point for angle1\n");
    bu_vls_printf(vp, "%s", " adc anchorpoint_a2  [# # #]   set or get anchor point for angle2\n");
    bu_vls_printf(vp, "%s", " adc anchorpoint_dst [# # #]   set or get anchor point for tick distance\n");
    bu_vls_printf(vp, "%s", " adc -i                        any of the above appropriate commands will interpret parameters as increments\n");
    bu_vls_printf(vp, "%s", " adc reset                     reset angles, location, and tick distance\n");
    bu_vls_printf(vp, "%s", " adc help                      prints this help message\n");
}


/*
 * Note - this needs to be rewritten to accept keyword/value pairs so
 * that multiple attributes can be set with a single command call.
 */
int
ged_adc_core(struct ged *gedp,
	int argc,
	const char *argv[])
{
    char *command;
    char *parameter;
    char **argp = (char **)argv;
    double scanval;
    point_t user_pt;		/* Value(s) provided by user */
    point_t scaled_pos;
    int incr_flag;
    int i;

    GED_CHECK_VIEW(gedp, BRLCAD_ERROR);
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    fastf_t view_scale = bv_scale_get(view);
    fastf_t adc_scale = (gedp->dbip) ? view_scale * gedp->dbip->dbi_base2local : view_scale;
    double sval = (gedp->dbip) ? gedp->dbip->dbi_local2base : 1.0;
    mat_t model2view;
    mat_t view2model;
    bv_model2view_get(model2view, view);
    bv_view2model_get(view2model, view);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (argc < 2 || 6 < argc) {
	adc_usage(gedp->ged_result_str, argv[0]);
	return BRLCAD_ERROR;
    }

    command = (char *)argv[0];

    if (BU_STR_EQUAL(argv[1], "-i")) {
	if (argc < 5) {
	    bu_vls_printf(gedp->ged_result_str, "%s: -i option specified without an op-val pair", command);
	    return BRLCAD_ERROR;
	}

	incr_flag = 1;
	parameter = (char *)argv[2];
	argc -= 3;
	argp += 3;
    } else {
	incr_flag = 0;
	parameter = (char *)argv[1];
	argc -= 2;
	argp += 2;
    }

    for (i = 0; i < argc; ++i) {
	if (sscanf(argp[i], "%lf", &scanval) != 1) {
	    adc_usage(gedp->ged_result_str, command);
	    return BRLCAD_ERROR;
	}
	user_pt[i] = scanval;
    }

    struct bv_adc_state adc;
    if (!bv_adc_state_get(&adc, view))
	return BRLCAD_ERROR;
#define ADC_COMMIT_RETURN(_ret) \
    do { \
	bv_adc_state_set(view, &adc); \
	return (_ret); \
    } while (0)

    if (BU_STR_EQUAL(parameter, "draw")) {
	if (argc == 0) {
	    bu_vls_printf(gedp->ged_result_str, "%d", adc.draw);
	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 1) {
	    i = (int)user_pt[X];
	    if (adc_draw_set(gedp, view_ctx, &adc, i) != BRLCAD_OK)
		return BRLCAD_ERROR;

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s draw' command accepts 0 or 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "a1")) {
	if (argc == 0) {
	    bu_vls_printf(gedp->ged_result_str, "%g", adc.a1);
	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 1) {
	    if (!adc.anchor_a1) {
		if (incr_flag)
		    adc.a1 += user_pt[0];
		else
		    adc.a1 = user_pt[0];

		adc.dv_a1 = (1.0 - (adc.a1 / 45.0)) * BV_VIEW_MAX;
	    }

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s a1' command accepts only 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "a2")) {
	if (argc == 0) {
	    bu_vls_printf(gedp->ged_result_str, "%g", adc.a2);
	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 1) {
	    if (!adc.anchor_a2) {
		if (incr_flag)
		    adc.a2 += user_pt[0];
		else
		    adc.a2 = user_pt[0];

		adc.dv_a2 = (1.0 - (adc.a2 / 45.0)) * BV_VIEW_MAX;
	    }

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s a2' command accepts only 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "dst")) {
	if (argc == 0) {
	    bu_vls_printf(gedp->ged_result_str, "%g", adc.dst * adc_scale);
	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 1) {
	    if (!adc.anchor_dst) {
		if (incr_flag)
		    adc.dst += user_pt[0] / adc_scale;
		else
		    adc.dst = user_pt[0] / adc_scale;

		adc.dv_dist = (adc.dst / M_SQRT1_2 - 1.0) * BV_VIEW_MAX;
	    }

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s dst' command accepts 0 or 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "odst")) {
	if (argc == 0) {
	    bu_vls_printf(gedp->ged_result_str, "%d", adc.dv_dist);
	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 1) {
	    if (!adc.anchor_dst) {
		if (incr_flag)
		    adc.dv_dist += user_pt[0];
		else
		    adc.dv_dist = user_pt[0];

		adc.dst = (adc.dv_dist * BV_INV_VIEW + 1.0) * M_SQRT1_2;
	    }

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s odst' command accepts 0 or 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "dh")) {
	if (argc == 1) {
	    if (!adc.anchor_pos) {
		adc.pos_grid[X] += user_pt[0] / adc_scale;
		bv_adc_grid_to_view(&adc, view2model, BV_VIEW_MAX);
		MAT4X3PNT(adc.pos_model, view2model, adc.pos_view);
	    }

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s dh' command requires 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "dv")) {
	if (argc == 1) {
	    if (!adc.anchor_pos) {
		adc.pos_grid[Y] += user_pt[0] / adc_scale;
		bv_adc_grid_to_view(&adc, view2model, BV_VIEW_MAX);
		MAT4X3PNT(adc.pos_model, view2model, adc.pos_view);
	    }

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s dv' command requires 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "hv")) {
	if (argc == 0) {
	    bu_vls_printf(gedp->ged_result_str, "%g %g",
			  adc.pos_grid[X] * adc_scale,
			  adc.pos_grid[Y] * adc_scale);
	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 2) {
	    if (!adc.anchor_pos) {
		if (incr_flag) {
		    adc.pos_grid[X] += user_pt[X] / adc_scale;
		    adc.pos_grid[Y] += user_pt[Y] / adc_scale;
		} else {
		    adc.pos_grid[X] = user_pt[X] / adc_scale;
		    adc.pos_grid[Y] = user_pt[Y] / adc_scale;
		}

		adc.pos_grid[Z] = 0.0;
		bv_adc_grid_to_view(&adc, view2model, BV_VIEW_MAX);
		MAT4X3PNT(adc.pos_model, view2model, adc.pos_model);
	    }

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s hv' command requires 0 or 2 arguments\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "dx")) {
	if (argc == 1) {
	    if (!adc.anchor_pos) {
		double pval = (gedp->dbip) ? (user_pt[0] * gedp->dbip->dbi_local2base) : user_pt[0];
		adc.pos_model[X] += pval;
		bv_adc_model_to_view(&adc, model2view, BV_VIEW_MAX);
		bv_adc_view_to_grid(&adc, model2view);
	    }

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s dx' command requires 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "dy")) {
	if (argc == 1) {
	    if (!adc.anchor_pos) {
		double pval = (gedp->dbip) ? (user_pt[0] * gedp->dbip->dbi_local2base) : user_pt[0];
		adc.pos_model[Y] += pval;
		bv_adc_model_to_view(&adc, model2view, BV_VIEW_MAX);
		bv_adc_view_to_grid(&adc, model2view);
	    }

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s dy' command requires 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "dz")) {
	if (argc == 1) {
	    if (!adc.anchor_pos) {
		double pval = (gedp->dbip) ? (user_pt[0] * gedp->dbip->dbi_local2base) : user_pt[0];
		adc.pos_model[Z] += pval;
		bv_adc_model_to_view(&adc, model2view, BV_VIEW_MAX);
		bv_adc_view_to_grid(&adc, model2view);
	    }

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s dz' command requires 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "xyz")) {
	if (argc == 0) {
	    VSCALE(scaled_pos, adc.pos_model, sval);
	    bu_vls_printf(gedp->ged_result_str, "%g %g %g", V3ARGS(scaled_pos));
	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 3) {
	    VSCALE(user_pt, user_pt, sval);

	    if (incr_flag) {
		VADD2(adc.pos_model, adc.pos_model, user_pt);
	    } else {
		VMOVE(adc.pos_model, user_pt);
	    }

	    bv_adc_model_to_view(&adc, model2view, BV_VIEW_MAX);
	    bv_adc_view_to_grid(&adc, model2view);

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s xyz' command requires 0 or 3 arguments\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "x")) {
	if (argc == 0) {
	    bu_vls_printf(gedp->ged_result_str, "%d", adc.dv_x);
	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 1) {
	    if (!adc.anchor_pos) {
		if (incr_flag) {
		    adc.dv_x += user_pt[0];
		} else {
		    adc.dv_x = user_pt[0];
		}

		adc.pos_view[X] = adc.dv_x * BV_INV_VIEW;
		adc.pos_view[Y] = adc.dv_y * BV_INV_VIEW;
		bv_adc_view_to_grid(&adc, model2view);
		MAT4X3PNT(adc.pos_model, view2model, adc.pos_view);
	    }

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s x' command requires 0 or 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "y")) {
	if (argc == 0) {
	    bu_vls_printf(gedp->ged_result_str, "%d", adc.dv_y);
	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 1) {
	    if (!adc.anchor_pos) {
		if (incr_flag) {
		    adc.dv_y += user_pt[0];
		} else {
		    adc.dv_y = user_pt[0];
		}

		adc.pos_view[X] = adc.dv_x * BV_INV_VIEW;
		adc.pos_view[Y] = adc.dv_y * BV_INV_VIEW;
		bv_adc_view_to_grid(&adc, model2view);
		MAT4X3PNT(adc.pos_model, view2model, adc.pos_view);
	    }

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s y' command requires 0 or 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "anchor_pos")) {
	if (argc == 0) {
	    bu_vls_printf(gedp->ged_result_str, "%d", adc.anchor_pos);
	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 1) {
	    i = (int)user_pt[X];

	    if (i < 0 || 2 < i) {
		bu_vls_printf(gedp->ged_result_str, "The '%d anchor_pos' parameter accepts values of 0, 1, or 2.", i);
		return BRLCAD_ERROR;
	    }

	    adc.anchor_pos = i;
	    ged_calc_adc_pos(view_ctx, &adc);

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s anchor_pos' command accepts 0 or 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "anchor_a1")) {
	if (argc == 0) {
	    bu_vls_printf(gedp->ged_result_str, "%d", adc.anchor_a1);
	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 1) {
	    i = (int)user_pt[X];

	    if (i)
		adc.anchor_a1 = 1;
	    else
		adc.anchor_a1 = 0;

	    ged_calc_adc_a1(view_ctx, &adc);

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s anchor_a1' command accepts 0 or 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "anchorpoint_a1")) {
	if (argc == 0) {
	    VSCALE(scaled_pos, adc.anchor_pt_a1, sval);
	    bu_vls_printf(gedp->ged_result_str, "%g %g %g", V3ARGS(scaled_pos));

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 3) {
	    VSCALE(user_pt, user_pt, sval);

	    if (incr_flag) {
		VADD2(adc.anchor_pt_a1, adc.anchor_pt_a1, user_pt);
	    } else {
		VMOVE(adc.anchor_pt_a1, user_pt);
	    }

	    ged_calc_adc_a1(view_ctx, &adc);

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s anchorpoint_a1' command accepts 0 or 3 arguments\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "anchor_a2")) {
	if (argc == 0) {
	    bu_vls_printf(gedp->ged_result_str, "%d", adc.anchor_a2);

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 1) {
	    i = (int)user_pt[X];

	    if (i)
		adc.anchor_a2 = 1;
	    else
		adc.anchor_a2 = 0;

	    ged_calc_adc_a2(view_ctx, &adc);

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s anchor_a2' command accepts 0 or 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "anchorpoint_a2")) {
	if (argc == 0) {
	    VSCALE(scaled_pos, adc.anchor_pt_a2, sval);

	    bu_vls_printf(gedp->ged_result_str, "%g %g %g", V3ARGS(scaled_pos));

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 3) {
	    VSCALE(user_pt, user_pt, sval);

	    if (incr_flag) {
		VADD2(adc.anchor_pt_a2, adc.anchor_pt_a2, user_pt);
	    } else {
		VMOVE(adc.anchor_pt_a2, user_pt);
	    }

	    ged_calc_adc_a2(view_ctx, &adc);

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s anchorpoint_a2' command accepts 0 or 3 arguments\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "anchor_dst")) {
	if (argc == 0) {
	    bu_vls_printf(gedp->ged_result_str, "%d", adc.anchor_dst);

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 1) {
	    i = (int)user_pt[X];

	    if (i) {
		adc.anchor_dst = 1;
	    } else
		adc.anchor_dst = 0;

	    ged_calc_adc_dst(view_ctx, &adc);

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s anchor_dst' command accepts 0 or 1 argument\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "anchorpoint_dst")) {
	if (argc == 0) {
	    VSCALE(scaled_pos, adc.anchor_pt_dst, sval);
	    bu_vls_printf(gedp->ged_result_str, "%g %g %g", V3ARGS(scaled_pos));

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	} else if (argc == 3) {
	    VSCALE(user_pt, user_pt, sval);

	    if (incr_flag) {
		VADD2(adc.anchor_pt_dst, adc.anchor_pt_dst, user_pt);
	    } else {
		VMOVE(adc.anchor_pt_dst, user_pt);
	    }

	    ged_calc_adc_dst(view_ctx, &adc);

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s anchorpoint_dst' command accepts 0 or 3 arguments\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "reset")) {
	if (argc == 0) {
	    bv_adc_reset(&adc, view2model, model2view);

	    ADC_COMMIT_RETURN(BRLCAD_OK);
	}

	bu_vls_printf(gedp->ged_result_str, "The '%s reset' command accepts no arguments\n", command);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(parameter, "vars")) {
	adc_vls_print(view_ctx, &adc, sval, gedp->ged_result_str);
	ADC_COMMIT_RETURN(BRLCAD_OK);
    }

    if (BU_STR_EQUAL(parameter, "help")) {
	adc_usage(gedp->ged_result_str, command);
	return GED_HELP;
    }

    bu_vls_printf(gedp->ged_result_str, "%s: unrecognized command '%s'\n", command, parameter);
    adc_usage(gedp->ged_result_str, command);

    return BRLCAD_ERROR;
}


void
ged_calc_adc_pos(const struct ged_view_context *view_ctx, struct bv_adc_state *adc)
{
    mat_t model2view;
    mat_t view2model;
    const struct bv *view = bv_context_view_const((const struct bv_context *)view_ctx);
    bv_model2view_get(model2view, view);
    bv_view2model_get(view2model, view);

    if (adc->anchor_pos == 1) {
	bv_adc_model_to_view(adc, model2view, BV_VIEW_MAX);
	bv_adc_view_to_grid(adc, model2view);
    } else if (adc->anchor_pos == 2) {
	bv_adc_grid_to_view(adc, view2model, BV_VIEW_MAX);
	MAT4X3PNT(adc->pos_model, view2model, adc->pos_view);
    } else {
	bv_adc_view_to_grid(adc, model2view);
	MAT4X3PNT(adc->pos_model, view2model, adc->pos_view);
    }
}


void
ged_calc_adc_a1(const struct ged_view_context *view_ctx, struct bv_adc_state *adc)
{
    if (adc->anchor_a1) {
	fastf_t dx, dy;
	mat_t model2view;
	point_t view_pt;

	bv_model2view_get(model2view,
		bv_context_view_const((const struct bv_context *)view_ctx));
	MAT4X3PNT(view_pt, model2view, adc->anchor_pt_a1);
	dx = view_pt[X] * BV_VIEW_MAX - adc->dv_x;
	dy = view_pt[Y] * BV_VIEW_MAX - adc->dv_y;

	if (!ZERO(dx) || !ZERO(dy)) {
	    adc->a1 = RAD2DEG*atan2(dy, dx);
	    adc->dv_a1 = (1.0 - (adc->a1 / 45.0)) * BV_VIEW_MAX;
	}
    }
}


void
ged_calc_adc_a2(const struct ged_view_context *view_ctx, struct bv_adc_state *adc)
{
    if (adc->anchor_a2) {
	fastf_t dx, dy;
	mat_t model2view;
	point_t view_pt;

	bv_model2view_get(model2view,
		bv_context_view_const((const struct bv_context *)view_ctx));
	MAT4X3PNT(view_pt, model2view, adc->anchor_pt_a2);
	dx = view_pt[X] * BV_VIEW_MAX - adc->dv_x;
	dy = view_pt[Y] * BV_VIEW_MAX - adc->dv_y;

	if (!ZERO(dx) || !ZERO(dy)) {
	    adc->a2 = RAD2DEG*atan2(dy, dx);
	    adc->dv_a2 = (1.0 - (adc->a2 / 45.0)) * BV_VIEW_MAX;
	}
    }
}


void
ged_calc_adc_dst(const struct ged_view_context *view_ctx, struct bv_adc_state *adc)
{
    if (adc->anchor_dst) {
	fastf_t dist;
	fastf_t dx, dy;
	mat_t model2view;
	point_t view_pt;

	bv_model2view_get(model2view,
		bv_context_view_const((const struct bv_context *)view_ctx));
	MAT4X3PNT(view_pt, model2view, adc->anchor_pt_dst);

	dx = view_pt[X] * BV_VIEW_MAX - adc->dv_x;
	dy = view_pt[Y] * BV_VIEW_MAX - adc->dv_y;
	dist = sqrt(dx * dx + dy * dy);
	adc->dst = dist * BV_INV_VIEW;
	adc->dv_dist = (dist / M_SQRT1_2) - BV_VIEW_MAX;
    } else
	adc->dst = (adc->dv_dist * BV_INV_VIEW + 1.0) * M_SQRT1_2;
}


#include "../include/plugin.h"

#define GED_ADC_COMMANDS(X, XID) \
    X(adc, ged_adc_core, GED_CMD_DEFAULT) \

GED_DECLARE_COMMAND_SET(GED_ADC_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_adc", 1, GED_ADC_COMMANDS)

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
