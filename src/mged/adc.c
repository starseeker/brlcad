/*                           A D C . C
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
/** @file mged/adc.c
 *
 */

#include "common.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "bu/vls.h"
#include "vmath.h"
#include "ged.h"
#include "ged/view.h"
#include "rt/view.h"
#include "./mged.h"
#include "./mged_display.h"


static char adc_syntax1[] = "\
 adc			toggle display of angle/distance cursor\n\
 adc vars		print a list of all variables (i.e. var = val)\n\
 adc draw [0|1]		set or get the draw parameter\n\
 adc a1 [#]		set or get angle1\n\
 adc a2 [#]		set or get angle2\n\
 adc dst [#]		set or get radius (distance) of tick\n\
 adc odst [#]		set or get radius (distance) of tick (+-2047)\n\
 adc hv [# #]		set or get position (grid coordinates)\n\
 adc xyz [# # #]	set or get position (model coordinates)\n\
";

static char adc_syntax2[] = "\
 adc x [#]		set or get horizontal position (+-2047)\n\
 adc y [#]		set or get vertical position (+-2047)\n\
 adc dh #		add to horizontal position (grid coordinates)\n\
 adc dv #		add to vertical position (grid coordinates)\n\
 adc dx #		add to X position (model coordinates)\n\
 adc dy #		add to Y position (model coordinates)\n\
 adc dz #		add to Z position (model coordinates)\n\
";

static char adc_syntax3[] = "\
 adc anchor_pos	[0|1]	anchor ADC to current position in model coordinates\n\
 adc anchor_a1	[0|1]	anchor angle1 to go through anchorpoint_a1\n\
 adc anchor_a2	[0|1]	anchor angle2 to go through anchorpoint_a2\n\
 adc anchor_dst	[0|1]	anchor tick distance to go through anchorpoint_dst\n\
 adc anchorpoint_a1 [# # #]	set or get anchor point for angle1\n\
 adc anchorpoint_a2 [# # #]	set or get anchor point for angle2\n\
 adc anchorpoint_dst [# # #]	set or get anchor point for tick distance\n\
";

static char adc_syntax4[] = "\
 adc -i			any of the above appropriate commands will interpret parameters as increments\n\
 adc reset		reset angles, location, and tick distance\n\
 adc help		prints this help message\n\
";


void
adc_set_dirty_flag(struct mged_state *s)
{

    for (size_t i = 0; i < BU_PTBL_LEN(&active_display_set); i++) {
	struct mged_display *m_dmp = (struct mged_display *)BU_PTBL_GET(&active_display_set, i);
	if (mged_display_view_settings_shared(m_dmp, s->mged_curr_display)) {
	    mged_display_repaint_request(m_dmp, MGED_REPAINT_DEVICE_SETTING);
	}
    }
}


void
adc_set_scroll(struct mged_state *s)
{
    struct mged_display *save_m_dmp = s->mged_curr_display;

    for (size_t i = 0; i < BU_PTBL_LEN(&active_display_set); i++) {
	struct mged_display *m_dmp = (struct mged_display *)BU_PTBL_GET(&active_display_set, i);
	if (mged_display_view_settings_shared(m_dmp, save_m_dmp)) {
	    mged_current_display_set(s, m_dmp);
	    set_scroll(s);
	    mged_display_repaint_request(s->mged_curr_display, MGED_REPAINT_DEVICE_SETTING);
	}
    }

    mged_current_display_set(s, save_m_dmp);
}


static fastf_t
adc_view_local_scale(struct mged_state *s)
{
    void *view_ctx = view_state->vs_gvp;
    return bv_scale_get(mged_view_context_view_const(view_ctx)) *
	s->dbip->dbi_base2local;
}


static void
adc_model_To_adc_view(struct mged_state *s, struct bv_adc_state *adc)
{
    mat_t model2view;
    void *view_ctx = view_state->vs_gvp;

    bv_model2view_get(model2view, mged_view_context_view_const(view_ctx));
    bv_adc_model_to_view(adc, model2view, RT_VIEW_MAX);
}


static void
adc_grid_To_adc_view(struct mged_state *s, struct bv_adc_state *adc)
{
    mat_t model2view;
    void *view_ctx = view_state->vs_gvp;

    bv_model2view_get(model2view, mged_view_context_view_const(view_ctx));
    bv_adc_grid_to_view(adc, model2view, RT_VIEW_MAX);
}


static void
adc_view_To_adc_grid(struct mged_state *s, struct bv_adc_state *adc)
{
    mat_t model2view;
    void *view_ctx = view_state->vs_gvp;

    bv_model2view_get(model2view, mged_view_context_view_const(view_ctx));
    bv_adc_view_to_grid(adc, model2view);
}


static void
calc_adc_pos(struct mged_state *s, struct bv_adc_state *adc)
{
    mat_t view2model;
    void *view_ctx = view_state->vs_gvp;

    if (adc->anchor_pos == 1) {
	adc_model_To_adc_view(s, adc);
	adc_view_To_adc_grid(s, adc);
    } else if (adc->anchor_pos == 2) {
	adc_grid_To_adc_view(s, adc);
	bv_view2model_get(view2model, mged_view_context_view_const(view_ctx));
	MAT4X3PNT(adc->pos_model, view2model, adc->pos_view);
    } else {
	adc_view_To_adc_grid(s, adc);
	bv_view2model_get(view2model, mged_view_context_view_const(view_ctx));
	MAT4X3PNT(adc->pos_model, view2model, adc->pos_view);
    }
}


static void
calc_adc_a1(struct mged_state *s, struct bv_adc_state *adc)
{
    if (adc->anchor_a1) {
	fastf_t dx, dy;
	point_t view_pt;
	mat_t model2view;
	void *view_ctx = view_state->vs_gvp;

	bv_model2view_get(model2view, mged_view_context_view_const(view_ctx));
	MAT4X3PNT(view_pt, model2view, adc->anchor_pt_a1);
	dx = view_pt[X] * RT_VIEW_MAX - adc->dv_x;
	dy = view_pt[Y] * RT_VIEW_MAX - adc->dv_y;

	if (!ZERO(dx) || !ZERO(dy)) {
	    adc->a1 = RAD2DEG*atan2(dy, dx);
	    adc->dv_a1 = (1.0 - (adc->a1 / 45.0)) * RT_VIEW_MAX;
	}
    }
}


static void
calc_adc_a2(struct mged_state *s, struct bv_adc_state *adc)
{
    if (adc->anchor_a2) {
	fastf_t dx, dy;
	point_t view_pt;
	mat_t model2view;
	void *view_ctx = view_state->vs_gvp;

	bv_model2view_get(model2view, mged_view_context_view_const(view_ctx));
	MAT4X3PNT(view_pt, model2view, adc->anchor_pt_a2);
	dx = view_pt[X] * RT_VIEW_MAX - adc->dv_x;
	dy = view_pt[Y] * RT_VIEW_MAX - adc->dv_y;

	if (!ZERO(dx) || !ZERO(dy)) {
	    adc->a2 = RAD2DEG*atan2(dy, dx);
	    adc->dv_a2 = (1.0 - (adc->a2 / 45.0)) * RT_VIEW_MAX;
	}
    }
}


static void
calc_adc_dst(struct mged_state *s, struct bv_adc_state *adc)
{
    if (adc->anchor_dst) {
	fastf_t dist;
	fastf_t dx, dy;
	point_t view_pt;
	mat_t model2view;
	void *view_ctx = view_state->vs_gvp;

	bv_model2view_get(model2view, mged_view_context_view_const(view_ctx));
	MAT4X3PNT(view_pt, model2view, adc->anchor_pt_dst);

	dx = view_pt[X] * RT_VIEW_MAX - adc->dv_x;
	dy = view_pt[Y] * RT_VIEW_MAX - adc->dv_y;
	dist = sqrt(dx * dx + dy * dy);
	adc->dst = dist * RT_INV_VIEW;
	adc->dv_dist = (dist / M_SQRT1_2) - RT_VIEW_MAX;
    } else
	adc->dst = (adc->dv_dist * RT_INV_VIEW + 1.0) * M_SQRT1_2;
}


void
mged_adc_state_refresh(struct mged_state *s)
{
    struct bv_adc_state adc_record;
    struct bv_adc_state *adc = &adc_record;

    if (!mged_display_adc_state_get(s->mged_curr_display, adc))
	return;

    calc_adc_pos(s, adc);
    calc_adc_a1(s, adc);
    calc_adc_a2(s, adc);
    calc_adc_dst(s, adc);
    mged_display_adc_state_set(s->mged_curr_display, adc);
}


static void
mged_adc_reset(struct mged_state *s, struct bv_adc_state *adc)
{
    mat_t model2view;
    mat_t view2model;
    void *view_ctx = view_state->vs_gvp;

    bv_model2view_get(model2view, mged_view_context_view_const(view_ctx));
    bv_view2model_get(view2model, mged_view_context_view_const(view_ctx));
    bv_adc_reset(adc, view2model, model2view);
}


static void
adc_print_vars(struct mged_state *s, struct bv_adc_state *adc)
{
    struct bu_vls vls = BU_VLS_INIT_ZERO;
    fastf_t view_local_scale = adc_view_local_scale(s);

    bu_vls_printf(&vls, "draw = %d\n", adc->draw);
    bu_vls_printf(&vls, "a1 = %.15e\n", adc->a1);
    bu_vls_printf(&vls, "a2 = %.15e\n", adc->a2);
    bu_vls_printf(&vls, "dst = %.15e\n", adc->dst * view_local_scale);
    bu_vls_printf(&vls, "odst = %d\n", adc->dv_dist);
    bu_vls_printf(&vls, "hv = %.15e %.15e\n",
		  adc->pos_grid[X] * view_local_scale,
		  adc->pos_grid[Y] * view_local_scale);
    bu_vls_printf(&vls, "xyz = %.15e %.15e %.15e\n",
		  adc->pos_model[X] * s->dbip->dbi_base2local,
		  adc->pos_model[Y] * s->dbip->dbi_base2local,
		  adc->pos_model[Z] * s->dbip->dbi_base2local);
    bu_vls_printf(&vls, "x = %d\n", adc->dv_x);
    bu_vls_printf(&vls, "y = %d\n", adc->dv_y);
    bu_vls_printf(&vls, "anchor_pos = %d\n", adc->anchor_pos);
    bu_vls_printf(&vls, "anchor_a1 = %d\n", adc->anchor_a1);
    bu_vls_printf(&vls, "anchor_a2 = %d\n", adc->anchor_a2);
    bu_vls_printf(&vls, "anchor_dst = %d\n", adc->anchor_dst);
    bu_vls_printf(&vls, "anchorpoint_a1 = %.15e %.15e %.15e\n",
		  adc->anchor_pt_a1[X] * s->dbip->dbi_base2local,
		  adc->anchor_pt_a1[Y] * s->dbip->dbi_base2local,
		  adc->anchor_pt_a1[Z] * s->dbip->dbi_base2local);
    bu_vls_printf(&vls, "anchorpoint_a2 = %.15e %.15e %.15e\n",
		  adc->anchor_pt_a2[X] * s->dbip->dbi_base2local,
		  adc->anchor_pt_a2[Y] * s->dbip->dbi_base2local,
		  adc->anchor_pt_a2[Z] * s->dbip->dbi_base2local);
    bu_vls_printf(&vls, "anchorpoint_dst = %.15e %.15e %.15e\n",
		  adc->anchor_pt_dst[X] * s->dbip->dbi_base2local,
		  adc->anchor_pt_dst[Y] * s->dbip->dbi_base2local,
		  adc->anchor_pt_dst[Z] * s->dbip->dbi_base2local);
    Tcl_AppendResult(s->interp, bu_vls_addr(&vls), (char *)NULL);
    bu_vls_free(&vls);
}


int
f_adc (
    ClientData clientData,
    Tcl_Interp *interp,
    int argc,
    const char *argv[])
{
    struct cmdtab *ctp = (struct cmdtab *)clientData;
    MGED_CK_CMD(ctp);
    struct mged_state *s = ctp->s;

    struct bu_vls vls = BU_VLS_INIT_ZERO;
    const char *parameter;
    const char **argp = argv;
    struct bv_adc_state adc_record;
    struct bv_adc_state *adc = &adc_record;
    point_t user_pt;		/* Value(s) provided by user */
    point_t scaled_pos;
    fastf_t view_local_scale;
    mat_t view2model;
    int incr_flag;
    int i;
    void *view_ctx;

    CHECK_DBI_NULL;

    view_ctx = view_state->vs_gvp;

#define ADC_RETURN(_ret) do { \
	mged_display_adc_state_set(s->mged_curr_display, adc); \
	return (_ret); \
    } while (0)

    if (!mged_display_adc_state_get(s->mged_curr_display, adc))
	return TCL_ERROR;

    view_local_scale = adc_view_local_scale(s);
    bv_view2model_get(view2model, mged_view_context_view_const(view_ctx));

    if (6 < argc) {
	bu_vls_printf(&vls, "help adc");
	Tcl_Eval(interp, bu_vls_addr(&vls));
	bu_vls_free(&vls);

	ADC_RETURN(TCL_ERROR);
    }

    if (argc == 1) {
	if (adc->draw)
	    adc->draw = 0;
	else
	    adc->draw = 1;

	if (adc_auto) {
	    mged_adc_reset(s, adc);
	    adc_auto = 0;
	}

	adc_set_scroll(s);

	ADC_RETURN(TCL_OK);
    }

    if (BU_STR_EQUAL(argv[1], "-i")) {
	if (argc < 4) {
	    bu_vls_printf(&vls, "adc: -i option specified without an op-val pair");
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_ERROR);
	}

	incr_flag = 1;
	parameter = argv[2];
	argc -= 3;
	argp += 3;
    } else {
	incr_flag = 0;
	parameter = argv[1];
	argc -= 2;
	argp += 2;
    }

    for (i = 0; i < argc; ++i)
	user_pt[i] = atof(argp[i]);

    if (BU_STR_EQUAL(parameter, "draw")) {
	if (argc == 0) {
	    bu_vls_printf(&vls, "%d", adc->draw);
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 1) {
	    i = (int)user_pt[X];

	    if (!mged_display_adc_visibility_set(s->mged_curr_display, i))
		ADC_RETURN(TCL_ERROR);
	    adc->draw = i ? 1 : 0;

	    adc_set_scroll(s);

	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc draw' command accepts 0 or 1 argument\n", (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "a1")) {
	if (argc == 0) {
	    bu_vls_printf(&vls, "%.15e", adc->a1);
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 1) {
	    if (!adc->anchor_a1) {
		if (incr_flag)
		    adc->a1 += user_pt[0];
		else
		    adc->a1 = user_pt[0];

		adc->dv_a1 = (1.0 - (adc->a1 / 45.0)) * RT_VIEW_MAX;
		adc_set_dirty_flag(s);
	    }

	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc a1' command accepts only 1 argument\n", (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "a2")) {
	if (argc == 0) {
	    bu_vls_printf(&vls, "%.15e", adc->a2);
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 1) {
	    if (!adc->anchor_a2) {
		if (incr_flag)
		    adc->a2 += user_pt[0];
		else
		    adc->a2 = user_pt[0];

		adc->dv_a2 = (1.0 - (adc->a2 / 45.0)) * RT_VIEW_MAX;
		adc_set_dirty_flag(s);
	    }

	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc a2' command accepts 0 or 1 argument\n", (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "dst")) {
	if (argc == 0) {
	    bu_vls_printf(&vls, "%.15e", adc->dst * view_local_scale);
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 1) {
	    if (!adc->anchor_dst) {
		if (incr_flag)
		    adc->dst += user_pt[0] / view_local_scale;
		else
		    adc->dst = user_pt[0] / view_local_scale;

		adc->dv_dist = (adc->dst / M_SQRT1_2 - 1.0) * RT_VIEW_MAX;

		adc_set_dirty_flag(s);
	    }

	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc dst' command accepts 0 or 1 argument\n", (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "odst")) {
	if (argc == 0) {
	    bu_vls_printf(&vls, "%d", adc->dv_dist);
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 1) {
	    if (!adc->anchor_dst) {
		if (incr_flag)
		    adc->dv_dist += user_pt[0];
		else
		    adc->dv_dist = user_pt[0];

		adc->dst = (adc->dv_dist * RT_INV_VIEW + 1.0) * M_SQRT1_2;
		adc_set_dirty_flag(s);
	    }

	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc odst' command accepts 0 or 1 argument\n", (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "dh")) {
	if (argc == 1) {
	    if (!adc->anchor_pos) {
		adc->pos_grid[X] += user_pt[0] / view_local_scale;
		adc_grid_To_adc_view(s, adc);
		MAT4X3PNT(adc->pos_model, view2model, adc->pos_view);

		adc_set_dirty_flag(s);
	    }

	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc dh' command requires 1 argument\n", (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "dv")) {
	if (argc == 1) {
	    if (!adc->anchor_pos) {
		adc->pos_grid[Y] += user_pt[0] / view_local_scale;
		adc_grid_To_adc_view(s, adc);
		MAT4X3PNT(adc->pos_model, view2model, adc->pos_view);

		adc_set_dirty_flag(s);
	    }

	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc dv' command requires 1 argument\n", (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "hv")) {
	if (argc == 0) {
	    bu_vls_printf(&vls, "%.15e %.15e",
			  adc->pos_grid[X] * view_local_scale,
			  adc->pos_grid[Y] * view_local_scale);
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 2) {
	    if (!adc->anchor_pos) {
		if (incr_flag) {
		    adc->pos_grid[X] += user_pt[X] / view_local_scale;
		    adc->pos_grid[Y] += user_pt[Y] / view_local_scale;
		} else {
		    adc->pos_grid[X] = user_pt[X] / view_local_scale;
		    adc->pos_grid[Y] = user_pt[Y] / view_local_scale;
		}

		adc->pos_grid[Z] = 0.0;
		adc_grid_To_adc_view(s, adc);
		MAT4X3PNT(adc->pos_model, view2model, adc->pos_model);

		adc_set_dirty_flag(s);
	    }

	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc hv' command requires 0 or 2 arguments\n", (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "dx")) {
	if (argc == 1) {
	    if (!adc->anchor_pos) {
		adc->pos_model[X] += user_pt[0] * s->dbip->dbi_local2base;
		adc_model_To_adc_view(s, adc);
		adc_view_To_adc_grid(s, adc);

		adc_set_dirty_flag(s);
	    }

	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc dx' command requires 1 argument\n", (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "dy")) {
	if (argc == 1) {
	    if (!adc->anchor_pos) {
		adc->pos_model[Y] += user_pt[0] * s->dbip->dbi_local2base;
		adc_model_To_adc_view(s, adc);
		adc_view_To_adc_grid(s, adc);

		adc_set_dirty_flag(s);
	    }

	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc dy' command requires 1 argument\n", (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "dz")) {
	if (argc == 1) {
	    if (!adc->anchor_pos) {
		adc->pos_model[Z] += user_pt[0] * s->dbip->dbi_local2base;
		adc_model_To_adc_view(s, adc);
		adc_view_To_adc_grid(s, adc);

		adc_set_dirty_flag(s);
	    }

	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc dz' command requires 1 argument\n", (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "xyz")) {
	if (argc == 0) {
	    VSCALE(scaled_pos, adc->pos_model, s->dbip->dbi_base2local);

	    bu_vls_printf(&vls, "%.15e %.15e %.15e", V3ARGS(scaled_pos));
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 3) {
	    VSCALE(user_pt, user_pt, s->dbip->dbi_local2base);

	    if (incr_flag) {
		VADD2(adc->pos_model, adc->pos_model, user_pt);
	    } else {
		VMOVE(adc->pos_model, user_pt);
	    }

	    adc_model_To_adc_view(s, adc);
	    adc_view_To_adc_grid(s, adc);

	    adc_set_dirty_flag(s);

	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc xyz' command requires 0 or 3 arguments\n", (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "x")) {
	if (argc == 0) {
	    bu_vls_printf(&vls, "%d", adc->dv_x);
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 1) {
	    if (!adc->anchor_pos) {
		if (incr_flag) {
		    adc->dv_x += user_pt[0];
		} else {
		    adc->dv_x = user_pt[0];
		}

		adc->pos_view[X] = adc->dv_x * RT_INV_VIEW;
		adc->pos_view[Y] = adc->dv_y * RT_INV_VIEW;
		adc_view_To_adc_grid(s, adc);
		MAT4X3PNT(adc->pos_model, view2model, adc->pos_view);

		adc_set_dirty_flag(s);
	    }

	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc x' command requires 0 or 1 argument\n", (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "y")) {
	if (argc == 0) {
	    bu_vls_printf(&vls, "%d", adc->dv_y);
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 1) {
	    if (!adc->anchor_pos) {
		if (incr_flag) {
		    adc->dv_y += user_pt[0];
		} else {
		    adc->dv_y = user_pt[0];
		}

		adc->pos_view[X] = adc->dv_x * RT_INV_VIEW;
		adc->pos_view[Y] = adc->dv_y * RT_INV_VIEW;
		adc_view_To_adc_grid(s, adc);
		MAT4X3PNT(adc->pos_model, view2model, adc->pos_view);

		adc_set_dirty_flag(s);
	    }

	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc y' command requires 0 or 1 argument\n", (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "anchor_pos")) {
	if (argc == 0) {
	    bu_vls_printf(&vls, "%d", adc->anchor_pos);
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 1) {
	    i = (int)user_pt[X];

	    if (i < 0 || 2 < i) {
		Tcl_AppendResult(interp, "The 'adc anchor_pos parameter accepts values of 0, 1, or 2.",
				 (char *)NULL);
		ADC_RETURN(TCL_ERROR);
	    }

	    adc->anchor_pos = i;

	    calc_adc_pos(s, adc);
	    adc_set_dirty_flag(s);
	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc anchor_pos' command accepts 0 or 1 argument\n",
			 (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "anchor_a1")) {
	if (argc == 0) {
	    bu_vls_printf(&vls, "%d", adc->anchor_a1);
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 1) {
	    i = (int)user_pt[X];

	    if (i)
		adc->anchor_a1 = 1;
	    else
		adc->anchor_a1 = 0;

	    calc_adc_a1(s, adc);
	    adc_set_dirty_flag(s);
	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc anchor_a1' command accepts 0 or 1 argument\n",
			 (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "anchorpoint_a1")) {
	if (argc == 0) {
	    VSCALE(scaled_pos, adc->anchor_pt_a1, s->dbip->dbi_base2local);

	    bu_vls_printf(&vls, "%.15e %.15e %.15e", V3ARGS(scaled_pos));
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 3) {
	    VSCALE(user_pt, user_pt, s->dbip->dbi_local2base);

	    if (incr_flag) {
		VADD2(adc->anchor_pt_a1, adc->anchor_pt_a1, user_pt);
	    } else {
		VMOVE(adc->anchor_pt_a1, user_pt);
	    }

	    calc_adc_a1(s, adc);
	    adc_set_dirty_flag(s);
	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc anchorpoint_a1' command accepts 0 or 3 arguments\n",
			 (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "anchor_a2")) {
	if (argc == 0) {
	    bu_vls_printf(&vls, "%d", adc->anchor_a2);
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 1) {
	    i = (int)user_pt[X];

	    if (i)
		adc->anchor_a2 = 1;
	    else
		adc->anchor_a2 = 0;

	    calc_adc_a2(s, adc);
	    adc_set_dirty_flag(s);
	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc anchor_a2' command accepts 0 or 1 argument\n",
			 (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "anchorpoint_a2")) {
	if (argc == 0) {
	    VSCALE(scaled_pos, adc->anchor_pt_a2, s->dbip->dbi_base2local);

	    bu_vls_printf(&vls, "%.15e %.15e %.15e", V3ARGS(scaled_pos));
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 3) {
	    VSCALE(user_pt, user_pt, s->dbip->dbi_local2base);

	    if (incr_flag) {
		VADD2(adc->anchor_pt_a2, adc->anchor_pt_a2, user_pt);
	    } else {
		VMOVE(adc->anchor_pt_a2, user_pt);
	    }

	    calc_adc_a2(s, adc);
	    adc_set_dirty_flag(s);
	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc anchorpoint_a2' command accepts 0 or 3 arguments\n",
			 (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "anchor_dst")) {
	if (argc == 0) {
	    bu_vls_printf(&vls, "%d", adc->anchor_dst);
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 1) {
	    i = (int)user_pt[X];

	    if (i) {
		adc->anchor_dst = 1;
	    } else
		adc->anchor_dst = 0;

	    calc_adc_dst(s, adc);
	    adc_set_dirty_flag(s);
	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc anchor_dst' command accepts 0 or 1 argument\n",
			 (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "anchorpoint_dst")) {
	if (argc == 0) {
	    VSCALE(scaled_pos, adc->anchor_pt_dst, s->dbip->dbi_base2local);

	    bu_vls_printf(&vls, "%.15e %.15e %.15e", V3ARGS(scaled_pos));
	    Tcl_AppendResult(interp, bu_vls_addr(&vls), (char *)NULL);
	    bu_vls_free(&vls);

	    ADC_RETURN(TCL_OK);
	} else if (argc == 3) {
	    VSCALE(user_pt, user_pt, s->dbip->dbi_local2base);

	    if (incr_flag) {
		VADD2(adc->anchor_pt_dst, adc->anchor_pt_dst, user_pt);
	    } else {
		VMOVE(adc->anchor_pt_dst, user_pt);
	    }

	    calc_adc_dst(s, adc);
	    adc_set_dirty_flag(s);
	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc anchorpoint_dst' command accepts 0 or 3 arguments\n",
			 (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "reset")) {
	if (argc == 0) {
	    mged_adc_reset(s, adc);

	    adc_set_dirty_flag(s);
	    ADC_RETURN(TCL_OK);
	}

	Tcl_AppendResult(interp, "The 'adc reset' command accepts no arguments\n", (char *)NULL);
	ADC_RETURN(TCL_ERROR);
    }

    if (BU_STR_EQUAL(parameter, "vars")) {
	adc_print_vars(s, adc);
	ADC_RETURN(TCL_OK);
    }

    if (BU_STR_EQUAL(parameter, "help")) {
	Tcl_AppendResult(interp, "Usage:\n", adc_syntax1, adc_syntax2, adc_syntax3, adc_syntax4, (char *)NULL);
	ADC_RETURN(TCL_OK);
    }

    Tcl_AppendResult(interp, "ADC: unrecognized command: '",
		     argv[1], "'\nUsage:\n", adc_syntax1, adc_syntax2, adc_syntax3, adc_syntax4, (char *)NULL);
    ADC_RETURN(TCL_ERROR);
#undef ADC_RETURN
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
