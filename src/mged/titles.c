/*                        T I T L E S . C
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
/** @file mged/titles.c
 *
 */

#include "common.h"

#include <math.h>
#include <string.h>

#include "bv.h"
#include "vmath.h"
#include "bu/units.h"
#include "bn.h"
#include "ged.h"
#include "ged/view.h"
#include "rt/view.h"

#include "./mged.h"
#include "./sedit.h"
#include "./mged_display.h"
#include "./menu.h"
#include "./hud.h"

#define USE_OLD_MENUS 0

char *state_str[] = {
    "-ZOT-",
    "VIEWING",
    "SOL PICK",
    "SOL EDIT",
    "OBJ PICK",
    "OBJ PATH",
    "OBJ EDIT",
    "VERTPICK",
    "UNKNOWN",
};


// FIXME: Global
extern mat_t perspective_mat;  /* defined in dozoom.c */

/*
 * Prepare the numerical display of the currently edited solid/object.
 */
void
create_text_overlay(struct mged_state *s, struct bu_vls *vp)
{
    struct directory *dp;
    struct ged_draw_shape_record hrec;
    int have_highlight = mged_highlight_shape_record(s, &hrec);

    BU_CK_VLS(vp);

    /* Preserve line boundaries for retained HUD label publication. */

    /* print solid info at top of screen
     * Check if the highlighted shape still exists or it has been killed
     * before Accept was clicked.
     */
    if (MEDIT(s) && MEDIT(s)->edit_flag >= 0 && have_highlight && hrec.fullpath && hrec.fullpath->fp_len > 0) {
	dp = DB_FULL_PATH_GET(hrec.fullpath, hrec.fullpath->fp_len - 1);
	if (dp) {

	    bu_vls_strcat(vp, "** SOLID -- ");
	    bu_vls_strcat(vp, dp->d_namep);
	    bu_vls_strcat(vp, ": ");

	    vls_solid(s, vp, MEDIT(s), bn_mat_identity);

	    if (hrec.fullpath->fp_len > 1) {
		bu_vls_strcat(vp, "\n** PATH --  ");
		db_path_to_vls(vp, hrec.fullpath);
		bu_vls_strcat(vp, ": ");

		/* print the evaluated (path) solid parameters */
		vls_solid(s, vp, MEDIT(s), MEDIT(s)->e_mat);
	    }
	}
    }

    /* display path info for object editing also */
    if (s->global_editing_state == ST_O_EDIT && have_highlight && hrec.fullpath) {
	bu_vls_strcat(vp, "** PATH --  ");
	db_path_to_vls(vp, hrec.fullpath);
	bu_vls_strcat(vp, ": ");

	/* print the evaluated (path) solid parameters */
	if (MEDIT(s) && !hrec.evaluated_region) {
	    mat_t new_mat;
	    /* NOT an evaluated region */
	    /* object edit option selected */
	    bn_mat_mul(new_mat, MEDIT(s)->model_changes, MEDIT(s)->e_mat);

	    vls_solid(s, vp, MEDIT(s), new_mat);
	}
    }
    {
	char *start;
	char *p;
	int imax = 0;
	int i = 0;
	int j;
	struct bu_vls vls = BU_VLS_INIT_ZERO;

	start = bu_vls_addr(vp);
	/*
	 * Some display managers don't handle TABs properly, so
	 * we replace any TABs with spaces. Also, look for the
	 * maximum line length.
	 */
	for (p = start; *p != '\0'; ++p) {
	    if (*p == '\t')
		*p = ' ';
	    else if (*p == '\n') {
		if (i > imax)
		    imax = i;
		i = 0;
	    } else
		++i;
	}

	if (i > imax)
	    imax = i;

	/* Prep string for use with Tcl/Tk */
	++imax;
	i = 0;
	for (p = start; *p != '\0'; ++p) {
	    if (*p == '\n') {
		for (j = 0; j < imax - i; ++j)
		    bu_vls_putc(&vls, ' ');

		bu_vls_putc(&vls, *p);
		i = 0;
	    } else {
		bu_vls_putc(&vls, *p);
		++i;
	    }
	}

	Tcl_SetVar(s->interp, "edit_info", bu_vls_addr(&vls), TCL_GLOBAL_ONLY);
	bu_vls_free(&vls);
    }
}


/*
 * Add a vls string to the retained graphics-area text overlay.
 *
 * Set up for character output.  For the best generality, we
 * don't assume that the display can process a CRLF sequence,
 * so each line becomes a separate retained label.
 */
void
screen_vls(
	struct mged_state *s,
	struct mged_hud_builder *hud,
	int xbase,
	int ybase,
	struct bu_vls *vp)
{
    char *start;
    char *end;
    int y;

    BU_CK_VLS(vp);
    y = ybase;

    mged_hud_color_set(hud, color_scheme->cs_edit_info);

    start = bu_vls_addr(vp);
    while (*start != '\0') {
	if ((end = strchr(start, '\n')) == NULL) return;

	*end = '\0';

	(void)mged_hud_label_add(hud, start, GED2PM1(xbase), GED2PM1(y),
		0.0, 0);
	start = end+1;
	y += TEXT0_DY;
    }
}


#define MGED_STATUS_MIN_FONT_SIZE 8.0
#define MGED_STATUS_GLYPH_WIDTH_EM 0.62

/*
 * SoHUDLabel uses a proportional font, but an average glyph width of 0.62 em
 * is a conservative fit estimate for MGED's numeric status strings.  Keep the
 * preferred faceplate font when it fits and reduce it only as much as needed.
 */
static fastf_t
mged_status_font_fit(const char *text, int view_width, fastf_t preferred)
{
    if (!text || !text[0] || view_width <= 0 || preferred <= 0.0)
	return preferred;

    const size_t glyph_count = strlen(text);
    const fastf_t usable_width = view_width > 8 ?
	(fastf_t)(view_width - 8) : (fastf_t)view_width;
    const fastf_t fit = usable_width /
	((fastf_t)glyph_count * MGED_STATUS_GLYPH_WIDTH_EM);
    return fit < preferred ? fit : preferred;
}


static fastf_t
mged_status_font_clamped(const char *text, int view_width, fastf_t preferred)
{
    fastf_t fit = mged_status_font_fit(text, view_width, preferred);
    return fit < MGED_STATUS_MIN_FONT_SIZE ?
	MGED_STATUS_MIN_FONT_SIZE : fit;
}


/*
 * Produce titles, etc., on the screen.
 * NOTE that this routine depends on being called after the retained scene
 * refresh.
 */
void
dotitles(struct mged_state *s, struct bu_vls *overlay_vls)
{
    size_t i = 0;
    struct ged_draw_shape_record hrec;
    int have_highlight = mged_highlight_shape_record(s, &hrec);
    int highlighted_legacy_eval = have_highlight ? hrec.evaluated_region : 0;

    /* for menu computations */
    int x = 0;
    int y = 0;

    int yloc = 0;
    int xloc = 0;
    int scroll_ybot = 0;
    struct bu_vls vls = BU_VLS_INIT_ZERO;

    char cent_x[80] = {0};
    char cent_y[80] = {0};
    char cent_z[80] = {0};
    char size[80] = {0};
    char ang_x[80] = {0};
    char ang_y[80] = {0};
    char ang_z[80] = {0};

    int ss_line_not_drawn=1; /* true if the second status line has not been drawn */
    vect_t temp = VINIT_ZERO;
    fastf_t tmp_val = 0.0;
    mat_t view_center;
    vect_t view_aet = VINIT_ZERO;
    fastf_t view_perspective = 0.0;
    fastf_t view_scale = 1.0;
    fastf_t view_size = 1.0;

    if (s->dbip == DBI_NULL)
	return;

    struct ged_view_context *view_ctx = view_state->vs_gvp;
    const struct bv *view = bv_context_view_const((const struct bv_context *)view_ctx);
    struct bv_params_state params = BV_PARAMS_STATE_INIT;
    (void)bv_params_state_get(&params, view);
    struct mged_hud_builder hud;
    mged_hud_builder_init(&hud, view_ctx, "_faceplate/mged");
    mged_hud_font_size_set(&hud, params.font_size > 0 ? params.font_size : 20);
    mged_hud_line_style_set(&hud, mged_variables->mv_linewidth, 0);
    const int view_width = bv_width_get(view);
    const int view_height = bv_height_get(view);
    const fastf_t view_aspect = (view_width > 0 && view_height > 0) ?
	(fastf_t)view_width / (fastf_t)view_height : 1.0;
    const fastf_t preferred_status_font = params.font_size > 0 ?
	params.font_size : 20.0;
    fastf_t status_font_size = preferred_status_font;
    bv_center_mat_get(view_center, view);
    bv_aet_get(view_aet, view);
    view_perspective = bv_perspective_get(view);
    view_scale = bv_scale_get(view);
    view_size = bv_size_get(view);

    /* Set the Tcl variables to the appropriate values. */

    if (have_highlight && hrec.fullpath) {
	struct bu_vls path_lhs = BU_VLS_INIT_ZERO;
	struct bu_vls path_rhs = BU_VLS_INIT_ZERO;
	struct directory *dp;
	const struct db_full_path *dbfp = hrec.fullpath;

	if (!dbfp) {
	    goto done;
	}
	RT_CK_FULL_PATH(dbfp);

	for (i = 0; i < (size_t)highlight_path_pos; i++) {
	    if ((size_t)i < (size_t)dbfp->fp_len) {
		dp = DB_FULL_PATH_GET(dbfp, i);
		if (dp && dp->d_namep) {
		    bu_vls_printf(&path_lhs, "/%s", dp->d_namep);
		}
	    }
	}
	for (; i < dbfp->fp_len; i++) {
	    dp = DB_FULL_PATH_GET(dbfp, i);
	    if (dp && dp->d_namep) {
		bu_vls_printf(&path_rhs, "/%s", dp->d_namep);
	    }
	}

	bu_vls_printf(&vls, "%s(path_lhs)", MGED_DISPLAY_VAR);
	Tcl_SetVar(s->interp, bu_vls_addr(&vls), bu_vls_addr(&path_lhs), TCL_GLOBAL_ONLY);
	bu_vls_trunc(&vls, 0);
	bu_vls_printf(&vls, "%s(path_rhs)", MGED_DISPLAY_VAR);
	Tcl_SetVar(s->interp, bu_vls_addr(&vls), bu_vls_addr(&path_rhs), TCL_GLOBAL_ONLY);
	bu_vls_free(&path_rhs);
	bu_vls_free(&path_lhs);
    } else {
	bu_vls_printf(&vls, "%s(path_lhs)", MGED_DISPLAY_VAR);
	Tcl_SetVar(s->interp, bu_vls_addr(&vls), "", TCL_GLOBAL_ONLY);
	bu_vls_trunc(&vls, 0);
	bu_vls_printf(&vls, "%s(path_rhs)", MGED_DISPLAY_VAR);
	Tcl_SetVar(s->interp, bu_vls_addr(&vls), "", TCL_GLOBAL_ONLY);
    }

    /* take some care here to avoid buffer overrun */
    tmp_val = -view_center[MDX]*s->dbip->dbi_base2local;
    if (fabs(tmp_val) < 10e70) {
	sprintf(cent_x, "%.3f", tmp_val);
    } else {
	sprintf(cent_x, "%.3g", tmp_val);
    }
    tmp_val = -view_center[MDY]*s->dbip->dbi_base2local;
    if (fabs(tmp_val) < 10e70) {
	sprintf(cent_y, "%.3f", tmp_val);
    } else {
	sprintf(cent_y, "%.3g", tmp_val);
    }
    tmp_val = -view_center[MDZ]*s->dbip->dbi_base2local;
    if (fabs(tmp_val) < 10e70) {
	sprintf(cent_z, "%.3f", tmp_val);
    } else {
	sprintf(cent_z, "%.3g", tmp_val);
    }
    bu_vls_trunc(&vls, 0);
    bu_vls_printf(&vls, "cent=(%s %s %s)", cent_x, cent_y, cent_z);
    Tcl_SetVar(s->interp, bu_vls_addr(&s->mged_curr_display->display_center_name),
	       bu_vls_addr(&vls), TCL_GLOBAL_ONLY);

    tmp_val = view_size*s->dbip->dbi_base2local;
    if (fabs(tmp_val) < 10e70) {
	sprintf(size, "sz=%.3f", tmp_val);
    } else {
	sprintf(size, "sz=%.3g", tmp_val);
    }
    Tcl_SetVar(s->interp, bu_vls_addr(&s->mged_curr_display->display_size_name),
	       size, TCL_GLOBAL_ONLY);

    bu_vls_trunc(&vls, 0);
    bu_vls_printf(&vls, "%s(units)", MGED_DISPLAY_VAR);
    Tcl_SetVar(s->interp, bu_vls_addr(&vls),
	       (char *)bu_units_string(s->dbip->dbi_local2base), TCL_GLOBAL_ONLY);

    bu_vls_trunc(&vls, 0);
    bu_vls_printf(&vls, "az=%3.2f  el=%3.2f  tw=%3.2f", V3ARGS(view_aet));
    Tcl_SetVar(s->interp, bu_vls_addr(&s->mged_curr_display->display_aet_name),
	       bu_vls_addr(&vls), TCL_GLOBAL_ONLY);

    sprintf(ang_x, "%.2f", view_state->k.rot_view[X]);
    sprintf(ang_y, "%.2f", view_state->k.rot_view[Y]);
    sprintf(ang_z, "%.2f", view_state->k.rot_view[Z]);

    bu_vls_trunc(&vls, 0);
    bu_vls_printf(&vls, "ang=(%s %s %s)", ang_x, ang_y, ang_z);
    Tcl_SetVar(s->interp, bu_vls_addr(&s->mged_curr_display->display_ang_name),
	       bu_vls_addr(&vls), TCL_GLOBAL_ONLY);

    /* Label the vertices of the edited solid */
    if (MEDIT(s) &&
	(s->global_editing_state == ST_S_EDIT ||
	 s->global_editing_state == ST_O_EDIT) &&
	(MEDIT(s)->edit_flag >= 0 ||
	 (s->global_editing_state == ST_O_EDIT && !highlighted_legacy_eval))) {
	mat_t xform;
	struct rt_point_labels pl[8+1];
	point_t lines[2*4];	/* up to 4 lines to draw */
	int num_lines=0;

	if (view_perspective <= 0)
	    bn_mat_mul(xform, view_state->vs_model2objview, MEDIT(s)->e_mat);
	else {
	    mat_t tmat;

	    bn_mat_mul(tmat, view_state->vs_model2objview, MEDIT(s)->e_mat);
	    bn_mat_mul(xform, perspective_mat, tmat);
	}

	label_edited_solid(s, &num_lines, lines,  pl, 8+1, xform, &MEDIT(s)->es_int);

	mged_hud_color_set(&hud, color_scheme->cs_geo_label);
	for (i=0; i<(size_t)num_lines; i++)
	    (void)mged_hud_line_add(&hud,
			    GED2PM1(((int)(lines[i*2][X]*RT_VIEW_MAX))),
			    GED2PM1(((int)(lines[i*2][Y]*RT_VIEW_MAX)) * view_aspect),
			    GED2PM1(((int)(lines[i*2+1][X]*RT_VIEW_MAX))),
			    GED2PM1(((int)(lines[i*2+1][Y]*RT_VIEW_MAX)) * view_aspect));
	for (i=0; i<8+1; i++) {
	    if (pl[i].str[0] == '\0') break;
	    (void)mged_hud_label_add(&hud, pl[i].str,
			      GED2PM1(((int)(pl[i].pt[X]*RT_VIEW_MAX))+15),
			      GED2PM1(((int)(pl[i].pt[Y]*RT_VIEW_MAX))+15),
			      0.0, 0);
	}
    }

    if (mged_variables->mv_faceplate) {
	/* Line across the bottom, above two bottom status lines */
	mged_hud_color_set(&hud, color_scheme->cs_other_line);
	(void)mged_hud_line_add(&hud,
			GED2PM1((int)RT_VIEW_MIN), GED2PM1(TITLE_YBASE-TEXT1_DY),
			GED2PM1((int)RT_VIEW_MAX), GED2PM1(TITLE_YBASE-TEXT1_DY));

	if (mged_variables->mv_orig_gui) {
	    /* Enclose window in decorative box.  Mostly for alignment. */
	    (void)mged_hud_line_add(&hud,
			    GED2PM1((int)RT_VIEW_MIN), GED2PM1((int)RT_VIEW_MIN),
			    GED2PM1((int)RT_VIEW_MAX), GED2PM1((int)RT_VIEW_MIN));
	    (void)mged_hud_line_add(&hud,
			    GED2PM1((int)RT_VIEW_MAX), GED2PM1((int)RT_VIEW_MIN),
			    GED2PM1((int)RT_VIEW_MAX), GED2PM1((int)RT_VIEW_MAX));
	    (void)mged_hud_line_add(&hud,
			    GED2PM1((int)RT_VIEW_MAX), GED2PM1((int)RT_VIEW_MAX),
			    GED2PM1((int)RT_VIEW_MIN), GED2PM1((int)RT_VIEW_MAX));
	    (void)mged_hud_line_add(&hud,
			    GED2PM1((int)RT_VIEW_MIN), GED2PM1((int)RT_VIEW_MAX),
			    GED2PM1((int)RT_VIEW_MIN), GED2PM1((int)RT_VIEW_MIN));

	    /* Display scroll bars */
	    scroll_ybot = scroll_display(s, &hud, SCROLLY);
	    y = MENUY;
	    x = MENUX;

	    /* Display state and local unit in upper left corner, boxed */
	    mged_hud_color_set(&hud, color_scheme->cs_state_text1);
	    (void)mged_hud_label_add(&hud, state_str[s->global_editing_state],
		    GED2PM1(MENUX), GED2PM1(MENUY - MENU_DY),
		    0.0, 0);
	} else {
	    scroll_ybot = SCROLLY;
	    x = (int)RT_VIEW_MIN + 20;
	    y = (int)RT_VIEW_MAX+TEXT0_DY;
	}

	/*
	 * Print information about object illuminated
	 */
	if (have_highlight && hrec.fullpath &&
	    (s->global_editing_state == ST_O_PATH || s->global_editing_state==ST_O_PICK || s->global_editing_state==ST_S_PICK)) {

	    for (i=0; i < hrec.fullpath->fp_len; i++) {
		if (i == (size_t)highlight_path_pos  &&  s->global_editing_state == ST_O_PATH) {
		    mged_hud_color_set(&hud, color_scheme->cs_state_text1);
		    (void)mged_hud_label_add(&hud, "[MATRIX]",
			    GED2PM1(x), GED2PM1(y), 0.0, 0);
		    y += MENU_DY;
		}
		mged_hud_color_set(&hud, color_scheme->cs_state_text2);
		(void)mged_hud_label_add(&hud,
				  DB_FULL_PATH_GET(hrec.fullpath, i)->d_namep,
				  GED2PM1(x), GED2PM1(y), 0.0, 0);
		y += MENU_DY;
	    }
	}

	if (mged_variables->mv_orig_gui) {
	    mged_hud_color_set(&hud, color_scheme->cs_other_line);
	    (void)mged_hud_line_add(&hud,
			    GED2PM1(MENUXLIM), GED2PM1(y),
			    GED2PM1(MENUXLIM), GED2PM1((int)RT_VIEW_MAX));	/* vert. */
	    /*
	     * The top of the menu (if any) begins at the Y value specified.
	     */
	    mmenu_display(s, &hud, y);

	    /* print parameter locations on screen */
	    if (s->global_editing_state == ST_O_EDIT && highlighted_legacy_eval) {
		/* region is a processed region */
		MAT4X3PNT(temp, view_state->vs_model2objview, MEDIT(s)->e_keypoint);
		xloc = (int)(temp[X]*RT_VIEW_MAX);
		yloc = (int)(temp[Y]*RT_VIEW_MAX);
		mged_hud_color_set(&hud, color_scheme->cs_edit_info);
		(void)mged_hud_line_add(&hud,
				GED2PM1(xloc-TEXT0_DY), GED2PM1(yloc+TEXT0_DY),
				GED2PM1(xloc+TEXT0_DY), GED2PM1(yloc-TEXT0_DY));
		(void)mged_hud_line_add(&hud,
				GED2PM1(xloc-TEXT0_DY), GED2PM1(yloc-TEXT0_DY),
				GED2PM1(xloc+TEXT0_DY), GED2PM1(yloc+TEXT0_DY));
		(void)mged_hud_line_add(&hud,
				GED2PM1(xloc+TEXT0_DY), GED2PM1(yloc+TEXT0_DY),
				GED2PM1(xloc-TEXT0_DY), GED2PM1(yloc+TEXT0_DY));
		(void)mged_hud_line_add(&hud,
				GED2PM1(xloc+TEXT0_DY), GED2PM1(yloc-TEXT0_DY),
				GED2PM1(xloc-TEXT0_DY), GED2PM1(yloc-TEXT0_DY));
		(void)mged_hud_line_add(&hud,
				GED2PM1(xloc+TEXT0_DY), GED2PM1(yloc+TEXT0_DY),
				GED2PM1(xloc+TEXT0_DY), GED2PM1(yloc-TEXT0_DY));
		(void)mged_hud_line_add(&hud,
				GED2PM1(xloc-TEXT0_DY), GED2PM1(yloc+TEXT0_DY),
				GED2PM1(xloc-TEXT0_DY), GED2PM1(yloc-TEXT0_DY));
	    }
	}

	/*
	 * Prepare the numerical display of the currently edited solid/object.
	 */
	/* create_text_overlay(s, &vls); */
	if (mged_variables->mv_orig_gui) {
	    screen_vls(s, &hud, SOLID_XBASE, scroll_ybot+TEXT0_DY, overlay_vls);
	} else {
	    screen_vls(s, &hud, x, y, overlay_vls);
	}

	/*
	 * General status information on first status line
	 */
	bu_vls_trunc(&vls, 0);
	bu_vls_printf(&vls,
		      " cent=(%s, %s, %s), %s %s, ", cent_x, cent_y, cent_z,
		      size, bu_units_string(s->dbip->dbi_local2base));
	bu_vls_printf(&vls, "az=%3.2f el=%3.2f tw=%3.2f ang=(%s, %s, %s)", V3ARGS(view_aet),
		      ang_x, ang_y, ang_z);
	status_font_size = mged_status_font_fit(bu_vls_cstr(&vls),
		view_width, preferred_status_font);
	if (status_font_size < MGED_STATUS_MIN_FONT_SIZE) {
	    /* Rotation deltas are the least important telemetry at narrow widths. */
	    bu_vls_sprintf(&vls,
		    " cent=(%s, %s, %s), %s %s, az=%3.2f el=%3.2f tw=%3.2f",
		    cent_x, cent_y, cent_z, size,
		    bu_units_string(s->dbip->dbi_local2base), V3ARGS(view_aet));
	    status_font_size = mged_status_font_fit(bu_vls_cstr(&vls),
		    view_width, preferred_status_font);
	}
	if (status_font_size < MGED_STATUS_MIN_FONT_SIZE) {
	    /* Preserve scale and orientation when the center no longer fits. */
	    bu_vls_sprintf(&vls, " %s %s, az=%3.2f el=%3.2f tw=%3.2f",
		    size, bu_units_string(s->dbip->dbi_local2base),
		    V3ARGS(view_aet));
	    status_font_size = mged_status_font_fit(bu_vls_cstr(&vls),
		    view_width, preferred_status_font);
	}
	if (status_font_size < MGED_STATUS_MIN_FONT_SIZE) {
	    /* A readable view size is still useful in an extremely narrow view. */
	    bu_vls_sprintf(&vls, " %s %s", size,
		    bu_units_string(s->dbip->dbi_local2base));
	}
	status_font_size = mged_status_font_clamped(bu_vls_cstr(&vls),
		view_width, preferred_status_font);
	mged_hud_color_set(&hud, color_scheme->cs_status_text1);
	(void)mged_hud_label_add(&hud, bu_vls_addr(&vls),
		GED2PM1(TITLE_XBASE), GED2PM1(TITLE_YBASE),
		status_font_size, 0);
    } /* if faceplate !0 */

    /*
     * Second status line
     */

    /* Priorities for what to display:
     * 1.  adc info
     * 2.  keypoint
     * 3.  illuminated path
     *
     * This way the adc info will be displayed during editing
     */

    struct bv_adc_state adc = {0};
    (void)mged_display_adc_state_get(s->mged_curr_display, &adc);
    if (adc.draw) {
	fastf_t f;

	f = view_scale * s->dbip->dbi_base2local;
	/* Angle/Distance cursor */
	bu_vls_trunc(&vls, 0);
	bu_vls_printf(&vls,
		      " curs:  a1=%.1f,  a2=%.1f,  dst=%.3f,  cent=(%.3f, %.3f),  delta=(%.3f, %.3f)",
		      adc.a1, adc.a2,
		      adc.dst * f,
		      adc.pos_grid[X] * f, adc.pos_grid[Y] * f,
		      adc.pos_view[X] * f, adc.pos_view[Y] * f);
	if (mged_variables->mv_faceplate) {
	    mged_hud_color_set(&hud, color_scheme->cs_status_text2);
	    (void)mged_hud_label_add(&hud, bu_vls_addr(&vls),
		    GED2PM1(TITLE_XBASE), GED2PM1(TITLE_YBASE + TEXT1_DY),
		    mged_status_font_clamped(bu_vls_cstr(&vls), view_width,
			status_font_size), 0);
	}
	Tcl_SetVar(s->interp, bu_vls_addr(&s->mged_curr_display->display_adc_name),
		   bu_vls_addr(&vls), TCL_GLOBAL_ONLY);
	ss_line_not_drawn = 0;
    } else {
	Tcl_SetVar(s->interp, bu_vls_addr(&s->mged_curr_display->display_adc_name), "", TCL_GLOBAL_ONLY);
    }

    if (s->global_editing_state == ST_S_EDIT || s->global_editing_state == ST_O_EDIT) {
	struct bu_vls kp_vls = BU_VLS_INIT_ZERO;

	bu_vls_printf(&kp_vls,
		      " Keypoint: %s %s: (%g, %g, %g)",
		      OBJ[MEDIT(s)->es_int.idb_type].ft_name+3,	/* Skip ID_ */
		      MEDIT(s)->e_keytag,
		      MEDIT(s)->e_keypoint[X] * s->dbip->dbi_base2local,
		      MEDIT(s)->e_keypoint[Y] * s->dbip->dbi_base2local,
		      MEDIT(s)->e_keypoint[Z] * s->dbip->dbi_base2local);
	if (mged_variables->mv_faceplate && ss_line_not_drawn) {
	    mged_hud_color_set(&hud, color_scheme->cs_status_text2);
	    (void)mged_hud_label_add(&hud, bu_vls_addr(&kp_vls),
		    GED2PM1(TITLE_XBASE), GED2PM1(TITLE_YBASE + TEXT1_DY),
		    mged_status_font_clamped(bu_vls_cstr(&kp_vls), view_width,
			status_font_size), 0);
	    ss_line_not_drawn = 0;
	}

	bu_vls_trunc(&vls, 0);
	bu_vls_printf(&vls, "%s(keypoint)", MGED_DISPLAY_VAR);
	Tcl_SetVar(s->interp, bu_vls_addr(&vls), bu_vls_addr(&kp_vls), TCL_GLOBAL_ONLY);

	bu_vls_free(&kp_vls);
    } else {
	bu_vls_trunc(&vls, 0);
	bu_vls_printf(&vls, "%s(keypoint)", MGED_DISPLAY_VAR);
	Tcl_SetVar(s->interp, bu_vls_addr(&vls), "", TCL_GLOBAL_ONLY);
    }

    if (have_highlight && hrec.fullpath) {
	if (mged_variables->mv_faceplate && ss_line_not_drawn) {
	    bu_vls_trunc(&vls, 0);

	    /* Illuminated path */
	    bu_vls_strcat(&vls, " Path: ");
	    for (i=0; i < hrec.fullpath->fp_len; i++) {
		if (i == (size_t)highlight_path_pos  &&
		    (s->global_editing_state == ST_O_PATH || s->global_editing_state == ST_O_EDIT))
		    bu_vls_strcat(&vls, "/__MATRIX__");
		bu_vls_printf(&vls, "/%s",
			      DB_FULL_PATH_GET(hrec.fullpath, i)->d_namep);
	    }
	    mged_hud_color_set(&hud, color_scheme->cs_status_text2);
	    (void)mged_hud_label_add(&hud, bu_vls_addr(&vls),
		    GED2PM1(TITLE_XBASE), GED2PM1(TITLE_YBASE + TEXT1_DY),
		    mged_status_font_clamped(bu_vls_cstr(&vls), view_width,
			status_font_size), 0);

	    ss_line_not_drawn = 0;
	}
    }

    bu_vls_trunc(&vls, 0);
    bu_vls_printf(&vls, "%.2f fps", 1/frametime);
    if (mged_variables->mv_faceplate && ss_line_not_drawn) {
	mged_hud_color_set(&hud, color_scheme->cs_status_text2);
	(void)mged_hud_label_add(&hud, bu_vls_addr(&vls),
		GED2PM1(TITLE_XBASE), GED2PM1(TITLE_YBASE + TEXT1_DY),
		mged_status_font_clamped(bu_vls_cstr(&vls), view_width,
		    status_font_size), 0);
    }
    Tcl_SetVar(s->interp, bu_vls_addr(&s->mged_curr_display->display_fps_name),
	       bu_vls_addr(&vls), TCL_GLOBAL_ONLY);

done:
    (void)mged_hud_builder_publish(&hud);
    mged_hud_builder_free(&hud);
    bu_vls_free(&vls);
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
