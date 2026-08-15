/*                         S H A R E . C
 * BRL-CAD
 *
 * Copyright (c) 1998-2026 United States Government as represented by
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
/** @file mged/share.c
 *
 * Description -
 * Routines for sharing resources among display managers.
 *
 */

#include "common.h"

#include <math.h>
#include <string.h>

#include "vmath.h"
#include "bn.h"

#include "./mged.h"
#include "./mged_display.h"

#define RESOURCE_TYPE_ADC		0
#define RESOURCE_TYPE_AXES		1
#define RESOURCE_TYPE_COLOR_SCHEMES	2
#define RESOURCE_TYPE_GRID		3
#define RESOURCE_TYPE_MENU		4
#define RESOURCE_TYPE_MGED_VARIABLES	5
#define RESOURCE_TYPE_RUBBER_BAND	6
#define RESOURCE_TYPE_VIEW		7

#define SHARE_RESOURCE(uflag, str, resource, rc, dlp1, dlp2, vls, error_msg) \
    do { \
	if (uflag) { \
	    struct str *strp; \
\
	    if (dlp1->resource->rc > 1) {   /* must be sharing this resource */ \
		--dlp1->resource->rc; \
		strp = dlp1->resource; \
		BU_ALLOC(dlp1->resource, struct str); \
		*dlp1->resource = *strp;        /* struct copy */ \
		dlp1->resource->rc = 1; \
	    } \
	} else { \
	    /* must not be sharing this resource */ \
	    if (dlp1->resource != dlp2->resource) { \
		if (!--dlp2->resource->rc) \
		    bu_free((void *)dlp2->resource, error_msg); \
\
		dlp2->resource = dlp1->resource; \
		++dlp1->resource->rc; \
	    } \
	} \
    } while (0)


/* These vparse arrays are defined in set.c and cmd.c.  They are exposed as
 * extern declarations here because share.c needs to iterate over them when
 * sharing/unsharing resources between display managers.  Moving them to a
 * header is a future clean-up task; for now the explicit externs suffice. */
extern struct bu_structparse axes_vparse[];
extern struct bu_structparse color_scheme_vparse[];
extern struct bu_structparse grid_vparse[];
extern struct bu_structparse rubber_band_vparse[];
extern struct bu_structparse mged_vparse[];

void free_all_resources(struct mged_display *dlp);

static void
copy_faceplate_record(int uflag, struct mged_display *src, struct mged_display *dst, int is_adc)
{
    if (uflag || !src || !dst)
	return;

    if (is_adc) {
	struct bv_adc_state adc;
	if (mged_display_adc_state_get(src, &adc))
	    mged_display_adc_state_set(dst, &adc);
    } else {
	struct bv_grid_state grid;
	if (mged_display_grid_state_get(src, &grid))
	    mged_display_grid_state_set(dst, &grid);
    }

    mged_display_repaint_request(src, MGED_REPAINT_DEVICE_SETTING);
    mged_display_repaint_request(dst, MGED_REPAINT_DEVICE_SETTING);
}

/*
 * SYNOPSIS
 *	share [-u] res p1 [p2]
 *
 * DESCRIPTION
 *	Provides a mechanism to (un)share resources among display managers.
 *	Currently, there are nine different resources that can be shared.
 *	They are:
 *		ADC AXES COLOR_SCHEMES GRID MENU MGED_VARIABLES RUBBER_BAND VIEW
 *
 * EXAMPLES
 *	share res_type p1 p2	--->	causes 'p1' to share its resource of type 'res_type' with 'p2'
 *	share -u res_type p	--->	causes 'p' to no longer share resource of type 'res_type'
 */
int
f_share(ClientData clientData, Tcl_Interp *interpreter, int argc, const char *argv[])
{
    struct cmdtab *ctp = (struct cmdtab *)clientData;
    MGED_CK_CMD(ctp);

    int uflag = 0;		/* unshare flag */
    struct mged_display *dlp1 = MGED_DISPLAY_NULL;
    struct mged_display *dlp2 = MGED_DISPLAY_NULL;
    struct bu_vls vls = BU_VLS_INIT_ZERO;

    if (argc != 4) {
	bu_vls_printf(&vls, "helpdevel share");
	Tcl_Eval(interpreter, bu_vls_addr(&vls));

	bu_vls_free(&vls);
	return TCL_ERROR;
    }

    if (argv[1][0] == '-' && argv[1][1] == 'u') {
	uflag = 1;
	--argc;
	++argv;
    }

    for (size_t di = 0; di < BU_PTBL_LEN(&active_display_set); di++) {
	struct mged_display *m_dmp = (struct mged_display *)BU_PTBL_GET(&active_display_set, di);
	struct bu_vls *pname = mged_display_pathname_vls(m_dmp);
	if (BU_STR_EQUAL(argv[2], bu_vls_cstr(pname))) {
	    dlp1 = m_dmp;
	    break;
	}
    }

    if (dlp1 == MGED_DISPLAY_NULL) {
	Tcl_AppendResult(interpreter, "share: unrecognized path name - ", argv[2], "\n", (char *)NULL);
	bu_vls_free(&vls);
	return TCL_ERROR;
    }

    if (!uflag) {
	for (size_t di = 0; di < BU_PTBL_LEN(&active_display_set); di++) {
	    struct mged_display *m_dmp = (struct mged_display *)BU_PTBL_GET(&active_display_set, di);
	    struct bu_vls *pname = mged_display_pathname_vls(m_dmp);
	    if (BU_STR_EQUAL(argv[3], bu_vls_cstr(pname))) {
		dlp2 = m_dmp;
		break;
	    }
	}

	if (dlp2 == MGED_DISPLAY_NULL) {
	    Tcl_AppendResult(interpreter, "share: unrecognized path name - ", argv[3], "\n", (char *)NULL);
	    bu_vls_free(&vls);
	    return TCL_ERROR;
	}

	/* same graphics view */
	if (dlp1 == dlp2) {
	    bu_vls_free(&vls);
	    return TCL_OK;
	}
    }

    switch (argv[1][0]) {
	case 'a':
	case 'A':
	    if (argv[1][1] == 'd' || argv[1][1] == 'D')
		copy_faceplate_record(uflag, dlp1, dlp2, 1);
	    else if (argv[1][1] == 'x' || argv[1][1] == 'X') {
		SHARE_RESOURCE(uflag, _axes_state, display_axes_state, ax_rc, dlp1, dlp2, vls, "share: axes_state");
		if (!uflag)
		    dlp2->display_axes_state_dirty = 1;
	    }
	    else {
		bu_vls_printf(&vls, "share: resource type '%s' unknown\n", argv[1]);
		Tcl_AppendResult(interpreter, bu_vls_addr(&vls), (char *)NULL);

		bu_vls_free(&vls);
		return TCL_ERROR;
	    }
	    break;
	case 'c':
	case 'C':
	    SHARE_RESOURCE(uflag, _color_scheme, display_color_scheme, cs_rc, dlp1, dlp2, vls, "share: color_scheme");
	    if (!uflag) {
		dlp2->display_color_scheme_dirty = 1;
		dlp2->display_axes_state_dirty = 1;
		dlp2->display_adc_style_dirty = 1;
	    }
	    break;
	case 'g':
	case 'G':
	    copy_faceplate_record(uflag, dlp1, dlp2, 0);
	    break;
	case 'm':
	case 'M':
	    SHARE_RESOURCE(uflag, _menu_state, display_menu_state, ms_rc, dlp1, dlp2, vls, "share: menu_state");
	    break;
	case 'r':
	case 'R':
	    SHARE_RESOURCE(uflag, _rubber_band, display_rubber_band, rb_rc, dlp1, dlp2, vls, "share: rubber_band");
	    break;
	case 'v':
	case 'V':
	    if ((argv[1][1] == 'a' || argv[1][1] == 'A') &&
		(argv[1][2] == 'r' || argv[1][2] == 'R')) {
		SHARE_RESOURCE(uflag, _mged_variables, display_variables, mv_rc, dlp1, dlp2, vls, "share: mged_variables");
		if (!uflag)
		    dlp2->display_framebuffer_state_dirty = 1;
	    }
	    else if (argv[1][1] == 'i' || argv[1][1] == 'I') {
		if (!uflag) {
		    /* free dlp2's view_state resources if currently not sharing */
		    if (dlp2->display_view_state->vs_rc == 1)
			view_ring_destroy(dlp2);
		}

		SHARE_RESOURCE(uflag, _view_state, display_view_state, vs_rc, dlp1, dlp2, vls, "share: view_state");

		if (uflag) {
		    struct _view_state *ovsp;
		    ovsp = dlp1->display_view_state;

		    /* initialize dlp1's view_state */
		    if (ovsp != dlp1->display_view_state)
			view_ring_init(dlp1->display_view_state, ovsp);
		}
	    } else {
		bu_vls_printf(&vls, "share: resource type '%s' unknown\n", argv[1]);
		Tcl_AppendResult(interpreter, bu_vls_addr(&vls), (char *)NULL);

		bu_vls_free(&vls);
		return TCL_ERROR;
	    }

	    break;
	default:
	    bu_vls_printf(&vls, "share: resource type '%s' unknown\n", argv[1]);
	    Tcl_AppendResult(interpreter, bu_vls_addr(&vls), (char *)NULL);

	    bu_vls_free(&vls);
	    return TCL_ERROR;
    }

    if (!uflag) {
	mged_display_repaint_request(dlp2, MGED_REPAINT_DEVICE_SETTING);
    }

    bu_vls_free(&vls);
    return TCL_OK;
}


/*
 * SYNOPSIS
 *	rset [res_type [res [vals]]]
 *
 * DESCRIPTION
 *	Provides a mechanism to set resource values for some resource.
 *
 * EXAMPLES
 *	rset c bg 0 0 50	--->	sets the background color to dark blue
 */
// TODO - is e_type actually used here?
int
f_rset (ClientData clientData, Tcl_Interp *interpreter, int argc, const char *argv[])
{
    struct cmdtab *ctp = (struct cmdtab *)clientData;
    MGED_CK_CMD(ctp);
    struct mged_state *s = ctp->s;

    struct bv_grid_state grid = {0};
    struct bu_vls vls = BU_VLS_INIT_ZERO;
    (void)mged_display_grid_state_get(s->mged_curr_display, &grid);

    /* print values for all resources */
    if (argc == 1) {
	mged_vls_struct_parse(s, &vls, "Axes, res_type - ax", axes_vparse,
			      (const char *)axes_state, argc, argv);
	bu_vls_printf(&vls, "\n");
	mged_vls_struct_parse(s, &vls, "Color Schemes, res_type - c", color_scheme_vparse,
			      (const char *)color_scheme, argc, argv);
	bu_vls_printf(&vls, "\n");
	mged_vls_struct_parse(s, &vls, "Grid, res_type - g", grid_vparse,
			      (const char *)&grid, argc, argv);
	bu_vls_printf(&vls, "\n");
	mged_vls_struct_parse(s, &vls, "Rubber Band, res_type - r", rubber_band_vparse,
			      (const char *)rubber_band, argc, argv);
	bu_vls_printf(&vls, "\n");
	mged_vls_struct_parse(s, &vls, "MGED Variables, res_type - var", mged_vparse,
			      (const char *)mged_variables, argc, argv);

	Tcl_AppendResult(interpreter, bu_vls_addr(&vls), (char *)NULL);
	bu_vls_free(&vls);

	return TCL_OK;
    }

    switch (argv[1][0]) {
	case 'a':
	case 'A':
	    if (argv[1][1] == 'd' || argv[1][1] == 'D')
		bu_vls_printf(&vls, "rset: use the adc command for the 'adc' resource");
	    else if (argv[1][1] == 'x' || argv[1][1] == 'X')
		mged_vls_struct_parse(s, &vls, "Axes", axes_vparse,
				      (const char *)axes_state, argc-1, argv+1);
	    else {
		bu_vls_printf(&vls, "rset: resource type '%s' unknown\n", argv[1]);
		Tcl_AppendResult(interpreter, bu_vls_addr(&vls), (char *)NULL);

		bu_vls_free(&vls);
		return TCL_ERROR;
	    }
	    break;
	case 'c':
	case 'C':
	    mged_vls_struct_parse(s, &vls, "Color Schemes", color_scheme_vparse,
			  (const char *)color_scheme, argc-1, argv+1);
	    break;
	case 'g':
	case 'G':
	    mged_vls_struct_parse(s, &vls, "Grid", grid_vparse,
				  (const char *)&grid, argc-1, argv+1);
	    mged_display_grid_state_set(s->mged_curr_display, &grid);
	    break;
	case 'r':
	case 'R':
	    mged_vls_struct_parse(s, &vls, "Rubber Band", rubber_band_vparse,
				  (const char *)rubber_band, argc-1, argv+1);
	    break;
	case 'v':
	case 'V':
	    if ((argv[1][1] == 'a' || argv[1][1] == 'A') &&
		(argv[1][2] == 'r' || argv[1][2] == 'R'))
		mged_vls_struct_parse(s, &vls, "mged variables", mged_vparse,
				      (const char *)mged_variables, argc-1, argv+1);
	    else if (argv[1][1] == 'i' || argv[1][1] == 'I')
		bu_vls_printf(&vls, "rset: no support available for the 'view' resource");
	    else {
		bu_vls_printf(&vls, "rset: resource type '%s' unknown\n", argv[1]);
		Tcl_AppendResult(interpreter, bu_vls_addr(&vls), (char *)NULL);

		bu_vls_free(&vls);
		return TCL_ERROR;
	    }

	    break;
	default:
	    bu_vls_printf(&vls, "rset: resource type '%s' unknown\n", argv[1]);
	    Tcl_AppendResult(interpreter, bu_vls_addr(&vls), (char *)NULL);

	    bu_vls_free(&vls);
	    return TCL_ERROR;
    }

    Tcl_AppendResult(interpreter, bu_vls_addr(&vls), (char *)NULL);
    bu_vls_free(&vls);

    return TCL_OK;
}


/*
 * dlp1 takes control of dlp2's resources. dlp2 is
 * probably on its way out (i.e. being destroyed).
 */
void
usurp_all_resources(struct mged_display *dlp1, struct mged_display *dlp2)
{
    free_all_resources(dlp1);
    dlp1->display_view_state = dlp2->display_view_state;
    dlp1->display_menu_state = dlp2->display_menu_state;
    dlp1->display_rubber_band = dlp2->display_rubber_band;
    dlp1->display_variables = dlp2->display_variables;
    dlp1->display_color_scheme = dlp2->display_color_scheme;
    dlp1->display_color_scheme_dirty = 1;
    dlp1->display_axes_state = dlp2->display_axes_state;
    dlp1->display_axes_state_dirty = 1;
    dlp1->display_adc_style_dirty = 1;
    dlp1->display_framebuffer_state_dirty = 1;

    /* sanity */
    dlp2->display_view_state = (struct _view_state *)NULL;
    dlp2->display_menu_state = (struct _menu_state *)NULL;
    dlp2->display_rubber_band = (struct _rubber_band *)NULL;
    dlp2->display_variables = (struct _mged_variables *)NULL;
    dlp2->display_color_scheme = (struct _color_scheme *)NULL;
    dlp2->display_axes_state = (struct _axes_state *)NULL;
}


/*
 * - decrement the reference count of all resources
 * - free all resources that are not being used
 */
void
free_all_resources(struct mged_display *dlp)
{
    if (!--dlp->display_view_state->vs_rc) {
	view_ring_destroy(dlp);
	bu_free((void *)dlp->display_view_state, "free_all_resources: view_state");
    }

    if (!--dlp->display_menu_state->ms_rc)
	bu_free((void *)dlp->display_menu_state, "free_all_resources: menu_state");

    if (!--dlp->display_rubber_band->rb_rc)
	bu_free((void *)dlp->display_rubber_band, "free_all_resources: rubber_band");

    if (!--dlp->display_variables->mv_rc)
	bu_free((void *)dlp->display_variables, "free_all_resources: mged_variables");

    if (!--dlp->display_color_scheme->cs_rc)
	bu_free((void *)dlp->display_color_scheme, "free_all_resources: color_scheme");

    if (!--dlp->display_axes_state->ax_rc)
	bu_free((void *)dlp->display_axes_state, "free_all_resources: axes_state");
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
