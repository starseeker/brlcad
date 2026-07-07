/*                       A R R O W S . C
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
/** @file libtclcad/view/arrows.c
 *
 */
/** @} */

#include "common.h"
#include "bu/units.h"
#include "ged.h"
#include "ged/view.h"
#include "tclcad.h"

/* Private headers */
#include "ged/draw.h"
#include "../tclcad_private.h"
#include "../view/view.h"

/* The "view get" introspection path (getters in to_data_arrows_func) recovers
 * values through typed GED draw-view data-arrow facades, making the draw view
 * the canonical read source for Tcl introspection.
 *
 * Setters no longer write gv_tcl at all; they request typed data-arrow facade
 * operations for color, line_width, tip_length, tip_width, draw, and points. */

int
go_data_arrows(Tcl_Interp *interp,
	       struct ged *gedp,
	       void *draw_view_ctx,
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

    ret = to_data_arrows_func(interp, gedp, draw_view_ctx, argc, argv);
    to_refresh_suppress_all_end(current_top);
    if (ret & BRLCAD_ERROR)
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);

    return ret;
}


int
to_data_arrows(struct ged *gedp,
	       int argc,
	       const char *argv[],
	       ged_func_ptr UNUSED(func),
	       const char *usage,
	       int UNUSED(maxargs))
{
    void *view_ctx;
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

    /* shift the command name to argv[1] before calling to_data_arrows_func */
    argv[1] = argv[0];
    ret = to_data_arrows_func(current_top->to_interp, gedp, view_ctx, argc-1, argv+1);
    if (ret == BRLCAD_ERROR)
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);

    return ret;
}


int
to_data_arrows_func(Tcl_Interp *interp,
		    struct ged *gedp,
		    void *view_ctx,
		    int argc,
		    const char *argv[])
{
    /* T3: BSG object name is the only per-variant state needed here. */
    const char *bsg_name = (argv[0][0] == 's') ? "_tcl_sdata_arrows" : "_tcl_data_arrows";

    if (BU_STR_EQUAL(argv[1], "draw")) {
	if (argc == 2) {
	    bu_vls_printf(gedp->ged_result_str, "%d",
			  ged_draw_view_context_data_arrows_draw_get(view_ctx, bsg_name) ? 1 : 0);
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int i;

	    if (bu_sscanf(argv[2], "%d", &i) != 1)
		goto bad;

	    ged_draw_view_context_data_arrows_draw_set(view_ctx, bsg_name, i ? 1 : 0);
	    /* If no BSG object exists and draw=1 is requested, nothing to show
	     * yet (no points have been set); silently no-op. */

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "color")) {
	if (argc == 2) {
	    struct ged_draw_view_feature_style style = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	    if (ged_draw_view_context_data_arrows_style_get(view_ctx, bsg_name, &style) && style.color_valid) {
		bu_vls_printf(gedp->ged_result_str, "%d %d %d",
			      (int)style.color[0], (int)style.color[1], (int)style.color[2]);
	    } else {
		bu_vls_printf(gedp->ged_result_str, "0 0 0");
	    }
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

	    ged_draw_view_context_data_arrows_color_set(view_ctx, bsg_name, r, g, b);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "line_width")) {
	if (argc == 2) {
	    struct ged_draw_view_feature_style style = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	    if (ged_draw_view_context_data_arrows_style_get(view_ctx, bsg_name, &style))
		bu_vls_printf(gedp->ged_result_str, "%d", style.line_width);
	    else
		bu_vls_printf(gedp->ged_result_str, "0");
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int line_width;

	    if (bu_sscanf(argv[2], "%d", &line_width) != 1)
		goto bad;

	    ged_draw_view_context_data_arrows_line_width_set(view_ctx, bsg_name, line_width);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "points")) {
	register int i;

	if (argc == 2) {
	    point_t *pts = NULL;
	    size_t npts = 0;
	    if (ged_draw_view_context_data_arrows_points_copy(view_ctx, bsg_name, &pts, &npts)) {
		for (size_t _j = 0; _j < npts; _j++)
		    bu_vls_printf(gedp->ged_result_str, " {%lf %lf %lf} ", V3ARGS(pts[_j]));
		bu_free(pts, "GED draw view feature points copy");
	    }
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int ac;
	    const char **av;

	    if (Tcl_SplitList(interp, argv[2], &ac, &av) != TCL_OK) {
		bu_vls_printf(gedp->ged_result_str, "%s", Tcl_GetStringResult(interp));
		return BRLCAD_ERROR;
	    }

	    if (ac % 2) {
		bu_vls_printf(gedp->ged_result_str, "%s: must be an even number of points", argv[0]);
		Tcl_Free((char *)av);
		return BRLCAD_ERROR;
	    }

	    /* Save style from existing object before replacing it. */
	    struct ged_draw_view_feature_style saved_style = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	    ged_draw_view_context_data_arrows_style_get(view_ctx, bsg_name, &saved_style);

	    /* Clear out: remove old draw-view feature. */
	    if (ac < 2) {
		ged_draw_view_context_data_arrows_points_replace(view_ctx,
			bsg_name, NULL, 0, &saved_style);
		Tcl_Free((char *)av);
		to_refresh_view(view_ctx);
		return BRLCAD_OK;
	    }

	    /* Parse points into temporary local array. */
	    point_t *pts = (point_t *)bu_calloc(ac, sizeof(point_t), "arrow points");
	    for (i = 0; i < ac; ++i) {
		double scan[ELEMENTS_PER_VECT];

		if (bu_sscanf(av[i], "%lf %lf %lf", &scan[X], &scan[Y], &scan[Z]) != 3) {
		    bu_vls_printf(gedp->ged_result_str, "bad data point - %s\n", av[i]);
		    bu_free(pts, "arrow points");
		    Tcl_Free((char *)av);
		    return BRLCAD_ERROR;
		}
		VMOVE(pts[i], scan);
	    }

	    /* Rebuild draw-view data arrows from new points, preserving style. */
	    (void)ged_draw_view_context_data_arrows_points_replace(view_ctx, bsg_name, pts, (size_t)ac,
		    &saved_style);
	    bu_free(pts, "arrow points");
	    Tcl_Free((char *)av);
	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}
    }

    if (BU_STR_EQUAL(argv[1], "tip_length")) {
	if (argc == 2) {
	    fastf_t tip_length = 0.0;
	    if (ged_draw_view_context_data_arrows_tip_get(view_ctx, bsg_name, &tip_length, NULL)) {
		bu_vls_printf(gedp->ged_result_str, "%d", (int)tip_length);
	    } else {
		bu_vls_printf(gedp->ged_result_str, "0");
	    }
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int tip_length;

	    if (bu_sscanf(argv[2], "%d", &tip_length) != 1)
		goto bad;

	    fastf_t tip_width = 0.0;
	    if (ged_draw_view_context_data_arrows_tip_get(view_ctx, bsg_name, NULL, &tip_width))
		ged_draw_view_context_data_arrows_tip_set(view_ctx, bsg_name, (fastf_t)tip_length, tip_width);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "tip_width")) {
	if (argc == 2) {
	    fastf_t tip_width = 0.0;
	    if (ged_draw_view_context_data_arrows_tip_get(view_ctx, bsg_name, NULL, &tip_width)) {
		bu_vls_printf(gedp->ged_result_str, "%d", (int)tip_width);
	    } else {
		bu_vls_printf(gedp->ged_result_str, "0");
	    }
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int tip_width;

	    if (bu_sscanf(argv[2], "%d", &tip_width) != 1)
		goto bad;

	    fastf_t tip_length = 0.0;
	    if (ged_draw_view_context_data_arrows_tip_get(view_ctx, bsg_name, &tip_length, NULL))
		ged_draw_view_context_data_arrows_tip_set(view_ctx, bsg_name, tip_length, (fastf_t)tip_width);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

bad:
    return BRLCAD_ERROR;
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
