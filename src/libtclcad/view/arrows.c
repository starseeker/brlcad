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
#include "ged/view_feature_batch.h"
#include "tclcad.h"

/* Private headers */
#include "ged/draw.h"
#include "../tclcad_private.h"
#include "../view/view.h"

int
tclcad_data_arrows_publish(struct ged_view_context *view_ctx,
	const char *name, const struct tclcad_data_arrow_state *state)
{
    if (!view_ctx || !name || !state)
	return 0;

    struct ged_view_feature_batch_desc desc = GED_VIEW_FEATURE_BATCH_DESC_INIT;
    desc.owner_id = "tclcad-arrows";
    desc.owner_role = "tcl-overlay";
    desc.overlay_class = GED_VIEW_FEATURE_OVERLAY_CLASS_TCL_OVERLAY;
    desc.local = 1;
    struct ged_view_feature_batch *batch =
	ged_view_feature_batch_begin(view_ctx, &desc);
    if (!batch)
	return 0;

    struct ged_view_feature_style style = GED_VIEW_FEATURE_STYLE_INIT;
    style.visible = state->gdas_draw ? 1 : 0;
    style.selectable = 1;
    style.color_valid = 1;
    style.color[0] = (unsigned char)state->gdas_color[0];
    style.color[1] = (unsigned char)state->gdas_color[1];
    style.color[2] = (unsigned char)state->gdas_color[2];
    style.line_width = state->gdas_line_width;
    style.arrow = 1;
    style.arrow_tip_length = (fastf_t)state->gdas_tip_length;
    style.arrow_tip_width = (fastf_t)state->gdas_tip_width;
    const point_t *points = state->gdas_draw && state->gdas_num_points > 0 ?
	(const point_t *)state->gdas_points : NULL;
    const size_t count = points ? (size_t)state->gdas_num_points : 0;
    if (!ged_view_feature_batch_arrow_replace(batch, name, points, count,
	    &style)) {
	ged_view_feature_batch_abort(batch);
	return 0;
    }
    return ged_view_feature_batch_commit(batch);
}

int
go_data_arrows(Tcl_Interp *interp,
	       struct ged *gedp,
	       struct ged_view_context *draw_view_ctx,
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
    struct ged_view_context *view_ctx;
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
		    struct ged_view_context *view_ctx,
		    int argc,
		    const char *argv[])
{
    const int staged = argv[0][0] == 's';
    const char *feature_name = staged ? "_tcl_sdata_arrows" : "_tcl_data_arrows";
    tclcad_view_state *view_state =
	tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!view_state)
	return BRLCAD_ERROR;
    struct tclcad_data_arrow_state *state = staged ?
	&view_state->gv_sdata_arrows : &view_state->gv_data_arrows;

    if (BU_STR_EQUAL(argv[1], "draw")) {
	if (argc == 2) {
	    bu_vls_printf(gedp->ged_result_str, "%d", state->gdas_draw);
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int i;

	    if (bu_sscanf(argv[2], "%d", &i) != 1)
		goto bad;

	    state->gdas_draw = i ? 1 : 0;
	    (void)tclcad_data_arrows_publish(view_ctx, feature_name, state);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "color")) {
	if (argc == 2) {
	    bu_vls_printf(gedp->ged_result_str, "%d %d %d",
		state->gdas_color[0], state->gdas_color[1],
		state->gdas_color[2]);
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

	    VSET(state->gdas_color, r, g, b);
	    (void)tclcad_data_arrows_publish(view_ctx, feature_name, state);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "line_width")) {
	if (argc == 2) {
	    bu_vls_printf(gedp->ged_result_str, "%d", state->gdas_line_width);
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int line_width;

	    if (bu_sscanf(argv[2], "%d", &line_width) != 1)
		goto bad;

	    state->gdas_line_width = line_width;
	    (void)tclcad_data_arrows_publish(view_ctx, feature_name, state);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "points")) {
	register int i;

	if (argc == 2) {
	    for (int j = 0; j < state->gdas_num_points; j++)
		bu_vls_printf(gedp->ged_result_str, " {%lf %lf %lf} ",
		    V3ARGS(state->gdas_points[j]));
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

	    if (state->gdas_points) {
		bu_free(state->gdas_points, "TclCAD arrow points");
		state->gdas_points = NULL;
	    }
	    state->gdas_num_points = 0;

	    if (ac < 2) {
		(void)tclcad_data_arrows_publish(view_ctx, feature_name, state);
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

	    state->gdas_points = pts;
	    state->gdas_num_points = ac;
	    (void)tclcad_data_arrows_publish(view_ctx, feature_name, state);
	    Tcl_Free((char *)av);
	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}
    }

    if (BU_STR_EQUAL(argv[1], "tip_length")) {
	if (argc == 2) {
	    bu_vls_printf(gedp->ged_result_str, "%d", state->gdas_tip_length);
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int tip_length;

	    if (bu_sscanf(argv[2], "%d", &tip_length) != 1)
		goto bad;

	    state->gdas_tip_length = tip_length;
	    (void)tclcad_data_arrows_publish(view_ctx, feature_name, state);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "tip_width")) {
	if (argc == 2) {
	    bu_vls_printf(gedp->ged_result_str, "%d", state->gdas_tip_width);
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int tip_width;

	    if (bu_sscanf(argv[2], "%d", &tip_width) != 1)
		goto bad;

	    state->gdas_tip_width = tip_width;
	    (void)tclcad_data_arrows_publish(view_ctx, feature_name, state);

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
