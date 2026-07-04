/*                           A X E S . C
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
/** @file libtclcad/view/axes.c
 *
 */
/** @} */

#include "common.h"
#include "bu/units.h"
#include "dm.h"
#include "ged.h"
#include "ged/view.h"
#include "rt/view.h"
#include "tclcad.h"

/* Private headers */
#include "ged/draw.h"
#include "../tclcad_private.h"
#include "../view/view.h"

/* Phase T3 (drawing_stack_modernization): data-axes getters and setters now
 * operate through typed GED draw-view data-axis facades.  gv_tcl is no longer
 * read or written by this path; the draw-view feature is the canonical store.
 * BVDAS_DEFAULT_DM_WIDTH is the pixel-width fallback for dm_width in sf calcs. */
#define BVDAS_DEFAULT_DM_WIDTH 512  /* fallback pixel width when no DM is attached */

static fastf_t
_tclcad_data_axes_display_scale(void *view_ctx)
{
    fastf_t dm_width = (fastf_t)BVDAS_DEFAULT_DM_WIDTH;
    struct dm *dmp = (struct dm *)ged_view_context_display_manager_get(view_ctx);
    if (dmp) {
	int width = dm_get_width(dmp);
	if (width > 0)
	    dm_width = (fastf_t)width;
    }

    struct rt_view_info view_info = RT_VIEW_INFO_INIT;
    ged_view_context_info_get(&view_info, view_ctx);
    return view_info.size / dm_width;
}

int
to_axes(struct ged *gedp,
	void *view_ctx,
	struct rt_view_axes_state *gasp,
	int argc,
	const char *argv[],
	const char *usage)
{

    if (BU_STR_EQUAL(argv[2], "draw")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->draw);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int i;

	    if (bu_sscanf(argv[3], "%d", &i) != 1)
		goto bad;

	    if (i)
		gasp->draw = 1;
	    else
		gasp->draw = 0;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "axes_size")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%lf", gasp->axes_size);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    double size; /* must be double for scanf */

	    if (bu_sscanf(argv[3], "%lf", &size) != 1)
		goto bad;

	    gasp->axes_size = size;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "axes_pos")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%lf %lf %lf",
			  V3ARGS(gasp->axes_pos));
	    return BRLCAD_OK;
	}

	if (argc == 6) {
	    double x, y, z; /* must be double for scanf */

	    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
		bu_sscanf(argv[4], "%lf", &y) != 1 ||
		bu_sscanf(argv[5], "%lf", &z) != 1)
		goto bad;

	    VSET(gasp->axes_pos, x, y, z);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "axes_color")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d %d %d",
			  V3ARGS(gasp->axes_color));
	    return BRLCAD_OK;
	}

	if (argc == 6) {
	    int r, g, b;

	    /* set background color */
	    if (bu_sscanf(argv[3], "%d", &r) != 1 ||
		bu_sscanf(argv[4], "%d", &g) != 1 ||
		bu_sscanf(argv[5], "%d", &b) != 1)
		goto bad;

	    /* validate color */
	    if (r < 0 || 255 < r ||
		g < 0 || 255 < g ||
		b < 0 || 255 < b)
		goto bad;

	    VSET(gasp->axes_color, r, g, b);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "label_color")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d %d %d",
			  V3ARGS(gasp->label_color));
	    return BRLCAD_OK;
	}

	if (argc == 6) {
	    int r, g, b;

	    /* set background color */
	    if (bu_sscanf(argv[3], "%d", &r) != 1 ||
		bu_sscanf(argv[4], "%d", &g) != 1 ||
		bu_sscanf(argv[5], "%d", &b) != 1)
		goto bad;

	    /* validate color */
	    if (r < 0 || 255 < r ||
		g < 0 || 255 < g ||
		b < 0 || 255 < b)
		goto bad;

	    VSET(gasp->label_color, r, g, b);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "line_width")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->line_width);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int line_width;

	    if (bu_sscanf(argv[3], "%d", &line_width) != 1)
		goto bad;

	    gasp->line_width = line_width;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "pos_only")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->pos_only);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int i;

	    if (bu_sscanf(argv[3], "%d", &i) != 1)
		goto bad;

	    if (i)
		gasp->pos_only = 1;
	    else
		gasp->pos_only = 0;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "tick_color")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d %d %d",
			  V3ARGS(gasp->tick_color));
	    return BRLCAD_OK;
	}

	if (argc == 6) {
	    int r, g, b;

	    /* set background color */
	    if (bu_sscanf(argv[3], "%d", &r) != 1 ||
		bu_sscanf(argv[4], "%d", &g) != 1 ||
		bu_sscanf(argv[5], "%d", &b) != 1)
		goto bad;

	    /* validate color */
	    if (r < 0 || 255 < r ||
		g < 0 || 255 < g ||
		b < 0 || 255 < b)
		goto bad;

	    VSET(gasp->tick_color, r, g, b);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "tick_enable")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->tick_enabled);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int i;

	    if (bu_sscanf(argv[3], "%d", &i) != 1)
		goto bad;

	    if (i)
		gasp->tick_enabled = 1;
	    else
		gasp->tick_enabled = 0;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "tick_interval")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%f", gasp->tick_interval);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int tick_interval;

	    if (bu_sscanf(argv[3], "%d", &tick_interval) != 1)
		goto bad;

	    gasp->tick_interval = tick_interval;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "tick_length")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->tick_length);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int tick_length;

	    if (bu_sscanf(argv[3], "%d", &tick_length) != 1)
		goto bad;

	    gasp->tick_length = tick_length;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "tick_major_color")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d %d %d",
			  V3ARGS(gasp->tick_major_color));
	    return BRLCAD_OK;
	}

	if (argc == 6) {
	    int r, g, b;

	    /* set background color */
	    if (bu_sscanf(argv[3], "%d", &r) != 1 ||
		bu_sscanf(argv[4], "%d", &g) != 1 ||
		bu_sscanf(argv[5], "%d", &b) != 1)
		goto bad;

	    /* validate color */
	    if (r < 0 || 255 < r ||
		g < 0 || 255 < g ||
		b < 0 || 255 < b)
		goto bad;

	    VSET(gasp->tick_major_color, r, g, b);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "tick_major_length")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->tick_major_length);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int tick_major_length;

	    if (bu_sscanf(argv[3], "%d", &tick_major_length) != 1)
		goto bad;

	    gasp->tick_major_length = tick_major_length;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "ticks_per_major")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->ticks_per_major);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int ticks_per_major;

	    if (bu_sscanf(argv[3], "%d", &ticks_per_major) != 1)
		goto bad;

	    gasp->ticks_per_major = ticks_per_major;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "tick_threshold")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->tick_threshold);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int tick_threshold;

	    if (bu_sscanf(argv[3], "%d", &tick_threshold) != 1)
		goto bad;

	    if (tick_threshold < 1)
		tick_threshold = 1;

	    gasp->tick_threshold = tick_threshold;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "triple_color")) {
	if (argc == 3) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gasp->triple_color);
	    return BRLCAD_OK;
	}

	if (argc == 4) {
	    int i;

	    if (bu_sscanf(argv[3], "%d", &i) != 1)
		goto bad;

	    if (i)
		gasp->triple_color = 1;
	    else
		gasp->triple_color = 0;

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

bad:
    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
    return BRLCAD_ERROR;
}

int
go_data_axes(Tcl_Interp *interp,
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

    ret = to_data_axes_func(interp, gedp, draw_view_ctx, argc, argv);
    to_refresh_suppress_all_end(current_top);
    if (ret & BRLCAD_ERROR)
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);

    return ret;
}

int
to_data_axes(struct ged *gedp,
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

    /* shift the command name to argv[1] before calling to_data_axes_func */
    argv[1] = argv[0];
    ret = to_data_axes_func(current_top->to_interp, gedp, view_ctx, argc-1, argv+1);
    if (ret == BRLCAD_ERROR)
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);

    return ret;
}

int
to_data_axes_func(Tcl_Interp *interp,
		  struct ged *gedp,
		  void *view_ctx,
		  int argc,
		  const char *argv[])
{
    /* T3: BSG object name is the only per-variant state needed here. */
    const char *bsg_name = (argv[0][0] == 's') ? "_tcl_sdata_axes" : "_tcl_data_axes";

    if (BU_STR_EQUAL(argv[1], "draw")) {
	if (argc == 2) {
	    bu_vls_printf(gedp->ged_result_str, "%d",
			  ged_draw_view_context_data_axes_draw_get(view_ctx, bsg_name) ? 1 : 0);
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int i;

	    if (bu_sscanf(argv[2], "%d", &i) != 1)
		goto bad;

	    ged_draw_view_context_data_axes_draw_set(view_ctx, bsg_name, i ? 1 : 0);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "color")) {
	if (argc == 2) {
	    struct ged_draw_view_feature_style style = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	    if (ged_draw_view_context_data_axes_style_get(view_ctx, bsg_name, &style) && style.color_valid) {
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

	    ged_draw_view_context_data_axes_color_set(view_ctx, bsg_name, r, g, b);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "line_width")) {
	if (argc == 2) {
	    struct ged_draw_view_feature_style style = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	    if (ged_draw_view_context_data_axes_style_get(view_ctx, bsg_name, &style))
		bu_vls_printf(gedp->ged_result_str, "%d", style.line_width);
	    else
		bu_vls_printf(gedp->ged_result_str, "0");
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int line_width;

	    if (bu_sscanf(argv[2], "%d", &line_width) != 1)
		goto bad;

	    ged_draw_view_context_data_axes_line_width_set(view_ctx, bsg_name, line_width);

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "size")) {
	if (argc == 2) {
	    fastf_t size = 0.0;
	    fastf_t sf = _tclcad_data_axes_display_scale(view_ctx);
	    ged_draw_view_context_data_axes_size_get(view_ctx, bsg_name, sf,
		    &size);
	    bu_vls_printf(gedp->ged_result_str, "%lf", size);
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    double size; /* must be double for scanf */

	    if (bu_sscanf(argv[2], "%lf", &size) != 1)
		goto bad;

	    /* Extract current centers and rebuild with new halfAxesSize. */
	    struct ged_draw_view_feature_style saved_style = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	    ged_draw_view_context_data_axes_style_get(view_ctx, bsg_name,
		    &saved_style);

	    point_t *cpts = NULL;
	    size_t ncpts = 0;
	    (void)ged_draw_view_context_data_axes_centers_copy(view_ctx,
		    bsg_name, &cpts, &ncpts);

	    fastf_t sf = _tclcad_data_axes_display_scale(view_ctx);
	    fastf_t half = (fastf_t)size * 0.5f * sf;

	    if (cpts && ncpts)
		(void)ged_draw_view_context_data_axes_centers_replace(view_ctx, bsg_name, cpts, ncpts,
			half, &saved_style);
	    if (cpts)
		bu_free(cpts, "GED draw view axes centers copy");

	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "points")) {
	register int i;

	if (argc == 2) {
	    point_t *cpts = NULL;
	    size_t ncpts = 0;
	    if (ged_draw_view_context_data_axes_centers_copy(view_ctx, bsg_name, &cpts, &ncpts)) {
		for (size_t j = 0; j < ncpts; ++j)
		    bu_vls_printf(gedp->ged_result_str, " {%lf %lf %lf} ", V3ARGS(cpts[j]));
		if (cpts)
		    bu_free(cpts, "GED draw view axes centers copy");
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

	    /* Save style and size from existing object before replacing it. */
	    struct ged_draw_view_feature_style saved_style = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	    ged_draw_view_context_data_axes_style_get(view_ctx, bsg_name,
		    &saved_style);

	    /* Recover halfAxesSize from existing object (use default 1.0 if none). */
	    fastf_t half = 1.0;
	    ged_draw_view_context_data_axes_half_size_get(view_ctx, bsg_name,
		    &half);

	    /* Clear out: remove old GED draw-view feature. */
	    if (ac < 1) {
		ged_draw_view_context_data_axes_centers_replace(view_ctx,
			bsg_name, NULL, 0, half, &saved_style);
		to_refresh_view(view_ctx);
		Tcl_Free((char *)av);
		return BRLCAD_OK;
	    }

	    /* Parse new center points into temporary local array. */
	    point_t *pts = (point_t *)bu_calloc(ac, sizeof(point_t), "axes points");
	    for (i = 0; i < ac; ++i) {
		double scan[3];

		if (bu_sscanf(av[i], "%lf %lf %lf", &scan[X], &scan[Y], &scan[Z]) != 3) {
		    bu_vls_printf(gedp->ged_result_str, "bad data point - %s\n", av[i]);
		    bu_free(pts, "axes points");
		    Tcl_Free((char *)av);
		    return BRLCAD_ERROR;
		}
		VMOVE(pts[i], scan);
	    }

	    /* Rebuild draw-view data axes from new centers, preserving style. */
	    (void)ged_draw_view_context_data_axes_centers_replace(view_ctx, bsg_name, pts, (size_t)ac,
		    half, &saved_style);
	    bu_free(pts, "axes points");
	    Tcl_Free((char *)av);
	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}
    }

bad:
    return BRLCAD_ERROR;
}

int
to_model_axes(struct ged *gedp,
	      int argc,
	      const char *argv[],
	      ged_func_ptr UNUSED(func),
	      const char *usage,
	      int UNUSED(maxargs))
{
    void *view_ctx;

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

    struct rt_view_axes_state axes;
    if (!ged_view_context_model_axes_state_get(&axes, view_ctx))
	return BRLCAD_ERROR;
    int ret = to_axes(gedp, view_ctx, &axes, argc, argv, usage);
    if (ret == BRLCAD_OK)
	ged_view_context_model_axes_state_set(view_ctx, &axes);
    return ret;
}

int
go_view_axes(struct ged *gedp,
	     void *draw_view_ctx,
	     int argc,
	     const char *argv[],
	     const char *usage)
{
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

    struct rt_view_axes_state axes;
    if (!ged_view_context_view_axes_state_get(&axes, draw_view_ctx))
	return BRLCAD_ERROR;
    int ret = to_axes(gedp, draw_view_ctx, &axes, argc, argv, usage);
    if (ret == BRLCAD_OK)
	ged_view_context_view_axes_state_set(draw_view_ctx, &axes);
    return ret;
}


int
to_view_axes(struct ged *gedp,
	     int argc,
	     const char *argv[],
	     ged_func_ptr UNUSED(func),
	     const char *usage,
	     int UNUSED(maxargs))
{
    void *view_ctx;

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

    struct rt_view_axes_state axes;
    if (!ged_view_context_view_axes_state_get(&axes, view_ctx))
	return BRLCAD_ERROR;
    int ret = to_axes(gedp, view_ctx, &axes, argc, argv, usage);
    if (ret == BRLCAD_OK)
	ged_view_context_view_axes_state_set(view_ctx, &axes);
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
