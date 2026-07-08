/*                      D A T A _ L I N E S . C
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
/** @file libged/data_lines.c
 *
 * Logic for drawing arbitrary lines not associated with geometry.
 *
 * The data-line feature (_tcl_data_lines / _tcl_sdata_lines) is managed by
 * GED draw-view data-line facades.  The command keeps the active view opaque
 * and does not manipulate feature-store internals directly.
 *
 * Usage example (Archer / QGED):
 *
 * Archer> sdata_lines points {{0 -1000 0} {0 1000 0} {100 -1000 0} {100 1000 0} {-1000 10 0} {1000 10 0}}
 * Archer> sdata_lines draw 1
 * Archer> sdata_lines line_width 100
 * Archer> sdata_lines color 255 0 0
 *
 * Note that the active GED view must be set to the correct display manager
 * before calling this command to put the output in the correct display
 * manager.
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "bu/cmd.h"
#include "bu/malloc.h"
#include "bu/opt.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "bv.h"
#include "../ged_private.h"
#include "./ged_view.h"

struct view_dlines_state {
    struct ged *gedp;
    void *view_ctx;
    const char *feature_name;
};

static int
_view_dlines_cmd_draw(void *bs, int argc, const char **argv)
{
    struct view_dlines_state *vs = (struct view_dlines_state *)bs;
    struct ged *gedp = vs->gedp;
    void *view_ctx = vs->view_ctx;

    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "%d",
		ged_draw_view_context_data_lines_draw_get(view_ctx,
		    vs->feature_name));
	return BRLCAD_OK;
    }

    if (argc == 2) {
	int i;

	if (bu_sscanf(argv[1], "%d", &i) != 1) return BRLCAD_ERROR;

	(void)ged_draw_view_context_data_lines_draw_set(view_ctx,
		vs->feature_name, i);
	/* draw=1 is a no-op here; use "points" to create/re-enable. */

	ged_refresh_cb(gedp);

	return BRLCAD_OK;
    }

    return BRLCAD_ERROR;
}

static int
_view_dlines_cmd_snap(void *bs, int argc, const char **argv)
{
    struct view_dlines_state *vs = (struct view_dlines_state *)bs;
    struct ged *gedp = vs->gedp;
    void *view_ctx = vs->view_ctx;
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "%d", bv_snap_lines_get(view));
	return BRLCAD_OK;
    }

    if (argc == 2) {
	int i;

	if (bu_sscanf(argv[1], "%d", &i) != 1) return BRLCAD_ERROR;

	bv_snap_lines_set(view, i);

	return BRLCAD_OK;
    }

    return BRLCAD_ERROR;
}

static int
_view_dlines_cmd_color(void *bs, int argc, const char **argv)
{
    struct view_dlines_state *vs = (struct view_dlines_state *)bs;
    struct ged *gedp = vs->gedp;
    void *view_ctx = vs->view_ctx;

    if (argc == 1) {
	struct ged_draw_view_line_style style = {{0, 0, 0}, 0};
	if (ged_draw_view_context_data_lines_style_get(view_ctx,
		vs->feature_name, &style)) {
	    bu_vls_printf(gedp->ged_result_str, "%d %d %d",
			  style.color[0], style.color[1], style.color[2]);
	} else
	    bu_vls_printf(gedp->ged_result_str, "0 0 0");
	return BRLCAD_OK;
    }

    if (argc == 4) {
	int r, g, b;

	/* set background color */
	if (bu_sscanf(argv[1], "%d", &r) != 1 ||
		bu_sscanf(argv[2], "%d", &g) != 1 ||
		bu_sscanf(argv[3], "%d", &b) != 1)
	    return BRLCAD_ERROR;

	/* validate color */
	if (r < 0 || 255 < r ||
		g < 0 || 255 < g ||
		b < 0 || 255 < b)
	    return BRLCAD_ERROR;

	(void)ged_draw_view_context_data_lines_color_set(view_ctx,
		vs->feature_name, r, g, b);

	ged_refresh_cb(gedp);

	return BRLCAD_OK;
    }

    return BRLCAD_ERROR;
}

static int
_view_dlines_cmd_line_width(void *bs, int argc, const char **argv)
{
    struct view_dlines_state *vs = (struct view_dlines_state *)bs;
    struct ged *gedp = vs->gedp;
    void *view_ctx = vs->view_ctx;

    if (argc == 1) {
	struct ged_draw_view_line_style style = {{0, 0, 0}, 0};
	if (ged_draw_view_context_data_lines_style_get(view_ctx,
		vs->feature_name, &style))
	    bu_vls_printf(gedp->ged_result_str, "%d", style.line_width);
	else
	    bu_vls_printf(gedp->ged_result_str, "0");
	return BRLCAD_OK;
    }

    if (argc == 2) {
	int line_width;

	if (bu_sscanf(argv[1], "%d", &line_width) != 1)
	    return BRLCAD_ERROR;

	(void)ged_draw_view_context_data_lines_line_width_set(view_ctx,
		vs->feature_name, line_width);

	ged_refresh_cb(gedp);

	return BRLCAD_OK;
    }

    return BRLCAD_ERROR;
}

static int
_view_dlines_cmd_points(void *bs, int argc, const char **argv)
{
    struct view_dlines_state *vs = (struct view_dlines_state *)bs;
    struct ged *gedp = vs->gedp;
    void *view_ctx = vs->view_ctx;
    int i;

    if (argc == 1) {
	point_t *points = NULL;
	size_t point_count = 0;
	if (ged_draw_view_context_data_lines_points_copy(view_ctx,
		vs->feature_name, &points, &point_count)) {
	    for (size_t j = 0; j < point_count; j++)
		bu_vls_printf(gedp->ged_result_str, " {%lf %lf %lf} ", V3ARGS(points[j]));
	    if (points)
		bu_free(points, "GED draw view line points copy");
	}
	return BRLCAD_OK;
    }

    if (argc == 2) {
	int ac;
	const char **av;

	if (bu_argv_from_tcl_list(argv[1], &ac, &av)) {
	    bu_vls_printf(gedp->ged_result_str, "failed to parse list");
	    return BRLCAD_ERROR;
	}

	if (ac % 2) {
	    bu_vls_printf(gedp->ged_result_str, "%s: must be an even number of points", argv[0]);
	    bu_free((char *)av, "av");
	    return BRLCAD_ERROR;
	}

	/* Preserve existing style when replacing the points payload. */
	struct ged_draw_view_line_style saved_style = {{255, 255, 0}, 0};
	(void)ged_draw_view_context_data_lines_style_get(view_ctx,
		vs->feature_name, &saved_style);

	if (ac < 2) {
	    (void)ged_draw_view_context_data_lines_points_replace(view_ctx,
		    vs->feature_name, NULL, 0, NULL);
	    ged_refresh_cb(gedp);
	    bu_free((char *)av, "av");
	    return BRLCAD_OK;
	}

	point_t *pts = (point_t *)bu_calloc(ac, sizeof(point_t), "data points");
	for (i = 0; i < ac; ++i) {
	    double scan[3];

	    if (bu_sscanf(av[i], "%lf %lf %lf", &scan[X], &scan[Y], &scan[Z]) != 3) {
		bu_vls_printf(gedp->ged_result_str, "bad data point - %s\n", av[i]);
		bu_free((void *)pts, "data points");
		ged_refresh_cb(gedp);
		bu_free((char *)av, "av");
		return BRLCAD_ERROR;
	    }
	    VMOVE(pts[i], scan);
	}

	if (!ged_draw_view_context_data_lines_points_replace(view_ctx, vs->feature_name,
		(const point_t *)pts, (size_t)ac, &saved_style)) {
	    bu_free((void *)pts, "data points");
	    ged_refresh_cb(gedp);
	    bu_free((char *)av, "av");
	    return BRLCAD_ERROR;
	}
	bu_free((void *)pts, "data points");

	ged_refresh_cb(gedp);
	bu_free((char *)av, "av");
	return BRLCAD_OK;
    }

    return BRLCAD_ERROR;
}

const struct bu_cmdtab _view_dline_cmds[] = {
    { "draw",       _view_dlines_cmd_draw},
    { "color",      _view_dlines_cmd_color},
    { "line_width", _view_dlines_cmd_line_width},
    { "points",     _view_dlines_cmd_points},
    { "snap",       _view_dlines_cmd_snap},
    { (char *)NULL,      NULL}
};

int
ged_view_data_lines(struct ged *gedp, int argc, const char *argv[])
{
    struct view_dlines_state vs;
    const char *usage = "[s]data_lines [subcommand]";

    vs.gedp = gedp;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (6 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    GED_CHECK_VIEW(gedp, BRLCAD_ERROR);
    vs.view_ctx = ged_view_active_ctx(gedp);

    if (argv[0][0] == 's') {
	vs.feature_name = "_tcl_sdata_lines";
    } else {
	vs.feature_name = "_tcl_data_lines";
    }

    argc--;argv++;

    if (bu_cmd_valid(_view_dline_cmds, argv[0]) != BRLCAD_OK) {
	bu_vls_printf(gedp->ged_result_str, "invalid subcommand: %s", argv[1]);
	return BRLCAD_ERROR;
    }

    int ret;
    if (bu_cmd(_view_dline_cmds, argc, argv, 0, (void *)&vs, &ret) == BRLCAD_OK) {
	return BRLCAD_OK;
    }

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
