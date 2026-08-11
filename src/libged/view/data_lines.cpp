/*                    D A T A _ L I N E S . C P P
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
/** @file libged/view/data_lines.cpp
 *
 * Logic for drawing arbitrary lines not associated with geometry.
 *
 * Command state is owned by the GED view host and published as a view-local
 * retained feature batch.  Renderer realization is never queried to answer a
 * command getter.
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
#include <vector>

#include "bu/cmd.h"
#include "bu/malloc.h"
#include "bu/opt.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "bv.h"
#include "ged/view_feature_batch.h"
#include "../ged_private.h"
#include "../ged_view_data_line_private.h"
#include "./ged_view.h"

struct view_dlines_state {
    struct ged *gedp;
    struct ged_view_context *view_ctx;
    const char *feature_name;
    int staged;
};

static bool
_view_dlines_publish(const struct view_dlines_state *vs,
    const struct ged_view_data_line_state *state)
{
    if (!vs || !vs->view_ctx || !vs->feature_name || !state)
	return false;
    struct ged_view_feature_batch_desc desc = GED_VIEW_FEATURE_BATCH_DESC_INIT;
    desc.owner_id = "ged-data-lines";
    desc.owner_role = "user-overlay";
    desc.overlay_class = GED_VIEW_FEATURE_OVERLAY_CLASS_USER_ANNOTATION;
    desc.lifecycle = GED_VIEW_FEATURE_LIFECYCLE_PERSISTENT;
    desc.overlay_order = GED_VIEW_FEATURE_OVERLAY_ORDER_MODEL;
    desc.local = 1;
    struct ged_view_feature_batch *batch =
	ged_view_feature_batch_begin(vs->view_ctx, &desc);
    if (!batch)
	return false;
    std::vector<int> commands(state->point_count, GED_DRAW_VIEW_LINE_DRAW);
    for (size_t i = 0; i + 1 < commands.size(); i += 2) {
	commands[i] = GED_DRAW_VIEW_LINE_MOVE;
	commands[i + 1] = GED_DRAW_VIEW_LINE_DRAW;
    }
    struct ged_view_feature_style style = GED_VIEW_FEATURE_STYLE_INIT;
    style.visible = state->draw ? 1 : 0;
    style.selectable = 1;
    style.color_valid = 1;
    style.color[0] = (unsigned char)state->color[0];
    style.color[1] = (unsigned char)state->color[1];
    style.color[2] = (unsigned char)state->color[2];
    style.line_width = state->line_width;
    const point_t *points = state->draw && state->point_count >= 2 ?
	(const point_t *)state->points : NULL;
    const size_t point_count = points ? state->point_count : 0;
    if (!ged_view_feature_batch_line_set_replace(batch, vs->feature_name,
	    points, point_count ? commands.data() : NULL, point_count, &style)) {
	ged_view_feature_batch_abort(batch);
	return false;
    }
    return ged_view_feature_batch_commit(batch) != 0;
}

static bool
_view_dlines_update(const struct view_dlines_state *vs,
    const struct ged_view_data_line_state *state)
{
    return vs && state && ged_view_data_line_state_replace(vs->view_ctx,
	vs->staged, state) && _view_dlines_publish(vs, state);
}

extern "C" GED_EXPORT int
ged_view_data_line_state_publish(struct ged_view_context *view_ctx, int staged)
{
    if (!view_ctx)
	return 0;
    struct ged_view_data_line_state state = GED_VIEW_DATA_LINE_STATE_INIT;
    if (!ged_view_data_line_state_get(view_ctx, staged, &state))
	return 0;
    struct view_dlines_state vs = {};
    vs.view_ctx = view_ctx;
    vs.feature_name = staged ? "_tcl_sdata_lines" : "_tcl_data_lines";
    vs.staged = staged ? 1 : 0;
    const int ret = _view_dlines_publish(&vs, &state) ? 1 : 0;
    ged_view_data_line_state_clear(&state);
    return ret;
}

static int
_view_dlines_cmd_draw(void *bs, int argc, const char **argv)
{
    struct view_dlines_state *vs = (struct view_dlines_state *)bs;
    struct ged *gedp = vs->gedp;

    if (argc == 1) {
	struct ged_view_data_line_state state = GED_VIEW_DATA_LINE_STATE_INIT;
	if (!ged_view_data_line_state_get(vs->view_ctx, vs->staged, &state))
	    return BRLCAD_ERROR;
	bu_vls_printf(gedp->ged_result_str, "%d", state.draw);
	ged_view_data_line_state_clear(&state);
	return BRLCAD_OK;
    }

    if (argc == 2) {
	int i;

	if (bu_sscanf(argv[1], "%d", &i) != 1) return BRLCAD_ERROR;

	struct ged_view_data_line_state state = GED_VIEW_DATA_LINE_STATE_INIT;
	if (!ged_view_data_line_state_get(vs->view_ctx, vs->staged, &state))
	    return BRLCAD_ERROR;
	state.draw = i ? 1 : 0;
	const bool updated = _view_dlines_update(vs, &state);
	ged_view_data_line_state_clear(&state);
	if (!updated)
	    return BRLCAD_ERROR;

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
    struct ged_view_context *view_ctx = vs->view_ctx;
    struct bv *view = bv_context_view(ged_view_context_bv(view_ctx));
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

    if (argc == 1) {
	struct ged_view_data_line_state state = GED_VIEW_DATA_LINE_STATE_INIT;
	if (!ged_view_data_line_state_get(vs->view_ctx, vs->staged, &state))
	    return BRLCAD_ERROR;
	bu_vls_printf(gedp->ged_result_str, "%d %d %d", V3ARGS(state.color));
	ged_view_data_line_state_clear(&state);
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

	struct ged_view_data_line_state state = GED_VIEW_DATA_LINE_STATE_INIT;
	if (!ged_view_data_line_state_get(vs->view_ctx, vs->staged, &state))
	    return BRLCAD_ERROR;
	VSET(state.color, r, g, b);
	const bool updated = _view_dlines_update(vs, &state);
	ged_view_data_line_state_clear(&state);
	if (!updated)
	    return BRLCAD_ERROR;

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

    if (argc == 1) {
	struct ged_view_data_line_state state = GED_VIEW_DATA_LINE_STATE_INIT;
	if (!ged_view_data_line_state_get(vs->view_ctx, vs->staged, &state))
	    return BRLCAD_ERROR;
	bu_vls_printf(gedp->ged_result_str, "%d", state.line_width);
	ged_view_data_line_state_clear(&state);
	return BRLCAD_OK;
    }

    if (argc == 2) {
	int line_width;

	if (bu_sscanf(argv[1], "%d", &line_width) != 1)
	    return BRLCAD_ERROR;

	struct ged_view_data_line_state state = GED_VIEW_DATA_LINE_STATE_INIT;
	if (!ged_view_data_line_state_get(vs->view_ctx, vs->staged, &state))
	    return BRLCAD_ERROR;
	state.line_width = line_width;
	const bool updated = _view_dlines_update(vs, &state);
	ged_view_data_line_state_clear(&state);
	if (!updated)
	    return BRLCAD_ERROR;

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
    int i;

    if (argc == 1) {
	struct ged_view_data_line_state state = GED_VIEW_DATA_LINE_STATE_INIT;
	if (!ged_view_data_line_state_get(vs->view_ctx, vs->staged, &state))
	    return BRLCAD_ERROR;
	for (size_t j = 0; j < state.point_count; j++)
	    bu_vls_printf(gedp->ged_result_str, " {%lf %lf %lf} ",
		V3ARGS(state.points[j]));
	ged_view_data_line_state_clear(&state);
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

	struct ged_view_data_line_state state = GED_VIEW_DATA_LINE_STATE_INIT;
	if (!ged_view_data_line_state_get(vs->view_ctx, vs->staged, &state)) {
	    bu_free((char *)av, "av");
	    return BRLCAD_ERROR;
	}

	if (ac < 2) {
	    if (state.points)
		bu_free(state.points, "GED view data-line state points");
	    state.points = NULL;
	    state.point_count = 0;
	    const bool updated = _view_dlines_update(vs, &state);
	    ged_view_data_line_state_clear(&state);
	    ged_refresh_cb(gedp);
	    bu_free((char *)av, "av");
	    return updated ? BRLCAD_OK : BRLCAD_ERROR;
	}

	point_t *points = (point_t *)bu_calloc((size_t)ac, sizeof(point_t),
	    "GED data-lines points");
	for (i = 0; i < ac; ++i) {
	    double scan[3];

	    if (bu_sscanf(av[i], "%lf %lf %lf", &scan[X], &scan[Y], &scan[Z]) != 3) {
		bu_vls_printf(gedp->ged_result_str, "bad data point - %s\n", av[i]);
		bu_free(points, "GED data-lines points");
		ged_view_data_line_state_clear(&state);
		ged_refresh_cb(gedp);
		bu_free((char *)av, "av");
		return BRLCAD_ERROR;
	    }
	    VSET(points[i], scan[X], scan[Y], scan[Z]);
	}

	if (state.points)
	    bu_free(state.points, "GED view data-line state points");
	state.points = points;
	state.point_count = (size_t)ac;
	state.draw = 1;
	if (!_view_dlines_update(vs, &state)) {
	    ged_view_data_line_state_clear(&state);
	    ged_refresh_cb(gedp);
	    bu_free((char *)av, "av");
	    return BRLCAD_ERROR;
	}
	ged_view_data_line_state_clear(&state);
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
	vs.staged = 1;
    } else {
	vs.feature_name = "_tcl_data_lines";
	vs.staged = 0;
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
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
