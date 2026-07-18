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
 * The data-line feature (_tcl_data_lines / _tcl_sdata_lines) is stored as a
 * view-local BObol line set.
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

#include "BObol/BViewController.h"
#include "bu/cmd.h"
#include "bu/malloc.h"
#include "bu/opt.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "bv.h"
#include "../ged_bobol_private.hpp"
#include "../ged_private.h"
#include "./ged_view.h"

struct view_dlines_state {
    struct ged *gedp;
    struct ged_view_context *view_ctx;
    const char *feature_name;
};

static BObolFeatureOwner
_view_dlines_owner(const struct view_dlines_state *vs)
{
    return ged_bobol_view_feature_owner(vs ? vs->view_ctx : nullptr);
}

static BObolViewController *
_view_dlines_controller(const struct view_dlines_state *vs)
{
    return vs ? ged_bobol_view_controller(vs->view_ctx) : nullptr;
}

static BObolFeatureHandle
_view_dlines_handle(const struct view_dlines_state *vs)
{
    BObolViewController *controller = _view_dlines_controller(vs);
    if (!controller || !vs || !vs->feature_name)
	return BObolFeatureHandle();
    const BObolFeatureOwner owner = _view_dlines_owner(vs);
    return controller->features().findOwned(vs->feature_name,
	BOBOL_FEATURE_SCOPE_LOCAL, &owner);
}

static BObolOverlayInfo
_view_dlines_overlay(const struct view_dlines_state *vs)
{
    BObolOverlayInfo overlay;
    overlay.isOverlay = TRUE;
    overlay.ownerToken = vs ? vs->view_ctx : nullptr;
    overlay.role = BObolOverlayRole::Model;
    overlay.overlayClass = BObolOverlayClass::TclOverlay;
    overlay.lifecycle = BObolOverlayLifecycle::PerCommand;
    overlay.order = BObolOverlayOrder::PostTransparent;
    overlay.sourcePath = (vs && vs->feature_name) ? vs->feature_name : "";
    return overlay;
}

static bool
_view_dlines_style(const struct view_dlines_state *vs,
    BObolFeatureStyle *style)
{
    BObolViewController *controller = _view_dlines_controller(vs);
    const BObolFeatureHandle handle = _view_dlines_handle(vs);
    return controller && handle.isValid() && style &&
	controller->features().style(handle, *style);
}

static unsigned char
_view_dlines_color_channel(float value)
{
    if (value <= 0.0f)
	return 0;
    if (value >= 1.0f)
	return 255;
    return static_cast<unsigned char>(value * 255.0f + 0.5f);
}

static bool
_view_dlines_replace(struct view_dlines_state *vs,
    const std::vector<SbVec3f> &points, const BObolFeatureStyle &style)
{
    BObolViewController *controller = _view_dlines_controller(vs);
    if (!controller || !vs || !vs->feature_name)
	return false;

    const BObolFeatureOwner owner = _view_dlines_owner(vs);
    (void)controller->features().removeOwned(vs->feature_name,
	BOBOL_FEATURE_SCOPE_LOCAL, &owner);
    if (points.size() < 2)
	return true;

    std::vector<int32_t> commands(points.size(),
	static_cast<int32_t>(BObolLineCommand::Draw));
    for (size_t i = 0; i + 1 < commands.size(); i += 2) {
	commands[i] = static_cast<int32_t>(BObolLineCommand::Move);
	commands[i + 1] = static_cast<int32_t>(BObolLineCommand::Draw);
    }
    const BObolFeatureHandle handle = controller->features().publishLineSet(
	vs->feature_name, BObolFeatureScope::Local, points, commands, &style,
	&owner);
    return handle.isValid() && controller->features().setOverlayInfo(handle,
	_view_dlines_overlay(vs));
}

static int
_view_dlines_cmd_draw(void *bs, int argc, const char **argv)
{
    struct view_dlines_state *vs = (struct view_dlines_state *)bs;
    struct ged *gedp = vs->gedp;

    if (argc == 1) {
	BObolFeatureStyle style;
	const int visible = _view_dlines_style(vs, &style) &&
	    style.hasVisible && style.visible;
	bu_vls_printf(gedp->ged_result_str, "%d", visible);
	return BRLCAD_OK;
    }

    if (argc == 2) {
	int i;

	if (bu_sscanf(argv[1], "%d", &i) != 1) return BRLCAD_ERROR;

	if (!i) {
	    BObolViewController *controller = _view_dlines_controller(vs);
	    const BObolFeatureOwner owner = _view_dlines_owner(vs);
	    if (controller)
		(void)controller->features().removeOwned(vs->feature_name,
		    BOBOL_FEATURE_SCOPE_LOCAL, &owner);
	}
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
	BObolFeatureStyle style;
	if (_view_dlines_style(vs, &style) && style.hasColor) {
	    bu_vls_printf(gedp->ged_result_str, "%d %d %d",
		_view_dlines_color_channel(style.color[0]),
		_view_dlines_color_channel(style.color[1]),
		_view_dlines_color_channel(style.color[2]));
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

	BObolViewController *controller = _view_dlines_controller(vs);
	const BObolFeatureHandle handle = _view_dlines_handle(vs);
	if (controller && handle.isValid())
	    (void)controller->features().setColor(handle,
		SbColor(static_cast<float>(r) / 255.0f,
		    static_cast<float>(g) / 255.0f,
		    static_cast<float>(b) / 255.0f));

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
	BObolFeatureStyle style;
	if (_view_dlines_style(vs, &style) && style.hasLineWidth)
	    bu_vls_printf(gedp->ged_result_str, "%d", style.lineWidth);
	else
	    bu_vls_printf(gedp->ged_result_str, "0");
	return BRLCAD_OK;
    }

    if (argc == 2) {
	int line_width;

	if (bu_sscanf(argv[1], "%d", &line_width) != 1)
	    return BRLCAD_ERROR;

	BObolViewController *controller = _view_dlines_controller(vs);
	const BObolFeatureHandle handle = _view_dlines_handle(vs);
	if (controller && handle.isValid())
	    (void)controller->features().setLineWidth(handle, line_width);

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
	BObolViewController *controller = _view_dlines_controller(vs);
	const BObolFeatureHandle handle = _view_dlines_handle(vs);
	std::vector<SbVec3f> points;
	if (controller && handle.isValid() &&
	    controller->features().points(handle, points)) {
	    for (const SbVec3f &point : points)
		bu_vls_printf(gedp->ged_result_str, " {%lf %lf %lf} ",
		    static_cast<double>(point[0]),
		    static_cast<double>(point[1]),
		    static_cast<double>(point[2]));
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
	BObolFeatureStyle saved_style;
	if (!_view_dlines_style(vs, &saved_style)) {
	    saved_style.hasVisible = TRUE;
	    saved_style.visible = TRUE;
	    saved_style.hasColor = TRUE;
	    saved_style.color = SbColor(1.0f, 1.0f, 0.0f);
	    saved_style.hasLineWidth = TRUE;
	    saved_style.lineWidth = 0;
	}

	if (ac < 2) {
	    (void)_view_dlines_replace(vs, std::vector<SbVec3f>(),
		saved_style);
	    ged_refresh_cb(gedp);
	    bu_free((char *)av, "av");
	    return BRLCAD_OK;
	}

	std::vector<SbVec3f> points;
	points.reserve(static_cast<size_t>(ac));
	for (i = 0; i < ac; ++i) {
	    double scan[3];

	    if (bu_sscanf(av[i], "%lf %lf %lf", &scan[X], &scan[Y], &scan[Z]) != 3) {
		bu_vls_printf(gedp->ged_result_str, "bad data point - %s\n", av[i]);
		ged_refresh_cb(gedp);
		bu_free((char *)av, "av");
		return BRLCAD_ERROR;
	    }
	    points.emplace_back(static_cast<float>(scan[X]),
		static_cast<float>(scan[Y]), static_cast<float>(scan[Z]));
	}

	if (!_view_dlines_replace(vs, points, saved_style)) {
	    ged_refresh_cb(gedp);
	    bu_free((char *)av, "av");
	    return BRLCAD_ERROR;
	}
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
    if (!ged_view_context_display_endpoint_get(vs.view_ctx) &&
	!ged_view_context_display_endpoint_ensure(vs.view_ctx)) {
	bu_vls_printf(gedp->ged_result_str,
	    "unable to create display endpoint for current view");
	return BRLCAD_ERROR;
    }

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
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
