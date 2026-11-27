/*                       A U T O V I E W . C P P
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
/** @file libged/autoview.cpp
 *
 * The autoview command.
 *
 */

#include "common.h"

#include <cstdlib>
#include <algorithm>
#include <string>
#include <vector>

#include "BObol/BDatabaseSource.h"
#include "BObol/BMeasureAction.h"
#include "BObol/BSceneController.h"
#include "BObol/BViewController.h"
#include <Inventor/SoViewport.h>
#include "bu/opt.h"
#include "bv.h"
#include "ged/draw.h"
#include "ged/scene.h"
#include "../ged_bobol_private.hpp"
#include "../ged_draw_private.h"
#include "../ged_private.h"

/* Return 1 (and set *v) if the entire string parses as a number. */
static int
_autoview_arg_is_num(const char *s, double *v)
{
    char *endptr = NULL;
    double d;

    if (!s || *s == '\0')
	return 0;

    d = strtod(s, &endptr);
    if (endptr == s || *endptr != '\0')
	return 0;

    if (v)
	*v = d;
    return 1;
}

static int
_autoview_bounds_finite(const point_t min, const point_t max)
{
    return min && max && isfinite(min[X]) && isfinite(min[Y]) &&
	isfinite(min[Z]) && isfinite(max[X]) && isfinite(max[Y]) &&
	isfinite(max[Z]) && min[X] <= max[X] && min[Y] <= max[Y] &&
	min[Z] <= max[Z];
}

static int
_autoview_named_database_bounds(struct ged *gedp, int argc,
	const char *argv[], point_t min, point_t max)
{
    if (!gedp || !gedp->dbip || argc <= 0 || !argv)
	return 0;

    if (rt_obj_bounds(gedp->ged_result_str, gedp->dbip, argc, argv, 0,
	    min, max) == BRLCAD_OK && _autoview_bounds_finite(min, max))
	return 1;

    /* A combination containing an unbounded halfspace can make its raw
     * aggregate bounds non-finite.  Retain every finite requested object and
     * fall back to the finite, non-subtractive members of an affected
     * combination.  This is the same framing policy used by retained scene
     * sources and avoids poisoning the camera with NaN values. */
    VSETALL(min, INFINITY);
    VSETALL(max, -INFINITY);
    int have_bounds = 0;
    for (int i = 0; i < argc; i++) {
	point_t path_min;
	point_t path_max;
	struct bu_vls messages = BU_VLS_INIT_ZERO;
	const char *path = argv[i];
	int bounded = rt_obj_bounds(&messages, gedp->dbip, 1, &path, 0,
	    path_min, path_max) == BRLCAD_OK &&
	    _autoview_bounds_finite(path_min, path_max);
	bu_vls_free(&messages);
	if (!bounded)
	    bounded = ged_database_path_member_autoview_bounds(gedp, path,
		&path_min, &path_max);
	if (!bounded || !_autoview_bounds_finite(path_min, path_max))
	    continue;
	VMIN(min, path_min);
	VMAX(max, path_max);
	have_bounds = 1;
    }
    return have_bounds;
}

static int
_autoview_semantic_draw_bounds(struct ged *gedp,
	struct ged_view_context *view_ctx, vect_t *min, vect_t *max)
{
    if (!gedp || !gedp->dbip || !min || !max)
	return 0;

    struct bu_vls listed = BU_VLS_INIT_ZERO;
    size_t path_count = ged_scene_paths_append(gedp, view_ctx,
	GED_SCENE_DRAW_DEFAULT, GED_SCENE_PATHS_DRAW_INTENTS, &listed);
    if (!path_count && view_ctx)
	path_count = ged_scene_paths_append(gedp, NULL,
	    GED_SCENE_DRAW_DEFAULT, GED_SCENE_PATHS_DRAW_INTENTS, &listed);
    if (!path_count || !bu_vls_strlen(&listed)) {
	bu_vls_free(&listed);
	return 0;
    }

    std::vector<std::string> paths;
    const char *begin = bu_vls_cstr(&listed);
    const char *cursor = begin;
    while (true) {
	if (*cursor != '\n' && *cursor != '\r' && *cursor != '\0') {
	    cursor++;
	    continue;
	}
	if (cursor > begin)
	    paths.emplace_back(begin, static_cast<size_t>(cursor - begin));
	if (*cursor == '\0')
	    break;
	begin = ++cursor;
    }

    std::vector<const char *> path_args;
    path_args.reserve(paths.size());
    for (const std::string &path : paths)
	path_args.push_back(path.c_str());

    struct bu_vls messages = BU_VLS_INIT_ZERO;
    const int have_bounds = !path_args.empty() &&
	rt_obj_bounds(&messages, gedp->dbip,
	    static_cast<int>(path_args.size()), path_args.data(), 0, *min,
	    *max) == BRLCAD_OK && _autoview_bounds_finite(*min, *max);
    bu_vls_free(&messages);
    bu_vls_free(&listed);
    return have_bounds;
}

static int
_autoview_bobol_view_bounds(struct ged_view_context *view_ctx,
	point_t min, point_t max)
{
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!controller || !controller->getViewport() ||
	!controller->getViewport()->getRoot())
	return 0;

    SoBRLMeasureAction measure;
    measure.setGeometryPolicy(SoBRLMeasureAction::DISPLAY_LEVEL);
    measure.apply(controller->getViewport()->getRoot());
    const SbBox3f &bounds = measure.getBounds();
    if (bounds.isEmpty())
	return 0;

    const SbVec3f bmin = bounds.getMin();
    const SbVec3f bmax = bounds.getMax();
    VSET(min, bmin[X], bmin[Y], bmin[Z]);
    VSET(max, bmax[X], bmax[Y], bmax[Z]);
    return 1;
}

struct autoview_request {
    struct ged *gedp;
    struct ged_view_context *view_ctx;
    fastf_t factor;
    bool named_objects;
    int all_view_objs;
    point_t min;
    point_t max;
};

static int
_autoview_on_owner(void *data)
{
    autoview_request &request = *static_cast<autoview_request *>(data);
    if (!request.named_objects) {
	if (request.all_view_objs) {
	    if (!_autoview_bobol_view_bounds(request.view_ctx, request.min,
		    request.max))
		return BRLCAD_OK;
	} else {
	    /* A progressive provider owns the fit until it publishes final
	     * bounds or the user moves the camera.  Let it use the bounds it
	     * already has before considering expensive database fallback. */
	    if (ged_draw_obol_progressive_autoview_follow(request.gedp,
		    request.view_ctx, request.factor))
		return BRLCAD_OK;

	    /* A null presentation still needs finite database bounds to
	     * initialize the camera for callers such as rtwizard. */
	    int empty = 1;
	    if ((!ged_draw_obol_scene_database_autoview_bounds(request.gedp,
			&request.min, &request.max, &empty, 1) || empty) &&
		!_autoview_semantic_draw_bounds(request.gedp, request.view_ctx,
		    &request.min, &request.max))
		return BRLCAD_OK;
	}
    }

    bv_autoview_bounds(bv_context_view((struct bv_context *)request.view_ctx),
	request.factor, request.min, request.max);
    return BRLCAD_OK;
}

/*
 * Auto-adjust the view so that geometry is framed within the view.  By
 * default all displayed geometry is framed, but if a list of objects
 * (or full paths) is supplied the view is instead centered and sized to
 * frame only those objects.
 *
 * Usage:
 * autoview [options] [object ...]
 *
 * TODO - like draw2, this needs to be aware of whether we're using local or shared
 * grp sets - if we're adaptive plotting, for example.
 *
 */
extern "C" int
ged_autoview2_core(struct ged *gedp, int argc, const char *argv[])
{
    static const char *usage = "[options] [object ...]";
    struct bu_vls cvls = BU_VLS_INIT_ZERO;

    /* default, 0.5 model scale == 2.0 view factor */
    fastf_t factor = BV_AUTOVIEW_SCALE_DEFAULT;

    GED_CHECK_DRAWABLE(gedp, BRLCAD_ERROR);
    GED_CHECK_VIEW(gedp, BRLCAD_ERROR);
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);


    int all_view_objs = 0;
    int print_help = 0;
    fastf_t scale = -1.0;
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);

    struct bu_opt_desc d[5];
    BU_OPT(d[0], "h", "help",      "",        NULL,     &print_help, "Print help and exit");
    BU_OPT(d[1], "",   "all-objs", "",        NULL,  &all_view_objs, "Bound all non-faceplate view features");
    BU_OPT(d[2], "s", "scale",  "#", &bu_opt_fastf_t,         &scale, "Set view scale (model scale relative to view size)");
    BU_OPT(d[3], "V", "view",  "name", &bu_opt_vls,           &cvls, "Specify view to adjust");
    BU_OPT_NULL(d[4]);

    argc-=(argc>0); argv+=(argc>0); /* skip command name argv[0] */

    int opt_ret = bu_opt_parse(NULL, argc, argv, d);

    if (print_help) {
	_ged_cmd_help(gedp, usage, d);
	return BRLCAD_OK;
    }

    argc = opt_ret;

    if (bu_vls_strlen(&cvls)) {
	view_ctx = ged_view_find_ctx(gedp, bu_vls_cstr(&cvls));
	if (!view_ctx) {
	    bu_vls_printf(gedp->ged_result_str, "Specified view %s not found\n", bu_vls_cstr(&cvls));
	    bu_vls_free(&cvls);
	    return BRLCAD_ERROR;
	}
    }
    bu_vls_free(&cvls);


    /* Backward compatibility: a lone leading numeric argument historically
     * set the view scale.  Honor that only when -s/--scale was not supplied
     * and the token does not name an existing object, so an object with a
     * numeric name still wins. */
    if (scale < 0.0 && argc > 0) {
	double sval;
	if (_autoview_arg_is_num(argv[0], &sval) &&
	    (!gedp->dbip || db_lookup(gedp->dbip, argv[0], LOOKUP_QUIET) == RT_DIR_NULL)) {
	    scale = sval;
	    argc--;
	    argv++;
	}
    }

    if (scale > 0.0)
	factor = 1.0 / scale;

    autoview_request request = {gedp, view_ctx, factor, argc > 0,
	all_view_objs, VINIT_ZERO, VINIT_ZERO};
    if (request.named_objects && !_autoview_named_database_bounds(gedp,
		argc, argv, request.min, request.max)) {
	bu_vls_printf(gedp->ged_result_str,
	    "Unable to determine finite autoview bounds for the requested objects");
	return BRLCAD_ERROR;
    }

    /* Hosts can execute GED commands on a worker.  The retained scene and
     * progressive provider share camera state with the view's owner thread. */
    return ged_exec_owner_invoke(gedp, _autoview_on_owner, &request);
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
