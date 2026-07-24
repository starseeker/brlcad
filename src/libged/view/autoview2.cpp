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

#include "BObol/BDatabaseSource.h"
#include "BObol/BMeasureAction.h"
#include "BObol/BSceneController.h"
#include "BObol/BViewController.h"
#include <Inventor/SoViewport.h>
#include "bu/opt.h"
#include "bv.h"
#include "ged/draw.h"
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
_autoview_bobol_database_bounds(struct ged *gedp, vect_t *min, vect_t *max)
{
    int empty = 1;
    /* Explicit autoview command: allow the precise per-child member-bounds
     * fallback when no cheap source bounds exist yet.  (Once geometry is
     * realized the cheap source bounds are used regardless.) */
    return ged_draw_obol_scene_database_autoview_bounds(gedp, min, max,
	&empty, 1) && !empty;
}

static int
_autoview_obol_database_scene(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	fastf_t factor,
	int all_view_objs)
{
    if (all_view_objs)
	return 0;

    vect_t min, max;
    if (!_autoview_bobol_database_bounds(gedp, &min, &max))
	return 0;

    bv_autoview_bounds(bv_context_view((struct bv_context *)view_ctx),
	    factor, min, max);
    return 1;
}

static int
_autoview_bobol_view_scene(struct ged_view_context *view_ctx,
	fastf_t factor)
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
    vect_t min, max;
    VSET(min, bmin[0], bmin[1], bmin[2]);
    VSET(max, bmax[0], bmax[1], bmax[2]);
    bv_autoview_bounds(bv_context_view((struct bv_context *)view_ctx),
	factor, min, max);
    return 1;
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

    if (argc > 0) {
	/* Frame only the named objects.  Bound them directly from the
	 * database (they need not be displayed). */
	point_t min, max;
	if (rt_obj_bounds(gedp->ged_result_str, gedp->dbip, argc, argv, 0, min, max) != BRLCAD_OK)
	    return BRLCAD_ERROR;
	bv_autoview_bounds(bv_context_view((struct bv_context *)view_ctx),
		factor, min, max);
    } else {
	if (all_view_objs)
	    (void)_autoview_bobol_view_scene(view_ctx, factor);
	else
	    (void)_autoview_obol_database_scene(gedp, view_ctx, factor, 0);
	/* A deferred Obol root starts with a proxy and gains its full bounds as
	 * realization completes.  Continue fitting this command's view until
	 * that work settles, unless a later view change supersedes it. */
	if (!all_view_objs)
	    (void)ged_draw_obol_progressive_autoview_follow(gedp, view_ctx,
		    factor);
    }

    return BRLCAD_OK;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
