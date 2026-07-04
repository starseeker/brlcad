/*                      G O B J S . C P P
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
/** @file libged/view/gobjs.c
 *
 * Create and manage transient view features that are used to
 * interactively edit geometry parameters for database solids
 * and combs.  These edits change only the locally stored info,
 * without out disturbing the original object in the .g, until
 * the changes are specifically written back to that object
 * with an explicit command.
 *
 */

#include "common.h"

#include "bu/cmd.h"
#include "bu/str.h"
#include "bu/vls.h"

#include "ged/draw.h"
#include "ged/view.h"
#include "./ged_view.h"

int
_gobjs_cmd_create(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    void *view_ctx = gd->cv;
    const char *usage_string = "view gobjs name create";
    const char *purpose_string = "create an editing view obj from a database solid/comb";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (argc != 2) {
	bu_vls_printf(gedp->ged_result_str, "view gobjs create g_obj_name view_obj_name\n");
	return BRLCAD_ERROR;
    }
    gd->vobj = argv[0];

    if (!ged_draw_view_context_gobject_create(gedp, view_ctx, gd->vobj,
	    argv[1], gedp->ged_result_str))
	return BRLCAD_ERROR;

    // TODO - set the object callbacks

    // TODO - in principle, we should be sharing a lot of logic with the edit command -
    // the only difference is we won't be unpacking and writing the rt_db_internal.

    return BRLCAD_OK;
}

int
_gobjs_cmd_delete(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view gobjs name delete";
    const char *purpose_string = "delete view feature";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    if (!gd->vobj) {
	bu_vls_trunc(gedp->ged_result_str, 0);
	bu_vls_printf(gedp->ged_result_str, "No view feature named \n");
	return BRLCAD_ERROR;
    }

    const char *nargv[5] = {"view", "obj", "remove", gd->vobj, NULL};
    return _view_cmd_objs(bs, 4, nargv);
}

extern "C" int
_view_cmd_gobjs(void *bs, int argc, const char **argv)
{
    const char *usage_string = "view [options] gobjs [create|del] ...";
    const char *purpose_string = "deprecated alias for 'view obj -g <dbpath> ...'";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    if (!gd || !gd->gedp)
	return BRLCAD_ERROR;
    struct ged *gedp = gd->gedp;

    argc--; argv++; /* skip gobjs */
    if (argc <= 0) {
	const char *nargv[3] = {"view", "obj", "list"};
	return _view_cmd_objs(bs, 3, nargv);
    }

    if (BU_STR_EQUAL(argv[0], "create")) {
	if (argc != 3) {
	    bu_vls_printf(gedp->ged_result_str, "Usage: view gobjs create <dbpath> <name>\n");
	    return BRLCAD_ERROR;
	}
	const char *nargv[7] = {"view", "obj", "-g", argv[1], "create", argv[2], NULL};
	return _view_cmd_objs(bs, 6, nargv);
    }

    if (BU_STR_EQUAL(argv[0], "del") || BU_STR_EQUAL(argv[0], "delete") || BU_STR_EQUAL(argv[0], "remove")) {
	if (argc != 2) {
	    bu_vls_printf(gedp->ged_result_str, "Usage: view gobjs del <name>\n");
	    return BRLCAD_ERROR;
	}
	const char *nargv[5] = {"view", "obj", "remove", argv[1], NULL};
	return _view_cmd_objs(bs, 4, nargv);
    }

    bu_vls_printf(gedp->ged_result_str, "view gobjs is deprecated - use 'view obj -g <dbpath> create <name>' or 'view obj remove <name>'\n");
    return BRLCAD_ERROR;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
