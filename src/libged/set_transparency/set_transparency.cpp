/*                         S E T _ T R A N S P A R E N C Y . C
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
/** @file libged/set_transparency.c
 *
 * The set_transparency command.
 *
 */

#include "common.h"

#include "ged/draw.h"
#include "../ged_private.h"

/*
 * Set the transparency of the specified object
 *
 * Usage:
 * set_transparency obj tr
 *
 */
extern "C" int
ged_set_transparency_core(struct ged *gedp, int argc, const char *argv[])
{
    /* intentionally double for scan */
    double transparency;

    static const char *usage = "node tval";

    GED_CHECK_DRAWABLE(gedp, BRLCAD_ERROR);
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }


    if (argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    if (sscanf(argv[2], "%lf", &transparency) != 1) {
	bu_vls_printf(gedp->ged_result_str, "dgo_set_transparency: bad transparency - %s\n", argv[2]);
	return BRLCAD_ERROR;
    }

    struct ged_scene_path_request request;
    ged_scene_path_request_init(&request);
    request.path = argv[1];
    request.match = GED_SCENE_PATH_MATCH_SUBTREE;
    if (ged_scene_opacity_set(gedp, &request, 1.0 - transparency, NULL) !=
	GED_SCENE_OK) {
	bu_vls_printf(gedp->ged_result_str,
		"opacity must be between 0 and 1\n");
	return BRLCAD_ERROR;
    }

    return BRLCAD_OK;
}


#include "../include/plugin.h"

#define GED_SET_TRANSPARENCY_COMMANDS(X, XID) \
    X(set_transparency, ged_set_transparency_core, GED_CMD_DEFAULT) \

GED_DECLARE_COMMAND_SET(GED_SET_TRANSPARENCY_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_set_transparency", 1, GED_SET_TRANSPARENCY_COMMANDS)

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
