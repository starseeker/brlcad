/*                         P E R S P E C T I V E . C
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
/** @file libged/perspective.c
 *
 * The perspective command.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "bv.h"

#include "../ged_private.h"


int
ged_perspective_core(struct ged *gedp, int argc, const char *argv[])
{
    /* intentionally double for scan */
    double perspective;
    static const char *usage = "[angle]";

    GED_CHECK_VIEW(gedp, BRLCAD_ERROR);
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    void *view_ctx = ged_view_active_ctx(gedp);
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* get the perspective angle */
    if (argc == 1) {
	struct bobol_endpoint_property_value value =
	    BOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	if (ged_view_context_display_property_get(view_ctx, "view.perspective",
		&value) != BOBOL_ENDPOINT_PROPERTY_OK) {
	    bu_vls_printf(gedp->ged_result_str,
		"active view has no Obol perspective policy");
	    return BRLCAD_ERROR;
	}
	bu_vls_printf(gedp->ged_result_str, "%g", value.double_value);
	return BRLCAD_OK;
    }

    /* set perspective angle */
    if (argc == 2) {
	if (sscanf(argv[1], "%lf", &perspective) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "bad perspective angle - %s", argv[1]);
	    return BRLCAD_ERROR;
	}

	struct bobol_endpoint_property_value value =
	    BOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	value.type = BOBOL_ENDPOINT_PROPERTY_DOUBLE;
	value.double_value = perspective;
	if (ged_view_context_display_property_set(view_ctx, "view.perspective",
		&value) != BOBOL_ENDPOINT_PROPERTY_OK) {
	    bu_vls_printf(gedp->ged_result_str,
		"bad perspective angle - %s", argv[1]);
	    return BRLCAD_ERROR;
	}

	return BRLCAD_OK;
    }

    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
    return BRLCAD_ERROR;
}


#include "../include/plugin.h"

#define GED_PERSPECTIVE_COMMANDS(X, XID) \
    X(perspective, ged_perspective_core, GED_CMD_DEFAULT) \

GED_DECLARE_COMMAND_SET(GED_PERSPECTIVE_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_perspective", 1, GED_PERSPECTIVE_COMMANDS)

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
