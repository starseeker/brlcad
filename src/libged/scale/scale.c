/*                         S C A L E . C
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
/** @file libged/scale.c
 *
 * The scale command.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "rt/view.h"
#include "rt/view_legacy_bsg.h"

#include "../ged_private.h"


int
ged_scale_core(struct ged *gedp, int argc, const char *argv[])
{
    int ret;
    fastf_t sf1;
    fastf_t sf2;
    fastf_t sf3;

    if ((ret = ged_scale_args(gedp, argc, argv, &sf1, &sf2, &sf3)) != BRLCAD_OK)
	return ret;

    if (argc != 2) {
	bu_vls_printf(gedp->ged_result_str, "Can not scale xyz independently on a view.");
	return BRLCAD_ERROR;
    }

    if (sf1 <= SMALL_FASTF || INFINITY < sf1)
	return BRLCAD_OK;

    struct bsg_view *v = (struct bsg_view *)ged_view_active_ctx(gedp);

    /* scale the view */
    fastf_t view_scale = rt_view_scale_from_bsg(v) * sf1;
    if (view_scale < RT_VIEW_MIN_SIZE)
	view_scale = RT_VIEW_MIN_SIZE;

    rt_view_scale_set_bsg(v, view_scale);
    rt_view_update_bsg(v);

    return BRLCAD_OK;
}


#include "../include/plugin.h"

#define GED_SCALE_COMMANDS(X, XID) \
    X(scale, ged_scale_core, GED_CMD_DEFAULT) \
    X(sca, ged_scale_core, GED_CMD_DEFAULT) \

GED_DECLARE_COMMAND_SET(GED_SCALE_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_scale", 1, GED_SCALE_COMMANDS)

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
