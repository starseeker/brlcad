/*                         Q U A T . C
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
/** @file libged/quat.c
 *
 * The quat command.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "rt/view_legacy_bsg.h"

#include "../ged_private.h"
#include "./ged_view.h"


int
ged_quat_core(struct ged *gedp, int argc, const char *argv[])
{
    quat_t quat;
    double scan[4];
    static const char *usage = "a b c d";

    GED_CHECK_VIEW(gedp, BRLCAD_ERROR);
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    struct bsg_view *v = (struct bsg_view *)ged_view_active_ctx(gedp);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* return Viewrot as a quaternion */
    if (argc == 1) {
	rt_view_orientation_quat_from_bsg(quat, v);
	bu_vls_printf(gedp->ged_result_str, "%.12g %.12g %.12g %.12g", V4ARGS(quat));
	return BRLCAD_OK;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: view %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    /* Set the view orientation given a quaternion */
    if (sscanf(argv[1], "%lf", &scan[0]) != 1
	|| sscanf(argv[2], "%lf", &scan[1]) != 1
	|| sscanf(argv[3], "%lf", &scan[2]) != 1
	|| sscanf(argv[4], "%lf", &scan[3]) != 1)
    {
	bu_vls_printf(gedp->ged_result_str, "view %s: bad value detected - %s %s %s %s",
		      argv[0], argv[1], argv[2], argv[3], argv[4]);
	return BRLCAD_ERROR;
    }
    HMOVE(quat, scan);

    mat_t rotation;
    quat_quat2mat(rotation, quat);
    rt_view_rotation_set_bsg(v, rotation);
    rt_view_update_bsg(v);

    return BRLCAD_OK;
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
