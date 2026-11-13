/*                         L O O K A T . C
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
/** @file libged/lookat.c
 *
 * The lookat command.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "bv.h"

#include "../ged_private.h"


int
ged_lookat_core(struct ged *gedp, int argc, const char *argv[])
{
    point_t look;
    point_t eye;
    point_t tmp;
    point_t new_center;
    vect_t dir;
    vect_t view_aet;
    mat_t view2model;
    fastf_t view_scale;
    fastf_t new_az, new_el;
    double scan[3];
    static const char *usage = "x y z";
    double lbval = (gedp->dbip) ? gedp->dbip->dbi_local2base : 1.0;

    GED_CHECK_VIEW(gedp, BRLCAD_ERROR);
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 2 && argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    if (argc == 2) {
	if (bn_decode_vect(look, argv[1]) != 3) {
	    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	    return BRLCAD_ERROR;
	}
    } else {
	if (sscanf(argv[1], "%lf", &scan[X]) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "ged_lookat_core: bad X value - %s\n", argv[1]);
	    return BRLCAD_ERROR;
	}

	if (sscanf(argv[2], "%lf", &scan[Y]) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "ged_lookat_core: bad Y value - %s\n", argv[2]);
	    return BRLCAD_ERROR;
	}

	if (sscanf(argv[3], "%lf", &scan[Z]) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "ged_lookat_core: bad Z value - %s\n", argv[3]);
	    return BRLCAD_ERROR;
	}

	/* convert from double to fastf_t */
	VMOVE(look, scan);
    }

    VSCALE(look, look, lbval);

    VSET(tmp, 0.0, 0.0, 1.0);
    bv_view2model_get(view2model, view);
    MAT4X3PNT(eye, view2model, tmp);

    VSUB2(dir, eye, look);
    VUNITIZE(dir);
    bn_ae_vec(&new_az, &new_el, dir);

    bv_aet_get(view_aet, view);
    VSET(view_aet, new_az, new_el, view_aet[Z]);
    bv_aet_set(view, view_aet);

    view_scale = bv_scale_get(view);
    VJOIN1(new_center, eye, -view_scale, dir);
    bv_center_set(view, new_center);

    ged_view_context_update(view_ctx);

    return BRLCAD_OK;
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
