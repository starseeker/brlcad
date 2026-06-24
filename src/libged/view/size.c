/*                         S I Z E . C
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
/** @file libged/size.c
 *
 * The size command.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "rt/view.h"

#include "../ged_private.h"
#include "./ged_view.h"


int
ged_size_core(struct ged *gedp, int argc, const char *argv[])
{
    /* intentionally double for scan */
    double size;
    static const char *usage = "[s]";

    GED_CHECK_VIEW(gedp, BRLCAD_ERROR);
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    void *view_ctx = ged_view_active_ctx(gedp);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* get view size */
    if (argc == 1) {
	fastf_t view_size = ged_view_context_size_get(view_ctx);
	if (gedp->dbip) {
	    bu_vls_printf(gedp->ged_result_str, "%g",
		    view_size * gedp->dbip->dbi_base2local);
	} else {
	    bu_vls_printf(gedp->ged_result_str, "%g", view_size);
	}
	return BRLCAD_OK;
    }

    /* set view size */
    if (argc == 2) {
	if (sscanf(argv[1], "%lf", &size) != 1
	    || size <= 0
	    || ZERO(size))
	{
	    bu_vls_printf(gedp->ged_result_str, "bad size - %s", argv[1]);
	    return BRLCAD_ERROR;
	}

	fastf_t view_size = (gedp->dbip) ? gedp->dbip->dbi_local2base * size : size;
	if (view_size < RT_VIEW_MIN_SIZE)
	    view_size = RT_VIEW_MIN_SIZE;
	ged_view_context_size_set(view_ctx, view_size);
	ged_view_context_update(view_ctx);

	return BRLCAD_OK;
    }

    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
    return BRLCAD_ERROR;
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
