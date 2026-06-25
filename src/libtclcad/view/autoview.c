/*                        A U T O V I E W . C
 * BRL-CAD
 *
 * Copyright (c) 2000-2026 United States Government as represented by
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
/** @addtogroup libtclcad */
/** @{ */
/** @file libtclcad/view/autoview.c
 *
 */
/** @} */

#include "common.h"
#include "bu/units.h"
#include "ged.h"
#include "ged/view.h"
#include "tclcad.h"

/* Private headers */
#include "../tclcad_private.h"
#include "../view/view.h"

void
to_autoview_view(void *view_ctx, const char *scale)
{
    int ret;
    const char *av[3];

    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(view_ctx);
    if (!tvd)
	return;

    ged_view_active_ctx_set(tvd->gedp, view_ctx);
    av[0] = "autoview";
    av[1] = scale;
    av[2] = NULL;

    if (scale)
	ret = ged_exec_autoview(tvd->gedp, 2, (const char **)av);
    else
	ret = ged_exec_autoview(tvd->gedp, 1, (const char **)av);

    if (ret == BRLCAD_OK) {
	if (0 < bu_vls_strlen(&tvd->gdv_callback)) {
	    Tcl_Eval(current_top->to_interp, bu_vls_addr(&tvd->gdv_callback));
	}

	to_refresh_view(view_ctx);
    }
}

int
to_autoview(struct ged *gedp,
	    int argc,
	    const char *argv[],
	    ged_func_ptr UNUSED(func),
	    const char *usage,
	    int UNUSED(maxargs))
{
    void *view_ctx;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (argc > 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s [scale]", argv[0], usage);
	return BRLCAD_ERROR;
    }

    view_ctx = ged_view_find_ctx(gedp, argv[1]);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (argc > 2)
	to_autoview_view(view_ctx, argv[2]);
    else
	to_autoview_view(view_ctx, NULL);

    return BRLCAD_OK;
}


void
to_autoview_all_views(struct tclcad_obj *top)
{
    struct bu_ptbl *views = ged_view_set_views_ctx(top->to_gedp);
    for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	void *view_ctx = BU_PTBL_GET(views, i);
	to_autoview_view(view_ctx, NULL);
    }
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
