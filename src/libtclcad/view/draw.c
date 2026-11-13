/*                       D R A W . C
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
/** @file libtclcad/view/draw.c
 *
 */
/** @} */

#include "common.h"
#include "bv.h"
#include "ged.h"
#include "ged/scene.h"
#include "ged/draw.h"
#include "ged/view.h"
#include "tclcad.h"

/* Private headers */
#include "../tclcad_private.h"

#include "../view/view.h"



void
go_draw(struct ged_view_context *view_ctx)
{
    if (!view_ctx || !current_top || !current_top->to_gedp)
	return;

    struct ged_scene_redraw_request request;
    ged_scene_redraw_request_init(&request);
    request.view = view_ctx;
    (void)ged_scene_redraw(current_top->to_gedp, &request, NULL);
}


int
to_edit_redraw(struct ged *gedp,
	       int argc,
	       const char *argv[])
{
    if (argc != 2)
	return BRLCAD_ERROR;

    if (!gedp || !argv[1] || !argv[1][0])
	return BRLCAD_ERROR;

    /* Redraw is a semantic retained-source operation.  Supplying the edited
     * path lets the scene reducer match every owning draw intent without
     * exporting renderer records or re-executing draw commands. */
    const char *paths[] = {argv[1]};
    struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
    size_t vi;
    for (vi = 0; vi < BU_PTBL_LEN(views); vi++) {
	struct ged_view_context *view_ctx =
	    (struct ged_view_context *)BU_PTBL_GET(views, vi);
	struct ged_scene_redraw_request request;
	ged_scene_redraw_request_init(&request);
	request.view = view_ctx;
	request.paths = paths;
	request.path_count = 1;
	(void)ged_scene_redraw(gedp, &request, NULL);
    }

    to_refresh_all_views(current_top);
    return BRLCAD_OK;
}

int
to_redraw(struct ged *gedp,
	  int argc,
	  const char *argv[],
	  ged_func_ptr UNUSED(func),
	  const char *usage,
	  int UNUSED(maxargs))
{
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 2) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    return to_edit_redraw(gedp, argc, argv);
}

int
to_blast(struct ged *gedp,
	 int argc,
	 const char *argv[],
	 ged_func_ptr UNUSED(func),
	 const char *UNUSED(usage),
	 int UNUSED(maxargs))
{
    int ret;

    ret = ged_exec(gedp, argc, argv);

    if (ret != BRLCAD_OK)
	return ret;

    to_autoview_all_views(current_top);

    return ret;
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
