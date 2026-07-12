/*                       R E F R E S H . C
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
/** @file libtclcad/views/refresh.c
 *
 */
/** @} */

#include "common.h"
#include "brlobol/display_endpoint.h"
#include "bv.h"
#include "ged.h"
#include "ged/draw_obol.h"
#include "ged/view.h"
#include "tclcad.h"

/* Private headers */
#include "../tclcad_private.h"
#include "../view/view.h"

void
go_refresh_draw(struct ged *gedp, void *draw_view_ctx, int UNUSED(restore_zbuffer))
{
    if (!gedp || !draw_view_ctx)
	return;

    brlobol_display_endpoint_t *endpoint =
	ged_view_context_display_endpoint_get(draw_view_ctx);
    if (!endpoint)
	return;

    go_draw(draw_view_ctx);
    (void)ged_draw_obol_framebuffer_present(gedp);
    if (!brlobol_display_endpoint_view_sync(endpoint, draw_view_ctx))
	return;
    (void)brlobol_display_endpoint_request_frame(endpoint, "TclCAD refresh");
}

void
go_refresh(struct ged *gedp, void *draw_view_ctx)
{
    go_refresh_draw(gedp, draw_view_ctx, 0);
}

void
to_refresh_view(void *view_ctx)
{

    if (current_top == NULL || view_ctx == NULL)
	return;

    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    if (!view || !bv_refresh_enabled_get(view) ||
	    bv_refresh_suppressed_get(view))
	return;

    bv_refresh_request(view, GED_VIEW_REFRESH_ALL);
    if (to_is_viewable(view_ctx)) {
	(void)bv_refresh_consume(view);
	go_refresh(current_top->to_gedp, view_ctx);
	bv_refresh_complete(view);
    }
}

void
to_refresh_all_views(struct tclcad_obj *top)
{
    void *view_ctx;

    struct bu_ptbl *views = ged_view_set_views_ctx(top->to_gedp);
    for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	view_ctx = BU_PTBL_GET(views, i);
	to_refresh_view(view_ctx);
    }
}

int
to_refresh_all_enabled(struct tclcad_obj *top)
{
    void *view_ctx;

    if (!top)
	return 1;

    struct bu_ptbl *views = ged_view_set_views_ctx(top->to_gedp);
    for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	view_ctx = BU_PTBL_GET(views, i);
	const struct bv *view =
	    bv_context_view_const((const struct bv_context *)view_ctx);
	if (!bv_refresh_enabled_get(view))
	    return 0;
    }

    return 1;
}

void
to_refresh_all_set_enabled(struct tclcad_obj *top, int enabled)
{
    void *view_ctx;

    if (!top)
	return;

    struct bu_ptbl *views = ged_view_set_views_ctx(top->to_gedp);
    for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	view_ctx = BU_PTBL_GET(views, i);
	struct bv *view = bv_context_view((struct bv_context *)view_ctx);
	bv_refresh_enabled_set(view, enabled);
    }
}

void
to_refresh_suppress_all_begin(struct tclcad_obj *top)
{
    void *view_ctx;

    if (!top)
	return;

    struct bu_ptbl *views = ged_view_set_views_ctx(top->to_gedp);
    for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	view_ctx = BU_PTBL_GET(views, i);
	struct bv *view = bv_context_view((struct bv_context *)view_ctx);
	bv_refresh_suppress_begin(view);
    }
}

void
to_refresh_suppress_all_end(struct tclcad_obj *top)
{
    void *view_ctx;

    if (!top)
	return;

    struct bu_ptbl *views = ged_view_set_views_ctx(top->to_gedp);
    for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	view_ctx = BU_PTBL_GET(views, i);
	struct bv *view = bv_context_view((struct bv_context *)view_ctx);
	bv_refresh_suppress_end(view);
    }
}

int
to_refresh(struct ged *gedp,
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

    return to_handle_refresh(gedp, argv[1]);
}


int
to_refresh_all(struct ged *gedp,
	       int argc,
	       const char *argv[],
	       ged_func_ptr UNUSED(func),
	       const char *UNUSED(usage),
	       int UNUSED(maxargs))
{
    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s", argv[0]);
	return BRLCAD_ERROR;
    }

    to_refresh_all_views(current_top);

    return BRLCAD_OK;
}


int
to_refresh_on(struct ged *gedp,
	      int argc,
	      const char *argv[],
	      ged_func_ptr UNUSED(func),
	      const char *UNUSED(usage),
	      int UNUSED(maxargs))
{
    int on;
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (2 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s", argv[0]);
	return BRLCAD_ERROR;
    }

    /* Get refresh_on state */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "%d", to_refresh_all_enabled(current_top));
	return BRLCAD_OK;
    }

    /* Set refresh_on state */
    if (bu_sscanf(argv[1], "%d", &on) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s", argv[0]);
	return BRLCAD_ERROR;
    }

    to_refresh_all_set_enabled(current_top, on);

    return BRLCAD_OK;
}

int
to_handle_refresh(struct ged *gedp,
		  const char *name)
{
    void *view_ctx = ged_view_find_ctx(gedp, name);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", name);
	return BRLCAD_ERROR;
    }

    to_refresh_view(view_ctx);

    return BRLCAD_OK;
}


void
to_refresh_handler(void *clientdata)
{
    /* Possibly do more here */

    to_refresh_view(clientdata);
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
