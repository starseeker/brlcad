/*                            F B . C
 * BRL-CAD
 *
 * Copyright (c) 1997-2026 United States Government as represented by
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
/**
 *
 * A framebuffer object contains the attributes and
 * methods for controlling framebuffers.
 *
 */
/** @} */

#include "common.h"

#include "ged/display_obol_private.h"

#include <stdlib.h>
#include <string.h>

#ifdef HAVE_STRINGS_H
# include <strings.h>
#endif

#include "tcl.h"
#include "bu/cmd.h"
#include "bu/color.h"
#include "bu/getopt.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "BObol/BDisplayEndpoint.h"
#include "ged/display.h"
#include "ged/view.h"
#include "imgstream/fbserv.h"
#include "tclcad.h"
#include "./tclcad_private.h"
#include "./view/view.h"

void
to_fbs_callback(void *clientData)
{
    to_refresh_view((struct ged_view_context *)clientData);
}

static void
tclcad_fbs_error(Tcl_Interp *interp, const char *message)
{
    if (!interp)
	return;

    Tcl_Obj *obj = Tcl_GetObjResult(interp);
    if (Tcl_IsShared(obj))
	obj = Tcl_DuplicateObj(obj);
    Tcl_AppendStringsToObj(obj, message, (char *)NULL);
    Tcl_SetObjResult(interp, obj);
}


/* Bind the one GED session stream to the requested view's Obol endpoint.  Factory
 * instances such as TkObolEndpointHost are intentionally opaque and are not
 * BObolWindowHost instances.  The libged bridge owns its generic image
 * host and publishes that stream through the endpoint capture provider.
 * TclCAD owns only notifier callbacks. */
static int
tclcad_obol_framebuffer_bind(struct ged_view_context *view_ctx, Tcl_Interp *interp)
{
    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(view_ctx);
    if (!tvd || !tvd->gedp)
	return TCL_ERROR;

    bobol_display_endpoint_t *endpoint =
	ged_view_context_obol_endpoint_get(view_ctx);
    if (!endpoint)
	return TCL_OK;

    if (ged_view_framebuffer_backend_ensure(tvd->gedp,
	view_ctx) != BRLCAD_OK) {
	tclcad_fbs_error(interp,
	    "openfb: unable to attach the Obol framebuffer stream\n");
	return TCL_ERROR;
    }

    struct fbserv_obj *fbsp = tvd->gedp->ged_fbs;
    if (!fbsp)
	return TCL_ERROR;
    fbsp->fbs_callback = to_fbs_callback;
    fbsp->fbs_clientData = view_ctx;
    fbsp->fbs_interp = interp;

    /* The libged Obol bridge owns its IPC transport.  In particular, its
     * worker keeps a synchronous external renderer from blocking the Tcl
     * event loop that will later present the updated image. */
    return TCL_OK;
}


int
to_close_fbs(struct ged_view_context *view_ctx)
{
    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(view_ctx);
    if (!tvd || !tvd->gedp)
	return TCL_ERROR;

    struct fbserv_obj *fbsp = tvd->gedp->ged_fbs;
    if (!fbsp || fbsp->fbs_clientData != view_ctx)
	return TCL_OK;

    ged_view_framebuffer_release(tvd->gedp);
    fbsp->fbs_callback = NULL;
    fbsp->fbs_clientData = NULL;
    fbsp->fbs_interp = NULL;

    return TCL_OK;
}


/*
 * Open or rebind the view's endpoint-owned Obol framebuffer stream.
 */
int
to_open_fbs(struct ged_view_context *view_ctx, Tcl_Interp *interp)
{
    return tclcad_obol_framebuffer_bind(view_ctx, interp);
}



int
to_set_fb_mode(struct ged *gedp,
	       int argc,
	       const char *argv[],
	       ged_func_ptr UNUSED(func),
	       const char *usage,
	       int UNUSED(maxargs))
{
    int mode;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (3 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    struct ged_view_context *view_ctx = ged_view_find_ctx(gedp, argv[1]);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (!tclcad_view_data_from_view_ctx(view_ctx))
	return BRLCAD_ERROR;

    /* TclCAD's public command ordering differs from libbv's historical
     * ordering: 0=off, 1=underlay, 2=interlay, 3=overlay.  Keep the
     * endpoint as the rendering authority and translate only at this API. */
    if (argc == 2) {
	struct bv_display_property_value value =
	    BV_DISPLAY_PROPERTY_VALUE_INIT;
	if (ged_view_context_display_property_get(view_ctx,
		"composition.framebuffer.mode", &value) !=
	    BV_DISPLAY_PROPERTY_OK || !value.string_value) {
	    bu_vls_printf(gedp->ged_result_str,
		"View has no Obol framebuffer composition policy");
	    return BRLCAD_ERROR;
	}
	if (BU_STR_EQUAL(value.string_value, "underlay"))
	    mode = TCLCAD_OBJ_FB_MODE_UNDERLAY;
	else if (BU_STR_EQUAL(value.string_value, "interlay"))
	    mode = TCLCAD_OBJ_FB_MODE_INTERLAY;
	else if (BU_STR_EQUAL(value.string_value, "overlay"))
	    mode = TCLCAD_OBJ_FB_MODE_OVERLAY;
	else if (BU_STR_EQUAL(value.string_value, "off"))
	    mode = TCLCAD_OBJ_FB_MODE_OFF;
	else {
	    bu_vls_printf(gedp->ged_result_str,
		"View has an invalid Obol framebuffer composition policy");
	    return BRLCAD_ERROR;
	}
	bu_vls_printf(gedp->ged_result_str, "%d", mode);
	return BRLCAD_OK;
    }

    /* Set fb mode */
    if (bu_sscanf(argv[2], "%d", &mode) != 1) {
	bu_vls_printf(gedp->ged_result_str, "set_fb_mode: bad value - %s\n", argv[2]);
	return BRLCAD_ERROR;
    }

    if (mode < 0)
	mode = 0;
    else if (TCLCAD_OBJ_FB_MODE_OVERLAY < mode)
	mode = TCLCAD_OBJ_FB_MODE_OVERLAY;

    const char *composition = mode == TCLCAD_OBJ_FB_MODE_UNDERLAY ?
	"underlay" : mode == TCLCAD_OBJ_FB_MODE_INTERLAY ? "interlay" :
	mode == TCLCAD_OBJ_FB_MODE_OVERLAY ? "overlay" : "off";
    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    value.type = BV_DISPLAY_PROPERTY_ENUM;
    value.string_value = composition;
    if (ged_view_context_display_property_set(view_ctx,
	    "composition.framebuffer.mode", &value) !=
	BV_DISPLAY_PROPERTY_OK) {
	bu_vls_printf(gedp->ged_result_str,
	    "Unable to set Obol framebuffer composition policy");
	return BRLCAD_ERROR;
    }

    to_refresh_view(view_ctx);

    return BRLCAD_OK;
}


int
to_listen(struct ged *gedp,
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

    if (3 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    struct ged_view_context *view_ctx = ged_view_find_ctx(gedp, argv[1]);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (!tclcad_view_data_from_view_ctx(view_ctx) ||
	to_open_fbs(view_ctx, (Tcl_Interp *)gedp->ged_interp) != TCL_OK ||
	!gedp->ged_fbs ||
	!fbs_framebuffer_backend_installed(gedp->ged_fbs)) {
	bu_vls_printf(gedp->ged_result_str,
	    "%s listen: Obol framebuffer stream is unavailable\n", argv[0]);
	return BRLCAD_ERROR;
    }
    struct fbserv_obj *fbsp = gedp->ged_fbs;

    /* return the port number */
    if (argc == 2) {
	const char *ipc_env = fbs_ipc_child_addr_env(fbsp);
	if (ipc_env) {
	    const char *eq = strchr(ipc_env, '=');
	    const char *addr = eq ? eq + 1 : ipc_env;
	    bu_vls_printf(gedp->ged_result_str, "%s", addr);
	} else {
	    bu_vls_printf(gedp->ged_result_str, "%d", fbs_listener_port(fbsp));
	}
	return BRLCAD_OK;
    }

    if (argc == 3) {
	if (BU_STR_EQUAL(argv[2], "ipc")) {
	    if (fbs_can_close(fbsp))
		fbs_close(fbsp);
	    if (tclcad_listen_ipc(fbsp, (Tcl_Interp *)gedp->ged_interp) != BRLCAD_OK) {
		bu_vls_printf(gedp->ged_result_str, "listen: failed to start IPC listener\n");
		return BRLCAD_ERROR;
	    }
	    {
		const char *ipc_env = fbs_ipc_child_addr_env(fbsp);
		if (ipc_env) {
		    const char *eq = strchr(ipc_env, '=');
		    const char *addr = eq ? eq + 1 : ipc_env;
		    bu_vls_printf(gedp->ged_result_str, "%s", addr);
		}
	    }
	    return BRLCAD_OK;
	}

	int port;

	if (bu_sscanf(argv[2], "%d", &port) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "listen: bad value - %s\n", argv[2]);
	    return BRLCAD_ERROR;
	}

	if (port >= 0) {
	    if (fbs_open(fbsp, port) != BRLCAD_OK) {
		bu_vls_printf(gedp->ged_result_str,
		    "listen: failed to open framebuffer listener\n");
		return BRLCAD_ERROR;
	    }
	} else {
	    fbs_close(fbsp);
	}
	{
	    const char *ipc_env = fbs_ipc_child_addr_env(fbsp);
	    if (ipc_env) {
		const char *eq = strchr(ipc_env, '=');
		const char *addr = eq ? eq + 1 : ipc_env;
		bu_vls_printf(gedp->ged_result_str, "%s", addr);
	    } else {
		bu_vls_printf(gedp->ged_result_str, "%d", fbs_listener_port(fbsp));
	    }
	}
	return BRLCAD_OK;
    }

    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
    return BRLCAD_ERROR;
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
