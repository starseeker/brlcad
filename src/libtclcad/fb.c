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
/** @addtogroup libstruct fb */
/** @{ */
/**
 *
 * A framebuffer object contains the attributes and
 * methods for controlling framebuffers.
 *
 */
/** @} */

#include "common.h"

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
#include "../libdm/include/private.h"
#include "dm.h"
#include "dm/fbserv_legacy.h"
#include "ged/view.h"
#include "tclcad.h"
#include "./tclcad_private.h"
#include "./view/view.h"

void
to_fbs_callback(void *clientData)
{
    to_refresh_view(clientData);
}


int
to_close_fbs(void *view_ctx)
{
    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(view_ctx);
    if (!tvd)
	return TCL_ERROR;

    struct fb *fbp = fbs_legacy_framebuffer(&tvd->gdv_fbs);
    if (fbp == FB_NULL)
	return TCL_OK;

    fb_flush(fbp);
    fb_close_existing(fbp);
    dm_fbserv_set_framebuffer(&tvd->gdv_fbs, FB_NULL);

    return TCL_OK;
}


/*
 * Open/activate the display managers framebuffer.
 */
int
to_open_fbs(void *view_ctx, Tcl_Interp *interp)
{
    /* already open */
    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(view_ctx);
    if (!tvd)
	return TCL_ERROR;

    if (fbs_legacy_framebuffer(&tvd->gdv_fbs) != FB_NULL)
	return TCL_OK;

    dm_fbserv_set_framebuffer(&tvd->gdv_fbs, FB_NULL);

    if (fbs_legacy_framebuffer(&tvd->gdv_fbs) == FB_NULL) {
	Tcl_Obj *obj;

	obj = Tcl_GetObjResult(interp);
	if (Tcl_IsShared(obj))
	    obj = Tcl_DuplicateObj(obj);

	Tcl_AppendStringsToObj(obj, "openfb: failed to allocate framebuffer memory\n", (char *)NULL);

	Tcl_SetObjResult(interp, obj);
	return TCL_ERROR;
    }

    return TCL_OK;
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

    void *view_ctx = ged_view_find_ctx(gedp, argv[1]);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(view_ctx);
    if (!tvd)
	return BRLCAD_ERROR;

    /* Get fb mode */
    if (argc == 2) {
	bu_vls_printf(gedp->ged_result_str, "%d", tvd->gdv_fbs.fbs_mode);
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

    tvd->gdv_fbs.fbs_mode = mode;
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

    void *view_ctx = ged_view_find_ctx(gedp, argv[1]);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(view_ctx);
    if (!tvd)
	return BRLCAD_ERROR;

    if (fbs_legacy_framebuffer(&tvd->gdv_fbs) == FB_NULL) {
	bu_vls_printf(gedp->ged_result_str, "%s listen: framebuffer not open!\n", argv[0]);
	return BRLCAD_ERROR;
    }

    /* return the port number */
    if (argc == 2) {
	const char *ipc_env = fbs_ipc_child_addr_env(&tvd->gdv_fbs);
	if (ipc_env) {
	    const char *eq = strchr(ipc_env, '=');
	    const char *addr = eq ? eq + 1 : ipc_env;
	    bu_vls_printf(gedp->ged_result_str, "%s", addr);
	} else {
	    bu_vls_printf(gedp->ged_result_str, "%d", fbs_listener_port(&tvd->gdv_fbs));
	}
	return BRLCAD_OK;
    }

    if (argc == 3) {
	if (BU_STR_EQUAL(argv[2], "ipc")) {
	    if (fbs_listener_port(&tvd->gdv_fbs) >= 0)
		fbs_close(&tvd->gdv_fbs);
	    if (tclcad_listen_ipc(&tvd->gdv_fbs, (Tcl_Interp *)gedp->ged_interp) != BRLCAD_OK) {
		bu_vls_printf(gedp->ged_result_str, "listen: failed to start IPC listener\n");
		return BRLCAD_ERROR;
	    }
	    {
		const char *ipc_env = fbs_ipc_child_addr_env(&tvd->gdv_fbs);
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
	    // Set up fbo_fbs callbacks, then call fbs_open
	    tclcad_fbserv_set_transport(&tvd->gdv_fbs);
	    fbs_open(&tvd->gdv_fbs, port);
	} else {
	    fbs_close(&tvd->gdv_fbs);
	}
	{
	    const char *ipc_env = fbs_ipc_child_addr_env(&tvd->gdv_fbs);
	    if (ipc_env) {
		const char *eq = strchr(ipc_env, '=');
		const char *addr = eq ? eq + 1 : ipc_env;
		bu_vls_printf(gedp->ged_result_str, "%s", addr);
	    } else {
		bu_vls_printf(gedp->ged_result_str, "%d", fbs_listener_port(&tvd->gdv_fbs));
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
