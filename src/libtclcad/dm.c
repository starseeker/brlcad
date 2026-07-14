/*                            D M . C
 * BRL-CAD
 *
 * Copyright (c) 2004-2026 United States Government as represented by
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
/**
 *
 * Tcl initialization for the Tk Obol host and image utilities.
 *
 */

#include "common.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifdef HAVE_STRINGS_H
# include <strings.h>
#endif

#include "tcl.h"
#include "vmath.h"
#include "imgstream/fb_compat.h"
#include "bu/cmd.h"
#include "bu/log.h"
#include "bu/vls.h"
#include "tclcad.h"
/* private headers */
#include "brlcad_version.h"
#include "./tclcad_private.h"

#ifdef HAVE_TKOBOL_HOST
extern int BrlcadTkObolHost_Init(Tcl_Interp *interp);
#endif

/**
 * Hook function wrapper to the image-size Tcl command.
 */
int
fb_cmd_common_file_size(ClientData clientData, int argc, const char **argv)
{
    Tcl_Interp *interp = (Tcl_Interp *)clientData;
    size_t width, height;
    int pixel_size = 3;

    if (argc != 2 && argc != 3) {
	bu_log("wrong #args: should be \" fileName [#bytes/pixel]\"");
	return TCL_ERROR;
    }

    if (argc >= 3) {
	pixel_size = atoi(argv[2]);
    }

    if (imgstream_image_file_size(&width, &height, argv[1], pixel_size) > 0) {
	struct bu_vls vls = BU_VLS_INIT_ZERO;
	bu_vls_printf(&vls, "%lu %lu", (unsigned long)width, (unsigned long)height);
	Tcl_SetObjResult(interp,
			 Tcl_NewStringObj(bu_vls_addr(&vls), bu_vls_strlen(&vls)));
	bu_vls_free(&vls);
	return TCL_OK;
    }

    /* Signal error */
    char *zerr = bu_strdup("0 0");
    Tcl_SetResult(interp, zerr, TCL_STATIC);
    bu_free(zerr, "zerr");
    return TCL_OK;
}

static int
wrapper_func(ClientData data, Tcl_Interp *interp, int argc, const char *argv[])
{
    struct bu_cmdtab *ctp = (struct bu_cmdtab *)data;

    return ctp->ct_func(interp, argc, argv);
}

static void
register_cmds(Tcl_Interp *interp, struct bu_cmdtab *cmds)
{
    struct bu_cmdtab *ctp = NULL;

    for (ctp = cmds; ctp->ct_name != (char *)NULL; ctp++) {
	(void)Tcl_CreateCommand(interp, ctp->ct_name, wrapper_func, (ClientData)ctp, (Tcl_CmdDeleteProc *)NULL);
    }
}


TCLCAD_EXPORT int
Dm_Init(Tcl_Interp *interp)
{
    static struct bu_cmdtab cmdtab[] = {
	{"fb_common_file_size",	 fb_cmd_common_file_size},
	{(const char *)NULL, BU_CMD_NULL}
    };

    /* register commands */
    register_cmds(interp, cmdtab);

#ifdef HAVE_TKOBOL_HOST
    if (!tclcad_obol_host_factories_register())
	return TCL_ERROR;
    if (BrlcadTkObolHost_Init(interp) != TCL_OK)
	return TCL_ERROR;
#endif


    Tcl_PkgProvide(interp,  "Dm", brlcad_version());

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
