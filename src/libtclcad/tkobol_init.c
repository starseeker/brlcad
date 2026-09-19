/*                  T K O B O L _ I N I T . C
 * BRL-CAD
 *
 * Copyright (c) 2004-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/**
 * Tcl initialization for the Tk Obol host and image utilities.
 */

#include "common.h"

#include <stdlib.h>

#include "tcl.h"
#include "imgstream/fb_compat.h"
#include "bu/cmd.h"
#include "bu/log.h"
#include "bu/vls.h"
#include "tclcad.h"
#include "brlcad_version.h"
#include "./tclcad_private.h"

#ifdef HAVE_TKOBOL_HOST
extern int BrlcadTkObolHost_Init(Tcl_Interp *interp);
#endif

static int
tclcad_image_file_size(ClientData clientData, int argc, const char **argv)
{
    Tcl_Interp *interp = (Tcl_Interp *)clientData;
    size_t width, height;
    int pixel_size = 3;

    if (argc != 2 && argc != 3) {
	bu_log("wrong #args: should be \" fileName [#bytes/pixel]\"");
	return TCL_ERROR;
    }

    if (argc >= 3)
	pixel_size = atoi(argv[2]);

    if (imgstream_image_file_size(&width, &height, argv[1], pixel_size) > 0) {
	struct bu_vls vls = BU_VLS_INIT_ZERO;
	bu_vls_printf(&vls, "%lu %lu", (unsigned long)width,
		      (unsigned long)height);
	Tcl_SetObjResult(interp, Tcl_NewStringObj(bu_vls_addr(&vls),
		bu_vls_strlen(&vls)));
	bu_vls_free(&vls);
	return TCL_OK;
    }

    Tcl_SetObjResult(interp, Tcl_NewStringObj("0 0", 3));
    return TCL_OK;
}

static int
tclcad_command_wrapper(ClientData data, Tcl_Interp *interp, int argc,
	const char *argv[])
{
    struct bu_cmdtab *command = (struct bu_cmdtab *)data;
    return command->ct_func(interp, argc, argv);
}

static void
tclcad_register_commands(Tcl_Interp *interp, struct bu_cmdtab *commands)
{
    for (struct bu_cmdtab *command = commands;
	 command->ct_name != (char *)NULL; command++) {
	(void)Tcl_CreateCommand(interp, command->ct_name,
	    tclcad_command_wrapper, (ClientData)command,
	    (Tcl_CmdDeleteProc *)NULL);
    }
}

TCLCAD_EXPORT int
tclcad_tkobol_init(Tcl_Interp *interp)
{
    static struct bu_cmdtab commands[] = {
	{"fb_common_file_size", tclcad_image_file_size},
	{(const char *)NULL, BU_CMD_NULL}
    };

    tclcad_register_commands(interp, commands);

#ifdef HAVE_TKOBOL_HOST
    if (!tclcad_obol_host_factories_register())
	return TCL_ERROR;
    if (BrlcadTkObolHost_Init(interp) != TCL_OK)
	return TCL_ERROR;
#endif

    Tcl_PkgProvide(interp, "BrlcadTkObol", brlcad_version());
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
