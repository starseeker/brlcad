/*                       F A C E D E F . C
 * BRL-CAD
 *
 * Copyright (c) 1986-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file mged/facedef.c
 *
 */

#include "common.h"
#include <string.h>
#include "ged/edit.h"
#include "rt/primitives/arb8.h"
#include "./mged.h"

int
f_facedef(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{
    struct cmdtab *ctp = (struct cmdtab *)clientData;
    MGED_CK_CMD(ctp);
    struct mged_state *s = ctp->s;
    if (s->global_editing_state != ST_S_EDIT || MEDIT(s)->es_int.idb_type != ID_ARB8) {
        Tcl_AppendResult(interp, "facedef requires an active ARB primitive edit", (char *)NULL);
        return TCL_ERROR;
    }
    if (argc < 3) {
        Tcl_AppendResult(interp, MORE_ARGS_STR,
            argc < 2 ? "Enter face: " : "Enter definition form (a, b, c, d, q): ", (char *)NULL);
        return TCL_ERROR;
    }
    bu_vls_trunc(s->gedp->ged_result_str, 0);
    int ret = ged_edit_session_command(s->gedp, argc, argv);
    if (ret != BRLCAD_OK) {
        const char *message = bu_vls_cstr(s->gedp->ged_result_str);
        if (strstr(message, "incorrect number of definition parameters"))
            Tcl_AppendResult(interp, MORE_ARGS_STR, "Enter remaining face definition parameters: ", (char *)NULL);
        else
            Tcl_AppendResult(interp, message, (char *)NULL);
        return TCL_ERROR;
    }
    replot_editing_solid(0, NULL, s, NULL);
    s->update_views = 1;
    return TCL_OK;
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
