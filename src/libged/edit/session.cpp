/* BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
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

#include "common.h"
#include "ged/edit.h"
#include "../include/plugin.h"

#define GED_EDIT_SESSION_COMMANDS(X, XID) \
    X(oed, ged_session_oed, GED_CMD_DEFAULT) \
    X(oed_apply, ged_session_oed_apply, GED_CMD_DEFAULT) \
    X(oed_reset, ged_session_oed_reset, GED_CMD_DEFAULT) \
    X(orot, ged_session_orot, GED_CMD_DEFAULT) \
    X(translate, ged_session_translate, GED_CMD_DEFAULT) \
    X(rotobj, ged_session_rotobj, GED_CMD_DEFAULT) \
    X(qorot, ged_session_qorot, GED_CMD_DEFAULT) \
    X(M, ged_session_M, GED_CMD_DEFAULT) \
    X(extrude, ged_session_extrude, GED_CMD_DEFAULT) \
    X(mirface, ged_session_mirface, GED_CMD_DEFAULT) \
    X(permute, ged_session_permute, GED_CMD_DEFAULT) \
    X(edgedir, ged_session_edgedir, GED_CMD_DEFAULT) \
    X(facedef, ged_session_facedef, GED_CMD_DEFAULT) \
    X(rset, ged_session_rset, GED_CMD_DEFAULT) \
    X(knob, ged_session_knob, GED_CMD_DEFAULT) \
    X(mmenu_get, ged_session_mmenu_get, GED_CMD_DEFAULT) \
    X(sed, ged_session_sed, GED_CMD_DEFAULT) \
    X(accept, ged_session_accept, GED_CMD_DEFAULT) \
    X(reject, ged_session_reject, GED_CMD_DEFAULT) \
    X(sed_apply, ged_session_sed_apply, GED_CMD_DEFAULT) \
    X(sed_reset, ged_session_sed_reset, GED_CMD_DEFAULT) \
    X(press, ged_session_press, GED_CMD_DEFAULT) \
    X(status, ged_session_status, GED_CMD_DEFAULT) \
    X(get_sed, ged_session_get_sed, GED_CMD_DEFAULT) \
    X(get_sed_menus, ged_session_get_sed_menus, GED_CMD_DEFAULT) \
    X(p, ged_session_p, GED_CMD_DEFAULT) \

/* Command identity is also used for alias detection by the registry. */
#define GED_EDIT_SESSION_WRAPPER(name, function, options) \
    static int function(struct ged *gedp, int argc, const char **argv) \
    { return ged_edit_session_command(gedp, argc, argv); }
GED_EDIT_SESSION_COMMANDS(GED_EDIT_SESSION_WRAPPER, GED_EDIT_SESSION_WRAPPER)
#undef GED_EDIT_SESSION_WRAPPER

GED_DECLARE_COMMAND_SET(GED_EDIT_SESSION_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_edit_session", 1, GED_EDIT_SESSION_COMMANDS)

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
