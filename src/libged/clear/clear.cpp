/*                       C L E A R . C P P
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
/** @file libged/clear.cpp
 *
 * The clear command.
 *
 */

#include "common.h"

#include "../ged_private.h"

// clear is a no-op at the libged level - its purpose is to support ged_exec
// callback hooks that call any application callback functions assigned to this
// command.
extern "C" int
ged_clear_core(struct ged *gedp, int argc, const char **argv)
{
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);
    bu_vls_trunc(gedp->ged_result_str, 0);
    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", argv[0]);
	return GED_HELP;
    }
    return BRLCAD_OK;
}

#include "../include/plugin.h"

#define GED_CLEAR_COMMANDS(X, XID) \
    X(clear, ged_clear_core, GED_CMD_DEFAULT) \

GED_DECLARE_COMMAND_SET(GED_CLEAR_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_clear", 1, GED_CLEAR_COMMANDS)

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
