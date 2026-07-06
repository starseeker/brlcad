/*                         D R A W . C
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
/** @file libged/draw.c
 *
 * Draw command registration and compatibility shims.
 *
 * The implementation is intentionally small: normal draw/redraw command
 * semantics are owned by draw2.cpp and evaluated wire display is owned by
 * bigE.c's Obol/librt-backed provider path.  The historical _ged_drawtrees()
 * walker is no longer a live fallback for migrated draw behavior.
 */

#include "common.h"

#include "ged/draw.h"
#include "ged.h"

#include "../bsg_ged_draw_private.h"
#include "./ged_draw.h"


static int ged_draw_test_force_face_set_failure = 0;


void
ged_draw_test_force_primitive_face_set_failure(int enable)
{
    ged_draw_test_force_face_set_failure = enable ? 1 : 0;
}


int
ged_draw_test_primitive_face_set_failure_enabled(void)
{
    return ged_draw_test_force_face_set_failure;
}


extern int ged_draw2_core(struct ged *gedp, int argc, const char *argv[]);

int
ged_draw_core(struct ged *gedp, int argc, const char *argv[])
{
    return ged_draw2_core(gedp, argc, argv);
}


int
ged_ev_core(struct ged *gedp, int argc, const char *argv[])
{
    return ged_eval_wire_display_core(gedp, argc, argv, 1);
}


extern int ged_redraw2_core(struct ged *gedp, int argc, const char *argv[]);

int
ged_redraw_core(struct ged *gedp, int argc, const char *argv[])
{
    return ged_redraw2_core(gedp, argc, argv);
}


#include "../include/plugin.h"

#define GED_DRAW_COMMANDS(X, XID) \
    X(draw, ged_draw_core, GED_CMD_DEFAULT) \
    XID(E_cmd, "E", ged_E_core, GED_CMD_DEFAULT) \
    X(e, ged_draw_core, GED_CMD_DEFAULT) \
    X(ev, ged_ev_core, GED_CMD_DEFAULT) \
    X(redraw, ged_redraw_core, GED_CMD_DEFAULT) \
    X(loadview, ged_loadview_core, GED_CMD_DEFAULT) \
    X(preview, ged_preview_core, GED_CMD_DEFAULT)

GED_DECLARE_COMMAND_SET(GED_DRAW_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_draw", 1, GED_DRAW_COMMANDS)

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
