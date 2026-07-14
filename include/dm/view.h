/*                         V I E W . H
 * BRL-CAD
 *
 * Copyright (c) 1993-2026 United States Government as represented by
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
/** @addtogroup libdm */
/** @{ */
/** @file dm/view.h
 *
 * dm routines related to view management (migrated from libtclcad.) These may
 * ultimately take a different form or move elsewhere - the immediate idea here
 * is to extract the key logic from libtclcad for use in non-Tcl environments.
 *
 */

#include "common.h"

#include <stdint.h>

#include "vmath.h"

#include "dm/defines.h"

#ifndef DM_VIEW_H
#define DM_VIEW_H

__BEGIN_DECLS

struct dm_path_edit_params {
    int edit_mode;
    double dx;
    double dy;
    mat_t edit_mat;
};

#define DM_VIEW_REFRESH_VIEW        0x00000001u
#define DM_VIEW_REFRESH_DRAW        0x00000002u
#define DM_VIEW_REFRESH_EDIT        0x00000004u
#define DM_VIEW_REFRESH_FRAMEBUFFER 0x00000008u
#define DM_VIEW_REFRESH_OVERLAY     0x00000010u
#define DM_VIEW_REFRESH_FORCE       0x80000000u
#define DM_VIEW_REFRESH_ALL         0xffffffffu

DM_EXPORT extern int dm_view_context_width_get(const void *view_ctx);
DM_EXPORT extern int dm_view_context_height_get(const void *view_ctx);
DM_EXPORT extern int dm_view_context_dimensions_set(void *view_ctx, int width, int height);
DM_EXPORT extern int dm_view_context_refresh_request(void *view_ctx, uint32_t flags);

DM_EXPORT extern void dm_draw_faceplate(void *view_ctx);

/* As a temporary measure, require client codes to specifically ask to enable
 * the bits that require librt in the headers if they're not going to be
 * calling them.  Not ideal, but pulling in rt also pulls in openNURBS, which
 * can have significant implications. */
#ifdef DM_WITH_RT
#include "rt/wdb.h"
#endif /* DM_WITH_RT */

/* Stripped down form of dm_draw_viewobjs that does just what's needed for the new setup */
DM_EXPORT extern void dm_draw_objs(void *view_ctx);

__END_DECLS

#endif /* DM_VIEW_H */

/** @} */
/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
