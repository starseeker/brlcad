/*                    G E D _ D R A W _ R O O T . C
 * BRL-CAD
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
/** @file libged/ged_draw_scene_root.c
 *
 * Scene-root and view-root management for GED draw state.
 */

#include "common.h"

#include <stdlib.h>
#include <string.h>

#include "bu/ptbl.h"
#include "bu/str.h"
#include "bu/color.h"
#include "bu/hash.h"
#include "bg/plot3.h"
#include "bg/clip.h"
#include "bv.h"

#include "ged.h"
#include "ged/draw.h"
#include "rt/view.h"
#include "./ged_draw_private.h"
#include "./ged_private.h"

void *
ged_draw_active_view_ctx(struct ged *gedp)
{
    return ged_view_active_ctx(gedp);
}


void
ged_draw_active_view_ctx_set(struct ged *gedp, void *view_ctx)
{
    ged_view_active_ctx_set(gedp, view_ctx);
}


/*
 * Lazily create (on first draw) and return the per-GED retained scene root.
 * Also attaches the root to each GED view handle so the render path can
 * traverse the shared scene without exposing the pointer through the public
 * view struct.
 */
static void *
_sg_root(struct ged *gedp)
{
    if (!gedp)
	return NULL;
    return ged_draw_source_root_attach_view_contexts(gedp,
	    ged_draw_active_view_ctx(gedp), ged_view_set_views_ctx(gedp));
}



void *
ged_draw_view_context_scene_root(void *view_ctx)
{
    return ged_draw_scene_handle_context(ged_view_context_scene_root_ref(view_ctx));
}


int
ged_draw_view_context_scene_attached(void *view_ctx)
{
    return ged_view_context_scene_attached(view_ctx);
}


uint64_t
ged_draw_view_context_frame_revision(void *view_ctx)
{
    return bv_frame_revision_get(
	       bv_context_view_const((const struct bv_context *)view_ctx));
}


uint64_t
ged_draw_view_context_bump_frame_revision(void *view_ctx)
{
    return bv_frame_revision_bump(bv_context_view((struct bv_context *)view_ctx));
}


static void *
ged_draw_ensure_root(struct ged *gedp)
{
    if (!gedp)
        return NULL;
    return _sg_root(gedp);
}


int
ged_draw_ensure_root_attached(struct ged *gedp)
{
    return ged_draw_ensure_root(gedp) ? 1 : 0;
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
