/*                B S G _ G E D _ D R A W _ R O O T . C
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
/** @file libged/b_s_g___g_e_d___d_r_a_w___r_o_o_t_._c.c
 *
 * Retained scene-root and view-root compatibility bridge.
 */

#include "common.h"

#include <stdlib.h>
#include <string.h>

#include "bu/ptbl.h"
#include "bu/str.h"
#include "bu/color.h"
#include "bu/hash.h"
#include "bsg/appearance.h"
#include "bsg/defines.h"
#include "bsg/database_source.h"
#include "bsg/draw_intent.h"
#include "bsg/draw_set.h"
#include "bsg/draw_source.h"
#include "bsg/field.h"
#include "bsg/geometry.h"
#include "bsg/material.h"
#include "bsg/payload.h"
#include "bg/plot3.h"
#include "bsg/scene_object.h"
#include "bsg/selection.h"
#include "bsg/view_state.h"
#include "bg/clip.h"

#include "ged.h"
#include "ged/draw.h"
#include "rt/view.h"
#include "./bsg_ged_draw_private.h"
#include "./ged_private.h"

static void
_sg_view_context_scene_root_set(void *view_ctx, rt_view_scene_ref root)
{
    (void)rt_view_context_scene_root_ref_attach(view_ctx, root);
}


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


bsg_scene_ref
ged_draw_scene_root_ref(struct ged *gedp)
{
    return ged_draw_scene_ref_from_rt_view_ref(ged_scene_root_rt_ref(gedp));
}


static int
_sg_foreach_shape(bsg_scene_ref ref,
		  int (*cb)(bsg_scene_ref, void *),
		  void *userdata)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 1;

    if (ged_draw_scene_ref_is_shape(ref)) {
	if (ged_draw_shape_state_get_scene_ref(ref))
	    return (*cb)(ref, userdata);
	return 1;
    }

    for (size_t i = 0; i < ged_draw_scene_ref_child_count(ref); i++) {
	if (!_sg_foreach_shape(ged_draw_scene_ref_child_at(ref, i), cb, userdata))
	    return 0;
    }
    return 1;
}


int
ged_draw_scene_root_foreach_shape(struct ged *gedp,
				  int skip_overlay_groups,
				  int (*cb)(bsg_scene_ref, void *),
				  void *userdata)
{
    if (!gedp || !cb)
	return 1;

    bsg_scene_ref root_ref = ged_draw_scene_root_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return 1;

    for (size_t i = 0; i < ged_draw_scene_ref_child_count(root_ref); i++) {
	bsg_scene_ref child_ref = ged_draw_scene_ref_child_at(root_ref, i);
	if (skip_overlay_groups && ged_draw_group_scene_ref_is_overlay(child_ref))
	    continue;
	if (!_sg_foreach_shape(child_ref, cb, userdata))
	    return 0;
    }
    return 1;
}


static int
_sg_foreach_group(bsg_scene_ref ref,
		  int (*cb)(bsg_scene_ref, void *),
		  void *userdata)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 1;
    if (ged_draw_scene_ref_is_group(ref) && !(*cb)(ref, userdata))
	return 0;
    for (size_t i = 0; i < ged_draw_scene_ref_child_count(ref); i++) {
	if (!_sg_foreach_group(ged_draw_scene_ref_child_at(ref, i), cb, userdata))
	    return 0;
    }
    return 1;
}


int
ged_draw_scene_root_foreach_group(struct ged *gedp,
				  int (*cb)(bsg_scene_ref, void *),
				  void *userdata)
{
    if (!gedp || !cb)
	return 1;

    bsg_scene_ref root_ref = ged_draw_scene_root_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return 1;

    for (size_t i = 0; i < ged_draw_scene_ref_child_count(root_ref); i++) {
	if (!_sg_foreach_group(ged_draw_scene_ref_child_at(root_ref, i), cb, userdata))
	    return 0;
    }
    return 1;
}


int
ged_draw_scene_root_subtree_bounds(struct ged *gedp,
				   vect_t *min,
				   vect_t *max,
				   int pflag)
{
    if (!gedp)
	return 1;

    bsg_scene_ref root_ref = ged_draw_scene_root_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return 1;

    return ged_draw_scene_ref_subtree_bounds(root_ref, min, max, pflag);
}


bsg_scene_ref
ged_draw_scene_root_last_group_ref(struct ged *gedp)
{
    if (!gedp)
	return bsg_scene_ref_null();

    bsg_scene_ref root_ref = ged_draw_scene_root_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return bsg_scene_ref_null();

    size_t child_count = ged_draw_scene_ref_child_count(root_ref);
    if (child_count == 0)
	return bsg_scene_ref_null();

    return ged_draw_scene_ref_child_at(root_ref, child_count - 1);
}


int
ged_draw_scene_root_has_groups(struct ged *gedp)
{
    if (!gedp)
	return 0;

    bsg_scene_ref root_ref = ged_draw_scene_root_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return 0;

    return ged_draw_scene_ref_child_count(root_ref) > 0 ? 1 : 0;
}


int
ged_draw_scene_root_mutate(struct ged *gedp,
			   int (*cb)(bsg_scene_ref, void *),
			   void *userdata)
{
    if (!gedp || !cb)
	return 0;

    bsg_scene_ref root_ref = ged_draw_scene_root_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return 0;

    return (*cb)(root_ref, userdata);
}


/*
 * Lazily create (on first draw) and return the per-GED retained scene root.
 * Also attaches the root to each GED view handle so the render path can
 * traverse the shared scene without exposing the pointer through the public
 * view struct.
 */
static bsg_scene_ref
_sg_root(struct ged *gedp)
{
    struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
    rt_view_scene_ref root_rt_ref = ged_scene_root_rt_ref(gedp);
    if (!rt_view_scene_ref_is_null(root_rt_ref)) {
	void *active_view_ctx = ged_draw_active_view_ctx(gedp);
	if (active_view_ctx)
	    _sg_view_context_scene_root_set(active_view_ctx, root_rt_ref);
	if (views) {
	    for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
		void *view_ctx = BU_PTBL_GET(views, i);
		if (!view_ctx)
		    continue;
		_sg_view_context_scene_root_set(view_ctx, root_rt_ref);
	    }
	}
	return ged_draw_scene_ref_from_rt_view_ref(root_rt_ref);
    }

    void *view_ctx = ged_draw_active_view_ctx(gedp);
    if (!view_ctx)
        return bsg_scene_ref_null();

    bsg_scene_ref root_ref = ged_draw_view_context_group_create(view_ctx, "_draw_root");
    if (ged_draw_scene_ref_is_null(root_ref))
        return bsg_scene_ref_null();

    ged_draw_scene_ref_set_visible(root_ref, 1);

    root_rt_ref = ged_draw_scene_ref_to_rt_view_ref(root_ref);
    ged_scene_root_rt_ref_set(gedp, root_rt_ref);

    /* Store draw-tree bookkeeping on the root so freeing helpers can bump
     * gd_draw_rev and recycle nodes without carrying a ged pointer. */
    gedp->i->ged_gdp->bsg_ctx.draw_rev = &gedp->i->ged_gdp->gd_draw_rev;
    gedp->i->ged_gdp->bsg_ctx.fso =
	rt_view_set_context_recycle_pool(ged_view_set_ctx(gedp));
    gedp->i->ged_gdp->bsg_ctx.owner_data = gedp->i->ged_gdp;
    ged_draw_scene_ref_set_draw_context(root_ref, &gedp->i->ged_gdp->bsg_ctx);

    /* Register in all views so the BSG render loop can traverse the shared
     * retained scene directly without reading gv_objs.  No per-frame
     * bsg_scene_root_sync rebuild is needed. */
    _sg_view_context_scene_root_set(view_ctx, root_rt_ref);
    if (views) {
	for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	    void *view_iter_ctx = BU_PTBL_GET(views, i);
	    if (!view_iter_ctx)
		continue;
	    _sg_view_context_scene_root_set(view_iter_ctx, root_rt_ref);
	}
    }

    return root_ref;
}



rt_view_scene_ref
ged_draw_view_context_scene_root_rt_ref(void *view_ctx)
{
    return rt_view_context_scene_root_ref(view_ctx);
}


void *
ged_draw_view_context_scene_root(void *view_ctx)
{
    return rt_view_scene_ref_context(ged_draw_view_context_scene_root_rt_ref(view_ctx));
}


int
ged_draw_view_context_scene_attached(void *view_ctx)
{
    return rt_view_context_scene_attached(view_ctx);
}


uint64_t
ged_draw_view_context_frame_revision(void *view_ctx)
{
    return rt_view_context_frame_revision_get(view_ctx);
}


uint64_t
ged_draw_view_context_bump_frame_revision(void *view_ctx)
{
    return rt_view_context_frame_revision_bump(view_ctx);
}


bsg_scene_ref
ged_draw_ensure_root(struct ged *gedp)
{
    if (!gedp)
        return bsg_scene_ref_null();
    return _sg_root(gedp);
}


int
ged_draw_ensure_root_attached(struct ged *gedp)
{
    return !ged_draw_scene_ref_is_null(ged_draw_ensure_root(gedp));
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
