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
#include "bg/plot3.h"
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


static rt_view_scene_ref
_ged_draw_scene_ref_to_rt(bsg_scene_ref ref)
{
    rt_view_scene_ref rt_ref = RT_VIEW_SCENE_REF_NULL_INIT;
    rt_ref.opaque = ged_draw_scene_ref_context(ref);
    return rt_ref;
}


rt_view_scene_ref
ged_draw_scene_ref_to_rt_view_ref(bsg_scene_ref ref)
{
    return _ged_draw_scene_ref_to_rt(ref);
}


bsg_scene_ref
ged_draw_scene_ref_from_rt_view_ref(rt_view_scene_ref ref)
{
    return ged_draw_scene_ref_from_context(ref.opaque);
}


bsg_scene_ref
ged_draw_scene_root_ref(struct ged *gedp)
{
    return ged_draw_scene_ref_from_rt_view_ref(ged_scene_root_rt_ref(gedp));
}


struct _sg_ref_iter_ctx {
    int has_record_ancestor;
    int (*cb)(bsg_scene_ref, void *);
    void *userdata;
};


struct _sg_root_shape_iter_ctx {
    int skip_overlay_groups;
    int (*cb)(bsg_scene_ref, void *);
    void *userdata;
};


static int _sg_foreach_shape(bsg_scene_ref ref,
			     int has_record_ancestor,
			     int (*cb)(bsg_scene_ref, void *),
			     void *userdata);
static int _sg_foreach_group(bsg_scene_ref ref,
			     int (*cb)(bsg_scene_ref, void *),
			     void *userdata);


static int
_sg_foreach_shape_child_cb(bsg_scene_ref child_ref, void *userdata)
{
    struct _sg_ref_iter_ctx *ctx = (struct _sg_ref_iter_ctx *)userdata;

    return ctx ? _sg_foreach_shape(child_ref, ctx->has_record_ancestor,
	    ctx->cb, ctx->userdata) : 1;
}


static int
_sg_foreach_shape(bsg_scene_ref ref,
		  int has_record_ancestor,
		  int (*cb)(bsg_scene_ref, void *),
		  void *userdata)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 1;

    struct ged_draw_shape_record_summary shape_summary;
    int has_record = ged_draw_scene_context_shape_record_summary(
	    ged_draw_scene_ref_context(ref), &shape_summary);
    if (has_record) {
	int has_path = shape_summary.fullpath &&
	    shape_summary.fullpath->fp_len > 0;
	if ((has_path || !has_record_ancestor) && !(*cb)(ref, userdata))
	    return 0;
    }

    struct _sg_ref_iter_ctx args;
    args.has_record_ancestor = has_record_ancestor || has_record;
    args.cb = cb;
    args.userdata = userdata;
    int ret = ged_draw_scene_ref_foreach_child(ref, _sg_foreach_shape_child_cb,
	    &args);
    return ret;
}


static int
_sg_root_foreach_shape_child_cb(bsg_scene_ref child_ref, void *userdata)
{
    struct _sg_root_shape_iter_ctx *ctx =
	(struct _sg_root_shape_iter_ctx *)userdata;

    if (!ctx)
	return 1;

    if (ctx->skip_overlay_groups) {
	struct ged_draw_group_record_summary group_summary;
	if (ged_draw_scene_context_group_record_summary(
		    ged_draw_scene_ref_context(child_ref), &group_summary) &&
		group_summary.is_overlay)
	    return 1;
    }

    return _sg_foreach_shape(child_ref, 0, ctx->cb, ctx->userdata);
}


static int
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

    struct _sg_root_shape_iter_ctx ctx;
    ctx.skip_overlay_groups = skip_overlay_groups;
    ctx.cb = cb;
    ctx.userdata = userdata;
    return ged_draw_scene_ref_foreach_child(root_ref,
	    _sg_root_foreach_shape_child_cb, &ctx);
}


struct _sg_root_shape_ref_iter_ctx {
    struct ged *gedp;
    ged_draw_shape_ref_index_cb cb;
    void *userdata;
};


static int
_sg_root_foreach_shape_ref_cb(bsg_scene_ref ref, void *userdata)
{
    struct _sg_root_shape_ref_iter_ctx *ctx =
	(struct _sg_root_shape_ref_iter_ctx *)userdata;
    if (!ctx || !ctx->cb)
	return 1;

    ged_draw_shape_ref shape_ref =
	ged_draw_shape_ref_from_scene_ref(ctx->gedp, ref);
    if (ged_draw_shape_ref_is_null(shape_ref))
	return 1;
    return ctx->cb(shape_ref, ctx->userdata);
}


int
ged_draw_scene_root_foreach_shape_ref(struct ged *gedp,
				      int skip_overlay_groups,
				      ged_draw_shape_ref_index_cb cb,
				      void *userdata)
{
    if (!gedp || !cb)
	return 1;

    struct _sg_root_shape_ref_iter_ctx ctx;
    ctx.gedp = gedp;
    ctx.cb = cb;
    ctx.userdata = userdata;
    return ged_draw_scene_root_foreach_shape(gedp, skip_overlay_groups,
	    _sg_root_foreach_shape_ref_cb, &ctx);
}


static int
_sg_foreach_group_child_cb(bsg_scene_ref child_ref, void *userdata)
{
    struct _sg_ref_iter_ctx *ctx = (struct _sg_ref_iter_ctx *)userdata;

    return ctx ? _sg_foreach_group(child_ref, ctx->cb, ctx->userdata) : 1;
}


static int
_sg_foreach_group(bsg_scene_ref ref,
		  int (*cb)(bsg_scene_ref, void *),
		  void *userdata)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 1;
    struct ged_draw_scene_tree_summary summary;
    if (ged_draw_scene_context_tree_summary(ged_draw_scene_ref_context(ref),
		&summary) &&
	    summary.is_group && !(*cb)(ref, userdata))
	return 0;
    struct _sg_ref_iter_ctx args;
    args.has_record_ancestor = 0;
    args.cb = cb;
    args.userdata = userdata;
    int ret = ged_draw_scene_ref_foreach_child(ref, _sg_foreach_group_child_cb,
	    &args);
    return ret;
}


static int
_sg_root_foreach_group_child_cb(bsg_scene_ref child_ref, void *userdata)
{
    struct _sg_ref_iter_ctx *ctx = (struct _sg_ref_iter_ctx *)userdata;

    return ctx ? _sg_foreach_group(child_ref, ctx->cb, ctx->userdata) : 1;
}


static int
ged_draw_scene_root_foreach_group(struct ged *gedp,
				  int (*cb)(bsg_scene_ref, void *),
				  void *userdata)
{
    if (!gedp || !cb)
	return 1;

    bsg_scene_ref root_ref = ged_draw_scene_root_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return 1;

    struct _sg_ref_iter_ctx ctx;
    ctx.has_record_ancestor = 0;
    ctx.cb = cb;
    ctx.userdata = userdata;
    return ged_draw_scene_ref_foreach_child(root_ref,
	    _sg_root_foreach_group_child_cb, &ctx);
}


struct _sg_root_group_ref_iter_ctx {
    struct ged *gedp;
    ged_draw_group_ref_index_cb cb;
    void *userdata;
};


static int
_sg_root_foreach_group_ref_cb(bsg_scene_ref ref, void *userdata)
{
    struct _sg_root_group_ref_iter_ctx *ctx =
	(struct _sg_root_group_ref_iter_ctx *)userdata;
    if (!ctx || !ctx->cb)
	return 1;

    ged_draw_group_ref group_ref =
	ged_draw_group_ref_from_scene_ref(ctx->gedp, ref);
    if (ged_draw_group_ref_is_null(group_ref))
	return 1;
    return ctx->cb(group_ref, ctx->userdata);
}


int
ged_draw_scene_root_foreach_group_ref(struct ged *gedp,
				      ged_draw_group_ref_index_cb cb,
				      void *userdata)
{
    if (!gedp || !cb)
	return 1;

    struct _sg_root_group_ref_iter_ctx ctx;
    ctx.gedp = gedp;
    ctx.cb = cb;
    ctx.userdata = userdata;
    return ged_draw_scene_root_foreach_group(gedp,
	    _sg_root_foreach_group_ref_cb, &ctx);
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

    return ged_draw_scene_context_subtree_bounds(
	    ged_draw_scene_ref_context(root_ref), min, max, pflag);
}


static int
_sg_root_has_group_cb(bsg_scene_ref UNUSED(group_ref), void *userdata)
{
    int *has_group = (int *)userdata;
    if (has_group)
	*has_group = 1;
    return 0;
}


int
ged_draw_scene_root_has_groups(struct ged *gedp)
{
    if (!gedp)
	return 0;

    int has_group = 0;
    (void)ged_draw_scene_root_foreach_group(gedp, _sg_root_has_group_cb,
	    &has_group);
    return has_group;
}


int
ged_draw_scene_root_clear_all_scope_children(struct ged *gedp)
{
    if (!gedp)
	return 0;

    int cleared = 0;
    struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
    if (views) {
	for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	    void *view_ctx = (void *)BU_PTBL_GET(views, i);
	    if (!view_ctx || !rt_view_context_is_independent(view_ctx))
		continue;
	    bsg_scene_ref scope_ref = ged_draw_scene_ref_active_scope(gedp,
		    view_ctx, 0, 0);
	    cleared += ged_draw_scene_context_clear_scope_children(
		    ged_draw_scene_ref_context(scope_ref));
	}
    }

    bsg_scene_ref root_ref = ged_draw_scene_root_ref(gedp);
    cleared += ged_draw_scene_context_clear_scope_children(
	    ged_draw_scene_ref_context(root_ref));
    return cleared;
}


int
ged_draw_scene_root_erase_path(struct ged *gedp, const char *path)
{
    if (!gedp || !path)
	return 0;

    bsg_scene_ref root_ref = ged_draw_scene_root_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return 0;

    return ged_draw_scene_ref_erase_path_at_base(gedp, path, root_ref,
	    NULL, -1, 0);
}


int
ged_draw_scene_root_erase_path_prefix(struct ged *gedp, const char *path)
{
    if (!gedp || !path)
	return 0;

    bsg_scene_ref root_ref = ged_draw_scene_root_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return 0;

    return ged_draw_scene_ref_erase_path_prefix_at_base(gedp, path, root_ref,
	    NULL, -1, 0);
}


int
ged_draw_scene_root_erase_groups_by_name(struct ged *gedp, const char *name)
{
    if (!gedp || !name)
	return 0;

    bsg_scene_ref root_ref = ged_draw_scene_root_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return 0;

    return ged_draw_scene_ref_erase_groups_by_name(root_ref, name);
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
        return ged_draw_scene_ref_null();

    bsg_scene_ref root_ref = ged_draw_view_context_group_create(view_ctx, "_draw_root");
    if (ged_draw_scene_ref_is_null(root_ref))
        return ged_draw_scene_ref_null();

    ged_draw_scene_context_set_visible(ged_draw_scene_ref_context(root_ref), 1);

    root_rt_ref = ged_draw_scene_ref_to_rt_view_ref(root_ref);
    ged_scene_root_rt_ref_set(gedp, root_rt_ref);

    ged_draw_scene_context_attach_draw_bookkeeping(gedp,
	    ged_draw_scene_ref_context(root_ref));

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
        return ged_draw_scene_ref_null();
    return _sg_root(gedp);
}


int
ged_draw_ensure_root_attached(struct ged *gedp)
{
    return !ged_draw_scene_ref_is_null(ged_draw_ensure_root(gedp));
}


bsg_scene_ref
ged_draw_scene_ref_active_scope(struct ged *gedp, void *view_ctx, int create,
				int fallback_root)
{
    if (!gedp)
	return ged_draw_scene_ref_null();

    if (view_ctx && rt_view_context_is_independent(view_ctx)) {
	bsg_scene_ref root_ref = ged_draw_scene_ref_null();
	if (create || fallback_root) {
	    root_ref = create ? ged_draw_ensure_root(gedp) :
		ged_draw_scene_ref_from_rt_view_ref(ged_scene_root_rt_ref(gedp));
	    if (create && ged_draw_scene_ref_is_null(root_ref))
		return ged_draw_scene_ref_null();
	}

	bsg_scene_ref scope_ref = ged_draw_scene_ref_from_rt_view_ref(
		rt_view_context_independent_scope_ref(view_ctx, create));
	if (!ged_draw_scene_ref_is_null(scope_ref))
	    return scope_ref;
	return fallback_root ? root_ref : ged_draw_scene_ref_null();
    }

    return create ? ged_draw_ensure_root(gedp) :
	ged_draw_scene_ref_from_rt_view_ref(ged_scene_root_rt_ref(gedp));
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
