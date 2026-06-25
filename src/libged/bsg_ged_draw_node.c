/*             B S G _ G E D _ D R A W _ N O D E . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file libged/bsg_ged_draw_node.c
 *
 * Remaining scene-level helpers for GED draw overlays and release.
 */

#include "common.h"

#include "bsg/draw_set.h"
#include "bsg/scene_builder.h"
#include "ged/draw.h"
#include "./ged_private.h"
#include "./bsg_ged_draw_private.h"


static void
_ged_draw_scene_ref_release_data_recurse(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return;

    for (size_t i = 0; i < ged_draw_scene_ref_child_count(ref); i++)
	_ged_draw_scene_ref_release_data_recurse(
		ged_draw_scene_ref_child_at(ref, i));
    ged_draw_scene_ref_highlight_free_cb(ref);
}


void
ged_draw_scene_ref_release(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return;
    _ged_draw_scene_ref_release_data_recurse(ref);
    bsg_scene_ref_destroy(ref);
}


int
ged_draw_scene_ref_detach(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 0;

    bsg_scene_ref parent_ref = ged_draw_scene_ref_parent(ref);
    if (ged_draw_scene_ref_is_null(parent_ref))
	return 0;

    bsg_scene_remove_child(parent_ref, ref);
    return 1;
}


int
ged_draw_scene_ref_append_child(bsg_scene_ref parent_ref, bsg_scene_ref child_ref)
{
    if (ged_draw_scene_ref_is_null(parent_ref) ||
	    ged_draw_scene_ref_is_null(child_ref))
	return 0;

    bsg_scene_append_child(parent_ref, child_ref);
    bsg_scene_bbox_invalidate(parent_ref);
    return 1;
}


bsg_scene_ref
ged_draw_view_context_group_child_ensure(void *view_ctx,
					 bsg_scene_ref parent_ref,
					 const char *name,
					 void *dp_hint)
{
    if (ged_draw_scene_ref_is_null(parent_ref) || !name)
	return bsg_scene_ref_null();

    bsg_scene_ref existing = bsg_scene_group_find_child(parent_ref, name);
    if (!ged_draw_scene_ref_is_null(existing))
	return existing;

    if (!view_ctx)
	return bsg_scene_ref_null();

    bsg_scene_ref child_ref =
	bsg_scene_group_ensure_child(parent_ref, (struct bsg_view *)view_ctx,
		name, dp_hint);
    if (ged_draw_scene_ref_is_null(child_ref))
	return bsg_scene_ref_null();

    bsg_scene_bump_rev(parent_ref);
    return child_ref;
}


bsg_scene_ref
ged_draw_view_context_group_create(void *view_ctx, const char *name)
{
    if (!view_ctx || !name)
	return bsg_scene_ref_null();

    return bsg_scene_group_create((struct bsg_view *)view_ctx, name);
}


void
ged_draw_scene_ref_erase_nested_subpath(
	bsg_scene_ref parent_ref,
	const char * const *comp_names,
	size_t comp_count,
	size_t depth_start,
	int (*match_fn)(bsg_scene_ref shape_ref, void *match_ctx),
	void *match_ctx)
{
    if (ged_draw_scene_ref_is_null(parent_ref) || !comp_names ||
	    comp_count == 0 || depth_start >= comp_count)
	return;

    bsg_scene_erase_nested_subpath(parent_ref, comp_names, comp_count,
	    depth_start, match_fn, match_ctx);
}


void
ged_draw_scene_ref_free_group(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return;

    bsg_scene_free_group(ref);
}


void
ged_draw_scene_ref_free_group_contents(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return;

    bsg_scene_free_group_contents(ref);
}


/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
