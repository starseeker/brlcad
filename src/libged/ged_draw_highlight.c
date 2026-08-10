/*                    G E D _ D R A W _ H I G H L I G H T . C
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
/** @file libged/ged_draw_highlight.c
 *
 * Highlighted-shape tracking and highlight-state mutation.
 */

#include "common.h"

#include "ged.h"
#include "ged/draw.h"
#include "./ged_private.h"
#include "./ged_draw_private.h"

/* ------------------------------------------------------------------ */
/* highlighted draw-ref tracker                                       */
/* ------------------------------------------------------------------ */
static ged_draw_shape_ref
_sg_highlighted_shape_ref(const struct ged_drawable *gdp)
{
    ged_draw_shape_ref ref = GED_DRAW_SHAPE_REF_NULL;
    if (!gdp || !gdp->gd_highlight_token)
	return ref;
    ref.token = gdp->gd_highlight_token;
    ref.scene_revision = gdp->gd_highlight_scene_rev;
    return ref;
}


void
ged_draw_highlighted_shape_ref_invalidate(struct ged *gedp)
{
    if (!gedp || !gedp->i || !gedp->i->ged_gdp)
	return;
    struct ged_drawable *gdp = gedp->i->ged_gdp;
    if (!gdp->gd_highlight_token && !gdp->gd_highlight_scene_rev)
	return;
    gdp->gd_highlight_token = 0;
    gdp->gd_highlight_scene_rev = 0;
    gdp->gd_highlight_rev++;
}


/*
 * Register @p ref as the single currently highlighted shape.
 * Clears the highlight state of any previously registered shape first.
 * Passing a null ref clears the previously registered shape and deregisters
 * the single-highlight ref.  Operations that set multiple highlighted records
 * use ged_draw_highlighted_shape_ref_invalidate() instead so their highlighted
 * records are not undone while the O(N) clear fallback remains safe.
 *
 * Highlight identity is a draw ref.  Node highlight bits are updated here
 * only as a derived compatibility surface for older callers that still ask
 * appearance directly.
 */
static int
_sg_set_highlighted_shape_ref(struct ged *gedp, ged_draw_shape_ref ref)
{
    if (!gedp || !gedp->i || !gedp->i->ged_gdp)
	return 0;

    struct ged_drawable *gdp = gedp->i->ged_gdp;
    ged_draw_shape_ref old_ref = _sg_highlighted_shape_ref(gdp);

    if (old_ref.token == ref.token &&
	    old_ref.scene_revision == ref.scene_revision) {
	return !ged_draw_shape_ref_is_null(ref);
    }

    (void)ged_draw_shape_ref_set_highlighted(gedp, old_ref, 0);
    int set = !ged_draw_shape_ref_is_null(ref) &&
	ged_draw_shape_ref_set_highlighted(gedp, ref, 1);

    gdp->gd_highlight_token = set ? ref.token : 0;
    gdp->gd_highlight_scene_rev = set ? ref.scene_revision : 0;

    /* Every transition is itself a highlight-state change. */
    gdp->gd_highlight_rev++;
    return set || ged_draw_shape_ref_is_null(ref);
}


void
ged_draw_set_highlighted_shape_ref(struct ged *gedp, ged_draw_shape_ref ref)
{
    if (!gedp)
	return;
    _sg_set_highlighted_shape_ref(gedp, ref);
}


int
ged_draw_shape_set_highlighted(struct ged *gedp, ged_draw_shape_ref ref, int highlighted)
{
    if (!gedp || ged_draw_shape_ref_is_null(ref))
	return 0;

    if (!highlighted) {
	if (!ged_draw_shape_ref_set_highlighted(gedp, ref, 0))
	    return 0;
	struct ged_drawable *gdp = gedp->i ? gedp->i->ged_gdp : NULL;
	if (gdp && _sg_highlighted_shape_ref(gdp).token == ref.token)
	    (void)_sg_set_highlighted_shape_ref(gedp,
		    GED_DRAW_SHAPE_REF_NULL);
    } else {
	if (!_sg_set_highlighted_shape_ref(gedp, ref))
	    return 0;
    }
    return 1;
}


/**
 * Return the highlight-state revision counter.  Bumped on every transition
 * of the highlighted draw ref.
 *
 * Cache a snapshot, then compare against a later live read to detect
 * "highlight may have changed since I last looked" without calling
 * ged_draw_highlighted_shape repeatedly.
 */
uint64_t
ged_draw_highlight_revision(const struct ged *gedp)
{
    if (!gedp || !gedp->i || !gedp->i->ged_gdp)
        return 0;
    return gedp->i->ged_gdp->gd_highlight_rev;
}


uint64_t
ged_draw_highlight_revision_advance(struct ged *gedp)
{
    if (!gedp || !gedp->i || !gedp->i->ged_gdp)
	return 0;
    return ++gedp->i->ged_gdp->gd_highlight_rev;
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
