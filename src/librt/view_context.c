/*                     V I E W _ C O N T E X T . C
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
/** @file librt/view_context.c
 *
 * Neutral RT view-context helpers independent of any concrete view backend.
 */

#include "common.h"

#include "rt/view.h"
#include "./view_context_private.h"


int
rt_view_context_is_valid(const void *ctx)
{
    return (_rt_view_context_native_is(ctx) || rt_view_context_is_bsg(ctx)) ?
	   1 : 0;
}


int
rt_view_context_is_retained(const void *ctx)
{
    return rt_view_context_is_bsg(ctx);
}


rt_view_scene_ref
rt_view_scene_ref_null(void)
{
    rt_view_scene_ref ref = RT_VIEW_SCENE_REF_NULL_INIT;
    return ref;
}


rt_view_scene_ref
rt_view_scene_ref_make(void *opaque, unsigned int backend)
{
    rt_view_scene_ref ref = RT_VIEW_SCENE_REF_NULL_INIT;
    ref.opaque = opaque;
    ref.backend = opaque ? backend : RT_VIEW_SCENE_BACKEND_NONE;
    return ref;
}


int
rt_view_scene_ref_is_null(rt_view_scene_ref ref)
{
    return ref.opaque ? 0 : 1;
}


int
rt_view_scene_ref_equal(rt_view_scene_ref a, rt_view_scene_ref b)
{
    if (rt_view_scene_ref_is_null(a) || rt_view_scene_ref_is_null(b))
	return rt_view_scene_ref_is_null(a) && rt_view_scene_ref_is_null(b);
    if (a.backend != RT_VIEW_SCENE_BACKEND_NONE &&
	b.backend != RT_VIEW_SCENE_BACKEND_NONE &&
	a.backend != b.backend)
	return 0;
    return a.opaque == b.opaque;
}


void *
rt_view_scene_ref_context(rt_view_scene_ref ref)
{
    return ref.opaque;
}


unsigned int
rt_view_scene_ref_backend(rt_view_scene_ref ref)
{
    return ref.opaque ? ref.backend : RT_VIEW_SCENE_BACKEND_NONE;
}
