/*         G E D _ D R A W _ S C E N E _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged_draw_scene_private.h
 *
 * Private, lightweight handles for libged draw-scene objects.  The handles
 * keep concrete retained-render storage out of renderer-neutral public GED
 * records while allowing the draw registry and its backend adapter to share
 * identity without copying geometry.
 */

#ifndef GED_DRAW_SCENE_PRIVATE_H
#define GED_DRAW_SCENE_PRIVATE_H

#include "common.h"

typedef enum ged_draw_scene_backend {
    GED_DRAW_SCENE_BACKEND_NONE = 0,
    GED_DRAW_SCENE_BACKEND_OBOL = 1
} ged_draw_scene_backend;

typedef struct ged_draw_scene_handle {
    void *opaque;
    unsigned int backend;
} ged_draw_scene_handle;

#define GED_DRAW_SCENE_HANDLE_NULL_INIT { NULL, GED_DRAW_SCENE_BACKEND_NONE }

static inline __attribute__((unused)) ged_draw_scene_handle
ged_draw_scene_handle_null(void)
{
    ged_draw_scene_handle ref = GED_DRAW_SCENE_HANDLE_NULL_INIT;
    return ref;
}

static inline __attribute__((unused)) ged_draw_scene_handle
ged_draw_scene_handle_make(void *opaque, unsigned int backend)
{
    ged_draw_scene_handle ref = GED_DRAW_SCENE_HANDLE_NULL_INIT;
    ref.opaque = opaque;
    ref.backend = opaque ? backend : (unsigned int)GED_DRAW_SCENE_BACKEND_NONE;
    return ref;
}

static inline __attribute__((unused)) int
ged_draw_scene_handle_is_null(ged_draw_scene_handle ref)
{
    return ref.opaque ? 0 : 1;
}

static inline __attribute__((unused)) int
ged_draw_scene_handle_equal(ged_draw_scene_handle a,
			    ged_draw_scene_handle b)
{
    if (ged_draw_scene_handle_is_null(a) ||
	    ged_draw_scene_handle_is_null(b))
	return ged_draw_scene_handle_is_null(a) &&
	    ged_draw_scene_handle_is_null(b);
    if (a.backend != GED_DRAW_SCENE_BACKEND_NONE &&
	    b.backend != GED_DRAW_SCENE_BACKEND_NONE &&
	    a.backend != b.backend)
	return 0;
    return a.opaque == b.opaque;
}

static inline __attribute__((unused)) void *
ged_draw_scene_handle_context(ged_draw_scene_handle ref)
{
    return ref.opaque;
}

static inline __attribute__((unused)) unsigned int
ged_draw_scene_handle_backend(ged_draw_scene_handle ref)
{
    return ref.opaque ? ref.backend : (unsigned int)GED_DRAW_SCENE_BACKEND_NONE;
}

#endif /* GED_DRAW_SCENE_PRIVATE_H */
