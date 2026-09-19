/*            D I S P L A Y _ O B O L _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file ged/display_obol_private.h
 *
 * Uninstalled adapter between GED's renderer-neutral display handle and the
 * current Obol retained-render implementation.
 *
 * Only libged and BRL-CAD's in-tree graphical hosts may include this file.
 * Command implementations and external GED clients must use ged/display.h.
 * The conversions are representation-only and add no allocation, refcount,
 * virtual dispatch, or per-frame work.
 */

#ifndef GED_DISPLAY_OBOL_PRIVATE_H
#define GED_DISPLAY_OBOL_PRIVATE_H

#include "BObol/BDisplayEndpoint.h"
#include "ged/display.h"

__BEGIN_DECLS

struct ged_display_endpoint;
typedef struct ged_display_endpoint ged_display_endpoint_t;

GED_EXPORT ged_display_endpoint_t *
ged_view_context_display_endpoint_get(const struct ged_view_context *view);

GED_EXPORT int
ged_view_context_display_endpoint_set(struct ged_view_context *view,
	ged_display_endpoint_t *endpoint, int take_ownership);

static inline bobol_display_endpoint_t *
ged_view_context_obol_endpoint_get(const struct ged_view_context *view)
{
    return (bobol_display_endpoint_t *)
	ged_view_context_display_endpoint_get(view);
}

static inline int
ged_view_context_obol_endpoint_set(struct ged_view_context *view,
	bobol_display_endpoint_t *endpoint, int take_ownership)
{
    return ged_view_context_display_endpoint_set(
	view, (ged_display_endpoint_t *)endpoint, take_ownership);
}

/** In-tree bridge for binding an optional borrowed renderer window host. */
GED_EXPORT int
ged_view_framebuffer_backend_install(struct ged *gedp,
	struct ged_view_context *view_ctx,
	void *window_host,
	int width,
	int height,
	int present_on_flush);

__END_DECLS

#endif /* GED_DISPLAY_OBOL_PRIVATE_H */
