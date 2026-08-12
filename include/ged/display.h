/*                     D I S P L A Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file ged/display.h
 *
 * Renderer-neutral presentation operations for GED views.
 *
 * These functions express host-visible semantics only.  The active retained
 * renderer, framebuffer bridge, and scene implementation are private to GED.
 */

#ifndef GED_DISPLAY_H
#define GED_DISPLAY_H


#include <stddef.h>

#include "bv/display.h"
#include "ged/defines.h"

__BEGIN_DECLS

struct imgstream_fb;

/** Ensure the hosted view has a GED-owned display attachment. */
GED_EXPORT int
ged_view_context_display_endpoint_ensure(struct ged_view_context *view);

/** Copy a renderer-neutral presentation property from the view policy. */
GED_EXPORT int
ged_view_context_display_property_get(const struct ged_view_context *view,
	const char *name, struct bv_display_property_value *value);

/**
 * Set a renderer-neutral presentation property.
 *
 * GED retains the policy while an endpoint is absent, so headless clients may
 * configure a view before it is presented.
 */
GED_EXPORT int
ged_view_context_display_property_set(struct ged_view_context *view,
	const char *name, const struct bv_display_property_value *value);

typedef int (*ged_view_framebuffer_operation_t)(
	struct imgstream_fb *fb,
	void *userdata);

/** Present pending framebuffer-stream state on the active view owner thread. */
GED_EXPORT int
ged_view_framebuffer_present(struct ged *gedp);

/** Ensure a framebuffer backend exists for one view. */
GED_EXPORT int
ged_view_framebuffer_backend_ensure(struct ged *gedp,
	struct ged_view_context *view_ctx);

/** Run a serialized operation against the active view's image stream. */
GED_EXPORT int
ged_view_framebuffer_apply(struct ged *gedp,
	struct ged_view_context *view_ctx,
	ged_view_framebuffer_operation_t operation,
	void *userdata,
	int publish);

/**
 * Capture the presented view into caller-owned RGB/RGBA storage.
 *
 * Returns 1 on capture, 0 when no graphical display is attached, and -1 on
 * error.  The caller releases the returned allocation with bu_free().
 */
GED_EXPORT int
ged_view_display_image_capture(struct ged *gedp,
	struct ged_view_context *view_ctx,
	unsigned char **image,
	int flip,
	int alpha);

/** Release the GED-owned framebuffer bridge, if active. */
GED_EXPORT void
ged_view_framebuffer_release(struct ged *gedp);

/** Synchronize semantic faceplate/HUD state into the active presentation. */
GED_EXPORT int
ged_view_faceplate_sync(struct ged *gedp,
	struct ged_view_context *view_ctx);

/** Synchronize semantic lighting policy into the active presentation. */
GED_EXPORT int
ged_view_lighting_sync(struct ged *gedp,
	struct ged_view_context *view_ctx);

/** Synchronize semantic normal/shading policy into the active presentation. */
GED_EXPORT int
ged_view_shading_sync(struct ged *gedp,
	struct ged_view_context *view_ctx);

__END_DECLS

#endif /* GED_DISPLAY_H */
