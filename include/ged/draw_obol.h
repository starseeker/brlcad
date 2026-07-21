/*                    D R A W _ O B O L . H
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
/** @file ged/draw_obol.h
 *
 * C integration operations for endpoint-owned framebuffer presentation.
 * Drawing-aware libged commands use libBObol's typed C++ services directly;
 * these declarations do not define or own a second draw session.
 */

#ifndef GED_DRAW_OBOL_H
#define GED_DRAW_OBOL_H

#include "common.h"

#include <stddef.h>
#include <stdint.h>

#include "ged/defines.h"
#include "vmath.h"

__BEGIN_DECLS

struct rt_db_internal;
struct bg_tess_tol;
struct bn_tol;
struct imgstream_fb;
struct bobol_display_endpoint;

typedef int (*ged_draw_obol_framebuffer_operation_t)(
	struct imgstream_fb *fb,
	void *userdata);

/**
 * Present pending Obol framebuffer stream state into the attached view scene.
 *
 * Applications should call this from their view/update thread before render
 * traversal.  fbserv producer threads only update image-stream state; they do
 * not mutate Coin/Obol scene graph fields directly.
 */
GED_EXPORT int
ged_draw_obol_framebuffer_present(struct ged *gedp);

/**
 * Install libged's Obol/imgstream fbserv backend for @p view_ctx.
 *
 * This installs the backend operation table used by fbserv packet handlers.
 * If no transport is registered yet, libged installs its default descriptor
 * transport so the backend is safe to close.  Application hosts that own
 * toolkit-specific socket/notifier integration may replace it afterward.
 *
 * @p window_host is a borrowed BObolWindowHost pointer passed as opaque
 * storage so C callers do not depend on C++ types.  Passing NULL first uses
 * the host bound to the view's display endpoint and otherwise falls back to
 * libged's generic window host.  Non-positive dimensions are derived from the
 * view context or attached display endpoint.  @p present_on_flush should only
 * be set by callers that manage presentation on the configured host's
 * scene/update thread.
 */
GED_EXPORT int
ged_draw_obol_framebuffer_backend_install_for_view(struct ged *gedp,
	struct ged_view_context *view_ctx,
	void *window_host,
	int width,
	int height,
	int present_on_flush);

/**
 * Install the Obol framebuffer backend using the active endpoint's existing
 * host, falling back to libged's generic image host when none is attached.
 * Libimgstream provides the default descriptor transport; toolkit shells use
 * this form and may replace only that transport afterward.
 */
GED_EXPORT int
ged_draw_obol_framebuffer_backend_ensure_for_view(struct ged *gedp,
	struct ged_view_context *view_ctx);

/**
 * Run an operation against the active view's Obol-owned image stream.
 *
 * Libged ensures the framebuffer bridge and serializes @p operation with
 * fbserv transport access.  When @p publish is non-zero, a successful
 * operation is flushed into the attached Obol endpoint and schedules a view
 * refresh.  The framebuffer pointer is borrowed only for the callback.
 */
GED_EXPORT int
ged_draw_obol_framebuffer_apply_for_view(struct ged *gedp,
	struct ged_view_context *view_ctx,
	ged_draw_obol_framebuffer_operation_t operation,
	void *userdata,
	int publish);

/**
 * Render an Obol-backed GED view into a caller-owned RGB/RGBA image buffer.
 *
 * Returns 1 when @p view_ctx has a display endpoint and image data was written
 * to @p image, 0 when no endpoint is attached, and -1 on error.  The returned
 * buffer is owned by the caller and should be released with bu_free().  There
 * is no alternate retained-scene readback path.
 */
GED_EXPORT int
ged_draw_obol_view_display_image(struct ged *gedp,
	struct ged_view_context *view_ctx,
	unsigned char **image,
	int flip,
	int alpha);

/**
 * Release the libged-owned Obol framebuffer bridge, if one is active.
 *
 * Call this before destroying the display endpoint hosting the framebuffer.
 */
GED_EXPORT void
ged_draw_obol_framebuffer_release(struct ged *gedp);

/**
 * Synchronize faceplate/HUD state from a GED view context into that view's
 * Obol feature store.
 *
 * This is intentionally per-view: faceplate records describe one camera/view,
 * not shared model geometry.  Missing Obol-backed endpoints are a no-op so
 * callers may invoke this unconditionally from render/update paths.
 */
GED_EXPORT int
ged_draw_obol_faceplate_sync(struct ged *gedp,
					  struct ged_view_context *view_ctx);

/**
 * Synchronize per-view lighting state (bv_lighting_state: headlight enable,
 * camera-tracking, in-scene lights) into that view's live Obol controller.
 * Missing Obol-backed controllers are a no-op (headless views keep the bv
 * state), so callers may invoke this unconditionally.
 */
GED_EXPORT int
ged_draw_obol_lighting_sync(struct ged *gedp,
					  struct ged_view_context *view_ctx);

__END_DECLS

#endif /* GED_DRAW_OBOL_H */
