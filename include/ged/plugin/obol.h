/*                 G E D / P L U G I N / O B O L . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged/plugin/obol.h
 *
 * Deliberately backend-specific GED plugin extension.
 *
 * Ordinary GED clients should use the renderer-neutral scene and view-feature
 * APIs.  Command plugins that intentionally depend on Obol/Coin may include
 * this header to reach the active endpoint and, in C++, its view controller
 * and retained scene.  Such plugins are coupled to libBObol by design.
 */

#ifndef GED_PLUGIN_OBOL_H
#define GED_PLUGIN_OBOL_H

#include "BObol/BDisplayEndpoint.h"
#include "ged/view_types.h"

__BEGIN_DECLS

/**
 * Return the borrowed active Obol endpoint for @p view_ctx, or NULL.
 *
 * The result is valid only until endpoint replacement or view teardown and
 * may be used only on the view's owner thread.  A plugin must remove callbacks
 * and scene nodes it installed before it unloads.  The endpoint's controller
 * is available with bobol_display_endpoint_controller(); C++ plugins may use
 * ged_plugin_obol_view_controller() below.
 */
GED_EXPORT extern bobol_display_endpoint_t *
ged_plugin_obol_endpoint_get(const struct ged_view_context *view_ctx);

__END_DECLS

#ifdef __cplusplus
class BObolViewController;

/**
 * Return the borrowed active Obol view controller for @p view_ctx, or NULL.
 *
 * This is the unrestricted extension point.  Prefer publishing custom nodes
 * through BObolFeatureStore::publishCustomNode() so normal feature ownership,
 * visibility, picking metadata, and teardown remain available.  Direct scene
 * mutation is permitted but becomes the plugin's responsibility and must not
 * replace or clear GED's managed CAD roots.
 */
static inline BObolViewController *
ged_plugin_obol_view_controller(const struct ged_view_context *view_ctx)
{
    bobol_display_endpoint_t *endpoint =
	ged_plugin_obol_endpoint_get(view_ctx);
    return endpoint ? static_cast<BObolViewController *>(
	bobol_display_endpoint_controller(endpoint)) : NULL;
}
#endif

#endif /* GED_PLUGIN_OBOL_H */
