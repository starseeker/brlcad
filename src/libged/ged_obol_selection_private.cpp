/*        G E D _ O B O L _ S E L E C T I O N _ P R I V A T E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged_obol_selection_private.cpp
 *
 * Obol-backed view selection adapter.
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BADC.h"
#include "BObol/BDrawCache.h"
#include "BObol/BGrid.h"
#include "BObol/BInit.h"
#include "BObol/BImageSource.h"
#include "BObol/BDisplayEndpoint.h"
#include "BObol/BLodRealization.h"
#include "BObol/BLodService.h"
#include "BObol/BMeshShape.h"
#include "BObol/BSceneController.h"
#include "BObol/BSceneGroup.h"
#include "BObol/BSnapAction.h"
#include "BObol/BSourceRealization.h"
#include "BObol/BViewportImage.h"
#include "BObol/BVListShape.h"
#include "BObol/BViewAttachment.h"
#include "BObol/BViewController.h"
#include "BObol/BViewQuery.h"
#include "BObol/BViewStore.h"
#include "bg/line_layer.h"
#include "bg/plane.h"
#include "bg/polygon.h"
#include "bv.h"
#include "bu/hash.h"
#include "bu/malloc.h"
#include "bu/parallel.h"
#include "bu/path.h"
#include "bu/str.h"
#include "bu/datetime.h"
#include "bu/units.h"
#include "bu/vls.h"
#include "ged.h"
#include "ged/db_index.h"
#include "ged/draw.h"
#include "ged/display.h"
#include "ged/view.h"
#include "icv.h"
#include "rt/db5.h"
#include "rt/db_fullpath.h"
#include "rt/search.h"
#include "rt/view.h"
#include "vmath.h"

#include "./ged_bobol_private.hpp"
#include "./draw_obol_bridge_private.hpp"
#include "./ged_draw_view_private.h"
#include "./ged_private.h"

#include <algorithm>
#include <Inventor/SbVec2f.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Obol/cad/SoCADAssembly.h>
#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <atomic>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

static const char *
ged_obol_skip_leading_slash(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}


static BObolViewController *
ged_obol_view_controller_for_context(struct ged_view_context *view_ctx)
{
    return ged_bobol_view_controller(view_ctx);
}


static BObolViewController *
ged_obol_view_controller_ensure_for_context(struct ged_view_context *view_ctx,
	int sync_current_scene)
{
    struct ged *gedp = view_ctx ?
	static_cast<struct ged *>(ged_view_context_user_data_get(view_ctx)) : NULL;
    if (!gedp)
	return NULL;
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (controller)
	return controller;
    if (!ged_draw_obol_render_endpoint_ensure_for_view(gedp, view_ctx,
	    sync_current_scene))
	return NULL;
    return ged_bobol_view_controller(view_ctx);
}


static BObolFeatureOwner
ged_obol_feature_owner(struct ged_view_context *view_ctx, int local)
{
    BObolFeatureOwner owner;
    if (!local)
	return owner;
    owner.ownerToken = view_ctx;
    owner.ownerId = ged_obol_view_scope_name(view_ctx).c_str();
    owner.ownerRole = "view";
    return owner;
}


static int
ged_obol_selection_kind_from_ged(int kind)
{
    switch (kind) {
	case -1:
	    return BOBOL_SELECTION_ALL;
	case 0:
	    return BOBOL_SELECTION_SELECTED_PATH;
	case 1:
	    return BOBOL_SELECTION_HIGHLIGHTED_REF;
	default:
	    return INT_MIN;
    }
}

static int
ged_draw_obol_selection_available_impl(struct ged_view_context *view_ctx)
{
    return ged_obol_view_controller_ensure_for_context(view_ctx, 1) ? 1 : 0;
}

static size_t
ged_draw_obol_selection_count_impl(struct ged_view_context *view_ctx)
{
    BObolViewController *controller =
	ged_obol_view_controller_ensure_for_context(view_ctx, 1);
    if (!controller)
	return 0;
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    return controller->selection().count(&owner, BOBOL_SELECTION_ALL);
}

static int
ged_draw_obol_selection_clear_impl(struct ged_view_context *view_ctx)
{
    BObolViewController *controller =
	ged_obol_view_controller_ensure_for_context(view_ctx, 1);
    if (!controller)
	return 0;
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    controller->selection().clear(&owner, BOBOL_SELECTION_ALL);
    return 1;
}

extern "C" int
ged_draw_obol_view_context_selection_available(struct ged_view_context *view_ctx)
{
    return ged_draw_obol_selection_available_impl(view_ctx);
}

extern "C" size_t
ged_draw_obol_view_context_selection_count(struct ged_view_context *view_ctx)
{
    return ged_draw_obol_selection_count_impl(view_ctx);
}

struct ged_obol_selection_foreach_context {
    struct ged_view_context *view_ctx;
    ged_view_selection_visit_cb cb;
    void *data;
    int ok;
};

static int
ged_obol_selection_foreach_callback(const SbString &path, void *data)
{
    struct ged_obol_selection_foreach_context *ctx =
	    static_cast<struct ged_obol_selection_foreach_context *>(data);
    const char *path_str = path.getString();
    if (!ctx || !ctx->cb)
	return 0;
    if (!path_str || !path_str[0])
	return 1;
    if (!ctx->cb(ctx->view_ctx, path_str, ctx->data)) {
	ctx->ok = 0;
	return 0;
    }
    return 1;
}

extern "C" int
ged_draw_obol_view_context_selection_path_foreach(
    struct ged_view_context *view_ctx,
    ged_view_selection_visit_cb cb,
    void *data)
{
    if (!cb)
	return 0;
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return 0;

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    struct ged_obol_selection_foreach_context ctx;
    ctx.view_ctx = view_ctx;
    ctx.cb = cb;
    ctx.data = data;
    ctx.ok = 1;
    controller->selection().visitPaths(ged_obol_selection_foreach_callback,
				       &ctx, &owner, BOBOL_SELECTION_ALL);
    return ctx.ok;
}

extern "C" int
ged_draw_obol_view_context_selection_clear(struct ged_view_context *view_ctx)
{
    return ged_draw_obol_selection_clear_impl(view_ctx);
}

extern "C" int
ged_draw_obol_view_context_selection_contains_path(
    struct ged_view_context *view_ctx,
    int kind,
    const char *path)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !path || !path[0])
	return 0;

    int obol_kind = ged_obol_selection_kind_from_ged(kind);
    if (obol_kind == INT_MIN)
	return 0;

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    return controller->selection().containsPath(
	       ged_obol_skip_leading_slash(path), obol_kind, &owner) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_selection_add_path(
    struct ged_view_context *view_ctx,
    int kind,
    const char *path)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !path || !path[0])
	return 0;

    int obol_kind = ged_obol_selection_kind_from_ged(kind);
    if (obol_kind == INT_MIN || obol_kind == BOBOL_SELECTION_ALL)
	return 0;

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    return controller->selection().addPath(ged_obol_skip_leading_slash(path),
					   obol_kind, &owner) ? 1 : 0;
}


extern "C" int
ged_draw_obol_view_context_selection_remove_path(
    struct ged_view_context *view_ctx,
    int kind,
    const char *path)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !path || !path[0])
	return 0;

    int obol_kind = ged_obol_selection_kind_from_ged(kind);
    if (obol_kind == INT_MIN || obol_kind == BOBOL_SELECTION_ALL)
	return 0;

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    return controller->selection().removePath(
	ged_obol_skip_leading_slash(path), obol_kind, &owner) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_selection_set_path(
    struct ged_view_context *view_ctx,
    int kind,
    const char *path)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return 0;

    int obol_kind = ged_obol_selection_kind_from_ged(kind);
    if (obol_kind == INT_MIN)
	return 0;

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    if (obol_kind == BOBOL_SELECTION_ALL) {
	controller->selection().clear(&owner, BOBOL_SELECTION_ALL);
	return path && path[0] ? 0 : 1;
    }

    return controller->selection().setPath(
	       ged_obol_skip_leading_slash(path), obol_kind, &owner) ? 1 : 0;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
