/*                         Q G S C E N E S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgSceneSync.cpp */

#include "common.h"

#include "ged/display_obol_private.h"

#include "QgSceneSyncPrivate.h"

#include "ged/scene.h"
#include "ged/view.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgView.h"

int
qg_scene_bind(struct ged *gedp, QgView *display)
{
    if (!gedp || !display)
	return 0;
    struct bobol_display_endpoint *endpoint = display->displayEndpoint();
    struct ged_view_context *view_ctx =
	ged_view_context_from_bv(display->viewContext());
    if (!endpoint || !view_ctx)
	return 0;
    if (ged_view_context_obol_endpoint_get(view_ctx) == endpoint)
	return 1;
    return ged_view_context_obol_endpoint_set(view_ctx, endpoint, 0);
}

int
qg_scene_delta_accepts_view(const struct ged_scene_delta *delta,
	QgView *display)
{
    if (!display)
	return 0;
    struct ged_view_context *delta_view = ged_scene_delta_view(delta);
    return !delta_view ||
	ged_view_context_from_bv(display->viewContext()) == delta_view;
}

int
qg_scene_delta_has_view(const struct ged_scene_delta *delta)
{
    return ged_scene_delta_view(delta) ? 1 : 0;
}

int
qg_scene_delta_notify(struct ged *gedp, const struct ged_scene_delta *delta,
	QgView *display)
{
    if (!gedp || !delta || !display || !qg_scene_bind(gedp, display))
	return 0;

    const int changed =
	ged_scene_delta_revision_after(delta) !=
	    ged_scene_delta_revision_before(delta) ||
	ged_scene_delta_group_count(delta) > 0 ||
	ged_scene_delta_shape_count(delta) > 0;
    if (changed) {
	/* A scene delta changes the pixels represented by the retained scene; it
	 * is not a camera refresh.  This distinction matters now that the canvas
	 * may replay an immutable completed framebuffer: DRAWN requests one new
	 * presentation, while REFRESH alone is allowed to reuse the completed
	 * frame when camera synchronization finds no transform change. */
	display->need_update(QG_VIEW_DRAWN);
    }
    return changed;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
