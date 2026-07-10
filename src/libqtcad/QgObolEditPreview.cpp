/*            Q G O B O L E D I T P R E V I E W . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolEditPreview.cpp */

#include "common.h"

#include "QgObolEditPreviewPrivate.h"

#include "brlobol/view_controller.h"
#include "ged/draw.h"
#include "QgLegacyViewContext.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgView.h"

static int
qg_obol_edit_preview_command_to_ged(int32_t command)
{
    switch (command) {
	case QG_OBOL_EDIT_PREVIEW_DRAW:
	    return GED_DRAW_VIEW_LINE_DRAW;
	case QG_OBOL_EDIT_PREVIEW_POINT:
	    return GED_DRAW_VIEW_LINE_POINT_DRAW;
	case QG_OBOL_EDIT_PREVIEW_MOVE:
	default:
	    return GED_DRAW_VIEW_LINE_MOVE;
    }
}

static ged_draw_view_feature_ref
qg_obol_edit_preview_feature(QgView *display,
			     const char *previewId,
			     const char *source_path)
{
    if (!display || !previewId || !previewId[0])
	return GED_DRAW_VIEW_FEATURE_REF_NULL;

    qg_legacy_view *view = display->view();
    void *view_ctx = qg_legacy_view_to_context(view);
    if (!view_ctx)
	return GED_DRAW_VIEW_FEATURE_REF_NULL;

    return ged_draw_view_context_feature_overlay_ensure(view_ctx,
	    previewId, display, NULL, NULL, source_path);
}

int
qg_obol_edit_preview_update(QgView *display,
	const char *previewId,
	const char *identity,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	uint32_t sourceRevision,
	uint32_t inputsRevision)
{
    return qg_obol_edit_preview_update_with_intent(display, previewId,
	    identity, NULL, NULL, points, commands, count, sourceRevision,
	    inputsRevision);
}

int
qg_obol_edit_preview_update_with_intent(QgView *display,
	const char *previewId,
	const char *identity,
	const char *editIntentId,
	const char *editIntentRole,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	uint32_t sourceRevision,
	uint32_t inputsRevision)
{
    if (!display || !previewId || !previewId[0] || !points || !commands ||
	    count <= 0)
	return 0;

    ged_draw_view_feature_ref ref =
	qg_obol_edit_preview_feature(display, previewId, identity);
    if (!ged_draw_view_feature_ref_is_null(ref)) {
	point_t *ged_points = (point_t *)bu_calloc((size_t)count,
		sizeof(point_t),
		"qtcad edit preview GED points");
	int *ged_cmds = (int *)bu_calloc((size_t)count, sizeof(int),
		"qtcad edit preview GED commands");
	for (int i = 0; i < count; i++) {
	    ged_points[i][X] = points[i][0];
	    ged_points[i][Y] = points[i][1];
	    ged_points[i][Z] = points[i][2];
	    ged_cmds[i] = qg_obol_edit_preview_command_to_ged(commands[i]);
	}

	int ret = ged_draw_view_feature_edit_preview_replace(ref,
		identity, editIntentId, editIntentRole,
		ged_points, ged_cmds, (size_t)count,
		sourceRevision, inputsRevision);
	bu_free(ged_points, "qtcad edit preview GED points");
	bu_free(ged_cmds, "qtcad edit preview GED commands");
	if (!ret)
	    return 0;

	display->need_update(QG_VIEW_DRAWN);
	return 1;
    }

    BRLObolViewController *obol = display->obolViewController();
    if (!obol)
	return 0;

    int ret = obol->replaceEditPreviewWithIntent(previewId, identity,
	    editIntentId, editIntentRole, points, commands, count,
	    sourceRevision, inputsRevision);
    if (ret < 0)
	return 0;

    display->need_update(QG_VIEW_DRAWN);
    return 1;
}

int
qg_obol_edit_preview_clear(QgView *display, const char *previewId)
{
    if (!display || !previewId || !previewId[0])
	return 0;

    qg_legacy_view *view = display->view();
    void *view_ctx = qg_legacy_view_to_context(view);
    if (view_ctx && ged_draw_view_context_feature_remove(view_ctx, previewId)) {
	display->need_update(QG_VIEW_DRAWN);
	return 1;
    }

    BRLObolViewController *obol = display->obolViewController();
    if (!obol)
	return 0;
    if (obol->removeEditPreview(previewId) < 0)
	return 0;

    display->need_update(QG_VIEW_DRAWN);
    return 1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
