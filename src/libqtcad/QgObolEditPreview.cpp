/*            Q G O B O L E D I T P R E V I E W . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolEditPreview.cpp */

#include "common.h"

#include "qtcad/QgObolEditPreview.h"

#include "brlobol/view_controller.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgView.h"

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

    BRLObolViewController *obol = display->obolViewController();
    if (!obol)
	return 0;

    int ret = obol->replaceEditPreviewWithIntent(previewId, identity,
	    editIntentId, editIntentRole, points, commands, count,
	    sourceRevision, inputsRevision);
    if (ret < 0)
	return 0;

    display->need_update(QG_VIEW_DRAWN);
    return ret;
}

int
qg_obol_edit_preview_clear(QgView *display, const char *previewId)
{
    if (!display || !previewId || !previewId[0])
	return 0;

    BRLObolViewController *obol = display->obolViewController();
    if (!obol)
	return 0;

    int ret = obol->removeEditPreview(previewId);
    if (ret < 0)
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
