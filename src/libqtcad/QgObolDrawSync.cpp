/*              Q G O B O L D R A W S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolDrawSync.cpp */

#include "common.h"

#include "QgObolDrawSyncPrivate.h"
#include "QgObolSelectionSyncPrivate.h"

#include "BObol/BSceneController.h"
#include "BObol/BViewController.h"
#include "ged/draw.h"
#include "ged/draw_obol.h"
#include "ged/view.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgView.h"

int
qg_obol_display_accepts_ged_draw_transaction_view(
	const struct ged_draw_transaction *txn,
	QgView *display)
{
    if (!display)
	return 0;
    if (!txn || !txn->view)
	return 1;
    return display->viewContext() == txn->view;
}

int
qg_obol_ged_draw_transaction_has_view(
	const struct ged_draw_transaction *txn)
{
    return (txn && txn->view) ? 1 : 0;
}

int
qg_obol_sync_ged_draw_transaction(struct ged *gedp,
	const struct ged_draw_transaction *txn,
	const struct ged_draw_transaction_result *result,
	QgView *display)
{
    if (!gedp || !txn || (result && result->status < 0) || !display)
	return 0;

    BObolViewController *obol = display->obolViewController();
    if (!obol)
	return 0;

    struct ged_draw_transaction sync_txn = *txn;
    if (!sync_txn.view)
    sync_txn.view = display->viewContext();

    if (!ged_draw_obol_controller_attach_for_view(gedp, sync_txn.view,
	    obol, 1))
	return 0;

    if (!ged_view_context_is_independent(sync_txn.view) &&
	ged_draw_obol_controller(gedp) != obol) {
	(void)qg_obol_sync_selection_state_if_active(gedp, display, nullptr);
	display->need_update(QG_VIEW_REFRESH);
	return 1;
    }

    const int transaction_changed = result && result->status >= 0 &&
	(result->scene_revision_after != result->scene_revision_before ||
	 result->affected_groups > 0 || result->affected_shapes > 0 ||
	 result->stale_count > 0 || result->redrawn_count > 0);

    /* libged's database-source draw provider has already published a
     * successful draw to the primary endpoint before observers run. */
    if (sync_txn.kind == GED_DRAW_TXN_DRAW && result && result->status > 0 &&
	ged_draw_obol_controller(gedp) == obol) {
	(void)qg_obol_sync_selection_state_if_active(gedp, display, nullptr);
	display->need_update(QG_VIEW_REFRESH);
	return 1;
    }

    SoBRLSceneController *scene = obol->getSceneController();
    const int changed = ged_draw_obol_scene_sync_transaction(gedp, &sync_txn,
	    result, scene);
    const int selection_changed =
	qg_obol_sync_selection_state_if_active(gedp, display, nullptr);
    if (changed || transaction_changed || selection_changed)
	display->need_update(QG_VIEW_REFRESH);
    return changed || transaction_changed || selection_changed;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
