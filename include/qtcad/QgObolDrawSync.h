/*                  Q G O B O L D R A W S Y N C . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file qtcad/QgObolDrawSync.h */

#ifndef QGOBOLDRAWSYNC_H
#define QGOBOLDRAWSYNC_H

#include "qtcad/QgLegacyView.h"

class QgView;
struct ged;

/**
 * Synchronize one reconciled GED draw transaction into a QgView's Obol
 * controller.
 *
 * This is a transitional bridge for qtcad/qged while GED draw state is still
 * retained in BSG.  The function updates Obol database-source nodes and
 * realizes pending Obol geometry so the QgView can render/capture the migrated
 * scene.  Returns non-zero when the Obol scene changed.
 */
QTCAD_EXPORT int qg_obol_sync_draw_transaction(struct ged *gedp,
	const qg_legacy_view_draw_transaction *txn,
	const qg_legacy_view_draw_transaction_result *result,
	QgView *display);

#endif /* QGOBOLDRAWSYNC_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
