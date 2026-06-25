/*           Q G O B O L D R A W S Y N C P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolDrawSyncPrivate.h
 *
 * Private qtcad/qged bridge for Obol draw synchronization.  Draw-sync callers
 * use neutral GED transactions; qg legacy draw transaction sync is retired
 * with the former installed public draw-sync header.
 */

#ifndef QGOBOLDRAWSYNCPRIVATE_H
#define QGOBOLDRAWSYNCPRIVATE_H

#include "qtcad/defines.h"

class QgView;
struct ged;
struct ged_draw_transaction;
struct ged_draw_transaction_result;

QTCAD_EXPORT int qg_obol_display_accepts_ged_draw_transaction_view(
	const struct ged_draw_transaction *txn,
	QgView *display);

QTCAD_EXPORT int qg_obol_ged_draw_transaction_has_view(
	const struct ged_draw_transaction *txn);

QTCAD_EXPORT int qg_obol_sync_ged_draw_transaction(
	struct ged *gedp,
	const struct ged_draw_transaction *txn,
	const struct ged_draw_transaction_result *result,
	QgView *display);

#endif /* QGOBOLDRAWSYNCPRIVATE_H */

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
