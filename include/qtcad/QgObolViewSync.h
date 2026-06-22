/*                Q G O B O L V I E W S Y N C . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file qtcad/QgObolViewSync.h */

#ifndef QGOBOLVIEWSYNC_H
#define QGOBOLVIEWSYNC_H

#include "qtcad/QgLegacyView.h"

class QgView;
struct ged;

/**
 * Return non-zero when @p display should accept a sync publication scoped to
 * @p gedp's active view.  A null GED pointer or null GED active view is treated
 * as unscoped and accepted, matching the transitional callers.
 */
QTCAD_EXPORT int qg_obol_display_accepts_ged_active_view(struct ged *gedp,
	QgView *display);

/**
 * Return non-zero when @p display should accept a sync publication scoped to
 * @p txn's view.  A null transaction pointer or null transaction view is
 * treated as unscoped and accepted.
 */
QTCAD_EXPORT int qg_obol_display_accepts_draw_transaction_view(
	const qg_legacy_view_draw_transaction *txn, QgView *display);

#endif /* QGOBOLVIEWSYNC_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
