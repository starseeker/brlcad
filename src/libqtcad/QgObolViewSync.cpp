/*              Q G O B O L V I E W S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolViewSync.cpp */

#include "common.h"

#include "qtcad/QgObolViewSync.h"

#include "ged.h"
#include "ged/bsg_ged_draw.h"
#include "qtcad/QgLegacyViewBsg.h"
#include "qtcad/QgView.h"

int
qg_obol_display_accepts_ged_active_view(struct ged *gedp, QgView *display)
{
    if (!display)
	return 0;

    if (!gedp || !gedp->ged_gvp)
	return 1;

    return qg_legacy_view_to_bsg(display->view()) == gedp->ged_gvp;
}

int
qg_obol_display_accepts_draw_transaction_view(
	const struct ged_draw_transaction *txn, QgView *display)
{
    if (!display)
	return 0;

    if (!txn || !txn->view)
	return 1;

    return qg_legacy_view_to_bsg(display->view()) == txn->view;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
