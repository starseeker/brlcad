/*              Q G O B O L V I E W S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolViewSync.cpp */

#include "common.h"

#include "QgObolViewSyncPrivate.h"

#include "qtcad/QgLegacyView.h"
#include "qtcad/QgView.h"

int
qg_obol_display_accepts_ged_active_view(struct ged *gedp, QgView *display)
{
    if (!display)
	return 0;

    qg_legacy_view *active_view = qg_legacy_view_ged_active_get(gedp);
    if (!active_view)
	return 1;

    return display->view() == active_view;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
