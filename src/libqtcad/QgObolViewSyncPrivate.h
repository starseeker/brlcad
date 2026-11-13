/*            Q G O B O L V I E W S Y N C P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolViewSyncPrivate.h
 *
 * Private qtcad/qged bridge for active-view scoped Obol sync publications.
 */

#ifndef QGOBOLVIEWSYNCPRIVATE_H
#define QGOBOLVIEWSYNCPRIVATE_H

#include "qtcad/defines.h"

class QgView;
struct ged;

QTCAD_EXPORT int qg_obol_display_accepts_ged_active_view(struct ged *gedp,
	QgView *display);

#endif /* QGOBOLVIEWSYNCPRIVATE_H */

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
