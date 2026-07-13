/*      Q G O B O L S E L E C T I O N S Y N C P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolSelectionSyncPrivate.h
 *
 * Private qtcad/qged bridge for GED selection state synchronization.
 */

#ifndef QGOBOLSELECTIONSYNCPRIVATE_H
#define QGOBOLSELECTIONSYNCPRIVATE_H

#include "qtcad/defines.h"

class QgView;
struct ged;

QTCAD_EXPORT int qg_obol_sync_selection_state(struct ged *gedp,
	QgView *display,
	const char *setName);

QTCAD_EXPORT int qg_obol_sync_selection_state_if_active(struct ged *gedp,
	QgView *display,
	const char *setName);

#endif /* QGOBOLSELECTIONSYNCPRIVATE_H */

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
