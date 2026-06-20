/*          Q G O B O L S E L E C T I O N S Y N C . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file qtcad/QgObolSelectionSync.h */

#ifndef QGOBOLSELECTIONSYNC_H
#define QGOBOLSELECTIONSYNC_H

#include "qtcad/defines.h"

class QgView;
struct ged;

/**
 * Synchronize GED semantic selection state into a QgView's Obol realized
 * database-source geometry.  Returns non-zero when any Obol selection field
 * changed.
 */
QTCAD_EXPORT int qg_obol_sync_selection_state(struct ged *gedp,
	QgView *display,
	const char *setName);

#endif /* QGOBOLSELECTIONSYNC_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
