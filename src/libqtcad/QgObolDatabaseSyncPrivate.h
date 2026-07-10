/*        Q G O B O L D A T A B A S E S Y N C P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolDatabaseSyncPrivate.h
 *
 * Private qtcad bridge for synchronizing GED draw database sources into the
 * Obol view controller.  Public qtcad callers use the neutral GED draw sync
 * path in QgObolDrawSyncPrivate.h instead of this support layer directly.
 */

#ifndef QGOBOLDATABASESYNCPRIVATE_H
#define QGOBOLDATABASESYNCPRIVATE_H

#include "qtcad/defines.h"

#include <stdint.h>

class QgView;
struct db_i;

enum QgObolDatabaseDrawMode {
    QG_OBOL_DATABASE_WIREFRAME = 0,
    QG_OBOL_DATABASE_SHADED = 1
};

QTCAD_EXPORT int qg_obol_sync_database_sources(struct db_i *dbip,
	const char * const *paths,
	int path_count,
	int draw_mode,
	uint32_t source_revision,
	QgView *display);

QTCAD_EXPORT int qg_obol_remove_database_sources(const char * const *paths,
	int path_count,
	QgView *display);

QTCAD_EXPORT int qg_obol_clear_database_sources(QgView *display);

#endif /* QGOBOLDATABASESYNCPRIVATE_H */

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
