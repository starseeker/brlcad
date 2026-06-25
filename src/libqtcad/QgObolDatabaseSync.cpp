/*            Q G O B O L D A T A B A S E S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolDatabaseSync.cpp */

#include "common.h"

#include "QgObolDatabaseSyncPrivate.h"

#include "brlobol/database_source.h"
#include "brlobol/view_controller.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgView.h"

static int
qg_obol_database_draw_mode(int draw_mode)
{
    if (draw_mode == QG_OBOL_DATABASE_SHADED)
	return SoBRLDatabaseSource::SHADED;
    return SoBRLDatabaseSource::WIREFRAME;
}

static BRLObolViewController *
qg_obol_database_controller(QgView *display)
{
    if (!display)
	return NULL;
    return display->obolViewController();
}

int
qg_obol_sync_database_sources(struct db_i *dbip,
	const char * const *paths,
	int path_count,
	int draw_mode,
	uint32_t source_revision,
	QgView *display)
{
    if (!dbip || !paths || path_count <= 0)
	return 0;

    BRLObolViewController *obol = qg_obol_database_controller(display);
    if (!obol)
	return 0;

    int changed = 0;
    const int obol_draw_mode = qg_obol_database_draw_mode(draw_mode);
    for (int i = 0; i < path_count; i++) {
	const char *path = paths[i];
	if (!path || !path[0])
	    continue;
	if (obol->replaceDatabaseSource(path, dbip, obol_draw_mode,
		source_revision) > 0)
	    changed = 1;
    }

    if (changed) {
	obol->realizePending();
	display->need_update(QG_VIEW_REFRESH);
    }
    return changed;
}

int
qg_obol_remove_database_sources(const char * const *paths,
	int path_count,
	QgView *display)
{
    if (!paths || path_count <= 0)
	return 0;

    BRLObolViewController *obol = qg_obol_database_controller(display);
    if (!obol)
	return 0;

    int changed = 0;
    for (int i = 0; i < path_count; i++) {
	const char *path = paths[i];
	if (!path || !path[0])
	    continue;
	if (obol->removeDatabaseSource(path) > 0)
	    changed = 1;
    }

    if (changed)
	display->need_update(QG_VIEW_REFRESH);
    return changed;
}

int
qg_obol_clear_database_sources(QgView *display)
{
    BRLObolViewController *obol = qg_obol_database_controller(display);
    if (!obol)
	return 0;

    int removed = obol->clearDatabaseSources();
    if (removed > 0) {
	display->need_update(QG_VIEW_REFRESH);
	return 1;
    }
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
