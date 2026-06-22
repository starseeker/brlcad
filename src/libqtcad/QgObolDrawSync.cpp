/*              Q G O B O L D R A W S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolDrawSync.cpp */

#include "common.h"

#include "qtcad/QgObolDrawSync.h"

#include "brlobol/view_controller.h"
#include "bu/vls.h"
#include "qtcad/QgLegacyView.h"
#include "qtcad/QgObolDatabaseSync.h"
#include "qtcad/QgView.h"

#include <algorithm>
#include <stdint.h>
#include <string>
#include <vector>

static void
qg_obol_append_unique_path(std::vector<std::string> &paths, const char *path)
{
    if (!path || !path[0])
	return;
    std::string spath(path);
    if (std::find(paths.begin(), paths.end(), spath) == paths.end())
	paths.push_back(spath);
}

static void
qg_obol_append_unique_paths_from_words(std::vector<std::string> &paths,
	const char *words)
{
    if (!words || !words[0])
	return;

    std::string names(words);
    size_t pos = 0;
    while (pos < names.size()) {
	pos = names.find_first_not_of(" \t\r\n", pos);
	if (pos == std::string::npos)
	    break;
	size_t end = names.find_first_of(" \t\r\n", pos);
	std::string path = names.substr(pos,
		(end == std::string::npos) ? std::string::npos : end - pos);
	qg_obol_append_unique_path(paths, path.c_str());
	if (end == std::string::npos)
	    break;
	pos = end + 1;
    }
}

static std::vector<std::string>
qg_obol_transaction_paths(const qg_legacy_view_draw_transaction *txn,
	const qg_legacy_view_draw_transaction_result *result)
{
    std::vector<std::string> paths;
    const char *result_names = qg_legacy_view_draw_result_names(result);
    if (result_names) {
	qg_obol_append_unique_paths_from_words(paths, result_names);
	if (!paths.empty())
	    return paths;
    }

    const char *path = qg_legacy_view_draw_transaction_path(txn);
    if (path)
	qg_obol_append_unique_path(paths, path);

    int path_count = qg_legacy_view_draw_transaction_path_count(txn);
    for (int i = 0; i < path_count; i++)
	qg_obol_append_unique_path(paths,
		qg_legacy_view_draw_transaction_path_at(txn, i));

    return paths;
}

static uint32_t
qg_obol_source_revision(const qg_legacy_view_draw_transaction_result *result)
{
    uint64_t revision =
	qg_legacy_view_draw_result_scene_revision_after(result);
    if (!revision)
	return 0;

    uint32_t folded = static_cast<uint32_t>(revision ^ (revision >> 32));
    return folded ? folded : 1;
}

static int
qg_obol_draw_mode_from_legacy(int mode)
{
    if (mode == QG_LEGACY_VIEW_DRAW_MODE_SHADED ||
	    mode == QG_LEGACY_VIEW_DRAW_MODE_SHADED_BOTS)
	return QG_OBOL_DATABASE_SHADED;
    return QG_OBOL_DATABASE_WIREFRAME;
}

static int
qg_obol_transaction_draw_mode(struct ged *gedp,
	const qg_legacy_view_draw_transaction *txn)
{
    return qg_obol_draw_mode_from_legacy(
	    qg_legacy_view_draw_transaction_mode(gedp, txn));
}

static int
qg_obol_drawn_path_mode(struct ged *gedp, qg_legacy_view *view,
	const char *path)
{
    if (qg_legacy_view_draw_path_state(gedp, view, path,
	    QG_LEGACY_VIEW_DRAW_MODE_SHADED_BOTS) > 0 ||
	    qg_legacy_view_draw_path_state(gedp, view, path,
		QG_LEGACY_VIEW_DRAW_MODE_SHADED) > 0)
	return QG_OBOL_DATABASE_SHADED;
    return QG_OBOL_DATABASE_WIREFRAME;
}

static int
qg_obol_replace_paths(struct db_i *dbip,
	const std::vector<std::string> &paths,
	int drawMode,
	uint32_t sourceRevision,
	QgView *display)
{
    if (!dbip || paths.empty())
	return 0;

    std::vector<const char *> cpaths;
    cpaths.reserve(paths.size());
    for (const std::string &path : paths)
	cpaths.push_back(path.c_str());

    return qg_obol_sync_database_sources(dbip, cpaths.data(),
	static_cast<int>(cpaths.size()), drawMode, sourceRevision, display);
}

static int
qg_obol_remove_paths(const std::vector<std::string> &paths,
	QgView *display)
{
    if (paths.empty())
	return 0;

    std::vector<const char *> cpaths;
    cpaths.reserve(paths.size());
    for (const std::string &path : paths)
	cpaths.push_back(path.c_str());

    return qg_obol_remove_database_sources(cpaths.data(),
	static_cast<int>(cpaths.size()), display);
}

static int
qg_obol_sync_full_scene(struct ged *gedp,
	qg_legacy_view *view,
	uint32_t sourceRevision,
	QgView *display)
{
    if (!display)
	return 0;

    int changed = qg_obol_clear_database_sources(display);
    struct db_i *dbip = qg_legacy_view_ged_database(gedp);
    if (!dbip)
	return changed;

    struct bu_vls listed = BU_VLS_INIT_ZERO;
    (void)qg_legacy_view_draw_list_paths(gedp, view, -1, 0, &listed);

    std::vector<std::string> paths;
    qg_obol_append_unique_paths_from_words(paths, bu_vls_cstr(&listed));
    bu_vls_free(&listed);

    for (const std::string &path : paths) {
	int drawMode = qg_obol_drawn_path_mode(gedp, view, path.c_str());
	const char *cpath = path.c_str();
	if (qg_obol_sync_database_sources(dbip, &cpath, 1,
		drawMode, sourceRevision, display) > 0)
	    changed = 1;
    }

    return changed;
}

int
qg_obol_sync_draw_transaction(struct ged *gedp,
	const qg_legacy_view_draw_transaction *txn,
	const qg_legacy_view_draw_transaction_result *result,
	QgView *display)
{
    if (!txn || qg_legacy_view_draw_result_status(result) < 0 || !display)
	return 0;

    BRLObolViewController *obol = display->obolViewController();
    if (!obol)
	return 0;

    qg_legacy_view *view = qg_legacy_view_draw_transaction_view_get(txn);
    if (!view)
	view = display->view();
    uint32_t sourceRevision = qg_obol_source_revision(result);
    int changed = 0;

    switch (qg_legacy_view_draw_transaction_kind_get(txn)) {
	case QG_LEGACY_VIEW_DRAW_TXN_DRAW:
	{
	    std::vector<std::string> paths =
		qg_obol_transaction_paths(txn, result);
	    if (paths.empty()) {
		changed = qg_obol_sync_full_scene(gedp, view,
			sourceRevision, display);
	    } else {
		changed = qg_obol_replace_paths(
			qg_legacy_view_ged_database(gedp), paths,
			qg_obol_transaction_draw_mode(gedp, txn),
			sourceRevision, display);
	    }
	    break;
	}
	case QG_LEGACY_VIEW_DRAW_TXN_ERASE:
	{
	    std::vector<std::string> paths =
		qg_obol_transaction_paths(txn, result);
	    if (paths.empty())
		changed = qg_obol_sync_full_scene(gedp, view,
			sourceRevision, display);
	    else
		changed = qg_obol_remove_paths(paths, display);
	    break;
	}
	case QG_LEGACY_VIEW_DRAW_TXN_CLEAR:
	case QG_LEGACY_VIEW_DRAW_TXN_CLEAR_SCOPE:
	case QG_LEGACY_VIEW_DRAW_TXN_TEARDOWN:
	    changed = qg_obol_clear_database_sources(display);
	    break;
	case QG_LEGACY_VIEW_DRAW_TXN_VISIBILITY:
	case QG_LEGACY_VIEW_DRAW_TXN_MATERIAL_CHANGED:
	case QG_LEGACY_VIEW_DRAW_TXN_REFRESH_MATERIAL_COLORS:
	case QG_LEGACY_VIEW_DRAW_TXN_STALE_SOURCE:
	case QG_LEGACY_VIEW_DRAW_TXN_SOURCE_UPDATED:
	case QG_LEGACY_VIEW_DRAW_TXN_SOURCE_RENAMED:
	case QG_LEGACY_VIEW_DRAW_TXN_SOURCE_REFERENCES_REMOVED:
	case QG_LEGACY_VIEW_DRAW_TXN_ERASE_PREFIX:
	case QG_LEGACY_VIEW_DRAW_TXN_REDRAW:
	    changed = qg_obol_sync_full_scene(gedp, view,
		    sourceRevision, display);
	    break;
	case QG_LEGACY_VIEW_DRAW_TXN_DEFAULT_DRAW_MODE:
	case QG_LEGACY_VIEW_DRAW_TXN_TRANSPARENCY:
	case QG_LEGACY_VIEW_DRAW_TXN_HIGHLIGHT:
	case QG_LEGACY_VIEW_DRAW_TXN_NONE:
	default:
	    break;
    }

    return changed;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
