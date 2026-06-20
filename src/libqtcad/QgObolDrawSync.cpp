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
#include "ged/bsg_ged_draw.h"
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
qg_obol_transaction_paths(const struct ged_draw_transaction *txn,
	const struct ged_draw_transaction_result *result)
{
    std::vector<std::string> paths;
    if (result && bu_vls_strlen(&result->names)) {
	qg_obol_append_unique_paths_from_words(paths, bu_vls_cstr(&result->names));
	if (!paths.empty())
	    return paths;
    }

    if (txn && txn->path && txn->path[0])
	qg_obol_append_unique_path(paths, txn->path);

    if (txn && txn->paths && txn->path_count > 0) {
	for (int i = 0; i < txn->path_count; i++)
	    qg_obol_append_unique_path(paths, txn->paths[i]);
    }

    return paths;
}

static uint32_t
qg_obol_source_revision(const struct ged_draw_transaction_result *result)
{
    if (!result || !result->scene_revision_after)
	return 0;

    uint64_t revision = result->scene_revision_after;
    uint32_t folded = static_cast<uint32_t>(revision ^ (revision >> 32));
    return folded ? folded : 1;
}

static int
qg_obol_draw_mode_from_bsg(int mode)
{
    if (mode == BSG_DRAW_MODE_SHADED || mode == BSG_DRAW_MODE_SHADED_BOTS)
	return QG_OBOL_DATABASE_SHADED;
    return QG_OBOL_DATABASE_WIREFRAME;
}

static int
qg_obol_transaction_draw_mode(struct ged *gedp,
	const struct ged_draw_transaction *txn)
{
    int mode = -1;
    if (txn && txn->appearance)
	mode = txn->appearance->draw_mode;
    else if (txn && txn->kind == GED_DRAW_TXN_DRAW && txn->mode >= 0)
	mode = txn->mode;
    if (mode < 0 && gedp)
	mode = ged_draw_default_mode(gedp);
    return qg_obol_draw_mode_from_bsg(mode);
}

static int
qg_obol_drawn_path_mode(struct ged *gedp, struct bsg_view *view,
	const char *path)
{
    if (ged_draw_path_state(gedp, view, path, BSG_DRAW_MODE_SHADED_BOTS) > 0 ||
	    ged_draw_path_state(gedp, view, path, BSG_DRAW_MODE_SHADED) > 0)
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
	struct bsg_view *view,
	uint32_t sourceRevision,
	QgView *display)
{
    if (!display)
	return 0;

    int changed = qg_obol_clear_database_sources(display);
    if (!gedp || !gedp->dbip)
	return changed;

    struct bu_vls listed = BU_VLS_INIT_ZERO;
    (void)ged_draw_list_paths(gedp, view, -1, 0, &listed);

    std::vector<std::string> paths;
    qg_obol_append_unique_paths_from_words(paths, bu_vls_cstr(&listed));
    bu_vls_free(&listed);

    for (const std::string &path : paths) {
	int drawMode = qg_obol_drawn_path_mode(gedp, view, path.c_str());
	const char *cpath = path.c_str();
	if (qg_obol_sync_database_sources(gedp->dbip, &cpath, 1,
		drawMode, sourceRevision, display) > 0)
	    changed = 1;
    }

    return changed;
}

int
qg_obol_sync_draw_transaction(struct ged *gedp,
	const struct ged_draw_transaction *txn,
	const struct ged_draw_transaction_result *result,
	QgView *display)
{
    if (!txn || (result && result->status < 0) || !display)
	return 0;

    BRLObolViewController *obol = display->obolViewController();
    if (!obol)
	return 0;

    struct bsg_view *view = txn->view ? txn->view : display->view();
    uint32_t sourceRevision = qg_obol_source_revision(result);
    int changed = 0;

    switch (txn->kind) {
	case GED_DRAW_TXN_DRAW:
	{
	    std::vector<std::string> paths =
		qg_obol_transaction_paths(txn, result);
	    if (paths.empty()) {
		changed = qg_obol_sync_full_scene(gedp, view,
			sourceRevision, display);
	    } else {
		changed = qg_obol_replace_paths(gedp ? gedp->dbip : NULL, paths,
			qg_obol_transaction_draw_mode(gedp, txn),
			sourceRevision, display);
	    }
	    break;
	}
	case GED_DRAW_TXN_ERASE:
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
	case GED_DRAW_TXN_CLEAR:
	case GED_DRAW_TXN_CLEAR_SCOPE:
	case GED_DRAW_TXN_TEARDOWN:
	    changed = qg_obol_clear_database_sources(display);
	    break;
	case GED_DRAW_TXN_VISIBILITY:
	case GED_DRAW_TXN_MATERIAL_CHANGED:
	case GED_DRAW_TXN_REFRESH_MATERIAL_COLORS:
	case GED_DRAW_TXN_STALE_SOURCE:
	case GED_DRAW_TXN_SOURCE_UPDATED:
	case GED_DRAW_TXN_SOURCE_RENAMED:
	case GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED:
	case GED_DRAW_TXN_ERASE_PREFIX:
	case GED_DRAW_TXN_REDRAW:
	    changed = qg_obol_sync_full_scene(gedp, view,
		    sourceRevision, display);
	    break;
	case GED_DRAW_TXN_DEFAULT_DRAW_MODE:
	case GED_DRAW_TXN_TRANSPARENCY:
	case GED_DRAW_TXN_HIGHLIGHT:
	case GED_DRAW_TXN_NONE:
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
