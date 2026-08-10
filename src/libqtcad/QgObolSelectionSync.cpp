/*        Q G O B O L S E L E C T I O N S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolSelectionSync.cpp */

#include "common.h"

#include "QgObolSelectionSyncPrivate.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BMeshShape.h"
#include "BObol/BViewController.h"
#include "BObol/BVListShape.h"
#include "bu/vls.h"
#include "ged/selection_state.h"
#include "QgObolViewSyncPrivate.h"
#include "qtcad/QgView.h"

#include <Inventor/nodes/SoGroup.h>

#include <string>
#include <unordered_set>
#include <vector>

static const char *
qg_obol_skip_leading_slash(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

static std::string
qg_obol_normalized_path(const char *path)
{
    return std::string(qg_obol_skip_leading_slash(path));
}

static bool
qg_obol_path_equal(const char *a, const char *b)
{
    const std::string na = qg_obol_normalized_path(a);
    const std::string nb = qg_obol_normalized_path(b);
    return !na.empty() && !nb.empty() && na == nb;
}

static bool
qg_obol_path_is_or_contains(const char *candidate, const char *selected)
{
    const std::string cpath = qg_obol_normalized_path(candidate);
    const std::string spath = qg_obol_normalized_path(selected);
    if (cpath.empty() || spath.empty())
	return false;
    if (cpath == spath)
	return true;
    if (cpath.size() <= spath.size())
	return false;
    return cpath.compare(0, spath.size(), spath) == 0 &&
	cpath[spath.size()] == '/';
}

static std::vector<std::string>
qg_obol_selection_paths(struct ged *gedp, const char *setName)
{
    std::vector<std::string> ret;
    if (!ged_selection_state_available(gedp))
	return ret;

    struct bu_vls paths = BU_VLS_INIT_ZERO;
    if (!ged_selection_list_paths(gedp, setName, &paths)) {
	bu_vls_free(&paths);
	return ret;
    }

    const char *pstr = bu_vls_cstr(&paths);
    const char *start = pstr;
    for (const char *c = pstr; c && *c; c++) {
	if (*c != '\n')
	    continue;
	if (c > start)
	    ret.push_back(std::string(start, static_cast<size_t>(c - start)));
	start = c + 1;
    }
    if (start && *start)
	ret.push_back(std::string(start));
    bu_vls_free(&paths);

    return ret;
}

static bool
qg_obol_source_selected(SoBRLDatabaseSource *source,
	const std::vector<std::string> &selectedPaths)
{
    if (!source)
	return false;

    const char *sourcePath = source->path.getValue().getString();
    for (const std::string &spath : selectedPaths) {
	if (qg_obol_path_equal(sourcePath, spath.c_str()))
	    return true;
    }

    return false;
}

static bool
qg_obol_realized_path_selected(const char *realizedPath,
	const std::vector<std::string> &selectedPaths)
{
    for (const std::string &spath : selectedPaths) {
	if (qg_obol_path_is_or_contains(realizedPath, spath.c_str()))
	    return true;
    }
    return false;
}

static int
qg_obol_set_shape_selection(SoBRLVListShape *shape, bool selected)
{
    if (!shape || shape->selected.getValue() == (selected ? TRUE : FALSE))
	return 0;

    shape->selected = selected ? TRUE : FALSE;
    return 1;
}

static int
qg_obol_set_mesh_selection(SoBRLMeshShape *mesh, bool selected)
{
    if (!mesh || mesh->selected.getValue() == (selected ? TRUE : FALSE))
	return 0;

    mesh->selected = selected ? TRUE : FALSE;
    return 1;
}

static void
qg_obol_collect_render_sources(SoNode *node,
	std::unordered_set<SoBRLDatabaseSource *> &sources)
{
    if (!node)
	return;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoBRLDatabaseSource *source =
	    static_cast<SoBRLDatabaseSource *>(node);
	sources.insert(source);
	return;
    }
    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++)
	qg_obol_collect_render_sources(group->getChild(i), sources);
}

static int
qg_obol_sync_selection_state_impl(struct ged *gedp,
	QgView *display,
	const char *setName,
	bool require_active_selection)
{
    if (!gedp || !display)
	return 0;

    if (!qg_obol_display_accepts_ged_active_view(gedp, display))
	return 0;

    BObolViewController *obol = display->obolViewController();
    if (!obol)
	return 0;

    const std::vector<std::string> selectedPaths =
	qg_obol_selection_paths(gedp, setName);
    if (require_active_selection && selectedPaths.empty())
	return 0;
    std::vector<SbString> compactSelectedPaths;
    compactSelectedPaths.reserve(selectedPaths.size());
    for (const std::string &path : selectedPaths)
	compactSelectedPaths.push_back(SbString(path.c_str()));

    std::unordered_set<SoBRLDatabaseSource *> sources;
    qg_obol_collect_render_sources(obol->getRenderSceneRoot(), sources);
    if (sources.empty())
	qg_obol_collect_render_sources(obol->getSceneRoot(), sources);

    int changed = 0;
    for (SoBRLDatabaseSource *source : sources) {
	if (!source)
	    continue;

	const bool wholeSourceSelected =
	    qg_obol_source_selected(source, selectedPaths);
	/*
	 * Store the semantic path frontier even before a compact index exists.
	 * A cold progressive source may publish its 50k occurrence registry
	 * after the user selects a tree row; installCompactInstanceIndex()
	 * reapplies this retained frontier to the new entries.  The source's
	 * sorted path index also makes this proportional to the changed
	 * selection, rather than scanning every occurrence on the GUI thread.
	 */
	changed |= source->syncCompactInstanceSelectedPaths(
	    compactSelectedPaths);
	if (source->hasCompactInstanceIndex())
	    continue;

	for (int j = 0; j < source->getRealizedShapeCount(); j++) {
	    SoBRLVListShape *shape = source->getRealizedShape(j);
	    if (!shape)
		continue;
	    bool selected = wholeSourceSelected ||
		qg_obol_realized_path_selected(
			shape->sourcePath.getValue().getString(),
			selectedPaths);
	    changed |= qg_obol_set_shape_selection(shape, selected);
	}

	for (int j = 0; j < source->getRealizedMeshCount(); j++) {
	    SoBRLMeshShape *mesh = source->getRealizedMesh(j);
	    if (!mesh)
		continue;
	    bool selected = wholeSourceSelected ||
		qg_obol_realized_path_selected(
			mesh->sourcePath.getValue().getString(),
			selectedPaths);
	    changed |= qg_obol_set_mesh_selection(mesh, selected);
	}
    }

    if (changed) {
	obol->requestRender("selection-state");
	display->need_update(QG_VIEW_REFRESH | QG_VIEW_SELECT);
    }
    return changed;
}

int
qg_obol_sync_selection_state(struct ged *gedp,
	QgView *display,
	const char *setName)
{
    return qg_obol_sync_selection_state_impl(gedp, display, setName, false);
}

int
qg_obol_sync_selection_state_if_active(struct ged *gedp,
	QgView *display,
	const char *setName)
{
    return qg_obol_sync_selection_state_impl(gedp, display, setName, true);
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
