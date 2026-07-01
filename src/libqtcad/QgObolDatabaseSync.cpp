/*            Q G O B O L D A T A B A S E S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolDatabaseSync.cpp */

#include "common.h"

#include "QgLegacyViewContext.h"
#include "QgObolDatabaseSyncPrivate.h"

#include "brlobol/database_source.h"
#include "brlobol/view_controller.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgView.h"
#include "rt/view.h"

#include <stdio.h>
#include <string.h>

#include <string>
#include <vector>

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

static const char *
qg_obol_skip_leading_slash(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

static int
qg_obol_path_equal(const char *a, const char *b)
{
    if (!a || !b)
	return 0;
    if (strcmp(a, b) == 0)
	return 1;
    return strcmp(qg_obol_skip_leading_slash(a),
	    qg_obol_skip_leading_slash(b)) == 0;
}

static void *
qg_obol_database_view_context(QgView *display)
{
    if (!display || !display->view())
	return NULL;
    return qg_legacy_view_to_context(display->view());
}

static int
qg_obol_database_view_is_independent(QgView *display)
{
    return rt_view_context_is_independent(
	    qg_obol_database_view_context(display));
}

static std::string
qg_obol_database_view_scope_name(QgView *display)
{
    void *view_ctx = qg_obol_database_view_context(display);
    if (!view_ctx)
	return std::string("shared");

    const char *name = rt_view_context_name_get(view_ctx);
    if (name && name[0])
	return std::string(name);

    char fallback[64] = {0};
    snprintf(fallback, sizeof(fallback), "%p", view_ctx);
    return std::string(fallback);
}

static std::string
qg_obol_database_source_instance_key(QgView *display, const char *path)
{
    if (!path || !path[0])
	return std::string();

    if (!qg_obol_database_view_is_independent(display))
	return std::string(path);

    std::string key("ged-view:");
    key += qg_obol_database_view_scope_name(display);
    key += ":";
    key += qg_obol_skip_leading_slash(path);
    return key;
}

static std::string
qg_obol_database_source_instance_prefix(QgView *display)
{
    if (!qg_obol_database_view_is_independent(display))
	return std::string();

    std::string prefix("ged-view:");
    prefix += qg_obol_database_view_scope_name(display);
    prefix += ":";
    return prefix;
}

static int
qg_obol_database_source_instance_in_display_scope(
	const BRLObolDatabaseSourceSummary &summary,
	QgView *display)
{
    if (!qg_obol_database_view_is_independent(display))
	return summary.instanceKey.getLength() == 0 ||
	    qg_obol_path_equal(summary.instanceKey.getString(),
		summary.path.getString());

    const std::string prefix =
	qg_obol_database_source_instance_prefix(display);
    const char *key = summary.instanceKey.getString();
    return key && strncmp(key, prefix.c_str(), prefix.size()) == 0;
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
	const std::string instance_key =
	    qg_obol_database_source_instance_key(display, path);
	if (instance_key.empty())
	    continue;
	if (obol->replaceDatabaseSourceInstance(instance_key.c_str(), path,
		dbip, obol_draw_mode, source_revision) > 0)
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
	const std::string instance_key =
	    qg_obol_database_source_instance_key(display, path);
	if (instance_key.empty())
	    continue;
	if (obol->removeDatabaseSourceInstance(instance_key.c_str()) > 0)
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

    std::vector<std::string> instance_keys;
    for (int i = 0; i < obol->getDatabaseSourceCount(); i++) {
	BRLObolDatabaseSourceSummary summary;
	if (!obol->getDatabaseSourceSummary(i, summary) || !summary.valid ||
		!qg_obol_database_source_instance_in_display_scope(summary,
		    display))
	    continue;

	const char *instance_key = summary.instanceKey.getString();
	if (instance_key && instance_key[0])
	    instance_keys.push_back(std::string(instance_key));
	else if (summary.path.getLength() > 0)
	    instance_keys.push_back(std::string(summary.path.getString()));
    }

    int removed = 0;
    for (std::vector<std::string>::const_iterator it = instance_keys.begin();
	    it != instance_keys.end(); ++it) {
	if (obol->removeDatabaseSourceInstance(it->c_str()) > 0)
	    removed = 1;
    }

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
