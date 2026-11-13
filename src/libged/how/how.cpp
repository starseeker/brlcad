/*                         H O W . C
 * BRL-CAD
 *
 * Copyright (c) 2008-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file libged/how.c
 *
 * The how command.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <string.h>

#include "BObol/BDatabaseSource.h"
#include "BObol/BSceneController.h"
#include "bu/cmd.h"
#include "bu/str.h"
#include "../ged_bobol_private.hpp"
#include "../ged_private.h"

static const char *
how_skip_leading_slash(const char *path)
{
    while (path && *path == '/')
	path++;
    return path;
}


static bool
how_path_matches(const char *candidate, const char *query)
{
    candidate = how_skip_leading_slash(candidate);
    query = how_skip_leading_slash(query);
    if (!candidate || !query)
	return false;
    if (BU_STR_EQUAL(candidate, query))
	return true;
    const size_t n = strlen(query);
    return strlen(candidate) > n && !bu_strncmp(candidate, query, n) &&
	candidate[n] == '/';
}


static bool
how_bobol_source_find(struct ged *gedp, struct ged_view_context *view_ctx,
	const char *path, BObolDatabaseSourceSummary &match)
{
    BObolSceneController *scene = ged_bobol_scene(gedp);
    if (!scene || !path)
	return false;

    bool have_fallback = false;
    for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !ged_bobol_source_in_view(view_ctx, summary) || !summary.visible)
	    continue;
	const char *group_path = summary.parentGroupPath.getString();
	const char *source_path = summary.path.getString();
	if (!how_path_matches(group_path, path) &&
	    !how_path_matches(source_path, path))
	    continue;

	if ((group_path && BU_STR_EQUAL(how_skip_leading_slash(group_path),
		how_skip_leading_slash(path))) ||
	    (source_path && BU_STR_EQUAL(how_skip_leading_slash(source_path),
		how_skip_leading_slash(path)))) {
	    match = summary;
	    return true;
	}
	if (!have_fallback) {
	    match = summary;
	    have_fallback = true;
	}
    }
    return have_fallback;
}

/*
 * Returns "how" an object is being displayed.
 *
 * Usage:
 * how [-b] object
 *
 */
extern "C" int
ged_how_core(struct ged *gedp, int argc, const char *argv[])
{
    int both = 0;
    int dmode = 0;
    fastf_t transparency = 0.0;
    static const char *usage = "[-b] object";
    const char *obj_arg = NULL;

    GED_CHECK_DATABASE_OPEN(gedp, BRLCAD_ERROR);
    GED_CHECK_DRAWABLE(gedp, BRLCAD_ERROR);
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (3 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    if (argc == 3 &&
	argv[1][0] == '-' &&
	argv[1][1] == 'b') {
	both = 1;
    } else {
	if (argc != 2) {
	    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	    return BRLCAD_ERROR;
	}
    }

    obj_arg = both ? argv[2] : argv[1];
    BObolDatabaseSourceSummary match;
    if (!how_bobol_source_find(gedp, ged_view_active_ctx(gedp), obj_arg,
	    match))
	goto not_found;

    dmode = match.representationMode >= 0 ?
	match.representationMode : GED_DRAW_MODE_WIRE;
    transparency = match.transparency;
    if (dmode == _GED_HIDDEN_LINE) {
	if (both)
	    bu_vls_printf(gedp->ged_result_str, "%d 1", _GED_HIDDEN_LINE);
	else
	    bu_vls_printf(gedp->ged_result_str, "%d", _GED_HIDDEN_LINE);
    } else {
	if (both)
	    bu_vls_printf(gedp->ged_result_str, "%d %g", dmode, transparency);
	else
	    bu_vls_printf(gedp->ged_result_str, "%d", dmode);
    }
    return BRLCAD_OK;

not_found:
    bu_vls_printf(gedp->ged_result_str, "-1");

    return BRLCAD_OK;
}

#include "../include/plugin.h"

#define GED_HOW_COMMANDS(X, XID) \
    X(how, ged_how_core, GED_CMD_DEFAULT) \

GED_DECLARE_COMMAND_SET(GED_HOW_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_how", 1, GED_HOW_COMMANDS)

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
