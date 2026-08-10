/*                         W H O . C
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
/** @file libged/who.c
 *
 * The who command.
 *
 */

#include <string.h>

#include <set>
#include <string>

#include "BObol/BDatabaseSource.h"
#include "BObol/BSceneController.h"
#include "BObol/BViewController.h"
#include "ged.h"
#include "../ged_bobol_private.hpp"
#include "../ged_private.h"

extern "C" int ged_who2_core(struct ged *gedp, int argc, const char **argv);
extern "C" int ged_who_solids_core(struct ged *gedp, int argc, const char **argv);

static void
who_append_real_paths(struct ged *gedp, struct ged_view_context *view_ctx)
{
    BObolSceneController *scene = ged_bobol_scene(gedp);
    if (!scene)
	return;

    std::set<std::string> paths;
    for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !ged_bobol_source_in_view(view_ctx, summary) || !summary.visible)
	    continue;
	const char *path = summary.parentGroupPath.getLength() &&
	    !BU_STR_EQUAL(summary.parentGroupPath.getString(), "/") ?
	    summary.parentGroupPath.getString() : summary.path.getString();
	while (path && *path == '/')
	    path++;
	if (path && path[0])
	    paths.insert(path);
    }
    for (const std::string &path : paths)
	bu_vls_printf(gedp->ged_result_str, "%s ", path.c_str());
}


struct who_overlay_ctx {
    std::set<std::string> *names;
};

static int
who_overlay_record_cb(const BObolFeatureRecord &record, void *ud)
{
    struct who_overlay_ctx *ctx = static_cast<struct who_overlay_ctx *>(ud);
    if (!ctx || !ctx->names || !record.overlay.isOverlay ||
	(record.style.hasVisible && !record.style.visible) ||
	record.name.getLength() == 0)
	return 1;
    ctx->names->insert(record.name.getString());
    return 1;
}


static void
who_append_overlay_paths(struct ged *gedp, struct ged_view_context *view_ctx)
{
    std::set<std::string> names;
    struct who_overlay_ctx ctx = {&names};
    BObolViewController *local = ged_bobol_view_controller(view_ctx);
    if (local) {
	const BObolFeatureOwner owner = ged_bobol_view_feature_owner(view_ctx);
	local->features().visitRecords(who_overlay_record_cb, &ctx,
	    BOBOL_FEATURE_SCOPE_LOCAL, &owner);
	local->features().visitRecords(who_overlay_record_cb, &ctx,
	    BOBOL_FEATURE_SCOPE_SHARED, nullptr);
    }
    BObolViewController *shared = ged_bobol_shared_view_controller(gedp);
    if (shared && shared != local)
	shared->features().visitRecords(who_overlay_record_cb, &ctx,
	    BOBOL_FEATURE_SCOPE_SHARED, nullptr);

    for (const std::string &name : names)
	bu_vls_printf(gedp->ged_result_str, "%s ", name.c_str());
}


static int
who_legacy_selector(const char *arg)
{
    if (!arg)
	return 0;

    return BU_STR_EQUAL(arg, "r") || BU_STR_EQUAL(arg, "real") ||
	BU_STR_EQUAL(arg, "p") || BU_STR_EQUAL(arg, "phony") ||
	BU_STR_EQUAL(arg, "b") || BU_STR_EQUAL(arg, "both");
}

/*
 * List the objects currently prepped for drawing
 *
 * Usage:
 * who [r(eal)|p(hony)|b(oth)]
 * who solids [level]
 *
 */
extern "C" int
ged_who_core(struct ged *gedp, int argc, const char *argv[])
{
    if (argc > 1 && (BU_STR_EQUAL(argv[1], "solids") || BU_STR_EQUAL(argv[1], "report")))
	return ged_who_solids_core(gedp, argc, argv);

    if (argc != 2 || !who_legacy_selector(argv[1]))
	return ged_who2_core(gedp, argc, argv);


    GED_CHECK_DRAWABLE(gedp, BRLCAD_ERROR);
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    int list_real = 1;
    int list_overlays = 0;
    switch (argv[1][0]) {
	case 'b':
	    list_real = 1;
	    list_overlays = 1;
	    break;
	case 'p':
	    list_real = 0;
	    list_overlays = 1;
	    break;
	case 'r':
	default:
	    list_real = 1;
	    list_overlays = 0;
	    break;
    }

    if (list_real)
	who_append_real_paths(gedp, ged_view_active_ctx(gedp));

    if (list_overlays)
	who_append_overlay_paths(gedp, ged_view_active_ctx(gedp));

    return BRLCAD_OK;
}


#include "../include/plugin.h"

#define GED_WHO_COMMANDS(X, XID) \
    X(who, ged_who_core, GED_CMD_DEFAULT) \

GED_DECLARE_COMMAND_SET(GED_WHO_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_who", 1, GED_WHO_COMMANDS)

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
