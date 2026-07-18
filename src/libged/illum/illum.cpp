/*                       I L L U M . C P P
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
/** @file libged/illum/illum.cpp
 *
 * The illum command.
 *
 */

#include "common.h"
#include <string.h>
#include <vector>

#include "BObol/BExportAction.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include <Inventor/SoViewport.h>
#include "bv.h"
#include "ged.h"
#include "ged/draw.h"
#include "ged/view.h"
#include "../ged_bobol_private.hpp"
#include "../ged_private.h"

/* Callback data for labelvert */
struct labelvert_data {
    struct directory *dp;
    struct db_i *dbip;
    double base2local;
    std::vector<BObolLabel> labels;
};

#define LABELVERT_FEATURE_NAME "_LABELVERT_ffffff"

static int
labelvert_record_matches(struct db_i *dbip,
			 const char *path,
			 struct directory *dp)
{
    if (!dbip || !path || !path[0] || !dp)
	return 0;

    struct db_full_path fp;
    db_full_path_init(&fp);
    if (db_string_to_path(&fp, dbip, path) < 0) {
	db_free_full_path(&fp);
	return 0;
    }
    int ret = db_full_path_search(&fp, dp);
    db_free_full_path(&fp);
    return ret;
}

static int
labelvert_append_label(struct labelvert_data *lvd, const point_t pt)
{
    char label[256];

    if (!lvd)
	return 0;

    snprintf(label, sizeof(label), " %g, %g, %g",
	    pt[0] * lvd->base2local,
	    pt[1] * lvd->base2local,
	    pt[2] * lvd->base2local);
    BObolLabel item;
    item.text = label;
    item.point = SbVec3f(static_cast<float>(pt[X]),
	static_cast<float>(pt[Y]), static_cast<float>(pt[Z]));
    item.hasColor = TRUE;
    item.color = SbColor(1.0f, 1.0f, 1.0f);
    lvd->labels.push_back(item);
    return 1;
}

/* Usage:  labelvert solid(s) */
int
ged_labelvert_core(struct ged *gedp, int argc, const char *argv[])
{
    int i;
    struct labelvert_data lvd;
    static const char *usage = "object(s) - label vertices of wireframes of objects";

    if (!gedp || !gedp->dbip)
	return BRLCAD_ERROR;

    if (argc < 2) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, ": no current view set");
	return BRLCAD_ERROR;
    }

    lvd.dp = nullptr;
    lvd.dbip = gedp->dbip;
    lvd.base2local = gedp->dbip->dbi_base2local;

    BObolViewController *endpoint = ged_bobol_view_controller(view_ctx);
    if (!endpoint || !endpoint->getViewport() ||
	!endpoint->getViewport()->getRoot()) {
	bu_vls_printf(gedp->ged_result_str,
	    "labelvert: current view has no BObol display endpoint");
	return BRLCAD_ERROR;
    }
    SoBRLExportAction export_action;
    export_action.setGeometryPolicy(SoBRLExportAction::DISPLAY_LEVEL);
    export_action.apply(endpoint->getViewport()->getRoot());
    std::vector<SoBRLExportAction::ObjectRecord> records;
    export_action.collectObjectRecords(records,
	SoBRLExportAction::QUERY_VISIBLE_ONLY |
	SoBRLExportAction::QUERY_DATABASE_OBJECTS);

    for (i=1; i<argc; i++) {
	struct directory *dp;
	if ((dp = db_lookup(gedp->dbip, argv[i], LOOKUP_NOISY)) == RT_DIR_NULL)
	    continue;
	/* Find displayed uses of this database object. */
	lvd.dp = dp;
	for (const SoBRLExportAction::ObjectRecord &record : records) {
	    if (!record.databaseIntent ||
		!labelvert_record_matches(lvd.dbip,
		    record.path.getString(), lvd.dp))
		continue;
	    const auto append = [&lvd](const SbVec3f &point) {
		point_t model_point;
		VSET(model_point, point[0], point[1], point[2]);
		(void)labelvert_append_label(&lvd, model_point);
	    };
	    for (int index : record.lineIndices) {
		const SoBRLExportAction::LineRecord &line =
		    export_action.getLine(index);
		append(line.a);
		append(line.b);
	    }
	    for (int index : record.pointIndices)
		append(export_action.getPoint(index).point);
	    for (int index : record.triangleIndices) {
		const SoBRLExportAction::TriangleRecord &triangle =
		    export_action.getTriangle(index);
		append(triangle.a);
		append(triangle.b);
		append(triangle.c);
	    }
	}
    }

    BObolViewController *controller =
	ged_bobol_shared_view_controller(gedp);
    BObolFeatureHandle existing = controller ?
	controller->features().find(LABELVERT_FEATURE_NAME,
	    BOBOL_FEATURE_SCOPE_SHARED) : BObolFeatureHandle();
    if (controller && existing.isValid())
	(void)controller->features().remove(existing);
    const bool published = controller && (lvd.labels.empty() ||
	controller->features().publishLabels(LABELVERT_FEATURE_NAME,
	    BObolFeatureScope::Shared, lvd.labels).isValid());
    if (!published) {
	bu_vls_printf(gedp->ged_result_str, "failed to create labelvert feature\n");
	return BRLCAD_ERROR;
    }

    (void)bv_context_refresh_request(ged_view_context_bv(view_ctx),
	GED_VIEW_REFRESH_DRAW);
    return BRLCAD_OK;
}


/*
 * Illuminate/highlight database object
 *
 * Usage:
 * illum [-n] obj
 *
 */
int
ged_illum_core(struct ged *gedp, int argc, const char *argv[])
{
    int illum = 1;
    int changed = 0;
    struct ged_draw_transaction txn;
    struct ged_draw_transaction_result result;
    static const char *usage = "[-n] obj";

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

    if (argc == 3) {
	if (argv[1][0] == '-' && argv[1][1] == 'n')
	    illum = 0;
	else
	    goto bad;

	--argc;
	++argv;
    }

    if (argc != 2)
	goto bad;

    txn = ged_draw_transaction_make_value(GED_DRAW_TXN_HIGHLIGHT,
	    argv[1], (double)illum);
    ged_draw_transaction_result_init(&result);
    changed = ged_draw_apply_transaction(gedp, &txn, &result);
    ged_draw_transaction_result_free(&result);
    if (!changed) {
	bu_vls_printf(gedp->ged_result_str, "illum: %s not found", argv[1]);
	return BRLCAD_ERROR;
    }

    return BRLCAD_OK;

bad:
    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
    return BRLCAD_ERROR;
}

#include "../include/plugin.h"

#define GED_ILLUM_COMMANDS(X, XID) \
    X(illum, ged_illum_core, GED_CMD_DEFAULT) \
    X(labelvert, ged_labelvert_core, GED_CMD_DEFAULT) \

GED_DECLARE_COMMAND_SET(GED_ILLUM_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_illum", 1, GED_ILLUM_COMMANDS)

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
