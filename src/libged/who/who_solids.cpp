/*                    W H O _ S O L I D S . C P P
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
/** @file libged/who_solids.cpp
 *
 * Shared detailed solid reporting support for who/solid_report/x.
 *
 */

#include "common.h"

#include <algorithm>
#include <cstring>
#include <set>
#include <string>
#include <vector>

extern "C" {
#include "bu/opt.h"
#include "ged.h"
#include "ged/draw.h"
#include "rt/db_fullpath.h"
#include "rt/directory.h"
#include "rt/tree.h"
}

#include "BObol/BExportAction.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BSceneController.h"
#include "BObol/BViewController.h"
#include <Inventor/SoViewport.h>
#include "../alphanum.h"
#include "../ged_bobol_private.hpp"
#include "../ged_private.h"

struct who_solids_record {
    std::string path;
    int highlighted;
    unsigned char color[3];
    point_t bounds_center;
    fastf_t bounds_radius;
    size_t vlist_structure_count;
    size_t vlist_point_count;
    std::string geometry_report;
};

static void
who_solids_usage(struct ged *gedp, const char *cmd_name, int subcmd_usage)
{
    const char *subcmd = "  who solids [-V view] [-m #] [level]\n";

    if (subcmd_usage) {
	bu_vls_printf(gedp->ged_result_str,
		"Usage:\n"
		"%s"
		"\n"
		"level <= -2 : print primitive names\n"
		"level == -1 : print full paths\n"
		"level == 0  : print full paths and ILLUM state\n"
		"level == 1  : add center, size, region, and color info\n"
		"level == 2  : add vector list counts\n"
		"level >= 3  : add vector list commands\n",
		subcmd);
	return;
    }

    bu_vls_printf(gedp->ged_result_str,
	    "Usage:\n"
	    "  %s [level]\n"
	    "%s"
	    "\n"
	    "level <= -2 : print primitive names\n"
	    "level == -1 : print full paths\n"
	    "level == 0  : print full paths and ILLUM state\n"
	    "level == 1  : add center, size, region, and color info\n"
	    "level == 2  : add vector list counts\n"
	    "level >= 3  : add vector list commands\n",
	    cmd_name,
	    subcmd);
}


static int
who_solids_level(const char *arg, int *lvl)
{
    if (!arg || !lvl)
	return 0;

    const char *av[1] = {arg};
    if (bu_opt_int(NULL, 1, av, (void *)lvl) != 1)
	return 0;

    if (*lvl > 3)
	*lvl = 3;

    return 1;
}


static int
who_solids_path_state(struct db_i *dbip, const char *path, struct db_tree_state *ts)
{
    struct db_full_path fp;
    int ret = 0;

    if (!dbip || !path || !ts)
	return 0;

    db_full_path_init(&fp);
    RT_DBTS_INIT(ts);
    db_init_db_tree_state(ts, dbip);

    ret = (db_follow_path_for_state(ts, &fp, path, LOOKUP_QUIET) >= 0);

    db_free_full_path(&fp);
    if (!ret)
	db_free_db_tree_state(ts);

    return ret;
}


static void
who_solids_print_record(const struct who_solids_record &rec,
			struct db_i *dbip,
			int lvl,
			struct bu_vls *vls)
{
    struct db_tree_state ts = RT_DBTS_INIT_ZERO;
    unsigned char basecolor[3];
    int region_id = -1;
    int dflag = 0;
    int cflag = 0;
    const char *path = NULL;
    const char *leaf = NULL;

    if (!dbip || !vls)
	return;

    path = rec.path.c_str();
    if (lvl <= -2) {
	if (path && *path) {
	    leaf = strrchr(path, '/');
	    leaf = leaf ? leaf + 1 : path;
	}
	if (leaf)
	    bu_vls_printf(vls, "%s ", leaf);
	return;
    }

    bu_vls_printf(vls, "%s", path ? path : "");
    if ((lvl != -1) && rec.highlighted)
	bu_vls_printf(vls, " ILLUM");
    bu_vls_printf(vls, "\n");

    if (lvl <= 0)
	return;

    basecolor[0] = rec.color[0];
    basecolor[1] = rec.color[1];
    basecolor[2] = rec.color[2];
    if (who_solids_path_state(dbip, path, &ts)) {
	region_id = ts.ts_regionid;
	if (ts.ts_mater.ma_color_valid) {
	    basecolor[0] = ts.ts_mater.ma_color[0] * 255.0;
	    basecolor[1] = ts.ts_mater.ma_color[1] * 255.0;
	    basecolor[2] = ts.ts_mater.ma_color[2] * 255.0;
	} else {
	    dflag = 1;
	}
	db_free_db_tree_state(&ts);
    }

    cflag = (basecolor[0] != rec.color[0] ||
	     basecolor[1] != rec.color[1] ||
	     basecolor[2] != rec.color[2]);

    bu_vls_printf(vls, "  cent=(%.3f, %.3f, %.3f) sz=%g ",
		  rec.bounds_center[X] * dbip->dbi_base2local,
		  rec.bounds_center[Y] * dbip->dbi_base2local,
		  rec.bounds_center[Z] * dbip->dbi_base2local,
		  rec.bounds_radius * dbip->dbi_base2local);
    bu_vls_printf(vls, "reg=%d\n", region_id);
    bu_vls_printf(vls, "  basecolor=(%d, %d, %d) color=(%d, %d, %d)%s%s\n",
		  basecolor[0],
		  basecolor[1],
		  basecolor[2],
		  rec.color[0],
		  rec.color[1],
		  rec.color[2],
		  dflag ? " D" : "",
		  cflag ? " C" : "");

    if (lvl <= 1)
	return;

    if (lvl > 2)
	bu_vls_printf(vls, "%s", rec.geometry_report.c_str());

    bu_vls_printf(vls, "  %zu vlist structures, %zu pts\n",
		  rec.vlist_structure_count,
		  rec.vlist_point_count);
    bu_vls_printf(vls, "  %zu pts (via semantic export)\n",
		  rec.vlist_point_count);
}

static unsigned char
who_solids_color_channel(float value)
{
    if (value <= 0.0f)
	return 0;
    if (value >= 1.0f)
	return 255;
    return static_cast<unsigned char>(value * 255.0f + 0.5f);
}

static void
who_solids_collect_record(const SoBRLExportAction::ObjectRecord &rec,
	const SoBRLExportAction &export_action, int lvl,
	std::vector<struct who_solids_record> &objects)
{
    if (rec.path.getLength() == 0 || rec.nonDatabaseSource ||
	rec.overlayIntent || rec.hudIntent || rec.localSource || rec.sharedSource)
	return;

    struct who_solids_record out;
    out.path = rec.path.getString();
    out.highlighted = rec.highlighted;
    SbColor color = rec.color;
    if (!rec.lineIndices.empty()) {
	const SoBRLExportAction::LineRecord &line =
	    export_action.getLine(rec.lineIndices.front());
	color = line.colorOverride ? line.color :
	    (line.materialColorValid ? line.materialColor : line.color);
    } else if (!rec.pointIndices.empty()) {
	const SoBRLExportAction::PointRecord &point =
	    export_action.getPoint(rec.pointIndices.front());
	color = point.colorOverride ? point.color :
	    (point.materialColorValid ? point.materialColor : point.color);
    } else if (!rec.triangleIndices.empty()) {
	const SoBRLExportAction::TriangleRecord &triangle =
	    export_action.getTriangle(rec.triangleIndices.front());
	color = triangle.colorOverride ? triangle.color :
	    (triangle.materialColorValid ? triangle.materialColor :
	     triangle.color);
    }
    out.color[0] = who_solids_color_channel(color[0]);
    out.color[1] = who_solids_color_channel(color[1]);
    out.color[2] = who_solids_color_channel(color[2]);
    VSETALL(out.bounds_center, 0.0);
    out.bounds_radius = 0.0;
    if (!rec.bounds.isEmpty()) {
	const SbVec3f center = rec.bounds.getCenter();
	const SbVec3f diagonal = rec.bounds.getMax() - rec.bounds.getMin();
	VSET(out.bounds_center, center[0], center[1], center[2]);
	out.bounds_radius = 0.5 * diagonal.length();
    }
    out.vlist_structure_count = rec.lineIndices.size() +
	rec.pointIndices.size() + rec.triangleIndices.size();
    out.vlist_point_count = rec.lineIndices.size() * 2 +
	rec.pointIndices.size() + rec.triangleIndices.size() * 3;

    if (lvl > 2) {
	struct bu_vls report = BU_VLS_INIT_ZERO;
	for (int index : rec.lineIndices) {
	    const SoBRLExportAction::LineRecord &line =
		export_action.getLine(index);
	    bu_vls_printf(&report,
		"  line move (%g, %g, %g)\n  line draw (%g, %g, %g)\n",
		line.a[0], line.a[1], line.a[2],
		line.b[0], line.b[1], line.b[2]);
	}
	for (int index : rec.pointIndices) {
	    const SbVec3f &point = export_action.getPoint(index).point;
	    bu_vls_printf(&report, "  point draw (%g, %g, %g)\n",
		point[0], point[1], point[2]);
	}
	for (int index : rec.triangleIndices) {
	    const SoBRLExportAction::TriangleRecord &triangle =
		export_action.getTriangle(index);
	    const SbVec3f points[3] = {triangle.a, triangle.b, triangle.c};
	    for (size_t i = 0; i < 3; i++)
		bu_vls_printf(&report,
		    "  indexed face %zu (%g, %g, %g)\n", i,
		    points[i][0], points[i][1], points[i][2]);
	}
	out.geometry_report = bu_vls_cstr(&report);
	bu_vls_free(&report);
    }
    objects.push_back(out);
}

static void
who_solids_print_view(struct ged *gedp, struct ged_view_context *v,
	struct db_i *dbip,
	int mode, int lvl, struct bu_vls *vls)
{
    if (!v)
	return;

    BObolViewController *controller = ged_bobol_view_controller(v);
    SoBRLExportAction export_action;
    std::vector<SoBRLExportAction::ObjectRecord> objects;
    if (controller && controller->getViewport() &&
	controller->getViewport()->getRoot()) {
	(void)controller->realizePending();
	export_action.setGeometryPolicy(SoBRLExportAction::DISPLAY_LEVEL);
	export_action.apply(controller->getViewport()->getRoot());
	export_action.collectObjectRecords(objects,
	    SoBRLExportAction::QUERY_VISIBLE_ONLY, nullptr, mode);
    }
    std::vector<struct who_solids_record> records;
    std::set<std::string> recorded_instances;
    for (const SoBRLExportAction::ObjectRecord &object : objects)
	if (!object.ownerSourceInstanceKey.getLength() ||
	    recorded_instances.insert(
		object.ownerSourceInstanceKey.getString()).second)
	    who_solids_collect_record(object, export_action, lvl, records);

    /* A source may be valid and visible before it has published display-level
     * geometry.  Source inventory is the authoritative database-object list;
     * the semantic exporter above supplements it with realized detail. */
    std::vector<SoBRLDatabaseSource *> sources = controller ?
	controller->getRenderDatabaseSources() :
	std::vector<SoBRLDatabaseSource *>();
    if (sources.empty()) {
	BObolSceneController *scene = ged_bobol_scene(gedp);
	if (scene) {
	    for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
		SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
		if (source)
		    sources.push_back(source);
	    }
	}
    }
    for (SoBRLDatabaseSource *source : sources) {
	BObolDatabaseSourceSummary summary;
	if (!source || !source->getSummary(summary) || !summary.valid ||
	    !summary.visible || (mode >= 0 && summary.drawMode != mode) ||
	    (summary.instanceKey.getLength() &&
	     recorded_instances.find(summary.instanceKey.getString()) !=
		recorded_instances.end()))
	    continue;

	struct who_solids_record out;
	out.path = summary.path.getString();
	out.highlighted = summary.highlighted ? 1 : 0;
	const SbColor color = summary.colorOverride ? summary.color :
	    (summary.materialColorValid ? summary.materialColor : summary.color);
	out.color[0] = who_solids_color_channel(color[0]);
	out.color[1] = who_solids_color_channel(color[1]);
	out.color[2] = who_solids_color_channel(color[2]);
	VSETALL(out.bounds_center, 0.0);
	out.bounds_radius = 0.0;
	if (summary.sourceBoundsValid && !summary.sourceBounds.isEmpty()) {
	    const SbVec3f center = summary.sourceBounds.getCenter();
	    const SbVec3f diagonal = summary.sourceBounds.getMax() -
		summary.sourceBounds.getMin();
	    VSET(out.bounds_center, center[0], center[1], center[2]);
	    out.bounds_radius = 0.5 * diagonal.length();
	}
	out.vlist_structure_count = 0;
	out.vlist_point_count = 0;
	for (int i = 0; i < source->getRealizedShapeSummaryCount(); i++) {
	    BObolRealizedShapeSummary shape;
	    if (!source->getRealizedShapeSummary(i, shape) || !shape.valid)
		continue;
	    out.vlist_structure_count += static_cast<size_t>(
		shape.segmentCount + shape.pointPrimitiveCount +
		shape.triangleCount);
	    out.vlist_point_count += static_cast<size_t>(shape.pointCount);
	}
	records.push_back(out);
    }

    std::sort(records.begin(), records.end(),
	      [](const struct who_solids_record &a,
		 const struct who_solids_record &b) {
		  return alphanum_impl(a.path.c_str(), b.path.c_str(), NULL) < 0;
	      });

    for (size_t i = 0; i < records.size(); i++)
	who_solids_print_record(records[i], dbip, lvl, vls);
}


static int
who_solids_impl(struct ged *gedp, int argc, const char *argv[], int subcmd_usage)
{
    int print_help = 0;
    int lvl = 0;
    int mode = -1;
    struct bu_vls cvls = BU_VLS_INIT_ZERO;
    const char *cmd_name = argv[0];
    struct bu_opt_desc d[5];
    BU_OPT(d[0], "h", "help", "",     NULL,        &print_help, "Print help and exit");
    BU_OPT(d[1], "?", "",     "",     NULL,        &print_help, "");
    BU_OPT(d[2], "V", "view", "name", &bu_opt_vls, &cvls,       "Specify view to report");
    BU_OPT(d[3], "m", "mode", "#",    &bu_opt_int, &mode,       "Only report objects drawn in the specified drawing mode");
    BU_OPT_NULL(d[4]);

    GED_CHECK_DATABASE_OPEN(gedp, BRLCAD_ERROR);
    GED_CHECK_DRAWABLE(gedp, BRLCAD_ERROR);
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    bu_vls_trunc(gedp->ged_result_str, 0);

    if (subcmd_usage) {
	if (argc < 2) {
	    who_solids_usage(gedp, cmd_name, 1);
	    bu_vls_free(&cvls);
	    return BRLCAD_ERROR;
	}
	argc--;
	argv++;
	cmd_name = argv[0];
	argc--;
	argv++;
    } else {
	argc--;
	argv++;
    }

    int opt_ret = bu_opt_parse(NULL, argc, argv, d);
    if (opt_ret < 0) {
	who_solids_usage(gedp, cmd_name, subcmd_usage);
	bu_vls_free(&cvls);
	return BRLCAD_ERROR;
    }

    argc = opt_ret;

    if (print_help) {
	who_solids_usage(gedp, subcmd_usage ? "who" : cmd_name, subcmd_usage);
	bu_vls_free(&cvls);
	return BRLCAD_OK;
    }

    if (argc > 1 || (argc == 1 && !who_solids_level(argv[0], &lvl))) {
	who_solids_usage(gedp, subcmd_usage ? "who" : cmd_name, subcmd_usage);
	bu_vls_free(&cvls);
	return BRLCAD_ERROR;
    }

    struct ged_view_context *v = ged_view_active_ctx(gedp);
    if (bu_vls_strlen(&cvls)) {
	v = ged_view_find_ctx(gedp, bu_vls_cstr(&cvls));
	if (!v) {
	    bu_vls_printf(gedp->ged_result_str, "Specified view %s not found\n", bu_vls_cstr(&cvls));
	    bu_vls_free(&cvls);
	    return BRLCAD_ERROR;
	}
    }
    bu_vls_free(&cvls);

    if (!v) {
	bu_vls_printf(gedp->ged_result_str, "No view specified and no current view defined in GED");
	return BRLCAD_ERROR;
    }

    who_solids_print_view(gedp, v, gedp->dbip, mode, lvl,
	gedp->ged_result_str);
    return BRLCAD_OK;
}


extern "C" int
ged_solid_report_shared_core(struct ged *gedp, int argc, const char *argv[])
{
    return who_solids_impl(gedp, argc, argv, 0);
}


extern "C" int
ged_who_solids_core(struct ged *gedp, int argc, const char *argv[])
{
    return who_solids_impl(gedp, argc, argv, 1);
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
