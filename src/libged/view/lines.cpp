/*                      L I N E S . C P P
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
/** @file libged/view/lines.cpp
 *
 * Commands for view lines.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <vector>

#include "BObol/BViewController.h"
#include "bu/cmd.h"
#include "bu/color.h"
#include "bu/opt.h"
#include "bu/vls.h"

#include "../ged_private.h"
#include "../ged_bobol_private.hpp"
#include "./ged_view.h"

static std::vector<SbVec3f>
_line_points(const point_t *points, size_t point_count)
{
    std::vector<SbVec3f> out;
    out.reserve(point_count);
    for (size_t i = 0; points && i < point_count; i++)
	out.emplace_back(static_cast<float>(points[i][X]),
	    static_cast<float>(points[i][Y]),
	    static_cast<float>(points[i][Z]));
    return out;
}

static std::vector<int32_t>
_line_commands(size_t point_count)
{
    std::vector<int32_t> commands(point_count,
	static_cast<int32_t>(BObolLineCommand::Draw));
    if (!commands.empty())
	commands[0] = static_cast<int32_t>(BObolLineCommand::Move);
    return commands;
}

static BObolOverlayInfo
_line_overlay(struct ged_view_context *view_ctx, const char *name)
{
    BObolOverlayInfo overlay;
    overlay.isOverlay = TRUE;
    overlay.ownerToken = view_ctx;
    overlay.role = BObolOverlayRole::Model;
    overlay.overlayClass = BObolOverlayClass::UserAnnotation;
    overlay.lifecycle = BObolOverlayLifecycle::Persistent;
    overlay.order = BObolOverlayOrder::Model;
    overlay.sourcePath = name ? name : "";
    return overlay;
}

static int
_line_parse_points(struct ged *gedp, int argc, const char **argv, point_t **points_out, size_t *point_count_out, const char *usage_string)
{
    if (points_out)
	*points_out = NULL;
    if (point_count_out)
	*point_count_out = 0;
    if (!points_out || !point_count_out)
	return BRLCAD_ERROR;
    if (argc <= 0 || argc % 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    size_t point_count = (size_t)argc / 3;
    point_t *points = (point_t *)bu_calloc(point_count, sizeof(point_t),
	    "view annotation line points");
    for (size_t i = 0; i < point_count; i++) {
	for (int axis = 0; axis < 3; axis++) {
	    const char **av = &argv[(i * 3) + (size_t)axis];
	    if (bu_opt_fastf_t(NULL, 1, av, (void *)&(points[i][axis])) != 1) {
		bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", *av);
		bu_free(points, "view annotation line points");
		return BRLCAD_ERROR;
	    }
	}
    }

    *points_out = points;
    *point_count_out = point_count;
    return BRLCAD_OK;
}

int
_line_cmd_create(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view annotation line create <name> x y z [x y z ...]";
    const char *purpose_string = "create a polyline from one or more model-space points";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    point_t *points = NULL;
    size_t point_count = 0;
    if (_line_parse_points(gedp, argc, argv, &points, &point_count,
	    usage_string) != BRLCAD_OK)
	return BRLCAD_ERROR;

    int ret = 0;
    if (ged_bobol_feature_find(gd->gedp, gd->cv,
	    gd->vobj).handle.isValid()) {
	bu_vls_printf(gedp->ged_result_str,
	    "View feature named %s already exists\n", gd->vobj);
    } else {
	BObolViewController *controller = ged_bobol_feature_controller(
	    gd->gedp, gd->cv, gd->local_obj);
	if (controller) {
	    BObolFeatureOwner owner = ged_bobol_view_feature_owner(gd->cv);
	    const BObolFeatureScope scope = gd->local_obj ?
		BObolFeatureScope::Local : BObolFeatureScope::Shared;
	    BObolFeatureStyle style;
	    style.hasColor = TRUE;
	    style.color = SbColor(1.0f, 0.0f, 0.0f);
	    BObolFeatureHandle handle = controller->features().publishLineSet(
		gd->vobj, scope, _line_points(points, point_count),
		_line_commands(point_count), &style,
		gd->local_obj ? &owner : nullptr);
	    if (handle.isValid()) {
		(void)controller->features().setOverlayInfo(handle,
		    _line_overlay(gd->cv, gd->vobj));
		ret = 1;
	    }
	}
	if (!ret)
	    bu_vls_printf(gedp->ged_result_str, "Failed to create %s\n",
		gd->vobj);
    }
    bu_free(points, "view annotation line points");
    return ret ? BRLCAD_OK : BRLCAD_ERROR;
}

int
_line_cmd_append(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view annotation line append <name> x y z [x y z ...]";
    const char *purpose_string = "append one or more model-space points to a polyline";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    point_t *points = NULL;
    size_t point_count = 0;
    if (_line_parse_points(gedp, argc, argv, &points, &point_count,
	    usage_string) != BRLCAD_OK)
	return BRLCAD_ERROR;

    struct ged_bobol_feature_binding binding =
	ged_bobol_feature_find(gd->gedp, gd->cv, gd->vobj);
    int ret = binding.controller && binding.handle.isValid();
    if (!ret) {
	bu_vls_printf(gedp->ged_result_str,
	    "no view feature named %s\n", gd->vobj);
    } else {
	for (size_t i = 0; i < point_count; i++) {
	    if (!binding.controller->features().appendLinePoint(binding.handle,
		SbVec3f(static_cast<float>(points[i][X]),
		    static_cast<float>(points[i][Y]),
		    static_cast<float>(points[i][Z])))) {
		bu_vls_printf(gedp->ged_result_str,
		    "Failed to append point %zu to %s\n", i, gd->vobj);
		ret = 0;
		break;
	    }
	}
    }
    bu_free(points, "view annotation line points");
    return ret ? BRLCAD_OK : BRLCAD_ERROR;
}

int
_line_cmd_remove(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view annotation line remove <name> <zero-based-index> [count]";
    const char *purpose_string = "remove one or more points from a polyline";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (argc < 1 || argc > 2) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    int start = -1;
    if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&start) != 1 || start < 0) {
	bu_vls_printf(gedp->ged_result_str, "Invalid zero-based point index %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    int remove_count_arg = 1;
    if (argc == 2 &&
	    (bu_opt_int(NULL, 1, (const char **)&argv[1], (void *)&remove_count_arg) != 1 ||
	     remove_count_arg <= 0)) {
	bu_vls_printf(gedp->ged_result_str, "Invalid point count %s\n", argv[1]);
	return BRLCAD_ERROR;
    }

    struct ged_bobol_feature_binding binding =
	ged_bobol_feature_find(gd->gedp, gd->cv, gd->vobj);
    if (!binding.controller || !binding.handle.isValid()) {
	bu_vls_printf(gedp->ged_result_str,
	    "no view feature named %s\n", gd->vobj);
	return BRLCAD_ERROR;
    }

    BObolFeatureRecord record;
    if (!binding.controller->features().record(binding.handle, record)) {
	bu_vls_printf(gedp->ged_result_str,
	    "Failed to read points for %s\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    const size_t remove_start = static_cast<size_t>(start);
    const size_t remove_count = static_cast<size_t>(remove_count_arg);
    if (remove_start >= record.points.size() ||
	remove_count > record.points.size() - remove_start) {
	bu_vls_printf(gedp->ged_result_str,
	    "Point range [%zu, %zu) is outside %s point count %zu\n",
	    remove_start, remove_start + remove_count, gd->vobj,
	    record.points.size());
	return BRLCAD_ERROR;
    }

    record.points.erase(record.points.begin() +
	static_cast<std::vector<SbVec3f>::difference_type>(remove_start),
	record.points.begin() + static_cast<std::vector<SbVec3f>::difference_type>(
	    remove_start + remove_count));
    if (record.points.empty()) {
	if (!binding.controller->features().remove(binding.handle)) {
	    bu_vls_printf(gedp->ged_result_str, "Failed to update %s\n",
		gd->vobj);
	    return BRLCAD_ERROR;
	}
	return BRLCAD_OK;
    }

    BObolFeatureHandle updated = binding.controller->features().publishLineSet(
	record.name, record.scope, record.points,
	_line_commands(record.points.size()), &record.style,
	record.scope == BObolFeatureScope::Local ? &record.owner : nullptr);
    if (!updated.isValid()) {
	bu_vls_printf(gedp->ged_result_str, "Failed to update %s\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    if (record.overlay.isOverlay)
	(void)binding.controller->features().setOverlayInfo(updated,
	    record.overlay);
    return BRLCAD_OK;
}

int
_line_cmd_clear(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view annotation line clear <name>";
    const char *purpose_string = "remove all points from a polyline";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (argc != 0) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    struct ged_bobol_feature_binding binding =
	ged_bobol_feature_find(gd->gedp, gd->cv, gd->vobj);
    if (!binding.controller || !binding.handle.isValid()) {
	bu_vls_printf(gedp->ged_result_str,
	    "no view feature named %s\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    if (!binding.controller->features().remove(binding.handle)) {
	bu_vls_printf(gedp->ged_result_str, "Failed to clear %s\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

const struct bu_cmdtab _line_cmds[] = {
    { "create",          _line_cmd_create},
    { "append",          _line_cmd_append},
    { "remove",          _line_cmd_remove},
    { "clear",           _line_cmd_clear},
    { (char *)NULL,      NULL}
};

int
_view_cmd_lines(void *bs, int argc, const char **argv)
{
    int help = 0;
    int s_version = 0;
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;

    const char *usage_string = "view annotation line [options] [args]";
    const char *purpose_string = "manipulate view lines";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    if (!gd->cv) {
	bu_vls_printf(gedp->ged_result_str, ": no view specified or current in GED");
	return BRLCAD_ERROR;
    }


    // We know we're the lines command - start processing args
    argc--; argv++;

    // See if we have any high level options set
    struct bu_opt_desc d[3];
    BU_OPT(d[0], "h", "help",  "",  NULL,  &help,      "Print help");
    BU_OPT(d[1], "s", "",      "",  NULL,  &s_version, "Work with S version of data");
    BU_OPT_NULL(d[2]);

    gd->gopts = d;

    // High level options are only defined prior to the subcommand
    int cmd_pos = -1;
    for (int i = 0; i < argc; i++) {
	if (bu_cmd_valid(_line_cmds, argv[i]) == BRLCAD_OK) {
	    cmd_pos = i;
	    break;
	}
    }

    int acnt = (cmd_pos >= 0) ? cmd_pos : argc;
    (void)bu_opt_parse(NULL, acnt, argv, d);

    return _ged_subcmd_exec(gedp, d, _line_cmds, "view annotation line", "[options] subcommand [args]", gd, argc, argv, help, cmd_pos);
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
