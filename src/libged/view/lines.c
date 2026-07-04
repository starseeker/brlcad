/*                        L I N E S . C
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
/** @file libged/view/lines.c
 *
 * Commands for view lines.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "bu/cmd.h"
#include "bu/color.h"
#include "bu/opt.h"
#include "bu/vls.h"

#include "../ged_private.h"
#include "./ged_view.h"

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

static int *
_line_commands(size_t point_count)
{
    if (!point_count)
	return NULL;

    int *cmds = (int *)bu_calloc(point_count, sizeof(int),
	    "view annotation line commands");
    cmds[0] = GED_DRAW_VIEW_LINE_MOVE;
    for (size_t i = 1; i < point_count; i++)
	cmds[i] = GED_DRAW_VIEW_LINE_DRAW;
    return cmds;
}

static int
_line_replace_points(struct _ged_view_info *gd, const point_t *points, size_t point_count, const struct ged_draw_view_feature_style *style)
{
    int *cmds = _line_commands(point_count);
    int ret = ged_draw_view_context_lines_replace(gd->cv, gd->vobj,
	    gd->local_obj, points, cmds, point_count, style);
    if (cmds)
	bu_free(cmds, "view annotation line commands");
    return ret;
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

    if (ged_draw_view_context_feature_exists(gd->cv, gd->vobj)) {
        bu_vls_printf(gedp->ged_result_str, "View feature named %s already exists\n", gd->vobj);
        return BRLCAD_ERROR;
    }

    point_t *points = NULL;
    size_t point_count = 0;
    if (_line_parse_points(gedp, argc, argv, &points, &point_count,
	    usage_string) != BRLCAD_OK)
	return BRLCAD_ERROR;

    if (!ged_draw_view_context_lines_create_model_annotation(gd->cv, gd->vobj,
	    gd->local_obj, points[0])) {
	bu_vls_printf(gedp->ged_result_str, "Failed to create %s\n", gd->vobj);
	bu_free(points, "view annotation line points");
	return BRLCAD_ERROR;
    }

    for (size_t i = 1; i < point_count; i++) {
	if (!ged_draw_view_context_lines_append_point(gd->cv, gd->vobj,
		points[i])) {
	    (void)_line_replace_points(gd, NULL, 0, NULL);
	    bu_vls_printf(gedp->ged_result_str,
		    "Failed to append point %zu to %s\n", i, gd->vobj);
	    bu_free(points, "view annotation line points");
	    return BRLCAD_ERROR;
	}
    }

    bu_free(points, "view annotation line points");
    return BRLCAD_OK;
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

    if (!ged_draw_view_context_feature_exists(gd->cv, gd->vobj)) {
        bu_vls_printf(gedp->ged_result_str, "no view feature named %s\n", gd->vobj);
        return BRLCAD_ERROR;
    }

    point_t *points = NULL;
    size_t point_count = 0;
    if (_line_parse_points(gedp, argc, argv, &points, &point_count,
	    usage_string) != BRLCAD_OK)
	return BRLCAD_ERROR;

    for (size_t i = 0; i < point_count; i++) {
	if (!ged_draw_view_context_lines_append_point(gd->cv, gd->vobj,
		points[i])) {
	    bu_vls_printf(gedp->ged_result_str,
		    "Failed to append point %zu to %s\n", i, gd->vobj);
	    bu_free(points, "view annotation line points");
	    return BRLCAD_ERROR;
	}
    }

    bu_free(points, "view annotation line points");
    return BRLCAD_OK;
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

    if (!ged_draw_view_context_feature_exists(gd->cv, gd->vobj)) {
        bu_vls_printf(gedp->ged_result_str, "no view feature named %s\n", gd->vobj);
        return BRLCAD_ERROR;
    }

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

    point_t *points = NULL;
    size_t point_count = 0;
    if (!ged_draw_view_context_feature_points_copy(gd->cv, gd->vobj,
	    &points, &point_count)) {
	bu_vls_printf(gedp->ged_result_str, "Failed to read points for %s\n", gd->vobj);
	return BRLCAD_ERROR;
    }

    size_t remove_start = (size_t)start;
    size_t remove_count = (size_t)remove_count_arg;
    if (remove_start >= point_count || remove_count > point_count - remove_start) {
	bu_vls_printf(gedp->ged_result_str,
		"Point range [%zu, %zu) is outside %s point count %zu\n",
		remove_start, remove_start + remove_count, gd->vobj, point_count);
	if (points)
	    bu_free(points, "view annotation line copied points");
	return BRLCAD_ERROR;
    }

    struct ged_draw_view_feature_style style = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    int have_style = ged_draw_view_context_feature_style_get(gd->cv, gd->vobj,
	    &style);
    size_t new_count = point_count - remove_count;
    point_t *new_points = NULL;
    if (new_count) {
	new_points = (point_t *)bu_calloc(new_count, sizeof(point_t),
		"view annotation line rebuilt points");
	size_t out = 0;
	for (size_t i = 0; i < point_count; i++) {
	    if (i >= remove_start && i < remove_start + remove_count)
		continue;
	    VMOVE(new_points[out], points[i]);
	    out++;
	}
    }

    int ret = _line_replace_points(gd, (const point_t *)new_points,
	    new_count, have_style ? &style : NULL);

    if (new_points)
	bu_free(new_points, "view annotation line rebuilt points");
    if (points)
	bu_free(points, "view annotation line copied points");

    if (!ret) {
	bu_vls_printf(gedp->ged_result_str, "Failed to update %s\n", gd->vobj);
	return BRLCAD_ERROR;
    }

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

    if (!ged_draw_view_context_feature_exists(gd->cv, gd->vobj)) {
        bu_vls_printf(gedp->ged_result_str, "no view feature named %s\n", gd->vobj);
        return BRLCAD_ERROR;
    }

    if (!_line_replace_points(gd, NULL, 0, NULL)) {
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
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
