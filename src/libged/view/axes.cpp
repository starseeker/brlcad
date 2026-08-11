/*                      A X E S . C P P
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
/** @file libged/view/axes.cpp
 *
 * Commands for scene data axes (the faceplate axes for showing view XYZ
 * coordinate systems is handled separately - these are axes defined as data in
 * the 3D scene.)
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

struct ged_view_annotation_axes_state {
    point_t position;
    fastf_t size;
    int line_width;
    int color[3];
};

static int
_axes_publish(struct _ged_view_info *gd,
    const struct ged_view_annotation_axes_state &axes, int local)
{
    struct ged_view_feature_batch_desc desc = GED_VIEW_FEATURE_BATCH_DESC_INIT;
    desc.owner_id = "ged-view-annotation-axes";
    desc.owner_role = "user-annotation";
    desc.overlay_class = GED_VIEW_FEATURE_OVERLAY_CLASS_USER_ANNOTATION;
    desc.lifecycle = GED_VIEW_FEATURE_LIFECYCLE_PERSISTENT;
    desc.overlay_order = GED_VIEW_FEATURE_OVERLAY_ORDER_MODEL;
    desc.local = local;
    struct ged_view_feature_batch *batch =
	ged_view_feature_batch_begin(gd->cv, &desc);
    if (!batch)
	return 0;

    struct ged_view_feature_style style = GED_VIEW_FEATURE_STYLE_INIT;
    style.visible = 1;
    style.color_valid = 1;
    for (int i = 0; i < 3; i++)
	style.color[i] = static_cast<unsigned char>(
	    axes.color[i] < 0 ? 0 : axes.color[i] > 255 ? 255 : axes.color[i]);
    style.line_width = axes.line_width;
    if (!ged_view_feature_batch_axes_replace(batch, gd->vobj,
	    &axes.position, 1, axes.size, &style)) {
	ged_view_feature_batch_abort(batch);
	return 0;
    }
    return ged_view_feature_batch_commit(batch);
}

static int
_axes_state(struct _ged_view_info *gd,
    struct ged_view_annotation_axes_state *axes,
    struct ged_view_feature_summary *summary = nullptr)
{
    struct ged_view_feature_summary current = GED_VIEW_FEATURE_SUMMARY_INIT;
    if (!ged_view_feature_get_summary(gd->cv, gd->vobj, &current) ||
	!current.exists || current.kind != GED_VIEW_FEATURE_KIND_AXES) {
	bu_vls_printf(gd->gedp->ged_result_str,
	    "View feature named %s does not exist\n", gd->vobj);
	return 0;
    }

    memset(axes, 0, sizeof(*axes));
    if (!ged_view_feature_axes_copy(gd->cv, gd->vobj, 0,
	    axes->position, &axes->size)) {
	bu_vls_printf(gd->gedp->ged_result_str,
	    "View feature %s has no axes state\n", gd->vobj);
	return 0;
    }

    struct ged_view_feature_style style = GED_VIEW_FEATURE_STYLE_INIT;
    axes->line_width = 1;
    VSET(axes->color, 255, 255, 255);
    if (ged_view_feature_style_get(gd->cv, gd->vobj, &style)) {
	if (style.color_valid)
	    VSET(axes->color, style.color[0], style.color[1], style.color[2]);
	if (style.line_width > 0)
	    axes->line_width = style.line_width;
    }
    if (summary)
	*summary = current;
    return 1;
}

static int
_axes_replace(struct _ged_view_info *gd,
    const struct ged_view_annotation_axes_state &axes)
{
    struct ged_view_annotation_axes_state current;
    struct ged_view_feature_summary summary = GED_VIEW_FEATURE_SUMMARY_INIT;
    if (!_axes_state(gd, &current, &summary))
	return 0;
    if (!_axes_publish(gd, axes,
	    summary.scope == GED_VIEW_FEATURE_SCOPE_LOCAL)) {
	bu_vls_printf(gd->gedp->ged_result_str,
	    "Failed to set axes state for %s\n", gd->vobj);
	return 0;
    }
    return 1;
}

int
_axes_cmd_create(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view annotation axes create <name> x y z";
    const char *purpose_string = "define data axes at point x,y,z";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    point_t p;
    for (int i = 0; i < 3; i++) {
	if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[i], (void *)&(p[i])) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[i]);
	    return BRLCAD_ERROR;
	}
    }

    struct ged_view_annotation_axes_state l;
    memset(&l, 0, sizeof(l));
    VMOVE(l.position, p);
    l.line_width = 1;
    l.size = 10;
    VSET(l.color, 255, 255, 0);
    if (ged_view_feature_exists(gd->cv, gd->vobj)) {
	bu_vls_printf(gedp->ged_result_str,
	    "View feature named %s already exists\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    if (!_axes_publish(gd, l, gd->local_obj)) {
	bu_vls_printf(gedp->ged_result_str,
	    "Failed to set axes state for %s\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

int
_axes_cmd_pos(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view annotation axes pos <name> [x y z]";
    const char *purpose_string = "adjust axes position";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct ged_view_annotation_axes_state a;
    if (!_axes_state(gd, &a))
	return BRLCAD_ERROR;
    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%f %f %f\n", V3ARGS(a.position));
	return BRLCAD_OK;
    }
    if (argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    point_t p;
    for (int i = 0; i < 3; i++) {
	if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[i], (void *)&(p[i])) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[i]);
	    return BRLCAD_ERROR;
	}
    }

    VMOVE(a.position, p);
    return _axes_replace(gd, a) ?
	BRLCAD_OK : BRLCAD_ERROR;
}

int
_axes_cmd_size(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view annotation axes size <name> [#]";
    const char *purpose_string = "adjust axes size";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct ged_view_annotation_axes_state a;
    if (!_axes_state(gd, &a))
	return BRLCAD_ERROR;
     if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%f\n", a.size);
	return BRLCAD_OK;
    }

    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    fastf_t val;
    if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[0], (void *)&val) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    a.size = val;
    return _axes_replace(gd, a) ?
	BRLCAD_OK : BRLCAD_ERROR;
}

int
_axes_cmd_linewidth(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view annotation axes line_width <name> [#]";
    const char *purpose_string = "adjust axes line width";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct ged_view_annotation_axes_state a;
    if (!_axes_state(gd, &a))
	return BRLCAD_ERROR;
     if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", a.line_width);
	return BRLCAD_OK;
    }

    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int val;
    if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&val) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    if (val < 1) {
	bu_vls_printf(gedp->ged_result_str, "Smallest supported value is 1\n");
	return BRLCAD_ERROR;
    }

    a.line_width = val;
    return _axes_replace(gd, a) ?
	BRLCAD_OK : BRLCAD_ERROR;
}

int
_axes_cmd_axes_color(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view annotation axes axes_color <name> [r/g/b]";
    const char *purpose_string = "get/set color of axes";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct ged_view_annotation_axes_state a;
    if (!_axes_state(gd, &a))
	return BRLCAD_ERROR;
     if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d %d %d\n", a.color[0], a.color[1], a.color[2]);
	return BRLCAD_OK;
    }

    // For color need either 1 or 3 non-subcommand args
    if (argc != 1 && argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    struct bu_color c;
    int opt_ret = bu_opt_color(NULL, argc, (const char **)argv, (void *)&c);
    if (opt_ret != 1 && opt_ret != 3) {
	bu_vls_printf(gedp->ged_result_str, "Invalid color specifier\n");
	return BRLCAD_ERROR;
    }

    bu_color_to_rgb_ints(&c, &a.color[0], &a.color[1], &a.color[2]);
    return _axes_replace(gd, a) ?
	BRLCAD_OK : BRLCAD_ERROR;
}


const struct bu_cmdtab _axes_cmds[] = {
    { "create",            _axes_cmd_create},
    { "pos",               _axes_cmd_pos},
    { "size",              _axes_cmd_size},
    { "line_width",        _axes_cmd_linewidth},
    { "axes_color",        _axes_cmd_axes_color},
    { (char *)NULL,      NULL}
};

int
_view_cmd_axes(void *bs, int argc, const char **argv)
{
    int help = 0;
    int s_version = 0;
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;

    const char *usage_string = "view annotation axes [options] [args]";
    const char *purpose_string = "create/manipulate view axes";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    if (!gd->cv) {
	bu_vls_printf(gedp->ged_result_str, ": no view specified or current");
	return BRLCAD_ERROR;
    }


    // We know we're the axes command - start processing args
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
	if (bu_cmd_valid(_axes_cmds, argv[i]) == BRLCAD_OK) {
	    cmd_pos = i;
	    break;
	}
    }

    int acnt = (cmd_pos >= 0) ? cmd_pos : argc;
    (void)bu_opt_parse(NULL, acnt, argv, d);

    return _ged_subcmd_exec(gedp, d, _axes_cmds, "view annotation axes", "[options] subcommand [args]", gd, argc, argv, help, cmd_pos);
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
