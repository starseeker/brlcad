/*                      L A B E L S . C P P
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
/** @file libged/view/labels.cpp
 *
 * Commands for view labels.
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
#include "bv.h"

#include "../ged_private.h"
#include "../ged_bobol_private.hpp"
#include "./ged_view.h"

int
_label_cmd_create(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view annotation label create <name> text x y [z] [px py pz]";
    const char *purpose_string = "start a label at point x,y,[z], possibly targeting point px,py,pz";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (argc != 3 && argc != 4 && argc != 6 && argc != 7) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    point_t p;
    if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[1], (void *)&(p[0])) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[1]);
	return BRLCAD_ERROR;
    }
    if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[2], (void *)&(p[1])) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[2]);
	return BRLCAD_ERROR;
    }

    if (argc == 4 || argc == 7) {
	if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[3], (void *)&(p[2])) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[3]);
	    return BRLCAD_ERROR;
	}
    } else {
	fastf_t fx, fy;
	const struct bv *view =
	    bv_context_view_const(ged_view_context_bv_const(gd->cv));
	if (!bv_screen_to_view(&fx, &fy, view,
		(int)p[0], (int)p[1])) {
	    return BRLCAD_ERROR;
	}
	p[0] = fx;
	p[1] = fy;
	p[2] = 0;
	point_t tp;
	VMOVE(tp, p);
	mat_t view2model;
	bv_view2model_get(view2model, view);
	MAT4X3PNT(p, view2model, tp);
    }
    point_t target;
    if (argc == 6) {
	if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[3], (void *)&(target[0])) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[3]);
	    return BRLCAD_ERROR;
	}
	if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[4], (void *)&(target[1])) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[4]);
	    return BRLCAD_ERROR;
	}
	if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[5], (void *)&(target[2])) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[5]);
	    return BRLCAD_ERROR;
	}
    }

    if (argc == 7) {
	if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[4], (void *)&(target[0])) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[4]);
	    return BRLCAD_ERROR;
	}
	if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[5], (void *)&(target[1])) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[5]);
	    return BRLCAD_ERROR;
	}
	if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[6], (void *)&(target[2])) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[6]);
	    return BRLCAD_ERROR;
	}
    }

    if (ged_bobol_feature_find(gd->gedp, gd->cv, gd->vobj).handle.isValid()) {
	bu_vls_printf(gedp->ged_result_str,
	    "View feature named %s already exists\n", gd->vobj);
	return BRLCAD_ERROR;
    }

    BObolViewController *controller = ged_bobol_feature_controller(
	gd->gedp, gd->cv, gd->local_obj);
    if (!controller) {
	bu_vls_printf(gedp->ged_result_str, "Failed to create %s\n",
	    gd->vobj);
	return BRLCAD_ERROR;
    }

    BObolLabel label;
    label.text = argv[0];
    label.point = SbVec3f(static_cast<float>(p[X]),
	static_cast<float>(p[Y]), static_cast<float>(p[Z]));
    if (argc == 6 || argc == 7) {
	label.hasLeader = TRUE;
	label.target = SbVec3f(static_cast<float>(target[X]),
	    static_cast<float>(target[Y]), static_cast<float>(target[Z]));
    }
    std::vector<BObolLabel> labels(1, label);
    BObolFeatureStyle style;
    style.hasColor = TRUE;
    style.color = SbColor(1.0f, 1.0f, 0.0f);
    BObolFeatureOwner owner = ged_bobol_view_feature_owner(gd->cv);
    const BObolFeatureScope scope = gd->local_obj ?
	BObolFeatureScope::Local : BObolFeatureScope::Shared;
    const BObolFeatureOwner *owner_ptr = gd->local_obj ? &owner : nullptr;
    if (!controller->features().publishLabels(gd->vobj, scope, labels,
	    &style, owner_ptr).isValid()) {
	bu_vls_printf(gedp->ged_result_str, "Failed to create %s\n",
	    gd->vobj);
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

const struct bu_cmdtab _label_cmds[] = {
    { "create",          _label_cmd_create},
    { (char *)NULL,      NULL}
};

int
_view_cmd_labels(void *bs, int argc, const char **argv)
{
    int help = 0;
    int s_version = 0;
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;

    const char *usage_string = "view annotation label [options] [args]";
    const char *purpose_string = "create/manipulate view labels";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    if (!gd->cv) {
	bu_vls_printf(gedp->ged_result_str, ": no view specified or current in GED");
	return BRLCAD_ERROR;
    }


    // We know we're the labels command - start processing args
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
	if (bu_cmd_valid(_label_cmds, argv[i]) == BRLCAD_OK) {
	    cmd_pos = i;
	    break;
	}
    }

    int acnt = (cmd_pos >= 0) ? cmd_pos : argc;
    (void)bu_opt_parse(NULL, acnt, argv, d);

    return _ged_subcmd_exec(gedp, d, _label_cmds, "view annotation label", "[options] subcommand [args]", gd, argc, argv, help, cmd_pos);
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
