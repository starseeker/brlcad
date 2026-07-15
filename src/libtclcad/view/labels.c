/*                       L A B E L S . C
 * BRL-CAD
 *
 * Copyright (c) 2000-2026 United States Government as represented by
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
/** @addtogroup libtclcad */
/** @{ */
/** @file libtclcad/view/labels.c
 *
 */
/** @} */

#include "common.h"
#include "bu/units.h"
#include "ged.h"
#include "tclcad.h"

/* Private headers */
#include "ged/draw.h"
#include "../tclcad_private.h"
#include "../view/view.h"

/* Keep draw-view label features in sync with the TclCAD per-view data-labels
 * state so the current renderer draws labels through retained features.
 *
 * The draw, color and labels getters in to_data_labels_func recover values
 * through typed GED data-label facades instead of TclCAD per-view data
 * directly.  The size getter still uses TclCAD per-view data because font size
 * is not stored in the retained draw-view child features yet.
 *
 * TclCAD per-view data continues to be written by setters here because label
 * parsing still builds a TclCAD label-state input record.  This is now
 * internal state used solely for feature sync. */

static int
_tclcad_data_labels_sync_draw_view(void *view_ctx,
				   tclcad_label_state *gdlsp,
				   const char *name)
{
    if (!view_ctx || !gdlsp || !name)
	return 0;

    if (!gdlsp->gdls_draw || gdlsp->gdls_num_labels < 1)
	return ged_draw_view_context_data_labels_replace(view_ctx, name, 0,
		NULL, 0);

    size_t label_count = (size_t)gdlsp->gdls_num_labels;
    struct ged_draw_view_label_data *labels =
	(struct ged_draw_view_label_data *)bu_calloc(label_count,
		sizeof(struct ged_draw_view_label_data), "TclCAD data labels");

    for (size_t i = 0; i < label_count; i++) {
	struct ged_draw_view_label_data init = GED_DRAW_VIEW_LABEL_DATA_INIT;
	labels[i] = init;
	labels[i].text = gdlsp->gdls_labels[i];
	VMOVE(labels[i].point, gdlsp->gdls_points[i]);
	labels[i].color_valid = 1;
	labels[i].color[0] = (unsigned char)gdlsp->gdls_color[0];
	labels[i].color[1] = (unsigned char)gdlsp->gdls_color[1];
	labels[i].color[2] = (unsigned char)gdlsp->gdls_color[2];
	labels[i].font_size = gdlsp->gdls_size;
    }

    int ret = ged_draw_view_context_data_labels_replace(view_ctx, name, 1,
	    labels, label_count);
    bu_free(labels, "TclCAD data labels");
    return ret;
}

int
go_data_labels(Tcl_Interp *interp,
	       struct ged *gedp,
	       void *draw_view_ctx,
	       int argc,
	       const char *argv[],
	       const char *usage)
{
    int ret;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 2 || 5 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }
    to_refresh_suppress_all_begin(current_top);

    ret = to_data_labels_func(interp, gedp, draw_view_ctx, argc, argv);
    to_refresh_suppress_all_end(current_top);
    if (ret & BRLCAD_ERROR)
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);

    return ret;
}


int
to_data_labels(struct ged *gedp,
	       int argc,
	       const char *argv[],
	       ged_func_ptr UNUSED(func),
	       const char *usage,
	       int UNUSED(maxargs))
{
    void *view_ctx;
    int ret;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 3 || 6 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    view_ctx = ged_view_find_ctx(gedp, argv[1]);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* shift the command name to argv[1] before calling to_data_labels_func */
    argv[1] = argv[0];
    ret = to_data_labels_func(current_top->to_interp, gedp, view_ctx, argc-1, argv+1);
    if (ret == BRLCAD_ERROR)
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);

    return ret;
}


int
to_data_labels_func(Tcl_Interp *interp,
		    struct ged *gedp,
		    void *view_ctx,
		    int argc,
		    const char *argv[])
{
    tclcad_label_state *gdlsp =
	tclcad_view_label_state_from_view_ctx(view_ctx, argv[0][0] == 's');
    if (!gdlsp)
	return BRLCAD_ERROR;

    const char *feature_name = (argv[0][0] == 's') ? "_tcl_sdata_labels" : "_tcl_data_labels";

    if (BU_STR_EQUAL(argv[1], "draw")) {
	if (argc == 2) {
	    bu_vls_printf(gedp->ged_result_str, "%d",
			  ged_draw_view_context_data_labels_draw_get(view_ctx, feature_name));
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int i;

	    if (bu_sscanf(argv[2], "%d", &i) != 1)
		goto bad;

	    if (i)
		gdlsp->gdls_draw = 1;
	    else
		gdlsp->gdls_draw = 0;

	    (void)_tclcad_data_labels_sync_draw_view(view_ctx, gdlsp, feature_name);
	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "color")) {
	if (argc == 2) {
	    unsigned char rgb[3] = {0, 0, 0};
	    if (ged_draw_view_context_data_labels_color_get(view_ctx, feature_name, rgb)) {
		bu_vls_printf(gedp->ged_result_str, "%d %d %d",
			      (int)rgb[0], (int)rgb[1], (int)rgb[2]);
	    } else {
		bu_vls_printf(gedp->ged_result_str, "0 0 0");
	    }
	    return BRLCAD_OK;
	}

	if (argc == 5) {
	    int r, g, b;

	    /* set background color */
	    if (bu_sscanf(argv[2], "%d", &r) != 1 ||
		bu_sscanf(argv[3], "%d", &g) != 1 ||
		bu_sscanf(argv[4], "%d", &b) != 1)
		goto bad;

	    /* validate color */
	    if (r < 0 || 255 < r ||
		g < 0 || 255 < g ||
		b < 0 || 255 < b)
		goto bad;

	    VSET(gdlsp->gdls_color, r, g, b);

	    (void)_tclcad_data_labels_sync_draw_view(view_ctx, gdlsp, feature_name);
	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[1], "labels")) {
	register int i;

	/* { {{label this} {0 0 0}} {{label that} {100 100 100}} }*/

	if (argc == 2) {
	    size_t _child_cnt = ged_draw_view_context_data_labels_count(view_ctx, feature_name);
	    if (_child_cnt > 0) {
		for (size_t _k = 0; _k < _child_cnt; _k++) {
		    struct bu_vls text = BU_VLS_INIT_ZERO;
		    point_t pt;
		    if (!ged_draw_view_context_data_labels_copy(view_ctx, feature_name, _k, &text, pt, NULL)) {
			bu_vls_free(&text);
			continue;
		    }
		    bu_vls_printf(gedp->ged_result_str, "{{%s}", bu_vls_cstr(&text));
		    bu_vls_printf(gedp->ged_result_str, " {%lf %lf %lf}} ", V3ARGS(pt));
		    bu_vls_free(&text);
		}
	    }
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int ac;
	    const char **av;

	    if (Tcl_SplitList(interp, argv[2], &ac, &av) != TCL_OK) {
		bu_vls_printf(gedp->ged_result_str, "%s", Tcl_GetStringResult(interp));
		return BRLCAD_ERROR;
	    }

	    if (gdlsp->gdls_num_labels) {
		bu_argv_free(gdlsp->gdls_num_labels, gdlsp->gdls_labels);
		bu_free((void *)gdlsp->gdls_points, "data points");
		gdlsp->gdls_labels = (char **)0;
		gdlsp->gdls_points = (point_t *)0;
		gdlsp->gdls_num_labels = 0;
	    }

	    /* Clear out data points */
	    if (ac < 1) {
		Tcl_Free((char *)av);
		(void)_tclcad_data_labels_sync_draw_view(view_ctx, gdlsp, feature_name);
		to_refresh_view(view_ctx);
		return BRLCAD_OK;
	    }

	    gdlsp->gdls_num_labels = ac;
	    gdlsp->gdls_labels = (char **)bu_calloc(ac, sizeof(char *), "data labels");
	    gdlsp->gdls_points = (point_t *)bu_calloc(ac, sizeof(point_t), "data points");
	    for (i = 0; i < ac; ++i) {
		int sub_ac;
		const char **sub_av;
		double scan[ELEMENTS_PER_VECT];

		if (Tcl_SplitList(interp, av[i], &sub_ac, &sub_av) != TCL_OK) {
		    /*XXX Need a macro for the following lines. Do something similar for the rest. */
		    bu_free((void *)gdlsp->gdls_labels, "data labels");
		    bu_free((void *)gdlsp->gdls_points, "data points");
		    gdlsp->gdls_labels = (char **)0;
		    gdlsp->gdls_points = (point_t *)0;
		    gdlsp->gdls_num_labels = 0;

		    bu_vls_printf(gedp->ged_result_str, "%s", Tcl_GetStringResult(interp));
		    Tcl_Free((char *)av);
		    to_refresh_view(view_ctx);
		    return BRLCAD_ERROR;
		}

		if (sub_ac != 2) {
		    /*XXX Need a macro for the following lines. Do something similar for the rest. */
		    bu_free((void *)gdlsp->gdls_labels, "data labels");
		    bu_free((void *)gdlsp->gdls_points, "data points");
		    gdlsp->gdls_labels = (char **)0;
		    gdlsp->gdls_points = (point_t *)0;
		    gdlsp->gdls_num_labels = 0;

		    bu_vls_printf(gedp->ged_result_str, "Each list element must contain a label and a point (i.e. {{some label} {0 0 0}})");
		    Tcl_Free((char *)sub_av);
		    Tcl_Free((char *)av);
		    to_refresh_view(view_ctx);
		    return BRLCAD_ERROR;
		}

		if (bu_sscanf(sub_av[1], "%lf %lf %lf", &scan[X], &scan[Y], &scan[Z]) != 3) {
		    bu_vls_printf(gedp->ged_result_str, "bad data point - %s\n", sub_av[1]);

		    /*XXX Need a macro for the following lines. Do something similar for the rest. */
		    bu_free((void *)gdlsp->gdls_labels, "data labels");
		    bu_free((void *)gdlsp->gdls_points, "data points");
		    gdlsp->gdls_labels = (char **)0;
		    gdlsp->gdls_points = (point_t *)0;
		    gdlsp->gdls_num_labels = 0;

		    Tcl_Free((char *)sub_av);
		    Tcl_Free((char *)av);
		    to_refresh_view(view_ctx);
		    return BRLCAD_ERROR;
		}
		/* convert double to fastf_t */
		VMOVE(gdlsp->gdls_points[i], scan);

		gdlsp->gdls_labels[i] = bu_strdup(sub_av[0]);
		Tcl_Free((char *)sub_av);
	    }

	    Tcl_Free((char *)av);
	    (void)_tclcad_data_labels_sync_draw_view(view_ctx, gdlsp, feature_name);
	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}
    }

    if (BU_STR_EQUAL(argv[1], "size")) {
	if (argc == 2) {
	    bu_vls_printf(gedp->ged_result_str, "%d", gdlsp->gdls_size);
	    return BRLCAD_OK;
	}

	if (argc == 3) {
	    int size;

	    if (bu_sscanf(argv[2], "%d", &size) != 1)
		goto bad;

	    gdlsp->gdls_size = size;

	    (void)_tclcad_data_labels_sync_draw_view(view_ctx, gdlsp, feature_name);
	    to_refresh_view(view_ctx);
	    return BRLCAD_OK;
	}

	goto bad;
    }


bad:
    return BRLCAD_ERROR;
}

int
to_prim_label(struct ged *gedp,
	      int argc,
	      const char *argv[],
	      ged_func_ptr UNUSED(func),
	      const char *UNUSED(usage),
	      int UNUSED(maxargs))
{
    register int i;
    struct tclcad_ged_data *tgd = (struct tclcad_ged_data *)current_top->to_gedp->u_data;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* Free the previous list of primitives scheduled for labeling */
    if (tgd->go_prim_label_list_size) {
	for (i = 0; i < tgd->go_prim_label_list_size; ++i)
	    bu_vls_free(&tgd->go_prim_label_list[i]);
	bu_free(tgd->go_prim_label_list, "prim_label");
	tgd->go_prim_label_list = NULL;
	tgd->go_prim_label_list_size = 0;
    }

    /* Set the list of primitives scheduled for labeling */
    tgd->go_prim_label_list_size = argc - 1;
    if (tgd->go_prim_label_list_size < 1)
	return BRLCAD_OK;

    tgd->go_prim_label_list = (struct bu_vls *)bu_calloc(tgd->go_prim_label_list_size,
									  sizeof(struct bu_vls), "prim_label");
    for (i = 0; i < tgd->go_prim_label_list_size; ++i) {
	bu_vls_init(&tgd->go_prim_label_list[i]);
	bu_vls_printf(&tgd->go_prim_label_list[i], "%s", argv[i+1]);
    }

    return BRLCAD_OK;
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
