/*                         P L O T . C
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
/** @file libged/plot.cpp
 *
 * The plot command.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "bn.h"
#include "bg/plot3.h"
#include "bg/clip.h"
#include "bv.h"
#include "rt/view.h"

#include "BObol/BExportAction.h"
#include "BObol/BViewController.h"
#include <Inventor/SoViewport.h>

#include "../ged_private.h"
#include "../ged_bobol_private.hpp"

#if defined(HAVE_POPEN) && !defined(HAVE_DECL_POPEN) && !defined(popen)
extern FILE *popen(const char *command, const char *type);
#endif
#if defined(HAVE_POPEN) && !defined(HAVE_POPEN_DECL) && !defined(pclose)
extern int pclose(FILE *stream);
#endif

/* Callback data for plot operations */
struct plot_data {
    FILE *fp;
    mat_t model2view;
    int floating;
    mat_t center;
    fastf_t scale;
    int Three_D;
    int Z_clip;
    int Dashing;
    vect_t clipmin;
    vect_t clipmax;
};

static int
plot_color_component(float value)
{
    if (value <= 0.0f)
	return 0;
    if (value >= 1.0f)
	return 255;
    return static_cast<int>(value * 255.0f + 0.5f);
}

static void
plot_line(const SoBRLExportAction::LineRecord &line, struct plot_data *pd)
{
    if (!pd)
	return;

    point_t a = {line.a[0], line.a[1], line.a[2]};
    point_t b = {line.b[0], line.b[1], line.b[2]};
    const int red = plot_color_component(line.color[0]);
    const int green = plot_color_component(line.color[1]);
    const int blue = plot_color_component(line.color[2]);

    if (pd->Dashing != line.lineStyle) {
	pl_linmod(pd->fp, line.lineStyle ? "dotdashed" : "solid");
	pd->Dashing = line.lineStyle;
    }

    if (pd->floating) {
	pl_color(pd->fp, red, green, blue);
	pdv_3move(pd->fp, a);
	pdv_3cont(pd->fp, b);
	return;
    }

    vect_t start;
    vect_t fin;
    MAT4X3PNT(start, pd->model2view, a);
    MAT4X3PNT(fin, pd->model2view, b);
    if (bg_ray_vclip(start, fin, pd->clipmin, pd->clipmax) == 0)
	return;

    if (pd->Three_D) {
	pl_color(pd->fp, red, green, blue);
	pl_3line(pd->fp,
		(int)(start[X] * RT_VIEW_MAX),
		(int)(start[Y] * RT_VIEW_MAX),
		(int)(start[Z] * RT_VIEW_MAX),
		(int)(fin[X] * RT_VIEW_MAX),
		(int)(fin[Y] * RT_VIEW_MAX),
		(int)(fin[Z] * RT_VIEW_MAX));
    } else {
	pl_line(pd->fp,
		(int)(start[X] * RT_VIEW_MAX),
		(int)(start[Y] * RT_VIEW_MAX),
		(int)(fin[X] * RT_VIEW_MAX),
		(int)(fin[Y] * RT_VIEW_MAX));
    }
}

static void
plot_visible_lines(struct ged_view_context *view_ctx, struct plot_data *pd)
{
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!controller || !controller->getViewport() ||
	!controller->getViewport()->getRoot())
	return;

    SoBRLExportAction export_action;
    export_action.setGeometryPolicy(SoBRLExportAction::DISPLAY_LEVEL);
    export_action.apply(controller->getViewport()->getRoot());

    std::vector<SoBRLExportAction::ObjectRecord> records;
    export_action.collectObjectRecords(records,
	SoBRLExportAction::QUERY_VISIBLE_ONLY);
    for (const SoBRLExportAction::ObjectRecord &record : records) {
	for (int line_index : record.lineIndices)
	    plot_line(export_action.getLine(line_index), pd);
    }
}

void
dl_plot(struct ged_view_context *view_ctx, FILE *fp, mat_t model2view, int floating, mat_t center, fastf_t scale, int Three_D, int Z_clip)
{
    struct plot_data pd;

    pd.fp = fp;
    MAT_COPY(pd.model2view, model2view);
    pd.floating = floating;
    MAT_COPY(pd.center, center);
    pd.scale = scale;
    pd.Three_D = Three_D;
    pd.Z_clip = Z_clip;
    pd.Dashing = 0;

    if (floating) {
	pd_3space(fp,
		  -center[MDX] - scale,
		  -center[MDY] - scale,
		  -center[MDZ] - scale,
		  -center[MDX] + scale,
		  -center[MDY] + scale,
		  -center[MDZ] + scale);
	pl_linmod(fp, "solid");
	plot_visible_lines(view_ctx, &pd);
	return;
    }

    /*
     * Integer output version, either 2-D or 3-D.
     * Viewing region is from -1.0 to +1.0
     * which is mapped to integer space -2048 to +2048 for plotting.
     * Compute the clipping bounds of the screen in view space.
     */
    pd.clipmin[X] = -1.0;
    pd.clipmax[X] =  1.0;
    pd.clipmin[Y] = -1.0;
    pd.clipmax[Y] =  1.0;
    if (Z_clip) {
	pd.clipmin[Z] = -1.0;
	pd.clipmax[Z] =  1.0;
    } else {
	pd.clipmin[Z] = -1.0e20;
	pd.clipmax[Z] =  1.0e20;
    }

    if (Three_D)
	pl_3space(fp, (int)RT_VIEW_MIN, (int)RT_VIEW_MIN, (int)RT_VIEW_MIN, (int)RT_VIEW_MAX, (int)RT_VIEW_MAX, (int)RT_VIEW_MAX);
    else
	pl_space(fp, (int)RT_VIEW_MIN, (int)RT_VIEW_MIN, (int)RT_VIEW_MAX, (int)RT_VIEW_MAX);
    pl_erase(fp);
    pl_linmod(fp, "solid");
    plot_visible_lines(view_ctx, &pd);
}


/*
 * plot file [opts]
 * potential options might include:
 * grid, 3d w/color, |filter, infinite Z
 */
int
ged_plot_core(struct ged *gedp, int argc, const char *argv[])
{
    FILE *fp;
    int Three_D;			/* 0=2-D -vs- 1=3-D */
    int Z_clip;			/* Z clipping */
    int floating;			/* 3-D floating point plot */
    int is_pipe = 0;
    mat_t center;
    mat_t model2view;
    fastf_t scale;
    struct ged_view_context *view_ctx;
    static const char *plot_usage = "file [2|3] [f] [g] [z]";

    GED_CHECK_DRAWABLE(gedp, BRLCAD_ERROR);
    GED_CHECK_VIEW(gedp, BRLCAD_ERROR);
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], plot_usage);
	return GED_HELP;
    }

    /* Process any options */
    Three_D = 1;				/* 3-D w/color, by default */
    Z_clip = 0;				/* NO Z clipping, by default*/
    floating = 0;
    while (argv[1] != (char *)0 && argv[1][0] == '-') {
	switch (argv[1][1]) {
	    case 'f':
		floating = 1;
		break;
	    case '3':
		Three_D = 1;
		break;
	    case '2':
		Three_D = 0;		/* 2-D, for portability */
		break;
	    case 'g':
		/* do grid */
		bu_vls_printf(gedp->ged_result_str, "%s: grid unimplemented\n", argv[0]);
		break;
	    case 'z':
	    case 'Z':
		/* Enable Z clipping */
		bu_vls_printf(gedp->ged_result_str, "%s: Clipped in Z to viewing cube\n", argv[0]);
		Z_clip = 1;
		break;
	    default:
		bu_vls_printf(gedp->ged_result_str, "%s: bad PLOT option %s\n", argv[0], argv[1]);
		break;
	}
	argv++;
    }
    if (argv[1] == (char *)0) {
	bu_vls_printf(gedp->ged_result_str, "%s: no filename or filter specified\n", argv[0]);
	return BRLCAD_ERROR;
    }
    if (argv[1][0] == '|') {
	struct bu_vls str = BU_VLS_INIT_ZERO;
	bu_vls_strcpy(&str, &argv[1][1]);
	while ((++argv)[1] != (char *)0) {
	    bu_vls_strcat(&str, " ");
	    bu_vls_strcat(&str, argv[1]);
	}
	fp = popen(bu_vls_addr(&str), "wb");
	if (fp == NULL) {
	    perror(bu_vls_addr(&str));
	    return BRLCAD_ERROR;
	}

	bu_vls_printf(gedp->ged_result_str, "piped to %s\n", bu_vls_addr(&str));
	bu_vls_free(&str);
	is_pipe = 1;
    } else {
	fp = fopen(argv[1], "wb");
	if (fp == NULL) {
	    perror(argv[1]);
	    return BRLCAD_ERROR;
	}

	bu_vls_printf(gedp->ged_result_str, "plot stored in %s\n", argv[1]);
	is_pipe = 0;
    }

    view_ctx = ged_view_active_ctx(gedp);
    const struct bv *view = bv_context_view_const((const struct bv_context *)view_ctx);
    bv_model2view_get(model2view, view);
    bv_center_mat_get(center, view);
    scale = bv_scale_get(view);
    dl_plot(view_ctx, fp, model2view, floating, center, scale, Three_D, Z_clip);

    if (is_pipe)
	(void)pclose(fp);
    else
	(void)fclose(fp);

    return BRLCAD_OK;
}


#include "../include/plugin.h"

#define GED_PLOT_COMMANDS(X, XID) \
    X(plot, ged_plot_core, GED_CMD_DEFAULT) \

GED_DECLARE_COMMAND_SET(GED_PLOT_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_plot", 1, GED_PLOT_COMMANDS)

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
