/*                     F A C E P L A T E . C
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
/** @file libtclcad/faceplate.c
 *
 */
/** @} */

#include "common.h"
#include "bu/units.h"
#include "ged.h"
#include "rt/view_legacy_bsg.h"
#include "tclcad.h"

/* Private headers */
#include "../tclcad_private.h"
#include "../view/view.h"

int
to_faceplate(struct ged *gedp,
	     int argc,
	     const char *argv[],
	     ged_func_ptr UNUSED(func),
	     const char *usage,
	     int UNUSED(maxargs))
{
    int i;
    void *view_ctx;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 4 || 7 < argc)
	goto bad;

    view_ctx = ged_view_find_ctx(gedp, argv[1]);
    if (!view_ctx) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(argv[2], "center_dot")) {
	struct rt_view_other_state center_dot;
	if (!rt_view_context_center_dot_state_from_bsg(&center_dot, view_ctx))
	    goto bad;
	if (BU_STR_EQUAL(argv[3], "draw")) {
	    if (argc == 4) {
		bu_vls_printf(gedp->ged_result_str, "%d", center_dot.gos_draw);
		return BRLCAD_OK;
	    } else if (argc == 5) {
		if (bu_sscanf(argv[4], "%d", &i) != 1)
		    goto bad;

		if (i)
		    center_dot.gos_draw = 1;
		else
		    center_dot.gos_draw = 0;

		rt_view_context_center_dot_state_set_bsg(view_ctx, &center_dot);
		to_refresh_view(view_ctx);
		return BRLCAD_OK;
	    }
	}

	if (BU_STR_EQUAL(argv[3], "color")) {
	    if (argc == 4) {
		bu_vls_printf(gedp->ged_result_str, "%d %d %d", V3ARGS(center_dot.gos_line_color));
		return BRLCAD_OK;
	    } else if (argc == 7) {
		int r, g, b;

		if (bu_sscanf(argv[4], "%d", &r) != 1 ||
		    bu_sscanf(argv[5], "%d", &g) != 1 ||
		    bu_sscanf(argv[6], "%d", &b) != 1)
		    goto bad;

		VSET(center_dot.gos_line_color, r, g, b);
		rt_view_context_center_dot_state_set_bsg(view_ctx, &center_dot);
		to_refresh_view(view_ctx);
		return BRLCAD_OK;
	    }
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "prim_labels")) {
	struct rt_view_other_state prim_labels;
	if (!tclcad_view_prim_labels_state_from_view_ctx(&prim_labels, view_ctx))
	    goto bad;

	if (BU_STR_EQUAL(argv[3], "draw")) {
	    if (argc == 4) {
		bu_vls_printf(gedp->ged_result_str, "%d", prim_labels.gos_draw);
		return BRLCAD_OK;
	    } else if (argc == 5) {
		if (bu_sscanf(argv[4], "%d", &i) != 1)
		    goto bad;

		if (i)
		    prim_labels.gos_draw = 1;
		else
		    prim_labels.gos_draw = 0;

		if (!tclcad_view_prim_labels_state_set(view_ctx, &prim_labels))
		    goto bad;
		to_refresh_view(view_ctx);
		return BRLCAD_OK;
	    }
	}

	if (BU_STR_EQUAL(argv[3], "color")) {
	    if (argc == 4) {
		bu_vls_printf(gedp->ged_result_str, "%d %d %d", V3ARGS(prim_labels.gos_text_color));
		return BRLCAD_OK;
	    } else if (argc == 7) {
		int r, g, b;

		if (bu_sscanf(argv[4], "%d", &r) != 1 ||
		    bu_sscanf(argv[5], "%d", &g) != 1 ||
		    bu_sscanf(argv[6], "%d", &b) != 1)
		    goto bad;

		VSET(prim_labels.gos_text_color, r, g, b);
		if (!tclcad_view_prim_labels_state_set(view_ctx, &prim_labels))
		    goto bad;
		to_refresh_view(view_ctx);
		return BRLCAD_OK;
	    }
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "view_params")) {
	struct rt_view_params_state params;
	if (!rt_view_context_params_state_from_bsg(&params, view_ctx))
	    goto bad;
	if (BU_STR_EQUAL(argv[3], "draw")) {
	    if (argc == 4) {
		bu_vls_printf(gedp->ged_result_str, "%d", params.draw);
		return BRLCAD_OK;
	    } else if (argc == 5) {
		if (bu_sscanf(argv[4], "%d", &i) != 1)
		    goto bad;

		if (i)
		    params.draw = 1;
		else
		    params.draw = 0;

		rt_view_context_params_state_set_bsg(view_ctx, &params);
		to_refresh_view(view_ctx);
		return BRLCAD_OK;
	    }
	}

	if (BU_STR_EQUAL(argv[3], "color")) {
	    if (argc == 4) {
		bu_vls_printf(gedp->ged_result_str, "%d %d %d", V3ARGS(params.color));
		return BRLCAD_OK;
	    } else if (argc == 7) {
		int r, g, b;

		if (bu_sscanf(argv[4], "%d", &r) != 1 ||
		    bu_sscanf(argv[5], "%d", &g) != 1 ||
		    bu_sscanf(argv[6], "%d", &b) != 1)
		    goto bad;

		VSET(params.color, r, g, b);
		rt_view_context_params_state_set_bsg(view_ctx, &params);
		to_refresh_view(view_ctx);
		return BRLCAD_OK;
	    }
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "view_scale")) {
	struct rt_view_other_state scale_state;
	if (!rt_view_context_scale_overlay_state_from_bsg(&scale_state, view_ctx))
	    goto bad;
	if (BU_STR_EQUAL(argv[3], "draw")) {
	    if (argc == 4) {
		bu_vls_printf(gedp->ged_result_str, "%d", scale_state.gos_draw);
		return BRLCAD_OK;
	    } else if (argc == 5) {
		if (bu_sscanf(argv[4], "%d", &i) != 1)
		    goto bad;

		if (i)
		    scale_state.gos_draw = 1;
		else
		    scale_state.gos_draw = 0;

		rt_view_context_scale_overlay_state_set_bsg(view_ctx, &scale_state);
		to_refresh_view(view_ctx);
		return BRLCAD_OK;
	    }
	}

	if (BU_STR_EQUAL(argv[3], "color")) {
	    if (argc == 4) {
		bu_vls_printf(gedp->ged_result_str, "%d %d %d", V3ARGS(scale_state.gos_line_color));
		return BRLCAD_OK;
	    } else if (argc == 7) {
		int r, g, b;

		if (bu_sscanf(argv[4], "%d", &r) != 1 ||
		    bu_sscanf(argv[5], "%d", &g) != 1 ||
		    bu_sscanf(argv[6], "%d", &b) != 1)
		    goto bad;

		VSET(scale_state.gos_line_color, r, g, b);
		rt_view_context_scale_overlay_state_set_bsg(view_ctx, &scale_state);
		to_refresh_view(view_ctx);
		return BRLCAD_OK;
	    }
	}

	goto bad;
    }

bad:
    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
    return BRLCAD_ERROR;
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
