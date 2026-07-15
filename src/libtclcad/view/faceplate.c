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
#include "bv.h"
#include "bu/units.h"
#include "brlobol/display_endpoint.h"
#include "ged.h"
#include "ged/view.h"
#include "rt/view.h"
#include "tclcad.h"

/* Private headers */
#include "../tclcad_private.h"
#include "../view/view.h"

static int
tclcad_faceplate_color_endpoint_set(void *view_ctx, const char *property_name,
	int r, int g, int b)
{
    if (!view_ctx || !property_name || r < 0 || r > 255 || b < 0 ||
	b > 255 || g < 0 || g > 255)
	return BRLOBOL_ENDPOINT_PROPERTY_INVALID;

    struct brlobol_endpoint_property_value value =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    value.type = BRLOBOL_ENDPOINT_PROPERTY_COLOR3;
    value.color3[0] = r / 255.0;
    value.color3[1] = g / 255.0;
    value.color3[2] = b / 255.0;
    return ged_view_context_display_property_set(view_ctx, property_name,
	&value);
}

static int
tclcad_faceplate_bool_endpoint_set(void *view_ctx, const char *property_name,
	int enabled)
{
    if (!view_ctx || !property_name || enabled < 0 || enabled > 1)
	return BRLOBOL_ENDPOINT_PROPERTY_INVALID;

    struct brlobol_endpoint_property_value value =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    value.type = BRLOBOL_ENDPOINT_PROPERTY_BOOL;
    value.bool_value = enabled;
    return ged_view_context_display_property_set(view_ctx, property_name,
	&value);
}

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
    struct bv *view = NULL;

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
    view = bv_context_view((struct bv_context *)view_ctx);
    if (!view)
	goto bad;

    if (BU_STR_EQUAL(argv[2], "center_dot")) {
	struct bv_other_state center_dot;
	if (!bv_center_dot_state_get(&center_dot, view))
	    goto bad;
	if (BU_STR_EQUAL(argv[3], "draw")) {
	    if (argc == 4) {
		bu_vls_printf(gedp->ged_result_str, "%d", center_dot.gos_draw);
		return BRLCAD_OK;
	    } else if (argc == 5) {
		if (bu_sscanf(argv[4], "%d", &i) != 1)
		    goto bad;

		if (ged_view_context_display_endpoint_get(view_ctx)) {
		    if (tclcad_faceplate_bool_endpoint_set(view_ctx,
			"view.faceplate.center_dot.visible", i) !=
			BRLOBOL_ENDPOINT_PROPERTY_OK)
			goto bad;
		} else {
		    center_dot.gos_draw = i ? 1 : 0;
		    bv_center_dot_state_set(view, &center_dot);
		}
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

		if (ged_view_context_display_endpoint_get(view_ctx)) {
		    if (tclcad_faceplate_color_endpoint_set(view_ctx,
			"view.faceplate.center_dot.color", r, g, b) !=
		BRLOBOL_ENDPOINT_PROPERTY_OK)
			goto bad;
		} else {
		    if (r < 0 || r > 255 || g < 0 || g > 255 || b < 0 || b > 255)
			goto bad;
		    VSET(center_dot.gos_line_color, r, g, b);
		    bv_center_dot_state_set(view, &center_dot);
		}
		to_refresh_view(view_ctx);
		return BRLCAD_OK;
	    }
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "prim_labels")) {
	struct bv_other_state prim_labels;
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
	struct bv_params_state params;
	if (!bv_params_state_get(&params, view))
	    goto bad;
	if (BU_STR_EQUAL(argv[3], "draw")) {
	    if (argc == 4) {
		bu_vls_printf(gedp->ged_result_str, "%d", params.draw);
		return BRLCAD_OK;
	    } else if (argc == 5) {
		if (bu_sscanf(argv[4], "%d", &i) != 1)
		    goto bad;

		if (ged_view_context_display_endpoint_get(view_ctx)) {
		    if (tclcad_faceplate_bool_endpoint_set(view_ctx,
			"view.faceplate.params.visible", i) !=
			BRLOBOL_ENDPOINT_PROPERTY_OK)
			goto bad;
		} else {
		    params.draw = i ? 1 : 0;
		    bv_params_state_set(view, &params);
		}
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

		if (ged_view_context_display_endpoint_get(view_ctx)) {
		    if (tclcad_faceplate_color_endpoint_set(view_ctx,
			"view.faceplate.params.color", r, g, b) !=
		BRLOBOL_ENDPOINT_PROPERTY_OK)
			goto bad;
		} else {
		    if (r < 0 || r > 255 || g < 0 || g > 255 || b < 0 || b > 255)
			goto bad;
		    VSET(params.color, r, g, b);
		    bv_params_state_set(view, &params);
		}
		to_refresh_view(view_ctx);
		return BRLCAD_OK;
	    }
	}

	goto bad;
    }

    if (BU_STR_EQUAL(argv[2], "view_scale")) {
	struct bv_other_state scale_state;
	if (!bv_scale_overlay_state_get(&scale_state, view))
	    goto bad;
	if (BU_STR_EQUAL(argv[3], "draw")) {
	    if (argc == 4) {
		bu_vls_printf(gedp->ged_result_str, "%d", scale_state.gos_draw);
		return BRLCAD_OK;
	    } else if (argc == 5) {
		if (bu_sscanf(argv[4], "%d", &i) != 1)
		    goto bad;

		if (ged_view_context_display_endpoint_get(view_ctx)) {
		    if (tclcad_faceplate_bool_endpoint_set(view_ctx,
			"view.faceplate.scale.visible", i) !=
			BRLOBOL_ENDPOINT_PROPERTY_OK)
			goto bad;
		} else {
		    scale_state.gos_draw = i ? 1 : 0;
		    bv_scale_overlay_state_set(view, &scale_state);
		}
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

		if (ged_view_context_display_endpoint_get(view_ctx)) {
		    if (tclcad_faceplate_color_endpoint_set(view_ctx,
			"view.faceplate.scale.color", r, g, b) !=
		BRLOBOL_ENDPOINT_PROPERTY_OK)
			goto bad;
		} else {
		    if (r < 0 || r > 255 || g < 0 || g > 255 || b < 0 || b > 255)
			goto bad;
		    VSET(scale_state.gos_line_color, r, g, b);
		    bv_scale_overlay_state_set(view, &scale_state);
		}
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
