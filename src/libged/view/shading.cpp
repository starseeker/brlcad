/*                    S H A D I N G . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file libged/view/shading.cpp
 *
 * Per-view mesh normal presentation policy.  This policy changes only the
 * renderer representation; retained PoP topology and resident prefixes remain
 * reusable.
 */

#include "common.h"

#include <string.h>

#include "bu/opt.h"
#include "bu/vls.h"
#include "bv.h"
#include "ged/draw_obol.h"

#include "../ged_bobol_private.hpp"
#include "../ged_private.h"
#include "./ged_view.h"

static const char *
normal_style_name(int style)
{
    switch (style) {
	case BV_NORMAL_FLAT:
	    return "flat";
	case BV_NORMAL_SMOOTH:
	    return "smooth";
	default:
	    return "authored";
    }
}

static int
normal_style_parse(const char *value, int *style)
{
    if (!value || !style)
	return 0;
    if (BU_STR_EQUAL(value, "authored") || BU_STR_EQUAL(value, "auto")) {
	*style = BV_NORMAL_AUTHORED;
	return 1;
    }
    if (BU_STR_EQUAL(value, "flat")) {
	*style = BV_NORMAL_FLAT;
	return 1;
    }
    if (BU_STR_EQUAL(value, "smooth")) {
	*style = BV_NORMAL_SMOOTH;
	return 1;
    }
    return 0;
}

static int
shading_commit(struct ged *gedp, struct ged_view_context *view_ctx,
	struct bv *view, const struct bv_shading_state *shading)
{
    if (!bv_shading_state_set(view, shading))
	return BRLCAD_ERROR;
    (void)ged_draw_obol_shading_sync(gedp, view_ctx);
    (void)bv_refresh_request(view, GED_VIEW_REFRESH_VIEW);
    ged_refresh_cb(gedp);
    return BRLCAD_OK;
}

int
ged_shading_core(struct ged *gedp, int argc, const char *argv[])
{
    if (!gedp || !argc || !argv)
	return BRLCAD_ERROR;
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* consume the "shading" command word */
    argc--;
    argv++;

    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    struct bv *view = view_ctx ?
	bv_context_view((struct bv_context *)view_ctx) : NULL;
    struct bv_shading_state shading;
    if (!view || !bv_shading_state_get(&shading, view)) {
	bu_vls_printf(gedp->ged_result_str, "no current view set");
	return BRLCAD_ERROR;
    }

    if (!argc) {
	bu_vls_printf(gedp->ged_result_str, "normals %s\ncrease %g",
	    normal_style_name(shading.normal_style),
	    shading.normal_crease_angle);
	BObolViewController *controller =
	    ged_bobol_view_controller(view_ctx);
	if (controller) {
	    bu_vls_printf(gedp->ged_result_str,
		"\nrenderer_normals %s\nrenderer_crease %g",
		normal_style_name(
		    static_cast<int>(controller->getNormalStyle())),
		controller->getNormalCreaseAngle());
	}
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(argv[0], "normals")) {
	if (argc == 1) {
	    bu_vls_printf(gedp->ged_result_str, "%s",
		normal_style_name(shading.normal_style));
	    return BRLCAD_OK;
	}
	if (argc != 2 ||
	    !normal_style_parse(argv[1], &shading.normal_style)) {
	    bu_vls_printf(gedp->ged_result_str,
		"usage: view shading normals authored|flat|smooth");
	    return BRLCAD_ERROR;
	}
	return shading_commit(gedp, view_ctx, view, &shading);
    }

    if (BU_STR_EQUAL(argv[0], "crease")) {
	if (argc == 1) {
	    bu_vls_printf(gedp->ged_result_str, "%g",
		shading.normal_crease_angle);
	    return BRLCAD_OK;
	}
	fastf_t angle = 0.0;
	if (argc != 2 ||
	    bu_opt_fastf_t(NULL, 1, &argv[1], (void *)&angle) != 1 ||
	    angle < 0.0 || angle > 180.0) {
	    bu_vls_printf(gedp->ged_result_str,
		"crease angle must be between 0 and 180 degrees");
	    return BRLCAD_ERROR;
	}
	shading.normal_crease_angle = angle;
	return shading_commit(gedp, view_ctx, view, &shading);
    }

    bu_vls_printf(gedp->ged_result_str,
	"usage: view shading [normals authored|flat|smooth] [crease degrees]");
    return BRLCAD_ERROR;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
