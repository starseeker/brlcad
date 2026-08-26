/*                       C U T T I N G . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libged/view/cutting.cpp
 *
 * World-space cutting-plane controls.  This is intentionally distinct from
 * view.zclip: navigation changes camera clipping, never this model-space
 * section plane.
 */

#include "common.h"

#include <cerrno>
#include <cmath>
#include <cstdlib>

#include "bu/str.h"
#include "bu/vls.h"

#include "BObol/BViewController.h"
#include "ged/draw.h"

#include "../ged_bobol_private.hpp"
#include "../ged_private.h"
#include "./ged_view.h"

static int
cutting_number(struct ged *gedp, const char *text, double *value)
{
    if (!text || !value)
	return 0;
    char *end = NULL;
    errno = 0;
    const double parsed = strtod(text, &end);
    if (errno || end == text || *end || !std::isfinite(parsed)) {
	bu_vls_printf(gedp->ged_result_str,
	    "invalid finite coordinate: %s\n", text);
	return 0;
    }
    *value = parsed;
    return 1;
}

static int
cutting_vector(struct ged *gedp, const char *const *argv, SbVec3f *value)
{
    double coordinates[3] = {};
    for (size_t i = 0; i < 3; i++) {
	if (!cutting_number(gedp, argv[i], &coordinates[i]))
	    return 0;
    }
    *value = SbVec3f(static_cast<float>(coordinates[0]),
	static_cast<float>(coordinates[1]), static_cast<float>(coordinates[2]));
    return 1;
}

static SbVec3f
cutting_plane_point(const SbPlane &plane)
{
    return plane.getNormal() * plane.getDistanceFromOrigin();
}

int
ged_cutting_core(struct ged *gedp, int argc, const char *argv[])
{
    const char *usage =
	"view cutting [enable [0|1]|origin x y z|normal x y z]";
    if (!gedp || !argv || argc < 1)
	return BRLCAD_ERROR;
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!controller) {
	bu_vls_printf(gedp->ged_result_str,
	    "active view has no Obol cutting-plane controller\n");
	return BRLCAD_ERROR;
    }

    bu_vls_trunc(gedp->ged_result_str, 0);
    argc--;
    argv++;
    if (argc == 0) {
	const SbPlane plane = controller->getCuttingPlane();
	const SbVec3f point = cutting_plane_point(plane);
	const SbVec3f normal = plane.getNormal();
	bu_vls_printf(gedp->ged_result_str,
	    "enable %d\norigin %.17g %.17g %.17g\nnormal %.17g %.17g %.17g",
	    controller->isCuttingPlaneEnabled() ? 1 : 0,
	    static_cast<double>(point[0]), static_cast<double>(point[1]),
	    static_cast<double>(point[2]), static_cast<double>(normal[0]),
	    static_cast<double>(normal[1]), static_cast<double>(normal[2]));
	return BRLCAD_OK;
    }
    if (argc == 2 && BU_STR_EQUAL(argv[0], "enable")) {
	if (BU_STR_EQUAL(argv[1], "0"))
	    controller->setCuttingPlaneEnabled(FALSE);
	else if (BU_STR_EQUAL(argv[1], "1"))
	    controller->setCuttingPlaneEnabled(TRUE);
	else {
	    bu_vls_printf(gedp->ged_result_str,
		"usage: %s\n", usage);
	    return BRLCAD_ERROR;
	}
    } else if (argc == 4 &&
	(BU_STR_EQUAL(argv[0], "origin") || BU_STR_EQUAL(argv[0], "normal"))) {
	SbVec3f value;
	if (!cutting_vector(gedp, argv + 1, &value))
	    return BRLCAD_ERROR;
	const SbPlane current = controller->getCuttingPlane();
	const SbVec3f point = BU_STR_EQUAL(argv[0], "origin") ? value :
	    cutting_plane_point(current);
	SbVec3f normal = BU_STR_EQUAL(argv[0], "normal") ? value :
	    current.getNormal();
	if (normal.length() <= 1.0e-6f) {
	    bu_vls_printf(gedp->ged_result_str,
		"cutting-plane normal must be non-zero\n");
	    return BRLCAD_ERROR;
	}
	normal.normalize();
	if (!controller->setCuttingPlane(SbPlane(normal, point)))
	    return BRLCAD_ERROR;
    } else {
	bu_vls_printf(gedp->ged_result_str, "usage: %s\n", usage);
	return BRLCAD_ERROR;
    }

    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    if (view)
	(void)bv_refresh_request(view, GED_VIEW_REFRESH_VIEW);
    ged_refresh_cb(gedp);
    return BRLCAD_OK;
}

