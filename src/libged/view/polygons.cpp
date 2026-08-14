/*                    P O L Y G O N S . C P P
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
/** @file libged/view/polygons.cpp
 *
 * Commands for view polygons.
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
#include "bg/plane.h"
#include "bg/polygon_types.h"
#include "bv.h"
#include "rt/geom.h"

#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"

#include "../ged_bobol_private.hpp"
#include "../ged_private.h"
#include "./ged_view.h"


using _poly_binding = ged_bobol_polygon_binding;

static _poly_binding
_poly_find(struct _ged_view_info *gd, const char *name, bool local_only)
{
    return gd ? ged_bobol_polygon_find(gd->gedp, gd->cv, name,
	local_only) : _poly_binding();
}

static _poly_binding
_poly_binding_get(struct _ged_view_info *gd)
{
    return _poly_find(gd, gd ? gd->vobj : nullptr,
	gd ? gd->local_obj != 0 : false);
}

static BObolViewController *
_poly_scope_controller(struct _ged_view_info *gd)
{
    if (!gd || !gd->gedp || !gd->cv)
	return nullptr;
    return gd->local_obj ? ged_bobol_view_controller(gd->cv) :
	ged_bobol_shared_view_controller(gd->gedp);
}

static bool
_poly_record(struct _ged_view_info *gd, BObolPolygonRecord *record)
{
    const _poly_binding binding = _poly_binding_get(gd);
    return binding.controller && record &&
	binding.controller->polygons().record(binding.handle, *record);
}

static bool
_poly_exists(struct _ged_view_info *gd)
{
    return _poly_binding_get(gd).handle.isValid();
}

static BObolPolygonType
_poly_type(int type)
{
    switch (type) {
	case GED_VIEW_POLYGON_CIRCLE:
	    return BObolPolygonType::Circle;
	case GED_VIEW_POLYGON_ELLIPSE:
	    return BObolPolygonType::Ellipse;
	case GED_VIEW_POLYGON_RECTANGLE:
	    return BObolPolygonType::Rectangle;
	case GED_VIEW_POLYGON_SQUARE:
	    return BObolPolygonType::Square;
	default:
	    return BObolPolygonType::General;
    }
}

static BObolPolygonUpdate
_poly_update_type(int op)
{
    switch (op) {
	case GED_VIEW_POLYGON_UPDATE_PROPS_ONLY:
	    return BObolPolygonUpdate::PropsOnly;
	case GED_VIEW_POLYGON_UPDATE_PT_SELECT:
	    return BObolPolygonUpdate::PointSelect;
	case GED_VIEW_POLYGON_UPDATE_PT_SELECT_CLEAR:
	    return BObolPolygonUpdate::PointSelectClear;
	case GED_VIEW_POLYGON_UPDATE_PT_MOVE:
	    return BObolPolygonUpdate::PointMove;
	case GED_VIEW_POLYGON_UPDATE_PT_APPEND:
	    return BObolPolygonUpdate::PointAppend;
	default:
	    return BObolPolygonUpdate::Default;
    }
}

static bool
_poly_update(struct _ged_view_info *gd, int op)
{
    const _poly_binding binding = _poly_binding_get(gd);
    return binding.controller && binding.controller->polygons().update(
	binding.handle, _poly_update_type(op));
}

static bool
_poly_update_screen(struct _ged_view_info *gd, int x, int y, int op)
{
    const _poly_binding binding = _poly_binding_get(gd);
    if (!binding.controller || !gd || !gd->cv)
	return false;

    point_t model_point = VINIT_ZERO;
    if (!bv_screen_to_model(model_point,
	bv_context_view_const(ged_view_context_bv_const(gd->cv)),
	static_cast<fastf_t>(x), static_cast<fastf_t>(y)))
	return false;

    return binding.controller->polygons().updateModelPoint(binding.handle,
	SbVec3f(static_cast<float>(model_point[X]),
	    static_cast<float>(model_point[Y]),
	    static_cast<float>(model_point[Z])), _poly_update_type(op));
}

static void
_poly_project_point(point_t projected, plane_t *view_plane,
    struct ged_view_context *view_ctx, const point_t input)
{
    VMOVE(projected, input);
    if (!view_plane)
	return;
    HSET(*view_plane, 0.0, 0.0, 1.0, input[Z]);
    const struct bv *view = view_ctx ?
	bv_context_view_const(ged_view_context_bv_const(view_ctx)) : nullptr;
    if (!view || bv_plane_get(view_plane, view) != 0)
	return;

    fastf_t fx = 0.0;
    fastf_t fy = 0.0;
    point_t input_copy;
    point_t plane_point;
    VMOVE(input_copy, input);
    bg_plane_closest_pt(&fx, &fy, view_plane, &input_copy);
    bg_plane_pt_at(&plane_point, view_plane, fx, fy);
    VMOVE(projected, plane_point);
}

static SbColor
_poly_color(const struct bu_color *color)
{
    fastf_t rgb[3] = {0.0, 0.0, 0.0};
    if (color)
	bu_color_to_rgb_floats(color, rgb);
    return SbColor(static_cast<float>(rgb[0]),
	static_cast<float>(rgb[1]), static_cast<float>(rgb[2]));
}

static bool
_poly_color(struct bu_color *color, const SbColor &rgb)
{
    if (!color)
	return false;
    fastf_t values[3] = {
	static_cast<fastf_t>(rgb[0]),
	static_cast<fastf_t>(rgb[1]),
	static_cast<fastf_t>(rgb[2])
    };
    return bu_color_from_rgb_floats(color, values) != 0;
}

int
_poly_cmd_create(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view polygon create <name> x y [circ|ell|rect|sq]";
    const char *purpose_string = "create polygon";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (_poly_exists(gd)) {
	bu_vls_printf(gedp->ged_result_str, "View feature named %s already exists\n", gd->vobj);
	return BRLCAD_ERROR;
    }

    if (argc != 2 && argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int x,y;
    if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&x) != 1 || x < 0) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }
    if (bu_opt_int(NULL, 1, (const char **)&argv[1], (void *)&y) != 1 || y < 0) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[1]);
	return BRLCAD_ERROR;
    }

    point_t sp;
    if (!bv_screen_to_model(sp,
	    bv_context_view_const(ged_view_context_bv_const(gd->cv)),
	    (fastf_t)x, (fastf_t)y)) {
	bu_vls_printf(gedp->ged_result_str, "Failed to calculate screen point\n");
	return BRLCAD_ERROR;
    }

    int type = GED_VIEW_POLYGON_GENERAL;
    if (argc == 3) {
	if (BU_STR_EQUAL(argv[2], "circ") || BU_STR_EQUAL(argv[2], "circle"))
	    type = GED_VIEW_POLYGON_CIRCLE;
	if (BU_STR_EQUAL(argv[2], "ell") || BU_STR_EQUAL(argv[2], "ellipse"))
	    type = GED_VIEW_POLYGON_ELLIPSE;
	if (BU_STR_EQUAL(argv[2], "rect") || BU_STR_EQUAL(argv[2], "rectangle"))
	    type = GED_VIEW_POLYGON_RECTANGLE;
	if (BU_STR_EQUAL(argv[2], "sq") || BU_STR_EQUAL(argv[2], "square"))
	    type = GED_VIEW_POLYGON_SQUARE;
	if (type == GED_VIEW_POLYGON_GENERAL) {
	    bu_vls_printf(gedp->ged_result_str, "Unknown polygon type %s\n", argv[2]);
	    return BRLCAD_ERROR;
	}
    }

    BObolViewController *controller = _poly_scope_controller(gd);
    point_t origin = VINIT_ZERO;
    plane_t view_plane = HINIT_ZERO;
    _poly_project_point(origin, &view_plane, gd->cv, sp);
    const BObolPolygonHandle handle = controller ?
	controller->polygons().create(gd->vobj,
	    gd->local_obj ? BObolFeatureScope::Local :
	    BObolFeatureScope::Shared,
	    _poly_type(type),
	    SbVec3f(static_cast<float>(origin[X]),
		static_cast<float>(origin[Y]),
		static_cast<float>(origin[Z])), view_plane, 0.0f) :
	BObolPolygonHandle();
    if (!handle.isValid()) {
	bu_vls_printf(gedp->ged_result_str, "Failed to create %s\n", gd->vobj);
	return BRLCAD_ERROR;
    }

    return BRLCAD_OK;
}

int
_poly_cmd_select(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view polygon select <name> [contour] x y";
    const char *purpose_string = "select polygon point closest to point x,y";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (!_poly_exists(gd)) {
	bu_vls_printf(gedp->ged_result_str, "View feature %s does not exist\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    BObolPolygonRecord rec;
    if (!_poly_record(gd, &rec)) {
	bu_vls_printf(gedp->ged_result_str, "Specified object is not a view polygon.\n");
	return BRLCAD_ERROR;
    }

    if (rec.type != BObolPolygonType::General) {
	bu_vls_printf(gedp->ged_result_str, "Point selection is only supported for general polygons - specified object defines a constrained shape\n");
	return BRLCAD_ERROR;
    }

    if (argc != 2 && argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int ioffset = 0;
    int contour_ind = -1;
    if (argc == 3) {
	ioffset++;
	if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&contour_ind) != 1 || contour_ind < 0) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	    return BRLCAD_ERROR;
	}
    } else {
	contour_ind = 0;
    }
    int x,y;
    if (bu_opt_int(NULL, 1, (const char **)&argv[ioffset], (void *)&x) != 1 || x < 0) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[ioffset]);
	return BRLCAD_ERROR;
    }
    if (bu_opt_int(NULL, 1, (const char **)&argv[ioffset+1], (void *)&y) != 1 || y < 0) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[ioffset+1]);
	return BRLCAD_ERROR;
    }

    const _poly_binding binding = _poly_binding_get(gd);
    if (!binding.controller || !binding.controller->polygons().setCurrent(
	    binding.handle, contour_ind, -1))
	return BRLCAD_ERROR;
    if (!_poly_update_screen(gd, x, y, GED_VIEW_POLYGON_UPDATE_PT_SELECT))
	return BRLCAD_ERROR;

    return BRLCAD_OK;
}


int
_poly_cmd_append(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view polygon append <name> [contour] x y";
    const char *purpose_string = "append point to polygon";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (!_poly_exists(gd)) {
	bu_vls_printf(gedp->ged_result_str, "View feature %s does not exist\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    BObolPolygonRecord rec;
    if (!_poly_record(gd, &rec)) {
	bu_vls_printf(gedp->ged_result_str, "Specified object is not a view polygon.\n");
	return BRLCAD_ERROR;
    }

    if (rec.type != BObolPolygonType::General) {
	bu_vls_printf(gedp->ged_result_str, "Point appending is only supported for general polygons - specified object defines a constrained shape\n");
	return BRLCAD_ERROR;
    }

    if (argc != 2 && argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int ioffset = 0;
    int contour_ind = -1;
    if (argc == 3) {
	ioffset++;
	if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&contour_ind) != 1 || contour_ind < 0) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	    return BRLCAD_ERROR;
	}
    } else {
	contour_ind = 0;
    }
    int x,y;
    if (bu_opt_int(NULL, 1, (const char **)&argv[ioffset], (void *)&x) != 1 || x < 0) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[ioffset]);
	return BRLCAD_ERROR;
    }
    if (bu_opt_int(NULL, 1, (const char **)&argv[ioffset+1], (void *)&y) != 1 || y < 0) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[ioffset+1]);
	return BRLCAD_ERROR;
    }

    const _poly_binding binding = _poly_binding_get(gd);
    if (!binding.controller || !binding.controller->polygons().setCurrent(
	    binding.handle, contour_ind, -1))
	return BRLCAD_ERROR;
    if (!_poly_update_screen(gd, x, y, GED_VIEW_POLYGON_UPDATE_PT_APPEND))
	return BRLCAD_ERROR;

    return BRLCAD_OK;
}

int
_poly_cmd_move(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view polygon move <name> x y";
    const char *purpose_string = "move selected polygon point to x,y";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (!_poly_exists(gd)) {
	bu_vls_printf(gedp->ged_result_str, "View feature %s does not exist\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    BObolPolygonRecord rec;
    if (!_poly_record(gd, &rec)) {
	bu_vls_printf(gedp->ged_result_str, "Specified object is not a view polygon.\n");
	return BRLCAD_ERROR;
    }

    if (rec.type != BObolPolygonType::General) {
	bu_vls_printf(gedp->ged_result_str, "Individual point movement is only supported for general polygons - specified object defines a constrained shape.  Use \"update\" to adjust constrained shapes.\n");
	return BRLCAD_ERROR;
    }

    if (argc != 2) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int x,y;
    if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&x) != 1 || x < 0) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }
    if (bu_opt_int(NULL, 1, (const char **)&argv[1], (void *)&y) != 1 || y < 0) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[1]);
	return BRLCAD_ERROR;
    }

    if (!_poly_update_screen(gd, x, y, GED_VIEW_POLYGON_UPDATE_PT_MOVE))
	return BRLCAD_ERROR;

    return BRLCAD_OK;
}

int
_poly_cmd_clear(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view polygon clear <name>";
    const char *purpose_string = "clear all modification flags";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (!_poly_exists(gd)) {
	bu_vls_printf(gedp->ged_result_str, "View feature %s does not exist\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    const _poly_binding binding = _poly_binding_get(gd);
    if (!binding.controller || !binding.controller->polygons().setCurrent(
	    binding.handle, 0, -1)) {
	bu_vls_printf(gedp->ged_result_str, "Specified object is not a view polygon.\n");
	return BRLCAD_ERROR;
    }

    _poly_update(gd, GED_VIEW_POLYGON_UPDATE_DEFAULT);

    return BRLCAD_OK;
}

int
_poly_cmd_close(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view polygon close <name> [ind]";
    const char *purpose_string = "contour -> polygon";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (!_poly_exists(gd)) {
	bu_vls_printf(gedp->ged_result_str, "View feature %s does not exist\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    BObolPolygonRecord rec;
    if (!_poly_record(gd, &rec)) {
	bu_vls_printf(gedp->ged_result_str, "Specified object is not a view polygon.\n");
	return BRLCAD_ERROR;
    }

    if (rec.type != BObolPolygonType::General) {
	return BRLCAD_OK;
    }

    int ind = -1;
    if (argc == 1) {
	if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&ind) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	    return BRLCAD_ERROR;
	}
	if (ind < 0 || ind >= (int)rec.contourCount) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	    return BRLCAD_ERROR;
	}
    } else if (argc > 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    const _poly_binding binding = _poly_binding_get(gd);
    if (!binding.controller)
	return BRLCAD_ERROR;
    if (ind < 0) {
	if (!binding.controller->polygons().setAllContoursOpen(
		binding.handle, FALSE))
	    return BRLCAD_ERROR;
    } else if (!binding.controller->polygons().setContourOpen(
	    binding.handle, ind, FALSE)) {
	return BRLCAD_ERROR;
    }

    return BRLCAD_OK;
}

int
_poly_cmd_open(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view polygon open <name> [ind]";
    const char *purpose_string = "polygon -> contour";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (!_poly_exists(gd)) {
	bu_vls_printf(gedp->ged_result_str, "View feature %s does not exist\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    BObolPolygonRecord rec;
    if (!_poly_record(gd, &rec)) {
	bu_vls_printf(gedp->ged_result_str, "Specified object is not a view polygon.\n");
	return BRLCAD_ERROR;
    }

    if (rec.type != BObolPolygonType::General) {
	bu_vls_printf(gedp->ged_result_str, "Constrained polygon shapes are always closed.\n");
	return BRLCAD_ERROR;
    }

    int ind = -1;
    if (argc == 1) {
	if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&ind) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	    return BRLCAD_ERROR;
	}
	if (ind < 0 || ind >= (int)rec.contourCount) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	    return BRLCAD_ERROR;
	}
    } else if (argc > 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    const _poly_binding binding = _poly_binding_get(gd);
    if (!binding.controller)
	return BRLCAD_ERROR;
    if (ind < 0) {
	if (!binding.controller->polygons().setAllContoursOpen(
		binding.handle, TRUE))
	    return BRLCAD_ERROR;
    } else if (!binding.controller->polygons().setContourOpen(
	    binding.handle, ind, TRUE)) {
	return BRLCAD_ERROR;
    }

    return BRLCAD_OK;
}

int
_poly_cmd_area(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view polygon area <name>";
    const char *purpose_string = "report polygon area";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (!_poly_exists(gd)) {
	bu_vls_printf(gedp->ged_result_str, "View feature %s does not exist\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    const _poly_binding binding = _poly_binding_get(gd);
    const struct bv *view = gd->cv ?
	bv_context_view_const(ged_view_context_bv_const(gd->cv)) : nullptr;
    if (!binding.controller || !view) {
	bu_vls_printf(gedp->ged_result_str, "Specified object is not a view polygon.\n");
	return BRLCAD_ERROR;
    }
    const fastf_t area = binding.controller->polygons().area(
	binding.handle, bv_scale_get(view));

    if (gedp->dbip) {
	bu_vls_printf(gedp->ged_result_str, "%g", area * gedp->dbip->dbi_base2local);
    } else {
	bu_vls_printf(gedp->ged_result_str, "%g", area);
    }
    return BRLCAD_OK;
}

int
_poly_cmd_overlap(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view polygon overlap <obj1> <obj2>";
    const char *purpose_string = "report if two polygons overlap";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (!_poly_exists(gd)) {
	bu_vls_printf(gedp->ged_result_str, "View feature %s does not exist\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    if (!argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    const _poly_binding other = _poly_find(gd, argv[0], false);
    if (!other.handle.isValid()) {
	bu_vls_printf(gedp->ged_result_str, "View feature %s does not exist\n", argv[0]);
	return BRLCAD_ERROR;
    }

    // Have two polygons.  Check for overlaps, using the origin plane of the
    // obj1 polygon.
    const _poly_binding binding = _poly_binding_get(gd);
    const BObolPolygonHandle other_handle = binding.controller ?
	binding.controller->polygons().find(argv[0]) : BObolPolygonHandle();
    const struct bv *view = gd->cv ?
	bv_context_view_const(ged_view_context_bv_const(gd->cv)) : nullptr;
    if (!binding.controller || !other_handle.isValid() || !view) {
	bu_vls_printf(gedp->ged_result_str, "%s is not a view polygon.\n", argv[0]);
	return BRLCAD_ERROR;
    }
    const struct bn_tol tol = BN_TOL_INIT_TOL;
    const int ovlp = binding.controller->polygons().overlaps(
	binding.handle, other_handle, tol,
	bv_scale_get(view)) ? 1 : 0;

    bu_vls_printf(gedp->ged_result_str, "%d", ovlp);

    return BRLCAD_OK;
}

int
_poly_cmd_import(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view polygon import <name> <sketchname>";
    const char *purpose_string = "import polygon from sketch";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (_poly_exists(gd)) {
	bu_vls_printf(gedp->ged_result_str, "View feature named %s already exists\n", gd->vobj);
	return BRLCAD_ERROR;
    }

    if (!argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    if (!gedp->dbip) {
	bu_vls_printf(gedp->ged_result_str, "no database open\n");
	return BRLCAD_ERROR;
    }

    // Begin import
    struct directory *dp = db_lookup(gedp->dbip, argv[0], LOOKUP_QUIET);
    if (dp == RT_DIR_NULL) {
	return BRLCAD_ERROR;
    }

    BObolViewController *controller = _poly_scope_controller(gd);
    const BObolPolygonHandle handle = controller ?
	controller->polygons().importSketch(gd->vobj,
	    gd->local_obj ? BObolFeatureScope::Local :
	    BObolFeatureScope::Shared, gedp->dbip, dp) :
	BObolPolygonHandle();
    if (!handle.isValid()) {
	bu_vls_printf(gedp->ged_result_str, "Failed to create %s\n", gd->vobj);
	return BRLCAD_ERROR;
    }

    return BRLCAD_OK;
}

int
_poly_cmd_export(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view polygon export <name> <sketchname>";
    const char *purpose_string = "export polygon to sketch";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (!_poly_exists(gd)) {
	bu_vls_printf(gedp->ged_result_str, "View feature %s does not exist\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    BObolPolygonRecord record;
    if (!_poly_record(gd, &record)) {
	bu_vls_printf(gedp->ged_result_str, "Specified object is not a view polygon.\n");
	return BRLCAD_ERROR;
    }

    if (!argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    if (!gedp->dbip) {
	bu_vls_printf(gedp->ged_result_str, "no database open\n");
	return BRLCAD_ERROR;
    }

    GED_CHECK_EXISTS(gedp, argv[0], LOOKUP_QUIET, BRLCAD_ERROR);

    const _poly_binding binding = _poly_binding_get(gd);
    if (!binding.controller || !binding.controller->polygons().exportSketch(
	    binding.handle, gedp->dbip, argv[0])) {
	bu_vls_printf(gedp->ged_result_str, "Failed to create sketch.\n");
	return BRLCAD_ERROR;
    }

    return BRLCAD_OK;
}

int
_poly_cmd_fill(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view polygon fill <name> [dx dy spacing]";
    const char *purpose_string = "use lines to visualize polygon interior";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (!_poly_exists(gd)) {
	bu_vls_printf(gedp->ged_result_str, "View feature %s does not exist\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    BObolPolygonRecord rec;
    if (!_poly_record(gd, &rec)) {
	bu_vls_printf(gedp->ged_result_str, "Specified object is not a view polygon.\n");
	return BRLCAD_ERROR;
    }

    const _poly_binding binding = _poly_binding_get(gd);
    if (!binding.controller)
	return BRLCAD_ERROR;

    if (argc == 1 && BU_STR_EQUAL(argv[0], "0")) {
	if (!binding.controller->polygons().setFill(binding.handle, FALSE,
		rec.fillSlope, rec.fillSpacing))
	    return BRLCAD_ERROR;
	return BRLCAD_OK;
    }

    if (argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }
    vect2d_t vdir;
    fastf_t vdelta;
    if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[0], (void *)&vdir[0]) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }
    if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[1], (void *)&vdir[1]) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[1]);
	return BRLCAD_ERROR;
    }
    if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[2], (void *)&vdelta) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[2]);
	return BRLCAD_ERROR;
    }

    if (!binding.controller->polygons().setFill(binding.handle, TRUE,
	    SbVec2f(static_cast<float>(vdir[0]),
		static_cast<float>(vdir[1])), static_cast<float>(vdelta)))
	return BRLCAD_ERROR;

    return BRLCAD_OK;
}

int
_poly_cmd_fill_color(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view polygon fill_color <name> [r/g/b]";
    const char *purpose_string = "customize fill lines color (if fill is enabled)";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (!_poly_exists(gd)) {
	bu_vls_printf(gedp->ged_result_str, "View feature %s does not exist\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    BObolPolygonRecord record;
    if (!_poly_record(gd, &record)) {
	bu_vls_printf(gedp->ged_result_str, "Specified object is not a view polygon.\n");
	return BRLCAD_ERROR;
    }

    if (!argc) {
	unsigned char frgb[3];
	struct bu_color fill_color;
	if (!_poly_color(&fill_color, record.fillColor))
	    return BRLCAD_ERROR;
	bu_color_to_rgb_chars(&fill_color, (unsigned char *)frgb);

	bu_vls_printf(gedp->ged_result_str, "%d/%d/%d\n", frgb[0], frgb[1], frgb[2]);

	return BRLCAD_OK;
    }

    struct bu_color fill_color;
    if (bu_opt_color(NULL, 1, (const char **)&argv[0], (void *)&fill_color) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    const _poly_binding binding = _poly_binding_get(gd);
    if (!binding.controller || !binding.controller->polygons().setFillColor(
	    binding.handle, _poly_color(&fill_color)))
	return BRLCAD_ERROR;

    return BRLCAD_OK;
}

int
_poly_cmd_csg(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view polygon csg <obj1> <u|-|+> <obj2>";
    const char *purpose_string = "replace obj1's polygon with the result of obj2 u/-/+ obj1";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (!_poly_exists(gd)) {
	bu_vls_printf(gedp->ged_result_str, "View feature %s does not exist\n", gd->vobj);
	return BRLCAD_ERROR;
    }
    BObolPolygonRecord record;
    if (!_poly_record(gd, &record)) {
	bu_vls_printf(gedp->ged_result_str, "Specified object is not a view polygon.\n");
	return BRLCAD_ERROR;
    }

    if (argc != 2) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    enum bg_polygon_boolean_op op = BG_POLYGON_BOOLEAN_UNION;
    char c = argv[0][0];
    switch (c) {
	case 'u':
	    op = BG_POLYGON_BOOLEAN_UNION;
	    break;
	case '-':
	    op = BG_POLYGON_BOOLEAN_DIFFERENCE;
	    break;
	case '+':
	    op = BG_POLYGON_BOOLEAN_INTERSECTION;
	    break;
	default:
	    bu_vls_printf(gedp->ged_result_str, "Invalid boolean operator \"%s\"\n", argv[0]);
	    break;
    }

    const _poly_binding other = _poly_find(gd, argv[1], false);
    if (!other.handle.isValid()) {
	bu_vls_printf(gedp->ged_result_str, "View feature %s does not exist\n", argv[1]);
	return BRLCAD_ERROR;
    }

    const _poly_binding target = _poly_binding_get(gd);
    const BObolPolygonHandle stencil = target.controller ?
	target.controller->polygons().find(argv[1]) : BObolPolygonHandle();
    if (!target.controller || !stencil.isValid() ||
	    !target.controller->polygons().csg(target.handle, stencil, op))
	return BRLCAD_ERROR;

    return BRLCAD_OK;
}


const struct bu_cmdtab _poly_cmds[] = {
    { "append",          _poly_cmd_append},
    { "area",            _poly_cmd_area},
    { "clear",           _poly_cmd_clear},
    { "close",           _poly_cmd_close},
    { "create",          _poly_cmd_create},
    { "csg",             _poly_cmd_csg},
    { "export",          _poly_cmd_export},
    { "fill",            _poly_cmd_fill},
    { "fill_color",      _poly_cmd_fill_color},
    { "import",          _poly_cmd_import},
    { "move",            _poly_cmd_move},
    { "open",            _poly_cmd_open},
    { "overlap",         _poly_cmd_overlap},
    { "select",          _poly_cmd_select},
    { (char *)NULL,      NULL}
};

int
_view_cmd_polygons(void *bs, int argc, const char **argv)
{
    int help = 0;
    int s_version = 0;
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;

    const char *usage_string = "view polygon <verb> <name> [options] [args]";
    const char *purpose_string = "manipulate view polygons";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    if (!gd->cv) {
	bu_vls_printf(gedp->ged_result_str, ": no view current in GED");
	return BRLCAD_ERROR;
    }


    // We know we're the polygons command - start processing args
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
	if (bu_cmd_valid(_poly_cmds, argv[i]) == BRLCAD_OK) {
	    cmd_pos = i;
	    break;
	}
    }

    int acnt = (cmd_pos >= 0) ? cmd_pos : argc;
    (void)bu_opt_parse(NULL, acnt, argv, d);

    if (help) {
	_ged_cmd_help(gedp, usage_string, d);
	return GED_HELP;
    }

    return _ged_subcmd_exec(gedp, d, _poly_cmds, "view polygon <name>", "[options] subcommand [args]", gd, argc, argv, help, cmd_pos);
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
