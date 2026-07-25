/*                G E O M E T R Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2023-2026 United States Government as represented by
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
/** @file libged/brep/geometry.cpp
 *
 * NURBS geometry editing support for brep objects
 *
 */

#include "common.h"

#include "wdb.h"
#include "../ged_private.h"
#include "./ged_brep.h"

using namespace brlcad;

struct _ged_brep_igeo {
    struct _ged_brep_info *gb;
    struct bu_vls *vls;
    const struct bu_cmdtab *cmds;
};

static int
_brep_geo_msgs(void *bs, int argc, const char **argv, const char *us, const char *ps)
{
    struct _ged_brep_igeo *gb = (struct _ged_brep_igeo *)bs;
    if (argc == 2 && BU_STR_EQUAL(argv[1], HELPFLAG)) {
	bu_vls_printf(gb->vls, "%s\n%s\n", us, ps);
	return 1;
    }
    if (argc == 2 && BU_STR_EQUAL(argv[1], PURPOSEFLAG)) {
	bu_vls_printf(gb->vls, "%s\n", ps);
	return 1;
    }
    return 0;
}

static int
_brep_cmd_geo_vertex_create(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo v_create <x> <y> <z>";
    const char *purpose_string = "create a new vertex";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    argc--;argv++;
    if (argc < 3) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }
    ON_3dPoint position(atof(argv[0]), atof(argv[1]), atof(argv[2]));
    int vertex = brep_vertex_create(b_ip->brep, position);
    // Make the new one
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "create vertex! id = %d", vertex);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_vertex_remove(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo v_remove <v_id>";
    const char *purpose_string = "remove a vertex";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    argc--;argv++;
    if (argc < 1) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int v_id = atoi(argv[0]);
    bool res = brep_vertex_remove(b_ip->brep, v_id);
    if (!res) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to remove vertex %s\n", argv[0]);
	return BRLCAD_ERROR;
    }
    // Make the new one
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "remove vertex %d", v_id);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_curve2d_create_line(void *bs, int argc, const char **argv)
{

    const char *usage_string = "brep [options] <objname> geo c2_create_line <from_x> <from_y> <to_x> <to_y>";
    const char *purpose_string = "create a 2D parameter space geometric line";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;

    argc--;argv++;

    if (argc != 4) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    ON_2dPoint from(atof(argv[0]), atof(argv[1]));
    ON_2dPoint to(atof(argv[2]), atof(argv[3]));

    // Create a 2d line
    int curve_id = brep_curve2d_make_line(b_ip->brep, from, to);

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "create C2 curve! id = %d", curve_id);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_curve2d_remove(void *bs, int argc, const char **argv)
{

    const char *usage_string = "brep [options] <objname> geo remove <curve_id>";
    const char *purpose_string = "remove a 2D parameter space geometric curve";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;

    argc--;argv++;

    if (argc != 1) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }
    
    int curve_id = atoi(argv[0]);

    // Create a 2d line
    bool res = brep_curve2d_remove(b_ip->brep, curve_id);
    if (!res) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to remove curve %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "remove C2 curve %d", curve_id);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_curve3d_create(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo c3_create <x> <y> <z>";
    const char *purpose_string = "create a new NURBS curve";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;

    argc--;argv++;

    ON_3dPoint position(0, 0, 0);
    if (argc != 0 && argc != 3) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }
    if (argc == 3) {
	position = ON_3dPoint(atof(argv[0]), atof(argv[1]), atof(argv[2]));
    }
    // Create a template nurbs curve
    int curve_id = brep_curve_make(b_ip->brep, position);

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "create C3 curve! id = %d", curve_id);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_curve3d_in(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo c3_in <is_rational> <order> <cv_count> <cv1_x> <cv1_y> <cv1_z> <cv_w>(if rational) ...";
    const char *purpose_string = "create a new NURBS curve given detailed description";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 4) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    bool is_rational = atoi(argv[1]);
    int order = atoi(argv[2]);
    int cv_count = atoi(argv[3]);

    if (order <= 0 || cv_count <= 0 || cv_count < order) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid order or cv_count\n");
	return BRLCAD_ERROR;
    }

    if (argc < 4 + cv_count * (3 + (is_rational ? 1 : 0))) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments, you need to input %d more args about control vertices\n", 4 + cv_count * (3 + (is_rational ? 1 : 0)) - argc);
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;

    argc--;argv++;

    std::vector<ON_4dPoint> cv;
    for (int i = 0; i < cv_count; i++) {
	ON_4dPoint p;
	p.x = atof(argv[3 + i * (3 + (is_rational ? 1 : 0))]);
	p.y = atof(argv[3 + i * (3 + (is_rational ? 1 : 0)) + 1]);
	p.z = atof(argv[3 + i * (3 + (is_rational ? 1 : 0)) + 2]);
	if (is_rational) {
	    p.w = atof(argv[3 + i * (3 + (is_rational ? 1 : 0)) + 3]);
	} else {
	    p.w = 1.0;
	}
	cv.push_back(p);
    }
    int curve_id = brep_curve_in(b_ip->brep, is_rational, order, cv_count, cv);
    if (curve_id < 0) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to create curve\n");
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "create C3 curve! id = %d", curve_id);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_curve3d_interp(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo c3_interp <cv_count> <cv1_x> <cv1_y> <cv1_z> ...";
    const char *purpose_string = "create a new NURBS curve interpolating given control vertices";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 2) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    int cv_count = atoi(argv[1]);

    if (cv_count <= 0) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid cv_count\n");
	return BRLCAD_ERROR;
    }

    if (argc < 2 + cv_count * 3) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments, you need to input %d more args about control vertices\n", 2 + cv_count * 3 - argc);
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;

    argc--;argv++;

    std::vector<ON_3dPoint> points;
    for (int i = 0; i < cv_count; i++) {
	ON_3dPoint p;
	p.x = atof(argv[1 + i * 3]);
	p.y = atof(argv[1 + i * 3 + 1]);
	p.z = atof(argv[1 + i * 3 + 2]);
	points.push_back(p);
    }
    int curve_id = brep_curve_interpCrv(b_ip->brep, points);
    if (curve_id < 0) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to create curve\n");
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "create C3 curve! id = %d", curve_id);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_curve3d_copy(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo c3_copy <curve_id>";
    const char *purpose_string = "copy a NURBS curve";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 2) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    int curve_id = atoi(argv[1]);

    if (curve_id < 0) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid curve_id\n");
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;

    int res = brep_curve_copy(b_ip->brep, curve_id);
    if (!res) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to copy curve\n");
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "successful copy C3 curve! new curve id = %d", res);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_curve3d_remove(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo c3_remove <curve_id>";
    const char *purpose_string = "remove a NURBS curve";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 2) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    int curve_id = atoi(argv[1]);

    if (curve_id < 0) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid curve_id\n");
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;

    bool res = brep_curve_remove(b_ip->brep, curve_id);
    if (!res) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to remove curve\n");
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "successful remove C3 curve! id = %d", curve_id);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_curve3d_move(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo c3_move  <curve_id> <x> <y> <z>";
    const char *purpose_string = "move a NURBS curve to a specified position";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    argc--;argv++;
    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 4) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    ON_3dPoint p = ON_3dPoint(atof(argv[1]), atof(argv[2]), atof(argv[3]));
    bool flag = brep_curve_move(b_ip->brep, atoi(argv[0]), p);
    if (!flag) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to move curve %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_curve3d_set_cv(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo c3_set_cv <curve_id> <CV_id> <x> <y> <z> [<w>]";
    const char *purpose_string = "set the control vertex of a NURBS curve";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;

    if (argc < 6) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;

    ON_4dPoint p = ON_4dPoint(atof(argv[3]), atof(argv[4]), atof(argv[5]), argc == 7 ? atof(argv[6]) : 1.0);
    bool flag = brep_curve_set_cv(b_ip->brep, atoi(argv[1]), atoi(argv[2]), p);
    if (!flag) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to move control vertex %s of curve %s\n", argv[2], argv[1]);
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_curve3d_flip(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo c3_flip <curve_id>";
    const char *purpose_string = "Flip the direction of a NURBS curve";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 2) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    int curve_id = 0;
    try {
	curve_id = std::stoi(argv[1]);
    } catch (const std::exception &e) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid curve id\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", e.what());
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;

    bool flag = brep_curve_reverse(b_ip->brep, curve_id);
    if (!flag) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to reverse curve %s\n", argv[1]);
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_curve3d_insert_knot(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo c3_insert_knot <curve_id> <knot_value> <multiplicity>";
    const char *purpose_string = "Insert a knot into a NURBS curve";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }
    argc--;argv++;
    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 3) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    int curve_id = 0;
    try {
	curve_id = std::stoi(argv[0]);
    } catch (const std::exception &e) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid curve id\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", e.what());
	return BRLCAD_ERROR;
    }

    double knot_value = 0;
    try {
	knot_value = std::stod(argv[1]);
    } catch (const std::exception &e) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid knot value\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", e.what());
	return BRLCAD_ERROR;
    }

    int multiplicity = 0;
    try {
	multiplicity = std::stoi(argv[2]);
    } catch (const std::exception &e) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid multiplicity\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", e.what());
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    bool flag = brep_curve_insert_knot(b_ip->brep, curve_id, knot_value, multiplicity);
    if (!flag) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to insert knot into curve %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_curve3d_trim(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo c3_trim <curve_id> <start_param> <end_param>";
    const char *purpose_string = "trim a NURBS curve using start and end parameters";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }
    argc--;argv++;
    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 3) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    int curve_id = 0;
    try {
	curve_id = std::stoi(argv[0]);
    } catch (const std::exception &e) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid curve id\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", e.what());
	return BRLCAD_ERROR;
    }

    double start_param = 0;
    try {
	start_param = std::stod(argv[1]);
    } catch (const std::exception &e) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid start parameter\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", e.what());
	return BRLCAD_ERROR;
    }

    double end_param = 0;
    try {
	end_param = std::stod(argv[2]);
    } catch (const std::exception &e) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid end parameter\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", e.what());
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    bool flag = brep_curve_trim(b_ip->brep, curve_id, start_param, end_param);
    if (!flag) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to trim curve %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_curve3d_split(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo c3_split <curve_id> <param>";
    const char *purpose_string = "split a NURBS curve into two at a parameter";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }
    argc--;argv++;
    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 2) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    int curve_id = 0;
    try {
	curve_id = std::stoi(argv[0]);
    } catch (const std::exception &e) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid curve id\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", e.what());
	return BRLCAD_ERROR;
    }

    double param = 0;
    try {
	param = std::stod(argv[1]);
    } catch (const std::exception &e) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid parameter\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", e.what());
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    bool flag = brep_curve_split(b_ip->brep, curve_id, param);
    if (!flag) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to split curve %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "split curve %s at parameter %s. Old curve removed.\n", argv[0], argv[1]);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_curve3d_join(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo c3_join <curve_id_1> <curve_id_2>";
    const char *purpose_string = "join end of curve 1 to start of curve 2";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }
    argc--;argv++;
    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 2) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    int curve_id_1 = 0;
    int curve_id_2 = 0;
    try {
	curve_id_1 = std::stoi(argv[0]);
	curve_id_2 = std::stoi(argv[1]);
    } catch (const std::exception &e) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid curve id\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", e.what());
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    int flag = brep_curve_join(b_ip->brep, curve_id_1, curve_id_2);
    if (flag < 0) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to join curve %s and curve %s\n", argv[0], argv[1]);
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "joined curve id %d, old curves deleted.\n", flag);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_surface_create(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo s_create <x> <y> <z>";
    const char *purpose_string = "create a new NURBS surface";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    argc--;argv++;
    ON_3dPoint position(0, 0, 0);
    if (argc >= 3) {
	position = ON_3dPoint(atof(argv[0]), atof(argv[1]), atof(argv[2]));
    }
    int surf_id = brep_surface_make(b_ip->brep, position);

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "create surface! id = %d", surf_id);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_surface_interp(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo s_interp <cv_count_x> <cv_1_1_x> <cv_1_1_y> <cv_1_1_z> <cv_1_2_x> <cv_1_2_y> <cv_1_2_z> ...";
    const char *purpose_string = "create a new NURBS surface";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 3) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    argc--;argv++;
    // Load and check the number of arguments
    int cv_count_x = atoi(argv[0]);
    int cv_count_y = atoi(argv[1]);
    if (cv_count_x <= 2 || cv_count_y <= 2) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid cv_count, cv_count >= 3\n");
	return BRLCAD_ERROR;
    }
    
    if (argc < 2 + (cv_count_x *cv_count_y) * 3) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments, you need to input %d more args about control vertices\n", 2 + (cv_count_x * cv_count_y) * 3 - argc);
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;

    std::vector<ON_3dPoint> points;
    for (int i = 0; i < (cv_count_x * cv_count_y); i++) {
	ON_3dPoint p;
	p.x = atof(argv[2 + i * 3]);
	p.y = atof(argv[2 + i * 3 + 1]);
	p.z = atof(argv[2 + i * 3 + 2]);
	points.push_back(p);
    }
    int surface_id = brep_surface_interpCrv(b_ip->brep, cv_count_x, cv_count_y, points);
    if(surface_id < 0) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to create surface\n");
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "create surface! id = %d", surface_id);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_surface_copy(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo s_copy <surface_id>";
    const char *purpose_string = "copy a NURBS surface";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 2) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    int surface_id = atoi(argv[1]);

    if (surface_id < 0) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "invalid surface_id\n");
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    int res = brep_surface_copy(b_ip->brep, surface_id);
    if (res < 0) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, ": failed to copy surface\n");
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "successful copy surface! new surface id = %d", res);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_surface_birail(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo s_birail <curve_id_1> <curve_id_2>";
    const char *purpose_string = "create a new NURBS surface using two curves";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    argc--;argv++;
    if (argc != 2) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, " not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int curve_id_1 = atoi(argv[0]);
    int curve_id_2 = atoi(argv[1]);
    int surf_id = brep_surface_create_ruled(b_ip->brep, curve_id_1, curve_id_2);
    if (surf_id < 0) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to create surface\n");
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "create surface! id = %d", surf_id);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_surface_remove(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo s_remove <surface_id>";
    const char *purpose_string = "remove a NURBS surface";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    argc--;argv++;
    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 1) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, " not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    bool flag = brep_surface_remove(b_ip->brep, atoi(argv[0]));
    if (!flag) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, ": failed to remove surface %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_surface_move(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo s_move <surface_id> <x> <y> <z>";
    const char *purpose_string = "move a NURBS surface to a specified position";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    argc--;argv++;
    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 4) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, " not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    ON_3dPoint p = ON_3dPoint(atof(argv[1]), atof(argv[2]), atof(argv[3]));
    bool flag = brep_surface_move(b_ip->brep, atoi(argv[0]), p);
    if (!flag) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, ": failed to move surface %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_surface_set_cv(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo s_set_cv <surface_id> <cv_id_u> <cv_id_v> <x> <y> <z> [<w>]";
    const char *purpose_string = "set a control vertex of a NURBS surface to a specified position";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    argc--;argv++;
    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 6) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, " not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    ON_4dPoint p = ON_4dPoint(atof(argv[3]), atof(argv[4]), atof(argv[5]), argc == 7 ? atof(argv[6]) : 1);
    bool flag = brep_surface_set_cv(b_ip->brep, atoi(argv[0]), atoi(argv[1]), atoi(argv[2]), p);
    if (!flag) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, ": failed to move surface cv \n");
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_surface_trim(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo s_trim <surface_id> <dir> <start_param> <end_param>";
    const char *purpose_string = "trim a NURBS surface using start and end parameters.";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }
    argc--;argv++;
    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 4) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, " not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    int surface_id = atoi(argv[0]);
    int dir = atoi(argv[1]);
    double start_param = atof(argv[2]);
    double end_param = atof(argv[3]);
    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    bool flag = brep_surface_trim(b_ip->brep, surface_id, dir, start_param, end_param);
    if (!flag) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, ": failed to trim surface %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_surface_split(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo s_split <surface_id> <dir> <param>";
    const char *purpose_string = "split a NURBS surface into two given a parameter value.";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }
    argc--;argv++;
    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc < 3) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, " not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    int surface_id = atoi(argv[0]);
    int dir = atoi(argv[1]);
    double param = atof(argv[2]);
    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    bool flag = brep_surface_split(b_ip->brep, surface_id, dir, param);
    if (!flag) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, ": failed to split the surface %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_surface_tensor_product(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo s_tensor <curve_id_1> <curve_id_2>";
    const char *purpose_string = "create a new NURBS surface by extruding the first curve along the second curve.";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    argc--;argv++;
    if (argc != 2) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, " not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int curve_id_1 = atoi(argv[0]);
    int curve_id_2 = atoi(argv[1]);
    int surf_id = brep_surface_tensor_product(b_ip->brep, curve_id_1, curve_id_2);
    if (surf_id < 0) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, ": failed to create surface\n");
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "create surface! id = %d", surf_id);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_surface_revolution(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo s_revolution <curve_id> <start_x> <start_y> <start_z> <end_x> <end_y> <end_z> [<angle>]";
    const char *purpose_string = "create a new NURBS surface by rotating a curve around an axis by an angle.";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    argc--;argv++;
    if (argc < 7) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, " not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }
    int curve_id_1 = atoi(argv[0]);
    ON_3dPoint line_start(atof(argv[1]), atof(argv[2]), atof(argv[3]));
    ON_3dPoint line_end(atof(argv[4]), atof(argv[5]), atof(argv[6]));
    double angle = 2 * ON_PI;
    if(argc == 8) {
	angle = atof(argv[7]);
    }
    int surf_id = brep_surface_revolution(b_ip->brep, curve_id_1, line_start, line_end, angle);
    if (surf_id < 0) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to create surface\n");
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "create surface! id = %d", surf_id);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_surface_extract_vertex(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo s_ext_v <surface_id> <u> <v>";
    const char *purpose_string = "extract a vertex from a NURBS surface";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    argc--;argv++;
    
    if (argc < 3) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }
    
    int surface_id = atoi(argv[0]);
    double u = atof(argv[1]);
    double v = atof(argv[2]);
    int vertex_id = brep_surface_extract_vertex(b_ip->brep, surface_id, u, v);

    if(vertex_id < 0) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to create vertex\n");
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "create vertex! id = %d", vertex_id);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_surface_extract_curve(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo s_ext_c3 <surface_id> <dir> <param>";
    const char *purpose_string = "extract a curve from a NURBS surface";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    struct rt_brep_internal *b_ip = (struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    argc--;argv++;
    
    if (argc < 3) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "not enough arguments\n");
	bu_vls_printf(gib->gb->gedp->ged_result_str, "%s\n", usage_string);
	return BRLCAD_ERROR;
    }
    
    int surface_id = atoi(argv[0]);
    int dir = atoi(argv[1]);
    double param = atof(argv[2]);
    int curve_id = brep_surface_extract_curve(b_ip->brep, surface_id, dir, param);

    if(curve_id < 0) {
	bu_vls_printf(gib->gb->gedp->ged_result_str, "failed to create curve\n");
	return BRLCAD_ERROR;
    }

    // Update object in database
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK) {
	return BRLCAD_ERROR;
    }
    bu_vls_printf(gib->gb->gedp->ged_result_str, "create curve! id = %d", curve_id);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_nurbs_refine(void *bs, int argc, const char **argv)
{
    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    struct rt_brep_internal *b_ip =
	(struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    const char *cmd = argc > 0 ? argv[0] : "";
    const char *usage_string = NULL;
    const char *purpose_string = NULL;
    bool changed = false;

    if (BU_STR_EQUAL(cmd, "c3_set_weight")) {
	usage_string = "brep [options] <objname> geo c3_set_weight <curve_id> <cv_id> <weight>";
	purpose_string = "set a NURBS curve control vertex weight";
    } else if (BU_STR_EQUAL(cmd, "c3_clamp")) {
	usage_string = "brep [options] <objname> geo c3_clamp <curve_id> <end:0|1|2>";
	purpose_string = "clamp a NURBS curve start, end, or both without changing its locus";
    } else if (BU_STR_EQUAL(cmd, "c3_elevate")) {
	usage_string = "brep [options] <objname> geo c3_elevate <curve_id> <desired_degree>";
	purpose_string = "elevate a NURBS curve degree without changing its locus";
    } else if (BU_STR_EQUAL(cmd, "s_set_weight")) {
	usage_string = "brep [options] <objname> geo s_set_weight <surface_id> <cv_u> <cv_v> <weight>";
	purpose_string = "set an unreferenced NURBS surface control vertex weight";
    } else if (BU_STR_EQUAL(cmd, "s_reverse")) {
	usage_string = "brep [options] <objname> geo s_reverse <surface_id> <dir:0|1>";
	purpose_string = "reverse an unreferenced NURBS surface parameter direction";
    } else if (BU_STR_EQUAL(cmd, "s_transpose")) {
	usage_string = "brep [options] <objname> geo s_transpose <surface_id>";
	purpose_string = "transpose an unreferenced NURBS surface";
    } else if (BU_STR_EQUAL(cmd, "s_insert_knot")) {
	usage_string = "brep [options] <objname> geo s_insert_knot <surface_id> <dir:0|1> <knot> <multiplicity>";
	purpose_string = "insert a NURBS surface knot without changing its locus";
    } else if (BU_STR_EQUAL(cmd, "s_clamp")) {
	usage_string = "brep [options] <objname> geo s_clamp <surface_id> <dir:0|1> <end:0|1|2>";
	purpose_string = "clamp a NURBS surface direction without changing its locus";
    } else if (BU_STR_EQUAL(cmd, "s_elevate")) {
	usage_string = "brep [options] <objname> geo s_elevate <surface_id> <dir:0|1> <desired_degree>";
	purpose_string = "elevate a NURBS surface degree without changing its locus";
    } else {
	bu_vls_printf(gib->vls, "unknown NURBS refinement command\n");
	return BRLCAD_ERROR;
    }

    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    if (BU_STR_EQUAL(cmd, "c3_set_weight") && argc == 4) {
	changed = brep_curve_set_weight(b_ip->brep, atoi(argv[1]),
		atoi(argv[2]), atof(argv[3]));
    } else if (BU_STR_EQUAL(cmd, "c3_clamp") && argc == 3) {
	changed = brep_curve_clamp(b_ip->brep, atoi(argv[1]), atoi(argv[2]));
    } else if (BU_STR_EQUAL(cmd, "c3_elevate") && argc == 3) {
	changed = brep_curve_elevate_degree(b_ip->brep, atoi(argv[1]),
		atoi(argv[2]));
    } else if (BU_STR_EQUAL(cmd, "s_set_weight") && argc == 5) {
	changed = brep_surface_set_weight(b_ip->brep, atoi(argv[1]),
		atoi(argv[2]), atoi(argv[3]), atof(argv[4]));
    } else if (BU_STR_EQUAL(cmd, "s_reverse") && argc == 3) {
	changed = brep_surface_reverse(b_ip->brep, atoi(argv[1]), atoi(argv[2]));
    } else if (BU_STR_EQUAL(cmd, "s_transpose") && argc == 2) {
	changed = brep_surface_transpose(b_ip->brep, atoi(argv[1]));
    } else if (BU_STR_EQUAL(cmd, "s_insert_knot") && argc == 5) {
	changed = brep_surface_insert_knot(b_ip->brep, atoi(argv[1]),
		atoi(argv[2]), atof(argv[3]), atoi(argv[4]));
    } else if (BU_STR_EQUAL(cmd, "s_clamp") && argc == 4) {
	changed = brep_surface_clamp(b_ip->brep, atoi(argv[1]),
		atoi(argv[2]), atoi(argv[3]));
    } else if (BU_STR_EQUAL(cmd, "s_elevate") && argc == 4) {
	changed = brep_surface_elevate_degree(b_ip->brep, atoi(argv[1]),
		atoi(argv[2]), atoi(argv[3]));
    } else {
	bu_vls_printf(gib->vls, "invalid arguments\n%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    if (!changed) {
	bu_vls_printf(gib->vls,
		"%s failed or was rejected to preserve B-rep references\n", cmd);
	return BRLCAD_ERROR;
    }
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK)
	return BRLCAD_ERROR;
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_face_edit(void *bs, int argc, const char **argv)
{
    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    struct rt_brep_internal *b_ip =
	(struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    const char *cmd = argc > 0 ? argv[0] : "";
    const char *usage_string = NULL;
    const char *purpose_string = NULL;
    bool changed = false;

    if (BU_STR_EQUAL(cmd, "f_cv_safe")) {
	usage_string = "brep [options] <objname> geo f_cv_safe <face_id> <cv_u> <cv_v>";
	purpose_string = "report whether a face CV can be changed without moving a trim";
    } else if (BU_STR_EQUAL(cmd, "f_cv_editable")) {
	usage_string = "brep [options] <objname> geo f_cv_editable <face_id> <cv_u> <cv_v>";
	purpose_string = "report whether a face CV can be translated with current C0 coupling support";
    } else if (BU_STR_EQUAL(cmd, "f_cv_constraint")) {
	usage_string = "brep [options] <objname> geo f_cv_constraint <face_id> <cv_u> <cv_v>";
	purpose_string = "report trim influence and constrained-edit eligibility for a face CV";
    } else if (BU_STR_EQUAL(cmd, "f_set_cv")) {
	usage_string = "brep [options] <objname> geo f_set_cv <face_id> <cv_u> <cv_v> <x> <y> <z> [<w>]";
	purpose_string = "set a face CV using supported topology-preserving constraints";
    } else if (BU_STR_EQUAL(cmd, "f_move_cv")) {
	usage_string = "brep [options] <objname> geo f_move_cv <face_id> <cv_u> <cv_v> <dx> <dy> <dz>";
	purpose_string = "translate a face CV using supported topology-preserving constraints";
    } else if (BU_STR_EQUAL(cmd, "f_set_weight")) {
	usage_string = "brep [options] <objname> geo f_set_weight <face_id> <cv_u> <cv_v> <weight>";
	purpose_string = "set a topology-safe face CV weight, isolating a shared surface";
    } else if (BU_STR_EQUAL(cmd, "f_reverse_param")) {
	usage_string = "brep [options] <objname> geo f_reverse_param <face_id> <dir:0|1>";
	purpose_string = "reverse a face parameter direction and update its trims";
    } else if (BU_STR_EQUAL(cmd, "f_transpose")) {
	usage_string = "brep [options] <objname> geo f_transpose <face_id>";
	purpose_string = "transpose a face surface and its parameter-space trims";
    } else {
	bu_vls_printf(gib->vls, "unknown face editing command\n");
	return BRLCAD_ERROR;
    }

    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    if (BU_STR_EQUAL(cmd, "f_cv_safe") && argc == 4) {
	const bool safe = brep_face_cv_is_topology_safe(b_ip->brep,
		atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
	bu_vls_printf(gib->vls, "%d", safe ? 1 : 0);
	return BRLCAD_OK;
    } else if (BU_STR_EQUAL(cmd, "f_cv_editable") && argc == 4) {
	const bool editable = brep_face_cv_can_translate(b_ip->brep,
		atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
	bu_vls_printf(gib->vls, "%d", editable ? 1 : 0);
	return BRLCAD_OK;
    } else if (BU_STR_EQUAL(cmd, "f_cv_constraint") && argc == 4) {
	struct brep_face_cv_constraint status;
	if (!brep_face_cv_constraint_status(b_ip->brep,
		atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), &status)) {
	    bu_vls_printf(gib->vls, "invalid face or CV");
	    return BRLCAD_ERROR;
	}
	bu_vls_printf(gib->vls,
		"class %s trims %d edges %d adjacent_faces %d natural %d "
		"isoparametric %d general %d topology_safe %d can_translate %d",
		brep_cv_constraint_type_name(status.classification),
		status.trim_count, status.edge_count, status.other_face_count,
		status.natural_trim_count, status.isoparametric_trim_count,
		status.general_trim_count, status.topology_safe ? 1 : 0,
		status.can_translate ? 1 : 0);
	return BRLCAD_OK;
    } else if (BU_STR_EQUAL(cmd, "f_set_cv") && (argc == 7 || argc == 8)) {
	const int face_id = atoi(argv[1]);
	const int cv_u = atoi(argv[2]);
	const int cv_v = atoi(argv[3]);
	double weight = argc == 8 ? atof(argv[7]) : 1.0;
	if (argc == 7 && face_id >= 0
		&& face_id < b_ip->brep->m_F.Count()) {
	    const ON_NurbsSurface *surface =
		dynamic_cast<const ON_NurbsSurface *>(
		    b_ip->brep->m_F[face_id].SurfaceOf());
	    ON_4dPoint current;
	    if (surface && surface->GetCV(cv_u, cv_v,
			ON::euclidean_rational, &current.x))
		weight = current.w;
	}
	ON_4dPoint cv(atof(argv[4]), atof(argv[5]), atof(argv[6]), weight);
	changed = brep_face_set_cv_constrained(b_ip->brep, face_id,
		cv_u, cv_v, cv);
    } else if (BU_STR_EQUAL(cmd, "f_move_cv") && argc == 7) {
	ON_3dVector delta(atof(argv[4]), atof(argv[5]), atof(argv[6]));
	changed = brep_face_translate_cv_constrained(b_ip->brep, atoi(argv[1]),
		atoi(argv[2]), atoi(argv[3]), delta);
    } else if (BU_STR_EQUAL(cmd, "f_set_weight") && argc == 5) {
	changed = brep_face_set_weight(b_ip->brep, atoi(argv[1]),
		atoi(argv[2]), atoi(argv[3]), atof(argv[4]));
    } else if (BU_STR_EQUAL(cmd, "f_reverse_param") && argc == 3) {
	changed = brep_face_reverse_parameter(b_ip->brep, atoi(argv[1]),
		atoi(argv[2]));
    } else if (BU_STR_EQUAL(cmd, "f_transpose") && argc == 2) {
	changed = brep_face_transpose(b_ip->brep, atoi(argv[1]));
    } else {
	bu_vls_printf(gib->vls, "invalid arguments\n%s\n", usage_string);
	return BRLCAD_ERROR;
    }

    if (!changed) {
	bu_vls_printf(gib->vls,
		"%s failed or was rejected to preserve B-rep topology\n", cmd);
	return BRLCAD_ERROR;
    }
    if (brep_write_modified(gib->gb, b_ip->brep) != BRLCAD_OK)
	return BRLCAD_ERROR;
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_boundary_summary(void *bs, int argc, const char **argv)
{
    const char *usage_string =
	"brep [options] <objname> geo boundary_summary";
    const char *purpose_string =
	"classify B-rep edges for constrained C0 boundary editing";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc != 1) {
	bu_vls_printf(gib->vls, "invalid arguments\n%s\n", usage_string);
	return BRLCAD_ERROR;
    }
    struct rt_brep_internal *b_ip =
	(struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    int counts[BREP_EDGE_CONSTRAINT_SPECIAL + 1] = {0};
    int invalid = 0;
    for (int edge = 0; edge < b_ip->brep->m_E.Count(); ++edge) {
	const int classification =
	    brep_edge_constraint_type(b_ip->brep, edge);
	if (classification < 0 || classification > BREP_EDGE_CONSTRAINT_SPECIAL)
	    ++invalid;
	else
	    ++counts[classification];
    }
    bu_vls_printf(gib->vls,
	    "edges %d open %d natural %d isoparametric %d "
	    "one_isoparametric %d general %d special %d invalid %d",
	    b_ip->brep->m_E.Count(),
	    counts[BREP_EDGE_CONSTRAINT_OPEN],
	    counts[BREP_EDGE_CONSTRAINT_NATURAL],
	    counts[BREP_EDGE_CONSTRAINT_ISOPARAMETRIC],
	    counts[BREP_EDGE_CONSTRAINT_ONE_ISOPARAMETRIC],
	    counts[BREP_EDGE_CONSTRAINT_GENERAL],
	    counts[BREP_EDGE_CONSTRAINT_SPECIAL], invalid);
    return BRLCAD_OK;
}

static int
_brep_cmd_geo_validate(void *bs, int argc, const char **argv)
{
    const char *usage_string = "brep [options] <objname> geo validate";
    const char *purpose_string = "run the openNURBS B-rep validity checks";
    if (_brep_geo_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    struct _ged_brep_igeo *gib = (struct _ged_brep_igeo *)bs;
    if (argc != 1) {
	bu_vls_printf(gib->vls, "invalid arguments\n%s\n", usage_string);
	return BRLCAD_ERROR;
    }
    struct rt_brep_internal *b_ip =
	(struct rt_brep_internal *)gib->gb->intern.idb_ptr;
    const bool valid = brep_is_valid(b_ip->brep);
    bu_vls_printf(gib->vls, "%s", valid ? "valid" : "invalid");
    return valid ? BRLCAD_OK : BRLCAD_ERROR;
}

static void
_brep_geo_help(struct _ged_brep_igeo *bs, int argc, const char **argv)
{
    struct _ged_brep_igeo *gb = (struct _ged_brep_igeo *)bs;
    if (!argc || !argv) {
	bu_vls_printf(gb->vls, "brep [options] <objname> geo <subcommand> [args]\n");
	bu_vls_printf(gb->vls, "Available subcommands:\n");
	const struct bu_cmdtab *ctp = NULL;
	int ret;
	const char *helpflag[2];
	helpflag[1] = PURPOSEFLAG;
	for (ctp = gb->cmds; ctp->ct_name != (char *)NULL; ctp++) {
	    bu_vls_printf(gb->vls, "  %s\t\t\t", ctp->ct_name);
	    helpflag[0] = ctp->ct_name;
	    bu_cmd(gb->cmds, 2, helpflag, 0, (void *)gb, &ret);
	}
    }
    else {
	int ret;
	const char *helpflag[2];
	helpflag[0] = argv[0];
	helpflag[1] = HELPFLAG;
	bu_cmd(gb->cmds, 2, helpflag, 0, (void *)gb, &ret);
    }
}

const struct bu_cmdtab _brep_geo_cmds[] = {
    { "v_create",               _brep_cmd_geo_vertex_create},
    { "v_remove",               _brep_cmd_geo_vertex_remove},
    { "c2_create_line",         _brep_cmd_geo_curve2d_create_line},
    { "c2_remove",              _brep_cmd_geo_curve2d_remove},
    { "c3_create",              _brep_cmd_geo_curve3d_create},
    { "c3_in",                  _brep_cmd_geo_curve3d_in},
    { "c3_interp",              _brep_cmd_geo_curve3d_interp},
    { "c3_copy",                _brep_cmd_geo_curve3d_copy},
    { "c3_remove",              _brep_cmd_geo_curve3d_remove},
    { "c3_move",                _brep_cmd_geo_curve3d_move},
    { "c3_set_cv",              _brep_cmd_geo_curve3d_set_cv},
    { "c3_set_weight",          _brep_cmd_geo_nurbs_refine},
    { "c3_flip",                _brep_cmd_geo_curve3d_flip},
    { "c3_insert_knot",         _brep_cmd_geo_curve3d_insert_knot},
    { "c3_clamp",               _brep_cmd_geo_nurbs_refine},
    { "c3_elevate",             _brep_cmd_geo_nurbs_refine},
    { "c3_trim",                _brep_cmd_geo_curve3d_trim},
    { "c3_split",               _brep_cmd_geo_curve3d_split},
    { "c3_join",                _brep_cmd_geo_curve3d_join},
    { "s_create",               _brep_cmd_geo_surface_create},
    { "s_interp",               _brep_cmd_geo_surface_interp},
    { "s_copy",                 _brep_cmd_geo_surface_copy},
    { "s_birail",               _brep_cmd_geo_surface_birail},
    { "s_remove",               _brep_cmd_geo_surface_remove},
    { "s_move",                 _brep_cmd_geo_surface_move},
    { "s_set_cv",               _brep_cmd_geo_surface_set_cv},
    { "s_set_weight",           _brep_cmd_geo_nurbs_refine},
    { "s_reverse",              _brep_cmd_geo_nurbs_refine},
    { "s_transpose",            _brep_cmd_geo_nurbs_refine},
    { "s_insert_knot",          _brep_cmd_geo_nurbs_refine},
    { "s_clamp",                _brep_cmd_geo_nurbs_refine},
    { "s_elevate",              _brep_cmd_geo_nurbs_refine},
    { "s_trim",                 _brep_cmd_geo_surface_trim},
    { "s_split",                _brep_cmd_geo_surface_split},
    { "s_tensor",               _brep_cmd_geo_surface_tensor_product},
    { "s_revolution",           _brep_cmd_geo_surface_revolution},
    { "s_ext_v",                _brep_cmd_geo_surface_extract_vertex},
    { "s_ext_c3",               _brep_cmd_geo_surface_extract_curve},
    { "f_cv_safe",              _brep_cmd_geo_face_edit},
    { "f_cv_editable",          _brep_cmd_geo_face_edit},
    { "f_cv_constraint",        _brep_cmd_geo_face_edit},
    { "f_set_cv",               _brep_cmd_geo_face_edit},
    { "f_move_cv",              _brep_cmd_geo_face_edit},
    { "f_set_weight",           _brep_cmd_geo_face_edit},
    { "f_reverse_param",        _brep_cmd_geo_face_edit},
    { "f_transpose",            _brep_cmd_geo_face_edit},
    { "boundary_summary",       _brep_cmd_geo_boundary_summary},
    { "validate",               _brep_cmd_geo_validate},
    { (char *)NULL,             NULL}
};

int brep_geo(struct _ged_brep_info *gb, int argc, const char **argv)
{
    struct _ged_brep_igeo gib;
    gib.gb = gb;
    gib.vls = gb->gedp->ged_result_str;
    gib.cmds = _brep_geo_cmds;

    const ON_Brep *brep = ((struct rt_brep_internal *)(gb->intern.idb_ptr))->brep;
    if (brep == NULL) {
	bu_vls_printf(gib.vls, "not a brep object\n");
	return BRLCAD_ERROR;
    }

    if (!argc) {
	_brep_geo_help(&gib, 0, NULL);
	return BRLCAD_OK;
    }

    // TODO: add help flag
    if (argc > 1 && BU_STR_EQUAL(argv[1], HELPFLAG)) {
	argc--;argv++;
	argc--;argv++;
	_brep_geo_help(&gib, argc, argv);
	return BRLCAD_OK;
    }

    // Must have valid subcommand to process
    if (bu_cmd_valid(_brep_geo_cmds, argv[0]) != BRLCAD_OK) {
	bu_vls_printf(gib.vls, "invalid subcommand \"%s\" specified\n", argv[0]);
	_brep_geo_help(&gib, 0, NULL);
	return BRLCAD_ERROR;
    }

    int ret;
    if (bu_cmd(_brep_geo_cmds, argc, argv, 0, (void *)&gib, &ret) == BRLCAD_OK) {
	return ret;
    }
    return BRLCAD_ERROR;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
