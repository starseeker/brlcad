/*                    P O L Y _ S K E T C H . C
 * BRL-CAD
 *
 * Copyright (c) 2013-2026 United States Government as represented by
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

#include "common.h"

#include <string.h>

#include "vmath.h"
#include "bu/app.h"
#include "bu/file.h"
#include "bg/plane.h"
#include "bg/polygon.h"
#include "raytrace.h"
#include "rt/primitives/sketch.h"

#include "../primitives/sketch/polygons_private.h"

/* Compare model-space points using the standard BN distance tolerance. */
static int
points_within_tolerance(point_t a, point_t b)
{
    return NEAR_ZERO(DIST_PNT_PNT(a, b), BN_TOL_DIST);
}

static void
compare_polygon_data(const struct bg_polygon *op, const struct bg_polygon *rp, const char *msg)
{
    if (!op || !rp)
	bu_exit(EXIT_FAILURE, "%s: missing polygon data\n", msg);

    if (op->num_contours != rp->num_contours)
	bu_exit(EXIT_FAILURE, "%s: contour count changed\n", msg);

    for (size_t i = 0; i < op->num_contours; i++) {
	if ((op->hole && op->hole[i] ? 1 : 0) !=
	    (rp->hole && rp->hole[i] ? 1 : 0))
	    bu_exit(EXIT_FAILURE, "%s: contour %zu hole state changed\n", msg, i);
	if ((op->contour[i].open ? 1 : 0) !=
	    (rp->contour[i].open ? 1 : 0))
	    bu_exit(EXIT_FAILURE, "%s: contour %zu open state changed\n", msg, i);
	size_t pcnt = op->contour[i].num_points;
	if (pcnt != rp->contour[i].num_points)
	    bu_exit(EXIT_FAILURE, "%s: contour point count changed\n", msg);

	int contour_match = 0;
	for (size_t offset = 0; offset < pcnt; offset++) {
	    int match = 1;
	    for (size_t j = 0; j < pcnt; j++) {
		size_t rj = (offset + j) % pcnt;
		if (!points_within_tolerance(op->contour[i].point[j], rp->contour[i].point[rj])) {
		    match = 0;
		    break;
		}
	    }
	    if (match) {
		contour_match = 1;
		break;
	    }

	    match = 1;
	    for (size_t j = 0; j < pcnt; j++) {
		size_t rj = (offset + pcnt - j) % pcnt;
		if (!points_within_tolerance(op->contour[i].point[j], rp->contour[i].point[rj])) {
		    match = 0;
		    break;
		}
	    }
	    if (match) {
		contour_match = 1;
		break;
	    }
	}

	if (!contour_match)
	    bu_exit(EXIT_FAILURE, "%s: contour %zu points moved\n", msg, i);
    }
}

static void
compare_rt_polygon_data(const struct rt_sketch_polygon *orig, const struct rt_sketch_polygon *rt, const char *msg)
{
    compare_polygon_data(rt_sketch_polygon_bg_polygon(orig),
	    rt_sketch_polygon_bg_polygon(rt), msg);
}

static void
test_non_origin_plane_roundtrip(void)
{
    char ofile[MAXPATHLEN];
    bu_dir(ofile, MAXPATHLEN, BU_DIR_CURR, "poly_sketch_roundtrip.g", NULL);
    struct rt_wdb *wfp = wdb_fopen_v(ofile, 5);
    if (!wfp)
	bu_exit(EXIT_FAILURE, "Failed to create output database %s\n", ofile);

    struct rt_sketch_polygon p;
    memset(&p, 0, sizeof(p));
    p.type = RT_SKETCH_POLYGON_GENERAL;
    p.curr_contour_i = -1;
    p.curr_point_i = -1;
    p.polygon.num_contours = 1;
    p.polygon.hole = (int *)bu_calloc(1, sizeof(int), "gp_hole");
    p.polygon.contour = (struct bg_poly_contour *)bu_calloc(1, sizeof(struct bg_poly_contour), "gp_contour");
    p.polygon.contour[0].num_points = 3;
    p.polygon.contour[0].point = (point_t *)bu_calloc(3, sizeof(point_t), "gpc_point");

    VSET(p.polygon.contour[0].point[0], 11.0, 21.0, 30.0);
    VSET(p.polygon.contour[0].point[1], 15.0, 21.0, 30.0);
    VSET(p.polygon.contour[0].point[2], 11.0, 26.0, 30.0);

    /* Z=30 gives a non-zero plane origin, verifying exported axes are relative
     * to the sketch origin rather than world origin. */
    point_t plane_pt = {0.0, 0.0, 30.0};
    vect_t plane_n = {0.0, 0.0, 1.0};
    bg_plane_pt_nrml(&p.vp, plane_pt, plane_n);

    struct directory *odp = db_sketch_polygon_to_sketch(wfp->dbip, "roundtrip.s", &p, NULL);
    if (odp == RT_DIR_NULL)
	bu_exit(EXIT_FAILURE, "Failed to write non-origin polygon to output database %s\n", ofile);

    struct rt_sketch_polygon *rtobj = db_sketch_to_polygon("roundtrip_rt", wfp->dbip, odp);
    if (!rtobj)
	bu_exit(EXIT_FAILURE, "Failed to recreate non-origin polygon from sketch\n");

    compare_rt_polygon_data(&p, rtobj, "non-origin plane polygon roundtrip");

    rt_sketch_polygon_destroy(rtobj);
    bg_polygon_clear(&p.polygon);
    db_close(wfp->dbip);
    bu_file_delete(ofile);
}


static void
test_contour_topology_roundtrip(void)
{
    char ofile[MAXPATHLEN];
    bu_dir(ofile, MAXPATHLEN, BU_DIR_CURR,
	    "poly_sketch_topology.g", NULL);
    struct rt_wdb *wfp = wdb_fopen_v(ofile, 5);
    if (!wfp)
	bu_exit(EXIT_FAILURE, "Failed to create output database %s\n", ofile);

    struct rt_sketch_polygon p;
    memset(&p, 0, sizeof(p));
    p.type = RT_SKETCH_POLYGON_GENERAL;
    p.curr_contour_i = -1;
    p.curr_point_i = -1;
    p.polygon.num_contours = 3;
    p.polygon.hole = (int *)bu_calloc(3, sizeof(int), "gp_hole");
    p.polygon.contour = (struct bg_poly_contour *)bu_calloc(3,
	    sizeof(struct bg_poly_contour), "gp_contour");
    p.polygon.hole[1] = 1;

    p.polygon.contour[0].num_points = 4;
    p.polygon.contour[0].point = (point_t *)bu_calloc(4,
	    sizeof(point_t), "outer points");
    VSET(p.polygon.contour[0].point[0], 0.0, 0.0, 0.0);
    VSET(p.polygon.contour[0].point[1], 10.0, 0.0, 0.0);
    VSET(p.polygon.contour[0].point[2], 10.0, 10.0, 0.0);
    VSET(p.polygon.contour[0].point[3], 0.0, 10.0, 0.0);

    p.polygon.contour[1].num_points = 4;
    p.polygon.contour[1].point = (point_t *)bu_calloc(4,
	    sizeof(point_t), "hole points");
    VSET(p.polygon.contour[1].point[0], 2.0, 2.0, 0.0);
    VSET(p.polygon.contour[1].point[1], 2.0, 4.0, 0.0);
    VSET(p.polygon.contour[1].point[2], 4.0, 4.0, 0.0);
    VSET(p.polygon.contour[1].point[3], 4.0, 2.0, 0.0);

    p.polygon.contour[2].num_points = 3;
    p.polygon.contour[2].open = 1;
    p.polygon.contour[2].point = (point_t *)bu_calloc(3,
	    sizeof(point_t), "open points");
    VSET(p.polygon.contour[2].point[0], 20.0, 0.0, 0.0);
    VSET(p.polygon.contour[2].point[1], 22.0, 1.0, 0.0);
    VSET(p.polygon.contour[2].point[2], 24.0, 0.0, 0.0);

    point_t plane_pt = VINIT_ZERO;
    vect_t plane_n = {0.0, 0.0, 1.0};
    bg_plane_pt_nrml(&p.vp, plane_pt, plane_n);

    struct directory *dp = db_sketch_polygon_to_sketch(wfp->dbip,
	    "topology.s", &p, NULL);
    if (dp == RT_DIR_NULL)
	bu_exit(EXIT_FAILURE, "Failed to write topology polygon\n");
    struct rt_sketch_polygon *roundtrip = db_sketch_to_polygon(
	    "topology_rt", wfp->dbip, dp);
    if (!roundtrip)
	bu_exit(EXIT_FAILURE, "Failed to read topology polygon\n");
    compare_rt_polygon_data(&p, roundtrip,
	    "hole and open contour roundtrip");

    rt_sketch_polygon_destroy(roundtrip);
    bg_polygon_clear(&p.polygon);
    db_close(wfp->dbip);
    bu_file_delete(ofile);
}


static void
write_test_poly_database(const char *ofile)
{
    struct rt_wdb *wfp = wdb_fopen_v(ofile, 5);
    if (!wfp)
	bu_exit(EXIT_FAILURE, "Failed to create input database %s\n", ofile);

    struct rt_sketch_polygon p;
    memset(&p, 0, sizeof(p));
    p.type = RT_SKETCH_POLYGON_GENERAL;
    p.curr_contour_i = -1;
    p.curr_point_i = -1;
    p.polygon.num_contours = 1;
    p.polygon.hole = (int *)bu_calloc(1, sizeof(int), "gp_hole");
    p.polygon.contour = (struct bg_poly_contour *)bu_calloc(1, sizeof(struct bg_poly_contour), "gp_contour");
    p.polygon.contour[0].num_points = 4;
    p.polygon.contour[0].point = (point_t *)bu_calloc(4, sizeof(point_t), "gpc_point");

    VSET(p.polygon.contour[0].point[0], 0.0, 0.0, 0.0);
    VSET(p.polygon.contour[0].point[1], 2.0, 0.0, 0.0);
    VSET(p.polygon.contour[0].point[2], 2.0, 2.0, 0.0);
    VSET(p.polygon.contour[0].point[3], 0.0, 2.0, 0.0);

    point_t plane_pt = {0.0, 0.0, 0.0};
    vect_t plane_n = {0.0, 0.0, 1.0};
    bg_plane_pt_nrml(&p.vp, plane_pt, plane_n);

    struct directory *odp = db_sketch_polygon_to_sketch(wfp->dbip, "poly.s", &p, NULL);
    if (odp == RT_DIR_NULL)
	bu_exit(EXIT_FAILURE, "Failed to write polygon to input database %s\n", ofile);

    bg_polygon_clear(&p.polygon);
    db_close(wfp->dbip);
}


int
main(int argc, char *argv[])
{
    struct db_i *dbip;
    struct directory *dp;
    const char *input_file = NULL;
    char generated_file[MAXPATHLEN] = {0};
    int generated_input = 0;

    bu_setprogname(argv[0]);

    if (argc == 1) {
	bu_dir(generated_file, MAXPATHLEN, BU_DIR_CURR,
		"poly_sketch_input.g", NULL);
	write_test_poly_database(generated_file);
	input_file = generated_file;
	generated_input = 1;
    } else if (argc == 2) {
	input_file = argv[1];
    } else {
	bu_exit(EXIT_FAILURE, "Usage: %s [file.g]", argv[0]);
    }

    // First, open the database and make sure we have poly.s

    dbip = db_open(input_file, DB_OPEN_READONLY);
    if (dbip == DBI_NULL) {
	bu_exit(EXIT_FAILURE, "ERROR: Unable to read from %s\n", input_file);
    }

    if (db_dirbuild(dbip) < 0) {
	bu_exit(EXIT_FAILURE, "ERROR: Unable to read from %s\n", input_file);
    }

    db_update_nref(dbip);

    dp = db_lookup(dbip, "poly.s", LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
	bu_exit(EXIT_FAILURE, "ERROR: Unable to look up object poly.s\n");

    struct rt_sketch_polygon *pobj = db_sketch_to_polygon("poly", dbip, dp);
    if (!pobj)
	bu_exit(EXIT_FAILURE, "Failed to create polygon from poly.s\n");

    char ofile[MAXPATHLEN];
    bu_dir(ofile, MAXPATHLEN, BU_DIR_CURR, "poly_sketch_out.g", NULL);
    struct rt_wdb *wfp = wdb_fopen_v(ofile, 5);
    if (!wfp)
	bu_exit(EXIT_FAILURE, "Failed to create output database %s\n", ofile);

    struct directory *odp = db_sketch_polygon_to_sketch(wfp->dbip, "poly.s", pobj, NULL);
    if (odp == RT_DIR_NULL)
	bu_exit(EXIT_FAILURE, "Failed to write polygon to output database %s\n", ofile);

    struct rt_sketch_polygon *opobj = db_sketch_to_polygon("poly_out", wfp->dbip, odp);
    if (!opobj)
	bu_exit(EXIT_FAILURE, "Failed to create polygon from exported poly.s\n");

    compare_rt_polygon_data(pobj, opobj, "imported sketch polygon roundtrip");

    rt_sketch_polygon_destroy(opobj);

    /* Updating a synchronized polygon must replace the existing sketch in
     * place.  Creation deliberately still rejects an existing name. */
    struct rt_sketch_polygon_data updated;
    rt_sketch_polygon_data_init(&updated);
    updated.type = pobj->type;
    updated.fill_flag = pobj->fill_flag;
    V2MOVE(updated.fill_dir, pobj->fill_dir);
    updated.fill_delta = pobj->fill_delta;
    BU_COLOR_CPY(&updated.fill_color, &pobj->fill_color);
    VMOVE(updated.origin_point, pobj->origin_point);
    HMOVE(updated.vp, pobj->vp);
    updated.vZ = pobj->vZ;
    updated.have_edge_color = pobj->have_edge_color;
    BU_COLOR_CPY(&updated.edge_color, &pobj->edge_color);
    if (bg_polygon_copy(&updated.polygon,
	    rt_sketch_polygon_bg_polygon(pobj)))
	bu_exit(EXIT_FAILURE, "Failed to copy synchronized polygon data\n");
    updated.polygon.contour[0].point[1][X] += 3.0;

    if (db_sketch_polygon_data_to_sketch(wfp->dbip, "poly.s", &updated) !=
	    RT_DIR_NULL)
	bu_exit(EXIT_FAILURE, "Sketch creation overwrote an existing object\n");
    struct directory *updated_dp = db_sketch_polygon_data_update_sketch(
	    wfp->dbip, "poly.s", &updated);
    if (updated_dp == RT_DIR_NULL || updated_dp != odp)
	bu_exit(EXIT_FAILURE, "Failed to update existing synchronized sketch\n");
    struct rt_sketch_polygon *updated_obj = db_sketch_to_polygon(
	    "poly_updated", wfp->dbip, updated_dp);
    if (!updated_obj)
	bu_exit(EXIT_FAILURE, "Failed to read updated synchronized sketch\n");
    compare_polygon_data(&updated.polygon,
	    rt_sketch_polygon_bg_polygon(updated_obj),
	    "synchronized sketch update");
    rt_sketch_polygon_destroy(updated_obj);
    rt_sketch_polygon_data_free(&updated);

    rt_sketch_polygon_destroy(pobj);
    db_close(dbip);
    db_close(wfp->dbip);
    bu_file_delete(ofile);
    if (generated_input)
	bu_file_delete(generated_file);

    test_non_origin_plane_roundtrip();
    test_contour_topology_roundtrip();

    return 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
