/*                      P O L Y G O N _ O P . C
 *
 * BRL-CAD
 *
 * Copyright (c) 2015-2026 United States Government as represented by
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
/** @file polygon_op.c
 *
 * Test polygon boolean clipping routines
 *
 */

#include "common.h"

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "bu.h"
#include "vmath.h"
#include "bg.h"
#include "bg/plot3.h"

static int
_bg_polygon_diff(struct bg_polygon *p1, struct bg_polygon *p2)
{
    if (!p1 && !p2) return 0;
    if ((!p1 && p2) || (p1 && !p2)) return 1;
    if (p1->num_contours != p2->num_contours) return 1;
    if (p1->hole && p2->hole) {
	for (size_t i = 0; i < p1->num_contours; i++) {
	    if (p1->hole[i] != p2->hole[i]) {
		return 1;
	    }
	}
    }
    for (size_t i = 0; i < p1->num_contours; i++) {
	struct bg_poly_contour *c1 = &(p1->contour[i]);
	struct bg_poly_contour *c2 = &(p2->contour[i]);

	if (c1->num_points != c2->num_points) {
	    return 1;
	}

	// Clipper may return the points with a different starting point.
	// To handle this, make an initial pass through the points to
	// find the first matching point in p2, and use that offset when
	// checking subsequent points
	size_t offset = 0;
	for (size_t j = 0; j < c1->num_points; j++) {
	    if (DIST_PNT_PNT_SQ(c1->point[0], c2->point[j]) > VUNITIZE_TOL) {
		continue;
	    } else {
		// Distance match - we have our point
		offset = j;
		break;
	    }
	}

	// Have alignment, check the points
	for (size_t j = 0; j < c1->num_points; j++) {
	    size_t p2_ind = ((offset + j) >= c1->num_points) ? (j - offset) : (offset + j);
	    if (DIST_PNT_PNT_SQ(c1->point[j], c2->point[p2_ind]) > VUNITIZE_TOL) {
		return 1;
	    }
	}

    }
    return 0;
}

static int
test_polygon_value_operations(void)
{
    int failures = 0;
    plane_t plane = {0.0, 0.0, 1.0, 0.0};
    point_t first = {0.0, 0.0, 0.0};
    point_t opposite = {4.0, 2.0, 0.0};
    struct bg_polygon rectangle = BG_POLYGON_INIT_ZERO;
    if (bg_polygon_make_rectangle(&rectangle, plane, first, opposite, 0) ||
	rectangle.num_contours != 1 ||
	bg_polygon_point_count(&rectangle) != 4 ||
	rectangle.contour[0].open)
	failures++;

    const fastf_t area = bg_polygon_area(&rectangle, CLIPPER_MAX,
	plane, 1.0);
    if (!NEAR_EQUAL(area, 8.0, 1.0e-6))
	failures++;

    struct bg_polygon copy = BG_POLYGON_INIT_ZERO;
    if (bg_polygon_copy(&copy, &rectangle) ||
	_bg_polygon_diff(&copy, &rectangle))
	failures++;

    vect_t delta = {3.0, -2.0, 1.0};
    point_t translated_first;
    VADD2(translated_first, rectangle.contour[0].point[0], delta);
    if (bg_polygon_translate(&copy, delta))
	failures++;
    size_t contour = SIZE_MAX;
    size_t point = SIZE_MAX;
    fastf_t distance_sq = INFINITY;
    if (bg_polygon_closest_point(&copy, translated_first, &contour,
	    &point, &distance_sq) || contour != 0 || point != 0 ||
	    !NEAR_ZERO(distance_sq, SMALL_FASTF))
	failures++;

    point_t appended = {8.0, 8.0, 8.0};
    if (bg_polygon_append_point(&copy, 0, appended) ||
	bg_polygon_point_count(&copy) != 5 ||
	bg_polygon_remove_point(&copy, 0, 4) ||
	bg_polygon_point_count(&copy) != 4 ||
	bg_polygon_contour_open_set(&copy, 0, 1) ||
	!copy.contour[0].open ||
	bg_polygon_contours_open_set(&copy, 0) || copy.contour[0].open)
	failures++;

    struct bg_polygon removable = BG_POLYGON_INIT_ZERO;
    removable.num_contours = 2;
    removable.hole = (int *)bu_calloc(2, sizeof(int),
	"removable polygon holes");
    removable.contour = (struct bg_poly_contour *)bu_calloc(2,
	sizeof(struct bg_poly_contour), "removable polygon contours");
    removable.hole[1] = 1;
    for (size_t i = 0; i < 2; i++) {
	removable.contour[i].num_points = 1;
	removable.contour[i].point = (point_t *)bu_calloc(1,
	    sizeof(point_t), "removable polygon point");
	VSET(removable.contour[i].point[0], (fastf_t)i + 1.0, 0.0, 0.0);
    }
    if (bg_polygon_remove_point(&removable, 0, 0) ||
	removable.num_contours != 1 || !removable.hole ||
	removable.hole[0] != 1 ||
	!NEAR_EQUAL(removable.contour[0].point[0][X], 2.0, SMALL_FASTF) ||
	bg_polygon_remove_point(&removable, 0, 0) ||
	removable.num_contours || removable.contour || removable.hole)
	failures++;
    bg_polygon_clear(&removable);

    struct bg_polygon ellipse = BG_POLYGON_INIT_ZERO;
    if (bg_polygon_make_ellipse(&ellipse, plane, first, opposite, 0, 32) ||
	bg_polygon_point_count(&ellipse) != 32) {
	bu_log("ellipse construction contract failure: count=%zu\n",
	    bg_polygon_point_count(&ellipse));
	failures++;
	}

    struct bg_polygon circle = BG_POLYGON_INIT_ZERO;
    if (bg_polygon_make_ellipse(&circle, plane, first, opposite, 1, 32) ||
	bg_polygon_point_count(&circle) != 32 ||
	!NEAR_EQUAL(DIST_PNT_PNT(circle.contour[0].point[0], first),
	    sqrt(20.0), 1.0e-6)) {
	bu_log("circle construction contract failure: count=%zu radius=%g\n",
	    bg_polygon_point_count(&circle),
	    DIST_PNT_PNT(circle.contour[0].point[0], first));
	failures++;
	}

    vect2d_t slope = {1.0, 0.0};
    struct bg_polygon hatch = BG_POLYGON_INIT_ZERO;
    if (bg_polygon_hatch(&hatch, &rectangle, plane, slope, 0.5) ||
	!hatch.num_contours || !bg_polygon_point_count(&hatch))
	failures++;
    for (size_t i = 0; i < hatch.num_contours; i++) {
	if (!hatch.contour[i].open)
	    failures++;
    }

    struct bg_polygon moved = BG_POLYGON_INIT_ZERO;
    if (bg_polygon_move(&moved, &ellipse) || ellipse.num_contours ||
	!moved.num_contours)
	failures++;

    /* Triangulation must reject empty/malformed values without leaving stale
     * output pointers, and must treat multiple positive contours as distinct
     * components rather than silently reclassifying all but the first as
     * holes. */
    int *faces = (int *)(uintptr_t)1;
    int face_count = 7;
    point_t *vertices = (point_t *)(uintptr_t)1;
    int vertex_count = 9;
    struct bg_polygon empty = BG_POLYGON_INIT_ZERO;
    if (bg_polygon_triangulate(&faces, &face_count, &vertices,
	    &vertex_count, &empty, TRI_EAR_CLIPPING) == 0 || faces ||
	    vertices || face_count || vertex_count) {
	bu_log("empty polygon triangulation contract failure\n");
	failures++;
	}

    struct bg_polygon malformed = BG_POLYGON_INIT_ZERO;
    malformed.num_contours = 1;
    malformed.contour = (struct bg_poly_contour *)bu_calloc(1,
	    sizeof(struct bg_poly_contour), "malformed polygon contour");
    malformed.contour[0].num_points = 3;
    if (bg_polygon_triangulate(&faces, &face_count, &vertices,
	    &vertex_count, &malformed, TRI_EAR_CLIPPING) == 0 || faces ||
	    vertices || face_count || vertex_count) {
	bu_log("malformed polygon triangulation contract failure\n");
	failures++;
	}
    bg_polygon_clear(&malformed);

    struct bg_polygon components = BG_POLYGON_INIT_ZERO;
    components.num_contours = 2;
    components.hole = (int *)bu_calloc(2, sizeof(int),
	    "component polygon holes");
    components.contour = (struct bg_poly_contour *)bu_calloc(2,
	    sizeof(struct bg_poly_contour), "component polygon contours");
    for (size_t i = 0; i < 2; i++) {
	components.contour[i].num_points = 4;
	components.contour[i].point = (point_t *)bu_calloc(4,
	    sizeof(point_t), "component polygon points");
	const fastf_t x = (fastf_t)i * 4.0;
	VSET(components.contour[i].point[0], x, 0.0, 0.0);
	VSET(components.contour[i].point[1], x + 2.0, 0.0, 0.0);
	VSET(components.contour[i].point[2], x + 2.0, 2.0, 0.0);
	VSET(components.contour[i].point[3], x, 2.0, 0.0);
    }
    if (bg_polygon_triangulate(&faces, &face_count, &vertices,
	    &vertex_count, &components, TRI_EAR_CLIPPING) ||
	    face_count != 4 || vertex_count < 8 || !faces || !vertices) {
	bu_log("multi-component triangulation failure: %d faces, %d vertices\n",
	    face_count, vertex_count);
	failures++;
	}
    if (faces)
	bu_free(faces, "component triangulation faces");
    if (vertices)
	bu_free(vertices, "component triangulation vertices");
    bg_polygon_clear(&components);

    /* Overlap classification must honor arbitrary contour nesting.  This
     * polygon is an outer square, a hole, and a filled island in that hole. */
    struct bg_polygon nested = BG_POLYGON_INIT_ZERO;
    nested.num_contours = 3;
    nested.hole = (int *)bu_calloc(3, sizeof(int), "nested polygon holes");
    nested.contour = (struct bg_poly_contour *)bu_calloc(3,
	    sizeof(struct bg_poly_contour), "nested polygon contours");
    nested.hole[1] = 1;
    const fastf_t square_min[3] = {0.0, 2.0, 4.0};
    const fastf_t square_max[3] = {10.0, 8.0, 6.0};
    for (size_t i = 0; i < 3; i++) {
	nested.contour[i].num_points = 4;
	nested.contour[i].point = (point_t *)bu_calloc(4, sizeof(point_t),
	    "nested polygon points");
	VSET(nested.contour[i].point[0], square_min[i], square_min[i], 0.0);
	VSET(nested.contour[i].point[1], square_max[i], square_min[i], 0.0);
	VSET(nested.contour[i].point[2], square_max[i], square_max[i], 0.0);
	VSET(nested.contour[i].point[3], square_min[i], square_max[i], 0.0);
    }
    struct bn_tol tolerance = BN_TOL_INIT_TOL;
    struct bg_polygon probe = BG_POLYGON_INIT_ZERO;
    point_t probe_min = {4.5, 4.5, 0.0};
    point_t probe_max = {5.5, 5.5, 0.0};
    if (bg_polygon_make_rectangle(&probe, plane, probe_min, probe_max, 0) ||
	!bg_polygon_overlaps(&nested, &probe, plane, &tolerance, CLIPPER_MAX)) {
	bu_log("nested polygon island overlap failure\n");
	failures++;
    }
    bg_polygon_clear(&probe);
    VSET(probe_min, 3.0, 3.0, 0.0);
    VSET(probe_max, 3.5, 3.5, 0.0);
    if (bg_polygon_make_rectangle(&probe, plane, probe_min, probe_max, 0) ||
	bg_polygon_overlaps(&nested, &probe, plane, &tolerance, CLIPPER_MAX)) {
	bu_log("nested polygon hole overlap failure\n");
	failures++;
    }
    bg_polygon_clear(&probe);
    bg_polygon_clear(&nested);

    bg_polygon_clear(&rectangle);
    bg_polygon_clear(&rectangle);
    bg_polygon_clear(&copy);
    bg_polygon_clear(&ellipse);
    bg_polygon_clear(&hatch);
    bg_polygon_clear(&moved);
    bg_polygon_clear(&circle);
    return failures;
}

int
main(int argc, const char **argv)
{
    int ret = 0;

    bu_setprogname(argv[0]);

    int plot_files = 0;
    if (argc > 1 && BU_STR_EQUAL(argv[1], "--plot")) {
	plot_files = 1;
    }

    ret += test_polygon_value_operations();

    /* Polygon 1 */
    struct bg_polygon p1 = BG_POLYGON_INIT_ZERO;
    p1.num_contours = 1;
    p1.contour = (struct bg_poly_contour *)bu_calloc(1, sizeof(struct bg_poly_contour), "c1 points");
    p1.contour[0].num_points = 4;
    p1.contour[0].point = (point_t *)bu_calloc(4, sizeof(point_t), "p1.contour[0] points");
    VSET(p1.contour[0].point[0], -3, -3, 0);
    VSET(p1.contour[0].point[1], 3, -3, 0);
    VSET(p1.contour[0].point[2], 3, 3, 0);
    VSET(p1.contour[0].point[3], -3, 3, 0);
    if (plot_files) {
	bg_polygon_plot("p1.plot3", (const point_t *)p1.contour[0].point, p1.contour[0].num_points, 255, 0, 0);
    }

    /* Polygon 2 */
    struct bg_polygon p2 = BG_POLYGON_INIT_ZERO;
    p2.num_contours = 1;
    p2.contour = (struct bg_poly_contour *)bu_calloc(1, sizeof(struct bg_poly_contour), "c1 points");
    p2.contour[0].num_points = 4;
    p2.contour[0].point = (point_t *)bu_calloc(4, sizeof(point_t), "p2.contour[0] points");
    VSET(p2.contour[0].point[0], 0, 0, 0);
    VSET(p2.contour[0].point[1], 4, 0, 0);
    VSET(p2.contour[0].point[2], 4, 5, 0);
    VSET(p2.contour[0].point[3], 0, 5, 0);
    if (plot_files) {
	bg_polygon_plot("p2.plot3", (const point_t *)p2.contour[0].point, p2.contour[0].num_points, 0, 255, 0);
    }

    /* Union expected result */
    struct bg_polygon union_expected = BG_POLYGON_INIT_ZERO;
    union_expected.num_contours = 1;
    union_expected.contour = (struct bg_poly_contour *)bu_calloc(1, sizeof(struct bg_poly_contour), "c1 points");
    union_expected.contour[0].num_points = 8;
    union_expected.contour[0].point = (point_t *)bu_calloc(8, sizeof(point_t), "union_expected.contour[0] points");
    VSET(union_expected.contour[0].point[0],  3,  0,  0);
    VSET(union_expected.contour[0].point[1],  4,  0,  0);
    VSET(union_expected.contour[0].point[2],  4,  5,  0);
    VSET(union_expected.contour[0].point[3],  0,  5,  0);
    VSET(union_expected.contour[0].point[4],  0,  3,  0);
    VSET(union_expected.contour[0].point[5], -3,  3,  0);
    VSET(union_expected.contour[0].point[6], -3, -3,  0);
    VSET(union_expected.contour[0].point[7],  3, -3,  0);

    /* Calculate union and compare it with the expected result */
    struct bg_polygon ur = BG_POLYGON_INIT_ZERO;
    ret += bg_polygon_boolean(&ur, BG_POLYGON_BOOLEAN_UNION,
	&p1, &p2, 1.0, NULL);
    if (plot_files) {
	bg_polygon_plot("ur.plot3", (const point_t *)ur.contour[0].point,
	    ur.contour[0].num_points, 0, 0, 255);
    }
    ret += _bg_polygon_diff(&ur, &union_expected);


    /* Difference expected result */
    struct bg_polygon difference_expected = BG_POLYGON_INIT_ZERO;
    difference_expected.num_contours = 1;
    difference_expected.contour = (struct bg_poly_contour *)bu_calloc(1, sizeof(struct bg_poly_contour), "c1 points");
    difference_expected.contour[0].num_points = 6;
    difference_expected.contour[0].point = (point_t *)bu_calloc(6, sizeof(point_t), "difference_expected.contour[0] points");
    VSET(difference_expected.contour[0].point[0],  3,  0,  0);
    VSET(difference_expected.contour[0].point[1],  0,  0,  0);
    VSET(difference_expected.contour[0].point[2],  0,  3,  0);
    VSET(difference_expected.contour[0].point[3], -3,  3,  0);
    VSET(difference_expected.contour[0].point[4], -3, -3,  0);
    VSET(difference_expected.contour[0].point[5],  3, -3,  0);

    /* Calculate difference and compare it with the expected result */
    struct bg_polygon dr = BG_POLYGON_INIT_ZERO;
    ret += bg_polygon_boolean(&dr, BG_POLYGON_BOOLEAN_DIFFERENCE,
	&p1, &p2, 1.0, NULL);
    if (plot_files) {
	bg_polygon_plot("dr.plot3", (const point_t *)dr.contour[0].point,
	    dr.contour[0].num_points, 0, 0, 255);
    }
    ret += _bg_polygon_diff(&dr, &difference_expected);

    /* Intersection expected result */
    struct bg_polygon intersection_expected = BG_POLYGON_INIT_ZERO;
    intersection_expected.num_contours = 1;
    intersection_expected.contour = (struct bg_poly_contour *)bu_calloc(1, sizeof(struct bg_poly_contour), "c1 points");
    intersection_expected.contour[0].num_points = 4;
    intersection_expected.contour[0].point = (point_t *)bu_calloc(4, sizeof(point_t), "intersection_expected.contour[0] points");
    VSET(intersection_expected.contour[0].point[0],  3,  3,  0);
    VSET(intersection_expected.contour[0].point[1],  0,  3,  0);
    VSET(intersection_expected.contour[0].point[2],  0,  0,  0);
    VSET(intersection_expected.contour[0].point[3],  3,  0,  0);

    /* Calculate intersection and compare it with the expected result */
    struct bg_polygon ir = BG_POLYGON_INIT_ZERO;
    ret += bg_polygon_boolean(&ir, BG_POLYGON_BOOLEAN_INTERSECTION,
	&p1, &p2, 1.0, NULL);
    if (plot_files) {
	bg_polygon_plot("ir.plot3", (const point_t *)ir.contour[0].point,
	    ir.contour[0].num_points, 0, 0, 255);
    }
    ret += _bg_polygon_diff(&ir, &intersection_expected);

    /* Note - clipper doesn't yet support Xor */


    /* Next, test a polygon full contained within another polygon. */

    /* Polygon 3 */
    struct bg_polygon p3 = BG_POLYGON_INIT_ZERO;
    p3.num_contours = 1;
    p3.contour = (struct bg_poly_contour *)bu_calloc(1, sizeof(struct bg_poly_contour), "c1 points");
    p3.contour[0].num_points = 4;
    p3.contour[0].point = (point_t *)bu_calloc(4, sizeof(point_t), "p3.contour[0] points");
    VSET(p3.contour[0].point[0], -3, -3, 0);
    VSET(p3.contour[0].point[1], 3, -3, 0);
    VSET(p3.contour[0].point[2], 3, 3, 0);
    VSET(p3.contour[0].point[3], -3, 3, 0);
    if (plot_files) {
	bg_polygon_plot("p3.plot3", (const point_t *)p3.contour[0].point, p3.contour[0].num_points, 255, 0, 0);
    }

    /* Polygon 4 */
    struct bg_polygon p4 = BG_POLYGON_INIT_ZERO;
    p4.num_contours = 1;
    p4.contour = (struct bg_poly_contour *)bu_calloc(1, sizeof(struct bg_poly_contour), "c1 points");
    p4.contour[0].num_points = 4;
    p4.contour[0].point = (point_t *)bu_calloc(4, sizeof(point_t), "p4.contour[0] points");
    VSET(p4.contour[0].point[0], -2, -2, 0);
    VSET(p4.contour[0].point[1], 2, -2, 0);
    VSET(p4.contour[0].point[2], 2, 2, 0);
    VSET(p4.contour[0].point[3], -2, 2, 0);
    if (plot_files) {
	bg_polygon_plot("p4.plot3", (const point_t *)p4.contour[0].point, p4.contour[0].num_points, 0, 255, 0);
    }

    /* Calculate union and compare it with the expected result */
    struct bg_polygon ucr = BG_POLYGON_INIT_ZERO;
    ret += bg_polygon_boolean(&ucr, BG_POLYGON_BOOLEAN_UNION,
	&p3, &p4, 1.0, NULL);
    if (plot_files) {
	bg_polygon_plot("ucr.plot3", (const point_t *)ucr.contour[0].point,
	    ucr.contour[0].num_points, 0, 0, 255);
    }
    ret += _bg_polygon_diff(&ucr, &p3);

    /* Calculate intersection and compare it with the expected result */
    struct bg_polygon icr = BG_POLYGON_INIT_ZERO;
    ret += bg_polygon_boolean(&icr, BG_POLYGON_BOOLEAN_INTERSECTION,
	&p3, &p4, 1.0, NULL);
    if (plot_files) {
	bg_polygon_plot("icr.plot3", (const point_t *)icr.contour[0].point,
	    icr.contour[0].num_points, 0, 0, 255);
    }
    ret += _bg_polygon_diff(&icr, &p4);

    /* NOTE: The difference case is the first one that will exercise the creation of a hole - this changes
     * the contour topology of the polygon to multiple-contour. */
    struct bg_polygon de2 = BG_POLYGON_INIT_ZERO;
    de2.num_contours = 2;
    de2.hole = (int *)bu_calloc(de2.num_contours, sizeof(int), "hole");
    de2.contour = (struct bg_poly_contour *)bu_calloc(2, sizeof(struct bg_poly_contour), "points");
    de2.hole[0] = 0;
    de2.contour[0].num_points = 4;
    de2.contour[0].point = (point_t *)bu_calloc(4, sizeof(point_t), "de2.contour[0] points");
    VSET(de2.contour[0].point[0], -3, -3, 0);
    VSET(de2.contour[0].point[1], 3, -3, 0);
    VSET(de2.contour[0].point[2], 3, 3, 0);
    VSET(de2.contour[0].point[3], -3, 3, 0);
    // Compared to p4, the order of the points will be reversed as it is now defining a hole
    de2.hole[1] = 1;
    de2.contour[1].num_points = 4;
    de2.contour[1].point = (point_t *)bu_calloc(4, sizeof(point_t), "de2.contour[1] points");
    VSET(de2.contour[1].point[0], -2, -2, 0);
    VSET(de2.contour[1].point[1], -2, 2, 0);
    VSET(de2.contour[1].point[2], 2, 2, 0);
    VSET(de2.contour[1].point[3], 2, -2, 0);

    /* Calculate difference and compare it with the expected result */
    struct bg_polygon dcr = BG_POLYGON_INIT_ZERO;
    ret += bg_polygon_boolean(&dcr, BG_POLYGON_BOOLEAN_DIFFERENCE,
	&p3, &p4, 1.0, NULL);
    if (plot_files) {
	bg_polygon_plot("dcr_1.plot3", (const point_t *)dcr.contour[0].point,
	    dcr.contour[0].num_points, 0, 0, 255);
	bg_polygon_plot("dcr_2.plot3", (const point_t *)dcr.contour[1].point,
	    dcr.contour[0].num_points, 0, 0, 255);
    }
    ret += _bg_polygon_diff(&dcr, &de2);

    /* Calculate difference the opposite way - this should be a null return, as
     * p4 is inside p3 so subtracting p3 from it leaves no geometry */
    struct bg_polygon dcr2 = BG_POLYGON_INIT_ZERO;
    ret += bg_polygon_boolean(&dcr2, BG_POLYGON_BOOLEAN_DIFFERENCE,
	&p4, &p3, 1.0, NULL);
    ret += dcr2.num_contours ? 1 : 0;

    bg_polygon_clear(&p1);
    bg_polygon_clear(&p2);
    bg_polygon_clear(&union_expected);
    bg_polygon_clear(&ur);
    bg_polygon_clear(&difference_expected);
    bg_polygon_clear(&dr);
    bg_polygon_clear(&intersection_expected);
    bg_polygon_clear(&ir);
    bg_polygon_clear(&p3);
    bg_polygon_clear(&p4);
    bg_polygon_clear(&ucr);
    bg_polygon_clear(&icr);
    bg_polygon_clear(&de2);
    bg_polygon_clear(&dcr);
    bg_polygon_clear(&dcr2);

    return ret;
}


/** @} */
/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
