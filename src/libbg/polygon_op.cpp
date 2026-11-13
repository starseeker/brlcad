/*                 P O L Y G O N _ O P . C P P
 * BRL-CAD
 *
 * Copyright (c) 2020-2026 United States Government as represented by
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
/** @file polygon_op.cpp
 *
 * Routines for operating on polygons.
 *
 */

#include "common.h"

#include "clipper.hpp"

#include "vmath.h"
#include "bu/log.h"
#include "bu/malloc.h"
#include "bn/mat.h"
#include "bg/plane.h"
#include "bg/defines.h"
extern "C" {
#include "bg/polygon.h"
}

fastf_t
bg_polygon_area(const struct bg_polygon *polygon, fastf_t coordinate_scale,
	const plane_t plane, fastf_t view_scale)
{
    ClipperLib::Path poly;
    fastf_t area = 0.0;

    if (!polygon || (polygon->num_contours && !polygon->contour) ||
	NEAR_ZERO(coordinate_scale, SMALL_FASTF))
	return 0.0;

    plane_t local_plane = HINIT_ZERO;
    if (plane)
	HMOVE(local_plane, plane);
    for (size_t j = 0; j < polygon->num_contours; ++j) {
	const size_t n = polygon->contour[j].num_points;
	if (n && !polygon->contour[j].point)
	    return 0.0;
	poly.resize(n);
	for (size_t k = 0; k < n; k++) {
	    fastf_t fx = polygon->contour[j].point[k][X];
	    fastf_t fy = polygon->contour[j].point[k][Y];
	    if (plane)
		(void)bg_plane_closest_pt(&fx, &fy, &local_plane,
		    (point_t *)&polygon->contour[j].point[k]);

	    poly[k].X = (ClipperLib::long64)(fx * coordinate_scale);
	    poly[k].Y = (ClipperLib::long64)(fy * coordinate_scale);
	}

	const fastf_t contour_area = fabs((fastf_t)ClipperLib::Area(poly));
	area += (polygon->hole && polygon->hole[j]) ?
	    -contour_area : contour_area;
    }

    const fastf_t area_scale = 1.0 /
	(coordinate_scale * coordinate_scale) * view_scale * view_scale;

    return area * area_scale;
}


extern "C" int
bg_polygon_overlaps(const struct bg_polygon *a,
	const struct bg_polygon *b, const plane_t plane,
	const struct bn_tol *tol, fastf_t coordinate_scale)
{
    if (!a || !b || !tol ||
	NEAR_ZERO(coordinate_scale, SMALL_FASTF))
	return 0;

    /* Use the boolean topology engine for contour nesting and holes.  The
     * former bespoke segment/containment path handled only one nesting level
     * and disagreed with polygon booleans for islands inside holes. */
    struct bg_polygon intersection = BG_POLYGON_INIT_ZERO;
    if (bg_polygon_boolean(&intersection,
	    BG_POLYGON_BOOLEAN_INTERSECTION, a, b, coordinate_scale, plane))
	return 0;

    const fastf_t area = fabs(bg_polygon_area(&intersection,
	coordinate_scale, plane, 1.0));
    bg_polygon_clear(&intersection);
    return area > tol->dist_sq;
}


static int
load_polygon(ClipperLib::Clipper &clipper, ClipperLib::PolyType ptype,
	const struct bg_polygon *polygon, fastf_t coordinate_scale,
	const plane_t plane)
{
    ClipperLib::Path curr_poly;
    plane_t local_plane = HINIT_ZERO;
    if (!polygon || (polygon->num_contours && !polygon->contour))
	return 1;
    if (plane)
	HMOVE(local_plane, plane);

    for (size_t j = 0; j < polygon->num_contours; ++j) {
	const struct bg_poly_contour *contour = &polygon->contour[j];
	const size_t n = contour->num_points;
	if (!n)
	    continue;
	if (!contour->point)
	    return 1;
	curr_poly.resize(n);
	for (size_t k = 0; k < n; k++) {
	    fastf_t fx = contour->point[k][X];
	    fastf_t fy = contour->point[k][Y];
	    if (plane)
		bg_plane_closest_pt(&fx, &fy, &local_plane,
		    (point_t *)&contour->point[k]);

	    curr_poly[k].X = (ClipperLib::long64)(fx * coordinate_scale);
	    curr_poly[k].Y = (ClipperLib::long64)(fy * coordinate_scale);
	}

	try {
	    if (!clipper.AddPath(curr_poly, ptype, !contour->open))
		return 1;
	} catch (...) {
	    return 1;
	}
    }

    return 0;
}

static int
load_polygons(ClipperLib::Clipper &clipper, ClipperLib::PolyType ptype,
	const struct bg_polygons *polygons, fastf_t coordinate_scale,
	const plane_t plane)
{
    if (!polygons || (polygons->num_polygons && !polygons->polygon))
	return 1;
    for (size_t i = 0; i < polygons->num_polygons; ++i) {
	if (load_polygon(clipper, ptype, &polygons->polygon[i],
		coordinate_scale, plane))
	    return 1;
    }

    return 0;
}

static int
extract_polygon(struct bg_polygon *result,
	ClipperLib::PolyTree &clipper_polytree, fastf_t inverse_scale,
	const plane_t plane)
{
    if (!result)
	return 1;
    struct bg_polygon extracted = BG_POLYGON_INIT_ZERO;
    const size_t num_contours = clipper_polytree.Total();
    if (!num_contours)
	return bg_polygon_move(result, &extracted);

    extracted.num_contours = num_contours;
    extracted.hole = (int *)bu_calloc(num_contours, sizeof(int), "hole");
    extracted.contour = (struct bg_poly_contour *)bu_calloc(num_contours,
	    sizeof(struct bg_poly_contour), "contour");
    plane_t local_plane = HINIT_ZERO;
    if (plane)
	HMOVE(local_plane, plane);

    ClipperLib::PolyNode *polynode = clipper_polytree.GetFirst();
    size_t n = 0;
    while (polynode) {
	ClipperLib::Path &path = polynode->Contour;

	extracted.hole[n] = polynode->IsHole();
	extracted.contour[n].num_points = path.size();
	extracted.contour[n].open = polynode->IsOpen();
	extracted.contour[n].point = (point_t *)bu_calloc(path.size(),
	    sizeof(point_t), "point");

	for (size_t j = 0; j < path.size(); ++j) {
	    const fastf_t x = (fastf_t)path[j].X * inverse_scale;
	    const fastf_t y = (fastf_t)path[j].Y * inverse_scale;
	    if (plane) {
		bg_plane_pt_at(&extracted.contour[n].point[j],
		    &local_plane, x, y);
	    } else {
		VSET(extracted.contour[n].point[j], x, y, 0);
	    }
	}

	++n;
	polynode = polynode->GetNext();
    }

    extracted.num_contours = n;
    return bg_polygon_move(result, &extracted);
}

static int
clipper_operation(ClipperLib::ClipType *clipper_op,
	enum bg_polygon_boolean_op op)
{
    switch (op) {
    case BG_POLYGON_BOOLEAN_INTERSECTION:
	*clipper_op = ClipperLib::ctIntersection;
	return 0;
    case BG_POLYGON_BOOLEAN_UNION:
	*clipper_op = ClipperLib::ctUnion;
	return 0;
    case BG_POLYGON_BOOLEAN_DIFFERENCE:
	*clipper_op = ClipperLib::ctDifference;
	return 0;
    case BG_POLYGON_BOOLEAN_XOR:
	*clipper_op = ClipperLib::ctXor;
	return 0;
    default:
	return 1;
    }
}

int
bg_polygon_boolean(struct bg_polygon *result,
	enum bg_polygon_boolean_op op, const struct bg_polygon *subject,
	const struct bg_polygon *clip, fastf_t coordinate_scale,
	const plane_t plane)
{
    if (!result || !subject || !clip ||
	NEAR_ZERO(coordinate_scale, SMALL_FASTF))
	return 1;
    ClipperLib::ClipType clipper_op;
    if (clipper_operation(&clipper_op, op))
	return 1;

    ClipperLib::Clipper clipper_engine;
    ClipperLib::PolyTree clipped;
    if (load_polygon(clipper_engine, ClipperLib::ptSubject, subject,
	    coordinate_scale, plane) ||
	load_polygon(clipper_engine, ClipperLib::ptClip, clip,
	    coordinate_scale, plane))
	return 1;
    try {
	if (!clipper_engine.Execute(clipper_op, clipped,
		ClipperLib::pftEvenOdd, ClipperLib::pftEvenOdd))
	    return 1;
    } catch (...) {
	return 1;
    }
    return extract_polygon(result, clipped, 1.0 / coordinate_scale, plane);
}

int
bg_polygons_boolean(struct bg_polygon *result,
	enum bg_polygon_boolean_op op, const struct bg_polygons *subject,
	const struct bg_polygons *clip, fastf_t coordinate_scale,
	const plane_t plane)
{
    if (!result || !subject || !clip ||
	NEAR_ZERO(coordinate_scale, SMALL_FASTF))
	return 1;
    ClipperLib::ClipType clipper_op;
    if (clipper_operation(&clipper_op, op))
	return 1;

    ClipperLib::Clipper clipper_engine;
    ClipperLib::PolyTree clipped;
    if (load_polygons(clipper_engine, ClipperLib::ptSubject, subject,
	    coordinate_scale, plane) ||
	load_polygons(clipper_engine, ClipperLib::ptClip, clip,
	    coordinate_scale, plane))
	return 1;
    try {
	if (!clipper_engine.Execute(clipper_op, clipped,
		ClipperLib::pftEvenOdd, ClipperLib::pftEvenOdd))
	    return 1;
    } catch (...) {
	return 1;
    }
    return extract_polygon(result, clipped, 1.0 / coordinate_scale, plane);
}

int
bg_polygon_hatch(struct bg_polygon *result,
	const struct bg_polygon *polygon, const plane_t plane,
	const vect2d_t slope, fastf_t spacing)
{
    if (!result || !polygon || !plane || !slope ||
	polygon->num_contours < 1 || !polygon->contour ||
	fabs(spacing) < BN_TOL_DIST)
	return 1;

    plane_t local_plane;
    HMOVE(local_plane, plane);
    vect2d_t direction;
    V2MOVE(direction, slope);
    if (MAG2SQ(direction) < SMALL_FASTF)
	V2SET(direction, 1.0, 0.0);
    V2UNITIZE(direction);
    if (direction[X] < 0.0 ||
	(NEAR_ZERO(direction[X], SMALL_FASTF) && direction[Y] < 0.0)) {
	direction[X] = -direction[X];
	direction[Y] = -direction[Y];
    }

    struct bg_polygon mask = BG_POLYGON_INIT_ZERO;
    if (bg_polygon_copy(&mask, polygon))
	return 1;
    vect2d_t minimum = {MAX_FASTF, MAX_FASTF};
    vect2d_t maximum = {-MAX_FASTF, -MAX_FASTF};
    size_t projected_count = 0;
    for (size_t i = 0; i < mask.num_contours; i++) {
	for (size_t j = 0; j < mask.contour[i].num_points; j++) {
	    fastf_t x = 0.0, y = 0.0;
	    bg_plane_closest_pt(&x, &y, &local_plane,
		&mask.contour[i].point[j]);
	    VSET(mask.contour[i].point[j], x, y, 0.0);
	    vect2d_t projected = {x, y};
	    V2MINMAX(minimum, maximum, projected);
	    projected_count++;
	}
    }
    if (!projected_count) {
	bg_polygon_clear(&mask);
	return 1;
    }

    const fastf_t diagonal = DIST_PNT2_PNT2(maximum, minimum);
    const fastf_t line_spacing = fabs(spacing);
    const size_t each_side = (size_t)ceil(0.5 * diagonal / line_spacing);
    const size_t line_count = each_side * 2 + 1;
    struct bg_polygon lines = BG_POLYGON_INIT_ZERO;
    lines.num_contours = line_count;
    lines.hole = (int *)bu_calloc(line_count, sizeof(int), "hatch holes");
    lines.contour = (struct bg_poly_contour *)bu_calloc(line_count,
	    sizeof(struct bg_poly_contour), "hatch contours");

    vect2d_t center, perpendicular;
    center[X] = 0.5 * (minimum[X] + maximum[X]);
    center[Y] = 0.5 * (minimum[Y] + maximum[Y]);
    V2SET(perpendicular, -direction[Y], direction[X]);
    const fastf_t half_length = 0.55 * diagonal + line_spacing;
    for (size_t i = 0; i < line_count; i++) {
	const fastf_t offset = ((fastf_t)i - (fastf_t)each_side) *
	    line_spacing;
	vect2d_t line_center, start, end;
	V2JOIN1(line_center, center, offset, perpendicular);
	V2JOIN1(start, line_center, -half_length, direction);
	V2JOIN1(end, line_center, half_length, direction);
	lines.contour[i].num_points = 2;
	lines.contour[i].open = 1;
	lines.contour[i].point = (point_t *)bu_calloc(2, sizeof(point_t),
	    "hatch line points");
	VSET(lines.contour[i].point[0], start[X], start[Y], 0.0);
	VSET(lines.contour[i].point[1], end[X], end[Y], 0.0);
    }

    struct bg_polygon clipped = BG_POLYGON_INIT_ZERO;
    const int boolean_status = bg_polygon_boolean(&clipped,
	BG_POLYGON_BOOLEAN_INTERSECTION, &lines, &mask, CLIPPER_MAX, NULL);
    bg_polygon_clear(&lines);
    bg_polygon_clear(&mask);
    if (boolean_status) {
	bg_polygon_clear(&clipped);
	return 1;
    }

    for (size_t i = 0; i < clipped.num_contours; i++) {
	clipped.contour[i].open = 1;
	for (size_t j = 0; j < clipped.contour[i].num_points; j++) {
	    const fastf_t x = clipped.contour[i].point[j][X];
	    const fastf_t y = clipped.contour[i].point[j][Y];
	    bg_plane_pt_at(&clipped.contour[i].point[j], &local_plane, x, y);
	}
    }
    return bg_polygon_move(result, &clipped);
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
