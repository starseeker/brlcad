/*                     P O L Y G O N . C
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

#include <bio.h>
#include <string.h>

#include "bu/malloc.h"
#include "bu/sort.h"
#include "bg/plane.h"
#ifndef PLOT_PREFIX_STR
#  define PLOT_PREFIX_STR bg_plot3_
#endif
#include "bg/plot3.h"
#include "bn/tol.h"
#include "bg/polygon.h"

void
bg_polygon_init(struct bg_polygon *polygon)
{
    if (!polygon)
	return;
    polygon->num_contours = 0;
    polygon->hole = NULL;
    polygon->contour = NULL;
}

void
bg_polygon_clear(struct bg_polygon *polygon)
{
    if (!polygon)
	return;

    if (polygon->contour) {
	for (size_t j = 0; j < polygon->num_contours; ++j) {
	    if (polygon->contour[j].point)
		bu_free((void *)polygon->contour[j].point, "contour points");
	}
	bu_free((void *)polygon->contour, "contour");
    }
    if (polygon->hole)
	bu_free((void *)polygon->hole, "hole");
    bg_polygon_init(polygon);
}

void
bg_polygons_init(struct bg_polygons *polygons)
{
    if (!polygons)
	return;
    polygons->num_polygons = 0;
    polygons->polygon = NULL;
}

void
bg_polygons_clear(struct bg_polygons *polygons)
{
    if (!polygons)
	return;


    if (polygons->polygon) {
	for (size_t i = 0; i < polygons->num_polygons; ++i)
	    bg_polygon_clear(&polygons->polygon[i]);
	bu_free((void *)polygons->polygon, "polygons");
    }
    bg_polygons_init(polygons);
}

void
bg_polygon_view_bbox(point2d_t *bmin, point2d_t *bmax,
	const struct bg_polygon *p, const mat_t model2view)
{
    if (!bmin || !bmax || !p)
	return;

    // Initialize
    V2SET(*bmin, INFINITY, INFINITY);
    V2SET(*bmax, -INFINITY, -INFINITY);

    if (!p->num_contours || !p->contour)
	return;

    // NOTE:  Holes don't define positive area, so their points are not
    // considered for the bbox dimensions even if they are outside the positive
    // contours.  ONLY considering positive contour points.
    for (size_t i = 0; i < p->num_contours; i++) {
	const struct bg_poly_contour *c = &p->contour[i];
	if (!c->num_points)
	    continue;
	for (size_t j = 0; j < c->num_points; j++) {
	    point_t vpoint;
	    MAT4X3PNT(vpoint, model2view, c->point[j]);
	    point2d_t v2d;
	    v2d[0] = vpoint[0];
	    v2d[1] = vpoint[1];
	    V2MINMAX(*bmin, *bmax, v2d);
	}
    }
}

int
bg_3d_polygon_area(fastf_t *area, size_t npts, const point_t *pts)
{
    size_t i;
    vect_t v1, v2, tmp, tot = VINIT_ZERO;
    plane_t plane_eqn;
    struct bn_tol tol;

    if (!pts || !area || npts < 3)
	return 1;
    BN_TOL_INIT(&tol);
    tol.dist_sq = BN_TOL_DIST * BN_TOL_DIST;
    if (bg_make_plane_3pnts(plane_eqn, pts[0], pts[1], pts[2], &tol) == -1)
	return 1;

    switch (npts) {
	case 3:
	    /* Triangular Face - for triangular face T:V0, V1, V2,
	     * area = 0.5 * [(V2 - V0) x (V1 - V0)] */
	    VSUB2(v1, pts[1], pts[0]);
	    VSUB2(v2, pts[2], pts[0]);
	    VCROSS(tot, v2, v1);
	    break;
	case 4:
	    /* Quadrilateral Face - for planar quadrilateral
	     * Q:V0, V1, V2, V3 with unit normal N,
	     * area = N/2 ⋅ [(V2 - V0) x (V3 - V1)] */
	    VSUB2(v1, pts[2], pts[0]);
	    VSUB2(v2, pts[3], pts[1]);
	    VCROSS(tot, v2, v1);
	    break;
	default:
	    /* N-Sided Face - compute area using Green's Theorem */
	    for (i = 0; i < npts; i++) {
		VCROSS(tmp, pts[i], pts[i + 1 == npts ? 0 : i + 1]);
		VADD2(tot, tot, tmp);
	    }
	    break;
    }
    *area = fabs(VDOT(plane_eqn, tot)) * 0.5;
    return 0;
}


int
bg_3d_polygon_centroid(point_t *cent, size_t npts, const point_t *pts)
{
    size_t i;
    fastf_t x_0 = 0.0;
    fastf_t x_1 = 0.0;
    fastf_t y_0 = 0.0;
    fastf_t y_1 = 0.0;
    fastf_t z_0 = 0.0;
    fastf_t z_1 = 0.0;
    fastf_t a = 0.0;
    fastf_t signedArea = 0.0;

    if (!pts || !cent || npts < 3)
	return 1;
    /* Calculate Centroid projection for face for x-y-plane */
    for (i = 0; i < npts-1; i++) {
	x_0 = pts[i][0];
	y_0 = pts[i][1];
	x_1 = pts[i+1][0];
	y_1 = pts[i+1][1];
	a = x_0 *y_1 - x_1*y_0;
	signedArea += a;
	*cent[0] += (x_0 + x_1)*a;
	*cent[1] += (y_0 + y_1)*a;
    }
    x_0 = pts[i][0];
    y_0 = pts[i][1];
    x_1 = pts[0][0];
    y_1 = pts[0][1];
    a = x_0 *y_1 - x_1*y_0;
    signedArea += a;
    *cent[0] += (x_0 + x_1)*a;
    *cent[1] += (y_0 + y_1)*a;

    signedArea *= 0.5;
    *cent[0] /= (6.0*signedArea);
    *cent[1] /= (6.0*signedArea);

    /* calculate Centroid projection for face for x-z-plane */

    signedArea = 0.0;
    for (i = 0; i < npts-1; i++) {
	x_0 = pts[i][0];
	z_0 = pts[i][2];
	x_1 = pts[i+1][0];
	z_1 = pts[i+1][2];
	a = x_0 *z_1 - x_1*z_0;
	signedArea += a;
	*cent[2] += (z_0 + z_1)*a;
    }
    x_0 = pts[i][0];
    z_0 = pts[i][2];
    x_1 = pts[0][0];
    z_1 = pts[0][2];
    a = x_0 *z_1 - x_1*z_0;
    signedArea += a;
    *cent[2] += (z_0 + z_1)*a;

    signedArea *= 0.5;
    *cent[2] /= (6.0*signedArea);
    return 0;
}


int
bg_3d_polygon_make_pnts_planes(size_t *npts, point_t **pts, size_t neqs, const plane_t *eqs)
{
    size_t i, j, k, l;
    if (!npts || !pts || neqs < 4 || !eqs)
	return 1;
    /* find all vertices */
    for (i = 0; i < neqs - 2; i++) {
	for (j = i + 1; j < neqs - 1; j++) {
	    for (k = j + 1; k < neqs; k++) {
		point_t pt;
		int keep_point = 1;
		if (bg_make_pnt_3planes(pt, eqs[i], eqs[j], eqs[k]) < 0)
		    continue;
		/* discard pt if it is outside the polyhedron */
		for (l = 0; l < neqs; l++) {
		    if (l == i || l == j || l == k)
			continue;
		    if (DIST_PNT_PLANE(pt, eqs[l]) > BN_TOL_DIST) {
			keep_point = 0;
			break;
		    }
		}
		/* found a good point, add it to each of the intersecting faces */
		if (keep_point) {
		    VMOVE(pts[i][npts[i]], (pt)); npts[i]++;
		    VMOVE(pts[j][npts[j]], (pt)); npts[j]++;
		    VMOVE(pts[k][npts[k]], (pt)); npts[k]++;
		}
	    }
	}
    }
    return 0;
}


static int
sort_ccw_3d(const void *x, const void *y, void *cmp)
{
    vect_t tmp;
    VCROSS(tmp, ((fastf_t *)x), ((fastf_t *)y));
    return VDOT(*((point_t *)cmp), tmp);
}


int
bg_3d_polygon_sort_ccw(size_t npts, point_t *pts, plane_t cmp)
{
    size_t i;
    point_t centroid;

    if (!pts || npts < 3)
	return 1;

    /* Compute the centroid of all points.  The sort_ccw_3d comparator
     * measures angles using cross products of the raw position vectors,
     * which is equivalent to sorting by angle around the *origin*.  For
     * faces not centred at the origin the resulting order can be wrong
     * (self-intersecting polygon), so translate all points to be centred
     * at the origin before sorting and translate back afterwards.        */
    VSETALL(centroid, 0.0);
    for (i = 0; i < npts; i++)
	VADD2(centroid, centroid, pts[i]);
    VSCALE(centroid, centroid, 1.0 / (double)npts);

    for (i = 0; i < npts; i++)
	VSUB2(pts[i], pts[i], centroid);

    bu_sort(pts, npts, sizeof(point_t), sort_ccw_3d, &cmp);

    for (i = 0; i < npts; i++)
	VADD2(pts[i], pts[i], centroid);

    return 0;
}

int
bg_polygon_direction(size_t npts, const point2d_t *pts, const int *pt_indices)
{
    size_t i;
    double sum = 0;
    const int *pt_order = NULL;
    int *tmp_pt_order = NULL;
    /* If no array of indices into pts is supplied, construct a
     * temporary version based on the point order in the array */
    if (pt_indices) pt_order = pt_indices;
    if (!pt_order) {
	tmp_pt_order = (int *)bu_calloc(npts, sizeof(size_t), "temp ordering array");
	for (i = 0; i < npts; i++)
	    tmp_pt_order[i] = (int)i;
	pt_order = (const int *)tmp_pt_order;
    }

    /* Conduct the actual CCW test */
    for (i = 0; i < npts; i++) {
	if (i + 1 == npts) {
	    sum += (pts[pt_order[0]][0] - pts[pt_order[i]][0]) * (pts[pt_order[0]][1] + pts[pt_order[i]][1]);
	} else {
	    sum += (pts[pt_order[i+1]][0] - pts[pt_order[i]][0]) * (pts[pt_order[i+1]][1] + pts[pt_order[i]][1]);
	}
    }

    /* clean up and evaluate results */
    bu_free(tmp_pt_order, "free tmp_pt_order");
    if (NEAR_ZERO(sum, SMALL_FASTF))
	return 0;
    return (sum > 0) ? BG_CW : BG_CCW;
}


int
bg_polygon_copy(struct bg_polygon *destination,
	const struct bg_polygon *source)
{
    if (!destination || !source)
	return 1;
    if (destination == source)
	return 0;

    struct bg_polygon copy = BG_POLYGON_INIT_ZERO;
    if (source->num_contours && !source->contour)
	return 1;

    copy.num_contours = source->num_contours;
    if (!source->num_contours) {
	bg_polygon_clear(destination);
	return 0;
    }

    copy.hole = (int *)bu_calloc(source->num_contours, sizeof(int), "hole");
    copy.contour = (struct bg_poly_contour *)bu_calloc(source->num_contours,
	    sizeof(struct bg_poly_contour), "contour");
    for (size_t i = 0; i < source->num_contours; i++)
	copy.hole[i] = source->hole ? source->hole[i] : 0;
    for (size_t i = 0; i < source->num_contours; i++) {
	if (source->contour[i].num_points && !source->contour[i].point) {
	    bg_polygon_clear(&copy);
	    return 1;
	}
	copy.contour[i].num_points = source->contour[i].num_points;
	copy.contour[i].open = source->contour[i].open;
	if (!source->contour[i].num_points)
	    continue;
	copy.contour[i].point = (point_t *)bu_calloc(
		source->contour[i].num_points, sizeof(point_t), "point");
	for (size_t j = 0; j < source->contour[i].num_points; j++) {
	    VMOVE(copy.contour[i].point[j], source->contour[i].point[j]);
	}
    }

    bg_polygon_clear(destination);
    *destination = copy;
    return 0;
}

int
bg_polygon_move(struct bg_polygon *destination, struct bg_polygon *source)
{
    if (!destination || !source)
	return 1;
    if (destination == source)
	return 0;
    bg_polygon_clear(destination);
    *destination = *source;
    bg_polygon_init(source);
    return 0;
}

size_t
bg_polygon_point_count(const struct bg_polygon *polygon)
{
    size_t count = 0;
    if (!polygon || !polygon->contour)
	return 0;
    for (size_t i = 0; i < polygon->num_contours; i++)
	count += polygon->contour[i].num_points;
    return count;
}

int
bg_polygon_append_point(struct bg_polygon *polygon, size_t contour,
	const point_t point)
{
    if (!polygon || !point || contour >= polygon->num_contours ||
	    !polygon->contour)
	return 1;

    struct bg_poly_contour *c = &polygon->contour[contour];
    const size_t count = c->num_points;
    c->point = (point_t *)bu_realloc(c->point,
	    (count + 1) * sizeof(point_t), "polygon contour points");
    VMOVE(c->point[count], point);
    c->num_points++;
    return 0;
}

int
bg_polygon_remove_point(struct bg_polygon *polygon, size_t contour,
	size_t point_index)
{
    if (!polygon || !polygon->contour || contour >= polygon->num_contours)
	return 1;

    struct bg_poly_contour *c = &polygon->contour[contour];
    if (!c->point || point_index >= c->num_points)
	return 1;

    if (c->num_points > 1) {
	if (point_index + 1 < c->num_points)
	    memmove(&c->point[point_index], &c->point[point_index + 1],
		    (c->num_points - point_index - 1) * sizeof(point_t));
	c->num_points--;
	c->point = (point_t *)bu_realloc(c->point,
	    c->num_points * sizeof(point_t), "polygon contour points");
	return 0;
    }

    bu_free(c->point, "polygon contour points");
    if (contour + 1 < polygon->num_contours) {
	memmove(&polygon->contour[contour], &polygon->contour[contour + 1],
	    (polygon->num_contours - contour - 1) *
	    sizeof(struct bg_poly_contour));
	if (polygon->hole)
	    memmove(&polygon->hole[contour], &polygon->hole[contour + 1],
		    (polygon->num_contours - contour - 1) * sizeof(int));
    }
    polygon->num_contours--;
    if (!polygon->num_contours) {
	bu_free(polygon->contour, "polygon contours");
	polygon->contour = NULL;
	if (polygon->hole) {
	    bu_free(polygon->hole, "polygon holes");
	    polygon->hole = NULL;
	}
	return 0;
    }

    polygon->contour = (struct bg_poly_contour *)bu_realloc(
	polygon->contour,
	polygon->num_contours * sizeof(struct bg_poly_contour),
	"polygon contours");
    if (polygon->hole)
	polygon->hole = (int *)bu_realloc(polygon->hole,
	    polygon->num_contours * sizeof(int), "polygon holes");
    return 0;
}

int
bg_polygon_translate(struct bg_polygon *polygon, const vect_t delta)
{
    if (!polygon || !delta || (polygon->num_contours && !polygon->contour))
	return 1;
    for (size_t i = 0; i < polygon->num_contours; i++) {
	struct bg_poly_contour *c = &polygon->contour[i];
	if (c->num_points && !c->point)
	    return 1;
	for (size_t j = 0; j < c->num_points; j++)
	    VADD2(c->point[j], c->point[j], delta);
    }
    return 0;
}

int
bg_polygon_contour_open_set(struct bg_polygon *polygon, size_t contour,
	int open)
{
    if (!polygon || !polygon->contour || contour >= polygon->num_contours)
	return 1;
    polygon->contour[contour].open = open ? 1 : 0;
    return 0;
}

int
bg_polygon_contours_open_set(struct bg_polygon *polygon, int open)
{
    if (!polygon || (polygon->num_contours && !polygon->contour))
	return 1;
    for (size_t i = 0; i < polygon->num_contours; i++)
	polygon->contour[i].open = open ? 1 : 0;
    return 0;
}

int
bg_polygon_closest_point(const struct bg_polygon *polygon,
	const point_t point, size_t *contour, size_t *point_index,
	fastf_t *distance_sq)
{
    if (!polygon || !point || !polygon->contour)
	return 1;

    fastf_t best = INFINITY;
    size_t best_contour = 0;
    size_t best_point = 0;
    int found = 0;
    for (size_t i = 0; i < polygon->num_contours; i++) {
	const struct bg_poly_contour *c = &polygon->contour[i];
	if (c->num_points && !c->point)
	    return 1;
	for (size_t j = 0; j < c->num_points; j++) {
	    const fastf_t d = DIST_PNT_PNT_SQ(c->point[j], point);
	    if (!found || d < best) {
		best = d;
		best_contour = i;
		best_point = j;
		found = 1;
	    }
	}
    }
    if (!found)
	return 1;
    if (contour)
	*contour = best_contour;
    if (point_index)
	*point_index = best_point;
    if (distance_sq)
	*distance_sq = best;
    return 0;
}

static int
bg_polygon_prepare_single_contour(struct bg_polygon *polygon, size_t count)
{
    if (!polygon || !count)
	return 1;
    struct bg_polygon replacement = BG_POLYGON_INIT_ZERO;
    replacement.num_contours = 1;
    replacement.hole = (int *)bu_calloc(1, sizeof(int), "polygon hole");
    replacement.contour = (struct bg_poly_contour *)bu_calloc(1,
	    sizeof(struct bg_poly_contour), "polygon contour");
    replacement.contour[0].num_points = count;
    replacement.contour[0].point = (point_t *)bu_calloc(count,
	    sizeof(point_t), "polygon points");
    bg_polygon_clear(polygon);
    *polygon = replacement;
    return 0;
}

int
bg_polygon_make_rectangle(struct bg_polygon *polygon, const plane_t plane,
	const point_t first_corner, const point_t opposite_corner, int square)
{
    if (!polygon || !plane || !first_corner || !opposite_corner)
	return 1;

    plane_t local_plane;
    HMOVE(local_plane, plane);
    fastf_t x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    if (bg_plane_closest_pt(&x1, &y1, &local_plane,
	    (point_t *)first_corner) ||
	    bg_plane_closest_pt(&x2, &y2, &local_plane,
	    (point_t *)opposite_corner))
	return 1;
    if (square) {
	const fastf_t dx = x2 - x1;
	const fastf_t dy = y2 - y1;
	const fastf_t side = fmax(fabs(dx), fabs(dy));
	x2 = dx < 0.0 ? x1 - side : x1 + side;
	y2 = dy < 0.0 ? y1 - side : y1 + side;
    }

    struct bg_polygon replacement = BG_POLYGON_INIT_ZERO;
    if (bg_polygon_prepare_single_contour(&replacement, 4))
	return 1;
    replacement.contour[0].open = 0;
    (void)bg_plane_pt_at(&replacement.contour[0].point[0], &local_plane,
	    x1, y1);
    (void)bg_plane_pt_at(&replacement.contour[0].point[1], &local_plane,
	    x1, y2);
    (void)bg_plane_pt_at(&replacement.contour[0].point[2], &local_plane,
	    x2, y2);
    (void)bg_plane_pt_at(&replacement.contour[0].point[3], &local_plane,
	    x2, y1);
    return bg_polygon_move(polygon, &replacement);
}

int
bg_polygon_make_ellipse(struct bg_polygon *polygon, const plane_t plane,
	const point_t center, const point_t radius_point, int circle,
	size_t segment_count)
{
    if (!polygon || !plane || !center || !radius_point || segment_count < 3)
	return 1;

    plane_t local_plane;
    HMOVE(local_plane, plane);
    fastf_t cx = 0.0, cy = 0.0, px = 0.0, py = 0.0;
    if (bg_plane_closest_pt(&cx, &cy, &local_plane, (point_t *)center) ||
	    bg_plane_closest_pt(&px, &py, &local_plane,
	    (point_t *)radius_point))
	return 1;
    fastf_t rx = px - cx;
    fastf_t ry = py - cy;
    if (circle) {
	const fastf_t r = sqrt(rx * rx + ry * ry);
	rx = r;
	ry = r;
    }

    struct bg_polygon replacement = BG_POLYGON_INIT_ZERO;
    if (bg_polygon_prepare_single_contour(&replacement, segment_count))
	return 1;
    replacement.contour[0].open = 0;
    const double two_pi = 6.283185307179586476925286766559;
    for (size_t i = 0; i < segment_count; i++) {
	const double angle = two_pi * (double)i / (double)segment_count;
	(void)bg_plane_pt_at(&replacement.contour[0].point[i], &local_plane,
		cx + cos(angle) * rx, cy + sin(angle) * ry);
    }
    return bg_polygon_move(polygon, &replacement);
}

void
bg_polygon_plot_2d(const char *filename, const point2d_t *pnts, int npnts, int r, int g, int b)
{
    point_t bnp;
    FILE* plot_file = fopen(filename, "wb");
    pl_color(plot_file, r, g, b);

    VSET(bnp, pnts[0][X], pnts[0][Y], 0);
    pdv_3move(plot_file, bnp);

    for (int i = 1; i < npnts; i++) {
	VSET(bnp, pnts[i][X], pnts[i][Y], 0);
	pdv_3cont(plot_file, bnp);
    }

    VSET(bnp, pnts[0][X], pnts[0][Y], 0);
    pdv_3cont(plot_file, bnp);

    fclose(plot_file);
}

void
bg_polygon_plot(const char *filename, const point_t *pnts, int npnts, int r, int g, int b)
{
    point_t bnp;
    FILE* plot_file = fopen(filename, "wb");
    pl_color(plot_file, r, g, b);

    VSET(bnp, pnts[0][X], pnts[0][Y], 0);
    pdv_3move(plot_file, bnp);

    for (int i = 1; i < npnts; i++) {
	VSET(bnp, pnts[i][X], pnts[i][Y], pnts[i][Z]);
	pdv_3cont(plot_file, bnp);
    }

    VSET(bnp, pnts[0][X], pnts[0][Y], pnts[0][Z]);
    pdv_3cont(plot_file, bnp);

    fclose(plot_file);
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
