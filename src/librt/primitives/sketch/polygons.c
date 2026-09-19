/*                      P O L Y G O N S . C
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
/** @file polygons.c
 *
 * Sketch routines for view polygons.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "bu/avs.h"
#include "bu/cmd.h"
#include "bu/color.h"
#include "bu/opt.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "bg/plane.h"
#include "bg/polygon.h"
#include "rt/defines.h"
#include "rt/directory.h"
#include "rt/db_attr.h"
#include "rt/db_internal.h"
#include "rt/db5.h"
#include "rt/db_instance.h"
#include "rt/db_io.h"
#include "rt/functab.h"
#include "rt/geom.h"
#include "rt/primitives/sketch.h"

#include "./polygons_private.h"

struct segment_node {
    struct bu_list l;
    int reverse;
    int used;
    void *segment;
};

struct contour_node {
    struct bu_list l;
    struct bu_list head;
};

static void
rt_sketch_polygon_set_default_color(struct bu_color *color, double r, double g, double b)
{
    if (!color)
	return;

    color->buc_rgb[RED] = r;
    color->buc_rgb[GRN] = g;
    color->buc_rgb[BLU] = b;
    color->buc_rgb[ALP] = 0.0;
}

void
rt_sketch_polygon_data_init(struct rt_sketch_polygon_data *poly)
{
    if (!poly)
	return;

    memset(poly, 0, sizeof(*poly));
    poly->type = RT_SKETCH_POLYGON_GENERAL;
    V2SET(poly->fill_dir, 1.0, 0.0);
    poly->fill_delta = 1.0;
    rt_sketch_polygon_set_default_color(&poly->fill_color, 1.0, 1.0, 1.0);
    rt_sketch_polygon_set_default_color(&poly->edge_color, 1.0, 1.0, 1.0);
    HSET(poly->vp, 0.0, 0.0, 1.0, 0.0);
    poly->polygon = (struct bg_polygon)BG_POLYGON_INIT_ZERO;
}

void
rt_sketch_polygon_data_free(struct rt_sketch_polygon_data *poly)
{
    if (!poly)
	return;

    bg_polygon_clear(&poly->polygon);
    rt_sketch_polygon_data_init(poly);
}

int
rt_sketch_polygon_data_copy(struct rt_sketch_polygon_data *dest,
	const struct rt_sketch_polygon_data *src)
{
    struct rt_sketch_polygon_data copy;

    if (!dest || !src)
	return -1;

    if (dest == src)
	return 0;

    copy = *src;
    copy.polygon = (struct bg_polygon)BG_POLYGON_INIT_ZERO;
    if (bg_polygon_copy(&copy.polygon, &src->polygon))
	return -1;

    bg_polygon_clear(&dest->polygon);
    *dest = copy;
    return 0;
}

static int
rt_sketch_polygon_to_data(struct rt_sketch_polygon_data *data,
	const struct rt_sketch_polygon *poly)
{
    if (!data || !poly)
	return -1;

    rt_sketch_polygon_data_init(data);
    data->type = poly->type;
    data->fill_flag = poly->fill_flag;
    V2MOVE(data->fill_dir, poly->fill_dir);
    data->fill_delta = poly->fill_delta;
    BU_COLOR_CPY(&data->fill_color, &poly->fill_color);
    VMOVE(data->origin_point, poly->origin_point);
    HMOVE(data->vp, poly->vp);
    data->vZ = poly->vZ;
    data->have_edge_color = poly->have_edge_color;
    BU_COLOR_CPY(&data->edge_color, &poly->edge_color);
    return bg_polygon_copy(&data->polygon, &poly->polygon);
}

static void
rt_sketch_polygon_from_data(struct rt_sketch_polygon *poly,
	const struct rt_sketch_polygon_data *data)
{
    if (!poly || !data)
	return;

    memset(poly, 0, sizeof(*poly));
    poly->type = data->type;
    poly->fill_flag = data->fill_flag;
    V2MOVE(poly->fill_dir, data->fill_dir);
    poly->fill_delta = data->fill_delta;
    BU_COLOR_CPY(&poly->fill_color, &data->fill_color);
    poly->curr_contour_i = -1;
    poly->curr_point_i = -1;
    VMOVE(poly->origin_point, data->origin_point);
    HMOVE(poly->vp, data->vp);
    poly->vZ = data->vZ;
    poly->polygon = data->polygon;
    poly->have_edge_color = data->have_edge_color;
    BU_COLOR_CPY(&poly->edge_color, &data->edge_color);
}

static struct rt_sketch_polygon *
db_sketch_to_polygon_internal(const char *UNUSED(sname), struct db_i *dbip, struct directory *dp)
{
    if (!dbip || !dp)
	return NULL;

    // Begin import
    size_t ncontours = 0;
    struct bu_list HeadSegmentNodes;
    struct bu_list HeadContourNodes;
    struct segment_node *all_segment_nodes;
    struct segment_node *curr_snode;
    struct contour_node *curr_cnode;

    struct rt_db_internal intern;
    struct rt_sketch_internal *sketch_ip;
    mat_t mat;
    MAT_IDN(mat);
    if (rt_db_get_internal(&intern, dp, dbip, mat) < 0) {
	return NULL;
    }
    sketch_ip = (struct rt_sketch_internal *)intern.idb_ptr;
    RT_SKETCH_CK_MAGIC(sketch_ip);

    if (sketch_ip->vert_count < 3 || sketch_ip->curve.count < 1) {
	rt_db_free_internal(&intern);
	return NULL;
    }

    // Have a sketch - create an empty polygon
    struct rt_sketch_polygon *poly;
    BU_GET(poly, struct rt_sketch_polygon);
    memset(poly, 0, sizeof(*poly));
    poly->type = RT_SKETCH_POLYGON_GENERAL;
    V2SET(poly->fill_dir, 1.0, 0.0);
    poly->fill_delta = 1.0;
    rt_sketch_polygon_set_default_color(&poly->fill_color, 1.0, 1.0, 1.0);
    rt_sketch_polygon_set_default_color(&poly->edge_color, 1.0, 1.0, 1.0);
    poly->curr_contour_i = -1;
    poly->curr_point_i = -1;
    VMOVE(poly->origin_point, sketch_ip->V);

    /* Start translating the sketch info into a polygon */
    all_segment_nodes = (struct segment_node *)bu_calloc(sketch_ip->curve.count, sizeof(struct segment_node), "all_segment_nodes");

    BU_LIST_INIT(&HeadSegmentNodes);
    BU_LIST_INIT(&HeadContourNodes);
    for (size_t n = 0; n < sketch_ip->curve.count; ++n) {
	all_segment_nodes[n].segment = sketch_ip->curve.segment[n];
	all_segment_nodes[n].reverse = sketch_ip->curve.reverse[n];
	BU_LIST_INSERT(&HeadSegmentNodes, &all_segment_nodes[n].l);
    }
    curr_cnode = (struct contour_node *)0;
    while (BU_LIST_NON_EMPTY(&HeadSegmentNodes)) {
	struct segment_node *unused_snode = BU_LIST_FIRST(segment_node, &HeadSegmentNodes);
	uint32_t *magic = (uint32_t *)unused_snode->segment;
	struct line_seg *unused_lsg;

	BU_LIST_DEQUEUE(&unused_snode->l);

	/* For the moment, skipping everything except line segments */
	if (*magic != CURVE_LSEG_MAGIC)
	    continue;

	unused_lsg = (struct line_seg *)unused_snode->segment;
	if (unused_snode->reverse) {
	    int tmp = unused_lsg->start;
	    unused_lsg->start = unused_lsg->end;
	    unused_lsg->end = tmp;
	}

	/* Find a contour to add the unused segment to. */
	for (BU_LIST_FOR(curr_cnode, contour_node, &HeadContourNodes)) {
	    for (BU_LIST_FOR(curr_snode, segment_node, &curr_cnode->head)) {
		struct line_seg *curr_lsg = (struct line_seg *)curr_snode->segment;

		if (unused_lsg->start == curr_lsg->end) {
		    unused_snode->used = 1;
		    BU_LIST_APPEND(&curr_snode->l, &unused_snode->l);
		    goto end;
		}

		if (unused_lsg->end == curr_lsg->start) {
		    unused_snode->used = 1;
		    BU_LIST_INSERT(&curr_snode->l, &unused_snode->l);
		    goto end;
		}
	    }
	}

end:
	if (!unused_snode->used) {
	    ++ncontours;
	    BU_ALLOC(curr_cnode, struct contour_node);
	    BU_LIST_INSERT(&HeadContourNodes, &curr_cnode->l);
	    BU_LIST_INIT(&curr_cnode->head);
	    BU_LIST_INSERT(&curr_cnode->head, &unused_snode->l);
	}
    }

    poly->polygon.num_contours = ncontours;
    poly->polygon.hole = (int *)bu_calloc(ncontours, sizeof(int), "gp_hole");
    poly->polygon.contour = (struct bg_poly_contour *)bu_calloc(ncontours, sizeof(struct bg_poly_contour), "gp_contour");

    size_t j = 0;
    fastf_t dmax = 0.0;
    while (BU_LIST_NON_EMPTY(&HeadContourNodes)) {
	size_t k = 0;
	size_t nsegments = 0;
	struct line_seg *curr_lsg = NULL;

	curr_cnode = BU_LIST_FIRST(contour_node, &HeadContourNodes);
	BU_LIST_DEQUEUE(&curr_cnode->l);

	/* Count the number of segments in this contour */
	for (BU_LIST_FOR(curr_snode, segment_node, &curr_cnode->head))
	    ++nsegments;

	struct segment_node *first_snode =
	    BU_LIST_FIRST(segment_node, &curr_cnode->head);
	struct segment_node *last_snode =
	    BU_LIST_LAST(segment_node, &curr_cnode->head);
	struct line_seg *first_lsg = (struct line_seg *)first_snode->segment;
	struct line_seg *last_lsg = (struct line_seg *)last_snode->segment;
	const int open = first_lsg->start != last_lsg->end;
	const size_t npoints = nsegments + (open ? 1 : 0);

	poly->polygon.contour[j].num_points = npoints;
	poly->polygon.contour[j].open = open;
	poly->polygon.contour[j].point = (point_t *)bu_calloc(npoints, sizeof(point_t), "gpc_point");

	while (BU_LIST_NON_EMPTY(&curr_cnode->head)) {
	    curr_snode = BU_LIST_FIRST(segment_node, &curr_cnode->head);
	    BU_LIST_DEQUEUE(&curr_snode->l);

	    curr_lsg = (struct line_seg *)curr_snode->segment;

	    /* Convert from UV space to model space */
	    VJOIN2(poly->polygon.contour[j].point[k], sketch_ip->V,
		    sketch_ip->verts[curr_lsg->start][0], sketch_ip->u_vec,
		    sketch_ip->verts[curr_lsg->start][1], sketch_ip->v_vec);
	    fastf_t dtmp = DIST_PNT_PNT(sketch_ip->V, poly->polygon.contour[j].point[k]);
	    if (dtmp > dmax)
		dmax = dtmp;
	    ++k;
	}
	if (open) {
	    VJOIN2(poly->polygon.contour[j].point[k], sketch_ip->V,
		    sketch_ip->verts[last_lsg->end][0], sketch_ip->u_vec,
		    sketch_ip->verts[last_lsg->end][1], sketch_ip->v_vec);
	    fastf_t dtmp = DIST_PNT_PNT(sketch_ip->V,
		    poly->polygon.contour[j].point[k]);
	    if (dtmp > dmax)
		dmax = dtmp;
	}

	/* free contour node */
	bu_free((void *)curr_cnode, "curr_cnode");

	++j;
    }

    /* Clean up */
    bu_free((void *)all_segment_nodes, "all_segment_nodes");

    /* Unlike interactive sketch creation, the plane of an imported sketch comes from the
     * sketch parameters. */
    vect_t pn;
    VCROSS(pn, sketch_ip->u_vec, sketch_ip->v_vec);
    bg_plane_pt_nrml(&poly->vp, sketch_ip->V, pn);

    // check attributes for visual properties
    struct bu_attribute_value_set lavs;
    bu_avs_init_empty(&lavs);
    if (!db5_get_attributes(dbip, &lavs, dp)) {
	const char *val = NULL;
	// Check for various polygon properties
	val = bu_avs_get(&lavs, "POLYGON_EDGE_COLOR");
	if (val) {
	    struct bu_color bc;
	    if (bu_opt_color(NULL, 1, (const char **)&val, (void *)&bc) == 1) {
		BU_COLOR_CPY(&poly->edge_color, &bc);
		poly->have_edge_color = 1;
	    }
	}
	val = bu_avs_get(&lavs, "POLYGON_FILL_COLOR");
	if (val) {
	    bu_opt_color(NULL, 1, (const char **)&val, (void *)&poly->fill_color);
	}
	val = bu_avs_get(&lavs, "POLYGON_FILL");
	if (val && BU_STR_EQUAL(val, "1")) {
	    poly->fill_flag = 1;
	}
	val = bu_avs_get(&lavs, "POLYGON_FILL_SLOPE_X");
	if (val) {
	    bu_opt_fastf_t(NULL, 1, (const char **)&val, (void *)&poly->fill_dir[0]);
	}
	val = bu_avs_get(&lavs, "POLYGON_FILL_SLOPE_Y");
	if (val) {
	    bu_opt_fastf_t(NULL, 1, (const char **)&val, (void *)&poly->fill_dir[1]);
	}
	val = bu_avs_get(&lavs, "POLYGON_FILL_DELTA");
	if (val) {
	    bu_opt_fastf_t(NULL, 1, (const char **)&val, (void *)&poly->fill_delta);
	}
	val = bu_avs_get(&lavs, "POLYGON_TYPE");
	if (val && BU_STR_EQUAL(val, "CIRCLE")) {
	    poly->type = RT_SKETCH_POLYGON_CIRCLE;
	}
	if (val && BU_STR_EQUAL(val, "ELLIPSE")) {
	    poly->type = RT_SKETCH_POLYGON_ELLIPSE;
	}
	if (val && BU_STR_EQUAL(val, "RECTANGLE")) {
	    poly->type = RT_SKETCH_POLYGON_RECTANGLE;
	}
	if (val && BU_STR_EQUAL(val, "SQUARE")) {
	    poly->type = RT_SKETCH_POLYGON_SQUARE;
	}
	if (val && BU_STR_EQUAL(val, "GENERAL")) {
	    poly->type = RT_SKETCH_POLYGON_GENERAL;
	}
	val = bu_avs_get(&lavs, "POLYGON_CONTOUR_HOLES");
	if (val) {
	    const char *cp = val;
	    for (size_t i = 0; i < poly->polygon.num_contours; i++) {
		while (*cp && (isspace((unsigned char)*cp) || *cp == ','))
		    cp++;
		if (!*cp)
		    break;
		char *endp = NULL;
		long hval = strtol(cp, &endp, 10);
		if (endp == cp)
		    break;
		poly->polygon.hole[i] = hval ? 1 : 0;
		cp = endp;
	    }
	}
    }
    bu_avs_free(&lavs);

    rt_db_free_internal(&intern);
    return poly;
}

struct rt_sketch_polygon *
db_sketch_to_polygon(const char *sname, struct db_i *dbip, struct directory *dp)
{
    return db_sketch_to_polygon_internal(sname, dbip, dp);
}

int
db_sketch_to_polygon_data(struct rt_sketch_polygon_data *data,
	const char *sname,
	struct db_i *dbip,
	struct directory *dp)
{
    if (!data)
	return -1;

    struct rt_sketch_polygon *poly = db_sketch_to_polygon_internal(sname,
	    dbip, dp);
    if (!poly)
	return -1;

    rt_sketch_polygon_data_free(data);
    if (rt_sketch_polygon_to_data(data, poly)) {
	rt_sketch_polygon_destroy(poly);
	return -1;
    }
    rt_sketch_polygon_destroy(poly);
    return 0;
}

const struct bg_polygon *
rt_sketch_polygon_bg_polygon(const struct rt_sketch_polygon *poly)
{
    return poly ? &poly->polygon : NULL;
}

void
rt_sketch_polygon_destroy(struct rt_sketch_polygon *poly)
{
    if (!poly)
	return;
    bg_polygon_clear(&poly->polygon);
    BU_PUT(poly, struct rt_sketch_polygon);
}

static struct directory *
db_polygon_data_to_sketch(struct db_i *dbip, const char *sname,
	const struct rt_sketch_polygon *poly, const unsigned char edge_rgb[3],
	int update)
{
    if (!dbip || !sname || !sname[0] || !poly)
	return NULL;

    struct directory *dp = db_lookup(dbip, sname, LOOKUP_QUIET);
    if (!update && dp != RT_DIR_NULL) {
	bu_log("Object %s already exists\n", sname);
	return NULL;
    }
    if (update && dp == RT_DIR_NULL) {
	bu_log("Sketch %s does not exist\n", sname);
	return NULL;
    }
    if (update && (dp->d_flags & RT_DIR_COMB || dp->d_minor_type != ID_SKETCH)) {
	bu_log("Object %s is not a sketch\n", sname);
	return NULL;
    }

    size_t num_verts = 0;
    size_t num_segments = 0;
    struct rt_db_internal internal;
    struct rt_sketch_internal *sketch_ip;
    struct line_seg *lsg;
    plane_t vp;
    HMOVE(vp, poly->vp);

    for (size_t j = 0; j < poly->polygon.num_contours; ++j) {
	const size_t npoints = poly->polygon.contour[j].num_points;
	num_verts += npoints;
	if (npoints > 1)
	    num_segments += npoints -
		(poly->polygon.contour[j].open ? 1 : 0);
    }

    if (num_verts < 3 || !num_segments) {
	return NULL;
    }

    RT_DB_INTERNAL_INIT(&internal);
    internal.idb_major_type = DB5_MAJORTYPE_BRLCAD;
    internal.idb_type = ID_SKETCH;
    internal.idb_meth = &OBJ[ID_SKETCH];

    BU_ALLOC(internal.idb_ptr, struct rt_sketch_internal);
    sketch_ip = (struct rt_sketch_internal *)internal.idb_ptr;
    sketch_ip->magic = RT_SKETCH_INTERNAL_MAGIC;
    sketch_ip->vert_count = num_verts;
    sketch_ip->verts = (point2d_t *)bu_calloc(sketch_ip->vert_count, sizeof(point2d_t), "sketch_ip->verts");
    sketch_ip->curve.count = num_segments;
    sketch_ip->curve.reverse = (int *)bu_calloc(sketch_ip->curve.count, sizeof(int), "sketch_ip->curve.reverse");
    sketch_ip->curve.segment = (void **)bu_calloc(sketch_ip->curve.count, sizeof(void *), "sketch_ip->curve.segment");


    /* Plane origin is sketch origin */
    bg_plane_pt_at(&sketch_ip->V, &vp, 0, 0);
    point_t u_end, v_end;
    bg_plane_pt_at(&u_end, &vp, 1, 0);
    bg_plane_pt_at(&v_end, &vp, 0, 1);
    VSUB2(sketch_ip->u_vec, u_end, sketch_ip->V);
    VSUB2(sketch_ip->v_vec, v_end, sketch_ip->V);

    int n = 0;
    int s = 0;
    for (size_t j = 0; j < poly->polygon.num_contours; ++j) {
	size_t cstart = n;
	size_t k = 0;
	for (k = 0; k < poly->polygon.contour[j].num_points; ++k) {
	    bg_plane_closest_pt(&sketch_ip->verts[n][0], &sketch_ip->verts[n][1], &vp, &poly->polygon.contour[j].point[k]);

	    if (k) {
		BU_ALLOC(lsg, struct line_seg);
		sketch_ip->curve.segment[s++] = (void *)lsg;
		lsg->magic = CURVE_LSEG_MAGIC;
		lsg->start = n-1;
		lsg->end = n;
	    }

	    ++n;
	}

	if (k > 1 && !poly->polygon.contour[j].open) {
	    BU_ALLOC(lsg, struct line_seg);
	    sketch_ip->curve.segment[s++] = (void *)lsg;
	    lsg->magic = CURVE_LSEG_MAGIC;
	    lsg->start = n-1;
	    lsg->end = cstart;
	}
    }


    int created = 0;
    if (dp == RT_DIR_NULL) {
	dp = db_diradd(dbip, sname, RT_DIR_PHONY_ADDR, 0, RT_DIR_SOLID,
		(void *)&internal.idb_type);
	if (dp == RT_DIR_NULL) {
	    rt_db_free_internal(&internal);
	    return NULL;
	}
	created = 1;
    }

    if (rt_db_put_internal(dp, dbip, &internal) < 0) {
	if (created)
	    db_dirdelete(dbip, dp);
	return NULL;
    }

    // write attributes to save visual properties

    struct bu_attribute_value_set lavs;
    bu_avs_init_empty(&lavs);
    if (!db5_get_attributes(dbip, &lavs, dp)) {
	struct bu_vls val = BU_VLS_INIT_ZERO;
	unsigned char prgb[3] = {255, 255, 0};
	if (edge_rgb) {
	    prgb[0] = edge_rgb[0];
	    prgb[1] = edge_rgb[1];
	    prgb[2] = edge_rgb[2];
	} else if (poly->have_edge_color) {
	    bu_color_to_rgb_chars(&poly->edge_color, prgb);
	}
	bu_vls_sprintf(&val, "%d/%d/%d", prgb[0], prgb[1], prgb[2]);
	bu_avs_add(&lavs, "POLYGON_EDGE_COLOR", bu_vls_cstr(&val));
	unsigned char rgb[3];
	bu_color_to_rgb_chars(&poly->fill_color, rgb);
	bu_vls_sprintf(&val, "%d/%d/%d", rgb[0], rgb[1], rgb[2]);
	bu_avs_add(&lavs, "POLYGON_FILL_COLOR", bu_vls_cstr(&val));
	bu_vls_sprintf(&val, "%d", poly->fill_flag);
	bu_avs_add(&lavs, "POLYGON_FILL", bu_vls_cstr(&val));
	bu_vls_sprintf(&val, "%g", poly->fill_dir[0]);
	bu_avs_add(&lavs, "POLYGON_FILL_SLOPE_X", bu_vls_cstr(&val));
	bu_vls_sprintf(&val, "%g", poly->fill_dir[1]);
	bu_avs_add(&lavs, "POLYGON_FILL_SLOPE_Y", bu_vls_cstr(&val));
	bu_vls_sprintf(&val, "%g", poly->fill_delta);
	bu_avs_add(&lavs, "POLYGON_FILL_DELTA", bu_vls_cstr(&val));
	switch (poly->type) {
	    case RT_SKETCH_POLYGON_CIRCLE:
		bu_vls_sprintf(&val, "CIRCLE");
		break;
	    case RT_SKETCH_POLYGON_ELLIPSE:
		bu_vls_sprintf(&val, "ELLIPSE");
		break;
	    case RT_SKETCH_POLYGON_RECTANGLE:
		bu_vls_sprintf(&val, "RECTANGLE");
		break;
	    case RT_SKETCH_POLYGON_SQUARE:
		bu_vls_sprintf(&val, "SQUARE");
		break;
	    default:
		bu_vls_sprintf(&val, "GENERAL");
		break;
	}
	bu_avs_add(&lavs, "POLYGON_TYPE", bu_vls_cstr(&val));
	bu_vls_trunc(&val, 0);
	for (size_t j = 0; j < poly->polygon.num_contours; j++) {
	    bu_vls_printf(&val, "%s%d", j ? "," : "",
		    (poly->polygon.hole && poly->polygon.hole[j]) ? 1 : 0);
	}
	bu_avs_add(&lavs, "POLYGON_CONTOUR_HOLES", bu_vls_cstr(&val));
    }
    db5_update_attributes(dp, &lavs, dbip);
    bu_avs_free(&lavs);

    return dp;
}

struct directory *
db_sketch_polygon_to_sketch(struct db_i *dbip, const char *sname, const struct rt_sketch_polygon *poly, const unsigned char edge_rgb[3])
{
    return db_polygon_data_to_sketch(dbip, sname, poly, edge_rgb, 0);
}

struct directory *
db_sketch_polygon_data_to_sketch(struct db_i *dbip,
	const char *sname,
	const struct rt_sketch_polygon_data *data)
{
    if (!data)
	return NULL;

    struct rt_sketch_polygon poly;
    rt_sketch_polygon_from_data(&poly, data);
    return db_polygon_data_to_sketch(dbip, sname, &poly, NULL, 0);
}

struct directory *
db_sketch_polygon_data_update_sketch(struct db_i *dbip,
	const char *sname,
	const struct rt_sketch_polygon_data *data)
{
    if (!data)
	return NULL;

    struct rt_sketch_polygon poly;
    rt_sketch_polygon_from_data(&poly, data);
    return db_polygon_data_to_sketch(dbip, sname, &poly, NULL, 1);
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
