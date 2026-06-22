/*              P O L Y G O N S _ L E G A C Y _ B S G . C
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
/** @file polygons_legacy_bsg.c
 *
 * Transitional BSG view-polygon adapters for RT sketch polygon conversion.
 */

#include "common.h"

#include <string.h>

#include "vmath.h"
#include "bu/avs.h"
#include "bu/color.h"
#include "bu/opt.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "bsg/polygon.h"
#include "bsg/view_state.h"
#include "rt/db_attr.h"
#include "rt/db5.h"
#include "rt/primitives/sketch_legacy_bsg.h"

#include "./polygons_private.h"

static int
rt_sketch_polygon_type_from_bsg(int type)
{
    switch (type) {
	case BSG_POLYGON_CIRCLE:
	    return RT_SKETCH_POLYGON_CIRCLE;
	case BSG_POLYGON_ELLIPSE:
	    return RT_SKETCH_POLYGON_ELLIPSE;
	case BSG_POLYGON_RECTANGLE:
	    return RT_SKETCH_POLYGON_RECTANGLE;
	case BSG_POLYGON_SQUARE:
	    return RT_SKETCH_POLYGON_SQUARE;
	default:
	    return RT_SKETCH_POLYGON_GENERAL;
    }
}

static int
rt_sketch_polygon_type_to_bsg(int type)
{
    switch (type) {
	case RT_SKETCH_POLYGON_CIRCLE:
	    return BSG_POLYGON_CIRCLE;
	case RT_SKETCH_POLYGON_ELLIPSE:
	    return BSG_POLYGON_ELLIPSE;
	case RT_SKETCH_POLYGON_RECTANGLE:
	    return BSG_POLYGON_RECTANGLE;
	case RT_SKETCH_POLYGON_SQUARE:
	    return BSG_POLYGON_SQUARE;
	default:
	    return BSG_POLYGON_GENERAL;
    }
}

static bsg_polygon_ref
rt_sketch_polygon_ref_to_bsg(rt_view_polygon_ref ref)
{
    bsg_polygon_ref bsg_ref = { ref.token, ref.revision };
    return bsg_ref;
}

static rt_view_polygon_ref
rt_sketch_polygon_ref_from_bsg(bsg_polygon_ref ref)
{
    rt_view_polygon_ref rt_ref = { ref.token, ref.revision };
    return rt_ref;
}

static void
rt_sketch_polygon_to_bsg(struct bsg_polygon *bp, const struct rt_sketch_polygon *poly)
{
    memset(bp, 0, sizeof(*bp));
    if (!poly)
	return;

    bp->type = rt_sketch_polygon_type_to_bsg(poly->type);
    bp->fill_flag = poly->fill_flag;
    V2MOVE(bp->fill_dir, poly->fill_dir);
    bp->fill_delta = poly->fill_delta;
    BU_COLOR_CPY(&bp->fill_color, &poly->fill_color);
    bp->curr_contour_i = poly->curr_contour_i;
    bp->curr_point_i = poly->curr_point_i;
    VMOVE(bp->origin_point, poly->origin_point);
    HMOVE(bp->vp, poly->vp);
    bp->vZ = poly->vZ;
    bp->polygon = poly->polygon;
    bp->u_data = poly->u_data;
}

static void
rt_sketch_polygon_from_bsg(struct rt_sketch_polygon *poly, const struct bsg_polygon *bp)
{
    memset(poly, 0, sizeof(*poly));
    if (!bp)
	return;

    poly->type = rt_sketch_polygon_type_from_bsg(bp->type);
    poly->fill_flag = bp->fill_flag;
    V2MOVE(poly->fill_dir, bp->fill_dir);
    poly->fill_delta = bp->fill_delta;
    BU_COLOR_CPY(&poly->fill_color, &bp->fill_color);
    poly->curr_contour_i = bp->curr_contour_i;
    poly->curr_point_i = bp->curr_point_i;
    VMOVE(poly->origin_point, bp->origin_point);
    HMOVE(poly->vp, bp->vp);
    poly->vZ = bp->vZ;
    poly->polygon = bp->polygon;
    poly->u_data = bp->u_data;
}

static void
db_sketch_view_attrs_apply(struct db_i *dbip, struct directory *dp, struct bsg_view *sv)
{
    if (!dbip || !dp || !sv)
	return;

    int have_view = 1;
    struct bu_attribute_value_set lavs;
    bu_avs_init_empty(&lavs);
    if (db5_get_attributes(dbip, &lavs, dp)) {
	bu_avs_free(&lavs);
	return;
    }

    const char *val = bu_avs_get(&lavs, "VIEWSCALE");
    if (val) {
	bu_opt_fastf_t(NULL, 1, (const char **)&val, (void *)&sv->gv_scale);
    } else {
	have_view = 0;
    }

    if (have_view) {
	val = bu_avs_get(&lavs, "ROTATION");
	if (val) {
	    quat_t quat;
	    char *av[5] = {NULL};
	    char *lp = bu_strdup(val);
	    if (bu_argv_from_string(av, 4, lp) != 4) {
		have_view = 0;
	    } else {
		bu_opt_fastf_t(NULL, 1, (const char **)&av[0], (void *)&quat[0]);
		bu_opt_fastf_t(NULL, 1, (const char **)&av[1], (void *)&quat[1]);
		bu_opt_fastf_t(NULL, 1, (const char **)&av[2], (void *)&quat[2]);
		bu_opt_fastf_t(NULL, 1, (const char **)&av[3], (void *)&quat[3]);
		quat_quat2mat(sv->gv_rotation, quat);
	    }
	    bu_free(lp, "val cpy");
	} else {
	    have_view = 0;
	}
    }

    if (have_view) {
	val = bu_avs_get(&lavs, "CENTER");
	if (val) {
	    quat_t quat;
	    char *av[5] = {NULL};
	    char *lp = bu_strdup(val);
	    if (bu_argv_from_string(av, 4, lp) != 4) {
		have_view = 0;
	    } else {
		bu_opt_fastf_t(NULL, 1, (const char **)&av[0], (void *)&quat[0]);
		bu_opt_fastf_t(NULL, 1, (const char **)&av[1], (void *)&quat[1]);
		bu_opt_fastf_t(NULL, 1, (const char **)&av[2], (void *)&quat[2]);
		bu_opt_fastf_t(NULL, 1, (const char **)&av[3], (void *)&quat[3]);
		quat_quat2mat(sv->gv_center, quat);
	    }
	    bu_free(lp, "val cpy");
	}
    }

    bu_avs_free(&lavs);
}

rt_view_polygon_ref
db_sketch_to_view_polygon_ref(const char *sname, struct db_i *dbip, struct directory *dp, void *view_ctx)
{
    return db_sketch_to_view_polygon_scoped_ref(sname, dbip, dp, view_ctx, 0);
}

rt_view_polygon_ref
db_sketch_to_view_polygon_scoped_ref(const char *sname, struct db_i *dbip, struct directory *dp, void *view_ctx, int local)
{
    struct bsg_view *sv = (struct bsg_view *)view_ctx;

    if (!sv)
	return RT_VIEW_POLYGON_REF_NULL;

    struct rt_sketch_polygon *poly = db_sketch_to_polygon(sname, dbip, dp);
    if (!poly)
	return RT_VIEW_POLYGON_REF_NULL;

    db_sketch_view_attrs_apply(dbip, dp, sv);

    struct bsg_polygon bp;
    rt_sketch_polygon_to_bsg(&bp, poly);
    bsg_polygon_ref ref = bsg_polygon_ref_create_from_data(sv, sname, local, &bp);
    if (!bsg_polygon_ref_is_null(ref) && poly->have_edge_color) {
	(void)bsg_polygon_set_visual(ref, &poly->edge_color, NULL,
		poly->fill_dir[0], poly->fill_dir[1], poly->fill_delta, poly->vZ,
		poly->fill_flag);
    }

    rt_sketch_polygon_destroy(poly);
    return rt_sketch_polygon_ref_from_bsg(ref);
}

struct directory *
db_view_polygon_ref_to_sketch(struct db_i *dbip, const char *sname, rt_view_polygon_ref ref)
{
    struct bsg_polygon_record rec;
    bsg_polygon_ref bsg_ref = rt_sketch_polygon_ref_to_bsg(ref);

    if (!bsg_polygon_record_get(bsg_ref, &rec))
	return NULL;

    struct rt_sketch_polygon poly;
    rt_sketch_polygon_from_bsg(&poly, bsg_polygon_data(bsg_ref));
    return db_sketch_polygon_to_sketch(dbip, sname, &poly, rec.edge_color);
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
