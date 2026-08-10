/*                        S K E T C H . H
 * BRL-CAD
 *
 * Copyright (c) 1993-2026 United States Government as represented by
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
/** @addtogroup rt_sketch */
/** @{ */
/** @file rt/primitives/sketch.h */

#ifndef RT_PRIMITIVES_SKETCH_H
#define RT_PRIMITIVES_SKETCH_H

#include "common.h"
#include "vmath.h"
#include "bg/polygon_types.h"
#include "bu/color.h"
#include "bu/list.h"
#include "bu/vls.h"
#include "bn/tol.h"
#include "rt/defines.h"
#include "rt/directory.h"
#include "rt/db_instance.h"

__BEGIN_DECLS

struct bg_tess_tol;
struct rt_db_internal;
struct rt_primitive_lod_realization;
struct rt_sketch_polygon;

#define RT_SKETCH_POLYGON_GENERAL   0
#define RT_SKETCH_POLYGON_CIRCLE    1
#define RT_SKETCH_POLYGON_ELLIPSE   2
#define RT_SKETCH_POLYGON_RECTANGLE 3
#define RT_SKETCH_POLYGON_SQUARE    4

/**
 * Public value representation of view-polygon data stored in a SKETCH.
 *
 * The bg_polygon memory is caller-owned.  Initialize with
 * rt_sketch_polygon_data_init, release with rt_sketch_polygon_data_free, and
 * use rt_sketch_polygon_data_copy for deep copies.
 */
struct rt_sketch_polygon_data {
    int type;
    int fill_flag;
    vect2d_t fill_dir;
    fastf_t fill_delta;
    struct bu_color fill_color;
    point_t origin_point;
    plane_t vp;
    fastf_t vZ;
    struct bg_polygon polygon;

    int have_edge_color;
    struct bu_color edge_color;
};

/* SKETCH specific editing info */
struct rt_sketch_edit {
    int curr_vert;  /* index of the currently selected vertex (-1 = none) */
    int curr_seg;   /* index of the currently selected curve segment (-1 = none) */
    /* Mouse-proximity pick: ft_edit_xy stores cursor here; ft_edit reads it */
    point_t v_pos;      /* view-space cursor position set by ft_edit_xy */
    int v_pos_valid;    /* non-zero when v_pos holds a pending proximity query */
};

RT_EXPORT extern int rt_check_curve(const struct rt_curve *crv,
				    const struct rt_sketch_internal *skt,
				    int noisy);

RT_EXPORT extern void rt_curve_reverse_segment(uint32_t *lng);
RT_EXPORT extern void rt_curve_order_segments(struct rt_curve *crv);

RT_EXPORT extern void rt_copy_curve(struct rt_curve *crv_out,
				    const struct rt_curve *crv_in);

RT_EXPORT extern void rt_curve_free(struct rt_curve *crv);
RT_EXPORT extern void rt_copy_curve(struct rt_curve *crv_out,
				    const struct rt_curve *crv_in);
RT_EXPORT extern struct rt_sketch_internal *rt_copy_sketch(const struct rt_sketch_internal *sketch_ip);

RT_EXPORT extern int rt_sketch_wireframe_line_set(struct rt_primitive_lod_realization *realization, struct rt_db_internal *ip, const struct bg_tess_tol *ttol);

RT_EXPORT extern void rt_sketch_polygon_data_init(struct rt_sketch_polygon_data *poly);

RT_EXPORT extern void rt_sketch_polygon_data_free(struct rt_sketch_polygon_data *poly);

RT_EXPORT extern int rt_sketch_polygon_data_copy(struct rt_sketch_polygon_data *dest,
						 const struct rt_sketch_polygon_data *src);

RT_EXPORT extern int
db_sketch_to_polygon_data(struct rt_sketch_polygon_data *poly,
			  const char *sname,
			  struct db_i *dbip,
			  struct directory *dp);

RT_EXPORT extern struct directory *
db_sketch_polygon_data_to_sketch(struct db_i *dbip,
				 const char *sname,
				 const struct rt_sketch_polygon_data *poly);

RT_EXPORT extern struct directory *
db_sketch_polygon_to_sketch(struct db_i *dbip, const char *sname, const struct rt_sketch_polygon *poly, const unsigned char edge_rgb[3]);

RT_EXPORT extern struct rt_sketch_polygon *
db_sketch_to_polygon(const char *sname, struct db_i *dbip, struct directory *dp);

RT_EXPORT extern const struct bg_polygon *
rt_sketch_polygon_bg_polygon(const struct rt_sketch_polygon *poly);

RT_EXPORT extern void
rt_sketch_polygon_destroy(struct rt_sketch_polygon *poly);

__END_DECLS

/** @} */

#endif /* RT_PRIMITIVES_SKETCH_H */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
