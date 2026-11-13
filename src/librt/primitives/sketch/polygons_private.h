/*              P O L Y G O N S _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
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
/** @file primitives/sketch/polygons_private.h
 *
 * Private neutral sketch polygon state shared by the RT implementation and
 * transitional legacy adapters.
 */

#ifndef LIBRT_PRIMITIVES_SKETCH_POLYGONS_PRIVATE_H
#define LIBRT_PRIMITIVES_SKETCH_POLYGONS_PRIVATE_H

#include "common.h"

#include "rt/primitives/sketch.h"

__BEGIN_DECLS

struct rt_sketch_polygon {
    int                 type;
    int                 fill_flag;
    vect2d_t            fill_dir;
    fastf_t             fill_delta;
    struct bu_color     fill_color;
    long                curr_contour_i;
    long                curr_point_i;
    point_t             origin_point;
    plane_t             vp;
    fastf_t             vZ;
    struct bg_polygon   polygon;
    void                *u_data;

    int                 have_edge_color;
    struct bu_color     edge_color;
};

__END_DECLS

#endif /* LIBRT_PRIMITIVES_SKETCH_POLYGONS_PRIVATE_H */
