/*                 P O L Y G O N _ T Y P E S. H
 * BRL-CAD
 *
 * Copyright (c) 2004-2026 United States Government as represented by
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

/*----------------------------------------------------------------------*/
/* @file polygon_types.h */
/** @addtogroup bg_polygon */
/** @{ */

/**
 *  @brief Functions for working with polygons
 */

#ifndef BG_POLYGON_TYPES_H
#define BG_POLYGON_TYPES_H

#include <stddef.h>

#include "common.h"
#include "vmath.h"

__BEGIN_DECLS

/** Boolean operation applied to one or more planar polygons. */
enum bg_polygon_boolean_op {
    BG_POLYGON_BOOLEAN_NONE = 0,
    BG_POLYGON_BOOLEAN_UNION,
    BG_POLYGON_BOOLEAN_DIFFERENCE,
    BG_POLYGON_BOOLEAN_INTERSECTION,
    BG_POLYGON_BOOLEAN_XOR
};

struct bg_poly_contour {
    size_t    num_points;
    point_t   *point;               /* in model coordinates */
    int       open;                 /* 0 = closed, 1 = open */
};

struct bg_polygon {
    size_t                  num_contours;
    int                     *hole;
    struct bg_poly_contour  *contour;
};

#define BG_POLYGON_INIT_ZERO {0, NULL, NULL}

struct bg_polygons {
    size_t            num_polygons;
    struct bg_polygon *polygon;
};

#define BG_POLYGONS_INIT_ZERO {0, NULL}

__END_DECLS

#endif  /* BG_POLYGON_TYPES_H */
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
