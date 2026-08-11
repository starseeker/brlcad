/*             V I E W _ P O L Y G O N _ T Y P E S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged/view_polygon_types.h
 *
 * Value types for view-owned polygon editing.
 */

#ifndef GED_VIEW_POLYGON_TYPES_H
#define GED_VIEW_POLYGON_TYPES_H

#include <stddef.h>
#include <stdint.h>

#include "vmath.h"
#include "bg/polygon_types.h"
#include "bu/color.h"
#include "ged/defines.h"

__BEGIN_DECLS

enum ged_view_polygon_type {
    GED_VIEW_POLYGON_GENERAL = 0,
    GED_VIEW_POLYGON_CIRCLE = 1,
    GED_VIEW_POLYGON_ELLIPSE = 2,
    GED_VIEW_POLYGON_RECTANGLE = 3,
    GED_VIEW_POLYGON_SQUARE = 4
};

enum ged_view_polygon_update {
    GED_VIEW_POLYGON_UPDATE_DEFAULT = 0,
    GED_VIEW_POLYGON_UPDATE_PROPS_ONLY = 1,
    GED_VIEW_POLYGON_UPDATE_PT_SELECT = 2,
    GED_VIEW_POLYGON_UPDATE_PT_SELECT_CLEAR = 3,
    GED_VIEW_POLYGON_UPDATE_PT_MOVE = 4,
    GED_VIEW_POLYGON_UPDATE_PT_APPEND = 5
};

/** Opaque owner/id/generation reference to a managed view polygon. */
typedef struct ged_view_polygon_ref {
    uint64_t owner;
    uint64_t id;
    uint64_t generation;
} ged_view_polygon_ref;

#define GED_VIEW_POLYGON_REF_NULL_INIT {0, 0, 0}

#ifdef __cplusplus
#  define GED_VIEW_POLYGON_REF_NULL ged_view_polygon_ref{0, 0, 0}
#else
#  define GED_VIEW_POLYGON_REF_NULL ((ged_view_polygon_ref){0, 0, 0})
#endif

struct ged_view_polygon_record {
    ged_view_polygon_ref ref;
    const char *name;
    enum ged_view_polygon_type type;
    int fill_flag;
    vect2d_t fill_dir;
    fastf_t fill_delta;
    struct bu_color fill_color;
    unsigned char edge_color[3];
    long curr_contour_i;
    long curr_point_i;
    int first_contour_open;
    size_t contour_count;
    size_t point_count;
    point_t origin_point;
    plane_t vp;
    fastf_t vZ;
    void *user_data;
};

typedef int (*ged_view_polygon_record_cb)(
    ged_view_polygon_ref ref,
    const struct ged_view_polygon_record *record,
    void *data);

__END_DECLS

#endif /* GED_VIEW_POLYGON_TYPES_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
