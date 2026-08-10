/* S C E N E _ R E C O R D _ T Y P E S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged/scene_record_types.h
 *
 * Lightweight types used by the legacy realized-record query surface.  New
 * semantic scene clients normally need only ged/scene_types.h.
 */

#ifndef GED_SCENE_RECORD_TYPES_H
#define GED_SCENE_RECORD_TYPES_H

#include <stddef.h>
#include <stdint.h>

#include "vmath.h"

struct ged_draw_shape_record;
struct ged_draw_group_record;
struct ged_draw_view_db_object_record;
struct ged_draw_shape_candidate;
struct ged_draw_view_record_query;
struct ged_draw_view_rendered_object_summary_s;

typedef struct ged_draw_shape_ref {
    uintptr_t token;
    uint64_t scene_revision;
} ged_draw_shape_ref;

typedef struct ged_draw_group_ref {
    uintptr_t token;
    uint64_t scene_revision;
} ged_draw_group_ref;

#define GED_DRAW_SHAPE_REF_NULL_INIT {0, 0}
#define GED_DRAW_GROUP_REF_NULL_INIT {0, 0}

#ifdef __cplusplus
#  define GED_DRAW_SHAPE_REF_NULL ged_draw_shape_ref{0, 0}
#  define GED_DRAW_GROUP_REF_NULL ged_draw_group_ref{0, 0}
#else
#  define GED_DRAW_SHAPE_REF_NULL ((ged_draw_shape_ref){0, 0})
#  define GED_DRAW_GROUP_REF_NULL ((ged_draw_group_ref){0, 0})
#endif

typedef int (*ged_draw_view_db_object_record_cb)(
    const struct ged_draw_view_db_object_record *record,
    void *client_data);

typedef int (*ged_draw_view_segment_cb)(const point_t start,
	const point_t end, void *client_data);

typedef int (*ged_draw_view_point_cb)(const point_t point,
	void *client_data);

typedef int (*ged_draw_shape_ref_index_cb)(ged_draw_shape_ref ref,
	void *client_data);

typedef int (*ged_draw_shape_candidate_cb)(
    const struct ged_draw_shape_candidate *candidate,
    void *client_data);

typedef struct ged_draw_view_rendered_object_summary_s
    ged_draw_view_rendered_object_summary_t;

#endif /* GED_SCENE_RECORD_TYPES_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
