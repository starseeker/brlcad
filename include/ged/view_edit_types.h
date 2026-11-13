/*                 V I E W _ E D I T _ T Y P E S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged/view_edit_types.h
 *
 * Renderer-neutral, view-owned edit-preview values.  References contain no
 * pointer and do not retain their owning view.
 */

#ifndef GED_VIEW_EDIT_TYPES_H
#define GED_VIEW_EDIT_TYPES_H

#include <stddef.h>
#include <stdint.h>

#include "vmath.h"
#include "ged/defines.h"

struct bg_tess_tol;
struct bn_tol;
struct db_i;
struct rt_db_internal;

__BEGIN_DECLS

typedef struct ged_view_edit_ref {
    uint64_t owner;
    uint64_t id;
    uint64_t generation;
} ged_view_edit_ref;

#define GED_VIEW_EDIT_REF_NULL_INIT {0, 0, 0}
#ifdef __cplusplus
#  define GED_VIEW_EDIT_REF_NULL ged_view_edit_ref{0, 0, 0}
#else
#  define GED_VIEW_EDIT_REF_NULL ((ged_view_edit_ref){0, 0, 0})
#endif

enum ged_view_edit_geometry_family {
    GED_VIEW_EDIT_GEOMETRY_UNKNOWN = 0,
    GED_VIEW_EDIT_GEOMETRY_TRANSIENT_PREVIEW = 1
};

enum ged_view_edit_preview_event {
    GED_VIEW_EDIT_PREVIEW_BEGIN = 0,
    GED_VIEW_EDIT_PREVIEW_UPDATE,
    GED_VIEW_EDIT_PREVIEW_COMMIT,
    GED_VIEW_EDIT_PREVIEW_CANCEL,
    GED_VIEW_EDIT_PREVIEW_REPLACE_SOURCE,
    GED_VIEW_EDIT_PREVIEW_DISCARD
};

struct ged_view_edit_transaction {
    enum ged_view_edit_preview_event event;
    ged_view_edit_ref feature;
    const char *feature_name;
    const char *source_path;
    const char *edit_intent_id;
    const char *edit_intent_role;
    const point_t *points;
    const int *commands;
    size_t point_count;
    struct db_i *dbip;
    struct rt_db_internal *internal;
    const fastf_t *matrix;
    const struct bg_tess_tol *ttol;
    const struct bn_tol *tol;
    uint32_t source_revision;
    uint32_t inputs_revision;
    int color_valid;
    unsigned char color[3];
};

__END_DECLS

#endif /* GED_VIEW_EDIT_TYPES_H */
