/*                S E L E C T I O N _ T Y P E S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged/selection_types.h
 *
 * Renderer-neutral semantic selection value and callback types.
 */

#ifndef GED_SELECTION_TYPES_H
#define GED_SELECTION_TYPES_H

#include <stddef.h>

#include "vmath.h"
#include "ged/defines.h"

struct bu_vls;
struct ged_pick_result;
struct ged_view_context;

__BEGIN_DECLS

enum ged_selection_snap_kind {
    GED_SELECTION_SNAP_GRID = 1,
    GED_SELECTION_SNAP_ENDPOINT = 2
};

/** Stable hash collections derived from one semantic selection set. */
enum ged_selection_hash_kind {
    GED_SELECTION_HASH_SELECTED_PATH = 1,
    GED_SELECTION_HASH_ACTIVE_PATH,
    GED_SELECTION_HASH_ACTIVE_PARENT_PATH,
    GED_SELECTION_HASH_IMMEDIATE_PARENT_OBJECT,
    GED_SELECTION_HASH_GRAND_PARENT_OBJECT
};

typedef int (*ged_view_selection_visit_cb)(
    struct ged_view_context *view_ctx,
    const char *path,
    void *data);

typedef void (*ged_view_selection_path_cb)(const char *path, void *data);

__END_DECLS

#endif /* GED_SELECTION_TYPES_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
