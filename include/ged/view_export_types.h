/*                  V I E W _ E X P O R T _ T Y P E S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged/view_export_types.h
 *
 * Renderer-neutral values used to export the current presented view.
 */

#ifndef GED_VIEW_EXPORT_TYPES_H
#define GED_VIEW_EXPORT_TYPES_H

#include <stddef.h>

#include "vmath.h"

struct ged_view_context;

/** Aggregate properties of visible database presentation geometry. */
struct ged_view_database_export_summary {
    size_t object_count;
    int has_unsupported_subtraction;
};

typedef int (*ged_view_export_segment_func_t)(
    const point_t start, const point_t end, void *client_data);

#endif /* GED_VIEW_EXPORT_TYPES_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
