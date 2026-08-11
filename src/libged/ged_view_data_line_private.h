/*         G E D _ V I E W _ D A T A _ L I N E _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged_view_data_line_private.h
 *
 * In-tree command-state contract shared by libged's legacy data_lines command
 * and libtclcad's point-move adapter.  This is not a retained-render API: the
 * producer owns these values and publishes copies through feature batches.
 */

#ifndef LIBGED_GED_VIEW_DATA_LINE_PRIVATE_H
#define LIBGED_GED_VIEW_DATA_LINE_PRIVATE_H

#include "common.h"

#include "ged/defines.h"
#include "vmath.h"

__BEGIN_DECLS

struct ged_view_data_line_state {
    int draw;
    int color[3];
    int line_width;
    point_t *points;
    size_t point_count;
};

#define GED_VIEW_DATA_LINE_STATE_INIT {0, {255, 255, 0}, 0, NULL, 0}

GED_EXPORT extern void ged_view_data_line_state_clear(
	struct ged_view_data_line_state *state);
GED_EXPORT extern int ged_view_data_line_state_get(
	const struct ged_view_context *view_ctx, int staged,
	struct ged_view_data_line_state *state);
GED_EXPORT extern int ged_view_data_line_state_replace(
	struct ged_view_context *view_ctx, int staged,
	const struct ged_view_data_line_state *state);
GED_EXPORT extern int ged_view_data_line_state_publish(
	struct ged_view_context *view_ctx, int staged);

__END_DECLS

#endif /* LIBGED_GED_VIEW_DATA_LINE_PRIVATE_H */
