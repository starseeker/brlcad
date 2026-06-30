/*          B S G _ G E D _ D R A W _ V I E W _ P R I V A T E . H
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
/** @file bsg_ged_draw_view_private.h
 *
 * Narrow BSG-backed view adapter surface for libged internals.
 *
 * Public view feature and polygon declarations live in ged/draw.h.  This header
 * only carries private selection routing and transaction hooks still needed by
 * the transitional view adapter.
 */

#ifndef LIBGED_BSG_GED_DRAW_VIEW_PRIVATE_H
#define LIBGED_BSG_GED_DRAW_VIEW_PRIVATE_H

#include "common.h"

#include "ged/draw.h"

__BEGIN_DECLS

enum ged_draw_view_selection_kind {
    GED_DRAW_VIEW_SELECTION_ALL = -1,
    GED_DRAW_VIEW_SELECTION_SELECTED_PATH = 0,
    GED_DRAW_VIEW_SELECTION_HIGHLIGHTED_REF = 1
};

GED_EXPORT extern int ged_draw_view_context_selection_available(void *view_ctx);
GED_EXPORT extern int ged_draw_view_context_selection_contains_path(
	void *view_ctx,
	enum ged_draw_view_selection_kind kind,
	const char *path);
GED_EXPORT extern int ged_draw_view_context_selection_add_path(
	void *view_ctx,
	enum ged_draw_view_selection_kind kind,
	const char *path);
GED_EXPORT extern int ged_draw_view_context_selection_set_path(
	void *view_ctx,
	enum ged_draw_view_selection_kind kind,
	const char *path);
GED_EXPORT extern void ged_draw_view_context_lod_bounds_callback_set(void *view_ctx);

__END_DECLS

#endif /* LIBGED_BSG_GED_DRAW_VIEW_PRIVATE_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
