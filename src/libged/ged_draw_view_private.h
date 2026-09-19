/*              G E D _ D R A W _ V I E W _ P R I V A T E . H
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
/** @file ged_draw_view_private.h
 *
 * Private view adapter surface for libged internals.
 *
 * Public view feature and polygon declarations live in ged/draw.h.  This header
 * only carries private selection routing and transaction hooks needed by GED's
 * view adapter.
 */

#ifndef LIBGED_GED_DRAW_VIEW_PRIVATE_H
#define LIBGED_GED_DRAW_VIEW_PRIVATE_H

#include "common.h"

#include "ged/draw.h"

#ifdef __cplusplus
class BObolViewAttachment;

/* Value-handle owner identities.  The low bit records local/shared scope;
 * the remaining bits identify a live GED view host record. */
extern uint64_t ged_view_context_reference_owner(
	const struct ged_view_context *view_ctx, int local);
extern BObolViewAttachment *ged_view_context_obol_attachment(
	const struct ged_view_context *view_ctx);
extern int ged_view_context_obol_attachment_bind(
	struct ged_view_context *view_ctx,
	BObolViewAttachment *attachment);
extern int ged_view_context_obol_attachment_unbind(
	struct ged_view_context *view_ctx,
	BObolViewAttachment *attachment);
#endif

__BEGIN_DECLS

enum ged_selection_kind {
    GED_SELECTION_ALL = -1,
    GED_SELECTION_SELECTED_PATH = 0,
    GED_SELECTION_HIGHLIGHTED_REF = 1
};

extern int ged_selection_available(struct ged_view_context *view_ctx);
extern int ged_selection_contains_path(
	struct ged_view_context *view_ctx,
	enum ged_selection_kind kind,
	const char *path);
extern int ged_selection_add_path(
	struct ged_view_context *view_ctx,
	enum ged_selection_kind kind,
	const char *path);
extern int ged_selection_remove_path(
	struct ged_view_context *view_ctx,
	enum ged_selection_kind kind,
	const char *path);
extern int ged_selection_apply_path_delta(
	struct ged_view_context *view_ctx,
	enum ged_selection_kind kind,
	const char *const *added_paths,
	size_t added_count,
	const char *const *removed_paths,
	size_t removed_count);
extern int ged_selection_set_path(
	struct ged_view_context *view_ctx,
	enum ged_selection_kind kind,
	const char *path);
extern void ged_draw_source_lod_bounds_callback_set(struct ged_view_context *view_ctx);
extern int ged_draw_source_lod_bounds_callback_is(const struct ged_view_context *view_ctx);

__END_DECLS

#endif /* LIBGED_GED_DRAW_VIEW_PRIVATE_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
