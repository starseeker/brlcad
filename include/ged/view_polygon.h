/* V I E W _ P O L Y G O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged/view_polygon.h
 *
 * View-owned polygon editing.  Every operation that resolves a reference
 * takes its owning view explicitly.  References from another view, or stale
 * references from an earlier store generation, fail without changing state.
 */

#ifndef GED_VIEW_POLYGON_H
#define GED_VIEW_POLYGON_H

#include "ged/view_polygon_types.h"

__BEGIN_DECLS

GED_EXPORT extern int
ged_view_polygon_ref_is_null(ged_view_polygon_ref ref);

GED_EXPORT extern ged_view_polygon_ref
ged_view_polygon_find(struct ged_view_context *view_ctx, const char *name);

GED_EXPORT extern ged_view_polygon_ref
ged_view_polygon_find_scoped(struct ged_view_context *view_ctx,
	const char *name, int local_only);

GED_EXPORT extern ged_view_polygon_ref
ged_view_polygon_select(struct ged_view_context *view_ctx,
	const point_t model_point);

GED_EXPORT extern ged_view_polygon_ref
ged_view_polygon_create(struct ged_view_context *view_ctx, const char *name,
	int local, enum ged_view_polygon_type type, const point_t model_point);

GED_EXPORT extern ged_view_polygon_ref
ged_view_polygon_dup(struct ged_view_context *view_ctx, const char *name,
	const char *new_name);

GED_EXPORT extern ged_view_polygon_ref
ged_view_polygon_import_sketch(struct ged_view_context *view_ctx,
	const char *name, struct db_i *dbip, struct directory *dp, int local);

GED_EXPORT extern void
ged_view_polygon_visit_records(struct ged_view_context *view_ctx,
	ged_view_polygon_record_cb callback, void *data);

GED_EXPORT extern size_t
ged_view_polygon_snap_count(struct ged_view_context *view_ctx,
	ged_view_polygon_ref exclude);

GED_EXPORT extern int
ged_view_polygon_clear_point_selection(struct ged_view_context *view_ctx);

GED_EXPORT extern int
ged_view_polygon_snap_exclude_set(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);

GED_EXPORT extern struct directory *
ged_view_polygon_export_sketch(struct ged_view_context *view_ctx,
	struct db_i *dbip, const char *name, ged_view_polygon_ref ref);

GED_EXPORT extern int
ged_view_polygon_record_get(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, struct ged_view_polygon_record *record);

GED_EXPORT extern int
ged_view_polygon_has_data(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);

GED_EXPORT extern int
ged_view_polygon_update(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, enum ged_view_polygon_update op);

GED_EXPORT extern int
ged_view_polygon_update_screen_pt(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, int x, int y, enum ged_view_polygon_update op);

GED_EXPORT extern int
ged_view_polygon_move(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const point_t current_point,
	const point_t previous_point);

GED_EXPORT extern int
ged_view_polygon_set_name(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const char *name);

GED_EXPORT extern int
ged_view_polygon_set_visual(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const struct bu_color *edge_color,
	const struct bu_color *fill_color, fastf_t fill_slope_x,
	fastf_t fill_slope_y, fastf_t fill_density, fastf_t vZ, int fill_flag);

GED_EXPORT extern int
ged_view_polygon_set_current(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, long contour_i, long point_i);

GED_EXPORT extern int
ged_view_polygon_set_contour_open(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, long contour_i, int open);

GED_EXPORT extern int
ged_view_polygon_set_all_contours_open(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, int open);

GED_EXPORT extern int
ged_view_polygon_close(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);

GED_EXPORT extern int
ged_view_polygon_clear_selected_point(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);

GED_EXPORT extern int
ged_view_polygon_remove(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);

GED_EXPORT extern void *
ged_view_polygon_user_data(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);

GED_EXPORT extern int
ged_view_polygon_user_data_set(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, void *user_data);

GED_EXPORT extern int
ged_view_polygon_area(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, fastf_t *area);

GED_EXPORT extern int
ged_view_polygon_overlap(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const char *other_name,
	const struct bn_tol *tol, int *overlap);

GED_EXPORT extern int
ged_view_polygon_set_fill(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, int fill_flag, fastf_t fill_slope_x,
	fastf_t fill_slope_y, fastf_t fill_density);

GED_EXPORT extern int
ged_view_polygon_fill_color_get(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, struct bu_color *fill_color);

GED_EXPORT extern int
ged_view_polygon_fill_color_set(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const struct bu_color *fill_color);

GED_EXPORT extern int
ged_view_polygon_csg_name(struct ged_view_context *view_ctx,
	ged_view_polygon_ref target, const char *other_name,
	enum bg_polygon_boolean_op op);

GED_EXPORT extern int
ged_view_polygon_csg(struct ged_view_context *view_ctx,
	ged_view_polygon_ref target, ged_view_polygon_ref stencil,
	enum bg_polygon_boolean_op op);

__END_DECLS

#endif /* GED_VIEW_POLYGON_H */
