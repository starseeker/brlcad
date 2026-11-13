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

/** Return non-zero if @p ref is the null polygon reference. */
GED_EXPORT extern int
ged_view_polygon_ref_is_null(ged_view_polygon_ref ref);

/** Find the polygon named @p name in the effective view scope. */
GED_EXPORT extern ged_view_polygon_ref
ged_view_polygon_find(struct ged_view_context *view_ctx, const char *name);

/** Find a named polygon, optionally restricting the search to local state. */
GED_EXPORT extern ged_view_polygon_ref
ged_view_polygon_find_scoped(struct ged_view_context *view_ctx,
	const char *name, int local_only);

/** Find the polygon control point nearest @p model_point. */
GED_EXPORT extern ged_view_polygon_ref
ged_view_polygon_select(struct ged_view_context *view_ctx,
	const point_t model_point);

/** Create a view-owned polygon and return its generation-checked reference. */
GED_EXPORT extern ged_view_polygon_ref
ged_view_polygon_create(struct ged_view_context *view_ctx, const char *name,
	int local, enum ged_view_polygon_type type, const point_t model_point);

/** Duplicate a named polygon under @p new_name. */
GED_EXPORT extern ged_view_polygon_ref
ged_view_polygon_dup(struct ged_view_context *view_ctx, const char *name,
	const char *new_name);

/** Import @p dp as a view-owned polygon backed by sketch geometry. */
GED_EXPORT extern ged_view_polygon_ref
ged_view_polygon_import_sketch(struct ged_view_context *view_ctx,
	const char *name, struct db_i *dbip, struct directory *dp, int local);

/** Visit immutable snapshots of all polygons in the effective view scope. */
GED_EXPORT extern void
ged_view_polygon_visit_records(struct ged_view_context *view_ctx,
	ged_view_polygon_record_cb callback, void *data);

/** Return the number of polygon snap candidates other than @p exclude. */
GED_EXPORT extern size_t
ged_view_polygon_snap_count(struct ged_view_context *view_ctx,
	ged_view_polygon_ref exclude);

/** Clear selected control points from all polygons in the view. */
GED_EXPORT extern int
ged_view_polygon_clear_point_selection(struct ged_view_context *view_ctx);

/** Clear whole-polygon selection from all polygons in the view. */
GED_EXPORT extern int
ged_view_polygon_clear_selection(struct ged_view_context *view_ctx);

/** Exclude @p ref from subsequent polygon snap candidate queries. */
GED_EXPORT extern int
ged_view_polygon_snap_exclude_set(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);

/** Export @p ref to a newly created sketch directory entry. */
GED_EXPORT extern struct directory *
ged_view_polygon_export_sketch(struct ged_view_context *view_ctx,
	struct db_i *dbip, const char *name, ged_view_polygon_ref ref);

/** Replace the linked existing sketch with the current polygon state. */
GED_EXPORT extern struct directory *
ged_view_polygon_update_sketch(struct ged_view_context *view_ctx,
	struct db_i *dbip, ged_view_polygon_ref ref);

/** Set the stable database sketch name synchronized with @p ref. */
GED_EXPORT extern int
ged_view_polygon_sketch_name_set(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const char *name);

/**
 * Synchronize the database sketch linked to @p ref with its current polygon
 * state and publish the corresponding GED database event.
 *
 * Return -1 on error, 0 when the polygon has no linked sketch, 1 when a
 * missing sketch was created, and 2 when an existing sketch was updated.
 */
GED_EXPORT extern int
ged_view_polygon_sync_sketch(struct ged *gedp,
	struct ged_view_context *view_ctx, ged_view_polygon_ref ref);

/** Copy the current immutable state of @p ref into @p record. */
GED_EXPORT extern int
ged_view_polygon_record_get(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, struct ged_view_polygon_record *record);

/** Return non-zero when @p ref contains usable polygon geometry. */
GED_EXPORT extern int
ged_view_polygon_has_data(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);

/** Apply a non-positional polygon edit operation. */
GED_EXPORT extern int
ged_view_polygon_update(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, enum ged_view_polygon_update op);

/** Apply a screen-positioned polygon edit operation. */
GED_EXPORT extern int
ged_view_polygon_update_screen_pt(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, int x, int y, enum ged_view_polygon_update op);

/** Apply a model-positioned polygon edit operation. */
GED_EXPORT extern int
ged_view_polygon_update_model_pt(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const point_t model_point,
	enum ged_view_polygon_update op);

/** Move @p ref by the model-space delta between two sample points. */
GED_EXPORT extern int
ged_view_polygon_move(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const point_t current_point,
	const point_t previous_point);

/** Rename @p ref within its owning view scope. */
GED_EXPORT extern int
ged_view_polygon_set_name(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const char *name);

/** Replace the complete visual presentation of @p ref. */
GED_EXPORT extern int
ged_view_polygon_set_visual(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const struct bu_color *edge_color,
	const struct bu_color *fill_color, fastf_t fill_slope_x,
	fastf_t fill_slope_y, fastf_t fill_density, fastf_t vZ, int fill_flag);

/** Set the active contour and point used by interactive editing. */
GED_EXPORT extern int
ged_view_polygon_set_current(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, long contour_i, long point_i);

/** Set whole-polygon selection presentation for @p ref. */
GED_EXPORT extern int
ged_view_polygon_set_selected(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, int selected);

/** Set whether one contour is open. */
GED_EXPORT extern int
ged_view_polygon_set_contour_open(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, long contour_i, int open);

/** Set the open state of every contour in @p ref. */
GED_EXPORT extern int
ged_view_polygon_set_all_contours_open(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, int open);

/** Close the active contour in @p ref. */
GED_EXPORT extern int
ged_view_polygon_close(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);

/** Clear the selected control point in @p ref. */
GED_EXPORT extern int
ged_view_polygon_clear_selected_point(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);

/** Remove @p ref from its owning view. */
GED_EXPORT extern int
ged_view_polygon_remove(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);

/** Return the borrowed application data associated with @p ref. */
GED_EXPORT extern void *
ged_view_polygon_user_data(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);

/** Associate borrowed application data with @p ref. */
GED_EXPORT extern int
ged_view_polygon_user_data_set(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, void *user_data);

/** Calculate the model-space area of @p ref. */
GED_EXPORT extern int
ged_view_polygon_area(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, fastf_t *area);

/** Test @p ref for overlap with the named polygon. */
GED_EXPORT extern int
ged_view_polygon_overlap(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const char *other_name,
	const struct bn_tol *tol, int *overlap);

/** Set fill visibility and hatch parameters for @p ref. */
GED_EXPORT extern int
ged_view_polygon_set_fill(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, int fill_flag, fastf_t fill_slope_x,
	fastf_t fill_slope_y, fastf_t fill_density);

/** Copy the fill color of @p ref into @p fill_color. */
GED_EXPORT extern int
ged_view_polygon_fill_color_get(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, struct bu_color *fill_color);

/** Set the fill color of @p ref. */
GED_EXPORT extern int
ged_view_polygon_fill_color_set(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref, const struct bu_color *fill_color);

/** Apply a polygon Boolean using a named stencil polygon. */
GED_EXPORT extern int
ged_view_polygon_csg_name(struct ged_view_context *view_ctx,
	ged_view_polygon_ref target, const char *other_name,
	enum bg_polygon_boolean_op op);

/** Apply a polygon Boolean between two owner-checked references. */
GED_EXPORT extern int
ged_view_polygon_csg(struct ged_view_context *view_ctx,
	ged_view_polygon_ref target, ged_view_polygon_ref stencil,
	enum bg_polygon_boolean_op op);

__END_DECLS

#endif /* GED_VIEW_POLYGON_H */
