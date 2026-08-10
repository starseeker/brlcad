/* P O L Y G O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged/polygon.h */

#ifndef GED_POLYGON_H
#define GED_POLYGON_H

#include "ged/polygon_types.h"

__BEGIN_DECLS

GED_EXPORT extern int
ged_polygon_ref_is_null(ged_polygon_ref ref);

GED_EXPORT extern ged_polygon_ref
ged_polygon_find(struct ged_view_context *view_ctx,
				   const char *name);

GED_EXPORT extern ged_polygon_ref
ged_polygon_find_scoped(struct ged_view_context *view_ctx,
	const char *name,
	int local_only);

GED_EXPORT extern ged_polygon_ref
ged_polygon_select(struct ged_view_context *view_ctx,
				     const point_t model_point);

GED_EXPORT extern ged_polygon_ref
ged_polygon_create(struct ged_view_context *view_ctx,
				     const char *name,
				     int local,
				     int type,
				     const point_t screen_point);

GED_EXPORT extern ged_polygon_ref
ged_polygon_dup(struct ged_view_context *view_ctx,
				  const char *name,
				  const char *new_name);

GED_EXPORT extern ged_polygon_ref
ged_polygon_import_sketch(const char *name,
	struct db_i *dbip,
	struct directory *dp,
	struct ged_view_context *view_ctx,
	int local);

GED_EXPORT extern void
ged_polygon_visit_records(
	struct ged_view_context *view_ctx,
	ged_polygon_record_cb callback,
	void *data);

GED_EXPORT extern size_t
ged_polygon_snap_count(struct ged_view_context *view_ctx,
	ged_polygon_ref exclude);

GED_EXPORT extern int
ged_polygon_clear_point_selection(struct ged_view_context *view_ctx);

GED_EXPORT extern int
ged_polygon_snap_exclude_set(
	struct ged_view_context *view_ctx,
	ged_polygon_ref ref);

GED_EXPORT extern struct directory *
ged_polygon_export_sketch(struct db_i *dbip,
				    const char *name,
				    ged_polygon_ref ref);

GED_EXPORT extern int
ged_polygon_record_get(ged_polygon_ref ref,
				 struct ged_polygon_record *record);

GED_EXPORT extern int
ged_polygon_has_data(ged_polygon_ref ref);

GED_EXPORT extern int
ged_polygon_update(ged_polygon_ref ref,
				     struct ged_view_context *view_ctx,
				     int op);

GED_EXPORT extern int
ged_polygon_update_screen_pt(ged_polygon_ref ref,
	struct ged_view_context *view_ctx,
	int x,
	int y,
	int op);

GED_EXPORT extern int
ged_polygon_move(ged_polygon_ref ref,
			   point_t *current_point,
			   point_t *previous_point);

GED_EXPORT extern int
ged_polygon_set_name(ged_polygon_ref ref,
			       const char *name);

GED_EXPORT extern int
ged_polygon_set_visual(ged_polygon_ref ref,
				 const struct bu_color *edge_color,
				 const struct bu_color *fill_color,
				 fastf_t fill_slope_x,
				 fastf_t fill_slope_y,
				 fastf_t fill_density,
				 fastf_t vZ,
				 int fill_flag);

GED_EXPORT extern int
ged_polygon_set_current(ged_polygon_ref ref,
				  long contour_i,
				  long point_i);

GED_EXPORT extern int
ged_polygon_set_contour_open(ged_polygon_ref ref,
				       long contour_i,
				       int open);

GED_EXPORT extern int
ged_polygon_set_all_contours_open(ged_polygon_ref ref,
	int open);

GED_EXPORT extern int
ged_polygon_close(ged_polygon_ref ref);

GED_EXPORT extern int
ged_polygon_clear_selected_point(ged_polygon_ref ref);

GED_EXPORT extern int
ged_polygon_remove(ged_polygon_ref ref);

GED_EXPORT extern void *
ged_polygon_user_data(ged_polygon_ref ref);

GED_EXPORT extern int
ged_polygon_user_data_set(ged_polygon_ref ref,
				    void *user_data);

GED_EXPORT extern int
ged_polygon_area(ged_polygon_ref ref,
				   struct ged_view_context *view_ctx,
				   fastf_t *area);

GED_EXPORT extern int
ged_polygon_overlap(ged_polygon_ref ref,
				      struct ged_view_context *view_ctx,
				      const char *other_name,
				      const struct bn_tol *tol,
				      int *overlap);

GED_EXPORT extern int
ged_polygon_set_fill(ged_polygon_ref ref,
			       int fill_flag,
			       fastf_t fill_slope_x,
			       fastf_t fill_slope_y,
			       fastf_t fill_density);

GED_EXPORT extern int
ged_polygon_fill_color_get(ged_polygon_ref ref,
				     struct bu_color *fill_color);

GED_EXPORT extern int
ged_polygon_fill_color_set(ged_polygon_ref ref,
				     const struct bu_color *fill_color);

GED_EXPORT extern int
ged_polygon_csg_name(ged_polygon_ref target,
				  struct ged_view_context *view_ctx,
				  const char *other_name,
				  bg_clip_t op);

GED_EXPORT extern int
ged_polygon_csg(ged_polygon_ref target,
			  ged_polygon_ref stencil,
			  bg_clip_t op);

__END_DECLS

#endif /* GED_POLYGON_H */
