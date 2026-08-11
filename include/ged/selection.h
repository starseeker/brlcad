/* S E L E C T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged/selection.h */

#ifndef GED_SELECTION_H
#define GED_SELECTION_H

#include "ged/selection_types.h"
#include "ged/scene_types.h"

__BEGIN_DECLS

GED_EXPORT extern size_t
ged_view_selection_count(struct ged_view_context *view_ctx);

GED_EXPORT extern int
ged_view_selection_visit(
    struct ged_view_context *view_ctx,
    ged_view_selection_visit_cb cb,
    void *data);

GED_EXPORT extern int
ged_view_selection_clear(struct ged_view_context *view_ctx);

/** Add one exact semantic occurrence to the selected-path set. */
GED_EXPORT extern int
ged_view_selection_add_occurrence(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    ged_scene_occurrence_ref occurrence,
    struct ged_view_context **selection_view_ctx,
    struct bu_vls *path);

/** Replace the transient highlighted-reference selection with an occurrence. */
GED_EXPORT extern int
ged_view_selection_highlight_occurrence_set(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    ged_scene_occurrence_ref occurrence);

GED_EXPORT extern int
ged_view_selection_set_pick(
    struct ged_view_context *view_ctx,
    const struct ged_pick_result *result,
    ged_view_selection_path_cb cb,
    void *data);

GED_EXPORT extern struct ged_pick_result *
ged_pick_result_create(void);

GED_EXPORT extern void
ged_pick_result_free(struct ged_pick_result *result);

GED_EXPORT extern size_t
ged_pick_result_count(const struct ged_pick_result *result);

GED_EXPORT extern int
ged_pick_result_path(const struct ged_pick_result *result,
			  size_t index,
			  struct bu_vls *path_out);

GED_EXPORT extern fastf_t
ged_pick_result_hit_dist(const struct ged_pick_result *result,
			      size_t index);

struct ged_pick_detail {
    uint32_t source_id;
    int primitive_kind;
    int primitive_index;
    int material_id;
    int face_vertex_index[3];
    int nearest_face_vertex_index;
    point_t model_point;
    int model_point_valid;
};

#define GED_PICK_DETAIL_INIT { 0, 0, -1, 0, {-1, -1, -1}, -1, \
    VINIT_ZERO, 0 }

GED_EXPORT extern int
ged_pick_result_detail(const struct ged_pick_result *result,
	size_t index, struct ged_pick_detail *detail);

GED_EXPORT extern int
ged_pick_result_append_detail(struct ged_pick_result *result,
	const char *path, fastf_t hit_dist,
	const struct ged_pick_detail *detail);

GED_EXPORT extern int
ged_pick_result_append_path(struct ged_pick_result *result,
				 const char *path,
				 fastf_t hit_dist);

GED_EXPORT extern int
ged_pick_result_append_copy(struct ged_pick_result *dest,
				 const struct ged_pick_result *src,
				 size_t index,
				 fastf_t hit_dist);

GED_EXPORT extern struct ged_pick_result *
ged_pick_result_filter_first(const struct ged_pick_result *src);

GED_EXPORT extern struct ged_pick_result *
ged_pick_point(struct ged_view_context *view_ctx,
				 int x,
				 int y,
				 int first_only);

GED_EXPORT extern struct ged_pick_result *
ged_pick_nearest(struct ged_view_context *view_ctx,
				   int x,
				   int y);

GED_EXPORT extern struct ged_pick_result *
ged_pick_rect(struct ged_view_context *view_ctx,
				int x0,
				int y0,
				int x1,
				int y1);

GED_EXPORT extern struct ged_pick_result *
ged_pick_semantic_path(struct ged *gedp,
					 struct ged_view_context *view_ctx,
					 const char *path_pattern);

GED_EXPORT extern int
ged_view_selection_snap(struct ged_view_context *view_ctx,
	const point_t sample,
	enum ged_selection_snap_kind kind,
	point_t candidate);


__END_DECLS

#endif /* GED_SELECTION_H */
