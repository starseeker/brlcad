/*       S C E N E _ R E C O R D _ A P I _ I N T E R N A L . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged_scene_record_api_private.h
 *
 * Private realization-record queries used by the semantic reducer and Obol
 * adapter.  Application code must use ged/scene.h and ged/view_export.h.
 */

#ifndef GED_SCENE_RECORD_API_PRIVATE_H
#define GED_SCENE_RECORD_API_PRIVATE_H

#include "common.h"

#include "./ged_draw_types_private.h"
#include "ged/scene.h"
#include "./ged_scene_record_private.h"

__BEGIN_DECLS

extern int
ged_view_selection_add_shape(struct ged *gedp,
	struct ged_view_context *view_ctx, ged_draw_shape_ref ref,
	struct ged_view_context **selection_view_ctx, struct bu_vls *path);

extern int
ged_view_selection_set_highlight(struct ged *gedp,
	struct ged_view_context *view_ctx, ged_draw_shape_ref ref);

extern int
ged_draw_shape_ref_is_null(ged_draw_shape_ref ref);

extern int
ged_draw_group_ref_is_null(ged_draw_group_ref ref);

extern int
ged_draw_shape_record_get(struct ged *gedp,
			  ged_draw_shape_ref ref,
			  struct ged_draw_shape_record *out);


extern int
ged_draw_group_record_get(struct ged *gedp,
			  ged_draw_group_ref ref,
			  struct ged_draw_group_record *out);

extern void
ged_draw_foreach_group_record(struct ged *gedp,
			      int (*cb)(const struct ged_draw_group_record *rec,
					void *userdata),
			      void *userdata);

extern void
ged_draw_foreach_shape_record(struct ged *gedp,
			      int (*cb)(const struct ged_draw_shape_record *rec,
					void *userdata),
			      void *userdata);

extern void
ged_draw_foreach_view_record_query(
    struct ged_view_context *view_ctx,
    const struct ged_draw_view_record_query *query,
    ged_draw_view_db_object_record_cb cb,
    void *userdata);

extern void
ged_draw_foreach_view_db_object_record(struct ged_view_context *view_ctx,
				       ged_draw_view_db_object_record_cb cb,
				       void *userdata);

extern void
ged_draw_foreach_visible_view_db_object_record(struct ged_view_context *view_ctx,
	ged_draw_view_db_object_record_cb cb,
	void *userdata);

extern void
ged_draw_foreach_visible_view_db_object_record_mode(
    struct ged_view_context *view_ctx,
    int draw_mode,
    ged_draw_view_db_object_record_cb cb,
    void *userdata);

extern void
ged_draw_foreach_visible_view_record(struct ged_view_context *view_ctx,
				     ged_draw_view_db_object_record_cb cb,
				     void *userdata);

extern int
ged_draw_view_db_object_record_foreach_segment(
    const struct ged_draw_view_db_object_record *rec,
    ged_draw_view_segment_cb cb,
    void *userdata);

extern int
ged_draw_view_db_object_record_foreach_point(
    const struct ged_draw_view_db_object_record *rec,
    ged_draw_view_point_cb cb,
    void *userdata);

extern int
ged_draw_view_db_object_record_has_segments(
    const struct ged_draw_view_db_object_record *rec);

extern void
ged_draw_view_db_object_record_geometry_report(
    const struct ged_draw_view_db_object_record *rec,
    struct bu_vls *out);

extern int
ged_draw_view_db_object_record_annotation_summary(
    const struct ged_draw_view_db_object_record *rec,
    size_t point_index,
    struct ged_draw_view_annotation_summary *out);

extern int
ged_draw_view_db_object_record_line_summary(
    const struct ged_draw_view_db_object_record *rec,
    struct ged_draw_view_line_summary *out);

extern int
ged_draw_view_db_object_record_line_point_at(
    const struct ged_draw_view_db_object_record *rec,
    size_t index,
    point_t out);

extern int
ged_draw_view_db_object_record_line_command_at(
    const struct ged_draw_view_db_object_record *rec,
    size_t index,
    int *out);

extern int
ged_draw_view_db_object_record_surface_summary(
    const struct ged_draw_view_db_object_record *rec,
    struct ged_draw_view_surface_summary *out);

extern int
ged_draw_view_db_object_record_surface_index_at(
    const struct ged_draw_view_db_object_record *rec,
    size_t index,
    int *out);

extern int
ged_draw_view_rendered_object_summary(
    struct ged_view_context *view_ctx,
    uint64_t cache_identity,
    ged_draw_view_rendered_object_summary_t *out);

extern ged_draw_shape_ref
ged_draw_first_shape_ref(struct ged *gedp);

extern ged_draw_shape_ref
ged_draw_shape_ref_at(struct ged *gedp, int idx);

extern int
ged_draw_shape_ref_index(struct ged *gedp, ged_draw_shape_ref ref);

extern int
ged_draw_shape_ref_index_for_component(struct ged *gedp,
				       const char *path,
				       ged_draw_shape_ref_index_cb cb,
				       void *userdata);

extern int
ged_draw_shape_ref_index_for_path_hash(struct ged *gedp,
				       unsigned long long path_hash,
				       ged_draw_shape_ref_index_cb cb,
				       void *userdata);

extern void
ged_draw_foreach_visible_shape_candidate(struct ged *gedp,
					 ged_draw_shape_candidate_cb cb,
					 void *userdata);

extern ged_draw_shape_ref
ged_draw_shape_ref_for_candidate(
    struct ged *gedp,
    const struct ged_draw_shape_candidate *candidate);

extern ged_draw_shape_ref
ged_draw_advance_shape_ref(struct ged *gedp, ged_draw_shape_ref ref, int delta);

extern ged_draw_group_ref
ged_draw_group_ref_of_shape(struct ged *gedp, ged_draw_shape_ref ref);

extern struct ged_view_context *
ged_draw_shape_ref_view_context(struct ged *gedp, ged_draw_shape_ref ref);

extern int
ged_draw_shape_ref_get_color(struct ged *gedp, ged_draw_shape_ref ref, unsigned char rgb[3]);

extern int
ged_draw_shape_ref_display_summary(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   struct ged_draw_scene_display_summary *out);

extern int
ged_draw_shape_ref_material_summary(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    struct ged_draw_shape_material_summary *out);

extern int
ged_draw_shape_ref_last_point(struct ged *gedp,
			      ged_draw_shape_ref ref,
			      point_t out);

extern int
ged_draw_shape_ref_line_summary(struct ged *gedp,
				ged_draw_shape_ref ref,
				struct ged_draw_view_line_summary *out);

extern int
ged_draw_shape_ref_line_point_at(struct ged *gedp,
				 ged_draw_shape_ref ref,
				 size_t index,
				 point_t out);

extern int
ged_draw_shape_ref_line_command_at(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   size_t index,
				   int *out);

extern int
ged_draw_shape_ref_geometry_summary(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    struct ged_draw_shape_geometry_summary *out);

extern uint64_t
ged_scene_frame_revision(struct ged_view_context *view_ctx);

extern uint64_t
ged_scene_bump_frame_revision(struct ged_view_context *view_ctx);


/**
 * Return non-zero when @p rec belongs to the drawable database scene for
 * @p v.  Shared draw groups match non-independent views; independent views
 * match only groups in their private view scope.  Passing NULL for @p v
 * queries only shared draw groups.
 *
 * This helper answers view-scope membership only.  Callers that list renderable
 * database paths must still check record visibility, overlay status, and draw
 * mode as appropriate.
 */
extern int
ged_draw_group_record_in_view(const struct ged_draw_group_record *rec,
			      struct ged_view_context *view_ctx);

/* Private mutation hooks retained only for low-level adapter diagnostics.
 * Interactive clients publish retained edit previews instead. */
extern int ged_scene_internal_shape_translate_geometry(
    struct ged *gedp, ged_draw_shape_ref ref, const vect_t translation);
extern int ged_scene_internal_shape_set_center(
    struct ged *gedp, ged_draw_shape_ref ref, const point_t center);
extern int ged_scene_internal_shape_geometry_clear(
    struct ged *gedp, ged_draw_shape_ref ref);
extern int ged_scene_internal_shape_bounds_update(
    struct ged *gedp, ged_draw_shape_ref ref, int *bad_command);

__END_DECLS

#endif /* GED_SCENE_RECORD_API_PRIVATE_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
