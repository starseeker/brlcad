/* S C E N E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged/scene.h */

#ifndef GED_SCENE_H
#define GED_SCENE_H

#include "common.h"
#include "ged/draw_types.h"

__BEGIN_DECLS

GED_EXPORT extern struct ged_draw_transaction
ged_draw_transaction_make(ged_draw_transaction_kind kind, const char *path);

GED_EXPORT extern struct ged_draw_transaction
ged_draw_transaction_make_value(ged_draw_transaction_kind kind,
				const char *path,
				double value);

GED_EXPORT extern void
ged_draw_transaction_result_init(struct ged_draw_transaction_result *result);

GED_EXPORT extern void
ged_draw_transaction_result_free(struct ged_draw_transaction_result *result);

GED_EXPORT extern int
ged_draw_apply_transaction(struct ged *gedp,
			   const struct ged_draw_transaction *txn,
			   struct ged_draw_transaction_result *result);

/**
 * Subscribe to post-transaction draw-state events.
 *
 * Observers are called synchronously after a successful draw transaction has
 * reconciled libged draw state and populated a stable
 * ged_draw_transaction_result.  The returned token may be passed to
 * ged_draw_observer_remove(); 0 indicates registration failure.
 */
GED_EXPORT extern ged_draw_observer_token
ged_draw_observer_add(struct ged *gedp,
		      ged_draw_transaction_observer_func_t func,
		      void *client_data);

/**
 * Remove a draw-state observer.  Safe to call from inside an observer callback.
 * Returns 1 when a live observer was removed, 0 otherwise.
 */
GED_EXPORT extern int
ged_draw_observer_remove(struct ged *gedp, ged_draw_observer_token token);


GED_EXPORT extern int
ged_draw_scene_available(struct ged *gedp);

/**
 * Structured-path convenience wrappers for common command operations.
 */
GED_EXPORT extern int
ged_draw_apply_erase_path(struct ged *gedp,
			  const struct db_full_path *dfp);

GED_EXPORT extern int
ged_draw_apply_erase_path_prefix(struct ged *gedp,
				 const struct db_full_path *dfp);

/**
 * Return the current default draw mode used by draw commands that do not
 * specify an explicit render style override.
 */
GED_EXPORT extern int
ged_draw_default_mode(const struct ged *gedp);

/**
 * Erase from @p gedp's drawn scene the entry whose path equals @p dfp.
 * When a scene group's path is a strict ancestor of the erase path, the
 * matching sub-tree is removed without disturbing sibling sub-groups.
 */
GED_EXPORT extern void
ged_draw_erase_path(struct ged *gedp,
		    const struct db_full_path *dfp);

/**
 * Erase from @p gedp's drawn-object set every scene object whose path
 * contains @p name as one of its directory components.
 *
 * This is used when a database object is renamed or deleted and all drawn
 * references to that name must be removed.
 */
GED_EXPORT extern void
ged_draw_erase_name(struct ged *gedp, const char *name);

/**
 * Erase every drawn scene object whose path has @p dfp as a prefix subset.
 */
GED_EXPORT extern void
ged_draw_erase_path_prefix(struct ged *gedp,
			   const struct db_full_path *dfp);

/**
 * Compute the axis-aligned bounding box of the visible drawn scene in
 * @p gedp's active view set.  When @p pflag is zero, non-database overlay
 * records are excluded.
 * Returns 1 if the result is empty (no contributing objects), 0 otherwise.
 *
 * The bounds are returned in model coordinates.
 */
GED_EXPORT extern int
ged_draw_bounds(struct ged *gedp, vect_t *min, vect_t *max,
		int pflag);

/**
 * Set the draw-scene highlight state for every drawn shape.  A non-zero
 * @p highlighted value asserts highlight state; zero clears it.
 */
GED_EXPORT extern void
ged_draw_set_highlight_state(struct ged *gedp, int highlighted);

/**
 * Set highlight state on all draw shapes whose database full path begins with
 * @p path[0..path_pos].  This is the semantic operation behind MGED matrix
 * path picking, where multiple rendered shapes can share the same editable
 * path prefix.
 *
 * Returns the number of matching shape records.
 */
GED_EXPORT extern int
ged_draw_set_highlighted_path_prefix(struct ged *gedp,
				     const struct db_full_path *path,
				     size_t path_pos,
				     int highlighted);

/**
 * Return the highlight-state revision counter.  Bumped on every transition
 * of the highlighted draw ref.  Cache a snapshot and compare to detect
 * "highlight may have changed since I last looked" cheaply.
 */
GED_EXPORT extern uint64_t
ged_draw_highlight_revision(const struct ged *gedp);

/**
 * Return the material revision counter.  The counter is incremented by
 * ged_draw_bump_material_revision() whenever the material/color table changes.
 * ged_draw_refresh_material_colors() does NOT increment the counter; it only
 * stamps per-shape color revisions.  Callers that cache per-shape color
 * data can compare a saved snapshot to detect when a recolor sweep is needed.
 */
GED_EXPORT extern uint64_t
ged_draw_material_revision(const struct ged *gedp);

/**
 * Bump the material revision counter.
 *
 * Call this after any operation that changes the effective material or
 * color table (e.g. 'color', 'mater', 'rmater', 'edmater') so that the
 * next ged_draw_refresh_material_colors() call recolors shapes that are
 * now stale.
 */
GED_EXPORT extern void
ged_draw_bump_material_revision(struct ged *gedp);

/**
 * Refresh per-object base color from the dbip's region/material table
 * (mater_struct chain) for every drawn scene object whose color stamp is
 * stale relative to the material revision counter.  Shapes that
 * were already colored since the last material-change event are skipped.
 *
 * Callers must invoke ged_draw_bump_material_revision() before this function
 * to signal that a material change occurred; without that bump, successive
 * calls will skip all already-stamped shapes.
 */
GED_EXPORT extern void
ged_draw_refresh_material_colors(struct ged *gedp);

/**
 * Return the current draw-scene structural revision as an unsigned long long.
 * Returns 0 when the scene is unavailable or after a clear operation.
 *
 * Prefer ged_draw_scene_revision() for new code that only needs change
 * detection.
 */
GED_EXPORT extern unsigned long long
ged_draw_scene_hash(struct ged *gedp);

/**
 * Return the raw structural revision counter for @p gedp's draw tree.
 * The counter is incremented on every structural mutation (group/shape
 * addition or removal) and reset to 0 by ged_draw_clear().  Returns
 * 0 when no objects have been drawn since the last zap (or ever).
 *
 * Callers that need to detect "drawn set changed" cheaply should compare
 * snapshots of this value rather than recomputing a content hash.
 */
GED_EXPORT extern uint64_t
ged_draw_scene_revision(struct ged *gedp);

GED_EXPORT extern int
ged_draw_shape_ref_is_null(ged_draw_shape_ref ref);

GED_EXPORT extern int
ged_draw_group_ref_is_null(ged_draw_group_ref ref);

GED_EXPORT extern int
ged_scene_node_ref_is_null(ged_scene_node_ref ref);

GED_EXPORT extern int
ged_scene_node_ref_equal(ged_scene_node_ref a, ged_scene_node_ref b);

GED_EXPORT extern int
ged_draw_shape_record_get(struct ged *gedp,
			  ged_draw_shape_ref ref,
			  struct ged_draw_shape_record *out);


GED_EXPORT extern int
ged_draw_group_record_get(struct ged *gedp,
			  ged_draw_group_ref ref,
			  struct ged_draw_group_record *out);

GED_EXPORT extern void
ged_draw_foreach_group_record(struct ged *gedp,
			      int (*cb)(const struct ged_draw_group_record *rec,
					void *userdata),
			      void *userdata);

GED_EXPORT extern void
ged_draw_foreach_shape_record(struct ged *gedp,
			      int (*cb)(const struct ged_draw_shape_record *rec,
					void *userdata),
			      void *userdata);

GED_EXPORT extern void
ged_draw_foreach_view_record_query(
    struct ged_view_context *view_ctx,
    const struct ged_draw_view_record_query *query,
    ged_draw_view_db_object_record_cb cb,
    void *userdata);

GED_EXPORT extern void
ged_draw_foreach_view_db_object_record(struct ged_view_context *view_ctx,
				       ged_draw_view_db_object_record_cb cb,
				       void *userdata);

GED_EXPORT extern void
ged_draw_foreach_visible_view_db_object_record(struct ged_view_context *view_ctx,
	ged_draw_view_db_object_record_cb cb,
	void *userdata);

GED_EXPORT extern void
ged_draw_foreach_visible_view_db_object_record_mode(
    struct ged_view_context *view_ctx,
    int draw_mode,
    ged_draw_view_db_object_record_cb cb,
    void *userdata);

GED_EXPORT extern void
ged_draw_foreach_visible_view_record(struct ged_view_context *view_ctx,
				     ged_draw_view_db_object_record_cb cb,
				     void *userdata);

GED_EXPORT extern int
ged_draw_view_db_object_record_foreach_segment(
    const struct ged_draw_view_db_object_record *rec,
    ged_draw_view_segment_cb cb,
    void *userdata);

GED_EXPORT extern int
ged_draw_view_db_object_record_foreach_point(
    const struct ged_draw_view_db_object_record *rec,
    ged_draw_view_point_cb cb,
    void *userdata);

GED_EXPORT extern int
ged_draw_view_db_object_record_has_segments(
    const struct ged_draw_view_db_object_record *rec);

GED_EXPORT extern void
ged_draw_view_db_object_record_geometry_report(
    const struct ged_draw_view_db_object_record *rec,
    struct bu_vls *out);

GED_EXPORT extern int
ged_draw_view_db_object_record_annotation_summary(
    const struct ged_draw_view_db_object_record *rec,
    size_t point_index,
    struct ged_draw_view_annotation_summary *out);

GED_EXPORT extern int
ged_draw_view_db_object_record_line_summary(
    const struct ged_draw_view_db_object_record *rec,
    struct ged_draw_view_line_summary *out);

GED_EXPORT extern int
ged_draw_view_db_object_record_line_point_at(
    const struct ged_draw_view_db_object_record *rec,
    size_t index,
    point_t out);

GED_EXPORT extern int
ged_draw_view_db_object_record_line_command_at(
    const struct ged_draw_view_db_object_record *rec,
    size_t index,
    int *out);

GED_EXPORT extern int
ged_draw_view_db_object_record_surface_summary(
    const struct ged_draw_view_db_object_record *rec,
    struct ged_draw_view_surface_summary *out);

GED_EXPORT extern int
ged_draw_view_db_object_record_surface_index_at(
    const struct ged_draw_view_db_object_record *rec,
    size_t index,
    int *out);

GED_EXPORT extern int
ged_draw_view_rendered_object_summary(
    struct ged_view_context *view_ctx,
    uint64_t cache_identity,
    ged_draw_view_rendered_object_summary_t *out);

GED_EXPORT extern ged_draw_shape_ref
ged_draw_first_shape_ref(struct ged *gedp);

GED_EXPORT extern ged_draw_shape_ref
ged_draw_shape_ref_at(struct ged *gedp, int idx);

GED_EXPORT extern int
ged_draw_shape_ref_index(struct ged *gedp, ged_draw_shape_ref ref);

GED_EXPORT extern int
ged_draw_shape_ref_index_for_component(struct ged *gedp,
				       const char *path,
				       ged_draw_shape_ref_index_cb cb,
				       void *userdata);

GED_EXPORT extern int
ged_draw_shape_ref_index_for_path_hash(struct ged *gedp,
				       unsigned long long path_hash,
				       ged_draw_shape_ref_index_cb cb,
				       void *userdata);

GED_EXPORT extern void
ged_draw_foreach_visible_shape_candidate(struct ged *gedp,
					 ged_draw_shape_candidate_cb cb,
					 void *userdata);

GED_EXPORT extern ged_draw_shape_ref
ged_draw_shape_ref_for_candidate(
    struct ged *gedp,
    const struct ged_draw_shape_candidate *candidate);

GED_EXPORT extern ged_draw_shape_ref
ged_draw_advance_shape_ref(struct ged *gedp, ged_draw_shape_ref ref, int delta);

GED_EXPORT extern ged_draw_group_ref
ged_draw_group_ref_of_shape(struct ged *gedp, ged_draw_shape_ref ref);

GED_EXPORT extern ged_scene_node_ref
ged_scene_first_node(struct ged *gedp);

GED_EXPORT extern ged_scene_node_ref
ged_scene_shape_node(struct ged *gedp, ged_draw_shape_ref ref);

GED_EXPORT extern ged_scene_node_ref
ged_scene_shape_cache_node(struct ged *gedp, ged_draw_shape_ref ref);

GED_EXPORT extern struct ged_view_context *
ged_draw_shape_ref_view_context(struct ged *gedp, ged_draw_shape_ref ref);

GED_EXPORT extern ged_scene_node_ref
ged_scene_group_node(struct ged *gedp, ged_draw_group_ref ref);

GED_EXPORT extern ged_draw_shape_ref
ged_scene_shape_ref(struct ged *gedp, ged_scene_node_ref node);

GED_EXPORT extern void
ged_draw_set_highlighted_shape_ref(struct ged *gedp, ged_draw_shape_ref ref);

GED_EXPORT extern int
ged_draw_shape_set_highlighted(struct ged *gedp, ged_draw_shape_ref ref, int highlighted);

GED_EXPORT extern int
ged_draw_shape_ref_get_color(struct ged *gedp, ged_draw_shape_ref ref, unsigned char rgb[3]);

GED_EXPORT extern int
ged_draw_shape_ref_set_color(struct ged *gedp, ged_draw_shape_ref ref, const unsigned char rgb[3]);

GED_EXPORT extern int
ged_draw_shape_ref_set_visible(struct ged *gedp,
			       ged_draw_shape_ref ref,
			       int visible);

GED_EXPORT extern int
ged_draw_shape_ref_display_summary(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   struct ged_draw_scene_display_summary *out);

GED_EXPORT extern int
ged_draw_shape_ref_material_summary(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    struct ged_draw_shape_material_summary *out);

GED_EXPORT extern int
ged_draw_shape_ref_set_material_color(struct ged *gedp,
				      ged_draw_shape_ref ref,
				      const unsigned char rgb[3]);

GED_EXPORT extern int
ged_draw_shape_ref_set_evaluated_region(struct ged *gedp,
					ged_draw_shape_ref ref,
					int evaluated_region);

GED_EXPORT extern int
ged_draw_shape_ref_last_point(struct ged *gedp,
			      ged_draw_shape_ref ref,
			      point_t out);

GED_EXPORT extern int
ged_draw_shape_ref_line_summary(struct ged *gedp,
				ged_draw_shape_ref ref,
				struct ged_draw_view_line_summary *out);

GED_EXPORT extern int
ged_draw_shape_ref_line_point_at(struct ged *gedp,
				 ged_draw_shape_ref ref,
				 size_t index,
				 point_t out);

GED_EXPORT extern int
ged_draw_shape_ref_line_command_at(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   size_t index,
				   int *out);

GED_EXPORT extern int
ged_draw_shape_ref_geometry_summary(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    struct ged_draw_shape_geometry_summary *out);

GED_EXPORT extern int
ged_scene_node_line_summary(struct ged *gedp, ged_scene_node_ref node,
				    struct ged_draw_view_line_summary *out);

GED_EXPORT extern int
ged_scene_node_line_point_at(struct ged *gedp, ged_scene_node_ref node,
				     size_t index,
				     point_t out);

GED_EXPORT extern int
ged_scene_node_line_command_at(struct ged *gedp, ged_scene_node_ref node,
				       size_t index,
				       int *out);

GED_EXPORT extern int
ged_scene_node_geometry_summary(struct ged *gedp, ged_scene_node_ref node,
					struct ged_draw_shape_geometry_summary *out);

GED_EXPORT extern int
ged_scene_node_has_state(struct ged *gedp, ged_scene_node_ref node);

GED_EXPORT extern ged_scene_node_ref
ged_scene_node_source(struct ged *gedp, ged_scene_node_ref node);

GED_EXPORT extern int
ged_annotation_hud_sync(struct ged_view_context *view_ctx);

struct ged_scene_export_report {
    int export_record_found;
    int render_item_found;
    int backend_node_found;
    int export_render_consistent;
    int export_backend_consistent;
};

#define GED_SCENE_EXPORT_CONSISTENCY_INIT { 0, 0, 0, 0, 0 }

GED_EXPORT extern int
ged_scene_check_export(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char *drawn_prefix,
    struct ged_scene_export_report *summary);

GED_EXPORT extern uint64_t
ged_scene_frame_revision(struct ged_view_context *view_ctx);

GED_EXPORT extern uint64_t
ged_scene_bump_frame_revision(struct ged_view_context *view_ctx);


GED_EXPORT extern const struct db_full_path *
ged_scene_node_fullpath(struct ged *gedp, ged_scene_node_ref node);

GED_EXPORT extern int
ged_scene_group_dbpath(struct ged *gedp,
			      ged_scene_node_ref group,
			      struct db_full_path *out);

GED_EXPORT extern int
ged_scene_group_is_overlay(struct ged *gedp, ged_scene_node_ref group);

GED_EXPORT extern int
ged_scene_node_display_summary(struct ged *gedp, ged_scene_node_ref node,
				       struct ged_draw_scene_display_summary *out);


GED_EXPORT extern int
ged_scene_node_tree_summary(struct ged *gedp, ged_scene_node_ref node,
				    struct ged_draw_scene_tree_summary *out);

GED_EXPORT extern ged_scene_node_ref
ged_scene_node_child_at(struct ged *gedp, ged_scene_node_ref node, size_t index);

GED_EXPORT extern ged_scene_node_ref
ged_scene_node_parent(struct ged *gedp, ged_scene_node_ref node);

GED_EXPORT extern const char *
ged_scene_node_name(struct ged *gedp, ged_scene_node_ref node);

GED_EXPORT extern int
ged_scene_node_subtree_bounds(struct ged *gedp, ged_scene_node_ref node,
				      vect_t *min,
				      vect_t *max,
				      int include_overlays);

GED_EXPORT extern int
ged_draw_shape_ref_translate_geometry(struct ged *gedp,
				      ged_draw_shape_ref ref,
				      const vect_t xlate);

GED_EXPORT extern int
ged_draw_shape_ref_set_center(struct ged *gedp,
			      ged_draw_shape_ref ref,
			      const point_t center);

GED_EXPORT extern int
ged_draw_shape_ref_geometry_clear(struct ged *gedp,
				  ged_draw_shape_ref ref);

GED_EXPORT extern int
ged_draw_shape_ref_update_bounds_from_geometry(struct ged *gedp,
	ged_draw_shape_ref ref,
	int *bad_cmd);

GED_EXPORT extern ged_draw_shape_ref
ged_draw_highlighted_shape_ref(const struct ged *gedp);

GED_EXPORT extern ged_draw_shape_ref
ged_draw_highlight_shape_ref_by_name(struct ged *gedp, const char *name);


/**
 * Returns 1 if @p gedp has at least one drawn shape, 0 otherwise.
 */
GED_EXPORT extern int
ged_draw_has_shapes(struct ged *gedp);

/**
 * Report whether @p path is represented in the semantic draw scene.
 *
 * Return values match qged/MGED drawn-state conventions:
 *   - 0: path is not drawn
 *   - 1: path is fully drawn
 *   - 2: path is partially drawn
 *
 * The query is answered from draw records and indexes when available.  When
 * hierarchy coverage must be distinguished from a partial draw and no direct
 * index answer exists, the slow path walks only @p path's database subtree and
 * compares its leaf paths to active drawn shapes.  Pass @p mode < 0 to query
 * all draw modes.
 */
GED_EXPORT extern int
ged_draw_path_state(struct ged *gedp,
		    struct ged_view_context *view_ctx,
		    const char *path,
		    int mode);

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
GED_EXPORT extern int
ged_draw_group_record_in_view(const struct ged_draw_group_record *rec,
			      struct ged_view_context *view_ctx);

/**
 * Append command-visible draw paths to @p result and return the number
 * appended.
 *
 * When @p expanded is zero, paths are draw-intent/group paths suitable for
 * command-level listings such as "who" and render-script replay.  When
 * @p expanded is non-zero, paths are realized leaf shape paths.  Overlay,
 * hidden, empty, and out-of-view records are omitted.  Pass @p mode < 0 to list
 * all draw modes.
 */
GED_EXPORT extern size_t
ged_draw_list_paths(struct ged *gedp,
		    struct ged_view_context *view_ctx,
		    int mode,
		    int expanded,
		    struct bu_vls *result);

/**
 * Return non-zero if command-visible drawn database paths exist for @p v and
 * @p mode.  This uses the same visibility/view/overlay filtering as
 * ged_draw_list_paths(..., expanded=0) without formatting or sorting paths.
 */
GED_EXPORT extern int
ged_draw_has_paths(struct ged *gedp,
		   struct ged_view_context *view_ctx,
		   int mode);

/**
 * Returns the first drawn scene object in display order, or NULL if
 * none are drawn.
 */
GED_EXPORT extern int
ged_draw_shape_count(struct ged *gedp);

/**
 * Erase all scene groups from @p gedp's drawn-object set and free associated
 * draw records and scene objects.
 * This is the "zap" operation.
 */
GED_EXPORT extern void
ged_draw_clear(struct ged *gedp);

/**
 * Erase database draw groups in one draw scope.  Passing NULL for @p v clears
 * shared database draw groups.  Passing an independent view clears that view's
 * independent database draw scope.
 *
 * Returns the number of top-level draw groups removed.
 */
GED_EXPORT extern int
ged_draw_clear_view(struct ged *gedp, struct ged_view_context *view_ctx);

/**
 * Return 1 if any scene groups exist (i.e., something is drawn), 0 if the
 * display is empty.
 */
GED_EXPORT extern int
ged_draw_has_groups(struct ged *gedp);

__END_DECLS

#endif /* GED_SCENE_H */
