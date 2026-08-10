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

#include "ged/defines.h"
#include "ged/scene_types.h"

struct bu_vls;

__BEGIN_DECLS

/** Initialize a database scene style with renderer-neutral defaults. */
GED_EXPORT extern void
ged_scene_style_init(struct ged_scene_style *style);

/** Initialize an automatic realization policy. */
GED_EXPORT extern void
ged_scene_realization_policy_init(
    struct ged_scene_realization_policy *policy);

/** Initialize a typed database draw request. */
GED_EXPORT extern void
ged_scene_draw_request_init(struct ged_scene_draw_request *request);

/** Initialize a typed exact-path erase request. */
GED_EXPORT extern void
ged_scene_erase_request_init(struct ged_scene_erase_request *request);

/** Initialize a typed path-scoped scene request. */
GED_EXPORT extern void
ged_scene_path_request_init(struct ged_scene_path_request *request);

/** Initialize a typed whole-scene redraw request. */
GED_EXPORT extern void
ged_scene_redraw_request_init(struct ged_scene_redraw_request *request);

/** Initialize a typed canonical-scene clear request. */
GED_EXPORT extern void
ged_scene_clear_request_init(struct ged_scene_clear_request *request);

/** Initialize a typed exact-occurrence edit-scope request. */
GED_EXPORT extern void
ged_scene_edit_request_init(struct ged_scene_edit_request *request);

/** Allocate a reusable opaque semantic scene result. */
GED_EXPORT extern struct ged_scene_result *
ged_scene_result_create(void);

/** Reset a semantic scene result while retaining its allocation. */
GED_EXPORT extern void
ged_scene_result_clear(struct ged_scene_result *result);

/** Destroy a semantic scene result created by ged_scene_result_create(). */
GED_EXPORT extern void
ged_scene_result_destroy(struct ged_scene_result *result);

/** Return the operation status stored in a semantic scene result. */
GED_EXPORT extern enum ged_scene_status
ged_scene_result_status(const struct ged_scene_result *result);

/** Return non-zero when the operation changed semantic or presentation state. */
GED_EXPORT extern int
ged_scene_result_changed(const struct ged_scene_result *result);

/** Return the number of affected database paths recorded by the operation. */
GED_EXPORT extern size_t
ged_scene_result_path_count(const struct ged_scene_result *result);

/** Return a callback-lifetime path at @p index, or NULL when out of range. */
GED_EXPORT extern const char *
ged_scene_result_path_at(const struct ged_scene_result *result, size_t index);

/** Return the number of affected semantic draw groups. */
GED_EXPORT extern size_t
ged_scene_result_group_count(const struct ged_scene_result *result);

/** Return the number of affected realized shapes reported by the backend. */
GED_EXPORT extern size_t
ged_scene_result_shape_count(const struct ged_scene_result *result);

/** Return the number of edit-scope conflicts reported by an operation. */
GED_EXPORT extern size_t
ged_scene_result_conflict_count(const struct ged_scene_result *result);

/** Return the scene revision before the operation. */
GED_EXPORT extern uint64_t
ged_scene_result_revision_before(const struct ged_scene_result *result);

/** Return the scene revision after the operation. */
GED_EXPORT extern uint64_t
ged_scene_result_revision_after(const struct ged_scene_result *result);

/** Return a borrowed diagnostic string, or an empty string when none exists. */
GED_EXPORT extern const char *
ged_scene_result_diagnostic(const struct ged_scene_result *result);

/** Apply a typed draw request to the semantic scene. */
GED_EXPORT extern enum ged_scene_status
ged_scene_draw(struct ged *gedp,
	       const struct ged_scene_draw_request *request,
	       struct ged_scene_result *result);

/** Apply a typed erase request to the semantic scene. */
GED_EXPORT extern enum ged_scene_status
ged_scene_erase(struct ged *gedp,
		const struct ged_scene_erase_request *request,
		struct ged_scene_result *result);

/** Rebuild retained sources selected by a typed redraw request. */
GED_EXPORT extern enum ged_scene_status
ged_scene_redraw(struct ged *gedp,
		 const struct ged_scene_redraw_request *request,
		 struct ged_scene_result *result);

/** Set visibility for a path-scoped portion of the semantic scene. */
GED_EXPORT extern enum ged_scene_status
ged_scene_visibility_set(struct ged *gedp,
			 const struct ged_scene_path_request *request,
			 int visible,
			 struct ged_scene_result *result);

/** Set opacity in the inclusive range zero through one for a scene path. */
GED_EXPORT extern enum ged_scene_status
ged_scene_opacity_set(struct ged *gedp,
		      const struct ged_scene_path_request *request,
		      double opacity,
		      struct ged_scene_result *result);

/** Set legacy explicit highlight state for a scene path. */
GED_EXPORT extern enum ged_scene_status
ged_scene_highlight_set(struct ged *gedp,
			const struct ged_scene_path_request *request,
			int highlighted,
			struct ged_scene_result *result);

/** Clear all transient occurrence highlights in the semantic scene. */
GED_EXPORT extern enum ged_scene_status
ged_scene_highlights_clear(struct ged *gedp,
			   struct ged_scene_result *result);

/** Clear canonical or view-scoped database draw intent. */
GED_EXPORT extern enum ged_scene_status
ged_scene_clear(struct ged *gedp,
		const struct ged_scene_clear_request *request,
		struct ged_scene_result *result);

/** Set the default presentation mode used by later draw requests. */
GED_EXPORT extern enum ged_scene_status
ged_scene_default_draw_mode_set(struct ged *gedp,
				enum ged_scene_draw_mode mode,
				struct ged_scene_result *result);

/** Notify the scene that database material colors changed and reapply them. */
GED_EXPORT extern enum ged_scene_status
ged_scene_materials_changed(struct ged *gedp,
			    struct ged_scene_result *result);

/** Return the current database-material presentation revision. */
GED_EXPORT extern uint64_t
ged_scene_material_revision(const struct ged *gedp);

/** Return the current semantic scene revision, or zero if unavailable. */
GED_EXPORT extern uint64_t
ged_scene_revision(const struct ged *gedp);

/** Return non-zero when the GED context can retain or display a scene. */
GED_EXPORT extern int
ged_scene_available(const struct ged *gedp);

/** Return the default draw mode used by requests that do not override it. */
GED_EXPORT extern enum ged_scene_draw_mode
ged_scene_default_draw_mode(const struct ged *gedp);

/** Return non-zero when at least one semantic draw path is present. */
GED_EXPORT extern int
ged_scene_has_paths(struct ged *gedp,
		    struct ged_view_context *view,
		    enum ged_scene_draw_mode mode);

/** Report whether a database path is not drawn, fully drawn, or partial. */
GED_EXPORT extern enum ged_scene_path_state
ged_scene_path_state_get(struct ged *gedp,
			 struct ged_view_context *view,
			 const char *path,
			 enum ged_scene_draw_mode mode);

/**
 * Append newline-delimited semantic paths and return the number appended.
 * GED_SCENE_DRAW_DEFAULT selects all draw modes for this query.
 */
GED_EXPORT extern size_t
ged_scene_paths_append(struct ged *gedp,
		       struct ged_view_context *view,
		       enum ged_scene_draw_mode mode,
		       enum ged_scene_path_listing listing,
		       struct bu_vls *result);

/** Return the number of currently realized semantic occurrences. */
GED_EXPORT extern size_t
ged_scene_occurrence_count(struct ged *gedp);

/** Return non-zero when at least one realized occurrence is present. */
GED_EXPORT extern int
ged_scene_has_occurrences(struct ged *gedp);

/**
 * Compute model-space scene bounds.  Returns zero on success and one when no
 * object contributes to the requested bounds scope.
 */
GED_EXPORT extern int
ged_scene_bounds(struct ged *gedp, vect_t *min, vect_t *max,
		 enum ged_scene_bounds_scope scope);

/** Return the revision of the occurrence-highlight presentation state. */
GED_EXPORT extern uint64_t
ged_scene_highlight_revision(const struct ged *gedp);

/**
 * Acquire an independently addressable scope for interactive editing.
 *
 * The returned reference is owned by @p gedp and must be closed exactly once.
 * Occurrence paths are reported through @p result when it is non-NULL.
 */
GED_EXPORT extern enum ged_scene_status
ged_scene_edit_acquire(struct ged *gedp,
		       const struct ged_scene_edit_request *request,
		       ged_scene_edit_scope_ref *scope,
		       struct ged_scene_result *result);

/** Close a semantic edit scope and report any intervening draw conflicts. */
GED_EXPORT extern enum ged_scene_status
ged_scene_edit_release(struct ged *gedp,
		       ged_scene_edit_scope_ref scope,
		       enum ged_scene_edit_outcome outcome,
		       struct ged_scene_result *result);

/** Return non-zero when an edit-scope reference is null. */
GED_EXPORT extern int
ged_scene_edit_scope_ref_is_null(ged_scene_edit_scope_ref scope);

/** Return the semantic operation represented by a committed delta. */
GED_EXPORT extern enum ged_scene_delta_kind
ged_scene_delta_kind(const struct ged_scene_delta *delta);

/** Return the view scope of a committed delta, or NULL for all views. */
GED_EXPORT extern struct ged_view_context *
ged_scene_delta_view(const struct ged_scene_delta *delta);

/** Return the number of structured affected paths in a committed delta. */
GED_EXPORT extern size_t
ged_scene_delta_path_count(const struct ged_scene_delta *delta);

/** Return a callback-lifetime affected path, or NULL when out of range. */
GED_EXPORT extern const char *
ged_scene_delta_path_at(const struct ged_scene_delta *delta, size_t index);

/** Return non-zero when a delta has no narrower affected-path scope. */
GED_EXPORT extern int
ged_scene_delta_affects_all(const struct ged_scene_delta *delta);

/** Return the number of semantic draw groups changed by a delta. */
GED_EXPORT extern size_t
ged_scene_delta_group_count(const struct ged_scene_delta *delta);

/** Return the number of realized shapes changed by a delta. */
GED_EXPORT extern size_t
ged_scene_delta_shape_count(const struct ged_scene_delta *delta);

/** Return non-zero when a delta changes presentation but not draw intent. */
GED_EXPORT extern int
ged_scene_delta_is_presentation_only(const struct ged_scene_delta *delta);

/** Return the semantic scene revision before a committed delta. */
GED_EXPORT extern uint64_t
ged_scene_delta_revision_before(const struct ged_scene_delta *delta);

/** Return the semantic scene revision after a committed delta. */
GED_EXPORT extern uint64_t
ged_scene_delta_revision_after(const struct ged_scene_delta *delta);

/** Return a callback-lifetime diagnostic string, or an empty string. */
GED_EXPORT extern const char *
ged_scene_delta_diagnostic(const struct ged_scene_delta *delta);

/**
 * Subscribe to committed semantic scene deltas.
 *
 * Observers run synchronously after state reconciliation.  Registration
 * returns zero on failure.  Callbacks must copy any values they retain.
 */
GED_EXPORT extern ged_scene_observer_token
ged_scene_observer_add(struct ged *gedp,
		       ged_scene_observer_func_t callback,
		       void *client_data);

/** Remove a semantic scene observer, returning non-zero if it was live. */
GED_EXPORT extern int
ged_scene_observer_remove(struct ged *gedp, ged_scene_observer_token token);

__END_DECLS

#endif /* GED_SCENE_H */
