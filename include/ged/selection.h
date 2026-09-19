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

#include "ged/db_index.h"
#include "ged/selection_types.h"
#include "ged/scene_types.h"

struct bu_vls;

__BEGIN_DECLS

/** Return non-zero when @p gedp owns an initialized selection service. */
GED_EXPORT extern int
ged_selection_state_available(struct ged *gedp);

/**
 * Begin a nestable semantic-selection mutation batch.
 *
 * Selection mutations made inside the outermost batch are accumulated and
 * presented to attached views and the rendering backend once, when
 * ged_selection_batch_end() is called.
 */
GED_EXPORT extern int
ged_selection_batch_begin(struct ged *gedp);

/** End a selection batch and present its accumulated sparse changes. */
GED_EXPORT extern int
ged_selection_batch_end(struct ged *gedp);

/** Return the number of explicitly selected paths in @p set_name. */
GED_EXPORT extern size_t
ged_selection_count(struct ged *gedp, const char *set_name);

/** Append one non-default selection-set name per line to @p out. */
GED_EXPORT extern int
ged_selection_list_set_names(struct ged *gedp, struct bu_vls *out);

/**
 * Append selected paths from sets matching @p set_pattern to @p out.
 *
 * A NULL or empty pattern identifies the default set.  The return value is
 * the number of matching sets.
 */
GED_EXPORT extern int
ged_selection_list_paths(struct ged *gedp, const char *set_pattern,
	struct bu_vls *out);

/** Return the number of sets matching @p set_pattern. */
GED_EXPORT extern int
ged_selection_set_match_count(struct ged *gedp, const char *set_pattern);

/** Return a deterministic hash of the selected paths in @p set_name. */
GED_EXPORT extern unsigned long long
ged_selection_state_hash(struct ged *gedp, const char *set_name);

/**
 * Copy a renderer-neutral selection identity collection.
 *
 * The return value is the total number of identifiers in @p kind.  If
 * @p hashes is non-NULL, at most @p capacity identifiers are copied.  A
 * caller may therefore query the required capacity with a NULL buffer and
 * then copy the collection without parsing or resolving path strings.
 */
GED_EXPORT extern size_t
ged_selection_hashes(struct ged *gedp, const char *set_name,
	enum ged_selection_hash_kind kind, ged_db_index_id *hashes,
	size_t capacity);

/** Clear one exact semantic selection set. */
GED_EXPORT extern int
ged_selection_clear(struct ged *gedp, const char *set_name);

/** Clear every semantic selection set matching @p set_pattern. */
GED_EXPORT extern int
ged_selection_clear_matching(struct ged *gedp, const char *set_pattern);

/** Select one canonical database path in an exact selection set. */
GED_EXPORT extern int
ged_selection_select_path(struct ged *gedp, const char *set_name,
	const char *path, int recompute_hierarchy);

/** Atomically add multiple canonical database paths to one selection set. */
GED_EXPORT extern int
ged_selection_select_paths(struct ged *gedp, const char *set_name,
	const char *const *paths, size_t path_count, int recompute_hierarchy);

/** Select one path in the single set matching @p set_pattern. */
GED_EXPORT extern int
ged_selection_select_path_matching(struct ged *gedp,
	const char *set_pattern, const char *path, int recompute_hierarchy);

/** Deselect one canonical database path from an exact selection set. */
GED_EXPORT extern int
ged_selection_deselect_path(struct ged *gedp, const char *set_name,
	const char *path, int recompute_hierarchy);

/** Atomically remove multiple paths and their selected descendants. */
GED_EXPORT extern int
ged_selection_deselect_paths(struct ged *gedp, const char *set_name,
	const char *const *paths, size_t path_count, int recompute_hierarchy);

/** Deselect one path from the single set matching @p set_pattern. */
GED_EXPORT extern int
ged_selection_deselect_path_matching(struct ged *gedp,
	const char *set_pattern, const char *path, int recompute_hierarchy);

/** Select one path represented by stable database-index identifiers. */
GED_EXPORT extern int
ged_selection_select_path_ids(struct ged *gedp, const char *set_name,
	const ged_db_index_id *path, size_t path_len, int recompute_hierarchy);

/** Deselect one path represented by stable database-index identifiers. */
GED_EXPORT extern int
ged_selection_deselect_path_ids(struct ged *gedp, const char *set_name,
	const ged_db_index_id *path, size_t path_len, int recompute_hierarchy);

/** Invalidate and lazily recompute hierarchy indexes for one selection set. */
GED_EXPORT extern int
ged_selection_recompute(struct ged *gedp, const char *set_name);

/** Invalidate hierarchy indexes for sets matching @p set_pattern. */
GED_EXPORT extern int
ged_selection_recompute_matching(struct ged *gedp,
	const char *set_pattern);

/** Expand selected combination paths to their terminal descendants. */
GED_EXPORT extern int
ged_selection_expand_matching(struct ged *gedp, const char *set_pattern);

/** Collapse complete sibling selections to their common parent paths. */
GED_EXPORT extern int
ged_selection_collapse_matching(struct ged *gedp, const char *set_pattern);

/** Return non-zero when @p path_hash is explicitly selected. */
GED_EXPORT extern int
ged_selection_is_path_selected(struct ged *gedp, const char *set_name,
	ged_db_index_id path_hash);

/** Return non-zero when @p path_hash is selected or below a selected path. */
GED_EXPORT extern int
ged_selection_is_path_active(struct ged *gedp, const char *set_name,
	ged_db_index_id path_hash);

/** Return non-zero when @p path_hash is an ancestor of a selected path. */
GED_EXPORT extern int
ged_selection_is_path_active_parent(struct ged *gedp,
	const char *set_name, ged_db_index_id path_hash);

/** Return non-zero when @p object_id is any parent of a selected object. */
GED_EXPORT extern int
ged_selection_is_object_parent(struct ged *gedp, const char *set_name,
	ged_db_index_id object_id);

/** Return non-zero when @p object_id immediately parents a selection. */
GED_EXPORT extern int
ged_selection_is_object_immediate_parent(struct ged *gedp,
	const char *set_name, ged_db_index_id object_id);

/** Return non-zero when @p object_id is a non-immediate selection ancestor. */
GED_EXPORT extern int
ged_selection_is_object_grand_parent(struct ged *gedp,
	const char *set_name, ged_db_index_id object_id);

/** Return the number of semantic paths selected in @p view_ctx. */
GED_EXPORT extern size_t
ged_view_selection_count(struct ged_view_context *view_ctx);

/** Visit the selected semantic paths without exposing renderer state. */
GED_EXPORT extern int
ged_view_selection_visit(
    struct ged_view_context *view_ctx,
    ged_view_selection_visit_cb cb,
    void *data);

/** Clear the semantic selection set owned by @p view_ctx. */
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

/** Replace the selection with paths represented by @p result. */
GED_EXPORT extern int
ged_view_selection_set_pick(
    struct ged_view_context *view_ctx,
    const struct ged_pick_result *result,
    ged_view_selection_path_cb cb,
    void *data);

/** Allocate an empty picking result owned by the caller. */
GED_EXPORT extern struct ged_pick_result *
ged_pick_result_create(void);

/** Release a picking result returned by this API. */
GED_EXPORT extern void
ged_pick_result_free(struct ged_pick_result *result);

/** Return the number of ordered hits in @p result. */
GED_EXPORT extern size_t
ged_pick_result_count(const struct ged_pick_result *result);

/** Copy one hit's semantic path into caller-owned @p path_out. */
GED_EXPORT extern int
ged_pick_result_path(const struct ged_pick_result *result,
			  size_t index,
			  struct bu_vls *path_out);

/** Return one hit's model-space ray distance, or a negative value on error. */
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

/** Initialize caller-owned primitive-pick detail storage. */
GED_EXPORT extern void
ged_pick_detail_init(struct ged_pick_detail *detail);

/** Return initialized primitive-pick detail storage. */
GED_HEADER_INLINE struct ged_pick_detail
ged_pick_detail_default(void)
{
    struct ged_pick_detail detail;
    ged_pick_detail_init(&detail);
    return detail;
}

/** Copy the renderer-neutral primitive detail for one hit. */
GED_EXPORT extern int
ged_pick_result_detail(const struct ged_pick_result *result,
	size_t index, struct ged_pick_detail *detail);

/** Append a semantic hit with optional primitive detail. */
GED_EXPORT extern int
ged_pick_result_append_detail(struct ged_pick_result *result,
	const char *path, fastf_t hit_dist,
	const struct ged_pick_detail *detail);

/** Append a semantic path hit without primitive detail. */
GED_EXPORT extern int
ged_pick_result_append_path(struct ged_pick_result *result,
				 const char *path,
				 fastf_t hit_dist);

/** Copy one hit from @p src into @p dest with a replacement distance. */
GED_EXPORT extern int
ged_pick_result_append_copy(struct ged_pick_result *dest,
				 const struct ged_pick_result *src,
				 size_t index,
				 fastf_t hit_dist);

/** Allocate a result containing only the first hit from @p src. */
GED_EXPORT extern struct ged_pick_result *
ged_pick_result_filter_first(const struct ged_pick_result *src);

/** Pick semantic paths at a screen point. */
GED_EXPORT extern struct ged_pick_result *
ged_pick_point(struct ged_view_context *view_ctx,
				 int x,
				 int y,
				 int first_only);

/** Pick the nearest semantic path at a screen point. */
GED_EXPORT extern struct ged_pick_result *
ged_pick_nearest(struct ged_view_context *view_ctx,
				   int x,
				   int y);

/** Pick semantic paths intersecting a screen-space rectangle. */
GED_EXPORT extern struct ged_pick_result *
ged_pick_rect(struct ged_view_context *view_ctx,
				int x0,
				int y0,
				int x1,
				int y1);

/** Resolve a semantic path expression without a graphical renderer pick. */
GED_EXPORT extern struct ged_pick_result *
ged_pick_semantic_path(struct ged *gedp,
					 struct ged_view_context *view_ctx,
					 const char *path_pattern);

/** Find the requested model-space snap candidate nearest @p sample. */
GED_EXPORT extern int
ged_view_selection_snap(struct ged_view_context *view_ctx,
	const point_t sample,
	enum ged_selection_snap_kind kind,
	point_t candidate);


__END_DECLS

#endif /* GED_SELECTION_H */
