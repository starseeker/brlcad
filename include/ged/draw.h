/*                      D R A W . H
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
/** @addtogroup ged_view
 *
 * @brief
 * GED draw-scene API.
 *
 * Modern callers should treat the draw set as semantic records:
 *   - ged_draw_shape_ref / ged_draw_group_ref are stable handles for the
 *     current draw-scene revision.
 *   - ged_draw_shape_record / ged_draw_group_record expose path identity,
 *     draw source/style, visibility, and highlight state without requiring
 *     callers to inspect BSG node layout.
 *
 * The BSG tree remains the private implementation used by the GED/BSG bridge.
 * Installed GED callers should not need BSG node pointers to query or update
 * the draw scene.
 *
 * Drawn database objects are keyed by struct db_full_path.  String path
 * formatting is available for logging/UI, but graph mutation APIs use the
 * structured path form.
 *
 * C++ callers that need to mirror this semantic draw state into an
 * Obol/libbrlobol view controller should use ged/draw_obol.h.
 */
/** @{ */
/* @file ged/draw.h */

#ifndef GED_DRAW_H
#define GED_DRAW_H

#include "common.h"

#include "vmath.h"
#include "bg/polygon_types.h"
#include "bu/color.h"
#include "bu/list.h"
#include "bu/ptbl.h"
#include "bu/vls.h"
#include "ged/defines.h"

/* Forward declaration to keep the include surface narrow.  Callers that
 * actually use the db_full_path variants will already be including
 * "raytrace.h" or "rt/db_fullpath.h" anyway. */
struct db_full_path;
struct db_i;
struct bn_tol;
struct bg_tess_tol;
struct bg_line_layer_builder;
struct directory;
struct rt_view_info;
struct rt_view_lod_policy;

typedef enum ged_draw_stale_reason {
    GED_DRAW_STALE_NONE = 0,
    GED_DRAW_STALE_SOURCE_CHANGED,
    GED_DRAW_STALE_VIEW_INPUT_CHANGED,
    GED_DRAW_STALE_SETTINGS_CHANGED,
    GED_DRAW_STALE_FORCED,
    GED_DRAW_STALE_UPDATE_FAILED
} ged_draw_stale_reason;

typedef enum ged_draw_mode {
    GED_DRAW_MODE_WIRE         = 0,
    GED_DRAW_MODE_SHADED_BOTS  = 1,
    GED_DRAW_MODE_SHADED       = 2,
    GED_DRAW_MODE_EVAL_WIRE    = 3,
    GED_DRAW_MODE_HIDDEN_LINE  = 4,
    GED_DRAW_MODE_EVAL_POINTS  = 5
} ged_draw_mode;

struct ged_draw_appearance_settings {
    int draw_mode;
    int mixed_modes;
    fastf_t transparency;
    int color_override;
    unsigned char color[3];
    int s_line_width;
    fastf_t s_arrow_tip_length;
    fastf_t s_arrow_tip_width;
    int draw_solid_lines_only;
    int draw_non_subtract_only;
    int strict_fallback;
};

#define GED_DRAW_APPEARANCE_SETTINGS_INIT {GED_DRAW_MODE_WIRE, 0, 1.0, 0, {255, 0, 0}, 1, 0.0, 0.0, 0, 0, 0}

typedef enum ged_draw_transaction_kind {
    GED_DRAW_TXN_NONE = 0,
    GED_DRAW_TXN_DRAW,
    GED_DRAW_TXN_ERASE,
    GED_DRAW_TXN_ERASE_PREFIX,
    GED_DRAW_TXN_REDRAW,
    GED_DRAW_TXN_VISIBILITY,
    GED_DRAW_TXN_HIGHLIGHT,
    GED_DRAW_TXN_MATERIAL_CHANGED,
    GED_DRAW_TXN_REFRESH_MATERIAL_COLORS,
    GED_DRAW_TXN_TRANSPARENCY,
    GED_DRAW_TXN_DEFAULT_DRAW_MODE,
    GED_DRAW_TXN_STALE_SOURCE,
    GED_DRAW_TXN_SOURCE_UPDATED,
    GED_DRAW_TXN_SOURCE_RENAMED,
    GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED,
    GED_DRAW_TXN_CLEAR,
    GED_DRAW_TXN_CLEAR_SCOPE,
    GED_DRAW_TXN_TEARDOWN
} ged_draw_transaction_kind;

struct ged_draw_transaction {
    ged_draw_transaction_kind kind;
    const char *path;                 /**< borrowed command-facing path/name */
    const char *new_path;             /**< optional borrowed replacement path/name */
    const char **paths;               /**< optional borrowed draw path array */
    int path_count;                   /**< entries in @p paths */
    const struct db_full_path *dfp;   /**< optional structured path target */
    void *view;                       /**< optional borrowed view context for view-aware transactions */
    int mode;                         /**< optional draw-mode filter; <0 means all modes */
    const void *appearance;           /**< optional borrowed ged_draw_appearance_settings */
    int autoview;                     /**< draw transaction may update view bounds */
    double value;                     /**< visibility/highlight/transparency/mode */
    ged_draw_stale_reason stale_reason;
    int removed;
    int redraw;
};

struct ged_draw_transaction_result {
    int status;                       /**< >=0 success/count, <0 error */
    int affected_groups;
    int affected_shapes;
    int stale_count;
    int redrawn_count;
    int error_count;
    uint64_t scene_revision_before;
    uint64_t scene_revision_after;
    struct bu_vls names;              /**< affected user-facing target names */
    struct bu_vls errors;             /**< diagnostic details */
};

typedef uintptr_t ged_draw_observer_token;
typedef void (*ged_draw_transaction_observer_func_t)(
	struct ged *gedp,
	const struct ged_draw_transaction *txn,
	const struct ged_draw_transaction_result *result,
	void *client_data);

__BEGIN_DECLS

typedef struct ged_draw_shape_ref {
    uintptr_t token;
    uint64_t scene_revision;
} ged_draw_shape_ref;

typedef struct ged_draw_group_ref {
    uintptr_t token;
    uint64_t scene_revision;
} ged_draw_group_ref;

struct ged_draw_view_line_style {
    int color[3];
    int line_width;
};

struct ged_draw_view_feature_style {
    int visible;
    int color_valid;
    unsigned char color[3];
    int line_width;
    int line_style;
    int arrow;
    fastf_t arrow_tip_length;
    fastf_t arrow_tip_width;
};

#define GED_DRAW_VIEW_FEATURE_STYLE_INIT { -1, 0, {0, 0, 0}, -1, -1, -1, -1.0, -1.0 }

struct ged_draw_view_label_data {
    const char *text;
    point_t point;
    int color_valid;
    unsigned char color[3];
    int line_flag;
    point_t target;
    int anchor;
    int arrow;
};

#define GED_DRAW_VIEW_LABEL_DATA_INIT { NULL, VINIT_ZERO, 0, {0, 0, 0}, 0, VINIT_ZERO, 0, 0 }

struct ged_draw_view_axes_state {
    point_t position;
    fastf_t size;
    int line_width;
    int color[3];
};

struct ged_draw_view_line_layer_data {
    const char *name;
    const point_t *points;
    const int *commands;
    size_t point_count;
    struct ged_draw_view_feature_style style;
};

#define GED_DRAW_VIEW_LINE_LAYER_DATA_INIT { NULL, NULL, NULL, 0, GED_DRAW_VIEW_FEATURE_STYLE_INIT }

enum ged_draw_view_snap_kind {
    GED_DRAW_VIEW_SNAP_GRID = 1,
    GED_DRAW_VIEW_SNAP_ENDPOINT = 2
};

typedef struct rt_view_lod_policy ged_draw_view_lod_policy;

typedef int (*ged_draw_view_feature_depth_cb)(fastf_t depth, void *data);
typedef int (*ged_draw_view_context_selection_path_cb)(void *view_ctx,
						       const char *path,
						       void *data);

enum ged_draw_view_polygon_type {
    GED_DRAW_VIEW_POLYGON_GENERAL = 0,
    GED_DRAW_VIEW_POLYGON_CIRCLE = 1,
    GED_DRAW_VIEW_POLYGON_ELLIPSE = 2,
    GED_DRAW_VIEW_POLYGON_RECTANGLE = 3,
    GED_DRAW_VIEW_POLYGON_SQUARE = 4
};

enum ged_draw_view_polygon_update {
    GED_DRAW_VIEW_POLYGON_UPDATE_DEFAULT = 0,
    GED_DRAW_VIEW_POLYGON_UPDATE_PROPS_ONLY = 1,
    GED_DRAW_VIEW_POLYGON_UPDATE_PT_SELECT = 2,
    GED_DRAW_VIEW_POLYGON_UPDATE_PT_SELECT_CLEAR = 3,
    GED_DRAW_VIEW_POLYGON_UPDATE_PT_MOVE = 4,
    GED_DRAW_VIEW_POLYGON_UPDATE_PT_APPEND = 5
};

typedef struct ged_draw_view_polygon_ref {
    uintptr_t token;
    uint64_t revision;
} ged_draw_view_polygon_ref;

#define GED_DRAW_VIEW_POLYGON_REF_NULL_INIT {0, 0}

#ifdef __cplusplus
#  define GED_DRAW_VIEW_POLYGON_REF_NULL ged_draw_view_polygon_ref{0, 0}
#else
#  define GED_DRAW_VIEW_POLYGON_REF_NULL ((ged_draw_view_polygon_ref){0, 0})
#endif

struct ged_draw_view_polygon_record {
    ged_draw_view_polygon_ref ref;
    const char *name;
    int type;
    int fill_flag;
    vect2d_t fill_dir;
    fastf_t fill_delta;
    struct bu_color fill_color;
    unsigned char edge_color[3];
    long curr_contour_i;
    long curr_point_i;
    int first_contour_open;
    size_t contour_count;
    size_t point_count;
    point_t origin_point;
    plane_t vp;
    fastf_t vZ;
    void *user_data;
};

#define GED_DRAW_SHAPE_REF_NULL_INIT {0, 0}
#define GED_DRAW_GROUP_REF_NULL_INIT {0, 0}

#ifdef __cplusplus
#  define GED_DRAW_SHAPE_REF_NULL ged_draw_shape_ref{0, 0}
#  define GED_DRAW_GROUP_REF_NULL ged_draw_group_ref{0, 0}
#else
#  define GED_DRAW_SHAPE_REF_NULL ((ged_draw_shape_ref){0, 0})
#  define GED_DRAW_GROUP_REF_NULL ((ged_draw_group_ref){0, 0})
#endif

struct ged_draw_shape_record {
    ged_draw_shape_ref ref;
    ged_draw_group_ref group;
    const struct db_full_path *fullpath; /**< borrowed; valid while draw scene is unchanged */
    const char *display_name;            /**< borrowed */
    const char *leaf_name;               /**< borrowed */
    unsigned long long path_hash;
    int visible;
    int highlighted;
    int selected;
    int evaluated_region;
    int stale;
    const char *stale_reason;
    uint64_t source_revision;
    uint64_t inputs_revision;
    uint64_t realized_source_revision;
    uint64_t realized_inputs_revision;
    uint64_t drawn_revision;
    fastf_t transparency;
    int draw_mode;
    int line_width;
    point_t center;
};

struct ged_draw_shape_source_snapshot {
    struct db_i *dbip;
    const struct db_full_path *fullpath; /**< borrowed; valid while draw scene is unchanged */
    struct directory *leaf_dp;           /**< borrowed leaf directory, if known */
    const char *name;                    /**< borrowed fallback display name */
    const struct bn_tol *tol;            /**< borrowed draw source tolerance */
    const struct bg_tess_tol *ttol;      /**< borrowed draw source tessellation tolerance */
};

struct ged_draw_group_record {
    ged_draw_group_ref ref;
    const char *path;                    /**< borrowed */
    const struct db_full_path *fullpath; /**< borrowed; valid while draw scene is unchanged */
    void *view;                          /**< borrowed view context; NULL for view-independent groups */
    int draw_mode;
    fastf_t transparency;
    int visible;
    int is_overlay;
    int is_view_scope;
    int in_view_scope;
    int is_view_source;
    int is_local_source;
    int shape_count;
};

struct ged_draw_view_db_object_record {
    const char *path;                    /**< borrowed; valid during callback */
    const char *type_name;               /**< borrowed static string */
    const char *geometry_name;           /**< borrowed static string */
    int draw_mode;
    fastf_t transparency;
    int evaluated_region;
    int is_database_source;
    int non_database_source;
    int is_database_intent;
    int is_local_source;
    int is_view_source;
    int highlighted;
    int selected;
    int visible;
    int line_style;
    unsigned char color[3];
    mat_t model_mat;
    point_t bounds_center;
    fastf_t bounds_radius;
    int has_bounds;
    size_t vlist_structure_count;
    size_t vlist_point_count;
    uint64_t cache_identity;
    uint64_t source_identity;
    uintptr_t detail_token;              /**< private; valid during callback */
};

struct ged_draw_view_annotation_summary {
    int valid;
    int display_space;
    size_t point_count;
    size_t segment_count;
    int has_point;
    point_t point;
    int line_segment_valid;
    int line_start;
    int line_end;
    int text_segment_valid;
    int text_ref_point;
    const char *text;                    /**< borrowed; valid during callback */
    uint64_t cache_identity;
    uint64_t source_identity;
};

enum ged_draw_view_line_command {
    GED_DRAW_VIEW_LINE_MOVE = 0,
    GED_DRAW_VIEW_LINE_DRAW = 1
};

struct ged_draw_view_line_summary {
    int valid;
    size_t point_count;
    uint64_t cache_identity;
    uint64_t source_identity;
};

struct ged_draw_shape_geometry_summary {
    int valid;
    const char *geometry_name;           /**< borrowed static string */
    size_t point_count;
    size_t index_count;
};

struct ged_draw_scene_display_summary {
    int valid;
    int is_database_source;
    int has_draw_intent;
    const char *intent_path;             /**< borrowed; valid while scene is unchanged */
    int intent_draw_mode;
    int visible;
    int highlighted;
    int line_style;
    int line_width;
    fastf_t transparency;
    int draw_mode;
    int material_valid;
    unsigned char material_color[3];
};

struct ged_draw_database_source_summary {
    int valid;
    int is_database_source;
    int has_state;
    int stale;
    const char *database_path;          /**< borrowed; valid while scene is unchanged */
    const char *stale_reason;           /**< borrowed static string */
    uint64_t source_revision;
    uint64_t inputs_revision;
    uint64_t realized_source_revision;
    uint64_t realized_inputs_revision;
    uint64_t realization_identity;
};

struct ged_draw_shape_material_summary {
    int valid;
    uint64_t material_revision;
    unsigned char material_color[3];
};

struct ged_draw_scene_tree_summary {
    int valid;
    int is_group;
    int is_shape;
    int has_parent;
    const char *name;                    /**< borrowed; valid while scene is unchanged */
    const struct db_full_path *fullpath; /**< borrowed; valid while scene is unchanged */
    int draw_tree_depth;
    size_t child_count;
};

struct ged_draw_view_surface_summary {
    int valid;
    size_t point_count;
    size_t normal_count;
    size_t index_count;
    size_t face_count;
    int normals_per_index;
    int material_valid;
    int material_draw_mode;
    fastf_t material_transparency;
    int material_highlighted;
    unsigned char material_color[3];
    uint64_t cache_identity;
    uint64_t source_identity;
};

struct ged_draw_view_rendered_object_summary {
    int found;
    int is_indexed_face_set;
    int is_annotation;
    size_t surface_point_count;
    size_t surface_normal_count;
    size_t surface_index_count;
    size_t surface_face_count;
    int surface_normals_per_index;
    int surface_material_valid;
    size_t annotation_segment_count;
    uint64_t source_identity;
    int material_draw_mode;
    fastf_t material_transparency;
    int selection_visible;
    int selection_highlighted;
};

struct ged_draw_view_feature_summary {
    int exists;
    int is_overlay;
    int is_label;
    int is_transient_preview;
    int visible;
    unsigned char color[3];
    size_t child_count;
    size_t geometry_command_count;
};

#define GED_DRAW_VIEW_FEATURE_SUMMARY_INIT { 0, 0, 0, 0, 0, {0, 0, 0}, 0, 0 }

typedef int (*ged_draw_view_db_object_record_cb)(
	const struct ged_draw_view_db_object_record *rec,
	void *userdata);

typedef int (*ged_draw_view_segment_cb)(const point_t a,
					const point_t b,
					void *userdata);

typedef int (*ged_draw_view_point_cb)(const point_t pt,
				      void *userdata);

typedef int (*ged_draw_shape_ref_index_cb)(ged_draw_shape_ref ref,
					   void *userdata);

#define GED_DRAW_VIEW_RECORD_QUERY_VIEW_OBJECTS 0x01u
#define GED_DRAW_VIEW_RECORD_QUERY_DB_OBJECTS   0x02u
#define GED_DRAW_VIEW_RECORD_QUERY_LOCAL_ONLY   0x04u

struct ged_draw_view_record_query {
    unsigned int flags;
    const char *glob;
    int draw_mode;                       /**< negative means any draw mode */
};

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
ged_draw_mark_database_change(struct ged *gedp,
			      const char *path,
			      ged_draw_stale_reason reason);

GED_EXPORT extern int
ged_draw_apply_database_update(struct ged *gedp,
			       const char *path,
			       int removed,
			       int redraw);

GED_EXPORT extern const char *
ged_draw_stale_reason_name(ged_draw_stale_reason reason);

/**
 * Return non-zero when @p gedp has an initialized draw scene.
 */
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
ged_draw_shape_record_get(struct ged *gedp,
			  ged_draw_shape_ref ref,
			  struct ged_draw_shape_record *out);

GED_EXPORT extern int
ged_draw_shape_ref_source_snapshot(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	struct ged_draw_shape_source_snapshot *out);

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
	void *view_ctx,
	const struct ged_draw_view_record_query *query,
	ged_draw_view_db_object_record_cb cb,
	void *userdata);

GED_EXPORT extern void
ged_draw_foreach_view_db_object_record(void *view_ctx,
				       ged_draw_view_db_object_record_cb cb,
				       void *userdata);

GED_EXPORT extern void
ged_draw_foreach_visible_view_db_object_record(void *view_ctx,
					      ged_draw_view_db_object_record_cb cb,
					      void *userdata);

GED_EXPORT extern void
ged_draw_foreach_visible_view_db_object_record_mode(
	void *view_ctx,
	int draw_mode,
	ged_draw_view_db_object_record_cb cb,
	void *userdata);

GED_EXPORT extern void
ged_draw_foreach_visible_view_record(void *view_ctx,
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
	void *view_ctx,
	uint64_t cache_identity,
	struct ged_draw_view_rendered_object_summary *out);

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

GED_EXPORT extern ged_draw_shape_ref
ged_draw_advance_shape_ref(struct ged *gedp, ged_draw_shape_ref ref, int delta);

GED_EXPORT extern ged_draw_group_ref
ged_draw_group_ref_of_shape(struct ged *gedp, ged_draw_shape_ref ref);

GED_EXPORT extern void *
ged_draw_first_shape_context(struct ged *gedp);

GED_EXPORT extern void *
ged_draw_shape_ref_context(struct ged *gedp, ged_draw_shape_ref ref);

GED_EXPORT extern void *
ged_draw_shape_ref_cache_context(struct ged *gedp, ged_draw_shape_ref ref);

GED_EXPORT extern void *
ged_draw_shape_ref_view_context(struct ged *gedp, ged_draw_shape_ref ref);

GED_EXPORT extern void *
ged_draw_group_ref_context(struct ged *gedp, ged_draw_group_ref ref);

GED_EXPORT extern ged_draw_shape_ref
ged_draw_shape_ref_from_context(struct ged *gedp, void *shape_ctx);

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
ged_draw_shape_context_line_summary(void *shape_ctx,
				    struct ged_draw_view_line_summary *out);

GED_EXPORT extern int
ged_draw_shape_context_line_point_at(void *shape_ctx,
				     size_t index,
				     point_t out);

GED_EXPORT extern int
ged_draw_shape_context_line_command_at(void *shape_ctx,
				       size_t index,
				       int *out);

GED_EXPORT extern int
ged_draw_shape_context_geometry_summary(void *shape_ctx,
					struct ged_draw_shape_geometry_summary *out);

GED_EXPORT extern int
ged_draw_shape_context_has_state(void *shape_ctx);

GED_EXPORT extern void *
ged_draw_shape_context_source(void *shape_ctx);

GED_EXPORT extern void *
ged_draw_view_context_scene_root(void *view_ctx);

GED_EXPORT extern int
ged_draw_view_context_scene_attached(void *view_ctx);

GED_EXPORT extern int
ged_draw_view_context_hud_sync(void *view_ctx);

GED_EXPORT extern size_t
ged_draw_view_context_selection_count(void *view_ctx);

GED_EXPORT extern int
ged_draw_view_context_selection_path_foreach(
	void *view_ctx,
	ged_draw_view_context_selection_path_cb cb,
	void *data);

GED_EXPORT extern int
ged_draw_view_context_selection_clear(void *view_ctx);

GED_EXPORT extern int
ged_draw_view_selection_add_shape_ref_context(
	struct ged *gedp,
	void *view_ctx,
	ged_draw_shape_ref ref,
	void **selection_view_ctx,
	struct bu_vls *path);

GED_EXPORT extern uint64_t
ged_draw_view_context_frame_revision(void *view_ctx);

GED_EXPORT extern uint64_t
ged_draw_view_context_bump_frame_revision(void *view_ctx);

GED_EXPORT extern int
ged_draw_view_context_feature_remove(void *view_ctx,
				     const char *name);

GED_EXPORT extern int
ged_draw_view_context_features_remove_prefix(void *view_ctx,
					     const char *prefix);

GED_EXPORT extern int
ged_draw_view_context_feature_exists(void *view_ctx,
				     const char *name);

GED_EXPORT extern int
ged_draw_view_context_feature_visible(void *view_ctx,
				      const char *name);

GED_EXPORT extern int
ged_draw_view_context_feature_visible_set(void *view_ctx,
					  const char *name,
					  int visible);

GED_EXPORT extern int
ged_draw_view_context_feature_style_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_feature_style *style);

GED_EXPORT extern int
ged_draw_view_context_feature_style_apply(
	void *view_ctx,
	const char *name,
	const struct ged_draw_view_feature_style *style,
	int recursive);

GED_EXPORT extern int
ged_draw_view_context_feature_realize(void *view_ctx,
				      const char *name,
				      int recursive);

GED_EXPORT extern int
ged_draw_view_context_gobject_create(
	struct ged *gedp,
	void *view_ctx,
	const char *db_path,
	const char *gobject_name,
	struct bu_vls *result);

GED_EXPORT extern int
ged_draw_view_context_feature_depth(void *view_ctx,
				    const char *name,
				    int mode,
				    fastf_t *depth);

GED_EXPORT extern int
ged_draw_view_context_feature_depth_foreach(
	void *view_ctx,
	int mode,
	ged_draw_view_feature_depth_cb cb,
	void *data);

GED_EXPORT extern int
ged_draw_view_context_snap_first_candidate(void *view_ctx,
					   const point_t sample,
					   enum ged_draw_view_snap_kind kind,
					   point_t candidate);

GED_EXPORT extern void
ged_draw_view_context_info_get(struct rt_view_info *view_info,
			       const void *view_ctx);

GED_EXPORT extern int
ged_draw_view_context_lod_policy_get(
	ged_draw_view_lod_policy *policy,
	const void *view_ctx);

GED_EXPORT extern int
ged_draw_view_context_lod_policy_apply(
	void *view_ctx,
	const ged_draw_view_lod_policy *policy);

GED_EXPORT extern int
ged_draw_view_context_lod_policy_apply_bot_threshold(
	void *view_ctx,
	const ged_draw_view_lod_policy *policy,
	size_t bot_threshold);

GED_EXPORT extern int
ged_draw_view_context_indexed_face_set_replace(
	void *view_ctx,
	const char *name,
	int local,
	const point_t *points,
	size_t point_count,
	const vect_t *normals,
	size_t normal_count,
	const int *indices,
	size_t index_count,
	const struct ged_draw_view_feature_style *style);

GED_EXPORT extern int
ged_draw_view_context_lines_replace(
	void *view_ctx,
	const char *name,
	int local,
	const point_t *points,
	const int *cmds,
	size_t point_count,
	const struct ged_draw_view_feature_style *style);

GED_EXPORT extern int
ged_draw_view_context_lines_create_model_annotation(
	void *view_ctx,
	const char *name,
	int local,
	const point_t point);

GED_EXPORT extern int
ged_draw_view_context_lines_append_point(void *view_ctx,
					 const char *name,
					 const point_t point);

GED_EXPORT extern int
ged_draw_view_context_line_layer_builder_replace(
	void *view_ctx,
	const char *name,
	int local,
	const struct bg_line_layer_builder *builder);

GED_EXPORT extern int
ged_draw_view_context_line_layers_replace(
	void *view_ctx,
	const char *name,
	int local,
	const struct ged_draw_view_line_layer_data *layers,
	size_t layer_count,
	const struct ged_draw_view_feature_style *style);

GED_EXPORT extern int
ged_draw_view_context_tcl_polygons_replace(
	void *view_ctx,
	const char *name,
	const point_t *points,
	const int *cmds,
	size_t point_count,
	const struct ged_draw_view_feature_style *style);

GED_EXPORT extern int
ged_draw_view_context_tcl_labels_replace(
	void *view_ctx,
	const char *name,
	int draw,
	const struct ged_draw_view_label_data *labels,
	size_t label_count);

GED_EXPORT extern int
ged_draw_view_context_labels_replace(
	void *view_ctx,
	const char *name,
	int local,
	const struct ged_draw_view_label_data *labels,
	size_t label_count);

GED_EXPORT extern int
ged_draw_view_context_label_create(void *view_ctx,
				   const char *name,
				   int local,
				   const char *text,
				   const point_t point,
				   const point_t target,
				   int has_target);

GED_EXPORT extern size_t
ged_draw_view_context_label_count(void *view_ctx,
				  const char *name);

GED_EXPORT extern int
ged_draw_view_context_label_copy(void *view_ctx,
				 const char *name,
				 size_t index,
				 struct bu_vls *text,
				 point_t point,
				 unsigned char rgb[3]);

GED_EXPORT extern int
ged_draw_view_context_label_point_set(void *view_ctx,
				      const char *name,
				      size_t index,
				      const point_t point);

GED_EXPORT extern int
ged_draw_view_context_line_color_set(void *view_ctx,
				     const char *name,
				     int r,
				     int g,
				     int b);

GED_EXPORT extern int
ged_draw_view_context_line_width_set(void *view_ctx,
				     const char *name,
				     int line_width);

GED_EXPORT extern int
ged_draw_view_context_line_style_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_line_style *style);

GED_EXPORT extern int
ged_draw_view_context_feature_points_copy(void *view_ctx,
					  const char *name,
					  point_t **points,
					  size_t *point_count);

GED_EXPORT extern int
ged_draw_view_context_feature_line_command_at(
	void *view_ctx,
	const char *name,
	size_t index,
	int *out);

GED_EXPORT extern int
ged_draw_view_context_lines_points_copy(
	void *view_ctx,
	const char *name,
	point_t **points,
	size_t *point_count);

GED_EXPORT extern int
ged_draw_view_context_tcl_lines_replace(
	void *view_ctx,
	const char *name,
	const point_t *points,
	size_t point_count,
	const struct ged_draw_view_line_style *style);

GED_EXPORT extern int
ged_draw_view_context_arrow_tip_get(void *view_ctx,
				    const char *name,
				    fastf_t *tip_length,
				    fastf_t *tip_width);

GED_EXPORT extern int
ged_draw_view_context_arrow_tip_set(void *view_ctx,
				    const char *name,
				    fastf_t tip_length,
				    fastf_t tip_width);

GED_EXPORT extern int
ged_draw_view_context_tcl_arrows_replace(
	void *view_ctx,
	const char *name,
	const point_t *points,
	size_t point_count,
	const struct ged_draw_view_feature_style *style);

GED_EXPORT extern int
ged_draw_view_context_feature_axes_centers_copy(
	void *view_ctx,
	const char *name,
	point_t **centers,
	size_t *center_count);

GED_EXPORT extern int
ged_draw_view_context_tcl_axes_replace(
	void *view_ctx,
	const char *name,
	const point_t *centers,
	size_t center_count,
	fastf_t half_axes_size,
	const struct ged_draw_view_feature_style *style);

GED_EXPORT extern int
ged_draw_view_context_axes_create(
	void *view_ctx,
	const char *name,
	int local,
	const struct ged_draw_view_axes_state *state);

GED_EXPORT extern int
ged_draw_view_context_axes_state_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_axes_state *state);

GED_EXPORT extern int
ged_draw_view_context_axes_state_replace(
	void *view_ctx,
	const char *name,
	const struct ged_draw_view_axes_state *state);

GED_EXPORT extern int
ged_draw_view_polygon_ref_is_null(ged_draw_view_polygon_ref ref);

GED_EXPORT extern ged_draw_view_polygon_ref
ged_draw_view_context_polygon_find(void *view_ctx,
				   const char *name);

GED_EXPORT extern ged_draw_view_polygon_ref
ged_draw_view_context_polygon_find_scoped(void *view_ctx,
					  const char *name,
					  int local_only);

GED_EXPORT extern ged_draw_view_polygon_ref
ged_draw_view_context_polygon_create(void *view_ctx,
				     const char *name,
				     int local,
				     int type,
				     const point_t screen_point);

GED_EXPORT extern ged_draw_view_polygon_ref
ged_draw_view_context_polygon_import_sketch(const char *name,
					    struct db_i *dbip,
					    struct directory *dp,
					    void *view_ctx,
					    int local);

GED_EXPORT extern int
ged_draw_view_polygon_export_sketch(struct db_i *dbip,
				    const char *name,
				    ged_draw_view_polygon_ref ref);

GED_EXPORT extern int
ged_draw_view_polygon_record_get(ged_draw_view_polygon_ref ref,
				 struct ged_draw_view_polygon_record *record);

GED_EXPORT extern int
ged_draw_view_polygon_has_data(ged_draw_view_polygon_ref ref);

GED_EXPORT extern int
ged_draw_view_context_polygon_update(ged_draw_view_polygon_ref ref,
				     void *view_ctx,
				     int op);

GED_EXPORT extern int
ged_draw_view_context_polygon_update_screen_pt(ged_draw_view_polygon_ref ref,
					       void *view_ctx,
					       int x,
					       int y,
					       int op);

GED_EXPORT extern int
ged_draw_view_polygon_set_current(ged_draw_view_polygon_ref ref,
				  long contour_i,
				  long point_i);

GED_EXPORT extern int
ged_draw_view_polygon_set_contour_open(ged_draw_view_polygon_ref ref,
				       long contour_i,
				       int open);

GED_EXPORT extern int
ged_draw_view_polygon_set_all_contours_open(ged_draw_view_polygon_ref ref,
					    int open);

GED_EXPORT extern int
ged_draw_view_context_polygon_area(ged_draw_view_polygon_ref ref,
				   void *view_ctx,
				   fastf_t *area);

GED_EXPORT extern int
ged_draw_view_context_polygon_overlap(ged_draw_view_polygon_ref ref,
				      void *view_ctx,
				      const char *other_name,
				      const struct bn_tol *tol,
				      int *overlap);

GED_EXPORT extern int
ged_draw_view_polygon_set_fill(ged_draw_view_polygon_ref ref,
			       int fill_flag,
			       fastf_t fill_slope_x,
			       fastf_t fill_slope_y,
			       fastf_t fill_density);

GED_EXPORT extern int
ged_draw_view_polygon_fill_color_get(ged_draw_view_polygon_ref ref,
				     struct bu_color *fill_color);

GED_EXPORT extern int
ged_draw_view_polygon_fill_color_set(ged_draw_view_polygon_ref ref,
				     const struct bu_color *fill_color);

GED_EXPORT extern int
ged_draw_view_context_polygon_csg(ged_draw_view_polygon_ref target,
				  void *view_ctx,
				  const char *other_name,
				  bg_clip_t op);

GED_EXPORT extern int
ged_draw_view_context_feature_summary(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_feature_summary *summary);

GED_EXPORT extern const struct db_full_path *
ged_draw_scene_context_fullpath(void *scene_ctx);

GED_EXPORT extern int
ged_draw_group_context_dbpath(struct ged *gedp,
			      void *group_ctx,
			      struct db_full_path *out);

GED_EXPORT extern int
ged_draw_group_context_is_overlay(void *group_ctx);

GED_EXPORT extern int
ged_draw_scene_context_display_summary(void *scene_ctx,
				       struct ged_draw_scene_display_summary *out);

GED_EXPORT extern int
ged_draw_scene_context_source_summary(void *scene_ctx,
				      struct ged_draw_database_source_summary *out);

GED_EXPORT extern int
ged_draw_scene_context_tree_summary(void *scene_ctx,
				    struct ged_draw_scene_tree_summary *out);

GED_EXPORT extern void *
ged_draw_scene_context_child_at(void *scene_ctx, size_t index);

GED_EXPORT extern void *
ged_draw_scene_context_parent(void *scene_ctx);

GED_EXPORT extern const char *
ged_draw_scene_context_name(void *scene_ctx);

GED_EXPORT extern int
ged_draw_scene_context_subtree_bounds(void *scene_ctx,
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

GED_EXPORT extern int
ged_draw_shape_ref_publish_line_set(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    const point_t *points,
				    const int *commands,
				    size_t point_count);

GED_EXPORT extern int
ged_draw_shape_ref_publish_point_set(struct ged *gedp,
				     ged_draw_shape_ref ref,
				     const point_t *points,
				     size_t point_count);

GED_EXPORT extern int
ged_draw_shape_ref_publish_indexed_face_set(struct ged *gedp,
					    ged_draw_shape_ref ref,
					    const point_t *points,
					    size_t point_count,
					    const vect_t *normals,
					    size_t normal_count,
					    const int *indices,
					    size_t index_count);

GED_EXPORT extern int
ged_draw_shape_ref_publish_primitive_wireframe(struct ged *gedp,
					       ged_draw_shape_ref ref,
					       struct rt_db_internal *ip,
					       const struct bg_tess_tol *ttol,
					       const struct bn_tol *tol,
					       void *view_ctx,
					       int adaptive);

GED_EXPORT extern int
ged_draw_shape_ref_eval_wireframe(struct ged *gedp,
				  ged_draw_shape_ref ref);

GED_EXPORT extern int
ged_draw_shape_ref_eval_points(struct ged *gedp,
			       ged_draw_shape_ref ref);

GED_EXPORT extern ged_draw_shape_ref
ged_draw_highlighted_shape_ref(const struct ged *gedp);

GED_EXPORT extern ged_draw_shape_ref
ged_draw_highlight_shape_ref_by_name(struct ged *gedp, const char *name);

GED_EXPORT extern int
ged_draw_view_selection_set_highlighted_shape_ref(struct ged *gedp,
						  void *view_ctx,
						  ged_draw_shape_ref ref);

/**
 * Returns 1 if @p gedp has at least one drawn shape, 0 otherwise.
 */
GED_EXPORT extern int
ged_draw_has_shapes(struct ged *gedp);

/**
 * Report whether @p path is represented in the retained draw scene.
 *
 * Return values match qged/MGED drawn-state conventions:
 *   - 0: path is not drawn
 *   - 1: path is fully drawn
 *   - 2: path is partially drawn
 *
 * The query is answered from retained draw records.  When hierarchy coverage
 * must be distinguished from a partial draw, the no-index fallback walks only
 * @p path's database subtree and compares its leaf paths to active retained
 * shapes.  Pass @p mode < 0 to query all draw modes.
 */
GED_EXPORT extern int
ged_draw_path_state(struct ged *gedp,
		    void *view_ctx,
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
			      void *view_ctx);

/**
 * Append retained draw paths to @p result and return the number appended.
 *
 * When @p expanded is zero, paths are draw-intent/group paths suitable for
 * command-level listings such as "who" and render-script replay.  When
 * @p expanded is non-zero, paths are realized leaf shape paths.  Overlay,
 * hidden, empty, and out-of-view records are omitted.  Pass @p mode < 0 to list
 * all draw modes.
 */
GED_EXPORT extern size_t
ged_draw_list_paths(struct ged *gedp,
		    void *view_ctx,
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
		   void *view_ctx,
		   int mode);

/**
 * Returns the first drawn scene object in display order, or NULL if
 * none are drawn.
 */
GED_EXPORT extern int
ged_draw_shape_count(struct ged *gedp);

/**
 * Erase all scene groups from @p gedp's drawn-object set, destroying vlists
 * and freeing all associated scene objects.
 * This is the "zap" operation.
 */
GED_EXPORT extern void
ged_draw_clear(struct ged *gedp);

/**
 * Erase database draw groups in one draw scope.  Passing NULL for @p v clears
 * shared database draw groups.  Passing an independent view clears that view's
 * independent database draw scope.
 *
 * Returns the number of top-level retained draw groups removed.
 */
GED_EXPORT extern int
ged_draw_clear_view(struct ged *gedp, void *view_ctx);

/**
 * Return 1 if any scene groups exist (i.e., something is drawn), 0 if the
 * display is empty.
 */
GED_EXPORT extern int
ged_draw_has_groups(struct ged *gedp);

__END_DECLS

#endif /* GED_DRAW_H */

/** @} */
/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
