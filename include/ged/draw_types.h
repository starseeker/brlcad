/*                      D R A W _ T Y P E S . H
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
 *     callers to inspect the underlying scene-graph implementation.
 *
 * libged owns the semantic draw records and mirrors them into Obol/libBObol
 * scene controllers.  Installed GED callers should use this API rather than
 * direct scene-graph node pointers to query or update the draw scene.
 *
 * Drawn database objects are keyed by struct db_full_path.  String path
 * formatting is available for logging/UI, but graph mutation APIs use the
 * structured path form.
 *
 * C++ callers that need to mirror this semantic draw state into an
 * Obol/libBObol view controller should use ged/display.h.
 */
/** @{ */
/* @file ged/draw_types.h */

#ifndef GED_DRAW_TYPES_H
#define GED_DRAW_TYPES_H

#include "common.h"

#include <stddef.h>

#include "vmath.h"
#include "bg/polygon_types.h"
#include "bu/color.h"
#include "bu/list.h"
#include "bu/ptbl.h"
#include "bu/vls.h"
#include "bv/view.h"
#include "ged/defines.h"
#include "ged/draw_scene.h"

/* Forward declaration to keep the include surface narrow.  Callers that
 * actually use the db_full_path variants will already be including
 * "raytrace.h" or "rt/db_fullpath.h" anyway. */
struct db_full_path;
struct db_i;
struct bn_tol;
struct bg_tess_tol;
struct bg_line_layer_builder;
struct directory;
struct ged;
struct ged_pick_result;
struct bv_view_info;
struct rt_db_internal;
struct bg_tess_tol;
struct bn_tol;

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
    /* Legacy draw-option opacity: 0 is clear and 1 is opaque. */
    fastf_t transparency;
    int color_override;
    unsigned char color[3];
    int s_line_width;
    fastf_t s_arrow_tip_length;
    fastf_t s_arrow_tip_width;
    int draw_solid_lines_only;
    int draw_non_subtract_only;
    int strict_fallback;
    int defer_leaf_expansion;        /**< Obol database-source draw may publish the requested root first. */
};

#define GED_DRAW_APPEARANCE_SETTINGS_INIT {GED_DRAW_MODE_WIRE, 0, 1.0, 0, {255, 0, 0}, 1, 0.0, 0.0, 0, 0, 0, 0}

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
    struct ged_view_context *view;    /**< optional borrowed view for view-aware transactions */
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
    int presentation_only;            /**< retained source data changed only by instance visibility/style */
    int progressive_data_complete;     /**< retained root already has complete leaf/source request registry */
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

/**
 * Opaque reference to a temporary draw-frontier promotion.
 *
 * Promotions are owned by one ged instance.  They do not retain that owner,
 * its database, a view, or any renderer object.  Structural draw changes may
 * make an automatic collapse conflict, but do not make the value unsafe to
 * pass to ged_draw_release_promotion().
 */
typedef struct ged_draw_promotion_ref {
    uint64_t owner;
    uint64_t id;
    uint64_t generation;
} ged_draw_promotion_ref;

#define GED_DRAW_PROMOTION_REF_NULL_INIT {0, 0, 0}
#ifdef __cplusplus
#  define GED_DRAW_PROMOTION_REF_NULL ged_draw_promotion_ref{0, 0, 0}
#else
#  define GED_DRAW_PROMOTION_REF_NULL ((ged_draw_promotion_ref){0, 0, 0})
#endif

enum ged_draw_promotion_scope {
    /** Promote only the path occurrence named by the request. */
    GED_DRAW_PROMOTE_EXACT_OCCURRENCE = 0,
    /** Promote every drawn occurrence of the target path's leaf object. */
    GED_DRAW_PROMOTE_ALL_OBJECT_OCCURRENCES = 1
};

enum ged_draw_promotion_outcome {
    GED_DRAW_PROMOTION_CANCEL = 0,
    GED_DRAW_PROMOTION_COMMIT = 1
};

struct ged_draw_promotion_request {
    const char *path;                 /**< borrowed command-facing path */
    const struct db_full_path *dfp;   /**< optional structured path */
    struct ged_view_context *view;    /**< optional view scope */
    int mode;                         /**< draw mode filter; <0 means all */
    int scope;                        /**< enum ged_draw_promotion_scope */
    int auto_draw;                    /**< draw an otherwise-undrawn target */
    const char *role;                 /**< optional borrowed diagnostic role */
};

#define GED_DRAW_PROMOTION_REQUEST_INIT { \
    NULL, NULL, NULL, -1, GED_DRAW_PROMOTE_EXACT_OCCURRENCE, 1, NULL \
}

struct ged_draw_promotion_result {
    int status;
    ged_draw_promotion_ref promotion;
    size_t occurrence_count;
    size_t replaced_root_count;
    size_t replacement_path_count;
    size_t conflict_count;
    uint64_t scene_revision_before;
    uint64_t scene_revision_after;
    struct bu_vls paths;              /**< promoted occurrence paths */
    struct bu_vls errors;
};

__BEGIN_DECLS

typedef struct ged_draw_shape_ref {
    uintptr_t token;
    uint64_t scene_revision;
} ged_draw_shape_ref;

typedef struct ged_draw_group_ref {
    uintptr_t token;
    uint64_t scene_revision;
} ged_draw_group_ref;

/**
 * Opaque value reference to a node in libged's semantic scene inventory.
 *
 * The owner, id, and generation fields are identifiers, not recoverable
 * pointers.  A reference does not retain its GED owner and all operations
 * require that owner explicitly.  Structural scene changes invalidate a
 * reference conservatively; stale or cross-owner references are rejected
 * before any backing record is accessed.
 */
typedef struct ged_scene_node_ref {
    uint64_t owner;
    uint64_t id;
    uint64_t generation;
    uint64_t scene_revision;
} ged_scene_node_ref;

#define GED_SCENE_NODE_REF_NULL_INIT {0, 0, 0, 0}
#ifdef __cplusplus
#  define GED_SCENE_NODE_REF_NULL ged_scene_node_ref{0, 0, 0, 0}
#else
#  define GED_SCENE_NODE_REF_NULL ((ged_scene_node_ref){0, 0, 0, 0})
#endif

struct ged_draw_view_line_style {
    int color[3];
    int line_width;
};

struct ged_view_feature_style {
    int visible;
    int selectable;
    int color_valid;
    unsigned char color[3];
    int line_width;
    int line_style;
    int arrow;
    fastf_t arrow_tip_length;
    fastf_t arrow_tip_width;
};

#define GED_VIEW_FEATURE_STYLE_INIT { -1, -1, 0, {0, 0, 0}, -1, -1, -1, -1.0, -1.0 }

/**
 * Opaque value reference to a managed view feature.
 *
 * The values identify the issuing view owner, the object, and the lifetime of
 * the backing store; none is a recoverable pointer.  References do not retain
 * their owner.  Operations reject a reference after feature removal or
 * recreation, endpoint/controller replacement, view teardown, or store
 * replacement.  Resolve and mutation calls belong on the view owner thread.
 * Callers must copy and compare the complete value and must not interpret its
 * fields individually.
 */
typedef struct ged_view_feature_ref {
    uint64_t owner;
    uint64_t id;
    uint64_t generation;
} ged_view_feature_ref;

#define GED_VIEW_FEATURE_REF_NULL_INIT {0, 0, 0}
#ifdef __cplusplus
#  define GED_VIEW_FEATURE_REF_NULL ged_view_feature_ref{0, 0, 0}
#else
#  define GED_VIEW_FEATURE_REF_NULL ((ged_view_feature_ref){0, 0, 0})
#endif

enum ged_view_feature_family {
    GED_VIEW_FEATURE_UNKNOWN = 0,
    GED_VIEW_FEATURE_TRANSIENT_PREVIEW = 1
};

enum ged_view_edit_preview_event {
    GED_VIEW_EDIT_PREVIEW_BEGIN = 0,
    GED_VIEW_EDIT_PREVIEW_UPDATE,
    GED_VIEW_EDIT_PREVIEW_COMMIT,
    GED_VIEW_EDIT_PREVIEW_CANCEL,
    GED_VIEW_EDIT_PREVIEW_REPLACE_SOURCE,
    GED_VIEW_EDIT_PREVIEW_DISCARD
};

struct ged_view_edit_transaction {
    enum ged_view_edit_preview_event event;
    ged_view_feature_ref feature;
    const char *feature_name;
    const void *owner;
    const char *source_path;
    const char *edit_intent_id;
    const char *edit_intent_role;
    const point_t *points;
    const int *commands;
    size_t point_count;
    struct db_i *dbip;
    struct rt_db_internal *internal;
    const fastf_t *matrix;
    const struct bg_tess_tol *ttol;
    const struct bn_tol *tol;
    uint32_t source_revision;
    uint32_t inputs_revision;
    int color_valid;
    unsigned char color[3];
};

#define GED_VIEW_EDIT_TRANSACTION_INIT { \
    GED_VIEW_EDIT_PREVIEW_UPDATE, GED_VIEW_FEATURE_REF_NULL, \
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, \
    NULL, 0, 0, 0, {255, 255, 255} }

struct ged_view_feature_label {
    const char *text;
    point_t point;
    int color_valid;
    unsigned char color[3];
    fastf_t font_size;
};

struct ged_annotation_label {
    const char *text;
    point_t point;
    int color_valid;
    unsigned char color[3];
    int line_flag;
    point_t target;
    int anchor;
    int arrow;
    fastf_t font_size;
};

#define GED_ANNOTATION_LABEL_INIT { NULL, VINIT_ZERO, 0, {0, 0, 0}, 0, VINIT_ZERO, 0, 0, 0.0 }

struct ged_annotation_axes {
    point_t position;
    fastf_t size;
    int line_width;
    int color[3];
};

struct ged_annotation_line_layer {
    const char *name;
    const point_t *points;
    const int *commands;
    size_t point_count;
    struct ged_view_feature_style style;
};

#define GED_ANNOTATION_LINE_LAYER_INIT { NULL, NULL, NULL, 0, GED_VIEW_FEATURE_STYLE_INIT }

struct ged_result_scene;

enum ged_result_status {
    GED_RESULT_NONE = 0,
    GED_RESULT_ACCEPTED,
    GED_RESULT_UPDATED,
    GED_RESULT_REMOVED,
    GED_RESULT_FAILED
};

struct ged_result_event {
    int status;
    const char *feature_name;
    const char *command;
    const char *diagnostic;
    uint64_t feature_id;
    uint64_t feature_revision;
};

#define GED_RESULT_INIT { GED_RESULT_NONE, NULL, NULL, NULL, 0, 0 }

typedef void (*ged_result_callback)(
    const struct ged_result_event *result,
    void *data);

struct ged_result_custom_node_request {
    const char *feature_name;
    const char *owner_id;
    const char *owner_role;
    uint64_t generation;
    int local;
};

#define GED_RESULT_SCENE_CUSTOM_NODE_REQUEST_INIT { NULL, NULL, NULL, 0, 0 }

typedef void *(*ged_result_custom_node_cb)(
    const struct ged_result_custom_node_request *request,
    void *data);

struct ged_result_desc {
    const char *owner_id;
    const char *owner_role;
    const char *run_id;
    uint64_t generation;
    int local;
    ged_result_callback result_cb;
    void *result_cb_data;
};

#define GED_RESULT_SCENE_DESC_INIT { NULL, NULL, NULL, 0, 0, NULL, NULL }

struct ged_result_metadata {
    const char *key;
    const char *value;
};

#define GED_RESULT_SCENE_METADATA_INIT { NULL, NULL }

enum ged_selection_snap_kind {
    GED_SELECTION_SNAP_GRID = 1,
    GED_SELECTION_SNAP_ENDPOINT = 2
};

typedef struct bv_lod_policy ged_draw_source_lod_policy;

typedef int (*ged_view_feature_depth_cb)(fastf_t depth, void *data);
typedef int (*ged_view_selection_visit_cb)(struct ged_view_context *view_ctx,
	const char *path,
	void *data);
typedef void (*ged_view_selection_path_cb)(const char *path, void *data);

enum ged_polygon_type {
    GED_POLYGON_GENERAL = 0,
    GED_POLYGON_CIRCLE = 1,
    GED_POLYGON_ELLIPSE = 2,
    GED_POLYGON_RECTANGLE = 3,
    GED_POLYGON_SQUARE = 4
};

enum ged_polygon_update {
    GED_POLYGON_UPDATE_DEFAULT = 0,
    GED_POLYGON_UPDATE_PROPS_ONLY = 1,
    GED_POLYGON_UPDATE_PT_SELECT = 2,
    GED_POLYGON_UPDATE_PT_SELECT_CLEAR = 3,
    GED_POLYGON_UPDATE_PT_MOVE = 4,
    GED_POLYGON_UPDATE_PT_APPEND = 5
};

/**
 * Opaque value reference to a managed view polygon.
 *
 * It has the same lifetime and owner-thread rules as
 * ged_view_feature_ref.  In particular, it contains no pointer and does
 * not keep a view, endpoint, controller, store, or polygon alive.
 */
typedef struct ged_polygon_ref {
    uint64_t owner;
    uint64_t id;
    uint64_t generation;
} ged_polygon_ref;

#define GED_POLYGON_REF_NULL_INIT {0, 0, 0}

#ifdef __cplusplus
#  define GED_POLYGON_REF_NULL ged_polygon_ref{0, 0, 0}
#else
#  define GED_POLYGON_REF_NULL ((ged_polygon_ref){0, 0, 0})
#endif

struct ged_polygon_record {
    ged_polygon_ref ref;
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

typedef int (*ged_polygon_record_cb)(
	ged_polygon_ref ref,
	const struct ged_polygon_record *record,
	void *data);

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
    const char *name;                    /**< borrowed display name */
    const struct bn_tol *tol;            /**< borrowed draw source tolerance */
    const struct bg_tess_tol *ttol;      /**< borrowed draw source tessellation tolerance */
};

struct ged_draw_group_record {
    ged_draw_group_ref ref;
    const char *path;                    /**< borrowed */
    const struct db_full_path *fullpath; /**< borrowed; valid while draw scene is unchanged */
    struct ged_view_context *view;       /**< borrowed; NULL for view-independent groups */
    struct ged_draw_appearance_settings appearance;
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
    GED_DRAW_VIEW_LINE_DRAW = 1,
    GED_DRAW_VIEW_LINE_POINT_DRAW = 12
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
    int selected;
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

struct ged_draw_view_rendered_object_summary_s {
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
typedef struct ged_draw_view_rendered_object_summary_s ged_draw_view_rendered_object_summary_t;

enum ged_view_feature_kind {
    GED_VIEW_FEATURE_KIND_UNKNOWN = 0,
    GED_VIEW_FEATURE_KIND_LINES,
    GED_VIEW_FEATURE_KIND_INDEXED_LINES,
    GED_VIEW_FEATURE_KIND_POINTS,
    GED_VIEW_FEATURE_KIND_LABELS,
    GED_VIEW_FEATURE_KIND_ARROW,
    GED_VIEW_FEATURE_KIND_AXES,
    GED_VIEW_FEATURE_KIND_LINE_LAYER,
    GED_VIEW_FEATURE_KIND_EDIT_PREVIEW,
    GED_VIEW_FEATURE_KIND_INDEXED_FACE_SET,
    GED_VIEW_FEATURE_KIND_POLYGON_OVERLAY,
    GED_VIEW_FEATURE_KIND_HUD_LABEL,
    GED_VIEW_FEATURE_KIND_CUSTOM_NODE
};

enum ged_view_feature_scope {
    GED_VIEW_FEATURE_SCOPE_UNKNOWN = 0,
    GED_VIEW_FEATURE_SCOPE_SHARED,
    GED_VIEW_FEATURE_SCOPE_LOCAL
};

enum ged_view_feature_overlay_class {
    GED_VIEW_FEATURE_OVERLAY_CLASS_UNKNOWN = 0,
    GED_VIEW_FEATURE_OVERLAY_CLASS_NONE,
    GED_VIEW_FEATURE_OVERLAY_CLASS_FACEPLATE,
    GED_VIEW_FEATURE_OVERLAY_CLASS_EDIT_HANDLE,
    GED_VIEW_FEATURE_OVERLAY_CLASS_MEASURE,
    GED_VIEW_FEATURE_OVERLAY_CLASS_SELECTION_RUBBER_BAND,
    GED_VIEW_FEATURE_OVERLAY_CLASS_SNAP_GUIDE,
    GED_VIEW_FEATURE_OVERLAY_CLASS_COMMAND_RESULT,
    GED_VIEW_FEATURE_OVERLAY_CLASS_DIAGNOSTIC,
    GED_VIEW_FEATURE_OVERLAY_CLASS_TCL_OVERLAY,
    GED_VIEW_FEATURE_OVERLAY_CLASS_POLYGON_EDIT,
    GED_VIEW_FEATURE_OVERLAY_CLASS_SKETCH_EDIT,
    GED_VIEW_FEATURE_OVERLAY_CLASS_USER_ANNOTATION
};

enum ged_view_feature_lifecycle {
    GED_VIEW_FEATURE_LIFECYCLE_UNKNOWN = 0,
    GED_VIEW_FEATURE_LIFECYCLE_NONE,
    GED_VIEW_FEATURE_LIFECYCLE_PERSISTENT,
    GED_VIEW_FEATURE_LIFECYCLE_PER_FRAME,
    GED_VIEW_FEATURE_LIFECYCLE_PER_COMMAND,
    GED_VIEW_FEATURE_LIFECYCLE_PER_TOOL,
    GED_VIEW_FEATURE_LIFECYCLE_PER_VIEW,
    GED_VIEW_FEATURE_LIFECYCLE_SHARED_VIEW_SET,
    GED_VIEW_FEATURE_LIFECYCLE_AUTO_REMOVE_ON_SOURCE
};

#define GED_VIEW_FEATURE_OWNER_ID_MAX 64
#define GED_VIEW_FEATURE_OWNER_ROLE_MAX 64

struct ged_view_feature_summary {
    int exists;
    int is_overlay;
    int is_label;
    int is_transient_preview;
    int is_command_result;
    int visible;
    int kind;
    int scope;
    int overlay_class;
    int lifecycle;
    unsigned char color[3];
    size_t child_count;
    size_t geometry_command_count;
    size_t metadata_count;
    size_t primitive_metadata_count;
    size_t selected_primitive_count;
    size_t highlighted_primitive_count;
    uint64_t owner_generation;
    char owner_id[GED_VIEW_FEATURE_OWNER_ID_MAX];
    char owner_role[GED_VIEW_FEATURE_OWNER_ROLE_MAX];
};

#define GED_VIEW_FEATURE_SUMMARY_INIT { 0, 0, 0, 0, 0, 0, GED_VIEW_FEATURE_KIND_UNKNOWN, GED_VIEW_FEATURE_SCOPE_UNKNOWN, GED_VIEW_FEATURE_OVERLAY_CLASS_UNKNOWN, GED_VIEW_FEATURE_LIFECYCLE_UNKNOWN, {0, 0, 0}, 0, 0, 0, 0, 0, 0, 0, {'\0'}, {'\0'} }

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

struct ged_draw_shape_candidate {
    const char *path;
    const char *instance_key;
    int draw_mode;
};

typedef int (*ged_draw_shape_candidate_cb)(
    const struct ged_draw_shape_candidate *candidate,
    void *userdata);

#define GED_DRAW_VIEW_RECORD_QUERY_VIEW_OBJECTS 0x01u
#define GED_DRAW_VIEW_RECORD_QUERY_DB_OBJECTS   0x02u
#define GED_DRAW_VIEW_RECORD_QUERY_LOCAL_ONLY   0x04u

struct ged_draw_view_record_query {
    unsigned int flags;
    const char *glob;
    int draw_mode;                       /**< negative means any draw mode */
};


__END_DECLS

#endif /* GED_DRAW_TYPES_H */
