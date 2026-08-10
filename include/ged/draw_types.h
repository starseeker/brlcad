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
 * libged owns the semantic draw records and sends committed changes to a
 * private retained-render backend.  Installed GED callers use this API rather
 * than renderer scene-graph nodes to query or update the draw scene.
 *
 * Drawn database objects are keyed by struct db_full_path.  String path
 * formatting is available for logging/UI, but graph mutation APIs use the
 * structured path form.
 *
 * Hosts that attach a display backend use the renderer-neutral operations in
 * ged/display.h.
 */
/** @{ */
/* @file ged/draw_types.h */

#ifndef GED_DRAW_TYPES_H
#define GED_DRAW_TYPES_H

#include <stddef.h>

#include "vmath.h"
#include "bg/polygon_types.h"
#include "bu/color.h"
#include "bu/list.h"
#include "bu/ptbl.h"
#include "bu/vls.h"
#include "bv/view.h"
#include "ged/defines.h"
#include "ged/polygon_types.h"
#include "ged/selection_types.h"
#include "ged/scene_record_types.h"
#include "ged/view_feature_types.h"

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
    int defer_leaf_expansion;        /**< Retain the requested root while leaves realize progressively. */
};

#define GED_DRAW_APPEARANCE_SETTINGS_INIT {GED_DRAW_MODE_WIRE, 0, 1.0, 0, {255, 0, 0}, 1, 0.0, 0.0, 0, 0, 0, 0}

__BEGIN_DECLS

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

struct ged_draw_shape_candidate {
    const char *path;
    const char *instance_key;
    int draw_mode;
};

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
