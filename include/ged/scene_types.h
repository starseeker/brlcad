/* S C E N E _ T Y P E S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 2.1.
 */
/** @file ged/scene_types.h
 *
 * Renderer-neutral semantic scene request types.
 */

#ifndef GED_SCENE_TYPES_H
#define GED_SCENE_TYPES_H

#include <stddef.h>
#include <stdint.h>

#include "vmath.h"

struct ged;
struct db_full_path;
struct ged_scene_delta;
struct ged_scene_result;
struct ged_view_context;

/** Status returned by semantic scene operations. */
enum ged_scene_status {
    GED_SCENE_OK = 0,
    GED_SCENE_INVALID,
    GED_SCENE_NOT_FOUND,
    GED_SCENE_ERROR
};

/** Database geometry presentation mode. */
enum ged_scene_draw_mode {
    GED_SCENE_DRAW_DEFAULT = -1,
    GED_SCENE_DRAW_WIRE = 0,
    GED_SCENE_DRAW_SHADED_BOTS,
    GED_SCENE_DRAW_SHADED,
    GED_SCENE_DRAW_EVALUATED_WIRE,
    GED_SCENE_DRAW_HIDDEN_LINE,
    GED_SCENE_DRAW_EVALUATED_POINTS
};

/** Coverage of a database path by the current semantic scene. */
enum ged_scene_path_state {
    GED_SCENE_PATH_NOT_DRAWN = 0,
    GED_SCENE_PATH_DRAWN,
    GED_SCENE_PATH_PARTIALLY_DRAWN
};

/** Path identity requested from ged_scene_paths_append(). */
enum ged_scene_path_listing {
    GED_SCENE_PATHS_DRAW_INTENTS = 0,
    GED_SCENE_PATHS_REALIZED_OCCURRENCES
};

/**
 * Stable identity of one database occurrence represented by the scene.
 *
 * The fields are opaque.  References remain valid across unrelated scene
 * commits and become invalid only when their occurrence is retired or their
 * owning GED context is destroyed.  Clients must not compare or interpret
 * individual fields; use ged_scene_occurrence_ref_equal().
 */
typedef struct ged_scene_occurrence_ref {
    uintptr_t owner;
    uint64_t id;
    uint64_t generation;
} ged_scene_occurrence_ref;

#define GED_SCENE_OCCURRENCE_REF_NULL_INIT {0, 0, 0}
#ifdef __cplusplus
#  define GED_SCENE_OCCURRENCE_REF_NULL ged_scene_occurrence_ref{0, 0, 0}
#else
#  define GED_SCENE_OCCURRENCE_REF_NULL ((ged_scene_occurrence_ref){0, 0, 0})
#endif

/**
 * Callback-lifetime semantic snapshot of a realized database occurrence.
 *
 * Strings and @p fullpath are borrowed and must not be retained.  This record
 * deliberately contains semantic identity and presentation state only; mesh,
 * vlist, renderer-node, cache, and LoD storage are backend-private.
 */
struct ged_scene_occurrence_info {
    ged_scene_occurrence_ref ref;
    const struct db_full_path *fullpath;
    const char *path;
    const char *leaf_name;
    unsigned long long path_hash;
    struct ged_view_context *view;
    enum ged_scene_draw_mode draw_mode;
    double opacity;
    int visible;
    int highlighted;
    int selected;
    int evaluated_region;
    int line_width;
    point_t center;
};

/** Lightweight visible-occurrence candidate used by interactive navigation. */
struct ged_scene_occurrence_candidate {
    const char *path;
    const char *instance_key;
    enum ged_scene_draw_mode draw_mode;
};

typedef int (*ged_scene_occurrence_func_t)(
    const struct ged_scene_occurrence_info *occurrence,
    void *client_data);

typedef int (*ged_scene_occurrence_candidate_func_t)(
    const struct ged_scene_occurrence_candidate *candidate,
    void *client_data);

/** Geometry classes included in a scene-bounds query. */
enum ged_scene_bounds_scope {
    GED_SCENE_BOUNDS_DATABASE = 0,
    GED_SCENE_BOUNDS_ALL
};

/** How a path-scoped presentation request matches database occurrences. */
enum ged_scene_path_match {
    GED_SCENE_PATH_MATCH_EXACT = 0,
    GED_SCENE_PATH_MATCH_SUBTREE,
    GED_SCENE_PATH_MATCH_OBJECT
};

/** Renderer realization policy for a draw request. */
enum ged_scene_realization_mode {
    GED_SCENE_REALIZE_AUTO = 0,
    GED_SCENE_REALIZE_EAGER,
    GED_SCENE_REALIZE_PROGRESSIVE
};

/** Flags modifying database scene presentation. */
enum ged_scene_style_flag {
    GED_SCENE_STYLE_SOLID_LINES_ONLY = 1u << 0,
    GED_SCENE_STYLE_NON_SUBTRACT_ONLY = 1u << 1
};

/** Renderer-neutral database scene style. */
struct ged_scene_style {
    enum ged_scene_draw_mode draw_mode;
    double opacity;
    unsigned char color[3];
    int color_override;
    int line_width;
    int mixed_modes;
    unsigned flags;
};

/** Controls how database geometry is initially realized. */
struct ged_scene_realization_policy {
    enum ged_scene_realization_mode mode;
    int strict;
};

/** Typed request for drawing one or more database paths. */
struct ged_scene_draw_request {
    struct ged_view_context *view;
    const char *const *paths;
    size_t path_count;
    struct ged_scene_style style;
    struct ged_scene_realization_policy realization;
    int autoview;
};

/** Matching rule for an erase operation. */
enum ged_scene_erase_match {
    GED_SCENE_ERASE_EXACT = 0,
    GED_SCENE_ERASE_SUBTREE
};

/** Typed request for erasing a database path from the semantic scene. */
struct ged_scene_erase_request {
    struct ged_view_context *view;
    const char *path;
    enum ged_scene_erase_match match;
    enum ged_scene_draw_mode draw_mode;
};

/** Typed path scope shared by visibility, opacity, and highlight updates. */
struct ged_scene_path_request {
    struct ged_view_context *view;
    const char *path;
    enum ged_scene_draw_mode draw_mode;
    enum ged_scene_path_match match;
};

/** Typed request to rebuild all or selected retained draw sources. */
struct ged_scene_redraw_request {
    struct ged_view_context *view;
    const char *const *paths;
    size_t path_count;
};

/** Semantic scene storage selected by a clear request. */
enum ged_scene_clear_scope {
    GED_SCENE_CLEAR_ALL = 0,
    GED_SCENE_CLEAR_VIEW
};

/** Typed request to clear canonical or view-scoped database draw intent. */
struct ged_scene_clear_request {
    struct ged_view_context *view;
    enum ged_scene_clear_scope scope;
};

/** Opaque, generation-checked reference to an active semantic edit scope. */
typedef struct ged_scene_edit_scope_ref {
    uint64_t owner;
    uint64_t id;
    uint64_t generation;
} ged_scene_edit_scope_ref;

#define GED_SCENE_EDIT_SCOPE_REF_NULL_INIT {0, 0, 0}
#ifdef __cplusplus
#  define GED_SCENE_EDIT_SCOPE_REF_NULL ged_scene_edit_scope_ref{0, 0, 0}
#else
#  define GED_SCENE_EDIT_SCOPE_REF_NULL ((ged_scene_edit_scope_ref){0, 0, 0})
#endif

/** Database occurrences made independently addressable by an edit scope. */
enum ged_scene_edit_occurrence_scope {
    GED_SCENE_EDIT_EXACT_OCCURRENCE = 0,
    GED_SCENE_EDIT_ALL_DRAWN_OCCURRENCES
};

/** How an edit scope is closed. */
enum ged_scene_edit_outcome {
    GED_SCENE_EDIT_CANCEL = 0,
    GED_SCENE_EDIT_COMMIT
};

/** Request an independently addressable semantic scene scope for editing. */
struct ged_scene_edit_request {
    struct ged_view_context *view;
    const char *path;
    enum ged_scene_draw_mode draw_mode;
    enum ged_scene_edit_occurrence_scope occurrences;
    int draw_if_absent;
    const char *purpose;
};

/** Semantic change described by a committed scene delta. */
enum ged_scene_delta_kind {
    GED_SCENE_DELTA_UNKNOWN = 0,
    GED_SCENE_DELTA_DRAW,
    GED_SCENE_DELTA_ERASE,
    GED_SCENE_DELTA_REDRAW,
    GED_SCENE_DELTA_VISIBILITY,
    GED_SCENE_DELTA_HIGHLIGHT,
    GED_SCENE_DELTA_STYLE,
    GED_SCENE_DELTA_EDIT_SCOPE,
    GED_SCENE_DELTA_SOURCE,
    GED_SCENE_DELTA_CLEAR,
    GED_SCENE_DELTA_TEARDOWN
};

/** Opaque token identifying a semantic scene observer registration. */
typedef uintptr_t ged_scene_observer_token;

/**
 * Callback invoked synchronously after a semantic scene commit.
 *
 * The delta and all values borrowed from it are valid only for the duration
 * of the callback.  A callback may remove its own observer registration.
 */
typedef void (*ged_scene_observer_func_t)(
    struct ged *gedp,
    const struct ged_scene_delta *delta,
    void *client_data);

#endif /* GED_SCENE_TYPES_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
