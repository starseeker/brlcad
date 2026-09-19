/* S C E N E _ R E C O R D _ T Y P E S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged_scene_record_private.h
 *
 * Private renderer-adapter records used while libged realizes semantic scene
 * state.  These values are never part of the installed GED API.
 */

#ifndef GED_SCENE_RECORD_PRIVATE_H
#define GED_SCENE_RECORD_PRIVATE_H

#include <stddef.h>
#include <stdint.h>

#include "vmath.h"
#include "./ged_draw_types_private.h"

struct bg_tess_tol;
struct bn_tol;
struct db_full_path;
struct db_i;
struct directory;
struct ged_view_context;

typedef struct ged_draw_shape_ref {
    uintptr_t token;
    uint64_t scene_revision;
} ged_draw_shape_ref;

typedef struct ged_draw_group_ref {
    uintptr_t token;
    uint64_t scene_revision;
} ged_draw_group_ref;

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
    const struct db_full_path *fullpath;
    const char *display_name;
    const char *leaf_name;
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
    const struct db_full_path *fullpath;
    struct directory *leaf_dp;
    const char *name;
    const struct bn_tol *tol;
    const struct bg_tess_tol *ttol;
};

struct ged_draw_group_record {
    ged_draw_group_ref ref;
    const char *path;
    const struct db_full_path *fullpath;
    struct ged_view_context *view;
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
    const char *path;
    const char *type_name;
    const char *geometry_name;
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
    uintptr_t detail_token;
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
    const char *text;
    uint64_t cache_identity;
    uint64_t source_identity;
};

struct ged_draw_view_line_summary {
    int valid;
    size_t point_count;
    uint64_t cache_identity;
    uint64_t source_identity;
};

struct ged_draw_shape_geometry_summary {
    int valid;
    const char *geometry_name;
    size_t point_count;
    size_t index_count;
};

struct ged_draw_scene_display_summary {
    int valid;
    int is_database_source;
    int has_draw_intent;
    const char *intent_path;
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
    const char *database_path;
    const char *stale_reason;
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
    const char *name;
    const struct db_full_path *fullpath;
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
    int draw_mode;
};

typedef int (*ged_draw_view_db_object_record_cb)(
    const struct ged_draw_view_db_object_record *record,
    void *client_data);

typedef int (*ged_draw_view_segment_cb)(const point_t start,
	const point_t end, void *client_data);

typedef int (*ged_draw_view_point_cb)(const point_t point,
	void *client_data);

typedef int (*ged_draw_shape_ref_index_cb)(ged_draw_shape_ref ref,
	void *client_data);

typedef int (*ged_draw_shape_candidate_cb)(
    const struct ged_draw_shape_candidate *candidate,
    void *client_data);

typedef struct ged_draw_view_rendered_object_summary_s
    ged_draw_view_rendered_object_summary_t;

#endif /* GED_SCENE_RECORD_PRIVATE_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
