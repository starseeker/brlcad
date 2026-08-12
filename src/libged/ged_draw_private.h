/*                  G E D _ D R A W _ P R I V A T E . H
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
/** @file ged_draw_private.h
 *
 * Private node-facing draw-scene bridge for libged implementation files.
 *
 * Installed GED drawing APIs are record/ref/transaction based.  These helpers
 * are intentionally not installed: public callers use semantic records and
 * refs while these helpers connect libged to its Obol scene controller.
 */

#ifndef LIBGED_GED_DRAW_PRIVATE_H
#define LIBGED_GED_DRAW_PRIVATE_H

#include "common.h"

#include <stdint.h>

#include "bu/vls.h"
#include "vmath.h"

#include "BObol/BMeshLodCache.h"
#include "ged/draw.h"
#include "./ged_draw_source_private.h"
#include "ged/draw_scene.h"
#include "ged/display.h"
#include "ged/scene.h"
#include "ged/selection.h"
#include "rt/db_fullpath.h"
#include "rt/view.h"

#include "./ged_scene_reducer_private.h"

__BEGIN_DECLS

struct bu_ptbl;
struct bu_vls;
struct bg_tess_tol;
struct bn_tol;
struct db_i;
struct directory;
struct rt_db_internal;
struct rt_bot_internal;
struct rt_brep_internal;
struct rt_pg_internal;
struct db_tree_state;
struct model;
struct nmgregion;
struct bobol_display_endpoint;
struct BObolDrawMetadataRecord;
struct ged_draw_obol_database_source_record;

extern void ged_view_feature_info_get(struct bv_view_info *view_info,
	const struct ged_view_context *view_ctx);

extern uint64_t ged_draw_material_revision(const struct ged *gedp);
extern void ged_draw_bump_material_revision(struct ged *gedp);
extern uint64_t ged_draw_scene_revision(struct ged *gedp);
extern int ged_draw_scene_available(struct ged *gedp);
extern int ged_draw_default_mode(const struct ged *gedp);
extern int ged_draw_bounds(struct ged *gedp, vect_t *min, vect_t *max,
	int include_overlays);
extern uint64_t ged_draw_highlight_revision(const struct ged *gedp);
extern uint64_t ged_draw_highlight_revision_advance(struct ged *gedp);
extern void ged_draw_set_highlighted_shape_ref(struct ged *gedp,
	ged_draw_shape_ref ref);
extern int ged_draw_shape_set_highlighted(struct ged *gedp,
	ged_draw_shape_ref ref, int highlighted);
extern int ged_draw_has_shapes(struct ged *gedp);
extern int ged_draw_path_state(struct ged *gedp,
	struct ged_view_context *view_ctx, const char *path, int mode);
extern size_t ged_draw_list_paths(struct ged *gedp,
	struct ged_view_context *view_ctx, int mode, int expanded,
	struct bu_vls *result);
extern int ged_draw_has_paths(struct ged *gedp,
	struct ged_view_context *view_ctx, int mode);
extern int ged_draw_shape_count(struct ged *gedp);
extern int ged_draw_has_groups(struct ged *gedp);
extern void ged_draw_clear(struct ged *gedp);
extern int ged_draw_clear_view(struct ged *gedp,
	struct ged_view_context *view_ctx);

/* Obol traversal identity used only inside the backend adapter and its
 * white-box tests.  Semantic clients must use draw paths and typed scene
 * deltas rather than renderer-node handles. */
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

extern int ged_scene_node_ref_is_null(ged_scene_node_ref ref);
extern int ged_scene_node_ref_equal(ged_scene_node_ref a,
	ged_scene_node_ref b);
extern ged_scene_node_ref ged_scene_first_node(struct ged *gedp);
extern ged_scene_node_ref ged_scene_shape_node(struct ged *gedp,
	ged_draw_shape_ref ref);
extern ged_scene_node_ref ged_scene_shape_cache_node(struct ged *gedp,
	ged_draw_shape_ref ref);
extern ged_scene_node_ref ged_scene_group_node(struct ged *gedp,
	ged_draw_group_ref ref);
extern ged_draw_shape_ref ged_scene_shape_ref(struct ged *gedp,
	ged_scene_node_ref node);
extern int ged_scene_node_line_summary(struct ged *gedp,
	ged_scene_node_ref node, struct ged_draw_view_line_summary *out);
extern int ged_scene_node_line_point_at(struct ged *gedp,
	ged_scene_node_ref node, size_t index, point_t out);
extern int ged_scene_node_line_command_at(struct ged *gedp,
	ged_scene_node_ref node, size_t index, int *out);
extern int ged_scene_node_geometry_summary(struct ged *gedp,
	ged_scene_node_ref node, struct ged_draw_shape_geometry_summary *out);
extern int ged_scene_node_has_state(struct ged *gedp,
	ged_scene_node_ref node);
extern ged_scene_node_ref ged_scene_node_source(struct ged *gedp,
	ged_scene_node_ref node);
extern const struct db_full_path *ged_scene_node_fullpath(struct ged *gedp,
	ged_scene_node_ref node);
extern int ged_scene_group_dbpath(struct ged *gedp,
	ged_scene_node_ref group, struct db_full_path *out);
extern int ged_scene_group_is_overlay(struct ged *gedp,
	ged_scene_node_ref group);
extern int ged_scene_node_display_summary(struct ged *gedp,
	ged_scene_node_ref node, struct ged_draw_scene_display_summary *out);
extern int ged_scene_node_tree_summary(struct ged *gedp,
	ged_scene_node_ref node, struct ged_draw_scene_tree_summary *out);
extern ged_scene_node_ref ged_scene_node_child_at(struct ged *gedp,
	ged_scene_node_ref node, size_t index);
extern ged_scene_node_ref ged_scene_node_parent(struct ged *gedp,
	ged_scene_node_ref node);
extern const char *ged_scene_node_name(struct ged *gedp,
	ged_scene_node_ref node);
extern int ged_scene_node_subtree_bounds(struct ged *gedp,
	ged_scene_node_ref node, vect_t *min, vect_t *max,
	int include_overlays);
extern int ged_draw_source_node_summary(struct ged *gedp,
	ged_scene_node_ref node,
	struct ged_draw_database_source_summary *out);

typedef int (*ged_draw_group_ref_index_cb)(ged_draw_group_ref ref,
					   void *userdata);
typedef int (*ged_draw_obol_group_path_cb)(
	struct ged *gedp,
	const char *path,
	void *userdata);
typedef int (*ged_draw_obol_database_source_path_cb)(
	struct ged *gedp,
	const char *path,
	void *userdata);
typedef int (*ged_draw_obol_database_source_record_cb)(
	struct ged *gedp,
	const struct ged_draw_obol_database_source_record *record,
	void *userdata);

enum ged_draw_overlay_geometry_kind {
    GED_DRAW_OVERLAY_GEOMETRY_NONE = 0,
    GED_DRAW_OVERLAY_GEOMETRY_LINE_SET,
    GED_DRAW_OVERLAY_GEOMETRY_POINT_SET,
    GED_DRAW_OVERLAY_GEOMETRY_INDEXED_FACE_SET
};

struct ged_draw_obol_realization_policy_summary {
    int valid;
    int role_flags;
    int view_dependent;
    int csg_lod_enabled;
    int mesh_lod_enabled;
    fastf_t view_scale;
    fastf_t lod_scale;
    uint64_t bot_threshold;
    fastf_t curve_scale;
    fastf_t point_scale;
};

enum {
    GED_DRAW_OBOL_DATABASE_SOURCE_MATERIAL_INHERIT = 0,
    GED_DRAW_OBOL_DATABASE_SOURCE_MATERIAL_DATABASE = 1
};

enum {
    GED_DRAW_OBOL_DATABASE_SOURCE_REALIZATION_STALE = 0,
    GED_DRAW_OBOL_DATABASE_SOURCE_REALIZATION_CURRENT = 1
};

/* Renderer-adapter service records.  These are deliberately private: public
 * integrations use the view-lod command and semantic draw/source APIs. */
struct ged_draw_obol_lod_service_status_s {
    int attached;
    int running;
    int auto_submit;
    size_t worker_count;
    size_t in_flight;
    size_t pending_tasks;
    size_t queued_results;
    size_t queued_cache_writes;
    size_t delayed_tasks;
    unsigned int last_visited_mesh_count;
    unsigned int last_submitted_task_count;
    unsigned int last_updated_cut_count;
    unsigned int last_skipped_mesh_count;
    size_t last_result_count;
    unsigned int last_matched_result_count;
    unsigned int last_applied_result_count;
    unsigned int last_rejected_result_count;
    unsigned int last_unmatched_result_count;
    size_t active_mesh_payloads;
    size_t active_aabb_proxy_payloads;
    size_t active_obb_proxy_payloads;
    char last_diagnostics[2048];
};
typedef struct ged_draw_obol_lod_service_status_s ged_draw_obol_lod_service_status_t;

struct ged_draw_obol_source_expansion_status {
    size_t child_count;
    size_t considered;
    size_t expanded;
    size_t existing;
    size_t skipped_non_union;
    size_t skipped_duplicate_instance;
    size_t expanded_non_union;
    size_t expanded_duplicate_instance;
    size_t skipped_invalid;
    size_t remaining;
    size_t proxy_published;
    size_t metadata_applied;
    size_t comb_sources;
    size_t leaf_sources;
};

struct ged_draw_obol_source_prewarm_status {
    size_t child_count;
    size_t considered;
    size_t submitted;
    size_t already_cached;
    size_t skipped_non_union;
    size_t skipped_duplicate_instance;
    size_t shared_request;
    size_t non_union_children;
    size_t duplicate_instances;
    size_t skipped_invalid;
    size_t remaining;
    size_t comb_sources;
    size_t leaf_sources;
};

struct ged_draw_obol_database_source_record {
    int valid;
    const char *database_path;
    const char *source_path;
    const char *instance_key;
    const char *owner_group_path;
    int visible;
    int selected;
    int highlighted;
    fastf_t transparency;
    int draw_mode;
    int material_policy;
    uint64_t source_revision;
    uint64_t inputs_revision;
    uint64_t realized_source_revision;
    uint64_t realized_inputs_revision;
    ged_draw_stale_reason stale_reason;
    uint64_t realization_identity;
    int realization_status;
    int realization_role_flags;
    int realization_view_dependent;
    int realization_csg_lod_enabled;
    int realization_mesh_lod_enabled;
    fastf_t realization_view_scale;
    fastf_t realization_lod_scale;
    uint64_t realization_bot_threshold;
    fastf_t realization_curve_scale;
    fastf_t realization_point_scale;
};

struct ged_draw_obol_database_source_runtime {
    int valid;
    const char *database_path;
    struct db_i *dbip;
    fastf_t tessellation_abs_tol;
    fastf_t tessellation_rel_tol;
    fastf_t tessellation_norm_tol;
    uint64_t lod_bot_threshold;
    int draw_size_valid;
    fastf_t draw_size;
    struct BObolMeshLod *mesh_lod;
    int mesh_lod_bounds_valid;
    point_t mesh_lod_bmin;
    point_t mesh_lod_bmax;
};

struct ged_draw_obol_draw_state_summary {
    int valid;
    int draw_mode_valid;
    int draw_mode;
    int line_style;
    int draw_mat_valid;
    mat_t draw_mat;
};

struct ged_draw_obol_scene_context_info {
    char *path;
    char *instance_key;
    char *name;
    int node_kind;
    int is_group;
    int is_shape;
    int is_database_source;
    int has_parent;
    int draw_tree_depth;
    size_t child_count;
    int draw_mode_valid;
    int draw_mode;
};

struct ged_draw_overlay_geometry {
    enum ged_draw_overlay_geometry_kind kind;
    const point_t *points;
    size_t point_count;
    const int *commands;
    size_t command_count;
    const vect_t *normals;
    size_t normal_count;
    const int *indices;
    size_t index_count;
};

/* Obol adapter operations.  None of these declarations is installed; public
 * callers use endpoints, semantic source APIs, or the `view lod` command. */
extern int ged_draw_obol_progressive_available(
	struct ged *gedp,
	struct ged_view_context *view_ctx);
extern int ged_draw_obol_view_lod_enabled(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	int draw_mode);
extern void ged_draw_obol_framebuffer_endpoint_detach(
	struct ged *gedp,
	struct bobol_display_endpoint *endpoint);
extern int ged_draw_obol_lod_service_start(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	size_t worker_count);
extern int ged_draw_obol_lod_service_stop(
	struct ged *gedp,
	struct ged_view_context *view_ctx);
extern int ged_draw_obol_lod_service_poll(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	size_t max_results,
	ged_draw_obol_lod_service_status_t *status);
extern int ged_draw_obol_lod_service_status(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	ged_draw_obol_lod_service_status_t *status);
extern int ged_draw_obol_view_lod_policy_changed(
	struct ged *gedp,
	struct ged_view_context *view_ctx);
extern size_t ged_draw_obol_lod_service_prewarm(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	int argc,
	const char * const *argv,
	ged_draw_obol_lod_service_status_t *status);
extern int ged_draw_obol_database_source_expand_children(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	const char *path,
	int ged_draw_mode,
	size_t max_children,
	struct ged_draw_obol_source_expansion_status *status);
extern size_t
ged_draw_obol_database_source_prewarm_child_aabb_proxies(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	const char *path,
	int ged_draw_mode,
	size_t max_children,
	struct ged_draw_obol_source_prewarm_status *status);
extern size_t
ged_draw_obol_database_source_prewarm_visible_child_aabb_proxies(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	const char *root_path,
	int ged_draw_mode,
	size_t max_sources,
	size_t max_children_per_source,
	struct ged_draw_obol_source_prewarm_status *status);
extern int
ged_draw_obol_database_source_expand_visible_children(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	const char *root_path,
	int ged_draw_mode,
	size_t max_sources,
	size_t max_children_per_source,
	struct ged_draw_obol_source_expansion_status *status);
extern int ged_draw_obol_database_source_count(
	struct ged *gedp,
	int skip_overlay_groups,
	size_t *out);
extern int ged_draw_obol_database_source_publish_line_set_for_path(
	struct ged *gedp,
	const char *path,
	const point_t *points,
	const int *commands,
	size_t point_count);
extern int ged_draw_obol_database_source_publish_point_set_for_path(
	struct ged *gedp,
	const char *path,
	const point_t *points,
	size_t point_count);
extern int
ged_draw_obol_database_source_publish_indexed_face_set_for_path(
	struct ged *gedp,
	const char *path,
	const point_t *points,
	size_t point_count,
	const vect_t *normals,
	size_t normal_count,
	const int *indices,
	size_t index_count);
extern int
ged_draw_obol_database_source_publish_lod_indexed_face_set_for_path(
	struct ged *gedp,
	const char *path,
	const point_t *points,
	size_t point_count,
	const vect_t *normals,
	size_t normal_count,
	const int *indices,
	size_t index_count);
extern int
ged_draw_obol_database_source_publish_primitive_wireframe_for_path(
	struct ged *gedp,
	const char *path,
	struct rt_db_internal *ip,
	const struct bg_tess_tol *ttol,
	const struct bn_tol *tol);

extern void ged_draw_overlay_erase_name_context(struct ged *gedp,
							   struct ged_view_context *view_ctx,
							   const char *name);
extern size_t ged_draw_obol_view_context_clear(
	struct ged_view_context *view_ctx,
	int flags);
extern int ged_draw_obol_view_context_selection_available(
	struct ged_view_context *view_ctx);
extern size_t ged_draw_obol_view_context_selection_count(
	struct ged_view_context *view_ctx);
extern int ged_draw_obol_view_context_selection_path_foreach(
	struct ged_view_context *view_ctx,
	ged_view_selection_visit_cb cb,
	void *data);
extern int ged_draw_obol_view_context_selection_clear(
	struct ged_view_context *view_ctx);
extern int ged_draw_obol_view_context_selection_contains_path(
	struct ged_view_context *view_ctx,
	int kind,
	const char *path);
extern int ged_draw_obol_view_context_selection_add_path(
    struct ged_view_context *view_ctx,
    int kind,
    const char *path);
extern int ged_draw_obol_view_context_selection_remove_path(
    struct ged_view_context *view_ctx,
    int kind,
    const char *path);
extern int ged_draw_obol_view_context_selection_set_path(
	struct ged_view_context *view_ctx,
	int kind,
	const char *path);
extern struct ged_pick_result *ged_draw_obol_view_context_pick_point(
	struct ged_view_context *view_ctx,
	int x,
	int y,
	int first_only);
extern struct ged_pick_result *ged_draw_obol_view_context_pick_nearest(
	struct ged_view_context *view_ctx,
	int x,
	int y);
extern struct ged_pick_result *ged_draw_obol_view_context_pick_rect(
	struct ged_view_context *view_ctx,
	int x0,
	int y0,
	int x1,
	int y1);
extern int ged_draw_obol_view_context_snap_first_candidate(
	struct ged_view_context *view_ctx,
	const point_t sample,
	enum ged_selection_snap_kind kind,
	point_t candidate);
extern int ged_draw_obol_view_feature_ref_is_null(
	ged_view_edit_ref ref);
extern int ged_draw_obol_view_feature_remove_ref(
	struct ged_view_context *view_ctx,
	ged_view_edit_ref ref);
extern int ged_draw_obol_view_context_edit_preview_publish_event(
	struct ged_view_context *view_ctx,
	ged_view_edit_ref feature,
	enum ged_view_edit_preview_event event,
	const char *source_path);
extern ged_view_edit_ref
ged_draw_obol_view_context_feature_overlay_ensure(
	struct ged_view_context *view_ctx,
	const char *name,
	const char *source_path);
extern ged_view_edit_ref
ged_draw_obol_view_context_feature_label_ensure(
	struct ged_view_context *view_ctx,
	const char *name);
extern int ged_draw_obol_view_feature_set_visible(
	struct ged_view_context *view_ctx,
	ged_view_edit_ref ref,
	int visible);
extern int ged_draw_obol_view_feature_set_color(
	struct ged_view_context *view_ctx,
	ged_view_edit_ref ref,
	int r,
	int g,
	int b);
extern int ged_draw_obol_view_feature_touch(
	struct ged_view_context *view_ctx,
	ged_view_edit_ref ref);
extern int ged_draw_obol_view_feature_labels_replace(
	struct ged_view_context *view_ctx,
	ged_view_edit_ref ref,
	const struct ged_view_feature_label *labels,
	size_t label_count);
extern int ged_draw_obol_view_feature_points_replace(
    struct ged_view_context *view_ctx,
    ged_view_edit_ref ref,
    enum ged_view_edit_geometry_family family,
    const point_t *points,
    const int *cmds,
    size_t point_count);
extern int ged_draw_obol_view_feature_edit_preview_replace(
    struct ged_view_context *view_ctx,
    ged_view_edit_ref ref,
    const char *source_path,
    const char *edit_intent_id,
    const char *edit_intent_role,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    uint32_t source_revision,
    uint32_t inputs_revision);
extern int ged_draw_obol_view_feature_clear_geometry(
    struct ged_view_context *view_ctx,
    ged_view_edit_ref ref);
extern int ged_draw_overlay_geometry_insert_context(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	const char *name,
	const struct ged_draw_overlay_geometry *geometry,
	const unsigned char basecolor[3],
	fastf_t transparency,
	int draw_mode,
	int csoltab,
	ged_draw_shape_ref *out);

struct ged_draw_shape_record_summary {
    const struct db_full_path *fullpath;
    const char *display_name;
    const char *leaf_name;
    ged_draw_group_ref owning_group_ref;
    unsigned long long path_hash;
    uint64_t source_revision;
    uint64_t inputs_revision;
    uint64_t realized_source_revision;
    uint64_t realized_inputs_revision;
    int stale;
    const char *stale_reason;
    int visible;
    int selected;
    int highlighted;
    int evaluated_region;
    uint64_t drawn_revision;
    fastf_t transparency;
    int draw_mode;
    int line_width;
    point_t center;
};

struct ged_draw_group_record_summary {
    const char *path;
    struct ged_view_context *view_ctx;
    int draw_mode;
    fastf_t transparency;
    int visible;
    int is_overlay;
    int is_view_scope;
    int in_view_scope;
    int is_view_source;
    int is_local_source;
};

enum ged_draw_view_export_query_flag {
    GED_DRAW_VIEW_EXPORT_QUERY_VISIBLE_ONLY = 0x0001u,
    GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS = 0x0002u,
    GED_DRAW_VIEW_EXPORT_QUERY_VIEW_OBJECTS = 0x0004u,
    GED_DRAW_VIEW_EXPORT_QUERY_LOCAL_ONLY = 0x0008u
};

enum ged_draw_view_export_render_flag {
    GED_DRAW_VIEW_EXPORT_RENDER_VISIBLE_ONLY = 0x01u,
    GED_DRAW_VIEW_EXPORT_RENDER_PAYLOAD_PREPARE = 0x02u
};

#define GED_DRAW_VIEW_EXPORT_DRAW_MODE_ANY (-1)

enum ged_draw_view_export_geometry_kind {
    GED_DRAW_VIEW_EXPORT_GEOMETRY_NONE = 0,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_MESH,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_LINE_SET,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_POINT_SET,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_FACE_SET,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_TEXT,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_IMAGE,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_OVERLAY,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_HUD,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_ANNOTATION,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_CSG_PROXY,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_BREP_PROXY,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_EDIT_PREVIEW,
    GED_DRAW_VIEW_EXPORT_GEOMETRY_EXTERNAL_PROXY
};

#define GED_DRAW_VIEW_LINE_POINT_DRAW 12

#define GED_DRAW_OBOL_ANNOTATION_SEGMENT_NONE 0
#define GED_DRAW_OBOL_ANNOTATION_SEGMENT_LINE 1
#define GED_DRAW_OBOL_ANNOTATION_SEGMENT_TEXT 2

struct ged_draw_obol_annotation_segment {
    int kind;
    int line_start;
    int line_end;
    int text_ref_point;
    const char *text;
};

struct ged_draw_obol_annotation_summary {
    int valid;
    size_t point_count;
    size_t segment_count;
    int line_segment_valid;
    int line_start;
    int line_end;
    int text_segment_valid;
    int text_ref_point;
    const char *text;
    uint64_t cache_identity;
    uint64_t source_identity;
};

enum ged_draw_obol_view_feature_kind {
    GED_DRAW_OBOL_VIEW_FEATURE_KIND_UNKNOWN = 0,
    GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINES,
    GED_DRAW_OBOL_VIEW_FEATURE_KIND_INDEXED_LINES,
    GED_DRAW_OBOL_VIEW_FEATURE_KIND_POINTS,
    GED_DRAW_OBOL_VIEW_FEATURE_KIND_LABELS,
    GED_DRAW_OBOL_VIEW_FEATURE_KIND_ARROW,
    GED_DRAW_OBOL_VIEW_FEATURE_KIND_AXES,
    GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINE_LAYER,
    GED_DRAW_OBOL_VIEW_FEATURE_KIND_EDIT_PREVIEW,
    GED_DRAW_OBOL_VIEW_FEATURE_KIND_INDEXED_FACE_SET,
    GED_DRAW_OBOL_VIEW_FEATURE_KIND_POLYGON_OVERLAY,
    GED_DRAW_OBOL_VIEW_FEATURE_KIND_HUD_LABEL,
    GED_DRAW_OBOL_VIEW_FEATURE_KIND_CUSTOM_NODE
};

struct ged_draw_obol_view_feature_record {
    const char *name;
    enum ged_draw_obol_view_feature_kind kind;
    int local;
    int visible;
    int realized;
    unsigned char color[3];
    int line_style;
    int line_width;
    size_t point_count;
    size_t command_count;
    size_t index_count;
    size_t normal_count;
    size_t label_count;
    size_t axes_center_count;
    size_t child_count;
    const char *line_layer_parent_name;
    size_t line_layer_index;
};

typedef int (*ged_draw_obol_view_feature_record_cb)(
	const struct ged_draw_obol_view_feature_record *rec,
	void *userdata);

struct ged_draw_view_export_detail {
    struct bu_vls path;
    struct ged_draw_view_db_object_record record;
    enum ged_draw_view_export_geometry_kind geometry_kind;
    struct {
	point_t *points;
	int *commands;
	int *indices;
	size_t point_count;
	size_t command_count;
	size_t index_count;
    } arrays;
    struct {
	point_t *points;
	int *indices;
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
    } surface;
    struct {
	int display_space;
	point_t *points;
	size_t point_count;
	size_t segment_count;
	int line_segment_valid;
	int line_start;
	int line_end;
	int text_segment_valid;
	int text_ref_point;
	char *text;
    } annotation;
};

extern void *ged_draw_source_root_attach_view_contexts(
	struct ged *gedp,
	struct ged_view_context *active_view_ctx,
	struct bu_ptbl *views);
extern int ged_draw_obol_view_context_feature_store_active(
	struct ged_view_context *view_ctx);
extern int ged_draw_obol_view_context_hud_axes_replace(
	struct ged_view_context *view_ctx,
	const char *name,
	const struct bv_axes_state *axes,
	const mat_t rotation);
extern int ged_draw_obol_view_context_feature_records_foreach(
	struct ged_view_context *view_ctx,
	unsigned int query_flags,
	const char *glob,
	ged_draw_obol_view_feature_record_cb cb,
	void *userdata);
extern ged_view_polygon_ref
ged_draw_obol_view_context_polygon_find(
	struct ged_view_context *view_ctx,
	const char *name);
extern ged_view_polygon_ref
ged_draw_obol_view_context_polygon_find_scoped(
	struct ged_view_context *view_ctx,
	const char *name,
	int local_only);
extern ged_view_polygon_ref
ged_draw_obol_view_context_polygon_select(
	struct ged_view_context *view_ctx,
	const point_t model_point);
extern ged_view_polygon_ref
ged_draw_obol_view_context_polygon_create(
	struct ged_view_context *view_ctx,
	const char *name,
	int local,
	enum ged_view_polygon_type type,
	const point_t screen_point);
extern ged_view_polygon_ref
ged_draw_obol_view_context_polygon_dup(
	struct ged_view_context *view_ctx,
	const char *name,
	const char *new_name);
extern ged_view_polygon_ref
ged_draw_obol_view_context_polygon_import_sketch(
	const char *name,
	struct db_i *dbip,
	struct directory *dp,
	struct ged_view_context *view_ctx,
	int local);
extern void ged_draw_obol_view_context_polygon_visit_records(
	struct ged_view_context *view_ctx,
	ged_view_polygon_record_cb callback,
	void *data);
extern size_t ged_draw_obol_view_context_polygon_snap_count(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref exclude);
extern int ged_draw_obol_view_context_polygon_clear_point_selection(
	struct ged_view_context *view_ctx);
extern int ged_draw_obol_view_context_polygon_snap_exclude_set(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);
extern struct directory *ged_draw_obol_view_polygon_export_sketch(
	struct ged_view_context *view_ctx,
	struct db_i *dbip,
	const char *name,
	ged_view_polygon_ref ref);
extern int ged_draw_obol_view_polygon_record_get(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref,
	struct ged_view_polygon_record *record);
extern int ged_draw_obol_view_context_polygon_update(
	ged_view_polygon_ref ref,
	struct ged_view_context *view_ctx,
	int op);
extern int ged_draw_obol_view_context_polygon_update_screen_pt(
	ged_view_polygon_ref ref,
	struct ged_view_context *view_ctx,
	int x,
	int y,
	int op);
extern int ged_draw_obol_view_polygon_move(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref,
	const point_t current_point,
	const point_t previous_point);
extern int ged_draw_obol_view_polygon_set_name(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref,
	const char *name);
extern int ged_draw_obol_view_polygon_set_visual(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref,
	const struct bu_color *edge_color,
	const struct bu_color *fill_color,
	fastf_t fill_slope_x,
	fastf_t fill_slope_y,
	fastf_t fill_density,
	fastf_t vZ,
	int fill_flag);
extern int ged_draw_obol_view_context_polygon_set_current(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref,
	long contour_i,
	long point_i);
extern int ged_draw_obol_view_context_polygon_set_contour_open(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref,
	long contour_i,
	int open);
extern int ged_draw_obol_view_polygon_set_all_contours_open(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref,
	int open);
extern int ged_draw_obol_view_polygon_close(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);
extern int ged_draw_obol_view_polygon_clear_selected_point(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);
extern int ged_draw_obol_view_polygon_remove(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);
extern void *ged_draw_obol_view_polygon_user_data(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref);
extern int ged_draw_obol_view_polygon_user_data_set(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref,
	void *user_data);
extern int ged_draw_obol_view_context_polygon_area(
	ged_view_polygon_ref ref,
	struct ged_view_context *view_ctx,
	fastf_t *area);
extern int ged_draw_obol_view_context_polygon_overlap(
	ged_view_polygon_ref ref,
	struct ged_view_context *view_ctx,
	const char *other_name,
	const struct bn_tol *tol,
	int *overlap);
extern int ged_draw_obol_view_polygon_csg(
	struct ged_view_context *view_ctx,
	ged_view_polygon_ref target,
	ged_view_polygon_ref stencil,
	enum bg_polygon_boolean_op op);
extern size_t ged_draw_obol_view_context_label_count(
	struct ged_view_context *view_ctx,
	const char *name);
extern int ged_draw_obol_view_context_label_copy(
	struct ged_view_context *view_ctx,
	const char *name,
	size_t index,
	struct bu_vls *text,
	point_t point,
	unsigned char rgb[3]);
extern int ged_draw_obol_view_context_feature_points_copy(
	struct ged_view_context *view_ctx,
	const char *name,
	point_t **points,
	size_t *point_count);
extern int ged_draw_obol_view_context_feature_indices_copy(
	struct ged_view_context *view_ctx,
	const char *name,
	int **indices,
	size_t *index_count);
extern int ged_draw_obol_view_context_feature_line_command_at(
	struct ged_view_context *view_ctx,
	const char *name,
	size_t index,
	int *out);
extern int ged_draw_obol_view_context_feature_layer_points_copy(
	struct ged_view_context *view_ctx,
	const char *name,
	size_t layer_index,
	point_t **points,
	size_t *point_count);
extern int ged_draw_obol_view_context_feature_layer_line_command_at(
	struct ged_view_context *view_ctx,
	const char *name,
	size_t layer_index,
	size_t point_index,
	int *out);
extern int ged_draw_obol_overlay_erase_name_context(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	const char *name);
extern int ged_draw_obol_overlay_geometry_insert_context(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	const char *name,
	const struct ged_draw_overlay_geometry *geometry,
	const unsigned char basecolor[3],
	fastf_t transparency,
	int ged_draw_mode,
	char **shape_path_out);
extern int ged_draw_source_root_foreach_shape_ref(
	struct ged *gedp,
	int skip_overlay_groups,
	ged_draw_shape_ref_index_cb cb,
	void *userdata);
extern int ged_draw_source_root_foreach_group_ref(
	struct ged *gedp,
	ged_draw_group_ref_index_cb cb,
	void *userdata);
extern int ged_draw_source_root_subtree_bounds(
	struct ged *gedp,
	vect_t *min,
	vect_t *max,
	int include_overlays);
extern int ged_draw_source_root_has_groups(struct ged *gedp);
extern int ged_draw_source_root_clear_all_scope_children(
	struct ged *gedp);
extern int ged_draw_source_erase_path_at_root(
	struct ged *gedp,
	const char *path);
extern int ged_draw_source_erase_path_prefix_at_root(
	struct ged *gedp,
	const char *path);
extern int ged_draw_source_erase_groups_by_name_at_root(
	struct ged *gedp,
	const char *name);
extern int ged_draw_source_erase_path_in_active_scope(
	struct ged *gedp,
	const char *path,
	struct ged_view_context *view_ctx,
	int mode);
extern int ged_draw_source_erase_path_prefix_in_active_scope(
	struct ged *gedp,
	const char *path,
	struct ged_view_context *view_ctx,
	int mode);
extern int ged_draw_source_erase_component_name_in_active_scope(
	struct ged *gedp,
	const char *name,
	struct ged_view_context *view_ctx,
	int mode,
	int nonroot_only);
extern int ged_draw_source_clear_db_groups_in_scope(
	struct ged *gedp,
	struct ged_view_context *view_ctx);
extern void ged_draw_log(int level, const char *fmt, ...) _BU_ATTR_PRINTF23;
extern void ged_draw_test_force_primitive_face_set_failure(int enable);
extern int ged_draw_test_primitive_face_set_failure_enabled(void);
extern void ged_draw_view_context_foreach_export_record(
	struct ged_view_context *view_ctx,
	unsigned int query_flags,
	unsigned int render_flags,
	const char *glob,
	int draw_mode,
	ged_draw_view_db_object_record_cb cb,
	void *userdata);

/* Private typed-NMG draw style bits.  The values match the historical plot
 * flags so command option plumbing can stay simple while the API names stay
 * out of the legacy vlist/plot namespace. */
enum ged_draw_nmg_style {
    GED_DRAW_NMG_STYLE_VECTOR = 0,
    GED_DRAW_NMG_STYLE_POLYGON = 1,
    GED_DRAW_NMG_STYLE_VISUALIZE_NORMALS = 2,
    GED_DRAW_NMG_STYLE_USE_VU_NORMALS = 4,
    GED_DRAW_NMG_STYLE_NO_SURFACES = 8
};

extern int ged_draw_shape_ref_refresh_material_color(struct ged *gedp,
								ged_draw_shape_ref ref,
								struct db_i *dbip,
								uint64_t mater_rev);
extern int ged_draw_registry_shape_ref_set_indexed_fullpath(
	struct ged *gedp,
	ged_draw_shape_ref shape_ref,
	const struct db_full_path *path);
extern int ged_draw_registry_group_ref_set_indexed_fullpath(
	struct ged *gedp,
	ged_draw_group_ref group_ref,
	const struct db_full_path *path);
extern const char *ged_draw_registry_shape_ref_semantic_path(
	struct ged *gedp,
	ged_draw_shape_ref shape_ref);
extern const char *ged_draw_registry_shape_ref_cached_semantic_path(
	struct ged *gedp,
	ged_draw_shape_ref shape_ref);
extern char *ged_draw_registry_shape_ref_cached_fullpath_dup(
	struct ged *gedp,
	ged_draw_shape_ref shape_ref);
extern const char *ged_draw_registry_group_ref_semantic_path(
	struct ged *gedp,
	ged_draw_group_ref group_ref);
extern int ged_draw_registry_shape_ref_apply_region(
	struct ged *gedp,
	ged_draw_shape_ref shape_ref,
	int region_id,
	int aircode,
	int los,
	int material_id);
extern const char *ged_draw_dbpath_skip_lead_slash(const char *s);

extern struct ged_view_context *ged_draw_active_view_ctx(struct ged *gedp);
extern void ged_draw_active_view_ctx_set(struct ged *gedp,
						    struct ged_view_context *view_ctx);
extern int ged_draw_ensure_root_attached(struct ged *gedp);

extern int ged_draw_overlay_geometry_insert(struct ged *gedp,
						       const char *name,
						       const struct ged_draw_overlay_geometry *geometry,
						       long int rgb,
						       fastf_t transparency,
						       int draw_mode,
						       int csoltab,
						       ged_draw_shape_ref *out);
extern void ged_draw_overlay_erase_name(struct ged *gedp,
						   const char *name);
extern int ged_draw_redraw_group_ref(struct ged *gedp,
						ged_draw_group_ref ref,
						int skip_subtractions);
extern int ged_draw_group_ref_redraw_wireframe(struct ged *gedp,
							  ged_draw_group_ref ref,
							  struct db_i *dbip,
							  const struct bn_tol *tol,
							  const struct bg_tess_tol *ttol,
							  struct ged_view_context *view_ctx,
							  int skip_subtractions);
extern int ged_draw_shape_ref_redraw_wireframe(struct ged *gedp,
							  ged_draw_shape_ref ref,
							  struct db_i *dbip,
							  const struct bn_tol *tol,
							  const struct bg_tess_tol *ttol,
							  struct ged_view_context *view_ctx,
							  int skip_subtractions);
extern ged_draw_group_ref ged_draw_group_ref_lookup_or_create(struct ged *gedp,
									 const struct db_full_path *dfp);

extern int ged_draw_erase_path_string(struct ged *gedp,
						 const char *path);
extern void ged_draw_erase_name(struct ged *gedp, const char *name);
extern int ged_draw_erase_path_prefix_string(struct ged *gedp,
								const char *path);
extern int ged_draw_erase_path_string_scoped(struct ged *gedp,
							const char *path,
							struct ged_view_context *view_ctx,
							int mode);
extern int ged_draw_erase_path_prefix_string_scoped(struct ged *gedp,
							       const char *path,
							       struct ged_view_context *view_ctx,
							       int mode);
extern int ged_draw_erase_nonroot_component_string_scoped(struct ged *gedp,
								     const char *name,
								     struct ged_view_context *view_ctx,
								     int mode);
extern int ged_draw_erase_component_string_scoped(struct ged *gedp,
							     const char *name,
							     struct ged_view_context *view_ctx,
							     int mode);
extern int ged_draw_apply_draw_inner(
	struct ged *gedp, const struct ged_draw_transaction *txn,
	const char *path, struct ged_draw_transaction_result *result);

extern void ged_draw_highlighted_shape_ref_invalidate(struct ged *gedp);

extern ged_draw_scene_handle ged_view_context_scene_root_ref(
    const struct ged_view_context *view);
extern int ged_view_context_scene_root_ref_attach(
    struct ged_view_context *view, ged_draw_scene_handle root_ref);

extern ged_draw_shape_ref ged_draw_registry_shape_ref_from_source_ref(
	struct ged *gedp,
	ged_draw_scene_handle scene_ref);
extern ged_draw_group_ref ged_draw_registry_group_ref_from_source_ref(
	struct ged *gedp,
	ged_draw_scene_handle scene_ref);
extern ged_draw_scene_handle ged_draw_registry_shape_ref_scene_handle(
	struct ged *gedp,
	ged_draw_shape_ref ref);
extern ged_draw_scene_handle ged_draw_registry_shape_ref_cache_scene_handle(
	struct ged *gedp,
	ged_draw_shape_ref ref);
extern ged_draw_scene_handle ged_draw_registry_group_ref_scene_handle(
	struct ged *gedp,
	ged_draw_group_ref ref);
extern int ged_draw_shape_ref_record_summary(struct ged *gedp,
							ged_draw_shape_ref ref,
							struct ged_draw_shape_record_summary *out);
extern int ged_draw_shape_ref_owning_group_record_summary(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	struct ged_draw_group_record_summary *out);
extern int ged_draw_shape_ref_database_source_summary(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	struct ged_draw_database_source_summary *out);
extern int ged_draw_obol_shape_ref_index_for_path_hash(
	struct ged *gedp,
	unsigned long long path_hash,
	ged_draw_shape_ref_index_cb cb,
	void *userdata);
extern int ged_draw_obol_scene_database_bounds(
	struct ged *gedp,
	vect_t *min,
	vect_t *max,
	int *empty_out);
extern int ged_draw_obol_scene_database_autoview_bounds(
	struct ged *gedp,
	vect_t *min,
	vect_t *max,
	int *empty_out,
	int allow_member_bounds);
extern int ged_database_path_member_autoview_bounds(
	struct ged *gedp,
	const char *path,
	vect_t *min,
	vect_t *max);
extern int ged_draw_obol_progressive_autoview_follow(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	fastf_t factor);
extern int ged_draw_obol_scene_root_child_count(
	struct ged *gedp,
	size_t *out);
extern void ged_draw_obol_context_tokens_free(struct ged *gedp);
extern void ged_draw_obol_scene_context_info_free(
	struct ged_draw_obol_scene_context_info *info);
extern int ged_draw_obol_scene_context_info_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_obol_scene_context_info *out);
extern int ged_draw_obol_scene_child_context_info_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	struct ged_draw_obol_scene_context_info *out);
extern int ged_draw_obol_database_source_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_database_source_summary *out);
extern int ged_draw_obol_database_source_paths_foreach(
	struct ged *gedp,
	int skip_overlay_groups,
	ged_draw_obol_database_source_path_cb cb,
	void *userdata);
extern int ged_draw_obol_database_source_records_foreach(
	struct ged *gedp,
	int skip_overlay_groups,
	ged_draw_obol_database_source_record_cb cb,
	void *userdata);
extern int ged_draw_obol_visible_database_source_records_foreach_fast(
	struct ged *gedp,
	ged_draw_obol_database_source_record_cb cb,
	void *userdata);
extern int ged_draw_obol_shape_paths_foreach(
	struct ged *gedp,
	int skip_overlay_groups,
	ged_draw_obol_database_source_path_cb cb,
	void *userdata);
extern int ged_draw_obol_group_database_source_paths_foreach(
	struct ged *gedp,
	const char *group_path,
	ged_draw_obol_database_source_path_cb cb,
	void *userdata);
extern int ged_draw_obol_group_database_source_records_foreach(
	struct ged *gedp,
	const char *group_path,
	ged_draw_obol_database_source_record_cb cb,
	void *userdata);
extern int ged_draw_obol_database_source_owner_group_path_for_path(
	struct ged *gedp,
	const char *path,
	struct bu_vls *out);
extern int ged_draw_obol_database_source_display_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_scene_display_summary *out);
extern int ged_draw_obol_database_source_bounds_for_path(
	struct ged *gedp,
	const char *path,
	vect_t *min,
	vect_t *max,
	int include_overlays,
	int *empty_out);
extern int ged_draw_obol_group_display_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_scene_display_summary *out);
extern int ged_draw_obol_group_subtree_bounds_for_path(
	struct ged *gedp,
	const char *path,
	vect_t *min,
	vect_t *max,
	int include_overlays,
	int *empty_out);
extern int ged_draw_obol_database_source_last_point_for_path(
	struct ged *gedp,
	const char *path,
	point_t out);
extern int ged_draw_obol_database_source_line_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_view_line_summary *out);
extern int ged_draw_obol_database_source_line_point_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	point_t out);
extern int ged_draw_obol_database_source_line_command_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	int *out);
extern int ged_draw_obol_database_source_line_data_copy_for_path(
	struct ged *gedp,
	const char *path,
	point_t **points,
	int **commands,
	size_t *point_count);
extern int ged_draw_obol_database_source_surface_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_view_surface_summary *out);
extern int ged_draw_obol_database_source_surface_point_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	point_t out);
extern int ged_draw_obol_database_source_surface_index_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	int *out);
extern int ged_draw_obol_shape_display_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_scene_display_summary *out);
extern int ged_draw_obol_shape_geometry_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_shape_geometry_summary *out);
extern int ged_draw_obol_shape_surface_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_view_surface_summary *out);
extern int ged_draw_obol_shape_surface_point_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	point_t out);
extern int ged_draw_obol_shape_surface_index_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	int *out);
extern int ged_draw_obol_shape_line_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_view_line_summary *out);
extern int ged_draw_obol_shape_line_point_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	point_t out);
extern int ged_draw_obol_shape_line_command_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	int *out);
extern int ged_draw_obol_database_source_translate_vlist_for_path(
	struct ged *gedp,
	const char *path,
	const vect_t xlate);
extern int ged_draw_obol_database_source_clear_vlist_for_path(
	struct ged *gedp,
	const char *path);
extern int
ged_draw_obol_database_source_publish_annotation_line_set_for_path(
	struct ged *gedp,
	const char *path,
	const point_t *points,
	const int *commands,
	size_t point_count);
extern int
ged_draw_obol_database_source_publish_annotation_record_for_path(
	struct ged *gedp,
	const char *path,
	const point_t base_point,
	const point_t *annotation_points,
	size_t annotation_point_count,
	const struct ged_draw_obol_annotation_segment *segments,
	size_t segment_count,
	const point_t *line_points,
	const int *line_commands,
	size_t line_point_count);
extern int
ged_draw_obol_database_source_annotation_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_obol_annotation_summary *out);
extern int
ged_draw_obol_database_source_annotation_point_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	point_t out);
extern int
ged_draw_obol_database_source_publish_auxiliary_line_set_for_path(
	struct ged *gedp,
	const char *path,
	const char *name,
	const point_t *points,
	const int *commands,
	size_t point_count,
	const struct ged_draw_scene_display_summary *display_state);
extern int
ged_draw_obol_database_source_publish_auxiliary_source_line_set_for_path(
	struct ged *gedp,
	const char *path,
	const char *source_path,
	const char *display_name,
	const point_t *points,
	const int *commands,
	size_t point_count,
	const struct ged_draw_scene_display_summary *display_state);
extern int
ged_draw_obol_database_source_clear_auxiliary_shapes_for_path(
	struct ged *gedp,
	const char *path);
extern int ged_draw_obol_database_source_clear_mesh_for_path(
	struct ged *gedp,
	const char *path);
extern int ged_draw_obol_local_shape_publish_line_set_for_path(
	struct ged *gedp,
	const char *group_path,
	const char *shape_path,
	const char *display_name,
	const point_t *points,
	const int *commands,
	size_t point_count,
	const struct ged_draw_scene_display_summary *display_state);
extern int ged_draw_obol_local_shape_set_record_role_for_path(
	struct ged *gedp,
	const char *shape_path,
	const char *record_role);
extern int
ged_draw_obol_local_shape_publish_indexed_face_set_for_path(
	struct ged *gedp,
	const char *group_path,
	const char *shape_path,
	const char *display_name,
	const point_t *points,
	size_t point_count,
	const vect_t *normals,
	size_t normal_count,
	const int *indices,
	size_t index_count,
	const struct ged_draw_scene_display_summary *display_state);
extern int ged_draw_obol_database_source_set_vlist_center_for_path(
	struct ged *gedp,
	const char *path,
	const point_t center);
extern int
ged_draw_obol_database_source_update_vlist_bounds_for_path(
	struct ged *gedp,
	const char *path);
extern int ged_draw_obol_database_source_set_placement_for_path(
	struct ged *gedp,
	const char *path,
	int draw_mat_valid,
	const mat_t draw_mat,
	int draw_center_valid,
	const point_t draw_center,
	int draw_size_valid,
	fastf_t draw_size);
extern int ged_draw_obol_database_source_set_bounds_for_path(
	struct ged *gedp,
	const char *path,
	const point_t bmin,
	const point_t bmax);
extern int ged_draw_obol_database_source_set_display_name_for_path(
	struct ged *gedp,
	const char *path,
	const char *display_name);
extern int ged_draw_obol_database_source_geometry_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_shape_geometry_summary *out);
extern int ged_draw_obol_database_source_geometry_summary_for_path_mode(
	struct ged *gedp,
	const char *path,
	int draw_mode_valid,
	int ged_draw_mode,
	struct ged_draw_shape_geometry_summary *out);
extern int ged_draw_obol_database_source_material_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_shape_material_summary *out);
extern int ged_draw_obol_database_source_refresh_material_color_for_path(
	struct ged *gedp,
	const char *path,
	struct db_i *dbip,
	uint64_t material_revision);
extern int ged_draw_obol_database_source_evaluated_region_for_path(
	struct ged *gedp,
	const char *path,
	int *out);
extern int ged_draw_obol_database_source_set_evaluated_region_for_path(
	struct ged *gedp,
	const char *path,
	int evaluated_region);
extern int ged_draw_obol_database_source_set_region_metadata_for_path(
	struct ged *gedp,
	const char *path,
	int region_id,
	int aircode,
	int los,
	int material_id);
extern int ged_draw_obol_database_source_apply_draw_metadata_for_path(
	struct ged *gedp,
	const char *path,
	const struct BObolDrawMetadataRecord *record);
extern int ged_draw_obol_database_source_update_display_for_path(
	struct ged *gedp,
	const char *path,
	int visible_valid,
	int visible,
	int selected_valid,
	int selected,
	int highlighted_valid,
	int highlighted,
	int draw_mode_valid,
	int ged_draw_mode,
	int line_style_valid,
	int line_style,
	int line_width_valid,
	int line_width,
	int transparency_valid,
	fastf_t transparency,
	int color_valid,
	const unsigned char color[3],
	int material_color_valid,
	const unsigned char material_color[3],
	int material_revision_valid,
	uint64_t material_revision);
extern int ged_draw_obol_database_source_set_selected_for_instance_key(
	struct ged *gedp,
	const char *instance_key,
	int selected);
extern int ged_draw_obol_database_sources_sync_selected_paths(
    struct ged *gedp,
    const char *const *paths,
    size_t path_count);
extern int ged_draw_obol_database_sources_apply_selection_delta(
    struct ged *gedp,
    const char *const *added_paths,
    size_t added_count,
    const char *const *removed_paths,
    size_t removed_count,
    const char *const *selected_paths,
    size_t selected_count);
extern int ged_draw_obol_shape_update_display_for_path(
	struct ged *gedp,
	const char *path,
	int visible_valid,
	int visible,
	int selected_valid,
	int selected,
	int highlighted_valid,
	int highlighted,
	int draw_mode_valid,
	int ged_draw_mode,
	int line_style_valid,
	int line_style,
	int line_width_valid,
	int line_width,
	int transparency_valid,
	fastf_t transparency,
	int color_valid,
	const unsigned char color[3],
	int material_color_valid,
	const unsigned char material_color[3],
	int material_revision_valid,
	uint64_t material_revision);
extern int ged_draw_obol_database_source_update_appearance_for_path(
	struct ged *gedp,
	const char *path,
	int line_width_valid,
	int line_width,
	int transparency_valid,
	fastf_t transparency,
	int color_override_valid,
	int color_override,
	int color_valid,
	const unsigned char color[3]);
extern int ged_draw_obol_database_source_mark_stale_for_path(
	struct ged *gedp,
	const char *path,
	ged_draw_stale_reason reason);
extern int ged_draw_obol_database_source_record_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_obol_database_source_record *out);
extern int ged_draw_obol_database_source_runtime_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_obol_database_source_runtime *out);
extern int ged_draw_obol_database_source_set_mesh_lod_for_path(
	struct ged *gedp,
	const char *path,
	struct BObolMeshLod *lod);
extern int ged_draw_obol_database_source_set_mesh_lod_bounds_for_path(
	struct ged *gedp,
	const char *path,
	const point_t bmin,
	const point_t bmax);
extern int ged_draw_obol_database_source_apply_record_for_path(
	struct ged *gedp,
	const char *path,
	const struct ged_draw_obol_database_source_record *record);
extern int ged_draw_obol_database_source_draw_state_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_obol_draw_state_summary *out);
extern int
ged_draw_obol_database_source_set_realization_for_path(
	struct ged *gedp,
	const char *path,
	int current,
	int failed,
	uint64_t realized_source_revision,
	uint64_t realized_inputs_revision,
	ged_draw_stale_reason reason);
extern int
ged_draw_obol_database_source_realization_policy_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_obol_realization_policy_summary *out);
extern int
ged_draw_obol_database_source_set_realization_roles_for_path(
	struct ged *gedp,
	const char *path,
	int role_flags);
extern int
ged_draw_obol_database_source_set_realization_view_policy_for_path(
	struct ged *gedp,
	const char *path,
	int view_dependent,
	int csg_lod_enabled,
	int mesh_lod_enabled,
	fastf_t view_scale,
	fastf_t lod_scale,
	int view_width,
	int view_height,
	uint64_t bot_threshold,
	fastf_t curve_scale,
	fastf_t point_scale);
extern int
ged_draw_obol_database_source_realize_for_path(struct ged *gedp,
	const char *path);
extern int
ged_draw_obol_database_source_realize_pending(struct ged *gedp);

int
ged_draw_obol_database_sources_redraw(struct ged *gedp,
	struct ged_view_context *view_ctx,
	const char *path,
	int draw_mode);
extern int ged_draw_obol_controller_attach_opaque_for_view(
	struct ged *gedp, struct ged_view_context *view_ctx, void *controller,
	int sync_current_scene);
extern int ged_draw_obol_render_endpoint_ensure_for_view(
	struct ged *gedp, struct ged_view_context *view_ctx, int sync_current_scene);
extern void ged_draw_obol_controller_detach_for_view(
	struct ged *gedp, struct ged_view_context *view_ctx);
extern int ged_draw_obol_scene_controller_ensure_owned(
	struct ged *gedp,
	int sync_current_scene);
extern int ged_draw_obol_scene_controller_attached(
	struct ged *gedp);
extern int ged_draw_obol_scene_controller_full_synced(
    struct ged *gedp);
extern int ged_draw_obol_scene_controller_owned_internal(
    struct ged *gedp);
extern int ged_draw_obol_scene_sync_attached_transaction(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result);
extern int ged_draw_obol_backend_apply_transaction(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result);
extern uint64_t ged_scene_revision_advance(struct ged *gedp);
extern void ged_scene_delta_dispatch_internal(
    struct ged *gedp,
    enum ged_scene_delta_kind kind,
    struct ged_view_context *view_ctx,
    const char *const *paths,
    size_t path_count,
    int affects_all,
    size_t group_count,
    size_t shape_count,
    int presentation_only,
    uint64_t revision_before,
    uint64_t revision_after,
    const char *diagnostic);
extern int ged_draw_obol_database_source_ensure_for_path(
	struct ged *gedp,
	const char *path,
	struct db_i *dbip,
	int ged_draw_mode,
	uint64_t source_revision);
extern int ged_draw_obol_database_source_ensure_for_path_with_placement(
	struct ged *gedp,
	const char *path,
	struct db_i *dbip,
	int ged_draw_mode,
	uint64_t source_revision,
	int draw_mat_valid,
	const mat_t draw_mat,
	int draw_center_valid,
	const point_t draw_center,
	int draw_size_valid,
	fastf_t draw_size);
extern int ged_draw_obol_database_source_rename_for_path(
	struct ged *gedp,
	const char *path,
	const char *new_path,
	uint64_t source_revision);
extern int ged_draw_obol_database_source_move_to_group_for_path(
	struct ged *gedp,
	const char *source_path,
	const char *group_path);
extern int ged_draw_obol_database_source_remove_for_path(
	struct ged *gedp,
	const char *path);
extern int
ged_draw_obol_database_sources_remove_for_path_prefix(
	struct ged *gedp,
	const char *path_prefix);
extern int
ged_draw_obol_database_sources_remove_for_path_prefix_in_scope(
	struct ged *gedp,
	const char *path_prefix,
	struct ged_view_context *view_ctx);
extern int
ged_draw_obol_database_sources_remove_for_path_prefix_in_scope_mode(
	struct ged *gedp,
	const char *path_prefix,
	struct ged_view_context *view_ctx,
	int mode);
extern int
ged_draw_obol_active_database_sources_remove_for_path_prefix(
	struct ged *gedp,
	const char *path_prefix);
extern int
ged_draw_obol_active_database_sources_remove_for_path_prefix_mode(
	struct ged *gedp,
	const char *path_prefix,
	int mode);
extern int
ged_draw_obol_database_sources_remove_for_component_name(
	struct ged *gedp,
	const char *name,
	int nonroot_only,
	int mode);
extern int ged_draw_obol_database_sources_clear(struct ged *gedp);
extern int ged_draw_obol_database_sources_clear_in_scope(
	struct ged *gedp,
	struct ged_view_context *view_ctx);
extern int ged_draw_obol_scene_clear(struct ged *gedp);
extern int ged_draw_obol_active_scene_clear(struct ged *gedp);
extern int ged_draw_obol_groups_remove_for_component_name(
	struct ged *gedp,
	const char *name);
extern int ged_draw_obol_groups_remove_for_path_prefix(
	struct ged *gedp,
	const char *path_prefix);
extern int ged_draw_obol_group_remove_for_path(
	struct ged *gedp,
	const char *path);
extern int ged_draw_obol_group_clear_for_path(
	struct ged *gedp,
	const char *path);
extern int ged_draw_obol_group_rename_for_path(
	struct ged *gedp,
	const char *path,
	const char *new_path);
int ged_draw_group_erase_subpath_for_path(struct ged *gedp,
	const char *parent_path, const char *subpath);
extern int ged_draw_obol_group_update_display_for_path(
	struct ged *gedp,
	const char *path,
	int visible_valid,
	int visible);
extern int ged_draw_obol_group_ensure_for_path(
	struct ged *gedp,
	const char *path,
	const char *intent_path,
	int mode,
	int overlay);
extern int ged_draw_obol_group_record_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_group_record_summary *out);
extern int ged_draw_obol_group_paths_foreach(
	struct ged *gedp,
	int skip_overlay_groups,
	ged_draw_obol_group_path_cb cb,
	void *userdata);
extern int ged_draw_obol_group_shape_count_for_path(
	struct ged *gedp,
	const char *path,
	int *out);
extern int ged_draw_obol_group_descendant_group_count_for_path(
	struct ged *gedp,
	const char *path,
	size_t *out);
extern int ged_draw_obol_group_child_count_for_path(
	struct ged *gedp,
	const char *path,
	size_t *out);
extern int ged_draw_obol_group_update_appearance_for_path(
	struct ged *gedp,
	const char *path,
	const struct ged_draw_appearance_settings *settings);
extern int ged_draw_obol_group_appearance_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_appearance_settings *settings);
extern int ged_draw_obol_group_update_draw_intent_for_path(
	struct ged *gedp,
	const char *path,
	const char *intent_path,
	int mode_valid,
	int mode,
	int overlay_valid,
	int overlay);
extern int ged_draw_group_ref_record_summary(struct ged *gedp,
							 ged_draw_group_ref ref,
							 struct ged_draw_group_record_summary *out);
extern int ged_draw_group_ref_tree_summary(struct ged *gedp,
						      ged_draw_group_ref ref,
						      struct ged_draw_scene_tree_summary *out);
extern int ged_draw_group_ref_shape_count(struct ged *gedp,
						     ged_draw_group_ref ref,
						     int *out);

extern int ged_draw_group_ref_index_for_component(struct ged *gedp,
							    const char *path,
							    ged_draw_group_ref_index_cb cb,
							    void *userdata);
extern int ged_draw_obol_shape_ref_index_for_component(
	struct ged *gedp,
	const char *path,
	ged_draw_shape_ref_index_cb cb,
	void *userdata);
extern int ged_draw_obol_group_ref_index_for_component(
	struct ged *gedp,
	const char *path,
	ged_draw_group_ref_index_cb cb,
	void *userdata);
struct ged_draw_index_stats {
    uint64_t shape_component_queries;
    uint64_t shape_component_candidates;
    uint64_t group_component_queries;
    uint64_t group_component_candidates;
    uint64_t path_queries;
    uint64_t path_candidates;
    uint64_t slow_path_shape_scans;
    uint64_t slow_path_group_scans;
};

extern void ged_draw_index_stats_get(struct ged *gedp,
						struct ged_draw_index_stats *stats);
extern void ged_draw_index_stats_reset(struct ged *gedp);

extern int ged_draw_group_ref_set_dbpath(struct ged *gedp,
						    ged_draw_group_ref ref,
						    const struct db_full_path *new_dfp);
extern int ged_draw_group_ref_set_mode(struct ged *gedp,
						  ged_draw_group_ref ref,
						  int mode);
extern int ged_draw_group_ref_set_appearance_settings(struct ged *gedp,
								ged_draw_group_ref ref,
								const struct ged_draw_appearance_settings *settings);
extern int ged_draw_group_ref_appearance_settings(struct ged *gedp,
							     ged_draw_group_ref ref,
							     struct ged_draw_appearance_settings *settings);
extern int ged_draw_group_ref_set_visible(struct ged *gedp,
						     ged_draw_group_ref ref,
						     int visible);
extern int ged_draw_shape_ref_set_visible(struct ged *gedp,
	ged_draw_shape_ref ref, int visible);
extern int ged_draw_shape_ref_set_color(struct ged *gedp,
	ged_draw_shape_ref ref, const unsigned char rgb[3]);
extern int ged_draw_shape_ref_set_material_color(struct ged *gedp,
	ged_draw_shape_ref ref, const unsigned char rgb[3]);
extern int ged_draw_shape_ref_set_evaluated_region(struct ged *gedp,
	ged_draw_shape_ref ref, int evaluated_region);
extern int ged_draw_shape_ref_set_highlighted(struct ged *gedp,
							 ged_draw_shape_ref ref,
							 int highlighted);
extern int ged_draw_shape_ref_set_selected(struct ged *gedp,
						      ged_draw_shape_ref ref,
						      int selected);
extern int ged_draw_shape_ref_set_transparency(struct ged *gedp,
							  ged_draw_shape_ref ref,
							  fastf_t transparency);
extern int ged_draw_shape_ref_mark_database_source_changed(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	ged_draw_stale_reason reason);
extern int ged_draw_shape_ref_lod_ensure(struct ged *gedp,
						    ged_draw_shape_ref ref,
						    struct ged_view_context *first_view_ctx,
						    struct ged_view_context **view_ctxs,
						    size_t view_ctx_count);
__END_DECLS

#ifdef __cplusplus
class BObolViewController;
class BObolSceneController;

/* Private C++ integration seam.  Scene ownership is per GED instance; view
 * controller identity is defined solely by the view's display endpoint.  The
 * returned pointers are borrowed and remain valid only while their owner is
 * attached. */
BObolSceneController *ged_bobol_scene(struct ged *gedp);
BObolViewController *ged_bobol_view_controller(
    struct ged_view_context *view_ctx);

void ged_draw_obol_scene_controller_detach(struct ged *gedp);
BObolSceneController *ged_draw_obol_scene_controller(
    struct ged *gedp);
int ged_draw_obol_scene_controller_owned(struct ged *gedp);
int ged_draw_obol_scene_sync_transaction(
    struct ged *gedp, const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    BObolSceneController *controller = NULL);
int ged_draw_obol_scene_sync_full_scene(
    struct ged *gedp, struct ged_view_context *view_ctx, uint32_t source_revision = 0,
    BObolSceneController *controller = NULL);
BObolViewController *ged_draw_obol_controller(struct ged *gedp);
#endif

#endif /* LIBGED_GED_DRAW_PRIVATE_H */
