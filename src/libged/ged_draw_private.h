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

#include "brlobol/mesh_lod_cache.h"
#include "ged/draw.h"
#include "ged/draw_obol.h"
#include "rt/db_fullpath.h"
#include "rt/view.h"

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
struct BRLObolDrawMetadataRecord;
struct ged_draw_obol_database_source_record;

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
    struct BRLObolMeshLod *mesh_lod;
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

GED_EXPORT extern void ged_draw_overlay_erase_name_context(struct ged *gedp,
							   void *view_ctx,
							   const char *name);
GED_EXPORT extern size_t ged_draw_obol_view_context_clear(
	void *view_ctx,
	int flags);
GED_EXPORT extern int ged_draw_obol_view_context_selection_available(
	void *view_ctx);
GED_EXPORT extern size_t ged_draw_obol_view_context_selection_count(
	void *view_ctx);
GED_EXPORT extern int ged_draw_obol_view_context_selection_path_foreach(
	void *view_ctx,
	ged_draw_view_context_selection_path_cb cb,
	void *data);
GED_EXPORT extern int ged_draw_obol_view_context_selection_clear(
	void *view_ctx);
GED_EXPORT extern int ged_draw_obol_view_context_selection_contains_path(
	void *view_ctx,
	int kind,
	const char *path);
GED_EXPORT extern int ged_draw_obol_view_context_selection_add_path(
	void *view_ctx,
	int kind,
	const char *path);
GED_EXPORT extern int ged_draw_obol_view_context_selection_set_path(
	void *view_ctx,
	int kind,
	const char *path);
GED_EXPORT extern struct ged_draw_pick_result *ged_draw_obol_view_context_pick_point(
	void *view_ctx,
	int x,
	int y,
	int first_only);
GED_EXPORT extern struct ged_draw_pick_result *ged_draw_obol_view_context_pick_nearest(
	void *view_ctx,
	int x,
	int y);
GED_EXPORT extern struct ged_draw_pick_result *ged_draw_obol_view_context_pick_rect(
	void *view_ctx,
	int x0,
	int y0,
	int x1,
	int y1);
GED_EXPORT extern int ged_draw_obol_view_context_snap_first_candidate(
	void *view_ctx,
	const point_t sample,
	enum ged_draw_view_snap_kind kind,
	point_t candidate);
GED_EXPORT extern int ged_draw_obol_view_feature_ref_is_null(
	ged_draw_view_feature_ref ref);
GED_EXPORT extern int ged_draw_obol_view_context_edit_preview_publish_event(
	void *view_ctx,
	ged_draw_view_feature_ref feature,
	enum ged_draw_view_edit_preview_event event,
	const char *source_path);
GED_EXPORT extern ged_draw_view_feature_ref
ged_draw_obol_view_context_feature_overlay_ensure(
	void *view_ctx,
	const char *name,
	const void *owner,
	const char *source_path);
GED_EXPORT extern ged_draw_view_feature_ref
ged_draw_obol_view_context_feature_label_ensure(
	void *view_ctx,
	const char *name,
	const void *owner);
GED_EXPORT extern int ged_draw_obol_view_feature_set_context(
	ged_draw_view_feature_ref ref,
	void *view_ctx);
GED_EXPORT extern int ged_draw_obol_view_feature_set_visible(
	ged_draw_view_feature_ref ref,
	int visible);
GED_EXPORT extern int ged_draw_obol_view_feature_set_color(
	ged_draw_view_feature_ref ref,
	int r,
	int g,
	int b);
GED_EXPORT extern int ged_draw_obol_view_feature_touch(
	ged_draw_view_feature_ref ref);
GED_EXPORT extern int ged_draw_obol_view_feature_labels_replace(
	ged_draw_view_feature_ref ref,
	const struct ged_draw_view_feature_label *labels,
	size_t label_count);
GED_EXPORT extern int ged_draw_obol_view_feature_points_replace(
    ged_draw_view_feature_ref ref,
    enum ged_draw_view_feature_family family,
    const point_t *points,
    const int *cmds,
    size_t point_count);
GED_EXPORT extern int ged_draw_obol_view_feature_edit_preview_replace(
    ged_draw_view_feature_ref ref,
    const char *source_path,
    const char *edit_intent_id,
    const char *edit_intent_role,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    uint32_t source_revision,
    uint32_t inputs_revision);
GED_EXPORT extern int ged_draw_obol_view_feature_clear_geometry(
    ged_draw_view_feature_ref ref);
GED_EXPORT extern int ged_draw_overlay_geometry_insert_context(
	struct ged *gedp,
	void *view_ctx,
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
    void *view_ctx;
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

GED_EXPORT extern int ged_draw_brep_mesh_lod_detail_setup(struct BRLObolMeshLod *lod,
							  struct db_i *dbip,
							  struct directory *dp,
							  const struct bg_tess_tol *ttol,
							  const struct bn_tol *tol);
GED_EXPORT extern int ged_draw_brep_mesh_lod_cache_prepare(struct BRLObolMeshLod **lod,
							   point_t bmin,
							   point_t bmax,
							   int *bounds_valid,
							   struct db_i *dbip,
							   struct directory *dp,
							   const struct bg_tess_tol *ttol,
							   const struct bn_tol *tol);
GED_EXPORT extern void *ged_draw_source_root_attach_view_contexts(
	struct ged *gedp,
	void *active_view_ctx,
	struct bu_ptbl *views);
GED_EXPORT extern int ged_draw_obol_view_context_feature_store_active(
	void *view_ctx);
GED_EXPORT extern int ged_draw_obol_view_context_hud_axes_replace(
	void *view_ctx,
	const char *name,
	const struct bv_axes_state *axes,
	const mat_t rotation);
GED_EXPORT extern int ged_draw_obol_view_context_hud_lines_replace(
	void *view_ctx,
	const char *name,
	const point_t *points,
	const int *cmds,
	size_t point_count,
	const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_obol_view_context_hud_labels_replace(
	void *view_ctx,
	const char *name,
	const struct ged_draw_view_label_data *labels,
	size_t label_count,
	const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_obol_view_context_hud_line_layers_replace(
	void *view_ctx,
	const char *name,
	const struct ged_draw_view_line_layer_data *layers,
	size_t layer_count,
	const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_obol_view_context_feature_exists(
	void *view_ctx,
	const char *name);
GED_EXPORT extern int ged_draw_obol_view_context_feature_remove(
	void *view_ctx,
	const char *name);
GED_EXPORT extern size_t ged_draw_obol_view_context_features_remove_prefix(
	void *view_ctx,
	const char *prefix);
GED_EXPORT extern int ged_draw_obol_view_context_feature_summary(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_feature_summary *summary);
GED_EXPORT extern size_t ged_draw_obol_view_context_feature_metadata_count(
	void *view_ctx,
	const char *name);
GED_EXPORT extern int ged_draw_obol_view_context_feature_metadata_copy(
	void *view_ctx,
	const char *name,
	size_t index,
	struct bu_vls *key,
	struct bu_vls *value);
GED_EXPORT extern size_t
ged_draw_obol_view_context_feature_primitive_metadata_count(
	void *view_ctx,
	const char *name,
	int primitive);
GED_EXPORT extern int
ged_draw_obol_view_context_feature_primitive_metadata_copy(
	void *view_ctx,
	const char *name,
	int primitive,
	size_t index,
	struct bu_vls *key,
	struct bu_vls *value);
GED_EXPORT extern int
ged_draw_obol_view_context_feature_pick_primitive_resolve(
	void *view_ctx,
	const char *picked_feature_name,
	int picked_primitive,
	int select,
	int highlight,
	struct bu_vls *feature_name,
	int *feature_primitive);
GED_EXPORT extern int
ged_draw_obol_view_context_feature_selected_primitives_replace(
	void *view_ctx,
	const char *name,
	const int *primitives,
	size_t primitive_count);
GED_EXPORT extern int
ged_draw_obol_view_context_feature_highlighted_primitives_replace(
	void *view_ctx,
	const char *name,
	const int *primitives,
	size_t primitive_count);
GED_EXPORT extern size_t
ged_draw_obol_view_context_feature_selected_primitive_count(
	void *view_ctx,
	const char *name);
GED_EXPORT extern size_t
ged_draw_obol_view_context_feature_highlighted_primitive_count(
	void *view_ctx,
	const char *name);
GED_EXPORT extern int
ged_draw_obol_view_context_feature_selected_primitive_at(
	void *view_ctx,
	const char *name,
	size_t index,
	int *primitive);
GED_EXPORT extern int
ged_draw_obol_view_context_feature_highlighted_primitive_at(
	void *view_ctx,
	const char *name,
	size_t index,
	int *primitive);
GED_EXPORT extern int ged_draw_obol_view_context_feature_visible(
	void *view_ctx,
	const char *name);
GED_EXPORT extern int ged_draw_obol_view_context_feature_visible_set(
	void *view_ctx,
	const char *name,
	int visible);
GED_EXPORT extern int ged_draw_obol_view_context_feature_style_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_obol_view_context_feature_style_apply(
	void *view_ctx,
	const char *name,
	const struct ged_draw_view_feature_style *style,
	int recursive);
GED_EXPORT extern int ged_draw_obol_view_context_feature_realize(
	void *view_ctx,
	const char *name,
	int recursive);
GED_EXPORT extern int ged_draw_obol_view_context_feature_depth(
	void *view_ctx,
	const char *name,
	int mode,
	fastf_t *depth);
GED_EXPORT extern int ged_draw_obol_view_context_feature_depth_foreach(
	void *view_ctx,
	int mode,
	ged_draw_view_feature_depth_cb cb,
	void *data);
GED_EXPORT extern int ged_draw_obol_view_context_feature_records_foreach(
	void *view_ctx,
	unsigned int query_flags,
	const char *glob,
	ged_draw_obol_view_feature_record_cb cb,
	void *userdata);
GED_EXPORT extern int ged_draw_obol_view_context_indexed_face_set_replace(
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
GED_EXPORT extern int ged_draw_obol_view_context_lines_replace(
	void *view_ctx,
	const char *name,
	int local,
	const point_t *points,
	const int *cmds,
	size_t point_count,
	const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_obol_view_context_tcl_polygons_replace(
	void *view_ctx,
	const char *name,
	const point_t *points,
	const int *cmds,
	size_t point_count,
	const struct ged_draw_view_feature_style *style);
GED_EXPORT extern ged_draw_view_polygon_ref
ged_draw_obol_view_context_polygon_find(
	void *view_ctx,
	const char *name);
GED_EXPORT extern ged_draw_view_polygon_ref
ged_draw_obol_view_context_polygon_find_scoped(
	void *view_ctx,
	const char *name,
	int local_only);
GED_EXPORT extern ged_draw_view_polygon_ref
ged_draw_obol_view_context_polygon_select(
	void *view_ctx,
	const point_t model_point);
GED_EXPORT extern ged_draw_view_polygon_ref
ged_draw_obol_view_context_polygon_create(
	void *view_ctx,
	const char *name,
	int local,
	int type,
	const point_t screen_point);
GED_EXPORT extern ged_draw_view_polygon_ref
ged_draw_obol_view_context_polygon_dup(
	void *view_ctx,
	const char *name,
	const char *new_name);
GED_EXPORT extern ged_draw_view_polygon_ref
ged_draw_obol_view_context_polygon_import_sketch(
	const char *name,
	struct db_i *dbip,
	struct directory *dp,
	void *view_ctx,
	int local);
GED_EXPORT extern void ged_draw_obol_view_context_polygon_visit_records(
	void *view_ctx,
	ged_draw_view_polygon_record_cb callback,
	void *data);
GED_EXPORT extern size_t ged_draw_obol_view_context_polygon_snap_count(
	void *view_ctx,
	ged_draw_view_polygon_ref exclude);
GED_EXPORT extern int ged_draw_obol_view_context_polygon_clear_point_selection(
	void *view_ctx);
GED_EXPORT extern int ged_draw_obol_view_context_polygon_snap_exclude_set(
	void *view_ctx,
	ged_draw_view_polygon_ref ref);
GED_EXPORT extern struct directory *ged_draw_obol_view_polygon_export_sketch(
	struct db_i *dbip,
	const char *name,
	ged_draw_view_polygon_ref ref);
GED_EXPORT extern int ged_draw_obol_view_polygon_record_get(
	ged_draw_view_polygon_ref ref,
	struct ged_draw_view_polygon_record *record);
GED_EXPORT extern int ged_draw_obol_view_context_polygon_update(
	ged_draw_view_polygon_ref ref,
	void *view_ctx,
	int op);
GED_EXPORT extern int ged_draw_obol_view_context_polygon_update_screen_pt(
	ged_draw_view_polygon_ref ref,
	void *view_ctx,
	int x,
	int y,
	int op);
GED_EXPORT extern int ged_draw_obol_view_polygon_move(
	ged_draw_view_polygon_ref ref,
	point_t *current_point,
	point_t *previous_point);
GED_EXPORT extern int ged_draw_obol_view_polygon_set_name(
	ged_draw_view_polygon_ref ref,
	const char *name);
GED_EXPORT extern int ged_draw_obol_view_polygon_set_visual(
	ged_draw_view_polygon_ref ref,
	const struct bu_color *edge_color,
	const struct bu_color *fill_color,
	fastf_t fill_slope_x,
	fastf_t fill_slope_y,
	fastf_t fill_density,
	fastf_t vZ,
	int fill_flag);
GED_EXPORT extern int ged_draw_obol_view_context_polygon_set_current(
	ged_draw_view_polygon_ref ref,
	long contour_i,
	long point_i);
GED_EXPORT extern int ged_draw_obol_view_context_polygon_set_contour_open(
	ged_draw_view_polygon_ref ref,
	long contour_i,
	int open);
GED_EXPORT extern int ged_draw_obol_view_polygon_set_all_contours_open(
	ged_draw_view_polygon_ref ref,
	int open);
GED_EXPORT extern int ged_draw_obol_view_polygon_close(
	ged_draw_view_polygon_ref ref);
GED_EXPORT extern int ged_draw_obol_view_polygon_clear_selected_point(
	ged_draw_view_polygon_ref ref);
GED_EXPORT extern int ged_draw_obol_view_polygon_remove(
	ged_draw_view_polygon_ref ref);
GED_EXPORT extern void *ged_draw_obol_view_polygon_user_data(
	ged_draw_view_polygon_ref ref);
GED_EXPORT extern int ged_draw_obol_view_polygon_user_data_set(
	ged_draw_view_polygon_ref ref,
	void *user_data);
GED_EXPORT extern int ged_draw_obol_view_context_polygon_area(
	ged_draw_view_polygon_ref ref,
	void *view_ctx,
	fastf_t *area);
GED_EXPORT extern int ged_draw_obol_view_context_polygon_overlap(
	ged_draw_view_polygon_ref ref,
	void *view_ctx,
	const char *other_name,
	const struct bn_tol *tol,
	int *overlap);
GED_EXPORT extern int ged_draw_obol_view_polygon_csg(
	ged_draw_view_polygon_ref target,
	ged_draw_view_polygon_ref stencil,
	bg_clip_t op);
GED_EXPORT extern int ged_draw_obol_view_context_line_layer_builder_replace(
	void *view_ctx,
	const char *name,
	int local,
	const struct bg_line_layer_builder *builder);
GED_EXPORT extern int ged_draw_obol_view_context_diagnostic_line_layer_builder_replace(
	void *view_ctx,
	const char *name,
	const struct bg_line_layer_builder *builder);
GED_EXPORT extern int ged_draw_obol_view_context_line_layers_replace(
	void *view_ctx,
	const char *name,
	int local,
	const struct ged_draw_view_line_layer_data *layers,
	size_t layer_count,
	const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_obol_view_context_lines_create_model_annotation(
	void *view_ctx,
	const char *name,
	int local,
	const point_t point);
GED_EXPORT extern int ged_draw_obol_view_context_lines_append_point(
	void *view_ctx,
	const char *name,
	const point_t point);
GED_EXPORT extern int ged_draw_obol_view_context_label_create(
	void *view_ctx,
	const char *name,
	int local,
	const char *text,
	const point_t point,
	const point_t target,
	int has_target);
GED_EXPORT extern int ged_draw_obol_view_context_labels_replace(
	void *view_ctx,
	const char *name,
	int local,
	const struct ged_draw_view_label_data *labels,
	size_t label_count);
GED_EXPORT extern int ged_draw_obol_view_context_tcl_labels_replace(
	void *view_ctx,
	const char *name,
	int draw,
	const struct ged_draw_view_label_data *labels,
	size_t label_count);
GED_EXPORT extern size_t ged_draw_obol_view_context_label_count(
	void *view_ctx,
	const char *name);
GED_EXPORT extern int ged_draw_obol_view_context_label_copy(
	void *view_ctx,
	const char *name,
	size_t index,
	struct bu_vls *text,
	point_t point,
	unsigned char rgb[3]);
GED_EXPORT extern int ged_draw_obol_view_context_label_point_set(
	void *view_ctx,
	const char *name,
	size_t index,
	const point_t point);
GED_EXPORT extern int ged_draw_obol_view_context_line_style_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_line_style *style);
GED_EXPORT extern int ged_draw_obol_view_context_line_color_set(
	void *view_ctx,
	const char *name,
	int r,
	int g,
	int b);
GED_EXPORT extern int ged_draw_obol_view_context_line_width_set(
	void *view_ctx,
	const char *name,
	int line_width);
GED_EXPORT extern int ged_draw_obol_view_context_feature_points_copy(
	void *view_ctx,
	const char *name,
	point_t **points,
	size_t *point_count);
GED_EXPORT extern int ged_draw_obol_view_context_feature_indices_copy(
	void *view_ctx,
	const char *name,
	int **indices,
	size_t *index_count);
GED_EXPORT extern int ged_draw_obol_view_context_feature_line_command_at(
	void *view_ctx,
	const char *name,
	size_t index,
	int *out);
GED_EXPORT extern int ged_draw_obol_view_context_feature_layer_points_copy(
	void *view_ctx,
	const char *name,
	size_t layer_index,
	point_t **points,
	size_t *point_count);
GED_EXPORT extern int ged_draw_obol_view_context_feature_layer_line_command_at(
	void *view_ctx,
	const char *name,
	size_t layer_index,
	size_t point_index,
	int *out);
GED_EXPORT extern int ged_draw_obol_view_context_tcl_lines_replace(
	void *view_ctx,
	const char *name,
	const point_t *points,
	size_t point_count,
	const struct ged_draw_view_line_style *style);
GED_EXPORT extern int ged_draw_obol_view_context_arrow_tip_get(
	void *view_ctx,
	const char *name,
	fastf_t *tip_length,
	fastf_t *tip_width);
GED_EXPORT extern int ged_draw_obol_view_context_arrow_tip_set(
	void *view_ctx,
	const char *name,
	fastf_t tip_length,
	fastf_t tip_width);
GED_EXPORT extern int ged_draw_obol_view_context_tcl_arrows_replace(
	void *view_ctx,
	const char *name,
	point_t *points,
	size_t point_count,
	const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_obol_view_context_feature_axes_centers_copy(
	void *view_ctx,
	const char *name,
	point_t **centers,
	size_t *center_count);
GED_EXPORT extern int ged_draw_obol_view_context_tcl_axes_replace(
	void *view_ctx,
	const char *name,
	point_t *centers,
	size_t center_count,
	fastf_t half_axes_size,
	const struct ged_draw_view_feature_style *style);
GED_EXPORT extern int ged_draw_obol_view_context_axes_create(
	void *view_ctx,
	const char *name,
	int local,
	const struct ged_draw_view_axes_state *state);
GED_EXPORT extern int ged_draw_obol_view_context_axes_state_get(
	void *view_ctx,
	const char *name,
	struct ged_draw_view_axes_state *state);
GED_EXPORT extern int ged_draw_obol_view_context_axes_state_replace(
	void *view_ctx,
	const char *name,
	const struct ged_draw_view_axes_state *state);
GED_EXPORT extern int ged_draw_obol_overlay_erase_name_context(
	struct ged *gedp,
	void *view_ctx,
	const char *name);
GED_EXPORT extern int ged_draw_obol_overlay_geometry_insert_context(
	struct ged *gedp,
	void *view_ctx,
	const char *name,
	const struct ged_draw_overlay_geometry *geometry,
	const unsigned char basecolor[3],
	fastf_t transparency,
	int ged_draw_mode,
	char **shape_path_out);
GED_EXPORT extern int ged_draw_source_root_foreach_shape_ref(
	struct ged *gedp,
	int skip_overlay_groups,
	ged_draw_shape_ref_index_cb cb,
	void *userdata);
GED_EXPORT extern int ged_draw_source_root_foreach_group_ref(
	struct ged *gedp,
	ged_draw_group_ref_index_cb cb,
	void *userdata);
GED_EXPORT extern int ged_draw_source_root_subtree_bounds(
	struct ged *gedp,
	vect_t *min,
	vect_t *max,
	int include_overlays);
GED_EXPORT extern int ged_draw_source_root_has_groups(struct ged *gedp);
GED_EXPORT extern int ged_draw_source_root_clear_all_scope_children(
	struct ged *gedp);
GED_EXPORT extern int ged_draw_source_erase_path_at_root(
	struct ged *gedp,
	const char *path);
GED_EXPORT extern int ged_draw_source_erase_path_prefix_at_root(
	struct ged *gedp,
	const char *path);
GED_EXPORT extern int ged_draw_source_erase_groups_by_name_at_root(
	struct ged *gedp,
	const char *name);
GED_EXPORT extern int ged_draw_source_erase_path_in_active_scope(
	struct ged *gedp,
	const char *path,
	void *view_ctx,
	int mode);
GED_EXPORT extern int ged_draw_source_erase_path_prefix_in_active_scope(
	struct ged *gedp,
	const char *path,
	void *view_ctx,
	int mode);
GED_EXPORT extern int ged_draw_source_erase_component_name_in_active_scope(
	struct ged *gedp,
	const char *name,
	void *view_ctx,
	int mode,
	int nonroot_only);
GED_EXPORT extern int ged_draw_source_clear_db_groups_in_scope(
	struct ged *gedp,
	void *view_ctx);
GED_EXPORT extern void ged_draw_log(int level, const char *fmt, ...) _BU_ATTR_PRINTF23;
GED_EXPORT extern void ged_draw_test_force_primitive_face_set_failure(int enable);
GED_EXPORT extern int ged_draw_test_primitive_face_set_failure_enabled(void);
GED_EXPORT extern void ged_draw_view_context_foreach_export_record(
	void *view_ctx,
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

GED_EXPORT extern int ged_draw_shape_ref_refresh_material_color(struct ged *gedp,
								ged_draw_shape_ref ref,
								struct db_i *dbip,
								uint64_t mater_rev);
GED_EXPORT extern int ged_draw_registry_shape_ref_set_indexed_fullpath(
	struct ged *gedp,
	ged_draw_shape_ref shape_ref,
	const struct db_full_path *path);
GED_EXPORT extern int ged_draw_registry_group_ref_set_indexed_fullpath(
	struct ged *gedp,
	ged_draw_group_ref group_ref,
	const struct db_full_path *path);
GED_EXPORT extern const char *ged_draw_registry_shape_ref_semantic_path(
	struct ged *gedp,
	ged_draw_shape_ref shape_ref);
GED_EXPORT extern const char *ged_draw_registry_shape_ref_cached_semantic_path(
	struct ged *gedp,
	ged_draw_shape_ref shape_ref);
GED_EXPORT extern char *ged_draw_registry_shape_ref_cached_fullpath_dup(
	struct ged *gedp,
	ged_draw_shape_ref shape_ref);
GED_EXPORT extern const char *ged_draw_registry_group_ref_semantic_path(
	struct ged *gedp,
	ged_draw_group_ref group_ref);
GED_EXPORT extern int ged_draw_registry_shape_ref_apply_region(
	struct ged *gedp,
	ged_draw_shape_ref shape_ref,
	int region_id,
	int aircode,
	int los,
	int material_id);
GED_EXPORT extern const char *ged_draw_dbpath_skip_lead_slash(const char *s);

GED_EXPORT extern void *ged_draw_active_view_ctx(struct ged *gedp);
GED_EXPORT extern void ged_draw_active_view_ctx_set(struct ged *gedp,
						    void *view_ctx);
GED_EXPORT extern int ged_draw_ensure_root_attached(struct ged *gedp);

GED_EXPORT extern int ged_draw_overlay_geometry_insert(struct ged *gedp,
						       const char *name,
						       const struct ged_draw_overlay_geometry *geometry,
						       long int rgb,
						       fastf_t transparency,
						       int draw_mode,
						       int csoltab,
						       ged_draw_shape_ref *out);
GED_EXPORT extern void ged_draw_overlay_erase_name(struct ged *gedp,
						   const char *name);
GED_EXPORT extern int ged_draw_redraw_group_ref(struct ged *gedp,
						ged_draw_group_ref ref,
						int skip_subtractions);
GED_EXPORT extern int ged_draw_group_ref_redraw_wireframe(struct ged *gedp,
							  ged_draw_group_ref ref,
							  struct db_i *dbip,
							  const struct bn_tol *tol,
							  const struct bg_tess_tol *ttol,
							  void *view_ctx,
							  int skip_subtractions);
GED_EXPORT extern int ged_draw_shape_ref_redraw_wireframe(struct ged *gedp,
							  ged_draw_shape_ref ref,
							  struct db_i *dbip,
							  const struct bn_tol *tol,
							  const struct bg_tess_tol *ttol,
							  void *view_ctx,
							  int skip_subtractions);
GED_EXPORT extern ged_draw_group_ref ged_draw_group_ref_lookup_or_create(struct ged *gedp,
									 const struct db_full_path *dfp);

GED_EXPORT extern int ged_draw_erase_path_string(struct ged *gedp,
							 const char *path);
GED_EXPORT extern int ged_draw_erase_path_prefix_string(struct ged *gedp,
								const char *path);
GED_EXPORT extern int ged_draw_erase_path_string_scoped(struct ged *gedp,
							const char *path,
							void *view_ctx,
							int mode);
GED_EXPORT extern int ged_draw_erase_path_prefix_string_scoped(struct ged *gedp,
							       const char *path,
							       void *view_ctx,
							       int mode);
GED_EXPORT extern int ged_draw_erase_nonroot_component_string_scoped(struct ged *gedp,
								     const char *name,
								     void *view_ctx,
								     int mode);
GED_EXPORT extern int ged_draw_erase_component_string_scoped(struct ged *gedp,
							     const char *name,
							     void *view_ctx,
							     int mode);

GED_EXPORT extern void ged_draw_highlighted_shape_ref_invalidate(struct ged *gedp);

GED_EXPORT extern ged_draw_shape_ref ged_draw_registry_shape_ref_from_source_ref(
	struct ged *gedp,
	ged_draw_scene_handle scene_ref);
GED_EXPORT extern ged_draw_group_ref ged_draw_registry_group_ref_from_source_ref(
	struct ged *gedp,
	ged_draw_scene_handle scene_ref);
GED_EXPORT extern ged_draw_scene_handle ged_draw_registry_shape_ref_scene_handle(
	struct ged *gedp,
	ged_draw_shape_ref ref);
GED_EXPORT extern ged_draw_scene_handle ged_draw_registry_shape_ref_cache_scene_handle(
	struct ged *gedp,
	ged_draw_shape_ref ref);
GED_EXPORT extern ged_draw_scene_handle ged_draw_registry_group_ref_scene_handle(
	struct ged *gedp,
	ged_draw_group_ref ref);
GED_EXPORT extern int ged_draw_shape_ref_record_summary(struct ged *gedp,
							ged_draw_shape_ref ref,
							struct ged_draw_shape_record_summary *out);
GED_EXPORT extern int ged_draw_shape_ref_owning_group_record_summary(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	struct ged_draw_group_record_summary *out);
GED_EXPORT extern int ged_draw_shape_ref_database_source_summary(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	struct ged_draw_database_source_summary *out);
GED_EXPORT extern int ged_draw_obol_shape_ref_index_for_path_hash(
	struct ged *gedp,
	unsigned long long path_hash,
	ged_draw_shape_ref_index_cb cb,
	void *userdata);
GED_EXPORT extern int ged_draw_obol_scene_database_bounds(
	struct ged *gedp,
	vect_t *min,
	vect_t *max,
	int *empty_out);
GED_EXPORT extern int ged_draw_obol_scene_database_autoview_bounds(
	struct ged *gedp,
	vect_t *min,
	vect_t *max,
	int *empty_out);
GED_EXPORT extern void ged_draw_obol_preserved_sources_free(
	struct ged *gedp);
GED_EXPORT extern int ged_draw_obol_highlight_state_set(
	struct ged *gedp,
	int highlighted);
GED_EXPORT extern int ged_draw_obol_scene_root_child_count(
	struct ged *gedp,
	size_t *out);
GED_EXPORT extern void ged_draw_obol_context_tokens_free(struct ged *gedp);
GED_EXPORT extern void ged_draw_obol_scene_context_info_free(
	struct ged_draw_obol_scene_context_info *info);
GED_EXPORT extern int ged_draw_obol_scene_context_info_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_obol_scene_context_info *out);
GED_EXPORT extern int ged_draw_obol_scene_child_context_info_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	struct ged_draw_obol_scene_context_info *out);
GED_EXPORT extern int ged_draw_obol_database_source_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_database_source_summary *out);
GED_EXPORT extern int ged_draw_obol_database_source_publication_begin(
	struct ged *gedp,
	void *view_ctx,
	int publication_draw_mode);
GED_EXPORT extern void ged_draw_obol_database_source_publication_appearance_set(
	struct ged *gedp,
	const struct ged_draw_appearance_settings *settings);
GED_EXPORT extern int ged_draw_obol_database_source_publication_appearance_active(
	struct ged *gedp);
GED_EXPORT extern void ged_draw_obol_database_source_publication_group_set(
	struct ged *gedp,
	const char *group_path);
GED_EXPORT extern void ged_draw_obol_database_source_publication_end(
	struct ged *gedp);
GED_EXPORT extern int ged_draw_obol_apply_draw_paths(
	struct ged *gedp,
	void *view_ctx,
	const char **paths,
	int path_count,
	const struct ged_draw_appearance_settings *settings,
	struct ged_draw_transaction_result *result);
GED_EXPORT extern int ged_draw_obol_database_source_paths_foreach(
	struct ged *gedp,
	int skip_overlay_groups,
	ged_draw_obol_database_source_path_cb cb,
	void *userdata);
GED_EXPORT extern int ged_draw_obol_database_source_records_foreach(
	struct ged *gedp,
	int skip_overlay_groups,
	ged_draw_obol_database_source_record_cb cb,
	void *userdata);
GED_EXPORT extern int ged_draw_obol_visible_database_source_records_foreach_fast(
	struct ged *gedp,
	ged_draw_obol_database_source_record_cb cb,
	void *userdata);
GED_EXPORT extern int ged_draw_obol_shape_paths_foreach(
	struct ged *gedp,
	int skip_overlay_groups,
	ged_draw_obol_database_source_path_cb cb,
	void *userdata);
GED_EXPORT extern int ged_draw_obol_group_database_source_paths_foreach(
	struct ged *gedp,
	const char *group_path,
	ged_draw_obol_database_source_path_cb cb,
	void *userdata);
GED_EXPORT extern int ged_draw_obol_group_database_source_records_foreach(
	struct ged *gedp,
	const char *group_path,
	ged_draw_obol_database_source_record_cb cb,
	void *userdata);
GED_EXPORT extern int ged_draw_obol_database_source_owner_group_path_for_path(
	struct ged *gedp,
	const char *path,
	struct bu_vls *out);
GED_EXPORT extern int ged_draw_obol_database_source_display_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_scene_display_summary *out);
GED_EXPORT extern int ged_draw_obol_database_source_bounds_for_path(
	struct ged *gedp,
	const char *path,
	vect_t *min,
	vect_t *max,
	int include_overlays,
	int *empty_out);
GED_EXPORT extern int ged_draw_obol_group_display_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_scene_display_summary *out);
GED_EXPORT extern int ged_draw_obol_group_subtree_bounds_for_path(
	struct ged *gedp,
	const char *path,
	vect_t *min,
	vect_t *max,
	int include_overlays,
	int *empty_out);
GED_EXPORT extern int ged_draw_obol_database_source_last_point_for_path(
	struct ged *gedp,
	const char *path,
	point_t out);
GED_EXPORT extern int ged_draw_obol_database_source_line_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_view_line_summary *out);
GED_EXPORT extern int ged_draw_obol_database_source_line_point_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	point_t out);
GED_EXPORT extern int ged_draw_obol_database_source_line_command_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	int *out);
GED_EXPORT extern int ged_draw_obol_database_source_surface_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_view_surface_summary *out);
GED_EXPORT extern int ged_draw_obol_database_source_surface_point_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	point_t out);
GED_EXPORT extern int ged_draw_obol_database_source_surface_index_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	int *out);
GED_EXPORT extern int ged_draw_obol_shape_display_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_scene_display_summary *out);
GED_EXPORT extern int ged_draw_obol_shape_geometry_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_shape_geometry_summary *out);
GED_EXPORT extern int ged_draw_obol_shape_surface_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_view_surface_summary *out);
GED_EXPORT extern int ged_draw_obol_shape_surface_point_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	point_t out);
GED_EXPORT extern int ged_draw_obol_shape_surface_index_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	int *out);
GED_EXPORT extern int ged_draw_obol_shape_line_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_view_line_summary *out);
GED_EXPORT extern int ged_draw_obol_shape_line_point_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	point_t out);
GED_EXPORT extern int ged_draw_obol_shape_line_command_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	int *out);
GED_EXPORT extern int ged_draw_obol_database_source_translate_vlist_for_path(
	struct ged *gedp,
	const char *path,
	const vect_t xlate);
GED_EXPORT extern int ged_draw_obol_database_source_clear_vlist_for_path(
	struct ged *gedp,
	const char *path);
GED_EXPORT extern int
ged_draw_obol_database_source_publish_annotation_line_set_for_path(
	struct ged *gedp,
	const char *path,
	const point_t *points,
	const int *commands,
	size_t point_count);
GED_EXPORT extern int
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
GED_EXPORT extern int
ged_draw_obol_database_source_annotation_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_obol_annotation_summary *out);
GED_EXPORT extern int
ged_draw_obol_database_source_annotation_point_at_for_path(
	struct ged *gedp,
	const char *path,
	size_t index,
	point_t out);
GED_EXPORT extern int
ged_draw_obol_database_source_publish_auxiliary_line_set_for_path(
	struct ged *gedp,
	const char *path,
	const char *name,
	const point_t *points,
	const int *commands,
	size_t point_count,
	const struct ged_draw_scene_display_summary *display_state);
GED_EXPORT extern int
ged_draw_obol_database_source_publish_auxiliary_source_line_set_for_path(
	struct ged *gedp,
	const char *path,
	const char *source_path,
	const char *display_name,
	const point_t *points,
	const int *commands,
	size_t point_count,
	const struct ged_draw_scene_display_summary *display_state);
GED_EXPORT extern int
ged_draw_obol_database_source_clear_auxiliary_shapes_for_path(
	struct ged *gedp,
	const char *path);
GED_EXPORT extern int ged_draw_obol_database_source_clear_mesh_for_path(
	struct ged *gedp,
	const char *path);
GED_EXPORT extern int ged_draw_obol_local_shape_publish_line_set_for_path(
	struct ged *gedp,
	const char *group_path,
	const char *shape_path,
	const char *display_name,
	const point_t *points,
	const int *commands,
	size_t point_count,
	const struct ged_draw_scene_display_summary *display_state);
GED_EXPORT extern int ged_draw_obol_local_shape_set_record_role_for_path(
	struct ged *gedp,
	const char *shape_path,
	const char *record_role);
GED_EXPORT extern int
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
GED_EXPORT extern int ged_draw_obol_database_source_set_vlist_center_for_path(
	struct ged *gedp,
	const char *path,
	const point_t center);
GED_EXPORT extern int
ged_draw_obol_database_source_update_vlist_bounds_for_path(
	struct ged *gedp,
	const char *path);
GED_EXPORT extern int ged_draw_obol_database_source_set_placement_for_path(
	struct ged *gedp,
	const char *path,
	int draw_mat_valid,
	const mat_t draw_mat,
	int draw_center_valid,
	const point_t draw_center,
	int draw_size_valid,
	fastf_t draw_size);
GED_EXPORT extern int ged_draw_obol_database_source_set_bounds_for_path(
	struct ged *gedp,
	const char *path,
	const point_t bmin,
	const point_t bmax);
GED_EXPORT extern int ged_draw_obol_database_source_set_display_name_for_path(
	struct ged *gedp,
	const char *path,
	const char *display_name);
GED_EXPORT extern int ged_draw_obol_database_source_geometry_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_shape_geometry_summary *out);
GED_EXPORT extern int ged_draw_obol_database_source_geometry_summary_for_path_mode(
	struct ged *gedp,
	const char *path,
	int draw_mode_valid,
	int ged_draw_mode,
	struct ged_draw_shape_geometry_summary *out);
GED_EXPORT extern int ged_draw_obol_database_source_material_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_shape_material_summary *out);
GED_EXPORT extern int ged_draw_obol_database_source_refresh_material_color_for_path(
	struct ged *gedp,
	const char *path,
	struct db_i *dbip,
	uint64_t material_revision);
GED_EXPORT extern int ged_draw_obol_database_sources_refresh_material_colors(
	struct ged *gedp,
	struct db_i *dbip,
	uint64_t material_revision);
GED_EXPORT extern int ged_draw_obol_database_source_evaluated_region_for_path(
	struct ged *gedp,
	const char *path,
	int *out);
GED_EXPORT extern int ged_draw_obol_database_source_set_evaluated_region_for_path(
	struct ged *gedp,
	const char *path,
	int evaluated_region);
GED_EXPORT extern int ged_draw_obol_database_source_set_region_metadata_for_path(
	struct ged *gedp,
	const char *path,
	int region_id,
	int aircode,
	int los,
	int material_id);
GED_EXPORT extern int ged_draw_obol_database_source_apply_draw_metadata_for_path(
	struct ged *gedp,
	const char *path,
	const struct BRLObolDrawMetadataRecord *record);
GED_EXPORT extern int ged_draw_obol_database_source_update_display_for_path(
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
GED_EXPORT extern int ged_draw_obol_database_source_set_selected_for_instance_key(
	struct ged *gedp,
	const char *instance_key,
	int selected);
GED_EXPORT extern int ged_draw_obol_database_sources_sync_selected_paths(
	struct ged *gedp,
	const char *const *paths,
	size_t path_count);
GED_EXPORT extern int ged_draw_obol_shape_update_display_for_path(
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
GED_EXPORT extern int ged_draw_obol_database_source_update_appearance_for_path(
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
GED_EXPORT extern int ged_draw_obol_database_source_mark_stale_for_path(
	struct ged *gedp,
	const char *path,
	ged_draw_stale_reason reason);
GED_EXPORT extern int ged_draw_obol_database_source_record_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_obol_database_source_record *out);
GED_EXPORT extern int ged_draw_obol_database_source_runtime_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_obol_database_source_runtime *out);
GED_EXPORT extern int ged_draw_obol_database_source_set_mesh_lod_for_path(
	struct ged *gedp,
	const char *path,
	struct BRLObolMeshLod *lod);
GED_EXPORT extern int ged_draw_obol_database_source_set_mesh_lod_bounds_for_path(
	struct ged *gedp,
	const char *path,
	const point_t bmin,
	const point_t bmax);
GED_EXPORT extern int ged_draw_obol_database_source_apply_record_for_path(
	struct ged *gedp,
	const char *path,
	const struct ged_draw_obol_database_source_record *record);
GED_EXPORT extern int ged_draw_obol_database_source_draw_state_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_obol_draw_state_summary *out);
GED_EXPORT extern int
ged_draw_obol_database_source_set_realization_for_path(
	struct ged *gedp,
	const char *path,
	int current,
	int failed,
	uint64_t realized_source_revision,
	uint64_t realized_inputs_revision,
	ged_draw_stale_reason reason);
GED_EXPORT extern int
ged_draw_obol_database_source_realization_policy_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_obol_realization_policy_summary *out);
GED_EXPORT extern int
ged_draw_obol_database_source_set_realization_roles_for_path(
	struct ged *gedp,
	const char *path,
	int role_flags);
GED_EXPORT extern int
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
GED_EXPORT extern int
ged_draw_obol_database_source_realize_for_path(struct ged *gedp,
	const char *path);
GED_EXPORT extern int
ged_draw_obol_database_source_realize_pending(struct ged *gedp);

int
ged_draw_obol_database_sources_redraw(struct ged *gedp,
	void *view_ctx,
	const char *path,
	int draw_mode);
GED_EXPORT extern int ged_draw_obol_scene_controller_ensure_owned(
	struct ged *gedp,
	int sync_current_scene);
GED_EXPORT extern int ged_draw_obol_scene_controller_attached(
	struct ged *gedp);
GED_EXPORT extern int ged_draw_obol_scene_controller_full_synced(
	struct ged *gedp);
GED_EXPORT extern int ged_draw_obol_scene_sync_attached_transaction(
	struct ged *gedp,
	const struct ged_draw_transaction *txn,
	const struct ged_draw_transaction_result *result);
GED_EXPORT extern int ged_draw_obol_database_source_ensure_for_path(
	struct ged *gedp,
	const char *path,
	struct db_i *dbip,
	int ged_draw_mode,
	uint64_t source_revision);
GED_EXPORT extern int ged_draw_obol_database_source_ensure_for_path_with_placement(
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
GED_EXPORT extern int ged_draw_obol_database_source_rename_for_path(
	struct ged *gedp,
	const char *path,
	const char *new_path,
	uint64_t source_revision);
GED_EXPORT extern int ged_draw_obol_database_source_move_to_group_for_path(
	struct ged *gedp,
	const char *source_path,
	const char *group_path);
GED_EXPORT extern int ged_draw_obol_database_source_remove_for_path(
	struct ged *gedp,
	const char *path);
GED_EXPORT extern int
ged_draw_obol_database_sources_remove_for_path_prefix(
	struct ged *gedp,
	const char *path_prefix);
GED_EXPORT extern int
ged_draw_obol_database_sources_remove_for_path_prefix_in_scope(
	struct ged *gedp,
	const char *path_prefix,
	void *view_ctx);
GED_EXPORT extern int
ged_draw_obol_database_sources_remove_for_path_prefix_in_scope_mode(
	struct ged *gedp,
	const char *path_prefix,
	void *view_ctx,
	int mode);
GED_EXPORT extern int
ged_draw_obol_active_database_sources_remove_for_path_prefix(
	struct ged *gedp,
	const char *path_prefix);
GED_EXPORT extern int
ged_draw_obol_active_database_sources_remove_for_path_prefix_mode(
	struct ged *gedp,
	const char *path_prefix,
	int mode);
GED_EXPORT extern int
ged_draw_obol_database_sources_remove_for_component_name(
	struct ged *gedp,
	const char *name,
	int nonroot_only,
	int mode);
GED_EXPORT extern int ged_draw_obol_database_sources_clear(struct ged *gedp);
GED_EXPORT extern int ged_draw_obol_database_sources_clear_in_scope(
	struct ged *gedp,
	void *view_ctx);
GED_EXPORT extern int ged_draw_obol_scene_clear(struct ged *gedp);
GED_EXPORT extern int ged_draw_obol_active_scene_clear(struct ged *gedp);
GED_EXPORT extern int ged_draw_obol_groups_remove_for_component_name(
	struct ged *gedp,
	const char *name);
GED_EXPORT extern int ged_draw_obol_groups_remove_for_path_prefix(
	struct ged *gedp,
	const char *path_prefix);
GED_EXPORT extern int ged_draw_obol_group_remove_for_path(
	struct ged *gedp,
	const char *path);
GED_EXPORT extern int ged_draw_obol_group_clear_for_path(
	struct ged *gedp,
	const char *path);
GED_EXPORT extern int ged_draw_obol_group_rename_for_path(
	struct ged *gedp,
	const char *path,
	const char *new_path);
GED_EXPORT extern int ged_draw_obol_group_erase_subpath_for_path(
	struct ged *gedp,
	const char *parent_path,
	const char *subpath);
GED_EXPORT extern int ged_draw_obol_group_update_display_for_path(
	struct ged *gedp,
	const char *path,
	int visible_valid,
	int visible);
GED_EXPORT extern int ged_draw_obol_group_ensure_for_path(
	struct ged *gedp,
	const char *path,
	const char *intent_path,
	int mode,
	int overlay);
GED_EXPORT extern int ged_draw_obol_group_record_summary_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_group_record_summary *out);
GED_EXPORT extern int ged_draw_obol_group_paths_foreach(
	struct ged *gedp,
	int skip_overlay_groups,
	ged_draw_obol_group_path_cb cb,
	void *userdata);
GED_EXPORT extern int ged_draw_obol_group_shape_count_for_path(
	struct ged *gedp,
	const char *path,
	int *out);
GED_EXPORT extern int ged_draw_obol_group_descendant_group_count_for_path(
	struct ged *gedp,
	const char *path,
	size_t *out);
GED_EXPORT extern int ged_draw_obol_group_child_count_for_path(
	struct ged *gedp,
	const char *path,
	size_t *out);
GED_EXPORT extern int ged_draw_obol_group_update_appearance_for_path(
	struct ged *gedp,
	const char *path,
	const struct ged_draw_appearance_settings *settings);
GED_EXPORT extern int ged_draw_obol_group_appearance_for_path(
	struct ged *gedp,
	const char *path,
	struct ged_draw_appearance_settings *settings);
GED_EXPORT extern int ged_draw_obol_group_update_draw_intent_for_path(
	struct ged *gedp,
	const char *path,
	const char *intent_path,
	int mode_valid,
	int mode,
	int overlay_valid,
	int overlay);
GED_EXPORT extern int ged_draw_group_ref_record_summary(struct ged *gedp,
							 ged_draw_group_ref ref,
							 struct ged_draw_group_record_summary *out);
GED_EXPORT extern int ged_draw_group_ref_tree_summary(struct ged *gedp,
						      ged_draw_group_ref ref,
						      struct ged_draw_scene_tree_summary *out);
GED_EXPORT extern int ged_draw_group_ref_shape_count(struct ged *gedp,
						     ged_draw_group_ref ref,
						     int *out);

GED_EXPORT extern int ged_draw_group_ref_index_for_component(struct ged *gedp,
							    const char *path,
							    ged_draw_group_ref_index_cb cb,
							    void *userdata);
GED_EXPORT extern int ged_draw_obol_shape_ref_index_for_component(
	struct ged *gedp,
	const char *path,
	ged_draw_shape_ref_index_cb cb,
	void *userdata);
GED_EXPORT extern int ged_draw_obol_group_ref_index_for_component(
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

GED_EXPORT extern void ged_draw_index_stats_get(struct ged *gedp,
						struct ged_draw_index_stats *stats);
GED_EXPORT extern void ged_draw_index_stats_reset(struct ged *gedp);

GED_EXPORT extern int ged_draw_group_ref_set_dbpath(struct ged *gedp,
						    ged_draw_group_ref ref,
						    const struct db_full_path *new_dfp);
GED_EXPORT extern int ged_draw_group_ref_set_mode(struct ged *gedp,
						  ged_draw_group_ref ref,
						  int mode);
GED_EXPORT extern int ged_draw_group_ref_set_appearance_settings(struct ged *gedp,
								ged_draw_group_ref ref,
								const struct ged_draw_appearance_settings *settings);
GED_EXPORT extern int ged_draw_group_ref_appearance_settings(struct ged *gedp,
							     ged_draw_group_ref ref,
							     struct ged_draw_appearance_settings *settings);
GED_EXPORT extern int ged_draw_group_ref_set_visible(struct ged *gedp,
						     ged_draw_group_ref ref,
						     int visible);
GED_EXPORT extern int ged_draw_shape_ref_set_visible(struct ged *gedp,
						     ged_draw_shape_ref ref,
						     int visible);
GED_EXPORT extern int ged_draw_shape_ref_get_color(struct ged *gedp,
					   ged_draw_shape_ref ref,
					   unsigned char rgb[3]);
GED_EXPORT extern int ged_draw_shape_ref_set_color(struct ged *gedp,
					   ged_draw_shape_ref ref,
					   const unsigned char rgb[3]);
GED_EXPORT extern int ged_draw_shape_ref_set_highlighted(struct ged *gedp,
							 ged_draw_shape_ref ref,
							 int highlighted);
GED_EXPORT extern int ged_draw_shape_ref_set_selected(struct ged *gedp,
						      ged_draw_shape_ref ref,
						      int selected);
GED_EXPORT extern int ged_draw_shape_ref_set_transparency(struct ged *gedp,
							  ged_draw_shape_ref ref,
							  fastf_t transparency);
GED_EXPORT extern int ged_draw_shape_ref_mark_database_source_changed(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	ged_draw_stale_reason reason);
GED_EXPORT extern int ged_draw_shape_ref_lod_ensure(struct ged *gedp,
						    ged_draw_shape_ref ref,
						    void *first_view_ctx,
						    void **view_ctxs,
						    size_t view_ctx_count);
GED_EXPORT extern void *ged_draw_shape_ref_view_context(struct ged *gedp,
							ged_draw_shape_ref ref);
GED_EXPORT extern int ged_draw_view_selection_add_shape_ref_context(
	struct ged *gedp,
	void *view_ctx,
	ged_draw_shape_ref ref,
	void **selection_view_ctx,
    struct bu_vls *path);
GED_EXPORT extern void ged_draw_registry_source_ref_highlight_free(
	ged_draw_scene_handle scene_ref);

__END_DECLS

#endif /* LIBGED_GED_DRAW_PRIVATE_H */
