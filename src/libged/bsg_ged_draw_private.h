/*              B S G _ G E D _ D R A W _ P R I V A T E . H
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
/** @file bsg_ged_draw_private.h
 *
 * Private node-facing draw-scene bridge for libged implementation files.
 *
 * Installed GED drawing APIs are record/ref/transaction based.  These helpers
 * are intentionally not installed: they preserve the current BSG-backed draw
 * implementation while public callers migrate to semantic records and refs.
 */

#ifndef LIBGED_BSG_GED_DRAW_PRIVATE_H
#define LIBGED_BSG_GED_DRAW_PRIVATE_H

#include "common.h"

#include <stdint.h>

#include "bu/vls.h"
#include "vmath.h"

#include "bsg/geometry.h"
#include "bsg/scene_builder.h"
#include "./bsg_ged_draw_view_private.h"
#include "ged/draw.h"
#include "rt/db_fullpath.h"
#include "rt/view.h"

__BEGIN_DECLS

struct bsg_view;
struct bsg_view_set;
struct bu_ptbl;
struct bu_vls;
struct bg_tess_tol;
struct bn_tol;
struct db_i;
struct directory;
struct rt_mesh_lod;
struct rt_db_internal;
struct rt_bot_internal;
struct rt_brep_internal;
struct rt_pg_internal;
struct db_tree_state;
struct ged_draw_source_state;
struct model;
struct nmgregion;

typedef enum {
    GED_DRAW_SHAPE_USER_DATA_NONE = 0,
    GED_DRAW_SHAPE_USER_DATA_RT_DB_INTERNAL
} ged_draw_shape_user_data_kind;

typedef struct ged_draw_shape_draft ged_draw_shape_draft;

typedef struct ged_draw_shape_state {
    struct db_full_path s_fullpath;
    struct directory *leaf_dp;
    char *display_name;
    unsigned long long path_hash;
    int region_id;
    int aircode;
    int los;
    int material_id;
    uint64_t source_revision;
    uint64_t inputs_revision;
    uint64_t realized_source_revision;
    uint64_t realized_inputs_revision;
    ged_draw_stale_reason stale_reason;
    void (*invalidate)(bsg_scene_ref, void *, void *);
    void *invalidate_data;
    struct ged *gedp;
    void *u_data;
    ged_draw_shape_user_data_kind u_data_kind;
    bsg_scene_ref source_ref;
    struct ged_draw_source_state *source_data;
    size_t geometry_command_count;
    uint64_t geometry_revision;
} ged_draw_shape_state;

enum ged_draw_overlay_geometry_kind {
    GED_DRAW_OVERLAY_GEOMETRY_NONE = 0,
    GED_DRAW_OVERLAY_GEOMETRY_LINE_SET,
    GED_DRAW_OVERLAY_GEOMETRY_POINT_SET,
    GED_DRAW_OVERLAY_GEOMETRY_INDEXED_FACE_SET
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

GED_EXPORT extern bsg_scene_ref ged_draw_view_context_overlay_geometry_create(
	void *view_ctx,
	const char *name,
	enum ged_draw_overlay_geometry_kind kind);
GED_EXPORT extern bsg_scene_ref ged_draw_view_context_overlay_scene_find(void *view_ctx,
									 const char *name);
GED_EXPORT extern void ged_draw_view_context_overlay_name_erase(void *view_ctx,
								const char *name);
GED_EXPORT extern int ged_draw_view_context_overlay_scene_append(void *view_ctx,
								 bsg_scene_ref scene);
GED_EXPORT extern int ged_draw_view_overlay_command_result_owner_set(bsg_scene_ref scene,
								     const void *owner,
								     const char *source_path);
GED_EXPORT extern bsg_scene_ref ged_draw_view_context_overlay_create(void *view_ctx,
								     const char *name);
GED_EXPORT extern int ged_draw_view_context_overlay_internal_create(
	struct ged *gedp,
	void *view_ctx,
	const char *name,
	struct db_full_path *fp,
	struct rt_db_internal **ip,
	bsg_scene_ref *out);

struct ged_draw_source_state {
    struct db_i *dbip;
    struct db_full_path *fp;
    const struct bn_tol *tol;
    const struct bg_tess_tol *ttol;
    struct rt_mesh_lod *rt_mesh_lod;
    point_t mesh_lod_bmin;
    point_t mesh_lod_bmax;
    int mesh_lod_bounds_valid;
    uint64_t source_revision;
    uint64_t inputs_revision;
    uint64_t realized_source_revision;
    uint64_t realized_inputs_revision;
    ged_draw_stale_reason stale_reason;
};

typedef enum {
    GED_DRAW_DATABASE_SOURCE_MATERIAL_INHERIT = 0,
    GED_DRAW_DATABASE_SOURCE_MATERIAL_DATABASE = 1
} ged_draw_database_source_material_policy;

typedef enum {
    GED_DRAW_DATABASE_SOURCE_REALIZATION_STALE = 0,
    GED_DRAW_DATABASE_SOURCE_REALIZATION_CURRENT = 1
} ged_draw_database_source_realization_status;

typedef enum {
    GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_NONE = 0,
    GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_CSG = 1,
    GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_MESH = 2
} ged_draw_database_source_realization_role;

struct ged_draw_database_source_record {
    const char *database_path;
    int draw_mode;
    ged_draw_database_source_material_policy material_policy;
    uint64_t source_revision;
    uint64_t inputs_revision;
    uint64_t realized_source_revision;
    uint64_t realized_inputs_revision;
    ged_draw_stale_reason stale_reason;
    uint64_t realization_identity;
    ged_draw_database_source_realization_status realization_status;
    int realization_role_flags;
    int realization_view_dependent;
    fastf_t realization_view_scale;
    uint64_t realization_bot_threshold;
    fastf_t realization_curve_scale;
    fastf_t realization_point_scale;
};

struct ged_draw_source_state_record {
    struct db_i *dbip;
    const struct db_full_path *fullpath;
    const struct bn_tol *tol;
    const struct bg_tess_tol *ttol;
    ged_draw_stale_reason stale_reason;
};

struct ged_draw_shape_record_summary {
    const struct db_full_path *fullpath;
    const char *display_name;
    const char *leaf_name;
    bsg_scene_ref owning_group_ref;
    unsigned long long path_hash;
    uint64_t source_revision;
    uint64_t inputs_revision;
    uint64_t realized_source_revision;
    uint64_t realized_inputs_revision;
    int stale;
    const char *stale_reason;
    int visible;
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

struct ged_draw_scene_draw_state_summary {
    const char *name;
    int draw_mode;
    int line_style;
    int pipeline_candidate;
    int draw_mat_valid;
    mat_t draw_mat;
    int bounds_valid;
    point_t bounds_min;
    point_t bounds_max;
    void *view_ctx;
};

struct ged_draw_source_runtime_summary {
    struct db_i *dbip;
    const struct db_full_path *fullpath;
    struct directory *leaf_dp;
    const char *name;
    const struct bn_tol *tol;
    const struct bg_tess_tol *ttol;
    struct rt_mesh_lod *rt_mesh_lod;
    int mesh_lod_bounds_valid;
    point_t mesh_lod_bmin;
    point_t mesh_lod_bmax;
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

GED_EXPORT extern void ged_draw_source_data_free(void *data);
GED_EXPORT extern void ged_draw_scene_ref_realization_set_view_context_policy(bsg_scene_ref ref,
									      const void *view_ctx);
GED_EXPORT extern int ged_draw_brep_mesh_lod_detail_setup(struct rt_mesh_lod *lod,
							  struct db_i *dbip,
							  struct directory *dp,
							  const struct bg_tess_tol *ttol,
							  const struct bn_tol *tol);
GED_EXPORT extern int ged_draw_brep_mesh_lod_cache_prepare(struct rt_mesh_lod **lod,
							   point_t bmin,
							   point_t bmax,
							   int *bounds_valid,
							   struct db_i *dbip,
							   struct directory *dp,
							   const struct bg_tess_tol *ttol,
							   const struct bn_tol *tol);
GED_EXPORT extern int ged_draw_scene_ref_detach(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_append_child(bsg_scene_ref parent_ref,
						      bsg_scene_ref child_ref);
GED_EXPORT extern bsg_scene_ref ged_draw_view_context_group_child_ensure(
	void *view_ctx,
	bsg_scene_ref parent_ref,
	const char *name,
	void *dp_hint);
GED_EXPORT extern bsg_scene_ref ged_draw_view_context_group_create(
	void *view_ctx,
	const char *name);
GED_EXPORT extern bsg_scene_ref ged_draw_view_context_database_source_create(
	void *view_ctx,
	const char *name);
GED_EXPORT extern bsg_scene_ref ged_draw_view_context_geometry_create(
	void *view_ctx,
	const char *name);
GED_EXPORT extern int ged_draw_scene_ref_is_database_source(
	bsg_scene_ref ref);
GED_EXPORT extern void ged_draw_scene_ref_erase_nested_subpath(
	bsg_scene_ref parent_ref,
	const char * const *comp_names,
	size_t comp_count,
	size_t depth_start,
	int (*match_fn)(bsg_scene_ref shape_ref, void *match_ctx),
	void *match_ctx);
GED_EXPORT extern void ged_draw_scene_ref_free_group(bsg_scene_ref ref);
GED_EXPORT extern void ged_draw_scene_ref_free_group_contents(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_attach_source_state_record(
	bsg_scene_ref ref,
	const struct ged_draw_source_state_record *record);
GED_EXPORT extern int ged_draw_scene_ref_database_source_record_get(
	bsg_scene_ref ref,
	struct ged_draw_database_source_record *record);
GED_EXPORT extern int ged_draw_scene_ref_database_source_record_apply(
	bsg_scene_ref ref,
	const struct ged_draw_database_source_record *record);
GED_EXPORT extern void ged_draw_log(int level, const char *fmt, ...) _BU_ATTR_PRINTF23;
GED_EXPORT extern void ged_draw_test_force_primitive_face_set_failure(int enable);
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

GED_EXPORT extern ged_draw_shape_draft *ged_draw_shape_draft_create_context(struct ged *gedp,
									    void *view_ctx,
									    int registered);
GED_EXPORT extern void ged_draw_shape_draft_destroy(ged_draw_shape_draft *draft);
GED_EXPORT extern int ged_draw_scene_ref_create_draft_pair(struct ged *gedp,
							   void *view_ctx,
							   bsg_scene_ref *source_ref,
							   bsg_scene_ref *shape_ref);
GED_EXPORT extern int ged_draw_shape_draft_publish_line_set(ged_draw_shape_draft *draft,
							    const point_t *points,
							    const int *commands,
							    size_t point_count);
GED_EXPORT extern int ged_draw_shape_draft_publish_primitive_face_set(ged_draw_shape_draft *draft,
								      struct rt_db_internal *ip,
								      const struct bg_tess_tol *ttol,
								      const struct bn_tol *tol,
								      const struct rt_view_info *view_info);
GED_EXPORT extern int ged_draw_shape_draft_publish_bot_wireframe_line_set(ged_draw_shape_draft *draft,
									  const struct rt_bot_internal *bot);
GED_EXPORT extern int ged_draw_shape_draft_publish_brep_wireframe_line_set(ged_draw_shape_draft *draft,
									   const struct rt_brep_internal *bi,
									   const struct bn_tol *tol);
GED_EXPORT extern int ged_draw_shape_draft_publish_poly_wireframe_line_set(ged_draw_shape_draft *draft,
									   const struct rt_pg_internal *pg);
GED_EXPORT extern int ged_draw_shape_draft_publish_nmg_region(ged_draw_shape_draft *draft,
							      const struct nmgregion *r,
							      int style);
GED_EXPORT extern int ged_draw_shape_draft_apply_known_bounds(ged_draw_shape_draft *draft,
								      const point_t min,
								      const point_t max);
GED_EXPORT extern int ged_draw_shape_draft_publish_primitive_wireframe(ged_draw_shape_draft *draft,
								       struct rt_db_internal *ip,
								       const struct bg_tess_tol *ttol,
								       const struct bn_tol *tol,
								       void *view_ctx,
								       int adaptive);
GED_EXPORT extern int ged_draw_shape_draft_apply_tree_result_state(
	ged_draw_shape_draft *draft,
	int dashed,
	int has_region,
	int region_id,
	int aircode,
	int los,
	int material_id,
	const struct ged_draw_appearance_settings *settings);
GED_EXPORT extern int ged_draw_shape_draft_apply_tree_legacy_color(
	ged_draw_shape_draft *draft,
	const unsigned char wireframe_color_override[3],
	const struct db_tree_state *tsp);
GED_EXPORT extern void color_soltab_scene_ref(struct db_i *dbip, bsg_scene_ref shape_ref);
GED_EXPORT extern int ged_draw_scene_ref_refresh_material_color(
	struct db_i *dbip,
	bsg_scene_ref shape_ref,
	uint64_t mater_rev);
GED_EXPORT extern int ged_draw_shape_draft_apply_path_source_state(
	ged_draw_shape_draft *draft,
	struct db_i *dbip,
	const struct db_full_path *pathp,
	const struct bn_tol *tol,
	const struct bg_tess_tol *ttol,
	int has_draw_mat,
	const mat_t draw_mat,
	const char *display_name);
GED_EXPORT extern int ged_draw_shape_draft_apply_late_display_state(
	ged_draw_shape_draft *draft,
	const struct db_full_path *path,
	int line_style,
	const struct ged_draw_appearance_settings *settings,
	const unsigned char material_rgb[3],
	int highlighted);
GED_EXPORT extern int ged_draw_shape_draft_apply_evaluated_path_display(
	ged_draw_shape_draft *draft,
	const struct ged_draw_appearance_settings *settings,
	const unsigned char material_rgb[3],
	int has_transform,
	const mat_t transform,
	int is_subtraction);
GED_EXPORT extern int ged_draw_shape_draft_apply_database_leaf_display(
	ged_draw_shape_draft *draft,
	const struct ged_draw_appearance_settings *settings,
	int bool_op,
	const unsigned char material_rgb[3],
	int has_draw_size,
	fastf_t draw_size);
GED_EXPORT extern ged_draw_shape_ref ged_draw_create_evaluated_path_shape_ref(
	struct ged *gedp,
	void *view_ctx,
	const char *path,
	const struct ged_draw_appearance_settings *settings);
GED_EXPORT extern int ged_draw_append_tree_shape_to_group(
	struct ged *gedp,
	void *view_ctx,
	ged_draw_group_ref group_ref,
	const struct ged_draw_appearance_settings *settings,
	const struct db_full_path *pathp,
	struct db_tree_state *tsp,
	struct rt_db_internal *ip);
GED_EXPORT extern int ged_draw_add_tree_line_set_to_group(
	struct ged *gedp,
	void *view_ctx,
	ged_draw_group_ref group_ref,
	const struct ged_draw_appearance_settings *settings,
	int dashflag,
	const struct db_full_path *pathp,
	struct db_tree_state *tsp,
	const unsigned char wireframe_color_override[3],
	const point_t *points,
	const int *commands,
	size_t point_count);
GED_EXPORT extern int ged_draw_add_tree_nmg_region_to_group(
	struct ged *gedp,
	void *view_ctx,
	ged_draw_group_ref group_ref,
	const struct ged_draw_appearance_settings *settings,
	int dashflag,
	const struct db_full_path *pathp,
	struct db_tree_state *tsp,
	const unsigned char wireframe_color_override[3],
	const struct nmgregion *r,
	int style);
GED_EXPORT extern int ged_draw_add_tree_primitive_face_set_to_group(
	struct ged *gedp,
	void *view_ctx,
	ged_draw_group_ref group_ref,
	const struct ged_draw_appearance_settings *settings,
	int dashflag,
	const struct db_full_path *pathp,
	struct db_tree_state *tsp,
	const unsigned char wireframe_color_override[3],
	const struct rt_db_internal *ip,
	int force_failure);
GED_EXPORT extern int ged_draw_add_tree_primitive_wireframe_to_group(
	struct ged *gedp,
	void *view_ctx,
	ged_draw_group_ref group_ref,
	const struct ged_draw_appearance_settings *settings,
	int dashflag,
	const struct db_full_path *pathp,
	struct db_tree_state *tsp,
	const unsigned char wireframe_color_override[3],
	const struct rt_db_internal *ip);
GED_EXPORT extern ged_draw_shape_ref ged_draw_shape_draft_commit_to_group(ged_draw_shape_draft *draft,
									  ged_draw_group_ref group_ref);
GED_EXPORT extern int ged_draw_shape_draft_commit_database_leaf_to_scene_ref(
	ged_draw_shape_draft *draft,
	bsg_scene_ref parent_ref);

GED_EXPORT extern bsg_scene_ref ged_draw_scene_ref_null(void);
GED_EXPORT extern bsg_scene_ref ged_draw_scene_ref_from_context(void *scene_ctx);
GED_EXPORT extern int ged_draw_scene_ref_is_null(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_equal(bsg_scene_ref a, bsg_scene_ref b);
GED_EXPORT extern bsg_scene_ref ged_draw_scene_ref_parent(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_foreach_child(
	bsg_scene_ref ref,
	int (*cb)(bsg_scene_ref child_ref, void *userdata),
	void *userdata);
GED_EXPORT extern void *ged_draw_scene_ref_parent_context(bsg_scene_ref ref);
GED_EXPORT extern void *ged_draw_scene_ref_registry_owner(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_display_summary(
	bsg_scene_ref ref,
	struct ged_draw_scene_display_summary *out);
GED_EXPORT extern int ged_draw_scene_ref_source_summary(
	bsg_scene_ref ref,
	struct ged_draw_database_source_summary *out);
GED_EXPORT extern bsg_scene_ref ged_draw_scene_ref_source_owner(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_tree_summary(
	bsg_scene_ref ref,
	struct ged_draw_scene_tree_summary *out);
GED_EXPORT extern int ged_draw_scene_ref_material_summary(
	bsg_scene_ref ref,
	struct ged_draw_shape_material_summary *out);
GED_EXPORT extern int ged_draw_scene_ref_shape_record_summary(
	bsg_scene_ref ref,
	struct ged_draw_shape_record_summary *out);
GED_EXPORT extern int ged_draw_scene_ref_group_record_summary(
	bsg_scene_ref ref,
	struct ged_draw_group_record_summary *out);
GED_EXPORT extern int ged_draw_scene_ref_draw_state_summary(
	bsg_scene_ref ref,
	struct ged_draw_scene_draw_state_summary *out);
GED_EXPORT extern int ged_draw_scene_ref_source_snapshot(
	bsg_scene_ref ref,
	struct ged_draw_shape_source_snapshot *out);
GED_EXPORT extern void ged_draw_scene_ref_realize_dispatch(
	bsg_scene_ref ref,
	void *view_ctx);
GED_EXPORT extern int ged_draw_scene_ref_commit_database_leaf_draft(
	bsg_scene_ref parent_ref,
	struct ged *gedp,
	void *view_ctx,
	struct db_i *dbip,
	const struct db_full_path *path,
	const mat_t draw_mat,
	const struct bn_tol *tol,
	const struct bg_tess_tol *ttol,
	const struct ged_draw_appearance_settings *settings,
	int bool_op,
	const unsigned char rgb[3],
	int has_draw_size,
	fastf_t draw_size);
GED_EXPORT extern int ged_draw_scene_ref_geometry_summary(
	bsg_scene_ref ref,
	struct ged_draw_shape_geometry_summary *out);
GED_EXPORT extern int ged_draw_scene_ref_last_point(bsg_scene_ref ref,
						    point_t out);
GED_EXPORT extern int ged_draw_scene_ref_line_summary(
	bsg_scene_ref ref,
	struct ged_draw_view_line_summary *out);
GED_EXPORT extern int ged_draw_scene_ref_line_point_at(
	bsg_scene_ref ref,
	size_t index,
	point_t out);
GED_EXPORT extern int ged_draw_scene_ref_line_command_at(
	bsg_scene_ref ref,
	size_t index,
	int *out);
GED_EXPORT extern int ged_draw_scene_ref_translate_geometry(
	bsg_scene_ref ref,
	const vect_t xlate);
GED_EXPORT extern int ged_draw_scene_ref_geometry_clear(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_publish_bot_wireframe_line_set(bsg_scene_ref ref,
									const struct rt_bot_internal *bot);
GED_EXPORT extern int ged_draw_scene_ref_publish_brep_wireframe_line_set(bsg_scene_ref ref,
									 const struct rt_brep_internal *bi,
									 const struct bn_tol *tol);
GED_EXPORT extern int ged_draw_scene_ref_publish_bspline_wireframe_line_set(bsg_scene_ref ref,
									    struct rt_db_internal *ip,
									    const struct bn_tol *tol);
GED_EXPORT extern int ged_draw_scene_ref_publish_indexed_face_set(bsg_scene_ref ref,
								  const point_t *points,
								  size_t point_count,
								  const vect_t *normals,
								  size_t normal_count,
								  const int *indices,
								  size_t index_count);
GED_EXPORT extern int ged_draw_scene_ref_publish_line_set(bsg_scene_ref ref,
							  const point_t *points,
							  const int *commands,
							  size_t point_count);
GED_EXPORT extern int ged_draw_scene_ref_publish_point_set(bsg_scene_ref ref,
							   const point_t *points,
							   size_t point_count);
GED_EXPORT extern int ged_draw_scene_ref_publish_primitive_face_set(bsg_scene_ref ref,
								    struct rt_db_internal *ip,
								    const struct bg_tess_tol *ttol,
								    const struct bn_tol *tol,
								    const struct rt_view_info *view_info);
GED_EXPORT extern int ged_draw_scene_ref_publish_poly_wireframe_line_set(bsg_scene_ref ref,
									 const struct rt_pg_internal *pg);
GED_EXPORT extern int ged_draw_scene_ref_geometry_publish_nmg_region(bsg_scene_ref ref,
								     const struct nmgregion *r,
								     int style);
GED_EXPORT extern int ged_draw_scene_ref_update_bounds_from_geometry(bsg_scene_ref ref,
								     int *bad_cmd);
GED_EXPORT extern void ged_draw_scene_ref_set_draw_center(bsg_scene_ref ref,
							  const point_t center);
GED_EXPORT extern int ged_draw_scene_ref_set_bounds(bsg_scene_ref ref,
						    const point_t min,
						    const point_t max,
						    int valid);
GED_EXPORT extern int ged_draw_scene_ref_subtree_bounds(bsg_scene_ref ref,
							vect_t *min,
							vect_t *max,
							int include_overlays);
GED_EXPORT extern int ged_draw_scene_ref_update_bounds_context(bsg_scene_ref ref,
							       void *view_ctx);
GED_EXPORT extern int ged_draw_scene_ref_sync_aux_display_state(bsg_scene_ref source_ref,
								bsg_scene_ref shape_ref);
GED_EXPORT extern int ged_draw_scene_ref_set_name(bsg_scene_ref ref,
						  const char *name);
GED_EXPORT extern int ged_draw_scene_ref_set_visible(bsg_scene_ref ref,
						     int visible);
GED_EXPORT extern int ged_draw_scene_root_attach_draw_bookkeeping(struct ged *gedp,
								  bsg_scene_ref root_ref);
GED_EXPORT extern int ged_draw_scene_ref_set_transparency(bsg_scene_ref ref,
							  fastf_t transparency);
GED_EXPORT extern int ged_draw_scene_ref_set_color(bsg_scene_ref ref,
						   const unsigned char rgb[3]);
GED_EXPORT extern int ged_draw_scene_ref_set_highlighted(bsg_scene_ref ref,
							 int highlighted);
GED_EXPORT extern int ged_draw_scene_ref_set_transform(bsg_scene_ref ref,
						       const mat_t mat);
GED_EXPORT extern int ged_draw_scene_ref_set_draw_mat(bsg_scene_ref ref,
						      const mat_t mat);
GED_EXPORT extern int ged_draw_scene_ref_set_draw_size(bsg_scene_ref ref,
						       fastf_t size);
GED_EXPORT extern int ged_draw_scene_ref_set_line_style(bsg_scene_ref ref,
							int dashed);
GED_EXPORT extern int ged_draw_scene_ref_set_line_width(bsg_scene_ref ref,
							int line_width);
GED_EXPORT extern void ged_draw_scene_ref_set_draw_mode(bsg_scene_ref ref, int draw_mode);
GED_EXPORT extern int ged_draw_scene_ref_ensure_draw_intent(bsg_scene_ref ref,
							    const char *path,
							    int draw_mode);
GED_EXPORT extern int ged_draw_scene_ref_set_draw_intent_path(bsg_scene_ref ref,
							      const char *path,
							      int fallback_draw_mode);
GED_EXPORT extern int ged_draw_scene_ref_set_draw_intent_mode(bsg_scene_ref ref,
							      int draw_mode,
							      const char *fallback_path);
GED_EXPORT extern int ged_draw_scene_ref_set_draw_intent_appearance_settings(
	bsg_scene_ref ref,
	const struct ged_draw_appearance_settings *settings,
	const char *fallback_path);
GED_EXPORT extern int ged_draw_scene_ref_draw_intent_appearance_settings(
	bsg_scene_ref ref,
	struct ged_draw_appearance_settings *settings);
GED_EXPORT extern void ged_draw_scene_ref_apply_qray_work_flag(bsg_scene_ref ref, int wflag);
GED_EXPORT extern int ged_draw_scene_ref_set_legacy_color_info(
	bsg_scene_ref ref,
	const unsigned char basecolor[3],
	int user_color,
	int default_color);
GED_EXPORT extern int ged_draw_scene_ref_set_legacy_eval_flag(bsg_scene_ref ref,
							      int evaluated_region);
GED_EXPORT extern int ged_draw_scene_ref_set_legacy_uses_default_color(bsg_scene_ref ref,
								       int default_color);
GED_EXPORT extern int ged_draw_scene_ref_set_legacy_region_id(bsg_scene_ref ref,
							      int region_id);
GED_EXPORT extern int ged_draw_shape_draft_apply_db_object_marker(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_bump_appearance_revision(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_apply_overlay_geometry_attributes(bsg_scene_ref ref,
									   const unsigned char rgb[3],
									   fastf_t transparency,
									   int draw_mode);
GED_EXPORT extern int ged_draw_scene_ref_apply_display_settings(
	bsg_scene_ref ref,
	const struct ged_draw_appearance_settings *settings);
GED_EXPORT extern void ged_draw_scene_ref_set_material_rgb(bsg_scene_ref ref,
							   const unsigned char rgb[3]);
GED_EXPORT extern int ged_draw_scene_ref_realize_view_inputs_changed(
	bsg_scene_ref ref,
	void *view_ctx);
GED_EXPORT extern void ged_draw_scene_ref_invalidate(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_ensure_registry_entry(struct ged *gedp,
							       bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_apply_path_state(struct ged *gedp,
							  bsg_scene_ref ref,
							  const struct db_full_path *path);
GED_EXPORT extern int ged_draw_shape_draft_apply_registry_region(struct ged *gedp,
								 bsg_scene_ref ref,
								 int region_id,
								 int aircode,
								 int los,
								 int material_id);
GED_EXPORT extern int ged_draw_scene_ref_publish_primitive_wireframe(bsg_scene_ref ref,
								     struct rt_db_internal *ip,
								     const struct bg_tess_tol *ttol,
								     const struct bn_tol *tol,
								     void *view_ctx,
								     const struct rt_view_info *view_info,
								     int adaptive);
GED_EXPORT extern void ged_draw_scene_ref_realize(bsg_scene_ref ref, void *view_ctx);
GED_EXPORT extern const char *ged_draw_dbpath_skip_lead_slash(const char *s);

GED_EXPORT extern ged_draw_shape_state *ged_draw_shape_state_get_scene_ref(bsg_scene_ref ref);
GED_EXPORT extern ged_draw_shape_state *ged_draw_shape_state_ensure_scene_ref(struct ged *gedp,
									     bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_set_indexed_fullpath(struct ged *gedp,
							      bsg_scene_ref ref,
							      const struct db_full_path *path);

GED_EXPORT extern void *ged_draw_active_view_ctx(struct ged *gedp);
GED_EXPORT extern void ged_draw_active_view_ctx_set(struct ged *gedp,
						    void *view_ctx);
GED_EXPORT extern int ged_draw_ensure_root_attached(struct ged *gedp);
GED_EXPORT extern bsg_scene_ref ged_draw_ensure_root(struct ged *gedp);
GED_EXPORT extern bsg_scene_ref ged_draw_scene_root_ref(struct ged *gedp);
GED_EXPORT extern int ged_draw_scene_root_foreach_shape(struct ged *gedp,
							int skip_overlay_groups,
							int (*cb)(bsg_scene_ref ref, void *userdata),
							void *userdata);
GED_EXPORT extern int ged_draw_scene_ref_foreach_shape(bsg_scene_ref ref,
						       int (*cb)(bsg_scene_ref ref, void *userdata),
						       void *userdata);
GED_EXPORT extern int ged_draw_scene_root_foreach_group(struct ged *gedp,
							int (*cb)(bsg_scene_ref ref, void *userdata),
							void *userdata);
GED_EXPORT extern int ged_draw_scene_root_subtree_bounds(struct ged *gedp,
							 vect_t *min,
							 vect_t *max,
							 int pflag);
GED_EXPORT extern int ged_draw_scene_root_mutate(struct ged *gedp,
						 int (*cb)(bsg_scene_ref root_ref,
							   void *userdata),
						 void *userdata);

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
GED_EXPORT extern void ged_draw_scene_ref_release(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_redraw_group_ref(struct ged *gedp,
						ged_draw_group_ref ref,
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

GED_EXPORT extern ged_draw_shape_ref ged_draw_shape_ref_from_scene_ref(struct ged *gedp,
								       bsg_scene_ref ref);
GED_EXPORT extern ged_draw_group_ref ged_draw_group_ref_from_scene_ref(struct ged *gedp,
								       bsg_scene_ref ref);
GED_EXPORT extern bsg_scene_ref ged_draw_scene_ref_from_rt_view_ref(rt_view_scene_ref ref);
GED_EXPORT extern rt_view_scene_ref ged_draw_scene_ref_to_rt_view_ref(bsg_scene_ref ref);
GED_EXPORT extern bsg_scene_ref ged_draw_scene_ref_active_scope(struct ged *gedp,
								void *view_ctx,
								int create,
								int fallback_root);
GED_EXPORT extern void *ged_draw_scene_ref_context(bsg_scene_ref ref);
GED_EXPORT extern bsg_scene_ref ged_draw_registry_shape_scene_ref(struct ged *gedp,
								  ged_draw_shape_ref ref);
GED_EXPORT extern bsg_scene_ref ged_draw_registry_group_scene_ref(struct ged *gedp,
								  ged_draw_group_ref ref);
GED_EXPORT extern bsg_scene_ref ged_draw_shape_scene_ref_from_cache_ref(struct ged *gedp,
									ged_draw_shape_ref ref);

typedef int (*ged_draw_scene_ref_index_cb)(bsg_scene_ref ref, void *userdata);

GED_EXPORT extern int ged_draw_shape_index_for_component(struct ged *gedp,
							const char *path,
							ged_draw_scene_ref_index_cb cb,
							void *userdata);
GED_EXPORT extern int ged_draw_group_index_for_component(struct ged *gedp,
							const char *path,
							ged_draw_scene_ref_index_cb cb,
							void *userdata);
GED_EXPORT extern int ged_draw_shape_index_for_path_hash(struct ged *gedp,
							 unsigned long long path_hash,
							 ged_draw_scene_ref_index_cb cb,
							 void *userdata);
struct ged_draw_index_stats {
    uint64_t shape_component_queries;
    uint64_t shape_component_candidates;
    uint64_t group_component_queries;
    uint64_t group_component_candidates;
    uint64_t path_queries;
    uint64_t path_candidates;
    uint64_t fallback_shape_scans;
    uint64_t fallback_group_scans;
};

GED_EXPORT extern void ged_draw_index_stats_get(struct ged *gedp,
						struct ged_draw_index_stats *stats);
GED_EXPORT extern void ged_draw_index_stats_reset(struct ged *gedp);

GED_EXPORT extern bsg_scene_ref ged_draw_first_shape_scene_ref(struct ged *gedp);
GED_EXPORT extern int ged_draw_group_ref_set_dbpath(struct ged *gedp,
						    ged_draw_group_ref ref,
						    const struct db_full_path *new_dfp);
GED_EXPORT extern int ged_draw_group_ref_set_mode(struct ged *gedp,
						  ged_draw_group_ref ref,
						  int mode);
GED_EXPORT extern int ged_draw_group_ref_set_appearance_settings(struct ged *gedp,
								ged_draw_group_ref ref,
								const struct ged_draw_appearance_settings *settings);
GED_EXPORT extern int ged_draw_shape_ref_apply_qray_work_flag(struct ged *gedp,
							      ged_draw_shape_ref ref,
							      int wflag);
GED_EXPORT extern int ged_draw_shape_ref_release(struct ged *gedp,
						 ged_draw_shape_ref ref);
GED_EXPORT extern int ged_draw_shape_ref_realize_context(struct ged *gedp,
							 ged_draw_shape_ref ref,
							 void *view_ctx);
GED_EXPORT extern int ged_draw_shape_ref_set_visible(struct ged *gedp,
						     ged_draw_shape_ref ref,
						     int visible);
GED_EXPORT extern int ged_draw_shape_ref_get_color(struct ged *gedp,
					   ged_draw_shape_ref ref,
					   unsigned char rgb[3]);
GED_EXPORT extern int ged_draw_shape_ref_set_color(struct ged *gedp,
					   ged_draw_shape_ref ref,
					   const unsigned char rgb[3]);
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
GED_EXPORT extern int ged_draw_group_ref_append_scene_ref(struct ged *gedp,
							  ged_draw_group_ref ref,
							  bsg_scene_ref shape_ref);

GED_EXPORT extern void ged_draw_scene_ref_highlight_free_cb(bsg_scene_ref ref);

__END_DECLS

#endif /* LIBGED_BSG_GED_DRAW_PRIVATE_H */
