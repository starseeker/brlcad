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

typedef enum {
    GED_DRAW_DATABASE_EVENT_MODIFIED = 0,
    GED_DRAW_DATABASE_EVENT_REMOVED = 1,
    GED_DRAW_DATABASE_EVENT_RENAMED = 2
} ged_draw_database_event_kind;

struct ged_draw_database_event_record {
    const char *path;
    const char *old_path;
    ged_draw_database_event_kind kind;
};

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

GED_EXPORT extern void ged_draw_source_data_free(void *data);
GED_EXPORT extern void ged_draw_scene_ref_realization_set_view_context_policy(bsg_scene_ref ref,
									      const void *view_ctx);
GED_EXPORT extern int ged_draw_brep_mesh_lod_detail_setup(struct rt_mesh_lod *lod,
							  struct db_i *dbip,
							  struct directory *dp,
							  const struct bg_tess_tol *ttol,
							  const struct bn_tol *tol);
GED_EXPORT extern bsg_scene_ref ged_draw_shape_source_ref(bsg_scene_ref ref);
GED_EXPORT extern void ged_draw_scene_ref_set_source_ref(bsg_scene_ref ref,
							 bsg_scene_ref source_ref);
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
GED_EXPORT extern int ged_draw_scene_ref_has_database_source(
	bsg_scene_ref ref);
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
GED_EXPORT extern int ged_draw_scene_ref_source_ensure(bsg_scene_ref ref);
GED_EXPORT extern void ged_draw_scene_ref_database_source_sync(
	bsg_scene_ref ref,
	const struct ged_draw_source_state *source_data,
	const ged_draw_shape_state *shape_data);
GED_EXPORT extern int ged_draw_database_source_record_has_state(
	const struct ged_draw_database_source_record *record);
GED_EXPORT extern int ged_draw_database_source_record_is_stale(
	const struct ged_draw_database_source_record *record);
GED_EXPORT extern int ged_draw_scene_ref_database_source_record_get(
	bsg_scene_ref ref,
	struct ged_draw_database_source_record *record);
GED_EXPORT extern int ged_draw_scene_ref_database_source_record_apply(
	bsg_scene_ref ref,
	const struct ged_draw_database_source_record *record);
GED_EXPORT extern void ged_draw_log(int level, const char *fmt, ...) _BU_ATTR_PRINTF23;
GED_EXPORT extern void ged_draw_scene_ref_mark_realization_stale(bsg_scene_ref ref,
								 uint64_t source_revision,
								 uint64_t inputs_revision,
								 ged_draw_stale_reason reason);
GED_EXPORT extern void ged_draw_scene_ref_mark_realization_result(bsg_scene_ref ref,
								  uint64_t source_revision,
								  uint64_t inputs_revision,
								  int failed);
GED_EXPORT extern int ged_draw_source_primitive_has_lod_realize(const struct rt_db_internal *ip);
GED_EXPORT extern void ged_draw_test_force_primitive_face_set_failure(int enable);

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
GED_EXPORT extern int ged_draw_shape_draft_set_fullpath(ged_draw_shape_draft *draft,
							const struct db_full_path *path);
GED_EXPORT extern int ged_draw_shape_draft_set_region(ged_draw_shape_draft *draft,
						      int region_id,
						      int aircode,
						      int los,
						      int material_id);
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
GED_EXPORT extern int ged_draw_shape_draft_clear_geometry(ged_draw_shape_draft *draft);
GED_EXPORT extern int ged_draw_shape_draft_update_scene_bounds(ged_draw_shape_draft *draft);
GED_EXPORT extern int ged_draw_shape_draft_update_bounds_from_geometry(ged_draw_shape_draft *draft,
								       int *bad_cmd);
GED_EXPORT extern int ged_draw_shape_draft_set_bounds_from_minmax(ged_draw_shape_draft *draft,
								  const point_t min,
								  const point_t max,
								  int set_node_bounds);
GED_EXPORT extern int ged_draw_shape_draft_set_center_size(ged_draw_shape_draft *draft,
							   const point_t center,
							   fastf_t size);
GED_EXPORT extern int ged_draw_shape_draft_publish_primitive_wireframe(ged_draw_shape_draft *draft,
								       struct rt_db_internal *ip,
								       const struct bg_tess_tol *ttol,
								       const struct bn_tol *tol,
								       void *view_ctx,
								       int adaptive);
GED_EXPORT extern int ged_draw_shape_draft_set_highlighted(ged_draw_shape_draft *draft,
							   int highlighted);
GED_EXPORT extern int ged_draw_shape_draft_set_line_style(ged_draw_shape_draft *draft,
							  int dashed);
GED_EXPORT extern int ged_draw_shape_draft_set_line_width(ged_draw_shape_draft *draft,
							  int line_width);
GED_EXPORT extern int ged_draw_shape_draft_set_legacy_color_info(ged_draw_shape_draft *draft,
								 const unsigned char basecolor[3],
								 int user_color,
								 int default_color);
GED_EXPORT extern int ged_draw_shape_draft_set_legacy_uses_default_color(ged_draw_shape_draft *draft,
									 int default_color);
GED_EXPORT extern int ged_draw_shape_draft_set_legacy_eval_flag(ged_draw_shape_draft *draft,
								int eflag);
GED_EXPORT extern int ged_draw_shape_draft_set_legacy_region_id(ged_draw_shape_draft *draft,
								int region_id);
GED_EXPORT extern int ged_draw_shape_draft_color_soltab(ged_draw_shape_draft *draft,
							struct db_i *dbip);
GED_EXPORT extern void color_soltab_scene_ref(struct db_i *dbip, bsg_scene_ref shape_ref);
GED_EXPORT extern int ged_draw_shape_draft_set_name(ged_draw_shape_draft *draft,
						    const char *name);
GED_EXPORT extern int ged_draw_shape_draft_set_source_state(ged_draw_shape_draft *draft,
							    struct ged_draw_source_state *data);
GED_EXPORT extern int ged_draw_shape_draft_mark_db_object(ged_draw_shape_draft *draft);
GED_EXPORT extern int ged_draw_shape_draft_apply_appearance_settings(ged_draw_shape_draft *draft,
								     const struct ged_draw_appearance_settings *settings);
GED_EXPORT extern int ged_draw_shape_draft_bump_appearance_revision(ged_draw_shape_draft *draft);
GED_EXPORT extern int ged_draw_shape_draft_set_material_rgb(ged_draw_shape_draft *draft,
							    unsigned char r,
							    unsigned char g,
							    unsigned char b);
GED_EXPORT extern int ged_draw_shape_draft_set_transform(ged_draw_shape_draft *draft,
							 const mat_t mat);
GED_EXPORT extern int ged_draw_shape_draft_set_bounds(ged_draw_shape_draft *draft,
						      const point_t min,
						      const point_t max);
GED_EXPORT extern int ged_draw_shape_draft_set_draw_size(ged_draw_shape_draft *draft,
							 fastf_t size);
GED_EXPORT extern int ged_draw_shape_draft_set_visible(ged_draw_shape_draft *draft,
						       int visible);
GED_EXPORT extern int ged_draw_shape_draft_set_transparency(ged_draw_shape_draft *draft,
							    fastf_t transparency);
GED_EXPORT extern int ged_draw_shape_draft_set_draw_mode(ged_draw_shape_draft *draft,
							 int draw_mode);
GED_EXPORT extern int ged_draw_shape_draft_set_draw_mat(ged_draw_shape_draft *draft,
							const mat_t mat);
GED_EXPORT extern ged_draw_shape_ref ged_draw_shape_draft_commit_to_group(ged_draw_shape_draft *draft,
									  ged_draw_group_ref group_ref);
GED_EXPORT extern ged_draw_shape_ref ged_draw_shape_draft_commit_to_last_group(ged_draw_shape_draft *draft);
GED_EXPORT extern bsg_scene_ref ged_draw_shape_draft_commit_to_scene_ref(ged_draw_shape_draft *draft,
									 bsg_scene_ref parent_ref);

GED_EXPORT extern const char *ged_draw_scene_ref_name(bsg_scene_ref ref);
GED_EXPORT extern const struct db_full_path *ged_draw_scene_ref_fullpath(bsg_scene_ref ref);
GED_EXPORT extern struct directory *ged_draw_scene_ref_leaf_dp(bsg_scene_ref ref);
GED_EXPORT extern ged_draw_shape_state *ged_draw_scene_ref_shape_state(bsg_scene_ref ref);
GED_EXPORT extern struct ged_draw_source_state *ged_draw_scene_ref_source_data(bsg_scene_ref ref);
GED_EXPORT extern void ged_draw_scene_ref_source_data_set(bsg_scene_ref ref,
							  struct ged_draw_source_state *data);
GED_EXPORT extern struct bu_list *ged_draw_scene_ref_geometry_pool(bsg_scene_ref ref);
GED_EXPORT extern bsg_scene_ref ged_draw_scene_ref_null(void);
GED_EXPORT extern bsg_scene_ref ged_draw_scene_ref_from_context(void *scene_ctx);
GED_EXPORT extern int ged_draw_scene_ref_is_null(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_equal(bsg_scene_ref a, bsg_scene_ref b);
GED_EXPORT extern int ged_draw_scene_ref_is_group(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_is_shape(bsg_scene_ref ref);
GED_EXPORT extern bsg_scene_ref ged_draw_scene_ref_parent(bsg_scene_ref ref);
GED_EXPORT extern size_t ged_draw_scene_ref_child_count(bsg_scene_ref ref);
GED_EXPORT extern bsg_scene_ref ged_draw_scene_ref_child_at(bsg_scene_ref ref,
							    size_t idx);
GED_EXPORT extern int ged_draw_scene_ref_visible(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_highlighted(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_legacy_eval_flag(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_legacy_region_id(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_legacy_basecolor(bsg_scene_ref ref,
							  unsigned char rgb[3]);
GED_EXPORT extern int ged_draw_scene_ref_legacy_user_color(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_legacy_default_color(bsg_scene_ref ref);
GED_EXPORT extern uint64_t ged_draw_scene_ref_drawn_revision(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_changed(bsg_scene_ref ref);
GED_EXPORT extern fastf_t ged_draw_scene_ref_transparency(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_line_width(bsg_scene_ref ref);
GED_EXPORT extern fastf_t ged_draw_scene_ref_draw_size(bsg_scene_ref ref);
GED_EXPORT extern void ged_draw_scene_ref_draw_center(bsg_scene_ref ref,
						      point_t center);
GED_EXPORT extern int ged_draw_scene_ref_color(bsg_scene_ref ref,
					       unsigned char rgb[3]);
GED_EXPORT extern int ged_draw_scene_ref_transform(bsg_scene_ref ref,
						   mat_t mat);
GED_EXPORT extern int ged_draw_scene_ref_is_view_scope(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_is_view_source(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_is_local_source(bsg_scene_ref ref);
GED_EXPORT extern void *ged_draw_scene_ref_draw_context_owner(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_display_summary(
	bsg_scene_ref ref,
	struct ged_draw_scene_display_summary *out);
GED_EXPORT extern int ged_draw_scene_ref_source_summary(
	bsg_scene_ref ref,
	struct ged_draw_database_source_summary *out);
GED_EXPORT extern int ged_draw_scene_ref_tree_summary(
	bsg_scene_ref ref,
	struct ged_draw_scene_tree_summary *out);
GED_EXPORT extern int ged_draw_scene_ref_material_summary(
	bsg_scene_ref ref,
	struct ged_draw_shape_material_summary *out);
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
GED_EXPORT extern int ged_draw_scene_ref_update_indexed_face_set(bsg_scene_ref ref,
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
GED_EXPORT extern int ged_draw_scene_ref_geometry_publish_nmg_model(bsg_scene_ref ref,
								    const struct model *m,
								    int style);
GED_EXPORT extern int ged_draw_scene_ref_update_bounds_from_geometry(bsg_scene_ref ref,
								     int *bad_cmd);
GED_EXPORT extern void ged_draw_scene_ref_set_draw_center(bsg_scene_ref ref,
							  const point_t center);
GED_EXPORT extern int ged_draw_scene_ref_bounds(bsg_scene_ref ref,
						point_t min,
						point_t max);
GED_EXPORT extern int ged_draw_scene_ref_set_bounds(bsg_scene_ref ref,
						    const point_t min,
						    const point_t max,
						    int valid);
GED_EXPORT extern int ged_draw_scene_ref_set_bounds_from_minmax(bsg_scene_ref ref,
								const point_t min,
								const point_t max,
								int set_scene_bounds);
GED_EXPORT extern int ged_draw_scene_ref_subtree_bounds(bsg_scene_ref ref,
							vect_t *min,
							vect_t *max,
							int include_overlays);
GED_EXPORT extern int ged_draw_scene_ref_draw_tree_depth(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_update_bounds_context(bsg_scene_ref ref,
							       void *view_ctx);
GED_EXPORT extern void *ged_draw_scene_ref_view_context(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_copy_aux_display_state(bsg_scene_ref dst,
								bsg_scene_ref src);
GED_EXPORT extern int ged_draw_scene_ref_set_name(bsg_scene_ref ref,
						  const char *name);
GED_EXPORT extern int ged_draw_scene_ref_set_visible(bsg_scene_ref ref,
						     int visible);
GED_EXPORT extern int ged_draw_scene_ref_set_draw_context(bsg_scene_ref ref,
							  void *draw_ctx);
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
GED_EXPORT extern int ged_draw_scene_ref_draw_mat(bsg_scene_ref ref, mat_t mat);
GED_EXPORT extern int ged_draw_scene_ref_draw_mode(bsg_scene_ref ref);
GED_EXPORT extern void ged_draw_scene_ref_set_draw_mode(bsg_scene_ref ref, int draw_mode);
GED_EXPORT extern int ged_draw_scene_ref_set_changed(bsg_scene_ref ref,
						     int changed);
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
GED_EXPORT extern int ged_draw_scene_ref_draw_intent_is_overlay(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_revalidate_draw_intent(
	bsg_scene_ref ref,
	const struct ged_draw_database_event_record *event);
GED_EXPORT extern int ged_draw_scene_ref_draw_intent_appearance_settings(
	bsg_scene_ref ref,
	struct ged_draw_appearance_settings *settings);
GED_EXPORT extern const char *ged_draw_scene_ref_draw_intent_path(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_draw_intent_mode(bsg_scene_ref ref,
							 int *draw_mode);
GED_EXPORT extern int ged_draw_scene_ref_copy_draw_intent(bsg_scene_ref dst,
							  bsg_scene_ref src);
GED_EXPORT extern int ged_draw_scene_ref_line_style(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_strict_fallback(bsg_scene_ref ref);
GED_EXPORT extern void ged_draw_scene_ref_set_work_flag(bsg_scene_ref ref, int wflag);
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
GED_EXPORT extern int ged_draw_scene_ref_mark_db_object(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_set_non_database_source(
	bsg_scene_ref ref,
	int non_database_source);
GED_EXPORT extern int ged_draw_scene_ref_bump_changed(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_apply_overlay_geometry_attributes(bsg_scene_ref ref,
									   const unsigned char rgb[3],
									   fastf_t transparency,
									   int draw_mode);
GED_EXPORT extern int ged_draw_scene_ref_apply_appearance_settings(
	bsg_scene_ref ref,
	const struct ged_draw_appearance_settings *settings);
GED_EXPORT extern void ged_draw_scene_ref_set_material_rgb(bsg_scene_ref ref,
							   const unsigned char rgb[3]);
GED_EXPORT extern uint64_t ged_draw_scene_ref_material_revision(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_set_material_revision(bsg_scene_ref ref,
							       uint64_t mater_rev);
GED_EXPORT extern int ged_draw_scene_ref_realization_current(bsg_scene_ref ref);
GED_EXPORT extern void ged_draw_scene_ref_realization_set_current(bsg_scene_ref ref,
								  int current);
GED_EXPORT extern void ged_draw_scene_ref_realization_set_roles(bsg_scene_ref ref,
								int csg_obj,
								int mesh_obj);
GED_EXPORT extern fastf_t ged_draw_scene_ref_realization_view_scale(bsg_scene_ref ref);
GED_EXPORT extern fastf_t ged_draw_scene_ref_realization_curve_scale(bsg_scene_ref ref);
GED_EXPORT extern fastf_t ged_draw_scene_ref_realization_point_scale(bsg_scene_ref ref);
GED_EXPORT extern void ged_draw_scene_ref_realization_set_view_policy(bsg_scene_ref ref,
								      int view_dependent,
								      fastf_t view_scale,
								      size_t bot_threshold,
								      fastf_t curve_scale,
								      fastf_t point_scale);
GED_EXPORT extern int ged_draw_scene_ref_realization_pipeline_candidate(bsg_scene_ref ref);
GED_EXPORT extern void ged_draw_scene_ref_realization_reset(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_mark_view_inputs_changed(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_mark_realized_current(bsg_scene_ref ref);
GED_EXPORT extern void ged_draw_scene_ref_invalidate(bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_prepare(struct ged *gedp,
						 bsg_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_set_fullpath(struct ged *gedp,
						      bsg_scene_ref ref,
						      const struct db_full_path *path);
GED_EXPORT extern int ged_draw_scene_ref_set_region(struct ged *gedp,
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
GED_EXPORT extern void ged_draw_shape_state_release_scene_ref(bsg_scene_ref ref);
GED_EXPORT extern ged_draw_shape_state *ged_draw_shape_ref_set_fullpath(bsg_scene_ref ref,
						    struct ged *gedp,
						    const struct db_full_path *path);
GED_EXPORT extern void ged_draw_shape_state_set_fullpath(ged_draw_shape_state *data,
							const struct db_full_path *path);
GED_EXPORT extern void ged_draw_scene_ref_index_remove(struct ged *gedp,
						       bsg_scene_ref ref);
GED_EXPORT extern void ged_draw_scene_ref_index_add(struct ged *gedp,
						    bsg_scene_ref ref);
GED_EXPORT extern void ged_draw_shape_state_set_region(ged_draw_shape_state *data,
						      int region_id,
						      int aircode,
						      int los,
						      int material_id);

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
GED_EXPORT extern int ged_draw_scene_root_foreach_group(struct ged *gedp,
							int (*cb)(bsg_scene_ref ref, void *userdata),
							void *userdata);
GED_EXPORT extern int ged_draw_scene_root_subtree_bounds(struct ged *gedp,
							 vect_t *min,
							 vect_t *max,
							 int pflag);
GED_EXPORT extern bsg_scene_ref ged_draw_scene_root_last_group_ref(struct ged *gedp);
GED_EXPORT extern int ged_draw_scene_root_has_groups(struct ged *gedp);
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
GED_EXPORT extern void ged_draw_revalidate_db_event(struct ged *gedp,
						    const struct ged_draw_database_event_record *ev);
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
GED_EXPORT extern size_t ged_draw_scene_ref_child_count_rt(rt_view_scene_ref ref);
GED_EXPORT extern int ged_draw_scene_ref_revalidate_draw_intent_rt(rt_view_scene_ref ref,
								   const struct ged_draw_database_event_record *event);
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
GED_EXPORT extern int ged_draw_group_index_for_path_hash(struct ged *gedp,
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
GED_EXPORT extern const char *ged_draw_group_scene_ref_path(bsg_scene_ref group_ref);
GED_EXPORT extern int ged_draw_group_scene_ref_mode(bsg_scene_ref group_ref);
GED_EXPORT extern int ged_draw_group_scene_ref_dbpath(struct ged *gedp,
						      bsg_scene_ref group_ref,
						      struct db_full_path *out);
GED_EXPORT extern int ged_draw_group_ref_set_dbpath(struct ged *gedp,
						    ged_draw_group_ref ref,
						    const struct db_full_path *new_dfp);
GED_EXPORT extern int ged_draw_group_ref_set_mode(struct ged *gedp,
						  ged_draw_group_ref ref,
						  int mode);
GED_EXPORT extern int ged_draw_group_ref_set_appearance_settings(struct ged *gedp,
								ged_draw_group_ref ref,
								const struct ged_draw_appearance_settings *settings);
GED_EXPORT extern int ged_draw_shape_ref_set_flag(struct ged *gedp,
						  ged_draw_shape_ref ref,
						  int flag);
GED_EXPORT extern int ged_draw_shape_ref_set_work_flag(struct ged *gedp,
						       ged_draw_shape_ref ref,
						       int wflag);
GED_EXPORT extern int ged_draw_shape_ref_line_style(struct ged *gedp,
						    ged_draw_shape_ref ref);
GED_EXPORT extern int ged_draw_shape_ref_release(struct ged *gedp,
						 ged_draw_shape_ref ref);
GED_EXPORT extern int ged_draw_shape_ref_realize_context(struct ged *gedp,
							 ged_draw_shape_ref ref,
							 void *view_ctx);
GED_EXPORT extern int ged_draw_shape_ref_reset_node(struct ged *gedp,
						    ged_draw_shape_ref ref);
GED_EXPORT extern int ged_draw_shape_ref_set_visible(struct ged *gedp,
						     ged_draw_shape_ref ref,
						     int visible);
GED_EXPORT extern int ged_draw_shape_ref_get_color(struct ged *gedp,
					   ged_draw_shape_ref ref,
					   unsigned char rgb[3]);
GED_EXPORT extern int ged_draw_shape_ref_set_color(struct ged *gedp,
					   ged_draw_shape_ref ref,
					   const unsigned char rgb[3]);
GED_EXPORT extern int ged_draw_shape_ref_refresh_material(struct ged *gedp,
							  ged_draw_shape_ref ref,
							  unsigned char rgb[3]);
GED_EXPORT extern int ged_draw_shape_ref_lod_ensure(struct ged *gedp,
						    ged_draw_shape_ref ref,
						    void *first_view_ctx,
						    void **view_ctxs,
						    size_t view_ctx_count);
GED_EXPORT extern int ged_draw_shape_ref_pipeline_candidate(struct ged *gedp,
							    ged_draw_shape_ref ref);
GED_EXPORT extern void *ged_draw_shape_ref_view_context(struct ged *gedp,
							ged_draw_shape_ref ref);
GED_EXPORT extern int ged_draw_view_selection_add_shape_ref_context(
	struct ged *gedp,
	void *view_ctx,
	ged_draw_shape_ref ref,
	void **selection_view_ctx,
	struct bu_vls *path);
GED_EXPORT extern void ged_draw_append_scene_ref_to_last_group(struct ged *gedp,
							       bsg_scene_ref shape_ref);
GED_EXPORT extern int ged_draw_group_scene_ref_is_overlay(bsg_scene_ref group_ref);
GED_EXPORT extern int ged_draw_group_ref_append_scene_ref(struct ged *gedp,
							  ged_draw_group_ref ref,
							  bsg_scene_ref shape_ref);

GED_EXPORT extern void ged_draw_scene_ref_highlight_free_cb(bsg_scene_ref ref);

__END_DECLS

#endif /* LIBGED_BSG_GED_DRAW_PRIVATE_H */
