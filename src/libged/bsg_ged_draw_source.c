/*              B S G _ G E D _ D R A W _ S O U R C E . C
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this file; see the file named COPYING for more information.
 */
/** @file bsg_ged_draw_source.c
 *
 * GED draw-source metadata helpers backed by public BSG database-source fields.
 */

#include "common.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "bg/plane.h"
#include "bg/sat.h"
#include "bn/mat.h"
#include "bu/hash.h"
#include "bu/log.h"
#include "bu/malloc.h"
#include "bu/parallel.h"
#include "bu/path.h"
#include "bu/str.h"
#include "bsg/appearance.h"
#include "bsg/database_source.h"
#include "bsg/draw_ctx.h"
#include "bsg/draw_set.h"
#include "bsg/draw_source.h"
#include "bsg/draw_intent.h"
#include "bsg/export.h"
#include "bsg/field.h"
#include "bsg/geometry.h"
#include "bsg/material.h"
#include "bsg/render.h"
#include "bsg/scene_builder.h"
#include "ged/draw.h"
#include "ged/db_index.h"
#include "ged/view.h"
#include "nmg.h"
#include "rt/db_instance.h"
#include "rt/db_io.h"
#include "rt/func.h"
#include "rt/functab.h"
#include "rt/global.h"
#include "rt/geom.h"
#include "rt/vlist.h"
#include "rt/view.h"
#include "../librt/librt_private.h"
#include "rt/primitives/arb8.h"
#include "rt/primitives/brep.h"
#include "rt/primitives/bspline.h"
#include "rt/primitives/datum.h"
#include "rt/primitives/dsp.h"
#include "rt/primitives/ebm.h"
#include "rt/primitives/ell.h"
#include "rt/primitives/ehy.h"
#include "rt/primitives/epa.h"
#include "rt/primitives/eto.h"
#include "rt/primitives/extrude.h"
#include "rt/primitives/hf.h"
#include "rt/primitives/hrt.h"
#include "rt/primitives/pipe.h"
#include "rt/primitives/rhc.h"
#include "rt/primitives/rpc.h"
#include "rt/primitives/revolve.h"
#include "rt/primitives/sketch.h"
#include "rt/primitives/superell.h"
#include "rt/primitives/tgc.h"
#include "rt/primitives/tor.h"
#include "rt/primitives/vol.h"
#include "rt/tree.h"
#include "./ged_private.h"
#include "./bsg_ged_draw_private.h"
#include "./bsg_ged_draw_shape_state_private.h"
#include "./bsg_ged_draw_view_private.h"


struct ged_draw_source_state;
struct ged_draw_source_state_record;
struct ged_draw_source_runtime_summary;
struct ged_draw_database_source_record;
struct ged_draw_scene_draw_state_summary;

static void _ged_draw_scene_ref_geometry_set_command_count(bsg_scene_ref ref, size_t vlen);
static bsg_scene_ref ged_draw_view_context_database_source_create(void *view_ctx,
								  const char *name);
static bsg_scene_ref ged_draw_view_context_geometry_create(void *view_ctx,
							   const char *name);
static int ged_draw_scene_ref_detach(bsg_scene_ref ref);
static int ged_draw_scene_ref_append_child(bsg_scene_ref parent_ref,
					   bsg_scene_ref child_ref);
static int ged_draw_scene_ref_foreach_child(
	bsg_scene_ref ref,
	int (*cb)(bsg_scene_ref child_ref, void *userdata),
	void *userdata);
static bsg_scene_ref ged_draw_view_context_group_child_ensure(
	void *view_ctx,
	bsg_scene_ref parent_ref,
	const char *name,
	void *dp_hint);
static int ged_draw_scene_ref_is_database_source(bsg_scene_ref ref);
static void ged_draw_scene_ref_release(bsg_scene_ref ref);
static bsg_scene_ref ged_draw_scene_ref_null(void);
static bsg_scene_ref ged_draw_scene_ref_from_context(void *scene_ctx);
static void *ged_draw_scene_ref_context(bsg_scene_ref ref);
static rt_view_scene_ref ged_draw_scene_ref_to_rt_view_ref(bsg_scene_ref ref);
static bsg_scene_ref ged_draw_scene_ref_from_rt_view_ref(rt_view_scene_ref ref);
static int ged_draw_scene_ref_is_null(bsg_scene_ref ref);
static int ged_draw_scene_ref_equal(bsg_scene_ref a, bsg_scene_ref b);
static ged_draw_shape_state *_ged_draw_shape_state_get_scene_ref(
	bsg_scene_ref ref);
static ged_draw_shape_ref _ged_draw_shape_ref_from_scene_ref(
	struct ged *gedp,
	bsg_scene_ref ref);
static ged_draw_group_ref _ged_draw_group_ref_from_scene_ref(
	struct ged *gedp,
	bsg_scene_ref ref);
static bsg_scene_ref _ged_draw_shape_ref_scene_ref(struct ged *gedp,
							   ged_draw_shape_ref ref);
static bsg_scene_ref _ged_draw_shape_ref_runtime_scene_ref(struct ged *gedp,
								   ged_draw_shape_ref ref);
static int _ged_draw_shape_ref_release(struct ged *gedp, ged_draw_shape_ref ref);
static bsg_scene_ref _ged_draw_group_ref_scene_ref(struct ged *gedp,
							   ged_draw_group_ref ref);
static int ged_draw_scene_ref_obol_group_ensure_path(
	bsg_scene_ref ref,
	const char *preferred_path,
	int mode,
	int overlay);
static int ged_draw_shape_ref_lod_ensure_obol(struct ged *gedp,
					      ged_draw_shape_ref ref,
					      void *first_view_ctx);
static int ged_draw_obol_scene_controller_is_owned(struct ged *gedp);
typedef int (*ged_draw_shape_ref_obol_path_cb)(struct ged *gedp,
					       const char *path,
					       void *userdata);
typedef int (*ged_draw_group_ref_obol_path_cb)(struct ged *gedp,
					       const char *path,
					       void *userdata);
typedef int (*ged_draw_scene_ref_obol_source_path_cb)(struct ged *gedp,
						      const char *path,
						      void *data);
static int ged_draw_scene_ref_obol_source_path_apply(
	bsg_scene_ref ref,
	ged_draw_scene_ref_obol_source_path_cb cb,
	void *data);
static int ged_draw_scene_ref_obol_source_path_apply_ensure(
	bsg_scene_ref ref,
	ged_draw_scene_ref_obol_source_path_cb cb,
	void *data);
static int _ged_draw_scene_ref_obol_source_snapshot_path_apply(
	bsg_scene_ref ref,
	ged_draw_scene_ref_obol_source_path_cb cb,
	void *data);
static int ged_draw_obol_database_source_needs_primary_shape_proxy(
	struct ged *gedp,
	const char *path);
static int ged_draw_obol_database_source_primary_shape_info(
	struct ged *gedp,
	const char *path,
	int draw_tree_depth,
	struct ged_draw_obol_scene_context_info *out);
static int ged_draw_obol_database_source_container_info(
	struct ged *gedp,
	const char *path,
	int draw_tree_depth,
	struct ged_draw_obol_scene_context_info *out);
static char *_ged_draw_scene_ref_owner_path_string(bsg_scene_ref ref);
static int _ged_draw_scene_ref_obol_ensure_source(bsg_scene_ref primary_ref);
static int _ged_draw_scene_ref_obol_ensure_owner_source(
	bsg_scene_ref primary_ref);
static ged_draw_group_ref _ged_draw_active_scope_group_ref(struct ged *gedp,
							   void *view_ctx,
							   int create,
							   int fallback_root);
static size_t ged_draw_scene_ref_child_count(bsg_scene_ref ref);
static bsg_scene_ref ged_draw_scene_ref_child_at(bsg_scene_ref ref, size_t idx);
static bsg_scene_ref ged_draw_scene_ref_shape_owning_group(bsg_scene_ref ref);
static int ged_draw_scene_ref_is_group(bsg_scene_ref ref);
static int ged_draw_scene_ref_is_view_scope(bsg_scene_ref ref);
static int ged_draw_scene_ref_draw_tree_depth(bsg_scene_ref ref);
static int ged_draw_scene_ref_draw_mode(bsg_scene_ref ref);
static int ged_draw_scene_ref_visible(bsg_scene_ref ref);
static void ged_draw_scene_ref_set_draw_center(bsg_scene_ref ref,
					       const point_t center);
static int ged_draw_scene_ref_set_bounds(bsg_scene_ref ref,
					 const point_t min,
					 const point_t max,
					 int valid);
static int ged_draw_scene_ref_set_visible(bsg_scene_ref ref, int visible);
static int ged_draw_scene_ref_set_transparency(bsg_scene_ref ref,
					       fastf_t transparency);
static int ged_draw_scene_ref_set_color(bsg_scene_ref ref,
					const unsigned char rgb[3]);
static int ged_draw_scene_ref_set_highlighted(bsg_scene_ref ref,
					      int highlighted);
static int ged_draw_scene_ref_set_transform(bsg_scene_ref ref,
					    const mat_t mat);
static int ged_draw_scene_ref_set_draw_mat(bsg_scene_ref ref,
					   const mat_t mat);
static int ged_draw_scene_ref_set_draw_size(bsg_scene_ref ref, fastf_t size);
static int ged_draw_scene_ref_set_line_style(bsg_scene_ref ref, int dashed);
static int ged_draw_scene_ref_set_line_width(bsg_scene_ref ref, int line_width);
static void ged_draw_scene_ref_set_draw_mode(bsg_scene_ref ref, int draw_mode);
static int ged_draw_scene_ref_ensure_draw_intent(bsg_scene_ref ref,
						 const char *path,
						 int draw_mode);
static int ged_draw_scene_ref_set_draw_intent_path(bsg_scene_ref ref,
						   const char *path,
						   int fallback_draw_mode);
static int ged_draw_scene_ref_set_draw_intent_mode(bsg_scene_ref ref,
						   int draw_mode,
						   const char *fallback_path);
static int ged_draw_scene_ref_draw_intent_appearance_settings(
	bsg_scene_ref ref,
	struct ged_draw_appearance_settings *settings);
static int ged_draw_scene_ref_set_legacy_color_info(
	bsg_scene_ref ref,
	const unsigned char basecolor[3],
	int user_color,
	int default_color);
static int ged_draw_scene_ref_set_legacy_eval_flag(bsg_scene_ref ref,
						   int evaluated_region);
static int ged_draw_scene_ref_legacy_region_id(bsg_scene_ref ref);
static int ged_draw_scene_ref_set_legacy_region_id(bsg_scene_ref ref,
						   int region_id);
static void ged_draw_scene_ref_set_material_rgb(bsg_scene_ref ref,
						const unsigned char rgb[3]);
static void ged_draw_scene_ref_invalidate(bsg_scene_ref ref);
static int ged_draw_scene_ref_apply_path_state(struct ged *gedp,
					       bsg_scene_ref ref,
					       const struct db_full_path *path);
static bsg_scene_ref ged_draw_scene_ref_parent(bsg_scene_ref ref);
static bsg_scene_ref ged_draw_scene_ref_source_owner(bsg_scene_ref ref);
static int ged_draw_scene_ref_erase_nested_db_subpath(
	bsg_scene_ref parent_ref,
	const struct db_full_path *subpath,
	size_t depth_start);
static int ged_draw_scene_ref_erase_matching_group_path_or_nested(
	struct ged *gedp,
	const char *path,
	bsg_scene_ref base_ref,
	const struct db_full_path *subpath,
	void *view_ctx,
	int mode,
	int apply_scope);
static int ged_draw_scene_ref_erase_groups_by_db_subpath(
	struct ged *gedp,
	bsg_scene_ref base_ref,
	const struct db_full_path *subpath,
	void *view_ctx,
	int mode,
	int apply_scope);
static int ged_draw_scene_ref_erase_subgroups_by_name(
	bsg_scene_ref parent_ref,
	const char *name);
static int ged_draw_scene_ref_erase_groups_by_path_prefix_string(
	struct ged *gedp,
	const char *path,
	bsg_scene_ref base_ref,
	void *view_ctx,
	int mode);
static int ged_draw_scene_ref_erase_shapes_by_component_name(
	struct ged *gedp,
	const char *name,
	bsg_scene_ref prune_base_ref,
	void *view_ctx,
	int mode,
	int nonroot_only);
static int ged_draw_scene_ref_erase_shapes_by_db_subpath(
	struct ged *gedp,
	const struct db_full_path *subpath,
	bsg_scene_ref prune_base_ref,
	void *view_ctx,
	int mode);
static int ged_draw_scene_ref_erase_shapes_by_path_prefix_string(
	struct ged *gedp,
	const char *path,
	bsg_scene_ref prune_base_ref,
	void *view_ctx,
	int mode);
static void ged_draw_scene_ref_free_group(bsg_scene_ref ref);
static void ged_draw_scene_ref_free_group_contents(bsg_scene_ref ref);
static int ged_draw_scene_ref_prune_empty_groups(struct ged *gedp,
							 bsg_scene_ref parent_ref,
							 void *view_ctx,
							 int mode);
static void ged_draw_obol_owner_structural_revision_bump(struct ged *gedp);
static int ged_draw_scene_ref_database_source_record_get(
	bsg_scene_ref ref,
	struct ged_draw_database_source_record *record);
static int ged_draw_scene_ref_database_source_record_apply(
	bsg_scene_ref ref,
	const struct ged_draw_database_source_record *record);
static void ged_draw_scene_ref_mark_database_source_redraw_result(
	bsg_scene_ref ref,
	int failed);
static int ged_draw_scene_ref_database_source_summary(
	bsg_scene_ref ref,
	struct ged_draw_database_source_summary *out);
static int ged_draw_scene_ref_draw_state_summary(
	bsg_scene_ref ref,
	struct ged_draw_scene_draw_state_summary *out);
static int ged_draw_scene_ref_tree_summary(
	bsg_scene_ref ref,
	struct ged_draw_scene_tree_summary *out);
static int ged_draw_scene_ref_shape_record_summary(
	bsg_scene_ref ref,
	struct ged_draw_shape_record_summary *out);
static int ged_draw_scene_ref_group_record_summary(
	bsg_scene_ref ref,
	struct ged_draw_group_record_summary *out);

struct _ged_draw_obol_update_display_ctx {
    int visible_valid;
    int visible;
    int highlighted_valid;
    int highlighted;
    int draw_mode_valid;
    int draw_mode;
    int line_style_valid;
    int line_style;
    int line_width_valid;
    int line_width;
    int transparency_valid;
    fastf_t transparency;
    int color_valid;
    const unsigned char *color;
    int material_color_valid;
    const unsigned char *material_color;
    int material_revision_valid;
    uint64_t material_revision;
};
static int _ged_draw_obol_update_display_cb(struct ged *gedp,
					    const char *path,
					    void *userdata);
struct _ged_draw_obol_set_evaluated_region_ctx {
    int evaluated_region;
};
static int _ged_draw_obol_set_evaluated_region_cb(struct ged *gedp,
						  const char *path,
						  void *userdata);
struct _ged_draw_obol_set_region_metadata_ctx {
    int region_id;
    int aircode;
    int los;
    int material_id;
};
static int _ged_draw_obol_set_region_metadata_cb(struct ged *gedp,
						 const char *path,
						 void *userdata);
typedef int (*ged_draw_scene_rt_ref_child_cb)(rt_view_scene_ref child_ref,
					      void *userdata);
static int ged_draw_scene_ref_erase_component_name_at_base(
	struct ged *gedp,
	const char *name,
	bsg_scene_ref base_ref,
	void *view_ctx,
	int mode,
	int nonroot_only);
static int ged_draw_scene_rt_ref_foreach_child(
	rt_view_scene_ref scene_ref,
	ged_draw_scene_rt_ref_child_cb cb,
	void *userdata);
static int ged_draw_scene_rt_ref_display_summary(
	rt_view_scene_ref scene_ref,
	struct ged_draw_scene_display_summary *out);
static int ged_draw_scene_ref_display_summary(
	bsg_scene_ref ref,
	struct ged_draw_scene_display_summary *out);
static int ged_draw_scene_rt_ref_source_summary(
	rt_view_scene_ref scene_ref,
	struct ged_draw_database_source_summary *out);
static int ged_draw_scene_rt_ref_tree_summary(
	rt_view_scene_ref scene_ref,
	struct ged_draw_scene_tree_summary *out);
static ged_draw_shape_state *ged_draw_scene_rt_ref_registry_state(
	rt_view_scene_ref scene_ref);
static struct ged *ged_draw_scene_rt_ref_owner_gedp(
	rt_view_scene_ref scene_ref);
static const char *ged_draw_scene_rt_ref_semantic_path(
	rt_view_scene_ref scene_ref);
static const struct db_full_path *ged_draw_scene_rt_ref_fullpath(
	rt_view_scene_ref scene_ref);
static const char *ged_draw_scene_rt_ref_name(rt_view_scene_ref scene_ref);
static int ged_draw_scene_rt_ref_subtree_bounds(
	rt_view_scene_ref scene_ref,
	vect_t *min,
	vect_t *max,
	int include_overlays);
static int ged_draw_scene_ref_geometry_clear(bsg_scene_ref ref);
static int ged_draw_scene_ref_publish_bot_wireframe_line_set(
	bsg_scene_ref ref,
	const struct rt_bot_internal *bot);
static int ged_draw_scene_ref_publish_brep_wireframe_line_set(
	bsg_scene_ref ref,
	const struct rt_brep_internal *bi,
	const struct bn_tol *tol,
	int *obol_published);
static int ged_draw_scene_ref_publish_bspline_wireframe_line_set(
	bsg_scene_ref ref,
	struct rt_db_internal *ip,
	const struct bn_tol *tol,
	int *obol_published);
static int ged_draw_scene_ref_publish_line_set(
	bsg_scene_ref ref,
	const point_t *points,
	const int *commands,
	size_t point_count);
static int ged_draw_scene_ref_publish_poly_wireframe_line_set(
	bsg_scene_ref ref,
	const struct rt_pg_internal *pg);
static int ged_draw_scene_ref_publish_primitive_wireframe(
	bsg_scene_ref ref,
	struct rt_db_internal *ip,
	const struct bg_tess_tol *ttol,
	const struct bn_tol *tol,
	void *view_ctx,
	const struct rt_view_info *view_info,
	int adaptive);
static int ged_draw_scene_ref_publish_obol_submodel_wireframe(
	bsg_scene_ref ref,
	struct rt_db_internal *ip,
	const struct bg_tess_tol *ttol,
	const struct bn_tol *tol);
static int ged_draw_scene_ref_publish_obol_primitive_wireframe(
	bsg_scene_ref ref,
	struct rt_db_internal *ip,
	const struct bg_tess_tol *ttol,
	const struct bn_tol *tol);
static int ged_draw_scene_ref_update_bounds_from_geometry(
	bsg_scene_ref ref,
	int *bad_cmd);
static int ged_draw_scene_ref_update_bounds_context(
	bsg_scene_ref ref,
	void *view_ctx);
static int ged_draw_rt_internal_payload_valid(const struct rt_db_internal *ip);
static int ged_draw_scene_ref_geometry_publish_nmg_model(bsg_scene_ref ref,
							 const struct model *m,
							 int style);
static bsg_scene_ref ged_draw_scene_ref_geometry_query_target(
	bsg_scene_ref ref);
static int _ged_draw_geometry_line_command_from_bsg(int command);
static bsg_scene_ref ged_draw_scene_ref_geometry_publication_target(
	bsg_scene_ref ref);
static int ged_draw_scene_ref_bounds(bsg_scene_ref ref, point_t min, point_t max);
static int ged_draw_scene_ref_subtree_bounds(bsg_scene_ref ref,
					     vect_t *min,
					     vect_t *max,
					     int include_overlays);
static int ged_draw_scene_ref_set_changed(bsg_scene_ref ref, int changed);
static int ged_draw_scene_ref_copy_draw_intent(bsg_scene_ref dst,
					       bsg_scene_ref src);
static int ged_draw_scene_ref_set_non_database_source(
	bsg_scene_ref ref,
	int non_database_source);
static int ged_draw_obol_database_source_bot_mesh_lod_realize(
	struct ged *gedp,
	const char *path,
	struct db_i *dbip,
	struct directory *dp,
	void *view_ctx);
static int ged_draw_obol_database_source_brep_mesh_lod_realize(
	struct ged *gedp,
	const char *path,
	struct db_i *dbip,
	struct directory *dp,
	const struct bg_tess_tol *ttol,
	const struct bn_tol *tol,
	void *view_ctx);
static int ged_draw_obol_database_source_adaptive_wireframe_realize(
	struct ged *gedp,
	const char *path,
	struct db_i *dbip,
	const struct db_full_path *fullpath,
	const mat_t draw_mat,
	const struct bn_tol *tol,
	void *view_ctx,
	fastf_t draw_size);
static void ged_draw_primitive_realization_line_set_free(
	struct rt_primitive_lod_realization *realization);
static int ged_draw_obol_database_source_publish_realization_line_set(
	struct ged *gedp,
	const char *path,
	struct rt_primitive_lod_realization *realization);
static int ged_draw_scene_ref_set_material_revision(bsg_scene_ref ref,
	uint64_t mater_rev);
static void *ged_draw_scene_ref_view_context(bsg_scene_ref ref);
static int ged_draw_scene_ref_draw_mat(bsg_scene_ref ref, mat_t mat);
static fastf_t ged_draw_scene_ref_draw_size(bsg_scene_ref ref);
static void ged_draw_scene_ref_draw_center(bsg_scene_ref ref, point_t center);
static int ged_draw_scene_ref_draw_intent_is_overlay(bsg_scene_ref ref);
static void _ged_draw_appearance_to_neutral(
	struct ged_draw_appearance_settings *out,
	const struct bsg_appearance_settings *in);
static const char *ged_draw_scene_ref_draw_intent_path(bsg_scene_ref ref);
static int ged_draw_scene_ref_draw_intent_mode(bsg_scene_ref ref,
					       int *draw_mode);
static int ged_draw_scene_ref_line_style(bsg_scene_ref ref);
static int ged_draw_scene_ref_realization_pipeline_candidate(bsg_scene_ref ref);
static struct ged *_ged_draw_scene_ref_owner_gedp(bsg_scene_ref ref);
static void ged_draw_scene_ref_remember_owner_gedp(bsg_scene_ref ref,
						   struct ged *gedp);
static int ged_draw_scene_ref_obol_set_realization(
	bsg_scene_ref ref,
	int current,
	int failed,
	uint64_t realized_source_revision,
	uint64_t realized_inputs_revision,
	ged_draw_stale_reason reason);
static int ged_draw_scene_ref_obol_realization_policy_summary(
	bsg_scene_ref ref,
	struct ged_draw_obol_realization_policy_summary *out);
static void _ged_draw_scene_ref_copy_display_state(bsg_scene_ref dst, bsg_scene_ref src);
static int _ged_draw_scene_ref_apply_db_object_marker(bsg_scene_ref ref);
static bsg_scene_ref _ged_draw_scene_ref_create_child_shape(bsg_scene_ref primary_ref,
							    const char *source_name,
							    const char *shape_name,
							    bsg_scene_ref *child_source_out);
static void _ged_draw_scene_ref_clear_child_sources(bsg_scene_ref primary_ref);
static int ged_draw_scene_ref_publish_submodel_wireframe_children(bsg_scene_ref ref,
								  struct rt_db_internal *ip,
								  const struct bg_tess_tol *ttol,
								  const struct bn_tol *tol,
								  void *view_ctx);
enum {
    GED_DRAW_NMG_CNURB_INTERIOR_SAMPLES = 10,
    GED_DRAW_NMG_CNURB_SAMPLE_POINTS = GED_DRAW_NMG_CNURB_INTERIOR_SAMPLES + 2,
    GED_DRAW_NMG_SNURB_INTERIOR_SAMPLES = 10
};


struct ged_draw_source_state {
    struct db_i *dbip;
    struct db_full_path *fp;
    const struct bn_tol *tol;
    const struct bg_tess_tol *ttol;
    struct BRLObolMeshLod *mesh_lod;
    point_t mesh_lod_bmin;
    point_t mesh_lod_bmax;
    int mesh_lod_bounds_valid;
    uint64_t source_revision;
    uint64_t inputs_revision;
    uint64_t realized_source_revision;
    uint64_t realized_inputs_revision;
    ged_draw_stale_reason stale_reason;
};


struct ged_draw_source_state_record {
    struct db_i *dbip;
    const struct db_full_path *fullpath;
    const struct bn_tol *tol;
    const struct bg_tess_tol *ttol;
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


struct ged_draw_source_runtime_summary {
    struct db_i *dbip;
    const struct db_full_path *fullpath;
    struct directory *leaf_dp;
    const char *name;
    const struct bn_tol *tol;
    const struct bg_tess_tol *ttol;
    struct BRLObolMeshLod *mesh_lod;
    int mesh_lod_bounds_valid;
    point_t mesh_lod_bmin;
    point_t mesh_lod_bmax;
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


static void
_ged_draw_view_export_detail_init(struct ged_draw_view_export_detail *detail)
{
    if (!detail)
	return;
    memset(detail, 0, sizeof(*detail));
    BU_VLS_INIT(&detail->path);
    MAT_IDN(detail->record.model_mat);
}


static void
_ged_draw_view_export_detail_free(struct ged_draw_view_export_detail *detail)
{
    if (!detail)
	return;
    bu_vls_free(&detail->path);
    if (detail->arrays.points)
	bu_free(detail->arrays.points, "ged draw view export array points");
    if (detail->arrays.commands)
	bu_free(detail->arrays.commands, "ged draw view export array commands");
    if (detail->arrays.indices)
	bu_free(detail->arrays.indices, "ged draw view export array indices");
    if (detail->surface.points)
	bu_free(detail->surface.points, "ged draw view export surface points");
    if (detail->surface.indices)
	bu_free(detail->surface.indices, "ged draw view export surface indices");
    if (detail->annotation.points)
	bu_free(detail->annotation.points,
		"ged draw view export annotation points");
    if (detail->annotation.text)
	bu_free(detail->annotation.text, "ged draw view export annotation text");
}


static int
_ged_draw_view_export_points_copy(point_t **dst,
				  const point_t *src,
				  size_t count,
				  const char *label)
{
    if (!dst)
	return 0;
    *dst = NULL;
    if (!count)
	return 1;
    if (!src)
	return 1;

    *dst = (point_t *)bu_calloc(count, sizeof(point_t), label);
    for (size_t i = 0; i < count; i++)
	VMOVE((*dst)[i], src[i]);
    return 1;
}


static int
_ged_draw_view_export_ints_copy(int **dst,
				const int *src,
				size_t count,
				const char *label)
{
    if (!dst)
	return 0;
    *dst = NULL;
    if (!count)
	return 1;
    if (!src)
	return 1;

    *dst = (int *)bu_calloc(count, sizeof(int), label);
    memcpy(*dst, src, count * sizeof(int));
    return 1;
}


static enum ged_draw_view_export_geometry_kind
_ged_draw_view_export_geometry_kind_from_bsg(bsg_render_geometry_kind kind)
{
    switch (kind) {
	case BSG_RENDER_GEOMETRY_MESH:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_MESH;
	case BSG_RENDER_GEOMETRY_LINE_SET:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET;
	case BSG_RENDER_GEOMETRY_INDEXED_LINE_SET:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_LINE_SET;
	case BSG_RENDER_GEOMETRY_POINT_SET:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_POINT_SET;
	case BSG_RENDER_GEOMETRY_INDEXED_FACE_SET:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_FACE_SET;
	case BSG_RENDER_GEOMETRY_TEXT:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_TEXT;
	case BSG_RENDER_GEOMETRY_IMAGE:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_IMAGE;
	case BSG_RENDER_GEOMETRY_OVERLAY:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_OVERLAY;
	case BSG_RENDER_GEOMETRY_HUD:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_HUD;
	case BSG_RENDER_GEOMETRY_ANNOTATION:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_ANNOTATION;
	case BSG_RENDER_GEOMETRY_CSG_PROXY:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_CSG_PROXY;
	case BSG_RENDER_GEOMETRY_BREP_PROXY:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_BREP_PROXY;
	case BSG_RENDER_GEOMETRY_EDIT_PREVIEW:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_EDIT_PREVIEW;
	case BSG_RENDER_GEOMETRY_EXTERNAL_PROXY:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_EXTERNAL_PROXY;
	case BSG_RENDER_GEOMETRY_NONE:
	default:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_NONE;
    }
}


static const char *
_ged_draw_view_export_type_name_from_bsg(
	bsg_render_geometry_role geometry_role)
{
    switch (geometry_role) {
	case BSG_RENDER_GEOMETRY_ROLE_AXES_WIDGET:
	    return "axes";
	case BSG_RENDER_GEOMETRY_ROLE_LINE_SET:
	    return "line";
	case BSG_RENDER_GEOMETRY_ROLE_TEXT_LABEL:
	    return "label";
	case BSG_RENDER_GEOMETRY_ROLE_POLYGON_REGION:
	    return "polygon";
	case BSG_RENDER_GEOMETRY_ROLE_DATABASE_OBJECT:
	    return "gobj";
	default:
	    return "object";
    }
}


static const char *
_ged_draw_view_export_geometry_name_from_bsg(
	bsg_render_geometry_kind geometry_kind)
{
    switch (geometry_kind) {
	case BSG_RENDER_GEOMETRY_MESH:
	    return "mesh";
	case BSG_RENDER_GEOMETRY_LINE_SET:
	    return "line-set";
	case BSG_RENDER_GEOMETRY_INDEXED_LINE_SET:
	    return "indexed-line-set";
	case BSG_RENDER_GEOMETRY_POINT_SET:
	    return "point-set";
	case BSG_RENDER_GEOMETRY_INDEXED_FACE_SET:
	    return "indexed-face-set";
	case BSG_RENDER_GEOMETRY_TEXT:
	    return "text";
	case BSG_RENDER_GEOMETRY_IMAGE:
	    return "image";
	case BSG_RENDER_GEOMETRY_OVERLAY:
	    return "overlay";
	case BSG_RENDER_GEOMETRY_HUD:
	    return "hud";
	case BSG_RENDER_GEOMETRY_ANNOTATION:
	    return "annotation";
	case BSG_RENDER_GEOMETRY_CSG_PROXY:
	    return "csg-proxy";
	case BSG_RENDER_GEOMETRY_BREP_PROXY:
	    return "brep-proxy";
	case BSG_RENDER_GEOMETRY_EDIT_PREVIEW:
	    return "edit-preview";
	case BSG_RENDER_GEOMETRY_EXTERNAL_PROXY:
	    return "external-proxy";
	case BSG_RENDER_GEOMETRY_NONE:
	default:
	    return "none";
    }
}


static int
_ged_draw_view_line_command_from_bsg(int command)
{
    switch (command) {
	case BSG_GEOMETRY_LINE_MOVE:
	    return GED_DRAW_VIEW_LINE_MOVE;
	case BSG_GEOMETRY_LINE_DRAW:
	    return GED_DRAW_VIEW_LINE_DRAW;
	case BSG_GEOMETRY_POINT_DRAW:
	    return GED_DRAW_VIEW_LINE_POINT_DRAW;
	default:
	    return command;
    }
}


static int
_ged_draw_view_export_detail_arrays_from_bsg(
	struct ged_draw_view_export_detail *detail,
	const struct bsg_render_geometry *geometry)
{
    if (!detail || !geometry)
	return 0;

    detail->arrays.point_count = geometry->arrays.point_count;
    detail->arrays.command_count = geometry->arrays.command_count;
    detail->arrays.index_count = geometry->arrays.index_count;
    if (!_ged_draw_view_export_points_copy(&detail->arrays.points,
	    (const point_t *)geometry->arrays.points,
	    geometry->arrays.point_count,
	    "ged draw view export array points"))
	return 0;
    if (!_ged_draw_view_export_ints_copy(&detail->arrays.commands,
	    geometry->arrays.commands, geometry->arrays.command_count,
	    "ged draw view export array commands"))
	return 0;
    for (size_t i = 0; i < detail->arrays.command_count; i++)
	detail->arrays.commands[i] =
	    _ged_draw_view_line_command_from_bsg(detail->arrays.commands[i]);
    if (!_ged_draw_view_export_ints_copy(&detail->arrays.indices,
	    geometry->arrays.indices, geometry->arrays.index_count,
	    "ged draw view export array indices"))
	return 0;

    return 1;
}


static int
_ged_draw_view_export_detail_surface_from_bsg(
	struct ged_draw_view_export_detail *detail,
	const struct bsg_render_geometry *geometry)
{
    if (!detail || !geometry)
	return 0;

    detail->surface.point_count = geometry->surface.point_count;
    detail->surface.normal_count = geometry->surface.normal_count;
    detail->surface.index_count = geometry->surface.index_count;
    detail->surface.face_count = geometry->surface.face_count;
    detail->surface.normals_per_index =
	(geometry->surface.normal_kind == BSG_RENDER_SURFACE_NORMALS_PER_INDEX);
    detail->surface.material_valid = geometry->surface.material_valid;
    detail->surface.material_draw_mode =
	geometry->surface.material.draw_mode;
    detail->surface.material_transparency =
	geometry->surface.material.transparency;
    detail->surface.material_highlighted =
	geometry->surface.material.highlighted;
    detail->surface.material_color[0] =
	geometry->surface.material.color[0];
    detail->surface.material_color[1] =
	geometry->surface.material.color[1];
    detail->surface.material_color[2] =
	geometry->surface.material.color[2];

    if (!_ged_draw_view_export_points_copy(&detail->surface.points,
	    (const point_t *)geometry->surface.points,
	    geometry->surface.point_count,
	    "ged draw view export surface points"))
	return 0;
    if (!_ged_draw_view_export_ints_copy(&detail->surface.indices,
	    geometry->surface.indices, geometry->surface.index_count,
	    "ged draw view export surface indices"))
	return 0;

    return 1;
}


static int
_ged_draw_view_export_detail_annotation_from_bsg(
	struct ged_draw_view_export_detail *detail,
	const struct bsg_render_geometry *geometry)
{
    if (!detail || !geometry)
	return 0;

    detail->annotation.display_space =
	(geometry->annotation.space == BSG_ANNOTATION_SPACE_DISPLAY);
    detail->annotation.point_count = geometry->annotation.point_count;
    detail->annotation.segment_count = geometry->annotation.segment_count;
    if (!_ged_draw_view_export_points_copy(&detail->annotation.points,
	    (const point_t *)geometry->annotation.points,
	    geometry->annotation.point_count,
	    "ged draw view export annotation points"))
	return 0;

    const struct bsg_annotation_segment *segments =
	geometry->annotation.segments;
    if (!segments)
	return 1;

    for (size_t i = 0; i < geometry->annotation.segment_count; i++) {
	const struct bsg_annotation_segment *seg = &segments[i];
	if (!detail->annotation.line_segment_valid &&
		seg->kind == BSG_ANNOTATION_SEGMENT_LINE) {
	    detail->annotation.line_segment_valid = 1;
	    detail->annotation.line_start = seg->data.line.start;
	    detail->annotation.line_end = seg->data.line.end;
	}
	if (!detail->annotation.text_segment_valid &&
		seg->kind == BSG_ANNOTATION_SEGMENT_TEXT) {
	    detail->annotation.text_segment_valid = 1;
	    detail->annotation.text_ref_point = seg->data.text.ref_pt;
	    if (seg->data.text.text)
		detail->annotation.text = bu_strdup(seg->data.text.text);
	}
    }

    return 1;
}


static int
_ged_draw_view_export_detail_geometry_from_bsg(
	struct ged_draw_view_export_detail *detail,
	const struct bsg_export_record *export_rec)
{
    if (!detail || !export_rec)
	return 0;

    detail->geometry_kind =
	_ged_draw_view_export_geometry_kind_from_bsg(
	    export_rec->geometry.kind);
    switch (detail->geometry_kind) {
	case GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET:
	case GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_LINE_SET:
	case GED_DRAW_VIEW_EXPORT_GEOMETRY_POINT_SET:
	    return _ged_draw_view_export_detail_arrays_from_bsg(detail,
		    &export_rec->geometry);
	case GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_FACE_SET:
	    return _ged_draw_view_export_detail_surface_from_bsg(detail,
		    &export_rec->geometry);
	case GED_DRAW_VIEW_EXPORT_GEOMETRY_ANNOTATION:
	    return _ged_draw_view_export_detail_annotation_from_bsg(detail,
		    &export_rec->geometry);
	default:
	    return 1;
    }
}


static int
_ged_draw_view_export_detail_from_bsg(
	struct ged_draw_view_export_detail *detail,
	const struct bsg_export_record *export_rec)
{
    if (!detail || !export_rec)
	return 0;

    _ged_draw_view_export_detail_init(detail);
    bu_vls_strcpy(&detail->path, bu_vls_cstr(&export_rec->path));

    struct ged_draw_view_db_object_record *rec = &detail->record;
    rec->path = bu_vls_cstr(&detail->path);
    rec->type_name =
	_ged_draw_view_export_type_name_from_bsg(
	    export_rec->source.geometry_role);
    rec->geometry_name =
	_ged_draw_view_export_geometry_name_from_bsg(
	    export_rec->geometry.kind);
    rec->draw_mode = export_rec->draw_mode;
    rec->transparency = export_rec->transparency;
    rec->evaluated_region = export_rec->evaluated_region;
    rec->is_database_source =
	(export_rec->source.scope == BSG_RENDER_SOURCE_SCOPE_DATABASE);
    rec->non_database_source = export_rec->non_database_source;
    rec->is_database_intent =
	(export_rec->source.draw_intent == BSG_RENDER_DRAW_INTENT_DATABASE);
    rec->is_local_source =
	(export_rec->source.scope == BSG_RENDER_SOURCE_SCOPE_VIEW_LOCAL);
    rec->is_view_source =
	(export_rec->source.scope == BSG_RENDER_SOURCE_SCOPE_VIEW_SHARED ||
	 export_rec->source.scope == BSG_RENDER_SOURCE_SCOPE_VIEW_LOCAL ||
	 export_rec->source.draw_intent == BSG_RENDER_DRAW_INTENT_OVERLAY ||
	 export_rec->source.draw_intent == BSG_RENDER_DRAW_INTENT_HUD);
    rec->highlighted = export_rec->highlighted;
    rec->selected = export_rec->selected;
    rec->visible = export_rec->visible;
    rec->line_style = export_rec->line_style;
    rec->color[0] = export_rec->color[0];
    rec->color[1] = export_rec->color[1];
    rec->color[2] = export_rec->color[2];
    MAT_COPY(rec->model_mat, export_rec->model_mat);
    VMOVE(rec->bounds_center, export_rec->bounds_center);
    rec->bounds_radius = export_rec->bounds_radius;
    rec->has_bounds = export_rec->has_bounds;
    rec->vlist_structure_count = export_rec->vlist_structure_count;
    rec->vlist_point_count = export_rec->vlist_point_count;
    rec->cache_identity = export_rec->cache_identity;
    rec->source_identity = export_rec->source.source_id;
    rec->detail_token = (uintptr_t)detail;

    if (!_ged_draw_view_export_detail_geometry_from_bsg(detail, export_rec)) {
	_ged_draw_view_export_detail_free(detail);
	return 0;
    }

    return 1;
}

static uint64_t
_ged_draw_view_export_hash_cstr(const char *str)
{
    if (!str || !str[0])
	return 0;
    return bu_data_hash(str, strlen(str) * sizeof(char));
}


static int
_ged_draw_view_export_glob_match(const char *glob,
				 const char *record_path,
				 const char *shape_path)
{
    if (!glob || !glob[0])
	return 1;

    if (record_path && bu_path_match(glob, record_path, 0) == 0)
	return 1;
    if (record_path && record_path[0] == '/' &&
	    bu_path_match(glob, record_path + 1, 0) == 0)
	return 1;
    if (shape_path && bu_path_match(glob, shape_path, 0) == 0)
	return 1;
    if (shape_path && shape_path[0] == '/' &&
	    bu_path_match(glob, shape_path + 1, 0) == 0)
	return 1;

    return 0;
}


static int
_ged_draw_view_export_obol_display_matches(
	const struct ged_draw_scene_display_summary *display,
	unsigned int query_flags,
	unsigned int render_flags,
	const char *glob,
	int draw_mode,
	const char *record_path,
	const char *shape_path)
{
    if (!display || !display->valid)
	return 0;

    if (((query_flags & GED_DRAW_VIEW_EXPORT_QUERY_VISIBLE_ONLY) ||
		(render_flags & GED_DRAW_VIEW_EXPORT_RENDER_VISIBLE_ONLY)) &&
	    !display->visible)
	return 0;

    if (draw_mode != GED_DRAW_VIEW_EXPORT_DRAW_MODE_ANY &&
	    display->draw_mode != draw_mode)
	return 0;

    const int is_database_source = display->is_database_source ? 1 : 0;
    const int is_view_source = is_database_source ? 0 : 1;
    const int wants_database =
	(query_flags & GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS) ? 1 : 0;
    const int wants_view =
	(query_flags & GED_DRAW_VIEW_EXPORT_QUERY_VIEW_OBJECTS) ? 1 : 0;

    if ((query_flags & GED_DRAW_VIEW_EXPORT_QUERY_LOCAL_ONLY) &&
	    !is_view_source)
	return 0;

    if ((wants_database || wants_view) &&
	    ((is_database_source && !wants_database) ||
	     (is_view_source && !wants_view)))
	return 0;

    return _ged_draw_view_export_glob_match(glob, record_path, shape_path);
}


static const char *
_ged_draw_view_export_type_name_from_obol_geometry(const char *geometry_name)
{
    if (geometry_name && BU_STR_EQUAL(geometry_name, "line-set"))
	return "line";
    if (geometry_name && BU_STR_EQUAL(geometry_name, "indexed-face-set"))
	return "polygon";
    if (geometry_name && BU_STR_EQUAL(geometry_name, "annotation"))
	return "annotation";
    return "object";
}


static void
_ged_draw_view_export_record_bounds_from_points(
	struct ged_draw_view_db_object_record *rec,
	const point_t *points,
	size_t point_count)
{
    if (!rec || !points || !point_count)
	return;

    point_t bmin;
    point_t bmax;
    VMOVE(bmin, points[0]);
    VMOVE(bmax, points[0]);
    for (size_t i = 1; i < point_count; i++)
	VMINMAX(bmin, bmax, points[i]);

    rec->bounds_center[X] = 0.5 * (bmin[X] + bmax[X]);
    rec->bounds_center[Y] = 0.5 * (bmin[Y] + bmax[Y]);
    rec->bounds_center[Z] = 0.5 * (bmin[Z] + bmax[Z]);
    rec->bounds_radius = 0.0;
    for (size_t i = 0; i < point_count; i++) {
	vect_t delta;
	VSUB2(delta, points[i], rec->bounds_center);
	const fastf_t radius = MAGNITUDE(delta);
	if (radius > rec->bounds_radius)
	    rec->bounds_radius = radius;
    }
    rec->has_bounds = 1;
}


static int
_ged_draw_view_export_detail_line_from_obol(
	struct ged_draw_view_export_detail *detail,
	struct ged *gedp,
	const char *shape_path,
	int database_source)
{
    if (!detail || !gedp || !shape_path)
	return 0;

    struct ged_draw_view_line_summary line_summary;
    int valid = database_source ?
	ged_draw_obol_database_source_line_summary_for_path(gedp,
		shape_path, &line_summary) :
	ged_draw_obol_shape_line_summary_for_path(gedp, shape_path,
		&line_summary);
    if (!valid || !line_summary.valid)
	return 0;

    detail->geometry_kind = GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET;
    detail->arrays.point_count = line_summary.point_count;
    detail->arrays.command_count = line_summary.point_count;
    if (line_summary.point_count) {
	detail->arrays.points = (point_t *)bu_calloc(line_summary.point_count,
		sizeof(point_t), "GED Obol export line points");
	detail->arrays.commands = (int *)bu_calloc(line_summary.point_count,
		sizeof(int), "GED Obol export line commands");
    }

    size_t structure_count = 0;
    for (size_t i = 0; i < line_summary.point_count; i++) {
	int point_valid = database_source ?
	    ged_draw_obol_database_source_line_point_at_for_path(gedp,
		    shape_path, i, detail->arrays.points[i]) :
	    ged_draw_obol_shape_line_point_at_for_path(gedp, shape_path, i,
		    detail->arrays.points[i]);
	if (!point_valid)
	    return 0;
	int command_valid = database_source ?
	    ged_draw_obol_database_source_line_command_at_for_path(gedp,
		    shape_path, i, &detail->arrays.commands[i]) :
	    ged_draw_obol_shape_line_command_at_for_path(gedp, shape_path, i,
		    &detail->arrays.commands[i]);
	if (!command_valid)
	    detail->arrays.commands[i] =
		(i % 2) ? GED_DRAW_VIEW_LINE_DRAW : GED_DRAW_VIEW_LINE_MOVE;
	if (detail->arrays.commands[i] != GED_DRAW_VIEW_LINE_MOVE)
	    structure_count++;
    }

    detail->record.vlist_structure_count = structure_count;
    detail->record.vlist_point_count = line_summary.point_count;
    _ged_draw_view_export_record_bounds_from_points(&detail->record,
	    detail->arrays.points, detail->arrays.point_count);
    return 1;
}


static int
_ged_draw_view_export_detail_surface_from_obol(
	struct ged_draw_view_export_detail *detail,
	struct ged *gedp,
	const char *shape_path,
	int database_source)
{
    if (!detail || !gedp || !shape_path)
	return 0;

    struct ged_draw_view_surface_summary surface_summary;
    int valid = database_source ?
	ged_draw_obol_database_source_surface_summary_for_path(gedp,
		shape_path, &surface_summary) :
	ged_draw_obol_shape_surface_summary_for_path(gedp, shape_path,
		&surface_summary);
    if (!valid || !surface_summary.valid)
	return 0;

    detail->geometry_kind = GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_FACE_SET;
    detail->surface.point_count = surface_summary.point_count;
    detail->surface.normal_count = surface_summary.normal_count;
    detail->surface.index_count = surface_summary.index_count;
    detail->surface.face_count = surface_summary.face_count;
    detail->surface.normals_per_index = surface_summary.normals_per_index;
    detail->surface.material_valid = surface_summary.material_valid;
    detail->surface.material_draw_mode = surface_summary.material_draw_mode;
    detail->surface.material_transparency =
	surface_summary.material_transparency;
    detail->surface.material_highlighted =
	surface_summary.material_highlighted;
    detail->surface.material_color[0] = surface_summary.material_color[0];
    detail->surface.material_color[1] = surface_summary.material_color[1];
    detail->surface.material_color[2] = surface_summary.material_color[2];

    if (surface_summary.point_count) {
	detail->surface.points = (point_t *)bu_calloc(
		surface_summary.point_count, sizeof(point_t),
		"GED Obol export surface points");
    }
    if (surface_summary.index_count) {
	detail->surface.indices = (int *)bu_calloc(
		surface_summary.index_count, sizeof(int),
		"GED Obol export surface indices");
    }

    for (size_t i = 0; i < surface_summary.point_count; i++) {
	int point_valid = database_source ?
	    ged_draw_obol_database_source_surface_point_at_for_path(gedp,
		    shape_path, i, detail->surface.points[i]) :
	    ged_draw_obol_shape_surface_point_at_for_path(gedp, shape_path,
		    i, detail->surface.points[i]);
	if (!point_valid)
	    return 0;
    }
    for (size_t i = 0; i < surface_summary.index_count; i++) {
	int index_valid = database_source ?
	    ged_draw_obol_database_source_surface_index_at_for_path(gedp,
		    shape_path, i, &detail->surface.indices[i]) :
	    ged_draw_obol_shape_surface_index_at_for_path(gedp, shape_path,
		    i, &detail->surface.indices[i]);
	if (!index_valid)
	    return 0;
    }

    if (surface_summary.cache_identity)
	detail->record.cache_identity = surface_summary.cache_identity;
    if (surface_summary.source_identity)
	detail->record.source_identity = surface_summary.source_identity;
    _ged_draw_view_export_record_bounds_from_points(&detail->record,
	    detail->surface.points, detail->surface.point_count);
    return 1;
}


static int
_ged_draw_view_export_detail_annotation_from_obol(
	struct ged_draw_view_export_detail *detail,
	struct ged *gedp,
	const char *shape_path)
{
    if (!detail || !gedp || !shape_path)
	return 0;

    struct ged_draw_obol_annotation_summary summary;
    memset(&summary, 0, sizeof(summary));
    if (!ged_draw_obol_database_source_annotation_summary_for_path(gedp,
	    shape_path, &summary) || !summary.valid)
	return 0;

    detail->geometry_kind = GED_DRAW_VIEW_EXPORT_GEOMETRY_ANNOTATION;
    detail->annotation.display_space = 1;
    detail->annotation.point_count = summary.point_count;
    detail->annotation.segment_count = summary.segment_count;
    detail->annotation.line_segment_valid = summary.line_segment_valid;
    detail->annotation.line_start = summary.line_start;
    detail->annotation.line_end = summary.line_end;
    detail->annotation.text_segment_valid = summary.text_segment_valid;
    detail->annotation.text_ref_point = summary.text_ref_point;
    if (summary.text)
	detail->annotation.text = bu_strdup(summary.text);
    if (summary.cache_identity)
	detail->record.cache_identity = summary.cache_identity;
    if (summary.source_identity)
	detail->record.source_identity = summary.source_identity;

    if (summary.point_count) {
	detail->annotation.points = (point_t *)bu_calloc(
		summary.point_count, sizeof(point_t),
		"GED Obol export annotation points");
    }

    for (size_t i = 0; i < summary.point_count; i++) {
	if (!ged_draw_obol_database_source_annotation_point_at_for_path(
		gedp, shape_path, i, detail->annotation.points[i]))
	    return 0;
    }

    _ged_draw_view_export_record_bounds_from_points(&detail->record,
	    detail->annotation.points, detail->annotation.point_count);
    return 1;
}


static int
_ged_draw_view_export_detail_geometry_from_obol(
	struct ged_draw_view_export_detail *detail,
	struct ged *gedp,
	const char *shape_path,
	const char *geometry_name,
	int database_source)
{
    if (!detail || !geometry_name)
	return 0;

    if (BU_STR_EQUAL(geometry_name, "line-set"))
	return _ged_draw_view_export_detail_line_from_obol(detail, gedp,
		shape_path, database_source);
    if (BU_STR_EQUAL(geometry_name, "indexed-face-set"))
	return _ged_draw_view_export_detail_surface_from_obol(detail, gedp,
		shape_path, database_source);
    if (BU_STR_EQUAL(geometry_name, "annotation"))
	return database_source ?
	    _ged_draw_view_export_detail_annotation_from_obol(detail, gedp,
		    shape_path) : 0;
    return 0;
}


static int
_ged_draw_view_export_detail_from_obol_shape(
	struct ged_draw_view_export_detail *detail,
	struct ged *gedp,
	void *view_ctx,
	const char *shape_path,
	unsigned int query_flags,
	unsigned int render_flags,
	const char *glob,
	int draw_mode)
{
    if (!detail || !gedp || !shape_path || !shape_path[0])
	return 0;

    _ged_draw_view_export_detail_init(detail);

    struct ged_draw_scene_display_summary display;
    struct ged_draw_shape_geometry_summary geometry;
    memset(&display, 0, sizeof(display));
    memset(&geometry, 0, sizeof(geometry));
    int database_source = 0;
    const char *geometry_path = shape_path;
    int valid = ged_draw_obol_shape_display_summary_for_path(gedp,
	    shape_path, &display);
    if (valid && display.valid) {
	valid = ged_draw_obol_shape_geometry_summary_for_path(gedp,
		shape_path, &geometry);
    }
    if (valid && display.valid && geometry.valid && geometry.geometry_name &&
	    display.is_database_source && display.intent_path &&
	    display.intent_path[0]) {
	struct ged_draw_scene_display_summary source_display;
	struct ged_draw_shape_geometry_summary source_geometry;
	memset(&source_display, 0, sizeof(source_display));
	memset(&source_geometry, 0, sizeof(source_geometry));
	const char *source_path = display.intent_path;
	int source_valid =
	    ged_draw_obol_database_source_display_summary_for_path(gedp,
		source_path, &source_display) &&
	    source_display.valid &&
	    ged_draw_obol_database_source_geometry_summary_for_path(gedp,
		source_path, &source_geometry) &&
	    source_geometry.valid && source_geometry.geometry_name;
	if (!source_valid && !BU_STR_EQUAL(source_path, shape_path)) {
	    memset(&source_display, 0, sizeof(source_display));
	    memset(&source_geometry, 0, sizeof(source_geometry));
	    source_path = shape_path;
	    source_valid =
		ged_draw_obol_database_source_display_summary_for_path(gedp,
		    source_path, &source_display) &&
		source_display.valid &&
		ged_draw_obol_database_source_geometry_summary_for_path(gedp,
		    source_path, &source_geometry) &&
		source_geometry.valid && source_geometry.geometry_name;
	}
	if (source_valid) {
	    display = source_display;
	    geometry = source_geometry;
	    database_source = 1;
	    geometry_path = source_path;
	}
    }
    if (!valid || !display.valid || !geometry.valid ||
	    !geometry.geometry_name) {
	memset(&display, 0, sizeof(display));
	memset(&geometry, 0, sizeof(geometry));
	if (!ged_draw_obol_database_source_display_summary_for_path(gedp,
		shape_path, &display) || !display.valid ||
		!ged_draw_obol_database_source_geometry_summary_for_path(gedp,
		    shape_path, &geometry) || !geometry.valid ||
		!geometry.geometry_name) {
	    _ged_draw_view_export_detail_free(detail);
	    return 0;
	}
	database_source = 1;
	geometry_path = shape_path;
    }

    if (!geometry.geometry_name) {
	_ged_draw_view_export_detail_free(detail);
	return 0;
    }

    const char *record_path =
	(display.intent_path && display.intent_path[0]) ?
	display.intent_path : shape_path;
    if (!_ged_draw_view_export_obol_display_matches(&display, query_flags,
	    render_flags, glob, draw_mode, record_path, shape_path)) {
	_ged_draw_view_export_detail_free(detail);
	return 0;
    }

    bu_vls_strcpy(&detail->path, record_path);

    struct ged_draw_view_db_object_record *rec = &detail->record;
    rec->path = bu_vls_cstr(&detail->path);
    rec->type_name =
	_ged_draw_view_export_type_name_from_obol_geometry(
		geometry.geometry_name);
    rec->geometry_name = geometry.geometry_name;
    rec->draw_mode = display.draw_mode;
    rec->transparency = display.transparency;
    rec->is_database_source = display.is_database_source ? 1 : 0;
    rec->non_database_source = rec->is_database_source ? 0 : 1;
    rec->is_database_intent =
	(rec->is_database_source && display.has_draw_intent) ? 1 : 0;
    rec->is_local_source = rec->is_database_source ? 0 : 1;
    rec->is_view_source = rec->is_local_source ? 1 : 0;
    rec->highlighted = display.highlighted;
    rec->selected = 0;
    if (view_ctx) {
	rec->selected = ged_draw_view_context_selection_contains_path(
		view_ctx, GED_DRAW_VIEW_SELECTION_SELECTED_PATH,
		record_path) ? 1 : 0;
	if (!rec->selected && shape_path &&
		(!record_path || !BU_STR_EQUAL(record_path, shape_path))) {
	    rec->selected = ged_draw_view_context_selection_contains_path(
		    view_ctx, GED_DRAW_VIEW_SELECTION_SELECTED_PATH,
		    shape_path) ? 1 : 0;
	}
    }
    rec->visible = display.visible;
    rec->line_style = display.line_style;
    rec->color[0] = display.material_color[0];
    rec->color[1] = display.material_color[1];
    rec->color[2] = display.material_color[2];
    MAT_IDN(rec->model_mat);
    rec->cache_identity = _ged_draw_view_export_hash_cstr(shape_path);
    rec->source_identity = _ged_draw_view_export_hash_cstr(record_path);
    rec->detail_token = (uintptr_t)detail;

    if (!_ged_draw_view_export_detail_geometry_from_obol(detail, gedp,
	    geometry_path, geometry.geometry_name, database_source)) {
	_ged_draw_view_export_detail_free(detail);
	return 0;
    }

    if (!rec->cache_identity)
	rec->cache_identity = _ged_draw_view_export_hash_cstr(shape_path);
    if (!rec->source_identity)
	rec->source_identity = _ged_draw_view_export_hash_cstr(record_path);
    if (!rec->source_identity)
	rec->source_identity = rec->cache_identity;

    return 1;
}


struct ged_draw_obol_view_export_ctx {
    struct ged *gedp;
    void *view_ctx;
    unsigned int query_flags;
    unsigned int render_flags;
    const char *glob;
    int draw_mode;
    ged_draw_view_db_object_record_cb cb;
    void *userdata;
    int keep_going;
};


static void
_ged_draw_view_export_visit_obol_info(
	struct ged_draw_obol_view_export_ctx *ctx,
	const struct ged_draw_obol_scene_context_info *info,
	const char *parent_path);


static void
_ged_draw_view_export_visit_obol_children(
	struct ged_draw_obol_view_export_ctx *ctx,
	const struct ged_draw_obol_scene_context_info *info)
{
    if (!ctx || !ctx->keep_going || !info || !info->path)
	return;

    for (size_t i = 0; ctx->keep_going && i < info->child_count; i++) {
	struct ged_draw_obol_scene_context_info child_info;
	memset(&child_info, 0, sizeof(child_info));
	if (!ged_draw_obol_scene_child_context_info_for_path(ctx->gedp,
		info->path, i, &child_info))
	    continue;
	_ged_draw_view_export_visit_obol_info(ctx, &child_info, info->path);
	ged_draw_obol_scene_context_info_free(&child_info);
    }
}


static void
_ged_draw_view_export_visit_obol_info(
	struct ged_draw_obol_view_export_ctx *ctx,
	const struct ged_draw_obol_scene_context_info *info,
	const char *parent_path)
{
    if (!ctx || !ctx->keep_going || !info || !info->path)
	return;

    if (info->is_shape || info->is_database_source) {
	struct ged_draw_view_export_detail detail;
	if (_ged_draw_view_export_detail_from_obol_shape(&detail, ctx->gedp,
		ctx->view_ctx, info->path, ctx->query_flags,
		ctx->render_flags, ctx->glob, ctx->draw_mode)) {
	    ctx->keep_going = ctx->cb(&detail.record, ctx->userdata);
	    _ged_draw_view_export_detail_free(&detail);
	}
	if (!ctx->keep_going)
	    return;
    }

    if (parent_path && BU_STR_EQUAL(info->path, parent_path))
	return;

    _ged_draw_view_export_visit_obol_children(ctx, info);
}


static const char *
_ged_draw_view_export_type_name_from_obol_feature(
	const struct ged_draw_obol_view_feature_record *feature)
{
    if (!feature)
	return "object";

    switch (feature->kind) {
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_LABELS:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_HUD_LABEL:
	    return "label";
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_AXES:
	    return "axes";
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_INDEXED_FACE_SET:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_POLYGON_OVERLAY:
	    return "polygon";
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_POINTS:
	    return "point";
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINES:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_INDEXED_LINES:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_ARROW:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINE_LAYER:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_EDIT_PREVIEW:
	    return "line";
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_CUSTOM_NODE:
	    return "object";
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_UNKNOWN:
	default:
	    return "object";
    }
}


static const char *
_ged_draw_view_export_geometry_name_from_obol_feature(
	const struct ged_draw_obol_view_feature_record *feature)
{
    if (!feature)
	return NULL;

    switch (feature->kind) {
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_INDEXED_FACE_SET:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_POLYGON_OVERLAY:
	    return "indexed-face-set";
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_POINTS:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_AXES:
	    return "point-set";
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_LABELS:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_HUD_LABEL:
	    return "annotation";
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINES:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_INDEXED_LINES:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_ARROW:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINE_LAYER:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_EDIT_PREVIEW:
	    return "line-set";
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_CUSTOM_NODE:
	    return NULL;
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_UNKNOWN:
	default:
	    return NULL;
    }
}


static enum ged_draw_view_export_geometry_kind
_ged_draw_view_export_geometry_kind_from_obol_feature(
	const struct ged_draw_obol_view_feature_record *feature)
{
    if (!feature)
	return GED_DRAW_VIEW_EXPORT_GEOMETRY_NONE;

    switch (feature->kind) {
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_INDEXED_FACE_SET:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_POLYGON_OVERLAY:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_FACE_SET;
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_POINTS:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_AXES:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_POINT_SET;
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_LABELS:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_HUD_LABEL:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_ANNOTATION;
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINES:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_INDEXED_LINES:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_ARROW:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINE_LAYER:
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_EDIT_PREVIEW:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET;
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_CUSTOM_NODE:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_NONE;
	case GED_DRAW_OBOL_VIEW_FEATURE_KIND_UNKNOWN:
	default:
	    return GED_DRAW_VIEW_EXPORT_GEOMETRY_NONE;
    }
}


static size_t
_ged_draw_view_export_face_count_from_indices(const int *indices,
					      size_t index_count)
{
    if (!indices || !index_count)
	return 0;

    size_t faces = 0;
    size_t verts = 0;
    int have_separator = 0;
    for (size_t i = 0; i < index_count; i++) {
	if (indices[i] < 0) {
	    have_separator = 1;
	    if (verts >= 3)
		faces++;
	    verts = 0;
	} else {
	    verts++;
	}
    }
    if (have_separator) {
	if (verts >= 3)
	    faces++;
	return faces;
    }

    return index_count / 3;
}


static int
_ged_draw_view_export_detail_geometry_from_obol_feature(
	struct ged_draw_view_export_detail *detail,
	void *view_ctx,
	const struct ged_draw_obol_view_feature_record *feature)
{
    if (!detail || !view_ctx || !feature || !feature->name)
	return 0;

    detail->geometry_kind =
	_ged_draw_view_export_geometry_kind_from_obol_feature(feature);
    if (detail->geometry_kind == GED_DRAW_VIEW_EXPORT_GEOMETRY_NONE)
	return 1;

    if (detail->geometry_kind == GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET ||
	    detail->geometry_kind == GED_DRAW_VIEW_EXPORT_GEOMETRY_POINT_SET) {
	if (feature->line_layer_parent_name) {
	    if (!ged_draw_obol_view_context_feature_layer_points_copy(
		    view_ctx, feature->line_layer_parent_name,
		    feature->line_layer_index, &detail->arrays.points,
		    &detail->arrays.point_count))
		return 0;
	} else {
	    if (!ged_draw_obol_view_context_feature_points_copy(view_ctx,
		    feature->name, &detail->arrays.points,
		    &detail->arrays.point_count))
		return 0;
	}
	detail->arrays.command_count = detail->arrays.point_count;
	if (detail->arrays.command_count) {
	    detail->arrays.commands = (int *)bu_calloc(
		    detail->arrays.command_count, sizeof(int),
		    "GED Obol feature export line commands");
	}
	size_t structure_count = 0;
	for (size_t i = 0; i < detail->arrays.command_count; i++) {
	    if (detail->geometry_kind == GED_DRAW_VIEW_EXPORT_GEOMETRY_POINT_SET)
		detail->arrays.commands[i] = GED_DRAW_VIEW_LINE_POINT_DRAW;
	    else if (feature->line_layer_parent_name &&
		    !ged_draw_obol_view_context_feature_layer_line_command_at(
		    view_ctx, feature->line_layer_parent_name,
		    feature->line_layer_index, i, &detail->arrays.commands[i]))
		detail->arrays.commands[i] =
		    (i % 2) ? GED_DRAW_VIEW_LINE_DRAW :
		    GED_DRAW_VIEW_LINE_MOVE;
	    else if (!feature->line_layer_parent_name &&
		    !ged_draw_obol_view_context_feature_line_command_at(
		    view_ctx, feature->name, i, &detail->arrays.commands[i]))
		detail->arrays.commands[i] =
		    (i % 2) ? GED_DRAW_VIEW_LINE_DRAW :
		    GED_DRAW_VIEW_LINE_MOVE;
	    if (detail->arrays.commands[i] != GED_DRAW_VIEW_LINE_MOVE)
		structure_count++;
	}
	detail->record.vlist_structure_count = structure_count;
	detail->record.vlist_point_count = detail->arrays.point_count;
	_ged_draw_view_export_record_bounds_from_points(&detail->record,
		detail->arrays.points, detail->arrays.point_count);
	return 1;
    }

    if (detail->geometry_kind ==
	    GED_DRAW_VIEW_EXPORT_GEOMETRY_INDEXED_FACE_SET) {
	if (!ged_draw_obol_view_context_feature_points_copy(view_ctx,
		feature->name, &detail->surface.points,
		&detail->surface.point_count))
	    return 0;
	if (!ged_draw_obol_view_context_feature_indices_copy(view_ctx,
		feature->name, &detail->surface.indices,
		&detail->surface.index_count))
	    return 0;
	detail->surface.normal_count = feature->normal_count;
	detail->surface.face_count =
	    _ged_draw_view_export_face_count_from_indices(
		    detail->surface.indices, detail->surface.index_count);
	detail->surface.material_valid = 1;
	detail->surface.material_draw_mode = GED_DRAW_MODE_SHADED;
	detail->surface.material_transparency = 0.0;
	detail->surface.material_highlighted = 0;
	detail->surface.material_color[0] = feature->color[0];
	detail->surface.material_color[1] = feature->color[1];
	detail->surface.material_color[2] = feature->color[2];
	detail->record.vlist_structure_count = detail->surface.face_count;
	detail->record.vlist_point_count = detail->surface.point_count;
	_ged_draw_view_export_record_bounds_from_points(&detail->record,
		detail->surface.points, detail->surface.point_count);
	return 1;
    }

    return 1;
}


static int
_ged_draw_view_export_detail_from_obol_feature(
	struct ged_draw_view_export_detail *detail,
	void *view_ctx,
	const struct ged_draw_obol_view_feature_record *feature)
{
    if (!detail || !view_ctx || !feature || !feature->name ||
	    !feature->name[0])
	return 0;

    _ged_draw_view_export_detail_init(detail);
    bu_vls_strcpy(&detail->path, feature->name);

    struct ged_draw_view_db_object_record *rec = &detail->record;
    rec->path = bu_vls_cstr(&detail->path);
    rec->type_name =
	_ged_draw_view_export_type_name_from_obol_feature(feature);
    rec->geometry_name =
	_ged_draw_view_export_geometry_name_from_obol_feature(feature);
    rec->draw_mode =
	(feature->kind == GED_DRAW_OBOL_VIEW_FEATURE_KIND_INDEXED_FACE_SET ||
	 feature->kind == GED_DRAW_OBOL_VIEW_FEATURE_KIND_POLYGON_OVERLAY) ?
	GED_DRAW_MODE_SHADED : GED_DRAW_MODE_WIRE;
    rec->transparency = 0.0;
    rec->evaluated_region = 0;
    rec->is_database_source = 0;
    rec->non_database_source = 1;
    rec->is_database_intent = 0;
    rec->is_local_source = feature->local;
    rec->is_view_source = 1;
    rec->highlighted = 0;
    rec->selected = ged_draw_view_context_selection_contains_path(view_ctx,
	    GED_DRAW_VIEW_SELECTION_SELECTED_PATH, feature->name) ? 1 : 0;
    rec->visible = feature->visible;
    rec->line_style = feature->line_style;
    rec->color[0] = feature->color[0];
    rec->color[1] = feature->color[1];
    rec->color[2] = feature->color[2];
    MAT_IDN(rec->model_mat);
    rec->cache_identity = _ged_draw_view_export_hash_cstr(feature->name);
    rec->source_identity = rec->cache_identity;
    rec->detail_token = (uintptr_t)detail;

    if (!_ged_draw_view_export_detail_geometry_from_obol_feature(detail,
	    view_ctx, feature)) {
	_ged_draw_view_export_detail_free(detail);
	return 0;
    }

    return 1;
}


struct ged_draw_obol_feature_export_ctx {
    void *view_ctx;
    ged_draw_view_db_object_record_cb cb;
    void *userdata;
    int keep_going;
};


struct ged_draw_obol_polygon_export_ctx {
    void *view_ctx;
    unsigned int query_flags;
    const char *glob;
    int draw_mode;
    ged_draw_view_db_object_record_cb cb;
    void *userdata;
    int keep_going;
};


static int
_ged_draw_view_export_polygon_record_cb(
	rt_view_polygon_ref UNUSED(ref),
	const struct rt_view_polygon_record *polygon,
	void *userdata)
{
    struct ged_draw_obol_polygon_export_ctx *ctx =
	(struct ged_draw_obol_polygon_export_ctx *)userdata;
    if (!ctx || !ctx->keep_going || !ctx->cb || !polygon ||
	    !polygon->name || !polygon->name[0])
	return 0;

    if ((ctx->query_flags & GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS) &&
	    !(ctx->query_flags & GED_DRAW_VIEW_EXPORT_QUERY_VIEW_OBJECTS))
	return 0;
    if (ctx->draw_mode != GED_DRAW_VIEW_EXPORT_DRAW_MODE_ANY &&
	    ctx->draw_mode != GED_DRAW_MODE_WIRE)
	return 1;
    if (!_ged_draw_view_export_glob_match(ctx->glob, polygon->name, NULL))
	return 1;

    struct ged_draw_view_export_detail detail;
    _ged_draw_view_export_detail_init(&detail);
    bu_vls_strcpy(&detail.path, polygon->name);
    detail.geometry_kind = GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET;

    struct ged_draw_view_db_object_record *rec = &detail.record;
    rec->path = bu_vls_cstr(&detail.path);
    rec->type_name = "polygon";
    rec->geometry_name = "line-set";
    rec->draw_mode = GED_DRAW_MODE_WIRE;
    rec->transparency = 0.0;
    rec->evaluated_region = 0;
    rec->is_database_source = 0;
    rec->non_database_source = 1;
    rec->is_database_intent = 0;
    rec->is_local_source =
	(ctx->query_flags & GED_DRAW_VIEW_EXPORT_QUERY_LOCAL_ONLY) ? 1 : 0;
    rec->is_view_source = 1;
    rec->highlighted = 0;
    rec->selected = ged_draw_view_context_selection_contains_path(
	    ctx->view_ctx, GED_DRAW_VIEW_SELECTION_SELECTED_PATH,
	    polygon->name) ? 1 : 0;
    rec->visible = 1;
    rec->line_style = 0;
    rec->color[0] = polygon->edge_color[0];
    rec->color[1] = polygon->edge_color[1];
    rec->color[2] = polygon->edge_color[2];
    MAT_IDN(rec->model_mat);
    rec->cache_identity = _ged_draw_view_export_hash_cstr(polygon->name);
    rec->source_identity = rec->cache_identity;
    rec->vlist_structure_count = polygon->contour_count;
    rec->vlist_point_count = polygon->point_count;
    VMOVE(rec->bounds_center, polygon->origin_point);
    rec->bounds_radius = 0.0;
    rec->has_bounds = 1;
    rec->detail_token = (uintptr_t)&detail;

    ctx->keep_going = ctx->cb(rec, ctx->userdata);
    _ged_draw_view_export_detail_free(&detail);
    return ctx->keep_going;
}


static void
_ged_draw_view_export_polygon_records_from_rt(
	void *view_ctx,
	unsigned int query_flags,
	const char *glob,
	int draw_mode,
	ged_draw_view_db_object_record_cb cb,
	void *userdata,
	int *keep_going)
{
    if (!keep_going || !*keep_going || !view_ctx || !cb)
	return;

    const int wants_db =
	(query_flags & GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS) ? 1 : 0;
    const int wants_view =
	(query_flags & GED_DRAW_VIEW_EXPORT_QUERY_VIEW_OBJECTS) ? 1 : 0;
    if (wants_db && !wants_view)
	return;

    struct ged_draw_obol_polygon_export_ctx ctx;
    ctx.view_ctx = view_ctx;
    ctx.query_flags = query_flags;
    ctx.glob = glob;
    ctx.draw_mode = draw_mode;
    ctx.cb = cb;
    ctx.userdata = userdata;
    ctx.keep_going = 1;
    rt_view_context_polygon_visit_records(view_ctx,
	    _ged_draw_view_export_polygon_record_cb, &ctx);
    *keep_going = ctx.keep_going;
}


static int
_ged_draw_view_export_feature_record_cb(
	const struct ged_draw_obol_view_feature_record *feature,
	void *userdata)
{
    struct ged_draw_obol_feature_export_ctx *ctx =
	(struct ged_draw_obol_feature_export_ctx *)userdata;
    if (!ctx || !ctx->keep_going || !ctx->cb || !feature)
	return 0;

    struct ged_draw_view_export_detail detail;
    if (_ged_draw_view_export_detail_from_obol_feature(&detail,
	    ctx->view_ctx, feature)) {
	ctx->keep_going = ctx->cb(&detail.record, ctx->userdata);
	_ged_draw_view_export_detail_free(&detail);
    }
    return ctx->keep_going;
}


static void
_ged_draw_view_export_feature_records_from_obol(
	void *view_ctx,
	unsigned int query_flags,
	const char *glob,
	ged_draw_view_db_object_record_cb cb,
	void *userdata,
	int *keep_going)
{
    if (!keep_going || !*keep_going || !view_ctx || !cb)
	return;

    struct ged_draw_obol_feature_export_ctx ctx;
    ctx.view_ctx = view_ctx;
    ctx.cb = cb;
    ctx.userdata = userdata;
    ctx.keep_going = 1;
    ged_draw_obol_view_context_feature_records_foreach(view_ctx,
	    query_flags, glob, _ged_draw_view_export_feature_record_cb, &ctx);
    *keep_going = ctx.keep_going;
}


static void
_ged_draw_view_export_records_from_obol(
	struct ged *gedp,
	void *view_ctx,
	unsigned int query_flags,
	unsigned int render_flags,
	const char *glob,
	int draw_mode,
	ged_draw_view_db_object_record_cb cb,
	void *userdata,
	int *keep_going)
{
    if (!keep_going || !*keep_going || !gedp || !cb)
	return;

    struct ged_draw_obol_scene_context_info root_info;
    memset(&root_info, 0, sizeof(root_info));
    if (!ged_draw_obol_scene_context_info_for_path(gedp, "/", &root_info))
	return;

    struct ged_draw_obol_view_export_ctx ctx;
    ctx.gedp = gedp;
    ctx.view_ctx = view_ctx;
    ctx.query_flags = query_flags;
    ctx.render_flags = render_flags;
    ctx.glob = glob;
    ctx.draw_mode = draw_mode;
    ctx.cb = cb;
    ctx.userdata = userdata;
    ctx.keep_going = 1;

    _ged_draw_view_export_visit_obol_info(&ctx, &root_info, NULL);
    *keep_going = ctx.keep_going;
    ged_draw_obol_scene_context_info_free(&root_info);

    _ged_draw_view_export_polygon_records_from_rt(view_ctx, query_flags,
	    glob, draw_mode, cb, userdata, keep_going);
    _ged_draw_view_export_feature_records_from_obol(view_ctx, query_flags,
	    glob, cb, userdata, keep_going);
}


static unsigned int
_ged_draw_view_export_query_flags_to_bsg(unsigned int flags)
{
    unsigned int bsg_flags = 0;
    if (flags & GED_DRAW_VIEW_EXPORT_QUERY_VISIBLE_ONLY)
	bsg_flags |= BSG_EXPORT_QUERY_VISIBLE_ONLY;
    if (flags & GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS)
	bsg_flags |= BSG_EXPORT_QUERY_DB_OBJECTS;
    if (flags & GED_DRAW_VIEW_EXPORT_QUERY_VIEW_OBJECTS)
	bsg_flags |= BSG_EXPORT_QUERY_VIEW_OBJECTS;
    if (flags & GED_DRAW_VIEW_EXPORT_QUERY_LOCAL_ONLY)
	bsg_flags |= BSG_EXPORT_QUERY_LOCAL_ONLY;
    return bsg_flags;
}


static unsigned int
_ged_draw_view_export_render_flags_to_bsg(unsigned int flags)
{
    unsigned int bsg_flags = 0;
    if (flags & GED_DRAW_VIEW_EXPORT_RENDER_VISIBLE_ONLY)
	bsg_flags |= BSG_RENDER_FLAG_VISIBLE_ONLY;
    if (flags & GED_DRAW_VIEW_EXPORT_RENDER_PAYLOAD_PREPARE)
	bsg_flags |= BSG_RENDER_FLAG_PAYLOAD_PREPARE;
    return bsg_flags;
}


void
ged_draw_view_context_foreach_export_record(
	void *view_ctx,
	unsigned int query_flags,
	unsigned int render_flags,
	const char *glob,
	int draw_mode,
	ged_draw_view_db_object_record_cb cb,
	void *userdata)
{
    struct bsg_view *v = (struct bsg_view *)view_ctx;

    if (!v || !cb)
	return;

    struct ged *gedp = (struct ged *)ged_view_context_user_data_get(view_ctx);
    int keep_going = 1;
    if (ged_draw_obol_scene_controller_is_owned(gedp)) {
	_ged_draw_view_export_records_from_obol(gedp, view_ctx, query_flags,
		render_flags, glob, draw_mode, cb, userdata, &keep_going);
	return;
    }

    struct bsg_export_request request;
    bsg_export_request_init(&request, v);
    request.query_flags = _ged_draw_view_export_query_flags_to_bsg(query_flags);
    request.render_flags =
	_ged_draw_view_export_render_flags_to_bsg(render_flags);
    request.glob = glob;
    if (draw_mode != GED_DRAW_VIEW_EXPORT_DRAW_MODE_ANY)
	request.draw_mode = draw_mode;

    struct bsg_export_result *result = bsg_export_query(&request);
    if (!result)
	return;

    for (size_t i = 0; i < bsg_export_result_count(result); i++) {
	const struct bsg_export_record *export_rec =
	    bsg_export_result_get(result, i);
	if (!export_rec)
	    continue;

	struct ged_draw_view_export_detail detail;
	if (!_ged_draw_view_export_detail_from_bsg(&detail, export_rec))
	    continue;

	keep_going = cb(&detail.record, userdata);
	_ged_draw_view_export_detail_free(&detail);
	if (!keep_going)
	    break;
    }

    bsg_export_result_free(result);
}

static int
ged_draw_rt_path_matches_prefix(const char *path, const char *prefix)
{
    size_t n;

    if (!path || !path[0] || !prefix || !prefix[0])
	return 0;

    while (*path == '/')
	path++;
    while (*prefix == '/')
	prefix++;
    n = strlen(prefix);
    return BU_STR_EQUAL(path, prefix) ||
	(strlen(path) > n && bu_strncmp(path, prefix, n) == 0 &&
	 path[n] == '/');
}

static const char *
ged_draw_rt_path_skip_leading_slash(const char *path)
{
    if (!path)
	return NULL;
    while (*path == '/')
	path++;
    return path;
}

static int
ged_draw_rt_path_matches_pattern(const char *path, const char *pattern)
{
    if (ged_draw_rt_path_matches_prefix(path, pattern))
	return 1;
    path = ged_draw_rt_path_skip_leading_slash(path);
    pattern = ged_draw_rt_path_skip_leading_slash(pattern);
    return (path && pattern && !bu_path_match(pattern, path, 0)) ? 1 : 0;
}

struct ged_draw_rt_record_match_ctx {
    const char *pattern;
    int found;
};

static int
ged_draw_rt_record_match_cb(
	const struct ged_draw_view_db_object_record *rec,
	void *userdata)
{
    struct ged_draw_rt_record_match_ctx *ctx =
	(struct ged_draw_rt_record_match_ctx *)userdata;
    if (!ctx || !rec || !rec->path)
	return 1;
    if (rec->visible && rec->is_database_source &&
	    ged_draw_rt_path_matches_pattern(rec->path, ctx->pattern)) {
	ctx->found = 1;
	return 0;
    }
    return 1;
}

struct ged_draw_rt_pick_ctx {
    void *view_ctx;
    const char *pattern;
    void *result;
    int count;
};

static int
ged_draw_rt_pick_record_cb(
	const struct ged_draw_view_db_object_record *rec,
	void *userdata)
{
    struct ged_draw_rt_pick_ctx *ctx =
	(struct ged_draw_rt_pick_ctx *)userdata;
    if (!ctx || !ctx->result || !rec || !rec->path)
	return 1;
    if (!rec->visible || !rec->is_database_source ||
	    !ged_draw_rt_path_matches_pattern(rec->path, ctx->pattern))
	return 1;

    const char *source_path = ged_draw_rt_path_matches_prefix(rec->path,
	    ctx->pattern) ? ged_draw_rt_path_skip_leading_slash(ctx->pattern) :
	ged_draw_rt_path_skip_leading_slash(rec->path);
    if (rt_view_pick_result_context_append_path(ctx->result, ctx->view_ctx,
	    0, 0, source_path, 0.0))
	ctx->count++;
    return 1;
}

struct ged_draw_rt_source_match_ctx {
    const char *prefix;
    int found;
};

static int
ged_draw_rt_source_match_cb(struct ged *UNUSED(gedp),
			    const char *path,
			    void *userdata)
{
    struct ged_draw_rt_source_match_ctx *ctx =
	(struct ged_draw_rt_source_match_ctx *)userdata;
    if (!ctx || !path)
	return 1;
    if (ged_draw_rt_path_matches_prefix(path, ctx->prefix)) {
	ctx->found = 1;
	return 0;
    }
    return 1;
}

static void *
ged_draw_rt_obol_pick_semantic_path(
	void *view_ctx,
	const char *path_pattern,
	void *data)
{
    struct ged *gedp = (struct ged *)data;
    if (!view_ctx || !path_pattern || !path_pattern[0] ||
	    !ged_draw_obol_scene_controller_is_owned(gedp))
	return NULL;

    void *result = rt_view_pick_result_context_create();
    if (!result)
	return NULL;

    struct ged_draw_rt_pick_ctx ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.view_ctx = view_ctx;
    ctx.pattern = path_pattern;
    ctx.result = result;
    ged_draw_view_context_foreach_export_record(view_ctx,
	    GED_DRAW_VIEW_EXPORT_QUERY_VISIBLE_ONLY |
	    GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS,
	    GED_DRAW_VIEW_EXPORT_RENDER_VISIBLE_ONLY |
	    GED_DRAW_VIEW_EXPORT_RENDER_PAYLOAD_PREPARE,
	    NULL, GED_DRAW_VIEW_EXPORT_DRAW_MODE_ANY,
	    ged_draw_rt_pick_record_cb, &ctx);

    if (ctx.count > 0)
	return result;

    rt_view_pick_result_context_free(result);
    return NULL;
}

static int
ged_draw_rt_obol_render_export_consistency(
	void *view_ctx,
	const char *drawn_prefix,
	struct rt_view_render_export_consistency *summary,
	void *data)
{
    struct ged *gedp = (struct ged *)data;
    if (summary)
	memset(summary, 0, sizeof(*summary));
    if (!view_ctx || !drawn_prefix || !drawn_prefix[0] || !summary ||
	    !ged_draw_obol_scene_controller_is_owned(gedp))
	return 0;

    struct ged_draw_rt_record_match_ctx record_ctx;
    memset(&record_ctx, 0, sizeof(record_ctx));
    record_ctx.pattern = drawn_prefix;
    ged_draw_view_context_foreach_export_record(view_ctx,
	    GED_DRAW_VIEW_EXPORT_QUERY_VISIBLE_ONLY |
	    GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS,
	    GED_DRAW_VIEW_EXPORT_RENDER_VISIBLE_ONLY |
	    GED_DRAW_VIEW_EXPORT_RENDER_PAYLOAD_PREPARE,
	    NULL, GED_DRAW_VIEW_EXPORT_DRAW_MODE_ANY,
	    ged_draw_rt_record_match_cb, &record_ctx);

    struct ged_draw_rt_source_match_ctx source_ctx;
    memset(&source_ctx, 0, sizeof(source_ctx));
    source_ctx.prefix = drawn_prefix;
    (void)ged_draw_obol_database_source_paths_foreach(gedp, 1,
	    ged_draw_rt_source_match_cb, &source_ctx);

    summary->export_record_found = record_ctx.found;
    summary->render_item_found = source_ctx.found;
    summary->backend_node_found = source_ctx.found;
    summary->export_render_consistent = record_ctx.found && source_ctx.found;
    summary->export_backend_consistent = record_ctx.found && source_ctx.found;
    return 1;
}

int
ged_draw_view_context_obol_scene_adapter_attach(
	struct ged *gedp,
	void *view_ctx)
{
    if (!gedp || !view_ctx)
	return 0;

    struct rt_view_context_scene_adapter adapter;
    memset(&adapter, 0, sizeof(adapter));
    adapter.pick_semantic_path = ged_draw_rt_obol_pick_semantic_path;
    adapter.render_export_consistency =
	ged_draw_rt_obol_render_export_consistency;
    adapter.data = gedp;
    return rt_view_context_scene_adapter_set(view_ctx, &adapter);
}


static const char *
ged_draw_scene_ref_name(bsg_scene_ref ref)
{
    return bsg_scene_name(ref);
}


static const struct db_full_path *
ged_draw_scene_ref_fullpath(bsg_scene_ref ref)
{
    ged_draw_shape_state *bd = _ged_draw_shape_state_get_scene_ref(ref);
    if (!bd || bd->s_fullpath.fp_len <= 0)
        return NULL;
    return &bd->s_fullpath;
}


static struct directory *
ged_draw_scene_ref_leaf_dp(bsg_scene_ref ref)
{
    ged_draw_shape_state *bd = _ged_draw_shape_state_get_scene_ref(ref);
    if (!bd)
	return RT_DIR_NULL;
    if (bd->leaf_dp)
	return bd->leaf_dp;
    if (bd->s_fullpath.fp_len > 0)
	return DB_FULL_PATH_CUR_DIR(&bd->s_fullpath);
    return RT_DIR_NULL;
}


static ged_draw_shape_state *
ged_draw_scene_ref_shape_state(bsg_scene_ref ref)
{
    return _ged_draw_shape_state_get_scene_ref(ref);
}


static void
_ged_draw_scene_ref_set_source_ref(bsg_scene_ref ref, bsg_scene_ref source_ref)
{
    ged_draw_shape_state *shape_data = ged_draw_scene_ref_shape_state(ref);
    if (shape_data)
	shape_data->source_ref = ged_draw_scene_ref_to_rt_view_ref(source_ref);
}


static bsg_scene_ref
ged_draw_scene_ref_null(void)
{
    return bsg_scene_ref_null();
}


static bsg_scene_ref
ged_draw_scene_ref_from_context(void *scene_ctx)
{
    bsg_scene_ref ref = bsg_scene_ref_null();
    ref.opaque = scene_ctx;
    return ref;
}


static void *
ged_draw_scene_ref_context(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? NULL : ref.opaque;
}

#define GED_DRAW_OBOL_CONTEXT_TOKEN_MAGIC 0x47444f626f6c4354ULL

struct ged_draw_obol_context_token {
    uint64_t magic;
    uint64_t token;
    struct ged *gedp;
    char *path;
    char *name;
    void *parent_ctx;
    int node_kind;
    int is_group;
    int is_shape;
    int is_database_source;
    int has_parent;
    int draw_tree_depth;
    size_t child_count;
    int draw_mode_valid;
    int draw_mode;
    struct db_full_path fullpath;
    int fullpath_valid;
    struct bn_tol obol_snapshot_tol;
    struct bg_tess_tol obol_snapshot_ttol;
};


static struct ged_draw_obol_context_token *
ged_draw_obol_context_from_scene_ctx(void *scene_ctx)
{
    struct ged_draw_obol_context_token *token =
	(struct ged_draw_obol_context_token *)scene_ctx;
    if (!token || token->magic != GED_DRAW_OBOL_CONTEXT_TOKEN_MAGIC)
	return NULL;
    return token;
}


static struct ged_draw_obol_context_token *
ged_draw_obol_context_from_rt_ref(rt_view_scene_ref scene_ref)
{
    if (rt_view_scene_ref_backend(scene_ref) != RT_VIEW_SCENE_BACKEND_OBOL)
	return NULL;
    return ged_draw_obol_context_from_scene_ctx(
	    rt_view_scene_ref_context(scene_ref));
}


static ged_draw_group_ref ged_draw_obol_group_ref_for_path(struct ged *gedp,
							   const char *path);


static ged_draw_group_ref
ged_draw_obol_context_group_ref(struct ged *gedp,
				struct ged_draw_obol_context_token *token)
{
    if (!gedp || !token || !token->is_group || token->is_database_source)
	return GED_DRAW_GROUP_REF_NULL;

    rt_view_scene_ref group_scene_ref =
	rt_view_scene_ref_make((void *)token, RT_VIEW_SCENE_BACKEND_OBOL);
    ged_draw_group_ref ref =
	ged_draw_registry_group_ref_from_source_ref(gedp, group_scene_ref);
    ref.scene_revision = 0;
    if (!ged_draw_group_ref_is_null(ref) && token->fullpath_valid)
	(void)ged_draw_registry_group_ref_set_indexed_fullpath(gedp, ref,
		&token->fullpath);
    return ref;
}


static ged_draw_group_ref
ged_draw_obol_top_group_ref_for_path(struct ged *gedp, const char *path)
{
    if (!gedp || !path || !path[0])
	return GED_DRAW_GROUP_REF_NULL;

    const char *start = ged_draw_dbpath_skip_lead_slash(path);
    if (!start || !start[0])
	return GED_DRAW_GROUP_REF_NULL;

    const char *slash = strchr(start, '/');
    if (!slash)
	return ged_draw_obol_group_ref_for_path(gedp, start);

    struct bu_vls top = BU_VLS_INIT_ZERO;
    bu_vls_strncpy(&top, start, (size_t)(slash - start));
    ged_draw_group_ref ref = ged_draw_obol_group_ref_for_path(gedp,
	    bu_vls_cstr(&top));
    bu_vls_free(&top);
    return ref;
}


static ged_draw_group_ref
ged_draw_obol_top_group_ref_for_fullpath(struct ged *gedp,
					 const struct db_full_path *fullpath)
{
    if (!gedp || !fullpath || fullpath->fp_len <= 0 ||
	    !fullpath->fp_names[0] || !fullpath->fp_names[0]->d_namep)
	return GED_DRAW_GROUP_REF_NULL;
    return ged_draw_obol_group_ref_for_path(gedp,
	    fullpath->fp_names[0]->d_namep);
}


static ged_draw_group_ref
ged_draw_obol_shape_owner_group_ref(struct ged *gedp,
				    struct ged_draw_obol_context_token *token)
{
    if (!gedp || !token)
	return GED_DRAW_GROUP_REF_NULL;

    struct ged_draw_obol_context_token *top_group = NULL;
    for (struct ged_draw_obol_context_token *ancestor =
	    ged_draw_obol_context_from_scene_ctx(token->parent_ctx);
	    ancestor;
	    ancestor = ged_draw_obol_context_from_scene_ctx(
		    ancestor->parent_ctx)) {
	if (ancestor->is_database_source && ancestor->path) {
	    ged_draw_group_ref ref = ancestor->fullpath_valid ?
		ged_draw_obol_top_group_ref_for_fullpath(gedp,
			&ancestor->fullpath) :
		ged_draw_obol_top_group_ref_for_path(gedp, ancestor->path);
	    if (!ged_draw_group_ref_is_null(ref))
		return ref;
	}
	if (ancestor->is_group && !ancestor->is_database_source &&
		ancestor->path && !BU_STR_EQUAL(ancestor->path, "/"))
	    top_group = ancestor;
    }

    return ged_draw_obol_context_group_ref(gedp, top_group);
}


static void
ged_draw_obol_context_token_free(struct ged_draw_obol_context_token *token)
{
    if (!token)
	return;
    token->magic = 0;
    if (token->path)
	bu_free(token->path, "GED Obol context token path");
    if (token->name)
	bu_free(token->name, "GED Obol context token name");
    if (token->fullpath_valid)
	db_free_full_path(&token->fullpath);
    BU_PUT(token, struct ged_draw_obol_context_token);
}


static char *
ged_draw_obol_db_path_without_instance_suffixes(const char *path)
{
    if (!path)
	return NULL;

    char *lookup_path = (char *)bu_malloc(strlen(path) + 1,
	    "GED Obol DB lookup path");
    char *out = lookup_path;
    for (const char *cp = path; *cp; cp++) {
	if (*cp == '@' && cp[1] >= '0' && cp[1] <= '9') {
	    while (cp[1] && cp[1] != '/')
		cp++;
	    continue;
	}
	*out++ = *cp;
    }
    *out = '\0';
    return lookup_path;
}


static int
ged_draw_obol_context_db_path_components_exist(
	struct db_i *dbip,
	const char *path)
{
    if (!dbip || !path || !path[0])
	return 0;

    const char *cp = path;
    while (*cp == '/')
	cp++;

    int found_component = 0;
    while (*cp) {
	while (*cp == '/')
	    cp++;
	if (!*cp)
	    break;

	const char *slash = strchr(cp, '/');
	size_t len = slash ? (size_t)(slash - cp) : strlen(cp);
	if (len == 0) {
	    cp = slash ? slash + 1 : cp + len;
	    continue;
	}

	char *component = (char *)bu_malloc(len + 1,
		"GED Obol context DB path component");
	memcpy(component, cp, len);
	component[len] = '\0';
	struct directory *dp = db_lookup(dbip, component, LOOKUP_QUIET);
	bu_free(component, "GED Obol context DB path component");
	if (dp == RT_DIR_NULL)
	    return 0;
	found_component = 1;

	if (!slash)
	    break;
	cp = slash + 1;
    }

    return found_component;
}


void
ged_draw_obol_context_tokens_free(struct ged *gedp)
{
    if (!gedp || !gedp->i || !gedp->i->ged_gdp ||
	    !gedp->i->ged_gdp->gd_obol_context_tokens_init)
	return;

    struct ged_drawable *gdp = gedp->i->ged_gdp;
    for (size_t i = 0; i < BU_PTBL_LEN(&gdp->gd_obol_context_tokens); i++) {
	struct ged_draw_obol_context_token *token =
	    (struct ged_draw_obol_context_token *)BU_PTBL_GET(
		    &gdp->gd_obol_context_tokens, i);
	ged_draw_obol_context_token_free(token);
    }
    bu_ptbl_free(&gdp->gd_obol_context_tokens);
    gdp->gd_obol_context_tokens_init = 0;
    gdp->gd_obol_next_context_token = 1;
}


static void
ged_draw_obol_context_token_fullpath_init(
	struct ged_draw_obol_context_token *token)
{
    if (!token)
	return;
    db_full_path_init(&token->fullpath);
    token->fullpath_valid = 0;
    if (!token->gedp || !token->gedp->dbip || !token->path ||
	    !token->path[0] || BU_STR_EQUAL(token->path, "/"))
	return;
    char *db_lookup_path =
	ged_draw_obol_db_path_without_instance_suffixes(token->path);
    if (!db_lookup_path || !db_lookup_path[0]) {
	if (db_lookup_path)
	    bu_free(db_lookup_path, "GED Obol DB lookup path");
	return;
    }
    if (!ged_draw_obol_context_db_path_components_exist(
	    token->gedp->dbip, db_lookup_path)) {
	bu_free(db_lookup_path, "GED Obol DB lookup path");
	return;
    }

    if (db_string_to_path(&token->fullpath, token->gedp->dbip,
	    db_lookup_path) == 0)
	token->fullpath_valid = 1;
    else {
	db_free_full_path(&token->fullpath);
	db_full_path_init(&token->fullpath);
    }
    bu_free(db_lookup_path, "GED Obol DB lookup path");
}


static void *
ged_draw_obol_context_token_create(
	struct ged *gedp,
	const struct ged_draw_obol_scene_context_info *info,
	void *parent_ctx)
{
    if (!gedp || !gedp->i || !gedp->i->ged_gdp || !info || !info->path)
	return NULL;

    struct ged_drawable *gdp = gedp->i->ged_gdp;
    if (!gdp->gd_obol_context_tokens_init) {
	BU_PTBL_INIT(&gdp->gd_obol_context_tokens);
	gdp->gd_obol_context_tokens_init = 1;
	gdp->gd_obol_next_context_token = 1;
    }

    struct ged_draw_obol_context_token *token;
    BU_GET(token, struct ged_draw_obol_context_token);
    memset(token, 0, sizeof(*token));
    token->magic = GED_DRAW_OBOL_CONTEXT_TOKEN_MAGIC;
    token->token = gdp->gd_obol_next_context_token++;
    token->gedp = gedp;
    token->path = bu_strdup(info->path);
    token->name = bu_strdup(info->name ? info->name : info->path);
    token->parent_ctx = parent_ctx;
    token->node_kind = info->node_kind;
    token->is_group = info->is_group;
    token->is_shape = info->is_shape;
    token->is_database_source = info->is_database_source;
    token->has_parent = info->has_parent;
    token->draw_tree_depth = info->draw_tree_depth;
    token->child_count = info->child_count;
    token->draw_mode_valid = info->draw_mode_valid;
    token->draw_mode = info->draw_mode;
    BN_TOL_INIT_SET_TOL(&token->obol_snapshot_tol);
    BG_TESS_TOL_INIT_SET_TOL(&token->obol_snapshot_ttol);
    ged_draw_obol_context_token_fullpath_init(token);
    bu_ptbl_ins(&gdp->gd_obol_context_tokens, (long *)token);
    return (void *)token;
}


static void *
ged_draw_obol_context_token_find(
	struct ged *gedp,
	const struct ged_draw_obol_scene_context_info *info,
	void *parent_ctx)
{
    if (!gedp || !gedp->i || !gedp->i->ged_gdp || !info || !info->path ||
	    !gedp->i->ged_gdp->gd_obol_context_tokens_init)
	return NULL;

    struct ged_drawable *gdp = gedp->i->ged_gdp;
    for (size_t i = 0; i < BU_PTBL_LEN(&gdp->gd_obol_context_tokens); i++) {
	struct ged_draw_obol_context_token *token =
	    (struct ged_draw_obol_context_token *)BU_PTBL_GET(
		    &gdp->gd_obol_context_tokens, i);
	if (!token || token->magic != GED_DRAW_OBOL_CONTEXT_TOKEN_MAGIC ||
		!token->path || !BU_STR_EQUAL(token->path, info->path) ||
		token->node_kind != info->node_kind ||
		token->is_group != info->is_group ||
		token->is_shape != info->is_shape ||
		token->is_database_source != info->is_database_source ||
		token->draw_mode_valid != info->draw_mode_valid ||
		(token->draw_mode_valid &&
		 token->draw_mode != info->draw_mode))
	    continue;

	if (parent_ctx && !token->parent_ctx)
	    token->parent_ctx = parent_ctx;
	token->has_parent = info->has_parent;
	token->draw_tree_depth = info->draw_tree_depth;
	token->child_count = info->child_count;
	token->draw_mode_valid = info->draw_mode_valid;
	token->draw_mode = info->draw_mode;
	if (!token->fullpath_valid)
	    ged_draw_obol_context_token_fullpath_init(token);
	return (void *)token;
    }

    return NULL;
}


static void *
ged_draw_obol_context_token_for_path_mode(
	struct ged *gedp,
	const char *path,
	void *parent_ctx,
	int draw_mode_valid,
	int draw_mode)
{
    if (!gedp || !path || !path[0])
	return NULL;

    struct ged_draw_obol_scene_context_info info;
    memset(&info, 0, sizeof(info));
    if (!ged_draw_obol_scene_context_info_for_path(gedp, path, &info))
	return NULL;
    info.draw_mode_valid = draw_mode_valid ? 1 : 0;
    info.draw_mode = draw_mode;

    void *ctx = ged_draw_obol_context_token_find(gedp, &info, parent_ctx);
    if (!ctx)
	ctx = ged_draw_obol_context_token_create(gedp, &info, parent_ctx);
    ged_draw_obol_scene_context_info_free(&info);
    return ctx;
}


static void *
ged_draw_obol_context_token_for_path(
	struct ged *gedp,
	const char *path,
	void *parent_ctx)
{
    return ged_draw_obol_context_token_for_path_mode(gedp, path, parent_ctx,
	    0, GED_DRAW_MODE_WIRE);
}


static char *
ged_draw_obol_context_parent_path(const char *path)
{
    if (!path || !path[0] || BU_STR_EQUAL(path, "/"))
	return NULL;

    const char *start = path;
    while (*start == '/')
	start++;
    if (!start[0])
	return NULL;

    size_t len = strlen(start);
    while (len > 0 && start[len - 1] == '/')
	len--;
    if (len == 0)
	return NULL;

    const char *slash = NULL;
    for (size_t i = 0; i < len; i++) {
	if (start[i] == '/')
	    slash = &start[i];
    }

    if (!slash)
	return bu_strdup("/");

    size_t parent_len = (size_t)(slash - start);
    if (parent_len == 0)
	return bu_strdup("/");

    char *parent = (char *)bu_malloc(parent_len + 1,
	    "GED Obol context parent path");
    memcpy(parent, start, parent_len);
    parent[parent_len] = '\0';
    return parent;
}


static void *
ged_draw_obol_context_parent_token_for_path(
	struct ged *gedp,
	const char *path)
{
    char *parent_path = ged_draw_obol_context_parent_path(path);
    if (!parent_path)
	return NULL;

    void *grandparent_ctx =
	ged_draw_obol_context_parent_token_for_path(gedp, parent_path);
    void *parent_ctx = ged_draw_obol_context_token_for_path(gedp,
	    parent_path, grandparent_ctx);
    bu_free(parent_path, "GED Obol context parent path");
    return parent_ctx;
}


static void *ged_draw_obol_group_context_parent_token_for_path(
	struct ged *gedp,
	const char *path);
static void *ged_draw_obol_group_context_token_for_path(
	struct ged *gedp,
	const char *path,
	void *parent_ctx);


static int
ged_draw_obol_path_equal_normalized(const char *a, const char *b)
{
    if (!a || !b)
	return 0;

    a = ged_draw_dbpath_skip_lead_slash(a);
    b = ged_draw_dbpath_skip_lead_slash(b);

    size_t alen = strlen(a);
    size_t blen = strlen(b);
    while (alen > 0 && a[alen - 1] == '/')
	alen--;
    while (blen > 0 && b[blen - 1] == '/')
	blen--;

    return alen == blen && bu_strncmp(a, b, alen) == 0;
}


static void *
ged_draw_obol_database_source_parent_token_for_path(
	struct ged *gedp,
	const char *path)
{
    struct bu_vls owner_group_path = BU_VLS_INIT_ZERO;
    if (ged_draw_obol_database_source_owner_group_path_for_path(gedp,
	    path, &owner_group_path) && bu_vls_strlen(&owner_group_path) > 0) {
	const char *owner_path = bu_vls_cstr(&owner_group_path);
	if (BU_STR_EQUAL(owner_path, "/")) {
	    bu_vls_free(&owner_group_path);
	    return ged_draw_obol_context_token_for_path(gedp, "/", NULL);
	}
	if (!ged_draw_obol_path_equal_normalized(owner_path, path)) {
	    void *grandparent_ctx =
		ged_draw_obol_group_context_parent_token_for_path(gedp,
			owner_path);
	    void *parent_ctx = ged_draw_obol_group_context_token_for_path(
		    gedp, owner_path, grandparent_ctx);
	    bu_vls_free(&owner_group_path);
	    return parent_ctx;
	}
    }
    bu_vls_free(&owner_group_path);

    return ged_draw_obol_context_parent_token_for_path(gedp, path);
}


static void *
ged_draw_obol_context_token_with_parent_for_path_mode(
	struct ged *gedp,
	const char *path,
	int draw_mode_valid,
	int draw_mode)
{
    void *parent_ctx =
	ged_draw_obol_database_source_parent_token_for_path(gedp, path);
    return ged_draw_obol_context_token_for_path_mode(gedp, path, parent_ctx,
	    draw_mode_valid, draw_mode);
}


static void *
ged_draw_obol_context_token_with_parent_for_path(
	struct ged *gedp,
	const char *path)
{
    return ged_draw_obol_context_token_with_parent_for_path_mode(gedp, path,
	    0, GED_DRAW_MODE_WIRE);
}


static char *
ged_draw_obol_context_leaf_name_dup(const char *path)
{
    if (!path || !path[0])
	return bu_strdup("");
    if (BU_STR_EQUAL(path, "/"))
	return bu_strdup("/");

    const char *start = path;
    while (*start == '/')
	start++;
    const char *slash = strrchr(start, '/');
    return bu_strdup((slash && slash[1]) ? slash + 1 : start);
}


static int
ged_draw_obol_group_context_depth(const char *path)
{
    if (!path || !path[0] || BU_STR_EQUAL(path, "/"))
	return 0;

    int depth = 0;
    const char *cp = path;
    while (*cp) {
	while (*cp == '/')
	    cp++;
	if (!*cp)
	    break;
	depth++;
	while (*cp && *cp != '/')
	    cp++;
    }
    return depth;
}


static void *ged_draw_obol_group_context_parent_token_for_path(
	struct ged *gedp,
	const char *path);


static void *
ged_draw_obol_group_context_token_for_path(
	struct ged *gedp,
	const char *path,
	void *parent_ctx)
{
    if (!gedp || !path || !path[0])
	return NULL;

    struct ged_draw_group_record_summary group_summary;
    memset(&group_summary, 0, sizeof(group_summary));
    if (!ged_draw_obol_group_record_summary_for_path(gedp, path,
	    &group_summary))
	return NULL;

    struct ged_draw_obol_scene_context_info info;
    memset(&info, 0, sizeof(info));
    info.path = bu_strdup(path);
    info.name = ged_draw_obol_context_leaf_name_dup(path);
    info.node_kind = 1;
    info.is_group = 1;
    info.is_shape = 0;
    info.is_database_source = 0;
    info.has_parent = BU_STR_EQUAL(path, "/") ? 0 : 1;
    info.draw_tree_depth = ged_draw_obol_group_context_depth(path);

    size_t child_count = 0;
    if (ged_draw_obol_group_child_count_for_path(gedp, path, &child_count))
	info.child_count = child_count;

    void *ctx = ged_draw_obol_context_token_find(gedp, &info, parent_ctx);
    if (!ctx)
	ctx = ged_draw_obol_context_token_create(gedp, &info, parent_ctx);
    ged_draw_obol_scene_context_info_free(&info);
    return ctx;
}


static void *
ged_draw_obol_group_context_parent_token_for_path(
	struct ged *gedp,
	const char *path)
{
    char *parent_path = ged_draw_obol_context_parent_path(path);
    if (!parent_path)
	return NULL;

    void *parent_ctx = NULL;
    if (BU_STR_EQUAL(parent_path, "/")) {
	parent_ctx = ged_draw_obol_context_token_for_path(gedp, parent_path,
		NULL);
    } else {
	void *grandparent_ctx =
	    ged_draw_obol_group_context_parent_token_for_path(gedp,
		    parent_path);
	parent_ctx = ged_draw_obol_group_context_token_for_path(gedp,
		parent_path, grandparent_ctx);
    }

    bu_free(parent_path, "GED Obol group context parent path");
    return parent_ctx;
}


static void *
ged_draw_obol_group_context_token_with_parent_for_path(
	struct ged *gedp,
	const char *path)
{
    void *parent_ctx =
	ged_draw_obol_group_context_parent_token_for_path(gedp, path);
    return ged_draw_obol_group_context_token_for_path(gedp, path,
	    parent_ctx);
}


static ged_draw_group_ref
ged_draw_obol_group_ref_for_path(struct ged *gedp, const char *path)
{
    void *scene_ctx = ged_draw_obol_group_context_token_with_parent_for_path(
	    gedp, path);
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_ctx(scene_ctx);
    if (!token || !token->is_group || token->is_database_source)
	return GED_DRAW_GROUP_REF_NULL;

    rt_view_scene_ref scene_ref = rt_view_scene_ref_make(scene_ctx,
	    RT_VIEW_SCENE_BACKEND_OBOL);
    ged_draw_group_ref ref =
	ged_draw_registry_group_ref_from_source_ref(gedp, scene_ref);
    if (!ged_draw_group_ref_is_null(ref) && token->fullpath_valid)
	(void)ged_draw_registry_group_ref_set_indexed_fullpath(gedp, ref,
		&token->fullpath);
    return ref;
}


static int
ged_draw_obol_context_token_tree_summary(
	struct ged_draw_obol_context_token *token,
	const struct db_full_path *fallback_fullpath,
	struct ged_draw_scene_tree_summary *out)
{
    if (!token || !out)
	return 0;

    memset(out, 0, sizeof(*out));
    out->valid = 1;
    out->is_group = token->is_group;
    out->is_shape = token->is_shape;
    out->has_parent = token->has_parent;
    out->name = token->name;
    out->fullpath = token->fullpath_valid ? &token->fullpath :
	fallback_fullpath;
    out->draw_tree_depth = token->draw_tree_depth;
    out->child_count = token->child_count;

    size_t obol_child_count = 0;
    if (ged_draw_obol_group_child_count_for_path(token->gedp,
	    token->path, &obol_child_count))
	out->child_count = obol_child_count;
    if (token->is_database_source &&
	    ged_draw_obol_database_source_needs_primary_shape_proxy(
		token->gedp, token->path))
	out->child_count++;
    return 1;
}


static int
ged_draw_obol_group_tree_summary_for_path(
	struct ged *gedp,
	const char *path,
	const struct db_full_path *fallback_fullpath,
	struct ged_draw_scene_tree_summary *out)
{
    void *scene_ctx = ged_draw_obol_group_context_token_with_parent_for_path(
	    gedp, path);
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_ctx(scene_ctx);
    if (!token || !token->is_group || token->is_database_source)
	return 0;

    return ged_draw_obol_context_token_tree_summary(token,
	    fallback_fullpath, out);
}


static int
ged_draw_obol_path_equal(const char *a, const char *b)
{
    if (!a || !b)
	return 0;
    if (BU_STR_EQUAL(a, b))
	return 1;
    return BU_STR_EQUAL(ged_draw_dbpath_skip_lead_slash(a),
	    ged_draw_dbpath_skip_lead_slash(b));
}


static int
ged_draw_obol_database_source_has_primary_shape_child(
	struct ged *gedp,
	const char *path)
{
    struct ged_draw_obol_scene_context_info info;
    memset(&info, 0, sizeof(info));
    int has_primary = 0;

    if (!gedp || !path || !path[0] ||
	    !ged_draw_obol_scene_child_context_info_for_path(gedp, path, 0,
		&info))
	return 0;

    has_primary = info.is_shape && info.is_database_source &&
	ged_draw_obol_path_equal(info.path, path);
    ged_draw_obol_scene_context_info_free(&info);
    return has_primary;
}


static int
ged_draw_obol_database_source_needs_primary_shape_proxy(
	struct ged *gedp,
	const char *path)
{
    struct ged_draw_database_source_summary source_summary;
    if (!gedp || !path || !path[0])
	return 0;

    memset(&source_summary, 0, sizeof(source_summary));
    if (!ged_draw_obol_database_source_summary_for_path(gedp, path,
	    &source_summary) || !source_summary.valid)
	return 0;

    return ged_draw_obol_database_source_has_primary_shape_child(gedp, path) ?
	0 : 1;
}


static int
ged_draw_obol_database_source_primary_shape_info(
	struct ged *gedp,
	const char *path,
	int draw_tree_depth,
	struct ged_draw_obol_scene_context_info *out)
{
    if (!out)
	return 0;
    memset(out, 0, sizeof(*out));
    if (!gedp || !path || !path[0] ||
	    !ged_draw_obol_database_source_needs_primary_shape_proxy(gedp,
		path))
	return 0;

    const char *leaf = strrchr(path, '/');
    leaf = (leaf && leaf[1]) ? leaf + 1 : ged_draw_dbpath_skip_lead_slash(path);
    out->path = bu_strdup(path);
    out->name = bu_strdup((leaf && leaf[0]) ? leaf : path);
    out->node_kind = 0;
    out->is_group = 0;
    out->is_shape = 1;
    out->is_database_source = 1;
    out->has_parent = 1;
    out->draw_tree_depth = draw_tree_depth;
    out->child_count = 0;
    return 1;
}


static int
ged_draw_obol_database_source_container_info(
	struct ged *gedp,
	const char *path,
	int draw_tree_depth,
	struct ged_draw_obol_scene_context_info *out)
{
    struct ged_draw_database_source_summary source_summary;
    if (!out)
	return 0;
    memset(out, 0, sizeof(*out));
    if (!gedp || !path || !path[0])
	return 0;

    memset(&source_summary, 0, sizeof(source_summary));
    if (!ged_draw_obol_database_source_summary_for_path(gedp, path,
	    &source_summary) || !source_summary.valid)
	return 0;

    const char *leaf = strrchr(path, '/');
    leaf = (leaf && leaf[1]) ? leaf + 1 : ged_draw_dbpath_skip_lead_slash(path);
    out->path = bu_strdup(path);
    out->name = bu_strdup((leaf && leaf[0]) ? leaf : path);
    out->node_kind = 0;
    out->is_group = 0;
    out->is_shape = 0;
    out->is_database_source = 1;
    out->has_parent = 1;
    out->draw_tree_depth = draw_tree_depth;

    size_t child_count = 0;
    while (1) {
	struct ged_draw_obol_scene_context_info child_info;
	memset(&child_info, 0, sizeof(child_info));
	if (!ged_draw_obol_scene_child_context_info_for_path(gedp, path,
		child_count, &child_info))
	    break;
	ged_draw_obol_scene_context_info_free(&child_info);
	child_count++;
    }
    out->child_count = child_count;
    return 1;
}


static rt_view_scene_ref
ged_draw_scene_ref_to_rt_view_ref(bsg_scene_ref ref)
{
    return rt_view_scene_ref_make(ged_draw_scene_ref_context(ref),
	    RT_VIEW_SCENE_BACKEND_BSG);
}


static bsg_scene_ref
ged_draw_scene_ref_from_rt_view_ref(rt_view_scene_ref ref)
{
    if (rt_view_scene_ref_backend(ref) != RT_VIEW_SCENE_BACKEND_NONE &&
	    rt_view_scene_ref_backend(ref) != RT_VIEW_SCENE_BACKEND_BSG)
	return ged_draw_scene_ref_null();

    return ged_draw_scene_ref_from_context(rt_view_scene_ref_context(ref));
}


static rt_view_scene_ref
ged_draw_scene_context_rt_ref(void *scene_ctx)
{
    if (ged_draw_obol_context_from_scene_ctx(scene_ctx))
	return rt_view_scene_ref_make(scene_ctx, RT_VIEW_SCENE_BACKEND_OBOL);

    return ged_draw_scene_ref_to_rt_view_ref(
	    ged_draw_scene_ref_from_context(scene_ctx));
}


static int
ged_draw_scene_ref_is_null(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref);
}


static int
ged_draw_scene_ref_equal(bsg_scene_ref a, bsg_scene_ref b)
{
    return bsg_scene_ref_equal(a, b);
}


static ged_draw_shape_state *
_ged_draw_shape_state_get_scene_ref(bsg_scene_ref ref)
{
    return ged_draw_shape_state_get_scene_rt_ref(
	    ged_draw_scene_ref_to_rt_view_ref(ref));
}


static ged_draw_shape_ref
_ged_draw_shape_ref_from_scene_ref(struct ged *gedp, bsg_scene_ref ref)
{
    return ged_draw_registry_shape_ref_from_source_ref(gedp,
	    ged_draw_scene_ref_to_rt_view_ref(ref));
}


static ged_draw_group_ref
_ged_draw_group_ref_from_scene_ref(struct ged *gedp, bsg_scene_ref ref)
{
    return ged_draw_registry_group_ref_from_source_ref(gedp,
	    ged_draw_scene_ref_to_rt_view_ref(ref));
}


static bsg_scene_ref
_ged_draw_shape_ref_scene_ref(struct ged *gedp, ged_draw_shape_ref ref)
{
    return ged_draw_scene_ref_from_rt_view_ref(
	    ged_draw_registry_shape_ref_rt_ref(gedp, ref));
}


static bsg_scene_ref
_ged_draw_shape_ref_runtime_scene_ref(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	shape_ref = ged_draw_scene_ref_from_rt_view_ref(
		ged_draw_registry_shape_ref_cache_rt_ref(gedp, ref));
    return shape_ref;
}


static struct ged_draw_obol_context_token *
ged_draw_shape_ref_obol_token(struct ged *gedp, ged_draw_shape_ref ref)
{
    rt_view_scene_ref scene_ref = ged_draw_registry_shape_ref_rt_ref(gedp,
	    ref);
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_rt_ref(scene_ref);
    if (token)
	return token;

    scene_ref = ged_draw_registry_shape_ref_cache_rt_ref(gedp, ref);
    return ged_draw_obol_context_from_rt_ref(scene_ref);
}


static struct ged_draw_obol_context_token *
ged_draw_group_ref_obol_token(struct ged *gedp, ged_draw_group_ref ref)
{
    rt_view_scene_ref scene_ref = ged_draw_registry_group_ref_rt_ref(gedp,
	    ref);
    return ged_draw_obol_context_from_rt_ref(scene_ref);
}


static bsg_scene_ref
_ged_draw_group_ref_scene_ref(struct ged *gedp, ged_draw_group_ref ref)
{
    return ged_draw_scene_ref_from_rt_view_ref(
	    ged_draw_registry_group_ref_rt_ref(gedp, ref));
}


static int
_ged_draw_shape_ref_try_obol_path(struct ged *gedp,
				  const char *path,
				  ged_draw_shape_ref_obol_path_cb cb,
				  void *userdata)
{
    if (!path || !path[0] || !cb)
	return 0;
    return (*cb)(gedp, path, userdata);
}


static int
_ged_draw_shape_ref_try_obol_paths(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   bsg_scene_ref shape_ref,
				   ged_draw_shape_ref_obol_path_cb cb,
				   void *userdata,
				   const char *free_label)
{
    const char *semantic_path =
	ged_draw_registry_shape_ref_semantic_path(gedp, ref);
    if (_ged_draw_shape_ref_try_obol_path(gedp, semantic_path, cb,
	    userdata))
	return 1;

    const char *cached_semantic_path =
	ged_draw_registry_shape_ref_cached_semantic_path(gedp, ref);
    if (cached_semantic_path &&
	    (!semantic_path || !BU_STR_EQUAL(cached_semantic_path,
		semantic_path)) &&
	    _ged_draw_shape_ref_try_obol_path(gedp, cached_semantic_path, cb,
		userdata))
	return 1;

    struct ged_draw_obol_context_token *token =
	ged_draw_shape_ref_obol_token(gedp, ref);
    if (token && token->path &&
	    (!semantic_path || !BU_STR_EQUAL(token->path, semantic_path)) &&
	    _ged_draw_shape_ref_try_obol_path(gedp, token->path, cb,
		userdata))
	return 1;

    char *cached_fullpath =
	ged_draw_registry_shape_ref_cached_fullpath_dup(gedp, ref);
    if (cached_fullpath) {
	int matched = ((!semantic_path || !BU_STR_EQUAL(cached_fullpath,
		    semantic_path)) &&
		(!cached_semantic_path || !BU_STR_EQUAL(cached_fullpath,
		    cached_semantic_path)) &&
		(!token || !token->path || !BU_STR_EQUAL(cached_fullpath,
		    token->path)) &&
		_ged_draw_shape_ref_try_obol_path(gedp, cached_fullpath, cb,
		    userdata));
	bu_free(cached_fullpath, free_label ? free_label :
		"GED Obol cached shape path");
	if (matched)
	    return 1;
    }

    struct ged_draw_shape_record_summary shape_summary = {0};
    if (!ged_draw_scene_ref_shape_record_summary(shape_ref, &shape_summary))
	return 0;

    if (shape_summary.fullpath) {
	char *path = db_path_to_string(shape_summary.fullpath);
	if (path) {
	    int matched = _ged_draw_shape_ref_try_obol_path(gedp, path, cb,
		    userdata);
	    bu_free(path, free_label ? free_label : "GED Obol path");
	    if (matched)
		return 1;
	}
    }

    if (shape_summary.display_name &&
	    (!semantic_path || !BU_STR_EQUAL(shape_summary.display_name,
		semantic_path)) &&
	    _ged_draw_shape_ref_try_obol_path(gedp, shape_summary.display_name,
		cb, userdata))
	return 1;

    return _ged_draw_shape_ref_try_obol_path(gedp, shape_summary.leaf_name,
	    cb, userdata);
}


static int
_ged_draw_shape_ref_try_obol_paths_lazy(struct ged *gedp,
					ged_draw_shape_ref ref,
					bsg_scene_ref *shape_ref_out,
					ged_draw_shape_ref_obol_path_cb cb,
					void *userdata,
					const char *free_label)
{
    bsg_scene_ref shape_ref = ged_draw_scene_ref_null();
    int obol_updated = _ged_draw_shape_ref_try_obol_paths(gedp, ref,
	    shape_ref, cb, userdata, free_label);
    if (shape_ref_out)
	*shape_ref_out = shape_ref;
    return obol_updated;
}


struct ged_draw_shape_ref_obol_context_ctx {
    void *context;
};


static int
_ged_draw_shape_ref_obol_context_cb(struct ged *gedp,
				    const char *path,
				    void *userdata)
{
    struct ged_draw_shape_ref_obol_context_ctx *ctx =
	(struct ged_draw_shape_ref_obol_context_ctx *)userdata;
    if (!ctx)
	return 0;

    ctx->context = ged_draw_obol_context_token_with_parent_for_path(gedp,
	    path);
    return ctx->context ? 1 : 0;
}


static void *
_ged_draw_shape_ref_obol_context(struct ged *gedp,
				 ged_draw_shape_ref ref)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_shape_ref_obol_token(gedp, ref);
    if (token)
	return (void *)token;

    struct ged_draw_shape_ref_obol_context_ctx ctx;
    ctx.context = NULL;

    bsg_scene_ref shape_ref = ged_draw_scene_ref_from_rt_view_ref(
	    ged_draw_registry_shape_ref_cache_rt_ref(gedp, ref));
    (void)_ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_shape_ref_obol_context_cb, &ctx,
	    "GED Obol shape-ref context path");
    return ctx.context;
}


static int
_ged_draw_group_ref_try_obol_path(struct ged *gedp,
				  const char *path,
				  ged_draw_group_ref_obol_path_cb cb,
				  void *userdata)
{
    if (!path || !path[0] || !cb)
	return 0;
    return (*cb)(gedp, path, userdata);
}


static int
_ged_draw_group_ref_try_obol_paths(struct ged *gedp,
				   ged_draw_group_ref ref,
				   bsg_scene_ref group_ref,
				   ged_draw_group_ref_obol_path_cb cb,
				   void *userdata,
				   const char *free_label)
{
    const char *semantic_path =
	ged_draw_registry_group_ref_semantic_path(gedp, ref);
    if (_ged_draw_group_ref_try_obol_path(gedp, semantic_path, cb,
	    userdata))
	return 1;

    struct ged_draw_obol_context_token *token =
	ged_draw_group_ref_obol_token(gedp, ref);
    if (token && token->path &&
	    (!semantic_path || !BU_STR_EQUAL(token->path, semantic_path)) &&
	    _ged_draw_group_ref_try_obol_path(gedp, token->path, cb,
		userdata))
	return 1;

    struct ged_draw_scene_tree_summary tree_summary = {0};
    if (ged_draw_scene_ref_tree_summary(group_ref, &tree_summary) &&
	    tree_summary.fullpath) {
	char *path = db_path_to_string(tree_summary.fullpath);
	if (path) {
	    int matched = _ged_draw_group_ref_try_obol_path(gedp, path, cb,
		    userdata);
	    bu_free(path, free_label ? free_label : "GED Obol group path");
	    if (matched)
		return 1;
	}
    }

    struct ged_draw_group_record_summary group_summary = {0};
    if (!ged_draw_scene_ref_group_record_summary(group_ref, &group_summary))
	return 0;
    if (group_summary.path &&
	    (!semantic_path || !BU_STR_EQUAL(group_summary.path,
		semantic_path)) &&
	    _ged_draw_group_ref_try_obol_path(gedp, group_summary.path, cb,
		userdata))
	return 1;

    return 0;
}


static int
_ged_draw_group_ref_try_obol_paths_ensure(
	struct ged *gedp,
	ged_draw_group_ref ref,
	bsg_scene_ref group_ref,
	ged_draw_group_ref_obol_path_cb cb,
	void *userdata,
	const char *free_label,
	const char *preferred_path,
	int mode,
	int overlay)
{
    if (_ged_draw_group_ref_try_obol_paths(gedp, ref, group_ref, cb,
	    userdata, free_label))
	return 1;

    if (ged_draw_scene_ref_is_null(group_ref))
	group_ref = _ged_draw_group_ref_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(group_ref) ||
	    !ged_draw_scene_ref_obol_group_ensure_path(group_ref,
		preferred_path, mode, overlay))
	return 0;

    return _ged_draw_group_ref_try_obol_paths(gedp, ref, group_ref, cb,
	    userdata, free_label);
}


struct ged_draw_group_ref_obol_context_ctx {
    void *context;
};


static int
ged_draw_group_ref_obol_context_cb(struct ged *gedp,
				   const char *path,
				   void *userdata)
{
    struct ged_draw_group_ref_obol_context_ctx *ctx =
	(struct ged_draw_group_ref_obol_context_ctx *)userdata;
    if (!ctx)
	return 0;

    ctx->context = ged_draw_obol_group_context_token_with_parent_for_path(
	    gedp, path);
    return ctx->context ? 1 : 0;
}


static void *
ged_draw_group_ref_obol_context(struct ged *gedp, ged_draw_group_ref ref)
{
    if (!gedp || ged_draw_group_ref_is_null(ref))
	return NULL;

    struct ged_draw_group_ref_obol_context_ctx ctx;
    ctx.context = NULL;

    bsg_scene_ref group_ref = _ged_draw_group_ref_scene_ref(gedp, ref);
    (void)_ged_draw_group_ref_try_obol_paths(gedp, ref, group_ref,
	    ged_draw_group_ref_obol_context_cb, &ctx,
	    "GED Obol group-ref context path");
    return ctx.context;
}


static int
ged_draw_scene_ref_obol_group_path_apply(bsg_scene_ref ref,
					 ged_draw_group_ref_obol_path_cb cb,
					 void *userdata,
					 const char *free_label)
{
    if (ged_draw_scene_ref_is_null(ref) || !cb)
	return 0;

    rt_view_scene_ref scene_ref = ged_draw_scene_ref_to_rt_view_ref(ref);
    struct ged *gedp = ged_draw_scene_rt_ref_owner_gedp(scene_ref);
    if (!gedp)
	return 0;

    const char *semantic_path = ged_draw_scene_rt_ref_semantic_path(scene_ref);
    if (_ged_draw_group_ref_try_obol_path(gedp, semantic_path, cb,
	    userdata))
	return 1;

    ged_draw_group_ref group_token =
	_ged_draw_group_ref_from_scene_ref(gedp, ref);
    const char *registry_path =
	ged_draw_registry_group_ref_semantic_path(gedp, group_token);
    if (registry_path &&
	    (!semantic_path || !BU_STR_EQUAL(registry_path,
		semantic_path)) &&
	    _ged_draw_group_ref_try_obol_path(gedp, registry_path, cb,
		userdata))
	return 1;

    struct ged_draw_scene_tree_summary tree_summary;
    if (ged_draw_scene_ref_tree_summary(ref, &tree_summary)) {
	if (tree_summary.fullpath) {
	    char *path = db_path_to_string(tree_summary.fullpath);
	    if (path) {
		int matched = _ged_draw_group_ref_try_obol_path(gedp, path,
			cb, userdata);
		bu_free(path, free_label ? free_label : "GED Obol group path");
		if (matched)
		    return 1;
	    }
	}
	if (tree_summary.name &&
		(!semantic_path || !BU_STR_EQUAL(tree_summary.name,
		    semantic_path)) &&
		(!registry_path || !BU_STR_EQUAL(tree_summary.name,
		    registry_path)) &&
		_ged_draw_group_ref_try_obol_path(gedp, tree_summary.name,
		    cb, userdata))
	    return 1;
    }

    const char *intent_path = ged_draw_scene_ref_draw_intent_path(ref);
    if (intent_path &&
	    (!semantic_path || !BU_STR_EQUAL(intent_path, semantic_path)) &&
	    (!registry_path || !BU_STR_EQUAL(intent_path, registry_path)) &&
	    _ged_draw_group_ref_try_obol_path(gedp, intent_path, cb,
		userdata))
	return 1;

    return 0;
}


static int
ged_draw_scene_ref_obol_group_ensure_path(
	bsg_scene_ref ref,
	const char *preferred_path,
	int mode,
	int overlay)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 0;

    struct ged *gedp = _ged_draw_scene_ref_owner_gedp(ref);
    if (!gedp)
	return 0;

    const char *ensure_path = preferred_path;
    char *fullpath_path = NULL;

    rt_view_scene_ref scene_ref = ged_draw_scene_ref_to_rt_view_ref(ref);
    const char *semantic_path = ged_draw_scene_rt_ref_semantic_path(scene_ref);
    if ((!ensure_path || !ensure_path[0]) && semantic_path && semantic_path[0])
	ensure_path = semantic_path;

    ged_draw_group_ref group_token = _ged_draw_group_ref_from_scene_ref(gedp,
	    ref);
    const char *registry_path =
	ged_draw_registry_group_ref_semantic_path(gedp, group_token);
    if ((!ensure_path || !ensure_path[0]) && registry_path &&
	    registry_path[0])
	ensure_path = registry_path;

    struct ged_draw_scene_tree_summary tree_summary;
    if (ged_draw_scene_ref_tree_summary(ref, &tree_summary)) {
	if ((!ensure_path || !ensure_path[0]) && tree_summary.fullpath) {
	    fullpath_path = db_path_to_string(tree_summary.fullpath);
	    if (fullpath_path && fullpath_path[0])
		ensure_path = fullpath_path;
	}
	if ((!ensure_path || !ensure_path[0]) && tree_summary.name &&
		tree_summary.name[0])
	    ensure_path = tree_summary.name;
    }

    const char *intent_path = ged_draw_scene_ref_draw_intent_path(ref);
    if ((!ensure_path || !ensure_path[0]) && intent_path && intent_path[0])
	ensure_path = intent_path;

    int ensured = 0;
    if (ensure_path && ensure_path[0])
	ensured = ged_draw_obol_group_ensure_for_path(gedp, ensure_path,
		ensure_path, mode, overlay);

    if (fullpath_path)
	bu_free(fullpath_path, "GED Obol group ensure path");
    return ensured;
}


static int
ged_draw_scene_ref_obol_group_path_apply_ensure(
	bsg_scene_ref ref,
	ged_draw_group_ref_obol_path_cb cb,
	void *userdata,
	const char *free_label,
	const char *preferred_path,
	int mode,
	int overlay)
{
    if (ged_draw_scene_ref_obol_group_path_apply(ref, cb, userdata,
	    free_label))
	return 1;

    if (!ged_draw_scene_ref_obol_group_ensure_path(ref, preferred_path, mode,
	    overlay))
	return 0;

    return ged_draw_scene_ref_obol_group_path_apply(ref, cb, userdata,
	    free_label);
}


struct _ged_draw_obol_display_summary_ctx {
    struct ged_draw_scene_display_summary *out;
};

static int
_ged_draw_obol_display_summary_cb(struct ged *gedp,
				  const char *path,
				  void *userdata)
{
    struct _ged_draw_obol_display_summary_ctx *ctx =
	(struct _ged_draw_obol_display_summary_ctx *)userdata;
    if (!ctx)
	return 0;
    if (ged_draw_obol_database_source_display_summary_for_path(
	    gedp, path, ctx->out))
	return 1;
    return ged_draw_obol_shape_display_summary_for_path(gedp, path,
	    ctx->out);
}


struct _ged_draw_obol_source_summary_ctx {
    struct ged_draw_database_source_summary *out;
};

static int
_ged_draw_obol_source_summary_cb(struct ged *gedp,
				 const char *path,
				 void *userdata)
{
    struct _ged_draw_obol_source_summary_ctx *ctx =
	(struct _ged_draw_obol_source_summary_ctx *)userdata;
    return ctx ? ged_draw_obol_database_source_summary_for_path(
	    gedp, path, ctx->out) : 0;
}


struct _ged_draw_obol_group_record_ctx {
    struct ged_draw_group_record_summary *out;
};

static int
_ged_draw_obol_group_record_cb(struct ged *gedp,
			       const char *path,
			       void *userdata)
{
    struct _ged_draw_obol_group_record_ctx *ctx =
	(struct _ged_draw_obol_group_record_ctx *)userdata;
    return ctx ? ged_draw_obol_group_record_summary_for_path(gedp, path,
	    ctx->out) : 0;
}


static int
ged_draw_scene_ref_obol_group_record_summary(
	bsg_scene_ref ref,
	struct ged_draw_group_record_summary *out)
{
    if (!out)
	return 0;

    struct _ged_draw_obol_group_record_ctx ctx = {out};
    return ged_draw_scene_ref_obol_group_path_apply(ref,
	    _ged_draw_obol_group_record_cb, &ctx,
	    "GED Obol source-adapter group record path");
}


static int
ged_draw_scene_ref_obol_group_remove_cb(struct ged *gedp,
					const char *path,
					void *UNUSED(userdata))
{
    return ged_draw_obol_group_remove_for_path(gedp, path);
}


static int
ged_draw_scene_ref_obol_group_remove(bsg_scene_ref ref)
{
    return ged_draw_scene_ref_obol_group_path_apply(ref,
	    ged_draw_scene_ref_obol_group_remove_cb, NULL,
	    "GED Obol source-adapter group remove path");
}


static int
ged_draw_scene_ref_obol_group_clear_cb(struct ged *gedp,
				       const char *path,
				       void *UNUSED(userdata))
{
    return ged_draw_obol_group_clear_for_path(gedp, path);
}


static int
ged_draw_scene_ref_obol_group_clear(bsg_scene_ref ref)
{
    return ged_draw_scene_ref_obol_group_path_apply(ref,
	    ged_draw_scene_ref_obol_group_clear_cb, NULL,
	    "GED Obol source-adapter group clear path");
}


struct _ged_draw_obol_group_erase_subpath_ctx {
    const char *subpath;
};


static int
ged_draw_scene_ref_obol_group_erase_subpath_cb(struct ged *gedp,
					       const char *path,
					       void *userdata)
{
    struct _ged_draw_obol_group_erase_subpath_ctx *ctx =
	(struct _ged_draw_obol_group_erase_subpath_ctx *)userdata;
    return ctx && ctx->subpath ?
	ged_draw_obol_group_erase_subpath_for_path(gedp, path,
		ctx->subpath) : 0;
}


static int
ged_draw_scene_ref_obol_group_erase_subpath(bsg_scene_ref ref,
					    const char *subpath)
{
    if (!subpath || !subpath[0])
	return 0;

    struct _ged_draw_obol_group_erase_subpath_ctx ctx = {subpath};
    return ged_draw_scene_ref_obol_group_path_apply(ref,
	    ged_draw_scene_ref_obol_group_erase_subpath_cb, &ctx,
	    "GED Obol source-adapter group subpath erase path");
}


struct _ged_draw_obol_group_shape_count_ctx {
    int *out;
};

static int
_ged_draw_obol_group_shape_count_cb(struct ged *gedp,
				    const char *path,
				    void *userdata)
{
    struct _ged_draw_obol_group_shape_count_ctx *ctx =
	(struct _ged_draw_obol_group_shape_count_ctx *)userdata;
    return ctx ? ged_draw_obol_group_shape_count_for_path(gedp, path,
	    ctx->out) : 0;
}


struct _ged_draw_obol_group_child_count_ctx {
    size_t *out;
};

static int
_ged_draw_obol_group_child_count_cb(struct ged *gedp,
				    const char *path,
				    void *userdata)
{
    struct _ged_draw_obol_group_child_count_ctx *ctx =
	(struct _ged_draw_obol_group_child_count_ctx *)userdata;
    return ctx ? ged_draw_obol_group_child_count_for_path(gedp, path,
	    ctx->out) : 0;
}


struct _ged_draw_obol_group_tree_summary_ctx {
    struct ged_draw_scene_tree_summary *out;
    const struct db_full_path *fallback_fullpath;
};

static int
_ged_draw_obol_group_tree_summary_cb(struct ged *gedp,
				     const char *path,
				     void *userdata)
{
    struct _ged_draw_obol_group_tree_summary_ctx *ctx =
	(struct _ged_draw_obol_group_tree_summary_ctx *)userdata;
    return ctx ? ged_draw_obol_group_tree_summary_for_path(gedp, path,
	    ctx->fallback_fullpath, ctx->out) : 0;
}


struct _ged_draw_obol_group_visible_ctx {
    int visible_valid;
    int visible;
};

static int
_ged_draw_obol_group_visible_cb(struct ged *gedp,
				const char *path,
				void *userdata)
{
    struct _ged_draw_obol_group_visible_ctx *ctx =
	(struct _ged_draw_obol_group_visible_ctx *)userdata;
    return ctx ? ged_draw_obol_group_update_display_for_path(gedp, path,
	    ctx->visible_valid, ctx->visible) : 0;
}


struct _ged_draw_obol_group_draw_intent_ctx {
    const char *intent_path;
    int mode_valid;
    int mode;
    int overlay_valid;
    int overlay;
};

static int
_ged_draw_obol_group_draw_intent_cb(struct ged *gedp,
				    const char *path,
				    void *userdata)
{
    struct _ged_draw_obol_group_draw_intent_ctx *ctx =
	(struct _ged_draw_obol_group_draw_intent_ctx *)userdata;
    return ctx ? ged_draw_obol_group_update_draw_intent_for_path(gedp,
	    path, ctx->intent_path, ctx->mode_valid, ctx->mode,
	    ctx->overlay_valid, ctx->overlay) : 0;
}


struct _ged_draw_obol_group_appearance_ctx {
    const struct ged_draw_appearance_settings *settings;
};

static int
_ged_draw_obol_group_update_appearance_cb(struct ged *gedp,
					  const char *path,
					  void *userdata)
{
    struct _ged_draw_obol_group_appearance_ctx *ctx =
	(struct _ged_draw_obol_group_appearance_ctx *)userdata;
    return ctx ? ged_draw_obol_group_update_appearance_for_path(gedp, path,
	    ctx->settings) : 0;
}


struct _ged_draw_obol_group_appearance_read_ctx {
    struct ged_draw_appearance_settings *settings;
};

static int
_ged_draw_obol_group_appearance_read_cb(struct ged *gedp,
					const char *path,
					void *userdata)
{
    struct _ged_draw_obol_group_appearance_read_ctx *ctx =
	(struct _ged_draw_obol_group_appearance_read_ctx *)userdata;
    return ctx ? ged_draw_obol_group_appearance_for_path(gedp, path,
	    ctx->settings) : 0;
}


struct _ged_draw_obol_material_summary_ctx {
    struct ged_draw_shape_material_summary *out;
};

static int
_ged_draw_obol_material_summary_cb(struct ged *gedp,
				   const char *path,
				   void *userdata)
{
    struct _ged_draw_obol_material_summary_ctx *ctx =
	(struct _ged_draw_obol_material_summary_ctx *)userdata;
    return ctx ? ged_draw_obol_database_source_material_summary_for_path(
	    gedp, path, ctx->out) : 0;
}


struct _ged_draw_obol_evaluated_region_ctx {
    int *out;
};

static int
_ged_draw_obol_evaluated_region_cb(struct ged *gedp,
				   const char *path,
				   void *userdata)
{
    struct _ged_draw_obol_evaluated_region_ctx *ctx =
	(struct _ged_draw_obol_evaluated_region_ctx *)userdata;
    return ctx ? ged_draw_obol_database_source_evaluated_region_for_path(
	    gedp, path, ctx->out) : 0;
}


static int
_ged_draw_obol_set_evaluated_region_cb(struct ged *gedp,
				       const char *path,
				       void *userdata)
{
    struct _ged_draw_obol_set_evaluated_region_ctx *ctx =
	(struct _ged_draw_obol_set_evaluated_region_ctx *)userdata;
    return ctx ? ged_draw_obol_database_source_set_evaluated_region_for_path(
	    gedp, path, ctx->evaluated_region) : 0;
}


static int
_ged_draw_obol_set_region_metadata_cb(struct ged *gedp,
				      const char *path,
				      void *userdata)
{
    struct _ged_draw_obol_set_region_metadata_ctx *ctx =
	(struct _ged_draw_obol_set_region_metadata_ctx *)userdata;
    return ctx ? ged_draw_obol_database_source_set_region_metadata_for_path(
	    gedp, path, ctx->region_id, ctx->aircode, ctx->los,
	    ctx->material_id) : 0;
}


struct _ged_draw_obol_line_point_ctx {
    size_t index;
    fastf_t *out;
};

static int
_ged_draw_obol_last_point_cb(struct ged *gedp,
			     const char *path,
			     void *userdata)
{
    fastf_t *out = (fastf_t *)userdata;
    return out ? ged_draw_obol_database_source_last_point_for_path(
	    gedp, path, out) : 0;
}


struct _ged_draw_obol_line_summary_ctx {
    struct ged_draw_view_line_summary *out;
};

static int
_ged_draw_obol_line_summary_cb(struct ged *gedp,
			       const char *path,
			       void *userdata)
{
    struct _ged_draw_obol_line_summary_ctx *ctx =
	(struct _ged_draw_obol_line_summary_ctx *)userdata;
    if (!ctx)
	return 0;
    if (ged_draw_obol_database_source_line_summary_for_path(gedp, path,
	    ctx->out))
	return 1;
    return ged_draw_obol_shape_line_summary_for_path(gedp, path,
	    ctx->out);
}


static int
_ged_draw_obol_line_point_cb(struct ged *gedp,
			     const char *path,
			     void *userdata)
{
    struct _ged_draw_obol_line_point_ctx *ctx =
	(struct _ged_draw_obol_line_point_ctx *)userdata;
    if (!ctx)
	return 0;
    if (ged_draw_obol_database_source_line_point_at_for_path(gedp, path,
	    ctx->index, ctx->out))
	return 1;
    return ged_draw_obol_shape_line_point_at_for_path(gedp, path,
	    ctx->index, ctx->out);
}


struct _ged_draw_obol_line_command_ctx {
    size_t index;
    int *out;
};

static int
_ged_draw_obol_line_command_cb(struct ged *gedp,
			       const char *path,
			       void *userdata)
{
    struct _ged_draw_obol_line_command_ctx *ctx =
	(struct _ged_draw_obol_line_command_ctx *)userdata;
    if (!ctx)
	return 0;
    if (ged_draw_obol_database_source_line_command_at_for_path(gedp, path,
	    ctx->index, ctx->out))
	return 1;
    return ged_draw_obol_shape_line_command_at_for_path(gedp, path,
	    ctx->index, ctx->out);
}


struct _ged_draw_obol_geometry_summary_ctx {
    struct ged_draw_shape_geometry_summary *out;
};

static int
_ged_draw_obol_geometry_summary_cb(struct ged *gedp,
				   const char *path,
				   void *userdata)
{
    struct _ged_draw_obol_geometry_summary_ctx *ctx =
	(struct _ged_draw_obol_geometry_summary_ctx *)userdata;
    if (!ctx)
	return 0;
    if (ged_draw_obol_database_source_geometry_summary_for_path(gedp, path,
	    ctx->out))
	return 1;
    return ged_draw_obol_shape_geometry_summary_for_path(gedp, path,
	    ctx->out);
}


static int
_ged_draw_obol_clear_vlist_cb(struct ged *gedp,
			      const char *path,
			      void *UNUSED(userdata))
{
    return ged_draw_obol_database_source_clear_vlist_for_path(gedp, path);
}


static int
_ged_draw_obol_clear_mesh_cb(struct ged *gedp,
			     const char *path,
			     void *UNUSED(userdata))
{
    return ged_draw_obol_database_source_clear_mesh_for_path(gedp, path);
}


static int
_ged_draw_obol_clear_auxiliary_shapes_cb(struct ged *gedp,
					 const char *path,
					 void *UNUSED(userdata))
{
    return ged_draw_obol_database_source_clear_auxiliary_shapes_for_path(gedp,
	    path);
}


struct _ged_draw_obol_publish_line_ctx {
    const point_t *points;
    const int *commands;
    size_t point_count;
};

static int
_ged_draw_obol_publish_line_cb(struct ged *gedp,
			       const char *path,
			       void *userdata)
{
    struct _ged_draw_obol_publish_line_ctx *ctx =
	(struct _ged_draw_obol_publish_line_ctx *)userdata;
    return ctx ? ged_draw_obol_database_source_publish_line_set_for_path(
	    gedp, path, ctx->points, ctx->commands, ctx->point_count) : 0;
}


static int
_ged_draw_obol_publish_annotation_line_cb(struct ged *gedp,
					  const char *path,
					  void *userdata)
{
    struct _ged_draw_obol_publish_line_ctx *ctx =
	(struct _ged_draw_obol_publish_line_ctx *)userdata;
    return ctx ?
	ged_draw_obol_database_source_publish_annotation_line_set_for_path(
		gedp, path, ctx->points, ctx->commands,
		ctx->point_count) : 0;
}


struct _ged_draw_obol_publish_annotation_record_ctx {
    const point_t base_point;
    const point_t *annotation_points;
    size_t annotation_point_count;
    const struct ged_draw_obol_annotation_segment *segments;
    size_t segment_count;
    const point_t *line_points;
    const int *line_commands;
    size_t line_point_count;
};

static int
_ged_draw_obol_publish_annotation_record_cb(struct ged *gedp,
					    const char *path,
					    void *userdata)
{
    struct _ged_draw_obol_publish_annotation_record_ctx *ctx =
	(struct _ged_draw_obol_publish_annotation_record_ctx *)userdata;
    return ctx ?
	ged_draw_obol_database_source_publish_annotation_record_for_path(
		gedp, path, ctx->base_point, ctx->annotation_points,
		ctx->annotation_point_count, ctx->segments,
		ctx->segment_count, ctx->line_points, ctx->line_commands,
		ctx->line_point_count) : 0;
}


struct _ged_draw_obol_publish_auxiliary_line_ctx {
    const char *name;
    const point_t *points;
    const int *commands;
    size_t point_count;
    const struct ged_draw_scene_display_summary *display_state;
};

static int
_ged_draw_obol_publish_auxiliary_line_cb(struct ged *gedp,
					 const char *path,
					 void *userdata)
{
    struct _ged_draw_obol_publish_auxiliary_line_ctx *ctx =
	(struct _ged_draw_obol_publish_auxiliary_line_ctx *)userdata;
    return ctx ? ged_draw_obol_database_source_publish_auxiliary_line_set_for_path(
	    gedp, path, ctx->name, ctx->points, ctx->commands,
	    ctx->point_count, ctx->display_state) : 0;
}


struct _ged_draw_obol_placement_ctx {
    bsg_scene_ref ref;
};

static int
_ged_draw_obol_sync_placement_cb(struct ged *gedp,
				 const char *path,
				 void *userdata)
{
    struct _ged_draw_obol_placement_ctx *ctx =
	(struct _ged_draw_obol_placement_ctx *)userdata;
    if (!ctx || ged_draw_scene_ref_is_null(ctx->ref))
	return 0;

    mat_t draw_mat;
    point_t draw_center;
    MAT_IDN(draw_mat);
    VSET(draw_center, 0.0, 0.0, 0.0);

    const int draw_mat_valid =
	ged_draw_scene_ref_draw_mat(ctx->ref, draw_mat) ? 1 : 0;
    const fastf_t draw_size = ged_draw_scene_ref_draw_size(ctx->ref);
    const int bounds_valid = draw_size > 0.0 ? 1 : 0;
    if (bounds_valid)
	ged_draw_scene_ref_draw_center(ctx->ref, draw_center);

    return ged_draw_obol_database_source_set_placement_for_path(gedp,
	    path, draw_mat_valid, draw_mat, bounds_valid, draw_center,
	    bounds_valid, draw_size);
}


static int
ged_draw_shape_ref_obol_sync_source_placement(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	bsg_scene_ref shape_ref)
{
    struct _ged_draw_obol_placement_ctx ctx = {shape_ref};
    return _ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_sync_placement_cb, &ctx,
	    "GED Obol source placement path");
}


static int
ged_draw_scene_ref_obol_sync_source_placement(bsg_scene_ref ref)
{
    struct _ged_draw_obol_placement_ctx ctx = {ref};
    return ged_draw_scene_ref_obol_source_path_apply(ref,
	    _ged_draw_obol_sync_placement_cb, &ctx);
}


struct _ged_draw_obol_center_ctx {
    point_t center;
};

static int
_ged_draw_obol_set_center_cb(struct ged *gedp,
			     const char *path,
			     void *userdata)
{
    struct _ged_draw_obol_center_ctx *ctx =
	(struct _ged_draw_obol_center_ctx *)userdata;
    return ctx ? ged_draw_obol_database_source_set_placement_for_path(gedp,
	    path, 0, NULL, 1, ctx->center, 0, 0.0) : 0;
}


static int
_ged_draw_obol_update_bounds_cb(struct ged *gedp,
				const char *path,
				void *UNUSED(userdata))
{
    return ged_draw_obol_database_source_update_vlist_bounds_for_path(gedp,
	    path);
}


static int
_ged_draw_obol_update_display_cb(struct ged *gedp,
				 const char *path,
				 void *userdata)
{
    struct _ged_draw_obol_update_display_ctx *ctx =
	(struct _ged_draw_obol_update_display_ctx *)userdata;
    if (!ctx)
	return 0;
    if (ged_draw_obol_database_source_update_display_for_path(
	    gedp, path,
	    ctx->visible_valid, ctx->visible,
	    ctx->highlighted_valid, ctx->highlighted,
	    ctx->draw_mode_valid, ctx->draw_mode,
	    ctx->line_style_valid, ctx->line_style,
	    ctx->line_width_valid, ctx->line_width,
	    ctx->transparency_valid, ctx->transparency,
	    ctx->color_valid, ctx->color,
	    ctx->material_color_valid, ctx->material_color,
	    ctx->material_revision_valid, ctx->material_revision))
	return 1;
    return ged_draw_obol_shape_update_display_for_path(gedp, path,
	    ctx->visible_valid, ctx->visible,
	    ctx->highlighted_valid, ctx->highlighted,
	    ctx->draw_mode_valid, ctx->draw_mode,
	    ctx->line_style_valid, ctx->line_style,
	    ctx->line_width_valid, ctx->line_width,
	    ctx->transparency_valid, ctx->transparency,
	    ctx->color_valid, ctx->color,
	    ctx->material_color_valid, ctx->material_color,
	    ctx->material_revision_valid, ctx->material_revision);
}


static int
ged_draw_scene_ref_obol_update_display(
	bsg_scene_ref ref,
	int visible_valid,
	int visible,
	int highlighted_valid,
	int highlighted,
	int draw_mode_valid,
	int draw_mode,
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
	uint64_t material_revision)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 0;

    /* Group-only updates still use their dedicated public paths elsewhere. */
    if (ged_draw_scene_ref_is_group(ref) &&
	    !ged_draw_scene_ref_is_database_source(ref))
	return 0;

    struct _ged_draw_obol_update_display_ctx ctx = {
	visible_valid,
	visible,
	highlighted_valid,
	highlighted,
	draw_mode_valid,
	draw_mode,
	line_style_valid,
	line_style,
	line_width_valid,
	line_width,
	transparency_valid,
	transparency,
	color_valid,
	color,
	material_color_valid,
	material_color,
	material_revision_valid,
	material_revision
    };
    return ged_draw_scene_ref_obol_source_path_apply_ensure(ref,
	    (ged_draw_scene_ref_obol_source_path_cb)_ged_draw_obol_update_display_cb,
	    &ctx);
}


struct _ged_draw_obol_placement_patch_ctx {
    int draw_mat_valid;
    const fastf_t *draw_mat;
    int draw_center_valid;
    const fastf_t *draw_center;
    int draw_size_valid;
    fastf_t draw_size;
};


static int
_ged_draw_obol_set_placement_cb(struct ged *gedp,
				const char *path,
				void *userdata)
{
    struct _ged_draw_obol_placement_patch_ctx *ctx =
	(struct _ged_draw_obol_placement_patch_ctx *)userdata;
    return ctx ? ged_draw_obol_database_source_set_placement_for_path(gedp,
	    path, ctx->draw_mat_valid, ctx->draw_mat,
	    ctx->draw_center_valid, ctx->draw_center,
	    ctx->draw_size_valid, ctx->draw_size) : 0;
}


static int
ged_draw_scene_ref_obol_set_placement(bsg_scene_ref ref,
				      int draw_mat_valid,
				      const mat_t draw_mat,
				      int draw_center_valid,
				      const point_t draw_center,
				      int draw_size_valid,
				      fastf_t draw_size)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 0;

    if (ged_draw_scene_ref_is_group(ref) &&
	    !ged_draw_scene_ref_is_database_source(ref))
	return 0;

    struct _ged_draw_obol_placement_patch_ctx ctx = {
	draw_mat_valid,
	draw_mat,
	draw_center_valid,
	draw_center,
	draw_size_valid,
	draw_size
    };
    return ged_draw_scene_ref_obol_source_path_apply_ensure(ref,
	    _ged_draw_obol_set_placement_cb, &ctx);
}


struct _ged_draw_obol_bounds_ctx {
    point_t bmin;
    point_t bmax;
};


static int
_ged_draw_obol_set_bounds_cb(struct ged *gedp,
			     const char *path,
			     void *userdata)
{
    struct _ged_draw_obol_bounds_ctx *ctx =
	(struct _ged_draw_obol_bounds_ctx *)userdata;
    return ctx ? ged_draw_obol_database_source_set_bounds_for_path(gedp,
	    path, ctx->bmin, ctx->bmax) : 0;
}


static int
ged_draw_scene_ref_obol_set_bounds(bsg_scene_ref ref,
				   const point_t bmin,
				   const point_t bmax)
{
    if (ged_draw_scene_ref_is_null(ref) || !bmin || !bmax)
	return 0;

    if (ged_draw_scene_ref_is_group(ref) &&
	    !ged_draw_scene_ref_is_database_source(ref))
	return 0;

    struct _ged_draw_obol_bounds_ctx ctx;
    VMOVE(ctx.bmin, bmin);
    VMOVE(ctx.bmax, bmax);
    return ged_draw_scene_ref_obol_source_path_apply_ensure(ref,
	    _ged_draw_obol_set_bounds_cb, &ctx);
}


struct _ged_draw_obol_stale_ctx {
    ged_draw_stale_reason reason;
};

static int
_ged_draw_obol_mark_stale_cb(struct ged *gedp,
			     const char *path,
			     void *userdata)
{
    struct _ged_draw_obol_stale_ctx *ctx =
	(struct _ged_draw_obol_stale_ctx *)userdata;
    return ctx ? ged_draw_obol_database_source_mark_stale_for_path(gedp,
	    path, ctx->reason) : 0;
}


static int
_ged_draw_obol_remove_source_cb(struct ged *gedp,
				const char *path,
				void *UNUSED(userdata))
{
    return ged_draw_obol_database_source_remove_for_path(gedp, path);
}


struct _ged_draw_obol_translate_ctx {
    vect_t xlate;
};

static int
_ged_draw_obol_translate_cb(struct ged *gedp,
			    const char *path,
			    void *userdata)
{
    struct _ged_draw_obol_translate_ctx *ctx =
	(struct _ged_draw_obol_translate_ctx *)userdata;
    return ctx ? ged_draw_obol_database_source_translate_vlist_for_path(
	    gedp, path, ctx->xlate) : 0;
}


static bsg_scene_ref
_ged_draw_source_root_scene_ref(struct ged *gedp)
{
    return _ged_draw_group_ref_scene_ref(gedp, ged_scene_root_group_ref(gedp));
}


static ged_draw_group_ref
_ged_draw_active_scope_group_ref(struct ged *gedp,
				 void *view_ctx,
				 int create,
				 int fallback_root)
{
    ged_draw_group_ref null_ref = GED_DRAW_GROUP_REF_NULL;
    if (!gedp)
	return null_ref;

    if (view_ctx && rt_view_context_is_independent(view_ctx)) {
	ged_draw_group_ref root_ref = null_ref;
	if (create || fallback_root) {
	    if (create && !ged_draw_ensure_root_attached(gedp))
		return null_ref;
	    root_ref = ged_scene_root_group_ref(gedp);
	}

	rt_view_scene_ref scope_rt_ref =
	    rt_view_context_independent_scope_ref(view_ctx, create);
	bsg_scene_ref scope_ref = ged_draw_scene_ref_from_rt_view_ref(scope_rt_ref);
	if (!ged_draw_scene_ref_is_null(scope_ref)) {
	    ged_draw_group_ref group_ref =
		_ged_draw_group_ref_from_scene_ref(gedp, scope_ref);
	    if (!ged_draw_group_ref_is_null(group_ref)) {
		group_ref.scene_revision = 0;
		return group_ref;
	    }
	}
	return fallback_root ? root_ref : null_ref;
    }

    if (create && !ged_draw_ensure_root_attached(gedp))
	return null_ref;
    return ged_scene_root_group_ref(gedp);
}


static bsg_scene_ref
_ged_draw_active_scope_scene_ref(struct ged *gedp,
				 void *view_ctx,
				 int create,
				 int fallback_root)
{
    return _ged_draw_group_ref_scene_ref(gedp,
	    _ged_draw_active_scope_group_ref(gedp, view_ctx, create,
		fallback_root));
}


static int
ged_draw_scene_ref_is_group(bsg_scene_ref ref)
{
    return !bsg_scene_ref_is_null(ref) &&
	bsg_scene_ref_type(ref) == BSG_SCENE_ELEMENT_GROUP;
}


static int
ged_draw_scene_ref_is_shape(bsg_scene_ref ref)
{
    return !bsg_scene_ref_is_null(ref) &&
	bsg_scene_ref_type(ref) == BSG_SCENE_ELEMENT_SHAPE;
}


static bsg_scene_ref
ged_draw_scene_ref_parent(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? bsg_scene_ref_null() :
	bsg_scene_parent(ref);
}


static bsg_scene_ref
ged_draw_scene_ref_source_owner(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return bsg_scene_ref_null();
    if (bsg_database_source_ref_is_container(
		bsg_database_source_ref_from_scene(ref)))
	return ref;

    ged_draw_shape_state *shape_data = _ged_draw_shape_state_get_scene_ref(ref);
    if (shape_data && !rt_view_scene_ref_is_null(shape_data->source_ref)) {
	bsg_scene_ref source_ref =
	    ged_draw_scene_ref_from_rt_view_ref(shape_data->source_ref);
	if (!ged_draw_scene_ref_is_null(source_ref))
	    return source_ref;
    }

    bsg_scene_ref parent_ref = ged_draw_scene_ref_parent(ref);
    while (!ged_draw_scene_ref_is_null(parent_ref)) {
	if (bsg_database_source_ref_is_container(
		    bsg_database_source_ref_from_scene(parent_ref)))
	    return parent_ref;
	parent_ref = ged_draw_scene_ref_parent(parent_ref);
    }

    return ref;
}


static bsg_scene_ref
ged_draw_view_context_database_source_create(void *view_ctx, const char *name)
{
    struct bsg_view *v = (struct bsg_view *)view_ctx;
    bsg_database_source_ref source =
	bsg_database_source_ref_create(v, name ? name : "_db_source");
    return bsg_database_source_ref_as_scene(source);
}


static int
ged_draw_scene_ref_is_database_source(bsg_scene_ref ref)
{
    bsg_database_source_ref source =
	bsg_database_source_ref_from_scene(ref);
    if (bsg_database_source_ref_is_null(source))
	return 0;
    return bsg_scene_is_database_source(ref) ||
	bsg_database_source_ref_is_container(source);
}


static bsg_database_source_stale_reason
ged_draw_source_stale_reason_to_bsg(ged_draw_stale_reason reason)
{
    switch (reason) {
	case GED_DRAW_STALE_SOURCE_CHANGED:
	    return BSG_DATABASE_SOURCE_STALE_SOURCE_CHANGED;
	case GED_DRAW_STALE_VIEW_INPUT_CHANGED:
	    return BSG_DATABASE_SOURCE_STALE_VIEW_INPUT_CHANGED;
	case GED_DRAW_STALE_SETTINGS_CHANGED:
	    return BSG_DATABASE_SOURCE_STALE_SETTINGS_CHANGED;
	case GED_DRAW_STALE_FORCED:
	    return BSG_DATABASE_SOURCE_STALE_FORCED;
	case GED_DRAW_STALE_UPDATE_FAILED:
	    return BSG_DATABASE_SOURCE_STALE_UPDATE_FAILED;
	case GED_DRAW_STALE_NONE:
	default:
	    return BSG_DATABASE_SOURCE_STALE_NONE;
    }
}


static ged_draw_stale_reason
ged_draw_source_stale_reason_from_bsg(bsg_database_source_stale_reason reason)
{
    switch (reason) {
	case BSG_DATABASE_SOURCE_STALE_SOURCE_CHANGED:
	    return GED_DRAW_STALE_SOURCE_CHANGED;
	case BSG_DATABASE_SOURCE_STALE_VIEW_INPUT_CHANGED:
	    return GED_DRAW_STALE_VIEW_INPUT_CHANGED;
	case BSG_DATABASE_SOURCE_STALE_SETTINGS_CHANGED:
	    return GED_DRAW_STALE_SETTINGS_CHANGED;
	case BSG_DATABASE_SOURCE_STALE_FORCED:
	    return GED_DRAW_STALE_FORCED;
	case BSG_DATABASE_SOURCE_STALE_UPDATE_FAILED:
	    return GED_DRAW_STALE_UPDATE_FAILED;
	case BSG_DATABASE_SOURCE_STALE_NONE:
	default:
	    return GED_DRAW_STALE_NONE;
    }
}


static bsg_database_source_material_policy
ged_draw_source_material_policy_to_bsg(
	ged_draw_database_source_material_policy policy)
{
    switch (policy) {
	case GED_DRAW_DATABASE_SOURCE_MATERIAL_DATABASE:
	    return BSG_DATABASE_SOURCE_MATERIAL_DATABASE;
	case GED_DRAW_DATABASE_SOURCE_MATERIAL_INHERIT:
	default:
	    return BSG_DATABASE_SOURCE_MATERIAL_INHERIT;
    }
}


static ged_draw_database_source_material_policy
ged_draw_source_material_policy_from_bsg(
	bsg_database_source_material_policy policy)
{
    switch (policy) {
	case BSG_DATABASE_SOURCE_MATERIAL_DATABASE:
	    return GED_DRAW_DATABASE_SOURCE_MATERIAL_DATABASE;
	case BSG_DATABASE_SOURCE_MATERIAL_INHERIT:
	default:
	    return GED_DRAW_DATABASE_SOURCE_MATERIAL_INHERIT;
    }
}


static bsg_database_source_realization_status
ged_draw_source_realization_status_to_bsg(
	ged_draw_database_source_realization_status status)
{
    switch (status) {
	case GED_DRAW_DATABASE_SOURCE_REALIZATION_CURRENT:
	    return BSG_DATABASE_SOURCE_REALIZATION_CURRENT;
	case GED_DRAW_DATABASE_SOURCE_REALIZATION_STALE:
	default:
	    return BSG_DATABASE_SOURCE_REALIZATION_STALE;
    }
}


static ged_draw_database_source_realization_status
ged_draw_source_realization_status_from_bsg(
	bsg_database_source_realization_status status)
{
    switch (status) {
	case BSG_DATABASE_SOURCE_REALIZATION_CURRENT:
	    return GED_DRAW_DATABASE_SOURCE_REALIZATION_CURRENT;
	case BSG_DATABASE_SOURCE_REALIZATION_STALE:
	default:
	    return GED_DRAW_DATABASE_SOURCE_REALIZATION_STALE;
    }
}


static int
ged_draw_source_realization_roles_to_bsg(int roles)
{
    int out = BSG_DATABASE_SOURCE_REALIZATION_ROLE_NONE;
    if (roles & GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_CSG)
	out |= BSG_DATABASE_SOURCE_REALIZATION_ROLE_CSG;
    if (roles & GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_MESH)
	out |= BSG_DATABASE_SOURCE_REALIZATION_ROLE_MESH;
    return out;
}


static int
ged_draw_source_realization_roles_from_bsg(int roles)
{
    int out = GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_NONE;
    if (roles & BSG_DATABASE_SOURCE_REALIZATION_ROLE_CSG)
	out |= GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_CSG;
    if (roles & BSG_DATABASE_SOURCE_REALIZATION_ROLE_MESH)
	out |= GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_MESH;
    return out;
}


static void
ged_draw_database_source_record_from_obol(
	struct ged_draw_database_source_record *record,
	const struct ged_draw_obol_database_source_record *obol_record)
{
    if (!record || !obol_record || !obol_record->valid)
	return;

    memset(record, 0, sizeof(*record));
    record->database_path = obol_record->database_path;
    record->draw_mode = obol_record->draw_mode;
    record->material_policy =
	obol_record->material_policy ==
	GED_DRAW_OBOL_DATABASE_SOURCE_MATERIAL_DATABASE ?
	GED_DRAW_DATABASE_SOURCE_MATERIAL_DATABASE :
	GED_DRAW_DATABASE_SOURCE_MATERIAL_INHERIT;
    record->source_revision = obol_record->source_revision;
    record->inputs_revision = obol_record->inputs_revision;
    record->realized_source_revision =
	obol_record->realized_source_revision;
    record->realized_inputs_revision =
	obol_record->realized_inputs_revision;
    record->stale_reason = obol_record->stale_reason;
    record->realization_identity = obol_record->realization_identity;
    record->realization_status =
	obol_record->realization_status ==
	GED_DRAW_OBOL_DATABASE_SOURCE_REALIZATION_CURRENT ?
	GED_DRAW_DATABASE_SOURCE_REALIZATION_CURRENT :
	GED_DRAW_DATABASE_SOURCE_REALIZATION_STALE;
    record->realization_role_flags =
	obol_record->realization_role_flags;
    record->realization_view_dependent =
	obol_record->realization_view_dependent;
    record->realization_view_scale =
	obol_record->realization_view_scale;
    record->realization_bot_threshold =
	obol_record->realization_bot_threshold;
    record->realization_curve_scale =
	obol_record->realization_curve_scale;
    record->realization_point_scale =
	obol_record->realization_point_scale;
}


static void
ged_draw_database_source_record_to_obol(
	struct ged_draw_obol_database_source_record *obol_record,
	const struct ged_draw_database_source_record *record)
{
    if (!obol_record || !record)
	return;

    memset(obol_record, 0, sizeof(*obol_record));
    obol_record->valid = 1;
    obol_record->database_path = record->database_path;
    obol_record->draw_mode = record->draw_mode;
    obol_record->material_policy =
	record->material_policy ==
	GED_DRAW_DATABASE_SOURCE_MATERIAL_DATABASE ?
	GED_DRAW_OBOL_DATABASE_SOURCE_MATERIAL_DATABASE :
	GED_DRAW_OBOL_DATABASE_SOURCE_MATERIAL_INHERIT;
    obol_record->source_revision = record->source_revision;
    obol_record->inputs_revision = record->inputs_revision;
    obol_record->realized_source_revision =
	record->realized_source_revision;
    obol_record->realized_inputs_revision =
	record->realized_inputs_revision;
    obol_record->stale_reason = record->stale_reason;
    obol_record->realization_identity = record->realization_identity;
    obol_record->realization_status =
	record->realization_status ==
	GED_DRAW_DATABASE_SOURCE_REALIZATION_CURRENT ?
	GED_DRAW_OBOL_DATABASE_SOURCE_REALIZATION_CURRENT :
	GED_DRAW_OBOL_DATABASE_SOURCE_REALIZATION_STALE;
    obol_record->realization_role_flags =
	record->realization_role_flags;
    obol_record->realization_view_dependent =
	record->realization_view_dependent;
    obol_record->realization_view_scale =
	record->realization_view_scale;
    obol_record->realization_bot_threshold =
	record->realization_bot_threshold;
    obol_record->realization_curve_scale =
	record->realization_curve_scale;
    obol_record->realization_point_scale =
	record->realization_point_scale;
}


static int
ged_draw_scene_ref_obol_database_source_record_get_cb(
	struct ged *gedp,
	const char *path,
	void *data)
{
    return ged_draw_obol_database_source_record_for_path(gedp, path,
	    (struct ged_draw_obol_database_source_record *)data);
}


static int
ged_draw_scene_ref_obol_database_source_record_get(
	bsg_scene_ref ref,
	struct ged_draw_database_source_record *record)
{
    struct ged_draw_obol_database_source_record obol_record;
    memset(&obol_record, 0, sizeof(obol_record));
    if (!ged_draw_scene_ref_obol_source_path_apply(ref,
	    ged_draw_scene_ref_obol_database_source_record_get_cb,
	    &obol_record))
	return 0;

    ged_draw_database_source_record_from_obol(record, &obol_record);
    return obol_record.valid ? 1 : 0;
}


static int
ged_draw_scene_ref_obol_database_source_record_apply_cb(
	struct ged *gedp,
	const char *path,
	void *data)
{
    return ged_draw_obol_database_source_apply_record_for_path(gedp, path,
	    (const struct ged_draw_obol_database_source_record *)data);
}


static int
ged_draw_scene_ref_obol_database_source_record_apply(
	bsg_scene_ref ref,
	const struct ged_draw_database_source_record *record)
{
    struct ged_draw_obol_database_source_record obol_record;
    ged_draw_database_source_record_to_obol(&obol_record, record);
    return ged_draw_scene_ref_obol_source_path_apply(ref,
	    ged_draw_scene_ref_obol_database_source_record_apply_cb,
	    &obol_record);
}


static int
ged_draw_scene_ref_obol_database_source_runtime_cb(
	struct ged *gedp,
	const char *path,
	void *data)
{
    return ged_draw_obol_database_source_runtime_for_path(gedp, path,
	    (struct ged_draw_obol_database_source_runtime *)data);
}


static int
ged_draw_scene_ref_obol_database_source_runtime(
	bsg_scene_ref ref,
	struct ged_draw_obol_database_source_runtime *out)
{
    if (!out)
	return 0;
    memset(out, 0, sizeof(*out));
    return ged_draw_scene_ref_obol_source_path_apply(ref,
	    ged_draw_scene_ref_obol_database_source_runtime_cb, out);
}


static int
ged_draw_scene_ref_obol_draw_state_cb(
	struct ged *gedp,
	const char *path,
	void *data)
{
    return ged_draw_obol_database_source_draw_state_for_path(gedp, path,
	    (struct ged_draw_obol_draw_state_summary *)data);
}


static int
ged_draw_scene_ref_obol_draw_state_summary(
	bsg_scene_ref ref,
	struct ged_draw_obol_draw_state_summary *out)
{
    if (!out)
	return 0;
    memset(out, 0, sizeof(*out));
    return ged_draw_scene_ref_obol_source_path_apply(ref,
	    ged_draw_scene_ref_obol_draw_state_cb, out);
}


static int
ged_draw_scene_ref_database_source_record_get(
	bsg_scene_ref ref,
	struct ged_draw_database_source_record *record)
{
    if (!record)
	return 0;

    memset(record, 0, sizeof(*record));
    if (ged_draw_scene_ref_obol_database_source_record_get(ref, record))
	return 1;

    bsg_database_source_ref source =
	bsg_database_source_ref_from_scene(ref);
    if (bsg_database_source_ref_is_null(source))
	return 0;

    struct bsg_database_source_record bsg_record;
    memset(&bsg_record, 0, sizeof(bsg_record));
    if (!bsg_database_source_record_get(source, &bsg_record))
	return 0;

    record->database_path = bsg_record.database_path;
    record->draw_mode = (int)bsg_record.draw_mode;
    record->material_policy =
	ged_draw_source_material_policy_from_bsg(
		bsg_record.material_policy);
    record->source_revision = bsg_record.source_revision;
    record->inputs_revision = bsg_record.inputs_revision;
    record->realized_source_revision =
	bsg_record.realized_source_revision;
    record->realized_inputs_revision =
	bsg_record.realized_inputs_revision;
    record->stale_reason =
	ged_draw_source_stale_reason_from_bsg(bsg_record.stale_reason);
    record->realization_identity = bsg_record.realization_identity;
    record->realization_status =
	ged_draw_source_realization_status_from_bsg(
		bsg_record.realization_status);
    record->realization_role_flags =
	ged_draw_source_realization_roles_from_bsg(
		bsg_record.realization_role_flags);
    record->realization_view_dependent =
	bsg_record.realization_view_dependent;
    record->realization_view_scale =
	(fastf_t)bsg_record.realization_view_scale;
    record->realization_bot_threshold =
	bsg_record.realization_bot_threshold;
    record->realization_curve_scale =
	(fastf_t)bsg_record.realization_curve_scale;
    record->realization_point_scale =
	(fastf_t)bsg_record.realization_point_scale;
    return 1;
}


static int
ged_draw_scene_ref_database_source_record_apply(
	bsg_scene_ref ref,
	const struct ged_draw_database_source_record *record)
{
    if (!record)
	return 0;

    int obol_updated =
	ged_draw_scene_ref_obol_database_source_record_apply(ref, record);
    if (obol_updated)
	return 1;

    bsg_database_source_ref source =
	bsg_database_source_ref_from_scene(ref);
    if (bsg_database_source_ref_is_null(source))
	return 0;

    struct bsg_database_source_record bsg_record;
    memset(&bsg_record, 0, sizeof(bsg_record));
    bsg_record.database_path = record->database_path;
    bsg_record.draw_mode = (bsg_draw_mode)record->draw_mode;
    bsg_record.material_policy =
	ged_draw_source_material_policy_to_bsg(record->material_policy);
    bsg_record.source_revision = record->source_revision;
    bsg_record.inputs_revision = record->inputs_revision;
    bsg_record.realized_source_revision =
	record->realized_source_revision;
    bsg_record.realized_inputs_revision =
	record->realized_inputs_revision;
    bsg_record.stale_reason =
	ged_draw_source_stale_reason_to_bsg(record->stale_reason);
    bsg_record.realization_identity = record->realization_identity;
    bsg_record.realization_status =
	ged_draw_source_realization_status_to_bsg(
		record->realization_status);
    bsg_record.realization_role_flags =
	ged_draw_source_realization_roles_to_bsg(
		record->realization_role_flags);
    bsg_record.realization_view_dependent =
	record->realization_view_dependent;
    bsg_record.realization_view_scale =
	(double)record->realization_view_scale;
    bsg_record.realization_bot_threshold =
	record->realization_bot_threshold;
    bsg_record.realization_curve_scale =
	(double)record->realization_curve_scale;
    bsg_record.realization_point_scale =
	(double)record->realization_point_scale;

    /* Obol owns database-source record mutation; retained BSG record apply is
     * no-owned-scene legacy adapter compatibility only. */
    return bsg_database_source_record_apply(source, &bsg_record);
}


static bsg_scene_ref
ged_draw_scene_ref_shape_owning_group(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return ged_draw_scene_ref_null();

    bsg_scene_ref group_ref = ged_draw_scene_ref_parent(ref);
    while (!ged_draw_scene_ref_is_null(group_ref)) {
	bsg_scene_ref parent_ref = ged_draw_scene_ref_parent(group_ref);
	if (ged_draw_scene_ref_is_null(parent_ref))
	    return group_ref;
	if (ged_draw_scene_ref_is_view_scope(parent_ref))
	    return group_ref;
	bsg_scene_ref grandparent_ref = ged_draw_scene_ref_parent(parent_ref);
	if (ged_draw_scene_ref_is_null(grandparent_ref))
	    return group_ref;
	group_ref = parent_ref;
    }
    return ged_draw_scene_ref_null();
}


static size_t
ged_draw_scene_ref_child_count(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? 0 : bsg_scene_child_count(ref);
}


static bsg_scene_ref
ged_draw_scene_ref_child_at(bsg_scene_ref ref, size_t idx)
{
    return bsg_scene_ref_is_null(ref) ? bsg_scene_ref_null() :
	bsg_scene_child_at(ref, idx);
}


static int
ged_draw_scene_ref_foreach_child(bsg_scene_ref ref,
				 int (*cb)(bsg_scene_ref child_ref, void *userdata),
				 void *userdata)
{
    if (ged_draw_scene_ref_is_null(ref) || !cb)
	return 1;

    size_t child_count = ged_draw_scene_ref_child_count(ref);
    for (size_t i = 0; i < child_count; i++) {
	if (!(*cb)(ged_draw_scene_ref_child_at(ref, i), userdata))
	    return 0;
    }

    return 1;
}


struct _ged_draw_source_shape_iter_ctx {
    int has_record_ancestor;
    int (*cb)(bsg_scene_ref, void *);
    void *userdata;
};


static int ged_draw_source_scene_ref_foreach_shape_impl(
	bsg_scene_ref ref,
	int has_record_ancestor,
	int (*cb)(bsg_scene_ref, void *),
	void *userdata);


static int
ged_draw_source_scene_ref_foreach_shape_child_cb(bsg_scene_ref child_ref,
						 void *userdata)
{
    struct _ged_draw_source_shape_iter_ctx *ctx =
	(struct _ged_draw_source_shape_iter_ctx *)userdata;
    return ctx ? ged_draw_source_scene_ref_foreach_shape_impl(child_ref,
	    ctx->has_record_ancestor, ctx->cb, ctx->userdata) : 1;
}


static int
ged_draw_source_scene_ref_foreach_shape_impl(
	bsg_scene_ref ref,
	int has_record_ancestor,
	int (*cb)(bsg_scene_ref, void *),
	void *userdata)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 1;

    struct ged_draw_shape_record_summary shape_summary;
    int has_record = ged_draw_scene_ref_shape_record_summary(ref,
	    &shape_summary);
    if (has_record) {
	int has_path = shape_summary.fullpath &&
	    shape_summary.fullpath->fp_len > 0;
	if ((has_path || !has_record_ancestor) && !(*cb)(ref, userdata))
	    return 0;
    }

    struct _ged_draw_source_shape_iter_ctx args;
    args.has_record_ancestor = has_record_ancestor || has_record;
    args.cb = cb;
    args.userdata = userdata;
    return ged_draw_scene_ref_foreach_child(ref,
	    ged_draw_source_scene_ref_foreach_shape_child_cb, &args);
}


static int
ged_draw_source_scene_ref_foreach_shape(bsg_scene_ref ref,
					int (*cb)(bsg_scene_ref, void *),
					void *userdata)
{
    if (!cb)
	return 1;
    return ged_draw_source_scene_ref_foreach_shape_impl(ref, 0, cb, userdata);
}


struct ged_draw_source_root_shape_ref_ctx {
    struct ged *gedp;
    ged_draw_shape_ref_index_cb cb;
    void *userdata;
};


static int
_ged_draw_source_root_shape_ref_cb(bsg_scene_ref shape_ref,
				   void *userdata)
{
    struct ged_draw_source_root_shape_ref_ctx *ctx =
	(struct ged_draw_source_root_shape_ref_ctx *)userdata;
    if (!ctx || !ctx->cb)
	return 1;

    ged_draw_shape_ref ref = _ged_draw_shape_ref_from_scene_ref(ctx->gedp,
	    shape_ref);
    if (ged_draw_shape_ref_is_null(ref))
	return 1;
    return ctx->cb(ref, ctx->userdata);
}


static ged_draw_shape_ref
ged_draw_obol_shape_ref_for_path_mode(struct ged *gedp,
				      const char *path,
				      int require_database_source,
				      int draw_mode_valid,
				      int draw_mode)
{
    if (!gedp || !path || !path[0])
	return GED_DRAW_SHAPE_REF_NULL;

    void *scene_ctx = ged_draw_obol_context_token_with_parent_for_path_mode(
	    gedp, path, draw_mode_valid, draw_mode);
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_ctx(scene_ctx);
    if (!token || (require_database_source && !token->is_database_source) ||
	    (!require_database_source && !token->is_database_source &&
	     !token->is_shape))
	return GED_DRAW_SHAPE_REF_NULL;

    rt_view_scene_ref scene_ref = rt_view_scene_ref_make(scene_ctx,
	    RT_VIEW_SCENE_BACKEND_OBOL);
    ged_draw_shape_ref ref =
	ged_draw_registry_shape_ref_from_source_ref(gedp, scene_ref);
    if (ged_draw_shape_ref_is_null(ref))
	return GED_DRAW_SHAPE_REF_NULL;

    if (token->fullpath_valid)
	(void)ged_draw_registry_shape_ref_set_indexed_fullpath(gedp, ref,
		&token->fullpath);

    return ref;
}


static ged_draw_shape_ref
ged_draw_obol_shape_ref_for_path(struct ged *gedp,
				 const char *path,
				 int require_database_source)
{
    return ged_draw_obol_shape_ref_for_path_mode(gedp, path,
	    require_database_source, 0, GED_DRAW_MODE_WIRE);
}


static ged_draw_shape_ref
ged_draw_obol_shape_ref_for_database_source_record(
	struct ged *gedp,
	const struct ged_draw_obol_database_source_record *record)
{
    if (!record || !record->database_path || !record->database_path[0])
	return GED_DRAW_SHAPE_REF_NULL;
    return ged_draw_obol_shape_ref_for_path_mode(gedp, record->database_path,
	    1, record->valid, record->draw_mode);
}


static int
_ged_draw_source_root_obol_shape_path_cb(struct ged *gedp,
					 const char *path,
					 void *userdata)
{
    struct ged_draw_source_root_shape_ref_ctx *ctx =
	(struct ged_draw_source_root_shape_ref_ctx *)userdata;
    if (!gedp || !path || !path[0] || !ctx || !ctx->cb)
	return 1;

    ged_draw_shape_ref ref =
	ged_draw_obol_shape_ref_for_path(gedp, path, 0);
    if (ged_draw_shape_ref_is_null(ref))
	return 1;

    return ctx->cb(ref, ctx->userdata);
}


static int
_ged_draw_source_root_obol_shape_record_cb(
	struct ged *gedp,
	const struct ged_draw_obol_database_source_record *record,
	void *userdata)
{
    struct ged_draw_source_root_shape_ref_ctx *ctx =
	(struct ged_draw_source_root_shape_ref_ctx *)userdata;
    if (!gedp || !record || !record->database_path ||
	    !record->database_path[0] || !ctx || !ctx->cb)
	return 1;

    ged_draw_shape_ref ref =
	ged_draw_obol_shape_ref_for_database_source_record(gedp, record);
    if (ged_draw_shape_ref_is_null(ref))
	return 1;

    return ctx->cb(ref, ctx->userdata);
}


struct ged_draw_source_root_shape_child_ctx {
    struct ged_draw_source_root_shape_ref_ctx shape_ctx;
    int skip_overlay_groups;
};


static int
_ged_draw_source_root_shape_child_cb(bsg_scene_ref child_ref,
				     void *userdata)
{
    struct ged_draw_source_root_shape_child_ctx *ctx =
	(struct ged_draw_source_root_shape_child_ctx *)userdata;
    if (!ctx)
	return 1;

    if (ctx->skip_overlay_groups) {
	struct ged_draw_group_record_summary group_summary;
	if (ged_draw_scene_ref_group_record_summary(child_ref,
		&group_summary) &&
		group_summary.is_overlay)
	    return 1;
    }

    return ged_draw_source_scene_ref_foreach_shape(child_ref,
	    _ged_draw_source_root_shape_ref_cb, &ctx->shape_ctx);
}


int
ged_draw_source_root_foreach_shape_ref(struct ged *gedp,
				       int skip_overlay_groups,
				       ged_draw_shape_ref_index_cb cb,
				       void *userdata)
{
    if (!gedp || !cb)
	return 1;

    struct ged_draw_source_root_shape_ref_ctx obol_ctx;
    obol_ctx.gedp = gedp;
    obol_ctx.cb = cb;
    obol_ctx.userdata = userdata;
    if (ged_draw_obol_scene_controller_full_synced(gedp)) {
	int obol_status = ged_draw_obol_database_source_records_foreach(gedp,
		skip_overlay_groups, _ged_draw_source_root_obol_shape_record_cb,
		&obol_ctx);
	if (obol_status == 0)
	    return 0;
	int shape_status = ged_draw_obol_shape_paths_foreach(gedp,
		skip_overlay_groups, _ged_draw_source_root_obol_shape_path_cb,
		&obol_ctx);
	if (shape_status >= 0)
	    return shape_status;
	if (obol_status >= 0)
	    return obol_status;
    }

    bsg_scene_ref root_ref = _ged_draw_source_root_scene_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return 1;

    struct ged_draw_source_root_shape_child_ctx ctx;
    ctx.shape_ctx.gedp = gedp;
    ctx.shape_ctx.cb = cb;
    ctx.shape_ctx.userdata = userdata;
    ctx.skip_overlay_groups = skip_overlay_groups;
    return ged_draw_scene_ref_foreach_child(root_ref,
	    _ged_draw_source_root_shape_child_cb, &ctx);
}


struct ged_draw_source_root_group_ref_ctx {
    struct ged *gedp;
    ged_draw_group_ref_index_cb cb;
    void *userdata;
};


static int _ged_draw_source_root_foreach_group_scene(
	struct ged_draw_source_root_group_ref_ctx *ctx,
	bsg_scene_ref scene_ref);


static int
_ged_draw_source_root_obol_group_path_cb(struct ged *gedp,
					 const char *path,
					 void *userdata)
{
    struct ged_draw_source_root_group_ref_ctx *ctx =
	(struct ged_draw_source_root_group_ref_ctx *)userdata;
    if (!ctx || !path || !path[0])
	return 1;

    ged_draw_group_ref ref = ged_draw_obol_group_ref_for_path(gedp, path);
    if (ged_draw_group_ref_is_null(ref))
	return 1;

    return ctx->cb(ref, ctx->userdata);
}


static int
_ged_draw_source_root_group_child_cb(bsg_scene_ref child_ref,
				     void *userdata)
{
    struct ged_draw_source_root_group_ref_ctx *ctx =
	(struct ged_draw_source_root_group_ref_ctx *)userdata;
    return ctx ? _ged_draw_source_root_foreach_group_scene(ctx,
	    child_ref) : 1;
}


static int
_ged_draw_source_root_foreach_group_scene(
	struct ged_draw_source_root_group_ref_ctx *ctx,
	bsg_scene_ref scene_ref)
{
    if (!ctx || ged_draw_scene_ref_is_null(scene_ref))
	return 1;

    if (ged_draw_scene_ref_is_group(scene_ref)) {
	ged_draw_group_ref ref = _ged_draw_group_ref_from_scene_ref(ctx->gedp,
		scene_ref);
	if (!ged_draw_group_ref_is_null(ref) && !ctx->cb(ref, ctx->userdata))
	    return 0;
    }

    return ged_draw_scene_ref_foreach_child(scene_ref,
	    _ged_draw_source_root_group_child_cb, ctx);
}


int
ged_draw_source_root_foreach_group_ref(struct ged *gedp,
				       ged_draw_group_ref_index_cb cb,
				       void *userdata)
{
    if (!gedp || !cb)
	return 1;

    struct ged_draw_source_root_group_ref_ctx obol_ctx;
    obol_ctx.gedp = gedp;
    obol_ctx.cb = cb;
    obol_ctx.userdata = userdata;
    int obol_status = ged_draw_obol_scene_controller_full_synced(gedp) ?
	ged_draw_obol_group_paths_foreach(gedp, 0,
	    _ged_draw_source_root_obol_group_path_cb, &obol_ctx) : -1;
    if (obol_status >= 0)
	return obol_status;

    bsg_scene_ref root_ref = _ged_draw_source_root_scene_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return 1;

    struct ged_draw_source_root_group_ref_ctx ctx;
    ctx.gedp = gedp;
    ctx.cb = cb;
    ctx.userdata = userdata;
    return ged_draw_scene_ref_foreach_child(root_ref,
	    _ged_draw_source_root_group_child_cb, &ctx);
}


static int
_ged_draw_source_root_has_group_cb(ged_draw_group_ref UNUSED(group_ref),
				   void *userdata)
{
    int *has_group = (int *)userdata;
    if (has_group)
	*has_group = 1;
    return 0;
}


int
ged_draw_source_root_has_groups(struct ged *gedp)
{
    if (!gedp)
	return 0;

    size_t obol_group_count = 0;
    if (ged_draw_obol_group_descendant_group_count_for_path(gedp, "/",
	    &obol_group_count))
	return obol_group_count > 0 ? 1 : 0;

    int has_group = 0;
    (void)ged_draw_source_root_foreach_group_ref(gedp,
	    _ged_draw_source_root_has_group_cb, &has_group);
    return has_group;
}


int
ged_draw_source_root_subtree_bounds(struct ged *gedp,
				    vect_t *min,
				    vect_t *max,
				    int include_overlays)
{
    if (!gedp)
	return 1;

    if (!include_overlays) {
	int obol_empty = 1;
	if (ged_draw_obol_scene_controller_full_synced(gedp) &&
		ged_draw_obol_scene_database_bounds(gedp, min, max,
		&obol_empty))
	    return obol_empty;
    }

    bsg_scene_ref root_ref = _ged_draw_source_root_scene_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return 1;

    return ged_draw_scene_ref_subtree_bounds(root_ref, min, max,
	    include_overlays);
}


static void _ged_draw_scene_ref_release_data_recurse(bsg_scene_ref ref);


static int
_ged_draw_scene_ref_release_data_recurse_cb(bsg_scene_ref child_ref,
					    void *UNUSED(userdata))
{
    _ged_draw_scene_ref_release_data_recurse(child_ref);
    return 1;
}


static void
_ged_draw_scene_ref_release_data_recurse(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return;

    (void)ged_draw_scene_ref_foreach_child(ref,
	    _ged_draw_scene_ref_release_data_recurse_cb, NULL);
    ged_draw_registry_source_ref_highlight_free(
	    ged_draw_scene_ref_to_rt_view_ref(ref));
}


static void
ged_draw_scene_ref_release(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return;
    _ged_draw_scene_ref_release_data_recurse(ref);
    bsg_scene_ref_destroy(ref);
}


static int
ged_draw_scene_ref_release_source_owner(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 0;

    bsg_scene_ref release_ref = ged_draw_scene_ref_source_owner(ref);
    if (ged_draw_scene_ref_is_null(release_ref))
	release_ref = ref;
    (void)ged_draw_scene_ref_detach(release_ref);
    ged_draw_scene_ref_release(release_ref);
    return 1;
}


static int
ged_draw_scene_ref_detach(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 0;

    bsg_scene_ref parent_ref = ged_draw_scene_ref_parent(ref);
    if (ged_draw_scene_ref_is_null(parent_ref))
	return 0;

    bsg_scene_remove_child(parent_ref, ref);
    return 1;
}


static int
ged_draw_scene_ref_append_child(bsg_scene_ref parent_ref, bsg_scene_ref child_ref)
{
    if (ged_draw_scene_ref_is_null(parent_ref) ||
	    ged_draw_scene_ref_is_null(child_ref))
	return 0;

    bsg_scene_append_child(parent_ref, child_ref);
    bsg_scene_bbox_invalidate(parent_ref);
    return 1;
}


static bsg_scene_ref
ged_draw_view_context_group_child_ensure(void *view_ctx,
					 bsg_scene_ref parent_ref,
					 const char *name,
					 void *dp_hint)
{
    if (ged_draw_scene_ref_is_null(parent_ref) || !name)
	return bsg_scene_ref_null();

    bsg_scene_ref existing = bsg_scene_group_find_child(parent_ref, name);
    if (!ged_draw_scene_ref_is_null(existing))
	return existing;

    if (!view_ctx)
	return bsg_scene_ref_null();

    bsg_scene_ref child_ref =
	bsg_scene_group_ensure_child(parent_ref, (struct bsg_view *)view_ctx,
		name, dp_hint);
    if (ged_draw_scene_ref_is_null(child_ref))
	return bsg_scene_ref_null();

    bsg_scene_bump_rev(parent_ref);
    return child_ref;
}


static rt_view_scene_ref
ged_draw_view_context_group_create_ref(void *view_ctx, const char *name)
{
    if (!view_ctx || !name)
	return rt_view_scene_ref_null();

    return ged_draw_scene_ref_to_rt_view_ref(
	    bsg_scene_group_create((struct bsg_view *)view_ctx, name));
}


static void
ged_draw_scene_ref_erase_nested_subpath(
	bsg_scene_ref parent_ref,
	const char * const *comp_names,
	size_t comp_count,
	size_t depth_start,
	int (*match_fn)(bsg_scene_ref shape_ref, void *match_ctx),
	void *match_ctx)
{
    if (ged_draw_scene_ref_is_null(parent_ref) || !comp_names ||
	    comp_count == 0 || depth_start >= comp_count)
	return;

    bsg_scene_erase_nested_subpath(parent_ref, comp_names, comp_count,
	    depth_start, match_fn, match_ctx);
}


static char *
ged_draw_db_full_path_range_string(const struct db_full_path *fp,
				   size_t first_idx,
				   size_t end_idx)
{
    if (!fp || first_idx >= fp->fp_len || end_idx <= first_idx)
	return NULL;

    if (end_idx > fp->fp_len)
	end_idx = fp->fp_len;

    struct bu_vls path = BU_VLS_INIT_ZERO;
    for (size_t i = first_idx; i < end_idx; i++) {
	struct directory *dp = fp->fp_names[i];
	if (!dp || !dp->d_namep || !dp->d_namep[0])
	    continue;
	if (bu_vls_strlen(&path))
	    bu_vls_putc(&path, '/');
	bu_vls_strcat(&path, dp->d_namep);
    }

    char *ret = NULL;
    if (bu_vls_strlen(&path))
	ret = bu_strdup(bu_vls_cstr(&path));
    bu_vls_free(&path);
    return ret;
}


static int
ged_draw_obol_group_erase_nested_db_full_path(struct ged *gedp,
					      const struct db_full_path *fp)
{
    if (!gedp || !fp || fp->fp_len < 2)
	return 0;

    char *parent_path = ged_draw_db_full_path_range_string(fp, 0,
	    fp->fp_len - 1);
    char *subpath = ged_draw_db_full_path_range_string(fp, fp->fp_len - 1,
	    fp->fp_len);
    int erased = parent_path && subpath ?
	ged_draw_obol_group_erase_subpath_for_path(gedp, parent_path,
		subpath) : 0;

    if (parent_path)
	bu_free(parent_path, "GED Obol nested group parent path");
    if (subpath)
	bu_free(subpath, "GED Obol nested group subpath");

    return erased;
}


static int
ged_draw_scene_ref_nested_subpath_match_cb(bsg_scene_ref shape_ref,
					   void *match_ctx)
{
    if (ged_draw_scene_ref_is_null(shape_ref) || !match_ctx)
	return 0;

    struct ged_draw_shape_record_summary shape_summary;
    if (!ged_draw_scene_ref_shape_record_summary(shape_ref, &shape_summary) ||
	    !shape_summary.fullpath)
	return 0;

    const struct db_full_path *subpath =
	(const struct db_full_path *)match_ctx;
    return db_full_path_match_top(subpath, shape_summary.fullpath);
}


static int
ged_draw_scene_ref_erase_nested_db_subpath(bsg_scene_ref parent_ref,
					   const struct db_full_path *subpath,
					   size_t depth_start)
{
    if (ged_draw_scene_ref_is_null(parent_ref) || !subpath ||
	    subpath->fp_len == 0 || depth_start >= subpath->fp_len)
	return 0;

    const char **names = (const char **)bu_malloc(
	    sizeof(const char *) * subpath->fp_len,
	    "draw source nested subpath names");
    for (size_t i = 0; i < subpath->fp_len; i++)
	names[i] = subpath->fp_names[i]->d_namep;

    char *obol_subpath = ged_draw_db_full_path_range_string(subpath,
	    depth_start, subpath->fp_len);
    if (obol_subpath) {
	(void)ged_draw_scene_ref_obol_group_erase_subpath(parent_ref,
		obol_subpath);
	bu_free(obol_subpath, "GED Obol nested group subpath");
    }

    ged_draw_scene_ref_erase_nested_subpath(parent_ref, names,
	    subpath->fp_len, depth_start,
	    ged_draw_scene_ref_nested_subpath_match_cb, (void *)subpath);

    bu_free(names, "draw source nested subpath names");

    struct ged_draw_scene_tree_summary parent_summary;
    if (ged_draw_scene_ref_tree_summary(parent_ref, &parent_summary) &&
	    parent_summary.child_count == 0) {
	ged_draw_scene_ref_free_group(parent_ref);
	return 2;
    }

    return 1;
}


static void
ged_draw_scene_ref_free_group(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return;

    (void)ged_draw_scene_ref_obol_group_remove(ref);
    bsg_scene_free_group(ref);
}


static void
ged_draw_scene_ref_free_group_contents(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return;

    (void)ged_draw_scene_ref_obol_group_clear(ref);
    bsg_scene_free_group_contents(ref);
}


struct _ged_draw_clear_scope_snapshot {
    bsg_scene_ref *refs;
    size_t len;
    size_t cap;
};


static void
_ged_draw_clear_scope_snapshot_append(struct _ged_draw_clear_scope_snapshot *snap,
				      bsg_scene_ref ref)
{
    if (!snap || ged_draw_scene_ref_is_null(ref))
	return;

    if (snap->len == snap->cap) {
	size_t new_cap = snap->cap ? snap->cap * 2 : 16;
	snap->refs = (bsg_scene_ref *)bu_realloc(snap->refs,
		new_cap * sizeof(bsg_scene_ref),
		"draw source clear scope snapshot");
	snap->cap = new_cap;
    }

    snap->refs[snap->len++] = ref;
}


static void
_ged_draw_clear_scope_snapshot_free(struct _ged_draw_clear_scope_snapshot *snap)
{
    if (!snap)
	return;
    if (snap->refs)
	bu_free(snap->refs, "draw source clear scope snapshot");
    snap->refs = NULL;
    snap->len = 0;
    snap->cap = 0;
}


struct _ged_draw_erase_subgroups_by_name_ctx {
    struct _ged_draw_clear_scope_snapshot *snap;
};


static int
_ged_draw_erase_subgroups_by_name_child_cb(bsg_scene_ref child,
					   void *userdata)
{
    struct _ged_draw_erase_subgroups_by_name_ctx *ctx =
	(struct _ged_draw_erase_subgroups_by_name_ctx *)userdata;
    if (!ctx || ged_draw_scene_ref_is_null(child))
	return 1;

    struct ged_draw_scene_tree_summary child_summary;
    if (!ged_draw_scene_ref_tree_summary(child, &child_summary) ||
	    !child_summary.is_group)
	return 1;

    _ged_draw_clear_scope_snapshot_append(ctx->snap, child);
    return 1;
}


static int
ged_draw_scene_ref_erase_subgroups_by_name(bsg_scene_ref parent_ref,
					   const char *name)
{
    if (ged_draw_scene_ref_is_null(parent_ref) || !name)
	return 0;

    struct _ged_draw_clear_scope_snapshot snap = {NULL, 0, 0};
    struct _ged_draw_erase_subgroups_by_name_ctx ctx;
    ctx.snap = &snap;
    (void)ged_draw_scene_ref_foreach_child(parent_ref,
	    _ged_draw_erase_subgroups_by_name_child_cb, &ctx);

    int erased = 0;
    for (size_t i = 0; i < snap.len; i++) {
	bsg_scene_ref child_ref = snap.refs[i];
	struct ged_draw_scene_tree_summary child_summary;
	if (ged_draw_scene_ref_tree_summary(child_ref, &child_summary) &&
		child_summary.name && BU_STR_EQUAL(child_summary.name, name)) {
	    ged_draw_scene_ref_free_group(child_ref);
	    erased++;
	} else {
	    erased += ged_draw_scene_ref_erase_subgroups_by_name(child_ref,
		    name);
	}
    }

    _ged_draw_clear_scope_snapshot_free(&snap);
    return erased;
}


static int
_ged_draw_scene_path_has_component_name(const char *path,
					const char *name,
					size_t first_idx)
{
    if (!path || !name)
	return 0;

    name = ged_draw_dbpath_skip_lead_slash(name);
    size_t namelen = strlen(name);
    if (!namelen)
	return 0;

    size_t idx = 0;
    const char *p = path;
    while (*p) {
	while (*p == '/')
	    p++;
	if (!*p)
	    break;
	const char *slash = strchr(p, '/');
	size_t len = slash ? (size_t)(slash - p) : strlen(p);
	if (idx >= first_idx && len == namelen &&
		bu_strncmp(p, name, len) == 0)
	    return 1;
	if (!slash)
	    break;
	p = slash + 1;
	idx++;
    }

    return 0;
}


static int
_ged_draw_fullpath_has_component_name(const struct db_full_path *fp,
				      const char *name,
				      size_t first_idx)
{
    if (!fp || !name)
	return 0;

    name = ged_draw_dbpath_skip_lead_slash(name);
    size_t namelen = strlen(name);
    if (!namelen)
	return 0;

    for (size_t i = first_idx; i < fp->fp_len; i++) {
	struct directory *dp = fp->fp_names[i];
	if (dp && strlen(dp->d_namep) == namelen &&
		bu_strncmp(dp->d_namep, name, namelen) == 0)
	    return 1;
    }

    return 0;
}


static void
ged_draw_obol_owner_structural_revision_bump(struct ged *gedp)
{
    if (gedp && gedp->i && gedp->i->ged_gdp)
	gedp->i->ged_gdp->gd_draw_rev++;
}


struct ged_draw_obol_shape_component_index_ctx {
    struct ged *gedp;
    const char *path;
    ged_draw_shape_ref_index_cb cb;
    void *userdata;
    int count;
};


static int
ged_draw_obol_shape_component_index_record_cb(
	struct ged *gedp,
	const struct ged_draw_obol_database_source_record *record,
	void *userdata)
{
    struct ged_draw_obol_shape_component_index_ctx *ctx =
	(struct ged_draw_obol_shape_component_index_ctx *)userdata;
    const char *source_path = record ? record->database_path : NULL;
    if (!ctx || !source_path || !source_path[0])
	return 1;
    if (!_ged_draw_scene_path_has_component_name(source_path, ctx->path, 0))
	return 1;

    ctx->count++;
    if (!ctx->cb)
	return 1;

    ged_draw_shape_ref ref =
	ged_draw_obol_shape_ref_for_database_source_record(gedp, record);
    if (ged_draw_shape_ref_is_null(ref))
	return 1;

    return ctx->cb(ref, ctx->userdata);
}


int
ged_draw_obol_shape_ref_index_for_component(
	struct ged *gedp,
	const char *path,
	ged_draw_shape_ref_index_cb cb,
	void *userdata)
{
    if (!gedp || !path || !path[0])
	return -1;

    struct ged_draw_obol_shape_component_index_ctx ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.gedp = gedp;
    ctx.path = path;
    ctx.cb = cb;
    ctx.userdata = userdata;
    int status = ged_draw_obol_scene_controller_full_synced(gedp) ?
	ged_draw_obol_database_source_records_foreach(gedp, 0,
	    ged_draw_obol_shape_component_index_record_cb, &ctx) : -1;
    return status >= 0 ? ctx.count : -1;
}


static unsigned long long
ged_draw_obol_database_source_path_hash(struct ged *gedp, const char *path)
{
    if (!gedp || !gedp->dbip || !path || !path[0])
	return 0;

    struct db_full_path dfp;
    db_full_path_init(&dfp);
    if (db_string_to_path(&dfp, gedp->dbip, path) != 0) {
	db_free_full_path(&dfp);
	return 0;
    }
    if (dfp.fp_len <= 0) {
	db_free_full_path(&dfp);
	return 0;
    }

    ged_db_index_id *components =
	(ged_db_index_id *)bu_calloc(dfp.fp_len, sizeof(ged_db_index_id),
		"GED Obol shape path hash components");
    for (size_t i = 0; i < dfp.fp_len; i++) {
	const struct directory *dp = dfp.fp_names[i];
	const char *name = (dp && dp->d_namep) ? dp->d_namep : "";
	components[i] = bu_data_hash(name, strlen(name) * sizeof(char));
    }

    unsigned long long hash = ged_db_index_path_hash(gedp, components,
	    dfp.fp_len, 0);
    bu_free(components, "GED Obol shape path hash components");
    db_free_full_path(&dfp);
    return hash;
}


struct ged_draw_obol_shape_path_hash_index_ctx {
    struct ged *gedp;
    unsigned long long path_hash;
    ged_draw_shape_ref_index_cb cb;
    void *userdata;
    int count;
};


static int
ged_draw_obol_shape_path_hash_index_record_cb(
	struct ged *gedp,
	const struct ged_draw_obol_database_source_record *record,
	void *userdata)
{
    struct ged_draw_obol_shape_path_hash_index_ctx *ctx =
	(struct ged_draw_obol_shape_path_hash_index_ctx *)userdata;
    const char *source_path = record ? record->database_path : NULL;
    if (!ctx || !source_path || !source_path[0])
	return 1;

    unsigned long long source_hash =
	ged_draw_obol_database_source_path_hash(gedp, source_path);
    if (!source_hash || source_hash != ctx->path_hash)
	return 1;

    ctx->count++;
    if (!ctx->cb)
	return 1;

    ged_draw_shape_ref ref =
	ged_draw_obol_shape_ref_for_database_source_record(gedp, record);
    if (ged_draw_shape_ref_is_null(ref))
	return 1;

    return ctx->cb(ref, ctx->userdata);
}


int
ged_draw_obol_shape_ref_index_for_path_hash(
	struct ged *gedp,
	unsigned long long path_hash,
	ged_draw_shape_ref_index_cb cb,
	void *userdata)
{
    if (!gedp || !path_hash)
	return -1;

    struct ged_draw_obol_shape_path_hash_index_ctx ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.gedp = gedp;
    ctx.path_hash = path_hash;
    ctx.cb = cb;
    ctx.userdata = userdata;
    int status = ged_draw_obol_scene_controller_full_synced(gedp) ?
	ged_draw_obol_database_source_records_foreach(gedp, 0,
	    ged_draw_obol_shape_path_hash_index_record_cb, &ctx) : -1;
    return status >= 0 ? ctx.count : -1;
}


struct ged_draw_obol_group_component_index_ctx {
    struct ged *gedp;
    const char *path;
    ged_draw_group_ref_index_cb cb;
    void *userdata;
    int count;
};


static int
ged_draw_obol_group_component_index_path_cb(struct ged *gedp,
					    const char *group_path,
					    void *userdata)
{
    struct ged_draw_obol_group_component_index_ctx *ctx =
	(struct ged_draw_obol_group_component_index_ctx *)userdata;
    if (!ctx || !group_path || !group_path[0])
	return 1;
    if (!_ged_draw_scene_path_has_component_name(group_path, ctx->path, 0))
	return 1;

    ctx->count++;
    if (!ctx->cb)
	return 1;

    ged_draw_group_ref ref = ged_draw_obol_group_ref_for_path(gedp,
	    group_path);
    if (ged_draw_group_ref_is_null(ref))
	return 1;

    return ctx->cb(ref, ctx->userdata);
}


int
ged_draw_obol_group_ref_index_for_component(
	struct ged *gedp,
	const char *path,
	ged_draw_group_ref_index_cb cb,
	void *userdata)
{
    if (!gedp || !path || !path[0])
	return -1;

    struct ged_draw_obol_group_component_index_ctx ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.gedp = gedp;
    ctx.path = path;
    ctx.cb = cb;
    ctx.userdata = userdata;
    int status = ged_draw_obol_scene_controller_full_synced(gedp) ?
	ged_draw_obol_group_paths_foreach(gedp, 0,
	    ged_draw_obol_group_component_index_path_cb, &ctx) : -1;
    return status >= 0 ? ctx.count : -1;
}


static int
_ged_draw_clear_scope_child_cb(bsg_scene_ref child, void *userdata)
{
    _ged_draw_clear_scope_snapshot_append(
	    (struct _ged_draw_clear_scope_snapshot *)userdata, child);
    return 1;
}


static int
ged_draw_scene_ref_clear_scope_children(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 0;

    struct _ged_draw_clear_scope_snapshot snap = {NULL, 0, 0};
    (void)ged_draw_scene_ref_foreach_child(ref, _ged_draw_clear_scope_child_cb,
	    &snap);

    int cleared = 0;
    for (size_t i = 0; i < snap.len; i++) {
	bsg_scene_ref child_ref = snap.refs[i];
	(void)ged_draw_scene_ref_obol_group_remove(child_ref);
	ged_draw_scene_ref_free_group_contents(child_ref);
	(void)ged_draw_scene_ref_detach(child_ref);
	ged_draw_scene_ref_release(child_ref);
	cleared++;
    }

    _ged_draw_clear_scope_snapshot_free(&snap);
    return cleared;
}


int
ged_draw_source_root_clear_all_scope_children(struct ged *gedp)
{
    if (!gedp)
	return 0;

    int obol_cleared = ged_draw_obol_scene_controller_full_synced(gedp) ?
	ged_draw_obol_active_scene_clear(gedp) : 0;
    if (obol_cleared) {
	/* Obol owns source-root clear; retained clear is fallback. */
	ged_draw_obol_owner_structural_revision_bump(gedp);
	return 1;
    }
    if (!ged_draw_obol_scene_controller_full_synced(gedp))
	(void)ged_draw_obol_active_scene_clear(gedp);

    int cleared = 0;
    struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
    if (views) {
	for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	    void *view_ctx = (void *)BU_PTBL_GET(views, i);
	    if (!view_ctx || !rt_view_context_is_independent(view_ctx))
		continue;
	    bsg_scene_ref scope_ref =
		_ged_draw_active_scope_scene_ref(gedp, view_ctx, 0, 0);
	    cleared += ged_draw_scene_ref_clear_scope_children(scope_ref);
	}
    }

    cleared += ged_draw_scene_ref_clear_scope_children(
	    _ged_draw_source_root_scene_ref(gedp));
    return cleared || obol_cleared;
}


struct _ged_draw_scene_rt_ref_foreach_child_ctx {
    ged_draw_scene_rt_ref_child_cb cb;
    void *userdata;
};


static int
_ged_draw_scene_rt_ref_foreach_child_cb(bsg_scene_ref child_ref,
					void *userdata)
{
    struct _ged_draw_scene_rt_ref_foreach_child_ctx *ctx =
	(struct _ged_draw_scene_rt_ref_foreach_child_ctx *)userdata;
    if (!ctx || !ctx->cb)
	return 1;

    return (*ctx->cb)(ged_draw_scene_ref_to_rt_view_ref(child_ref),
	    ctx->userdata);
}


static int
ged_draw_scene_rt_ref_foreach_child(rt_view_scene_ref scene_ref,
				    ged_draw_scene_rt_ref_child_cb cb,
				    void *userdata)
{
    bsg_scene_ref ref = ged_draw_scene_ref_from_rt_view_ref(scene_ref);
    if (ged_draw_scene_ref_is_null(ref) || !cb)
	return 1;

    struct _ged_draw_scene_rt_ref_foreach_child_ctx ctx;
    ctx.cb = cb;
    ctx.userdata = userdata;
    return ged_draw_scene_ref_foreach_child(
	    ref, _ged_draw_scene_rt_ref_foreach_child_cb, &ctx);
}


static int
ged_draw_scene_ref_group_matches_erase_scope(struct ged *gedp,
					     bsg_scene_ref group_ref,
					     void *view_ctx,
					     int mode)
{
    if (!gedp || ged_draw_scene_ref_is_null(group_ref))
	return 0;

    struct ged_draw_group_record_summary group_summary;
    if (!ged_draw_scene_ref_group_record_summary(group_ref, &group_summary))
	return 0;
    if (group_summary.is_overlay)
	return 0;
    if (!view_ctx) {
	if (group_summary.in_view_scope)
	    return 0;
    } else if (rt_view_context_is_independent(view_ctx)) {
	if (!group_summary.in_view_scope || group_summary.view_ctx != view_ctx)
	    return 0;
    } else if (group_summary.in_view_scope && group_summary.view_ctx != view_ctx) {
	return 0;
    }
    if (mode >= 0 && (int)group_summary.draw_mode != mode)
	return 0;
    return 1;
}


static int
ged_draw_scene_ref_group_matches_optional_scope(struct ged *gedp,
						bsg_scene_ref group_ref,
						void *view_ctx,
						int mode)
{
    if (!gedp || ged_draw_scene_ref_is_null(group_ref))
	return 0;

    struct ged_draw_group_record_summary group_summary;
    if (!ged_draw_scene_ref_group_record_summary(group_ref, &group_summary))
	return 0;
    if (group_summary.is_overlay)
	return 0;
    if (view_ctx) {
	if (rt_view_context_is_independent(view_ctx)) {
	    if (!group_summary.in_view_scope || group_summary.view_ctx != view_ctx)
		return 0;
	} else if (group_summary.in_view_scope && group_summary.view_ctx != view_ctx) {
	    return 0;
	}
    }
    if (mode >= 0 && (int)group_summary.draw_mode != mode)
	return 0;
    return 1;
}


static size_t
_ged_draw_scene_path_norm_len(const char *path)
{
    if (!path)
	return 0;
    size_t len = strlen(path);
    while (len && path[len - 1] == '/')
	len--;
    return len;
}


static int
_ged_draw_scene_path_string_is_prefix(const char *prefix, const char *path)
{
    if (!prefix || !path)
	return 0;

    prefix = ged_draw_dbpath_skip_lead_slash(prefix);
    path = ged_draw_dbpath_skip_lead_slash(path);

    size_t plen = _ged_draw_scene_path_norm_len(prefix);
    size_t path_len = _ged_draw_scene_path_norm_len(path);
    if (!plen || path_len < plen)
	return 0;
    if (bu_strncmp(prefix, path, plen) != 0)
	return 0;
    return path_len == plen || path[plen] == '/';
}


struct _ged_draw_prune_groups_ctx {
    struct _ged_draw_clear_scope_snapshot *snap;
};


static int
_ged_draw_prune_groups_child_cb(bsg_scene_ref child, void *userdata)
{
    struct _ged_draw_prune_groups_ctx *ctx =
	(struct _ged_draw_prune_groups_ctx *)userdata;
    if (!ctx || ged_draw_scene_ref_is_null(child))
	return 1;

    struct ged_draw_scene_tree_summary child_summary;
    if (!ged_draw_scene_ref_tree_summary(child, &child_summary) ||
	    !child_summary.is_group)
	return 1;

    struct ged_draw_group_record_summary group_summary;
    if (ged_draw_scene_ref_group_record_summary(child, &group_summary) &&
	    group_summary.is_overlay)
	return 1;

    _ged_draw_clear_scope_snapshot_append(ctx->snap, child);
    return 1;
}


static int
_ged_draw_scene_ref_erase_groups_by_name(bsg_scene_ref parent_ref,
					 const char *name)
{
    if (ged_draw_scene_ref_is_null(parent_ref) || !name)
	return 0;

    struct _ged_draw_clear_scope_snapshot snap = {NULL, 0, 0};
    struct _ged_draw_prune_groups_ctx ctx;
    ctx.snap = &snap;
    (void)ged_draw_scene_ref_foreach_child(parent_ref,
	    _ged_draw_prune_groups_child_cb, &ctx);

    int erased = 0;
    for (size_t i = 0; i < snap.len; i++) {
	bsg_scene_ref group_ref = snap.refs[i];
	struct ged_draw_group_record_summary group_summary;
	if (!ged_draw_scene_ref_group_record_summary(group_ref,
		&group_summary) || !group_summary.path)
	    continue;

	if (_ged_draw_scene_path_has_component_name(group_summary.path, name,
		0)) {
	    ged_draw_scene_ref_free_group(group_ref);
	    erased++;
	} else {
	    erased += ged_draw_scene_ref_erase_subgroups_by_name(group_ref,
		    name);
	}
    }

    _ged_draw_clear_scope_snapshot_free(&snap);
    return erased;
}


int
ged_draw_source_erase_groups_by_name_at_root(struct ged *gedp,
					     const char *name)
{
    if (!gedp || !name)
	return 0;

    int obol_erased =
	ged_draw_obol_groups_remove_for_component_name(gedp, name);
    if (obol_erased) {
	/* Obol owns group-name erase; retained group erase is fallback. */
	ged_draw_obol_owner_structural_revision_bump(gedp);
	return 1;
    }

    bsg_scene_ref root_ref = _ged_draw_source_root_scene_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return obol_erased;

    int bsg_erased = _ged_draw_scene_ref_erase_groups_by_name(root_ref,
	    name);
    return bsg_erased || obol_erased;
}


static int
_ged_draw_group_child_cb(bsg_scene_ref child, void *userdata)
{
    struct _ged_draw_clear_scope_snapshot *snap =
	(struct _ged_draw_clear_scope_snapshot *)userdata;
    if (!snap || ged_draw_scene_ref_is_null(child))
	return 1;

    struct ged_draw_scene_tree_summary child_summary;
    if (!ged_draw_scene_ref_tree_summary(child, &child_summary) ||
	    !child_summary.is_group)
	return 1;

    _ged_draw_clear_scope_snapshot_append(snap, child);
    return 1;
}


static int
ged_draw_scene_ref_erase_matching_group_path_or_nested(
	struct ged *gedp,
	const char *path,
	bsg_scene_ref base_ref,
	const struct db_full_path *subpath,
	void *view_ctx,
	int mode,
	int apply_scope)
{
    if (!gedp || !path || ged_draw_scene_ref_is_null(base_ref))
	return 0;

    struct _ged_draw_clear_scope_snapshot snap = {NULL, 0, 0};
    struct _ged_draw_prune_groups_ctx ctx;
    ctx.snap = &snap;
    (void)ged_draw_scene_ref_foreach_child(base_ref,
	    _ged_draw_prune_groups_child_cb, &ctx);

    int matched = 0;
    for (size_t i = 0; i < snap.len; i++) {
	bsg_scene_ref group_ref = snap.refs[i];
	if (apply_scope && !ged_draw_scene_ref_group_matches_erase_scope(gedp,
		group_ref, view_ctx, mode))
	    continue;

	struct ged_draw_group_record_summary group_summary;
	if (!ged_draw_scene_ref_group_record_summary(group_ref,
		&group_summary) || !group_summary.path)
	    continue;

	if (BU_STR_EQUAL(path, group_summary.path)) {
	    ged_draw_scene_ref_free_group(group_ref);
	    matched = 1;
	    break;
	}

	if (!subpath)
	    continue;

	struct db_full_path group_path;
	if (db_string_to_path(&group_path, gedp->dbip, group_summary.path) != 0)
	    continue;

	int is_ancestor = db_full_path_match_top(&group_path, subpath);
	size_t ancestor_depth = group_path.fp_len;
	db_free_full_path(&group_path);

	if (is_ancestor && ancestor_depth < subpath->fp_len) {
	    (void)ged_draw_scene_ref_erase_nested_db_subpath(group_ref,
		    subpath, ancestor_depth);
	    matched = 1;
	    break;
	}
    }

    _ged_draw_clear_scope_snapshot_free(&snap);
    return matched;
}


static int
ged_draw_scene_ref_erase_groups_by_db_subpath(struct ged *gedp,
					      bsg_scene_ref base_ref,
					      const struct db_full_path *subpath,
					      void *view_ctx,
					      int mode,
					      int apply_scope)
{
    if (!gedp || ged_draw_scene_ref_is_null(base_ref) || !subpath)
	return 0;

    int erased_total = 0;
    int restart;
    do {
	restart = 0;
	struct _ged_draw_clear_scope_snapshot snap = {NULL, 0, 0};
	(void)ged_draw_scene_ref_foreach_child(base_ref,
		_ged_draw_group_child_cb, &snap);

	for (size_t i = 0; i < snap.len; i++) {
	    bsg_scene_ref group_ref = snap.refs[i];
	    if (apply_scope && !ged_draw_scene_ref_group_matches_erase_scope(gedp,
		    group_ref, view_ctx, mode))
		continue;

	    struct ged_draw_group_record_summary group_summary;
	    if (!ged_draw_scene_ref_group_record_summary(group_ref,
		    &group_summary) || !group_summary.path)
		continue;

	    struct db_full_path group_path;
	    if (db_string_to_path(&group_path, gedp->dbip, group_summary.path) != 0)
		continue;

	    if (db_full_path_subset(&group_path, subpath, 0)) {
		db_free_full_path(&group_path);
		ged_draw_scene_ref_free_group(group_ref);
		erased_total++;
		restart = 1;
		break;
	    }

	    if (db_full_path_match_top(&group_path, subpath) &&
		    group_path.fp_len < subpath->fp_len) {
		size_t depth = group_path.fp_len;
		uint64_t rev0 = gedp->i->ged_gdp->gd_draw_rev;
		db_free_full_path(&group_path);
		int erased = ged_draw_scene_ref_erase_nested_db_subpath(group_ref,
			subpath, depth);
		erased_total += erased;
		if (erased > 1) {
		    restart = 1;
		    break;
		}
		if (gedp->i->ged_gdp->gd_draw_rev != rev0) {
		    restart = 1;
		    break;
		}
		continue;
	    }

	    db_free_full_path(&group_path);
	}

	_ged_draw_clear_scope_snapshot_free(&snap);
    } while (restart);

    return erased_total;
}


int
ged_draw_source_clear_db_groups_in_scope(struct ged *gedp, void *view_ctx)
{
    if (!gedp)
	return 0;

    int obol_cleared = ged_draw_obol_database_sources_clear_in_scope(gedp,
	    view_ctx);
    if (obol_cleared) {
	/* Obol owns scoped database-group clear; retained clear is fallback. */
	ged_draw_obol_owner_structural_revision_bump(gedp);
	return 1;
    }

    bsg_scene_ref base_ref =
	_ged_draw_active_scope_scene_ref(gedp, view_ctx, 0, 0);
    if (ged_draw_scene_ref_is_null(base_ref))
	return obol_cleared;

    struct _ged_draw_clear_scope_snapshot snap = {NULL, 0, 0};
    struct _ged_draw_prune_groups_ctx ctx;
    ctx.snap = &snap;
    (void)ged_draw_scene_ref_foreach_child(base_ref,
	    _ged_draw_prune_groups_child_cb, &ctx);

    int cleared = 0;
    for (size_t i = 0; i < snap.len; i++) {
	bsg_scene_ref group_ref = snap.refs[i];
	if (!ged_draw_scene_ref_group_matches_erase_scope(gedp, group_ref,
		view_ctx, -1))
	    continue;
	ged_draw_scene_ref_free_group(group_ref);
	cleared++;
    }

    _ged_draw_clear_scope_snapshot_free(&snap);
    return cleared || obol_cleared;
}


static int
ged_draw_scene_ref_erase_groups_by_path_prefix_string(struct ged *gedp,
						      const char *path,
						      bsg_scene_ref base_ref,
						      void *view_ctx,
						      int mode)
{
    if (!gedp || !path || ged_draw_scene_ref_is_null(base_ref))
	return 0;

    int erased = 0;
    int restart;
    do {
	restart = 0;
	struct _ged_draw_clear_scope_snapshot snap = {NULL, 0, 0};
	struct _ged_draw_prune_groups_ctx ctx;
	ctx.snap = &snap;
	(void)ged_draw_scene_ref_foreach_child(base_ref,
		_ged_draw_prune_groups_child_cb, &ctx);

	for (size_t i = 0; i < snap.len; i++) {
	    bsg_scene_ref group_ref = snap.refs[i];
	    if (!ged_draw_scene_ref_group_matches_optional_scope(gedp,
		    group_ref, view_ctx, mode))
		continue;

	    struct ged_draw_group_record_summary group_summary;
	    if (!ged_draw_scene_ref_group_record_summary(group_ref,
		    &group_summary) || !group_summary.path)
		continue;

	    if (_ged_draw_scene_path_string_is_prefix(path, group_summary.path)) {
		ged_draw_scene_ref_free_group(group_ref);
		erased++;
		restart = 1;
		break;
	    }
	}

	_ged_draw_clear_scope_snapshot_free(&snap);
    } while (restart);

    return erased;
}


static int
_ged_draw_scene_ref_erase_path_at_base(struct ged *gedp,
				       const char *path,
				       bsg_scene_ref base_ref,
				       void *view_ctx,
				       int mode,
				       int apply_scope)
{
    if (!gedp || !path || ged_draw_scene_ref_is_null(base_ref))
	return 0;

    struct db_full_path subpath;
    int found_subpath = (db_string_to_path(&subpath, gedp->dbip, path) == 0);
    uint64_t erase_rev0 = gedp->i->ged_gdp->gd_draw_rev;
    int obol_group_erased = 0;

    if (mode < 0 && (!view_ctx || !rt_view_context_is_independent(view_ctx))) {
	if (found_subpath)
	    obol_group_erased =
		ged_draw_obol_group_erase_nested_db_full_path(gedp,
			&subpath);
	if (!obol_group_erased)
	    obol_group_erased = ged_draw_obol_group_remove_for_path(gedp,
		    path);
    }

    int erased = ged_draw_scene_ref_erase_matching_group_path_or_nested(gedp,
	    path, base_ref, found_subpath ? &subpath : NULL, view_ctx, mode,
	    apply_scope);
    if (found_subpath) {
	if (gedp->i->ged_gdp->gd_draw_rev == erase_rev0)
	    erased += ged_draw_scene_ref_erase_shapes_by_db_subpath(gedp,
		    &subpath, base_ref, view_ctx, mode);
	db_free_full_path(&subpath);
    }

    return erased || obol_group_erased;
}


static int
ged_draw_obol_erase_exact_path_owner(struct ged *gedp, const char *path,
				     void *view_ctx,
				     int mode)
{
    if (!gedp || !path || !path[0])
	return 0;

    int independent_scope = view_ctx && rt_view_context_is_independent(view_ctx);
    int obol_erased = independent_scope ?
	ged_draw_obol_database_sources_remove_for_path_prefix_in_scope_mode(
		gedp, path, view_ctx, mode) :
	ged_draw_obol_active_database_sources_remove_for_path_prefix_mode(
		gedp, path, mode);

    if (independent_scope) {
	if (obol_erased) {
	    ged_draw_obol_owner_structural_revision_bump(gedp);
	    return 1;
	}
	return 0;
    }

    if (mode >= 0) {
	if (obol_erased) {
	    ged_draw_obol_owner_structural_revision_bump(gedp);
	    return 1;
	}
	return 0;
    }

    if (obol_erased) {
	ged_draw_obol_owner_structural_revision_bump(gedp);
	return 1;
    }

    int obol_group_erased = ged_draw_obol_group_remove_for_path(gedp, path);
    if (obol_group_erased) {
	ged_draw_obol_owner_structural_revision_bump(gedp);
	return 1;
    }

    struct db_full_path subpath;
    int found_subpath = (db_string_to_path(&subpath, gedp->dbip, path) == 0);
    if (found_subpath) {
	obol_group_erased =
	    ged_draw_obol_group_erase_nested_db_full_path(gedp, &subpath);
	db_free_full_path(&subpath);
    }

    if (obol_erased || obol_group_erased) {
	ged_draw_obol_owner_structural_revision_bump(gedp);
	return 1;
    }

    return 0;
}


static int
ged_draw_obol_erase_path_prefix_owner(struct ged *gedp, const char *path,
				      void *view_ctx,
				      int mode)
{
    if (!gedp || !path || !path[0])
	return 0;

    int obol_erased =
	ged_draw_obol_database_sources_remove_for_path_prefix_in_scope_mode(
		gedp, path, view_ctx, mode);
    if (view_ctx && rt_view_context_is_independent(view_ctx)) {
	if (obol_erased) {
	    ged_draw_obol_owner_structural_revision_bump(gedp);
	    return 1;
	}
	return 0;
    }

    if (mode >= 0) {
	if (obol_erased) {
	    ged_draw_obol_owner_structural_revision_bump(gedp);
	    return 1;
	}
	return 0;
    }

    int obol_group_erased =
	ged_draw_obol_groups_remove_for_path_prefix(gedp, path);

    if (obol_erased || obol_group_erased) {
	ged_draw_obol_owner_structural_revision_bump(gedp);
	return 1;
    }

    return 0;
}


int
ged_draw_source_erase_path_at_root(struct ged *gedp, const char *path)
{
    if (!gedp || !path)
	return 0;

    int obol_erased = ged_draw_obol_erase_exact_path_owner(gedp, path,
	    NULL, -1);
    /* Obol owns exact-path source/group erase; retained erase is fallback. */

    bsg_scene_ref root_ref = _ged_draw_source_root_scene_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return obol_erased;

    int bsg_erased = _ged_draw_scene_ref_erase_path_at_base(gedp, path,
	    root_ref, NULL, -1, 0);
    return bsg_erased || obol_erased;
}


static int
_ged_draw_scene_ref_erase_path_prefix_at_base(struct ged *gedp,
					      const char *path,
					      bsg_scene_ref base_ref,
					      void *view_ctx,
					      int mode,
					      int apply_scope)
{
    if (!gedp || !path || ged_draw_scene_ref_is_null(base_ref))
	return 0;

    int erased = ged_draw_scene_ref_erase_groups_by_path_prefix_string(gedp,
	    path, base_ref, view_ctx, mode);
    erased += ged_draw_scene_ref_erase_shapes_by_path_prefix_string(gedp,
	    path, base_ref, view_ctx, mode);
    if (erased)
	return erased;

    if (ged_db_index_available(gedp) &&
	    !ged_db_index_path_resolve(gedp, path, NULL, 0))
	return 0;

    struct db_full_path subpath;
    if (db_string_to_path(&subpath, gedp->dbip, path) != 0)
	return 0;

    uint64_t erase_rev0 = gedp->i->ged_gdp->gd_draw_rev;
    erased += ged_draw_scene_ref_erase_groups_by_db_subpath(gedp, base_ref,
	    &subpath, view_ctx, mode, apply_scope);
    if (gedp->i->ged_gdp->gd_draw_rev == erase_rev0)
	erased += ged_draw_scene_ref_erase_shapes_by_db_subpath(gedp,
		&subpath, base_ref, view_ctx, mode);

    db_free_full_path(&subpath);
    return erased;
}


int
ged_draw_source_erase_path_prefix_at_root(struct ged *gedp, const char *path)
{
    if (!gedp || !path)
	return 0;

    int obol_erased = ged_draw_obol_erase_path_prefix_owner(gedp, path,
	    NULL, -1);
    /* Obol owns path-prefix source/group erase; retained erase is fallback. */

    bsg_scene_ref root_ref = _ged_draw_source_root_scene_ref(gedp);
    if (ged_draw_scene_ref_is_null(root_ref))
	return obol_erased;

    int bsg_erased = _ged_draw_scene_ref_erase_path_prefix_at_base(gedp, path,
	    root_ref, NULL, -1, 0);
    return bsg_erased || obol_erased;
}


int
ged_draw_source_erase_path_in_active_scope(struct ged *gedp,
					   const char *path,
					   void *view_ctx,
					   int mode)
{
    if (!gedp || !path)
	return 0;

    int obol_erased = ged_draw_obol_erase_exact_path_owner(gedp, path,
	    view_ctx, mode);
    if (obol_erased)
	return 1;

    bsg_scene_ref base_ref =
	_ged_draw_active_scope_scene_ref(gedp, view_ctx, 0, 0);
    int bsg_erased = _ged_draw_scene_ref_erase_path_at_base(gedp, path,
	    base_ref, view_ctx, mode, 1);
    return bsg_erased || obol_erased;
}


int
ged_draw_source_erase_path_prefix_in_active_scope(struct ged *gedp,
						  const char *path,
						  void *view_ctx,
						  int mode)
{
    if (!gedp || !path)
	return 0;

    int obol_erased = ged_draw_obol_erase_path_prefix_owner(gedp, path,
	    view_ctx, mode);

    bsg_scene_ref base_ref =
	_ged_draw_active_scope_scene_ref(gedp, view_ctx, 0, 0);
    int bsg_erased = _ged_draw_scene_ref_erase_path_prefix_at_base(gedp, path,
	    base_ref, view_ctx, mode, 1);
    return bsg_erased || obol_erased;
}


struct _ged_draw_erase_path_record_ctx {
    struct ged *gedp;
    const struct db_full_path *subpath;
    void *view_ctx;
    int mode;
    struct bu_ptbl shape_refs;
};


struct _ged_draw_erase_path_string_record_ctx {
    struct ged *gedp;
    const char *path;
    void *view_ctx;
    int mode;
    struct bu_ptbl shape_refs;
};


struct _ged_draw_erase_component_record_ctx {
    struct ged *gedp;
    const char *name;
    void *view_ctx;
    int mode;
    int nonroot_only;
    struct bu_ptbl shape_refs;
};


static int
_ged_draw_scene_ref_erase_path_shape_cb(const struct ged_draw_shape_record *rec,
					void *userdata)
{
    struct _ged_draw_erase_path_record_ctx *ctx =
	(struct _ged_draw_erase_path_record_ctx *)userdata;
    if (!ctx || !rec || !rec->fullpath || !ctx->subpath)
	return 1;
    if (ctx->mode >= 0 && (int)rec->draw_mode != ctx->mode)
	return 1;
    if (!db_full_path_match_top(ctx->subpath, rec->fullpath))
	return 1;

    struct ged_draw_group_record grec;
    if (!ged_draw_group_record_get(ctx->gedp, rec->group, &grec))
	return 1;
    if (grec.is_overlay)
	return 1;
    if (ctx->view_ctx &&
	    !ged_draw_group_record_in_view(&grec, ctx->view_ctx))
	return 1;

    ged_draw_shape_ref *ref = (ged_draw_shape_ref *)bu_malloc(
	    sizeof(ged_draw_shape_ref), "path erase shape ref");
    *ref = rec->ref;
    bu_ptbl_ins(&ctx->shape_refs, (long *)ref);
    return 1;
}


static int
_ged_draw_scene_ref_erase_path_string_shape_cb(
	const struct ged_draw_shape_record *rec,
	void *userdata)
{
    struct _ged_draw_erase_path_string_record_ctx *ctx =
	(struct _ged_draw_erase_path_string_record_ctx *)userdata;
    if (!ctx || !rec || !rec->fullpath || !ctx->path)
	return 1;
    if (ctx->mode >= 0 && (int)rec->draw_mode != ctx->mode)
	return 1;

    char *rec_path = db_path_to_string(rec->fullpath);
    if (!rec_path)
	return 1;
    int matches = _ged_draw_scene_path_string_is_prefix(ctx->path, rec_path);
    bu_free(rec_path, "draw erase record path string");
    if (!matches)
	return 1;

    struct ged_draw_group_record grec;
    if (!ged_draw_group_record_get(ctx->gedp, rec->group, &grec))
	return 1;
    if (grec.is_overlay)
	return 1;
    if (ctx->view_ctx &&
	    !ged_draw_group_record_in_view(&grec, ctx->view_ctx))
	return 1;

    ged_draw_shape_ref *ref = (ged_draw_shape_ref *)bu_malloc(
	    sizeof(ged_draw_shape_ref), "path string erase shape ref");
    *ref = rec->ref;
    bu_ptbl_ins(&ctx->shape_refs, (long *)ref);
    return 1;
}


static int
_ged_draw_scene_ref_erase_component_shape_cb(
	const struct ged_draw_shape_record *rec,
	void *userdata)
{
    struct _ged_draw_erase_component_record_ctx *ctx =
	(struct _ged_draw_erase_component_record_ctx *)userdata;
    if (!ctx || !rec || !rec->fullpath)
	return 1;
    if (ctx->mode >= 0 && (int)rec->draw_mode != ctx->mode)
	return 1;

    struct ged_draw_group_record grec;
    if (!ged_draw_group_record_get(ctx->gedp, rec->group, &grec))
	return 1;
    if (grec.is_overlay)
	return 1;
    if (ctx->view_ctx &&
	    !ged_draw_group_record_in_view(&grec, ctx->view_ctx))
	return 1;

    size_t first_idx = ctx->nonroot_only ? 1 : 0;
    if (rec->display_name &&
	    !_ged_draw_scene_path_has_component_name(rec->display_name,
		ctx->name, first_idx))
	return 1;
    if (!rec->display_name &&
	    !_ged_draw_fullpath_has_component_name(rec->fullpath, ctx->name,
		first_idx))
	return 1;

    ged_draw_shape_ref *ref = (ged_draw_shape_ref *)bu_malloc(
	    sizeof(ged_draw_shape_ref), "component erase shape ref");
    *ref = rec->ref;
    bu_ptbl_ins(&ctx->shape_refs, (long *)ref);
    return 1;
}


static int
_ged_draw_release_shape_refs_and_prune(struct ged *gedp,
				       struct bu_ptbl *shape_refs,
				       bsg_scene_ref prune_base_ref,
				       void *view_ctx,
				       int mode,
				       const char *free_label)
{
    int released = 0;
    for (size_t i = 0; i < BU_PTBL_LEN(shape_refs); i++) {
	ged_draw_shape_ref *ref =
	    (ged_draw_shape_ref *)BU_PTBL_GET(shape_refs, i);
	if (ref) {
	    ref->scene_revision = 0;
	    if (_ged_draw_shape_ref_release(gedp, *ref))
		released++;
	    bu_free(ref, free_label);
	}
    }

    if (released)
	(void)ged_draw_scene_ref_prune_empty_groups(gedp, prune_base_ref,
		view_ctx, mode);

    return released;
}


static int
ged_draw_scene_ref_erase_shapes_by_component_name(
	struct ged *gedp,
	const char *name,
	bsg_scene_ref prune_base_ref,
	void *view_ctx,
	int mode,
	int nonroot_only)
{
    if (!gedp || !name)
	return 0;

    struct _ged_draw_erase_component_record_ctx ctx;
    ctx.gedp = gedp;
    ctx.name = ged_draw_dbpath_skip_lead_slash(name);
    ctx.view_ctx = view_ctx;
    ctx.mode = mode;
    ctx.nonroot_only = nonroot_only ? 1 : 0;
    bu_ptbl_init(&ctx.shape_refs, 8, "component erase shape refs");

    ged_draw_foreach_shape_record(gedp,
	    _ged_draw_scene_ref_erase_component_shape_cb, &ctx);
    int released = _ged_draw_release_shape_refs_and_prune(gedp,
	    &ctx.shape_refs, prune_base_ref, view_ctx, mode,
	    "component erase shape ref");
    bu_ptbl_free(&ctx.shape_refs);
    return released;
}


static int
ged_draw_obol_erase_component_owner(struct ged *gedp,
				    const char *name,
				    void *view_ctx,
				    int mode,
				    int nonroot_only)
{
    if (!gedp || !name || !name[0])
	return 0;

    if (view_ctx && rt_view_context_is_independent(view_ctx))
	return 0;

    int obol_erased =
	ged_draw_obol_database_sources_remove_for_component_name(gedp,
		name, nonroot_only, mode);
    if (obol_erased) {
	ged_draw_obol_owner_structural_revision_bump(gedp);
	return 1;
    }

    return 0;
}


int
ged_draw_source_erase_component_name_in_active_scope(struct ged *gedp,
						    const char *name,
						    void *view_ctx,
						    int mode,
						    int nonroot_only)
{
    if (!gedp || !name)
	return 0;

    int obol_erased = ged_draw_obol_erase_component_owner(gedp, name,
	    view_ctx, mode, nonroot_only);
    if (obol_erased) {
	/* Obol owns shared-scope component erase; retained erase is fallback. */
	return 1;
    }

    bsg_scene_ref base_ref =
	_ged_draw_active_scope_scene_ref(gedp, view_ctx, 0, 0);
    int bsg_erased = ged_draw_scene_ref_erase_component_name_at_base(gedp, name,
	    base_ref, view_ctx, mode, nonroot_only);
    return bsg_erased || obol_erased;
}


static int
ged_draw_scene_ref_erase_component_name_at_base(struct ged *gedp,
						const char *name,
						bsg_scene_ref base_ref,
						void *view_ctx,
						int mode,
						int nonroot_only)
{
    return ged_draw_scene_ref_erase_shapes_by_component_name(gedp, name,
	    base_ref, view_ctx, mode, nonroot_only);
}


static int
ged_draw_scene_ref_erase_shapes_by_db_subpath(
	struct ged *gedp,
	const struct db_full_path *subpath,
	bsg_scene_ref prune_base_ref,
	void *view_ctx,
	int mode)
{
    if (!gedp || !subpath)
	return 0;

    struct _ged_draw_erase_path_record_ctx ctx;
    ctx.gedp = gedp;
    ctx.subpath = subpath;
    ctx.view_ctx = view_ctx;
    ctx.mode = mode;
    bu_ptbl_init(&ctx.shape_refs, 8, "path erase shape refs");

    ged_draw_foreach_shape_record(gedp,
	    _ged_draw_scene_ref_erase_path_shape_cb, &ctx);
    int released = _ged_draw_release_shape_refs_and_prune(gedp,
	    &ctx.shape_refs, prune_base_ref, view_ctx, mode,
	    "path erase shape ref");
    bu_ptbl_free(&ctx.shape_refs);
    return released;
}


static int
ged_draw_scene_ref_erase_shapes_by_path_prefix_string(
	struct ged *gedp,
	const char *path,
	bsg_scene_ref prune_base_ref,
	void *view_ctx,
	int mode)
{
    if (!gedp || !path)
	return 0;

    struct _ged_draw_erase_path_string_record_ctx ctx;
    ctx.gedp = gedp;
    ctx.path = ged_draw_dbpath_skip_lead_slash(path);
    ctx.view_ctx = view_ctx;
    ctx.mode = mode;
    bu_ptbl_init(&ctx.shape_refs, 8, "path string erase shape refs");

    ged_draw_foreach_shape_record(gedp,
	    _ged_draw_scene_ref_erase_path_string_shape_cb, &ctx);
    int released = _ged_draw_release_shape_refs_and_prune(gedp,
	    &ctx.shape_refs, prune_base_ref, view_ctx, mode,
	    "path string erase shape ref");
    bu_ptbl_free(&ctx.shape_refs);
    return released;
}


static int
ged_draw_scene_ref_prune_empty_groups(struct ged *gedp,
				      bsg_scene_ref parent_ref,
				      void *view_ctx,
				      int mode)
{
    if (!gedp || ged_draw_scene_ref_is_null(parent_ref))
	return 0;

    int pruned = 0;
    struct _ged_draw_clear_scope_snapshot snap = {NULL, 0, 0};
    struct _ged_draw_prune_groups_ctx ctx;
    ctx.snap = &snap;
    (void)ged_draw_scene_ref_foreach_child(parent_ref,
	    _ged_draw_prune_groups_child_cb, &ctx);

    for (size_t i = 0; i < snap.len; i++) {
	bsg_scene_ref group_ref = snap.refs[i];
	pruned += ged_draw_scene_ref_prune_empty_groups(gedp, group_ref,
		view_ctx, mode);
	struct ged_draw_scene_tree_summary group_summary;
	if (ged_draw_scene_ref_tree_summary(group_ref, &group_summary) &&
		group_summary.child_count == 0 &&
		ged_draw_scene_ref_group_matches_optional_scope(gedp, group_ref,
		    view_ctx, mode)) {
	    ged_draw_scene_ref_free_group(group_ref);
	    pruned++;
	}
    }

    _ged_draw_clear_scope_snapshot_free(&snap);
    return pruned;
}


static int
ged_draw_scene_ref_visible(bsg_scene_ref ref)
{
    return !bsg_scene_ref_is_null(ref) && bsg_scene_visible(ref);
}


static int
ged_draw_scene_ref_highlighted(bsg_scene_ref ref)
{
    return !bsg_scene_ref_is_null(ref) && bsg_scene_highlighted(ref);
}


static int
ged_draw_scene_ref_legacy_eval_flag(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? 0 : bsg_scene_legacy_eval_flag(ref);
}


static int
ged_draw_scene_ref_legacy_region_id(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? 0 : bsg_scene_legacy_region_id(ref);
}


static int
ged_draw_scene_ref_legacy_basecolor(bsg_scene_ref ref, unsigned char rgb[3])
{
    if (!rgb)
	return 0;
    rgb[0] = 0;
    rgb[1] = 0;
    rgb[2] = 0;
    if (bsg_scene_ref_is_null(ref))
	return 0;
    bsg_scene_legacy_basecolor(ref, rgb);
    return 1;
}


static int
ged_draw_scene_ref_legacy_user_color(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? 0 : bsg_scene_legacy_user_color(ref);
}


static int
ged_draw_scene_ref_legacy_default_color(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? 0 :
	bsg_scene_legacy_default_color(ref);
}


static int
ged_draw_database_source_record_summary_has_state(
	const struct ged_draw_database_source_record *record)
{
    if (!record)
	return 0;
    return (record->database_path && record->database_path[0]) ||
	record->source_revision != 0 ||
	record->inputs_revision != 0 ||
	record->realized_source_revision != 0 ||
	record->realized_inputs_revision != 0 ||
	record->stale_reason != GED_DRAW_STALE_NONE ||
	record->realization_identity != 0;
}


static int
ged_draw_database_source_record_summary_is_stale(
	const struct ged_draw_database_source_record *record)
{
    if (!record)
	return 0;
    return record->stale_reason != GED_DRAW_STALE_NONE ||
	record->source_revision != record->realized_source_revision ||
	record->inputs_revision != record->realized_inputs_revision;
}


static uint64_t
ged_draw_scene_ref_drawn_revision(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? 0 : bsg_scene_drawn_rev(ref);
}


static int
ged_draw_scene_ref_changed(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? 0 : bsg_scene_changed(ref);
}


static fastf_t
ged_draw_scene_ref_transparency(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? 0.0 : bsg_scene_transparency(ref);
}


static int
ged_draw_scene_ref_line_width(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? 0 : bsg_scene_line_width(ref);
}


static fastf_t
ged_draw_scene_ref_draw_size(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? 0.0 : bsg_scene_draw_size(ref);
}


static void
ged_draw_scene_ref_draw_center(bsg_scene_ref ref, point_t center)
{
    if (bsg_scene_ref_is_null(ref) || !center)
	return;

    bsg_scene_draw_center(ref, center);
}


static int
ged_draw_scene_ref_color(bsg_scene_ref ref, unsigned char rgb[3])
{
    if (bsg_scene_ref_is_null(ref) || !rgb)
	return 0;

    bsg_scene_color(ref, &rgb[0], &rgb[1], &rgb[2]);
    return 1;
}


static int
ged_draw_scene_ref_transform(bsg_scene_ref ref, mat_t mat)
{
    if (bsg_scene_ref_is_null(ref) || !mat)
	return 0;

    bsg_scene_transform(ref, mat);
    return 1;
}


static int
ged_draw_scene_ref_is_view_scope(bsg_scene_ref ref)
{
    return !bsg_scene_ref_is_null(ref) && bsg_scene_is_view_scope(ref);
}


static int
ged_draw_scene_ref_is_view_source(bsg_scene_ref ref)
{
    return !bsg_scene_ref_is_null(ref) && bsg_scene_is_view_source(ref);
}


static int
ged_draw_scene_ref_is_local_source(bsg_scene_ref ref)
{
    return !bsg_scene_ref_is_null(ref) && bsg_scene_is_local_source(ref);
}


static int
ged_draw_scene_ref_display_summary(
	bsg_scene_ref ref,
	struct ged_draw_scene_display_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (bsg_scene_ref_is_null(ref))
	return 0;

    struct _ged_draw_obol_display_summary_ctx ctx = {out};
    if (ged_draw_scene_ref_obol_source_path_apply(ref,
	    _ged_draw_obol_display_summary_cb, &ctx))
	return 1;

    memset(out, 0, sizeof(*out));
    const struct bsg_draw_intent *intent = bsg_scene_draw_intent(ref);
    out->valid = 1;
    out->is_database_source = ged_draw_scene_ref_is_database_source(ref);
    out->has_draw_intent = intent ? 1 : 0;
    out->intent_path = intent ? bsg_draw_intent_path(intent) : NULL;
    out->intent_draw_mode = intent ? bsg_draw_intent_mode(intent) : -1;
    out->visible = ged_draw_scene_ref_visible(ref);
    out->highlighted = ged_draw_scene_ref_highlighted(ref);
    out->line_style = ged_draw_scene_ref_line_style(ref);
    out->line_width = ged_draw_scene_ref_line_width(ref);
    out->transparency = ged_draw_scene_ref_transparency(ref);
    out->draw_mode = ged_draw_scene_ref_draw_mode(ref);
    bsg_scene_material_get_rgb(ref,
	    &out->material_color[0],
	    &out->material_color[1],
	    &out->material_color[2]);
    out->material_valid = 1;
    return 1;
}


int
ged_draw_scene_context_display_summary(void *scene_ctx,
				       struct ged_draw_scene_display_summary *out)
{
    return ged_draw_scene_rt_ref_display_summary(
	    ged_draw_scene_context_rt_ref(scene_ctx), out);
}


static int
ged_draw_scene_rt_ref_display_summary(
	rt_view_scene_ref scene_ref,
	struct ged_draw_scene_display_summary *out)
{
    bsg_scene_ref ref = ged_draw_scene_ref_from_rt_view_ref(scene_ref);
    struct ged *gedp = ged_draw_scene_rt_ref_owner_gedp(scene_ref);
    const char *path = ged_draw_scene_rt_ref_semantic_path(scene_ref);
    if (gedp && path) {
	struct ged_draw_obol_context_token *token =
	    ged_draw_obol_context_from_rt_ref(scene_ref);
	struct ged_draw_scene_tree_summary tree_summary;
	int is_group = token ?
	    (token->is_group && !token->is_database_source) : 0;
	if (!token && ged_draw_scene_ref_tree_summary(ref, &tree_summary))
	    is_group = tree_summary.is_group;
	if (is_group) {
	    if (ged_draw_obol_group_display_summary_for_path(gedp, path, out))
		return 1;
	} else if (token && token->is_shape) {
	    if (ged_draw_obol_shape_display_summary_for_path(gedp, path,
		    out))
		return 1;
	    if (ged_draw_obol_database_source_display_summary_for_path(gedp,
		    path, out))
		return 1;
	} else {
	    if (ged_draw_obol_database_source_display_summary_for_path(gedp,
		    path, out))
		return 1;
	}
    }

    return ged_draw_scene_ref_display_summary(ref, out);
}


static int
ged_draw_scene_ref_source_summary(
	bsg_scene_ref ref,
	struct ged_draw_database_source_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (bsg_scene_ref_is_null(ref))
	return 0;

    struct ged_draw_database_source_record record;
    memset(&record, 0, sizeof(record));
    if (!ged_draw_scene_ref_database_source_record_get(ref, &record))
	return 0;

    out->valid = 1;
    out->is_database_source = ged_draw_scene_ref_is_database_source(ref);
    out->has_state =
	ged_draw_database_source_record_summary_has_state(&record);
    out->stale =
	ged_draw_database_source_record_summary_is_stale(&record);
    out->database_path = record.database_path;
    out->source_revision = record.source_revision;
    out->inputs_revision = record.inputs_revision;
    out->realized_source_revision = record.realized_source_revision;
    out->realized_inputs_revision = record.realized_inputs_revision;
    out->realization_identity = record.realization_identity;
    if (record.stale_reason != GED_DRAW_STALE_NONE)
	out->stale_reason = ged_draw_stale_reason_name(record.stale_reason);
    else if (record.source_revision != record.realized_source_revision)
	out->stale_reason =
	    ged_draw_stale_reason_name(GED_DRAW_STALE_SOURCE_CHANGED);
    else if (record.inputs_revision != record.realized_inputs_revision)
	out->stale_reason =
	    ged_draw_stale_reason_name(GED_DRAW_STALE_VIEW_INPUT_CHANGED);
    else
	out->stale_reason = ged_draw_stale_reason_name(GED_DRAW_STALE_NONE);

    return 1;
}


static int
ged_draw_scene_ref_database_source_summary(
	bsg_scene_ref ref,
	struct ged_draw_database_source_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (bsg_scene_ref_is_null(ref))
	return 0;

    return ged_draw_scene_ref_source_summary(
	    ged_draw_scene_ref_source_owner(ref), out);
}


static void *
ged_draw_scene_ref_source_context(bsg_scene_ref ref)
{
    if (bsg_scene_ref_is_null(ref))
	return NULL;

    return ged_draw_scene_ref_context(ged_draw_scene_ref_source_owner(ref));
}


int
ged_draw_shape_context_has_state(void *shape_ctx)
{
    if (!shape_ctx)
	return 0;

    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_ctx(shape_ctx);
    if (token)
	return token->is_database_source || token->is_shape;

    struct ged_draw_shape_record_summary summary;
    return ged_draw_scene_ref_shape_record_summary(
	    ged_draw_scene_ref_from_context(shape_ctx), &summary);
}


void *
ged_draw_shape_context_source(void *shape_ctx)
{
    if (!shape_ctx)
	return NULL;

    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_ctx(shape_ctx);
    if (token) {
	if (token->is_database_source)
	    return shape_ctx;
	if (token->parent_ctx) {
	    struct ged_draw_obol_context_token *parent =
		ged_draw_obol_context_from_scene_ctx(token->parent_ctx);
	    if (parent && parent->is_database_source)
		return token->parent_ctx;
	}
	return NULL;
    }

    return ged_draw_scene_ref_source_context(
	    ged_draw_scene_ref_from_context(shape_ctx));
}


void *
ged_draw_first_shape_context(struct ged *gedp)
{
    return ged_draw_shape_ref_context(gedp, ged_draw_first_shape_ref(gedp));
}


void *
ged_draw_shape_ref_context(struct ged *gedp, ged_draw_shape_ref ref)
{
    void *obol_ctx = _ged_draw_shape_ref_obol_context(gedp, ref);
    if (obol_ctx)
	return obol_ctx;

    return rt_view_scene_ref_context(
	    ged_draw_registry_shape_ref_rt_ref(gedp, ref));
}


void *
ged_draw_shape_ref_cache_context(struct ged *gedp, ged_draw_shape_ref ref)
{
    return rt_view_scene_ref_context(
	    ged_draw_registry_shape_ref_cache_rt_ref(gedp, ref));
}


int
ged_draw_shape_ref_record_summary(struct ged *gedp,
				  ged_draw_shape_ref ref,
				  struct ged_draw_shape_record_summary *out)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_shape_ref_obol_token(gedp, ref);
    if (token && token->is_database_source && token->path &&
	    token->path[0] && out) {
	memset(out, 0, sizeof(*out));
	int obol_summary_scope = token->draw_mode_valid ?
	    ged_draw_obol_database_source_publication_begin(gedp, NULL,
		    token->draw_mode) : 0;
	out->fullpath = token->fullpath_valid ? &token->fullpath : NULL;
	out->display_name = token->path;
	out->leaf_name = token->name;
	if (token->fullpath_valid && token->fullpath.fp_len > 0) {
	    unsigned long long *components =
		(unsigned long long *)bu_calloc(token->fullpath.fp_len,
			sizeof(unsigned long long),
			"GED Obol shape record path hash components");
	    for (size_t i = 0; i < token->fullpath.fp_len; i++) {
		const struct directory *dp = token->fullpath.fp_names[i];
		const char *name = (dp && dp->d_namep) ? dp->d_namep : "";
		components[i] = bu_data_hash(name,
			strlen(name) * sizeof(char));
	    }
	    out->path_hash = bu_data_hash(components,
		    token->fullpath.fp_len * sizeof(unsigned long long));
	    bu_free(components, "GED Obol shape record path hash components");
	}

	out->owning_group_ref = token->fullpath_valid ?
	    ged_draw_obol_top_group_ref_for_fullpath(gedp, &token->fullpath) :
	    ged_draw_obol_top_group_ref_for_path(gedp, token->path);
	if (ged_draw_group_ref_is_null(out->owning_group_ref)) {
	    struct bu_vls owner_group_path = BU_VLS_INIT_ZERO;
	    if (ged_draw_obol_database_source_owner_group_path_for_path(gedp,
		    token->path, &owner_group_path)) {
		out->owning_group_ref = ged_draw_obol_group_ref_for_path(gedp,
			bu_vls_cstr(&owner_group_path));
	    }
	    bu_vls_free(&owner_group_path);
	}

	struct ged_draw_database_source_summary source_summary;
	if (ged_draw_obol_database_source_summary_for_path(gedp,
		token->path, &source_summary) && source_summary.has_state) {
	    out->source_revision = source_summary.source_revision;
	    out->inputs_revision = source_summary.inputs_revision;
	    out->realized_source_revision =
		source_summary.realized_source_revision;
	    out->realized_inputs_revision =
		source_summary.realized_inputs_revision;
	    out->stale = source_summary.stale;
	    out->stale_reason = source_summary.stale_reason;
	}
	if (!out->stale_reason)
	    out->stale_reason = ged_draw_stale_reason_name(
		    GED_DRAW_STALE_NONE);
	out->visible = 1;
	out->draw_mode = token->draw_mode_valid ? token->draw_mode :
	    GED_DRAW_MODE_WIRE;
	out->line_width = 1;

	struct ged_draw_scene_display_summary display_summary;
	if (ged_draw_obol_database_source_display_summary_for_path(gedp,
		token->path, &display_summary) && display_summary.valid) {
	    out->visible = display_summary.visible;
	    out->highlighted = display_summary.highlighted;
	    out->transparency = display_summary.transparency;
	    out->draw_mode = display_summary.draw_mode;
	    out->line_width = display_summary.line_width;
	}

	int evaluated_region = 0;
	if (ged_draw_obol_database_source_evaluated_region_for_path(gedp,
		token->path, &evaluated_region))
	    out->evaluated_region = evaluated_region;

	vect_t bounds_min, bounds_max;
	int empty = 1;
	if (ged_draw_obol_database_source_bounds_for_path(gedp, token->path,
		&bounds_min, &bounds_max, 1, &empty) && !empty) {
	    VADD2SCALE(out->center, bounds_min, bounds_max, 0.5);
	}
	if (obol_summary_scope)
	    ged_draw_obol_database_source_publication_end(gedp);
	return 1;
    }

    if (token && token->is_shape && !token->is_database_source &&
	    token->path && token->path[0] && out) {
	memset(out, 0, sizeof(*out));
	out->fullpath = token->fullpath_valid ? &token->fullpath : NULL;
	out->display_name = token->path;
	out->leaf_name = token->name;
	out->path_hash = bu_data_hash(token->path,
		strlen(token->path) * sizeof(char));

	out->owning_group_ref =
	    ged_draw_obol_shape_owner_group_ref(gedp, token);

	if (!out->stale_reason)
	    out->stale_reason = ged_draw_stale_reason_name(
		    GED_DRAW_STALE_NONE);
	out->visible = 1;
	out->draw_mode = GED_DRAW_MODE_WIRE;
	out->line_width = 1;

	struct ged_draw_scene_display_summary display_summary;
	if (ged_draw_obol_shape_display_summary_for_path(gedp, token->path,
		&display_summary) && display_summary.valid) {
	    out->visible = display_summary.visible;
	    out->highlighted = display_summary.highlighted;
	    out->transparency = display_summary.transparency;
	    out->draw_mode = display_summary.draw_mode;
	    out->line_width = display_summary.line_width;
	}
	out->drawn_revision = ged_draw_scene_revision(gedp);

	struct ged_draw_view_line_summary line_summary;
	if (ged_draw_obol_shape_line_summary_for_path(gedp, token->path,
		&line_summary) && line_summary.valid &&
		line_summary.point_count > 0) {
	    point_t first = VINIT_ZERO;
	    point_t last = VINIT_ZERO;
	    if (ged_draw_obol_shape_line_point_at_for_path(gedp, token->path,
		    0, first) &&
		    ged_draw_obol_shape_line_point_at_for_path(gedp,
			token->path, line_summary.point_count - 1, last))
		VADD2SCALE(out->center, first, last, 0.5);
	}
	return 1;
    }

    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    if (!ged_draw_scene_ref_shape_record_summary(shape_ref, out))
	return 0;

    int obol_evaluated_region = 0;
    struct _ged_draw_obol_evaluated_region_ctx ctx = {
	&obol_evaluated_region
    };
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_evaluated_region_cb, &ctx,
	    "GED Obol evaluated-region read path"))
	out->evaluated_region = obol_evaluated_region;
    return 1;
}


int
ged_draw_shape_ref_display_summary(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   struct ged_draw_scene_display_summary *out)
{
    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    struct _ged_draw_obol_display_summary_ctx ctx = {out};
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_display_summary_cb, &ctx,
	    "GED Obol display summary path"))
	return 1;

    return ged_draw_scene_ref_display_summary(shape_ref, out);
}


int
ged_draw_shape_ref_owning_group_record_summary(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	struct ged_draw_group_record_summary *out)
{
    struct ged_draw_shape_record_summary shape_summary;
    if (!ged_draw_shape_ref_record_summary(gedp, ref, &shape_summary))
	return 0;
    return ged_draw_group_ref_record_summary(
	    gedp, shape_summary.owning_group_ref, out);
}


int
ged_draw_shape_ref_database_source_summary(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	struct ged_draw_database_source_summary *out)
{
    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    struct _ged_draw_obol_source_summary_ctx ctx = {out};
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_source_summary_cb, &ctx,
	    "GED Obol database source summary path"))
	return 1;

    return ged_draw_scene_ref_database_source_summary(shape_ref, out);
}


ged_draw_group_ref
ged_draw_group_ref_of_shape(struct ged *gedp, ged_draw_shape_ref ref)
{
    struct ged_draw_shape_record_summary shape_summary;
    if (!ged_draw_shape_ref_record_summary(gedp, ref, &shape_summary))
	return GED_DRAW_GROUP_REF_NULL;
    return shape_summary.owning_group_ref;
}


void *
ged_draw_group_ref_context(struct ged *gedp, ged_draw_group_ref ref)
{
    void *obol_ctx = ged_draw_group_ref_obol_context(gedp, ref);
    if (obol_ctx)
	return obol_ctx;

    return rt_view_scene_ref_context(
	    ged_draw_registry_group_ref_rt_ref(gedp, ref));
}


int
ged_draw_group_ref_record_summary(struct ged *gedp,
				  ged_draw_group_ref ref,
				  struct ged_draw_group_record_summary *out)
{
    bsg_scene_ref group_ref = _ged_draw_group_ref_scene_ref(gedp, ref);
    struct ged_draw_obol_context_token *token =
	ged_draw_group_ref_obol_token(gedp, ref);
    if ((token || ged_draw_scene_ref_is_null(group_ref)) && out) {
	memset(out, 0, sizeof(*out));
	struct _ged_draw_obol_group_record_ctx ctx = {out};
	if (_ged_draw_group_ref_try_obol_paths(gedp, ref, group_ref,
		_ged_draw_obol_group_record_cb, &ctx,
		"GED Obol group record path"))
	    return 1;
    }

    if (!ged_draw_scene_ref_group_record_summary(group_ref, out))
	return 0;

    struct _ged_draw_obol_group_record_ctx ctx = {out};
    (void)_ged_draw_group_ref_try_obol_paths(gedp, ref, group_ref,
	    _ged_draw_obol_group_record_cb, &ctx,
	    "GED Obol group record path");
    return 1;
}


int
ged_draw_group_ref_tree_summary(struct ged *gedp,
				ged_draw_group_ref ref,
				struct ged_draw_scene_tree_summary *out)
{
    if (!out)
	return 0;
    memset(out, 0, sizeof(*out));

    bsg_scene_ref group_ref = _ged_draw_group_ref_scene_ref(gedp, ref);
    struct ged_draw_scene_tree_summary retained_summary;
    memset(&retained_summary, 0, sizeof(retained_summary));
    int have_retained_summary =
	ged_draw_scene_ref_tree_summary(group_ref, &retained_summary);
    struct _ged_draw_obol_group_tree_summary_ctx ctx = {
	out,
	have_retained_summary ? retained_summary.fullpath : NULL
    };
    if (_ged_draw_group_ref_try_obol_paths(gedp, ref, group_ref,
	    _ged_draw_obol_group_tree_summary_cb, &ctx,
	    "GED Obol group tree-summary path"))
	return 1;

    if (!have_retained_summary)
	return 0;

    *out = retained_summary;
    size_t obol_child_count = 0;
    struct _ged_draw_obol_group_child_count_ctx count_ctx = {
	&obol_child_count
    };
    if (_ged_draw_group_ref_try_obol_paths(gedp, ref, group_ref,
	    _ged_draw_obol_group_child_count_cb, &count_ctx,
	    "GED Obol group child-count path"))
	out->child_count = obol_child_count;
    return 1;
}


static int
ged_draw_group_ref_count_shape_cb(bsg_scene_ref UNUSED(ref), void *ud)
{
    int *count = (int *)ud;
    if (count)
	(*count)++;
    return 1;
}


int
ged_draw_group_ref_shape_count(struct ged *gedp,
			       ged_draw_group_ref ref,
			       int *out)
{
    if (out)
	*out = 0;
    if (!out)
	return 0;

    bsg_scene_ref group_ref = _ged_draw_group_ref_scene_ref(gedp, ref);

    struct _ged_draw_obol_group_shape_count_ctx ctx = {out};
    if (_ged_draw_group_ref_try_obol_paths(gedp, ref, group_ref,
	    _ged_draw_obol_group_shape_count_cb, &ctx,
	    "GED Obol group shape-count path"))
	return 1;

    if (ged_draw_scene_ref_is_null(group_ref))
	return 0;

    ged_draw_source_scene_ref_foreach_shape(group_ref,
	    ged_draw_group_ref_count_shape_cb, out);
    return 1;
}


static bsg_scene_ref
ged_draw_group_ref_find_or_create_child_group(struct ged *gedp,
					      bsg_scene_ref parent_ref,
					      const char *comp_name)
{
    if (!gedp || ged_draw_scene_ref_is_null(parent_ref) || !comp_name)
	return ged_draw_scene_ref_null();

    struct directory *dp = RT_DIR_NULL;
    if (gedp->dbip)
	dp = db_lookup(gedp->dbip, comp_name, LOOKUP_QUIET);

    return ged_draw_view_context_group_child_ensure(
	    ged_draw_active_view_ctx(gedp), parent_ref, comp_name,
	    (void *)dp);
}


static bsg_scene_ref
ged_draw_group_ref_add_path(struct ged *gedp, const char *name)
{
    if (!gedp || !name || !gedp->dbip)
	return ged_draw_scene_ref_null();

    void *view_ctx = ged_draw_active_view_ctx(gedp);
    bsg_scene_ref base_ref =
	_ged_draw_active_scope_scene_ref(gedp, view_ctx, 1, 1);
    if (ged_draw_scene_ref_is_null(base_ref))
	return ged_draw_scene_ref_null();

    struct db_full_path pathcomp;
    if (db_string_to_path(&pathcomp, gedp->dbip, name) != 0) {
	const char *cp = strrchr(name, '/');
	cp = cp ? cp + 1 : name;
	struct directory *dp = db_lookup(gedp->dbip, cp, LOOKUP_NOISY);
	if (dp == RT_DIR_NULL)
	    return ged_draw_scene_ref_null();
	return ged_draw_group_ref_find_or_create_child_group(gedp, base_ref,
		cp);
    }

    if (pathcomp.fp_len == 0) {
	db_free_full_path(&pathcomp);
	return ged_draw_scene_ref_null();
    }

    bsg_scene_ref cur_ref = base_ref;
    for (size_t i = 0; i < pathcomp.fp_len; i++) {
	const char *comp = pathcomp.fp_names[i]->d_namep;
	bsg_scene_ref child_ref =
	    ged_draw_group_ref_find_or_create_child_group(gedp, cur_ref, comp);
	if (ged_draw_scene_ref_is_null(child_ref)) {
	    db_free_full_path(&pathcomp);
	    return ged_draw_scene_ref_null();
	}
	cur_ref = child_ref;
    }

    db_free_full_path(&pathcomp);
    return cur_ref;
}


static bsg_scene_ref
ged_draw_group_ref_lookup_or_create_scene_ref(
	struct ged *gedp,
	const struct db_full_path *dfp)
{
    if (!gedp || !dfp)
	return ged_draw_scene_ref_null();

    char *s = db_path_to_string(dfp);
    if (!s)
	return ged_draw_scene_ref_null();

    const char *path = ged_draw_dbpath_skip_lead_slash(s);

    bsg_scene_ref group_ref = ged_draw_group_ref_add_path(gedp, path);
    if (!ged_draw_scene_ref_is_null(group_ref))
	(void)ged_draw_scene_ref_ensure_draw_intent(group_ref, path,
		GED_DRAW_MODE_WIRE);

    if (!ged_draw_scene_ref_is_null(group_ref))
	(void)ged_draw_scene_ref_apply_path_state(gedp, group_ref, dfp);

    bu_free(s, "draw group path string");
    return group_ref;
}


static int
ged_draw_obol_scene_controller_is_owned(struct ged *gedp)
{
    return (gedp && gedp->i && gedp->i->ged_gdp &&
	    gedp->i->ged_gdp->gd_obol_scene_controller &&
	    (gedp->i->ged_gdp->gd_obol_scene_controller_owned ||
	     gedp->i->ged_gdp->gd_obol_controller_owned) &&
	    gedp->i->ged_gdp->gd_obol_scene_controller_full_sync) ? 1 : 0;
}


ged_draw_group_ref
ged_draw_group_ref_lookup_or_create(struct ged *gedp,
				    const struct db_full_path *dfp)
{
    if (!gedp || !dfp)
	return GED_DRAW_GROUP_REF_NULL;

    char *s = db_path_to_string(dfp);
    if (!s)
	return GED_DRAW_GROUP_REF_NULL;

    const char *path = ged_draw_dbpath_skip_lead_slash(s);
    if (ged_draw_obol_scene_controller_is_owned(gedp)) {
	if (!ged_draw_obol_group_ensure_for_path(gedp, path, path,
		GED_DRAW_MODE_WIRE, 0)) {
	    bu_free(s, "draw group path string");
	    return GED_DRAW_GROUP_REF_NULL;
	}

	ged_draw_group_ref obol_ref = ged_draw_obol_group_ref_for_path(gedp,
		path);
	if (!ged_draw_group_ref_is_null(obol_ref))
	    (void)ged_draw_registry_group_ref_set_indexed_fullpath(gedp,
		    obol_ref, dfp);
	bu_free(s, "draw group path string");
	return obol_ref;
    }

    bu_free(s, "draw group path string");

    bsg_scene_ref group_ref =
	ged_draw_group_ref_lookup_or_create_scene_ref(gedp, dfp);
    if (ged_draw_scene_ref_is_null(group_ref))
	return GED_DRAW_GROUP_REF_NULL;

    return _ged_draw_group_ref_from_scene_ref(gedp, group_ref);
}


int
ged_draw_group_ref_set_dbpath(struct ged *gedp,
			      ged_draw_group_ref ref,
			      const struct db_full_path *new_dfp)
{
    if (!gedp || !new_dfp)
	return 0;

    char *s = db_path_to_string(new_dfp);
    if (!s)
	return 0;

    const char *path = ged_draw_dbpath_skip_lead_slash(s);
    struct ged_draw_group_record_summary group_summary;
    const char *old_path = ged_draw_registry_group_ref_semantic_path(gedp,
	    ref);
    if (ged_draw_group_ref_record_summary(gedp, ref, &group_summary)) {
	if (!old_path)
	    old_path = group_summary.path;
    }

    if (!old_path) {
	struct ged_draw_obol_context_token *token =
	    ged_draw_group_ref_obol_token(gedp, ref);
	if (token && token->path)
	    old_path = token->path;
    }

    int obol_renamed = ged_draw_obol_group_rename_for_path(gedp,
	    old_path ? old_path : path, path);
    int obol_updated = obol_renamed;
    if (!obol_renamed)
	obol_updated = ged_draw_obol_group_update_draw_intent_for_path(gedp,
	    old_path ? old_path : path, path, 0, 0, 0, 0);

    if (obol_updated) {
	(void)ged_draw_registry_group_ref_set_indexed_fullpath(gedp, ref,
		new_dfp);
	bu_free(s, "ged_draw_group_ref_set_dbpath: path string");
	return 1;
    }

    bu_free(s, "ged_draw_group_ref_set_dbpath: path string");
    /* Obol owns public group DB-path mutation; retained fallback is retired. */
    return 0;
}


int
ged_draw_group_ref_set_mode(struct ged *gedp,
			    ged_draw_group_ref ref,
			    int mode)
{
    bsg_scene_ref group_ref = _ged_draw_group_ref_scene_ref(gedp, ref);
    struct _ged_draw_obol_group_draw_intent_ctx ctx = {
	NULL,
	1,
	mode,
	0,
	0
    };
    return _ged_draw_group_ref_try_obol_paths_ensure(gedp, ref,
	    group_ref, _ged_draw_obol_group_draw_intent_cb, &ctx,
	    "GED Obol group mode path", NULL, mode, 0);
}


int
ged_draw_group_ref_set_appearance_settings(struct ged *gedp,
					   ged_draw_group_ref ref,
					   const struct ged_draw_appearance_settings *settings)
{
    if (!settings)
	return 0;

    bsg_scene_ref group_ref = _ged_draw_group_ref_scene_ref(gedp, ref);
    struct _ged_draw_obol_group_appearance_ctx ctx = {settings};
    return _ged_draw_group_ref_try_obol_paths_ensure(gedp, ref,
	    group_ref, _ged_draw_obol_group_update_appearance_cb, &ctx,
	    "GED Obol group appearance path", NULL, settings->draw_mode, 0);
}


int
ged_draw_group_ref_appearance_settings(struct ged *gedp,
				       ged_draw_group_ref ref,
				       struct ged_draw_appearance_settings *settings)
{
    if (!settings)
	return 0;

    bsg_scene_ref group_ref = _ged_draw_group_ref_scene_ref(gedp, ref);

    const char *path = NULL;
    struct ged_draw_group_record_summary group_record;
    if (!ged_draw_scene_ref_is_null(group_ref) &&
	    ged_draw_scene_ref_group_record_summary(group_ref, &group_record))
	path = group_record.path;
    if (!path) {
	struct ged_draw_scene_tree_summary group_summary;
	if (!ged_draw_scene_ref_is_null(group_ref) &&
		ged_draw_scene_ref_tree_summary(group_ref, &group_summary))
	    path = group_summary.name;
    }
    struct _ged_draw_obol_group_appearance_read_ctx ctx = {settings};
    if (_ged_draw_group_ref_try_obol_paths(gedp, ref, group_ref,
	    _ged_draw_obol_group_appearance_read_cb, &ctx,
	    "GED Obol group appearance read path"))
	return 1;

    if (ged_draw_scene_ref_is_null(group_ref))
	return 0;

    return ged_draw_scene_ref_draw_intent_appearance_settings(group_ref,
	    settings);
}


ged_draw_shape_ref
ged_draw_shape_ref_from_context(struct ged *gedp, void *shape_ctx)
{
    return ged_draw_registry_shape_ref_from_source_ref(gedp,
	    ged_draw_scene_context_rt_ref(shape_ctx));
}


static ged_draw_shape_state *
ged_draw_scene_rt_ref_registry_state(rt_view_scene_ref scene_ref)
{
    if (rt_view_scene_ref_is_null(scene_ref))
	return NULL;
    if (ged_draw_obol_context_from_rt_ref(scene_ref))
	return NULL;
    return ged_draw_shape_state_get_scene_rt_ref(scene_ref);
}


static struct ged *
ged_draw_scene_rt_ref_owner_gedp(rt_view_scene_ref scene_ref)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_rt_ref(scene_ref);
    if (token)
	return token->gedp;

    ged_draw_shape_state *state =
	ged_draw_scene_rt_ref_registry_state(scene_ref);
    if (state)
	return state->gedp;

    bsg_scene_ref ref = ged_draw_scene_ref_from_rt_view_ref(scene_ref);
    struct bsg_draw_ctx *draw_ctx = bsg_scene_draw_ctx(ref);
    return draw_ctx ? (struct ged *)draw_ctx->owner_data : NULL;
}


static int
ged_draw_scene_rt_ref_is_source_root(struct ged *gedp,
				     rt_view_scene_ref scene_ref)
{
    if (!gedp || rt_view_scene_ref_backend(scene_ref) !=
	    RT_VIEW_SCENE_BACKEND_BSG)
	return 0;

    bsg_scene_ref ref = ged_draw_scene_ref_from_rt_view_ref(scene_ref);
    bsg_scene_ref root_ref = _ged_draw_source_root_scene_ref(gedp);
    if (ged_draw_scene_ref_is_null(ref))
	return 0;
    if (!ged_draw_scene_ref_is_null(root_ref) &&
	    ged_draw_scene_ref_equal(ref, root_ref))
	return 1;

    struct bsg_draw_ctx *draw_ctx = bsg_scene_draw_ctx(ref);
    return draw_ctx && draw_ctx->owner_data == gedp &&
	ged_draw_scene_ref_is_null(ged_draw_scene_ref_parent(ref));
}


static const char *
ged_draw_scene_rt_ref_semantic_path(rt_view_scene_ref scene_ref)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_rt_ref(scene_ref);
    if (token)
	return (token->path && token->path[0]) ? token->path : NULL;

    ged_draw_shape_state *state =
	ged_draw_scene_rt_ref_registry_state(scene_ref);
    return (state && state->display_name && state->display_name[0]) ?
	state->display_name : NULL;
}


static int
ged_draw_shape_context_obol_path_apply(
	void *shape_ctx,
	ged_draw_shape_ref_obol_path_cb cb,
	void *userdata,
	const char *free_label)
{
    if (!shape_ctx || !cb)
	return 0;

    rt_view_scene_ref scene_ref = ged_draw_scene_context_rt_ref(shape_ctx);
    struct ged *gedp = ged_draw_scene_rt_ref_owner_gedp(scene_ref);
    if (!gedp)
	return 0;

    bsg_scene_ref bsg_ref = ged_draw_scene_ref_from_rt_view_ref(scene_ref);
    if (!ged_draw_scene_ref_is_null(bsg_ref)) {
	ged_draw_shape_ref shape_ref =
	    _ged_draw_shape_ref_from_scene_ref(gedp, bsg_ref);
	if (_ged_draw_shape_ref_try_obol_paths(gedp, shape_ref, bsg_ref, cb,
		userdata, free_label))
	    return 1;
    }

    return _ged_draw_shape_ref_try_obol_path(gedp,
	    ged_draw_scene_rt_ref_semantic_path(scene_ref), cb, userdata);
}


const struct db_full_path *
ged_draw_scene_context_fullpath(void *scene_ctx)
{
    return ged_draw_scene_rt_ref_fullpath(
	    ged_draw_scene_context_rt_ref(scene_ctx));
}


static const struct db_full_path *
ged_draw_scene_rt_ref_fullpath(rt_view_scene_ref scene_ref)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_rt_ref(scene_ref);
    if (token)
	return token->fullpath_valid ? &token->fullpath : NULL;

    ged_draw_shape_state *state =
	ged_draw_scene_rt_ref_registry_state(scene_ref);
    if (state && state->s_fullpath.fp_len > 0)
	return &state->s_fullpath;

    struct ged_draw_scene_tree_summary summary;
    if (!ged_draw_scene_ref_tree_summary(
		ged_draw_scene_ref_from_rt_view_ref(scene_ref), &summary))
	return NULL;
    return summary.fullpath;
}


int
ged_draw_group_context_dbpath(struct ged *gedp,
			      void *group_ctx,
			      struct db_full_path *out)
{
    (void)gedp;
    if (!group_ctx || !out)
	return -1;

    if (ged_draw_group_context_is_overlay(group_ctx))
	return -1;

    if (ged_draw_obol_context_from_scene_ctx(group_ctx)) {
	const struct db_full_path *stored =
	    ged_draw_scene_context_fullpath(group_ctx);
	if (!stored || stored->fp_len <= 0)
	    return -1;
	db_dup_full_path(out, stored);
	return 0;
    }

    bsg_scene_ref group_ref = ged_draw_scene_ref_from_context(group_ctx);
    struct ged_draw_group_record_summary group_summary;
    if (ged_draw_scene_ref_group_record_summary(group_ref, &group_summary) &&
	    group_summary.is_overlay)
	return -1;

    const struct db_full_path *stored = ged_draw_scene_context_fullpath(group_ctx);
    if (!stored || stored->fp_len <= 0)
	return -1;
    db_dup_full_path(out, stored);
    return 0;
}


int
ged_draw_group_context_is_overlay(void *group_ctx)
{
    if (!group_ctx)
	return 0;

    rt_view_scene_ref scene_ref = ged_draw_scene_context_rt_ref(group_ctx);
    struct ged *gedp = ged_draw_scene_rt_ref_owner_gedp(scene_ref);
    const char *path = ged_draw_scene_rt_ref_semantic_path(scene_ref);
    if (gedp && path) {
	struct ged_draw_group_record_summary obol_summary;
	memset(&obol_summary, 0, sizeof(obol_summary));
	if (ged_draw_obol_group_record_summary_for_path(gedp, path,
		&obol_summary))
	    return obol_summary.is_overlay;
    }

    struct ged_draw_group_record_summary group_summary;
    if (!ged_draw_scene_ref_group_record_summary(
	    ged_draw_scene_ref_from_context(group_ctx), &group_summary))
	return 0;
    return group_summary.is_overlay;
}


int
ged_draw_scene_context_source_summary(void *scene_ctx,
				      struct ged_draw_database_source_summary *out)
{
    return ged_draw_scene_rt_ref_source_summary(
	    ged_draw_scene_context_rt_ref(scene_ctx), out);
}


static int
ged_draw_scene_rt_ref_source_summary(
	rt_view_scene_ref scene_ref,
	struct ged_draw_database_source_summary *out)
{
    struct ged *gedp = ged_draw_scene_rt_ref_owner_gedp(scene_ref);
    const char *path = ged_draw_scene_rt_ref_semantic_path(scene_ref);
    if (gedp && path &&
	    ged_draw_obol_database_source_summary_for_path(gedp, path, out))
	return 1;

    return ged_draw_scene_ref_database_source_summary(
	    ged_draw_scene_ref_from_rt_view_ref(scene_ref), out);
}


int
ged_draw_scene_context_tree_summary(void *scene_ctx,
				    struct ged_draw_scene_tree_summary *out)
{
    return ged_draw_scene_rt_ref_tree_summary(
	    ged_draw_scene_context_rt_ref(scene_ctx), out);
}


static int
ged_draw_scene_rt_ref_tree_summary(
	rt_view_scene_ref scene_ref,
	struct ged_draw_scene_tree_summary *out)
{
    if (!out)
	return 0;

    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_rt_ref(scene_ref);
    if (token)
	return ged_draw_obol_context_token_tree_summary(token, NULL, out);

    bsg_scene_ref ref = ged_draw_scene_ref_from_rt_view_ref(scene_ref);
    struct ged_draw_scene_tree_summary retained_summary;
    memset(&retained_summary, 0, sizeof(retained_summary));
    int have_retained_summary =
	ged_draw_scene_ref_tree_summary(ref, &retained_summary);

    ged_draw_shape_state *state =
	ged_draw_scene_rt_ref_registry_state(scene_ref);
    const struct db_full_path *fallback_fullpath =
	have_retained_summary ? retained_summary.fullpath : NULL;
    if (state && state->s_fullpath.fp_len > 0)
	fallback_fullpath = &state->s_fullpath;

    struct ged *gedp = ged_draw_scene_rt_ref_owner_gedp(scene_ref);
    if (gedp && ged_draw_scene_rt_ref_is_source_root(gedp, scene_ref) &&
	    ged_draw_obol_scene_controller_full_synced(gedp)) {
	size_t obol_child_count = 0;
	if (ged_draw_obol_scene_root_child_count(gedp, &obol_child_count)) {
	    if (have_retained_summary)
		*out = retained_summary;
	    else {
		memset(out, 0, sizeof(*out));
		out->valid = 1;
		out->is_group = 1;
		out->is_shape = 0;
		out->has_parent = 0;
		out->name = "_draw_root";
	    }
	    out->child_count = obol_child_count;
	    return 1;
	}
    }

    const char *path = ged_draw_scene_rt_ref_semantic_path(scene_ref);
    if (gedp && path) {
	void *obol_ctx = ged_draw_obol_context_token_with_parent_for_path(
		gedp, path);
	token = ged_draw_obol_context_from_scene_ctx(obol_ctx);
	if (ged_draw_obol_context_token_tree_summary(token,
		fallback_fullpath, out)) {
	    if (out->child_count == 0 && have_retained_summary &&
		    retained_summary.child_count > 0)
		out->child_count = retained_summary.child_count;
	    return 1;
	}
    }

    if (!have_retained_summary)
	return 0;

    *out = retained_summary;
    if (out->is_group) {
	size_t obol_child_count = 0;
	if (gedp && path &&
		ged_draw_obol_group_child_count_for_path(gedp, path,
		    &obol_child_count))
	    out->child_count = obol_child_count;
    }

    if (state && state->s_fullpath.fp_len > 0)
	out->fullpath = &state->s_fullpath;
    return 1;
}


struct ged_draw_context_child_at_ctx {
    size_t index;
    size_t cur;
    rt_view_scene_ref ref;
};


static int
_ged_draw_context_child_at_cb(rt_view_scene_ref child_ref, void *userdata)
{
    struct ged_draw_context_child_at_ctx *ctx =
	(struct ged_draw_context_child_at_ctx *)userdata;
    if (!ctx)
	return 0;
    if (ctx->cur == ctx->index) {
	ctx->ref = child_ref;
	return 0;
    }
    ctx->cur++;
    return 1;
}


static rt_view_scene_ref
ged_draw_scene_rt_ref_child_at(rt_view_scene_ref scene_ref, size_t index)
{
    struct ged_draw_context_child_at_ctx child_ctx;
    child_ctx.index = index;
    child_ctx.cur = 0;
    child_ctx.ref = rt_view_scene_ref_null();

    (void)ged_draw_scene_rt_ref_foreach_child(scene_ref,
	    _ged_draw_context_child_at_cb,
	    &child_ctx);
    return child_ctx.ref;
}


static void *
ged_draw_obol_scene_context_child_at(
	rt_view_scene_ref scene_ref,
	void *parent_ctx,
	size_t index)
{
    struct ged *gedp = ged_draw_scene_rt_ref_owner_gedp(scene_ref);
    struct ged_draw_obol_scene_context_info info;
    memset(&info, 0, sizeof(info));

    if (gedp && ged_draw_scene_rt_ref_is_source_root(gedp, scene_ref) &&
	    ged_draw_obol_scene_controller_full_synced(gedp)) {
	if (index > (size_t)INT_MAX)
	    return NULL;
	if (!ged_draw_obol_scene_child_context_info_for_path(gedp, "/",
		index, &info))
	    return NULL;
	void *child_ctx = ged_draw_obol_context_token_create(gedp, &info,
		parent_ctx);
	ged_draw_obol_scene_context_info_free(&info);
	return child_ctx;
    }

    const char *path = ged_draw_scene_rt_ref_semantic_path(scene_ref);
    if (!gedp || !path)
	return NULL;

    size_t child_index = index;
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_rt_ref(scene_ref);
    bsg_scene_ref ref = ged_draw_scene_ref_from_rt_view_ref(scene_ref);
    int source_context = token ? token->is_database_source :
	ged_draw_scene_ref_is_database_source(ref);
    if (!source_context && token && token->is_group &&
	    !token->is_database_source && index == 0 &&
	    ged_draw_obol_database_source_container_info(gedp, path,
		token->draw_tree_depth + 1, &info)) {
	void *child_ctx = ged_draw_obol_context_token_create(gedp, &info,
		parent_ctx);
	ged_draw_obol_scene_context_info_free(&info);
	return child_ctx;
    }

    int need_primary_proxy = source_context ?
	ged_draw_obol_database_source_needs_primary_shape_proxy(gedp, path) : 0;
    if (need_primary_proxy) {
	struct ged_draw_scene_tree_summary parent_summary;
	memset(&parent_summary, 0, sizeof(parent_summary));
	if (index == 0 &&
		ged_draw_scene_rt_ref_tree_summary(scene_ref,
		    &parent_summary) &&
		ged_draw_obol_database_source_primary_shape_info(gedp, path,
		    parent_summary.draw_tree_depth + 1, &info)) {
	    void *child_ctx = ged_draw_obol_context_token_create(gedp, &info,
		    parent_ctx);
	    ged_draw_obol_scene_context_info_free(&info);
	    return child_ctx;
	}
	child_index--;
    }

    if (!ged_draw_obol_scene_child_context_info_for_path(gedp, path,
	    child_index, &info))
	return NULL;

    void *child_ctx = ged_draw_obol_context_token_create(gedp, &info,
	    parent_ctx);
    ged_draw_obol_scene_context_info_free(&info);
    return child_ctx;
}


void *
ged_draw_scene_context_child_at(void *scene_ctx, size_t index)
{
    rt_view_scene_ref scene_ref = ged_draw_scene_context_rt_ref(scene_ctx);
    void *obol_child = ged_draw_obol_scene_context_child_at(scene_ref,
	    scene_ctx, index);
    if (obol_child)
	return obol_child;

    return rt_view_scene_ref_context(ged_draw_scene_rt_ref_child_at(
		scene_ref, index));
}


static rt_view_scene_ref
ged_draw_scene_rt_ref_parent(rt_view_scene_ref scene_ref)
{
    bsg_scene_ref parent_ref = ged_draw_scene_ref_parent(
	    ged_draw_scene_ref_from_rt_view_ref(scene_ref));
    return ged_draw_scene_ref_to_rt_view_ref(parent_ref);
}


static void *
ged_draw_obol_scene_context_parent(rt_view_scene_ref scene_ref)
{
    struct ged *gedp = ged_draw_scene_rt_ref_owner_gedp(scene_ref);
    const char *path = ged_draw_scene_rt_ref_semantic_path(scene_ref);
    bsg_scene_ref ref = ged_draw_scene_ref_from_rt_view_ref(scene_ref);

    if (!gedp || !path || !ged_draw_scene_ref_is_database_source(ref))
	return NULL;

    return ged_draw_obol_database_source_parent_token_for_path(gedp, path);
}


void *
ged_draw_scene_context_parent(void *scene_ctx)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_ctx(scene_ctx);
    if (token)
	return token->parent_ctx;

    rt_view_scene_ref scene_ref = ged_draw_scene_context_rt_ref(scene_ctx);
    void *obol_parent = ged_draw_obol_scene_context_parent(scene_ref);
    if (obol_parent)
	return obol_parent;

    return rt_view_scene_ref_context(ged_draw_scene_rt_ref_parent(
		scene_ref));
}


const char *
ged_draw_scene_context_name(void *scene_ctx)
{
    return ged_draw_scene_rt_ref_name(ged_draw_scene_context_rt_ref(scene_ctx));
}


static const char *
ged_draw_scene_rt_ref_name(rt_view_scene_ref scene_ref)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_rt_ref(scene_ref);
    if (token)
	return token->name ? token->name : token->path;

    struct ged_draw_scene_tree_summary summary;
    if (ged_draw_scene_ref_tree_summary(
		ged_draw_scene_ref_from_rt_view_ref(scene_ref), &summary) &&
	    summary.name)
	return summary.name;

    ged_draw_shape_state *state =
	ged_draw_scene_rt_ref_registry_state(scene_ref);
    if (state && state->leaf_dp && state->leaf_dp->d_namep)
	return state->leaf_dp->d_namep;
    return state ? state->display_name : NULL;
}


int
ged_draw_scene_context_subtree_bounds(void *scene_ctx,
				      vect_t *min,
				      vect_t *max,
				      int include_overlays)
{
    return ged_draw_scene_rt_ref_subtree_bounds(
	    ged_draw_scene_context_rt_ref(scene_ctx), min, max,
	    include_overlays);
}


static int
ged_draw_scene_rt_ref_subtree_bounds(rt_view_scene_ref scene_ref,
				     vect_t *min,
				     vect_t *max,
				     int include_overlays)
{
    if (!min || !max)
	return 1;

    struct ged *gedp = ged_draw_scene_rt_ref_owner_gedp(scene_ref);
    const char *path = ged_draw_scene_rt_ref_semantic_path(scene_ref);
    if (gedp && path) {
	bsg_scene_ref ref = ged_draw_scene_ref_from_rt_view_ref(scene_ref);
	struct ged_draw_scene_tree_summary tree_summary;
	struct ged_draw_obol_context_token *token =
	    ged_draw_obol_context_from_rt_ref(scene_ref);
	int empty = 1;
	int is_group_context = token ?
	    (token->is_group && !token->is_database_source) : 0;
	if (!token && ged_draw_scene_ref_tree_summary(ref, &tree_summary))
	    is_group_context = tree_summary.is_group;
	if (is_group_context) {
	    if (ged_draw_obol_group_subtree_bounds_for_path(gedp, path,
		    min, max, include_overlays, &empty))
		return empty;
	} else if (ged_draw_obol_database_source_bounds_for_path(gedp,
		path, min, max, include_overlays, &empty)) {
	    return empty;
	}
    }

    return ged_draw_scene_ref_subtree_bounds(
	    ged_draw_scene_ref_from_rt_view_ref(scene_ref), min, max,
	    include_overlays);
}


static int
ged_draw_scene_ref_tree_summary(
	bsg_scene_ref ref,
	struct ged_draw_scene_tree_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (bsg_scene_ref_is_null(ref))
	return 0;

    out->valid = 1;
    out->is_group = ged_draw_scene_ref_is_group(ref);
    out->is_shape = ged_draw_scene_ref_is_shape(ref);
    out->has_parent =
	ged_draw_scene_ref_is_null(ged_draw_scene_ref_parent(ref)) ? 0 : 1;
    out->name = ged_draw_scene_ref_name(ref);
    out->fullpath = ged_draw_scene_ref_fullpath(ref);
    out->draw_tree_depth = ged_draw_scene_ref_draw_tree_depth(ref);
    out->child_count = ged_draw_scene_ref_child_count(ref);
    return 1;
}


static int
ged_draw_scene_ref_material_summary(
	bsg_scene_ref ref,
	struct ged_draw_shape_material_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (bsg_scene_ref_is_null(ref))
	return 0;

    struct _ged_draw_obol_material_summary_ctx ctx = {out};
    if (ged_draw_scene_ref_obol_source_path_apply(ref,
	    _ged_draw_obol_material_summary_cb, &ctx))
	return 1;

    memset(out, 0, sizeof(*out));
    out->valid = 1;
    out->material_revision = (uint64_t)bsg_scene_material_revision(ref);
    bsg_scene_material_get_rgb(ref,
	    &out->material_color[0],
	    &out->material_color[1],
	    &out->material_color[2]);
    return 1;
}


int
ged_draw_shape_ref_material_summary(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    struct ged_draw_shape_material_summary *out)
{
    if (!out)
	return 0;

    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	shape_ref = ged_draw_scene_ref_from_rt_view_ref(
		ged_draw_registry_shape_ref_cache_rt_ref(gedp, ref));
    struct _ged_draw_obol_material_summary_ctx ctx = {out};
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_material_summary_cb, &ctx,
	    "GED Obol material summary path"))
	return 1;

    memset(out, 0, sizeof(*out));
    return ged_draw_scene_ref_material_summary(shape_ref, out);
}


static int
ged_draw_scene_ref_shape_record_summary(
	bsg_scene_ref ref,
	struct ged_draw_shape_record_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!ged_draw_scene_ref_is_shape(ref))
	return 0;

    ged_draw_shape_state *shape_data = _ged_draw_shape_state_get_scene_ref(ref);
    if (!shape_data)
	return 0;

    out->fullpath = &shape_data->s_fullpath;
    out->display_name = shape_data->display_name;
    out->owning_group_ref = _ged_draw_group_ref_from_scene_ref(
	    shape_data->gedp, ged_draw_scene_ref_shape_owning_group(ref));
    out->path_hash = shape_data->path_hash;
    out->source_revision = shape_data->source_revision;
    out->inputs_revision = shape_data->inputs_revision;
    out->realized_source_revision = shape_data->realized_source_revision;
    out->realized_inputs_revision = shape_data->realized_inputs_revision;
    out->stale = (shape_data->stale_reason != GED_DRAW_STALE_NONE ||
	    shape_data->source_revision != shape_data->realized_source_revision ||
	    shape_data->inputs_revision != shape_data->realized_inputs_revision);
    out->visible = ged_draw_scene_ref_visible(ref);
    out->highlighted = ged_draw_scene_ref_highlighted(ref);
    out->evaluated_region = ged_draw_scene_ref_legacy_eval_flag(ref);
    out->drawn_revision = ged_draw_scene_ref_drawn_revision(ref);
    out->transparency = ged_draw_scene_ref_transparency(ref);
    out->draw_mode = ged_draw_scene_ref_draw_mode(ref);
    out->line_width = ged_draw_scene_ref_line_width(ref);
    ged_draw_scene_ref_draw_center(ref, out->center);

    if (shape_data->stale_reason != GED_DRAW_STALE_NONE)
	out->stale_reason = ged_draw_stale_reason_name(shape_data->stale_reason);
    else if (shape_data->source_revision != shape_data->realized_source_revision)
	out->stale_reason =
	    ged_draw_stale_reason_name(GED_DRAW_STALE_SOURCE_CHANGED);
    else if (shape_data->inputs_revision != shape_data->realized_inputs_revision)
	out->stale_reason =
	    ged_draw_stale_reason_name(GED_DRAW_STALE_VIEW_INPUT_CHANGED);
    else
	out->stale_reason = ged_draw_stale_reason_name(GED_DRAW_STALE_NONE);

    if (shape_data->leaf_dp)
	out->leaf_name = shape_data->leaf_dp->d_namep;
    else if (shape_data->s_fullpath.fp_len > 0 &&
	    shape_data->s_fullpath.fp_names[shape_data->s_fullpath.fp_len - 1])
	out->leaf_name =
	    shape_data->s_fullpath.fp_names[shape_data->s_fullpath.fp_len - 1]->d_namep;

    return 1;
}


static int
_ged_draw_scene_ref_in_view_scope(bsg_scene_ref ref)
{
    bsg_scene_ref cur = ref;

    while (!ged_draw_scene_ref_is_null(cur)) {
	if (ged_draw_scene_ref_is_view_scope(cur))
	    return 1;
	cur = ged_draw_scene_ref_parent(cur);
    }

    return 0;
}


static int _ged_draw_scene_ref_first_shape_record_child_cb(
	bsg_scene_ref child_ref,
	void *userdata);


static bsg_scene_ref
_ged_draw_scene_ref_first_shape_record(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return ged_draw_scene_ref_null();

    if (_ged_draw_shape_state_get_scene_ref(ref))
	return ref;

    bsg_scene_ref shape_ref = ged_draw_scene_ref_null();
    (void)ged_draw_scene_ref_foreach_child(ref,
	    _ged_draw_scene_ref_first_shape_record_child_cb, &shape_ref);
    return shape_ref;
}


static int
_ged_draw_scene_ref_first_shape_record_child_cb(bsg_scene_ref child_ref,
						void *userdata)
{
    bsg_scene_ref *shape_ref = (bsg_scene_ref *)userdata;

    if (!shape_ref)
	return 1;
    *shape_ref = _ged_draw_scene_ref_first_shape_record(child_ref);
    return ged_draw_scene_ref_is_null(*shape_ref) ? 1 : 0;
}


static int
ged_draw_scene_ref_group_record_summary(
	bsg_scene_ref ref,
	struct ged_draw_group_record_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (ged_draw_scene_ref_is_null(ref) || !ged_draw_scene_ref_is_group(ref))
	return 0;

    out->path = ged_draw_scene_ref_draw_intent_path(ref);
    if (!ged_draw_scene_ref_draw_intent_mode(ref, &out->draw_mode)) {
	bsg_scene_ref shape_ref = _ged_draw_scene_ref_first_shape_record(ref);
	if (!ged_draw_scene_ref_is_null(shape_ref))
	    out->draw_mode = ged_draw_scene_ref_draw_mode(shape_ref);
	else
	    out->draw_mode = GED_DRAW_MODE_WIRE;
    }

    bsg_scene_ref first_shape = _ged_draw_scene_ref_first_shape_record(ref);
    if (!ged_draw_scene_ref_is_null(first_shape))
	out->transparency = ged_draw_scene_ref_transparency(first_shape);
    out->visible = ged_draw_scene_ref_visible(ref);
    out->is_overlay = ged_draw_scene_ref_draw_intent_is_overlay(ref);
    out->is_view_scope = ged_draw_scene_ref_is_view_scope(ref);
    out->in_view_scope = _ged_draw_scene_ref_in_view_scope(ref);
    out->is_view_source = ged_draw_scene_ref_is_view_source(ref);
    out->is_local_source = ged_draw_scene_ref_is_local_source(ref);
    out->view_ctx = ged_draw_scene_ref_view_context(ref);
    (void)ged_draw_scene_ref_obol_group_record_summary(ref, out);
    return 1;
}


static int
ged_draw_scene_ref_draw_state_summary(
	bsg_scene_ref ref,
	struct ged_draw_scene_draw_state_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (ged_draw_scene_ref_is_null(ref))
	return 0;

    out->name = ged_draw_scene_ref_name(ref);
    out->draw_mode = ged_draw_scene_ref_draw_mode(ref);
    out->line_style = ged_draw_scene_ref_line_style(ref);
    out->pipeline_candidate =
	ged_draw_scene_ref_realization_pipeline_candidate(ref) ? 1 : 0;
    out->draw_mat_valid = ged_draw_scene_ref_draw_mat(ref, out->draw_mat);
    out->bounds_valid = ged_draw_scene_ref_bounds(ref, out->bounds_min,
	    out->bounds_max);
    out->view_ctx = ged_draw_scene_ref_view_context(ref);

    struct ged_draw_obol_draw_state_summary obol_summary;
    memset(&obol_summary, 0, sizeof(obol_summary));
    if (ged_draw_scene_ref_obol_draw_state_summary(ref, &obol_summary) &&
	    obol_summary.valid) {
	if (obol_summary.draw_mode_valid)
	    out->draw_mode = obol_summary.draw_mode;
	out->line_style = obol_summary.line_style;
	if (obol_summary.draw_mat_valid) {
	    out->draw_mat_valid = 1;
	    MAT_COPY(out->draw_mat, obol_summary.draw_mat);
	}
    }

    return 1;
}


static int
ged_draw_scene_ref_source_snapshot(
	bsg_scene_ref ref,
	struct ged_draw_shape_source_snapshot *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    ged_draw_shape_state *shape_data = _ged_draw_shape_state_get_scene_ref(ref);
    struct ged_draw_source_state *source =
	shape_data ? (struct ged_draw_source_state *)shape_data->source_data : NULL;

    struct ged_draw_obol_database_source_runtime obol_runtime;
    memset(&obol_runtime, 0, sizeof(obol_runtime));
    if (shape_data &&
	    ged_draw_scene_ref_obol_database_source_runtime(ref,
		&obol_runtime) &&
	    obol_runtime.valid &&
	    obol_runtime.dbip) {
	if (shape_data->s_fullpath.fp_len <= 0 &&
		obol_runtime.database_path &&
		obol_runtime.database_path[0]) {
	    const char *path = ged_draw_dbpath_skip_lead_slash(
		    obol_runtime.database_path);
	    if (path && path[0] &&
		    db_string_to_path(&shape_data->s_fullpath,
			obol_runtime.dbip, path) == 0)
		shape_data->leaf_dp =
		    DB_FULL_PATH_CUR_DIR(&shape_data->s_fullpath);
	}

	if (source && source->tol)
	    shape_data->obol_snapshot_tol = *source->tol;
	else
	    BN_TOL_INIT_SET_TOL(&shape_data->obol_snapshot_tol);

	if (source && source->ttol)
	    shape_data->obol_snapshot_ttol = *source->ttol;
	else
	    BG_TESS_TOL_INIT_SET_TOL(&shape_data->obol_snapshot_ttol);
	if (obol_runtime.tessellation_abs_tol >= 0.0)
	    shape_data->obol_snapshot_ttol.abs =
		obol_runtime.tessellation_abs_tol;
	if (obol_runtime.tessellation_rel_tol >= 0.0)
	    shape_data->obol_snapshot_ttol.rel =
		obol_runtime.tessellation_rel_tol;
	if (obol_runtime.tessellation_norm_tol >= 0.0)
	    shape_data->obol_snapshot_ttol.norm =
		obol_runtime.tessellation_norm_tol;
	shape_data->obol_snapshot_valid = 1;

	out->dbip = obol_runtime.dbip;
	out->fullpath = shape_data->s_fullpath.fp_len > 0 ?
	    &shape_data->s_fullpath :
	    ged_draw_scene_ref_fullpath(ref);
	out->leaf_dp = shape_data->s_fullpath.fp_len > 0 ?
	    DB_FULL_PATH_CUR_DIR(&shape_data->s_fullpath) :
	    ged_draw_scene_ref_leaf_dp(ref);
	out->name = shape_data->display_name ?
	    shape_data->display_name :
	    (obol_runtime.database_path && obol_runtime.database_path[0] ?
	     obol_runtime.database_path : ged_draw_scene_ref_name(ref));
	out->tol = &shape_data->obol_snapshot_tol;
	out->ttol = &shape_data->obol_snapshot_ttol;
	return 1;
    }

    if (!source)
	return 0;

    out->dbip = source->dbip;
    out->fullpath = ged_draw_scene_ref_fullpath(ref);
    out->leaf_dp = ged_draw_scene_ref_leaf_dp(ref);
    out->name = ged_draw_scene_ref_name(ref);
    out->tol = source->tol;
    out->ttol = source->ttol;
    return 1;
}


static int
ged_draw_shape_ref_obol_clear_vlist(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    bsg_scene_ref shape_ref)
{
    return _ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_clear_vlist_cb, NULL,
	    "GED Obol clear VLIST path");
}


static int
ged_draw_shape_ref_obol_clear_mesh(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   bsg_scene_ref shape_ref)
{
    return _ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_clear_mesh_cb, NULL,
	    "GED Obol clear mesh path");
}


static int
ged_draw_shape_ref_obol_clear_auxiliary_shapes(struct ged *gedp,
					       ged_draw_shape_ref ref,
					       bsg_scene_ref shape_ref)
{
    return _ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_clear_auxiliary_shapes_cb, NULL,
	    "GED Obol clear auxiliary VLIST path");
}


static int
ged_draw_shape_ref_obol_publish_line_set(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	bsg_scene_ref shape_ref,
	const point_t *points,
	const int *commands,
	size_t point_count)
{
    struct _ged_draw_obol_publish_line_ctx ctx = {
	points,
	commands,
	point_count
    };
    (void)ged_draw_shape_ref_obol_sync_source_placement(gedp, ref,
	    shape_ref);
    return _ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_publish_line_cb, &ctx,
	    "GED Obol publish VLIST path");
}


static int
ged_draw_scene_ref_obol_publish_annotation_line_set(
	bsg_scene_ref ref,
	const point_t *points,
	const int *commands,
	size_t point_count)
{
    ged_draw_shape_state *shape_data = _ged_draw_shape_state_get_scene_ref(ref);
    if (!shape_data || !shape_data->gedp)
	return 0;

    ged_draw_shape_ref shape_ref =
	_ged_draw_shape_ref_from_scene_ref(shape_data->gedp, ref);
    (void)_ged_draw_scene_ref_obol_ensure_owner_source(ref);
    struct _ged_draw_obol_publish_line_ctx ctx = {
	points,
	commands,
	point_count
    };
    (void)ged_draw_shape_ref_obol_sync_source_placement(shape_data->gedp,
	    shape_ref, ref);
    return _ged_draw_shape_ref_try_obol_paths(shape_data->gedp, shape_ref, ref,
	    _ged_draw_obol_publish_annotation_line_cb, &ctx,
	    "GED Obol publish annotation VLIST path");
}


static int
ged_draw_scene_ref_obol_publish_auxiliary_line_set(
	bsg_scene_ref ref,
	const char *name,
	const point_t *points,
	const int *commands,
	size_t point_count)
{
    ged_draw_shape_state *shape_data = _ged_draw_shape_state_get_scene_ref(ref);
    if (!shape_data || !shape_data->gedp || !name || !name[0])
	return 0;

    ged_draw_shape_ref shape_ref =
	_ged_draw_shape_ref_from_scene_ref(shape_data->gedp, ref);
    (void)_ged_draw_scene_ref_obol_ensure_owner_source(ref);
    struct ged_draw_scene_display_summary display_state;
    memset(&display_state, 0, sizeof(display_state));
    (void)ged_draw_scene_ref_display_summary(ref, &display_state);
    struct _ged_draw_obol_publish_auxiliary_line_ctx ctx = {
	name,
	points,
	commands,
	point_count,
	display_state.valid ? &display_state : NULL
    };
    (void)ged_draw_shape_ref_obol_sync_source_placement(shape_data->gedp,
	    shape_ref, ref);
    if (_ged_draw_shape_ref_try_obol_paths(shape_data->gedp, shape_ref, ref,
	    _ged_draw_obol_publish_auxiliary_line_cb, &ctx,
	    "GED Obol publish auxiliary VLIST path"))
	return 1;

    return _ged_draw_scene_ref_obol_source_snapshot_path_apply(ref,
	    _ged_draw_obol_publish_auxiliary_line_cb, &ctx);
}


static int
ged_draw_shape_ref_obol_update_vlist_bounds(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	bsg_scene_ref shape_ref)
{
    return _ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_update_bounds_cb, NULL,
	    "GED Obol update VLIST bounds path");
}


int
ged_draw_shape_ref_set_center(struct ged *gedp, ged_draw_shape_ref ref,
			      const point_t center)
{
    if (!center)
	return 0;
    struct _ged_draw_obol_center_ctx ctx;
    VMOVE(ctx.center, center);
    /* Obol owns public shape center mutation; retained fallback is retired. */
    return _ged_draw_shape_ref_try_obol_paths_lazy(gedp, ref, NULL,
	    _ged_draw_obol_set_center_cb, &ctx,
	    "GED Obol set VLIST center path");
}


int
ged_draw_shape_ref_geometry_clear(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    int obol_cleared = 0;
    obol_cleared = ged_draw_shape_ref_obol_clear_vlist(gedp, ref,
	    shape_ref);
    obol_cleared =
	ged_draw_shape_ref_obol_clear_mesh(gedp, ref, shape_ref) ||
	obol_cleared;
    obol_cleared =
	ged_draw_shape_ref_obol_clear_auxiliary_shapes(gedp, ref, shape_ref) ||
	obol_cleared;
    if (obol_cleared)
	return 1;
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;
    return ged_draw_scene_ref_geometry_clear(shape_ref);
}


int
ged_draw_shape_ref_update_bounds_from_geometry(struct ged *gedp,
					       ged_draw_shape_ref ref,
					       int *bad_cmd)
{
    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    if (bad_cmd)
	*bad_cmd = 0;
    if (ged_draw_shape_ref_obol_update_vlist_bounds(gedp, ref, shape_ref))
	return 1;
    return ged_draw_scene_ref_update_bounds_from_geometry(shape_ref, bad_cmd);
}


int
ged_draw_shape_ref_source_snapshot(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   struct ged_draw_shape_source_snapshot *out)
{
    if (!out)
	return 0;
    memset(out, 0, sizeof(*out));

    struct ged_draw_obol_context_token *token =
	ged_draw_shape_ref_obol_token(gedp, ref);
    if (token && token->is_database_source && token->path &&
	    token->path[0]) {
	struct ged_draw_obol_database_source_runtime runtime;
	memset(&runtime, 0, sizeof(runtime));
	if (!ged_draw_obol_database_source_runtime_for_path(gedp,
		token->path, &runtime) || !runtime.valid || !runtime.dbip)
	    return 0;
	if (!token->fullpath_valid)
	    ged_draw_obol_context_token_fullpath_init(token);
	if (!token->fullpath_valid)
	    return 0;

	BN_TOL_INIT_SET_TOL(&token->obol_snapshot_tol);
	BG_TESS_TOL_INIT_SET_TOL(&token->obol_snapshot_ttol);
	if (runtime.tessellation_abs_tol >= 0.0)
	    token->obol_snapshot_ttol.abs = runtime.tessellation_abs_tol;
	if (runtime.tessellation_rel_tol >= 0.0)
	    token->obol_snapshot_ttol.rel = runtime.tessellation_rel_tol;
	if (runtime.tessellation_norm_tol >= 0.0)
	    token->obol_snapshot_ttol.norm = runtime.tessellation_norm_tol;

	out->dbip = runtime.dbip;
	out->fullpath = &token->fullpath;
	out->leaf_dp = DB_FULL_PATH_CUR_DIR(&token->fullpath);
	out->name = token->path;
	out->tol = &token->obol_snapshot_tol;
	out->ttol = &token->obol_snapshot_ttol;
	return 1;
    }

    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	shape_ref = ged_draw_scene_ref_from_rt_view_ref(
		ged_draw_registry_shape_ref_cache_rt_ref(gedp, ref));
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;

    return ged_draw_scene_ref_source_snapshot(shape_ref, out);
}


static int
ged_draw_rt_vlist_to_ged_line_command(int rt_cmd, int *ged_cmd,
				      int *duplicate)
{
    if (ged_cmd)
	*ged_cmd = -1;
    if (duplicate)
	*duplicate = 0;

    switch (rt_cmd) {
	case RT_VLIST_LINE_MOVE:
	case RT_VLIST_POLY_MOVE:
	case RT_VLIST_TRI_MOVE:
	    if (ged_cmd)
		*ged_cmd = GED_DRAW_VIEW_LINE_MOVE;
	    return 1;
	case RT_VLIST_LINE_DRAW:
	case RT_VLIST_POLY_DRAW:
	case RT_VLIST_POLY_END:
	case RT_VLIST_TRI_DRAW:
	case RT_VLIST_TRI_END:
	    if (ged_cmd)
		*ged_cmd = GED_DRAW_VIEW_LINE_DRAW;
	    return 1;
	case RT_VLIST_POINT_DRAW:
	    if (ged_cmd)
		*ged_cmd = GED_DRAW_VIEW_LINE_MOVE;
	    if (duplicate)
		*duplicate = 1;
	    return 1;
	default:
	    return 0;
    }
}


static int
ged_draw_rt_vlist_to_ged_line_arrays(struct bu_list *vhead,
				     point_t **points_out,
				     int **commands_out,
				     size_t *point_count_out)
{
    if (points_out)
	*points_out = NULL;
    if (commands_out)
	*commands_out = NULL;
    if (point_count_out)
	*point_count_out = 0;
    if (!vhead || !points_out || !commands_out || !point_count_out)
	return 0;

    size_t point_count = 0;
    rt_vlist *vp = NULL;
    BU_LIST_EACH(vhead, vp, rt_vlist) {
	for (size_t i = 0; i < vp->nused; i++) {
	    int duplicate = 0;
	    if (!ged_draw_rt_vlist_to_ged_line_command(vp->cmd[i], NULL,
		    &duplicate))
		continue;
	    if (point_count >= SIZE_MAX - (duplicate ? 2 : 1))
		return 0;
	    point_count += duplicate ? 2 : 1;
	}
    }

    point_t *points = NULL;
    int *commands = NULL;
    if (point_count) {
	points = (point_t *)bu_calloc(point_count, sizeof(point_t),
		"GED Obol RT VLIST line points");
	commands = (int *)bu_calloc(point_count, sizeof(int),
		"GED Obol RT VLIST line commands");
    }

    size_t idx = 0;
    BU_LIST_EACH(vhead, vp, rt_vlist) {
	for (size_t i = 0; i < vp->nused; i++) {
	    int command = -1;
	    int duplicate = 0;
	    if (!ged_draw_rt_vlist_to_ged_line_command(vp->cmd[i], &command,
		    &duplicate))
		continue;
	    VMOVE(points[idx], vp->pt[i]);
	    commands[idx++] = command;
	    if (duplicate) {
		VMOVE(points[idx], vp->pt[i]);
		commands[idx++] = GED_DRAW_VIEW_LINE_DRAW;
	    }
	}
    }

    *points_out = points;
    *commands_out = commands;
    *point_count_out = point_count;
    return 1;
}


struct ged_draw_obol_submodel_publish_ctx {
    struct ged *gedp;
    const char *source_path;
    const struct ged_draw_scene_display_summary *display_state;
    size_t child_count;
    int failed;
};


static union tree *
ged_draw_obol_submodel_wireframe_leaf(struct db_tree_state *tsp,
				      const struct db_full_path *pathp,
				      struct directory *dp,
				      void *client_data)
{
    struct ged_draw_obol_submodel_publish_ctx *ctx =
	(struct ged_draw_obol_submodel_publish_ctx *)client_data;
    char *path_name = NULL;
    const char *name = "submodel_leaf";
    struct bu_list vhead;
    point_t *points = NULL;
    int *commands = NULL;
    size_t point_count = 0;
    union tree *curtree = TREE_NULL;
    struct rt_db_internal local_ip;
    int have_local_ip = 0;

    if (!ctx || ctx->failed || !tsp || !dp)
	return TREE_NULL;
    RT_DB_INTERNAL_INIT(&local_ip);

    if (pathp && pathp->fp_len > 0) {
	path_name = db_path_to_string(pathp);
	if (path_name && path_name[0])
	    name = path_name;
    } else if (dp->d_namep) {
	name = dp->d_namep;
    }

    if (!tsp->ts_dbip ||
	    rt_db_get_internal(&local_ip, dp, tsp->ts_dbip, NULL) < 0) {
	ctx->failed = 1;
	goto cleanup;
    }
    have_local_ip = 1;

    BU_LIST_INIT(&vhead);
    if (!ged_draw_rt_internal_payload_valid(&local_ip)) {
	ctx->failed = 1;
	goto cleanup;
    }
    if (rt_obj_plot(&vhead, &local_ip, tsp->ts_ttol, tsp->ts_tol) < 0 ||
	    !ged_draw_rt_vlist_to_ged_line_arrays(&vhead, &points,
		&commands, &point_count)) {
	ctx->failed = 1;
	goto cleanup;
    }

    for (size_t i = 0; i < point_count; i++) {
	point_t transformed;
	MAT4X3PNT(transformed, tsp->ts_mat, points[i]);
	VMOVE(points[i], transformed);
    }

    if (!ged_draw_obol_database_source_publish_auxiliary_line_set_for_path(
	    ctx->gedp, ctx->source_path, name,
	    (const point_t *)points, commands, point_count,
	    ctx->display_state)) {
	ctx->failed = 1;
	goto cleanup;
    }

    ctx->child_count++;
    BU_GET(curtree, union tree);
    RT_TREE_INIT(curtree);
    curtree->tr_op = OP_NOP;

cleanup:
    if (points)
	bu_free(points, "GED Obol RT VLIST line points");
    if (commands)
	bu_free(commands, "GED Obol RT VLIST line commands");
    RT_FREE_VLIST(&rt_vlfree, &vhead);
    if (have_local_ip)
	rt_db_free_internal(&local_ip);
    if (path_name)
	bu_free(path_name, "GED Obol submodel leaf path string");
    return curtree;
}


static int
ged_draw_obol_database_source_publish_submodel_wireframe_for_path(
	struct ged *gedp,
	const char *source_path,
	struct rt_db_internal *ip,
	const struct bg_tess_tol *ttol,
	const struct bn_tol *tol)
{
    if (!gedp || !source_path || !source_path[0] || !ip ||
	    ip->idb_type != ID_SUBMODEL || !ip->idb_ptr)
	return 0;

    struct ged_draw_database_source_summary source_summary;
    memset(&source_summary, 0, sizeof(source_summary));
    if (!ged_draw_obol_database_source_summary_for_path(gedp, source_path,
	    &source_summary) || !source_summary.valid)
	return 0;

    struct rt_submodel_internal *sip =
	(struct rt_submodel_internal *)ip->idb_ptr;
    RT_SUBMODEL_CK_MAGIC(sip);

    struct db_i *dbip = DBI_NULL;
    int close_db = 0;
    if (bu_vls_strlen(&sip->file) != 0) {
	dbip = db_open(bu_vls_addr(&sip->file), DB_OPEN_READONLY);
	if (dbip == DBI_NULL)
	    return -1;
	close_db = 1;
	if (!db_is_directory_non_empty(dbip) && db_dirbuild(dbip) < 0) {
	    db_close(dbip);
	    return -1;
	}
    } else {
	RT_CK_DBI(sip->dbip);
	dbip = (struct db_i *)sip->dbip;
    }

    struct bn_tol local_tol;
    const struct bn_tol *use_tol = tol;
    if (!use_tol) {
	BN_TOL_INIT_SET_TOL(&local_tol);
	use_tol = &local_tol;
    }

    struct bg_tess_tol local_ttol;
    const struct bg_tess_tol *use_ttol = ttol;
    if (!use_ttol) {
	BG_TESS_TOL_INIT_SET_TOL(&local_ttol);
	use_ttol = &local_ttol;
    }

    (void)ged_draw_obol_database_source_clear_vlist_for_path(gedp,
	    source_path);
    (void)ged_draw_obol_database_source_clear_auxiliary_shapes_for_path(gedp,
	    source_path);

    struct ged_draw_scene_display_summary display_state;
    memset(&display_state, 0, sizeof(display_state));
    int have_display_state =
	ged_draw_obol_database_source_display_summary_for_path(gedp,
		source_path, &display_state) && display_state.valid;

    struct db_tree_state state;
    RT_DBTS_INIT(&state);
    state.ts_dbip = dbip;
    state.ts_ttol = use_ttol;
    state.ts_tol = use_tol;
    MAT_COPY(state.ts_mat, sip->root2leaf);

    struct ged_draw_obol_submodel_publish_ctx ctx;
    ctx.gedp = gedp;
    ctx.source_path = source_path;
    ctx.display_state = have_display_state ? &display_state : NULL;
    ctx.child_count = 0;
    ctx.failed = 0;

    const char *argv[2];
    argv[0] = bu_vls_addr(&sip->treetop);
    argv[1] = NULL;
    int ret = db_walk_tree_leaf_instances(dbip, 1, argv, 1, &state, 0, NULL,
	    ged_draw_obol_submodel_wireframe_leaf, (void *)&ctx);

    if (close_db)
	db_close(dbip);
    if (ret < 0 || ctx.failed)
	return -1;
    return ctx.child_count ? 1 : 0;
}


static int
ged_draw_rt_internal_payload_valid(const struct rt_db_internal *ip)
{
    if (!ip || !ip->idb_ptr || !ip->idb_meth)
	return 0;
    if (ip->idb_meth->magic != RT_FUNCTAB_MAGIC)
	return 0;
    if (!ip->idb_meth->ft_internal_magic)
	return 1;
    uint32_t actual = *((const uint32_t *)ip->idb_ptr);
    if (actual != ip->idb_meth->ft_internal_magic) {
	bu_log("ged draw invalid internal payload for %s: actual=0x%08x expected=0x%08x\n",
		ip->idb_meth->ft_name[0] ? ip->idb_meth->ft_name : "unknown",
		(unsigned int)actual,
		(unsigned int)ip->idb_meth->ft_internal_magic);
	return 0;
    }
    return 1;
}


static int
ged_draw_scene_ref_publish_obol_submodel_wireframe(
	bsg_scene_ref ref,
	struct rt_db_internal *ip,
	const struct bg_tess_tol *ttol,
	const struct bn_tol *tol)
{
    struct ged *gedp = _ged_draw_scene_ref_owner_gedp(ref);
    if (!gedp)
	return 0;

    (void)_ged_draw_scene_ref_obol_ensure_owner_source(ref);

    char *path = _ged_draw_scene_ref_owner_path_string(ref);
    if (path) {
	int result =
	    ged_draw_obol_database_source_publish_submodel_wireframe_for_path(
		    gedp, path, ip, ttol, tol);
	bu_free(path, "GED Obol submodel source path");
	if (result)
	    return result;
    }

    struct ged_draw_shape_source_snapshot source;
    memset(&source, 0, sizeof(source));
    if (!ged_draw_scene_ref_source_snapshot(ref, &source))
	return 0;

    return (source.name && source.name[0]) ?
	ged_draw_obol_database_source_publish_submodel_wireframe_for_path(
		gedp, source.name, ip, ttol, tol) : 0;
}


static int
ged_draw_scene_ref_publish_obol_primitive_wireframe(
	bsg_scene_ref ref,
	struct rt_db_internal *ip,
	const struct bg_tess_tol *ttol,
	const struct bn_tol *tol)
{
    struct ged *gedp = _ged_draw_scene_ref_owner_gedp(ref);
    if (!gedp)
	return 0;

    char *path = _ged_draw_scene_ref_owner_path_string(ref);
    if (path) {
	(void)ged_draw_scene_ref_obol_sync_source_placement(ref);
	int result =
	    ged_draw_obol_database_source_publish_primitive_wireframe_for_path(
		gedp, path, ip, ttol, tol);
	bu_free(path, "GED Obol primitive source path");
	if (result)
	    return result;
    }

    struct ged_draw_shape_source_snapshot source;
    memset(&source, 0, sizeof(source));
    if (!ged_draw_scene_ref_source_snapshot(ref, &source))
	return 0;

    (void)ged_draw_scene_ref_obol_sync_source_placement(ref);
    return (source.name && source.name[0]) ?
	ged_draw_obol_database_source_publish_primitive_wireframe_for_path(
	    gedp, source.name, ip, ttol, tol) : 0;
}


static int
ged_draw_shape_ref_redraw_wireframe_obol(struct ged *gedp,
					 ged_draw_shape_ref ref,
					 struct db_i *dbip,
					 const struct bn_tol *tol,
					 const struct bg_tess_tol *ttol,
					 void *UNUSED(view_ctx),
					 int skip_subtractions)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_shape_ref_obol_token(gedp, ref);
    if (!token || !token->is_database_source || !token->path ||
	    !token->path[0])
	return 0;

    struct ged_draw_obol_draw_state_summary draw_state;
    memset(&draw_state, 0, sizeof(draw_state));
    if (!ged_draw_obol_database_source_draw_state_for_path(gedp,
	    token->path, &draw_state) || !draw_state.valid ||
	    draw_state.draw_mode != GED_DRAW_MODE_WIRE)
	return 0;
    if (skip_subtractions && draw_state.line_style)
	return 0;
    if (!draw_state.draw_mat_valid)
	MAT_IDN(draw_state.draw_mat);

    struct ged_draw_obol_database_source_runtime runtime;
    memset(&runtime, 0, sizeof(runtime));
    if (!ged_draw_obol_database_source_runtime_for_path(gedp, token->path,
	    &runtime) || !runtime.valid)
	return -1;
    if (!runtime.dbip)
	runtime.dbip = dbip;
    if (!runtime.dbip)
	return -1;

    struct db_full_path local_fullpath;
    db_full_path_init(&local_fullpath);
    int free_local_fullpath = 0;
    const struct db_full_path *fullpath =
	token->fullpath_valid ? &token->fullpath : NULL;
    if (!fullpath && runtime.database_path && runtime.database_path[0]) {
	const char *path = ged_draw_dbpath_skip_lead_slash(
		runtime.database_path);
	if (path && path[0] &&
		db_string_to_path(&local_fullpath, runtime.dbip, path) == 0) {
	    fullpath = &local_fullpath;
	    free_local_fullpath = 1;
	}
    }
    if (!fullpath || fullpath->fp_len <= 0) {
	if (free_local_fullpath)
	    db_free_full_path(&local_fullpath);
	return -1;
    }

    const struct bn_tol *use_tol = tol;
    struct bn_tol local_tol;
    if (!use_tol) {
	BN_TOL_INIT_SET_TOL(&local_tol);
	use_tol = &local_tol;
    }

    const struct bg_tess_tol *use_ttol = ttol;
    struct bg_tess_tol local_ttol;
    if (!use_ttol) {
	BG_TESS_TOL_INIT_SET_TOL(&local_ttol);
	if (runtime.tessellation_abs_tol >= 0.0)
	    local_ttol.abs = runtime.tessellation_abs_tol;
	if (runtime.tessellation_rel_tol >= 0.0)
	    local_ttol.rel = runtime.tessellation_rel_tol;
	if (runtime.tessellation_norm_tol >= 0.0)
	    local_ttol.norm = runtime.tessellation_norm_tol;
	use_ttol = &local_ttol;
    }

    struct directory *dp = NULL;
    struct rt_db_internal local_ip;
    int have_local_ip = 0;
    int ret = -1;
    RT_DB_INTERNAL_INIT(&local_ip);
    if (!ged_draw_obol_database_source_set_placement_for_path(gedp,
	    token->path, 1, draw_state.draw_mat, 0, NULL, 0, 0.0))
	goto cleanup_fullpath;

    dp = DB_FULL_PATH_CUR_DIR(fullpath);
    if (!dp || !(dp->d_flags & RT_DIR_SOLID))
	goto cleanup_fullpath;

    if (rt_db_get_internal(&local_ip, dp, runtime.dbip, NULL) < 0)
	goto cleanup_fullpath;
    have_local_ip = 1;

    ret = ged_draw_obol_database_source_publish_primitive_wireframe_for_path(
	    gedp, token->path, &local_ip, use_ttol, use_tol) > 0 ? 0 : -1;

cleanup_fullpath:
    if (have_local_ip)
	rt_db_free_internal(&local_ip);
    if (free_local_fullpath)
	db_free_full_path(&local_fullpath);

    struct ged_draw_database_source_summary source_summary;
    memset(&source_summary, 0, sizeof(source_summary));
    if (ged_draw_obol_database_source_summary_for_path(gedp, token->path,
	    &source_summary) && source_summary.valid) {
	const int failed = ret < 0;
	(void)ged_draw_obol_database_source_set_realization_for_path(gedp,
		token->path, failed ? 0 : 1, failed ? 1 : 0,
		failed ? source_summary.realized_source_revision :
		    source_summary.source_revision,
		failed ? source_summary.realized_inputs_revision :
		    source_summary.inputs_revision,
		failed ? GED_DRAW_STALE_UPDATE_FAILED : GED_DRAW_STALE_NONE);
    }

    return ret;
}


static int
ged_draw_scene_ref_wireframe_shape(bsg_scene_ref shape_ref, void *view_ctx,
				   struct db_i *dbip,
				   const struct bn_tol *tol,
				   const struct bg_tess_tol *ttol)
{
    struct ged_draw_shape_record_summary shape_summary;
    if (!ged_draw_scene_ref_shape_record_summary(shape_ref, &shape_summary))
	return 0;
    const struct db_full_path *fp = shape_summary.fullpath;
    if (!fp || fp->fp_len <= 0)
	return 0;

    struct ged_draw_scene_draw_state_summary draw_state = {0};
    if (!ged_draw_scene_ref_draw_state_summary(shape_ref, &draw_state) ||
	    !draw_state.draw_mat_valid)
	return -1;

    struct rt_db_internal dbintern;
    RT_DB_INTERNAL_INIT(&dbintern);
    struct rt_db_internal *ip = &dbintern;
    int ret = -1;
    int get_ret = rt_db_get_internal(ip, DB_FULL_PATH_CUR_DIR(fp),
	    dbip, draw_state.draw_mat);
    if (get_ret < 0)
	return -1;

    ged_draw_view_lod_policy lod_policy;
    rt_view_context_lod_policy_get(&lod_policy, view_ctx);
    if (view_ctx && lod_policy.csg_enabled) {
	struct rt_view_info view_info;
	rt_view_context_info_get(&view_info, view_ctx);
	ret = ged_draw_scene_ref_publish_primitive_wireframe(shape_ref, ip,
		ttol, tol, view_ctx, &view_info, 1);
    }
    if (ret < 0)
	ret = ged_draw_scene_ref_publish_primitive_wireframe(shape_ref, ip,
		ttol, tol, view_ctx, NULL, 0);

    rt_db_free_internal(ip);

    if (ret < 0) {
	if (DB_FULL_PATH_CUR_DIR(fp)) {
	    bu_log("%s: plot failure\n",
		   DB_FULL_PATH_CUR_DIR(fp)->d_namep);
	}
	return -1;
    }

    return 0;
}


int
ged_draw_shape_ref_redraw_wireframe(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    struct db_i *dbip,
				    const struct bn_tol *tol,
				    const struct bg_tess_tol *ttol,
				    void *view_ctx,
				    int skip_subtractions)
{
    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return ged_draw_shape_ref_redraw_wireframe_obol(gedp, ref, dbip,
		tol, ttol, view_ctx, skip_subtractions);
    ged_draw_scene_ref_remember_owner_gedp(shape_ref, gedp);

    struct ged_draw_scene_draw_state_summary draw_state = {0};
    if (!ged_draw_scene_ref_draw_state_summary(shape_ref, &draw_state) ||
	    draw_state.draw_mode != GED_DRAW_MODE_WIRE)
	return 0;
    if (skip_subtractions && draw_state.line_style)
	return 0;

    struct ged_draw_shape_record_summary shape_summary;
    if (!ged_draw_scene_ref_shape_record_summary(shape_ref, &shape_summary))
	return 0;
    const struct db_full_path *fp = shape_summary.fullpath;
    if (!fp || fp->fp_len <= 0)
	return 0;

    struct directory *dp = DB_FULL_PATH_CUR_DIR(fp);
    if (dp && (dp->d_flags & RT_DIR_COMB))
	return 0;

    ged_draw_scene_ref_geometry_clear(shape_ref);
    int ret = ged_draw_scene_ref_wireframe_shape(shape_ref, view_ctx, dbip,
	    tol, ttol);
    ged_draw_scene_ref_mark_database_source_redraw_result(shape_ref,
	    ret < 0);
    return ret;
}


struct ged_draw_group_redraw_wireframe_ctx {
    struct ged *gedp;
    struct db_i *dbip;
    const struct bn_tol *tol;
    const struct bg_tess_tol *ttol;
    void *view_ctx;
    int skip_subtractions;
    int ret;
};


static int
ged_draw_group_ref_redraw_wireframe_shape_cb(bsg_scene_ref shape_ref,
					     void *userdata)
{
    struct ged_draw_group_redraw_wireframe_ctx *ctx =
	(struct ged_draw_group_redraw_wireframe_ctx *)userdata;
    if (!ctx || ged_draw_scene_ref_is_null(shape_ref))
	return 1;

    ged_draw_shape_ref ref =
	_ged_draw_shape_ref_from_scene_ref(ctx->gedp, shape_ref);
    ctx->ret += ged_draw_shape_ref_redraw_wireframe(ctx->gedp, ref,
	    ctx->dbip, ctx->tol, ctx->ttol, ctx->view_ctx,
	    ctx->skip_subtractions);

    return 1;
}


static int
ged_draw_group_ref_redraw_wireframe_obol_path_cb(struct ged *gedp,
						 const char *path,
						 void *userdata)
{
    struct ged_draw_group_redraw_wireframe_ctx *ctx =
	(struct ged_draw_group_redraw_wireframe_ctx *)userdata;
    if (!ctx || !path || !path[0])
	return 1;

    ctx->ret += ged_draw_obol_database_source_realize_for_path(gedp, path);
    return 1;
}


int
ged_draw_group_ref_redraw_wireframe(struct ged *gedp,
				    ged_draw_group_ref ref,
				    struct db_i *dbip,
				    const struct bn_tol *tol,
				    const struct bg_tess_tol *ttol,
				    void *view_ctx,
				    int skip_subtractions)
{
    struct ged_draw_group_redraw_wireframe_ctx ctx;
    ctx.gedp = gedp;
    ctx.dbip = dbip;
    ctx.tol = tol;
    ctx.ttol = ttol;
    ctx.view_ctx = view_ctx;
    ctx.skip_subtractions = skip_subtractions;
    ctx.ret = 0;

    struct ged_draw_group_record_summary group_summary;
    memset(&group_summary, 0, sizeof(group_summary));
    if (ged_draw_group_ref_record_summary(gedp, ref, &group_summary) &&
	    group_summary.path) {
	int obol_status = ged_draw_obol_group_database_source_paths_foreach(
		gedp, group_summary.path,
		ged_draw_group_ref_redraw_wireframe_obol_path_cb, &ctx);
	if (obol_status >= 0)
	    return ctx.ret;
    }

    bsg_scene_ref group_ref = _ged_draw_group_ref_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(group_ref))
	return -1;

    ged_draw_source_scene_ref_foreach_shape(group_ref,
	    ged_draw_group_ref_redraw_wireframe_shape_cb, &ctx);

    return ctx.ret;
}


static int
_ged_draw_shape_ref_release(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    int obol_removed = _ged_draw_shape_ref_try_obol_paths(gedp, ref,
	    shape_ref, _ged_draw_obol_remove_source_cb, NULL,
	    "GED Obol source release path");
    if (obol_removed) {
	ged_draw_obol_owner_structural_revision_bump(gedp);
	return 1;
    }

    /* Obol owns public shape release; retained source-owner release is fallback. */
    return ged_draw_scene_ref_release_source_owner(shape_ref);
}


static int
ged_draw_scene_rt_ref_set_visible(rt_view_scene_ref scene_ref, int visible)
{
    return ged_draw_scene_ref_set_visible(
	    ged_draw_scene_ref_from_rt_view_ref(scene_ref), visible);
}


static int
ged_draw_shape_ref_obol_update_display_lazy(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	bsg_scene_ref *shape_ref_out,
	int visible_valid,
	int visible,
	int highlighted_valid,
	int highlighted,
	int draw_mode_valid,
	int draw_mode,
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
	uint64_t material_revision)
{
    struct _ged_draw_obol_update_display_ctx ctx = {
	visible_valid,
	visible,
	highlighted_valid,
	highlighted,
	draw_mode_valid,
	draw_mode,
	line_style_valid,
	line_style,
	line_width_valid,
	line_width,
	transparency_valid,
	transparency,
	color_valid,
	color,
	material_color_valid,
	material_color,
	material_revision_valid,
	material_revision
    };
    return _ged_draw_shape_ref_try_obol_paths_lazy(gedp, ref, shape_ref_out,
	    _ged_draw_obol_update_display_cb, &ctx,
	    "GED Obol display mutation path");
}


int
ged_draw_shape_ref_set_visible(struct ged *gedp, ged_draw_shape_ref ref,
			       int visible)
{
    int obol_updated = ged_draw_shape_ref_obol_update_display_lazy(gedp,
	    ref, NULL, 1, visible, 0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0.0, 0, NULL, 0, NULL, 0, 0);
    /* Obol owns public shape visibility mutation; retained fallback is retired. */
    return obol_updated;
}


int
ged_draw_group_ref_set_visible(struct ged *gedp, ged_draw_group_ref ref,
			       int visible)
{
    bsg_scene_ref group_ref = _ged_draw_group_ref_scene_ref(gedp, ref);
    struct _ged_draw_obol_group_visible_ctx ctx = {1, visible};
    return _ged_draw_group_ref_try_obol_paths_ensure(gedp, ref,
	    group_ref, _ged_draw_obol_group_visible_cb, &ctx,
	    "GED Obol group visible path", NULL, GED_DRAW_MODE_WIRE, 0);
}


int
ged_draw_shape_ref_get_color(struct ged *gedp, ged_draw_shape_ref ref,
			     unsigned char rgb[3])
{
    if (!rgb)
	return 0;

    struct ged_draw_shape_material_summary material;
    if (!ged_draw_shape_ref_material_summary(gedp, ref, &material) ||
	    !material.valid)
	return 0;

    rgb[0] = material.material_color[0];
    rgb[1] = material.material_color[1];
    rgb[2] = material.material_color[2];
    return 1;
}


int
ged_draw_shape_ref_set_color(struct ged *gedp, ged_draw_shape_ref ref,
			     const unsigned char rgb[3])
{
    int obol_updated = ged_draw_shape_ref_obol_update_display_lazy(gedp,
	    ref, NULL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0.0, 1, rgb, 0, NULL, 0, 0);
    /* Obol owns public shape color mutation; retained fallback is retired. */
    return obol_updated;
}


int
ged_draw_shape_ref_set_highlighted(struct ged *gedp, ged_draw_shape_ref ref,
				   int highlighted)
{
    int obol_updated = ged_draw_shape_ref_obol_update_display_lazy(gedp,
	    ref, NULL, 0, 0, 1, highlighted, 0, 0, 0, 0, 0, 0,
	    0, 0.0, 0, NULL, 0, NULL, 0, 0);
    /* Obol owns public shape highlight mutation; retained fallback is retired. */
    return obol_updated;
}


int
ged_draw_shape_ref_set_transparency(struct ged *gedp, ged_draw_shape_ref ref,
				    fastf_t transparency)
{
    int obol_updated = ged_draw_shape_ref_obol_update_display_lazy(gedp,
	    ref, NULL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	    1, transparency, 0, NULL, 0, NULL, 0, 0);
    /* Obol owns public shape transparency mutation; retained fallback is retired. */
    return obol_updated;
}


int
ged_draw_shape_ref_mark_database_source_changed(struct ged *gedp,
						ged_draw_shape_ref ref,
						ged_draw_stale_reason reason)
{
    struct _ged_draw_obol_stale_ctx ctx = {reason};
    int obol_marked = _ged_draw_shape_ref_try_obol_paths_lazy(gedp, ref,
	    NULL, _ged_draw_obol_mark_stale_cb, &ctx,
	    "GED Obol stale source path");
    /* Obol owns public database-source stale mutation; retained fallback is retired. */
    return obol_marked;
}


int
ged_draw_shape_ref_set_material_color(struct ged *gedp,
				      ged_draw_shape_ref ref,
				      const unsigned char rgb[3])
{
    if (!rgb)
	return 0;
    int obol_updated = ged_draw_shape_ref_obol_update_display_lazy(gedp, ref,
	    NULL,
	    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0.0, 0, NULL, 1, rgb, 0, 0);
    /* Obol owns public shape material-color mutation; retained fallback is retired. */
    return obol_updated;
}


int
ged_draw_shape_ref_refresh_material_color(struct ged *gedp,
					  ged_draw_shape_ref ref,
					  struct db_i *dbip,
					  uint64_t mater_rev)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_shape_ref_obol_token(gedp, ref);
    if (token && token->is_database_source && token->path &&
	    token->path[0])
	return ged_draw_obol_database_source_refresh_material_color_for_path(
		gedp, token->path, dbip, mater_rev);

    return 0;
}


int
ged_draw_shape_ref_set_evaluated_region(struct ged *gedp,
					ged_draw_shape_ref ref,
					int evaluated_region)
{
    struct _ged_draw_obol_set_evaluated_region_ctx ctx = {
	evaluated_region
    };
    int obol_updated = _ged_draw_shape_ref_try_obol_paths_lazy(gedp, ref,
	    NULL, _ged_draw_obol_set_evaluated_region_cb, &ctx,
	    "GED Obol evaluated-region mutation path");
    /* Obol owns public evaluated-region mutation; retained fallback is retired. */
    return obol_updated;
}


static int
ged_draw_shape_ref_lod_ensure_obol(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   void *first_view_ctx)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_shape_ref_obol_token(gedp, ref);
    if (!token || !token->is_database_source || !token->path ||
	    !token->path[0] || !first_view_ctx)
	return 0;

    ged_draw_view_lod_policy policy;
    rt_view_context_lod_policy_get(&policy, first_view_ctx);

    struct ged_draw_obol_draw_state_summary draw_state;
    memset(&draw_state, 0, sizeof(draw_state));
    if (!ged_draw_obol_database_source_draw_state_for_path(gedp,
	    token->path, &draw_state) || !draw_state.valid)
	return 1;

    struct ged_draw_obol_realization_policy_summary existing_policy;
    memset(&existing_policy, 0, sizeof(existing_policy));
    int candidate =
	ged_draw_obol_database_source_realization_policy_for_path(gedp,
	    token->path, &existing_policy) &&
	existing_policy.valid &&
	(existing_policy.role_flags &
	 (GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_CSG |
	  GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_MESH));

    struct db_full_path local_fullpath;
    db_full_path_init(&local_fullpath);
    int free_local_fullpath = 0;
    const struct db_full_path *fullpath =
	token->fullpath_valid ? &token->fullpath : NULL;
    struct ged_draw_obol_database_source_runtime runtime;
    memset(&runtime, 0, sizeof(runtime));
    (void)ged_draw_obol_database_source_runtime_for_path(gedp, token->path,
	    &runtime);
    if (!fullpath && runtime.valid && runtime.dbip && runtime.database_path &&
	    runtime.database_path[0]) {
	const char *path = ged_draw_dbpath_skip_lead_slash(
		runtime.database_path);
	if (path && path[0] &&
		db_string_to_path(&local_fullpath, runtime.dbip, path) == 0) {
	    fullpath = &local_fullpath;
	    free_local_fullpath = 1;
	}
    }

    struct directory *dp = (fullpath && fullpath->fp_len > 0) ?
	DB_FULL_PATH_CUR_DIR(fullpath) : NULL;
    const int mode = draw_state.draw_mode;
    int roles = GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_NONE;
    if (policy.csg_enabled && mode == GED_DRAW_MODE_WIRE)
	roles |= GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_CSG;
    if (dp && policy.mesh_enabled &&
	    ((dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT &&
	      (mode == GED_DRAW_MODE_WIRE ||
	       mode == GED_DRAW_MODE_SHADED_BOTS)) ||
	     (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP &&
	      (mode == GED_DRAW_MODE_SHADED_BOTS ||
	       mode == GED_DRAW_MODE_SHADED))))
	roles |= GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_MESH;

    if (!candidate && roles == GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_NONE) {
	if (free_local_fullpath)
	    db_free_full_path(&local_fullpath);
	return 1;
    }

    if (!ged_draw_obol_database_source_set_realization_roles_for_path(gedp,
	    token->path, roles)) {
	if (free_local_fullpath)
	    db_free_full_path(&local_fullpath);
	return 0;
    }

    int policy_updated =
	ged_draw_obol_database_source_set_realization_view_policy_for_path(
	    gedp,
	    token->path,
	    (policy.csg_enabled || policy.mesh_enabled) ? 1 : 0,
	    rt_view_context_scale_get(first_view_ctx),
	    (uint64_t)policy.bot_threshold,
	    policy.curve_scale,
	    policy.point_scale);
    if (!policy_updated) {
	if (free_local_fullpath)
	    db_free_full_path(&local_fullpath);
	return 0;
    }

    if ((roles & GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_CSG) &&
	    !(roles & GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_MESH) &&
	    dp && !(dp->d_flags & RT_DIR_COMB)) {
	struct db_i *realize_dbip = runtime.valid && runtime.dbip ?
	    runtime.dbip : gedp->dbip;
	struct bn_tol local_tol = BN_TOL_INIT_TOL;
	mat_t draw_mat;
	if (draw_state.draw_mat_valid)
	    MAT_COPY(draw_mat, draw_state.draw_mat);
	else
	    MAT_IDN(draw_mat);
	(void)ged_draw_obol_database_source_adaptive_wireframe_realize(gedp,
		token->path, realize_dbip, fullpath, draw_mat, &local_tol,
		first_view_ctx, runtime.draw_size_valid ? runtime.draw_size :
		    0.0);
    }

    if ((roles & GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_MESH) && dp &&
	    dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT) {
	struct db_i *realize_dbip = runtime.valid && runtime.dbip ?
	    runtime.dbip : gedp->dbip;
	(void)ged_draw_obol_database_source_bot_mesh_lod_realize(gedp,
		token->path, realize_dbip, dp, first_view_ctx);
    }
    if ((roles & GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_MESH) && dp &&
	    dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP) {
	struct db_i *realize_dbip = runtime.valid && runtime.dbip ?
	    runtime.dbip : gedp->dbip;
	struct bg_tess_tol local_ttol;
	BG_TESS_TOL_INIT_SET_TOL(&local_ttol);
	if (runtime.tessellation_abs_tol >= 0.0)
	    local_ttol.abs = runtime.tessellation_abs_tol;
	if (runtime.tessellation_rel_tol >= 0.0)
	    local_ttol.rel = runtime.tessellation_rel_tol;
	if (runtime.tessellation_norm_tol >= 0.0)
	    local_ttol.norm = runtime.tessellation_norm_tol;
	struct bn_tol local_tol = BN_TOL_INIT_TOL;
	(void)ged_draw_obol_database_source_brep_mesh_lod_realize(gedp,
		token->path, realize_dbip, dp, &local_ttol, &local_tol,
		first_view_ctx);
    }

    if (free_local_fullpath)
	db_free_full_path(&local_fullpath);
    return 1;
}


int
ged_draw_shape_ref_lod_ensure(struct ged *gedp, ged_draw_shape_ref ref,
			      void *first_view_ctx,
			      void **view_ctxs,
			      size_t view_ctx_count)
{
    if (!gedp || ged_draw_shape_ref_is_null(ref))
	return 0;

    struct ged_draw_obol_context_token *token =
	ged_draw_shape_ref_obol_token(gedp, ref);
    if (!token || !token->is_database_source || !token->path ||
	    !token->path[0])
	return 0;

    int ensured = 0;
    for (size_t i = 0; i < view_ctx_count; i++) {
	if (view_ctxs && view_ctxs[i] &&
		ged_draw_shape_ref_lod_ensure_obol(gedp, ref, view_ctxs[i]))
	    ensured = 1;
    }
    if (!ensured && first_view_ctx)
	ensured = ged_draw_shape_ref_lod_ensure_obol(gedp, ref,
		first_view_ctx);

    return ensured;
}


void *
ged_draw_shape_ref_view_context(struct ged *gedp, ged_draw_shape_ref ref)
{
    if (ged_draw_shape_ref_obol_token(gedp, ref))
	return ged_draw_active_view_ctx(gedp);

    bsg_scene_ref shape_ref = _ged_draw_shape_ref_runtime_scene_ref(gedp,
	    ref);
    struct ged_draw_scene_draw_state_summary draw_state = {0};
    if (!ged_draw_scene_ref_draw_state_summary(shape_ref, &draw_state))
	return NULL;
    return draw_state.view_ctx;
}


void
ged_draw_overlay_erase_name_context(struct ged *gedp, void *view_ctx,
				    const char *name)
{
    if (!view_ctx || !name)
	return;

    if (ged_draw_obol_view_context_feature_store_active(view_ctx))
	(void)ged_draw_obol_overlay_erase_name_context(gedp, view_ctx, name);
}


int
ged_draw_overlay_geometry_insert_context(
	struct ged *gedp,
	void *view_ctx,
	const char *name,
	const struct ged_draw_overlay_geometry *geometry,
	const unsigned char basecolor[3],
	fastf_t transparency,
	int draw_mode,
	int csoltab,
	ged_draw_shape_ref *out)
{
    if (out)
	*out = GED_DRAW_SHAPE_REF_NULL;
    if (!gedp || !view_ctx || !name || !geometry || !basecolor)
	return -1;

    if (!ged_draw_obol_view_context_feature_store_active(view_ctx))
	return 0;

    char *obol_shape_path = NULL;
    int ret = ged_draw_obol_overlay_geometry_insert_context(gedp,
	    view_ctx, name, geometry, basecolor, transparency, draw_mode,
	    &obol_shape_path);
    if (!ret)
	return -1;
    if (out && obol_shape_path)
	*out = ged_draw_obol_shape_ref_for_path(gedp, obol_shape_path, 0);
    if (obol_shape_path)
	bu_free(obol_shape_path, "GED Obol overlay shape path");
    ged_draw_obol_owner_structural_revision_bump(gedp);
    (void)csoltab;
    return (!out || !ged_draw_shape_ref_is_null(*out)) ? 0 : -1;
}


typedef int (*ged_draw_mesh_lod_publish_cb)(
	void *data,
	const point_t *points,
	size_t point_count,
	const vect_t *normals,
	size_t normal_count,
	const int *indices,
	size_t index_count);


static int
ged_draw_mesh_lod_publish_current(
	struct BRLObolMeshLod *lod,
	ged_draw_mesh_lod_publish_cb cb,
	void *cb_data)
{
    struct BRLObolMeshLodData lod_data;
    if (!lod || !cb || !brlobol_mesh_lod_data_get(lod, &lod_data))
	return 0;

    const int *faces = lod_data.faces;
    const point_t *points = lod_data.points;
    size_t face_count = lod_data.face_count;
    size_t point_count = lod_data.point_count;
    size_t valid_faces = 0;
    for (size_t i = 0; i < face_count; i++) {
	int valid = 1;
	for (int j = 0; j < 3; j++) {
	    int idx = faces[3*i + j];
	    if (idx < 0 || (size_t)idx >= point_count) {
		valid = 0;
		break;
	    }
	}
	if (valid)
	    valid_faces++;
    }
    if (!valid_faces)
	return 0;

    size_t index_count = valid_faces * 4;
    size_t normal_count = valid_faces * 3;
    int *indices = (int *)bu_calloc(index_count, sizeof(int), "mesh lod indexed-face indices");
    vect_t *normals = (vect_t *)bu_calloc(normal_count, sizeof(vect_t), "mesh lod indexed-face normals");
    const point_t *normal_points =
	(lod_data.points_orig && lod_data.point_orig_count >= point_count) ?
	lod_data.points_orig : points;
    size_t out = 0;
    size_t normal_out = 0;

    for (size_t i = 0; i < face_count; i++) {
	int bad_face = 0;
	for (int j = 0; j < 3; j++) {
	    int idx = faces[3*i + j];
	    if (idx < 0 || (size_t)idx >= point_count) {
		bad_face = 1;
		break;
	    }
	}
	if (bad_face)
	    continue;

	vect_t face_normal, ab, ac;
	VSUB2(ab, normal_points[faces[3*i + 0]], normal_points[faces[3*i + 1]]);
	VSUB2(ac, normal_points[faces[3*i + 0]], normal_points[faces[3*i + 2]]);
	VCROSS(face_normal, ab, ac);
	VUNITIZE(face_normal);

	for (int j = 0; j < 3; j++) {
	    indices[out] = faces[3*i + j];
	    if (lod_data.normals && lod_data.normal_count >= (3*i + j + 1)) {
		VMOVE(normals[normal_out], lod_data.normals[3*i + j]);
		if (ZERO(MAGNITUDE(normals[normal_out])))
		    VMOVE(normals[normal_out], face_normal);
	    } else {
		VMOVE(normals[normal_out], face_normal);
	    }
	    normal_out++;
	    out++;
	}
	indices[out++] = -1;
    }

    int ret = (*cb)(cb_data, points, point_count, normals, normal_count,
	    indices, index_count);
    bu_free(normals, "mesh lod indexed-face normals");
    bu_free(indices, "mesh lod indexed-face indices");
    return ret;
}


struct ged_draw_mesh_lod_obol_publish {
    struct ged *gedp;
    const char *path;
};


static int
ged_draw_mesh_lod_obol_publish_cb(
	void *data,
	const point_t *points,
	size_t point_count,
	const vect_t *normals,
	size_t normal_count,
	const int *indices,
	size_t index_count)
{
    struct ged_draw_mesh_lod_obol_publish *ctx =
	(struct ged_draw_mesh_lod_obol_publish *)data;
    return ctx && ctx->gedp && ctx->path ?
	ged_draw_obol_database_source_publish_lod_indexed_face_set_for_path(
		ctx->gedp, ctx->path, points, point_count, normals,
		normal_count, indices, index_count) : 0;
}


static int
ged_draw_obol_database_source_publish_current_mesh_lod(
	struct ged *gedp,
	const char *path,
	struct BRLObolMeshLod *lod)
{
    struct ged_draw_mesh_lod_obol_publish ctx;
    ctx.gedp = gedp;
    ctx.path = path;
    return ged_draw_mesh_lod_publish_current(lod,
	    ged_draw_mesh_lod_obol_publish_cb, &ctx);
}


static int
ged_draw_obol_database_source_adaptive_wireframe_realize(
	struct ged *gedp,
	const char *path,
	struct db_i *dbip,
	const struct db_full_path *fullpath,
	const mat_t draw_mat,
	const struct bn_tol *tol,
	void *view_ctx,
	fastf_t draw_size)
{
    if (!gedp || !path || !path[0] || !dbip || !fullpath ||
	    fullpath->fp_len <= 0 || !draw_mat || !view_ctx)
	return 0;

    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, DB_FULL_PATH_CUR_DIR(fullpath), dbip,
	    draw_mat) < 0)
	return 0;

    int published = 0;
    if (intern.idb_meth && intern.idb_meth->ft_lod_realize) {
	struct rt_primitive_lod_realization realization;
	struct rt_view_info view_info = RT_VIEW_INFO_INIT;
	memset(&realization, 0, sizeof(realization));
	ged_draw_view_context_info_get(&view_info, view_ctx);
	if (intern.idb_meth->ft_lod_realize(&realization, &intern, tol,
		&view_info, draw_size) >= 0 && realization.has_line_set) {
	    published =
		ged_draw_obol_database_source_publish_realization_line_set(
			gedp, path, &realization);
	} else {
	    ged_draw_primitive_realization_line_set_free(&realization);
	}
    }

    rt_db_free_internal(&intern);

    if (published) {
	struct ged_draw_obol_database_source_record record;
	memset(&record, 0, sizeof(record));
	if (ged_draw_obol_database_source_record_for_path(gedp, path,
		&record)) {
	    (void)ged_draw_obol_database_source_set_realization_for_path(gedp,
		    path, 1, 0, record.source_revision,
		    record.inputs_revision, GED_DRAW_STALE_NONE);
	}
    }

    return published;
}


static int
ged_draw_obol_database_source_bot_mesh_lod_realize(
	struct ged *gedp,
	const char *path,
	struct db_i *dbip,
	struct directory *dp,
	void *view_ctx)
{
    if (!gedp || !path || !path[0] || !dbip || !dp || !view_ctx)
	return 0;

    struct ged_draw_obol_database_source_runtime runtime;
    memset(&runtime, 0, sizeof(runtime));
    if (!ged_draw_obol_database_source_runtime_for_path(gedp, path,
	    &runtime) || !runtime.valid)
	return 0;

    struct BRLObolMeshLod *mesh_lod = runtime.mesh_lod;
    if (!mesh_lod) {
	struct BRLObolMeshLodCacheStatus status = BRLOBOL_MESH_LOD_CACHE_STATUS_INIT;
	if (brlobol_mesh_lod_cache_status(dbip, dp->d_namep, &status) != BRLCAD_OK)
	    return 0;
	if (!status.has_cache_key || !status.has_cached_payload ||
		status.stale_cache_entry) {
	    if (brlobol_mesh_lod_cache_refresh(dbip, dp->d_namep, &status) != BRLCAD_OK)
		return 0;
	}
	if (!status.has_cache_key || !status.has_cached_payload)
	    return 0;

	mesh_lod = brlobol_mesh_lod_get(dbip, dp->d_namep);
	if (!mesh_lod) {
	    if (brlobol_mesh_lod_cache_refresh(dbip, dp->d_namep, &status) != BRLCAD_OK ||
		    !status.has_cache_key || !status.has_cached_payload)
		return 0;
	    mesh_lod = brlobol_mesh_lod_get(dbip, dp->d_namep);
	}
	if (!mesh_lod)
	    return 0;

	if (!ged_draw_obol_database_source_set_mesh_lod_for_path(gedp, path,
		mesh_lod)) {
	    brlobol_mesh_lod_destroy(mesh_lod);
	    return 0;
	}
    }

    struct rt_view_info view_info = RT_VIEW_INFO_INIT;
    ged_draw_view_context_info_get(&view_info, view_ctx);
    int level = brlobol_mesh_lod_load_view(mesh_lod, &view_info, 0);
    if (level < 0)
	bu_log("Error loading info for initial Obol LoD view\n");

    struct BRLObolMeshLodInfo info = BRLOBOL_MESH_LOD_INFO_INIT;
    if (brlobol_mesh_lod_info_get(mesh_lod, &info))
	(void)ged_draw_obol_database_source_set_mesh_lod_bounds_for_path(gedp,
		path, info.bmin, info.bmax);

    int published = ged_draw_obol_database_source_publish_current_mesh_lod(
	    gedp, path, mesh_lod);
    if (published) {
	struct ged_draw_obol_database_source_record record;
	memset(&record, 0, sizeof(record));
	if (ged_draw_obol_database_source_record_for_path(gedp, path,
		&record)) {
	    (void)ged_draw_obol_database_source_set_realization_for_path(gedp,
		    path, 1, 0, record.source_revision,
		    record.inputs_revision, GED_DRAW_STALE_NONE);
	}
    }
    return published;
}


static int
ged_draw_obol_database_source_brep_mesh_lod_realize(
	struct ged *gedp,
	const char *path,
	struct db_i *dbip,
	struct directory *dp,
	const struct bg_tess_tol *ttol,
	const struct bn_tol *tol,
	void *view_ctx)
{
    if (!gedp || !path || !path[0] || !dbip || !dp || !view_ctx)
	return 0;

    struct ged_draw_obol_database_source_runtime runtime;
    memset(&runtime, 0, sizeof(runtime));
    if (!ged_draw_obol_database_source_runtime_for_path(gedp, path,
	    &runtime) || !runtime.valid)
	return 0;

    struct BRLObolMeshLod *mesh_lod = runtime.mesh_lod;
    if (!mesh_lod) {
	point_t bmin, bmax;
	int bounds_valid = 0;
	if (!ged_draw_brep_mesh_lod_cache_prepare(&mesh_lod, bmin, bmax,
		&bounds_valid, dbip, dp, ttol, tol))
	    return 0;
	if (!mesh_lod)
	    return 0;

	if (!ged_draw_brep_mesh_lod_detail_setup(mesh_lod, dbip, dp, ttol,
		tol)) {
	    brlobol_mesh_lod_destroy(mesh_lod);
	    return 0;
	}

	if (!ged_draw_obol_database_source_set_mesh_lod_for_path(gedp, path,
		mesh_lod)) {
	    brlobol_mesh_lod_destroy(mesh_lod);
	    return 0;
	}

	if (bounds_valid)
	    (void)ged_draw_obol_database_source_set_mesh_lod_bounds_for_path(
		    gedp, path, bmin, bmax);
    }

    struct rt_view_info view_info = RT_VIEW_INFO_INIT;
    ged_draw_view_context_info_get(&view_info, view_ctx);
    int level = brlobol_mesh_lod_load_view(mesh_lod, &view_info, 0);
    if (level < 0)
	bu_log("Error loading info for initial Obol BREP LoD view\n");

    struct BRLObolMeshLodInfo info = BRLOBOL_MESH_LOD_INFO_INIT;
    if (brlobol_mesh_lod_info_get(mesh_lod, &info))
	(void)ged_draw_obol_database_source_set_mesh_lod_bounds_for_path(gedp,
		path, info.bmin, info.bmax);

    int published = ged_draw_obol_database_source_publish_current_mesh_lod(
	    gedp, path, mesh_lod);
    if (published) {
	struct ged_draw_obol_database_source_record record;
	memset(&record, 0, sizeof(record));
	if (ged_draw_obol_database_source_record_for_path(gedp, path,
		&record)) {
	    (void)ged_draw_obol_database_source_set_realization_for_path(gedp,
		    path, 1, 0, record.source_revision,
		    record.inputs_revision, GED_DRAW_STALE_NONE);
	}
    }
    return published;
}


static struct ged *
_ged_draw_scene_ref_owner_gedp(bsg_scene_ref ref)
{
    ged_draw_shape_state *shape_data = _ged_draw_shape_state_get_scene_ref(ref);
    return shape_data ? shape_data->gedp : NULL;
}


static void
ged_draw_scene_ref_remember_owner_gedp(bsg_scene_ref ref, struct ged *gedp)
{
    ged_draw_shape_state *shape_data = _ged_draw_shape_state_get_scene_ref(ref);
    if (shape_data && !shape_data->gedp)
	shape_data->gedp = gedp;
}


static int
ged_draw_scene_ref_obol_source_path_apply(
	bsg_scene_ref ref,
	ged_draw_scene_ref_obol_source_path_cb cb,
	void *data)
{
    struct ged *gedp = _ged_draw_scene_ref_owner_gedp(ref);
    if (!gedp || !cb)
	return 0;

    bsg_scene_ref source_ref = ged_draw_scene_ref_source_owner(ref);
    const char *source_path = ged_draw_scene_rt_ref_semantic_path(
	    ged_draw_scene_ref_to_rt_view_ref(source_ref));
    if (source_path && source_path[0] && (*cb)(gedp, source_path, data))
	return 1;

    source_path = ged_draw_scene_ref_draw_intent_path(source_ref);
    if (source_path && source_path[0] && (*cb)(gedp, source_path, data))
	return 1;

    source_path = ged_draw_scene_ref_draw_intent_path(ref);
    if (source_path && source_path[0] && (*cb)(gedp, source_path, data))
	return 1;

    const struct db_full_path *fp = ged_draw_scene_ref_fullpath(ref);
    if (fp) {
	char *path = db_path_to_string(fp);
	if (path) {
	    int matched = (*cb)(gedp, path, data);
	    bu_free(path, "GED Obol source fullpath");
	    if (matched)
		return 1;
	}
    }

    ged_draw_shape_ref shape_ref = _ged_draw_shape_ref_from_scene_ref(gedp,
	    ref);
    return _ged_draw_shape_ref_try_obol_paths(gedp, shape_ref, ref,
	    (ged_draw_shape_ref_obol_path_cb)cb, data,
	    "GED Obol source path");
}


static int
ged_draw_scene_ref_obol_source_path_apply_ensure(
	bsg_scene_ref ref,
	ged_draw_scene_ref_obol_source_path_cb cb,
	void *data)
{
    if (ged_draw_scene_ref_obol_source_path_apply(ref, cb, data))
	return 1;

    if (!_ged_draw_scene_ref_obol_ensure_source(ref))
	return 0;

    return ged_draw_scene_ref_obol_source_path_apply(ref, cb, data);
}


struct ged_draw_scene_ref_obol_realization_update {
    int current;
    int failed;
    uint64_t realized_source_revision;
    uint64_t realized_inputs_revision;
    ged_draw_stale_reason reason;
};


static int
ged_draw_scene_ref_obol_set_realization_cb(
	struct ged *gedp,
	const char *path,
	void *data)
{
    struct ged_draw_scene_ref_obol_realization_update *update =
	(struct ged_draw_scene_ref_obol_realization_update *)data;
    return update ?
	ged_draw_obol_database_source_set_realization_for_path(gedp, path,
		update->current, update->failed,
		update->realized_source_revision,
		update->realized_inputs_revision,
		update->reason) : 0;
}


static int
ged_draw_scene_ref_obol_set_realization(
	bsg_scene_ref ref,
	int current,
	int failed,
	uint64_t realized_source_revision,
	uint64_t realized_inputs_revision,
	ged_draw_stale_reason reason)
{
    struct ged_draw_scene_ref_obol_realization_update update;
    update.current = current;
    update.failed = failed;
    update.realized_source_revision = realized_source_revision;
    update.realized_inputs_revision = realized_inputs_revision;
    update.reason = reason;

    return ged_draw_scene_ref_obol_source_path_apply(ref,
	    ged_draw_scene_ref_obol_set_realization_cb, &update);
}


static int
ged_draw_scene_ref_obol_realization_policy_summary_cb(
	struct ged *gedp,
	const char *path,
	void *data)
{
    return ged_draw_obol_database_source_realization_policy_for_path(gedp,
	    path, (struct ged_draw_obol_realization_policy_summary *)data);
}


static int
ged_draw_scene_ref_obol_realization_policy_summary(
	bsg_scene_ref ref,
	struct ged_draw_obol_realization_policy_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    return ged_draw_scene_ref_obol_source_path_apply(ref,
	    ged_draw_scene_ref_obol_realization_policy_summary_cb, out);
}


static int
ged_draw_scene_ref_geometry_clear(bsg_scene_ref ref)
{
    ged_draw_shape_state *shape_data = _ged_draw_shape_state_get_scene_ref(ref);
    int obol_cleared = 0;
    if (!shape_data)
	return 0;
    if (shape_data->gedp) {
	ged_draw_shape_ref shape_ref =
	    _ged_draw_shape_ref_from_scene_ref(shape_data->gedp, ref);
	obol_cleared = ged_draw_shape_ref_obol_clear_vlist(shape_data->gedp,
		shape_ref, ref);
	obol_cleared =
	    ged_draw_shape_ref_obol_clear_mesh(shape_data->gedp, shape_ref,
		    ref) ||
	    obol_cleared;
	obol_cleared =
	    ged_draw_shape_ref_obol_clear_auxiliary_shapes(shape_data->gedp,
		    shape_ref, ref) ||
	    obol_cleared;
    }

    int bsg_cleared = 0;
    if (!obol_cleared) {
	bsg_scene_ref geometry_ref =
	    ged_draw_scene_ref_geometry_publication_target(ref);
	if (!bsg_scene_ref_is_null(geometry_ref)) {
	    (void)bsg_geometry_ref_clear(
		    bsg_scene_ref_as_geometry(geometry_ref));
	    bsg_cleared = 1;
	}
    }
    if (!bsg_cleared && !obol_cleared)
	return 0;

    /* Obol owns source-adapter geometry clear; BSG is fallback. */
    shape_data->geometry_command_count = 0;
    shape_data->geometry_revision++;
    bsg_scene_invalidate(ref);
    return 1;
}


static void
_ged_draw_bounds_include_point(point_t bmin, point_t bmax, const point_t pt)
{
    VMINMAX(bmin, bmax, pt);
}


static void
_ged_draw_bounds_include_point_marker(point_t bmin, point_t bmax, const point_t pt)
{
    point_t marker_min;
    point_t marker_max;

    VSET(marker_min, pt[X] - 1.0, pt[Y] - 1.0, pt[Z] - 1.0);
    VSET(marker_max, pt[X] + 1.0, pt[Y] + 1.0, pt[Z] + 1.0);
    VMIN(bmin, marker_min);
    VMAX(bmax, marker_max);
}


static int
_ged_draw_line_geometry_bounds(bsg_geometry_ref geometry,
			       point_t bmin,
			       point_t bmax,
			       size_t *length,
			       int *bad_cmd)
{
    bsg_field_ref coords = bsg_geometry_ref_coordinates_field(geometry);
    bsg_field_ref commands = bsg_geometry_ref_primitive_sets_field(geometry);
    size_t count = bsg_field_multi_count(coords);

    for (size_t i = 0; i < count; i++) {
	point_t pt;
	int cmd = (i % 2) ? BSG_GEOMETRY_LINE_DRAW : BSG_GEOMETRY_LINE_MOVE;

	if (!bsg_field_multi_point_at(coords, i, pt))
	    return 0;
	(void)bsg_field_multi_int_at(commands, i, &cmd);

	switch (cmd) {
	    case BSG_GEOMETRY_LINE_MOVE:
	    case BSG_GEOMETRY_LINE_DRAW:
		_ged_draw_bounds_include_point(bmin, bmax, pt);
		break;
	    case BSG_GEOMETRY_POINT_DRAW:
		_ged_draw_bounds_include_point_marker(bmin, bmax, pt);
		break;
	    default:
		if (bad_cmd)
		    *bad_cmd = cmd;
		return 0;
	}
    }

    if (length)
	*length = count;
    return count ? 1 : 0;
}


static int
_ged_draw_point_geometry_bounds(bsg_geometry_ref geometry,
				point_t bmin,
				point_t bmax,
				size_t *length)
{
    bsg_field_ref coords = bsg_geometry_ref_coordinates_field(geometry);
    size_t count = bsg_field_multi_count(coords);

    for (size_t i = 0; i < count; i++) {
	point_t pt;
	if (!bsg_field_multi_point_at(coords, i, pt))
	    return 0;
	_ged_draw_bounds_include_point_marker(bmin, bmax, pt);
    }

    if (length)
	*length = count;
    return count ? 1 : 0;
}


static int
_ged_draw_indexed_geometry_bounds(bsg_geometry_ref geometry,
				  point_t bmin,
				  point_t bmax,
				  size_t *length)
{
    bsg_field_ref coords = bsg_geometry_ref_coordinates_field(geometry);
    bsg_field_ref indices = bsg_geometry_ref_indices_field(geometry);
    size_t coord_count = bsg_field_multi_count(coords);
    size_t index_count = bsg_field_multi_count(indices);
    size_t used = 0;

    for (size_t i = 0; i < index_count; i++) {
	int idx = -1;
	point_t pt;

	if (!bsg_field_multi_int_at(indices, i, &idx))
	    return 0;
	if (idx < 0)
	    continue;
	if ((size_t)idx >= coord_count)
	    return 0;
	if (!bsg_field_multi_point_at(coords, (size_t)idx, pt))
	    return 0;
	_ged_draw_bounds_include_point(bmin, bmax, pt);
	used++;
    }

    if (length)
	*length = index_count;
    return used ? 1 : 0;
}


static int
_ged_draw_geometry_line_command_from_bsg(int command)
{
    switch (command) {
	case BSG_GEOMETRY_LINE_MOVE:
	    return GED_DRAW_VIEW_LINE_MOVE;
	case BSG_GEOMETRY_LINE_DRAW:
	    return GED_DRAW_VIEW_LINE_DRAW;
	default:
	    return command;
    }
}


static const char *
_ged_draw_geometry_kind_name(bsg_geometry_node_kind kind)
{
    switch (kind) {
	case BSG_GEOMETRY_NODE_LINE_SET:
	    return "line-set";
	case BSG_GEOMETRY_NODE_INDEXED_LINE_SET:
	    return "indexed-line-set";
	case BSG_GEOMETRY_NODE_POINT_SET:
	    return "point-set";
	case BSG_GEOMETRY_NODE_INDEXED_FACE_SET:
	    return "indexed-face-set";
	case BSG_GEOMETRY_NODE_MESH:
	    return "mesh";
	case BSG_GEOMETRY_NODE_TEXT:
	    return "text";
	case BSG_GEOMETRY_NODE_IMAGE:
	    return "image";
	case BSG_GEOMETRY_NODE_FRAMEBUFFER:
	    return "framebuffer";
	case BSG_GEOMETRY_NODE_OVERLAY:
	    return "overlay";
	case BSG_GEOMETRY_NODE_HUD:
	    return "hud";
	case BSG_GEOMETRY_NODE_ANNOTATION:
	    return "annotation";
	case BSG_GEOMETRY_NODE_CSG_PROXY:
	    return "csg-proxy";
	case BSG_GEOMETRY_NODE_BREP_PROXY:
	    return "brep-proxy";
	case BSG_GEOMETRY_NODE_EDIT_PREVIEW:
	    return "edit-preview";
	case BSG_GEOMETRY_NODE_NONE:
	default:
	    return "none";
    }
}


static int
_ged_draw_geometry_last_face_point(bsg_geometry_ref geometry, point_t out)
{
    bsg_field_ref coords = bsg_geometry_ref_coordinates_field(geometry);
    bsg_field_ref indices = bsg_geometry_ref_indices_field(geometry);
    size_t index_count = bsg_field_multi_count(indices);
    int first = -1;
    int vert = 0;
    int have = 0;

    for (size_t i = 0; i < index_count; i++) {
	int idx = -1;
	if (!bsg_field_multi_int_at(indices, i, &idx))
	    continue;
	if (idx < 0) {
	    if (first >= 0 && vert >= 3 &&
		    bsg_field_multi_point_at(coords, (size_t)first, out))
		have = 1;
	    first = -1;
	    vert = 0;
	    continue;
	}
	if (first < 0)
	    first = idx;
	vert++;
    }

    if (first >= 0 && vert >= 3 &&
	    bsg_field_multi_point_at(coords, (size_t)first, out))
	have = 1;

    return have;
}


static int
_ged_draw_scene_ref_line_geometry(bsg_scene_ref ref, bsg_geometry_ref *out)
{
    if (out)
	memset(out, 0, sizeof(*out));
    if (bsg_scene_ref_is_null(ref) || !out)
	return 0;

    bsg_scene_ref geometry_ref = ged_draw_scene_ref_geometry_query_target(ref);
    if (bsg_scene_ref_is_null(geometry_ref))
	return 0;

    bsg_geometry_ref geometry = bsg_scene_ref_as_geometry(geometry_ref);
    if (bsg_geometry_ref_kind(geometry) != BSG_GEOMETRY_NODE_LINE_SET)
	return 0;

    *out = geometry;
    return 1;
}


static int
_ged_draw_scene_ref_geometry_bounds(bsg_scene_ref ref,
				    point_t bmin,
				    point_t bmax,
				    size_t *length,
				    int *bad_cmd)
{
    if (length)
	*length = 0;
    if (bad_cmd)
	*bad_cmd = 0;

    bsg_scene_ref geometry_ref = ged_draw_scene_ref_geometry_query_target(ref);
    if (bsg_scene_ref_is_null(geometry_ref))
	return 0;

    bsg_geometry_ref geometry = bsg_scene_ref_as_geometry(geometry_ref);
    switch (bsg_geometry_ref_kind(geometry)) {
	case BSG_GEOMETRY_NODE_LINE_SET:
	    return _ged_draw_line_geometry_bounds(geometry, bmin, bmax,
		    length, bad_cmd);
	case BSG_GEOMETRY_NODE_POINT_SET:
	    return _ged_draw_point_geometry_bounds(geometry, bmin, bmax, length);
	case BSG_GEOMETRY_NODE_INDEXED_LINE_SET:
	case BSG_GEOMETRY_NODE_INDEXED_FACE_SET:
	    return _ged_draw_indexed_geometry_bounds(geometry, bmin, bmax, length);
	default:
	    return 0;
    }
}


static int
ged_draw_scene_ref_geometry_summary(bsg_scene_ref ref,
				    struct ged_draw_shape_geometry_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (bsg_scene_ref_is_null(ref))
	return 0;

    bsg_scene_ref geometry_ref = ged_draw_scene_ref_geometry_query_target(ref);
    if (bsg_scene_ref_is_null(geometry_ref))
	return 0;

    bsg_geometry_ref geometry = bsg_scene_ref_as_geometry(geometry_ref);
    bsg_geometry_node_kind kind = bsg_geometry_ref_kind(geometry);
    if (kind == BSG_GEOMETRY_NODE_NONE)
	return 0;

    out->valid = 1;
    out->geometry_name = _ged_draw_geometry_kind_name(kind);
    out->point_count =
	bsg_field_multi_count(bsg_geometry_ref_coordinates_field(geometry));
    out->index_count =
	bsg_field_multi_count(bsg_geometry_ref_indices_field(geometry));
    return 1;
}


static int
ged_draw_scene_ref_last_point(bsg_scene_ref ref, point_t out)
{
    if (bsg_scene_ref_is_null(ref) || !out)
	return 0;

    bsg_scene_ref geometry_ref = ged_draw_scene_ref_geometry_query_target(ref);
    if (bsg_scene_ref_is_null(geometry_ref))
	return 0;

    bsg_geometry_ref geometry = bsg_scene_ref_as_geometry(geometry_ref);
    bsg_geometry_node_kind kind = bsg_geometry_ref_kind(geometry);
    if (kind == BSG_GEOMETRY_NODE_NONE)
	return 0;
    if (kind == BSG_GEOMETRY_NODE_INDEXED_FACE_SET)
	return _ged_draw_geometry_last_face_point(geometry, out);

    bsg_field_ref coords = bsg_geometry_ref_coordinates_field(geometry);
    size_t count = bsg_field_multi_count(coords);
    if (!count)
	return 0;
    return bsg_field_multi_point_at(coords, count - 1, out);
}


static int
ged_draw_scene_ref_line_summary(bsg_scene_ref ref,
				struct ged_draw_view_line_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    bsg_geometry_ref geometry;
    if (!_ged_draw_scene_ref_line_geometry(ref, &geometry))
	return 0;

    out->valid = 1;
    out->point_count =
	bsg_field_multi_count(bsg_geometry_ref_coordinates_field(geometry));
    return 1;
}


static int
ged_draw_scene_ref_line_point_at(bsg_scene_ref ref,
				 size_t index,
				 point_t out)
{
    if (!out)
	return 0;

    bsg_geometry_ref geometry;
    if (!_ged_draw_scene_ref_line_geometry(ref, &geometry))
	return 0;

    return bsg_field_multi_point_at(
	    bsg_geometry_ref_coordinates_field(geometry), index, out);
}


static int
ged_draw_scene_ref_line_command_at(bsg_scene_ref ref,
				   size_t index,
				   int *out)
{
    if (!out)
	return 0;

    bsg_geometry_ref geometry;
    if (!_ged_draw_scene_ref_line_geometry(ref, &geometry))
	return 0;

    int command = -1;
    if (!bsg_field_multi_int_at(
	    bsg_geometry_ref_primitive_sets_field(geometry), index, &command))
	return 0;
    *out = _ged_draw_geometry_line_command_from_bsg(command);
    return 1;
}


static int
_ged_draw_translate_line_set(bsg_geometry_ref geometry, const vect_t xlate)
{
    bsg_field_ref coords = bsg_geometry_ref_coordinates_field(geometry);
    bsg_field_ref primitives = bsg_geometry_ref_primitive_sets_field(geometry);
    size_t count = bsg_field_multi_count(coords);
    point_t *points = NULL;
    int *commands = NULL;
    int ok = 0;

    if (!count)
	return 0;

    points = (point_t *)bu_calloc(count, sizeof(point_t),
	    "translated line-set points");
    commands = (int *)bu_calloc(count, sizeof(int),
	    "translated line-set commands");
    for (size_t i = 0; i < count; i++) {
	point_t pt;
	int cmd = BSG_GEOMETRY_LINE_MOVE;
	if (!bsg_field_multi_point_at(coords, i, pt))
	    goto cleanup;
	(void)bsg_field_multi_int_at(primitives, i, &cmd);
	VADD2(points[i], pt, xlate);
	commands[i] = cmd;
    }

    ok = bsg_geometry_ref_set_line_set(geometry, (const point_t *)points,
	    commands, count);

cleanup:
    if (points)
	bu_free(points, "translated line-set points");
    if (commands)
	bu_free(commands, "translated line-set commands");
    return ok;
}


static int
_ged_draw_translate_point_set(bsg_geometry_ref geometry, const vect_t xlate)
{
    bsg_field_ref coords = bsg_geometry_ref_coordinates_field(geometry);
    size_t count = bsg_field_multi_count(coords);
    point_t *points = NULL;
    int ok = 0;

    if (!count)
	return 0;

    points = (point_t *)bu_calloc(count, sizeof(point_t),
	    "translated point-set points");
    for (size_t i = 0; i < count; i++) {
	point_t pt;
	if (!bsg_field_multi_point_at(coords, i, pt))
	    goto cleanup;
	VADD2(points[i], pt, xlate);
    }

    ok = bsg_geometry_ref_set_point_set(geometry, (const point_t *)points, count);

cleanup:
    if (points)
	bu_free(points, "translated point-set points");
    return ok;
}


static int
_ged_draw_translate_face_set(bsg_geometry_ref geometry, const vect_t xlate)
{
    bsg_field_ref coords = bsg_geometry_ref_coordinates_field(geometry);
    bsg_field_ref normals_field = bsg_geometry_ref_normals_field(geometry);
    bsg_field_ref indices_field = bsg_geometry_ref_indices_field(geometry);
    size_t point_count = bsg_field_multi_count(coords);
    size_t normal_count = bsg_field_multi_count(normals_field);
    size_t index_count = bsg_field_multi_count(indices_field);
    point_t *points = NULL;
    vect_t *normals = NULL;
    int *indices = NULL;
    int ok = 0;

    if (!point_count || !index_count)
	return 0;

    points = (point_t *)bu_calloc(point_count, sizeof(point_t),
	    "translated face-set points");
    normals = normal_count ? (vect_t *)bu_calloc(normal_count, sizeof(vect_t),
	    "translated face-set normals") : NULL;
    indices = (int *)bu_calloc(index_count, sizeof(int),
	    "translated face-set indices");

    for (size_t i = 0; i < point_count; i++) {
	point_t pt;
	if (!bsg_field_multi_point_at(coords, i, pt))
	    goto cleanup;
	VADD2(points[i], pt, xlate);
    }
    for (size_t i = 0; i < normal_count; i++) {
	vect_t normal;
	if (!bsg_field_multi_point_at(normals_field, i, normal))
	    goto cleanup;
	VMOVE(normals[i], normal);
    }
    for (size_t i = 0; i < index_count; i++) {
	if (!bsg_field_multi_int_at(indices_field, i, &indices[i]))
	    goto cleanup;
    }

    ok = bsg_geometry_ref_set_indexed_face_set(geometry,
	    (const point_t *)points, point_count,
	    (const vect_t *)normals, normal_count,
	    indices, index_count);

cleanup:
    if (points)
	bu_free(points, "translated face-set points");
    if (normals)
	bu_free(normals, "translated face-set normals");
    if (indices)
	bu_free(indices, "translated face-set indices");
    return ok;
}


static int
ged_draw_scene_ref_translate_geometry(bsg_scene_ref ref, const vect_t xlate)
{
    if (bsg_scene_ref_is_null(ref) || !xlate)
	return 0;

    bsg_geometry_ref geometry = bsg_scene_ref_as_geometry(ref);
    int ok = 0;
    switch (bsg_geometry_ref_kind(geometry)) {
	case BSG_GEOMETRY_NODE_LINE_SET:
	    ok = _ged_draw_translate_line_set(geometry, xlate);
	    break;
	case BSG_GEOMETRY_NODE_POINT_SET:
	    ok = _ged_draw_translate_point_set(geometry, xlate);
	    break;
	case BSG_GEOMETRY_NODE_INDEXED_FACE_SET:
	    ok = _ged_draw_translate_face_set(geometry, xlate);
	    break;
	default:
	    return 0;
    }

    if (ok) {
	ged_draw_shape_state *shape_data = ged_draw_scene_ref_shape_state(ref);
	if (shape_data)
	    shape_data->geometry_revision++;
	ged_draw_scene_ref_invalidate(ref);
    }
    return ok;
}


int
ged_draw_shape_ref_last_point(struct ged *gedp, ged_draw_shape_ref ref, point_t out)
{
    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_last_point_cb, out,
	    "GED Obol line last-point path"))
	return 1;

    return ged_draw_scene_ref_last_point(shape_ref, out);
}


int
ged_draw_shape_ref_line_summary(struct ged *gedp,
				ged_draw_shape_ref ref,
				struct ged_draw_view_line_summary *out)
{
    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    struct _ged_draw_obol_line_summary_ctx ctx = {out};
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_line_summary_cb, &ctx,
	    "GED Obol line summary path"))
	return 1;

    return ged_draw_scene_ref_line_summary(shape_ref, out);
}


int
ged_draw_shape_ref_line_point_at(struct ged *gedp,
				 ged_draw_shape_ref ref,
				 size_t index,
				 point_t out)
{
    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    struct _ged_draw_obol_line_point_ctx ctx;
    ctx.index = index;
    ctx.out = out;
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_line_point_cb, &ctx,
	    "GED Obol line point path"))
	return 1;

    return ged_draw_scene_ref_line_point_at(shape_ref, index, out);
}


int
ged_draw_shape_ref_line_command_at(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   size_t index,
				   int *out)
{
    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    struct _ged_draw_obol_line_command_ctx ctx = {index, out};
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_line_command_cb, &ctx,
	    "GED Obol line command path"))
	return 1;

    return ged_draw_scene_ref_line_command_at(shape_ref, index, out);
}


int
ged_draw_shape_ref_geometry_summary(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    struct ged_draw_shape_geometry_summary *out)
{
    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    struct _ged_draw_obol_geometry_summary_ctx ctx = {out};
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_geometry_summary_cb, &ctx,
	    "GED Obol geometry summary path"))
	return 1;

    return ged_draw_scene_ref_geometry_summary(shape_ref, out);
}


int
ged_draw_shape_context_line_summary(void *shape_ctx,
				    struct ged_draw_view_line_summary *out)
{
    struct _ged_draw_obol_line_summary_ctx ctx = {out};
    if (ged_draw_shape_context_obol_path_apply(shape_ctx,
	    _ged_draw_obol_line_summary_cb, &ctx,
	    "GED Obol context line summary path"))
	return 1;

    return ged_draw_scene_ref_line_summary(
	    ged_draw_scene_ref_from_context(shape_ctx), out);
}


int
ged_draw_shape_context_line_point_at(void *shape_ctx,
				     size_t index,
				     point_t out)
{
    struct _ged_draw_obol_line_point_ctx ctx;
    ctx.index = index;
    ctx.out = out;
    if (ged_draw_shape_context_obol_path_apply(shape_ctx,
	    _ged_draw_obol_line_point_cb, &ctx,
	    "GED Obol context line point path"))
	return 1;

    return ged_draw_scene_ref_line_point_at(
	    ged_draw_scene_ref_from_context(shape_ctx), index, out);
}


int
ged_draw_shape_context_line_command_at(void *shape_ctx,
				       size_t index,
				       int *out)
{
    struct _ged_draw_obol_line_command_ctx ctx = {index, out};
    if (ged_draw_shape_context_obol_path_apply(shape_ctx,
	    _ged_draw_obol_line_command_cb, &ctx,
	    "GED Obol context line command path"))
	return 1;

    return ged_draw_scene_ref_line_command_at(
	    ged_draw_scene_ref_from_context(shape_ctx), index, out);
}


int
ged_draw_shape_context_geometry_summary(void *shape_ctx,
					struct ged_draw_shape_geometry_summary *out)
{
    struct _ged_draw_obol_geometry_summary_ctx ctx = {out};
    if (ged_draw_shape_context_obol_path_apply(shape_ctx,
	    _ged_draw_obol_geometry_summary_cb, &ctx,
	    "GED Obol context geometry summary path"))
	return 1;

    return ged_draw_scene_ref_geometry_summary(
	    ged_draw_scene_ref_from_context(shape_ctx), out);
}


int
ged_draw_shape_ref_translate_geometry(struct ged *gedp, ged_draw_shape_ref ref,
				      const vect_t xlate)
{
    bsg_scene_ref shape_ref = _ged_draw_shape_ref_scene_ref(gedp, ref);
    struct _ged_draw_obol_translate_ctx ctx;
    VMOVE(ctx.xlate, xlate);
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref, shape_ref,
	    _ged_draw_obol_translate_cb, &ctx,
	    "GED Obol translate VLIST path"))
	return 1;

    return ged_draw_scene_ref_translate_geometry(shape_ref, xlate);
}


static int
ged_draw_scene_ref_update_bounds_from_geometry(bsg_scene_ref ref, int *bad_cmd)
{
    if (bad_cmd)
	*bad_cmd = 0;
    if (bsg_scene_ref_is_null(ref))
	return 0;

    point_t bmin, bmax;
    size_t length = 0;
    VSET(bmin, INFINITY, INFINITY, INFINITY);
    VSET(bmax, -INFINITY, -INFINITY, -INFINITY);

    if (!_ged_draw_scene_ref_geometry_bounds(ref, bmin, bmax, &length, bad_cmd))
	return 0;

    _ged_draw_scene_ref_geometry_set_command_count(ref, length);

    vect_t center;
    center[X] = (bmin[X] + bmax[X]) * 0.5;
    center[Y] = (bmin[Y] + bmax[Y]) * 0.5;
    center[Z] = (bmin[Z] + bmax[Z]) * 0.5;
    ged_draw_scene_ref_set_draw_center(ref, center);

    fastf_t size = bmax[X] - bmin[X];
    V_MAX(size, bmax[Y] - bmin[Y]);
    V_MAX(size, bmax[Z] - bmin[Z]);
    (void)ged_draw_scene_ref_set_draw_size(ref, size);
    ged_draw_scene_ref_set_non_database_source(ref, 0);
    return 1;
}


static void
ged_draw_scene_ref_set_draw_center(bsg_scene_ref ref, const point_t center)
{
    if (ged_draw_scene_ref_obol_set_placement(ref,
	    0, NULL,
	    1, center,
	    0, 0.0))
	return;

    bsg_scene_set_draw_center(ref, center);
}


static int
ged_draw_scene_ref_bounds(bsg_scene_ref ref, point_t min, point_t max)
{
    if (bsg_scene_ref_is_null(ref) || !min || !max)
	return 0;

    bsg_scene_bounds(ref, min, max);
    return 1;
}


static int
ged_draw_scene_ref_set_bounds(bsg_scene_ref ref,
			      const point_t min,
			      const point_t max,
			      int valid)
{
    if (bsg_scene_ref_is_null(ref) || !min || !max)
	return 0;

    if (valid && ged_draw_scene_ref_obol_set_bounds(ref, min, max))
	return 1;

    bsg_scene_set_bounds(ref, min, max, valid);
    return 1;
}


static int
ged_draw_scene_ref_subtree_bounds(bsg_scene_ref ref,
				  vect_t *min,
				  vect_t *max,
				  int include_overlays)
{
    if (bsg_scene_ref_is_null(ref) || !min || !max)
	return 1;

    return bsg_scene_subtree_bbox(ref, min, max, include_overlays);
}


static int
ged_draw_scene_ref_draw_tree_depth(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? 0 : bsg_scene_draw_tree_depth(ref);
}


static int
ged_draw_scene_ref_update_bounds_context(bsg_scene_ref ref, void *view_ctx)
{
    return bsg_scene_update_bounds(ref, (struct bsg_view *)view_ctx);
}


static void *
ged_draw_scene_ref_view_context(bsg_scene_ref ref)
{
    return (void *)bsg_scene_view(ref);
}


static void *
ged_draw_scene_ref_publication_view_context(bsg_scene_ref ref)
{
    void *view_ctx = ged_draw_scene_ref_view_context(ref);
    if (view_ctx)
	return view_ctx;

    ged_draw_shape_state *shape_data = _ged_draw_shape_state_get_scene_ref(ref);
    if (!shape_data || !shape_data->gedp)
	return NULL;

    return ged_draw_active_view_ctx(shape_data->gedp);
}


static bsg_scene_ref
ged_draw_view_context_geometry_create(void *view_ctx, const char *name)
{
    if (!view_ctx)
	return bsg_scene_ref_null();

    bsg_geometry_ref geometry_ref =
	bsg_geometry_ref_create((struct bsg_view *)view_ctx,
		name ? name : "geometry");
    return bsg_geometry_ref_as_scene(geometry_ref);
}


struct ged_draw_geometry_child_find_ctx {
    bsg_scene_ref geometry_ref;
    int require_geometry;
};


static int
_ged_draw_scene_ref_find_geometry_child_cb(bsg_scene_ref child_ref,
					   void *userdata)
{
    struct ged_draw_geometry_child_find_ctx *ctx =
	(struct ged_draw_geometry_child_find_ctx *)userdata;

    if (!ctx)
	return 0;
    if (ged_draw_scene_ref_is_shape(child_ref) &&
	    (!ctx->require_geometry ||
	     bsg_geometry_ref_kind(bsg_scene_ref_as_geometry(child_ref)) !=
	     BSG_GEOMETRY_NODE_NONE)) {
	ctx->geometry_ref = child_ref;
	return 0;
    }
    return 1;
}


static bsg_scene_ref
ged_draw_scene_ref_geometry_query_target(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return ged_draw_scene_ref_null();

    if (bsg_geometry_ref_kind(bsg_scene_ref_as_geometry(ref)) !=
	    BSG_GEOMETRY_NODE_NONE)
	return ref;

    struct ged_draw_geometry_child_find_ctx ctx;
    ctx.geometry_ref = ged_draw_scene_ref_null();
    ctx.require_geometry = 1;
    (void)ged_draw_scene_ref_foreach_child(ref,
	    _ged_draw_scene_ref_find_geometry_child_cb, &ctx);
    return ctx.geometry_ref;
}


static bsg_scene_ref
ged_draw_scene_ref_geometry_publication_target(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return ged_draw_scene_ref_null();
    if (ged_draw_scene_ref_is_shape(ref))
	return ref;

    struct ged_draw_geometry_child_find_ctx ctx;
    ctx.geometry_ref = ged_draw_scene_ref_null();
    ctx.require_geometry = 0;
    (void)ged_draw_scene_ref_foreach_child(ref,
	    _ged_draw_scene_ref_find_geometry_child_cb, &ctx);
    if (!ged_draw_scene_ref_is_null(ctx.geometry_ref))
	return ctx.geometry_ref;

    ged_draw_shape_state *shape_data = _ged_draw_shape_state_get_scene_ref(ref);
    if (!shape_data || !shape_data->gedp)
	return ged_draw_scene_ref_null();

    void *view_ctx = ged_draw_scene_ref_publication_view_context(ref);
    if (!view_ctx)
	return ged_draw_scene_ref_null();

    bsg_scene_ref geometry_ref =
	ged_draw_view_context_geometry_create(view_ctx,
		ged_draw_scene_ref_name(ref) ? ged_draw_scene_ref_name(ref) :
		"geometry");
    if (ged_draw_scene_ref_is_null(geometry_ref))
	return ged_draw_scene_ref_null();

    if (!ged_draw_scene_ref_append_child(ref, geometry_ref)) {
	ged_draw_scene_ref_release(geometry_ref);
	return ged_draw_scene_ref_null();
    }
    _ged_draw_scene_ref_copy_display_state(geometry_ref, ref);
    (void)_ged_draw_scene_ref_apply_db_object_marker(geometry_ref);
    (void)ged_draw_scene_ref_set_non_database_source(geometry_ref, 0);
    return geometry_ref;
}


static int
ged_draw_scene_ref_set_visible(bsg_scene_ref ref, int visible)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    if (ged_draw_scene_ref_obol_update_display(ref,
	    1, visible,
	    0, 0,
	    0, 0,
	    0, 0,
	    0, 0,
	    0, 0.0,
	    0, NULL,
	    0, NULL,
	    0, 0))
	return 1;

    bsg_scene_set_visible(ref, visible ? 1 : 0);
    return 1;
}


static int
ged_draw_scene_ref_set_transparency(bsg_scene_ref ref, fastf_t transparency)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    if (ged_draw_scene_ref_obol_update_display(ref,
	    0, 0,
	    0, 0,
	    0, 0,
	    0, 0,
	    0, 0,
	    1, transparency,
	    0, NULL,
	    0, NULL,
	    0, 0))
	return 1;

    bsg_scene_set_transparency(ref, transparency);
    return 1;
}


static int
ged_draw_scene_ref_set_color(bsg_scene_ref ref, const unsigned char rgb[3])
{
    if (bsg_scene_ref_is_null(ref) || !rgb)
	return 0;

    if (ged_draw_scene_ref_obol_update_display(ref,
	    0, 0,
	    0, 0,
	    0, 0,
	    0, 0,
	    0, 0,
	    0, 0.0,
	    1, rgb,
	    0, NULL,
	    0, 0))
	return 1;

    bsg_scene_set_color(ref, rgb[0], rgb[1], rgb[2]);
    return 1;
}


static int
ged_draw_scene_ref_set_highlighted(bsg_scene_ref ref, int highlighted)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    if (ged_draw_scene_ref_obol_update_display(ref,
	    0, 0,
	    1, highlighted,
	    0, 0,
	    0, 0,
	    0, 0,
	    0, 0.0,
	    0, NULL,
	    0, NULL,
	    0, 0))
	return 1;

    bsg_scene_set_highlighted(ref, highlighted ? 1 : 0);
    return 1;
}


static int
ged_draw_scene_ref_set_transform(bsg_scene_ref ref, const mat_t mat)
{
    if (bsg_scene_ref_is_null(ref) || !mat)
	return 0;

    bsg_scene_set_transform(ref, mat);
    return 1;
}


static int
ged_draw_scene_ref_set_draw_mat(bsg_scene_ref ref, const mat_t mat)
{
    if (bsg_scene_ref_is_null(ref) || !mat)
	return 0;

    if (ged_draw_scene_ref_obol_set_placement(ref,
	    1, mat,
	    0, NULL,
	    0, 0.0))
	return 1;

    bsg_scene_set_draw_mat(ref, mat);
    return 1;
}


static int
ged_draw_scene_ref_set_draw_size(bsg_scene_ref ref, fastf_t size)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    if (ged_draw_scene_ref_obol_set_placement(ref,
	    0, NULL,
	    0, NULL,
	    1, size))
	return 1;

    bsg_scene_set_draw_size(ref, size);
    return 1;
}


static int
ged_draw_scene_ref_set_line_style(bsg_scene_ref ref, int dashed)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    if (ged_draw_scene_ref_obol_update_display(ref,
	    0, 0,
	    0, 0,
	    0, 0,
	    1, dashed,
	    0, 0,
	    0, 0.0,
	    0, NULL,
	    0, NULL,
	    0, 0))
	return 1;

    bsg_scene_set_line_style(ref, dashed);
    return 1;
}


static int
ged_draw_scene_ref_set_line_width(bsg_scene_ref ref, int line_width)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    if (ged_draw_scene_ref_obol_update_display(ref,
	    0, 0,
	    0, 0,
	    0, 0,
	    0, 0,
	    1, line_width,
	    0, 0.0,
	    0, NULL,
	    0, NULL,
	    0, 0))
	return 1;

    bsg_scene_set_line_width(ref, line_width);
    return 1;
}


static int
ged_draw_scene_ref_draw_mat(bsg_scene_ref ref, mat_t mat)
{
    if (bsg_scene_ref_is_null(ref) || !mat)
	return 0;
    bsg_scene_draw_mat(ref, mat);
    return 1;
}


static int
ged_draw_scene_ref_draw_mode(bsg_scene_ref ref)
{
    return bsg_scene_dmode(ref);
}


static void
ged_draw_scene_ref_set_draw_mode(bsg_scene_ref ref, int draw_mode)
{
    if (ged_draw_scene_ref_obol_update_display(ref,
	    0, 0,
	    0, 0,
	    1, draw_mode,
	    0, 0,
	    0, 0,
	    0, 0.0,
	    0, NULL,
	    0, NULL,
	    0, 0))
	return;

    bsg_scene_set_dmode(ref, draw_mode);
    struct bsg_draw_intent *di = bsg_scene_draw_intent(ref);
    if (di)
	bsg_draw_intent_set_mode(di, (bsg_draw_mode)draw_mode);
}


static int
ged_draw_scene_ref_set_changed(bsg_scene_ref ref, int changed)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    bsg_scene_set_changed(ref, changed);
    return 1;
}


static int
ged_draw_scene_ref_ensure_draw_intent(bsg_scene_ref ref,
				      const char *path,
				      int draw_mode)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    if (ged_draw_scene_ref_is_group(ref) &&
	    !ged_draw_scene_ref_is_database_source(ref)) {
	struct _ged_draw_obol_group_draw_intent_ctx ctx = {
	    path,
	    1,
	    draw_mode,
	    0,
	    0
	};
	return ged_draw_scene_ref_obol_group_path_apply_ensure(ref,
		_ged_draw_obol_group_draw_intent_cb, &ctx,
		"GED Obol group ensure draw-intent path",
		path, draw_mode, 0);
    }

    if (bsg_scene_draw_intent(ref))
	return 1;

    struct bsg_draw_intent *di =
	bsg_draw_intent_create(path, (bsg_draw_mode)draw_mode);
    if (!di)
	return 0;

    bsg_scene_set_draw_intent(ref, di);
    return 1;
}


static int
ged_draw_scene_ref_set_draw_intent_path(bsg_scene_ref ref,
					const char *path,
					int fallback_draw_mode)
{
    if (bsg_scene_ref_is_null(ref) || !path)
	return 0;

    if (ged_draw_scene_ref_is_group(ref) &&
	    !ged_draw_scene_ref_is_database_source(ref)) {
	struct _ged_draw_obol_group_draw_intent_ctx ctx = {
	    path,
	    0,
	    0,
	    0,
	    0
	};
	return ged_draw_scene_ref_obol_group_path_apply_ensure(ref,
		_ged_draw_obol_group_draw_intent_cb, &ctx,
		"GED Obol group draw-intent path",
		path, fallback_draw_mode,
		ged_draw_scene_ref_draw_intent_is_overlay(ref));
    }

    struct bsg_draw_intent *di = bsg_scene_draw_intent(ref);
    if (di) {
	bsg_draw_intent_set_path(di, path);
	return 1;
    }

    return ged_draw_scene_ref_ensure_draw_intent(ref, path, fallback_draw_mode);
}


static int
ged_draw_scene_ref_apply_path_state(struct ged *gedp,
				    bsg_scene_ref ref,
				    const struct db_full_path *path)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 0;

    int indexed = 0;
    if (ged_draw_scene_ref_is_group(ref)) {
	ged_draw_group_ref group_ref =
	    _ged_draw_group_ref_from_scene_ref(gedp, ref);
	indexed = ged_draw_registry_group_ref_set_indexed_fullpath(gedp,
		group_ref, path);
    } else {
	ged_draw_shape_ref shape_ref =
	    _ged_draw_shape_ref_from_scene_ref(gedp, ref);
	indexed = ged_draw_registry_shape_ref_set_indexed_fullpath(gedp,
		shape_ref, path);
    }
    if (!indexed)
	return 0;

    if (path && path->fp_len > 0) {
	char *path_name = db_path_to_string(path);
	if (path_name) {
	    const char *semantic_path = path_name;
	    while (*semantic_path == '/')
		semantic_path++;
	    (void)ged_draw_scene_ref_set_draw_intent_path(ref,
		    semantic_path, GED_DRAW_MODE_WIRE);
	    bu_free(path_name, "ged draw shape intent path string");
	}
    }

    return 1;
}


static int
ged_draw_scene_ref_set_draw_intent_mode(bsg_scene_ref ref,
					int draw_mode,
					const char *fallback_path)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    if (ged_draw_scene_ref_is_group(ref) &&
	    !ged_draw_scene_ref_is_database_source(ref)) {
	struct _ged_draw_obol_group_draw_intent_ctx ctx = {
	    fallback_path,
	    1,
	    draw_mode,
	    0,
	    0
	};
	return ged_draw_scene_ref_obol_group_path_apply_ensure(ref,
		_ged_draw_obol_group_draw_intent_cb, &ctx,
		"GED Obol group draw-intent mode",
		fallback_path, draw_mode,
		ged_draw_scene_ref_draw_intent_is_overlay(ref));
    }

    struct bsg_draw_intent *di = bsg_scene_draw_intent(ref);
    if (di) {
	bsg_draw_intent_set_mode(di, (bsg_draw_mode)draw_mode);
	return 1;
    }

    return ged_draw_scene_ref_ensure_draw_intent(ref, fallback_path, draw_mode);
}


static int
ged_draw_scene_ref_draw_intent_is_overlay(bsg_scene_ref ref)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    const struct bsg_draw_intent *di = bsg_scene_draw_intent(ref);
    return bsg_draw_intent_is_overlay(di);
}


static int
ged_draw_scene_ref_draw_intent_appearance(bsg_scene_ref ref,
					  struct bsg_appearance_settings *settings)
{
    if (bsg_scene_ref_is_null(ref) || !settings)
	return 0;

    const struct bsg_draw_intent *di = bsg_scene_draw_intent(ref);
    return bsg_draw_intent_appearance(di, settings);
}


static void
_ged_draw_appearance_to_neutral(struct ged_draw_appearance_settings *out,
				const struct bsg_appearance_settings *in)
{
    if (!out)
	return;

    struct ged_draw_appearance_settings settings =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    if (in) {
	settings.draw_mode = in->draw_mode;
	settings.mixed_modes = in->mixed_modes;
	settings.transparency = in->transparency;
	settings.color_override = in->color_override;
	settings.color[0] = in->color[0];
	settings.color[1] = in->color[1];
	settings.color[2] = in->color[2];
	settings.s_line_width = in->s_line_width;
	settings.s_arrow_tip_length = in->s_arrow_tip_length;
	settings.s_arrow_tip_width = in->s_arrow_tip_width;
	settings.draw_solid_lines_only = in->draw_solid_lines_only;
	settings.draw_non_subtract_only = in->draw_non_subtract_only;
	settings.strict_fallback = in->strict_fallback;
    }

    *out = settings;
}


static int
ged_draw_scene_ref_draw_intent_appearance_settings(
	bsg_scene_ref ref,
	struct ged_draw_appearance_settings *settings)
{
    if (!settings)
	return 0;

    *settings = (struct ged_draw_appearance_settings)
	GED_DRAW_APPEARANCE_SETTINGS_INIT;

    struct bsg_appearance_settings bsg_settings = BSG_APPEARANCE_SETTINGS_INIT;
    if (!ged_draw_scene_ref_draw_intent_appearance(ref, &bsg_settings))
	return 0;

    _ged_draw_appearance_to_neutral(settings, &bsg_settings);
    return 1;
}


static const char *
ged_draw_scene_ref_draw_intent_path(bsg_scene_ref ref)
{
    if (bsg_scene_ref_is_null(ref))
	return NULL;

    const struct bsg_draw_intent *di = bsg_scene_draw_intent(ref);
    return bsg_draw_intent_path(di);
}


static int
ged_draw_scene_ref_draw_intent_mode(bsg_scene_ref ref, int *draw_mode)
{
    if (draw_mode)
	*draw_mode = -1;
    if (bsg_scene_ref_is_null(ref) || !draw_mode)
	return 0;

    const struct bsg_draw_intent *di = bsg_scene_draw_intent(ref);
    if (!di)
	return 0;

    *draw_mode = bsg_draw_intent_mode(di);
    return 1;
}


static int
ged_draw_scene_ref_copy_draw_intent(bsg_scene_ref dst, bsg_scene_ref src)
{
    if (bsg_scene_ref_is_null(dst) || bsg_scene_ref_is_null(src))
	return 0;

    const struct bsg_draw_intent *di = bsg_scene_draw_intent(src);
    const char *path = bsg_draw_intent_path(di);
    if (!path || !path[0])
	return 0;

    struct bsg_draw_intent *copy =
	bsg_draw_intent_create(path, bsg_draw_intent_mode(di));
    if (!copy)
	return 0;

    struct bsg_appearance_settings appearance;
    copy->di_lod = bsg_draw_intent_lod(di);
    copy->di_mixed = di->di_mixed;
    if (bsg_draw_intent_appearance(di, &appearance))
	bsg_draw_intent_set_appearance(copy, &appearance);
    bsg_scene_set_draw_intent(dst, copy);
    return 1;
}


static int
ged_draw_scene_ref_line_style(bsg_scene_ref ref)
{
    return bsg_scene_line_style(ref);
}


static int
ged_draw_scene_ref_set_legacy_color_info(bsg_scene_ref ref,
					 const unsigned char basecolor[3],
					 int user_color,
					 int default_color)
{
    if (bsg_scene_ref_is_null(ref) || !basecolor)
	return 0;

    bsg_scene_set_legacy_color_info(ref, basecolor, user_color, default_color);
    return 1;
}


static int
ged_draw_scene_ref_set_legacy_eval_flag(bsg_scene_ref ref,
					int evaluated_region)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    struct _ged_draw_obol_set_evaluated_region_ctx ctx = {
	evaluated_region ? 1 : 0
    };
    if (ged_draw_scene_ref_obol_source_path_apply_ensure(ref,
	    _ged_draw_obol_set_evaluated_region_cb, &ctx))
	return 1;

    bsg_scene_set_legacy_eval_flag(ref, evaluated_region ? 1 : 0);
    return 1;
}


static int
ged_draw_scene_ref_set_legacy_region_id(bsg_scene_ref ref, int region_id)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    struct _ged_draw_obol_set_region_metadata_ctx ctx = {
	region_id,
	0,
	0,
	0
    };
    if (ged_draw_scene_ref_obol_source_path_apply_ensure(ref,
	    _ged_draw_obol_set_region_metadata_cb, &ctx))
	return 1;

    bsg_scene_set_legacy_region_id(ref, region_id);
    return 1;
}


static int
_ged_draw_scene_ref_apply_db_object_marker(bsg_scene_ref ref)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    bsg_scene_mark_db_object(ref);
    return 1;
}


static int
ged_draw_scene_ref_set_non_database_source(bsg_scene_ref ref,
					   int non_database_source)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    bsg_scene_set_non_database_source(ref, non_database_source ? 1 : 0);
    return 1;
}


static struct bsg_draw_ctx *
ged_drawable_draw_bookkeeping_ctx(struct ged_drawable *gdp)
{
    if (!gdp)
	return NULL;

    if (!gdp->gd_draw_bookkeeping_ctx) {
	struct bsg_draw_ctx *ctx;
	BU_GET(ctx, struct bsg_draw_ctx);
	memset(ctx, 0, sizeof(*ctx));
	gdp->gd_draw_bookkeeping_ctx = (void *)ctx;
    }

    return (struct bsg_draw_ctx *)gdp->gd_draw_bookkeeping_ctx;
}


void
ged_draw_scene_context_free_draw_bookkeeping(struct ged *gedp)
{
    if (!gedp || !gedp->i || !gedp->i->ged_gdp ||
	    !gedp->i->ged_gdp->gd_draw_bookkeeping_ctx)
	return;

    struct bsg_draw_ctx *ctx =
	(struct bsg_draw_ctx *)gedp->i->ged_gdp->gd_draw_bookkeeping_ctx;
    gedp->i->ged_gdp->gd_draw_bookkeeping_ctx = NULL;
    BU_PUT(ctx, struct bsg_draw_ctx);
}


static int
ged_draw_scene_rt_ref_attach_draw_bookkeeping(struct ged *gedp,
					      rt_view_scene_ref scene_ref)
{
    bsg_scene_ref root_ref = ged_draw_scene_ref_from_rt_view_ref(scene_ref);
    if (!gedp || !gedp->i || !gedp->i->ged_gdp ||
	    bsg_scene_ref_is_null(root_ref))
	return 0;

    struct ged_drawable *gdp = gedp->i->ged_gdp;
    struct bsg_draw_ctx *ctx = ged_drawable_draw_bookkeeping_ctx(gdp);
    if (!ctx)
	return 0;

    ctx->draw_rev = &gdp->gd_draw_rev;
    ctx->fso = rt_view_set_context_recycle_pool(
	    ged_view_set_ctx(gedp));
    ctx->owner_data = gedp;
    bsg_scene_set_draw_ctx(root_ref, ctx);
    return 1;
}


ged_draw_group_ref
ged_draw_source_root_create_group_ref(struct ged *gedp, void *view_ctx)
{
    if (!gedp || !view_ctx)
	return GED_DRAW_GROUP_REF_NULL;

    rt_view_scene_ref root_rt_ref =
	ged_draw_view_context_group_create_ref(view_ctx, "_draw_root");
    if (rt_view_scene_ref_is_null(root_rt_ref))
	return GED_DRAW_GROUP_REF_NULL;

    ged_draw_group_ref root_group_ref =
	ged_draw_registry_group_ref_from_source_ref(gedp, root_rt_ref);
    if (ged_draw_group_ref_is_null(root_group_ref)) {
	ged_draw_scene_ref_release(
		ged_draw_scene_ref_from_rt_view_ref(root_rt_ref));
	return GED_DRAW_GROUP_REF_NULL;
    }
    root_group_ref.scene_revision = 0;
    (void)ged_draw_scene_rt_ref_set_visible(root_rt_ref, 1);
    ged_scene_root_group_ref_set(gedp, root_group_ref);
    (void)ged_draw_scene_rt_ref_attach_draw_bookkeeping(gedp, root_rt_ref);

    return root_group_ref;
}


static int
ged_draw_source_view_context_scene_root_set(void *view_ctx,
					   rt_view_scene_ref root)
{
    return rt_view_context_scene_root_ref_attach(view_ctx, root);
}


static rt_view_scene_ref
ged_draw_obol_source_root_rt_ref(struct ged *gedp)
{
    if (!gedp || !ged_draw_obol_scene_controller_attached(gedp))
	return rt_view_scene_ref_null();

    void *root_ctx = ged_draw_obol_context_token_for_path(gedp, "/", NULL);
    if (!root_ctx)
	return rt_view_scene_ref_null();

    return rt_view_scene_ref_make(root_ctx, RT_VIEW_SCENE_BACKEND_OBOL);
}


void *
ged_draw_source_root_attach_view_contexts(struct ged *gedp,
					  void *active_view_ctx,
					  struct bu_ptbl *views)
{
    if (!gedp)
	return NULL;

    rt_view_scene_ref obol_root_rt_ref =
	ged_draw_obol_source_root_rt_ref(gedp);
    ged_draw_group_ref root_group_ref = ged_scene_root_group_ref(gedp);
    rt_view_scene_ref root_rt_ref =
	ged_draw_registry_group_ref_rt_ref(gedp, root_group_ref);
    if (rt_view_scene_ref_is_null(root_rt_ref)) {
	if (!active_view_ctx && rt_view_scene_ref_is_null(obol_root_rt_ref))
	    return NULL;

	if (active_view_ctx && rt_view_scene_ref_is_null(obol_root_rt_ref)) {
	    root_group_ref = ged_draw_source_root_create_group_ref(gedp,
		    active_view_ctx);
	    root_rt_ref = ged_draw_registry_group_ref_rt_ref(gedp,
		    root_group_ref);
	    if (rt_view_scene_ref_is_null(root_rt_ref) &&
		    rt_view_scene_ref_is_null(obol_root_rt_ref))
		return NULL;
	}
    }

    rt_view_scene_ref attachable_root = !rt_view_scene_ref_is_null(
	    root_rt_ref) ? root_rt_ref : obol_root_rt_ref;
    rt_view_scene_ref fallback_root = !rt_view_scene_ref_is_null(
	    obol_root_rt_ref) ? obol_root_rt_ref : root_rt_ref;
    int attached = 0;

    if (active_view_ctx)
	attached += ged_draw_source_view_context_scene_root_set(
		active_view_ctx, attachable_root);
    if (views) {
	for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	    void *view_ctx = BU_PTBL_GET(views, i);
	    if (!view_ctx)
		continue;
	    attached += ged_draw_source_view_context_scene_root_set(view_ctx,
		    attachable_root);
	}
    }

    if (attached > 0)
	return rt_view_scene_ref_context(attachable_root);

    if (rt_view_scene_ref_is_null(fallback_root))
	return NULL;
    return rt_view_scene_ref_context(fallback_root);
}


static void
ged_draw_scene_ref_set_material_rgb(bsg_scene_ref ref,
				    const unsigned char rgb[3])
{
    if (bsg_scene_ref_is_null(ref) || !rgb)
	return;

    if (ged_draw_scene_ref_obol_update_display(ref,
	    0, 0,
	    0, 0,
	    0, 0,
	    0, 0,
	    0, 0,
	    0, 0.0,
	    0, NULL,
	    1, rgb,
	    0, 0))
	return;

    bsg_scene_material_set_rgb(ref, rgb[0], rgb[1], rgb[2]);
}


static int
ged_draw_scene_ref_set_material_revision(bsg_scene_ref ref,
					 uint64_t mater_rev)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    if (ged_draw_scene_ref_obol_update_display(ref,
	    0, 0,
	    0, 0,
	    0, 0,
	    0, 0,
	    0, 0,
	    0, 0.0,
	    0, NULL,
	    0, NULL,
	    1, mater_rev))
	return 1;

    bsg_scene_material_set_revision(ref, (uint32_t)mater_rev);
    return 1;
}


static void
ged_draw_scene_ref_mark_database_source_redraw_result(bsg_scene_ref ref,
						      int failed)
{
    if (ged_draw_scene_ref_is_null(ref))
	return;

    bsg_scene_ref source_ref = ged_draw_scene_ref_source_owner(ref);
    struct ged_draw_database_source_record record;
    memset(&record, 0, sizeof(record));
    if (!ged_draw_scene_ref_database_source_record_get(source_ref, &record))
	return;

    record.stale_reason = failed ? GED_DRAW_STALE_UPDATE_FAILED :
	GED_DRAW_STALE_NONE;
    if (!failed) {
	record.realized_source_revision = record.source_revision;
	record.realized_inputs_revision = record.inputs_revision;
    }

    (void)ged_draw_scene_ref_database_source_record_apply(source_ref,
	    &record);
    if (!ged_draw_scene_ref_equal(source_ref, ref))
	(void)ged_draw_scene_ref_database_source_record_apply(ref, &record);
    (void)ged_draw_scene_ref_obol_set_realization(ref, failed ? 0 : 1,
	    failed ? 1 : 0,
	    record.realized_source_revision,
	    record.realized_inputs_revision,
	    record.stale_reason);
}


static int
ged_draw_scene_ref_realization_pipeline_candidate(bsg_scene_ref ref)
{
    struct ged_draw_obol_realization_policy_summary obol_summary;
    memset(&obol_summary, 0, sizeof(obol_summary));
    if (ged_draw_scene_ref_obol_realization_policy_summary(ref,
	    &obol_summary) && obol_summary.valid)
	return (obol_summary.role_flags &
		(GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_CSG |
		 GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_MESH)) ? 1 : 0;

    struct ged_draw_database_source_record record;
    memset(&record, 0, sizeof(record));
    if (!ged_draw_scene_ref_database_source_record_get(
	    ged_draw_scene_ref_source_owner(ref), &record))
	return 0;
    return (record.realization_role_flags &
	    (GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_CSG |
	     GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_MESH)) ? 1 : 0;
}


static void
ged_draw_scene_ref_invalidate(bsg_scene_ref ref)
{
    bsg_scene_invalidate(ref);
}


static int
_ged_draw_scene_ref_publish_adaptive_wireframe(bsg_scene_ref ref,
					       struct rt_db_internal *ip,
					       const struct bn_tol *tol,
					       const struct rt_view_info *view_info)
{
    if (!ip || !ip->idb_meth)
	return -1;
    if (bsg_scene_ref_is_null(ref) || !ip->idb_meth->ft_lod_realize)
	return -1;

    struct rt_primitive_lod_realization realization;
    struct rt_view_info default_view_info = RT_VIEW_INFO_INIT;
    memset(&realization, 0, sizeof(realization));
    if (!view_info)
	view_info = &default_view_info;
    int ret = ip->idb_meth->ft_lod_realize(&realization, ip, tol, view_info,
	    ged_draw_scene_ref_draw_size(ref));
    if (ret < 0) {
	if (realization.line_points)
	    bu_free(realization.line_points, "primitive LoD line-set points");
	if (realization.line_commands)
	    bu_free(realization.line_commands, "primitive LoD line-set commands");
	return ret;
    }
    if (!realization.has_line_set) {
	if (realization.line_points)
	    bu_free(realization.line_points, "primitive LoD line-set points");
	if (realization.line_commands)
	    bu_free(realization.line_commands, "primitive LoD line-set commands");
	return -1;
    }
    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    if (ok) {
	ged_draw_shape_state *shape_data =
	    _ged_draw_shape_state_get_scene_ref(ref);
	if (shape_data && shape_data->gedp && shape_data->gedp->i &&
		shape_data->gedp->i->ged_gdp)
	    shape_data->gedp->i->ged_gdp->
		gd_draw_retained_primitive_wireframe_publications++;
    }
    if (realization.line_points)
	bu_free(realization.line_points, "primitive LoD line-set points");
    if (realization.line_commands)
	bu_free(realization.line_commands, "primitive LoD line-set commands");
    return ok ? ret : -1;
}


static void
ged_draw_primitive_realization_line_set_free(struct rt_primitive_lod_realization *realization)
{
    if (!realization)
	return;
    if (realization->line_points)
	bu_free(realization->line_points, "primitive realization line-set points");
    if (realization->line_commands)
	bu_free(realization->line_commands, "primitive realization line-set commands");
    realization->line_points = NULL;
    realization->line_commands = NULL;
    realization->line_count = 0;
    realization->line_capacity = 0;
    realization->has_line_set = 0;
}


static int
ged_draw_obol_database_source_publish_realization_line_set(
	struct ged *gedp,
	const char *path,
	struct rt_primitive_lod_realization *realization)
{
    int ok = 0;

    if (!gedp || !path || !path[0] || !realization ||
	    !realization->has_line_set)
	return 0;

    ok = realization->line_count ?
	ged_draw_obol_database_source_publish_line_set_for_path(gedp, path,
		(const point_t *)realization->line_points,
		realization->line_commands, realization->line_count) :
	ged_draw_obol_database_source_clear_vlist_for_path(gedp, path);

    ged_draw_primitive_realization_line_set_free(realization);
    return ok;
}


struct ged_draw_scene_ref_obol_realization_line_set {
    struct rt_primitive_lod_realization *realization;
};


static int
ged_draw_scene_ref_obol_realization_line_set_cb(
	struct ged *gedp,
	const char *path,
	void *data)
{
    struct ged_draw_scene_ref_obol_realization_line_set *ctx =
	(struct ged_draw_scene_ref_obol_realization_line_set *)data;

    if (!ctx || !ctx->realization || !ctx->realization->has_line_set)
	return 0;

    return ctx->realization->line_count ?
	ged_draw_obol_database_source_publish_line_set_for_path(gedp, path,
		(const point_t *)ctx->realization->line_points,
		ctx->realization->line_commands,
		ctx->realization->line_count) :
	ged_draw_obol_database_source_clear_vlist_for_path(gedp, path);
}


static int
ged_draw_scene_ref_obol_publish_realization_line_set(
	bsg_scene_ref ref,
	struct rt_primitive_lod_realization *realization)
{
    if (!realization || !realization->has_line_set)
	return 0;

    struct ged_draw_scene_ref_obol_realization_line_set ctx;
    ctx.realization = realization;
    (void)ged_draw_scene_ref_obol_sync_source_placement(ref);
    if (!ged_draw_scene_ref_obol_source_path_apply(ref,
	    ged_draw_scene_ref_obol_realization_line_set_cb, &ctx))
	return 0;

    ged_draw_primitive_realization_line_set_free(realization);
    return 1;
}


static int
ged_draw_scene_ref_publish_realization_line_set(bsg_scene_ref ref,
						struct rt_primitive_lod_realization *realization)
{
    int ok = 0;

    if (!realization || !realization->has_line_set)
	return 0;

    ok = realization->line_count ?
	ged_draw_scene_ref_publish_line_set(ref,
		(const point_t *)realization->line_points,
		realization->line_commands, realization->line_count) :
	ged_draw_scene_ref_geometry_clear(ref);

    ged_draw_primitive_realization_line_set_free(realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_brep_wireframe_line_set(bsg_scene_ref ref,
						   const struct rt_brep_internal *bi,
						   const struct bn_tol *tol,
						   int *obol_published)
{
    struct rt_primitive_lod_realization realization = { 0 };
    struct rt_db_internal brep_ip;
    int ret = 0;

    if (!bi)
	return 0;

    RT_BREP_CK_MAGIC(bi);
    if (!bi->brep)
	return 0;

    RT_DB_INTERNAL_INIT(&brep_ip);
    brep_ip.idb_major_type = DB5_MAJORTYPE_BRLCAD;
    brep_ip.idb_type = ID_BREP;
    brep_ip.idb_meth = &OBJ[ID_BREP];
    brep_ip.idb_ptr = (void *)bi;

    ret = rt_brep_wireframe_line_set(&realization, &brep_ip, tol);
    if (ret < 0) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }
    if (!realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    if (ged_draw_scene_ref_obol_publish_realization_line_set(ref,
	    &realization)) {
	if (obol_published)
	    *obol_published = 1;
	return 1;
    }

    return ged_draw_scene_ref_publish_realization_line_set(ref, &realization);
}


static int
ged_draw_scene_ref_publish_bspline_wireframe_line_set(bsg_scene_ref ref,
						      struct rt_db_internal *ip,
						      const struct bn_tol *tol,
						      int *obol_published)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = 0;

    if (!ip || ip->idb_type != ID_BSPLINE || !ip->idb_ptr)
	return 0;

    ret = rt_nurb_wireframe_line_set(&realization, ip, tol);
    if (ret < 0) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }
    if (!realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    if (ged_draw_scene_ref_obol_publish_realization_line_set(ref,
	    &realization)) {
	if (obol_published)
	    *obol_published = 1;
	return 1;
    }

    return ged_draw_scene_ref_publish_realization_line_set(ref, &realization);
}


struct ged_draw_line_set_buffer {
    point_t *points;
    int *commands;
    size_t count;
    size_t capacity;
};


static void
ged_draw_line_set_buffer_free(struct ged_draw_line_set_buffer *buffer)
{
    if (!buffer)
	return;
    if (buffer->points)
	bu_free(buffer->points, "GED draw dynamic line-set points");
    if (buffer->commands)
	bu_free(buffer->commands, "GED draw dynamic line-set commands");
    buffer->points = NULL;
    buffer->commands = NULL;
    buffer->count = 0;
    buffer->capacity = 0;
}


static void
ged_draw_line_set_buffer_reserve(struct ged_draw_line_set_buffer *buffer,
				 size_t count)
{
    size_t capacity;

    if (!buffer || count <= buffer->capacity)
	return;

    capacity = buffer->capacity ? buffer->capacity : 64;
    while (capacity < count)
	capacity *= 2;

    buffer->points = (point_t *)bu_realloc(buffer->points,
	    capacity * sizeof(point_t), "GED draw dynamic line-set points");
    buffer->commands = (int *)bu_realloc(buffer->commands,
	    capacity * sizeof(int), "GED draw dynamic line-set commands");
    buffer->capacity = capacity;
}


static void
ged_draw_line_set_buffer_append(struct ged_draw_line_set_buffer *buffer,
				const point_t pt,
				int command)
{
    if (!buffer)
	return;
    ged_draw_line_set_buffer_reserve(buffer, buffer->count + 1);
    VMOVE(buffer->points[buffer->count], pt);
    buffer->commands[buffer->count] = command;
    buffer->count++;
}


static int
ged_draw_scene_ref_publish_grip_wireframe_line_set(bsg_scene_ref ref,
						   const struct rt_grip_internal *grip)
{
    vect_t normal;
    point_t center;
    vect_t xbase;
    vect_t ybase;
    vect_t x_1;
    vect_t x_2;
    vect_t y_1;
    vect_t y_2;
    vect_t tip;
    point_t points[11];
    int commands[11];
    size_t idx = 0;

    if (!grip)
	return 0;
    RT_GRIP_CK_MAGIC(grip);

    VMOVE(center, grip->center);
    VMOVE(normal, grip->normal);

    bn_vec_perp(xbase, normal);
    VCROSS(ybase, xbase, normal);

    VUNITIZE(xbase);
    VUNITIZE(ybase);
    VSCALE(xbase, xbase, grip->mag / 4.0);
    VSCALE(ybase, ybase, grip->mag / 4.0);

    VADD2(x_1, center, xbase);
    VSUB2(x_2, center, xbase);
    VADD2(y_1, center, ybase);
    VSUB2(y_2, center, ybase);

    VMOVE(points[idx], x_1);
    commands[idx++] = BSG_GEOMETRY_LINE_MOVE;
    VMOVE(points[idx], y_1);
    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
    VMOVE(points[idx], x_2);
    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
    VMOVE(points[idx], y_2);
    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
    VMOVE(points[idx], x_1);
    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;

    VSCALE(tip, normal, grip->mag);
    VADD2(tip, center, tip);

    VMOVE(points[idx], x_1);
    commands[idx++] = BSG_GEOMETRY_LINE_MOVE;
    VMOVE(points[idx], tip);
    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
    VMOVE(points[idx], x_2);
    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
    VMOVE(points[idx], y_1);
    commands[idx++] = BSG_GEOMETRY_LINE_MOVE;
    VMOVE(points[idx], tip);
    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
    VMOVE(points[idx], y_2);
    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;

    return ged_draw_scene_ref_publish_line_set(ref, (const point_t *)points,
	    commands, idx);
}


static void
ged_draw_line_set_append_16_point_ring(point_t *points,
				       int *commands,
				       size_t *idx,
				       const fastf_t *ring)
{
    VMOVE(points[*idx], &ring[15 * ELEMENTS_PER_VECT]);
    commands[(*idx)++] = BSG_GEOMETRY_LINE_MOVE;
    for (int i = 0; i < 16; i++) {
	VMOVE(points[*idx], &ring[i * ELEMENTS_PER_VECT]);
	commands[(*idx)++] = BSG_GEOMETRY_LINE_DRAW;
    }
}


static void
ged_draw_line_set_append_16_point_connectors(point_t *points,
					     int *commands,
					     size_t *idx,
					     const fastf_t *top,
					     const fastf_t *bottom)
{
    for (int i = 0; i < 16; i += 4) {
	VMOVE(points[*idx], &top[i * ELEMENTS_PER_VECT]);
	commands[(*idx)++] = BSG_GEOMETRY_LINE_MOVE;
	VMOVE(points[*idx], &bottom[i * ELEMENTS_PER_VECT]);
	commands[(*idx)++] = BSG_GEOMETRY_LINE_DRAW;
    }
}


static int
ged_draw_scene_ref_publish_half_wireframe_line_set(bsg_scene_ref ref,
						   const struct rt_half_internal *half)
{
    plane_t eqn;
    vect_t cent;
    vect_t xbase;
    vect_t ybase;
    vect_t x_1;
    vect_t x_2;
    vect_t y_1;
    vect_t y_2;
    vect_t tip;
    point_t points[10];
    int commands[10];
    size_t idx = 0;

    if (!half)
	return 0;
    RT_HALF_CK_MAGIC(half);

    HMOVE(eqn, half->eqn);
    VSCALE(cent, eqn, eqn[W]);

    bn_vec_perp(xbase, &eqn[0]);
    VCROSS(ybase, xbase, eqn);

    VUNITIZE(xbase);
    VUNITIZE(ybase);
    VSCALE(xbase, xbase, 1000);
    VSCALE(ybase, ybase, 1000);

    VADD2(x_1, cent, xbase);
    VSUB2(x_2, cent, xbase);
    VADD2(y_1, cent, ybase);
    VSUB2(y_2, cent, ybase);

    VMOVE(points[idx], x_1);
    commands[idx++] = BSG_GEOMETRY_LINE_MOVE;
    VMOVE(points[idx], x_2);
    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
    VMOVE(points[idx], y_1);
    commands[idx++] = BSG_GEOMETRY_LINE_MOVE;
    VMOVE(points[idx], y_2);
    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
    VMOVE(points[idx], x_2);
    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
    VMOVE(points[idx], y_1);
    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
    VMOVE(points[idx], x_1);
    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
    VMOVE(points[idx], y_2);
    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;

    VSCALE(tip, eqn, 500);
    VADD2(tip, cent, tip);
    VMOVE(points[idx], cent);
    commands[idx++] = BSG_GEOMETRY_LINE_MOVE;
    VMOVE(points[idx], tip);
    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;

    return ged_draw_scene_ref_publish_line_set(ref, (const point_t *)points,
	    commands, idx);
}


static int
ged_draw_scene_ref_publish_arb_wireframe_line_set(bsg_scene_ref ref,
						  const struct rt_arb_internal *arb)
{
    static const int contours[4][4] = {
	{0, 1, 2, 3},
	{4, 0, 3, 7},
	{5, 4, 7, 6},
	{1, 5, 6, 2}
    };
    point_t points[16];
    int commands[16];
    size_t idx = 0;

    if (!arb)
	return 0;
    RT_ARB_CK_MAGIC(arb);

    for (int i = 0; i < 4; i++) {
	for (int j = 0; j < 4; j++) {
	    VMOVE(points[idx], arb->pt[contours[i][j]]);
	    commands[idx++] = (j == 0) ? BSG_GEOMETRY_LINE_MOVE :
		BSG_GEOMETRY_LINE_DRAW;
	}
    }

    return ged_draw_scene_ref_publish_line_set(ref, (const point_t *)points,
	    commands, idx);
}


static int
ged_draw_scene_ref_publish_arbn_wireframe_line_set(bsg_scene_ref ref,
						   const struct rt_arbn_internal *aip,
						   const struct bn_tol *tol)
{
    struct ged_draw_line_set_buffer buffer = { NULL, NULL, 0, 0 };
    int ok;

    if (!aip || !tol)
	return 0;
    RT_ARBN_CK_MAGIC(aip);

    for (size_t i = 0; i < aip->neqn - 1; i++) {
	for (size_t j = i + 1; j < aip->neqn; j++) {
	    double dot;
	    int point_count;
	    point_t a;
	    point_t b;
	    vect_t dist;

	    VSETALL(a, 0);
	    VSETALL(b, 0);

	    dot = VDOT(aip->eqn[i], aip->eqn[j]);
	    if (BN_VECT_ARE_PARALLEL(dot, tol))
		continue;

	    point_count = 0;
	    for (size_t k = 0; k < aip->neqn; k++) {
		size_t m;
		point_t pt;
		size_t next_k = 0;

		if (k == i || k == j)
		    continue;
		if (bg_make_pnt_3planes(pt, aip->eqn[i], aip->eqn[j], aip->eqn[k]) < 0)
		    continue;

		for (m = 0; m < aip->neqn; m++) {
		    if (i == m || j == m || k == m)
			continue;
		    if (VDOT(pt, aip->eqn[m])-aip->eqn[m][3] > tol->dist) {
			next_k = 1;
			break;
		    }
		}

		if (next_k != 0)
		    continue;

		if (point_count <= 0) {
		    ged_draw_line_set_buffer_append(&buffer, pt,
			    BSG_GEOMETRY_LINE_MOVE);
		    VMOVE(a, pt);
		} else if (point_count == 1) {
		    VSUB2(dist, pt, a);
		    if (MAGSQ(dist) < tol->dist_sq)
			continue;
		    ged_draw_line_set_buffer_append(&buffer, pt,
			    BSG_GEOMETRY_LINE_DRAW);
		    VMOVE(b, pt);
		} else {
		    VSUB2(dist, pt, a);
		    if (MAGSQ(dist) < tol->dist_sq)
			continue;
		    VSUB2(dist, pt, b);
		    if (MAGSQ(dist) < tol->dist_sq)
			continue;
		    bu_log("ARBN typed wireframe error, point_count=%d (>2) on edge %zu/%zu, non-convex\n",
			    point_count + 1, i, j);
		    VPRINT(" a", a);
		    VPRINT(" b", b);
		    VPRINT("pt", pt);
		    ged_draw_line_set_buffer_append(&buffer, pt,
			    BSG_GEOMETRY_LINE_DRAW);
		}
		point_count++;
	    }
	}
    }

    ok = ged_draw_scene_ref_publish_line_set(ref, (const point_t *)buffer.points,
	    buffer.commands, buffer.count);
    ged_draw_line_set_buffer_free(&buffer);
    return ok;
}


static int
ged_draw_scene_ref_publish_ars_wireframe_line_set(bsg_scene_ref ref,
						  const struct rt_ars_internal *ars)
{
    struct ged_draw_line_set_buffer buffer = { NULL, NULL, 0, 0 };
    int ok;

    if (!ars)
	return 0;
    RT_ARS_CK_MAGIC(ars);

    for (size_t i = 0; i < ars->ncurves; i++) {
	fastf_t *v1 = ars->curves[i];

	ged_draw_line_set_buffer_append(&buffer, v1,
		BSG_GEOMETRY_LINE_MOVE);
	v1 += ELEMENTS_PER_VECT;
	for (size_t j = 1; j <= ars->pts_per_curve; j++, v1 += ELEMENTS_PER_VECT)
	    ged_draw_line_set_buffer_append(&buffer, v1,
		    BSG_GEOMETRY_LINE_DRAW);
    }

    for (size_t i = 0; i < ars->pts_per_curve; i++) {
	ged_draw_line_set_buffer_append(&buffer,
		&ars->curves[0][i * ELEMENTS_PER_VECT],
		BSG_GEOMETRY_LINE_MOVE);
	for (size_t j = 1; j < ars->ncurves; j++)
	    ged_draw_line_set_buffer_append(&buffer,
		    &ars->curves[j][i * ELEMENTS_PER_VECT],
		    BSG_GEOMETRY_LINE_DRAW);
    }

    ok = ged_draw_scene_ref_publish_line_set(ref, (const point_t *)buffer.points,
	    buffer.commands, buffer.count);
    ged_draw_line_set_buffer_free(&buffer);
    return ok;
}


static int
ged_draw_scene_ref_publish_tgc_wireframe_line_set(bsg_scene_ref ref,
						  const struct rt_tgc_internal *tgc)
{
    fastf_t top[16 * ELEMENTS_PER_VECT];
    fastf_t bottom[16 * ELEMENTS_PER_VECT];
    point_t v;
    point_t top_center;
    vect_t a;
    vect_t b;
    vect_t c;
    vect_t d;
    vect_t h;
    point_t points[42];
    int commands[42];
    size_t idx = 0;

    if (!tgc)
	return 0;
    RT_TGC_CK_MAGIC(tgc);

    VMOVE(v, tgc->v);
    VMOVE(a, tgc->a);
    VMOVE(b, tgc->b);
    VMOVE(c, tgc->c);
    VMOVE(d, tgc->d);
    VMOVE(h, tgc->h);

    rt_ell_16pnts(bottom, v, a, b);
    VADD2(top_center, v, h);
    rt_ell_16pnts(top, top_center, c, d);

    ged_draw_line_set_append_16_point_ring(points, commands, &idx, top);
    ged_draw_line_set_append_16_point_ring(points, commands, &idx, bottom);
    ged_draw_line_set_append_16_point_connectors(points, commands, &idx,
	    top, bottom);

    return ged_draw_scene_ref_publish_line_set(ref, (const point_t *)points,
	    commands, idx);
}


static int
ged_draw_scene_ref_publish_cline_wireframe_line_set(bsg_scene_ref ref,
						    const struct rt_cline_internal *cline)
{
    fastf_t top[16 * ELEMENTS_PER_VECT];
    fastf_t bottom[16 * ELEMENTS_PER_VECT];
    point_t center;
    vect_t h;
    point_t top_center;
    vect_t unit_a;
    vect_t unit_b;
    vect_t a;
    vect_t b;
    fastf_t inner_radius;
    point_t points[84];
    int commands[84];
    size_t idx = 0;

    if (!cline)
	return 0;
    RT_CLINE_CK_MAGIC(cline);

    VMOVE(center, cline->v);
    VMOVE(h, cline->h);
    VADD2(top_center, center, h);
    bn_vec_ortho(unit_a, h);
    VCROSS(unit_b, unit_a, h);
    VUNITIZE(unit_b);
    VSCALE(a, unit_a, cline->radius);
    VSCALE(b, unit_b, cline->radius);

    rt_ell_16pnts(bottom, center, a, b);
    rt_ell_16pnts(top, top_center, a, b);
    ged_draw_line_set_append_16_point_ring(points, commands, &idx, top);
    ged_draw_line_set_append_16_point_ring(points, commands, &idx, bottom);
    ged_draw_line_set_append_16_point_connectors(points, commands, &idx,
	    top, bottom);

    if (cline->thickness > 0.0 && cline->thickness < cline->radius) {
	inner_radius = cline->radius - cline->thickness;
	VSCALE(a, unit_a, inner_radius);
	VSCALE(b, unit_b, inner_radius);

	rt_ell_16pnts(bottom, center, a, b);
	rt_ell_16pnts(top, top_center, a, b);
	ged_draw_line_set_append_16_point_ring(points, commands, &idx, top);
	ged_draw_line_set_append_16_point_ring(points, commands, &idx, bottom);
	ged_draw_line_set_append_16_point_connectors(points, commands, &idx,
		top, bottom);
    }

    return ged_draw_scene_ref_publish_line_set(ref, (const point_t *)points,
	    commands, idx);
}


static void
ged_draw_part_hemisphere(point_t ov[13],
			 const fastf_t *v,
			 const fastf_t *a,
			 const fastf_t *b,
			 const fastf_t *h)
{
    VADD2(ov[12], v, h);

    VADD2(ov[0], v, a);
    VJOIN2(ov[1], v, M_SQRT1_2, a, M_SQRT1_2, b);
    VADD2(ov[2], v, b);
    VJOIN2(ov[3], v, -M_SQRT1_2, a, M_SQRT1_2, b);
    VSUB2(ov[4], v, a);
    VJOIN2(ov[5], v, -M_SQRT1_2, a, -M_SQRT1_2, b);
    VSUB2(ov[6], v, b);
    VJOIN2(ov[7], v, M_SQRT1_2, a, -M_SQRT1_2, b);

    VJOIN2(ov[8], v, M_SQRT1_2, a, M_SQRT1_2, h);
    VJOIN2(ov[10], v, -M_SQRT1_2, a, M_SQRT1_2, h);

    VJOIN2(ov[9], v, M_SQRT1_2, b, M_SQRT1_2, h);
    VJOIN2(ov[11], v, -M_SQRT1_2, b, M_SQRT1_2, h);
}


static void
ged_draw_line_set_append_point(point_t *points,
			       int *commands,
			       size_t *idx,
			       const point_t pt,
			       int command)
{
    VMOVE(points[*idx], pt);
    commands[(*idx)++] = command;
}


static void
ged_draw_line_set_append_part_hemisphere_outline(point_t *points,
						 int *commands,
						 size_t *idx,
						 const point_t hemi[13])
{
    ged_draw_line_set_append_point(points, commands, idx, hemi[0],
	    BSG_GEOMETRY_LINE_MOVE);
    for (int i = 7; i >= 0; i--)
	ged_draw_line_set_append_point(points, commands, idx, hemi[i],
		BSG_GEOMETRY_LINE_DRAW);
    ged_draw_line_set_append_point(points, commands, idx, hemi[8],
	    BSG_GEOMETRY_LINE_DRAW);
    ged_draw_line_set_append_point(points, commands, idx, hemi[12],
	    BSG_GEOMETRY_LINE_DRAW);
    ged_draw_line_set_append_point(points, commands, idx, hemi[10],
	    BSG_GEOMETRY_LINE_DRAW);
    ged_draw_line_set_append_point(points, commands, idx, hemi[4],
	    BSG_GEOMETRY_LINE_DRAW);
    ged_draw_line_set_append_point(points, commands, idx, hemi[2],
	    BSG_GEOMETRY_LINE_MOVE);
    ged_draw_line_set_append_point(points, commands, idx, hemi[9],
	    BSG_GEOMETRY_LINE_DRAW);
    ged_draw_line_set_append_point(points, commands, idx, hemi[12],
	    BSG_GEOMETRY_LINE_DRAW);
    ged_draw_line_set_append_point(points, commands, idx, hemi[11],
	    BSG_GEOMETRY_LINE_DRAW);
    ged_draw_line_set_append_point(points, commands, idx, hemi[6],
	    BSG_GEOMETRY_LINE_DRAW);
}


static int
ged_draw_scene_ref_publish_part_wireframe_line_set(bsg_scene_ref ref,
						   const struct rt_part_internal *part)
{
    point_t sphere_rim[16];
    point_t vhemi[13];
    point_t hhemi[13];
    point_t v;
    point_t tail;
    vect_t a;
    vect_t b;
    vect_t c;
    vect_t as;
    vect_t bs;
    vect_t hs;
    vect_t hunit;
    point_t points[51];
    int commands[51];
    size_t idx = 0;

    if (!part)
	return 0;
    RT_PART_CK_MAGIC(part);
    VMOVE(v, part->part_V);

    if (part->part_type == RT_PARTICLE_TYPE_SPHERE) {
	VSET(a, part->part_vrad, 0, 0);
	VSET(b, 0, part->part_vrad, 0);
	VSET(c, 0, 0, part->part_vrad);

	rt_ell_16pnts(&sphere_rim[0][X], v, a, b);
	ged_draw_line_set_append_16_point_ring(points, commands, &idx,
		&sphere_rim[0][X]);

	rt_ell_16pnts(&sphere_rim[0][X], v, b, c);
	ged_draw_line_set_append_16_point_ring(points, commands, &idx,
		&sphere_rim[0][X]);

	rt_ell_16pnts(&sphere_rim[0][X], v, a, c);
	ged_draw_line_set_append_16_point_ring(points, commands, &idx,
		&sphere_rim[0][X]);

	return ged_draw_scene_ref_publish_line_set(ref, (const point_t *)points,
		commands, idx);
    }

    VMOVE(hunit, part->part_H);
    VUNITIZE(hunit);
    bn_vec_perp(a, hunit);
    VUNITIZE(a);
    VCROSS(b, hunit, a);
    VUNITIZE(b);

    VSCALE(as, a, part->part_vrad);
    VSCALE(bs, b, part->part_vrad);
    VSCALE(hs, hunit, -part->part_vrad);
    ged_draw_part_hemisphere(vhemi, v, as, bs, hs);

    VADD2(tail, v, part->part_H);
    VSCALE(as, a, part->part_hrad);
    VSCALE(bs, b, part->part_hrad);
    VSCALE(hs, hunit, part->part_hrad);
    ged_draw_part_hemisphere(hhemi, tail, as, bs, hs);

    ged_draw_line_set_append_part_hemisphere_outline(points, commands, &idx,
	    (const point_t *)vhemi);
    ged_draw_line_set_append_part_hemisphere_outline(points, commands, &idx,
	    (const point_t *)hhemi);

    for (int i = 0; i <= 6; i += 2) {
	ged_draw_line_set_append_point(points, commands, &idx, vhemi[i],
		BSG_GEOMETRY_LINE_MOVE);
	ged_draw_line_set_append_point(points, commands, &idx, hhemi[i],
		BSG_GEOMETRY_LINE_DRAW);
    }

    return ged_draw_scene_ref_publish_line_set(ref, (const point_t *)points,
	    commands, idx);
}


static int
ged_draw_scene_ref_publish_joint_wireframe_line_set(bsg_scene_ref ref,
						    const struct rt_joint_internal *joint)
{
    fastf_t top[16 * ELEMENTS_PER_VECT];
    fastf_t middle[16 * ELEMENTS_PER_VECT];
    fastf_t bottom[16 * ELEMENTS_PER_VECT];
    point_t location;
    point_t a = {5, 0, 0};
    point_t b = {0, 5, 0};
    point_t c = {0, 0, 5};
    point_t points[51];
    int commands[51];
    size_t idx = 0;

    if (!joint)
	return 0;
    RT_JOINT_CK_MAGIC(joint);

    VMOVE(location, joint->location);
    rt_ell_16pnts(top, location, a, b);
    rt_ell_16pnts(bottom, location, b, c);
    rt_ell_16pnts(middle, location, a, c);

    ged_draw_line_set_append_16_point_ring(points, commands, &idx, top);
    ged_draw_line_set_append_16_point_ring(points, commands, &idx, bottom);
    ged_draw_line_set_append_16_point_ring(points, commands, &idx, middle);

    return ged_draw_scene_ref_publish_line_set(ref, (const point_t *)points,
	    commands, idx);
}


static void
ged_draw_line_set_append_sphere_rings(point_t *points,
				      int *commands,
				      size_t *idx,
				      const point_t center,
				      fastf_t radius)
{
    fastf_t top[16 * ELEMENTS_PER_VECT];
    fastf_t middle[16 * ELEMENTS_PER_VECT];
    fastf_t bottom[16 * ELEMENTS_PER_VECT];
    point_t v;
    vect_t a;
    vect_t b;
    vect_t c;

    VMOVE(v, center);
    VSET(a, radius, 0, 0);
    VSET(b, 0, radius, 0);
    VSET(c, 0, 0, radius);

    rt_ell_16pnts(top, v, a, b);
    rt_ell_16pnts(bottom, v, b, c);
    rt_ell_16pnts(middle, v, a, c);

    ged_draw_line_set_append_16_point_ring(points, commands, idx, top);
    ged_draw_line_set_append_16_point_ring(points, commands, idx, bottom);
    ged_draw_line_set_append_16_point_ring(points, commands, idx, middle);
}


static int
ged_draw_scene_ref_publish_metaball_wireframe_line_set(bsg_scene_ref ref,
						       const struct rt_metaball_internal *metaball)
{
    const struct wdb_metaball_pnt *mbpt;
    point_t *points = NULL;
    int *commands = NULL;
    size_t point_count = 0;
    size_t idx = 0;
    int ok;

    if (!metaball)
	return 0;
    RT_METABALL_CK_MAGIC(metaball);

    for (BU_LIST_FOR(mbpt, wdb_metaball_pnt, &metaball->metaball_ctrl_head))
	point_count++;

    if (point_count) {
	points = (point_t *)bu_calloc(point_count * 51, sizeof(point_t),
		"GED draw metaball wireframe points");
	commands = (int *)bu_calloc(point_count * 51, sizeof(int),
		"GED draw metaball wireframe commands");

	for (BU_LIST_FOR(mbpt, wdb_metaball_pnt, &metaball->metaball_ctrl_head)) {
	    point_t center;
	    VMOVE(center, mbpt->coord);
	    ged_draw_line_set_append_sphere_rings(points, commands, &idx,
		    center, mbpt->field_strength / metaball->threshold);
	}
    }

    ok = ged_draw_scene_ref_publish_line_set(ref, (const point_t *)points,
	    commands, idx);

    if (points)
	bu_free(points, "GED draw metaball wireframe points");
    if (commands)
	bu_free(commands, "GED draw metaball wireframe commands");

    return ok;
}


static int
ged_draw_scene_ref_publish_hyp_wireframe_line_set(bsg_scene_ref ref,
						  const struct rt_hyp_internal *hyp_in)
{
    vect_t major_axis[8];
    vect_t minor_axis[8];
    vect_t height_axis[7];
    vect_t bunit;
    vect_t ell[16];
    vect_t ribs[16][7];
    point_t hyp_v;
    vect_t hyp_h;
    vect_t hyp_au;
    fastf_t hyp_r1;
    fastf_t hyp_r2;
    fastf_t hyp_c;
    fastf_t scale;
    fastf_t cos22_5 = 0.9238795325112867385;
    fastf_t cos67_5 = 0.3826834323650898373;
    point_t points[231];
    int commands[231];
    size_t idx = 0;

    if (!hyp_in)
	return 0;
    RT_HYP_CK_MAGIC(hyp_in);

    hyp_r1 = hyp_in->hyp_bnr * MAGNITUDE(hyp_in->hyp_A);
    hyp_r2 = hyp_in->hyp_bnr * hyp_in->hyp_b;
    hyp_c = sqrt(4 * MAGSQ(hyp_in->hyp_A) / MAGSQ(hyp_in->hyp_Hi) *
	    (1 - hyp_in->hyp_bnr * hyp_in->hyp_bnr));
    VSCALE(hyp_h, hyp_in->hyp_Hi, 0.5);
    VADD2(hyp_v, hyp_in->hyp_Vi, hyp_h);
    VMOVE(hyp_au, hyp_in->hyp_A);
    VUNITIZE(hyp_au);

    VCROSS(bunit, hyp_h, hyp_au);
    VUNITIZE(bunit);

    VMOVE(height_axis[0], hyp_h);
    VSCALE(height_axis[1], height_axis[0], 0.5);
    VSCALE(height_axis[2], height_axis[0], 0.25);
    VSETALL(height_axis[3], 0);
    VREVERSE(height_axis[4], height_axis[2]);
    VREVERSE(height_axis[5], height_axis[1]);
    VREVERSE(height_axis[6], height_axis[0]);

    for (int i = 0; i < 7; i++) {
	scale = sqrt(MAGSQ(height_axis[i]) * (hyp_c * hyp_c) /
		(hyp_r1 * hyp_r1) + 1);

	VSCALE(major_axis[0], hyp_au, hyp_r1 * scale);
	VSCALE(major_axis[1], major_axis[0], cos22_5);
	VSCALE(major_axis[2], major_axis[0], M_SQRT1_2);
	VSCALE(major_axis[3], major_axis[0], cos67_5);
	VREVERSE(major_axis[4], major_axis[3]);
	VREVERSE(major_axis[5], major_axis[2]);
	VREVERSE(major_axis[6], major_axis[1]);
	VREVERSE(major_axis[7], major_axis[0]);

	VSCALE(minor_axis[0], bunit, hyp_r2 * scale);
	VSCALE(minor_axis[1], minor_axis[0], cos22_5);
	VSCALE(minor_axis[2], minor_axis[0], M_SQRT1_2);
	VSCALE(minor_axis[3], minor_axis[0], cos67_5);
	VREVERSE(minor_axis[4], minor_axis[3]);
	VREVERSE(minor_axis[5], minor_axis[2]);
	VREVERSE(minor_axis[6], minor_axis[1]);
	VREVERSE(minor_axis[7], minor_axis[0]);

	VADD3(ell[ 0], hyp_v, height_axis[i], major_axis[0]);
	VADD4(ell[ 1], hyp_v, height_axis[i], major_axis[1], minor_axis[3]);
	VADD4(ell[ 2], hyp_v, height_axis[i], major_axis[2], minor_axis[2]);
	VADD4(ell[ 3], hyp_v, height_axis[i], major_axis[3], minor_axis[1]);
	VADD3(ell[ 4], hyp_v, height_axis[i], minor_axis[0]);
	VADD4(ell[ 5], hyp_v, height_axis[i], major_axis[4], minor_axis[1]);
	VADD4(ell[ 6], hyp_v, height_axis[i], major_axis[5], minor_axis[2]);
	VADD4(ell[ 7], hyp_v, height_axis[i], major_axis[6], minor_axis[3]);
	VADD3(ell[ 8], hyp_v, height_axis[i], major_axis[7]);
	VADD4(ell[ 9], hyp_v, height_axis[i], major_axis[6], minor_axis[4]);
	VADD4(ell[10], hyp_v, height_axis[i], major_axis[5], minor_axis[5]);
	VADD4(ell[11], hyp_v, height_axis[i], major_axis[4], minor_axis[6]);
	VADD3(ell[12], hyp_v, height_axis[i], minor_axis[7]);
	VADD4(ell[13], hyp_v, height_axis[i], major_axis[3], minor_axis[6]);
	VADD4(ell[14], hyp_v, height_axis[i], major_axis[2], minor_axis[5]);
	VADD4(ell[15], hyp_v, height_axis[i], major_axis[1], minor_axis[4]);

	ged_draw_line_set_append_point(points, commands, &idx, ell[15],
		BSG_GEOMETRY_LINE_MOVE);
	for (int j = 0; j < 16; j++) {
	    ged_draw_line_set_append_point(points, commands, &idx, ell[j],
		    BSG_GEOMETRY_LINE_DRAW);
	    VMOVE(ribs[j][i], ell[j]);
	}
    }

    for (int i = 0; i < 16; i++) {
	ged_draw_line_set_append_point(points, commands, &idx, ribs[i][0],
		BSG_GEOMETRY_LINE_MOVE);
	for (int j = 1; j < 7; j++)
	    ged_draw_line_set_append_point(points, commands, &idx, ribs[i][j],
		    BSG_GEOMETRY_LINE_DRAW);
    }

    return ged_draw_scene_ref_publish_line_set(ref, (const point_t *)points,
	    commands, idx);
}


static int
ged_draw_scene_ref_publish_pnts_wireframe_line_set(bsg_scene_ref ref,
						   const struct rt_pnts_internal *pnts)
{
    point_t *points = NULL;
    int *commands = NULL;
    size_t point_count;
    size_t command_count;
    size_t idx = 0;
    int ok;

    if (!pnts)
	return 0;
    RT_PNTS_CK_MAGIC(pnts);

    if (pnts->count == 0)
	return ged_draw_scene_ref_publish_line_set(ref, NULL, NULL, 0);

    point_count = (size_t)pnts->count;
    command_count = (pnts->scale > 0) ? point_count * 51 : point_count * 6;
    points = (point_t *)bu_calloc(command_count, sizeof(point_t),
	    "GED draw PNTS wireframe points");
    commands = (int *)bu_calloc(command_count, sizeof(int),
	    "GED draw PNTS wireframe commands");

    struct pnt *point = (struct pnt *)pnts->point;
    struct bu_list *head = &point->l;
    if (pnts->scale > 0) {
	for (BU_LIST_FOR(point, pnt, head)) {
	    point_t center;
	    VMOVE(center, point->v);
	    ged_draw_line_set_append_sphere_rings(points, commands, &idx,
		    center, pnts->scale);
	}
    } else {
	for (BU_LIST_FOR(point, pnt, head)) {
	    point_t a;
	    point_t b;
	    double vcoord = 1;
	    double hcoord = 1;

	    VSET(a, point->v[X] - hcoord, point->v[Y], point->v[Z]);
	    VSET(b, point->v[X] + hcoord, point->v[Y], point->v[Z]);
	    ged_draw_line_set_append_point(points, commands, &idx, a,
		    BSG_GEOMETRY_LINE_MOVE);
	    ged_draw_line_set_append_point(points, commands, &idx, b,
		    BSG_GEOMETRY_LINE_DRAW);

	    VSET(a, point->v[X], point->v[Y] - hcoord, point->v[Z]);
	    VSET(b, point->v[X], point->v[Y] + hcoord, point->v[Z]);
	    ged_draw_line_set_append_point(points, commands, &idx, a,
		    BSG_GEOMETRY_LINE_MOVE);
	    ged_draw_line_set_append_point(points, commands, &idx, b,
		    BSG_GEOMETRY_LINE_DRAW);

	    VSET(a, point->v[X], point->v[Y], point->v[Z] - vcoord);
	    VSET(b, point->v[X], point->v[Y], point->v[Z] + vcoord);
	    ged_draw_line_set_append_point(points, commands, &idx, a,
		    BSG_GEOMETRY_LINE_MOVE);
	    ged_draw_line_set_append_point(points, commands, &idx, b,
		    BSG_GEOMETRY_LINE_DRAW);
	}
    }

    ok = ged_draw_scene_ref_publish_line_set(ref, (const point_t *)points,
	    commands, idx);

    bu_free(points, "GED draw PNTS wireframe points");
    bu_free(commands, "GED draw PNTS wireframe commands");

    return ok;
}


static int
ged_draw_scene_ref_publish_script_wireframe_line_set(bsg_scene_ref ref,
						     const struct rt_script_internal *script)
{
    if (!script)
	return 0;
    RT_SCRIPT_CK_MAGIC(script);

    if (bu_vls_addr(&script->s_type))
	bu_log("Script data not found or not specified\n");

    return ged_draw_scene_ref_publish_line_set(ref, NULL, NULL, 0);
}


static int
ged_draw_scene_ref_publish_pipe_wireframe_line_set(bsg_scene_ref ref,
						   struct rt_db_internal *ip)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_pipe_wireframe_line_set(&realization, ip);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_tor_wireframe_line_set(bsg_scene_ref ref,
						  struct rt_db_internal *ip,
						  const struct bg_tess_tol *ttol)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_tor_wireframe_line_set(&realization, ip, ttol);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_rpc_wireframe_line_set(bsg_scene_ref ref,
						  struct rt_db_internal *ip,
						  const struct bg_tess_tol *ttol)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_rpc_wireframe_line_set(&realization, ip, ttol);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_rhc_wireframe_line_set(bsg_scene_ref ref,
						  struct rt_db_internal *ip,
						  const struct bg_tess_tol *ttol)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_rhc_wireframe_line_set(&realization, ip, ttol);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_epa_wireframe_line_set(bsg_scene_ref ref,
						  struct rt_db_internal *ip,
						  const struct bg_tess_tol *ttol)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_epa_wireframe_line_set(&realization, ip, ttol);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_ehy_wireframe_line_set(bsg_scene_ref ref,
						  struct rt_db_internal *ip,
						  const struct bg_tess_tol *ttol)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_ehy_wireframe_line_set(&realization, ip, ttol);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_eto_wireframe_line_set(bsg_scene_ref ref,
						  struct rt_db_internal *ip,
						  const struct bg_tess_tol *ttol)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_eto_wireframe_line_set(&realization, ip, ttol);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_superell_wireframe_line_set(bsg_scene_ref ref,
						       struct rt_db_internal *ip)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_superell_wireframe_line_set(&realization, ip);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_hrt_wireframe_line_set(bsg_scene_ref ref,
						  struct rt_db_internal *ip,
						  const struct bg_tess_tol *ttol)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_hrt_wireframe_line_set(&realization, ip, ttol);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_datum_wireframe_line_set(bsg_scene_ref ref,
						    struct rt_db_internal *ip)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_datum_wireframe_line_set(&realization, ip);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_sketch_wireframe_line_set(bsg_scene_ref ref,
						     struct rt_db_internal *ip,
						     const struct bg_tess_tol *ttol)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_sketch_wireframe_line_set(&realization, ip, ttol);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_extrude_wireframe_line_set(bsg_scene_ref ref,
						      struct rt_db_internal *ip,
						      const struct bg_tess_tol *ttol)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_extrude_wireframe_line_set(&realization, ip, ttol);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_revolve_wireframe_line_set(bsg_scene_ref ref,
						      struct rt_db_internal *ip,
						      const struct bg_tess_tol *ttol)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_revolve_wireframe_line_set(&realization, ip, ttol);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_ebm_wireframe_line_set(bsg_scene_ref ref,
						  struct rt_db_internal *ip)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_ebm_wireframe_line_set(&realization, ip);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_vol_wireframe_line_set(bsg_scene_ref ref,
						  struct rt_db_internal *ip)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_vol_wireframe_line_set(&realization, ip);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_dsp_wireframe_line_set(bsg_scene_ref ref,
						  struct rt_db_internal *ip,
						  const struct bg_tess_tol *ttol)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_dsp_wireframe_line_set(&realization, ip, ttol);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_hf_wireframe_line_set(bsg_scene_ref ref,
						 struct rt_db_internal *ip,
						 const struct bg_tess_tol *ttol)
{
    struct rt_primitive_lod_realization realization = { 0 };
    int ret = rt_hf_wireframe_line_set(&realization, ip, ttol);

    if (ret < 0 || !realization.has_line_set) {
	ged_draw_primitive_realization_line_set_free(&realization);
	return 0;
    }

    int ok = ged_draw_scene_ref_publish_line_set(ref,
	    (const point_t *)realization.line_points,
	    realization.line_commands, realization.line_count);
    ged_draw_primitive_realization_line_set_free(&realization);
    return ok;
}


static int
ged_draw_scene_ref_publish_nmg_wireframe_line_set(bsg_scene_ref ref,
						  struct rt_db_internal *ip)
{
    if (!ip || ip->idb_type != ID_NMG || !ip->idb_ptr)
	return 0;

    return ged_draw_scene_ref_geometry_publish_nmg_model(ref,
	    (const struct model *)ip->idb_ptr, GED_DRAW_NMG_STYLE_VECTOR);
}


static int
ged_draw_scene_ref_publish_ell_wireframe_line_set(bsg_scene_ref ref,
						  const struct rt_ell_internal *ell)
{
    fastf_t top[16 * ELEMENTS_PER_VECT];
    fastf_t middle[16 * ELEMENTS_PER_VECT];
    fastf_t bottom[16 * ELEMENTS_PER_VECT];
    point_t v;
    vect_t a;
    vect_t b;
    vect_t c;
    point_t points[51];
    int commands[51];
    size_t idx = 0;

    if (!ell)
	return 0;
    RT_ELL_CK_MAGIC(ell);

    VMOVE(v, ell->v);
    VMOVE(a, ell->a);
    VMOVE(b, ell->b);
    VMOVE(c, ell->c);

    rt_ell_16pnts(top, v, a, b);
    rt_ell_16pnts(bottom, v, b, c);
    rt_ell_16pnts(middle, v, a, c);

    ged_draw_line_set_append_16_point_ring(points, commands, &idx, top);
    ged_draw_line_set_append_16_point_ring(points, commands, &idx, bottom);
    ged_draw_line_set_append_16_point_ring(points, commands, &idx, middle);

    return ged_draw_scene_ref_publish_line_set(ref, (const point_t *)points,
	    commands, idx);
}

static int
ged_draw_scene_ref_publish_annot_record(bsg_scene_ref ref,
					struct rt_db_internal *ip)
{
    if (!ip || ip->idb_type != ID_ANNOT || !ip->idb_ptr)
	return 0;

    ged_draw_shape_state *shape_data = _ged_draw_shape_state_get_scene_ref(ref);
    if (!shape_data)
	return 0;

    struct rt_annot_internal *ann = (struct rt_annot_internal *)ip->idb_ptr;
    RT_ANNOT_CK_MAGIC(ann);

    point_t *points = NULL;
    point_t *obol_line_points = NULL;
    int *obol_line_commands = NULL;
    size_t obol_line_count = 0;
    struct bsg_annotation_segment *segments = NULL;
    struct ged_draw_obol_annotation_segment *obol_segments = NULL;
    struct bu_vls summary = BU_VLS_INIT_ZERO;
    mat_t model_mat, display_mat;
    int ok = 0;

    MAT_IDN(model_mat);
    MAT_IDN(display_mat);

    if (ann->vert_count) {
	points = (point_t *)bu_calloc(ann->vert_count, sizeof(point_t),
		"GED ANNOT points");
	for (size_t i = 0; i < ann->vert_count; i++)
	    VSET(points[i], ann->verts[i][X], ann->verts[i][Y], 0.0);
    }

    if (ann->ant.count && ann->vert_count && ann->ant.count <= SIZE_MAX / 2) {
	obol_line_points = (point_t *)bu_calloc(ann->ant.count * 2,
		sizeof(point_t), "GED ANNOT Obol line points");
	obol_line_commands = (int *)bu_calloc(ann->ant.count * 2,
		sizeof(int), "GED ANNOT Obol line commands");
    }

    if (ann->ant.count) {
	segments = (struct bsg_annotation_segment *)bu_calloc(ann->ant.count,
		sizeof(struct bsg_annotation_segment), "GED ANNOT segments");
	obol_segments = (struct ged_draw_obol_annotation_segment *)bu_calloc(
		ann->ant.count, sizeof(struct ged_draw_obol_annotation_segment),
		"GED ANNOT Obol segments");
	for (size_t i = 0; i < ann->ant.count; i++) {
	    uint32_t *lng = ann->ant.segments ? (uint32_t *)ann->ant.segments[i] : NULL;
	    if (!lng)
		continue;
	    segments[i].reverse = ann->ant.reverse ? ann->ant.reverse[i] : 0;
	    switch (*lng) {
		case CURVE_LSEG_MAGIC: {
		    const struct line_seg *lsg = (const struct line_seg *)lng;
		    segments[i].kind = BSG_ANNOTATION_SEGMENT_LINE;
		    segments[i].data.line.start = lsg->start;
		    segments[i].data.line.end = lsg->end;
		    obol_segments[i].kind =
			GED_DRAW_OBOL_ANNOTATION_SEGMENT_LINE;
		    obol_segments[i].line_start = lsg->start;
		    obol_segments[i].line_end = lsg->end;
		    if (obol_line_points &&
			    lsg->start >= 0 && lsg->end >= 0 &&
			    (size_t)lsg->start < ann->vert_count &&
			    (size_t)lsg->end < ann->vert_count) {
			VSET(obol_line_points[obol_line_count],
				ann->V[X] + ann->verts[lsg->start][X],
				ann->V[Y] + ann->verts[lsg->start][Y],
				ann->V[Z]);
			obol_line_commands[obol_line_count++] =
			    GED_DRAW_VIEW_LINE_MOVE;
			VSET(obol_line_points[obol_line_count],
				ann->V[X] + ann->verts[lsg->end][X],
				ann->V[Y] + ann->verts[lsg->end][Y],
				ann->V[Z]);
			obol_line_commands[obol_line_count++] =
			    GED_DRAW_VIEW_LINE_DRAW;
		    }
		    break;
		}
		case CURVE_CARC_MAGIC: {
		    const struct carc_seg *csg = (const struct carc_seg *)lng;
		    segments[i].kind = BSG_ANNOTATION_SEGMENT_CARC;
		    segments[i].data.carc.start = csg->start;
		    segments[i].data.carc.end = csg->end;
		    segments[i].data.carc.radius = csg->radius;
		    segments[i].data.carc.center_is_left = csg->center_is_left;
		    segments[i].data.carc.orientation = csg->orientation;
		    segments[i].data.carc.center = csg->center;
		    break;
		}
		case CURVE_NURB_MAGIC: {
		    const struct nurb_seg *nsg = (const struct nurb_seg *)lng;
		    segments[i].kind = BSG_ANNOTATION_SEGMENT_NURB;
		    segments[i].data.nurb.order = nsg->order;
		    segments[i].data.nurb.point_type = nsg->pt_type;
		    segments[i].data.nurb.knot_count = (nsg->k.k_size > 0) ?
			(size_t)nsg->k.k_size : 0;
		    segments[i].data.nurb.knots = nsg->k.knots;
		    segments[i].data.nurb.control_point_count = (nsg->c_size > 0) ?
			(size_t)nsg->c_size : 0;
		    segments[i].data.nurb.control_points = nsg->ctl_points;
		    segments[i].data.nurb.weights = nsg->weights;
		    break;
		}
		case CURVE_BEZIER_MAGIC: {
		    const struct bezier_seg *bsg = (const struct bezier_seg *)lng;
		    segments[i].kind = BSG_ANNOTATION_SEGMENT_BEZIER;
		    segments[i].data.bezier.degree = bsg->degree;
		    segments[i].data.bezier.control_point_count = (bsg->degree >= 0) ?
			(size_t)bsg->degree + 1 : 0;
		    segments[i].data.bezier.control_points = bsg->ctl_points;
		    break;
		}
		case ANN_TSEG_MAGIC: {
		    const struct txt_seg *tsg = (const struct txt_seg *)lng;
		    const char *label = BU_VLS_IS_INITIALIZED(&tsg->label) ?
			bu_vls_cstr(&tsg->label) : "";
		    segments[i].kind = BSG_ANNOTATION_SEGMENT_TEXT;
		    segments[i].data.text.ref_pt = tsg->ref_pt;
		    segments[i].data.text.relative_position = tsg->rel_pos;
		    segments[i].data.text.text = (char *)label;
		    segments[i].data.text.size = tsg->txt_size;
		    segments[i].data.text.rotation = tsg->txt_rot_angle;
		    obol_segments[i].kind =
			GED_DRAW_OBOL_ANNOTATION_SEGMENT_TEXT;
		    obol_segments[i].text_ref_point = tsg->ref_pt;
		    obol_segments[i].text = label;
		    if (label[0]) {
			if (bu_vls_strlen(&summary))
			    bu_vls_strcat(&summary, " ");
			bu_vls_strcat(&summary, label);
		    }
		    break;
		}
		default:
		    break;
	    }
	}
    }

    if (!bu_vls_strlen(&summary))
	bu_vls_strcpy(&summary, "annotation");

    bsg_scene_ref geometry_ref =
	ged_draw_scene_ref_geometry_publication_target(ref);
    if (bsg_scene_ref_is_null(geometry_ref))
	goto cleanup;

    if (obol_segments) {
	struct _ged_draw_obol_publish_annotation_record_ctx ctx = {
	    {ann->V[X], ann->V[Y], ann->V[Z]},
	    (const point_t *)points,
	    ann->vert_count,
	    obol_segments,
	    ann->ant.count,
	    (const point_t *)obol_line_points,
	    obol_line_commands,
	    obol_line_count
	};
	(void)ged_draw_scene_ref_obol_source_path_apply_ensure(ref,
		_ged_draw_obol_publish_annotation_record_cb, &ctx);
    } else if (obol_line_count) {
	ged_draw_scene_ref_obol_publish_annotation_line_set(ref,
		(const point_t *)obol_line_points, obol_line_commands,
		obol_line_count);
    }

    ok = bsg_geometry_ref_set_annotation(bsg_scene_ref_as_geometry(geometry_ref),
	    bu_vls_cstr(&summary), BSG_ANNOTATION_SPACE_DISPLAY, ann->V,
	    model_mat, display_mat, (const point_t *)points, ann->vert_count,
	    segments, ann->ant.count);
    if (ok) {
	shape_data->geometry_command_count = ann->ant.count;
	shape_data->geometry_revision++;
	bsg_scene_invalidate(ref);
    }

cleanup:
    if (points)
	bu_free(points, "GED ANNOT points");
    if (obol_line_points)
	bu_free(obol_line_points, "GED ANNOT Obol line points");
    if (obol_line_commands)
	bu_free(obol_line_commands, "GED ANNOT Obol line commands");
    if (segments)
	bu_free(segments, "GED ANNOT segments");
    if (obol_segments)
	bu_free(obol_segments, "GED ANNOT Obol segments");
    bu_vls_free(&summary);
    return ok;
}


static int
_ged_draw_scene_ref_publish_direct_primitive_wireframe(bsg_scene_ref ref,
						       struct rt_db_internal *ip,
						       const struct bg_tess_tol *ttol,
						       const struct bn_tol *tol,
						       void *view_ctx)
{
    int ok = 0;
    int obol_published = 0;

    if (!ip)
	return 0;
    if (!ged_draw_rt_internal_payload_valid(ip))
	return -1;

    switch (ip->idb_type) {
	case ID_BOT:
	    ok = ged_draw_scene_ref_publish_bot_wireframe_line_set(ref,
		    (const struct rt_bot_internal *)ip->idb_ptr);
	    break;
	case ID_POLY:
	    ok = ged_draw_scene_ref_publish_poly_wireframe_line_set(ref,
		    (const struct rt_pg_internal *)ip->idb_ptr);
	    break;
	case ID_BREP:
	    ok = ged_draw_scene_ref_publish_brep_wireframe_line_set(ref,
		    (const struct rt_brep_internal *)ip->idb_ptr, tol,
		    &obol_published);
	    break;
	case ID_BSPLINE:
	    ok = ged_draw_scene_ref_publish_bspline_wireframe_line_set(ref, ip,
		    tol, &obol_published);
	    break;
	case ID_TOR:
	    ok = ged_draw_scene_ref_publish_tor_wireframe_line_set(ref, ip, ttol);
	    break;
	case ID_RPC:
	    ok = ged_draw_scene_ref_publish_rpc_wireframe_line_set(ref, ip, ttol);
	    break;
	case ID_RHC:
	    ok = ged_draw_scene_ref_publish_rhc_wireframe_line_set(ref, ip, ttol);
	    break;
	case ID_EPA:
	    ok = ged_draw_scene_ref_publish_epa_wireframe_line_set(ref, ip, ttol);
	    break;
	case ID_EHY:
	    ok = ged_draw_scene_ref_publish_ehy_wireframe_line_set(ref, ip, ttol);
	    break;
	case ID_ETO:
	    ok = ged_draw_scene_ref_publish_eto_wireframe_line_set(ref, ip, ttol);
	    break;
	case ID_SUPERELL:
	    ok = ged_draw_scene_ref_publish_superell_wireframe_line_set(ref, ip);
	    break;
	case ID_HRT:
	    ok = ged_draw_scene_ref_publish_hrt_wireframe_line_set(ref, ip, ttol);
	    break;
	case ID_DATUM:
	    ok = ged_draw_scene_ref_publish_datum_wireframe_line_set(ref, ip);
	    break;
	case ID_SKETCH:
	    ok = ged_draw_scene_ref_publish_sketch_wireframe_line_set(ref, ip,
		    ttol);
	    break;
	case ID_EXTRUDE:
	    ok = ged_draw_scene_ref_publish_extrude_wireframe_line_set(ref, ip,
		    ttol);
	    break;
	case ID_REVOLVE:
	    ok = ged_draw_scene_ref_publish_revolve_wireframe_line_set(ref, ip,
		    ttol);
	    break;
	case ID_EBM:
	    ok = ged_draw_scene_ref_publish_ebm_wireframe_line_set(ref, ip);
	    break;
	case ID_VOL:
	    ok = ged_draw_scene_ref_publish_vol_wireframe_line_set(ref, ip);
	    break;
	case ID_DSP:
	    ok = ged_draw_scene_ref_publish_dsp_wireframe_line_set(ref, ip,
		    ttol);
	    break;
	case ID_HF:
	    ok = ged_draw_scene_ref_publish_hf_wireframe_line_set(ref, ip,
		    ttol);
	    break;
	case ID_NMG:
	    ok = ged_draw_scene_ref_publish_nmg_wireframe_line_set(ref, ip);
	    break;
	case ID_SUBMODEL:
	    ok = ged_draw_scene_ref_publish_obol_submodel_wireframe(ref, ip,
		    ttol, tol);
	    if (ok < 0)
		return -1;
	    if (!ok)
		ok = ged_draw_scene_ref_publish_submodel_wireframe_children(ref,
			ip, ttol, tol, view_ctx);
	    break;
	case ID_ANNOT:
	    ok = ged_draw_scene_ref_publish_annot_record(ref, ip);
	    break;
	case ID_ELL:
	case ID_SPH:
	    ok = ged_draw_scene_ref_publish_ell_wireframe_line_set(ref,
		    (const struct rt_ell_internal *)ip->idb_ptr);
	    break;
	case ID_ARB8:
	    ok = ged_draw_scene_ref_publish_arb_wireframe_line_set(ref,
		    (const struct rt_arb_internal *)ip->idb_ptr);
	    break;
	case ID_ARBN:
	    ok = ged_draw_scene_ref_publish_arbn_wireframe_line_set(ref,
		    (const struct rt_arbn_internal *)ip->idb_ptr, tol);
	    break;
	case ID_ARS:
	    ok = ged_draw_scene_ref_publish_ars_wireframe_line_set(ref,
		    (const struct rt_ars_internal *)ip->idb_ptr);
	    break;
	case ID_GRIP:
	    ok = ged_draw_scene_ref_publish_grip_wireframe_line_set(ref,
		    (const struct rt_grip_internal *)ip->idb_ptr);
	    break;
	case ID_HALF:
	    ok = ged_draw_scene_ref_publish_half_wireframe_line_set(ref,
		    (const struct rt_half_internal *)ip->idb_ptr);
	    break;
	case ID_CLINE:
	    ok = ged_draw_scene_ref_publish_cline_wireframe_line_set(ref,
		    (const struct rt_cline_internal *)ip->idb_ptr);
	    break;
	case ID_PARTICLE:
	    ok = ged_draw_scene_ref_publish_part_wireframe_line_set(ref,
		    (const struct rt_part_internal *)ip->idb_ptr);
	    break;
	case ID_JOINT:
	    ok = ged_draw_scene_ref_publish_joint_wireframe_line_set(ref,
		    (const struct rt_joint_internal *)ip->idb_ptr);
	    break;
	case ID_METABALL:
	    ok = ged_draw_scene_ref_publish_metaball_wireframe_line_set(ref,
		    (const struct rt_metaball_internal *)ip->idb_ptr);
	    break;
	case ID_HYP:
	    ok = ged_draw_scene_ref_publish_hyp_wireframe_line_set(ref,
		    (const struct rt_hyp_internal *)ip->idb_ptr);
	    break;
	case ID_PNTS:
	    ok = ged_draw_scene_ref_publish_pnts_wireframe_line_set(ref,
		    (const struct rt_pnts_internal *)ip->idb_ptr);
	    break;
	case ID_PIPE:
	    ok = ged_draw_scene_ref_publish_pipe_wireframe_line_set(ref, ip);
	    break;
	case ID_SCRIPT:
	    ok = ged_draw_scene_ref_publish_script_wireframe_line_set(ref,
		    (const struct rt_script_internal *)ip->idb_ptr);
	    break;
	case ID_TGC:
	case ID_REC:
	    ok = ged_draw_scene_ref_publish_tgc_wireframe_line_set(ref,
		    (const struct rt_tgc_internal *)ip->idb_ptr);
	    break;
	default:
	    return 0;
    }

    if (ok && !obol_published) {
	ged_draw_shape_state *shape_data =
	    _ged_draw_shape_state_get_scene_ref(ref);
	if (shape_data && shape_data->gedp && shape_data->gedp->i &&
		shape_data->gedp->i->ged_gdp)
	    shape_data->gedp->i->ged_gdp->
		gd_draw_retained_primitive_wireframe_publications++;
    }
    return ok ? 1 : -1;
}


static int
ged_draw_scene_ref_publish_primitive_wireframe(bsg_scene_ref ref,
					       struct rt_db_internal *ip,
					       const struct bg_tess_tol *ttol,
					       const struct bn_tol *tol,
					       void *view_ctx,
					       const struct rt_view_info *view_info,
					       int adaptive)
{
    if (!adaptive && ip && (ip->idb_type == ID_BREP ||
	    ip->idb_type == ID_BSPLINE)) {
	int direct = _ged_draw_scene_ref_publish_direct_primitive_wireframe(
		ref, ip, ttol, tol, view_ctx);
	if (direct > 0)
	    return 0;
	if (direct < 0)
	    return -1;
    }

    if (!adaptive || (ip && (ip->idb_type == ID_SUBMODEL ||
		    ip->idb_type == ID_ANNOT))) {
	int obol_direct = ged_draw_scene_ref_publish_obol_primitive_wireframe(
		ref, ip, ttol, tol);
	if (obol_direct > 0)
	    return 0;
	if (obol_direct < 0)
	    return -1;

	int direct = _ged_draw_scene_ref_publish_direct_primitive_wireframe(
		ref, ip, ttol, tol, view_ctx);
	if (direct > 0)
	    return 0;
	if (direct < 0)
	    return -1;
    }

    if (adaptive)
	return _ged_draw_scene_ref_publish_adaptive_wireframe(ref, ip, tol,
		view_info);

    return -1;
}


static int
ged_draw_scene_ref_publish_line_set(bsg_scene_ref ref,
				    const point_t *points,
				    const int *commands,
				    size_t point_count)
{
    ged_draw_shape_state *shape_data = _ged_draw_shape_state_get_scene_ref(ref);
    int obol_published = 0;

    if (!shape_data)
	return 0;
    if (point_count && !points)
	return 0;
    if (shape_data->gedp) {
	ged_draw_shape_ref shape_ref =
	    _ged_draw_shape_ref_from_scene_ref(shape_data->gedp, ref);
	obol_published = ged_draw_shape_ref_obol_publish_line_set(
		shape_data->gedp, shape_ref, ref, points, commands,
		point_count);
	if (!obol_published) {
	    struct _ged_draw_obol_publish_line_ctx ctx = {
		points,
		commands,
		point_count
	    };
	    (void)ged_draw_scene_ref_obol_sync_source_placement(ref);
	    obol_published = ged_draw_scene_ref_obol_source_path_apply(ref,
		    _ged_draw_obol_publish_line_cb, &ctx);
	    if (!obol_published) {
		(void)_ged_draw_scene_ref_obol_ensure_source(ref);
		obol_published =
		    _ged_draw_scene_ref_obol_source_snapshot_path_apply(ref,
			    _ged_draw_obol_publish_line_cb, &ctx);
	    }
	}
    }

    int bsg_published = 0;
    if (!obol_published) {
	bsg_scene_ref geometry_ref =
	    ged_draw_scene_ref_geometry_publication_target(ref);
	if (!bsg_scene_ref_is_null(geometry_ref)) {
	    bsg_geometry_ref geometry = bsg_scene_ref_as_geometry(geometry_ref);
	    bsg_published = bsg_geometry_ref_set_line_set(geometry, points,
		    commands, point_count);
	}
    }
    if (!bsg_published && !obol_published)
	return 0;

    /* Obol owns source-adapter line-set publication; BSG is fallback. */
    shape_data->geometry_command_count = point_count;
    shape_data->geometry_revision++;
    bsg_scene_invalidate(ref);
    return 1;
}


static int
_ged_draw_bot_face_order(const struct rt_bot_internal *bot,
			 size_t face_idx,
			 int vertex_order[3],
			 int normal_order[3])
{
    const int *face;

    if (!bot || !bot->faces || !vertex_order || !normal_order)
	return 0;

    face = &bot->faces[face_idx * 3];
    for (int i = 0; i < 3; i++) {
	if (face[i] < 0 || (size_t)face[i] >= bot->num_vertices)
	    return 0;
    }

    vertex_order[0] = face[0];
    normal_order[0] = 0;
    if (bot->orientation == RT_BOT_CW) {
	vertex_order[1] = face[2];
	vertex_order[2] = face[1];
	normal_order[1] = 2;
	normal_order[2] = 1;
    } else {
	vertex_order[1] = face[1];
	vertex_order[2] = face[2];
	normal_order[1] = 1;
	normal_order[2] = 2;
    }

    return 1;
}


static int
ged_draw_scene_ref_publish_bot_wireframe_line_set(bsg_scene_ref ref,
						  const struct rt_bot_internal *bot)
{
    point_t *points = NULL;
    int *commands = NULL;
    size_t valid_faces = 0;
    size_t idx = 0;
    int ok;

    if (!bot)
	return 0;

    RT_BOT_CK_MAGIC(bot);
    if (!bot->vertices || !bot->faces || !bot->num_vertices || !bot->num_faces)
	return ged_draw_scene_ref_geometry_clear(ref);

    for (size_t i = 0; i < bot->num_faces; i++) {
	int vertex_order[3];
	int normal_order[3];
	if (_ged_draw_bot_face_order(bot, i, vertex_order, normal_order))
	    valid_faces++;
    }

    if (!valid_faces)
	return ged_draw_scene_ref_geometry_clear(ref);
    if (valid_faces > ((size_t)-1) / 4)
	return 0;

    size_t point_count = valid_faces * 4;
    points = (point_t *)bu_calloc(point_count, sizeof(point_t),
	    "GED draw BoT wireframe line-set points");
    commands = (int *)bu_calloc(point_count, sizeof(int),
	    "GED draw BoT wireframe line-set commands");

    for (size_t i = 0; i < bot->num_faces; i++) {
	int vertex_order[3];
	int normal_order[3];

	if (!_ged_draw_bot_face_order(bot, i, vertex_order, normal_order))
	    continue;

	VMOVE(points[idx], &bot->vertices[(size_t)vertex_order[0] * 3]);
	commands[idx++] = BSG_GEOMETRY_LINE_MOVE;
	VMOVE(points[idx], &bot->vertices[(size_t)vertex_order[1] * 3]);
	commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
	VMOVE(points[idx], &bot->vertices[(size_t)vertex_order[2] * 3]);
	commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
	VMOVE(points[idx], &bot->vertices[(size_t)vertex_order[0] * 3]);
	commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
    }

    ok = ged_draw_scene_ref_publish_line_set(ref, (const point_t *)points,
	    commands, idx);

    bu_free(points, "GED draw BoT wireframe line-set points");
    bu_free(commands, "GED draw BoT wireframe line-set commands");
    return ok;
}


static int
_ged_draw_pg_face_valid(const struct rt_pg_face_internal *face)
{
    return face && face->verts && face->npts >= 3;
}


static int
ged_draw_scene_ref_publish_poly_wireframe_line_set(bsg_scene_ref ref,
						   const struct rt_pg_internal *pg)
{
    point_t *points = NULL;
    int *commands = NULL;
    size_t point_count = 0;
    size_t idx = 0;
    int ok;

    if (!pg)
	return 0;

    RT_PG_CK_MAGIC(pg);
    if (!pg->poly || !pg->npoly)
	return ged_draw_scene_ref_geometry_clear(ref);

    for (size_t i = 0; i < pg->npoly; i++) {
	const struct rt_pg_face_internal *face = &pg->poly[i];
	if (!_ged_draw_pg_face_valid(face))
	    continue;
	if (point_count > ((size_t)-1) - (face->npts + 1))
	    return 0;
	point_count += face->npts + 1;
    }

    if (!point_count)
	return ged_draw_scene_ref_geometry_clear(ref);

    points = (point_t *)bu_calloc(point_count, sizeof(point_t),
	    "GED draw POLY wireframe line-set points");
    commands = (int *)bu_calloc(point_count, sizeof(int),
	    "GED draw POLY wireframe line-set commands");

    for (size_t i = 0; i < pg->npoly; i++) {
	const struct rt_pg_face_internal *face = &pg->poly[i];

	if (!_ged_draw_pg_face_valid(face))
	    continue;

	VMOVE(points[idx], &face->verts[0]);
	commands[idx++] = BSG_GEOMETRY_LINE_MOVE;
	for (size_t j = 1; j < face->npts; j++) {
	    VMOVE(points[idx], &face->verts[3 * j]);
	    commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
	}
	VMOVE(points[idx], &face->verts[0]);
	commands[idx++] = BSG_GEOMETRY_LINE_DRAW;
    }

    ok = ged_draw_scene_ref_publish_line_set(ref, (const point_t *)points,
	    commands, idx);

    bu_free(points, "GED draw POLY wireframe line-set points");
    bu_free(commands, "GED draw POLY wireframe line-set commands");
    return ok;
}


static int
_ged_draw_nmg_vertexuse_coord(const struct vertexuse *vu, point_t point)
{
    if (!vu || !point)
	return 0;

    NMG_CK_VERTEXUSE(vu);
    NMG_CK_VERTEX(vu->v_p);
    if (!vu->v_p->vg_p)
	return 0;
    NMG_CK_VERTEX_G(vu->v_p->vg_p);
    VMOVE(point, vu->v_p->vg_p->coord);
    return 1;
}


static int
_ged_draw_nmg_edgeuse_has_cnurb(const struct edgeuse *eu)
{
    if (!eu)
	return 0;
    NMG_CK_EDGEUSE(eu);
    return (eu->g.magic_p && *eu->g.magic_p == NMG_EDGE_G_CNURB_MAGIC) ? 1 : 0;
}


static const struct face_g_snurb *
_ged_draw_nmg_faceuse_snurb(const struct faceuse *fu)
{
    if (!fu || !fu->f_p || !fu->f_p->g.magic_p)
	return NULL;
    if (*fu->f_p->g.magic_p != NMG_FACE_G_SNURB_MAGIC)
	return NULL;
    NMG_CK_FACE_G_SNURB(fu->f_p->g.snurb_p);
    return fu->f_p->g.snurb_p;
}


static fastf_t
_ged_draw_nmg_point_dist_sq(const point_t a, const point_t b)
{
    vect_t diff;

    VSUB2(diff, a, b);
    return MAGSQ(diff);
}


static int
_ged_draw_nmg_cnurb_linear_uv(const struct edgeuse *eu,
			      point_t start_uv,
			      point_t end_uv)
{
    const struct vertexuse *start_vu;
    const struct vertexuse *end_vu;

    if (!eu || !start_uv || !end_uv || !eu->eumate_p)
	return 0;

    start_vu = eu->vu_p;
    end_vu = eu->eumate_p->vu_p;
    if (!start_vu || !end_vu ||
	    !start_vu->a.magic_p || !end_vu->a.magic_p ||
	    *start_vu->a.magic_p != NMG_VERTEXUSE_A_CNURB_MAGIC ||
	    *end_vu->a.magic_p != NMG_VERTEXUSE_A_CNURB_MAGIC)
	return 0;

    VMOVE(start_uv, start_vu->a.cnurb_p->param);
    VMOVE(end_uv, end_vu->a.cnurb_p->param);
    return 1;
}


static int
_ged_draw_nmg_cnurb_eval_point(const struct edge_g_cnurb *cnurb,
			       const struct face_g_snurb *snurb,
			       fastf_t t,
			       point_t point)
{
    hpoint_t uvw = HINIT_ZERO;
    hpoint_t xyz = HINIT_ZERO;
    int coords;

    if (!cnurb || !point || cnurb->order <= 0)
	return 0;

    NMG_CK_EDGE_G_CNURB(cnurb);
    coords = RT_NURB_EXTRACT_COORDS(cnurb->pt_type);
    if (coords < 2 || coords > 4)
	return 0;

    nmg_nurb_c_eval(cnurb, t, uvw);
    if (RT_NURB_IS_PT_RATIONAL(cnurb->pt_type)) {
	int widx = coords - 1;
	fastf_t weight = uvw[widx];

	if (ZERO(weight))
	    return 0;
	for (int i = 0; i < widx; i++)
	    uvw[i] /= weight;
    }

    if (snurb) {
	int scoords;

	NMG_CK_FACE_G_SNURB(snurb);
	nmg_nurb_s_eval(snurb, uvw[0], uvw[1], xyz);
	scoords = RT_NURB_EXTRACT_COORDS(snurb->pt_type);
	if (RT_NURB_IS_PT_RATIONAL(snurb->pt_type)) {
	    int widx = scoords - 1;
	    fastf_t weight = xyz[widx];

	    if (widx < 1 || ZERO(weight))
		return 0;
	    for (int i = 0; i < widx; i++)
		xyz[i] /= weight;
	}
	VSET(point, xyz[X], xyz[Y], xyz[Z]);
    } else {
	VSET(point, uvw[X], uvw[Y], (coords > 2) ? uvw[Z] : 0.0);
    }

    return 1;
}


static int
_ged_draw_nmg_cnurb_sample_count(const struct edgeuse *eu,
				 const struct face_g_snurb *snurb,
				 size_t *point_count)
{
    const struct edge_g_cnurb *cnurb;

    if (point_count)
	*point_count = 0;
    if (!_ged_draw_nmg_edgeuse_has_cnurb(eu))
	return 0;

    cnurb = eu->g.cnurb_p;
    NMG_CK_EDGE_G_CNURB(cnurb);
    if (cnurb->order <= 0 && !snurb) {
	if (point_count)
	    *point_count = 2;
	return 1;
    }
    if (cnurb->order > 0 &&
	    (cnurb->order < 2 || !cnurb->k.knots ||
	     cnurb->k.k_size <= cnurb->order ||
	     cnurb->c_size < 2))
	return 0;

    if (point_count)
	*point_count = GED_DRAW_NMG_CNURB_SAMPLE_POINTS;
    return 1;
}


static int
_ged_draw_nmg_sample_cnurb_edge(const struct edgeuse *eu,
				const struct face_g_snurb *snurb,
				point_t *samples,
				size_t sample_count)
{
    const struct edge_g_cnurb *cnurb;
    point_t start = VINIT_ZERO;
    point_t end = VINIT_ZERO;

    if (!_ged_draw_nmg_edgeuse_has_cnurb(eu) || !samples || sample_count < 2)
	return 0;
    if (!_ged_draw_nmg_vertexuse_coord(eu->vu_p, start) ||
	    !_ged_draw_nmg_vertexuse_coord(eu->eumate_p->vu_p, end))
	return 0;

    cnurb = eu->g.cnurb_p;
    NMG_CK_EDGE_G_CNURB(cnurb);

    VMOVE(samples[0], start);
    VMOVE(samples[sample_count - 1], end);

    if (cnurb->order <= 0) {
	if (snurb) {
	    point_t start_uv = VINIT_ZERO;
	    point_t end_uv = VINIT_ZERO;

	    if (!_ged_draw_nmg_cnurb_linear_uv(eu, start_uv, end_uv))
		return 0;
	    for (size_t i = 1; i + 1 < sample_count; i++) {
		fastf_t t = (fastf_t)i / (fastf_t)(sample_count - 1);
		point_t uvw;

		VBLEND2(uvw, 1.0 - t, start_uv, t, end_uv);
		nmg_eval_linear_trim_curve(snurb, uvw, samples[i]);
	    }
	} else {
	    for (size_t i = 1; i + 1 < sample_count; i++) {
		fastf_t t = (fastf_t)i / (fastf_t)(sample_count - 1);

		VBLEND2(samples[i], 1.0 - t, start, t, end);
	    }
	}
	return 1;
    }

    {
	fastf_t t0 = cnurb->k.knots[cnurb->order - 1];
	fastf_t t1 = cnurb->k.knots[cnurb->k.k_size - cnurb->order + 1];
	point_t cstart = VINIT_ZERO;
	point_t cend = VINIT_ZERO;
	int reverse = 0;

	if (ZERO(t1 - t0))
	    return 0;
	if (!_ged_draw_nmg_cnurb_eval_point(cnurb, snurb, t0, cstart) ||
		!_ged_draw_nmg_cnurb_eval_point(cnurb, snurb, t1, cend))
	    return 0;
	reverse = (_ged_draw_nmg_point_dist_sq(cstart, end) +
		_ged_draw_nmg_point_dist_sq(cend, start) <
		_ged_draw_nmg_point_dist_sq(cstart, start) +
		_ged_draw_nmg_point_dist_sq(cend, end));

	for (size_t i = 1; i + 1 < sample_count; i++) {
	    fastf_t f = (fastf_t)i / (fastf_t)(sample_count - 1);
	    fastf_t t = reverse ? (t1 - (t1 - t0) * f) : (t0 + (t1 - t0) * f);

	    if (!_ged_draw_nmg_cnurb_eval_point(cnurb, snurb, t, samples[i]))
		return 0;
	}
    }

    return 1;
}


static int
_ged_draw_nmg_edge_line_count(const struct edgeuse *eu,
			      const struct face_g_snurb *snurb,
			      size_t *point_count)
{
    point_t point;
    point_t mate_point;

    if (point_count)
	*point_count = 0;
    if (!eu || !eu->eumate_p)
	return 0;

    NMG_CK_EDGEUSE(eu);
    if (_ged_draw_nmg_edgeuse_has_cnurb(eu))
	return _ged_draw_nmg_cnurb_sample_count(eu, snurb, point_count);

    if (!_ged_draw_nmg_vertexuse_coord(eu->vu_p, point) ||
	    !_ged_draw_nmg_vertexuse_coord(eu->eumate_p->vu_p, mate_point))
	return 0;
    if (point_count)
	*point_count = 2;
    return 1;
}


static void
_ged_draw_nmg_snurb_grid_free(struct face_g_snurb *row_refined,
			      struct face_g_snurb *grid,
			      struct knot_vector *tkv1,
			      struct knot_vector *tkv2,
			      struct knot_vector *tau1,
			      struct knot_vector *tau2)
{
    if (grid)
	nmg_nurb_free_snurb(grid);
    if (row_refined)
	nmg_nurb_free_snurb(row_refined);
    if (tau1 && tau1->knots)
	bu_free((char *)tau1->knots, "GED draw NMG snurb tau1 knots");
    if (tau2 && tau2->knots)
	bu_free((char *)tau2->knots, "GED draw NMG snurb tau2 knots");
    if (tkv1 && tkv1->knots)
	bu_free((char *)tkv1->knots, "GED draw NMG snurb tkv1 knots");
    if (tkv2 && tkv2->knots)
	bu_free((char *)tkv2->knots, "GED draw NMG snurb tkv2 knots");
}


static struct face_g_snurb *
_ged_draw_nmg_snurb_grid_refine(const struct face_g_snurb *fg,
				struct face_g_snurb **row_refined,
				struct knot_vector *tkv1,
				struct knot_vector *tkv2,
				struct knot_vector *tau1,
				struct knot_vector *tau2)
{
    struct face_g_snurb *grid = NULL;

    if (row_refined)
	*row_refined = NULL;
    if (!fg || !row_refined || !tkv1 || !tkv2 || !tau1 || !tau2)
	return NULL;

    NMG_CK_FACE_G_SNURB(fg);
    nmg_nurb_kvgen(tkv1, fg->u.knots[0],
	    fg->u.knots[fg->u.k_size - 1],
	    GED_DRAW_NMG_SNURB_INTERIOR_SAMPLES);
    nmg_nurb_kvgen(tkv2, fg->v.knots[0],
	    fg->v.knots[fg->v.k_size - 1],
	    GED_DRAW_NMG_SNURB_INTERIOR_SAMPLES);
    nmg_nurb_kvmerge(tau1, tkv1, &fg->u);
    nmg_nurb_kvmerge(tau2, tkv2, &fg->v);

    *row_refined = nmg_nurb_s_refine(fg, RT_NURB_SPLIT_COL, tau2);
    if (!*row_refined)
	return NULL;
    NMG_CK_SNURB(*row_refined);

    grid = nmg_nurb_s_refine(*row_refined, RT_NURB_SPLIT_ROW, tau1);
    if (!grid)
	return NULL;
    NMG_CK_SNURB(grid);
    return grid;
}


static int
_ged_draw_nmg_snurb_grid_line_count(const struct face_g_snurb *fg,
				    size_t *point_count)
{
    struct knot_vector tkv1 = {0, 0, NULL};
    struct knot_vector tkv2 = {0, 0, NULL};
    struct knot_vector tau1 = {0, 0, NULL};
    struct knot_vector tau2 = {0, 0, NULL};
    struct face_g_snurb *row_refined = NULL;
    struct face_g_snurb *grid = NULL;
    int ok = 0;

    if (point_count)
	*point_count = 0;
    if (!fg)
	return 0;

    grid = _ged_draw_nmg_snurb_grid_refine(fg, &row_refined, &tkv1, &tkv2,
	    &tau1, &tau2);
    if (!grid)
	goto cleanup;
    if (grid->s_size[0] <= 0 || grid->s_size[1] <= 0)
	goto cleanup;

    if (point_count)
	*point_count = (size_t)grid->s_size[0] * (size_t)grid->s_size[1] * 2;
    ok = 1;

cleanup:
    _ged_draw_nmg_snurb_grid_free(row_refined, grid, &tkv1, &tkv2, &tau1,
	    &tau2);
    return ok;
}


static int
_ged_draw_nmg_loop_line_count(const struct loopuse *lu,
			      const struct face_g_snurb *snurb,
			      size_t *point_count)
{
    const struct edgeuse *eu;
    size_t points = 0;
    int first = 1;

    if (point_count)
	*point_count = 0;
    if (!lu)
	return 0;

    NMG_CK_LOOPUSE(lu);
    if (BU_LIST_FIRST_MAGIC(&lu->down_hd) == NMG_VERTEXUSE_MAGIC) {
	point_t point;
	const struct vertexuse *vu = BU_LIST_FIRST(vertexuse, &lu->down_hd);
	if (!_ged_draw_nmg_vertexuse_coord(vu, point))
	    return 0;
	if (point_count)
	    *point_count = 2;
	return 1;
    }
    if (BU_LIST_FIRST_MAGIC(&lu->down_hd) != NMG_EDGEUSE_MAGIC)
	return 0;

    for (BU_LIST_FOR(eu, edgeuse, &lu->down_hd)) {
	size_t edge_points = 0;

	if (!_ged_draw_nmg_edge_line_count(eu, snurb, &edge_points))
	    return 0;
	if (edge_points < 2)
	    return 0;
	points += first ? edge_points : edge_points - 1;
	first = 0;
    }

    if (!points)
	return 0;
    if (point_count)
	*point_count = points;
    return 1;
}


static int
_ged_draw_nmg_loop_normal_line_count(const struct loopuse *lu, size_t *point_count)
{
    const struct edgeuse *eu;
    size_t vertices = 0;
    size_t vertex_normals = 0;

    if (point_count)
	*point_count = 0;
    if (!lu)
	return 0;

    NMG_CK_LOOPUSE(lu);
    if (BU_LIST_FIRST_MAGIC(&lu->down_hd) == NMG_VERTEXUSE_MAGIC)
	return 1;
    if (BU_LIST_FIRST_MAGIC(&lu->down_hd) != NMG_EDGEUSE_MAGIC)
	return 0;

    for (BU_LIST_FOR(eu, edgeuse, &lu->down_hd)) {
	point_t point;
	const struct vertexuse *vu = eu->vu_p;
	if (_ged_draw_nmg_edgeuse_has_cnurb(eu))
	    return 0;
	if (!_ged_draw_nmg_vertexuse_coord(vu, point))
	    return 0;
	if (vu->a.magic_p && *vu->a.magic_p == NMG_VERTEXUSE_A_PLANE_MAGIC)
	    vertex_normals++;
	vertices++;
    }

    if (vertices > 2 && point_count)
	*point_count = 2 + 2 * vertex_normals;
    return 1;
}


static int
_ged_draw_nmg_append_line(point_t *points,
			  int *commands,
			  size_t point_count,
			  size_t *point_idx,
			  const point_t point,
			  int command)
{
    if (!points || !commands || !point_idx || !point || *point_idx >= point_count)
	return 0;

    VMOVE(points[*point_idx], point);
    commands[*point_idx] = command;
    (*point_idx)++;
    return 1;
}


static int
_ged_draw_nmg_append_snurb_grid_lines(point_t *points,
				      int *commands,
				      size_t point_count,
				      size_t *point_idx,
				      const struct face_g_snurb *fg)
{
    struct knot_vector tkv1 = {0, 0, NULL};
    struct knot_vector tkv2 = {0, 0, NULL};
    struct knot_vector tau1 = {0, 0, NULL};
    struct knot_vector tau2 = {0, 0, NULL};
    struct face_g_snurb *row_refined = NULL;
    struct face_g_snurb *grid = NULL;
    int coords = 0;
    int ok = 0;

    if (!points || !commands || !point_idx || !fg)
	return 0;

    grid = _ged_draw_nmg_snurb_grid_refine(fg, &row_refined, &tkv1, &tkv2,
	    &tau1, &tau2);
    if (!grid)
	goto cleanup;
    if (grid->s_size[0] <= 0 || grid->s_size[1] <= 0)
	goto cleanup;

    coords = RT_NURB_EXTRACT_COORDS(grid->pt_type);
    if (coords < 3)
	goto cleanup;

    if (RT_NURB_IS_PT_RATIONAL(grid->pt_type)) {
	fastf_t *vp = grid->ctl_points;
	for (int i = 0; i < grid->s_size[0] * grid->s_size[1]; i++) {
	    fastf_t w = vp[3];
	    if (NEAR_ZERO(w, SMALL_FASTF))
		goto cleanup;
	    fastf_t one_over_w = 1.0 / w;
	    vp[0] *= one_over_w;
	    vp[1] *= one_over_w;
	    vp[2] *= one_over_w;
	    vp[3] *= one_over_w;
	    vp += coords;
	}
    }

    {
	fastf_t *vp = grid->ctl_points;
	for (int row = 0; row < grid->s_size[0]; row++) {
	    point_t point;
	    VSET(point, vp[0], vp[1], vp[2]);
	    if (!_ged_draw_nmg_append_line(points, commands, point_count,
			point_idx, point, BSG_GEOMETRY_LINE_MOVE))
		goto cleanup;
	    vp += coords;
	    for (int col = 1; col < grid->s_size[1]; col++) {
		VSET(point, vp[0], vp[1], vp[2]);
		if (!_ged_draw_nmg_append_line(points, commands, point_count,
			    point_idx, point, BSG_GEOMETRY_LINE_DRAW))
		    goto cleanup;
		vp += coords;
	    }
	}
    }

    for (int col = 0; col < grid->s_size[1]; col++) {
	int stride = grid->s_size[1] * coords;
	fastf_t *vp = &grid->ctl_points[col * coords];
	point_t point;

	VSET(point, vp[0], vp[1], vp[2]);
	if (!_ged_draw_nmg_append_line(points, commands, point_count,
		point_idx, point, BSG_GEOMETRY_LINE_MOVE))
	    goto cleanup;
	vp += stride;
	for (int row = 1; row < grid->s_size[0]; row++) {
	    VSET(point, vp[0], vp[1], vp[2]);
	    if (!_ged_draw_nmg_append_line(points, commands, point_count,
		    point_idx, point, BSG_GEOMETRY_LINE_DRAW))
		goto cleanup;
	    vp += stride;
	}
    }

    ok = 1;

cleanup:
    _ged_draw_nmg_snurb_grid_free(row_refined, grid, &tkv1, &tkv2, &tau1,
	    &tau2);
    return ok;
}


static int
_ged_draw_nmg_append_loop_normal_lines(point_t *points,
				       int *commands,
				       size_t point_count,
				       size_t *point_idx,
				       const struct loopuse *lu,
				       const vect_t normal)
{
    const struct edgeuse *eu;
    point_t centroid = VINIT_ZERO;
    point_t first = VINIT_ZERO;
    size_t vertices = 0;
    int is_first = 1;

    if (!lu)
	return 0;

    NMG_CK_LOOPUSE(lu);
    if (BU_LIST_FIRST_MAGIC(&lu->down_hd) == NMG_VERTEXUSE_MAGIC)
	return 1;
    if (BU_LIST_FIRST_MAGIC(&lu->down_hd) != NMG_EDGEUSE_MAGIC)
	return 0;

    for (BU_LIST_FOR(eu, edgeuse, &lu->down_hd)) {
	point_t point;
	if (_ged_draw_nmg_edgeuse_has_cnurb(eu))
	    return 0;
	if (!_ged_draw_nmg_vertexuse_coord(eu->vu_p, point))
	    return 0;
	if (is_first) {
	    VMOVE(first, point);
	    is_first = 0;
	}
	VADD2(centroid, centroid, point);
	vertices++;
    }

    if (vertices <= 2)
	return 1;

    {
	double scale;
	vect_t tocent;
	point_t tip;

	VSCALE(centroid, centroid, 1.0 / (double)vertices);
	VSUB2(tocent, first, centroid);
	scale = MAGNITUDE(tocent) * 0.5;

	if (!_ged_draw_nmg_append_line(points, commands, point_count, point_idx,
		    centroid, BSG_GEOMETRY_LINE_MOVE))
	    return 0;
	VJOIN1(tip, centroid, scale, normal);
	if (!_ged_draw_nmg_append_line(points, commands, point_count, point_idx,
		    tip, BSG_GEOMETRY_LINE_DRAW))
	    return 0;

	for (BU_LIST_FOR(eu, edgeuse, &lu->down_hd)) {
	    const struct vertexuse *vu = eu->vu_p;
	    point_t point;
	    if (!vu->a.magic_p || *vu->a.magic_p != NMG_VERTEXUSE_A_PLANE_MAGIC)
		continue;
	    if (!_ged_draw_nmg_vertexuse_coord(vu, point))
		return 0;
	    if (!_ged_draw_nmg_append_line(points, commands, point_count,
			point_idx, point, BSG_GEOMETRY_LINE_MOVE))
		return 0;
	    VJOIN1(tip, point, scale, vu->a.plane_p->N);
	    if (!_ged_draw_nmg_append_line(points, commands, point_count,
			point_idx, tip, BSG_GEOMETRY_LINE_DRAW))
		return 0;
	}
    }

    return 1;
}


static int
_ged_draw_nmg_append_loop_lines(point_t *points,
				int *commands,
				size_t point_count,
				size_t *point_idx,
				const struct loopuse *lu,
				const struct face_g_snurb *snurb)
{
    const struct edgeuse *eu;
    int is_first = 1;

    if (!lu)
	return 0;

    NMG_CK_LOOPUSE(lu);
    if (BU_LIST_FIRST_MAGIC(&lu->down_hd) == NMG_VERTEXUSE_MAGIC) {
	point_t point;
	const struct vertexuse *vu = BU_LIST_FIRST(vertexuse, &lu->down_hd);
	if (!_ged_draw_nmg_vertexuse_coord(vu, point))
	    return 0;
	if (!_ged_draw_nmg_append_line(points, commands, point_count, point_idx,
		    point, BSG_GEOMETRY_LINE_MOVE))
	    return 0;
	return _ged_draw_nmg_append_line(points, commands, point_count, point_idx,
		point, BSG_GEOMETRY_LINE_DRAW);
    }
    if (BU_LIST_FIRST_MAGIC(&lu->down_hd) != NMG_EDGEUSE_MAGIC)
	return 0;

    for (BU_LIST_FOR(eu, edgeuse, &lu->down_hd)) {
	if (_ged_draw_nmg_edgeuse_has_cnurb(eu)) {
	    point_t samples[GED_DRAW_NMG_CNURB_SAMPLE_POINTS];
	    size_t sample_count = 0;
	    size_t start = is_first ? 0 : 1;

	    if (!_ged_draw_nmg_cnurb_sample_count(eu, snurb, &sample_count) ||
		    sample_count > GED_DRAW_NMG_CNURB_SAMPLE_POINTS ||
		    !_ged_draw_nmg_sample_cnurb_edge(eu, snurb, samples,
			sample_count))
		return 0;
	    for (size_t i = start; i < sample_count; i++) {
		if (!_ged_draw_nmg_append_line(points, commands, point_count,
			    point_idx, samples[i],
			    (is_first && i == 0) ? BSG_GEOMETRY_LINE_MOVE :
			    BSG_GEOMETRY_LINE_DRAW))
		    return 0;
	    }
	    is_first = 0;
	} else {
	    point_t point;
	    point_t mate_point;

	    if (!_ged_draw_nmg_vertexuse_coord(eu->vu_p, point) ||
		    !_ged_draw_nmg_vertexuse_coord(eu->eumate_p->vu_p,
			mate_point))
		return 0;
	    if (is_first) {
		if (!_ged_draw_nmg_append_line(points, commands, point_count,
			    point_idx, point, BSG_GEOMETRY_LINE_MOVE))
		    return 0;
		is_first = 0;
	    }
	    if (!_ged_draw_nmg_append_line(points, commands, point_count,
			point_idx, mate_point, BSG_GEOMETRY_LINE_DRAW))
		return 0;
	}
    }

    return is_first ? 0 : 1;
}


static int
_ged_draw_nmg_wire_edge_line_count(const struct bu_list *eu_hd, size_t *point_count)
{
    const struct edgeuse *eu;
    size_t points = 0;

    if (point_count)
	*point_count = 0;
    if (!eu_hd)
	return 0;

    for (BU_LIST_FOR(eu, edgeuse, eu_hd)) {
	size_t edge_points = 0;

	NMG_CK_EDGEUSE(eu);
	if (!_ged_draw_nmg_edge_line_count(eu,
		    _ged_draw_nmg_faceuse_snurb(nmg_find_fu_of_eu(eu)),
		    &edge_points))
	    return 0;
	points += edge_points;
    }

    if (point_count)
	*point_count = points;
    return 1;
}


static int
_ged_draw_nmg_append_wire_edge_lines(point_t *points,
				     int *commands,
				     size_t point_count,
				     size_t *point_idx,
				     const struct bu_list *eu_hd)
{
    const struct edgeuse *eu;

    if (!eu_hd)
	return 0;

    for (BU_LIST_FOR(eu, edgeuse, eu_hd)) {
	const struct face_g_snurb *snurb =
	    _ged_draw_nmg_faceuse_snurb(nmg_find_fu_of_eu(eu));

	NMG_CK_EDGEUSE(eu);
	if (_ged_draw_nmg_edgeuse_has_cnurb(eu)) {
	    point_t samples[GED_DRAW_NMG_CNURB_SAMPLE_POINTS];
	    size_t sample_count = 0;

	    if (!_ged_draw_nmg_cnurb_sample_count(eu, snurb, &sample_count) ||
		    sample_count > GED_DRAW_NMG_CNURB_SAMPLE_POINTS ||
		    !_ged_draw_nmg_sample_cnurb_edge(eu, snurb, samples,
			sample_count))
		return 0;
	    for (size_t i = 0; i < sample_count; i++) {
		if (!_ged_draw_nmg_append_line(points, commands, point_count,
			    point_idx, samples[i],
			    (i == 0) ? BSG_GEOMETRY_LINE_MOVE :
			    BSG_GEOMETRY_LINE_DRAW))
		    return 0;
	    }
	} else {
	    point_t point;
	    point_t mate_point;

	    if (!_ged_draw_nmg_vertexuse_coord(eu->vu_p, point) ||
		    !_ged_draw_nmg_vertexuse_coord(eu->eumate_p->vu_p,
			mate_point))
		return 0;
	    if (!_ged_draw_nmg_append_line(points, commands, point_count,
			point_idx, point, BSG_GEOMETRY_LINE_MOVE))
		return 0;
	    if (!_ged_draw_nmg_append_line(points, commands, point_count,
			point_idx, mate_point, BSG_GEOMETRY_LINE_DRAW))
		return 0;
	}
    }

    return 1;
}


static int
_ged_draw_nmg_surface_style(int style)
{
    return ((style & GED_DRAW_NMG_STYLE_POLYGON) &&
	    !(style & GED_DRAW_NMG_STYLE_NO_SURFACES));
}


static void
_ged_draw_scene_ref_copy_display_state(bsg_scene_ref dst, bsg_scene_ref src)
{
    unsigned char rgb[3] = {0, 0, 0};
    unsigned char basecolor[3] = {0, 0, 0};
    mat_t mat;
    point_t bmin, bmax;
    vect_t center = VINIT_ZERO;
    struct ged_draw_shape_material_summary material;

    if (ged_draw_scene_ref_is_null(dst) || ged_draw_scene_ref_is_null(src))
	return;

    (void)ged_draw_scene_ref_set_visible(dst, ged_draw_scene_ref_visible(src));
    if (ged_draw_scene_ref_color(src, rgb))
	(void)ged_draw_scene_ref_set_color(dst, rgb);
    if (ged_draw_scene_ref_material_summary(src, &material)) {
	ged_draw_scene_ref_set_material_rgb(dst, material.material_color);
	(void)ged_draw_scene_ref_set_material_revision(dst,
		material.material_revision);
    }

    MAT_IDN(mat);
    if (ged_draw_scene_ref_transform(src, mat))
	(void)ged_draw_scene_ref_set_transform(dst, mat);
    MAT_IDN(mat);
    if (ged_draw_scene_ref_draw_mat(src, mat))
	(void)ged_draw_scene_ref_set_draw_mat(dst, mat);

    if (ged_draw_scene_ref_bounds(src, bmin, bmax))
	(void)ged_draw_scene_ref_set_bounds(dst, bmin, bmax, 1);
    ged_draw_scene_ref_draw_center(src, center);
    ged_draw_scene_ref_set_draw_center(dst, center);
    (void)ged_draw_scene_ref_set_draw_size(dst,
	    ged_draw_scene_ref_draw_size(src));

    (void)ged_draw_scene_ref_set_line_style(dst,
	    ged_draw_scene_ref_line_style(src));
    (void)ged_draw_scene_ref_set_line_width(dst,
	    ged_draw_scene_ref_line_width(src));
    (void)ged_draw_scene_ref_legacy_basecolor(src, basecolor);
    (void)ged_draw_scene_ref_set_legacy_color_info(dst, basecolor,
	    ged_draw_scene_ref_legacy_user_color(src),
	    ged_draw_scene_ref_legacy_default_color(src));
    (void)ged_draw_scene_ref_set_legacy_eval_flag(dst,
	    ged_draw_scene_ref_legacy_eval_flag(src));
    (void)ged_draw_scene_ref_set_legacy_region_id(dst,
	    ged_draw_scene_ref_legacy_region_id(src));
    (void)ged_draw_scene_ref_set_highlighted(dst,
	    ged_draw_scene_ref_highlighted(src));
    ged_draw_scene_ref_set_draw_mode(dst, ged_draw_scene_ref_draw_mode(src));
    (void)ged_draw_scene_ref_set_transparency(dst,
	    ged_draw_scene_ref_transparency(src));
    (void)ged_draw_scene_ref_set_changed(dst,
	    ged_draw_scene_ref_changed(src));

    (void)ged_draw_scene_ref_copy_draw_intent(dst, src);
}


static void
_ged_draw_scene_ref_copy_display_state_local(bsg_scene_ref dst, bsg_scene_ref src)
{
    unsigned char rgb[3] = {0, 0, 0};
    unsigned char basecolor[3] = {0, 0, 0};
    mat_t mat;
    point_t bmin, bmax;
    vect_t center = VINIT_ZERO;
    if (ged_draw_scene_ref_is_null(dst) || ged_draw_scene_ref_is_null(src))
	return;

    bsg_scene_set_visible(dst, ged_draw_scene_ref_visible(src));
    if (ged_draw_scene_ref_color(src, rgb))
	bsg_scene_set_color(dst, rgb[0], rgb[1], rgb[2]);

    MAT_IDN(mat);
    if (ged_draw_scene_ref_transform(src, mat))
	bsg_scene_set_transform(dst, mat);
    MAT_IDN(mat);
    if (ged_draw_scene_ref_draw_mat(src, mat))
	bsg_scene_set_draw_mat(dst, mat);

    if (ged_draw_scene_ref_bounds(src, bmin, bmax))
	bsg_scene_set_bounds(dst, bmin, bmax, 1);
    ged_draw_scene_ref_draw_center(src, center);
    bsg_scene_set_draw_center(dst, center);
    bsg_scene_set_draw_size(dst, ged_draw_scene_ref_draw_size(src));

    bsg_scene_set_line_style(dst, ged_draw_scene_ref_line_style(src));
    bsg_scene_set_line_width(dst, ged_draw_scene_ref_line_width(src));
    (void)ged_draw_scene_ref_legacy_basecolor(src, basecolor);
    bsg_scene_set_legacy_color_info(dst, basecolor,
	    ged_draw_scene_ref_legacy_user_color(src),
	    ged_draw_scene_ref_legacy_default_color(src));
    bsg_scene_set_legacy_eval_flag(dst,
	    ged_draw_scene_ref_legacy_eval_flag(src));
    bsg_scene_set_legacy_region_id(dst,
	    ged_draw_scene_ref_legacy_region_id(src));
    bsg_scene_set_highlighted(dst, ged_draw_scene_ref_highlighted(src));
    bsg_scene_set_dmode(dst, ged_draw_scene_ref_draw_mode(src));
    bsg_scene_set_transparency(dst, ged_draw_scene_ref_transparency(src));
    bsg_scene_set_changed(dst, ged_draw_scene_ref_changed(src));

    (void)ged_draw_scene_ref_copy_draw_intent(dst, src);
}


static void
_ged_draw_scene_ref_set_semantic_path(bsg_scene_ref ref, const char *path)
{
    if (bsg_scene_ref_is_null(ref) || !path || !path[0])
	return;

    const char *semantic_path = path;
    while (*semantic_path == '/')
	semantic_path++;
    if (!semantic_path[0])
	semantic_path = path;

    if (ged_draw_scene_ref_set_draw_intent_path(ref, semantic_path,
	    GED_DRAW_MODE_WIRE))
	(void)ged_draw_scene_ref_set_draw_intent_mode(ref,
		GED_DRAW_MODE_WIRE, semantic_path);
}


static bsg_scene_ref
_ged_draw_scene_ref_create_child_shape(bsg_scene_ref primary_ref,
				       const char *source_name,
				       const char *shape_name,
				       bsg_scene_ref *child_source_out)
{
    if (child_source_out)
	*child_source_out = bsg_scene_ref_null();

    if (bsg_scene_ref_is_null(primary_ref))
	return bsg_scene_ref_null();

    ged_draw_shape_state *shape_data =
	ged_draw_scene_ref_shape_state(primary_ref);
    bsg_scene_ref parent_source_ref = ged_draw_scene_ref_source_owner(primary_ref);
    void *view_ctx = ged_draw_scene_ref_publication_view_context(primary_ref);
    if (!view_ctx)
	view_ctx = ged_draw_scene_ref_publication_view_context(parent_source_ref);

    if (!shape_data || !shape_data->gedp || !view_ctx ||
	    bsg_scene_ref_is_null(parent_source_ref))
	return bsg_scene_ref_null();

    bsg_scene_ref child_source_ref =
	ged_draw_view_context_database_source_create(view_ctx,
		source_name ? source_name : "_submodel_source");
    if (bsg_scene_ref_is_null(child_source_ref))
	return bsg_scene_ref_null();

    bsg_scene_ref child_shape_ref =
	ged_draw_view_context_geometry_create(view_ctx,
		shape_name ? shape_name : "submodel_leaf");
    if (bsg_scene_ref_is_null(child_shape_ref)) {
	ged_draw_scene_ref_release(child_source_ref);
	return bsg_scene_ref_null();
    }

    if (!ged_draw_scene_ref_append_child(child_source_ref, child_shape_ref)) {
	ged_draw_scene_ref_release(child_shape_ref);
	ged_draw_scene_ref_release(child_source_ref);
	return bsg_scene_ref_null();
    }
    if (ged_draw_shape_ref_is_null(
	    _ged_draw_shape_ref_from_scene_ref(shape_data->gedp,
		child_shape_ref))) {
	ged_draw_scene_ref_release(child_source_ref);
	return bsg_scene_ref_null();
    }

    _ged_draw_scene_ref_set_source_ref(child_shape_ref, child_source_ref);
    ged_draw_scene_ref_geometry_clear(child_shape_ref);

    _ged_draw_scene_ref_copy_display_state_local(child_source_ref,
	    primary_ref);
    _ged_draw_scene_ref_copy_display_state_local(child_shape_ref,
	    primary_ref);
    _ged_draw_scene_ref_apply_db_object_marker(child_shape_ref);
    ged_draw_scene_ref_set_non_database_source(child_source_ref, 0);
    ged_draw_scene_ref_set_non_database_source(child_shape_ref, 0);

    if (!ged_draw_scene_ref_append_child(parent_source_ref, child_source_ref)) {
	ged_draw_scene_ref_release(child_source_ref);
	return bsg_scene_ref_null();
    }
    if (shape_data->gedp && shape_data->gedp->i &&
	    shape_data->gedp->i->ged_gdp)
	shape_data->gedp->i->ged_gdp->gd_draw_retained_child_source_creations++;
    if (child_source_out)
	*child_source_out = child_source_ref;
    return child_shape_ref;
}


static void
_ged_draw_scene_ref_clear_child_sources(bsg_scene_ref primary_ref)
{
    bsg_scene_ref source_ref = ged_draw_scene_ref_source_owner(primary_ref);
    if (ged_draw_scene_ref_is_null(primary_ref) ||
	    ged_draw_scene_ref_is_null(source_ref))
	return;

    size_t child_count = ged_draw_scene_ref_child_count(source_ref);
    while (child_count > 0) {
	child_count--;
	bsg_scene_ref child_ref =
	    ged_draw_scene_ref_child_at(source_ref, child_count);
	if (ged_draw_scene_ref_equal(child_ref, primary_ref))
	    continue;
	(void)ged_draw_scene_ref_detach(child_ref);
	ged_draw_scene_ref_release(child_ref);
    }
    ged_draw_scene_ref_invalidate(source_ref);
}


static int
_ged_draw_scene_ref_obol_source_snapshot_path_apply(
	bsg_scene_ref ref,
	ged_draw_scene_ref_obol_source_path_cb cb,
	void *data)
{
    if (ged_draw_scene_ref_obol_source_path_apply(ref, cb, data))
	return 1;

    struct ged *gedp = _ged_draw_scene_ref_owner_gedp(ref);
    if (!gedp || !cb)
	return 0;

    struct ged_draw_shape_source_snapshot source;
    memset(&source, 0, sizeof(source));
    if (!ged_draw_scene_ref_source_snapshot(ref, &source))
	return 0;

    if (source.fullpath) {
	char *path = db_path_to_string(source.fullpath);
	if (path) {
	    int matched = (*cb)(gedp, path, data);
	    bu_free(path, "GED Obol source snapshot fullpath");
	    if (matched)
		return 1;
	}
    }

    return (source.name && source.name[0]) ?
	(*cb)(gedp, source.name, data) : 0;
}


static char *
_ged_draw_scene_ref_owner_path_string(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return NULL;

    const struct db_full_path *fp = ged_draw_scene_ref_fullpath(ref);
    if (fp && fp->fp_len > 0)
	return db_path_to_string(fp);

    const char *intent_path = ged_draw_scene_ref_draw_intent_path(ref);
    if (intent_path && intent_path[0])
	return bu_strdup(intent_path);

    rt_view_scene_ref scene_ref = ged_draw_scene_ref_to_rt_view_ref(ref);
    const char *semantic_path = ged_draw_scene_rt_ref_semantic_path(scene_ref);
    if (semantic_path && semantic_path[0])
	return bu_strdup(semantic_path);

    struct ged_draw_shape_source_snapshot source;
    memset(&source, 0, sizeof(source));
    if (!ged_draw_scene_ref_source_snapshot(ref, &source))
	return NULL;

    if (source.fullpath && source.fullpath->fp_len > 0)
	return db_path_to_string(source.fullpath);

    return (source.name && source.name[0]) ? bu_strdup(source.name) : NULL;
}


static int
_ged_draw_scene_ref_obol_ensure_source(bsg_scene_ref primary_ref)
{
    struct ged_draw_shape_source_snapshot source;
    struct ged_draw_scene_display_summary display;

    if (ged_draw_scene_ref_is_null(primary_ref))
	return 0;

    struct ged *gedp = _ged_draw_scene_ref_owner_gedp(primary_ref);
    if (!gedp || !ged_draw_obol_scene_controller_ensure_owned(gedp, 1))
	return 0;

    memset(&source, 0, sizeof(source));
    if (!ged_draw_scene_ref_source_snapshot(primary_ref, &source) ||
	    !source.dbip)
	return 0;

    memset(&display, 0, sizeof(display));
    int draw_mode = GED_DRAW_MODE_WIRE;
    if (ged_draw_scene_ref_display_summary(primary_ref, &display) &&
	    display.valid && display.draw_mode >= 0)
	draw_mode = display.draw_mode;

    uint64_t source_revision = 0;
    struct ged_draw_shape_record_summary record_summary;
    memset(&record_summary, 0, sizeof(record_summary));
    if (ged_draw_scene_ref_shape_record_summary(primary_ref,
	    &record_summary))
	source_revision = record_summary.source_revision;

    char *path = _ged_draw_scene_ref_owner_path_string(primary_ref);
    if (path) {
	int matched = ged_draw_obol_database_source_ensure_for_path(
		gedp, path, source.dbip, draw_mode, source_revision);
	bu_free(path, "GED Obol source ensure path");
	if (matched)
	    return 1;
    }

    return (source.name && source.name[0]) ?
	ged_draw_obol_database_source_ensure_for_path(
		gedp, source.name, source.dbip, draw_mode,
		source_revision) : 0;
}


static int
_ged_draw_scene_ref_obol_ensure_owner_source(bsg_scene_ref primary_ref)
{
    if (_ged_draw_scene_ref_obol_ensure_source(primary_ref))
	return 1;

    struct ged *gedp = _ged_draw_scene_ref_owner_gedp(primary_ref);
    if (!gedp || !gedp->dbip)
	return 0;

    char *path = _ged_draw_scene_ref_owner_path_string(primary_ref);
    if (!path)
	return 0;

    int draw_mode = ged_draw_scene_ref_draw_mode(primary_ref);
    if (draw_mode < 0)
	draw_mode = GED_DRAW_MODE_WIRE;
    int matched = ged_draw_obol_database_source_ensure_for_path(gedp,
	    path, gedp->dbip, draw_mode, 0);
    bu_free(path, "GED Obol owner source ensure path");
    return matched;
}


static int
_ged_draw_scene_ref_publish_obol_auxiliary_from_line_ref(
	bsg_scene_ref primary_ref,
	bsg_scene_ref line_ref,
	const char *name)
{
    if (ged_draw_scene_ref_is_null(primary_ref) ||
	    ged_draw_scene_ref_is_null(line_ref) || !name || !name[0])
	return 0;

    struct ged_draw_view_line_summary summary;
    memset(&summary, 0, sizeof(summary));
    struct _ged_draw_obol_line_summary_ctx obol_summary_ctx = {&summary};
    int line_from_obol =
	_ged_draw_scene_ref_obol_source_snapshot_path_apply(line_ref,
		_ged_draw_obol_line_summary_cb, &obol_summary_ctx);
    if (!line_from_obol &&
	    !ged_draw_scene_ref_line_summary(line_ref, &summary))
	return 0;
    if (!summary.valid || summary.point_count == 0)
	return 0;

    point_t *points = (point_t *)bu_calloc(summary.point_count,
	    sizeof(point_t), "GED Obol auxiliary submodel points");
    int *commands = (int *)bu_calloc(summary.point_count, sizeof(int),
	    "GED Obol auxiliary submodel commands");
    int ok = 0;

    for (size_t i = 0; i < summary.point_count; i++) {
	int point_ok = 0;
	int command_ok = 0;
	if (line_from_obol) {
	    struct _ged_draw_obol_line_point_ctx point_ctx;
	    point_ctx.index = i;
	    point_ctx.out = points[i];
	    point_ok =
		_ged_draw_scene_ref_obol_source_snapshot_path_apply(line_ref,
			_ged_draw_obol_line_point_cb, &point_ctx);

	    struct _ged_draw_obol_line_command_ctx command_ctx = {i,
		&commands[i]};
	    command_ok =
		_ged_draw_scene_ref_obol_source_snapshot_path_apply(line_ref,
			_ged_draw_obol_line_command_cb, &command_ctx);
	} else {
	    point_ok = ged_draw_scene_ref_line_point_at(line_ref, i,
		    points[i]);
	    command_ok = ged_draw_scene_ref_line_command_at(line_ref, i,
		    &commands[i]);
	}
	if (!point_ok || !command_ok)
	    goto cleanup;
    }

    (void)_ged_draw_scene_ref_obol_ensure_source(primary_ref);
    ok = ged_draw_scene_ref_obol_publish_auxiliary_line_set(primary_ref,
	    name, (const point_t *)points, commands, summary.point_count);

cleanup:
    bu_free(points, "GED Obol auxiliary submodel points");
    bu_free(commands, "GED Obol auxiliary submodel commands");
    return ok;
}


struct ged_draw_submodel_publish_ctx {
    bsg_scene_ref parent_ref;
    const struct bg_tess_tol *ttol;
    const struct bn_tol *tol;
    void *view_ctx;
    size_t child_count;
    int failed;
};


static union tree *
_ged_draw_submodel_wireframe_leaf(struct db_tree_state *tsp,
				  const struct db_full_path *pathp,
				  struct rt_db_internal *ip,
				  void *client_data)
{
    struct ged_draw_submodel_publish_ctx *ctx =
	(struct ged_draw_submodel_publish_ctx *)client_data;
    char *path_name = NULL;
    const char *name = "submodel_leaf";
    bsg_scene_ref child_source_ref = bsg_scene_ref_null();
    int child_obol_source_preexisting = 0;

    if (!ctx || ctx->failed || !tsp || !ip)
	return TREE_NULL;

    if (pathp && pathp->fp_len > 0) {
	path_name = db_path_to_string(pathp);
	if (path_name && path_name[0])
	    name = path_name;
    } else if (ip->idb_meth) {
	name = ip->idb_meth->ft_name;
    }

    struct ged *gedp = _ged_draw_scene_ref_owner_gedp(ctx->parent_ref);
    if (gedp && name && name[0]) {
	struct ged_draw_database_source_summary source_summary;
	memset(&source_summary, 0, sizeof(source_summary));
	child_obol_source_preexisting =
	    ged_draw_obol_database_source_summary_for_path(gedp, name,
		    &source_summary) &&
	    source_summary.valid;
    }

    bsg_scene_ref child_shape_ref =
	_ged_draw_scene_ref_create_child_shape(ctx->parent_ref, name, name,
		&child_source_ref);
    if (bsg_scene_ref_is_null(child_shape_ref)) {
	ctx->failed = 1;
	if (path_name)
	    bu_free(path_name, "submodel leaf path string");
	return TREE_NULL;
    }

    _ged_draw_scene_ref_set_semantic_path(child_source_ref, name);
    _ged_draw_scene_ref_set_semantic_path(child_shape_ref, name);

    int ret = ged_draw_scene_ref_publish_primitive_wireframe(child_shape_ref,
	    ip, ctx->ttol ? ctx->ttol : tsp->ts_ttol,
	    ctx->tol ? ctx->tol : tsp->ts_tol, ctx->view_ctx, NULL, 0);
    if (ret < 0) {
	bu_log("ged draw submodel leaf %s (%s) wireframe failure\n",
		name, ip->idb_meth ? ip->idb_meth->ft_name : "unknown");
	if (!ged_draw_scene_ref_is_null(child_source_ref)) {
	    (void)ged_draw_scene_ref_detach(child_source_ref);
	    ged_draw_scene_ref_release(child_source_ref);
	}
	ctx->failed = 1;
	if (path_name)
	    bu_free(path_name, "submodel leaf path string");
	return TREE_NULL;
    }

    int obol_aux_published =
	_ged_draw_scene_ref_publish_obol_auxiliary_from_line_ref(
		ctx->parent_ref, child_shape_ref, name);
    if (obol_aux_published) {
	/* Obol auxiliary VLISTs own published submodel leaves. */
	(void)ged_draw_scene_ref_detach(child_source_ref);
	ged_draw_scene_ref_release(child_source_ref);
	child_source_ref = bsg_scene_ref_null();
	if (gedp && name && name[0] && !child_obol_source_preexisting)
	    (void)ged_draw_obol_database_source_remove_for_path(gedp, name);
    } else {
	(void)ged_draw_scene_ref_update_bounds_from_geometry(child_shape_ref,
		NULL);
	if (!ged_draw_scene_ref_is_null(child_source_ref))
	    (void)ged_draw_scene_ref_update_bounds_context(child_source_ref,
		ctx->view_ctx);
    }
    ctx->child_count++;

    if (path_name)
	bu_free(path_name, "submodel leaf path string");

    union tree *curtree;
    BU_GET(curtree, union tree);
    RT_TREE_INIT(curtree);
    curtree->tr_op = OP_NOP;
    return curtree;
}


static int
ged_draw_scene_ref_publish_submodel_wireframe_children(bsg_scene_ref ref,
						       struct rt_db_internal *ip,
						       const struct bg_tess_tol *ttol,
						       const struct bn_tol *tol,
						       void *view_ctx)
{
    if (bsg_scene_ref_is_null(ref) || !ip || ip->idb_type != ID_SUBMODEL ||
	    !ip->idb_ptr)
	return 0;

    struct rt_submodel_internal *sip =
	(struct rt_submodel_internal *)ip->idb_ptr;
    RT_SUBMODEL_CK_MAGIC(sip);

    _ged_draw_scene_ref_clear_child_sources(ref);
    ged_draw_scene_ref_geometry_clear(ref);

    struct db_i *dbip = DBI_NULL;
    int close_db = 0;
    if (bu_vls_strlen(&sip->file) != 0) {
	dbip = db_open(bu_vls_addr(&sip->file), DB_OPEN_READONLY);
	if (dbip == DBI_NULL) {
	    bu_log("Cannot open geometry database file (%s) to store plot\n",
		    bu_vls_addr(&sip->file));
	    return 0;
	}
	close_db = 1;
	if (!db_is_directory_non_empty(dbip) && db_dirbuild(dbip) < 0) {
	    bu_log("ged draw submodel db_dirbuild() failure\n");
	    db_close(dbip);
	    return 0;
	}
    } else {
	RT_CK_DBI(sip->dbip);
	dbip = (struct db_i *)sip->dbip;
    }

    struct db_tree_state state;
    RT_DBTS_INIT(&state);
    state.ts_dbip = dbip;
    state.ts_ttol = ttol;
    state.ts_tol = tol;
    MAT_COPY(state.ts_mat, sip->root2leaf);

    struct ged_draw_submodel_publish_ctx ctx;
    ctx.parent_ref = ref;
    ctx.ttol = ttol;
    ctx.tol = tol;
    ctx.view_ctx = view_ctx ? view_ctx : ged_draw_scene_ref_view_context(ref);
    ctx.child_count = 0;
    ctx.failed = 0;

    const char *argv[2];
    argv[0] = bu_vls_addr(&sip->treetop);
    argv[1] = NULL;
    int ret = db_walk_tree(dbip, 1, argv, 1, &state, 0, NULL,
	    _ged_draw_submodel_wireframe_leaf, (void *)&ctx);

    if (close_db)
	db_close(dbip);

    if (ret < 0 || ctx.failed) {
	bu_log("ged draw submodel db_walk_tree(%s) failure\n",
		bu_vls_addr(&sip->treetop));
	return 0;
    }

    bsg_scene_ref source_ref = ged_draw_scene_ref_source_owner(ref);
    if (!ged_draw_scene_ref_is_null(source_ref))
	(void)ged_draw_scene_ref_update_bounds_context(source_ref, ctx.view_ctx);
    (void)ged_draw_scene_ref_update_bounds_context(ref, ctx.view_ctx);
    return 1;
}


static int
_ged_draw_nmg_line_set_stats(const struct nmgregion *r,
			     int style,
			     size_t *point_count)
{
    const struct shell *s;
    const int surface_style = _ged_draw_nmg_surface_style(style);
    size_t points = 0;

    if (point_count)
	*point_count = 0;
    if (!r)
	return 0;

    if (surface_style && (style & GED_DRAW_NMG_STYLE_VISUALIZE_NORMALS))
	return 0;

    NMG_CK_REGION(r);
    for (BU_LIST_FOR(s, shell, &r->s_hd)) {
	const struct faceuse *fu;
	const struct loopuse *lu;
	size_t shell_points = 0;

	NMG_CK_SHELL(s);
	for (BU_LIST_FOR(fu, faceuse, &s->fu_hd)) {
	    const struct face_g_snurb *snurb = NULL;

	    if (fu->orientation != OT_SAME)
		continue;
	    if (surface_style)
		return 0;

	    NMG_CK_FACEUSE(fu);
	    NMG_CK_FACE(fu->f_p);
	    if (!fu->f_p->g.magic_p)
		return 0;
	    if (*fu->f_p->g.magic_p == NMG_FACE_G_SNURB_MAGIC) {
		if (style & GED_DRAW_NMG_STYLE_VISUALIZE_NORMALS)
		    return 0;
		snurb = _ged_draw_nmg_faceuse_snurb(fu);
		if (!(style & GED_DRAW_NMG_STYLE_NO_SURFACES)) {
		    size_t snurb_points = 0;
		    if (!_ged_draw_nmg_snurb_grid_line_count(snurb,
				&snurb_points))
			return 0;
		    shell_points += snurb_points;
		}
	    } else if (*fu->f_p->g.magic_p == NMG_FACE_G_PLANE_MAGIC) {
		NMG_CK_FACE_G_PLANE(fu->f_p->g.plane_p);
	    } else {
		return 0;
	    }

	    for (BU_LIST_FOR(lu, loopuse, &fu->lu_hd)) {
		size_t loop_points = 0;
		if (!_ged_draw_nmg_loop_line_count(lu, snurb, &loop_points))
		    return 0;
		shell_points += loop_points;
		if (style & GED_DRAW_NMG_STYLE_VISUALIZE_NORMALS) {
		    size_t normal_points = 0;
		    if (!_ged_draw_nmg_loop_normal_line_count(lu, &normal_points))
			return 0;
		    shell_points += normal_points;
		}
	    }
	}

	for (BU_LIST_FOR(lu, loopuse, &s->lu_hd)) {
	    size_t loop_points = 0;
	    if (!_ged_draw_nmg_loop_line_count(lu, NULL, &loop_points))
		return 0;
	    shell_points += loop_points;
	}

	{
	    size_t edge_points = 0;
	    if (!_ged_draw_nmg_wire_edge_line_count(&s->eu_hd, &edge_points))
		return 0;
	    shell_points += edge_points;
	}

	if (s->vu_p) {
	    point_t point;
	    if (!_ged_draw_nmg_vertexuse_coord(s->vu_p, point))
		return 0;
	    shell_points += 2;
	}

	points += shell_points;
    }

    if (!points)
	return 0;
    if (point_count)
	*point_count = points;
    return 1;
}


static int
_ged_draw_nmg_append_region_line_set(const struct nmgregion *r,
				     int style,
				     point_t *points,
				     int *commands,
				     size_t point_count,
				     size_t *point_idx)
{
    const struct shell *s;
    const int surface_style = _ged_draw_nmg_surface_style(style);

    if (!r || !points || !commands || !point_idx)
	return 0;

    for (BU_LIST_FOR(s, shell, &r->s_hd)) {
	const struct faceuse *fu;
	const struct loopuse *lu;

	for (BU_LIST_FOR(fu, faceuse, &s->fu_hd)) {
	    const struct face_g_snurb *snurb = NULL;
	    vect_t face_normal = VINIT_ZERO;

	    if (fu->orientation != OT_SAME)
		continue;
	    if (surface_style)
		return 0;
	    snurb = _ged_draw_nmg_faceuse_snurb(fu);
	    if (!snurb) {
		NMG_GET_FU_NORMAL(face_normal, fu);
	    } else if (style & GED_DRAW_NMG_STYLE_VISUALIZE_NORMALS) {
		return 0;
	    } else if (!(style & GED_DRAW_NMG_STYLE_NO_SURFACES)) {
		if (!_ged_draw_nmg_append_snurb_grid_lines(points, commands,
			    point_count, point_idx, snurb))
		    return 0;
	    }
	    for (BU_LIST_FOR(lu, loopuse, &fu->lu_hd)) {
		if (!_ged_draw_nmg_append_loop_lines(points, commands,
			    point_count, point_idx, lu, snurb))
		    return 0;
		if ((style & GED_DRAW_NMG_STYLE_VISUALIZE_NORMALS) &&
			!_ged_draw_nmg_append_loop_normal_lines(points, commands,
			    point_count, point_idx, lu, face_normal))
		    return 0;
	    }
	}

	for (BU_LIST_FOR(lu, loopuse, &s->lu_hd)) {
	    if (!_ged_draw_nmg_append_loop_lines(points, commands, point_count,
			point_idx, lu, NULL))
		return 0;
	}

	if (!_ged_draw_nmg_append_wire_edge_lines(points, commands, point_count,
		    point_idx, &s->eu_hd))
	    return 0;

	if (s->vu_p) {
	    point_t point;
	    if (!_ged_draw_nmg_vertexuse_coord(s->vu_p, point))
		return 0;
	    if (!_ged_draw_nmg_append_line(points, commands, point_count,
			point_idx, point, BSG_GEOMETRY_LINE_MOVE))
		return 0;
	    if (!_ged_draw_nmg_append_line(points, commands, point_count,
			point_idx, point, BSG_GEOMETRY_LINE_DRAW))
		return 0;
	}
    }

    return 1;
}


static int
_ged_draw_nmg_model_line_set_stats(const struct model *m,
				   int style,
				   size_t *point_count)
{
    const struct nmgregion *r;
    size_t points = 0;
    size_t regions = 0;

    if (point_count)
	*point_count = 0;
    if (!m)
	return 0;

    NMG_CK_MODEL(m);
    for (BU_LIST_FOR(r, nmgregion, &m->r_hd)) {
	size_t region_points = 0;
	NMG_CK_REGION(r);
	if (!_ged_draw_nmg_line_set_stats(r, style, &region_points))
	    return 0;
	points += region_points;
	regions++;
    }

    if (!regions || !points)
	return 0;
    if (point_count)
	*point_count = points;
    return 1;
}


static int
_ged_draw_nmg_fill_model_line_set(const struct model *m,
				  int style,
				  point_t *points,
				  int *commands,
				  size_t point_count)
{
    const struct nmgregion *r;
    size_t point_idx = 0;

    if (!m || !points || !commands)
	return 0;

    NMG_CK_MODEL(m);
    for (BU_LIST_FOR(r, nmgregion, &m->r_hd)) {
	NMG_CK_REGION(r);
	if (!_ged_draw_nmg_append_region_line_set(r, style, points, commands,
		point_count, &point_idx))
	    return 0;
    }

    return (point_idx == point_count) ? 1 : 0;
}


static int
ged_draw_scene_ref_geometry_publish_nmg_model(bsg_scene_ref ref,
					     const struct model *m,
					     int style)
{
    point_t *points = NULL;
    int *commands = NULL;
    size_t point_count = 0;
    int ok = 0;

    if (bsg_scene_ref_is_null(ref) || !m)
	return 0;

    if (!_ged_draw_nmg_model_line_set_stats(m, style, &point_count))
	return 0;

    points = (point_t *)bu_calloc(point_count, sizeof(point_t),
	    "GED draw NMG model line-set points");
    commands = (int *)bu_calloc(point_count, sizeof(int),
	    "GED draw NMG model line-set commands");

    if (!_ged_draw_nmg_fill_model_line_set(m, style, points, commands,
		point_count))
	goto cleanup;

    ok = ged_draw_scene_ref_publish_line_set(ref, (const point_t *)points,
	    commands, point_count);

cleanup:
    if (points)
	bu_free(points, "GED draw NMG model line-set points");
    if (commands)
	bu_free(commands, "GED draw NMG model line-set commands");
    return ok;
}


static void
_ged_draw_scene_ref_geometry_set_command_count(bsg_scene_ref ref, size_t vlen)
{
    ged_draw_shape_state *shape_data =
	_ged_draw_shape_state_get_scene_ref(ref);
    if (!shape_data)
	return;
    shape_data->geometry_command_count = vlen;
    shape_data->geometry_revision++;
}



/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
