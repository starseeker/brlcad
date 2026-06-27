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
#include <string.h>

#include "bg/plane.h"
#include "bg/sat.h"
#include "bu/log.h"
#include "bu/malloc.h"
#include "bu/parallel.h"
#include "bu/str.h"
#include "bsg/appearance.h"
#include "bsg/database_source.h"
#include "bsg/draw_ctx.h"
#include "bsg/draw_set.h"
#include "bsg/draw_source.h"
#include "bsg/draw_intent.h"
#include "bsg/export.h"
#include "bsg/field.h"
#include "bsg/feature.h"
#include "bsg/geometry.h"
#include "bsg/material.h"
#include "bsg/overlay.h"
#include "bsg/render.h"
#include "bsg/scene_builder.h"
#include "ged/draw.h"
#include "ged/db_index.h"
#include "ged/view.h"
#include "nmg.h"
#include "rt/db_instance.h"
#include "rt/db_io.h"
#include "rt/functab.h"
#include "rt/geom.h"
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


static void _ged_draw_scene_ref_geometry_set_command_count(bsg_scene_ref ref, size_t vlen);
static bsg_scene_ref ged_draw_view_context_overlay_geometry_create(
	void *view_ctx,
	const char *name,
	enum ged_draw_overlay_geometry_kind kind);
static bsg_scene_ref ged_draw_view_context_overlay_scene_find(void *view_ctx,
							      const char *name);
static void ged_draw_view_context_overlay_name_erase(void *view_ctx,
						     const char *name);
static int ged_draw_view_context_overlay_scene_append(void *view_ctx,
						      bsg_scene_ref scene);
static int ged_draw_view_overlay_command_result_owner_set(bsg_scene_ref scene,
							  const void *owner,
							  const char *source_path);
static bsg_scene_ref ged_draw_view_context_overlay_create(void *view_ctx,
							  const char *name);
static bsg_scene_ref ged_draw_view_context_database_source_create(void *view_ctx,
								  const char *name);
static bsg_scene_ref ged_draw_view_context_geometry_create(void *view_ctx,
							   const char *name);
static int ged_draw_scene_ref_detach(bsg_scene_ref ref);
static int ged_draw_scene_ref_append_child(bsg_scene_ref parent_ref,
					   bsg_scene_ref child_ref);
static int ged_draw_scene_ref_append_source_owner_to_group(
	bsg_scene_ref group_ref,
	bsg_scene_ref ref);
static bsg_scene_ref ged_draw_view_context_group_child_ensure(
	void *view_ctx,
	bsg_scene_ref parent_ref,
	const char *name,
	void *dp_hint);
static int ged_draw_scene_ref_is_database_source(bsg_scene_ref ref);
static void ged_draw_scene_ref_release(bsg_scene_ref ref);
static int ged_draw_scene_ref_create_draft_pair(struct ged *gedp,
						void *view_ctx,
						bsg_scene_ref *source_out,
						bsg_scene_ref *shape_out);
static int ged_draw_group_ref_append_scene_ref(struct ged *gedp,
					       ged_draw_group_ref ref,
					       bsg_scene_ref shape_ref);
static size_t ged_draw_scene_ref_child_count(bsg_scene_ref ref);
static bsg_scene_ref ged_draw_scene_ref_child_at(bsg_scene_ref ref, size_t idx);
static bsg_scene_ref ged_draw_scene_ref_shape_owning_group(bsg_scene_ref ref);
static int ged_draw_scene_ref_draw_tree_depth(bsg_scene_ref ref);
static int ged_draw_scene_ref_draw_mode(bsg_scene_ref ref);
static void ged_draw_scene_ref_realization_set_view_context_policy(
	bsg_scene_ref ref,
	const void *view_ctx);
static void ged_draw_scene_ref_set_draw_center(bsg_scene_ref ref,
					       const point_t center);
static int ged_draw_scene_ref_set_bounds(bsg_scene_ref ref,
					 const point_t min,
					 const point_t max,
					 int valid);
static int ged_draw_scene_ref_sync_aux_display_state(bsg_scene_ref source_ref,
						     bsg_scene_ref shape_ref);
static int ged_draw_scene_ref_set_name(bsg_scene_ref ref, const char *name);
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
static int ged_draw_scene_ref_set_draw_intent_appearance_settings(
	bsg_scene_ref ref,
	const struct ged_draw_appearance_settings *settings,
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
static int ged_draw_scene_ref_set_legacy_uses_default_color(
	bsg_scene_ref ref,
	int default_color);
static int ged_draw_scene_ref_set_legacy_region_id(bsg_scene_ref ref,
						   int region_id);
static int ged_draw_shape_draft_apply_db_object_marker(bsg_scene_ref ref);
static int ged_draw_scene_ref_bump_appearance_revision(bsg_scene_ref ref);
static int ged_draw_scene_ref_apply_overlay_geometry_attributes(
	bsg_scene_ref ref,
	const unsigned char rgb[3],
	fastf_t transparency,
	int draw_mode);
static int ged_draw_scene_ref_apply_display_settings(
	bsg_scene_ref ref,
	const struct ged_draw_appearance_settings *settings);
static void ged_draw_scene_ref_set_material_rgb(bsg_scene_ref ref,
						const unsigned char rgb[3]);
static int ged_draw_scene_ref_realize_view_inputs_changed(bsg_scene_ref ref,
							  void *view_ctx);
static void ged_draw_scene_ref_invalidate(bsg_scene_ref ref);
static int ged_draw_scene_ref_apply_path_state(struct ged *gedp,
					       bsg_scene_ref ref,
					       const struct db_full_path *path);
static void ged_draw_scene_ref_realize(bsg_scene_ref ref, void *view_ctx);
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
static int ged_draw_scene_ref_attach_source_state_record(
	bsg_scene_ref ref,
	const struct ged_draw_source_state_record *record);
static int ged_draw_scene_ref_database_source_record_get(
	bsg_scene_ref ref,
	struct ged_draw_database_source_record *record);
static int ged_draw_scene_ref_database_source_record_apply(
	bsg_scene_ref ref,
	const struct ged_draw_database_source_record *record);
static int ged_draw_scene_ref_mark_database_source_changed(
	bsg_scene_ref ref,
	ged_draw_stale_reason reason);
static void ged_draw_scene_ref_mark_database_source_redraw_result(
	bsg_scene_ref ref,
	int failed);
static void color_soltab_scene_ref(struct db_i *dbip, bsg_scene_ref shape_ref);
static int ged_draw_scene_ref_refresh_material_color(
	struct db_i *dbip,
	bsg_scene_ref shape_ref,
	uint64_t mater_rev);
static int ged_draw_shape_draft_commit_database_leaf_to_scene_ref(
	ged_draw_shape_draft *draft,
	bsg_scene_ref parent_ref);
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
static void ged_draw_scene_ref_realize_dispatch(
	bsg_scene_ref ref,
	void *view_ctx);
static int ged_draw_scene_ref_source_runtime_summary(
	bsg_scene_ref ref,
	struct ged_draw_source_runtime_summary *out);
static int ged_draw_source_primitive_has_lod_realize(const struct rt_db_internal *ip);
static int ged_draw_scene_ref_geometry_clear(bsg_scene_ref ref);
static int ged_draw_scene_ref_publish_bot_wireframe_line_set(
	bsg_scene_ref ref,
	const struct rt_bot_internal *bot);
static int ged_draw_scene_ref_publish_brep_wireframe_line_set(
	bsg_scene_ref ref,
	const struct rt_brep_internal *bi,
	const struct bn_tol *tol);
static int ged_draw_scene_ref_publish_bspline_wireframe_line_set(
	bsg_scene_ref ref,
	struct rt_db_internal *ip,
	const struct bn_tol *tol);
static int ged_draw_scene_ref_publish_indexed_face_set(
	bsg_scene_ref ref,
	const point_t *points,
	size_t point_count,
	const vect_t *normals,
	size_t normal_count,
	const int *indices,
	size_t index_count);
static int ged_draw_scene_ref_publish_line_set(
	bsg_scene_ref ref,
	const point_t *points,
	const int *commands,
	size_t point_count);
static int ged_draw_scene_ref_publish_point_set(
	bsg_scene_ref ref,
	const point_t *points,
	size_t point_count);
static int ged_draw_scene_ref_publish_primitive_face_set(
	bsg_scene_ref ref,
	struct rt_db_internal *ip,
	const struct bg_tess_tol *ttol,
	const struct bn_tol *tol,
	const struct rt_view_info *view_info);
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
static int ged_draw_scene_ref_geometry_publish_nmg_region(
	bsg_scene_ref ref,
	const struct nmgregion *r,
	int style);
static int ged_draw_scene_ref_update_bounds_from_geometry(
	bsg_scene_ref ref,
	int *bad_cmd);
static int ged_draw_scene_ref_update_bounds_context(
	bsg_scene_ref ref,
	void *view_ctx);
static int ged_draw_scene_ref_update_indexed_face_set(bsg_scene_ref ref,
						      const point_t *points,
						      size_t point_count,
						      const vect_t *normals,
						      size_t normal_count,
						      const int *indices,
						      size_t index_count);
static int ged_draw_scene_ref_publish_current_face_set(bsg_scene_ref ref,
						       struct rt_db_internal *ip,
						       const struct bg_tess_tol *ttol,
						       const struct bn_tol *tol,
						       const struct rt_view_info *view_info);
static int ged_draw_scene_ref_publish_current_wireframe(bsg_scene_ref ref,
							struct rt_db_internal *ip,
							const struct bg_tess_tol *ttol,
							const struct bn_tol *tol,
							void *view_ctx,
							const struct rt_view_info *view_info,
							int adaptive);
static int ged_draw_scene_ref_publish_current_nmg_region(bsg_scene_ref ref,
							 const struct nmgregion *r,
							 int style);
static int ged_draw_scene_ref_geometry_publish_nmg_model(bsg_scene_ref ref,
							 const struct model *m,
							 int style);
static bsg_scene_ref ged_draw_scene_ref_geometry_query_target(
	bsg_scene_ref ref);
static bsg_scene_ref ged_draw_scene_ref_geometry_publication_target(
	bsg_scene_ref ref);
static int ged_draw_scene_ref_bounds(bsg_scene_ref ref, point_t min, point_t max);
static int ged_draw_scene_ref_subtree_bounds(bsg_scene_ref ref,
					     vect_t *min,
					     vect_t *max,
					     int include_overlays);
static void ged_draw_scene_ref_realization_finish_view(
	bsg_scene_ref ref,
	void *view_ctx,
	void *scene_view_ctx);
static void ged_draw_scene_ref_realization_finish_current_view(
	bsg_scene_ref ref,
	void *view_ctx,
	void *scene_view_ctx);
static int ged_draw_scene_ref_set_changed(bsg_scene_ref ref, int changed);
static int ged_draw_scene_ref_copy_draw_intent(bsg_scene_ref dst,
					       bsg_scene_ref src);
static int ged_draw_scene_ref_set_non_database_source(
	bsg_scene_ref ref,
	int non_database_source);
static void ged_draw_scene_ref_realization_set_current(bsg_scene_ref ref,
						       int current);
static void ged_draw_scene_ref_realization_mark_current(bsg_scene_ref ref);
static void ged_draw_scene_ref_realization_clear_stale(bsg_scene_ref ref);
static void ged_draw_scene_ref_realization_set_roles(bsg_scene_ref ref,
						     int csg_obj,
						     int mesh_obj);
static void ged_draw_scene_ref_realization_prepare_wireframe(bsg_scene_ref ref);
static void ged_draw_scene_ref_realization_prepare_mesh(bsg_scene_ref ref);
static void ged_draw_scene_ref_realization_prepare_surface(bsg_scene_ref ref);
static int _ged_draw_scene_ref_mark_view_inputs_changed(bsg_scene_ref ref);
static int ged_draw_scene_ref_realization_prepare_view_redraw(
	bsg_scene_ref ref,
	void *view_ctx);
static int ged_draw_scene_ref_publish_current_mesh_lod(bsg_scene_ref ref);
static uint64_t ged_draw_scene_ref_material_revision(bsg_scene_ref ref);
static int ged_draw_scene_ref_set_material_revision(bsg_scene_ref ref,
						    uint64_t mater_rev);
static void *ged_draw_scene_ref_view_context(bsg_scene_ref ref);
static int ged_draw_scene_ref_draw_mat(bsg_scene_ref ref, mat_t mat);
static int ged_draw_scene_ref_draw_intent_is_overlay(bsg_scene_ref ref);
static const char *ged_draw_scene_ref_draw_intent_path(bsg_scene_ref ref);
static int ged_draw_scene_ref_draw_intent_mode(bsg_scene_ref ref,
					       int *draw_mode);
static int ged_draw_scene_ref_line_style(bsg_scene_ref ref);
static int ged_draw_scene_ref_strict_fallback(bsg_scene_ref ref);
static int ged_draw_scene_ref_realization_current(bsg_scene_ref ref);
static fastf_t ged_draw_scene_ref_realization_view_scale(bsg_scene_ref ref);
static fastf_t ged_draw_scene_ref_realization_curve_scale(bsg_scene_ref ref);
static fastf_t ged_draw_scene_ref_realization_point_scale(bsg_scene_ref ref);
static int ged_draw_scene_ref_realization_pipeline_candidate(bsg_scene_ref ref);
static struct ged *_ged_draw_scene_ref_owner_gedp(bsg_scene_ref ref);
static int _ged_draw_scene_ref_user_data_set(
	bsg_scene_ref ref,
	void *data,
	ged_draw_shape_user_data_kind kind);
static void _ged_draw_scene_ref_copy_display_state(bsg_scene_ref dst, bsg_scene_ref src);
static int _ged_draw_scene_ref_apply_aux_display_state(bsg_scene_ref dst, bsg_scene_ref src);
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
								  struct bsg_view *v);

enum {
    GED_DRAW_NMG_CNURB_INTERIOR_SAMPLES = 10,
    GED_DRAW_NMG_CNURB_SAMPLE_POINTS = GED_DRAW_NMG_CNURB_INTERIOR_SAMPLES + 2,
    GED_DRAW_NMG_SNURB_INTERIOR_SAMPLES = 10
};


struct ged_draw_shape_draft {
    struct ged *gedp;
    void *view_ctx;
    bsg_scene_ref source_ref;
    bsg_scene_ref shape_ref;
    int committed;
};


static void
_ged_draw_shape_draft_destroy_refs(bsg_scene_ref source_ref, bsg_scene_ref shape_ref)
{
    if (!ged_draw_scene_ref_is_null(source_ref) &&
	    !ged_draw_scene_ref_is_null(shape_ref))
	(void)ged_draw_scene_ref_detach(shape_ref);
    if (!ged_draw_scene_ref_is_null(shape_ref))
	ged_draw_scene_ref_release(shape_ref);
    if (!ged_draw_scene_ref_is_null(source_ref))
	ged_draw_scene_ref_release(source_ref);
}


static void
_ged_draw_shape_draft_sync_aux_geometry(ged_draw_shape_draft *draft)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->source_ref) ||
	    ged_draw_scene_ref_is_null(draft->shape_ref))
	return;

    (void)ged_draw_scene_ref_sync_aux_display_state(draft->source_ref,
	    draft->shape_ref);
}


ged_draw_shape_draft *
ged_draw_shape_draft_create_context(struct ged *gedp, void *view_ctx, int registered)
{
    bsg_scene_ref source_ref = ged_draw_scene_ref_null();
    bsg_scene_ref shape_ref = ged_draw_scene_ref_null();
    if (!ged_draw_scene_ref_create_draft_pair(gedp, view_ctx, &source_ref,
	    &shape_ref))
	return NULL;

    ged_draw_shape_draft *draft;
    BU_GET(draft, ged_draw_shape_draft);
    draft->gedp = gedp;
    draft->view_ctx = view_ctx;
    draft->source_ref = source_ref;
    draft->shape_ref = shape_ref;
    draft->committed = 0;
    if (registered && !ged_draw_shape_draft_apply_db_object_marker(shape_ref)) {
	ged_draw_shape_draft_destroy(draft);
	return NULL;
    }
    return draft;
}


void
ged_draw_shape_draft_destroy(ged_draw_shape_draft *draft)
{
    if (!draft)
	return;
    if (!draft->committed)
	_ged_draw_shape_draft_destroy_refs(draft->source_ref, draft->shape_ref);
    draft->source_ref = ged_draw_scene_ref_null();
    draft->shape_ref = ged_draw_scene_ref_null();
    BU_PUT(draft, ged_draw_shape_draft);
}


static int
_ged_draw_shape_draft_apply_path_state(ged_draw_shape_draft *draft,
				       const struct db_full_path *path)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_apply_path_state(draft->gedp, draft->shape_ref,
	    path);
}


int
ged_draw_shape_draft_publish_line_set(ged_draw_shape_draft *draft,
				      const point_t *points,
				      const int *commands,
				      size_t point_count)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_publish_line_set(draft->shape_ref, points,
	    commands, point_count);
}


int
ged_draw_shape_draft_publish_bot_wireframe_line_set(ged_draw_shape_draft *draft,
						    const struct rt_bot_internal *bot)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !bot)
	return 0;
    return ged_draw_scene_ref_publish_bot_wireframe_line_set(draft->shape_ref,
	    bot);
}


int
ged_draw_shape_draft_publish_primitive_face_set(ged_draw_shape_draft *draft,
						struct rt_db_internal *ip,
						const struct bg_tess_tol *ttol,
						const struct bn_tol *tol,
						const struct rt_view_info *view_info)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !ip)
	return 0;
    return ged_draw_scene_ref_publish_primitive_face_set(draft->shape_ref, ip,
	    ttol, tol, view_info);
}


int
ged_draw_shape_draft_publish_brep_wireframe_line_set(ged_draw_shape_draft *draft,
						     const struct rt_brep_internal *bi,
						     const struct bn_tol *tol)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !bi)
	return 0;
    return ged_draw_scene_ref_publish_brep_wireframe_line_set(draft->shape_ref,
	    bi, tol);
}


int
ged_draw_shape_draft_publish_poly_wireframe_line_set(ged_draw_shape_draft *draft,
						     const struct rt_pg_internal *pg)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !pg)
	return 0;
    return ged_draw_scene_ref_publish_poly_wireframe_line_set(draft->shape_ref,
	    pg);
}


int
ged_draw_shape_draft_publish_nmg_region(ged_draw_shape_draft *draft,
					const struct nmgregion *r,
					int style)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !r)
	return 0;
    return ged_draw_scene_ref_geometry_publish_nmg_region(draft->shape_ref,
	    r, style);
}


static int
_ged_draw_shape_draft_refresh_scene_bounds(ged_draw_shape_draft *draft)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    return ged_draw_scene_ref_update_bounds_context(draft->shape_ref,
	    draft->view_ctx);
}


int
ged_draw_shape_draft_apply_known_bounds(ged_draw_shape_draft *draft,
					const point_t min,
					const point_t max)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;
    if (!min || !max)
	return 0;

    vect_t center;
    center[X] = (min[X] + max[X]) * 0.5;
    center[Y] = (min[Y] + max[Y]) * 0.5;
    center[Z] = (min[Z] + max[Z]) * 0.5;
    ged_draw_scene_ref_set_draw_center(draft->shape_ref, center);

    fastf_t size = max[X] - min[X];
    V_MAX(size, max[Y] - min[Y]);
    V_MAX(size, max[Z] - min[Z]);
    if (!ged_draw_scene_ref_set_draw_size(draft->shape_ref, size))
	return 0;

    return ged_draw_scene_ref_set_bounds(draft->shape_ref, min, max, 1);
}


int
ged_draw_shape_draft_publish_primitive_wireframe(ged_draw_shape_draft *draft,
						 struct rt_db_internal *ip,
						 const struct bg_tess_tol *ttol,
						 const struct bn_tol *tol,
						 void *view_ctx,
						 int adaptive)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return -1;
    struct rt_view_info view_info;
    const struct rt_view_info *view_infop = NULL;
    if (adaptive && view_ctx) {
	rt_view_context_info_get(&view_info, view_ctx);
	view_infop = &view_info;
    }
    int status = ged_draw_scene_ref_publish_primitive_wireframe(draft->shape_ref, ip,
	    ttol, tol, view_ctx, view_infop, adaptive);
    if (status >= 0 && !adaptive)
	(void)ged_draw_scene_ref_update_bounds_from_geometry(draft->shape_ref,
		NULL);
    return status;
}


int
ged_draw_shape_draft_apply_tree_result_state(ged_draw_shape_draft *draft,
					     int dashed,
					     int has_region,
					     int region_id,
					     int aircode,
					     int los,
					     int material_id,
					     const struct ged_draw_appearance_settings *settings)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !settings)
	return 0;
    (void)_ged_draw_shape_draft_refresh_scene_bounds(draft);
    if (!ged_draw_scene_ref_set_highlighted(draft->shape_ref, 0))
	return 0;
    if (!ged_draw_scene_ref_set_line_style(draft->shape_ref, dashed))
	return 0;
    if (!ged_draw_scene_ref_set_legacy_eval_flag(draft->shape_ref, 0))
	return 0;
    if (!has_region)
	return ged_draw_scene_ref_apply_display_settings(draft->shape_ref,
		settings);
    if (!ged_draw_scene_ref_set_legacy_region_id(draft->shape_ref,
	    region_id))
	return 0;
    if (!ged_draw_shape_context_apply_registry_region(draft->gedp,
	    ged_draw_scene_ref_context(draft->shape_ref), region_id, aircode,
	    los, material_id))
	return 0;
    return ged_draw_scene_ref_apply_display_settings(draft->shape_ref,
	    settings);
}


static int
_ged_draw_shape_draft_apply_name(ged_draw_shape_draft *draft, const char *name)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !name)
	return 0;
    if (!ged_draw_scene_ref_set_name(draft->shape_ref, name))
	return 0;
    if (!ged_draw_scene_ref_is_null(draft->source_ref)) {
	struct bu_vls source_name = BU_VLS_INIT_ZERO;
	bu_vls_sprintf(&source_name, "%s:source", name);
	(void)ged_draw_scene_ref_set_name(draft->source_ref,
		bu_vls_cstr(&source_name));
	bu_vls_free(&source_name);
    }
    return 1;
}


int
ged_draw_shape_draft_apply_tree_legacy_color(
	ged_draw_shape_draft *draft,
	const unsigned char wireframe_color_override[3],
	const struct db_tree_state *tsp)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;

    unsigned char bcolor[3] = {255, 0, 0};
    int user_color = 0;
    int default_color = 0;

    if (wireframe_color_override) {
	user_color = 1;
	bcolor[RED] = wireframe_color_override[RED];
	bcolor[GRN] = wireframe_color_override[GRN];
	bcolor[BLU] = wireframe_color_override[BLU];
    } else if (tsp) {
	if (tsp->ts_mater.ma_color_valid) {
	    bcolor[RED] = tsp->ts_mater.ma_color[RED] * 255.0;
	    bcolor[GRN] = tsp->ts_mater.ma_color[GRN] * 255.0;
	    bcolor[BLU] = tsp->ts_mater.ma_color[BLU] * 255.0;
	} else {
	    default_color = 1;
	}
    }

    if (!ged_draw_scene_ref_set_legacy_color_info(draft->shape_ref,
	    bcolor, user_color, default_color))
	return 0;
    color_soltab_scene_ref(tsp ? tsp->ts_dbip : NULL, draft->shape_ref);
    return 1;
}


static int
_ged_draw_shape_draft_attach_source_inputs(
	ged_draw_shape_draft *draft,
	struct db_i *dbip,
	const struct db_full_path *pathp,
	const struct bn_tol *tol,
	const struct bg_tess_tol *ttol)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return 0;

    struct ged_draw_source_state_record record;
    memset(&record, 0, sizeof(record));
    record.dbip = dbip;
    record.fullpath = pathp;
    record.tol = tol;
    record.ttol = ttol;
    record.stale_reason = GED_DRAW_STALE_NONE;
    return ged_draw_scene_ref_attach_source_state_record(draft->shape_ref,
	    &record);
}


int
ged_draw_shape_draft_apply_path_source_state(
	ged_draw_shape_draft *draft,
	struct db_i *dbip,
	const struct db_full_path *pathp,
	const struct bn_tol *tol,
	const struct bg_tess_tol *ttol,
	int has_draw_mat,
	const mat_t draw_mat,
	const char *display_name)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !pathp)
	return 0;
    if (display_name && !_ged_draw_shape_draft_apply_name(draft, display_name))
	return 0;
    if (!_ged_draw_shape_draft_apply_path_state(draft, pathp))
	return 0;
    if (!_ged_draw_shape_draft_attach_source_inputs(draft, dbip, pathp, tol,
	    ttol))
	return 0;
    if (has_draw_mat &&
	    (!draw_mat || !ged_draw_scene_ref_set_draw_mat(draft->shape_ref,
		draw_mat)))
	return 0;
    return 1;
}


static int
_ged_draw_shape_draft_apply_material_rgb(ged_draw_shape_draft *draft,
					 const unsigned char rgb[3])
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !rgb)
	return 0;
    if (!ged_draw_scene_ref_set_color(draft->shape_ref, rgb))
	return 0;
    ged_draw_scene_ref_set_material_rgb(draft->shape_ref, rgb);
    return 1;
}


static int
_ged_draw_shape_draft_apply_material_rgb_revision(ged_draw_shape_draft *draft,
						  const unsigned char rgb[3])
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) || !rgb)
	return 0;
    if (!_ged_draw_shape_draft_apply_material_rgb(draft, rgb))
	return 0;
    return ged_draw_scene_ref_bump_appearance_revision(draft->shape_ref);
}


int
ged_draw_shape_draft_apply_late_display_state(
	ged_draw_shape_draft *draft,
	const struct db_full_path *path,
	int line_style,
	const struct ged_draw_appearance_settings *settings,
	const unsigned char material_rgb[3],
	int highlighted)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) ||
	    !path || !settings || !material_rgb)
	return 0;
    if (!_ged_draw_shape_draft_refresh_scene_bounds(draft))
	return 0;
    if (!_ged_draw_shape_draft_apply_path_state(draft, path))
	return 0;
    if (!ged_draw_scene_ref_set_line_style(draft->shape_ref, line_style))
	return 0;
    if (!ged_draw_scene_ref_apply_display_settings(draft->shape_ref,
	    settings))
	return 0;
    if (!_ged_draw_shape_draft_apply_material_rgb(draft, material_rgb))
	return 0;
    return ged_draw_scene_ref_set_highlighted(draft->shape_ref, highlighted);
}


int
ged_draw_shape_draft_apply_evaluated_path_display(
	ged_draw_shape_draft *draft,
	const struct ged_draw_appearance_settings *settings,
	const unsigned char material_rgb[3],
	int has_transform,
	const mat_t transform,
	int is_subtraction)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) ||
	    !settings || !material_rgb)
	return 0;

    if (!_ged_draw_shape_draft_apply_material_rgb(draft, material_rgb))
	return 0;
    if (has_transform && transform &&
	    !ged_draw_scene_ref_set_transform(draft->shape_ref, transform))
	return 0;
    if (!settings->draw_solid_lines_only &&
	    !ged_draw_scene_ref_set_line_style(draft->shape_ref,
		is_subtraction ? 1 : 0))
	return 0;
    if (settings->draw_non_subtract_only && is_subtraction &&
	    !ged_draw_scene_ref_set_visible(draft->shape_ref, 0))
	return 0;
    return ged_draw_scene_ref_apply_display_settings(draft->shape_ref,
	    settings);
}


int
ged_draw_shape_draft_apply_database_leaf_display(
	ged_draw_shape_draft *draft,
	const struct ged_draw_appearance_settings *settings,
	int bool_op,
	const unsigned char material_rgb[3],
	int has_draw_size,
	fastf_t draw_size)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) ||
	    !material_rgb)
	return 0;
    if (settings &&
	    !ged_draw_scene_ref_apply_display_settings(draft->shape_ref,
		settings))
	return 0;
    if ((!settings || !settings->draw_solid_lines_only) &&
	    !ged_draw_scene_ref_set_line_style(draft->shape_ref,
		(bool_op == 4) ? 1 : 0))
	return 0;
    if (!_ged_draw_shape_draft_apply_material_rgb_revision(draft,
	    material_rgb))
	return 0;
    if (has_draw_size &&
	    !ged_draw_scene_ref_set_draw_size(draft->shape_ref, draw_size))
	return 0;
    return 1;
}


ged_draw_shape_ref
ged_draw_shape_draft_commit_to_group(ged_draw_shape_draft *draft,
				     ged_draw_group_ref group_ref)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref))
	return GED_DRAW_SHAPE_REF_NULL;

    _ged_draw_shape_draft_sync_aux_geometry(draft);

    if (!ged_draw_group_ref_append_scene_ref(draft->gedp, group_ref, draft->shape_ref)) {
	ged_draw_shape_draft_destroy(draft);
	return GED_DRAW_SHAPE_REF_NULL;
    }
    ged_draw_shape_ref ref =
	ged_draw_shape_ref_from_scene_ref(draft->gedp, draft->shape_ref);
    if (ged_draw_shape_ref_is_null(ref)) {
	ged_draw_shape_draft_destroy(draft);
	return GED_DRAW_SHAPE_REF_NULL;
    }
    ref.scene_revision = 0;

    draft->committed = 1;
    ged_draw_shape_draft_destroy(draft);
    return ref;
}


static int
ged_draw_shape_draft_commit_database_leaf_to_scene_ref(
	ged_draw_shape_draft *draft,
	bsg_scene_ref parent_ref)
{
    if (!draft || ged_draw_scene_ref_is_null(draft->shape_ref) ||
	    ged_draw_scene_ref_is_null(parent_ref))
	return 0;

    _ged_draw_shape_draft_sync_aux_geometry(draft);

    if (!ged_draw_scene_ref_append_child(parent_ref, draft->source_ref)) {
	ged_draw_shape_draft_destroy(draft);
	return 0;
    }
    draft->committed = 1;
    ged_draw_shape_draft_destroy(draft);
    return 1;
}


static int
ged_draw_shape_draft_commit_tree_result(
	ged_draw_shape_draft *draft,
	struct ged *gedp,
	ged_draw_group_ref group_ref,
	const struct ged_draw_appearance_settings *settings,
	int dashflag,
	const struct db_full_path *pathp,
	struct db_tree_state *tsp,
	const unsigned char wireframe_color_override[3])
{
    if (!draft || !gedp || !settings)
	return 0;

    (void)gedp;
    (void)pathp;

    if (!ged_draw_shape_draft_apply_tree_legacy_color(draft,
	    wireframe_color_override, tsp))
	return 0;
    if (!ged_draw_shape_draft_apply_tree_result_state(draft, dashflag,
	    tsp ? 1 : 0, tsp ? tsp->ts_regionid : 0,
	    tsp ? tsp->ts_aircode : 0, tsp ? tsp->ts_los : 0,
	    tsp ? tsp->ts_gmater : 0, settings))
	return 0;

    bu_semaphore_acquire(RT_SEM_MODEL);
    (void)ged_draw_shape_draft_commit_to_group(draft, group_ref);
    bu_semaphore_release(RT_SEM_MODEL);

    return 1;
}


static int
ged_draw_shape_draft_prepare_tree_source_state(
	ged_draw_shape_draft *draft,
	struct ged *gedp,
	const struct db_full_path *pathp,
	struct db_tree_state *tsp)
{
    if (!draft || !gedp || !pathp)
	return 0;

    return ged_draw_shape_draft_apply_path_source_state(draft,
	    (tsp && tsp->ts_dbip) ? tsp->ts_dbip : gedp->dbip,
	    pathp, tsp ? tsp->ts_tol : NULL, tsp ? tsp->ts_ttol : NULL,
	    tsp ? 1 : 0, tsp ? tsp->ts_mat : NULL, NULL);
}


static int
_ged_draw_path_is_subtraction(struct ged *gedp,
			      const struct db_full_path *path)
{
    if (!gedp || !gedp->dbip || !path)
	return 0;

    int op = db_fp_op(path, gedp->dbip, 0);
    return (op == DB_OP_SUBTRACT || op == OP_SUBTRACT) ? 1 : 0;
}


ged_draw_shape_ref
ged_draw_create_evaluated_path_shape_ref(
	struct ged *gedp,
	void *view_ctx,
	const char *path,
	const struct ged_draw_appearance_settings *settings)
{
    if (!gedp || !gedp->dbip || !view_ctx || !path || !path[0] || !settings)
	return GED_DRAW_SHAPE_REF_NULL;

    struct db_full_path dfp;
    db_full_path_init(&dfp);
    if (db_string_to_path(&dfp, gedp->dbip, path) != 0)
	return GED_DRAW_SHAPE_REF_NULL;

    ged_draw_group_ref group_ref =
	ged_draw_group_ref_lookup_or_create(gedp, &dfp);
    if (ged_draw_group_ref_is_null(group_ref)) {
	db_free_full_path(&dfp);
	return GED_DRAW_SHAPE_REF_NULL;
    }
    ged_draw_group_ref_set_appearance_settings(gedp, group_ref, settings);

    ged_draw_shape_draft *draft = ged_draw_shape_draft_create_context(gedp,
	    view_ctx, 1);
    if (!draft) {
	db_free_full_path(&dfp);
	return GED_DRAW_SHAPE_REF_NULL;
    }

    char *name = db_path_to_string(&dfp);
    ged_draw_shape_draft_apply_path_source_state(draft, gedp->dbip, &dfp,
	    NULL, NULL, 0, NULL,
	    name ? ged_draw_dbpath_skip_lead_slash(name) : NULL);

    struct bu_color c;
    unsigned char color[3] = {255, 0, 0};
    db_full_path_color(&c, &dfp, gedp->dbip);
    bu_color_to_rgb_chars(&c, color);
    if (settings->color_override) {
	color[0] = settings->color[0];
	color[1] = settings->color[1];
	color[2] = settings->color[2];
    }

    mat_t node_mat;
    MAT_IDN(node_mat);
    int has_node_mat = db_path_to_mat(gedp->dbip, &dfp, node_mat, 0) ? 1 : 0;

    point_t bmin, bmax;
    VSETALL(bmin, INFINITY);
    VSETALL(bmax, -INFINITY);
    const char *bounds_path = name ? ged_draw_dbpath_skip_lead_slash(name) :
	path;
    const char *bounds_argv[1] = {bounds_path};
    if (rt_obj_bounds(gedp->ged_result_str, gedp->dbip, 1, bounds_argv, 0,
	    bmin, bmax) != BRLCAD_ERROR) {
	ged_draw_shape_draft_apply_known_bounds(draft, bmin, bmax);
    }

    int is_subtraction = _ged_draw_path_is_subtraction(gedp, &dfp);
    ged_draw_shape_draft_apply_evaluated_path_display(draft, settings,
	    color, has_node_mat, node_mat, is_subtraction);

    ged_draw_shape_ref ref =
	ged_draw_shape_draft_commit_to_group(draft, group_ref);
    if (name)
	bu_free(name, "draw evaluated path name");
    db_free_full_path(&dfp);
    return ref;
}


int
ged_draw_append_tree_shape_to_group(
	struct ged *gedp,
	void *view_ctx,
	ged_draw_group_ref group_ref,
	const struct ged_draw_appearance_settings *settings,
	const struct db_full_path *pathp,
	struct db_tree_state *tsp,
	struct rt_db_internal *ip)
{
    if (!gedp || !view_ctx || ged_draw_group_ref_is_null(group_ref) ||
	    !settings || !pathp || !tsp || !ip || !ip->idb_meth)
	return 0;

    point_t min, max;
    VSETALL(min, INFINITY);
    VSETALL(max, -INFINITY);

    ged_draw_shape_draft *draft = ged_draw_shape_draft_create_context(gedp,
	    view_ctx, 1);
    if (!draft)
	return 0;

    if (!ged_draw_shape_draft_prepare_tree_source_state(draft, gedp, pathp,
	    tsp)) {
	ged_draw_shape_draft_destroy(draft);
	return 0;
    }

    int have_bounds = 0;
    if (ip->idb_meth->ft_bbox &&
	    ip->idb_meth->ft_bbox(ip, &min, &max, tsp->ts_tol) >= 0) {
	ged_draw_shape_draft_apply_known_bounds(draft, min, max);
	have_bounds = 1;
    }

    if (!have_bounds) {
	int plot_status = ged_draw_shape_draft_publish_primitive_wireframe(draft,
		ip, tsp->ts_ttol, tsp->ts_tol, NULL, 0);
	if (plot_status < 0) {
	    if (pathp && DB_FULL_PATH_CUR_DIR(pathp)) {
		bu_log("%s: plot failure\n",
			DB_FULL_PATH_CUR_DIR(pathp)->d_namep);
	    } else {
		bu_log("plot failure - invalid path\n");
	    }
	    ged_draw_shape_draft_destroy(draft);
	    return 0;
	}
    }

    int dashed = settings->draw_solid_lines_only ? 0 :
	((tsp->ts_sofar & (TS_SOFAR_MINUS|TS_SOFAR_INTER)) ? 1 : 0);
    const unsigned char *tree_color_override = NULL;
    unsigned char obj_color[3];
    struct db_tree_state color_state;
    const struct db_tree_state *commit_tsp = tsp;

    if (ip->idb_type == ID_GRIP) {
	if (settings->color_override) {
	    tree_color_override = settings->color;
	} else {
	    color_state = *tsp;
	    color_state.ts_mater.ma_color[RED] = 0;
	    color_state.ts_mater.ma_color[GRN] = 128;
	    color_state.ts_mater.ma_color[BLU] = 128;
	    commit_tsp = &color_state;
	}
    } else if (settings->color_override) {
	tree_color_override = settings->color;
    } else {
	const char *attr_color = bu_avs_get(&ip->idb_avs,
		db5_standard_attribute(ATTR_COLOR));
	if (attr_color) {
	    int i;
	    int color[3];
	    int color_cnt = sscanf(attr_color, "%3i%*c%3i%*c%3i", color+0,
		    color+1, color+2);
	    if (color_cnt == 3 && color[0] >= 0 && color[1] >= 0 &&
		    color[2] >= 0) {
		for (i = 0; i < 3; i++) {
		    if (color[i] > 255)
			color[i] = 255;
		}
		obj_color[RED] = (unsigned char)color[RED];
		obj_color[GRN] = (unsigned char)color[GRN];
		obj_color[BLU] = (unsigned char)color[BLU];
		tree_color_override = obj_color;
	    }
	}
    }

    return ged_draw_shape_draft_commit_tree_result(draft, gedp, group_ref,
	    settings, dashed, pathp, (struct db_tree_state *)commit_tsp,
	    tree_color_override);
}


static ged_draw_shape_draft *
_ged_draw_tree_draft_create(struct ged *gedp, void *view_ctx)
{
    if (!gedp || !view_ctx)
	return NULL;

    ged_draw_shape_draft *draft = ged_draw_shape_draft_create_context(gedp,
	    view_ctx, 1);
    if (!draft)
	return NULL;

    return draft;
}


int
ged_draw_add_tree_line_set_to_group(
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
	size_t point_count)
{
    ged_draw_shape_draft *draft = _ged_draw_tree_draft_create(gedp, view_ctx);
    if (!draft)
	return 0;

    if (!ged_draw_shape_draft_prepare_tree_source_state(draft, gedp, pathp,
	    tsp)) {
	ged_draw_shape_draft_destroy(draft);
	return 0;
    }

    if (!ged_draw_shape_draft_publish_line_set(draft, points, commands,
	    point_count)) {
	ged_draw_shape_draft_destroy(draft);
	return 0;
    }

    return ged_draw_shape_draft_commit_tree_result(draft, gedp, group_ref,
	    settings, dashflag, pathp, tsp, wireframe_color_override);
}


int
ged_draw_add_tree_nmg_region_to_group(
	struct ged *gedp,
	void *view_ctx,
	ged_draw_group_ref group_ref,
	const struct ged_draw_appearance_settings *settings,
	int dashflag,
	const struct db_full_path *pathp,
	struct db_tree_state *tsp,
	const unsigned char wireframe_color_override[3],
	const struct nmgregion *r,
	int style)
{
    if (!r)
	return 0;

    ged_draw_shape_draft *draft = _ged_draw_tree_draft_create(gedp, view_ctx);
    if (!draft)
	return 0;

    if (!ged_draw_shape_draft_prepare_tree_source_state(draft, gedp, pathp,
	    tsp)) {
	ged_draw_shape_draft_destroy(draft);
	return 0;
    }

    if (!ged_draw_shape_draft_publish_nmg_region(draft, r, style)) {
	ged_draw_shape_draft_destroy(draft);
	return 0;
    }

    return ged_draw_shape_draft_commit_tree_result(draft, gedp, group_ref,
	    settings, dashflag, pathp, tsp, wireframe_color_override);
}


int
ged_draw_add_tree_primitive_face_set_to_group(
	struct ged *gedp,
	void *view_ctx,
	ged_draw_group_ref group_ref,
	const struct ged_draw_appearance_settings *settings,
	int dashflag,
	const struct db_full_path *pathp,
	struct db_tree_state *tsp,
	const unsigned char wireframe_color_override[3],
	const struct rt_db_internal *ip,
	int force_failure)
{
    if (!ip)
	return 0;

    ged_draw_shape_draft *draft = _ged_draw_tree_draft_create(gedp, view_ctx);
    if (!draft)
	return 0;

    if (!ged_draw_shape_draft_prepare_tree_source_state(draft, gedp, pathp,
	    tsp)) {
	ged_draw_shape_draft_destroy(draft);
	return 0;
    }

    int ok = 0;
    if (!force_failure && ip->idb_meth && ip->idb_meth->ft_indexed_face_set) {
	struct rt_view_info view_info;
	ged_draw_view_context_info_get(&view_info, view_ctx);
	ok = ged_draw_shape_draft_publish_primitive_face_set(draft,
		(struct rt_db_internal *)ip, tsp ? tsp->ts_ttol : NULL,
		tsp ? tsp->ts_tol : NULL, &view_info);
    }

    if (!ok) {
	ged_draw_shape_draft_destroy(draft);
	return 0;
    }

    return ged_draw_shape_draft_commit_tree_result(draft, gedp, group_ref,
	    settings, dashflag, pathp, tsp, wireframe_color_override);
}


int
ged_draw_add_tree_primitive_wireframe_to_group(
	struct ged *gedp,
	void *view_ctx,
	ged_draw_group_ref group_ref,
	const struct ged_draw_appearance_settings *settings,
	int dashflag,
	const struct db_full_path *pathp,
	struct db_tree_state *tsp,
	const unsigned char wireframe_color_override[3],
	const struct rt_db_internal *ip)
{
    if (!ip)
	return 0;

    ged_draw_shape_draft *draft = _ged_draw_tree_draft_create(gedp, view_ctx);
    if (!draft)
	return 0;

    if (!ged_draw_shape_draft_prepare_tree_source_state(draft, gedp, pathp,
	    tsp)) {
	ged_draw_shape_draft_destroy(draft);
	return 0;
    }

    int ok = 0;
    switch (ip->idb_type) {
	case ID_BOT:
	    {
		const struct rt_bot_internal *bot =
		    (const struct rt_bot_internal *)ip->idb_ptr;
		RT_BOT_CK_MAGIC(bot);
		ok = ged_draw_shape_draft_publish_bot_wireframe_line_set(draft,
			bot);
	    }
	    break;
	case ID_POLY:
	    {
		const struct rt_pg_internal *pg =
		    (const struct rt_pg_internal *)ip->idb_ptr;
		RT_PG_CK_MAGIC(pg);
		ok = ged_draw_shape_draft_publish_poly_wireframe_line_set(draft,
			pg);
	    }
	    break;
	case ID_BREP:
	    {
		const struct rt_brep_internal *bi =
		    (const struct rt_brep_internal *)ip->idb_ptr;
		RT_BREP_CK_MAGIC(bi);
		ok = ged_draw_shape_draft_publish_brep_wireframe_line_set(draft,
			bi, tsp ? tsp->ts_tol : NULL);
	    }
	    break;
	default:
	    break;
    }

    if (!ok) {
	ged_draw_shape_draft_destroy(draft);
	return 0;
    }

    return ged_draw_shape_draft_commit_tree_result(draft, gedp, group_ref,
	    settings, dashflag, pathp, tsp, wireframe_color_override);
}


static const char *
_ged_draw_source_view_context_name(const void *view_ctx)
{
    if (!view_ctx)
	return "NULL";

    const char *name = ged_view_context_name_get(view_ctx);
    return name ? name : "";
}


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

	int keep_going = cb(&detail.record, userdata);
	_ged_draw_view_export_detail_free(&detail);
	if (!keep_going)
	    break;
    }

    bsg_export_result_free(result);
}


static const char *
ged_draw_scene_ref_name(bsg_scene_ref ref)
{
    return bsg_scene_name(ref);
}


static const struct db_full_path *
ged_draw_scene_ref_fullpath(bsg_scene_ref ref)
{
    ged_draw_shape_state *bd = ged_draw_shape_state_get_scene_ref(ref);
    if (!bd || bd->s_fullpath.fp_len <= 0)
        return NULL;
    return &bd->s_fullpath;
}


static struct directory *
ged_draw_scene_ref_leaf_dp(bsg_scene_ref ref)
{
    ged_draw_shape_state *bd = ged_draw_shape_state_get_scene_ref(ref);
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
    return ged_draw_shape_state_get_scene_ref(ref);
}


static void
_ged_draw_scene_ref_set_source_ref(bsg_scene_ref ref, bsg_scene_ref source_ref)
{
    ged_draw_shape_state *shape_data = ged_draw_scene_ref_shape_state(ref);
    if (shape_data)
	shape_data->source_ref = source_ref;
}


bsg_scene_ref
ged_draw_scene_ref_null(void)
{
    return bsg_scene_ref_null();
}


bsg_scene_ref
ged_draw_scene_ref_from_context(void *scene_ctx)
{
    bsg_scene_ref ref = bsg_scene_ref_null();
    ref.opaque = scene_ctx;
    return ref;
}


void *
ged_draw_scene_ref_context(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? NULL : ref.opaque;
}


int
ged_draw_scene_ref_is_null(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref);
}


int
ged_draw_scene_ref_equal(bsg_scene_ref a, bsg_scene_ref b)
{
    return bsg_scene_ref_equal(a, b);
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

    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
    if (shape_data && !ged_draw_scene_ref_is_null(shape_data->source_ref))
	return shape_data->source_ref;

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


static int
ged_draw_scene_ref_database_source_record_get(
	bsg_scene_ref ref,
	struct ged_draw_database_source_record *record)
{
    if (!record)
	return 0;

    memset(record, 0, sizeof(*record));
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

    return bsg_database_source_record_apply(source, &bsg_record);
}


struct ged_draw_source_revision_record {
    uint64_t source_revision;
    uint64_t inputs_revision;
    uint64_t realized_source_revision;
    uint64_t realized_inputs_revision;
    ged_draw_stale_reason stale_reason;
};


static void
ged_draw_scene_ref_database_source_sync(bsg_scene_ref ref,
					const struct ged_draw_source_state *source_data,
					const ged_draw_shape_state *shape_data);


static void
ged_draw_source_revision_record_from_state(
	struct ged_draw_source_revision_record *record,
	const struct ged_draw_source_state *source_data,
	const ged_draw_shape_state *shape_data)
{
    if (!record)
	return;

    memset(record, 0, sizeof(*record));
    if (shape_data) {
	record->source_revision = shape_data->source_revision;
	record->inputs_revision = shape_data->inputs_revision;
	record->realized_source_revision =
	    shape_data->realized_source_revision;
	record->realized_inputs_revision =
	    shape_data->realized_inputs_revision;
	record->stale_reason = shape_data->stale_reason;
	return;
    }

    if (source_data) {
	record->source_revision = source_data->source_revision;
	record->inputs_revision = source_data->inputs_revision;
	record->realized_source_revision =
	    source_data->realized_source_revision;
	record->realized_inputs_revision =
	    source_data->realized_inputs_revision;
	record->stale_reason = source_data->stale_reason;
    }
}


static void
ged_draw_scene_ref_database_source_sync(bsg_scene_ref ref,
					const struct ged_draw_source_state *source_data,
					const ged_draw_shape_state *shape_data)
{
    if (ged_draw_scene_ref_is_null(ref))
	return;

    const struct db_full_path *fp = NULL;
    char *path = NULL;
    bsg_scene_ref source_ref = ged_draw_scene_ref_source_owner(ref);
    struct ged_draw_database_source_record record = {0};

    (void)ged_draw_scene_ref_database_source_record_get(source_ref, &record);

    if (shape_data)
	fp = &shape_data->s_fullpath;
    else if (source_data)
	fp = source_data->fp;

    if (fp)
	path = db_path_to_string(fp);

    record.database_path = path ? path : "";
    struct ged_draw_scene_draw_state_summary draw_summary;
    if (ged_draw_scene_ref_draw_state_summary(ref, &draw_summary))
	record.draw_mode = draw_summary.draw_mode;
    record.material_policy = fp ?
	GED_DRAW_DATABASE_SOURCE_MATERIAL_DATABASE :
	GED_DRAW_DATABASE_SOURCE_MATERIAL_INHERIT;

    if (shape_data || source_data) {
	struct ged_draw_source_revision_record revisions;
	ged_draw_source_revision_record_from_state(&revisions, source_data,
		shape_data);
	record.source_revision = revisions.source_revision;
	record.inputs_revision = revisions.inputs_revision;
	record.realized_source_revision = revisions.realized_source_revision;
	record.realized_inputs_revision = revisions.realized_inputs_revision;
	record.stale_reason = revisions.stale_reason;
	if (shape_data)
	    record.realization_identity = shape_data->path_hash;
    }

    (void)ged_draw_scene_ref_database_source_record_apply(source_ref,
	    &record);
    if (!ged_draw_scene_ref_equal(source_ref, ref))
	(void)ged_draw_scene_ref_database_source_record_apply(ref, &record);

    if (path)
	bu_free(path, "db_path_to_string");
}


void
ged_draw_source_data_free(void *data)
{
    struct ged_draw_source_state *d = (struct ged_draw_source_state *)data;
    if (!d)
	return;
    if (d->fp) {
	db_free_full_path(d->fp);
	bu_free(d->fp, "ged draw source full path");
	d->fp = NULL;
    }
    if (d->rt_mesh_lod) {
	rt_mesh_lod_destroy(d->rt_mesh_lod);
	d->rt_mesh_lod = NULL;
    }
    BU_PUT(d, struct ged_draw_source_state);
}


static void
ged_draw_scene_ref_source_data_set(bsg_scene_ref ref,
				   struct ged_draw_source_state *data)
{
    if (ged_draw_scene_ref_is_null(ref) || !data)
	return;
    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
    if (!shape_data) {
	ged_draw_source_data_free(data);
	return;
    }
    if (shape_data->source_data && shape_data->source_data != data)
	ged_draw_source_data_free(shape_data->source_data);
    shape_data->source_data = data;
    ged_draw_scene_ref_database_source_sync(ref, data, shape_data);
}


static int
ged_draw_scene_ref_attach_source_state_record(
	bsg_scene_ref ref,
	const struct ged_draw_source_state_record *record)
{
    if (ged_draw_scene_ref_is_null(ref) || !record || !record->dbip ||
	    !record->fullpath || record->fullpath->fp_len <= 0)
	return 0;

    struct ged_draw_source_state *data;
    BU_GET(data, struct ged_draw_source_state);
    memset(data, 0, sizeof(*data));

    BU_GET(data->fp, struct db_full_path);
    db_full_path_init(data->fp);
    db_dup_full_path(data->fp, record->fullpath);

    data->dbip = record->dbip;
    data->tol = record->tol;
    data->ttol = record->ttol;
    if (!data->tol || !data->ttol) {
	struct rt_wdb *wdbp = wdb_dbopen(data->dbip, RT_WDB_TYPE_DB_DEFAULT);
	if (wdbp) {
	    if (!data->tol)
		data->tol = &wdbp->wdb_tol;
	    if (!data->ttol)
		data->ttol = &wdbp->wdb_ttol;
	}
    }
    data->stale_reason = record->stale_reason;

    ged_draw_scene_ref_source_data_set(ref, data);
    return ged_draw_shape_state_get_scene_ref(ref) ? 1 : 0;
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


int
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
    ged_draw_scene_context_highlight_free_cb(ged_draw_scene_ref_context(ref));
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


static void *
ged_draw_scene_ref_parent_context(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return NULL;

    bsg_scene_ref parent_ref = ged_draw_scene_ref_parent(ref);
    if (ged_draw_scene_ref_is_null(parent_ref))
	return NULL;

    return ged_draw_scene_ref_context(parent_ref);
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


static int
ged_draw_scene_ref_append_source_owner_to_group(bsg_scene_ref group_ref,
						bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(group_ref) ||
	    ged_draw_scene_ref_is_null(ref))
	return 0;

    bsg_scene_ref attach_ref = ged_draw_scene_ref_source_owner(ref);
    if (ged_draw_scene_ref_is_null(attach_ref))
	attach_ref = ref;

    return ged_draw_scene_ref_append_child(group_ref, attach_ref);
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


bsg_scene_ref
ged_draw_view_context_group_create(void *view_ctx, const char *name)
{
    if (!view_ctx || !name)
	return bsg_scene_ref_null();

    return bsg_scene_group_create((struct bsg_view *)view_ctx, name);
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

    bsg_scene_free_group(ref);
}


static void
ged_draw_scene_ref_free_group_contents(bsg_scene_ref ref)
{
    if (ged_draw_scene_ref_is_null(ref))
	return;

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
	ged_draw_scene_ref_free_group_contents(child_ref);
	(void)ged_draw_scene_ref_detach(child_ref);
	ged_draw_scene_ref_release(child_ref);
	cleared++;
    }

    _ged_draw_clear_scope_snapshot_free(&snap);
    return cleared;
}


int
ged_draw_scene_context_clear_scope_children(void *scene_ctx)
{
    if (!scene_ctx)
	return 0;

    return ged_draw_scene_ref_clear_scope_children(
	    ged_draw_scene_ref_from_context(scene_ctx));
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


int
ged_draw_scene_ref_erase_groups_by_name(bsg_scene_ref parent_ref,
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
ged_draw_scene_ref_clear_db_groups_in_scope(struct ged *gedp, void *view_ctx)
{
    if (!gedp)
	return 0;

    bsg_scene_ref base_ref = ged_draw_scene_ref_active_scope(gedp, view_ctx,
	    0, 0);
    if (ged_draw_scene_ref_is_null(base_ref))
	return 0;

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
    return cleared;
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


int
ged_draw_scene_ref_erase_path_at_base(struct ged *gedp,
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

    int erased = ged_draw_scene_ref_erase_matching_group_path_or_nested(gedp,
	    path, base_ref, found_subpath ? &subpath : NULL, view_ctx, mode,
	    apply_scope);
    if (found_subpath) {
	if (gedp->i->ged_gdp->gd_draw_rev == erase_rev0)
	    erased += ged_draw_scene_ref_erase_shapes_by_db_subpath(gedp,
		    &subpath, base_ref, view_ctx, mode);
	db_free_full_path(&subpath);
    }

    return erased;
}


int
ged_draw_scene_ref_erase_path_prefix_at_base(struct ged *gedp,
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
ged_draw_scene_ref_erase_path_in_active_scope(struct ged *gedp,
					      const char *path,
					      void *view_ctx,
					      int mode)
{
    if (!gedp || !path)
	return 0;

    bsg_scene_ref base_ref = ged_draw_scene_ref_active_scope(gedp, view_ctx,
	    0, 0);
    return ged_draw_scene_ref_erase_path_at_base(gedp, path, base_ref,
	    view_ctx, mode, 1);
}


int
ged_draw_scene_ref_erase_path_prefix_in_active_scope(struct ged *gedp,
						     const char *path,
						     void *view_ctx,
						     int mode)
{
    if (!gedp || !path)
	return 0;

    bsg_scene_ref base_ref = ged_draw_scene_ref_active_scope(gedp, view_ctx,
	    0, 0);
    return ged_draw_scene_ref_erase_path_prefix_at_base(gedp, path, base_ref,
	    view_ctx, mode, 1);
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
	    if (ged_draw_shape_ref_release(gedp, *ref))
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


int
ged_draw_scene_ref_erase_component_name_in_active_scope(struct ged *gedp,
						       const char *name,
						       void *view_ctx,
						       int mode,
						       int nonroot_only)
{
    if (!gedp || !name)
	return 0;

    bsg_scene_ref base_ref = ged_draw_scene_ref_active_scope(gedp, view_ctx,
	    0, 0);
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


static void
color_soltab_scene_ref(struct db_i *dbip, bsg_scene_ref shape_ref)
{
    const struct mater *mp;
    unsigned char basecolor[3] = {0, 0, 0};
    int region_id = ged_draw_scene_ref_legacy_region_id(shape_ref);

    (void)ged_draw_scene_ref_set_legacy_uses_default_color(shape_ref, 0);
    (void)ged_draw_scene_ref_legacy_basecolor(shape_ref, basecolor);

    if (ged_draw_scene_ref_legacy_user_color(shape_ref)) {
	ged_draw_scene_ref_set_material_rgb(shape_ref, basecolor);
	return;
    }

    if (dbip) {
	for (mp = db_mater_head(dbip); mp != MATER_NULL; mp = mp->mt_forw) {
	    if (region_id <= mp->mt_high &&
		    region_id >= mp->mt_low) {
		unsigned char mater_color[3] = {
		    (unsigned char)mp->mt_r,
		    (unsigned char)mp->mt_g,
		    (unsigned char)mp->mt_b
		};
		ged_draw_scene_ref_set_material_rgb(shape_ref, mater_color);
		return;
	    }
	}
    }

    if (ged_draw_scene_ref_legacy_default_color(shape_ref))
	(void)ged_draw_scene_ref_set_legacy_uses_default_color(shape_ref, 1);

    ged_draw_scene_ref_set_material_rgb(shape_ref, basecolor);
}


static int
ged_draw_scene_ref_refresh_material_color(
	struct db_i *dbip,
	bsg_scene_ref shape_ref,
	uint64_t mater_rev)
{
    if (ged_draw_scene_ref_material_revision(shape_ref) == mater_rev)
	return 1;

    color_soltab_scene_ref(dbip, shape_ref);
    return ged_draw_scene_ref_set_material_revision(shape_ref, mater_rev);
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


void *
ged_draw_scene_context_registry_owner(void *scene_ctx)
{
    bsg_scene_ref ref = ged_draw_scene_ref_from_context(scene_ctx);
    if (bsg_scene_ref_is_null(ref))
	return NULL;

    bsg_scene_ref root_ref = ref;
    bsg_scene_ref parent_ref = bsg_scene_parent(ref);
    while (!bsg_scene_ref_is_null(parent_ref)) {
	root_ref = parent_ref;
	parent_ref = bsg_scene_parent(root_ref);
    }

    struct bsg_draw_ctx *ctx = bsg_scene_draw_ctx(root_ref);
    return ctx ? ctx->owner_data : NULL;
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
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!scene_ctx)
	return 0;

    return ged_draw_scene_ref_display_summary(
	    ged_draw_scene_ref_from_context(scene_ctx), out);
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

    struct ged_draw_shape_record_summary summary;
    return ged_draw_scene_ref_shape_record_summary(
	    ged_draw_scene_ref_from_context(shape_ctx), &summary);
}


void *
ged_draw_shape_context_source(void *shape_ctx)
{
    if (!shape_ctx)
	return NULL;

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
    return ged_draw_scene_ref_context(
	    ged_draw_registry_shape_scene_ref(gedp, ref));
}


void *
ged_draw_shape_ref_cache_context(struct ged *gedp, ged_draw_shape_ref ref)
{
    return ged_draw_scene_ref_context(
	    ged_draw_shape_scene_ref_from_cache_ref(gedp, ref));
}


int
ged_draw_shape_ref_record_summary(struct ged *gedp,
				  ged_draw_shape_ref ref,
				  struct ged_draw_shape_record_summary *out)
{
    return ged_draw_scene_ref_shape_record_summary(
	    ged_draw_registry_shape_scene_ref(gedp, ref), out);
}


int
ged_draw_shape_ref_owning_group_record_summary(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	struct ged_draw_group_record_summary *out)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    struct ged_draw_shape_record_summary shape_summary;
    if (!ged_draw_scene_ref_shape_record_summary(shape_ref, &shape_summary))
	return 0;
    return ged_draw_scene_ref_group_record_summary(
	    shape_summary.owning_group_ref, out);
}


int
ged_draw_shape_ref_database_source_summary(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	struct ged_draw_database_source_summary *out)
{
    return ged_draw_scene_ref_database_source_summary(
	    ged_draw_registry_shape_scene_ref(gedp, ref), out);
}


ged_draw_group_ref
ged_draw_group_ref_of_shape(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    struct ged_draw_shape_record_summary shape_summary;
    if (!ged_draw_scene_ref_shape_record_summary(shape_ref, &shape_summary))
	return GED_DRAW_GROUP_REF_NULL;
    return ged_draw_group_ref_from_scene_ref(gedp,
	    shape_summary.owning_group_ref);
}


void *
ged_draw_group_ref_context(struct ged *gedp, ged_draw_group_ref ref)
{
    return ged_draw_scene_ref_context(
	    ged_draw_registry_group_scene_ref(gedp, ref));
}


int
ged_draw_group_ref_record_summary(struct ged *gedp,
				  ged_draw_group_ref ref,
				  struct ged_draw_group_record_summary *out)
{
    return ged_draw_scene_ref_group_record_summary(
	    ged_draw_registry_group_scene_ref(gedp, ref), out);
}


int
ged_draw_group_ref_tree_summary(struct ged *gedp,
				ged_draw_group_ref ref,
				struct ged_draw_scene_tree_summary *out)
{
    return ged_draw_scene_ref_tree_summary(
	    ged_draw_registry_group_scene_ref(gedp, ref), out);
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

    bsg_scene_ref group_ref = ged_draw_registry_group_scene_ref(gedp, ref);
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
    bsg_scene_ref base_ref = ged_draw_scene_ref_active_scope(gedp, view_ctx,
	    1, 1);
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
    if (!ged_draw_scene_ref_is_null(group_ref)) {
	(void)ged_draw_scene_ref_ensure_draw_intent(group_ref, path,
		GED_DRAW_MODE_WIRE);
    }

    if (!ged_draw_scene_ref_is_null(group_ref))
	(void)ged_draw_scene_ref_apply_path_state(gedp, group_ref, dfp);

    bu_free(s, "draw group path string");
    return group_ref;
}


ged_draw_group_ref
ged_draw_group_ref_lookup_or_create(struct ged *gedp,
				    const struct db_full_path *dfp)
{
    bsg_scene_ref group_ref =
	ged_draw_group_ref_lookup_or_create_scene_ref(gedp, dfp);
    if (ged_draw_scene_ref_is_null(group_ref))
	return GED_DRAW_GROUP_REF_NULL;

    return ged_draw_group_ref_from_scene_ref(gedp, group_ref);
}


int
ged_draw_group_ref_set_dbpath(struct ged *gedp,
			      ged_draw_group_ref ref,
			      const struct db_full_path *new_dfp)
{
    if (!gedp || !new_dfp)
	return 0;

    bsg_scene_ref group_ref = ged_draw_registry_group_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(group_ref))
	return 0;

    char *s = db_path_to_string(new_dfp);
    if (!s)
	return 0;

    const char *path = ged_draw_dbpath_skip_lead_slash(s);
    ged_draw_scene_ref_set_name(group_ref, path);

    int mode = GED_DRAW_MODE_WIRE;
    struct ged_draw_group_record_summary group_summary;
    if (ged_draw_scene_ref_group_record_summary(group_ref, &group_summary))
	mode = group_summary.draw_mode;
    (void)ged_draw_scene_ref_set_draw_intent_path(group_ref, path, mode);

    (void)ged_draw_scene_ref_apply_path_state(gedp, group_ref, new_dfp);
    bu_free(s, "ged_draw_group_ref_set_dbpath: path string");
    return 1;
}


int
ged_draw_group_ref_set_mode(struct ged *gedp,
			    ged_draw_group_ref ref,
			    int mode)
{
    bsg_scene_ref group_ref = ged_draw_registry_group_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(group_ref))
	return 0;

    const char *path = NULL;
    struct ged_draw_group_record_summary group_record;
    if (ged_draw_scene_ref_group_record_summary(group_ref, &group_record))
	path = group_record.path;
    if (!path) {
	struct ged_draw_scene_tree_summary group_summary;
	if (ged_draw_scene_ref_tree_summary(group_ref, &group_summary))
	    path = group_summary.name;
    }
    (void)ged_draw_scene_ref_set_draw_intent_mode(group_ref, mode, path);
    return 1;
}


int
ged_draw_group_ref_set_appearance_settings(struct ged *gedp,
					   ged_draw_group_ref ref,
					   const struct ged_draw_appearance_settings *settings)
{
    if (!settings)
	return 0;

    bsg_scene_ref group_ref = ged_draw_registry_group_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(group_ref))
	return 0;

    const char *path = NULL;
    struct ged_draw_scene_tree_summary group_summary;
    if (ged_draw_scene_ref_tree_summary(group_ref, &group_summary))
	path = group_summary.name;
    return ged_draw_scene_ref_set_draw_intent_appearance_settings(group_ref,
	    settings, path);
}


int
ged_draw_group_ref_appearance_settings(struct ged *gedp,
				       ged_draw_group_ref ref,
				       struct ged_draw_appearance_settings *settings)
{
    if (!settings)
	return 0;

    bsg_scene_ref group_ref = ged_draw_registry_group_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(group_ref))
	return 0;

    return ged_draw_scene_ref_draw_intent_appearance_settings(group_ref,
	    settings);
}


static int
ged_draw_group_ref_append_scene_ref(struct ged *gedp,
				    ged_draw_group_ref ref,
				    bsg_scene_ref shape_ref)
{
    /* Appending a single draw command's leaf shapes bumps the draw-scene
     * revision after the first insert.  Keep resolving this private mutation
     * handle by token for the rest of the batch; public record lookups still
     * use revision-stamped refs and reject stale handles. */
    ref.scene_revision = 0;
    bsg_scene_ref group_ref = ged_draw_registry_group_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(group_ref) ||
	    ged_draw_scene_ref_is_null(shape_ref))
	return 0;

    struct ged_draw_shape_record_summary shape_summary;
    int fp_len =
	ged_draw_scene_ref_shape_record_summary(shape_ref, &shape_summary) &&
	shape_summary.fullpath ? (int)shape_summary.fullpath->fp_len : 0;

    if (!gedp || fp_len == 0)
	return ged_draw_scene_ref_append_source_owner_to_group(group_ref,
		shape_ref);

    struct ged_draw_scene_tree_summary group_tree;
    int group_depth =
	ged_draw_scene_ref_tree_summary(group_ref, &group_tree) ?
	group_tree.draw_tree_depth : 0;
    bsg_scene_ref cur_ref = group_ref;
    for (int d = group_depth; d < fp_len - 1; d++) {
	if (!shape_summary.fullpath->fp_names[d])
	    break;
	const char *comp = shape_summary.fullpath->fp_names[d]->d_namep;
	bsg_scene_ref child_ref =
	    ged_draw_group_ref_find_or_create_child_group(gedp, cur_ref, comp);
	if (ged_draw_scene_ref_is_null(child_ref))
	    break;
	cur_ref = child_ref;
    }

    return ged_draw_scene_ref_append_source_owner_to_group(cur_ref,
	    shape_ref);
}


ged_draw_shape_ref
ged_draw_shape_ref_from_context(struct ged *gedp, void *shape_ctx)
{
    return ged_draw_shape_ref_from_scene_ref(gedp,
	    ged_draw_scene_ref_from_context(shape_ctx));
}


const struct db_full_path *
ged_draw_scene_context_fullpath(void *scene_ctx)
{
    if (!scene_ctx)
	return NULL;

    struct ged_draw_scene_tree_summary summary;
    if (!ged_draw_scene_ref_tree_summary(
		ged_draw_scene_ref_from_context(scene_ctx), &summary))
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
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!scene_ctx)
	return 0;

    bsg_scene_ref scene_ref = ged_draw_scene_ref_from_context(scene_ctx);
    return ged_draw_scene_ref_database_source_summary(scene_ref, out);
}


int
ged_draw_scene_context_tree_summary(void *scene_ctx,
				    struct ged_draw_scene_tree_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!scene_ctx)
	return 0;

    bsg_scene_ref scene_ref = ged_draw_scene_ref_from_context(scene_ctx);
    return ged_draw_scene_ref_tree_summary(scene_ref, out);
}


struct ged_draw_context_child_at_ctx {
    size_t index;
    size_t cur;
    void *ctx;
};


static int
_ged_draw_context_child_at_cb(bsg_scene_ref child_ref, void *userdata)
{
    struct ged_draw_context_child_at_ctx *ctx =
	(struct ged_draw_context_child_at_ctx *)userdata;
    if (!ctx)
	return 0;
    if (ctx->cur == ctx->index) {
	ctx->ctx = ged_draw_scene_ref_context(child_ref);
	return 0;
    }
    ctx->cur++;
    return 1;
}


void *
ged_draw_scene_context_child_at(void *scene_ctx, size_t index)
{
    if (!scene_ctx)
	return NULL;

    struct ged_draw_context_child_at_ctx child_ctx;
    child_ctx.index = index;
    child_ctx.cur = 0;
    child_ctx.ctx = NULL;

    (void)ged_draw_scene_ref_foreach_child(
	    ged_draw_scene_ref_from_context(scene_ctx),
	    _ged_draw_context_child_at_cb,
	    &child_ctx);
    return child_ctx.ctx;
}


void *
ged_draw_scene_context_parent(void *scene_ctx)
{
    if (!scene_ctx)
	return NULL;

    return ged_draw_scene_ref_parent_context(
	    ged_draw_scene_ref_from_context(scene_ctx));
}


const char *
ged_draw_scene_context_name(void *scene_ctx)
{
    if (!scene_ctx)
	return NULL;

    struct ged_draw_scene_tree_summary summary;
    if (!ged_draw_scene_ref_tree_summary(
		ged_draw_scene_ref_from_context(scene_ctx), &summary))
	return NULL;

    return summary.name;
}


int
ged_draw_scene_context_subtree_bounds(void *scene_ctx,
				      vect_t *min,
				      vect_t *max,
				      int include_overlays)
{
    if (!scene_ctx || !min || !max)
	return 1;

    bsg_scene_ref scene_ref = ged_draw_scene_ref_from_context(scene_ctx);
    if (ged_draw_scene_ref_is_null(scene_ref))
	return 1;

    return ged_draw_scene_ref_subtree_bounds(scene_ref, min, max, include_overlays);
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

    memset(out, 0, sizeof(*out));
    return ged_draw_scene_ref_material_summary(
	    ged_draw_registry_shape_scene_ref(gedp, ref), out);
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

    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
    if (!shape_data)
	return 0;

    out->fullpath = &shape_data->s_fullpath;
    out->display_name = shape_data->display_name;
    out->owning_group_ref = ged_draw_scene_ref_shape_owning_group(ref);
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


int
ged_draw_scene_context_shape_record_summary(
	void *shape_ctx,
	struct ged_draw_shape_record_summary *out)
{
    return ged_draw_scene_ref_shape_record_summary(
	    ged_draw_scene_ref_from_context(shape_ctx), out);
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

    if (ged_draw_shape_state_get_scene_ref(ref))
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
    return 1;
}


int
ged_draw_scene_context_group_record_summary(
	void *group_ctx,
	struct ged_draw_group_record_summary *out)
{
    return ged_draw_scene_ref_group_record_summary(
	    ged_draw_scene_ref_from_context(group_ctx), out);
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
    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
    struct ged_draw_source_state *source =
	shape_data ? shape_data->source_data : NULL;
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


int
ged_draw_shape_ref_set_center(struct ged *gedp, ged_draw_shape_ref ref,
			      const point_t center)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref) || !center)
	return 0;
    ged_draw_scene_ref_set_draw_center(shape_ref, center);
    return 1;
}


int
ged_draw_shape_ref_geometry_clear(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;
    return ged_draw_scene_ref_geometry_clear(shape_ref);
}


int
ged_draw_shape_ref_update_bounds_from_geometry(struct ged *gedp,
					       ged_draw_shape_ref ref,
					       int *bad_cmd)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
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

    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;

    return ged_draw_scene_ref_source_snapshot(shape_ref, out);
}


int
ged_draw_shape_ref_publish_line_set(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    const point_t *points,
				    const int *commands,
				    size_t point_count)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;
    return ged_draw_scene_ref_publish_line_set(shape_ref, points, commands,
	    point_count);
}


int
ged_draw_shape_ref_publish_indexed_face_set(struct ged *gedp,
					    ged_draw_shape_ref ref,
					    const point_t *points,
					    size_t point_count,
					    const vect_t *normals,
					    size_t normal_count,
					    const int *indices,
					    size_t index_count)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;
    return ged_draw_scene_ref_publish_indexed_face_set(shape_ref, points,
	    point_count, normals, normal_count, indices, index_count);
}


int
ged_draw_shape_ref_publish_primitive_wireframe(struct ged *gedp,
					       ged_draw_shape_ref ref,
					       struct rt_db_internal *ip,
					       const struct bg_tess_tol *ttol,
					       const struct bn_tol *tol,
					       void *view_ctx,
					       int adaptive)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return -1;
    struct rt_view_info view_info;
    const struct rt_view_info *view_infop = NULL;
    if (adaptive && view_ctx) {
	rt_view_context_info_get(&view_info, view_ctx);
	view_infop = &view_info;
    }
    return ged_draw_scene_ref_publish_primitive_wireframe(shape_ref, ip, ttol,
	    tol, view_ctx, view_infop, adaptive);
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
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;

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
	ged_draw_shape_ref_from_scene_ref(ctx->gedp, shape_ref);
    ctx->ret += ged_draw_shape_ref_redraw_wireframe(ctx->gedp, ref,
	    ctx->dbip, ctx->tol, ctx->ttol, ctx->view_ctx,
	    ctx->skip_subtractions);

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
    bsg_scene_ref group_ref = ged_draw_registry_group_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(group_ref))
	return -1;

    struct ged_draw_group_redraw_wireframe_ctx ctx;
    ctx.gedp = gedp;
    ctx.dbip = dbip;
    ctx.tol = tol;
    ctx.ttol = ttol;
    ctx.view_ctx = view_ctx;
    ctx.skip_subtractions = skip_subtractions;
    ctx.ret = 0;

    ged_draw_source_scene_ref_foreach_shape(group_ref,
	    ged_draw_group_ref_redraw_wireframe_shape_cb, &ctx);

    return ctx.ret;
}


int
ged_draw_shape_ref_release(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    return ged_draw_scene_ref_release_source_owner(shape_ref);
}


int
ged_draw_shape_ref_realize_context(struct ged *gedp, ged_draw_shape_ref ref,
				   void *view_ctx)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;
    ged_draw_scene_ref_realize(shape_ref, view_ctx);
    return 1;
}


int
ged_draw_scene_context_set_visible(void *scene_ctx, int visible)
{
    if (!scene_ctx)
	return 0;

    return ged_draw_scene_ref_set_visible(
	    ged_draw_scene_ref_from_context(scene_ctx), visible);
}


int
ged_draw_shape_ref_set_visible(struct ged *gedp, ged_draw_shape_ref ref,
			       int visible)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    return ged_draw_scene_ref_set_visible(shape_ref, visible);
}


int
ged_draw_group_ref_set_visible(struct ged *gedp, ged_draw_group_ref ref,
			       int visible)
{
    bsg_scene_ref group_ref = ged_draw_registry_group_scene_ref(gedp, ref);
    return ged_draw_scene_ref_set_visible(group_ref, visible);
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
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    return ged_draw_scene_ref_set_color(shape_ref, rgb);
}


int
ged_draw_shape_ref_set_highlighted(struct ged *gedp, ged_draw_shape_ref ref,
				   int highlighted)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    return ged_draw_scene_ref_set_highlighted(shape_ref, highlighted);
}


int
ged_draw_shape_ref_set_transparency(struct ged *gedp, ged_draw_shape_ref ref,
				    fastf_t transparency)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    return ged_draw_scene_ref_set_transparency(shape_ref, transparency);
}


int
ged_draw_shape_ref_mark_database_source_changed(struct ged *gedp,
						ged_draw_shape_ref ref,
						ged_draw_stale_reason reason)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    return ged_draw_scene_ref_mark_database_source_changed(shape_ref, reason);
}


int
ged_draw_shape_ref_set_material_color(struct ged *gedp,
				      ged_draw_shape_ref ref,
				      const unsigned char rgb[3])
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref) || !rgb)
	return 0;
    ged_draw_scene_ref_set_material_rgb(shape_ref, rgb);
    return 1;
}


int
ged_draw_shape_ref_refresh_material_color(struct ged *gedp,
					  ged_draw_shape_ref ref,
					  struct db_i *dbip,
					  uint64_t mater_rev)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;
    return ged_draw_scene_ref_refresh_material_color(dbip, shape_ref,
	    mater_rev);
}


int
ged_draw_shape_ref_set_evaluated_region(struct ged *gedp,
					ged_draw_shape_ref ref,
					int evaluated_region)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    return ged_draw_scene_ref_set_legacy_eval_flag(shape_ref,
	    evaluated_region);
}


int
ged_draw_shape_ref_lod_ensure(struct ged *gedp, ged_draw_shape_ref ref,
			      void *first_view_ctx,
			      void **view_ctxs,
			      size_t view_ctx_count)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref) || !first_view_ctx)
	return 0;

    ged_draw_view_lod_policy policy;
    rt_view_context_lod_policy_get(&policy, first_view_ctx);
    struct ged_draw_scene_draw_state_summary draw_state = {0};
    if (!ged_draw_scene_ref_draw_state_summary(shape_ref, &draw_state))
	return 1;
    int candidate = draw_state.pipeline_candidate;
    struct ged_draw_database_source_summary source_summary;
    int source_backed =
	ged_draw_scene_ref_database_source_summary(shape_ref, &source_summary) &&
	source_summary.has_state;
    if (!candidate && !source_backed)
	return 1;
    if (!candidate) {
	struct ged_draw_shape_source_snapshot source = {0};
	struct directory *dp = NULL;
	if (ged_draw_scene_ref_source_snapshot(shape_ref, &source))
	    dp = source.leaf_dp;
	int mode = draw_state.draw_mode;
	int mesh_ref = dp && policy.mesh_enabled &&
	    ((dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT &&
	      (mode == 0 || mode == 1)) ||
	     (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP &&
	      mode == 1));
	if (!mesh_ref && !(policy.csg_enabled && mode == 0))
	    return 1;
    }
    int realized = 0;
    for (size_t i = 0; i < view_ctx_count; i++) {
	if (view_ctxs && view_ctxs[i]) {
	    (void)ged_draw_scene_ref_realize_view_inputs_changed(shape_ref,
		    view_ctxs[i]);
	    realized = 1;
	}
    }
    if (!realized) {
	(void)ged_draw_scene_ref_realize_view_inputs_changed(shape_ref,
		first_view_ctx);
    }
    return 1;
}


void *
ged_draw_shape_ref_view_context(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    struct ged_draw_scene_draw_state_summary draw_state = {0};
    if (!ged_draw_scene_ref_draw_state_summary(shape_ref, &draw_state))
	return NULL;
    return draw_state.view_ctx;
}


static int
ged_draw_scene_ref_source_runtime_summary(
	bsg_scene_ref ref,
	struct ged_draw_source_runtime_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
    struct ged_draw_source_state *source =
	shape_data ? shape_data->source_data : NULL;
    if (!source)
	return 0;

    out->dbip = source->dbip;
    out->fullpath = ged_draw_scene_ref_fullpath(ref);
    out->leaf_dp = ged_draw_scene_ref_leaf_dp(ref);
    out->name = ged_draw_scene_ref_name(ref);
    out->tol = source->tol;
    out->ttol = source->ttol;
    out->rt_mesh_lod = source->rt_mesh_lod;
    out->mesh_lod_bounds_valid = source->mesh_lod_bounds_valid;
    VMOVE(out->mesh_lod_bmin, source->mesh_lod_bmin);
    VMOVE(out->mesh_lod_bmax, source->mesh_lod_bmax);
    return 1;
}


static int
ged_draw_scene_ref_source_mesh_lod_set(
	bsg_scene_ref ref,
	struct rt_mesh_lod *lod)
{
    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
    struct ged_draw_source_state *source =
	shape_data ? shape_data->source_data : NULL;
    if (!source)
	return 0;

    source->rt_mesh_lod = lod;
    return 1;
}


static int
ged_draw_scene_ref_source_mesh_lod_bounds_set(
	bsg_scene_ref ref,
	const point_t bmin,
	const point_t bmax)
{
    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
    struct ged_draw_source_state *source =
	shape_data ? shape_data->source_data : NULL;
    if (!source)
	return 0;

    VMOVE(source->mesh_lod_bmin, bmin);
    VMOVE(source->mesh_lod_bmax, bmax);
    source->mesh_lod_bounds_valid = 1;
    return 1;
}


static int
_ged_draw_scene_ref_set_transformed_bounds(
	bsg_scene_ref ref,
	const point_t src_bmin,
	const point_t src_bmax)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 0;

    struct ged_draw_scene_draw_state_summary draw_state;
    memset(&draw_state, 0, sizeof(draw_state));
    if (!ged_draw_scene_ref_draw_state_summary(ref, &draw_state) ||
	    !draw_state.draw_mat_valid)
	return 0;

    point_t corners[8] = {
	{src_bmin[X], src_bmin[Y], src_bmin[Z]},
	{src_bmin[X], src_bmin[Y], src_bmax[Z]},
	{src_bmin[X], src_bmax[Y], src_bmin[Z]},
	{src_bmin[X], src_bmax[Y], src_bmax[Z]},
	{src_bmax[X], src_bmin[Y], src_bmin[Z]},
	{src_bmax[X], src_bmin[Y], src_bmax[Z]},
	{src_bmax[X], src_bmax[Y], src_bmin[Z]},
	{src_bmax[X], src_bmax[Y], src_bmax[Z]}
    };

    point_t bmin, bmax, tp;
    MAT4X3PNT(tp, draw_state.draw_mat, corners[0]);
    VMOVE(bmin, tp);
    VMOVE(bmax, tp);
    for (int i = 1; i < 8; i++) {
	MAT4X3PNT(tp, draw_state.draw_mat, corners[i]);
	VMIN(bmin, tp);
	VMAX(bmax, tp);
    }

    return ged_draw_scene_ref_set_bounds(ref, bmin, bmax, 1);
}


static int
ged_draw_scene_ref_source_mesh_lod_bounds_restore(bsg_scene_ref ref)
{
    struct ged_draw_source_runtime_summary source;
    memset(&source, 0, sizeof(source));
    if (ged_draw_scene_ref_is_null(ref) ||
	    !ged_draw_scene_ref_source_runtime_summary(ref, &source) ||
	    !source.mesh_lod_bounds_valid)
	return 0;

    return _ged_draw_scene_ref_set_transformed_bounds(ref,
	    source.mesh_lod_bmin, source.mesh_lod_bmax);
}


static int
ged_draw_scene_ref_source_mesh_lod_load_view(
	bsg_scene_ref ref,
	void *view_ctx)
{
    struct ged_draw_source_runtime_summary source;
    memset(&source, 0, sizeof(source));
    if (ged_draw_scene_ref_is_null(ref) ||
	    !ged_draw_scene_ref_source_runtime_summary(ref, &source) ||
	    !source.rt_mesh_lod)
	return -1;

    return rt_mesh_lod_load_view_scene_ref(source.rt_mesh_lod,
	    ged_draw_scene_ref_to_rt_view_ref(ref), view_ctx, 0);
}


static int
ged_draw_scene_ref_source_bot_mesh_lod_realize(
	bsg_scene_ref ref,
	void *view_ctx)
{
    if (ged_draw_scene_ref_is_null(ref) || !view_ctx)
	return 0;

    ged_draw_scene_ref_realization_prepare_mesh(ref);

    struct ged_draw_source_runtime_summary source;
    memset(&source, 0, sizeof(source));
    if (!ged_draw_scene_ref_source_runtime_summary(ref, &source))
	return 0;
    ged_draw_log(1, "bot_lod_mesh_realize %s[%s]", source.name,
	    _ged_draw_source_view_context_name(view_ctx));

    if (!source.rt_mesh_lod) {
	struct db_i *dbip = source.dbip;
	struct directory *dp = source.leaf_dp;

	if (!dbip || !dp)
	    return 0;

	struct rt_mesh_lod_cache_status status = RT_MESH_LOD_CACHE_STATUS_INIT;
	if (db_mesh_lod_status(dbip, dp->d_namep, &status) != BRLCAD_OK)
	    return 0;
	if (!status.has_cache_key || !status.has_cached_payload || status.stale_cache_entry) {
	    if (db_mesh_lod_refresh(dbip, dp->d_namep, &status) != BRLCAD_OK)
		return 0;
	}
	if (!status.has_cache_key || !status.has_cached_payload)
	    return 0;

	struct rt_mesh_lod *rt_lod = db_mesh_lod_get(dbip, dp->d_namep);
	if (!rt_lod) {
	    if (db_mesh_lod_refresh(dbip, dp->d_namep, &status) != BRLCAD_OK ||
		    !status.has_cache_key || !status.has_cached_payload)
		return 0;
	    rt_lod = db_mesh_lod_get(dbip, dp->d_namep);
	}
	if (!rt_lod)
	    return 0;

	ged_draw_scene_ref_source_mesh_lod_set(ref, rt_lod);

	int level = ged_draw_scene_ref_source_mesh_lod_load_view(ref, view_ctx);
	if (level < 0)
	    bu_log("Error loading info for initial LoD view\n");

	struct rt_mesh_lod_info info = RT_MESH_LOD_INFO_INIT;
	if (rt_mesh_lod_info_get(rt_lod, &info))
	    ged_draw_scene_ref_source_mesh_lod_bounds_set(ref, info.bmin,
		    info.bmax);

	(void)ged_draw_scene_ref_source_mesh_lod_bounds_restore(ref);
    }

    (void)ged_draw_scene_ref_source_mesh_lod_load_view(ref, view_ctx);
    return ged_draw_scene_ref_publish_current_mesh_lod(ref);
}


static int
ged_draw_scene_ref_source_brep_mesh_lod_realize(
	bsg_scene_ref ref,
	void *view_ctx)
{
    if (ged_draw_scene_ref_is_null(ref) || !view_ctx)
	return 0;

    ged_draw_scene_ref_realization_prepare_mesh(ref);

    struct ged_draw_source_runtime_summary source;
    memset(&source, 0, sizeof(source));
    if (!ged_draw_scene_ref_source_runtime_summary(ref, &source))
	return 0;
    ged_draw_log(1, "brep_lod_mesh_realize %s[%s]", source.name,
	    _ged_draw_source_view_context_name(view_ctx));

    if (!source.rt_mesh_lod) {
	struct db_i *dbip = source.dbip;
	struct directory *dp = source.leaf_dp;

	if (!dbip || !dp)
	    return 0;

	const struct bn_tol *tol = source.tol;
	const struct bg_tess_tol *ttol = source.ttol;
	struct rt_mesh_lod *rt_lod = NULL;
	point_t bmin, bmax;
	int bounds_valid = 0;
	if (!ged_draw_brep_mesh_lod_cache_prepare(&rt_lod, bmin, bmax,
		&bounds_valid, dbip, dp, ttol, tol))
	    return 0;
	if (!rt_lod)
	    return 0;

	ged_draw_scene_ref_source_mesh_lod_set(ref, rt_lod);
	if (bounds_valid)
	    ged_draw_scene_ref_source_mesh_lod_bounds_set(ref, bmin, bmax);

	(void)ged_draw_scene_ref_source_mesh_lod_bounds_restore(ref);

	if (!ged_draw_brep_mesh_lod_detail_setup(rt_lod, dbip, dp, ttol, tol))
	    return 0;

	int level = ged_draw_scene_ref_source_mesh_lod_load_view(ref, view_ctx);
	if (level < 0)
	    bu_log("Error loading info for initial LoD view\n");
    }

    (void)ged_draw_scene_ref_source_mesh_lod_load_view(ref, view_ctx);
    return ged_draw_scene_ref_publish_current_mesh_lod(ref);
}


static int
ged_draw_scene_ref_source_adaptive_wireframe_update(
	bsg_scene_ref ref,
	void *view_ctx,
	int force)
{
    if (ged_draw_scene_ref_is_null(ref) || !view_ctx)
	return 0;

    ged_draw_view_lod_policy policy;
    ged_draw_view_context_lod_policy_get(&policy, view_ctx);
    if (!policy.csg_enabled)
	return 0;

    struct ged_draw_scene_draw_state_summary draw_state;
    memset(&draw_state, 0, sizeof(draw_state));
    if (!ged_draw_scene_ref_draw_state_summary(ref, &draw_state))
	return 0;

    ged_draw_scene_ref_realization_prepare_wireframe(ref);

    if (!draw_state.bounds_valid)
	return 0;
    if (!(ged_view_context_perspective_get(view_ctx) > SMALL_FASTF)) {
	point_t obb_center = VINIT_ZERO;
	vect_t obb_extent1 = VINIT_ZERO;
	vect_t obb_extent2 = VINIT_ZERO;
	vect_t obb_extent3 = VINIT_ZERO;
	if (!rt_view_context_obb_get(view_ctx, obb_center, obb_extent1,
		obb_extent2, obb_extent3) ||
		!bg_sat_aabb_obb(draw_state.bounds_min,
		    draw_state.bounds_max, obb_center, obb_extent1,
		    obb_extent2, obb_extent3))
	    return 0;
    }

    int rework = force ? 1 : 0;
    if (!rework && !NEAR_EQUAL(ged_draw_scene_ref_realization_curve_scale(ref),
	    policy.curve_scale, SMALL_FASTF))
	rework = 1;
    if (!rework && !NEAR_EQUAL(ged_draw_scene_ref_realization_point_scale(ref),
	    policy.point_scale, SMALL_FASTF))
	rework = 1;
    if (!rework) {
	fastf_t view_scale = ged_draw_scene_ref_realization_view_scale(ref);
	fastf_t current_view_scale = ged_view_context_scale_get(view_ctx);
	fastf_t delta = view_scale * 0.1/view_scale;
	if (!NEAR_EQUAL(view_scale, current_view_scale, delta))
	    rework = 1;
    }
    if (!rework)
	return 0;

    struct ged_draw_shape_source_snapshot source;
    memset(&source, 0, sizeof(source));
    if (!ged_draw_scene_ref_source_snapshot(ref, &source))
	return 0;
    (void)ged_draw_scene_ref_realization_prepare_view_redraw(ref, view_ctx);

    struct directory *dp = source.leaf_dp;
    struct db_i *dbip = source.dbip;
    struct rt_db_internal dbintern;
    RT_DB_INTERNAL_INIT(&dbintern);
    struct rt_db_internal *ip = &dbintern;
    if (!draw_state.draw_mat_valid)
	return 0;
    int ret = rt_db_get_internal(ip, dp, dbip, draw_state.draw_mat);
    if (ret < 0)
	return 0;

    struct rt_view_info view_info;
    ged_draw_view_context_info_get(&view_info, view_ctx);
    (void)ged_draw_scene_ref_publish_current_wireframe(ref, ip, NULL,
	    source.tol, view_ctx, &view_info, 1);

    rt_db_free_internal(ip);

    return 1;
}


static int
ged_draw_scene_ref_source_wireframe_realize(
	bsg_scene_ref ref,
	void *view_ctx,
	struct rt_db_internal *ip)
{
    struct ged_draw_shape_source_snapshot source;
    memset(&source, 0, sizeof(source));
    if (!ged_draw_scene_ref_source_snapshot(ref, &source))
	return 0;

    ged_draw_scene_ref_realization_prepare_wireframe(ref);

    ged_draw_view_lod_policy policy;
    ged_draw_view_context_lod_policy_get(&policy, view_ctx);

    if (!view_ctx || !policy.csg_enabled) {
	struct ged_draw_scene_draw_state_summary draw_state;
	memset(&draw_state, 0, sizeof(draw_state));
	if (!ged_draw_scene_ref_draw_state_summary(ref, &draw_state))
	    return 0;
	return ged_draw_scene_ref_publish_current_wireframe(ref, ip, source.ttol,
		source.tol, draw_state.view_ctx, NULL, 0);
    }

    if (ged_draw_source_primitive_has_lod_realize(ip))
	return ged_draw_scene_ref_source_adaptive_wireframe_update(ref, view_ctx,
		1);

    struct ged_draw_scene_draw_state_summary draw_state;
    memset(&draw_state, 0, sizeof(draw_state));
    if (!ged_draw_scene_ref_draw_state_summary(ref, &draw_state))
	return 0;
    return ged_draw_scene_ref_publish_current_wireframe(ref, ip, source.ttol,
	    source.tol, draw_state.view_ctx, NULL, 0);
}


static int
ged_draw_scene_ref_source_face_set_realize(
	bsg_scene_ref ref,
	void *view_ctx,
	struct rt_db_internal *ip,
	const char *name)
{
    if (!ip || !ip->idb_meth || !ip->idb_meth->ft_indexed_face_set)
	return 0;

    struct ged_draw_shape_source_snapshot source;
    memset(&source, 0, sizeof(source));
    if (!ged_draw_scene_ref_source_snapshot(ref, &source))
	return 1;

    struct rt_view_info view_info;
    ged_draw_view_context_info_get(&view_info, view_ctx);
    if (!ged_draw_scene_ref_publish_current_face_set(ref, ip, source.ttol,
	    source.tol, &view_info)) {
	bu_log("ERROR(%s): %s shaded face-set publication failed\n",
		name ? name : "<unknown>", ip->idb_meth->ft_label);
    }

    return 1;
}


static int
_ged_draw_scene_ref_source_tessellation_fallback(
	bsg_scene_ref ref,
	void *view_ctx,
	struct rt_db_internal *ip,
	const char *name,
	const char *mode_name)
{
    if (ged_draw_scene_ref_strict_fallback(ref)) {
	bu_log("ERROR(%s): %s tessellation failed; geometry cleared\n",
		name ? name : "<unknown>", mode_name ? mode_name : "draw");
	ged_draw_scene_ref_realization_clear_stale(ref);
	return 0;
    }

    bu_log("WARNING(%s): %s tessellation failed; falling back to wireframe\n",
	    name ? name : "<unknown>", mode_name ? mode_name : "draw");
    return ged_draw_scene_ref_source_wireframe_realize(ref, view_ctx, ip);
}


static int
ged_draw_scene_ref_source_tessellate_nmg_realize(
	bsg_scene_ref ref,
	void *view_ctx,
	struct rt_db_internal *ip,
	const char *name,
	const char *mode_name)
{
    struct ged_draw_shape_source_snapshot source;
    memset(&source, 0, sizeof(source));
    if (!ged_draw_scene_ref_source_snapshot(ref, &source))
	return 0;

    struct directory *dp = source.leaf_dp;
    const struct bn_tol *tol = source.tol;
    const struct bg_tess_tol *ttol = source.ttol;
    RT_CK_DB_INTERNAL(ip);
    RT_CK_DIR(dp);
    BN_CK_TOL(tol);
    BG_CK_TESS_TOL(ttol);
    if (!ip->idb_meth || !ip->idb_meth->ft_tessellate) {
	bu_log("ERROR(%s): tessellation support not available\n", dp->d_namep);
	return _ged_draw_scene_ref_source_tessellation_fallback(ref, view_ctx,
		ip, name ? name : dp->d_namep, mode_name);
    }

    struct model *m = nmg_mm();
    struct nmgregion *r = (struct nmgregion *)NULL;
    if (ip->idb_meth->ft_tessellate(&r, m, ip, ttol, tol) < 0) {
	bu_log("ERROR(%s): tessellation failure\n", dp->d_namep);
	nmg_km(m);
	return _ged_draw_scene_ref_source_tessellation_fallback(ref, view_ctx,
		ip, name ? name : dp->d_namep, mode_name);
    }

    NMG_CK_REGION(r);
    if (!ged_draw_scene_ref_publish_current_nmg_region(ref, r,
	    GED_DRAW_NMG_STYLE_POLYGON)) {
	nmg_km(m);
	return _ged_draw_scene_ref_source_tessellation_fallback(ref, view_ctx,
		ip, name ? name : dp->d_namep, mode_name);
    }
    nmg_km(m);

    return 1;
}


static int
ged_draw_scene_ref_source_standard_realize(
	bsg_scene_ref ref,
	void *view_ctx)
{
    struct ged_draw_scene_draw_state_summary draw_state;
    memset(&draw_state, 0, sizeof(draw_state));
    if (!ged_draw_scene_ref_draw_state_summary(ref, &draw_state))
	return 0;

    struct ged_draw_source_runtime_summary source;
    memset(&source, 0, sizeof(source));
    if (!ged_draw_scene_ref_source_runtime_summary(ref, &source))
	return 0;

    struct db_i *dbip = source.dbip;
    const struct db_full_path *fp = source.fullpath;
    if (fp && fp->fp_len <= 0)
	return 0;
    struct directory *dp = source.leaf_dp;
    if (!dp)
	return 0;

    struct rt_db_internal dbintern;
    RT_DB_INTERNAL_INIT(&dbintern);
    struct rt_db_internal *ip = &dbintern;
    if (!draw_state.draw_mat_valid)
	return 0;
    int ret = rt_db_get_internal(ip, dp, dbip, draw_state.draw_mat);
    if (ret < 0)
	return 0;

    if (ip->idb_major_type != DB5_MAJORTYPE_BRLCAD) {
	(void)ged_draw_scene_ref_source_wireframe_realize(ref, view_ctx, ip);
	goto geom_done;
    }

    if (ip->idb_minor_type == DB5_MINORTYPE_BRLCAD_PIPE) {
	(void)ged_draw_scene_ref_source_wireframe_realize(ref, view_ctx, ip);
	goto geom_done;
    }

    if (draw_state.draw_mode > 0 &&
	    ged_draw_scene_ref_source_face_set_realize(ref, view_ctx, ip,
		dp->d_namep)) {
	goto geom_done;
    }

    switch (draw_state.draw_mode) {
	case GED_DRAW_MODE_WIRE:
	case GED_DRAW_MODE_SHADED_BOTS:
	    (void)ged_draw_scene_ref_source_wireframe_realize(ref, view_ctx, ip);
	    ged_draw_scene_ref_set_draw_mode(ref, GED_DRAW_MODE_WIRE);
	    break;
	case GED_DRAW_MODE_SHADED:
	    (void)ged_draw_scene_ref_source_tessellate_nmg_realize(ref,
		    view_ctx, ip, dp->d_namep, "shaded");
	    break;
	case GED_DRAW_MODE_EVAL_WIRE:
	    bu_log("Error - got too deep into _scene_obj_draw routine with evaluated wireframe mode\n");
	    goto cleanup;
	case GED_DRAW_MODE_HIDDEN_LINE:
	    (void)ged_draw_scene_ref_source_tessellate_nmg_realize(ref,
		    view_ctx, ip, dp->d_namep, "hidden-line");
	    break;
	case GED_DRAW_MODE_EVAL_POINTS:
	    bu_log("Error - got too deep into _scene_obj_draw routine with evaluated sampled-points mode\n");
	    goto cleanup;
	default:
	    (void)ged_draw_scene_ref_source_wireframe_realize(ref, view_ctx, ip);
	    break;
    }

geom_done:
    ged_draw_scene_ref_realization_finish_view(ref, view_ctx,
	    draw_state.view_ctx);

cleanup:
    rt_db_free_internal(&dbintern);
    return 1;
}


static int
ged_draw_scene_ref_source_realize(
	bsg_scene_ref ref,
	void *view_ctx)
{
    if (!ged_draw_scene_ref_is_shape(ref))
	return 0;

    struct ged_draw_scene_draw_state_summary draw_state;
    memset(&draw_state, 0, sizeof(draw_state));
    if (!ged_draw_scene_ref_draw_state_summary(ref, &draw_state))
	return 0;

    struct ged_draw_source_runtime_summary source;
    memset(&source, 0, sizeof(source));
    if (!ged_draw_scene_ref_source_runtime_summary(ref, &source))
	return 0;

    const struct db_full_path *fp = source.fullpath;
    if (fp && fp->fp_len <= 0)
	return 0;

    struct directory *dp = source.leaf_dp;
    if (!dp)
	return 0;

    ged_draw_view_lod_policy policy;
    ged_draw_view_context_lod_policy_get(&policy, view_ctx);

    if (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT && view_ctx &&
	    policy.mesh_enabled &&
	    (draw_state.draw_mode == GED_DRAW_MODE_WIRE ||
	     draw_state.draw_mode == GED_DRAW_MODE_SHADED_BOTS))
	return ged_draw_scene_ref_source_bot_mesh_lod_realize(ref, view_ctx);

    if (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP && view_ctx &&
	    policy.mesh_enabled &&
	    draw_state.draw_mode == GED_DRAW_MODE_SHADED_BOTS)
	return ged_draw_scene_ref_source_brep_mesh_lod_realize(ref, view_ctx);

    return ged_draw_scene_ref_source_standard_realize(ref, view_ctx);
}


struct ged_draw_realize_children_ctx {
    void *view_ctx;
};


static int
_ged_draw_realize_child_cb(bsg_scene_ref child_ref, void *userdata)
{
    struct ged_draw_realize_children_ctx *ctx =
	(struct ged_draw_realize_children_ctx *)userdata;

    ged_draw_scene_ref_realize(child_ref, ctx ? ctx->view_ctx : NULL);
    return 1;
}


static int
ged_draw_scene_ref_container_realize(
	bsg_scene_ref ref,
	void *view_ctx)
{
    struct ged_draw_source_runtime_summary source;
    memset(&source, 0, sizeof(source));
    if (ged_draw_scene_ref_is_shape(ref) &&
	    ged_draw_scene_ref_source_runtime_summary(ref, &source))
	return 0;

    struct ged_draw_realize_children_ctx ctx = { view_ctx };
    (void)ged_draw_scene_ref_foreach_child(ref, _ged_draw_realize_child_cb,
	    &ctx);
    return 1;
}


static int
ged_draw_scene_ref_evaluated_realize(
	bsg_scene_ref ref,
	void *view_ctx)
{
    struct ged_draw_scene_draw_state_summary draw_state;
    memset(&draw_state, 0, sizeof(draw_state));
    if (!ged_draw_scene_ref_draw_state_summary(ref, &draw_state))
	return 0;

    if (draw_state.draw_mode != GED_DRAW_MODE_EVAL_WIRE &&
	    draw_state.draw_mode != GED_DRAW_MODE_EVAL_POINTS)
	return 0;

    struct ged *gedp = _ged_draw_scene_ref_owner_gedp(ref);
    ged_draw_shape_ref shape_ref = gedp ?
	ged_draw_shape_ref_from_scene_ref(gedp, ref) : GED_DRAW_SHAPE_REF_NULL;

    if (draw_state.draw_mode == GED_DRAW_MODE_EVAL_WIRE)
	ged_draw_shape_ref_eval_wireframe(gedp, shape_ref);
    else
	ged_draw_shape_ref_eval_points(gedp, shape_ref);

    ged_draw_scene_ref_realization_finish_current_view(ref, view_ctx, NULL);
    return 1;
}


static void
ged_draw_scene_ref_realize_dispatch(
	bsg_scene_ref ref,
	void *view_ctx)
{
    if (ged_draw_scene_ref_realization_current(ref) && !view_ctx)
	return;

    struct ged_draw_scene_draw_state_summary draw_state;
    memset(&draw_state, 0, sizeof(draw_state));
    if (!ged_draw_scene_ref_draw_state_summary(ref, &draw_state))
	return;
    ged_draw_log(1, "draw_scene %s[%s]", draw_state.name,
	    _ged_draw_source_view_context_name(view_ctx));

    ged_draw_view_lod_policy policy;
    ged_draw_view_context_lod_policy_get(&policy, view_ctx);
    if (view_ctx && !policy.csg_enabled && !policy.mesh_enabled) {
	ged_draw_scene_ref_realize(ref, NULL);
	return;
    }

    if (ged_draw_scene_ref_container_realize(ref, view_ctx))
	return;

    if (ged_draw_scene_ref_evaluated_realize(ref, view_ctx))
	return;

    (void)ged_draw_scene_ref_source_realize(ref, view_ctx);
}


static int
_ged_draw_scene_ref_overlay_internal_set(
	struct ged *gedp,
	bsg_scene_ref ref,
	struct db_full_path *fp,
	struct rt_db_internal **ip)
{
    if (ged_draw_scene_ref_is_null(ref) || !ip || !*ip)
	return 0;

    if (!ged_draw_scene_ref_apply_path_state(gedp, ref, fp))
	return 0;

    if (!_ged_draw_scene_ref_user_data_set(ref, (void *)*ip,
	    GED_DRAW_SHAPE_USER_DATA_RT_DB_INTERNAL))
	return 0;

    *ip = NULL;
    return 1;
}


static void
ged_draw_scene_ref_realize(bsg_scene_ref ref, void *view_ctx)
{
    ged_draw_scene_ref_realize_dispatch(ref, view_ctx);
}


static bsg_scene_ref
ged_draw_view_context_overlay_scene_find(void *view_ctx, const char *name)
{
    if (!view_ctx || !name)
	return bsg_scene_ref_null();

    return bsg_overlay_find_scene((struct bsg_view *)view_ctx, name);
}


static void
ged_draw_view_context_overlay_name_erase(void *view_ctx, const char *name)
{
    if (!view_ctx || !name)
	return;

    bsg_overlay_erase_name((struct bsg_view *)view_ctx, name);
}


static int
ged_draw_view_context_overlay_scene_append(void *view_ctx, bsg_scene_ref scene)
{
    if (!view_ctx || ged_draw_scene_ref_is_null(scene))
	return 0;

    return bsg_overlay_append_scene((struct bsg_view *)view_ctx, scene);
}


static int
ged_draw_view_overlay_command_result_owner_set(bsg_scene_ref scene,
					       const void *owner,
					       const char *source_path)
{
    if (ged_draw_scene_ref_is_null(scene))
	return 0;

    return bsg_overlay_register_scene_owner(scene, owner,
	    BSG_OVERLAY_ROLE_MODEL,
	    BSG_OVERLAY_CLASS_COMMAND_RESULT,
	    BSG_OVERLAY_LC_PERSISTENT,
	    BSG_OVERLAY_ORDER_MODEL,
	    source_path, 0);
}


static bsg_scene_ref
ged_draw_view_context_overlay_geometry_create(
	void *view_ctx,
	const char *name,
	enum ged_draw_overlay_geometry_kind kind)
{
    struct bsg_view *view = (struct bsg_view *)view_ctx;
    if (!view || !name)
	return bsg_scene_ref_null();

    switch (kind) {
	case GED_DRAW_OVERLAY_GEOMETRY_LINE_SET:
	    return bsg_geometry_ref_as_scene(
		    bsg_line_set_ref_as_geometry(bsg_line_set_ref_create(view, name)));
	case GED_DRAW_OVERLAY_GEOMETRY_POINT_SET:
	    return bsg_geometry_ref_as_scene(
		    bsg_point_set_ref_as_geometry(bsg_point_set_ref_create(view, name)));
	case GED_DRAW_OVERLAY_GEOMETRY_INDEXED_FACE_SET:
	    return bsg_geometry_ref_as_scene(
		    bsg_indexed_face_set_ref_as_geometry(
			bsg_indexed_face_set_ref_create(view, name)));
	default:
	    return bsg_scene_ref_null();
    }
}


static bsg_scene_ref
ged_draw_view_context_overlay_create(void *view_ctx, const char *name)
{
    if (!view_ctx || !name)
	return bsg_scene_ref_null();

    bsg_feature_ref ref = bsg_feature_create_overlay(
	    (struct bsg_view *)view_ctx, name, 0);
    if (bsg_feature_ref_is_null(ref))
	return bsg_scene_ref_null();

    return bsg_feature_ref_as_scene(ref);
}


static int
_ged_draw_view_context_overlay_internal_create(
	struct ged *gedp,
	void *view_ctx,
	const char *name,
	struct db_full_path *fp,
	struct rt_db_internal **ip,
	bsg_scene_ref *out)
{
    if (out)
	*out = ged_draw_scene_ref_null();

    bsg_scene_ref ref = ged_draw_view_context_overlay_create(view_ctx, name);
    if (ged_draw_scene_ref_is_null(ref))
	return 0;

    if (!_ged_draw_scene_ref_overlay_internal_set(gedp, ref, fp, ip)) {
	ged_draw_scene_ref_release(ref);
	return 0;
    }

    if (out)
	*out = ref;
    return 1;
}


int
ged_draw_view_context_overlay_internal_create_context(
	struct ged *gedp,
	void *view_ctx,
	const char *name,
	struct db_full_path *fp,
	struct rt_db_internal **ip,
	void **out_ctx)
{
    if (out_ctx)
	*out_ctx = NULL;

    bsg_scene_ref ref = ged_draw_scene_ref_null();
    if (!_ged_draw_view_context_overlay_internal_create(gedp, view_ctx, name,
	    fp, ip, &ref))
	return 0;

    if (out_ctx)
	*out_ctx = ged_draw_scene_ref_context(ref);
    return 1;
}


void
ged_draw_overlay_erase_name_context(struct ged *gedp, void *view_ctx,
				    const char *name)
{
    if (!view_ctx || !name)
	return;

    bsg_scene_ref ref = ged_draw_view_context_overlay_scene_find(view_ctx, name);
    if (!ged_draw_scene_ref_is_null(ref))
	ged_draw_scene_context_highlight_free_cb(ged_draw_scene_ref_context(ref));
    ged_draw_view_context_overlay_name_erase(view_ctx, name);
    (void)gedp;
}


static int
_ged_overlay_publish_geometry(bsg_scene_ref ref,
			      const struct ged_draw_overlay_geometry *geometry)
{
    if (ged_draw_scene_ref_is_null(ref) || !geometry)
	return 0;

    switch (geometry->kind) {
	case GED_DRAW_OVERLAY_GEOMETRY_LINE_SET:
	    if (geometry->point_count != geometry->command_count)
		return 0;
	    return ged_draw_scene_ref_publish_line_set(ref,
		    geometry->points, geometry->commands, geometry->point_count);
	case GED_DRAW_OVERLAY_GEOMETRY_POINT_SET:
	    return ged_draw_scene_ref_publish_point_set(ref,
		    geometry->points, geometry->point_count);
	case GED_DRAW_OVERLAY_GEOMETRY_INDEXED_FACE_SET:
	    return ged_draw_scene_ref_publish_indexed_face_set(ref,
		    geometry->points, geometry->point_count,
		    geometry->normals, geometry->normal_count,
		    geometry->indices, geometry->index_count);
	default:
	    return 0;
    }
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

    bsg_scene_ref overlay_scene =
	ged_draw_view_context_overlay_geometry_create(view_ctx, name,
		geometry->kind);
    if (ged_draw_scene_ref_is_null(overlay_scene))
	return -1;

    (void)ged_draw_view_overlay_command_result_owner_set(overlay_scene,
	    gedp, name);
    (void)ged_draw_scene_ref_set_name(overlay_scene, name);

    if (!ged_draw_scene_context_ensure_registry_entry(gedp,
	    ged_draw_scene_ref_context(overlay_scene))) {
	ged_draw_scene_ref_release(overlay_scene);
	return -1;
    }

    if (!_ged_overlay_publish_geometry(overlay_scene, geometry)) {
	ged_draw_scene_ref_release(overlay_scene);
	return -1;
    }
    ged_draw_scene_ref_update_bounds_context(overlay_scene, view_ctx);

    if (!ged_draw_view_context_overlay_scene_append(view_ctx, overlay_scene)) {
	ged_draw_scene_ref_release(overlay_scene);
	return -1;
    }

    (void)ged_draw_scene_ref_apply_overlay_geometry_attributes(overlay_scene,
	    basecolor, transparency, draw_mode);

    if (csoltab && gedp->dbip)
	color_soltab_scene_ref(gedp->dbip, overlay_scene);

    if (out)
	*out = ged_draw_shape_ref_from_scene_ref(gedp, overlay_scene);

    return 0;
}


static int
ged_draw_scene_ref_commit_database_leaf_draft(
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
	fastf_t draw_size)
{
    if (ged_draw_scene_ref_is_null(parent_ref) || !gedp || !dbip || !path ||
	    !draw_mat || !rgb)
	return 0;

    ged_draw_shape_draft *draft = ged_draw_shape_draft_create_context(gedp,
	    view_ctx, 1);
    if (!draft)
	return 0;

    char *name = db_path_to_string(path);
    ged_draw_shape_draft_apply_path_source_state(draft, dbip, path, tol,
	    ttol, 1, draw_mat, name);
    bu_free(name, "path string");
    ged_draw_shape_draft_apply_database_leaf_display(draft, settings,
	    bool_op, rgb, has_draw_size, draw_size);

    return ged_draw_shape_draft_commit_database_leaf_to_scene_ref(draft,
	    parent_ref);
}


int
ged_draw_scene_context_commit_database_leaf_draft(
	void *parent_ctx,
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
	fastf_t draw_size)
{
    if (!parent_ctx)
	return 0;

    return ged_draw_scene_ref_commit_database_leaf_draft(
	    ged_draw_scene_ref_from_context(parent_ctx), gedp, view_ctx, dbip,
	    path, draw_mat, tol, ttol, settings, bool_op, rgb, has_draw_size,
	    draw_size);
}


static int
ged_draw_scene_ref_publish_current_mesh_lod(bsg_scene_ref ref)
{
    struct rt_mesh_lod_data data;
    struct ged_draw_source_runtime_summary source;
    memset(&source, 0, sizeof(source));
    if (ged_draw_scene_ref_is_null(ref) ||
	    !ged_draw_scene_ref_source_runtime_summary(ref, &source) ||
	    !rt_mesh_lod_data_get(source.rt_mesh_lod, &data))
	return 0;

    const int *faces = data.faces;
    const point_t *points = data.points;
    size_t face_count = data.face_count;
    size_t point_count = data.point_count;
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
    const point_t *normal_points = (data.points_orig && data.point_orig_count >= point_count) ?
	data.points_orig : points;
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
	    if (data.normals && data.normal_count >= (3*i + j + 1)) {
		VMOVE(normals[normal_out], data.normals[3*i + j]);
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

    int ret = ged_draw_scene_ref_update_indexed_face_set(ref,
	    points, point_count, normals, normal_count, indices, index_count);
    bu_free(normals, "mesh lod indexed-face normals");
    bu_free(indices, "mesh lod indexed-face indices");
    if (ret) {
	(void)ged_draw_scene_ref_source_mesh_lod_bounds_restore(ref);
	ged_draw_scene_ref_invalidate(ref);
    }
    return ret;
}


static struct ged *
_ged_draw_scene_ref_owner_gedp(bsg_scene_ref ref)
{
    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
    return shape_data ? shape_data->gedp : NULL;
}


static int
_ged_draw_scene_ref_user_data_set(
	bsg_scene_ref ref,
	void *data,
	ged_draw_shape_user_data_kind kind)
{
    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
    if (!shape_data)
	return 0;
    shape_data->u_data = data;
    shape_data->u_data_kind = kind;
    return 1;
}


static int
ged_draw_scene_ref_geometry_clear(bsg_scene_ref ref)
{
    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
    if (!shape_data)
	return 0;
    bsg_scene_ref geometry_ref =
	ged_draw_scene_ref_geometry_publication_target(ref);
    (void)bsg_geometry_ref_clear(bsg_scene_ref_as_geometry(geometry_ref));
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
    return ged_draw_scene_ref_last_point(
	    ged_draw_registry_shape_scene_ref(gedp, ref), out);
}


int
ged_draw_shape_ref_line_summary(struct ged *gedp,
				ged_draw_shape_ref ref,
				struct ged_draw_view_line_summary *out)
{
    return ged_draw_scene_ref_line_summary(
	    ged_draw_registry_shape_scene_ref(gedp, ref), out);
}


int
ged_draw_shape_ref_line_point_at(struct ged *gedp,
				 ged_draw_shape_ref ref,
				 size_t index,
				 point_t out)
{
    return ged_draw_scene_ref_line_point_at(
	    ged_draw_registry_shape_scene_ref(gedp, ref), index, out);
}


int
ged_draw_shape_ref_line_command_at(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   size_t index,
				   int *out)
{
    return ged_draw_scene_ref_line_command_at(
	    ged_draw_registry_shape_scene_ref(gedp, ref), index, out);
}


int
ged_draw_shape_context_line_summary(void *shape_ctx,
				    struct ged_draw_view_line_summary *out)
{
    return ged_draw_scene_ref_line_summary(
	    ged_draw_scene_ref_from_context(shape_ctx), out);
}


int
ged_draw_shape_context_line_point_at(void *shape_ctx,
				     size_t index,
				     point_t out)
{
    return ged_draw_scene_ref_line_point_at(
	    ged_draw_scene_ref_from_context(shape_ctx), index, out);
}


int
ged_draw_shape_context_line_command_at(void *shape_ctx,
				       size_t index,
				       int *out)
{
    return ged_draw_scene_ref_line_command_at(
	    ged_draw_scene_ref_from_context(shape_ctx), index, out);
}


int
ged_draw_shape_context_geometry_summary(void *shape_ctx,
					struct ged_draw_shape_geometry_summary *out)
{
    return ged_draw_scene_ref_geometry_summary(
	    ged_draw_scene_ref_from_context(shape_ctx), out);
}


int
ged_draw_shape_ref_translate_geometry(struct ged *gedp, ged_draw_shape_ref ref,
				      const vect_t xlate)
{
    return ged_draw_scene_ref_translate_geometry(
	    ged_draw_registry_shape_scene_ref(gedp, ref), xlate);
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
    bsg_scene_set_draw_center(ref, center);

    fastf_t size = bmax[X] - bmin[X];
    V_MAX(size, bmax[Y] - bmin[Y]);
    V_MAX(size, bmax[Z] - bmin[Z]);
    bsg_scene_set_draw_size(ref, size);
    ged_draw_scene_ref_set_non_database_source(ref, 0);
    return 1;
}


static void
ged_draw_scene_ref_set_draw_center(bsg_scene_ref ref, const point_t center)
{
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


static void
ged_draw_scene_ref_realization_finish_view(bsg_scene_ref ref,
					   void *view_ctx,
					   void *scene_view_ctx)
{
    (void)ged_draw_scene_ref_update_bounds_context(ref, view_ctx);
    if (scene_view_ctx)
	ged_draw_scene_ref_realization_set_view_context_policy(ref,
		scene_view_ctx);
}


static void
ged_draw_scene_ref_realization_finish_current_view(bsg_scene_ref ref,
						   void *view_ctx,
						   void *scene_view_ctx)
{
    ged_draw_scene_ref_realization_finish_view(ref, view_ctx, scene_view_ctx);
    ged_draw_scene_ref_realization_set_current(ref, 1);
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

    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
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


static void
_ged_draw_scene_ref_destroy_pair(bsg_scene_ref source_ref, bsg_scene_ref shape_ref)
{
    if (!bsg_scene_ref_is_null(source_ref) &&
	    !bsg_scene_ref_is_null(shape_ref))
	(void)ged_draw_scene_ref_detach(shape_ref);
    if (!bsg_scene_ref_is_null(shape_ref))
	ged_draw_scene_ref_release(shape_ref);
    if (!bsg_scene_ref_is_null(source_ref))
	ged_draw_scene_ref_release(source_ref);
}


static int
ged_draw_scene_ref_create_draft_pair(struct ged *gedp,
				     void *view_ctx,
				     bsg_scene_ref *source_out,
				     bsg_scene_ref *shape_out)
{
    if (source_out)
	*source_out = bsg_scene_ref_null();
    if (shape_out)
	*shape_out = bsg_scene_ref_null();

    if (!gedp || !view_ctx || !source_out || !shape_out)
	return 0;
    if (!ged_draw_ensure_root_attached(gedp))
	return 0;

    bsg_scene_ref source_scene =
	ged_draw_view_context_database_source_create(view_ctx, "_db_source");
    if (bsg_scene_ref_is_null(source_scene))
	return 0;

    bsg_scene_ref shape_scene =
	ged_draw_view_context_geometry_create(view_ctx, "geometry");
    if (bsg_scene_ref_is_null(shape_scene)) {
	_ged_draw_scene_ref_destroy_pair(source_scene, shape_scene);
	return 0;
    }
    if (!ged_draw_scene_ref_append_child(source_scene, shape_scene)) {
	_ged_draw_scene_ref_destroy_pair(source_scene, shape_scene);
	return 0;
    }

    if (!ged_draw_scene_context_ensure_registry_entry(gedp,
	    ged_draw_scene_ref_context(shape_scene))) {
	_ged_draw_scene_ref_destroy_pair(source_scene, shape_scene);
	return 0;
    }
    _ged_draw_scene_ref_set_source_ref(shape_scene, source_scene);
    ged_draw_scene_ref_geometry_clear(shape_scene);

    *source_out = source_scene;
    *shape_out = shape_scene;
    return 1;
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

    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
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
_ged_draw_scene_ref_apply_aux_display_state(bsg_scene_ref dst, bsg_scene_ref src)
{
    unsigned char rgb[3] = {0, 0, 0};
    unsigned char basecolor[3] = {0, 0, 0};
    mat_t mat;
    struct ged_draw_shape_material_summary material;

    if (ged_draw_scene_ref_is_null(dst) || ged_draw_scene_ref_is_null(src))
	return 0;

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
    ged_draw_scene_ref_set_non_database_source(dst, 0);
    return 1;
}


struct ged_draw_aux_display_sync_ctx {
    bsg_scene_ref shape_ref;
};


static int
_ged_draw_scene_ref_sync_aux_display_cb(bsg_scene_ref child_ref,
					void *userdata)
{
    struct ged_draw_aux_display_sync_ctx *ctx =
	(struct ged_draw_aux_display_sync_ctx *)userdata;

    if (!ctx || ged_draw_scene_ref_equal(child_ref, ctx->shape_ref))
	return 1;

    (void)_ged_draw_scene_ref_apply_aux_display_state(child_ref,
	    ctx->shape_ref);
    return 1;
}


static int
ged_draw_scene_ref_sync_aux_display_state(bsg_scene_ref source_ref,
					  bsg_scene_ref shape_ref)
{
    if (ged_draw_scene_ref_is_null(source_ref) ||
	    ged_draw_scene_ref_is_null(shape_ref))
	return 0;

    struct ged_draw_aux_display_sync_ctx ctx;
    ctx.shape_ref = shape_ref;
    return ged_draw_scene_ref_foreach_child(source_ref,
	    _ged_draw_scene_ref_sync_aux_display_cb, &ctx);
}


static int
ged_draw_scene_ref_set_name(bsg_scene_ref ref, const char *name)
{
    if (bsg_scene_ref_is_null(ref) || !name)
	return 0;

    bsg_scene_set_name(ref, name);
    return 1;
}


static int
ged_draw_scene_ref_set_visible(bsg_scene_ref ref, int visible)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    bsg_scene_set_visible(ref, visible ? 1 : 0);
    return 1;
}


static int
ged_draw_scene_ref_set_transparency(bsg_scene_ref ref, fastf_t transparency)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    bsg_scene_set_transparency(ref, transparency);
    return 1;
}


static int
ged_draw_scene_ref_set_color(bsg_scene_ref ref, const unsigned char rgb[3])
{
    if (bsg_scene_ref_is_null(ref) || !rgb)
	return 0;

    bsg_scene_set_color(ref, rgb[0], rgb[1], rgb[2]);
    return 1;
}


static int
ged_draw_scene_ref_set_highlighted(bsg_scene_ref ref, int highlighted)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

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

    bsg_scene_set_draw_mat(ref, mat);
    return 1;
}


static int
ged_draw_scene_ref_set_draw_size(bsg_scene_ref ref, fastf_t size)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    bsg_scene_set_draw_size(ref, size);
    return 1;
}


static int
ged_draw_scene_ref_set_line_style(bsg_scene_ref ref, int dashed)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    bsg_scene_set_line_style(ref, dashed);
    return 1;
}


static int
ged_draw_scene_ref_set_line_width(bsg_scene_ref ref, int line_width)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

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

    if (!ged_draw_scene_context_set_indexed_fullpath(gedp,
	    ged_draw_scene_ref_context(ref), path))
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

    struct bsg_draw_intent *di = bsg_scene_draw_intent(ref);
    if (di) {
	bsg_draw_intent_set_mode(di, (bsg_draw_mode)draw_mode);
	return 1;
    }

    return ged_draw_scene_ref_ensure_draw_intent(ref, fallback_path, draw_mode);
}


static int
ged_draw_scene_ref_set_draw_intent_appearance(
	bsg_scene_ref ref,
	const struct bsg_appearance_settings *settings,
	const char *fallback_path)
{
    if (bsg_scene_ref_is_null(ref) || !settings)
	return 0;

    struct bsg_draw_intent *di = bsg_scene_draw_intent(ref);
    if (!di) {
	if (!ged_draw_scene_ref_ensure_draw_intent(ref, fallback_path,
		settings->draw_mode))
	    return 0;
	di = bsg_scene_draw_intent(ref);
    }

    if (!di)
	return 0;

    bsg_draw_intent_set_appearance(di, settings);
    return 1;
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
ged_draw_bsg_appearance_from_neutral(struct bsg_appearance_settings *out,
				     const struct ged_draw_appearance_settings *in)
{
    if (!out)
	return;

    struct bsg_appearance_settings settings = BSG_APPEARANCE_SETTINGS_INIT;
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
ged_draw_scene_ref_set_draw_intent_appearance_settings(
	bsg_scene_ref ref,
	const struct ged_draw_appearance_settings *settings,
	const char *fallback_path)
{
    if (!settings)
	return 0;

    struct bsg_appearance_settings bsg_settings = BSG_APPEARANCE_SETTINGS_INIT;
    ged_draw_bsg_appearance_from_neutral(&bsg_settings, settings);
    return ged_draw_scene_ref_set_draw_intent_appearance(ref, &bsg_settings,
	    fallback_path);
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
ged_draw_scene_ref_strict_fallback(bsg_scene_ref ref)
{
    return bsg_scene_strict_fallback(ref);
}


static void
ged_draw_scene_ref_apply_qray_work_flag(bsg_scene_ref ref, int wflag)
{
    bsg_scene_set_work_flag(ref, wflag);
}


int
ged_draw_shape_ref_apply_qray_work_flag(struct ged *gedp, ged_draw_shape_ref ref, int wflag)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (ged_draw_scene_ref_is_null(shape_ref))
	return 0;
    ged_draw_scene_ref_apply_qray_work_flag(shape_ref, wflag);
    return 1;
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

    bsg_scene_set_legacy_eval_flag(ref, evaluated_region ? 1 : 0);
    return 1;
}


static int
ged_draw_scene_ref_set_legacy_uses_default_color(bsg_scene_ref ref,
						 int default_color)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    bsg_scene_set_legacy_uses_default_color(ref, default_color ? 1 : 0);
    return 1;
}


static int
ged_draw_scene_ref_set_legacy_region_id(bsg_scene_ref ref, int region_id)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

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
ged_draw_shape_draft_apply_db_object_marker(bsg_scene_ref ref)
{
    return _ged_draw_scene_ref_apply_db_object_marker(ref);
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


static int
ged_draw_scene_ref_bump_appearance_revision(bsg_scene_ref ref)
{
    return ged_draw_scene_ref_set_changed(ref,
	    ged_draw_scene_ref_changed(ref) + 1);
}


static int
ged_draw_scene_ref_apply_overlay_geometry_attributes(bsg_scene_ref ref,
						     const unsigned char rgb[3],
						     fastf_t transparency,
						     int draw_mode)
{
    if (bsg_scene_ref_is_null(ref) || !rgb)
	return 0;

    bsg_scene_set_highlighted(ref, 0);
    bsg_scene_set_line_style(ref, 0);
    bsg_scene_set_legacy_eval_flag(ref, 1);
    bsg_scene_set_legacy_color_info(ref, rgb, 0, 0);
    bsg_scene_set_color(ref, rgb[0], rgb[1], rgb[2]);
    bsg_scene_set_legacy_region_id(ref, 0);
    bsg_scene_set_work_flag(ref, 0);
    bsg_scene_set_transparency(ref, transparency);
    bsg_scene_set_dmode(ref, draw_mode);
    return 1;
}


static int
ged_draw_scene_ref_apply_settings(bsg_scene_ref ref,
				  const struct bsg_appearance_settings *settings)
{
    if (bsg_scene_ref_is_null(ref) || !settings)
	return 0;

    (void)bsg_scene_apply_appearance_settings(ref, settings);
    return 1;
}


static int
ged_draw_scene_ref_apply_display_settings(
	bsg_scene_ref ref,
	const struct ged_draw_appearance_settings *settings)
{
    if (bsg_scene_ref_is_null(ref) || !settings)
	return 0;

    if (!ged_draw_scene_ref_set_transparency(ref, settings->transparency))
	return 0;
    ged_draw_scene_ref_set_draw_mode(ref, settings->draw_mode);

    struct bsg_appearance_settings bsg_settings = BSG_APPEARANCE_SETTINGS_INIT;
    ged_draw_bsg_appearance_from_neutral(&bsg_settings, settings);
    return ged_draw_scene_ref_apply_settings(ref, &bsg_settings);
}


int
ged_draw_scene_context_attach_draw_bookkeeping(struct ged *gedp,
					       void *scene_ctx)
{
    bsg_scene_ref root_ref = ged_draw_scene_ref_from_context(scene_ctx);
    if (!gedp || !gedp->i || !gedp->i->ged_gdp ||
	    bsg_scene_ref_is_null(root_ref))
	return 0;

    struct ged_drawable *gdp = gedp->i->ged_gdp;
    gdp->bsg_ctx.draw_rev = &gdp->gd_draw_rev;
    gdp->bsg_ctx.fso = rt_view_set_context_recycle_pool(
	    ged_view_set_ctx(gedp));
    gdp->bsg_ctx.owner_data = gdp;
    bsg_scene_set_draw_ctx(root_ref, &gdp->bsg_ctx);
    return 1;
}


static void
ged_draw_scene_ref_set_material_rgb(bsg_scene_ref ref,
				    const unsigned char rgb[3])
{
    if (bsg_scene_ref_is_null(ref) || !rgb)
	return;
    bsg_scene_material_set_rgb(ref, rgb[0], rgb[1], rgb[2]);
}


static uint64_t
ged_draw_scene_ref_material_revision(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? 0 :
	(uint64_t)bsg_scene_material_revision(ref);
}


static int
ged_draw_scene_ref_set_material_revision(bsg_scene_ref ref,
					 uint64_t mater_rev)
{
    if (bsg_scene_ref_is_null(ref))
	return 0;

    bsg_scene_material_set_revision(ref, (uint32_t)mater_rev);
    return 1;
}


static int
ged_draw_scene_ref_realization_current(bsg_scene_ref ref)
{
    struct ged_draw_database_source_record record;
    memset(&record, 0, sizeof(record));
    if (!ged_draw_scene_ref_database_source_record_get(
	    ged_draw_scene_ref_source_owner(ref), &record))
	return 0;
    return record.realization_status ==
	GED_DRAW_DATABASE_SOURCE_REALIZATION_CURRENT;
}


static void
ged_draw_scene_ref_realization_set_current(bsg_scene_ref ref, int current)
{
    bsg_scene_ref source_ref = ged_draw_scene_ref_source_owner(ref);
    struct ged_draw_database_source_record record;
    memset(&record, 0, sizeof(record));
    (void)ged_draw_scene_ref_database_source_record_get(source_ref, &record);
    record.realization_status = current ?
	GED_DRAW_DATABASE_SOURCE_REALIZATION_CURRENT :
	GED_DRAW_DATABASE_SOURCE_REALIZATION_STALE;
    (void)ged_draw_scene_ref_database_source_record_apply(source_ref,
	    &record);
    if (!bsg_scene_ref_equal(source_ref, ref))
	(void)ged_draw_scene_ref_database_source_record_apply(ref, &record);
}


static void
ged_draw_scene_ref_realization_mark_current(bsg_scene_ref ref)
{
    bsg_scene_ref source_ref = ged_draw_scene_ref_source_owner(ref);
    struct ged_draw_database_source_record record;
    memset(&record, 0, sizeof(record));
    if (!ged_draw_scene_ref_database_source_record_get(source_ref, &record))
	return;

    record.realized_source_revision = record.source_revision;
    record.realized_inputs_revision = record.inputs_revision;
    record.stale_reason = GED_DRAW_STALE_NONE;
    (void)ged_draw_scene_ref_database_source_record_apply(source_ref,
	    &record);
    if (!bsg_scene_ref_equal(source_ref, ref))
	(void)ged_draw_scene_ref_database_source_record_apply(ref, &record);
}


static int
ged_draw_scene_ref_mark_database_source_changed(bsg_scene_ref ref,
						ged_draw_stale_reason reason)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 0;

    bsg_scene_ref source_ref = ged_draw_scene_ref_source_owner(ref);
    struct ged_draw_database_source_record record;
    memset(&record, 0, sizeof(record));
    (void)ged_draw_scene_ref_database_source_record_get(source_ref, &record);
    record.source_revision++;
    if (record.source_revision == 0)
	record.source_revision = 1;
    record.stale_reason = reason ? reason : GED_DRAW_STALE_SOURCE_CHANGED;
    (void)ged_draw_scene_ref_database_source_record_apply(source_ref,
	    &record);
    if (!ged_draw_scene_ref_equal(source_ref, ref))
	(void)ged_draw_scene_ref_database_source_record_apply(ref, &record);
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
}


static void
ged_draw_scene_ref_realization_clear_stale(bsg_scene_ref ref)
{
    ged_draw_scene_ref_geometry_clear(ref);
    ged_draw_scene_ref_realization_set_current(ref, 0);
}


static void
ged_draw_scene_ref_realization_set_roles(bsg_scene_ref ref, int csg_obj, int mesh_obj)
{
    int flags = GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_NONE;
    if (csg_obj)
	flags |= GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_CSG;
    if (mesh_obj)
	flags |= GED_DRAW_DATABASE_SOURCE_REALIZATION_ROLE_MESH;
    bsg_scene_ref source_ref = ged_draw_scene_ref_source_owner(ref);
    struct ged_draw_database_source_record record;
    memset(&record, 0, sizeof(record));
    (void)ged_draw_scene_ref_database_source_record_get(source_ref, &record);
    record.realization_role_flags = flags;
    (void)ged_draw_scene_ref_database_source_record_apply(source_ref,
	    &record);
    if (!bsg_scene_ref_equal(source_ref, ref))
	(void)ged_draw_scene_ref_database_source_record_apply(ref, &record);
}


static void
ged_draw_scene_ref_realization_prepare_wireframe(bsg_scene_ref ref)
{
    ged_draw_scene_ref_realization_set_roles(ref, 1, 0);
}


static void
ged_draw_scene_ref_realization_prepare_mesh(bsg_scene_ref ref)
{
    ged_draw_scene_ref_realization_set_roles(ref, 0, 1);
}


static void
ged_draw_scene_ref_realization_prepare_surface(bsg_scene_ref ref)
{
    ged_draw_scene_ref_realization_set_roles(ref, 0, 0);
}


static int
ged_draw_scene_ref_realization_prepare_view_redraw(bsg_scene_ref ref,
						   void *view_ctx)
{
    if (!view_ctx)
	return 0;
    (void)_ged_draw_scene_ref_mark_view_inputs_changed(ref);
    ged_draw_scene_ref_realization_set_view_context_policy(ref, view_ctx);
    return 1;
}


static int
_ged_draw_scene_ref_mark_view_inputs_changed(bsg_scene_ref ref)
{
    bsg_scene_ref source_ref = ged_draw_scene_ref_source_owner(ref);
    struct ged_draw_database_source_record record;
    memset(&record, 0, sizeof(record));
    if (!ged_draw_scene_ref_database_source_record_get(source_ref, &record))
	return 0;

    record.inputs_revision++;
    record.stale_reason = GED_DRAW_STALE_VIEW_INPUT_CHANGED;
    (void)ged_draw_scene_ref_database_source_record_apply(source_ref,
	    &record);
    if (!bsg_scene_ref_equal(source_ref, ref))
	(void)ged_draw_scene_ref_database_source_record_apply(ref, &record);
    return 1;
}


static int
ged_draw_scene_ref_realize_view_inputs_changed(bsg_scene_ref ref,
					       void *view_ctx)
{
    if (ged_draw_scene_ref_is_null(ref))
	return 0;

    (void)_ged_draw_scene_ref_mark_view_inputs_changed(ref);
    ged_draw_scene_ref_invalidate(ref);
    ged_draw_scene_ref_realize(ref, view_ctx);
    return 1;
}


static fastf_t
ged_draw_scene_ref_realization_view_scale(bsg_scene_ref ref)
{
    struct ged_draw_database_source_record record;
    memset(&record, 0, sizeof(record));
    if (!ged_draw_scene_ref_database_source_record_get(
	    ged_draw_scene_ref_source_owner(ref), &record))
	return 0.0;
    return record.realization_view_scale;
}


static fastf_t
ged_draw_scene_ref_realization_curve_scale(bsg_scene_ref ref)
{
    struct ged_draw_database_source_record record;
    memset(&record, 0, sizeof(record));
    if (!ged_draw_scene_ref_database_source_record_get(
	    ged_draw_scene_ref_source_owner(ref), &record))
	return 0.0;
    return record.realization_curve_scale;
}


static fastf_t
ged_draw_scene_ref_realization_point_scale(bsg_scene_ref ref)
{
    struct ged_draw_database_source_record record;
    memset(&record, 0, sizeof(record));
    if (!ged_draw_scene_ref_database_source_record_get(
	    ged_draw_scene_ref_source_owner(ref), &record))
	return 0.0;
    return record.realization_point_scale;
}


static void
ged_draw_scene_ref_realization_set_view_context_policy(bsg_scene_ref ref,
						       const void *view_ctx)
{
    if (!view_ctx)
	return;

    ged_draw_view_lod_policy policy;
    rt_view_context_lod_policy_get(&policy, view_ctx);

    bsg_scene_ref source_ref = ged_draw_scene_ref_source_owner(ref);
    struct ged_draw_database_source_record record;
    memset(&record, 0, sizeof(record));
    (void)ged_draw_scene_ref_database_source_record_get(source_ref, &record);
    record.realization_view_dependent = policy.csg_enabled ? 1 : 0;
    record.realization_view_scale = rt_view_context_scale_get(view_ctx);
    record.realization_bot_threshold = (uint64_t)policy.bot_threshold;
    record.realization_curve_scale = policy.curve_scale;
    record.realization_point_scale = policy.point_scale;
    (void)ged_draw_scene_ref_database_source_record_apply(source_ref,
	    &record);
    if (!bsg_scene_ref_equal(source_ref, ref))
	(void)ged_draw_scene_ref_database_source_record_apply(ref, &record);
}


static int
ged_draw_scene_ref_realization_pipeline_candidate(bsg_scene_ref ref)
{
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
						   const struct bn_tol *tol)
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

    return ged_draw_scene_ref_publish_realization_line_set(ref, &realization);
}


static int
ged_draw_scene_ref_publish_bspline_wireframe_line_set(bsg_scene_ref ref,
						      struct rt_db_internal *ip,
						      const struct bn_tol *tol)
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

    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
    if (!shape_data)
	return 0;

    struct rt_annot_internal *ann = (struct rt_annot_internal *)ip->idb_ptr;
    RT_ANNOT_CK_MAGIC(ann);

    point_t *points = NULL;
    struct bsg_annotation_segment *segments = NULL;
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

    if (ann->ant.count) {
	segments = (struct bsg_annotation_segment *)bu_calloc(ann->ant.count,
		sizeof(struct bsg_annotation_segment), "GED ANNOT segments");
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
    if (segments)
	bu_free(segments, "GED ANNOT segments");
    bu_vls_free(&summary);
    return ok;
}


static int
_ged_draw_scene_ref_publish_direct_primitive_wireframe(bsg_scene_ref ref,
						       struct rt_db_internal *ip,
						       const struct bg_tess_tol *ttol,
						       const struct bn_tol *tol,
						       struct bsg_view *v)
{
    int ok = 0;

    if (!ip)
	return 0;

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
		    (const struct rt_brep_internal *)ip->idb_ptr, tol);
	    break;
	case ID_BSPLINE:
	    ok = ged_draw_scene_ref_publish_bspline_wireframe_line_set(ref, ip,
		    tol);
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
	    ok = ged_draw_scene_ref_publish_submodel_wireframe_children(ref, ip,
		    ttol, tol, v);
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
    struct bsg_view *v = (struct bsg_view *)view_ctx;

    if (!adaptive || (ip && (ip->idb_type == ID_SUBMODEL ||
		    ip->idb_type == ID_ANNOT))) {
	int direct = _ged_draw_scene_ref_publish_direct_primitive_wireframe(
		ref, ip, ttol, tol, v);
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
ged_draw_scene_ref_publish_current_wireframe(bsg_scene_ref ref,
					     struct rt_db_internal *ip,
					     const struct bg_tess_tol *ttol,
					     const struct bn_tol *tol,
					     void *view_ctx,
					     const struct rt_view_info *view_info,
					     int adaptive)
{
    if (ged_draw_scene_ref_is_null(ref) || !ip)
	return 0;

    if (!adaptive)
	ged_draw_scene_ref_geometry_clear(ref);

    if (ged_draw_scene_ref_publish_primitive_wireframe(ref, ip, ttol, tol,
	    view_ctx, view_info, adaptive) < 0)
	return 0;

    if (adaptive) {
	ged_draw_scene_ref_update_bounds_context(ref, view_ctx);
	ged_draw_scene_ref_realization_mark_current(ref);
	ged_draw_scene_ref_invalidate(ref);
    } else {
	ged_draw_scene_ref_realization_set_current(ref, 1);
    }

    return 1;
}


static int
ged_draw_scene_ref_publish_line_set(bsg_scene_ref ref,
				    const point_t *points,
				    const int *commands,
				    size_t point_count)
{
    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);

    if (!shape_data)
	return 0;
    if (point_count && !points)
	return 0;
    bsg_scene_ref geometry_ref =
	ged_draw_scene_ref_geometry_publication_target(ref);
    if (bsg_scene_ref_is_null(geometry_ref))
	return 0;
    bsg_geometry_ref geometry = bsg_scene_ref_as_geometry(geometry_ref);
    if (!bsg_geometry_ref_set_line_set(geometry, points, commands,
	    point_count))
	return 0;

    shape_data->geometry_command_count = point_count;
    shape_data->geometry_revision++;
    bsg_scene_invalidate(ref);
    return 1;
}


static int
ged_draw_scene_ref_publish_point_set(bsg_scene_ref ref,
				     const point_t *points,
				     size_t point_count)
{
    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);

    if (!shape_data)
	return 0;
    if (point_count && !points)
	return 0;
    bsg_scene_ref geometry_ref =
	ged_draw_scene_ref_geometry_publication_target(ref);
    if (!bsg_geometry_ref_set_point_set(bsg_scene_ref_as_geometry(geometry_ref),
	    points, point_count))
	return 0;

    shape_data->geometry_command_count = point_count;
    shape_data->geometry_revision++;
    bsg_scene_invalidate(ref);
    return 1;
}


static int
ged_draw_scene_ref_publish_indexed_face_set(bsg_scene_ref ref,
					    const point_t *points,
					    size_t point_count,
					    const vect_t *normals,
					    size_t normal_count,
					    const int *indices,
					    size_t index_count)
{
    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);

    if (!shape_data)
	return 0;
    if (!point_count || !index_count)
	return ged_draw_scene_ref_geometry_clear(ref);
    if (!points || !indices)
	return 0;

    bsg_scene_ref geometry_ref =
	ged_draw_scene_ref_geometry_publication_target(ref);
    int ok = bsg_geometry_ref_set_indexed_face_set(bsg_scene_ref_as_geometry(geometry_ref),
	    points, point_count, normals, normal_count, indices, index_count);
    if (!ok)
	return 0;

    shape_data->geometry_command_count = point_count;
    shape_data->geometry_revision++;
    bsg_scene_invalidate(ref);
    return 1;
}


static int
ged_draw_scene_ref_update_indexed_face_set(bsg_scene_ref ref,
					   const point_t *points,
					   size_t point_count,
					   const vect_t *normals,
					   size_t normal_count,
					   const int *indices,
					   size_t index_count)
{
    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);

    if (!shape_data)
	return 0;
    if (!point_count || !index_count)
	return ged_draw_scene_ref_geometry_clear(ref);
    if (!points || !indices)
	return 0;

    bsg_scene_ref geometry_ref =
	ged_draw_scene_ref_geometry_publication_target(ref);
    int ok = bsg_geometry_ref_update_indexed_face_set(bsg_scene_ref_as_geometry(geometry_ref),
	    points, point_count, normals, normal_count, indices, index_count);
    if (!ok)
	return 0;

    shape_data->geometry_command_count = point_count;
    shape_data->geometry_revision++;
    bsg_scene_invalidate(ref);
    return 1;
}


static void
_ged_draw_primitive_indexed_face_set_free(struct rt_primitive_indexed_face_set *face_set)
{
    if (!face_set)
	return;
    if (face_set->points)
	bu_free(face_set->points, "primitive indexed-face points");
    if (face_set->normals)
	bu_free(face_set->normals, "primitive indexed-face normals");
    if (face_set->indices)
	bu_free(face_set->indices, "primitive indexed-face indices");
    memset(face_set, 0, sizeof(*face_set));
}


static int
ged_draw_scene_ref_publish_primitive_face_set(bsg_scene_ref ref,
					     struct rt_db_internal *ip,
					     const struct bg_tess_tol *ttol,
					     const struct bn_tol *tol,
					     const struct rt_view_info *view_info)
{
    struct rt_primitive_indexed_face_set face_set;
    struct rt_view_info default_view_info = RT_VIEW_INFO_INIT;
    int ret;
    int ok;

    memset(&face_set, 0, sizeof(face_set));

    if (bsg_scene_ref_is_null(ref) || !ip || !ip->idb_meth ||
	    !ip->idb_meth->ft_indexed_face_set)
	return 0;

    if (!view_info)
	view_info = &default_view_info;
    ret = ip->idb_meth->ft_indexed_face_set(&face_set, ip, ttol, tol,
	    view_info);
    if (ret != BRLCAD_OK) {
	_ged_draw_primitive_indexed_face_set_free(&face_set);
	return 0;
    }

    ok = face_set.point_count || face_set.index_count ?
	ged_draw_scene_ref_publish_indexed_face_set(ref,
		(const point_t *)face_set.points, face_set.point_count,
		(const vect_t *)face_set.normals, face_set.normal_count,
		face_set.indices, face_set.index_count) :
	ged_draw_scene_ref_geometry_clear(ref);

    _ged_draw_primitive_indexed_face_set_free(&face_set);
    return ok;
}


static int
ged_draw_scene_ref_publish_current_face_set(bsg_scene_ref ref,
					    struct rt_db_internal *ip,
					    const struct bg_tess_tol *ttol,
					    const struct bn_tol *tol,
					    const struct rt_view_info *view_info)
{
    ged_draw_scene_ref_realization_prepare_surface(ref);
    return ged_draw_scene_ref_publish_primitive_face_set(ref, ip, ttol, tol,
	    view_info);
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
    if (!ged_draw_scene_context_ensure_registry_entry(shape_data->gedp,
	    ged_draw_scene_ref_context(child_shape_ref))) {
	ged_draw_scene_ref_release(child_source_ref);
	return bsg_scene_ref_null();
    }

    _ged_draw_scene_ref_set_source_ref(child_shape_ref, child_source_ref);
    ged_draw_scene_ref_geometry_clear(child_shape_ref);

    _ged_draw_scene_ref_copy_display_state(child_source_ref, primary_ref);
    _ged_draw_scene_ref_copy_display_state(child_shape_ref, primary_ref);
    _ged_draw_scene_ref_apply_db_object_marker(child_shape_ref);
    ged_draw_scene_ref_set_non_database_source(child_source_ref, 0);
    ged_draw_scene_ref_set_non_database_source(child_shape_ref, 0);

    if (!ged_draw_scene_ref_append_child(parent_source_ref, child_source_ref)) {
	ged_draw_scene_ref_release(child_source_ref);
	return bsg_scene_ref_null();
    }
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


struct ged_draw_submodel_publish_ctx {
    bsg_scene_ref parent_ref;
    const struct bg_tess_tol *ttol;
    const struct bn_tol *tol;
    struct bsg_view *view;
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

    if (!ctx || ctx->failed || !tsp || !ip)
	return TREE_NULL;

    if (pathp && pathp->fp_len > 0) {
	path_name = db_path_to_string(pathp);
	if (path_name && path_name[0])
	    name = path_name;
    } else if (ip->idb_meth) {
	name = ip->idb_meth->ft_name;
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
	    ctx->tol ? ctx->tol : tsp->ts_tol, ctx->view, NULL, 0);
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

    (void)ged_draw_scene_ref_update_bounds_from_geometry(child_shape_ref, NULL);
    if (!ged_draw_scene_ref_is_null(child_source_ref))
	(void)ged_draw_scene_ref_update_bounds_context(child_source_ref,
		ctx->view);
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
						       struct bsg_view *v)
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
    ctx.view = v ? v : (struct bsg_view *)ged_draw_scene_ref_view_context(ref);
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
	(void)ged_draw_scene_ref_update_bounds_context(source_ref, ctx.view);
    (void)ged_draw_scene_ref_update_bounds_context(ref, ctx.view);
    return 1;
}


static bsg_scene_ref
_ged_draw_scene_ref_create_aux_geometry(bsg_scene_ref primary_ref,
					const char *name)
{
    bsg_scene_ref source_ref = ged_draw_scene_ref_source_owner(primary_ref);
    bsg_scene_ref aux_ref;
    void *view_ctx = ged_draw_scene_ref_publication_view_context(primary_ref);

    if (ged_draw_scene_ref_is_null(primary_ref) ||
	    ged_draw_scene_ref_is_null(source_ref) || !view_ctx)
	return bsg_scene_ref_null();

    aux_ref = ged_draw_view_context_geometry_create(view_ctx, name);
    if (ged_draw_scene_ref_is_null(aux_ref))
	return bsg_scene_ref_null();

    _ged_draw_scene_ref_copy_display_state(aux_ref, primary_ref);
    _ged_draw_scene_ref_apply_db_object_marker(aux_ref);
    ged_draw_scene_ref_set_non_database_source(aux_ref, 0);
    if (!ged_draw_scene_ref_append_child(source_ref, aux_ref)) {
	ged_draw_scene_ref_release(aux_ref);
	return bsg_scene_ref_null();
    }
    return aux_ref;
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
_ged_draw_nmg_normal_line_stats(const struct nmgregion *r,
				int style,
				size_t *point_count)
{
    const struct shell *s;
    size_t points = 0;

    if (point_count)
	*point_count = 0;
    if (!r || !_ged_draw_nmg_surface_style(style) ||
	    !(style & GED_DRAW_NMG_STYLE_VISUALIZE_NORMALS))
	return 0;

    NMG_CK_REGION(r);
    for (BU_LIST_FOR(s, shell, &r->s_hd)) {
	const struct faceuse *fu;
	NMG_CK_SHELL(s);
	for (BU_LIST_FOR(fu, faceuse, &s->fu_hd)) {
	    const struct loopuse *lu;
	    if (fu->orientation != OT_SAME)
		continue;
	    NMG_CK_FACEUSE(fu);
	    NMG_CK_FACE(fu->f_p);
	    if (!fu->f_p->g.magic_p ||
		    *fu->f_p->g.magic_p == NMG_FACE_G_SNURB_MAGIC)
		return 0;
	    NMG_CK_FACE_G_PLANE(fu->f_p->g.plane_p);

	    for (BU_LIST_FOR(lu, loopuse, &fu->lu_hd)) {
		size_t loop_points = 0;
		if (!_ged_draw_nmg_loop_normal_line_count(lu, &loop_points))
		    return 0;
		points += loop_points;
	    }
	}
    }

    if (!points)
	return 0;
    if (point_count)
	*point_count = points;
    return 1;
}


static int
_ged_draw_nmg_wire_remainder_line_stats(const struct nmgregion *r,
					size_t *point_count)
{
    const struct shell *s;
    size_t points = 0;

    if (point_count)
	*point_count = 0;
    if (!r)
	return 0;

    NMG_CK_REGION(r);
    for (BU_LIST_FOR(s, shell, &r->s_hd)) {
	const struct loopuse *lu;
	size_t edge_points = 0;

	NMG_CK_SHELL(s);
	for (BU_LIST_FOR(lu, loopuse, &s->lu_hd)) {
	    size_t loop_points = 0;
	    if (!_ged_draw_nmg_loop_line_count(lu, NULL, &loop_points))
		return 0;
	    points += loop_points;
	}

	if (!_ged_draw_nmg_wire_edge_line_count(&s->eu_hd, &edge_points))
	    return 0;
	points += edge_points;

	if (s->vu_p) {
	    point_t point;
	    if (!_ged_draw_nmg_vertexuse_coord(s->vu_p, point))
		return 0;
	    points += 2;
	}
    }

    if (point_count)
	*point_count = points;
    return 1;
}


static int
_ged_draw_nmg_fill_normal_line_set(const struct nmgregion *r,
				   int style,
				   point_t *points,
				   int *commands,
				   size_t point_count)
{
    const struct shell *s;
    size_t point_idx = 0;

    if (!r || !points || !commands || !_ged_draw_nmg_surface_style(style) ||
	    !(style & GED_DRAW_NMG_STYLE_VISUALIZE_NORMALS))
	return 0;

    for (BU_LIST_FOR(s, shell, &r->s_hd)) {
	const struct faceuse *fu;
	for (BU_LIST_FOR(fu, faceuse, &s->fu_hd)) {
	    const struct loopuse *lu;
	    vect_t face_normal = VINIT_ZERO;

	    if (fu->orientation != OT_SAME)
		continue;
	    NMG_GET_FU_NORMAL(face_normal, fu);
	    for (BU_LIST_FOR(lu, loopuse, &fu->lu_hd)) {
		if (!_ged_draw_nmg_append_loop_normal_lines(points, commands,
			    point_count, &point_idx, lu, face_normal))
		    return 0;
	    }
	}
    }

    return (point_idx == point_count) ? 1 : 0;
}


static int
_ged_draw_nmg_fill_wire_remainder_line_set(const struct nmgregion *r,
					   point_t *points,
					   int *commands,
					   size_t point_count)
{
    const struct shell *s;
    size_t point_idx = 0;

    if (!r || !points || !commands)
	return 0;

    for (BU_LIST_FOR(s, shell, &r->s_hd)) {
	const struct loopuse *lu;

	for (BU_LIST_FOR(lu, loopuse, &s->lu_hd)) {
	    if (!_ged_draw_nmg_append_loop_lines(points, commands, point_count,
			&point_idx, lu, NULL))
		return 0;
	}

	if (!_ged_draw_nmg_append_wire_edge_lines(points, commands, point_count,
		    &point_idx, &s->eu_hd))
	    return 0;

	if (s->vu_p) {
	    point_t point;
	    if (!_ged_draw_nmg_vertexuse_coord(s->vu_p, point))
		return 0;
	    if (!_ged_draw_nmg_append_line(points, commands, point_count,
			&point_idx, point, BSG_GEOMETRY_LINE_MOVE))
		return 0;
	    if (!_ged_draw_nmg_append_line(points, commands, point_count,
			&point_idx, point, BSG_GEOMETRY_LINE_DRAW))
		return 0;
	}
    }

    return (point_idx == point_count) ? 1 : 0;
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
_ged_draw_nmg_fill_line_set(const struct nmgregion *r,
			    int style,
			    point_t *points,
			    int *commands,
			    size_t point_count)
{
    size_t point_idx = 0;

    if (!_ged_draw_nmg_append_region_line_set(r, style, points, commands,
	    point_count, &point_idx))
	return 0;
    return (point_idx == point_count) ? 1 : 0;
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
_ged_draw_nmg_loop_vertex_count(const struct loopuse *lu, size_t *count)
{
    const struct edgeuse *eu;
    size_t n = 0;

    if (count)
	*count = 0;
    if (!lu)
	return 0;

    NMG_CK_LOOPUSE(lu);
    if (BU_LIST_FIRST_MAGIC(&lu->down_hd) != NMG_EDGEUSE_MAGIC)
	return 0;

    for (BU_LIST_FOR(eu, edgeuse, &lu->down_hd)) {
	NMG_CK_EDGEUSE(eu);
	NMG_CK_VERTEXUSE(eu->vu_p);
	NMG_CK_VERTEX(eu->vu_p->v_p);
	if (!eu->vu_p->v_p->vg_p)
	    return 0;
	NMG_CK_VERTEX_G(eu->vu_p->v_p->vg_p);
	n++;
    }

    if (n < 3)
	return 0;
    if (count)
	*count = n;
    return 1;
}


static int
_ged_draw_nmg_indexed_face_stats(const struct nmgregion *r,
				 int style,
				 size_t *point_count,
				 size_t *index_count)
{
    const struct shell *s;
    size_t points = 0;
    size_t indices = 0;

    if (point_count)
	*point_count = 0;
    if (index_count)
	*index_count = 0;
    if (!r)
	return 0;

    if (!_ged_draw_nmg_surface_style(style))
	return 0;

    NMG_CK_REGION(r);
    for (BU_LIST_FOR(s, shell, &r->s_hd)) {
	const struct faceuse *fu;

	NMG_CK_SHELL(s);
	for (BU_LIST_FOR(fu, faceuse, &s->fu_hd)) {
	    const struct loopuse *lu;

	    NMG_CK_FACEUSE(fu);
	    if (fu->orientation != OT_SAME)
		continue;
	    NMG_CK_FACE(fu->f_p);
	    if (!fu->f_p->g.magic_p ||
		    *fu->f_p->g.magic_p == NMG_FACE_G_SNURB_MAGIC)
		return 0;
	    NMG_CK_FACE_G_PLANE(fu->f_p->g.plane_p);

	    for (BU_LIST_FOR(lu, loopuse, &fu->lu_hd)) {
		size_t n = 0;
		if (!_ged_draw_nmg_loop_vertex_count(lu, &n))
		    return 0;
		points += n;
		indices += n + 1;
	    }
	}
    }

    if (!points || !indices)
	return 0;
    if (point_count)
	*point_count = points;
    if (index_count)
	*index_count = indices;
    return 1;
}


static int
_ged_draw_nmg_fill_indexed_face_set(const struct nmgregion *r,
				    int style,
				    point_t *points,
				    vect_t *normals,
				    int *indices,
				    size_t point_count,
				    size_t index_count)
{
    const struct shell *s;
    size_t point_idx = 0;
    size_t index_idx = 0;

    if (!r || !points || !normals || !indices)
	return 0;

    for (BU_LIST_FOR(s, shell, &r->s_hd)) {
	const struct faceuse *fu;

	for (BU_LIST_FOR(fu, faceuse, &s->fu_hd)) {
	    const struct loopuse *lu;
	    vect_t face_normal = VINIT_ZERO;

	    if (fu->orientation != OT_SAME)
		continue;
	    NMG_GET_FU_NORMAL(face_normal, fu);

	    for (BU_LIST_FOR(lu, loopuse, &fu->lu_hd)) {
		const struct edgeuse *eu;

		for (BU_LIST_FOR(eu, edgeuse, &lu->down_hd)) {
		    const struct vertexuse *vu = eu->vu_p;
		    vect_t normal = VINIT_ZERO;

		    if (point_idx >= point_count || index_idx >= index_count)
			return 0;

		    VMOVE(points[point_idx], vu->v_p->vg_p->coord);
		    if ((style & GED_DRAW_NMG_STYLE_USE_VU_NORMALS) &&
			    vu->a.magic_p &&
			    *vu->a.magic_p == NMG_VERTEXUSE_A_PLANE_MAGIC) {
			VMOVE(normal, vu->a.plane_p->N);
		    } else {
			VMOVE(normal, face_normal);
		    }
		    VMOVE(normals[point_idx], normal);
		    indices[index_idx++] = (int)point_idx++;
		}

		if (index_idx >= index_count)
		    return 0;
		indices[index_idx++] = -1;
	    }
	}
    }

    return (point_idx == point_count && index_idx == index_count) ? 1 : 0;
}


static int
_ged_draw_scene_ref_publish_nmg_wire_remainders(bsg_scene_ref primary_ref,
						const struct nmgregion *r)
{
    bsg_scene_ref wire_ref = bsg_scene_ref_null();
    point_t *points = NULL;
    int *commands = NULL;
    size_t point_count = 0;
    int ok = 0;

    if (!_ged_draw_nmg_wire_remainder_line_stats(r, &point_count))
	return 0;
    if (!point_count)
	return 1;

    points = (point_t *)bu_calloc(point_count, sizeof(point_t),
	    "GED draw NMG mixed-wire points");
    commands = (int *)bu_calloc(point_count, sizeof(int),
	    "GED draw NMG mixed-wire commands");

    if (!_ged_draw_nmg_fill_wire_remainder_line_set(r, points, commands,
		point_count))
	goto cleanup;

    wire_ref = _ged_draw_scene_ref_create_aux_geometry(primary_ref,
	    "nmg_mixed_wire");
    if (bsg_scene_ref_is_null(wire_ref))
	goto cleanup;

    ok = bsg_geometry_ref_set_line_set(bsg_scene_ref_as_geometry(wire_ref),
	    (const point_t *)points, commands, point_count);
    if (!ok)
	goto cleanup;

    (void)ged_draw_scene_ref_update_bounds_from_geometry(wire_ref, NULL);
    ged_draw_scene_ref_invalidate(wire_ref);

cleanup:
    if (!ok && !ged_draw_scene_ref_is_null(wire_ref)) {
	(void)ged_draw_scene_ref_detach(wire_ref);
	ged_draw_scene_ref_release(wire_ref);
    }
    if (points)
	bu_free(points, "GED draw NMG mixed-wire points");
    if (commands)
	bu_free(commands, "GED draw NMG mixed-wire commands");
    return ok;
}


static int
_ged_draw_scene_ref_publish_nmg_normal_lines(bsg_scene_ref primary_ref,
					     const struct nmgregion *r,
					     int style)
{
    bsg_scene_ref normal_ref = bsg_scene_ref_null();
    point_t *points = NULL;
    int *commands = NULL;
    size_t point_count = 0;
    int ok = 0;

    if (!_ged_draw_nmg_normal_line_stats(r, style, &point_count))
	return 0;

    points = (point_t *)bu_calloc(point_count, sizeof(point_t),
	    "GED draw NMG normal-line points");
    commands = (int *)bu_calloc(point_count, sizeof(int),
	    "GED draw NMG normal-line commands");

    if (!_ged_draw_nmg_fill_normal_line_set(r, style, points, commands,
		point_count))
	goto cleanup;

    normal_ref = _ged_draw_scene_ref_create_aux_geometry(primary_ref,
	    "nmg_surface_normals");
    if (bsg_scene_ref_is_null(normal_ref))
	goto cleanup;

    ok = bsg_geometry_ref_set_line_set(bsg_scene_ref_as_geometry(normal_ref),
	    (const point_t *)points, commands, point_count);
    if (!ok)
	goto cleanup;

    (void)ged_draw_scene_ref_update_bounds_from_geometry(normal_ref, NULL);
    ged_draw_scene_ref_invalidate(normal_ref);

cleanup:
    if (!ok && !ged_draw_scene_ref_is_null(normal_ref)) {
	(void)ged_draw_scene_ref_detach(normal_ref);
	ged_draw_scene_ref_release(normal_ref);
    }
    if (points)
	bu_free(points, "GED draw NMG normal-line points");
    if (commands)
	bu_free(commands, "GED draw NMG normal-line commands");
    return ok;
}


static int
ged_draw_scene_ref_geometry_publish_nmg_region(bsg_scene_ref ref,
					       const struct nmgregion *r,
					       int style)
{
    ged_draw_shape_state *shape_data = ged_draw_shape_state_get_scene_ref(ref);
    point_t *points = NULL;
    int *commands = NULL;
    vect_t *normals = NULL;
    int *indices = NULL;
    size_t point_count = 0;
    size_t index_count = 0;
    int ok = 0;

    if (!shape_data || !r)
	return 0;

    if (_ged_draw_nmg_indexed_face_stats(r, style, &point_count, &index_count)) {
	points = (point_t *)bu_calloc(point_count, sizeof(point_t),
		"GED draw NMG indexed-face points");
	normals = (vect_t *)bu_calloc(point_count, sizeof(vect_t),
		"GED draw NMG indexed-face normals");
	indices = (int *)bu_calloc(index_count, sizeof(int),
		"GED draw NMG indexed-face indices");

	if (!_ged_draw_nmg_fill_indexed_face_set(r, style, points, normals, indices,
		    point_count, index_count))
	    goto cleanup;

	ok = bsg_geometry_ref_set_indexed_face_set(bsg_scene_ref_as_geometry(ref),
		(const point_t *)points, point_count,
		(const vect_t *)normals, point_count,
		indices, index_count);
	if (!ok)
	    goto cleanup;

	shape_data->geometry_command_count = point_count;
	shape_data->geometry_revision++;
	bsg_scene_invalidate(ref);
	if ((style & GED_DRAW_NMG_STYLE_VISUALIZE_NORMALS) &&
		!_ged_draw_scene_ref_publish_nmg_normal_lines(ref, r, style)) {
	    ok = 0;
	    goto cleanup;
	}
	if (!_ged_draw_scene_ref_publish_nmg_wire_remainders(ref, r)) {
	    ok = 0;
	    goto cleanup;
	}
	goto cleanup;
    }

    if (!_ged_draw_nmg_line_set_stats(r, style, &point_count))
	return 0;

    points = (point_t *)bu_calloc(point_count, sizeof(point_t),
	    "GED draw NMG line-set points");
    commands = (int *)bu_calloc(point_count, sizeof(int),
	    "GED draw NMG line-set commands");

    if (!_ged_draw_nmg_fill_line_set(r, style, points, commands, point_count))
	goto cleanup;

    ok = bsg_geometry_ref_set_line_set(bsg_scene_ref_as_geometry(ref),
	    (const point_t *)points, commands, point_count);
    if (!ok)
	goto cleanup;

    shape_data->geometry_command_count = point_count;
    shape_data->geometry_revision++;
    bsg_scene_invalidate(ref);

cleanup:
    if (points)
	bu_free(points, "GED draw NMG geometry points");
    if (commands)
	bu_free(commands, "GED draw NMG line-set commands");
    if (normals)
	bu_free(normals, "GED draw NMG indexed-face normals");
    if (indices)
	bu_free(indices, "GED draw NMG indexed-face indices");
    return ok;
}


static int
ged_draw_scene_ref_publish_current_nmg_region(bsg_scene_ref ref,
					     const struct nmgregion *r,
					     int style)
{
    if (!ged_draw_scene_ref_geometry_publish_nmg_region(ref, r, style))
	return 0;
    ged_draw_scene_ref_realization_set_current(ref, 1);
    return 1;
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
	ged_draw_shape_state_get_scene_ref(ref);
    if (!shape_data)
	return;
    shape_data->geometry_command_count = vlen;
    shape_data->geometry_revision++;
}


static int
ged_draw_source_primitive_has_lod_realize(const struct rt_db_internal *ip)
{
    return (ip && ip->idb_meth && ip->idb_meth->ft_lod_realize) ? 1 : 0;
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
