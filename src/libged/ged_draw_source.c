/*                  G E D _ D R A W _ S O U R C E . C
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
/** @file ged_draw_source.c
 *
 * GED draw-source metadata helpers backed by semantic draw records and Obol
 * scene-controller state, with legacy backend adapters confined to local
 * implementation details.
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
#include "bv.h"
#include "ged/draw.h"
#include "ged/db_index.h"
#include "ged/view.h"
#include "nmg.h"
#include "rt/db_instance.h"
#include "rt/db_io.h"
#include "rt/calc.h"
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
#include "./ged_draw_private.h"
#include "./ged_draw_shape_state_private.h"
#include "./ged_draw_view_private.h"


struct ged_draw_source_state;
struct ged_draw_source_state_record;
struct ged_draw_source_runtime_summary;
struct ged_draw_database_source_record;
struct ged_draw_scene_draw_state_summary;

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
static void ged_draw_obol_owner_structural_revision_bump(struct ged *gedp);

struct _ged_draw_obol_update_display_ctx {
    int visible_valid;
    int visible;
    int selected_valid;
    int selected;
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
typedef int (*ged_draw_scene_handle_child_cb)(ged_draw_scene_handle child_ref,
					      void *userdata);
static int ged_draw_scene_handle_foreach_child(
	ged_draw_scene_handle scene_ref,
	ged_draw_scene_handle_child_cb cb,
	void *userdata);
static void *ged_draw_obol_scene_context_child_at(
	ged_draw_scene_handle scene_ref,
	void *parent_ctx,
	size_t index);
static int ged_draw_source_scene_handle_foreach_shape(
	ged_draw_scene_handle scene_ref,
	ged_draw_scene_handle_child_cb cb,
	void *userdata);
static int ged_draw_scene_handle_display_summary(
	ged_draw_scene_handle scene_ref,
	struct ged_draw_scene_display_summary *out);
static int ged_draw_scene_handle_source_summary(
	ged_draw_scene_handle scene_ref,
	struct ged_draw_database_source_summary *out);
static int ged_draw_scene_handle_tree_summary(
	ged_draw_scene_handle scene_ref,
	struct ged_draw_scene_tree_summary *out);
static struct ged *ged_draw_scene_handle_owner_gedp(
	ged_draw_scene_handle scene_ref);
static const char *ged_draw_scene_handle_semantic_path(
	ged_draw_scene_handle scene_ref);
static const struct db_full_path *ged_draw_scene_handle_fullpath(
	ged_draw_scene_handle scene_ref);
static const char *ged_draw_scene_handle_name(ged_draw_scene_handle scene_ref);
static int ged_draw_scene_handle_subtree_bounds(
	ged_draw_scene_handle scene_ref,
	vect_t *min,
	vect_t *max,
	int include_overlays);
static ged_draw_scene_handle _ged_draw_source_root_scene_handle(struct ged *gedp);
static int ged_draw_rt_internal_payload_valid(const struct rt_db_internal *ip);
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
    const struct bn_tol *tol,
    void *view_ctx,
    fastf_t draw_size);
static void ged_draw_primitive_realization_line_set_free(
	struct rt_primitive_lod_realization *realization);
static int ged_draw_obol_database_source_publish_realization_line_set(
	struct ged *gedp,
	const char *path,
	struct rt_primitive_lod_realization *realization);
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
    int realization_csg_lod_enabled;
    int realization_mesh_lod_enabled;
    fastf_t realization_view_scale;
    fastf_t realization_lod_scale;
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
	point_t *points,
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
_ged_draw_view_export_base_instance_key(struct bu_vls *out,
					const char *instance_key)
{
    if (!out)
	return;
    bu_vls_trunc(out, 0);
    if (!instance_key || !instance_key[0])
	return;

    const char mode_marker[] = ":ged-draw-mode:";
    const char *marker = NULL;
    const char *candidate = instance_key;
    while ((candidate = strstr(candidate, mode_marker)) != NULL) {
	marker = candidate;
	candidate++;
    }
    if (marker)
	bu_vls_strncpy(out, instance_key, (size_t)(marker - instance_key));
    else
	bu_vls_strcpy(out, instance_key);
}


static int
_ged_draw_view_export_path_equal(const char *a, const char *b)
{
    if (!a || !b)
	return 0;
    while (*a == '/')
	a++;
    while (*b == '/')
	b++;
    return BU_STR_EQUAL(a, b);
}


static int
_ged_draw_view_export_obol_source_record_in_view(
	const struct ged_draw_obol_database_source_record *record,
	void *view_ctx)
{
    if (!record || !record->database_path || !record->database_path[0])
	return 0;

    const char *instance_key = record->instance_key;
    if (!view_ctx || !ged_view_context_is_independent(view_ctx)) {
	if (!instance_key || !instance_key[0])
	    return 1;

	struct bu_vls base_key = BU_VLS_INIT_ZERO;
	_ged_draw_view_export_base_instance_key(&base_key, instance_key);
	const char *base = bu_vls_cstr(&base_key);
	int ret = 0;
	if (!base || !base[0])
	    ret = 1;
	else if (bu_strncmp(base, "ged-view:", 9) == 0)
	    ret = 0;
	else if (bu_strncmp(base, "brlcad-direct:", 14) == 0)
	    ret = 1;
	else
	    ret = _ged_draw_view_export_path_equal(base,
		    record->database_path);
	bu_vls_free(&base_key);
	return ret;
    }

    const char *view_name = bv_name_get(
			       bv_context_view_const((const struct bv_context *)view_ctx));
    char fallback[64] = {0};
    if (!view_name || !view_name[0]) {
	snprintf(fallback, sizeof(fallback), "%p", view_ctx);
	view_name = fallback;
    }

    struct bu_vls base_key = BU_VLS_INIT_ZERO;
    struct bu_vls prefix = BU_VLS_INIT_ZERO;
    _ged_draw_view_export_base_instance_key(&base_key, instance_key);
    bu_vls_printf(&prefix, "ged-view:%s:", view_name);
    int ret = (bu_vls_strlen(&base_key) > 0 &&
	    bu_strncmp(bu_vls_cstr(&base_key), bu_vls_cstr(&prefix),
		bu_vls_strlen(&prefix)) == 0) ? 1 : 0;
    bu_vls_free(&prefix);
    bu_vls_free(&base_key);
    return ret;
}


static int
_ged_draw_view_export_obol_source_record_matches(
	const struct ged_draw_obol_view_export_ctx *ctx,
	const struct ged_draw_obol_database_source_record *record)
{
    if (!ctx || !record || !record->valid || !record->database_path ||
	    !record->database_path[0])
	return 0;

    if (!_ged_draw_view_export_obol_source_record_in_view(record,
	    ctx->view_ctx))
	return 0;

    if (((ctx->query_flags & GED_DRAW_VIEW_EXPORT_QUERY_VISIBLE_ONLY) ||
		(ctx->render_flags & GED_DRAW_VIEW_EXPORT_RENDER_VISIBLE_ONLY)) &&
	    !record->visible)
	return 0;

    if (ctx->draw_mode != GED_DRAW_VIEW_EXPORT_DRAW_MODE_ANY &&
	    record->draw_mode != ctx->draw_mode)
	return 0;

    const int wants_database =
	(ctx->query_flags & GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS) ? 1 : 0;
    const int wants_view =
	(ctx->query_flags & GED_DRAW_VIEW_EXPORT_QUERY_VIEW_OBJECTS) ? 1 : 0;
    if ((ctx->query_flags & GED_DRAW_VIEW_EXPORT_QUERY_LOCAL_ONLY) ||
	    ((wants_database || wants_view) && !wants_database))
	return 0;

    return _ged_draw_view_export_glob_match(ctx->glob,
	    record->database_path, NULL);
}


static int
_ged_draw_view_export_detail_from_obol_source_record(
	struct ged_draw_view_export_detail *detail,
	const struct ged_draw_obol_view_export_ctx *ctx,
	const struct ged_draw_obol_database_source_record *source_record)
{
    if (!detail || !ctx || !ctx->gedp || !source_record ||
	    !source_record->database_path || !source_record->database_path[0])
	return 0;

    _ged_draw_view_export_detail_init(detail);

    struct ged_draw_shape_geometry_summary geometry;
    memset(&geometry, 0, sizeof(geometry));
    int geometry_valid =
	ged_draw_obol_database_source_geometry_summary_for_path_mode(
	    ctx->gedp, source_record->database_path, 1,
	    source_record->draw_mode, &geometry);
    if (!geometry_valid || !geometry.valid || !geometry.geometry_name) {
	memset(&geometry, 0, sizeof(geometry));
	geometry_valid =
	    ged_draw_obol_database_source_geometry_summary_for_path(
		ctx->gedp, source_record->database_path, &geometry);
    }
    if (!geometry_valid || !geometry.valid || !geometry.geometry_name) {
	_ged_draw_view_export_detail_free(detail);
	return 0;
    }

    bu_vls_strcpy(&detail->path, source_record->database_path);

    struct ged_draw_view_db_object_record *rec = &detail->record;
    rec->path = bu_vls_cstr(&detail->path);
    rec->type_name =
	_ged_draw_view_export_type_name_from_obol_geometry(
		geometry.geometry_name);
    rec->geometry_name = geometry.geometry_name;
    rec->draw_mode = source_record->draw_mode;
    rec->transparency = source_record->transparency;
    rec->is_database_source = 1;
    rec->non_database_source = 0;
    rec->is_database_intent = 1;
    rec->is_local_source = 0;
    rec->is_view_source = 0;
    rec->highlighted = source_record->highlighted;
    rec->selected = source_record->selected;
    if (ctx->view_ctx && !rec->selected)
	rec->selected = ged_draw_view_context_selection_contains_path(
		ctx->view_ctx, GED_DRAW_VIEW_SELECTION_SELECTED_PATH,
		source_record->database_path) ? 1 : 0;
    rec->visible = source_record->visible;
    rec->line_style = 0;

    struct ged_draw_scene_display_summary display;
    memset(&display, 0, sizeof(display));
    if (ged_draw_obol_database_source_display_summary_for_path(ctx->gedp,
	    source_record->database_path, &display) && display.valid &&
	    display.material_valid) {
	rec->color[0] = display.material_color[0];
	rec->color[1] = display.material_color[1];
	rec->color[2] = display.material_color[2];
    }

    MAT_IDN(rec->model_mat);
    rec->cache_identity =
	_ged_draw_view_export_hash_cstr(source_record->instance_key);
    if (!rec->cache_identity)
	rec->cache_identity =
	    _ged_draw_view_export_hash_cstr(source_record->database_path);
    rec->source_identity =
	_ged_draw_view_export_hash_cstr(source_record->database_path);
    if (!rec->source_identity)
	rec->source_identity = rec->cache_identity;
    rec->detail_token = (uintptr_t)detail;

    if (!_ged_draw_view_export_detail_geometry_from_obol(detail, ctx->gedp,
	    source_record->database_path, geometry.geometry_name, 1)) {
	_ged_draw_view_export_detail_free(detail);
	return 0;
    }

    return 1;
}


static int
_ged_draw_view_export_obol_source_record_cb(
	struct ged *UNUSED(gedp),
	const struct ged_draw_obol_database_source_record *record,
	void *userdata)
{
    struct ged_draw_obol_view_export_ctx *ctx =
	(struct ged_draw_obol_view_export_ctx *)userdata;
    if (!ctx || !ctx->keep_going)
	return 0;
    if (!_ged_draw_view_export_obol_source_record_matches(ctx, record))
	return 1;

    struct ged_draw_view_export_detail detail;
    if (_ged_draw_view_export_detail_from_obol_source_record(&detail, ctx,
	    record)) {
	ctx->keep_going = ctx->cb(&detail.record, ctx->userdata);
	_ged_draw_view_export_detail_free(&detail);
    }

    return ctx->keep_going;
}


static void
_ged_draw_view_export_source_records_from_obol(
	struct ged_draw_obol_view_export_ctx *ctx,
	int *keep_going)
{
    if (!ctx || !keep_going || !*keep_going)
	return;

    ctx->keep_going = *keep_going;
    if (ged_draw_obol_database_source_records_foreach(ctx->gedp, 1,
	    _ged_draw_view_export_obol_source_record_cb, ctx) >= 0)
	*keep_going = ctx->keep_going;
}


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
	ged_draw_view_polygon_ref UNUSED(ref),
	const struct ged_draw_view_polygon_record *polygon,
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
    ged_draw_view_context_polygon_visit_records(view_ctx,
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

    _ged_draw_view_export_source_records_from_obol(&ctx, keep_going);
    if (!*keep_going) {
	ged_draw_obol_scene_context_info_free(&root_info);
	return;
    }
    ctx.keep_going = 1;

    _ged_draw_view_export_visit_obol_info(&ctx, &root_info, NULL);
    *keep_going = ctx.keep_going;
    ged_draw_obol_scene_context_info_free(&root_info);

    _ged_draw_view_export_polygon_records_from_rt(view_ctx, query_flags,
	    glob, draw_mode, cb, userdata, keep_going);
    _ged_draw_view_export_feature_records_from_obol(view_ctx, query_flags,
	    glob, cb, userdata, keep_going);
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
    if (!view_ctx || !cb)
	return;

    struct ged *gedp = (struct ged *)ged_view_context_user_data_get(view_ctx);
    if (!ged_draw_obol_scene_controller_is_owned(gedp))
	return;

    int keep_going = 1;
    _ged_draw_view_export_records_from_obol(gedp, view_ctx, query_flags,
	render_flags, glob, draw_mode, cb, userdata, &keep_going);
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
    struct ged_draw_pick_result *result;
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
    if (ged_draw_pick_result_append_path(ctx->result, source_path, 0.0))
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

struct ged_draw_pick_result *
ged_draw_view_context_pick_semantic_path(
	struct ged *gedp,
	void *view_ctx,
	const char *path_pattern)
{
    if (!view_ctx || !path_pattern || !path_pattern[0] ||
	    !ged_draw_obol_scene_controller_is_owned(gedp))
	return NULL;

    struct ged_draw_pick_result *result = ged_draw_pick_result_create();
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

    ged_draw_pick_result_free(result);
    return NULL;
}

int
ged_draw_view_context_render_export_consistency(
	struct ged *gedp,
	void *view_ctx,
	const char *drawn_prefix,
	struct ged_draw_render_export_consistency *summary)
{
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

#define GED_DRAW_OBOL_CONTEXT_TOKEN_MAGIC 0x47444f626f6c4354ULL

struct ged_draw_obol_context_token {
    uint64_t magic;
    uint64_t token;
    struct ged *gedp;
    char *path;
    char *instance_key;
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
ged_draw_obol_context_from_scene_handle(ged_draw_scene_handle scene_ref)
{
    if (ged_draw_scene_handle_backend(scene_ref) != GED_DRAW_SCENE_BACKEND_OBOL)
	return NULL;
    return ged_draw_obol_context_from_scene_ctx(
	    ged_draw_scene_handle_context(scene_ref));
}


static ged_draw_group_ref ged_draw_obol_group_ref_for_path(struct ged *gedp,
							   const char *path);


static ged_draw_group_ref
ged_draw_obol_context_group_ref(struct ged *gedp,
				struct ged_draw_obol_context_token *token)
{
    if (!gedp || !token || !token->is_group || token->is_database_source)
	return GED_DRAW_GROUP_REF_NULL;

    ged_draw_scene_handle group_scene_ref =
	ged_draw_scene_handle_make((void *)token, GED_DRAW_SCENE_BACKEND_OBOL);
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
    if (token->instance_key)
	bu_free(token->instance_key, "GED Obol context token instance key");
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
    token->instance_key = info->instance_key ?
			  bu_strdup(info->instance_key) : NULL;
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


static int
ged_draw_obol_context_token_key_equal(const char *a, const char *b)
{
    if ((!a || !a[0]) && (!b || !b[0]))
	return 1;
    if (!a || !b)
	return 0;
    return BU_STR_EQUAL(a, b);
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
		!ged_draw_obol_context_token_key_equal(token->instance_key,
		    info->instance_key) ||
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
ged_draw_obol_context_token_for_path_mode_instance(
	struct ged *gedp,
	const char *path,
	void *parent_ctx,
	int draw_mode_valid,
	int draw_mode,
	const char *instance_key)
{
    if (!gedp || !path || !path[0])
	return NULL;

    struct ged_draw_obol_scene_context_info info;
    memset(&info, 0, sizeof(info));
    if (!ged_draw_obol_scene_context_info_for_path(gedp, path, &info))
	return NULL;
    info.draw_mode_valid = draw_mode_valid ? 1 : 0;
    info.draw_mode = draw_mode;
    if (instance_key && instance_key[0])
	info.instance_key = bu_strdup(instance_key);

    void *ctx = ged_draw_obol_context_token_find(gedp, &info, parent_ctx);
    if (!ctx)
	ctx = ged_draw_obol_context_token_create(gedp, &info, parent_ctx);
    ged_draw_obol_scene_context_info_free(&info);
    return ctx;
}


static void *
ged_draw_obol_context_token_for_path_mode(
	struct ged *gedp,
	const char *path,
	void *parent_ctx,
	int draw_mode_valid,
	int draw_mode)
{
    return ged_draw_obol_context_token_for_path_mode_instance(gedp, path,
	    parent_ctx, draw_mode_valid, draw_mode, NULL);
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
ged_draw_obol_context_token_with_parent_for_path_mode_instance(
	struct ged *gedp,
	const char *path,
	int draw_mode_valid,
	int draw_mode,
	const char *instance_key)
{
    void *parent_ctx =
	ged_draw_obol_database_source_parent_token_for_path(gedp, path);
    return ged_draw_obol_context_token_for_path_mode_instance(gedp, path,
	    parent_ctx, draw_mode_valid, draw_mode, instance_key);
}


static void *
ged_draw_obol_context_token_with_parent_for_path_mode(
	struct ged *gedp,
	const char *path,
	int draw_mode_valid,
	int draw_mode)
{
    return ged_draw_obol_context_token_with_parent_for_path_mode_instance(
	    gedp, path, draw_mode_valid, draw_mode, NULL);
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

    ged_draw_scene_handle scene_ref = ged_draw_scene_handle_make(scene_ctx,
	    GED_DRAW_SCENE_BACKEND_OBOL);
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


static ged_draw_scene_handle
ged_draw_scene_context_handle(void *scene_ctx)
{
    if (ged_draw_obol_context_from_scene_ctx(scene_ctx))
	return ged_draw_scene_handle_make(scene_ctx, GED_DRAW_SCENE_BACKEND_OBOL);
    return ged_draw_scene_handle_null();
}


static struct ged_draw_obol_context_token *
ged_draw_shape_ref_obol_token(struct ged *gedp, ged_draw_shape_ref ref)
{
    ged_draw_scene_handle scene_ref = ged_draw_registry_shape_ref_scene_handle(gedp,
	    ref);
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_handle(scene_ref);
    if (token)
	return token;

    scene_ref = ged_draw_registry_shape_ref_cache_scene_handle(gedp, ref);
    return ged_draw_obol_context_from_scene_handle(scene_ref);
}


static struct ged_draw_obol_context_token *
ged_draw_group_ref_obol_token(struct ged *gedp, ged_draw_group_ref ref)
{
    ged_draw_scene_handle scene_ref = ged_draw_registry_group_ref_scene_handle(gedp,
	    ref);
    return ged_draw_obol_context_from_scene_handle(scene_ref);
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
				   ged_draw_shape_ref_obol_path_cb cb,
				   void *userdata,
				   const char *free_label)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_shape_ref_obol_token(gedp, ref);
    if (token && token->instance_key && token->instance_key[0] &&
	    _ged_draw_shape_ref_try_obol_path(gedp, token->instance_key, cb,
		userdata))
	return 1;

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
    (void)free_label;
    return 0;
}


static int
_ged_draw_shape_ref_try_obol_paths_lazy(struct ged *gedp,
					ged_draw_shape_ref ref,
					ged_draw_shape_ref_obol_path_cb cb,
					void *userdata,
					const char *free_label)
{
    return _ged_draw_shape_ref_try_obol_paths(gedp, ref, cb, userdata,
	    free_label);
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

    (void)_ged_draw_shape_ref_try_obol_paths(gedp, ref,
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
    (void)free_label;
    return 0;
}


static int
_ged_draw_group_ref_try_obol_paths_ensure(
	struct ged *gedp,
	ged_draw_group_ref ref,
	ged_draw_group_ref_obol_path_cb cb,
	void *userdata,
	const char *free_label,
	const char *preferred_path,
	int mode,
	int overlay)
{

    if (_ged_draw_group_ref_try_obol_paths(gedp, ref, cb,
	    userdata, free_label))
	return 1;

    if (!preferred_path || !preferred_path[0] ||
	    !ged_draw_obol_group_ensure_for_path(gedp, preferred_path,
		preferred_path, mode, overlay))
	return 0;

    return _ged_draw_group_ref_try_obol_path(gedp, preferred_path, cb,
	    userdata);
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

    (void)_ged_draw_group_ref_try_obol_paths(gedp, ref,
	    ged_draw_group_ref_obol_context_cb, &ctx,
	    "GED Obol group-ref context path");
    return ctx.context;
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
    int draw_mode_valid;
    int draw_mode;
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
    if (ged_draw_obol_database_source_geometry_summary_for_path_mode(gedp,
	    path, ctx->draw_mode_valid, ctx->draw_mode, ctx->out))
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
    int updated = ged_draw_obol_database_source_update_display_for_path(
	    gedp, path,
	    ctx->visible_valid, ctx->visible,
	    ctx->selected_valid, ctx->selected,
	    ctx->highlighted_valid, ctx->highlighted,
	    ctx->draw_mode_valid, ctx->draw_mode,
	    ctx->line_style_valid, ctx->line_style,
	    ctx->line_width_valid, ctx->line_width,
	    ctx->transparency_valid, ctx->transparency,
	    ctx->color_valid, ctx->color,
	    ctx->material_color_valid, ctx->material_color,
	    ctx->material_revision_valid, ctx->material_revision);
    if (updated) {
	struct ged_draw_obol_database_source_record record;
	memset(&record, 0, sizeof(record));
	if (ged_draw_obol_database_source_record_for_path(gedp, path,
		&record) && record.database_path && record.database_path[0] &&
		!BU_STR_EQUAL(record.database_path, path)) {
	    (void)ged_draw_obol_database_source_update_display_for_path(
		gedp, record.database_path,
		ctx->visible_valid, ctx->visible,
		ctx->selected_valid, ctx->selected,
		ctx->highlighted_valid, ctx->highlighted,
		ctx->draw_mode_valid, ctx->draw_mode,
		ctx->line_style_valid, ctx->line_style,
		ctx->line_width_valid, ctx->line_width,
		ctx->transparency_valid, ctx->transparency,
		ctx->color_valid, ctx->color,
		ctx->material_color_valid, ctx->material_color,
		ctx->material_revision_valid, ctx->material_revision);
	}
	return 1;
    }
    return ged_draw_obol_shape_update_display_for_path(gedp, path,
	    ctx->visible_valid, ctx->visible,
	    ctx->selected_valid, ctx->selected,
	    ctx->highlighted_valid, ctx->highlighted,
	    ctx->draw_mode_valid, ctx->draw_mode,
	    ctx->line_style_valid, ctx->line_style,
	    ctx->line_width_valid, ctx->line_width,
	    ctx->transparency_valid, ctx->transparency,
	    ctx->color_valid, ctx->color,
	    ctx->material_color_valid, ctx->material_color,
	    ctx->material_revision_valid, ctx->material_revision);
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


static ged_draw_scene_handle
_ged_draw_source_root_scene_handle(struct ged *gedp)
{
    if (!gedp)
	return ged_draw_scene_handle_null();

    return ged_draw_registry_group_ref_scene_handle(gedp,
	    ged_scene_root_group_ref(gedp));
}


struct _ged_draw_source_rt_shape_iter_ctx {
    int has_record_ancestor;
    ged_draw_scene_handle_child_cb cb;
    void *userdata;
};


static int ged_draw_source_scene_handle_foreach_shape_impl(
	ged_draw_scene_handle scene_ref,
	int has_record_ancestor,
	ged_draw_scene_handle_child_cb cb,
	void *userdata);


static int
ged_draw_source_scene_handle_foreach_shape_child_cb(
	ged_draw_scene_handle child_ref,
	void *userdata)
{
    struct _ged_draw_source_rt_shape_iter_ctx *ctx =
	(struct _ged_draw_source_rt_shape_iter_ctx *)userdata;
    return ctx ? ged_draw_source_scene_handle_foreach_shape_impl(child_ref,
	    ctx->has_record_ancestor, ctx->cb, ctx->userdata) : 1;
}


static int
ged_draw_source_scene_handle_foreach_shape_impl(
	ged_draw_scene_handle scene_ref,
	int has_record_ancestor,
	ged_draw_scene_handle_child_cb cb,
	void *userdata)
{
    if (ged_draw_scene_handle_is_null(scene_ref))
	return 1;

    struct ged_draw_scene_tree_summary tree_summary;
    memset(&tree_summary, 0, sizeof(tree_summary));
    int has_record = ged_draw_scene_handle_tree_summary(scene_ref,
	    &tree_summary) && tree_summary.is_shape;
    if (has_record) {
	int has_path = tree_summary.fullpath &&
	    tree_summary.fullpath->fp_len > 0;
	if ((has_path || !has_record_ancestor) && !(*cb)(scene_ref, userdata))
	    return 0;
    }

    struct _ged_draw_source_rt_shape_iter_ctx args;
    args.has_record_ancestor = has_record_ancestor || has_record;
    args.cb = cb;
    args.userdata = userdata;
    return ged_draw_scene_handle_foreach_child(scene_ref,
	    ged_draw_source_scene_handle_foreach_shape_child_cb, &args);
}


static int
ged_draw_source_scene_handle_foreach_shape(
	ged_draw_scene_handle scene_ref,
	ged_draw_scene_handle_child_cb cb,
	void *userdata)
{
    if (!cb)
	return 1;
    return ged_draw_source_scene_handle_foreach_shape_impl(scene_ref, 0,
	    cb, userdata);
}


struct ged_draw_source_root_shape_ref_ctx {
    struct ged *gedp;
    ged_draw_shape_ref_index_cb cb;
    void *userdata;
};


static int
_ged_draw_source_root_shape_ref_cb(ged_draw_scene_handle scene_ref,
				   void *userdata)
{
    struct ged_draw_source_root_shape_ref_ctx *ctx =
	(struct ged_draw_source_root_shape_ref_ctx *)userdata;
    if (!ctx || !ctx->cb)
	return 1;

    ged_draw_shape_ref ref = ged_draw_registry_shape_ref_from_source_ref(
	    ctx->gedp, scene_ref);
    if (ged_draw_shape_ref_is_null(ref))
	return 1;
    const struct db_full_path *fullpath =
	ged_draw_scene_handle_fullpath(scene_ref);
    if (fullpath && fullpath->fp_len > 0)
	(void)ged_draw_registry_shape_ref_set_indexed_fullpath(ctx->gedp,
		ref, fullpath);
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

    ged_draw_scene_handle scene_ref = ged_draw_scene_handle_make(scene_ctx,
	    GED_DRAW_SCENE_BACKEND_OBOL);
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

    void *scene_ctx =
	ged_draw_obol_context_token_with_parent_for_path_mode_instance(gedp,
	    record->database_path, record->valid, record->draw_mode,
	    record->instance_key);
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_ctx(scene_ctx);
    if (!token || !token->is_database_source)
	return GED_DRAW_SHAPE_REF_NULL;

    ged_draw_scene_handle scene_ref = ged_draw_scene_handle_make(scene_ctx,
	    GED_DRAW_SCENE_BACKEND_OBOL);
    ged_draw_shape_ref ref =
	ged_draw_registry_shape_ref_from_source_ref(gedp, scene_ref);
    if (ged_draw_shape_ref_is_null(ref))
	return GED_DRAW_SHAPE_REF_NULL;

    if (token->fullpath_valid)
	(void)ged_draw_registry_shape_ref_set_indexed_fullpath(gedp, ref,
		&token->fullpath);

    return ref;
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
_ged_draw_source_root_shape_child_cb(ged_draw_scene_handle child_ref,
				     void *userdata)
{
    struct ged_draw_source_root_shape_child_ctx *ctx =
	(struct ged_draw_source_root_shape_child_ctx *)userdata;
    if (!ctx)
	return 1;

    if (ctx->skip_overlay_groups) {
	struct ged_draw_group_record_summary group_summary;
	struct ged_draw_obol_context_token *token =
	    ged_draw_obol_context_from_scene_handle(child_ref);
	if (token && token->gedp && token->path &&
		ged_draw_obol_group_record_summary_for_path(token->gedp,
		    token->path, &group_summary) && group_summary.is_overlay)
	    return 1;
    }

    return ged_draw_source_scene_handle_foreach_shape(child_ref,
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

    ged_draw_scene_handle root_ref = _ged_draw_source_root_scene_handle(gedp);
    if (ged_draw_scene_handle_is_null(root_ref))
	return 1;

    struct ged_draw_source_root_shape_child_ctx ctx;
    ctx.shape_ctx.gedp = gedp;
    ctx.shape_ctx.cb = cb;
    ctx.shape_ctx.userdata = userdata;
    ctx.skip_overlay_groups = skip_overlay_groups;
    return ged_draw_scene_handle_foreach_child(root_ref,
	    _ged_draw_source_root_shape_child_cb, &ctx);
}


struct ged_draw_visible_shape_candidate_ctx {
    ged_draw_shape_candidate_cb cb;
    void *userdata;
};


static int
_ged_draw_visible_shape_candidate_cb(
	struct ged *UNUSED(gedp),
	const struct ged_draw_obol_database_source_record *record,
	void *userdata)
{
    struct ged_draw_visible_shape_candidate_ctx *ctx =
	(struct ged_draw_visible_shape_candidate_ctx *)userdata;
    if (!ctx || !ctx->cb || !record || !record->database_path)
	return 0;

    struct ged_draw_shape_candidate candidate;
    candidate.path = record->database_path;
    candidate.instance_key = record->instance_key;
    candidate.draw_mode = record->draw_mode;
    return ctx->cb(&candidate, ctx->userdata);
}


void
ged_draw_foreach_visible_shape_candidate(struct ged *gedp,
					 ged_draw_shape_candidate_cb cb,
					 void *userdata)
{
    if (!gedp || !cb || !ged_draw_obol_scene_controller_full_synced(gedp))
	return;

    struct ged_draw_visible_shape_candidate_ctx ctx;
    ctx.cb = cb;
    ctx.userdata = userdata;
    (void)ged_draw_obol_visible_database_source_records_foreach_fast(gedp,
	    _ged_draw_visible_shape_candidate_cb, &ctx);
}


ged_draw_shape_ref
ged_draw_shape_ref_for_candidate(
    struct ged *gedp,
    const struct ged_draw_shape_candidate *candidate)
{
    if (!gedp || !candidate || !candidate->path || !candidate->path[0])
	return GED_DRAW_SHAPE_REF_NULL;

    struct ged_draw_obol_database_source_record record;
    memset(&record, 0, sizeof(record));
    record.valid = 1;
    record.database_path = candidate->path;
    record.instance_key = candidate->instance_key;
    record.draw_mode = candidate->draw_mode;
    return ged_draw_obol_shape_ref_for_database_source_record(gedp, &record);
}


struct ged_draw_source_root_group_ref_ctx {
    struct ged *gedp;
    ged_draw_group_ref_index_cb cb;
    void *userdata;
};


static int _ged_draw_source_root_foreach_group_scene(
	struct ged_draw_source_root_group_ref_ctx *ctx,
	ged_draw_scene_handle scene_ref);


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
_ged_draw_source_root_group_child_cb(ged_draw_scene_handle child_ref,
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
	ged_draw_scene_handle scene_ref)
{
    if (!ctx || ged_draw_scene_handle_is_null(scene_ref))
	return 1;

    struct ged_draw_scene_tree_summary tree_summary;
    memset(&tree_summary, 0, sizeof(tree_summary));
    if (ged_draw_scene_handle_tree_summary(scene_ref, &tree_summary) &&
	    tree_summary.is_group) {
	ged_draw_group_ref ref =
	    ged_draw_registry_group_ref_from_source_ref(ctx->gedp, scene_ref);
	if (!ged_draw_group_ref_is_null(ref) && !ctx->cb(ref, ctx->userdata))
	    return 0;
    }

    return ged_draw_scene_handle_foreach_child(scene_ref,
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

    ged_draw_scene_handle root_ref = _ged_draw_source_root_scene_handle(gedp);
    if (ged_draw_scene_handle_is_null(root_ref))
	return 1;

    struct ged_draw_source_root_group_ref_ctx ctx;
    ctx.gedp = gedp;
    ctx.cb = cb;
    ctx.userdata = userdata;
    return ged_draw_scene_handle_foreach_child(root_ref,
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

    ged_draw_scene_handle root_ref = _ged_draw_source_root_scene_handle(gedp);
    if (ged_draw_scene_handle_is_null(root_ref))
	return 1;

    return ged_draw_scene_handle_subtree_bounds(root_ref, min, max,
	    include_overlays);
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


int
ged_draw_source_root_clear_all_scope_children(struct ged *gedp)
{
    if (!gedp)
	return 0;

    int obol_cleared = ged_draw_obol_active_scene_clear(gedp);
    if (obol_cleared) {
	ged_draw_obol_owner_structural_revision_bump(gedp);
	return 1;
    }
    return 0;
}


static int
ged_draw_scene_handle_foreach_child(ged_draw_scene_handle scene_ref,
				    ged_draw_scene_handle_child_cb cb,
				    void *userdata)
{
    if (ged_draw_scene_handle_is_null(scene_ref) || !cb)
	return 1;

    struct ged_draw_scene_tree_summary summary;
    if (!ged_draw_scene_handle_tree_summary(scene_ref, &summary))
	return 1;

    void *parent_ctx = ged_draw_scene_handle_context(scene_ref);
    for (size_t i = 0; i < summary.child_count; i++) {
	void *child_ctx = ged_draw_obol_scene_context_child_at(scene_ref,
		parent_ctx, i);
	if (!child_ctx)
	    continue;
	ged_draw_scene_handle child_ref = ged_draw_scene_handle_make(child_ctx,
		GED_DRAW_SCENE_BACKEND_OBOL);
	if (!(*cb)(child_ref, userdata))
	    return 0;
    }
    return 1;
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
	ged_draw_obol_owner_structural_revision_bump(gedp);
	return 1;
    }
    return 0;
}


int
ged_draw_source_clear_db_groups_in_scope(struct ged *gedp, void *view_ctx)
{
    if (!gedp)
	return 0;

    int obol_cleared = ged_draw_obol_database_sources_clear_in_scope(gedp,
	    view_ctx);
    if (obol_cleared) {
	ged_draw_obol_owner_structural_revision_bump(gedp);
	return 1;
    }
    return 0;
}


static int
ged_draw_obol_erase_exact_path_owner(struct ged *gedp, const char *path,
				     void *view_ctx,
				     int mode)
{
    if (!gedp || !path || !path[0])
	return 0;

    int independent_scope = view_ctx && ged_view_context_is_independent(view_ctx);
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
    if (view_ctx && ged_view_context_is_independent(view_ctx)) {
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

    return ged_draw_obol_erase_exact_path_owner(gedp, path, NULL, -1);
}


int
ged_draw_source_erase_path_prefix_at_root(struct ged *gedp, const char *path)
{
    if (!gedp || !path)
	return 0;

    return ged_draw_obol_erase_path_prefix_owner(gedp, path, NULL, -1);
}


int
ged_draw_source_erase_path_in_active_scope(struct ged *gedp,
					   const char *path,
					   void *view_ctx,
					   int mode)
{
    if (!gedp || !path)
	return 0;

    return ged_draw_obol_erase_exact_path_owner(gedp, path, view_ctx, mode);
}


int
ged_draw_source_erase_path_prefix_in_active_scope(struct ged *gedp,
						  const char *path,
						  void *view_ctx,
						  int mode)
{
    if (!gedp || !path)
	return 0;

    return ged_draw_obol_erase_path_prefix_owner(gedp, path, view_ctx, mode);
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

    if (view_ctx && ged_view_context_is_independent(view_ctx))
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
    return obol_erased ? 1 : 0;
}

int
ged_draw_scene_context_display_summary(void *scene_ctx,
				       struct ged_draw_scene_display_summary *out)
{
    return ged_draw_scene_handle_display_summary(
	    ged_draw_scene_context_handle(scene_ctx), out);
}


static int
ged_draw_scene_handle_display_summary(
	ged_draw_scene_handle scene_ref,
	struct ged_draw_scene_display_summary *out)
{
    struct ged *gedp = ged_draw_scene_handle_owner_gedp(scene_ref);
    const char *path = ged_draw_scene_handle_semantic_path(scene_ref);
    if (gedp && path) {
	struct ged_draw_obol_context_token *token =
	    ged_draw_obol_context_from_scene_handle(scene_ref);
	int is_group = token && token->is_group &&
	    !token->is_database_source;
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
    if (out)
	memset(out, 0, sizeof(*out));
    return 0;
}


int
ged_draw_shape_context_has_state(void *shape_ctx)
{
    if (!shape_ctx)
	return 0;

    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_ctx(shape_ctx);
    return token && (token->is_database_source || token->is_shape);
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
    return NULL;
}


void *
ged_draw_first_shape_context(struct ged *gedp)
{
    return ged_draw_shape_ref_context(gedp, ged_draw_first_shape_ref(gedp));
}


void *
ged_draw_shape_ref_context(struct ged *gedp, ged_draw_shape_ref ref)
{
    return _ged_draw_shape_ref_obol_context(gedp, ref);
}


void *
ged_draw_shape_ref_cache_context(struct ged *gedp, ged_draw_shape_ref ref)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_handle(
	    ged_draw_registry_shape_ref_cache_scene_handle(gedp, ref));
    return token ? (void *)token : NULL;
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
	    out->selected = display_summary.selected;
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
	    out->selected = display_summary.selected;
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

    return 0;
}


int
ged_draw_shape_ref_display_summary(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   struct ged_draw_scene_display_summary *out)
{
    struct _ged_draw_obol_display_summary_ctx ctx = {out};
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref,
	    _ged_draw_obol_display_summary_cb, &ctx,
	    "GED Obol display summary path"))
	return 1;

    if (out)
	memset(out, 0, sizeof(*out));
    return 0;
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
    struct _ged_draw_obol_source_summary_ctx ctx = {out};
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref,
	    _ged_draw_obol_source_summary_cb, &ctx,
	    "GED Obol database source summary path"))
	return 1;

    if (out)
	memset(out, 0, sizeof(*out));
    return 0;
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
    return ged_draw_group_ref_obol_context(gedp, ref);
}


int
ged_draw_group_ref_record_summary(struct ged *gedp,
				  ged_draw_group_ref ref,
				  struct ged_draw_group_record_summary *out)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_group_ref_obol_token(gedp, ref);
    if ((token || ged_draw_group_ref_is_null(ref)) && out) {
	memset(out, 0, sizeof(*out));
	struct _ged_draw_obol_group_record_ctx ctx = {out};
	if (_ged_draw_group_ref_try_obol_paths(gedp, ref,
		_ged_draw_obol_group_record_cb, &ctx,
		"GED Obol group record path"))
	    return 1;
    }

    if (out)
	memset(out, 0, sizeof(*out));
    return 0;
}


int
ged_draw_group_ref_tree_summary(struct ged *gedp,
				ged_draw_group_ref ref,
				struct ged_draw_scene_tree_summary *out)
{
    if (!out)
	return 0;
    memset(out, 0, sizeof(*out));

    struct _ged_draw_obol_group_tree_summary_ctx ctx = {
	out,
	NULL
    };
    if (_ged_draw_group_ref_try_obol_paths(gedp, ref,
	    _ged_draw_obol_group_tree_summary_cb, &ctx,
	    "GED Obol group tree-summary path"))
	return 1;

    return 0;
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

    struct _ged_draw_obol_group_shape_count_ctx ctx = {out};
    if (_ged_draw_group_ref_try_obol_paths(gedp, ref,
	    _ged_draw_obol_group_shape_count_cb, &ctx,
	    "GED Obol group shape-count path"))
	return 1;

    return 0;
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
    if (!ged_draw_obol_scene_controller_ensure_owned(gedp, 1)) {
	bu_free(s, "draw group path string");
	return GED_DRAW_GROUP_REF_NULL;
    }

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
    /* Obol owns public group DB-path mutation. */
    return 0;
}


int
ged_draw_group_ref_set_mode(struct ged *gedp,
			    ged_draw_group_ref ref,
			    int mode)
{
    struct _ged_draw_obol_group_draw_intent_ctx ctx = {
	NULL,
	1,
	mode,
	0,
	0
    };
    return _ged_draw_group_ref_try_obol_paths_ensure(gedp, ref,
	    _ged_draw_obol_group_draw_intent_cb, &ctx,
	    "GED Obol group mode path", NULL, mode, 0);
}


int
ged_draw_group_ref_set_appearance_settings(struct ged *gedp,
					   ged_draw_group_ref ref,
					   const struct ged_draw_appearance_settings *settings)
{
    if (!settings)
	return 0;

    struct _ged_draw_obol_group_appearance_ctx ctx = {settings};
    return _ged_draw_group_ref_try_obol_paths_ensure(gedp, ref,
	    _ged_draw_obol_group_update_appearance_cb, &ctx,
	    "GED Obol group appearance path", NULL, settings->draw_mode, 0);
}


int
ged_draw_group_ref_appearance_settings(struct ged *gedp,
				       ged_draw_group_ref ref,
				       struct ged_draw_appearance_settings *settings)
{
    if (!settings)
	return 0;

    struct _ged_draw_obol_group_appearance_read_ctx ctx = {settings};
    if (_ged_draw_group_ref_try_obol_paths(gedp, ref,
	    _ged_draw_obol_group_appearance_read_cb, &ctx,
	    "GED Obol group appearance read path"))
	return 1;

    return 0;
}


ged_draw_shape_ref
ged_draw_shape_ref_from_context(struct ged *gedp, void *shape_ctx)
{
    return ged_draw_registry_shape_ref_from_source_ref(gedp,
	    ged_draw_scene_context_handle(shape_ctx));
}


static struct ged *
ged_draw_scene_handle_owner_gedp(ged_draw_scene_handle scene_ref)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_handle(scene_ref);
    if (token)
	return token->gedp;
    return NULL;
}


static int
ged_draw_scene_handle_is_source_root(struct ged *gedp,
				     ged_draw_scene_handle scene_ref)
{
    if (!gedp || ged_draw_scene_handle_is_null(scene_ref))
	return 0;

    ged_draw_scene_handle root_ref = _ged_draw_source_root_scene_handle(gedp);
    if (!ged_draw_scene_handle_is_null(root_ref) &&
	    ged_draw_scene_handle_equal(scene_ref, root_ref))
	return 1;

    return 0;
}


static const char *
ged_draw_scene_handle_semantic_path(ged_draw_scene_handle scene_ref)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_handle(scene_ref);
    if (token)
	return (token->path && token->path[0]) ? token->path : NULL;
    return NULL;
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

    ged_draw_scene_handle scene_ref = ged_draw_scene_context_handle(shape_ctx);
    struct ged *gedp = ged_draw_scene_handle_owner_gedp(scene_ref);
    if (!gedp)
	return 0;

    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_handle(scene_ref);
    if (token && token->instance_key && token->instance_key[0] &&
	    _ged_draw_shape_ref_try_obol_path(gedp, token->instance_key, cb,
		userdata))
	return 1;

    (void)free_label;
    return _ged_draw_shape_ref_try_obol_path(gedp,
	    ged_draw_scene_handle_semantic_path(scene_ref), cb, userdata);
}


const struct db_full_path *
ged_draw_scene_context_fullpath(void *scene_ctx)
{
    return ged_draw_scene_handle_fullpath(
	    ged_draw_scene_context_handle(scene_ctx));
}


static const struct db_full_path *
ged_draw_scene_handle_fullpath(ged_draw_scene_handle scene_ref)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_handle(scene_ref);
    if (token)
	return token->fullpath_valid ? &token->fullpath : NULL;
    return NULL;
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

    if (!ged_draw_obol_context_from_scene_ctx(group_ctx))
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

    ged_draw_scene_handle scene_ref = ged_draw_scene_context_handle(group_ctx);
    struct ged *gedp = ged_draw_scene_handle_owner_gedp(scene_ref);
    const char *path = ged_draw_scene_handle_semantic_path(scene_ref);
    if (gedp && path) {
	struct ged_draw_group_record_summary obol_summary;
	memset(&obol_summary, 0, sizeof(obol_summary));
	if (ged_draw_obol_group_record_summary_for_path(gedp, path,
		&obol_summary))
	    return obol_summary.is_overlay;
    }
    return 0;
}


int
ged_draw_scene_context_source_summary(void *scene_ctx,
				      struct ged_draw_database_source_summary *out)
{
    return ged_draw_scene_handle_source_summary(
	    ged_draw_scene_context_handle(scene_ctx), out);
}


static int
ged_draw_scene_handle_source_summary(
	ged_draw_scene_handle scene_ref,
	struct ged_draw_database_source_summary *out)
{
    struct ged *gedp = ged_draw_scene_handle_owner_gedp(scene_ref);
    const char *path = ged_draw_scene_handle_semantic_path(scene_ref);
    if (gedp && path &&
	    ged_draw_obol_database_source_summary_for_path(gedp, path, out))
	return 1;
    if (out)
	memset(out, 0, sizeof(*out));
    return 0;
}


int
ged_draw_scene_context_tree_summary(void *scene_ctx,
				    struct ged_draw_scene_tree_summary *out)
{
    return ged_draw_scene_handle_tree_summary(
	    ged_draw_scene_context_handle(scene_ctx), out);
}


static int
ged_draw_scene_handle_tree_summary(
	ged_draw_scene_handle scene_ref,
	struct ged_draw_scene_tree_summary *out)
{
    if (!out)
	return 0;

    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_handle(scene_ref);
    if (token)
	return ged_draw_obol_context_token_tree_summary(token, NULL, out);
    memset(out, 0, sizeof(*out));
    return 0;
}


static void *
ged_draw_obol_scene_context_child_at(
	ged_draw_scene_handle scene_ref,
	void *parent_ctx,
	size_t index)
{
    struct ged *gedp = ged_draw_scene_handle_owner_gedp(scene_ref);
    struct ged_draw_obol_scene_context_info info;
    memset(&info, 0, sizeof(info));

    if (gedp && ged_draw_scene_handle_is_source_root(gedp, scene_ref) &&
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

    const char *path = ged_draw_scene_handle_semantic_path(scene_ref);
    if (!gedp || !path)
	return NULL;

    size_t child_index = index;
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_handle(scene_ref);
    int source_context = token && token->is_database_source;
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
		ged_draw_scene_handle_tree_summary(scene_ref,
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
    ged_draw_scene_handle scene_ref = ged_draw_scene_context_handle(scene_ctx);
    void *obol_child = ged_draw_obol_scene_context_child_at(scene_ref,
	    scene_ctx, index);
    return obol_child;
}


void *
ged_draw_scene_context_parent(void *scene_ctx)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_ctx(scene_ctx);
    if (token)
	return token->parent_ctx;
    return NULL;
}


const char *
ged_draw_scene_context_name(void *scene_ctx)
{
    return ged_draw_scene_handle_name(ged_draw_scene_context_handle(scene_ctx));
}


static const char *
ged_draw_scene_handle_name(ged_draw_scene_handle scene_ref)
{
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_handle(scene_ref);
    if (token)
	return token->name ? token->name : token->path;
    return NULL;
}


int
ged_draw_scene_context_subtree_bounds(void *scene_ctx,
				      vect_t *min,
				      vect_t *max,
				      int include_overlays)
{
    return ged_draw_scene_handle_subtree_bounds(
	    ged_draw_scene_context_handle(scene_ctx), min, max,
	    include_overlays);
}


static int
ged_draw_scene_handle_subtree_bounds(ged_draw_scene_handle scene_ref,
				     vect_t *min,
				     vect_t *max,
				     int include_overlays)
{
    if (!min || !max)
	return 1;

    struct ged *gedp = ged_draw_scene_handle_owner_gedp(scene_ref);
    const char *path = ged_draw_scene_handle_semantic_path(scene_ref);
    if (gedp && path) {
	struct ged_draw_obol_context_token *token =
	    ged_draw_obol_context_from_scene_handle(scene_ref);
	int empty = 1;
	int is_group_context = token && token->is_group &&
	    !token->is_database_source;
	if (is_group_context) {
	    if (ged_draw_obol_group_subtree_bounds_for_path(gedp, path,
		    min, max, include_overlays, &empty))
		return empty;
	} else if (ged_draw_obol_database_source_bounds_for_path(gedp,
		path, min, max, include_overlays, &empty)) {
	    return empty;
	}
    }
    return 1;
}


int
ged_draw_shape_ref_material_summary(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    struct ged_draw_shape_material_summary *out)
{
    if (!out)
	return 0;

    struct _ged_draw_obol_material_summary_ctx ctx = {out};
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref,
	    _ged_draw_obol_material_summary_cb, &ctx,
	    "GED Obol material summary path"))
	return 1;

    memset(out, 0, sizeof(*out));
    return 0;
}


static int
ged_draw_shape_ref_obol_clear_vlist(struct ged *gedp,
				    ged_draw_shape_ref ref)
{
    return _ged_draw_shape_ref_try_obol_paths(gedp, ref,
	    _ged_draw_obol_clear_vlist_cb, NULL,
	    "GED Obol clear VLIST path");
}


static int
ged_draw_shape_ref_obol_clear_mesh(struct ged *gedp,
				   ged_draw_shape_ref ref)
{
    return _ged_draw_shape_ref_try_obol_paths(gedp, ref,
	    _ged_draw_obol_clear_mesh_cb, NULL,
	    "GED Obol clear mesh path");
}


static int
ged_draw_shape_ref_obol_clear_auxiliary_shapes(struct ged *gedp,
					       ged_draw_shape_ref ref)
{
    return _ged_draw_shape_ref_try_obol_paths(gedp, ref,
	    _ged_draw_obol_clear_auxiliary_shapes_cb, NULL,
	    "GED Obol clear auxiliary VLIST path");
}


static int
ged_draw_shape_ref_obol_update_vlist_bounds(
	struct ged *gedp,
	ged_draw_shape_ref ref)
{
    return _ged_draw_shape_ref_try_obol_paths(gedp, ref,
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
    /* Obol owns public shape center mutation. */
    return _ged_draw_shape_ref_try_obol_paths_lazy(gedp, ref,
	    _ged_draw_obol_set_center_cb, &ctx,
	    "GED Obol set VLIST center path");
}


int
ged_draw_shape_ref_geometry_clear(struct ged *gedp, ged_draw_shape_ref ref)
{
    int obol_cleared = 0;
    obol_cleared = ged_draw_shape_ref_obol_clear_vlist(gedp, ref);
    obol_cleared =
	ged_draw_shape_ref_obol_clear_mesh(gedp, ref) ||
	obol_cleared;
    obol_cleared =
	ged_draw_shape_ref_obol_clear_auxiliary_shapes(gedp, ref) ||
	obol_cleared;
    return obol_cleared ? 1 : 0;
}


int
ged_draw_shape_ref_update_bounds_from_geometry(struct ged *gedp,
					       ged_draw_shape_ref ref,
					       int *bad_cmd)
{
    if (bad_cmd)
	*bad_cmd = 0;
    return ged_draw_shape_ref_obol_update_vlist_bounds(gedp, ref);
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

    return 0;
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


int
ged_draw_view_feature_primitive_wireframe_replace(
    ged_draw_view_feature_ref ref,
    struct db_i *dbip,
    struct rt_db_internal *ip,
    const mat_t mat,
    const struct bg_tess_tol *ttol,
    const struct bn_tol *tol)
{
    struct bu_list vhead;
    struct rt_db_internal local_ip;
    struct rt_db_internal *plot_ip = ip;
    point_t *points = NULL;
    int *commands = NULL;
    size_t point_count = 0;
    int have_local_ip = 0;
    int ret = 0;

    RT_DB_INTERNAL_INIT(&local_ip);
    BU_LIST_INIT(&vhead);

    if (ged_draw_view_feature_ref_is_null(ref) || !ip)
	return 0;

    RT_CK_DB_INTERNAL(ip);
    if (!ged_draw_rt_internal_payload_valid(ip))
	goto cleanup;

    if (mat) {
	if (rt_matrix_transform(&local_ip, mat, ip, 0, dbip) < 0)
	    goto cleanup;
	have_local_ip = 1;
	plot_ip = &local_ip;
	if (!ged_draw_rt_internal_payload_valid(plot_ip))
	    goto cleanup;
    }

    if (rt_obj_plot(&vhead, plot_ip, ttol, tol) < 0 ||
	    !ged_draw_rt_vlist_to_ged_line_arrays(&vhead, &points,
		&commands, &point_count))
	goto cleanup;

    ret = ged_draw_view_feature_points_replace(ref,
	    GED_DRAW_VIEW_FEATURE_TRANSIENT_PREVIEW,
	    (const point_t *)points, commands, point_count);

cleanup:
    if (points)
	bu_free(points, "GED edit-preview RT VLIST line points");
    if (commands)
	bu_free(commands, "GED edit-preview RT VLIST line commands");
    RT_FREE_VLIST(&rt_vlfree, &vhead);
    if (have_local_ip)
	rt_db_free_internal(&local_ip);
    return ret;
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


int
ged_draw_shape_ref_redraw_wireframe(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    struct db_i *dbip,
				    const struct bn_tol *tol,
				    const struct bg_tess_tol *ttol,
				    void *view_ctx,
				    int skip_subtractions)
{
    return ged_draw_shape_ref_redraw_wireframe_obol(gedp, ref, dbip,
	    tol, ttol, view_ctx, skip_subtractions);
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
    return -1;
}


static int
ged_draw_shape_ref_obol_update_display_lazy(
	struct ged *gedp,
	ged_draw_shape_ref ref,
	int visible_valid,
	int visible,
	int selected_valid,
	int selected,
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
	selected_valid,
	selected,
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
    return _ged_draw_shape_ref_try_obol_paths_lazy(gedp, ref,
	    _ged_draw_obol_update_display_cb, &ctx,
	    "GED Obol display mutation path");
}


int
ged_draw_shape_ref_set_visible(struct ged *gedp, ged_draw_shape_ref ref,
			       int visible)
{
    int obol_updated = ged_draw_shape_ref_obol_update_display_lazy(gedp,
	    ref, 1, visible, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0.0, 0, NULL, 0, NULL, 0, 0);
    /* Obol owns public shape visibility mutation. */
    return obol_updated;
}


int
ged_draw_group_ref_set_visible(struct ged *gedp, ged_draw_group_ref ref,
			       int visible)
{
    struct _ged_draw_obol_group_visible_ctx ctx = {1, visible};
    return _ged_draw_group_ref_try_obol_paths_ensure(gedp, ref,
	    _ged_draw_obol_group_visible_cb, &ctx,
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
	    ref, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0.0, 1, rgb, 0, NULL, 0, 0);
    /* Obol owns public shape color mutation. */
    return obol_updated;
}


int
ged_draw_shape_ref_set_highlighted(struct ged *gedp, ged_draw_shape_ref ref,
				   int highlighted)
{
    int obol_updated = ged_draw_shape_ref_obol_update_display_lazy(gedp,
	    ref, 0, 0, 0, 0, 1, highlighted, 0, 0, 0, 0,
	    0, 0, 0, 0.0, 0, NULL, 0, NULL, 0, 0);
    /* Obol owns public shape highlight mutation. */
    return obol_updated;
}


int
ged_draw_shape_ref_set_selected(struct ged *gedp, ged_draw_shape_ref ref,
				int selected)
{
    int obol_updated = ged_draw_shape_ref_obol_update_display_lazy(gedp,
	    ref, 0, 0, 1, selected, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0.0, 0, NULL, 0, NULL, 0, 0);
    /* Obol owns public shape selection mutation. */
    return obol_updated;
}


int
ged_draw_shape_ref_set_transparency(struct ged *gedp, ged_draw_shape_ref ref,
				    fastf_t transparency)
{
    int obol_updated = ged_draw_shape_ref_obol_update_display_lazy(gedp,
	    ref, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	    1, transparency, 0, NULL, 0, NULL, 0, 0);
    /* Obol owns public shape transparency mutation. */
    return obol_updated;
}


int
ged_draw_shape_ref_mark_database_source_changed(struct ged *gedp,
						ged_draw_shape_ref ref,
						ged_draw_stale_reason reason)
{
    struct _ged_draw_obol_stale_ctx ctx = {reason};
    int obol_marked = _ged_draw_shape_ref_try_obol_paths_lazy(gedp, ref,
	    _ged_draw_obol_mark_stale_cb, &ctx,
	    "GED Obol stale source path");
    /* Obol owns public database-source stale mutation. */
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
	    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0.0, 0, NULL, 1, rgb, 0, 0);
    /* Obol owns public shape material-color mutation. */
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
	    _ged_draw_obol_set_evaluated_region_cb, &ctx,
	    "GED Obol evaluated-region mutation path");
    /* Obol owns public evaluated-region mutation. */
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
    ged_draw_view_context_lod_policy_get(&policy, first_view_ctx);

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
	    policy.csg_enabled ? 1 : 0,
	    policy.mesh_enabled ? 1 : 0,
	    bv_scale_get(bv_context_view_const((const struct bv_context *)first_view_ctx)),
	    policy.scale,
	    bv_width_get(bv_context_view_const((const struct bv_context *)first_view_ctx)),
	    bv_height_get(bv_context_view_const((const struct bv_context *)first_view_ctx)),
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
	(void)ged_draw_obol_database_source_adaptive_wireframe_realize(gedp,
		token->path, realize_dbip, fullpath, &local_tol,
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
    return NULL;
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

    int ret = (*cb)(cb_data, (const point_t *)points, point_count,
	    (const vect_t *)normals, normal_count, (const int *)indices,
	    index_count);
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
	const struct bn_tol *tol,
	void *view_ctx,
	fastf_t draw_size)
{
    if (!gedp || !path || !path[0] || !dbip || !fullpath ||
	    fullpath->fp_len <= 0 || !view_ctx)
	return 0;

    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, DB_FULL_PATH_CUR_DIR(fullpath), dbip,
	    NULL) < 0)
	return 0;

    int published = 0;
    if (intern.idb_meth && intern.idb_meth->ft_lod_realize) {
	struct rt_primitive_lod_realization realization;
	struct bv_view_info view_info = BV_VIEW_INFO_INIT;
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

    struct bv_view_info view_info = BV_VIEW_INFO_INIT;
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

    struct bv_view_info view_info = BV_VIEW_INFO_INIT;
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


int
ged_draw_shape_ref_last_point(struct ged *gedp, ged_draw_shape_ref ref, point_t out)
{
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref,
	    _ged_draw_obol_last_point_cb, out,
	    "GED Obol line last-point path"))
	return 1;

    return 0;
}


int
ged_draw_shape_ref_line_summary(struct ged *gedp,
				ged_draw_shape_ref ref,
				struct ged_draw_view_line_summary *out)
{
    struct _ged_draw_obol_line_summary_ctx ctx = {out};
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref,
	    _ged_draw_obol_line_summary_cb, &ctx,
	    "GED Obol line summary path"))
	return 1;

    return 0;
}


int
ged_draw_shape_ref_line_point_at(struct ged *gedp,
				 ged_draw_shape_ref ref,
				 size_t index,
				 point_t out)
{
    struct _ged_draw_obol_line_point_ctx ctx;
    ctx.index = index;
    ctx.out = out;
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref,
	    _ged_draw_obol_line_point_cb, &ctx,
	    "GED Obol line point path"))
	return 1;

    return 0;
}


int
ged_draw_shape_ref_line_command_at(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   size_t index,
				   int *out)
{
    struct _ged_draw_obol_line_command_ctx ctx = {index, out};
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref,
	    _ged_draw_obol_line_command_cb, &ctx,
	    "GED Obol line command path"))
	return 1;

    return 0;
}


int
ged_draw_shape_ref_geometry_summary(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    struct ged_draw_shape_geometry_summary *out)
{
    struct _ged_draw_obol_geometry_summary_ctx ctx = {out, 0,
						      GED_DRAW_MODE_WIRE};
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_handle(
	    ged_draw_registry_shape_ref_scene_handle(gedp, ref));
    if (token && token->draw_mode_valid) {
	ctx.draw_mode_valid = 1;
	ctx.draw_mode = token->draw_mode;
    }
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref,
	    _ged_draw_obol_geometry_summary_cb, &ctx,
	    "GED Obol geometry summary path"))
	return 1;

    return 0;
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

    return 0;
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

    return 0;
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

    return 0;
}


int
ged_draw_shape_context_geometry_summary(void *shape_ctx,
					struct ged_draw_shape_geometry_summary *out)
{
    struct _ged_draw_obol_geometry_summary_ctx ctx = {out, 0,
						      GED_DRAW_MODE_WIRE};
    struct ged_draw_obol_context_token *token =
	ged_draw_obol_context_from_scene_ctx(shape_ctx);
    if (token && token->draw_mode_valid) {
	ctx.draw_mode_valid = 1;
	ctx.draw_mode = token->draw_mode;
    }
    if (ged_draw_shape_context_obol_path_apply(shape_ctx,
	    _ged_draw_obol_geometry_summary_cb, &ctx,
	    "GED Obol context geometry summary path"))
	return 1;

    return 0;
}


int
ged_draw_shape_ref_translate_geometry(struct ged *gedp, ged_draw_shape_ref ref,
				      const vect_t xlate)
{
    struct _ged_draw_obol_translate_ctx ctx;
    VMOVE(ctx.xlate, xlate);
    if (_ged_draw_shape_ref_try_obol_paths(gedp, ref,
	    _ged_draw_obol_translate_cb, &ctx,
	    "GED Obol translate VLIST path"))
	return 1;

    return 0;
}


static int
ged_draw_source_view_context_scene_root_set(void *view_ctx,
					   ged_draw_scene_handle root)
{
    return ged_view_context_scene_root_ref_attach(view_ctx, root);
}


static ged_draw_scene_handle
ged_draw_obol_source_root_scene_handle(struct ged *gedp)
{
    if (!gedp || !ged_draw_obol_scene_controller_attached(gedp))
	return ged_draw_scene_handle_null();

    void *root_ctx = ged_draw_obol_context_token_for_path(gedp, "/", NULL);
    if (!root_ctx)
	return ged_draw_scene_handle_null();

    return ged_draw_scene_handle_make(root_ctx, GED_DRAW_SCENE_BACKEND_OBOL);
}


void *
ged_draw_source_root_attach_view_contexts(struct ged *gedp,
					  void *active_view_ctx,
					  struct bu_ptbl *views)
{
    if (!gedp)
	return NULL;

    ged_draw_scene_handle root_scene_handle =
	ged_draw_obol_source_root_scene_handle(gedp);
    if (ged_draw_scene_handle_is_null(root_scene_handle))
	return NULL;

    int attached = 0;

    if (active_view_ctx)
	attached += ged_draw_source_view_context_scene_root_set(
		active_view_ctx, root_scene_handle);
    if (views) {
	for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	    void *view_ctx = BU_PTBL_GET(views, i);
	    if (!view_ctx)
		continue;
	    attached += ged_draw_source_view_context_scene_root_set(view_ctx,
		    root_scene_handle);
	}
    }

    if (attached > 0)
	return ged_draw_scene_handle_context(root_scene_handle);

    return NULL;
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
