/*                B S G _ G E D _ D R A W _ R E F S . C
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
/** @file libged/b_s_g___g_e_d___d_r_a_w___r_e_f_s_._c.c
 *
 * Draw shape-ref mutation compatibility bridge.
 */

#include "common.h"

#include <stdlib.h>
#include <string.h>

#include "bu/malloc.h"
#include "bu/ptbl.h"
#include "bu/str.h"
#include "bu/color.h"
#include "bu/hash.h"
#include "bsg/appearance.h"
#include "bsg/defines.h"
#include "bsg/database_source.h"
#include "bsg/draw_ctx.h"
#include "bsg/draw_intent.h"
#include "bsg/draw_set.h"
#include "bsg/draw_source.h"
#include "bsg/field.h"
#include "bsg/geometry.h"
#include "bsg/material.h"
#include "bsg/payload.h"
#include "bg/plot3.h"
#include "bsg/scene_builder.h"
#include "bsg/scene_object.h"
#include "bsg/view_set.h"
#include "bsg/view_state.h"
#include "bg/clip.h"

#include "ged.h"
#include "ged/draw.h"
#include "rt/view.h"
#include "./ged_private.h"
#include "./bsg_ged_draw_private.h"

static bsg_scene_ref
_ged_draw_scene_ref_from_context(void *scene_ctx)
{
    bsg_scene_ref ref = BSG_SCENE_REF_NULL_INIT;
    ref.opaque = scene_ctx;
    return ref;
}


static rt_view_scene_ref
_ged_draw_scene_ref_to_rt(bsg_scene_ref ref)
{
    rt_view_scene_ref rt_ref = RT_VIEW_SCENE_REF_NULL_INIT;
    rt_ref.opaque = ref.opaque;
    return rt_ref;
}


rt_view_scene_ref
ged_draw_scene_ref_to_rt_view_ref(bsg_scene_ref ref)
{
    return _ged_draw_scene_ref_to_rt(ref);
}


bsg_scene_ref
ged_draw_scene_ref_from_rt_view_ref(rt_view_scene_ref ref)
{
    bsg_scene_ref bsg_ref = BSG_SCENE_REF_NULL_INIT;
    bsg_ref.opaque = ref.opaque;
    return bsg_ref;
}


static void *
_ged_draw_context_from_scene_ref(bsg_scene_ref ref)
{
    return bsg_scene_ref_is_null(ref) ? NULL : ref.opaque;
}


void *
ged_draw_scene_ref_context(bsg_scene_ref ref)
{
    return _ged_draw_context_from_scene_ref(ref);
}


int
ged_draw_shape_ref_set_flag(struct ged *gedp, ged_draw_shape_ref ref, int flag)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref))
	return 0;
    bsg_scene_set_visible(shape_ref, flag ? 1 : 0);
    return 1;
}


int
ged_draw_shape_ref_set_work_flag(struct ged *gedp, ged_draw_shape_ref ref, int wflag)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref))
	return 0;
    ged_draw_scene_ref_set_work_flag(shape_ref, wflag);
    return 1;
}


int
ged_draw_shape_ref_line_style(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref))
	return -1;
    return ged_draw_scene_ref_line_style(shape_ref);
}


static void
_ged_draw_shape_ref_geometry_changed(bsg_scene_ref shape_ref)
{
    ged_draw_shape_state *shape_data = ged_draw_scene_ref_shape_state(shape_ref);
    if (shape_data)
	shape_data->geometry_revision++;
    bsg_scene_invalidate(shape_ref);
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
	    if (first >= 0 && vert >= 3 && bsg_field_multi_point_at(coords, (size_t)first, out))
		have = 1;
	    first = -1;
	    vert = 0;
	    continue;
	}
	if (first < 0)
	    first = idx;
	vert++;
    }

    if (first >= 0 && vert >= 3 && bsg_field_multi_point_at(coords, (size_t)first, out))
	have = 1;

    return have;
}


void *
ged_draw_first_shape_context(struct ged *gedp)
{
    return _ged_draw_context_from_scene_ref(ged_draw_first_shape_scene_ref(gedp));
}


void *
ged_draw_shape_ref_context(struct ged *gedp, ged_draw_shape_ref ref)
{
    return _ged_draw_context_from_scene_ref(
	    ged_draw_registry_shape_scene_ref(gedp, ref));
}


void *
ged_draw_shape_ref_cache_context(struct ged *gedp, ged_draw_shape_ref ref)
{
    return _ged_draw_context_from_scene_ref(
	    ged_draw_shape_scene_ref_from_cache_ref(gedp, ref));
}


void *
ged_draw_group_ref_context(struct ged *gedp, ged_draw_group_ref ref)
{
    return _ged_draw_context_from_scene_ref(
	    ged_draw_registry_group_scene_ref(gedp, ref));
}


ged_draw_shape_ref
ged_draw_shape_ref_from_context(struct ged *gedp, void *shape_ctx)
{
    return ged_draw_shape_ref_from_scene_ref(gedp,
	    _ged_draw_scene_ref_from_context(shape_ctx));
}


int
ged_draw_shape_ref_last_point(struct ged *gedp, ged_draw_shape_ref ref, point_t out)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref) || !out)
	return 0;

    bsg_geometry_ref geometry = bsg_scene_ref_as_geometry(shape_ref);
    if (bsg_node_ref_is_null(geometry.shape.node))
	return 0;

    if (bsg_geometry_ref_kind(geometry) == BSG_GEOMETRY_NODE_INDEXED_FACE_SET)
	return _ged_draw_geometry_last_face_point(geometry, out);

    bsg_field_ref coords = bsg_geometry_ref_coordinates_field(geometry);
    size_t count = bsg_field_multi_count(coords);
    if (!count)
	return 0;
    return bsg_field_multi_point_at(coords, count - 1, out);
}


static int
_ged_draw_shape_ref_line_command_from_bsg(int command)
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
_ged_draw_scene_ref_line_geometry(bsg_scene_ref shape_ref,
				  bsg_geometry_ref *out)
{
    if (out)
	memset(out, 0, sizeof(*out));
    if (bsg_scene_ref_is_null(shape_ref) || !out)
	return 0;

    bsg_geometry_ref geometry = bsg_scene_ref_as_geometry(shape_ref);
    if (bsg_node_ref_is_null(geometry.shape.node) ||
	    bsg_geometry_ref_kind(geometry) != BSG_GEOMETRY_NODE_LINE_SET)
	return 0;

    *out = geometry;
    return 1;
}


static int
_ged_draw_shape_ref_line_geometry(struct ged *gedp,
				  ged_draw_shape_ref ref,
				  bsg_geometry_ref *out)
{
    if (!gedp || ged_draw_shape_ref_is_null(ref) || !out) {
	if (out)
	    memset(out, 0, sizeof(*out));
	return 0;
    }

    return _ged_draw_scene_ref_line_geometry(
	    ged_draw_registry_shape_scene_ref(gedp, ref), out);
}


int
ged_draw_shape_ref_line_summary(struct ged *gedp,
				ged_draw_shape_ref ref,
				struct ged_draw_view_line_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    bsg_geometry_ref geometry;
    if (!_ged_draw_shape_ref_line_geometry(gedp, ref, &geometry))
	return 0;

    out->valid = 1;
    out->point_count =
	bsg_field_multi_count(bsg_geometry_ref_coordinates_field(geometry));
    return 1;
}


int
ged_draw_shape_ref_line_point_at(struct ged *gedp,
				 ged_draw_shape_ref ref,
				 size_t index,
				 point_t out)
{
    if (!out)
	return 0;

    bsg_geometry_ref geometry;
    if (!_ged_draw_shape_ref_line_geometry(gedp, ref, &geometry))
	return 0;

    return bsg_field_multi_point_at(
	    bsg_geometry_ref_coordinates_field(geometry), index, out);
}


int
ged_draw_shape_ref_line_command_at(struct ged *gedp,
				   ged_draw_shape_ref ref,
				   size_t index,
				   int *out)
{
    if (!out)
	return 0;

    bsg_geometry_ref geometry;
    if (!_ged_draw_shape_ref_line_geometry(gedp, ref, &geometry))
	return 0;

    int command = -1;
    if (!bsg_field_multi_int_at(
	    bsg_geometry_ref_primitive_sets_field(geometry), index, &command))
	return 0;
    *out = _ged_draw_shape_ref_line_command_from_bsg(command);
    return 1;
}


static int
_ged_draw_shape_context_line_geometry(void *shape_ctx, bsg_geometry_ref *out)
{
    if (out)
	memset(out, 0, sizeof(*out));
    if (!shape_ctx || !out)
	return 0;

    bsg_scene_ref shape_ref = BSG_SCENE_REF_NULL_INIT;
    shape_ref.opaque = shape_ctx;
    return _ged_draw_scene_ref_line_geometry(shape_ref, out);
}


int
ged_draw_shape_context_line_summary(void *shape_ctx,
				    struct ged_draw_view_line_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    bsg_geometry_ref geometry;
    if (!_ged_draw_shape_context_line_geometry(shape_ctx, &geometry))
	return 0;

    out->valid = 1;
    out->point_count =
	bsg_field_multi_count(bsg_geometry_ref_coordinates_field(geometry));
    return 1;
}


int
ged_draw_shape_context_line_point_at(void *shape_ctx,
				     size_t index,
				     point_t out)
{
    if (!out)
	return 0;

    bsg_geometry_ref geometry;
    if (!_ged_draw_shape_context_line_geometry(shape_ctx, &geometry))
	return 0;

    return bsg_field_multi_point_at(
	    bsg_geometry_ref_coordinates_field(geometry), index, out);
}


int
ged_draw_shape_context_line_command_at(void *shape_ctx,
				       size_t index,
				       int *out)
{
    if (!out)
	return 0;

    bsg_geometry_ref geometry;
    if (!_ged_draw_shape_context_line_geometry(shape_ctx, &geometry))
	return 0;

    int command = -1;
    if (!bsg_field_multi_int_at(
	    bsg_geometry_ref_primitive_sets_field(geometry), index, &command))
	return 0;
    *out = _ged_draw_shape_ref_line_command_from_bsg(command);
    return 1;
}


int
ged_draw_shape_context_geometry_summary(void *shape_ctx,
					struct ged_draw_shape_geometry_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!shape_ctx)
	return 0;

    bsg_scene_ref shape_ref = BSG_SCENE_REF_NULL_INIT;
    shape_ref.opaque = shape_ctx;
    bsg_geometry_ref geometry = bsg_scene_ref_as_geometry(shape_ref);
    if (bsg_node_ref_is_null(geometry.shape.node))
	return 0;

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


int
ged_draw_shape_context_has_state(void *shape_ctx)
{
    if (!shape_ctx)
	return 0;

    return ged_draw_scene_ref_shape_state(
	    _ged_draw_scene_ref_from_context(shape_ctx)) ? 1 : 0;
}


void *
ged_draw_shape_context_source(void *shape_ctx)
{
    if (!shape_ctx)
	return NULL;

    return _ged_draw_context_from_scene_ref(
	    ged_draw_shape_source_ref(
		_ged_draw_scene_ref_from_context(shape_ctx)));
}


const struct db_full_path *
ged_draw_scene_context_fullpath(void *scene_ctx)
{
    if (!scene_ctx)
	return NULL;

    return ged_draw_scene_ref_fullpath(
	    _ged_draw_scene_ref_from_context(scene_ctx));
}


int
ged_draw_group_context_dbpath(struct ged *gedp,
			      void *group_ctx,
			      struct db_full_path *out)
{
    if (!group_ctx)
	return -1;

    return ged_draw_group_scene_ref_dbpath(gedp,
	    _ged_draw_scene_ref_from_context(group_ctx), out);
}


int
ged_draw_group_context_is_overlay(void *group_ctx)
{
    if (!group_ctx)
	return 0;

    return ged_draw_group_scene_ref_is_overlay(
	    _ged_draw_scene_ref_from_context(group_ctx));
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

    bsg_scene_ref scene_ref = BSG_SCENE_REF_NULL_INIT;
    scene_ref.opaque = scene_ctx;
    if (bsg_scene_ref_is_null(scene_ref))
	return 0;

    const struct bsg_draw_intent *intent = bsg_scene_draw_intent(scene_ref);
    bsg_database_source_ref source_ref =
	bsg_database_source_ref_from_scene(scene_ref);
    out->valid = 1;
    out->is_database_source =
	bsg_scene_is_database_source(scene_ref) ||
	bsg_database_source_ref_is_container(source_ref);
    out->has_draw_intent = intent ? 1 : 0;
    out->intent_path = intent ? bsg_draw_intent_path(intent) : NULL;
    out->intent_draw_mode = intent ? bsg_draw_intent_mode(intent) : -1;
    out->visible = bsg_scene_visible(scene_ref);
    out->highlighted = bsg_scene_highlighted(scene_ref);
    out->line_style = bsg_scene_line_style(scene_ref);
    out->line_width = bsg_scene_line_width(scene_ref);
    out->transparency = bsg_scene_transparency(scene_ref);
    out->draw_mode = bsg_scene_dmode(scene_ref);
    bsg_scene_material_get_rgb(scene_ref,
	    &out->material_color[0],
	    &out->material_color[1],
	    &out->material_color[2]);
    out->material_valid = 1;
    return 1;
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

    bsg_scene_ref scene_ref = BSG_SCENE_REF_NULL_INIT;
    scene_ref.opaque = scene_ctx;
    if (bsg_scene_ref_is_null(scene_ref))
	return 0;

    bsg_database_source_ref source_ref =
	bsg_database_source_ref_from_scene(scene_ref);
    if (bsg_database_source_ref_is_null(source_ref))
	return 0;

    struct bsg_database_source_record record;
    memset(&record, 0, sizeof(record));
    if (!bsg_database_source_record_get(source_ref, &record))
	return 0;

    out->valid = 1;
    out->is_database_source =
	bsg_scene_is_database_source(scene_ref) ||
	bsg_database_source_ref_is_container(source_ref);
    out->has_state = ged_draw_database_source_record_has_state(&record);
    out->stale = bsg_database_source_ref_is_stale(source_ref);
    out->database_path = record.database_path;
    out->source_revision = record.source_revision;
    out->inputs_revision = record.inputs_revision;
    out->realized_source_revision = record.realized_source_revision;
    out->realized_inputs_revision = record.realized_inputs_revision;
    out->realization_identity = record.realization_identity;
    if (record.stale_reason != BSG_DATABASE_SOURCE_STALE_NONE)
	out->stale_reason =
	    ged_draw_database_source_stale_reason_name(record.stale_reason);
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


int
ged_draw_scene_context_tree_summary(void *scene_ctx,
				    struct ged_draw_scene_tree_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!scene_ctx)
	return 0;

    bsg_scene_ref scene_ref = BSG_SCENE_REF_NULL_INIT;
    scene_ref.opaque = scene_ctx;
    if (bsg_scene_ref_is_null(scene_ref))
	return 0;

    enum bsg_scene_element_type scene_type = bsg_scene_ref_type(scene_ref);
    out->valid = 1;
    out->is_group = (scene_type == BSG_SCENE_ELEMENT_GROUP) ? 1 : 0;
    out->is_shape = (scene_type == BSG_SCENE_ELEMENT_SHAPE) ? 1 : 0;
    out->has_parent = bsg_scene_ref_is_null(bsg_scene_parent(scene_ref)) ? 0 : 1;
    out->draw_tree_depth = bsg_scene_draw_tree_depth(scene_ref);
    out->child_count = bsg_scene_child_count(scene_ref);
    return 1;
}


void *
ged_draw_scene_context_child_at(void *scene_ctx, size_t index)
{
    if (!scene_ctx)
	return NULL;

    bsg_scene_ref scene_ref = BSG_SCENE_REF_NULL_INIT;
    scene_ref.opaque = scene_ctx;
    if (bsg_scene_ref_is_null(scene_ref))
	return NULL;

    bsg_scene_ref child_ref = bsg_scene_child_at(scene_ref, index);
    if (bsg_scene_ref_is_null(child_ref))
	return NULL;

    return child_ref.opaque;
}


void *
ged_draw_scene_context_parent(void *scene_ctx)
{
    if (!scene_ctx)
	return NULL;

    bsg_scene_ref scene_ref = BSG_SCENE_REF_NULL_INIT;
    scene_ref.opaque = scene_ctx;
    if (bsg_scene_ref_is_null(scene_ref))
	return NULL;

    bsg_scene_ref parent_ref = bsg_scene_parent(scene_ref);
    if (bsg_scene_ref_is_null(parent_ref))
	return NULL;

    return parent_ref.opaque;
}


const char *
ged_draw_scene_context_name(void *scene_ctx)
{
    if (!scene_ctx)
	return NULL;

    bsg_scene_ref scene_ref = BSG_SCENE_REF_NULL_INIT;
    scene_ref.opaque = scene_ctx;
    if (bsg_scene_ref_is_null(scene_ref))
	return NULL;

    return bsg_scene_name(scene_ref);
}


int
ged_draw_scene_context_subtree_bounds(void *scene_ctx,
				      vect_t *min,
				      vect_t *max,
				      int include_overlays)
{
    if (!scene_ctx || !min || !max)
	return 1;

    bsg_scene_ref scene_ref = BSG_SCENE_REF_NULL_INIT;
    scene_ref.opaque = scene_ctx;
    if (bsg_scene_ref_is_null(scene_ref))
	return 1;

    return bsg_scene_subtree_bbox(scene_ref, min, max, include_overlays);
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

    points = (point_t *)bu_calloc(count, sizeof(point_t), "translated line-set points");
    commands = (int *)bu_calloc(count, sizeof(int), "translated line-set commands");
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

    points = (point_t *)bu_calloc(count, sizeof(point_t), "translated point-set points");
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


int
ged_draw_shape_ref_translate_geometry(struct ged *gedp, ged_draw_shape_ref ref,
				      const vect_t xlate)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref) || !xlate)
	return 0;

    bsg_geometry_ref geometry = bsg_scene_ref_as_geometry(shape_ref);
    if (bsg_node_ref_is_null(geometry.shape.node))
	return 0;

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
    if (ok)
	_ged_draw_shape_ref_geometry_changed(shape_ref);
    return ok;
}


int
ged_draw_shape_ref_set_center(struct ged *gedp, ged_draw_shape_ref ref,
			      const point_t center)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref) || !center)
	return 0;
    ged_draw_scene_ref_set_draw_center(shape_ref, center);
    return 1;
}


int
ged_draw_shape_ref_geometry_clear(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref))
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
    if (bsg_scene_ref_is_null(shape_ref))
	return 0;

    struct ged_draw_source_state *source = ged_draw_scene_ref_source_data(shape_ref);
    if (!source)
	return 0;

    out->dbip = source->dbip;
    out->fullpath = ged_draw_scene_ref_fullpath(shape_ref);
    out->leaf_dp = ged_draw_scene_ref_leaf_dp(shape_ref);
    out->name = ged_draw_scene_ref_name(shape_ref);
    out->tol = source->tol;
    out->ttol = source->ttol;
    return 1;
}


int
ged_draw_shape_ref_publish_line_set(struct ged *gedp,
				    ged_draw_shape_ref ref,
				    const point_t *points,
				    const int *commands,
				    size_t point_count)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref))
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
    if (bsg_scene_ref_is_null(shape_ref))
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
    struct bsg_view *v = (struct bsg_view *)view_ctx;

    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref))
	return -1;
    struct rt_view_info view_info;
    if (adaptive)
	rt_view_context_info_get(&view_info, v);
    return ged_draw_scene_ref_publish_primitive_wireframe(shape_ref, ip, ttol,
	    tol, v, adaptive ? &view_info : NULL, adaptive);
}


int
ged_draw_shape_ref_release(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref))
	return 0;

    bsg_scene_ref release_ref = ged_draw_shape_source_ref(shape_ref);
    if (bsg_scene_ref_is_null(release_ref))
	release_ref = shape_ref;
    bsg_scene_ref parent_ref = bsg_scene_parent(release_ref);
    if (!bsg_scene_ref_is_null(parent_ref))
	bsg_scene_remove_child(parent_ref, release_ref);
    bsg_scene_ref_destroy(release_ref);
    return 1;
}


int
ged_draw_shape_ref_realize_context(struct ged *gedp, ged_draw_shape_ref ref,
				   void *view_ctx)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref))
	return 0;
    ged_draw_scene_ref_realize(shape_ref, view_ctx);
    return 1;
}


int
ged_draw_shape_ref_set_view(struct ged *gedp, ged_draw_shape_ref ref,
			    struct bsg_view *v)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref))
	return 0;
    bsg_scene_set_view(shape_ref, v);
    return 1;
}


int
ged_draw_shape_ref_reset_node(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref))
	return 0;
    ged_draw_scene_ref_geometry_clear(shape_ref);
    rt_mesh_lod_free_scene_ref(_ged_draw_scene_ref_to_rt(shape_ref));
    ged_draw_scene_ref_realization_reset(shape_ref);
    bsg_scene_invalidate(shape_ref);
    return 1;
}


int
ged_draw_shape_ref_set_visible(struct ged *gedp, ged_draw_shape_ref ref,
			       int visible)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref))
	return 0;
    bsg_scene_set_visible(shape_ref, visible ? 1 : 0);
    return 1;
}


int
ged_draw_shape_ref_get_color(struct ged *gedp, ged_draw_shape_ref ref,
			     unsigned char rgb[3])
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref) || !rgb)
	return 0;
    bsg_scene_color(shape_ref, &rgb[0], &rgb[1], &rgb[2]);
    return 1;
}


int
ged_draw_shape_ref_set_color(struct ged *gedp, ged_draw_shape_ref ref,
			     const unsigned char rgb[3])
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref) || !rgb)
	return 0;
    bsg_scene_set_color(shape_ref, rgb[0], rgb[1], rgb[2]);
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
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref))
	return 0;

    out->valid = 1;
    out->material_revision = (uint64_t)bsg_scene_material_revision(shape_ref);
    bsg_scene_material_get_rgb(shape_ref,
	    &out->material_color[0],
	    &out->material_color[1],
	    &out->material_color[2]);
    return 1;
}


int
ged_draw_shape_ref_set_material_color(struct ged *gedp,
				      ged_draw_shape_ref ref,
				      const unsigned char rgb[3])
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref) || !rgb)
	return 0;
    bsg_scene_material_set_rgb(shape_ref, rgb[0], rgb[1], rgb[2]);
    return 1;
}


int
ged_draw_shape_ref_set_evaluated_region(struct ged *gedp,
					ged_draw_shape_ref ref,
					int evaluated_region)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref))
	return 0;
    bsg_scene_set_legacy_eval_flag(shape_ref, evaluated_region ? 1 : 0);
    return 1;
}


int
ged_draw_shape_ref_sync_settings(struct ged *gedp, ged_draw_shape_ref ref,
				 struct bsg_appearance_settings *settings,
				 int current_mode,
				 int *changed)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (changed)
	*changed = 0;
    if (bsg_scene_ref_is_null(shape_ref))
	return 0;
    if (!settings || settings->draw_mode != current_mode)
	return 1;

    if (ged_draw_scene_ref_line_style(shape_ref) && settings->draw_non_subtract_only) {
	if (bsg_scene_visible(shape_ref)) {
	    bsg_scene_set_visible(shape_ref, 0);
	    if (changed)
		*changed = 1;
	}
    } else {
	if (!bsg_scene_visible(shape_ref)) {
	    bsg_scene_set_visible(shape_ref, 1);
	    if (changed)
		*changed = 1;
	}
    }

    if (ged_draw_scene_ref_apply_settings(shape_ref, settings) && changed)
	*changed = 1;
    return 1;
}


int
ged_draw_shape_ref_refresh_material(struct ged *gedp, ged_draw_shape_ref ref,
				    unsigned char rgb[3])
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref) || !rgb)
	return 0;
    ged_draw_scene_ref_set_material_rgb(shape_ref, rgb);
    return 1;
}


int
ged_draw_shape_ref_adaptive_sync(struct ged *gedp, ged_draw_shape_ref ref,
				 struct bsg_view **views,
				 size_t view_count,
				 int *changed)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);

    if (changed)
	*changed = 0;
    if (bsg_scene_ref_is_null(shape_ref))
	return 0;
    if (!ged_draw_scene_ref_realization_pipeline_candidate(shape_ref))
	return 1;

    for (size_t i = 0; i < view_count; i++) {
	struct bsg_view *v = views ? views[i] : NULL;
	if (!v)
	    continue;
	ged_draw_scene_ref_mark_view_inputs_changed(shape_ref);
	bsg_scene_invalidate(shape_ref);
	if (changed)
	    *changed = 1;
	break;
    }
    return 1;
}


int
ged_draw_shape_ref_lod_ensure(struct ged *gedp, ged_draw_shape_ref ref,
			      struct bsg_view *first_view,
			      struct bsg_view **views,
			      size_t view_count)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref) || !first_view)
	return 0;

    ged_draw_view_lod_policy policy;
    rt_view_context_lod_policy_get(&policy, first_view);
    int candidate = ged_draw_scene_ref_realization_pipeline_candidate(shape_ref);
    int source_backed = ged_draw_scene_ref_source_data(shape_ref) ? 1 : 0;
    if (!candidate && !source_backed)
	return 1;
    if (!candidate) {
	struct directory *dp = ged_draw_scene_ref_leaf_dp(shape_ref);
	int mode = ged_draw_scene_ref_draw_mode(shape_ref);
	int mesh_ref = dp && policy.mesh_enabled &&
	    ((dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT &&
	      (mode == 0 || mode == 1)) ||
	     (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BREP &&
	      mode == 1));
	if (!mesh_ref && !(policy.csg_enabled && mode == 0))
	    return 1;
    }
    int realized = 0;
    for (size_t i = 0; i < view_count; i++) {
	if (views && views[i]) {
	    ged_draw_scene_ref_mark_view_inputs_changed(shape_ref);
	    bsg_scene_invalidate(shape_ref);
	    ged_draw_scene_ref_realize(shape_ref, views[i]);
	    realized = 1;
	}
    }
    if (!realized) {
	ged_draw_scene_ref_mark_view_inputs_changed(shape_ref);
	bsg_scene_invalidate(shape_ref);
	ged_draw_scene_ref_realize(shape_ref, first_view);
    }
    return 1;
}


int
ged_draw_shape_ref_pipeline_candidate(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    if (bsg_scene_ref_is_null(shape_ref))
	return 0;
    return ged_draw_scene_ref_realization_pipeline_candidate(shape_ref) ? 1 : 0;
}


void *
ged_draw_shape_ref_view_context(struct ged *gedp, ged_draw_shape_ref ref)
{
    bsg_scene_ref shape_ref = ged_draw_registry_shape_scene_ref(gedp, ref);
    return bsg_scene_ref_is_null(shape_ref) ? NULL : bsg_scene_view(shape_ref);
}


/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
