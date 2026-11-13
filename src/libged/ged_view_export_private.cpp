/*          G E D _ V I E W _ E X P O R T _ P R I V A T E . C P P
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
/** @file ged_view_export_private.cpp
 *
 * Renderer-neutral view export and semantic path picking adapters.  Heavy
 * retained geometry remains owned by the rendering backend; this unit copies
 * only data explicitly requested by an export consumer.
 */

#include "common.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "BObol/BDatabaseSource.h"
#include "BObol/BSceneController.h"
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
#include "./ged_scene_record_api_private.h"
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
#include "rt/primitives/arb8.h"
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
#include "./ged_bobol_private.hpp"
#include "./ged_private.h"
#include "./ged_draw_private.h"
#include "./ged_draw_shape_state_private.h"
#include "./ged_draw_view_private.h"


static int
ged_draw_obol_scene_controller_is_owned(struct ged *gedp)
{
    return (ged_draw_obol_scene_controller_owned_internal(gedp) &&
	    ged_draw_obol_scene_controller_full_synced(gedp)) ? 1 : 0;
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

    if (database_source) {
	point_t *compact_points = NULL;
	int *compact_commands = NULL;
	size_t compact_count = 0;
	if (ged_draw_obol_database_source_line_data_copy_for_path(gedp,
		shape_path, &compact_points, &compact_commands,
		&compact_count)) {
	    detail->geometry_kind = GED_DRAW_VIEW_EXPORT_GEOMETRY_LINE_SET;
	    detail->arrays.points = compact_points;
	    detail->arrays.commands = compact_commands;
	    detail->arrays.point_count = compact_count;
	    detail->arrays.command_count = compact_count;
	    size_t structure_count = 0;
	    for (size_t i = 0; i < compact_count; i++) {
		if (compact_commands[i] != GED_DRAW_VIEW_LINE_MOVE)
		    structure_count++;
	    }
	    detail->record.vlist_structure_count = structure_count;
	    detail->record.vlist_point_count = compact_count;
	    _ged_draw_view_export_record_bounds_from_points(&detail->record,
		compact_points, compact_count);
	    return 1;
	}
    }

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
	struct ged_view_context *view_ctx,
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
	rec->selected = ged_selection_contains_path(
		view_ctx, GED_SELECTION_SELECTED_PATH,
		record_path) ? 1 : 0;
	if (!rec->selected && shape_path &&
		(!record_path || !BU_STR_EQUAL(record_path, shape_path))) {
	    rec->selected = ged_selection_contains_path(
		    view_ctx, GED_SELECTION_SELECTED_PATH,
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
    struct ged_view_context *view_ctx;
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
	struct ged_view_context *view_ctx)
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
	snprintf(fallback, sizeof(fallback), "%p", (void *)view_ctx);
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
	rec->selected = ged_selection_contains_path(
		ctx->view_ctx, GED_SELECTION_SELECTED_PATH,
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
	struct ged_view_context *view_ctx,
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
	struct ged_view_context *view_ctx,
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
    rec->selected = ged_selection_contains_path(view_ctx,
	    GED_SELECTION_SELECTED_PATH, feature->name) ? 1 : 0;
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
    struct ged_view_context *view_ctx;
    ged_draw_view_db_object_record_cb cb;
    void *userdata;
    int keep_going;
};


struct ged_draw_obol_polygon_export_ctx {
    struct ged_view_context *view_ctx;
    unsigned int query_flags;
    const char *glob;
    int draw_mode;
    ged_draw_view_db_object_record_cb cb;
    void *userdata;
    int keep_going;
};


static int
_ged_draw_view_export_polygon_record_cb(
	ged_view_polygon_ref UNUSED(ref),
	const struct ged_view_polygon_record *polygon,
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
    rec->selected = ged_selection_contains_path(
	    ctx->view_ctx, GED_SELECTION_SELECTED_PATH,
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
	struct ged_view_context *view_ctx,
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
    ged_view_polygon_visit_records(view_ctx,
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
	struct ged_view_context *view_ctx,
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
	struct ged_view_context *view_ctx,
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
	struct ged_view_context *view_ctx,
	unsigned int query_flags,
	unsigned int render_flags,
	const char *glob,
	int draw_mode,
	ged_draw_view_db_object_record_cb cb,
	void *userdata)
{
    if (!view_ctx || !cb)
	return;

    struct ged *gedp = ged_view_context_owner(view_ctx);
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

struct ged_draw_rt_pick_ctx {
    struct ged_view_context *view_ctx;
    const char *pattern;
    struct ged_pick_result *result;
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
    if (ged_pick_result_append_path(ctx->result, source_path, 0.0))
	ctx->count++;
    return 1;
}

struct ged_pick_result *
ged_pick_semantic_path(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	const char *path_pattern)
{
    if (!view_ctx || !path_pattern || !path_pattern[0] ||
	    !ged_draw_obol_scene_controller_is_owned(gedp))
	return NULL;

    struct ged_pick_result *result = ged_pick_result_create();
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

    ged_pick_result_free(result);
    return NULL;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
