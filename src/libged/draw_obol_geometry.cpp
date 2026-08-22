/*           D R A W _ O B O L _ G E O M E T R Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file draw_obol_geometry.cpp
 *
 * Retained-geometry publication, inspection, and placement at the private
 * GED/Obol boundary.  Transaction reduction and progressive realization do
 * not belong in this unit.
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BADC.h"
#include "BObol/BDrawCache.h"
#include "BObol/BGrid.h"
#include "BObol/BInit.h"
#include "BObol/BImageSource.h"
#include "BObol/BDisplayEndpoint.h"
#include "BObol/BLodRealization.h"
#include "BObol/BLodService.h"
#include "BObol/BMeshShape.h"
#include "BObol/BSceneController.h"
#include "BObol/BSceneGroup.h"
#include "BObol/BSnapAction.h"
#include "BObol/BSourceRealization.h"
#include "BObol/BViewportImage.h"
#include "BObol/BVListShape.h"
#include "BObol/BViewAttachment.h"
#include "BObol/BViewController.h"
#include "BObol/BViewQuery.h"
#include "BObol/BViewStore.h"
#include "bg/line_layer.h"
#include "bg/plane.h"
#include "bg/polygon.h"
#include "bv.h"
#include "bu/hash.h"
#include "bu/malloc.h"
#include "bu/parallel.h"
#include "bu/path.h"
#include "bu/str.h"
#include "bu/datetime.h"
#include "bu/units.h"
#include "bu/vls.h"
#include "ged.h"
#include "ged/db_index.h"
#include "ged/draw.h"
#include "ged/display.h"
#include "ged/selection.h"
#include "ged/view.h"
#include "icv.h"
#include "rt/db5.h"
#include "rt/db_fullpath.h"
#include "rt/search.h"
#include "rt/view.h"
#include "vmath.h"

#include "./ged_bobol_private.hpp"
#include "./draw_obol_bridge_private.hpp"
#include "./draw_obol_overlay_private.hpp"
#include "./ged_draw_view_private.h"
#include "./ged_private.h"

#include <algorithm>
#include <Inventor/SbVec2f.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Obol/cad/SoCADAssembly.h>
#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <atomic>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

static const char *ged_obol_leaf_name_from_path(const char *path);

static SoBRLVListShape *
ged_obol_owned_vlist_shape_for_source(SoBRLDatabaseSource *source,
				      const char *fallback_path,
				      int create)
{
    if (!source)
	return NULL;

    SoBRLVListShape *fallback = NULL;
    const int count = source->getRealizedShapeCount();
    for (int i = 0; i < count; i++) {
	SoBRLVListShape *shape = source->getRealizedShape(i);
	if (!shape)
	    continue;
	if (!fallback)
	    fallback = shape;
	const SoBRLVListShape *geom = shape->getGeometrySource();
	if (geom->point.getNum() > 0 || geom->command.getNum() > 0)
	    return shape;
    }

    if (fallback || !create)
	return fallback;

    const char *source_path = source->path.getValue().getString();
    if (!source_path || !source_path[0])
	source_path = fallback_path ? fallback_path : "";
    const char *source_name = ged_obol_leaf_name_from_path(source_path);

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->sourcePath = source_path;
    shape->sourceName = source_name;
    shape->sourceType = "line-set";
    shape->sourceId = source->realizedRevision.getValue();
    shape->displayName = source_name;
    shape->geometryName = source_name;
    shape->sourceIdentity = source_path;
    shape->cacheIdentity = source_path;
    shape->databaseIntent = TRUE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = FALSE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = FALSE;
    shape->drawMode = ged_obol_lod_draw_mode_from_database_source(source);
    shape->recordRole = "database";
    shape->geometryKind = "line-set";
    shape->visible = source->visible.getValue();
    shape->highlighted = source->highlighted.getValue();
    shape->lineStyle = source->lineStyle.getValue();
    shape->lineWidth = source->lineWidth.getValue();
    shape->transparency = source->transparency.getValue();
    shape->hiddenLine = shape->drawMode.getValue() ==
			BOBOL_LOD_DRAW_HIDDEN_LINE ? TRUE : FALSE;
    shape->colorOverride = source->colorOverride.getValue();
    shape->color = source->color.getValue();
    shape->materialColorValid = source->materialColorValid.getValue();
    shape->materialColor = source->materialColor.getValue();
    shape->materialRevision = source->materialRevision.getValue();
    shape->drawMatrixValid = source->drawMatrixValid.getValue();
    shape->drawMatrix = source->drawMatrix.getValue();
    shape->drawCenterValid = source->drawCenterValid.getValue();
    shape->drawCenter = source->drawCenter.getValue();
    shape->drawSizeValid = source->drawSizeValid.getValue();
    shape->drawSize = source->drawSize.getValue();
    source->addChild(shape);
    return shape;
}

static SoBRLVListShape *
ged_obol_owned_vlist_shape_for_path(struct ged *gedp, const char *path)
{
    return ged_obol_owned_vlist_shape_for_source(
	       ged_obol_owned_database_source_for_path(gedp, path), path, 0);
}


int
ged_obol_vlist_shape_is_annotation(SoBRLVListShape *shape)
{
    if (!shape)
	return 0;

    const char *kind = shape->geometryKind.getValue().getString();
    if (kind && BU_STR_EQUAL(kind, "annotation"))
	return 1;

    const char *source_type = shape->sourceType.getValue().getString();
    return (source_type && BU_STR_EQUAL(source_type, "annotation")) ? 1 : 0;
}


int
ged_obol_vlist_shape_has_annotation_record(SoBRLVListShape *shape)
{
    if (!shape)
	return 0;
    const SoBRLVListShape *geom = shape->getGeometrySource();
    return (geom->annotationPoint.getNum() > 0 ||
	    geom->annotationSegmentKind.getNum() > 0) ? 1 : 0;
}


SoBRLVListShape *
ged_obol_owned_annotation_vlist_shape_for_source(
    SoBRLDatabaseSource *source,
    const char *fallback_path)
{
    if (!source)
	return NULL;

    SoBRLVListShape *annotation_shape = NULL;
    const int count = source->getRealizedShapeCount();
    for (int i = 0; i < count; i++) {
	SoBRLVListShape *shape = source->getRealizedShape(i);
	if (!ged_obol_vlist_shape_is_annotation(shape))
	    continue;
	if (ged_obol_vlist_shape_has_annotation_record(shape))
	    return shape;
	if (!annotation_shape)
	    annotation_shape = shape;
    }

    return annotation_shape ? annotation_shape :
	   ged_obol_owned_vlist_shape_for_source(source, fallback_path, 0);
}


static SoBRLVListShape *
ged_obol_owned_annotation_vlist_shape_for_path(struct ged *gedp,
	const char *path)
{
    return ged_obol_owned_annotation_vlist_shape_for_source(
	       ged_obol_owned_database_source_for_path(gedp, path), path);
}

static const char *
ged_obol_leaf_name_from_path(const char *path)
{
    if (!path || !path[0])
	return "";
    const char *leaf = strrchr(path, '/');
    return (leaf && leaf[1]) ? leaf + 1 : path;
}

static SoBRLMeshShape *
ged_obol_owned_mesh_shape_for_source(SoBRLDatabaseSource *source,
				     const char *fallback_path,
				     int create)
{
    if (!source)
	return NULL;

    SoBRLMeshShape *shape = source->getRealizedMesh();
    if (shape || !create)
	return shape;

    const char *source_path = source->path.getValue().getString();
    if (!source_path || !source_path[0])
	source_path = fallback_path ? fallback_path : "";
    const char *source_name = ged_obol_leaf_name_from_path(source_path);

    shape = new SoBRLMeshShape;
    shape->sourcePath = source_path;
    shape->sourceName = source_name;
    shape->sourceType = "indexed-face-set";
    shape->sourceId = source->realizedRevision.getValue();
    shape->displayName = source_name;
    shape->geometryName = source_name;
    shape->sourceIdentity = source_path;
    shape->cacheIdentity = source_path;
    shape->databaseIntent = TRUE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = FALSE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = FALSE;
    shape->drawMode = ged_obol_lod_draw_mode_from_database_source(source);
    shape->recordRole = "database";
    shape->geometryKind = "surface";
    shape->visible = source->visible.getValue();
    shape->highlighted = source->highlighted.getValue();
    shape->lineStyle = source->lineStyle.getValue();
    shape->lineWidth = source->lineWidth.getValue();
    shape->transparency = source->transparency.getValue();
    shape->hiddenLine = shape->drawMode.getValue() ==
			BOBOL_LOD_DRAW_HIDDEN_LINE ? TRUE : FALSE;
    shape->colorOverride = source->colorOverride.getValue();
    shape->color = source->color.getValue();
    shape->materialColorValid = source->materialColorValid.getValue();
    shape->materialColor = source->materialColor.getValue();
    shape->materialRevision = source->materialRevision.getValue();
    shape->drawMatrixValid = source->drawMatrixValid.getValue();
    shape->drawMatrix = source->drawMatrix.getValue();
    shape->drawCenterValid = source->drawCenterValid.getValue();
    shape->drawCenter = source->drawCenter.getValue();
    shape->drawSizeValid = source->drawSizeValid.getValue();
    shape->drawSize = source->drawSize.getValue();
    source->addChild(shape);
    return shape;
}

static SoBRLMeshShape *
ged_obol_owned_mesh_shape_for_path(struct ged *gedp, const char *path,
				   int create)
{
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    return ged_obol_owned_mesh_shape_for_source(source, path, create);
}

static uint64_t
ged_obol_hash_sb_string(const SbString &value)
{
    const char *str = value.getString();
    if (!str || !str[0])
	return 0;
    return bu_data_hash(str, strlen(str) * sizeof(char));
}

static uint64_t
ged_obol_hash_cstr(const char *str)
{
    if (!str || !str[0])
	return 0;
    return bu_data_hash(str, strlen(str) * sizeof(char));
}

static int
ged_obol_vlist_command_to_ged(int command)
{
    if (command == SoBRLVListShape::MOVE)
	return GED_DRAW_VIEW_LINE_MOVE;
    if (command == SoBRLVListShape::DRAW)
	return GED_DRAW_VIEW_LINE_DRAW;
    if (command == SoBRLVListShape::POINT)
	return GED_DRAW_VIEW_LINE_POINT_DRAW;
    return command;
}

int32_t
ged_obol_vlist_command_from_ged(int command, size_t index)
{
    if (command == GED_DRAW_VIEW_LINE_MOVE)
	return SoBRLVListShape::MOVE;
    if (command == GED_DRAW_VIEW_LINE_DRAW)
	return SoBRLVListShape::DRAW;
    if (command == GED_DRAW_VIEW_LINE_POINT_DRAW)
	return SoBRLVListShape::POINT;
    if (command < 0 && (index % 2) == 0)
	return SoBRLVListShape::MOVE;
    if (command < 0)
	return SoBRLVListShape::DRAW;
    return -1;
}

static void
ged_obol_vlist_shape_set_precise_points(SoBRLVListShape *shape,
					const point_t *points,
					size_t point_count)
{
    if (!shape || !points || point_count > static_cast<size_t>(INT_MAX))
	return;

    std::vector<double> precise_points;
    precise_points.reserve(point_count * 3);
    for (size_t i = 0; i < point_count; i++) {
	precise_points.push_back(points[i][0]);
	precise_points.push_back(points[i][1]);
	precise_points.push_back(points[i][2]);
    }

    shape->setPrecisePoints(precise_points.empty() ? NULL :
			    precise_points.data(), static_cast<int>(point_count));
}

template <typename ShapeT>
static void
ged_obol_local_shape_apply_common_state(
    ShapeT *shape,
    const char *shape_path,
    const char *display_name,
    const char *source_type,
    const char *geometry_kind,
    const struct ged_draw_scene_display_summary *display_state)
{
    if (!shape || !shape_path)
	return;

    const char *leaf = ged_obol_leaf_name_from_path(shape_path);
    const char *name = (display_name && display_name[0]) ? display_name : leaf;
    const int local_draw_mode = (display_state && display_state->valid) ?
				display_state->draw_mode : GED_DRAW_MODE_WIRE;

    shape->sourcePath = shape_path;
    shape->sourceName = leaf;
    shape->sourceType = source_type ? source_type : "local";
    shape->sourceId = static_cast<uint32_t>(ged_obol_hash_cstr(shape_path));
    shape->displayName = name ? name : leaf;
    shape->geometryName = name ? name : leaf;
    shape->sourceIdentity = shape_path;
    shape->cacheIdentity = shape_path;
    const char *intent_path = (display_state && display_state->valid &&
			       display_state->intent_path && display_state->intent_path[0]) ?
			      display_state->intent_path : "";
    shape->ownerSourcePath = intent_path;
    shape->databaseIntent = intent_path[0] ? TRUE : FALSE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = TRUE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = TRUE;
    shape->drawMode = ged_obol_lod_draw_mode_from_ged(local_draw_mode);
    shape->recordRole = "local";
    shape->geometryKind = geometry_kind ? geometry_kind : "local";
    shape->visible = (!display_state || !display_state->valid ||
		      display_state->visible) ? TRUE : FALSE;
    shape->highlighted = (display_state && display_state->valid &&
			  display_state->highlighted) ? TRUE : FALSE;
    shape->lineStyle = (display_state && display_state->valid) ?
		       display_state->line_style : 0;
    shape->lineWidth = (display_state && display_state->valid) ?
		       display_state->line_width : 0;
    shape->transparency = (display_state && display_state->valid) ?
			  static_cast<float>(display_state->transparency) : 0.0f;
    shape->colorOverride = FALSE;
    shape->color = SbColor(1.0f, 1.0f, 1.0f);
    shape->materialColorValid = (display_state && display_state->valid &&
				 display_state->material_valid) ? TRUE : FALSE;
    if (display_state && display_state->valid &&
	display_state->material_valid) {
	const SbColor material =
	    ged_obol_color_from_rgb(display_state->material_color);
	shape->materialColor = material;
	shape->color = material;
	shape->colorOverride = TRUE;
    }
    shape->materialRevision = 0;
}

static SoBRLVListShape *
ged_obol_local_vlist_shape_for_path(
    BObolSceneController *scene,
    const char *group_path,
    const char *shape_path)
{
    if (!scene || !group_path || !group_path[0] ||
	!shape_path || !shape_path[0])
	return NULL;

    if (!scene->ensureGroup(group_path))
	return NULL;

    SoNode *node = scene->findShape(shape_path);
    if (node && !node->isOfType(SoBRLVListShape::getClassTypeId())) {
	(void)scene->removeShape(shape_path);
	node = NULL;
    }
    if (node) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	(void)scene->moveShapeToGroup(shape_path, group_path);
	return shape;
    }

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->sourcePath = shape_path;
    const int appended = scene->appendChildToGroup(group_path, shape);
    if (appended < 0)
	return NULL;
    return shape;
}

static SoBRLMeshShape *
ged_obol_local_mesh_shape_for_path(
    BObolSceneController *scene,
    const char *group_path,
    const char *shape_path)
{
    if (!scene || !group_path || !group_path[0] ||
	!shape_path || !shape_path[0])
	return NULL;

    if (!scene->ensureGroup(group_path))
	return NULL;

    SoNode *node = scene->findShape(shape_path);
    if (node && !node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	(void)scene->removeShape(shape_path);
	node = NULL;
    }
    if (node) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	(void)scene->moveShapeToGroup(shape_path, group_path);
	return shape;
    }

    SoBRLMeshShape *shape = new SoBRLMeshShape;
    shape->sourcePath = shape_path;
    const int appended = scene->appendChildToGroup(group_path, shape);
    if (appended < 0)
	return NULL;
    return shape;
}

extern "C" int
ged_draw_obol_local_shape_publish_line_set_for_path(
    struct ged *gedp,
    const char *group_path,
    const char *shape_path,
    const char *display_name,
    const point_t *points,
    const int *commands,
    size_t point_count,
    const struct ged_draw_scene_display_summary *display_state)
{
    if (!group_path || !group_path[0] || !shape_path || !shape_path[0] ||
	point_count > static_cast<size_t>(INT_MAX))
	return 0;
    if (point_count && !points)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    SoBRLVListShape *shape =
	ged_obol_local_vlist_shape_for_path(scene, group_path, shape_path);
    if (!shape)
	return 0;

    std::vector<SbVec3f> obol_points;
    std::vector<int32_t> obol_commands;
    obol_points.reserve(point_count);
    obol_commands.reserve(point_count);
    for (size_t i = 0; i < point_count; i++) {
	const int command = commands ? commands[i] : -1;
	const int32_t obol_command =
	    ged_obol_vlist_command_from_ged(command, i);
	if (obol_command < 0)
	    return 0;

	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
	obol_commands.push_back(obol_command);
    }

    ged_obol_local_shape_apply_common_state(shape, shape_path, display_name,
					    "local-line-set", "line", display_state);
    shape->setLineSet(obol_points.empty() ? NULL : obol_points.data(),
		      obol_commands.empty() ? NULL : obol_commands.data(),
		      static_cast<int>(point_count));
    ged_obol_vlist_shape_set_precise_points(shape, points, point_count);
    return 1;
}

extern "C" int
ged_draw_obol_local_shape_set_record_role_for_path(
    struct ged *gedp,
    const char *shape_path,
    const char *record_role)
{
    if (!shape_path || !shape_path[0] || !record_role)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    SoNode *node = scene->findShape(shape_path);
    if (!node)
	return 0;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	shape->recordRole = record_role;
	return 1;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	shape->recordRole = record_role;
	return 1;
    }

    return 0;
}

extern "C" int
ged_draw_obol_shape_display_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_scene_display_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node)
	return 0;
    const char *record_role = ged_obol_shape_node_record_role(node);
    if (record_role && (BU_STR_EQUAL(record_role, "view-feature") ||
			BU_STR_EQUAL(record_role, "view-polygon")))
	return 0;

    out->valid = 1;
    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	const char *intent_path = shape->ownerSourcePath.getValue().getString();
	out->is_database_source = shape->databaseIntent.getValue() ? 1 : 0;
	out->has_draw_intent = 1;
	out->intent_path = (intent_path && intent_path[0]) ? intent_path : path;
	out->intent_draw_mode =
	    ged_obol_lod_draw_mode_to_ged(shape->drawMode.getValue());
	out->visible = shape->visible.getValue() ? 1 : 0;
	out->selected = shape->selected.getValue() ? 1 : 0;
	out->highlighted = shape->highlighted.getValue() ? 1 : 0;
	out->line_style = shape->lineStyle.getValue();
	out->line_width = shape->lineWidth.getValue();
	out->transparency = shape->transparency.getValue();
	out->draw_mode =
	    ged_obol_lod_draw_mode_to_ged(shape->drawMode.getValue());
	const SbColor material = shape->materialColorValid.getValue() ?
				 shape->materialColor.getValue() : shape->color.getValue();
	ged_obol_rgb_from_color(material, out->material_color);
	out->material_valid = 1;
	return 1;
    }
    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	const char *intent_path = shape->ownerSourcePath.getValue().getString();
	out->is_database_source = shape->databaseIntent.getValue() ? 1 : 0;
	out->has_draw_intent = 1;
	out->intent_path = (intent_path && intent_path[0]) ? intent_path : path;
	out->intent_draw_mode =
	    ged_obol_lod_draw_mode_to_ged(shape->drawMode.getValue());
	out->visible = shape->visible.getValue() ? 1 : 0;
	out->selected = shape->selected.getValue() ? 1 : 0;
	out->highlighted = shape->highlighted.getValue() ? 1 : 0;
	out->line_style = shape->lineStyle.getValue();
	out->line_width = shape->lineWidth.getValue();
	out->transparency = shape->transparency.getValue();
	out->draw_mode =
	    ged_obol_lod_draw_mode_to_ged(shape->drawMode.getValue());
	const SbColor material = shape->materialColorValid.getValue() ?
				 shape->materialColor.getValue() : shape->color.getValue();
	ged_obol_rgb_from_color(material, out->material_color);
	out->material_valid = 1;
	return 1;
    }

    memset(out, 0, sizeof(*out));
    return 0;
}

extern "C" int
ged_draw_obol_shape_geometry_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_shape_geometry_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node)
	return 0;

    out->valid = 1;
    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	const SoBRLVListShape *geom = shape->getGeometrySource();
	const char *kind = shape->geometryKind.getValue().getString();
	out->geometry_name = (kind && BU_STR_EQUAL(kind, "annotation")) ?
			     "annotation" : "line-set";
	out->point_count = static_cast<size_t>(geom->point.getNum());
	out->index_count = 0;
	return 1;
    }
    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	const SoBRLMeshShape *geom = shape->getGeometrySource();
	out->geometry_name = "indexed-face-set";
	out->point_count = static_cast<size_t>(geom->point.getNum());
	out->index_count =
	    static_cast<size_t>(shape->getTriangleCount()) * 4;
	return 1;
    }

    memset(out, 0, sizeof(*out));
    return 0;
}

extern "C" int
ged_draw_obol_shape_surface_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_view_surface_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node || !node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return 0;

    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
    const size_t triangle_count =
	static_cast<size_t>(shape->getTriangleCount());
    const SoBRLMeshShape *geom = shape->getGeometrySource();

    out->valid = 1;
    out->point_count = static_cast<size_t>(geom->point.getNum());
    out->normal_count = triangle_count * 3;
    out->index_count = triangle_count * 4;
    out->face_count = triangle_count;
    out->normals_per_index = 1;
    out->material_valid = 1;
    out->material_draw_mode =
	ged_obol_lod_draw_mode_to_ged(shape->drawMode.getValue());
    out->material_transparency = shape->transparency.getValue();
    out->material_highlighted = shape->highlighted.getValue() ? 1 : 0;
    const SbColor material = shape->materialColorValid.getValue() ?
			     shape->materialColor.getValue() : shape->color.getValue();
    ged_obol_rgb_from_color(material, out->material_color);
    out->cache_identity = ged_obol_hash_cstr(path);
    const char *source_path = shape->sourcePath.getValue().getString();
    if (!source_path || !source_path[0])
	source_path = shape->ownerSourcePath.getValue().getString();
    out->source_identity = ged_obol_hash_cstr(source_path);
    if (!out->source_identity)
	out->source_identity = out->cache_identity;
    return 1;
}

extern "C" int
ged_draw_obol_shape_surface_point_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    point_t out)
{
    if (!out)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node || !node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return 0;

    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
    const SoBRLMeshShape *geom = shape->getGeometrySource();
    if (index >= static_cast<size_t>(geom->point.getNum()))
	return 0;

    const SbVec3f &point = geom->point[static_cast<int>(index)];
    out[0] = point[0];
    out[1] = point[1];
    out[2] = point[2];
    return 1;
}

extern "C" int
ged_draw_obol_shape_surface_index_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    int *out)
{
    if (!out)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node || !node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return 0;

    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
    const SoBRLMeshShape *geom = shape->getGeometrySource();
    const size_t triangle_count =
	static_cast<size_t>(shape->getTriangleCount());
    if (index >= triangle_count * 4)
	return 0;

    const size_t face_offset = index % 4;
    if (face_offset == 3) {
	*out = -1;
	return 1;
    }

    const size_t coord_index = (index / 4) * 3 + face_offset;
    if (coord_index >= static_cast<size_t>(geom->coordIndex.getNum()))
	return 0;
    *out = geom->coordIndex[static_cast<int>(coord_index)];
    return 1;
}

extern "C" int
ged_draw_obol_shape_line_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_view_line_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node || !node->isOfType(SoBRLVListShape::getClassTypeId()))
	return 0;

    SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
    const SoBRLVListShape *geom = shape->getGeometrySource();
    out->valid = 1;
    out->point_count = static_cast<size_t>(geom->point.getNum());
    return 1;
}

extern "C" int
ged_draw_obol_shape_line_point_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    point_t out)
{
    if (!out)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node || !node->isOfType(SoBRLVListShape::getClassTypeId()))
	return 0;

    SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
    const SoBRLVListShape *geom = shape->getGeometrySource();
    if (index >= static_cast<size_t>(geom->point.getNum()))
	return 0;

    double precise[3];
    if (shape->getPrecisePoint(static_cast<int>(index), precise)) {
	out[0] = precise[0];
	out[1] = precise[1];
	out[2] = precise[2];
	return 1;
    }

    const SbVec3f &point = geom->point[static_cast<int>(index)];
    out[0] = point[0];
    out[1] = point[1];
    out[2] = point[2];
    return 1;
}

extern "C" int
ged_draw_obol_shape_line_command_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    int *out)
{
    if (!out)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node || !node->isOfType(SoBRLVListShape::getClassTypeId()))
	return 0;

    SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
    const SoBRLVListShape *geom = shape->getGeometrySource();
    if (index >= static_cast<size_t>(geom->command.getNum()))
	return 0;

    *out = ged_obol_vlist_command_to_ged(
	       geom->command[static_cast<int>(index)]);
    return 1;
}

extern "C" int
ged_draw_obol_database_source_last_point_for_path(
    struct ged *gedp,
    const char *path,
    point_t out)
{
    if (!out)
	return 0;

    SoBRLVListShape *shape =
	ged_obol_owned_annotation_vlist_shape_for_path(gedp, path);
    const SoBRLVListShape *geom = shape ? shape->getGeometrySource() : NULL;
    if (!geom || geom->point.getNum() <= 0)
	return 0;

    const SbVec3f &point = geom->point[geom->point.getNum() - 1];
    out[0] = point[0];
    out[1] = point[1];
    out[2] = point[2];
    return 1;
}

extern "C" int
ged_draw_obol_database_source_line_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_view_line_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    SoBRLVListShape *shape =
	ged_obol_owned_annotation_vlist_shape_for_source(source, path);

    out->valid = 1;
    if (!shape)
	return 1;

    const SoBRLVListShape *geom = shape->getGeometrySource();
    out->point_count = static_cast<size_t>(geom->point.getNum());
    out->cache_identity =
	ged_obol_hash_sb_string(shape->cacheIdentity.getValue());
    out->source_identity =
	ged_obol_hash_sb_string(shape->sourceIdentity.getValue());
    return 1;
}

extern "C" int
ged_draw_obol_database_source_line_point_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    point_t out)
{
    if (!out)
	return 0;

    SoBRLVListShape *shape =
	ged_obol_owned_annotation_vlist_shape_for_path(gedp, path);
    const SoBRLVListShape *geom = shape ? shape->getGeometrySource() : NULL;
    if (!geom || index >= static_cast<size_t>(geom->point.getNum()))
	return 0;

    double precise[3];
    if (shape->getPrecisePoint(static_cast<int>(index), precise)) {
	out[0] = precise[0];
	out[1] = precise[1];
	out[2] = precise[2];
	return 1;
    }

    const SbVec3f &point = geom->point[static_cast<int>(index)];
    out[0] = point[0];
    out[1] = point[1];
    out[2] = point[2];
    return 1;
}

extern "C" int
ged_draw_obol_database_source_line_command_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    int *out)
{
    if (!out)
	return 0;

    SoBRLVListShape *shape = ged_obol_owned_vlist_shape_for_path(gedp, path);
    const SoBRLVListShape *geom = shape ? shape->getGeometrySource() : NULL;
    if (!geom || index >= static_cast<size_t>(geom->command.getNum()))
	return 0;

    *out = ged_obol_vlist_command_to_ged(
	       geom->command[static_cast<int>(index)]);
    return 1;
}

extern "C" int
ged_draw_obol_database_source_line_data_copy_for_path(
    struct ged *gedp,
    const char *path,
    point_t **points,
    int **commands,
    size_t *point_count)
{
    if (points)
	*points = NULL;
    if (commands)
	*commands = NULL;
    if (point_count)
	*point_count = 0;
    if (!points || !commands || !point_count)
	return 0;

    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source || !source->hasCompactInstanceIndex())
	return 0;

    std::vector<SbVec3f> compactPoints;
    std::vector<int32_t> compactCommands;
    if (!source->copyCompactWireGeometry(compactPoints, compactCommands) ||
	compactPoints.size() != compactCommands.size())
	return 0;

    point_t *copiedPoints = (point_t *)bu_calloc(compactPoints.size(),
	sizeof(point_t), "GED Obol compact export points");
    int *copiedCommands = (int *)bu_calloc(compactCommands.size(),
	sizeof(int), "GED Obol compact export commands");
    for (size_t i = 0; i < compactPoints.size(); i++) {
	copiedPoints[i][X] = compactPoints[i][0];
	copiedPoints[i][Y] = compactPoints[i][1];
	copiedPoints[i][Z] = compactPoints[i][2];
	copiedCommands[i] = compactCommands[i] == 0 ?
	    GED_DRAW_VIEW_LINE_MOVE : GED_DRAW_VIEW_LINE_DRAW;
    }

    *points = copiedPoints;
    *commands = copiedCommands;
    *point_count = compactPoints.size();
    return 1;
}

extern "C" int
ged_draw_obol_database_source_surface_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_view_surface_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    SoBRLMeshShape *shape = ged_obol_owned_mesh_shape_for_source(source,
			    path, 0);
    if (!shape) {
	out->valid = 1;
	return 1;
    }

    int material_draw_mode =
	ged_obol_lod_draw_mode_to_ged(shape->drawMode.getValue());
    BObolDatabaseSourceSummary source_summary;
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (ged_obol_database_source_controller_summary_for_source(scene,
	    source, source_summary) && source_summary.valid)
	material_draw_mode =
	    ged_obol_database_source_exact_draw_mode_to_ged(gedp,
		source_summary, source);

    const size_t triangle_count =
	static_cast<size_t>(shape->getTriangleCount());
    const SoBRLMeshShape *geom = shape->getGeometrySource();
    out->valid = 1;
    out->point_count = static_cast<size_t>(geom->point.getNum());
    out->normal_count = triangle_count * 3;
    out->index_count = triangle_count * 4;
    out->face_count = triangle_count;
    out->normals_per_index = 1;
    out->material_valid = 1;
    out->material_draw_mode = material_draw_mode;
    out->material_transparency = shape->transparency.getValue();
    out->material_highlighted = shape->highlighted.getValue() ? 1 : 0;
    const SbColor material = shape->materialColorValid.getValue() ?
			     shape->materialColor.getValue() : shape->color.getValue();
    ged_obol_rgb_from_color(material, out->material_color);
    out->cache_identity =
	ged_obol_hash_sb_string(shape->cacheIdentity.getValue());
    out->source_identity =
	ged_obol_hash_sb_string(shape->sourceIdentity.getValue());
    return 1;
}

extern "C" int
ged_draw_obol_database_source_surface_point_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    point_t out)
{
    if (!out)
	return 0;

    SoBRLMeshShape *shape = ged_obol_owned_mesh_shape_for_path(gedp, path, 0);
    const SoBRLMeshShape *geom = shape ? shape->getGeometrySource() : NULL;
    if (!geom || index >= static_cast<size_t>(geom->point.getNum()))
	return 0;

    const SbVec3f &point = geom->point[static_cast<int>(index)];
    out[0] = point[0];
    out[1] = point[1];
    out[2] = point[2];
    return 1;
}

extern "C" int
ged_draw_obol_database_source_surface_index_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    int *out)
{
    if (!out)
	return 0;

    SoBRLMeshShape *shape = ged_obol_owned_mesh_shape_for_path(gedp, path, 0);
    if (!shape)
	return 0;

    const SoBRLMeshShape *geom = shape->getGeometrySource();
    const size_t triangle_count =
	static_cast<size_t>(shape->getTriangleCount());
    if (index >= triangle_count * 4)
	return 0;

    const size_t face_offset = index % 4;
    if (face_offset == 3) {
	*out = -1;
	return 1;
    }

    const size_t coord_index = (index / 4) * 3 + face_offset;
    if (coord_index >= static_cast<size_t>(geom->coordIndex.getNum()))
	return 0;
    *out = geom->coordIndex[static_cast<int>(coord_index)];
    return 1;
}

extern "C" int
ged_draw_obol_database_source_translate_vlist_for_path(
    struct ged *gedp,
    const char *path,
    const vect_t xlate)
{
    if (!xlate)
	return 0;

    SoBRLVListShape *shape =
	ged_obol_owned_annotation_vlist_shape_for_path(gedp, path);
    if (!shape)
	return 0;

    return shape->translatePoints(SbVec3f(
				      static_cast<float>(xlate[0]),
				      static_cast<float>(xlate[1]),
				      static_cast<float>(xlate[2]))) ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_clear_vlist_for_path(
    struct ged *gedp,
    const char *path)
{
    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    const int cleared =
	scene->clearDatabaseSourceInstanceExternalPrimaryGeometry(
	    source_instance_key.c_str());
    return cleared >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_publish_primitive_wireframe_for_path(
    struct ged *gedp,
    const char *path,
    struct rt_db_internal *ip,
    const struct bg_tess_tol *ttol,
    const struct bn_tol *tol)
{
    if (!gedp || !path || !path[0] || !ip)
	return 0;

    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    return scene->publishDatabaseSourceInstancePrimitiveWireframe(
	       source_instance_key.c_str(), ip, ttol, tol);
}

extern "C" int
ged_draw_obol_database_source_publish_line_set_for_path(
    struct ged *gedp,
    const char *path,
    const point_t *points,
    const int *commands,
    size_t point_count)
{
    if (point_count > static_cast<size_t>(INT_MAX))
	return 0;
    if (point_count && !points)
	return 0;

    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    std::vector<SbVec3f> obol_points;
    std::vector<int32_t> obol_commands;
    std::vector<double> precise_points;
    obol_points.reserve(point_count);
    obol_commands.reserve(point_count);
    precise_points.reserve(point_count * 3);
    for (size_t i = 0; i < point_count; i++) {
	const int command = commands ? commands[i] : -1;
	const int32_t obol_command =
	    ged_obol_vlist_command_from_ged(command, i);
	if (obol_command < 0)
	    return 0;

	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
	obol_commands.push_back(obol_command);
	precise_points.push_back(points[i][0]);
	precise_points.push_back(points[i][1]);
	precise_points.push_back(points[i][2]);
    }

    BObolExternalLineSet line_set;
    line_set.points = obol_points.empty() ? NULL : obol_points.data();
    line_set.commands = obol_commands.empty() ? NULL : obol_commands.data();
    line_set.precisePoints = precise_points.empty() ? NULL :
			     precise_points.data();
    line_set.count = static_cast<int>(point_count);
    const int published =
	scene->publishDatabaseSourceInstanceExternalLineSet(
	    source_instance_key.c_str(), line_set);
    return published > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_publish_annotation_line_set_for_path(
    struct ged *gedp,
    const char *path,
    const point_t *points,
    const int *commands,
    size_t point_count)
{
    if (point_count > static_cast<size_t>(INT_MAX))
	return 0;
    if (point_count && !points)
	return 0;

    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    std::vector<SbVec3f> obol_points;
    std::vector<int32_t> obol_commands;
    std::vector<double> precise_points;
    obol_points.reserve(point_count);
    obol_commands.reserve(point_count);
    precise_points.reserve(point_count * 3);
    for (size_t i = 0; i < point_count; i++) {
	const int command = commands ? commands[i] : -1;
	const int32_t obol_command =
	    ged_obol_vlist_command_from_ged(command, i);
	if (obol_command < 0)
	    return 0;

	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
	obol_commands.push_back(obol_command);
	precise_points.push_back(points[i][0]);
	precise_points.push_back(points[i][1]);
	precise_points.push_back(points[i][2]);
    }

    BObolExternalLineSet line_set;
    line_set.points = obol_points.empty() ? NULL : obol_points.data();
    line_set.commands = obol_commands.empty() ? NULL : obol_commands.data();
    line_set.precisePoints = precise_points.empty() ? NULL :
			     precise_points.data();
    line_set.count = static_cast<int>(point_count);
    line_set.sourceType = "annotation";
    line_set.geometryKind = "annotation";
    const int published =
	scene->publishDatabaseSourceInstanceExternalLineSet(
	    source_instance_key.c_str(), line_set);
    return published > 0 ? 1 : 0;
}

extern "C" int
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
    size_t line_point_count)
{
    if ((annotation_point_count && !annotation_points) ||
	(segment_count && !segments) ||
	segment_count > static_cast<size_t>(INT_MAX) ||
	annotation_point_count > static_cast<size_t>(INT_MAX) ||
	line_point_count > static_cast<size_t>(INT_MAX) ||
	(line_point_count && !line_points))
	return 0;

    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    std::vector<SbVec3f> obol_line_points;
    std::vector<int32_t> obol_line_commands;
    std::vector<double> precise_line_points;
    obol_line_points.reserve(line_point_count);
    obol_line_commands.reserve(line_point_count);
    precise_line_points.reserve(line_point_count * 3);
    for (size_t i = 0; i < line_point_count; i++) {
	const int command = line_commands ? line_commands[i] : -1;
	const int32_t obol_command =
	    ged_obol_vlist_command_from_ged(command, i);
	if (obol_command < 0)
	    return 0;

	obol_line_points.push_back(SbVec3f(
				       static_cast<float>(line_points[i][0]),
				       static_cast<float>(line_points[i][1]),
				       static_cast<float>(line_points[i][2])));
	obol_line_commands.push_back(obol_command);
	precise_line_points.push_back(line_points[i][0]);
	precise_line_points.push_back(line_points[i][1]);
	precise_line_points.push_back(line_points[i][2]);
    }

    std::vector<SbVec3f> obol_annotation_points;
    std::vector<double> precise_annotation_points;
    obol_annotation_points.reserve(annotation_point_count);
    precise_annotation_points.reserve(annotation_point_count * 3);
    if (annotation_point_count && annotation_points) {
	for (size_t i = 0; i < annotation_point_count; i++) {
	    obol_annotation_points.push_back(SbVec3f(
						 static_cast<float>(annotation_points[i][0]),
						 static_cast<float>(annotation_points[i][1]),
						 static_cast<float>(annotation_points[i][2])));
	    precise_annotation_points.push_back(annotation_points[i][0]);
	    precise_annotation_points.push_back(annotation_points[i][1]);
	    precise_annotation_points.push_back(annotation_points[i][2]);
	}
    }

    std::vector<BObolExternalAnnotationSegment> obol_segments;
    obol_segments.reserve(segment_count);
    for (size_t i = 0; i < segment_count; i++) {
	const struct ged_draw_obol_annotation_segment *seg = &segments[i];
	BObolExternalAnnotationSegment obol_seg;
	if (seg->kind == GED_DRAW_OBOL_ANNOTATION_SEGMENT_LINE)
	    obol_seg.kind = BObolExternalAnnotationSegment::SEGMENT_LINE;
	else if (seg->kind == GED_DRAW_OBOL_ANNOTATION_SEGMENT_TEXT)
	    obol_seg.kind = BObolExternalAnnotationSegment::SEGMENT_TEXT;
	else
	    obol_seg.kind = BObolExternalAnnotationSegment::SEGMENT_NONE;
	obol_seg.lineStart = seg->line_start;
	obol_seg.lineEnd = seg->line_end;
	obol_seg.textRefPoint = seg->text_ref_point;
	obol_seg.text = seg->text;
	obol_segments.push_back(obol_seg);
    }

    BObolExternalAnnotation annotation;
    annotation.basePoint = base_point ? SbVec3f(
			       static_cast<float>(base_point[0]),
			       static_cast<float>(base_point[1]),
			       static_cast<float>(base_point[2])) :
			   SbVec3f(0.0f, 0.0f, 0.0f);
    annotation.linePoints = obol_line_points.empty() ? NULL :
			    obol_line_points.data();
    annotation.lineCommands = obol_line_commands.empty() ? NULL :
			      obol_line_commands.data();
    annotation.preciseLinePoints = precise_line_points.empty() ? NULL :
				   precise_line_points.data();
    annotation.linePointCount = static_cast<int>(line_point_count);
    annotation.annotationPoints = obol_annotation_points.empty() ? NULL :
				  obol_annotation_points.data();
    annotation.preciseAnnotationPoints =
	precise_annotation_points.empty() ? NULL :
	precise_annotation_points.data();
    annotation.annotationPointCount =
	static_cast<int>(annotation_point_count);
    annotation.segments = obol_segments.empty() ? NULL : obol_segments.data();
    annotation.segmentCount = static_cast<int>(segment_count);
    const int published =
	scene->publishDatabaseSourceInstanceExternalAnnotation(
	    source_instance_key.c_str(), annotation);
    return published > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_annotation_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_obol_annotation_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    SoBRLVListShape *shape =
	ged_obol_owned_annotation_vlist_shape_for_path(gedp, path);
    if (!shape)
	return 0;

    const char *kind = shape->geometryKind.getValue().getString();
    if (!kind || !BU_STR_EQUAL(kind, "annotation"))
	return 0;

    const SoBRLVListShape *geom = shape->getGeometrySource();
    out->valid = 1;
    out->point_count = static_cast<size_t>(geom->annotationPoint.getNum());
    out->segment_count =
	static_cast<size_t>(geom->annotationSegmentKind.getNum());
    out->cache_identity =
	ged_obol_hash_sb_string(shape->cacheIdentity.getValue());
    out->source_identity =
	ged_obol_hash_sb_string(shape->sourceIdentity.getValue());

    for (int i = 0; i < geom->annotationSegmentKind.getNum(); i++) {
	const int kind_value = geom->annotationSegmentKind[i];
	if (!out->line_segment_valid &&
	    kind_value == GED_DRAW_OBOL_ANNOTATION_SEGMENT_LINE) {
	    out->line_segment_valid = 1;
	    out->line_start = (i < geom->annotationSegmentStart.getNum()) ?
			      geom->annotationSegmentStart[i] : 0;
	    out->line_end = (i < geom->annotationSegmentEnd.getNum()) ?
			    geom->annotationSegmentEnd[i] : 0;
	}
	if (!out->text_segment_valid &&
	    kind_value == GED_DRAW_OBOL_ANNOTATION_SEGMENT_TEXT) {
	    out->text_segment_valid = 1;
	    out->text_ref_point = (i < geom->annotationTextRefPoint.getNum()) ?
				  geom->annotationTextRefPoint[i] : 0;
	    out->text = (i < geom->annotationText.getNum()) ?
			geom->annotationText[i].getString() : "";
	}
    }

    return 1;
}

extern "C" int
ged_draw_obol_database_source_annotation_point_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    point_t out)
{
    if (!out)
	return 0;

    SoBRLVListShape *shape =
	ged_obol_owned_annotation_vlist_shape_for_path(gedp, path);
    const SoBRLVListShape *geom = shape ? shape->getGeometrySource() : NULL;
    if (!geom || index >= static_cast<size_t>(geom->annotationPoint.getNum()))
	return 0;

    double precise[3];
    if (shape->getPreciseAnnotationPoint(static_cast<int>(index), precise)) {
	out[0] = precise[0];
	out[1] = precise[1];
	out[2] = precise[2];
	return 1;
    }

    const SbVec3f &point = geom->annotationPoint[static_cast<int>(index)];
    out[0] = point[0];
    out[1] = point[1];
    out[2] = point[2];
    return 1;
}


extern "C" int
ged_draw_obol_database_source_publish_auxiliary_line_set_for_path(
    struct ged *gedp,
    const char *path,
    const char *name,
    const point_t *points,
    const int *commands,
    size_t point_count,
    const struct ged_draw_scene_display_summary *display_state)
{
    if (!name || !name[0] || point_count > static_cast<size_t>(INT_MAX))
	return 0;
    if (point_count && !points)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    std::vector<SbVec3f> obol_points;
    std::vector<int32_t> obol_commands;
    obol_points.reserve(point_count);
    obol_commands.reserve(point_count);
    for (size_t i = 0; i < point_count; i++) {
	const int command = commands ? commands[i] : -1;
	const int32_t obol_command =
	    ged_obol_vlist_command_from_ged(command, i);
	if (obol_command < 0)
	    return 0;

	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
	obol_commands.push_back(obol_command);
    }

    BObolAuxiliaryLineSetDisplayState obol_display;
    const BObolAuxiliaryLineSetDisplayState *obol_display_ptr = NULL;
    if (display_state && display_state->valid) {
	obol_display.valid = TRUE;
	obol_display.drawMode =
	    ged_obol_lod_draw_mode_from_ged(display_state->draw_mode);
	obol_display.visible = display_state->visible ? TRUE : FALSE;
	obol_display.highlighted = display_state->highlighted ? TRUE : FALSE;
	obol_display.lineStyle = display_state->line_style;
	obol_display.lineWidth = display_state->line_width;
	obol_display.transparency =
	    static_cast<float>(display_state->transparency);
	obol_display.materialColorValid =
	    display_state->material_valid ? TRUE : FALSE;
	if (display_state->material_valid)
	    obol_display.materialColor =
		ged_obol_color_from_rgb(display_state->material_color);
	obol_display.materialRevision = 0;
	obol_display_ptr = &obol_display;
    }

    int changed_any = 0;
    for (const std::string &source_instance_key : instance_keys) {
	const int changed =
	    scene->publishDatabaseSourceInstanceAuxiliaryLineSet(
		source_instance_key.c_str(),
		name,
		obol_points.empty() ? NULL : obol_points.data(),
		obol_commands.empty() ? NULL : obol_commands.data(),
		static_cast<int>(point_count),
		obol_display_ptr);
	if (changed < 0)
	    return 0;
	if (changed > 0)
	    changed_any = 1;
    }
    return changed_any;
}

extern "C" int
ged_draw_obol_database_source_publish_auxiliary_source_line_set_for_path(
    struct ged *gedp,
    const char *path,
    const char *source_path,
    const char *display_name,
    const point_t *points,
    const int *commands,
    size_t point_count,
    const struct ged_draw_scene_display_summary *display_state)
{
    if (!source_path || !source_path[0] ||
	point_count > static_cast<size_t>(INT_MAX))
	return 0;
    if (point_count && !points)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> owner_instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (owner_instance_keys.empty())
	return 0;

    std::vector<SbVec3f> obol_points;
    std::vector<int32_t> obol_commands;
    obol_points.reserve(point_count);
    obol_commands.reserve(point_count);
    for (size_t i = 0; i < point_count; i++) {
	const int command = commands ? commands[i] : -1;
	const int32_t obol_command =
	    ged_obol_vlist_command_from_ged(command, i);
	if (obol_command < 0)
	    return 0;

	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
	obol_commands.push_back(obol_command);
    }

    BObolAuxiliaryLineSetDisplayState obol_display;
    const BObolAuxiliaryLineSetDisplayState *obol_display_ptr = NULL;
    if (display_state && display_state->valid) {
	obol_display.valid = TRUE;
	obol_display.drawMode =
	    ged_obol_lod_draw_mode_from_ged(display_state->draw_mode);
	obol_display.visible = display_state->visible ? TRUE : FALSE;
	obol_display.highlighted = display_state->highlighted ? TRUE : FALSE;
	obol_display.lineStyle = display_state->line_style;
	obol_display.lineWidth = display_state->line_width;
	obol_display.transparency =
	    static_cast<float>(display_state->transparency);
	obol_display.materialColorValid =
	    display_state->material_valid ? TRUE : FALSE;
	if (display_state->material_valid)
	    obol_display.materialColor =
		ged_obol_color_from_rgb(display_state->material_color);
	obol_display.materialRevision = 0;
	obol_display_ptr = &obol_display;
    }

    int changed_any = 0;
    for (const std::string &owner_instance_key : owner_instance_keys) {
	const int changed =
	    scene->publishDatabaseSourceInstanceAuxiliarySourceLineSet(
		owner_instance_key.c_str(),
		source_path,
		display_name ? display_name : source_path,
		obol_points.empty() ? NULL : obol_points.data(),
		obol_commands.empty() ? NULL : obol_commands.data(),
		static_cast<int>(point_count),
		obol_display_ptr);
	if (changed < 0)
	    return 0;
	if (changed > 0)
	    changed_any = 1;
    }
    if (!changed_any)
	return 0;

    SoBRLDatabaseSource *source = scene->findDatabaseSource(source_path);
    SoBRLVListShape *shape = source ? source->getRealizedShape(0) : NULL;
    if (shape)
	ged_obol_vlist_shape_set_precise_points(shape, points, point_count);
    return 1;
}

extern "C" int
ged_draw_obol_database_source_clear_auxiliary_shapes_for_path(
    struct ged *gedp,
    const char *path)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	const int cleared = scene->clearDatabaseSourceInstanceAuxiliaryShapes(
				source_instance_key.c_str());
	if (cleared < 0)
	    return 0;
	if (cleared > 0)
	    applied = 1;
    }
    return applied;
}

extern "C" int
ged_draw_obol_database_source_publish_point_set_for_path(
    struct ged *gedp,
    const char *path,
    const point_t *points,
    size_t point_count)
{
    if (point_count > static_cast<size_t>(INT_MAX))
	return 0;
    if (point_count && !points)
	return 0;

    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    std::vector<SbVec3f> obol_points;
    std::vector<double> precise_points;
    obol_points.reserve(point_count);
    precise_points.reserve(point_count * 3);
    for (size_t i = 0; i < point_count; i++) {
	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
	precise_points.push_back(points[i][0]);
	precise_points.push_back(points[i][1]);
	precise_points.push_back(points[i][2]);
    }

    BObolExternalPointSet point_set;
    point_set.points = obol_points.empty() ? NULL : obol_points.data();
    point_set.precisePoints = precise_points.empty() ? NULL :
			      precise_points.data();
    point_set.count = static_cast<int>(point_count);
    const int published =
	scene->publishDatabaseSourceInstanceExternalPointSet(
	    source_instance_key.c_str(), point_set);
    return published > 0 ? 1 : 0;
}

static int
ged_obol_indexed_face_finish(std::vector<int32_t> &face,
			     size_t point_count,
			     std::vector<int32_t> &triangles,
			     size_t *face_count,
			     unsigned int *face_stamp,
			     std::vector<unsigned int> &seen)
{
    if (face.size() < 3)
	return 0;

    for (size_t i = 1; i + 1 < face.size(); i++) {
	triangles.push_back(face[0]);
	triangles.push_back(face[i]);
	triangles.push_back(face[i + 1]);
    }

    face.clear();
    if (face_count)
	(*face_count)++;
    if (face_stamp && seen.size() == point_count) {
	if (*face_stamp == UINT_MAX) {
	    std::fill(seen.begin(), seen.end(), 0);
	    *face_stamp = 1;
	} else {
	    (*face_stamp)++;
	}
    }
    return 1;
}

static int
ged_obol_indexed_faces_to_triangles(const int *indices,
				    size_t index_count,
				    size_t point_count,
				    std::vector<int32_t> &triangles,
				    size_t *face_count_out,
				    size_t *vertex_index_count_out)
{
    if (!indices || !index_count || !point_count ||
	point_count > static_cast<size_t>(INT_MAX) ||
	index_count > static_cast<size_t>(INT_MAX))
	return 0;

    size_t face_count = 0;
    size_t vertex_index_count = 0;
    unsigned int face_stamp = 1;
    std::vector<unsigned int> seen(point_count, 0);
    std::vector<int32_t> face;

    for (size_t i = 0; i < index_count; i++) {
	const int idx = indices[i];
	if (idx < 0) {
	    if (idx != -1 || !ged_obol_indexed_face_finish(face,
		    point_count, triangles, &face_count, &face_stamp, seen))
		return 0;
	    continue;
	}

	if (static_cast<size_t>(idx) >= point_count)
	    return 0;
	if (seen[static_cast<size_t>(idx)] == face_stamp)
	    return 0;
	seen[static_cast<size_t>(idx)] = face_stamp;
	vertex_index_count++;
	face.push_back(static_cast<int32_t>(idx));
    }

    if (!face.empty() && !ged_obol_indexed_face_finish(face,
	    point_count, triangles, &face_count, &face_stamp, seen))
	return 0;
    if (!face_count || triangles.empty())
	return 0;

    if (face_count_out)
	*face_count_out = face_count;
    if (vertex_index_count_out)
	*vertex_index_count_out = vertex_index_count;
    return 1;
}

static int
ged_obol_indexed_face_triangle_normals(const vect_t *normals,
	size_t normal_count, const int *indices, size_t index_count,
	size_t point_count, size_t face_count, size_t vertex_index_count,
	const std::vector<int32_t> &triangles,
	std::vector<SbVec3f> &triangle_normals)
{
    triangle_normals.clear();
    if (!normal_count)
	return 1;
    if (!normals || !indices || !index_count)
	return 0;

    enum normal_binding {
	NORMAL_PER_CORNER,
	NORMAL_PER_RAW_INDEX,
	NORMAL_PER_POINT,
	NORMAL_PER_FACE
    };
    normal_binding binding;
    if (normal_count == vertex_index_count)
	binding = NORMAL_PER_CORNER;
    else if (normal_count == index_count)
	binding = NORMAL_PER_RAW_INDEX;
    else if (normal_count == point_count)
	binding = NORMAL_PER_POINT;
    else if (normal_count == face_count)
	binding = NORMAL_PER_FACE;
    else
	return 0;

    std::vector<int> face_indices;
    std::vector<SbVec3f> face_normals;
    size_t corner_index = 0;
    size_t face_index = 0;
    auto append_face = [&]() -> int {
	if (face_indices.size() < 3 || face_normals.size() != face_indices.size())
	    return 0;
	for (size_t i = 1; i + 1 < face_indices.size(); i++) {
	    triangle_normals.push_back(face_normals[0]);
	    triangle_normals.push_back(face_normals[i]);
	    triangle_normals.push_back(face_normals[i + 1]);
	}
	face_indices.clear();
	face_normals.clear();
	face_index++;
	return 1;
    };

    for (size_t i = 0; i < index_count; i++) {
	const int point_index = indices[i];
	if (point_index < 0) {
	    if (point_index != -1 || !append_face())
		return 0;
	    continue;
	}
	if (static_cast<size_t>(point_index) >= point_count)
	    return 0;
	size_t normal_index = 0;
	switch (binding) {
	    case NORMAL_PER_CORNER:
		normal_index = corner_index;
		break;
	    case NORMAL_PER_RAW_INDEX:
		normal_index = i;
		break;
	    case NORMAL_PER_POINT:
		normal_index = static_cast<size_t>(point_index);
		break;
	    case NORMAL_PER_FACE:
		normal_index = face_index;
		break;
	}
	if (normal_index >= normal_count)
	    return 0;
	face_indices.push_back(point_index);
	face_normals.push_back(SbVec3f(
	    static_cast<float>(normals[normal_index][X]),
	    static_cast<float>(normals[normal_index][Y]),
	    static_cast<float>(normals[normal_index][Z])));
	corner_index++;
    }
    if (!face_indices.empty() && !append_face())
	return 0;
    return face_index == face_count &&
	triangle_normals.size() == triangles.size();
}

extern "C" int
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
    const struct ged_draw_scene_display_summary *display_state)
{
    if (!group_path || !group_path[0] || !shape_path || !shape_path[0])
	return 0;
    if (!point_count || !index_count)
	return 0;
    if (!points || !indices || (normal_count && !normals) ||
	point_count > static_cast<size_t>(INT_MAX))
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<int32_t> triangles;
    size_t face_count = 0;
    size_t vertex_index_count = 0;
    if (!ged_obol_indexed_faces_to_triangles(indices, index_count,
	    point_count, triangles, &face_count, &vertex_index_count))
	return 0;
    if (normal_count && normal_count != vertex_index_count &&
	normal_count != index_count &&
	normal_count != point_count && normal_count != face_count)
	return 0;
    if (triangles.size() > static_cast<size_t>(INT_MAX))
	return 0;

    SoBRLMeshShape *shape =
	ged_obol_local_mesh_shape_for_path(scene, group_path, shape_path);
    if (!shape)
	return 0;

    std::vector<SbVec3f> obol_points;
    obol_points.reserve(point_count);
    for (size_t i = 0; i < point_count; i++) {
	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
    }

    std::vector<SbVec3f> triangle_normals;
    if (!ged_obol_indexed_face_triangle_normals(normals, normal_count,
	indices, index_count, point_count, face_count, vertex_index_count,
	triangles, triangle_normals))
	return 0;

    ged_obol_local_shape_apply_common_state(shape, shape_path, display_name,
					    "local-indexed-face-set", "surface", display_state);
    shape->setIndexedTriangles(obol_points.data(),
			       static_cast<int>(obol_points.size()),
			       triangles.data(),
			       static_cast<int>(triangles.size()),
			       triangle_normals.empty() ? NULL : triangle_normals.data(),
			       static_cast<int>(triangle_normals.size()));
    return 1;
}

extern "C" int
ged_draw_obol_database_source_clear_mesh_for_path(
    struct ged *gedp,
    const char *path)
{
    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    const int cleared =
	scene->clearDatabaseSourceInstanceExternalPrimaryGeometry(
	    source_instance_key.c_str());
    return cleared >= 0 ? 1 : 0;
}

static int
ged_obol_database_source_publish_indexed_face_set_for_path(
    struct ged_bobol_publication_context *publication,
    struct ged *gedp,
    const char *path,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count,
    int lod_backed)
{
    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key, publication))
	return 0;

    if (!point_count || !index_count) {
	BObolExternalTriangleMesh empty_mesh;
	const int published =
	    scene->publishDatabaseSourceInstanceExternalTriangleMesh(
		source_instance_key.c_str(), empty_mesh);
	return published > 0 ? 1 : 0;
    }

    if (!points || !indices || (normal_count && !normals) ||
	point_count > static_cast<size_t>(INT_MAX))
	return 0;

    std::vector<int32_t> triangles;
    size_t face_count = 0;
    size_t vertex_index_count = 0;
    if (!ged_obol_indexed_faces_to_triangles(indices, index_count,
	    point_count, triangles, &face_count, &vertex_index_count))
	return 0;
    if (normal_count && normal_count != vertex_index_count &&
	normal_count != index_count &&
	normal_count != point_count && normal_count != face_count)
	return 0;
    if (triangles.size() > static_cast<size_t>(INT_MAX))
	return 0;

    std::vector<SbVec3f> obol_points;
    obol_points.reserve(point_count);
    for (size_t i = 0; i < point_count; i++) {
	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
    }

    BObolExternalTriangleMesh triangle_mesh;
    triangle_mesh.points = obol_points.empty() ? NULL : obol_points.data();
    triangle_mesh.pointCount = static_cast<int>(obol_points.size());
    triangle_mesh.indices = triangles.empty() ? NULL : triangles.data();
    triangle_mesh.indexCount = static_cast<int>(triangles.size());
    std::vector<SbVec3f> triangle_normals;
    if (!ged_obol_indexed_face_triangle_normals(normals, normal_count,
	indices, index_count, point_count, face_count, vertex_index_count,
	triangles, triangle_normals))
	return 0;
    triangle_mesh.normals = triangle_normals.empty() ? NULL :
	triangle_normals.data();
    triangle_mesh.normalCount = static_cast<int>(triangle_normals.size());
    triangle_mesh.lodBacked = lod_backed ? TRUE : FALSE;
    const int published =
	scene->publishDatabaseSourceInstanceExternalTriangleMesh(
	    source_instance_key.c_str(), triangle_mesh);
    return published > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_publish_indexed_face_set_for_path(
    struct ged *gedp,
    const char *path,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count)
{
	return ged_obol_database_source_publish_indexed_face_set_for_path(NULL,
	    gedp, path, points, point_count, normals, normal_count, indices,
	    index_count, 0);
}

int
ged_bobol_database_source_publish_indexed_face_set_for_path(
    struct ged_bobol_publication_context *publication,
    const char *path,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count)
{
    if (!publication || !publication->active || !publication->gedp)
	return 0;
    return ged_obol_database_source_publish_indexed_face_set_for_path(
	publication, publication->gedp, path, points, point_count, normals,
	normal_count, indices, index_count, 0);
}

extern "C" int
ged_draw_obol_database_source_publish_lod_indexed_face_set_for_path(
    struct ged *gedp,
    const char *path,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count)
{
	return ged_obol_database_source_publish_indexed_face_set_for_path(NULL,
	    gedp, path, points, point_count, normals, normal_count, indices,
	    index_count, 1);
}

extern "C" int
ged_draw_obol_database_source_set_vlist_center_for_path(
    struct ged *gedp,
    const char *path,
    const point_t center)
{
    if (!center)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *source =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	if (!source)
	    return 0;

	BObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid)
	    return 0;

	const int changed = scene->setDatabaseSourceInstancePlacementState(
				source_instance_key.c_str(),
				summary.drawMatrixValid,
				summary.drawMatrix,
				TRUE,
				SbVec3f(
				    static_cast<float>(center[0]),
				    static_cast<float>(center[1]),
				    static_cast<float>(center[2])),
				summary.drawSizeValid,
				summary.drawSize);
	if (changed < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}

extern "C" int
ged_draw_obol_database_source_update_vlist_bounds_for_path(
    struct ged *gedp,
    const char *path)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *source =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	if (!source)
	    return 0;

	SoBRLVListShape *shape =
	    ged_obol_owned_vlist_shape_for_source(source, path, 0);
	if (!shape)
	    return 0;

	if (!shape->updateDrawBoundsFromPoints())
	    return 0;

	BObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid) {
	    applied = 1;
	    continue;
	}

	SbBool nextDrawMatrixValid = summary.drawMatrixValid ? TRUE : FALSE;
	SbMatrix nextDrawMatrix = summary.drawMatrix;
	if (!nextDrawMatrixValid && shape->drawCenterValid.getValue() &&
	    shape->drawSizeValid.getValue()) {
	    nextDrawMatrixValid = TRUE;
	    nextDrawMatrix = SbMatrix::identity();
	}

	(void)scene->setDatabaseSourceInstancePlacementState(
	    source_instance_key.c_str(),
	    nextDrawMatrixValid,
	    nextDrawMatrix,
	    shape->drawCenterValid.getValue(),
	    shape->drawCenter.getValue(),
	    shape->drawSizeValid.getValue(),
	    shape->drawSize.getValue());
	applied = 1;
    }
    return applied;
}

extern "C" int
ged_draw_obol_database_source_set_placement_for_path(
    struct ged *gedp,
    const char *path,
    int draw_mat_valid,
    const mat_t draw_mat,
    int draw_center_valid,
    const point_t draw_center,
    int draw_size_valid,
    fastf_t draw_size)
{
    if ((draw_mat_valid && !draw_mat) ||
	(draw_center_valid && !draw_center))
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *source =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	if (!source)
	    return 0;

	BObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid)
	    return 0;

	SbBool nextDrawMatrixValid = draw_mat_valid ?
				     TRUE : summary.drawMatrixValid;
	SbMatrix nextDrawMatrix = draw_mat_valid ?
				  ged_obol_sbmatrix_from_mat(draw_mat) :
				  summary.drawMatrix;
	SbBool nextDrawCenterValid = draw_center_valid ?
				     TRUE : summary.drawCenterValid;
	SbVec3f nextDrawCenter = draw_center_valid ?
				 SbVec3f(static_cast<float>(draw_center[0]),
					 static_cast<float>(draw_center[1]),
					 static_cast<float>(draw_center[2])) :
				 summary.drawCenter;
	SbBool nextDrawSizeValid = draw_size_valid ?
				   TRUE : summary.drawSizeValid;
	float nextDrawSize = draw_size_valid ?
			     static_cast<float>(draw_size) : summary.drawSize;

	const int changed = scene->setDatabaseSourceInstancePlacementState(
				source_instance_key.c_str(),
				nextDrawMatrixValid, nextDrawMatrix,
				nextDrawCenterValid,
				nextDrawCenter, nextDrawSizeValid,
				nextDrawSize);
	if (changed < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}
