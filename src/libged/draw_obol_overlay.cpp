/*              D R A W _ O B O L _ O V E R L A Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file draw_obol_overlay.cpp
 *
 * Obol overlay publication and GED/Inventor value-boundary conversions.
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include "BObol/BMeshShape.h"
#include "BObol/BSceneController.h"
#include "BObol/BVListShape.h"
#include "BObol/BViewController.h"
#include "bu/malloc.h"
#include "ged.h"
#include "vmath.h"

#include "./draw_obol_bridge_private.hpp"
#include "./draw_obol_overlay_private.hpp"
#include "./ged_bobol_private.hpp"
#include "./ged_draw_private.h"

#include <Inventor/SbMatrix.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/nodes/SoNode.h>
#include <math.h>
#include <stdint.h>
#include <string>
#include <vector>

static unsigned char
ged_obol_rgb_channel(float value)
{
    if (value <= 0.0f)
	return 0;
    if (value >= 1.0f)
	return 255;
    return static_cast<unsigned char>(value * 255.0f + 0.5f);
}

void
ged_obol_rgb_from_color(const SbColor &color, unsigned char rgb[3])
{
    rgb[0] = ged_obol_rgb_channel(color.getValue()[0]);
    rgb[1] = ged_obol_rgb_channel(color.getValue()[1]);
    rgb[2] = ged_obol_rgb_channel(color.getValue()[2]);
}

void
ged_obol_mat_from_sbmatrix(const SbMatrix &matrix, mat_t mat)
{
    const SbMat &m = matrix.getValue();
    mat[0] = m[0][0];
    mat[4] = m[0][1];
    mat[8] = m[0][2];
    mat[12] = m[0][3];
    mat[1] = m[1][0];
    mat[5] = m[1][1];
    mat[9] = m[1][2];
    mat[13] = m[1][3];
    mat[2] = m[2][0];
    mat[6] = m[2][1];
    mat[10] = m[2][2];
    mat[14] = m[2][3];
    mat[3] = m[3][0];
    mat[7] = m[3][1];
    mat[11] = m[3][2];
    mat[15] = m[3][3];
}

SbMatrix
ged_obol_sbmatrix_from_mat(const mat_t mat)
{
    return SbMatrix(
	       static_cast<float>(mat[0]), static_cast<float>(mat[4]),
	       static_cast<float>(mat[8]), static_cast<float>(mat[12]),
	       static_cast<float>(mat[1]), static_cast<float>(mat[5]),
	       static_cast<float>(mat[9]), static_cast<float>(mat[13]),
	       static_cast<float>(mat[2]), static_cast<float>(mat[6]),
	       static_cast<float>(mat[10]), static_cast<float>(mat[14]),
	       static_cast<float>(mat[3]), static_cast<float>(mat[7]),
	       static_cast<float>(mat[11]), static_cast<float>(mat[15]));
}

SbColor
ged_obol_color_from_rgb(const unsigned char rgb[3])
{
    return SbColor(
	       static_cast<float>(rgb[0]) / 255.0f,
	       static_cast<float>(rgb[1]) / 255.0f,
	       static_cast<float>(rgb[2]) / 255.0f);
}

/* ged_draw_appearance_settings predates the retained scene and calls its
 * opacity field "transparency".  Keep that command-facing contract at this
 * boundary; Obol display state uses conventional transparency. */
float
ged_obol_transparency_from_appearance_opacity(fastf_t opacity)
{
    if (opacity <= 0.0)
	return 1.0f;
    if (opacity >= 1.0)
	return 0.0f;
    return 1.0f - static_cast<float>(opacity);
}

fastf_t
ged_obol_appearance_opacity_from_transparency(float transparency)
{
    if (transparency <= 0.0f)
	return 1.0;
    if (transparency >= 1.0f)
	return 0.0;
    return static_cast<fastf_t>(1.0f - transparency);
}

fastf_t
ged_obol_reported_transparency(float transparency)
{
    const double scaled = (double)transparency * 1000000.0;
    return (fastf_t)(floor(scaled + 0.5) / 1000000.0);
}

static BObolViewController *
ged_obol_view_controller_for_context(struct ged_view_context *view_ctx)
{
    return ged_bobol_view_controller(view_ctx);
}

static BObolViewController *
ged_obol_view_controller_ensure_for_context(struct ged_view_context *view_ctx,
	int sync_current_scene)
{
    struct ged *gedp = view_ctx ?
			   ged_view_context_owner(view_ctx) :
		       NULL;
    if (!gedp)
	return NULL;
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (controller)
	return controller;
    if (!ged_draw_obol_render_endpoint_ensure_for_view(gedp, view_ctx,
	    sync_current_scene))
	return NULL;
    return ged_bobol_view_controller(view_ctx);
}

static BObolViewController *
ged_obol_shared_view_controller_ensure_for_context(struct ged_view_context *view_ctx,
	int sync_current_scene)
{
    struct ged *gedp = view_ctx ?
			   ged_view_context_owner(view_ctx) :
		       NULL;
    if (!gedp)
	return NULL;
    BObolViewController *controller = ged_draw_obol_controller(gedp);
    if (controller)
	return controller;
    if (!ged_draw_obol_scene_controller_ensure_owned(gedp, sync_current_scene))
	return NULL;
    return ged_draw_obol_controller(gedp);
}

BObolViewController *
ged_obol_view_controller_for_scope(struct ged_view_context *view_ctx,
				   int local,
				   int sync_current_scene)
{
    return local ?
	   ged_obol_view_controller_ensure_for_context(view_ctx,
		   sync_current_scene) :
	   ged_obol_shared_view_controller_ensure_for_context(view_ctx,
		   sync_current_scene);
}


static std::string
ged_obol_overlay_group_path_for_name(const char *name)
{
    std::string path("_overlays");
    if (name && name[0]) {
	path += "/";
	path += name;
    }
    return path;
}

static std::string
ged_obol_overlay_shape_path_for_name(const char *name)
{
    std::string path = ged_obol_overlay_group_path_for_name(name);
    path += "/geometry";
    return path;
}

static void
ged_obol_overlay_group_state_set(BObolSceneController *scene,
				 const char *path,
				 const unsigned char *rgb,
				 fastf_t transparency,
				 int draw_mode)
{
    if (!scene || !path || !path[0])
	return;

    const SbColor color = rgb ? ged_obol_color_from_rgb(rgb) :
			  SbColor(1.0f, 1.0f, 1.0f);
    const int obol_draw_mode = draw_mode ?
			       ged_obol_lod_draw_mode_from_ged(draw_mode) :
			       BOBOL_LOD_DRAW_DIAGNOSTIC;
    (void)scene->setGroupDrawIntent(path, path, obol_draw_mode,
				    BOBOL_LOD_DRAW_WIRE, TRUE, 0);
    (void)scene->setGroupDisplayState(path, TRUE, FALSE, FALSE, 0, 0,
				      static_cast<float>(transparency), rgb ? TRUE : FALSE, color,
				      rgb ? TRUE : FALSE, color, 0);
}

template <typename ShapeT>
static void
ged_obol_overlay_shape_common_set(ShapeT *shape,
				  const char *name,
				  const std::string &shape_path,
				  const unsigned char basecolor[3],
				  fastf_t transparency,
				  int draw_mode,
				  const char *source_type,
				  const char *geometry_kind)
{
    if (!shape)
	return;

    const char *display_name = name && name[0] ? name : "overlay";
    shape->sourcePath = shape_path.c_str();
    shape->sourceName = display_name;
    shape->sourceType = source_type;
    shape->displayName = display_name;
    shape->geometryName = display_name;
    shape->sourceIdentity = shape_path.c_str();
    shape->cacheIdentity = shape_path.c_str();
    shape->ownerSourcePath = "_overlays";
    shape->databaseIntent = FALSE;
    shape->overlayIntent = TRUE;
    shape->hudIntent = FALSE;
    shape->localSource = TRUE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = TRUE;
    shape->drawMode = draw_mode ?
		      ged_obol_lod_draw_mode_from_ged(draw_mode) :
		      BOBOL_LOD_DRAW_DIAGNOSTIC;
    shape->recordRole = "overlay";
    shape->geometryKind = geometry_kind;
    shape->visible = TRUE;
    shape->selectable = TRUE;
    shape->transparency = static_cast<float>(transparency);
    shape->colorOverride = TRUE;
    shape->color = ged_obol_color_from_rgb(basecolor);
    shape->materialColorValid = TRUE;
    shape->materialColor = ged_obol_color_from_rgb(basecolor);
}

static std::vector<SbVec3f>
ged_obol_overlay_points(const struct ged_draw_overlay_geometry *geometry)
{
    std::vector<SbVec3f> points;
    if (!geometry || !geometry->points || !geometry->point_count)
	return points;

    points.reserve(geometry->point_count);
    for (size_t i = 0; i < geometry->point_count; i++)
	points.push_back(SbVec3f(
			     static_cast<float>(geometry->points[i][X]),
			     static_cast<float>(geometry->points[i][Y]),
			     static_cast<float>(geometry->points[i][Z])));
    return points;
}

static SoBRLVListShape *
ged_obol_overlay_vlist_shape_create(const char *name,
				    const std::string &shape_path,
				    const struct ged_draw_overlay_geometry *geometry,
				    const unsigned char basecolor[3],
				    fastf_t transparency,
				    int draw_mode)
{
    if (!geometry)
	return NULL;

    std::vector<SbVec3f> points = ged_obol_overlay_points(geometry);
    if (points.empty())
	return NULL;

    std::vector<int32_t> commands;
    const char *source_type = "line-set";
    const char *geometry_kind = "line";
    if (geometry->kind == GED_DRAW_OVERLAY_GEOMETRY_POINT_SET) {
	commands.assign(points.size(),
			static_cast<int32_t>(BObolLineCommand::Point));
	source_type = "point-set";
	geometry_kind = "point";
    } else {
	if (geometry->kind != GED_DRAW_OVERLAY_GEOMETRY_LINE_SET ||
	    !geometry->commands ||
	    geometry->command_count != geometry->point_count)
	    return NULL;
	commands.reserve(geometry->command_count);
	for (size_t i = 0; i < geometry->command_count; i++)
	    commands.push_back(ged_obol_line_command_from_ged(
				   geometry->commands[i], i));
    }

    SoBRLVListShape *shape = new SoBRLVListShape;
    ged_obol_overlay_shape_common_set(shape, name, shape_path, basecolor,
				      transparency, draw_mode, source_type, geometry_kind);
    shape->setLineSet(points.data(), commands.data(),
		      static_cast<int>(points.size()));
    return shape;
}

static void
ged_obol_overlay_append_face_triangles(const std::vector<int32_t> &face,
				       size_t point_count,
				       std::vector<int32_t> &triangles)
{
    if (face.size() < 3)
	return;

    std::vector<int32_t> clean_face = face;
    if (clean_face.size() > 3 && clean_face.front() == clean_face.back())
	clean_face.pop_back();
    if (clean_face.size() < 3)
	return;

    for (size_t i = 0; i < clean_face.size(); i++) {
	if (clean_face[i] < 0 ||
	    static_cast<size_t>(clean_face[i]) >= point_count)
	    return;
    }

    for (size_t i = 1; i + 1 < clean_face.size(); i++) {
	triangles.push_back(clean_face[0]);
	triangles.push_back(clean_face[i]);
	triangles.push_back(clean_face[i + 1]);
    }
}

static std::vector<int32_t>
ged_obol_overlay_triangles(const struct ged_draw_overlay_geometry *geometry)
{
    std::vector<int32_t> triangles;
    if (!geometry || !geometry->indices || !geometry->index_count ||
	!geometry->point_count)
	return triangles;

    std::vector<int32_t> face;
    for (size_t i = 0; i < geometry->index_count; i++) {
	const int32_t idx = static_cast<int32_t>(geometry->indices[i]);
	if (idx < 0) {
	    ged_obol_overlay_append_face_triangles(face,
						   geometry->point_count, triangles);
	    face.clear();
	} else {
	    face.push_back(idx);
	}
    }
    ged_obol_overlay_append_face_triangles(face, geometry->point_count,
					   triangles);
    return triangles;
}

static SoBRLMeshShape *
ged_obol_overlay_mesh_shape_create(const char *name,
				   const std::string &shape_path,
				   const struct ged_draw_overlay_geometry *geometry,
				   const unsigned char basecolor[3],
				   fastf_t transparency,
				   int draw_mode)
{
    if (!geometry || geometry->kind != GED_DRAW_OVERLAY_GEOMETRY_INDEXED_FACE_SET)
	return NULL;

    std::vector<SbVec3f> points = ged_obol_overlay_points(geometry);
    std::vector<int32_t> triangles = ged_obol_overlay_triangles(geometry);
    if (points.empty() || triangles.empty())
	return NULL;

    SoBRLMeshShape *shape = new SoBRLMeshShape;
    ged_obol_overlay_shape_common_set(shape, name, shape_path, basecolor,
				      transparency, draw_mode, "indexed-face-set", "surface");
    shape->setIndexedTriangles(points.data(), static_cast<int>(points.size()),
			       triangles.data(), static_cast<int>(triangles.size()));
    return shape;
}

extern "C" int
ged_draw_obol_overlay_erase_name_context(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char *name)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name || !name[0])
	return 0;

    BObolSceneController *scene = controller->getSceneController();
    if (!scene)
	return 0;

    const std::string overlay_path = ged_obol_overlay_group_path_for_name(name);
    const int removed = scene->removeGroup(overlay_path.c_str());
    if (scene->getGroupChildCount("_overlays") == 0)
	(void)scene->removeGroup("_overlays");
    if (removed > 0) {
	scene->realizePending();
	controller->requestLodCapacityRender("overlay-remove");
    }
    (void)gedp;
    return removed > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_overlay_geometry_insert_context(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char *name,
    const struct ged_draw_overlay_geometry *geometry,
    const unsigned char basecolor[3],
    fastf_t transparency,
    int draw_mode,
    char **shape_path_out)
{
    if (shape_path_out)
	*shape_path_out = NULL;
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name || !name[0] || !geometry || !basecolor)
	return 0;

    BObolSceneController *scene = controller->getSceneController();
    if (!scene)
	return 0;

    const std::string overlay_path = ged_obol_overlay_group_path_for_name(name);
    const std::string shape_path = ged_obol_overlay_shape_path_for_name(name);
    (void)scene->removeGroup(overlay_path.c_str());
    if (!scene->ensureGroup("_overlays") ||
	!scene->ensureGroup(overlay_path.c_str()))
	return 0;

    ged_obol_overlay_group_state_set(scene, "_overlays", NULL, 0.0, 0);
    ged_obol_overlay_group_state_set(scene, overlay_path.c_str(), basecolor,
				     transparency, draw_mode);

    SoNode *shape = NULL;
    if (geometry->kind == GED_DRAW_OVERLAY_GEOMETRY_INDEXED_FACE_SET)
	shape = ged_obol_overlay_mesh_shape_create(name, shape_path, geometry,
		basecolor, transparency, draw_mode);
    else
	shape = ged_obol_overlay_vlist_shape_create(name, shape_path, geometry,
		basecolor, transparency, draw_mode);
    if (!shape) {
	(void)scene->removeGroup(overlay_path.c_str());
	if (scene->getGroupChildCount("_overlays") == 0)
	    (void)scene->removeGroup("_overlays");
	return 0;
    }

    if (scene->appendChildToGroup(overlay_path.c_str(), shape) < 0) {
	shape->ref();
	shape->unref();
	(void)scene->removeGroup(overlay_path.c_str());
	if (scene->getGroupChildCount("_overlays") == 0)
	    (void)scene->removeGroup("_overlays");
	return 0;
    }

    if (shape_path_out)
	*shape_path_out = bu_strdup(shape_path.c_str());
    controller->requestLodCapacityRender("overlay-insert");
    (void)gedp;
    return 1;
}

