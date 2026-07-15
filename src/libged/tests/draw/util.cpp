/*                       U T I L . C P P
 * BRL-CAD
 *
 * Copyright (c) 2018-2026 United States Government as represented by
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
/** @file util.cpp
 *
 * Utility testing for drawing routines
 *
 */

#include "common.h"

#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <stdio.h>
#include <inttypes.h>
#include <fstream>
#include <thread>

#include <bu.h>
#include <brlobol/init.h>
#include <brlobol/database_source.h>
#include <brlobol/mesh_lod_cache.h>
#include <brlobol/scene_controller.h>
#include <brlobol/vlist_shape.h>
#include <brlobol/view_controller.h>
#include <icv.h>
#include <Inventor/SbColor.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/SoDB.h>
#include <Inventor/SoOffscreenRenderer.h>
#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoMatrixTransform.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <rt/view.h>
#include "view_test_util.h"
#include <ged.h>
#include <ged/draw.h>
#include <ged/draw_obol.h>
#include <ged/event_txn.h>

// In order to handle changes to .g geometry contents, we need to defined
// callbacks for the librt hooks that will update the working data structures.
// In Qt we have libqtcad handle this, but as we are not using a QgModel we
// need to do it ourselves.
extern "C" void
ged_changed_callback(struct db_i *dbip, struct directory *dp, int mode, void *u_data)
{
    struct ged *gedp = (struct ged *)u_data;
    const char *name = (dp) ? dp->d_namep : NULL;

    // Need to invalidate any LoD caches associated with this dp
    if (dbip && dp && dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT)
	(void)brlobol_mesh_lod_cache_invalidate(dbip, dp->d_namep, NULL);

    switch (mode) {
	case 0:
	    ged_event_notify_object_modified(gedp, name, 1, NULL);
	    break;
	case 1:
	    ged_event_notify_object_added(gedp, name, NULL);
	    break;
	case 2:
	    ged_event_notify_object_removed(gedp, name, NULL);
	    break;
	default:
	    bu_log("changed callback mode error: %d\n", mode);
    }
}

static int
draw_test_obol_capture_enabled(void)
{
    const char *value = std::getenv("GED_TEST_DRAW_OBOL_CAPTURE");
    return value ? bu_str_true(value) : 1;
}

static int
draw_test_obol_debug_enabled(void)
{
    const char *value = std::getenv("GED_TEST_DRAW_OBOL_DEBUG");
    return value ? bu_str_true(value) : 0;
}

static int
draw_test_compare_debug_enabled(void)
{
    const char *value = std::getenv("GED_TEST_DRAW_COMPARE_DEBUG");
    return value ? bu_str_true(value) : 0;
}

static void *
draw_test_active_view_ctx(struct ged *gedp)
{
    if (!gedp)
	return NULL;

    struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
    if (views && BU_PTBL_LEN(views) > 0)
	return BU_PTBL_GET(views, 0);

    return ged_view_active_ctx(gedp);
}

static SoDB::ContextManager *
draw_test_obol_context_manager(void)
{
    static SoDB::ContextManager *manager = SoDB::createOSMesaContextManager();
    return manager;
}

static int
draw_test_sync_obol_camera(BRLObolViewController *controller, void *view_ctx)
{
    if (!controller || !view_ctx)
	return 0;

    return controller->syncCameraFromViewContext(view_ctx) ? 1 : 0;
}

extern "C" int
draw_test_obol_view_init(struct ged *gedp, void *view_ctx, int width, int height)
{
    if (!gedp || !view_ctx || width <= 0 || height <= 0)
	return BRLCAD_ERROR;

    bv_dimensions_set(DRAW_TEST_BV(view_ctx), width, height);
    if (gedp->dbip) {
	bv_unit_conversion_set(DRAW_TEST_BV(view_ctx), gedp->dbip->dbi_local2base,
		gedp->dbip->dbi_base2local);
    }
    if (!ged_draw_obol_render_endpoint_ensure_for_view(gedp, view_ctx, 1))
	return BRLCAD_ERROR;
    BRLObolViewController *controller =
	(BRLObolViewController *)ged_draw_obol_controller_opaque_for_view(view_ctx);
    if (!controller)
	return BRLCAD_ERROR;

    controller->setViewportSize((unsigned int)width, (unsigned int)height);
    return draw_test_sync_obol_camera(controller, view_ctx) ?
	BRLCAD_OK : BRLCAD_ERROR;
}

extern "C" int
draw_test_obol_progressive_drain(struct ged *gedp, void *view_ctx,
	unsigned int max_attempts, unsigned int sleep_milliseconds)
{
    if (!gedp || !view_ctx || !max_attempts)
	return 0;

    BRLObolViewController *controller =
	static_cast<BRLObolViewController *>(
	    ged_draw_obol_controller_opaque_for_view(view_ctx));
    if (!controller)
	return 0;

    for (unsigned int attempt = 0; attempt < max_attempts; attempt++) {
	BRLObolProgressiveStatus status;
	(void)controller->advanceProgressiveWork(NULL, &status);
	if (!status.hasMore)
	    return 1;
	if (sleep_milliseconds) {
	    std::this_thread::sleep_for(std::chrono::milliseconds(
		sleep_milliseconds));
	}
    }
    return 0;
}

static void
draw_test_obol_debug_dump_node(int id, int source_index, SoNode *node,
			       int depth, int child_index)
{
    if (!node || depth > 3)
	return;

    const char *type_name = node->getTypeId().getName().getString();
    bu_log("draw-obol-debug[%03d]: source[%d] node depth=%d child=%d ptr=%p type=%s isGroup=%d isMatrix=%d\n",
	   id, source_index, depth, child_index, (void *)node, type_name,
	   node->isOfType(SoGroup::getClassTypeId()) ? 1 : 0,
	   node->isOfType(SoMatrixTransform::getClassTypeId()) ? 1 : 0);

    if (node->isOfType(SoMatrixTransform::getClassTypeId())) {
	SoMatrixTransform *transform = static_cast<SoMatrixTransform *>(node);
	SbMatrix matrix = transform->matrix.getValue();
	const SbMat &m = matrix;
	bu_log("draw-obol-debug[%03d]: source[%d] node depth=%d matrix rows=[%.9g %.9g %.9g %.9g] [%.9g %.9g %.9g %.9g] [%.9g %.9g %.9g %.9g] [%.9g %.9g %.9g %.9g]\n",
	       id, source_index, depth,
	       m[0][0], m[0][1], m[0][2], m[0][3],
	       m[1][0], m[1][1], m[1][2], m[1][3],
	       m[2][0], m[2][1], m[2][2], m[2][3],
	       m[3][0], m[3][1], m[3][2], m[3][3]);
    }

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;

    SoGroup *group = static_cast<SoGroup *>(node);
    const int limit = group->getNumChildren() < 4 ?
		      group->getNumChildren() : 4;
    for (int i = 0; i < limit; i++)
	draw_test_obol_debug_dump_node(id, source_index, group->getChild(i),
				       depth + 1, i);
}

static int
draw_test_write_rgb_png(const char *filename, const unsigned char *buffer,
			int width, int height)
{
    if (!filename || !buffer || width <= 0 || height <= 0)
	return BRLCAD_ERROR;

    icv_image_t *img = icv_create(static_cast<size_t>(width),
					static_cast<size_t>(height), ICV_COLOR_SPACE_RGB);
    if (!img)
	return BRLCAD_ERROR;

    for (int y = 0; y < height; y++) {
	for (int x = 0; x < width; x++) {
	    const size_t src = (static_cast<size_t>(y) *
				static_cast<size_t>(width) + static_cast<size_t>(x)) * 3;
	    const size_t dst = (static_cast<size_t>(y) *
				static_cast<size_t>(width) + static_cast<size_t>(x)) * 3;
	    img->data[dst] = static_cast<double>(buffer[src]) / 255.0;
	    img->data[dst + 1] = static_cast<double>(buffer[src + 1]) / 255.0;
	    img->data[dst + 2] = static_cast<double>(buffer[src + 2]) / 255.0;
	}
    }

    const int ret = icv_write(img, filename, BU_MIME_IMAGE_PNG);
    icv_destroy(img);
    return ret;
}

static void
draw_test_obol_debug_dump(struct ged *gedp, int id,
			  BRLObolViewController *controller, void *view_ctx)
{
    if (!draw_test_obol_debug_enabled() || !gedp || !controller || !view_ctx)
	return;

    const char *debug_id = std::getenv("GED_TEST_DRAW_OBOL_DEBUG_ID");
    if (debug_id && debug_id[0] && std::atoi(debug_id) != id)
	return;

    point_t center = VINIT_ZERO;
    mat_t center_mat;
    if (bv_center_mat_get(center_mat, DRAW_TEST_BV_CONST(view_ctx)))
	MAT_DELTAS_GET_NEG(center, center_mat);

    point_t eye = VINIT_ZERO;
    (void)bv_eye_pos_get(eye, DRAW_TEST_BV_CONST(view_ctx));

    vect_t aet = VINIT_ZERO;
    (void)bv_aet_get(aet, DRAW_TEST_BV_CONST(view_ctx));

    bu_log("draw-obol-debug[%03d]: view=%s size=%g scale=%g dims=%dx%d center=(%.17g %.17g %.17g) eye=(%.17g %.17g %.17g) aet=(%.17g %.17g %.17g) perspective=%g\n",
	   id,
	   bv_name_get(DRAW_TEST_BV_CONST(view_ctx)) ?
	   bv_name_get(DRAW_TEST_BV_CONST(view_ctx)) : "(null)",
	   bv_size_get(DRAW_TEST_BV_CONST(view_ctx)),
	   bv_scale_get(DRAW_TEST_BV_CONST(view_ctx)),
	   bv_width_get(DRAW_TEST_BV_CONST(view_ctx)),
	   bv_height_get(DRAW_TEST_BV_CONST(view_ctx)),
	   center[X], center[Y], center[Z],
	   eye[X], eye[Y], eye[Z],
	   aet[X], aet[Y], aet[Z],
	   bv_perspective_get(DRAW_TEST_BV_CONST(view_ctx)));

    SoCamera *camera = controller->getCamera();
    if (camera) {
	const SbVec3f &pos = camera->position.getValue();
	bu_log("draw-obol-debug[%03d]: camera=%s pos=(%.9g %.9g %.9g) aspect=%.9g focal=%.9g near=%.9g far=%.9g\n",
	       id, camera->getTypeId().getName().getString(),
	       pos[0], pos[1], pos[2],
	       camera->aspectRatio.getValue(),
	       camera->focalDistance.getValue(),
	       camera->nearDistance.getValue(),
	       camera->farDistance.getValue());
	if (camera->isOfType(SoOrthographicCamera::getClassTypeId())) {
	    SoOrthographicCamera *ortho =
		static_cast<SoOrthographicCamera *>(camera);
	    bu_log("draw-obol-debug[%03d]: ortho-height=%.9g\n",
		   id, ortho->height.getValue());
	}
    }

    /* Endpoint-only controllers render the GED-owned shared root rather than
     * their empty private scene.  Report the scene that supplies autoview
     * bounds and database geometry. */
    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	scene = controller->getSceneController();
    if (!scene)
	return;

    BRLObolSceneSummary scene_summary;
    if (scene->getSceneSummary(scene_summary) && scene_summary.valid) {
	bu_log("draw-obol-debug[%03d]: scene rootChildren=%d dbSources=%d nonDbRootChildren=%d visited=%u realized=%u failed=%u structural=%" PRIu64 " frame=%" PRIu64 "\n",
	       id,
	       scene_summary.rootChildCount,
	       scene_summary.databaseSourceCount,
	       scene_summary.nonDatabaseRootChildCount,
	       scene_summary.lastVisitedSourceCount,
	       scene_summary.lastRealizedSourceCount,
	       scene_summary.lastFailedSourceCount,
	       scene_summary.structuralRevision,
	       scene_summary.frameRevision);
    }

    SbBox3f bounds;
    if (scene->getSceneSubtreeBounds("/", TRUE, bounds) && !bounds.isEmpty()) {
	const SbVec3f bmin = bounds.getMin();
	const SbVec3f bmax = bounds.getMax();
	bu_log("draw-obol-debug[%03d]: subtree-bounds / overlays=true min=(%.9g %.9g %.9g) max=(%.9g %.9g %.9g)\n",
	       id, bmin[0], bmin[1], bmin[2], bmax[0], bmax[1], bmax[2]);
    }
    if (scene->getSceneSubtreeBounds("/", FALSE, bounds) && !bounds.isEmpty()) {
	const SbVec3f bmin = bounds.getMin();
	const SbVec3f bmax = bounds.getMax();
	bu_log("draw-obol-debug[%03d]: subtree-bounds / overlays=false min=(%.9g %.9g %.9g) max=(%.9g %.9g %.9g)\n",
	       id, bmin[0], bmin[1], bmin[2], bmax[0], bmax[1], bmax[2]);
    }

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary source;
	if (!scene->getDatabaseSourceSummary(i, source) || !source.valid)
	    continue;
	SoBRLDatabaseSource *source_node = scene->getDatabaseSource(i);

	bu_log("draw-obol-debug[%03d]: source[%d] ptr=%p path=%s instance=%s repr=%s reprMode=%d mode=%d visible=%d selected=%d highlighted=%d colorOverride=%d color=(%.9g %.9g %.9g) materialValid=%d material=(%.9g %.9g %.9g) stale=%d status=%d shapes=%d meshes=%d drawMatrix=%d drawCenter=%d drawSize=%d sourceBounds=%d children=%d\n",
	       id, i,
	       (void *)source_node,
	       source.path.getString(),
	       source.instanceKey.getString(),
	       source.representationKey.getString(),
	       source.representationMode,
	       source.drawMode,
	       source.visible ? 1 : 0,
	       source.selected ? 1 : 0,
	       source.highlighted ? 1 : 0,
	       source.colorOverride ? 1 : 0,
	       source.color[0],
	       source.color[1],
	       source.color[2],
	       source.materialColorValid ? 1 : 0,
	       source.materialColor[0],
	       source.materialColor[1],
	       source.materialColor[2],
	       source.stale ? 1 : 0,
	       source.realizationStatus,
	       source.realizedShapeCount,
	       source.realizedMeshCount,
	       source.drawMatrixValid ? 1 : 0,
	       source.drawCenterValid ? 1 : 0,
	       source.drawSizeValid ? 1 : 0,
	       source.sourceBoundsValid ? 1 : 0,
	       source_node ? source_node->getNumChildren() : -1);
	if (source_node) {
	    for (int child_index = 0;
		 child_index < source_node->getNumChildren() &&
		 child_index < 4;
		 child_index++) {
		SoNode *child = source_node->getChild(child_index);
		const char *type_name = child ?
					child->getTypeId().getName().getString() : "<null>";
		bu_log("draw-obol-debug[%03d]: source[%d] child[%d]=%p type=%s isGroup=%d isMatrix=%d\n",
		       id, i, child_index, (void *)child, type_name,
		       child && child->isOfType(SoGroup::getClassTypeId()) ? 1 : 0,
		       child && child->isOfType(SoMatrixTransform::getClassTypeId()) ? 1 : 0);
		draw_test_obol_debug_dump_node(id, i, child, 1, child_index);
	    }
	}
	if (source.drawMatrixValid) {
	    const SbMat &m = source.drawMatrix.getValue();
	    bu_log("draw-obol-debug[%03d]: source[%d] drawMatrix rows=[%.9g %.9g %.9g %.9g] [%.9g %.9g %.9g %.9g] [%.9g %.9g %.9g %.9g] [%.9g %.9g %.9g %.9g]\n",
		   id, i,
		   m[0][0], m[0][1], m[0][2], m[0][3],
		   m[1][0], m[1][1], m[1][2], m[1][3],
		   m[2][0], m[2][1], m[2][2], m[2][3],
		   m[3][0], m[3][1], m[3][2], m[3][3]);
	}
	if (source.sourceBoundsValid && !source.sourceBounds.isEmpty()) {
	    const SbVec3f bmin = source.sourceBounds.getMin();
	    const SbVec3f bmax = source.sourceBounds.getMax();
	    bu_log("draw-obol-debug[%03d]: source[%d] sourceBounds min=(%.9g %.9g %.9g) max=(%.9g %.9g %.9g)\n",
		   id, i, bmin[0], bmin[1], bmin[2],
		   bmax[0], bmax[1], bmax[2]);
	}
	if (source_node) {
	    SoBRLVListShape *shape = source_node->getRealizedShape(0);
	    if (shape) {
		SbBox3f shape_bounds;
		shape_bounds.makeEmpty();
		for (int point_index = 0; point_index < shape->point.getNum();
		     point_index++)
		    shape_bounds.extendBy(shape->point[point_index]);
		if (!shape_bounds.isEmpty()) {
		    const SbVec3f bmin = shape_bounds.getMin();
		    const SbVec3f bmax = shape_bounds.getMax();
		    bu_log("draw-obol-debug[%03d]: source[%d] direct-shape ptr=%p selected=%d highlighted=%d points=%d bounds min=(%.9g %.9g %.9g) max=(%.9g %.9g %.9g)\n",
			   id, i, (void *)shape,
			   shape->selected.getValue() ? 1 : 0,
			   shape->highlighted.getValue() ? 1 : 0,
			   shape->point.getNum(),
			   bmin[0], bmin[1], bmin[2],
			   bmax[0], bmax[1], bmax[2]);
		}
	    }
	}
    }

    const int shape_count = scene->getRealizedShapeSummaryCount();
    const int shape_limit = shape_count < 32 ? shape_count : 32;
    bu_log("draw-obol-debug[%03d]: realized-shape-count=%d\n",
	   id, shape_count);
    for (int i = 0; i < shape_limit; i++) {
	BRLObolRealizedShapeSummary shape;
	if (!scene->getRealizedShapeSummary(i, shape) || !shape.valid)
	    continue;

	bu_log("draw-obol-debug[%03d]: shape[%d] path=%s owner=%s kind=%d drawMode=%d visible=%d selected=%d highlighted=%d points=%d commands=%d segments=%d tris=%d bounds=%d materialValid=%d material=(%.9g %.9g %.9g) region=%d air=%d gmater=%d los=%d shader=%s\n",
	       id, i,
	       shape.path.getString(),
	       shape.ownerSourcePath.getString(),
	       shape.shapeKind,
	       shape.drawMode,
	       shape.visible ? 1 : 0,
	       shape.selected ? 1 : 0,
	       shape.highlighted ? 1 : 0,
	       shape.pointCount,
	       shape.commandCount,
	       shape.segmentCount,
	       shape.triangleCount,
	       shape.boundsValid ? 1 : 0,
	       shape.materialColorValid ? 1 : 0,
	       shape.materialColor[0],
	       shape.materialColor[1],
	       shape.materialColor[2],
	       shape.regionId,
	       shape.airCode,
	       shape.materialId,
	       shape.los,
	       shape.materialShader.getString());
	if (shape.boundsValid && !shape.bounds.isEmpty()) {
	    const SbVec3f bmin = shape.bounds.getMin();
	    const SbVec3f bmax = shape.bounds.getMax();
	    bu_log("draw-obol-debug[%03d]: shape[%d] bounds min=(%.9g %.9g %.9g) max=(%.9g %.9g %.9g)\n",
		   id, i, bmin[0], bmin[1], bmin[2],
		   bmax[0], bmax[1], bmax[2]);
	}
    }
}

static int
draw_test_obol_screengrab_impl(struct ged *gedp, void *view_ctx, int id,
			       const char *filename)
{
    SoDB::ContextManager *manager = draw_test_obol_context_manager();
    if (!manager) {
	bu_log("Obol draw-test capture requested, but OSMesa context manager is unavailable\n");
	return -1;
    }
    brlobol_init(manager);

    void *v = view_ctx ? view_ctx : draw_test_active_view_ctx(gedp);
    if (!v)
	return -1;

    BRLObolViewController *controller =
	(BRLObolViewController *)ged_draw_obol_controller_opaque_for_view(v);
    if (!controller) {
	(void)ged_draw_obol_scene_controller_ensure(gedp, 1);
	controller = ged_draw_obol_controller(gedp);
    }
    if (!controller)
	return -1;

    int width = bv_width_get(DRAW_TEST_BV_CONST(v));
    int height = bv_height_get(DRAW_TEST_BV_CONST(v));
    if (width <= 0 || height <= 0)
	return -1;

    controller->setViewportSize(static_cast<unsigned int>(width),
				static_cast<unsigned int>(height));
    if (!draw_test_sync_obol_camera(controller, v))
	return -1;
    (void)ged_draw_obol_view_context_faceplate_sync(gedp, v);
    (void)controller->realizePending();
    (void)ged_selection_draw_sync(gedp, NULL);
    draw_test_obol_debug_dump(gedp, id, controller, v);

    const SbViewportRegion &region = controller->getViewportRegion();
    SbVec2s size = region.getViewportSizePixels();
    if (size[0] <= 0 || size[1] <= 0 || !controller->getViewport())
	return -1;

    SoOffscreenRenderer renderer(manager, region);
    renderer.setComponents(SoOffscreenRenderer::RGB);
    renderer.setBackgroundColor(SbColor(0.0f, 0.0f, 0.0f));
    /* Capture the composed presentation graph, including framebuffer layers.
     * Rendering only the viewport scene drops retained underlay/overlay nodes. */
    if (!renderer.render(controller->getRenderRoot()))
	return -1;

    const unsigned char *buffer = renderer.getBuffer();
    if (!buffer)
	return -1;

    return draw_test_write_rgb_png(filename, buffer, size[0], size[1]) == BRLCAD_OK ?
	   1 : -1;
}

extern "C" int
draw_test_obol_screengrab(struct ged *gedp, int id, const char *filename)
{
    return draw_test_obol_screengrab_impl(gedp, NULL, id, filename);
}

static int
draw_test_obol_screengrab_if_enabled(struct ged *gedp, int id,
				     const char *filename)
{
    if (!draw_test_obol_capture_enabled())
	return 0;

    return draw_test_obol_screengrab_impl(gedp, NULL, id, filename);
}

extern "C" int
draw_test_obol_screengrab_view_if_enabled(struct ged *gedp, void *view_ctx,
	int id, const char *filename)
{
    if (!draw_test_obol_capture_enabled())
	return 0;

    return draw_test_obol_screengrab_impl(gedp, view_ctx, id, filename);
}

extern "C" long
draw_test_count_nonblack_pixels(const char *filename)
{
    icv_image_t *img = icv_read(filename, BU_MIME_IMAGE_PNG, 0, 0);
    if (!img)
	return -1;

    long cnt = 0;
    size_t npix = (size_t)img->width * img->height;
    double *d = img->data;
    for (size_t i = 0; i < npix; i++) {
	if (d[3*i] > 0.001 || d[3*i+1] > 0.001 || d[3*i+2] > 0.001)
	    cnt++;
    }

    icv_destroy(img);
    return cnt;
}

extern "C" int
draw_test_images_differ(const char *a, const char *b, int offmany_threshold)
{
    icv_image_t *ia = icv_read(a, BU_MIME_IMAGE_PNG, 0, 0);
    icv_image_t *ib = icv_read(b, BU_MIME_IMAGE_PNG, 0, 0);
    if (!ia || !ib) {
	if (ia)
	    icv_destroy(ia);
	if (ib)
	    icv_destroy(ib);
	return -1;
    }

    int match = 0;
    int off1 = 0;
    int offmany = 0;
    int diff = icv_diff(&match, &off1, &offmany, ia, ib);
    icv_destroy(ia);
    icv_destroy(ib);

    if (!diff)
	return 0;

    return (offmany >= offmany_threshold) ? 1 : 0;
}

extern "C" void
dm_refresh(struct ged *gedp)
{
    void *v = draw_test_active_view_ctx(gedp);
    if (!v)
	return;
    struct ged_draw_transaction txn =
	ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
    txn.view = v;
    ged_draw_apply_transaction(gedp, &txn, NULL);

    BRLObolViewController *controller =
	(BRLObolViewController *)ged_draw_obol_controller_opaque_for_view(v);
    if (controller)
	(void)controller->realizePending();
}

extern "C" void
scene_clear(struct ged *gedp)
{
    const char *s_av[1] = {"Z"};
    ged_exec_Z(gedp, 1, s_av);
    dm_refresh(gedp);
}

extern "C" int
unpack_apng(const char *src_dir, const char *apng_name, const char *out_dir, const char *prefix)
{
    // copy empty.png
    struct bu_vls empty_src = BU_VLS_INIT_ZERO;
    struct bu_vls empty_dst = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&empty_src, "%s/empty.png", src_dir);
    bu_vls_sprintf(&empty_dst, "%s/empty.png", out_dir);
    {
        std::ifstream empty_in(bu_vls_cstr(&empty_src), std::ios::binary);
        std::ofstream empty_out(bu_vls_cstr(&empty_dst), std::ios::binary);
        if (empty_in && empty_out) {
            empty_out << empty_in.rdbuf();
        }
    }
    bu_vls_free(&empty_src);
    bu_vls_free(&empty_dst);

    // read apng
    struct bu_vls apng_path = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&apng_path, "%s/%s", src_dir, apng_name);

    icv_anim_t *anim = icv_anim_read(bu_vls_cstr(&apng_path));
    if (!anim) {
        bu_log("Failed to read APNG: %s\n", bu_vls_cstr(&apng_path));
        bu_vls_free(&apng_path);
        return -1;
    }
    bu_vls_free(&apng_path);

    size_t num = icv_anim_num_frames(anim);
    for (size_t i = 0; i < num; ++i) {
        icv_image_t *img = icv_anim_get_frame(anim, i);
        if (img) {
            struct bu_vls out_path = BU_VLS_INIT_ZERO;
            bu_vls_sprintf(&out_path, "%s/%s%03zu_ctrl.png", out_dir, prefix, i + 1);
            if (icv_write(img, bu_vls_cstr(&out_path), BU_MIME_IMAGE_PNG) != 0) {
                bu_log("Failed to write extracted frame: %s\n", bu_vls_cstr(&out_path));
                icv_destroy(img);
                bu_vls_free(&out_path);
                icv_anim_destroy(anim);
                return -1;
            }
            icv_destroy(img);
            bu_vls_free(&out_path);
        }
    }

    icv_anim_destroy(anim);
    return 0;
}


extern "C" int
img_cmp(int id, struct ged *gedp, const char *cdir, bool clear_scene, bool clear_image, int soft_fail, fastf_t approximate_check, const char *clear_root, const char *img_root)
{
    icv_image_t *ctrl, *timg;
    struct bu_vls tname = BU_VLS_INIT_ZERO;
    struct bu_vls cname = BU_VLS_INIT_ZERO;
    if (id <= 0) {
	bu_vls_sprintf(&tname, "%s.png", clear_root);
	bu_vls_sprintf(&cname, "%s/empty.png", cdir);
    } else {
	bu_vls_sprintf(&tname, "%s%03d.png", img_root, id);
	bu_vls_sprintf(&cname, "%s/%s%03d_ctrl.png", cdir, img_root, id);
    }

    dm_refresh(gedp);
    int obol_capture = draw_test_obol_screengrab_if_enabled(gedp, id,
		       bu_vls_cstr(&tname));
    if (obol_capture != 1)
	bu_file_delete(bu_vls_cstr(&tname));

    timg = icv_read(bu_vls_cstr(&tname), BU_MIME_IMAGE_PNG, 0, 0);
    if (!timg) {
	if (soft_fail) {
	    bu_log("Failed to read %s\n", bu_vls_cstr(&tname));
	    if (clear_scene)
		scene_clear(gedp);
	    bu_vls_free(&tname);
	    return BRLCAD_ERROR;
	}
	bu_exit(EXIT_FAILURE, "failed to read %s\n", bu_vls_cstr(&tname));
    }
    ctrl = icv_read(bu_vls_cstr(&cname), BU_MIME_IMAGE_PNG, 0, 0);
    if (!ctrl) {
	if (soft_fail) {
	    bu_log("Failed to read %s\n", bu_vls_cstr(&cname));
	    if (clear_scene)
		scene_clear(gedp);
	    bu_vls_free(&tname);
	    bu_vls_free(&cname);
	    return BRLCAD_ERROR;
	}
	bu_exit(EXIT_FAILURE, "failed to read %s\n", bu_vls_cstr(&cname));
    }
    if (draw_test_compare_debug_enabled()) {
	bu_log("draw-compare-debug[%d]: target=%s ctrl=%s target-size=%zux%zu ctrl-size=%zux%zu target-cspace=%d ctrl-cspace=%d\n",
	       id, bu_vls_cstr(&tname), bu_vls_cstr(&cname),
	       timg->width, timg->height, ctrl->width, ctrl->height,
	       timg->color_space, ctrl->color_space);
    }
    bu_vls_free(&cname);
    int matching_cnt = 0;
    int off_by_1_cnt = 0;
    int off_by_many_cnt = 0;
    int iret = icv_diff(&matching_cnt, &off_by_1_cnt, &off_by_many_cnt, ctrl, timg);
    if (iret) {
	if (approximate_check > 0) {
	    // First, if we're allowing approximate and all we have are off by one errors,
	    // allow it.
	    if (!off_by_many_cnt) {
		bu_log("%d approximate matching enabled, no off by many - passing.  %d matching, %d off by 1\n", id, matching_cnt, off_by_1_cnt);
		iret = 0;
		// We don't need it as a pass/fail, but report the SSIM score.
		fastf_t pret = icv_adiff(ctrl, timg, ICV_DIFF_SSIM);
		bu_log("icv_adiff SSIM metric(%d): %g \n", id, pret);
	    } else {
		// We have off-by-many pixels, so compare structural similarity.
		fastf_t pret = icv_adiff(ctrl, timg, ICV_DIFF_SSIM);
		// SSIM ranges from -1 to 1, with 1 indicating identical images.
		// approximate_check specifies the minimum accepted similarity.
		bu_log("icv_adiff SSIM metric(%d): %g \n", id, pret);
		if (pret > approximate_check) {
		    iret = 0;
		}
	    }
	}

	if (iret) {

	    bu_log("%d %s diff failed.  %d matching, %d off by 1, %d off by many\n", id, img_root, matching_cnt, off_by_1_cnt, off_by_many_cnt);

	    // Generate a diff image for debugging work
	    icv_image_t *diff_img = icv_diffimg(ctrl, timg);
	    if (diff_img) {
		struct bu_vls diff_name = BU_VLS_INIT_ZERO;
		bu_vls_sprintf(&diff_name, "%s_%d_diff.png", img_root, id);
		icv_write(diff_img, bu_vls_cstr(&diff_name), BU_MIME_IMAGE_PNG);
		bu_vls_free(&diff_name);
	    }

	    // If we're in soft fail mode, we're not keeping the image unless
	    // the user requested it. In hard fail, we leave the last image for
	    // inspection.
	    if (soft_fail && clear_image)
		bu_file_delete(bu_vls_cstr(&tname));
#if 0
	    // Dump an ascii rendering of the difference image to the log.  In
	    // scenarios such as CI systems, where we don't have a way to
	    // readily inspect the difference images, this output can still
	    // give us at least some idea of what is going on as long as the
	    // user can get their terminal font small enough.  We deliberately
	    // don't do any resizing of the image - we want to leave diff
	    // pixels intact, so we accept the terminal aspect ratio distortion
	    // as a trade-off for pixel rendering accuracy.
	    if (diff_img) {
		struct icv_ascii_art_params iparams = ICV_ASCII_ART_PARAMS_DEFAULT;
		iparams.output_color = 1;
		char *aart = icv_ascii_art(diff_img, &iparams);
		bu_log("%s\n", aart);
	    }
#endif
	    // Done with images
	    bu_vls_free(&tname);
	    icv_destroy(ctrl);
	    icv_destroy(timg);
	    icv_destroy(diff_img);

	    // If we're in soft fail, return
	    if (soft_fail) {
		if (clear_scene)
		    scene_clear(gedp);
		return BRLCAD_ERROR;
	    }

	    // Hard fail
	    bu_exit(EXIT_FAILURE, "Diff failure found, test aborted.\n");

	}
    }

    // Clean up
    if (clear_image)
	bu_file_delete(bu_vls_cstr(&tname));
    bu_vls_free(&tname);
    icv_destroy(ctrl);
    icv_destroy(timg);

    // If we're supposed to clear the scene, do so now.
    if (clear_scene)
	scene_clear(gedp);

    return BRLCAD_OK;
}

extern "C" int
img_not_empty(int id, struct ged *gedp, const char *cdir, bool clear_scene, bool clear_image, int soft_fail, const char *UNUSED(clear_root), const char *img_root)
{
    icv_image_t *ctrl, *timg;
    struct bu_vls tname = BU_VLS_INIT_ZERO;
    struct bu_vls cname = BU_VLS_INIT_ZERO;

    bu_vls_sprintf(&tname, "%s%03d_semantic.png", img_root, id);
    bu_vls_sprintf(&cname, "%s/empty.png", cdir);

    dm_refresh(gedp);
    int obol_capture = draw_test_obol_screengrab_if_enabled(gedp, id,
		       bu_vls_cstr(&tname));
    if (obol_capture != 1)
	bu_file_delete(bu_vls_cstr(&tname));

    timg = icv_read(bu_vls_cstr(&tname), BU_MIME_IMAGE_PNG, 0, 0);
    ctrl = icv_read(bu_vls_cstr(&cname), BU_MIME_IMAGE_PNG, 0, 0);
    bu_vls_free(&cname);

    if (!timg || !ctrl) {
	if (timg) icv_destroy(timg);
	if (ctrl) icv_destroy(ctrl);
	if (clear_image)
	    bu_file_delete(bu_vls_cstr(&tname));
	bu_vls_free(&tname);
	if (clear_scene)
	    scene_clear(gedp);
	if (soft_fail)
	    return BRLCAD_ERROR;
	bu_exit(EXIT_FAILURE, "failed to read semantic image comparison inputs\n");
    }

    int matching_cnt = 0;
    int off_by_1_cnt = 0;
    int off_by_many_cnt = 0;
    int diff = icv_diff(&matching_cnt, &off_by_1_cnt, &off_by_many_cnt, ctrl, timg);
    icv_destroy(ctrl);
    icv_destroy(timg);

    if (clear_image)
	bu_file_delete(bu_vls_cstr(&tname));
    bu_vls_free(&tname);

    if (!diff) {
	bu_log("%d %s semantic non-empty check failed: image matches empty scene\n", id, img_root);
	if (clear_scene)
	    scene_clear(gedp);
	if (soft_fail)
	    return BRLCAD_ERROR;
	bu_exit(EXIT_FAILURE, "Semantic empty-image failure found, test aborted.\n");
    }

    bu_log("%d %s semantic non-empty check passed against empty scene.  %d empty pixels matching, %d nearly-empty pixels, %d visibly drawn pixels\n",
	   id, img_root, matching_cnt, off_by_1_cnt, off_by_many_cnt);

    if (clear_scene)
	scene_clear(gedp);

    return BRLCAD_OK;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
