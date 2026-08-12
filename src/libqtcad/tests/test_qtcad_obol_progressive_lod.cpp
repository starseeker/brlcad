/*     T E S T _ Q T C A D _ O B O L _ P R O G R E S S I V E _ L O D . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include "bv.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BDrawCache.h"
#include "BObol/BLodService.h"
#include "BObol/BMeshShape.h"
#include "BObol/BSceneController.h"
#include "BObol/BVListShape.h"
#include "BObol/BViewLod.h"
#include "BObol/BViewController.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/path.h"
#include "bu/process.h"
#include "bu/str.h"
#include "bu/datetime.h"
#include "ged.h"
#include "ged/db_index.h"
#include "ged/draw.h"
#include "ged/display.h"
#include "ged/event.h"
#include "qtcad/QgView.h"

#include <Inventor/SoViewport.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>

#include <QApplication>
#include <QCoreApplication>
#include <QImage>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctype.h>
#include <fstream>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

struct progressive_lod_case {
    const char *name;
    const char *sourceDb;
    const char *facetizeRoot;
    const char *drawTarget;
    int useFacetize;
    unsigned int minVisitedMeshCount;
    size_t minActivePayloadCount;
    int aabbTimeoutMs;
    int obbTimeoutMs;
    int meshTimeoutMs;
    int copyCount;
    int startupOnly;
    int startupRefine;
    int startupExpand;
    int startupAutoExpand;
    int startupWireLod;
    int shadedLod;
    int slowLodDelays;
    int diagnosticStages;
    int startupExpandPasses;
};

struct progressive_lod_timings {
    double copySeconds;
    double openSeconds;
    double viewSetupSeconds;
    double cacheClearSeconds;
    double tolSeconds;
    double facetizeSeconds;
    double repeatTargetSeconds;
    double lodPolicySeconds;
    double drawApplySeconds;
    double syncSeconds;
    double viewAllSeconds;
    double firstReadbackSeconds;
    double firstVisualSeconds;
    double serviceStartSeconds;
    double submitSeconds;
    double aabbSeconds;
    double obbSeconds;
    double meshSeconds;
};

static void
sanitize_case_name(char *dst, size_t dst_len, const char *src)
{
    size_t out = 0;

    if (!dst || dst_len == 0)
	return;

    if (!src || !src[0])
	src = "case";

    for (const char *p = src; *p && out + 1 < dst_len; p++) {
	unsigned char ch = (unsigned char)*p;
	dst[out++] = isalnum(ch) ? (char)ch : '_';
    }

    if (out == 0 && dst_len > 1)
	dst[out++] = 'x';

    dst[out] = '\0';
}

static int
parse_positive_int(const char *str, int fallback)
{
    if (!str || !str[0])
	return fallback;

    char *end = NULL;
    long value = strtol(str, &end, 10);
    if (end == str || value <= 0 || value > INT_MAX)
	return fallback;

    return (int)value;
}

static size_t
parse_positive_size(const char *str, size_t fallback)
{
    if (!str || !str[0])
	return fallback;

    char *end = NULL;
    unsigned long value = strtoul(str, &end, 10);
    if (end == str || value == 0)
	return fallback;

    return (size_t)value;
}

static double
elapsed_seconds(int64_t start)
{
    return (double)(bu_gettime() - start) / 1000000.0;
}

static int
phase_logging_enabled(const struct progressive_lod_case &testCase)
{
    const char *env = getenv("QTCAD_OBOL_PROGRESSIVE_LOD_PHASE_LOG");
    return testCase.startupOnly || testCase.startupRefine ||
	   testCase.startupWireLod || testCase.shadedLod ||
	   testCase.diagnosticStages ||
	   (env && BU_STR_EQUAL(env, "1"));
}

static void
record_progressive_lod_phase(double &slot,
			     int64_t phaseStart,
			     int64_t totalStart,
			     const struct progressive_lod_case &testCase,
			     const char *phase)
{
    slot = elapsed_seconds(phaseStart);
    if (phase_logging_enabled(testCase)) {
	fprintf(stderr,
		"qtcad_progressive_lod_phase case=%s phase=%s seconds=%.3f total=%.3f\n",
		testCase.name, phase, slot, elapsed_seconds(totalStart));
	fflush(stderr);
    }
}

static void
print_progressive_lod_timings(FILE *fp,
			      const struct progressive_lod_timings &timings)
{
    fprintf(fp, " timing_copy=%.3f timing_open=%.3f timing_view_setup=%.3f",
	    timings.copySeconds, timings.openSeconds,
	    timings.viewSetupSeconds);
    fprintf(fp, " timing_cache_clear=%.3f timing_tol=%.3f",
	    timings.cacheClearSeconds, timings.tolSeconds);
    fprintf(fp, " timing_facetize=%.3f timing_repeat=%.3f",
	    timings.facetizeSeconds, timings.repeatTargetSeconds);
    fprintf(fp, " timing_lod_policy=%.3f timing_draw_apply=%.3f",
	    timings.lodPolicySeconds, timings.drawApplySeconds);
    fprintf(fp, " timing_sync=%.3f timing_view_all=%.3f",
	    timings.syncSeconds, timings.viewAllSeconds);
    fprintf(fp, " timing_readback0=%.3f timing_first_visual=%.3f",
	    timings.firstReadbackSeconds, timings.firstVisualSeconds);
    fprintf(fp, " timing_service_start=%.3f timing_submit=%.3f",
	    timings.serviceStartSeconds, timings.submitSeconds);
    fprintf(fp, " timing_aabb=%.3f timing_obb=%.3f timing_mesh=%.3f",
	    timings.aabbSeconds, timings.obbSeconds, timings.meshSeconds);
}


static int
qtcad_obol_render_database_source_count(struct ged *gedp,
					BObolViewController *controller)
{
    (void)gedp;
    if (!controller)
	return 0;
    BObolSceneController renderScene(controller->getRenderSceneRoot());
    return renderScene.getDatabaseSourceCount();
}

static int
qtcad_obol_path_depth(const char *path)
{
    int depth = 0;
    int in_component = 0;

    if (!path)
	return 0;

    while (*path) {
	if (*path == '/') {
	    in_component = 0;
	} else if (!in_component) {
	    depth++;
	    in_component = 1;
	}
	path++;
    }

    return depth;
}

static int
qtcad_obol_scene_database_source_max_depth(
    BObolSceneController *scene,
    const char *root_path)
{
    int max_depth = 0;
    if (!scene || !root_path || !root_path[0])
	return 0;

    const char *root = root_path;
    while (*root == '/')
	root++;
    const size_t root_len = strlen(root);

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;
	const char *path = summary.path.getString();
	while (path && *path == '/')
	    path++;
	if (!path || bu_strncmp(path, root, root_len) != 0 ||
	    (path[root_len] != '\0' && path[root_len] != '/'))
	    continue;
	int depth = qtcad_obol_path_depth(path);
	if (depth > max_depth)
	    max_depth = depth;
    }

    return max_depth;
}

static int
qtcad_obol_render_database_source_max_depth(
    struct ged *gedp,
    BObolViewController *controller,
    const char *root_path)
{
    (void)gedp;
    if (!controller)
	return 0;
    BObolSceneController renderScene(controller->getRenderSceneRoot());
    return qtcad_obol_scene_database_source_max_depth(&renderScene,
	root_path);
}


static int
qtcad_obol_scene_realized_database_geometry_count(BObolSceneController *scene)
{
    int count = 0;
    if (!scene)
	return 0;

    const int sourceCount = scene->getDatabaseSourceCount();
	for (int i = 0; i < sourceCount; i++) {
	    SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	    if (!source)
		continue;
	    if (source->hasRealizedWireGeometry())
		count++;
	    if (source->hasRealizedMeshGeometry())
		count++;
    }

    return count;
}


static int
qtcad_obol_render_realized_database_geometry_count(
    struct ged *gedp,
    BObolViewController *controller)
{
    (void)gedp;
    if (!controller)
	return 0;
    BObolSceneController renderScene(controller->getRenderSceneRoot());
    return qtcad_obol_scene_realized_database_geometry_count(&renderScene);
}

static int
qtcad_obol_render_compact_instance_count(
    struct ged *gedp,
    BObolViewController *controller)
{
    (void)gedp;
    if (!controller)
	return 0;
    BObolSceneController renderScene(controller->getRenderSceneRoot());
    int count = 0;
    for (int i = 0; i < renderScene.getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = renderScene.getDatabaseSource(i);
	if (source && source->isCompactOccurrenceRegistry())
	    count += source->getCompactInstanceCount();
    }
    return count;
}

static size_t
qtcad_obol_render_visible_compact_mesh_demand_count(
    BObolViewController *controller, int includeSubtract,
    int uniqueOccurrenceKeys = 0, int inFrustumOnly = 0)
{
    if (!controller)
	return 0;
    BObolSceneController renderScene(controller->getRenderSceneRoot());
    size_t count = 0;
    std::unordered_set<std::string> occurrenceKeys;
    SbViewVolume viewVolume;
    if (inFrustumOnly) {
	const SbVec2s viewport =
	    controller->getViewportRegion().getViewportSizePixels();
	const float aspect = viewport[1] > 0 ?
	    static_cast<float>(viewport[0]) /
		static_cast<float>(viewport[1]) : 1.0f;
	viewVolume = controller->getCamera()->getViewVolume(aspect);
    }
    for (int sourceIndex = 0;
	 sourceIndex < renderScene.getDatabaseSourceCount(); sourceIndex++) {
	SoBRLDatabaseSource *source =
	    renderScene.getDatabaseSource(sourceIndex);
	if (!source || !source->isCompactOccurrenceRegistry())
	    continue;
	for (int instanceIndex = 0;
	     instanceIndex < source->getCompactInstanceCount();
	     instanceIndex++) {
	    BObolCompactInstanceHandle handle;
	    BObolCompactInstanceSummary summary;
	    if (!source->getCompactInstanceHandle(instanceIndex, handle) ||
		!source->getCompactInstanceSummary(handle, summary) ||
		!summary.valid || !summary.visible)
		continue;
	    if (!includeSubtract &&
		summary.booleanOperation ==
		    SoBRLDatabaseSource::BOOLEAN_SUBTRACT)
		continue;
	    if (inFrustumOnly) {
		SbBox3f worldBounds = summary.localBounds;
		worldBounds.transform(summary.localToSource);
		if (worldBounds.isEmpty() ||
		    !viewVolume.intersect(worldBounds))
		    continue;
	    }
	    /* A resident base mesh is already richer than a box and may be below
	     * the configured PoP threshold.  Count only occurrences that actually
	     * advertise source-backed/view-LoD demand. */
	    if (summary.sourceMeshRequestValid || summary.lodBacked) {
		if (uniqueOccurrenceKeys)
		    occurrenceKeys.insert(
			summary.sourceInstanceKey.getString());
		else
		    count++;
	    }
	}
    }
    return uniqueOccurrenceKeys ? occurrenceKeys.size() : count;
}


static int
qtcad_obol_color_matches_metadata(const SbColor &color,
				  const unsigned char expected[3])
{
    return std::fabs(static_cast<double>(color[0]) -
		     static_cast<double>(expected[0]) / 255.0) < 1.0e-5 &&
	   std::fabs(static_cast<double>(color[1]) -
		     static_cast<double>(expected[1]) / 255.0) < 1.0e-5 &&
	   std::fabs(static_cast<double>(color[2]) -
		     static_cast<double>(expected[2]) / 255.0) < 1.0e-5;
}


template <typename ShapeT>
static int
qtcad_obol_shape_matches_draw_metadata(
    ShapeT *shape,
    const struct BObolDrawMetadataRecord *record)
{
    if (!shape || !record)
	return 0;

    if (record->hasRegionId &&
	shape->regionId.getValue() != record->regionId)
	return 0;
    if (record->hasAircode &&
	shape->airCode.getValue() != record->aircode)
	return 0;
    if (record->hasLos && shape->los.getValue() != record->los)
	return 0;
    if (record->hasMaterialId &&
	shape->materialId.getValue() != record->materialId)
	return 0;
    if (record->hasColor &&
	(!shape->materialColorValid.getValue() ||
	 !qtcad_obol_color_matches_metadata(
	     shape->materialColor.getValue(), record->color)))
	return 0;
    if (record->hasShader &&
	!BU_STR_EQUAL(shape->materialShader.getValue().getString(),
		      record->shader))
	return 0;

    return 1;
}


static int
qtcad_obol_source_matches_draw_metadata(
    SoBRLDatabaseSource *source,
    const struct BObolDrawMetadataRecord *record)
{
    if (!source || !record)
	return 0;

    BObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid ||
	!summary.databaseMetadataValid)
	return 0;

    if (record->hasRegionId &&
	summary.databaseRegionId != record->regionId)
	return 0;
    if (record->hasAircode &&
	summary.databaseAirCode != record->aircode)
	return 0;
    if (record->hasLos && summary.databaseLos != record->los)
	return 0;
    if (record->hasMaterialId &&
	summary.databaseMaterialId != record->materialId)
	return 0;
    if (record->hasColor &&
	(!summary.databaseMaterialColorValid ||
	 !qtcad_obol_color_matches_metadata(
	     summary.databaseMaterialColor, record->color)))
	return 0;
    if (record->hasShader &&
	!BU_STR_EQUAL(summary.databaseMaterialShader.getString(),
		      record->shader))
	return 0;

    return 1;
}


static int
qtcad_obol_scene_has_draw_metadata(
    BObolSceneController *scene,
    const struct BObolDrawMetadataRecord *record)
{
    if (!scene || !record)
	return 0;

    const int sourceCount = scene->getDatabaseSourceCount();
    for (int i = 0; i < sourceCount; i++) {
	SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	if (!source)
	    continue;
	if (qtcad_obol_source_matches_draw_metadata(source, record))
	    return 1;
	const int vlistCount = source->getRealizedShapeCount();
	for (int j = 0; j < vlistCount; j++) {
	    if (qtcad_obol_shape_matches_draw_metadata(
		    source->getRealizedShape(j), record))
		return 1;
	}
	const int meshCount = source->getRealizedMeshCount();
	for (int j = 0; j < meshCount; j++) {
	    if (qtcad_obol_shape_matches_draw_metadata(
		    source->getRealizedMesh(j), record))
		return 1;
	}
    }

    return 0;
}


static int
qtcad_obol_render_has_draw_metadata(
    struct ged *gedp,
    BObolViewController *controller,
    const struct BObolDrawMetadataRecord *record)
{
    (void)gedp;
    if (!controller)
	return 0;
    BObolSceneController renderScene(controller->getRenderSceneRoot());
    return qtcad_obol_scene_has_draw_metadata(&renderScene, record);
}


static int
qtcad_obol_autoview_refresh(struct ged *gedp,
			    QgView &view,
			    BObolViewController *controller,
			    const char *render_reason)
{
    const char *autoview_cmd[1] = {"autoview"};
    if (ged_exec_autoview(gedp, 1, autoview_cmd) != BRLCAD_OK)
	return 0;

    view.need_update(QG_VIEW_REFRESH);
    if (controller)
	controller->requestRender(render_reason ? render_reason :
				  "progressive-lod-autoview");
    QCoreApplication::processEvents();
    return 1;
}

static void
qtcad_obol_request_view_frame(QgView &view,
	BObolViewController *controller, const char *render_reason)
{
    view.need_update(QG_VIEW_REFRESH);
    if (controller)
	controller->requestRender(render_reason ? render_reason :
		"progressive-lod-frame");
    QCoreApplication::processEvents();
}

/* A hidden software canvas has no Qt paint delivery.  Keep ordinary image
 * readback observational, but let tests explicitly emulate the presentation
 * contract when a progressive admission gate is waiting for a rendered-frame
 * timing sample.  This mirrors QgSW::paintEvent: advance/publish, traverse,
 * report timing, and consume only the request that was actually presented. */
static void
qtcad_obol_present_view_frame(QgView &view,
	BObolViewController *controller, QImage &image)
{
    image = QImage();
    if (!controller)
	return;

    (void)controller->advanceProgressiveWork(NULL, NULL);
    const uint64_t started = controller->beginRenderTiming();
    view.get_viewport_image(image);
    controller->completeRenderTiming(started);
    if (!image.isNull() && controller->isRenderRequested())
	(void)controller->consumeRenderRequest(NULL);
}

static int
qtcad_obol_scene_autoview_refresh(QgView &view,
	BObolViewController *controller, const char *render_reason)
{
    struct ged_view_context *view_ctx =
	ged_view_context_from_bv(view.viewContext());
    SoNode *root = controller ? controller->getRenderSceneRoot() : NULL;
    if (!root && controller)
	root = controller->getSceneRoot();
    if (!view_ctx || !root)
	return 0;

    SoGetBoundingBoxAction bboxAction(controller->getViewportRegion());
    bboxAction.apply(root);
    const SbBox3f bounds = bboxAction.getBoundingBox();
    if (bounds.isEmpty())
	return 0;

    point_t bmin;
    point_t bmax;
    const SbVec3f min = bounds.getMin();
    const SbVec3f max = bounds.getMax();
    VSET(bmin, min[0], min[1], min[2]);
    VSET(bmax, max[0], max[1], max[2]);
	if (!bv_autoview_bounds(bv_context_view(ged_view_context_bv(view_ctx)),
	BV_AUTOVIEW_SCALE_DEFAULT, bmin, bmax))
	return 0;
	bv_context_update(ged_view_context_bv(view_ctx),
	BV_CONTEXT_CHANGED_VIEW);
    qtcad_obol_request_view_frame(view, controller, render_reason);
    return controller->syncCameraFromViewContext(
	view_ctx) ? 1 : 0;
}


static void
qtcad_obol_print_scene_source_diagnostics(const char *label,
	BObolSceneController *scene)
{
    if (!label || !scene)
	return;

    const int sourceCount = scene->getDatabaseSourceCount();
    fprintf(stderr,
	    " qtcad_progressive_lod_diag_source_set label=%s count=%d",
	    label, sourceCount);
    for (int i = 0; i < sourceCount; i++) {
	SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	if (!source)
	    continue;
	SbBox3f sourceBounds;
	const SbBool sourceBoundsValid = source->getSourceBounds(sourceBounds);
	fprintf(stderr,
		" source[%d]={ptr=%p path=%s instance=%s representation=%d draw_mode=%d role_flags=%d view_dependent=%d mesh_lod=%d bot_threshold=%u children=%d realized_shapes=%d realized_meshes=%d bounds_valid=%d",
		i, (void *)source, source->path.getValue().getString(),
		source->instanceKey.getValue().getString(),
		source->representationMode.getValue(),
		source->drawMode.getValue(),
		source->realizationRoleFlags.getValue(),
		source->realizationViewDependent.getValue() ? 1 : 0,
		source->realizationMeshLodEnabled.getValue() ? 1 : 0,
		source->lodBotThreshold.getValue(),
		source->getNumChildren(), source->getRealizedShapeCount(),
		source->getRealizedMeshCount(), sourceBoundsValid ? 1 : 0);
	if (source->getRealizedShapeCount() > 0) {
	    SoBRLVListShape *shape = source->getRealizedShape(0);
	    if (shape) {
		SbColor color = shape->materialColor.getValue();
		fprintf(stderr,
			" first_vlist_region=%d first_vlist_air=%d first_vlist_los=%d first_vlist_material=%d first_vlist_color_valid=%d first_vlist_color=(%.9g %.9g %.9g) first_vlist_shader=%s",
			shape->regionId.getValue(), shape->airCode.getValue(),
			shape->los.getValue(), shape->materialId.getValue(),
			shape->materialColorValid.getValue() ? 1 : 0,
			color[0], color[1], color[2],
			shape->materialShader.getValue().getString());
	    }
	}
	if (source->getRealizedMeshCount() > 0) {
	    SoBRLMeshShape *mesh = source->getRealizedMesh(0);
	    if (mesh) {
		SbColor color = mesh->materialColor.getValue();
		fprintf(stderr,
			" first_mesh_region=%d first_mesh_air=%d first_mesh_los=%d first_mesh_material=%d first_mesh_color_valid=%d first_mesh_color=(%.9g %.9g %.9g) first_mesh_shader=%s",
			mesh->regionId.getValue(), mesh->airCode.getValue(),
			mesh->los.getValue(), mesh->materialId.getValue(),
			mesh->materialColorValid.getValue() ? 1 : 0,
			color[0], color[1], color[2],
			mesh->materialShader.getValue().getString());
	    }
	}
	if (sourceBoundsValid) {
	    fprintf(stderr,
		    " bounds_min=(%.9g %.9g %.9g) bounds_max=(%.9g %.9g %.9g)",
		    sourceBounds.getMin()[0], sourceBounds.getMin()[1],
		    sourceBounds.getMin()[2], sourceBounds.getMax()[0],
		    sourceBounds.getMax()[1], sourceBounds.getMax()[2]);
	}
	fprintf(stderr, "}");
    }
    fprintf(stderr, "\n");
}


static void
qtcad_obol_print_render_diagnostics(struct ged *gedp,
				    const struct progressive_lod_case &testCase,
				    BObolViewController *controller,
				    const char *stage)
{
    (void)gedp;
    const char *env = getenv("QTCAD_OBOL_PROGRESSIVE_LOD_DIAGNOSTICS");
    if (!env || !BU_STR_EQUAL(env, "1") || !controller)
	return;

    SoNode *root = controller->getRenderSceneRoot();
    if (!root)
	root = controller->getSceneRoot();

    SbBox3f bbox;
    bbox.makeEmpty();
    if (root) {
	SoGetBoundingBoxAction bboxAction(controller->getViewportRegion());
	bboxAction.apply(root);
	bbox = bboxAction.getBoundingBox();
    }

    SoCamera *camera = controller->getCamera();
    if (camera) {
	const SbVec3f position = camera->position.getValue();
	fprintf(stderr,
		"qtcad_progressive_lod_diag case=%s stage=%s bbox_empty=%d bbox_min=(%.9g %.9g %.9g) bbox_max=(%.9g %.9g %.9g) camera_pos=(%.9g %.9g %.9g) near=%.9g far=%.9g focal=%.9g",
		testCase.name, stage ? stage : "",
		bbox.isEmpty() ? 1 : 0,
		bbox.getMin()[0], bbox.getMin()[1], bbox.getMin()[2],
		bbox.getMax()[0], bbox.getMax()[1], bbox.getMax()[2],
		position[0], position[1], position[2],
		camera->nearDistance.getValue(),
		camera->farDistance.getValue(),
		camera->focalDistance.getValue());
	if (camera->isOfType(SoOrthographicCamera::getClassTypeId())) {
	    SoOrthographicCamera *orthographic =
		static_cast<SoOrthographicCamera *>(camera);
	    fprintf(stderr, " ortho_height=%.9g", orthographic->height.getValue());
	} else if (camera->isOfType(SoPerspectiveCamera::getClassTypeId())) {
	    SoPerspectiveCamera *perspective =
		static_cast<SoPerspectiveCamera *>(camera);
	    fprintf(stderr, " perspective_height_angle=%.9g",
		    perspective->heightAngle.getValue());
	}
	fprintf(stderr, "\n");
    } else {
	fprintf(stderr,
		"qtcad_progressive_lod_diag case=%s stage=%s bbox_empty=%d bbox_min=(%.9g %.9g %.9g) bbox_max=(%.9g %.9g %.9g) camera=null\n",
		testCase.name, stage ? stage : "",
		bbox.isEmpty() ? 1 : 0,
		bbox.getMin()[0], bbox.getMin()[1], bbox.getMin()[2],
		bbox.getMax()[0], bbox.getMax()[1], bbox.getMax()[2]);
    }
    BObolSceneController renderScene(controller->getRenderSceneRoot());
    qtcad_obol_print_scene_source_diagnostics("render", &renderScene);
    fflush(stderr);
}


static void
set_default_env(const char *name, const char *value)
{
    if (!getenv(name))
	bu_setenv(name, value, 1);
}

static int
copy_file(const char *src, const char *dst)
{
    std::ifstream in(src, std::ios::binary);
    std::ofstream out(dst, std::ios::binary);
    if (!in || !out)
	return 0;
    out << in.rdbuf();
    return in.good() || in.eof();
}

static int
lit_pixel_count(const QImage &image)
{
    QImage rgba = image.convertToFormat(QImage::Format_RGBA8888);
    int count = 0;

    for (int y = 0; y < rgba.height(); y++) {
	const unsigned char *line = rgba.constScanLine(y);
	for (int x = 0; x < rgba.width(); x++) {
	    const unsigned char *p = line + x * 4;
	    if (p[0] > 32 || p[1] > 32 || p[2] > 32)
		count++;
	}
    }

    return count;
}

struct lit_pixel_bounds {
    int valid;
    int lit;
    int minX;
    int minY;
    int maxX;
    int maxY;
    double centerX;
    double centerY;
};

static int
lit_pixel_bounds_get(const QImage &image, struct lit_pixel_bounds *bounds)
{
    if (!bounds)
	return 0;

    memset(bounds, 0, sizeof(*bounds));
    if (image.isNull())
	return 0;

    QImage rgba = image.convertToFormat(QImage::Format_RGBA8888);
    bounds->minX = rgba.width();
    bounds->minY = rgba.height();
    bounds->maxX = -1;
    bounds->maxY = -1;

    for (int y = 0; y < rgba.height(); y++) {
	const unsigned char *line = rgba.constScanLine(y);
	for (int x = 0; x < rgba.width(); x++) {
	    const unsigned char *p = line + x * 4;
	    if (p[0] <= 32 && p[1] <= 32 && p[2] <= 32)
		continue;
	    bounds->lit++;
	    if (x < bounds->minX)
		bounds->minX = x;
	    if (y < bounds->minY)
		bounds->minY = y;
	    if (x > bounds->maxX)
		bounds->maxX = x;
	    if (y > bounds->maxY)
		bounds->maxY = y;
	}
    }

    if (bounds->lit <= 0 || bounds->maxX < bounds->minX ||
	bounds->maxY < bounds->minY)
	return 0;

    bounds->valid = 1;
    bounds->centerX = (static_cast<double>(bounds->minX) +
		       static_cast<double>(bounds->maxX)) * 0.5;
    bounds->centerY = (static_cast<double>(bounds->minY) +
		       static_cast<double>(bounds->maxY)) * 0.5;
    return 1;
}

static int
lit_pixel_bounds_centered(const QImage &image,
			  double max_center_offset_fraction,
			  struct lit_pixel_bounds *bounds)
{
    struct lit_pixel_bounds local_bounds;
    if (!bounds)
	bounds = &local_bounds;
    if (!lit_pixel_bounds_get(image, bounds))
	return 0;

    const double expectedX =
	(static_cast<double>(image.width()) - 1.0) * 0.5;
    const double expectedY =
	(static_cast<double>(image.height()) - 1.0) * 0.5;
    const double allowedX =
	static_cast<double>(image.width()) * max_center_offset_fraction;
    const double allowedY =
	static_cast<double>(image.height()) * max_center_offset_fraction;
    return std::fabs(bounds->centerX - expectedX) <= allowedX &&
	   std::fabs(bounds->centerY - expectedY) <= allowedY;
}

static int
image_byte_diff(const QImage &a, const QImage &b)
{
    if (a.isNull() || b.isNull() ||
	a.width() != b.width() ||
	a.height() != b.height())
	return -1;

    QImage ar = a.convertToFormat(QImage::Format_RGBA8888);
    QImage br = b.convertToFormat(QImage::Format_RGBA8888);
    int diff = 0;
    for (int y = 0; y < ar.height(); y++) {
	const unsigned char *ap = ar.constScanLine(y);
	const unsigned char *bp = br.constScanLine(y);
	for (int i = 0; i < ar.width() * 4; i++) {
	    if (ap[i] != bp[i])
		diff++;
	}
    }
    return diff;
}

static int
accumulate_realized_mesh_count(SoNode *node)
{
    int count = 0;

    if (!node)
	return 0;

    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return 1;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoBRLDatabaseSource *source =
	    static_cast<SoBRLDatabaseSource *>(node);
	if (source->isCompactOccurrenceRegistry()) {
	    int compactCount = 0;
	    for (int i = 0; i < source->getCompactInstanceCount(); i++) {
		BObolCompactInstanceHandle handle;
		BObolCompactInstanceSummary summary;
		if (source->getCompactInstanceHandle(i, handle) &&
		    source->getCompactInstanceSummary(handle, summary) &&
		    summary.valid && summary.meshGeometry)
		    compactCount++;
	    }
	    return compactCount;
	}
	return source->getRealizedMeshCount();
    }

    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    count += accumulate_realized_mesh_count(group->getChild(i));
    }

    return count;
}

static int
realized_mesh_count(BObolViewController *controller)
{
    return controller ?
	accumulate_realized_mesh_count(controller->getRenderSceneRoot()) : 0;
}

/* The renderer's service/frontier records are private to libged.  This
 * integration test consumes the supported `view lod` key/value contract and
 * keeps only the fields needed for its diagnostics. */
struct qtcad_obol_lod_service_status {
    int auto_submit;
    size_t in_flight;
    size_t pending_tasks;
    size_t queued_results;
    size_t queued_cache_writes;
    size_t delayed_tasks;
    unsigned int last_submitted_task_count;
    size_t active_aabb_proxy_payloads;
};

struct qtcad_obol_source_prewarm_status {
    size_t child_count;
    size_t considered;
    size_t submitted;
    size_t already_cached;
    size_t skipped_non_union;
    size_t skipped_duplicate_instance;
    size_t shared_request;
    size_t non_union_children;
    size_t duplicate_instances;
    size_t skipped_invalid;
    size_t remaining;
    size_t comb_sources;
    size_t leaf_sources;
};

struct qtcad_obol_source_expansion_status {
    size_t child_count;
    size_t considered;
    size_t expanded;
    size_t existing;
    size_t skipped_non_union;
    size_t skipped_duplicate_instance;
    size_t expanded_non_union;
    size_t expanded_duplicate_instance;
    size_t skipped_invalid;
    size_t remaining;
    size_t proxy_published;
    size_t metadata_applied;
    size_t comb_sources;
    size_t leaf_sources;
};

static size_t
qtcad_obol_result_value(const char *result, const char *key)
{
    if (!result || !key)
	return 0;
    const char *value = strstr(result, key);
    if (!value)
	return 0;
    value += strlen(key);
    if (*value != '=' && *value != ':')
	return 0;
    value++;
    while (*value == ' ')
	value++;
    return (size_t)strtoull(value, NULL, 10);
}

static int
qtcad_obol_lod_service_status_get(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    struct qtcad_obol_lod_service_status *status)
{
    (void)view_ctx;
    if (!gedp || !status)
	return 0;
    const char *av[5] = {"view", "lod", "service", "status", NULL};
    if (ged_exec_view(gedp, 4, av) != BRLCAD_OK)
	return 0;
    const char *result = bu_vls_addr(gedp->ged_result_str);
    memset(status, 0, sizeof(*status));
    status->auto_submit = (int)
	qtcad_obol_result_value(result, "auto_submit");
    status->in_flight = qtcad_obol_result_value(result, "in_flight");
    status->pending_tasks = qtcad_obol_result_value(result, "pending_tasks");
    status->queued_results = qtcad_obol_result_value(result, "queued_results");
    status->queued_cache_writes =
	qtcad_obol_result_value(result, "queued_cache_writes");
    status->delayed_tasks = qtcad_obol_result_value(result, "delayed_tasks");
    status->last_submitted_task_count = (unsigned int)
	qtcad_obol_result_value(result, "last_submitted_tasks");
    status->active_aabb_proxy_payloads =
	qtcad_obol_result_value(result, "active_lod_aabb_proxies");
    return 1;
}

static size_t
qtcad_obol_frontier_prewarm(
    struct ged *gedp,
    const char *path,
    int draw_mode,
    size_t max_sources,
    size_t max_children,
    struct qtcad_obol_source_prewarm_status *status)
{
    char mode[32] = {0};
    char sources[32] = {0};
    char children[32] = {0};
    snprintf(mode, sizeof(mode), "%d", draw_mode);
    snprintf(sources, sizeof(sources), "%zu", max_sources);
    snprintf(children, sizeof(children), "%zu", max_children);
    const char *av[9] = {"view", "lod", "frontier", "prewarm", path,
	mode, sources, children, NULL};
    if (!status || ged_exec_view(gedp, 8, av) != BRLCAD_OK)
	return 0;
    const char *result = bu_vls_addr(gedp->ged_result_str);
    memset(status, 0, sizeof(*status));
#define QTCAD_PREWARM_VALUE(_field) \
    status->_field = qtcad_obol_result_value(result, #_field)
    QTCAD_PREWARM_VALUE(child_count);
    QTCAD_PREWARM_VALUE(considered);
    QTCAD_PREWARM_VALUE(submitted);
    QTCAD_PREWARM_VALUE(already_cached);
    QTCAD_PREWARM_VALUE(skipped_non_union);
    QTCAD_PREWARM_VALUE(skipped_duplicate_instance);
    QTCAD_PREWARM_VALUE(shared_request);
    QTCAD_PREWARM_VALUE(non_union_children);
    QTCAD_PREWARM_VALUE(duplicate_instances);
    QTCAD_PREWARM_VALUE(skipped_invalid);
    QTCAD_PREWARM_VALUE(remaining);
    QTCAD_PREWARM_VALUE(comb_sources);
    QTCAD_PREWARM_VALUE(leaf_sources);
#undef QTCAD_PREWARM_VALUE
    return qtcad_obol_result_value(result, "result");
}

static int
qtcad_obol_frontier_expand(
    struct ged *gedp,
    const char *path,
    int draw_mode,
    size_t max_sources,
    size_t max_children,
    struct qtcad_obol_source_expansion_status *status)
{
    char mode[32] = {0};
    char sources[32] = {0};
    char children[32] = {0};
    snprintf(mode, sizeof(mode), "%d", draw_mode);
    snprintf(sources, sizeof(sources), "%zu", max_sources);
    snprintf(children, sizeof(children), "%zu", max_children);
    const char *av[9] = {"view", "lod", "frontier", "expand", path,
	mode, sources, children, NULL};
    if (!status || ged_exec_view(gedp, 8, av) != BRLCAD_OK)
	return 0;
    const char *result = bu_vls_addr(gedp->ged_result_str);
    memset(status, 0, sizeof(*status));
#define QTCAD_EXPAND_VALUE(_field) \
    status->_field = qtcad_obol_result_value(result, #_field)
    QTCAD_EXPAND_VALUE(child_count);
    QTCAD_EXPAND_VALUE(considered);
    QTCAD_EXPAND_VALUE(expanded);
    QTCAD_EXPAND_VALUE(existing);
    QTCAD_EXPAND_VALUE(skipped_non_union);
    QTCAD_EXPAND_VALUE(skipped_duplicate_instance);
    QTCAD_EXPAND_VALUE(expanded_non_union);
    QTCAD_EXPAND_VALUE(expanded_duplicate_instance);
    QTCAD_EXPAND_VALUE(skipped_invalid);
    QTCAD_EXPAND_VALUE(remaining);
    QTCAD_EXPAND_VALUE(proxy_published);
    QTCAD_EXPAND_VALUE(metadata_applied);
    QTCAD_EXPAND_VALUE(comb_sources);
    QTCAD_EXPAND_VALUE(leaf_sources);
#undef QTCAD_EXPAND_VALUE
    return (int)qtcad_obol_result_value(result, "result");
}

static int
qtcad_obol_wait_for_lod_service_idle(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    int timeout_ms,
    qtcad_obol_lod_service_status *status)
{
    qtcad_obol_lod_service_status local_status;
    if (!status)
	status = &local_status;

    for (int elapsed = 0; elapsed <= timeout_ms; elapsed += 25) {
	/* Cache prewarm runs entirely in the service; pumping the view here can
	 * re-enter rendering before the bounded cache request has settled. */
	memset(status, 0, sizeof(*status));
	if (!qtcad_obol_lod_service_status_get(gedp, view_ctx, status))
	    return 0;
	if (status->in_flight == 0 &&
	    status->pending_tasks == 0 &&
	    status->queued_cache_writes == 0 &&
	    status->delayed_tasks == 0)
	    return 1;
	std::this_thread::sleep_for(std::chrono::milliseconds(25));
    }

    return 0;
}

static void
qtcad_obol_prewarm_status_accumulate(
    struct qtcad_obol_source_prewarm_status *dst,
    const struct qtcad_obol_source_prewarm_status *src)
{
    if (!dst || !src)
	return;
    dst->child_count += src->child_count;
    dst->considered += src->considered;
    dst->submitted += src->submitted;
    dst->already_cached += src->already_cached;
    dst->skipped_non_union += src->skipped_non_union;
    dst->skipped_duplicate_instance += src->skipped_duplicate_instance;
    dst->skipped_invalid += src->skipped_invalid;
    dst->remaining += src->remaining;
    dst->comb_sources += src->comb_sources;
    dst->leaf_sources += src->leaf_sources;
}

static void
qtcad_obol_expansion_status_accumulate(
    struct qtcad_obol_source_expansion_status *dst,
    const struct qtcad_obol_source_expansion_status *src)
{
    if (!dst || !src)
	return;
    dst->child_count += src->child_count;
    dst->considered += src->considered;
    dst->expanded += src->expanded;
    dst->existing += src->existing;
    dst->skipped_non_union += src->skipped_non_union;
    dst->skipped_duplicate_instance += src->skipped_duplicate_instance;
    dst->skipped_invalid += src->skipped_invalid;
    dst->remaining += src->remaining;
    dst->proxy_published += src->proxy_published;
    dst->metadata_applied += src->metadata_applied;
    dst->comb_sources += src->comb_sources;
    dst->leaf_sources += src->leaf_sources;
}

static int
process_until_proxy_kind(BObolViewController *controller,
			 BObolLodService &service,
			 int expected_kind,
			 size_t min_active_payload_count,
			 size_t &applied_total,
			 int timeout_ms,
			 const char *label)
{
    for (int elapsed = 0; elapsed <= timeout_ms; elapsed += 25) {
	QCoreApplication::processEvents();
	(void)controller->processPendingLodResults(1);
	applied_total += controller->getLastLodAppliedResultCount();
	if (controller->getLastLodRejectedResultCount() != 0) {
	    fprintf(stderr, "qtcad Obol progressive LoD rejected %s results: %s\n",
		    label, controller->getLastLodDiagnostics().getString());
	    return 0;
	}
	if (controller->getActiveLodProxyPayloadCount(expected_kind) >=
	    min_active_payload_count)
	    return 1;
	if (service.inFlightCount() == 0 &&
	    service.queuedResultCountForDiagnostics() == 0 &&
	    service.pendingTaskCountForDiagnostics() == 0)
	    break;
	std::this_thread::sleep_for(std::chrono::milliseconds(25));
    }

    fprintf(stderr, "qtcad Obol progressive LoD did not reach %s proxy kind %d with %zu active payloads, current=%d active=%zu applied=%zu pending=%zu queued=%zu in_flight=%zu\n",
	    label, expected_kind, min_active_payload_count,
	    0,
	    controller->getActiveLodProxyPayloadCount(expected_kind), applied_total,
	    service.pendingTaskCountForDiagnostics(),
	    service.queuedResultCountForDiagnostics(),
	    service.inFlightCount());
    return 0;
}

static int
process_until_mesh(BObolViewController *controller,
		   BObolLodService &service,
		   size_t min_active_payload_count,
		   size_t &applied_total,
		   int timeout_ms)
{
    for (int elapsed = 0; elapsed <= timeout_ms; elapsed += 25) {
	QCoreApplication::processEvents();
	(void)controller->processPendingLodResults(8);
	applied_total += controller->getLastLodAppliedResultCount();
	if (controller->getLastLodRejectedResultCount() != 0) {
	    fprintf(stderr, "qtcad Obol progressive LoD rejected mesh results: %s\n",
		    controller->getLastLodDiagnostics().getString());
	    return 0;
	}
	if (controller->getActiveLodMeshPayloadCount() >=
	    min_active_payload_count)
	    return 1;
	if (service.inFlightCount() == 0 &&
	    service.queuedResultCountForDiagnostics() == 0 &&
	    service.pendingTaskCountForDiagnostics() == 0)
	    break;
	std::this_thread::sleep_for(std::chrono::milliseconds(25));
    }

    fprintf(stderr, "qtcad Obol progressive LoD did not reach %zu active mesh payloads, active=%zu applied=%zu pending=%zu queued=%zu in_flight=%zu\n",
	    min_active_payload_count, controller->getActiveLodMeshPayloadCount(),
	    applied_total, service.pendingTaskCountForDiagnostics(),
	    service.queuedResultCountForDiagnostics(),
	    service.inFlightCount());
    return 0;
}

static int
run_ged_command(struct ged *gedp, int argc, const char **argv,
		int (*func)(struct ged *, int, const char **),
		const char *label)
{
    int ret = (*func)(gedp, argc, argv);
    if (ret == BRLCAD_OK)
	return 1;

    fprintf(stderr, "%s failed: %s\n", label,
	    bu_vls_cstr(gedp->ged_result_str));
    return 0;
}

static int
create_repeated_draw_target(struct ged *gedp,
			    const char *sourceTarget,
			    int copyCount,
			    char *drawTarget,
			    size_t drawTargetLen)
{
    std::vector<std::string> copyNames;
    std::vector<const char *> groupArgv;

    if (!gedp || !sourceTarget || !drawTarget || drawTargetLen == 0)
	return 0;

    if (copyCount <= 1) {
	snprintf(drawTarget, drawTargetLen, "%s", sourceTarget);
	return 1;
    }

    snprintf(drawTarget, drawTargetLen, "lod_stress_group");
    copyNames.reserve((size_t)copyCount - 1);
    groupArgv.reserve((size_t)copyCount + 2);
    groupArgv.push_back("g");
    groupArgv.push_back(drawTarget);
    groupArgv.push_back(sourceTarget);

    for (int i = 1; i < copyCount; i++) {
	char name[MAXPATHLEN] = {0};
	snprintf(name, sizeof(name), "lod_stress_%03d.bot", i);
	copyNames.push_back(name);

	const char *copyArgv[4] = {"cp", sourceTarget,
				   copyNames.back().c_str(), NULL
				  };
	if (!run_ged_command(gedp, 3, copyArgv, ged_exec_cp,
			     "stress copy")) {
	    return 0;
	}
	groupArgv.push_back(copyNames.back().c_str());
    }

    if (!run_ged_command(gedp, (int)groupArgv.size(), groupArgv.data(),
			 ged_exec_g, "stress group"))
	return 0;

    ged_db_index_refresh(gedp);
    ged_event_notify_batch_rebuild(gedp, NULL);
    return 1;
}

static int
run_progressive_lod_case(const struct progressive_lod_case &testCase)
{
    int64_t total_start = bu_gettime();
    int64_t phase_start = 0;
    char case_name[MAXPATHLEN] = {0};
    std::string tmp_db;
    std::string cache_leaf;
    char cache_dir[MAXPATHLEN] = {0};
    char active_draw_target[MAXPATHLEN] = {0};
    struct ged *gedp = NULL;
    BObolViewController *controller = NULL;
    QgView view(NULL, QgViewType::SW);
    BObolLodService service;
    int service_started = 0;
    int endpoint_attached = 0;
    struct ged_scene_result *draw_result = NULL;
    struct progressive_lod_timings timings = {};
    const int startup_deferred =
	testCase.startupOnly || testCase.startupRefine ||
	testCase.startupExpand || testCase.startupAutoExpand ||
	testCase.startupWireLod;
    const int production_auto_expand =
	testCase.startupAutoExpand || testCase.startupWireLod;
    struct BObolDrawMetadataRecord startup_metadata;
    bobol_draw_metadata_record_init(&startup_metadata);

    auto cleanup = [&]() {
	if (controller) {
	    controller->setLodAutoSubmit(FALSE);
	    if (controller->getLodService() == &service)
		controller->setLodService(NULL);
	}
	if (service_started)
	    service.stop();
	ged_scene_result_destroy(draw_result);
	if (endpoint_attached)
	    (void)ged_view_context_obol_endpoint_set(
		ged_view_context_from_bv(view.viewContext()),
		NULL, 0);
	if (gedp)
	    ged_close(gedp);
	bu_file_delete(tmp_db.c_str());
	bu_dirclear(cache_dir);
    };

    sanitize_case_name(case_name, sizeof(case_name), testCase.name);
    tmp_db = "qtcad_obol_progressive_lod_";
    tmp_db += case_name;
    tmp_db += "_";
    tmp_db += std::to_string(bu_pid());
    tmp_db += ".g";

    cache_leaf = "qtcad_obol_progressive_lod_";
    cache_leaf += case_name;
    cache_leaf += "_";
    cache_leaf += std::to_string(bu_pid());
    cache_leaf += "_cache";
    bu_dir(cache_dir, MAXPATHLEN, BU_DIR_CURR, cache_leaf.c_str(), NULL);
    bu_dirclear(cache_dir);
    bu_mkdir(cache_dir);
    bu_setenv("BU_DIR_CACHE", cache_dir, 1);
    /*
     * Normal delivery coalesces completed stages so first presentation uses
     * the best available payload.  Delays are exclusively for the explicit
     * diagnostic test that verifies each publication boundary.
     */
    if (testCase.diagnosticStages) {
	/* An explicit diagnostic request must be independent of its caller. */
	bu_setenv("BOBOL_LOD_OBB_TASK_DELAY_MS",
		  testCase.slowLodDelays ? "350" : "75", 1);
	bu_setenv("BOBOL_LOD_TASK_DELAY_MS",
		  testCase.slowLodDelays ? "700" : "150", 1);
    } else {
	set_default_env("BOBOL_LOD_OBB_TASK_DELAY_MS", "0");
	set_default_env("BOBOL_LOD_TASK_DELAY_MS", "0");
    }

    phase_start = bu_gettime();
    bu_file_delete(tmp_db.c_str());
    if (!copy_file(testCase.sourceDb, tmp_db.c_str())) {
	fprintf(stderr, "unable to copy %s to %s\n", testCase.sourceDb, tmp_db.c_str());
	cleanup();
	return 0;
    }
    record_progressive_lod_phase(timings.copySeconds, phase_start,
				 total_start, testCase, "copy");

    phase_start = bu_gettime();
    gedp = ged_open("db", tmp_db.c_str(), 1);
    if (!gedp) {
	fprintf(stderr, "failed to open %s\n", tmp_db.c_str());
	cleanup();
	return 0;
    }
    record_progressive_lod_phase(timings.openSeconds, phase_start,
				 total_start, testCase, "open");

    phase_start = bu_gettime();
    view.resize(512, 512);
	struct ged_view_context *ged_view_ctx =
	    ged_view_context_from_bv(view.viewContext());
	ged_view_active_ctx_set(gedp, ged_view_ctx);
	if (!ged_view_set_context_add(ged_view_set_ctx(gedp), ged_view_ctx) ||
	    !ged_view_context_host_attach(gedp, ged_view_ctx)) {
	    fprintf(stderr, "failed to register progressive LoD view context\n");
	    cleanup();
	    return 0;
	}
	(void)bv_unit_conversion_set(bv_context_view(
	    view.viewContext()),
	    gedp->dbip->dbi_local2base, gedp->dbip->dbi_base2local);
    record_progressive_lod_phase(timings.viewSetupSeconds, phase_start,
				 total_start, testCase, "view_setup");

    {
	const char *av[5] = {"view", "lod", "cache", "clear", NULL};
	phase_start = bu_gettime();
	if (!run_ged_command(gedp, 4, av, ged_exec_view,
			     "view lod cache clear")) {
	    cleanup();
	    return 0;
	}
	record_progressive_lod_phase(timings.cacheClearSeconds, phase_start,
				     total_start, testCase, "cache_clear");
    }
    {
	const char *av[4] = {"tol", "rel", "0.0002", NULL};
	phase_start = bu_gettime();
	if (!run_ged_command(gedp, 3, av, ged_exec_tol, "tol rel")) {
	    cleanup();
	    return 0;
	}
	record_progressive_lod_phase(timings.tolSeconds, phase_start,
				     total_start, testCase, "tol");
    }
    if (testCase.useFacetize) {
	const char *av[5] = {"facetize", "-r", testCase.facetizeRoot,
			     testCase.drawTarget, NULL
			    };
	phase_start = bu_gettime();
	if (!run_ged_command(gedp, 4, av, ged_exec_facetize,
			     "facetize")) {
	    cleanup();
	    return 0;
	}
	ged_db_index_refresh(gedp);
	ged_event_notify_batch_rebuild(gedp, NULL);
	record_progressive_lod_phase(timings.facetizeSeconds, phase_start,
				     total_start, testCase, "facetize");
    }
    phase_start = bu_gettime();
    if (!create_repeated_draw_target(gedp, testCase.drawTarget,
				     testCase.copyCount, active_draw_target,
				     sizeof(active_draw_target))) {
	cleanup();
	return 0;
    }
    record_progressive_lod_phase(timings.repeatTargetSeconds, phase_start,
				 total_start, testCase, "repeat_target");
    phase_start = bu_gettime();
    {
	const char *av[5] = {"view", "lod", "mesh", "1", NULL};
	if (!run_ged_command(gedp, 4, av, ged_exec_view,
			     "view lod mesh 1")) {
	    cleanup();
	    return 0;
	}
    }
    {
	const char *av[5] = {"view", "lod", "scale", "0.8", NULL};
	if (!run_ged_command(gedp, 4, av, ged_exec_view,
			     "view lod scale")) {
	    cleanup();
	    return 0;
	}
    }
    record_progressive_lod_phase(timings.lodPolicySeconds, phase_start,
				 total_start, testCase, "lod_policy");

    if (startup_deferred) {
	startup_metadata.directoryFound = 1;
	startup_metadata.isComb = 1;
	startup_metadata.hasRegionId = 1;
	startup_metadata.regionId = 123456;
	startup_metadata.hasAircode = 1;
	startup_metadata.aircode = 7;
	startup_metadata.hasLos = 1;
	startup_metadata.los = 89;
	startup_metadata.hasMaterialId = 1;
	startup_metadata.materialId = 654321;
	startup_metadata.hasShader = 1;
	bu_strlcpy(startup_metadata.shader, "cache-startup-shader",
		   sizeof(startup_metadata.shader));
	if (bobol_draw_path_metadata_cache_store(gedp->dbip,
		active_draw_target, &startup_metadata, NULL) != BRLCAD_OK) {
	    fprintf(stderr,
		    "qtcad progressive LoD startup metadata seed failed for %s\n",
		    active_draw_target);
	    cleanup();
	    return 0;
	}
    }

    controller = view.obolViewController();
    if (!controller) {
	fprintf(stderr, "QgView did not provide an Obol controller\n");
	cleanup();
	return 0;
    }
    controller->clearDatabaseSources();
    if (!ged_view_context_obol_endpoint_set(
	    ged_view_context_from_bv(view.viewContext()),
	    view.displayEndpoint(), 0)) {
	fprintf(stderr, "failed to attach qtcad display endpoint for GED draw\n");
	cleanup();
	return 0;
    }
    endpoint_attached = 1;
    /* The production automatic fixture mirrors `e -m1`: begin with compact
     * leaf boxes, shade BoTs through view-aware LoD payloads, and use ordinary
     * geometry for the remaining leaves. */
    enum ged_scene_draw_mode startup_draw_mode = testCase.startupAutoExpand ?
			   GED_SCENE_DRAW_SHADED_BOTS :
			   ((testCase.startupOnly || testCase.startupExpand ||
			     testCase.startupWireLod) ?
			    GED_SCENE_DRAW_WIRE : GED_SCENE_DRAW_SHADED);

    const char *draw_path = active_draw_target;
    struct ged_scene_draw_request draw_request;
    ged_scene_draw_request_init(&draw_request);
    draw_request.view = ged_view_context_from_bv(view.viewContext());
    draw_request.paths = &draw_path;
    draw_request.path_count = 1;
    draw_request.style.draw_mode = startup_draw_mode;
    draw_request.realization.mode = startup_deferred ?
	GED_SCENE_REALIZE_PROGRESSIVE : GED_SCENE_REALIZE_EAGER;
    draw_request.autoview = production_auto_expand ? 1 : 0;

    int sync_changed = 0;
    int render_source_count = 0;
    int draw_ret = 0;
    if (testCase.shadedLod) {
	BObolSceneController *scene = controller->getSceneController();
	if (!scene) {
	    fprintf(stderr,
		    "qtcad Obol progressive LoD shaded-lod has no scene controller\n");
	    cleanup();
	    return 0;
	}
	controller->clearDatabaseSources();
	controller->clearViewLodState();
	const uint32_t source_revision = 1;
	const char *root_instance_key = "shaded_lod_root";
	phase_start = bu_gettime();
	if (scene->replaceDatabaseSourceInstanceRepresentation(
		root_instance_key, active_draw_target,
		root_instance_key,
		SoBRLDatabaseSource::REPRESENTATION_SHADED_BOTS,
		gedp->dbip, SoBRLDatabaseSource::SHADED,
		source_revision) < 0 ||
	    scene->setDatabaseSourceInstanceRealizationRoleFlags(
		root_instance_key,
		SoBRLDatabaseSource::REALIZATION_ROLE_MESH) < 0 ||
	    scene->setDatabaseSourceInstanceRealizationViewPolicy(
		root_instance_key, TRUE, FALSE, TRUE, 1.0f, 1.0f,
		512, 512, 1, 1.0f, 1.0f) < 0 ||
	    !controller->realizePending()) {
	    record_progressive_lod_phase(timings.drawApplySeconds,
					 phase_start, total_start,
					 testCase, "obol_shaded_realize");
	    fprintf(stderr,
		    "qtcad Obol progressive LoD shaded-lod Obol realization failed: diagnostics=%s\n",
		    controller->getLastDiagnostics().getString());
	    cleanup();
	    return 0;
	}
	record_progressive_lod_phase(timings.drawApplySeconds,
				     phase_start, total_start, testCase,
				     "obol_shaded_realize");
	render_source_count = controller->getDatabaseSourceCount();
    } else {
	draw_result = ged_scene_result_create();
	phase_start = bu_gettime();
	draw_ret = ged_scene_draw(gedp, &draw_request, draw_result) ==
	    GED_SCENE_OK ? 1 : -1;
	record_progressive_lod_phase(timings.drawApplySeconds, phase_start,
				     total_start, testCase, "draw_apply");
	render_source_count =
	    qtcad_obol_render_database_source_count(gedp, controller);
    }
    if (!testCase.shadedLod) {
	phase_start = bu_gettime();
	sync_changed = ged_scene_result_changed(draw_result);
	if (sync_changed)
	    view.need_update(QG_VIEW_REFRESH);
	record_progressive_lod_phase(timings.syncSeconds, phase_start,
				     total_start, testCase, "obol_sync");
	render_source_count =
	    qtcad_obol_render_database_source_count(gedp, controller);
    } else if (phase_logging_enabled(testCase)) {
	fprintf(stderr,
		"qtcad_progressive_lod_phase case=%s phase=obol_sync seconds=0.000 total=%.3f skipped=already_synced\n",
		testCase.name, elapsed_seconds(total_start));
	fflush(stderr);
    }
    if (draw_ret < 0 || render_source_count <= 0) {
	fprintf(stderr, "qtcad progressive LoD draw/sync failed: draw_ret=%d sync_changed=%d source_count=%d errors=%s\n",
		draw_ret, sync_changed, render_source_count,
		ged_scene_result_diagnostic(draw_result));
	cleanup();
	return 0;
    }
    int startup_realized_geometry_count = -1;
    if (startup_deferred) {
	startup_realized_geometry_count =
	    qtcad_obol_render_realized_database_geometry_count(gedp,
		controller);
	if (!qtcad_obol_render_has_draw_metadata(gedp, controller,
		&startup_metadata)) {
	    fprintf(stderr,
		    "qtcad progressive LoD startup proxy did not publish cached path metadata\n");
	    BObolSceneController renderScene(
		controller->getRenderSceneRoot());
	    qtcad_obol_print_scene_source_diagnostics("render",
		&renderScene);
	    cleanup();
	    return 0;
	}
    }
    ged_scene_result_destroy(draw_result);
    draw_result = NULL;
    qtcad_obol_print_render_diagnostics(gedp, testCase, controller,
					"before-view-all");

    phase_start = bu_gettime();
    /*
     * A production automatic draw has already autoviewed the best bound
     * available in ged_draw_apply_transaction() and armed the progressive
     * follow operation.  Re-autoviewing the same partial scene here advances
     * the bv frame revision and correctly tells that follow operation that a
     * user changed the view.  On a cold draw this left the test framed around
     * the first small proxy while authoritative leaf bounds arrived elsewhere
     * in model space, so frustum-aware LoD submission rejected every leaf.
     *
     * Preserve the transaction's view revision in that case and only mirror
     * its current bv state into the hidden canvas camera.  Other test modes do
     * not arm progressive follow and still need the explicit scene autoview.
     */
    const int initial_autoview = production_auto_expand ?
	(controller->syncCameraFromViewContext(
	    ged_view_context_from_bv(view.viewContext())) ? 1 : 0) :
	qtcad_obol_scene_autoview_refresh(view, controller,
	    "progressive-lod-view-all");
    /* A cold asynchronous draw may not have a publishable bound in its first
     * UI tick.  Its transaction already armed progressive autoview, which will
     * apply as soon as a region box or streamed leaf arrives. */
    if (!initial_autoview && !production_auto_expand) {
	fprintf(stderr, "qtcad progressive LoD autoview failed: case=%s\n",
	    testCase.name);
	cleanup();
	return 0;
    }
    record_progressive_lod_phase(timings.viewAllSeconds, phase_start,
				 total_start, testCase, "view_all");
    qtcad_obol_print_render_diagnostics(gedp, testCase, controller,
					"after-view-all");

    int startup_auto_source_count_before = -1;
    int startup_auto_depth_before = -1;
    int startup_auto_geometry_before = -1;
    int startup_auto_compact_before = -1;
    if (production_auto_expand) {
	startup_auto_source_count_before =
	    qtcad_obol_render_database_source_count(gedp, controller);
	startup_auto_depth_before =
	    qtcad_obol_render_database_source_max_depth(gedp, controller,
		active_draw_target);
	startup_auto_geometry_before =
	    qtcad_obol_render_realized_database_geometry_count(gedp,
		controller);
	startup_auto_compact_before =
	    qtcad_obol_render_compact_instance_count(gedp, controller);
    }

    QImage frame0;
    phase_start = bu_gettime();
    qtcad_obol_request_view_frame(view, controller,
	"progressive-lod-initial-frame");
    qtcad_obol_present_view_frame(view, controller, frame0);
    record_progressive_lod_phase(timings.firstReadbackSeconds, phase_start,
				 total_start, testCase, "readback0");
    timings.firstVisualSeconds = elapsed_seconds(total_start);

    if (startup_deferred) {
	int lit0 = lit_pixel_count(frame0);
	const int startup_has_first_visual =
	    startup_realized_geometry_count > 0 && !frame0.isNull() &&
	    lit0 >= 20;
	if (!startup_has_first_visual && testCase.startupOnly) {
	    fprintf(stderr,
		    "qtcad_progressive_lod_startup case=%s lit=%d sources=%d geometry=%d",
		    testCase.name, lit0, render_source_count,
		    startup_realized_geometry_count);
	    print_progressive_lod_timings(stderr, timings);
	    fprintf(stderr, "\n");
	    cleanup();
	    return 1;
	}
	if (!startup_has_first_visual && !testCase.startupRefine &&
	    !testCase.startupExpand &&
	    !testCase.startupAutoExpand && !testCase.startupWireLod) {
	    fprintf(stderr,
		    "qtcad Obol progressive LoD startup image check failed: case=%s lit=%d",
		    testCase.name, lit0);
	    print_progressive_lod_timings(stderr, timings);
	    fprintf(stderr, "\n");
	    cleanup();
	    return 0;
	}
	fprintf(stderr,
		"qtcad_progressive_lod_startup case=%s lit=%d sources=%d geometry=%d",
		testCase.name, lit0, render_source_count,
		startup_realized_geometry_count);
	print_progressive_lod_timings(stderr, timings);
	fprintf(stderr, "\n");
	if (testCase.startupOnly) {
	    cleanup();
	    return 1;
	}
	if (production_auto_expand) {
	    int source_count_before = startup_auto_source_count_before;
	    int max_depth_before = startup_auto_depth_before;
	    int geometry_before = startup_auto_geometry_before;
	    int compact_before = startup_auto_compact_before;
	    int source_count_after = source_count_before;
	    int max_depth_after = max_depth_before;
	    int realized_after = geometry_before;
	    int compact_after = compact_before;
	    size_t lod_mesh_after =
		controller->getActiveLodMeshPayloadCount();
	    size_t visible_mesh_demand_after =
		qtcad_obol_render_visible_compact_mesh_demand_count(
		    controller, 0);
	    size_t all_visible_mesh_demand_after =
		qtcad_obol_render_visible_compact_mesh_demand_count(
		    controller, 1);
	    size_t unique_visible_mesh_demand_after =
		qtcad_obol_render_visible_compact_mesh_demand_count(
		    controller, 0, 1);
	    size_t projected_mesh_demand_after =
		qtcad_obol_render_visible_compact_mesh_demand_count(
		    controller, 0, 1, 1);
	    int lit_auto = 0;
	    int diff_auto = 0;
	    QImage frame_auto;
	    BObolProgressiveStatus progressive_status;
	    qtcad_obol_lod_service_status service_status;
	    memset(&service_status, 0, sizeof(service_status));
	    const int settle_timeout_ms =
		std::max(testCase.meshTimeoutMs, 90000);
	    const int64_t settle_deadline =
		bu_gettime() + static_cast<int64_t>(settle_timeout_ms) * 1000;
	    int64_t last_progress_time = bu_gettime();
	    int64_t last_feedback_frame_time = 0;
	    size_t last_progress_mesh_count = lod_mesh_after;
	    int last_progress_compact_count = compact_after;
	    size_t feedback_frame_count = 0;
	    while (bu_gettime() < settle_deadline) {
		QCoreApplication::processEvents();
		(void)controller->advanceProgressiveWork(NULL,
		    &progressive_status);
		source_count_after =
		    qtcad_obol_render_database_source_count(gedp,
			controller);
		max_depth_after =
		    qtcad_obol_render_database_source_max_depth(gedp,
			controller, active_draw_target);
		realized_after =
		    qtcad_obol_render_realized_database_geometry_count(gedp,
			controller);
		compact_after =
		    qtcad_obol_render_compact_instance_count(gedp, controller);
		lod_mesh_after = controller->getActiveLodMeshPayloadCount();
		visible_mesh_demand_after =
		    qtcad_obol_render_visible_compact_mesh_demand_count(
			controller, 0);
		all_visible_mesh_demand_after =
		    qtcad_obol_render_visible_compact_mesh_demand_count(
			controller, 1);
		unique_visible_mesh_demand_after =
		    qtcad_obol_render_visible_compact_mesh_demand_count(
			controller, 0, 1);
		projected_mesh_demand_after =
		    qtcad_obol_render_visible_compact_mesh_demand_count(
			controller, 0, 1, 1);
		(void)qtcad_obol_lod_service_status_get(gedp,
		    ged_view_context_from_bv(view.viewContext()),
		    &service_status);
		const int64_t settle_now = bu_gettime();
		if (lod_mesh_after != last_progress_mesh_count ||
		    compact_after != last_progress_compact_count) {
		    last_progress_mesh_count = lod_mesh_after;
		    last_progress_compact_count = compact_after;
		    last_progress_time = settle_now;
		}
		const bool service_idle =
		    service_status.pending_tasks == 0 &&
		    service_status.in_flight == 0 &&
		    service_status.queued_results == 0 &&
		    service_status.queued_cache_writes == 0 &&
		    !controller->hasPendingLodResults() &&
		    !controller->hasPendingLodSubmissions();
		const bool otherwise_settled =
		    compact_after > compact_before &&
		    projected_mesh_demand_after > 0 &&
		    lod_mesh_after >= projected_mesh_demand_after &&
		    realized_after >= geometry_before &&
		    service_idle;
		/* A hidden software canvas has no automatic paint delivery.
		 * The production controller deliberately waits for presentation
		 * feedback before admitting a richer cut.  Present only when that
		 * gate is explicit, or when an otherwise idle provider has made no
		 * structural progress for 100 ms.  This exercises the real contract
		 * without turning every 5 ms polling iteration into an OSMesa
		 * render/readback. */
		const bool progress_stalled =
		    settle_now - last_progress_time >= 100000;
		const bool feedback_rate_limited =
		    last_feedback_frame_time > 0 &&
		    settle_now - last_feedback_frame_time < 100000;
		const bool presentation_gate =
		    controller->hasPendingLodRefinementFrame() ||
		    controller->isLodInteractionActive();
		/*
		 * A production canvas presents while camera input is active.  A
		 * hidden test canvas must do the same even when the bounded worker
		 * queue is non-idle: the controller may intentionally retain
		 * richer immutable results until that presentation completes.
		 * Requiring service_idle here made the harness wait forever on
		 * the very result queue whose publication was waiting for a frame.
		 */
		if (progressive_status.hasMore &&
		    !feedback_rate_limited &&
		    (controller->hasPendingLodRefinementFrame() ||
		     (progress_stalled &&
		      (service_idle || presentation_gate)))) {
		    QImage feedback_frame;
		    qtcad_obol_request_view_frame(view, controller,
			"progressive-lod-feedback-frame");
		    qtcad_obol_present_view_frame(view, controller,
			feedback_frame);
		    last_feedback_frame_time = bu_gettime();
		    last_progress_time = last_feedback_frame_time;
		    feedback_frame_count++;
		}
		if (otherwise_settled && !progressive_status.hasMore)
		    break;
		std::this_thread::sleep_for(std::chrono::milliseconds(5));
	    }
	    /* Read the framebuffer once after the structural/service state has
	     * settled.  Repeated readbacks each take substantially longer than a
	     * progression pump and used to turn this into an artificial
	     * throttling test rather than a production-path convergence test. */
	    view.get_viewport_image(frame_auto);
	    lit_auto = lit_pixel_count(frame_auto);
	    diff_auto = image_byte_diff(frame0, frame_auto);
	    (void)qtcad_obol_lod_service_status_get(gedp,
		ged_view_context_from_bv(view.viewContext()), &service_status);
	    if (!service_status.auto_submit) {
		fprintf(stderr,
			"qtcad Obol progressive LoD startup did not enable the shared automatic LoD pipeline: case=%s auto_submit=%d last_submitted=%u\n",
			testCase.name, service_status.auto_submit,
			service_status.last_submitted_task_count);
		cleanup();
		return 0;
	    }
	    if (compact_after <= compact_before ||
		projected_mesh_demand_after == 0 ||
		lod_mesh_after < projected_mesh_demand_after ||
		progressive_status.hasMore ||
		realized_after < geometry_before ||
		frame_auto.isNull() || lit_auto < 20 ||
		(diff_auto <= 0 && lit0 < 20)) {
		fprintf(stderr,
			"qtcad Obol progressive LoD startup-auto-expand failed: case=%s sources=%d->%d depth=%d->%d geometry=%d->%d compact=%d->%d lod_mesh=%zu/%zu projected (%zu drawable, %zu unique, %zu all) lit=%d diff=%d feedback_frames=%zu last_submitted=%u active_aabb=%zu pending=%zu in_flight=%zu cache_writes=%zu progressive={submitted=%zu cached=%zu expanded=%zu existing=%zu remaining=%zu proxies=%zu changed=%d more=%d}",
			testCase.name, source_count_before, source_count_after,
			max_depth_before, max_depth_after,
			geometry_before, realized_after,
			compact_before, compact_after,
			lod_mesh_after,
			projected_mesh_demand_after,
			visible_mesh_demand_after,
			unique_visible_mesh_demand_after,
			all_visible_mesh_demand_after,
			lit_auto, diff_auto,
			feedback_frame_count,
			service_status.last_submitted_task_count,
			service_status.active_aabb_proxy_payloads,
			service_status.pending_tasks,
			service_status.in_flight,
			service_status.queued_cache_writes,
			progressive_status.submitted,
			progressive_status.alreadyCached,
			progressive_status.expanded,
			progressive_status.existing,
			progressive_status.remaining,
			progressive_status.proxyPublished,
			progressive_status.changed,
			progressive_status.hasMore);
		print_progressive_lod_timings(stderr, timings);
		fprintf(stderr, "\n");
		qtcad_obol_print_render_diagnostics(gedp, testCase, controller,
		    "startup-auto-expand-failure");
		cleanup();
		return 0;
	    }
	    fprintf(stderr,
		    "qtcad_progressive_lod_startup_auto_expand case=%s sources=%d->%d depth=%d->%d geometry=%d->%d compact=%d->%d lod_mesh=%zu lit=%d diff=%d feedback_frames=%zu last_submitted=%u active_aabb=%zu pending=%zu in_flight=%zu cache_writes=%zu",
		    testCase.name, source_count_before, source_count_after,
		    max_depth_before, max_depth_after,
		    geometry_before, realized_after,
		    compact_before, compact_after,
		    lod_mesh_after,
		    lit_auto, diff_auto,
		    feedback_frame_count,
		    service_status.last_submitted_task_count,
		    service_status.active_aabb_proxy_payloads,
		    service_status.pending_tasks, service_status.in_flight,
		    service_status.queued_cache_writes);
	    print_progressive_lod_timings(stderr, timings);
	    fprintf(stderr, "\n");
	    cleanup();
	    return 1;
	}
	if (testCase.startupExpand) {
	    if (testCase.startupExpand) {
		/* This fixture exercises cache-backed source expansion, not mesh
		 * LoD submission.  Bound its worker pool and leave auto-submit off
		 * so its completion condition is solely the requested prewarm. */
		if (!controller->ensureManagedLodService(4)) {
		    fprintf(stderr,
			    "qtcad Obol progressive LoD startup-expand could not start its prewarm service: case=%s\n",
			    testCase.name);
		    cleanup();
		    return 0;
		}
		controller->setLodAutoSubmit(FALSE);
	    }
	    struct qtcad_obol_source_expansion_status expansion_status;
	    struct qtcad_obol_source_prewarm_status prewarm_status;
	    qtcad_obol_lod_service_status service_status;
	    memset(&expansion_status, 0, sizeof(expansion_status));
	    memset(&prewarm_status, 0, sizeof(prewarm_status));
	    memset(&service_status, 0, sizeof(service_status));
	    int source_count_before =
		qtcad_obol_render_database_source_count(gedp, controller);
	    int max_depth_before =
		qtcad_obol_render_database_source_max_depth(gedp, controller,
		    active_draw_target);
	    double prewarm_seconds = 0.0;
	    double expand_seconds = 0.0;
	    size_t total_prewarm_submitted = 0;
	    int expanded = 0;
	    for (int pass = 0; pass < testCase.startupExpandPasses; pass++) {
		struct qtcad_obol_source_prewarm_status pass_prewarm_status;
		struct qtcad_obol_source_expansion_status pass_expansion_status;
		memset(&pass_prewarm_status, 0, sizeof(pass_prewarm_status));
		memset(&pass_expansion_status, 0, sizeof(pass_expansion_status));

		double pass_prewarm_seconds = 0.0;
		phase_start = bu_gettime();
		size_t prewarm_submitted =
		    qtcad_obol_frontier_prewarm(
			gedp, active_draw_target,
			(int)draw_request.style.draw_mode, 1,
			testCase.minVisitedMeshCount, &pass_prewarm_status);
		total_prewarm_submitted += prewarm_submitted;
		if (!qtcad_obol_wait_for_lod_service_idle(gedp,
			ged_view_context_from_bv(view.viewContext()), 10000,
			&service_status)) {
		    record_progressive_lod_phase(pass_prewarm_seconds,
						 phase_start, total_start,
						 testCase,
						 "prewarm_child_aabb_cache");
		    prewarm_seconds += pass_prewarm_seconds;
		    fprintf(stderr,
			    "qtcad Obol progressive LoD startup-expand timed out waiting for child AABB prewarm: case=%s pass=%d submitted=%zu status_submitted=%zu children=%zu considered=%zu pending=%zu in_flight=%zu delayed=%zu cache_writes=%zu",
			    testCase.name, pass + 1, prewarm_submitted,
			    pass_prewarm_status.submitted,
			    pass_prewarm_status.child_count,
			    pass_prewarm_status.considered,
			    service_status.pending_tasks,
			    service_status.in_flight,
			    service_status.delayed_tasks,
			    service_status.queued_cache_writes);
		    print_progressive_lod_timings(stderr, timings);
		    fprintf(stderr, "\n");
		    cleanup();
		    return 0;
		}
		record_progressive_lod_phase(pass_prewarm_seconds,
					     phase_start, total_start,
					     testCase,
					     "prewarm_child_aabb_cache");
		prewarm_seconds += pass_prewarm_seconds;
		qtcad_obol_prewarm_status_accumulate(&prewarm_status,
						     &pass_prewarm_status);
		if (prewarm_submitted == 0 &&
		    pass_prewarm_status.already_cached == 0 &&
		    pass_prewarm_status.comb_sources == 0) {
		    fprintf(stderr,
			    "qtcad Obol progressive LoD startup-expand failed to submit visible child-frontier AABB prewarm: case=%s pass=%d children=%zu considered=%zu combs=%zu already_cached=%zu skipped_non_union=%zu skipped_duplicate=%zu skipped_invalid=%zu",
			    testCase.name, pass + 1,
			    pass_prewarm_status.child_count,
			    pass_prewarm_status.considered,
			    pass_prewarm_status.comb_sources,
			    pass_prewarm_status.already_cached,
			    pass_prewarm_status.skipped_non_union,
			    pass_prewarm_status.skipped_duplicate_instance,
			    pass_prewarm_status.skipped_invalid);
		    print_progressive_lod_timings(stderr, timings);
		    fprintf(stderr, "\n");
		    cleanup();
		    return 0;
		}

		double pass_expand_seconds = 0.0;
		phase_start = bu_gettime();
		int pass_expanded =
		    qtcad_obol_frontier_expand(
			gedp, active_draw_target,
			(int)draw_request.style.draw_mode, 1,
			testCase.minVisitedMeshCount, &pass_expansion_status);
		record_progressive_lod_phase(pass_expand_seconds, phase_start,
					     total_start, testCase,
					     "expand_children");
		expand_seconds += pass_expand_seconds;
		if (pass_expanded)
		    expanded = 1;
		qtcad_obol_expansion_status_accumulate(&expansion_status,
						       &pass_expansion_status);
	    }
	    QImage frame_expand;
	    int source_count_after =
		qtcad_obol_render_database_source_count(gedp, controller);
	    int max_depth_after =
		qtcad_obol_render_database_source_max_depth(gedp, controller,
		    active_draw_target);
	    if (!qtcad_obol_autoview_refresh(gedp, view, controller,
					     "startup-expand-autoview")) {
		fprintf(stderr,
			"qtcad Obol progressive LoD startup-expand failed to autoview expanded proxy bounds: case=%s",
			testCase.name);
		print_progressive_lod_timings(stderr, timings);
		fprintf(stderr, "\n");
		cleanup();
		return 0;
	    }
	    QCoreApplication::processEvents();
	    view.get_viewport_image(frame_expand);
	    int lit_expand = lit_pixel_count(frame_expand);
	    int diff_expand = image_byte_diff(frame0, frame_expand);
	    struct lit_pixel_bounds expand_bounds;
	    int centered_expand =
		lit_pixel_bounds_centered(frame_expand, 0.30, &expand_bounds);
	    if (!expanded || expansion_status.expanded == 0 ||
		source_count_after <= source_count_before ||
		(testCase.startupExpandPasses > 1 &&
		 max_depth_after <= max_depth_before + 1) ||
		expansion_status.proxy_published == 0 ||
		frame_expand.isNull() || lit_expand < 20 ||
		diff_expand <= 0) {
		fprintf(stderr,
			"qtcad Obol progressive LoD startup-expand failed: case=%s passes=%d expanded=%d status_expanded=%zu children=%zu considered=%zu remaining=%zu prewarm_submitted=%zu prewarm_ready=%zu proxies=%zu metadata=%zu sources=%d->%d depth=%d->%d lit=%d diff=%d centered=%d lit_bbox=(%d,%d)-(%d,%d) center=(%.1f,%.1f)",
			testCase.name, testCase.startupExpandPasses,
			expanded, expansion_status.expanded,
			expansion_status.child_count,
			expansion_status.considered,
			expansion_status.remaining,
			prewarm_status.submitted,
			prewarm_status.already_cached,
			expansion_status.proxy_published,
			expansion_status.metadata_applied,
			source_count_before, source_count_after,
			max_depth_before, max_depth_after, lit_expand, diff_expand,
			centered_expand, expand_bounds.minX, expand_bounds.minY,
			expand_bounds.maxX, expand_bounds.maxY,
			expand_bounds.centerX, expand_bounds.centerY);
		print_progressive_lod_timings(stderr, timings);
		fprintf(stderr, "\n");
		cleanup();
		return 0;
	    }
	    fprintf(stderr,
		    "qtcad_progressive_lod_startup_expand case=%s passes=%d expanded=%zu children=%zu considered=%zu existing=%zu remaining=%zu prewarm_submitted=%zu prewarm_ready=%zu total_prewarm_submitted=%zu prewarm_seconds=%.3f proxies=%zu metadata=%zu sources=%d->%d depth=%d->%d lit=%d diff=%d expand_seconds=%.3f",
		    testCase.name, testCase.startupExpandPasses,
		    expansion_status.expanded,
		    expansion_status.child_count, expansion_status.considered,
		    expansion_status.existing, expansion_status.remaining,
		    prewarm_status.submitted, prewarm_status.already_cached,
		    total_prewarm_submitted, prewarm_seconds,
		    expansion_status.proxy_published,
		    expansion_status.metadata_applied,
		    source_count_before, source_count_after,
		    max_depth_before, max_depth_after, lit_expand, diff_expand,
		    expand_seconds);
	    print_progressive_lod_timings(stderr, timings);
	    fprintf(stderr, "\n");
	    if (testCase.startupExpand) {
		cleanup();
		return 1;
	    }
	}

	    struct ged_scene_draw_request refine_request = draw_request;
	    refine_request.realization.mode = GED_SCENE_REALIZE_EAGER;
	    refine_request.autoview = 0;

	    draw_result = ged_scene_result_create();
	    double refine_draw_seconds = 0.0;
	    phase_start = bu_gettime();
	    int refine_draw_ret = ged_scene_draw(gedp, &refine_request,
		draw_result) == GED_SCENE_OK ? 1 : -1;
	    record_progressive_lod_phase(refine_draw_seconds, phase_start,
					 total_start, testCase,
					 "refine_draw_apply");

	    double refine_sync_seconds = 0.0;
	    phase_start = bu_gettime();
	    if (ged_scene_result_changed(draw_result))
		view.need_update(QG_VIEW_REFRESH);
	    record_progressive_lod_phase(refine_sync_seconds, phase_start,
		    total_start, testCase, "refine_obol_sync");
	    int refine_source_count =
		qtcad_obol_render_database_source_count(gedp, controller);
	    int refine_mesh_count = realized_mesh_count(controller);

	    if (refine_draw_ret < 0 || refine_source_count <= 0 ||
		refine_mesh_count <= 0) {
		fprintf(stderr,
			"qtcad progressive LoD startup-refine draw failed: draw_ret=%d source_count=%d mesh_count=%d errors=%s\n",
			refine_draw_ret, refine_source_count, refine_mesh_count,
			ged_scene_result_diagnostic(draw_result));
		cleanup();
		return 0;
	    }

	    ged_scene_result_destroy(draw_result);
	    draw_result = NULL;
	    qtcad_obol_print_render_diagnostics(gedp, testCase, controller,
						"after-refine-draw");
    }

    phase_start = bu_gettime();
    if (!service.start(4, TRUE)) {
	fprintf(stderr, "unable to start qtcad Obol LoD service\n");
	cleanup();
	return 0;
    }
    service_started = 1;
    controller->setLodService(&service);
    record_progressive_lod_phase(timings.serviceStartSeconds, phase_start,
				 total_start, testCase, "service_start");

    int source_mesh_count = realized_mesh_count(controller);
    if (source_mesh_count <= 0) {
	fprintf(stderr, "qtcad Obol progressive LoD could not find realized mesh geometry\n");
	cleanup();
	return 0;
	}

    int64_t submit_start = bu_gettime();
    int submitted = controller->submitLodRequestsIfNeeded();
    record_progressive_lod_phase(timings.submitSeconds, submit_start,
				 total_start, testCase, "submit");
    unsigned int visited_mesh_count = controller->getLastLodVisitedMeshCount();
    unsigned int submitted_task_count =
	controller->getLastLodSubmittedTaskCount();

    size_t max_active_tasks = 0;
    size_t max_queued_results = 0;
    size_t max_queued_cache_writes = 0;
    service.getQueueLimits(max_active_tasks, max_queued_results,
	max_queued_cache_writes);
    const size_t submission_limit = std::min(max_active_tasks,
	max_queued_results);
    const uint64_t rejected_tasks =
	service.rejectedTaskCountForDiagnostics();
    if (submitted <= 0 || visited_mesh_count == 0 ||
	submitted_task_count != static_cast<unsigned int>(submitted) ||
	(submission_limit && submitted_task_count > submission_limit) ||
	(visited_mesh_count * 3 > submitted_task_count &&
	 rejected_tasks == 0)) {
	fprintf(stderr, "qtcad Obol progressive LoD submitted unexpected work: submitted=%d source_meshes=%d visited=%u tasks=%u diagnostics=%s\n",
		submitted, source_mesh_count, visited_mesh_count,
		submitted_task_count,
		controller->getLastLodDiagnostics().getString());
	cleanup();
	return 0;
    }
    if (visited_mesh_count < testCase.minVisitedMeshCount) {
	fprintf(stderr, "qtcad Obol progressive LoD did not exercise multiple LoD candidates: source_meshes=%d visited=%u tasks=%u diagnostics=%s\n",
		source_mesh_count, visited_mesh_count, submitted_task_count,
		controller->getLastLodDiagnostics().getString());
	cleanup();
	return 0;
    }

    size_t applied_total = 0;
    size_t active_aabb_payload_count = 0;
    size_t active_obb_payload_count = 0;
    QImage frame1;
    QImage frame2;
    if (testCase.diagnosticStages) {
	if (!process_until_proxy_kind(controller, service,
			      BOBOL_LOD_PROXY_AABB,
			      testCase.minActivePayloadCount, applied_total,
			      testCase.aabbTimeoutMs, "AABB")) {
	    cleanup();
	    return 0;
	}
	record_progressive_lod_phase(timings.aabbSeconds, submit_start,
				     total_start, testCase, "aabb");
	active_aabb_payload_count =
	    controller->getActiveLodProxyPayloadCount(BOBOL_LOD_PROXY_AABB);
	qtcad_obol_request_view_frame(view, controller,
	    "progressive-lod-aabb-frame");
	view.get_viewport_image(frame1);

	if (!process_until_proxy_kind(controller, service,
			      BOBOL_LOD_PROXY_OBB,
			      testCase.minActivePayloadCount, applied_total,
			      testCase.obbTimeoutMs, "OBB")) {
	    cleanup();
	    return 0;
	}
	record_progressive_lod_phase(timings.obbSeconds, submit_start,
				     total_start, testCase, "obb");
	active_obb_payload_count =
	    controller->getActiveLodProxyPayloadCount(BOBOL_LOD_PROXY_OBB);
	qtcad_obol_request_view_frame(view, controller,
	    "progressive-lod-obb-frame");
	view.get_viewport_image(frame2);
    } else {
	/* Coalesced delivery need not expose transient proxy stages. */
	frame1 = frame0;
	frame2 = frame0;
    }

	if (!process_until_mesh(controller, service,
				testCase.minActivePayloadCount, applied_total,
				testCase.meshTimeoutMs)) {
	cleanup();
	return 0;
    }
    record_progressive_lod_phase(timings.meshSeconds, submit_start,
				 total_start, testCase, "mesh");
    double final_phase_seconds = 0.0;
    phase_start = bu_gettime();
    size_t active_mesh_payload_count =
	controller->getActiveLodMeshPayloadCount();
    record_progressive_lod_phase(final_phase_seconds, phase_start,
				 total_start, testCase,
				 "active_mesh_payloads");
    phase_start = bu_gettime();
    QImage frame3;
    qtcad_obol_request_view_frame(view, controller,
	"progressive-lod-mesh-frame");
    view.get_viewport_image(frame3);
    record_progressive_lod_phase(final_phase_seconds, phase_start,
				 total_start, testCase, "readback3");

    if (phase_logging_enabled(testCase)) {
	fprintf(stderr,
		"qtcad_progressive_lod_service_stop_state case=%s pending=%zu queued_results=%zu queued_cache_writes=%zu in_flight=%zu delayed=%zu active_mesh_payloads=%zu\n",
		testCase.name,
		service.pendingTaskCountForDiagnostics(),
		service.queuedResultCountForDiagnostics(),
		service.queuedCacheWriteCountForDiagnostics(),
		service.inFlightCount(),
		service.delayedTaskCountForDiagnostics(),
		active_mesh_payload_count);
	fflush(stderr);
    }

    phase_start = bu_gettime();
    controller->setLodAutoSubmit(FALSE);
    controller->setLodService(NULL);
    service.stop();
    service_started = 0;
    record_progressive_lod_phase(final_phase_seconds, phase_start,
				 total_start, testCase, "service_stop");

    if (applied_total == 0) {
	fprintf(stderr, "qtcad Obol progressive LoD applied no results\n");
	cleanup();
	return 0;
    }

    phase_start = bu_gettime();
    int lit0 = lit_pixel_count(frame0);
    int lit1 = lit_pixel_count(frame1);
    int lit2 = lit_pixel_count(frame2);
    int lit3 = lit_pixel_count(frame3);
    int diff01 = image_byte_diff(frame0, frame1);
    int diff12 = image_byte_diff(frame1, frame2);
    int diff23 = image_byte_diff(frame2, frame3);
    int diff03 = image_byte_diff(frame0, frame3);
    record_progressive_lod_phase(final_phase_seconds, phase_start,
				 total_start, testCase, "image_metrics");
    const int any_frame_difference =
	diff01 > 0 || diff12 > 0 || diff23 > 0;
    /* Non-diagnostic delivery is allowed to coalesce directly to the final
     * payload.  Requiring a raster difference there turns the intended
     * final-before-first-frame fast path into a timing-dependent failure. */
    const int image_progression_ok = !testCase.diagnosticStages ||
	(diff03 > 0 && any_frame_difference);
    if (frame0.isNull() || frame1.isNull() || frame2.isNull() ||
	frame3.isNull() || (lit0 < 20 && !testCase.startupRefine) ||
	lit1 < 20 || lit2 < 20 ||
	lit3 < 20 || !image_progression_ok) {
	fprintf(stderr, "qtcad Obol progressive LoD image check failed: case=%s lit=%d,%d,%d,%d diff=%d,%d,%d,%d applied=%zu source_meshes=%d visited=%u active_payloads=%zu,%zu,%zu",
		testCase.name,
		lit0, lit1, lit2, lit3, diff01, diff12, diff23, diff03,
		applied_total, source_mesh_count, visited_mesh_count,
		active_aabb_payload_count, active_obb_payload_count,
		active_mesh_payload_count);
	print_progressive_lod_timings(stderr, timings);
	fprintf(stderr, "\n");
	cleanup();
	return 0;
    }

    fprintf(stderr, "qtcad_progressive_lod case=%s diff01=%d diff12=%d diff23=%d diff03=%d lit=%d,%d,%d,%d applied=%zu source_meshes=%d visited=%u active_payloads=%zu,%zu,%zu",
	    testCase.name, diff01, diff12, diff23, diff03, lit0, lit1, lit2,
	    lit3, applied_total, source_mesh_count, visited_mesh_count,
	    active_aabb_payload_count, active_obb_payload_count,
	    active_mesh_payload_count);
    print_progressive_lod_timings(stderr, timings);
    fprintf(stderr, "\n");

    cleanup();
    return 1;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    bu_setenv("QT_QPA_PLATFORM", "offscreen", 1);
    bu_setenv("QTCAD_OBOL_FORCE_OSMESA", "1", 1);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    struct progressive_lod_case testCase;
    testCase.name = "moss";
    testCase.sourceDb = NULL;
    testCase.facetizeRoot = "all.g";
    testCase.drawTarget = "all.bot";
    testCase.useFacetize = 1;
    testCase.minVisitedMeshCount = 2;
    testCase.minActivePayloadCount = 2;
    testCase.aabbTimeoutMs = 1000;
    testCase.obbTimeoutMs = 1500;
    testCase.meshTimeoutMs = 3000;
    testCase.copyCount = 1;
    testCase.startupOnly = 0;
    testCase.startupRefine = 0;
    testCase.startupExpand = 0;
    testCase.startupAutoExpand = 0;
    testCase.startupWireLod = 0;
    testCase.shadedLod = 0;
    testCase.slowLodDelays = 0;
    testCase.diagnosticStages = 0;
    testCase.startupExpandPasses = 1;

    if (argc == 2) {
	testCase.sourceDb = argv[1];
    } else {
	int arg = 1;
	int slow = 0;
	while (arg < argc && argv[arg] && argv[arg][0] == '-') {
	    if (BU_STR_EQUAL(argv[arg], "--slow")) {
		slow = 1;
		arg++;
		continue;
	    }
	    if (BU_STR_EQUAL(argv[arg], "--diagnostic-stages")) {
		testCase.diagnosticStages = 1;
		arg++;
		continue;
	    }
	    if (BU_STR_EQUAL(argv[arg], "--startup-only")) {
		testCase.startupOnly = 1;
		arg++;
		continue;
	    }
	    if (BU_STR_EQUAL(argv[arg], "--startup-refine")) {
		testCase.startupRefine = 1;
		arg++;
		continue;
	    }
	    if (BU_STR_EQUAL(argv[arg], "--startup-expand")) {
		testCase.startupExpand = 1;
		arg++;
		continue;
	    }
	    if (BU_STR_EQUAL(argv[arg], "--startup-auto-expand")) {
		testCase.startupAutoExpand = 1;
		arg++;
		continue;
	    }
	    if (BU_STR_EQUAL(argv[arg], "--startup-wire-lod")) {
		testCase.startupWireLod = 1;
		arg++;
		continue;
	    }
	    if (BU_STR_EQUAL(argv[arg], "--shaded-lod")) {
		testCase.shadedLod = 1;
		arg++;
		continue;
	    }
	    if (BU_STR_EQUAL(argv[arg], "--startup-expand-passes")) {
		if (arg + 1 >= argc)
		    FAIL("--startup-expand-passes requires a positive integer");
		testCase.startupExpandPasses =
		    parse_positive_int(argv[arg + 1], 1);
		arg += 2;
		continue;
	    }
	    break;
	}

	if ((testCase.startupOnly ? 1 : 0) +
	    (testCase.startupRefine ? 1 : 0) +
	    (testCase.startupExpand ? 1 : 0) +
	    (testCase.startupAutoExpand ? 1 : 0) +
	    (testCase.startupWireLod ? 1 : 0) > 1)
	    FAIL("--startup-only, --startup-refine, --startup-expand, --startup-auto-expand, and --startup-wire-lod are mutually exclusive");
	if (testCase.shadedLod &&
	    (testCase.startupOnly || testCase.startupRefine ||
	     testCase.startupExpand || testCase.startupAutoExpand ||
	     testCase.startupWireLod))
	    FAIL("--shaded-lod is a separate direct shaded mode test");

	const char *run_slow_tests = getenv("BRLCAD_RUN_SLOW_TESTS");
	if (slow && (!run_slow_tests ||
		     !BU_STR_EQUAL(run_slow_tests, "1"))) {
	    fprintf(stderr, "SKIP: set BRLCAD_RUN_SLOW_TESTS=1 to run slow qtcad Obol progressive LoD case\n");
	    return 123;
	}

	if (argc - arg < 4 || argc - arg > 7)
	    FAIL("usage: test_qtcad_obol_progressive_lod <moss.g> OR test_qtcad_obol_progressive_lod [--slow] [--diagnostic-stages] [--startup-only|--startup-refine|--startup-expand|--startup-auto-expand|--startup-wire-lod] [--shaded-lod] [--startup-expand-passes N] <db> <facetize-root|-> <draw-target> <facetize|direct> [min-visited] [min-active-payloads] [copy-count]");

	if (testCase.startupOnly)
	    testCase.name = "startup";
	else if (testCase.startupRefine)
	    testCase.name = "startup_refine";
	else if (testCase.startupExpand)
	    testCase.name = "startup_expand";
	else if (testCase.startupAutoExpand)
	    testCase.name = "startup_auto_expand";
	else if (testCase.startupWireLod)
	    testCase.name = "startup_wire_lod";
	else if (testCase.shadedLod)
	    testCase.name = "shaded_lod";
	else
	    testCase.name = slow ? "slow" : "custom";
	testCase.sourceDb = argv[arg++];
	testCase.facetizeRoot = argv[arg++];
	testCase.drawTarget = argv[arg++];
	const char *mode = argv[arg++];
	if (BU_STR_EQUAL(mode, "facetize"))
	    testCase.useFacetize = 1;
	else if (BU_STR_EQUAL(mode, "direct"))
	    testCase.useFacetize = 0;
	else
	    FAIL("progressive LoD mode must be 'facetize' or 'direct'");

	if (testCase.useFacetize &&
	    (!testCase.facetizeRoot || BU_STR_EQUAL(testCase.facetizeRoot, "-")))
	    FAIL("facetize mode requires a facetize root");

	testCase.minVisitedMeshCount = (unsigned int)parse_positive_int(
					   arg < argc ? argv[arg++] : NULL, slow ? 10 : 2);
	testCase.minActivePayloadCount = parse_positive_size(
					     arg < argc ? argv[arg++] : NULL, 2);
	testCase.copyCount = parse_positive_int(
				 arg < argc ? argv[arg++] : NULL, 1);
	if (slow) {
	    testCase.slowLodDelays = 1;
	    testCase.diagnosticStages = 1;
	    testCase.aabbTimeoutMs = 30000;
	    testCase.obbTimeoutMs = 45000;
	    testCase.meshTimeoutMs = 90000;
	}
    }

    QApplication app(argc, argv);

    if (!run_progressive_lod_case(testCase))
	FAIL("qtcad Obol progressive LoD workflow should pass");

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
