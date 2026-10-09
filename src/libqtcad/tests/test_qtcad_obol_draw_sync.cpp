/*          T E S T _ Q T C A D _ O B O L _ D R A W _ S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bv.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BLodService.h"
#include "BObol/BLodRealization.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "BObol/BSceneController.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/malloc.h"
#include "bu/process.h"
#include "bu/str.h"
#include "ged.h"
#include "ged/scene.h"
#include "ged/display_obol_private.h"
#include "ged/draw.h"
#include "ged/event.h"
#include "ged/view.h"
#include "QgSceneSyncPrivate.h"
#include "QgCanvasState.h"
#include "qtcad/QgView.h"
#include "qtcad_obol_test_presentation.h"
#include "raytrace.h"
#include "rt/db_fullpath.h"
#include "wdb.h"

#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoGroup.h>

#include <QApplication>
#include <QCoreApplication>
#include <QImage>
#include <QPainter>

#include <algorithm>
#include <chrono>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include <string>
#include <thread>
#include <vector>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

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

static int
nonblack_pixel_count(const QImage &image)
{
    QImage rgba = image.convertToFormat(QImage::Format_RGBA8888);
    int count = 0;
    for (int y = 0; y < rgba.height(); y++) {
	const unsigned char *line = rgba.constScanLine(y);
	for (int x = 0; x < rgba.width(); x++) {
	    const unsigned char *p = line + x * 4;
	    if (p[0] > 2 || p[1] > 2 || p[2] > 2)
		count++;
	}
    }
    return count;
}

static int
image_byte_diff(const QImage &a, const QImage &b)
{
    QImage ar = a.convertToFormat(QImage::Format_RGBA8888);
    QImage br = b.convertToFormat(QImage::Format_RGBA8888);
    if (ar.size() != br.size())
	return -1;

    int changed = 0;
    for (int y = 0; y < ar.height(); y++) {
	const unsigned char *ap = ar.constScanLine(y);
	const unsigned char *bp = br.constScanLine(y);
	for (int x = 0; x < ar.bytesPerLine(); x++)
	    changed += ap[x] != bp[x] ? 1 : 0;
    }
    return changed;
}

static int
camera_positions_differ(const SbVec3f &a, const SbVec3f &b, float min_delta)
{
    return fabsf(a[0] - b[0]) > min_delta ||
	fabsf(a[1] - b[1]) > min_delta ||
	fabsf(a[2] - b[2]) > min_delta;
}

/* The direct camera is placed one focal distance behind the BRL-CAD view
 * center.  Position therefore includes orientation-dependent depth and is
 * not a valid assertion for a recenter operation. */
static SbVec3f
camera_focal_point(const SoCamera *camera)
{
    if (!camera)
	return SbVec3f(0.0f, 0.0f, 0.0f);
    SbVec3f forward(0.0f, 0.0f, -1.0f);
    camera->orientation.getValue().multVec(forward, forward);
    return camera->position.getValue() +
	forward * camera->focalDistance.getValue();
}

static int
make_draw_sync_db(const char *dbpath)
{
    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp)
	return 0;

    point_t bmin = {-1.0, -1.0, -1.0};
    point_t bmax = { 1.0,  1.0,  1.0};
    point_t center = {5.0, 0.0, 0.0};

    static constexpr int mesh_grid_size = 12;
    static constexpr fastf_t mesh_height_step = 0.05;
    const int mesh_vertex_count =
	(mesh_grid_size + 1) * (mesh_grid_size + 1);
    const int mesh_face_count = mesh_grid_size * mesh_grid_size * 2;
    std::vector<fastf_t> mesh_vertices(mesh_vertex_count * 3, 0.0);
    std::vector<int> mesh_faces(mesh_face_count * 3, 0);
    for (int y = 0; y <= mesh_grid_size; y++) {
	for (int x = 0; x <= mesh_grid_size; x++) {
	    const int vertex = y * (mesh_grid_size + 1) + x;
	    mesh_vertices[3 * vertex + X] = static_cast<fastf_t>(x);
	    mesh_vertices[3 * vertex + Y] = static_cast<fastf_t>(y);
	    mesh_vertices[3 * vertex + Z] =
		static_cast<fastf_t>((x + y) % 3) * mesh_height_step;
	}
    }
    for (int y = 0; y < mesh_grid_size; y++) {
	for (int x = 0; x < mesh_grid_size; x++) {
	    const int cell = y * mesh_grid_size + x;
	    const int v0 = y * (mesh_grid_size + 1) + x;
	    const int v1 = v0 + 1;
	    const int v2 = v0 + mesh_grid_size + 1;
	    const int v3 = v2 + 1;
	    mesh_faces[6 * cell + 0] = v0;
	    mesh_faces[6 * cell + 1] = v1;
	    mesh_faces[6 * cell + 2] = v3;
	    mesh_faces[6 * cell + 3] = v0;
	    mesh_faces[6 * cell + 4] = v3;
	    mesh_faces[6 * cell + 5] = v2;
	}
    }

    int ret = mk_rpp(wdbp, "box.s", bmin, bmax) == 0 &&
	mk_sph(wdbp, "ball.s", center, 1.0) == 0 &&
	mk_bot(wdbp, "flow.bot", RT_BOT_SURFACE, RT_BOT_CCW, 0,
	    mesh_vertex_count, mesh_face_count, mesh_vertices.data(),
	    mesh_faces.data(), NULL, NULL) == 0;
    struct wmember pair;
    BU_LIST_INIT(&pair.l);
    ret = ret &&
	mk_addmember("box.s", &pair.l, NULL, WMOP_UNION) != NULL &&
	mk_addmember("ball.s", &pair.l, NULL, WMOP_UNION) != NULL &&
	mk_lcomb(wdbp, "pair.c", &pair, 0, NULL, NULL, NULL, 0) == 0;
    wdb_close(wdbp);
    return ret;
}

static const char *
test_skip_leading_slash(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

static int
test_path_equal(const char *a, const char *b)
{
    if (!a || !b)
	return 0;
    if (BU_STR_EQUAL(a, b))
	return 1;
    return BU_STR_EQUAL(test_skip_leading_slash(a),
	    test_skip_leading_slash(b));
}

static int
test_shape_record_matches_path(const struct ged_scene_occurrence_info *record,
	const char *path)
{
    if (!record || !path || !path[0])
	return 0;
    if (test_path_equal(record->path, path) ||
	    test_path_equal(record->leaf_name, path))
	return 1;

    if (!record->fullpath || record->fullpath->fp_len <= 0)
	return 0;

    char *recordPath = db_path_to_string(record->fullpath);
    if (!recordPath)
	return 0;
    int matched = test_path_equal(recordPath, path);
    bu_free(recordPath, "qtcad Obol draw-sync test record path");
    return matched;
}

struct test_shape_record_path_context {
    const char *path;
    struct ged_scene_occurrence_info *out;
    int found;
};

static int
test_shape_record_by_path_cb(const struct ged_scene_occurrence_info *record,
	void *userdata)
{
    struct test_shape_record_path_context *ctx =
	(struct test_shape_record_path_context *)userdata;
    if (!ctx || !record || !test_shape_record_matches_path(record,
	    ctx->path))
	return 1;

    if (ctx->out)
	*ctx->out = *record;
    ctx->found = 1;
    return 0;
}

static int
test_shape_record_by_path(struct ged *gedp,
	const char *path,
	struct ged_scene_occurrence_info *out)
{
    if (!gedp || !path || !path[0])
	return 0;

    struct test_shape_record_path_context ctx;
    ctx.path = path;
    ctx.out = out;
    ctx.found = 0;
    ged_scene_occurrences_visit(gedp, test_shape_record_by_path_cb, &ctx);
    return ctx.found;
}

static SoBRLDatabaseSource *
render_source(BObolViewController *controller, int index);

static void
collect_render_sources(SoNode *node,
	std::vector<SoBRLDatabaseSource *> &sources)
{
    if (!node)
	return;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoBRLDatabaseSource *source =
	    static_cast<SoBRLDatabaseSource *>(node);
	if (std::find(sources.begin(), sources.end(), source) == sources.end())
	    sources.push_back(source);
	return;
    }
    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;
    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++)
	collect_render_sources(group->getChild(i), sources);
}

static std::vector<SoBRLDatabaseSource *>
render_sources(BObolViewController *controller)
{
    std::vector<SoBRLDatabaseSource *> sources;
    if (controller)
	collect_render_sources(controller->getRenderSceneRoot(), sources);
    if (sources.empty() && controller)
	collect_render_sources(controller->getSceneRoot(), sources);
    return sources;
}

static int
render_source_count(BObolViewController *controller)
{
    return static_cast<int>(render_sources(controller).size());
}

static int
visible_compact_occurrence_count(SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;
    int count = 0;
    for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	BObolCompactOccurrence occurrence;
	if (source->getCompactOccurrence(i, occurrence) &&
		occurrence.summary.visible)
	    count++;
    }
    return count;
}

static int
compact_occurrence_for_path(SoBRLDatabaseSource *source, const char *path,
	BObolCompactOccurrence *out)
{
    if (!source || !path)
	return 0;
    for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	BObolCompactOccurrence occurrence;
	if (!source->getCompactOccurrence(i, occurrence) ||
	    !test_path_equal(occurrence.summary.path.getString(), path))
	    continue;
	if (out)
	    *out = occurrence;
	return 1;
    }
    return 0;
}

static SoBRLDatabaseSource *
render_source(BObolViewController *controller, int index)
{
    std::vector<SoBRLDatabaseSource *> sources = render_sources(controller);
    return index >= 0 && static_cast<size_t>(index) < sources.size() ?
	sources[static_cast<size_t>(index)] : NULL;
}

static SoBRLDatabaseSource *
source_for_path(BObolViewController *controller, const char *path)
{
    if (!controller || !path)
	return NULL;
    for (SoBRLDatabaseSource *source : render_sources(controller)) {
	if (source && test_path_equal(source->path.getValue().getString(),
		path))
	    return source;
    }
    return NULL;
}

static SoBRLDatabaseSource *
source_for_path_mode(BObolViewController *controller,
	const char *path,
	int draw_mode)
{
    if (!controller || !path)
	return NULL;
    for (SoBRLDatabaseSource *source : render_sources(controller)) {
	if (source && test_path_equal(source->path.getValue().getString(),
		path) && source->drawMode.getValue() == draw_mode)
	    return source;
    }
    return NULL;
}

static int
scene_display_summary_by_path(BObolViewController *controller,
	const char *path,
	int nodeKind,
	BObolSceneDisplaySummary *out)
{
    if (!controller || !controller->getSceneController() || !path)
	return 0;

    const BObolSceneController *scene = controller->getSceneController();
    for (int i = 0; i < scene->getSceneDisplaySummaryCount(); i++) {
	BObolSceneDisplaySummary summary;
	if (!scene->getSceneDisplaySummary(i, summary) || !summary.valid)
	    continue;
	if (summary.nodeKind == nodeKind &&
		test_path_equal(summary.path.getString(), path)) {
	    if (out)
		*out = summary;
	    return 1;
	}
    }
    return 0;
}

struct draw_observer_sync_context {
    QgView *view;
    int calls;
    int changed;
};

static void
test_draw_observer(struct ged *gedp,
	const struct ged_scene_delta *delta,
	void *client_data)
{
    struct draw_observer_sync_context *ctx =
	static_cast<struct draw_observer_sync_context *>(client_data);
    if (!ctx || !ctx->view)
	return;
    ctx->calls++;
    ctx->changed += qg_scene_delta_notify(gedp, delta, ctx->view);
}

static int
scene_result_matches(enum ged_scene_status status,
	struct ged_scene_result *result,
	int expect_changed)
{
    const int changed = ged_scene_result_changed(result);
    const int matches = status == GED_SCENE_OK &&
	(expect_changed ? changed != 0 : changed == 0);
    ged_scene_result_destroy(result);
    return matches;
}

static int
pending_hud_provider(BObolViewController *, void *,
    const BObolProgressiveOptions *, BObolProgressiveStatus *status)
{
    status->hasMore = 1;
    status->sourcePreparationTotalUnits = 1;
    return 0;
}

static BObolLodConvergenceStatus
active_lod_overlay_status(int phase)
{
    BObolLodConvergenceStatus status;
    status.hasLodState = TRUE;
    status.phase = phase;
    status.outcome = BOBOL_LOD_PRESENTATION_ACTIVE;
    status.terminal = FALSE;
    status.viewReady = FALSE;
    status.fraction = 0.25f;
    status.episodeRevision = 1;
    status.viewRevision = 1;
    status.policyRevision = 1;
    return status;
}

static int
check_native_lod_overlay_state(void)
{
    BObolLodConvergenceStatus status =
	active_lod_overlay_status(BOBOL_LOD_CONVERGENCE_PREPARING);
    status.proxyReasons.sourcePreparationOccurrenceCount = 1250;
    status.proxyReasons.reasonMask =
	BOBOL_LOD_PROXY_REASON_SOURCE_PREPARATION;
    status.presentedStructuralBoxCount = 1250;
    status.activeFaces = 96;
    status.temporaryCoverageOccurrenceCount = 2;
    status.sourcePreparationPending = TRUE;
    status.sourcePreparationCompletedUnits = 1079;
    status.sourcePreparationTotalUnits = 1082;
    status.episode.elapsedMilliseconds = 1500;
    status.progressEstimateAvailable = TRUE;
    status.estimatedFraction = 0.42f;
    QgLodProgressOverlayState overlay =
	qgcanvas_lod_progress_overlay_state(status);
    if (!overlay.visible || !overlay.determinate || overlay.percent != 42 ||
	overlay.title !=
	    QStringLiteral("Preparing geometry") ||
	!overlay.detail.contains(QStringLiteral("1.3k temporary boxes")) ||
	!overlay.detail.contains(
	    QStringLiteral("3 source preparations remaining")) ||
	!overlay.detail.contains(
	    QStringLiteral("2 temporary coverage previews")) ||
	!overlay.detail.contains(QStringLiteral("first mesh pending")) ||
	!overlay.detail.contains(QStringLiteral("1.5 s elapsed")))
	FAIL("native LoD overlay did not describe cold source preparation");

    status = active_lod_overlay_status(BOBOL_LOD_CONVERGENCE_DISCOVERING);
    status.proxyReasons.visibilityPlanningOccurrenceCount = 4;
    status.proxyReasons.reasonMask =
	BOBOL_LOD_PROXY_REASON_VISIBILITY_PLANNING;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (!overlay.visible ||
	overlay.title != QStringLiteral("Preparing geometry"))
	FAIL("native LoD overlay did not describe visibility planning");

    status = active_lod_overlay_status(BOBOL_LOD_CONVERGENCE_REFINING);
    status.proxyReasons.geometryPreparationOccurrenceCount = 3;
    status.proxyReasons.reasonMask =
	BOBOL_LOD_PROXY_REASON_GEOMETRY_PREPARATION;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (!overlay.visible ||
	overlay.title != QStringLiteral("Preparing geometry"))
	FAIL("native LoD overlay did not describe geometry preparation");

    status.producerStageMask =
	(1u << (BOBOL_LOD_PRODUCER_STAGE_SOURCE_HASHING - 1)) |
	(1u << (BOBOL_LOD_PRODUCER_STAGE_FACE_CLASSIFICATION - 1));
    status.producerStage = BOBOL_LOD_PRODUCER_STAGE_SOURCE_HASHING;
    status.producerStageTaskCount = 2;
    status.activeProducerCount = 3;
    status.producerStageTaskCounts[
	BOBOL_LOD_PRODUCER_STAGE_SOURCE_HASHING] = 2;
    status.producerStageTaskCounts[
	BOBOL_LOD_PRODUCER_STAGE_FACE_CLASSIFICATION] = 1;
    status.producerStageCompletedUnits = 25;
    status.producerStageTotalUnits = 100;
    status.oldestPendingTaskAgeMicroseconds = 1500000;
    status.maximumProducerQueueWaitMicroseconds = 250000;
    status.maximumProducerElapsedMicroseconds = 2400000;
    status.producerStageElapsedMicroseconds = 1100000;
    status.activeProducerSourceFaceCount = 1200000;
    status.activeProducerSourceByteCount = 32 * 1024 * 1024;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (!overlay.visible ||
	overlay.title != QStringLiteral("Preparing geometry") ||
	!overlay.detail.contains(QStringLiteral("oldest queued 1.5 s")) ||
	!overlay.detail.contains(QStringLiteral("current stage 1.1 s")) ||
	!overlay.detail.contains(QStringLiteral("active 2.4 s total")) ||
	!overlay.detail.contains(QStringLiteral("queued 250 ms before start")) ||
	!overlay.detail.contains(QStringLiteral("source 1.2M faces, 32.0 MiB")) ||
	!overlay.detail.contains(QStringLiteral("2 hashing, 1 classifying")) ||
	!overlay.detail.contains(QStringLiteral("hashing 25% (25/100)")))
	FAIL("native LoD overlay did not expose producer latency and source size");

    status.producerStageMask |=
	1u << (BOBOL_LOD_PRODUCER_STAGE_COVERAGE_PREVIEW - 1);
    status.producerStage = BOBOL_LOD_PRODUCER_STAGE_COVERAGE_PREVIEW;
    status.producerStageTaskCount = 1;
    status.producerStageTaskCounts[
	BOBOL_LOD_PRODUCER_STAGE_COVERAGE_PREVIEW] = 1;
    status.producerStageCompletedUnits = 750;
    status.producerStageTotalUnits = 1000;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (!overlay.visible ||
	overlay.title != QStringLiteral("Preparing geometry") ||
	!overlay.detail.contains(QStringLiteral("1 sampling coverage")) ||
	!overlay.detail.contains(
	    QStringLiteral("sampling coverage 75% (750/1.0k)")))
	FAIL("native LoD overlay did not expose cold coverage sampling");

    status.producerStageMask |=
	1u << (BOBOL_LOD_PRODUCER_STAGE_ASSET_SERIALIZATION - 1);
    status.producerStage = BOBOL_LOD_PRODUCER_STAGE_ASSET_SERIALIZATION;
    status.producerStageTaskCount = 1;
    status.producerStageTaskCounts[
	BOBOL_LOD_PRODUCER_STAGE_ASSET_SERIALIZATION] = 1;
    status.producerStageCompletedUnits = 0;
    status.producerStageTotalUnits = 0;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (!overlay.visible ||
	overlay.title != QStringLiteral("Preparing geometry") ||
	!overlay.detail.contains(QStringLiteral("1 waiting for shared asset")))
	FAIL("native LoD overlay hid resident-asset serialization");

    status = active_lod_overlay_status(BOBOL_LOD_CONVERGENCE_PREPARING);
    status.pendingTasks = 2;
    status.oldestPendingTaskAgeMicroseconds = 1750000;
    status.proxyReasons.sourcePreparationOccurrenceCount = 2;
    status.proxyReasons.reasonMask =
	BOBOL_LOD_PROXY_REASON_SOURCE_PREPARATION;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (!overlay.visible ||
	overlay.title != QStringLiteral("Preparing geometry") ||
	!overlay.detail.contains(QStringLiteral("oldest queued 1.8 s")))
	FAIL("native LoD overlay hid queued cold-start work");

    status.runnableQueuedTasks = 1;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (overlay.title != QStringLiteral("Preparing geometry") ||
	!overlay.detail.contains(QStringLiteral("1 ready for worker")))
	FAIL("native LoD overlay hid runnable queued work");

    status.pendingTasks = 0;
    status.runnableQueuedTasks = 0;
    status.cpuAdmissionWaitingTasks = 1;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (overlay.title != QStringLiteral("Preparing geometry") ||
	!overlay.detail.contains(QStringLiteral("1 waiting for CPU")))
	FAIL("native LoD overlay hid process CPU admission wait");

    status.pendingTasks = 0;
    status.cpuAdmissionWaitingTasks = 0;
    status.dependencyBlockedTasks = 0;
    status.transientMemoryBlockedTasks = 2;
    status.taskSubmissionCapacityBlocked = TRUE;
    status.resultSubmissionCapacityBlocked = TRUE;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (overlay.title != QStringLiteral("Preparing geometry") ||
	!overlay.detail.contains(QStringLiteral(
	    "2 waiting on transient memory")) ||
	!overlay.detail.contains(QStringLiteral("producer task slots full")) ||
	!overlay.detail.contains(QStringLiteral("result queue slots full")))
	FAIL("native LoD overlay hid memory/capacity-blocked work");

    status.pendingTasks = 2;
    status.transientMemoryBlockedTasks = 0;
    status.taskSubmissionCapacityBlocked = FALSE;
    status.resultSubmissionCapacityBlocked = FALSE;
    status.dependencyBlockedTasks = 2;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (overlay.title != QStringLiteral("Preparing geometry") ||
	!overlay.detail.contains(QStringLiteral("2 waiting on prerequisites")))
	FAIL("native LoD overlay hid dependency-blocked work");

    status = active_lod_overlay_status(BOBOL_LOD_CONVERGENCE_PREPARING);
    status.proxyReasons = BObolLodProxyReasonStatus();
    status.proxyReasons.rendererPreparationOccurrenceCount = 2;
    status.proxyReasons.reasonMask =
	BOBOL_LOD_PROXY_REASON_RENDERER_PREPARATION;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (!overlay.visible ||
	overlay.title != QStringLiteral("Preparing geometry"))
	FAIL("native LoD overlay did not describe renderer preparation");

    status = active_lod_overlay_status(BOBOL_LOD_CONVERGENCE_PREPARING);
    status.presentedStructuralBoxCount = 12;
    status.proxyReasons.sourcePreparationOccurrenceCount = 6;
    status.proxyReasons.frameBudgetOccurrenceCount = 5;
    status.proxyReasons.terminalFailureOccurrenceCount = 1;
    status.proxyReasons.reasonMask =
	BOBOL_LOD_PROXY_REASON_SOURCE_PREPARATION |
	BOBOL_LOD_PROXY_REASON_FRAME_BUDGET |
	BOBOL_LOD_PROXY_REASON_TERMINAL_FAILURE;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (!overlay.detail.contains(QStringLiteral("6 temporary boxes")) ||
	!overlay.detail.contains(QStringLiteral("5 budget-limited boxes")) ||
	!overlay.detail.contains(QStringLiteral("1 failed box")) ||
	overlay.detail.contains(QStringLiteral("12 box proxies")))
	FAIL("native LoD overlay conflated temporary and terminal boxes");

    status = active_lod_overlay_status(BOBOL_LOD_CONVERGENCE_REFINING);
    status.episode.firstMeshReached = TRUE;
    status.activeSourceFaces = 48;
    status.sourceMeshOccurrenceCount = 2;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (!overlay.visible ||
	overlay.title != QStringLiteral("Refining visible detail") ||
	!overlay.detail.contains(QStringLiteral("48 source triangles")) ||
	overlay.detail.contains(QStringLiteral("first mesh pending")))
	FAIL("native LoD overlay hid post-mesh settling");

    status = BObolLodConvergenceStatus();
    status.hasLodState = TRUE;
    status.phase = BOBOL_LOD_CONVERGENCE_IDLE;
    status.outcome = BOBOL_LOD_PRESENTATION_CONSTRAINED;
    status.terminal = TRUE;
    status.viewReady = TRUE;
    status.fraction = 1.0f;
    status.progressEstimateAvailable = TRUE;
    status.estimatedFraction = 1.0f;
    status.memoryLimited = TRUE;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (!overlay.visible || overlay.percent != 100 ||
	overlay.title != QStringLiteral("Detail limited by memory") ||
	!overlay.detail.contains(
	    QStringLiteral("best available under current budget")))
	FAIL("native LoD overlay hid a terminal memory limit");

    status.memoryLimited = FALSE;
    status.performanceLimited = TRUE;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (!overlay.visible ||
	overlay.title != QStringLiteral("Detail limited by frame budget") ||
	qgcanvas_lod_overlay_color(overlay) != QColor(112, 235, 135))
	FAIL("native LoD overlay did not mark a terminal frame limit complete");

    status.performanceLimited = FALSE;
    status.gpuMemoryPressure = TRUE;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (!overlay.visible ||
	overlay.title != QStringLiteral("Detail limited by memory"))
	FAIL("native LoD overlay hid terminal GPU memory pressure");

    status.gpuMemoryPressure = FALSE;
    status.viewReady = FALSE;
    status.terminalError = TRUE;
    status.outcome = BOBOL_LOD_PRESENTATION_ERROR;
    status.phase = BOBOL_LOD_CONVERGENCE_ERROR;
    status.failedSourceCount = 1;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (!overlay.visible ||
	overlay.title != QStringLiteral("Geometry preparation failed") ||
	!overlay.detail.contains(QStringLiteral("1 failed source")))
	FAIL("native LoD overlay hid a terminal geometry error");

    status = BObolLodConvergenceStatus();
    status.hasLodState = TRUE;
    status.phase = BOBOL_LOD_CONVERGENCE_IDLE;
    status.terminal = TRUE;
    status.viewReady = TRUE;
    status.fraction = 1.0f;
    overlay = qgcanvas_lod_progress_overlay_state(status);
    if (overlay.visible)
	FAIL("native LoD overlay did not clear after unconstrained readiness");

    return 0;
}

static int
check_native_lod_overlay_stability(void)
{
    QgLodProgressPresentationState presentation;
    BObolLodConvergenceStatus status =
	active_lod_overlay_status(BOBOL_LOD_CONVERGENCE_REFINING);
    status.episodeRevision = 7;
    status.viewRevision = 11;
    status.policyRevision = 13;
    status.episode.firstMeshReached = TRUE;
    status.progressEstimateAvailable = TRUE;
    status.estimatedFraction = 0.2f;
    status.estimatedRemainingMilliseconds = 9000;
    status.episode.elapsedMilliseconds = 100;

    QgLodProgressOverlayState overlay =
	qgcanvas_stabilize_lod_progress_overlay(presentation, status,
	    qgcanvas_lod_progress_overlay_state(status));
    if (!overlay.determinate || overlay.percent != 20 || overlay.etaVisible ||
	overlay.title != QStringLiteral("Refining visible detail") ||
	!overlay.detail.startsWith(
	    QStringLiteral("20% | estimating remaining time")))
	FAIL("native LoD overlay did not begin a stable progress episode");

    status.estimatedFraction = 0.35f;
    status.estimatedRemainingMilliseconds = 8200;
    status.episode.elapsedMilliseconds = 900;
    overlay = qgcanvas_stabilize_lod_progress_overlay(presentation, status,
	qgcanvas_lod_progress_overlay_state(status));
    status.estimatedFraction = 0.5f;
    status.estimatedRemainingMilliseconds = 7400;
    status.progressEstimateRefinementCycleBased = TRUE;
    status.estimatedRemainingRefinementCycles = 11;
    status.episode.elapsedMilliseconds = 1700;
    overlay = qgcanvas_stabilize_lod_progress_overlay(presentation, status,
	qgcanvas_lod_progress_overlay_state(status));
    if (!overlay.etaVisible || overlay.percent != 50 ||
	!overlay.detail.contains(
	    QStringLiteral("~11 refinement cycles, about 10 s remaining")))
	FAIL("native LoD overlay did not qualify a stable ETA");

    /* Internal admission-policy epochs are implementation details when the
     * controller keeps the user-visible episode identity. */
    status.policyRevision++;
    status.estimatedFraction = 0.3f;
    status.estimatedRemainingMilliseconds = 8000;
    status.episode.elapsedMilliseconds = 1750;
    overlay = qgcanvas_stabilize_lod_progress_overlay(presentation, status,
	qgcanvas_lod_progress_overlay_state(status));
    if (!overlay.determinate || overlay.percent != 50)
	FAIL("raw policy revision restarted the native progress episode");

    /* A lost cycle forecast is meaningful: frontier discovery or a flat frame
     * must return the card to indeterminate activity immediately. */
    status.progressEstimateAvailable = FALSE;
    status.progressEstimateRefinementCycleBased = FALSE;
    status.estimatedRemainingRefinementCycles = 0;
    status.estimatedFraction = 0.0f;
    status.estimatedRemainingMilliseconds = 0;
    status.fraction = 0.1f;
    status.activePayloadCount = 100;
    status.satisfiedPayloadCount = 75;
    status.retainedRenderCost = 950;
    status.renderCostBudget = 1000;
    status.episode.elapsedMilliseconds = 1800;
    overlay = qgcanvas_stabilize_lod_progress_overlay(presentation, status,
	qgcanvas_lod_progress_overlay_state(status));
    if (overlay.determinate || overlay.etaVisible ||
	!overlay.detail.startsWith(
	    QStringLiteral("25 visible items still refining")) ||
	!overlay.detail.contains(
	    QStringLiteral("render budget 95% allocated")))
	FAIL("native LoD overlay hid exact work after losing its forecast");

    /* A sustained unknown phase remains indeterminate. */
    status.episode.elapsedMilliseconds = 2600;
    overlay = qgcanvas_stabilize_lod_progress_overlay(presentation, status,
	qgcanvas_lod_progress_overlay_state(status));
    if (overlay.determinate || overlay.etaVisible ||
	!overlay.detail.startsWith(
	    QStringLiteral("25 visible items still refining")))
	FAIL("native LoD overlay did not retain exact work without an ETA");

    status.progressEstimateAvailable = TRUE;
    status.estimatedFraction = 0.3f;
    status.estimatedRemainingMilliseconds = 8000;
    status.episode.elapsedMilliseconds = 2700;
    overlay = qgcanvas_stabilize_lod_progress_overlay(presentation, status,
	qgcanvas_lod_progress_overlay_state(status));
    if (!overlay.determinate || overlay.percent != 30 || overlay.etaVisible)
	FAIL("native LoD overlay did not restart a newly qualified forecast");

    status.episodeRevision++;
    status.viewRevision++;
    status.episode.elapsedMilliseconds = 0;
    status.estimatedFraction = 0.15f;
    status.estimatedRemainingMilliseconds = 12000;
    overlay = qgcanvas_stabilize_lod_progress_overlay(presentation, status,
	qgcanvas_lod_progress_overlay_state(status));
    if (overlay.percent != 15 || overlay.etaVisible)
	FAIL("native LoD overlay did not reset for a new view episode");

    /* Reaching a visually stable frame is itself enough to enter the
     * reserved finalization band, even if the remaining certificate work is
     * temporarily unranked. */
    status.episode.stableViewReached = TRUE;
    status.episode.elapsedMilliseconds = 100;
    status.progressEstimateAvailable = FALSE;
    status.fraction = 0.6f;
    overlay = qgcanvas_stabilize_lod_progress_overlay(presentation, status,
	qgcanvas_lod_progress_overlay_state(status));
    if (!overlay.determinate || overlay.percent != 95 || overlay.etaVisible ||
	overlay.title != QStringLiteral("Finalizing view") ||
	!overlay.detail.contains(QStringLiteral(
	    "finalizing | 25 visible items still refining")))
	FAIL("native LoD overlay did not reserve its finalization band");

    status.episode.elapsedMilliseconds = 200;
    status.progressEstimateAvailable = TRUE;
	status.progressEstimateRefinementCycleBased = TRUE;
	status.estimatedRemainingRefinementCycles = 12;
	status.estimatedFraction = 0.96f;
	status.estimatedRemainingMilliseconds = 5000;
    overlay = qgcanvas_stabilize_lod_progress_overlay(presentation, status,
	qgcanvas_lod_progress_overlay_state(status));
    if (overlay.title != QStringLiteral("Finalizing view") ||
	overlay.percent != 96 || overlay.etaVisible ||
	!overlay.detail.contains(
	    QStringLiteral("finalizing | estimating remaining time")))
	FAIL("native LoD overlay did not begin qualifying its final-tail ETA");

    status.episode.elapsedMilliseconds = 1000;
    status.estimatedRemainingRefinementCycles = 10;
    status.estimatedFraction = 0.97f;
    status.estimatedRemainingMilliseconds = 4200;
    overlay = qgcanvas_stabilize_lod_progress_overlay(presentation, status,
	qgcanvas_lod_progress_overlay_state(status));
    if (overlay.percent != 97 || overlay.etaVisible)
	FAIL("native LoD overlay exposed an unqualified final-tail ETA");

    status.episode.elapsedMilliseconds = 1800;
    status.estimatedRemainingRefinementCycles = 8;
    status.estimatedFraction = 0.98f;
    status.estimatedRemainingMilliseconds = 3400;
    overlay = qgcanvas_stabilize_lod_progress_overlay(presentation, status,
	qgcanvas_lod_progress_overlay_state(status));
    if (overlay.percent != 98 || !overlay.etaVisible ||
	!overlay.detail.contains(QStringLiteral(
	    "finalizing | ~8 refinement cycles, under 5 s remaining")))
	FAIL("native LoD overlay did not publish a qualified final-tail ETA");

    return 0;
}

static int
check_native_lod_overlay_paint(void)
{
    QWidget widget;
    widget.resize(360, 120);

    QImage baseline(widget.size(), QImage::Format_ARGB32_Premultiplied);
    baseline.fill(Qt::black);
    QImage painted = baseline;

    BObolLodConvergenceStatus status =
	active_lod_overlay_status(BOBOL_LOD_CONVERGENCE_PREPARING);
    status.proxyReasons.geometryPreparationOccurrenceCount = 24;
    status.proxyReasons.reasonMask =
	BOBOL_LOD_PROXY_REASON_GEOMETRY_PREPARATION;
    status.presentedStructuralBoxCount = 24;
    status.progressEstimateAvailable = TRUE;
    status.estimatedFraction = 0.5f;

    QgCanvasState state;
    state.lod_progress_overlay = qgcanvas_lod_progress_overlay_state(status);
    state.lod_progress_overlay_dirty = true;
    {
	QPainter painter(&painted);
	qgcanvas_paint_lod_progress_overlay(state, &widget, painter);
    }
    if (state.lod_progress_overlay_dirty)
	FAIL("native LoD overlay paint did not consume dirty state");
    if (image_byte_diff(baseline, painted) <= 0)
	FAIL("native LoD overlay paint did not modify the framebuffer");
    if (painted.pixelColor(painted.width() - 1, painted.height() - 1) !=
	QColor(Qt::black))
	FAIL("native LoD overlay paint escaped its upper-left card bounds");

    const auto paintedBounds = [&baseline](const QImage &image) {
	QRect bounds;
	for (int y = 0; y < image.height(); ++y) {
	    for (int x = 0; x < image.width(); ++x) {
		if (image.pixelColor(x, y) != baseline.pixelColor(x, y))
		    bounds = bounds.united(QRect(x, y, 1, 1));
	    }
	}
	return bounds;
    };
    QImage shortText = baseline;
    state.lod_progress_overlay.title = QStringLiteral("Working");
    state.lod_progress_overlay.detail = QStringLiteral("50%");
    state.lod_progress_overlay_dirty = true;
    {
	QPainter painter(&shortText);
	qgcanvas_paint_lod_progress_overlay(state, &widget, painter);
    }
    if (paintedBounds(shortText).size() != paintedBounds(painted).size())
	FAIL("native LoD overlay card size changed with its text");

    const QColor accent = qgcanvas_lod_overlay_color(
	state.lod_progress_overlay);
    int accentPixels = 0;
    for (int y = 0; y < painted.height(); ++y) {
	for (int x = 0; x < painted.width(); ++x) {
	    if (painted.pixelColor(x, y) == accent)
		accentPixels++;
	}
    }
    if (accentPixels < 4)
	FAIL("native LoD overlay did not paint its progress indication");

    QImage hidden = baseline;
    state.lod_progress_overlay.visible = false;
    state.lod_progress_overlay_dirty = true;
    {
	QPainter painter(&hidden);
	qgcanvas_paint_lod_progress_overlay(state, &widget, painter);
    }
    if (state.lod_progress_overlay_dirty ||
	image_byte_diff(baseline, hidden) != 0)
	FAIL("hidden native LoD overlay modified the framebuffer");
    return 0;
}

static int
check_native_lod_overlay_repaint_priority(void)
{
    QWidget widget;
    QgCanvasState state;
    state.lod_progress_overlay_dirty = true;
    state.lod_progress_overlay.visible = true;
    state.lod_progress_overlay_last_request =
	std::chrono::steady_clock::now() + std::chrono::hours(1);
    const std::chrono::steady_clock::time_point throttled =
	state.lod_progress_overlay_last_request;
    qgcanvas_request_lod_overlay_repaint(state, &widget);
    if (state.lod_progress_overlay_last_request != throttled)
	FAIL("ordinary native overlay repaint bypassed its cadence limit");

    /* Hiding the last visible card is a terminal state transition, not a
     * periodic animation sample.  It must survive a preceding update inside
     * the cadence window after the progressive timer retires. */
    state.lod_progress_overlay.visible = false;
    qgcanvas_request_lod_overlay_repaint(state, &widget);
    if (state.lod_progress_overlay_last_request == throttled)
	FAIL("native overlay hide was suppressed by the repaint throttle");
    return 0;
}

static bool
render_hud_frame(BObolViewController *controller)
{
    unsigned char *pixels = NULL;
    const int result = controller->renderToImage(&pixels);
    const bool rendered = result == BRLCAD_OK && pixels;
    bu_free(pixels, "HUD publication regression frame");
    if (rendered)
	controller->noteFramePresented();
    return rendered;
}

static int
check_terminal_error_hud(QgView &view, SoBRLDatabaseSource *source)
{
    BObolViewController *controller = view.obolViewController();
    struct ged_view_context *view_ctx =
	ged_view_context_from_bv(view.viewContext());
    /* This hidden test canvas receives no resize event.  GED needs the same
     * physical viewport as the renderer to position its retained HUD. */
    const SbVec2s viewport_size =
	controller->getViewportRegion().getViewportSizePixels();
    (void)bv_context_dimensions_set(view.viewContext(),
	viewport_size[0], viewport_size[1]);
    ged_view_lod_policy original_policy = BV_LOD_POLICY_INIT;
    if (!ged_view_lod_policy_get(&original_policy, view_ctx))
	FAIL("HUD regression should read the view policy");
    struct ged *gedp = ged_view_context_owner(view_ctx);
    const char *healthy_paths[] = {"box.s"};
    struct ged_scene_draw_request healthy_draw;
    ged_scene_draw_request_init(&healthy_draw);
    healthy_draw.view = view_ctx;
    healthy_draw.paths = healthy_paths;
    healthy_draw.path_count = 1;
    healthy_draw.realization.mode = GED_SCENE_REALIZE_EAGER;
    healthy_draw.style.draw_mode = GED_SCENE_DRAW_SHADED;
    struct ged_scene_result *result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_draw(gedp, &healthy_draw, result),
	result, 1))
	FAIL("HUD regression should retain a healthy sibling beside the failed source");
    const SbBool original_auto_submit = controller->isLodAutoSubmitEnabled();
    ged_view_lod_policy policy = original_policy;
    policy.policy = BV_LOD_AUTO;
    policy.mesh_enabled = 1;
    policy.csg_enabled = 1;
    if (!ged_view_lod_policy_apply(view_ctx, &policy))
	FAIL("HUD regression should enable automatic LoD");
    const bool needs_service = controller->getLodService() == NULL;
    if (needs_service && !controller->ensureManagedLodService(1))
	FAIL("HUD regression should provide the service required by the planner");
    controller->setLodAutoSubmit(TRUE);
    source->realizationStatus = SoBRLDatabaseSource::FAILED;

    const uint64_t provider = controller->registerProgressiveProvider(
	pending_hud_provider, NULL);
    if (!provider)
	FAIL("HUD regression should retain an independent preparation provider");
    (void)controller->advanceProgressiveWork();
    /* Earlier camera tests may still own an unbracketed interaction.  Retire
     * that setup work before initializing the publication cache; the tested
     * provider transition below does not depend on the quiet-view clock. */
    const auto setup_deadline = std::chrono::steady_clock::now() +
	std::chrono::seconds(2);
    while (controller->isLodInteractionActive()) {
	if (controller->isRenderRequested() && !render_hud_frame(controller))
	    FAIL("HUD regression should present its initial view");
	(void)controller->advanceProgressiveWork();
	if (std::chrono::steady_clock::now() >= setup_deadline)
	    FAIL("HUD regression should begin with a quiet view");
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    /* This test deliberately validates renderer-neutral retained records.
     * Production QgView hosts select the native Qt overlay instead. */
    if (!ged_view_lod_progress_presentation_mode_set(view_ctx,
	    GED_VIEW_LOD_PROGRESS_PRESENTATION_RETAINED))
	FAIL("HUD regression should select retained presentation ownership");

    /* Use the production publication helper with its own cache, and omit
     * periodic samples.  No event-loop or debugger timing can accidentally
     * publish the terminal sample on behalf of the transition under test. */
    QgCanvasState state;
    state.v = view.viewContext();
    state.obol = controller;
    BObolLodConvergenceStatus status;
    controller->getLodConvergenceStatus(status);
    if (status.phase != BOBOL_LOD_CONVERGENCE_ERROR || status.terminal ||
	!status.sourcePreparationPending ||
	!qgcanvas_sync_obol_lod_progress(state, false))
	FAIL("an active provider should keep a published geometry error nonterminal");

    controller->unregisterProgressiveProvider(provider);
    /* Bound control/render effects for these two small realized sources by
     * steps, independently of wall-clock speed. */
    static constexpr unsigned int maximumRetirementSteps = 16;
    unsigned int steps = 0;
    for (; steps < maximumRetirementSteps; ++steps) {
	(void)controller->advanceProgressiveWork();
	if (controller->isRenderRequested() && !render_hud_frame(controller))
	    FAIL("HUD regression should render the remaining presentation");
	controller->getLodConvergenceStatus(status);
	if (status.terminalError)
	    break;
    }
    if (!status.terminalError || status.viewReady || status.fraction < 1.0f) {
	fprintf(stderr, "HUD retirement: steps=%u phase=%d terminal=%d "
	    "errors=%u fraction=%g owner=%d obligations=%u pending=%d render=%d\n",
	    steps, status.phase, status.terminal, status.failedSourceCount,
	    status.fraction, status.controlOwner, status.controlObligationMask,
	    controller->hasProgressiveWorkPending(), controller->isRenderRequested());
	FAIL("retiring the provider should reach an honest terminal error");
    }
    if (!qgcanvas_sync_obol_lod_progress(state, false))
	FAIL("terminal error must publish its final HUD without a periodic sample");

    BObolFeatureRecord fill;
    BObolFeatureRecord label;
    BObolFeatureRecord track;
    const BObolFeatureStore &features = controller->features();
    static constexpr float pixelTolerance = 1.0e-4f;
    if (!features.record(features.find("_faceplate/lod_progress_fill"), fill) ||
	!features.record(features.find("_faceplate/lod_progress_label"), label) ||
	!features.record(features.find("_faceplate/lod_progress_track"), track) ||
	fill.points.size() != 2 || track.points.size() != 2 ||
	label.labels.size() != 1 ||
	label.labels[0].text != "View incomplete  1 geometry error" ||
	fill.style.color != SbColor(1.0f, 90.0f / 255.0f, 80.0f / 255.0f) ||
	fill.points[1][1] <= fill.points[0][1] ||
	fabsf(track.points[0][1] - track.points[1][1]) > pixelTolerance ||
	fabsf(fill.points[1][1] - track.points[0][1]) > pixelTolerance) {
	fprintf(stderr, "HUD records: fill=%zu layers=%zu track=%zu labels=%zu\n",
	    fill.points.size(), fill.layers.size(), track.points.size(), label.labels.size());
	if (fill.points.size() == 2 && track.points.size() == 2 && !label.labels.empty())
	    fprintf(stderr, "HUD details: fillY=%g,%g trackY=%g,%g color=%g,%g,%g label=%s\n",
		fill.points[0][1], fill.points[1][1], track.points[0][1], track.points[1][1],
		fill.style.color[0], fill.style.color[1], fill.style.color[2],
		label.labels[0].text.getString());
	FAIL("terminal error should retain a full red bar, terminal cap and error label");
    }

    enum ged_view_lod_progress_presentation_mode presentationMode =
	GED_VIEW_LOD_PROGRESS_PRESENTATION_NONE;
    if (!ged_view_lod_progress_presentation_mode_get(&presentationMode,
	    view_ctx) || presentationMode !=
	    GED_VIEW_LOD_PROGRESS_PRESENTATION_RETAINED)
	FAIL("retained LoD HUD ownership should be observable");
    if (!ged_view_lod_progress_presentation_mode_set(view_ctx,
	    GED_VIEW_LOD_PROGRESS_PRESENTATION_NATIVE_HOST) ||
	features.exists("_faceplate/lod_progress_fill") ||
	features.exists("_faceplate/lod_progress_label") ||
	features.exists("_faceplate/lod_progress_track") ||
	ged_view_faceplate_sync(gedp, view_ctx) != BRLCAD_OK ||
	features.exists("_faceplate/lod_progress_fill") ||
	features.exists("_faceplate/lod_progress_label") ||
	features.exists("_faceplate/lod_progress_track"))
	FAIL("native LoD HUD ownership should atomically suppress retained records");
    if (!ged_view_lod_progress_presentation_mode_set(view_ctx,
	    GED_VIEW_LOD_PROGRESS_PRESENTATION_RETAINED) ||
	!features.exists("_faceplate/lod_progress_fill") ||
	!features.exists("_faceplate/lod_progress_label") ||
	!features.exists("_faceplate/lod_progress_track"))
	FAIL("retained LoD HUD ownership should republish all progress records");

    const uint64_t feature_revision = features.presentationRevision();
    if (qgcanvas_sync_obol_lod_progress(state, false) ||
	features.presentationRevision() != feature_revision)
	FAIL("a pending HUD frame must not reopen the error publication transition");
    if (!controller->isRenderRequested() || !render_hud_frame(controller))
	FAIL("the final error HUD should own a presentation frame");
    if (qgcanvas_sync_obol_lod_progress(state, false) ||
	features.presentationRevision() != feature_revision ||
	controller->isRenderRequested() || controller->hasProgressiveWorkPending())
	FAIL("the final error HUD should retire without another publication or repaint");
    printf("Terminal error HUD publishes its full bar and retires without repaint\n");

    source->realizationStatus = SoBRLDatabaseSource::REALIZED;
    if (!ged_view_lod_policy_apply(view_ctx, &original_policy))
	FAIL("HUD regression should restore the view policy");
    controller->setLodAutoSubmit(original_auto_submit);
    if (needs_service)
	controller->stopManagedLodService();
    struct ged_scene_erase_request healthy_erase;
    ged_scene_erase_request_init(&healthy_erase);
    healthy_erase.view = view_ctx;
    healthy_erase.path = "box.s";
    result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_erase(gedp, &healthy_erase, result),
	result, 1))
	FAIL("HUD regression should remove its healthy sibling fixture");
    return 0;
}

static int
production_flow_view_open(struct ged *gedp,
	std::unique_ptr<QgView> &view,
	struct draw_observer_sync_context &observer,
	ged_scene_observer_token &observer_token)
{
    static constexpr int viewport_width = 180;
    static constexpr int viewport_height = 140;

    view = std::make_unique<QgView>(nullptr, QgViewType::SW);
    view->resize(viewport_width, viewport_height);
    struct ged_view_context *view_ctx =
	ged_view_context_from_bv(view->viewContext());
    if (!view_ctx || !view->isValid() ||
	!ged_view_set_context_add(ged_view_set_ctx(gedp), view_ctx) ||
	!ged_view_context_host_attach(gedp, view_ctx))
	return 0;

    ged_view_active_ctx_set(gedp, view_ctx);
    (void)bv_unit_conversion_set(bv_context_view(view->viewContext()),
	gedp->dbip->dbi_local2base, gedp->dbip->dbi_base2local);
    if (!qg_scene_bind(gedp, view.get()))
	return 0;

    observer.view = view.get();
    observer.calls = 0;
    observer.changed = 0;
    observer_token = ged_scene_observer_add(gedp,
	test_draw_observer, &observer);
    return observer_token ? 1 : 0;
}

static void
production_flow_view_close(struct ged *gedp,
	std::unique_ptr<QgView> &view,
	ged_scene_observer_token &observer_token)
{
    if (observer_token) {
	(void)ged_scene_observer_remove(gedp, observer_token);
	observer_token = 0;
    }
    ged_view_active_ctx_set(gedp, NULL);
    view.reset();
}

static int
production_flow_mesh_policy_enable(struct ged *gedp)
{
    const char *lod_enable[3] = {"view", "lod", "1"};
    const char *mesh_enable[4] = {"view", "lod", "mesh", "1"};
    const char *bot_threshold[4] = {
	"view", "lod", "bot_threshold", "0"
    };
    return ged_exec_view(gedp, 3, lod_enable) == BRLCAD_OK &&
	ged_exec_view(gedp, 4, mesh_enable) == BRLCAD_OK &&
	ged_exec_view(gedp, 4, bot_threshold) == BRLCAD_OK;
}

static bool
production_flow_service_idle(const BObolViewController *controller)
{
    BObolLodService *service = controller ? controller->getLodService() : NULL;
    return service && service->workStatus().isIdle() &&
	!controller->hasPendingLodResults() &&
	!controller->hasPendingLodSubmissions();
}

static BObolLodResult
production_flow_delayed_result(const BObolLodRequest &request, void *)
{
    BObolLodResult result;
    result.request = request;
    result.resultKind = BOBOL_LOD_RESULT_AABB;
    result.qualityTier = request.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.terminal = TRUE;
    return result;
}

static int
exercise_production_graphical_lifecycle(struct ged *gedp)
{
    static constexpr const char *mesh_path = "flow.bot";
    static constexpr fastf_t loading_camera_scale = 18.0;
    static constexpr int task_delay_milliseconds = 250;
    static constexpr int delayed_task_timeout_seconds = 4;
    static constexpr int worker_retirement_timeout_seconds = 2;
    static constexpr int terminal_timeout_seconds = 8;
    static constexpr int minimum_lit_pixels = 10;

    std::unique_ptr<QgView> view;
    struct draw_observer_sync_context observer = {NULL, 0, 0};
    ged_scene_observer_token observer_token = 0;
    auto fail = [&](const char *message) {
	production_flow_view_close(gedp, view, observer_token);
	fprintf(stderr, "FAIL: %s\n", message);
	return 1;
    };

    if (!production_flow_view_open(gedp, view, observer, observer_token) ||
	!production_flow_mesh_policy_enable(gedp))
	return fail("graphical production flow should open its first LoD view");

    BObolViewController *controller = view->obolViewController();
    struct ged_view_context *view_ctx =
	ged_view_context_from_bv(view->viewContext());
    const char *draw_mesh[3] = {"draw", "-m1", mesh_path};
    if (!controller || ged_exec_draw(gedp, 3, draw_mesh) != BRLCAD_OK ||
	observer.calls <= 0 || observer.changed <= 0)
	return fail("graphical production flow should draw through the GED observer");
    if (!controller->syncCameraFromViewContext(view_ctx))
	return fail("graphical production flow should sync its initial camera");

    controller->requestLodCapacityRender("graphical-production-first-frame");
    view->need_update(QG_VIEW_REFRESH);
    QCoreApplication::processEvents();
    (void)controller->advanceProgressiveWork(NULL, NULL);
    QImage first_frame;
    (void)qtcad_obol_present_requested_frame(*view, controller,
	    first_frame);
    if (first_frame.isNull() ||
	lit_pixel_count(first_frame) < minimum_lit_pixels)
	return fail("graphical production flow should present a visible first frame");

    BObolLodService *service = controller->getLodService();
    /* The headless half of FLOW-01 closes during a real delayed mesh task.
     * Here an explicit service delay isolates the Qt ownership boundary from
     * view-planner policy: destroying QgView must stop its managed workers
     * regardless of whether this small model needs another visible cut. */
    BObolLodTask delayed_task;
    delayed_task.generation = controller->beginLodGeneration();
    delayed_task.request.databaseId = "db://qtcad-production-flow";
    delayed_task.request.databaseRevision = 1;
    delayed_task.request.sourceRevision = 1;
    delayed_task.request.objectPath = mesh_path;
    delayed_task.request.objectName = mesh_path;
    delayed_task.request.providerId = "qtcad-production-flow";
    delayed_task.request.providerVersion = "1";
    delayed_task.request.qualityTier = BOBOL_LOD_QUALITY_PROXY;
    delayed_task.realize = production_flow_delayed_result;
    delayed_task.debugDelayMilliseconds = task_delay_milliseconds;
    delayed_task.publishResult = FALSE;
    if (!service || !delayed_task.generation ||
	service->submit(delayed_task) == 0)
	return fail("graphical production flow should submit bounded delayed work");

    bool delayed_task_observed = false;
    const auto delayed_deadline = std::chrono::steady_clock::now() +
	std::chrono::seconds(delayed_task_timeout_seconds);
    do {
	QCoreApplication::processEvents();
	(void)controller->advanceProgressiveWork(NULL, NULL);
	service = controller->getLodService();
	if (service && service->delayedTaskCountForDiagnostics() > 0) {
	    delayed_task_observed = true;
	    break;
	}
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    } while (std::chrono::steady_clock::now() < delayed_deadline);

    SoBRLDatabaseSource *loading_source = source_for_path(controller,
	mesh_path);
    if (!delayed_task_observed || !loading_source)
	return fail("graphical production flow should expose delayed mesh work");

    struct bv *loading_view = bv_context_view(view->viewContext());
    const uint64_t camera_revision = bv_frame_revision_get(loading_view);
    bv_scale_set(loading_view, loading_camera_scale);
    view->need_update(QG_VIEW_REFRESH);
    QCoreApplication::processEvents();
    if (bv_frame_revision_get(loading_view) <= camera_revision ||
	!controller->syncCameraFromViewContext(view_ctx))
	return fail("graphical production flow should accept camera input while loading");

    const size_t closing_worker_count = controller->getManagedLodWorkerCount();
#if defined(__linux__)
    const size_t closing_thread_count =
	bu_file_list("/proc/self/task", "[0-9]*", NULL);
#endif
    if (closing_worker_count == 0)
	return fail("graphical production flow should close a worker-active view");
    production_flow_view_close(gedp, view, observer_token);
#if defined(__linux__)
    bool workers_retired = false;
    const auto worker_deadline = std::chrono::steady_clock::now() +
	std::chrono::seconds(worker_retirement_timeout_seconds);
    do {
	const size_t current_threads =
	    bu_file_list("/proc/self/task", "[0-9]*", NULL);
	workers_retired = current_threads + closing_worker_count <=
	    closing_thread_count;
	if (!workers_retired)
	    std::this_thread::sleep_for(std::chrono::milliseconds(1));
    } while (!workers_retired &&
	std::chrono::steady_clock::now() < worker_deadline);
    if (!workers_retired)
	return fail("destroying the graphical view should retire its LoD workers");
#endif

    if (!production_flow_view_open(gedp, view, observer, observer_token) ||
	!production_flow_mesh_policy_enable(gedp))
	return fail("graphical production flow should reopen its LoD view");

    controller = view->obolViewController();
    view_ctx = ged_view_context_from_bv(view->viewContext());
    struct bv *reopened_view = bv_context_view(view->viewContext());
    bv_scale_set(reopened_view, loading_camera_scale);
    view->need_update(QG_VIEW_REFRESH);
    QCoreApplication::processEvents();
    if (!controller || !controller->syncCameraFromViewContext(view_ctx) ||
	!source_for_path(controller, mesh_path))
	return fail("reopened graphical view should consume the retained GED scene");

    controller->requestLodCapacityRender("graphical-production-reopen-baseline");
    QImage baseline_frame;
    (void)qtcad_obol_present_requested_frame(*view, controller,
	    baseline_frame);
    if (baseline_frame.isNull())
	return fail("reopened graphical view should present its retained baseline");

    observer.calls = 0;
    observer.changed = 0;
    const char *redraw[1] = {"redraw"};
    if (ged_exec_redraw(gedp, 1, redraw) != BRLCAD_OK ||
	observer.calls <= 0 || observer.changed <= 0)
	return fail("reopened graphical view should redraw through the GED observer");

    BObolProgressiveOptions options;
    options.forceTerminalLodRefinement = TRUE;
    BObolProgressiveStatus progress;
    BObolLodConvergenceStatus convergence;
    QImage terminal_frame;
    const auto terminal_deadline = std::chrono::steady_clock::now() +
	std::chrono::seconds(terminal_timeout_seconds);
    do {
	QCoreApplication::processEvents();
	(void)controller->advanceProgressiveWork(&options, &progress);
	controller->requestLodCapacityRender(
	    "graphical-production-terminal-frame");
	QImage candidate;
	if (qtcad_obol_present_requested_frame(*view, controller, candidate))
	    terminal_frame = candidate;
	controller->getLodConvergenceStatus(convergence);
	if (convergence.terminal && production_flow_service_idle(controller))
	    break;
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    } while (std::chrono::steady_clock::now() < terminal_deadline);

    SoBRLDatabaseSource *terminal_source = source_for_path(controller,
	mesh_path);
    const int pixel_difference = image_byte_diff(baseline_frame,
	terminal_frame);
    if (!convergence.terminal || !convergence.viewReady ||
	convergence.terminalError || !production_flow_service_idle(controller) ||
	!terminal_source || !terminal_source->isCompactOccurrenceRegistry() ||
	controller->getActiveLodMeshPayloadCount() == 0 ||
	terminal_frame.isNull() ||
	lit_pixel_count(terminal_frame) < minimum_lit_pixels ||
	pixel_difference <= 0) {
	fprintf(stderr,
	    "graphical production terminal=%d ready=%d error=%d source=%p "
	    "compact=%d meshes=%zu lit=%d diff=%d service_idle=%d\n",
	    convergence.terminal, convergence.viewReady,
	    convergence.terminalError, static_cast<void *>(terminal_source),
	    terminal_source ? terminal_source->isCompactOccurrenceRegistry() : -1,
	    controller->getActiveLodMeshPayloadCount(),
	    terminal_frame.isNull() ? 0 : lit_pixel_count(terminal_frame),
	    pixel_difference, production_flow_service_idle(controller));
	return fail("reopened graphical view should reach terminal output and release transient work");
    }

    const char *erase_mesh[2] = {"erase", mesh_path};
    if (ged_exec_erase(gedp, 2, erase_mesh) != BRLCAD_OK ||
	source_for_path(controller, mesh_path))
	return fail("graphical production flow should erase its retained source");

    production_flow_view_close(gedp, view, observer_token);
    std::puts("PASS QgView production flow lifecycle: draw, camera input, close, reopen, terminal image, resource release");
    return 0;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    QApplication app(argc, argv);

    if (check_native_lod_overlay_state())
	return 1;
    if (check_native_lod_overlay_stability())
	return 1;
    if (check_native_lod_overlay_paint())
	return 1;
    if (check_native_lod_overlay_repaint_priority())
	return 1;

    const bool production_flow = argc > 1 && BU_STR_EQUAL(argv[1],
	"production-flow-lifecycle");
    char production_cache_leaf[64] = {0};
    char production_cache_dir[MAXPATHLEN] = {0};
    if (production_flow) {
	snprintf(production_cache_leaf, sizeof(production_cache_leaf),
	    "qtcad_obol_draw_sync_cache_%d", bu_pid());
	bu_dir(production_cache_dir, MAXPATHLEN, BU_DIR_CURR,
	    production_cache_leaf, NULL);
	bu_dirclear(production_cache_dir);
	bu_mkdir(production_cache_dir);
	bu_setenv("BU_DIR_CACHE", production_cache_dir, 1);
    }

    const char *dbpath = "qtcad_obol_draw_sync_tmp.g";
    if (!make_draw_sync_db(dbpath))
	FAIL("failed to create qtcad Obol draw-sync test database");

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open qtcad Obol draw-sync test database");

    if (production_flow) {
	const int result = exercise_production_graphical_lifecycle(gedp);
	ged_close(gedp);
	bu_file_delete(dbpath);
	bu_dirclear(production_cache_dir);
	return result;
    }

    QgView view(NULL, QgViewType::SW);
    view.resize(180, 140);

    /* A renderer attachment must consume the complete current semantic
     * snapshot, including visibility-frontier changes committed while the
     * scene was headless. */
    struct ged_view_context *headless_view = ged_view_context_create();
    if (!headless_view || !ged_view_context_host_attach(gedp, headless_view))
	FAIL("headless snapshot test should create an attached view context");
    ged_view_active_ctx_set(gedp, headless_view);
    const char *headless_paths[] = {"pair.c"};
    struct ged_scene_draw_request headless_draw;
    ged_scene_draw_request_init(&headless_draw);
    headless_draw.view = headless_view;
    headless_draw.paths = headless_paths;
    headless_draw.path_count = 1;
    headless_draw.realization.mode = GED_SCENE_REALIZE_EAGER;
    struct ged_scene_result *headless_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_draw(gedp, &headless_draw,
	    headless_result), headless_result, 1))
	FAIL("headless draw intent should commit without a renderer");
    struct ged_scene_path_request headless_path;
    ged_scene_path_request_init(&headless_path);
    headless_path.view = headless_view;
    headless_path.path = "box.s";
    headless_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_visibility_set(gedp,
	    &headless_path, 0, headless_result), headless_result, 0))
	FAIL("exact compact path matching must not fall back to a leaf basename");
    headless_path.path = "pair.c/ball.s";
    headless_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_opacity_set(gedp, &headless_path,
	    0.25, headless_result), headless_result, 1))
	FAIL("headless compact opacity should be retained before realization");
    headless_path.path = "ball.s";
    headless_path.match = GED_SCENE_PATH_MATCH_OBJECT;
    headless_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_highlight_set(gedp, &headless_path,
	    1, headless_result), headless_result, 1))
	FAIL("headless object highlight should be retained before realization");
    struct ged_scene_erase_request headless_erase;
    ged_scene_erase_request_init(&headless_erase);
    headless_erase.view = headless_view;
    headless_erase.path = "pair.c/box.s";
    headless_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_erase(gedp, &headless_erase,
	    headless_result), headless_result, 1))
	FAIL("headless subpath erase should commit a visibility frontier");
    BObolViewController *controller = view.obolViewController();
    if (!controller || render_source_count(controller) != 0)
	FAIL("headless semantic commits must not construct renderer sources");
    if (!ged_view_context_obol_endpoint_set(headless_view,
	    view.displayEndpoint(), 0))
	FAIL("delayed endpoint attachment should accept a semantic snapshot");
    SoBRLDatabaseSource *headless_source =
	source_for_path(controller, "pair.c");
    if (!headless_source || headless_source->getCompactInstanceCount() != 2 ||
	    visible_compact_occurrence_count(headless_source) != 1)
	FAIL("delayed endpoint snapshot should preserve the headless erase frontier");
    BObolCompactOccurrence headless_ball;
    if (!compact_occurrence_for_path(headless_source, "pair.c/ball.s",
	    &headless_ball) || !headless_ball.summary.highlighted ||
	    fabsf(headless_ball.summary.transparency - 0.75f) > 1.0e-6f) {
	fprintf(stderr, "headless ball found=%d highlighted=%d transparency=%g\n",
	    compact_occurrence_for_path(headless_source, "pair.c/ball.s", NULL),
	    headless_ball.summary.highlighted ? 1 : 0,
	    headless_ball.summary.transparency);
	FAIL("delayed endpoint snapshot should preserve compact highlight and opacity");
    }

    /* Addressing one compact occurrence must remain available to editing and
     * picking clients without realizing sibling scene records. */
    struct ged_scene_path_request occurrence_request;
    ged_scene_path_request_init(&occurrence_request);
    occurrence_request.view = headless_view;
    occurrence_request.path = "pair.c/ball.s";
    occurrence_request.match = GED_SCENE_PATH_MATCH_EXACT;
    ged_scene_occurrence_ref occurrence = ged_scene_occurrence_resolve(gedp,
	&occurrence_request);
    struct ged_scene_occurrence_info occurrence_info;
    if (ged_scene_occurrence_ref_is_null(occurrence) ||
	!ged_scene_occurrence_get(gedp, occurrence, &occurrence_info) ||
	!occurrence_info.path ||
	!BU_STR_EQUAL(test_skip_leading_slash(occurrence_info.path),
	    "pair.c/ball.s") || render_source_count(controller) != 1)
	FAIL("exact compact occurrence resolution should retain one source and one semantic path");

    occurrence_request.path = "pair.c";
    occurrence_request.match = GED_SCENE_PATH_MATCH_SUBTREE;
    occurrence = ged_scene_occurrence_resolve(gedp, &occurrence_request);
    if (ged_scene_occurrence_ref_is_null(occurrence) ||
	!ged_scene_occurrence_get(gedp, occurrence, &occurrence_info) ||
	!occurrence_info.path ||
	!BU_STR_EQUAL(test_skip_leading_slash(occurrence_info.path),
	    "pair.c/ball.s"))
	FAIL("compact subtree occurrence resolution should skip hidden children in deterministic order");

    occurrence_request.path = "pair.c/missing.s";
    occurrence_request.match = GED_SCENE_PATH_MATCH_EXACT;
    if (!ged_scene_occurrence_ref_is_null(
	    ged_scene_occurrence_resolve(gedp, &occurrence_request)))
	FAIL("missing compact paths must not manufacture occurrence records");

    ged_scene_path_request_init(&headless_path);
    headless_path.view = headless_view;
    headless_path.path = "pair.c/ball.s";
    headless_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_visibility_set(gedp,
	    &headless_path, 0, headless_result), headless_result, 1))
	FAIL("compact visibility should update retained scene intent");
    headless_source = source_for_path(controller, "pair.c");
    if (!headless_source || visible_compact_occurrence_count(headless_source))
	FAIL("attached endpoint should apply compact visibility without leaf records");
    if (!compact_occurrence_for_path(headless_source, "pair.c/ball.s",
	    &headless_ball) || !headless_ball.summary.highlighted ||
	    fabsf(headless_ball.summary.transparency - 0.75f) > 1.0e-6f)
	FAIL("compact visibility must preserve prior highlight and opacity state");
    (void)ged_view_context_obol_endpoint_set(headless_view, NULL, 0);
    struct ged_scene_clear_request headless_clear;
    ged_scene_clear_request_init(&headless_clear);
    headless_clear.view = headless_view;
    headless_clear.scope = GED_SCENE_CLEAR_VIEW;
    headless_result = ged_scene_result_create();
    if (ged_scene_clear(gedp, &headless_clear, headless_result) != GED_SCENE_OK)
	FAIL("headless snapshot test should clear its view-scoped draw intent");
    ged_scene_result_destroy(headless_result);
    ged_view_context_free(headless_view);
    controller->clearDatabaseSources();

    struct ged_view_context *ged_view_ctx =
	ged_view_context_from_bv(view.viewContext());
	ged_view_active_ctx_set(gedp, ged_view_ctx);
	(void)ged_view_context_host_attach(gedp, ged_view_ctx);
    controller->clearDatabaseSources();

    struct draw_observer_sync_context obs = {&view, 0, 0};
    ged_scene_observer_token observerToken =
	ged_scene_observer_add(gedp, test_draw_observer, &obs);
    if (!observerToken)
	FAIL("GED draw observer should register for qtcad Obol draw sync");

    const char *draw_cmd[2] = {"draw", "box.s"};
    if (ged_exec_draw(gedp, 2, draw_cmd) != BRLCAD_OK)
	FAIL("real GED draw command should succeed for observer sync");
    if (obs.calls <= 0 || obs.changed <= 0)
	FAIL("real GED draw command should notify and sync qtcad Obol");
    if (render_source_count(controller) != 1)
	FAIL("observer-synced GED draw should create one Obol database source");
    SoBRLDatabaseSource *observerSource = source_for_path(controller, "box.s");
    if (!observerSource ||
	    observerSource->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	    !observerSource->hasRealizedWireGeometry())
	FAIL("observer-synced GED draw should realize Obol wire geometry");

    controller->getViewport()->viewAll();
    controller->requestLodCapacityRender("observer-draw-visible");
    QCoreApplication::processEvents();
    QImage observerImage;
    view.get_viewport_image(observerImage);
    int observerLit = observerImage.isNull() ? 0 : lit_pixel_count(observerImage);
    if (observerImage.isNull() || observerLit < 10) {
	fprintf(stderr,
		"FAIL: observer-synced GED draw should be visible through qtcad capture (null=%d lit=%d)\n",
		observerImage.isNull() ? 1 : 0, observerLit);
	return 1;
    }

    SoCamera *camera = controller->getCamera();
    if (!camera)
	FAIL("qtcad Obol controller should expose a camera for view sync");
    point_t offcenter = {100.0, 0.0, 0.0};
    bv_center_set(bv_context_view(static_cast<struct bv_context *>(view.viewContext())), offcenter);
    bv_scale_set(bv_context_view(static_cast<struct bv_context *>(view.viewContext())), 250.0);
    view.need_update(QG_VIEW_REFRESH);
    SbVec3f offTargetCamera = camera_focal_point(camera);
    if (offTargetCamera[0] < 50.0f)
	FAIL("qtcad refresh should sync Obol camera from GED view state");

    const char *autoview_cmd[1] = {"autoview"};
    if (ged_exec_autoview(gedp, 1, autoview_cmd) != BRLCAD_OK)
	FAIL("real GED autoview command should succeed for qtcad Obol view sync");
    view.need_update(QG_VIEW_REFRESH);
    SbVec3f autoviewCamera = camera_focal_point(camera);
    if (!camera_positions_differ(offTargetCamera, autoviewCamera, 10.0f))
	FAIL("GED autoview should update qtcad Obol camera through view refresh");
    if (fabsf(autoviewCamera[0]) > 25.0f)
	FAIL("GED autoview should recenter qtcad Obol camera near drawn geometry");
    controller->requestLodCapacityRender("observer-autoview-visible");
    QCoreApplication::processEvents();
    QImage autoviewImage;
    view.get_viewport_image(autoviewImage);
    if (autoviewImage.isNull() || lit_pixel_count(autoviewImage) < 10)
	FAIL("GED-autoviewed Obol scene should stay visible through qtcad capture");

    obs.calls = 0;
    obs.changed = 0;
    const char *erase_cmd[2] = {"erase", "box.s"};
    if (ged_exec_erase(gedp, 2, erase_cmd) != BRLCAD_OK)
	FAIL("real GED erase command should succeed for observer sync");
    if (obs.calls <= 0 || obs.changed <= 0)
	FAIL("real GED erase command should notify and sync qtcad Obol");
    if (render_source_count(controller) != 0)
	FAIL("observer-synced GED erase should remove Obol database sources");

    /* A nested erase/redraw is a visibility-frontier edit of the retained
     * top-level source.  It deliberately does not rebuild that source, but it
     * must invalidate the completed framebuffer or Qt will replay obsolete
     * pixels forever.  Exercise the actual GED observer bridge used by qged. */
    const char *observer_pair_paths[] = {"pair.c"};
    struct ged_scene_draw_request observer_pair_draw;
    ged_scene_draw_request_init(&observer_pair_draw);
    observer_pair_draw.view = ged_view_ctx;
    observer_pair_draw.paths = observer_pair_paths;
    observer_pair_draw.path_count = 1;
    observer_pair_draw.realization.mode = GED_SCENE_REALIZE_EAGER;
    struct ged_scene_result *observer_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_draw(gedp, &observer_pair_draw,
	    observer_result), observer_result, 1))
	FAIL("observer retained-frontier setup draw should succeed");
    SoBRLDatabaseSource *pairSource = source_for_path(controller, "pair.c");
    if (!pairSource || visible_compact_occurrence_count(pairSource) != 2)
	FAIL("observer retained-frontier setup should expose both occurrences");
    const uint64_t pairMeshInventoryRevision =
	pairSource->getDisplayMeshLodRevision();

    controller->clearRenderRequest();
    uint64_t frontierSerial = controller->renderRequestSerialGet();
    obs.calls = 0;
    obs.changed = 0;
    struct ged_scene_erase_request observer_nested_erase;
    ged_scene_erase_request_init(&observer_nested_erase);
    observer_nested_erase.view = ged_view_ctx;
    observer_nested_erase.path = "pair.c/box.s";
    observer_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_erase(gedp,
	    &observer_nested_erase, observer_result), observer_result, 1))
	FAIL("observer nested erase should succeed");
    SbBool frontierCapacityRelevant = TRUE;
    if (obs.calls <= 0 || obs.changed <= 0 ||
	visible_compact_occurrence_count(pairSource) != 1 ||
	pairSource->getDisplayMeshLodRevision() !=
	    pairMeshInventoryRevision ||
	controller->renderRequestSerialGet() <= frontierSerial ||
	!controller->consumeRenderRequest(NULL,
	    &frontierCapacityRelevant) || frontierCapacityRelevant)
	FAIL("observer nested erase should request a presentation-only frame");

    frontierSerial = controller->renderRequestSerialGet();
    obs.calls = 0;
    obs.changed = 0;
    const char *observer_nested_paths[] = {"pair.c/box.s"};
    struct ged_scene_draw_request observer_nested_draw;
    ged_scene_draw_request_init(&observer_nested_draw);
    observer_nested_draw.view = ged_view_ctx;
    observer_nested_draw.paths = observer_nested_paths;
    observer_nested_draw.path_count = 1;
    observer_nested_draw.realization.mode = GED_SCENE_REALIZE_EAGER;
    observer_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_draw(gedp, &observer_nested_draw,
	    observer_result), observer_result, 1))
	FAIL("observer nested redraw should succeed");
    frontierCapacityRelevant = TRUE;
    if (obs.calls <= 0 || obs.changed <= 0 ||
	visible_compact_occurrence_count(pairSource) != 2 ||
	pairSource->getDisplayMeshLodRevision() !=
	    pairMeshInventoryRevision ||
	controller->renderRequestSerialGet() <= frontierSerial ||
	!controller->consumeRenderRequest(NULL,
	    &frontierCapacityRelevant) || frontierCapacityRelevant)
	FAIL("observer nested redraw should request a presentation-only frame");

    struct ged_scene_erase_request observer_pair_erase;
    ged_scene_erase_request_init(&observer_pair_erase);
    observer_pair_erase.view = ged_view_ctx;
    observer_pair_erase.path = "pair.c";
    observer_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_erase(gedp,
	    &observer_pair_erase, observer_result), observer_result, 1) ||
	render_source_count(controller) != 0)
	FAIL("observer retained-frontier test should cleanly erase its root");
    if (ged_scene_observer_remove(gedp, observerToken) != 1)
	FAIL("GED draw observer should unregister after qtcad Obol sync test");
    controller->clearDatabaseSources();

    const char *box_paths[] = {"box.s"};
    struct ged_scene_draw_request draw_box;
    ged_scene_draw_request_init(&draw_box);
    draw_box.view = ged_view_ctx;
    draw_box.paths = box_paths;
    draw_box.path_count = 1;
    draw_box.realization.mode = GED_SCENE_REALIZE_EAGER;
    draw_box.style.color_override = 1;
    draw_box.style.color[0] = 10;
    draw_box.style.color[1] = 20;
    draw_box.style.color[2] = 30;
    draw_box.style.line_width = 5;
    struct ged_scene_result *scene_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_draw(gedp, &draw_box, scene_result),
	    scene_result, 1))
	FAIL("GED draw should sync a wire Obol database source");
    if (render_source_count(controller) != 1)
	FAIL("Obol draw sync should create one database source");
    SoBRLDatabaseSource *source = render_source(controller, 0);
    if (!source ||
	    !test_path_equal(source->path.getValue().getString(), "box.s") ||
	    source->drawMode.getValue() != SoBRLDatabaseSource::WIREFRAME ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	    !source->hasRealizedWireGeometry())
	FAIL("wire Obol database source should preserve path and realized geometry");
    struct ged_scene_occurrence_info box_record;
    if (!test_shape_record_by_path(gedp, "box.s", &box_record))
	FAIL("GED draw should expose neutral source state for box.s");
    if ((source->visible.getValue() ? 1 : 0) != box_record.visible ||
	    (source->highlighted.getValue() ? 1 : 0) != box_record.highlighted ||
	    source->lineWidth.getValue() != box_record.line_width ||
	    fabs(source->transparency.getValue() - (1.0 - box_record.opacity)) >
	    1.0e-6)
	FAIL("Obol draw sync should seed source display state from semantic occurrence state");
    if (!source->materialColorValid.getValue())
	FAIL("Obol draw sync should publish a valid database material color");
    if (source_for_path(controller, "box.s") != source)
	FAIL("Obol draw sync should publish box.s as the authoritative source owner");

    controller->getViewport()->viewAll();
    controller->requestLodCapacityRender("draw-sync-visible");
    QCoreApplication::processEvents();
    QImage visibleImage;
    view.get_viewport_image(visibleImage);
    if (visibleImage.isNull() || nonblack_pixel_count(visibleImage) < 10)
	FAIL("Obol-synced GED draw should be visible through qtcad capture");

    if (ged_event_notify_object_modified(gedp, "box.s", 1, NULL) != GED_EVENT_OK)
	FAIL("GED stale-source event should refresh Obol source state");
    struct ged_scene_occurrence_info stale_record;
    if (!test_shape_record_by_path(gedp, "box.s", &stale_record))
	FAIL("GED stale-source event should retain the semantic occurrence");
    source = source_for_path(controller, "box.s");
    if (!source)
	FAIL("Obol stale-source sync should retain the box source");
    uint32_t expectedSourceRevision = source->sourceRevision.getValue();
    uint32_t expectedInputsRevision = source->inputsRevision.getValue();
    if (!expectedSourceRevision)
	FAIL("Obol stale-source event should advance a nonzero source revision");
    if (source->sourceRevision.getValue() != expectedSourceRevision ||
	    source->inputsRevision.getValue() != expectedInputsRevision ||
	    source->realizedSourceRevision.getValue() != expectedSourceRevision ||
	    source->realizedInputsRevision.getValue() != expectedInputsRevision ||
	    source->stale.getValue() ||
	    source->staleReason.getValue() != SoBRLDatabaseSource::STALE_NONE)
	FAIL("Obol stale-source sync should realize current GED source revisions");
    BObolRealizedShapeSummary shape_summary;
    if (source->getRealizedShapeSummaryCount() <= 0 ||
	    !source->getRealizedShapeSummary(0, shape_summary) ||
	    shape_summary.ownerSourceRevision != expectedSourceRevision ||
	    shape_summary.ownerInputsRevision != expectedInputsRevision ||
	    shape_summary.ownerSourceStale)
	FAIL("Obol realized summary should carry current GED source lineage");

    struct ged_scene_erase_request erase_box;
    ged_scene_erase_request_init(&erase_box);
    erase_box.view = ged_view_ctx;
    erase_box.path = "box.s";
    scene_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_erase(gedp, &erase_box, scene_result),
	    scene_result, 1))
	FAIL("GED erase should remove an Obol database source");
    if (render_source_count(controller) != 0)
	FAIL("Obol draw sync should remove erased database sources");
    if (scene_display_summary_by_path(controller, "box.s",
	    BObolSceneTreeSummary::NODE_GROUP, NULL))
	FAIL("Obol draw sync should prune empty GED draw groups after erase");

    const char *ball_paths[] = {"ball.s"};
    struct ged_scene_draw_request draw_ball;
    ged_scene_draw_request_init(&draw_ball);
    draw_ball.view = ged_view_ctx;
    draw_ball.paths = ball_paths;
    draw_ball.path_count = 1;
    draw_ball.realization.mode = GED_SCENE_REALIZE_EAGER;
    draw_ball.style.draw_mode = GED_SCENE_DRAW_SHADED;
    scene_result = ged_scene_result_create();
    int drew_ball = scene_result_matches(
	ged_scene_draw(gedp, &draw_ball, scene_result), scene_result, 1);
    if (!drew_ball)
	FAIL("GED shaded draw should sync a shaded Obol database source");
    source = render_source(controller, 0);
    if (!source ||
	    !test_path_equal(source->path.getValue().getString(), "ball.s") ||
	    source->drawMode.getValue() != SoBRLDatabaseSource::SHADED ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
	    !source->hasRealizedMeshGeometry())
	FAIL("shaded Obol database source should preserve draw mode and mesh geometry");

    if (check_terminal_error_hud(view, source))
	return 1;

    const char *paths[2] = {"box.s", "ball.s"};
    struct ged_scene_draw_request draw_both;
    ged_scene_draw_request_init(&draw_both);
    draw_both.view = ged_view_ctx;
    draw_both.paths = paths;
    draw_both.path_count = 2;
    draw_both.realization.mode = GED_SCENE_REALIZE_EAGER;
    draw_both.style.draw_mode = GED_SCENE_DRAW_WIRE;
    draw_both.style.mixed_modes = 1;
    scene_result = ged_scene_result_create();
    int drew_both = scene_result_matches(
	ged_scene_draw(gedp, &draw_both, scene_result), scene_result, 1);
    if (!drew_both)
	FAIL("multi-path GED draw should sync multiple Obol database sources");
    if (render_source_count(controller) != 3 ||
	    !source_for_path_mode(controller, "box.s",
		SoBRLDatabaseSource::WIREFRAME) ||
	    !source_for_path_mode(controller, "ball.s",
		SoBRLDatabaseSource::WIREFRAME) ||
	    !source_for_path_mode(controller, "ball.s",
		SoBRLDatabaseSource::SHADED))
	FAIL("multi-path Obol draw sync should retain shared and representation-specific database sources");

    SoBRLDatabaseSource *boxWireBeforeRedraw =
	source_for_path_mode(controller, "box.s", SoBRLDatabaseSource::WIREFRAME);
    SoBRLDatabaseSource *ballWireBeforeRedraw =
	source_for_path_mode(controller, "ball.s", SoBRLDatabaseSource::WIREFRAME);
    SoBRLDatabaseSource *ballShadedBeforeRedraw =
	source_for_path_mode(controller, "ball.s", SoBRLDatabaseSource::SHADED);
    struct ged_scene_redraw_request redraw_all;
    ged_scene_redraw_request_init(&redraw_all);
    redraw_all.view = ged_view_ctx;
    scene_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_redraw(gedp, &redraw_all,
	    scene_result), scene_result, 1))
	FAIL("GED redraw should invalidate the retained Obol render");
    if (render_source_count(controller) != 3 ||
	    source_for_path_mode(controller, "box.s",
		SoBRLDatabaseSource::WIREFRAME) != boxWireBeforeRedraw ||
	    source_for_path_mode(controller, "ball.s",
		SoBRLDatabaseSource::WIREFRAME) != ballWireBeforeRedraw ||
	    source_for_path_mode(controller, "ball.s",
		SoBRLDatabaseSource::SHADED) != ballShadedBeforeRedraw)
	FAIL("Obol redraw should preserve retained database sources");

    struct ged_scene_clear_request clear_all;
    ged_scene_clear_request_init(&clear_all);
    scene_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_clear(gedp, &clear_all,
	    scene_result), scene_result, 1))
	FAIL("GED clear should clear Obol database sources");
    if (render_source_count(controller) != 0)
	FAIL("Obol draw sync should clear all database sources");

    const char *nested_paths[] = {"pair.c/box.s"};
    struct ged_scene_draw_request draw_nested;
    ged_scene_draw_request_init(&draw_nested);
    draw_nested.view = ged_view_ctx;
    draw_nested.paths = nested_paths;
    draw_nested.path_count = 1;
    draw_nested.realization.mode = GED_SCENE_REALIZE_EAGER;
    scene_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_draw(gedp, &draw_nested,
	    scene_result), scene_result, 1))
	FAIL("nested GED draw should sync an Obol database source");
    if (render_source_count(controller) != 1 ||
	    !source_for_path(controller, "pair.c/box.s"))
	FAIL("nested Obol draw sync should retain one full-path database source");
    source = source_for_path(controller, "pair.c/box.s");
    if (!source ||
	    source->drawMode.getValue() != SoBRLDatabaseSource::WIREFRAME)
	FAIL("full-path Obol draw sync should publish the nested source owner");

    struct ged_scene_erase_request erase_nested;
    ged_scene_erase_request_init(&erase_nested);
    erase_nested.view = ged_view_ctx;
    erase_nested.path = "pair.c/box.s";
    scene_result = ged_scene_result_create();
    if (!scene_result_matches(ged_scene_erase(gedp, &erase_nested,
	    scene_result), scene_result, 1))
	FAIL("nested GED erase should remove an Obol database source");
    if (render_source_count(controller) != 0 ||
	    scene_display_summary_by_path(controller, "pair.c/box.s",
		BObolSceneTreeSummary::NODE_GROUP, NULL) ||
	    scene_display_summary_by_path(controller, "pair.c",
		BObolSceneTreeSummary::NODE_GROUP, NULL))
	FAIL("nested Obol erase should leave no synthetic GED group ancestors");

    ged_close(gedp);
    bu_file_delete(dbpath);
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
