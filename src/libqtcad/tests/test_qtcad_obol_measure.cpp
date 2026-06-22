/*          T E S T _ Q T C A D _ O B O L _ M E A S U R E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/database_source.h"
#include "brlobol/line_layer_overlay.h"
#include "brlobol/lod_mesh_shape.h"
#include "brlobol/lod_service.h"
#include "brlobol/measure_action.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/view_controller.h"
#include "bu/app.h"
#include "bu/color.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/str.h"
#include "ged.h"
#include "qtcad/QgLegacyView.h"
#include "qtcad/QgMeasureFilter.h"
#include "qtcad/QgObolDrawSync.h"
#include "qtcad/QgObolMeasure.h"
#include "qtcad/QgView.h"
#include "wdb.h"

#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoSeparator.h>

#include <QApplication>
#include <QMouseEvent>

#include <chrono>
#include <math.h>
#include <stdio.h>
#include <string>
#include <string.h>
#include <thread>
#include <vector>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

static int
make_measure_db(const char *dbpath)
{
    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp)
	return 0;

    point_t bmin = {-1.0, -1.0, -1.0};
    point_t bmax = { 1.0,  1.0,  1.0};
    int ret = mk_rpp(wdbp, "box.s", bmin, bmax) == 0;
    wdb_close(wdbp);
    return ret;
}

static int
apply_and_sync(struct ged *gedp,
	QgView *view,
	qg_legacy_view_draw_transaction *txn)
{
    qg_legacy_view_draw_transaction_result result;
    qg_legacy_view_draw_result_init(&result);
    int draw_ret = qg_legacy_view_draw_transaction_apply(gedp, txn, &result);
    int changed = qg_obol_sync_draw_transaction(gedp, txn, &result, view);
    qg_legacy_view_draw_result_free(&result);

    return draw_ret >= 0 && changed != 0;
}

static int
nearly_equal(float a, float b)
{
    return fabsf(a - b) < 0.001f;
}

static int
near_point(const SbVec3f &p, float x, float y, float z)
{
    return nearly_equal(p[0], x) && nearly_equal(p[1], y) &&
	nearly_equal(p[2], z);
}

static BRLObolLodResult
qtcad_measure_source_result(const BRLObolLodRequest &request,
	int resultKind,
	int qualityTier)
{
    BRLObolLodResult result;

    result.request = request;
    result.cacheKey = brlobol_lod_cache_key(request);
    result.resultKind = resultKind;
    result.qualityTier = qualityTier;
    result.providerStatus = BRLOBOL_LOD_PROVIDER_READY;
    result.counts.faceCount = 1;
    result.counts.pointCount = 3;
    result.terminal = TRUE;

    result.geometry.kind = BRLOBOL_LOD_GEOMETRY_OBOL_MESH;
    result.geometry.providerId = request.providerId;
    result.geometry.providerVersion = request.providerVersion;
    result.geometry.cacheKey = result.cacheKey;
    result.geometry.activeLevel =
	resultKind == BRLOBOL_LOD_RESULT_MESH ? 1 : -1;
    result.geometry.borrowed = FALSE;

    result.mesh.points.push_back(SbVec3f(-1.0f, -1.0f, 0.0f));
    result.mesh.points.push_back(SbVec3f(1.0f, -1.0f, 0.0f));
    result.mesh.points.push_back(SbVec3f(-1.0f, 1.0f, 0.0f));
    result.mesh.coordIndex.push_back(0);
    result.mesh.coordIndex.push_back(1);
    result.mesh.coordIndex.push_back(2);
    result.mesh.faceIndex.push_back(42);
    result.mesh.vertexIndex.push_back(10);
    result.mesh.vertexIndex.push_back(11);
    result.mesh.vertexIndex.push_back(12);
    for (const SbVec3f &point : result.mesh.points)
	result.bounds.extendBy(point);

    return result;
}

static BRLObolLodResult
qtcad_measure_source_task(const BRLObolLodRequest &request,
	void *UNUSED(userData))
{
    return qtcad_measure_source_result(request,
	    BRLOBOL_LOD_RESULT_FULL_DETAIL,
	    BRLOBOL_LOD_QUALITY_FULL_DETAIL);
}

static int
wait_for_qtcad_measure_source_result(BRLObolLodService &service)
{
    for (int i = 0; i < 400; i++) {
	if (service.inFlightCount() == 0 &&
		service.queuedResultCountForDiagnostics() >= 1)
	    return 0;
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    return 1;
}

static QMouseEvent
mouse_event(QEvent::Type type,
	int x,
	int y,
	Qt::MouseButton button,
	Qt::MouseButtons buttons)
{
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    return QMouseEvent(type, QPoint(x, y), button, buttons, Qt::NoModifier);
#else
    return QMouseEvent(type, QPointF(x, y), QPointF(x, y), button, buttons,
	    Qt::NoModifier);
#endif
}

static SoBRLLineLayerOverlay *
find_measure_overlay(BRLObolViewController *controller,
	const char *overlayId)
{
    if (!controller || !controller->getSceneRoot() ||
	    !controller->getSceneRoot()->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    SoGroup *root = static_cast<SoGroup *>(controller->getSceneRoot());
    for (int i = 0; i < root->getNumChildren(); i++) {
	SoNode *node = root->getChild(i);
	if (!node || !node->isOfType(SoBRLLineLayerOverlay::getClassTypeId()))
	    continue;
	SoBRLLineLayerOverlay *overlay =
	    static_cast<SoBRLLineLayerOverlay *>(node);
	if (BU_STR_EQUAL(overlay->overlayId.getValue().getString(), overlayId))
	    return overlay;
    }
    return NULL;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    QApplication app(argc, argv);

    const char *dbpath = "qtcad_obol_measure_tmp.g";
    if (!make_measure_db(dbpath))
	FAIL("failed to create qtcad Obol measure test database");

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open qtcad Obol measure test database");

    QgView view(NULL, QgView_SW);
    view.resize(180, 140);
    qg_legacy_view_ged_active_set(gedp, view.view());

    BRLObolViewController *controller = view.obolViewController();
    if (!controller)
	FAIL("QgView should expose an Obol controller");
    controller->clearDatabaseSources();
    controller->setViewportSize(180, 140);

    qg_legacy_view_draw_appearance *shaded_appearance =
	qg_legacy_view_draw_appearance_create(QG_LEGACY_VIEW_DRAW_MODE_SHADED);
    qg_legacy_view_draw_transaction draw_box;
    qg_legacy_view_draw_transaction_init(&draw_box,
	    QG_LEGACY_VIEW_DRAW_TXN_DRAW, "box.s");
    qg_legacy_view_draw_transaction_view_set(&draw_box, view.view());
    qg_legacy_view_draw_transaction_appearance_set(&draw_box,
	    shaded_appearance);
    int drew_box = apply_and_sync(gedp, &view, &draw_box);
    qg_legacy_view_draw_appearance_destroy(shaded_appearance);
    if (!drew_box)
	FAIL("GED shaded draw should sync box source into Obol");

    controller->getViewport()->viewAll();

    SbVec3f pickedPoint;
    std::string pickedPath;
    if (!qg_obol_measure_pick_point(&view, 90, 70, pickedPoint,
	    &pickedPath))
	FAIL("qtcad Obol measure helper should pick the centered box");
    if (pickedPath != "/box.s")
	FAIL("qtcad Obol measure helper should preserve picked path identity");
    if (view.legacyBackendInitialized())
	FAIL("qtcad Obol measure pick should not initialize the legacy display manager");

    SbVec3f guidePoints[2] = {
	pickedPoint,
	SbVec3f(pickedPoint[0] + 0.5f, pickedPoint[1], pickedPoint[2])
    };
    struct bu_color green = BU_COLOR_GREEN;
    if (!qg_obol_measure_update_overlay(&view, "tool:measurement",
	    guidePoints, 2, &green))
	FAIL("qtcad Obol measure helper should publish an Obol overlay");
    if (view.legacyBackendInitialized())
	FAIL("qtcad Obol measure overlay should not initialize the legacy display manager");

    SoBRLLineLayerOverlay *overlay =
	find_measure_overlay(controller, "tool:measurement");
    if (!overlay || overlay->getLayerShapeCount() != 1 ||
	    overlay->getPointCount() != 2 ||
	    overlay->selectable.getValue())
	FAIL("qtcad Obol measure overlay should be a non-selectable line layer");

    if (!qg_obol_measure_clear_overlay(&view, "tool:measurement") ||
	    find_measure_overlay(controller, "tool:measurement"))
	FAIL("qtcad Obol measure helper should clear its overlay");
    if (view.legacyBackendInitialized())
	FAIL("qtcad Obol measure clear should not initialize the legacy display manager");

    QMeasure3DFilter filter;
    filter.set_view_widget(&view);
    filter.update_color(&green);

    QMouseEvent press = mouse_event(QEvent::MouseButtonPress, 90, 70,
	    Qt::LeftButton, Qt::LeftButton);
    if (!filter.eventFilter(NULL, &press) || filter.mode != 1)
	FAIL("qtcad 3D measure filter should start with an Obol hit point");
    overlay = find_measure_overlay(controller, "tool:measurement");
    if (!overlay || overlay->getPointCount() != 1)
	FAIL("qtcad 3D measure filter should start with a single-point Obol overlay");
    if (view.legacyBackendInitialized())
	FAIL("qtcad 3D measure filter start should not initialize the legacy display manager");

    QMouseEvent move = mouse_event(QEvent::MouseMove, 115, 70,
	    Qt::NoButton, Qt::LeftButton);
    if (!filter.eventFilter(NULL, &move))
	FAIL("qtcad 3D measure filter should update the Obol guide line");
    if (view.legacyBackendInitialized())
	FAIL("qtcad 3D measure filter update should not initialize the legacy display manager");

    QMouseEvent release = mouse_event(QEvent::MouseButtonRelease, 115, 70,
	    Qt::LeftButton, Qt::LeftButton);
    if (!filter.eventFilter(NULL, &release) || filter.mode != 2)
	FAIL("qtcad 3D measure filter should finalize the first Obol measure segment");
    if (filter.length1() <= 0.0)
	FAIL("qtcad 3D Obol measure segment should report a positive length");

    overlay = find_measure_overlay(controller, "tool:measurement");
    if (!overlay || overlay->getPointCount() != 2)
	FAIL("qtcad 3D measure filter should publish the finalized segment as an Obol overlay");
    if (view.legacyBackendInitialized())
	FAIL("qtcad 3D measure filter finalize should not initialize the legacy display manager");

    QMouseEvent rightPress = mouse_event(QEvent::MouseButtonPress, 90, 70,
	    Qt::RightButton, Qt::RightButton);
    if (!filter.eventFilter(NULL, &rightPress) ||
	    find_measure_overlay(controller, "tool:measurement"))
	FAIL("qtcad 3D measure filter should clear the Obol overlay on right click");
    if (view.legacyBackendInitialized())
	FAIL("qtcad 3D measure filter clear should not initialize the legacy display manager");

    SoSeparator *lodRoot = new SoSeparator;
    lodRoot->ref();
    SoBRLDatabaseSource *lodDatabase = new SoBRLDatabaseSource;
    lodDatabase->configureDatabaseSource("lod-measure-db", gedp->dbip,
	    SoBRLDatabaseSource::WIREFRAME, 1);
    SoBRLLodMeshShape *lodMesh = new SoBRLLodMeshShape;
    SbVec3f lodPoints[3] = {
	SbVec3f(-1.0f, -1.0f, 0.0f),
	SbVec3f(1.0f, -1.0f, 0.0f),
	SbVec3f(-1.0f, 1.0f, 0.0f)
    };
    int32_t lodIndices[3] = {0, 1, 2};
    lodMesh->sourcePath = "/lod-measure.bot";
    lodMesh->sourceName = "lod-measure.bot";
    lodMesh->sourceType = "bot";
    lodMesh->sourceId = 9201;
    lodMesh->editIntentId = "edit::lod-measure/face";
    lodMesh->editIntentRole = "exact-measure";
    lodMesh->setIndexedTriangles(lodPoints, 3, lodIndices, 3);

    BRLObolLodRequest displayRequest;
    lodMesh->makeLodRequest(displayRequest,
	    "db://qtcad-obol-measure-test",
	    1,
	    1,
	    1,
	    BRLOBOL_LOD_DRAW_SHADED,
	    "rt_mesh_lod",
	    "rt-cache-v1",
	    BRLOBOL_LOD_QUALITY_FAST_DISPLAY);
    BRLObolLodResult displayResult = qtcad_measure_source_result(
	    displayRequest, BRLOBOL_LOD_RESULT_MESH,
	    BRLOBOL_LOD_QUALITY_FAST_DISPLAY);
    if (!lodMesh->applyStagedLodResult(displayResult, &displayRequest) ||
	    !lodMesh->isLodDisplayActive() ||
	    lodMesh->hasFullDetailMesh())
	FAIL("qtcad LoD measure fixture should have active display LoD without resident full detail");

    lodRoot->addChild(lodDatabase);
    lodRoot->addChild(lodMesh);
    controller->setSceneRoot(lodRoot);
    lodRoot->unref();

    SbVec3f query(-0.2f, -0.2f, 0.02f);
    SoBRLMeasureAction sourceRequestMeasure;
    sourceRequestMeasure.setGeometryPolicy(SoBRLMeasureAction::FULL_DETAIL);
    sourceRequestMeasure.setQueryPoint(query);
    sourceRequestMeasure.apply(lodRoot);
    BRLObolLodRequest sourceLodRequest;
    if (sourceRequestMeasure.getSourceBackedFullDetailRequestCount() != 1 ||
	    !sourceRequestMeasure.makeSourceBackedFullDetailLodRequest(0,
		sourceLodRequest))
	FAIL("qtcad LoD measure fixture should build a query-scoped source full-detail request");

    BRLObolLodService sourceService;
    if (!sourceService.start(1, TRUE))
	FAIL("qtcad LoD measure source service should start");
    controller->setLodService(&sourceService);

    BRLObolLodTask sourceTask;
    sourceTask.request = sourceLodRequest;
    sourceTask.realize = qtcad_measure_source_task;
    if (sourceService.submit(sourceTask) == 0) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad LoD measure source service should accept the ready task");
    }
    if (wait_for_qtcad_measure_source_result(sourceService)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad LoD measure source result should become ready");
    }

    QgObolMeasureGeometryRecord exactMeasure;
    if (!qg_obol_measure_geometry_full_detail(&view, &query, exactMeasure)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol measure should consume ready source-backed full detail");
    }
    if (exactMeasure.shapeCount != 1 ||
	    exactMeasure.triangleCount != 1 ||
	    exactMeasure.segmentCount != 0 ||
	    exactMeasure.surfaceArea < 1.9f ||
	    exactMeasure.surfaceArea > 2.1f ||
	    !exactMeasure.boundsValid) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol measure should report source-backed full-detail geometry metrics");
    }
    if (!near_point(exactMeasure.boundsMin, -1.0f, -1.0f, 0.0f) ||
	    !near_point(exactMeasure.boundsMax, 1.0f, 1.0f, 0.0f)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol measure should preserve source-backed bounds");
    }
    if (!exactMeasure.hasNearestPrimitive ||
	    exactMeasure.nearestPrimitiveKind != QgObolMeasureGeometryRecord::FACE ||
	    exactMeasure.nearestPrimitiveIndex != 42 ||
	    exactMeasure.nearestFaceVertexIndexA != 10 ||
	    exactMeasure.nearestFaceVertexIndexB != 11 ||
	    exactMeasure.nearestFaceVertexIndexC != 12 ||
	    exactMeasure.nearestFaceEdgeSlot != 1 ||
	    exactMeasure.nearestFaceEdgeVertexIndexA != 11 ||
	    exactMeasure.nearestFaceEdgeVertexIndexB != 12 ||
	    exactMeasure.nearestFaceVertexSlot != 0 ||
	    exactMeasure.nearestFaceVertexIndex != 10 ||
	    exactMeasure.nearestEditIntentId != "edit::lod-measure/face" ||
	    exactMeasure.nearestEditIntentRole != "exact-measure" ||
	    exactMeasure.nearestPath != "/lod-measure.bot" ||
	    !near_point(exactMeasure.nearestPoint, -0.2f, -0.2f, 0.0f)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol measure should preserve source-backed face, sub-entity identity, and edit intent");
    }
    if (sourceService.queuedResultCountForDiagnostics() != 0) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol measure should drain only its matching result");
    }
    if (exactMeasure.sourceFullDetailPending ||
	    exactMeasure.submittedSourceRequestCount != 0) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol measure should not report pending source detail after consuming a ready result");
    }

    controller->setExactFullDetailBudget(0, 2);
    QgObolMeasureGeometryRecord overBudgetMeasure;
    if (qg_obol_measure_geometry_full_detail(&view, &query,
	    overBudgetMeasure)) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol measure should not treat missing over-budget source detail as measured geometry");
    }
    if (!overBudgetMeasure.sourceFullDetailPending ||
	    overBudgetMeasure.submittedSourceRequestCount != 1) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol measure should report pending submitted source detail");
    }
    if (wait_for_qtcad_measure_source_result(sourceService)) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol measure over-budget source result should become ready");
    }
    std::vector<BRLObolLodResult> overBudgetResults;
    sourceService.drainResults(overBudgetResults);
    if (overBudgetResults.size() != 1 ||
	    overBudgetResults[0].providerStatus != BRLOBOL_LOD_PROVIDER_FALLBACK ||
	    strcmp(overBudgetResults[0].diagnostic.getString(),
		"RT source full-detail provider request exceeds full-detail limits") != 0 ||
	    overBudgetResults[0].mesh.isValid()) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol measure should pass controller full-detail budget to source provider");
    }
    controller->setExactFullDetailBudget(0, 0);

    if (view.legacyBackendInitialized()) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad source-backed exact Obol measure should not initialize the legacy display manager");
    }
    controller->setLodService(NULL);
    sourceService.stop();

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
