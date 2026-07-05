/*          T E S T _ Q T C A D _ O B O L _ E X P O R T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/database_source.h"
#include "brlobol/export_action.h"
#include "brlobol/lod_mesh_shape.h"
#include "brlobol/lod_service.h"
#include "brlobol/vlist_shape.h"
#include "brlobol/view_controller.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "ged.h"
#include "qtcad/QgObolExport.h"
#include "qtcad/QgView.h"
#include "wdb.h"

#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoSeparator.h>

#include <QApplication>

#include <chrono>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <thread>
#include <vector>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

static int
make_export_db(const char *dbpath)
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
qtcad_export_source_result(const BRLObolLodRequest &request,
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
qtcad_export_source_task(const BRLObolLodRequest &request,
	void *UNUSED(userData))
{
    return qtcad_export_source_result(request,
	    BRLOBOL_LOD_RESULT_FULL_DETAIL,
	    BRLOBOL_LOD_QUALITY_FULL_DETAIL);
}

static int
wait_for_qtcad_export_source_result(BRLObolLodService &service)
{
    for (int i = 0; i < 400; i++) {
	if (service.inFlightCount() == 0 &&
		service.queuedResultCountForDiagnostics() >= 1)
	    return 0;
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    return 1;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);

    QApplication app(argc, argv);

    const char *dbpath = "qtcad_obol_export_tmp.g";
    if (!make_export_db(dbpath))
	FAIL("failed to create qtcad Obol export test database");

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open qtcad Obol export test database");

    QgView view(NULL, QgView_SW);
    view.resize(180, 140);
    qg_legacy_view_ged_active_set(gedp, view.view());

    BRLObolViewController *controller = view.obolViewController();
    if (!controller)
	FAIL("QgView should expose an Obol controller");
    controller->clearDatabaseSources();
    controller->setViewportSize(180, 140);

    SoSeparator *lodRoot = new SoSeparator;
    lodRoot->ref();
    SoBRLDatabaseSource *lodDatabase = new SoBRLDatabaseSource;
    lodDatabase->configureDatabaseSource("lod-export-db", gedp->dbip,
	    SoBRLDatabaseSource::WIREFRAME, 1);
    SoBRLLodMeshShape *lodMesh = new SoBRLLodMeshShape;
    SbVec3f lodPoints[3] = {
	SbVec3f(-1.0f, -1.0f, 0.0f),
	SbVec3f(1.0f, -1.0f, 0.0f),
	SbVec3f(-1.0f, 1.0f, 0.0f)
    };
    int32_t lodIndices[3] = {0, 1, 2};
    lodMesh->sourcePath = "/lod-export.bot";
    lodMesh->sourceName = "lod-export.bot";
    lodMesh->sourceType = "bot";
    lodMesh->sourceId = 9301;
    lodMesh->displayName = "LoD Export Display";
    lodMesh->geometryName = "lod export geometry";
    lodMesh->cacheIdentity = "cache://qtcad/lod-export";
    lodMesh->sourceIdentity = "db://qtcad/lod-export.bot";
    lodMesh->databaseIntent = TRUE;
    lodMesh->sharedSource = FALSE;
    lodMesh->drawMode = BRLOBOL_LOD_DRAW_SHADED;
    lodMesh->recordRole = "database";
    lodMesh->geometryKind = "surface";
    lodMesh->editIntentId = "edit::lod-export/face";
    lodMesh->editIntentRole = "exact-export";
    lodMesh->setIndexedTriangles(lodPoints, 3, lodIndices, 3);

    SoBRLVListShape *viewOverlay = new SoBRLVListShape;
    SbVec3f overlayPoints[3] = {
	SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 0.0f, 0.0f),
	SbVec3f(2.0f, 0.0f, 0.0f)
    };
    int32_t overlayCommands[3] = {
	SoBRLVListShape::MOVE,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::POINT
    };
    viewOverlay->sourcePath = "overlay::qtcad-line-query";
    viewOverlay->sourceName = "qtcad-line-query";
    viewOverlay->sourceType = "overlay";
    viewOverlay->sourceId = 9302;
    viewOverlay->displayName = "Line Query Overlay";
    viewOverlay->geometryName = "line query geometry";
    viewOverlay->cacheIdentity = "cache://qtcad/line-query";
    viewOverlay->sourceIdentity = "view://qtcad/line-query";
    viewOverlay->overlayIntent = TRUE;
    viewOverlay->localSource = TRUE;
    viewOverlay->nonDatabaseSource = TRUE;
    viewOverlay->drawMode = BRLOBOL_LOD_DRAW_WIRE;
    viewOverlay->recordRole = "overlay";
    viewOverlay->visible = FALSE;
    viewOverlay->setLineSet(overlayPoints, overlayCommands, 3);

    BRLObolLodRequest displayRequest;
    lodMesh->makeLodRequest(displayRequest,
	    "db://qtcad-obol-export-test",
	    1,
	    1,
	    1,
	    BRLOBOL_LOD_DRAW_SHADED,
	    "brlobol_mesh_lod",
	    "brlobol-cache-v1",
	    BRLOBOL_LOD_QUALITY_FAST_DISPLAY);
    BRLObolLodResult displayResult = qtcad_export_source_result(
	    displayRequest, BRLOBOL_LOD_RESULT_MESH,
	    BRLOBOL_LOD_QUALITY_FAST_DISPLAY);
    if (!lodMesh->applyStagedLodResult(displayResult, &displayRequest) ||
	    !lodMesh->isLodDisplayActive() ||
	    lodMesh->hasFullDetailMesh())
	FAIL("qtcad LoD export fixture should have active display LoD without resident full detail");

    lodRoot->addChild(lodDatabase);
    lodRoot->addChild(lodMesh);
    lodRoot->addChild(viewOverlay);
    controller->setSceneRoot(lodRoot);
    lodRoot->unref();

    SoBRLExportAction sourceRequestExport;
    sourceRequestExport.setGeometryPolicy(SoBRLExportAction::FULL_DETAIL);
    sourceRequestExport.apply(controller->getViewport()->getRoot());
    BRLObolLodRequest sourceLodRequest;
    if (sourceRequestExport.getSourceBackedFullDetailRequestCount() != 1 ||
	    !sourceRequestExport.makeSourceBackedFullDetailLodRequest(0,
		sourceLodRequest))
	FAIL("qtcad LoD export fixture should build a source full-detail request");

    BRLObolLodService sourceService;
    if (!sourceService.start(1, TRUE))
	FAIL("qtcad LoD export source service should start");
    controller->setLodService(&sourceService);

    BRLObolLodTask sourceTask;
    sourceTask.request = sourceLodRequest;
    sourceTask.realize = qtcad_export_source_task;
    if (sourceService.submit(sourceTask) == 0) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad LoD export source service should accept the ready task");
    }
    if (wait_for_qtcad_export_source_result(sourceService)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad LoD export source result should become ready");
    }

    QgObolExportGeometryRecord exactExport;
    if (!qg_obol_export_geometry_full_detail(&view, exactExport)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol export should consume ready source-backed full detail");
    }
    if (exactExport.lineCount != 0 ||
	    exactExport.pointCount != 0 ||
	    exactExport.triangleCount != 1 ||
	    exactExport.triangles.size() != 1 ||
	    !exactExport.boundsValid) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol export should report source-backed full-detail geometry metrics");
    }
    if (!near_point(exactExport.boundsMin, -1.0f, -1.0f, 0.0f) ||
	    !near_point(exactExport.boundsMax, 1.0f, 1.0f, 0.0f)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol export should preserve source-backed bounds");
    }
    const QgObolExportTriangleRecord &triangle = exactExport.triangles[0];
    const uint64_t cacheIdentity =
	SoBRLExportAction::identityValue("cache://qtcad/lod-export");
    const uint64_t sourceIdentity =
	SoBRLExportAction::identityValue("db://qtcad/lod-export.bot");
    if (triangle.path != "/lod-export.bot" ||
	    triangle.sourceName != "lod-export.bot" ||
	    triangle.sourceType != "bot" ||
	    triangle.displayName != "LoD Export Display" ||
	    triangle.geometryName != "lod export geometry" ||
	    triangle.cacheIdentity != "cache://qtcad/lod-export" ||
	    triangle.sourceIdentity != "db://qtcad/lod-export.bot" ||
	    triangle.cacheIdentityValue != cacheIdentity ||
	    triangle.sourceIdentityValue != sourceIdentity ||
	    triangle.editIntentId != "edit::lod-export/face" ||
	    triangle.editIntentRole != "exact-export" ||
	    triangle.sourceId != 9301 ||
	    !triangle.databaseIntent ||
	    triangle.overlayIntent ||
	    triangle.hudIntent ||
	    triangle.localSource ||
	    triangle.sharedSource ||
	    triangle.nonDatabaseSource ||
	    triangle.drawMode != BRLOBOL_LOD_DRAW_SHADED ||
	    triangle.recordRole != "database" ||
	    triangle.geometryKind != "surface" ||
	    triangle.primitiveIndex != 42 ||
	    triangle.vertexIndexA != 10 ||
	    triangle.vertexIndexB != 11 ||
	    triangle.vertexIndexC != 12 ||
	    !near_point(triangle.a, -1.0f, -1.0f, 0.0f) ||
	    !near_point(triangle.b, 1.0f, -1.0f, 0.0f) ||
	    !near_point(triangle.c, -1.0f, 1.0f, 0.0f)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol export should preserve source-backed triangle identity and edit intent");
    }
    if (sourceService.queuedResultCountForDiagnostics() != 0) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol export should drain only its matching result");
    }
    if (exactExport.sourceFullDetailPending ||
	    exactExport.submittedSourceRequestCount != 0) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol export should not report pending source detail after consuming a ready result");
    }

    QgObolExportObjectQuery objectQuery;
    objectQuery.flags = QG_OBOL_EXPORT_QUERY_VISIBLE_ONLY |
	QG_OBOL_EXPORT_QUERY_DATABASE_OBJECTS;
    objectQuery.glob = "*lod-export.bot";
    objectQuery.drawMode = BRLOBOL_LOD_DRAW_SHADED;
    objectQuery.geometryPolicy = QG_OBOL_EXPORT_DISPLAY_LEVEL;
    QgObolExportObjectQueryResult objectExport;
    if (!qg_obol_export_object_records(&view, objectQuery, objectExport) ||
	    objectExport.objects.size() != 1 ||
	    objectExport.lineCount != 0 ||
	    objectExport.pointCount != 0 ||
	    objectExport.triangleCount != 1 ||
	    objectExport.sourceFullDetailPending ||
	    objectExport.submittedSourceRequestCount != 0) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad Obol object query should report display-level database object records");
    }
    const QgObolExportObjectRecord &object = objectExport.objects[0];
    if (object.path != "/lod-export.bot" ||
	    object.displayName != "LoD Export Display" ||
	    object.geometryName != "lod export geometry" ||
	    object.cacheIdentity != "cache://qtcad/lod-export" ||
	    object.sourceIdentity != "db://qtcad/lod-export.bot" ||
	    object.cacheIdentityValue != cacheIdentity ||
	    object.sourceIdentityValue != sourceIdentity ||
	    !object.databaseIntent ||
	    object.overlayIntent ||
	    object.localSource ||
	    object.sharedSource ||
	    object.nonDatabaseSource ||
	    object.drawMode != BRLOBOL_LOD_DRAW_SHADED ||
	    object.recordRole != "database" ||
	    object.geometryKind != "surface" ||
	    object.lineCount != 0 ||
	    object.pointPrimitiveCount != 0 ||
	    object.triangleCount != 1 ||
	    !object.boundsValid ||
	    !object.surfaceSummary.valid ||
	    object.surfaceSummary.pointCount != 3 ||
	    object.surfaceSummary.indexCount != 3 ||
	    object.surfaceSummary.faceCount != 1 ||
	    object.surfaceSummary.cacheIdentityValue != cacheIdentity ||
	    object.surfaceSummary.sourceIdentityValue != sourceIdentity ||
	    object.surfacePoints.size() != 3 ||
	    object.surfaceIndices.size() != 3 ||
	    object.surfaceIndices[0] != 0 ||
	    object.surfaceIndices[1] != 1 ||
	    object.surfaceIndices[2] != 2 ||
	    !near_point(object.surfacePoints[0], -1.0f, -1.0f, 0.0f) ||
	    !near_point(object.surfacePoints[1], 1.0f, -1.0f, 0.0f) ||
	    !near_point(object.surfacePoints[2], -1.0f, 1.0f, 0.0f)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad Obol object query should preserve neutral object metadata and compact surface detail");
    }
    QgObolExportObjectQuery localOnlyQuery;
    localOnlyQuery.flags = QG_OBOL_EXPORT_QUERY_LOCAL_ONLY;
    localOnlyQuery.geometryPolicy = QG_OBOL_EXPORT_DISPLAY_LEVEL;
    QgObolExportObjectQueryResult localOnlyExport;
    if (qg_obol_export_object_records(&view, localOnlyQuery,
	    localOnlyExport) ||
	    !localOnlyExport.objects.empty()) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad Obol object query should honor local-only filtering");
    }

    viewOverlay->visible = TRUE;
    QgObolExportObjectQuery viewObjectQuery;
    viewObjectQuery.flags = QG_OBOL_EXPORT_QUERY_VISIBLE_ONLY |
	QG_OBOL_EXPORT_QUERY_VIEW_OBJECTS |
	QG_OBOL_EXPORT_QUERY_LOCAL_ONLY;
    viewObjectQuery.glob = "overlay::*";
    viewObjectQuery.drawMode = BRLOBOL_LOD_DRAW_WIRE;
    viewObjectQuery.geometryPolicy = QG_OBOL_EXPORT_DISPLAY_LEVEL;
    QgObolExportObjectQueryResult viewObjectExport;
    const uint64_t overlayCacheIdentity =
	SoBRLExportAction::identityValue("cache://qtcad/line-query");
    const uint64_t overlaySourceIdentity =
	SoBRLExportAction::identityValue("view://qtcad/line-query");
    if (!qg_obol_export_object_records(&view, viewObjectQuery,
	    viewObjectExport) ||
	    viewObjectExport.objects.size() != 1 ||
	    viewObjectExport.lineCount != 1 ||
	    viewObjectExport.pointCount != 1 ||
	    viewObjectExport.triangleCount != 1) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad Obol object query should report view-local line and point records");
    }
    const QgObolExportObjectRecord &viewObject = viewObjectExport.objects[0];
    if (viewObject.path != "overlay::qtcad-line-query" ||
	    viewObject.displayName != "Line Query Overlay" ||
	    viewObject.geometryName != "line query geometry" ||
	    viewObject.cacheIdentity != "cache://qtcad/line-query" ||
	    viewObject.sourceIdentity != "view://qtcad/line-query" ||
	    viewObject.cacheIdentityValue != overlayCacheIdentity ||
	    viewObject.sourceIdentityValue != overlaySourceIdentity ||
	    viewObject.databaseIntent ||
	    !viewObject.overlayIntent ||
	    !viewObject.localSource ||
	    viewObject.sharedSource ||
	    !viewObject.nonDatabaseSource ||
	    viewObject.drawMode != BRLOBOL_LOD_DRAW_WIRE ||
	    viewObject.recordRole != "overlay" ||
	    viewObject.geometryKind != "mixed" ||
	    viewObject.lineCount != 1 ||
	    viewObject.pointPrimitiveCount != 1 ||
	    viewObject.triangleCount != 0 ||
	    !viewObject.lineSummary.valid ||
	    viewObject.lineSummary.pointCount != 2 ||
	    viewObject.lineSummary.segmentCount != 1 ||
	    viewObject.lineSummary.cacheIdentityValue != overlayCacheIdentity ||
	    viewObject.lineSummary.sourceIdentityValue != overlaySourceIdentity ||
	    viewObject.linePoints.size() != 2 ||
	    viewObject.lineCommands.size() != 2 ||
	    viewObject.lineCommands[0] != SoBRLExportAction::LINE_MOVE ||
	    viewObject.lineCommands[1] != SoBRLExportAction::LINE_DRAW ||
	    !near_point(viewObject.linePoints[0], 0.0f, 0.0f, 0.0f) ||
	    !near_point(viewObject.linePoints[1], 1.0f, 0.0f, 0.0f) ||
	    !viewObject.pointSummary.valid ||
	    viewObject.pointSummary.pointCount != 1 ||
	    viewObject.pointSummary.cacheIdentityValue != overlayCacheIdentity ||
	    viewObject.pointSummary.sourceIdentityValue != overlaySourceIdentity ||
	    viewObject.points.size() != 1 ||
	    !near_point(viewObject.points[0], 2.0f, 0.0f, 0.0f) ||
	    viewObject.surfaceSummary.valid ||
	    !viewObject.surfacePoints.empty() ||
	    !viewObject.surfaceIndices.empty()) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad Obol object query should preserve view-local line and point detail");
    }
    viewOverlay->visible = FALSE;

    controller->setExactFullDetailBudget(0, 2);
    QgObolExportGeometryRecord overBudgetExport;
    if (qg_obol_export_geometry_full_detail(&view, overBudgetExport)) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol export should not treat missing over-budget source detail as exported geometry");
    }
    if (!overBudgetExport.sourceFullDetailPending ||
	    overBudgetExport.submittedSourceRequestCount != 1) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol export should report pending submitted source detail");
    }
    if (wait_for_qtcad_export_source_result(sourceService)) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol export over-budget source result should become ready");
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
	FAIL("qtcad exact Obol export should pass controller full-detail budget to source provider");
    }
    controller->setExactFullDetailBudget(0, 0);

    if (view.legacyBackendInitialized()) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad source-backed exact Obol export should not initialize the legacy display manager");
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
