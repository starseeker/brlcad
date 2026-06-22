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
    lodMesh->editIntentId = "edit::lod-export/face";
    lodMesh->editIntentRole = "exact-export";
    lodMesh->setIndexedTriangles(lodPoints, 3, lodIndices, 3);

    BRLObolLodRequest displayRequest;
    lodMesh->makeLodRequest(displayRequest,
	    "db://qtcad-obol-export-test",
	    1,
	    1,
	    1,
	    BRLOBOL_LOD_DRAW_SHADED,
	    "rt_mesh_lod",
	    "rt-cache-v1",
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
    if (triangle.path != "/lod-export.bot" ||
	    triangle.sourceName != "lod-export.bot" ||
	    triangle.sourceType != "bot" ||
	    triangle.editIntentId != "edit::lod-export/face" ||
	    triangle.editIntentRole != "exact-export" ||
	    triangle.sourceId != 9301 ||
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
