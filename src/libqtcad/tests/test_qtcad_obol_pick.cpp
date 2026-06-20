/*              T E S T _ Q T C A D _ O B O L _ P I C K . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/database_source.h"
#include "brlobol/lod_mesh_shape.h"
#include "brlobol/lod_service.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/pick_detail.h"
#include "brlobol/view_controller.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "ged.h"
#include "ged/bsg_ged_draw.h"
#include "qtcad/QgObolDrawSync.h"
#include "qtcad/QgObolPick.h"
#include "qtcad/QgSelectFilter.h"
#include "qtcad/QgView.h"
#include "raytrace.h"
#include "wdb.h"

#include <Inventor/SoViewport.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoSeparator.h>

#include <QApplication>
#include <QMouseEvent>

#include <chrono>
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
make_pick_db(const char *dbpath)
{
    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp)
	return 0;

    point_t bmin = {-1.0, -1.0, -1.0};
    point_t bmax = { 1.0,  1.0,  1.0};
    int ret = mk_rpp(wdbp, "box.s", bmin, bmax) == 0;
    point_t center = {0.0, 0.0, 0.0};
    ret = ret && mk_sph(wdbp, "implicit.s", center, 0.75) == 0;

    struct wmember region;
    BU_LIST_INIT(&region.l);
    unsigned char color[3] = {32, 96, 192};
    ret = ret &&
	mk_addmember("implicit.s", &region.l, NULL, WMOP_UNION) != NULL &&
	mk_lrcomb(wdbp, "implicit.r", &region, 1, "plastic", "",
	    color, 77, 2, 33, 100, 0) == 0;
    wdb_close(wdbp);
    return ret;
}

static int
apply_and_sync(struct ged *gedp,
	QgView *view,
	struct ged_draw_transaction *txn)
{
    struct ged_draw_transaction_result result;
    ged_draw_transaction_result_init(&result);
    int draw_ret = ged_draw_apply_transaction(gedp, txn, &result);
    int changed = qg_obol_sync_draw_transaction(gedp, txn, &result, view);
    ged_draw_transaction_result_free(&result);

    return draw_ret >= 0 && changed != 0;
}

static BRLObolLodResult
qtcad_source_pick_result(const BRLObolLodRequest &request,
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
    result.geometry.activeLevel = resultKind == BRLOBOL_LOD_RESULT_MESH ? 1 : -1;
    result.geometry.borrowed = FALSE;

    result.mesh.points.push_back(SbVec3f(-1.0f, -1.0f, 0.0f));
    result.mesh.points.push_back(SbVec3f(1.0f, -1.0f, 0.0f));
    result.mesh.points.push_back(SbVec3f(-1.0f, 1.0f, 0.0f));
    result.mesh.coordIndex.push_back(0);
    result.mesh.coordIndex.push_back(1);
    result.mesh.coordIndex.push_back(2);
    for (const SbVec3f &point : result.mesh.points)
	result.bounds.extendBy(point);

    return result;
}

static BRLObolLodResult
qtcad_source_pick_task(const BRLObolLodRequest &request, void *UNUSED(userData))
{
    return qtcad_source_pick_result(request, BRLOBOL_LOD_RESULT_FULL_DETAIL,
	    BRLOBOL_LOD_QUALITY_FULL_DETAIL);
}

static int
wait_for_qtcad_source_pick_result(BRLObolLodService &service)
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

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    QApplication app(argc, argv);

    const char *dbpath = "qtcad_obol_pick_tmp.g";
    if (!make_pick_db(dbpath))
	FAIL("failed to create qtcad Obol pick test database");

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open qtcad Obol pick test database");

    QgView view(NULL, QgView_SW, NULL);
    view.resize(180, 140);
    gedp->ged_gvp = view.view();

    BRLObolViewController *controller = view.obolViewController();
    if (!controller)
	FAIL("QgView should expose an Obol controller");
    controller->clearDatabaseSources();
    controller->setViewportSize(180, 140);

    struct bsg_appearance_settings shaded_settings = BSG_APPEARANCE_SETTINGS_INIT;
    shaded_settings.draw_mode = BSG_DRAW_MODE_SHADED;
    struct ged_draw_transaction draw_box =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "box.s");
    draw_box.view = view.view();
    draw_box.appearance = &shaded_settings;
    if (!apply_and_sync(gedp, &view, &draw_box))
	FAIL("GED shaded draw should sync box source into Obol");

    controller->getViewport()->viewAll();

    std::vector<QgObolPickRecord> picks;
    if (qg_obol_pick_point(&view, 90, 70, 8.0f, false, picks) != 1)
	FAIL("qtcad Obol point pick should find the centered box");
    if (picks[0].path != "/box.s" ||
	    picks[0].primitiveKind != QgObolPickRecord::FACE ||
	    picks[0].primitiveIndex < 0 ||
	    picks[0].sourceName != "box.s" ||
	    picks[0].sourceType != "arb8")
	FAIL("qtcad Obol pick record should preserve BRL-CAD pick detail identity");
    if (view.dmp())
	FAIL("qtcad Obol point pick should not initialize the legacy display manager");

    QgSelectPntFilter filter;
    filter.set_view(view.view());
    filter.set_view_widget(&view);
    QMouseEvent release = mouse_event(QEvent::MouseButtonRelease, 90, 70,
	    Qt::LeftButton, Qt::LeftButton);
    if (!filter.eventFilter(NULL, &release))
	FAIL("qtcad select point filter should accept the release event");
    const std::vector<std::string> &paths = filter.selected_paths();
    if (paths.size() != 1 || paths[0] != "box.s")
	FAIL("qtcad select point filter should expose normalized Obol pick paths");
    if (view.dmp())
	FAIL("qtcad select point filter should not initialize the legacy display manager for Obol picks");

    std::vector<QgObolPickRecord> rectPicks;
    if (qg_obol_pick_rect(&view, 70, 50, 110, 90, 8.0f, false, rectPicks) <= 0)
	FAIL("qtcad Obol rectangle pick should find the centered box");
    bool foundRectPath = false;
    for (const QgObolPickRecord &pick : rectPicks) {
	if (pick.path == "/box.s") {
	    foundRectPath = true;
	    break;
	}
    }
    if (!foundRectPath)
	FAIL("qtcad Obol rectangle pick should preserve BRL-CAD path identity");
    if (view.dmp())
	FAIL("qtcad Obol rectangle pick should not initialize the legacy display manager");

    QgSelectBoxFilter boxFilter;
    boxFilter.set_view(view.view());
    boxFilter.set_view_widget(&view);
    QMouseEvent boxPress = mouse_event(QEvent::MouseButtonPress, 70, 50,
	    Qt::LeftButton, Qt::LeftButton);
    QMouseEvent boxMove = mouse_event(QEvent::MouseMove, 110, 90,
	    Qt::NoButton, Qt::LeftButton);
    QMouseEvent boxRelease = mouse_event(QEvent::MouseButtonRelease, 110, 90,
	    Qt::LeftButton, Qt::LeftButton);
    if (!boxFilter.eventFilter(NULL, &boxPress) ||
	    !boxFilter.eventFilter(NULL, &boxMove) ||
	    !boxFilter.eventFilter(NULL, &boxRelease))
	FAIL("qtcad select box filter should accept the Obol rectangle workflow");
    const std::vector<std::string> &boxPaths = boxFilter.selected_paths();
    if (boxPaths.size() != 1 || boxPaths[0] != "box.s")
	FAIL("qtcad select box filter should expose normalized Obol pick paths");
    if (view.dmp())
	FAIL("qtcad select box filter should not initialize the legacy display manager for Obol picks");

    QgSelectRayFilter rayFilter;
    rayFilter.dbip = gedp->dbip;
    rayFilter.set_view(view.view());
    rayFilter.set_view_widget(&view);
    QMouseEvent rayRelease = mouse_event(QEvent::MouseButtonRelease, 90, 70,
	    Qt::LeftButton, Qt::LeftButton);
    if (!rayFilter.eventFilter(NULL, &rayRelease))
	FAIL("qtcad select ray filter should accept the Obol-backed release event");
    const std::vector<std::string> &rayPaths = rayFilter.selected_paths();
    if (rayPaths.size() != 1 || rayPaths[0] != "box.s")
	FAIL("qtcad select ray filter should expose normalized Obol pick paths");
    if (view.dmp())
	FAIL("qtcad select ray filter should not initialize the legacy display manager for Obol picks");

    SoSeparator *lodRoot = new SoSeparator;
    lodRoot->ref();
    SoBRLDatabaseSource *lodDatabase = new SoBRLDatabaseSource;
    lodDatabase->configureDatabaseSource("lod-pick-db", gedp->dbip,
	    SoBRLDatabaseSource::WIREFRAME, 1);
    SoBRLLodMeshShape *lodMesh = new SoBRLLodMeshShape;
    SbVec3f lodPoints[3] = {
	SbVec3f(-1.0f, -1.0f, 0.0f),
	SbVec3f(1.0f, -1.0f, 0.0f),
	SbVec3f(-1.0f, 1.0f, 0.0f)
    };
    int32_t lodIndices[3] = {0, 1, 2};
    lodMesh->sourcePath = "/lod-pick.bot";
    lodMesh->sourceName = "lod-pick.bot";
    lodMesh->sourceType = "bot";
    lodMesh->sourceId = 9001;
    lodMesh->setIndexedTriangles(lodPoints, 3, lodIndices, 3);
    lodMesh->setPickGeometryPolicy(SoBRLMeshShape::PICK_FULL_DETAIL);

    BRLObolLodRequest displayRequest;
    lodMesh->makeLodRequest(displayRequest,
	    "db://qtcad-obol-pick-test",
	    1,
	    1,
	    1,
	    BRLOBOL_LOD_DRAW_SHADED,
	    "rt_mesh_lod",
	    "rt-cache-v1",
	    BRLOBOL_LOD_QUALITY_FAST_DISPLAY);
    BRLObolLodResult displayResult = qtcad_source_pick_result(displayRequest,
	    BRLOBOL_LOD_RESULT_MESH, BRLOBOL_LOD_QUALITY_FAST_DISPLAY);
    if (!lodMesh->applyStagedLodResult(displayResult, &displayRequest) ||
	    !lodMesh->isLodDisplayActive() ||
	    lodMesh->hasFullDetailMesh())
	FAIL("qtcad LoD pick fixture should have active display LoD without resident full detail");

    lodRoot->addChild(lodDatabase);
    lodRoot->addChild(lodMesh);
    controller->setSceneRoot(lodRoot);
    lodRoot->unref();
    controller->getViewport()->viewAll();

    BRLObolLodRequest sourceLodRequest;
    const SbViewportRegion &pickRegion = controller->getViewportRegion();
    SbVec2s pickSize = pickRegion.getViewportSizePixels();
    SoRayPickAction seededPickAction(pickRegion);
    seededPickAction.setPoint(SbVec2s(90,
	    static_cast<short>(pickSize[1] - 1 - 70)));
    seededPickAction.setRadius(8.0f);
    seededPickAction.setPickAll(FALSE);
    seededPickAction.apply(controller->getViewport()->getRoot());
    SoBRLSourceMeshPickAction seededSourcePickAction;
    seededSourcePickAction.setRay(seededPickAction.getLine().getPosition(),
	    seededPickAction.getLine().getDirection());
    seededSourcePickAction.apply(controller->getViewport()->getRoot());
    if (seededSourcePickAction.getSourceBackedFullDetailRequestCount() != 1 ||
	    !seededSourcePickAction.makeSourceBackedFullDetailLodRequest(0,
		sourceLodRequest))
	FAIL("qtcad LoD pick fixture should build a bounded source full-detail request");

    BRLObolLodService sourceService;
    if (!sourceService.start(1, TRUE))
	FAIL("qtcad LoD pick source service should start");
    controller->setLodService(&sourceService);

    BRLObolLodTask sourceTask;
    sourceTask.request = sourceLodRequest;
    sourceTask.realize = qtcad_source_pick_task;
    if (sourceService.submit(sourceTask) == 0) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad LoD pick source service should accept the ready task");
    }
    if (wait_for_qtcad_source_pick_result(sourceService)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad LoD pick source result should become ready");
    }

    std::vector<QgObolPickRecord> lodPicks;
    if (qg_obol_pick_point(&view, 90, 70, 8.0f, false, lodPicks) != 1) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad Obol point pick should consume ready source-backed full detail");
    }
    if (lodPicks[0].path != "/lod-pick.bot" ||
	    lodPicks[0].sourceName != "lod-pick.bot" ||
	    lodPicks[0].sourceType != "bot" ||
	    lodPicks[0].sourceId != 9001 ||
	    lodPicks[0].primitiveKind != QgObolPickRecord::FACE ||
	    lodPicks[0].primitiveIndex != 0 ||
	    lodPicks[0].distance <= 0.0f) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad Obol source-backed pick should preserve exact mesh identity");
    }
    if (sourceService.queuedResultCountForDiagnostics() != 0) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad Obol source-backed pick should drain only its matching result");
    }

    controller->setExactFullDetailBudget(0, 2);
    std::vector<QgObolPickRecord> overBudgetPicks;
    if (qg_obol_pick_point(&view, 90, 70, 8.0f, false,
	    overBudgetPicks) != 0) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad source-backed exact Obol pick should not use display LoD as over-budget exact geometry");
    }
    if (wait_for_qtcad_source_pick_result(sourceService)) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad source-backed exact Obol pick over-budget result should become ready");
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
	FAIL("qtcad source-backed exact Obol pick should pass controller full-detail budget to source provider");
    }
    controller->setExactFullDetailBudget(0, 0);

    if (view.dmp()) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad source-backed Obol pick should not initialize the legacy display manager");
    }
    controller->setLodService(NULL);
    sourceService.stop();

    SoSeparator *implicitRoot = new SoSeparator;
    implicitRoot->ref();
    SoBRLDatabaseSource *implicitDatabase = new SoBRLDatabaseSource;
    implicitDatabase->configureDatabaseSource("implicit.r", gedp->dbip,
	    SoBRLDatabaseSource::WIREFRAME, 3);
    implicitRoot->addChild(implicitDatabase);
    controller->setSceneRoot(implicitRoot);
    implicitRoot->unref();
    SoCamera *camera = controller->getCamera();
    if (camera)
	camera->position = SbVec3f(0.0f, 0.0f, 5.0f);

    std::vector<QgObolPickRecord> implicitPicks;
    if (qg_obol_pick_point(&view, 90, 70, 8.0f, false,
	    implicitPicks) != 1) {
	FAIL("qtcad Obol point pick should use librt exact implicit pick");
    }
    if (implicitPicks[0].path.find("implicit.r") == std::string::npos ||
	    implicitPicks[0].sourceName != "implicit.s" ||
	    implicitPicks[0].sourceType != "sph" ||
	    implicitPicks[0].primitiveKind != QgObolPickRecord::IMPLICIT_SOLID ||
	    implicitPicks[0].regionId != 77 ||
	    implicitPicks[0].airCode != 2 ||
	    implicitPicks[0].materialId != 33 ||
	    implicitPicks[0].los != 100 ||
	    implicitPicks[0].distance <= 0.0f) {
	FAIL("qtcad Obol librt exact implicit pick should preserve RT hit identity");
    }

    if (view.dmp())
	FAIL("qtcad librt exact implicit pick should not initialize the legacy display manager");

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
