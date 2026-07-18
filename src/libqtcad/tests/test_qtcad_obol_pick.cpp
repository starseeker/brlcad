/*              T E S T _ Q T C A D _ O B O L _ P I C K . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BLodMeshShape.h"
#include "BObol/BLodService.h"
#include "BObol/BMeshShape.h"
#include "BObol/BPickDetail.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "ged.h"
#include "ged/draw.h"
#include "QgObolDrawSyncPrivate.h"
#include "qtcad/QgObolPick.h"
#include "qtcad/QgSelectFilter.h"
#include "qtcad/QgView.h"
#include "raytrace.h"
#include "wdb.h"

#include <Inventor/SoViewport.h>
#include <Inventor/SbRotation.h>
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
    point_t far_bmin = {-1.0, -1.0, -3.0};
    point_t far_bmax = { 1.0,  1.0, -2.0};
    ret = ret && mk_rpp(wdbp, "farbox.s", far_bmin, far_bmax) == 0;
    point_t center = {0.0, 0.0, 0.0};
    ret = ret && mk_sph(wdbp, "implicit.s", center, 0.75) == 0;

    struct wmember region;
    BU_LIST_INIT(&region.l);
    unsigned char color[3] = {32, 96, 192};
    ret = ret &&
	mk_addmember("implicit.s", &region.l, NULL, WMOP_UNION) != NULL &&
	mk_lrcomb(wdbp, "implicit.r", &region, 1, "plastic", "",
	    color, 77, 2, 33, 100, 0) == 0;

    fastf_t lod_vertices[9] = {
	-1.0, -1.0, 0.0,
	 1.0, -1.0, 0.0,
	-1.0,  1.0, 0.0
    };
    int lod_faces[3] = {0, 1, 2};
    ret = ret && mk_bot(wdbp, "lod-pick.bot", RT_BOT_SOLID,
	    RT_BOT_UNORIENTED, 0, 3, 1, lod_vertices, lod_faces,
	    NULL, NULL) == 0;

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
    int changed = qg_obol_sync_ged_draw_transaction(gedp, txn, &result, view);
    ged_draw_transaction_result_free(&result);

    return draw_ret >= 0 && changed != 0;
}

static BObolLodResult
qtcad_source_pick_result(const BObolLodRequest &request,
	int resultKind,
	int qualityTier)
{
    BObolLodResult result;

    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = resultKind;
    result.qualityTier = qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.counts.faceCount = 1;
    result.counts.pointCount = 3;
    result.terminal = TRUE;

    result.geometry.kind = BOBOL_LOD_GEOMETRY_OBOL_MESH;
    result.geometry.providerId = request.providerId;
    result.geometry.providerVersion = request.providerVersion;
    result.geometry.cacheKey = result.cacheKey;
    result.geometry.activeLevel = resultKind == BOBOL_LOD_RESULT_MESH ? 1 : -1;
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

static BObolLodResult
qtcad_source_pick_task(const BObolLodRequest &request, void *UNUSED(userData))
{
    return qtcad_source_pick_result(request, BOBOL_LOD_RESULT_FULL_DETAIL,
	    BOBOL_LOD_QUALITY_FULL_DETAIL);
}

static int
wait_for_qtcad_source_pick_result(BObolLodService &service)
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
make_mouse_event(QEvent::Type type,
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

    QgView view(NULL, QgViewType::SW);
    view.resize(180, 140);
    struct ged_view_context *view_ctx =
	ged_view_context_from_bv(view.viewContext());
    ged_view_active_ctx_set(gedp, view_ctx);
    (void)ged_view_context_host_attach(gedp, view_ctx);

    BObolViewController *controller = view.obolViewController();
    if (!controller)
	FAIL("QgView should expose an Obol controller");
    controller->clearDatabaseSources();
    controller->setViewportSize(180, 140);

    struct ged_draw_appearance_settings shaded_appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    shaded_appearance.draw_mode = GED_DRAW_MODE_SHADED;
    struct ged_draw_transaction draw_box =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "box.s");
    draw_box.view = view_ctx;
    draw_box.appearance = &shaded_appearance;
    int drew_box = apply_and_sync(gedp, &view, &draw_box);
    if (!drew_box)
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

    QgSelectPntFilter filter;
    filter.set_view_widget(&view);
    QMouseEvent release = make_mouse_event(QEvent::MouseButtonRelease, 90, 70,
	    Qt::LeftButton, Qt::LeftButton);
    if (!filter.eventFilter(NULL, &release))
	FAIL("qtcad select point filter should accept the release event");
    const std::vector<std::string> &paths = filter.selected_paths();
    if (paths.size() != 1 || paths[0] != "box.s")
	FAIL("qtcad select point filter should expose normalized Obol pick paths");

    BObolFeatureStyle command_style;
    command_style.hasSelectable = TRUE;
    command_style.selectable = TRUE;
    command_style.hasColor = TRUE;
    command_style.color = SbColor(1.0f, 1.0f, 0.0f);
    BObolLineLayer command_layer;
    command_layer.name = "rtcheck::overlaps/yellow";
    command_layer.points.push_back(SbVec3f(-0.5f, 0.0f, 2.0f));
    command_layer.commands.push_back(static_cast<int32_t>(
	    BObolLineCommand::Move));
    command_layer.points.push_back(SbVec3f(0.5f, 0.0f, 2.0f));
    command_layer.commands.push_back(static_cast<int32_t>(
	    BObolLineCommand::Draw));
    std::vector<BObolLineLayer> command_layers;
    command_layers.push_back(command_layer);
    BObolFeatureHandle command_handle =
	controller->features().publishLineLayers("rtcheck::overlaps",
		BObolFeatureScope::Shared, command_layers, &command_style);
    std::vector<BObolFeatureMetadata> command_feature_metadata;
    BObolFeatureMetadata command_schema_metadata;
    command_schema_metadata.key = "result.schema";
    command_schema_metadata.value = "brlcad.rtcheck.overlap.v1";
    command_feature_metadata.push_back(command_schema_metadata);
    std::vector<BObolFeatureMetadata> command_primitive_metadata;
    BObolFeatureMetadata command_metadata;
    command_metadata.key = "overlap.objects";
    command_metadata.value = "box.s cone.s";
    command_primitive_metadata.push_back(command_metadata);
	if (!command_handle.isValid() ||
	    !controller->features().replaceMetadata(command_handle,
		command_feature_metadata) ||
	    !controller->features().replacePrimitiveMetadata(command_handle,
		0, command_primitive_metadata))
	FAIL("qtcad command-result pick fixture should publish primitive metadata");
    SoCamera *command_camera = controller->getCamera();
    if (command_camera) {
	command_camera->position = SbVec3f(0.0f, 0.0f, 5.0f);
	command_camera->orientation = SbRotation::identity();
    }

    std::vector<QgObolPickRecord> commandRayPicks;
    if (qg_obol_pick_ray(&view, SbVec3f(0.0f, 0.0f, 5.0f),
	    SbVec3f(0.0f, 0.0f, -1.0f), false, commandRayPicks) != 1)
	FAIL("qtcad command-result ray pick should hit the front line-layer child");
    if (commandRayPicks[0].path != "rtcheck::overlaps/yellow" ||
	    !commandRayPicks[0].featurePickResolved ||
	    commandRayPicks[0].featureName != "rtcheck::overlaps" ||
	    commandRayPicks[0].featurePrimitiveIndex != 0 ||
	    commandRayPicks[0].featureMetadata.size() != 1 ||
	    commandRayPicks[0].featureMetadata[0].first != "result.schema" ||
	    commandRayPicks[0].featureMetadata[0].second !=
		"brlcad.rtcheck.overlap.v1" ||
	    commandRayPicks[0].featurePrimitiveMetadata.size() != 1 ||
	    commandRayPicks[0].featurePrimitiveMetadata[0].first !=
		"overlap.objects" ||
	    commandRayPicks[0].featurePrimitiveMetadata[0].second !=
		"box.s cone.s")
	FAIL("qtcad command-result pick should resolve child layer hits to parent primitive metadata");
    if (qg_obol_pick_apply_feature_state(&view, commandRayPicks[0], true,
	    true) != 1)
	FAIL("qtcad command-result pick helper should apply parent primitive state");
    BObolFeatureRecord command_record;
    if (!controller->features().record(command_handle, command_record) ||
	    command_record.selectedPrimitives.size() != 1 ||
	    command_record.selectedPrimitives[0] != 0 ||
	    command_record.highlightedPrimitives.size() != 1 ||
	    command_record.highlightedPrimitives[0] != 0)
	FAIL("qtcad command-result pick helper should update parent primitive state");
    std::vector<int32_t> no_primitives;
    if (!controller->features().replaceSelectedPrimitives(command_handle,
	    no_primitives) ||
	    !controller->features().replaceHighlightedPrimitives(command_handle,
		no_primitives))
	FAIL("qtcad command-result pick fixture should clear primitive state");

    QgSelectPntFilter commandFilter;
    commandFilter.first_only = true;
    commandFilter.set_view_widget(&view);
    QMouseEvent commandRelease = make_mouse_event(QEvent::MouseButtonRelease, 90,
	    70, Qt::LeftButton, Qt::LeftButton);
    if (!commandFilter.eventFilter(NULL, &commandRelease))
	FAIL("qtcad select point filter should accept command-result picks");
    if (!commandFilter.selected_paths().empty())
	FAIL("qtcad command-result selection should not publish database paths");
    if (!controller->features().record(command_handle, command_record) ||
	    command_record.selectedPrimitives.size() != 1 ||
	    command_record.selectedPrimitives[0] != 0 ||
	    command_record.highlightedPrimitives.size() != 1 ||
	    command_record.highlightedPrimitives[0] != 0)
	FAIL("qtcad select point filter should apply command-result primitive state");
    controller->features().remove(command_handle);

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

    QgSelectBoxFilter boxFilter;
    boxFilter.set_view_widget(&view);
    QMouseEvent boxPress = make_mouse_event(QEvent::MouseButtonPress, 70, 50,
	    Qt::LeftButton, Qt::LeftButton);
    QMouseEvent boxMove = make_mouse_event(QEvent::MouseMove, 110, 90,
	    Qt::NoButton, Qt::LeftButton);
    QMouseEvent boxRelease = make_mouse_event(QEvent::MouseButtonRelease, 110, 90,
	    Qt::LeftButton, Qt::LeftButton);
    if (!boxFilter.eventFilter(NULL, &boxPress) ||
	    !boxFilter.eventFilter(NULL, &boxMove) ||
	    !boxFilter.eventFilter(NULL, &boxRelease))
	FAIL("qtcad select box filter should accept the Obol rectangle workflow");
    const std::vector<std::string> &boxPaths = boxFilter.selected_paths();
    if (boxPaths.size() != 1 || boxPaths[0] != "box.s")
	FAIL("qtcad select box filter should expose normalized Obol pick paths");

    QgSelectRayFilter rayFilter;
    rayFilter.dbip = gedp->dbip;
    rayFilter.set_view_widget(&view);
    QMouseEvent rayRelease = make_mouse_event(QEvent::MouseButtonRelease, 90, 70,
	    Qt::LeftButton, Qt::LeftButton);
    if (!rayFilter.eventFilter(NULL, &rayRelease))
	FAIL("qtcad select ray filter should accept the Obol-backed release event");
    const std::vector<std::string> &rayPaths = rayFilter.selected_paths();
    if (rayPaths.size() != 1 || rayPaths[0] != "box.s")
	FAIL("qtcad select ray filter should expose normalized Obol pick paths");

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
    lodMesh->editIntentId = "edit::lod-pick/face";
    lodMesh->editIntentRole = "exact-pick";
    lodMesh->setIndexedTriangles(lodPoints, 3, lodIndices, 3);
    lodMesh->setPickGeometryPolicy(SoBRLMeshShape::PICK_FULL_DETAIL);

    BObolLodRequest displayRequest;
    lodMesh->makeLodRequest(displayRequest,
	    "db://qtcad-obol-pick-test",
	    1,
	    1,
	    1,
	    BOBOL_LOD_DRAW_SHADED,
	    "bobol_mesh_lod",
	    "bobol-cache-v1",
	    BOBOL_LOD_QUALITY_FAST_DISPLAY);
    BObolLodResult displayResult = qtcad_source_pick_result(displayRequest,
	    BOBOL_LOD_RESULT_MESH, BOBOL_LOD_QUALITY_FAST_DISPLAY);
    if (!lodMesh->applyStagedLodResult(displayResult, &displayRequest) ||
	    !lodMesh->isLodDisplayActive() ||
	    lodMesh->hasFullDetailMesh())
	FAIL("qtcad LoD pick fixture should have active display LoD without resident full detail");

    lodRoot->addChild(lodDatabase);
    lodRoot->addChild(lodMesh);
    controller->setSceneRoot(lodRoot);
    lodRoot->unref();
    controller->getViewport()->viewAll();

    BObolLodRequest sourceLodRequest;
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

    BObolLodService sourceService;
    if (!sourceService.start(1, TRUE))
	FAIL("qtcad LoD pick source service should start");
    controller->setLodService(&sourceService);

    BObolLodTask sourceTask;
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
	    lodPicks[0].editIntentId != "edit::lod-pick/face" ||
	    lodPicks[0].editIntentRole != "exact-pick" ||
	    lodPicks[0].sourceId != 9001 ||
	    lodPicks[0].primitiveKind != QgObolPickRecord::FACE ||
	    lodPicks[0].primitiveIndex != 42 ||
	    lodPicks[0].faceVertexIndexA != 10 ||
	    lodPicks[0].faceVertexIndexB != 11 ||
	    lodPicks[0].faceVertexIndexC != 12 ||
	    lodPicks[0].nearestFaceEdgeSlot != 1 ||
	    lodPicks[0].nearestFaceEdgeVertexIndexA != 11 ||
	    lodPicks[0].nearestFaceEdgeVertexIndexB != 12 ||
	    lodPicks[0].nearestFaceVertexSlot != 0 ||
	    lodPicks[0].nearestFaceVertexIndex != 10 ||
	    lodPicks[0].distance <= 0.0f) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad Obol source-backed pick should preserve exact mesh, sub-entity identity, and edit intent");
    }
    if (sourceService.queuedResultCountForDiagnostics() != 0) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad Obol source-backed pick should drain only its matching result");
    }

    controller->setExactFullDetailBudget(0, 2);
    std::vector<QgObolPickRecord> overBudgetPicks;
    int overBudgetSubmittedCount = -1;
    if (qg_obol_pick_point(&view, 90, 70, 8.0f, false,
	    overBudgetPicks, &overBudgetSubmittedCount) != 0) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad source-backed exact Obol pick should not use display LoD as over-budget exact geometry");
    }
    if (overBudgetSubmittedCount != 1) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad source-backed exact Obol pick should report pending submitted source detail");
    }
    if (wait_for_qtcad_source_pick_result(sourceService)) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad source-backed exact Obol pick over-budget result should become ready");
    }
    std::vector<BObolLodResult> overBudgetResults;
    sourceService.drainResults(overBudgetResults);
    if (overBudgetResults.size() != 1 ||
	    overBudgetResults[0].providerStatus != BOBOL_LOD_PROVIDER_FALLBACK ||
	    bu_strcmp(overBudgetResults[0].diagnostic.getString(),
		"RT source full-detail provider request exceeds full-detail limits") != 0 ||
	    overBudgetResults[0].mesh.isValid()) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad source-backed exact Obol pick should pass controller full-detail budget to source provider");
    }

    QgSelectPntFilter pendingSourcePickFilter;
    pendingSourcePickFilter.set_view_widget(&view);
    QMouseEvent pendingSourcePickRelease = make_mouse_event(
	    QEvent::MouseButtonRelease, 90, 70, Qt::LeftButton,
	    Qt::LeftButton);
    if (!pendingSourcePickFilter.eventFilter(NULL, &pendingSourcePickRelease)) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad point select filter should consume pending source-backed exact pick");
    }
    if (!pendingSourcePickFilter.selected_paths().empty()) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad point select filter should not fall back to legacy selection while exact source pick is pending");
    }
    if (wait_for_qtcad_source_pick_result(sourceService)) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad point select filter pending source pick result should become ready");
    }
    {
	std::vector<BObolLodResult> pendingFilterResults;
	sourceService.drainResults(pendingFilterResults);
	if (pendingFilterResults.size() != 1 ||
		pendingFilterResults[0].providerStatus !=
		    BOBOL_LOD_PROVIDER_FALLBACK) {
	    controller->setExactFullDetailBudget(0, 0);
	    controller->setLodService(NULL);
	    sourceService.stop();
	    FAIL("qtcad point select filter should submit the pending source-backed exact pick request");
	}
    }

    QgSelectBoxFilter pendingSourceBoxFilter;
    pendingSourceBoxFilter.set_view_widget(&view);
    QMouseEvent pendingSourceBoxPress = make_mouse_event(
	    QEvent::MouseButtonPress, 70, 50, Qt::LeftButton,
	    Qt::LeftButton);
    QMouseEvent pendingSourceBoxMove = make_mouse_event(
	    QEvent::MouseMove, 110, 90, Qt::NoButton,
	    Qt::LeftButton);
    QMouseEvent pendingSourceBoxRelease = make_mouse_event(
	    QEvent::MouseButtonRelease, 110, 90, Qt::LeftButton,
	    Qt::LeftButton);
    if (!pendingSourceBoxFilter.eventFilter(NULL, &pendingSourceBoxPress) ||
	    !pendingSourceBoxFilter.eventFilter(NULL, &pendingSourceBoxMove) ||
	    !pendingSourceBoxFilter.eventFilter(NULL, &pendingSourceBoxRelease)) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad box select filter should consume pending source-backed exact pick");
    }
    if (!pendingSourceBoxFilter.selected_paths().empty()) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad box select filter should not fall back to legacy selection while exact source pick is pending");
    }
    if (wait_for_qtcad_source_pick_result(sourceService)) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad box select filter pending source pick result should become ready");
    }
    {
	std::vector<BObolLodResult> pendingBoxResults;
	sourceService.drainResults(pendingBoxResults);
	if (pendingBoxResults.empty()) {
	    controller->setExactFullDetailBudget(0, 0);
	    controller->setLodService(NULL);
	    sourceService.stop();
	    FAIL("qtcad box select filter should submit pending source-backed exact pick requests");
	}
	for (const BObolLodResult &result : pendingBoxResults) {
	    if (result.providerStatus != BOBOL_LOD_PROVIDER_FALLBACK) {
		controller->setExactFullDetailBudget(0, 0);
		controller->setLodService(NULL);
		sourceService.stop();
		FAIL("qtcad box select filter pending source pick should honor full-detail budget fallback");
	    }
	}
    }

    QgSelectRayFilter pendingSourceRayFilter;
    pendingSourceRayFilter.dbip = gedp->dbip;
    pendingSourceRayFilter.set_view_widget(&view);
    QMouseEvent pendingSourceRayRelease = make_mouse_event(
	    QEvent::MouseButtonRelease, 90, 70, Qt::LeftButton,
	    Qt::LeftButton);
    if (!pendingSourceRayFilter.eventFilter(NULL, &pendingSourceRayRelease)) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad ray select filter should consume pending source-backed exact pick");
    }
    if (!pendingSourceRayFilter.selected_paths().empty()) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad ray select filter should not fall back to legacy selection while exact source pick is pending");
    }
    if (wait_for_qtcad_source_pick_result(sourceService)) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad ray select filter pending source pick result should become ready");
    }
    {
	std::vector<BObolLodResult> pendingRayResults;
	sourceService.drainResults(pendingRayResults);
	if (pendingRayResults.empty()) {
	    controller->setExactFullDetailBudget(0, 0);
	    controller->setLodService(NULL);
	    sourceService.stop();
	    FAIL("qtcad ray select filter should submit the pending source-backed exact pick request");
	}
	for (const BObolLodResult &result : pendingRayResults) {
	    if (result.providerStatus != BOBOL_LOD_PROVIDER_FALLBACK) {
		controller->setExactFullDetailBudget(0, 0);
		controller->setLodService(NULL);
		sourceService.stop();
		FAIL("qtcad ray select filter pending source pick should honor full-detail budget fallback");
	    }
	}
    }
    controller->setExactFullDetailBudget(0, 0);

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
    if (controller->getRtPickCacheCount() != 1 ||
	    !controller->getRtPickCache(0) ||
	    !controller->getRtPickCache(0)->isReady())
	FAIL("qtcad Obol librt exact implicit pick should retain a controller RT pick cache");
    BObolRtPickCache *implicitRtPickCache = controller->getRtPickCache(0);

    std::vector<QgObolPickRecord> implicitRectPicks;
    if (qg_obol_pick_rect(&view, 80, 60, 100, 80, 8.0f, false,
	    implicitRectPicks) <= 0) {
	FAIL("qtcad Obol rectangle pick should reuse controller-cached librt exact implicit picks");
    }
    if (controller->getRtPickCacheCount() != 1 ||
	    controller->getRtPickCache(0) != implicitRtPickCache)
	FAIL("qtcad Obol rectangle pick should keep reusing the controller RT pick cache");
    if (implicitRectPicks[0].path.find("implicit.r") == std::string::npos ||
	    implicitRectPicks[0].sourceName != "implicit.s" ||
	    implicitRectPicks[0].primitiveKind !=
	    QgObolPickRecord::IMPLICIT_SOLID ||
	    implicitRectPicks[0].distance <= 0.0f) {
	FAIL("qtcad Obol cached rectangle librt pick should preserve RT hit identity");
    }

    std::vector<QgObolPickRecord> implicitRayPicks;
    int implicitRayPickCount = qg_obol_pick_ray(&view,
	    SbVec3f(0.0f, 0.0f, 5.0f),
	    SbVec3f(0.0f, 0.0f, -1.0f), false, implicitRayPicks);
    if (implicitRayPickCount != 1) {
	fprintf(stderr, "implicit ray pick count=%d cache_count=%d\n",
		implicitRayPickCount, controller->getRtPickCacheCount());
	FAIL("qtcad Obol ray pick should reuse controller-cached librt exact implicit picks");
    }
    if (controller->getRtPickCacheCount() != 1 ||
	    controller->getRtPickCache(0) != implicitRtPickCache)
	FAIL("qtcad Obol ray pick should keep reusing the controller RT pick cache");
    if (implicitRayPicks[0].path.find("implicit.r") == std::string::npos ||
	    implicitRayPicks[0].sourceName != "implicit.s" ||
	    implicitRayPicks[0].primitiveKind !=
	    QgObolPickRecord::IMPLICIT_SOLID ||
	    implicitRayPicks[0].distance <= 0.0f)
	FAIL("qtcad Obol cached ray librt pick should preserve RT hit identity");

    QgSelectRayFilter implicitRayFilter;
    implicitRayFilter.dbip = gedp->dbip;
    implicitRayFilter.set_view_widget(&view);
    QMouseEvent implicitRayRelease = make_mouse_event(QEvent::MouseButtonRelease,
	    90, 70, Qt::LeftButton, Qt::LeftButton);
    if (!implicitRayFilter.eventFilter(NULL, &implicitRayRelease))
	FAIL("qtcad select ray filter should accept the Obol explicit-ray implicit workflow");
    const std::vector<std::string> &implicitRayPaths =
	implicitRayFilter.selected_paths();
    if (implicitRayPaths.size() != 1 ||
	    implicitRayPaths[0].find("implicit.r") == std::string::npos)
	FAIL("qtcad select ray filter should expose Obol explicit-ray implicit paths");
    if (controller->getRtPickCacheCount() != 1 ||
	    controller->getRtPickCache(0) != implicitRtPickCache)
	FAIL("qtcad select ray filter should reuse the controller RT pick cache");

    QgSelectRayFilter noLegacyDbRayFilter;
    noLegacyDbRayFilter.set_view_widget(&view);
    QMouseEvent noLegacyDbRayRelease = make_mouse_event(QEvent::MouseButtonRelease,
	    90, 70, Qt::LeftButton, Qt::LeftButton);
    if (!noLegacyDbRayFilter.eventFilter(NULL, &noLegacyDbRayRelease))
	FAIL("qtcad select ray filter should run Obol explicit-ray selection without legacy dbip");
    const std::vector<std::string> &noLegacyDbRayPaths =
	noLegacyDbRayFilter.selected_paths();
    if (noLegacyDbRayPaths.size() != 1 ||
	    noLegacyDbRayPaths[0].find("implicit.r") == std::string::npos)
	FAIL("qtcad select ray filter without legacy dbip should expose Obol explicit-ray implicit paths");
    if (controller->getRtPickCacheCount() != 1 ||
	    controller->getRtPickCache(0) != implicitRtPickCache)
	FAIL("qtcad select ray filter without legacy dbip should reuse the controller RT pick cache");

    SoSeparator *mixedRoot = new SoSeparator;
    mixedRoot->ref();
    SoBRLDatabaseSource *farBoxDatabase = new SoBRLDatabaseSource;
    farBoxDatabase->configureDatabaseSource("farbox.s", gedp->dbip,
	    SoBRLDatabaseSource::SHADED, 6);
    if (!farBoxDatabase->realizeDatabaseMesh())
	FAIL("qtcad mixed pick fixture should realize the farther display box");
    mixedRoot->addChild(farBoxDatabase);
    SoBRLDatabaseSource *mixedImplicitDatabase = new SoBRLDatabaseSource;
    mixedImplicitDatabase->configureDatabaseSource("implicit.r", gedp->dbip,
	    SoBRLDatabaseSource::WIREFRAME, 7);
    mixedRoot->addChild(mixedImplicitDatabase);
    controller->setSceneRoot(mixedRoot);
    mixedRoot->unref();
    if (camera)
	camera->position = SbVec3f(0.0f, 0.0f, 5.0f);

    std::vector<QgObolPickRecord> mixedNearestPicks;
    if (qg_obol_pick_point(&view, 90, 70, 8.0f, false,
	    mixedNearestPicks) != 1)
	FAIL("qtcad mixed Obol/RT point pick should return one nearest hit");
    if (mixedNearestPicks[0].sourceName != "implicit.s" ||
	    mixedNearestPicks[0].primitiveKind !=
	    QgObolPickRecord::IMPLICIT_SOLID ||
	    mixedNearestPicks[0].distance <= 0.0f)
	FAIL("qtcad mixed Obol/RT point pick should choose the nearer librt hit");

    std::vector<QgObolPickRecord> mixedAllPicks;
    if (qg_obol_pick_point(&view, 90, 70, 8.0f, true,
	    mixedAllPicks) < 2)
	FAIL("qtcad mixed Obol/RT pick-all should keep display and librt hits");
    bool foundFarBox = false;
    for (const QgObolPickRecord &pick : mixedAllPicks) {
	if (pick.path == "/farbox.s" || pick.sourceName == "farbox.s") {
	    foundFarBox = true;
	    break;
	}
    }
    if (mixedAllPicks[0].sourceName != "implicit.s" || !foundFarBox)
	FAIL("qtcad mixed Obol/RT pick-all should order merged hits by distance");

    std::vector<QgObolPickRecord> mixedRayAllPicks;
    if (qg_obol_pick_ray(&view, SbVec3f(0.0f, 0.0f, 5.0f),
	    SbVec3f(0.0f, 0.0f, -1.0f), true, mixedRayAllPicks) < 2)
	FAIL("qtcad mixed Obol/RT ray pick-all should keep display and librt hits");
    bool foundRayFarBox = false;
    for (const QgObolPickRecord &pick : mixedRayAllPicks) {
	if (pick.path == "/farbox.s" || pick.sourceName == "farbox.s") {
	    foundRayFarBox = true;
	    break;
	}
    }
    if (mixedRayAllPicks[0].sourceName != "implicit.s" || !foundRayFarBox)
	FAIL("qtcad mixed Obol/RT ray pick-all should order merged hits by distance");

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
