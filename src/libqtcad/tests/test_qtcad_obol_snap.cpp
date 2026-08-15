/*              T E S T _ Q T C A D _ O B O L _ S N A P . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "bv.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BLodMeshShape.h"
#include "BObol/BLodService.h"
#include "BObol/BMeshShape.h"
#include "BObol/BSnapAction.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "ged.h"
#include "ged/draw.h"
#include "QgSceneSyncPrivate.h"
#include "qtcad/QgObolSnap.h"
#include "qtcad/QgPolyFilter.h"
#include "qtcad/QgView.h"
#include "qtcad/QgViewFilter.h"
#include "raytrace.h"
#include "wdb.h"

#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoSeparator.h>

#include <QApplication>
#include <QMouseEvent>

#include <algorithm>
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

class SnapProbeFilter : public QgViewFilter {
public:
    bool sync(QEvent *event)
    {
	return view_sync(event) != NULL;
    }
};

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

static int
make_snap_db(const char *dbpath)
{
    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp)
	return 0;

    point_t bmin = {-1.0, -1.0, -1.0};
    point_t bmax = { 1.0,  1.0,  1.0};
    int ret = mk_rpp(wdbp, "box.s", bmin, bmax) == 0;

    point_t omin = { 9.0,  9.0,  9.0};
    point_t omax = {11.0, 11.0, 11.0};
    ret = ret && mk_rpp(wdbp, "snap_only.s", omin, omax) == 0;

    fastf_t lod_vertices[9] = {
	-1.0, -1.0, 0.0,
	 1.0, -1.0, 0.0,
	-1.0,  1.0, 0.0
    };
    int lod_faces[3] = {0, 1, 2};
    ret = ret && mk_bot(wdbp, "lod-snap.bot", RT_BOT_SOLID,
	    RT_BOT_UNORIENTED, 0, 3, 1, lod_vertices, lod_faces,
	    NULL, NULL) == 0;

    wdb_close(wdbp);
    return ret;
}

static int
apply_and_sync(struct ged *gedp,
	QgView *view,
	struct ged_view_context *view_ctx,
	const char *path,
	enum ged_scene_draw_mode mode)
{
    if (!qg_scene_bind(gedp, view))
	return 0;
    struct ged_scene_draw_request request;
    ged_scene_draw_request_init(&request);
    request.view = view_ctx;
    request.paths = &path;
    request.path_count = 1;
    request.style.draw_mode = mode;
    request.realization.mode = GED_SCENE_REALIZE_EAGER;
    struct ged_scene_result *result = ged_scene_result_create();
    const enum ged_scene_status status = ged_scene_draw(gedp, &request,
	result);
    const int changed = ged_scene_result_changed(result);
    ged_scene_result_destroy(result);
    return status == GED_SCENE_OK && changed;
}

static BObolLodResult
qtcad_source_snap_result(const BObolLodRequest &request,
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
    result.geometry.activeCut =
	resultKind == BOBOL_LOD_RESULT_MESH ? 1 : -1;
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
qtcad_source_snap_task(const BObolLodRequest &request, void *UNUSED(userData))
{
    return qtcad_source_snap_result(request, BOBOL_LOD_RESULT_FULL_DETAIL,
	    BOBOL_LOD_QUALITY_FULL_DETAIL);
}

static int
wait_for_qtcad_source_snap_result(BObolLodService &service)
{
    for (int i = 0; i < 400; i++) {
	if (service.inFlightCount() == 0 &&
		service.queuedResultCountForDiagnostics() >= 1)
	    return 0;
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    return 1;
}

static int
qtcad_make_source_snap_request(BObolViewController *controller,
	const SbVec3f &query,
	float tolerance,
	uint32_t enabledKinds,
	BObolLodRequest &request)
{
    if (!controller || !controller->getViewport() ||
	    !controller->getViewport()->getRoot())
	return 0;

    SoBRLSnapAction snapAction;
    snapAction.setQueryPoint(query);
    snapAction.setTolerance(tolerance);
    snapAction.setEnabledKinds(enabledKinds);
    snapAction.setPriorityPolicy(SoBRLSnapAction::FEATURE_PRIORITY);
    snapAction.setGeometryPolicy(SoBRLSnapAction::FULL_DETAIL);
    snapAction.apply(controller->getViewport()->getRoot());
    return snapAction.getSourceBackedFullDetailRequestCount() == 1 &&
	snapAction.makeSourceBackedFullDetailLodRequest(0, request);
}

static int
qtcad_submit_source_snap_result(BObolLodService &service,
	const BObolLodRequest &request)
{
    BObolLodTask task;
    task.request = request;
    task.realize = qtcad_source_snap_task;
    if (service.submit(task) == 0)
	return 0;

    return wait_for_qtcad_source_snap_result(service) == 0;
}

static QMouseEvent
left_move_at(int x, int y)
{
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    return QMouseEvent(QEvent::MouseMove, QPoint(x, y),
	    Qt::NoButton, Qt::NoButton, Qt::NoModifier);
#else
    return QMouseEvent(QEvent::MouseMove, QPointF(x, y),
	    QPointF(x, y), Qt::NoButton, Qt::NoButton, Qt::NoModifier);
#endif
}

static QMouseEvent
left_press_at(int x, int y)
{
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    return QMouseEvent(QEvent::MouseButtonPress, QPoint(x, y),
	    Qt::LeftButton, Qt::LeftButton, Qt::NoModifier);
#else
    return QMouseEvent(QEvent::MouseButtonPress, QPointF(x, y),
	    QPointF(x, y), Qt::LeftButton, Qt::LeftButton, Qt::NoModifier);
#endif
}

static QMouseEvent
left_drag_at(int x, int y)
{
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    return QMouseEvent(QEvent::MouseMove, QPoint(x, y),
	    Qt::NoButton, Qt::LeftButton, Qt::NoModifier);
#else
    return QMouseEvent(QEvent::MouseMove, QPointF(x, y),
	    QPointF(x, y), Qt::NoButton, Qt::LeftButton, Qt::NoModifier);
#endif
}

static void
set_center_query(void *view, fastf_t x, fastf_t y, fastf_t z)
{
    mat_t view2model;
    mat_t model2view;
    mat_t rotation;
    (void)bv_context_dimensions_set(
	static_cast<struct bv_context *>(view), 200, 200);
    bv_size_set(bv_context_view(static_cast<struct bv_context *>(view)), 2.0);
    MAT_IDN(rotation);
    bv_rotation_set(bv_context_view(static_cast<struct bv_context *>(view)),
	rotation);
    MAT_IDN(view2model);
    MAT_DELTAS(view2model, x, y, z);
    bv_view2model_set(bv_context_view(static_cast<struct bv_context *>(view)), view2model);
    MAT_IDN(model2view);
    MAT_DELTAS(model2view, -x, -y, -z);
    bv_model2view_set(bv_context_view(static_cast<struct bv_context *>(view)), model2view);
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    QApplication app(argc, argv);

    const char *dbpath = "qtcad_obol_snap_tmp.g";
    if (!make_snap_db(dbpath))
	FAIL("failed to create qtcad Obol snap test database");

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open qtcad Obol snap test database");

    QgView view(NULL, QgViewType::SW);
    view.resize(200, 200);
    struct ged_view_context *view_ctx =
	ged_view_context_from_bv(view.viewContext());
    ged_view_active_ctx_set(gedp, view_ctx);
    (void)ged_view_context_host_attach(gedp, view_ctx);

    BObolViewController *controller = view.obolViewController();
    if (!controller)
	FAIL("QgView should expose an Obol controller");
    controller->clearDatabaseSources();
    controller->setViewportSize(200, 200);

    if (!apply_and_sync(gedp, &view, view_ctx, "box.s",
	GED_SCENE_DRAW_WIRE))
	FAIL("GED wire draw should sync box source into Obol");

    QgObolSnapRecord record;
    if (!qg_obol_snap_point(&view, SbVec3f(1.02f, 1.02f, 1.02f),
	    0.1f, QgObolSnapRecord::ENDPOINT, record))
	FAIL("qtcad Obol snap helper should find a box endpoint");
    if (record.path != "/box.s" ||
	    record.kind != QgObolSnapRecord::ENDPOINT ||
	    record.primitiveIndex < 0 ||
	    !near_point(record.point, 1.0f, 1.0f, 1.0f))
	FAIL("qtcad Obol snap helper should preserve endpoint identity");

    if (controller->replaceDatabaseSource("snap_only.s", gedp->dbip,
	    SoBRLDatabaseSource::WIREFRAME, 2) <= 0 ||
	    !controller->realizePending())
	FAIL("test should add an Obol-only snap source");

    set_center_query(view.viewContext(), 11.02, 11.02, 11.02);
    bv_snap_source_flags_set(bv_context_view(static_cast<struct bv_context *>(view.viewContext())), BV_SNAP_DB);
    bv_snap_lines_set(bv_context_view(static_cast<struct bv_context *>(view.viewContext())), 1);

    SnapProbeFilter filter;
    filter.set_view_widget(&view);
    QMouseEvent move = left_move_at(100, 100);
    if (!filter.sync(&move))
	FAIL("qtcad view filter should accept a snap probe mouse event");
    point_t snapped_point = VINIT_ZERO;
    bv_current_point_get(snapped_point, bv_context_view_const(static_cast<const struct bv_context *>(view.viewContext())));
    if (!nearly_equal((float)snapped_point[X], 11.0f) ||
	    !nearly_equal((float)snapped_point[Y], 11.0f) ||
	    !nearly_equal((float)snapped_point[Z], 11.0f))
	FAIL("qtcad view filter should refine database snapping through Obol");

    set_center_query(view.viewContext(), 11.02, 11.02, 11.02);
    bv_snap_source_flags_set(bv_context_view(static_cast<struct bv_context *>(view.viewContext())), BV_SNAP_VIEW);
    QMouseEvent viewScopedMove = left_move_at(100, 100);
    if (!filter.sync(&viewScopedMove))
	FAIL("qtcad view filter should accept a view-scoped snap probe event");
    bv_current_point_get(snapped_point, bv_context_view_const(static_cast<const struct bv_context *>(view.viewContext())));
    if (nearly_equal((float)snapped_point[X], 11.0f) &&
	    nearly_equal((float)snapped_point[Y], 11.0f) &&
	    nearly_equal((float)snapped_point[Z], 11.0f))
	FAIL("qtcad view-scoped snapping should not apply database refinement");

    /* Exercise polygon filters end-to-end.  The shared view filter owns raw
     * screen conversion and snapping; polygon filters must consume its model
     * point rather than converting the original Qt coordinates again. */
    set_center_query(view.viewContext(), 0.0, 0.0, 0.0);
    struct bv *bv_view = bv_context_view(
	static_cast<struct bv_context *>(view.viewContext()));
    bv_snap_source_flags_set(bv_view, BV_SNAP_VIEW);
    bv_snap_lines_set(bv_view, 0);
    struct bv_grid_state grid = BV_GRID_STATE_INIT;
    grid.snap = 1;
    grid.res_h = 0.5;
    grid.res_v = 0.5;
    if (!bv_grid_state_set(bv_view, &grid))
	FAIL("polygon snap test should configure a model grid");

    QgPolyCreateFilter gridCreate;
    gridCreate.set_view_widget(&view);
    gridCreate.ptype = QG_POLYGON_CIRCLE;
    QMouseEvent gridPress = left_press_at(106, 94);
    if (!gridCreate.eventFilter(NULL, &gridPress) ||
	ged_view_polygon_ref_is_null(gridCreate.polygon))
	FAIL("polygon create filter should accept a grid-snap press");
    struct ged_view_polygon_record polygonRecord;
    if (!ged_view_polygon_record_get(view_ctx, gridCreate.polygon,
	    &polygonRecord) ||
	!nearly_equal((float)polygonRecord.origin_point[X], 0.0f) ||
	!nearly_equal((float)polygonRecord.origin_point[Y], 0.0f))
	FAIL("polygon creation should use the snapped grid model point");
    if (!ged_view_polygon_set_selected(view_ctx, gridCreate.polygon, 1) ||
	!ged_view_polygon_record_get(view_ctx, gridCreate.polygon,
	    &polygonRecord) || !polygonRecord.selected)
	FAIL("polygon selection should be retained in the shared scene record");
    if (!ged_view_polygon_clear_selection(view_ctx) ||
	!ged_view_polygon_record_get(view_ctx, gridCreate.polygon,
	    &polygonRecord) || polygonRecord.selected)
	FAIL("polygon selection clearing should update the shared scene record");

    QMouseEvent gridDrag = left_drag_at(146, 100);
    if (!gridCreate.eventFilter(NULL, &gridDrag))
	FAIL("polygon create filter should accept a snapped shape drag");
    struct bg_polygon gridGeometry = BG_POLYGON_INIT_ZERO;
    if (!controller->polygons().copyGeometry(
	    BObolPolygonHandle(gridCreate.polygon.id, 1), &gridGeometry) ||
	gridGeometry.num_contours != 1 ||
	gridGeometry.contour[0].num_points < 4)
	FAIL("polygon grid-snap drag should produce constrained geometry");
    point_t dragPoint = VINIT_ZERO;
    if (!bv_current_point_get(dragPoint, bv_view) ||
	!nearly_equal((float)dragPoint[X], 0.5f) ||
	!nearly_equal((float)dragPoint[Y], 0.0f))
	FAIL("polygon shape drag should retain the snapped model point");
    ged_view_polygon_ref control = ged_view_polygon_create(view_ctx,
	"polygon-grid-control", 1, GED_VIEW_POLYGON_CIRCLE,
	polygonRecord.origin_point);
    struct bg_polygon controlGeometry = BG_POLYGON_INIT_ZERO;
    if (ged_view_polygon_ref_is_null(control) ||
	!ged_view_polygon_update_model_pt(view_ctx, control, dragPoint,
	    GED_VIEW_POLYGON_UPDATE_DEFAULT) ||
	!controller->polygons().copyGeometry(BObolPolygonHandle(control.id, 1),
	    &controlGeometry) ||
	controlGeometry.num_contours != gridGeometry.num_contours ||
	controlGeometry.contour[0].num_points !=
	    gridGeometry.contour[0].num_points)
	FAIL("polygon grid-snap control geometry should be available");
    for (size_t i = 0; i < gridGeometry.contour[0].num_points; i++) {
	if (!VNEAR_EQUAL(gridGeometry.contour[0].point[i],
		controlGeometry.contour[0].point[i], 0.0001))
	    FAIL("polygon update should use the snapped grid model point");
    }
    bg_polygon_clear(&gridGeometry);
    bg_polygon_clear(&controlGeometry);
    if (!ged_view_polygon_remove(view_ctx, control))
	FAIL("polygon grid-snap control should be removable");

    QgPolyMoveFilter gridMove;
    gridMove.set_view_widget(&view);
    gridMove.polygon = gridCreate.polygon;
    QMouseEvent movePress = left_press_at(100, 100);
    QMouseEvent moveDrag = left_drag_at(146, 100);
    const bool movePressHandled = gridMove.eventFilter(NULL, &movePress);
    point_t moveStart = VINIT_ZERO;
    if (!movePressHandled || !bv_current_point_get(moveStart, bv_view))
	FAIL("polygon move should retain its snapped start point");
    ged_view_polygon_ref moveControl = ged_view_polygon_create(view_ctx,
	"polygon-grid-move-control", 1, GED_VIEW_POLYGON_CIRCLE,
	polygonRecord.origin_point);
    if (ged_view_polygon_ref_is_null(moveControl) ||
	!ged_view_polygon_update_model_pt(view_ctx, moveControl, dragPoint,
	    GED_VIEW_POLYGON_UPDATE_DEFAULT) ||
	!ged_view_polygon_move(view_ctx, moveControl, dragPoint, moveStart))
	FAIL("polygon move control should accept snapped model points");
    struct ged_view_polygon_record moveControlRecord;
    if (!ged_view_polygon_record_get(view_ctx, moveControl,
	    &moveControlRecord))
	FAIL("polygon move control record should be available");
    const bool moveDragHandled = gridMove.eventFilter(NULL, &moveDrag);
    const int moveRecord = ged_view_polygon_record_get(view_ctx,
	gridCreate.polygon, &polygonRecord);
    if (!movePressHandled || !moveDragHandled || !moveRecord ||
	!VNEAR_EQUAL(polygonRecord.origin_point,
	    moveControlRecord.origin_point, 0.0001))
	FAIL("polygon move should use snapped current and previous model points");
    if (!ged_view_polygon_remove(view_ctx, moveControl))
	FAIL("polygon move control should be removable");
    if (!ged_view_polygon_remove(view_ctx, gridCreate.polygon))
	FAIL("polygon grid-snap fixture should be removable");

    grid.snap = 0;
    if (!bv_grid_state_set(bv_view, &grid))
	FAIL("polygon line-snap test should disable grid snapping");
    point_t lineStart = {-0.5, 0.0, 0.0};
    point_t lineEnd = {0.5, 0.0, 0.0};
    ged_view_polygon_ref snapLine = ged_view_polygon_create(view_ctx,
	"polygon-snap-line", 1, GED_VIEW_POLYGON_GENERAL, lineStart);
    if (ged_view_polygon_ref_is_null(snapLine) ||
	!ged_view_polygon_update_model_pt(view_ctx, snapLine, lineEnd,
	    GED_VIEW_POLYGON_UPDATE_PT_APPEND))
	FAIL("polygon line-snap fixture should be created through the public API");
    bv_snap_lines_set(bv_view, 1);
    ged_view_polygon_snap_exclude_set(view_ctx, GED_VIEW_POLYGON_REF_NULL);

    QgObolSnapRecord lineProbe;
    const int lineProbeFound = qg_obol_snap_point_filtered(&view,
	SbVec3f(0.2f, 0.05f, 0.0f), 0.1f,
	QgObolSnapRecord::ENDPOINT | QgObolSnapRecord::MIDPOINT |
	QgObolSnapRecord::LINE_NEAREST,
	SoBRLSnapAction::VIEW_SOURCES, lineProbe);
    if (!lineProbeFound || lineProbe.path != "polygon-snap-line" ||
	lineProbe.kind != QgObolSnapRecord::LINE_NEAREST ||
	!nearly_equal(lineProbe.point[0], 0.2f) ||
	!nearly_equal(lineProbe.point[1], 0.0f))
	FAIL("Obol source filtering should retain view-polygon line candidates");

    QgPolyCreateFilter lineCreate;
    lineCreate.set_view_widget(&view);
    lineCreate.ptype = QG_POLYGON_CIRCLE;
    QMouseEvent linePress = left_press_at(120, 95);
    if (!lineCreate.eventFilter(NULL, &linePress) ||
	ged_view_polygon_ref_is_null(lineCreate.polygon) ||
	!ged_view_polygon_record_get(view_ctx, lineCreate.polygon,
	    &polygonRecord) ||
	!nearly_equal((float)polygonRecord.origin_point[X], 0.2f) ||
	!nearly_equal((float)polygonRecord.origin_point[Y], 0.0f))
	FAIL("polygon creation should snap to view-polygon segments");
    if (!ged_view_polygon_remove(view_ctx, lineCreate.polygon) ||
	!ged_view_polygon_remove(view_ctx, snapLine))
	FAIL("polygon line-snap fixtures should be removable");
    bv_snap_lines_set(bv_view, 0);

    SoSeparator *lodRoot = new SoSeparator;
    lodRoot->ref();
    SoBRLDatabaseSource *lodDatabase = new SoBRLDatabaseSource;
    lodDatabase->configureDatabaseSource("lod-snap-db", gedp->dbip,
	    SoBRLDatabaseSource::WIREFRAME, 1);
    SoBRLLodMeshShape *lodMesh = new SoBRLLodMeshShape;
    SbVec3f lodPoints[3] = {
	SbVec3f(-1.0f, -1.0f, 0.0f),
	SbVec3f(1.0f, -1.0f, 0.0f),
	SbVec3f(-1.0f, 1.0f, 0.0f)
    };
    int32_t lodIndices[3] = {0, 1, 2};
    lodMesh->sourcePath = "/lod-snap.bot";
    lodMesh->sourceName = "lod-snap.bot";
    lodMesh->sourceType = "bot";
    lodMesh->sourceId = 9101;
    lodMesh->editIntentId = "edit::lod-snap/face";
    lodMesh->editIntentRole = "exact-snap";
    lodMesh->setIndexedTriangles(lodPoints, 3, lodIndices, 3);

    BObolLodRequest displayRequest;
    lodMesh->makeLodRequest(displayRequest,
	    "db://qtcad-obol-snap-test",
	    1,
	    1,
	    1,
	    BOBOL_LOD_DRAW_SHADED,
	    "bobol_mesh_lod",
	    "bobol-cache-v1",
	    BOBOL_LOD_QUALITY_FAST_DISPLAY);
    BObolLodResult displayResult = qtcad_source_snap_result(displayRequest,
	    BOBOL_LOD_RESULT_MESH, BOBOL_LOD_QUALITY_FAST_DISPLAY);
    if (!lodMesh->applyStagedLodResult(displayResult, &displayRequest) ||
	    !lodMesh->isLodDisplayActive() ||
	    lodMesh->hasFullDetailMesh())
	FAIL("qtcad LoD snap fixture should have active display LoD without resident full detail");

    lodRoot->addChild(lodDatabase);
    lodRoot->addChild(lodMesh);
    controller->setSceneRoot(lodRoot);
    lodRoot->unref();

    BObolLodRequest sourceLodRequest;
    if (!qtcad_make_source_snap_request(controller,
	    SbVec3f(-0.2f, -0.2f, 0.02f), 0.1f,
	    QgObolSnapRecord::FACE_NEAREST, sourceLodRequest))
	FAIL("qtcad LoD snap fixture should build a bounded source full-detail request");

    BObolLodService sourceService;
    if (!sourceService.start(1, TRUE))
	FAIL("qtcad LoD snap source service should start");
    controller->setLodService(&sourceService);

    if (!qtcad_submit_source_snap_result(sourceService, sourceLodRequest)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad LoD snap source result should become ready");
    }

    QgObolSnapRecord exactSnap;
    if (!qg_obol_snap_point_full_detail(&view,
	    SbVec3f(-0.2f, -0.2f, 0.02f), 0.1f,
	    QgObolSnapRecord::FACE_NEAREST, exactSnap)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol snap should consume ready source-backed full detail");
    }
    if (exactSnap.path != "/lod-snap.bot" ||
	    exactSnap.editIntentId != "edit::lod-snap/face" ||
	    exactSnap.editIntentRole != "exact-snap" ||
	    exactSnap.kind != QgObolSnapRecord::FACE_NEAREST ||
	    exactSnap.primitiveIndex != 42 ||
	    exactSnap.vertexIndex != -1 ||
	    exactSnap.edgeSlot != -1 ||
	    !near_point(exactSnap.point, -0.2f, -0.2f, 0.0f)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol snap should preserve source-backed face identity and edit intent");
    }
    if (sourceService.queuedResultCountForDiagnostics() != 0) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol snap should drain only its matching result");
    }

    BObolLodRequest vertexLodRequest;
    if (!qtcad_make_source_snap_request(controller,
	    SbVec3f(-0.96f, -0.96f, 0.0f), 0.1f,
	    QgObolSnapRecord::VERTEX, vertexLodRequest) ||
	    !qtcad_submit_source_snap_result(sourceService, vertexLodRequest)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol vertex snap source result should become ready");
    }
    QgObolSnapRecord vertexSnap;
    if (!qg_obol_snap_point_full_detail(&view,
	    SbVec3f(-0.96f, -0.96f, 0.0f), 0.1f,
	    QgObolSnapRecord::VERTEX, vertexSnap) ||
	    vertexSnap.path != "/lod-snap.bot" ||
	    vertexSnap.kind != QgObolSnapRecord::VERTEX ||
	    vertexSnap.primitiveIndex != 42 ||
	    vertexSnap.vertexIndex != 10 ||
	    vertexSnap.edgeSlot != -1 ||
	    vertexSnap.edgeVertexIndexA != -1 ||
	    vertexSnap.edgeVertexIndexB != -1 ||
	    !near_point(vertexSnap.point, -1.0f, -1.0f, 0.0f)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol snap should preserve source-backed vertex identity");
    }
    if (sourceService.queuedResultCountForDiagnostics() != 0) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol vertex snap should drain only its matching result");
    }

    BObolLodRequest edgeLodRequest;
    if (!qtcad_make_source_snap_request(controller,
	    SbVec3f(0.0f, -0.97f, 0.0f), 0.1f,
	    QgObolSnapRecord::EDGE_NEAREST, edgeLodRequest) ||
	    !qtcad_submit_source_snap_result(sourceService, edgeLodRequest)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol edge snap source result should become ready");
    }
    QgObolSnapRecord edgeSnap;
    if (!qg_obol_snap_point_full_detail(&view,
	    SbVec3f(0.0f, -0.97f, 0.0f), 0.1f,
	    QgObolSnapRecord::EDGE_NEAREST, edgeSnap) ||
	    edgeSnap.path != "/lod-snap.bot" ||
	    edgeSnap.kind != QgObolSnapRecord::EDGE_NEAREST ||
	    edgeSnap.primitiveIndex != 42 ||
	    edgeSnap.vertexIndex != -1 ||
	    edgeSnap.edgeSlot != 0 ||
	    edgeSnap.edgeVertexIndexA != 10 ||
	    edgeSnap.edgeVertexIndexB != 11 ||
	    !near_point(edgeSnap.point, 0.0f, -1.0f, 0.0f)) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol snap should preserve source-backed edge identity");
    }
    if (sourceService.queuedResultCountForDiagnostics() != 0) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol edge snap should drain only its matching result");
    }
    if (edgeSnap.sourceFullDetailPending ||
	    edgeSnap.submittedSourceRequestCount != 0) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol snap should not report pending source detail after consuming a ready result");
    }

    controller->setExactFullDetailBudget(0, 2);
    QgObolSnapRecord overBudgetSnap;
    if (qg_obol_snap_point_full_detail(&view,
	    SbVec3f(-0.2f, -0.2f, 0.02f), 0.1f,
	    QgObolSnapRecord::FACE_NEAREST, overBudgetSnap)) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol snap should not treat missing over-budget source detail as a snap candidate");
    }
    if (!overBudgetSnap.sourceFullDetailPending ||
	    overBudgetSnap.submittedSourceRequestCount != 1) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol snap should report pending submitted source detail");
    }
    if (wait_for_qtcad_source_snap_result(sourceService)) {
	controller->setExactFullDetailBudget(0, 0);
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad exact Obol snap over-budget source result should become ready");
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
	FAIL("qtcad exact Obol snap should pass controller full-detail budget to source provider");
    }
    controller->setExactFullDetailBudget(0, 0);

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
