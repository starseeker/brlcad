/*              T E S T _ Q T C A D _ O B O L _ S N A P . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

extern "C" {
#include "rt/view_legacy_bsg.h"
}

#include "brlobol/database_source.h"
#include "brlobol/lod_mesh_shape.h"
#include "brlobol/lod_service.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/snap_action.h"
#include "brlobol/view_controller.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "ged.h"
#include "ged/bsg_ged_draw.h"
#include "qtcad/QgLegacyViewBsg.h"
#include "qtcad/QgObolDrawSync.h"
#include "qtcad/QgObolSnap.h"
#include "qtcad/QgView.h"
#include "qtcad/QgViewFilter.h"
#include "raytrace.h"
#include "wdb.h"

#include <Inventor/SoViewport.h>
#include <Inventor/nodes/SoSeparator.h>

#include <QApplication>
#include <QMouseEvent>

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
qtcad_source_snap_result(const BRLObolLodRequest &request,
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
qtcad_source_snap_task(const BRLObolLodRequest &request, void *UNUSED(userData))
{
    return qtcad_source_snap_result(request, BRLOBOL_LOD_RESULT_FULL_DETAIL,
	    BRLOBOL_LOD_QUALITY_FULL_DETAIL);
}

static int
wait_for_qtcad_source_snap_result(BRLObolLodService &service)
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
qtcad_make_source_snap_request(BRLObolViewController *controller,
	const SbVec3f &query,
	float tolerance,
	uint32_t enabledKinds,
	BRLObolLodRequest &request)
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
qtcad_submit_source_snap_result(BRLObolLodService &service,
	const BRLObolLodRequest &request)
{
    BRLObolLodTask task;
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

static void
set_center_query(struct bsg_view *v, fastf_t x, fastf_t y, fastf_t z)
{
    mat_t view2model;
    rt_view_dimensions_set_bsg(v, 200, 200);
    rt_view_size_set_bsg(v, 2.0);
    MAT_IDN(view2model);
    MAT_DELTAS(view2model, x, y, z);
    rt_view_view2model_set_bsg(v, view2model);
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

    QgView view(NULL, QgView_SW, NULL);
    view.resize(200, 200);
    gedp->ged_gvp = qg_legacy_view_to_bsg(view.view());

    BRLObolViewController *controller = view.obolViewController();
    if (!controller)
	FAIL("QgView should expose an Obol controller");
    controller->clearDatabaseSources();
    controller->setViewportSize(200, 200);

    struct ged_draw_transaction draw_box =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "box.s");
    draw_box.view = qg_legacy_view_to_bsg(view.view());
    if (!apply_and_sync(gedp, &view, &draw_box))
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
    if (view.dmp())
	FAIL("qtcad Obol snap helper should not initialize the legacy display manager");

    if (controller->replaceDatabaseSource("snap_only.s", gedp->dbip,
	    SoBRLDatabaseSource::WIREFRAME, 2) <= 0 ||
	    !controller->realizePending())
	FAIL("test should add an Obol-only snap source");
    if (view.dmp())
	FAIL("qtcad Obol snap setup should not initialize the legacy display manager");

    struct bsg_view *bv = qg_legacy_view_to_bsg(view.view());
    set_center_query(bv, 11.02, 11.02, 11.02);
    rt_view_snap_source_flags_set_bsg(bv, RT_VIEW_SNAP_DB_BSG);
    rt_view_snap_lines_set_bsg(bv, 1);

    SnapProbeFilter filter;
    filter.set_view_widget(&view);
    QMouseEvent move = left_move_at(100, 100);
    if (!filter.sync(&move))
	FAIL("qtcad view filter should accept a snap probe mouse event");
    point_t snapped_point = VINIT_ZERO;
    rt_view_current_point_from_bsg(snapped_point, bv);
    if (!nearly_equal((float)snapped_point[X], 11.0f) ||
	    !nearly_equal((float)snapped_point[Y], 11.0f) ||
	    !nearly_equal((float)snapped_point[Z], 11.0f))
	FAIL("qtcad view filter should refine database snapping through Obol");
    if (view.dmp())
	FAIL("qtcad Obol snap refinement should not initialize the legacy display manager");

    set_center_query(bv, 11.02, 11.02, 11.02);
    rt_view_snap_source_flags_set_bsg(bv, RT_VIEW_SNAP_VIEW_BSG);
    QMouseEvent viewScopedMove = left_move_at(100, 100);
    if (!filter.sync(&viewScopedMove))
	FAIL("qtcad view filter should accept a view-scoped snap probe event");
    rt_view_current_point_from_bsg(snapped_point, bv);
    if (nearly_equal((float)snapped_point[X], 11.0f) &&
	    nearly_equal((float)snapped_point[Y], 11.0f) &&
	    nearly_equal((float)snapped_point[Z], 11.0f))
	FAIL("qtcad view filter should leave view-scoped snapping on the BSG path");
    if (view.dmp())
	FAIL("qtcad view-scoped snap fallback should not initialize the legacy display manager");

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

    BRLObolLodRequest displayRequest;
    lodMesh->makeLodRequest(displayRequest,
	    "db://qtcad-obol-snap-test",
	    1,
	    1,
	    1,
	    BRLOBOL_LOD_DRAW_SHADED,
	    "rt_mesh_lod",
	    "rt-cache-v1",
	    BRLOBOL_LOD_QUALITY_FAST_DISPLAY);
    BRLObolLodResult displayResult = qtcad_source_snap_result(displayRequest,
	    BRLOBOL_LOD_RESULT_MESH, BRLOBOL_LOD_QUALITY_FAST_DISPLAY);
    if (!lodMesh->applyStagedLodResult(displayResult, &displayRequest) ||
	    !lodMesh->isLodDisplayActive() ||
	    lodMesh->hasFullDetailMesh())
	FAIL("qtcad LoD snap fixture should have active display LoD without resident full detail");

    lodRoot->addChild(lodDatabase);
    lodRoot->addChild(lodMesh);
    controller->setSceneRoot(lodRoot);
    lodRoot->unref();

    BRLObolLodRequest sourceLodRequest;
    if (!qtcad_make_source_snap_request(controller,
	    SbVec3f(-0.2f, -0.2f, 0.02f), 0.1f,
	    QgObolSnapRecord::FACE_NEAREST, sourceLodRequest))
	FAIL("qtcad LoD snap fixture should build a bounded source full-detail request");

    BRLObolLodService sourceService;
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

    BRLObolLodRequest vertexLodRequest;
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

    BRLObolLodRequest edgeLodRequest;
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
	FAIL("qtcad exact Obol snap should pass controller full-detail budget to source provider");
    }
    controller->setExactFullDetailBudget(0, 0);

    if (view.dmp()) {
	controller->setLodService(NULL);
	sourceService.stop();
	FAIL("qtcad source-backed exact Obol snap should not initialize the legacy display manager");
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
