/*    T E S T _ Q T C A D _ O B O L _ S E L E C T I O N _ S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BMeshShape.h"
#include "BObol/BLodService.h"
#include "BObol/BPerformance.h"
#include "BObol/BViewController.h"
#include "BObol/BVListShape.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "ged.h"
#include "ged/draw.h"
#include "ged/selection_state.h"
#include "QgObolDrawSyncPrivate.h"
#include "QgObolSelectionSyncPrivate.h"
#include "qtcad/QgView.h"
#include "raytrace.h"
#include "wdb.h"

#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoSeparator.h>

#include <Obol/cad/SoCADAssembly.h>

#include <QApplication>

#include <stdio.h>
#include <string.h>
#include <chrono>
#include <thread>
#include <vector>

static const int repeated_occurrence_count = 1024;

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

static int
make_selection_sync_db(const char *dbpath)
{
    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp)
	return 0;

    point_t bmin = {-1.0, -1.0, -1.0};
    point_t bmax = { 1.0,  1.0,  1.0};
    point_t center = {5.0, 0.0, 0.0};
    point_t extra_center = {0.0, 5.0, 0.0};

    int ret = mk_rpp(wdbp, "box.s", bmin, bmax) == 0 &&
	mk_rpp(wdbp, "tile.s", bmin, bmax) == 0 &&
	mk_sph(wdbp, "ball.s", center, 1.0) == 0 &&
	mk_sph(wdbp, "extra.s", extra_center, 0.5) == 0;
    struct wmember pair;
    BU_LIST_INIT(&pair.l);
    ret = ret &&
	mk_addmember("box.s", &pair.l, NULL, WMOP_UNION) != NULL &&
	mk_addmember("ball.s", &pair.l, NULL, WMOP_UNION) != NULL &&
	mk_lcomb(wdbp, "pair.c", &pair, 0, NULL, NULL, NULL, 0) == 0;

    struct wmember repeated;
    BU_LIST_INIT(&repeated.l);
    for (int i = 0; ret && i < repeated_occurrence_count; i++) {
	mat_t transform;
	MAT_IDN(transform);
	MAT_DELTAS(transform,
	    static_cast<fastf_t>((i % 32) * 4),
	    static_cast<fastf_t>((i / 32) * 4), 0.0);
	ret = mk_addmember("tile.s", &repeated.l, transform,
	    WMOP_UNION) != NULL;
    }
    ret = ret && mk_lcomb(wdbp, "repeated.c", &repeated, 0, NULL,
	NULL, NULL, 0) == 0;
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

static SoBRLDatabaseSource *
find_source_in_node(SoNode *node, const char *path)
{
    if (!node || !path)
	return NULL;

    const char *normalized_path = path;
    while (*normalized_path == '/')
	normalized_path++;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoBRLDatabaseSource *source =
	    static_cast<SoBRLDatabaseSource *>(node);
	const char *source_path = source ? source->path.getValue().getString() : NULL;
	if (source_path) {
	    while (*source_path == '/')
		source_path++;
	    if (BU_STR_EQUAL(source_path, normalized_path))
		return source;
	}
    }

    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    SoBRLDatabaseSource *source = find_source_in_node(
		group->getChild(i), path);
	    if (source)
		return source;
	}
    }
    return NULL;
}

static int
compact_summary_for_path(SoBRLDatabaseSource *source, const char *path,
	BObolCompactInstanceHandle &handle,
	BObolCompactInstanceSummary &summary)
{
    if (!source || !path)
	return 0;
    for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	BObolCompactInstanceHandle candidate;
	BObolCompactInstanceSummary candidateSummary;
	if (!source->getCompactInstanceHandle(i, candidate) ||
	    !source->getCompactInstanceSummary(candidate, candidateSummary))
	    continue;
	const char *candidatePath = candidateSummary.path.getString();
	if (candidatePath && strstr(candidatePath, path)) {
	    handle = candidate;
	    summary = candidateSummary;
	    return 1;
	}
    }
    return 0;
}

static union tree *
find_combination_leaf(union tree *tp, const char *name)
{
    if (!tp || !name)
	return NULL;
    if (tp->tr_op == OP_DB_LEAF)
	return BU_STR_EQUAL(tp->tr_l.tl_name, name) ? tp : NULL;
    if (tp->tr_op != OP_UNION && tp->tr_op != OP_INTERSECT &&
	tp->tr_op != OP_SUBTRACT && tp->tr_op != OP_XOR)
	return NULL;
    union tree *found = find_combination_leaf(tp->tr_b.tb_left, name);
    return found ? found : find_combination_leaf(tp->tr_b.tb_right, name);
}

static union tree *
find_nth_combination_leaf(union tree *tp, const char *name, int target,
	int &current)
{
    if (!tp || !name || target < 0)
	return NULL;
    if (tp->tr_op == OP_DB_LEAF) {
	if (!BU_STR_EQUAL(tp->tr_l.tl_name, name))
	    return NULL;
	return current++ == target ? tp : NULL;
    }
    if (tp->tr_op != OP_UNION && tp->tr_op != OP_INTERSECT &&
	tp->tr_op != OP_SUBTRACT && tp->tr_op != OP_XOR)
	return NULL;
    union tree *found = find_nth_combination_leaf(tp->tr_b.tb_left, name,
	target, current);
    return found ? found : find_nth_combination_leaf(tp->tr_b.tb_right, name,
	target, current);
}

static SoBRLDatabaseSource *
find_source(BObolViewController *controller, const char *path)
{
    if (!controller)
	return NULL;
    SoBRLDatabaseSource *source = find_source_in_node(
	controller->getRenderSceneRoot(), path);
    return source ? source : find_source_in_node(controller->getSceneRoot(), path);
}

static int
source_has_selected_shapes(SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;

    if (source->hasCompactInstanceIndex()) {
	for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	    BObolCompactInstanceHandle handle;
	    BObolCompactInstanceSummary summary;
	    if (source->getCompactInstanceHandle(i, handle) &&
		source->getCompactInstanceSummary(handle, summary) &&
		summary.selected)
		return 1;
	}
	return 0;
    }

    for (int i = 0; i < source->getRealizedShapeCount(); i++) {
	SoBRLVListShape *shape = source->getRealizedShape(i);
	if (shape && shape->selected.getValue())
	    return 1;
    }
    for (int i = 0; i < source->getRealizedMeshCount(); i++) {
	SoBRLMeshShape *mesh = source->getRealizedMesh(i);
	if (mesh && mesh->selected.getValue())
	    return 1;
    }
    return 0;
}

static int
compact_selected_count(SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;
    int count = 0;
    for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	BObolCompactInstanceHandle handle;
	BObolCompactInstanceSummary summary;
	if (source->getCompactInstanceHandle(i, handle) &&
	    source->getCompactInstanceSummary(handle, summary) &&
	    summary.selected)
	    count++;
    }
    return count;
}

static int
compact_hidden_count(SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;
    int count = 0;
    for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	BObolCompactInstanceHandle handle;
	BObolCompactInstanceSummary summary;
	if (source->getCompactInstanceHandle(i, handle) &&
	    source->getCompactInstanceSummary(handle, summary) &&
	    !summary.visible)
	    count++;
    }
    return count;
}

static int
test_compact_lod_submission_batches(SoBRLDatabaseSource *source)
{
    if (!source || !source->getDatabase() ||
	source->getCompactInstanceCount() != 3)
	return 0;

    BObolLodService service;
    service.setQueueLimits(3, 3, 1);
    if (!service.start(2, FALSE))
	return 0;

    SoSeparator *root = new SoSeparator;
    root->ref();
    root->addChild(source);
    SoOrthographicCamera *camera = new SoOrthographicCamera;
    BObolViewController batchingController(root, camera);
    batchingController.setLodService(&service);
    const uint64_t generation = service.beginGeneration();
    int totalSubmitted = 0;
    int batchCount = 0;
    for (; batchCount < 4; batchCount++) {
	const int submitted = batchingController.submitLodRequests(&service,
	    generation, TRUE, 0);
	if (submitted != 2)
	    break;
	totalSubmitted += submitted;
	for (int i = 0; i < 400 && service.inFlightCount() != 0; i++)
	    std::this_thread::sleep_for(std::chrono::milliseconds(5));
	if (service.inFlightCount() != 0 ||
	    service.queuedResultCountForDiagnostics() == 0)
	    break;
	std::vector<BObolLodResult> drained;
	service.drainResults(drained);
	if (!batchingController.hasPendingLodSubmissions()) {
	    batchCount++;
	    break;
	}
	if (batchCount >= 2)
	    break;
    }
    batchingController.setLodService(NULL);
    service.stop();
    BObolLodService undersizedService;
    undersizedService.setQueueLimits(2, 2, 1);
    int undersizedRejected = 0;
    if (undersizedService.start(1, FALSE)) {
	SoOrthographicCamera *undersizedCamera = new SoOrthographicCamera;
	BObolViewController undersizedController(root, undersizedCamera);
	undersizedController.setLodService(&undersizedService);
	undersizedRejected =
	    undersizedController.submitLodRequests(&undersizedService,
		undersizedService.beginGeneration(), TRUE, 0) == 0 &&
	    strstr(undersizedController.getLastLodDiagnostics().getString(),
		"requires three result slots") != NULL;
	undersizedController.setLodService(NULL);
	undersizedService.stop();
    }
    root->unref();
    return totalSubmitted == 6 && batchCount == 3 && undersizedRejected ? 1 : 0;
}

static int
test_large_compact_lod_churn(SoBRLDatabaseSource *source)
{
    if (!source || !source->getDatabase() ||
	source->getCompactInstanceCount() <= 0 ||
	!source->hasRealizedMeshGeometry())
	return 0;

    BObolLodService service;
    service.setQueueLimits(96, 96, 96);
    if (!service.start(2, FALSE))
	return 0;

    SoSeparator *root = new SoSeparator;
    root->ref();
    root->addChild(source);
    SoOrthographicCamera *camera = new SoOrthographicCamera;
    camera->height = 80.0f;
    const uint32_t originalSourceRevision = source->sourceRevision.getValue();
    const uint32_t originalRealizedRevision =
	source->realizedSourceRevision.getValue();
    int passed = 0;
    {
	BObolViewController controller(root, camera);
	controller.setViewportSize(800, 600);
	controller.setLodService(&service);
	bu_setenv("BOBOL_LOD_AABB_TASK_DELAY_MS", "20", 1);
	bu_setenv("BOBOL_LOD_OBB_TASK_DELAY_MS", "20", 1);

	const int initialSubmitted = controller.submitLodRequestsIfNeeded();
	const int initialVisited = controller.getLastLodVisitedMeshCount();
	camera->height = 90.0f;
	const int cameraSubmitted = controller.submitLodRequestsIfNeeded();
	const int cameraVisited = controller.getLastLodVisitedMeshCount();
	source->sourceRevision = originalSourceRevision + 1;
	source->realizedSourceRevision = originalSourceRevision + 1;
	source->realizationStatus = SoBRLDatabaseSource::REALIZED;
	source->stale = FALSE;
	source->staleReason = SoBRLDatabaseSource::STALE_NONE;
	const int editSubmitted = controller.submitLodRequestsIfNeeded();
	const int editVisited = controller.getLastLodVisitedMeshCount();

	/*
	 * Camera and edit churn must not wait for every OBB refinement to finish.
	 * Drive completion as a view frame does, and require a usable proxy result
	 * within one second while the bounded remainder continues asynchronously.
	 */
	int processed = 0;
	unsigned int applied = 0;
	unsigned int rejected = 0;
	for (int i = 0; i < 500; i++) {
	    processed += static_cast<int>(
		controller.processPendingLodResults(8, 2000));
	    applied += controller.getLastLodAppliedResultCount();
	    rejected += controller.getLastLodRejectedResultCount();
	    if (applied > 0)
		break;
	    std::this_thread::sleep_for(std::chrono::milliseconds(2));
	}
	const size_t inFlight = service.inFlightCount();
	const size_t queued = service.queuedResultCountForDiagnostics();
	passed = initialSubmitted > 0 && cameraSubmitted > 0 &&
	    editSubmitted > 0 && initialSubmitted <= 96 &&
	    cameraSubmitted <= 96 && editSubmitted <= 96 &&
	    initialVisited <= 32 && cameraVisited <= 32 && editVisited <= 32 &&
	    processed > 0 && applied > 0 && rejected == 0 &&
	    inFlight <= 96 && queued <= 96;
	if (!passed)
	    fprintf(stderr, "large LoD churn: submit=%d/%d/%d visit=%d/%d/%d "
		"processed=%d applied=%u rejected=%u inflight=%zu queued=%zu\n",
		initialSubmitted, cameraSubmitted, editSubmitted,
		initialVisited, cameraVisited, editVisited, processed, applied,
		rejected, inFlight, queued);
	controller.setLodService(NULL);
    }
    bu_setenv("BOBOL_LOD_AABB_TASK_DELAY_MS", "0", 1);
    bu_setenv("BOBOL_LOD_OBB_TASK_DELAY_MS", "0", 1);
    source->sourceRevision = originalSourceRevision;
    source->realizedSourceRevision = originalRealizedRevision;
    service.stop();
    root->unref();
    return passed;
}

static int
test_compact_lod_scale(struct db_i *dbip)
{
    const char *runSlow = getenv("BRLCAD_RUN_SLOW_TESTS");
    if (!dbip)
	return 0;

    const int candidateCount = runSlow && runSlow[0] &&
	!BU_STR_EQUAL(runSlow, "0") ? 100000 : 16;
    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = "stress.c";
    source->drawMode = SoBRLDatabaseSource::SHADED;
    source->representationMode = SoBRLDatabaseSource::REPRESENTATION_SHADED;
    source->sourceRevision = 1;

    std::shared_ptr<Obol::PartGeometry> geometry =
	std::make_shared<Obol::PartGeometry>();
    geometry->shaded.emplace();
    geometry->shaded->positions = {
	SbVec3f(-1.0f, -1.0f, 0.0f),
	SbVec3f(1.0f, -1.0f, 0.0f),
	SbVec3f(0.0f, 1.0f, 0.0f)
    };
    geometry->shaded->indices = {0, 1, 2};
    geometry->shaded->bounds = SbBox3f(SbVec3f(-1.0f, -1.0f, 0.0f),
	SbVec3f(1.0f, 1.0f, 0.0f));

    std::vector<BObolCompactOccurrence> occurrences;
    occurrences.reserve(candidateCount);
    for (int i = 0; i < candidateCount; i++) {
	BObolCompactOccurrence occurrence;
	occurrence.geometry = geometry;
	occurrence.summary.valid = TRUE;
	occurrence.summary.path.sprintf("stress.c/tile.s@%d", i);
	occurrence.summary.sourceName = "tile.s";
	occurrence.summary.sourceType = "rpp";
	occurrence.summary.geometryKind = "mesh";
	occurrence.summary.shapeKind =
	    BObolRealizedShapeSummary::SHAPE_MESH;
	occurrence.summary.visible = TRUE;
	occurrence.summary.selectable = TRUE;
	occurrence.occurrenceIndex = static_cast<uint32_t>(i);
	occurrence.localTransform.setTranslate(SbVec3f(
	    static_cast<float>(i % 1000) * 3.0f,
	    static_cast<float>(i / 1000) * 3.0f, 0.0f));
	occurrences.push_back(std::move(occurrence));
    }
    const int installed = source->setCompactOccurrenceRegistry(occurrences);
    occurrences.clear();
    occurrences.shrink_to_fit();
    source->realizedSourceRevision = source->sourceRevision.getValue();
    source->realizedInputsRevision = source->inputsRevision.getValue();
    source->realizationStatus = SoBRLDatabaseSource::REALIZED;
    source->stale = FALSE;
    source->staleReason = SoBRLDatabaseSource::STALE_NONE;
    if (installed != candidateCount) {
	source->ref();
	source->unref();
	return 0;
    }
    return test_large_compact_lod_churn(source);
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    bu_setenv("LIBRT_USE_COMB_INSTANCE_SPECIFIERS", "1", 1);

    QApplication app(argc, argv);

    const char *dbpath = "qtcad_obol_selection_sync_tmp.g";
    if (!make_selection_sync_db(dbpath))
	FAIL("failed to create qtcad Obol selection-sync test database");

    struct ged *gedp = ged_open("db", dbpath, 1);
    if (!gedp)
	FAIL("failed to open qtcad Obol selection-sync test database");

    QgView view(NULL, QgViewType::SW);
    view.resize(180, 140);
	ged_view_active_ctx_set(gedp, view.viewContext());
	(void)ged_view_context_host_attach(gedp, view.viewContext());

    BObolViewController *controller = view.obolViewController();
    if (!controller)
	FAIL("QgView should expose an Obol controller");
    controller->clearDatabaseSources();

    struct ged_draw_transaction draw_box =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "box.s");
    draw_box.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &draw_box))
	FAIL("GED wire draw should sync box source into Obol");

    struct ged_draw_appearance_settings shaded_appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    shaded_appearance.draw_mode = GED_DRAW_MODE_SHADED;
    struct ged_draw_transaction draw_ball =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "ball.s");
    draw_ball.view = view.viewContext();
    draw_ball.appearance = &shaded_appearance;
    int drew_ball = apply_and_sync(gedp, &view, &draw_ball);
    if (!drew_ball)
	FAIL("GED shaded draw should sync ball source into Obol");

    SoBRLDatabaseSource *box = find_source(controller, "box.s");
    SoBRLDatabaseSource *ball = find_source(controller, "ball.s");
    if (!box || box->getCompactInstanceCount() <= 0)
	FAIL("box should have realized Obol wire geometry");
    if (!ball || ball->getCompactInstanceCount() <= 0)
	FAIL("ball should have realized Obol mesh geometry");

    if (qg_obol_sync_selection_state(gedp, &view, nullptr) != 0)
	FAIL("empty GED selection should not change unselected Obol geometry");
    if (qg_obol_sync_selection_state_if_active(gedp, &view, nullptr) != 0)
	FAIL("draw-triggered selection sync should skip an empty selection");
    if (source_has_selected_shapes(box) || source_has_selected_shapes(ball))
	FAIL("fresh Obol geometry should start unselected");

    if (!ged_selection_select_path(gedp, nullptr, "box.s", 1))
	FAIL("GED should select box path");
    if (!qg_obol_sync_selection_state(gedp, &view, nullptr))
	FAIL("box selection should update Obol selection fields");
    if (!source_has_selected_shapes(box) || source_has_selected_shapes(ball))
	FAIL("Obol selection should select only box geometry");
    if (qg_obol_sync_selection_state(gedp, &view, nullptr) != 0)
	FAIL("repeating selection sync should be stable");

    if (!ged_selection_select_path(gedp, nullptr, "ball.s", 1))
	FAIL("GED should select ball path");
    if (!qg_obol_sync_selection_state(gedp, &view, nullptr))
	FAIL("ball selection should update Obol mesh selection fields");
    if (!source_has_selected_shapes(box) || !source_has_selected_shapes(ball))
	FAIL("Obol selection should retain box and add ball geometry");

    if (!ged_selection_clear(gedp, nullptr))
	FAIL("GED should clear default selection");
    ged_selection_recompute(gedp, nullptr);
    if (!qg_obol_sync_selection_state(gedp, &view, nullptr))
	FAIL("clearing GED selection should clear Obol selection fields");
    if (source_has_selected_shapes(box) || source_has_selected_shapes(ball))
	FAIL("Obol selection fields should clear with GED selection");

    controller->clearDatabaseSources();
    struct ged_draw_transaction draw_pair =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "pair.c");
    draw_pair.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &draw_pair))
	FAIL("GED wire draw should publish a compact combination root");
    SoBRLDatabaseSource *pair = find_source(controller, "pair.c");
    if (!pair || !pair->hasCompactInstanceIndex() ||
	pair->getCompactInstanceCountForPath("pair.c", TRUE) != 2)
	FAIL("pair should retain two addressable compact occurrences");
    if (pair->getRealizedShapeCount() != 0 ||
	pair->getRealizedMeshCount() != 0 ||
	pair->getRealizedShapeSummaryCount() != 2)
	FAIL("aggregate occurrences should expose summaries without Coin shape carriers");
    BObolCompactInstanceHandle box_handle;
    BObolCompactInstanceHandle ball_handle;
    BObolCompactInstanceSummary box_initial;
    BObolCompactInstanceSummary ball_initial;
    if (!compact_summary_for_path(pair, "box.s", box_handle, box_initial) ||
	!compact_summary_for_path(pair, "ball.s", ball_handle, ball_initial))
	FAIL("aggregate revision test should resolve initial occurrences");
    if (!pair->hasRealizedWireGeometry() || pair->hasRealizedMeshGeometry())
	FAIL("aggregate geometry presence should not depend on Coin shape counts");
    const SbColor compact_metadata_color(0.2f, 0.3f, 0.4f);
    if (pair->setCompactInstanceMetadataForPath("box.s", FALSE,
	17, 3, 9, 75, TRUE, compact_metadata_color, SbString("plastic")) != 1)
	FAIL("leaf-addressed compact metadata should update one occurrence");
    BObolCompactInstanceSummary box_metadata;
    BObolCompactInstanceSummary ball_metadata;
    if (!pair->getCompactInstanceSummary(box_handle, box_metadata) ||
	!pair->getCompactInstanceSummary(ball_handle, ball_metadata) ||
	box_metadata.regionId != 17 || box_metadata.airCode != 3 ||
	box_metadata.materialId != 9 || box_metadata.los != 75 ||
	!box_metadata.materialColorValid ||
	box_metadata.materialColor != compact_metadata_color ||
	bu_strcmp(box_metadata.materialShader.getString(), "plastic") != 0 ||
	box_metadata.appearanceRevision <= box_initial.appearanceRevision ||
	box_metadata.geometryRevision != box_initial.geometryRevision ||
	ball_metadata.appearanceRevision != ball_initial.appearanceRevision ||
	ball_metadata.regionId != ball_initial.regionId)
	FAIL("compact metadata should remain occurrence-local and geometry-neutral");
    if (pair->setCompactInstanceMetadataForPath("box.s", FALSE,
	box_initial.regionId, box_initial.airCode, box_initial.materialId,
	box_initial.los, box_initial.materialColorValid,
	box_initial.materialColor, box_initial.materialShader) != 1)
	FAIL("compact metadata test should restore the initial occurrence state");
    if (!pair->getCompactInstanceSummary(box_handle, box_initial))
	FAIL("compact metadata restore should refresh the revision baseline");

    SoSeparator *independent_root = new SoSeparator;
    SoBRLSceneController independent_scene(independent_root);
    independent_scene.shareRealizationRepository(
	controller->getSceneController());
    BObolDatabaseSourceSummary pair_source_summary;
    if (!pair->getSummary(pair_source_summary) ||
	independent_scene.replaceDatabaseSourceInstanceRepresentation(
	    pair_source_summary.instanceKey.getString(),
	    pair_source_summary.path.getString(),
	    pair_source_summary.representationKey.getString(),
	    pair_source_summary.representationMode, gedp->dbip,
	    pair_source_summary.drawMode,
	    pair_source_summary.sourceRevision) < 0)
	FAIL("independent view should publish the shared aggregate source");
    bobol_performance_counters_set_enabled(1);
    bobol_performance_counters_reset();
    independent_scene.realizePending();
    struct BObolPerformanceCounters independent_counters;
    bobol_performance_counters_get(&independent_counters);
    bobol_performance_counters_set_enabled(0);
    SoBRLDatabaseSource *independent_pair =
	independent_scene.findDatabaseSourceInstance(
	    pair_source_summary.instanceKey.getString());
    BObolCompactInstanceHandle independent_box_handle;
    BObolCompactInstanceHandle independent_ball_handle;
    BObolCompactInstanceSummary independent_box;
    BObolCompactInstanceSummary independent_ball;
    if (!independent_pair || !independent_pair->isCompactOccurrenceRegistry() ||
	!compact_summary_for_path(independent_pair, "box.s",
	    independent_box_handle, independent_box) ||
	!compact_summary_for_path(independent_pair, "ball.s",
	    independent_ball_handle, independent_ball) ||
	independent_box.geometryIdentity != box_initial.geometryIdentity ||
	independent_ball.geometryIdentity != ball_initial.geometryIdentity ||
	independent_counters.wire_cache_hits < 2 ||
	independent_counters.plot_calls != 0)
	FAIL("independent view should reuse the shared geometry repository");

    if (independent_scene.setDatabaseSourceInstanceRealizationViewPolicy(
	    pair_source_summary.instanceKey.getString(), TRUE, FALSE, FALSE,
	    2.0f, 0.5f, 320, 240, 0, 1.0f, 1.0f) <= 0)
	FAIL("independent view should accept a distinct realization policy");
    independent_scene.realizePending();
    BObolCompactInstanceSummary primary_box_after_policy;
    BObolCompactInstanceSummary independent_box_after_policy;
    if (!pair->getCompactInstanceSummary(box_handle,
	    primary_box_after_policy) ||
	!independent_pair->getCompactInstanceSummary(independent_box_handle,
	    independent_box_after_policy) ||
	primary_box_after_policy.geometryIdentity != box_initial.geometryIdentity ||
	primary_box_after_policy.geometryRevision != box_initial.geometryRevision ||
	independent_box_after_policy.geometryIdentity !=
	    independent_box.geometryIdentity)
	FAIL("independent view policy should not mutate shared occurrence geometry");

    struct ged_draw_transaction draw_repeated =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "repeated.c");
    draw_repeated.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &draw_repeated))
	FAIL("large repeated draw should publish a compact occurrence registry");
    SoBRLDatabaseSource *repeated = find_source(controller, "repeated.c");
    if (!repeated || !repeated->isCompactOccurrenceRegistry() ||
	repeated->getCompactInstanceCount() != repeated_occurrence_count ||
	repeated->getRealizedShapeCount() != 0 ||
	repeated->getRealizedMeshCount() != 0)
	FAIL("large repeated draw should remain carrier-free and addressable");

    std::vector<BObolCompactInstanceHandle> repeated_handles;
    std::vector<BObolCompactInstanceSummary> repeated_before;
    repeated_handles.reserve(repeated_occurrence_count);
    repeated_before.reserve(repeated_occurrence_count);
    uint64_t repeated_geometry_identity = 0;
    for (int i = 0; i < repeated_occurrence_count; i++) {
	BObolCompactInstanceHandle handle;
	BObolCompactInstanceSummary summary;
	if (!repeated->getCompactInstanceHandle(i, handle) ||
	    !repeated->getCompactInstanceSummary(handle, summary) ||
	    !summary.wireGeometry || summary.meshGeometry ||
	    (repeated_geometry_identity &&
	     summary.geometryIdentity != repeated_geometry_identity))
	    FAIL("repeated occurrences should share one immutable wire part");
	repeated_geometry_identity = summary.geometryIdentity;
	repeated_handles.push_back(handle);
	repeated_before.push_back(summary);
    }

    BObolDatabaseSourceSummary repeated_source_summary;
    if (!repeated->getSummary(repeated_source_summary))
	FAIL("large repeated source should expose publication metadata");
    SoSeparator *third_root = new SoSeparator;
    SoBRLSceneController third_scene(third_root);
    third_scene.shareRealizationRepository(controller->getSceneController());
    if (independent_scene.replaceDatabaseSourceInstanceRepresentation(
	    repeated_source_summary.instanceKey.getString(),
	    repeated_source_summary.path.getString(),
	    repeated_source_summary.representationKey.getString(),
	    repeated_source_summary.representationMode, gedp->dbip,
	    repeated_source_summary.drawMode,
	    repeated_source_summary.sourceRevision) < 0 ||
	third_scene.replaceDatabaseSourceInstanceRepresentation(
	    repeated_source_summary.instanceKey.getString(),
	    repeated_source_summary.path.getString(),
	    repeated_source_summary.representationKey.getString(),
	    repeated_source_summary.representationMode, gedp->dbip,
	    repeated_source_summary.drawMode,
	    repeated_source_summary.sourceRevision) < 0)
	FAIL("additional views should publish the repeated aggregate source");
    if (independent_scene.setDatabaseSourceInstanceRealizationViewPolicy(
	    repeated_source_summary.instanceKey.getString(), TRUE, FALSE, FALSE,
	    1.5f, 0.25f, 640, 480, 0, 1.0f, 1.0f) <= 0 ||
	third_scene.setDatabaseSourceInstanceRealizationViewPolicy(
	    repeated_source_summary.instanceKey.getString(), TRUE, FALSE, FALSE,
	    0.75f, 1.0f, 1024, 768, 0, 1.0f, 1.0f) <= 0)
	FAIL("three views should accept independent realization policies");
    bobol_performance_counters_set_enabled(1);
    bobol_performance_counters_reset();
    independent_scene.realizePending();
    third_scene.realizePending();
    struct BObolPerformanceCounters repeated_view_counters;
    bobol_performance_counters_get(&repeated_view_counters);
    bobol_performance_counters_set_enabled(0);
    SoBRLDatabaseSource *independent_repeated =
	independent_scene.findDatabaseSourceInstance(
	    repeated_source_summary.instanceKey.getString());
    SoBRLDatabaseSource *third_repeated = third_scene.findDatabaseSourceInstance(
	repeated_source_summary.instanceKey.getString());
    if (!independent_repeated || !third_repeated ||
	independent_repeated->getCompactInstanceCount() !=
	    repeated_occurrence_count ||
	third_repeated->getCompactInstanceCount() != repeated_occurrence_count ||
	repeated_view_counters.wire_cache_misses != 0 ||
	repeated_view_counters.plot_calls != 0) {
	fprintf(stderr,
	    "repeated view reuse: independent=%p count=%d third=%p count=%d "
	    "hits=%" PRIu64 " misses=%" PRIu64 " plots=%" PRIu64 "\n",
	    static_cast<void *>(independent_repeated),
	    independent_repeated ? independent_repeated->getCompactInstanceCount() : -1,
	    static_cast<void *>(third_repeated),
	    third_repeated ? third_repeated->getCompactInstanceCount() : -1,
	    repeated_view_counters.wire_cache_hits,
	    repeated_view_counters.wire_cache_misses,
	    repeated_view_counters.plot_calls);
	FAIL("additional view policies should reuse the repeated geometry repository");
    }
    std::vector<BObolCompactInstanceHandle> repeated_view_handles[2];
    std::vector<BObolCompactInstanceSummary> repeated_view_before[2];
    for (int scene_index = 0; scene_index < 2; scene_index++) {
	SoBRLDatabaseSource *scene_source = scene_index ? third_repeated :
	    independent_repeated;
	repeated_view_handles[scene_index].reserve(repeated_occurrence_count);
	repeated_view_before[scene_index].reserve(repeated_occurrence_count);
	for (int i = 0; i < repeated_occurrence_count; i++) {
	    BObolCompactInstanceHandle handle;
	    BObolCompactInstanceSummary summary;
	    if (!scene_source->getCompactInstanceHandle(i, handle) ||
		!scene_source->getCompactInstanceSummary(handle, summary) ||
		summary.geometryIdentity != repeated_geometry_identity)
		FAIL("independent policies should preserve shared occurrence geometry");
	    repeated_view_handles[scene_index].push_back(handle);
	    repeated_view_before[scene_index].push_back(summary);
	}
    }

    struct directory *repeated_dp = db_lookup(gedp->dbip, "repeated.c",
	LOOKUP_QUIET);
    struct rt_db_internal repeated_internal;
    RT_DB_INTERNAL_INIT(&repeated_internal);
    if (!repeated_dp || rt_db_get_internal(&repeated_internal, repeated_dp,
	    gedp->dbip, NULL) < 0 ||
	repeated_internal.idb_type != ID_COMBINATION)
	FAIL("repeated transform edit should read the combination");
    struct rt_comb_internal *repeated_comb =
	static_cast<struct rt_comb_internal *>(repeated_internal.idb_ptr);
    int repeated_leaf_index = 0;
    union tree *changed_repeated_leaf = find_nth_combination_leaf(
	repeated_comb->tree, "tile.s", 637, repeated_leaf_index);
    if (!changed_repeated_leaf)
	FAIL("repeated transform edit should find the targeted occurrence");
    if (!changed_repeated_leaf->tr_l.tl_mat)
	changed_repeated_leaf->tr_l.tl_mat = static_cast<matp_t>(
	    bu_calloc(1, sizeof(mat_t), "repeated occurrence transform"));
    MAT_IDN(changed_repeated_leaf->tr_l.tl_mat);
    MAT_DELTAS(changed_repeated_leaf->tr_l.tl_mat, 73.0, 19.0, 5.0);
    if (rt_db_put_internal(repeated_dp, gedp->dbip, &repeated_internal) < 0)
	FAIL("repeated transform edit should update the combination");
    struct ged_draw_transaction transform_repeated =
	ged_draw_transaction_make(GED_DRAW_TXN_STALE_SOURCE, "repeated.c");
    transform_repeated.view = view.viewContext();
    transform_repeated.stale_reason = GED_DRAW_STALE_SOURCE_CHANGED;
    bobol_performance_counters_set_enabled(1);
    bobol_performance_counters_reset();
    if (!apply_and_sync(gedp, &view, &transform_repeated))
	FAIL("repeated transform edit should reconcile the retained aggregate");
    struct BObolPerformanceCounters repeated_edit_counters;
    bobol_performance_counters_get(&repeated_edit_counters);
    bobol_performance_counters_set_enabled(0);
    repeated = find_source(controller, "repeated.c");
    int changed_placements = 0;
    if (!repeated ||
	repeated->getCompactInstanceCount() != repeated_occurrence_count)
	FAIL("repeated transform edit should preserve the occurrence count");
    for (int i = 0; i < repeated_occurrence_count; i++) {
	BObolCompactInstanceSummary summary;
	if (!repeated->isCompactInstanceHandleValid(repeated_handles[i]) ||
	    !repeated->getCompactInstanceSummary(repeated_handles[i], summary) ||
	    summary.geometryIdentity != repeated_before[i].geometryIdentity ||
	    summary.geometryRevision != repeated_before[i].geometryRevision ||
	    summary.appearanceRevision != repeated_before[i].appearanceRevision ||
	    summary.visibilityRevision != repeated_before[i].visibilityRevision ||
	    summary.selectionRevision != repeated_before[i].selectionRevision)
	    FAIL("repeated transform edit should preserve handles and non-placement state");
	if (summary.placementRevision != repeated_before[i].placementRevision) {
	    if (summary.placementRevision <=
		repeated_before[i].placementRevision)
		FAIL("repeated placement revisions must advance monotonically");
	    changed_placements++;
	}
    }
    if (changed_placements != 1 || repeated_edit_counters.plot_calls != 0 ||
	repeated_edit_counters.wire_cache_misses != 0 ||
	repeated_edit_counters.source_replace_calls != 0 ||
	repeated_edit_counters.sources_visited != 0 ||
	repeated_edit_counters.sources_realized != 0 ||
	repeated_edit_counters.cad_compact_sources != 1)
	FAIL("localized repeated transform should change one placement without geometry work");

    BObolDatabaseSourceSummary repeated_after_summary;
    if (!repeated->getSummary(repeated_after_summary) ||
	independent_scene.replaceDatabaseSourceInstanceRepresentation(
	    repeated_after_summary.instanceKey.getString(),
	    repeated_after_summary.path.getString(),
	    repeated_after_summary.representationKey.getString(),
	    repeated_after_summary.representationMode, gedp->dbip,
	    repeated_after_summary.drawMode,
	    repeated_after_summary.sourceRevision) < 0 ||
	third_scene.replaceDatabaseSourceInstanceRepresentation(
	    repeated_after_summary.instanceKey.getString(),
	    repeated_after_summary.path.getString(),
	    repeated_after_summary.representationKey.getString(),
	    repeated_after_summary.representationMode, gedp->dbip,
	    repeated_after_summary.drawMode,
	    repeated_after_summary.sourceRevision) < 0)
	FAIL("database edit should republish the repeated source to all views");
    bobol_performance_counters_set_enabled(1);
    bobol_performance_counters_reset();
    independent_scene.realizePending();
    third_scene.realizePending();
    struct BObolPerformanceCounters repeated_republish_counters;
    bobol_performance_counters_get(&repeated_republish_counters);
    bobol_performance_counters_set_enabled(0);
    independent_repeated = independent_scene.findDatabaseSourceInstance(
	repeated_after_summary.instanceKey.getString());
    third_repeated = third_scene.findDatabaseSourceInstance(
	repeated_after_summary.instanceKey.getString());
    if (!independent_repeated || !third_repeated ||
	repeated_republish_counters.plot_calls != 0 ||
	repeated_republish_counters.wire_cache_misses != 0)
	FAIL("edited repeated source should reuse geometry in every view");
    for (int scene_index = 0; scene_index < 2; scene_index++) {
	SoBRLDatabaseSource *scene_source = scene_index ? third_repeated :
	    independent_repeated;
	int scene_changed_placements = 0;
	for (int i = 0; i < repeated_occurrence_count; i++) {
	    BObolCompactInstanceSummary summary;
	    if (!scene_source->isCompactInstanceHandleValid(
		    repeated_view_handles[scene_index][i]) ||
		!scene_source->getCompactInstanceSummary(
		    repeated_view_handles[scene_index][i], summary) ||
		summary.geometryIdentity != repeated_geometry_identity ||
		summary.geometryRevision !=
		    repeated_view_before[scene_index][i].geometryRevision ||
		summary.appearanceRevision !=
		    repeated_view_before[scene_index][i].appearanceRevision ||
		summary.visibilityRevision !=
		    repeated_view_before[scene_index][i].visibilityRevision ||
		summary.selectionRevision !=
		    repeated_view_before[scene_index][i].selectionRevision)
		FAIL("cross-view edit should preserve handles and geometry state");
	    if (summary.placementRevision !=
		repeated_view_before[scene_index][i].placementRevision)
		scene_changed_placements++;
	}
	if (scene_changed_placements != 1)
	    FAIL("cross-view edit should advance exactly one placement per scene");
    }

    struct ged_draw_transaction erase_repeated =
	ged_draw_transaction_make(GED_DRAW_TXN_ERASE, "repeated.c");
    erase_repeated.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &erase_repeated) ||
	find_source(controller, "repeated.c"))
	FAIL("top-level erase should release the primary repeated source");
    struct ged_draw_transaction redraw_shared_repeated =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "repeated.c");
    redraw_shared_repeated.view = view.viewContext();
    bobol_performance_counters_set_enabled(1);
    bobol_performance_counters_reset();
    if (!apply_and_sync(gedp, &view, &redraw_shared_repeated))
	FAIL("shared-resident repeated source should redraw");
    struct BObolPerformanceCounters shared_residency_counters;
    bobol_performance_counters_get(&shared_residency_counters);
    bobol_performance_counters_set_enabled(0);
    if (shared_residency_counters.plot_calls != 0 ||
	shared_residency_counters.wire_cache_misses != 0)
	FAIL("remaining view owners should retain shared repository geometry");

    struct ged_draw_transaction erase_repeated_again =
	ged_draw_transaction_make(GED_DRAW_TXN_ERASE, "repeated.c");
    erase_repeated_again.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &erase_repeated_again) ||
	independent_scene.removeDatabaseSourceInstance(
	    repeated_after_summary.instanceKey.getString()) <= 0 ||
	third_scene.removeDatabaseSourceInstance(
	    repeated_after_summary.instanceKey.getString()) <= 0)
	FAIL("last-owner test should erase the repeated source from every view");
    struct ged_draw_transaction redraw_evicted_repeated =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "repeated.c");
    redraw_evicted_repeated.view = view.viewContext();
    bobol_performance_counters_set_enabled(1);
    bobol_performance_counters_reset();
    if (!apply_and_sync(gedp, &view, &redraw_evicted_repeated))
	FAIL("last-owner-evicted repeated source should redraw");
    struct BObolPerformanceCounters evicted_residency_counters;
    bobol_performance_counters_get(&evicted_residency_counters);
    bobol_performance_counters_set_enabled(0);
    if (evicted_residency_counters.plot_calls != 1 ||
	evicted_residency_counters.wire_cache_misses != 1)
	FAIL("last owner release should evict one shared repeated geometry payload");

    bobol_performance_counters_set_enabled(1);
    bobol_performance_counters_reset();
    view.aet(35.0, 20.0, 5.0);
    struct BObolPerformanceCounters camera_only_counters;
    bobol_performance_counters_get(&camera_only_counters);
    bobol_performance_counters_set_enabled(0);
    if (camera_only_counters.source_replace_calls != 0 ||
	camera_only_counters.sources_visited != 0 ||
	camera_only_counters.sources_realized != 0 ||
	camera_only_counters.cad_compact_sources != 0 ||
	camera_only_counters.plot_calls != 0)
	FAIL("camera-only updates should not publish or realize database sources");

    struct ged_draw_transaction erase_repeated_for_lod =
	ged_draw_transaction_make(GED_DRAW_TXN_ERASE, "repeated.c");
    erase_repeated_for_lod.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &erase_repeated_for_lod))
	FAIL("large compact LoD test should replace the wire representation");
    struct ged_draw_transaction draw_repeated_shaded =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "repeated.c");
    draw_repeated_shaded.view = view.viewContext();
    draw_repeated_shaded.appearance = &shaded_appearance;
    if (!apply_and_sync(gedp, &view, &draw_repeated_shaded))
	FAIL("large compact LoD test should publish a shaded aggregate");
    SoBRLDatabaseSource *repeated_shaded =
	find_source(controller, "repeated.c");
    if (!test_large_compact_lod_churn(repeated_shaded))
	FAIL("large compact LoD scheduling should remain bounded under camera and edit churn");
    if (!test_compact_lod_scale(gedp->dbip))
	FAIL("typed compact LoD scheduling should remain bounded under churn");

    if (!ged_selection_select_path(gedp, nullptr, "pair.c/box.s", 1) ||
	!qg_obol_sync_selection_state(gedp, &view, nullptr))
	FAIL("nested GED selection should update the aggregate occurrence");
    if (pair->getCompactInstanceCountForPath("pair.c/box.s", FALSE) != 1)
	FAIL("aggregate registry should resolve the selected nested path");
    int selected_count = compact_selected_count(pair);
    if (selected_count != 1)
	FAIL("nested selection should select exactly one aggregate occurrence");
    BObolCompactInstanceSummary box_selected;
    BObolCompactInstanceSummary ball_unselected;
    if (!pair->getCompactInstanceSummary(box_handle, box_selected) ||
	!pair->getCompactInstanceSummary(ball_handle, ball_unselected) ||
	box_selected.selectionRevision <= box_initial.selectionRevision ||
	box_selected.geometryRevision != box_initial.geometryRevision ||
	box_selected.appearanceRevision != box_initial.appearanceRevision ||
	box_selected.placementRevision != box_initial.placementRevision ||
	box_selected.visibilityRevision != box_initial.visibilityRevision ||
	ball_unselected.selectionRevision != ball_initial.selectionRevision)
	FAIL("nested selection should advance only the selected occurrence's selection revision");

    struct ged_draw_transaction hide_box =
	ged_draw_transaction_make_value(GED_DRAW_TXN_VISIBILITY,
	    "pair.c/box.s", 0.0);
    hide_box.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &hide_box))
	FAIL("nested visibility transaction should address aggregate geometry");
    int hidden_count = compact_hidden_count(pair);
    if (hidden_count != 1)
	FAIL("nested visibility should hide exactly one aggregate occurrence");
    BObolCompactInstanceSummary box_hidden;
	if (!pair->getCompactInstanceSummary(box_handle, box_hidden))
	FAIL("nested visibility should retain the occurrence handle");
    if (box_hidden.visibilityRevision <= box_selected.visibilityRevision ||
	box_hidden.geometryRevision != box_selected.geometryRevision ||
	box_hidden.appearanceRevision != box_selected.appearanceRevision ||
	box_hidden.placementRevision != box_selected.placementRevision ||
	box_hidden.selectionRevision != box_selected.selectionRevision) {
	printf("visibility revisions before=%llu/%llu/%llu/%llu/%llu after=%llu/%llu/%llu/%llu/%llu\n",
	    (unsigned long long)box_selected.geometryRevision,
	    (unsigned long long)box_selected.appearanceRevision,
	    (unsigned long long)box_selected.placementRevision,
	    (unsigned long long)box_selected.visibilityRevision,
	    (unsigned long long)box_selected.selectionRevision,
	    (unsigned long long)box_hidden.geometryRevision,
	    (unsigned long long)box_hidden.appearanceRevision,
	    (unsigned long long)box_hidden.placementRevision,
	    (unsigned long long)box_hidden.visibilityRevision,
	    (unsigned long long)box_hidden.selectionRevision);
	FAIL("nested visibility should advance only the occurrence visibility revision");
    }

    struct ged_draw_transaction show_box =
	ged_draw_transaction_make_value(GED_DRAW_TXN_VISIBILITY,
	    "pair.c/box.s", 1.0);
    show_box.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &show_box))
	FAIL("nested visibility restore should address aggregate geometry");

    struct ged_draw_transaction erase_box =
	ged_draw_transaction_make(GED_DRAW_TXN_ERASE, "pair.c/box.s");
    erase_box.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &erase_box))
	FAIL("nested erase should address aggregate geometry");
    hidden_count = compact_hidden_count(pair);
    if (hidden_count != 1)
	FAIL("nested erase should hide exactly one aggregate occurrence");

    struct ged_draw_transaction redraw_box =
	ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, "pair.c/box.s");
    redraw_box.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &redraw_box))
	FAIL("nested redraw should restore aggregate geometry");
    hidden_count = compact_hidden_count(pair);
    if (hidden_count != 0)
	FAIL("nested redraw should restore the erased aggregate occurrence");

    BObolCompactInstanceSummary box_before;
    BObolCompactInstanceSummary ball_before;
    if (!compact_summary_for_path(pair, "box.s", box_handle, box_before) ||
	!compact_summary_for_path(pair, "ball.s", ball_handle, ball_before))
	FAIL("aggregate refresh test should resolve stable occurrence handles");
    struct ged_draw_transaction stale_pair_style =
	ged_draw_transaction_make(GED_DRAW_TXN_STALE_SOURCE, "pair.c");
    stale_pair_style.view = view.viewContext();
    stale_pair_style.stale_reason = GED_DRAW_STALE_SETTINGS_CHANGED;
    bobol_performance_counters_set_enabled(1);
    bobol_performance_counters_reset();
    if (!apply_and_sync(gedp, &view, &stale_pair_style))
	FAIL("style revision should republish aggregate instance state");
    struct BObolPerformanceCounters style_counters;
    bobol_performance_counters_get(&style_counters);
    bobol_performance_counters_set_enabled(0);
    pair = find_source(controller, "pair.c");
    BObolCompactInstanceSummary box_style_after;
    BObolCompactInstanceSummary ball_style_after;
    bool style_state_valid = pair &&
	pair->isCompactInstanceHandleValid(box_handle) &&
	pair->isCompactInstanceHandleValid(ball_handle) &&
	pair->getCompactInstanceSummary(box_handle, box_style_after) &&
	pair->getCompactInstanceSummary(ball_handle, ball_style_after);
    if (style_state_valid &&
	box_style_after.visibilityRevision != box_before.visibilityRevision)
	printf("style visibility revision before=%llu after=%llu\n",
	    (unsigned long long)box_before.visibilityRevision,
	    (unsigned long long)box_style_after.visibilityRevision);
    if (!style_state_valid ||
	box_style_after.geometryRevision != box_before.geometryRevision ||
	box_style_after.appearanceRevision != box_before.appearanceRevision ||
	box_style_after.placementRevision != box_before.placementRevision ||
	box_style_after.visibilityRevision != box_before.visibilityRevision ||
	box_style_after.selectionRevision != box_before.selectionRevision ||
	ball_style_after.geometryRevision != ball_before.geometryRevision ||
	ball_style_after.appearanceRevision != ball_before.appearanceRevision ||
	ball_style_after.placementRevision != ball_before.placementRevision ||
	ball_style_after.visibilityRevision != ball_before.visibilityRevision ||
	ball_style_after.selectionRevision != ball_before.selectionRevision ||
	style_counters.wire_cache_hits < 2 || style_counters.plot_calls != 0)
	FAIL("no-op style revision should preserve all retained occurrence channels");
    box_before = box_style_after;
    ball_before = ball_style_after;

    struct directory *box_dp = db_lookup(gedp->dbip, "box.s", LOOKUP_QUIET);
    struct rt_db_internal box_internal;
    RT_DB_INTERNAL_INIT(&box_internal);
    if (!box_dp || rt_db_get_internal(&box_internal, box_dp, gedp->dbip,
	    NULL) < 0 || box_internal.idb_type != ID_ARB8)
	FAIL("aggregate refresh test should read the box primitive");
    struct rt_arb_internal *arb =
	static_cast<struct rt_arb_internal *>(box_internal.idb_ptr);
    for (int i = 0; i < 8; i++)
	arb->pt[i][X] *= 2.0;
    if (rt_db_put_internal(box_dp, gedp->dbip, &box_internal) < 0)
	FAIL("aggregate refresh test should update the box primitive");

    struct ged_draw_transaction stale_box =
	ged_draw_transaction_make(GED_DRAW_TXN_STALE_SOURCE, "box.s");
    stale_box.view = view.viewContext();
    stale_box.stale_reason = GED_DRAW_STALE_SOURCE_CHANGED;
    if (!apply_and_sync(gedp, &view, &stale_box))
	FAIL("primitive edit should refresh the matching compact part");
    BObolCompactInstanceSummary box_after;
    BObolCompactInstanceSummary ball_after;
    if (!pair->isCompactInstanceHandleValid(box_handle) ||
	!pair->isCompactInstanceHandleValid(ball_handle) ||
	!pair->getCompactInstanceSummary(box_handle, box_after) ||
	!pair->getCompactInstanceSummary(ball_handle, ball_after))
	FAIL("primitive edit should preserve occurrence handles");
    const bool primitive_revision_mismatch =
	box_after.geometryRevision <= box_before.geometryRevision ||
	box_after.appearanceRevision != box_before.appearanceRevision ||
	box_after.placementRevision != box_before.placementRevision ||
	box_after.visibilityRevision != box_before.visibilityRevision ||
	box_after.selectionRevision != box_before.selectionRevision ||
	ball_after.geometryRevision != ball_before.geometryRevision ||
	ball_after.appearanceRevision != ball_before.appearanceRevision ||
	ball_after.placementRevision != ball_before.placementRevision ||
	ball_after.visibilityRevision != ball_before.visibilityRevision ||
	ball_after.selectionRevision != ball_before.selectionRevision;
    if (primitive_revision_mismatch)
	printf("primitive revisions box=%llu/%llu/%llu/%llu/%llu -> %llu/%llu/%llu/%llu/%llu ball=%llu/%llu/%llu/%llu/%llu -> %llu/%llu/%llu/%llu/%llu\n",
	    (unsigned long long)box_before.geometryRevision,
	    (unsigned long long)box_before.appearanceRevision,
	    (unsigned long long)box_before.placementRevision,
	    (unsigned long long)box_before.visibilityRevision,
	    (unsigned long long)box_before.selectionRevision,
	    (unsigned long long)box_after.geometryRevision,
	    (unsigned long long)box_after.appearanceRevision,
	    (unsigned long long)box_after.placementRevision,
	    (unsigned long long)box_after.visibilityRevision,
	    (unsigned long long)box_after.selectionRevision,
	    (unsigned long long)ball_before.geometryRevision,
	    (unsigned long long)ball_before.appearanceRevision,
	    (unsigned long long)ball_before.placementRevision,
	    (unsigned long long)ball_before.visibilityRevision,
	    (unsigned long long)ball_before.selectionRevision,
	    (unsigned long long)ball_after.geometryRevision,
	    (unsigned long long)ball_after.appearanceRevision,
	    (unsigned long long)ball_after.placementRevision,
	    (unsigned long long)ball_after.visibilityRevision,
	    (unsigned long long)ball_after.selectionRevision);
    if (pair->getCompactInstanceCount() != 2 ||
	pair->stale.getValue() ||
	box_after.geometryIdentity == box_before.geometryIdentity ||
	primitive_revision_mismatch ||
	ball_after.geometryIdentity != ball_before.geometryIdentity ||
	box_after.localBounds.isEmpty() ||
	box_after.localBounds.getMin()[X] > -1.9f ||
	box_after.localBounds.getMax()[X] < 1.9f)
	FAIL("primitive edit should replace only its shared part and preserve occurrence handles");

    struct directory *pair_dp = db_lookup(gedp->dbip, "pair.c", LOOKUP_QUIET);
    struct rt_db_internal pair_internal;
    RT_DB_INTERNAL_INIT(&pair_internal);
    if (!pair_dp || rt_db_get_internal(&pair_internal, pair_dp, gedp->dbip,
	    NULL) < 0 || pair_internal.idb_type != ID_COMBINATION)
	FAIL("aggregate structural refresh should read the combination");
    struct rt_comb_internal *comb =
	static_cast<struct rt_comb_internal *>(pair_internal.idb_ptr);
    union tree *extra_leaf = NULL;
    union tree *new_root = NULL;
    BU_GET(extra_leaf, union tree);
    RT_TREE_INIT(extra_leaf);
    extra_leaf->tr_l.tl_op = OP_DB_LEAF;
    extra_leaf->tr_l.tl_name = bu_strdup("extra.s");
    extra_leaf->tr_l.tl_mat = NULL;
    BU_GET(new_root, union tree);
    RT_TREE_INIT(new_root);
    new_root->tr_b.tb_op = OP_UNION;
    new_root->tr_b.tb_left = comb->tree;
    new_root->tr_b.tb_right = extra_leaf;
    comb->tree = new_root;
    if (rt_db_put_internal(pair_dp, gedp->dbip, &pair_internal) < 0)
	FAIL("aggregate structural refresh should append a combination member");
    struct ged_draw_transaction stale_pair =
	ged_draw_transaction_make(GED_DRAW_TXN_STALE_SOURCE, "pair.c");
    stale_pair.view = view.viewContext();
    stale_pair.stale_reason = GED_DRAW_STALE_SOURCE_CHANGED;
    bobol_performance_counters_set_enabled(1);
    bobol_performance_counters_reset();
    if (!apply_and_sync(gedp, &view, &stale_pair))
	FAIL("combination edit should rebuild the compact occurrence diff");
    struct BObolPerformanceCounters structural_counters;
    bobol_performance_counters_get(&structural_counters);
    bobol_performance_counters_set_enabled(0);
    if (structural_counters.wire_cache_hits < 2 ||
	structural_counters.plot_calls > 1) {
	fprintf(stderr, "structural cache counters: wire_hits=%" PRIu64
	    " plot_calls=%" PRIu64 "\n",
	    structural_counters.wire_cache_hits,
	    structural_counters.plot_calls);
	FAIL("structural diff should reuse unchanged retained wire parts");
    }
    pair = find_source(controller, "pair.c");
    BObolCompactInstanceHandle extra_handle;
    BObolCompactInstanceSummary extra_summary;
    if (!pair || !pair->isCompactOccurrenceRegistry() ||
	pair->getCompactInstanceCount() != 3 ||
	!pair->isCompactInstanceHandleValid(box_handle) ||
	!pair->isCompactInstanceHandleValid(ball_handle) ||
	!compact_summary_for_path(pair, "extra.s", extra_handle, extra_summary))
	FAIL("structural diff should preserve existing handles and add one occurrence");

    BObolCompactInstanceSummary box_structure_before;
    BObolCompactInstanceSummary ball_structure_before;
    BObolCompactInstanceSummary extra_structure_before;
    if (!pair->getCompactInstanceSummary(box_handle, box_structure_before) ||
	!pair->getCompactInstanceSummary(ball_handle, ball_structure_before) ||
	!pair->getCompactInstanceSummary(extra_handle, extra_structure_before))
	FAIL("transform diff should capture the retained occurrence state");

    RT_DB_INTERNAL_INIT(&pair_internal);
    if (rt_db_get_internal(&pair_internal, pair_dp, gedp->dbip, NULL) < 0 ||
	pair_internal.idb_type != ID_COMBINATION)
	FAIL("transform diff should read the combination");
    comb = static_cast<struct rt_comb_internal *>(pair_internal.idb_ptr);
    union tree *box_leaf = find_combination_leaf(comb->tree, "box.s");
    if (!box_leaf)
	FAIL("transform diff should find the box occurrence");
    if (!box_leaf->tr_l.tl_mat)
	box_leaf->tr_l.tl_mat = static_cast<matp_t>(
	    bu_calloc(1, sizeof(mat_t), "combination occurrence transform"));
    MAT_IDN(box_leaf->tr_l.tl_mat);
    MAT_DELTAS(box_leaf->tr_l.tl_mat, 3.0, 0.0, 0.0);
    if (rt_db_put_internal(pair_dp, gedp->dbip, &pair_internal) < 0)
	FAIL("transform diff should update the combination occurrence");

    struct ged_draw_transaction transform_pair =
	ged_draw_transaction_make(GED_DRAW_TXN_STALE_SOURCE, "pair.c");
    transform_pair.view = view.viewContext();
    transform_pair.stale_reason = GED_DRAW_STALE_SOURCE_CHANGED;
    if (!apply_and_sync(gedp, &view, &transform_pair))
	FAIL("combination transform should reconcile retained occurrences");
    pair = find_source(controller, "pair.c");
    BObolCompactInstanceSummary box_structure_after;
    BObolCompactInstanceSummary ball_structure_after;
    BObolCompactInstanceSummary extra_structure_after;
    if (!pair || pair->getCompactInstanceCount() != 3 ||
	!pair->getCompactInstanceSummary(box_handle, box_structure_after) ||
	!pair->getCompactInstanceSummary(ball_handle, ball_structure_after) ||
	!pair->getCompactInstanceSummary(extra_handle, extra_structure_after) ||
	box_structure_after.geometryIdentity !=
	    box_structure_before.geometryIdentity ||
	ball_structure_after.geometryIdentity !=
	    ball_structure_before.geometryIdentity ||
	extra_structure_after.geometryIdentity !=
	    extra_structure_before.geometryIdentity ||
	box_structure_after.localToSource.equals(
	    box_structure_before.localToSource, 0.000001f) ||
	!ball_structure_after.localToSource.equals(
	    ball_structure_before.localToSource, 0.000001f) ||
	!extra_structure_after.localToSource.equals(
	    extra_structure_before.localToSource, 0.000001f) ||
	box_structure_after.placementRevision <=
	    box_structure_before.placementRevision ||
	box_structure_after.geometryRevision !=
	    box_structure_before.geometryRevision ||
	box_structure_after.appearanceRevision !=
	    box_structure_before.appearanceRevision ||
	box_structure_after.visibilityRevision !=
	    box_structure_before.visibilityRevision ||
	box_structure_after.selectionRevision !=
	    box_structure_before.selectionRevision ||
	ball_structure_after.placementRevision !=
	    ball_structure_before.placementRevision ||
	extra_structure_after.placementRevision !=
	    extra_structure_before.placementRevision ||
	box_structure_after.selected != box_structure_before.selected)
	FAIL("transform diff should update only placement and preserve parts, handles, and display state");

    (void)ged_selection_clear(gedp, nullptr);
    ged_selection_recompute(gedp, nullptr);
    controller->clearDatabaseSources();
    struct ged_draw_appearance_settings pair_shaded_appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    pair_shaded_appearance.draw_mode = GED_DRAW_MODE_SHADED;
    struct ged_draw_transaction draw_pair_shaded =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "pair.c");
    draw_pair_shaded.view = view.viewContext();
    draw_pair_shaded.appearance = &pair_shaded_appearance;
    if (!apply_and_sync(gedp, &view, &draw_pair_shaded))
	FAIL("GED shaded draw should publish a compact combination root");
    pair = find_source(controller, "pair.c");
    if (!pair || !pair->isCompactOccurrenceRegistry() ||
	pair->getCompactInstanceCountForPath("pair.c", TRUE) != 3)
	FAIL("shaded pair should retain three addressable mesh occurrences");
    for (int i = 0; i < pair->getCompactInstanceCount(); i++) {
	BObolCompactInstanceHandle shaded_handle;
	BObolCompactInstanceSummary shaded_summary;
	if (!pair->getCompactInstanceHandle(i, shaded_handle) ||
	    !pair->getCompactInstanceSummary(shaded_handle, shaded_summary) ||
	    !shaded_summary.meshGeometry || shaded_summary.wireGeometry)
	    FAIL("shaded cache channels must not reuse normal wire payloads");
    }
    if (!test_compact_lod_submission_batches(pair))
	FAIL("compact LoD submission should resume across bounded queue batches");
    if (!ged_selection_select_path(gedp, nullptr, "pair.c/ball.s", 1) ||
	!qg_obol_sync_selection_state(gedp, &view, nullptr))
	FAIL("nested shaded selection should update the aggregate occurrence");
    selected_count = compact_selected_count(pair);
    if (selected_count != 1)
	FAIL("nested shaded selection should select exactly one mesh occurrence");

    BObolCompactInstanceHandle shaded_extra_handle;
    BObolCompactInstanceSummary shaded_box_before;
    BObolCompactInstanceSummary shaded_ball_before;
    BObolCompactInstanceSummary shaded_extra_before;
    if (!compact_summary_for_path(pair, "box.s", box_handle,
	    shaded_box_before) ||
	!compact_summary_for_path(pair, "ball.s", ball_handle,
	    shaded_ball_before) ||
	!compact_summary_for_path(pair, "extra.s", shaded_extra_handle,
	    shaded_extra_before))
	FAIL("shaded transform diff should capture retained mesh parts");
    RT_DB_INTERNAL_INIT(&pair_internal);
    if (rt_db_get_internal(&pair_internal, pair_dp, gedp->dbip, NULL) < 0 ||
	pair_internal.idb_type != ID_COMBINATION)
	FAIL("shaded transform diff should read the combination");
    comb = static_cast<struct rt_comb_internal *>(pair_internal.idb_ptr);
    extra_leaf = find_combination_leaf(comb->tree, "extra.s");
    if (!extra_leaf)
	FAIL("shaded transform diff should find the extra occurrence");
    if (!extra_leaf->tr_l.tl_mat)
	extra_leaf->tr_l.tl_mat = static_cast<matp_t>(
	    bu_calloc(1, sizeof(mat_t), "shaded occurrence transform"));
    MAT_IDN(extra_leaf->tr_l.tl_mat);
    MAT_DELTAS(extra_leaf->tr_l.tl_mat, 0.0, 0.0, 2.0);
    if (rt_db_put_internal(pair_dp, gedp->dbip, &pair_internal) < 0)
	FAIL("shaded transform diff should update the combination occurrence");
    struct ged_draw_transaction transform_shaded_pair =
	ged_draw_transaction_make(GED_DRAW_TXN_STALE_SOURCE, "pair.c");
    transform_shaded_pair.view = view.viewContext();
    transform_shaded_pair.stale_reason = GED_DRAW_STALE_SOURCE_CHANGED;
    bobol_performance_counters_set_enabled(1);
    bobol_performance_counters_reset();
    if (!apply_and_sync(gedp, &view, &transform_shaded_pair))
	FAIL("shaded transform should reconcile retained mesh occurrences");
    struct BObolPerformanceCounters shaded_structural_counters;
    bobol_performance_counters_get(&shaded_structural_counters);
    bobol_performance_counters_set_enabled(0);
    if (shaded_structural_counters.mesh_cache_hits < 3 ||
	shaded_structural_counters.mesh_cache_misses != 0 ||
	shaded_structural_counters.plot_calls != 0)
	FAIL("shaded transform should reuse every retained mesh part");
    pair = find_source(controller, "pair.c");
    BObolCompactInstanceSummary shaded_box_after;
    BObolCompactInstanceSummary shaded_ball_after;
    BObolCompactInstanceSummary shaded_extra_after;
    if (!pair || pair->getCompactInstanceCount() != 3 ||
	!pair->getCompactInstanceSummary(box_handle, shaded_box_after) ||
	!pair->getCompactInstanceSummary(ball_handle, shaded_ball_after) ||
	!pair->getCompactInstanceSummary(shaded_extra_handle,
	    shaded_extra_after) ||
	shaded_box_after.geometryIdentity != shaded_box_before.geometryIdentity ||
	shaded_ball_after.geometryIdentity != shaded_ball_before.geometryIdentity ||
	shaded_extra_after.geometryIdentity !=
	    shaded_extra_before.geometryIdentity ||
	!shaded_box_after.localToSource.equals(
	    shaded_box_before.localToSource, 0.000001f) ||
	!shaded_ball_after.localToSource.equals(
	    shaded_ball_before.localToSource, 0.000001f) ||
	shaded_extra_after.localToSource.equals(
	    shaded_extra_before.localToSource, 0.000001f) ||
	shaded_ball_after.selected != shaded_ball_before.selected)
	FAIL("shaded transform diff should preserve mesh parts, handles, siblings, and display state");

    if (!compact_summary_for_path(pair, "box.s", box_handle, box_before) ||
	!compact_summary_for_path(pair, "ball.s", ball_handle, ball_before))
	FAIL("shaded refresh test should resolve stable occurrence handles");
    struct directory *ball_dp = db_lookup(gedp->dbip, "ball.s", LOOKUP_QUIET);
    struct rt_db_internal ball_internal;
    RT_DB_INTERNAL_INIT(&ball_internal);
    if (!ball_dp || rt_db_get_internal(&ball_internal, ball_dp, gedp->dbip,
	    NULL) < 0 || ball_internal.idb_type != ID_ELL)
	FAIL("shaded refresh test should read the sphere primitive");
    struct rt_ell_internal *ell =
	static_cast<struct rt_ell_internal *>(ball_internal.idb_ptr);
    VSCALE(ell->a, ell->a, 1.5);
    VSCALE(ell->b, ell->b, 1.5);
    VSCALE(ell->c, ell->c, 1.5);
    if (rt_db_put_internal(ball_dp, gedp->dbip, &ball_internal) < 0)
	FAIL("shaded refresh test should update the sphere primitive");
    struct ged_draw_transaction stale_ball =
	ged_draw_transaction_make(GED_DRAW_TXN_STALE_SOURCE, "ball.s");
    stale_ball.view = view.viewContext();
    stale_ball.stale_reason = GED_DRAW_STALE_SOURCE_CHANGED;
    if (!apply_and_sync(gedp, &view, &stale_ball))
	FAIL("shaded primitive edit should refresh the matching compact part");
    if (!pair->isCompactInstanceHandleValid(box_handle) ||
	!pair->isCompactInstanceHandleValid(ball_handle) ||
	!pair->getCompactInstanceSummary(box_handle, box_after) ||
	!pair->getCompactInstanceSummary(ball_handle, ball_after) ||
	pair->getCompactInstanceCount() != 3 || pair->stale.getValue() ||
	box_after.geometryIdentity != box_before.geometryIdentity ||
	ball_after.geometryIdentity == ball_before.geometryIdentity)
	FAIL("shaded primitive edit should replace only its mesh part");

    (void)ged_selection_clear(gedp, nullptr);
    ged_selection_recompute(gedp, nullptr);
    controller->clearDatabaseSources();
    struct ged_draw_appearance_settings pair_hidden_appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
    pair_hidden_appearance.draw_mode = GED_DRAW_MODE_HIDDEN_LINE;
    struct ged_draw_transaction draw_pair_hidden =
	ged_draw_transaction_make(GED_DRAW_TXN_DRAW, "pair.c");
    draw_pair_hidden.view = view.viewContext();
    draw_pair_hidden.appearance = &pair_hidden_appearance;
    if (!apply_and_sync(gedp, &view, &draw_pair_hidden))
	FAIL("GED hidden-line draw should publish a compact combination root");
    pair = find_source(controller, "pair.c");
    if (!pair || !pair->isCompactOccurrenceRegistry() ||
	pair->representationMode.getValue() != GED_DRAW_MODE_HIDDEN_LINE ||
	pair->getCompactInstanceCountForPath("pair.c", TRUE) != 3 ||
	!pair->prepareCompiledAssembly() ||
	pair->getCompiledAssemblyInstanceCount() != 3)
	FAIL("hidden-line pair should compile three addressable mesh occurrences");
    if (!ged_selection_select_path(gedp, nullptr, "pair.c/box.s", 1) ||
	!qg_obol_sync_selection_state(gedp, &view, nullptr))
	FAIL("nested hidden-line selection should update the aggregate occurrence");
    selected_count = compact_selected_count(pair);
    if (selected_count != 1)
	FAIL("nested hidden-line selection should select exactly one mesh occurrence");

    BObolCompactInstanceHandle hidden_box_handle;
    BObolCompactInstanceHandle hidden_ball_handle;
    BObolCompactInstanceHandle hidden_extra_handle;
    BObolCompactInstanceSummary hidden_box_before;
    BObolCompactInstanceSummary hidden_ball_before;
    BObolCompactInstanceSummary hidden_extra_before;
    if (!compact_summary_for_path(pair, "box.s", hidden_box_handle,
	    hidden_box_before) ||
	!compact_summary_for_path(pair, "ball.s", hidden_ball_handle,
	    hidden_ball_before) ||
	!compact_summary_for_path(pair, "extra.s", hidden_extra_handle,
	    hidden_extra_before))
	FAIL("hidden-line removal diff should capture retained parts");
    RT_DB_INTERNAL_INIT(&pair_internal);
    if (rt_db_get_internal(&pair_internal, pair_dp, gedp->dbip, NULL) < 0 ||
	pair_internal.idb_type != ID_COMBINATION)
	FAIL("hidden-line removal diff should read the combination");
    comb = static_cast<struct rt_comb_internal *>(pair_internal.idb_ptr);
    if (db_tree_rm_dbleaf(&comb->tree, "extra.s", 0) != 0)
	FAIL("hidden-line removal diff should remove the extra occurrence");
    if (rt_db_put_internal(pair_dp, gedp->dbip, &pair_internal) < 0)
	FAIL("hidden-line removal diff should update the combination");
    struct ged_draw_transaction remove_hidden_member =
	ged_draw_transaction_make(GED_DRAW_TXN_STALE_SOURCE, "pair.c");
    remove_hidden_member.view = view.viewContext();
    remove_hidden_member.stale_reason = GED_DRAW_STALE_SOURCE_CHANGED;
    bobol_performance_counters_set_enabled(1);
    bobol_performance_counters_reset();
    if (!apply_and_sync(gedp, &view, &remove_hidden_member))
	FAIL("hidden-line member removal should reconcile retained occurrences");
    struct BObolPerformanceCounters hidden_remove_counters;
    bobol_performance_counters_get(&hidden_remove_counters);
    bobol_performance_counters_set_enabled(0);
    pair = find_source(controller, "pair.c");
    BObolCompactInstanceSummary hidden_box_after;
    BObolCompactInstanceSummary hidden_ball_after;
    if (!pair || pair->getCompactInstanceCount() != 2 ||
	!pair->getCompactInstanceSummary(hidden_box_handle, hidden_box_after) ||
	!pair->getCompactInstanceSummary(hidden_ball_handle, hidden_ball_after) ||
	pair->isCompactInstanceHandleValid(hidden_extra_handle) ||
	hidden_box_after.geometryIdentity != hidden_box_before.geometryIdentity ||
	hidden_ball_after.geometryIdentity != hidden_ball_before.geometryIdentity ||
	hidden_box_after.selected != hidden_box_before.selected ||
	hidden_remove_counters.mesh_cache_hits < 2 ||
	hidden_remove_counters.mesh_cache_misses != 0)
	FAIL("hidden-line removal should preserve surviving parts, handles, and display state without rebuilding meshes");

    if (db_rename(gedp->dbip, ball_dp, "orb.s") != 0)
	FAIL("hidden-line rename diff should rename the database primitive");
    RT_DB_INTERNAL_INIT(&pair_internal);
    if (rt_db_get_internal(&pair_internal, pair_dp, gedp->dbip, NULL) < 0 ||
	pair_internal.idb_type != ID_COMBINATION)
	FAIL("hidden-line rename diff should read the combination");
    comb = static_cast<struct rt_comb_internal *>(pair_internal.idb_ptr);
    union tree *orb_leaf = find_combination_leaf(comb->tree, "ball.s");
    if (!orb_leaf)
	FAIL("hidden-line rename diff should find the renamed occurrence");
    bu_free(orb_leaf->tr_l.tl_name, "old combination leaf name");
    orb_leaf->tr_l.tl_name = bu_strdup("orb.s");
    if (rt_db_put_internal(pair_dp, gedp->dbip, &pair_internal) < 0)
	FAIL("hidden-line rename diff should update the combination reference");
    struct ged_draw_transaction rename_hidden_member =
	ged_draw_transaction_make(GED_DRAW_TXN_SOURCE_RENAMED, "ball.s");
    rename_hidden_member.new_path = "orb.s";
    rename_hidden_member.view = view.viewContext();
    bobol_performance_counters_set_enabled(1);
    bobol_performance_counters_reset();
    if (!apply_and_sync(gedp, &view, &rename_hidden_member))
	FAIL("hidden-line member rename should reconcile retained occurrences");
    struct BObolPerformanceCounters hidden_rename_counters;
    bobol_performance_counters_get(&hidden_rename_counters);
    bobol_performance_counters_set_enabled(0);
    pair = find_source(controller, "pair.c");
    BObolCompactInstanceSummary hidden_orb_after;
    if (!pair || pair->getCompactInstanceCount() != 2 ||
	!pair->getCompactInstanceSummary(hidden_box_handle, hidden_box_after) ||
	!pair->getCompactInstanceSummary(hidden_ball_handle, hidden_orb_after) ||
	!BU_STR_EQUAL(hidden_orb_after.sourceName.getString(), "orb.s") ||
	hidden_orb_after.geometryIdentity != hidden_ball_before.geometryIdentity ||
	hidden_orb_after.selected != hidden_ball_before.selected ||
	hidden_rename_counters.mesh_cache_hits < 2 ||
	hidden_rename_counters.mesh_cache_misses > 1)
	FAIL("hidden-line rename should preserve aggregate parts, handles, and display state without rebuilding aggregate meshes");

    BObolCompactInstanceSummary material_box_before;
    BObolCompactInstanceSummary material_orb_before;
    if (!pair->getCompactInstanceSummary(hidden_box_handle,
	    material_box_before) ||
	!pair->getCompactInstanceSummary(hidden_ball_handle,
	    material_orb_before))
	FAIL("aggregate material update should capture retained occurrence state");
    RT_DB_INTERNAL_INIT(&pair_internal);
    if (rt_db_get_internal(&pair_internal, pair_dp, gedp->dbip, NULL) < 0 ||
	pair_internal.idb_type != ID_COMBINATION)
	FAIL("aggregate material update should read the combination");
    comb = static_cast<struct rt_comb_internal *>(pair_internal.idb_ptr);
    comb->rgb_valid = 1;
    comb->rgb[0] = 12;
    comb->rgb[1] = 200;
    comb->rgb[2] = 34;
    if (rt_db_put_internal(pair_dp, gedp->dbip, &pair_internal) < 0)
	FAIL("aggregate material update should write the combination color");
    struct ged_draw_transaction material_pair =
	ged_draw_transaction_make(GED_DRAW_TXN_MATERIAL_CHANGED, "pair.c");
    material_pair.view = view.viewContext();
    bobol_performance_counters_set_enabled(1);
    bobol_performance_counters_reset();
    if (!apply_and_sync(gedp, &view, &material_pair))
	FAIL("aggregate material transaction should advance retained material state");
    struct ged_draw_transaction refresh_pair_material =
	ged_draw_transaction_make(GED_DRAW_TXN_REFRESH_MATERIAL_COLORS,
	    "pair.c");
    refresh_pair_material.view = view.viewContext();
    if (!apply_and_sync(gedp, &view, &refresh_pair_material))
	FAIL("aggregate material refresh should update retained style");
    struct BObolPerformanceCounters material_counters;
    bobol_performance_counters_get(&material_counters);
    bobol_performance_counters_set_enabled(0);
    pair = find_source(controller, "pair.c");
    BObolCompactInstanceSummary material_box_after;
    BObolCompactInstanceSummary material_orb_after;
    BObolDatabaseSourceSummary material_source_after;
    if (!pair ||
	!pair->getCompactInstanceSummary(hidden_box_handle,
	    material_box_after) ||
	!pair->getCompactInstanceSummary(hidden_ball_handle,
	    material_orb_after) ||
	!pair->getSummary(material_source_after) ||
	!material_source_after.materialColorValid ||
	material_source_after.materialColor.getValue()[0] < 0.04f ||
	material_source_after.materialColor.getValue()[0] > 0.06f ||
	material_source_after.materialColor.getValue()[1] < 0.77f ||
	material_source_after.materialColor.getValue()[1] > 0.80f ||
	material_source_after.materialColor.getValue()[2] < 0.12f ||
	material_source_after.materialColor.getValue()[2] > 0.15f ||
	material_box_after.geometryIdentity != hidden_box_before.geometryIdentity ||
	material_orb_after.geometryIdentity != hidden_ball_before.geometryIdentity ||
	material_box_after.appearanceRevision <=
	    material_box_before.appearanceRevision ||
	material_orb_after.appearanceRevision <=
	    material_orb_before.appearanceRevision ||
	material_box_after.geometryRevision !=
	    material_box_before.geometryRevision ||
	material_orb_after.geometryRevision !=
	    material_orb_before.geometryRevision ||
	material_box_after.placementRevision !=
	    material_box_before.placementRevision ||
	material_orb_after.placementRevision !=
	    material_orb_before.placementRevision ||
	material_box_after.visibilityRevision !=
	    material_box_before.visibilityRevision ||
	material_orb_after.visibilityRevision !=
	    material_orb_before.visibilityRevision ||
	material_box_after.selectionRevision !=
	    material_box_before.selectionRevision ||
	material_orb_after.selectionRevision !=
	    material_orb_before.selectionRevision ||
	material_counters.mesh_realize_calls != 0)
	FAIL("aggregate material update should change only retained style state");

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
