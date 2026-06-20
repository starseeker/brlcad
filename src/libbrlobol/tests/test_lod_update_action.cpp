/*          T E S T _ L O D _ U P D A T E _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "wdb.h"

#include <Inventor/SbBox.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoSeparator.h>

#include <chrono>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <thread>
#include <vector>

static BRLObolLodRequest
make_request(const char *path, const char *name)
{
    BRLObolLodRequest request;

    request.databaseId = "db://lod-update-action-test";
    request.sourceRevision = 3;
    request.sourceContentHash = 0x123456;
    request.objectPath = path;
    request.objectName = name;
    request.viewRevision = 5;
    request.policyRevision = 7;
    request.drawMode = BRLOBOL_LOD_DRAW_SHADED;
    request.providerId = "lod-update-test";
    request.providerVersion = "1";
    request.qualityTier = BRLOBOL_LOD_QUALITY_ATTRIBUTES;
    request.bounds = SbBox3f(SbVec3f(0.0f, 0.0f, 0.0f),
	    SbVec3f(1.0f, 1.0f, 1.0f));

    return request;
}

static SoBRLMeshShape *
make_mesh(const char *path, const char *name)
{
    SoBRLMeshShape *mesh = new SoBRLMeshShape;
    SbVec3f points[3] = {
	SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 0.0f, 0.0f),
	SbVec3f(0.0f, 1.0f, 0.0f)
    };
    int32_t indices[3] = {0, 1, 2};

    mesh->sourcePath = path;
    mesh->sourceName = name;
    mesh->setIndexedTriangles(points, 3, indices, 3);
    return mesh;
}

static SoBRLLodMeshShape *
make_lod_mesh(const char *path, const char *name)
{
    SoBRLLodMeshShape *mesh = new SoBRLLodMeshShape;
    SbVec3f points[3] = {
	SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 0.0f, 0.0f),
	SbVec3f(0.0f, 1.0f, 0.0f)
    };
    int32_t indices[3] = {0, 1, 2};

    mesh->sourcePath = path;
    mesh->sourceName = name;
    mesh->setIndexedTriangles(points, 3, indices, 3);
    return mesh;
}

static BRLObolLodResult
attributes_result(const BRLObolLodRequest &request, const char *value)
{
    std::vector<BRLObolLodAttribute> attributes;
    BRLObolLodAttribute attribute;

    attribute.name = "display.intent";
    attribute.value = value;
    attributes.push_back(attribute);
    return brlobol_lod_attributes_result(request, attributes);
}

static BRLObolLodResult
mesh_payload_result(const BRLObolLodRequest &request)
{
    BRLObolLodResult result;

    result.request = request;
    result.cacheKey = brlobol_lod_cache_key(request);
    result.resultKind = BRLOBOL_LOD_RESULT_MESH;
    result.qualityTier = BRLOBOL_LOD_QUALITY_FAST_DISPLAY;
    result.providerStatus = BRLOBOL_LOD_PROVIDER_READY;
    result.geometry.kind = BRLOBOL_LOD_GEOMETRY_OBOL_MESH;
    result.geometry.providerId = request.providerId;
    result.geometry.providerVersion = request.providerVersion;
    result.geometry.cacheKey = result.cacheKey;
    result.geometry.activeLevel = 2;
    result.counts.faceCount = 2;
    result.counts.pointCount = 4;
    result.bounds = request.bounds;
    result.mesh.points.push_back(SbVec3f(0.0f, 0.0f, 0.0f));
    result.mesh.points.push_back(SbVec3f(2.0f, 0.0f, 0.0f));
    result.mesh.points.push_back(SbVec3f(0.0f, 2.0f, 0.0f));
    result.mesh.points.push_back(SbVec3f(0.0f, 0.0f, 2.0f));
    result.mesh.coordIndex.push_back(0);
    result.mesh.coordIndex.push_back(1);
    result.mesh.coordIndex.push_back(2);
    result.mesh.coordIndex.push_back(0);
    result.mesh.coordIndex.push_back(3);
    result.mesh.coordIndex.push_back(1);

    return result;
}

static BRLObolLodResult
aabb_task(const BRLObolLodRequest &request, void *UNUSED(userData))
{
    BRLObolLodCounts counts;

    counts.faceCount = 1;
    counts.pointCount = 3;
    return brlobol_lod_aabb_result(request, request.bounds, &counts);
}

static int
wait_for_service(BRLObolLodService &service)
{
    for (int i = 0; i < 400; i++) {
	if (service.inFlightCount() == 0 &&
		service.queuedResultCountForDiagnostics() == 1)
	    return 0;
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    printf("FAIL: LoD service did not produce update-action result\n");
    return 1;
}

static int
make_submit_test_db(char *dbpath, size_t dbpath_len, struct db_i **dbip_out)
{
    static const char *objname = "lod-submit.bot";
    fastf_t vertices[12] = {
	0.0, 0.0, 0.0,
	1.0, 0.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 0.0, 1.0
    };
    int faces[12] = {
	0, 1, 2,
	0, 3, 1,
	1, 3, 2,
	2, 3, 0
    };

    if (!dbpath || dbpath_len == 0 || !dbip_out)
	return 1;
    *dbip_out = NULL;

    FILE *fp = bu_temp_file(dbpath, dbpath_len);
    if (!fp) {
	printf("FAIL: LoD submit temp file\n");
	return 1;
    }
    fclose(fp);

    struct db_i *dbip = db_create(dbpath, 5);
    if (!dbip) {
	printf("FAIL: LoD submit db_create\n");
	bu_file_delete(dbpath);
	return 1;
    }

    struct rt_wdb *wdbp = wdb_dbopen(dbip, RT_WDB_TYPE_DB_DISK);
    if (!wdbp) {
	printf("FAIL: LoD submit wdb_dbopen\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    if (mk_bot(wdbp, objname, RT_BOT_SOLID, RT_BOT_UNORIENTED, 0,
	    4, 4, vertices, faces, NULL, NULL) != 0) {
	printf("FAIL: LoD submit mk_bot\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    *dbip_out = dbip;
    return 0;
}

static int
test_update_action_direct(void)
{
    SoSeparator *root = new SoSeparator;
    root->ref();

    SoBRLMeshShape *meshA = make_mesh("/mesh/a", "a");
    SoBRLMeshShape *meshB = make_mesh("/mesh/b", "b");
    SoBRLMeshShape *meshC = make_mesh("/mesh/c", "c");
    SoBRLMeshShape *meshD = make_mesh("/mesh/d", "d");
    root->addChild(meshA);
    root->addChild(meshB);
    root->addChild(meshC);
    root->addChild(meshD);

    BRLObolLodResult resultA =
	attributes_result(make_request("/mesh/a", "a"), "proxy");
    BRLObolLodResult resultB =
	attributes_result(make_request("/mesh/b", "b"), "stale");
    resultB.providerStatus = BRLOBOL_LOD_PROVIDER_STALE;
    resultB.stale = TRUE;
    resultB.diagnostic = "unit-test stale result";
    BRLObolLodResult resultC =
	attributes_result(make_request("/mesh/c", "c"), "cache-miss");
    resultC.providerStatus = BRLOBOL_LOD_PROVIDER_CACHE_MISS;
    resultC.diagnostic = "unit-test cache miss";
    BRLObolLodResult resultD =
	mesh_payload_result(make_request("/mesh/d", "d"));
    BRLObolLodResult resultMissing =
	attributes_result(make_request("/mesh/missing", "missing"), "missing");

    SoBRLLodUpdateAction update;
    update.addResult(resultA);
    update.addResult(resultB);
    update.addResult(resultC);
    update.addResult(resultD);
    update.addResult(resultMissing);
    update.apply(root);

    if (update.getResultCount() != 5 ||
	    update.getMatchedResultCount() != 4 ||
	    update.getAppliedResultCount() != 2 ||
	    update.getRejectedResultCount() != 2 ||
	    update.getUnmatchedResultCount() != 1 ||
	    update.getDiagnostics().getLength() == 0) {
	printf("FAIL: LoD update action direct counts\n");
	root->unref();
	return 1;
    }

    if (!meshA->lodStagedAvailable.getValue() ||
	    meshA->lodResultKind.getValue() != BRLOBOL_LOD_RESULT_ATTRIBUTES ||
	    meshA->lodAttributeValue.getNum() != 1 ||
	    strcmp(meshA->lodAttributeValue[0].getString(), "proxy") != 0) {
	printf("FAIL: LoD update action did not apply matched result\n");
	root->unref();
	return 1;
    }

    if (meshB->lodStagedAvailable.getValue() ||
	    meshB->lodProviderStatus.getValue() != BRLOBOL_LOD_PROVIDER_STALE ||
	    meshB->lodDiagnostic.getValue().getLength() == 0) {
	printf("FAIL: LoD update action did not reject stale result\n");
	root->unref();
	return 1;
    }

    if (meshC->lodStagedAvailable.getValue() ||
	    meshC->lodProviderStatus.getValue() != BRLOBOL_LOD_PROVIDER_CACHE_MISS ||
	    meshC->lodDiagnostic.getValue().getLength() == 0) {
	printf("FAIL: LoD update action did not reject cache-miss result\n");
	root->unref();
	return 1;
    }

    if (!meshD->lodStagedAvailable.getValue() ||
	    !meshD->lodAvailable.getValue() ||
	    meshD->lodResultKind.getValue() != BRLOBOL_LOD_RESULT_MESH ||
	    meshD->lodActiveLevel.getValue() != 2 ||
	    meshD->lodFaceCount.getValue() != 2 ||
	    meshD->lodPointCount.getValue() != 4 ||
	    meshD->point.getNum() != 4 ||
	    meshD->coordIndex.getNum() != 6 ||
	    meshD->getTriangleCount() != 2 ||
	    !meshD->hasFullDetailMesh() ||
	    meshD->getFullDetailTriangleCount() != 1) {
	printf("FAIL: LoD update action did not apply staged mesh payload\n");
	root->unref();
	return 1;
    }

    BRLObolLodRequest meshDRequest;
    meshD->makeLodRequest(meshDRequest,
	    "db://lod-update-action-test",
	    10,
	    11,
	    12,
	    BRLOBOL_LOD_DRAW_SHADED,
	    "rt_mesh_lod",
	    "rt-cache-v1",
	    BRLOBOL_LOD_QUALITY_FAST_DISPLAY);
    if (meshDRequest.sourceCounts.faceCount != 1 ||
	    meshDRequest.sourceCounts.pointCount != 3 ||
	    meshDRequest.bounds.isEmpty() ||
	    meshDRequest.bounds.getMax()[0] > 1.0f ||
	    meshDRequest.bounds.getMax()[1] > 1.0f) {
	printf("FAIL: LoD-backed mesh request did not preserve full-detail source identity\n");
	root->unref();
	return 1;
    }

    SoBRLExportAction exactExport;
    exactExport.apply(root);
    if (exactExport.getGeometryPolicy() != SoBRLExportAction::FULL_DETAIL ||
	    exactExport.getTriangleCount() != 4 ||
	    exactExport.getSkippedLodDisplayMeshCount() != 1) {
	printf("FAIL: default export did not use preserved full-detail mesh\n");
	root->unref();
	return 1;
    }

    SoBRLExportAction displayExport;
    displayExport.setGeometryPolicy(SoBRLExportAction::DISPLAY_LEVEL);
    displayExport.apply(root);
    if (displayExport.getGeometryPolicy() != SoBRLExportAction::DISPLAY_LEVEL ||
	    displayExport.getTriangleCount() != 5 ||
	    displayExport.getSkippedLodDisplayMeshCount() != 0) {
	printf("FAIL: display-level export did not include active display LoD mesh\n");
	root->unref();
	return 1;
    }

    SoBRLMeasureAction exactMeasure;
    exactMeasure.apply(root);
    if (exactMeasure.getGeometryPolicy() != SoBRLMeasureAction::FULL_DETAIL ||
	    exactMeasure.getTriangleCount() != 4 ||
	    exactMeasure.getSkippedLodDisplayMeshCount() != 1 ||
	    fabsf(exactMeasure.getSurfaceArea() - 2.0f) > 1.0e-5f) {
	printf("FAIL: default measure did not use preserved full-detail mesh\n");
	root->unref();
	return 1;
    }

    SoBRLMeasureAction displayMeasure;
    displayMeasure.setGeometryPolicy(SoBRLMeasureAction::DISPLAY_LEVEL);
    displayMeasure.apply(root);
    if (displayMeasure.getGeometryPolicy() != SoBRLMeasureAction::DISPLAY_LEVEL ||
	    displayMeasure.getTriangleCount() != 5 ||
	    displayMeasure.getSkippedLodDisplayMeshCount() != 0 ||
	    fabsf(displayMeasure.getSurfaceArea() - 5.5f) > 1.0e-5f) {
	printf("FAIL: display-level measure did not include active display LoD mesh\n");
	root->unref();
	return 1;
    }

    SoBRLSnapAction displaySnap;
    displaySnap.setEnabledKinds(SoBRLSnapAction::FACE_NEAREST);
    displaySnap.setQueryPoint(SbVec3f(1.5f, 0.2f, 0.0f));
    displaySnap.setTolerance(0.1f);
    displaySnap.apply(root);
    if (displaySnap.getGeometryPolicy() != SoBRLSnapAction::DISPLAY_LEVEL ||
	    !displaySnap.hasCandidate() ||
	    strcmp(displaySnap.getPath().getString(), "/mesh/d") != 0 ||
	    displaySnap.getSkippedLodDisplayMeshCount() != 0) {
	printf("FAIL: display-level snap did not use active display LoD mesh\n");
	root->unref();
	return 1;
    }

    SoBRLSnapAction exactSnap;
    exactSnap.setGeometryPolicy(SoBRLSnapAction::FULL_DETAIL);
    exactSnap.setEnabledKinds(SoBRLSnapAction::FACE_NEAREST);
    exactSnap.setQueryPoint(SbVec3f(1.5f, 0.2f, 0.0f));
    exactSnap.setTolerance(0.1f);
    exactSnap.apply(root);
    if (exactSnap.getGeometryPolicy() != SoBRLSnapAction::FULL_DETAIL ||
	    exactSnap.hasCandidate() ||
	    exactSnap.getSkippedLodDisplayMeshCount() != 1) {
	printf("FAIL: full-detail snap did not skip active display LoD mesh\n");
	root->unref();
	return 1;
    }

    meshA->selectable = FALSE;
    meshB->selectable = FALSE;
    meshC->selectable = FALSE;
    SoBRLSnapAction exactFullSnap;
    exactFullSnap.setGeometryPolicy(SoBRLSnapAction::FULL_DETAIL);
    exactFullSnap.setEnabledKinds(SoBRLSnapAction::FACE_NEAREST);
    exactFullSnap.setQueryPoint(SbVec3f(0.2f, 0.2f, 0.0f));
    exactFullSnap.setTolerance(0.1f);
    exactFullSnap.apply(root);
    if (!exactFullSnap.hasCandidate() ||
	    strcmp(exactFullSnap.getPath().getString(), "/mesh/d") != 0 ||
	    exactFullSnap.getSkippedLodDisplayMeshCount() != 1) {
	printf("FAIL: full-detail snap did not use preserved full-detail mesh\n");
	root->unref();
	return 1;
    }

    if (meshD->getPickGeometryPolicy() != SoBRLMeshShape::PICK_DISPLAY_LEVEL) {
	printf("FAIL: mesh pick policy should default to display-level geometry\n");
	root->unref();
	return 1;
    }

    SbViewportRegion pickViewport(100, 100);
    SoRayPickAction displayPick(pickViewport);
    displayPick.setRay(SbVec3f(1.5f, 0.2f, 5.0f),
	    SbVec3f(0.0f, 0.0f, -1.0f));
    displayPick.apply(root);
    const SoPickedPoint *pickedPoint = displayPick.getPickedPoint();
    const SoDetail *rawDetail = pickedPoint ? pickedPoint->getDetail(meshD) :
	NULL;
    const SoBRLPickDetail *pickDetail =
	(rawDetail && rawDetail->isOfType(SoBRLPickDetail::getClassTypeId())) ?
	static_cast<const SoBRLPickDetail *>(rawDetail) : NULL;
    if (!pickDetail ||
	    strcmp(pickDetail->getPath().getString(), "/mesh/d") != 0 ||
	    pickDetail->getPrimitiveKind() != SoBRLPickDetail::FACE) {
	printf("FAIL: display-level pick did not use active display LoD mesh\n");
	root->unref();
	return 1;
    }

    meshD->setPickGeometryPolicy(SoBRLMeshShape::PICK_FULL_DETAIL);
    SoRayPickAction exactMissPick(pickViewport);
    exactMissPick.setRay(SbVec3f(1.5f, 0.2f, 5.0f),
	    SbVec3f(0.0f, 0.0f, -1.0f));
    exactMissPick.apply(root);
    if (exactMissPick.getPickedPoint()) {
	printf("FAIL: full-detail pick did not skip active display LoD mesh\n");
	root->unref();
	return 1;
    }

    SoRayPickAction exactHitPick(pickViewport);
    exactHitPick.setRay(SbVec3f(0.2f, 0.2f, 5.0f),
	    SbVec3f(0.0f, 0.0f, -1.0f));
    exactHitPick.apply(root);
    pickedPoint = exactHitPick.getPickedPoint();
    rawDetail = pickedPoint ? pickedPoint->getDetail(meshD) : NULL;
    pickDetail =
	(rawDetail && rawDetail->isOfType(SoBRLPickDetail::getClassTypeId())) ?
	static_cast<const SoBRLPickDetail *>(rawDetail) : NULL;
    if (!pickDetail ||
	    strcmp(pickDetail->getPath().getString(), "/mesh/d") != 0 ||
	    pickDetail->getPrimitiveIndex() != 0 ||
	    pickDetail->getFaceVertexIndexA() != 0 ||
	    pickDetail->getFaceVertexIndexB() != 1 ||
	    pickDetail->getFaceVertexIndexC() != 2) {
	printf("FAIL: full-detail pick did not use preserved full-detail mesh\n");
	root->unref();
	return 1;
    }
    meshD->setPickGeometryPolicy(SoBRLMeshShape::PICK_DISPLAY_LEVEL);

    root->unref();
    return 0;
}

static int
test_update_action_service_drain(void)
{
    SoSeparator *root = new SoSeparator;
    root->ref();

    SoBRLMeshShape *mesh = make_mesh("/drained.bot", "drained.bot");
    root->addChild(mesh);

    BRLObolLodService service;
    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD update action service did not start\n");
	root->unref();
	return 1;
    }

    BRLObolLodTask task;
    task.generation = service.beginGeneration();
    task.request = make_request("drained.bot", "drained.bot");
    task.realize = aabb_task;
    if (service.submit(task) == 0) {
	printf("FAIL: LoD update action service did not accept task\n");
	service.stop();
	root->unref();
	return 1;
    }

    if (wait_for_service(service)) {
	service.stop();
	root->unref();
	return 1;
    }

    SoBRLLodUpdateAction update;
    if (update.drainService(service) != 1 ||
	    update.getResultCount() != 1) {
	printf("FAIL: LoD update action did not drain service result\n");
	service.stop();
	root->unref();
	return 1;
    }

    update.apply(root);
    if (update.getAppliedResultCount() != 1 ||
	    update.getUnmatchedResultCount() != 0 ||
	    !mesh->lodStagedAvailable.getValue() ||
	    mesh->lodResultKind.getValue() != BRLOBOL_LOD_RESULT_AABB ||
	    mesh->lodFaceCount.getValue() != 1 ||
	    mesh->lodPointCount.getValue() != 3) {
	printf("FAIL: LoD update action did not apply drained service result\n");
	service.stop();
	root->unref();
	return 1;
    }

    service.stop();
    root->unref();
    return 0;
}

static int
test_mesh_lod_request_and_view_info(void)
{
    SoBRLMeshShape *mesh = make_mesh("/mesh/request", "request");
    mesh->ref();
    mesh->sourceId = 44;
    mesh->lodPolicy = 9;
    if (mesh->isLodBackedMesh() || mesh->isLodDisplayActive() ||
	    !mesh->isLodPreserveFullDetailEnabled() ||
	    mesh->getPickGeometryPolicy() != SoBRLMeshShape::PICK_DISPLAY_LEVEL) {
	printf("FAIL: base mesh should not start in LoD-backed display mode\n");
	mesh->unref();
	return 1;
    }

    SoBRLLodMeshShape *lodMesh = make_lod_mesh("/mesh/lod-request",
	    "lod-request");
    lodMesh->ref();
    if (!lodMesh->isLodBackedMesh() || lodMesh->isLodDisplayActive() ||
	    lodMesh->isLodPreserveFullDetailEnabled() ||
	    lodMesh->getPickGeometryPolicy() != SoBRLMeshShape::PICK_DISPLAY_LEVEL) {
	printf("FAIL: LoD mesh shape did not initialize LoD-backed mode\n");
	lodMesh->unref();
	mesh->unref();
	return 1;
    }
    lodMesh->unref();

    BRLObolLodRequest request;
    mesh->makeLodRequest(request,
	    "db://request-test",
	    12,
	    34,
	    56,
	    BRLOBOL_LOD_DRAW_SHADED,
	    "rt_mesh_lod",
	    "rt-cache-v1",
	    BRLOBOL_LOD_QUALITY_FAST_DISPLAY);

    if (strcmp(request.databaseId.getString(), "db://request-test") != 0 ||
	    request.databaseRevision != 12 ||
	    request.sourceRevision != 44 ||
	    strcmp(request.objectPath.getString(), "/mesh/request") != 0 ||
	    strcmp(request.objectName.getString(), "request") != 0 ||
	    request.viewRevision != 34 ||
	    request.policyRevision != 56 ||
	    request.drawMode != BRLOBOL_LOD_DRAW_SHADED ||
	    request.lodPolicy != 9 ||
	    strcmp(request.providerId.getString(), "rt_mesh_lod") != 0 ||
	    request.qualityTier != BRLOBOL_LOD_QUALITY_FAST_DISPLAY ||
	    request.sourceCounts.faceCount != 1 ||
	    request.sourceCounts.pointCount != 3 ||
	    request.bounds.isEmpty()) {
	printf("FAIL: mesh LoD request helper did not preserve identity\n");
	mesh->unref();
	return 1;
    }

    BRLObolLodCacheKey firstKey = brlobol_lod_cache_key(request);
    BRLObolLodRequest sameRequest;
    mesh->makeLodRequest(sameRequest,
	    "db://request-test",
	    12,
	    34,
	    56,
	    BRLOBOL_LOD_DRAW_SHADED,
	    "rt_mesh_lod",
	    "rt-cache-v1",
	    BRLOBOL_LOD_QUALITY_FAST_DISPLAY);
    if (strcmp(firstKey.value.getString(),
	    brlobol_lod_cache_key(sameRequest).value.getString()) != 0) {
	printf("FAIL: mesh LoD request helper produced unstable cache key\n");
	mesh->unref();
	return 1;
    }

    SoSeparator *root = new SoSeparator;
    root->ref();
    SoOrthographicCamera *camera = new SoOrthographicCamera;
    camera->height = 20.0f;
    BRLObolViewController controller(root, camera);
    controller.setViewportSize(320, 240);

    struct rt_view_info info = RT_VIEW_INFO_INIT;
    if (!controller.getRtViewInfo(&info) ||
	    info.width != 320 ||
	    info.height != 240 ||
	    fabs(info.size - 20.0) > 1.0e-6) {
	printf("FAIL: view controller did not produce RT view info\n");
	root->unref();
	mesh->unref();
	return 1;
    }

    controller.setCamera(NULL);
    if (controller.getRtViewInfo(&info) ||
	    info.width != 320 ||
	    info.height != 240 ||
	    info.size <= 0.0) {
	printf("FAIL: view controller missing-camera RT view fallback failed\n");
	root->unref();
	mesh->unref();
	return 1;
    }

    root->unref();
    mesh->unref();
    return 0;
}

static int
test_mesh_lod_submit_action(void)
{
    char cache_dir[MAXPATHLEN] = {0};
    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;

    bu_dir(cache_dir, MAXPATHLEN, BU_DIR_CURR,
	    "brlobol_lod_submit_action_cache", NULL);
    bu_dirclear(cache_dir);
    bu_mkdir(cache_dir);
    bu_setenv("BU_DIR_CACHE", cache_dir, 1);

    if (make_submit_test_db(dbpath, sizeof(dbpath), &dbip)) {
	bu_dirclear(cache_dir);
	return 1;
    }

    SoSeparator *root = new SoSeparator;
    root->ref();
    SoBRLMeshShape *plainMesh = make_mesh("/plain.bot", "plain.bot");
    SoBRLMeshShape *mesh = make_lod_mesh("/lod-submit.bot", "lod-submit.bot");
    plainMesh->sourceId = 100;
    mesh->sourceId = 101;
    mesh->lodPolicy = 4;
    root->addChild(plainMesh);
    root->addChild(mesh);

    BRLObolLodService service;
    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD submit action service did not start\n");
	root->unref();
	db_mesh_lod_clear(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    struct rt_view_info view = RT_VIEW_INFO_INIT;
    view.width = 640;
    view.height = 480;
    view.size = 100.0;

    SoBRLMeshLodSubmitAction submit;
    submit.setService(&service);
    submit.setDatabase(dbip, "db://lod-submit-test", 2026);
    submit.setViewInfo(&view);
    submit.setGeneration(service.beginGeneration());
    submit.setRevisions(55, 66);
    submit.apply(root);

    if (submit.getVisitedMeshCount() != 2 ||
	    submit.getSubmittedTaskCount() != 1 ||
	    submit.getSkippedMeshCount() != 1 ||
	    submit.getDiagnostics().getLength() == 0) {
	printf("FAIL: LoD submit action did not submit expected mesh task\n");
	service.stop();
	root->unref();
	db_mesh_lod_clear(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    if (wait_for_service(service)) {
	service.stop();
	root->unref();
	db_mesh_lod_clear(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    SoBRLLodUpdateAction update;
    if (update.drainService(service) != 1 ||
	    update.getResultCount() != 1) {
	printf("FAIL: LoD submit action result drain failed\n");
	service.stop();
	root->unref();
	db_mesh_lod_clear(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    update.apply(root);
    int ret = 0;
    if (update.getAppliedResultCount() != 1 ||
	    update.getRejectedResultCount() != 0 ||
	    !mesh->lodAvailable.getValue() ||
	    !mesh->isLodDisplayActive() ||
	    mesh->lodResultKind.getValue() != BRLOBOL_LOD_RESULT_MESH ||
	    mesh->lodFaceCount.getValue() != 4 ||
	    mesh->lodPointCount.getValue() != 4 ||
	    mesh->getTriangleCount() != 4) {
	printf("FAIL: LoD submit action result did not update mesh payload\n");
	ret = 1;
    }
    if (!ret && mesh->hasFullDetailMesh()) {
	printf("FAIL: LoD-backed mesh retained full-detail payload after display LoD update\n");
	ret = 1;
    }
    if (!ret) {
	BRLObolLodRequest activeRequest;
	mesh->makeLodRequest(activeRequest,
		"db://lod-submit-test",
		2026,
		55,
		66,
		BRLOBOL_LOD_DRAW_SHADED,
		"rt_mesh_lod",
		"rt-cache-v1",
		BRLOBOL_LOD_QUALITY_FAST_DISPLAY);
	if (activeRequest.sourceCounts.faceCount != 1 ||
		activeRequest.sourceCounts.pointCount != 3 ||
		activeRequest.bounds.isEmpty() ||
		activeRequest.bounds.getMax()[0] > 1.0f ||
		activeRequest.bounds.getMax()[1] > 1.0f) {
	    printf("FAIL: LoD-backed mesh request did not keep source metrics without full-detail payload\n");
	    ret = 1;
	}
    }

    if (!ret) {
	int forcedLevel = mesh->lodActiveLevel.getValue();
	SoBRLMeshLodSubmitAction forcedSubmit;
	forcedSubmit.setService(&service);
	forcedSubmit.setDatabase(dbip, "db://lod-submit-test", 2026);
	forcedSubmit.setViewInfo(&view);
	forcedSubmit.setGeneration(service.beginGeneration());
	forcedSubmit.setRevisions(77, 88);
	forcedSubmit.setForcedLevel(forcedLevel);
	forcedSubmit.apply(root);

	if (!forcedSubmit.hasForcedLevel() ||
		forcedSubmit.getForcedLevel() != forcedLevel ||
		forcedSubmit.getSubmittedTaskCount() != 1) {
	    printf("FAIL: LoD submit action forced-level policy was not applied\n");
	    ret = 1;
	} else if (wait_for_service(service)) {
	    ret = 1;
	} else {
	    std::vector<BRLObolLodResult> forcedResults;
	    service.drainResults(forcedResults);
	    if (forcedResults.size() != 1 ||
		    forcedResults[0].geometry.activeLevel != forcedLevel ||
		    !forcedResults[0].mesh.isValid()) {
		printf("FAIL: LoD submit action forced-level result was not returned\n");
		ret = 1;
	    } else {
		SoBRLLodUpdateAction forcedUpdate;
		forcedUpdate.setResults(forcedResults);
		forcedUpdate.apply(root);
		if (forcedUpdate.getAppliedResultCount() != 1 ||
			mesh->lodActiveLevel.getValue() != forcedLevel) {
		    printf("FAIL: LoD submit action forced-level result was not applied\n");
		    ret = 1;
		}
	    }
	}
    }

    service.stop();
    root->unref();
    db_mesh_lod_clear(dbip);
    db_close(dbip);
    bu_file_delete(dbpath);
    bu_dirclear(cache_dir);
    return ret;
}

static int
test_view_controller_lod_submit_and_apply(void)
{
    char cache_dir[MAXPATHLEN] = {0};
    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;

    bu_dir(cache_dir, MAXPATHLEN, BU_DIR_CURR,
	    "brlobol_lod_view_controller_cache", NULL);
    bu_dirclear(cache_dir);
    bu_mkdir(cache_dir);
    bu_setenv("BU_DIR_CACHE", cache_dir, 1);

    if (make_submit_test_db(dbpath, sizeof(dbpath), &dbip)) {
	bu_dirclear(cache_dir);
	return 1;
    }

    SoSeparator *root = new SoSeparator;
    root->ref();
    SoOrthographicCamera *camera = new SoOrthographicCamera;
    camera->height = 80.0f;

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = "lod-submit.bot";
    source->lodBotThreshold = 1;
    source->sourceRevision = 303;
    source->realizedSourceRevision = 303;
    source->realizationStatus = SoBRLDatabaseSource::REALIZED;
    source->stale = FALSE;
    source->staleReason = SoBRLDatabaseSource::STALE_NONE;

    SoBRLMeshShape *mesh = make_lod_mesh("/lod-submit.bot", "lod-submit.bot");
    mesh->sourceId = 303;
    mesh->lodPolicy = 6;
    source->addChild(mesh);
    root->addChild(source);

    BRLObolLodService service;
    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD view-controller service did not start\n");
	root->unref();
	db_mesh_lod_clear(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    {
	BRLObolViewController controller(root, camera);
	controller.setLodService(&service);
	controller.setViewportSize(800, 600);
	controller.setLodPolicyRevision(99);
	controller.clearRenderRequest();

	uint64_t viewRevision = controller.getLodViewRevision();
	if (viewRevision == 0 || controller.getLodPolicyRevision() != 99) {
	    printf("FAIL: LoD view-controller revisions were not initialized\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	if (controller.submitLodRequests(NULL, service.beginGeneration()) != 1 ||
		controller.getLastLodVisitedMeshCount() != 1 ||
		controller.getLastLodSubmittedTaskCount() != 1 ||
		controller.getLastLodSkippedMeshCount() != 0 ||
		controller.getLastLodDiagnostics().getLength() != 0 ||
		!controller.isRenderRequested() ||
		strcmp(controller.getRenderReason().getString(),
		    "lod-submit") != 0) {
	    printf("FAIL: LoD view controller did not submit expected task\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	if (wait_for_service(service)) {
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	if (!controller.hasPendingLodResults() || !controller.isRenderRequested()) {
	    printf("FAIL: LoD view controller did not observe result-ready wakeup\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	controller.clearRenderRequest();
	if (!controller.hasPendingLodResults() ||
		controller.processPendingLodResults() != 1 ||
		controller.getLastLodResultCount() != 1 ||
		controller.getLastLodMatchedResultCount() != 1 ||
		controller.getLastLodAppliedResultCount() != 1 ||
		controller.getLastLodRejectedResultCount() != 0 ||
		controller.getLastLodUnmatchedResultCount() != 0 ||
		controller.getLastLodDiagnostics().getLength() != 0 ||
		!controller.isRenderRequested() ||
		strcmp(controller.getRenderReason().getString(),
		    "lod-result") != 0 ||
		controller.hasPendingLodResults()) {
	    printf("FAIL: LoD view controller did not apply service result\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	controller.clearRenderRequest();
	controller.setLodAutoSubmit(TRUE);
	if (!controller.isLodAutoSubmitEnabled() ||
		!controller.isRenderRequested() ||
		strcmp(controller.getRenderReason().getString(),
		    "lod-auto-submit") != 0) {
	    printf("FAIL: LoD view controller did not enable auto submit\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	controller.clearRenderRequest();
	if (controller.submitLodRequestsIfNeeded() != 1 ||
		controller.getLastLodSubmittedTaskCount() != 1 ||
		!controller.isRenderRequested() ||
		strcmp(controller.getRenderReason().getString(),
		    "lod-submit") != 0) {
	    printf("FAIL: LoD view controller did not auto-submit changed scene\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (controller.submitLodRequestsIfNeeded() != 0) {
	    printf("FAIL: LoD view controller duplicated unchanged auto-submit\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (wait_for_service(service)) {
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	controller.clearRenderRequest();
	if (controller.processPendingLodResults() != 1 ||
		controller.getLastLodAppliedResultCount() != 1 ||
		strcmp(controller.getRenderReason().getString(),
		    "lod-result") != 0) {
	    printf("FAIL: LoD view controller did not apply auto-submit result\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	source->lodBotThreshold = 2;
	source->realizationStatus = SoBRLDatabaseSource::REALIZED;
	source->stale = FALSE;
	source->staleReason = SoBRLDatabaseSource::STALE_NONE;
	controller.clearRenderRequest();
	if (controller.submitLodRequestsIfNeeded() != 1 ||
		controller.getLastLodSubmittedTaskCount() != 1 ||
		!controller.isRenderRequested() ||
		strcmp(controller.getRenderReason().getString(),
		    "lod-submit") != 0) {
	    printf("FAIL: LoD view controller did not auto-submit threshold policy change\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (wait_for_service(service)) {
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	controller.clearRenderRequest();
	if (controller.processPendingLodResults() != 1 ||
		controller.getLastLodAppliedResultCount() != 1 ||
		strcmp(controller.getRenderReason().getString(),
		    "lod-result") != 0) {
	    printf("FAIL: LoD view controller did not apply threshold policy result\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	uint64_t previousViewRevision = controller.getLodViewRevision();
	camera->height = 90.0f;
	controller.clearRenderRequest();
	if (controller.submitLodRequestsIfNeeded() != 1 ||
		controller.getLastLodSubmittedTaskCount() != 1 ||
		controller.getLodViewRevision() == previousViewRevision ||
		!controller.isRenderRequested() ||
		strcmp(controller.getRenderReason().getString(),
		    "lod-submit") != 0) {
	    printf("FAIL: LoD view controller did not auto-submit camera field change\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (wait_for_service(service)) {
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	controller.clearRenderRequest();
	if (controller.processPendingLodResults() != 1 ||
		controller.getLastLodAppliedResultCount() != 1 ||
		strcmp(controller.getRenderReason().getString(),
		    "lod-result") != 0) {
	    printf("FAIL: LoD view controller did not apply camera field result\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	int forcedLevel = mesh->lodActiveLevel.getValue();
	uint64_t previousPolicyRevision = controller.getLodPolicyRevision();
	controller.clearRenderRequest();
	controller.setLodForcedLevel(forcedLevel);
	if (!controller.hasLodForcedLevel() ||
		controller.getLodForcedLevel() != forcedLevel ||
		controller.getLodPolicyRevision() == previousPolicyRevision ||
		!controller.isRenderRequested() ||
		strcmp(controller.getRenderReason().getString(),
		    "lod-policy") != 0) {
	    printf("FAIL: LoD view controller did not record forced-level policy\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	controller.clearRenderRequest();
	if (controller.submitLodRequestsIfNeeded() != 1 ||
		controller.getLastLodSubmittedTaskCount() != 1 ||
		!controller.isRenderRequested() ||
		strcmp(controller.getRenderReason().getString(),
		    "lod-submit") != 0) {
	    printf("FAIL: LoD view controller did not auto-submit forced-level policy\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (wait_for_service(service)) {
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	controller.clearRenderRequest();
	if (controller.processPendingLodResults() != 1 ||
		controller.getLastLodAppliedResultCount() != 1 ||
		mesh->lodActiveLevel.getValue() != forcedLevel ||
		strcmp(controller.getRenderReason().getString(),
		    "lod-result") != 0) {
	    printf("FAIL: LoD view controller did not apply forced-level result\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	previousPolicyRevision = controller.getLodPolicyRevision();
	controller.clearRenderRequest();
	controller.clearLodForcedLevel();
	if (controller.hasLodForcedLevel() ||
		controller.getLodPolicyRevision() == previousPolicyRevision ||
		!controller.isRenderRequested() ||
		strcmp(controller.getRenderReason().getString(),
		    "lod-policy") != 0) {
	    printf("FAIL: LoD view controller did not clear forced-level policy\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	previousViewRevision = controller.getLodViewRevision();
	controller.setViewportSize(801, 600);
	if (controller.getLodViewRevision() == previousViewRevision) {
	    printf("FAIL: LoD view controller did not advance view revision\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	controller.clearRenderRequest();
	if (controller.submitLodRequestsIfNeeded() != 1 ||
		controller.getLastLodSubmittedTaskCount() != 1) {
	    printf("FAIL: LoD view controller did not auto-submit view change\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
    }

    int ret = 0;
    if (!mesh->lodAvailable.getValue() ||
	    mesh->lodResultKind.getValue() != BRLOBOL_LOD_RESULT_MESH ||
	    mesh->lodFaceCount.getValue() != 4 ||
	    mesh->lodPointCount.getValue() != 4 ||
	    mesh->getTriangleCount() != 4) {
	printf("FAIL: LoD view-controller result did not update mesh payload\n");
	ret = 1;
    }

    service.stop();
    root->unref();
    db_mesh_lod_clear(dbip);
    db_close(dbip);
    bu_file_delete(dbpath);
    bu_dirclear(cache_dir);
    return ret;
}

int
main(int argc, char **argv)
{
    if (argc != 1) {
	printf("Usage: %s\n", argv[0]);
	return 1;
    }
    bu_setprogname(argv[0]);

    brlobol_init(NULL);

    if (test_update_action_direct())
	return 1;
    if (test_update_action_service_drain())
	return 1;
    if (test_mesh_lod_request_and_view_info())
	return 1;
    if (test_mesh_lod_submit_action())
	return 1;
    if (test_view_controller_lod_submit_and_apply())
	return 1;

    return 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
