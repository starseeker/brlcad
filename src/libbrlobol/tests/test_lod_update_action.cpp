/*          T E S T _ L O D _ U P D A T E _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol.h"
#include "bv.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "rt/view.h"
#include "wdb.h"

#include <Inventor/SbBox.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbViewVolume.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoGroup.h>
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

static int
find_database_source_summary_by_instance(SoBRLSceneController &scene,
	const char *instanceKey,
	BRLObolDatabaseSourceSummary &summary)
{
    summary = BRLObolDatabaseSourceSummary();
    if (!instanceKey || !instanceKey[0])
	return 0;

    const int count = scene.getDatabaseSourceCount();
    for (int i = 0; i < count; i++) {
	BRLObolDatabaseSourceSummary candidate;
	if (!scene.getDatabaseSourceSummary(i, candidate) ||
	    !candidate.valid)
	    continue;
	if (strcmp(candidate.instanceKey.getString(), instanceKey) == 0) {
	    summary = candidate;
	    return 1;
	}
    }

    return 0;
}

static int
matrix_nearly_equal(const SbMatrix &a, const SbMatrix &b)
{
    return a.equals(b, 0.000001f);
}

static int
check_source_mesh_request(const BRLObolSourceMeshRequest &request,
			  const char *path,
			  const char *name,
			  uint32_t sourceId)
{
    if (strcmp(request.path.getString(), path) != 0 ||
	strcmp(request.sourceName.getString(), name) != 0 ||
	request.sourceId != sourceId ||
	request.faceCount != 1 ||
	request.pointCount != 3 ||
	request.bounds.isEmpty() ||
	request.bounds.getMax()[0] > 1.0f ||
	request.bounds.getMax()[1] > 1.0f) {
	return 1;
    }

    return 0;
}

static int
check_source_full_detail_lod_request(const BRLObolLodRequest &request,
				     const char *path,
				     const char *name)
{
    if (strcmp(request.objectPath.getString(), path) != 0 ||
	strcmp(request.objectName.getString(), name) != 0 ||
	strcmp(request.providerId.getString(), "rt_source_full_detail") != 0 ||
	strcmp(request.providerVersion.getString(), "direct-bot-v1") != 0 ||
	request.qualityTier != BRLOBOL_LOD_QUALITY_FULL_DETAIL ||
	request.sourceCounts.faceCount != 1 ||
	request.sourceCounts.pointCount != 3 ||
	request.bounds.isEmpty() ||
	request.bounds.getMax()[0] > 1.0f ||
	request.bounds.getMax()[1] > 1.0f) {
	return 1;
    }

    return 0;
}

static int
has_lod_provider_param(const BRLObolLodRequest &request, const char *name)
{
    for (size_t i = 0; i < request.providerParams.size(); i++) {
	if (strcmp(request.providerParams[i].name.getString(), name) == 0)
	    return 1;
    }

    return 0;
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
mesh_payload_task(const BRLObolLodRequest &request, void *UNUSED(userData))
{
    return mesh_payload_result(request);
}

static BRLObolLodResult
mesh_payload_variant_result(const BRLObolLodRequest &request,
			    int activeLevel,
			    int triangleCount)
{
    BRLObolLodResult result = mesh_payload_result(request);

    result.geometry.activeLevel = activeLevel;
    result.mesh.clear();
    result.mesh.points.push_back(SbVec3f(0.0f, 0.0f, 0.0f));
    result.mesh.points.push_back(SbVec3f(2.0f, 0.0f, 0.0f));
    result.mesh.points.push_back(SbVec3f(0.0f, 2.0f, 0.0f));
    result.mesh.points.push_back(SbVec3f(0.0f, 0.0f, 2.0f));
    result.mesh.points.push_back(SbVec3f(2.0f, 2.0f, 0.0f));
    static const int32_t indices[] = {
	0, 1, 2,
	0, 3, 1,
	1, 4, 2
    };
    const int clampedTriangleCount =
	triangleCount < 1 ? 1 : (triangleCount > 3 ? 3 : triangleCount);
    for (int i = 0; i < clampedTriangleCount * 3; i++)
	result.mesh.coordIndex.push_back(indices[i]);
    result.counts.faceCount = static_cast<uint64_t>(clampedTriangleCount);
    result.counts.pointCount = static_cast<uint64_t>(result.mesh.points.size());

    return result;
}

static BRLObolLodResult
mesh_payload_coarse_task(const BRLObolLodRequest &request, void *UNUSED(userData))
{
    return mesh_payload_variant_result(request, 1, 1);
}

static BRLObolLodResult
mesh_payload_refined_task(const BRLObolLodRequest &request, void *UNUSED(userData))
{
    return mesh_payload_variant_result(request, 2, 3);
}

static BRLObolLodResult
source_full_detail_result(const BRLObolLodRequest &request)
{
    BRLObolLodResult result;

    result.request = request;
    result.cacheKey = brlobol_lod_cache_key(request);
    result.resultKind = BRLOBOL_LOD_RESULT_FULL_DETAIL;
    result.qualityTier = BRLOBOL_LOD_QUALITY_FULL_DETAIL;
    result.providerStatus = BRLOBOL_LOD_PROVIDER_READY;
    result.geometry.kind = BRLOBOL_LOD_GEOMETRY_OBOL_MESH;
    result.geometry.providerId = request.providerId;
    result.geometry.providerVersion = request.providerVersion;
    result.geometry.cacheKey = result.cacheKey;
    result.geometry.activeLevel = -1;
    result.counts.faceCount = 1;
    result.counts.pointCount = 3;
    result.bounds = request.bounds;
    result.mesh.points.push_back(SbVec3f(0.0f, 0.0f, 0.0f));
    result.mesh.points.push_back(SbVec3f(1.0f, 0.0f, 0.0f));
    result.mesh.points.push_back(SbVec3f(0.0f, 1.0f, 0.0f));
    result.mesh.coordIndex.push_back(0);
    result.mesh.coordIndex.push_back(1);
    result.mesh.coordIndex.push_back(2);

    return result;
}

static BRLObolLodResult
source_full_detail_task(const BRLObolLodRequest &request, void *UNUSED(userData))
{
    return source_full_detail_result(request);
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
	    service.queuedResultCountForDiagnostics() >= 1)
	    return 0;
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    printf("FAIL: LoD service did not produce update-action result\n");
    return 1;
}

static int
wait_for_service_queued(BRLObolLodService &service,
			size_t minQueuedCount,
			SbBool requireDelayedTask,
			SbBool requireNoInFlight)
{
    for (int i = 0; i < 400; i++) {
	if (service.queuedResultCountForDiagnostics() >= minQueuedCount &&
	    (!requireDelayedTask ||
	     service.delayedTaskCountForDiagnostics() > 0) &&
	    (!requireNoInFlight || service.inFlightCount() == 0))
	    return 0;
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    printf("FAIL: LoD service did not reach expected queued result state\n");
    return 1;
}

static int
submit_source_full_detail_task(BRLObolLodService &service,
			       const BRLObolLodRequest &request)
{
    BRLObolLodTask task;

    task.generation = service.beginGeneration();
    task.request = request;
    task.realize = source_full_detail_task;
    if (service.submit(task) == 0 || wait_for_service(service))
	return 1;

    return 0;
}

static int
expected_view_lod_level(struct db_i *dbip,
			const char *name,
			const struct bv_view_info *view,
			int *level)
{
    if (!dbip || !name || !view || !level)
	return 1;

    struct BRLObolMeshLod *lod = brlobol_mesh_lod_get(dbip, name);
    if (!lod)
	return 1;

    struct BRLObolMeshLodInfo info = BRLOBOL_MESH_LOD_INFO_INIT;
    int ret = 0;
    if (brlobol_mesh_lod_load_view(lod, view, 0) < 0 ||
	!brlobol_mesh_lod_info_get(lod, &info) ||
	info.active_level < 0) {
	ret = 1;
    } else {
	*level = info.active_level;
    }

    brlobol_mesh_lod_destroy(lod);
    return ret;
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
make_rt_pick_test_db(char *dbpath, size_t dbpath_len, struct db_i **dbip_out)
{
    if (!dbpath || dbpath_len == 0 || !dbip_out)
	return 1;
    *dbip_out = NULL;

    FILE *fp = bu_temp_file(dbpath, dbpath_len);
    if (!fp) {
	printf("FAIL: RT exact pick temp file\n");
	return 1;
    }
    fclose(fp);

    struct db_i *dbip = db_create(dbpath, 5);
    if (!dbip) {
	printf("FAIL: RT exact pick db_create\n");
	bu_file_delete(dbpath);
	return 1;
    }

    struct rt_wdb *wdbp = wdb_dbopen(dbip, RT_WDB_TYPE_DB_DISK);
    if (!wdbp) {
	printf("FAIL: RT exact pick wdb_dbopen\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    point_t center = {0.0, 0.0, 0.0};
    if (mk_sph(wdbp, "implicit.s", center, 2.0) != 0) {
	printf("FAIL: RT exact pick mk_sph\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    struct wmember region;
    BU_LIST_INIT(&region.l);
    unsigned char color[3] = {32, 96, 192};
    if (!mk_addmember("implicit.s", &region.l, NULL, WMOP_UNION) ||
	mk_lrcomb(wdbp, "implicit.r", &region, 1, "plastic", "",
		  color, 77, 2, 33, 100, 0) != 0) {
	printf("FAIL: RT exact pick mk_lrcomb\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    *dbip_out = dbip;
    return 0;
}

static int
test_rt_exact_pick_provider(void)
{
    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;
    if (make_rt_pick_test_db(dbpath, sizeof(dbpath), &dbip))
	return 1;

    std::vector<SbString> paths;
    paths.push_back("implicit.r");

    BRLObolRtPickCache pickCache;
    if (!pickCache.prepare(dbip, paths) ||
	!pickCache.isReady() ||
	pickCache.getObjectPathCount() != 1) {
	printf("FAIL: RT exact pick provider did not prepare reusable pick cache\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    BRLObolRtPickResult pick;
    if (!pickCache.pickRay(pick,
			   SbVec3f(0.0f, 0.0f, 5.0f),
			   SbVec3f(0.0f, 0.0f, -1.0f)) ||
	!pick.hit ||
	fabsf(pick.distance - 3.0f) > 1.0e-5f ||
	fabsf(pick.point[2] - 2.0f) > 1.0e-5f ||
	pick.normal[2] < 0.9f ||
	pick.detail.getPrimitiveKind() != SoBRLPickDetail::IMPLICIT_SOLID ||
	pick.detail.getRegionId() != 77 ||
	pick.detail.getAirCode() != 2 ||
	pick.detail.getMaterialId() != 33 ||
	pick.detail.getLos() != 100 ||
	strcmp(pick.detail.getSourceName().getString(), "implicit.s") != 0 ||
	strcmp(pick.detail.getSourceType().getString(), "sph") != 0 ||
	!strstr(pick.detail.getPath().getString(), "implicit.r")) {
	printf("FAIL: RT exact pick provider did not return implicit comb hit identity\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    BRLObolRtPickResult miss;
    if (pickCache.pickRay(miss,
			  SbVec3f(10.0f, 10.0f, 5.0f),
			  SbVec3f(0.0f, 0.0f, -1.0f)) || miss.hit) {
	printf("FAIL: RT exact pick provider reported a miss ray as a hit\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    BRLObolRtPickResult wrapperPick;
    if (!brlobol_pick_rt_ray(wrapperPick, dbip, paths,
			     SbVec3f(0.0f, 0.0f, 5.0f),
			     SbVec3f(0.0f, 0.0f, -1.0f)) ||
	!wrapperPick.hit ||
	fabsf(wrapperPick.distance - pick.distance) > 1.0e-5f) {
	printf("FAIL: RT exact pick provider one-shot wrapper did not use cache-backed ray path\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    pickCache.clear();
    if (pickCache.isReady() || pickCache.getObjectPathCount() != 0) {
	printf("FAIL: RT exact pick provider reusable pick cache did not clear\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    db_close(dbip);
    bu_file_delete(dbpath);
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

    SoBRLLodUpdateAction missingViewStateUpdate;
    missingViewStateUpdate.addResult(resultD);
    missingViewStateUpdate.apply(root);
    if (missingViewStateUpdate.getMatchedResultCount() != 1 ||
	missingViewStateUpdate.getAppliedResultCount() != 0 ||
	missingViewStateUpdate.getRejectedResultCount() != 1 ||
	missingViewStateUpdate.getUnmatchedResultCount() != 0 ||
	missingViewStateUpdate.getDiagnostics().getLength() == 0 ||
	meshD->lodAvailable.getValue() ||
	meshD->isLodDisplayActive() ||
	meshD->point.getNum() != 3 ||
	meshD->coordIndex.getNum() != 3 ||
	meshD->getTriangleCount() != 1) {
	printf("FAIL: LoD update action without view state mutated mesh or did not reject result\n");
	root->unref();
	return 1;
    }

    BRLObolViewLodState viewState;
    SoBRLViewLodGroup *renderRoot = new SoBRLViewLodGroup;
    renderRoot->ref();
    renderRoot->setViewLodState(&viewState);
    renderRoot->addChild(root);

    SoBRLLodUpdateAction update;
    update.setViewLodState(&viewState);
    update.addResult(resultA);
    update.addResult(resultB);
    update.addResult(resultC);
    update.addResult(resultD);
    update.addResult(resultMissing);
    update.apply(root);

    const BRLObolViewLodState::MeshPayload *meshDPayload =
	viewState.findMesh(meshD);
    if (update.getResultCount() != 5 ||
	update.getMatchedResultCount() != 4 ||
	update.getAppliedResultCount() != 1 ||
	update.getRejectedResultCount() != 3 ||
	update.getUnmatchedResultCount() != 1 ||
	update.getDiagnostics().getLength() == 0 ||
	!meshDPayload ||
	meshDPayload->activeLevel != 2 ||
	meshDPayload->getTriangleCount() != 2) {
	printf("FAIL: LoD update action view-local counts or payload\n");
	renderRoot->unref();
	root->unref();
	return 1;
    }

    if (meshA->lodStagedAvailable.getValue() ||
	meshB->lodStagedAvailable.getValue() ||
	meshC->lodStagedAvailable.getValue() ||
	meshD->lodStagedAvailable.getValue() ||
	meshD->lodAvailable.getValue() ||
	meshD->isLodDisplayActive() ||
	meshD->point.getNum() != 3 ||
	meshD->coordIndex.getNum() != 3 ||
	meshD->getTriangleCount() != 1) {
	printf("FAIL: LoD update action view-local update mutated shared mesh fields\n");
	renderRoot->unref();
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
			  "brlobol_mesh_lod",
			  "brlobol-cache-v1",
			  BRLOBOL_LOD_QUALITY_FAST_DISPLAY);
    if (meshDRequest.sourceCounts.faceCount != 1 ||
	meshDRequest.sourceCounts.pointCount != 3 ||
	meshDRequest.bounds.isEmpty() ||
	meshDRequest.bounds.getMax()[0] > 1.0f ||
	meshDRequest.bounds.getMax()[1] > 1.0f) {
	printf("FAIL: LoD-backed mesh request did not preserve full-detail source identity\n");
	renderRoot->unref();
	root->unref();
	return 1;
    }

    SoBRLExportAction exactExport;
    exactExport.apply(renderRoot);
    if (exactExport.getGeometryPolicy() != SoBRLExportAction::FULL_DETAIL ||
	exactExport.getTriangleCount() != 4 ||
	exactExport.getSkippedLodDisplayMeshCount() != 1 ||
	exactExport.getSourceBackedFullDetailRequestCount() != 0) {
	printf("FAIL: default export did not use preserved full-detail mesh\n");
	renderRoot->unref();
	root->unref();
	return 1;
    }

    SoBRLExportAction displayExport;
    displayExport.setGeometryPolicy(SoBRLExportAction::DISPLAY_LEVEL);
    displayExport.apply(renderRoot);
    if (displayExport.getGeometryPolicy() != SoBRLExportAction::DISPLAY_LEVEL ||
	displayExport.getTriangleCount() != 5 ||
	displayExport.getSkippedLodDisplayMeshCount() != 0) {
	printf("FAIL: display-level export did not include active display LoD mesh\n");
	renderRoot->unref();
	root->unref();
	return 1;
    }

    SoBRLMeasureAction exactMeasure;
    exactMeasure.apply(renderRoot);
    if (exactMeasure.getGeometryPolicy() != SoBRLMeasureAction::FULL_DETAIL ||
	exactMeasure.getTriangleCount() != 4 ||
	exactMeasure.getSkippedLodDisplayMeshCount() != 1 ||
	exactMeasure.getSourceBackedFullDetailRequestCount() != 0 ||
	fabsf(exactMeasure.getSurfaceArea() - 2.0f) > 1.0e-5f) {
	printf("FAIL: default measure did not use preserved full-detail mesh\n");
	renderRoot->unref();
	root->unref();
	return 1;
    }

    SoBRLMeasureAction displayMeasure;
    displayMeasure.setGeometryPolicy(SoBRLMeasureAction::DISPLAY_LEVEL);
    displayMeasure.apply(renderRoot);
    if (displayMeasure.getGeometryPolicy() != SoBRLMeasureAction::DISPLAY_LEVEL ||
	displayMeasure.getTriangleCount() != 5 ||
	displayMeasure.getSkippedLodDisplayMeshCount() != 0 ||
	fabsf(displayMeasure.getSurfaceArea() - 5.5f) > 1.0e-5f) {
	printf("FAIL: display-level measure did not include active display LoD mesh\n");
	renderRoot->unref();
	root->unref();
	return 1;
    }

    SoBRLSnapAction displaySnap;
    displaySnap.setEnabledKinds(SoBRLSnapAction::FACE_NEAREST);
    displaySnap.setQueryPoint(SbVec3f(1.5f, 0.2f, 0.0f));
    displaySnap.setTolerance(0.1f);
    displaySnap.apply(renderRoot);
    if (displaySnap.getGeometryPolicy() != SoBRLSnapAction::DISPLAY_LEVEL ||
	!displaySnap.hasCandidate() ||
	strcmp(displaySnap.getPath().getString(), "/mesh/d") != 0 ||
	displaySnap.getSkippedLodDisplayMeshCount() != 0) {
	printf("FAIL: display-level snap did not use active display LoD mesh\n");
	renderRoot->unref();
	root->unref();
	return 1;
    }

    SoBRLSnapAction exactSnap;
    exactSnap.setGeometryPolicy(SoBRLSnapAction::FULL_DETAIL);
    exactSnap.setEnabledKinds(SoBRLSnapAction::FACE_NEAREST);
    exactSnap.setQueryPoint(SbVec3f(1.5f, 0.2f, 0.0f));
    exactSnap.setTolerance(0.1f);
    exactSnap.apply(renderRoot);
    if (exactSnap.getGeometryPolicy() != SoBRLSnapAction::FULL_DETAIL ||
	exactSnap.hasCandidate() ||
	exactSnap.getSkippedLodDisplayMeshCount() != 1 ||
	exactSnap.getSourceBackedFullDetailRequestCount() != 0) {
	printf("FAIL: full-detail snap did not skip active display LoD mesh\n");
	renderRoot->unref();
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
    exactFullSnap.apply(renderRoot);
    if (!exactFullSnap.hasCandidate() ||
	strcmp(exactFullSnap.getPath().getString(), "/mesh/d") != 0 ||
	exactFullSnap.getSkippedLodDisplayMeshCount() != 1 ||
	exactFullSnap.getSourceBackedFullDetailRequestCount() != 0) {
	printf("FAIL: full-detail snap did not use preserved full-detail mesh\n");
	renderRoot->unref();
	root->unref();
	return 1;
    }

    if (meshD->getPickGeometryPolicy() != SoBRLMeshShape::PICK_DISPLAY_LEVEL) {
	printf("FAIL: mesh pick policy should default to display-level geometry\n");
	renderRoot->unref();
	root->unref();
	return 1;
    }

    SbViewportRegion pickViewport(100, 100);
    SoRayPickAction displayPick(pickViewport);
    displayPick.setRay(SbVec3f(1.5f, 0.2f, 5.0f),
		       SbVec3f(0.0f, 0.0f, -1.0f));
    displayPick.apply(renderRoot);
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
	renderRoot->unref();
	root->unref();
	return 1;
    }

    meshD->setPickGeometryPolicy(SoBRLMeshShape::PICK_FULL_DETAIL);
    SoRayPickAction exactMissPick(pickViewport);
    exactMissPick.setRay(SbVec3f(1.5f, 0.2f, 5.0f),
			 SbVec3f(0.0f, 0.0f, -1.0f));
    exactMissPick.apply(renderRoot);
    if (exactMissPick.getPickedPoint()) {
	printf("FAIL: full-detail pick did not skip active display LoD mesh\n");
	renderRoot->unref();
	root->unref();
	return 1;
    }

    SoRayPickAction exactHitPick(pickViewport);
    exactHitPick.setRay(SbVec3f(0.2f, 0.2f, 5.0f),
			SbVec3f(0.0f, 0.0f, -1.0f));
    exactHitPick.apply(renderRoot);
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
	renderRoot->unref();
	root->unref();
	return 1;
    }
    meshD->setPickGeometryPolicy(SoBRLMeshShape::PICK_DISPLAY_LEVEL);

    renderRoot->unref();
    root->unref();
    return 0;
}

static int
test_scene_database_source_summary(void)
{
    SoSeparator *root = new SoSeparator;
    root->ref();

    SoSeparator *nonSourceChild = new SoSeparator;
    root->addChild(nonSourceChild);

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->path = "/summary/source";
    source->sourceRevision = 11;
    source->inputsRevision = 2;
    source->viewRevision = 3;
    source->highlighted = TRUE;
    source->lineStyle = 2;
    source->lineWidth = 3;
    source->transparency = 0.2f;
    source->materialColorValid = TRUE;
    source->materialColor = SbColor(0.1f, 0.2f, 0.3f);
    source->materialRevision = 55;
    root->addChild(source);

    SoBRLSceneController scene(root);
    const uint64_t initialStructuralRevision = scene.getStructuralRevision();
    const uint64_t initialFrameRevision = scene.getFrameRevision();
    BRLObolSceneSummary sceneSummary;
    BRLObolDatabaseSourceSummary summary;
    if (scene.getDatabaseSourceCount() != 1 ||
	initialStructuralRevision == 0 ||
	initialFrameRevision == 0 ||
	!scene.getSceneSummary(sceneSummary) ||
	!sceneSummary.valid ||
	!sceneSummary.hasRoot ||
	!sceneSummary.rootIsGroup ||
	sceneSummary.rootChildCount != 2 ||
	sceneSummary.databaseSourceCount != 1 ||
	sceneSummary.nonDatabaseRootChildCount != 1 ||
	sceneSummary.structuralRevision != initialStructuralRevision ||
	sceneSummary.frameRevision != initialFrameRevision ||
	scene.getDatabaseSource(0) != source ||
	scene.getDatabaseSource(1) != NULL ||
	!scene.getDatabaseSourceSummary(0, summary) ||
	!summary.valid ||
	strcmp(summary.path.getString(), "/summary/source") != 0 ||
	!summary.hasParent ||
	summary.drawTreeDepth != 1 ||
	strcmp(summary.parentGroupPath.getString(), "/") != 0 ||
	summary.sourceRevision != 11 ||
	summary.inputsRevision != 2 ||
	summary.viewRevision != 3 ||
	!summary.stale ||
	!(summary.staleReason & SoBRLDatabaseSource::STALE_SOURCE) ||
	summary.realizationStatus != SoBRLDatabaseSource::UNREALIZED ||
	summary.realizationIdentity.getLength() != 0 ||
	summary.realizedShapeCount != 0 ||
	summary.realizedMeshCount != 0) {
	printf("FAIL: scene controller should summarize initial database source state\n");
	root->unref();
	return 1;
    }

    {
	BRLObolViewController view(root, NULL);
	BRLObolDatabaseSourceSummary viewSummary;
	if (view.getDatabaseSourceCount() != 1 ||
	    view.getDatabaseSource(0) != source ||
	    !view.getDatabaseSourceSummary(0, viewSummary) ||
	    !viewSummary.hasParent ||
	    strcmp(viewSummary.parentGroupPath.getString(), "/") != 0) {
	    printf("FAIL: view controller should delegate database source lookup to scene controller\n");
	    root->unref();
	    return 1;
	}
    }

    if (!scene.realizePending() ||
	scene.getStructuralRevision() != initialStructuralRevision ||
	scene.getFrameRevision() == initialFrameRevision ||
	!scene.getSceneSummary(sceneSummary) ||
	sceneSummary.rootChildCount != 2 ||
	sceneSummary.databaseSourceCount != 1 ||
	sceneSummary.nonDatabaseRootChildCount != 1 ||
	sceneSummary.structuralRevision != scene.getStructuralRevision() ||
	sceneSummary.frameRevision != scene.getFrameRevision() ||
	sceneSummary.lastVisitedSourceCount != 1 ||
	sceneSummary.lastRealizedSourceCount != 1 ||
	sceneSummary.lastFailedSourceCount != 0 ||
	!scene.getDatabaseSourceSummary(0, summary) ||
	summary.stale ||
	summary.staleReason != SoBRLDatabaseSource::STALE_NONE ||
	summary.realizationStatus != SoBRLDatabaseSource::REALIZED ||
	summary.realizedSourceRevision != 11 ||
	summary.realizedInputsRevision != 2 ||
	summary.realizedViewRevision != 3 ||
	summary.realizedShapeCount != 1 ||
	summary.realizedMeshCount != 0 ||
	summary.realizationIdentity.getLength() == 0 ||
	!strstr(summary.realizationIdentity.getString(),
		"source=11;inputs=2;view=3")) {
	printf("FAIL: scene controller should summarize realized database source state\n");
	root->unref();
	return 1;
    }
    SbString firstIdentity = summary.realizationIdentity;
    const uint64_t realizedFrameRevision = scene.getFrameRevision();
    SoBRLVListShape *realizedShape = source->getRealizedShape();
    if (!realizedShape ||
	strcmp(realizedShape->ownerSourcePath.getValue().getString(),
	       "/summary/source") != 0 ||
	realizedShape->ownerSourceRevision.getValue() != 11 ||
	realizedShape->ownerInputsRevision.getValue() != 2 ||
	realizedShape->ownerViewRevision.getValue() != 3 ||
	realizedShape->ownerRealizedSourceRevision.getValue() != 11 ||
	realizedShape->ownerRealizedInputsRevision.getValue() != 2 ||
	realizedShape->ownerRealizedViewRevision.getValue() != 3 ||
	realizedShape->ownerRealizationStatus.getValue() !=
	SoBRLDatabaseSource::REALIZED ||
	strcmp(realizedShape->ownerRealizationIdentity.getValue().getString(),
	       firstIdentity.getString()) != 0 ||
	realizedShape->ownerSourceStale.getValue() ||
	realizedShape->ownerStaleReason.getValue() !=
	SoBRLDatabaseSource::STALE_NONE) {
	printf("FAIL: realized shape node should mirror source owner state after realization\n");
	root->unref();
	return 1;
    }

    source->inputsRevision = 4;
    if (!scene.getDatabaseSourceSummary(0, summary) ||
	scene.getFrameRevision() != realizedFrameRevision ||
	!summary.stale ||
	!(summary.staleReason & SoBRLDatabaseSource::STALE_INPUTS) ||
	summary.realizationStatus != SoBRLDatabaseSource::UNREALIZED ||
	strcmp(summary.realizationIdentity.getString(),
	       firstIdentity.getString()) != 0 ||
	!realizedShape->ownerSourceStale.getValue() ||
	!(realizedShape->ownerStaleReason.getValue() &
	  SoBRLDatabaseSource::STALE_INPUTS) ||
	realizedShape->ownerRealizationStatus.getValue() !=
	SoBRLDatabaseSource::UNREALIZED ||
	strcmp(realizedShape->ownerRealizationIdentity.getValue().getString(),
	       firstIdentity.getString()) != 0) {
	printf("FAIL: scene controller source summary should expose stale source state before refresh\n");
	root->unref();
	return 1;
    }

    if (!scene.realizePending() ||
	scene.getStructuralRevision() != initialStructuralRevision ||
	scene.getFrameRevision() == realizedFrameRevision ||
	!scene.getSceneSummary(sceneSummary) ||
	sceneSummary.rootChildCount != 2 ||
	sceneSummary.databaseSourceCount != 1 ||
	sceneSummary.structuralRevision != scene.getStructuralRevision() ||
	sceneSummary.frameRevision != scene.getFrameRevision() ||
	sceneSummary.lastVisitedSourceCount != 1 ||
	sceneSummary.lastRealizedSourceCount != 1 ||
	sceneSummary.lastFailedSourceCount != 0 ||
	!scene.getDatabaseSourceSummary(0, summary) ||
	summary.stale ||
	summary.realizedInputsRevision != 4 ||
	summary.realizationIdentity.getLength() == 0 ||
	strcmp(summary.realizationIdentity.getString(),
	       firstIdentity.getString()) == 0 ||
	!strstr(summary.realizationIdentity.getString(),
		"source=11;inputs=4;view=3")) {
	printf("FAIL: scene controller source summary should update realization identity after refresh\n");
	root->unref();
	return 1;
    }
    realizedShape = source->getRealizedShape();
    if (!realizedShape ||
	realizedShape->ownerInputsRevision.getValue() != 4 ||
	realizedShape->ownerRealizedInputsRevision.getValue() != 4 ||
	realizedShape->ownerRealizationStatus.getValue() !=
	SoBRLDatabaseSource::REALIZED ||
	strcmp(realizedShape->ownerRealizationIdentity.getValue().getString(),
	       summary.realizationIdentity.getString()) != 0 ||
	realizedShape->ownerSourceStale.getValue() ||
	realizedShape->ownerStaleReason.getValue() !=
	SoBRLDatabaseSource::STALE_NONE) {
	printf("FAIL: realized shape node should refresh mirrored source owner state\n");
	root->unref();
	return 1;
    }

    BRLObolRealizedShapeSummary shapeSummary;
    if (source->getRealizedShapeSummaryCount() != 1 ||
	scene.getRealizedShapeSummaryCount() != 1 ||
	!source->getRealizedShapeSummary(0, shapeSummary) ||
	!shapeSummary.valid ||
	shapeSummary.shapeKind != BRLObolRealizedShapeSummary::SHAPE_VLIST ||
	strcmp(shapeSummary.path.getString(), "/summary/source") != 0 ||
	shapeSummary.ownerSourceIndex != -1 ||
	strcmp(shapeSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	shapeSummary.ownerDrawMode != source->drawMode.getValue() ||
	shapeSummary.ownerSourceRevision != 11 ||
	shapeSummary.ownerInputsRevision != 4 ||
	shapeSummary.ownerViewRevision != 3 ||
	shapeSummary.ownerRealizedRevision != summary.realizedRevision ||
	shapeSummary.ownerRealizedSourceRevision != 11 ||
	shapeSummary.ownerRealizedInputsRevision != 4 ||
	shapeSummary.ownerRealizedViewRevision != 3 ||
	shapeSummary.ownerRealizationStatus !=
	SoBRLDatabaseSource::REALIZED ||
	shapeSummary.ownerRealizationDiagnostic.getLength() != 0 ||
	strcmp(shapeSummary.ownerRealizationIdentity.getString(),
	       summary.realizationIdentity.getString()) != 0 ||
	shapeSummary.ownerSourceStale ||
	shapeSummary.ownerStaleReason != SoBRLDatabaseSource::STALE_NONE ||
	strcmp(shapeSummary.sourceName.getString(), "summary/source") != 0 ||
	strcmp(shapeSummary.sourceType.getString(), "prototype") != 0 ||
	strcmp(shapeSummary.displayName.getString(), "summary/source") != 0 ||
	!shapeSummary.localSource ||
	!shapeSummary.nonDatabaseSource ||
	shapeSummary.databaseIntent ||
	shapeSummary.drawMode != BRLOBOL_LOD_DRAW_DIAGNOSTIC ||
	strcmp(shapeSummary.recordRole.getString(), "prototype") != 0 ||
	!shapeSummary.visible ||
	!shapeSummary.selectable ||
	shapeSummary.pointCount != 5 ||
	shapeSummary.commandCount != 5 ||
	shapeSummary.segmentCount != 4 ||
	shapeSummary.pointPrimitiveCount != 0 ||
	!shapeSummary.boundsValid) {
	printf("FAIL: database source should summarize realized vlist shape metadata\n");
	root->unref();
	return 1;
    }
    BRLObolRealizedShapeSummary sceneShapeSummary;
    if (!scene.getRealizedShapeSummary(0, sceneShapeSummary) ||
	!sceneShapeSummary.valid ||
	sceneShapeSummary.shapeKind !=
	BRLObolRealizedShapeSummary::SHAPE_VLIST ||
	strcmp(sceneShapeSummary.path.getString(),
	       shapeSummary.path.getString()) != 0 ||
	sceneShapeSummary.ownerSourceIndex != 0 ||
	strcmp(sceneShapeSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	strcmp(sceneShapeSummary.ownerRealizationIdentity.getString(),
	       summary.realizationIdentity.getString()) != 0 ||
	scene.getRealizedShapeSummary(1, sceneShapeSummary)) {
	printf("FAIL: scene controller should forward realized vlist shape summaries\n");
	root->unref();
	return 1;
    }

    SbVec3f auxPoints[2] = {
	SbVec3f(1.0f, 2.0f, 3.0f),
	SbVec3f(4.0f, 5.0f, 6.0f)
    };
    int32_t auxCommands[2] = {
	SoBRLVListShape::MOVE,
	SoBRLVListShape::DRAW
    };
    BRLObolAuxiliaryLineSetDisplayState auxDisplay;
    auxDisplay.valid = TRUE;
    auxDisplay.drawMode = BRLOBOL_LOD_DRAW_SHADED;
    auxDisplay.visible = FALSE;
    auxDisplay.highlighted = TRUE;
    auxDisplay.lineStyle = 7;
    auxDisplay.lineWidth = 13;
    auxDisplay.transparency = 0.625f;
    auxDisplay.materialColorValid = TRUE;
    auxDisplay.materialColor = SbColor(0.25f, 0.5f, 0.75f);
    auxDisplay.materialRevision = 88;
    if (scene.publishDatabaseSourceAuxiliaryLineSet("/summary/source",
	    "summary_aux", auxPoints, auxCommands, 2, &auxDisplay) != 1) {
	printf("FAIL: scene controller should publish auxiliary source line sets\n");
	root->unref();
	return 1;
    }
    SoBRLVListShape *auxShape =
	source->findAuxiliaryVListShape("summary_aux");
    if (!auxShape ||
	auxShape->point.getNum() != 2 ||
	auxShape->command.getNum() != 2 ||
	strcmp(auxShape->recordRole.getValue().getString(),
	       "auxiliary") != 0 ||
	source->getRealizedShapeSummaryCount() != 2 ||
	scene.getRealizedShapeSummaryCount() != 2 ||
	!source->getRealizedShapeSummary(1, shapeSummary) ||
	!shapeSummary.valid ||
	shapeSummary.shapeKind !=
	BRLObolRealizedShapeSummary::SHAPE_VLIST ||
	strcmp(shapeSummary.sourceName.getString(), "summary_aux") != 0 ||
	strcmp(shapeSummary.sourceType.getString(),
	       "auxiliary-line-set") != 0 ||
	strcmp(shapeSummary.geometryName.getString(), "summary_aux") != 0 ||
	strcmp(shapeSummary.recordRole.getString(), "auxiliary") != 0 ||
	!shapeSummary.databaseIntent ||
	shapeSummary.nonDatabaseSource ||
	shapeSummary.drawMode != BRLOBOL_LOD_DRAW_SHADED ||
	shapeSummary.visible ||
	!shapeSummary.highlighted ||
	shapeSummary.materialRevision != 88 ||
	!shapeSummary.materialColorValid ||
	fabsf(shapeSummary.materialColor[0] - 0.25f) > 1.0e-6f ||
	fabsf(shapeSummary.materialColor[1] - 0.5f) > 1.0e-6f ||
	fabsf(shapeSummary.materialColor[2] - 0.75f) > 1.0e-6f ||
	shapeSummary.pointCount != 2 ||
	shapeSummary.segmentCount != 1 ||
	strcmp(shapeSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0) {
	printf("FAIL: database source should summarize auxiliary line-set metadata\n");
	root->unref();
	return 1;
    }
    if (scene.markDatabaseSourceStale("/summary/source",
				      SoBRLDatabaseSource::STALE_DRAW) != 1 ||
	!scene.realizePending() ||
	!source->findAuxiliaryVListShape("summary_aux") ||
	source->getRealizedShapeSummaryCount() != 2 ||
	scene.getRealizedShapeSummaryCount() != 2) {
	printf("FAIL: database source realization should preserve auxiliary line-set ownership\n");
	root->unref();
	return 1;
    }
    if (scene.clearDatabaseSourceAuxiliaryShapes("/summary/source") != 1 ||
	source->findAuxiliaryVListShape("summary_aux") ||
	source->getRealizedShapeSummaryCount() != 1 ||
	scene.getRealizedShapeSummaryCount() != 1) {
	printf("FAIL: scene controller should clear auxiliary source line sets\n");
	root->unref();
	return 1;
    }

    SoBRLMeshShape *summaryMesh = make_mesh("/summary/mesh", "mesh");
    summaryMesh->sourceType = "bot";
    summaryMesh->sourceId = 44;
    summaryMesh->displayName = "Summary Mesh";
    summaryMesh->geometryName = "summary mesh geometry";
    summaryMesh->cacheIdentity = "cache://summary/mesh#44";
    summaryMesh->sourceIdentity = "db://summary/mesh";
    summaryMesh->databaseIntent = TRUE;
    summaryMesh->localSource = FALSE;
    summaryMesh->nonDatabaseSource = FALSE;
    summaryMesh->drawMode = BRLOBOL_LOD_DRAW_SHADED;
    summaryMesh->recordRole = "database";
    summaryMesh->geometryKind = "surface";
    summaryMesh->highlighted = TRUE;
    summaryMesh->lineStyle = 1;
    summaryMesh->lineWidth = 5;
    summaryMesh->transparency = 0.35f;
    summaryMesh->materialColorValid = TRUE;
    summaryMesh->materialColor = SbColor(0.25f, 0.5f, 0.75f);
    summaryMesh->materialRevision = 77;
    source->addChild(summaryMesh);
    if (source->getRealizedShapeSummaryCount() != 2 ||
	scene.getRealizedShapeSummaryCount() != 2 ||
	!source->getRealizedShapeSummary(1, shapeSummary) ||
	!shapeSummary.valid ||
	shapeSummary.shapeKind != BRLObolRealizedShapeSummary::SHAPE_MESH ||
	strcmp(shapeSummary.path.getString(), "/summary/mesh") != 0 ||
	shapeSummary.ownerSourceIndex != -1 ||
	strcmp(shapeSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	shapeSummary.ownerInputsRevision != 4 ||
	shapeSummary.ownerRealizationStatus !=
	SoBRLDatabaseSource::REALIZED ||
	strcmp(shapeSummary.sourceName.getString(), "mesh") != 0 ||
	strcmp(shapeSummary.sourceType.getString(), "bot") != 0 ||
	strcmp(shapeSummary.displayName.getString(), "Summary Mesh") != 0 ||
	strcmp(shapeSummary.cacheIdentity.getString(),
	       "cache://summary/mesh#44") != 0 ||
	strcmp(shapeSummary.sourceIdentity.getString(),
	       "db://summary/mesh") != 0 ||
	!shapeSummary.databaseIntent ||
	shapeSummary.localSource ||
	shapeSummary.nonDatabaseSource ||
	shapeSummary.drawMode != BRLOBOL_LOD_DRAW_SHADED ||
	strcmp(shapeSummary.recordRole.getString(), "database") != 0 ||
	strcmp(shapeSummary.geometryKind.getString(), "surface") != 0 ||
	shapeSummary.materialRevision != 77 ||
	shapeSummary.pointCount != 3 ||
	shapeSummary.indexCount != 3 ||
	shapeSummary.triangleCount != 1 ||
	!shapeSummary.boundsValid ||
	source->getRealizedShapeSummary(2, shapeSummary)) {
	printf("FAIL: database source should summarize realized mesh shape metadata\n");
	root->unref();
	return 1;
    }
    if (!scene.getRealizedShapeSummary(1, sceneShapeSummary) ||
	!sceneShapeSummary.valid ||
	sceneShapeSummary.shapeKind !=
	BRLObolRealizedShapeSummary::SHAPE_MESH ||
	strcmp(sceneShapeSummary.path.getString(), "/summary/mesh") != 0 ||
	sceneShapeSummary.ownerSourceIndex != 0 ||
	strcmp(sceneShapeSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	scene.getRealizedShapeSummary(2, sceneShapeSummary)) {
	printf("FAIL: scene controller should forward realized mesh shape summaries\n");
	root->unref();
	return 1;
    }

    SoBRLMaterialObject *summaryMaterial = new SoBRLMaterialObject;
    summaryMaterial->sourcePath = "/summary/material";
    summaryMaterial->sourceName = "mat";
    summaryMaterial->sourceType = "material";
    summaryMaterial->sourceId = 77;
    summaryMaterial->materialName = "matte";
    summaryMaterial->parentName = "base";
    summaryMaterial->materialSource = "database";
    summaryMaterial->addProperty("physical", "density", "7.8");
    summaryMaterial->addProperty("optical", "finish", "matte");
    source->addChild(summaryMaterial);

    BRLObolRealizedMaterialSummary materialSummary;
    SbString propertyGroup;
    SbString propertyName;
    SbString propertyValue;
    if (source->getRealizedMaterialSummaryCount() != 1 ||
	scene.getRealizedMaterialSummaryCount() != 1 ||
	!source->getRealizedMaterialSummary(0, materialSummary) ||
	!materialSummary.valid ||
	strcmp(materialSummary.sourcePath.getString(),
	       "/summary/material") != 0 ||
	strcmp(materialSummary.sourceName.getString(), "mat") != 0 ||
	strcmp(materialSummary.sourceType.getString(), "material") != 0 ||
	materialSummary.sourceId != 77 ||
	materialSummary.ownerSourceIndex != -1 ||
	strcmp(materialSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	materialSummary.ownerInputsRevision != 4 ||
	materialSummary.ownerRealizationStatus !=
	SoBRLDatabaseSource::REALIZED ||
	strcmp(materialSummary.ownerRealizationIdentity.getString(),
	       summary.realizationIdentity.getString()) != 0 ||
	strcmp(materialSummary.materialName.getString(), "matte") != 0 ||
	strcmp(materialSummary.parentName.getString(), "base") != 0 ||
	strcmp(materialSummary.materialSource.getString(), "database") != 0 ||
	materialSummary.propertyCount != 2 ||
	!source->getRealizedMaterialProperty(0, 0, propertyGroup,
		propertyName, propertyValue) ||
	strcmp(propertyGroup.getString(), "physical") != 0 ||
	strcmp(propertyName.getString(), "density") != 0 ||
	strcmp(propertyValue.getString(), "7.8") != 0 ||
	source->getRealizedMaterialProperty(0, 2, propertyGroup,
					    propertyName, propertyValue) ||
	source->getRealizedMaterialSummary(1, materialSummary)) {
	printf("FAIL: database source should summarize realized material objects\n");
	root->unref();
	return 1;
    }

    BRLObolRealizedMaterialSummary sceneMaterialSummary;
    if (!scene.getRealizedMaterialSummary(0, sceneMaterialSummary) ||
	!sceneMaterialSummary.valid ||
	sceneMaterialSummary.ownerSourceIndex != 0 ||
	strcmp(sceneMaterialSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	strcmp(sceneMaterialSummary.sourcePath.getString(),
	       "/summary/material") != 0 ||
	!scene.getRealizedMaterialProperty(0, 1, propertyGroup,
					   propertyName, propertyValue) ||
	strcmp(propertyGroup.getString(), "optical") != 0 ||
	strcmp(propertyName.getString(), "finish") != 0 ||
	strcmp(propertyValue.getString(), "matte") != 0 ||
	scene.getRealizedMaterialProperty(1, 0, propertyGroup,
					  propertyName, propertyValue) ||
	scene.getRealizedMaterialSummary(1, sceneMaterialSummary)) {
	printf("FAIL: scene controller should forward realized material object summaries\n");
	root->unref();
	return 1;
    }

    BRLObolSceneTreeSummary treeSummary;
    if (source->getRealizedTreeSummaryCount() != 4 ||
	!source->getRealizedTreeSummary(0, treeSummary) ||
	!treeSummary.valid ||
	treeSummary.nodeKind !=
	BRLObolSceneTreeSummary::NODE_DATABASE_SOURCE ||
	!treeSummary.isGroup ||
	treeSummary.isShape ||
	treeSummary.hasParent ||
	treeSummary.drawTreeDepth != 0 ||
	treeSummary.childCount != 3 ||
	treeSummary.ownerSourceIndex != -1 ||
	strcmp(treeSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	strcmp(treeSummary.path.getString(), "/summary/source") != 0 ||
	!source->getRealizedTreeSummary(1, treeSummary) ||
	!treeSummary.valid ||
	treeSummary.nodeKind != BRLObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	!treeSummary.isShape ||
	!treeSummary.hasParent ||
	treeSummary.drawTreeDepth != 1 ||
	treeSummary.childCount != 0 ||
	strcmp(treeSummary.path.getString(), "/summary/source") != 0 ||
	!source->getRealizedTreeSummary(2, treeSummary) ||
	treeSummary.nodeKind != BRLObolSceneTreeSummary::NODE_MESH_SHAPE ||
	strcmp(treeSummary.path.getString(), "/summary/mesh") != 0 ||
	!source->getRealizedTreeSummary(3, treeSummary) ||
	treeSummary.nodeKind !=
	BRLObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
	!treeSummary.isMaterialObject ||
	strcmp(treeSummary.path.getString(), "/summary/material") != 0 ||
	strcmp(treeSummary.displayName.getString(), "matte") != 0 ||
	source->getRealizedTreeSummary(4, treeSummary)) {
	printf("FAIL: database source should summarize realized source tree nodes\n");
	root->unref();
	return 1;
    }

    if (scene.getSceneTreeSummaryCount() != 6 ||
	!scene.getSceneTreeSummary(0, treeSummary) ||
	!treeSummary.valid ||
	treeSummary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
	!treeSummary.isGroup ||
	treeSummary.hasParent ||
	treeSummary.drawTreeDepth != 0 ||
	treeSummary.childCount != 2 ||
	!scene.getSceneTreeSummary(1, treeSummary) ||
	treeSummary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
	!treeSummary.hasParent ||
	treeSummary.ownerSourceIndex != -1 ||
	treeSummary.drawTreeDepth != 1 ||
	!scene.getSceneTreeSummary(2, treeSummary) ||
	treeSummary.nodeKind !=
	BRLObolSceneTreeSummary::NODE_DATABASE_SOURCE ||
	treeSummary.ownerSourceIndex != 0 ||
	strcmp(treeSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	treeSummary.drawTreeDepth != 1 ||
	!treeSummary.hasParent ||
	!scene.getSceneTreeSummary(3, treeSummary) ||
	treeSummary.nodeKind != BRLObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	treeSummary.ownerSourceIndex != 0 ||
	treeSummary.drawTreeDepth != 2 ||
	!scene.getSceneTreeSummary(4, treeSummary) ||
	treeSummary.nodeKind != BRLObolSceneTreeSummary::NODE_MESH_SHAPE ||
	strcmp(treeSummary.path.getString(), "/summary/mesh") != 0 ||
	!scene.getSceneTreeSummary(5, treeSummary) ||
	treeSummary.nodeKind !=
	BRLObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
	strcmp(treeSummary.path.getString(), "/summary/material") != 0 ||
	scene.getSceneTreeSummary(6, treeSummary)) {
	printf("FAIL: scene controller should summarize Obol scene tree nodes\n");
	root->unref();
	return 1;
    }

    BRLObolSceneDisplaySummary displaySummary;
    if (source->getRealizedDisplaySummaryCount() !=
	source->getRealizedTreeSummaryCount() ||
	!source->getRealizedDisplaySummary(0, displaySummary) ||
	!displaySummary.valid ||
	!displaySummary.isDatabaseSource ||
	!displaySummary.hasDrawIntent ||
	strcmp(displaySummary.intentPath.getString(),
	       "/summary/source") != 0 ||
	displaySummary.intentDrawMode != source->drawMode.getValue() ||
	displaySummary.drawMode != source->drawMode.getValue() ||
	!displaySummary.visible ||
	!displaySummary.highlighted ||
	displaySummary.lineStyle != 2 ||
	displaySummary.lineWidth != 3 ||
	displaySummary.transparency < 0.19 ||
	displaySummary.transparency > 0.21 ||
	!displaySummary.materialValid ||
	displaySummary.materialRevision != 55 ||
	displaySummary.materialColor[0] < 0.09f ||
	displaySummary.materialColor[0] > 0.11f ||
	displaySummary.materialColor[1] < 0.19f ||
	displaySummary.materialColor[1] > 0.21f ||
	displaySummary.materialColor[2] < 0.29f ||
	displaySummary.materialColor[2] > 0.31f ||
	displaySummary.ownerSourceIndex != -1 ||
	strcmp(displaySummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	!source->getRealizedDisplaySummary(1, displaySummary) ||
	displaySummary.nodeKind != BRLObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	!displaySummary.hasDrawIntent ||
	strcmp(displaySummary.intentPath.getString(),
	       "/summary/source") != 0 ||
	displaySummary.intentDrawMode != BRLOBOL_LOD_DRAW_DIAGNOSTIC ||
	!displaySummary.visible ||
	displaySummary.highlighted ||
	displaySummary.materialColor[0] < 0.99f ||
	!source->getRealizedDisplaySummary(2, displaySummary) ||
	displaySummary.nodeKind != BRLObolSceneTreeSummary::NODE_MESH_SHAPE ||
	strcmp(displaySummary.intentPath.getString(), "/summary/mesh") != 0 ||
	!displaySummary.highlighted ||
	displaySummary.lineStyle != 1 ||
	displaySummary.lineWidth != 5 ||
	displaySummary.transparency < 0.34 ||
	displaySummary.transparency > 0.36 ||
	displaySummary.materialRevision != 77 ||
	!displaySummary.materialValid ||
	displaySummary.materialColor[0] < 0.24f ||
	displaySummary.materialColor[0] > 0.26f ||
	displaySummary.materialColor[1] < 0.49f ||
	displaySummary.materialColor[1] > 0.51f ||
	displaySummary.materialColor[2] < 0.74f ||
	displaySummary.materialColor[2] > 0.76f ||
	!source->getRealizedDisplaySummary(3, displaySummary) ||
	displaySummary.nodeKind !=
	BRLObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
	displaySummary.hasDrawIntent ||
	strcmp(displaySummary.path.getString(), "/summary/material") != 0 ||
	source->getRealizedDisplaySummary(4, displaySummary)) {
	printf("FAIL: database source should summarize realized display state\n");
	root->unref();
	return 1;
    }

    if (scene.getSceneDisplaySummaryCount() !=
	scene.getSceneTreeSummaryCount() ||
	!scene.getSceneDisplaySummary(0, displaySummary) ||
	!displaySummary.valid ||
	displaySummary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
	displaySummary.hasDrawIntent ||
	displaySummary.ownerSourceIndex != -1 ||
	!scene.getSceneDisplaySummary(1, displaySummary) ||
	displaySummary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
	displaySummary.hasDrawIntent ||
	displaySummary.ownerSourceIndex != -1 ||
	!scene.getSceneDisplaySummary(2, displaySummary) ||
	!displaySummary.isDatabaseSource ||
	displaySummary.ownerSourceIndex != 0 ||
	strcmp(displaySummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	strcmp(displaySummary.intentPath.getString(),
	       "/summary/source") != 0 ||
	!displaySummary.visible ||
	!displaySummary.highlighted ||
	displaySummary.lineStyle != 2 ||
	displaySummary.lineWidth != 3 ||
	displaySummary.transparency < 0.19 ||
	displaySummary.transparency > 0.21 ||
	!displaySummary.materialValid ||
	displaySummary.materialRevision != 55 ||
	displaySummary.materialColor[0] < 0.09f ||
	displaySummary.materialColor[0] > 0.11f ||
	displaySummary.materialColor[2] < 0.29f ||
	displaySummary.materialColor[2] > 0.31f ||
	!scene.getSceneDisplaySummary(4, displaySummary) ||
	displaySummary.nodeKind != BRLObolSceneTreeSummary::NODE_MESH_SHAPE ||
	displaySummary.ownerSourceIndex != 0 ||
	!displaySummary.highlighted ||
	displaySummary.lineStyle != 1 ||
	displaySummary.lineWidth != 5 ||
	displaySummary.transparency < 0.34 ||
	displaySummary.transparency > 0.36 ||
	displaySummary.materialRevision != 77 ||
	displaySummary.materialColor[2] < 0.74f ||
	displaySummary.materialColor[2] > 0.76f ||
	!scene.getSceneDisplaySummary(5, displaySummary) ||
	displaySummary.nodeKind !=
	BRLObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
	strcmp(displaySummary.path.getString(), "/summary/material") != 0 ||
	scene.getSceneDisplaySummary(6, displaySummary)) {
	printf("FAIL: scene controller should summarize Obol scene display state\n");
	root->unref();
	return 1;
    }

    BRLObolSceneMaterialSummary materialStateSummary;
    if (source->getRealizedSceneMaterialSummaryCount() !=
	source->getRealizedTreeSummaryCount() ||
	!source->getRealizedSceneMaterialSummary(0,
		materialStateSummary) ||
	!materialStateSummary.valid ||
	materialStateSummary.materialValid ||
	materialStateSummary.ownerSourceIndex != -1 ||
	!source->getRealizedSceneMaterialSummary(1,
		materialStateSummary) ||
	!materialStateSummary.valid ||
	materialStateSummary.nodeKind !=
	BRLObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	!materialStateSummary.materialValid ||
	materialStateSummary.materialRevision != 0 ||
	materialStateSummary.materialColor[0] < 0.99f ||
	!source->getRealizedSceneMaterialSummary(2,
		materialStateSummary) ||
	materialStateSummary.nodeKind !=
	BRLObolSceneTreeSummary::NODE_MESH_SHAPE ||
	!materialStateSummary.materialValid ||
	materialStateSummary.materialRevision != 77 ||
	materialStateSummary.materialColor[0] < 0.24f ||
	materialStateSummary.materialColor[0] > 0.26f ||
	materialStateSummary.materialColor[2] < 0.74f ||
	materialStateSummary.materialColor[2] > 0.76f ||
	!source->getRealizedSceneMaterialSummary(3,
		materialStateSummary) ||
	materialStateSummary.nodeKind !=
	BRLObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
	materialStateSummary.materialValid ||
	source->getRealizedSceneMaterialSummary(4,
		materialStateSummary)) {
	printf("FAIL: database source should summarize realized shape material state\n");
	root->unref();
	return 1;
    }

    if (scene.getSceneMaterialSummaryCount() !=
	scene.getSceneTreeSummaryCount() ||
	!scene.getSceneMaterialSummary(0, materialStateSummary) ||
	!materialStateSummary.valid ||
	materialStateSummary.materialValid ||
	!scene.getSceneMaterialSummary(4, materialStateSummary) ||
	materialStateSummary.nodeKind !=
	BRLObolSceneTreeSummary::NODE_MESH_SHAPE ||
	materialStateSummary.ownerSourceIndex != 0 ||
	strcmp(materialStateSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	!materialStateSummary.materialValid ||
	materialStateSummary.materialRevision != 77 ||
	materialStateSummary.materialColor[1] < 0.49f ||
	materialStateSummary.materialColor[1] > 0.51f ||
	!scene.getSceneMaterialSummary(5, materialStateSummary) ||
	materialStateSummary.nodeKind !=
	BRLObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
	materialStateSummary.materialValid ||
	scene.getSceneMaterialSummary(6, materialStateSummary)) {
	printf("FAIL: scene controller should summarize Obol shape material state\n");
	root->unref();
	return 1;
    }

    BRLObolSceneBoundsSummary boundsSummary;
    if (source->getRealizedBoundsSummaryCount() !=
	source->getRealizedTreeSummaryCount() ||
	!source->getRealizedBoundsSummary(0, boundsSummary) ||
	!boundsSummary.valid ||
	boundsSummary.nodeKind !=
	BRLObolSceneTreeSummary::NODE_DATABASE_SOURCE ||
	!boundsSummary.boundsValid ||
	boundsSummary.ownerSourceIndex != -1 ||
	strcmp(boundsSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	boundsSummary.bounds.getMin()[0] > -1.74f ||
	boundsSummary.bounds.getMax()[0] < 1.74f ||
	boundsSummary.bounds.getMin()[1] > -1.74f ||
	boundsSummary.bounds.getMax()[1] < 1.74f ||
	!source->getRealizedBoundsSummary(1, boundsSummary) ||
	boundsSummary.nodeKind != BRLObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	!boundsSummary.boundsValid ||
	boundsSummary.bounds.getMin()[0] > -1.74f ||
	boundsSummary.bounds.getMax()[1] < 1.74f ||
	!source->getRealizedBoundsSummary(2, boundsSummary) ||
	boundsSummary.nodeKind != BRLObolSceneTreeSummary::NODE_MESH_SHAPE ||
	!boundsSummary.boundsValid ||
	boundsSummary.bounds.getMin()[0] < -0.01f ||
	boundsSummary.bounds.getMin()[0] > 0.01f ||
	boundsSummary.bounds.getMax()[0] < 0.99f ||
	boundsSummary.bounds.getMax()[1] < 0.99f ||
	!source->getRealizedBoundsSummary(3, boundsSummary) ||
	boundsSummary.nodeKind !=
	BRLObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
	boundsSummary.boundsValid ||
	source->getRealizedBoundsSummary(4, boundsSummary)) {
	printf("FAIL: database source should summarize realized subtree bounds\n");
	root->unref();
	return 1;
    }

    if (scene.getSceneBoundsSummaryCount() !=
	scene.getSceneTreeSummaryCount() ||
	!scene.getSceneBoundsSummary(0, boundsSummary) ||
	!boundsSummary.valid ||
	boundsSummary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
	!boundsSummary.boundsValid ||
	boundsSummary.ownerSourceIndex != -1 ||
	boundsSummary.bounds.getMin()[0] > -1.74f ||
	boundsSummary.bounds.getMax()[0] < 1.74f ||
	!scene.getSceneBoundsSummary(1, boundsSummary) ||
	boundsSummary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
	boundsSummary.boundsValid ||
	!scene.getSceneBoundsSummary(2, boundsSummary) ||
	boundsSummary.nodeKind !=
	BRLObolSceneTreeSummary::NODE_DATABASE_SOURCE ||
	boundsSummary.ownerSourceIndex != 0 ||
	strcmp(boundsSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	!boundsSummary.boundsValid ||
	!scene.getSceneBoundsSummary(4, boundsSummary) ||
	boundsSummary.nodeKind != BRLObolSceneTreeSummary::NODE_MESH_SHAPE ||
	boundsSummary.ownerSourceIndex != 0 ||
	boundsSummary.bounds.getMax()[0] < 0.99f ||
	boundsSummary.bounds.getMax()[1] < 0.99f ||
	!scene.getSceneBoundsSummary(5, boundsSummary) ||
	boundsSummary.nodeKind !=
	BRLObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
	boundsSummary.boundsValid ||
	scene.getSceneBoundsSummary(6, boundsSummary)) {
	printf("FAIL: scene controller should summarize Obol subtree bounds\n");
	root->unref();
	return 1;
    }

    const uint64_t beforeSourceDrawModeFrameRevision =
	scene.getFrameRevision();
    if (scene.setDatabaseSourceDrawMode("/summary/source",
					SoBRLDatabaseSource::SHADED) != 1 ||
	scene.getFrameRevision() <= beforeSourceDrawModeFrameRevision ||
	source->drawMode.getValue() != SoBRLDatabaseSource::SHADED ||
	summaryMesh->drawMode.getValue() != BRLOBOL_LOD_DRAW_SHADED) {
	printf("FAIL: scene controller should update source draw mode through the source owner\n");
	root->unref();
	return 1;
    }
    const uint64_t afterSourceShadedFrameRevision = scene.getFrameRevision();
    if (scene.setDatabaseSourceDrawMode("summary/source",
					SoBRLDatabaseSource::WIREFRAME) != 1 ||
	scene.getFrameRevision() <= afterSourceShadedFrameRevision ||
	source->drawMode.getValue() != SoBRLDatabaseSource::WIREFRAME ||
	summaryMesh->drawMode.getValue() != BRLOBOL_LOD_DRAW_WIRE ||
	!source->getRealizedShapeSummary(1, shapeSummary) ||
	shapeSummary.drawMode != BRLOBOL_LOD_DRAW_WIRE ||
	scene.setDatabaseSourceDrawMode("summary/source",
					SoBRLDatabaseSource::WIREFRAME) != 0 ||
	scene.setDatabaseSourceDrawMode("missing/source",
					SoBRLDatabaseSource::SHADED) != -1) {
	printf("FAIL: source draw-mode updates should sync database-intent realized shapes\n");
	root->unref();
	return 1;
    }

    const uint64_t beforeSourceMaterialPolicyFrameRevision =
	scene.getFrameRevision();
    if (scene.setDatabaseSourceMaterialPolicy("/summary/source",
	    SoBRLDatabaseSource::MATERIAL_DATABASE) != 1 ||
	scene.getFrameRevision() <= beforeSourceMaterialPolicyFrameRevision ||
	source->materialPolicy.getValue() !=
	SoBRLDatabaseSource::MATERIAL_DATABASE ||
	!source->getSummary(summary) ||
	summary.materialPolicy !=
	SoBRLDatabaseSource::MATERIAL_DATABASE) {
	printf("FAIL: scene controller should update source material policy through the source owner\n");
	root->unref();
	return 1;
    }
    const uint64_t afterSourceMaterialDatabaseFrameRevision =
	scene.getFrameRevision();
    if (scene.setDatabaseSourceMaterialPolicy("summary/source", 9999) != 1 ||
	scene.getFrameRevision() <=
	afterSourceMaterialDatabaseFrameRevision ||
	source->materialPolicy.getValue() !=
	SoBRLDatabaseSource::MATERIAL_INHERIT ||
	!source->getSummary(summary) ||
	summary.materialPolicy !=
	SoBRLDatabaseSource::MATERIAL_INHERIT ||
	scene.setDatabaseSourceMaterialPolicy("summary/source", -11) != 0 ||
	scene.setDatabaseSourceMaterialPolicy("missing/source",
		SoBRLDatabaseSource::MATERIAL_DATABASE) != -1) {
	printf("FAIL: source material policy updates should sanitize and report no-op/missing sources\n");
	root->unref();
	return 1;
    }

    root->unref();

    SoSeparator *groupRoot = new SoSeparator;
    groupRoot->ref();
    SoBRLSceneController groupScene(groupRoot);
    const uint64_t groupInitialStructuralRevision =
	groupScene.getStructuralRevision();
    const uint64_t groupInitialFrameRevision = groupScene.getFrameRevision();
    SoGroup *createdLeaf = groupScene.ensureGroup("assembly/leaf");
    const uint64_t afterGroupCreateStructuralRevision =
	groupScene.getStructuralRevision();
    const uint64_t afterGroupCreateFrameRevision =
	groupScene.getFrameRevision();
    SoBRLSceneGroup *retainedLeaf = NULL;
    if (createdLeaf &&
	createdLeaf->isOfType(SoBRLSceneGroup::getClassTypeId()))
	retainedLeaf = static_cast<SoBRLSceneGroup *>(createdLeaf);
    if (groupScene.findGroup("/") != groupRoot ||
	groupScene.findGroup(NULL) != NULL ||
	!createdLeaf ||
	!retainedLeaf ||
	strcmp(retainedLeaf->groupPath.getValue().getString(),
	       "assembly/leaf") != 0 ||
	groupScene.findGroup("assembly/leaf") != createdLeaf ||
	groupScene.ensureGroup("/assembly/leaf") != createdLeaf ||
	groupScene.getStructuralRevision() !=
	afterGroupCreateStructuralRevision ||
	groupScene.getFrameRevision() != afterGroupCreateFrameRevision ||
	groupScene.getGroupChildCount("/") != 1 ||
	groupScene.getGroupChildCount("assembly") != 1 ||
	groupScene.getGroupChildCount("missing") != -1 ||
	groupScene.getGroupDescendantGroupCount("/") != 2 ||
	groupScene.getGroupDescendantGroupCount("assembly") != 1 ||
	groupScene.getGroupDescendantGroupCount("assembly/leaf") != 0 ||
	groupScene.getGroupDescendantGroupCount("missing") != -1 ||
	groupScene.getStructuralRevision() <=
	groupInitialStructuralRevision ||
	groupScene.getFrameRevision() <= groupInitialFrameRevision) {
	printf("FAIL: scene controller should own slash-path group creation and lookup\n");
	groupRoot->unref();
	return 1;
    }

    if (groupScene.getSceneTreeSummaryCount() != 3 ||
	!groupScene.getSceneTreeSummary(0, treeSummary) ||
	treeSummary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
	strcmp(treeSummary.path.getString(), "/") != 0 ||
	treeSummary.childCount != 1 ||
	treeSummary.drawTreeDepth != 0 ||
	!groupScene.getSceneTreeSummary(1, treeSummary) ||
	treeSummary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
	strcmp(treeSummary.path.getString(), "assembly") != 0 ||
	treeSummary.childCount != 1 ||
	treeSummary.drawTreeDepth != 1 ||
	!groupScene.getSceneTreeSummary(2, treeSummary) ||
	treeSummary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
	strcmp(treeSummary.path.getString(), "assembly/leaf") != 0 ||
	treeSummary.childCount != 0 ||
	treeSummary.drawTreeDepth != 2 ||
	groupScene.getSceneTreeSummary(3, treeSummary)) {
	printf("FAIL: scene controller should summarize retained named group paths\n");
	groupRoot->unref();
	return 1;
    }

    if (!groupScene.getSceneDisplaySummary(2, displaySummary) ||
	displaySummary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
	strcmp(displaySummary.path.getString(), "assembly/leaf") != 0 ||
	!groupScene.getSceneBoundsSummary(2, boundsSummary) ||
	boundsSummary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
	strcmp(boundsSummary.path.getString(), "assembly/leaf") != 0) {
	printf("FAIL: scene controller display and bounds summaries should retain group paths\n");
	groupRoot->unref();
	return 1;
    }

    const uint64_t beforeGroupStateStructuralRevision =
	groupScene.getStructuralRevision();
    const uint64_t beforeGroupStateFrameRevision =
	groupScene.getFrameRevision();
    if (groupScene.setGroupDrawIntent("assembly/leaf",
				      "draw://assembly/leaf", BRLOBOL_LOD_DRAW_SHADED,
				      BRLOBOL_LOD_DRAW_WIRE, TRUE, 42) != 1 ||
	groupScene.setGroupDisplayState("assembly/leaf", FALSE, TRUE,
					TRUE, 6, 7, 0.45f, TRUE,
					SbColor(0.8f, 0.2f, 0.1f), TRUE,
					SbColor(0.3f, 0.4f, 0.5f), 123) != 1 ||
	groupScene.getStructuralRevision() !=
	beforeGroupStateStructuralRevision ||
	groupScene.getFrameRevision() <= beforeGroupStateFrameRevision ||
	!retainedLeaf->selected.getValue() ||
	!retainedLeaf->overlayIntent.getValue() ||
	retainedLeaf->fallbackDrawMode.getValue() !=
	BRLOBOL_LOD_DRAW_WIRE ||
	retainedLeaf->revalidationRevision.getValue() != 42) {
	printf("FAIL: scene controller should retain group draw and display state in Obol\n");
	groupRoot->unref();
	return 1;
    }
    const uint64_t afterGroupStateFrameRevision =
	groupScene.getFrameRevision();
    if (groupScene.setGroupDrawIntent("assembly/leaf",
				      "draw://assembly/leaf", BRLOBOL_LOD_DRAW_SHADED,
				      BRLOBOL_LOD_DRAW_WIRE, TRUE, 42) != 0 ||
	groupScene.setGroupDisplayState("assembly/leaf", FALSE, TRUE,
					TRUE, 6, 7, 0.45f, TRUE,
					SbColor(0.8f, 0.2f, 0.1f), TRUE,
					SbColor(0.3f, 0.4f, 0.5f), 123) != 0 ||
	groupScene.getFrameRevision() != afterGroupStateFrameRevision ||
	groupScene.setGroupDrawIntent("missing", "missing",
				      BRLOBOL_LOD_DRAW_WIRE, BRLOBOL_LOD_DRAW_WIRE,
				      FALSE, 0) != -1) {
	printf("FAIL: scene controller should treat retained group state reapplication as a no-op\n");
	groupRoot->unref();
	return 1;
    }

    if (!groupScene.getSceneDisplaySummary(2, displaySummary) ||
	displaySummary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
	strcmp(displaySummary.path.getString(), "assembly/leaf") != 0 ||
	!displaySummary.hasDrawIntent ||
	strcmp(displaySummary.intentPath.getString(),
	       "draw://assembly/leaf") != 0 ||
	displaySummary.intentDrawMode != BRLOBOL_LOD_DRAW_SHADED ||
	displaySummary.visible ||
	!displaySummary.highlighted ||
	displaySummary.lineStyle != 6 ||
	displaySummary.lineWidth != 7 ||
	fabs(displaySummary.transparency - 0.45) > 0.001 ||
	displaySummary.drawMode != BRLOBOL_LOD_DRAW_SHADED ||
	!displaySummary.materialValid ||
	displaySummary.materialRevision != 123 ||
	displaySummary.materialColor[0] < 0.29f ||
	displaySummary.materialColor[1] < 0.39f ||
	displaySummary.materialColor[2] < 0.49f ||
	!groupScene.getSceneMaterialSummary(2, materialStateSummary) ||
	materialStateSummary.materialValid) {
	printf("FAIL: scene controller summaries should expose retained group display state without making groups shape materials\n");
	groupRoot->unref();
	return 1;
    }

    SoBRLMeshShape *groupMesh = make_mesh("assembly/leaf/mesh", "mesh");
    const uint64_t beforeGroupAppendStructuralRevision =
	groupScene.getStructuralRevision();
    if (groupScene.appendChildToGroup("assembly/leaf", groupMesh) != 1 ||
	groupScene.appendChildToGroup("assembly/leaf", groupMesh) != 0 ||
	groupScene.getGroupChildCount("assembly/leaf") != 1 ||
	groupScene.getStructuralRevision() <=
	beforeGroupAppendStructuralRevision ||
	groupScene.appendChildToGroup("missing", groupMesh) != -1 ||
	groupScene.removeChildFromGroup("missing", groupMesh) != -1) {
	printf("FAIL: scene controller should append retained children to named groups\n");
	groupRoot->unref();
	return 1;
    }
    const uint64_t afterGroupAppendStructuralRevision =
	groupScene.getStructuralRevision();
    if (groupScene.getSceneTreeSummaryCount() != 4 ||
	!groupScene.getSceneTreeSummary(3, treeSummary) ||
	treeSummary.nodeKind != BRLObolSceneTreeSummary::NODE_MESH_SHAPE ||
	strcmp(treeSummary.path.getString(), "assembly/leaf/mesh") != 0 ||
	treeSummary.drawTreeDepth != 3 ||
	!groupScene.getSceneBoundsSummary(2, boundsSummary) ||
	!boundsSummary.boundsValid ||
	boundsSummary.bounds.getMax()[0] < 0.99f ||
	boundsSummary.bounds.getMax()[1] < 0.99f ||
	!groupScene.getSceneBoundsSummary(3, boundsSummary) ||
	boundsSummary.nodeKind != BRLObolSceneTreeSummary::NODE_MESH_SHAPE ||
	!boundsSummary.boundsValid) {
	printf("FAIL: scene controller retained group children should affect summaries and bounds\n");
	groupRoot->unref();
	return 1;
    }
    const uint64_t beforeShapeStateStructuralRevision =
	groupScene.getStructuralRevision();
    const uint64_t beforeShapeStateFrameRevision =
	groupScene.getFrameRevision();
    if (groupScene.findShape("assembly/leaf/mesh") != groupMesh ||
	groupScene.findShape("/assembly/leaf/mesh") != groupMesh ||
	groupScene.findShapeParent("assembly/leaf/mesh") != createdLeaf ||
	groupScene.findShape("missing") != NULL ||
	groupScene.setShapeDrawState("assembly/leaf/mesh",
				     BRLOBOL_LOD_DRAW_POINTS, TRUE, FALSE, FALSE) != 1 ||
	groupScene.setShapeDisplayState("assembly/leaf/mesh",
					FALSE, TRUE, TRUE, 8, 9, 0.55f, TRUE,
					SbColor(0.2f, 0.3f, 0.4f), TRUE,
					SbColor(0.5f, 0.6f, 0.7f), 222) != 1 ||
	groupScene.getStructuralRevision() !=
	beforeShapeStateStructuralRevision ||
	groupScene.getFrameRevision() <= beforeShapeStateFrameRevision) {
	printf("FAIL: scene controller should find and update retained shape state by path\n");
	groupRoot->unref();
	return 1;
    }
    const uint64_t afterShapeStateFrameRevision =
	groupScene.getFrameRevision();
    if (groupScene.setShapeDrawState("assembly/leaf/mesh",
				     BRLOBOL_LOD_DRAW_POINTS, TRUE, FALSE, FALSE) != 0 ||
	groupScene.setShapeDisplayState("assembly/leaf/mesh",
					FALSE, TRUE, TRUE, 8, 9, 0.55f, TRUE,
					SbColor(0.2f, 0.3f, 0.4f), TRUE,
					SbColor(0.5f, 0.6f, 0.7f), 222) != 0 ||
	groupScene.getFrameRevision() != afterShapeStateFrameRevision ||
	groupScene.setShapeDisplayState("missing", TRUE, FALSE,
					FALSE, 0, 0, 0.0f, FALSE,
					SbColor(1.0f, 1.0f, 1.0f), FALSE,
					SbColor(1.0f, 1.0f, 1.0f), 0) != -1) {
	printf("FAIL: scene controller should treat retained shape state reapplication as a no-op\n");
	groupRoot->unref();
	return 1;
    }
    if (!groupScene.getSceneDisplaySummary(3, displaySummary) ||
	displaySummary.nodeKind != BRLObolSceneTreeSummary::NODE_MESH_SHAPE ||
	strcmp(displaySummary.path.getString(), "assembly/leaf/mesh") != 0 ||
	displaySummary.visible ||
	!displaySummary.highlighted ||
	displaySummary.lineStyle != 8 ||
	displaySummary.lineWidth != 9 ||
	fabs(displaySummary.transparency - 0.55) > 0.001 ||
	displaySummary.drawMode != BRLOBOL_LOD_DRAW_POINTS ||
	!displaySummary.materialValid ||
	displaySummary.materialRevision != 222 ||
	displaySummary.materialColor[0] < 0.49f ||
	displaySummary.materialColor[1] < 0.59f ||
	displaySummary.materialColor[2] < 0.69f ||
	!groupScene.getSceneMaterialSummary(3, materialStateSummary) ||
	!materialStateSummary.materialValid ||
	materialStateSummary.materialRevision != 222) {
	printf("FAIL: scene controller summaries should expose retained shape state\n");
	groupRoot->unref();
	return 1;
    }
    const uint64_t beforeShapeSourceFrameRevision =
	groupScene.getFrameRevision();
    SbMatrix shapePlacementMatrix;
    shapePlacementMatrix.setTranslate(SbVec3f(2.0f, 3.0f, 4.0f));
    if (groupScene.setShapeSourceState("assembly/leaf/mesh",
				       "/owner/source", 10, 11, 12, 13, 14, 15, 16, 17,
				       "stale source", "realized::shape", TRUE, 18) != 1 ||
	groupScene.getFrameRevision() <= beforeShapeSourceFrameRevision ||
	strcmp(groupMesh->ownerSourcePath.getValue().getString(),
	       "/owner/source") != 0 ||
	groupMesh->ownerSourceRevision.getValue() != 10 ||
	groupMesh->ownerInputsRevision.getValue() != 11 ||
	groupMesh->ownerViewRevision.getValue() != 12 ||
	groupMesh->ownerRealizedRevision.getValue() != 13 ||
	groupMesh->ownerRealizedSourceRevision.getValue() != 14 ||
	groupMesh->ownerRealizedInputsRevision.getValue() != 15 ||
	groupMesh->ownerRealizedViewRevision.getValue() != 16 ||
	groupMesh->ownerRealizationStatus.getValue() != 17 ||
	strcmp(groupMesh->ownerRealizationDiagnostic.getValue().getString(),
	       "stale source") != 0 ||
	strcmp(groupMesh->ownerRealizationIdentity.getValue().getString(),
	       "realized::shape") != 0 ||
	!groupMesh->ownerSourceStale.getValue() ||
	groupMesh->ownerStaleReason.getValue() != 18) {
	printf("FAIL: scene controller should retain shape source revision state\n");
	groupRoot->unref();
	return 1;
    }
    const uint64_t beforeShapePlacementFrameRevision =
	groupScene.getFrameRevision();
    if (groupScene.setShapePlacementState("assembly/leaf/mesh",
					  TRUE, shapePlacementMatrix, TRUE,
					  SbVec3f(5.0f, 6.0f, 7.0f), TRUE, 12.5f) != 1 ||
	groupScene.getFrameRevision() <=
	beforeShapePlacementFrameRevision ||
	!groupMesh->drawMatrixValid.getValue() ||
	!matrix_nearly_equal(groupMesh->drawMatrix.getValue(),
			     shapePlacementMatrix) ||
	!groupMesh->drawCenterValid.getValue() ||
	groupMesh->drawCenter.getValue()[0] < 4.99f ||
	groupMesh->drawCenter.getValue()[1] < 5.99f ||
	groupMesh->drawCenter.getValue()[2] < 6.99f ||
	!groupMesh->drawSizeValid.getValue() ||
	fabs(groupMesh->drawSize.getValue() - 12.5f) > 0.001f) {
	printf("FAIL: scene controller should retain shape placement state\n");
	groupRoot->unref();
	return 1;
    }
    const uint64_t afterShapePlacementFrameRevision =
	groupScene.getFrameRevision();
    if (groupScene.setShapePlacementState("assembly/leaf/mesh",
					  TRUE, shapePlacementMatrix, TRUE,
					  SbVec3f(5.0f, 6.0f, 7.0f), TRUE, 12.5f) != 0 ||
	groupScene.getFrameRevision() !=
	afterShapePlacementFrameRevision ||
	groupScene.setShapePlacementState("missing", TRUE,
					  shapePlacementMatrix, TRUE, SbVec3f(5.0f, 6.0f, 7.0f),
					  TRUE, 12.5f) != -1 ||
	!groupScene.getSceneDisplaySummary(3, displaySummary) ||
	!displaySummary.drawMatrixValid ||
	!matrix_nearly_equal(displaySummary.drawMatrix,
			     shapePlacementMatrix) ||
	!displaySummary.drawCenterValid ||
	displaySummary.drawCenter[0] < 4.99f ||
	displaySummary.drawCenter[1] < 5.99f ||
	displaySummary.drawCenter[2] < 6.99f ||
	!displaySummary.drawSizeValid ||
	fabs(displaySummary.drawSize - 12.5f) > 0.001f) {
	printf("FAIL: scene controller should summarize retained shape placement state\n");
	groupRoot->unref();
	return 1;
    }
    if (groupScene.setShapeSourceState("assembly/leaf/mesh",
				       "/owner/source", 10, 11, 12, 13, 14, 15, 16, 17,
				       "stale source", "realized::shape", TRUE, 18) != 0 ||
	groupScene.getFrameRevision() != afterShapePlacementFrameRevision ||
	groupScene.setShapeSourceState("missing", "/owner/source",
				       10, 11, 12, 13, 14, 15, 16, 17, "stale source",
				       "realized::shape", TRUE, 18) != -1 ||
	!groupScene.getSceneTreeSummary(3, treeSummary) ||
	strcmp(treeSummary.ownerSourcePath.getString(),
	       "/owner/source") != 0 ||
	!groupScene.getSceneDisplaySummary(3, displaySummary) ||
	strcmp(displaySummary.ownerSourcePath.getString(),
	       "/owner/source") != 0 ||
	!groupScene.getSceneBoundsSummary(3, boundsSummary) ||
	strcmp(boundsSummary.ownerSourcePath.getString(),
	       "/owner/source") != 0) {
	printf("FAIL: scene controller should summarize retained shape owner source state\n");
	groupRoot->unref();
	return 1;
    }
    SoGroup *assemblyGroup = groupScene.findGroup("assembly");
    const uint64_t beforeShapeMoveStructuralRevision =
	groupScene.getStructuralRevision();
    if (!assemblyGroup ||
	groupScene.moveShapeToGroup("assembly/leaf/mesh",
				    "assembly/leaf") != 0 ||
	groupScene.moveShapeToGroup("assembly/leaf/mesh", "assembly") != 1 ||
	groupScene.findShapeParent("assembly/leaf/mesh") !=
	assemblyGroup ||
	groupScene.getGroupChildCount("assembly/leaf") != 0 ||
	groupScene.getGroupChildCount("assembly") != 2 ||
	groupScene.getStructuralRevision() <=
	beforeShapeMoveStructuralRevision ||
	groupScene.moveShapeToGroup("assembly/leaf/mesh",
				    "assembly/leaf") != 1 ||
	groupScene.findShapeParent("assembly/leaf/mesh") != createdLeaf ||
	groupScene.getGroupChildCount("assembly/leaf") != 1 ||
	groupScene.getGroupChildCount("assembly") != 1 ||
	groupScene.moveShapeToGroup("missing", "assembly") != 0 ||
	groupScene.moveShapeToGroup("assembly/leaf/mesh",
				    "missing") != -1) {
	printf("FAIL: scene controller should move retained shape ownership between groups\n");
	groupRoot->unref();
	return 1;
    }
    if (groupScene.removeChildFromGroup("assembly/leaf", groupMesh) != 1 ||
	groupScene.getGroupChildCount("assembly/leaf") != 0 ||
	groupScene.getStructuralRevision() <=
	afterGroupAppendStructuralRevision ||
	groupScene.getSceneTreeSummaryCount() != 3) {
	printf("FAIL: scene controller should remove retained children from named groups\n");
	groupRoot->unref();
	return 1;
    }
    SoBRLMeshShape *removeMesh = make_mesh("assembly/leaf/remove", "remove");
    const uint64_t beforeShapeRemoveStructuralRevision =
	groupScene.getStructuralRevision();
    if (groupScene.appendChildToGroup("assembly/leaf", removeMesh) != 1 ||
	groupScene.removeShape("/assembly/leaf/remove") != 1 ||
	groupScene.findShape("assembly/leaf/remove") != NULL ||
	groupScene.getGroupChildCount("assembly/leaf") != 0 ||
	groupScene.getStructuralRevision() <=
	beforeShapeRemoveStructuralRevision ||
	groupScene.removeShape("assembly/leaf/remove") != 0) {
	printf("FAIL: scene controller should remove retained shapes by path\n");
	groupRoot->unref();
	return 1;
    }

    SoGroup *createdChild = groupScene.ensureGroup("assembly/leaf/child");
    SoBRLSceneGroup *retainedChild = NULL;
    if (createdChild &&
	createdChild->isOfType(SoBRLSceneGroup::getClassTypeId()))
	retainedChild = static_cast<SoBRLSceneGroup *>(createdChild);
    const uint64_t afterGroupChildStructuralRevision =
	groupScene.getStructuralRevision();
    if (!createdChild ||
	!retainedChild ||
	groupScene.getGroupChildCount("assembly/leaf") != 1 ||
	groupScene.renameGroup("assembly/leaf", "renamed") != 1 ||
	groupScene.findGroup("assembly/leaf") != NULL ||
	groupScene.findGroup("assembly/renamed") != createdLeaf ||
	groupScene.findGroup("assembly/renamed/child") != createdChild ||
	strcmp(retainedLeaf->groupPath.getValue().getString(),
	       "assembly/renamed") != 0 ||
	strcmp(retainedChild->groupPath.getValue().getString(),
	       "assembly/renamed/child") != 0 ||
	groupScene.getStructuralRevision() <=
	afterGroupChildStructuralRevision ||
	!groupScene.getSceneTreeSummary(2, treeSummary) ||
	strcmp(treeSummary.path.getString(), "assembly/renamed") != 0 ||
	groupScene.renameGroup("assembly/renamed", "renamed") != 0 ||
	groupScene.renameGroup("assembly/renamed", "bad/name") != 0) {
	printf("FAIL: scene controller should rename retained group leaf paths in place\n");
	groupRoot->unref();
	return 1;
    }

    const uint64_t afterGroupRenameStructuralRevision =
	groupScene.getStructuralRevision();
    if (groupScene.eraseGroupSubpath("assembly", "renamed/child") != 1 ||
	groupScene.findGroup("assembly/renamed/child") != NULL ||
	groupScene.getGroupChildCount("assembly/renamed") != 0 ||
	groupScene.getStructuralRevision() <=
	afterGroupRenameStructuralRevision ||
	groupScene.eraseGroupSubpath("assembly", "renamed/child") != 0 ||
	groupScene.eraseGroupSubpath("missing", "renamed") != -1) {
	printf("FAIL: scene controller should erase nested retained group subpaths\n");
	groupRoot->unref();
	return 1;
    }
    createdChild = groupScene.ensureGroup("assembly/renamed/child");
    if (!createdChild) {
	printf("FAIL: scene controller should recreate nested retained group subpaths\n");
	groupRoot->unref();
	return 1;
    }
    if (groupScene.clearGroup("assembly/renamed") != 1 ||
	groupScene.getGroupChildCount("assembly/renamed") != 0 ||
	groupScene.clearGroup("assembly/renamed") != 0 ||
	groupScene.removeGroup("assembly/renamed") != 1 ||
	groupScene.findGroup("assembly/renamed") != NULL ||
	groupScene.getGroupChildCount("assembly") != 0 ||
	groupScene.removeGroup("assembly/renamed") != 0 ||
	groupScene.removeGroup("assembly") != 1 ||
	groupScene.getGroupChildCount("/") != 0 ||
	groupScene.clearGroup("/") != 0) {
	printf("FAIL: scene controller should clear and remove retained groups with revision accounting\n");
	groupRoot->unref();
	return 1;
    }
    {
	BRLObolViewController groupView(groupRoot);
	SbString renderReason;
	groupView.clearRenderRequest();
	SoGroup *viewLeaf = groupView.ensureGroup("view/leaf");
	if (!viewLeaf ||
	    groupView.findGroup("view/leaf") != viewLeaf ||
	    groupView.getGroupChildCount("view") != 1 ||
	    !groupView.consumeRenderRequest(&renderReason) ||
	    strcmp(renderReason.getString(), "scene-group") != 0) {
	    printf("FAIL: view controller should expose retained group creation and render invalidation\n");
	    groupRoot->unref();
	    return 1;
	}
	groupView.clearRenderRequest();
	if (groupView.ensureGroup("view/leaf") != viewLeaf ||
	    groupView.consumeRenderRequest(&renderReason)) {
	    printf("FAIL: view controller should not request render for no-op group ensure\n");
	    groupRoot->unref();
	    return 1;
	}
	if (viewLeaf->isOfType(SoBRLSceneGroup::getClassTypeId()) != TRUE ||
	    groupView.setGroupDrawIntent("view/leaf",
					 "draw://view/leaf", BRLOBOL_LOD_DRAW_POINTS,
					 BRLOBOL_LOD_DRAW_WIRE, FALSE, 5) != 1 ||
	    !groupView.consumeRenderRequest(&renderReason) ||
	    strcmp(renderReason.getString(), "scene-group") != 0) {
	    printf("FAIL: view controller should expose retained group draw intent state\n");
	    groupRoot->unref();
	    return 1;
	}
	groupView.clearRenderRequest();
	if (groupView.setGroupDrawIntent("view/leaf",
					 "draw://view/leaf", BRLOBOL_LOD_DRAW_POINTS,
					 BRLOBOL_LOD_DRAW_WIRE, FALSE, 5) != 0 ||
	    groupView.consumeRenderRequest(&renderReason)) {
	    printf("FAIL: view controller should not render for no-op retained group draw intent state\n");
	    groupRoot->unref();
	    return 1;
	}
	if (groupView.setGroupDisplayState("view/leaf", TRUE, TRUE,
					   FALSE, 1, 2, 0.25f, TRUE,
					   SbColor(0.1f, 0.2f, 0.3f), FALSE,
					   SbColor(1.0f, 1.0f, 1.0f), 77) != 1 ||
	    !groupView.consumeRenderRequest(&renderReason) ||
	    strcmp(renderReason.getString(), "scene-group") != 0 ||
	    groupView.setGroupDisplayState("missing", TRUE, FALSE,
					   FALSE, 0, 0, 0.0f, FALSE,
					   SbColor(1.0f, 1.0f, 1.0f), FALSE,
					   SbColor(1.0f, 1.0f, 1.0f), 0) != -1) {
	    printf("FAIL: view controller should expose retained group display state\n");
	    groupRoot->unref();
	    return 1;
	}
	groupView.clearRenderRequest();

	SoBRLMeshShape *viewMesh = make_mesh("view/leaf/mesh", "mesh");
	if (groupView.appendChildToGroup("view/leaf", viewMesh) != 1 ||
	    groupView.appendChildToGroup("view/leaf", viewMesh) != 0 ||
	    groupView.getGroupChildCount("view/leaf") != 1 ||
	    !groupView.consumeRenderRequest(&renderReason) ||
	    strcmp(renderReason.getString(), "scene-group") != 0) {
	    printf("FAIL: view controller should expose retained group child ownership\n");
	    groupRoot->unref();
	    return 1;
	}
	groupView.clearRenderRequest();
	if (groupView.findShape("view/leaf/mesh") != viewMesh ||
	    groupView.findShapeParent("view/leaf/mesh") != viewLeaf ||
	    groupView.setShapeDrawState("view/leaf/mesh",
					BRLOBOL_LOD_DRAW_DIAGNOSTIC, FALSE, TRUE, FALSE) != 1 ||
	    !groupView.consumeRenderRequest(&renderReason) ||
	    strcmp(renderReason.getString(), "scene-shape") != 0 ||
	    groupView.setShapeDrawState("view/leaf/mesh",
					BRLOBOL_LOD_DRAW_DIAGNOSTIC, FALSE, TRUE, FALSE) != 0) {
	    printf("FAIL: view controller should expose retained shape draw state\n");
	    groupRoot->unref();
	    return 1;
	}
	groupView.clearRenderRequest();
	if (groupView.setShapeDisplayState("view/leaf/mesh", TRUE, TRUE,
					   TRUE, 3, 4, 0.15f, TRUE,
					   SbColor(0.3f, 0.2f, 0.1f), TRUE,
					   SbColor(0.6f, 0.5f, 0.4f), 88) != 1 ||
	    !groupView.consumeRenderRequest(&renderReason) ||
	    strcmp(renderReason.getString(), "scene-shape") != 0 ||
	    groupView.setShapeDisplayState("missing", TRUE, FALSE,
					   FALSE, 0, 0, 0.0f, FALSE,
					   SbColor(1.0f, 1.0f, 1.0f), FALSE,
					   SbColor(1.0f, 1.0f, 1.0f), 0) != -1) {
	    printf("FAIL: view controller should expose retained shape display state\n");
	    groupRoot->unref();
	    return 1;
	}
	groupView.clearRenderRequest();
	SbMatrix viewPlacementMatrix;
	viewPlacementMatrix.setTranslate(SbVec3f(8.0f, 9.0f, 10.0f));
	if (groupView.setShapePlacementState("view/leaf/mesh",
					     TRUE, viewPlacementMatrix, TRUE,
					     SbVec3f(11.0f, 12.0f, 13.0f), TRUE, 14.5f) != 1 ||
	    !groupView.consumeRenderRequest(&renderReason) ||
	    strcmp(renderReason.getString(), "scene-shape") != 0 ||
	    groupView.setShapePlacementState("view/leaf/mesh",
					     TRUE, viewPlacementMatrix, TRUE,
					     SbVec3f(11.0f, 12.0f, 13.0f), TRUE, 14.5f) != 0 ||
	    groupView.setShapePlacementState("missing", TRUE,
					     viewPlacementMatrix, TRUE,
					     SbVec3f(11.0f, 12.0f, 13.0f), TRUE, 14.5f) != -1) {
	    printf("FAIL: view controller should expose retained shape placement state\n");
	    groupRoot->unref();
	    return 1;
	}
	groupView.clearRenderRequest();
	if (groupView.setShapeSourceState("view/leaf/mesh",
					  "view/source", 21, 22, 23, 24, 25, 26, 27, 28,
					  "view diagnostic", "view identity", TRUE, 29) != 1 ||
	    !groupView.consumeRenderRequest(&renderReason) ||
	    strcmp(renderReason.getString(), "scene-shape") != 0 ||
	    groupView.setShapeSourceState("view/leaf/mesh",
					  "view/source", 21, 22, 23, 24, 25, 26, 27, 28,
					  "view diagnostic", "view identity", TRUE, 29) != 0 ||
	    groupView.setShapeSourceState("missing", "view/source",
					  21, 22, 23, 24, 25, 26, 27, 28,
					  "view diagnostic", "view identity", TRUE, 29) != -1) {
	    printf("FAIL: view controller should expose retained shape source state\n");
	    groupRoot->unref();
	    return 1;
	}
	SoGroup *viewRootGroup = groupView.findGroup("view");
	groupView.clearRenderRequest();
	if (!viewRootGroup ||
	    groupView.moveShapeToGroup("view/leaf/mesh", "view") != 1 ||
	    groupView.findShapeParent("view/leaf/mesh") != viewRootGroup ||
	    !groupView.consumeRenderRequest(&renderReason) ||
	    strcmp(renderReason.getString(), "scene-shape") != 0 ||
	    groupView.moveShapeToGroup("view/leaf/mesh", "view/leaf") != 1 ||
	    groupView.findShapeParent("view/leaf/mesh") != viewLeaf ||
	    groupView.removeShape("view/leaf/mesh") != 1 ||
	    groupView.findShape("view/leaf/mesh") != NULL ||
	    groupView.getGroupChildCount("view/leaf") != 0) {
	    printf("FAIL: view controller should expose retained shape ownership changes\n");
	    groupRoot->unref();
	    return 1;
	}
	groupView.clearRenderRequest();
	if (groupView.ensureGroup("view/leaf/child") == NULL ||
	    groupView.renameGroup("view/leaf", "renamed") != 1 ||
	    groupView.findGroup("view/renamed/child") == NULL ||
	    groupView.eraseGroupSubpath("view", "renamed/child") != 1 ||
	    groupView.findGroup("view/renamed/child") != NULL ||
	    groupView.clearGroup("view/renamed") != 0 ||
	    groupView.removeGroup("view/renamed") != 1 ||
	    groupView.removeGroup("view") != 1 ||
	    groupView.getGroupChildCount("/") != 0 ||
	    !groupView.consumeRenderRequest(&renderReason) ||
	    strcmp(renderReason.getString(), "scene-group") != 0) {
	    printf("FAIL: view controller should expose retained group rename, nested erase, clear, and remove\n");
	    groupRoot->unref();
	    return 1;
	}
    }
    groupRoot->unref();

    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;
    if (make_submit_test_db(dbpath, sizeof(dbpath), &dbip))
	return 1;

    SoSeparator *ownedRoot = new SoSeparator;
    ownedRoot->ref();
    SoSeparator *ownedNonSourceChild = new SoSeparator;
    ownedRoot->addChild(ownedNonSourceChild);
    SoBRLSceneController ownedScene(ownedRoot);
    const uint64_t ownedInitialStructuralRevision =
	ownedScene.getStructuralRevision();
    const uint64_t ownedInitialFrameRevision = ownedScene.getFrameRevision();

    if (ownedScene.replaceDatabaseSource("lod-submit.bot", dbip,
					 SoBRLDatabaseSource::WIREFRAME, 21) != 1 ||
	ownedScene.getDatabaseSourceCount() != 1 ||
	ownedRoot->getNumChildren() != 2 ||
	ownedScene.getStructuralRevision() <= ownedInitialStructuralRevision ||
	ownedScene.getFrameRevision() <= ownedInitialFrameRevision ||
	!ownedScene.getSceneSummary(sceneSummary) ||
	sceneSummary.rootChildCount != 2 ||
	sceneSummary.databaseSourceCount != 1 ||
	sceneSummary.nonDatabaseRootChildCount != 1 ||
	sceneSummary.structuralRevision != ownedScene.getStructuralRevision() ||
	sceneSummary.frameRevision != ownedScene.getFrameRevision() ||
	!ownedScene.findDatabaseSource("/lod-submit.bot")) {
	printf("FAIL: scene controller should own database source replacement\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    const uint64_t afterAddStructuralRevision =
	ownedScene.getStructuralRevision();
    const uint64_t afterAddFrameRevision = ownedScene.getFrameRevision();

    SoBRLDatabaseSource *ownedSource =
	ownedScene.findDatabaseSource("lod-submit.bot");
    if (!ownedSource ||
	ownedSource->getDatabase() != dbip ||
	ownedSource->drawMode.getValue() != SoBRLDatabaseSource::WIREFRAME ||
	ownedSource->sourceRevision.getValue() != 21 ||
	!ownedScene.getDatabaseSourceSummary(0, summary) ||
	strcmp(summary.path.getString(), "lod-submit.bot") != 0 ||
	!summary.hasParent ||
	summary.drawTreeDepth != 1 ||
	strcmp(summary.parentGroupPath.getString(), "/") != 0 ||
	summary.drawMode != SoBRLDatabaseSource::WIREFRAME ||
	summary.sourceRevision != 21) {
	printf("FAIL: scene controller replacement should preserve source configuration\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    if (ownedScene.replaceDatabaseSource("/lod-submit.bot", dbip,
					 SoBRLDatabaseSource::SHADED, 22) != 1 ||
	ownedScene.getDatabaseSourceCount() != 1 ||
	ownedRoot->getNumChildren() != 2 ||
	ownedScene.getStructuralRevision() != afterAddStructuralRevision ||
	ownedScene.getFrameRevision() <= afterAddFrameRevision ||
	!ownedScene.getSceneSummary(sceneSummary) ||
	sceneSummary.rootChildCount != 2 ||
	sceneSummary.databaseSourceCount != 1 ||
	sceneSummary.nonDatabaseRootChildCount != 1 ||
	sceneSummary.structuralRevision != ownedScene.getStructuralRevision() ||
	sceneSummary.frameRevision != ownedScene.getFrameRevision() ||
	ownedScene.findDatabaseSource("lod-submit.bot") != ownedSource ||
	ownedSource->drawMode.getValue() != SoBRLDatabaseSource::SHADED ||
	ownedSource->sourceRevision.getValue() != 22) {
	printf("FAIL: scene controller should replace slash-equivalent database sources in place\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    if (ownedScene.moveDatabaseSourceToGroup("lod-submit.bot",
	    "draw/group") != 1 ||
	ownedScene.getDatabaseSourceCount() != 1 ||
	ownedScene.getDatabaseSource(0) != ownedSource ||
	ownedScene.findDatabaseSource("/lod-submit.bot") != ownedSource ||
	ownedRoot->getNumChildren() != 2 ||
	ownedScene.getGroupChildCount("draw/group") != 1 ||
	!ownedScene.getDatabaseSourceSummary(0, summary) ||
	!summary.hasParent ||
	summary.drawTreeDepth != 3 ||
	strcmp(summary.parentGroupPath.getString(), "draw/group") != 0 ||
	!ownedScene.getSceneSummary(sceneSummary) ||
	sceneSummary.databaseSourceCount != 1 ||
	sceneSummary.nonDatabaseRootChildCount != 2) {
	printf("FAIL: scene controller source summaries should expose retained group ownership\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    const uint64_t afterMoveStructuralRevision =
	ownedScene.getStructuralRevision();
    const uint64_t afterMoveFrameRevision = ownedScene.getFrameRevision();
    const SbColor sourceColor(0.5f, 0.25f, 0.125f);
    const SbColor sourceMaterial(0.25f, 0.5f, 0.75f);
    if (ownedScene.setDatabaseSourceState("lod-submit.bot", TRUE, 24, 44,
					  FALSE, FALSE, TRUE, 2, 6, 0.35f, TRUE, sourceColor, TRUE,
					  sourceMaterial, 77) != 1 ||
	ownedScene.getStructuralRevision() != afterMoveStructuralRevision ||
	ownedScene.getFrameRevision() <= afterMoveFrameRevision ||
	!ownedScene.getDatabaseSourceSummary(0, summary) ||
	summary.sourceRevision != 24 ||
	summary.inputsRevision != 44 ||
	summary.visible ||
	!summary.highlighted ||
	summary.lineStyle != 2 ||
	summary.lineWidth != 6 ||
	fabsf(summary.transparency - 0.35f) > 1.0e-6f ||
	!summary.colorOverride ||
	fabsf(summary.color[0] - sourceColor[0]) > 1.0e-6f ||
	fabsf(summary.color[1] - sourceColor[1]) > 1.0e-6f ||
	fabsf(summary.color[2] - sourceColor[2]) > 1.0e-6f ||
	!summary.materialColorValid ||
	fabsf(summary.materialColor[0] - sourceMaterial[0]) > 1.0e-6f ||
	fabsf(summary.materialColor[1] - sourceMaterial[1]) > 1.0e-6f ||
	fabsf(summary.materialColor[2] - sourceMaterial[2]) > 1.0e-6f ||
	summary.materialRevision != 77) {
	printf("FAIL: scene controller should own database source display/material state updates\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    const uint64_t afterSourceStateFrameRevision =
	ownedScene.getFrameRevision();
    BRLObolDatabaseSourceDisplayPatch sourcePatch;
    sourcePatch.lineWidthValid = TRUE;
    sourcePatch.lineWidth = 9;
    sourcePatch.transparencyValid = TRUE;
    sourcePatch.transparency = 0.62f;
    sourcePatch.colorOverrideValid = TRUE;
    sourcePatch.colorOverride = FALSE;
    sourcePatch.colorValid = TRUE;
    sourcePatch.color = SbColor(0.9f, 0.8f, 0.7f);
    if (ownedScene.setDatabaseSourceDisplayPatch("lod-submit.bot",
	    sourcePatch) != 1 ||
	ownedScene.getStructuralRevision() != afterMoveStructuralRevision ||
	ownedScene.getFrameRevision() <= afterSourceStateFrameRevision ||
	!ownedScene.getDatabaseSourceSummary(0, summary) ||
	summary.sourceRevision != 24 ||
	summary.inputsRevision != 44 ||
	summary.lineWidth != 9 ||
	fabsf(summary.transparency - 0.62f) > 1.0e-6f ||
	summary.colorOverride) {
	printf("FAIL: scene controller source display patch should update presentation without source revision changes\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    const uint64_t afterSourcePatchFrameRevision =
	ownedScene.getFrameRevision();
    if (ownedScene.setDatabaseSourceDisplayPatch("lod-submit.bot",
	    sourcePatch) != 0 ||
	ownedScene.getFrameRevision() != afterSourcePatchFrameRevision ||
	ownedScene.setDatabaseSourceDisplayPatch("missing.bot",
		sourcePatch) != -1) {
	printf("FAIL: scene controller source display patch should report no-op and missing updates\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    const uint64_t afterSourceNameBaseFrameRevision =
	ownedScene.getFrameRevision();
    if (ownedScene.setDatabaseSourceDisplayName("lod-submit.bot",
	    "LoD Submit Display") != 1 ||
	ownedScene.getStructuralRevision() != afterMoveStructuralRevision ||
	ownedScene.getFrameRevision() <= afterSourceNameBaseFrameRevision ||
	!ownedScene.getDatabaseSourceSummary(0, summary) ||
	strcmp(summary.displayName.getString(),
	       "LoD Submit Display") != 0 ||
	!ownedScene.getSceneTreeSummaryForPath("lod-submit.bot",
		treeSummary) ||
	strcmp(treeSummary.displayName.getString(),
	       "LoD Submit Display") != 0) {
	printf("FAIL: scene controller should own database source display name\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    const uint64_t afterSourceNameFrameRevision =
	ownedScene.getFrameRevision();
    if (ownedScene.setDatabaseSourceDisplayName("lod-submit.bot",
	    "LoD Submit Display") != 0 ||
	ownedScene.getFrameRevision() != afterSourceNameFrameRevision ||
	ownedScene.setDatabaseSourceDisplayName("missing.bot",
		"Missing Display") != -1) {
	printf("FAIL: scene controller source display name should report no-op and missing updates\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    const uint64_t afterSourceBoundsBaseFrameRevision =
	ownedScene.getFrameRevision();
    const SbVec3f ownedBoundsMin(-2.0f, -3.0f, -4.0f);
    const SbVec3f ownedBoundsMax(5.0f, 6.0f, 7.0f);
    SbBox3f ownedSourceBounds;
    if (ownedScene.setDatabaseSourceBoundsState("lod-submit.bot", TRUE,
	    ownedBoundsMin, ownedBoundsMax) != 1 ||
	ownedScene.getStructuralRevision() != afterMoveStructuralRevision ||
	ownedScene.getFrameRevision() <=
	afterSourceBoundsBaseFrameRevision ||
	!ownedScene.getDatabaseSourceSummary(0, summary) ||
	!summary.sourceBoundsValid ||
	fabsf(summary.sourceBounds.getMin()[0] + 2.0f) > 1.0e-6f ||
	fabsf(summary.sourceBounds.getMin()[1] + 3.0f) > 1.0e-6f ||
	fabsf(summary.sourceBounds.getMin()[2] + 4.0f) > 1.0e-6f ||
	fabsf(summary.sourceBounds.getMax()[0] - 5.0f) > 1.0e-6f ||
	fabsf(summary.sourceBounds.getMax()[1] - 6.0f) > 1.0e-6f ||
	fabsf(summary.sourceBounds.getMax()[2] - 7.0f) > 1.0e-6f ||
	!ownedScene.getSceneSubtreeBounds("lod-submit.bot", TRUE,
					  ownedSourceBounds) ||
	fabsf(ownedSourceBounds.getMin()[0] + 2.0f) > 1.0e-6f ||
	fabsf(ownedSourceBounds.getMax()[2] - 7.0f) > 1.0e-6f) {
	printf("FAIL: scene controller should own explicit database source bounds\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    const uint64_t afterSourceBoundsFrameRevision =
	ownedScene.getFrameRevision();
    if (ownedScene.setDatabaseSourceBoundsState("lod-submit.bot", TRUE,
	    ownedBoundsMin, ownedBoundsMax) != 0 ||
	ownedScene.getFrameRevision() != afterSourceBoundsFrameRevision ||
	ownedScene.setDatabaseSourceBoundsState("missing.bot", TRUE,
		ownedBoundsMin, ownedBoundsMax) != -1) {
	printf("FAIL: scene controller source bounds should report no-op and missing updates\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    if (ownedScene.setDatabaseSourceState("lod-submit.bot", FALSE, 99, 45,
					  FALSE, FALSE, TRUE, 2, 6, 0.35f, TRUE, sourceColor, TRUE,
					  sourceMaterial, 77) != 1 ||
	ownedScene.getStructuralRevision() != afterMoveStructuralRevision ||
	ownedScene.getFrameRevision() <= afterSourceStateFrameRevision ||
	!ownedScene.getDatabaseSourceSummary(0, summary) ||
	summary.sourceRevision != 24 ||
	summary.inputsRevision != 45) {
	printf("FAIL: scene controller source state should preserve source revision unless explicitly valid\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    const uint64_t afterRevisionPreserveFrameRevision =
	ownedScene.getFrameRevision();
    if (ownedScene.setDatabaseSourceState("lod-submit.bot", TRUE, 24, 45,
					  FALSE, FALSE, TRUE, 2, 6, 0.35f, TRUE, sourceColor, TRUE,
					  sourceMaterial, 77) != 0 ||
	ownedScene.getFrameRevision() != afterRevisionPreserveFrameRevision ||
	ownedScene.setDatabaseSourceState("missing.bot", TRUE, 1, 1,
					  TRUE, FALSE, FALSE, 0, 0, 0.0f, FALSE,
					  SbColor(1.0f, 1.0f, 1.0f), FALSE,
					  SbColor(1.0f, 1.0f, 1.0f), 0) != -1) {
	printf("FAIL: scene controller source state should report no-op and missing updates\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    const uint64_t afterMaterialPolicyFrameRevision =
	ownedScene.getFrameRevision();
    if (ownedScene.setDatabaseSourceMaterialPolicy("lod-submit.bot",
	    SoBRLDatabaseSource::MATERIAL_DATABASE) != 1 ||
	ownedScene.getStructuralRevision() != afterMoveStructuralRevision ||
	ownedScene.getFrameRevision() <= afterMaterialPolicyFrameRevision ||
	!ownedScene.getDatabaseSourceSummary(0, summary) ||
	summary.materialPolicy != SoBRLDatabaseSource::MATERIAL_DATABASE ||
	ownedScene.setDatabaseSourceMaterialPolicy("lod-submit.bot",
		SoBRLDatabaseSource::MATERIAL_DATABASE) != 0 ||
	ownedScene.setDatabaseSourceMaterialPolicy("missing.bot",
		SoBRLDatabaseSource::MATERIAL_DATABASE) != -1) {
	printf("FAIL: scene controller should own database source material policy\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    (void)ownedScene.realizePending();
    const uint64_t beforeStaleFrameRevision = ownedScene.getFrameRevision();
    if (ownedScene.markDatabaseSourceStale("lod-submit.bot",
					   SoBRLDatabaseSource::STALE_INPUTS) != 1 ||
	ownedScene.getFrameRevision() <= beforeStaleFrameRevision ||
	!ownedScene.getDatabaseSourceSummary(0, summary) ||
	!summary.stale ||
	!(summary.staleReason & SoBRLDatabaseSource::STALE_INPUTS) ||
	summary.realizationStatus != SoBRLDatabaseSource::UNREALIZED ||
	ownedScene.markDatabaseSourceStale("lod-submit.bot",
					   SoBRLDatabaseSource::STALE_INPUTS) != 0 ||
	ownedScene.markDatabaseSourceStale("missing.bot",
					   SoBRLDatabaseSource::STALE_INPUTS) != -1) {
	printf("FAIL: scene controller should own database source stale marking\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    if (ownedScene.setDatabaseSourceDrawMode("lod-submit.bot",
	    SoBRLDatabaseSource::WIREFRAME) < 0 ||
	ownedScene.setDatabaseSourceRealizationRoleFlags("lod-submit.bot",
		SoBRLDatabaseSource::REALIZATION_ROLE_MESH) < 0 ||
	!ownedScene.realizePending() ||
	!ownedScene.getDatabaseSourceSummary(0, summary) ||
	summary.stale ||
	summary.realizationStatus != SoBRLDatabaseSource::REALIZED ||
	summary.realizedMeshCount != 1 ||
	summary.realizedShapeCount != 0) {
	printf("FAIL: mesh realization role should realize database mesh even in wire draw mode\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    {
	SoBRLMeshShape *realizedMesh = ownedSource->getRealizedMesh();
	const SoBRLMeshShape *realizedMeshGeom = realizedMesh ?
	    realizedMesh->getGeometrySource() : NULL;
	if (!realizedMesh || !realizedMeshGeom ||
	    realizedMeshGeom->point.getNum() == 0 ||
	    realizedMeshGeom->coordIndex.getNum() == 0) {
	    printf("FAIL: mesh realization role should expose database mesh geometry\n");
	    ownedRoot->unref();
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    return 1;
	}
    }

    if (ownedScene.replaceDatabaseSource("lod-submit.bot", dbip,
					 SoBRLDatabaseSource::WIREFRAME, 25) != 1 ||
	ownedScene.getDatabaseSourceCount() != 1 ||
	ownedScene.getDatabaseSource(0) != ownedSource ||
	ownedScene.findDatabaseSource("lod-submit.bot") != ownedSource ||
	ownedSource->drawMode.getValue() != SoBRLDatabaseSource::WIREFRAME ||
	ownedSource->sourceRevision.getValue() != 25 ||
	!ownedScene.getDatabaseSourceSummary(0, summary) ||
	summary.drawTreeDepth != 3 ||
	strcmp(summary.parentGroupPath.getString(), "draw/group") != 0 ||
	ownedScene.getGroupChildCount("draw/group") != 1) {
	printf("FAIL: scene controller should replace nested database sources in place\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    if (ownedScene.moveDatabaseSourceToGroup("lod-submit.bot", "/") != 1 ||
	ownedScene.getDatabaseSourceCount() != 1 ||
	ownedScene.getDatabaseSource(0) != ownedSource ||
	ownedScene.getGroupChildCount("draw/group") != 0 ||
	!ownedScene.getDatabaseSourceSummary(0, summary) ||
	!summary.hasParent ||
	summary.drawTreeDepth != 1 ||
	strcmp(summary.parentGroupPath.getString(), "/") != 0 ||
	ownedScene.removeGroup("draw") != 1 ||
	ownedRoot->getNumChildren() != 2) {
	printf("FAIL: scene controller should move nested database sources back to the scene root\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    const uint64_t afterReplaceStructuralRevision =
	ownedScene.getStructuralRevision();
    const uint64_t afterReplaceFrameRevision = ownedScene.getFrameRevision();

    if (ownedScene.removeDatabaseSource("lod-submit.bot") != 1 ||
	ownedScene.getDatabaseSourceCount() != 0 ||
	ownedRoot->getNumChildren() != 1 ||
	ownedScene.getStructuralRevision() <= afterReplaceStructuralRevision ||
	ownedScene.getFrameRevision() <= afterReplaceFrameRevision ||
	!ownedScene.getSceneSummary(sceneSummary) ||
	sceneSummary.rootChildCount != 1 ||
	sceneSummary.databaseSourceCount != 0 ||
	sceneSummary.nonDatabaseRootChildCount != 1 ||
	ownedScene.removeDatabaseSource("lod-submit.bot") != 0) {
	printf("FAIL: scene controller should remove database sources while preserving non-source children\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    const uint64_t afterRemoveStructuralRevision =
	ownedScene.getStructuralRevision();
    const uint64_t afterRemoveFrameRevision = ownedScene.getFrameRevision();

    if (ownedScene.replaceDatabaseSource("lod-submit.bot", dbip,
					 SoBRLDatabaseSource::WIREFRAME, 23) != 1 ||
	ownedScene.replaceDatabaseSource("other-submit.bot", dbip,
					 SoBRLDatabaseSource::SHADED, 24) != 1 ||
	ownedScene.moveDatabaseSourceToGroup("other-submit.bot",
		"draw/group") != 1 ||
	ownedScene.getDatabaseSourceCount() != 2 ||
	ownedScene.clearDatabaseSources() != 2 ||
	ownedScene.getDatabaseSourceCount() != 0 ||
	ownedScene.getGroupChildCount("draw/group") != 0 ||
	ownedScene.removeGroup("draw") != 1 ||
	ownedRoot->getNumChildren() != 1 ||
	ownedScene.getStructuralRevision() <= afterRemoveStructuralRevision ||
	ownedScene.getFrameRevision() <= afterRemoveFrameRevision ||
	!ownedScene.getSceneSummary(sceneSummary) ||
	sceneSummary.rootChildCount != 1 ||
	sceneSummary.databaseSourceCount != 0 ||
	sceneSummary.nonDatabaseRootChildCount != 1) {
	printf("FAIL: scene controller should clear only database source children\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    const uint64_t afterClearStructuralRevision =
	ownedScene.getStructuralRevision();
    const uint64_t afterClearFrameRevision = ownedScene.getFrameRevision();
    if (ownedScene.clearDatabaseSources() != 0 ||
	ownedScene.getStructuralRevision() != afterClearStructuralRevision ||
	ownedScene.getFrameRevision() != afterClearFrameRevision ||
	!ownedScene.getSceneSummary(sceneSummary) ||
	sceneSummary.rootChildCount != 1 ||
	sceneSummary.databaseSourceCount != 0 ||
	sceneSummary.nonDatabaseRootChildCount != 1 ||
	sceneSummary.structuralRevision != afterClearStructuralRevision ||
	sceneSummary.frameRevision != afterClearFrameRevision) {
	printf("FAIL: scene controller no-op clear should not advance scene revisions\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    const char *sharedInstanceKey = "scope/shared/lod-submit.bot";
    const char *viewInstanceKey = "scope/view/V0/lod-submit.bot";
    if (ownedScene.replaceDatabaseSourceInstance(sharedInstanceKey,
	    "lod-submit.bot", dbip, SoBRLDatabaseSource::WIREFRAME,
	    31) != 1 ||
	ownedScene.replaceDatabaseSourceInstance(viewInstanceKey,
		"lod-submit.bot", dbip, SoBRLDatabaseSource::SHADED,
		32) != 1 ||
	ownedScene.getDatabaseSourceCount() != 2 ||
	!ownedScene.findDatabaseSourceInstance(sharedInstanceKey) ||
	!ownedScene.findDatabaseSourceInstance(viewInstanceKey)) {
	printf("FAIL: scene controller should own duplicate database paths by source instance key\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    BRLObolDatabaseSourceSummary sharedSummary;
    BRLObolDatabaseSourceSummary viewSummary;
    if (!find_database_source_summary_by_instance(ownedScene,
	    sharedInstanceKey, sharedSummary) ||
	!find_database_source_summary_by_instance(ownedScene,
		viewInstanceKey, viewSummary) ||
	strcmp(sharedSummary.path.getString(), "lod-submit.bot") != 0 ||
	strcmp(viewSummary.path.getString(), "lod-submit.bot") != 0 ||
	strcmp(sharedSummary.instanceKey.getString(),
	       sharedInstanceKey) != 0 ||
	strcmp(viewSummary.instanceKey.getString(), viewInstanceKey) != 0 ||
	sharedSummary.drawMode != SoBRLDatabaseSource::WIREFRAME ||
	viewSummary.drawMode != SoBRLDatabaseSource::SHADED) {
	printf("FAIL: scene controller source summaries should separate instance key from database path\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    if (ownedScene.moveDatabaseSourceInstanceToGroup(viewInstanceKey,
	    "views/V0") != 1 ||
	ownedScene.setDatabaseSourceInstanceState(sharedInstanceKey, TRUE,
		33, 44, FALSE, FALSE, TRUE, 3, 11, 0.5f, TRUE,
		SbColor(0.1f, 0.2f, 0.3f), FALSE,
		SbColor(1.0f, 1.0f, 1.0f), 0) != 1 ||
	!find_database_source_summary_by_instance(ownedScene,
		sharedInstanceKey, sharedSummary) ||
	!find_database_source_summary_by_instance(ownedScene,
		viewInstanceKey, viewSummary) ||
	sharedSummary.visible ||
	!sharedSummary.highlighted ||
	sharedSummary.lineWidth != 11 ||
	viewSummary.visible != TRUE ||
	viewSummary.highlighted != FALSE ||
	strcmp(viewSummary.parentGroupPath.getString(), "views/V0") != 0) {
	printf("FAIL: scene controller source instance state should mutate only the targeted owner\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    if (ownedScene.setDatabaseSourceInstanceRealizationRoleFlags(
	    sharedInstanceKey,
	    SoBRLDatabaseSource::REALIZATION_ROLE_MESH) < 0 ||
	ownedScene.setDatabaseSourceInstanceRealizationRoleFlags(
	    viewInstanceKey,
	    SoBRLDatabaseSource::REALIZATION_ROLE_MESH) < 0 ||
	!ownedScene.realizePending()) {
	printf("FAIL: scene controller should realize duplicate source instances independently\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    int sharedRealized = 0;
    int viewRealized = 0;
    const int realizedCount = ownedScene.getRealizedShapeSummaryCount();
    for (int i = 0; i < realizedCount; i++) {
	BRLObolRealizedShapeSummary instanceShapeSummary;
	if (!ownedScene.getRealizedShapeSummary(i, instanceShapeSummary) ||
	    !instanceShapeSummary.valid)
	    continue;
	const char *instanceShapePath =
	    instanceShapeSummary.path.getString();
	const int instancePathMatches =
	    strcmp(instanceShapePath, "lod-submit.bot") == 0 ||
	    strcmp(instanceShapePath, "/lod-submit.bot") == 0;
	if (strcmp(instanceShapeSummary.ownerSourceInstanceKey.getString(),
		   sharedInstanceKey) == 0 &&
	    instancePathMatches &&
	    strstr(instanceShapeSummary.sourceIdentity.getString(),
		   sharedInstanceKey))
	    sharedRealized = 1;
	if (strcmp(instanceShapeSummary.ownerSourceInstanceKey.getString(),
		   viewInstanceKey) == 0 &&
	    instancePathMatches &&
	    strstr(instanceShapeSummary.sourceIdentity.getString(),
		   viewInstanceKey))
	    viewRealized = 1;
    }
    if (!sharedRealized || !viewRealized) {
	printf("FAIL: realized source summaries should preserve scoped owner instance identity\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    if (ownedScene.removeDatabaseSourceInstance(sharedInstanceKey) != 1 ||
	ownedScene.getDatabaseSourceCount() != 1 ||
	ownedScene.findDatabaseSourceInstance(sharedInstanceKey) ||
	!ownedScene.findDatabaseSourceInstance(viewInstanceKey) ||
	!find_database_source_summary_by_instance(ownedScene,
		viewInstanceKey, viewSummary) ||
	strcmp(viewSummary.path.getString(), "lod-submit.bot") != 0 ||
	strcmp(viewSummary.parentGroupPath.getString(), "views/V0") != 0) {
	printf("FAIL: scene controller should remove only the targeted source instance\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    ownedRoot->unref();
    db_close(dbip);
    bu_file_delete(dbpath);
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
    task.realize = mesh_payload_task;
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

    BRLObolViewLodState viewState;
    SoBRLLodUpdateAction update;
    update.setViewLodState(&viewState);
    if (update.drainService(service) != 1 ||
	update.getResultCount() != 1) {
	printf("FAIL: LoD update action did not drain service result\n");
	service.stop();
	root->unref();
	return 1;
    }

    update.apply(root);
    const BRLObolViewLodState::MeshPayload *payload =
	viewState.findMesh(mesh);
    if (update.getAppliedResultCount() != 1 ||
	update.getUnmatchedResultCount() != 0 ||
	!payload ||
	payload->getTriangleCount() != 2 ||
	mesh->lodStagedAvailable.getValue() ||
	mesh->lodAvailable.getValue() ||
	mesh->isLodDisplayActive() ||
	mesh->getTriangleCount() != 1) {
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
			 "brlobol_mesh_lod",
			 "brlobol-cache-v1",
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
	strcmp(request.providerId.getString(), "brlobol_mesh_lod") != 0 ||
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
			 "brlobol_mesh_lod",
			 "brlobol-cache-v1",
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
    camera->height = 15.0f;
    BRLObolViewController controller(root, camera);
    controller.setViewportSize(320, 240);

    struct bv_view_info info = BV_VIEW_INFO_INIT;
    if (!controller.getViewInfo(&info) ||
	info.width != 320 ||
	info.height != 240 ||
	fabs(info.size - 20.0) > 1.0e-6) {
	printf("FAIL: view controller did not produce RT view info\n");
	root->unref();
	mesh->unref();
	return 1;
    }

    struct bv_context *viewCtx = bv_context_create();
    struct bv *view = bv_context_view(viewCtx);
    if (!viewCtx || !view ||
	!bv_dimensions_set(view, 400, 200) ||
	!bv_size_set(view, 80.0)) {
	printf("FAIL: could not prepare libbv view context for Obol camera sync\n");
	if (viewCtx)
	    bv_context_destroy(viewCtx);
	root->unref();
	mesh->unref();
	return 1;
    }
    bv_context_update(viewCtx, BV_CONTEXT_CHANGED_VIEW);
    controller.setViewportSize(400, 200);
    if (!controller.syncCameraFromViewContext(viewCtx) ||
	!controller.getCamera() ||
	!controller.getCamera()->isOfType(
	    SoOrthographicCamera::getClassTypeId())) {
	printf("FAIL: view controller did not sync libbv view camera\n");
	bv_context_destroy(viewCtx);
	root->unref();
	mesh->unref();
	return 1;
    }
    SoOrthographicCamera *syncedCamera =
	static_cast<SoOrthographicCamera *>(controller.getCamera());
    if (fabs(syncedCamera->height.getValue() - 40.0) > 1.0e-6 ||
	fabs(syncedCamera->aspectRatio.getValue() - 2.0) > 1.0e-6 ||
	!controller.getViewInfo(&info) ||
	info.width != 400 ||
	info.height != 200 ||
	fabs(info.size - 80.0) > 1.0e-6) {
	printf("FAIL: direct Obol camera did not preserve horizontal libbv view size\n");
	bv_context_destroy(viewCtx);
	root->unref();
	mesh->unref();
	return 1;
    }

    point_t center = {10.0, -5.0, 2.0};
    vect_t aet = {35.0, 25.0, 0.0};
    if (!bv_center_set(view, center) ||
	!bv_aet_set(view, aet) ||
	!bv_context_update(viewCtx, BV_CONTEXT_CHANGED_VIEW) ||
	!controller.syncCameraFromViewContext(viewCtx)) {
	printf("FAIL: could not prepare oriented libbv view context for Obol camera sync\n");
	bv_context_destroy(viewCtx);
	root->unref();
	mesh->unref();
	return 1;
    }

    mat_t model2view;
    MAT_IDN(model2view);
    if (!bv_model2view_get(model2view, view)) {
	printf("FAIL: could not query libbv model2view matrix\n");
	bv_context_destroy(viewCtx);
	root->unref();
	mesh->unref();
	return 1;
    }

    const point_t samples[] = {
	{10.0, -5.0, 2.0},
	{20.0, -5.0, 2.0},
	{10.0, 5.0, 2.0},
	{10.0, -5.0, 12.0},
	{0.0, -15.0, -8.0}
    };
    const SbViewVolume viewVolume = syncedCamera->getViewVolume(0.0f);
    for (size_t i = 0; i < sizeof(samples) / sizeof(samples[0]); i++) {
	point_t brlView;
	MAT4X3PNT(brlView, model2view, samples[i]);
	SbVec3f obolScreen;
	viewVolume.projectToScreen(SbVec3f(
				       static_cast<float>(samples[i][X]),
				       static_cast<float>(samples[i][Y]),
				       static_cast<float>(samples[i][Z])), obolScreen);
	const double expectedX = 0.5 * (brlView[X] + 1.0);
	const double expectedY = 0.5 * (brlView[Y] * 2.0 + 1.0);
	if (fabs(obolScreen[0] - expectedX) > 1.0e-5 ||
	    fabs(obolScreen[1] - expectedY) > 1.0e-5) {
	    printf("FAIL: Obol camera projection mismatch for sample %zu: "
		   "BRL=(%.9f, %.9f) Obol=(%.9f, %.9f)\n",
		   i, expectedX, expectedY,
		   static_cast<double>(obolScreen[0]),
		   static_cast<double>(obolScreen[1]));
	    bv_context_destroy(viewCtx);
	    root->unref();
	    mesh->unref();
	    return 1;
	}
    }
    bv_context_destroy(viewCtx);

    controller.setCamera(NULL);
    if (controller.getViewInfo(&info) ||
	info.width != 400 ||
	info.height != 200 ||
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
test_export_record_metadata(void)
{
    SoSeparator *root = new SoSeparator;
    root->ref();

    SoBRLMeshShape *mesh = make_mesh("/metadata/mesh.bot", "mesh.bot");
    mesh->sourceType = "bot";
    mesh->sourceId = 501;
    mesh->displayName = "Mesh Display";
    mesh->geometryName = "mesh geometry";
    mesh->cacheIdentity = "cache://metadata/mesh";
    mesh->sourceIdentity = "db://metadata/mesh.bot";
    mesh->databaseIntent = TRUE;
    mesh->sharedSource = FALSE;
    mesh->drawMode = BRLOBOL_LOD_DRAW_SHADED;
    mesh->recordRole = "database";
    mesh->geometryKind = "surface";
    root->addChild(mesh);

    SoBRLVListShape *vlist = new SoBRLVListShape;
    SbVec3f points[3] = {
	SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 0.0f, 0.0f),
	SbVec3f(2.0f, 0.0f, 0.0f)
    };
    int32_t commands[3] = {
	SoBRLVListShape::MOVE,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::POINT
    };
    vlist->sourcePath = "/metadata/curve.sketch";
    vlist->sourceName = "curve.sketch";
    vlist->sourceType = "sketch";
    vlist->sourceId = 502;
    vlist->displayName = "Curve Display";
    vlist->geometryName = "curve geometry";
    vlist->cacheIdentity = "cache://metadata/curve";
    vlist->sourceIdentity = "db://metadata/curve.sketch";
    vlist->databaseIntent = TRUE;
    vlist->sharedSource = FALSE;
    vlist->drawMode = BRLOBOL_LOD_DRAW_WIRE;
    vlist->recordRole = "database";
    vlist->setLineSet(points, commands, 3);
    root->addChild(vlist);

    SoBRLVListShape *vlistShaded = new SoBRLVListShape;
    vlistShaded->sourcePath = "/metadata/curve.sketch";
    vlistShaded->sourceName = "curve.sketch";
    vlistShaded->sourceType = "sketch";
    vlistShaded->sourceId = 502;
    vlistShaded->displayName = "Curve Display";
    vlistShaded->geometryName = "curve geometry";
    vlistShaded->cacheIdentity = "cache://metadata/curve";
    vlistShaded->sourceIdentity = "db://metadata/curve.sketch";
    vlistShaded->databaseIntent = TRUE;
    vlistShaded->sharedSource = FALSE;
    vlistShaded->drawMode = BRLOBOL_LOD_DRAW_SHADED;
    vlistShaded->recordRole = "database";
    vlistShaded->setLineSet(points, commands, 2);
    root->addChild(vlistShaded);

    SoBRLLodMeshShape *sourceBacked =
	make_lod_mesh("/metadata/source.bot", "source.bot");
    sourceBacked->sourceType = "bot";
    sourceBacked->sourceId = 503;
    sourceBacked->displayName = "Source Mesh Display";
    sourceBacked->geometryName = "source mesh geometry";
    sourceBacked->cacheIdentity = "cache://metadata/source";
    sourceBacked->sourceIdentity = "db://metadata/source.bot";
    sourceBacked->databaseIntent = TRUE;
    sourceBacked->sharedSource = FALSE;
    sourceBacked->drawMode = BRLOBOL_LOD_DRAW_SHADED;
    sourceBacked->recordRole = "database";
    sourceBacked->geometryKind = "surface";
    root->addChild(sourceBacked);

    SoBRLExportAction exportAction;
    exportAction.apply(root);
    if (exportAction.getTriangleCount() != 1 ||
	exportAction.getPointCount() != 1 ||
	exportAction.getLineCount() != 2 ||
	exportAction.getSourceBackedFullDetailRequestCount() != 1) {
	printf("FAIL: export metadata test did not collect expected records\n");
	root->unref();
	return 1;
    }

    const uint64_t meshCacheIdentity =
	SoBRLExportAction::identityValue("cache://metadata/mesh");
    const uint64_t meshSourceIdentity =
	SoBRLExportAction::identityValue(SbString("db://metadata/mesh.bot"));
    const uint64_t curveCacheIdentity =
	SoBRLExportAction::identityValue("cache://metadata/curve");
    const uint64_t curveSourceIdentity =
	SoBRLExportAction::identityValue("db://metadata/curve.sketch");
    const uint64_t sourceCacheIdentity =
	SoBRLExportAction::identityValue("cache://metadata/source");
    const uint64_t sourceSourceIdentity =
	SoBRLExportAction::identityValue("db://metadata/source.bot");
    if (SoBRLExportAction::identityValue("") != 0 ||
	meshCacheIdentity == 0 ||
	meshSourceIdentity == 0 ||
	meshCacheIdentity !=
	SoBRLExportAction::identityValue("cache://metadata/mesh") ||
	meshCacheIdentity == curveCacheIdentity ||
	meshSourceIdentity == curveSourceIdentity) {
	printf("FAIL: export identity helper did not produce stable numeric identities\n");
	root->unref();
	return 1;
    }

    const SoBRLExportAction::TriangleRecord &meshRecord =
	exportAction.getTriangle(0);
    if (strcmp(meshRecord.displayName.getString(), "Mesh Display") != 0 ||
	strcmp(meshRecord.geometryName.getString(), "mesh geometry") != 0 ||
	strcmp(meshRecord.cacheIdentity.getString(),
	       "cache://metadata/mesh") != 0 ||
	strcmp(meshRecord.sourceIdentity.getString(),
	       "db://metadata/mesh.bot") != 0 ||
	meshRecord.cacheIdentityValue != meshCacheIdentity ||
	meshRecord.sourceIdentityValue != meshSourceIdentity ||
	!meshRecord.databaseIntent ||
	meshRecord.overlayIntent ||
	meshRecord.localSource ||
	meshRecord.sharedSource ||
	meshRecord.nonDatabaseSource ||
	meshRecord.drawMode != BRLOBOL_LOD_DRAW_SHADED ||
	strcmp(meshRecord.recordRole.getString(), "database") != 0 ||
	strcmp(meshRecord.geometryKind.getString(), "surface") != 0) {
	printf("FAIL: mesh export record did not preserve neutral metadata\n");
	root->unref();
	return 1;
    }

    const SoBRLExportAction::LineRecord &lineRecord = exportAction.getLine(0);
    const SoBRLExportAction::PointRecord &pointRecord =
	exportAction.getPoint(0);
    if (strcmp(lineRecord.displayName.getString(), "Curve Display") != 0 ||
	strcmp(lineRecord.geometryName.getString(), "curve geometry") != 0 ||
	strcmp(lineRecord.geometryKind.getString(), "line") != 0 ||
	lineRecord.cacheIdentityValue != curveCacheIdentity ||
	lineRecord.sourceIdentityValue != curveSourceIdentity ||
	strcmp(pointRecord.displayName.getString(), "Curve Display") != 0 ||
	strcmp(pointRecord.geometryKind.getString(), "point") != 0 ||
	pointRecord.cacheIdentityValue != curveCacheIdentity ||
	pointRecord.sourceIdentityValue != curveSourceIdentity ||
	pointRecord.drawMode != BRLOBOL_LOD_DRAW_WIRE ||
	!pointRecord.databaseIntent ||
	pointRecord.sharedSource) {
	printf("FAIL: vlist export records did not preserve neutral metadata\n");
	root->unref();
	return 1;
    }

    if (exportAction.getObjectRecordCount() != 3) {
	printf("FAIL: export object records did not group initial objects\n");
	root->unref();
	return 1;
    }
    const SoBRLExportAction::ObjectRecord &meshObject =
	exportAction.getObjectRecord(0);
    const SoBRLExportAction::ObjectRecord &vlistObject =
	exportAction.getObjectRecord(1);
    const SoBRLExportAction::ObjectRecord &vlistShadedObject =
	exportAction.getObjectRecord(2);
    if (strcmp(meshObject.path.getString(), "/metadata/mesh.bot") != 0 ||
	meshObject.triangleIndices.size() != 1 ||
	!meshObject.lineIndices.empty() ||
	!meshObject.pointIndices.empty() ||
	!meshObject.databaseIntent ||
	meshObject.sharedSource ||
	meshObject.cacheIdentityValue != meshCacheIdentity ||
	meshObject.sourceIdentityValue != meshSourceIdentity ||
	strcmp(meshObject.geometryKind.getString(), "surface") != 0 ||
	meshObject.bounds.isEmpty()) {
	printf("FAIL: mesh object export record did not summarize triangle metadata\n");
	root->unref();
	return 1;
    }
    if (strcmp(vlistObject.path.getString(), "/metadata/curve.sketch") != 0 ||
	vlistObject.lineIndices.size() != 1 ||
	vlistObject.pointIndices.size() != 1 ||
	!vlistObject.triangleIndices.empty() ||
	strcmp(vlistObject.geometryKind.getString(), "mixed") != 0 ||
	vlistObject.cacheIdentityValue != curveCacheIdentity ||
	vlistObject.sourceIdentityValue != curveSourceIdentity ||
	vlistObject.drawMode != BRLOBOL_LOD_DRAW_WIRE) {
	printf("FAIL: vlist object export record did not group line and point records\n");
	root->unref();
	return 1;
    }
    if (strcmp(vlistShadedObject.path.getString(),
	       "/metadata/curve.sketch") != 0 ||
	vlistShadedObject.lineIndices.size() != 1 ||
	!vlistShadedObject.pointIndices.empty() ||
	!vlistShadedObject.triangleIndices.empty() ||
	vlistShadedObject.drawMode != BRLOBOL_LOD_DRAW_SHADED ||
	vlistShadedObject.cacheIdentityValue != curveCacheIdentity ||
	vlistShadedObject.sourceIdentityValue != curveSourceIdentity) {
	printf("FAIL: object export records should not merge different draw modes\n");
	root->unref();
	return 1;
    }
    SoBRLExportAction::ObjectSurfaceSummary surfaceSummary;
    int surfaceIndex = -1;
    if (!exportAction.getObjectRecordSurfaceSummary(meshObject,
	    surfaceSummary) ||
	!surfaceSummary.valid ||
	surfaceSummary.pointCount != 3 ||
	surfaceSummary.indexCount != 3 ||
	surfaceSummary.faceCount != 1 ||
	surfaceSummary.cacheIdentityValue != meshCacheIdentity ||
	surfaceSummary.sourceIdentityValue != meshSourceIdentity ||
	!exportAction.getObjectRecordSurfaceIndex(meshObject, 2,
		surfaceIndex) ||
	surfaceIndex != 2) {
	printf("FAIL: mesh object export record did not provide surface detail readback\n");
	root->unref();
	return 1;
    }

    SoBRLExportAction::ObjectLineSummary lineSummary;
    SoBRLExportAction::ObjectPointSummary pointSummary;
    SbVec3f detailPoint;
    int lineCommand = -1;
    if (!exportAction.getObjectRecordLineSummary(vlistObject,
	    lineSummary) ||
	lineSummary.pointCount != 2 ||
	lineSummary.segmentCount != 1 ||
	!exportAction.getObjectRecordLinePoint(vlistObject, 0,
		detailPoint) ||
	!NEAR_ZERO(detailPoint[0], SMALL_FASTF) ||
	!exportAction.getObjectRecordLineCommand(vlistObject, 0,
		lineCommand) ||
	lineCommand != SoBRLExportAction::LINE_MOVE ||
	!exportAction.getObjectRecordLineCommand(vlistObject, 1,
		lineCommand) ||
	lineCommand != SoBRLExportAction::LINE_DRAW ||
	!exportAction.getObjectRecordPointSummary(vlistObject,
		pointSummary) ||
	pointSummary.pointCount != 1 ||
	lineSummary.cacheIdentityValue != curveCacheIdentity ||
	lineSummary.sourceIdentityValue != curveSourceIdentity ||
	pointSummary.cacheIdentityValue != curveCacheIdentity ||
	pointSummary.sourceIdentityValue != curveSourceIdentity ||
	!exportAction.getObjectRecordPoint(vlistObject, 0,
					   detailPoint) ||
	!NEAR_EQUAL(detailPoint[0], 2.0f, SMALL_FASTF)) {
	printf("FAIL: vlist object export record did not provide line/point detail readback\n");
	root->unref();
	return 1;
    }
    std::vector<SoBRLExportAction::ObjectRecord> queriedObjects;
    if (exportAction.collectObjectRecords(queriedObjects,
					  SoBRLExportAction::QUERY_VISIBLE_ONLY) != 3 ||
	exportAction.collectObjectRecords(queriedObjects,
					  SoBRLExportAction::QUERY_DATABASE_OBJECTS) != 3 ||
	exportAction.collectObjectRecords(queriedObjects,
					  SoBRLExportAction::QUERY_VIEW_OBJECTS) != 0 ||
	exportAction.collectObjectRecords(queriedObjects,
					  SoBRLExportAction::QUERY_LOCAL_ONLY) != 0 ||
	exportAction.collectObjectRecords(queriedObjects, 0,
					  "overlay::*") != 0 ||
	exportAction.collectObjectRecords(queriedObjects,
					  SoBRLExportAction::QUERY_DATABASE_OBJECTS, NULL,
					  BRLOBOL_LOD_DRAW_WIRE) != 1 ||
	exportAction.collectObjectRecords(queriedObjects,
					  SoBRLExportAction::QUERY_DATABASE_OBJECTS, NULL,
					  BRLOBOL_LOD_DRAW_SHADED) != 2) {
	printf("FAIL: export object record filters did not match expected object sets\n");
	root->unref();
	return 1;
    }

    const BRLObolSourceMeshRequest &sourceRequest =
	exportAction.getSourceBackedFullDetailRequest(0);
    if (strcmp(sourceRequest.displayName.getString(),
	       "Source Mesh Display") != 0 ||
	strcmp(sourceRequest.geometryName.getString(),
	       "source mesh geometry") != 0 ||
	strcmp(sourceRequest.cacheIdentity.getString(),
	       "cache://metadata/source") != 0 ||
	strcmp(sourceRequest.sourceIdentity.getString(),
	       "db://metadata/source.bot") != 0 ||
	!sourceRequest.databaseIntent ||
	sourceRequest.sharedSource ||
	sourceRequest.drawMode != BRLOBOL_LOD_DRAW_SHADED ||
	strcmp(sourceRequest.recordRole.getString(), "database") != 0 ||
	strcmp(sourceRequest.geometryKind.getString(), "surface") != 0) {
	printf("FAIL: source-backed request did not preserve neutral metadata\n");
	root->unref();
	return 1;
    }

    BRLObolLodRequest sourceLodRequest;
    if (!exportAction.makeSourceBackedFullDetailLodRequest(0,
	    sourceLodRequest)) {
	printf("FAIL: export metadata test did not make source-backed LoD request\n");
	root->unref();
	return 1;
    }
    int beforeTriangleCount = exportAction.getTriangleCount();
    if (!exportAction.appendSourceBackedFullDetailResult(sourceRequest,
	    source_full_detail_result(sourceLodRequest)) ||
	exportAction.getTriangleCount() != beforeTriangleCount + 1) {
	printf("FAIL: export metadata test did not append source-backed full-detail result\n");
	root->unref();
	return 1;
    }
    const SoBRLExportAction::TriangleRecord &sourceRecord =
	exportAction.getTriangle(beforeTriangleCount);
    if (strcmp(sourceRecord.displayName.getString(),
	       "Source Mesh Display") != 0 ||
	strcmp(sourceRecord.cacheIdentity.getString(),
	       "cache://metadata/source") != 0 ||
	strcmp(sourceRecord.sourceIdentity.getString(),
	       "db://metadata/source.bot") != 0 ||
	sourceRecord.cacheIdentityValue != sourceCacheIdentity ||
	sourceRecord.sourceIdentityValue != sourceSourceIdentity ||
	!sourceRecord.databaseIntent ||
	sourceRecord.sharedSource ||
	sourceRecord.nonDatabaseSource ||
	sourceRecord.drawMode != BRLOBOL_LOD_DRAW_SHADED ||
	strcmp(sourceRecord.recordRole.getString(), "database") != 0 ||
	strcmp(sourceRecord.geometryKind.getString(), "surface") != 0) {
	printf("FAIL: source-backed full-detail export record did not preserve neutral metadata\n");
	root->unref();
	return 1;
    }

    if (exportAction.getObjectRecordCount() != 4) {
	printf("FAIL: export object records did not rebuild after source-backed full-detail append\n");
	root->unref();
	return 1;
    }
    const SoBRLExportAction::ObjectRecord &sourceObject =
	exportAction.getObjectRecord(3);
    if (strcmp(sourceObject.path.getString(), "/metadata/source.bot") != 0 ||
	sourceObject.triangleIndices.size() != 1 ||
	!sourceObject.databaseIntent ||
	sourceObject.sharedSource ||
	sourceObject.nonDatabaseSource ||
	strcmp(sourceObject.displayName.getString(),
	       "Source Mesh Display") != 0 ||
	sourceObject.cacheIdentityValue != sourceCacheIdentity ||
	sourceObject.sourceIdentityValue != sourceSourceIdentity ||
	strcmp(sourceObject.cacheIdentity.getString(),
	       "cache://metadata/source") != 0) {
	printf("FAIL: source-backed object export record did not preserve neutral metadata\n");
	root->unref();
	return 1;
    }
    if (!exportAction.getObjectRecordSurfaceSummary(sourceObject,
	    surfaceSummary) ||
	surfaceSummary.pointCount != 3 ||
	surfaceSummary.indexCount != 3 ||
	surfaceSummary.faceCount != 1 ||
	surfaceSummary.cacheIdentityValue != sourceCacheIdentity ||
	surfaceSummary.sourceIdentityValue != sourceSourceIdentity ||
	!exportAction.getObjectRecordSurfaceIndex(sourceObject, 1,
		surfaceIndex) ||
	surfaceIndex != 1) {
	printf("FAIL: source-backed object export record did not provide surface detail readback\n");
	root->unref();
	return 1;
    }
    if (exportAction.collectObjectRecords(queriedObjects,
					  SoBRLExportAction::QUERY_DATABASE_OBJECTS) != 4 ||
	exportAction.collectObjectRecords(queriedObjects, 0,
					  "*source.bot") != 1) {
	printf("FAIL: export object record filters did not include source-backed full-detail object\n");
	root->unref();
	return 1;
    }

    root->unref();
    return 0;
}

static int
test_mesh_residency_budget_action(void)
{
    SoSeparator *root = new SoSeparator;
    root->ref();

    SoBRLMeshShape *preservingMesh = make_mesh("/budget/preserve.bot",
				     "preserve.bot");
    SoBRLMeshShape *lodMesh = make_lod_mesh("/budget/lod.bot", "lod.bot");
    preservingMesh->sourceId = 121;
    lodMesh->sourceId = 122;
    root->addChild(preservingMesh);
    root->addChild(lodMesh);

    BRLObolLodRequest preservingRequest =
	make_request("/budget/preserve.bot", "preserve.bot");
    BRLObolLodRequest lodRequest =
	make_request("/budget/lod.bot", "lod.bot");
    BRLObolLodResult preservingResult =
	mesh_payload_result(preservingRequest);
    BRLObolLodResult lodResult =
	mesh_payload_result(lodRequest);
    if (!preservingMesh->applyStagedLodResult(preservingResult,
	    &preservingRequest) ||
	!lodMesh->applyStagedLodResult(lodResult, &lodRequest) ||
	!preservingMesh->hasFullDetailMesh() ||
	lodMesh->hasFullDetailMesh() ||
	!preservingMesh->isLodDisplayActive() ||
	!lodMesh->isLodDisplayActive()) {
	printf("FAIL: mesh residency budget fixture did not stage LoD meshes\n");
	root->unref();
	return 1;
    }

    size_t initialBytes =
	preservingMesh->estimateResidentMeshBytes() +
	lodMesh->estimateResidentMeshBytes();
    size_t fullDetailBytes = preservingMesh->estimateFullDetailMeshBytes();
    size_t displayBytesAfterFullDetail =
	preservingMesh->estimateDisplayMeshBytes() +
	lodMesh->estimateDisplayMeshBytes();
    if (initialBytes == 0 || fullDetailBytes == 0 ||
	displayBytesAfterFullDetail == 0 ||
	initialBytes != fullDetailBytes + displayBytesAfterFullDetail) {
	printf("FAIL: mesh residency budget fixture byte accounting failed\n");
	root->unref();
	return 1;
    }

    SoBRLMeshResidencyAction fullDetailBudget;
    fullDetailBudget.setMaxResidentMeshBytes(displayBytesAfterFullDetail);
    fullDetailBudget.setEvictDisplayPayloads(FALSE);
    fullDetailBudget.apply(root);
    if (fullDetailBudget.getVisitedMeshCount() != 2 ||
	fullDetailBudget.getInitialResidentMeshBytes() != initialBytes ||
	fullDetailBudget.getFinalResidentMeshBytes() !=
	displayBytesAfterFullDetail ||
	fullDetailBudget.getFreedFullDetailBytes() != fullDetailBytes ||
	fullDetailBudget.getFreedDisplayBytes() != 0 ||
	fullDetailBudget.getEvictedFullDetailMeshCount() != 1 ||
	fullDetailBudget.getEvictedDisplayMeshCount() != 0 ||
	preservingMesh->hasFullDetailMesh() ||
	!preservingMesh->isLodDisplayActive() ||
	!lodMesh->isLodDisplayActive()) {
	printf("FAIL: mesh residency budget did not evict preserved full detail first\n");
	root->unref();
	return 1;
    }

    int controllerBudgetRet = 0;
    {
	SoOrthographicCamera *camera = new SoOrthographicCamera;
	BRLObolViewController controller(root, camera);
	controller.clearRenderRequest();
	size_t freedBytes = controller.evictMeshPayloadsToBudget(0, TRUE);
	BRLObolSourceMeshRequest preservingSourceRequest;
	BRLObolSourceMeshRequest lodSourceRequest;
	if (freedBytes != displayBytesAfterFullDetail ||
	    controller.getLastMeshBudgetVisitedMeshCount() != 2 ||
	    controller.getLastMeshBudgetInitialResidentBytes() !=
	    displayBytesAfterFullDetail ||
	    controller.getLastMeshBudgetFinalResidentBytes() != 0 ||
	    controller.getLastMeshBudgetFreedResidentBytes() !=
	    displayBytesAfterFullDetail ||
	    controller.getLastMeshBudgetFreedFullDetailBytes() != 0 ||
	    controller.getLastMeshBudgetFreedDisplayBytes() !=
	    displayBytesAfterFullDetail ||
	    controller.getLastMeshBudgetEvictedFullDetailMeshCount() != 0 ||
	    controller.getLastMeshBudgetEvictedDisplayMeshCount() != 2 ||
	    !controller.isRenderRequested() ||
	    strcmp(controller.getRenderReason().getString(),
		   "lod-memory-budget") != 0 ||
	    preservingMesh->isLodDisplayActive() ||
	    lodMesh->isLodDisplayActive() ||
	    preservingMesh->getTriangleCount() != 0 ||
	    lodMesh->getTriangleCount() != 0 ||
	    preservingMesh->estimateResidentMeshBytes() != 0 ||
	    lodMesh->estimateResidentMeshBytes() != 0 ||
	    !preservingMesh->makeSourceMeshRequest(
		preservingSourceRequest) ||
	    !lodMesh->makeSourceMeshRequest(lodSourceRequest) ||
	    check_source_mesh_request(preservingSourceRequest,
				      "/budget/preserve.bot", "preserve.bot", 121) ||
	    check_source_mesh_request(lodSourceRequest,
				      "/budget/lod.bot", "lod.bot", 122)) {
	    printf("FAIL: controller mesh residency budget did not preserve source-backed exact identity\n");
	    controllerBudgetRet = 1;
	}
    }
    if (controllerBudgetRet) {
	root->unref();
	return 1;
    }

    root->unref();
    return 0;
}

static int
test_view_controller_mesh_residency_budget_auto(void)
{
    SoSeparator *root = new SoSeparator;
    root->ref();
    SoBRLMeshShape *mesh = make_mesh("/budget/auto.bot", "auto.bot");
    mesh->sourceId = 131;
    root->addChild(mesh);

    int ret = 0;
    {
	SoOrthographicCamera *camera = new SoOrthographicCamera;
	BRLObolViewController controller(root, camera);
	BRLObolLodService service;

	if (controller.hasMeshResidencyBudget()) {
	    printf("FAIL: controller mesh residency budget should default disabled\n");
	    ret = 1;
	}
	if (!ret) {
	    controller.setMeshResidencyBudget(0, TRUE);
	    if (!controller.hasMeshResidencyBudget() ||
		controller.getMaxResidentMeshBytes() != 0 ||
		!controller.isMeshResidencyDisplayEvictionEnabled()) {
		printf("FAIL: controller mesh residency budget policy was not stored\n");
		ret = 1;
	    }
	}
	if (!ret && !service.start(1, TRUE)) {
	    printf("FAIL: controller mesh residency budget service did not start\n");
	    ret = 1;
	}

	if (!ret) {
	    BRLObolLodTask task;
	    task.generation = service.beginGeneration();
	    task.request = make_request("/budget/auto.bot", "auto.bot");
	    task.request.viewRevision = controller.getLodViewRevision();
	    task.request.policyRevision = controller.getLodPolicyRevision();
	    task.realize = mesh_payload_task;
	    if (service.submit(task) == 0) {
		printf("FAIL: controller mesh residency budget service did not accept task\n");
		ret = 1;
	    }
	}
	if (!ret && wait_for_service(service))
	    ret = 1;

	if (!ret) {
	    controller.clearRenderRequest();
	    int applied = controller.applyLodResults(&service);
	    BRLObolSourceMeshRequest sourceRequest;
	    if (applied != 1 ||
		controller.getLastLodAppliedResultCount() != 1 ||
		controller.getLastMeshBudgetVisitedMeshCount() != 1 ||
		controller.getLastMeshBudgetInitialResidentBytes() == 0 ||
		controller.getLastMeshBudgetFinalResidentBytes() != 0 ||
		controller.getLastMeshBudgetEvictedFullDetailMeshCount() != 1 ||
		controller.getLastMeshBudgetEvictedDisplayMeshCount() != 1 ||
		!controller.isRenderRequested() ||
		strcmp(controller.getRenderReason().getString(),
		       "lod-memory-budget") != 0 ||
		mesh->isLodDisplayActive() ||
		mesh->hasFullDetailMesh() ||
		mesh->getTriangleCount() != 0 ||
		mesh->estimateResidentMeshBytes() != 0 ||
		!mesh->makeSourceMeshRequest(sourceRequest) ||
		check_source_mesh_request(sourceRequest,
					  "/budget/auto.bot", "auto.bot", 131)) {
		printf("FAIL: controller mesh residency budget did not auto-evict applied LoD payload\n");
		ret = 1;
	    }
	}
	if (!ret) {
	    controller.clearMeshResidencyBudget();
	    if (controller.hasMeshResidencyBudget() ||
		!controller.isMeshResidencyDisplayEvictionEnabled()) {
		printf("FAIL: controller mesh residency budget clear failed\n");
		ret = 1;
	    }
	}

	service.stop();
    }

    root->unref();
    return ret;
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
	brlobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    struct bv_view_info view = BV_VIEW_INFO_INIT;
    view.width = 640;
    view.height = 480;
    view.size = 100.0;

    BRLObolLodRequest activeDuplicateRequest;
    mesh->makeLodRequest(activeDuplicateRequest,
			 "db://lod-submit-test",
			 2026,
			 55,
			 66,
			 BRLOBOL_LOD_DRAW_SHADED,
			 "brlobol_mesh_lod",
			 "brlobol-cache-v1",
			 BRLOBOL_LOD_QUALITY_FAST_DISPLAY);
    BRLObolLodTask activeTask;
    activeTask.generation = service.beginGeneration();
    activeTask.request = activeDuplicateRequest;
    activeTask.realize = aabb_task;
    activeTask.debugDelayMilliseconds = 200;
    if (service.submit(activeTask) == 0) {
	printf("FAIL: LoD submit action active-duplicate fixture did not queue request\n");
	service.stop();
	root->unref();
	brlobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    SoBRLMeshLodSubmitAction activeDuplicateSubmit;
    activeDuplicateSubmit.setService(&service);
    activeDuplicateSubmit.setDatabase(dbip, "db://lod-submit-test", 2026);
    activeDuplicateSubmit.setViewInfo(&view);
    activeDuplicateSubmit.setGeneration(service.beginGeneration());
    activeDuplicateSubmit.setRevisions(55, 66);
    activeDuplicateSubmit.apply(root);
    if (activeDuplicateSubmit.getVisitedMeshCount() != 2 ||
	activeDuplicateSubmit.getSubmittedTaskCount() != 0 ||
	activeDuplicateSubmit.getSkippedMeshCount() != 2 ||
	strstr(activeDuplicateSubmit.getDiagnostics().getString(),
	       "current LoD request is already active") == NULL) {
	printf("FAIL: LoD submit action did not skip active duplicate view/policy request\n");
	service.stop();
	root->unref();
	brlobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }
    if (wait_for_service(service)) {
	service.stop();
	root->unref();
	brlobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }
    {
	std::vector<BRLObolLodResult> activeResults;
	service.drainResults(activeResults);
	if (activeResults.size() != 1 ||
	    !brlobol_lod_result_matches_request(activeResults[0],
						activeDuplicateRequest)) {
	    printf("FAIL: LoD submit action active-duplicate fixture did not preserve queued request\n");
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
    }

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
	brlobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    if (wait_for_service(service)) {
	service.stop();
	root->unref();
	brlobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    std::vector<BRLObolLodResult> viewPolicyResults;
    if (service.drainResults(viewPolicyResults) != 1 ||
	viewPolicyResults.size() != 1) {
	printf("FAIL: LoD submit action result drain failed\n");
	service.stop();
	root->unref();
	brlobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    int expectedViewLevel = -1;
    if (expected_view_lod_level(dbip, "lod-submit.bot", &view,
				&expectedViewLevel)) {
	printf("FAIL: LoD submit action view-policy active level helper failed\n");
	service.stop();
	root->unref();
	brlobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    if (viewPolicyResults[0].geometry.activeLevel != expectedViewLevel) {
	printf("FAIL: LoD submit action did not use view-policy active level\n");
	service.stop();
	root->unref();
	brlobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    BRLObolViewLodState viewState;
    SoBRLViewLodGroup *renderRoot = new SoBRLViewLodGroup;
    renderRoot->ref();
    renderRoot->setViewLodState(&viewState);
    renderRoot->addChild(root);

    SoBRLLodUpdateAction update;
    update.setViewLodState(&viewState);
    update.setResults(viewPolicyResults);
    if (update.getResultCount() != 1) {
	printf("FAIL: LoD submit action result drain failed\n");
	service.stop();
	root->unref();
	brlobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    update.apply(root);
    int ret = 0;
    const BRLObolViewLodState::MeshPayload *activePayload =
	viewState.findMesh(mesh);
    if (update.getAppliedResultCount() != 1 ||
	update.getRejectedResultCount() != 0 ||
	!activePayload ||
	activePayload->resultKind != BRLOBOL_LOD_RESULT_MESH ||
	activePayload->counts.faceCount != 4 ||
	activePayload->counts.pointCount != 4 ||
	activePayload->getTriangleCount() != 4 ||
	mesh->lodAvailable.getValue() ||
	mesh->isLodDisplayActive() ||
	mesh->getTriangleCount() != 1) {
	printf("FAIL: LoD submit action result did not update mesh payload\n");
	ret = 1;
    }
    if (!ret && mesh->hasFullDetailMesh()) {
	printf("FAIL: LoD-backed mesh retained unexpected full-detail payload\n");
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
			     "brlobol_mesh_lod",
			     "brlobol-cache-v1",
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
	SoBRLMeshLodSubmitAction duplicateSubmit;
	duplicateSubmit.setService(&service);
	duplicateSubmit.setViewLodState(&viewState);
	duplicateSubmit.setDatabase(dbip, "db://lod-submit-test", 2026);
	duplicateSubmit.setViewInfo(&view);
	duplicateSubmit.setGeneration(service.beginGeneration());
	duplicateSubmit.setRevisions(55, 66);
	duplicateSubmit.apply(root);
	if (duplicateSubmit.getVisitedMeshCount() != 2 ||
	    duplicateSubmit.getSubmittedTaskCount() != 0 ||
	    duplicateSubmit.getSkippedMeshCount() != 2 ||
	    strstr(duplicateSubmit.getDiagnostics().getString(),
		   "current LoD request is already resident") == NULL ||
	    service.queuedResultCountForDiagnostics() != 0) {
	    printf("FAIL: LoD submit action did not skip already-resident view/policy request\n");
	    ret = 1;
	}
    }

    if (!ret) {
	SoBRLExportAction exactExport;
	exactExport.apply(renderRoot);
	if (exactExport.getGeometryPolicy() != SoBRLExportAction::FULL_DETAIL ||
	    exactExport.getTriangleCount() != 1 ||
	    exactExport.getSkippedLodDisplayMeshCount() != 1 ||
	    exactExport.getSourceBackedFullDetailRequestCount() != 1 ||
	    check_source_mesh_request(
		exactExport.getSourceBackedFullDetailRequest(0),
		"/lod-submit.bot", "lod-submit.bot", 101)) {
	    printf("FAIL: exact export did not request source-backed full-detail LoD mesh\n");
	    ret = 1;
	} else {
	    BRLObolLodRequest sourceLodRequest;
	    if (!exactExport.makeSourceBackedFullDetailLodRequest(0,
		    sourceLodRequest) ||
		check_source_full_detail_lod_request(sourceLodRequest,
			"/lod-submit.bot", "lod-submit.bot")) {
		printf("FAIL: exact export source-backed request did not convert to RT full-detail LoD request\n");
		ret = 1;
	    } else {
		int beforeTriangleCount = exactExport.getTriangleCount();
		BRLObolLodResult sourceResult =
		    source_full_detail_result(sourceLodRequest);
		std::vector<BRLObolLodResult> sourceResults;
		sourceResults.push_back(sourceResult);
		if (exactExport.consumeSourceBackedFullDetailResults(
			sourceResults) != 1 ||
		    exactExport.getTriangleCount() != beforeTriangleCount + 1 ||
		    exactExport.getTriangle(beforeTriangleCount).sourceId != 101 ||
		    exactExport.getTriangle(beforeTriangleCount).vertexIndexA != 0 ||
		    exactExport.getTriangle(beforeTriangleCount).vertexIndexB != 1 ||
		    exactExport.getTriangle(beforeTriangleCount).vertexIndexC != 2) {
		    printf("FAIL: exact export did not consume source-backed full-detail LoD result\n");
		    ret = 1;
		}
		if (!ret) {
		    SoBRLExportAction controllerExport;
		    controllerExport.apply(renderRoot);
		    BRLObolLodRequest controllerExportRequest;
		    if (controllerExport.getSourceBackedFullDetailRequestCount() != 1 ||
			!controllerExport.makeSourceBackedFullDetailLodRequest(0,
				controllerExportRequest)) {
			printf("FAIL: controller source-backed exact export helper did not collect source request\n");
			ret = 1;
		    } else {
			BRLObolLodService controllerExportService;
			if (!controllerExportService.start(1, TRUE)) {
			    printf("FAIL: controller source-backed exact export helper service did not start\n");
			    ret = 1;
			} else {
			    if (submit_source_full_detail_task(
				    controllerExportService,
				    controllerExportRequest)) {
				ret = 1;
			    } else {
				BRLObolViewController exportController(root,
								       NULL);
				exportController.setLodService(
				    &controllerExportService);
				int submittedCount = -1;
				beforeTriangleCount =
				    controllerExport.getTriangleCount();
				if (exportController.consumeExportSourceFullDetail(
					controllerExport, 0, &submittedCount) != 1 ||
				    submittedCount != 0 ||
				    controllerExport.getTriangleCount() !=
				    beforeTriangleCount + 1 ||
				    controllerExport.getTriangle(
					beforeTriangleCount).sourceId != 101) {
				    printf("FAIL: controller source-backed exact export helper did not consume matching LoD result\n");
				    ret = 1;
				}
			    }
			    controllerExportService.stop();
			}
		    }
		}
		BRLObolLodService sourceService;
		if (!sourceService.start(1, TRUE)) {
		    printf("FAIL: exact export source-backed submit helper service did not start\n");
		    ret = 1;
		} else {
		    if (exactExport.submitSourceBackedFullDetailRequests(
			    &sourceService, sourceService.beginGeneration(),
			    dbip) != 1) {
			printf("FAIL: exact export did not submit source-backed full-detail LoD request\n");
			ret = 1;
		    } else if (wait_for_service(sourceService)) {
			ret = 1;
		    } else {
			std::vector<BRLObolLodResult> submittedResults;
			sourceService.drainResults(submittedResults);
			if (submittedResults.size() != 1 ||
			    submittedResults[0].providerStatus !=
			    BRLOBOL_LOD_PROVIDER_STALE) {
			    printf("FAIL: exact export source-backed submit helper did not publish stale source result\n");
			    ret = 1;
			}
		    }
		    sourceService.stop();
		}
	    }
	}
    }

    if (!ret) {
	SoBRLMeasureAction exactMeasure;
	exactMeasure.setQueryPoint(SbVec3f(0.25f, 0.25f, 0.0f));
	exactMeasure.apply(renderRoot);
	if (exactMeasure.getGeometryPolicy() != SoBRLMeasureAction::FULL_DETAIL ||
	    exactMeasure.getTriangleCount() != 1 ||
	    exactMeasure.getSkippedLodDisplayMeshCount() != 1 ||
	    exactMeasure.getSourceBackedFullDetailRequestCount() != 1) {
	    printf("FAIL: exact measure did not request source-backed full-detail LoD mesh\n");
	    ret = 1;
	} else {
	    const BRLObolSourceMeshRequest &measureSourceRequest =
		exactMeasure.getSourceBackedFullDetailRequest(0);
	    if (check_source_mesh_request(measureSourceRequest,
					  "/lod-submit.bot", "lod-submit.bot", 101)) {
		printf("FAIL: exact measure did not request source-backed full-detail LoD mesh\n");
		ret = 1;
	    }
	    if (!measureSourceRequest.queryBoundsValid ||
		measureSourceRequest.queryBounds.isEmpty() ||
		measureSourceRequest.queryToleranceValid) {
		printf("FAIL: exact measure source-backed request did not carry bounded query metadata\n");
		ret = 1;
	    }
	}
	if (!ret) {
	    BRLObolLodRequest sourceLodRequest;
	    if (!exactMeasure.makeSourceBackedFullDetailLodRequest(0,
		    sourceLodRequest) ||
		check_source_full_detail_lod_request(sourceLodRequest,
			"/lod-submit.bot", "lod-submit.bot") ||
		!has_lod_provider_param(sourceLodRequest,
					"source_query.space") ||
		!has_lod_provider_param(sourceLodRequest,
					"source_query.bounds") ||
		has_lod_provider_param(sourceLodRequest,
				       "source_query.tolerance")) {
		printf("FAIL: exact measure source-backed request did not convert to RT full-detail LoD request\n");
		ret = 1;
	    } else {
		int beforeTriangleCount = exactMeasure.getTriangleCount();
		float beforeArea = exactMeasure.getSurfaceArea();
		BRLObolLodResult sourceResult =
		    source_full_detail_result(sourceLodRequest);
		std::vector<BRLObolLodResult> sourceResults;
		sourceResults.push_back(sourceResult);
		if (exactMeasure.consumeSourceBackedFullDetailResults(
			sourceResults) != 1 ||
		    exactMeasure.getTriangleCount() != beforeTriangleCount + 1 ||
		    exactMeasure.getSurfaceArea() <= beforeArea) {
		    printf("FAIL: exact measure did not consume source-backed full-detail LoD result\n");
		    ret = 1;
		}
		if (!ret) {
		    SoBRLMeasureAction controllerMeasure;
		    controllerMeasure.setQueryPoint(SbVec3f(0.25f, 0.25f,
							    0.0f));
		    controllerMeasure.apply(renderRoot);
		    BRLObolLodRequest controllerMeasureRequest;
		    if (controllerMeasure.getSourceBackedFullDetailRequestCount() != 1 ||
			!controllerMeasure.makeSourceBackedFullDetailLodRequest(0,
				controllerMeasureRequest)) {
			printf("FAIL: controller source-backed exact measure helper did not collect source request\n");
			ret = 1;
		    } else {
			BRLObolLodService controllerMeasureService;
			if (!controllerMeasureService.start(1, TRUE)) {
			    printf("FAIL: controller source-backed exact measure helper service did not start\n");
			    ret = 1;
			} else {
			    if (submit_source_full_detail_task(
				    controllerMeasureService,
				    controllerMeasureRequest)) {
				ret = 1;
			    } else {
				BRLObolViewController measureController(root,
									NULL);
				measureController.setLodService(
				    &controllerMeasureService);
				int submittedCount = -1;
				beforeTriangleCount =
				    controllerMeasure.getTriangleCount();
				beforeArea = controllerMeasure.getSurfaceArea();
				if (measureController.consumeMeasureSourceFullDetail(
					controllerMeasure, 0, &submittedCount) != 1 ||
				    submittedCount != 0 ||
				    controllerMeasure.getTriangleCount() !=
				    beforeTriangleCount + 1 ||
				    controllerMeasure.getSurfaceArea() <=
				    beforeArea) {
				    printf("FAIL: controller source-backed exact measure helper did not consume matching LoD result\n");
				    ret = 1;
				}
			    }
			    controllerMeasureService.stop();
			}
		    }
		}
	    }
	}
    }

    if (!ret) {
	SoBRLMeasureAction limitedMeasureMiss;
	limitedMeasureMiss.setQueryPoint(SbVec3f(5.0f, 5.0f, 0.0f));
	limitedMeasureMiss.setQueryDistanceLimit(0.1f);
	limitedMeasureMiss.apply(renderRoot);
	if (!limitedMeasureMiss.hasQueryDistanceLimit() ||
	    fabsf(limitedMeasureMiss.getQueryDistanceLimit() - 0.1f) >
	    1.0e-6f ||
	    limitedMeasureMiss.hasNearestPrimitive()) {
	    printf("FAIL: exact measure query distance limit did not filter resident nearest primitives\n");
	    ret = 1;
	}
    }

    if (!ret) {
	SoBRLMeasureAction boundedMeasure;
	boundedMeasure.setQueryPoint(SbVec3f(0.25f, 0.25f, 0.0f));
	boundedMeasure.setQueryDistanceLimit(0.5f);
	boundedMeasure.apply(renderRoot);
	if (boundedMeasure.getSourceBackedFullDetailRequestCount() != 1) {
	    printf("FAIL: bounded exact measure did not collect source-backed request\n");
	    ret = 1;
	} else {
	    const BRLObolSourceMeshRequest &boundedMeasureRequest =
		boundedMeasure.getSourceBackedFullDetailRequest(0);
	    if (!boundedMeasureRequest.queryBoundsValid ||
		boundedMeasureRequest.queryBounds.isEmpty() ||
		!boundedMeasureRequest.queryToleranceValid ||
		boundedMeasureRequest.queryTolerance <= 0.0f) {
		printf("FAIL: bounded exact measure source-backed request did not carry explicit query tolerance\n");
		ret = 1;
	    } else {
		BRLObolLodRequest boundedMeasureLodRequest;
		if (!boundedMeasure.makeSourceBackedFullDetailLodRequest(0,
			boundedMeasureLodRequest) ||
		    !has_lod_provider_param(boundedMeasureLodRequest,
					    "source_query.space") ||
		    !has_lod_provider_param(boundedMeasureLodRequest,
					    "source_query.bounds") ||
		    !has_lod_provider_param(boundedMeasureLodRequest,
					    "source_query.tolerance")) {
		    printf("FAIL: bounded exact measure source-backed request did not convert explicit query tolerance\n");
		    ret = 1;
		}
	    }
	}
    }

    if (!ret) {
	SoBRLSnapAction exactSnap;
	exactSnap.setGeometryPolicy(SoBRLSnapAction::FULL_DETAIL);
	exactSnap.setEnabledKinds(SoBRLSnapAction::FACE_NEAREST);
	exactSnap.setQueryPoint(SbVec3f(1.5f, 0.2f, 0.0f));
	exactSnap.setTolerance(0.1f);
	exactSnap.apply(renderRoot);
	if (exactSnap.hasCandidate() ||
	    exactSnap.getSkippedLodDisplayMeshCount() != 1 ||
	    exactSnap.getSourceBackedFullDetailRequestCount() != 1) {
	    printf("FAIL: exact snap did not request source-backed full-detail LoD mesh\n");
	    ret = 1;
	} else {
	    const BRLObolSourceMeshRequest &snapSourceRequest =
		exactSnap.getSourceBackedFullDetailRequest(0);
	    if (check_source_mesh_request(snapSourceRequest, "/lod-submit.bot",
					  "lod-submit.bot", 101)) {
		printf("FAIL: exact snap did not request source-backed full-detail LoD mesh\n");
		ret = 1;
	    } else if (!snapSourceRequest.queryBoundsValid ||
		       snapSourceRequest.queryBounds.isEmpty() ||
		       !snapSourceRequest.queryToleranceValid ||
		       snapSourceRequest.queryTolerance <= 0.0f) {
		printf("FAIL: exact snap source-backed request did not carry bounded query metadata\n");
		ret = 1;
	    }
	}
	if (!ret) {
	    BRLObolLodRequest sourceLodRequest;
	    if (!exactSnap.makeSourceBackedFullDetailLodRequest(0,
		    sourceLodRequest) ||
		check_source_full_detail_lod_request(sourceLodRequest,
			"/lod-submit.bot", "lod-submit.bot") ||
		!has_lod_provider_param(sourceLodRequest,
					"source_query.space") ||
		!has_lod_provider_param(sourceLodRequest,
					"source_query.bounds") ||
		!has_lod_provider_param(sourceLodRequest,
					"source_query.tolerance")) {
		printf("FAIL: exact snap source-backed request did not convert to RT full-detail LoD request\n");
		ret = 1;
	    } else {
		BRLObolLodResult sourceResult =
		    source_full_detail_result(sourceLodRequest);
		std::vector<BRLObolLodResult> sourceResults;
		sourceResults.push_back(sourceResult);
		exactSnap.setQueryPoint(SbVec3f(0.2f, 0.2f, 0.0f));
		if (exactSnap.consumeSourceBackedFullDetailResults(
			sourceResults) != 1 ||
		    !exactSnap.hasCandidate() ||
		    exactSnap.getKind() != SoBRLSnapAction::FACE_NEAREST ||
		    exactSnap.getPrimitiveIndex() != 0 ||
		    strcmp(exactSnap.getPath().getString(),
			   "/lod-submit.bot") != 0) {
		    printf("FAIL: exact snap did not consume source-backed full-detail LoD result\n");
		    ret = 1;
		}
		if (!ret) {
		    SoBRLSnapAction controllerSnap;
		    controllerSnap.setGeometryPolicy(SoBRLSnapAction::FULL_DETAIL);
		    controllerSnap.setEnabledKinds(SoBRLSnapAction::FACE_NEAREST);
		    controllerSnap.setQueryPoint(SbVec3f(1.5f, 0.2f, 0.0f));
		    controllerSnap.setTolerance(0.1f);
		    controllerSnap.apply(renderRoot);
		    BRLObolLodRequest controllerSnapRequest;
		    if (controllerSnap.getSourceBackedFullDetailRequestCount() != 1 ||
			!controllerSnap.makeSourceBackedFullDetailLodRequest(0,
				controllerSnapRequest)) {
			printf("FAIL: controller source-backed exact snap helper did not collect source request\n");
			ret = 1;
		    } else {
			BRLObolLodService controllerSnapService;
			if (!controllerSnapService.start(1, TRUE)) {
			    printf("FAIL: controller source-backed exact snap helper service did not start\n");
			    ret = 1;
			} else {
			    if (submit_source_full_detail_task(
				    controllerSnapService,
				    controllerSnapRequest)) {
				ret = 1;
			    } else {
				BRLObolViewController snapController(root,
								     NULL);
				snapController.setLodService(
				    &controllerSnapService);
				int submittedCount = -1;
				controllerSnap.setQueryPoint(SbVec3f(0.2f,
								     0.2f, 0.0f));
				if (snapController.consumeSnapSourceFullDetail(
					controllerSnap, 0, &submittedCount) != 1 ||
				    submittedCount != 0 ||
				    !controllerSnap.hasCandidate() ||
				    controllerSnap.getKind() !=
				    SoBRLSnapAction::FACE_NEAREST ||
				    controllerSnap.getPrimitiveIndex() != 0 ||
				    strcmp(controllerSnap.getPath().getString(),
					   "/lod-submit.bot") != 0) {
				    printf("FAIL: controller source-backed exact snap helper did not consume matching LoD result\n");
				    ret = 1;
				}
			    }
			    controllerSnapService.stop();
			}
		    }
		}
		BRLObolSourceMeshPickResult sourcePick;
		if (!brlobol_pick_source_full_detail_result(sourcePick,
			exactSnap.getSourceBackedFullDetailRequest(0),
			sourceResult,
			SbVec3f(0.2f, 0.2f, 5.0f),
			SbVec3f(0.0f, 0.0f, -1.0f)) ||
		    !sourcePick.hit ||
		    strcmp(sourcePick.detail.getPath().getString(),
			   "/lod-submit.bot") != 0 ||
		    sourcePick.detail.getPrimitiveKind() !=
		    SoBRLPickDetail::FACE ||
		    sourcePick.detail.getPrimitiveIndex() != 0 ||
		    sourcePick.detail.getFaceVertexIndexA() != 0 ||
		    sourcePick.detail.getFaceVertexIndexB() != 1 ||
		    sourcePick.detail.getFaceVertexIndexC() != 2 ||
		    fabsf(sourcePick.distance - 5.0f) > 1.0e-5f) {
		    printf("FAIL: exact pick did not consume source-backed full-detail LoD result\n");
		    ret = 1;
		}
		SoBRLSourceMeshPickAction sourcePickAction;
		mesh->setPickGeometryPolicy(SoBRLMeshShape::PICK_FULL_DETAIL);
		sourcePickAction.setRay(SbVec3f(0.2f, 0.2f, 5.0f),
					SbVec3f(0.0f, 0.0f, -1.0f));
		sourcePickAction.apply(renderRoot);
		if (sourcePickAction.getVisitedMeshCount() != 2 ||
		    sourcePickAction.getSourceBackedFullDetailRequestCount() != 1) {
		    printf("FAIL: exact pick action did not collect and consume source-backed full-detail LoD result\n");
		    ret = 1;
		} else {
		    BRLObolLodRequest pickLodRequest;
		    BRLObolSourceMeshPickResult actionPick;
		    const BRLObolSourceMeshRequest &pickSourceRequest =
			sourcePickAction.getSourceBackedFullDetailRequest(0);
		    if (!pickSourceRequest.queryRayValid ||
			pickSourceRequest.queryRayDirection.length() <= 0.0f ||
			!sourcePickAction.makeSourceBackedFullDetailLodRequest(0,
				pickLodRequest) ||
			check_source_full_detail_lod_request(pickLodRequest,
				"/lod-submit.bot", "lod-submit.bot") ||
			!has_lod_provider_param(pickLodRequest,
						"source_query.space") ||
			!has_lod_provider_param(pickLodRequest,
						"source_query.ray.origin") ||
			!has_lod_provider_param(pickLodRequest,
						"source_query.ray.direction")) {
			printf("FAIL: exact pick action did not collect and consume source-backed full-detail LoD result\n");
			ret = 1;
		    } else {
			std::vector<BRLObolLodResult> pickSourceResults;
			pickSourceResults.push_back(
			    source_full_detail_result(pickLodRequest));
			if (sourcePickAction.consumeSourceBackedFullDetailResults(
				actionPick, pickSourceResults) != 1 ||
			    !actionPick.hit ||
			    strcmp(actionPick.detail.getPath().getString(),
				   "/lod-submit.bot") != 0 ||
			    actionPick.detail.getPrimitiveIndex() != 0) {
			    printf("FAIL: exact pick action did not collect and consume source-backed full-detail LoD result\n");
			    ret = 1;
			}
			if (!ret) {
			    BRLObolLodService controllerPickService;
			    if (!controllerPickService.start(1, TRUE)) {
				printf("FAIL: controller source-backed exact pick helper service did not start\n");
				ret = 1;
			    } else {
				BRLObolLodTask controllerPickTask;
				controllerPickTask.generation =
				    controllerPickService.beginGeneration();
				controllerPickTask.request = pickLodRequest;
				controllerPickTask.realize = source_full_detail_task;
				if (controllerPickService.submit(controllerPickTask) == 0 ||
				    wait_for_service(controllerPickService)) {
				    ret = 1;
				} else {
				    BRLObolViewController pickController(root, NULL);
				    pickController.setLodService(&controllerPickService);
				    BRLObolSourceMeshPickResult controllerPick;
				    int submittedCount = -1;
				    if (pickController.pickSourceMeshExactRay(
					    controllerPick,
					    SbVec3f(0.2f, 0.2f, 5.0f),
					    SbVec3f(0.0f, 0.0f, -1.0f),
					    0, &submittedCount) != 1 ||
					submittedCount != 0 ||
					!controllerPick.hit ||
					strcmp(controllerPick.detail.getPath().getString(),
					       "/lod-submit.bot") != 0 ||
					controllerPick.detail.getPrimitiveIndex() != 0) {
					printf("FAIL: controller source-backed exact pick helper did not consume matching LoD result\n");
					ret = 1;
				    }
				}
				controllerPickService.stop();
			    }
			}
		    }
		}
		if (!ret) {
		    BRLObolSourceMeshRequest mappedSourceRequest =
			exactSnap.getSourceBackedFullDetailRequest(0);
		    mappedSourceRequest.faceCount = 0;

		    BRLObolLodResult mappedSourceResult =
			source_full_detail_result(sourceLodRequest);
		    mappedSourceResult.mesh.faceIndex.push_back(7);
		    mappedSourceResult.mesh.vertexIndex.push_back(10);
		    mappedSourceResult.mesh.vertexIndex.push_back(11);
		    mappedSourceResult.mesh.vertexIndex.push_back(12);

		    SoBRLExportAction mappedExport;
		    if (!mappedExport.appendSourceBackedFullDetailResult(
			    mappedSourceRequest, mappedSourceResult) ||
			mappedExport.getTriangleCount() != 1 ||
			mappedExport.getTriangle(0).primitiveIndex != 7 ||
			mappedExport.getTriangle(0).vertexIndexA != 10 ||
			mappedExport.getTriangle(0).vertexIndexB != 11 ||
			mappedExport.getTriangle(0).vertexIndexC != 12) {
			printf("FAIL: exact export did not preserve source face and vertex index mapping\n");
			ret = 1;
		    }

		    SoBRLMeasureAction mappedMeasure;
		    mappedMeasure.setQueryPoint(SbVec3f(0.2f, 0.2f, 0.0f));
		    if (!mappedMeasure.consumeSourceBackedFullDetailResult(
			    mappedSourceRequest, mappedSourceResult) ||
			!mappedMeasure.hasNearestPrimitive() ||
			mappedMeasure.getNearestPrimitiveIndex() != 7 ||
			mappedMeasure.getNearestFaceVertexIndexA() != 10 ||
			mappedMeasure.getNearestFaceVertexIndexB() != 11 ||
			mappedMeasure.getNearestFaceVertexIndexC() != 12 ||
			mappedMeasure.getNearestFaceVertexIndex() != 10) {
			printf("FAIL: exact measure did not preserve source face and vertex index mapping\n");
			ret = 1;
		    }

		    BRLObolSourceMeshRequest compactMeasureRequest =
			mappedSourceRequest;
		    compactMeasureRequest.faceCount = 2;
		    compactMeasureRequest.pointCount = 6;
		    compactMeasureRequest.queryBoundsValid = 1;
		    compactMeasureRequest.queryBounds = SbBox3f(
							    SbVec3f(-0.1f, -0.1f, -0.1f),
							    SbVec3f(0.2f, 0.2f, 0.2f));
		    compactMeasureRequest.queryToleranceValid = 1;
		    compactMeasureRequest.queryTolerance = 0.5f;
		    BRLObolLodResult compactMeasureResult =
			source_full_detail_result(sourceLodRequest);
		    compactMeasureResult.mesh.faceIndex.push_back(1);
		    compactMeasureResult.mesh.vertexIndex.push_back(30);
		    compactMeasureResult.mesh.vertexIndex.push_back(31);
		    compactMeasureResult.mesh.vertexIndex.push_back(32);
		    SoBRLMeasureAction compactMeasure;
		    compactMeasure.setQueryPoint(SbVec3f(0.0f, 0.0f, 0.0f));
		    compactMeasure.setQueryDistanceLimit(0.5f);
		    if (!compactMeasure.consumeSourceBackedFullDetailResult(
			    compactMeasureRequest, compactMeasureResult) ||
			!compactMeasure.hasNearestPrimitive() ||
			compactMeasure.getNearestPrimitiveIndex() != 1 ||
			compactMeasure.getNearestFaceVertexIndexA() != 30) {
			printf("FAIL: exact measure did not accept compact bounded source subset with identity mapping\n");
			ret = 1;
		    }
		    BRLObolLodResult compactMeasureNoFaceMap =
			compactMeasureResult;
		    compactMeasureNoFaceMap.mesh.faceIndex.clear();
		    SoBRLMeasureAction rejectedCompactMeasureFaceMap;
		    if (rejectedCompactMeasureFaceMap.consumeSourceBackedFullDetailResult(
			    compactMeasureRequest, compactMeasureNoFaceMap)) {
			printf("FAIL: exact measure accepted compact bounded source subset without face index mapping\n");
			ret = 1;
		    }
		    BRLObolLodResult compactMeasureNoVertexMap =
			compactMeasureResult;
		    compactMeasureNoVertexMap.mesh.vertexIndex.clear();
		    SoBRLMeasureAction rejectedCompactMeasureVertexMap;
		    if (rejectedCompactMeasureVertexMap.consumeSourceBackedFullDetailResult(
			    compactMeasureRequest, compactMeasureNoVertexMap)) {
			printf("FAIL: exact measure accepted compact bounded source subset without vertex index mapping\n");
			ret = 1;
		    }

		    SoBRLSnapAction mappedSnap;
		    mappedSnap.setEnabledKinds(SoBRLSnapAction::FACE_NEAREST);
		    mappedSnap.setQueryPoint(SbVec3f(0.2f, 0.2f, 0.0f));
		    mappedSnap.setTolerance(0.5f);
		    if (!mappedSnap.consumeSourceBackedFullDetailResult(
			    mappedSourceRequest, mappedSourceResult) ||
			!mappedSnap.hasCandidate() ||
			mappedSnap.getPrimitiveIndex() != 7) {
			printf("FAIL: exact snap did not preserve source face index mapping\n");
			ret = 1;
		    }
		    SoBRLSnapAction mappedVertexSnap;
		    mappedVertexSnap.setEnabledKinds(SoBRLSnapAction::VERTEX);
		    mappedVertexSnap.setQueryPoint(SbVec3f(0.0f, 0.0f, 0.0f));
		    mappedVertexSnap.setTolerance(0.5f);
		    if (!mappedVertexSnap.consumeSourceBackedFullDetailResult(
			    mappedSourceRequest, mappedSourceResult) ||
			!mappedVertexSnap.hasCandidate() ||
			mappedVertexSnap.getPrimitiveIndex() != 7 ||
			mappedVertexSnap.getVertexIndex() != 10) {
			printf("FAIL: exact snap did not preserve source vertex index mapping\n");
			ret = 1;
		    }
		    BRLObolSourceMeshRequest compactSnapRequest =
			mappedSourceRequest;
		    compactSnapRequest.faceCount = 2;
		    compactSnapRequest.pointCount = 6;
		    compactSnapRequest.queryBoundsValid = 1;
		    compactSnapRequest.queryBounds = SbBox3f(
							 SbVec3f(-0.1f, -0.1f, -0.1f),
							 SbVec3f(0.2f, 0.2f, 0.2f));
		    compactSnapRequest.queryToleranceValid = 1;
		    compactSnapRequest.queryTolerance = 0.5f;
		    BRLObolLodResult compactSnapResult =
			source_full_detail_result(sourceLodRequest);
		    compactSnapResult.mesh.faceIndex.push_back(1);
		    compactSnapResult.mesh.vertexIndex.push_back(30);
		    compactSnapResult.mesh.vertexIndex.push_back(31);
		    compactSnapResult.mesh.vertexIndex.push_back(32);
		    SoBRLSnapAction compactSnap;
		    compactSnap.setEnabledKinds(SoBRLSnapAction::VERTEX);
		    compactSnap.setQueryPoint(SbVec3f(0.0f, 0.0f, 0.0f));
		    compactSnap.setTolerance(0.5f);
		    if (!compactSnap.consumeSourceBackedFullDetailResult(
			    compactSnapRequest, compactSnapResult) ||
			!compactSnap.hasCandidate() ||
			compactSnap.getPrimitiveIndex() != 1 ||
			compactSnap.getVertexIndex() != 30) {
			printf("FAIL: exact snap did not accept compact source vertex subset with vertex index mapping\n");
			ret = 1;
		    }
		    BRLObolLodResult compactSnapNoVertexMap =
			compactSnapResult;
		    compactSnapNoVertexMap.mesh.vertexIndex.clear();
		    SoBRLSnapAction rejectedCompactSnap;
		    if (rejectedCompactSnap.consumeSourceBackedFullDetailResult(
			    compactSnapRequest, compactSnapNoVertexMap)) {
			printf("FAIL: exact snap accepted compact source vertex subset without vertex index mapping\n");
			ret = 1;
		    }

		    BRLObolSourceMeshPickResult mappedPick;
		    if (!brlobol_pick_source_full_detail_result(mappedPick,
			    mappedSourceRequest, mappedSourceResult,
			    SbVec3f(0.2f, 0.2f, 5.0f),
			    SbVec3f(0.0f, 0.0f, -1.0f)) ||
			!mappedPick.hit ||
			mappedPick.detail.getPrimitiveIndex() != 7 ||
			mappedPick.detail.getFaceVertexIndexA() != 10 ||
			mappedPick.detail.getFaceVertexIndexB() != 11 ||
			mappedPick.detail.getFaceVertexIndexC() != 12 ||
			mappedPick.detail.getNearestFaceVertexIndex() != 10) {
			printf("FAIL: exact pick did not preserve source face and vertex index mapping\n");
			ret = 1;
		    }

		    BRLObolSourceMeshRequest rayScopedPickRequest =
			mappedSourceRequest;
		    rayScopedPickRequest.faceCount = 4;
		    rayScopedPickRequest.pointCount = 6;
		    rayScopedPickRequest.queryRayValid = 1;
		    rayScopedPickRequest.queryRayDirection =
			SbVec3f(0.0f, 0.0f, -1.0f);
		    BRLObolLodResult rayScopedPickResult =
			source_full_detail_result(sourceLodRequest);
		    rayScopedPickResult.mesh.faceIndex.push_back(2);
		    rayScopedPickResult.mesh.vertexIndex.push_back(20);
		    rayScopedPickResult.mesh.vertexIndex.push_back(21);
		    rayScopedPickResult.mesh.vertexIndex.push_back(22);
		    BRLObolSourceMeshPickResult rayScopedPick;
		    if (!brlobol_pick_source_full_detail_result(rayScopedPick,
			    rayScopedPickRequest, rayScopedPickResult,
			    SbVec3f(0.2f, 0.2f, 5.0f),
			    SbVec3f(0.0f, 0.0f, -1.0f)) ||
			!rayScopedPick.hit ||
			rayScopedPick.detail.getPrimitiveIndex() != 2 ||
			rayScopedPick.detail.getFaceVertexIndexA() != 20 ||
			rayScopedPick.detail.getFaceVertexIndexB() != 21 ||
			rayScopedPick.detail.getFaceVertexIndexC() != 22) {
			printf("FAIL: exact pick did not accept ray-scoped source face and vertex subset\n");
			ret = 1;
		    }
		    BRLObolLodResult rayScopedPickNoVertexMap =
			rayScopedPickResult;
		    rayScopedPickNoVertexMap.mesh.vertexIndex.clear();
		    if (brlobol_pick_source_full_detail_result(rayScopedPick,
			    rayScopedPickRequest, rayScopedPickNoVertexMap,
			    SbVec3f(0.2f, 0.2f, 5.0f),
			    SbVec3f(0.0f, 0.0f, -1.0f))) {
			printf("FAIL: exact pick accepted compact ray-scoped source point subset without vertex index mapping\n");
			ret = 1;
		    }
		    rayScopedPickRequest.queryRayValid = 0;
		    if (brlobol_pick_source_full_detail_result(rayScopedPick,
			    rayScopedPickRequest, rayScopedPickResult,
			    SbVec3f(0.2f, 0.2f, 5.0f),
			    SbVec3f(0.0f, 0.0f, -1.0f))) {
			printf("FAIL: exact pick accepted non-query source face subset\n");
			ret = 1;
		    }
		}
	    }
	}
    }

    if (!ret) {
	size_t displayBytes = viewState.estimateDisplayMeshBytes();
	unsigned int evictedMeshCount = 0;
	size_t freedBytes = viewState.evictDisplayMeshes(&evictedMeshCount);
	if (displayBytes == 0 ||
	    freedBytes != displayBytes ||
	    evictedMeshCount != 1 ||
	    viewState.estimateDisplayMeshBytes() != 0 ||
	    viewState.findMesh(mesh) ||
	    mesh->isLodDisplayActive() ||
	    mesh->lodAvailable.getValue() ||
	    mesh->getTriangleCount() != 1) {
	    printf("FAIL: view-local display eviction did not clear view payload without mutating mesh\n");
	    ret = 1;
	}
    }

    if (!ret) {
	SoBRLExportAction evictedDisplayExactExport;
	evictedDisplayExactExport.apply(renderRoot);
	if (evictedDisplayExactExport.getGeometryPolicy() !=
	    SoBRLExportAction::FULL_DETAIL ||
	    evictedDisplayExactExport.getTriangleCount() != 1 ||
	    evictedDisplayExactExport.getSkippedLodDisplayMeshCount() != 0 ||
	    evictedDisplayExactExport.getSourceBackedFullDetailRequestCount() != 1 ||
	    check_source_mesh_request(
		evictedDisplayExactExport.getSourceBackedFullDetailRequest(0),
		"/lod-submit.bot", "lod-submit.bot", 101)) {
	    printf("FAIL: evicted active display mesh did not keep source-backed exact export\n");
	    ret = 1;
	}
    }

    if (!ret) {
	SoBRLMeasureAction evictedDisplayExactMeasure;
	evictedDisplayExactMeasure.apply(renderRoot);
	if (evictedDisplayExactMeasure.getGeometryPolicy() !=
	    SoBRLMeasureAction::FULL_DETAIL ||
	    evictedDisplayExactMeasure.getTriangleCount() != 1 ||
	    evictedDisplayExactMeasure.getSkippedLodDisplayMeshCount() != 0 ||
	    evictedDisplayExactMeasure.getSourceBackedFullDetailRequestCount() != 1 ||
	    check_source_mesh_request(
		evictedDisplayExactMeasure.getSourceBackedFullDetailRequest(0),
		"/lod-submit.bot", "lod-submit.bot", 101)) {
	    printf("FAIL: evicted active display mesh did not keep source-backed exact measure\n");
	    ret = 1;
	}
    }

    if (!ret) {
	SoBRLSnapAction evictedDisplayExactSnap;
	evictedDisplayExactSnap.setGeometryPolicy(SoBRLSnapAction::FULL_DETAIL);
	evictedDisplayExactSnap.setEnabledKinds(SoBRLSnapAction::FACE_NEAREST);
	evictedDisplayExactSnap.setQueryPoint(SbVec3f(1.5f, 0.2f, 0.0f));
	evictedDisplayExactSnap.setTolerance(0.1f);
	evictedDisplayExactSnap.apply(renderRoot);
	if (evictedDisplayExactSnap.hasCandidate() ||
	    evictedDisplayExactSnap.getSkippedLodDisplayMeshCount() != 0 ||
	    evictedDisplayExactSnap.getSourceBackedFullDetailRequestCount() != 1 ||
	    check_source_mesh_request(
		evictedDisplayExactSnap.getSourceBackedFullDetailRequest(0),
		"/lod-submit.bot", "lod-submit.bot", 101)) {
	    printf("FAIL: evicted active display mesh did not keep source-backed exact snap\n");
	    ret = 1;
	}
    }

    if (!ret) {
	SoBRLSourceMeshPickAction evictedDisplayPick;
	BRLObolLodRequest evictedDisplayPickRequest;
	mesh->setPickGeometryPolicy(SoBRLMeshShape::PICK_FULL_DETAIL);
	evictedDisplayPick.setRay(SbVec3f(0.2f, 0.2f, 5.0f),
				  SbVec3f(0.0f, 0.0f, -1.0f));
	evictedDisplayPick.apply(renderRoot);
	if (evictedDisplayPick.getSourceBackedFullDetailRequestCount() != 1 ||
	    !evictedDisplayPick.makeSourceBackedFullDetailLodRequest(0,
		    evictedDisplayPickRequest) ||
	    check_source_full_detail_lod_request(evictedDisplayPickRequest,
		    "/lod-submit.bot", "lod-submit.bot")) {
	    printf("FAIL: evicted active display mesh did not keep source-backed exact pick\n");
	    ret = 1;
	}
	mesh->setPickGeometryPolicy(SoBRLMeshShape::PICK_DISPLAY_LEVEL);
    }

    if (!ret) {
	int forcedLevel = expectedViewLevel;
	SoBRLMeshLodSubmitAction forcedSubmit;
	forcedSubmit.setService(&service);
	forcedSubmit.setViewLodState(&viewState);
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
		forcedUpdate.setViewLodState(&viewState);
		forcedUpdate.setResults(forcedResults);
		forcedUpdate.apply(root);
		const BRLObolViewLodState::MeshPayload *forcedPayload =
		    viewState.findMesh(mesh);
		if (forcedUpdate.getAppliedResultCount() != 1 ||
		    !forcedPayload ||
		    forcedPayload->activeLevel != forcedLevel) {
		    printf("FAIL: LoD submit action forced-level result was not applied\n");
		    ret = 1;
		}
	    }
	}
    }

    service.stop();
    renderRoot->unref();
    root->unref();
    brlobol_mesh_lod_cache_clear_database(dbip);
    db_close(dbip);
    bu_file_delete(dbpath);
    bu_dirclear(cache_dir);
    return ret;
}

static int
test_view_controller_source_backed_multi_source_exact_submit(void)
{
    char wrong_dbpath[MAXPATHLEN] = {0};
    char right_dbpath[MAXPATHLEN] = {0};
    struct db_i *wrong_dbip = NULL;
    struct db_i *right_dbip = NULL;

    if (make_rt_pick_test_db(wrong_dbpath, sizeof(wrong_dbpath),
			     &wrong_dbip))
	return 1;
    if (make_submit_test_db(right_dbpath, sizeof(right_dbpath),
			    &right_dbip)) {
	db_close(wrong_dbip);
	bu_file_delete(wrong_dbpath);
	return 1;
    }

    SoSeparator *root = new SoSeparator;
    root->ref();

    SoBRLDatabaseSource *wrongSource = new SoBRLDatabaseSource;
    wrongSource->setDatabase(wrong_dbip);
    wrongSource->path = "implicit.r";
    wrongSource->sourceRevision = 111;

    SoBRLDatabaseSource *rightSource = new SoBRLDatabaseSource;
    rightSource->setDatabase(right_dbip);
    rightSource->path = "lod-submit.bot";
    rightSource->sourceRevision = 222;

    SoBRLLodMeshShape *mesh = make_lod_mesh("/lod-submit.bot",
					    "lod-submit.bot");
    mesh->sourceId = 222;
    BRLObolLodRequest displayRequest =
	make_request("/lod-submit.bot", "lod-submit.bot");
    BRLObolLodResult displayResult = mesh_payload_result(displayRequest);
    if (!mesh->applyStagedLodResult(displayResult, &displayRequest) ||
	!mesh->isLodDisplayActive()) {
	printf("FAIL: controller multi-source exact submit fixture did not stage LoD mesh\n");
	root->unref();
	brlobol_mesh_lod_cache_clear_database(right_dbip);
	db_close(right_dbip);
	db_close(wrong_dbip);
	bu_file_delete(right_dbpath);
	bu_file_delete(wrong_dbpath);
	return 1;
    }

    rightSource->addChild(mesh);
    root->addChild(wrongSource);
    root->addChild(rightSource);

    BRLObolLodService service;
    if (!service.start(1, TRUE)) {
	printf("FAIL: controller multi-source source-backed exact service did not start\n");
	root->unref();
	brlobol_mesh_lod_cache_clear_database(right_dbip);
	db_close(right_dbip);
	db_close(wrong_dbip);
	bu_file_delete(right_dbpath);
	bu_file_delete(wrong_dbpath);
	return 1;
    }

    int ret = 0;
    BRLObolViewController controller(root, NULL);
    controller.setLodService(&service);

    SoBRLExportAction exportAction;
    exportAction.apply(root);
    if (exportAction.getSourceBackedFullDetailRequestCount() != 1) {
	printf("FAIL: controller multi-source source-backed exact submit did not collect request\n");
	ret = 1;
    }

    int submittedCount = -1;
    if (!ret &&
	(controller.consumeExportSourceFullDetail(exportAction,
		service.beginGeneration(), &submittedCount) != 0 ||
	 submittedCount != 1)) {
	printf("FAIL: controller multi-source source-backed exact submit did not submit one request\n");
	ret = 1;
    }

    if (!ret && wait_for_service(service))
	ret = 1;

    if (!ret) {
	std::vector<BRLObolLodResult> submittedResults;
	service.drainResults(submittedResults);
	const char *rightDbId = right_dbip->dbi_filename ?
				right_dbip->dbi_filename : "";
	if (submittedResults.size() != 1 ||
	    submittedResults[0].providerStatus !=
	    BRLOBOL_LOD_PROVIDER_STALE ||
	    strcmp(submittedResults[0].request.databaseId.getString(),
		   rightDbId) != 0 ||
	    submittedResults[0].request.sourceRevision != 222 ||
	    strcmp(submittedResults[0].request.objectPath.getString(),
		   "/lod-submit.bot") != 0) {
	    printf("FAIL: controller multi-source source-backed exact submit did not use matching database source\n");
	    ret = 1;
	} else {
	    BRLObolLodTask task;
	    task.generation = service.beginGeneration();
	    task.request = submittedResults[0].request;
	    task.realize = source_full_detail_task;
	    if (service.submit(task) == 0 || wait_for_service(service)) {
		ret = 1;
	    } else {
		const int beforeTriangleCount =
		    exportAction.getTriangleCount();
		submittedCount = -1;
		if (controller.consumeExportSourceFullDetail(exportAction,
			service.beginGeneration(), &submittedCount) != 1 ||
		    submittedCount != 0 ||
		    exportAction.getTriangleCount() !=
		    beforeTriangleCount + 1 ||
		    exportAction.getTriangle(beforeTriangleCount).sourceId !=
		    222) {
		    printf("FAIL: controller multi-source source-backed exact submit did not consume matching database-scoped result\n");
		    ret = 1;
		}
	    }
	}
    }

    service.stop();
    root->unref();
    brlobol_mesh_lod_cache_clear_database(right_dbip);
    db_close(right_dbip);
    db_close(wrong_dbip);
    bu_file_delete(right_dbpath);
    bu_file_delete(wrong_dbpath);
    return ret;
}

static int
test_view_controller_source_backed_partial_ready_submit(void)
{
    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;

    if (make_submit_test_db(dbpath, sizeof(dbpath), &dbip))
	return 1;

    SoSeparator *root = new SoSeparator;
    root->ref();

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = "lod-submit.bot";
    source->sourceRevision = 333;

    SoBRLLodMeshShape *readyMesh = make_lod_mesh("/lod-submit.bot",
				   "lod-submit.bot");
    readyMesh->sourceId = 333;
    BRLObolLodRequest readyDisplayRequest =
	make_request("/lod-submit.bot", "lod-submit.bot");
    BRLObolLodResult readyDisplayResult =
	mesh_payload_result(readyDisplayRequest);
    if (!readyMesh->applyStagedLodResult(readyDisplayResult,
					 &readyDisplayRequest) ||
	!readyMesh->isLodDisplayActive()) {
	printf("FAIL: controller partial-ready exact fixture did not stage ready LoD mesh\n");
	root->unref();
	brlobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    SoBRLLodMeshShape *missingMesh = make_lod_mesh(
					 "/lod-submit-missing.bot", "lod-submit-missing.bot");
    missingMesh->sourceId = 333;
    BRLObolLodRequest missingDisplayRequest =
	make_request("/lod-submit-missing.bot", "lod-submit-missing.bot");
    BRLObolLodResult missingDisplayResult =
	mesh_payload_result(missingDisplayRequest);
    if (!missingMesh->applyStagedLodResult(missingDisplayResult,
					   &missingDisplayRequest) ||
	!missingMesh->isLodDisplayActive()) {
	printf("FAIL: controller partial-ready exact fixture did not stage missing LoD mesh\n");
	root->unref();
	brlobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    source->addChild(readyMesh);
    source->addChild(missingMesh);
    root->addChild(source);

    BRLObolLodService service;
    if (!service.start(1, TRUE)) {
	printf("FAIL: controller partial-ready source-backed exact service did not start\n");
	root->unref();
	brlobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    int ret = 0;
    SoBRLExportAction exportAction;
    exportAction.apply(root);
    if (exportAction.getSourceBackedFullDetailRequestCount() != 2) {
	printf("FAIL: controller partial-ready source-backed exact test did not collect two requests\n");
	ret = 1;
    }

    BRLObolLodRequest readySourceRequest;
    if (!ret && !exportAction.makeSourceBackedFullDetailLodRequest(0,
	    readySourceRequest)) {
	printf("FAIL: controller partial-ready source-backed exact test did not convert ready request\n");
	ret = 1;
    }

    if (!ret) {
	BRLObolLodTask readyTask;
	readyTask.generation = service.beginGeneration();
	readyTask.request = readySourceRequest;
	readyTask.realize = source_full_detail_task;
	if (service.submit(readyTask) == 0 || wait_for_service(service))
	    ret = 1;
    }

    if (!ret) {
	BRLObolViewController controller(root, NULL);
	controller.setLodService(&service);
	const int beforeTriangleCount = exportAction.getTriangleCount();
	int submittedCount = -1;
	if (controller.consumeExportSourceFullDetail(exportAction,
		service.beginGeneration(), &submittedCount) != 1 ||
	    submittedCount != 1 ||
	    exportAction.getTriangleCount() != beforeTriangleCount + 1 ||
	    exportAction.getTriangle(beforeTriangleCount).sourceId != 333) {
	    printf("FAIL: controller partial-ready exact helper did not consume ready result and submit missing request\n");
	    ret = 1;
	}
    }

    service.stop();

    if (!ret) {
	readyMesh->setPickGeometryPolicy(SoBRLMeshShape::PICK_FULL_DETAIL);
	missingMesh->setPickGeometryPolicy(SoBRLMeshShape::PICK_FULL_DETAIL);

	BRLObolLodService pickService;
	if (!pickService.start(1, TRUE)) {
	    printf("FAIL: controller partial-ready exact pick service did not start\n");
	    ret = 1;
	} else {
	    const SbVec3f rayOrigin(0.2f, 0.2f, 5.0f);
	    const SbVec3f rayDirection(0.0f, 0.0f, -1.0f);
	    SoBRLSourceMeshPickAction requestAction;
	    requestAction.setRay(rayOrigin, rayDirection);
	    requestAction.apply(root);

	    BRLObolLodRequest readyPickRequest;
	    if (requestAction.getSourceBackedFullDetailRequestCount() != 2 ||
		!requestAction.makeSourceBackedFullDetailLodRequest(0,
			readyPickRequest)) {
		printf("FAIL: controller partial-ready exact pick did not collect two source requests\n");
		ret = 1;
	    } else {
		BRLObolLodTask readyPickTask;
		readyPickTask.generation = pickService.beginGeneration();
		readyPickTask.request = readyPickRequest;
		readyPickTask.realize = source_full_detail_task;
		if (pickService.submit(readyPickTask) == 0 ||
		    wait_for_service(pickService)) {
		    ret = 1;
		} else {
		    BRLObolViewController pickController(root, NULL);
		    pickController.setLodService(&pickService);
		    BRLObolSourceMeshPickResult pick;
		    int submittedCount = -1;
		    if (pickController.pickSourceMeshExactRay(pick,
			    rayOrigin, rayDirection,
			    pickService.beginGeneration(),
			    &submittedCount) != 1 ||
			submittedCount != 1 ||
			!pick.hit ||
			strcmp(pick.detail.getPath().getString(),
			       "/lod-submit.bot") != 0 ||
			pick.detail.getPrimitiveIndex() != 0) {
			printf("FAIL: controller partial-ready exact pick helper did not consume ready result and submit missing request\n");
			ret = 1;
		    }
		}
	    }
	    pickService.stop();
	}
    }

    root->unref();
    brlobol_mesh_lod_cache_clear_database(dbip);
    db_close(dbip);
    bu_file_delete(dbpath);
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
	brlobol_mesh_lod_cache_clear_database(dbip);
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
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	if (controller.submitLodRequests(NULL, service.beginGeneration()) != 3 ||
	    controller.getLastLodVisitedMeshCount() != 1 ||
	    controller.getLastLodSubmittedTaskCount() != 3 ||
	    controller.getLastLodSkippedMeshCount() != 0 ||
	    controller.getLastLodDiagnostics().getLength() != 0 ||
	    !controller.isRenderRequested() ||
	    strcmp(controller.getRenderReason().getString(),
		   "lod-submit") != 0) {
	    printf("FAIL: LoD view controller did not submit expected task\n");
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	if (wait_for_service(service)) {
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	if (!controller.hasPendingLodResults() || !controller.isRenderRequested()) {
	    printf("FAIL: LoD view controller did not observe result-ready wakeup\n");
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	controller.clearRenderRequest();
	if (!controller.hasPendingLodResults() ||
	    controller.processPendingLodResults() != 3 ||
	    controller.getLastLodResultCount() != 3 ||
	    controller.getLastLodMatchedResultCount() != 3 ||
	    controller.getLastLodAppliedResultCount() != 3 ||
	    controller.getLastLodRejectedResultCount() != 0 ||
	    controller.getLastLodUnmatchedResultCount() != 0 ||
	    controller.getLastLodDiagnostics().getLength() != 0 ||
	    !controller.getViewLodState()->findMesh(mesh) ||
	    !controller.getViewLodState()->findProxy(mesh) ||
	    !controller.isRenderRequested() ||
	    strcmp(controller.getRenderReason().getString(),
		   "lod-result") != 0 ||
	    controller.hasPendingLodResults()) {
	    printf("FAIL: LoD view controller did not apply service result\n");
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
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
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	controller.clearRenderRequest();
	if (controller.submitLodRequestsIfNeeded() != 0 ||
	    controller.getLastLodSubmittedTaskCount() != 0 ||
	    controller.getLastLodSkippedMeshCount() != 1 ||
	    strstr(controller.getLastLodDiagnostics().getString(),
		   "current LoD request is already resident") == NULL ||
	    controller.isRenderRequested()) {
	    printf("FAIL: LoD view controller did not skip resident changed-scene LoD request\n");
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (controller.submitLodRequestsIfNeeded() != 0) {
	    printf("FAIL: LoD view controller duplicated unchanged auto-submit\n");
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (service.queuedResultCountForDiagnostics() != 0) {
	    printf("FAIL: LoD view controller queued duplicate resident request result\n");
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
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
	if (controller.submitLodRequestsIfNeeded() != 0 ||
	    controller.getLastLodSubmittedTaskCount() != 0 ||
	    controller.getLastLodSkippedMeshCount() != 1 ||
	    strstr(controller.getLastLodDiagnostics().getString(),
		   "current LoD request is already resident") == NULL ||
	    controller.isRenderRequested()) {
	    printf("FAIL: LoD view controller did not skip resident threshold policy request\n");
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (service.queuedResultCountForDiagnostics() != 0) {
	    printf("FAIL: LoD view controller queued duplicate threshold request result\n");
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	uint64_t previousViewRevision = controller.getLodViewRevision();
	camera->height = 90.0f;
	controller.clearRenderRequest();
	if (controller.submitLodRequestsIfNeeded() != 3 ||
	    controller.getLastLodSubmittedTaskCount() != 3 ||
	    controller.getLodViewRevision() == previousViewRevision ||
	    !controller.isRenderRequested() ||
	    strcmp(controller.getRenderReason().getString(),
		   "lod-submit") != 0) {
	    printf("FAIL: LoD view controller did not auto-submit camera field change\n");
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (wait_for_service(service)) {
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	controller.clearRenderRequest();
	if (controller.processPendingLodResults() != 3 ||
	    controller.getLastLodAppliedResultCount() != 3 ||
	    strcmp(controller.getRenderReason().getString(),
		   "lod-result") != 0) {
	    printf("FAIL: LoD view controller did not apply camera field result\n");
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	const BRLObolViewLodState::MeshPayload *activePayload =
	    controller.getViewLodState()->findMesh(mesh);
	int forcedLevel = activePayload ? activePayload->activeLevel : 0;
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
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	controller.clearRenderRequest();
	if (controller.submitLodRequestsIfNeeded() != 3 ||
	    controller.getLastLodSubmittedTaskCount() != 3 ||
	    !controller.isRenderRequested() ||
	    strcmp(controller.getRenderReason().getString(),
		   "lod-submit") != 0) {
	    printf("FAIL: LoD view controller did not auto-submit forced-level policy\n");
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (wait_for_service(service)) {
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	controller.clearRenderRequest();
	if (controller.processPendingLodResults() != 3 ||
	    controller.getLastLodAppliedResultCount() != 3 ||
	    !controller.getViewLodState()->findMesh(mesh) ||
	    controller.getViewLodState()->findMesh(mesh)->activeLevel !=
	    forcedLevel ||
	    strcmp(controller.getRenderReason().getString(),
		   "lod-result") != 0) {
	    printf("FAIL: LoD view controller did not apply forced-level result\n");
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
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
	    brlobol_mesh_lod_cache_clear_database(dbip);
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
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	controller.clearRenderRequest();
	if (controller.submitLodRequestsIfNeeded() != 3 ||
	    controller.getLastLodSubmittedTaskCount() != 3) {
	    printf("FAIL: LoD view controller did not auto-submit view change\n");
	    service.stop();
	    root->unref();
	    brlobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
    }

    int ret = 0;
    if (mesh->lodAvailable.getValue() ||
	mesh->isLodDisplayActive() ||
	mesh->point.getNum() != 3 ||
	mesh->coordIndex.getNum() != 3 ||
	mesh->getTriangleCount() != 1) {
	printf("FAIL: LoD view-controller result mutated shared mesh payload\n");
	ret = 1;
    }

    service.stop();
    root->unref();
    brlobol_mesh_lod_cache_clear_database(dbip);
    db_close(dbip);
    bu_file_delete(dbpath);
    bu_dirclear(cache_dir);
    return ret;
}

static int
test_view_controller_shared_lod_is_view_local(void)
{
    SoSeparator *shared = new SoSeparator;
    shared->ref();
    SoBRLLodMeshShape *mesh = make_lod_mesh("/shared/lod.bot",
					    "lod.bot");
    shared->addChild(mesh);

    int ret = 0;
    {
	BRLObolViewController viewA(shared, NULL);
	BRLObolViewController viewB(shared, NULL);

	BRLObolLodRequest request = make_request("/shared/lod.bot",
				    "lod.bot");
	BRLObolLodResult resultA =
	    mesh_payload_variant_result(request, 1, 2);
	BRLObolLodResult resultB =
	    mesh_payload_variant_result(request, 3, 3);

	SoBRLLodUpdateAction updateA;
	updateA.setViewLodState(viewA.getViewLodState());
	updateA.addResult(resultA);
	updateA.apply(viewA.getRenderSceneRoot());
	SoBRLLodUpdateAction updateB;
	updateB.setViewLodState(viewB.getViewLodState());
	updateB.addResult(resultB);
	updateB.apply(viewB.getRenderSceneRoot());

	if (updateA.getAppliedResultCount() != 1 ||
	    updateB.getAppliedResultCount() != 1 ||
	    !viewA.getViewLodState()->findMesh(mesh) ||
	    !viewB.getViewLodState()->findMesh(mesh)) {
	    printf("FAIL: shared LoD view-local update did not bind payloads\n");
	    ret = 1;
	}

	if (!ret) {
	    SoBRLExportAction exportA;
	    exportA.setGeometryPolicy(SoBRLExportAction::DISPLAY_LEVEL);
	    exportA.apply(viewA.getRenderSceneRoot());
	    SoBRLExportAction exportB;
	    exportB.setGeometryPolicy(SoBRLExportAction::DISPLAY_LEVEL);
	    exportB.apply(viewB.getRenderSceneRoot());
	    if (exportA.getTriangleCount() != 2 ||
		exportB.getTriangleCount() != 3) {
		printf("FAIL: shared LoD view-local exports did not stay independent\n");
		ret = 1;
	    }
	}

	if (!ret) {
	    SoBRLExportAction exactA;
	    exactA.apply(viewA.getRenderSceneRoot());
	    if (exactA.getTriangleCount() != 0 ||
		exactA.getSkippedLodDisplayMeshCount() != 1 ||
		exactA.getSourceBackedFullDetailRequestCount() != 1 ||
		check_source_mesh_request(
		    exactA.getSourceBackedFullDetailRequest(0),
		    "/shared/lod.bot", "lod.bot", 0)) {
		printf("FAIL: shared LoD view-local exact export did not request source mesh\n");
		ret = 1;
	    }
	}

	if (!ret &&
	    (mesh->isLodDisplayActive() ||
	     mesh->lodAvailable.getValue() ||
	     mesh->point.getNum() != 3 ||
	     mesh->coordIndex.getNum() != 3 ||
	     mesh->getTriangleCount() != 1)) {
	    printf("FAIL: shared LoD view-local update mutated shared mesh fields\n");
	    ret = 1;
	}
    }

    shared->unref();
    return ret;
}

static int
test_view_controller_progressive_lod_results(void)
{
    SoSeparator *root = new SoSeparator;
    root->ref();
    SoBRLLodMeshShape *mesh = make_lod_mesh("/progress/lod.bot", "lod.bot");
    root->addChild(mesh);

    BRLObolLodService service;
    if (!service.start(2, TRUE)) {
	printf("FAIL: progressive LoD service did not start\n");
	root->unref();
	return 1;
    }

    int ret = 0;
    {
	BRLObolViewController controller(root, NULL);
	controller.setLodService(&service);
	controller.setLodPolicyRevision(17);
	controller.clearRenderRequest();

	BRLObolLodRequest request = make_request("/progress/lod.bot", "lod.bot");
	request.viewRevision = controller.getLodViewRevision();
	request.policyRevision = controller.getLodPolicyRevision();

	uint64_t generation = service.beginGeneration();
	BRLObolLodTask coarseTask;
	coarseTask.generation = generation;
	coarseTask.request = request;
	coarseTask.realize = mesh_payload_coarse_task;
	BRLObolLodTask refinedTask;
	refinedTask.generation = generation;
	refinedTask.request = request;
	refinedTask.realize = mesh_payload_refined_task;
	refinedTask.debugDelayMilliseconds = 200;

	if (service.submit(coarseTask) == 0 ||
	    service.submit(refinedTask) == 0) {
	    printf("FAIL: progressive LoD service did not accept staged tasks\n");
	    ret = 1;
	}

	if (!ret && wait_for_service_queued(service, 1, TRUE, FALSE))
	    ret = 1;

	if (!ret) {
	    controller.clearRenderRequest();
	    if (controller.processPendingLodResults(1) != 1 ||
		controller.getLastLodAppliedResultCount() != 1 ||
		!controller.isRenderRequested() ||
		strcmp(controller.getRenderReason().getString(),
		       "lod-result") != 0) {
		printf("FAIL: progressive LoD controller did not apply coarse result\n");
		ret = 1;
	    }
	}

	if (!ret) {
	    const BRLObolViewLodState::MeshPayload *payload =
		controller.getViewLodState()->findMesh(mesh);
	    SoBRLExportAction displayExport;
	    displayExport.setGeometryPolicy(SoBRLExportAction::DISPLAY_LEVEL);
	    displayExport.apply(controller.getRenderSceneRoot());
	    if (!payload ||
		payload->activeLevel != 1 ||
		payload->getTriangleCount() != 1 ||
		displayExport.getTriangleCount() != 1 ||
		mesh->lodAvailable.getValue() ||
		mesh->isLodDisplayActive() ||
		mesh->getTriangleCount() != 1) {
		printf("FAIL: progressive LoD coarse result did not become visible view payload\n");
		ret = 1;
	    }
	}

	if (!ret && wait_for_service_queued(service, 1, FALSE, TRUE))
	    ret = 1;

	if (!ret) {
	    controller.clearRenderRequest();
	    if (controller.processPendingLodResults(1) != 1 ||
		controller.getLastLodAppliedResultCount() != 1 ||
		!controller.isRenderRequested() ||
		strcmp(controller.getRenderReason().getString(),
		       "lod-result") != 0) {
		printf("FAIL: progressive LoD controller did not apply refined result\n");
		ret = 1;
	    }
	}

	if (!ret) {
	    const BRLObolViewLodState::MeshPayload *payload =
		controller.getViewLodState()->findMesh(mesh);
	    SoBRLExportAction displayExport;
	    displayExport.setGeometryPolicy(SoBRLExportAction::DISPLAY_LEVEL);
	    displayExport.apply(controller.getRenderSceneRoot());
	    if (!payload ||
		payload->activeLevel != 2 ||
		payload->getTriangleCount() != 3 ||
		displayExport.getTriangleCount() != 3 ||
		service.queuedResultCountForDiagnostics() != 0 ||
		mesh->lodAvailable.getValue() ||
		mesh->isLodDisplayActive() ||
		mesh->getTriangleCount() != 1) {
		printf("FAIL: progressive LoD refined result did not replace coarse view payload\n");
		ret = 1;
	    }
	}
    }

    service.stop();
    root->unref();
    return ret;
}

static int
test_view_local_lod_only_pick(void)
{
    SoSeparator *root = new SoSeparator;
    root->ref();
    SoBRLLodMeshShape *mesh = make_lod_mesh("/pick/lod.bot", "lod.bot");
    mesh->sourceId = 515;
    root->addChild(mesh);

    if (mesh->evictDisplayMeshPreservingSourceMetrics() == 0 ||
	!mesh->isLodBackedMesh() ||
	mesh->getTriangleCount() != 0 ||
	mesh->point.getNum() != 0 ||
	mesh->coordIndex.getNum() != 0) {
	printf("FAIL: view-local LoD pick fixture did not evict shared mesh arrays\n");
	root->unref();
	return 1;
    }

    int ret = 0;
    {
	BRLObolViewController controller(root, NULL);
	BRLObolLodRequest request = make_request("/pick/lod.bot", "lod.bot");
	request.viewRevision = controller.getLodViewRevision();
	request.policyRevision = controller.getLodPolicyRevision();
	BRLObolLodResult result = mesh_payload_variant_result(request, 1, 1);

	SoBRLLodUpdateAction update;
	update.setViewLodState(controller.getViewLodState());
	update.addResult(result);
	update.apply(controller.getRenderSceneRoot());
	if (update.getAppliedResultCount() != 1 ||
	    !controller.getViewLodState()->findMesh(mesh)) {
	    printf("FAIL: view-local LoD pick fixture did not bind display payload\n");
	    ret = 1;
	}

	if (!ret) {
	    SbViewportRegion viewport(100, 100);
	    SoRayPickAction displayPick(viewport);
	    displayPick.setRay(SbVec3f(1.5f, 0.2f, 5.0f),
			       SbVec3f(0.0f, 0.0f, -1.0f));
	    displayPick.apply(controller.getRenderSceneRoot());
	    const SoPickedPoint *pickedPoint = displayPick.getPickedPoint();
	    const SoDetail *rawDetail = pickedPoint ?
					pickedPoint->getDetail(mesh) : NULL;
	    const SoBRLPickDetail *pickDetail =
		(rawDetail &&
		 rawDetail->isOfType(SoBRLPickDetail::getClassTypeId())) ?
		static_cast<const SoBRLPickDetail *>(rawDetail) : NULL;
	    if (!pickDetail ||
		strcmp(pickDetail->getPath().getString(), "/pick/lod.bot") != 0 ||
		pickDetail->getSourceId() != 515 ||
		pickDetail->getPrimitiveKind() != SoBRLPickDetail::FACE ||
		pickDetail->getPrimitiveIndex() != 0) {
		printf("FAIL: display pick did not use view-local LoD-only payload identity\n");
		ret = 1;
	    }
	}

	if (!ret) {
	    mesh->setPickGeometryPolicy(SoBRLMeshShape::PICK_FULL_DETAIL);
	    SbViewportRegion viewport(100, 100);
	    SoRayPickAction exactPick(viewport);
	    exactPick.setRay(SbVec3f(1.5f, 0.2f, 5.0f),
			     SbVec3f(0.0f, 0.0f, -1.0f));
	    exactPick.apply(controller.getRenderSceneRoot());
	    if (exactPick.getPickedPoint()) {
		printf("FAIL: full-detail pick incorrectly used view-local LoD-only payload\n");
		ret = 1;
	    }
	    mesh->setPickGeometryPolicy(SoBRLMeshShape::PICK_DISPLAY_LEVEL);
	}
    }

    root->unref();
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
    if (test_scene_database_source_summary())
	return 1;
    if (test_update_action_service_drain())
	return 1;
    if (test_mesh_lod_request_and_view_info())
	return 1;
    if (test_export_record_metadata())
	return 1;
    if (test_rt_exact_pick_provider())
	return 1;
    if (test_mesh_residency_budget_action())
	return 1;
    if (test_view_controller_mesh_residency_budget_auto())
	return 1;
    if (test_mesh_lod_submit_action())
	return 1;
    if (test_view_controller_source_backed_multi_source_exact_submit())
	return 1;
    if (test_view_controller_source_backed_partial_ready_submit())
	return 1;
    if (test_view_controller_shared_lod_is_view_local())
	return 1;
    if (test_view_controller_progressive_lod_results())
	return 1;
    if (test_view_local_lod_only_pick())
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
