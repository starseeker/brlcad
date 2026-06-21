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
#include <Inventor/actions/SoGetBoundingBoxAction.h>
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
		service.queuedResultCountForDiagnostics() == 1)
	    return 0;
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    printf("FAIL: LoD service did not produce update-action result\n");
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
	const struct rt_view_info *view,
	int *level)
{
    if (!dbip || !name || !view || !level)
	return 1;

    struct rt_mesh_lod *lod = db_mesh_lod_get(dbip, name);
    if (!lod)
	return 1;

    struct rt_mesh_lod_info info = RT_MESH_LOD_INFO_INIT;
    int ret = 0;
    if (rt_mesh_lod_load_view(lod, view, 0) < 0 ||
	    !rt_mesh_lod_info_get(lod, &info) ||
	    info.active_level < 0) {
	ret = 1;
    } else {
	*level = info.active_level;
    }

    rt_mesh_lod_destroy(lod);
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
	    exactExport.getSkippedLodDisplayMeshCount() != 1 ||
	    exactExport.getSourceBackedFullDetailRequestCount() != 0) {
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
	    exactMeasure.getSourceBackedFullDetailRequestCount() != 0 ||
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
	    exactSnap.getSkippedLodDisplayMeshCount() != 1 ||
	    exactSnap.getSourceBackedFullDetailRequestCount() != 0) {
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
	    exactFullSnap.getSkippedLodDisplayMeshCount() != 1 ||
	    exactFullSnap.getSourceBackedFullDetailRequestCount() != 0) {
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

    size_t displayBytes = meshD->estimateDisplayMeshBytes();
    size_t fullDetailBytes = meshD->estimateFullDetailMeshBytes();
    size_t residentBytes = meshD->estimateResidentMeshBytes();
    if (displayBytes == 0 || fullDetailBytes == 0 ||
	    residentBytes != displayBytes + fullDetailBytes) {
	printf("FAIL: resident LoD mesh byte accounting did not include display and full-detail payloads\n");
	root->unref();
	return 1;
    }

    size_t freedBytes = meshD->evictFullDetailMesh();
    BRLObolSourceMeshRequest evictedRequest;
    if (freedBytes != fullDetailBytes ||
	    meshD->hasFullDetailMesh() ||
	    meshD->estimateFullDetailMeshBytes() != 0 ||
	    meshD->estimateResidentMeshBytes() != displayBytes ||
	    !meshD->makeSourceMeshRequest(evictedRequest) ||
	    check_source_mesh_request(evictedRequest, "/mesh/d", "d", 0)) {
	printf("FAIL: full-detail eviction did not preserve source-backed request metadata\n");
	root->unref();
	return 1;
    }

    SoBRLExportAction evictedExactExport;
    evictedExactExport.apply(root);
    if (evictedExactExport.getGeometryPolicy() != SoBRLExportAction::FULL_DETAIL ||
	    evictedExactExport.getTriangleCount() != 3 ||
	    evictedExactExport.getSkippedLodDisplayMeshCount() != 1 ||
	    evictedExactExport.getSourceBackedFullDetailRequestCount() != 1 ||
	    check_source_mesh_request(
		evictedExactExport.getSourceBackedFullDetailRequest(0),
		"/mesh/d", "d", 0)) {
	printf("FAIL: evicted full-detail mesh did not fall back to source-backed exact export request\n");
	root->unref();
	return 1;
    }

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

    BRLObolLodRequest activeDuplicateRequest;
    mesh->makeLodRequest(activeDuplicateRequest,
	    "db://lod-submit-test",
	    2026,
	    55,
	    66,
	    BRLOBOL_LOD_DRAW_SHADED,
	    "rt_mesh_lod",
	    "rt-cache-v1",
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
	db_mesh_lod_clear(dbip);
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
    {
	std::vector<BRLObolLodResult> activeResults;
	service.drainResults(activeResults);
	if (activeResults.size() != 1 ||
		!brlobol_lod_result_matches_request(activeResults[0],
		    activeDuplicateRequest)) {
	    printf("FAIL: LoD submit action active-duplicate fixture did not preserve queued request\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
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

    std::vector<BRLObolLodResult> viewPolicyResults;
    if (service.drainResults(viewPolicyResults) != 1 ||
	    viewPolicyResults.size() != 1) {
	printf("FAIL: LoD submit action result drain failed\n");
	service.stop();
	root->unref();
	db_mesh_lod_clear(dbip);
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
	db_mesh_lod_clear(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    if (viewPolicyResults[0].geometry.activeLevel != expectedViewLevel) {
	printf("FAIL: LoD submit action did not use view-policy active level\n");
	service.stop();
	root->unref();
	db_mesh_lod_clear(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    SoBRLLodUpdateAction update;
    update.setResults(viewPolicyResults);
    if (update.getResultCount() != 1) {
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
	SoBRLMeshLodSubmitAction duplicateSubmit;
	duplicateSubmit.setService(&service);
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
	exactExport.apply(root);
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
		    controllerExport.apply(root);
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
	exactMeasure.apply(root);
	if (exactMeasure.getGeometryPolicy() != SoBRLMeasureAction::FULL_DETAIL ||
		exactMeasure.getTriangleCount() != 1 ||
		exactMeasure.getSkippedLodDisplayMeshCount() != 1 ||
		exactMeasure.getSourceBackedFullDetailRequestCount() != 1 ||
		check_source_mesh_request(
		    exactMeasure.getSourceBackedFullDetailRequest(0),
		    "/lod-submit.bot", "lod-submit.bot", 101)) {
	    printf("FAIL: exact measure did not request source-backed full-detail LoD mesh\n");
	    ret = 1;
	} else {
	    const BRLObolSourceMeshRequest &measureSourceRequest =
		exactMeasure.getSourceBackedFullDetailRequest(0);
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
		    controllerMeasure.apply(root);
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
	limitedMeasureMiss.apply(root);
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
	boundedMeasure.apply(root);
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
	exactSnap.apply(root);
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
		    controllerSnap.apply(root);
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
		sourcePickAction.apply(root);
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
	size_t displayBytes = mesh->estimateDisplayMeshBytes();
	size_t residentBytes = mesh->estimateResidentMeshBytes();
	size_t freedBytes = mesh->evictActiveDisplayMesh();
	BRLObolSourceMeshRequest evictedDisplayRequest;
	SoGetBoundingBoxAction bboxAction(SbViewportRegion(100, 100));
	bboxAction.apply(mesh);
	if (displayBytes == 0 ||
		residentBytes != displayBytes ||
		freedBytes != displayBytes ||
		mesh->isLodDisplayActive() ||
		mesh->getTriangleCount() != 0 ||
		mesh->estimateDisplayMeshBytes() != 0 ||
		mesh->estimateResidentMeshBytes() != 0 ||
		!mesh->makeSourceMeshRequest(evictedDisplayRequest) ||
		check_source_mesh_request(evictedDisplayRequest,
		    "/lod-submit.bot", "lod-submit.bot", 101) ||
		bboxAction.getBoundingBox().isEmpty() ||
		bboxAction.getBoundingBox().getMax()[0] > 1.0f ||
		bboxAction.getBoundingBox().getMax()[1] > 1.0f) {
	    printf("FAIL: active display eviction did not preserve source-backed LoD identity\n");
	    ret = 1;
	}
    }

    if (!ret) {
	SoBRLExportAction evictedDisplayExactExport;
	evictedDisplayExactExport.apply(root);
	if (evictedDisplayExactExport.getGeometryPolicy() !=
		SoBRLExportAction::FULL_DETAIL ||
		evictedDisplayExactExport.getTriangleCount() != 1 ||
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
	evictedDisplayExactMeasure.apply(root);
	if (evictedDisplayExactMeasure.getGeometryPolicy() !=
		SoBRLMeasureAction::FULL_DETAIL ||
		evictedDisplayExactMeasure.getTriangleCount() != 1 ||
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
	evictedDisplayExactSnap.apply(root);
	if (evictedDisplayExactSnap.hasCandidate() ||
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
	evictedDisplayPick.apply(root);
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
	db_mesh_lod_clear(right_dbip);
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
	db_mesh_lod_clear(right_dbip);
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
    db_mesh_lod_clear(right_dbip);
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
	db_mesh_lod_clear(dbip);
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
	db_mesh_lod_clear(dbip);
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
	db_mesh_lod_clear(dbip);
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
    db_mesh_lod_clear(dbip);
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
	if (controller.submitLodRequestsIfNeeded() != 0 ||
		controller.getLastLodSubmittedTaskCount() != 0 ||
		controller.getLastLodSkippedMeshCount() != 1 ||
		strstr(controller.getLastLodDiagnostics().getString(),
		    "current LoD request is already resident") == NULL ||
		controller.isRenderRequested()) {
	    printf("FAIL: LoD view controller did not skip resident changed-scene LoD request\n");
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
	if (service.queuedResultCountForDiagnostics() != 0) {
	    printf("FAIL: LoD view controller queued duplicate resident request result\n");
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
	if (controller.submitLodRequestsIfNeeded() != 0 ||
		controller.getLastLodSubmittedTaskCount() != 0 ||
		controller.getLastLodSkippedMeshCount() != 1 ||
		strstr(controller.getLastLodDiagnostics().getString(),
		    "current LoD request is already resident") == NULL ||
		controller.isRenderRequested()) {
	    printf("FAIL: LoD view controller did not skip resident threshold policy request\n");
	    service.stop();
	    root->unref();
	    db_mesh_lod_clear(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (service.queuedResultCountForDiagnostics() != 0) {
	    printf("FAIL: LoD view controller queued duplicate threshold request result\n");
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
