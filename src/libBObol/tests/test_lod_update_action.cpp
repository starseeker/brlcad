/*          T E S T _ L O D _ U P D A T E _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "BObol.h"
#include "BObol/BLodMeshShape.h"
#include "BObol/BLodRealization.h"
#include "BObol/BLodService.h"
#include "BObol/BLodUpdateAction.h"
#include "BObol/BMaterialObject.h"
#include "BObol/BMeshLodCache.h"
#include "BObol/BMeshLodSubmitAction.h"
#include "BObol/BMeshResidencyAction.h"
#include "BObol/BMeshShape.h"
#include "BObol/BVListShape.h"
#include "BObol/BViewLod.h"
#include "BObol/BViewAttachment.h"
#include "bv.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/time.h"
#include "rt/view.h"
#include "wdb.h"

#include <Obol/cad/SoCADAssembly.h>

#include <Inventor/SbBox.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbViewVolume.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoSeparator.h>

#include <atomic>
#include <chrono>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <string.h>
#include <thread>
#include <vector>

static BObolLodRequest
make_request(const char *path, const char *name)
{
    BObolLodRequest request;

    request.databaseId = "db://lod-update-action-test";
    request.sourceRevision = 3;
    request.sourceContentHash = 0x123456;
    request.objectPath = path;
    request.objectName = name;
    request.viewRevision = 5;
    request.policyRevision = 7;
    request.drawMode = BOBOL_LOD_DRAW_SHADED;
    request.providerId = "lod-update-test";
    request.providerVersion = "1";
    request.qualityTier = BOBOL_LOD_QUALITY_ATTRIBUTES;
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
test_compact_staged_source_lease(void)
{
    std::weak_ptr<const BObolStagedSourceMesh> weakSource;
    {
	BObolCompactOccurrenceStream stream;
	std::shared_ptr<struct rt_bot_internal> bot =
	    std::make_shared<struct rt_bot_internal>();
	std::shared_ptr<BObolStagedSourceMesh> staged =
	    std::make_shared<BObolStagedSourceMesh>();
	staged->owner = bot;
	staged->bot = bot.get();
	staged->assetName = "staged-test.bot";
	staged->sourceRevision = 42;
	staged->byteCount = 4096;
	weakSource = staged;

	BObolSourceMeshRequest request;
	request.stagedSource = staged;
	if (!stream.retainStagedSource(staged) ||
	    stream.stagedSourceByteCount() != staged->byteCount ||
	    request.stagedSource.expired()) {
	    printf("FAIL: compact staged source lease was not retained\n");
	    return 1;
	}

	/* The request itself is intentionally weak: the bounded stream window,
	 * not every copied occurrence record, owns the imported source arrays. */
	bot.reset();
	staged.reset();
	if (weakSource.expired() || request.stagedSource.expired()) {
	    printf("FAIL: compact staged source lease died before stream adoption\n");
	    return 1;
	}
    }
    if (!weakSource.expired()) {
	printf("FAIL: compact staged source lease survived stream destruction\n");
	return 1;
    }
    return 0;
}

static int
find_database_source_summary_by_instance(BObolSceneController &scene,
	const char *instanceKey,
	BObolDatabaseSourceSummary &summary)
{
    summary = BObolDatabaseSourceSummary();
    if (!instanceKey || !instanceKey[0])
	return 0;

    const int count = scene.getDatabaseSourceCount();
    for (int i = 0; i < count; i++) {
	BObolDatabaseSourceSummary candidate;
	if (!scene.getDatabaseSourceSummary(i, candidate) ||
	    !candidate.valid)
	    continue;
	if (bu_strcmp(candidate.instanceKey.getString(), instanceKey) == 0) {
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

static SoNode *
find_named_node(SoNode *node, const char *name)
{
    if (!node || !name)
	return NULL;
    if (bu_strcmp(node->getName().getString(), name) == 0)
	return node;
    if (!node->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *found = find_named_node(group->getChild(i), name);
	if (found)
	    return found;
    }
    return NULL;
}

static SoGroup *
find_node_parent(SoNode *node, SoNode *wanted)
{
    if (!node || !wanted || !node->isOfType(SoGroup::getClassTypeId()))
	return NULL;

    SoGroup *group = static_cast<SoGroup *>(node);
    if (group->findChild(wanted) >= 0)
	return group;
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoGroup *parent = find_node_parent(group->getChild(i), wanted);
	if (parent)
	    return parent;
    }
    return NULL;
}

static int
check_controller_light_order(BObolViewController &controller,
			     SoCamera *camera)
{
    SoNode *renderRoot = controller.getRenderRoot();
    SoNode *environment =
	find_named_node(renderRoot, "BObolRenderEnvironment");
    SoNode *headlight = find_named_node(renderRoot, "BObolHeadlight");
    SoNode *fill = find_named_node(renderRoot, "BObolStudioFill");
    SoNode *rim = find_named_node(renderRoot, "BObolStudioRim");
    SoNode *sceneLights = find_named_node(renderRoot, "BObolSceneLights");
    SoGroup *parent = find_node_parent(renderRoot, headlight);

    if (!renderRoot || !camera || !environment || !headlight || !fill || !rim ||
	!sceneLights || !parent ||
	!headlight->isOfType(SoDirectionalLight::getClassTypeId()) ||
	!fill->isOfType(SoDirectionalLight::getClassTypeId()) ||
	!rim->isOfType(SoDirectionalLight::getClassTypeId()) ||
	find_node_parent(renderRoot, camera) != parent ||
	find_node_parent(renderRoot, environment) != parent ||
	find_node_parent(renderRoot, fill) != parent ||
	find_node_parent(renderRoot, rim) != parent ||
	find_node_parent(renderRoot, sceneLights) != parent) {
	return 1;
    }

    const int environmentIndex = parent->findChild(environment);
    const int cameraIndex = parent->findChild(camera);
    const int headlightIndex = parent->findChild(headlight);
    const int fillIndex = parent->findChild(fill);
    const int rimIndex = parent->findChild(rim);
    const int sceneLightsIndex = parent->findChild(sceneLights);
    return environmentIndex != 0 ||
	cameraIndex <= environmentIndex ||
	headlightIndex != cameraIndex + 1 ||
	fillIndex != headlightIndex + 1 ||
	rimIndex != fillIndex + 1 ||
	sceneLightsIndex <= rimIndex;
}

static int
test_view_controller_fixed_gl_light_order(void)
{
    SoSeparator *root = new SoSeparator;
    root->ref();
    int ret = 0;
    {
	SoOrthographicCamera *camera = new SoOrthographicCamera;
	BObolViewController controller(root, camera);
	if (check_controller_light_order(controller, camera)) {
	    printf("FAIL: camera light rig does not traverse immediately after camera\n");
	    ret = 1;
	}
	std::vector<BObolSceneLightRealization> cameraLights;
	controller.getCameraLights(cameraLights);
	if (!ret && (controller.getLightingProfile() !=
		BObolViewController::LIGHTING_STUDIO ||
	    cameraLights.size() != 3u)) {
	    printf("FAIL: controller did not initialize the studio light rig\n");
	    ret = 1;
	}
	controller.setLightingProfile(BObolViewController::LIGHTING_MGED);
	controller.getCameraLights(cameraLights);
	if (!ret && (cameraLights.size() != 1u ||
	    std::fabs(controller.getLightingAmbientIntensity() - 0.3f) >
		1.0e-6f)) {
	    printf("FAIL: controller MGED light profile is not the historical rig\n");
	    ret = 1;
	}
	controller.setLightingProfile(BObolViewController::LIGHTING_STUDIO);
	controller.setLightingEnabled(FALSE);
	controller.getCameraLights(cameraLights);
	if (!ret && !cameraLights.empty()) {
	    printf("FAIL: disabled controller lighting still exported camera lights\n");
	    ret = 1;
	}
	controller.setLightingEnabled(TRUE);
	controller.getCameraLights(cameraLights);
	if (!ret && cameraLights.size() != 3u) {
	    printf("FAIL: re-enabled controller lighting did not restore studio rig\n");
	    ret = 1;
	}

	SoOrthographicCamera *replacement = new SoOrthographicCamera;
	controller.setCamera(replacement);
	if (!ret && check_controller_light_order(controller, replacement)) {
	    printf("FAIL: replacement camera broke fixed-GL light ordering\n");
	    ret = 1;
	}
    }
    root->unref();
    return ret;
}

static int
test_view_controller_presentation_deadline_contract(void)
{
    SoSeparator *root = new SoSeparator;
    root->ref();
    int ret = 0;
    {
	BObolViewController controller(root, NULL);
	if (controller.getInteractivePresentationFrameDeadline() !=
		40000000ULL ||
	    controller.getStablePresentationFrameDeadline() !=
		100000000ULL ||
	    controller.getCurrentPresentationFrameDeadline() !=
		100000000ULL) {
	    printf("FAIL: presentation deadline defaults\n");
	    ret = 1;
	}

	controller.setPresentationFrameDeadlines(5000000ULL, 9000000ULL);
	controller.setLodAutoSubmit(TRUE);
	controller.beginLodInteraction();
	if (!ret &&
		controller.getCurrentPresentationFrameDeadline() != 5000000ULL) {
	    printf("FAIL: interactive presentation deadline selection\n");
	    ret = 1;
	}
	controller.endLodInteraction();

	controller.clearRenderRequest();
	const size_t budgetBeforeInterruptedPresentation =
	    controller.getCurrentLodRenderCostBudget();
	controller.notePresentationRenderInterrupted(123000000ULL);
	if (!ret &&
		(controller.getInterruptedPresentationFrameCount() != 1 ||
		 controller.getLastInterruptedPresentationTimeNanoseconds() !=
		     123000000ULL ||
		 !controller.isRenderRequested() ||
		 bu_strcmp(controller.getRenderReason().getString(),
		     "render-deadline") != 0 ||
		 controller.getCurrentPresentationFrameDeadline() !=
		     7500000ULL ||
		 controller.getCurrentLodRenderCostBudget() !=
		     budgetBeforeInterruptedPresentation)) {
	    printf("FAIL: interrupted presentation did not retain the "
		   "progressive retry barrier and proven scene budget\n");
	    ret = 1;
	}
	controller.notePresentationRenderInterrupted(0);
	if (!ret && controller.getInterruptedPresentationFrameCount() != 1) {
	    printf("FAIL: zero-duration presentation counted as interrupted\n");
	    ret = 1;
	}
    }
    root->unref();
    return ret;
}

static int
check_source_mesh_request(const BObolSourceMeshRequest &request,
			  const char *path,
			  const char *name,
			  uint32_t sourceId)
{
    if (bu_strcmp(request.path.getString(), path) != 0 ||
	bu_strcmp(request.sourceName.getString(), name) != 0 ||
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
check_source_full_detail_lod_request(const BObolLodRequest &request,
				     const char *path,
				     const char *name)
{
    if (bu_strcmp(request.objectPath.getString(), path) != 0 ||
	bu_strcmp(request.objectName.getString(), name) != 0 ||
	bu_strcmp(request.providerId.getString(), "rt_source_full_detail") != 0 ||
	bu_strcmp(request.providerVersion.getString(), "direct-bot-v1") != 0 ||
	request.qualityTier != BOBOL_LOD_QUALITY_FULL_DETAIL ||
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
has_lod_provider_param(const BObolLodRequest &request, const char *name)
{
    for (size_t i = 0; i < request.providerParams.size(); i++) {
	if (bu_strcmp(request.providerParams[i].name.getString(), name) == 0)
	    return 1;
    }

    return 0;
}

static BObolLodResult
attributes_result(const BObolLodRequest &request, const char *value)
{
    std::vector<BObolLodAttribute> attributes;
    BObolLodAttribute attribute;

    attribute.name = "display.intent";
    attribute.value = value;
    attributes.push_back(attribute);
    return bobol_lod_attributes_result(request, attributes);
}

static BObolLodResult
mesh_payload_result(const BObolLodRequest &request)
{
    BObolLodResult result;

    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = BOBOL_LOD_RESULT_MESH;
    result.qualityTier = BOBOL_LOD_QUALITY_FAST_DISPLAY;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.geometry.kind = BOBOL_LOD_GEOMETRY_OBOL_MESH;
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

static BObolLodResult
mesh_payload_task(const BObolLodRequest &request, void *UNUSED(userData))
{
    return mesh_payload_result(request);
}

static BObolLodResult
aabb_payload_task(const BObolLodRequest &request, void *UNUSED(userData))
{
    return bobol_lod_aabb_result(request, request.bounds, NULL);
}

static BObolLodResult
mesh_payload_variant_result(const BObolLodRequest &request,
			    int activeLevel,
			    int triangleCount)
{
    BObolLodResult result = mesh_payload_result(request);

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

static BObolLodResult
mesh_payload_coarse_task(const BObolLodRequest &request, void *UNUSED(userData))
{
    return mesh_payload_variant_result(request, 1, 1);
}

static BObolLodResult
mesh_payload_refined_task(const BObolLodRequest &request, void *UNUSED(userData))
{
    return mesh_payload_variant_result(request, 2, 3);
}

static BObolLodResult
mesh_payload_requested_refined_task(const BObolLodRequest &request,
	void *UNUSED(userData))
{
    return mesh_payload_variant_result(request, request.requestedLevel, 3);
}

static BObolLodResult
copied_progressive_result_task(const BObolLodRequest &request, void *userData)
{
    const BObolLodResult *prototype =
	static_cast<const BObolLodResult *>(userData);
    if (!prototype)
	return BObolLodResult();
    BObolLodResult result = *prototype;
    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    return result;
}

static void
count_frame_request(void *userData, const char *UNUSED(reason))
{
    std::atomic<unsigned int> *count =
	static_cast<std::atomic<unsigned int> *>(userData);
    if (count)
	count->fetch_add(1, std::memory_order_relaxed);
}

struct large_payload_task_data {
    size_t pointCount;
    std::atomic<uintptr_t> pointStorage;
    std::atomic<uintptr_t> indexStorage;

    large_payload_task_data(void) :
	pointCount(0),
	pointStorage(0),
	indexStorage(0)
    {
    }
};

static BObolLodResult
large_mesh_payload_task(const BObolLodRequest &request, void *userData)
{
    large_payload_task_data *data =
	static_cast<large_payload_task_data *>(userData);
    BObolLodResult result = mesh_payload_result(request);
    result.mesh.clear();
    const size_t pointCount = data ? data->pointCount : 0;
    result.mesh.points.resize(pointCount, SbVec3f(0.0f, 0.0f, 0.0f));
    const size_t indexCount = pointCount >= 3 ?
	pointCount - (pointCount % 3) : 0;
    result.mesh.coordIndex.resize(indexCount, 0);
    result.counts.pointCount = result.mesh.points.size();
    result.counts.faceCount = result.mesh.coordIndex.size() / 3;
    if (data) {
	data->pointStorage.store(reinterpret_cast<uintptr_t>(
	    result.mesh.points.data()));
	data->indexStorage.store(reinterpret_cast<uintptr_t>(
	    result.mesh.coordIndex.data()));
    }
    return result;
}

static BObolLodResult
source_full_detail_result(const BObolLodRequest &request)
{
    BObolLodResult result;

    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = BOBOL_LOD_RESULT_FULL_DETAIL;
    result.qualityTier = BOBOL_LOD_QUALITY_FULL_DETAIL;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.geometry.kind = BOBOL_LOD_GEOMETRY_OBOL_MESH;
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

static BObolLodResult
source_full_detail_task(const BObolLodRequest &request, void *UNUSED(userData))
{
    return source_full_detail_result(request);
}

static BObolLodResult
aabb_task(const BObolLodRequest &request, void *UNUSED(userData))
{
    BObolLodCounts counts;

    counts.faceCount = 1;
    counts.pointCount = 3;
    return bobol_lod_aabb_result(request, request.bounds, &counts);
}

static int
wait_for_service(BObolLodService &service)
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
wait_for_service_queued(BObolLodService &service,
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
wait_for_service_delayed(BObolLodService &service)
{
    for (int i = 0; i < 400; i++) {
	if (service.delayedTaskCountForDiagnostics() > 0)
	    return 0;
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    printf("FAIL: LoD service did not start delayed source work\n");
    return 1;
}

static int
wait_for_service_idle(BObolLodService &service)
{
    for (int i = 0; i < 400; i++) {
	if (service.inFlightCount() == 0 &&
	    service.pendingTaskCountForDiagnostics() == 0 &&
	    service.queuedResultCountForDiagnostics() == 0 &&
	    service.queuedCacheWriteCountForDiagnostics() == 0 &&
	    service.delayedTaskCountForDiagnostics() == 0 &&
	    service.activeRequestCountForDiagnostics() == 0)
	    return 0;
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    printf("FAIL: LoD service did not drain cancelled source work\n");
    return 1;
}

static int
submit_source_full_detail_task(BObolLodService &service,
			       const BObolLodRequest &request)
{
    BObolLodTask task;

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

    struct BObolMeshLod *lod = bobol_mesh_lod_get(dbip, name);
    if (!lod)
	return 1;

    struct BObolMeshLodInfo info = BOBOL_MESH_LOD_INFO_INIT;
    int ret = 0;
    if (bobol_mesh_lod_load_view(lod, view, 0) < 0 ||
	!bobol_mesh_lod_info_get(lod, &info) ||
	info.active_level < 0) {
	ret = 1;
    } else {
	*level = info.active_level;
    }

    bobol_mesh_lod_destroy(lod);
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

    BObolRtPickCache pickCache;
    if (!pickCache.prepare(dbip, paths) ||
	!pickCache.isReady() ||
	pickCache.getObjectPathCount() != 1) {
	printf("FAIL: RT exact pick provider did not prepare reusable pick cache\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    BObolRtPickResult pick;
    if (!pickCache.pickRay(pick,
			   SbVec3f(0.0f, 0.0f, 5.0f),
			   SbVec3f(0.0f, 0.0f, -1.0f)) ||
	!pick.hit ||
	fabsf(pick.distance - 3.0f) > 1.0e-5f ||
	fabsf(pick.point[2] - 2.0f) > 1.0e-5f ||
	pick.normal[2] < 0.9f ||
	!pick.uvValid || pick.uv[0] < 0.0f || pick.uv[0] > 1.0f ||
	pick.uv[1] < 0.0f || pick.uv[1] > 1.0f ||
	pick.detail.getPrimitiveKind() != SoBRLPickDetail::IMPLICIT_SOLID ||
	pick.detail.getRegionId() != 77 ||
	pick.detail.getAirCode() != 2 ||
	pick.detail.getMaterialId() != 33 ||
	pick.detail.getLos() != 100 ||
	bu_strcmp(pick.detail.getSourceName().getString(), "implicit.s") != 0 ||
	bu_strcmp(pick.detail.getSourceType().getString(), "sph") != 0 ||
	!strstr(pick.detail.getPath().getString(), "implicit.r")) {
	printf("FAIL: RT exact pick provider did not return implicit comb hit identity\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    BObolRtPickResult exitPick;
    if (!pickCache.pickRay(exitPick,
	    SbVec3f(0.0f, 0.0f, 1.9999f),
	    SbVec3f(0.0f, 0.0f, -1.0f), NULL, 0, NULL, 1.0e-4f) ||
	!exitPick.hit || fabsf(exitPick.distance - 3.9999f) > 1.0e-4f ||
	fabsf(exitPick.point[2] + 2.0f) > 1.0e-4f ||
	exitPick.normal[2] > -0.9f ||
	!exitPick.uvValid || exitPick.uv[0] < 0.0f ||
	exitPick.uv[0] > 1.0f || exitPick.uv[1] < 0.0f ||
	exitPick.uv[1] > 1.0f ||
	exitPick.detail.getPrimitiveKind() != SoBRLPickDetail::IMPLICIT_SOLID) {
	printf("FAIL: RT exact pick minimum distance did not select solid exit boundary\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    const SbPlane sectionPlane(SbVec3f(0.0f, 0.0f, -1.0f),
	SbVec3f(0.0f, 0.0f, 0.0f));
    BObolRtPickResult sectionPick;
    if (!pickCache.pickRay(sectionPick,
	    SbVec3f(0.0f, 0.0f, 5.0f),
	    SbVec3f(0.0f, 0.0f, -1.0f), &sectionPlane, 1) ||
	!sectionPick.hit || fabsf(sectionPick.distance - 5.0f) > 1.0e-5f ||
	fabsf(sectionPick.point[2]) > 1.0e-5f ||
	sectionPick.normal[2] > -0.9f ||
	sectionPick.uvValid ||
	sectionPick.detail.getPrimitiveIndex() != -1) {
	printf("FAIL: RT exact pick did not realize a cut-plane section hit\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    const SbPlane removePlane(SbVec3f(0.0f, 0.0f, 1.0f),
	SbVec3f(0.0f, 0.0f, 3.0f));
    BObolRtPickResult clippedMiss;
    if (pickCache.pickRay(clippedMiss,
	    SbVec3f(0.0f, 0.0f, 5.0f),
	    SbVec3f(0.0f, 0.0f, -1.0f), &removePlane, 1) ||
	clippedMiss.hit) {
	printf("FAIL: RT exact pick selected geometry removed by cut plane\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    BObolRtPickResult miss;
    if (pickCache.pickRay(miss,
			  SbVec3f(10.0f, 10.0f, 5.0f),
			  SbVec3f(0.0f, 0.0f, -1.0f)) || miss.hit) {
	printf("FAIL: RT exact pick provider reported a miss ray as a hit\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    BObolRtPickResult wrapperPick;
    if (!bobol_pick_rt_ray(wrapperPick, dbip, paths,
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

    BObolLodResult resultA =
	attributes_result(make_request("/mesh/a", "a"), "proxy");
    BObolLodResult resultB =
	attributes_result(make_request("/mesh/b", "b"), "stale");
    resultB.providerStatus = BOBOL_LOD_PROVIDER_STALE;
    resultB.stale = TRUE;
    resultB.diagnostic = "unit-test stale result";
    BObolLodResult resultC =
	attributes_result(make_request("/mesh/c", "c"), "cache-miss");
    resultC.providerStatus = BOBOL_LOD_PROVIDER_CACHE_MISS;
    resultC.diagnostic = "unit-test cache miss";
    BObolLodResult resultD =
	mesh_payload_result(make_request("/mesh/d", "d"));
    BObolLodResult resultMissing =
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

    BObolViewLodState viewState;
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

    const BObolViewLodState::MeshPayload *meshDPayload =
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

    BObolLodRequest meshDRequest;
    meshD->makeLodRequest(meshDRequest,
			  "db://lod-update-action-test",
			  10,
			  11,
			  12,
			  BOBOL_LOD_DRAW_SHADED,
			  "bobol_mesh_lod",
			  "bobol-cache-v1",
			  BOBOL_LOD_QUALITY_FAST_DISPLAY);
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
	bu_strcmp(displaySnap.getPath().getString(), "/mesh/d") != 0 ||
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
	bu_strcmp(exactFullSnap.getPath().getString(), "/mesh/d") != 0 ||
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
	bu_strcmp(pickDetail->getPath().getString(), "/mesh/d") != 0 ||
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
	bu_strcmp(pickDetail->getPath().getString(), "/mesh/d") != 0 ||
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

    BObolSceneController scene(root);
    const uint64_t initialStructuralRevision = scene.getStructuralRevision();
    const uint64_t initialFrameRevision = scene.getFrameRevision();
    BObolSceneSummary sceneSummary;
    BObolDatabaseSourceSummary summary;
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
	scene.findDatabaseSourceRoutingId(
	    source->getCompactSourceRoutingId()) != source ||
	scene.findDatabaseSourceRoutingId(0) != NULL ||
	!scene.getDatabaseSourceSummary(0, summary) ||
	!summary.valid ||
	bu_strcmp(summary.path.getString(), "/summary/source") != 0 ||
	!summary.hasParent ||
	summary.drawTreeDepth != 1 ||
	bu_strcmp(summary.parentGroupPath.getString(), "/") != 0 ||
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
	BObolViewController view(root, NULL);
	BObolDatabaseSourceSummary viewSummary;
	if (view.getDatabaseSourceCount() != 1 ||
	    view.getDatabaseSource(0) != source ||
	    !view.getDatabaseSourceSummary(0, viewSummary) ||
	    !viewSummary.hasParent ||
	    bu_strcmp(viewSummary.parentGroupPath.getString(), "/") != 0) {
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
	bu_strcmp(realizedShape->ownerSourcePath.getValue().getString(),
	       "/summary/source") != 0 ||
	realizedShape->ownerSourceRevision.getValue() != 11 ||
	realizedShape->ownerInputsRevision.getValue() != 2 ||
	realizedShape->ownerViewRevision.getValue() != 3 ||
	realizedShape->ownerRealizedSourceRevision.getValue() != 11 ||
	realizedShape->ownerRealizedInputsRevision.getValue() != 2 ||
	realizedShape->ownerRealizedViewRevision.getValue() != 3 ||
	realizedShape->ownerRealizationStatus.getValue() !=
	SoBRLDatabaseSource::REALIZED ||
	bu_strcmp(realizedShape->ownerRealizationIdentity.getValue().getString(),
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
	bu_strcmp(summary.realizationIdentity.getString(),
	       firstIdentity.getString()) != 0 ||
	!realizedShape->ownerSourceStale.getValue() ||
	!(realizedShape->ownerStaleReason.getValue() &
	  SoBRLDatabaseSource::STALE_INPUTS) ||
	realizedShape->ownerRealizationStatus.getValue() !=
	SoBRLDatabaseSource::UNREALIZED ||
	bu_strcmp(realizedShape->ownerRealizationIdentity.getValue().getString(),
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
	bu_strcmp(summary.realizationIdentity.getString(),
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
	bu_strcmp(realizedShape->ownerRealizationIdentity.getValue().getString(),
	       summary.realizationIdentity.getString()) != 0 ||
	realizedShape->ownerSourceStale.getValue() ||
	realizedShape->ownerStaleReason.getValue() !=
	SoBRLDatabaseSource::STALE_NONE) {
	printf("FAIL: realized shape node should refresh mirrored source owner state\n");
	root->unref();
	return 1;
    }

    BObolRealizedShapeSummary shapeSummary;
    if (source->getRealizedShapeSummaryCount() != 1 ||
	scene.getRealizedShapeSummaryCount() != 1 ||
	!source->getRealizedShapeSummary(0, shapeSummary) ||
	!shapeSummary.valid ||
	shapeSummary.shapeKind != BObolRealizedShapeSummary::SHAPE_VLIST ||
	bu_strcmp(shapeSummary.path.getString(), "/summary/source") != 0 ||
	shapeSummary.ownerSourceIndex != -1 ||
	bu_strcmp(shapeSummary.ownerSourcePath.getString(),
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
	bu_strcmp(shapeSummary.ownerRealizationIdentity.getString(),
	       summary.realizationIdentity.getString()) != 0 ||
	shapeSummary.ownerSourceStale ||
	shapeSummary.ownerStaleReason != SoBRLDatabaseSource::STALE_NONE ||
	bu_strcmp(shapeSummary.sourceName.getString(), "summary/source") != 0 ||
	bu_strcmp(shapeSummary.sourceType.getString(), "prototype") != 0 ||
	bu_strcmp(shapeSummary.displayName.getString(), "summary/source") != 0 ||
	!shapeSummary.localSource ||
	!shapeSummary.nonDatabaseSource ||
	shapeSummary.databaseIntent ||
	shapeSummary.drawMode != BOBOL_LOD_DRAW_DIAGNOSTIC ||
	bu_strcmp(shapeSummary.recordRole.getString(), "prototype") != 0 ||
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
    BObolRealizedShapeSummary sceneShapeSummary;
    if (!scene.getRealizedShapeSummary(0, sceneShapeSummary) ||
	!sceneShapeSummary.valid ||
	sceneShapeSummary.shapeKind !=
	BObolRealizedShapeSummary::SHAPE_VLIST ||
	bu_strcmp(sceneShapeSummary.path.getString(),
	       shapeSummary.path.getString()) != 0 ||
	sceneShapeSummary.ownerSourceIndex != 0 ||
	bu_strcmp(sceneShapeSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	bu_strcmp(sceneShapeSummary.ownerRealizationIdentity.getString(),
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
    BObolAuxiliaryLineSetDisplayState auxDisplay;
    auxDisplay.valid = TRUE;
    auxDisplay.drawMode = BOBOL_LOD_DRAW_SHADED;
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
	bu_strcmp(auxShape->recordRole.getValue().getString(),
	       "auxiliary") != 0 ||
	source->getRealizedShapeSummaryCount() != 2 ||
	scene.getRealizedShapeSummaryCount() != 2 ||
	!source->getRealizedShapeSummary(1, shapeSummary) ||
	!shapeSummary.valid ||
	shapeSummary.shapeKind !=
	BObolRealizedShapeSummary::SHAPE_VLIST ||
	bu_strcmp(shapeSummary.sourceName.getString(), "summary_aux") != 0 ||
	bu_strcmp(shapeSummary.sourceType.getString(),
	       "auxiliary-line-set") != 0 ||
	bu_strcmp(shapeSummary.geometryName.getString(), "summary_aux") != 0 ||
	bu_strcmp(shapeSummary.recordRole.getString(), "auxiliary") != 0 ||
	!shapeSummary.databaseIntent ||
	shapeSummary.nonDatabaseSource ||
	shapeSummary.drawMode != BOBOL_LOD_DRAW_SHADED ||
	shapeSummary.visible ||
	!shapeSummary.highlighted ||
	shapeSummary.materialRevision != 88 ||
	!shapeSummary.materialColorValid ||
	fabsf(shapeSummary.materialColor[0] - 0.25f) > 1.0e-6f ||
	fabsf(shapeSummary.materialColor[1] - 0.5f) > 1.0e-6f ||
	fabsf(shapeSummary.materialColor[2] - 0.75f) > 1.0e-6f ||
	shapeSummary.pointCount != 2 ||
	shapeSummary.segmentCount != 1 ||
	bu_strcmp(shapeSummary.ownerSourcePath.getString(),
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
    summaryMesh->drawMode = BOBOL_LOD_DRAW_SHADED;
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
	shapeSummary.shapeKind != BObolRealizedShapeSummary::SHAPE_MESH ||
	bu_strcmp(shapeSummary.path.getString(), "/summary/mesh") != 0 ||
	shapeSummary.ownerSourceIndex != -1 ||
	bu_strcmp(shapeSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	shapeSummary.ownerInputsRevision != 4 ||
	shapeSummary.ownerRealizationStatus !=
	SoBRLDatabaseSource::REALIZED ||
	bu_strcmp(shapeSummary.sourceName.getString(), "mesh") != 0 ||
	bu_strcmp(shapeSummary.sourceType.getString(), "bot") != 0 ||
	bu_strcmp(shapeSummary.displayName.getString(), "Summary Mesh") != 0 ||
	bu_strcmp(shapeSummary.cacheIdentity.getString(),
	       "cache://summary/mesh#44") != 0 ||
	bu_strcmp(shapeSummary.sourceIdentity.getString(),
	       "db://summary/mesh") != 0 ||
	!shapeSummary.databaseIntent ||
	shapeSummary.localSource ||
	shapeSummary.nonDatabaseSource ||
	shapeSummary.drawMode != BOBOL_LOD_DRAW_SHADED ||
	bu_strcmp(shapeSummary.recordRole.getString(), "database") != 0 ||
	bu_strcmp(shapeSummary.geometryKind.getString(), "surface") != 0 ||
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
	BObolRealizedShapeSummary::SHAPE_MESH ||
	bu_strcmp(sceneShapeSummary.path.getString(), "/summary/mesh") != 0 ||
	sceneShapeSummary.ownerSourceIndex != 0 ||
	bu_strcmp(sceneShapeSummary.ownerSourcePath.getString(),
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

    BObolRealizedMaterialSummary materialSummary;
    SbString propertyGroup;
    SbString propertyName;
    SbString propertyValue;
    if (source->getRealizedMaterialSummaryCount() != 1 ||
	scene.getRealizedMaterialSummaryCount() != 1 ||
	!source->getRealizedMaterialSummary(0, materialSummary) ||
	!materialSummary.valid ||
	bu_strcmp(materialSummary.sourcePath.getString(),
	       "/summary/material") != 0 ||
	bu_strcmp(materialSummary.sourceName.getString(), "mat") != 0 ||
	bu_strcmp(materialSummary.sourceType.getString(), "material") != 0 ||
	materialSummary.sourceId != 77 ||
	materialSummary.ownerSourceIndex != -1 ||
	bu_strcmp(materialSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	materialSummary.ownerInputsRevision != 4 ||
	materialSummary.ownerRealizationStatus !=
	SoBRLDatabaseSource::REALIZED ||
	bu_strcmp(materialSummary.ownerRealizationIdentity.getString(),
	       summary.realizationIdentity.getString()) != 0 ||
	bu_strcmp(materialSummary.materialName.getString(), "matte") != 0 ||
	bu_strcmp(materialSummary.parentName.getString(), "base") != 0 ||
	bu_strcmp(materialSummary.materialSource.getString(), "database") != 0 ||
	materialSummary.propertyCount != 2 ||
	!source->getRealizedMaterialProperty(0, 0, propertyGroup,
		propertyName, propertyValue) ||
	bu_strcmp(propertyGroup.getString(), "physical") != 0 ||
	bu_strcmp(propertyName.getString(), "density") != 0 ||
	bu_strcmp(propertyValue.getString(), "7.8") != 0 ||
	source->getRealizedMaterialProperty(0, 2, propertyGroup,
					    propertyName, propertyValue) ||
	source->getRealizedMaterialSummary(1, materialSummary)) {
	printf("FAIL: database source should summarize realized material objects\n");
	root->unref();
	return 1;
    }

    BObolRealizedMaterialSummary sceneMaterialSummary;
    if (!scene.getRealizedMaterialSummary(0, sceneMaterialSummary) ||
	!sceneMaterialSummary.valid ||
	sceneMaterialSummary.ownerSourceIndex != 0 ||
	bu_strcmp(sceneMaterialSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	bu_strcmp(sceneMaterialSummary.sourcePath.getString(),
	       "/summary/material") != 0 ||
	!scene.getRealizedMaterialProperty(0, 1, propertyGroup,
					   propertyName, propertyValue) ||
	bu_strcmp(propertyGroup.getString(), "optical") != 0 ||
	bu_strcmp(propertyName.getString(), "finish") != 0 ||
	bu_strcmp(propertyValue.getString(), "matte") != 0 ||
	scene.getRealizedMaterialProperty(1, 0, propertyGroup,
					  propertyName, propertyValue) ||
	scene.getRealizedMaterialSummary(1, sceneMaterialSummary)) {
	printf("FAIL: scene controller should forward realized material object summaries\n");
	root->unref();
	return 1;
    }

    BObolSceneTreeSummary treeSummary;
    if (source->getRealizedTreeSummaryCount() != 4 ||
	!source->getRealizedTreeSummary(0, treeSummary) ||
	!treeSummary.valid ||
	treeSummary.nodeKind !=
	BObolSceneTreeSummary::NODE_DATABASE_SOURCE ||
	!treeSummary.isGroup ||
	treeSummary.isShape ||
	treeSummary.hasParent ||
	treeSummary.drawTreeDepth != 0 ||
	treeSummary.childCount != 3 ||
	treeSummary.ownerSourceIndex != -1 ||
	bu_strcmp(treeSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	bu_strcmp(treeSummary.path.getString(), "/summary/source") != 0 ||
	!source->getRealizedTreeSummary(1, treeSummary) ||
	!treeSummary.valid ||
	treeSummary.nodeKind != BObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	!treeSummary.isShape ||
	!treeSummary.hasParent ||
	treeSummary.drawTreeDepth != 1 ||
	treeSummary.childCount != 0 ||
	bu_strcmp(treeSummary.path.getString(), "/summary/source") != 0 ||
	!source->getRealizedTreeSummary(2, treeSummary) ||
	treeSummary.nodeKind != BObolSceneTreeSummary::NODE_MESH_SHAPE ||
	bu_strcmp(treeSummary.path.getString(), "/summary/mesh") != 0 ||
	!source->getRealizedTreeSummary(3, treeSummary) ||
	treeSummary.nodeKind !=
	BObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
	!treeSummary.isMaterialObject ||
	bu_strcmp(treeSummary.path.getString(), "/summary/material") != 0 ||
	bu_strcmp(treeSummary.displayName.getString(), "matte") != 0 ||
	source->getRealizedTreeSummary(4, treeSummary)) {
	printf("FAIL: database source should summarize realized source tree nodes\n");
	root->unref();
	return 1;
    }

    if (scene.getSceneTreeSummaryCount() != 6 ||
	!scene.getSceneTreeSummary(0, treeSummary) ||
	!treeSummary.valid ||
	treeSummary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	!treeSummary.isGroup ||
	treeSummary.hasParent ||
	treeSummary.drawTreeDepth != 0 ||
	treeSummary.childCount != 2 ||
	!scene.getSceneTreeSummary(1, treeSummary) ||
	treeSummary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	!treeSummary.hasParent ||
	treeSummary.ownerSourceIndex != -1 ||
	treeSummary.drawTreeDepth != 1 ||
	!scene.getSceneTreeSummary(2, treeSummary) ||
	treeSummary.nodeKind !=
	BObolSceneTreeSummary::NODE_DATABASE_SOURCE ||
	treeSummary.ownerSourceIndex != 0 ||
	bu_strcmp(treeSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	treeSummary.drawTreeDepth != 1 ||
	!treeSummary.hasParent ||
	!scene.getSceneTreeSummary(3, treeSummary) ||
	treeSummary.nodeKind != BObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	treeSummary.ownerSourceIndex != 0 ||
	treeSummary.drawTreeDepth != 2 ||
	!scene.getSceneTreeSummary(4, treeSummary) ||
	treeSummary.nodeKind != BObolSceneTreeSummary::NODE_MESH_SHAPE ||
	bu_strcmp(treeSummary.path.getString(), "/summary/mesh") != 0 ||
	!scene.getSceneTreeSummary(5, treeSummary) ||
	treeSummary.nodeKind !=
	BObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
	bu_strcmp(treeSummary.path.getString(), "/summary/material") != 0 ||
	scene.getSceneTreeSummary(6, treeSummary)) {
	printf("FAIL: scene controller should summarize Obol scene tree nodes\n");
	root->unref();
	return 1;
    }

    BObolSceneDisplaySummary displaySummary;
    if (source->getRealizedDisplaySummaryCount() !=
	source->getRealizedTreeSummaryCount() ||
	!source->getRealizedDisplaySummary(0, displaySummary) ||
	!displaySummary.valid ||
	!displaySummary.isDatabaseSource ||
	!displaySummary.hasDrawIntent ||
	bu_strcmp(displaySummary.intentPath.getString(),
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
	bu_strcmp(displaySummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	!source->getRealizedDisplaySummary(1, displaySummary) ||
	displaySummary.nodeKind != BObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	!displaySummary.hasDrawIntent ||
	bu_strcmp(displaySummary.intentPath.getString(),
	       "/summary/source") != 0 ||
	displaySummary.intentDrawMode != BOBOL_LOD_DRAW_DIAGNOSTIC ||
	!displaySummary.visible ||
	displaySummary.highlighted ||
	displaySummary.materialColor[0] < 0.99f ||
	!source->getRealizedDisplaySummary(2, displaySummary) ||
	displaySummary.nodeKind != BObolSceneTreeSummary::NODE_MESH_SHAPE ||
	bu_strcmp(displaySummary.intentPath.getString(), "/summary/mesh") != 0 ||
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
	BObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
	displaySummary.hasDrawIntent ||
	bu_strcmp(displaySummary.path.getString(), "/summary/material") != 0 ||
	source->getRealizedDisplaySummary(4, displaySummary)) {
	printf("FAIL: database source should summarize realized display state\n");
	root->unref();
	return 1;
    }

    if (scene.getSceneDisplaySummaryCount() !=
	scene.getSceneTreeSummaryCount() ||
	!scene.getSceneDisplaySummary(0, displaySummary) ||
	!displaySummary.valid ||
	displaySummary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	displaySummary.hasDrawIntent ||
	displaySummary.ownerSourceIndex != -1 ||
	!scene.getSceneDisplaySummary(1, displaySummary) ||
	displaySummary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	displaySummary.hasDrawIntent ||
	displaySummary.ownerSourceIndex != -1 ||
	!scene.getSceneDisplaySummary(2, displaySummary) ||
	!displaySummary.isDatabaseSource ||
	displaySummary.ownerSourceIndex != 0 ||
	bu_strcmp(displaySummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	bu_strcmp(displaySummary.intentPath.getString(),
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
	displaySummary.nodeKind != BObolSceneTreeSummary::NODE_MESH_SHAPE ||
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
	BObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
	bu_strcmp(displaySummary.path.getString(), "/summary/material") != 0 ||
	scene.getSceneDisplaySummary(6, displaySummary)) {
	printf("FAIL: scene controller should summarize Obol scene display state\n");
	root->unref();
	return 1;
    }

    BObolSceneMaterialSummary materialStateSummary;
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
	BObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	!materialStateSummary.materialValid ||
	materialStateSummary.materialRevision != 0 ||
	materialStateSummary.materialColor[0] < 0.99f ||
	!source->getRealizedSceneMaterialSummary(2,
		materialStateSummary) ||
	materialStateSummary.nodeKind !=
	BObolSceneTreeSummary::NODE_MESH_SHAPE ||
	!materialStateSummary.materialValid ||
	materialStateSummary.materialRevision != 77 ||
	materialStateSummary.materialColor[0] < 0.24f ||
	materialStateSummary.materialColor[0] > 0.26f ||
	materialStateSummary.materialColor[2] < 0.74f ||
	materialStateSummary.materialColor[2] > 0.76f ||
	!source->getRealizedSceneMaterialSummary(3,
		materialStateSummary) ||
	materialStateSummary.nodeKind !=
	BObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
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
	BObolSceneTreeSummary::NODE_MESH_SHAPE ||
	materialStateSummary.ownerSourceIndex != 0 ||
	bu_strcmp(materialStateSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	!materialStateSummary.materialValid ||
	materialStateSummary.materialRevision != 77 ||
	materialStateSummary.materialColor[1] < 0.49f ||
	materialStateSummary.materialColor[1] > 0.51f ||
	!scene.getSceneMaterialSummary(5, materialStateSummary) ||
	materialStateSummary.nodeKind !=
	BObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
	materialStateSummary.materialValid ||
	scene.getSceneMaterialSummary(6, materialStateSummary)) {
	printf("FAIL: scene controller should summarize Obol shape material state\n");
	root->unref();
	return 1;
    }

    BObolSceneBoundsSummary boundsSummary;
    if (source->getRealizedBoundsSummaryCount() !=
	source->getRealizedTreeSummaryCount() ||
	!source->getRealizedBoundsSummary(0, boundsSummary) ||
	!boundsSummary.valid ||
	boundsSummary.nodeKind !=
	BObolSceneTreeSummary::NODE_DATABASE_SOURCE ||
	!boundsSummary.boundsValid ||
	boundsSummary.ownerSourceIndex != -1 ||
	bu_strcmp(boundsSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	boundsSummary.bounds.getMin()[0] > -1.74f ||
	boundsSummary.bounds.getMax()[0] < 1.74f ||
	boundsSummary.bounds.getMin()[1] > -1.74f ||
	boundsSummary.bounds.getMax()[1] < 1.74f ||
	!source->getRealizedBoundsSummary(1, boundsSummary) ||
	boundsSummary.nodeKind != BObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	!boundsSummary.boundsValid ||
	boundsSummary.bounds.getMin()[0] > -1.74f ||
	boundsSummary.bounds.getMax()[1] < 1.74f ||
	!source->getRealizedBoundsSummary(2, boundsSummary) ||
	boundsSummary.nodeKind != BObolSceneTreeSummary::NODE_MESH_SHAPE ||
	!boundsSummary.boundsValid ||
	boundsSummary.bounds.getMin()[0] < -0.01f ||
	boundsSummary.bounds.getMin()[0] > 0.01f ||
	boundsSummary.bounds.getMax()[0] < 0.99f ||
	boundsSummary.bounds.getMax()[1] < 0.99f ||
	!source->getRealizedBoundsSummary(3, boundsSummary) ||
	boundsSummary.nodeKind !=
	BObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
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
	boundsSummary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	!boundsSummary.boundsValid ||
	boundsSummary.ownerSourceIndex != -1 ||
	boundsSummary.bounds.getMin()[0] > -1.74f ||
	boundsSummary.bounds.getMax()[0] < 1.74f ||
	!scene.getSceneBoundsSummary(1, boundsSummary) ||
	boundsSummary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	boundsSummary.boundsValid ||
	!scene.getSceneBoundsSummary(2, boundsSummary) ||
	boundsSummary.nodeKind !=
	BObolSceneTreeSummary::NODE_DATABASE_SOURCE ||
	boundsSummary.ownerSourceIndex != 0 ||
	bu_strcmp(boundsSummary.ownerSourcePath.getString(),
	       "/summary/source") != 0 ||
	!boundsSummary.boundsValid ||
	!scene.getSceneBoundsSummary(4, boundsSummary) ||
	boundsSummary.nodeKind != BObolSceneTreeSummary::NODE_MESH_SHAPE ||
	boundsSummary.ownerSourceIndex != 0 ||
	boundsSummary.bounds.getMax()[0] < 0.99f ||
	boundsSummary.bounds.getMax()[1] < 0.99f ||
	!scene.getSceneBoundsSummary(5, boundsSummary) ||
	boundsSummary.nodeKind !=
	BObolSceneTreeSummary::NODE_MATERIAL_OBJECT ||
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
	summaryMesh->drawMode.getValue() != BOBOL_LOD_DRAW_SHADED) {
	printf("FAIL: scene controller should update source draw mode through the source owner\n");
	root->unref();
	return 1;
    }
    const uint64_t afterSourceShadedFrameRevision = scene.getFrameRevision();
    if (scene.setDatabaseSourceDrawMode("summary/source",
					SoBRLDatabaseSource::WIREFRAME) != 1 ||
	scene.getFrameRevision() <= afterSourceShadedFrameRevision ||
	source->drawMode.getValue() != SoBRLDatabaseSource::WIREFRAME ||
	summaryMesh->drawMode.getValue() != BOBOL_LOD_DRAW_WIRE ||
	!source->getRealizedShapeSummary(1, shapeSummary) ||
	shapeSummary.drawMode != BOBOL_LOD_DRAW_WIRE ||
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
    BObolSceneController groupScene(groupRoot);
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
	bu_strcmp(retainedLeaf->groupPath.getValue().getString(),
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
	treeSummary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	bu_strcmp(treeSummary.path.getString(), "/") != 0 ||
	treeSummary.childCount != 1 ||
	treeSummary.drawTreeDepth != 0 ||
	!groupScene.getSceneTreeSummary(1, treeSummary) ||
	treeSummary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	bu_strcmp(treeSummary.path.getString(), "assembly") != 0 ||
	treeSummary.childCount != 1 ||
	treeSummary.drawTreeDepth != 1 ||
	!groupScene.getSceneTreeSummary(2, treeSummary) ||
	treeSummary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	bu_strcmp(treeSummary.path.getString(), "assembly/leaf") != 0 ||
	treeSummary.childCount != 0 ||
	treeSummary.drawTreeDepth != 2 ||
	groupScene.getSceneTreeSummary(3, treeSummary)) {
	printf("FAIL: scene controller should summarize retained named group paths\n");
	groupRoot->unref();
	return 1;
    }

    if (!groupScene.getSceneDisplaySummary(2, displaySummary) ||
	displaySummary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	bu_strcmp(displaySummary.path.getString(), "assembly/leaf") != 0 ||
	!groupScene.getSceneBoundsSummary(2, boundsSummary) ||
	boundsSummary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	bu_strcmp(boundsSummary.path.getString(), "assembly/leaf") != 0) {
	printf("FAIL: scene controller display and bounds summaries should retain group paths\n");
	groupRoot->unref();
	return 1;
    }

    const uint64_t beforeGroupStateStructuralRevision =
	groupScene.getStructuralRevision();
    const uint64_t beforeGroupStateFrameRevision =
	groupScene.getFrameRevision();
    if (groupScene.setGroupDrawIntent("assembly/leaf",
				      "draw://assembly/leaf", BOBOL_LOD_DRAW_SHADED,
				      BOBOL_LOD_DRAW_WIRE, TRUE, 42) != 1 ||
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
	BOBOL_LOD_DRAW_WIRE ||
	retainedLeaf->revalidationRevision.getValue() != 42) {
	printf("FAIL: scene controller should retain group draw and display state in Obol\n");
	groupRoot->unref();
	return 1;
    }
    const uint64_t afterGroupStateFrameRevision =
	groupScene.getFrameRevision();
    if (groupScene.setGroupDrawIntent("assembly/leaf",
				      "draw://assembly/leaf", BOBOL_LOD_DRAW_SHADED,
				      BOBOL_LOD_DRAW_WIRE, TRUE, 42) != 0 ||
	groupScene.setGroupDisplayState("assembly/leaf", FALSE, TRUE,
					TRUE, 6, 7, 0.45f, TRUE,
					SbColor(0.8f, 0.2f, 0.1f), TRUE,
					SbColor(0.3f, 0.4f, 0.5f), 123) != 0 ||
	groupScene.getFrameRevision() != afterGroupStateFrameRevision ||
	groupScene.setGroupDrawIntent("missing", "missing",
				      BOBOL_LOD_DRAW_WIRE, BOBOL_LOD_DRAW_WIRE,
				      FALSE, 0) != -1) {
	printf("FAIL: scene controller should treat retained group state reapplication as a no-op\n");
	groupRoot->unref();
	return 1;
    }

    if (!groupScene.getSceneDisplaySummary(2, displaySummary) ||
	displaySummary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	bu_strcmp(displaySummary.path.getString(), "assembly/leaf") != 0 ||
	!displaySummary.hasDrawIntent ||
	bu_strcmp(displaySummary.intentPath.getString(),
	       "draw://assembly/leaf") != 0 ||
	displaySummary.intentDrawMode != BOBOL_LOD_DRAW_SHADED ||
	displaySummary.visible ||
	!displaySummary.highlighted ||
	displaySummary.lineStyle != 6 ||
	displaySummary.lineWidth != 7 ||
	fabs(displaySummary.transparency - 0.45) > 0.001 ||
	displaySummary.drawMode != BOBOL_LOD_DRAW_SHADED ||
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
	treeSummary.nodeKind != BObolSceneTreeSummary::NODE_MESH_SHAPE ||
	bu_strcmp(treeSummary.path.getString(), "assembly/leaf/mesh") != 0 ||
	treeSummary.drawTreeDepth != 3 ||
	!groupScene.getSceneBoundsSummary(2, boundsSummary) ||
	!boundsSummary.boundsValid ||
	boundsSummary.bounds.getMax()[0] < 0.99f ||
	boundsSummary.bounds.getMax()[1] < 0.99f ||
	!groupScene.getSceneBoundsSummary(3, boundsSummary) ||
	boundsSummary.nodeKind != BObolSceneTreeSummary::NODE_MESH_SHAPE ||
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
				     BOBOL_LOD_DRAW_POINTS, TRUE, FALSE, FALSE) != 1 ||
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
				     BOBOL_LOD_DRAW_POINTS, TRUE, FALSE, FALSE) != 0 ||
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
	displaySummary.nodeKind != BObolSceneTreeSummary::NODE_MESH_SHAPE ||
	bu_strcmp(displaySummary.path.getString(), "assembly/leaf/mesh") != 0 ||
	displaySummary.visible ||
	!displaySummary.highlighted ||
	displaySummary.lineStyle != 8 ||
	displaySummary.lineWidth != 9 ||
	fabs(displaySummary.transparency - 0.55) > 0.001 ||
	displaySummary.drawMode != BOBOL_LOD_DRAW_POINTS ||
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
	bu_strcmp(groupMesh->ownerSourcePath.getValue().getString(),
	       "/owner/source") != 0 ||
	groupMesh->ownerSourceRevision.getValue() != 10 ||
	groupMesh->ownerInputsRevision.getValue() != 11 ||
	groupMesh->ownerViewRevision.getValue() != 12 ||
	groupMesh->ownerRealizedRevision.getValue() != 13 ||
	groupMesh->ownerRealizedSourceRevision.getValue() != 14 ||
	groupMesh->ownerRealizedInputsRevision.getValue() != 15 ||
	groupMesh->ownerRealizedViewRevision.getValue() != 16 ||
	groupMesh->ownerRealizationStatus.getValue() != 17 ||
	bu_strcmp(groupMesh->ownerRealizationDiagnostic.getValue().getString(),
	       "stale source") != 0 ||
	bu_strcmp(groupMesh->ownerRealizationIdentity.getValue().getString(),
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
	bu_strcmp(treeSummary.ownerSourcePath.getString(),
	       "/owner/source") != 0 ||
	!groupScene.getSceneDisplaySummary(3, displaySummary) ||
	bu_strcmp(displaySummary.ownerSourcePath.getString(),
	       "/owner/source") != 0 ||
	!groupScene.getSceneBoundsSummary(3, boundsSummary) ||
	bu_strcmp(boundsSummary.ownerSourcePath.getString(),
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
	bu_strcmp(retainedLeaf->groupPath.getValue().getString(),
	       "assembly/renamed") != 0 ||
	bu_strcmp(retainedChild->groupPath.getValue().getString(),
	       "assembly/renamed/child") != 0 ||
	groupScene.getStructuralRevision() <=
	afterGroupChildStructuralRevision ||
	!groupScene.getSceneTreeSummary(2, treeSummary) ||
	bu_strcmp(treeSummary.path.getString(), "assembly/renamed") != 0 ||
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
	BObolViewController groupView(groupRoot);
	SbString renderReason;
	groupView.clearRenderRequest();
	SoGroup *viewLeaf = groupView.ensureGroup("view/leaf");
	if (!viewLeaf ||
	    groupView.findGroup("view/leaf") != viewLeaf ||
	    groupView.getGroupChildCount("view") != 1 ||
	    !groupView.consumeRenderRequest(&renderReason) ||
	    bu_strcmp(renderReason.getString(), "scene-group") != 0) {
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
					 "draw://view/leaf", BOBOL_LOD_DRAW_POINTS,
					 BOBOL_LOD_DRAW_WIRE, FALSE, 5) != 1 ||
	    !groupView.consumeRenderRequest(&renderReason) ||
	    bu_strcmp(renderReason.getString(), "scene-group") != 0) {
	    printf("FAIL: view controller should expose retained group draw intent state\n");
	    groupRoot->unref();
	    return 1;
	}
	groupView.clearRenderRequest();
	if (groupView.setGroupDrawIntent("view/leaf",
					 "draw://view/leaf", BOBOL_LOD_DRAW_POINTS,
					 BOBOL_LOD_DRAW_WIRE, FALSE, 5) != 0 ||
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
	    bu_strcmp(renderReason.getString(), "scene-group") != 0 ||
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
	    bu_strcmp(renderReason.getString(), "scene-group") != 0) {
	    printf("FAIL: view controller should expose retained group child ownership\n");
	    groupRoot->unref();
	    return 1;
	}
	groupView.clearRenderRequest();
	if (groupView.findShape("view/leaf/mesh") != viewMesh ||
	    groupView.findShapeParent("view/leaf/mesh") != viewLeaf ||
	    groupView.setShapeDrawState("view/leaf/mesh",
					BOBOL_LOD_DRAW_DIAGNOSTIC, FALSE, TRUE, FALSE) != 1 ||
	    !groupView.consumeRenderRequest(&renderReason) ||
	    bu_strcmp(renderReason.getString(), "scene-shape") != 0 ||
	    groupView.setShapeDrawState("view/leaf/mesh",
					BOBOL_LOD_DRAW_DIAGNOSTIC, FALSE, TRUE, FALSE) != 0) {
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
	    bu_strcmp(renderReason.getString(), "scene-shape") != 0 ||
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
	    bu_strcmp(renderReason.getString(), "scene-shape") != 0 ||
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
	    bu_strcmp(renderReason.getString(), "scene-shape") != 0 ||
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
	    bu_strcmp(renderReason.getString(), "scene-shape") != 0 ||
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
	    bu_strcmp(renderReason.getString(), "scene-group") != 0) {
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
    BObolSceneController ownedScene(ownedRoot);
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
	bu_strcmp(summary.path.getString(), "lod-submit.bot") != 0 ||
	!summary.hasParent ||
	summary.drawTreeDepth != 1 ||
	bu_strcmp(summary.parentGroupPath.getString(), "/") != 0 ||
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
	bu_strcmp(summary.parentGroupPath.getString(), "draw/group") != 0 ||
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
    BObolDatabaseSourceDisplayPatch sourcePatch;
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
	bu_strcmp(summary.displayName.getString(),
	       "LoD Submit Display") != 0 ||
	!ownedScene.getSceneTreeSummaryForPath("lod-submit.bot",
		treeSummary) ||
	bu_strcmp(treeSummary.displayName.getString(),
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
	!ownedSource->hasCompactInstanceIndex() ||
	!ownedSource->hasRealizedMeshGeometry() ||
	summary.realizedMeshCount != 0 ||
	summary.realizedShapeCount != 0) {
	printf("FAIL: mesh realization role should realize database mesh even in wire draw mode\n");
	ownedRoot->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    {
	BObolRealizedShapeSummary realizedMeshSummary;
	if (!ownedSource->getRealizedShapeSummary(0, realizedMeshSummary) ||
	    !realizedMeshSummary.valid ||
	    realizedMeshSummary.shapeKind !=
		BObolRealizedShapeSummary::SHAPE_MESH ||
	    realizedMeshSummary.pointCount == 0 ||
	    realizedMeshSummary.triangleCount == 0) {
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
	bu_strcmp(summary.parentGroupPath.getString(), "draw/group") != 0 ||
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
	bu_strcmp(summary.parentGroupPath.getString(), "/") != 0 ||
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

    BObolDatabaseSourceSummary sharedSummary;
    BObolDatabaseSourceSummary viewSummary;
    if (!find_database_source_summary_by_instance(ownedScene,
	    sharedInstanceKey, sharedSummary) ||
	!find_database_source_summary_by_instance(ownedScene,
		viewInstanceKey, viewSummary) ||
	bu_strcmp(sharedSummary.path.getString(), "lod-submit.bot") != 0 ||
	bu_strcmp(viewSummary.path.getString(), "lod-submit.bot") != 0 ||
	bu_strcmp(sharedSummary.instanceKey.getString(),
	       sharedInstanceKey) != 0 ||
	bu_strcmp(viewSummary.instanceKey.getString(), viewInstanceKey) != 0 ||
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
	bu_strcmp(viewSummary.parentGroupPath.getString(), "views/V0") != 0) {
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
	BObolRealizedShapeSummary instanceShapeSummary;
	if (!ownedScene.getRealizedShapeSummary(i, instanceShapeSummary) ||
	    !instanceShapeSummary.valid)
	    continue;
	const char *instanceShapePath =
	    instanceShapeSummary.path.getString();
	const int instancePathMatches =
	    bu_strcmp(instanceShapePath, "lod-submit.bot") == 0 ||
	    bu_strcmp(instanceShapePath, "/lod-submit.bot") == 0;
	if (bu_strcmp(instanceShapeSummary.ownerSourceInstanceKey.getString(),
		   sharedInstanceKey) == 0 &&
	    instancePathMatches &&
	    strstr(instanceShapeSummary.sourceIdentity.getString(),
		   sharedInstanceKey))
	    sharedRealized = 1;
	if (bu_strcmp(instanceShapeSummary.ownerSourceInstanceKey.getString(),
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
	bu_strcmp(viewSummary.path.getString(), "lod-submit.bot") != 0 ||
	bu_strcmp(viewSummary.parentGroupPath.getString(), "views/V0") != 0) {
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

    BObolLodService service;
    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD update action service did not start\n");
	root->unref();
	return 1;
    }

    BObolLodTask task;
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

    BObolViewLodState viewState;
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
    const BObolViewLodState::MeshPayload *payload =
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

    BObolLodRequest request;
    mesh->makeLodRequest(request,
			 "db://request-test",
			 12,
			 34,
			 56,
			 BOBOL_LOD_DRAW_SHADED,
			 "bobol_mesh_lod",
			 BOBOL_MESH_LOD_PROVIDER_VERSION,
			 BOBOL_LOD_QUALITY_FAST_DISPLAY);

    if (bu_strcmp(request.databaseId.getString(), "db://request-test") != 0 ||
	request.databaseRevision != 12 ||
	request.sourceRevision != 44 ||
	bu_strcmp(request.objectPath.getString(), "/mesh/request") != 0 ||
	bu_strcmp(request.objectName.getString(), "request") != 0 ||
	request.viewRevision != 34 ||
	request.policyRevision != 56 ||
	request.drawMode != BOBOL_LOD_DRAW_SHADED ||
	request.lodPolicy != 9 ||
	bu_strcmp(request.providerId.getString(), "bobol_mesh_lod") != 0 ||
	request.qualityTier != BOBOL_LOD_QUALITY_FAST_DISPLAY ||
	request.sourceCounts.faceCount != 1 ||
	request.sourceCounts.pointCount != 3 ||
	request.bounds.isEmpty()) {
	printf("FAIL: mesh LoD request helper did not preserve identity\n");
	mesh->unref();
	return 1;
    }

    BObolLodCacheKey firstKey = bobol_lod_cache_key(request);
    BObolLodRequest sameRequest;
    mesh->makeLodRequest(sameRequest,
			 "db://request-test",
			 12,
			 34,
			 56,
			 BOBOL_LOD_DRAW_SHADED,
			 "bobol_mesh_lod",
			 BOBOL_MESH_LOD_PROVIDER_VERSION,
			 BOBOL_LOD_QUALITY_FAST_DISPLAY);
    if (bu_strcmp(firstKey.value.getString(),
	       bobol_lod_cache_key(sameRequest).value.getString()) != 0) {
	printf("FAIL: mesh LoD request helper produced unstable cache key\n");
	mesh->unref();
	return 1;
    }

    SoSeparator *root = new SoSeparator;
    root->ref();
    SoOrthographicCamera *camera = new SoOrthographicCamera;
    camera->height = 15.0f;
    BObolViewController controller(root, camera);
    controller.setViewportSize(320, 240);
    struct bv_lod_policy controllerPolicy;
    bv_lod_policy_init(&controllerPolicy);
    controllerPolicy.scale = 3.5;
    controller.getViewAttachment()->setLodPolicy(&controllerPolicy);

    struct bv_view_info info = BV_VIEW_INFO_INIT;
    if (!controller.getViewInfo(&info) ||
	info.width != 320 ||
	info.height != 240 ||
	fabs(info.size - 20.0) > 1.0e-6 ||
	fabs(info.lod.scale - 3.5) > 1.0e-6) {
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

    /* A trace of -1 is the numerically difficult matrix-to-quaternion case.
     * It is also an ordinary CAD view (ae 0 90), so verify the complete
     * projection rather than just accepting a finite quaternion. */
    VSET(aet, 0.0, 90.0, 0.0);
    if (!bv_aet_set(view, aet) ||
	!bv_context_update(viewCtx, BV_CONTEXT_CHANGED_VIEW) ||
	!controller.syncCameraFromViewContext(viewCtx) ||
	!bv_model2view_get(model2view, view)) {
	printf("FAIL: could not prepare gimbal-position Obol camera sync\n");
	bv_context_destroy(viewCtx);
	root->unref();
	mesh->unref();
	return 1;
    }
    const SbViewVolume topViewVolume = syncedCamera->getViewVolume(0.0f);
    for (size_t i = 0; i < sizeof(samples) / sizeof(samples[0]); i++) {
	point_t brlView;
	MAT4X3PNT(brlView, model2view, samples[i]);
	SbVec3f obolScreen;
	topViewVolume.projectToScreen(SbVec3f(
					  static_cast<float>(samples[i][X]),
					  static_cast<float>(samples[i][Y]),
					  static_cast<float>(samples[i][Z])),
				      obolScreen);
	const double expectedX = 0.5 * (brlView[X] + 1.0);
	const double expectedY = 0.5 * (brlView[Y] * 2.0 + 1.0);
	if (fabs(obolScreen[0] - expectedX) > 1.0e-5 ||
	    fabs(obolScreen[1] - expectedY) > 1.0e-5) {
	    printf("FAIL: gimbal-position Obol projection mismatch for sample "
		   "%zu: BRL=(%.9f, %.9f) Obol=(%.9f, %.9f)\n",
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
    mesh->drawMode = BOBOL_LOD_DRAW_SHADED;
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
    vlist->drawMode = BOBOL_LOD_DRAW_WIRE;
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
    vlistShaded->drawMode = BOBOL_LOD_DRAW_SHADED;
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
    sourceBacked->drawMode = BOBOL_LOD_DRAW_SHADED;
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
    if (bu_strcmp(meshRecord.displayName.getString(), "Mesh Display") != 0 ||
	bu_strcmp(meshRecord.geometryName.getString(), "mesh geometry") != 0 ||
	bu_strcmp(meshRecord.cacheIdentity.getString(),
	       "cache://metadata/mesh") != 0 ||
	bu_strcmp(meshRecord.sourceIdentity.getString(),
	       "db://metadata/mesh.bot") != 0 ||
	meshRecord.cacheIdentityValue != meshCacheIdentity ||
	meshRecord.sourceIdentityValue != meshSourceIdentity ||
	!meshRecord.databaseIntent ||
	meshRecord.overlayIntent ||
	meshRecord.localSource ||
	meshRecord.sharedSource ||
	meshRecord.nonDatabaseSource ||
	meshRecord.drawMode != BOBOL_LOD_DRAW_SHADED ||
	bu_strcmp(meshRecord.recordRole.getString(), "database") != 0 ||
	bu_strcmp(meshRecord.geometryKind.getString(), "surface") != 0) {
	printf("FAIL: mesh export record did not preserve neutral metadata\n");
	root->unref();
	return 1;
    }

    const SoBRLExportAction::LineRecord &lineRecord = exportAction.getLine(0);
    const SoBRLExportAction::PointRecord &pointRecord =
	exportAction.getPoint(0);
    if (bu_strcmp(lineRecord.displayName.getString(), "Curve Display") != 0 ||
	bu_strcmp(lineRecord.geometryName.getString(), "curve geometry") != 0 ||
	bu_strcmp(lineRecord.geometryKind.getString(), "line") != 0 ||
	lineRecord.cacheIdentityValue != curveCacheIdentity ||
	lineRecord.sourceIdentityValue != curveSourceIdentity ||
	bu_strcmp(pointRecord.displayName.getString(), "Curve Display") != 0 ||
	bu_strcmp(pointRecord.geometryKind.getString(), "point") != 0 ||
	pointRecord.cacheIdentityValue != curveCacheIdentity ||
	pointRecord.sourceIdentityValue != curveSourceIdentity ||
	pointRecord.drawMode != BOBOL_LOD_DRAW_WIRE ||
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
    if (bu_strcmp(meshObject.path.getString(), "/metadata/mesh.bot") != 0 ||
	meshObject.triangleIndices.size() != 1 ||
	!meshObject.lineIndices.empty() ||
	!meshObject.pointIndices.empty() ||
	!meshObject.databaseIntent ||
	meshObject.sharedSource ||
	meshObject.cacheIdentityValue != meshCacheIdentity ||
	meshObject.sourceIdentityValue != meshSourceIdentity ||
	bu_strcmp(meshObject.geometryKind.getString(), "surface") != 0 ||
	meshObject.bounds.isEmpty()) {
	printf("FAIL: mesh object export record did not summarize triangle metadata\n");
	root->unref();
	return 1;
    }
    if (bu_strcmp(vlistObject.path.getString(), "/metadata/curve.sketch") != 0 ||
	vlistObject.lineIndices.size() != 1 ||
	vlistObject.pointIndices.size() != 1 ||
	!vlistObject.triangleIndices.empty() ||
	bu_strcmp(vlistObject.geometryKind.getString(), "mixed") != 0 ||
	vlistObject.cacheIdentityValue != curveCacheIdentity ||
	vlistObject.sourceIdentityValue != curveSourceIdentity ||
	vlistObject.drawMode != BOBOL_LOD_DRAW_WIRE) {
	printf("FAIL: vlist object export record did not group line and point records\n");
	root->unref();
	return 1;
    }
    if (bu_strcmp(vlistShadedObject.path.getString(),
	       "/metadata/curve.sketch") != 0 ||
	vlistShadedObject.lineIndices.size() != 1 ||
	!vlistShadedObject.pointIndices.empty() ||
	!vlistShadedObject.triangleIndices.empty() ||
	vlistShadedObject.drawMode != BOBOL_LOD_DRAW_SHADED ||
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
					  BOBOL_LOD_DRAW_WIRE) != 1 ||
	exportAction.collectObjectRecords(queriedObjects,
					  SoBRLExportAction::QUERY_DATABASE_OBJECTS, NULL,
					  BOBOL_LOD_DRAW_SHADED) != 2) {
	printf("FAIL: export object record filters did not match expected object sets\n");
	root->unref();
	return 1;
    }

    const BObolSourceMeshRequest &sourceRequest =
	exportAction.getSourceBackedFullDetailRequest(0);
    if (bu_strcmp(sourceRequest.displayName.getString(),
	       "Source Mesh Display") != 0 ||
	bu_strcmp(sourceRequest.geometryName.getString(),
	       "source mesh geometry") != 0 ||
	bu_strcmp(sourceRequest.cacheIdentity.getString(),
	       "cache://metadata/source") != 0 ||
	bu_strcmp(sourceRequest.sourceIdentity.getString(),
	       "db://metadata/source.bot") != 0 ||
	!sourceRequest.databaseIntent ||
	sourceRequest.sharedSource ||
	sourceRequest.drawMode != BOBOL_LOD_DRAW_SHADED ||
	bu_strcmp(sourceRequest.recordRole.getString(), "database") != 0 ||
	bu_strcmp(sourceRequest.geometryKind.getString(), "surface") != 0) {
	printf("FAIL: source-backed request did not preserve neutral metadata\n");
	root->unref();
	return 1;
    }

    BObolLodRequest sourceLodRequest;
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
    if (bu_strcmp(sourceRecord.displayName.getString(),
	       "Source Mesh Display") != 0 ||
	bu_strcmp(sourceRecord.cacheIdentity.getString(),
	       "cache://metadata/source") != 0 ||
	bu_strcmp(sourceRecord.sourceIdentity.getString(),
	       "db://metadata/source.bot") != 0 ||
	sourceRecord.cacheIdentityValue != sourceCacheIdentity ||
	sourceRecord.sourceIdentityValue != sourceSourceIdentity ||
	!sourceRecord.databaseIntent ||
	sourceRecord.sharedSource ||
	sourceRecord.nonDatabaseSource ||
	sourceRecord.drawMode != BOBOL_LOD_DRAW_SHADED ||
	bu_strcmp(sourceRecord.recordRole.getString(), "database") != 0 ||
	bu_strcmp(sourceRecord.geometryKind.getString(), "surface") != 0) {
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
    if (bu_strcmp(sourceObject.path.getString(), "/metadata/source.bot") != 0 ||
	sourceObject.triangleIndices.size() != 1 ||
	!sourceObject.databaseIntent ||
	sourceObject.sharedSource ||
	sourceObject.nonDatabaseSource ||
	bu_strcmp(sourceObject.displayName.getString(),
	       "Source Mesh Display") != 0 ||
	sourceObject.cacheIdentityValue != sourceCacheIdentity ||
	sourceObject.sourceIdentityValue != sourceSourceIdentity ||
	bu_strcmp(sourceObject.cacheIdentity.getString(),
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

    BObolLodRequest preservingRequest =
	make_request("/budget/preserve.bot", "preserve.bot");
    BObolLodRequest lodRequest =
	make_request("/budget/lod.bot", "lod.bot");
    BObolLodResult preservingResult =
	mesh_payload_result(preservingRequest);
    BObolLodResult lodResult =
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
	BObolViewController controller(root, camera);
	controller.clearRenderRequest();
	size_t freedBytes = controller.evictMeshPayloadsToBudget(0, TRUE);
	BObolSourceMeshRequest preservingSourceRequest;
	BObolSourceMeshRequest lodSourceRequest;
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
	    bu_strcmp(controller.getRenderReason().getString(),
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
	BObolViewController controller(root, camera);
	BObolLodService service;
	uint64_t taskGeneration = 0;

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
	    BObolLodTask task;
	    task.generation = service.beginGeneration();
	    taskGeneration = task.generation;
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
	    int applied = controller.applyLodResults(
		&service, 0, 0, taskGeneration);
	    BObolSourceMeshRequest sourceRequest;
	    if (applied != 1 ||
		controller.getLastLodAppliedResultCount() != 1 ||
		controller.getLastMeshBudgetVisitedMeshCount() != 1 ||
		controller.getLastMeshBudgetInitialResidentBytes() == 0 ||
		controller.getLastMeshBudgetFinalResidentBytes() != 0 ||
		controller.getLastMeshBudgetEvictedFullDetailMeshCount() != 1 ||
		controller.getLastMeshBudgetEvictedDisplayMeshCount() != 1 ||
		!controller.isRenderRequested() ||
		bu_strcmp(controller.getRenderReason().getString(),
		       "lod-memory-budget") != 0 ||
		mesh->isLodDisplayActive() ||
		mesh->hasFullDetailMesh() ||
		mesh->getTriangleCount() != 0 ||
		mesh->estimateResidentMeshBytes() != 0 ||
		!mesh->makeSourceMeshRequest(sourceRequest) ||
		check_source_mesh_request(sourceRequest,
					  "/budget/auto.bot", "auto.bot", 131)) {
		printf("FAIL: controller mesh residency budget did not auto-evict "
		       "applied LoD payload (applied=%d matched=%u initial=%zu "
		       "final=%zu full=%u display=%u render=%d reason=%s "
		       "lod=%d exact=%d triangles=%d resident=%zu)\n",
		       applied, controller.getLastLodAppliedResultCount(),
		       controller.getLastMeshBudgetInitialResidentBytes(),
		       controller.getLastMeshBudgetFinalResidentBytes(),
		       controller.getLastMeshBudgetEvictedFullDetailMeshCount(),
		       controller.getLastMeshBudgetEvictedDisplayMeshCount(),
		       controller.isRenderRequested() ? 1 : 0,
		       controller.getRenderReason().getString(),
		       mesh->isLodDisplayActive() ? 1 : 0,
		       mesh->hasFullDetailMesh() ? 1 : 0,
		       mesh->getTriangleCount(),
		       mesh->estimateResidentMeshBytes());
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
	   "bobol_lod_submit_action_cache", NULL);
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

    BObolLodService service;
    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD submit action service did not start\n");
	root->unref();
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    struct bv_view_info view = BV_VIEW_INFO_INIT;
    view.width = 640;
    view.height = 480;
    view.size = 100.0;

    BObolLodRequest activeDuplicateRequest;
    mesh->makeLodRequest(activeDuplicateRequest,
			 "db://lod-submit-test",
			 2026,
			 55,
			 66,
			 BOBOL_LOD_DRAW_SHADED,
			 "bobol_mesh_lod",
			 BOBOL_MESH_LOD_PROVIDER_VERSION,
			 BOBOL_LOD_QUALITY_FAST_DISPLAY);
    BObolLodTask activeTask;
    activeTask.generation = service.beginGeneration();
    activeTask.request = activeDuplicateRequest;
    activeTask.realize = aabb_task;
    activeTask.debugDelayMilliseconds = 200;
    if (service.submit(activeTask) == 0) {
	printf("FAIL: LoD submit action active-duplicate fixture did not queue request\n");
	service.stop();
	root->unref();
	bobol_mesh_lod_cache_clear_database(dbip);
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
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }
    if (wait_for_service(service)) {
	service.stop();
	root->unref();
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }
    {
	std::vector<BObolLodResult> activeResults;
	service.drainResults(activeResults);
	if (activeResults.size() != 1 ||
	    !bobol_lod_result_matches_request(activeResults[0],
						activeDuplicateRequest)) {
	    printf("FAIL: LoD submit action active-duplicate fixture did not preserve queued request\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
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
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    if (wait_for_service(service)) {
	service.stop();
	root->unref();
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    std::vector<BObolLodResult> viewPolicyResults;
    if (service.drainResults(viewPolicyResults) != 1 ||
	viewPolicyResults.size() != 1) {
	printf("FAIL: LoD submit action result drain failed\n");
	service.stop();
	root->unref();
	bobol_mesh_lod_cache_clear_database(dbip);
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
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    /*
     * A newly resident asset publishes its minimum useful mesh first.  This
     * gives every leaf coverage before any one leaf consumes optional suffix
     * capacity.  The unchanged view request must then refine directly to the
     * view-policy target on its next pass.
     */
    if (!viewPolicyResults[0].progressiveMesh ||
	viewPolicyResults[0].geometry.activeLevel !=
	    viewPolicyResults[0].progressiveMesh->minimumLevel() ||
	viewPolicyResults[0].terminal) {
	printf("FAIL: LoD submit action did not publish coverage-first minimum\n");
	service.stop();
	root->unref();
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    SoBRLMeshLodSubmitAction refineSubmit;
    refineSubmit.setService(&service);
    refineSubmit.setDatabase(dbip, "db://lod-submit-test", 2026);
    refineSubmit.setViewInfo(&view);
    refineSubmit.setGeneration(service.beginGeneration());
    refineSubmit.setRevisions(55, 66);
    refineSubmit.apply(root);
    if (refineSubmit.getSubmittedTaskCount() != 1 ||
	wait_for_service(service)) {
	printf("FAIL: LoD submit action did not queue coverage refinement\n");
	service.stop();
	root->unref();
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }
    std::vector<BObolLodResult> refinedViewPolicyResults;
    if (service.drainResults(refinedViewPolicyResults) != 1 ||
	refinedViewPolicyResults.size() != 1 ||
	refinedViewPolicyResults[0].geometry.activeLevel != expectedViewLevel) {
	printf("FAIL: LoD submit action did not use view-policy refinement level\n");
	service.stop();
	root->unref();
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    /* Projected demand, not whole-view model size, selects interactive PoP
     * levels.  Settled demand must refine beyond motion demand, and an
     * off-frustum mesh must submit no payload at all. */
    SoOrthographicCamera *projectedCamera = new SoOrthographicCamera;
    projectedCamera->ref();
    projectedCamera->position = SbVec3f(0.5f, 0.5f, 5.0f);
    projectedCamera->height = 10.0f;
    const SbViewVolume projectedVolume =
	projectedCamera->getViewVolume(640.0f / 480.0f);
    int motionLevel = -1;
    int settledLevel = -1;
    for (int pass = 0; pass < 2; pass++) {
	SoBRLMeshLodSubmitAction projectedSubmit;
	projectedSubmit.setService(&service);
	projectedSubmit.setDatabase(dbip, "db://lod-submit-test", 2026);
	projectedSubmit.setViewInfo(&view);
	projectedSubmit.setViewVolume(&projectedVolume, pass ? 1.0f : 4.0f);
	projectedSubmit.setGeneration(service.beginGeneration());
	projectedSubmit.setRevisions(70 + pass, 80 + pass);
	projectedSubmit.apply(root);
	if (projectedSubmit.getSubmittedTaskCount() != 1 ||
	    wait_for_service(service)) {
	    printf("FAIL: projected LoD demand did not submit visible mesh\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    projectedCamera->unref();
	    return 1;
	}
	std::vector<BObolLodResult> projectedResults;
	if (service.drainResults(projectedResults) != 1 ||
	    projectedResults.size() != 1 ||
	    projectedResults[0].request.requestedLevel < 0) {
	    printf("FAIL: projected LoD demand result was not view-quantized\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    projectedCamera->unref();
	    return 1;
	}
	if (pass)
	    settledLevel = projectedResults[0].request.requestedLevel;
	else
	    motionLevel = projectedResults[0].request.requestedLevel;
    }
    if (settledLevel <= motionLevel) {
	printf("FAIL: settled projected LoD demand did not refine motion level\n");
	service.stop();
	root->unref();
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	projectedCamera->unref();
	return 1;
    }

    projectedCamera->position = SbVec3f(1000.0f, 0.0f, 5.0f);
    const SbViewVolume offscreenVolume =
	projectedCamera->getViewVolume(640.0f / 480.0f);
    SoBRLMeshLodSubmitAction offscreenSubmit;
    offscreenSubmit.setService(&service);
    offscreenSubmit.setDatabase(dbip, "db://lod-submit-test", 2026);
    offscreenSubmit.setViewInfo(&view);
    offscreenSubmit.setViewVolume(&offscreenVolume, 1.0f);
    offscreenSubmit.setGeneration(service.beginGeneration());
    offscreenSubmit.setRevisions(90, 91);
    offscreenSubmit.apply(root);
    if (offscreenSubmit.getSubmittedTaskCount() != 0 ||
	offscreenSubmit.getSkippedMeshCount() != 2) {
	printf("FAIL: off-frustum LoD demand submitted display geometry\n");
	service.stop();
	root->unref();
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	projectedCamera->unref();
	return 1;
    }
    projectedCamera->unref();

    BObolViewLodState viewState;
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
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    update.apply(root);
    int ret = 0;
    const BObolViewLodState::MeshPayload *activePayload =
	viewState.findMesh(mesh);
    if (update.getAppliedResultCount() != 1 ||
	update.getRejectedResultCount() != 0 ||
	!activePayload ||
	activePayload->resultKind != BOBOL_LOD_RESULT_MESH ||
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
	BObolLodRequest activeRequest;
	mesh->makeLodRequest(activeRequest,
			     "db://lod-submit-test",
			     2026,
			     55,
			     66,
			     BOBOL_LOD_DRAW_SHADED,
			     "bobol_mesh_lod",
			     BOBOL_MESH_LOD_PROVIDER_VERSION,
			     BOBOL_LOD_QUALITY_FAST_DISPLAY);
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
	    BObolLodRequest sourceLodRequest;
	    if (!exactExport.makeSourceBackedFullDetailLodRequest(0,
		    sourceLodRequest) ||
		check_source_full_detail_lod_request(sourceLodRequest,
			"/lod-submit.bot", "lod-submit.bot")) {
		printf("FAIL: exact export source-backed request did not convert to RT full-detail LoD request\n");
		ret = 1;
	    } else {
		int beforeTriangleCount = exactExport.getTriangleCount();
		BObolLodResult sourceResult =
		    source_full_detail_result(sourceLodRequest);
		std::vector<BObolLodResult> sourceResults;
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
		    BObolLodRequest controllerExportRequest;
		    if (controllerExport.getSourceBackedFullDetailRequestCount() != 1 ||
			!controllerExport.makeSourceBackedFullDetailLodRequest(0,
				controllerExportRequest)) {
			printf("FAIL: controller source-backed exact export helper did not collect source request\n");
			ret = 1;
		    } else {
			BObolLodService controllerExportService;
			if (!controllerExportService.start(1, TRUE)) {
			    printf("FAIL: controller source-backed exact export helper service did not start\n");
			    ret = 1;
			} else {
			    if (submit_source_full_detail_task(
				    controllerExportService,
				    controllerExportRequest)) {
				ret = 1;
			    } else {
				BObolViewController exportController(root,
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
		BObolLodService sourceService;
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
			std::vector<BObolLodResult> submittedResults;
			sourceService.drainResults(submittedResults);
			if (submittedResults.size() != 1 ||
			    submittedResults[0].providerStatus !=
			    BOBOL_LOD_PROVIDER_STALE) {
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
	    const BObolSourceMeshRequest &measureSourceRequest =
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
	    BObolLodRequest sourceLodRequest;
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
		BObolLodResult sourceResult =
		    source_full_detail_result(sourceLodRequest);
		std::vector<BObolLodResult> sourceResults;
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
		    BObolLodRequest controllerMeasureRequest;
		    if (controllerMeasure.getSourceBackedFullDetailRequestCount() != 1 ||
			!controllerMeasure.makeSourceBackedFullDetailLodRequest(0,
				controllerMeasureRequest)) {
			printf("FAIL: controller source-backed exact measure helper did not collect source request\n");
			ret = 1;
		    } else {
			BObolLodService controllerMeasureService;
			if (!controllerMeasureService.start(1, TRUE)) {
			    printf("FAIL: controller source-backed exact measure helper service did not start\n");
			    ret = 1;
			} else {
			    if (submit_source_full_detail_task(
				    controllerMeasureService,
				    controllerMeasureRequest)) {
				ret = 1;
			    } else {
				BObolViewController measureController(root,
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
	    const BObolSourceMeshRequest &boundedMeasureRequest =
		boundedMeasure.getSourceBackedFullDetailRequest(0);
	    if (!boundedMeasureRequest.queryBoundsValid ||
		boundedMeasureRequest.queryBounds.isEmpty() ||
		!boundedMeasureRequest.queryToleranceValid ||
		boundedMeasureRequest.queryTolerance <= 0.0f) {
		printf("FAIL: bounded exact measure source-backed request did not carry explicit query tolerance\n");
		ret = 1;
	    } else {
		BObolLodRequest boundedMeasureLodRequest;
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
	    const BObolSourceMeshRequest &snapSourceRequest =
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
	    BObolLodRequest sourceLodRequest;
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
		BObolLodResult sourceResult =
		    source_full_detail_result(sourceLodRequest);
		std::vector<BObolLodResult> sourceResults;
		sourceResults.push_back(sourceResult);
		exactSnap.setQueryPoint(SbVec3f(0.2f, 0.2f, 0.0f));
		if (exactSnap.consumeSourceBackedFullDetailResults(
			sourceResults) != 1 ||
		    !exactSnap.hasCandidate() ||
		    exactSnap.getKind() != SoBRLSnapAction::FACE_NEAREST ||
		    exactSnap.getPrimitiveIndex() != 0 ||
		    bu_strcmp(exactSnap.getPath().getString(),
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
		    BObolLodRequest controllerSnapRequest;
		    if (controllerSnap.getSourceBackedFullDetailRequestCount() != 1 ||
			!controllerSnap.makeSourceBackedFullDetailLodRequest(0,
				controllerSnapRequest)) {
			printf("FAIL: controller source-backed exact snap helper did not collect source request\n");
			ret = 1;
		    } else {
			BObolLodService controllerSnapService;
			if (!controllerSnapService.start(1, TRUE)) {
			    printf("FAIL: controller source-backed exact snap helper service did not start\n");
			    ret = 1;
			} else {
			    if (submit_source_full_detail_task(
				    controllerSnapService,
				    controllerSnapRequest)) {
				ret = 1;
			    } else {
				BObolViewController snapController(root,
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
				    bu_strcmp(controllerSnap.getPath().getString(),
					   "/lod-submit.bot") != 0) {
				    printf("FAIL: controller source-backed exact snap helper did not consume matching LoD result\n");
				    ret = 1;
				}
			    }
			    controllerSnapService.stop();
			}
		    }
		}
		BObolSourceMeshPickResult sourcePick;
		if (!bobol_pick_source_full_detail_result(sourcePick,
			exactSnap.getSourceBackedFullDetailRequest(0),
			sourceResult,
			SbVec3f(0.2f, 0.2f, 5.0f),
			SbVec3f(0.0f, 0.0f, -1.0f)) ||
		    !sourcePick.hit ||
		    bu_strcmp(sourcePick.detail.getPath().getString(),
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
		    BObolLodRequest pickLodRequest;
		    BObolSourceMeshPickResult actionPick;
		    const BObolSourceMeshRequest &pickSourceRequest =
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
			std::vector<BObolLodResult> pickSourceResults;
			pickSourceResults.push_back(
			    source_full_detail_result(pickLodRequest));
			if (sourcePickAction.consumeSourceBackedFullDetailResults(
				actionPick, pickSourceResults) != 1 ||
			    !actionPick.hit ||
			    bu_strcmp(actionPick.detail.getPath().getString(),
				   "/lod-submit.bot") != 0 ||
			    actionPick.detail.getPrimitiveIndex() != 0) {
			    printf("FAIL: exact pick action did not collect and consume source-backed full-detail LoD result\n");
			    ret = 1;
			}
			if (!ret) {
			    BObolLodService controllerPickService;
			    if (!controllerPickService.start(1, TRUE)) {
				printf("FAIL: controller source-backed exact pick helper service did not start\n");
				ret = 1;
			    } else {
				BObolLodTask controllerPickTask;
				controllerPickTask.generation =
				    controllerPickService.beginGeneration();
				controllerPickTask.request = pickLodRequest;
				controllerPickTask.realize = source_full_detail_task;
				if (controllerPickService.submit(controllerPickTask) == 0 ||
				    wait_for_service(controllerPickService)) {
				    ret = 1;
				} else {
				    BObolViewController pickController(root, NULL);
				    pickController.setLodService(&controllerPickService);
				    BObolSourceMeshPickResult controllerPick;
				    int submittedCount = -1;
				    if (pickController.pickSourceMeshExactRay(
					    controllerPick,
					    SbVec3f(0.2f, 0.2f, 5.0f),
					    SbVec3f(0.0f, 0.0f, -1.0f),
					    0, &submittedCount) != 1 ||
					submittedCount != 0 ||
					!controllerPick.hit ||
					bu_strcmp(controllerPick.detail.getPath().getString(),
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
		    BObolSourceMeshRequest mappedSourceRequest =
			exactSnap.getSourceBackedFullDetailRequest(0);
		    mappedSourceRequest.faceCount = 0;

		    BObolLodResult mappedSourceResult =
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

		    BObolSourceMeshRequest compactMeasureRequest =
			mappedSourceRequest;
		    compactMeasureRequest.faceCount = 2;
		    compactMeasureRequest.pointCount = 6;
		    compactMeasureRequest.queryBoundsValid = 1;
		    compactMeasureRequest.queryBounds = SbBox3f(
							    SbVec3f(-0.1f, -0.1f, -0.1f),
							    SbVec3f(0.2f, 0.2f, 0.2f));
		    compactMeasureRequest.queryToleranceValid = 1;
		    compactMeasureRequest.queryTolerance = 0.5f;
		    BObolLodResult compactMeasureResult =
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
		    BObolLodResult compactMeasureNoFaceMap =
			compactMeasureResult;
		    compactMeasureNoFaceMap.mesh.faceIndex.clear();
		    SoBRLMeasureAction rejectedCompactMeasureFaceMap;
		    if (rejectedCompactMeasureFaceMap.consumeSourceBackedFullDetailResult(
			    compactMeasureRequest, compactMeasureNoFaceMap)) {
			printf("FAIL: exact measure accepted compact bounded source subset without face index mapping\n");
			ret = 1;
		    }
		    BObolLodResult compactMeasureNoVertexMap =
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
		    BObolSourceMeshRequest compactSnapRequest =
			mappedSourceRequest;
		    compactSnapRequest.faceCount = 2;
		    compactSnapRequest.pointCount = 6;
		    compactSnapRequest.queryBoundsValid = 1;
		    compactSnapRequest.queryBounds = SbBox3f(
							 SbVec3f(-0.1f, -0.1f, -0.1f),
							 SbVec3f(0.2f, 0.2f, 0.2f));
		    compactSnapRequest.queryToleranceValid = 1;
		    compactSnapRequest.queryTolerance = 0.5f;
		    BObolLodResult compactSnapResult =
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
		    BObolLodResult compactSnapNoVertexMap =
			compactSnapResult;
		    compactSnapNoVertexMap.mesh.vertexIndex.clear();
		    SoBRLSnapAction rejectedCompactSnap;
		    if (rejectedCompactSnap.consumeSourceBackedFullDetailResult(
			    compactSnapRequest, compactSnapNoVertexMap)) {
			printf("FAIL: exact snap accepted compact source vertex subset without vertex index mapping\n");
			ret = 1;
		    }

		    BObolSourceMeshPickResult mappedPick;
		    if (!bobol_pick_source_full_detail_result(mappedPick,
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

		    BObolSourceMeshRequest rayScopedPickRequest =
			mappedSourceRequest;
		    rayScopedPickRequest.faceCount = 4;
		    rayScopedPickRequest.pointCount = 6;
		    rayScopedPickRequest.queryRayValid = 1;
		    rayScopedPickRequest.queryRayDirection =
			SbVec3f(0.0f, 0.0f, -1.0f);
		    BObolLodResult rayScopedPickResult =
			source_full_detail_result(sourceLodRequest);
		    rayScopedPickResult.mesh.faceIndex.push_back(2);
		    rayScopedPickResult.mesh.vertexIndex.push_back(20);
		    rayScopedPickResult.mesh.vertexIndex.push_back(21);
		    rayScopedPickResult.mesh.vertexIndex.push_back(22);
		    BObolSourceMeshPickResult rayScopedPick;
		    if (!bobol_pick_source_full_detail_result(rayScopedPick,
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
		    BObolLodResult rayScopedPickNoVertexMap =
			rayScopedPickResult;
		    rayScopedPickNoVertexMap.mesh.vertexIndex.clear();
		    if (bobol_pick_source_full_detail_result(rayScopedPick,
			    rayScopedPickRequest, rayScopedPickNoVertexMap,
			    SbVec3f(0.2f, 0.2f, 5.0f),
			    SbVec3f(0.0f, 0.0f, -1.0f))) {
			printf("FAIL: exact pick accepted compact ray-scoped source point subset without vertex index mapping\n");
			ret = 1;
		    }
		    rayScopedPickRequest.queryRayValid = 0;
		    if (bobol_pick_source_full_detail_result(rayScopedPick,
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
	BObolLodRequest evictedDisplayPickRequest;
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
	    std::vector<BObolLodResult> forcedResults;
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
		const BObolViewLodState::MeshPayload *forcedPayload =
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
    bobol_mesh_lod_cache_clear_database(dbip);
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
    BObolLodRequest displayRequest =
	make_request("/lod-submit.bot", "lod-submit.bot");
    BObolLodResult displayResult = mesh_payload_result(displayRequest);
    if (!mesh->applyStagedLodResult(displayResult, &displayRequest) ||
	!mesh->isLodDisplayActive()) {
	printf("FAIL: controller multi-source exact submit fixture did not stage LoD mesh\n");
	root->unref();
	bobol_mesh_lod_cache_clear_database(right_dbip);
	db_close(right_dbip);
	db_close(wrong_dbip);
	bu_file_delete(right_dbpath);
	bu_file_delete(wrong_dbpath);
	return 1;
    }

    rightSource->addChild(mesh);
    root->addChild(wrongSource);
    root->addChild(rightSource);

    BObolLodService service;
    if (!service.start(1, TRUE)) {
	printf("FAIL: controller multi-source source-backed exact service did not start\n");
	root->unref();
	bobol_mesh_lod_cache_clear_database(right_dbip);
	db_close(right_dbip);
	db_close(wrong_dbip);
	bu_file_delete(right_dbpath);
	bu_file_delete(wrong_dbpath);
	return 1;
    }

    int ret = 0;
    BObolViewController controller(root, NULL);
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
	std::vector<BObolLodResult> submittedResults;
	service.drainResults(submittedResults);
	const char *rightDbId = right_dbip->dbi_filename ?
				right_dbip->dbi_filename : "";
	if (submittedResults.size() != 1 ||
	    submittedResults[0].providerStatus !=
	    BOBOL_LOD_PROVIDER_STALE ||
	    bu_strcmp(submittedResults[0].request.databaseId.getString(),
		   rightDbId) != 0 ||
	    submittedResults[0].request.sourceRevision != 222 ||
	    bu_strcmp(submittedResults[0].request.objectPath.getString(),
		   "/lod-submit.bot") != 0) {
	    printf("FAIL: controller multi-source source-backed exact submit did not use matching database source\n");
	    ret = 1;
	} else {
	    BObolLodTask task;
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
    bobol_mesh_lod_cache_clear_database(right_dbip);
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
    BObolLodRequest readyDisplayRequest =
	make_request("/lod-submit.bot", "lod-submit.bot");
    BObolLodResult readyDisplayResult =
	mesh_payload_result(readyDisplayRequest);
    if (!readyMesh->applyStagedLodResult(readyDisplayResult,
					 &readyDisplayRequest) ||
	!readyMesh->isLodDisplayActive()) {
	printf("FAIL: controller partial-ready exact fixture did not stage ready LoD mesh\n");
	root->unref();
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    SoBRLLodMeshShape *missingMesh = make_lod_mesh(
					 "/lod-submit-missing.bot", "lod-submit-missing.bot");
    missingMesh->sourceId = 333;
    BObolLodRequest missingDisplayRequest =
	make_request("/lod-submit-missing.bot", "lod-submit-missing.bot");
    BObolLodResult missingDisplayResult =
	mesh_payload_result(missingDisplayRequest);
    if (!missingMesh->applyStagedLodResult(missingDisplayResult,
					   &missingDisplayRequest) ||
	!missingMesh->isLodDisplayActive()) {
	printf("FAIL: controller partial-ready exact fixture did not stage missing LoD mesh\n");
	root->unref();
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    source->addChild(readyMesh);
    source->addChild(missingMesh);
    root->addChild(source);

    BObolLodService service;
    if (!service.start(1, TRUE)) {
	printf("FAIL: controller partial-ready source-backed exact service did not start\n");
	root->unref();
	bobol_mesh_lod_cache_clear_database(dbip);
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

    BObolLodRequest readySourceRequest;
    if (!ret && !exportAction.makeSourceBackedFullDetailLodRequest(0,
	    readySourceRequest)) {
	printf("FAIL: controller partial-ready source-backed exact test did not convert ready request\n");
	ret = 1;
    }

    if (!ret) {
	BObolLodTask readyTask;
	readyTask.generation = service.beginGeneration();
	readyTask.request = readySourceRequest;
	readyTask.realize = source_full_detail_task;
	if (service.submit(readyTask) == 0 || wait_for_service(service))
	    ret = 1;
    }

    if (!ret) {
	BObolViewController controller(root, NULL);
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

	BObolLodService pickService;
	if (!pickService.start(1, TRUE)) {
	    printf("FAIL: controller partial-ready exact pick service did not start\n");
	    ret = 1;
	} else {
	    const SbVec3f rayOrigin(0.2f, 0.2f, 5.0f);
	    const SbVec3f rayDirection(0.0f, 0.0f, -1.0f);
	    SoBRLSourceMeshPickAction requestAction;
	    requestAction.setRay(rayOrigin, rayDirection);
	    requestAction.apply(root);

	    BObolLodRequest readyPickRequest;
	    if (requestAction.getSourceBackedFullDetailRequestCount() != 2 ||
		!requestAction.makeSourceBackedFullDetailLodRequest(0,
			readyPickRequest)) {
		printf("FAIL: controller partial-ready exact pick did not collect two source requests\n");
		ret = 1;
	    } else {
		BObolLodTask readyPickTask;
		readyPickTask.generation = pickService.beginGeneration();
		readyPickTask.request = readyPickRequest;
		readyPickTask.realize = source_full_detail_task;
		if (pickService.submit(readyPickTask) == 0 ||
		    wait_for_service(pickService)) {
		    ret = 1;
		} else {
		    BObolViewController pickController(root, NULL);
		    pickController.setLodService(&pickService);
		    BObolSourceMeshPickResult pick;
		    int submittedCount = -1;
		    if (pickController.pickSourceMeshExactRay(pick,
			    rayOrigin, rayDirection,
			    pickService.beginGeneration(),
			    &submittedCount) != 1 ||
			submittedCount != 1 ||
			!pick.hit ||
			bu_strcmp(pick.detail.getPath().getString(),
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
    bobol_mesh_lod_cache_clear_database(dbip);
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
	   "bobol_lod_view_controller_cache", NULL);
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

    BObolLodService service;
    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD view-controller service did not start\n");
	root->unref();
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    {
	std::atomic<unsigned int> frameRequestCount(0);
	BObolViewController controller(root, camera);
	controller.setFrameRequestCallback(
	    count_frame_request, &frameRequestCount);
	controller.setLodService(&service);
	controller.setViewportSize(800, 600);
	controller.setLodPolicyRevision(99);
	controller.clearRenderRequest();

	uint64_t viewRevision = controller.getLodViewRevision();
	if (viewRevision == 0 || controller.getLodPolicyRevision() != 99) {
	    printf("FAIL: LoD view-controller revisions were not initialized\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
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
	    controller.isRenderRequested()) {
	    printf("FAIL: LoD view controller did not submit expected task\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	if (wait_for_service(service)) {
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	if (!controller.hasPendingLodResults() ||
	    !controller.hasProgressiveWorkPending() ||
	    controller.isRenderRequested()) {
	    printf("FAIL: LoD view controller did not observe result-ready "
		   "wakeup (results=%d progressive=%d render=%d reason=%s)\n",
		   controller.hasPendingLodResults() ? 1 : 0,
		   controller.hasProgressiveWorkPending() ? 1 : 0,
		   controller.isRenderRequested() ? 1 : 0,
		   controller.getRenderReason().getString());
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
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
	    !controller.getViewLodState()->findMesh(mesh) ||
	    !controller.isRenderRequested() ||
	    bu_strcmp(controller.getRenderReason().getString(),
		   "lod-result") != 0 ||
	    controller.hasPendingLodResults()) {
	    printf("FAIL: LoD view controller did not apply service result\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	controller.clearRenderRequest();
	controller.setLodAutoSubmit(TRUE);
	if (!controller.isLodAutoSubmitEnabled() ||
	    !controller.isRenderRequested() ||
	    bu_strcmp(controller.getRenderReason().getString(),
		   "lod-auto-submit") != 0) {
	    printf("FAIL: LoD view controller did not enable auto submit\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	/*
	 * Present the coverage-first minimum, then let the normal progressive
	 * pump finish this one-mesh view.  Camera-retarget assertions below are
	 * stable-view tests; leaving the first publication pending would instead
	 * (correctly) require provider growth on the later zoom.
	 */
	BObolProgressiveOptions settleOptions;
	settleOptions.forceTerminalLodRefinement = TRUE;
	for (int settlePass = 0; settlePass < 16; ++settlePass) {
	    const BObolViewLodState::MeshPayload *settlePayload =
		controller.getViewLodState()->findMesh(mesh);
	    if (settlePayload &&
		settlePayload->activeLevel == settlePayload->requestedLevel &&
		service.pendingTaskCountForDiagnostics() == 0 &&
		service.queuedResultCountForDiagnostics() == 0)
		break;
	    const uint64_t presentedStarted = controller.beginRenderTiming();
	    controller.completeRenderTiming(presentedStarted);
	    controller.noteFramePresented();
	    controller.clearRenderRequest();
	    BObolProgressiveStatus settleStatus;
	    (void)controller.advanceProgressiveWork(
		&settleOptions, &settleStatus);
	    if (service.pendingTaskCountForDiagnostics() > 0 &&
		wait_for_service(service)) {
		printf("FAIL: LoD view controller could not finish initial "
		       "coverage refinement\n");
		service.stop();
		root->unref();
		bobol_mesh_lod_cache_clear_database(dbip);
		db_close(dbip);
		bu_file_delete(dbpath);
		bu_dirclear(cache_dir);
		return 1;
	    }
	}
	const BObolViewLodState::MeshPayload *settledPayload =
	    controller.getViewLodState()->findMesh(mesh);
	if (!settledPayload ||
	    settledPayload->activeLevel != settledPayload->requestedLevel) {
	    printf("FAIL: LoD view controller stranded its coverage-first "
		   "minimum (active=%d requested=%d tasks=%u cuts=%u "
		   "skipped=%u progressive=%d frame=%d pending=%zu "
		   "queued=%zu render=%d reason=%s diagnostic=%s)\n",
		   settledPayload ? settledPayload->activeLevel : -1,
		   settledPayload ? settledPayload->requestedLevel : -1,
		   controller.getLastLodSubmittedTaskCount(),
		   controller.getLastLodUpdatedCutCount(),
		   controller.getLastLodSkippedMeshCount(),
		   controller.hasProgressiveWorkPending() ? 1 : 0,
		   controller.hasPendingLodRefinementFrame() ? 1 : 0,
		   service.pendingTaskCountForDiagnostics(),
		   service.queuedResultCountForDiagnostics(),
		   controller.isRenderRequested() ? 1 : 0,
		   controller.getRenderReason().getString(),
		   controller.getLastLodDiagnostics().getString());
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	controller.clearRenderRequest();
	if (controller.submitLodRequestsIfNeeded() != 0 ||
	    controller.isRenderRequested()) {
	    printf("FAIL: LoD view controller did not skip resident "
		   "unchanged LoD request (visited/submitted/skipped/render="
		   "%u/%u/%u/%d diagnostic=%s)\n",
		   controller.getLastLodVisitedMeshCount(),
		   controller.getLastLodSubmittedTaskCount(),
		   controller.getLastLodSkippedMeshCount(),
		   controller.isRenderRequested() ? 1 : 0,
		   controller.getLastLodDiagnostics().getString());
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (controller.submitLodRequestsIfNeeded() != 0) {
	    printf("FAIL: LoD view controller duplicated unchanged auto-submit\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (service.queuedResultCountForDiagnostics() != 0) {
	    printf("FAIL: LoD view controller queued duplicate resident request result\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
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
	    controller.isRenderRequested()) {
	    printf("FAIL: LoD view controller did not skip resident threshold policy request\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (service.queuedResultCountForDiagnostics() != 0) {
	    printf("FAIL: LoD view controller queued duplicate threshold request result\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	uint64_t previousViewRevision = controller.getLodViewRevision();
	/* Cross a projected-level boundary.  A small camera change that remains
	 * within the resident level must not restart PoP loading. */
	/* Interactive demand deliberately targets 4 px error.  Zoom far enough
	 * to cross a cut under that responsiveness policy. */
	camera->height = 2.0f;
	controller.clearRenderRequest();
	if (controller.submitLodRequestsIfNeeded() != 0 ||
	    controller.getLastLodSubmittedTaskCount() != 0 ||
	    controller.getLastLodUpdatedCutCount() != 1 ||
	    controller.getLodViewRevision() == previousViewRevision ||
	    !controller.isRenderRequested() ||
	    bu_strcmp(controller.getRenderReason().getString(),
		   "lod-cut") != 0) {
	    const BObolViewLodState::MeshPayload *failedPayload =
		controller.getViewLodState()->findMesh(mesh);
	    printf("FAIL: LoD view controller did not retarget the resident "
		   "camera cut (tasks=%d cuts=%d skipped=%d view=%llu/%llu render=%d "
		   "active=%d requested=%d resident=%d max=%d "
		   "reason=%s diagnostic=%s)\n",
		   controller.getLastLodSubmittedTaskCount(),
		   controller.getLastLodUpdatedCutCount(),
		   controller.getLastLodSkippedMeshCount(),
		   (unsigned long long)controller.getLodViewRevision(),
		   (unsigned long long)previousViewRevision,
		   controller.isRenderRequested() ? 1 : 0,
		   failedPayload ? failedPayload->activeLevel : -1,
		   failedPayload ? failedPayload->requestedLevel : -1,
		   failedPayload ? failedPayload->residentLevel : -1,
		   failedPayload && failedPayload->progressiveMesh ?
		       failedPayload->progressiveMesh->maximumLevel() : -1,
		   controller.getRenderReason().getString(),
		   controller.getLastLodDiagnostics().getString());
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (service.pendingTaskCountForDiagnostics() != 0 ||
	    service.queuedResultCountForDiagnostics() != 0) {
	    printf("FAIL: resident camera cut incorrectly scheduled provider work\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	const BObolViewLodState::MeshPayload *activePayload =
	    controller.getViewLodState()->findMesh(mesh);
	if (activePayload && activePayload->progressiveMesh) {
	    const int savedActiveLevel = activePayload->activeLevel;
	    const int savedRequestedLevel = activePayload->requestedLevel;
	    const int residentLevel =
		activePayload->progressiveMesh->residentLevel();
	    const int pendingRequestedLevel = std::min(
		activePayload->progressiveMesh->maximumLevel(),
		residentLevel + 1);
	    if (pendingRequestedLevel > residentLevel) {
		std::vector<BObolLodResidentDemand> demands;
		if (!controller.getViewLodState()->retargetMeshPayload(
			activePayload, savedActiveLevel, pendingRequestedLevel,
			controller.getLodViewRevision(),
			controller.getLodPolicyRevision())) {
		    printf("FAIL: resident PoP cut could not retain a richer pending demand\n");
		    service.stop();
		    root->unref();
		    bobol_mesh_lod_cache_clear_database(dbip);
		    db_close(dbip);
		    bu_file_delete(dbpath);
		    bu_dirclear(cache_dir);
		    return 1;
		}
		controller.getViewLodState()->residentMeshDemands(demands);
		if (activePayload->activeLevel != savedActiveLevel ||
		    activePayload->requestedLevel != pendingRequestedLevel ||
		    demands.size() != 1 ||
		    demands[0].level != savedActiveLevel) {
		    printf("FAIL: resident demand did not follow the stable active PoP cut\n");
		    service.stop();
		    root->unref();
		    bobol_mesh_lod_cache_clear_database(dbip);
		    db_close(dbip);
		    bu_file_delete(dbpath);
		    bu_dirclear(cache_dir);
		    return 1;
		}
		if (!controller.getViewLodState()->retargetMeshPayload(
			activePayload, savedActiveLevel, savedRequestedLevel,
			controller.getLodViewRevision(),
			controller.getLodPolicyRevision())) {
		    printf("FAIL: resident PoP demand fixture could not restore its view target\n");
		    service.stop();
		    root->unref();
		    bobol_mesh_lod_cache_clear_database(dbip);
		    db_close(dbip);
		    bu_file_delete(dbpath);
		    bu_dirclear(cache_dir);
		    return 1;
		}
	    }
	}
	const int residentCameraLevel =
	    activePayload ? activePayload->activeLevel : -1;
	previousViewRevision = controller.getLodViewRevision();
	camera->height = 2.002f;
	controller.clearRenderRequest();
	/* The changed camera still needs a frame.  What must remain idle is the
	 * PoP provider: no new task, result, or payload replacement. */
	const int jitterSubmitted = controller.submitLodRequestsIfNeeded();
	if (jitterSubmitted != 0 ||
	    controller.getLastLodSubmittedTaskCount() != 0 ||
	    controller.getLastLodSkippedMeshCount() != 1 ||
	    controller.getLodViewRevision() == previousViewRevision ||
	    controller.getViewLodState()->findMesh(mesh) != activePayload ||
	    !activePayload ||
	    activePayload->activeLevel != residentCameraLevel ||
	    service.queuedResultCountForDiagnostics() != 0 ||
	    !controller.isRenderRequested()) {
	    printf("FAIL: sub-level camera jitter restarted resident PoP loading "
		   "(return=%d tasks=%d skipped=%d revision=%llu/%llu "
		   "payload=%p/%p level=%d/%d queued=%zu render=%d diagnostic=%s)\n",
		   jitterSubmitted,
		   controller.getLastLodSubmittedTaskCount(),
		   controller.getLastLodSkippedMeshCount(),
		   (unsigned long long)controller.getLodViewRevision(),
		   (unsigned long long)previousViewRevision,
		   (const void *)controller.getViewLodState()->findMesh(mesh),
		   (const void *)activePayload,
		   activePayload ? activePayload->activeLevel : -1,
		   residentCameraLevel,
		   service.queuedResultCountForDiagnostics(),
		   controller.isRenderRequested() ? 1 : 0,
		   controller.getLastLodDiagnostics().getString());
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	int forcedLevel = activePayload ? activePayload->activeLevel : 0;
	uint64_t previousPolicyRevision = controller.getLodPolicyRevision();
	controller.clearRenderRequest();
	controller.setLodForcedLevel(forcedLevel);
	if (!controller.hasLodForcedLevel() ||
	    controller.getLodForcedLevel() != forcedLevel ||
	    controller.getLodPolicyRevision() == previousPolicyRevision ||
	    !controller.isRenderRequested() ||
	    bu_strcmp(controller.getRenderReason().getString(),
		   "lod-policy") != 0) {
	    printf("FAIL: LoD view controller did not record forced-level policy\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	controller.clearRenderRequest();
	if (controller.submitLodRequestsIfNeeded() != 1 ||
	    controller.getLastLodSubmittedTaskCount() != 1 ||
	    controller.isRenderRequested()) {
	    printf("FAIL: LoD view controller did not auto-submit forced-level policy\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (wait_for_service(service)) {
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	controller.clearRenderRequest();
	if (controller.processPendingLodResults() != 1 ||
	    controller.getLastLodAppliedResultCount() != 1 ||
	    !controller.getViewLodState()->findMesh(mesh) ||
	    controller.getViewLodState()->findMesh(mesh)->activeLevel !=
	    forcedLevel ||
	    bu_strcmp(controller.getRenderReason().getString(),
		   "lod-result") != 0) {
	    printf("FAIL: LoD view controller did not apply forced-level result\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
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
	    bu_strcmp(controller.getRenderReason().getString(),
		   "lod-policy") != 0) {
	    printf("FAIL: LoD view controller did not clear forced-level policy\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
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
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	controller.clearRenderRequest();
	if (controller.submitLodRequestsIfNeeded() != 0 ||
	    controller.getLastLodSubmittedTaskCount() != 0 ||
	    service.pendingTaskCountForDiagnostics() != 0 ||
	    service.inFlightCount() != 0) {
	    printf("FAIL: LoD view controller reloaded a resident cut for a "
		   "sub-level viewport change\n");
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}

	/* A zoom is an active demand epoch, not pose-only frame reuse.  A richer
	 * suffix which finishes while the gesture remains bracketed must publish
	 * immediately in a bounded progressive pump; waiting for release leaves a
	 * visibly magnified coarse image throughout continuous wheel/trackpad
	 * input.  Start the owned generation before changing scale so the camera
	 * revision remains the result's authoritative epoch. */
	const uint64_t zoomGeneration = controller.beginLodGeneration();
	controller.beginLodInteraction();
	camera->height = 0.25f;
	controller.clearRenderRequest();
	(void)controller.submitLodRequestsIfNeeded();
	const BObolViewLodState::MeshPayload *zoomBefore =
	    controller.getViewLodState()->findMesh(mesh);
	const int zoomPublishedLevel = zoomBefore ?
	    std::max(zoomBefore->activeLevel + 1, 42) : 42;
	BObolLodRequest zoomRequest =
	    make_request("/lod-submit.bot", "lod-submit.bot");
	zoomRequest.sourceRevision = 303;
	zoomRequest.viewRevision = controller.getLodViewRevision();
	zoomRequest.policyRevision = controller.getLodPolicyRevision();
	zoomRequest.requestedLevel = zoomPublishedLevel;
	BObolLodTask zoomResultTask;
	zoomResultTask.generation = zoomGeneration;
	zoomResultTask.request = zoomRequest;
	zoomResultTask.realize = mesh_payload_requested_refined_task;
	zoomResultTask.debugDelayMilliseconds = 200;
	if (!zoomGeneration || service.submit(zoomResultTask) == 0) {
	    printf("FAIL: active zoom refinement result setup\n");
	    controller.endLodInteraction();
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	if (wait_for_service(service)) {
	    printf("FAIL: active zoom refinement did not finish under bracketed input\n");
	    controller.endLodInteraction();
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
	    db_close(dbip);
	    bu_file_delete(dbpath);
	    bu_dirclear(cache_dir);
	    return 1;
	}
	BObolProgressiveStatus zoomStatus;
	(void)controller.advanceProgressiveWork(NULL, &zoomStatus);
	const BObolViewLodState::MeshPayload *zoomAfter =
	    controller.getViewLodState()->findMesh(mesh);
	const unsigned int callbacksBeforePresentation =
	    frameRequestCount.load(std::memory_order_relaxed);
	const uint64_t zoomPresentationStarted = controller.beginRenderTiming();
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
	controller.completeRenderTiming(zoomPresentationStarted);
	const bool qualityResumeRequested =
	    frameRequestCount.load(std::memory_order_relaxed) >
		callbacksBeforePresentation;
	controller.endLodInteraction();
	if (zoomStatus.lodResultsProcessed == 0 ||
	    zoomStatus.lodResultsApplied == 0 ||
	    !zoomAfter || zoomAfter == zoomBefore ||
	    zoomAfter->activeLevel != zoomPublishedLevel ||
	    zoomAfter->getTriangleCount() != 3 ||
	    !qualityResumeRequested) {
	    printf("FAIL: scale-changing interaction held its richer result "
		   "until quiet (processed=%zu applied=%zu payload=%p/%p "
		   "level=%d triangles=%d resume=%d interaction=%d scale=%d "
		   "diagnostics=%s)\n",
		   zoomStatus.lodResultsProcessed,
		   zoomStatus.lodResultsApplied,
		   (const void *)zoomAfter, (const void *)zoomBefore,
		   zoomAfter ? zoomAfter->activeLevel : -1,
		   zoomAfter ? zoomAfter->getTriangleCount() : -1,
		   qualityResumeRequested ? 1 : 0,
		   controller.isLodInteractionActive() ? 1 : 0,
		   controller.isLodScaleChangingInteraction() ? 1 : 0,
		   controller.getLastLodDiagnostics().getString());
	    service.stop();
	    root->unref();
	    bobol_mesh_lod_cache_clear_database(dbip);
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
    bobol_mesh_lod_cache_clear_database(dbip);
    db_close(dbip);
    bu_file_delete(dbpath);
    bu_dirclear(cache_dir);
    return ret;
}

static int
test_view_controller_lod_source_churn(void)
{
    static const int churnCount = 64;
    static const uint32_t firstRevision = 1000;
    static const uint32_t finalRevision = 9000;
    char cache_dir[MAXPATHLEN] = {0};
    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;

    bu_dir(cache_dir, MAXPATHLEN, BU_DIR_CURR,
	   "bobol_lod_source_churn_cache", NULL);
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

    BObolLodService service;
    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD source-churn service did not start\n");
	root->unref();
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }
    service.setQueueLimits(4, 4, 4);

    bu_setenv("BOBOL_LOD_TASK_DELAY_MS", "40", 1);

    int ret = 0;
    {
	BObolViewController controller(root, camera);
	controller.setLodService(&service);
	controller.setViewportSize(800, 600);
	controller.setLodPolicyRevision(404);

	for (int i = 0; !ret && i < churnCount; i++) {
	    const uint32_t revision = firstRevision +
		static_cast<uint32_t>(i);
	    if (controller.replaceDatabaseSource("lod-submit.bot", dbip,
		    SoBRLDatabaseSource::SHADED, revision) != 1) {
		printf("FAIL: LoD source churn did not publish cycle %d\n", i);
		ret = 1;
		break;
	    }

	    SoBRLDatabaseSource *source = controller.getDatabaseSource(0);
	    SoBRLLodMeshShape *mesh =
		make_lod_mesh("/lod-submit.bot", "lod-submit.bot");
	    if (!source) {
		printf("FAIL: LoD source churn lost published cycle %d\n", i);
		ret = 1;
		break;
	    }
	    source->lodBotThreshold = 1;
	    source->realizedSourceRevision = revision;
	    source->realizedInputsRevision = source->inputsRevision.getValue();
	    source->realizationStatus = SoBRLDatabaseSource::REALIZED;
	    source->stale = FALSE;
	    source->staleReason = SoBRLDatabaseSource::STALE_NONE;
	    mesh->sourceId = revision;
	    mesh->lodPolicy = 6;
	    source->addChild(mesh);

	    BObolLodRequest residentRequest =
		make_request("/lod-submit.bot", "lod-submit.bot");
	    residentRequest.sourceRevision = revision;
	    BObolLodResult residentMesh =
		mesh_payload_result(residentRequest);
	    BObolLodCounts proxyCounts;
	    proxyCounts.faceCount = 1;
	    proxyCounts.pointCount = 3;
	    BObolLodResult residentProxy = bobol_lod_aabb_result(
		residentRequest, residentRequest.bounds, &proxyCounts);
	    if (!controller.getViewLodState()->applyMeshResult(mesh,
		    residentMesh) ||
		!controller.getViewLodState()->applyProxyResult(mesh,
		    residentProxy) ||
		controller.getActiveLodMeshPayloadCount() != 1 ||
		controller.getActiveLodProxyPayloadCount(
		    BOBOL_LOD_PROXY_AABB) != 1) {
		printf("FAIL: LoD source churn did not seed cycle %d view state\n",
		       i);
		ret = 1;
		break;
	    }

	    controller.clearRenderRequest();
	    if (controller.submitLodRequestsIfNeeded(TRUE, 1) != 1 ||
		controller.getLastLodSubmittedTaskCount() != 1 ||
		wait_for_service_delayed(service)) {
		printf("FAIL: LoD source churn did not start cycle %d work\n", i);
		ret = 1;
		break;
	    }
	    const uint64_t generation = service.currentGeneration();
	    if (generation == 0 ||
		controller.setDatabaseSourceState(
		    "lod-submit.bot", FALSE, revision,
		    source->inputsRevision.getValue(),
		    source->visible.getValue(), !source->selected.getValue(),
		    source->highlighted.getValue(), source->lineStyle.getValue(),
		    source->lineWidth.getValue(), source->transparency.getValue(),
		    source->colorOverride.getValue(), source->color.getValue(),
		    source->materialColorValid.getValue(),
		    source->materialColor.getValue(),
		    source->materialRevision.getValue()) != 1 ||
		service.isGenerationCancelled(generation) ||
		controller.getActiveLodMeshPayloadCount() != 1 ||
		controller.getActiveLodProxyPayloadCount(
		    BOBOL_LOD_PROXY_AABB) != 1 ||
		service.inFlightCount() > 4 ||
		service.pendingTaskCountForDiagnostics() > 4 ||
		service.queuedResultCountForDiagnostics() > 4 ||
		service.queuedCacheWriteCountForDiagnostics() > 4) {
		printf("FAIL: LoD source churn exceeded bounded cycle %d work\n", i);
		ret = 1;
		break;
	    }

	    int changed = 0;
	    switch (i % 4) {
		case 0:
		    changed = controller.markDatabaseSourceStale(
			"lod-submit.bot", SoBRLDatabaseSource::STALE_SOURCE);
		    break;
		case 1:
		    changed = controller.setDatabaseSourceState(
			"lod-submit.bot", TRUE, revision + 1,
			source->inputsRevision.getValue(),
			source->visible.getValue(), source->selected.getValue(),
			source->highlighted.getValue(), source->lineStyle.getValue(),
			source->lineWidth.getValue(), source->transparency.getValue(),
			source->colorOverride.getValue(), source->color.getValue(),
			source->materialColorValid.getValue(),
			source->materialColor.getValue(),
			source->materialRevision.getValue());
		    break;
		case 2:
		    changed = controller.setDatabaseSourceBoundsState(
			"lod-submit.bot", TRUE,
			SbVec3f(0.0f, 0.0f, 0.0f),
			SbVec3f(static_cast<float>(i + 2), 2.0f, 2.0f));
		    break;
		default:
		    changed = controller.removeDatabaseSource("lod-submit.bot");
		    break;
	    }

	    if (changed != 1 || !service.isGenerationCancelled(generation) ||
		controller.getActiveLodMeshPayloadCount() != 0 ||
		controller.getActiveLodProxyPayloadCount(
		    BOBOL_LOD_PROXY_NONE) != 0 ||
		controller.getActiveLodCadPayloadCount() != 0 ||
		wait_for_service_idle(service) ||
		controller.processPendingLodResults() != 0) {
		printf("FAIL: LoD source churn did not cancel and clear cycle %d\n",
		       i);
		ret = 1;
		break;
	    }

	    if (controller.getDatabaseSourceCount() > 0 &&
		controller.removeDatabaseSource("lod-submit.bot") != 1) {
		printf("FAIL: LoD source churn did not erase cycle %d source\n", i);
		ret = 1;
		break;
	    }
	    if (controller.getDatabaseSourceCount() != 0) {
		printf("FAIL: LoD source churn retained cycle %d source\n", i);
		ret = 1;
	    }
	}

	bu_setenv("BOBOL_LOD_TASK_DELAY_MS", "0", 1);

	if (!ret) {
	    if (controller.replaceDatabaseSource("lod-submit.bot", dbip,
		    SoBRLDatabaseSource::SHADED, finalRevision) != 1) {
		printf("FAIL: LoD source churn did not publish final source\n");
		ret = 1;
	    }
	}

	SoBRLDatabaseSource *finalSource = ret ? NULL :
	    controller.getDatabaseSource(0);
	SoBRLLodMeshShape *finalMesh = NULL;
	if (!ret && finalSource) {
	    finalSource->lodBotThreshold = 1;
	    finalSource->realizedSourceRevision = finalRevision;
	    finalSource->realizedInputsRevision =
		finalSource->inputsRevision.getValue();
	    finalSource->realizationStatus = SoBRLDatabaseSource::REALIZED;
	    finalSource->stale = FALSE;
	    finalSource->staleReason = SoBRLDatabaseSource::STALE_NONE;
	    finalMesh = make_lod_mesh("/lod-submit.bot", "lod-submit.bot");
	    finalMesh->sourceId = finalRevision;
	    finalMesh->lodPolicy = 6;
	    finalSource->addChild(finalMesh);

	    if (controller.submitLodRequestsIfNeeded(TRUE, 1) != 1 ||
		wait_for_service(service) ||
		controller.processPendingLodResults() != 1 ||
		controller.getLastLodMatchedResultCount() != 1 ||
		controller.getLastLodAppliedResultCount() != 1 ||
		controller.getLastLodRejectedResultCount() != 0 ||
		controller.getLastLodUnmatchedResultCount() != 0 ||
		controller.getActiveLodMeshPayloadCount() != 1 ||
		controller.getActiveLodProxyPayloadCount(
		    BOBOL_LOD_PROXY_AABB) != 0 ||
		controller.getActiveLodProxyPayloadCount(
		    BOBOL_LOD_PROXY_OBB) != 0 ||
		controller.getActiveLodCadPayloadCount() != 0 ||
		!controller.getViewLodState()->findMesh(finalMesh) ||
		controller.getDatabaseSourceCount() != 1 ||
		bu_strcmp(finalSource->path.getValue().getString(),
		    "lod-submit.bot") != 0 ||
		finalSource->sourceRevision.getValue() != finalRevision) {
		printf("FAIL: LoD source churn final generation was not authoritative\n");
		ret = 1;
	    }
	}
	if (!ret && (!finalSource || !finalMesh || wait_for_service_idle(service)))
	    ret = 1;
    }

    bu_setenv("BOBOL_LOD_TASK_DELAY_MS", "0", 1);
    service.stop();
    root->unref();
    bobol_mesh_lod_cache_clear_database(dbip);
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
	BObolViewController viewA(shared, NULL);
	BObolViewController viewB(shared, NULL);

	BObolLodRequest request = make_request("/shared/lod.bot",
				    "lod.bot");
	BObolLodResult resultA =
	    mesh_payload_variant_result(request, 1, 2);
	BObolLodResult resultB =
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

    BObolLodService service;
    if (!service.start(2, TRUE)) {
	printf("FAIL: progressive LoD service did not start\n");
	root->unref();
	return 1;
    }

    int ret = 0;
    {
	BObolViewController controller(root, NULL);
	controller.setLodService(&service);
	controller.setLodPolicyRevision(17);
	controller.clearRenderRequest();
	std::atomic<unsigned int> frameRequestCount(0);
	controller.setFrameRequestCallback(
	    count_frame_request, &frameRequestCount);

	BObolLodRequest request = make_request("/progress/lod.bot", "lod.bot");
	request.viewRevision = controller.getLodViewRevision();
	request.policyRevision = controller.getLodPolicyRevision();
	request.requestedLevel = 2;

	uint64_t generation = controller.beginLodGeneration();
	BObolLodTask coarseTask;
	coarseTask.generation = generation;
	coarseTask.request = request;
	coarseTask.realize = mesh_payload_coarse_task;
	BObolLodTask refinedTask;
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
	    const unsigned int callbacksBefore =
		frameRequestCount.load(std::memory_order_relaxed);
	    controller.clearRenderRequest();
	    if (controller.processPendingLodResults(1) != 1 ||
		controller.getLastLodAppliedResultCount() != 1 ||
		!controller.isRenderRequested() ||
		!controller.hasPendingLodRefinementFrame() ||
		frameRequestCount.load(std::memory_order_relaxed) <=
		    callbacksBefore ||
		bu_strcmp(controller.getRenderReason().getString(),
		       "lod-result") != 0) {
		printf("FAIL: progressive LoD controller did not apply and schedule the coarse result frame\n");
		ret = 1;
	    }
	}

	if (!ret) {
	    const BObolViewLodState::MeshPayload *payload =
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

	if (!ret) {
	    const uint64_t started = controller.beginRenderTiming();
	    std::this_thread::sleep_for(std::chrono::milliseconds(1));
	    controller.completeRenderTiming(started);
	    if (controller.hasPendingLodRefinementFrame() ||
		!controller.hasProgressiveWorkPending() ||
		!controller.isRenderRequested()) {
		printf("FAIL: completed coarse presentation did not release progressive LoD refinement\n");
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
		bu_strcmp(controller.getRenderReason().getString(),
		       "lod-result") != 0) {
		printf("FAIL: progressive LoD controller did not apply refined result\n");
		ret = 1;
	    }
	}

	if (!ret) {
	    const BObolViewLodState::MeshPayload *payload =
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
	controller.clearFrameRequestCallback(&frameRequestCount);
    }

    service.stop();
    root->unref();
    return ret;
}

static int
test_view_controller_large_result_transfer(void)
{
    SoSeparator *root = new SoSeparator;
    root->ref();
    SoBRLLodMeshShape *mesh = make_lod_mesh("/large/lod.bot", "lod.bot");
    root->addChild(mesh);

    BObolLodService service;
    if (!service.start(1, FALSE)) {
	root->unref();
	return 1;
    }

    int ret = 0;
    {
	BObolViewController controller(root, NULL);
	controller.setLodService(&service);
	controller.setLodPolicyRevision(23);
	BObolLodRequest request = make_request("/large/lod.bot", "lod.bot");
	request.viewRevision = controller.getLodViewRevision();
	request.policyRevision = controller.getLodPolicyRevision();

	large_payload_task_data data;
	data.pointCount = 1500000;
	BObolLodTask task;
	task.generation = controller.beginLodGeneration();
	task.request = request;
	task.realize = large_mesh_payload_task;
	task.realizeData = &data;
	if (service.submit(task) == 0 || wait_for_service(service)) {
	    ret = 1;
	} else {
	    const uintptr_t workerPoints = data.pointStorage.load();
	    const uintptr_t workerIndices = data.indexStorage.load();
	    const int64_t started = bu_gettime();
	    const size_t processed = controller.processPendingLodResults(1, 4000);
	    const int64_t elapsed = bu_gettime() - started;
	    const BObolViewLodState::MeshPayload *payload =
		controller.getViewLodState()->findMesh(mesh);
	    if (processed != 1 ||
		controller.getLastLodAppliedResultCount() != 1 || !payload ||
		payload->mesh.points.size() != data.pointCount ||
		workerPoints == 0 || workerIndices == 0 ||
		reinterpret_cast<uintptr_t>(payload->mesh.points.data()) !=
		    workerPoints ||
		reinterpret_cast<uintptr_t>(payload->mesh.coordIndex.data()) !=
		    workerIndices || elapsed < 0 || elapsed > 4000) {
		printf("FAIL: large LoD result transfer copied storage or exceeded bounded apply time (elapsed=%lld us)\n",
		    (long long)elapsed);
		ret = 1;
	    }
	}
	controller.setLodService(NULL);
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
	BObolViewController controller(root, NULL);
	BObolLodRequest request = make_request("/pick/lod.bot", "lod.bot");
	request.viewRevision = controller.getLodViewRevision();
	request.policyRevision = controller.getLodPolicyRevision();
	BObolLodResult result = mesh_payload_variant_result(request, 1, 1);

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
		bu_strcmp(pickDetail->getPath().getString(), "/pick/lod.bot") != 0 ||
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

static std::shared_ptr<const Obol::PartGeometry>
compact_duplicate_proxy_geometry(void)
{
    std::shared_ptr<Obol::PartGeometry> geometry(new Obol::PartGeometry);
    Obol::WireRep wire;
    wire.segmentPoints.push_back(SbVec3f(0.0f, 0.0f, 0.0f));
    wire.segmentPoints.push_back(SbVec3f(1.0f, 0.0f, 0.0f));
    wire.segmentIds.push_back(1);
    wire.bounds = SbBox3f(SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 0.0f, 0.0f));
    geometry->wire = std::move(wire);
    return geometry;
}

static std::shared_ptr<const Obol::PartGeometry>
compact_projected_proxy_geometry(void)
{
    std::shared_ptr<Obol::PartGeometry> geometry(new Obol::PartGeometry);
    Obol::WireRep wire;
    wire.segmentPoints.push_back(SbVec3f(0.0f, 0.0f, 0.0f));
    wire.segmentPoints.push_back(SbVec3f(1.0f, 1.0f, 1.0f));
    wire.segmentIds.push_back(1);
    wire.bounds = SbBox3f(SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 1.0f, 1.0f));
    geometry->wire = std::move(wire);
    geometry->subpixelProxyEligible = true;
    geometry->structuralProxy = true;
    return geometry;
}

static std::shared_ptr<const Obol::PartGeometry>
compact_native_triangle_geometry(void)
{
    std::shared_ptr<Obol::PartGeometry> geometry(new Obol::PartGeometry);
    Obol::TriMesh mesh;
    mesh.positions.push_back(SbVec3f(0.0f, 0.0f, 0.0f));
    mesh.positions.push_back(SbVec3f(1.0f, 0.0f, 0.0f));
    mesh.positions.push_back(SbVec3f(0.0f, 1.0f, 0.0f));
    mesh.indices.push_back(0);
    mesh.indices.push_back(1);
    mesh.indices.push_back(2);
    mesh.bounds = SbBox3f(SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 1.0f, 0.0f));
    geometry->shaded = std::move(mesh);
    return geometry;
}

static int
test_compact_many_leaf_scene_admission(void)
{
    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;
    if (make_submit_test_db(dbpath, sizeof(dbpath), &dbip))
	return 1;

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->ref();
    source->setDatabase(dbip);
    source->path = "many-leaf-root.c";
    source->instanceKey = "many-leaf-source";
    source->sourceRevision = 901;
    source->lodBotThreshold = 1;
    source->representationMode =
	SoBRLDatabaseSource::REPRESENTATION_SHADED;

    /*
     * Worker roles describe publication ownership, not the immutable
     * geometry route.  A warm external publication must not make a wire BoT
     * fall back to the legacy plotted-vlist walker and overwrite its PoP
     * request contract.
     */
    SoBRLDatabaseSource *wireContract = new SoBRLDatabaseSource;
    wireContract->ref();
    wireContract->drawMode = SoBRLDatabaseSource::WIREFRAME;
    wireContract->representationMode =
	SoBRLDatabaseSource::REPRESENTATION_WIRE;
    wireContract->realizationRoleFlags =
	SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL;
    wireContract->realizationViewDependent = TRUE;
    wireContract->realizationMeshLodEnabled = TRUE;
    wireContract->lodBotThreshold = 1;
    const bool wireUsesMesh = wireContract->usesMeshRealization() ? true : false;
    wireContract->unref();
    if (!wireUsesMesh) {
	printf("FAIL: external view-managed wire source lost its mesh "
	       "realization contract\n");
	source->unref();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    static const size_t leafCount = 128;
    std::vector<BObolCompactOccurrence> occurrences;
    occurrences.reserve(leafCount);
    for (size_t i = 0; i < leafCount; ++i) {
	BObolCompactOccurrence occurrence;
	occurrence.geometry = compact_projected_proxy_geometry();
	occurrence.summary.valid = TRUE;
	occurrence.summary.shapeKind =
	    BObolRealizedShapeSummary::SHAPE_VLIST;
	SbString path;
	path.sprintf("many-leaf-root.c/leaf-%zu.bot", i);
	occurrence.summary.path = path;
	occurrence.summary.sourceName = "lod-submit.bot";
	occurrence.summary.sourceType = "bot";
	occurrence.summary.sourceId = static_cast<uint32_t>(1000 + i);
	occurrence.summary.visible = TRUE;
	occurrence.summary.selectable = TRUE;
	occurrence.summary.selected =
	    i == leafCount - 1 ? TRUE : FALSE;
	occurrence.lodBacked = TRUE;
	occurrence.sourceMeshRequestValid = TRUE;
	occurrence.sourceMeshRequest.path = path;
	occurrence.sourceMeshRequest.sourceName = "lod-submit.bot";
	/* Every occurrence aliases one source asset, while its source count is
	 * deliberately large enough to exercise the provisional first-cut
	 * reservation rather than the tiny tetrahedron fixture. */
	occurrence.sourceMeshRequest.meshAssetPath = "lod-submit.bot";
	occurrence.sourceMeshRequest.meshAssetName = "lod-submit.bot";
	occurrence.sourceMeshRequest.faceCount = 10000;
	occurrence.sourceMeshRequest.pointCount = 5000;
	occurrence.sourceMeshRequest.bounds = SbBox3f(
	    SbVec3f(0.0f, 0.0f, 0.0f),
	    SbVec3f(1.0f, 1.0f, 1.0f));
	SbMatrix placement;
	placement.setTranslate(SbVec3f(
	    static_cast<float>(i % 16) * 1.5f,
	    static_cast<float>(i / 16) * 1.5f, 0.0f));
	occurrence.localTransform = placement;
	occurrences.push_back(occurrence);
    }

    int ret = 0;
    if (source->setCompactOccurrenceRegistry(occurrences) !=
	    static_cast<int>(leafCount)) {
	printf("FAIL: many-leaf compact registry setup\n");
	ret = 1;
    }
    if (!ret &&
	source->getDisplayMeshLodRequestCount() != leafCount) {
	printf("FAIL: many-leaf compact mesh-request contract count "
	       "(requests=%zu expected=%zu)\n",
	       source->getDisplayMeshLodRequestCount(), leafCount);
	ret = 1;
    }
    if (!ret) {
	/*
	 * Source-mesh identity is independent of view and presentation policy.
	 * A camera-policy invalidation may start a new detached realization, but
	 * it must not make already-published PoP requests disappear meanwhile.
	 * Wireframe has no shaded source fallback and exposed this as a permanent
	 * zero-progress submission loop on the first warm draw.
	 */
	source->markStale(SoBRLDatabaseSource::STALE_VIEW);
	if (!source->hasDisplayMeshLodRequests() ||
	    source->getDisplayMeshLodRequestCount() != leafCount) {
	    printf("FAIL: view staleness invalidated immutable compact "
		   "mesh-request contracts\n");
	    ret = 1;
	}
    }
    if (!ret) {
	std::vector<SbString> selectedPaths(1,
	    occurrences.back().summary.path);
	if (source->syncCompactInstanceSelectedPaths(selectedPaths) <= 0) {
	    printf("FAIL: many-leaf compact priority selection setup\n");
	    ret = 1;
	}
    }
    if (!ret) {
	const uint64_t beforeVisibility =
	    source->getDisplayMeshLodRevision();
	std::vector<SbString> paths(1, occurrences[7].summary.path);
	std::vector<SbBool> states(1, FALSE);
	std::vector<size_t> changedEntries;
	if (source->setCompactInstanceVisibilityOverrides(paths, states) <= 0 ||
	    !source->getDisplayMeshLodChangedEntries(
		beforeVisibility, changedEntries) ||
	    changedEntries.size() != 1 || changedEntries[0] != 7) {
	    printf("FAIL: compact visibility did not publish a selective "
		   "LoD demand delta\n");
	    ret = 1;
	} else {
	    const uint64_t hiddenRevision =
		source->getDisplayMeshLodRevision();
	    changedEntries.clear();
	    if (source->clearCompactInstanceVisibilityFrontier() <= 0 ||
		!source->getDisplayMeshLodChangedEntries(
		    hiddenRevision, changedEntries) ||
		changedEntries.size() != 1 || changedEntries[0] != 7) {
		printf("FAIL: compact visibility restore did not publish a "
		       "selective LoD demand delta\n");
		ret = 1;
	    }
	}
    }
    BObolLodService service;
    if (!ret && !service.start(1, TRUE)) {
	printf("FAIL: many-leaf LoD service setup\n");
	ret = 1;
    }

    struct bv_view_info view = BV_VIEW_INFO_INIT;
    view.width = 640;
    view.height = 480;
    view.size = 64.0;
    SoOrthographicCamera *camera = new SoOrthographicCamera;
    camera->ref();
    camera->position = SbVec3f(12.0f, 6.0f, 5.0f);
    camera->height = 32.0f;
    const SbViewVolume volume =
	camera->getViewVolume(640.0f / 480.0f);

    std::vector<size_t> pinnedPlan;
    BObolViewLodState sharedViewState;
    if (!ret) {
	SoBRLMeshLodSubmitAction firstWindow;
	firstWindow.setService(&service);
	firstWindow.setDatabase(dbip, "db://many-leaf-test", 2026);
	firstWindow.setViewInfo(&view);
	firstWindow.setViewVolume(&volume, 1.0f);
	firstWindow.setGeneration(service.beginGeneration());
	firstWindow.setRevisions(1, 1);
	const size_t provisionalCost = 832;
	firstWindow.setRefinementCostBudget(16 * provisionalCost);
	firstWindow.setCompactEntryRange(0, 16);
	firstWindow.apply(source);
	firstWindow.getCompactEntryPlan(pinnedPlan);
	if (firstWindow.getVisitedMeshCount() != 16 ||
	    firstWindow.getVisibleMeshCount() != 16 ||
	    firstWindow.getCoveredVisibleMeshCount() != 0 ||
	    firstWindow.getSubmittedTaskCount() != 16 ||
	    firstWindow.getRefinementCostBudgetUsed() !=
		16 * provisionalCost ||
	    firstWindow.getRefinementBudgetBlockedCount() != 0 ||
	    !firstWindow.hasDeferredCompactEntries() ||
	    firstWindow.getCompactEntryNext() != 16 ||
	    pinnedPlan.size() != leafCount ||
	    pinnedPlan.front() != leafCount - 1) {
	    printf("FAIL: many-leaf aggregate admission or priority plan "
		   "(visited=%u visible=%zu covered=%zu tasks=%u "
		   "used=%zu blocked=%u deferred=%d "
		   "next=%zu plan=%zu first=%zu)\n",
		   firstWindow.getVisitedMeshCount(),
		   firstWindow.getVisibleMeshCount(),
		   firstWindow.getCoveredVisibleMeshCount(),
		   firstWindow.getSubmittedTaskCount(),
		   firstWindow.getRefinementCostBudgetUsed(),
		   firstWindow.getRefinementBudgetBlockedCount(),
		   firstWindow.hasDeferredCompactEntries() ? 1 : 0,
		   firstWindow.getCompactEntryNext(), pinnedPlan.size(),
		   pinnedPlan.empty() ? SIZE_MAX : pinnedPlan.front());
	    ret = 1;
	}
    }

    if (!ret) {
	point_t points[4];
	VSET(points[0], 0.0, 0.0, 0.0);
	VSET(points[1], 1.0, 0.0, 0.0);
	VSET(points[2], 0.0, 1.0, 0.0);
	VSET(points[3], 0.0, 0.0, 1.0);
	int faces[12] = {
	    0, 1, 2, 0, 3, 1, 1, 3, 2, 2, 3, 0
	};
	struct BObolMeshLodData data = {};
	data.faces = faces;
	data.face_count = 4;
	data.points = points;
	data.point_count = 4;
	data.points_orig = points;
	data.point_orig_count = 4;
	VSET(data.bmin, 0.0, 0.0, 0.0);
	VSET(data.bmax, 1.0, 1.0, 1.0);
	struct BObolMeshLodHierarchyInfo hierarchy =
	    BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
	hierarchy.min_level = 0;
	hierarchy.max_level = 4;
	hierarchy.resident_level = 4;
	for (int level = 0; level <= 4; ++level) {
	    hierarchy.face_count[level] = 4;
	    hierarchy.point_count[level] = 4;
	}
	VSET(hierarchy.quantization_min, 0.0, 0.0, 0.0);
	VSET(hierarchy.quantization_max, 1.0, 1.0, 1.0);
	BObolLodProgressiveMeshPtr progressive(
	    new BObolLodProgressiveMesh);
	BObolCompactInstanceHandle handle;
	BObolCompactInstanceSummary summary;
	if (!progressive->update(data, hierarchy, 4, FALSE) ||
	    !source->getCompactInstanceHandle(
		static_cast<int>(pinnedPlan[0]), handle) ||
	    !source->getCompactInstanceSummary(handle, summary)) {
	    printf("FAIL: many-leaf shared asset fixture setup\n");
	    ret = 1;
	} else {
	    BObolLodRequest residentRequest;
	    residentRequest.databaseId = "db://many-leaf-test";
	    residentRequest.databaseRevision = 2026;
	    residentRequest.sourceRevision = source->sourceRevision.getValue();
	    residentRequest.objectPath = summary.meshAssetPath;
	    residentRequest.objectName = summary.meshAssetName;
	    residentRequest.occurrenceKey = summary.sourceInstanceKey;
	    residentRequest.viewRevision = 1;
	    residentRequest.policyRevision = 1;
	    residentRequest.drawMode = BOBOL_LOD_DRAW_SHADED;
	    residentRequest.providerId = "bobol_mesh_lod";
	    residentRequest.providerVersion =
		BOBOL_MESH_LOD_PROVIDER_VERSION;
	    residentRequest.qualityTier =
		BOBOL_LOD_QUALITY_FAST_DISPLAY;
	    residentRequest.requestedLevel = 4;
	    residentRequest.bounds = SbBox3f(
		SbVec3f(0.0f, 0.0f, 0.0f),
		SbVec3f(1.0f, 1.0f, 1.0f));
	    BObolLodResult residentResult;
	    residentResult.request = residentRequest;
	    residentResult.cacheKey =
		bobol_lod_cache_key(residentRequest);
	    residentResult.geometry.kind =
		BOBOL_LOD_GEOMETRY_MESH_LOD_CACHE;
	    residentResult.geometry.providerId =
		residentRequest.providerId;
	    residentResult.geometry.providerVersion =
		residentRequest.providerVersion;
	    residentResult.geometry.cacheKey =
		bobol_lod_asset_cache_key(residentRequest);
	    residentResult.geometry.activeLevel = 4;
	    residentResult.progressiveMesh = progressive;
	    residentResult.residentLevel = 4;
	    residentResult.resultKind = BOBOL_LOD_RESULT_MESH;
	    residentResult.qualityTier =
		BOBOL_LOD_QUALITY_FAST_DISPLAY;
	    residentResult.providerStatus = BOBOL_LOD_PROVIDER_READY;
	    residentResult.bounds = progressive->bounds();
	    residentResult.counts.faceCount = 4;
	    residentResult.counts.pointCount = 4;
	    residentResult.counts.originalPointCount = 4;
	    residentResult.terminal = TRUE;
	    if (!sharedViewState.applySourceResult(
		    source, residentResult)) {
		printf("FAIL: many-leaf shared asset binding setup\n");
		ret = 1;
	    }
	}
    }

    if (!ret) {
	SoBRLMeshLodSubmitAction sharedReuse;
	sharedReuse.setService(&service);
	sharedReuse.setDatabase(dbip, "db://many-leaf-test", 2026);
	sharedReuse.setViewInfo(&view);
	sharedReuse.setViewVolume(&volume, 1.0f);
	sharedReuse.setGeneration(service.beginGeneration());
	sharedReuse.setRevisions(1, 1);
	sharedReuse.setViewLodState(&sharedViewState);
	const size_t residentOccurrenceCost = 69;
	sharedReuse.setRefinementCostBudget(residentOccurrenceCost);
	sharedReuse.setCompactEntryPlan(pinnedPlan);
	sharedReuse.setCompactEntryRange(2, 1);
	sharedReuse.apply(source);
	if (sharedReuse.getVisitedMeshCount() != 1 ||
	    sharedReuse.getSubmittedTaskCount() != 0 ||
	    sharedReuse.getUpdatedCutCount() != 1 ||
	    sharedReuse.getRefinementCostBudgetUsed() !=
		residentOccurrenceCost ||
	    sharedViewState.cadMeshPayloadCount() != 2) {
	    printf("FAIL: many-leaf shared PoP asset was not bound directly "
		   "(visited=%u tasks=%u cuts=%u used=%zu payloads=%zu)\n",
		   sharedReuse.getVisitedMeshCount(),
		   sharedReuse.getSubmittedTaskCount(),
		   sharedReuse.getUpdatedCutCount(),
		   sharedReuse.getRefinementCostBudgetUsed(),
		   sharedViewState.cadMeshPayloadCount());
	    ret = 1;
	}
    }

    if (!ret) {
	SoBRLMeshLodSubmitAction secondWindow;
	secondWindow.setService(&service);
	secondWindow.setDatabase(dbip, "db://many-leaf-test", 2026);
	secondWindow.setViewInfo(&view);
	secondWindow.setViewVolume(&volume, 1.0f);
	secondWindow.setGeneration(service.beginGeneration());
	secondWindow.setRevisions(1, 1);
	secondWindow.setRefinementCostBudget(0);
	secondWindow.setCompactEntryPlan(pinnedPlan);
	secondWindow.setCompactEntryRange(16, 16);
	secondWindow.apply(source);
	std::vector<size_t> continuedPlan;
	secondWindow.getCompactEntryPlan(continuedPlan);
	if (secondWindow.getVisitedMeshCount() != 16 ||
	    secondWindow.getSubmittedTaskCount() != 0 ||
	    secondWindow.getUpdatedCutCount() != 0 ||
	    secondWindow.getRefinementBudgetBlockedCount() != 16 ||
	    secondWindow.getCompactEntryNext() != 32 ||
	    continuedPlan != pinnedPlan ||
	    sharedViewState.cadMeshPayloadCount() != 2) {
	    printf("FAIL: many-leaf pinned plan continuation "
		   "(visited=%u tasks=%u cuts=%u blocked=%u next=%zu "
		   "plan=%zu payloads=%zu)\n",
		   secondWindow.getVisitedMeshCount(),
		   secondWindow.getSubmittedTaskCount(),
		   secondWindow.getUpdatedCutCount(),
		   secondWindow.getRefinementBudgetBlockedCount(),
		   secondWindow.getCompactEntryNext(),
		   continuedPlan.size(),
		   sharedViewState.cadMeshPayloadCount());
	    ret = 1;
	}
    }

    if (!ret) {
	SoBRLMeshLodSubmitAction retainedAdmission;
	retainedAdmission.setService(&service);
	retainedAdmission.setDatabase(
	    dbip, "db://many-leaf-test", 2026);
	retainedAdmission.setViewInfo(&view);
	retainedAdmission.setViewVolume(&volume, 1.0f);
	retainedAdmission.setGeneration(service.beginGeneration());
	retainedAdmission.setRevisions(2, 2);
	retainedAdmission.setViewLodState(&sharedViewState);
	retainedAdmission.setRefinementCostBudget(0);
	const size_t residentOccurrenceCost = 69;
	retainedAdmission.setRetainedSceneCostBudget(
	    residentOccurrenceCost);
	retainedAdmission.setSubmissionTaskLimit(0);
	retainedAdmission.setCompactEntryPlan(pinnedPlan);
	/* Large-scene retained recovery is intentionally windowed.  The first
	 * slice accounts for its own occurrence without scanning/mutating every
	 * payload on the UI thread. */
	retainedAdmission.setCompactEntryRange(0, 1);
	retainedAdmission.apply(source);
	if (retainedAdmission.getRetainedSceneCostBudgetUsed() !=
		residentOccurrenceCost ||
	    retainedAdmission.getUpdatedCutCount() != 0 ||
	    sharedViewState.activeFaceCount() != 8 ||
	    sharedViewState.cadMeshPayloadCount() != 2) {
	    printf("FAIL: bounded retained recovery changed occurrences "
		   "outside its first window "
		   "(used=%zu cuts=%u faces=%zu payloads=%zu)\n",
		   retainedAdmission.getRetainedSceneCostBudgetUsed(),
		   retainedAdmission.getUpdatedCutCount(),
		   sharedViewState.activeFaceCount(),
		   sharedViewState.cadMeshPayloadCount());
	    ret = 1;
	}
    }

    if (!ret) {
	SoBRLMeshLodSubmitAction retainedAdmissionTail;
	retainedAdmissionTail.setService(&service);
	retainedAdmissionTail.setDatabase(
	    dbip, "db://many-leaf-test", 2026);
	retainedAdmissionTail.setViewInfo(&view);
	retainedAdmissionTail.setViewVolume(&volume, 1.0f);
	retainedAdmissionTail.setGeneration(service.beginGeneration());
	retainedAdmissionTail.setRevisions(2, 2);
	retainedAdmissionTail.setViewLodState(&sharedViewState);
	retainedAdmissionTail.setRefinementCostBudget(0);
	retainedAdmissionTail.setRetainedSceneCostBudget(0);
	retainedAdmissionTail.setPreserveMeshCoverage(TRUE);
	retainedAdmissionTail.setSubmissionTaskLimit(0);
	retainedAdmissionTail.setCompactEntryPlan(pinnedPlan);
	retainedAdmissionTail.setCompactEntryRange(1, SIZE_MAX);
	retainedAdmissionTail.apply(source);
	if (retainedAdmissionTail.getRetainedSceneCostBudgetUsed() !=
		(leafCount - 1) * 69 ||
	    retainedAdmissionTail.getUpdatedCutCount() != leafCount - 2 ||
	    retainedAdmissionTail.getRefinementBudgetBlockedCount() == 0 ||
	    sharedViewState.activeFaceCount() != leafCount * 4 ||
	    sharedViewState.cadMeshPayloadCount() != leafCount ||
	    sharedViewState.cadMeshPayloadCountForSource(source) !=
		leafCount) {
	    printf("FAIL: bounded retained recovery tail did not preserve "
		   "the minimum visible mesh floor "
		   "(used=%zu cuts=%u faces=%zu payloads=%zu "
		   "source_payloads=%zu)\n",
		   retainedAdmissionTail.getRetainedSceneCostBudgetUsed(),
		   retainedAdmissionTail.getUpdatedCutCount(),
		   sharedViewState.activeFaceCount(),
		   sharedViewState.cadMeshPayloadCount(),
		   sharedViewState.cadMeshPayloadCountForSource(source));
	    ret = 1;
	}
    }

    /* Re-entering a view with a resident shared asset keeps every minimum
     * mesh available.  The retained budget controls additional submitted
     * faces, but a richer PoP snap with the same face population is free and
     * should advance rather than preserving avoidable geometric error. */
    if (!ret) {
	SoBRLMeshLodSubmitAction coverageReadmission;
	coverageReadmission.setService(&service);
	coverageReadmission.setDatabase(
	    dbip, "db://many-leaf-test", 2026);
	coverageReadmission.setViewInfo(&view);
	coverageReadmission.setViewVolume(&volume, 1.0f);
	coverageReadmission.setGeneration(service.beginGeneration());
	coverageReadmission.setRevisions(3, 3);
	coverageReadmission.setViewLodState(&sharedViewState);
	coverageReadmission.setRefinementCostBudget(0);
	coverageReadmission.setRetainedSceneCostBudget(2 * 69);
	coverageReadmission.setPreserveMeshCoverage(TRUE);
	coverageReadmission.setSubmissionTaskLimit(0);
	coverageReadmission.setCompactEntryPlan(pinnedPlan);
	coverageReadmission.apply(source);
	if (coverageReadmission.getRetainedSceneCostBudgetUsed() !=
		leafCount * 69 ||
	    coverageReadmission.getUpdatedCutCount() != leafCount - 2 ||
	    sharedViewState.activeFaceCount() != leafCount * 4 ||
	    sharedViewState.cadMeshPayloadCount() != leafCount) {
	    size_t occurrenceBindings = 0;
	    size_t shadedBindings = 0;
	    size_t minimumBindings = 0;
	    size_t terminalBindings = 0;
	    for (size_t occurrenceIndex = 0;
		 occurrenceIndex < occurrences.size(); occurrenceIndex++) {
		BObolCompactInstanceHandle occurrenceHandle;
		BObolCompactInstanceSummary occurrenceSummary;
		if (!source->getCompactInstanceHandle(
			static_cast<int>(occurrenceIndex), occurrenceHandle) ||
		    !source->getCompactInstanceSummary(
			occurrenceHandle, occurrenceSummary))
		    continue;
		const BObolViewLodState::CadPayload *bound =
		    sharedViewState.findCadForOccurrence(
			source, occurrenceSummary.sourceInstanceKey);
		if (!bound)
		    continue;
		occurrenceBindings++;
		if (bound->drawMode == BOBOL_LOD_DRAW_SHADED)
		    shadedBindings++;
		if (bound->activeLevel == 0)
		    minimumBindings++;
		if (bound->activeLevel == 4)
		    terminalBindings++;
	    }
	    printf("FAIL: resident shared asset readmission did not preserve "
		   "coverage and advance zero-cost PoP cuts "
		   "(used=%zu cuts=%u faces=%zu payloads=%zu "
		   "occurrences=%zu shaded=%zu min=%zu terminal=%zu)\n",
		   coverageReadmission.getRetainedSceneCostBudgetUsed(),
		   coverageReadmission.getUpdatedCutCount(),
		   sharedViewState.activeFaceCount(),
		   sharedViewState.cadMeshPayloadCount(),
		   occurrenceBindings, shadedBindings, minimumBindings,
		   terminalBindings);
	    ret = 1;
	}
    }

    /*
     * An explicit retained-admission caller may cross the ordinary
     * one-mesh-per-visible-leaf floor (for example, a renderer emergency
     * policy).  The normal controller does not use this as a byte allocator.
     * When requested, the structural occurrence remains as the fallback and
     * the selected occurrence must win the bounded slot regardless of
     * unordered payload-map iteration.
     */
    if (!ret) {
	BObolCompactInstanceHandle selectedHandle;
	BObolCompactInstanceSummary selectedSummary;
	if (!source->getCompactInstanceHandle(
		static_cast<int>(leafCount - 1), selectedHandle) ||
	    !source->getCompactInstanceSummary(
		selectedHandle, selectedSummary)) {
	    printf("FAIL: retained-memory pressure selected fixture lookup\n");
	    ret = 1;
	} else {
	    SoBRLMeshLodSubmitAction memoryRecovery;
	    memoryRecovery.setService(&service);
	    memoryRecovery.setDatabase(
		dbip, "db://many-leaf-test", 2026);
	    memoryRecovery.setViewInfo(&view);
	    memoryRecovery.setViewVolume(&volume, 1.0f);
	    memoryRecovery.setGeneration(service.beginGeneration());
	    memoryRecovery.setRevisions(4, 4);
	    memoryRecovery.setViewLodState(&sharedViewState);
	    memoryRecovery.setRefinementCostBudget(0);
	    memoryRecovery.setRetainedSceneCostBudget(69);
	    memoryRecovery.setPreserveMeshCoverage(FALSE);
	    memoryRecovery.setSubmissionTaskLimit(0);
	    memoryRecovery.setCompactEntryPlan(pinnedPlan);
	    memoryRecovery.apply(source);
	    const BObolViewLodState::CadPayload *selectedPayload =
		sharedViewState.findCadForOccurrence(
		    source, selectedSummary.sourceInstanceKey);
	    if (memoryRecovery.getRetainedSceneCostBudgetUsed() != 69 ||
		memoryRecovery.getUpdatedCutCount() != leafCount - 1 ||
		memoryRecovery.getRefinementBudgetBlockedCount() == 0 ||
		sharedViewState.activeFaceCount() != 4 ||
		sharedViewState.cadMeshPayloadCount() != 1 ||
		!selectedPayload) {
		printf("FAIL: retained-memory pressure did not preserve the "
		       "selected occurrence and retire lower-priority payloads "
		       "(used=%zu cuts=%u blocked=%u faces=%zu payloads=%zu "
		       "selected=%d)\n",
		       memoryRecovery.getRetainedSceneCostBudgetUsed(),
		       memoryRecovery.getUpdatedCutCount(),
		       memoryRecovery.getRefinementBudgetBlockedCount(),
		       sharedViewState.activeFaceCount(),
		       sharedViewState.cadMeshPayloadCount(),
		       selectedPayload ? 1 : 0);
		ret = 1;
	    }
	}
    }

    /*
     * Active selection and erase frontiers must classify only the appended
     * tail of a streaming registry.  Besides protecting the semantics, this
     * guards against reintroducing the accumulated N-entry mask rebuild that
     * made batched cold realization quadratic.
     */
    if (!ret) {
	const std::vector<SbString> selectedPaths = {
	    SbString("many-leaf-root.c/late")
	};
	const std::vector<SbString> visibilityPaths = selectedPaths;
	const std::vector<SbBool> visibilityStates = {FALSE};
	if (source->syncCompactInstanceSelectedPaths(selectedPaths) <= 0 ||
	    source->setCompactInstanceVisibilityOverrides(
		visibilityPaths, visibilityStates) <= 0) {
	    printf("FAIL: streamed frontier fixture setup\n");
	    ret = 1;
	} else {
	    BObolCompactOccurrence selectedTail = occurrences.front();
	    selectedTail.summary.path =
		"many-leaf-root.c/late/selected.bot";
	    selectedTail.sourceMeshRequest.path =
		selectedTail.summary.path;
	    selectedTail.summary.sourceId = 9001;
	    selectedTail.occurrenceIndex =
		static_cast<uint32_t>(leafCount);
	    BObolCompactOccurrence ordinaryTail = occurrences.front();
	    ordinaryTail.summary.path =
		"many-leaf-root.c/other/ordinary.bot";
	    ordinaryTail.sourceMeshRequest.path =
		ordinaryTail.summary.path;
	    ordinaryTail.summary.sourceId = 9002;
	    ordinaryTail.occurrenceIndex =
		static_cast<uint32_t>(leafCount + 1);
	    const std::vector<BObolCompactOccurrence> tail = {
		selectedTail, ordinaryTail
	    };
	    BObolCompactInstanceHandle selectedHandle;
	    BObolCompactInstanceHandle ordinaryHandle;
	    BObolCompactInstanceSummary selectedSummary;
	    BObolCompactInstanceSummary ordinarySummary;
	    if (source->mergeCompactOccurrences(tail) != 2 ||
		source->getCompactInstanceCount() !=
		    static_cast<int>(leafCount + 2) ||
		!source->getCompactInstanceHandle(
		    static_cast<int>(leafCount), selectedHandle) ||
		!source->getCompactInstanceHandle(
		    static_cast<int>(leafCount + 1), ordinaryHandle) ||
		!source->getCompactInstanceSummary(
		    selectedHandle, selectedSummary) ||
		!source->getCompactInstanceSummary(
		    ordinaryHandle, ordinarySummary) ||
		!selectedSummary.selected || selectedSummary.visible ||
		ordinarySummary.selected || !ordinarySummary.visible) {
		printf("FAIL: streamed compact tail did not inherit active "
		       "selection/visibility frontiers sparsely\n");
		ret = 1;
	    }
	}
    }

    /* A coarse, visually dominant part must lead an otherwise equal retained
     * frontier.  Keep this fixture isolated from the shared-instance tests
     * above: real vehicle scenes mix tiny hardware with wheels, rotors, and
     * skins whose projected error needs a hard perceptual floor. */
    if (!ret) {
	SoBRLDatabaseSource *prioritySource = new SoBRLDatabaseSource;
	prioritySource->ref();
	prioritySource->setDatabase(dbip);
	prioritySource->path = "priority-root.c";
	prioritySource->instanceKey = "priority-source";
	prioritySource->sourceRevision = 902;
	prioritySource->lodBotThreshold = 1;
	prioritySource->representationMode =
	    SoBRLDatabaseSource::REPRESENTATION_SHADED;

	BObolCompactOccurrence small = occurrences.front();
	small.summary.path = "priority-root.c/small.bot";
	small.summary.selected = FALSE;
	small.summary.highlighted = FALSE;
	small.sourceMeshRequest.path = small.summary.path;
	small.sourceMeshRequest.meshAssetPath = "lod-submit.bot";
	small.localTransform.makeIdentity();
	small.sourceMeshRequest.meshAssetTransform.makeIdentity();
	BObolCompactOccurrence prominent = small;
	prominent.summary.path = "priority-root.c/prominent.bot";
	prominent.sourceMeshRequest.path = prominent.summary.path;
	prominent.sourceMeshRequest.meshAssetTransform.setScale(
	    SbVec3f(6.0f, 6.0f, 6.0f));
	const std::vector<BObolCompactOccurrence> priorityOccurrences = {
	    small, prominent
	};
	if (prioritySource->setCompactOccurrenceRegistry(
		priorityOccurrences) != 2) {
	    printf("FAIL: perceptual priority registry setup\n");
	    ret = 1;
	} else {
	    point_t priorityPoints[4];
	    VSET(priorityPoints[0], 0.0, 0.0, 0.0);
	    VSET(priorityPoints[1], 1.0, 0.0, 0.0);
	    VSET(priorityPoints[2], 0.0, 1.0, 0.0);
	    VSET(priorityPoints[3], 0.0, 0.0, 1.0);
	    int priorityFaces[12] = {
		0, 1, 2, 0, 3, 1, 1, 3, 2, 2, 3, 0
	    };
	    struct BObolMeshLodData priorityData = {};
	    priorityData.faces = priorityFaces;
	    priorityData.face_count = 4;
	    priorityData.points = priorityPoints;
	    priorityData.point_count = 4;
	    priorityData.points_orig = priorityPoints;
	    priorityData.point_orig_count = 4;
	    VSET(priorityData.bmin, 0.0, 0.0, 0.0);
	    VSET(priorityData.bmax, 1.0, 1.0, 1.0);
	    struct BObolMeshLodHierarchyInfo priorityHierarchy =
		BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
	    priorityHierarchy.min_level = 0;
	    priorityHierarchy.max_level = 1;
	    priorityHierarchy.resident_level = 1;
	    priorityHierarchy.face_count[0] = 1;
	    priorityHierarchy.face_count[1] = 4;
	    priorityHierarchy.point_count[0] = 3;
	    priorityHierarchy.point_count[1] = 4;
	    VSET(priorityHierarchy.quantization_min, 0.0, 0.0, 0.0);
	    VSET(priorityHierarchy.quantization_max, 1.0, 1.0, 1.0);
	    BObolLodProgressiveMeshPtr priorityMesh(
		new BObolLodProgressiveMesh);
	    BObolViewLodState priorityState;
	    bool priorityReady = priorityMesh->update(
		priorityData, priorityHierarchy, 1, FALSE) ? true : false;
	    for (int priorityIndex = 0;
		 priorityReady && priorityIndex < 2; priorityIndex++) {
		BObolCompactInstanceHandle priorityHandle;
		BObolCompactInstanceSummary prioritySummary;
		if (!prioritySource->getCompactInstanceHandle(
			priorityIndex, priorityHandle) ||
		    !prioritySource->getCompactInstanceSummary(
			priorityHandle, prioritySummary)) {
		    priorityReady = false;
		    break;
		}
		BObolLodResult priorityResult;
		priorityResult.request.databaseId = "db://priority-test";
		priorityResult.request.databaseRevision = 2026;
		priorityResult.request.sourceRevision =
		    prioritySource->sourceRevision.getValue();
		priorityResult.request.objectPath =
		    prioritySummary.meshAssetPath;
		priorityResult.request.objectName =
		    prioritySummary.meshAssetName;
		priorityResult.request.occurrenceKey =
		    prioritySummary.sourceInstanceKey;
		priorityResult.request.viewRevision = 6;
		priorityResult.request.policyRevision = 6;
		priorityResult.request.drawMode = BOBOL_LOD_DRAW_SHADED;
		priorityResult.request.providerId = "bobol_mesh_lod";
		priorityResult.request.providerVersion =
		    BOBOL_MESH_LOD_PROVIDER_VERSION;
		priorityResult.request.requestedLevel = 1;
		priorityResult.request.bounds = priorityMesh->bounds();
		priorityResult.cacheKey =
		    bobol_lod_cache_key(priorityResult.request);
		priorityResult.geometry.kind =
		    BOBOL_LOD_GEOMETRY_MESH_LOD_CACHE;
		priorityResult.geometry.activeLevel = 0;
		priorityResult.progressiveMesh = priorityMesh;
		priorityResult.residentLevel = 1;
		priorityResult.resultKind = BOBOL_LOD_RESULT_MESH;
		priorityResult.providerStatus = BOBOL_LOD_PROVIDER_READY;
		priorityResult.bounds = priorityMesh->bounds();
		priorityResult.counts.faceCount = 1;
		priorityResult.counts.pointCount = 3;
		priorityResult.counts.originalPointCount = 3;
		priorityResult.terminal = FALSE;
		priorityReady = priorityState.applySourceResult(
		    prioritySource, priorityResult) ? true : false;
	    }
	    if (!priorityReady) {
		printf("FAIL: perceptual priority retained payload setup\n");
		ret = 1;
	    } else {
		SoBRLMeshLodSubmitAction priorityAction;
		priorityAction.setService(&service);
		priorityAction.setDatabase(
		    dbip, "db://priority-test", 2026);
		priorityAction.setViewInfo(&view);
		priorityAction.setViewVolume(&volume, 1.0f);
		priorityAction.setGeneration(service.beginGeneration());
		priorityAction.setRevisions(6, 6);
		priorityAction.setViewLodState(&priorityState);
		priorityAction.setRefinementCostBudget(0);
		priorityAction.setSubmissionTaskLimit(0);
		priorityAction.apply(prioritySource);
		std::vector<size_t> priorityPlan;
		priorityAction.getCompactEntryPlan(priorityPlan);
		if (priorityPlan.size() != 2 || priorityPlan.front() != 1) {
		    printf("FAIL: prominent coarse occurrence did not lead the "
			   "perceptual frontier (plan=%zu first=%zu)\n",
			   priorityPlan.size(),
			   priorityPlan.empty() ? SIZE_MAX :
				priorityPlan.front());
		    ret = 1;
		}
	    }
	}
	prioritySource->unref();
    }

    camera->unref();
    service.stop();
    source->unref();
    bobol_mesh_lod_cache_clear_database(dbip);
    db_close(dbip);
    bu_file_delete(dbpath);
    return ret;
}

static int
test_compact_native_mesh_not_pop_submitted(void)
{
    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;
    if (make_submit_test_db(dbpath, sizeof(dbpath), &dbip))
	return 1;

    SoBRLViewLodGroup *root = new SoBRLViewLodGroup;
    root->ref();
    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = "native-analytic.tgc";
    source->instanceKey = "native-analytic-instance";
    source->lodBotThreshold = 1;
    source->representationMode = SoBRLDatabaseSource::REPRESENTATION_SHADED;
    root->addChild(source);

    BObolCompactOccurrence occurrence;
    occurrence.geometry = compact_native_triangle_geometry();
    occurrence.summary.valid = TRUE;
    occurrence.summary.shapeKind =
	BObolRealizedShapeSummary::SHAPE_MESH;
    occurrence.summary.path = "root.c/native-analytic.tgc";
    occurrence.summary.sourceName = "native-analytic.tgc";
    occurrence.summary.sourceType = "tgc";
    occurrence.summary.sourceId = 711;
    occurrence.summary.visible = TRUE;
    occurrence.summary.selectable = TRUE;
    occurrence.lodBacked = FALSE;
    occurrence.sourceMeshRequestValid = FALSE;
    std::vector<BObolCompactOccurrence> occurrences(1, occurrence);

    int ret = 0;
    if (source->setCompactOccurrenceRegistry(occurrences) != 1) {
	printf("FAIL: compact native mesh setup\n");
	ret = 1;
    }

    BObolLodService service;
    if (!ret && !service.start(1, TRUE)) {
	printf("FAIL: compact native mesh service did not start\n");
	ret = 1;
    }

    if (!ret) {
	SoBRLMeshLodSubmitAction submit;
	submit.setService(&service);
	submit.setDatabase(dbip, "db://compact-native-test", 2026);
	submit.setGeneration(service.beginGeneration());
	submit.setRevisions(71, 72);
	submit.apply(root);
	if (submit.getVisitedMeshCount() != 1 ||
	    submit.getSubmittedTaskCount() != 0 ||
	    submit.getSkippedMeshCount() != 1 ||
	    service.pendingTaskCountForDiagnostics() != 0 ||
	    service.queuedResultCountForDiagnostics() != 0) {
	    printf("FAIL: native analytic mesh was submitted to the BoT PoP provider\n");
	    ret = 1;
	}
    }

    service.stop();
    root->unref();
    bobol_mesh_lod_cache_clear_database(dbip);
    db_close(dbip);
    bu_file_delete(dbpath);
    return ret;
}

/*
 * A growing compact registry is one standing source, not a succession of new
 * scenes.  Exercise the subtle transition which production cold streams use:
 * an initial sub-window inventory completes, a large exact delta starts a
 * bounded scan, and another inventory batch arrives while that delta is
 * pending.  The pending delta must finish its remaining suffix rather than
 * restart at entry zero merely because displayMeshLodRevision advanced.
 */
static int
test_view_controller_compact_append_preserves_submission_cursor(void)
{
    static const size_t initialCount = 512;
    static const size_t deltaCount = 2500;
    static const size_t quietWave = 2048;
    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;
    if (make_submit_test_db(dbpath, sizeof(dbpath), &dbip))
	return 1;

    SoBRLViewLodGroup *root = new SoBRLViewLodGroup;
    root->ref();
    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = "streaming-root.c";
    source->instanceKey = "streaming-root-instance";
    source->lodBotThreshold = 1;
    source->sourceRevision = 2026;
    source->realizedSourceRevision = 2026;
    source->realizedInputsRevision = source->inputsRevision.getValue();
    source->realizedViewRevision = source->viewRevision.getValue();
    /* The source-wide authoritative worker is deliberately still running.
     * Its revision-validated compact contracts must already be plannable. */
    source->realizationStatus = SoBRLDatabaseSource::UNREALIZED;
    source->stale = FALSE;
    source->staleReason = SoBRLDatabaseSource::STALE_NONE;
    root->addChild(source);

    BObolCompactOccurrence prototype;
    prototype.geometry = compact_projected_proxy_geometry();
    prototype.summary.valid = TRUE;
    prototype.summary.shapeKind =
	BObolRealizedShapeSummary::SHAPE_VLIST;
    prototype.summary.sourceName = "lod-submit.bot";
    prototype.summary.sourceType = "bot";
    prototype.summary.geometryKind = "aabb";
    prototype.summary.sourceId = 2026;
    prototype.summary.visible = TRUE;
    prototype.summary.selectable = TRUE;
    prototype.lodBacked = TRUE;
    prototype.sourceMeshRequestValid = TRUE;
    prototype.sourceMeshRequest.sourceName = "lod-submit.bot";
    prototype.sourceMeshRequest.faceCount = 4;
    prototype.sourceMeshRequest.pointCount = 4;
    prototype.sourceMeshRequest.bounds = SbBox3f(
	SbVec3f(0.0f, 0.0f, 0.0f), SbVec3f(1.0f, 1.0f, 1.0f));
    prototype.sourceMeshRequest.meshAssetBounds =
	prototype.sourceMeshRequest.bounds;
    prototype.sourceMeshRequest.meshAssetTransform = SbMatrix::identity();
    /* Keep every fixture occurrence outside the camera.  The test exercises
     * only bounded cursor progression and therefore schedules no provider
     * work or cache writes. */
    prototype.localTransform.setTranslate(
	SbVec3f(100000.0f, 0.0f, 0.0f));

    std::vector<BObolCompactOccurrence> occurrences;
    occurrences.reserve(initialCount);
    for (size_t i = 0; i < initialCount; i++) {
	BObolCompactOccurrence occurrence = prototype;
	SbString path;
	path.sprintf("streaming-root.c/leaf-%zu.bot", i);
	occurrence.summary.path = path;
	occurrence.sourceMeshRequest.path = path;
	occurrence.occurrenceIndex = static_cast<uint32_t>(i);
	occurrences.push_back(std::move(occurrence));
    }

    int ret = 0;
    if (source->setCompactOccurrenceRegistry(occurrences) !=
	    static_cast<int>(initialCount)) {
	printf("FAIL: compact append cursor fixture registry setup\n");
	ret = 1;
    }
    BObolCompactLodInstanceSummary firstSummary;
    size_t firstEntryIndex = SIZE_MAX;
    if (!ret &&
	(!source->getCompactLodInstanceSummary(0, firstSummary) ||
	 firstSummary.sourceInstanceKey.getLength() == 0 ||
	 !source->getCompactInstanceIndex(
	     firstSummary.sourceInstanceKey.getString(), firstEntryIndex) ||
	 firstEntryIndex != 0)) {
	printf("FAIL: compact occurrence key did not resolve to its stable "
	       "source-local index\n");
	ret = 1;
    }

    BObolLodService service;
    if (!ret && !service.start(1, TRUE)) {
	printf("FAIL: compact append cursor fixture service start\n");
	ret = 1;
    }

    SoOrthographicCamera *camera = new SoOrthographicCamera;
    camera->ref();
    camera->height = 10.0f;
    if (!ret) {
	BObolViewController controller(root, camera);
	controller.setViewportSize(800, 600);
	controller.setLodService(&service);

	if (controller.submitLodRequestsIfNeeded() != 0 ||
	    controller.getLastLodVisitedMeshCount() != initialCount ||
	    controller.hasPendingLodSubmissions()) {
	    printf("FAIL: compact append cursor did not complete its initial "
		   "sub-window inventory (visited=%u pending=%d diagnostics=%s)\n",
		   controller.getLastLodVisitedMeshCount(),
		   controller.hasPendingLodSubmissions() ? 1 : 0,
		   controller.getLastLodDiagnostics().getString());
	    ret = 1;
	}

	std::vector<BObolCompactOccurrence> delta;
	delta.reserve(deltaCount);
	for (size_t i = 0; i < deltaCount; i++) {
	    BObolCompactOccurrence occurrence = prototype;
	    SbString path;
	    path.sprintf("streaming-root.c/delta-%zu.bot", i);
	    occurrence.summary.path = path;
	    occurrence.sourceMeshRequest.path = path;
	    occurrence.occurrenceIndex =
		static_cast<uint32_t>(initialCount + i);
	    delta.push_back(std::move(occurrence));
	}
	if (!ret && source->mergeCompactOccurrences(delta) !=
		static_cast<int>(deltaCount)) {
	    printf("FAIL: compact append cursor could not append large delta\n");
	    ret = 1;
	}
	size_t deltaVisited = 0;
	if (!ret) {
	    const int submitted = controller.submitLodRequestsIfNeeded();
	    deltaVisited = controller.getLastLodVisitedMeshCount();
	    if (submitted != 0 || deltaVisited == 0 ||
		deltaVisited > quietWave ||
		!controller.hasPendingLodSubmissions()) {
	    printf("FAIL: compact append cursor did not stop at the first "
		   "bounded delta wave (visited=%u pending=%d diagnostics=%s)\n",
		   controller.getLastLodVisitedMeshCount(),
		   controller.hasPendingLodSubmissions() ? 1 : 0,
		   controller.getLastLodDiagnostics().getString());
	    ret = 1;
	    }
	}

	BObolCompactOccurrence appended = prototype;
	appended.summary.path = "streaming-root.c/appended.bot";
	appended.sourceMeshRequest.path = appended.summary.path;
	appended.occurrenceIndex =
	    static_cast<uint32_t>(initialCount + deltaCount);
	if (!ret && source->mergeCompactOccurrences(
		std::vector<BObolCompactOccurrence>(1, appended)) != 1) {
	    printf("FAIL: compact append cursor could not append inventory\n");
	    ret = 1;
	}

	/* The exact remaining delta plus the one new tail entry proves the
	 * cursor continued and the pinned plan was extended in place.  The
	 * production action is bounded by both a nominal 2048-entry wave and a
	 * wall-clock deadline, so a loaded test host may need more than one
	 * call.  Assert lossless aggregate progress rather than assuming the
	 * timing limit can never win. */
	const size_t completeDelta = deltaCount + 1;
	for (size_t attempt = 0;
	     !ret && deltaVisited < completeDelta && attempt < 8;
	     ++attempt) {
	    if (controller.submitLodRequestsIfNeeded() != 0) {
		printf("FAIL: compact append inventory submitted unexpected "
		       "provider work\n");
		ret = 1;
		break;
	    }
	    const size_t visited =
		controller.getLastLodVisitedMeshCount();
	    if (visited == 0 || visited > completeDelta - deltaVisited) {
		printf("FAIL: compact append inventory restarted or stalled the "
		       "submission cursor (visited=%zu aggregate=%zu expected=%zu "
		       "pending=%d diagnostics=%s)\n",
		       visited, deltaVisited, completeDelta,
		       controller.hasPendingLodSubmissions() ? 1 : 0,
		       controller.getLastLodDiagnostics().getString());
		ret = 1;
		break;
	    }
	    deltaVisited += visited;
	}
	if (!ret &&
	    (deltaVisited != completeDelta ||
	     !controller.hasPendingLodSubmissions())) {
	    printf("FAIL: compact append inventory did not complete its exact "
		   "extended delta (visited=%zu expected=%zu pending=%d "
		   "diagnostics=%s)\n",
		   deltaVisited, completeDelta,
		   controller.hasPendingLodSubmissions() ? 1 : 0,
		   controller.getLastLodDiagnostics().getString());
	    ret = 1;
	}
	/*
	 * Completing the selective append plan must not make that subset the
	 * controller's authoritative coverage proof.  It is possible for older
	 * structural ranges to have missed an earlier scan while the registry
	 * was growing.  The delta is followed by one bounded full-inventory pass
	 * before the controller may report idle.
	 */
	size_t fullRescanVisited = 0;
	for (size_t attempt = 0;
	     !ret && controller.hasPendingLodSubmissions() && attempt < 8;
	     ++attempt) {
	    if (controller.submitLodRequestsIfNeeded() != 0) {
		printf("FAIL: compact append full coverage rescan submitted "
		       "unexpected provider work\n");
		ret = 1;
		break;
	    }
	    fullRescanVisited += controller.getLastLodVisitedMeshCount();
	}
	const size_t completeInventory =
	    initialCount + deltaCount + 1;
	if (!ret &&
	    (controller.hasPendingLodSubmissions() ||
	     fullRescanVisited != completeInventory)) {
	    printf("FAIL: compact append delta did not finish one complete "
		   "coverage rescan (visited=%zu expected=%zu pending=%d "
		   "diagnostics=%s)\n",
		   fullRescanVisited, completeInventory,
		   controller.hasPendingLodSubmissions() ? 1 : 0,
		   controller.getLastLodDiagnostics().getString());
	    ret = 1;
	}
	controller.setLodService(NULL);
    }

    camera->unref();
    service.stop();
    root->unref();
    bobol_mesh_lod_cache_clear_database(dbip);
    db_close(dbip);
    bu_file_delete(dbpath);
    return ret;
}

static int
test_compact_aabb_stream_upgrade(void)
{
    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->ref();
    source->path = "stream-root.c";
    source->instanceKey = "stream-root-instance";

    BObolCompactOccurrence proxy;
    proxy.geometry = compact_projected_proxy_geometry();
    proxy.summary.valid = TRUE;
    proxy.summary.shapeKind = BObolRealizedShapeSummary::SHAPE_MESH;
    proxy.summary.path = "stream-root.c/leaf.bot";
    proxy.summary.sourceName = "leaf.bot";
    proxy.summary.sourceType = "bot";
    proxy.summary.geometryKind = "aabb";
    proxy.summary.sourceId = 101;
    proxy.summary.visible = TRUE;
    proxy.summary.selectable = TRUE;
    proxy.lodBacked = FALSE;
    proxy.sourceMeshRequestValid = FALSE;
    proxy.sourceMeshRequest.path = proxy.summary.path;
    proxy.sourceMeshRequest.sourceName = proxy.summary.sourceName;
    proxy.sourceMeshRequest.faceCount = 4096;
    proxy.sourceMeshRequest.pointCount = 2048;
    proxy.sourceMeshRequest.bounds = SbBox3f(
	SbVec3f(0.0f, 0.0f, 0.0f), SbVec3f(1.0f, 1.0f, 1.0f));
    proxy.sourceMeshRequest.meshAssetBounds =
	proxy.sourceMeshRequest.bounds;
    proxy.geometryTransform.setScale(SbVec3f(2.0f, 3.0f, 4.0f));
    SbMatrix proxyGeometryTranslation;
    proxyGeometryTranslation.setTranslate(SbVec3f(1.0f, 2.0f, 3.0f));
    proxy.geometryTransform.multRight(proxyGeometryTranslation);
    proxy.localTransform.setTranslate(SbVec3f(10.0f, 20.0f, 30.0f));

    std::vector<BObolCompactOccurrence> initial(1, proxy);
    int ret = 0;
    BObolCompactInstanceHandle initialHandle;
    if (source->setCompactOccurrenceRegistry(initial) != 1 ||
	!source->getCompactInstanceHandle(0, initialHandle)) {
	printf("FAIL: compact AABB stream upgrade fixture setup\n");
	ret = 1;
    }

    BObolCompactOccurrence contractedProxy = proxy;
    const std::shared_ptr<const Obol::PartGeometry> contractedGeometry =
	compact_native_triangle_geometry();
    contractedProxy.geometry = contractedGeometry;
    contractedProxy.geometryTransform = SbMatrix::identity();
    contractedProxy.lodBacked = TRUE;
    contractedProxy.sourceMeshRequestValid = TRUE;
    contractedProxy.summary.sourceId = 102;
    if (!ret && source->mergeCompactOccurrences(
	    std::vector<BObolCompactOccurrence>(1, contractedProxy)) != 1) {
	printf("FAIL: compact streamed AABB did not accept its later PoP contract\n");
	ret = 1;
    }

    BObolCompactOccurrence contractedOccurrence;
    BObolCompactLodInstanceSummary lodSummary;
    SbBox3f proxyBounds;
    if (!ret &&
	(!source->getCompactOccurrence(0, contractedOccurrence) ||
	 !source->getCompactLodInstanceSummary(0, lodSummary) ||
	 !contractedOccurrence.lodBacked ||
	 !contractedOccurrence.sourceMeshRequestValid ||
	 contractedOccurrence.geometry != proxy.geometry ||
	 contractedOccurrence.geometry == contractedGeometry ||
	 !contractedOccurrence.geometryTransform.equals(
	     proxy.geometryTransform, 0.000001f) ||
	 !contractedOccurrence.localTransform.equals(
	     proxy.localTransform, 0.000001f) ||
	 !lodSummary.localToSource.equals(
	     proxy.localTransform, 0.000001f) ||
	 !lodSummary.localBounds.getMin().equals(
	     SbVec3f(0.0f, 0.0f, 0.0f), 0.0001f) ||
	 !lodSummary.localBounds.getMax().equals(
	     SbVec3f(1.0f, 1.0f, 1.0f), 0.0001f) ||
	 !source->getSourceBounds(proxyBounds) ||
	 !proxyBounds.getMin().equals(SbVec3f(11.0f, 22.0f, 33.0f),
	     0.0001f) ||
	 !proxyBounds.getMax().equals(SbVec3f(13.0f, 25.0f, 37.0f),
	     0.0001f) ||
	 !BU_STR_EQUAL(
	     contractedOccurrence.summary.geometryKind.getString(), "aabb") ||
	 source->getCompactInstanceCount() != 1)) {
	printf("FAIL: compact same-tier PoP contract/transform was not retained\n");
	ret = 1;
    }

    /*
     * Initialize the retained view presentation while the structural proxy
     * is current.  The later richer occurrence must update this existing
     * assembly, not merely the compact-index metadata.
     */
    BObolViewLodState presentationState;
    std::vector<const BObolViewLodState::CadPayload *> noPayloads;
    SoCADAssembly *retainedPresentation = NULL;
    std::vector<Obol::InstanceId> retainedIds;
    if (!ret) {
	(void)source->compactViewLodAssembly(
	    noPayloads, &presentationState);
	retainedPresentation =
	    presentationState.findCadPresentation(source);
	retainedIds = retainedPresentation ?
	    retainedPresentation->instanceIds() :
	    std::vector<Obol::InstanceId>();
	const std::optional<Obol::InstanceRecord> retainedRecord =
	    retainedIds.size() == 1 ?
	    retainedPresentation->getInstanceRecord(retainedIds[0]) :
	    std::optional<Obol::InstanceRecord>();
	const Obol::PartGeometry *retainedGeometry = retainedRecord ?
	    retainedPresentation->partGeometry(retainedRecord->part) : NULL;
	if (!retainedPresentation || retainedIds.size() != 1 ||
	    !retainedGeometry || !retainedGeometry->structuralProxy) {
	    printf("FAIL: compact AABB upgrade did not initialize its retained "
		   "structural presentation\n");
	    ret = 1;
	}
    }

    BObolCompactOccurrence mesh = contractedProxy;
    mesh.geometry = compact_native_triangle_geometry();
    mesh.geometryTransform = SbMatrix::identity();
    mesh.summary.geometryKind = "surface";
    mesh.summary.sourceId = 103;
    mesh.lodBacked = TRUE;
    if (!ret && source->mergeCompactOccurrences(
	    std::vector<BObolCompactOccurrence>(1, mesh)) != 1) {
	printf("FAIL: compact streamed mesh did not supersede its AABB tier\n");
	ret = 1;
    }

    BObolCompactInstanceHandle upgradedHandle;
    BObolCompactInstanceSummary upgradedSummary;
    BObolCompactOccurrence upgradedOccurrence;
    if (!ret &&
	(!source->getCompactInstanceHandle(0, upgradedHandle) ||
	 !source->getCompactInstanceSummary(upgradedHandle, upgradedSummary) ||
	 !source->getCompactOccurrence(0, upgradedOccurrence) ||
	 upgradedHandle.sourceNodeId != initialHandle.sourceNodeId ||
	 upgradedHandle.instanceWord0 != initialHandle.instanceWord0 ||
	 upgradedHandle.instanceWord1 != initialHandle.instanceWord1 ||
	 source->getCompactInstanceCount() != 1 ||
	 !upgradedSummary.meshGeometry || upgradedSummary.wireGeometry ||
	 !upgradedOccurrence.geometry ||
	 !upgradedOccurrence.geometry->shaded ||
	 !upgradedOccurrence.geometryTransform.equals(
	     SbMatrix::identity(), 0.000001f) ||
	 !upgradedOccurrence.localTransform.equals(
	     proxy.localTransform, 0.000001f) ||
	 !BU_STR_EQUAL(upgradedOccurrence.summary.geometryKind.getString(),
	     "surface"))) {
	printf("FAIL: compact AABB upgrade did not preserve identity and install mesh geometry\n");
	ret = 1;
    }
    if (!ret) {
	(void)source->compactViewLodAssembly(
	    noPayloads, &presentationState);
	SoCADAssembly *upgradedPresentation =
	    presentationState.findCadPresentation(source);
	const std::optional<Obol::InstanceRecord> upgradedRecord =
	    retainedIds.size() == 1 && upgradedPresentation ?
	    upgradedPresentation->getInstanceRecord(retainedIds[0]) :
	    std::optional<Obol::InstanceRecord>();
	const Obol::PartGeometry *presentedGeometry = upgradedRecord ?
	    upgradedPresentation->partGeometry(upgradedRecord->part) : NULL;
	if (upgradedPresentation != retainedPresentation ||
	    !presentedGeometry ||
	    presentedGeometry != upgradedOccurrence.geometry.get() ||
	    presentedGeometry->structuralProxy ||
	    !presentedGeometry->shaded) {
	    printf("FAIL: compact AABB metadata upgraded without replacing "
		   "the retained structural presentation\n");
	    ret = 1;
	}
    }

    if (!ret && source->mergeCompactOccurrences(initial) != 0) {
	printf("FAIL: late compact AABB displaced an installed mesh tier\n");
	ret = 1;
    }

    /*
     * Moving the camera while a detached realization finishes must not make
     * its source-valid native/PoP registry stale.  The view selects active
     * PoP prefixes after adoption; it is not part of geometry identity.
     */
    SoBRLDatabaseSource *detached = new SoBRLDatabaseSource;
    detached->ref();
    detached->path = source->path.getValue();
    detached->instanceKey = source->instanceKey.getValue();
    detached->sourceRevision = source->sourceRevision.getValue();
    detached->inputsRevision = source->inputsRevision.getValue();
    detached->drawMode = source->drawMode.getValue();
    detached->representationMode = source->representationMode.getValue();
    detached->viewRevision = 17;
    source->viewRevision = 18;
    source->markStale(SoBRLDatabaseSource::STALE_SOURCE);
    if (!ret && source->hasDisplayMeshLodRequests()) {
	printf("FAIL: source invalidation did not revoke the compact "
	       "mesh-request epoch\n");
	ret = 1;
    }
    if (!ret &&
	(detached->setCompactOccurrenceRegistry(
	     std::vector<BObolCompactOccurrence>(1, mesh)) != 1 ||
	 source->adoptDetachedCompactRealization(detached) != 1)) {
	printf("FAIL: camera change rejected source-valid detached compact geometry\n");
	ret = 1;
    }
    if (!ret &&
	(!source->getCompactInstanceHandle(0, upgradedHandle) ||
	 !source->getCompactInstanceSummary(upgradedHandle, upgradedSummary) ||
	 !upgradedSummary.meshGeometry ||
	 !source->hasDisplayMeshLodRequests() ||
	 source->getDisplayMeshLodRequestCount() != 1 ||
	 !BU_STR_EQUAL(upgradedSummary.geometryKind.getString(), "surface"))) {
	printf("FAIL: detached compact adoption lost mesh geometry or failed "
	       "to restore its mesh-request epoch\n");
	ret = 1;
    }
    detached->unref();

    /*
     * A cold stream publishes a synthetic whole-target extent before its
     * leaves.  Final authoritative adoption must retire that one retained
     * instance as a sparse visibility update.  In particular, the compact
     * index may already say "hidden" while a view still has the earlier
     * visible generation, so adoption must publish a new visibility revision
     * rather than treating the matching boolean as proof of presentation.
     */
    SoBRLDatabaseSource *retirementSource = new SoBRLDatabaseSource;
    retirementSource->ref();
    retirementSource->path = "retirement-root.c";
    retirementSource->instanceKey = "retirement-root-instance";

    BObolCompactOccurrence overview = proxy;
    overview.summary.path = "retirement-root.c/@draw-extent";
    overview.summary.sourceName = "retirement-root.c";
    overview.summary.sourceType = "proxy";
    overview.summary.geometryKind = "overview-aabb";
    overview.summary.recordRole = "lod-overview";
    overview.summary.selectable = FALSE;
    overview.lodBacked = FALSE;
    overview.sourceMeshRequestValid = FALSE;

    BObolCompactOccurrence retirementMesh = mesh;
    retirementMesh.summary.path = "retirement-root.c/leaf.bot";
    retirementMesh.summary.sourceName = "leaf.bot";
    retirementMesh.summary.recordRole = "";
    retirementMesh.summary.sourceId = 201;
    std::vector<BObolCompactOccurrence> streamed;
    streamed.push_back(overview);
    streamed.push_back(retirementMesh);
    if (!ret && retirementSource->setCompactOccurrenceRegistry(streamed) != 2) {
	printf("FAIL: compact overview retirement fixture setup\n");
	ret = 1;
    }

    BObolViewLodState retirementState;
    SoCADAssembly *retirementPresentation = NULL;
    Obol::InstanceId overviewId;
    Obol::InstanceId meshId;
    if (!ret) {
	(void)retirementSource->compactViewLodAssembly(
	    noPayloads, &retirementState);
	retirementPresentation =
	    retirementState.findCadPresentation(retirementSource);
	const std::vector<Obol::InstanceId> ids = retirementPresentation ?
	    retirementPresentation->instanceIds() :
	    std::vector<Obol::InstanceId>();
	for (const Obol::InstanceId &id : ids) {
	    const std::optional<Obol::InstanceRecord> record =
		retirementPresentation->getInstanceRecord(id);
	    const Obol::PartGeometry *geometry = record ?
		retirementPresentation->partGeometry(record->part) : NULL;
	    if (geometry && geometry->structuralProxy)
		overviewId = id;
	    else if (geometry)
		meshId = id;
	}
	if (!retirementPresentation || ids.size() != 2 ||
	    !overviewId.isValid() || !meshId.isValid() ||
	    retirementPresentation->isInstanceHidden(overviewId)) {
	    printf("FAIL: compact overview was not initially retained visible\n");
	    ret = 1;
	}
    }

    SoBRLDatabaseSource *retirementDetached = new SoBRLDatabaseSource;
    retirementDetached->ref();
    retirementDetached->path = retirementSource->path.getValue();
    retirementDetached->instanceKey =
	retirementSource->instanceKey.getValue();
    retirementDetached->sourceRevision =
	retirementSource->sourceRevision.getValue();
    retirementDetached->inputsRevision =
	retirementSource->inputsRevision.getValue();
    retirementDetached->drawMode = retirementSource->drawMode.getValue();
    retirementDetached->representationMode =
	retirementSource->representationMode.getValue();
    if (!ret &&
	(retirementDetached->setCompactOccurrenceRegistry(
	     std::vector<BObolCompactOccurrence>(1, retirementMesh)) != 1 ||
	 retirementSource->adoptDetachedCompactRealization(
	     retirementDetached, TRUE) != 1)) {
	printf("FAIL: compact authoritative overview retirement adoption\n");
	ret = 1;
    }
    if (!ret) {
	(void)retirementSource->compactViewLodAssembly(
	    noPayloads, &retirementState);
	SoCADAssembly *retiredPresentation =
	    retirementState.findCadPresentation(retirementSource);
	if (retiredPresentation != retirementPresentation ||
	    !retiredPresentation->isInstanceHidden(overviewId) ||
	    retiredPresentation->isInstanceHidden(meshId) ||
	    retiredPresentation->instanceCount() != 2) {
	    printf("FAIL: compact overview retirement did not sparsely hide "
		   "the retained extent\n");
	    ret = 1;
	}
    }
    /*
     * A root redraw republishes the user-facing draw path as visible.  That
     * state must restore real leaves without resurrecting the internal
     * overview record which authoritative adoption just retired.
     */
    if (!ret &&
	retirementSource->setCompactInstanceDisplayStateForPath(
	    "retirement-root.c", TRUE, 1, TRUE, 0, FALSE, 0, FALSE) < 0) {
	printf("FAIL: compact overview retirement redraw setup\n");
	ret = 1;
    }
    if (!ret) {
	(void)retirementSource->compactViewLodAssembly(
	    noPayloads, &retirementState);
	SoCADAssembly *redrawnPresentation =
	    retirementState.findCadPresentation(retirementSource);
	if (redrawnPresentation != retirementPresentation ||
	    !redrawnPresentation->isInstanceHidden(overviewId) ||
	    redrawnPresentation->isInstanceHidden(meshId)) {
	    printf("FAIL: root redraw resurrected a retired compact overview\n");
	    ret = 1;
	}
    }
    retirementDetached->unref();
    retirementSource->unref();

    source->unref();
    return ret;
}

static int
test_compact_mesh_lod_projection_and_mode_parity(void)
{
    char cacheDir[MAXPATHLEN] = {0};
    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;
    int ret = 0;

    bu_dir(cacheDir, MAXPATHLEN, BU_DIR_CURR,
	   "bobol_compact_lod_projection_cache", NULL);
    bu_dirclear(cacheDir);
    bu_mkdir(cacheDir);
    bu_setenv("BU_DIR_CACHE", cacheDir, 1);
    if (make_submit_test_db(dbpath, sizeof(dbpath), &dbip)) {
	bu_dirclear(cacheDir);
	return 1;
    }

    BObolViewLodState viewState;
    SoBRLViewLodGroup *root = new SoBRLViewLodGroup;
    root->ref();
    root->setViewLodState(&viewState);

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = "lod-submit.bot";
    source->instanceKey = "compact-projected-source";
    source->sourceRevision = 515;
    source->lodBotThreshold = 1;
    source->representationMode = SoBRLDatabaseSource::REPRESENTATION_SHADED;
    source->selectedColor = SbColor(0.25f, 0.75f, 0.125f);
    source->drawMatrixValid = TRUE;
    SbMatrix placement;
    placement.setTranslate(SbVec3f(100.0f, 0.0f, 0.0f));
    source->drawMatrix = placement;
    root->addChild(source);

    BObolCompactOccurrence occurrence;
    occurrence.geometry = compact_projected_proxy_geometry();
    occurrence.summary.valid = TRUE;
    occurrence.summary.shapeKind =
	BObolRealizedShapeSummary::SHAPE_VLIST;
    occurrence.summary.path = "root.c/lod-submit.bot";
    occurrence.summary.sourceName = "lod-submit.bot";
    occurrence.summary.sourceType = "bot";
    occurrence.summary.sourceId = 515;
    occurrence.summary.visible = TRUE;
    occurrence.summary.selectable = TRUE;
    occurrence.lodBacked = TRUE;
    occurrence.sourceMeshRequestValid = TRUE;
    occurrence.sourceMeshRequest.path = occurrence.summary.path;
    occurrence.sourceMeshRequest.sourceName = "lod-submit.bot";
    occurrence.sourceMeshRequest.faceCount = 4;
    occurrence.sourceMeshRequest.pointCount = 4;
    occurrence.sourceMeshRequest.bounds = SbBox3f(
	SbVec3f(0.0f, 0.0f, 0.0f), SbVec3f(1.0f, 1.0f, 1.0f));
    occurrence.sourceMeshRequest.meshAssetBounds =
	occurrence.sourceMeshRequest.bounds;
    occurrence.sourceMeshRequest.meshAssetTransform = SbMatrix::identity();
    /* Model cold coverage as a normalized proxy whose local coordinates are
     * intentionally different from the source BoT arrays.  The view-LoD
     * presentation must discard this scale when it swaps in the cached mesh. */
    occurrence.geometryTransform.setScale(SbVec3f(7.0f, 11.0f, 13.0f));
    std::vector<BObolCompactOccurrence> occurrences(1, occurrence);
    if (source->setCompactOccurrenceRegistry(occurrences) != 1) {
	printf("FAIL: compact projected LoD occurrence setup\n");
	ret = 1;
    }

    BObolLodService service;
    if (!ret && !service.start(1, TRUE)) {
	printf("FAIL: compact projected LoD service did not start\n");
	ret = 1;
    }

    SoOrthographicCamera *camera = new SoOrthographicCamera;
    camera->ref();
    camera->position = SbVec3f(100.5f, 0.5f, 5.0f);
    camera->height = 10.0f;
    const SbViewVolume volume = camera->getViewVolume(640.0f / 480.0f);
    struct bv_view_info view = BV_VIEW_INFO_INIT;
    view.width = 640;
    view.height = 480;
    view.size = 10.0;

    std::vector<BObolLodResult> results;
    if (!ret) {
	SoBRLMeshLodSubmitAction submit;
	submit.setService(&service);
	submit.setDatabase(dbip, "db://compact-projected-test", 2026);
	submit.setViewInfo(&view);
	submit.setViewVolume(&volume, 1.0f);
	submit.setGeneration(service.beginGeneration());
	submit.setRevisions(61, 62);
	submit.apply(root);
	if (submit.getVisitedMeshCount() != 1 ||
	    submit.getSubmittedTaskCount() != 1 ||
	    wait_for_service(service) ||
	    service.drainResults(results) != 1 || results.size() != 1 ||
	    results[0].request.projectedPixelDiameter < 47.0f ||
	    results[0].request.projectedPixelDiameter > 49.0f ||
	    results[0].request.projectedPixelArea < 2200.0f ||
	    results[0].request.projectedPixelArea > 2400.0f ||
	    results[0].request.projectedPixelPerimeter < 190.0f ||
	    results[0].request.projectedPixelPerimeter > 194.0f ||
	    results[0].request.requestedLevel != 6 ||
	    !results[0].progressiveMesh ||
	    results[0].geometry.activeLevel !=
		results[0].progressiveMesh->minimumLevel() ||
	    results[0].terminal) {
	    printf("FAIL: compact projected LoD did not apply source placement exactly once\n");
	    ret = 1;
	}
    }

    if (!ret) {
	SoBRLLodUpdateAction update;
	update.setViewLodState(&viewState);
	update.setResults(results);
	update.apply(root);
	const BObolViewLodState::CadPayload *payload =
	    viewState.findCadForResult(results[0]);
	if (update.getAppliedResultCount() != 1 || !payload ||
	    payload->activeLevel < 0 ||
	    payload->counts.faceCount != results[0].counts.faceCount) {
	    printf("FAIL: compact projected LoD payload was not installed\n");
	    ret = 1;
	}
    }

    /* Finish the coverage-first request before testing stable presentation
     * behavior.  An unchanged demand must jump from the installed minimum
     * directly to the requested cut, without walking nominal PoP levels. */
    if (!ret) {
	SoBRLMeshLodSubmitAction refine;
	refine.setService(&service);
	refine.setDatabase(dbip, "db://compact-projected-test", 2026);
	refine.setViewInfo(&view);
	refine.setViewVolume(&volume, 1.0f);
	refine.setGeneration(service.beginGeneration());
	refine.setRevisions(61, 62);
	refine.setViewLodState(&viewState);
	refine.apply(root);
	std::vector<BObolLodResult> refinedResults;
	const bool retainedCut =
	    refine.getSubmittedTaskCount() == 0 &&
	    refine.getUpdatedCutCount() == 1;
	if ((!retainedCut && refine.getSubmittedTaskCount() != 1) ||
	    (!retainedCut && wait_for_service(service))) {
	    printf("FAIL: compact projected LoD did not finish coverage "
		   "refinement (tasks=%u cuts=%u pending=%u results=%zu "
		   "active=%d requested=%d terminal=%d diagnostics=%s)\n",
		   refine.getSubmittedTaskCount(),
		   refine.getUpdatedCutCount(),
		   refine.getPendingRetainedRefinementCount(),
		   refinedResults.size(),
		   refinedResults.empty() ? -1 :
		       refinedResults[0].geometry.activeLevel,
		   refinedResults.empty() ? -1 :
		       refinedResults[0].request.requestedLevel,
		   refinedResults.empty() ? 0 :
		       (refinedResults[0].terminal ? 1 : 0),
		   refine.getDiagnostics().getString());
	    ret = 1;
	} else if (retainedCut) {
	    const BObolViewLodState::CadPayload *payload =
		viewState.findCadForResult(results[0]);
	    if (service.drainResults(refinedResults) != 0 ||
		!payload ||
		payload->activeLevel != payload->requestedLevel) {
		printf("FAIL: compact projected LoD retained refinement did "
		       "not reach its requested cut\n");
		ret = 1;
	    } else {
		/* Keep the fixture's result descriptor synchronized with the
		 * in-place retained cut used by the remaining presentation
		 * assertions. */
		results[0].geometry.activeLevel = payload->activeLevel;
		results[0].residentLevel = payload->residentLevel;
		results[0].counts = payload->counts;
		results[0].terminal = TRUE;
	    }
	} else {
	    if (service.drainResults(refinedResults) != 1 ||
		refinedResults.size() != 1 ||
		refinedResults[0].geometry.activeLevel !=
		    refinedResults[0].request.requestedLevel ||
		!refinedResults[0].terminal) {
		printf("FAIL: compact projected LoD provider refinement did "
		       "not reach its requested cut\n");
		ret = 1;
	    } else {
		SoBRLLodUpdateAction update;
		update.setViewLodState(&viewState);
		update.setResults(refinedResults);
		update.apply(root);
		if (update.getAppliedResultCount() != 1) {
		printf("FAIL: compact projected LoD refinement was not installed\n");
		ret = 1;
		} else {
		    results.swap(refinedResults);
		}
	    }
	}
    }

    /* Renderer-only motion limits must be reversible without another PoP
     * planning pass.  A deliberately slow measured frame lowers the moving
     * ceiling; ending a pose-only orthographic gesture restores the prior
     * stable ceiling immediately while the refinement debounce remains
     * active. */
    if (!ret) {
	BObolViewController restoreController(root, camera);
	restoreController.setLodAutoSubmit(TRUE);
	BObolViewLodState *restoreState =
	    restoreController.getViewLodState();
	if (!restoreState ||
	    !restoreState->applySourceResult(source, results[0])) {
	    printf("FAIL: pose-only stable presentation fixture apply\n");
	    ret = 1;
	} else {
	    const int stableCeiling =
		restoreController.getLodInteractiveProgressiveCeiling();
	    const uint64_t now = restoreController.beginRenderTiming();
	    restoreController.completeRenderTiming(
		now > 100000000ULL ? now - 100000000ULL : 1);
	    restoreController.beginLodInteraction();
	    const int motionCeiling =
		restoreController.getLodInteractiveProgressiveCeiling();
	    restoreController.endLodInteraction();
	    const int releaseCeiling =
		restoreController.getLodInteractiveProgressiveCeiling();
	    if (stableCeiling != -1 || motionCeiling < 0 ||
		releaseCeiling != stableCeiling ||
		restoreController.isLodGestureActive() ||
		!restoreController.isLodInteractionActive()) {
		printf("FAIL: pose-only release did not immediately restore "
		       "the resident stable presentation (stable=%d motion=%d "
		       "release=%d gesture=%d interactive=%d)\n",
		       stableCeiling, motionCeiling, releaseCeiling,
		       restoreController.isLodGestureActive() ? 1 : 0,
		       restoreController.isLodInteractionActive() ? 1 : 0);
		ret = 1;
	    }
	}
    }

    /* Face count alone is not an occurrence-identity witness.  Replacing a
     * payload during the gesture with an equal-cost payload must invalidate
     * the snapshot and keep the cautious release ceiling until the quiet
     * allocator verifies the new population. */
    if (!ret) {
	BObolViewController changedController(root, camera);
	changedController.setLodAutoSubmit(TRUE);
	BObolViewLodState *changedState =
	    changedController.getViewLodState();
	if (!changedState ||
	    !changedState->applySourceResult(source, results[0])) {
	    printf("FAIL: changed-population presentation fixture apply\n");
	    ret = 1;
	} else {
	    const uint64_t now = changedController.beginRenderTiming();
	    changedController.completeRenderTiming(
		now > 100000000ULL ? now - 100000000ULL : 1);
	    changedController.beginLodInteraction();
	    const int motionCeiling =
		changedController.getLodInteractiveProgressiveCeiling();
	    const BObolViewLodState::CadPayload *changedPayload =
		changedState->findCadForResult(results[0]);
	    if (!changedPayload ||
		!changedState->removeCadPayload(changedPayload) ||
		!changedState->applySourceResult(source, results[0])) {
		printf("FAIL: changed-population presentation mutation\n");
		ret = 1;
	    } else {
		changedController.endLodInteraction();
		const int releaseCeiling =
		    changedController.getLodInteractiveProgressiveCeiling();
		if (motionCeiling < 0 || releaseCeiling != motionCeiling) {
		    printf("FAIL: equal-cost population change incorrectly "
			   "restored a stale presentation (motion=%d release=%d)\n",
			   motionCeiling, releaseCeiling);
		    ret = 1;
		}
	    }
	}
    }

    /* A terminal provider failure must not turn a quiet view into an
     * unbounded submit/result/reject loop.  It covers the occurrence with
     * its structural fallback for this exact demand, while a new view epoch
     * remains eligible to retry and an eventual success clears the record. */
    if (!ret) {
	BObolViewLodState failedState;
	BObolLodResult failedResult;
	failedResult.request = results[0].request;
	failedResult.cacheKey = bobol_lod_cache_key(failedResult.request);
	failedResult.providerStatus = BOBOL_LOD_PROVIDER_CACHE_MISS;
	failedResult.terminal = TRUE;
	failedResult.diagnostic = "unit-test terminal cache failure";
	if (!failedState.applySourceResult(source, failedResult) ||
	    !failedState.hasCadOccurrenceTerminalFailure(
		source, failedResult.request)) {
	    printf("FAIL: compact terminal provider failure was not retained\n");
	    ret = 1;
	}
	if (!ret) {
	    SoBRLMeshLodSubmitAction failedSubmit;
	    failedSubmit.setService(&service);
	    failedSubmit.setDatabase(
		dbip, "db://compact-projected-test", 2026);
	    failedSubmit.setViewInfo(&view);
	    failedSubmit.setViewVolume(&volume, 1.0f);
	    failedSubmit.setGeneration(service.beginGeneration());
	    failedSubmit.setRevisions(61, 62);
	    failedSubmit.setViewLodState(&failedState);
	    failedSubmit.apply(root);
	    if (failedSubmit.getVisitedMeshCount() != 1 ||
		failedSubmit.getVisibleMeshCount() != 1 ||
		failedSubmit.getCoveredVisibleMeshCount() != 1 ||
		failedSubmit.getSubmittedTaskCount() != 0 ||
		failedSubmit.getSkippedMeshCount() != 1) {
		printf("FAIL: compact terminal provider failure was "
		       "resubmitted (visited=%u visible=%zu covered=%zu "
		       "tasks=%u skipped=%u diagnostics=%s)\n",
		       failedSubmit.getVisitedMeshCount(),
		       failedSubmit.getVisibleMeshCount(),
		       failedSubmit.getCoveredVisibleMeshCount(),
		       failedSubmit.getSubmittedTaskCount(),
		       failedSubmit.getSkippedMeshCount(),
		       failedSubmit.getDiagnostics().getString());
		ret = 1;
	    }
	}
	BObolLodRequest nextViewRequest = failedResult.request;
	nextViewRequest.viewRevision++;
	if (!ret && failedState.hasCadOccurrenceTerminalFailure(
		source, nextViewRequest)) {
	    printf("FAIL: compact terminal provider failure suppressed a new "
		   "view epoch\n");
	    ret = 1;
	}
	if (!ret &&
	    (!failedState.applySourceResult(source, results[0]) ||
	     failedState.hasCadOccurrenceTerminalFailure(
		 source, failedResult.request))) {
	    printf("FAIL: compact successful result did not clear its exact "
		   "terminal failure\n");
	    ret = 1;
	}
    }

    /* A retained hierarchy may make a very large population jump between
     * adjacent levels.  The submit action must reserve that growth before
     * changing the cut or scheduling a provider task, while preserving the
     * richer view demand for a later budget epoch. */
    if (!ret) {
	point_t points[6];
	int faces[12] = {
	    0, 1, 2,
	    0, 2, 3,
	    0, 3, 4,
	    0, 4, 5
	};
	for (int i = 0; i < 6; ++i)
	    VSET(points[i], (i == 1 || i == 2) ? 1.0 : 0.0,
		(i >= 2 && i <= 4) ? 1.0 : 0.0,
		i == 5 ? 1.0 : 0.0);
	struct BObolMeshLodData data = {};
	data.faces = faces;
	data.face_count = 4;
	data.points = points;
	data.point_count = 6;
	data.points_orig = points;
	data.point_orig_count = 6;
	VSET(data.bmin, 0.0, 0.0, 0.0);
	VSET(data.bmax, 1.0, 1.0, 1.0);
	struct BObolMeshLodHierarchyInfo hierarchy =
	    BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
	hierarchy.min_level = 0;
	hierarchy.max_level = 1;
	hierarchy.resident_level = 1;
	hierarchy.face_count[0] = 1;
	hierarchy.point_count[0] = 3;
	hierarchy.face_count[1] = 4;
	hierarchy.point_count[1] = 6;
	VSET(hierarchy.quantization_min, 0.0, 0.0, 0.0);
	VSET(hierarchy.quantization_max, 1.0, 1.0, 1.0);

	BObolLodProgressiveMeshPtr progressive(
	    new BObolLodProgressiveMesh);
	BObolLodProgressiveMeshPtr partialProgressive(
	    new BObolLodProgressiveMesh);
	struct BObolMeshLodData partialData = data;
	partialData.face_count = 1;
	partialData.point_count = 3;
	partialData.point_orig_count = 3;
	struct BObolMeshLodHierarchyInfo partialHierarchy = hierarchy;
	partialHierarchy.resident_level = 0;
	BObolLodResult budgetResult = results[0];
	if (!progressive->update(data, hierarchy, 1, FALSE) ||
	    !partialProgressive->update(
		partialData, partialHierarchy, 0, FALSE) ||
	    partialProgressive->faceCount(1) != 1 ||
	    partialProgressive->hierarchyFaceCount(1) != 4) {
	    printf("FAIL: retained refinement budget fixture setup\n");
	    ret = 1;
	} else {
	    budgetResult.progressiveMesh = partialProgressive;
	    budgetResult.mesh.clear();
	    budgetResult.geometry.activeLevel = 0;
	    budgetResult.residentLevel = 0;
	    budgetResult.request.requestedLevel = 1;
	    budgetResult.counts.faceCount = 1;
	    budgetResult.counts.pointCount = 3;
	    budgetResult.counts.originalPointCount = 3;
	    budgetResult.counts.normalCount = 0;
	    budgetResult.terminal = FALSE;

	    BObolViewLodState budgetState;
	    SoBRLLodUpdateAction budgetUpdate;
	    budgetUpdate.setViewLodState(&budgetState);
	    budgetUpdate.addResult(budgetResult);
	    budgetUpdate.apply(root);
	    const BObolViewLodState::CadPayload *payload =
		budgetState.findCadForResult(budgetResult);
	    std::vector<BObolCompactOccurrence> rekeyedOccurrences =
		occurrences;
	    rekeyedOccurrences[0].summary.path =
		"root.c/post-pca/lod-submit.bot";
	    rekeyedOccurrences[0].sourceMeshRequest.meshAssetPath =
		budgetResult.request.objectPath;
	    rekeyedOccurrences[0].sourceMeshRequest.meshAssetName =
		budgetResult.request.objectName;
	    if (source->setCompactOccurrenceRegistry(rekeyedOccurrences) < 0) {
		printf("FAIL: retained refinement rekey fixture setup\n");
		ret = 1;
	    }

	    /* A cumulative resident suffix remains useful when a newer wheel
	     * event advances the camera epoch before its worker result is drained.
	     * Rebase that result to the current occurrence demand instead of
	     * rejecting it, but preserve the provider's independently admitted
	     * presentation cut: residency is not permission to draw the whole
	     * newly loaded prefix. */
	    if (!ret) {
		BObolLodProgressiveMeshPtr rebaseProgressive(
		    new BObolLodProgressiveMesh);
		BObolLodResult rebaseSeed = budgetResult;
		if (!rebaseProgressive->update(
			partialData, partialHierarchy, 0, FALSE)) {
		    printf("FAIL: cumulative suffix rebase fixture setup\n");
		    ret = 1;
		} else {
		    rebaseSeed.progressiveMesh = rebaseProgressive;
		    rebaseSeed.mesh.clear();
		    BObolCompactInstanceHandle rebaseHandle;
		    BObolCompactInstanceSummary rebaseSummary;
		    if (!source->getCompactInstanceHandle(0, rebaseHandle) ||
			!source->getCompactInstanceSummary(
			    rebaseHandle, rebaseSummary)) {
			printf("FAIL: cumulative suffix rebase occurrence setup\n");
			ret = 1;
		    }
		    rebaseSeed.request.occurrenceKey =
			rebaseSummary.sourceInstanceKey;
		    rebaseSeed.request.objectPath =
			rebaseSummary.meshAssetPath;
		    rebaseSeed.request.objectName =
			rebaseSummary.meshAssetName;
		    rebaseSeed.geometry.activeLevel = 0;
		    rebaseSeed.residentLevel = 0;
		    rebaseSeed.request.requestedLevel = 1;
		    rebaseSeed.counts.faceCount = 1;
		    rebaseSeed.counts.pointCount = 3;
		    rebaseSeed.counts.originalPointCount = 3;
		    rebaseSeed.counts.normalCount = 0;
		    rebaseSeed.terminal = FALSE;

		    BObolViewController rebaseController(root, camera);
		    rebaseController.setLodService(&service);
		    rebaseController.setLodAutoSubmit(TRUE);
		    BObolViewLodState *rebaseState =
			rebaseController.getViewLodState();
		    rebaseSeed.request.viewRevision =
			rebaseController.getLodViewRevision();
		    rebaseSeed.request.policyRevision =
			rebaseController.getLodPolicyRevision();
		    rebaseSeed.request.sourceRoutingId =
			source->getCompactSourceRoutingId();
		    rebaseSeed.cacheKey =
			bobol_lod_cache_key(rebaseSeed.request);
		    const uint64_t staleView =
			rebaseSeed.request.viewRevision.value();
		    if (!rebaseState ||
			!rebaseState->applySourceResult(source, rebaseSeed)) {
			printf("FAIL: cumulative suffix rebase seed apply\n");
			ret = 1;
		    } else if (!rebaseProgressive->update(
			    data, hierarchy, 1, FALSE)) {
			printf("FAIL: cumulative suffix rebase growth\n");
			ret = 1;
		    } else {
			BObolLodResult completed = rebaseSeed;
			completed.residentLevel = 1;
			completed.geometry.activeLevel = 0;
			completed.counts.faceCount = 1;
			completed.counts.pointCount = 3;
			completed.counts.originalPointCount = 3;
			completed.terminal = FALSE;
			completed.preparedCadGeometry =
			    rebaseProgressive->prepareCadGeometry(
				BOBOL_LOD_DRAW_SHADED,
				&completed.preparedCadGeometryRevision);
			const uint64_t generation =
			    rebaseController.beginLodGeneration();
			BObolLodTask task;
			task.generation = generation;
			task.request = completed.request;
			task.realize = copied_progressive_result_task;
			task.realizeData = &completed;
			task.debugDelayMilliseconds = 50;
			if (!generation || service.submit(task) == 0) {
			    printf("FAIL: cumulative suffix rebase task submit\n");
			    ret = 1;
			} else {
			    rebaseController.setViewportSize(641, 480);
			    const BObolViewLodState::CadPayload *current =
				rebaseState->findCadForResult(rebaseSeed);
			    /* Deliberately leave the retained occurrence stamped with
			     * staleView.  Continuous wheel events can outrun both the
			     * worker result and occurrence metadata; immutable asset
			     * identity, not a current camera stamp, is the safe rebase
			     * witness. */
			    const int waited = wait_for_service(service);
			    rebaseController.clearProgressiveWorkPending();
			    const int applied = waited ? -1 :
				rebaseController.applyLodResults(
				    &service, 0, 0, generation);
			    if (!current || waited || applied != 1) {
				printf("FAIL: cumulative suffix current-view rebase "
				       "application (current=%p "
				       "wait=%d applied=%d processed=%zu matched=%u "
				       "rejected=%u unmatched=%u routing=%llu/%llu "
				       "diagnostics=%s)\n",
				       (const void *)current, waited, applied,
				       rebaseController.getLastLodResultCount(),
				       rebaseController.getLastLodMatchedResultCount(),
				       rebaseController.getLastLodRejectedResultCount(),
				       rebaseController.getLastLodUnmatchedResultCount(),
				       (unsigned long long)
					   completed.request.sourceRoutingId.value(),
				       (unsigned long long)
					   source->getCompactSourceRoutingId(),
				       rebaseController.getLastLodDiagnostics().getString());
				ret = 1;
			    } else {
				const BObolViewLodState::CadPayload *after =
				    rebaseState->findCadForResult(rebaseSeed);
				if (!after || after->activeLevel != 0 ||
				    after->residentLevel != 1 ||
				    after->requestedLevel != 1 ||
				    after->viewRevision !=
					rebaseController.getLodViewRevision() ||
				    after->viewRevision == staleView ||
				    after->counts.faceCount != 1 ||
				    !rebaseController.hasProgressiveWorkPending()) {
				    printf("FAIL: cumulative suffix was discarded "
					   "after a newer view (active=%d "
					   "resident=%d requested=%d "
					   "view=%llu/%llu stale=%llu "
					   "faces=%zu pending=%d diagnostics=%s)\n",
					   after ? after->activeLevel : -1,
					   after ? after->residentLevel : -1,
					   after ? after->requestedLevel : -1,
					   (unsigned long long)(after ?
					       after->viewRevision : 0),
					   (unsigned long long)
					       rebaseController.
						   getLodViewRevision(),
					   (unsigned long long)staleView,
					   after ? after->counts.faceCount : 0,
					   rebaseController.
					       hasProgressiveWorkPending() ? 1 : 0,
					   rebaseController.
					       getLastLodDiagnostics().getString());
				    ret = 1;
				}
			    }
			}
		    }
		}
	    }

	    SoOrthographicCamera *budgetCamera =
		new SoOrthographicCamera;
	    budgetCamera->ref();
	    budgetCamera->position = SbVec3f(100.5f, 0.5f, 5.0f);
	    budgetCamera->height = 300.0f;
	    const SbViewVolume budgetVolume =
		budgetCamera->getViewVolume(640.0f / 480.0f);

	    SoBRLMeshLodSubmitAction blocked;
	    blocked.setService(&service);
	    blocked.setDatabase(dbip, "db://compact-projected-test", 2026);
	    blocked.setViewInfo(&view);
	    blocked.setViewVolume(&budgetVolume, 1.0f);
	    blocked.setGeneration(service.beginGeneration());
	    blocked.setRevisions(63, 64);
	    blocked.setViewLodState(&budgetState);
	    blocked.setRefinementCostBudget(3);
	    blocked.apply(root);
	    payload = budgetState.findCadForResult(budgetResult);
	    BObolCompactInstanceHandle budgetHandle;
	    BObolCompactInstanceSummary budgetSummary;
	    (void)source->getCompactInstanceHandle(0, budgetHandle);
	    (void)source->getCompactInstanceSummary(
		budgetHandle, budgetSummary);
	    if (budgetUpdate.getAppliedResultCount() != 1 || !payload ||
		blocked.getSubmittedTaskCount() != 0 ||
		blocked.getUpdatedCutCount() != 0 ||
		blocked.getPendingRetainedRefinementCount() != 1 ||
		blocked.getRefinementBudgetBlockedCount() != 1 ||
		blocked.getRefinementCostBudgetUsed() != 0 ||
		payload->activeLevel != 0 || payload->requestedLevel != 1 ||
		budgetState.activeFaceCount() != 1) {
		printf("FAIL: retained refinement exceeded its aggregate face "
		       "budget (tasks=%u cuts=%u pending=%u blocked=%u used=%zu "
		       "active=%d requested=%d faces=%zu payload_path=%s "
		       "asset_path=%s diagnostics=%s)\n",
		       blocked.getSubmittedTaskCount(),
		       blocked.getUpdatedCutCount(),
		       blocked.getPendingRetainedRefinementCount(),
		       blocked.getRefinementBudgetBlockedCount(),
		       blocked.getRefinementCostBudgetUsed(),
		       payload ? payload->activeLevel : -1,
		       payload ? payload->requestedLevel : -1,
		       budgetState.activeFaceCount(),
		       payload ? payload->sourcePath.getString() : "",
		       budgetSummary.meshAssetPath.getString(),
		       blocked.getDiagnostics().getString());
		ret = 1;
	    }

	    if (!ret) {
		/* A hard resident-memory denial and a scene-cost denial are
		 * independent terminal witnesses.  Rechecking the same demand
		 * must not erase the resident admission epoch merely because the
		 * scene allocator also cannot afford its next transition. */
		BObolLodResult limitedResult = budgetResult;
		limitedResult.request.viewRevision = 81;
		limitedResult.request.policyRevision = 82;
		limitedResult.memoryLimited = TRUE;
		const uint64_t denialRevision =
		    service.residentMeshAdmissionRevision() + 1;
		limitedResult.residentAdmissionRevision = denialRevision;
		BObolViewLodState limitedState;
		if (!limitedState.applySourceResult(source, limitedResult)) {
		    printf("FAIL: memory-limited frontier fixture apply\n");
		    ret = 1;
		} else {
		    std::vector<SbString> blockedFrontier;
		    limitedState.unsatisfiedCadOccurrenceKeys(
			source, denialRevision, blockedFrontier);
		    std::vector<SbString> reopenedFrontier;
		    limitedState.unsatisfiedCadOccurrenceKeys(
			source, service.residentMeshAdmissionRevision(),
			reopenedFrontier);
		    SoBRLMeshLodSubmitAction repeatedDenial;
		    repeatedDenial.setService(&service);
		    repeatedDenial.setDatabase(
			dbip, "db://compact-projected-test", 2026);
		    repeatedDenial.setViewInfo(&view);
		    repeatedDenial.setViewVolume(&budgetVolume, 1.0f);
		    repeatedDenial.setGeneration(service.beginGeneration());
		    repeatedDenial.setRevisions(81, 82);
		    repeatedDenial.setViewLodState(&limitedState);
		    repeatedDenial.setRefinementCostBudget(3);
		    repeatedDenial.apply(root);
		    const BObolViewLodState::CadPayload *limitedPayload =
			limitedState.findCadForResult(limitedResult);
		    if (!blockedFrontier.empty() ||
			reopenedFrontier.size() != 1 || !limitedPayload ||
			repeatedDenial.getSubmittedTaskCount() != 0 ||
			repeatedDenial.getPendingRetainedRefinementCount() != 1 ||
			repeatedDenial.getRefinementBudgetBlockedCount() != 1 ||
			!limitedPayload->memoryLimited ||
			limitedPayload->residentAdmissionRevision !=
			    denialRevision) {
			printf("FAIL: repeated budget denial lost or rescanned its "
			       "memory-capacity witness (blocked=%zu reopened=%zu "
			       "tasks=%u pending=%u budget_blocked=%u "
			       "limited=%d admission=%llu/%llu)\n",
			       blockedFrontier.size(), reopenedFrontier.size(),
			       repeatedDenial.getSubmittedTaskCount(),
			       repeatedDenial.
				   getPendingRetainedRefinementCount(),
			       repeatedDenial.getRefinementBudgetBlockedCount(),
			       limitedPayload && limitedPayload->memoryLimited ?
				   1 : 0,
			       (unsigned long long)(limitedPayload ?
				   limitedPayload->residentAdmissionRevision : 0),
			       (unsigned long long)denialRevision);
			ret = 1;
		    }
		}
	    }

	    if (!ret) {
		/* A scale-quality experiment is scene-global even though one
		 * submit action may visit several occurrences.  Let the first
		 * complete marginal transition use one unit beyond the ordinary
		 * budget, then prove the second occurrence remains pending. */
		std::vector<BObolCompactOccurrence> trialOccurrences =
		    rekeyedOccurrences;
		BObolCompactOccurrence secondTrial = trialOccurrences.front();
		secondTrial.summary.path =
		    "root.c/post-pca-second/lod-submit.bot";
		trialOccurrences.push_back(secondTrial);
		if (source->setCompactOccurrenceRegistry(trialOccurrences) < 0) {
		    printf("FAIL: one-over-budget trial registry setup\n");
		    ret = 1;
		} else {
		    BObolViewLodState trialState;
		    std::vector<BObolLodResult> trialResults;
		    for (size_t trialIndex = 0;
			 trialIndex < trialOccurrences.size(); ++trialIndex) {
			BObolCompactInstanceHandle trialHandle;
			BObolCompactInstanceSummary trialSummary;
			if (!source->getCompactInstanceHandle(
				static_cast<int>(trialIndex), trialHandle) ||
			    !source->getCompactInstanceSummary(
				trialHandle, trialSummary)) {
			    ret = 1;
			    break;
			}
			BObolLodResult trialResult = budgetResult;
			trialResult.progressiveMesh = progressive;
			trialResult.geometry.activeLevel = 0;
			trialResult.residentLevel = 1;
			trialResult.request.requestedLevel = 1;
			trialResult.request.occurrenceKey =
			    trialSummary.sourceInstanceKey;
			trialResult.request.objectPath =
			    trialSummary.meshAssetPath;
			trialResult.request.objectName =
			    trialSummary.meshAssetName;
			if (!trialState.applySourceResult(source, trialResult)) {
			    ret = 1;
			    break;
			}
			trialResults.push_back(trialResult);
		    }
		    if (ret) {
			printf("FAIL: one-over-budget trial payload setup\n");
		    } else {
			SoBRLMeshLodSubmitAction trial;
			trial.setService(&service);
			trial.setDatabase(
			    dbip, "db://compact-projected-test", 2026);
			trial.setViewInfo(&view);
			trial.setViewVolume(&budgetVolume, 1.0f);
			trial.setGeneration(service.beginGeneration());
			trial.setRevisions(63, 65);
			trial.setViewLodState(&trialState);
			trial.setRefinementCostBudget(3);
			trial.setOneOverBudgetRefinementLimit(1);
			trial.apply(root);
			unsigned int coarseCount = 0;
			unsigned int richCount = 0;
			for (const BObolLodResult &trialResult : trialResults) {
			    const BObolViewLodState::CadPayload *trialPayload =
				trialState.findCadForResult(trialResult);
			    if (trialPayload && trialPayload->activeLevel == 0)
				coarseCount++;
			    if (trialPayload && trialPayload->activeLevel == 1)
				richCount++;
			}
			if (trial.getOneOverBudgetRefinementLimit() != 1 ||
			    !trial.getOneOverBudgetRefinementUsed() ||
			    trial.getUpdatedCutCount() != 1 ||
			    trial.getRefinementCostBudgetUsed() != 4 ||
			    trial.getRefinementBudgetBlockedCount() < 1 ||
			    trial.getPendingRetainedRefinementCount() != 1 ||
			    coarseCount != 1 || richCount != 1 ||
			    trialState.activeFaceCount() != 5) {
			    printf("FAIL: one-over-budget trial was not unique and "
				   "exactly charged (limit=%zu used=%d cuts=%u "
				   "cost=%zu blocked=%u pending=%u coarse=%u "
				   "rich=%u faces=%zu diagnostics=%s)\n",
				   trial.getOneOverBudgetRefinementLimit(),
				   trial.getOneOverBudgetRefinementUsed() ? 1 : 0,
				   trial.getUpdatedCutCount(),
				   trial.getRefinementCostBudgetUsed(),
				   trial.getRefinementBudgetBlockedCount(),
				   trial.getPendingRetainedRefinementCount(),
				   coarseCount, richCount,
				   trialState.activeFaceCount(),
				   trial.getDiagnostics().getString());
			    ret = 1;
			}
		    }
		}
		if (source->setCompactOccurrenceRegistry(rekeyedOccurrences) < 0) {
		    printf("FAIL: one-over-budget trial registry restore\n");
		    ret = 1;
		}
	    }

	    if (!ret) {
		BObolLodResult admittedResult = budgetResult;
		admittedResult.progressiveMesh = progressive;
		admittedResult.residentLevel = 1;
		BObolViewLodState admittedState;
		if (!admittedState.applySourceResult(source, admittedResult)) {
		    printf("FAIL: affordable retained refinement fixture apply\n");
		    ret = 1;
		}
		BObolViewLodState ceilingState;
		if (!ret &&
		    !ceilingState.applySourceResult(source, admittedResult)) {
		    printf("FAIL: presentation ceiling fixture apply\n");
		    ret = 1;
		}
		if (!ret) {
		    SoBRLMeshLodSubmitAction ceilingAction;
		    ceilingAction.setService(&service);
		    ceilingAction.setDatabase(
			dbip, "db://compact-projected-test", 2026);
		    ceilingAction.setViewInfo(&view);
		    ceilingAction.setViewVolume(&budgetVolume, 1.0f);
		    ceilingAction.setGeneration(service.beginGeneration());
		    ceilingAction.setRevisions(65, 66);
		    ceilingAction.setViewLodState(&ceilingState);
		    ceilingAction.setRefinementCostBudget(4);
		    ceilingAction.setRefinementLevelCeiling(0);
		    ceilingAction.apply(root);
		    const BObolViewLodState::CadPayload *ceilingPayload =
			ceilingState.findCadForResult(admittedResult);
		    if (!ceilingPayload ||
			ceilingAction.getRefinementLevelCeiling() != 0 ||
			ceilingAction.getUpdatedCutCount() != 0 ||
			ceilingAction.getRefinementCostBudgetUsed() != 0 ||
			ceilingPayload->activeLevel != 0 ||
			ceilingPayload->residentLevel != 1) {
			printf("FAIL: renderer ceiling allowed hidden occurrence "
			       "promotion (ceiling=%d cuts=%u used=%zu active=%d "
			       "resident=%d)\n",
			       ceilingAction.getRefinementLevelCeiling(),
			       ceilingAction.getUpdatedCutCount(),
			       ceilingAction.getRefinementCostBudgetUsed(),
			       ceilingPayload ? ceilingPayload->activeLevel : -1,
			       ceilingPayload ? ceilingPayload->residentLevel : -1);
			ret = 1;
		    }
		}
		SoBRLMeshLodSubmitAction admitted;
		admitted.setService(&service);
		admitted.setDatabase(dbip, "db://compact-projected-test",
		    2026);
		admitted.setViewInfo(&view);
		admitted.setViewVolume(&budgetVolume, 1.0f);
		admitted.setGeneration(service.beginGeneration());
		admitted.setRevisions(65, 66);
		admitted.setViewLodState(&admittedState);
		admitted.setRefinementCostBudget(4);
		admitted.apply(root);
		payload = admittedState.findCadForResult(admittedResult);
		if (!payload || admitted.getSubmittedTaskCount() != 0 ||
		    admitted.getUpdatedCutCount() != 1 ||
		    admitted.getPendingRetainedRefinementCount() != 0 ||
		    admitted.getRefinementBudgetBlockedCount() != 0 ||
		    admitted.getRefinementCostBudgetUsed() != 4 ||
		    payload->activeLevel != 1 ||
		    admittedState.activeFaceCount() != 4) {
		    printf("FAIL: retained refinement did not use an affordable "
			   "render-cost budget (tasks=%u cuts=%u pending=%u blocked=%u "
			   "used=%zu active=%d faces=%zu diagnostics=%s)\n",
			   admitted.getSubmittedTaskCount(),
			   admitted.getUpdatedCutCount(),
			   admitted.getPendingRetainedRefinementCount(),
			   admitted.getRefinementBudgetBlockedCount(),
			   admitted.getRefinementCostBudgetUsed(),
			   payload ? payload->activeLevel : -1,
			   admittedState.activeFaceCount(),
			   admitted.getDiagnostics().getString());
		    ret = 1;
		}
		if (!ret) {
		    /*
		     * A bounded coverage-recovery window may temporarily draw
		     * the minimum prefix, but that cut is not the view target.
		     * Keep the richer request on the sparse unsatisfied frontier
		     * so a following pass cannot falsely declare convergence.
		     */
		    SoBRLMeshLodSubmitAction recovery;
		    recovery.setService(&service);
		    recovery.setDatabase(
			dbip, "db://compact-projected-test", 2026);
		    recovery.setViewInfo(&view);
		    recovery.setViewVolume(&budgetVolume, 1.0f);
		    recovery.setGeneration(service.beginGeneration());
		    recovery.setRevisions(67, 68);
		    recovery.setViewLodState(&admittedState);
		    recovery.setAllowLevelDowngrade(TRUE);
		    recovery.setPreserveMeshCoverage(TRUE);
		    recovery.setRefinementCostBudget(0);
		    recovery.setRetainedSceneCostBudget(66);
		    recovery.setSubmissionTaskLimit(0);
		    recovery.setCompactEntryRange(0, 1);
		    recovery.apply(root);
		    payload = admittedState.findCadForResult(admittedResult);
		    std::vector<SbString> unsatisfied;
		    admittedState.unsatisfiedCadOccurrenceKeys(
			source, 0, unsatisfied);
		    if (!payload ||
			recovery.getSubmittedTaskCount() != 0 ||
			recovery.getUpdatedCutCount() != 1 ||
			recovery.getPendingRetainedRefinementCount() != 1 ||
			recovery.getRefinementBudgetBlockedCount() != 0 ||
			recovery.getRetainedSceneCostBudgetUsed() != 66 ||
			payload->activeLevel != 0 ||
			payload->requestedLevel != 1 ||
			unsatisfied.size() != 1 ||
			bu_strcmp(unsatisfied[0].getString(),
			    payload->sourceInstanceKey.getString()) != 0) {
			printf("FAIL: bounded retained recovery replaced the "
			       "view target with its temporary minimum cut "
			       "(tasks=%u cuts=%u pending=%u blocked=%u "
			       "used=%zu active=%d requested=%d "
			       "unsatisfied=%zu)\n",
			       recovery.getSubmittedTaskCount(),
			       recovery.getUpdatedCutCount(),
			       recovery.
				   getPendingRetainedRefinementCount(),
			       recovery.getRefinementBudgetBlockedCount(),
			       recovery.getRetainedSceneCostBudgetUsed(),
			       payload ? payload->activeLevel : -1,
			       payload ? payload->requestedLevel : -1,
			       unsatisfied.size());
			ret = 1;
		    }
		}
	    }
	    if (!ret) {
		/* A shared asset may retain several levels above this
		 * occurrence's active cut while the view requests a still-richer
		 * nonresident level.  Advance to the next affordable retained
		 * level instead of repeatedly trying (and failing) to charge the
		 * whole active-to-resident jump. */
		struct BObolMeshLodHierarchyInfo multiHierarchy = hierarchy;
		multiHierarchy.max_level = 3;
		multiHierarchy.resident_level = 2;
		multiHierarchy.face_count[0] = 1;
		multiHierarchy.face_count[1] = 2;
		multiHierarchy.face_count[2] = 4;
		multiHierarchy.face_count[3] = 8;
		multiHierarchy.point_count[0] = 3;
		multiHierarchy.point_count[1] = 4;
		multiHierarchy.point_count[2] = 6;
		multiHierarchy.point_count[3] = 8;
		BObolLodProgressiveMeshPtr multiProgressive(
		    new BObolLodProgressiveMesh);
		BObolLodResult multiResult = budgetResult;
		if (!multiProgressive->update(
			data, multiHierarchy, 2, FALSE)) {
		    printf("FAIL: multi-level retained refinement fixture setup\n");
		    ret = 1;
		} else {
		    multiResult.progressiveMesh = multiProgressive;
		    multiResult.geometry.activeLevel = 0;
		    multiResult.residentLevel = 2;
		    multiResult.request.requestedLevel = 3;
		    multiResult.counts.faceCount = 1;
		    multiResult.counts.pointCount = 3;
		    multiResult.counts.originalPointCount = 3;
		    multiResult.terminal = FALSE;
		    /* A perceptually ordered finite-budget visit is scored for
		     * one marginal population transition.  Even when the full
		     * level-2 cut fits, admitting it in one leap would invalidate
		     * global minimax ordering and allow one prominent occurrence
		     * to consume the allowance intended for its peers. */
		    if (!ret) {
			BObolLodResult fairResult = multiResult;
			fairResult.request.requestedLevel = 2;
			BObolViewLodState fairState;
			if (!fairState.applySourceResult(source, fairResult)) {
			    printf("FAIL: transition-limited refinement fixture apply\n");
			    ret = 1;
			} else {
			    SoBRLMeshLodSubmitAction fair;
			    fair.setService(&service);
			    fair.setDatabase(dbip,
				"db://compact-projected-test", 2026);
			    fair.setViewInfo(&view);
			    fair.setViewVolume(&budgetVolume, 0.5f);
			    fair.setGeneration(service.beginGeneration());
			    fair.setRevisions(73, 74);
			    fair.setViewLodState(&fairState);
			    fair.setRefinementCostBudget(3);
			    fair.setTransitionLimitedRefinement(TRUE);
			    fair.apply(root);
			    payload = fairState.findCadForResult(fairResult);
			    if (!payload || fair.getSubmittedTaskCount() != 0 ||
				fair.getUpdatedCutCount() != 1 ||
				fair.getPendingRetainedRefinementCount() != 1 ||
				fair.getRefinementCostBudgetUsed() != 1 ||
				payload->activeLevel != 1 ||
				payload->requestedLevel != 2 ||
				fairState.activeFaceCount() != 2) {
				printf("FAIL: perceptual refinement consumed more "
				       "than one populated transition (tasks=%u "
				       "cuts=%u pending=%u used=%zu active=%d "
				       "requested=%d faces=%zu diagnostics=%s)\n",
				       fair.getSubmittedTaskCount(),
				       fair.getUpdatedCutCount(),
				       fair.getPendingRetainedRefinementCount(),
				       fair.getRefinementCostBudgetUsed(),
				       payload ? payload->activeLevel : -1,
				       payload ? payload->requestedLevel : -1,
				       fairState.activeFaceCount(),
				       fair.getDiagnostics().getString());
				ret = 1;
			    }
			}
		    }
		    BObolViewLodState multiState;
		    if (!multiState.applySourceResult(source, multiResult)) {
			printf("FAIL: multi-level retained refinement apply\n");
			ret = 1;
		    } else {
			std::vector<SbString> unsatisfied;
			multiState.unsatisfiedCadOccurrenceKeys(
			    source, 0, unsatisfied);
			if (unsatisfied.size() != 1 ||
			    bu_strcmp(unsatisfied[0].getString(),
				multiResult.request.occurrenceKey.getString()) !=
				0) {
			    printf("FAIL: retained refinement frontier did not "
				   "retain its unsatisfied occurrence\n");
			    ret = 1;
			}
			SoBRLMeshLodSubmitAction incremental;
			incremental.setService(&service);
			incremental.setDatabase(dbip,
			    "db://compact-projected-test", 2026);
			incremental.setViewInfo(&view);
			incremental.setViewVolume(
			    &budgetVolume, 0.25f);
			incremental.setGeneration(
			    service.beginGeneration());
			incremental.setRevisions(69, 70);
			incremental.setViewLodState(&multiState);
		incremental.setRefinementCostBudget(1);
			incremental.apply(root);
			payload = multiState.findCadForResult(multiResult);
			if (!payload ||
			    incremental.getSubmittedTaskCount() != 0 ||
			    incremental.getUpdatedCutCount() != 1 ||
			    incremental.getPendingRetainedRefinementCount() != 1 ||
			    incremental.getRefinementBudgetBlockedCount() < 1 ||
		    incremental.getRefinementCostBudgetUsed() != 1 ||
			    payload->activeLevel != 1 ||
			    multiState.activeFaceCount() != 2) {
			    printf("FAIL: retained refinement did not take an "
				   "affordable intermediate cut (tasks=%u cuts=%u "
				   "pending=%u blocked=%u used=%zu active=%d "
				   "faces=%zu diagnostics=%s)\n",
				   incremental.getSubmittedTaskCount(),
				   incremental.getUpdatedCutCount(),
				   incremental.
				       getPendingRetainedRefinementCount(),
				   incremental.
				       getRefinementBudgetBlockedCount(),
			   incremental.getRefinementCostBudgetUsed(),
				   payload ? payload->activeLevel : -1,
				   multiState.activeFaceCount(),
				   incremental.getDiagnostics().getString());
			    ret = 1;
			}
			if (!ret) {
			    unsatisfied.clear();
			    multiState.unsatisfiedCadOccurrenceKeys(
				source, 0, unsatisfied);
			    if (unsatisfied.size() != 1 ||
				!multiState.retargetCadPayload(
				    payload, payload->activeLevel,
				    payload->activeLevel, 71, 72)) {
				printf("FAIL: retained refinement frontier "
				       "lost pending work before satisfaction\n");
				ret = 1;
			    } else {
				multiState.unsatisfiedCadOccurrenceKeys(
				    source, 0, unsatisfied);
				if (!unsatisfied.empty()) {
				    printf("FAIL: retained refinement frontier "
					   "kept a satisfied occurrence\n");
				    ret = 1;
				}
			    }
			}
		    }
		}
	    }
	    budgetCamera->unref();
	    if (source->setCompactOccurrenceRegistry(occurrences) < 0) {
		printf("FAIL: retained refinement fixture registry restore\n");
		ret = 1;
	    }
	}
    }

    /* Shaded and wireframe are two presentations of the same resident PoP
     * cut.  Switching modes must neither load another cut nor fall back to
     * the original box. */
    if (!ret) {
	std::vector<const BObolViewLodState::CadPayload *> payloads;
	viewState.findCadPayloads(source, payloads);
	(void)source->compactViewLodAssembly(payloads, &viewState);
	SoCADAssembly *assembly = viewState.findCadPresentation(source);
	std::vector<Obol::InstanceId> ids =
	    assembly ? assembly->instanceIds() :
	    std::vector<Obol::InstanceId>();
	std::optional<Obol::InstanceRecord> record =
	    ids.size() == 1 ? assembly->getInstanceRecord(ids[0]) :
	    std::optional<Obol::InstanceRecord>();
	const Obol::PartGeometry *geometry =
	    record ? assembly->partGeometry(record->part) : NULL;
	if (!assembly || ids.size() != 1 || !geometry ||
	    !results[0].preparedCadGeometry ||
	    geometry != results[0].preparedCadGeometry.get() ||
	    !geometry->shaded || geometry->wire ||
	    !geometry->shaded->isProgressive() || !record ||
	    !record->localToRoot.equals(placement, 0.000001f) ||
	    record->lodLevel != results[0].geometry.activeLevel ||
	    geometry->shaded->indexCountAtLevel(record->lodLevel) !=
		results[0].counts.faceCount * 3) {
	    printf("FAIL: shaded compact presentation did not retain the "
		   "worker-owned resident PoP geometry and asset transform "
		   "(assembly=%p instances=%zu geometry=%p prepared=%p "
		   "revision=%llu/%llu shaded=%d wire=%d progressive=%d "
		   "record=%d level=%d/%d prepared_range=%d..%d "
		   "indices=%zu/%zu/%llu)\n",
		   (void *)assembly, ids.size(), (const void *)geometry,
		   (const void *)results[0].preparedCadGeometry.get(),
		   (unsigned long long)results[0].preparedCadGeometryRevision,
		   (unsigned long long)(results[0].progressiveMesh ?
		       results[0].progressiveMesh->revision() : 0),
		   geometry && geometry->shaded ? 1 : 0,
		   geometry && geometry->wire ? 1 : 0,
		   geometry && geometry->shaded &&
		       geometry->shaded->isProgressive() ? 1 : 0,
		   record ? 1 : 0, record ? record->lodLevel : -1,
		   results[0].geometry.activeLevel,
		   results[0].preparedCadGeometry &&
		       results[0].preparedCadGeometry->shaded ?
		       results[0].preparedCadGeometry->shaded->
			   progressiveMinimumLevel : -1,
		   results[0].preparedCadGeometry &&
		       results[0].preparedCadGeometry->shaded ?
		       results[0].preparedCadGeometry->shaded->
			   progressiveResidentLevel : -1,
		   geometry && geometry->shaded && record ?
		       geometry->shaded->indexCountAtLevel(record->lodLevel) : 0,
		   results[0].preparedCadGeometry &&
		       results[0].preparedCadGeometry->shaded ?
		       results[0].preparedCadGeometry->shaded->
			   indexCountAtLevel(results[0].geometry.activeLevel) : 0,
		   (unsigned long long)results[0].counts.faceCount * 3);
	    ret = 1;
	}
    }

    /* Compact visibility and selection are presentation revisions, not PoP
     * payload revisions.  A view-specific LoD assembly must update those
     * instance sets in place even though its resident geometry and view-state
     * payload revision are unchanged. */
    if (!ret) {
	std::vector<const BObolViewLodState::CadPayload *> payloads;
	viewState.findCadPayloads(source, payloads);
	const BObolViewLodState::CadPayload *payloadBefore =
	    viewState.findCadForResult(results[0]);
	const size_t payloadCount = viewState.cadPayloadCount();
	(void)source->compactViewLodAssembly(payloads, &viewState);
	SoCADAssembly *assembly = viewState.findCadPresentation(source);
	const std::vector<Obol::InstanceId> ids =
	    assembly ? assembly->instanceIds() :
	    std::vector<Obol::InstanceId>();
	std::vector<SbString> paths = {
	    SbString("root.c/lod-submit.bot")
	};
	std::vector<SbBool> hidden = {FALSE};
	const int visibilityChanged =
	    source->setCompactInstanceVisibilityOverrides(paths, hidden);
	if (!assembly || ids.size() != 1 || visibilityChanged <= 0 ||
	    source->currentCompactViewLodAssembly(&viewState) ||
	    !source->compactViewLodAssembly(payloads, &viewState) ||
	    viewState.findCadPresentation(source) != assembly ||
	    !assembly->isInstanceHidden(ids[0]) ||
	    viewState.findCadForResult(results[0]) != payloadBefore ||
	    viewState.cadPayloadCount() != payloadCount) {
	    printf("FAIL: compact retained visibility did not update the "
		   "view-specific LoD assembly in place\n");
	    ret = 1;
	}
	if (!ret &&
	    (source->clearCompactInstanceVisibilityFrontier() <= 0 ||
	     !source->compactViewLodAssembly(payloads, &viewState) ||
	     viewState.findCadPresentation(source) != assembly ||
	     assembly->isInstanceHidden(ids[0]))) {
	    printf("FAIL: compact retained visibility did not restore the "
		   "view-specific LoD instance\n");
	    ret = 1;
	}
	if (!ret) {
	    BObolCompactInstanceHandle handle;
	    BObolCompactInstanceSummary summary;
	    const int selectionChanged =
		source->syncCompactInstanceSelectedPaths(paths);
	    const bool presentationWasStale =
		source->currentCompactViewLodAssembly(&viewState) == NULL;
	    auto *selectedAssembly =
		source->compactViewLodAssembly(payloads, &viewState);
	    const std::optional<Obol::InstanceRecord> selectedRecord =
		assembly->getInstanceRecord(ids[0]);
	    const bool selectedStyleCorrect = selectedRecord &&
		selectedRecord->style.hasColorOverride &&
		fabsf(selectedRecord->style.color[0] - 0.25f) < 1.0e-6f &&
		fabsf(selectedRecord->style.color[1] - 0.75f) < 1.0e-6f &&
		fabsf(selectedRecord->style.color[2] - 0.125f) < 1.0e-6f;
	    if (selectionChanged <= 0 ||
		!presentationWasStale ||
		!source->getCompactInstanceHandle(0, handle) ||
		!source->getCompactInstanceSummary(handle, summary) ||
		!summary.selected ||
		!selectedAssembly ||
		viewState.findCadPresentation(source) != assembly ||
		!source->currentCompactViewLodAssembly(&viewState) ||
		!selectedStyleCorrect ||
		viewState.findCadForResult(results[0]) != payloadBefore ||
		viewState.cadPayloadCount() != payloadCount) {
		printf("FAIL: compact retained selection did not refresh the "
		       "view-specific LoD assembly\n");
		ret = 1;
	    }
	}
	if (!ret) {
	    /*
	     * Erasing a selected subpath only hides its retained occurrence.
	     * Restoring that visibility must not silently replace the selected
	     * instance style with the normal style.  This is the exact sequence
	     * used by qged's hierarchy erase/draw workflow.
	     */
	    std::vector<SbBool> selectedHidden = {FALSE};
	    if (source->setCompactInstanceVisibilityOverrides(paths,
		    selectedHidden) <= 0 ||
		!source->compactViewLodAssembly(payloads, &viewState) ||
		!assembly->isInstanceHidden(ids[0])) {
		printf("FAIL: selected compact occurrence was not hidden in place\n");
		ret = 1;
	    }
	    if (!ret &&
		(source->clearCompactInstanceVisibilityFrontier() <= 0 ||
		 !source->compactViewLodAssembly(payloads, &viewState) ||
		 assembly->isInstanceHidden(ids[0]))) {
		printf("FAIL: selected compact occurrence was not restored in place\n");
		ret = 1;
	    }
	    const std::optional<Obol::InstanceRecord> restoredRecord =
		assembly->getInstanceRecord(ids[0]);
	    const bool restoredSelectedStyle = restoredRecord &&
		restoredRecord->style.hasColorOverride &&
		fabsf(restoredRecord->style.color[0] - 0.25f) < 1.0e-6f &&
		fabsf(restoredRecord->style.color[1] - 0.75f) < 1.0e-6f &&
		fabsf(restoredRecord->style.color[2] - 0.125f) < 1.0e-6f;
	    if (!ret && !restoredSelectedStyle) {
		printf("FAIL: selected compact occurrence lost its style across "
		       "hide/restore\n");
		ret = 1;
	    }
	}
	if (!ret) {
	    const std::vector<SbString> noPaths;
	    (void)source->syncCompactInstanceSelectedPaths(noPaths);
	    (void)source->compactViewLodAssembly(payloads, &viewState);
	    const std::optional<Obol::InstanceRecord> deselectedRecord =
		assembly->getInstanceRecord(ids[0]);
	    const bool stillSelectedStyle = deselectedRecord &&
		fabsf(deselectedRecord->style.color[0] - 0.25f) < 1.0e-6f &&
		fabsf(deselectedRecord->style.color[1] - 0.75f) < 1.0e-6f &&
		fabsf(deselectedRecord->style.color[2] - 0.125f) < 1.0e-6f;
	    if (!deselectedRecord || stillSelectedStyle) {
		printf("FAIL: compact retained deselection did not restore the "
		       "normal view-specific LoD style\n");
		ret = 1;
	    }
	}
    }

    /* A retained population prefix may already contain every point and face
     * required by finer coordinate-only cuts.  The Obol part must advertise
     * that drawable frontier, not clamp the instance back to the last level
     * that happened to add array entries. */
    if (!ret) {
	const BObolViewLodState::CadPayload *payload =
	    viewState.findCadForResult(results[0]);
	int drawableLevel = payload && payload->progressiveMesh ?
	    payload->progressiveMesh->residentLevel() : -1;
	if (payload && payload->progressiveMesh) {
	    for (int level = drawableLevel + 1;
		 level <= payload->progressiveMesh->maximumLevel(); ++level) {
		if (!payload->progressiveMesh->canDrawLevel(level))
		    break;
		drawableLevel = level;
	    }
	}
	if (!payload || !payload->progressiveMesh ||
	    drawableLevel <= payload->progressiveMesh->residentLevel()) {
	    printf("FAIL: compact coordinate-only LoD fixture has no drawable plateau\n");
	    ret = 1;
	} else if (!viewState.retargetCadPayload(payload, drawableLevel,
		drawableLevel, 63, 64)) {
	    printf("FAIL: compact coordinate-only LoD cut could not be retargeted\n");
	    ret = 1;
	} else {
	    std::vector<const BObolViewLodState::CadPayload *> payloads;
	    viewState.findCadPayloads(source, payloads);
	    (void)source->compactViewLodAssembly(payloads, &viewState);
	    SoCADAssembly *assembly = viewState.findCadPresentation(source);
	    std::vector<Obol::InstanceId> ids =
		assembly ? assembly->instanceIds() :
		std::vector<Obol::InstanceId>();
	    std::optional<Obol::InstanceRecord> record =
		ids.size() == 1 ? assembly->getInstanceRecord(ids[0]) :
		std::optional<Obol::InstanceRecord>();
	    const Obol::PartGeometry *geometry =
		record ? assembly->partGeometry(record->part) : NULL;
	    if (!record || record->lodLevel != drawableLevel || !geometry ||
		!geometry->shaded ||
		geometry->shaded->progressiveResidentLevel < drawableLevel) {
		printf("FAIL: Obol clamped a drawable coordinate-only PoP cut "
		       "to its population level\n");
		ret = 1;
	    }
	}
	if (!ret) {
	    const BObolViewLodState::CadPayload *richer =
		viewState.findCadForResult(results[0]);
	    if (!viewState.retargetCadPayload(richer,
		    results[0].geometry.activeLevel,
		    results[0].request.requestedLevel,
		    results[0].request.viewRevision.value(),
		    results[0].request.policyRevision.value())) {
		printf("FAIL: compact coordinate-only LoD fixture could not "
		       "restore its original cut\n");
		ret = 1;
	    }
	}
    }

    /* Normal selection is presentation policy, not LoD residency policy.
     * Switching it must invalidate only the renderer assembly while retaining
     * the exact installed PoP payload and its resident cut. */
    if (!ret) {
	std::vector<const BObolViewLodState::CadPayload *> payloads;
	viewState.findCadPayloads(source, payloads);
	const BObolViewLodState::CadPayload *before =
	    viewState.findCadForResult(results[0]);
	const size_t payloadCount = viewState.cadPayloadCount();
	const uint64_t cadRevision = viewState.cadRevision();
	SoCADAssembly *authoredAssembly =
	    viewState.findCadPresentation(source);
	std::vector<Obol::InstanceId> authoredIds =
	    authoredAssembly ? authoredAssembly->instanceIds() :
	    std::vector<Obol::InstanceId>();
	std::optional<Obol::InstanceRecord> authoredRecord =
	    authoredIds.size() == 1 ?
	    authoredAssembly->getInstanceRecord(authoredIds[0]) :
	    std::optional<Obol::InstanceRecord>();
	const Obol::PartGeometry *authoredGeometry = authoredRecord ?
	    authoredAssembly->partGeometry(authoredRecord->part) : NULL;
	std::vector<SbVec3f> authoredNormals =
	    authoredGeometry && authoredGeometry->shaded ?
	    authoredGeometry->shaded->normals : std::vector<SbVec3f>();
	const bool authoredStayedIndexed =
	    authoredGeometry && authoredGeometry->shaded &&
	    authoredGeometry->shaded->positions.size() ==
		results[0].counts.pointCount &&
	    authoredGeometry->shaded->indices.size() ==
		results[0].counts.faceCount * 3;
	/* No authored normals is a valid renderer contract.  The default
	 * presentation must preserve the shared PoP index topology and let Obol
	 * derive face normals, rather than expanding to three vertices per
	 * triangle. */
	if (!authoredAssembly || !before || !authoredNormals.empty() ||
	    !authoredStayedIndexed) {
	    printf("FAIL: normal-policy test requires a resident presentation\n");
	    ret = 1;
	} else {
	    viewState.setNormalStyle(BObolViewLodState::NORMAL_SMOOTH, 180.0f);
	    const BObolViewLodState::CadPayload *after =
		viewState.findCadForResult(results[0]);
	    if (viewState.findCadPresentation(source) ||
		after != before ||
		viewState.cadPayloadCount() != payloadCount ||
		viewState.cadRevision() != cadRevision ||
		viewState.getNormalStyle() !=
		    BObolViewLodState::NORMAL_SMOOTH ||
		fabsf(viewState.getNormalCreaseAngle() - 180.0f) > 0.0001f) {
		printf("FAIL: normal policy discarded or mutated resident PoP data\n");
		ret = 1;
	    }
	}
	if (!ret) {
	    (void)source->compactViewLodAssembly(payloads, &viewState);
	    SoCADAssembly *smoothAssembly =
		viewState.findCadPresentation(source);
	    std::vector<Obol::InstanceId> smoothIds =
		smoothAssembly ? smoothAssembly->instanceIds() :
		std::vector<Obol::InstanceId>();
	    std::optional<Obol::InstanceRecord> smoothRecord =
		smoothIds.size() == 1 ?
		smoothAssembly->getInstanceRecord(smoothIds[0]) :
		std::optional<Obol::InstanceRecord>();
	    const Obol::PartGeometry *smoothGeometry = smoothRecord ?
		smoothAssembly->partGeometry(smoothRecord->part) : NULL;
	    const std::vector<SbVec3f> smoothNormals =
		smoothGeometry && smoothGeometry->shaded ?
		smoothGeometry->shaded->normals : std::vector<SbVec3f>();
	    bool normalsChanged =
		smoothNormals.size() != authoredNormals.size();
	    for (size_t i = 0;
		!normalsChanged && i < smoothNormals.size(); ++i) {
		if ((smoothNormals[i] - authoredNormals[i]).sqrLength() >
		    1.0e-8f)
		    normalsChanged = true;
	    }
	    if (!smoothAssembly || smoothNormals.empty() ||
		!normalsChanged) {
		printf("FAIL: normal policy did not rebuild the renderer presentation\n");
		ret = 1;
	    }
	}
    }

    if (!ret) {
	source->representationMode =
	    SoBRLDatabaseSource::REPRESENTATION_WIRE;
	std::vector<const BObolViewLodState::CadPayload *> payloads;
	viewState.findCadPayloads(source, payloads);
	const BObolViewLodState::CadPayload *before =
	    viewState.findCadForResult(results[0]);
	(void)source->compactViewLodAssembly(payloads, &viewState);
	SoCADAssembly *assembly = viewState.findCadPresentation(source);
	std::vector<Obol::InstanceId> ids =
	    assembly ? assembly->instanceIds() :
	    std::vector<Obol::InstanceId>();
	std::optional<Obol::InstanceRecord> record =
	    ids.size() == 1 ? assembly->getInstanceRecord(ids[0]) :
	    std::optional<Obol::InstanceRecord>();
	const Obol::PartGeometry *geometry =
	    record ? assembly->partGeometry(record->part) : NULL;
	const BObolViewLodState::CadPayload *after =
	    viewState.findCadForResult(results[0]);
	if (!assembly || ids.size() != 1 || !geometry ||
	    geometry->shaded || !geometry->wire ||
	    !geometry->wire->isProgressive() || !record ||
	    record->lodLevel != results[0].geometry.activeLevel ||
	    geometry->wire->segmentCountAtLevel(record->lodLevel) !=
		results[0].counts.faceCount * 3 ||
	    before != after || !after ||
	    after->activeLevel != results[0].geometry.activeLevel ||
	    after->counts.faceCount != results[0].counts.faceCount) {
	    printf("FAIL: wire compact presentation did not preserve the shaded PoP cut\n");
	    ret = 1;
	}
    }

    if (!ret) {
	SoBRLMeshLodSubmitAction wireSubmit;
	wireSubmit.setService(&service);
	wireSubmit.setViewLodState(&viewState);
	wireSubmit.setDatabase(dbip, "db://compact-projected-test", 2026);
	wireSubmit.setViewInfo(&view);
	wireSubmit.setViewVolume(&volume, 1.0f);
	wireSubmit.setGeneration(service.beginGeneration());
	wireSubmit.setRevisions(63, 62);
	wireSubmit.apply(root);
	if (wireSubmit.getVisitedMeshCount() != 1 ||
	    wireSubmit.getSubmittedTaskCount() != 0 ||
	    wireSubmit.getSkippedMeshCount() != 1 ||
	    service.queuedResultCountForDiagnostics() != 0) {
	    printf("FAIL: wire mode restarted an already-resident PoP request\n");
	    ret = 1;
	}
    }

    if (!ret) {
	/* A visible piece of an object larger than the window still needs error
	 * control based on the object's full projected scale. */
	SbMatrix largePlacement;
	largePlacement.setScale(SbVec3f(100.0f, 100.0f, 1.0f));
	(void)source->setPlacementState(TRUE, largePlacement, FALSE,
	    SbVec3f(0.0f, 0.0f, 0.0f), FALSE, 0.0f);
	camera->position = SbVec3f(0.0f, 0.0f, 5.0f);
	const SbViewVolume partialVolume =
	    camera->getViewVolume(640.0f / 480.0f);
	SoBRLMeshLodSubmitAction partialSubmit;
	partialSubmit.setService(&service);
	partialSubmit.setDatabase(dbip, "db://compact-projected-test", 2026);
	partialSubmit.setViewInfo(&view);
	partialSubmit.setViewVolume(&partialVolume, 1.0f);
	partialSubmit.setGeneration(service.beginGeneration());
	partialSubmit.setRevisions(64, 62);
	partialSubmit.apply(root);
	std::vector<BObolLodResult> partialResults;
	if (partialSubmit.getSubmittedTaskCount() != 1 ||
	    wait_for_service(service) ||
	    service.drainResults(partialResults) != 1 ||
	    partialResults.size() != 1 ||
	    partialResults[0].request.projectedPixelDiameter < 4790.0f ||
	    partialResults[0].request.projectedPixelDiameter > 4810.0f ||
	    partialResults[0].request.requestedLevel != 13) {
	    printf("FAIL: partially visible mesh demand was clipped to viewport size\n");
	    ret = 1;
	}
    }

    camera->unref();
    service.stop();
    root->setViewLodState(NULL);
    root->unref();
    bobol_mesh_lod_cache_clear_database(dbip);
    db_close(dbip);
    bu_file_delete(dbpath);
    bu_dirclear(cacheDir);
    return ret;
}

static int
test_view_lod_mesh_eviction_preserves_proxy(void)
{
    SoBRLMeshShape *mesh = make_mesh("/eviction/mesh.bot", "mesh.bot");
    mesh->ref();
    BObolViewLodState viewState;
    BObolLodRequest meshRequest =
	make_request("/eviction/mesh.bot", "mesh.bot");
    BObolLodResult meshResult = mesh_payload_result(meshRequest);
    BObolLodRequest proxyRequest = meshRequest;
    proxyRequest.providerId = "rt_proxy_aabb";
    proxyRequest.providerVersion = "rt-proxy-v1";
    proxyRequest.qualityTier = BOBOL_LOD_QUALITY_PROXY;
    BObolLodResult proxyResult = bobol_lod_aabb_result(proxyRequest,
	proxyRequest.bounds, NULL);

    int ret = 0;
    if (!viewState.applyDisplayResult(mesh, meshResult) ||
	!viewState.applyDisplayResult(mesh, proxyResult) ||
	viewState.meshPayloadCount() != 1 ||
	viewState.proxyPayloadCount(BOBOL_LOD_PROXY_AABB) != 1) {
	printf("FAIL: mesh eviction fixture did not install display payloads\n");
	ret = 1;
    }

    unsigned int evicted = 0;
	const size_t freed = viewState.evictDisplayMeshPayloads(&evicted);
    if (!ret && (freed == 0 || evicted != 1 || viewState.findMesh(mesh) ||
	!viewState.findProxy(mesh) || viewState.meshPayloadCount() != 0 ||
	viewState.proxyPayloadCount(BOBOL_LOD_PROXY_AABB) != 1)) {
	printf("FAIL: mesh eviction discarded the retained proxy frontier\n");
	ret = 1;
    }

    mesh->unref();
    return ret;
}

static int
test_view_controller_compact_direct_result_route(void)
{
    SoSeparator *sceneRoot = new SoSeparator;
    sceneRoot->ref();
    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->path = "direct-route.c";
    source->instanceKey = "direct-route";
    source->sourceRevision = 17;
    sceneRoot->addChild(source);

    BObolCompactOccurrence occurrence;
    occurrence.geometry = compact_duplicate_proxy_geometry();
    occurrence.summary.valid = TRUE;
    occurrence.summary.shapeKind =
	BObolRealizedShapeSummary::SHAPE_VLIST;
    occurrence.summary.path = "direct-route.c/leaf.s";
    occurrence.summary.sourceName = "leaf.s";
    occurrence.summary.sourceType = "proxy";
    occurrence.summary.sourceId = 17;
    occurrence.summary.visible = TRUE;
    occurrence.summary.selectable = TRUE;
    occurrence.localTransform = SbMatrix::identity();
    occurrence.lodBacked = TRUE;
    occurrence.occurrenceIndex = 0;
    std::vector<BObolCompactOccurrence> occurrences(1, occurrence);
    if (source->setCompactOccurrenceRegistry(occurrences) != 1) {
	printf("FAIL: compact direct-route occurrence setup\n");
	sceneRoot->unref();
	return 1;
    }

    BObolCompactInstanceHandle handle;
    BObolCompactInstanceSummary summary;
    if (!source->getCompactInstanceHandle(0, handle) ||
	!source->getCompactInstanceSummary(handle, summary) ||
	summary.sourceInstanceKey.getLength() == 0) {
	printf("FAIL: compact direct-route occurrence identity\n");
	sceneRoot->unref();
	return 1;
    }

    SoSeparator *unrelatedRenderRoot = new SoSeparator;
    unrelatedRenderRoot->ref();
    BObolLodService service;
    int ret = 0;
    {
	BObolViewController controller(sceneRoot, NULL);
	/* The source is deliberately absent from the render graph.  A generic
	 * update action cannot discover it; only the scene-controller routing
	 * index can apply this result. */
	controller.setRenderSceneRoot(unrelatedRenderRoot);
	controller.setLodService(&service);
	controller.setLodPolicyRevision(71);
	if (!service.start(1, TRUE)) {
	    printf("FAIL: compact direct-route service did not start\n");
	    ret = 1;
	}

	BObolLodRequest request =
	    make_request("direct-route.c/leaf.s", "leaf.s");
	request.sourceRevision = source->sourceRevision.getValue();
	request.sourceRoutingId = source->getCompactSourceRoutingId();
	request.occurrenceKey = summary.sourceInstanceKey;
	request.viewRevision = controller.getLodViewRevision();
	request.policyRevision = controller.getLodPolicyRevision();
	request.drawMode = BOBOL_LOD_DRAW_WIRE;
	request.qualityTier = BOBOL_LOD_QUALITY_PROXY;

	if (!ret) {
	    BObolLodTask task;
	    task.generation = controller.beginLodGeneration();
	    task.request = request;
	    task.realize = aabb_payload_task;
	    if (service.submit(task) == 0 || wait_for_service(service)) {
		ret = 1;
	    } else if (controller.processPendingLodResults(1, 0) != 1) {
		printf("FAIL: compact direct-route result was not applied\n");
		ret = 1;
	    }
	}

	const BObolViewLodState::CadPayload *payload =
	    controller.getViewLodState()->findCadForOccurrence(
		source, summary.sourceInstanceKey);
	if (!ret && (controller.getLastLodMatchedResultCount() != 1 ||
	    controller.getLastLodAppliedResultCount() != 1 ||
	    controller.getLastLodRejectedResultCount() != 0 ||
	    controller.getLastLodUnmatchedResultCount() != 0 ||
	    !payload || payload->resultKind != BOBOL_LOD_RESULT_AABB ||
	    payload->proxy.kind != BOBOL_LOD_PROXY_AABB)) {
	    printf("FAIL: compact result did not use direct indexed routing "
		"(matched=%u applied=%u rejected=%u unmatched=%u payload=%d "
		"kind=%d proxy=%d diagnostics=%s)\n",
		controller.getLastLodMatchedResultCount(),
		controller.getLastLodAppliedResultCount(),
		controller.getLastLodRejectedResultCount(),
		controller.getLastLodUnmatchedResultCount(),
		payload ? 1 : 0,
		payload ? payload->resultKind : -1,
		payload ? payload->proxy.kind : -1,
		controller.getLastLodDiagnostics().getString());
	    ret = 1;
	}
	if (!ret) {
	    /*
	     * Refining/replacing one compact occurrence updates its one
	     * owner-thread metadata slot.  The payload allocation and direct
	     * source/occurrence lookup must remain stable; only immutable
	     * renderer geometry generations are exchanged.
	     */
	    const BObolViewLodState::CadPayload *payloadBefore = payload;
	    BObolLodTask replacement;
	    replacement.generation = controller.beginLodGeneration();
	    replacement.request = request;
	    replacement.realize = aabb_payload_task;
	    if (service.submit(replacement) == 0 ||
		wait_for_service(service) ||
		controller.processPendingLodResults(1, 0) != 1) {
		printf("FAIL: compact direct-route replacement was not applied\n");
		ret = 1;
	    } else {
		payload = controller.getViewLodState()->findCadForOccurrence(
		    source, summary.sourceInstanceKey);
		if (!payload || payload != payloadBefore ||
		    payload->resultKind != BOBOL_LOD_RESULT_AABB) {
		    printf("FAIL: compact direct-route replacement did not "
			   "reuse its authoritative payload slot\n");
		    ret = 1;
		}
	    }
	}
	if (!ret) {
	    /*
	     * A cold compact source has no preceding mesh.  Its first useful
	     * result is nevertheless only a nonterminal PoP prefix and must
	     * keep the render/refinement barrier armed.  A later camera event
	     * must not be required to make progress.
	     */
	    if (!controller.getViewLodState()->removeCadPayload(payload)) {
		printf("FAIL: compact partial-route fixture did not clear its proxy\n");
		ret = 1;
	    } else {
		controller.clearRenderRequest();
		request.requestedLevel = 2;
		BObolLodTask task;
		task.generation = controller.beginLodGeneration();
		task.request = request;
		task.realize = mesh_payload_coarse_task;
		if (service.submit(task) == 0 || wait_for_service(service)) {
		    ret = 1;
		} else if (controller.processPendingLodResults(1, 0) != 1 ||
		    controller.getLastLodAppliedResultCount() != 1 ||
		    !controller.isRenderRequested() ||
		    !controller.hasPendingLodRefinementFrame()) {
		    printf("FAIL: first compact nonterminal mesh did not "
			"schedule its presentation/refinement barrier\n");
		    ret = 1;
		} else {
		    /* A coalesced/consumed host request must not strand the
		     * barrier.  The next progressive pump owns the invariant
		     * and reasserts a presentation edge. */
		    controller.clearRenderRequest();
		    (void)controller.advanceProgressiveWork(NULL, NULL);
		    if (!controller.isRenderRequested() ||
			!controller.hasProgressiveWorkPending()) {
			printf("FAIL: pending compact refinement barrier did "
			       "not reassert its host frame request\n");
			ret = 1;
		    }
		}
	    }
	}
	controller.setLodService(NULL);
    }
    service.stop();
    unrelatedRenderRoot->unref();
    sceneRoot->unref();
    return ret;
}

static int
test_compact_occurrence_lod_identity(void)
{
    BObolViewLodState viewState;
    SoBRLViewLodGroup *root = new SoBRLViewLodGroup;
    root->ref();
    root->setViewLodState(&viewState);

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->path = "duplicate-root.c";
    source->instanceKey = "duplicate-root";
    root->addChild(source);

    std::shared_ptr<const Obol::PartGeometry> geometry =
	compact_duplicate_proxy_geometry();
    BObolCompactOccurrence first;
    first.geometry = geometry;
    first.summary.valid = TRUE;
    first.summary.shapeKind = BObolRealizedShapeSummary::SHAPE_VLIST;
    first.summary.path = "duplicate-root.c/duplicate.s";
    first.summary.sourceName = "duplicate.s";
    first.summary.sourceType = "proxy";
    first.summary.sourceId = 91;
    first.summary.visible = TRUE;
    first.summary.selectable = TRUE;
    first.localTransform = SbMatrix::identity();
    first.lodBacked = TRUE;
    first.occurrenceIndex = 0;

    BObolCompactOccurrence second = first;
    second.localTransform.setTranslate(SbVec3f(10.0f, 0.0f, 0.0f));
    second.occurrenceIndex = 1;
    std::vector<BObolCompactOccurrence> occurrences;
    occurrences.push_back(first);
    occurrences.push_back(second);
    if (source->setCompactOccurrenceRegistry(occurrences) != 2) {
	printf("FAIL: compact occurrence LoD identity setup\n");
	root->setViewLodState(NULL);
	root->unref();
	return 1;
    }

    BObolCompactInstanceHandle firstHandle;
    BObolCompactInstanceHandle secondHandle;
    BObolCompactInstanceSummary firstSummary;
    BObolCompactInstanceSummary secondSummary;
    if (!source->getCompactInstanceHandle(0, firstHandle) ||
	!source->getCompactInstanceHandle(1, secondHandle) ||
	!source->getCompactInstanceSummary(firstHandle, firstSummary) ||
	!source->getCompactInstanceSummary(secondHandle, secondSummary) ||
	firstSummary.sourceInstanceKey.getLength() == 0 ||
	secondSummary.sourceInstanceKey.getLength() == 0 ||
	bu_strcmp(firstSummary.sourceInstanceKey.getString(),
		secondSummary.sourceInstanceKey.getString()) == 0) {
	printf("FAIL: compact occurrence LoD identity handles\n");
	root->setViewLodState(NULL);
	root->unref();
	return 1;
    }

    BObolLodRequest firstRequest =
	make_request("duplicate-root.c/duplicate.s", "duplicate.s");
    firstRequest.drawMode = BOBOL_LOD_DRAW_WIRE;
    firstRequest.qualityTier = BOBOL_LOD_QUALITY_PROXY;
    firstRequest.occurrenceKey = firstSummary.sourceInstanceKey;
    BObolLodRequest secondRequest = firstRequest;
    secondRequest.occurrenceKey = secondSummary.sourceInstanceKey;
    const SbBox3f firstBounds(SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(2.0f, 2.0f, 2.0f));
    const SbBox3f secondBounds(SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(3.0f, 3.0f, 3.0f));
    BObolLodResult firstResult =
	bobol_lod_aabb_result(firstRequest, firstBounds, NULL);
    BObolLodResult secondResult =
	bobol_lod_aabb_result(secondRequest, secondBounds, NULL);

    SoBRLLodUpdateAction update;
    update.setViewLodState(&viewState);
    update.addResult(firstResult);
    update.addResult(secondResult);
    update.apply(root);

    int ret = 0;
    if (update.getMatchedResultCount() != 2 ||
	update.getAppliedResultCount() != 2 ||
	viewState.cadPayloadCount() != 2) {
	printf("FAIL: compact occurrence LoD results did not bind independently\n");
	ret = 1;
    }
    if (!ret) {
	/* The OBB must replace, rather than coexist with, the first
	 * occurrence's AABB.  A late coarse result must not displace it. */
	BObolLodRequest obbRequest = firstRequest;
	obbRequest.providerId = "rt_proxy_obb";
	obbRequest.providerVersion = "rt-proxy-v1";
	BObolLodProxy obb;
	obb.kind = BOBOL_LOD_PROXY_OBB;
	obb.bounds = firstBounds;
	obb.center = SbVec3f(1.0f, 1.0f, 1.0f);
	obb.axisX = SbVec3f(1.0f, 0.0f, 0.0f);
	obb.axisY = SbVec3f(0.0f, 1.0f, 0.0f);
	obb.axisZ = SbVec3f(0.0f, 0.0f, 1.0f);
	obb.halfExtents = SbVec3f(1.0f, 1.0f, 1.0f);
	BObolLodResult obbResult =
	    bobol_lod_proxy_result(obbRequest, obb, NULL);
	SoBRLLodUpdateAction obbUpdate;
	obbUpdate.setViewLodState(&viewState);
	obbUpdate.addResult(obbResult);
	obbUpdate.apply(root);
	const BObolViewLodState::CadPayload *obbPayload =
	    viewState.findCadForResult(obbResult);
	if (obbUpdate.getMatchedResultCount() != 1 ||
	    obbUpdate.getAppliedResultCount() != 1 ||
	    viewState.cadPayloadCount() != 2 || !obbPayload ||
	    obbPayload->proxy.kind != BOBOL_LOD_PROXY_OBB) {
	    printf("FAIL: compact OBB refinement did not replace its AABB\n");
	    ret = 1;
	}

	BObolLodRequest lateAabbRequest = firstRequest;
	lateAabbRequest.providerId = "rt_proxy_aabb";
	lateAabbRequest.providerVersion = "rt-proxy-v1";
	BObolLodResult lateAabb =
	    bobol_lod_aabb_result(lateAabbRequest, firstBounds, NULL);
	SoBRLLodUpdateAction lateUpdate;
	lateUpdate.setViewLodState(&viewState);
	lateUpdate.addResult(lateAabb);
	lateUpdate.apply(root);
	obbPayload = viewState.findCadForResult(obbResult);
	if (lateUpdate.getMatchedResultCount() != 1 ||
	    lateUpdate.getAppliedResultCount() != 1 ||
	    viewState.cadPayloadCount() != 2 || !obbPayload ||
	    obbPayload->proxy.kind != BOBOL_LOD_PROXY_OBB) {
	    printf("FAIL: late compact AABB displaced a better OBB\n");
	    ret = 1;
	}
    }
    if (!ret) {
	SbViewportRegion viewport(100, 100);
	SoGetBoundingBoxAction boundsAction(viewport);
	boundsAction.apply(root);
	const SbBox3f bounds = boundsAction.getBoundingBox();
	if (bounds.isEmpty() || fabs(bounds.getMin()[0]) > 0.0001f ||
	    fabs(bounds.getMax()[0] - 13.0f) > 0.0001f ||
	    fabs(bounds.getMax()[1] - 3.0f) > 0.0001f) {
	    printf("FAIL: compact occurrence LoD result aliased duplicate paths\n");
	    ret = 1;
	}
    }

    root->setViewLodState(NULL);
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

    bobol_init(NULL);

    if (test_compact_staged_source_lease())
	return 1;
    if (test_view_controller_fixed_gl_light_order())
	return 1;
    if (test_view_controller_presentation_deadline_contract())
	return 1;
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
    if (test_view_lod_mesh_eviction_preserves_proxy())
	return 1;
    if (test_mesh_lod_submit_action())
	return 1;
    if (test_compact_many_leaf_scene_admission())
	return 1;
    if (test_compact_native_mesh_not_pop_submitted())
	return 1;
    if (test_view_controller_compact_append_preserves_submission_cursor())
	return 1;
    if (test_compact_aabb_stream_upgrade())
	return 1;
    if (test_compact_mesh_lod_projection_and_mode_parity())
	return 1;
    if (test_view_controller_source_backed_multi_source_exact_submit())
	return 1;
    if (test_view_controller_source_backed_partial_ready_submit())
	return 1;
    if (test_view_controller_shared_lod_is_view_local())
	return 1;
    if (test_view_controller_progressive_lod_results())
	return 1;
    if (test_view_controller_large_result_transfer())
	return 1;
    if (test_view_local_lod_only_pick())
	return 1;
    if (test_view_controller_compact_direct_result_route())
	return 1;
    if (test_compact_occurrence_lod_identity())
	return 1;
    if (test_view_controller_lod_submit_and_apply())
	return 1;
    if (test_view_controller_lod_source_churn())
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
