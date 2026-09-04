/*                T E S T _ L O D _ S E R V I C E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "BObol.h"
#include "BObol/BLodRealization.h"
#include "BObol/BLodService.h"
#include "BObol/BMeshLodCache.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BLodService.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "wdb.h"

#include <Obol/cad/SoCADAssembly.h>

#include <Inventor/SbBox.h>
#include <Inventor/SbVec3f.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <limits>
#include <mutex>
#include <numeric>
#include <stdio.h>
#include <string.h>
#include <thread>
#include <vector>

static BObolLodRequest
make_request(const char *name)
{
    BObolLodRequest request;

    request.databaseId = "db://lod-service-test";
    request.databaseRevision = 9;
    request.sourceRevision = 17;
    request.sourceContentHash = 0x5577;
    request.objectPath = name;
    request.objectName = name;
    request.viewRevision = 23;
    request.policyRevision = 29;
    request.drawMode = BOBOL_LOD_DRAW_SHADED;
    request.lodPolicy = 3;
    request.providerId = "service-test";
    request.providerVersion = "1";
    request.qualityTier = BOBOL_LOD_QUALITY_FAST_DISPLAY;
    request.bounds = SbBox3f(SbVec3f(-1.0f, -1.0f, -1.0f),
			     SbVec3f(1.0f, 1.0f, 1.0f));
    request.sourceCounts.faceCount = 100;
    request.sourceCounts.pointCount = 101;

    return request;
}

struct ServiceTestContext {
    std::mutex mutex;
    std::vector<int> executionOrder;
    std::vector<int> cacheWriteOrder;
};

struct TaskData {
    ServiceTestContext *context;
    int value;
};

static BObolLodResult
ready_task(const BObolLodRequest &request, void *userData)
{
    TaskData *data = static_cast<TaskData *>(userData);

    {
	std::lock_guard<std::mutex> lock(data->context->mutex);
	data->context->executionOrder.push_back(data->value);
    }

    BObolLodResult result;
    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = BOBOL_LOD_RESULT_AABB;
    result.qualityTier = request.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.terminal = TRUE;
    result.counts.faceCount = (uint64_t)data->value;
    result.bounds = request.bounds;
    return result;
}

static void
cache_write(const BObolLodResult &result, void *userData)
{
    ServiceTestContext *context = static_cast<ServiceTestContext *>(userData);
    std::lock_guard<std::mutex> lock(context->mutex);
    context->cacheWriteOrder.push_back((int)result.counts.faceCount);
}

static int
wait_for_settled(BObolLodService &service, size_t expectedResults)
{
    for (int i = 0; i < 400; i++) {
	if (service.inFlightCount() == 0 &&
	    service.queuedCacheWriteCountForDiagnostics() == 0 &&
	    service.queuedResultCountForDiagnostics() >= expectedResults)
	    return 0;
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    printf("FAIL: LoD service did not settle: inFlight=%zu results=%zu cache=%zu\n",
	   service.inFlightCount(),
	   service.queuedResultCountForDiagnostics(),
	   service.queuedCacheWriteCountForDiagnostics());
    return 1;
}

static int
wait_for_resident_compaction(BObolLodService &service)
{
    for (int i = 0; i < 2000; i++) {
	if (service.
		pendingResidentMeshCompactionCountForDiagnostics() == 0)
	    return 0;
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    printf("FAIL: resident LoD compaction did not settle: pending=%zu\n",
	   service.pendingResidentMeshCompactionCountForDiagnostics());
    return 1;
}

static size_t
provider_param_count(const BObolLodRequest &request, const char *name)
{
    size_t count = 0;

    if (!name)
	return 0;

    for (size_t i = 0; i < request.providerParams.size(); i++) {
	if (bu_strcmp(request.providerParams[i].name.getString(), name) == 0)
	    count++;
    }

    return count;
}

static int
test_dependency_order_and_cache_write(void)
{
    BObolLodService service;
    ServiceTestContext context;
    TaskData firstData;
    TaskData secondData;

    firstData.context = &context;
    firstData.value = 1;
    secondData.context = &context;
    secondData.value = 2;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start\n");
	return 1;
    }
    if (!service.isRunning()) {
	printf("FAIL: LoD service did not report running\n");
	return 1;
    }

    uint64_t generation = service.beginGeneration();

    BObolLodTask first;
    first.generation = generation;
    first.request = make_request("/first.bot");
    first.realize = ready_task;
    first.realizeData = &firstData;
    first.writeCache = TRUE;
    first.cacheWrite = cache_write;
    first.cacheWriteData = &context;

    BObolLodTask second;
    second.generation = generation;
    second.request = make_request("/second.bot");
    second.realize = ready_task;
    second.realizeData = &secondData;
    second.writeCache = TRUE;
    second.cacheWrite = cache_write;
    second.cacheWriteData = &context;

    uint64_t firstId = service.submit(first);
    second.addDependency(firstId);
    uint64_t secondId = service.submit(second);
    if (firstId == 0 || secondId == 0 || firstId == secondId) {
	printf("FAIL: LoD service did not assign task ids\n");
	return 1;
    }

    if (wait_for_settled(service, 2))
	return 1;

    std::vector<BObolLodResult> results;
    if (service.drainResults(results) != 2 || results.size() != 2) {
	printf("FAIL: LoD service did not drain expected results\n");
	return 1;
    }
    if (service.queuedResultCountForDiagnostics() != 0) {
	printf("FAIL: LoD service retained drained results\n");
	return 1;
    }

    if (results[0].counts.faceCount != 1 ||
	results[1].counts.faceCount != 2 ||
	!bobol_lod_result_matches_request(results[0], first.request) ||
	!bobol_lod_result_matches_request(results[1], second.request)) {
	printf("FAIL: LoD service result ordering or request identity failed\n");
	return 1;
    }

    {
	std::lock_guard<std::mutex> lock(context.mutex);
	if (context.executionOrder.size() != 2 ||
	    context.executionOrder[0] != 1 ||
	    context.executionOrder[1] != 2 ||
	    context.cacheWriteOrder.size() != 2 ||
	    context.cacheWriteOrder[0] != 1 ||
	    context.cacheWriteOrder[1] != 2) {
	    printf("FAIL: LoD service dependency/cache-write order failed\n");
	    return 1;
	}
    }

    service.stop();
    if (service.isRunning()) {
	printf("FAIL: LoD service still reports running after stop\n");
	return 1;
    }

    return 0;
}

static int
test_filtered_result_drain(void)
{
    ServiceTestContext context;
    TaskData firstData;
    TaskData secondData;
    TaskData thirdData;
    BObolLodService service;

    firstData.context = &context;
    firstData.value = 1;
    secondData.context = &context;
    secondData.value = 2;
    thirdData.context = &context;
    thirdData.value = 3;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD filtered-drain service did not start\n");
	return 1;
    }

    uint64_t generation = service.beginGeneration();
    BObolLodTask first;
    first.generation = generation;
    first.request = make_request("/filtered-first.bot");
    first.realize = ready_task;
    first.realizeData = &firstData;

    BObolLodTask second;
    second.generation = generation;
    second.request = make_request("/filtered-second.bot");
    second.realize = ready_task;
    second.realizeData = &secondData;

    BObolLodTask third;
    third.generation = generation;
    third.request = make_request("/filtered-third.bot");
    third.realize = ready_task;
    third.realizeData = &thirdData;

    if (service.submit(first) == 0 ||
	service.submit(second) == 0 ||
	service.submit(third) == 0) {
	printf("FAIL: LoD filtered-drain service did not accept tasks\n");
	service.stop();
	return 1;
    }

    if (wait_for_settled(service, 3)) {
	service.stop();
	return 1;
    }

    std::vector<BObolLodRequest> requests;
    requests.push_back(second.request);
    std::vector<BObolLodResult> matched;
    if (service.drainMatchingResults(matched, requests) != 1 ||
	matched.size() != 1 ||
	!bobol_lod_result_matches_request(matched[0], second.request) ||
	service.queuedResultCountForDiagnostics() != 2) {
	printf("FAIL: LoD filtered result drain did not isolate requested result\n");
	service.stop();
	return 1;
    }

    std::vector<BObolLodResult> remaining;
    if (service.drainResults(remaining) != 2 ||
	remaining.size() != 2 ||
	!bobol_lod_result_matches_request(remaining[0], first.request) ||
	!bobol_lod_result_matches_request(remaining[1], third.request)) {
	printf("FAIL: LoD filtered result drain did not preserve unmatched queue order\n");
	service.stop();
	return 1;
    }

    service.stop();
    return 0;
}

static BObolLodResult
ready_mesh_task(const BObolLodRequest &request, void *userData)
{
    BObolLodResult result = ready_task(request, userData);
    result.resultKind = BOBOL_LOD_RESULT_MESH;
    result.mesh.points.push_back(SbVec3f(0.0f, 0.0f, 0.0f));
    result.mesh.points.push_back(SbVec3f(1.0f, 0.0f, 0.0f));
    result.mesh.points.push_back(SbVec3f(0.0f, 1.0f, 0.0f));
    result.mesh.coordIndex.push_back(0);
    result.mesh.coordIndex.push_back(1);
    result.mesh.coordIndex.push_back(2);
    return result;
}

static int
test_occurrence_result_coalescing(void)
{
    BObolLodService service;
    ServiceTestContext context;
    TaskData oldData{&context, 1};
    TaskData newData{&context, 2};
    TaskData meshData{&context, 3};
    TaskData routedData{&context, 4};
    if (!service.start(2, FALSE)) {
	printf("FAIL: LoD coalescing service did not start\n");
	return 1;
    }

    const uint64_t generation = service.beginGeneration();
    BObolLodTask oldTask;
    oldTask.generation = generation;
    oldTask.request = make_request("/coalesced.bot");
    oldTask.realize = ready_task;
    oldTask.realizeData = &oldData;
    oldTask.debugDelayMilliseconds = 80;

    BObolLodTask newTask = oldTask;
    newTask.request.viewRevision++;
    newTask.realizeData = &newData;
    newTask.debugDelayMilliseconds = 0;

    BObolLodTask meshTask = newTask;
    meshTask.realize = ready_mesh_task;
    meshTask.realizeData = &meshData;

    /*
     * A second live source may consume the same occurrence.  Its result must
     * neither be rejected as duplicate active work nor coalesced into the
     * first source's publication slot.
     */
    BObolLodTask routedTask = newTask;
    routedTask.request.sourceRoutingId = 91;
    routedTask.realizeData = &routedData;

    if (service.submit(oldTask) == 0 || service.submit(newTask) == 0 ||
	service.submit(meshTask) == 0 || service.submit(routedTask) == 0 ||
	wait_for_settled(service, 2)) {
	printf("FAIL: LoD service did not execute coalescing tasks\n");
	service.stop();
	return 1;
    }

    std::vector<BObolLodResult> results;
    service.drainResults(results);
    const BObolLodResult *aabb = NULL;
    const BObolLodResult *routedAabb = NULL;
    const BObolLodResult *mesh = NULL;
    for (const BObolLodResult &result : results) {
	if (result.resultKind == BOBOL_LOD_RESULT_AABB &&
	    result.request.sourceRoutingId == 91)
	    routedAabb = &result;
	else if (result.resultKind == BOBOL_LOD_RESULT_AABB)
	    aabb = &result;
	if (result.resultKind == BOBOL_LOD_RESULT_MESH)
	    mesh = &result;
    }
    if (results.size() != 3 || !aabb || !routedAabb || !mesh ||
	aabb->request.viewRevision != newTask.request.viewRevision ||
	aabb->counts.faceCount != 2 ||
	routedAabb->counts.faceCount != 4 ||
	mesh->counts.faceCount != 3 ||
	service.coalescedResultCountForDiagnostics() != 1 ||
	service.discardedStaleResultCountForDiagnostics() != 0) {
	printf("FAIL: LoD service did not retain newest per-occurrence stages\n");
	service.stop();
	return 1;
    }

    service.stop();
    return 0;
}

struct BlockingTaskData {
    std::mutex mutex;
    std::condition_variable cv;
    bool started;
    bool release;
    int calls;

    BlockingTaskData(void) : started(false), release(false), calls(0)
    {
    }
};

static BObolLodResult
blocking_task(const BObolLodRequest &request, void *userData)
{
    BlockingTaskData *data = static_cast<BlockingTaskData *>(userData);

    {
	std::unique_lock<std::mutex> lock(data->mutex);
	data->started = true;
	data->calls++;
	data->cv.notify_all();
	while (!data->release)
	    data->cv.wait(lock);
    }

    BObolLodResult result;
    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = BOBOL_LOD_RESULT_MESH;
    result.qualityTier = request.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.counts.faceCount = 7;
    return result;
}

static int
wait_for_blocking_task_started(BlockingTaskData &data)
{
    std::unique_lock<std::mutex> lock(data.mutex);
    if (!data.cv.wait_for(lock, std::chrono::seconds(2),
			  [&data] { return data.started; })) {
	printf("FAIL: LoD blocking task did not start\n");
	return 1;
    }
    return 0;
}

static int
test_generation_cancellation(void)
{
    BObolLodService service;
    BlockingTaskData blockingData;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for cancellation test\n");
	return 1;
    }

    uint64_t generation = service.beginGeneration();
    BObolLodTask task;
    task.generation = generation;
    task.request = make_request("/cancelled.bot");
    task.realize = blocking_task;
    task.realizeData = &blockingData;

    if (service.submit(task) == 0) {
	printf("FAIL: LoD service did not accept cancellable task\n");
	return 1;
    }

    if (wait_for_blocking_task_started(blockingData))
	return 1;

    service.cancelGeneration(generation);
    if (!service.isGenerationCancelled(generation)) {
	printf("FAIL: LoD service did not record generation cancellation\n");
	return 1;
    }

    {
	std::lock_guard<std::mutex> lock(blockingData.mutex);
	blockingData.release = true;
	blockingData.cv.notify_all();
    }

    for (int i = 0; i < 400 && service.inFlightCount() != 0; i++)
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    if (service.inFlightCount() != 0 ||
	service.queuedResultCountForDiagnostics() != 0 ||
	service.discardedStaleResultCountForDiagnostics() == 0) {
	printf("FAIL: LoD service published or retained a cancelled generation result\n");
	return 1;
    }

    if (blockingData.calls != 1) {
	printf("FAIL: LoD service cancellation callback count mismatch\n");
	return 1;
    }

    service.stop();
    return 0;
}

static int
test_destruction_waits_for_in_flight_task(void)
{
    BObolLodService *service = new BObolLodService;
    BlockingTaskData blockingData;
    std::atomic<bool> destroyStarted(false);
    std::atomic<bool> destroyReturned(false);

    if (!service->start(1, FALSE)) {
	printf("FAIL: LoD service did not start for in-flight destruction test\n");
	delete service;
	return 1;
    }

    BObolLodTask task;
    task.generation = service->beginGeneration();
    task.request = make_request("/destroy-in-flight.bot");
    task.realize = blocking_task;
    task.realizeData = &blockingData;
    if (service->submit(task) == 0) {
	printf("FAIL: LoD service rejected in-flight destruction task\n");
	delete service;
	return 1;
    }

    if (wait_for_blocking_task_started(blockingData)) {
	{
	    std::lock_guard<std::mutex> lock(blockingData.mutex);
	    blockingData.release = true;
	    blockingData.cv.notify_all();
	}
	delete service;
	return 1;
    }

    std::thread destroyer([&]() {
	destroyStarted.store(true, std::memory_order_release);
	delete service;
	destroyReturned.store(true, std::memory_order_release);
    });
    while (!destroyStarted.load(std::memory_order_acquire))
	std::this_thread::yield();
    std::this_thread::sleep_for(std::chrono::milliseconds(25));

    int ret = 0;
    if (destroyReturned.load(std::memory_order_acquire)) {
	printf("FAIL: LoD service destruction returned while a worker callback was active\n");
	ret = 1;
    }

    {
	std::lock_guard<std::mutex> lock(blockingData.mutex);
	blockingData.release = true;
	blockingData.cv.notify_all();
    }
    destroyer.join();

    if (!destroyReturned.load(std::memory_order_acquire) ||
	blockingData.calls != 1) {
	printf("FAIL: LoD service destruction did not drain the active worker callback\n");
	ret = 1;
    }

    return ret;
}

static int
test_stop_discards_undrained_state(void)
{
    BObolLodService service;
    ServiceTestContext context;
    TaskData data;
    data.context = &context;
    data.value = 19;

    if (!service.start(1, FALSE)) {
	printf("FAIL: LoD service did not start for teardown test\n");
	return 1;
    }

    BObolLodTask task;
    task.generation = service.beginGeneration();
    task.request = make_request("/teardown.bot");
    task.realize = ready_task;
    task.realizeData = &data;
    if (service.submit(task) == 0) {
	printf("FAIL: LoD service rejected teardown task\n");
	service.stop();
	return 1;
    }

    for (int i = 0; i < 400 &&
	 service.queuedResultCountForDiagnostics() == 0; i++)
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    if (service.queuedResultCountForDiagnostics() != 1) {
	printf("FAIL: LoD teardown task did not publish a result\n");
	service.stop();
	return 1;
    }

    service.stop();
    if (service.isRunning() || service.inFlightCount() != 0 ||
	service.pendingTaskCountForDiagnostics() != 0 ||
	service.queuedResultCountForDiagnostics() != 0 ||
	service.queuedCacheWriteCountForDiagnostics() != 0 ||
	service.activeRequestCountForDiagnostics() != 0 ||
	service.completedTaskCountForDiagnostics() != 0 ||
	service.cancelledGenerationCountForDiagnostics() != 0 ||
	service.currentGeneration() != 0) {
	printf("FAIL: LoD service retained state after stop\n");
	return 1;
    }

    return 0;
}

static BObolLodResult
superseded_task(const BObolLodRequest &request, void *UNUSED(userData))
{
    BObolLodRequest supersededRequest = request;
    supersededRequest.viewRevision++;

    BObolLodResult result;
    result.request = supersededRequest;
    result.cacheKey = bobol_lod_cache_key(supersededRequest);
    result.resultKind = BOBOL_LOD_RESULT_MESH;
    result.qualityTier = supersededRequest.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.counts.faceCount = 11;
    return result;
}

static int
test_superseded_result_rejection(void)
{
    BObolLodService service;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for superseded-result test\n");
	return 1;
    }

    BObolLodTask task;
    task.generation = service.beginGeneration();
    task.request = make_request("/stale.bot");
    task.realize = superseded_task;

    if (service.submit(task) == 0) {
	printf("FAIL: LoD service did not accept superseded-result task\n");
	return 1;
    }

    if (wait_for_settled(service, 1))
	return 1;

    std::vector<BObolLodResult> results;
    service.drainResults(results);
    if (results.size() != 1 ||
	results[0].providerStatus != BOBOL_LOD_PROVIDER_SUPERSEDED ||
	!results[0].stale ||
	!results[0].terminal ||
	results[0].diagnostic.getLength() == 0) {
	printf("FAIL: LoD service did not reject superseded result\n");
	return 1;
    }

    service.stop();
    return 0;
}

static BObolLodResult
staged_payload_task(const BObolLodRequest &request, void *UNUSED(userData))
{
    std::vector<BObolLodAttribute> attributes;
    BObolLodAttribute attribute;

    attribute.name = "display.intent";
    attribute.value = "proxy";
    attributes.push_back(attribute);

    return bobol_lod_attributes_result(request, attributes);
}

static int
test_staged_payload_delivery(void)
{
    BObolLodService service;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for staged-payload test\n");
	return 1;
    }

    BObolLodTask task;
    task.generation = service.beginGeneration();
    task.request = make_request("/attributes.bot");
    task.realize = staged_payload_task;

    if (service.submit(task) == 0) {
	printf("FAIL: LoD service did not accept staged-payload task\n");
	return 1;
    }

    if (wait_for_settled(service, 1))
	return 1;

    std::vector<BObolLodResult> results;
    service.drainResults(results);
    if (results.size() != 1 ||
	results[0].resultKind != BOBOL_LOD_RESULT_ATTRIBUTES ||
	results[0].qualityTier != BOBOL_LOD_QUALITY_ATTRIBUTES ||
	results[0].attributes.size() != 1 ||
	bu_strcmp(results[0].attributes[0].value.getString(), "proxy") != 0 ||
	!bobol_lod_result_matches_request(results[0], task.request)) {
	printf("FAIL: LoD service did not deliver staged payload result\n");
	return 1;
    }

    service.stop();
    return 0;
}

struct ResultReadyContext {
    std::mutex mutex;
    std::condition_variable cv;
    int wakeups;
    size_t maxQueuedCount;

    ResultReadyContext(void) : wakeups(0), maxQueuedCount(0)
    {
    }
};

static void
result_ready_cb(BObolLodService *service, void *userData)
{
    ResultReadyContext *context =
	static_cast<ResultReadyContext *>(userData);
    size_t queued = service ? service->queuedResultCountForDiagnostics() : 0;

    {
	std::lock_guard<std::mutex> lock(context->mutex);
	context->wakeups++;
	if (queued > context->maxQueuedCount)
	    context->maxQueuedCount = queued;
    }
    context->cv.notify_all();
}

static int
wait_for_wakeups(ResultReadyContext &context, int expected)
{
    std::unique_lock<std::mutex> lock(context.mutex);
    if (!context.cv.wait_for(lock, std::chrono::seconds(2),
			     [&context, expected] { return context.wakeups >= expected; })) {
	printf("FAIL: LoD result-ready wakeup count did not reach %d\n",
	       expected);
	return 1;
    }
    return 0;
}

static int
test_result_ready_subscription(void)
{
    BObolLodService service;
    ServiceTestContext serviceContext;
    ResultReadyContext readyContext;
    TaskData firstData{&serviceContext, 101};
    TaskData secondData{&serviceContext, 102};
    TaskData thirdData{&serviceContext, 103};
    TaskData fourthData{&serviceContext, 104};

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for result-ready test\n");
	return 1;
    }

    BObolLodSubscriberId subscriber =
	service.subscribeResultReady(result_ready_cb, &readyContext);
    if (subscriber == 0) {
	printf("FAIL: LoD service did not create result-ready subscription\n");
	service.stop();
	return 1;
    }

    uint64_t generation = service.beginGeneration();
    BObolLodTask first;
    first.generation = generation;
    first.request = make_request("/ready-first.bot");
    first.realize = ready_task;
    first.realizeData = &firstData;
    BObolLodTask second;
    second.generation = generation;
    second.request = make_request("/ready-second.bot");
    second.realize = ready_task;
    second.realizeData = &secondData;

    if (service.submit(first) == 0 || service.submit(second) == 0) {
	printf("FAIL: LoD service did not accept result-ready burst tasks\n");
	service.unsubscribeResultReady(subscriber);
	service.stop();
	return 1;
    }

    if (wait_for_settled(service, 2) || wait_for_wakeups(readyContext, 1)) {
	service.unsubscribeResultReady(subscriber);
	service.stop();
	return 1;
    }
    {
	std::lock_guard<std::mutex> lock(readyContext.mutex);
	if (readyContext.wakeups != 1 || readyContext.maxQueuedCount == 0) {
	    printf("FAIL: LoD result-ready wakeup did not coalesce burst\n");
	    service.unsubscribeResultReady(subscriber);
	    service.stop();
	    return 1;
	}
    }

    std::vector<BObolLodResult> results;
    if (service.drainResults(results) != 2 ||
	service.queuedResultCountForDiagnostics() != 0) {
	printf("FAIL: LoD result-ready test did not drain burst results\n");
	service.unsubscribeResultReady(subscriber);
	service.stop();
	return 1;
    }

    BObolLodTask third;
    third.generation = generation;
    third.request = make_request("/ready-third.bot");
    third.realize = ready_task;
    third.realizeData = &thirdData;
    if (service.submit(third) == 0 ||
	wait_for_settled(service, 1) ||
	wait_for_wakeups(readyContext, 2)) {
	printf("FAIL: LoD result-ready wakeup did not refire after drain\n");
	service.unsubscribeResultReady(subscriber);
	service.stop();
	return 1;
    }

    results.clear();
    if (service.drainResults(results) != 1) {
	printf("FAIL: LoD result-ready test did not drain third result\n");
	service.unsubscribeResultReady(subscriber);
	service.stop();
	return 1;
    }

    service.unsubscribeResultReady(subscriber);
    int wakeupsAfterUnsubscribe = 0;
    {
	std::lock_guard<std::mutex> lock(readyContext.mutex);
	wakeupsAfterUnsubscribe = readyContext.wakeups;
    }

    BObolLodTask fourth;
    fourth.generation = generation;
    fourth.request = make_request("/ready-fourth.bot");
    fourth.realize = ready_task;
    fourth.realizeData = &fourthData;
    if (service.submit(fourth) == 0 || wait_for_settled(service, 1)) {
	printf("FAIL: LoD result-ready test did not run post-unsubscribe task\n");
	service.stop();
	return 1;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(25));
    {
	std::lock_guard<std::mutex> lock(readyContext.mutex);
	if (readyContext.wakeups != wakeupsAfterUnsubscribe) {
	    printf("FAIL: LoD result-ready unsubscribe did not suppress callbacks\n");
	    service.stop();
	    return 1;
	}
    }

    service.stop();
    return 0;
}

struct SelfUnsubscribeContext {
    std::mutex mutex;
    std::condition_variable cv;
    BObolLodSubscriberId subscriber;
    int calls;

    SelfUnsubscribeContext(void) : subscriber(0), calls(0)
    {
    }
};

static void
self_unsubscribe_result_ready(BObolLodService *service, void *userData)
{
    SelfUnsubscribeContext *context =
	static_cast<SelfUnsubscribeContext *>(userData);
    if (!context || !service)
	return;

    {
	std::lock_guard<std::mutex> lock(context->mutex);
	context->calls++;
	context->cv.notify_all();
    }
    service->unsubscribeResultReady(context->subscriber);
}

static int
test_result_ready_self_unsubscribe(void)
{
    BObolLodService service;
    ServiceTestContext serviceContext;
    SelfUnsubscribeContext callbackContext;
    TaskData firstData{&serviceContext, 121};
    TaskData secondData{&serviceContext, 122};

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for self-unsubscribe test\n");
	return 1;
    }

    callbackContext.subscriber = service.subscribeResultReady(
	self_unsubscribe_result_ready, &callbackContext);
    if (callbackContext.subscriber == 0) {
	printf("FAIL: LoD service did not create self-unsubscribe subscription\n");
	service.stop();
	return 1;
    }

    const uint64_t generation = service.beginGeneration();
    BObolLodTask first;
    first.generation = generation;
    first.request = make_request("/ready-self-unsubscribe-first.bot");
    first.realize = ready_task;
    first.realizeData = &firstData;
    if (service.submit(first) == 0 || wait_for_settled(service, 1)) {
	printf("FAIL: LoD service did not deliver self-unsubscribe result\n");
	service.stop();
	return 1;
    }

    {
	std::unique_lock<std::mutex> lock(callbackContext.mutex);
	if (!callbackContext.cv.wait_for(lock, std::chrono::seconds(2),
				 [&callbackContext] {
				     return callbackContext.calls == 1;
				 })) {
	    printf("FAIL: LoD self-unsubscribe callback did not return\n");
	    service.stop();
	    return 1;
	}
    }

    std::vector<BObolLodResult> results;
    service.drainResults(results);
    BObolLodTask second;
    second.generation = generation;
    second.request = make_request("/ready-self-unsubscribe-second.bot");
    second.realize = ready_task;
    second.realizeData = &secondData;
    if (service.submit(second) == 0 || wait_for_settled(service, 1)) {
	printf("FAIL: LoD service did not run after self-unsubscribe\n");
	service.stop();
	return 1;
    }
    {
	std::lock_guard<std::mutex> lock(callbackContext.mutex);
	if (callbackContext.calls != 1) {
	    printf("FAIL: LoD self-unsubscribe callback was invoked again\n");
	    service.stop();
	    return 1;
	}
    }

    service.stop();
    return 0;
}

struct CrossUnsubscribeContext {
    std::mutex mutex;
    std::condition_variable cv;
    BObolLodSubscriberId targetSubscriber;
    int removerCalls;
    int targetCalls;
    int observerCalls;

    CrossUnsubscribeContext(void) : targetSubscriber(0), removerCalls(0),
	targetCalls(0), observerCalls(0)
    {
    }
};

static void
cross_unsubscribe_remover(BObolLodService *service, void *userData)
{
    CrossUnsubscribeContext *context =
	static_cast<CrossUnsubscribeContext *>(userData);
    if (!context || !service)
	return;

    service->unsubscribeResultReady(context->targetSubscriber);
    {
	std::lock_guard<std::mutex> lock(context->mutex);
	context->removerCalls++;
	context->cv.notify_all();
    }
}

static void
cross_unsubscribe_target(BObolLodService *UNUSED(service), void *userData)
{
    CrossUnsubscribeContext *context =
	static_cast<CrossUnsubscribeContext *>(userData);
    if (!context)
	return;

    std::lock_guard<std::mutex> lock(context->mutex);
    context->targetCalls++;
    context->cv.notify_all();
}

static void
cross_unsubscribe_observer(BObolLodService *UNUSED(service), void *userData)
{
    CrossUnsubscribeContext *context =
	static_cast<CrossUnsubscribeContext *>(userData);
    if (!context)
	return;

    std::lock_guard<std::mutex> lock(context->mutex);
    context->observerCalls++;
    context->cv.notify_all();
}

static int
test_result_ready_cross_unsubscribe(void)
{
    BObolLodService service;
    ServiceTestContext serviceContext;
    CrossUnsubscribeContext callbackContext;
    TaskData firstData{&serviceContext, 123};
    TaskData secondData{&serviceContext, 124};

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for cross-unsubscribe test\n");
	return 1;
    }

    if (service.subscribeResultReady(cross_unsubscribe_remover,
				     &callbackContext) == 0 ||
	(callbackContext.targetSubscriber = service.subscribeResultReady(
	    cross_unsubscribe_target, &callbackContext)) == 0 ||
	service.subscribeResultReady(cross_unsubscribe_observer,
				     &callbackContext) == 0) {
	printf("FAIL: LoD service did not create cross-unsubscribe subscriptions\n");
	service.stop();
	return 1;
    }

    const uint64_t generation = service.beginGeneration();
    BObolLodTask first;
    first.generation = generation;
    first.request = make_request("/ready-cross-unsubscribe-first.bot");
    first.realize = ready_task;
    first.realizeData = &firstData;
    if (service.submit(first) == 0 || wait_for_settled(service, 1)) {
	printf("FAIL: LoD service did not deliver cross-unsubscribe result\n");
	service.stop();
	return 1;
    }

    {
	std::lock_guard<std::mutex> lock(callbackContext.mutex);
	if (callbackContext.removerCalls != 1 || callbackContext.targetCalls != 0 ||
	    callbackContext.observerCalls != 1) {
	    printf("FAIL: LoD cross-unsubscribe did not suppress the reserved callback\n");
	    service.stop();
	    return 1;
	}
    }

    std::vector<BObolLodResult> results;
    service.drainResults(results);
    BObolLodTask second;
    second.generation = generation;
    second.request = make_request("/ready-cross-unsubscribe-second.bot");
    second.realize = ready_task;
    second.realizeData = &secondData;
    if (service.submit(second) == 0 || wait_for_settled(service, 1)) {
	printf("FAIL: LoD service did not continue after cross-unsubscribe\n");
	service.stop();
	return 1;
    }
    {
	std::lock_guard<std::mutex> lock(callbackContext.mutex);
	if (callbackContext.removerCalls != 2 || callbackContext.targetCalls != 0 ||
	    callbackContext.observerCalls != 2) {
	    printf("FAIL: LoD cross-unsubscribe changed surviving callback delivery\n");
	    service.stop();
	    return 1;
	}
    }

    service.stop();
    return 0;
}

struct CleanupMarker {
    std::atomic<int> *calls;
    std::atomic<int> *cleanups;
};

static BObolLodResult
cleanup_task(const BObolLodRequest &request, void *userData)
{
    CleanupMarker *marker = static_cast<CleanupMarker *>(userData);
    marker->calls->fetch_add(1);

    BObolLodResult result;
    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = BOBOL_LOD_RESULT_AABB;
    result.qualityTier = request.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.bounds = request.bounds;
    return result;
}

static void
cleanup_marker_free(void *userData)
{
    CleanupMarker *marker = static_cast<CleanupMarker *>(userData);
    marker->cleanups->fetch_add(1);
    delete marker;
}

static int
test_task_realize_data_cleanup(void)
{
    std::atomic<int> calls(0);
    std::atomic<int> cleanups(0);
    BObolLodService service;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for cleanup test\n");
	return 1;
    }

    BObolLodTask task;
    task.generation = service.beginGeneration();
    task.request = make_request("/cleanup.bot");
    task.realize = cleanup_task;
    task.realizeData = new CleanupMarker{&calls, &cleanups};
    task.realizeDataFree = cleanup_marker_free;

    if (service.submit(task) == 0) {
	printf("FAIL: LoD service did not accept cleanup task\n");
	cleanup_marker_free(task.realizeData);
	service.stop();
	return 1;
    }

    if (wait_for_settled(service, 1)) {
	service.stop();
	return 1;
    }

    std::vector<BObolLodResult> results;
    service.drainResults(results);
    service.stop();
    if (results.size() != 1 ||
	calls.load() != 1 ||
	cleanups.load() != 1) {
	printf("FAIL: LoD service did not clean completed task data\n");
	return 1;
    }

    std::atomic<int> pendingCalls(0);
    std::atomic<int> pendingCleanups(0);
    BObolLodService pendingService;
    if (!pendingService.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for pending-cleanup test\n");
	return 1;
    }

    BObolLodTask pendingTask;
    pendingTask.generation = pendingService.beginGeneration();
    pendingTask.request = make_request("/pending-cleanup.bot");
    pendingTask.realize = cleanup_task;
    pendingTask.realizeData =
	new CleanupMarker{&pendingCalls, &pendingCleanups};
    pendingTask.realizeDataFree = cleanup_marker_free;
    pendingTask.addDependency(999999);

    if (pendingService.submit(pendingTask) == 0) {
	printf("FAIL: LoD service did not accept pending-cleanup task\n");
	cleanup_marker_free(pendingTask.realizeData);
	pendingService.stop();
	return 1;
    }

    for (int i = 0; i < 100 &&
	 pendingService.pendingTaskCountForDiagnostics() != 1; i++)
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    pendingService.stop();
    if (pendingCalls.load() != 0 ||
	pendingCleanups.load() != 1 ||
	pendingService.inFlightCount() != 0) {
	printf("FAIL: LoD service did not clean pending task data on stop\n");
	return 1;
    }

    return 0;
}

static int
test_queue_limits_and_pending_cancellation(void)
{
    BObolLodService service;
    service.setQueueLimits(2, 1, 1);
    size_t maxActive = 0;
    size_t maxResults = 0;
    size_t maxCacheWrites = 0;
    service.getQueueLimits(maxActive, maxResults, maxCacheWrites);
    if (maxActive != 2 || maxResults != 1 || maxCacheWrites != 1) {
	printf("FAIL: LoD service queue limits did not round trip\n");
	return 1;
    }
    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for queue-limit test\n");
	return 1;
    }

    std::atomic<int> calls(0);
    std::atomic<int> cleanups(0);
    const uint64_t generation = service.beginGeneration();
    BObolLodTask first;
    first.generation = generation;
    first.request = make_request("/bounded-first.bot");
    first.realize = cleanup_task;
    first.realizeData = new CleanupMarker{&calls, &cleanups};
    first.realizeDataFree = cleanup_marker_free;
    first.publishResult = TRUE;
    first.addDependency(999998);
    if (service.submit(first) == 0) {
	cleanup_marker_free(first.realizeData);
	service.stop();
	printf("FAIL: LoD service rejected first bounded task\n");
	return 1;
    }

    BObolLodTask rejected = first;
    rejected.request = make_request("/bounded-rejected.bot");
    rejected.realizeData = new CleanupMarker{&calls, &cleanups};
    if (service.submit(rejected) != 0 ||
	service.rejectedTaskCountForDiagnostics() != 1) {
	cleanup_marker_free(rejected.realizeData);
	service.stop();
	printf("FAIL: LoD service did not enforce reserved result capacity\n");
	return 1;
    }
    cleanup_marker_free(rejected.realizeData);

    service.cancelGeneration(generation);
    if (!service.isGenerationCancelled(generation) ||
	service.currentGeneration() != 0 ||
	service.pendingTaskCountForDiagnostics() != 0 ||
	service.inFlightCount() != 0 || calls.load() != 0 ||
	cleanups.load() != 2) {
	service.stop();
	printf("FAIL: LoD generation cancellation did not immediately free pending work\n");
	return 1;
    }
    service.stop();
    return 0;
}

static int
test_large_pending_cancellation_and_generation_history(void)
{
    static const int requestCount = 10000;
    BObolLodService service;
    service.setQueueLimits(requestCount + 16, 1, 1);
    if (!service.start(1, FALSE)) {
	printf("FAIL: LoD scale service did not start\n");
	return 1;
    }

    const uint64_t generation = service.beginGeneration();
    for (int i = 0; i < requestCount; i++) {
	char name[96] = {0};
	snprintf(name, sizeof(name), "/scale/%d.bot", i);
	BObolLodTask task;
	task.generation = generation;
	task.request = make_request(name);
	task.realize = ready_task;
	task.publishResult = FALSE;
	task.addDependency(UINT64_MAX - 1);
	if (service.submitIfNotActive(task) == 0) {
	    printf("FAIL: LoD scale service rejected request %d\n", i);
	    service.stop();
	    return 1;
	}
    }
    if (service.pendingTaskCountForDiagnostics() != requestCount ||
	service.activeRequestCountForDiagnostics() != requestCount ||
	service.inFlightCount() != requestCount) {
	printf("FAIL: LoD scale service did not retain bounded pending state\n");
	service.stop();
	return 1;
    }

    service.cancelGeneration(generation);
    if (service.pendingTaskCountForDiagnostics() != 0 ||
	service.activeRequestCountForDiagnostics() != 0 ||
	service.inFlightCount() != 0 ||
	service.completedTaskCountForDiagnostics() != 0) {
	printf("FAIL: LoD scale cancellation did not release request bookkeeping\n");
	service.stop();
	return 1;
    }

    uint64_t newestGeneration = 0;
    for (int i = 0; i < 2048; i++) {
	newestGeneration = service.beginGeneration();
	service.cancelGeneration(newestGeneration);
    }
    if (!service.isGenerationCancelled(newestGeneration) ||
	service.cancelledGenerationCountForDiagnostics() > 1024) {
	printf("FAIL: LoD cancellation history is not bounded\n");
	service.stop();
	return 1;
    }

    service.stop();
    return 0;
}

struct DebugDelayTaskData {
    std::mutex mutex;
    int calls;
    SbBool overrideResultRequest;
    BObolLodRequest resultRequest;

    DebugDelayTaskData(void) : calls(0), overrideResultRequest(FALSE)
    {
    }
};

static BObolLodResult
debug_delay_task(const BObolLodRequest &request, void *userData)
{
    DebugDelayTaskData *data = static_cast<DebugDelayTaskData *>(userData);
    BObolLodRequest resultRequest = request;

    {
	std::lock_guard<std::mutex> lock(data->mutex);
	data->calls++;
	if (data->overrideResultRequest)
	    resultRequest = data->resultRequest;
    }

    BObolLodResult result;
    result.request = resultRequest;
    result.cacheKey = bobol_lod_cache_key(resultRequest);
    result.resultKind = BOBOL_LOD_RESULT_AABB;
    result.qualityTier = resultRequest.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.terminal = TRUE;
    return result;
}

static int
wait_for_debug_delay(BObolLodService &service)
{
    for (int i = 0; i < 200; i++) {
	if (service.delayedTaskCountForDiagnostics() == 1)
	    return 0;
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    printf("FAIL: LoD service did not report delayed task: delayed=%zu inFlight=%zu\n",
	   service.delayedTaskCountForDiagnostics(),
	   service.inFlightCount());
    return 1;
}

static int
test_debug_delay_cancellation(void)
{
    BObolLodService service;
    DebugDelayTaskData data;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for debug-delay test\n");
	return 1;
    }

    uint64_t generation = service.beginGeneration();
    BObolLodTask task;
    task.generation = generation;
    task.request = make_request("/debug-delay.bot");
    task.realize = debug_delay_task;
    task.realizeData = &data;
    task.debugDelayMilliseconds = 200;

    if (service.submit(task) == 0) {
	printf("FAIL: LoD service did not accept debug-delay task\n");
	return 1;
    }

    if (wait_for_debug_delay(service))
	return 1;

    service.cancelGeneration(generation);


    for (int i = 0; i < 400 && service.inFlightCount() != 0; i++)
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    if (service.inFlightCount() != 0 ||
	service.queuedResultCountForDiagnostics() != 0 ||
	service.discardedStaleResultCountForDiagnostics() == 0 ||
	service.delayedTaskCountForDiagnostics() != 0) {
	printf("FAIL: LoD service did not cancel delayed task cleanly\n");
	return 1;
    }

    {
	std::lock_guard<std::mutex> lock(data.mutex);
	if (data.calls != 0) {
	    printf("FAIL: LoD service invoked callback during cancelled debug delay\n");
	    return 1;
	}
    }

    service.stop();
    return 0;
}

struct BlockingCacheWriteData {
    std::mutex mutex;
    std::condition_variable cv;
    bool started;
    bool release;
    int calls;
    std::vector<int> values;

    BlockingCacheWriteData(void) : started(false), release(false), calls(0)
    {
    }
};

static void
blocking_cache_write(const BObolLodResult &result, void *userData)
{
    BlockingCacheWriteData *data =
	static_cast<BlockingCacheWriteData *>(userData);

    std::unique_lock<std::mutex> lock(data->mutex);
    data->started = true;
    data->calls++;
    data->values.push_back((int)result.counts.faceCount);
    data->cv.notify_all();
    while (!data->release)
	data->cv.wait(lock);
}

static int
wait_for_blocking_cache_write_started(BlockingCacheWriteData &data)
{
    std::unique_lock<std::mutex> lock(data.mutex);
    if (!data.cv.wait_for(lock, std::chrono::seconds(2),
			  [&data] { return data.started; })) {
	printf("FAIL: LoD blocking cache write did not start\n");
	return 1;
    }
    return 0;
}

static int
wait_for_cache_write_count(BObolLodService &service, size_t expected)
{
    for (int i = 0; i < 400; i++) {
	if (service.queuedCacheWriteCountForDiagnostics() >= expected)
	    return 0;
	std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    printf("FAIL: LoD cache-write queue did not reach %zu: cache=%zu\n",
	   expected, service.queuedCacheWriteCountForDiagnostics());
    return 1;
}

static int
test_cancelled_cache_write_not_persisted(void)
{
    BObolLodService service;
    ServiceTestContext context;
    BlockingCacheWriteData blockingWrite;
    TaskData firstData{&context, 401};
    TaskData cancelledData{&context, 402};
    TaskData finalData{&context, 403};

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for cache cancellation test\n");
	return 1;
    }

    const uint64_t cancelledGeneration = service.beginGeneration();
    BObolLodTask first;
    first.generation = cancelledGeneration;
    first.request = make_request("/cache-first.bot");
    first.realize = ready_task;
    first.realizeData = &firstData;
    first.publishResult = FALSE;
    first.writeCache = TRUE;
    first.cacheWrite = blocking_cache_write;
    first.cacheWriteData = &blockingWrite;
    if (service.submit(first) == 0 ||
	wait_for_blocking_cache_write_started(blockingWrite)) {
	printf("FAIL: LoD service did not start the cache-write blocker\n");
	service.stop();
	return 1;
    }

    BObolLodTask cancelled = first;
    cancelled.request = make_request("/cache-cancelled.bot");
    cancelled.realizeData = &cancelledData;
    cancelled.cacheWrite = cache_write;
    cancelled.cacheWriteData = &context;
    if (service.submit(cancelled) == 0 ||
	wait_for_cache_write_count(service, 2)) {
	printf("FAIL: LoD service did not queue the cancellable cache write\n");
	{
	    std::lock_guard<std::mutex> lock(blockingWrite.mutex);
	    blockingWrite.release = true;
	    blockingWrite.cv.notify_all();
	}
	service.stop();
	return 1;
    }

    service.cancelGeneration(cancelledGeneration);
    if (!service.isGenerationCancelled(cancelledGeneration) ||
	service.queuedCacheWriteCountForDiagnostics() != 1) {
	printf("FAIL: LoD generation cancellation retained queued cache work\n");
	{
	    std::lock_guard<std::mutex> lock(blockingWrite.mutex);
	    blockingWrite.release = true;
	    blockingWrite.cv.notify_all();
	}
	service.stop();
	return 1;
    }

    {
	std::lock_guard<std::mutex> lock(blockingWrite.mutex);
	blockingWrite.release = true;
	blockingWrite.cv.notify_all();
    }
    if (wait_for_settled(service, 0)) {
	service.stop();
	return 1;
    }
    {
	std::lock_guard<std::mutex> lock(blockingWrite.mutex);
	if (blockingWrite.calls != 1 || blockingWrite.values.size() != 1 ||
	    blockingWrite.values[0] != 401) {
	    printf("FAIL: LoD cancellation persisted queued stale cache work\n");
	    service.stop();
	    return 1;
	}
    }
    {
	std::lock_guard<std::mutex> lock(context.mutex);
	if (!context.cacheWriteOrder.empty()) {
	    printf("FAIL: LoD cancellation invoked the queued stale cache callback\n");
	    service.stop();
	    return 1;
	}
    }

    const uint64_t liveGeneration = service.beginGeneration();
    BObolLodTask finalTask;
    finalTask.generation = liveGeneration;
    finalTask.request = make_request("/cache-live.bot");
    finalTask.realize = ready_task;
    finalTask.realizeData = &finalData;
    finalTask.publishResult = FALSE;
    finalTask.writeCache = TRUE;
    finalTask.cacheWrite = cache_write;
    finalTask.cacheWriteData = &context;
    if (service.submit(finalTask) == 0 || wait_for_settled(service, 0)) {
	printf("FAIL: LoD service did not resume cache writes after cancellation\n");
	service.stop();
	return 1;
    }
    {
	std::lock_guard<std::mutex> lock(context.mutex);
	if (context.cacheWriteOrder.size() != 1 ||
	    context.cacheWriteOrder[0] != 403) {
	    printf("FAIL: LoD service cache writer did not isolate live generation\n");
	    service.stop();
	    return 1;
	}
    }

    service.stop();
    return 0;
}

static int
make_provider_test_db_version(char *dbpath, size_t dbpath_len,
	struct db_i **dbip_out, int databaseVersion)
{
    static const char *objname = "lod-provider.bot";
    static const char *tri_objname = "lod-two-tri.bot";
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
    fastf_t tri_vertices[18] = {
	0.0, 0.0, 0.0,
	1.0, 0.0, 0.0,
	0.0, 1.0, 0.0,
	10.0, 0.0, 0.0,
	11.0, 0.0, 0.0,
	10.0, 1.0, 0.0
    };
    int tri_faces[6] = {
	0, 1, 2,
	3, 4, 5
    };

    if (!dbpath || dbpath_len == 0 || !dbip_out)
	return 1;
    *dbip_out = NULL;

    FILE *fp = bu_temp_file(dbpath, dbpath_len);
    if (!fp) {
	printf("FAIL: LoD Obol mesh provider temp file\n");
	return 1;
    }
    fclose(fp);

    struct db_i *dbip = db_create(dbpath, databaseVersion);
    if (!dbip) {
	printf("FAIL: LoD Obol mesh provider db_create\n");
	bu_file_delete(dbpath);
	return 1;
    }

    struct rt_wdb *wdbp = wdb_dbopen(dbip, RT_WDB_TYPE_DB_DISK);
    if (!wdbp) {
	printf("FAIL: LoD Obol mesh provider wdb_dbopen\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    if (mk_bot(wdbp, objname, RT_BOT_SOLID, RT_BOT_UNORIENTED, 0,
	       4, 4, vertices, faces, NULL, NULL) != 0) {
	printf("FAIL: LoD Obol mesh provider mk_bot\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    if (mk_bot(wdbp, tri_objname, RT_BOT_SOLID, RT_BOT_UNORIENTED, 0,
	       6, 2, tri_vertices, tri_faces, NULL, NULL) != 0) {
	printf("FAIL: LoD Obol mesh provider mk_bot disjoint triangles\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    struct wmember snapshot_members;
    BU_LIST_INIT(&snapshot_members.l);
    if (!mk_addmember(objname, &snapshot_members.l, NULL, WMOP_UNION) ||
	mk_lcomb(wdbp, "lod-snapshot.c", &snapshot_members, 0, NULL, NULL,
	    NULL, 0) != 0) {
	printf("FAIL: LoD detached snapshot combination\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    *dbip_out = dbip;
    return 0;
}

static int
make_provider_test_db(char *dbpath, size_t dbpath_len, struct db_i **dbip_out)
{
    return make_provider_test_db_version(dbpath, dbpath_len, dbip_out, 5);
}

static int
test_database_lease_outlives_submitter(void)
{
    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;
    if (make_provider_test_db(dbpath, sizeof(dbpath), &dbip))
	return 1;

    BObolLodService service;
    if (!service.start(1, FALSE)) {
	db_close(dbip);
	bu_file_delete(dbpath);
	printf("FAIL: database-lease service setup\n");
	return 1;
    }

    BObolRtSourceFullDetailProvider *provider =
	new BObolRtSourceFullDetailProvider;
    if (!provider->setDatabase(dbip)) {
	delete provider;
	service.stop();
	db_close(dbip);
	bu_file_delete(dbpath);
	printf("FAIL: database lease acquisition\n");
	return 1;
    }

    BObolLodTask task;
    task.generation = service.beginGeneration();
    task.request = make_request("/lod-provider.bot");
    task.request.objectName = "lod-provider.bot";
    task.request.sourceCounts.faceCount = 4;
    task.request.sourceCounts.pointCount = 4;
    task.realize = bobol_rt_source_full_detail_provider_task;
    task.realizeData = provider;
    task.realizeDataFree = bobol_rt_source_full_detail_provider_free;
    task.debugDelayMilliseconds = 100;
    if (!service.submit(task)) {
	delete provider;
	service.stop();
	db_close(dbip);
	bu_file_delete(dbpath);
	printf("FAIL: database-lease task submission\n");
	return 1;
    }

    /* Release the application's only database use before the worker wakes. */
    db_close(dbip);
    dbip = NULL;

    if (wait_for_settled(service, 1)) {
	service.stop();
	bu_file_delete(dbpath);
	return 1;
    }
    std::vector<BObolLodResult> results;
    service.drainResults(results);
    service.stop();
    bu_file_delete(dbpath);
    if (results.size() != 1 ||
	results.front().providerStatus != BOBOL_LOD_PROVIDER_READY ||
	results.front().resultKind != BOBOL_LOD_RESULT_FULL_DETAIL ||
	!results.front().mesh.isValid()) {
	printf("FAIL: background provider lost its database lease\n");
	return 1;
    }
    return 0;
}

static int
check_source_full_detail_result(const BObolLodResult &result,
				const BObolLodRequest &request, const char *label)
{
    if (result.resultKind != BOBOL_LOD_RESULT_FULL_DETAIL ||
	result.qualityTier != BOBOL_LOD_QUALITY_FULL_DETAIL ||
	result.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	result.geometry.kind != BOBOL_LOD_GEOMETRY_OBOL_MESH ||
	!result.geometry.isValid() ||
	result.counts.faceCount != 4 ||
	result.counts.pointCount != 4 ||
	!result.mesh.isValid() ||
	result.mesh.points.size() != 4 ||
	result.mesh.coordIndex.size() != 12 ||
	result.mesh.faceIndex.size() != 4 ||
	result.mesh.faceIndex[0] != 0 ||
	result.mesh.faceIndex[1] != 1 ||
	result.mesh.faceIndex[2] != 2 ||
	result.mesh.faceIndex[3] != 3 ||
	!result.mesh.vertexIndex.empty() ||
	!bobol_lod_result_matches_request(result, request)) {
	printf("FAIL: LoD RT source full-detail provider %s\n", label);
	return 1;
    }

    return 0;
}

static int
test_active_request_duplicate_suppression(void)
{
    BObolLodService service;
    DebugDelayTaskData data;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for active-request test\n");
	return 1;
    }

    BObolLodTask task;
    task.generation = service.beginGeneration();
    task.request = make_request("/active-duplicate.bot");
    task.realize = debug_delay_task;
    task.realizeData = &data;
    task.debugDelayMilliseconds = 80;
    BObolLodRequest movedDemand = task.request;
    movedDemand.viewRevision.advance();
    movedDemand.policyRevision.advance();
    data.overrideResultRequest = TRUE;
    data.resultRequest = movedDemand;

    uint64_t firstId = service.submitIfNotActive(task);
    if (firstId == 0 || !service.hasActiveRequest(task.request)) {
	printf("FAIL: LoD service did not track active request identity\n");
	service.stop();
	return 1;
    }

    uint64_t duplicateId = service.submitIfNotActive(task);
    if (duplicateId != 0) {
	printf("FAIL: LoD service accepted duplicate active request\n");
	service.stop();
	return 1;
    }
    BObolLodTask richerTask = task;
    richerTask.request.requestedCut = task.request.requestedCut + 1;
    if (!service.hasActiveRequest(richerTask.request) ||
	service.submitIfNotActive(richerTask) != 0) {
	printf("FAIL: LoD service queued a superseding level for the same "
	       "serialized resident occurrence\n");
	service.stop();
	return 1;
    }
    if (!service.updateActiveRequestDemand(movedDemand)) {
	printf("FAIL: LoD service did not retain newer demand for active "
	       "asset producer\n");
	service.stop();
	return 1;
    }

    if (wait_for_settled(service, 1)) {
	service.stop();
	return 1;
    }
    if (!service.hasActiveRequest(task.request) ||
	service.submitIfNotActive(task) != 0) {
	printf("FAIL: LoD service released request identity before queued "
	       "result publication\n");
	service.stop();
	return 1;
    }

    std::vector<BObolLodResult> results;
    if (service.drainResults(results) != 1 || results.size() != 1) {
	printf("FAIL: LoD service duplicate suppression changed result count\n");
	service.stop();
	return 1;
    }
    if (!bobol_lod_result_matches_request(results[0], movedDemand) ||
	results[0].providerStatus != BOBOL_LOD_PROVIDER_READY) {
	printf("FAIL: LoD service rejected a result resolved against retained "
	       "latest producer demand\n");
	service.stop();
	return 1;
    }
    if (service.hasActiveRequest(task.request)) {
	printf("FAIL: LoD service retained request identity after result drain\n");
	service.stop();
	return 1;
    }
    if (service.updateActiveRequestDemand(movedDemand)) {
	printf("FAIL: LoD service accepted demand for an inactive producer\n");
	service.stop();
	return 1;
    }
    {
	std::lock_guard<std::mutex> lock(data.mutex);
	data.overrideResultRequest = FALSE;
    }

    uint64_t secondId = service.submitIfNotActive(task);
    if (secondId == 0 || secondId == firstId ||
	!service.hasActiveRequest(task.request)) {
	printf("FAIL: LoD service did not allow completed request resubmission\n");
	service.stop();
	return 1;
    }
    if (wait_for_settled(service, 1)) {
	service.stop();
	return 1;
    }

    results.clear();
    if (service.drainResults(results) != 1 || results.size() != 1) {
	printf("FAIL: LoD service completed-request resubmission did not publish result\n");
	service.stop();
	return 1;
    }

    {
	std::lock_guard<std::mutex> lock(data.mutex);
	if (data.calls != 2) {
	    printf("FAIL: LoD service duplicate active request executed unexpectedly\n");
	    service.stop();
	    return 1;
	}
    }

    /* A result slot may coalesce a newer content revision before the owner
     * drains it.  Duplicate ownership must follow the replacement request;
     * otherwise the old content key remains suppressed forever and the new
     * key appears available while its result is still queued. */
    BObolLodTask queuedOld = task;
    queuedOld.debugDelayMilliseconds = 0;
    if (service.submit(queuedOld) == 0 || wait_for_settled(service, 1)) {
	printf("FAIL: LoD service did not queue old coalescing fixture\n");
	service.stop();
	return 1;
    }
    BObolLodTask queuedNew = queuedOld;
    queuedNew.request.sourceRevision++;
    queuedNew.request.sourceContentHash++;
    if (service.submit(queuedNew) == 0 || wait_for_settled(service, 1) ||
	service.hasActiveRequest(queuedOld.request) ||
	!service.hasActiveRequest(queuedNew.request)) {
	printf("FAIL: LoD service did not transfer queued request identity "
	       "during result coalescing\n");
	service.stop();
	return 1;
    }
    results.clear();
    if (service.drainResults(results) != 1 || results.size() != 1 ||
	!bobol_lod_result_matches_request(results[0], queuedNew.request) ||
	service.hasActiveRequest(queuedOld.request) ||
	service.hasActiveRequest(queuedNew.request)) {
	printf("FAIL: LoD service retained coalesced request ownership after "
	       "publication\n");
	service.stop();
	return 1;
    }

    service.stop();
    return 0;
}

static int
test_asset_producer_coalescing(void)
{
    BObolLodService service;
    DebugDelayTaskData data;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for asset producer test\n");
	return 1;
    }

    BObolLodTask first;
    first.generation = service.beginGeneration();
    first.request = make_request("/shared-asset.bot");
    first.request.occurrenceKey = "/assembly/first/shared-asset.bot";
    first.request.coalesceAssetProducer = TRUE;
    first.realize = debug_delay_task;
    first.realizeData = &data;
    first.debugDelayMilliseconds = 80;

    BObolLodTask sibling = first;
    sibling.request.occurrenceKey = "/assembly/second/shared-asset.bot";
    std::vector<BObolLodTask> tasks{first, sibling};
    std::vector<uint64_t> taskIds;
    if (service.submitBatch(tasks, taskIds, TRUE) != 1 ||
	taskIds.size() != 2 || taskIds[0] == 0 || taskIds[1] != 0 ||
	!service.hasActiveRequest(sibling.request)) {
	printf("FAIL: LoD service did not coalesce sibling asset producers\n");
	service.stop();
	return 1;
    }

    if (wait_for_settled(service, 1)) {
	service.stop();
	return 1;
    }
    std::vector<BObolLodResult> results;
    if (service.drainResults(results) != 1 || results.size() != 1 ||
	!bobol_lod_result_matches_request(results.front(), first.request) ||
	service.hasActiveRequest(sibling.request)) {
	printf("FAIL: LoD service did not release coalesced producer ownership\n");
	service.stop();
	return 1;
    }

    const uint64_t siblingId = service.submitIfNotActive(sibling);
    if (!siblingId || wait_for_settled(service, 1)) {
	printf("FAIL: LoD service did not admit sibling after asset publication\n");
	service.stop();
	return 1;
    }
    results.clear();
    if (service.drainResults(results) != 1 || results.size() != 1 ||
	!bobol_lod_result_matches_request(results.front(), sibling.request)) {
	printf("FAIL: LoD service did not publish the later sibling result\n");
	service.stop();
	return 1;
    }
    {
	std::lock_guard<std::mutex> lock(data.mutex);
	if (data.calls != 2) {
	    printf("FAIL: LoD service executed duplicate sibling producer work\n");
	    service.stop();
	    return 1;
	}
    }

    service.stop();
    return 0;
}

struct SharedProducerLeaseState {
    std::mutex mutex;
    std::condition_variable cv;
    bool started = false;
    bool release = false;
    size_t calls = 0;
};

struct SharedProducerLeaseFixture {
    SharedProducerLeaseState state;
    BObolLodService service;

    ~SharedProducerLeaseFixture(void)
    {
	{
	    std::lock_guard<std::mutex> lock(state.mutex);
	    state.release = true;
	}
	state.cv.notify_all();
	service.stop();
    }
};

static BObolLodResult
shared_producer_lease_task(const BObolLodRequest &request, void *userData)
{
    SharedProducerLeaseState *state =
	static_cast<SharedProducerLeaseState *>(userData);
    {
	std::unique_lock<std::mutex> lock(state->mutex);
	state->started = true;
	++state->calls;
	state->cv.notify_all();
	state->cv.wait(lock, [state] { return state->release; });
    }

    BObolLodResult result;
    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = BOBOL_LOD_RESULT_AABB;
    result.qualityTier = request.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.terminal = TRUE;
    result.counts.faceCount = 17;
    result.bounds = request.bounds;
    return result;
}

static BObolLodTask
shared_producer_lease_task_for(SharedProducerLeaseFixture &fixture,
	uint64_t generation, const char *occurrence)
{
    BObolLodTask task;
    task.generation = generation;
    task.request = make_request("/shared-producer-lease.bot");
    task.request.occurrenceKey = occurrence;
    task.request.coalesceAssetProducer = TRUE;
    task.realize = shared_producer_lease_task;
    task.realizeData = &fixture.state;
    return task;
}

static int
wait_for_shared_producer_started(SharedProducerLeaseFixture &fixture,
	const char *scenario)
{
    std::unique_lock<std::mutex> lock(fixture.state.mutex);
    if (fixture.state.cv.wait_for(lock, std::chrono::seconds(2),
	    [&fixture] { return fixture.state.started; }))
	return 0;
    printf("FAIL: shared producer did not start (%s)\n", scenario);
    return 1;
}

static void
release_shared_producer(SharedProducerLeaseFixture &fixture)
{
    {
	std::lock_guard<std::mutex> lock(fixture.state.mutex);
	fixture.state.release = true;
    }
    fixture.state.cv.notify_all();
}

static int
check_shared_producer_clean(BObolLodService &service,
	const char *scenario)
{
    if (service.inFlightCount() == 0 &&
	service.queuedResultCountForDiagnostics() == 0 &&
	service.activeRequestCountForDiagnostics() == 0 &&
	service.sharedProducerCountForDiagnostics() == 0 &&
	service.sharedProducerLeaseCountForDiagnostics() == 0)
	return 0;
    printf("FAIL: shared producer retained ownership (%s): "
	   "in_flight=%zu results=%zu requests=%zu producers=%zu leases=%zu\n",
	   scenario, service.inFlightCount(),
	   service.queuedResultCountForDiagnostics(),
	   service.activeRequestCountForDiagnostics(),
	   service.sharedProducerCountForDiagnostics(),
	   service.sharedProducerLeaseCountForDiagnostics());
    return 1;
}

static int
check_shared_consumer_result(const std::vector<BObolLodResult> &results,
	uint64_t generation, const BObolLodRequest &request,
	int expectedStatus, const char *scenario)
{
    if (results.size() == 1 && results.front().generation == generation &&
	bobol_lod_result_matches_request(results.front(), request) &&
	results.front().providerStatus == expectedStatus)
	return 0;
    printf("FAIL: shared producer consumer result (%s): count=%zu\n",
	   scenario, results.size());
    return 1;
}

static int
test_shared_producer_two_generation_delivery(void)
{
    SharedProducerLeaseFixture fixture;
    if (!fixture.service.start(1, FALSE))
	return 1;
    const uint64_t firstGeneration = fixture.service.beginGeneration();
    const uint64_t secondGeneration = fixture.service.beginGeneration();
    BObolLodTask first = shared_producer_lease_task_for(
	fixture, firstGeneration, "/lease/first");
    BObolLodTask second = shared_producer_lease_task_for(
	fixture, secondGeneration, "/lease/second");
    if (!fixture.service.submitIfNotActive(first) ||
	fixture.service.submitIfNotActive(second) != 0 ||
	fixture.service.sharedProducerCountForDiagnostics() != 1 ||
	fixture.service.sharedProducerLeaseCountForDiagnostics() != 2 ||
	fixture.service.sharedProducerLeaseCountForGeneration(
	    firstGeneration) != 1 ||
	fixture.service.sharedProducerLeaseCountForGeneration(
	    secondGeneration) != 1) {
	printf("FAIL: two-generation shared producer lease acquisition\n");
	return 1;
    }
    release_shared_producer(fixture);
    if (wait_for_settled(fixture.service, 2))
	return 1;

    std::vector<BObolLodResult> firstResults;
    std::vector<BObolLodResult> secondResults;
    fixture.service.drainGenerationResults(firstResults, firstGeneration);
    fixture.service.drainGenerationResults(secondResults, secondGeneration);
    if (check_shared_consumer_result(firstResults, firstGeneration,
	    first.request, BOBOL_LOD_PROVIDER_READY, "producer delivery") ||
	check_shared_consumer_result(secondResults, secondGeneration,
	    second.request, BOBOL_LOD_PROVIDER_SUPERSEDED,
	    "secondary replay") ||
	fixture.state.calls != 1 ||
	check_shared_producer_clean(fixture.service,
	    "two-generation delivery"))
	return 1;
    return 0;
}

static int
test_shared_producer_first_consumer_cancellation(void)
{
    SharedProducerLeaseFixture fixture;
    if (!fixture.service.start(1, FALSE))
	return 1;
    const uint64_t firstGeneration = fixture.service.beginGeneration();
    const uint64_t secondGeneration = fixture.service.beginGeneration();
    BObolLodTask first = shared_producer_lease_task_for(
	fixture, firstGeneration, "/lease/cancelled-first");
    BObolLodTask second = shared_producer_lease_task_for(
	fixture, secondGeneration, "/lease/survivor");
    if (!fixture.service.submitIfNotActive(first) ||
	wait_for_shared_producer_started(fixture, "cancel first") ||
	!fixture.service.updateActiveRequestDemand(
	    second.request, secondGeneration))
	return 1;
    fixture.service.cancelGeneration(firstGeneration);
    if (!fixture.service.isGenerationCancelled(firstGeneration) ||
	fixture.service.isProducerCancelled(
	    firstGeneration, first.request) ||
	fixture.service.sharedProducerLeaseCountForDiagnostics() != 1) {
	printf("FAIL: surviving lease did not retain cancelled owner's producer\n");
	return 1;
    }
    release_shared_producer(fixture);
    if (wait_for_settled(fixture.service, 1))
	return 1;
    std::vector<BObolLodResult> cancelledResults;
    std::vector<BObolLodResult> survivorResults;
    fixture.service.drainGenerationResults(
	cancelledResults, firstGeneration);
    fixture.service.drainGenerationResults(
	survivorResults, secondGeneration);
    if (!cancelledResults.empty() ||
	check_shared_consumer_result(survivorResults, secondGeneration,
	    second.request, BOBOL_LOD_PROVIDER_SUPERSEDED,
	    "survivor replay") ||
	fixture.state.calls != 1 ||
	check_shared_producer_clean(fixture.service,
	    "first consumer cancellation"))
	return 1;
    return 0;
}

static int
test_shared_producer_late_join_during_build(void)
{
    SharedProducerLeaseFixture fixture;
    if (!fixture.service.start(1, FALSE))
	return 1;
    const uint64_t firstGeneration = fixture.service.beginGeneration();
    BObolLodTask first = shared_producer_lease_task_for(
	fixture, firstGeneration, "/lease/build-owner");
    if (!fixture.service.submitIfNotActive(first) ||
	wait_for_shared_producer_started(fixture, "late build join"))
	return 1;
    const uint64_t lateGeneration = fixture.service.beginGeneration();
    BObolLodTask late = shared_producer_lease_task_for(
	fixture, lateGeneration, "/lease/build-late");
    if (!fixture.service.updateActiveRequestDemand(
	    late.request, lateGeneration) ||
	fixture.service.sharedProducerLeaseCountForDiagnostics() != 2) {
	printf("FAIL: late consumer did not join building producer\n");
	return 1;
    }
    release_shared_producer(fixture);
    if (wait_for_settled(fixture.service, 2))
	return 1;
    std::vector<BObolLodResult> ownerResults;
    std::vector<BObolLodResult> lateResults;
    fixture.service.drainGenerationResults(ownerResults, firstGeneration);
    fixture.service.drainGenerationResults(lateResults, lateGeneration);
    if (check_shared_consumer_result(ownerResults, firstGeneration,
	    first.request, BOBOL_LOD_PROVIDER_READY, "build owner") ||
	check_shared_consumer_result(lateResults, lateGeneration,
	    late.request, BOBOL_LOD_PROVIDER_SUPERSEDED,
	    "late build replay") ||
	fixture.state.calls != 1 ||
	check_shared_producer_clean(fixture.service, "late build join"))
	return 1;
    return 0;
}

static int
test_shared_producer_late_join_while_result_queued(void)
{
    SharedProducerLeaseFixture fixture;
    if (!fixture.service.start(1, FALSE))
	return 1;
    const uint64_t firstGeneration = fixture.service.beginGeneration();
    BObolLodTask first = shared_producer_lease_task_for(
	fixture, firstGeneration, "/lease/result-owner");
    if (!fixture.service.submitIfNotActive(first) ||
	wait_for_shared_producer_started(fixture, "late result join"))
	return 1;
    release_shared_producer(fixture);
    if (wait_for_settled(fixture.service, 1))
	return 1;

    const uint64_t lateGeneration = fixture.service.beginGeneration();
    BObolLodTask late = shared_producer_lease_task_for(
	fixture, lateGeneration, "/lease/result-late");
    if (fixture.service.submitIfNotActive(late) != 0 ||
	fixture.service.queuedResultCountForDiagnostics() != 2 ||
	fixture.service.sharedProducerLeaseCountForDiagnostics() != 2) {
	printf("FAIL: late consumer did not lease queued shared result\n");
	return 1;
    }
    std::vector<BObolLodResult> lateResults;
    std::vector<BObolLodResult> ownerResults;
    fixture.service.drainGenerationResults(lateResults, lateGeneration);
    fixture.service.drainGenerationResults(ownerResults, firstGeneration);
    if (check_shared_consumer_result(lateResults, lateGeneration,
	    late.request, BOBOL_LOD_PROVIDER_SUPERSEDED,
	    "late queued-result replay") ||
	check_shared_consumer_result(ownerResults, firstGeneration,
	    first.request, BOBOL_LOD_PROVIDER_READY,
	    "queued-result owner") ||
	fixture.state.calls != 1 ||
	check_shared_producer_clean(fixture.service,
	    "late queued-result join"))
	return 1;
    return 0;
}

static int
test_shared_producer_final_consumer_cancellation_during_build(void)
{
    SharedProducerLeaseFixture fixture;
    if (!fixture.service.start(1, FALSE))
	return 1;
    const uint64_t generation = fixture.service.beginGeneration();
    BObolLodTask task = shared_producer_lease_task_for(
	fixture, generation, "/lease/final-build");
    if (!fixture.service.submitIfNotActive(task) ||
	wait_for_shared_producer_started(fixture, "final build cancel"))
	return 1;
    fixture.service.cancelGeneration(generation);
    if (!fixture.service.isProducerCancelled(generation, task.request) ||
	fixture.service.sharedProducerLeaseCountForDiagnostics() != 0) {
	printf("FAIL: final build lease cancellation did not cancel producer\n");
	return 1;
    }
    release_shared_producer(fixture);
    if (wait_for_settled(fixture.service, 0) ||
	fixture.state.calls != 1 ||
	check_shared_producer_clean(fixture.service,
	    "final cancellation during build"))
	return 1;
    return 0;
}

static int
test_shared_producer_final_consumer_cancellation_with_result(void)
{
    SharedProducerLeaseFixture fixture;
    if (!fixture.service.start(1, FALSE))
	return 1;
    const uint64_t generation = fixture.service.beginGeneration();
    BObolLodTask task = shared_producer_lease_task_for(
	fixture, generation, "/lease/final-result");
    if (!fixture.service.submitIfNotActive(task) ||
	wait_for_shared_producer_started(fixture, "final result cancel"))
	return 1;
    release_shared_producer(fixture);
    if (wait_for_settled(fixture.service, 1))
	return 1;
    fixture.service.cancelGeneration(generation);
    if (fixture.state.calls != 1 ||
	check_shared_producer_clean(fixture.service,
	    "final cancellation with queued result"))
	return 1;
    return 0;
}

static int
test_batch_submission(void)
{
    BObolLodService service;
    ServiceTestContext context;
    TaskData firstData{&context, 1};
    TaskData duplicateData{&context, 2};
    TaskData routedData{&context, 3};
    TaskData populationData{&context, 4};

    if (!service.start(2, TRUE)) {
	printf("FAIL: LoD service did not start for batch submission test\n");
	return 1;
    }

    const uint64_t generation = service.beginGeneration();
    BObolLodTask first;
    first.generation = generation;
    first.request = make_request("/batch.bot");
    first.request.occurrenceKey = "/batch.bot#0";
    first.realize = ready_task;
    first.realizeData = &firstData;

    BObolLodTask duplicate = first;
    duplicate.request.viewRevision++;
    duplicate.realizeData = &duplicateData;

    BObolLodTask routed = first;
    routed.request.sourceRoutingId = 91;
    routed.realizeData = &routedData;

    BObolLodTask repopulated = routed;
    repopulated.request.sourcePopulationEpoch = 2;
    repopulated.realizeData = &populationData;

    std::vector<BObolLodTask> tasks{
	first, duplicate, routed, repopulated};
    std::vector<uint64_t> taskIds;
    if (service.submitBatch(tasks, taskIds, TRUE) != 3 ||
	taskIds.size() != tasks.size() ||
	taskIds[0] == 0 || taskIds[1] != 0 || taskIds[2] == 0 ||
	taskIds[3] == 0) {
	printf("FAIL: LoD batch did not suppress only the duplicate route\n");
	service.stop();
	return 1;
    }

    if (wait_for_settled(service, 2)) {
	service.stop();
	return 1;
    }

    std::vector<BObolLodResult> results;
    if (service.drainResults(results) != 3 || results.size() != 3) {
	printf("FAIL: LoD batch did not publish distinct route generations\n");
	service.stop();
	return 1;
    }

    bool foundFirst = false;
    bool foundRouted = false;
    bool foundRepopulated = false;
    for (const BObolLodResult &result : results) {
	if (result.request.sourceRoutingId == 0 &&
	    result.counts.faceCount == 1)
	    foundFirst = true;
	if (result.request.sourceRoutingId == 91 &&
	    result.request.sourcePopulationEpoch == 0 &&
	    result.counts.faceCount == 3)
	    foundRouted = true;
	if (result.request.sourceRoutingId == 91 &&
	    result.request.sourcePopulationEpoch == 2 &&
	    result.counts.faceCount == 4)
	    foundRepopulated = true;
    }
    {
	std::lock_guard<std::mutex> lock(context.mutex);
	if (!foundFirst || !foundRouted || !foundRepopulated ||
	    context.executionOrder.size() != 3 ||
	    std::find(context.executionOrder.begin(),
		context.executionOrder.end(), 2) !=
		context.executionOrder.end()) {
	    printf("FAIL: LoD batch executed a suppressed duplicate\n");
	    service.stop();
	    return 1;
	}
    }

    service.stop();
    return 0;
}

static int
test_rt_source_full_detail_provider_task(void)
{
    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;

    if (make_provider_test_db(dbpath, sizeof(dbpath), &dbip))
	return 1;

    BObolRtSourceFullDetailProvider provider;
    provider.setDatabase(dbip);
    provider.validateSourceMetrics = TRUE;

    BObolLodRequest request = make_request("/lod-provider.bot");
    request.objectName = "lod-provider.bot";
    request.providerId = "rt_source_full_detail";
    request.providerVersion = "direct-bot-v1";
    request.qualityTier = BOBOL_LOD_QUALITY_FULL_DETAIL;
    request.sourceCounts.faceCount = 4;
    request.sourceCounts.pointCount = 4;

    int ret = 0;
    BObolLodResult directResult =
	bobol_rt_source_full_detail_provider_task(request, &provider);
    if (check_source_full_detail_result(directResult, request,
					"did not return direct BoT mesh result"))
	ret = 1;

    BObolLodRequest preparedRequest = request;
    preparedRequest.addProviderParam("display.prepared_cad_only", "1");
    BObolLodResult preparedResult =
	bobol_rt_source_full_detail_provider_task(preparedRequest, &provider);
    if (preparedResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	preparedResult.resultKind != BOBOL_LOD_RESULT_FULL_DETAIL ||
	!preparedResult.preparedCadGeometry ||
	!preparedResult.preparedCadGeometry->shaded ||
	preparedResult.mesh.isValid() ||
	preparedResult.preparedCadGeometryRevision == 0 ||
	preparedResult.counts.faceCount != 4 ||
	preparedResult.counts.pointCount == 0 ||
	preparedResult.counts.byteCount == 0 ||
	!bobol_lod_result_matches_request(preparedResult, preparedRequest)) {
	printf("FAIL: LoD RT source full-detail provider did not prepare terminal CAD geometry\n");
	ret = 1;
    }

    BObolLodRequest scopedRequest = request;
    scopedRequest.addProviderParam("source_query.space", "source_local");
    scopedRequest.addProviderParam("source_query.bounds",
				   "0.7 0.15 -0.05 0.8 0.25 0.05");
    scopedRequest.addProviderParam("source_query.tolerance", "0.05");
    BObolLodResult scopedResult =
	bobol_rt_source_full_detail_provider_task(scopedRequest, &provider);
    if (scopedResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	scopedResult.resultKind != BOBOL_LOD_RESULT_FULL_DETAIL ||
	scopedResult.mesh.points.size() != 4 ||
	scopedResult.mesh.coordIndex.size() != 6 ||
	scopedResult.mesh.faceIndex.size() != 2 ||
	scopedResult.mesh.faceIndex[0] != 0 ||
	scopedResult.mesh.faceIndex[1] != 2 ||
	scopedResult.mesh.vertexIndex.size() != 4 ||
	scopedResult.mesh.vertexIndex[0] != 0 ||
	scopedResult.mesh.vertexIndex[1] != 1 ||
	scopedResult.mesh.vertexIndex[2] != 2 ||
	scopedResult.mesh.vertexIndex[3] != 3 ||
	scopedResult.counts.faceCount != 2 ||
	!bobol_lod_result_matches_request(scopedResult, scopedRequest)) {
	printf("FAIL: LoD RT source full-detail provider query bounds did not reduce returned face payload\n");
	ret = 1;
    }

    BObolRtSourceFullDetailProvider scopedLimitProvider;
    scopedLimitProvider.setDatabase(dbip);
    scopedLimitProvider.validateSourceMetrics = TRUE;
    scopedLimitProvider.maxFullDetailFaceCount = 2;
    scopedLimitProvider.maxFullDetailPointCount = 4;
    BObolLodResult scopedLimitedResult =
	bobol_rt_source_full_detail_provider_task(scopedRequest,
	    &scopedLimitProvider);
    if (scopedLimitedResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	scopedLimitedResult.resultKind != BOBOL_LOD_RESULT_FULL_DETAIL ||
	scopedLimitedResult.mesh.points.size() != 4 ||
	scopedLimitedResult.mesh.coordIndex.size() != 6 ||
	scopedLimitedResult.mesh.faceIndex.size() != 2 ||
	scopedLimitedResult.counts.faceCount != 2 ||
	scopedLimitedResult.counts.pointCount != 4 ||
	!bobol_lod_result_matches_request(scopedLimitedResult,
					    scopedRequest)) {
	printf("FAIL: LoD RT source full-detail provider scoped bounds limit should apply after payload reduction\n");
	ret = 1;
    }

    BObolLodRequest emptyScopedRequest = request;
    emptyScopedRequest.addProviderParam("source_query.space", "source_local");
    emptyScopedRequest.addProviderParam("source_query.bounds",
					"5.0 5.0 5.0 5.1 5.1 5.1");
    emptyScopedRequest.addProviderParam("source_query.tolerance", "0.05");
    BObolLodResult emptyScopedResult =
	bobol_rt_source_full_detail_provider_task(emptyScopedRequest,
	    &provider);
    if (emptyScopedResult.providerStatus != BOBOL_LOD_PROVIDER_FALLBACK ||
	emptyScopedResult.resultKind != BOBOL_LOD_RESULT_NONE ||
	emptyScopedResult.counts.faceCount != 0 ||
	emptyScopedResult.counts.pointCount != 0 ||
	emptyScopedResult.mesh.isValid() ||
	bu_strcmp(emptyScopedResult.diagnostic.getString(),
	       "RT source full-detail provider scoped query matched no faces") != 0 ||
	!bobol_lod_result_matches_request(emptyScopedResult,
					    emptyScopedRequest)) {
	printf("FAIL: LoD RT source full-detail provider scoped bounds miss should not expand to whole-object payload\n");
	ret = 1;
    }

    BObolLodRequest wrongSpaceBoundsRequest = request;
    wrongSpaceBoundsRequest.addProviderParam("source_query.space", "world");
    wrongSpaceBoundsRequest.addProviderParam("source_query.bounds",
	    "0.7 0.15 -0.05 0.8 0.25 0.05");
    wrongSpaceBoundsRequest.addProviderParam("source_query.tolerance",
	    "0.05");
    BObolLodResult wrongSpaceBoundsResult =
	bobol_rt_source_full_detail_provider_task(wrongSpaceBoundsRequest,
	    &provider);
    if (wrongSpaceBoundsResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	wrongSpaceBoundsResult.resultKind !=
	BOBOL_LOD_RESULT_FULL_DETAIL ||
	wrongSpaceBoundsResult.mesh.points.size() != 4 ||
	wrongSpaceBoundsResult.mesh.coordIndex.size() != 12 ||
	wrongSpaceBoundsResult.counts.faceCount != 4 ||
	wrongSpaceBoundsResult.counts.pointCount != 4 ||
	!bobol_lod_result_matches_request(wrongSpaceBoundsResult,
					    wrongSpaceBoundsRequest)) {
	printf("FAIL: LoD RT source full-detail provider should ignore non-source-local bounds when reducing payloads\n");
	ret = 1;
    }

    BObolRtSourceFullDetailProvider wrongSpaceLimitProvider;
    wrongSpaceLimitProvider.setDatabase(dbip);
    wrongSpaceLimitProvider.validateSourceMetrics = TRUE;
    wrongSpaceLimitProvider.maxFullDetailFaceCount = 2;
    wrongSpaceLimitProvider.maxFullDetailPointCount = 4;
    BObolLodResult wrongSpaceLimitedResult =
	bobol_rt_source_full_detail_provider_task(wrongSpaceBoundsRequest,
	    &wrongSpaceLimitProvider);
    if (wrongSpaceLimitedResult.providerStatus !=
	BOBOL_LOD_PROVIDER_FALLBACK ||
	wrongSpaceLimitedResult.resultKind != BOBOL_LOD_RESULT_NONE ||
	wrongSpaceLimitedResult.mesh.isValid() ||
	bu_strcmp(wrongSpaceLimitedResult.diagnostic.getString(),
	       "RT source full-detail provider request exceeds full-detail limits") != 0 ||
	!bobol_lod_result_matches_request(wrongSpaceLimitedResult,
					    wrongSpaceBoundsRequest)) {
	printf("FAIL: LoD RT source full-detail provider should not bypass whole-object limits for non-source-local bounds\n");
	ret = 1;
    }

    BObolLodRequest malformedBoundsRequest = request;
    malformedBoundsRequest.addProviderParam("source_query.space",
					    "source_local");
    malformedBoundsRequest.addProviderParam("source_query.bounds",
					    "0.7 0.15 -0.05 0.8 0.25 0.05 trailing");
    malformedBoundsRequest.addProviderParam("source_query.tolerance",
					    "0.05");
    BObolLodResult malformedBoundsResult =
	bobol_rt_source_full_detail_provider_task(malformedBoundsRequest,
	    &provider);
    if (malformedBoundsResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	malformedBoundsResult.resultKind !=
	BOBOL_LOD_RESULT_FULL_DETAIL ||
	malformedBoundsResult.mesh.points.size() != 4 ||
	malformedBoundsResult.mesh.coordIndex.size() != 12 ||
	malformedBoundsResult.counts.faceCount != 4 ||
	malformedBoundsResult.counts.pointCount != 4 ||
	!bobol_lod_result_matches_request(malformedBoundsResult,
					    malformedBoundsRequest)) {
	printf("FAIL: LoD RT source full-detail provider should ignore malformed source-local bounds when reducing payloads\n");
	ret = 1;
    }

    BObolLodRequest malformedToleranceRequest = request;
    malformedToleranceRequest.addProviderParam("source_query.space",
	    "source_local");
    malformedToleranceRequest.addProviderParam("source_query.bounds",
	    "0.7 0.15 -0.05 0.8 0.25 0.05");
    malformedToleranceRequest.addProviderParam("source_query.tolerance",
	    "0.05 trailing");
    BObolLodResult malformedToleranceLimitedResult =
	bobol_rt_source_full_detail_provider_task(malformedToleranceRequest,
	    &wrongSpaceLimitProvider);
    if (malformedToleranceLimitedResult.providerStatus !=
	BOBOL_LOD_PROVIDER_FALLBACK ||
	malformedToleranceLimitedResult.resultKind !=
	BOBOL_LOD_RESULT_NONE ||
	malformedToleranceLimitedResult.mesh.isValid() ||
	bu_strcmp(malformedToleranceLimitedResult.diagnostic.getString(),
	       "RT source full-detail provider request exceeds full-detail limits") != 0 ||
	!bobol_lod_result_matches_request(
	    malformedToleranceLimitedResult, malformedToleranceRequest)) {
	printf("FAIL: LoD RT source full-detail provider should not bypass whole-object limits for malformed source-local tolerance\n");
	ret = 1;
    }

    BObolLodRequest duplicateSpaceBoundsRequest = request;
    duplicateSpaceBoundsRequest.addProviderParam("source_query.space",
	    "source_local");
    duplicateSpaceBoundsRequest.addProviderParam("source_query.space",
	    "source_local");
    duplicateSpaceBoundsRequest.addProviderParam("source_query.bounds",
	    "0.7 0.15 -0.05 0.8 0.25 0.05");
    duplicateSpaceBoundsRequest.addProviderParam("source_query.tolerance",
	    "0.05");
    BObolLodResult duplicateSpaceBoundsResult =
	bobol_rt_source_full_detail_provider_task(
	    duplicateSpaceBoundsRequest, &provider);
    if (duplicateSpaceBoundsResult.providerStatus !=
	BOBOL_LOD_PROVIDER_READY ||
	duplicateSpaceBoundsResult.resultKind !=
	BOBOL_LOD_RESULT_FULL_DETAIL ||
	duplicateSpaceBoundsResult.mesh.points.size() != 4 ||
	duplicateSpaceBoundsResult.mesh.coordIndex.size() != 12 ||
	duplicateSpaceBoundsResult.counts.faceCount != 4 ||
	duplicateSpaceBoundsResult.counts.pointCount != 4 ||
	!bobol_lod_result_matches_request(duplicateSpaceBoundsResult,
					    duplicateSpaceBoundsRequest)) {
	printf("FAIL: LoD RT source full-detail provider should ignore duplicate query-space params when reducing payloads\n");
	ret = 1;
    }

    BObolLodRequest duplicateBoundsRequest = request;
    duplicateBoundsRequest.addProviderParam("source_query.space",
					    "source_local");
    duplicateBoundsRequest.addProviderParam("source_query.bounds",
					    "0.7 0.15 -0.05 0.8 0.25 0.05");
    duplicateBoundsRequest.addProviderParam("source_query.bounds",
					    "0.7 0.15 -0.05 0.8 0.25 0.05");
    duplicateBoundsRequest.addProviderParam("source_query.tolerance",
					    "0.05");
    BObolLodResult duplicateBoundsLimitedResult =
	bobol_rt_source_full_detail_provider_task(duplicateBoundsRequest,
	    &wrongSpaceLimitProvider);
    if (duplicateBoundsLimitedResult.providerStatus !=
	BOBOL_LOD_PROVIDER_FALLBACK ||
	duplicateBoundsLimitedResult.resultKind !=
	BOBOL_LOD_RESULT_NONE ||
	duplicateBoundsLimitedResult.mesh.isValid() ||
	bu_strcmp(duplicateBoundsLimitedResult.diagnostic.getString(),
	       "RT source full-detail provider request exceeds full-detail limits") != 0 ||
	!bobol_lod_result_matches_request(duplicateBoundsLimitedResult,
					    duplicateBoundsRequest)) {
	printf("FAIL: LoD RT source full-detail provider should not bypass whole-object limits for duplicate bounds params\n");
	ret = 1;
    }

    BObolLodRequest mixedScopedRequest = request;
    mixedScopedRequest.addProviderParam("source_query.space",
					"source_local");
    mixedScopedRequest.addProviderParam("source_query.bounds",
					"0.7 0.15 -0.05 0.8 0.25 0.05");
    mixedScopedRequest.addProviderParam("source_query.tolerance", "0.05");
    mixedScopedRequest.addProviderParam("source_query.ray.origin",
					"0.2 0.2 5.0");
    mixedScopedRequest.addProviderParam("source_query.ray.direction",
					"0.0 0.0 -1.0");
    BObolLodResult mixedScopedResult =
	bobol_rt_source_full_detail_provider_task(mixedScopedRequest,
	    &provider);
    if (mixedScopedResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	mixedScopedResult.resultKind != BOBOL_LOD_RESULT_FULL_DETAIL ||
	mixedScopedResult.mesh.points.size() != 4 ||
	mixedScopedResult.mesh.coordIndex.size() != 12 ||
	mixedScopedResult.counts.faceCount != 4 ||
	mixedScopedResult.counts.pointCount != 4 ||
	!bobol_lod_result_matches_request(mixedScopedResult,
					    mixedScopedRequest)) {
	printf("FAIL: LoD RT source full-detail provider should ignore mixed scoped query kinds when reducing payloads\n");
	ret = 1;
    }

    BObolLodResult mixedScopedLimitedResult =
	bobol_rt_source_full_detail_provider_task(mixedScopedRequest,
	    &wrongSpaceLimitProvider);
    if (mixedScopedLimitedResult.providerStatus !=
	BOBOL_LOD_PROVIDER_FALLBACK ||
	mixedScopedLimitedResult.resultKind != BOBOL_LOD_RESULT_NONE ||
	mixedScopedLimitedResult.mesh.isValid() ||
	bu_strcmp(mixedScopedLimitedResult.diagnostic.getString(),
	       "RT source full-detail provider request exceeds full-detail limits") != 0 ||
	!bobol_lod_result_matches_request(mixedScopedLimitedResult,
					    mixedScopedRequest)) {
	printf("FAIL: LoD RT source full-detail provider should not bypass whole-object limits for mixed scoped query kinds\n");
	ret = 1;
    }

    BObolLodRequest rayRequest = request;
    rayRequest.addProviderParam("source_query.space", "source_local");
    rayRequest.addProviderParam("source_query.ray.origin",
				"0.2 0.2 5.0");
    rayRequest.addProviderParam("source_query.ray.direction",
				"0.0 0.0 -1.0");
    BObolLodResult rayResult =
	bobol_rt_source_full_detail_provider_task(rayRequest, &provider);
    if (rayResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	rayResult.resultKind != BOBOL_LOD_RESULT_FULL_DETAIL ||
	rayResult.mesh.points.size() != 4 ||
	rayResult.mesh.coordIndex.size() != 6 ||
	rayResult.mesh.faceIndex.size() != 2 ||
	rayResult.mesh.faceIndex[0] != 0 ||
	rayResult.mesh.faceIndex[1] != 2 ||
	rayResult.mesh.vertexIndex.size() != 4 ||
	rayResult.mesh.vertexIndex[0] != 0 ||
	rayResult.mesh.vertexIndex[1] != 1 ||
	rayResult.mesh.vertexIndex[2] != 2 ||
	rayResult.mesh.vertexIndex[3] != 3 ||
	rayResult.counts.faceCount != 2 ||
	!bobol_lod_result_matches_request(rayResult, rayRequest)) {
	printf("FAIL: LoD RT source full-detail provider query ray did not reduce returned face payload\n");
	ret = 1;
    }

    BObolLodRequest compactRayRequest = make_request("/lod-two-tri.bot");
    compactRayRequest.objectName = "lod-two-tri.bot";
    compactRayRequest.providerId = "rt_source_full_detail";
    compactRayRequest.providerVersion = "direct-bot-v1";
    compactRayRequest.qualityTier = BOBOL_LOD_QUALITY_FULL_DETAIL;
    compactRayRequest.sourceCounts.faceCount = 2;
    compactRayRequest.sourceCounts.pointCount = 6;
    compactRayRequest.addProviderParam("source_query.space", "source_local");
    compactRayRequest.addProviderParam("source_query.ray.origin",
				       "0.2 0.2 5.0");
    compactRayRequest.addProviderParam("source_query.ray.direction",
				       "0.0 0.0 -1.0");
    BObolLodResult compactRayResult =
	bobol_rt_source_full_detail_provider_task(compactRayRequest,
	    &provider);
    if (compactRayResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	compactRayResult.resultKind != BOBOL_LOD_RESULT_FULL_DETAIL ||
	compactRayResult.mesh.points.size() != 3 ||
	compactRayResult.mesh.coordIndex.size() != 3 ||
	compactRayResult.mesh.coordIndex[0] != 0 ||
	compactRayResult.mesh.coordIndex[1] != 1 ||
	compactRayResult.mesh.coordIndex[2] != 2 ||
	compactRayResult.mesh.faceIndex.size() != 1 ||
	compactRayResult.mesh.faceIndex[0] != 0 ||
	compactRayResult.mesh.vertexIndex.size() != 3 ||
	compactRayResult.mesh.vertexIndex[0] != 0 ||
	compactRayResult.mesh.vertexIndex[1] != 1 ||
	compactRayResult.mesh.vertexIndex[2] != 2 ||
	compactRayResult.counts.faceCount != 1 ||
	compactRayResult.counts.pointCount != 3 ||
	!bobol_lod_result_matches_request(compactRayResult,
					    compactRayRequest)) {
	printf("FAIL: LoD RT source full-detail provider query ray did not compact source vertex payload\n");
	ret = 1;
    }

    BObolRtSourceFullDetailProvider compactRayLimitProvider;
    compactRayLimitProvider.setDatabase(dbip);
    compactRayLimitProvider.validateSourceMetrics = TRUE;
    compactRayLimitProvider.maxFullDetailFaceCount = 1;
    compactRayLimitProvider.maxFullDetailPointCount = 3;
    BObolLodResult compactRayLimitedResult =
	bobol_rt_source_full_detail_provider_task(compactRayRequest,
	    &compactRayLimitProvider);
    if (compactRayLimitedResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	compactRayLimitedResult.resultKind != BOBOL_LOD_RESULT_FULL_DETAIL ||
	compactRayLimitedResult.mesh.points.size() != 3 ||
	compactRayLimitedResult.mesh.coordIndex.size() != 3 ||
	compactRayLimitedResult.mesh.faceIndex.size() != 1 ||
	compactRayLimitedResult.mesh.vertexIndex.size() != 3 ||
	compactRayLimitedResult.counts.faceCount != 1 ||
	compactRayLimitedResult.counts.pointCount != 3 ||
	!bobol_lod_result_matches_request(compactRayLimitedResult,
					    compactRayRequest)) {
	printf("FAIL: LoD RT source full-detail provider scoped ray limit should apply after payload reduction\n");
	ret = 1;
    }

    BObolLodRequest missRayRequest = compactRayRequest;
    missRayRequest.providerParams.clear();
    missRayRequest.addProviderParam("source_query.space", "source_local");
    missRayRequest.addProviderParam("source_query.ray.origin",
				    "5.0 5.0 5.0");
    missRayRequest.addProviderParam("source_query.ray.direction",
				    "0.0 0.0 -1.0");
    BObolLodResult missRayResult =
	bobol_rt_source_full_detail_provider_task(missRayRequest, &provider);
    if (missRayResult.providerStatus != BOBOL_LOD_PROVIDER_FALLBACK ||
	missRayResult.resultKind != BOBOL_LOD_RESULT_NONE ||
	missRayResult.counts.faceCount != 0 ||
	missRayResult.counts.pointCount != 0 ||
	missRayResult.mesh.isValid() ||
	bu_strcmp(missRayResult.diagnostic.getString(),
	       "RT source full-detail provider scoped query matched no faces") != 0 ||
	!bobol_lod_result_matches_request(missRayResult,
					    missRayRequest)) {
	printf("FAIL: LoD RT source full-detail provider scoped ray miss should not expand to whole-object payload\n");
	ret = 1;
    }

    BObolLodRequest wrongSpaceRayRequest = compactRayRequest;
    wrongSpaceRayRequest.providerParams.clear();
    wrongSpaceRayRequest.addProviderParam("source_query.space", "world");
    wrongSpaceRayRequest.addProviderParam("source_query.ray.origin",
					  "0.2 0.2 5.0");
    wrongSpaceRayRequest.addProviderParam("source_query.ray.direction",
					  "0.0 0.0 -1.0");
    BObolLodResult wrongSpaceRayResult =
	bobol_rt_source_full_detail_provider_task(wrongSpaceRayRequest,
	    &provider);
    if (wrongSpaceRayResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	wrongSpaceRayResult.resultKind != BOBOL_LOD_RESULT_FULL_DETAIL ||
	wrongSpaceRayResult.mesh.points.size() != 6 ||
	wrongSpaceRayResult.mesh.coordIndex.size() != 6 ||
	!wrongSpaceRayResult.mesh.vertexIndex.empty() ||
	wrongSpaceRayResult.counts.faceCount != 2 ||
	wrongSpaceRayResult.counts.pointCount != 6 ||
	!bobol_lod_result_matches_request(wrongSpaceRayResult,
					    wrongSpaceRayRequest)) {
	printf("FAIL: LoD RT source full-detail provider should ignore non-source-local rays when reducing payloads\n");
	ret = 1;
    }

    BObolLodRequest malformedRayRequest = compactRayRequest;
    malformedRayRequest.providerParams.clear();
    malformedRayRequest.addProviderParam("source_query.space",
					 "source_local");
    malformedRayRequest.addProviderParam("source_query.ray.origin",
					 "0.2 0.2 5.0");
    malformedRayRequest.addProviderParam("source_query.ray.direction",
					 "0.0 0.0 -1.0 trailing");
    BObolLodResult malformedRayResult =
	bobol_rt_source_full_detail_provider_task(malformedRayRequest,
	    &provider);
    if (malformedRayResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	malformedRayResult.resultKind != BOBOL_LOD_RESULT_FULL_DETAIL ||
	malformedRayResult.mesh.points.size() != 6 ||
	malformedRayResult.mesh.coordIndex.size() != 6 ||
	!malformedRayResult.mesh.vertexIndex.empty() ||
	malformedRayResult.counts.faceCount != 2 ||
	malformedRayResult.counts.pointCount != 6 ||
	!bobol_lod_result_matches_request(malformedRayResult,
					    malformedRayRequest)) {
	printf("FAIL: LoD RT source full-detail provider should ignore malformed source-local rays when reducing payloads\n");
	ret = 1;
    }

    BObolLodRequest duplicateRayRequest = compactRayRequest;
    duplicateRayRequest.providerParams.clear();
    duplicateRayRequest.addProviderParam("source_query.space",
					 "source_local");
    duplicateRayRequest.addProviderParam("source_query.ray.origin",
					 "0.2 0.2 5.0");
    duplicateRayRequest.addProviderParam("source_query.ray.direction",
					 "0.0 0.0 -1.0");
    duplicateRayRequest.addProviderParam("source_query.ray.direction",
					 "0.0 0.0 -1.0");
    BObolLodResult duplicateRayResult =
	bobol_rt_source_full_detail_provider_task(duplicateRayRequest,
	    &provider);
    if (duplicateRayResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	duplicateRayResult.resultKind != BOBOL_LOD_RESULT_FULL_DETAIL ||
	duplicateRayResult.mesh.points.size() != 6 ||
	duplicateRayResult.mesh.coordIndex.size() != 6 ||
	!duplicateRayResult.mesh.vertexIndex.empty() ||
	duplicateRayResult.counts.faceCount != 2 ||
	duplicateRayResult.counts.pointCount != 6 ||
	!bobol_lod_result_matches_request(duplicateRayResult,
					    duplicateRayRequest)) {
	printf("FAIL: LoD RT source full-detail provider should ignore duplicate ray params when reducing payloads\n");
	ret = 1;
    }

    BObolLodRequest measureHintRequest = request;
    measureHintRequest.addProviderParam("source_query.space", "source_local");
    measureHintRequest.addProviderParam("source_query.bounds",
					"0.7 0.15 -0.05 0.8 0.25 0.05");
    BObolLodResult measureHintResult =
	bobol_rt_source_full_detail_provider_task(measureHintRequest, &provider);
    if (measureHintResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	measureHintResult.mesh.coordIndex.size() != 12 ||
	measureHintResult.mesh.faceIndex.size() != 4 ||
	!measureHintResult.mesh.vertexIndex.empty() ||
	measureHintResult.counts.faceCount != 4 ||
	!bobol_lod_result_matches_request(measureHintResult,
					    measureHintRequest)) {
	printf("FAIL: LoD RT source full-detail provider should keep measure query hints whole-object without tolerance\n");
	ret = 1;
    }

    BObolLodRequest boundedMeasureRequest = request;
    boundedMeasureRequest.addProviderParam("source_query.space",
					   "source_local");
    boundedMeasureRequest.addProviderParam("source_query.bounds",
					   "0.7 0.15 -0.05 0.8 0.25 0.05");
    boundedMeasureRequest.addProviderParam("source_query.tolerance", "0.05");
    BObolLodResult boundedMeasureResult =
	bobol_rt_source_full_detail_provider_task(boundedMeasureRequest,
	    &provider);
    if (boundedMeasureResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	boundedMeasureResult.mesh.points.size() != 4 ||
	boundedMeasureResult.mesh.coordIndex.size() != 6 ||
	boundedMeasureResult.mesh.faceIndex.size() != 2 ||
	boundedMeasureResult.mesh.vertexIndex.size() != 4 ||
	boundedMeasureResult.counts.faceCount != 2 ||
	boundedMeasureResult.counts.pointCount != 4 ||
	!bobol_lod_result_matches_request(boundedMeasureResult,
					    boundedMeasureRequest)) {
	printf("FAIL: LoD RT source full-detail provider should reduce explicit bounded measure query payloads\n");
	ret = 1;
    }

    BObolLodService service;
    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD RT source full-detail provider service did not start\n");
	ret = 1;
    } else {
	BObolLodTask task;
	task.generation = service.beginGeneration();
	task.request = request;
	task.realize = bobol_rt_source_full_detail_provider_task;
	task.realizeData = &provider;

	if (service.submit(task) == 0) {
	    printf("FAIL: LoD RT source full-detail provider service did not accept task\n");
	    ret = 1;
	} else if (wait_for_settled(service, 1)) {
	    ret = 1;
	} else {
	    std::vector<BObolLodResult> results;
	    service.drainResults(results);
	    if (results.size() != 1 ||
		check_source_full_detail_result(results[0], request,
						"did not publish service BoT mesh result"))
		ret = 1;
	}
	service.stop();
    }

    BObolLodRequest staleRequest = request;
    staleRequest.sourceCounts.faceCount = 99;
    BObolLodResult staleResult =
	bobol_rt_source_full_detail_provider_task(staleRequest, &provider);
    if (staleResult.providerStatus != BOBOL_LOD_PROVIDER_STALE ||
	!staleResult.stale ||
	bu_strcmp(staleResult.diagnostic.getString(),
	       "RT source full-detail provider source metrics changed") != 0 ||
	staleResult.mesh.isValid()) {
	printf("FAIL: LoD RT source full-detail provider did not reject stale source metrics\n");
	ret = 1;
    }

    BObolSourceMeshRequest sourceRequest;
    sourceRequest.path = "/lod-provider.bot";
    sourceRequest.sourceName = "lod-provider.bot";
    sourceRequest.sourceType = "bot";
    sourceRequest.sourceId = 7001;
    sourceRequest.faceCount = 4;
    sourceRequest.pointCount = 4;
    sourceRequest.bounds = request.bounds;

    BObolLodRequest convertedRequest;
    if (!bobol_lod_rt_source_full_detail_request_from_source_mesh_request(
	    convertedRequest, sourceRequest, &request) ||
	bu_strcmp(convertedRequest.providerId.getString(),
	       "rt_source_full_detail") != 0 ||
	bu_strcmp(convertedRequest.providerVersion.getString(),
	       "direct-bot-v1") != 0 ||
	convertedRequest.qualityTier != BOBOL_LOD_QUALITY_FULL_DETAIL ||
	convertedRequest.sourceCounts.faceCount != 4 ||
	convertedRequest.sourceCounts.pointCount != 4 ||
	bu_strcmp(convertedRequest.objectPath.getString(),
	       "/lod-provider.bot") != 0 ||
	bu_strcmp(convertedRequest.objectName.getString(),
	       "lod-provider.bot") != 0) {
	printf("FAIL: LoD RT source full-detail helper did not convert source request\n");
	ret = 1;
    }

    /* Occurrence metadata must remain independent from canonical geometry
     * residency.  A transformed-copy occurrence reports its own path/name,
     * while both progressive and exact providers fetch the representative
     * asset in its own local coordinate system. */
    BObolSourceMeshRequest sharedAssetRequest = sourceRequest;
    sharedAssetRequest.path = "/assembly/copy.bot";
    sharedAssetRequest.sourceName = "copy.bot";
    sharedAssetRequest.meshAssetPath = "/assembly/original.bot";
    sharedAssetRequest.meshAssetName = "original.bot";
    sharedAssetRequest.meshAssetBounds = SbBox3f(
	SbVec3f(-2.0f, -3.0f, -4.0f), SbVec3f(2.0f, 3.0f, 4.0f));
    BObolLodRequest sharedAssetConverted;
    if (!bobol_lod_rt_source_full_detail_request_from_source_mesh_request(
	    sharedAssetConverted, sharedAssetRequest, &request) ||
	bu_strcmp(sharedAssetConverted.objectPath.getString(),
	    "/assembly/original.bot") != 0 ||
	bu_strcmp(sharedAssetConverted.objectName.getString(),
	    "original.bot") != 0 ||
	sharedAssetConverted.bounds != sharedAssetRequest.meshAssetBounds ||
	bu_strcmp(sharedAssetRequest.path.getString(),
	    "/assembly/copy.bot") != 0 ||
	bu_strcmp(sharedAssetRequest.sourceName.getString(), "copy.bot") != 0) {
	printf("FAIL: LoD source helper did not separate occurrence and canonical asset identity\n");
	ret = 1;
    }

    BObolSourceMeshRequest templatedScopedSourceRequest = sourceRequest;
    templatedScopedSourceRequest.queryBoundsValid = 1;
    templatedScopedSourceRequest.queryBounds = SbBox3f(
	    SbVec3f(0.7f, 0.15f, -0.05f),
	    SbVec3f(0.8f, 0.25f, 0.05f));
    templatedScopedSourceRequest.queryToleranceValid = 1;
    templatedScopedSourceRequest.queryTolerance = 0.05f;
    BObolLodRequest staleTemplateRequest = request;
    staleTemplateRequest.addProviderParam("source_query.space", "world");
    staleTemplateRequest.addProviderParam("source_query.bounds",
					  "9.0 9.0 9.0 10.0 10.0 10.0");
    staleTemplateRequest.addProviderParam("source_query.ray.origin",
					  "0.2 0.2 5.0");
    staleTemplateRequest.addProviderParam("source_query.ray.direction",
					  "0.0 0.0 -1.0");
    staleTemplateRequest.addProviderParam("provider-template", "kept");

    BObolLodRequest cleanTemplatedRequest;
    if (!bobol_lod_rt_source_full_detail_request_from_source_mesh_request(
	    cleanTemplatedRequest, templatedScopedSourceRequest,
	    &staleTemplateRequest) ||
	provider_param_count(cleanTemplatedRequest,
			     "source_query.space") != 1 ||
	provider_param_count(cleanTemplatedRequest,
			     "source_query.bounds") != 1 ||
	provider_param_count(cleanTemplatedRequest,
			     "source_query.tolerance") != 1 ||
	provider_param_count(cleanTemplatedRequest,
			     "source_query.ray.origin") != 0 ||
	provider_param_count(cleanTemplatedRequest,
			     "source_query.ray.direction") != 0 ||
	provider_param_count(cleanTemplatedRequest,
			     "provider-template") != 1) {
	printf("FAIL: LoD RT source full-detail helper should replace stale template query params\n");
	ret = 1;
    } else {
	BObolRtSourceFullDetailProvider templatedLimitProvider;
	templatedLimitProvider.setDatabase(dbip);
	templatedLimitProvider.validateSourceMetrics = TRUE;
	templatedLimitProvider.maxFullDetailFaceCount = 2;
	templatedLimitProvider.maxFullDetailPointCount = 4;
	BObolLodResult cleanTemplatedResult =
	    bobol_rt_source_full_detail_provider_task(
		cleanTemplatedRequest, &templatedLimitProvider);
	if (cleanTemplatedResult.providerStatus !=
	    BOBOL_LOD_PROVIDER_READY ||
	    cleanTemplatedResult.resultKind !=
	    BOBOL_LOD_RESULT_FULL_DETAIL ||
	    cleanTemplatedResult.mesh.coordIndex.size() != 6 ||
	    cleanTemplatedResult.mesh.faceIndex.size() != 2 ||
	    cleanTemplatedResult.counts.faceCount != 2 ||
	    cleanTemplatedResult.counts.pointCount != 4 ||
	    !bobol_lod_result_matches_request(cleanTemplatedResult,
						cleanTemplatedRequest)) {
	    printf("FAIL: LoD RT source full-detail helper should submit current source query params after template cleanup\n");
	    ret = 1;
	}
    }

    BObolLodService submitService;
    if (!submitService.start(1, TRUE)) {
	printf("FAIL: LoD RT source full-detail helper service did not start\n");
	ret = 1;
    } else {
	uint64_t taskId = bobol_lod_submit_rt_source_full_detail_request(
			      &submitService, submitService.beginGeneration(), sourceRequest,
			      dbip, &request, 10, 10);
	if (taskId == 0) {
	    printf("FAIL: LoD RT source full-detail helper did not submit source request\n");
	    ret = 1;
	} else if (wait_for_settled(submitService, 1)) {
	    ret = 1;
	} else {
	    std::vector<BObolLodResult> helperResults;
	    submitService.drainResults(helperResults);
	    if (helperResults.size() != 1 ||
		check_source_full_detail_result(helperResults[0],
						convertedRequest,
						"did not publish helper-submitted BoT mesh result"))
		ret = 1;
	}
	submitService.stop();
    }

    BObolRtSourceFullDetailProvider limitProvider;
    limitProvider.setDatabase(dbip);
    limitProvider.maxFullDetailFaceCount = 3;
    BObolLodResult limitedResult =
	bobol_rt_source_full_detail_provider_task(request, &limitProvider);
    if (limitedResult.providerStatus != BOBOL_LOD_PROVIDER_FALLBACK ||
	bu_strcmp(limitedResult.diagnostic.getString(),
	       "RT source full-detail provider request exceeds full-detail limits") != 0 ||
	limitedResult.mesh.isValid()) {
	printf("FAIL: LoD RT source full-detail provider did not refuse over-budget full detail\n");
	ret = 1;
    }

    BObolLodRequest missingRequest = request;
    missingRequest.objectPath = "/missing.bot";
    missingRequest.objectName = "missing.bot";
    BObolLodResult missingResult =
	bobol_rt_source_full_detail_provider_task(missingRequest, &provider);
    if (missingResult.providerStatus != BOBOL_LOD_PROVIDER_ERROR ||
	bu_strcmp(missingResult.diagnostic.getString(),
	       "RT source full-detail provider could not find source object") != 0 ||
	missingResult.mesh.isValid()) {
	printf("FAIL: LoD RT source full-detail provider did not report missing source object\n");
	ret = 1;
    }

    db_close(dbip);
    bu_file_delete(dbpath);
    return ret;
}

static int
test_rt_mesh_provider_task(void)
{
    char cache_dir[MAXPATHLEN] = {0};
    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;
    std::vector<fastf_t> cachedVertices;
    std::vector<int> cachedFaces;
    std::vector<fastf_t> cachedNormals;

    bu_dir(cache_dir, MAXPATHLEN, BU_DIR_CURR, "bobol_lod_service_cache", NULL);
    bu_dirclear(cache_dir);
    bu_mkdir(cache_dir);
    bu_setenv("BU_DIR_CACHE", cache_dir, 1);

    if (make_provider_test_db(dbpath, sizeof(dbpath), &dbip)) {
	bu_dirclear(cache_dir);
	return 1;
    }

    /* Exceed one private-page target so the service test exercises view
	* unions and page-selective compaction rather than only the legacy-shaped
	* one-prefix case. */
    const int cachedGrid = 190;
    cachedVertices.reserve((cachedGrid + 1) * (cachedGrid + 1) * 3);
    for (int y = 0; y <= cachedGrid; y++) {
	for (int x = 0; x <= cachedGrid; x++) {
	    cachedVertices.push_back((fastf_t)x);
	    cachedVertices.push_back((fastf_t)y);
	    cachedVertices.push_back(0.0);
	}
    }
    cachedFaces.reserve(cachedGrid * cachedGrid * 6);
    for (int y = 0; y < cachedGrid; y++) {
	for (int x = 0; x < cachedGrid; x++) {
	    int v0 = y * (cachedGrid + 1) + x;
	    int v1 = v0 + 1;
	    int v2 = v0 + cachedGrid + 1;
	    int v3 = v2 + 1;
	    cachedFaces.push_back(v0);
	    cachedFaces.push_back(v1);
	    cachedFaces.push_back(v2);
	    cachedFaces.push_back(v1);
	    cachedFaces.push_back(v3);
	    cachedFaces.push_back(v2);
	}
    }
    cachedNormals.reserve(cachedFaces.size() * 3);
    for (size_t i = 0; i < cachedFaces.size(); i++) {
	cachedNormals.push_back(0.0);
	cachedNormals.push_back(0.0);
	cachedNormals.push_back(1.0 + (fastf_t)i);
    }

    struct BObolMeshLodCacheStatus storeStatus =
	    BOBOL_MESH_LOD_CACHE_STATUS_INIT;
    if (bobol_mesh_lod_cache_store_mesh(dbip, "lod-provider.bot",
					  (const point_t *)cachedVertices.data(),
					  cachedVertices.size() / 3,
					  (const vect_t *)cachedNormals.data(), cachedFaces.data(),
					  cachedFaces.size() / 3, 0x12345678, 0,
					  &storeStatus) != BRLCAD_OK ||
	!storeStatus.has_cache_key ||
	!storeStatus.has_cached_payload) {
	printf("FAIL: LoD Obol mesh provider did not store cached mesh normals "
	       "(key=%d payload=%d generated=%d stale=%d)\n",
	       storeStatus.has_cache_key, storeStatus.has_cached_payload,
	       storeStatus.generated_cache_entry, storeStatus.stale_cache_entry);
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    BObolMeshLodProvider provider;
    provider.setDatabase(dbip);
    provider.refreshMissing = TRUE;

    BObolLodTask task;
    task.generation = 1;
    task.request = make_request("/lod-provider.bot");
    task.request.objectName = "lod-provider.bot";
    task.request.providerId = "bobol_mesh_lod";
    task.request.providerVersion = "bobol-cache-v1";
    task.realize = bobol_mesh_lod_provider_task;
    task.realizeData = &provider;

    BObolLodService service;
    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD Obol mesh provider service did not start\n");
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    if (service.submit(task) == 0) {
	printf("FAIL: LoD Obol mesh provider service did not accept task\n");
	service.stop();
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    if (wait_for_settled(service, 1)) {
	service.stop();
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    std::vector<BObolLodResult> results;
    service.drainResults(results);

    int ret = 0;
    if (results.size() != 1 ||
	results[0].resultKind != BOBOL_LOD_RESULT_MESH ||
	results[0].qualityTier != task.request.qualityTier ||
	results[0].providerStatus != BOBOL_LOD_PROVIDER_READY ||
	results[0].geometry.kind != BOBOL_LOD_GEOMETRY_MESH_LOD_CACHE ||
	results[0].geometry.providerToken == 0 ||
	!results[0].geometry.isValid() ||
	results[0].counts.faceCount == 0 ||
	results[0].counts.faceCount > cachedFaces.size() / 3 ||
	results[0].counts.pointCount == 0 ||
	results[0].counts.pointCount > cachedVertices.size() / 3 ||
	!results[0].mesh.isValid() ||
	results[0].mesh.points.size() != results[0].counts.pointCount ||
	results[0].mesh.coordIndex.size() !=
	    results[0].counts.faceCount * 3 ||
	results[0].request.requestedCut != task.request.requestedCut ||
	results[0].resolvedCut < 0 ||
	!bobol_lod_result_matches_request(results[0], task.request)) {
	printf("FAIL: LoD Obol mesh provider task did not return cached mesh result (results=%zu kind=%d quality=%d status=%d geometry=%d token=%llu valid=%d faces=%llu points=%llu mesh_valid=%d mesh_points=%zu indices=%zu match=%d cut=%d requested=%d resolved=%d)\n",
	    results.size(),
	    results.empty() ? -1 : results[0].resultKind,
	    results.empty() ? -1 : results[0].qualityTier,
	    results.empty() ? -1 : results[0].providerStatus,
	    results.empty() ? -1 : results[0].geometry.kind,
	    static_cast<unsigned long long>(results.empty() ? 0 :
		results[0].geometry.providerToken),
	    results.empty() ? 0 : results[0].geometry.isValid(),
	    static_cast<unsigned long long>(results.empty() ? 0 :
		results[0].counts.faceCount),
	    static_cast<unsigned long long>(results.empty() ? 0 :
		results[0].counts.pointCount),
	    results.empty() ? 0 : results[0].mesh.isValid(),
	    results.empty() ? 0 : results[0].mesh.points.size(),
	    results.empty() ? 0 : results[0].mesh.coordIndex.size(),
	    results.empty() ? 0 :
		bobol_lod_result_matches_request(results[0], task.request),
	    results.empty() ? -1 : results[0].geometry.activeCut,
	    task.request.requestedCut,
	    results.empty() ? -1 : results[0].resolvedCut);
	ret = 1;
    }

    /* The retained provider must publish a bounded first prefix and then make
     * monotonic, budget-bounded progress toward the unchanged view target.
     * Cheap adjacent cuts may be collapsed into one delivery; requiring a
     * publication per nominal cut creates a severe many-asset tail. */
    BObolMeshLodProvider stagedProvider;
    stagedProvider.service = &service;
    stagedProvider.setDatabase(dbip);
    stagedProvider.refreshMissing = FALSE;
    stagedProvider.initialRefinementCostBudget = 4;
    stagedProvider.refinementGrowthFactor = 2.0;
    BObolLodRequest stagedRequest = task.request;
    stagedRequest.requestedCut = 4;
    BObolLodResult firstStage =
	service.realizeResidentMeshLod(stagedRequest, stagedProvider);
    const SbBool firstStagePrepared =
	firstStage.progressiveMesh &&
	((firstStage.preparedCadGeometry &&
	  firstStage.preparedCadGeometryRevision ==
	      firstStage.progressiveMesh->revision()) ||
	 (!firstStage.presentationLayers.empty() &&
	  std::all_of(firstStage.presentationLayers.begin(),
	      firstStage.presentationLayers.end(),
	      [](const BObolLodPresentationLayer &layer) {
		  return layer.isValid() ? true : false;
	      }))) ? TRUE : FALSE;
    BObolLodResult secondStage =
	service.realizeResidentMeshLod(stagedRequest, stagedProvider);
    stagedProvider.progressiveDelivery = FALSE;
    BObolLodResult terminalStage =
	service.realizeResidentMeshLod(stagedRequest, stagedProvider);
    /* Keep the terminal prefix resident but tell the provider that the view
     * is currently presenting the first stage.  Refinement must grow from
     * that visible prefix, not jump directly to the resident capacity. */
    stagedProvider.progressiveDelivery = TRUE;
    stagedProvider.useCurrentDrawCut = TRUE;
    stagedProvider.currentDrawCut = firstStage.geometry.activeCut;
    BObolLodResult residentAheadStage =
	service.realizeResidentMeshLod(stagedRequest, stagedProvider);
    BObolMeshLodProvider directProvider;
    directProvider.service = &service;
    directProvider.setDatabase(dbip);
    directProvider.refreshMissing = FALSE;
    directProvider.useCurrentDrawCut = TRUE;
    directProvider.currentDrawCut = firstStage.geometry.activeCut;
    BObolLodResult directStage =
	service.realizeResidentMeshLod(stagedRequest, directProvider);
    /* Zoom residency and presentation are separate contracts.  A worker may
     * make the complete pixel-demanded prefix resident while publishing the
     * previously affordable draw cut; a later frame can retarget the same
     * immutable arrays without another cache read. */
    BObolMeshLodProvider prefetchedProvider = directProvider;
    prefetchedProvider.progressiveDelivery = FALSE;
    prefetchedProvider.useDeliveryCutLimit = TRUE;
    prefetchedProvider.deliveryCutLimit = stagedRequest.requestedCut;
    prefetchedProvider.usePresentationCutLimit = TRUE;
    prefetchedProvider.presentationCutLimit =
	firstStage.geometry.activeCut;
    BObolLodResult prefetchedStage =
	service.realizeResidentMeshLod(stagedRequest, prefetchedProvider);
    if (firstStage.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	firstStage.geometry.activeCut < 0 ||
	firstStage.geometry.activeCut >= stagedRequest.requestedCut ||
	firstStage.terminal ||
	secondStage.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	secondStage.geometry.activeCut <= firstStage.geometry.activeCut ||
	secondStage.geometry.activeCut != firstStage.geometry.activeCut + 1 ||
	secondStage.geometry.activeCut >= stagedRequest.requestedCut ||
	secondStage.terminal ||
	terminalStage.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	terminalStage.geometry.activeCut != stagedRequest.requestedCut ||
	!terminalStage.terminal ||
	residentAheadStage.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	residentAheadStage.geometry.activeCut <=
	    firstStage.geometry.activeCut ||
	residentAheadStage.geometry.activeCut !=
	    firstStage.geometry.activeCut + 1 ||
	residentAheadStage.geometry.activeCut >=
	    stagedRequest.requestedCut ||
	residentAheadStage.terminal ||
	directStage.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	directStage.geometry.activeCut != stagedRequest.requestedCut ||
	!directStage.terminal ||
	prefetchedStage.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	prefetchedStage.geometry.activeCut !=
	    firstStage.geometry.activeCut ||
	prefetchedStage.residentCut != stagedRequest.requestedCut ||
	prefetchedStage.terminal ||
	!firstStage.progressiveMesh ||
	!firstStagePrepared ||
	firstStage.progressiveMesh != secondStage.progressiveMesh ||
	firstStage.progressiveMesh != terminalStage.progressiveMesh ||
	firstStage.progressiveMesh != residentAheadStage.progressiveMesh ||
	firstStage.progressiveMesh != directStage.progressiveMesh ||
	firstStage.counts.faceCount > secondStage.counts.faceCount ||
	secondStage.counts.faceCount > terminalStage.counts.faceCount) {
	printf("FAIL: retained LoD provider did not stage large growth or "
	    "collapse a cheap target (levels=%d/%d/%d resident-ahead=%d "
	    "direct=%d prefetched=%d/%d terminal=%d/%d/%d/%d/%d/%d "
	    "prepared=%d "
	    "faces=%llu/%llu/%llu/%llu/%llu diagnostic=%s)\n",
	    firstStage.geometry.activeCut, secondStage.geometry.activeCut,
	    terminalStage.geometry.activeCut,
	    residentAheadStage.geometry.activeCut,
	    directStage.geometry.activeCut,
	    prefetchedStage.geometry.activeCut,
	    prefetchedStage.residentCut,
	    firstStage.terminal ? 1 : 0,
	    secondStage.terminal ? 1 : 0, terminalStage.terminal ? 1 : 0,
	    residentAheadStage.terminal ? 1 : 0,
	    directStage.terminal ? 1 : 0,
	    prefetchedStage.terminal ? 1 : 0,
	    firstStagePrepared ? 1 : 0,
	    static_cast<unsigned long long>(firstStage.counts.faceCount),
	    static_cast<unsigned long long>(secondStage.counts.faceCount),
	    static_cast<unsigned long long>(terminalStage.counts.faceCount),
	    static_cast<unsigned long long>(
		residentAheadStage.counts.faceCount),
	    static_cast<unsigned long long>(directStage.counts.faceCount),
	    firstStage.diagnostic.getString());
	ret = 1;
    }

    /* Two views of one large logical leaf retain the union of their private
	 * pages.  Establish both snapshots before starting this isolated worker so
	 * the first view's now-stale trim is also exercised deterministically. */
    {
	struct BObolMeshLod *unionLod =
	    bobol_mesh_lod_get(dbip, "lod-provider.bot");
	BObolMeshLodHierarchyInfo unionHierarchy =
	    BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
	BObolLodService unionService;
	BObolMeshLodProvider unionProvider;
	unionProvider.service = &unionService;
	unionProvider.setDatabase(dbip);
	unionProvider.refreshMissing = FALSE;
	unionProvider.progressiveDelivery = FALSE;
	unionProvider.atomicRepresentationHandoff = TRUE;
	unionProvider.shrinkAfterCopy = TRUE;
	BObolLodRequest unionRequest = stagedRequest;
	if (!unionLod ||
	    !bobol_mesh_lod_hierarchy_info_get(unionLod, &unionHierarchy) ||
	    unionHierarchy.chunk_count < 2) {
	    printf("FAIL: resident page-union fixture is not chunked\n");
	    ret = 1;
	} else {
	    uint32_t latePage = UINT32_MAX;
	    int latePageFirstCut = -1;
	    for (uint32_t chunk = 0; chunk < unionHierarchy.chunk_count;
		    ++chunk) {
		for (int cut = unionHierarchy.min_cut;
			cut <= unionHierarchy.max_cut; ++cut) {
		    if (!unionHierarchy.chunks[chunk].cuts[cut].face_count)
			continue;
		    if (cut > unionHierarchy.min_cut) {
			latePage = chunk;
			latePageFirstCut = cut;
		    }
		    break;
		}
		if (latePage != UINT32_MAX)
		    break;
	    }
	    if (latePage == UINT32_MAX) {
		printf("FAIL: resident page-union fixture has no delayed page "
		       "population\n");
		ret = 1;
	    } else {
		BObolLodCounts firstPageCounts;
		firstPageCounts.faceCount = unionHierarchy.chunks[latePage].
		    cuts[latePageFirstCut].face_count;
		firstPageCounts.pointCount = unionHierarchy.chunks[latePage].
		    cuts[latePageFirstCut].point_count;
		firstPageCounts.originalPointCount =
		    firstPageCounts.pointCount;
		firstPageCounts.normalCount = unionHierarchy.has_normals ?
		    firstPageCounts.pointCount : 0;
		BObolLodService latePageService;
		BObolMeshLodProvider latePageProvider;
		latePageProvider.service = &latePageService;
		latePageProvider.setDatabase(dbip);
		latePageProvider.refreshMissing = FALSE;
		latePageProvider.initialRefinementCostBudget =
		    bobol_lod_render_cost_units(
			firstPageCounts, BOBOL_LOD_DRAW_SHADED, 1);
		BObolLodRequest latePageRequest = stagedRequest;
		latePageRequest.requestedCut = unionHierarchy.max_cut;
		latePageRequest.requiredChunks = {latePage};
		latePageRequest.occurrenceKey = "late-page-occurrence";
		BObolLodResult latePageResult =
		    latePageService.realizeResidentMeshLod(
			latePageRequest, latePageProvider);
		if (latePageResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
		    latePageResult.geometry.activeCut < latePageFirstCut ||
		    latePageResult.counts.faceCount == 0 ||
		    latePageResult.presentationLayers.empty()) {
		    printf("FAIL: delayed visible page population became a "
			   "terminal empty delivery (page=%u first=%d active=%d "
			   "status=%d faces=%llu layers=%zu diagnostic=%s)\n",
			   latePage, latePageFirstCut,
			   latePageResult.geometry.activeCut,
			   latePageResult.providerStatus,
			   static_cast<unsigned long long>(
			       latePageResult.counts.faceCount),
			   latePageResult.presentationLayers.size(),
			   latePageResult.diagnostic.getString());
		    ret = 1;
		}
		latePageProvider.initialRefinementCostBudget = 1;
		latePageProvider.resetExisting = TRUE;
		BObolLodResult indivisiblePageResult =
		    latePageService.realizeResidentMeshLod(
			latePageRequest, latePageProvider);
		if (indivisiblePageResult.providerStatus !=
			BOBOL_LOD_PROVIDER_READY ||
		    indivisiblePageResult.geometry.activeCut <
			latePageFirstCut ||
		    indivisiblePageResult.counts.faceCount == 0 ||
		    indivisiblePageResult.presentationLayers.empty()) {
		    printf("FAIL: an underestimated first-page reservation "
			   "published an empty or failed spatial delivery "
			   "(page=%u first=%d active=%d status=%d faces=%llu "
			   "layers=%zu diagnostic=%s)\n",
			   latePage, latePageFirstCut,
			   indivisiblePageResult.geometry.activeCut,
			   indivisiblePageResult.providerStatus,
			   static_cast<unsigned long long>(
			       indivisiblePageResult.counts.faceCount),
			   indivisiblePageResult.presentationLayers.size(),
			   indivisiblePageResult.diagnostic.getString());
		    ret = 1;
		}
		BObolLodRequest emptyPageRequest = latePageRequest;
		emptyPageRequest.requestedCut = latePageFirstCut - 1;
		BObolLodResult emptyPageResult =
		    latePageService.realizeResidentMeshLod(
			emptyPageRequest, latePageProvider);
		if (emptyPageResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
		    !emptyPageResult.progressiveMesh ||
		    !emptyPageResult.progressiveMesh->isValid() ||
		    emptyPageResult.geometry.activeCut !=
			emptyPageRequest.requestedCut ||
		    emptyPageResult.counts.faceCount != 0 ||
		    !emptyPageResult.presentationLayers.empty() ||
		    !emptyPageResult.payloadIsConsistent() ||
		    !emptyPageResult.terminal) {
		    printf("FAIL: a valid empty spatial population was not "
			   "published as a terminal zero-draw mesh "
			   "(requested=%d active=%d status=%d faces=%llu "
			   "layers=%zu terminal=%d diagnostic=%s)\n",
			   emptyPageRequest.requestedCut,
			   emptyPageResult.geometry.activeCut,
			   emptyPageResult.providerStatus,
			   static_cast<unsigned long long>(
			       emptyPageResult.counts.faceCount),
			   emptyPageResult.presentationLayers.size(),
			   emptyPageResult.terminal ? 1 : 0,
			   emptyPageResult.diagnostic.getString());
		    ret = 1;
		}
	    }
	    unionRequest.requestedCut = unionHierarchy.max_cut;
	    unionRequest.requiredChunks.resize(unionHierarchy.chunk_count);
	    std::iota(unionRequest.requiredChunks.begin(),
		unionRequest.requiredChunks.end(), 0u);
	    BObolLodResult unionRich = unionService.realizeResidentMeshLod(
		unionRequest, unionProvider);
	    std::vector<uint32_t> initiallyResident;
	    if (unionRich.progressiveMesh)
		unionRich.progressiveMesh->residentChunkIds(initiallyResident);
	    BObolLodResidentDemand firstView;
	    firstView.assetKey = unionRich.geometry.cacheKey.value;
	    firstView.cut = unionHierarchy.max_cut;
	    firstView.channelMask = 2;
	    firstView.chunkIds = {0};
	    BObolLodResidentDemand secondView = firstView;
	    secondView.chunkIds = {1};
	    secondView.cut = std::max(
		unionHierarchy.chunks[1].min_cut,
		unionHierarchy.max_cut - 1);
	    BObolLodResidentDemand allView = firstView;
	    allView.chunkIds.resize(unionHierarchy.chunk_count);
	    std::iota(allView.chunkIds.begin(), allView.chunkIds.end(), 0u);
	    const SbBool started = unionService.start(1, FALSE);
	    const size_t allQueued =
		unionService.scheduleResidentMeshCompaction(
		    0xa002, 1, {allView});
	    const size_t firstQueued =
		unionService.scheduleResidentMeshCompaction(
		    0xa001, 1, {firstView});
	    const size_t secondQueued =
		unionService.scheduleResidentMeshCompaction(
		    0xa002, 2, {secondView});
	    const int unionWait = started ?
		wait_for_resident_compaction(unionService) : 1;
	    std::vector<uint32_t> retained;
	    if (unionRich.progressiveMesh)
		unionRich.progressiveMesh->residentChunkIds(retained);
	    std::vector<BObolLodChunkCut> retainedCuts;
	    if (unionRich.progressiveMesh)
		unionRich.progressiveMesh->residentChunkCuts(retainedCuts);
	    const bool retainedBoth =
		std::binary_search(retained.begin(), retained.end(), 0u) &&
		std::binary_search(retained.begin(), retained.end(), 1u);
	    const bool retainedIndependentCuts = retainedCuts.size() == 2 &&
		retainedCuts[0] == BObolLodChunkCut{0, unionHierarchy.max_cut} &&
		retainedCuts[1] == BObolLodChunkCut{1, secondView.cut};
	    const std::vector<uint32_t> unionRetained = retained;
	    unionService.releaseResidentMeshConsumer(0xa002);
	    const size_t trimQueued =
		unionService.scheduleResidentMeshCompaction(
		    0xa001, 2, {firstView});
	    const int trimWait = started ?
		wait_for_resident_compaction(unionService) : 1;
	    retained.clear();
	    if (unionRich.progressiveMesh)
		unionRich.progressiveMesh->residentChunkIds(retained);
	    std::vector<BObolLodChunkCut> coverageFloorCuts;
	    if (unionRich.progressiveMesh)
		unionRich.progressiveMesh->residentChunkCuts(coverageFloorCuts);
	    const bool retainedWarmPages =
		coverageFloorCuts.size() == unionHierarchy.chunk_count &&
		coverageFloorCuts[0] ==
		    BObolLodChunkCut{0, unionHierarchy.max_cut} &&
		coverageFloorCuts[1] ==
		    BObolLodChunkCut{1, secondView.cut};
	    /* Expanding the visible page set while zooming out must be an atomic
	     * transition from the incumbent presentation to the new view target.
	     * A stale bounded-delivery stamp may ask for the absolute minimum, but
	     * publishing that page population would make the whole occurrence flash
	     * as a few coarse slabs. */
	    BObolMeshLodProvider continuityProvider = unionProvider;
	    continuityProvider.atomicRepresentationHandoff = FALSE;
	    continuityProvider.useCurrentDrawCut = TRUE;
	    continuityProvider.currentDrawCut = unionHierarchy.max_cut;
	    continuityProvider.useDeliveryCutLimit = TRUE;
	    continuityProvider.deliveryCutLimit = unionHierarchy.min_cut;
	    continuityProvider.usePresentationCutLimit = TRUE;
	    continuityProvider.presentationCutLimit = unionHierarchy.min_cut;
	    BObolLodRequest continuityRequest = unionRequest;
	    continuityRequest.requiredChunks = {0, 1};
	    continuityRequest.requestedCut = std::max(
		unionHierarchy.min_cut, unionHierarchy.max_cut - 1);
	    BObolLodResult continuityResult =
		unionService.realizeResidentMeshLod(
		    continuityRequest, continuityProvider);
	    BObolLodResult restoredUnion =
		unionService.realizeResidentMeshLod(
		    unionRequest, unionProvider);
	    std::vector<uint32_t> restored;
	    if (restoredUnion.progressiveMesh)
		restoredUnion.progressiveMesh->residentChunkIds(restored);
	    BObolLodResidentDemand incompleteView = allView;
	    incompleteView.chunkIds.push_back(unionHierarchy.chunk_count);
	    const int beforeIncompleteCut = unionRich.progressiveMesh ?
		unionRich.progressiveMesh->residentCut() : -1;
	    const size_t incompleteQueued =
		unionService.scheduleResidentMeshCompaction(
		    0xa001, 3, {incompleteView});
	    const int incompleteWait = started ?
		wait_for_resident_compaction(unionService) : 1;
	    const bool incompletePreserved = unionRich.progressiveMesh &&
		unionRich.progressiveMesh->residentCut() == beforeIncompleteCut;
	    (void)unionService.scheduleResidentMeshCompaction(
		0xa001, 4, {firstView});
	    /* Healthy memory above deliberately preserves page 1 as a latency
	     * cache.  Once an explicit resident limit is exceeded, the same demand
	     * may compact that off-screen page to the documented coverage floor
	     * while preserving the demanded page at full quality. */
	    int pressureCoverageCut = std::max(
		unionHierarchy.min_cut,
		unionRich.progressiveMesh->cutForScreenError(128.0, 1.0));
	    pressureCoverageCut = std::min(
		unionHierarchy.max_cut, pressureCoverageCut);
	    while (pressureCoverageCut < unionHierarchy.max_cut &&
		!unionHierarchy.chunks[1].cuts[pressureCoverageCut].face_count)
		pressureCoverageCut++;
	    const int expectedPressurePageCut = std::min(
		unionHierarchy.max_cut, pressureCoverageCut);
	    unionService.setResidentMeshLimit(1);
	    const size_t pressureQueued =
		unionService.scheduleResidentMeshCompaction(
		    0xa001, 5, {firstView});
	    const int pressureWait = started ?
		wait_for_resident_compaction(unionService) : 1;
	    std::vector<BObolLodChunkCut> pressureCuts;
	    if (unionRich.progressiveMesh)
		unionRich.progressiveMesh->residentChunkCuts(pressureCuts);
	    const bool pressureTrimmed =
		pressureCuts.size() == unionHierarchy.chunk_count &&
		pressureCuts[0] ==
		    BObolLodChunkCut{0, unionHierarchy.max_cut} &&
		pressureCuts[1] ==
		    BObolLodChunkCut{1, expectedPressurePageCut};
	    if (unionRich.providerStatus != BOBOL_LOD_PROVIDER_READY ||
		!unionRich.progressiveMesh || allQueued || firstQueued ||
		!secondQueued || !started || unionWait || !retainedBoth ||
		!retainedIndependentCuts ||
		trimQueued ||
		trimWait || !retainedWarmPages ||
		continuityResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
		continuityResult.geometry.activeCut !=
		    continuityRequest.requestedCut ||
		continuityResult.residentCut <
		    continuityRequest.requestedCut ||
		restoredUnion.providerStatus != BOBOL_LOD_PROVIDER_READY ||
		restoredUnion.progressiveMesh != unionRich.progressiveMesh ||
		restored.size() != unionHierarchy.chunk_count ||
		incompleteQueued != 0 || incompleteWait ||
		!incompletePreserved ||
		pressureQueued !=
		    (expectedPressurePageCut < unionHierarchy.max_cut ? 1u : 0u) ||
		pressureWait || !pressureTrimmed) {
		printf("FAIL: multi-view resident page union/trim/continuity/restore "
		       "(queued=%zu/%zu/%zu trim=%zu pressure=%zu "
		       "incomplete=%zu/%d/%d started=%d wait=%d/%d/%d "
		       "both=%d cuts=%d first=%d continuity=%d/%d restored=%zu/%u "
		       "pressure-cut=%d/%d status=%d/%d/%d)\n",
		       allQueued, firstQueued, secondQueued, trimQueued,
		       pressureQueued, incompleteQueued, incompleteWait,
		       incompletePreserved ? 1 : 0,
		       started ? 1 : 0, unionWait, trimWait,
		       pressureWait,
		       retainedBoth ? 1 : 0,
		       retainedIndependentCuts ? 1 : 0,
		       retainedWarmPages ? 1 : 0,
		       continuityResult.geometry.activeCut,
		       continuityRequest.requestedCut,
		       restored.size(), unionHierarchy.chunk_count,
		       pressureCuts.size() > 1 ? pressureCuts[1].cut : -1,
		       expectedPressurePageCut,
		       unionRich.providerStatus,
		       continuityResult.providerStatus,
		       restoredUnion.providerStatus);
		printf("  pages initial:");
		for (uint32_t page : initiallyResident)
		    printf(" %u", page);
		printf("\n  pages union:");
		for (uint32_t page : unionRetained)
		    printf(" %u", page);
		printf("\n  pages trimmed:");
		for (uint32_t page : retained)
		    printf(" %u", page);
		printf("\n  trimmed cuts:");
		for (const BObolLodChunkCut &entry : coverageFloorCuts)
		    printf(" %u:%d", entry.chunkId, entry.cut);
		printf("\n  pages restored:");
		for (uint32_t page : restored)
		    printf(" %u", page);
		printf("\n  pressure cuts:");
		for (const BObolLodChunkCut &entry : pressureCuts)
		    printf(" %u:%d", entry.chunkId, entry.cut);
		printf("\n");
		ret = 1;
	    }
	    unionService.stop();
	}
	if (unionLod)
	    bobol_mesh_lod_destroy(unionLod);
    }

    /* A new view reopening a complete persistent hierarchy must publish the
     * exact prefix selected by its bounded delivery policy and discard the
     * cache reader's duplicate arrays. */
    BObolLodService warmFirstService;
    BObolMeshLodProvider warmFirstProvider;
    warmFirstProvider.service = &warmFirstService;
    warmFirstProvider.setDatabase(dbip);
    warmFirstProvider.refreshMissing = FALSE;
    warmFirstProvider.shrinkAfterCopy = TRUE;
    BObolLodResult warmFirst = warmFirstService.realizeResidentMeshLod(
	stagedRequest, warmFirstProvider);
    SbBool warmFirstPlanningComplete = FALSE;
    size_t warmFirstMaintenance = 0;
    if (warmFirst.progressiveMesh && warmFirst.geometry.cacheKey.isValid()) {
	BObolLodResidentDemand demand;
	demand.assetKey = warmFirst.geometry.cacheKey.value;
	demand.cut = warmFirst.residentCut;
	demand.channelMask = 2;
	demand.chunkIds = warmFirst.request.requiredChunks;
	std::vector<BObolLodResidentDemand> demands(1, demand);
	warmFirstMaintenance = warmFirstService.scheduleResidentMeshCompaction(
	    0x9876, 1, demands, &warmFirstPlanningComplete);
    }
    if (warmFirst.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	!warmFirst.progressiveMesh ||
	warmFirst.geometry.activeCut != stagedRequest.requestedCut ||
	warmFirst.residentCut != stagedRequest.requestedCut ||
	!warmFirst.terminal || !warmFirstPlanningComplete ||
	warmFirstMaintenance != 0) {
	printf("FAIL: warm first publication did not load once and release "
	       "duplicate backing (status=%d active=%d min=%d resident=%d "
	       "requested=%d terminal=%d planned=%d maintenance=%zu)\n",
	       warmFirst.providerStatus, warmFirst.geometry.activeCut,
	       warmFirst.progressiveMesh ?
		   warmFirst.progressiveMesh->minimumCut() : -1,
	       warmFirst.residentCut, stagedRequest.requestedCut,
	       warmFirst.terminal ? 1 : 0,
	       warmFirstPlanningComplete ? 1 : 0,
	       warmFirstMaintenance);
	ret = 1;
    }
    warmFirstService.releaseResidentMeshConsumer(0x9876);

    /* Replacing one adaptive representation with another is an atomic
     * handoff.  Unlike a cold first publication, a replacement must not put
     * its minimum prefix on screen while the requested prefix is already
     * available in the persistent hierarchy. */
    BObolLodService atomicHandoffService;
    BObolMeshLodProvider atomicHandoffProvider;
    atomicHandoffProvider.service = &atomicHandoffService;
    atomicHandoffProvider.setDatabase(dbip);
    atomicHandoffProvider.refreshMissing = FALSE;
    atomicHandoffProvider.atomicRepresentationHandoff = TRUE;
    BObolLodResult atomicHandoff =
	atomicHandoffService.realizeResidentMeshLod(
	    stagedRequest, atomicHandoffProvider);
    if (atomicHandoff.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	!atomicHandoff.progressiveMesh ||
	atomicHandoff.geometry.activeCut != stagedRequest.requestedCut ||
	atomicHandoff.residentCut != stagedRequest.requestedCut ||
	!atomicHandoff.terminal) {
	printf("FAIL: adaptive representation handoff published a transient "
	       "minimum prefix (status=%d active=%d resident=%d requested=%d "
	       "terminal=%d)\n",
	       atomicHandoff.providerStatus,
	       atomicHandoff.geometry.activeCut,
	       atomicHandoff.residentCut, stagedRequest.requestedCut,
	       atomicHandoff.terminal ? 1 : 0);
	ret = 1;
    }

    /* A stable demand snapshot protects the presented prefix.  Once the
     * asset disappears from every consumer snapshot, the service must keep
     * only its minimum useful prefix and later restore richer data from the
     * same on-disk cache into the same retained mesh identity. */
    if (terminalStage.progressiveMesh &&
	terminalStage.geometry.cacheKey.isValid() &&
	terminalStage.progressiveMesh->residentCut() >
	    terminalStage.progressiveMesh->minimumCut()) {
	BObolLodResidentDemand demand;
	demand.assetKey = terminalStage.geometry.cacheKey.value;
	demand.cut = terminalStage.geometry.activeCut;
	demand.channelMask = 2;
	demand.chunkIds = terminalStage.request.requiredChunks;
	std::vector<BObolLodResidentDemand> demanded(1, demand);
	const int richLevel = terminalStage.progressiveMesh->residentCut();
	std::vector<BObolLodChunkCut> richChunkCuts;
	terminalStage.progressiveMesh->residentChunkCuts(richChunkCuts);
	const size_t richBytes =
	    service.residentMeshBytesForDiagnostics();
	const size_t maintenanceQueued =
	    service.scheduleResidentMeshCompaction(
		0x1234, 1, demanded);
	uint64_t maintenancePlanRevision = 0;
	size_t maintenanceCandidates = SIZE_MAX;
	const SbBool maintenancePlanCurrent =
	    service.residentMeshCompactionPlanForDiagnostics(
		0x1234, &maintenancePlanRevision, &maintenanceCandidates);
	const int maintenanceWait =
	    wait_for_resident_compaction(service);
	std::vector<BObolLodResidentCompaction> maintenanceResults;
	service.drainResidentMeshCompactions(
	    0x1234, maintenanceResults);
	const size_t maintainedBytes =
	    service.residentMeshBytesForDiagnostics();
	if (maintenanceQueued != 0 || !maintenancePlanCurrent ||
	    maintenancePlanRevision != 1 || maintenanceCandidates != 0 ||
	    maintenanceWait ||
	    !maintenanceResults.empty() ||
	    terminalStage.progressiveMesh->residentCut() != richLevel ||
	    maintainedBytes != richBytes) {
	    printf("FAIL: stable resident PoP demand queued redundant backing "
		   "cleanup (queued=%zu level=%d/%d "
		   "bytes=%zu/%zu results=%zu plan=%d/%llu candidates=%zu)\n",
		   maintenanceQueued,
		   terminalStage.progressiveMesh->residentCut(), richLevel,
		   maintainedBytes, richBytes, maintenanceResults.size(),
		   maintenancePlanCurrent ? 1 : 0,
		   static_cast<unsigned long long>(maintenancePlanRevision),
		   maintenanceCandidates);
	    ret = 1;
	}

	std::vector<BObolLodResidentDemand> hidden;
	const size_t queued =
	    service.scheduleResidentMeshCompaction(0x1234, 2, hidden);
	const int compactWait = wait_for_resident_compaction(service);
	const size_t compactBytes =
	    service.residentMeshBytesForDiagnostics();
	std::vector<uint32_t> coverageFloorPages;
	terminalStage.progressiveMesh->residentChunkIds(coverageFloorPages);
	std::vector<BObolLodChunkCut> coverageFloorCuts;
	terminalStage.progressiveMesh->residentChunkCuts(coverageFloorCuts);
	const int coverageFloorLevel =
	    terminalStage.progressiveMesh->residentCut();
	if (queued != 1 || compactWait ||
	    coverageFloorPages != terminalStage.request.requiredChunks ||
	    coverageFloorCuts == richChunkCuts ||
	    compactBytes >= richBytes) {
	    printf("FAIL: hidden resident PoP prefix did not compact "
		   "(count=%zu level=%d min=%d pages=%zu/%zu bytes=%zu/%zu)\n",
		   queued,
		   coverageFloorLevel,
		   terminalStage.progressiveMesh->minimumCut(),
		   coverageFloorPages.size(),
		   terminalStage.request.requiredChunks.size(),
		   compactBytes, richBytes);
	    ret = 1;
	}

	const uint64_t loadsBefore =
	    service.residentMeshCacheLoadCountForDiagnostics();
	BObolMeshLodProvider reloadProvider = stagedProvider;
	reloadProvider.progressiveDelivery = FALSE;
	reloadProvider.useCurrentDrawCut = FALSE;
	BObolLodResult reloaded =
	    service.realizeResidentMeshLod(stagedRequest, reloadProvider);
	if (reloaded.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	    !reloaded.terminal ||
	    reloaded.geometry.activeCut != stagedRequest.requestedCut ||
	    reloaded.progressiveMesh != terminalStage.progressiveMesh ||
	    reloaded.progressiveMesh->residentCut() !=
		stagedRequest.requestedCut ||
	    service.residentMeshCacheLoadCountForDiagnostics() <=
		loadsBefore) {
	    printf("FAIL: compacted resident PoP prefix did not reload "
		   "from cache in place (status=%d terminal=%d active=%d "
		   "resident=%d requested=%d same=%d loads=%llu/%llu)\n",
		   reloaded.providerStatus, reloaded.terminal ? 1 : 0,
		   reloaded.geometry.activeCut,
		   reloaded.progressiveMesh ?
		       reloaded.progressiveMesh->residentCut() : -1,
		   stagedRequest.requestedCut,
		   reloaded.progressiveMesh == terminalStage.progressiveMesh ?
		       1 : 0,
		   static_cast<unsigned long long>(
		       service.residentMeshCacheLoadCountForDiagnostics()),
		       static_cast<unsigned long long>(loadsBefore));
	    ret = 1;
	}

	/* A quiet trim is planned from a complete old demand snapshot.  Renewed
	 * use may arrive before a worker takes that plan; it must invalidate the
	 * trim rather than let the old worker shorten the same shared progressive
	 * object after the richer request has already returned.  Queue before
	 * starting this isolated service to make that ordering deterministic. */
	BObolLodService staleCompactionService;
	BObolMeshLodProvider staleProvider = reloadProvider;
	staleProvider.service = &staleCompactionService;
	BObolLodResult staleMinimum =
	    staleCompactionService.realizeResidentMeshLod(
		stagedRequest, staleProvider);
	BObolLodResult staleRich =
	    staleCompactionService.realizeResidentMeshLod(
		stagedRequest, staleProvider);
	size_t staleQueued = 0;
	BObolLodResult renewedRich;
	if (staleRich.providerStatus == BOBOL_LOD_PROVIDER_READY &&
	    staleRich.progressiveMesh &&
	    staleRich.progressiveMesh->residentCut() >
		staleRich.progressiveMesh->minimumCut()) {
	    BObolLodResidentDemand staleDemand;
	    staleDemand.assetKey = staleRich.geometry.cacheKey.value;
	    staleDemand.cut = staleRich.progressiveMesh->minimumCut();
	    staleDemand.channelMask = 2;
	    staleDemand.chunkIds = staleRich.request.requiredChunks;
	    std::vector<BObolLodResidentDemand> staleDemands(1, staleDemand);
	    staleQueued = staleCompactionService.scheduleResidentMeshCompaction(
		0x5678, 1, staleDemands);
	    renewedRich = staleCompactionService.realizeResidentMeshLod(
		stagedRequest, staleProvider);
	}
	const int staleRichLevel = renewedRich.progressiveMesh ?
	    renewedRich.progressiveMesh->residentCut() : -1;
	const SbBool staleStarted =
	    staleCompactionService.start(1, FALSE);
	const int staleWait = staleStarted ?
	    wait_for_resident_compaction(staleCompactionService) : 1;
	std::vector<BObolLodResidentCompaction> staleCompletions;
	staleCompactionService.drainResidentMeshCompactions(
	    0x5678, staleCompletions);
	if (staleMinimum.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	    staleRich.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	    !staleRich.progressiveMesh || staleQueued != 1 ||
	    renewedRich.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	    renewedRich.progressiveMesh != staleRich.progressiveMesh ||
	    !staleStarted || staleWait || !staleCompletions.empty() ||
	    staleRich.progressiveMesh->residentCut() != staleRichLevel ||
	    staleRichLevel != stagedRequest.requestedCut) {
	    printf("FAIL: renewed resident use did not cancel a stale queued "
		   "compaction (queued=%zu started=%d wait=%d status=%d/%d "
		   "level=%d/%d results=%zu)\n",
		   staleQueued, staleStarted ? 1 : 0, staleWait,
		   staleRich.providerStatus, renewedRich.providerStatus,
		   staleRich.progressiveMesh ?
		       staleRich.progressiveMesh->residentCut() : -1,
		   staleRichLevel, staleCompletions.size());
	    ret = 1;
	}
	staleCompactionService.stop();

	/* A consumer-wide semantic edge must cancel a queued trim even when no
	 * provider request touches the asset.  Retained allocation can promote an
	 * already-resident cut without I/O, so per-asset use notification alone is
	 * not a complete stale-work guard. */
	BObolLodService invalidatedCompactionService;
	BObolMeshLodProvider invalidatedProvider = reloadProvider;
	invalidatedProvider.service = &invalidatedCompactionService;
	BObolLodResult invalidatedRich =
	    invalidatedCompactionService.realizeResidentMeshLod(
		stagedRequest, invalidatedProvider);
	size_t invalidatedQueued = 0;
	int invalidatedRichCut = -1;
	if (invalidatedRich.providerStatus == BOBOL_LOD_PROVIDER_READY &&
	    invalidatedRich.progressiveMesh) {
	    invalidatedRichCut = invalidatedRich.progressiveMesh->residentCut();
	    BObolLodResidentDemand oldDemand;
	    oldDemand.assetKey = invalidatedRich.geometry.cacheKey.value;
	    oldDemand.cut = invalidatedRich.progressiveMesh->minimumCut();
	    oldDemand.channelMask = 2;
	    oldDemand.chunkIds = invalidatedRich.request.requiredChunks;
	    invalidatedQueued =
		invalidatedCompactionService.scheduleResidentMeshCompaction(
		    0x5679, 17, {oldDemand});
	    invalidatedCompactionService.invalidateResidentMeshConsumer(
		0x5679);
	}
	const SbBool invalidatedStarted =
	    invalidatedCompactionService.start(1, FALSE);
	const int invalidatedWait = invalidatedStarted ?
	    wait_for_resident_compaction(invalidatedCompactionService) : 1;
	std::vector<BObolLodResidentCompaction> invalidatedCompletions;
	invalidatedCompactionService.drainResidentMeshCompactions(
	    0x5679, invalidatedCompletions);
	if (invalidatedRich.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	    !invalidatedRich.progressiveMesh || invalidatedQueued != 1 ||
	    !invalidatedStarted || invalidatedWait ||
	    !invalidatedCompletions.empty() ||
	    invalidatedRich.progressiveMesh->residentCut() !=
		invalidatedRichCut) {
	    printf("FAIL: consumer demand invalidation did not cancel a queued "
		   "resident compaction (queued=%zu started=%d wait=%d "
		   "cut=%d/%d results=%zu)\n",
		   invalidatedQueued, invalidatedStarted ? 1 : 0,
		   invalidatedWait,
		   invalidatedRich.progressiveMesh ?
		       invalidatedRich.progressiveMesh->residentCut() : -1,
		   invalidatedRichCut, invalidatedCompletions.size());
	    ret = 1;
	}
	invalidatedCompactionService.stop();

	if (reloaded.providerStatus == BOBOL_LOD_PROVIDER_READY &&
	    reloaded.progressiveMesh &&
	    reloaded.progressiveMesh->residentCut() >
		reloaded.progressiveMesh->minimumCut()) {
	    BObolLodResidentDemand compactDemand;
	    compactDemand.assetKey = reloaded.geometry.cacheKey.value;
	    compactDemand.cut =
		reloaded.progressiveMesh->minimumCut();
	    compactDemand.channelMask = 2;
	    compactDemand.chunkIds = reloaded.request.requiredChunks;
	    std::vector<BObolLodResidentDemand> compactDemands(
		1, compactDemand);
	    const int rendererRichCut = reloaded.progressiveMesh->residentCut();
	    const size_t rendererQueued =
		service.scheduleResidentMeshCompaction(
		    0x1234, 3, compactDemands);
	    const int waitResult =
		wait_for_resident_compaction(service);
	    std::vector<BObolLodResidentCompaction> completions;
	    service.drainResidentMeshCompactions(
		0x1234, completions);
	    std::vector<uint32_t> rendererFloorPages;
	    reloaded.progressiveMesh->residentChunkIds(rendererFloorPages);
	    if (rendererQueued != 1 || waitResult ||
		completions.size() != 1 ||
		completions[0].progressiveMesh != reloaded.progressiveMesh ||
		completions[0].residentCut !=
		    reloaded.progressiveMesh->residentCut() ||
		completions[0].residentCut > rendererRichCut ||
		rendererFloorPages != compactDemand.chunkIds ||
		completions[0].channelMask != 2 ||
		completions[0].preparedCadGeometry ||
		!completions[0].presentationLayers.empty()) {
		printf("FAIL: spatial resident compaction published "
		       "view-specific presentation geometry\n");
		ret = 1;
	    }
	}

	/*
	 * A hidden asset normally keeps its minimum prefix warm.  Under an
	 * explicit retained-memory limit, however, an asset absent from every
	 * consumer snapshot must be released completely and recreated from the
	 * persistent PoP cache if a later view needs it again.
	 */
	const BObolLodProgressiveMeshPtr priorEvictedMesh =
	    reloaded.progressiveMesh;
	const size_t bytesBeforeEviction =
	    service.residentMeshBytesForDiagnostics();
	const size_t assetsBeforeEviction =
	    service.residentMeshAssetCountForDiagnostics();
	const uint64_t evictionsBefore =
	    service.residentMeshEvictionCountForDiagnostics();
	service.setResidentMeshLimit(1);
	std::vector<BObolLodResidentDemand> noDemands;
	const size_t evictionQueued =
	    service.scheduleResidentMeshCompaction(
		0x1234, 4, noDemands);
	const int evictionWait =
	    wait_for_resident_compaction(service);
	const size_t bytesAfterEviction =
	    service.residentMeshBytesForDiagnostics();
	const size_t assetsAfterEviction =
	    service.residentMeshAssetCountForDiagnostics();
	BObolLodResult restoredAfterEviction =
	    service.realizeResidentMeshLod(stagedRequest, reloadProvider);
	if (!evictionQueued || evictionWait ||
	    service.getResidentMeshLimit() != 1 ||
	    bytesAfterEviction >= bytesBeforeEviction ||
	    assetsAfterEviction >= assetsBeforeEviction ||
	    service.residentMeshEvictionCountForDiagnostics() <=
		evictionsBefore ||
	    restoredAfterEviction.providerStatus !=
		BOBOL_LOD_PROVIDER_READY ||
	    !restoredAfterEviction.progressiveMesh ||
	    restoredAfterEviction.progressiveMesh ==
		priorEvictedMesh ||
	    restoredAfterEviction.residentAdmissionRevision == 0 ||
	    restoredAfterEviction.terminal ||
	    restoredAfterEviction.geometry.activeCut !=
		restoredAfterEviction.progressiveMesh->minimumCut() ||
	    service.reservedResidentMeshGrowthBytesForDiagnostics() != 0) {
	    printf("FAIL: retained-memory pressure did not evict and "
		   "restore an undemanded PoP asset "
		   "(queued=%zu wait=%d bytes=%zu/%zu assets=%zu/%zu "
		   "evictions=%llu/%llu restored=%d fresh=%d "
		   "limited=%d active=%d min=%d revision=%llu "
		   "reserved=%zu)\n",
		   evictionQueued, evictionWait,
		   bytesAfterEviction, bytesBeforeEviction,
		   assetsAfterEviction, assetsBeforeEviction,
		   static_cast<unsigned long long>(
		       service.residentMeshEvictionCountForDiagnostics()),
		   static_cast<unsigned long long>(evictionsBefore),
		   restoredAfterEviction.providerStatus,
		   restoredAfterEviction.progressiveMesh !=
		       priorEvictedMesh ? 1 : 0,
		   restoredAfterEviction.memoryLimited ? 1 : 0,
		   restoredAfterEviction.geometry.activeCut,
		   restoredAfterEviction.progressiveMesh ?
		       restoredAfterEviction.progressiveMesh->minimumCut() :
		       -1,
		   static_cast<unsigned long long>(
		       restoredAfterEviction.residentAdmissionRevision),
		   service.
		       reservedResidentMeshGrowthBytesForDiagnostics());
	    ret = 1;
	}

	/* The first publication is deliberately the minimum coverage mesh.
	 * Its next request encounters the byte denial; a repeated direct call
	 * must report that same decision without loading the cache prefix.
	 * Raising the target advances the capacity epoch and permits the exact
	 * same view demand to finish.  The submit action suppresses that repeated
	 * direct call in production. */
	const uint64_t loadsAfterMinimum =
	    service.residentMeshCacheLoadCountForDiagnostics();
	BObolLodResult deniedSuffix =
	    service.realizeResidentMeshLod(
		stagedRequest, reloadProvider);
	const uint64_t deniedRevision =
	    deniedSuffix.residentAdmissionRevision;
	const uint64_t loadsAfterDenial =
	    service.residentMeshCacheLoadCountForDiagnostics();
	BObolLodResult repeatedDenial =
	    service.realizeResidentMeshLod(
		stagedRequest, reloadProvider);
	const uint64_t loadsAfterRepeat =
	    service.residentMeshCacheLoadCountForDiagnostics();
	const uint64_t revisionBeforeRelax =
	    service.residentMeshAdmissionRevision();
	service.setResidentMeshLimit(SIZE_MAX);
	const uint64_t revisionAfterRelax =
	    service.residentMeshAdmissionRevision();
	BObolLodResult admittedAfterRelax =
	    service.realizeResidentMeshLod(
		stagedRequest, reloadProvider);
	/* A per-task working-set cap may select a poorer presentation from an
	 * already-richer immutable resident mesh.  There is no growth reservation
	 * in that case, but the terminal capacity witness still needs the current
	 * admission revision so an unchanged submit pass can suppress retries. */
	BObolMeshLodProvider transientProvider = reloadProvider;
	transientProvider.progressiveDelivery = FALSE;
	transientProvider.useDeliveryCutLimit = TRUE;
	transientProvider.deliveryCutLimit =
	    admittedAfterRelax.progressiveMesh ?
		admittedAfterRelax.progressiveMesh->minimumCut() : 0;
	transientProvider.transientMemoryLimited = TRUE;
	BObolLodResult transientLimited =
	    service.realizeResidentMeshLod(stagedRequest, transientProvider);
	BObolLodResult repeatedTransientLimited =
	    service.realizeResidentMeshLod(stagedRequest, transientProvider);
	const uint64_t transientRevision =
	    service.residentMeshAdmissionRevision();
	if (deniedSuffix.providerStatus !=
		BOBOL_LOD_PROVIDER_READY ||
	    !deniedSuffix.memoryLimited ||
	    deniedSuffix.residentAdmissionRevision == 0 ||
	    deniedSuffix.geometry.activeCut !=
		restoredAfterEviction.geometry.activeCut ||
	    loadsAfterDenial != loadsAfterMinimum ||
	    repeatedDenial.providerStatus !=
		BOBOL_LOD_PROVIDER_READY ||
	    !repeatedDenial.memoryLimited ||
	    repeatedDenial.residentAdmissionRevision != deniedRevision ||
	    repeatedDenial.geometry.activeCut !=
		restoredAfterEviction.geometry.activeCut ||
	    loadsAfterRepeat != loadsAfterMinimum ||
	    revisionAfterRelax == revisionBeforeRelax ||
	    admittedAfterRelax.providerStatus !=
		BOBOL_LOD_PROVIDER_READY ||
	    admittedAfterRelax.memoryLimited ||
	    !admittedAfterRelax.terminal ||
	    admittedAfterRelax.geometry.activeCut !=
		stagedRequest.requestedCut ||
	    transientLimited.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	    !transientLimited.memoryLimited || !transientLimited.terminal ||
	    transientLimited.residentAdmissionRevision == 0 ||
	    transientLimited.residentAdmissionRevision != transientRevision ||
	    transientLimited.geometry.activeCut !=
		transientProvider.deliveryCutLimit ||
	    repeatedTransientLimited.providerStatus !=
		BOBOL_LOD_PROVIDER_READY ||
	    !repeatedTransientLimited.memoryLimited ||
	    repeatedTransientLimited.residentAdmissionRevision !=
		transientRevision ||
	    repeatedTransientLimited.geometry.activeCut !=
		transientLimited.geometry.activeCut) {
	    printf("FAIL: resident byte admission did not suppress/release "
		   "optional suffix growth (denied=%d/%d repeated=%d/%d "
		   "levels=%d/%d "
		   "revision=%llu/%llu/%llu admitted=%d/%d "
		   "transient=%d/%d/%llu repeated=%d/%d/%llu)\n",
		   deniedSuffix.providerStatus,
		   deniedSuffix.memoryLimited ? 1 : 0,
		   repeatedDenial.providerStatus,
		   repeatedDenial.memoryLimited ? 1 : 0,
		   repeatedDenial.geometry.activeCut,
		   restoredAfterEviction.geometry.activeCut,
		   static_cast<unsigned long long>(deniedRevision),
		   static_cast<unsigned long long>(revisionBeforeRelax),
		   static_cast<unsigned long long>(revisionAfterRelax),
		   admittedAfterRelax.memoryLimited ? 1 : 0,
		   admittedAfterRelax.geometry.activeCut,
		   transientLimited.memoryLimited ? 1 : 0,
		   transientLimited.geometry.activeCut,
		   static_cast<unsigned long long>(
		       transientLimited.residentAdmissionRevision),
		   repeatedTransientLimited.memoryLimited ? 1 : 0,
		   repeatedTransientLimited.geometry.activeCut,
		   static_cast<unsigned long long>(
		       repeatedTransientLimited.residentAdmissionRevision));
	    ret = 1;
	}
    } else {
	printf("FAIL: retained LoD compaction fixture has no richer prefix\n");
	ret = 1;
    }

    service.stop();

    /* Analytic and BREP producers hand the service owned triangle arrays,
     * rather than a database BoT or a callback into a temporary internal.
     * Exercise that cold-cache path and then discard the staging owner: the
     * resident/persistent PoP asset must remain usable on its own. */
    struct GenericStagedOwner {
	std::vector<fastf_t> points;
	std::vector<fastf_t> normals;
	std::vector<int> faces;
    };
    std::shared_ptr<GenericStagedOwner> genericOwner =
	std::make_shared<GenericStagedOwner>();
    genericOwner->points = cachedVertices;
    genericOwner->normals = cachedNormals;
    genericOwner->faces = cachedFaces;
    std::shared_ptr<BObolStagedSourceMesh> genericStaged =
	std::make_shared<BObolStagedSourceMesh>();
    genericStaged->owner = genericOwner;
    genericStaged->points = reinterpret_cast<const point_t *>(
	genericOwner->points.data());
    genericStaged->normals = reinterpret_cast<const vect_t *>(
	genericOwner->normals.data());
    genericStaged->faces = genericOwner->faces.data();
    genericStaged->pointCount = genericOwner->points.size() / 3;
    genericStaged->faceCount = genericOwner->faces.size() / 3;
    genericStaged->contentKey = 0xb5e0b5e0ULL;
    genericStaged->assetName = "lod-two-tri.bot";
    genericStaged->sourceRevision = 17;

    BObolLodService genericService;
    BObolMeshLodProvider genericProvider;
    genericProvider.service = &genericService;
    genericProvider.setDatabase(dbip);
    genericProvider.refreshMissing = TRUE;
    genericProvider.progressiveDelivery = FALSE;
    genericProvider.useForcedCut = TRUE;
    genericProvider.forcedCut = 15;
    genericProvider.stagedSource = genericStaged;
    BObolLodRequest genericRequest = make_request("/lod-two-tri.bot");
    genericRequest.objectName = "lod-two-tri.bot";
    genericRequest.sourceContentHash = 0;
    genericRequest.sourceCounts.faceCount = genericStaged->faceCount;
    genericRequest.sourceCounts.pointCount = genericStaged->pointCount;
    BObolLodResult genericResult;
    BObolLodResult genericReload;
    struct BObolMeshLodCacheStatus genericStatus =
	BOBOL_MESH_LOD_CACHE_STATUS_INIT;
    if (genericService.start(1, TRUE)) {
	genericResult = genericService.realizeResidentMeshLod(
	    genericRequest, genericProvider);
	genericProvider.stagedSource.reset();
	genericStaged.reset();
	genericOwner.reset();
	genericReload = genericService.realizeResidentMeshLod(
	    genericRequest, genericProvider);
	(void)bobol_mesh_lod_cache_status(dbip, "lod-two-tri.bot",
	    &genericStatus);
	genericService.stop();
    }
    if (genericResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	!genericResult.progressiveMesh || !genericResult.geometry.isValid() ||
	!genericStatus.has_cache_key || !genericStatus.has_cached_payload ||
	genericReload.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	genericReload.progressiveMesh != genericResult.progressiveMesh ||
	!genericReload.geometry.isValid()) {
	printf("FAIL: owned generic staged mesh did not survive cold PoP "
	       "publication (status=%d reload=%d cache=%d/%d same=%d)\n",
	    genericResult.providerStatus, genericReload.providerStatus,
	    genericStatus.has_cache_key, genericStatus.has_cached_payload,
	    genericReload.progressiveMesh == genericResult.progressiveMesh ?
		1 : 0);
	ret = 1;
    }

    BObolMeshLodProvider cachedNormalProvider;
    cachedNormalProvider.setDatabase(dbip);
    cachedNormalProvider.useForcedCut = TRUE;
    cachedNormalProvider.forcedCut = 1;
    cachedNormalProvider.refreshMissing = FALSE;
    BObolLodResult cachedNormalResult =
	bobol_mesh_lod_provider_task(task.request, &cachedNormalProvider);
    SbBool sawSeededNormal = FALSE;
    for (size_t i = 0; i < cachedNormalResult.mesh.normals.size(); i++) {
	if (cachedNormalResult.mesh.normals[i][2] > 1.5f)
	    sawSeededNormal = TRUE;
    }
    if (cachedNormalResult.resultKind != BOBOL_LOD_RESULT_MESH ||
	cachedNormalResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	cachedNormalResult.geometry.activeCut != 1 ||
	cachedNormalResult.counts.normalCount == 0 ||
	cachedNormalResult.counts.normalCount !=
	cachedNormalResult.mesh.coordIndex.size() ||
	!cachedNormalResult.hasNormals ||
	!cachedNormalResult.mesh.isValid() ||
	cachedNormalResult.mesh.normals.size() !=
	cachedNormalResult.mesh.coordIndex.size() ||
	!sawSeededNormal ||
	!bobol_lod_result_matches_request(cachedNormalResult,
					    task.request)) {
	printf("FAIL: LoD Obol mesh provider did not return cached mesh normals "
	    "(level=%d normals=%llu has=%d payload=%zu indices=%zu seeded=%d)\n",
	    cachedNormalResult.geometry.activeCut,
	    static_cast<unsigned long long>(
		cachedNormalResult.counts.normalCount),
	    cachedNormalResult.hasNormals ? 1 : 0,
	    cachedNormalResult.mesh.normals.size(),
	    cachedNormalResult.mesh.coordIndex.size(),
	    sawSeededNormal ? 1 : 0);
	ret = 1;
    }

    BObolMeshLodProvider forcedProvider;
    forcedProvider.setDatabase(dbip);
    forcedProvider.useForcedCut = TRUE;
    forcedProvider.forcedCut = (results.size() == 1 &&
				  results[0].geometry.activeCut >= 0) ?
				 results[0].geometry.activeCut : 0;
    forcedProvider.refreshMissing = FALSE;
    forcedProvider.shrinkAfterCopy = TRUE;

    BObolLodRequest forcedRequest = task.request;
    forcedRequest.qualityTier = BOBOL_LOD_QUALITY_FULL_DETAIL;
    BObolLodResult forcedResult =
	bobol_mesh_lod_provider_task(forcedRequest, &forcedProvider);
    if (forcedResult.resultKind != BOBOL_LOD_RESULT_MESH ||
	forcedResult.qualityTier != BOBOL_LOD_QUALITY_FULL_DETAIL ||
	forcedResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	forcedResult.geometry.activeCut != forcedProvider.forcedCut ||
	!forcedResult.mesh.isValid()) {
	printf("FAIL: LoD Obol mesh provider forced-level task did not return "
	    "requested level mesh result (requested=%d active=%d status=%d "
	    "valid=%d)\n", forcedProvider.forcedCut,
	    forcedResult.geometry.activeCut, forcedResult.providerStatus,
	    forcedResult.mesh.isValid() ? 1 : 0);
	ret = 1;
    }

    struct BObolMeshLodCacheStatus invalidateStatus =
	    BOBOL_MESH_LOD_CACHE_STATUS_INIT;
    if (bobol_mesh_lod_cache_invalidate(dbip, "lod-provider.bot",
					  &invalidateStatus) != BRLCAD_OK ||
	!invalidateStatus.cleared_cache_entry ||
	!invalidateStatus.cleared_cache_key ||
	invalidateStatus.has_cache_key ||
	invalidateStatus.has_cached_payload) {
	printf("FAIL: LoD Obol mesh provider database invalidation status failed\n");
	ret = 1;
    }

    BObolMeshLodProvider staleProvider;
    staleProvider.setDatabase(dbip);
    staleProvider.refreshMissing = FALSE;

    BObolLodResult staleResult =
	bobol_mesh_lod_provider_task(task.request, &staleProvider);
    if (staleResult.providerStatus != BOBOL_LOD_PROVIDER_CACHE_MISS ||
	!staleResult.stale ||
	strstr(staleResult.diagnostic.getString(),
	       "Obol mesh LoD provider has no cache payload") == NULL ||
	staleResult.mesh.isValid()) {
	printf("FAIL: LoD Obol mesh provider did not report cache miss after database invalidation\n");
	ret = 1;
    }

    BObolMeshLodProvider refreshProvider;
    refreshProvider.setDatabase(dbip);
    refreshProvider.refreshMissing = TRUE;

    BObolLodResult refreshResult =
	bobol_mesh_lod_provider_task(task.request, &refreshProvider);
    if (refreshResult.resultKind != BOBOL_LOD_RESULT_MESH ||
	refreshResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	refreshResult.geometry.kind != BOBOL_LOD_GEOMETRY_MESH_LOD_CACHE ||
	!refreshResult.geometry.isValid() ||
	refreshResult.counts.faceCount != 4 ||
	refreshResult.counts.pointCount != 4 ||
	!refreshResult.mesh.isValid() ||
	!bobol_lod_result_matches_request(refreshResult, task.request)) {
	printf("FAIL: LoD Obol mesh provider did not refresh cache after database invalidation\n");
	ret = 1;
    }

    bobol_mesh_lod_cache_clear_database(dbip);
    db_close(dbip);
    bu_file_delete(dbpath);
    bu_dirclear(cache_dir);
    return ret;
}

static int
test_detached_database_snapshot_version(int databaseVersion)
{
    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;
    struct db_i *snapshotDb = NULL;
    struct db_i *postCloseSnapshotDb = NULL;
    SoBRLDatabaseSource *source = NULL;
    SoBRLDatabaseSource *detached = NULL;
    SoBRLDatabaseSource *postCloseDetached = NULL;
    SbString snapshotPath;
    SbString postCloseSnapshotPath;
    struct directory *dp = NULL;
    struct rt_db_internal intern;
    struct rt_bot_internal *bot = NULL;
    SbBox3f bounds;
    int workerRealized = 0;
    int postCloseWorkerRealized = 0;
    int mutationResult = 0;
    int ret = 1;

    RT_DB_INTERNAL_INIT(&intern);


    if (make_provider_test_db_version(dbpath, sizeof(dbpath), &dbip,
	    databaseVersion))
	return 1;

    source = new SoBRLDatabaseSource;
    source->ref();
    const char *sourcePath = databaseVersion == 4 ?
	"/lod-provider.bot" : "/lod-snapshot.c";
    source->configureDatabaseSourceInstance("snapshot-root",
	sourcePath, dbip, SoBRLDatabaseSource::WIREFRAME, 1);
    detached = source->createDetachedRealizationSource(&snapshotDb,
	&snapshotPath);
    if (!detached || !snapshotDb || snapshotDb == dbip) {
	printf("FAIL: detached database source did not create an independent snapshot\n");
	goto cleanup;
    }
    if (databaseVersion == 5 &&
	!db_lookup(snapshotDb, DB5_GLOBAL_OBJECT_NAME, LOOKUP_QUIET)) {
	printf("FAIL: detached database snapshot did not retain _GLOBAL metadata\n");
	goto cleanup;
    }
    if ((dbip->dbi_filename &&
	 (!snapshotDb->dbi_filename ||
	  bu_strcmp(snapshotDb->dbi_filename, dbip->dbi_filename) != 0)) ||
	(dbip->dbi_filepath &&
	 (!snapshotDb->dbi_filepath ||
	  bu_strcmp(snapshotDb->dbi_filepath[0], dbip->dbi_filepath[0]) != 0 ||
	  bu_strcmp(snapshotDb->dbi_filepath[1], dbip->dbi_filepath[1]) != 0))) {
	printf("FAIL: detached database snapshot did not retain file lookup context\n");
	goto cleanup;
    }

    dp = db_lookup(dbip, "lod-provider.bot", LOOKUP_QUIET);
    if (!dp || rt_db_get_internal(&intern, dp, dbip, NULL) < 0 ||
	intern.idb_type != ID_BOT || !intern.idb_ptr) {
	printf("FAIL: detached database snapshot mutation setup\n");
	if (intern.idb_ptr)
	    rt_db_free_internal(&intern);
	goto cleanup;
    }
    bot = static_cast<struct rt_bot_internal *>(intern.idb_ptr);
    for (size_t i = 0; i < bot->num_vertices; i++)
	bot->vertices[i * 3] += 100.0;
    {
	std::atomic<int> workerResult(0);
	std::thread worker([&]() {
	    workerResult.store(detached->realizeDatabaseWireframe() ? 1 : -1,
		std::memory_order_release);
	});
	mutationResult = rt_db_put_internal(dp, dbip, &intern);
	worker.join();
	workerRealized = workerResult.load(std::memory_order_acquire);
    }
    if (mutationResult < 0) {
	printf("FAIL: detached database snapshot live mutation\n");
	goto cleanup;
    }

    if (workerRealized != 1) {
	printf("FAIL: detached database snapshot worker realization\n");
	goto cleanup;
    }
    if (!detached->getSourceBounds(bounds) || bounds.isEmpty() ||
	bounds.getMin()[0] < -0.01f || bounds.getMax()[0] > 1.01f) {
	printf("FAIL: detached database snapshot observed post-snapshot mutation\n");
	goto cleanup;
    }

    postCloseDetached = source->createDetachedRealizationSource(
	&postCloseSnapshotDb, &postCloseSnapshotPath);
    if (!postCloseDetached || !postCloseSnapshotDb ||
	postCloseSnapshotDb == dbip) {
	printf("FAIL: detached database source did not snapshot for close test\n");
	goto cleanup;
    }

    source->unref();
    source = NULL;
    db_close(dbip);
    dbip = NULL;
    bu_file_delete(dbpath);
    dbpath[0] = '\0';

    {
	std::atomic<int> workerResult(0);
	std::thread worker([&]() {
	    workerResult.store(
		postCloseDetached->realizeDatabaseWireframe() ? 1 : -1,
		std::memory_order_release);
	});
	worker.join();
	postCloseWorkerRealized = workerResult.load(std::memory_order_acquire);
    }
    bounds.makeEmpty();
    if (postCloseWorkerRealized != 1 ||
	!postCloseDetached->getSourceBounds(bounds) || bounds.isEmpty() ||
	bounds.getMin()[0] < 99.99f || bounds.getMin()[0] > 100.01f ||
	bounds.getMax()[0] < 100.99f || bounds.getMax()[0] > 101.01f) {
	printf("FAIL: detached database snapshot did not survive live database close\n");
	goto cleanup;
    }

    ret = 0;
cleanup:
    if (postCloseDetached)
	postCloseDetached->unref();
    if (postCloseSnapshotDb)
	db_close(postCloseSnapshotDb);
    if (postCloseSnapshotPath.getLength() > 0)
	bu_file_delete(postCloseSnapshotPath.getString());
    if (detached)
	detached->unref();
    if (snapshotDb)
	db_close(snapshotDb);
    if (snapshotPath.getLength() > 0)
	bu_file_delete(snapshotPath.getString());
    if (source)
	source->unref();
    if (dbip)
	db_close(dbip);
    if (dbpath[0])
	bu_file_delete(dbpath);
    return ret;
}

static int
test_detached_database_snapshot(void)
{
    return test_detached_database_snapshot_version(5) ||
	test_detached_database_snapshot_version(4) ? 1 : 0;
}

static int
test_rt_obb_proxy_provider_request_bounds(void)
{
    {
	BObolLodRequest request = make_request("/proxy.bot");
	BObolRtProxyProvider provider;

	provider.proxyKind = BOBOL_LOD_PROXY_OBB;
	provider.useRequestBounds = TRUE;
	BObolLodResult result = bobol_rt_proxy_provider_task(request,
							  &provider);
	if (result.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	    result.resultKind != BOBOL_LOD_RESULT_PROXY ||
	    result.proxy.kind != BOBOL_LOD_PROXY_OBB ||
	    !result.proxy.isValid() ||
	    result.proxy.bounds != request.bounds ||
	    result.proxy.center != SbVec3f(0.0f, 0.0f, 0.0f) ||
	    result.proxy.halfExtents != SbVec3f(1.0f, 1.0f, 1.0f)) {
	    printf("FAIL: LoD OBB proxy provider did not preserve the OBB stage from request bounds\n");
	    return 1;
	}
    }

    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;
    if (make_provider_test_db(dbpath, sizeof(dbpath), &dbip))
	return 1;

    BObolLodRequest request = make_request("/lod-provider.bot");
    request.objectName = "lod-provider.bot";
    BObolRtProxyProvider provider;
    provider.proxyKind = BOBOL_LOD_PROXY_OBB;
    provider.useRequestBounds = FALSE;
    if (!provider.setDatabase(dbip)) {
	printf("FAIL: LoD OBB proxy provider database lease\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    const BObolLodResult result = bobol_rt_proxy_provider_task(request,
	&provider);
    const float xy = result.proxy.axisX.dot(result.proxy.axisY);
    const float xz = result.proxy.axisX.dot(result.proxy.axisZ);
    const float yz = result.proxy.axisY.dot(result.proxy.axisZ);
    constexpr float orthogonalityTolerance = 1.0e-5f;
    if (result.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	result.resultKind != BOBOL_LOD_RESULT_PROXY ||
	result.proxy.kind != BOBOL_LOD_PROXY_OBB ||
	!result.proxy.isValid() ||
	std::fabs(xy) > orthogonalityTolerance ||
	std::fabs(xz) > orthogonalityTolerance ||
	std::fabs(yz) > orthogonalityTolerance) {
	printf("FAIL: cached LoD OBB proxy lost canonical orthogonal axes\n");
	provider.clear();
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    provider.clear();
    db_close(dbip);
    bu_file_delete(dbpath);
    return 0;
}

struct WorkingSetTaskState {
    std::atomic<int> active{0};
    std::atomic<int> maximum{0};
};

static BObolLodResult
working_set_task(const BObolLodRequest &request, void *userData)
{
    WorkingSetTaskState *state =
	static_cast<WorkingSetTaskState *>(userData);
    const int active = state->active.fetch_add(1) + 1;
    int observed = state->maximum.load();
    while (active > observed &&
	!state->maximum.compare_exchange_weak(observed, active)) {
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(40));
    state->active.fetch_sub(1);

    BObolLodResult result;
    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = BOBOL_LOD_RESULT_AABB;
    result.qualityTier = request.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.terminal = TRUE;
    result.bounds = request.bounds;
    return result;
}

static int
test_working_set_admission(void)
{
    BObolLodService service;
    WorkingSetTaskState state;
    service.setWorkingSetLimit(0);
    if (service.getWorkingSetLimit() == 0) {
	printf("FAIL: automatic LoD working-set admission limit\n");
	return 1;
    }
    service.setWorkingSetLimit(100);
    if (service.getWorkingSetLimit() != 100 || !service.start(4, FALSE)) {
	printf("FAIL: LoD working-set admission setup\n");
	return 1;
    }

    for (int i = 0; i < 3; ++i) {
	BObolLodTask task;
	SbString name;
	name.sprintf("/working-set-%d.bot", i);
	task.request = make_request(name.getString());
	task.realize = working_set_task;
	task.realizeData = &state;
	task.estimatedWorkingSetBytes = 60;
	if (!service.submit(task)) {
	    printf("FAIL: LoD working-set task submission\n");
	    service.stop();
	    return 1;
	}
    }
    const int waitResult = wait_for_settled(service, 3);
    const int maximum = state.maximum.load();
    const size_t activeBytes =
	service.activeWorkingSetBytesForDiagnostics();
    const size_t peakBytes =
	service.peakWorkingSetBytesForDiagnostics();
    const size_t peakTasks =
	service.peakExecutingTaskCountForDiagnostics();
    if (waitResult || maximum != 1 || activeBytes != 0 ||
	peakBytes != 60 || peakTasks != 1) {
	printf("FAIL: byte-weighted LoD admission exceeded its working-set "
	       "limit (maximum=%d active_bytes=%zu peak_bytes=%zu "
	       "peak_tasks=%zu)\n",
	       maximum, activeBytes, peakBytes, peakTasks);
	service.stop();
	return 1;
    }

    WorkingSetTaskState rejected;
    BObolLodTask oversized;
    oversized.request = make_request("/working-set-oversized.bot");
    oversized.realize = working_set_task;
    oversized.realizeData = &rejected;
    oversized.estimatedWorkingSetBytes = 101;
    if (!service.submit(oversized) || wait_for_settled(service, 1) ||
	rejected.maximum.load() != 0 ||
	service.activeWorkingSetBytesForDiagnostics() != 0 ||
	service.peakWorkingSetBytesForDiagnostics() != 60) {
	printf("FAIL: service working-set governor entered oversized task\n");
	service.stop();
	return 1;
    }
    service.stop();
    return 0;
}

static int
test_process_working_set_admission(void)
{
    const size_t limit = bobol_lod_working_set_global_limit();
    if (limit < 4 || limit == SIZE_MAX ||
	bobol_lod_working_set_global_active_bytes() != 0 ||
	bobol_lod_working_set_global_active_tasks() != 0) {
	printf("FAIL: process LoD working-set governor setup\n");
	return 1;
    }

    const size_t reservation = limit - limit / 4;
    WorkingSetTaskState state;
    std::vector<std::thread> workers;
    for (int i = 0; i < 3; i++) {
	workers.push_back(std::thread([&]() {
	    if (!bobol_lod_working_set_acquire(reservation))
		return;
	    const int active = state.active.fetch_add(1) + 1;
	    int observed = state.maximum.load();
	    while (active > observed &&
		!state.maximum.compare_exchange_weak(observed, active)) {
	    }
	    std::this_thread::sleep_for(std::chrono::milliseconds(20));
	    state.active.fetch_sub(1);
	    bobol_lod_working_set_release(reservation);
	}));
    }
    for (std::thread &worker : workers)
	worker.join();

    if (state.maximum.load() != 1 ||
	bobol_lod_working_set_global_active_bytes() != 0 ||
	bobol_lod_working_set_global_active_tasks() != 0 ||
	bobol_lod_working_set_global_peak_bytes() < reservation ||
	bobol_lod_working_set_global_peak_tasks() < 1) {
	printf("FAIL: process LoD working-set governor admitted overlapping "
	       "large preparations\n");
	return 1;
    }
    return 0;
}

static int
test_process_working_set_rejection(void)
{
    const size_t limit = bobol_lod_working_set_global_limit();
    if (limit == SIZE_MAX)
	return 0;
    if (limit == SIZE_MAX - 1 ||
	bobol_lod_working_set_acquire(limit + 1) ||
	bobol_lod_working_set_global_active_bytes() != 0 ||
	bobol_lod_working_set_global_active_tasks() != 0) {
	printf("FAIL: process LoD governor admitted an oversized task\n");
	return 1;
    }
    return 0;
}

struct GenerationIsolationState {
    std::mutex mutex;
    std::condition_variable cv;
    bool entered = false;
    bool release = false;
    std::atomic<unsigned int> notifications{0};
};

static BObolLodResult
generation_blocking_task(const BObolLodRequest &request, void *userData)
{
    GenerationIsolationState *state =
	static_cast<GenerationIsolationState *>(userData);
    {
	std::unique_lock<std::mutex> lock(state->mutex);
	state->entered = true;
	state->cv.notify_all();
	state->cv.wait(lock, [state] { return state->release; });
    }

    BObolLodResult result;
    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = BOBOL_LOD_RESULT_AABB;
    result.qualityTier = request.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.terminal = TRUE;
    result.counts.faceCount = 101;
    result.bounds = request.bounds;
    return result;
}

static BObolLodResult
generation_ready_task(const BObolLodRequest &request, void *)
{
    BObolLodResult result;
    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = BOBOL_LOD_RESULT_AABB;
    result.qualityTier = request.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.terminal = TRUE;
    result.counts.faceCount = 202;
    result.bounds = request.bounds;
    return result;
}

static void
generation_result_ready(BObolLodService *, void *userData)
{
    GenerationIsolationState *state =
	static_cast<GenerationIsolationState *>(userData);
    state->notifications.fetch_add(1);
}

static int
test_generation_scoped_consumers(void)
{
    BObolLodService service;
    GenerationIsolationState state;
    if (!service.start(1, FALSE)) {
	printf("FAIL: generation-isolation service setup\n");
	return 1;
    }
    const BObolLodSubscriberId subscriber =
	service.subscribeResultReady(generation_result_ready, &state);
    const uint64_t firstGeneration = service.beginGeneration();
    const uint64_t secondGeneration = service.beginGeneration();

    BObolLodTask first;
    first.generation = firstGeneration;
    first.request = make_request("/generation-first.bot");
    first.realize = generation_blocking_task;
    first.realizeData = &state;
    if (!service.submit(first)) {
	printf("FAIL: first generation submission\n");
	return 1;
    }
    {
	std::unique_lock<std::mutex> lock(state.mutex);
	if (!state.cv.wait_for(lock, std::chrono::seconds(2),
		[&state] { return state.entered; })) {
	    printf("FAIL: first generation did not begin execution\n");
	    return 1;
	}
    }

    BObolLodTask second;
    second.generation = secondGeneration;
    second.request = make_request("/generation-second.bot");
    second.realize = generation_ready_task;
    if (!service.submit(second) ||
	service.activeTaskCountForGeneration(firstGeneration) != 1 ||
	service.executingTaskCountForGeneration(firstGeneration) != 1 ||
	service.pendingTaskCountForGeneration(firstGeneration) != 0 ||
	service.activeTaskCountForGeneration(secondGeneration) != 1 ||
	service.executingTaskCountForGeneration(secondGeneration) != 0 ||
	service.pendingTaskCountForGeneration(secondGeneration) != 1) {
	printf("FAIL: per-generation active/pending/executing counters\n");
	{
	    std::lock_guard<std::mutex> lock(state.mutex);
	    state.release = true;
	}
	state.cv.notify_all();
	return 1;
    }

    {
	std::lock_guard<std::mutex> lock(state.mutex);
	state.release = true;
    }
    state.cv.notify_all();
    if (wait_for_settled(service, 2))
	return 1;
    for (int i = 0; i < 200 && state.notifications.load() < 2; ++i)
	std::this_thread::sleep_for(std::chrono::milliseconds(5));

    if (service.activeTaskCountForGeneration(firstGeneration) != 0 ||
	service.activeTaskCountForGeneration(secondGeneration) != 0 ||
	service.queuedResultCountForGeneration(firstGeneration) != 1 ||
	service.queuedResultCountForGeneration(secondGeneration) != 1 ||
	state.notifications.load() != 2) {
	printf("FAIL: per-generation completion/result notification state "
	       "(notifications=%u)\n", state.notifications.load());
	return 1;
    }

    std::vector<BObolLodResult> secondResults;
    if (service.drainGenerationResults(
	    secondResults, secondGeneration) != 1 ||
	secondResults.size() != 1 ||
	secondResults.front().generation != secondGeneration ||
	secondResults.front().counts.faceCount != 202 ||
	service.queuedResultCountForGeneration(secondGeneration) != 0 ||
	service.queuedResultCountForGeneration(firstGeneration) != 1 ||
	service.queuedResultCountForDiagnostics() != 1) {
	printf("FAIL: draining one generation consumed another generation\n");
	return 1;
    }

    std::vector<BObolLodResult> firstResults;
    if (service.drainGenerationResults(firstResults, firstGeneration) != 1 ||
	firstResults.size() != 1 ||
	firstResults.front().generation != firstGeneration ||
	firstResults.front().counts.faceCount != 101 ||
	service.queuedResultCountForDiagnostics() != 0) {
	printf("FAIL: final generation-scoped drain\n");
	return 1;
    }

    service.unsubscribeResultReady(subscriber);
    service.stop();
    return 0;
}

struct IntermediateResultState {
    BObolLodService *service = NULL;
    uint64_t generation = 0;
    std::mutex mutex;
    std::condition_variable cv;
    bool entered = false;
    bool published = false;
    bool release = false;
};

static BObolLodResult
intermediate_result_task(const BObolLodRequest &request, void *userData)
{
    IntermediateResultState *state =
	static_cast<IntermediateResultState *>(userData);
    BObolLodResult preview;
    preview.request = request;
    preview.cacheKey = bobol_lod_cache_key(request);
    preview.resultKind = BOBOL_LOD_RESULT_AABB;
    preview.qualityTier = request.qualityTier;
    preview.providerStatus = BOBOL_LOD_PROVIDER_READY;
    preview.terminal = FALSE;
    preview.counts.faceCount = 11;
    preview.bounds = request.bounds;

    const bool published = state && state->service &&
	state->service->tryPublishIntermediateResult(
	    state->generation, std::move(preview));
    if (state) {
	std::unique_lock<std::mutex> lock(state->mutex);
	state->published = published;
	state->entered = true;
	state->cv.notify_all();
	state->cv.wait(lock, [state] { return state->release; });
    }

    BObolLodResult result;
    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = BOBOL_LOD_RESULT_AABB;
    result.qualityTier = request.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.terminal = TRUE;
    result.counts.faceCount = 22;
    result.bounds = request.bounds;
    return result;
}

static int
test_intermediate_result_lifecycle(void)
{
    BObolLodService service;
    IntermediateResultState state;
    if (!service.start(1, FALSE)) {
	printf("FAIL: intermediate-result service setup\n");
	return 1;
    }
    const uint64_t generation = service.beginGeneration();
    state.service = &service;
    state.generation = generation;

    BObolLodTask task;
    task.generation = generation;
    task.request = make_request("/intermediate.bot");
    task.realize = intermediate_result_task;
    task.realizeData = &state;
    if (!service.submit(task)) {
	printf("FAIL: intermediate-result submission\n");
	return 1;
    }
    {
	std::unique_lock<std::mutex> lock(state.mutex);
	if (!state.cv.wait_for(lock, std::chrono::seconds(2),
		[&state] { return state.entered; }) || !state.published) {
	    printf("FAIL: intermediate result was not published\n");
	    state.release = true;
	    lock.unlock();
	    state.cv.notify_all();
	    return 1;
	}
    }

    std::vector<BObolLodResult> preview;
    if (service.drainGenerationResults(preview, generation) != 1 ||
	preview.size() != 1 || preview.front().terminal ||
	preview.front().counts.faceCount != 11 ||
	service.activeTaskCountForGeneration(generation) != 1) {
	printf("FAIL: intermediate result did not preserve active task state\n");
	{
	    std::lock_guard<std::mutex> lock(state.mutex);
	    state.release = true;
	}
	state.cv.notify_all();
	return 1;
    }

    {
	std::lock_guard<std::mutex> lock(state.mutex);
	state.release = true;
    }
    state.cv.notify_all();
    if (wait_for_settled(service, 1))
	return 1;

    std::vector<BObolLodResult> finalResult;
    if (service.drainGenerationResults(finalResult, generation) != 1 ||
	finalResult.size() != 1 || !finalResult.front().terminal ||
	finalResult.front().counts.faceCount != 22 ||
	service.activeTaskCountForGeneration(generation) != 0 ||
	service.queuedResultCountForDiagnostics() != 0) {
	printf("FAIL: final result did not replace intermediate lifecycle\n");
	return 1;
    }

    service.stop();
    return 0;
}

static int
test_resident_memory_percentage_limit(void)
{
    BObolLodService service;
    const double maximumPercent =
	service.getMaximumResidentMeshAvailableMemoryPercent();
    if (std::fabs(maximumPercent - 50.0) >
	std::numeric_limits<double>::epsilon() ||
	std::fabs(service.getResidentMeshAvailableMemoryPercent()) >
	    std::numeric_limits<double>::epsilon() ||
	service.getResidentMeshAvailableMemoryBasisBytes() != 0) {
	printf("FAIL: resident memory percentage policy did not start in auto mode\n");
	return 1;
    }

    if (service.setResidentMeshAvailableMemoryPercent(0.0) ||
	service.setResidentMeshAvailableMemoryPercent(-1.0) ||
	service.setResidentMeshAvailableMemoryPercent(maximumPercent + 1.0) ||
	service.setResidentMeshAvailableMemoryPercent(
	    std::numeric_limits<double>::quiet_NaN()) ||
	service.setResidentMeshAvailableMemoryPercent(
	    std::numeric_limits<double>::infinity())) {
	printf("FAIL: resident memory percentage policy accepted an invalid value\n");
	return 1;
    }

    const double testPercent = 25.0;
    if (!service.setResidentMeshAvailableMemoryPercent(testPercent)) {
	printf("FAIL: resident memory percentage policy rejected a valid value\n");
	return 1;
    }
    const size_t basisBytes =
	service.getResidentMeshAvailableMemoryBasisBytes();
    const size_t expectedBytes = basisBytes / 4;
    if (std::fabs(service.getResidentMeshAvailableMemoryPercent() -
	    testPercent) > std::numeric_limits<double>::epsilon() ||
	basisBytes == 0 || expectedBytes == 0 ||
	service.getResidentMeshLimit() != expectedBytes) {
	printf("FAIL: resident memory percentage policy did not retain its snapshot\n");
	return 1;
    }

    const size_t directLimit = 12345;
    service.setResidentMeshLimit(directLimit);
    if (service.getResidentMeshLimit() != directLimit ||
	std::fabs(service.getResidentMeshAvailableMemoryPercent()) >
	    std::numeric_limits<double>::epsilon() ||
	service.getResidentMeshAvailableMemoryBasisBytes() != 0) {
	printf("FAIL: direct resident memory limit did not clear percentage state\n");
	return 1;
    }
    service.setResidentMeshLimit(0);
    if (service.getResidentMeshLimit() == 0 ||
	std::fabs(service.getResidentMeshAvailableMemoryPercent()) >
	    std::numeric_limits<double>::epsilon() ||
	service.getResidentMeshAvailableMemoryBasisBytes() != 0) {
	printf("FAIL: automatic resident memory policy was not restored\n");
	return 1;
    }
    return 0;
}

static int
test_managed_service_is_shared(void)
{
    BObolViewController first;
    BObolViewController second;
    BObolLodService *firstService = first.ensureManagedLodService(1);
    BObolLodService *secondService = second.ensureManagedLodService(2);
    if (!firstService || secondService != firstService ||
	first.getManagedLodWorkerCount() < 2 ||
	second.getManagedLodWorkerCount() < 2) {
	printf("FAIL: view controllers did not share and grow the managed service\n");
	return 1;
    }

    first.stopManagedLodService();
    if (first.getLodService() != NULL ||
	second.getLodService() != secondService || !secondService->isRunning()) {
	printf("FAIL: releasing one view stopped another view's managed service\n");
	return 1;
    }
    second.stopManagedLodService();
    return 0;
}

static int
test_presentation_layer_validation(void)
{
    BObolLodResult result;
    result.request = make_request("presentation-layer");
    result.resultKind = BOBOL_LOD_RESULT_MESH;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    BObolLodPresentationLayer coverage;
    coverage.layerKey = "coverage";
    const Obol::CadGeometryAdmission coverageAdmission =
	Obol::cadAdmitPartGeometry(Obol::PartGeometryBuilder());
    if (!coverageAdmission) {
	printf("FAIL: presentation layer geometry admission\n");
	return 1;
    }
    coverage.geometry = coverageAdmission.geometry.shared();
    coverage.geometryRevision = 1;
    coverage.coverage = TRUE;
    result.presentationLayers.push_back(coverage);
    result.canonicalizePayload();
    if (!result.payloadIsConsistent()) {
	printf("FAIL: valid presentation layer rejected\n");
	return 1;
    }

    result.presentationLayers.push_back(coverage);
    if (result.payloadIsConsistent()) {
	printf("FAIL: duplicate presentation layer accepted\n");
	return 1;
    }
    result.presentationLayers.pop_back();
    result.presentationLayers[0].geometry.reset();
    if (result.payloadIsConsistent()) {
	printf("FAIL: invalid presentation layer accepted\n");
	return 1;
    }
    return 0;
}

int
main(int argc, char **argv)
{
    if (argc != 1) {
	printf("Usage: %s\n", argv[0]);
	return 1;
    }
    bu_setprogname(argv[0]);

    char processCacheDir[MAXPATHLEN] = {0};
    bu_dir(processCacheDir, MAXPATHLEN, BU_DIR_CURR,
	"bobol_lod_service_process_cache", NULL);
    bu_dirclear(processCacheDir);
    bu_mkdir(processCacheDir);
    bu_setenv("BU_DIR_CACHE", processCacheDir, 1);
    char resolvedCacheDir[MAXPATHLEN] = {0};
    bu_dir(resolvedCacheDir, MAXPATHLEN, BU_DIR_CACHE, NULL);
    if (!BU_STR_EQUAL(processCacheDir, resolvedCacheDir)) {
	printf("FAIL: LoD service test cache isolation (%s != %s)\n",
	    processCacheDir, resolvedCacheDir);
	return 1;
    }
    bobol_init(NULL);

    const auto runIsolated = [&processCacheDir](int (*test)(void)) {
	bu_dirclear(processCacheDir);
	bu_mkdir(processCacheDir);
	bu_setenv("BU_DIR_CACHE", processCacheDir, 1);
	return test();
    };

    if (runIsolated(test_dependency_order_and_cache_write))
	return 1;
    if (runIsolated(test_filtered_result_drain))
	return 1;
    if (runIsolated(test_occurrence_result_coalescing))
	return 1;
    if (runIsolated(test_generation_cancellation))
	return 1;
    if (runIsolated(test_destruction_waits_for_in_flight_task))
	return 1;
    if (runIsolated(test_stop_discards_undrained_state))
	return 1;
    if (runIsolated(test_superseded_result_rejection))
	return 1;
    if (runIsolated(test_staged_payload_delivery))
	return 1;
    if (runIsolated(test_result_ready_subscription))
	return 1;
    if (runIsolated(test_result_ready_self_unsubscribe))
	return 1;
    if (runIsolated(test_result_ready_cross_unsubscribe))
	return 1;
    if (runIsolated(test_task_realize_data_cleanup))
	return 1;
    if (runIsolated(test_queue_limits_and_pending_cancellation))
	return 1;
    if (runIsolated(test_large_pending_cancellation_and_generation_history))
	return 1;
    if (runIsolated(test_debug_delay_cancellation))
	return 1;
    if (runIsolated(test_cancelled_cache_write_not_persisted))
	return 1;
    if (runIsolated(test_active_request_duplicate_suppression))
	return 1;
    if (runIsolated(test_asset_producer_coalescing))
	return 1;
    if (runIsolated(test_shared_producer_two_generation_delivery))
	return 1;
    if (runIsolated(test_shared_producer_first_consumer_cancellation))
	return 1;
    if (runIsolated(test_shared_producer_late_join_during_build))
	return 1;
    if (runIsolated(test_shared_producer_late_join_while_result_queued))
	return 1;
    if (runIsolated(
	    test_shared_producer_final_consumer_cancellation_during_build))
	return 1;
    if (runIsolated(
	    test_shared_producer_final_consumer_cancellation_with_result))
	return 1;
    if (runIsolated(test_batch_submission))
	return 1;
    if (runIsolated(test_rt_source_full_detail_provider_task))
	return 1;
    if (runIsolated(test_database_lease_outlives_submitter))
	return 1;
    if (runIsolated(test_rt_mesh_provider_task))
	return 1;
    if (runIsolated(test_detached_database_snapshot))
	return 1;
    if (runIsolated(test_rt_obb_proxy_provider_request_bounds))
	return 1;
    if (runIsolated(test_working_set_admission))
	return 1;
    if (runIsolated(test_process_working_set_admission))
	return 1;
    if (runIsolated(test_process_working_set_rejection))
	return 1;
    if (runIsolated(test_generation_scoped_consumers))
	return 1;
    if (runIsolated(test_intermediate_result_lifecycle))
	return 1;
    if (runIsolated(test_resident_memory_percentage_limit))
	return 1;
    if (runIsolated(test_managed_service_is_shared))
	return 1;
    if (runIsolated(test_presentation_layer_validation))
	return 1;

    bu_dirclear(processCacheDir);
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
