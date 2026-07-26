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

#include <Inventor/SbBox.h>
#include <Inventor/SbVec3f.h>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <mutex>
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

    if (service.submit(oldTask) == 0 || service.submit(newTask) == 0 ||
	service.submit(meshTask) == 0 || wait_for_settled(service, 2)) {
	printf("FAIL: LoD service did not execute coalescing tasks\n");
	service.stop();
	return 1;
    }

    std::vector<BObolLodResult> results;
    service.drainResults(results);
    const BObolLodResult *aabb = NULL;
    const BObolLodResult *mesh = NULL;
    for (const BObolLodResult &result : results) {
	if (result.resultKind == BOBOL_LOD_RESULT_AABB)
	    aabb = &result;
	if (result.resultKind == BOBOL_LOD_RESULT_MESH)
	    mesh = &result;
    }
    if (results.size() != 2 || !aabb || !mesh ||
	aabb->request.viewRevision != newTask.request.viewRevision ||
	aabb->counts.faceCount != 2 || mesh->counts.faceCount != 3 ||
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
stale_task(const BObolLodRequest &request, void *UNUSED(userData))
{
    BObolLodRequest stale = request;
    stale.viewRevision++;

    BObolLodResult result;
    result.request = stale;
    result.cacheKey = bobol_lod_cache_key(stale);
    result.resultKind = BOBOL_LOD_RESULT_MESH;
    result.qualityTier = stale.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.counts.faceCount = 11;
    return result;
}

static int
test_stale_result_rejection(void)
{
    BObolLodService service;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for stale-result test\n");
	return 1;
    }

    BObolLodTask task;
    task.generation = service.beginGeneration();
    task.request = make_request("/stale.bot");
    task.realize = stale_task;

    if (service.submit(task) == 0) {
	printf("FAIL: LoD service did not accept stale-result task\n");
	return 1;
    }

    if (wait_for_settled(service, 1))
	return 1;

    std::vector<BObolLodResult> results;
    service.drainResults(results);
    if (results.size() != 1 ||
	results[0].providerStatus != BOBOL_LOD_PROVIDER_STALE ||
	!results[0].stale ||
	!results[0].terminal ||
	results[0].diagnostic.getLength() == 0) {
	printf("FAIL: LoD service did not reject stale result\n");
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

    DebugDelayTaskData(void) : calls(0)
    {
    }
};

static BObolLodResult
debug_delay_task(const BObolLodRequest &request, void *userData)
{
    DebugDelayTaskData *data = static_cast<DebugDelayTaskData *>(userData);

    {
	std::lock_guard<std::mutex> lock(data->mutex);
	data->calls++;
    }

    BObolLodResult result;
    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = BOBOL_LOD_RESULT_AABB;
    result.qualityTier = request.qualityTier;
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

    if (wait_for_settled(service, 1)) {
	service.stop();
	return 1;
    }
    if (service.hasActiveRequest(task.request)) {
	printf("FAIL: LoD service retained active request after completion\n");
	service.stop();
	return 1;
    }

    std::vector<BObolLodResult> results;
    if (service.drainResults(results) != 1 || results.size() != 1) {
	printf("FAIL: LoD service duplicate suppression changed result count\n");
	service.stop();
	return 1;
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
    provider.dbip = dbip;
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
    scopedLimitProvider.dbip = dbip;
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
    wrongSpaceLimitProvider.dbip = dbip;
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
    compactRayLimitProvider.dbip = dbip;
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
	templatedLimitProvider.dbip = dbip;
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
    limitProvider.dbip = dbip;
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

    const int cachedGrid = 8;
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
	printf("FAIL: LoD Obol mesh provider did not store cached mesh normals\n");
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    BObolMeshLodProvider provider;
    provider.dbip = dbip;
    provider.useView = TRUE;
    provider.refreshMissing = TRUE;
    provider.view.size = 100.0;
    provider.view.width = 640;
    provider.view.height = 480;

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
	results[0].counts.faceCount <= 4 ||
	results[0].counts.faceCount > cachedFaces.size() / 3 ||
	results[0].counts.pointCount <= 4 ||
	results[0].counts.pointCount > cachedVertices.size() / 3 ||
	!results[0].mesh.isValid() ||
	results[0].mesh.points.size() != results[0].counts.pointCount ||
	results[0].mesh.coordIndex.size() !=
	    results[0].counts.faceCount * 3 ||
	!bobol_lod_result_matches_request(results[0], task.request)) {
	printf("FAIL: LoD Obol mesh provider task did not return cached mesh result (results=%zu kind=%d quality=%d status=%d geometry=%d token=%llu valid=%d faces=%llu points=%llu mesh_valid=%d mesh_points=%zu indices=%zu match=%d level=%d requested=%d)\n",
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
	    results.empty() ? -1 : results[0].geometry.activeLevel,
	    task.request.requestedLevel);
	ret = 1;
    }

    /* The retained provider must publish a bounded first prefix and then make
     * monotonic progress toward the unchanged view target.  Tiny budgets make
     * the behavior observable with this small fixture. */
    BObolMeshLodProvider stagedProvider;
    stagedProvider.service = &service;
    stagedProvider.dbip = dbip;
    stagedProvider.refreshMissing = FALSE;
    stagedProvider.initialRefinementFaceBudget = 2;
    stagedProvider.initialRefinementPointBudget = 4;
    stagedProvider.refinementGrowthFactor = 2.0;
    BObolLodRequest stagedRequest = task.request;
    stagedRequest.requestedLevel = 4;
    BObolLodResult firstStage =
	service.realizeResidentMeshLod(stagedRequest, stagedProvider);
    BObolLodResult secondStage =
	service.realizeResidentMeshLod(stagedRequest, stagedProvider);
    stagedProvider.progressiveDelivery = FALSE;
    BObolLodResult terminalStage =
	service.realizeResidentMeshLod(stagedRequest, stagedProvider);
    /* Keep the terminal prefix resident but tell the provider that the view
     * is currently presenting the first stage.  Refinement must grow from
     * that visible prefix, not jump directly to the resident capacity. */
    stagedProvider.progressiveDelivery = TRUE;
    stagedProvider.useCurrentDrawLevel = TRUE;
    stagedProvider.currentDrawLevel = firstStage.geometry.activeLevel;
    BObolLodResult residentAheadStage =
	service.realizeResidentMeshLod(stagedRequest, stagedProvider);
    if (firstStage.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	firstStage.geometry.activeLevel < 0 ||
	firstStage.geometry.activeLevel >= stagedRequest.requestedLevel ||
	firstStage.terminal ||
	secondStage.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	secondStage.geometry.activeLevel <= firstStage.geometry.activeLevel ||
	secondStage.geometry.activeLevel != firstStage.geometry.activeLevel + 1 ||
	secondStage.geometry.activeLevel >= stagedRequest.requestedLevel ||
	secondStage.terminal ||
	terminalStage.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	terminalStage.geometry.activeLevel != stagedRequest.requestedLevel ||
	!terminalStage.terminal ||
	residentAheadStage.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	residentAheadStage.geometry.activeLevel <=
	    firstStage.geometry.activeLevel ||
	residentAheadStage.geometry.activeLevel !=
	    firstStage.geometry.activeLevel + 1 ||
	residentAheadStage.geometry.activeLevel >=
	    stagedRequest.requestedLevel ||
	residentAheadStage.terminal ||
	!firstStage.progressiveMesh ||
	firstStage.progressiveMesh != secondStage.progressiveMesh ||
	firstStage.progressiveMesh != terminalStage.progressiveMesh ||
	firstStage.progressiveMesh != residentAheadStage.progressiveMesh ||
	firstStage.counts.faceCount > secondStage.counts.faceCount ||
	secondStage.counts.faceCount > terminalStage.counts.faceCount) {
	printf("FAIL: retained LoD provider did not stage population growth "
	    "(levels=%d/%d/%d resident-ahead=%d "
	    "terminal=%d/%d/%d/%d faces=%llu/%llu/%llu/%llu)\n",
	    firstStage.geometry.activeLevel, secondStage.geometry.activeLevel,
	    terminalStage.geometry.activeLevel,
	    residentAheadStage.geometry.activeLevel,
	    firstStage.terminal ? 1 : 0,
	    secondStage.terminal ? 1 : 0, terminalStage.terminal ? 1 : 0,
	    residentAheadStage.terminal ? 1 : 0,
	    static_cast<unsigned long long>(firstStage.counts.faceCount),
	    static_cast<unsigned long long>(secondStage.counts.faceCount),
	    static_cast<unsigned long long>(terminalStage.counts.faceCount),
	    static_cast<unsigned long long>(
		residentAheadStage.counts.faceCount));
	ret = 1;
    }

    service.stop();

    BObolMeshLodProvider cachedNormalProvider;
    cachedNormalProvider.dbip = dbip;
    cachedNormalProvider.useForcedLevel = TRUE;
    cachedNormalProvider.forcedLevel = 1;
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
	cachedNormalResult.geometry.activeLevel != 1 ||
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
	    cachedNormalResult.geometry.activeLevel,
	    static_cast<unsigned long long>(
		cachedNormalResult.counts.normalCount),
	    cachedNormalResult.hasNormals ? 1 : 0,
	    cachedNormalResult.mesh.normals.size(),
	    cachedNormalResult.mesh.coordIndex.size(),
	    sawSeededNormal ? 1 : 0);
	ret = 1;
    }

    BObolMeshLodProvider forcedProvider;
    forcedProvider.dbip = dbip;
    forcedProvider.useForcedLevel = TRUE;
    forcedProvider.forcedLevel = (results.size() == 1 &&
				  results[0].geometry.activeLevel >= 0) ?
				 results[0].geometry.activeLevel : 0;
    forcedProvider.refreshMissing = FALSE;
    forcedProvider.shrinkAfterCopy = TRUE;

    BObolLodRequest forcedRequest = task.request;
    forcedRequest.qualityTier = BOBOL_LOD_QUALITY_FULL_DETAIL;
    BObolLodResult forcedResult =
	bobol_mesh_lod_provider_task(forcedRequest, &forcedProvider);
    if (forcedResult.resultKind != BOBOL_LOD_RESULT_MESH ||
	forcedResult.qualityTier != BOBOL_LOD_QUALITY_FULL_DETAIL ||
	forcedResult.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	forcedResult.geometry.activeLevel != forcedProvider.forcedLevel ||
	!forcedResult.mesh.isValid()) {
	printf("FAIL: LoD Obol mesh provider forced-level task did not return "
	    "requested level mesh result (requested=%d active=%d status=%d "
	    "valid=%d)\n", forcedProvider.forcedLevel,
	    forcedResult.geometry.activeLevel, forcedResult.providerStatus,
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
    staleProvider.dbip = dbip;
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
    refreshProvider.dbip = dbip;
    refreshProvider.useView = TRUE;
    refreshProvider.refreshMissing = TRUE;
    refreshProvider.view = provider.view;

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
    bobol_init(NULL);

    if (test_dependency_order_and_cache_write())
	return 1;
    if (test_filtered_result_drain())
	return 1;
    if (test_occurrence_result_coalescing())
	return 1;
    if (test_generation_cancellation())
	return 1;
    if (test_destruction_waits_for_in_flight_task())
	return 1;
    if (test_stop_discards_undrained_state())
	return 1;
    if (test_stale_result_rejection())
	return 1;
    if (test_staged_payload_delivery())
	return 1;
    if (test_result_ready_subscription())
	return 1;
    if (test_result_ready_self_unsubscribe())
	return 1;
    if (test_result_ready_cross_unsubscribe())
	return 1;
    if (test_task_realize_data_cleanup())
	return 1;
    if (test_queue_limits_and_pending_cancellation())
	return 1;
    if (test_large_pending_cancellation_and_generation_history())
	return 1;
    if (test_debug_delay_cancellation())
	return 1;
    if (test_cancelled_cache_write_not_persisted())
	return 1;
    if (test_active_request_duplicate_suppression())
	return 1;
    if (test_rt_source_full_detail_provider_task())
	return 1;
    if (test_rt_mesh_provider_task())
	return 1;
    if (test_detached_database_snapshot())
	return 1;
    if (test_rt_obb_proxy_provider_request_bounds())
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
