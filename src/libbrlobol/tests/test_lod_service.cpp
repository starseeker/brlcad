/*                T E S T _ L O D _ S E R V I C E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/lod_service.h"
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

static BRLObolLodRequest
make_request(const char *name)
{
    BRLObolLodRequest request;

    request.databaseId = "db://lod-service-test";
    request.databaseRevision = 9;
    request.sourceRevision = 17;
    request.sourceContentHash = 0x5577;
    request.objectPath = name;
    request.objectName = name;
    request.viewRevision = 23;
    request.policyRevision = 29;
    request.drawMode = BRLOBOL_LOD_DRAW_SHADED;
    request.lodPolicy = 3;
    request.providerId = "service-test";
    request.providerVersion = "1";
    request.qualityTier = BRLOBOL_LOD_QUALITY_FAST_DISPLAY;
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

static BRLObolLodResult
ready_task(const BRLObolLodRequest &request, void *userData)
{
    TaskData *data = static_cast<TaskData *>(userData);

    {
	std::lock_guard<std::mutex> lock(data->context->mutex);
	data->context->executionOrder.push_back(data->value);
    }

    BRLObolLodResult result;
    result.request = request;
    result.cacheKey = brlobol_lod_cache_key(request);
    result.resultKind = BRLOBOL_LOD_RESULT_AABB;
    result.qualityTier = request.qualityTier;
    result.providerStatus = BRLOBOL_LOD_PROVIDER_READY;
    result.terminal = TRUE;
    result.counts.faceCount = (uint64_t)data->value;
    result.bounds = request.bounds;
    return result;
}

static void
cache_write(const BRLObolLodResult &result, void *userData)
{
    ServiceTestContext *context = static_cast<ServiceTestContext *>(userData);
    std::lock_guard<std::mutex> lock(context->mutex);
    context->cacheWriteOrder.push_back((int)result.counts.faceCount);
}

static int
wait_for_settled(BRLObolLodService &service, size_t expectedResults)
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
provider_param_count(const BRLObolLodRequest &request, const char *name)
{
    size_t count = 0;

    if (!name)
	return 0;

    for (size_t i = 0; i < request.providerParams.size(); i++) {
	if (strcmp(request.providerParams[i].name.getString(), name) == 0)
	    count++;
    }

    return count;
}

static int
test_dependency_order_and_cache_write(void)
{
    BRLObolLodService service;
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

    BRLObolLodTask first;
    first.generation = generation;
    first.request = make_request("/first.bot");
    first.realize = ready_task;
    first.realizeData = &firstData;
    first.writeCache = TRUE;
    first.cacheWrite = cache_write;
    first.cacheWriteData = &context;

    BRLObolLodTask second;
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

    std::vector<BRLObolLodResult> results;
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
	!brlobol_lod_result_matches_request(results[0], first.request) ||
	!brlobol_lod_result_matches_request(results[1], second.request)) {
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
    BRLObolLodService service;

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
    BRLObolLodTask first;
    first.generation = generation;
    first.request = make_request("/filtered-first.bot");
    first.realize = ready_task;
    first.realizeData = &firstData;

    BRLObolLodTask second;
    second.generation = generation;
    second.request = make_request("/filtered-second.bot");
    second.realize = ready_task;
    second.realizeData = &secondData;

    BRLObolLodTask third;
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

    std::vector<BRLObolLodRequest> requests;
    requests.push_back(second.request);
    std::vector<BRLObolLodResult> matched;
    if (service.drainMatchingResults(matched, requests) != 1 ||
	matched.size() != 1 ||
	!brlobol_lod_result_matches_request(matched[0], second.request) ||
	service.queuedResultCountForDiagnostics() != 2) {
	printf("FAIL: LoD filtered result drain did not isolate requested result\n");
	service.stop();
	return 1;
    }

    std::vector<BRLObolLodResult> remaining;
    if (service.drainResults(remaining) != 2 ||
	remaining.size() != 2 ||
	!brlobol_lod_result_matches_request(remaining[0], first.request) ||
	!brlobol_lod_result_matches_request(remaining[1], third.request)) {
	printf("FAIL: LoD filtered result drain did not preserve unmatched queue order\n");
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

static BRLObolLodResult
blocking_task(const BRLObolLodRequest &request, void *userData)
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

    BRLObolLodResult result;
    result.request = request;
    result.cacheKey = brlobol_lod_cache_key(request);
    result.resultKind = BRLOBOL_LOD_RESULT_MESH;
    result.qualityTier = request.qualityTier;
    result.providerStatus = BRLOBOL_LOD_PROVIDER_READY;
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
    BRLObolLodService service;
    BlockingTaskData blockingData;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for cancellation test\n");
	return 1;
    }

    uint64_t generation = service.beginGeneration();
    BRLObolLodTask task;
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

    if (wait_for_settled(service, 1))
	return 1;

    std::vector<BRLObolLodResult> results;
    service.drainResults(results);
    if (results.size() != 1 ||
	results[0].providerStatus != BRLOBOL_LOD_PROVIDER_CANCELLED ||
	!results[0].stale ||
	results[0].diagnostic.getLength() == 0 ||
	results[0].counts.faceCount != 0) {
	printf("FAIL: LoD service did not publish cancellation result\n");
	return 1;
    }

    if (blockingData.calls != 1) {
	printf("FAIL: LoD service cancellation callback count mismatch\n");
	return 1;
    }

    service.stop();
    return 0;
}

static BRLObolLodResult
stale_task(const BRLObolLodRequest &request, void *UNUSED(userData))
{
    BRLObolLodRequest stale = request;
    stale.viewRevision++;

    BRLObolLodResult result;
    result.request = stale;
    result.cacheKey = brlobol_lod_cache_key(stale);
    result.resultKind = BRLOBOL_LOD_RESULT_MESH;
    result.qualityTier = stale.qualityTier;
    result.providerStatus = BRLOBOL_LOD_PROVIDER_READY;
    result.counts.faceCount = 11;
    return result;
}

static int
test_stale_result_rejection(void)
{
    BRLObolLodService service;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for stale-result test\n");
	return 1;
    }

    BRLObolLodTask task;
    task.generation = service.beginGeneration();
    task.request = make_request("/stale.bot");
    task.realize = stale_task;

    if (service.submit(task) == 0) {
	printf("FAIL: LoD service did not accept stale-result task\n");
	return 1;
    }

    if (wait_for_settled(service, 1))
	return 1;

    std::vector<BRLObolLodResult> results;
    service.drainResults(results);
    if (results.size() != 1 ||
	results[0].providerStatus != BRLOBOL_LOD_PROVIDER_STALE ||
	!results[0].stale ||
	!results[0].terminal ||
	results[0].diagnostic.getLength() == 0) {
	printf("FAIL: LoD service did not reject stale result\n");
	return 1;
    }

    service.stop();
    return 0;
}

static BRLObolLodResult
staged_payload_task(const BRLObolLodRequest &request, void *UNUSED(userData))
{
    std::vector<BRLObolLodAttribute> attributes;
    BRLObolLodAttribute attribute;

    attribute.name = "display.intent";
    attribute.value = "proxy";
    attributes.push_back(attribute);

    return brlobol_lod_attributes_result(request, attributes);
}

static int
test_staged_payload_delivery(void)
{
    BRLObolLodService service;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for staged-payload test\n");
	return 1;
    }

    BRLObolLodTask task;
    task.generation = service.beginGeneration();
    task.request = make_request("/attributes.bot");
    task.realize = staged_payload_task;

    if (service.submit(task) == 0) {
	printf("FAIL: LoD service did not accept staged-payload task\n");
	return 1;
    }

    if (wait_for_settled(service, 1))
	return 1;

    std::vector<BRLObolLodResult> results;
    service.drainResults(results);
    if (results.size() != 1 ||
	results[0].resultKind != BRLOBOL_LOD_RESULT_ATTRIBUTES ||
	results[0].qualityTier != BRLOBOL_LOD_QUALITY_ATTRIBUTES ||
	results[0].attributes.size() != 1 ||
	strcmp(results[0].attributes[0].value.getString(), "proxy") != 0 ||
	!brlobol_lod_result_matches_request(results[0], task.request)) {
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
result_ready_cb(BRLObolLodService *service, void *userData)
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
    BRLObolLodService service;
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

    BRLObolLodSubscriberId subscriber =
	service.subscribeResultReady(result_ready_cb, &readyContext);
    if (subscriber == 0) {
	printf("FAIL: LoD service did not create result-ready subscription\n");
	service.stop();
	return 1;
    }

    uint64_t generation = service.beginGeneration();
    BRLObolLodTask first;
    first.generation = generation;
    first.request = make_request("/ready-first.bot");
    first.realize = ready_task;
    first.realizeData = &firstData;
    BRLObolLodTask second;
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

    std::vector<BRLObolLodResult> results;
    if (service.drainResults(results) != 2 ||
	service.queuedResultCountForDiagnostics() != 0) {
	printf("FAIL: LoD result-ready test did not drain burst results\n");
	service.unsubscribeResultReady(subscriber);
	service.stop();
	return 1;
    }

    BRLObolLodTask third;
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

    BRLObolLodTask fourth;
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

struct CleanupMarker {
    std::atomic<int> *calls;
    std::atomic<int> *cleanups;
};

static BRLObolLodResult
cleanup_task(const BRLObolLodRequest &request, void *userData)
{
    CleanupMarker *marker = static_cast<CleanupMarker *>(userData);
    marker->calls->fetch_add(1);

    BRLObolLodResult result;
    result.request = request;
    result.cacheKey = brlobol_lod_cache_key(request);
    result.resultKind = BRLOBOL_LOD_RESULT_AABB;
    result.qualityTier = request.qualityTier;
    result.providerStatus = BRLOBOL_LOD_PROVIDER_READY;
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
    BRLObolLodService service;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for cleanup test\n");
	return 1;
    }

    BRLObolLodTask task;
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

    std::vector<BRLObolLodResult> results;
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
    BRLObolLodService pendingService;
    if (!pendingService.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for pending-cleanup test\n");
	return 1;
    }

    BRLObolLodTask pendingTask;
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

struct DebugDelayTaskData {
    std::mutex mutex;
    int calls;

    DebugDelayTaskData(void) : calls(0)
    {
    }
};

static BRLObolLodResult
debug_delay_task(const BRLObolLodRequest &request, void *userData)
{
    DebugDelayTaskData *data = static_cast<DebugDelayTaskData *>(userData);

    {
	std::lock_guard<std::mutex> lock(data->mutex);
	data->calls++;
    }

    BRLObolLodResult result;
    result.request = request;
    result.cacheKey = brlobol_lod_cache_key(request);
    result.resultKind = BRLOBOL_LOD_RESULT_AABB;
    result.qualityTier = request.qualityTier;
    result.providerStatus = BRLOBOL_LOD_PROVIDER_READY;
    result.terminal = TRUE;
    return result;
}

static int
wait_for_debug_delay(BRLObolLodService &service)
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
    BRLObolLodService service;
    DebugDelayTaskData data;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for debug-delay test\n");
	return 1;
    }

    uint64_t generation = service.beginGeneration();
    BRLObolLodTask task;
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

    if (wait_for_settled(service, 1))
	return 1;

    std::vector<BRLObolLodResult> results;
    service.drainResults(results);
    if (results.size() != 1 ||
	results[0].providerStatus != BRLOBOL_LOD_PROVIDER_CANCELLED ||
	!results[0].stale ||
	results[0].diagnostic.getLength() == 0 ||
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

static int
make_provider_test_db(char *dbpath, size_t dbpath_len, struct db_i **dbip_out)
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
	printf("FAIL: LoD RT provider temp file\n");
	return 1;
    }
    fclose(fp);

    struct db_i *dbip = db_create(dbpath, 5);
    if (!dbip) {
	printf("FAIL: LoD RT provider db_create\n");
	bu_file_delete(dbpath);
	return 1;
    }

    struct rt_wdb *wdbp = wdb_dbopen(dbip, RT_WDB_TYPE_DB_DISK);
    if (!wdbp) {
	printf("FAIL: LoD RT provider wdb_dbopen\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    if (mk_bot(wdbp, objname, RT_BOT_SOLID, RT_BOT_UNORIENTED, 0,
	       4, 4, vertices, faces, NULL, NULL) != 0) {
	printf("FAIL: LoD RT provider mk_bot\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }
    if (mk_bot(wdbp, tri_objname, RT_BOT_SOLID, RT_BOT_UNORIENTED, 0,
	       6, 2, tri_vertices, tri_faces, NULL, NULL) != 0) {
	printf("FAIL: LoD RT provider mk_bot disjoint triangles\n");
	db_close(dbip);
	bu_file_delete(dbpath);
	return 1;
    }

    *dbip_out = dbip;
    return 0;
}

static int
check_source_full_detail_result(const BRLObolLodResult &result,
				const BRLObolLodRequest &request, const char *label)
{
    if (result.resultKind != BRLOBOL_LOD_RESULT_FULL_DETAIL ||
	result.qualityTier != BRLOBOL_LOD_QUALITY_FULL_DETAIL ||
	result.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	result.geometry.kind != BRLOBOL_LOD_GEOMETRY_OBOL_MESH ||
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
	!brlobol_lod_result_matches_request(result, request)) {
	printf("FAIL: LoD RT source full-detail provider %s\n", label);
	return 1;
    }

    return 0;
}

static int
test_active_request_duplicate_suppression(void)
{
    BRLObolLodService service;
    DebugDelayTaskData data;

    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD service did not start for active-request test\n");
	return 1;
    }

    BRLObolLodTask task;
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

    std::vector<BRLObolLodResult> results;
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

    BRLObolRtSourceFullDetailProvider provider;
    provider.dbip = dbip;
    provider.validateSourceMetrics = TRUE;

    BRLObolLodRequest request = make_request("/lod-provider.bot");
    request.objectName = "lod-provider.bot";
    request.providerId = "rt_source_full_detail";
    request.providerVersion = "direct-bot-v1";
    request.qualityTier = BRLOBOL_LOD_QUALITY_FULL_DETAIL;
    request.sourceCounts.faceCount = 4;
    request.sourceCounts.pointCount = 4;

    int ret = 0;
    BRLObolLodResult directResult =
	brlobol_rt_source_full_detail_provider_task(request, &provider);
    if (check_source_full_detail_result(directResult, request,
					"did not return direct BoT mesh result"))
	ret = 1;

    BRLObolLodRequest scopedRequest = request;
    scopedRequest.addProviderParam("source_query.space", "source_local");
    scopedRequest.addProviderParam("source_query.bounds",
				   "0.7 0.15 -0.05 0.8 0.25 0.05");
    scopedRequest.addProviderParam("source_query.tolerance", "0.05");
    BRLObolLodResult scopedResult =
	brlobol_rt_source_full_detail_provider_task(scopedRequest, &provider);
    if (scopedResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	scopedResult.resultKind != BRLOBOL_LOD_RESULT_FULL_DETAIL ||
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
	!brlobol_lod_result_matches_request(scopedResult, scopedRequest)) {
	printf("FAIL: LoD RT source full-detail provider query bounds did not reduce returned face payload\n");
	ret = 1;
    }

    BRLObolRtSourceFullDetailProvider scopedLimitProvider;
    scopedLimitProvider.dbip = dbip;
    scopedLimitProvider.validateSourceMetrics = TRUE;
    scopedLimitProvider.maxFullDetailFaceCount = 2;
    scopedLimitProvider.maxFullDetailPointCount = 4;
    BRLObolLodResult scopedLimitedResult =
	brlobol_rt_source_full_detail_provider_task(scopedRequest,
	    &scopedLimitProvider);
    if (scopedLimitedResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	scopedLimitedResult.resultKind != BRLOBOL_LOD_RESULT_FULL_DETAIL ||
	scopedLimitedResult.mesh.points.size() != 4 ||
	scopedLimitedResult.mesh.coordIndex.size() != 6 ||
	scopedLimitedResult.mesh.faceIndex.size() != 2 ||
	scopedLimitedResult.counts.faceCount != 2 ||
	scopedLimitedResult.counts.pointCount != 4 ||
	!brlobol_lod_result_matches_request(scopedLimitedResult,
					    scopedRequest)) {
	printf("FAIL: LoD RT source full-detail provider scoped bounds limit should apply after payload reduction\n");
	ret = 1;
    }

    BRLObolLodRequest emptyScopedRequest = request;
    emptyScopedRequest.addProviderParam("source_query.space", "source_local");
    emptyScopedRequest.addProviderParam("source_query.bounds",
					"5.0 5.0 5.0 5.1 5.1 5.1");
    emptyScopedRequest.addProviderParam("source_query.tolerance", "0.05");
    BRLObolLodResult emptyScopedResult =
	brlobol_rt_source_full_detail_provider_task(emptyScopedRequest,
	    &provider);
    if (emptyScopedResult.providerStatus != BRLOBOL_LOD_PROVIDER_FALLBACK ||
	emptyScopedResult.resultKind != BRLOBOL_LOD_RESULT_NONE ||
	emptyScopedResult.counts.faceCount != 0 ||
	emptyScopedResult.counts.pointCount != 0 ||
	emptyScopedResult.mesh.isValid() ||
	strcmp(emptyScopedResult.diagnostic.getString(),
	       "RT source full-detail provider scoped query matched no faces") != 0 ||
	!brlobol_lod_result_matches_request(emptyScopedResult,
					    emptyScopedRequest)) {
	printf("FAIL: LoD RT source full-detail provider scoped bounds miss should not expand to whole-object payload\n");
	ret = 1;
    }

    BRLObolLodRequest wrongSpaceBoundsRequest = request;
    wrongSpaceBoundsRequest.addProviderParam("source_query.space", "world");
    wrongSpaceBoundsRequest.addProviderParam("source_query.bounds",
	    "0.7 0.15 -0.05 0.8 0.25 0.05");
    wrongSpaceBoundsRequest.addProviderParam("source_query.tolerance",
	    "0.05");
    BRLObolLodResult wrongSpaceBoundsResult =
	brlobol_rt_source_full_detail_provider_task(wrongSpaceBoundsRequest,
	    &provider);
    if (wrongSpaceBoundsResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	wrongSpaceBoundsResult.resultKind !=
	BRLOBOL_LOD_RESULT_FULL_DETAIL ||
	wrongSpaceBoundsResult.mesh.points.size() != 4 ||
	wrongSpaceBoundsResult.mesh.coordIndex.size() != 12 ||
	wrongSpaceBoundsResult.counts.faceCount != 4 ||
	wrongSpaceBoundsResult.counts.pointCount != 4 ||
	!brlobol_lod_result_matches_request(wrongSpaceBoundsResult,
					    wrongSpaceBoundsRequest)) {
	printf("FAIL: LoD RT source full-detail provider should ignore non-source-local bounds when reducing payloads\n");
	ret = 1;
    }

    BRLObolRtSourceFullDetailProvider wrongSpaceLimitProvider;
    wrongSpaceLimitProvider.dbip = dbip;
    wrongSpaceLimitProvider.validateSourceMetrics = TRUE;
    wrongSpaceLimitProvider.maxFullDetailFaceCount = 2;
    wrongSpaceLimitProvider.maxFullDetailPointCount = 4;
    BRLObolLodResult wrongSpaceLimitedResult =
	brlobol_rt_source_full_detail_provider_task(wrongSpaceBoundsRequest,
	    &wrongSpaceLimitProvider);
    if (wrongSpaceLimitedResult.providerStatus !=
	BRLOBOL_LOD_PROVIDER_FALLBACK ||
	wrongSpaceLimitedResult.resultKind != BRLOBOL_LOD_RESULT_NONE ||
	wrongSpaceLimitedResult.mesh.isValid() ||
	strcmp(wrongSpaceLimitedResult.diagnostic.getString(),
	       "RT source full-detail provider request exceeds full-detail limits") != 0 ||
	!brlobol_lod_result_matches_request(wrongSpaceLimitedResult,
					    wrongSpaceBoundsRequest)) {
	printf("FAIL: LoD RT source full-detail provider should not bypass whole-object limits for non-source-local bounds\n");
	ret = 1;
    }

    BRLObolLodRequest malformedBoundsRequest = request;
    malformedBoundsRequest.addProviderParam("source_query.space",
					    "source_local");
    malformedBoundsRequest.addProviderParam("source_query.bounds",
					    "0.7 0.15 -0.05 0.8 0.25 0.05 trailing");
    malformedBoundsRequest.addProviderParam("source_query.tolerance",
					    "0.05");
    BRLObolLodResult malformedBoundsResult =
	brlobol_rt_source_full_detail_provider_task(malformedBoundsRequest,
	    &provider);
    if (malformedBoundsResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	malformedBoundsResult.resultKind !=
	BRLOBOL_LOD_RESULT_FULL_DETAIL ||
	malformedBoundsResult.mesh.points.size() != 4 ||
	malformedBoundsResult.mesh.coordIndex.size() != 12 ||
	malformedBoundsResult.counts.faceCount != 4 ||
	malformedBoundsResult.counts.pointCount != 4 ||
	!brlobol_lod_result_matches_request(malformedBoundsResult,
					    malformedBoundsRequest)) {
	printf("FAIL: LoD RT source full-detail provider should ignore malformed source-local bounds when reducing payloads\n");
	ret = 1;
    }

    BRLObolLodRequest malformedToleranceRequest = request;
    malformedToleranceRequest.addProviderParam("source_query.space",
	    "source_local");
    malformedToleranceRequest.addProviderParam("source_query.bounds",
	    "0.7 0.15 -0.05 0.8 0.25 0.05");
    malformedToleranceRequest.addProviderParam("source_query.tolerance",
	    "0.05 trailing");
    BRLObolLodResult malformedToleranceLimitedResult =
	brlobol_rt_source_full_detail_provider_task(malformedToleranceRequest,
	    &wrongSpaceLimitProvider);
    if (malformedToleranceLimitedResult.providerStatus !=
	BRLOBOL_LOD_PROVIDER_FALLBACK ||
	malformedToleranceLimitedResult.resultKind !=
	BRLOBOL_LOD_RESULT_NONE ||
	malformedToleranceLimitedResult.mesh.isValid() ||
	strcmp(malformedToleranceLimitedResult.diagnostic.getString(),
	       "RT source full-detail provider request exceeds full-detail limits") != 0 ||
	!brlobol_lod_result_matches_request(
	    malformedToleranceLimitedResult, malformedToleranceRequest)) {
	printf("FAIL: LoD RT source full-detail provider should not bypass whole-object limits for malformed source-local tolerance\n");
	ret = 1;
    }

    BRLObolLodRequest duplicateSpaceBoundsRequest = request;
    duplicateSpaceBoundsRequest.addProviderParam("source_query.space",
	    "source_local");
    duplicateSpaceBoundsRequest.addProviderParam("source_query.space",
	    "source_local");
    duplicateSpaceBoundsRequest.addProviderParam("source_query.bounds",
	    "0.7 0.15 -0.05 0.8 0.25 0.05");
    duplicateSpaceBoundsRequest.addProviderParam("source_query.tolerance",
	    "0.05");
    BRLObolLodResult duplicateSpaceBoundsResult =
	brlobol_rt_source_full_detail_provider_task(
	    duplicateSpaceBoundsRequest, &provider);
    if (duplicateSpaceBoundsResult.providerStatus !=
	BRLOBOL_LOD_PROVIDER_READY ||
	duplicateSpaceBoundsResult.resultKind !=
	BRLOBOL_LOD_RESULT_FULL_DETAIL ||
	duplicateSpaceBoundsResult.mesh.points.size() != 4 ||
	duplicateSpaceBoundsResult.mesh.coordIndex.size() != 12 ||
	duplicateSpaceBoundsResult.counts.faceCount != 4 ||
	duplicateSpaceBoundsResult.counts.pointCount != 4 ||
	!brlobol_lod_result_matches_request(duplicateSpaceBoundsResult,
					    duplicateSpaceBoundsRequest)) {
	printf("FAIL: LoD RT source full-detail provider should ignore duplicate query-space params when reducing payloads\n");
	ret = 1;
    }

    BRLObolLodRequest duplicateBoundsRequest = request;
    duplicateBoundsRequest.addProviderParam("source_query.space",
					    "source_local");
    duplicateBoundsRequest.addProviderParam("source_query.bounds",
					    "0.7 0.15 -0.05 0.8 0.25 0.05");
    duplicateBoundsRequest.addProviderParam("source_query.bounds",
					    "0.7 0.15 -0.05 0.8 0.25 0.05");
    duplicateBoundsRequest.addProviderParam("source_query.tolerance",
					    "0.05");
    BRLObolLodResult duplicateBoundsLimitedResult =
	brlobol_rt_source_full_detail_provider_task(duplicateBoundsRequest,
	    &wrongSpaceLimitProvider);
    if (duplicateBoundsLimitedResult.providerStatus !=
	BRLOBOL_LOD_PROVIDER_FALLBACK ||
	duplicateBoundsLimitedResult.resultKind !=
	BRLOBOL_LOD_RESULT_NONE ||
	duplicateBoundsLimitedResult.mesh.isValid() ||
	strcmp(duplicateBoundsLimitedResult.diagnostic.getString(),
	       "RT source full-detail provider request exceeds full-detail limits") != 0 ||
	!brlobol_lod_result_matches_request(duplicateBoundsLimitedResult,
					    duplicateBoundsRequest)) {
	printf("FAIL: LoD RT source full-detail provider should not bypass whole-object limits for duplicate bounds params\n");
	ret = 1;
    }

    BRLObolLodRequest mixedScopedRequest = request;
    mixedScopedRequest.addProviderParam("source_query.space",
					"source_local");
    mixedScopedRequest.addProviderParam("source_query.bounds",
					"0.7 0.15 -0.05 0.8 0.25 0.05");
    mixedScopedRequest.addProviderParam("source_query.tolerance", "0.05");
    mixedScopedRequest.addProviderParam("source_query.ray.origin",
					"0.2 0.2 5.0");
    mixedScopedRequest.addProviderParam("source_query.ray.direction",
					"0.0 0.0 -1.0");
    BRLObolLodResult mixedScopedResult =
	brlobol_rt_source_full_detail_provider_task(mixedScopedRequest,
	    &provider);
    if (mixedScopedResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	mixedScopedResult.resultKind != BRLOBOL_LOD_RESULT_FULL_DETAIL ||
	mixedScopedResult.mesh.points.size() != 4 ||
	mixedScopedResult.mesh.coordIndex.size() != 12 ||
	mixedScopedResult.counts.faceCount != 4 ||
	mixedScopedResult.counts.pointCount != 4 ||
	!brlobol_lod_result_matches_request(mixedScopedResult,
					    mixedScopedRequest)) {
	printf("FAIL: LoD RT source full-detail provider should ignore mixed scoped query kinds when reducing payloads\n");
	ret = 1;
    }

    BRLObolLodResult mixedScopedLimitedResult =
	brlobol_rt_source_full_detail_provider_task(mixedScopedRequest,
	    &wrongSpaceLimitProvider);
    if (mixedScopedLimitedResult.providerStatus !=
	BRLOBOL_LOD_PROVIDER_FALLBACK ||
	mixedScopedLimitedResult.resultKind != BRLOBOL_LOD_RESULT_NONE ||
	mixedScopedLimitedResult.mesh.isValid() ||
	strcmp(mixedScopedLimitedResult.diagnostic.getString(),
	       "RT source full-detail provider request exceeds full-detail limits") != 0 ||
	!brlobol_lod_result_matches_request(mixedScopedLimitedResult,
					    mixedScopedRequest)) {
	printf("FAIL: LoD RT source full-detail provider should not bypass whole-object limits for mixed scoped query kinds\n");
	ret = 1;
    }

    BRLObolLodRequest rayRequest = request;
    rayRequest.addProviderParam("source_query.space", "source_local");
    rayRequest.addProviderParam("source_query.ray.origin",
				"0.2 0.2 5.0");
    rayRequest.addProviderParam("source_query.ray.direction",
				"0.0 0.0 -1.0");
    BRLObolLodResult rayResult =
	brlobol_rt_source_full_detail_provider_task(rayRequest, &provider);
    if (rayResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	rayResult.resultKind != BRLOBOL_LOD_RESULT_FULL_DETAIL ||
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
	!brlobol_lod_result_matches_request(rayResult, rayRequest)) {
	printf("FAIL: LoD RT source full-detail provider query ray did not reduce returned face payload\n");
	ret = 1;
    }

    BRLObolLodRequest compactRayRequest = make_request("/lod-two-tri.bot");
    compactRayRequest.objectName = "lod-two-tri.bot";
    compactRayRequest.providerId = "rt_source_full_detail";
    compactRayRequest.providerVersion = "direct-bot-v1";
    compactRayRequest.qualityTier = BRLOBOL_LOD_QUALITY_FULL_DETAIL;
    compactRayRequest.sourceCounts.faceCount = 2;
    compactRayRequest.sourceCounts.pointCount = 6;
    compactRayRequest.addProviderParam("source_query.space", "source_local");
    compactRayRequest.addProviderParam("source_query.ray.origin",
				       "0.2 0.2 5.0");
    compactRayRequest.addProviderParam("source_query.ray.direction",
				       "0.0 0.0 -1.0");
    BRLObolLodResult compactRayResult =
	brlobol_rt_source_full_detail_provider_task(compactRayRequest,
	    &provider);
    if (compactRayResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	compactRayResult.resultKind != BRLOBOL_LOD_RESULT_FULL_DETAIL ||
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
	!brlobol_lod_result_matches_request(compactRayResult,
					    compactRayRequest)) {
	printf("FAIL: LoD RT source full-detail provider query ray did not compact source vertex payload\n");
	ret = 1;
    }

    BRLObolRtSourceFullDetailProvider compactRayLimitProvider;
    compactRayLimitProvider.dbip = dbip;
    compactRayLimitProvider.validateSourceMetrics = TRUE;
    compactRayLimitProvider.maxFullDetailFaceCount = 1;
    compactRayLimitProvider.maxFullDetailPointCount = 3;
    BRLObolLodResult compactRayLimitedResult =
	brlobol_rt_source_full_detail_provider_task(compactRayRequest,
	    &compactRayLimitProvider);
    if (compactRayLimitedResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	compactRayLimitedResult.resultKind != BRLOBOL_LOD_RESULT_FULL_DETAIL ||
	compactRayLimitedResult.mesh.points.size() != 3 ||
	compactRayLimitedResult.mesh.coordIndex.size() != 3 ||
	compactRayLimitedResult.mesh.faceIndex.size() != 1 ||
	compactRayLimitedResult.mesh.vertexIndex.size() != 3 ||
	compactRayLimitedResult.counts.faceCount != 1 ||
	compactRayLimitedResult.counts.pointCount != 3 ||
	!brlobol_lod_result_matches_request(compactRayLimitedResult,
					    compactRayRequest)) {
	printf("FAIL: LoD RT source full-detail provider scoped ray limit should apply after payload reduction\n");
	ret = 1;
    }

    BRLObolLodRequest missRayRequest = compactRayRequest;
    missRayRequest.providerParams.clear();
    missRayRequest.addProviderParam("source_query.space", "source_local");
    missRayRequest.addProviderParam("source_query.ray.origin",
				    "5.0 5.0 5.0");
    missRayRequest.addProviderParam("source_query.ray.direction",
				    "0.0 0.0 -1.0");
    BRLObolLodResult missRayResult =
	brlobol_rt_source_full_detail_provider_task(missRayRequest, &provider);
    if (missRayResult.providerStatus != BRLOBOL_LOD_PROVIDER_FALLBACK ||
	missRayResult.resultKind != BRLOBOL_LOD_RESULT_NONE ||
	missRayResult.counts.faceCount != 0 ||
	missRayResult.counts.pointCount != 0 ||
	missRayResult.mesh.isValid() ||
	strcmp(missRayResult.diagnostic.getString(),
	       "RT source full-detail provider scoped query matched no faces") != 0 ||
	!brlobol_lod_result_matches_request(missRayResult,
					    missRayRequest)) {
	printf("FAIL: LoD RT source full-detail provider scoped ray miss should not expand to whole-object payload\n");
	ret = 1;
    }

    BRLObolLodRequest wrongSpaceRayRequest = compactRayRequest;
    wrongSpaceRayRequest.providerParams.clear();
    wrongSpaceRayRequest.addProviderParam("source_query.space", "world");
    wrongSpaceRayRequest.addProviderParam("source_query.ray.origin",
					  "0.2 0.2 5.0");
    wrongSpaceRayRequest.addProviderParam("source_query.ray.direction",
					  "0.0 0.0 -1.0");
    BRLObolLodResult wrongSpaceRayResult =
	brlobol_rt_source_full_detail_provider_task(wrongSpaceRayRequest,
	    &provider);
    if (wrongSpaceRayResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	wrongSpaceRayResult.resultKind != BRLOBOL_LOD_RESULT_FULL_DETAIL ||
	wrongSpaceRayResult.mesh.points.size() != 6 ||
	wrongSpaceRayResult.mesh.coordIndex.size() != 6 ||
	!wrongSpaceRayResult.mesh.vertexIndex.empty() ||
	wrongSpaceRayResult.counts.faceCount != 2 ||
	wrongSpaceRayResult.counts.pointCount != 6 ||
	!brlobol_lod_result_matches_request(wrongSpaceRayResult,
					    wrongSpaceRayRequest)) {
	printf("FAIL: LoD RT source full-detail provider should ignore non-source-local rays when reducing payloads\n");
	ret = 1;
    }

    BRLObolLodRequest malformedRayRequest = compactRayRequest;
    malformedRayRequest.providerParams.clear();
    malformedRayRequest.addProviderParam("source_query.space",
					 "source_local");
    malformedRayRequest.addProviderParam("source_query.ray.origin",
					 "0.2 0.2 5.0");
    malformedRayRequest.addProviderParam("source_query.ray.direction",
					 "0.0 0.0 -1.0 trailing");
    BRLObolLodResult malformedRayResult =
	brlobol_rt_source_full_detail_provider_task(malformedRayRequest,
	    &provider);
    if (malformedRayResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	malformedRayResult.resultKind != BRLOBOL_LOD_RESULT_FULL_DETAIL ||
	malformedRayResult.mesh.points.size() != 6 ||
	malformedRayResult.mesh.coordIndex.size() != 6 ||
	!malformedRayResult.mesh.vertexIndex.empty() ||
	malformedRayResult.counts.faceCount != 2 ||
	malformedRayResult.counts.pointCount != 6 ||
	!brlobol_lod_result_matches_request(malformedRayResult,
					    malformedRayRequest)) {
	printf("FAIL: LoD RT source full-detail provider should ignore malformed source-local rays when reducing payloads\n");
	ret = 1;
    }

    BRLObolLodRequest duplicateRayRequest = compactRayRequest;
    duplicateRayRequest.providerParams.clear();
    duplicateRayRequest.addProviderParam("source_query.space",
					 "source_local");
    duplicateRayRequest.addProviderParam("source_query.ray.origin",
					 "0.2 0.2 5.0");
    duplicateRayRequest.addProviderParam("source_query.ray.direction",
					 "0.0 0.0 -1.0");
    duplicateRayRequest.addProviderParam("source_query.ray.direction",
					 "0.0 0.0 -1.0");
    BRLObolLodResult duplicateRayResult =
	brlobol_rt_source_full_detail_provider_task(duplicateRayRequest,
	    &provider);
    if (duplicateRayResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	duplicateRayResult.resultKind != BRLOBOL_LOD_RESULT_FULL_DETAIL ||
	duplicateRayResult.mesh.points.size() != 6 ||
	duplicateRayResult.mesh.coordIndex.size() != 6 ||
	!duplicateRayResult.mesh.vertexIndex.empty() ||
	duplicateRayResult.counts.faceCount != 2 ||
	duplicateRayResult.counts.pointCount != 6 ||
	!brlobol_lod_result_matches_request(duplicateRayResult,
					    duplicateRayRequest)) {
	printf("FAIL: LoD RT source full-detail provider should ignore duplicate ray params when reducing payloads\n");
	ret = 1;
    }

    BRLObolLodRequest measureHintRequest = request;
    measureHintRequest.addProviderParam("source_query.space", "source_local");
    measureHintRequest.addProviderParam("source_query.bounds",
					"0.7 0.15 -0.05 0.8 0.25 0.05");
    BRLObolLodResult measureHintResult =
	brlobol_rt_source_full_detail_provider_task(measureHintRequest, &provider);
    if (measureHintResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	measureHintResult.mesh.coordIndex.size() != 12 ||
	measureHintResult.mesh.faceIndex.size() != 4 ||
	!measureHintResult.mesh.vertexIndex.empty() ||
	measureHintResult.counts.faceCount != 4 ||
	!brlobol_lod_result_matches_request(measureHintResult,
					    measureHintRequest)) {
	printf("FAIL: LoD RT source full-detail provider should keep measure query hints whole-object without tolerance\n");
	ret = 1;
    }

    BRLObolLodRequest boundedMeasureRequest = request;
    boundedMeasureRequest.addProviderParam("source_query.space",
					   "source_local");
    boundedMeasureRequest.addProviderParam("source_query.bounds",
					   "0.7 0.15 -0.05 0.8 0.25 0.05");
    boundedMeasureRequest.addProviderParam("source_query.tolerance", "0.05");
    BRLObolLodResult boundedMeasureResult =
	brlobol_rt_source_full_detail_provider_task(boundedMeasureRequest,
	    &provider);
    if (boundedMeasureResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	boundedMeasureResult.mesh.points.size() != 4 ||
	boundedMeasureResult.mesh.coordIndex.size() != 6 ||
	boundedMeasureResult.mesh.faceIndex.size() != 2 ||
	boundedMeasureResult.mesh.vertexIndex.size() != 4 ||
	boundedMeasureResult.counts.faceCount != 2 ||
	boundedMeasureResult.counts.pointCount != 4 ||
	!brlobol_lod_result_matches_request(boundedMeasureResult,
					    boundedMeasureRequest)) {
	printf("FAIL: LoD RT source full-detail provider should reduce explicit bounded measure query payloads\n");
	ret = 1;
    }

    BRLObolLodService service;
    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD RT source full-detail provider service did not start\n");
	ret = 1;
    } else {
	BRLObolLodTask task;
	task.generation = service.beginGeneration();
	task.request = request;
	task.realize = brlobol_rt_source_full_detail_provider_task;
	task.realizeData = &provider;

	if (service.submit(task) == 0) {
	    printf("FAIL: LoD RT source full-detail provider service did not accept task\n");
	    ret = 1;
	} else if (wait_for_settled(service, 1)) {
	    ret = 1;
	} else {
	    std::vector<BRLObolLodResult> results;
	    service.drainResults(results);
	    if (results.size() != 1 ||
		check_source_full_detail_result(results[0], request,
						"did not publish service BoT mesh result"))
		ret = 1;
	}
	service.stop();
    }

    BRLObolLodRequest staleRequest = request;
    staleRequest.sourceCounts.faceCount = 99;
    BRLObolLodResult staleResult =
	brlobol_rt_source_full_detail_provider_task(staleRequest, &provider);
    if (staleResult.providerStatus != BRLOBOL_LOD_PROVIDER_STALE ||
	!staleResult.stale ||
	strcmp(staleResult.diagnostic.getString(),
	       "RT source full-detail provider source metrics changed") != 0 ||
	staleResult.mesh.isValid()) {
	printf("FAIL: LoD RT source full-detail provider did not reject stale source metrics\n");
	ret = 1;
    }

    BRLObolSourceMeshRequest sourceRequest;
    sourceRequest.path = "/lod-provider.bot";
    sourceRequest.sourceName = "lod-provider.bot";
    sourceRequest.sourceType = "bot";
    sourceRequest.sourceId = 7001;
    sourceRequest.faceCount = 4;
    sourceRequest.pointCount = 4;
    sourceRequest.bounds = request.bounds;

    BRLObolLodRequest convertedRequest;
    if (!brlobol_lod_rt_source_full_detail_request_from_source_mesh_request(
	    convertedRequest, sourceRequest, &request) ||
	strcmp(convertedRequest.providerId.getString(),
	       "rt_source_full_detail") != 0 ||
	strcmp(convertedRequest.providerVersion.getString(),
	       "direct-bot-v1") != 0 ||
	convertedRequest.qualityTier != BRLOBOL_LOD_QUALITY_FULL_DETAIL ||
	convertedRequest.sourceCounts.faceCount != 4 ||
	convertedRequest.sourceCounts.pointCount != 4 ||
	strcmp(convertedRequest.objectPath.getString(),
	       "/lod-provider.bot") != 0 ||
	strcmp(convertedRequest.objectName.getString(),
	       "lod-provider.bot") != 0) {
	printf("FAIL: LoD RT source full-detail helper did not convert source request\n");
	ret = 1;
    }

    BRLObolSourceMeshRequest templatedScopedSourceRequest = sourceRequest;
    templatedScopedSourceRequest.queryBoundsValid = 1;
    templatedScopedSourceRequest.queryBounds = SbBox3f(
	    SbVec3f(0.7f, 0.15f, -0.05f),
	    SbVec3f(0.8f, 0.25f, 0.05f));
    templatedScopedSourceRequest.queryToleranceValid = 1;
    templatedScopedSourceRequest.queryTolerance = 0.05f;
    BRLObolLodRequest staleTemplateRequest = request;
    staleTemplateRequest.addProviderParam("source_query.space", "world");
    staleTemplateRequest.addProviderParam("source_query.bounds",
					  "9.0 9.0 9.0 10.0 10.0 10.0");
    staleTemplateRequest.addProviderParam("source_query.ray.origin",
					  "0.2 0.2 5.0");
    staleTemplateRequest.addProviderParam("source_query.ray.direction",
					  "0.0 0.0 -1.0");
    staleTemplateRequest.addProviderParam("provider-template", "kept");

    BRLObolLodRequest cleanTemplatedRequest;
    if (!brlobol_lod_rt_source_full_detail_request_from_source_mesh_request(
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
	BRLObolRtSourceFullDetailProvider templatedLimitProvider;
	templatedLimitProvider.dbip = dbip;
	templatedLimitProvider.validateSourceMetrics = TRUE;
	templatedLimitProvider.maxFullDetailFaceCount = 2;
	templatedLimitProvider.maxFullDetailPointCount = 4;
	BRLObolLodResult cleanTemplatedResult =
	    brlobol_rt_source_full_detail_provider_task(
		cleanTemplatedRequest, &templatedLimitProvider);
	if (cleanTemplatedResult.providerStatus !=
	    BRLOBOL_LOD_PROVIDER_READY ||
	    cleanTemplatedResult.resultKind !=
	    BRLOBOL_LOD_RESULT_FULL_DETAIL ||
	    cleanTemplatedResult.mesh.coordIndex.size() != 6 ||
	    cleanTemplatedResult.mesh.faceIndex.size() != 2 ||
	    cleanTemplatedResult.counts.faceCount != 2 ||
	    cleanTemplatedResult.counts.pointCount != 4 ||
	    !brlobol_lod_result_matches_request(cleanTemplatedResult,
						cleanTemplatedRequest)) {
	    printf("FAIL: LoD RT source full-detail helper should submit current source query params after template cleanup\n");
	    ret = 1;
	}
    }

    BRLObolLodService submitService;
    if (!submitService.start(1, TRUE)) {
	printf("FAIL: LoD RT source full-detail helper service did not start\n");
	ret = 1;
    } else {
	uint64_t taskId = brlobol_lod_submit_rt_source_full_detail_request(
			      &submitService, submitService.beginGeneration(), sourceRequest,
			      dbip, &request, 10, 10);
	if (taskId == 0) {
	    printf("FAIL: LoD RT source full-detail helper did not submit source request\n");
	    ret = 1;
	} else if (wait_for_settled(submitService, 1)) {
	    ret = 1;
	} else {
	    std::vector<BRLObolLodResult> helperResults;
	    submitService.drainResults(helperResults);
	    if (helperResults.size() != 1 ||
		check_source_full_detail_result(helperResults[0],
						convertedRequest,
						"did not publish helper-submitted BoT mesh result"))
		ret = 1;
	}
	submitService.stop();
    }

    BRLObolRtSourceFullDetailProvider limitProvider;
    limitProvider.dbip = dbip;
    limitProvider.maxFullDetailFaceCount = 3;
    BRLObolLodResult limitedResult =
	brlobol_rt_source_full_detail_provider_task(request, &limitProvider);
    if (limitedResult.providerStatus != BRLOBOL_LOD_PROVIDER_FALLBACK ||
	strcmp(limitedResult.diagnostic.getString(),
	       "RT source full-detail provider request exceeds full-detail limits") != 0 ||
	limitedResult.mesh.isValid()) {
	printf("FAIL: LoD RT source full-detail provider did not refuse over-budget full detail\n");
	ret = 1;
    }

    BRLObolLodRequest missingRequest = request;
    missingRequest.objectPath = "/missing.bot";
    missingRequest.objectName = "missing.bot";
    BRLObolLodResult missingResult =
	brlobol_rt_source_full_detail_provider_task(missingRequest, &provider);
    if (missingResult.providerStatus != BRLOBOL_LOD_PROVIDER_ERROR ||
	strcmp(missingResult.diagnostic.getString(),
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

    bu_dir(cache_dir, MAXPATHLEN, BU_DIR_CURR, "brlobol_lod_service_cache", NULL);
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

    struct rt_mesh_lod_cache_status storeStatus =
	    RT_MESH_LOD_CACHE_STATUS_INIT;
    if (db_mesh_lod_store_mesh(dbip, "lod-provider.bot",
			       (const point_t *)cachedVertices.data(),
			       cachedVertices.size() / 3,
			       (const vect_t *)cachedNormals.data(), cachedFaces.data(),
			       cachedFaces.size() / 3, 0x12345678, 0.66,
			       &storeStatus) != BRLCAD_OK ||
	!storeStatus.has_cache_key ||
	!storeStatus.has_cached_payload) {
	printf("FAIL: LoD RT provider did not store cached mesh normals\n");
	db_mesh_lod_clear(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    BRLObolRtMeshLodProvider provider;
    provider.dbip = dbip;
    provider.useView = TRUE;
    provider.refreshMissing = TRUE;
    provider.view.size = 100.0;
    provider.view.width = 640;
    provider.view.height = 480;

    BRLObolLodTask task;
    task.generation = 1;
    task.request = make_request("/lod-provider.bot");
    task.request.objectName = "lod-provider.bot";
    task.request.providerId = "rt_mesh_lod";
    task.request.providerVersion = "rt-cache-v1";
    task.realize = brlobol_rt_mesh_lod_provider_task;
    task.realizeData = &provider;

    BRLObolLodService service;
    if (!service.start(1, TRUE)) {
	printf("FAIL: LoD RT provider service did not start\n");
	db_mesh_lod_clear(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    if (service.submit(task) == 0) {
	printf("FAIL: LoD RT provider service did not accept task\n");
	service.stop();
	db_mesh_lod_clear(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    if (wait_for_settled(service, 1)) {
	service.stop();
	db_mesh_lod_clear(dbip);
	db_close(dbip);
	bu_file_delete(dbpath);
	bu_dirclear(cache_dir);
	return 1;
    }

    std::vector<BRLObolLodResult> results;
    service.drainResults(results);
    service.stop();

    int ret = 0;
    if (results.size() != 1 ||
	results[0].resultKind != BRLOBOL_LOD_RESULT_MESH ||
	results[0].qualityTier != task.request.qualityTier ||
	results[0].providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	results[0].geometry.kind != BRLOBOL_LOD_GEOMETRY_RT_MESH_CACHE ||
	results[0].geometry.providerToken == 0 ||
	!results[0].geometry.isValid() ||
	results[0].counts.faceCount != 4 ||
	results[0].counts.pointCount != 4 ||
	!results[0].mesh.isValid() ||
	results[0].mesh.points.size() != 4 ||
	results[0].mesh.coordIndex.size() != 12 ||
	!brlobol_lod_result_matches_request(results[0], task.request)) {
	printf("FAIL: LoD RT provider task did not return cached mesh result\n");
	ret = 1;
    }

    BRLObolRtMeshLodProvider cachedNormalProvider;
    cachedNormalProvider.dbip = dbip;
    cachedNormalProvider.useForcedLevel = TRUE;
    cachedNormalProvider.forcedLevel = 1;
    cachedNormalProvider.refreshMissing = FALSE;
    BRLObolLodResult cachedNormalResult =
	brlobol_rt_mesh_lod_provider_task(task.request, &cachedNormalProvider);
    SbBool sawSeededNormal = FALSE;
    for (size_t i = 0; i < cachedNormalResult.mesh.normals.size(); i++) {
	if (cachedNormalResult.mesh.normals[i][2] > 1.5f)
	    sawSeededNormal = TRUE;
    }
    if (cachedNormalResult.resultKind != BRLOBOL_LOD_RESULT_MESH ||
	cachedNormalResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	cachedNormalResult.geometry.activeLevel != 1 ||
	cachedNormalResult.counts.normalCount == 0 ||
	cachedNormalResult.counts.normalCount !=
	cachedNormalResult.mesh.coordIndex.size() ||
	!cachedNormalResult.hasNormals ||
	!cachedNormalResult.mesh.isValid() ||
	cachedNormalResult.mesh.normals.size() !=
	cachedNormalResult.mesh.coordIndex.size() ||
	!sawSeededNormal ||
	!brlobol_lod_result_matches_request(cachedNormalResult,
					    task.request)) {
	printf("FAIL: LoD RT provider did not return cached mesh normals\n");
	ret = 1;
    }

    BRLObolRtMeshLodProvider forcedProvider;
    forcedProvider.dbip = dbip;
    forcedProvider.useForcedLevel = TRUE;
    forcedProvider.forcedLevel = (results.size() == 1 &&
				  results[0].geometry.activeLevel >= 0) ?
				 results[0].geometry.activeLevel : 0;
    forcedProvider.refreshMissing = FALSE;
    forcedProvider.shrinkAfterCopy = TRUE;

    BRLObolLodRequest forcedRequest = task.request;
    forcedRequest.qualityTier = BRLOBOL_LOD_QUALITY_FULL_DETAIL;
    BRLObolLodResult forcedResult =
	brlobol_rt_mesh_lod_provider_task(forcedRequest, &forcedProvider);
    if (forcedResult.resultKind != BRLOBOL_LOD_RESULT_MESH ||
	forcedResult.qualityTier != BRLOBOL_LOD_QUALITY_FULL_DETAIL ||
	forcedResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	forcedResult.geometry.activeLevel != forcedProvider.forcedLevel ||
	!forcedResult.mesh.isValid()) {
	printf("FAIL: LoD RT provider forced-level task did not return requested level mesh result\n");
	ret = 1;
    }

    struct rt_mesh_lod_cache_status invalidateStatus =
	    RT_MESH_LOD_CACHE_STATUS_INIT;
    if (db_mesh_lod_invalidate(dbip, "lod-provider.bot",
			       &invalidateStatus) != BRLCAD_OK ||
	!invalidateStatus.cleared_cache_entry ||
	!invalidateStatus.cleared_cache_key ||
	invalidateStatus.has_cache_key ||
	invalidateStatus.has_cached_payload) {
	printf("FAIL: LoD RT provider database invalidation status failed\n");
	ret = 1;
    }

    BRLObolRtMeshLodProvider staleProvider;
    staleProvider.dbip = dbip;
    staleProvider.refreshMissing = FALSE;

    BRLObolLodResult staleResult =
	brlobol_rt_mesh_lod_provider_task(task.request, &staleProvider);
    if (staleResult.providerStatus != BRLOBOL_LOD_PROVIDER_CACHE_MISS ||
	!staleResult.stale ||
	strcmp(staleResult.diagnostic.getString(),
	       "RT mesh LoD provider has no cache payload") != 0 ||
	staleResult.mesh.isValid()) {
	printf("FAIL: LoD RT provider did not report cache miss after database invalidation\n");
	ret = 1;
    }

    BRLObolRtMeshLodProvider refreshProvider;
    refreshProvider.dbip = dbip;
    refreshProvider.useView = TRUE;
    refreshProvider.refreshMissing = TRUE;
    refreshProvider.view = provider.view;

    BRLObolLodResult refreshResult =
	brlobol_rt_mesh_lod_provider_task(task.request, &refreshProvider);
    if (refreshResult.resultKind != BRLOBOL_LOD_RESULT_MESH ||
	refreshResult.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	refreshResult.geometry.kind != BRLOBOL_LOD_GEOMETRY_RT_MESH_CACHE ||
	!refreshResult.geometry.isValid() ||
	refreshResult.counts.faceCount != 4 ||
	refreshResult.counts.pointCount != 4 ||
	!refreshResult.mesh.isValid() ||
	!brlobol_lod_result_matches_request(refreshResult, task.request)) {
	printf("FAIL: LoD RT provider did not refresh cache after database invalidation\n");
	ret = 1;
    }

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

    if (test_dependency_order_and_cache_write())
	return 1;
    if (test_filtered_result_drain())
	return 1;
    if (test_generation_cancellation())
	return 1;
    if (test_stale_result_rejection())
	return 1;
    if (test_staged_payload_delivery())
	return 1;
    if (test_result_ready_subscription())
	return 1;
    if (test_task_realize_data_cleanup())
	return 1;
    if (test_debug_delay_cancellation())
	return 1;
    if (test_active_request_duplicate_suppression())
	return 1;
    if (test_rt_source_full_detail_provider_task())
	return 1;
    if (test_rt_mesh_provider_task())
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
