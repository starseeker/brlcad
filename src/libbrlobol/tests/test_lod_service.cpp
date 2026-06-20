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

    *dbip_out = dbip;
    return 0;
}

static int
test_rt_mesh_provider_task(void)
{
    char cache_dir[MAXPATHLEN] = {0};
    char dbpath[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;

    bu_dir(cache_dir, MAXPATHLEN, BU_DIR_CURR, "brlobol_lod_service_cache", NULL);
    bu_dirclear(cache_dir);
    bu_mkdir(cache_dir);
    bu_setenv("BU_DIR_CACHE", cache_dir, 1);

    if (make_provider_test_db(dbpath, sizeof(dbpath), &dbip)) {
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
