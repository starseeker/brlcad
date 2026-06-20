/*                L O D _ S E R V I C E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/lod_service.h"

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <mutex>
#include <set>
#include <string.h>
#include <thread>
#include <vector>

BRLObolRtMeshLodProvider::BRLObolRtMeshLodProvider(void)
{
    clear();
}

void
BRLObolRtMeshLodProvider::clear(void)
{
    dbip = NULL;
    rt_view_info_init(&view);
    useView = FALSE;
    refreshMissing = TRUE;
    useForcedLevel = FALSE;
    shrinkAfterCopy = TRUE;
    forcedLevel = 0;
    reset = 0;
}

BRLObolLodTask::BRLObolLodTask(void)
{
    clear();
}

void
BRLObolLodTask::clear(void)
{
    generation = 0;
    request.clear();
    dependencies.clear();
    realize = NULL;
    realizeData = NULL;
    realizeDataFree = NULL;
    cacheWrite = NULL;
    cacheWriteData = NULL;
    debugDelayMilliseconds = 0;
    publishResult = TRUE;
    writeCache = FALSE;
}

void
BRLObolLodTask::addDependency(uint64_t taskId)
{
    if (taskId != 0)
	dependencies.push_back(taskId);
}

static const char *
lod_request_object_name(const BRLObolLodRequest &request)
{
    const char *name = request.objectName.getString();
    if (name && name[0])
	return name;

    name = request.objectPath.getString();
    if (!name || !name[0])
	return NULL;

    const char *slash = strrchr(name, '/');
    if (slash && slash[1])
	return slash + 1;

    while (*name == '/')
	name++;
    return name[0] ? name : NULL;
}

static BRLObolLodResult
lod_provider_status_result(const BRLObolLodRequest &request, int status,
	const char *diagnostic)
{
    BRLObolLodResult result;

    result.request = request;
    result.cacheKey = brlobol_lod_cache_key(request);
    result.qualityTier = request.qualityTier;
    result.providerStatus = status;
    result.terminal = TRUE;
    result.diagnostic = diagnostic ? diagnostic : "";
    if (status == BRLOBOL_LOD_PROVIDER_CACHE_MISS ||
	    status == BRLOBOL_LOD_PROVIDER_STALE)
	result.stale = TRUE;

    return result;
}

BRLObolLodResult
brlobol_rt_mesh_lod_provider_task(const BRLObolLodRequest &request,
	void *userData)
{
    BRLObolRtMeshLodProvider *provider =
	static_cast<BRLObolRtMeshLodProvider *>(userData);
    if (!provider || !provider->dbip)
	return lod_provider_status_result(request, BRLOBOL_LOD_PROVIDER_ERROR,
	    "RT mesh LoD provider has no database");

    const char *name = lod_request_object_name(request);
    if (!name)
	return lod_provider_status_result(request, BRLOBOL_LOD_PROVIDER_ERROR,
	    "RT mesh LoD provider request has no object name");

    struct rt_mesh_lod_cache_status status = RT_MESH_LOD_CACHE_STATUS_INIT;
    if (db_mesh_lod_status(provider->dbip, name, &status) != BRLCAD_OK)
	return lod_provider_status_result(request, BRLOBOL_LOD_PROVIDER_ERROR,
	    "RT mesh LoD provider could not query cache status");

    if ((!status.has_cache_key || !status.has_cached_payload ||
		status.stale_cache_entry) && provider->refreshMissing) {
	if (db_mesh_lod_refresh(provider->dbip, name, &status) != BRLCAD_OK)
	    return lod_provider_status_result(request,
		BRLOBOL_LOD_PROVIDER_CACHE_MISS,
		"RT mesh LoD provider could not refresh cache entry");
    }

    struct rt_mesh_lod *lod = db_mesh_lod_get(provider->dbip, name);
    if (!lod)
	return lod_provider_status_result(request,
	    status.stale_cache_entry ? BRLOBOL_LOD_PROVIDER_STALE :
	    BRLOBOL_LOD_PROVIDER_CACHE_MISS,
	    "RT mesh LoD provider has no cache payload");

    int load_ret = provider->useForcedLevel ?
	rt_mesh_lod_load_level(lod, provider->forcedLevel, provider->reset) :
	(provider->useView ?
	    rt_mesh_lod_load_view(lod, &provider->view, provider->reset) :
	    rt_mesh_lod_load_view(lod, NULL, provider->reset));
    if (load_ret < 0) {
	rt_mesh_lod_destroy(lod);
	return lod_provider_status_result(request,
	    BRLOBOL_LOD_PROVIDER_CACHE_MISS,
	    "RT mesh LoD provider could not load a view level");
    }

    struct rt_mesh_lod_info info = RT_MESH_LOD_INFO_INIT;
    int have_info = rt_mesh_lod_info_get(lod, &info);
    if (!rt_mesh_lod_has_active_data(lod)) {
	BRLObolLodResult result =
	    brlobol_lod_result_from_rt_mesh_info(request, info, &status);
	rt_mesh_lod_destroy(lod);
	result.providerStatus = BRLOBOL_LOD_PROVIDER_CACHE_MISS;
	result.diagnostic = "RT mesh LoD provider loaded no active mesh data";
	return result;
    }
    if (!have_info) {
	rt_mesh_lod_destroy(lod);
	return lod_provider_status_result(request,
	    BRLOBOL_LOD_PROVIDER_CACHE_MISS,
	    "RT mesh LoD provider loaded no mesh metadata");
    }

    BRLObolLodResult result =
	brlobol_lod_result_from_rt_mesh_info(request, info, &status);
    if (result.providerStatus == BRLOBOL_LOD_PROVIDER_READY) {
	struct rt_mesh_lod_data data;
	if (!rt_mesh_lod_data_get(lod, &data) ||
		!brlobol_lod_mesh_payload_from_rt_mesh_data(result.mesh, data)) {
	    rt_mesh_lod_destroy(lod);
	    return lod_provider_status_result(request,
		BRLOBOL_LOD_PROVIDER_CACHE_MISS,
		"RT mesh LoD provider could not copy active mesh payload");
	}
	if (provider->shrinkAfterCopy)
	    rt_mesh_lod_memshrink(lod);
    }

    rt_mesh_lod_destroy(lod);
    return result;
}

void
brlobol_rt_mesh_lod_provider_free(void *userData)
{
    BRLObolRtMeshLodProvider *provider =
	static_cast<BRLObolRtMeshLodProvider *>(userData);
    delete provider;
}

struct BRLObolLodWorkItem {
    uint64_t id;
    BRLObolLodTask task;
};

struct BRLObolLodCacheWriteItem {
    BRLObolLodResult result;
    BRLObolLodCacheWriteProc write;
    void *writeData;
};

struct BRLObolLodSubscriber {
    BRLObolLodSubscriber(void) :
	id(0),
	callback(NULL),
	userData(NULL),
	active(FALSE),
	inFlight(0)
    {
    }

    BRLObolLodSubscriberId id;
    BRLObolLodResultReadyCB callback;
    void *userData;
    SbBool active;
    size_t inFlight;
};

struct BRLObolLodServicePrivate {
    explicit BRLObolLodServicePrivate(BRLObolLodService *newOwner) :
	owner(newOwner),
	running(FALSE),
	stopping(FALSE),
	cacheWriterStopping(FALSE),
	cacheWriterEnabled(FALSE),
	nextTaskId(1),
	nextSubscriberId(1),
	nextGeneration(0),
	activeGeneration(0),
	inFlight(0),
	cacheWriteInFlight(0),
	delayedTasks(0)
    {
    }

    BRLObolLodService *owner;
    mutable std::mutex mutex;
    std::condition_variable workerCv;
    std::condition_variable cacheWriterCv;
    std::condition_variable subscriberCv;
    std::vector<std::thread> workers;
    std::thread cacheWriter;
    SbBool running;
    SbBool stopping;
    SbBool cacheWriterStopping;
    SbBool cacheWriterEnabled;
    uint64_t nextTaskId;
    BRLObolLodSubscriberId nextSubscriberId;
    uint64_t nextGeneration;
    uint64_t activeGeneration;
    std::deque<BRLObolLodWorkItem> pending;
    std::deque<BRLObolLodResult> results;
    std::deque<BRLObolLodCacheWriteItem> cacheWrites;
    std::vector<BRLObolLodSubscriber> subscribers;
    std::set<uint64_t> completed;
    std::set<uint64_t> cancelledGenerations;
    size_t inFlight;
    size_t cacheWriteInFlight;
    size_t delayedTasks;
};

static SbBool
lod_generation_cancelled_unlocked(const BRLObolLodServicePrivate *p,
	uint64_t generation)
{
    return p->cancelledGenerations.find(generation) !=
	p->cancelledGenerations.end() ? TRUE : FALSE;
}

static SbBool
lod_generation_cancelled(BRLObolLodServicePrivate *p, uint64_t generation)
{
    std::lock_guard<std::mutex> lock(p->mutex);
    return lod_generation_cancelled_unlocked(p, generation);
}

static SbBool
lod_generation_cancelled_or_stopping(BRLObolLodServicePrivate *p,
	uint64_t generation)
{
    std::lock_guard<std::mutex> lock(p->mutex);
    return p->stopping ||
	lod_generation_cancelled_unlocked(p, generation) ? TRUE : FALSE;
}

static SbBool
lod_request_has_identity(const BRLObolLodRequest &request)
{
    return request.databaseId.getLength() > 0 ||
	request.objectPath.getLength() > 0 ||
	request.objectName.getLength() > 0 ||
	request.sourceRevision != 0 ||
	request.sourceContentHash != 0 ? TRUE : FALSE;
}

static BRLObolLodResult
lod_service_status_result(const BRLObolLodTask &task, int status,
	const char *diagnostic)
{
    BRLObolLodResult result;

    result.request = task.request;
    result.cacheKey = brlobol_lod_cache_key(task.request);
    result.qualityTier = task.request.qualityTier;
    result.providerStatus = status;
    result.terminal = TRUE;
    result.diagnostic = diagnostic ? diagnostic : "";
    if (status == BRLOBOL_LOD_PROVIDER_CANCELLED ||
	    status == BRLOBOL_LOD_PROVIDER_STALE)
	result.stale = TRUE;

    return result;
}

static void
lod_delayed_task_count_add(BRLObolLodServicePrivate *p, int delta)
{
    std::lock_guard<std::mutex> lock(p->mutex);

    if (delta > 0) {
	p->delayedTasks += (size_t)delta;
	return;
    }

    size_t decrement = (size_t)(-delta);
    p->delayedTasks = decrement > p->delayedTasks ?
	0 : p->delayedTasks - decrement;
}

static SbBool
lod_wait_for_debug_delay(BRLObolLodServicePrivate *p,
	const BRLObolLodTask &task)
{
    if (task.debugDelayMilliseconds == 0)
	return TRUE;

    lod_delayed_task_count_add(p, 1);

    uint32_t remaining = task.debugDelayMilliseconds;
    while (remaining > 0) {
	if (lod_generation_cancelled_or_stopping(p, task.generation)) {
	    lod_delayed_task_count_add(p, -1);
	    return FALSE;
	}

	uint32_t slice = remaining > 10 ? 10 : remaining;
	std::this_thread::sleep_for(std::chrono::milliseconds(slice));
	remaining -= slice;
    }

    lod_delayed_task_count_add(p, -1);
    return !lod_generation_cancelled_or_stopping(p, task.generation);
}

static void
lod_normalize_result(BRLObolLodResult &result, const BRLObolLodTask &task)
{
    if (!result.cacheKey.isValid()) {
	if (lod_request_has_identity(result.request)) {
	    result.cacheKey = brlobol_lod_cache_key(result.request);
	} else {
	    result.request = task.request;
	    result.cacheKey = brlobol_lod_cache_key(task.request);
	}
    }

    if (!lod_request_has_identity(result.request))
	result.request = task.request;

    if (!brlobol_lod_result_matches_request(result, task.request)) {
	result.providerStatus = BRLOBOL_LOD_PROVIDER_STALE;
	result.stale = TRUE;
	result.terminal = TRUE;
	if (result.diagnostic.getLength() == 0)
	    result.diagnostic = "LoD task returned stale request/result";
    }
}

static SbBool
lod_task_dependencies_ready(const BRLObolLodServicePrivate *p,
	const BRLObolLodTask &task)
{
    if (lod_generation_cancelled_unlocked(p, task.generation))
	return TRUE;

    for (size_t i = 0; i < task.dependencies.size(); i++) {
	if (p->completed.find(task.dependencies[i]) == p->completed.end())
	    return FALSE;
    }

    return TRUE;
}

static std::deque<BRLObolLodWorkItem>::iterator
lod_find_ready_task(BRLObolLodServicePrivate *p)
{
    std::deque<BRLObolLodWorkItem>::iterator it;

    for (it = p->pending.begin(); it != p->pending.end(); ++it) {
	if (lod_task_dependencies_ready(p, it->task))
	    return it;
    }

    return p->pending.end();
}

static BRLObolLodResult
lod_execute_task(BRLObolLodServicePrivate *p, const BRLObolLodTask &task)
{
    if (lod_generation_cancelled(p, task.generation))
	return lod_service_status_result(task, BRLOBOL_LOD_PROVIDER_CANCELLED,
	    "LoD task generation cancelled");

    if (!lod_wait_for_debug_delay(p, task))
	return lod_service_status_result(task, BRLOBOL_LOD_PROVIDER_CANCELLED,
	    "LoD task generation cancelled during debug delay");

    if (!task.realize)
	return lod_service_status_result(task, BRLOBOL_LOD_PROVIDER_ERROR,
	    "LoD task has no realization callback");

    BRLObolLodResult result = (*task.realize)(task.request, task.realizeData);

    if (lod_generation_cancelled(p, task.generation))
	return lod_service_status_result(task, BRLOBOL_LOD_PROVIDER_CANCELLED,
	    "LoD task generation cancelled");

    lod_normalize_result(result, task);
    return result;
}

static void
lod_task_free_realize_data(BRLObolLodTask &task)
{
    if (task.realizeDataFree && task.realizeData) {
	(*task.realizeDataFree)(task.realizeData);
	task.realizeData = NULL;
	task.realizeDataFree = NULL;
    }
}

struct BRLObolLodSubscriberCall {
    BRLObolLodSubscriberId id;
    BRLObolLodResultReadyCB callback;
    void *userData;
};

static std::vector<BRLObolLodSubscriberCall>
lod_collect_result_ready_callbacks(BRLObolLodServicePrivate *p)
{
    std::vector<BRLObolLodSubscriberCall> calls;
    std::lock_guard<std::mutex> lock(p->mutex);

    for (size_t i = 0; i < p->subscribers.size(); i++) {
	BRLObolLodSubscriber &subscriber = p->subscribers[i];
	if (!subscriber.active || !subscriber.callback)
	    continue;

	BRLObolLodSubscriberCall call;
	call.id = subscriber.id;
	call.callback = subscriber.callback;
	call.userData = subscriber.userData;
	subscriber.inFlight++;
	calls.push_back(call);
    }

    return calls;
}

static void
lod_complete_result_ready_callback(BRLObolLodServicePrivate *p,
	BRLObolLodSubscriberId id)
{
    std::lock_guard<std::mutex> lock(p->mutex);

    for (size_t i = 0; i < p->subscribers.size(); i++) {
	if (p->subscribers[i].id != id)
	    continue;
	if (p->subscribers[i].inFlight > 0)
	    p->subscribers[i].inFlight--;
	break;
    }

    p->subscriberCv.notify_all();
}

static void
lod_notify_result_ready(BRLObolLodServicePrivate *p)
{
    std::vector<BRLObolLodSubscriberCall> calls =
	lod_collect_result_ready_callbacks(p);

    for (size_t i = 0; i < calls.size(); i++) {
	if (calls[i].callback)
	    (*calls[i].callback)(p->owner, calls[i].userData);
	lod_complete_result_ready_callback(p, calls[i].id);
    }
}

static void
lod_finish_task(BRLObolLodServicePrivate *p, const BRLObolLodWorkItem &item,
	const BRLObolLodResult &result)
{
    SbBool notifyResultReady = FALSE;

    {
	std::lock_guard<std::mutex> lock(p->mutex);

	p->completed.insert(item.id);
	if (p->inFlight > 0)
	    p->inFlight--;

	if (item.task.publishResult) {
	    notifyResultReady = p->results.empty() ? TRUE : FALSE;
	    p->results.push_back(result);
	}

	if (p->cacheWriterEnabled && item.task.writeCache &&
		item.task.cacheWrite) {
	    BRLObolLodCacheWriteItem writeItem;
	    writeItem.result = result;
	    writeItem.write = item.task.cacheWrite;
	    writeItem.writeData = item.task.cacheWriteData;
	    p->cacheWrites.push_back(writeItem);
	}
    }

    p->workerCv.notify_all();
    p->cacheWriterCv.notify_one();
    if (notifyResultReady)
	lod_notify_result_ready(p);
}

static void
lod_worker_loop(BRLObolLodServicePrivate *p)
{
    for (;;) {
	BRLObolLodWorkItem item;

	{
	    std::unique_lock<std::mutex> lock(p->mutex);

	    while (!p->stopping && lod_find_ready_task(p) == p->pending.end())
		p->workerCv.wait(lock);

	    if (p->stopping)
		return;

	    std::deque<BRLObolLodWorkItem>::iterator ready =
		lod_find_ready_task(p);
	    if (ready == p->pending.end())
		continue;
	    item = *ready;
	    p->pending.erase(ready);
	}

	BRLObolLodResult result = lod_execute_task(p, item.task);
	lod_finish_task(p, item, result);
	lod_task_free_realize_data(item.task);
    }
}

static void
lod_cache_writer_loop(BRLObolLodServicePrivate *p)
{
    for (;;) {
	BRLObolLodCacheWriteItem item;

	{
	    std::unique_lock<std::mutex> lock(p->mutex);

	    while (!p->cacheWriterStopping && p->cacheWrites.empty())
		p->cacheWriterCv.wait(lock);

	    if (p->cacheWrites.empty() && p->cacheWriterStopping)
		return;

	    if (p->cacheWrites.empty())
		continue;

	    item = p->cacheWrites.front();
	    p->cacheWrites.pop_front();
	    p->cacheWriteInFlight++;
	}

	if (item.write)
	    (*item.write)(item.result, item.writeData);

	{
	    std::lock_guard<std::mutex> lock(p->mutex);
	    if (p->cacheWriteInFlight > 0)
		p->cacheWriteInFlight--;
	}
	p->cacheWriterCv.notify_all();
    }
}

BRLObolLodService::BRLObolLodService(void) :
    p(new BRLObolLodServicePrivate(this))
{
}

BRLObolLodService::~BRLObolLodService(void)
{
    this->stop();
    delete this->p;
    this->p = NULL;
}

SbBool
BRLObolLodService::start(size_t workerCount, SbBool startCacheWriter)
{
    if (workerCount == 0)
	workerCount = 1;

    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	if (this->p->running)
	    return TRUE;

	this->p->stopping = FALSE;
	this->p->cacheWriterStopping = FALSE;
	this->p->cacheWriterEnabled = startCacheWriter ? TRUE : FALSE;
	this->p->running = TRUE;
    }

    try {
	if (startCacheWriter)
	    this->p->cacheWriter =
		std::thread(lod_cache_writer_loop, this->p);

	for (size_t i = 0; i < workerCount; i++)
	    this->p->workers.push_back(std::thread(lod_worker_loop, this->p));
    } catch (...) {
	this->stop();
	return FALSE;
    }

    return TRUE;
}

void
BRLObolLodService::stop(void)
{
    std::vector<std::thread> workers;
    std::thread cacheWriter;
    std::deque<BRLObolLodWorkItem> pending;

    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	if (!this->p->running && this->p->workers.empty() &&
		!this->p->cacheWriter.joinable())
	    return;
	this->p->stopping = TRUE;
    }
    this->p->workerCv.notify_all();

    workers.swap(this->p->workers);
    for (size_t i = 0; i < workers.size(); i++) {
	if (workers[i].joinable())
	    workers[i].join();
    }

    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	pending.swap(this->p->pending);
	this->p->inFlight = 0;
	this->p->cacheWriterStopping = TRUE;
    }
    for (size_t i = 0; i < pending.size(); i++)
	lod_task_free_realize_data(pending[i].task);
    this->p->cacheWriterCv.notify_all();

    if (this->p->cacheWriter.joinable())
	cacheWriter.swap(this->p->cacheWriter);
    if (cacheWriter.joinable())
	cacheWriter.join();

    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	this->p->cacheWrites.clear();
	this->p->cacheWriteInFlight = 0;
	this->p->delayedTasks = 0;
	this->p->running = FALSE;
	this->p->stopping = FALSE;
	this->p->cacheWriterStopping = FALSE;
	this->p->cacheWriterEnabled = FALSE;
    }
}

SbBool
BRLObolLodService::isRunning(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->running;
}

uint64_t
BRLObolLodService::beginGeneration(void)
{
    std::lock_guard<std::mutex> lock(this->p->mutex);

    this->p->nextGeneration++;
    if (this->p->nextGeneration == 0)
	this->p->nextGeneration++;
    this->p->activeGeneration = this->p->nextGeneration;
    return this->p->activeGeneration;
}

uint64_t
BRLObolLodService::currentGeneration(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->activeGeneration;
}

void
BRLObolLodService::cancelGeneration(uint64_t generation)
{
    if (generation == 0)
	return;

    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	this->p->cancelledGenerations.insert(generation);
    }
    this->p->workerCv.notify_all();
}

SbBool
BRLObolLodService::isGenerationCancelled(uint64_t generation) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return lod_generation_cancelled_unlocked(this->p, generation);
}

uint64_t
BRLObolLodService::submit(const BRLObolLodTask &task)
{
    uint64_t id = 0;

    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	if (!this->p->running || this->p->stopping)
	    return 0;

	BRLObolLodWorkItem item;
	item.id = this->p->nextTaskId++;
	if (this->p->nextTaskId == 0)
	    this->p->nextTaskId++;
	item.task = task;
	if (item.task.generation == 0) {
	    if (this->p->activeGeneration == 0) {
		this->p->nextGeneration++;
		this->p->activeGeneration = this->p->nextGeneration;
	    }
	    item.task.generation = this->p->activeGeneration;
	}

	id = item.id;
	this->p->pending.push_back(item);
	this->p->inFlight++;
    }

    this->p->workerCv.notify_all();
    return id;
}

size_t
BRLObolLodService::drainResults(std::vector<BRLObolLodResult> &results,
	size_t maxResults)
{
    size_t count = 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);

    while (!this->p->results.empty() &&
	    (maxResults == 0 || count < maxResults)) {
	results.push_back(this->p->results.front());
	this->p->results.pop_front();
	count++;
    }

    return count;
}

BRLObolLodSubscriberId
BRLObolLodService::subscribeResultReady(BRLObolLodResultReadyCB callback,
	void *userData)
{
    if (!callback)
	return 0;

    std::lock_guard<std::mutex> lock(this->p->mutex);

    BRLObolLodSubscriber subscriber;
    subscriber.id = this->p->nextSubscriberId++;
    if (this->p->nextSubscriberId == 0)
	this->p->nextSubscriberId++;
    subscriber.callback = callback;
    subscriber.userData = userData;
    subscriber.active = TRUE;
    this->p->subscribers.push_back(subscriber);

    return subscriber.id;
}

void
BRLObolLodService::unsubscribeResultReady(BRLObolLodSubscriberId id)
{
    if (id == 0)
	return;

    std::unique_lock<std::mutex> lock(this->p->mutex);

    for (size_t i = 0; i < this->p->subscribers.size(); i++) {
	if (this->p->subscribers[i].id != id)
	    continue;

	this->p->subscribers[i].active = FALSE;
	this->p->subscriberCv.wait(lock, [this, id] {
	    for (size_t j = 0; j < this->p->subscribers.size(); j++) {
		if (this->p->subscribers[j].id == id)
		    return this->p->subscribers[j].inFlight == 0;
	    }
	    return true;
	});
	for (size_t j = 0; j < this->p->subscribers.size(); j++) {
	    if (this->p->subscribers[j].id == id) {
		this->p->subscribers.erase(
			this->p->subscribers.begin() + (long)j);
		break;
	    }
	}
	return;
    }
}

size_t
BRLObolLodService::inFlightCount(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->inFlight;
}

size_t
BRLObolLodService::pendingTaskCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->pending.size();
}

size_t
BRLObolLodService::queuedResultCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->results.size();
}

size_t
BRLObolLodService::queuedCacheWriteCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->cacheWrites.size() + this->p->cacheWriteInFlight;
}

size_t
BRLObolLodService::delayedTaskCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->delayedTasks;
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
