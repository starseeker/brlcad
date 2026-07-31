/*          S O U R C E _ R E A L I Z A T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file source_realization.cpp */

#include "common.h"

#include "BObol/BSourceRealization.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BLodService.h"
#include "bu/app.h"
#include "bu/file.h"
#include "bu/parallel.h"
#include "rt/db_io.h"

#include <Inventor/SbString.h>

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdlib>
#include <deque>
#include <mutex>
#include <thread>

struct SourceRealizationItem {
    SourceRealizationItem(void) :
	source(NULL), snapshotSourceDatabase(NULL), database(NULL),
	clientToken(0), estimatedWorkingSetBytes(0), drawMode(0),
	allowWireFallback(FALSE),
	probeWarmManifest(NULL), storeManifest(NULL), callbackData(NULL),
	state(BOBOL_SOURCE_REALIZATION_PENDING), warmManifest(false),
	manifestStored(false)
    {
    }

    ~SourceRealizationItem(void)
    {
	/*
	 * Detached sources borrow database.  Destroy the source first, matching
	 * the former libged ownership order, then close the independent handles.
	 */
	if (source)
	    source->unref();
	if (database)
	    db_close(database);
	if (snapshotSourceDatabase)
	    db_close(snapshotSourceDatabase);
	if (snapshotPath.getLength() > 0)
	    (void)bu_file_delete(snapshotPath.getString());
    }

    SoBRLDatabaseSource *source;
    struct db_i *snapshotSourceDatabase;
    struct db_i *database;
    SbString snapshotPath;
    std::shared_ptr<BObolCompactOccurrenceStream> stream;
    uint64_t clientToken;
    size_t estimatedWorkingSetBytes;
    int drawMode;
    SbBool allowWireFallback;
    BObolSourceWarmManifestProbe probeWarmManifest;
    BObolSourceManifestStore storeManifest;
    void *callbackData;
    std::atomic<int> state;
    std::atomic<bool> warmManifest;
    std::atomic<bool> manifestStored;
};

static bool
source_job_cancelled(const BObolSourceRealizationJobPrivate *job);

struct BObolSourceRealizationJobPrivate {
    BObolSourceRealizationJobPrivate(void) :
	state(BOBOL_SOURCE_REALIZATION_PENDING), cancelRequested(false),
	remaining(0), failed(false)
    {
    }

    std::atomic<int> state;
    std::atomic<bool> cancelRequested;
    std::atomic<size_t> remaining;
    std::atomic<bool> failed;
    std::vector<std::unique_ptr<SourceRealizationItem>> items;
};

struct SourceRealizationWork {
    std::shared_ptr<BObolSourceRealizationJobPrivate> job;
    size_t itemIndex = 0;
};

static size_t
source_work_estimated_bytes(const SourceRealizationWork &work)
{
    if (!work.job || work.itemIndex >= work.job->items.size() ||
	!work.job->items[work.itemIndex])
	return 1;
    const size_t estimate =
	work.job->items[work.itemIndex]->estimatedWorkingSetBytes;
    /*
     * An unknown root may still open a multi-gigabyte directory and build
     * hierarchy state before mesh-level accounting begins.
     */
    return estimate ? estimate : 256ULL * 1024ULL * 1024ULL;
}

static bool
source_job_cancelled(const BObolSourceRealizationJobPrivate *job)
{
    return !job || job->cancelRequested.load(std::memory_order_acquire);
}

static void
source_item_cancel(SourceRealizationItem *item)
{
    if (!item)
	return;
    if (item->stream)
	item->stream->requestCancel();
    item->state.store(BOBOL_SOURCE_REALIZATION_CANCELLED,
	std::memory_order_release);
}

static void
source_realize_item(const std::shared_ptr<BObolSourceRealizationJobPrivate> &job,
	size_t itemIndex)
{
    if (!job || itemIndex >= job->items.size())
	return;
    SourceRealizationItem *item = job->items[itemIndex].get();
    if (!item)
	return;
    if (source_job_cancelled(job.get())) {
	source_item_cancel(item);
	return;
    }

    item->state.store(BOBOL_SOURCE_REALIZATION_RUNNING,
	std::memory_order_release);
    SoBRLDatabaseSource *source = item->source;
    bool success = source != NULL;

    int warmManifest = 0;
    if (success && item->probeWarmManifest && item->stream &&
	item->snapshotSourceDatabase && !source_job_cancelled(job.get())) {
	warmManifest = item->probeWarmManifest(source,
	    item->snapshotSourceDatabase, item->drawMode,
	    source->sourceRevision.getValue(), item->stream.get(),
	    item->callbackData);
    }
    item->warmManifest.store(warmManifest == 2, std::memory_order_release);
    if (item->warmManifest.load(std::memory_order_acquire) &&
	item->snapshotSourceDatabase) {
	db_close(item->snapshotSourceDatabase);
	item->snapshotSourceDatabase = NULL;
    }

    if (success && !item->warmManifest && !item->database) {
	if (source_job_cancelled(job.get())) {
	    success = false;
	} else if (!item->snapshotSourceDatabase ||
	    !source->initializeDetachedRealizationDatabase(
		item->snapshotSourceDatabase, &item->database,
		&item->snapshotPath)) {
	    success = false;
	}
	if (item->snapshotSourceDatabase) {
	    db_close(item->snapshotSourceDatabase);
	    item->snapshotSourceDatabase = NULL;
	}
    }

    if (success && !source_job_cancelled(job.get())) {
	const bool mesh = source->usesMeshRealization() ? true : false;
	BObolCompactOccurrenceStream *stream =
	    item->stream ? item->stream.get() : NULL;
	SbBool realized = warmManifest == 2 ? TRUE :
	    (mesh ? source->realizeDatabaseMesh(stream) :
		source->realizeDatabaseWireframe(stream));
	if (!realized && mesh && item->allowWireFallback &&
	    !source_job_cancelled(job.get()))
	    realized = source->realizeDatabaseWireframe(stream);
	success = realized ? true : false;
	if (success && stream && !stream->hasCoverageBoundsComplete())
	    stream->setCoverageBoundsComplete(true);
	if (success &&
	    !item->warmManifest.load(std::memory_order_acquire) &&
	    item->storeManifest &&
	    !source_job_cancelled(job.get()))
	    item->manifestStored.store(
		item->storeManifest(item->database, source,
		    item->callbackData) ? true : false,
		std::memory_order_release);
    } else {
	success = false;
    }

    if (source_job_cancelled(job.get())) {
	source_item_cancel(item);
	return;
    }
    if (!success) {
	job->failed.store(true, std::memory_order_release);
	job->cancelRequested.store(true, std::memory_order_release);
	for (const std::unique_ptr<SourceRealizationItem> &sibling :
	     job->items) {
	    if (sibling && sibling->stream)
		sibling->stream->requestCancel();
	}
	item->state.store(BOBOL_SOURCE_REALIZATION_FAILED,
	    std::memory_order_release);
	return;
    }
    item->state.store(BOBOL_SOURCE_REALIZATION_COMPLETE,
	std::memory_order_release);
}

static size_t
source_realization_worker_count(void)
{
    size_t count = bu_avail_cpus();
    if (count < 1)
	count = 1;
    /*
     * Each task may open a large database directory and run internally
     * parallel PoP preparation.  A modest outer pool exposes independent-root
     * parallelism without multiplying transient memory by the CPU count.
     */
    count = std::min<size_t>(count, 8);
    const char *setting = std::getenv("BOBOL_SOURCE_REALIZATION_WORKERS");
    if (setting && setting[0]) {
	char *end = NULL;
	const unsigned long parsed = std::strtoul(setting, &end, 10);
	if (end && !end[0] && parsed > 0)
	    count = std::min<size_t>(parsed, 64);
    }
    return count;
}

struct BObolSourceRealizationCoordinatorPrivate {
    BObolSourceRealizationCoordinatorPrivate(void) :
	stopping(false), active(0), activeBytes(0),
	maxActiveBytes(bobol_lod_working_set_global_limit())
    {
    }

    std::mutex mutex;
    std::condition_variable cv;
    bool stopping;
    size_t active;
    size_t activeBytes;
    size_t maxActiveBytes;
    std::deque<SourceRealizationWork> queue;
    std::vector<std::thread> workers;
};

static void
source_realization_worker(BObolSourceRealizationCoordinatorPrivate *service)
{
    for (;;) {
	SourceRealizationWork work;
	size_t admittedBytes = 0;
	{
	    std::unique_lock<std::mutex> lock(service->mutex);
	    service->cv.wait(lock, [service]() {
		if (service->stopping)
		    return true;
		for (const SourceRealizationWork &candidate : service->queue) {
		    const size_t estimate =
			source_work_estimated_bytes(candidate);
		    if (!service->active ||
			service->maxActiveBytes == SIZE_MAX ||
			(estimate <= service->maxActiveBytes &&
			 service->activeBytes <=
			    service->maxActiveBytes - estimate))
			return true;
		}
		return false;
	    });
	    if (service->stopping && service->queue.empty())
		return;
	    auto selected = service->queue.end();
	    for (auto candidate = service->queue.begin();
		    candidate != service->queue.end(); ++candidate) {
		const size_t estimate =
		    source_work_estimated_bytes(*candidate);
		if (!service->active ||
		    service->maxActiveBytes == SIZE_MAX ||
		    (estimate <= service->maxActiveBytes &&
		     service->activeBytes <=
			service->maxActiveBytes - estimate)) {
		    selected = candidate;
		    admittedBytes = estimate;
		    break;
		}
	    }
	    if (selected == service->queue.end()) {
		if (service->stopping)
		    return;
		continue;
	    }
	    work = std::move(*selected);
	    service->queue.erase(selected);
	    service->active++;
	    service->activeBytes =
		admittedBytes > SIZE_MAX - service->activeBytes ?
		SIZE_MAX : service->activeBytes + admittedBytes;
	}

	source_realize_item(work.job, work.itemIndex);

	const size_t remaining = work.job ?
	    work.job->remaining.fetch_sub(1, std::memory_order_acq_rel) : 0;
	if (work.job && remaining == 1) {
	    const int terminal = work.job->failed.load(std::memory_order_acquire) ?
		BOBOL_SOURCE_REALIZATION_FAILED :
		(source_job_cancelled(work.job.get()) ?
		    BOBOL_SOURCE_REALIZATION_CANCELLED :
		    BOBOL_SOURCE_REALIZATION_COMPLETE);
	    /*
	     * A failed item cancels siblings, but failure remains the aggregate
	     * result.  Release-store publishes every item/source write to the
	     * owner-thread result pump.
	     */
	    work.job->state.store(terminal, std::memory_order_release);
	}

	{
	    std::lock_guard<std::mutex> lock(service->mutex);
	    service->active--;
	    service->activeBytes =
		admittedBytes >= service->activeBytes ?
		0 : service->activeBytes - admittedBytes;
	}
	service->cv.notify_all();
    }
}

BObolSourceRealizationRequest::BObolSourceRealizationRequest(void) :
    source(NULL), snapshotSourceDatabase(NULL), clientToken(0),
    estimatedWorkingSetBytes(0), drawMode(0),
    allowWireFallback(FALSE), probeWarmManifest(NULL), storeManifest(NULL),
    callbackData(NULL)
{
}

BObolSourceRealizationItemResult::BObolSourceRealizationItemResult(void) :
    source(NULL), clientToken(0), state(BOBOL_SOURCE_REALIZATION_PENDING),
    warmManifest(FALSE), manifestStored(FALSE)
{
}

BObolSourceRealizationJob::BObolSourceRealizationJob(
    const std::shared_ptr<BObolSourceRealizationJobPrivate> &state) :
    p(state)
{
}

BObolSourceRealizationJob::~BObolSourceRealizationJob(void)
{
    this->cancel();
}

void
BObolSourceRealizationJob::cancel(void)
{
    if (!this->p)
	return;
    this->p->cancelRequested.store(true, std::memory_order_release);
    for (const std::unique_ptr<SourceRealizationItem> &item : this->p->items) {
	if (item && item->stream)
	    item->stream->requestCancel();
    }
}

int
BObolSourceRealizationJob::state(void) const
{
    return this->p ? this->p->state.load(std::memory_order_acquire) :
	BOBOL_SOURCE_REALIZATION_CANCELLED;
}

SbBool
BObolSourceRealizationJob::isTerminal(void) const
{
    const int current = this->state();
    return current == BOBOL_SOURCE_REALIZATION_COMPLETE ||
	current == BOBOL_SOURCE_REALIZATION_FAILED ||
	current == BOBOL_SOURCE_REALIZATION_CANCELLED ? TRUE : FALSE;
}

size_t
BObolSourceRealizationJob::itemCount(void) const
{
    return this->p ? this->p->items.size() : 0;
}

SbBool
BObolSourceRealizationJob::itemResult(
    size_t index, BObolSourceRealizationItemResult &result) const
{
    result = BObolSourceRealizationItemResult();
    if (!this->p || index >= this->p->items.size() ||
	!this->p->items[index])
	return FALSE;
    const SourceRealizationItem *item = this->p->items[index].get();
    result.source = item->source;
    result.stream = item->stream;
    result.clientToken = item->clientToken;
    result.state = item->state.load(std::memory_order_acquire);
    result.warmManifest = item->warmManifest.load(std::memory_order_acquire) ?
	TRUE : FALSE;
    result.manifestStored =
	item->manifestStored.load(std::memory_order_acquire) ? TRUE : FALSE;
    return TRUE;
}

BObolSourceRealizationCoordinator &
BObolSourceRealizationCoordinator::global(void)
{
    static BObolSourceRealizationCoordinator coordinator;
    return coordinator;
}

BObolSourceRealizationCoordinator::BObolSourceRealizationCoordinator(void) :
    p(new BObolSourceRealizationCoordinatorPrivate)
{
    const size_t count = source_realization_worker_count();
    this->p->workers.reserve(count);
    for (size_t i = 0; i < count; i++)
	this->p->workers.push_back(
	    std::thread(source_realization_worker, this->p));
}

BObolSourceRealizationCoordinator::~BObolSourceRealizationCoordinator(void)
{
    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	this->p->stopping = true;
	for (SourceRealizationWork &work : this->p->queue) {
	    if (work.job)
		work.job->cancelRequested.store(true, std::memory_order_release);
	}
    }
    this->p->cv.notify_all();
    for (std::thread &worker : this->p->workers) {
	if (worker.joinable())
	    worker.join();
    }
    delete this->p;
}

std::shared_ptr<BObolSourceRealizationJob>
BObolSourceRealizationCoordinator::submit(
    std::vector<BObolSourceRealizationRequest> &requests)
{
    if (requests.empty())
	return std::shared_ptr<BObolSourceRealizationJob>();

    /*
     * Validate the complete batch before transferring any ownership.  A
     * malformed later item must not leave the caller with a half-consumed
     * request vector.
     */
    for (const BObolSourceRealizationRequest &request : requests) {
	if (!request.source || !request.snapshotSourceDatabase)
	    return std::shared_ptr<BObolSourceRealizationJob>();
    }

    std::shared_ptr<BObolSourceRealizationJobPrivate> job(
	new BObolSourceRealizationJobPrivate);
    job->items.reserve(requests.size());
    for (BObolSourceRealizationRequest &request : requests) {
	std::unique_ptr<SourceRealizationItem> item(
	    new SourceRealizationItem);
	item->source = request.source;
	item->snapshotSourceDatabase = request.snapshotSourceDatabase;
	item->stream = request.stream;
	item->clientToken = request.clientToken;
	item->estimatedWorkingSetBytes = request.estimatedWorkingSetBytes;
	item->drawMode = request.drawMode;
	item->allowWireFallback = request.allowWireFallback;
	item->probeWarmManifest = request.probeWarmManifest;
	item->storeManifest = request.storeManifest;
	item->callbackData = request.callbackData;
	request.source = NULL;
	request.snapshotSourceDatabase = NULL;
	job->items.push_back(std::move(item));
    }
    job->remaining.store(job->items.size(), std::memory_order_release);
    job->state.store(BOBOL_SOURCE_REALIZATION_RUNNING,
	std::memory_order_release);

    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	if (this->p->stopping)
	    return std::shared_ptr<BObolSourceRealizationJob>();
	for (size_t i = 0; i < job->items.size(); i++) {
	    SourceRealizationWork work;
	    work.job = job;
	    work.itemIndex = i;
	    this->p->queue.push_back(std::move(work));
	}
    }
    this->p->cv.notify_all();
    return std::shared_ptr<BObolSourceRealizationJob>(
	new BObolSourceRealizationJob(job));
}

size_t
BObolSourceRealizationCoordinator::workerCountForDiagnostics(void) const
{
    return this->p ? this->p->workers.size() : 0;
}

size_t
BObolSourceRealizationCoordinator::queuedItemCountForDiagnostics(void) const
{
    if (!this->p)
	return 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->queue.size();
}

size_t
BObolSourceRealizationCoordinator::activeItemCountForDiagnostics(void) const
{
    if (!this->p)
	return 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->active;
}

size_t
BObolSourceRealizationCoordinator::activeWorkingSetBytesForDiagnostics(
    void) const
{
    if (!this->p)
	return 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->activeBytes;
}

size_t
BObolSourceRealizationCoordinator::workingSetLimitBytesForDiagnostics(
    void) const
{
    if (!this->p)
	return 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->maxActiveBytes;
}
