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
#include "database_source_realization.h"
#include "draw_cache_private.h"
#include "parallel_budget_private.h"
#include "transaction_fault_private.h"
#include "bu/app.h"
#include "bu/file.h"
#include "bu/parallel.h"
#include "rt/db4.h"
#include "rt/db_io.h"

#include <Inventor/SbString.h>

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdlib>
#include <list>
#include <mutex>
#include <new>
#include <system_error>
#include <thread>

/* Internal cache lifetime barriers.  Construct the lazy cache registries
 * before the process-wide worker coordinator so reverse static destruction
 * joins every worker before either registry is released. */
void bobol_mesh_lod_cache_runtime_prepare(void);

struct BObolSourceRealizationItemPrivate {
    BObolSourceRealizationItemPrivate(void) :
	source(NULL), snapshotSourceDatabase(NULL), database(NULL),
	clientToken(0), estimatedWorkingSetBytes(0), drawMode(0),
	allowWireFallback(FALSE),
	probeWarmManifest(NULL), storeManifest(NULL),
	state(BOBOL_SOURCE_REALIZATION_PENDING), warmManifest(false),
	manifestStored(false)
    {
    }

    ~BObolSourceRealizationItemPrivate(void)
    {
	releaseRealization();
    }

    void releaseRealization(void)
    {
	/*
	 * Detached sources borrow database.  Destroy the source first, matching
	 * the former libged ownership order, then close the independent handles.
	 */
	if (source) {
	    source->unref();
	    source = NULL;
	}
	if (database) {
	    db_close(database);
	    database = NULL;
	}
	if (snapshotSourceDatabase) {
	    db_close(snapshotSourceDatabase);
	    snapshotSourceDatabase = NULL;
	}
	if (snapshotPath.getLength() > 0) {
	    (void)bu_file_delete(snapshotPath.getString());
	    snapshotPath.makeEmpty();
	}
    }

    void finish(int outcome)
    {
	/* Only the worker (or an unstarted item's queue owner) retires these
	 * resources. Cancellation requests cannot release an active callback's
	 * storage. Completed source/database pairs remain borrowed job results. */
	if (outcome != BOBOL_SOURCE_REALIZATION_COMPLETE) {
	    if (stream)
		(void)stream->cancelPublication();
	    releaseRealization();
	}
	callbackContext.reset();
	state.store(outcome, std::memory_order_release);
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
    std::shared_ptr<void> callbackContext;
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
    std::weak_ptr<BObolSourceRealizationCoordinatorPrivate> coordinator;
    std::vector<std::unique_ptr<BObolSourceRealizationItemPrivate>> items;
};

struct SourceRealizationWork {
    std::shared_ptr<BObolSourceRealizationJobPrivate> job;
    size_t itemIndex = 0;
    /* First-fit memory admission may bypass an older, temporarily too-large
     * root to keep otherwise idle workers useful.  Bound those bypasses so a
     * sustained stream of small roots cannot starve a vehicle-scale source. */
    size_t memoryBypassCount = 0;
};

static const size_t SOURCE_REALIZATION_MAX_MEMORY_BYPASSES = 8;
static const size_t SOURCE_REALIZATION_UNKNOWN_WORKING_SET_BYTES =
    256ULL * 1024ULL * 1024ULL;
static const size_t SOURCE_REALIZATION_LEAF_FIXED_WORKING_SET_BYTES =
    8ULL * 1024ULL * 1024ULL;
static const size_t SOURCE_REALIZATION_LEAF_SERIALIZED_COPY_COUNT = 3;

/*
 * A fixed reservation allowed concurrent giant-leaf imports even though each
 * immediately needed much more than the nominal 256 MiB.  Size direct imports
 * from the cheap directory record before opening them.  The caller separately
 * classifies combinations and bounded serialized coverage.
 */
static size_t
source_realization_leaf_working_set_bytes(const struct db_i *dbip,
	const struct directory *dp)
{
    if (!dp)
	return SOURCE_REALIZATION_UNKNOWN_WORKING_SET_BYTES;

    size_t encodedBytes = dp->d_len;
    if (dbip && db_version(dbip) < 5) {
	if (encodedBytes > SIZE_MAX / sizeof(union record))
	    return SIZE_MAX;
	encodedBytes *= sizeof(union record);
    }
    if (encodedBytes >
	(SIZE_MAX - SOURCE_REALIZATION_LEAF_FIXED_WORKING_SET_BYTES) /
	SOURCE_REALIZATION_LEAF_SERIALIZED_COPY_COUNT)
	return SIZE_MAX;
    return encodedBytes * SOURCE_REALIZATION_LEAF_SERIALIZED_COPY_COUNT +
	SOURCE_REALIZATION_LEAF_FIXED_WORKING_SET_BYTES;
}

static size_t
source_realization_request_working_set_bytes(
	const BObolSourceRealizationRequest &request, size_t capacityLimit)
{
    if (request.estimatedWorkingSetBytes)
	return request.estimatedWorkingSetBytes;
    if (!request.source || !request.snapshotSourceDatabase)
	return SOURCE_REALIZATION_UNKNOWN_WORKING_SET_BYTES;

    const char *path = request.source->path.getValue().getString();
    if (!path || !path[0])
	return SOURCE_REALIZATION_UNKNOWN_WORKING_SET_BYTES;

    struct db_full_path fullPath;
    db_full_path_init(&fullPath);
    const int pathValid = db_string_to_path(&fullPath,
	request.snapshotSourceDatabase, path) == 0;
    const struct directory *dp = pathValid ?
	DB_FULL_PATH_CUR_DIR(&fullPath) : NULL;
    const bool combination = dp && (dp->d_flags & RT_DIR_COMB);
    db_free_full_path(&fullPath);
    if (!dp)
	return SOURCE_REALIZATION_UNKNOWN_WORKING_SET_BYTES;
    /* Serialized BoT coverage has its own bounded scan/detail admission,
     * just like a combination's leaf census.  Reserve the complete source
     * allowance so directory setup runs alone, without charging a full
     * primitive import which this route expressly forbids. */
    if (combination || (request.stream &&
	bobol_database_source_uses_serialized_bot_coverage(
	    request.source, request.snapshotSourceDatabase)))
	return capacityLimit && capacityLimit != SIZE_MAX ? capacityLimit :
	    SOURCE_REALIZATION_UNKNOWN_WORKING_SET_BYTES;
    return source_realization_leaf_working_set_bytes(
	request.snapshotSourceDatabase, dp);
}

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
    return estimate ? estimate : SOURCE_REALIZATION_UNKNOWN_WORKING_SET_BYTES;
}

static bool
source_job_cancelled(const BObolSourceRealizationJobPrivate *job)
{
    return !job || job->cancelRequested.load(std::memory_order_acquire);
}

static bool
source_item_cancelled(const BObolSourceRealizationJobPrivate *job,
	const BObolSourceRealizationItemPrivate *item)
{
    return source_job_cancelled(job) || !item ||
	(item->stream && item->stream->isCancelled());
}

static void
source_job_request_cancel(BObolSourceRealizationJobPrivate *job)
{
    if (!job || job->cancelRequested.exchange(true, std::memory_order_acq_rel))
	return;
    /* One job may occupy many queue entries and worker slots.  Only its
     * first cancellation request needs to visit the immutable item vector. */
    for (const std::unique_ptr<BObolSourceRealizationItemPrivate> &item : job->items) {
	if (item && item->stream)
	    item->stream->requestCancel();
    }
    /* Requests without an occurrence stream still own queue resources. */
    bobol_source_realization_cancel_queued(job->coordinator);
}

static void
source_job_fail(const std::shared_ptr<BObolSourceRealizationJobPrivate> &job,
	size_t itemIndex)
{
    if (!job)
	return;
    job->failed.store(true, std::memory_order_release);
    source_job_request_cancel(job.get());
    if (itemIndex < job->items.size() && job->items[itemIndex])
	job->items[itemIndex]->finish(BOBOL_SOURCE_REALIZATION_FAILED);
}

static void
source_realize_item(const std::shared_ptr<BObolSourceRealizationJobPrivate> &job,
	size_t itemIndex)
{
    if (!job || itemIndex >= job->items.size())
	return;
    BObolSourceRealizationItemPrivate *item = job->items[itemIndex].get();
    if (!item)
	return;
    if (source_item_cancelled(job.get(), item)) {
	item->finish(BOBOL_SOURCE_REALIZATION_CANCELLED);
	return;
    }

    item->state.store(BOBOL_SOURCE_REALIZATION_RUNNING,
	std::memory_order_release);
    SoBRLDatabaseSource *source = item->source;
    bool success = source != NULL;

    int warmManifest = 0;
    if (success && item->probeWarmManifest && item->stream &&
	item->snapshotSourceDatabase && !source_item_cancelled(job.get(), item)) {
	warmManifest = item->probeWarmManifest(source,
	    item->snapshotSourceDatabase, item->drawMode,
	    source->sourceRevision.getValue(), item->stream.get(),
	    item->callbackContext.get());
    }
    item->warmManifest.store(warmManifest == 2, std::memory_order_release);
    if (item->warmManifest.load(std::memory_order_acquire) &&
	item->snapshotSourceDatabase) {
	db_close(item->snapshotSourceDatabase);
	item->snapshotSourceDatabase = NULL;
    }

    if (success && !item->warmManifest && !item->database) {
	if (source_item_cancelled(job.get(), item)) {
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

    if (success && !source_item_cancelled(job.get(), item)) {
	const bool mesh = source->usesMeshRealization() ? true : false;
	BObolCompactOccurrenceStream *stream =
	    item->stream ? item->stream.get() : NULL;
	SbBool realized = warmManifest == 2 ? TRUE :
	    bobol_database_source_construct_realization(source, mesh, stream);
	if (!realized && mesh && item->allowWireFallback &&
	    !source_item_cancelled(job.get(), item))
	    realized = bobol_database_source_construct_realization(source, FALSE, stream);
	success = realized ? true : false;
	if (success && stream && !source_item_cancelled(job.get(), item) &&
	    !stream->hasCoverageBoundsComplete() &&
	    source->hasExactSourceBounds()) {
	    SbBox3f finalBounds;
	    if (source->getSourceBounds(finalBounds) &&
		!finalBounds.isEmpty()) {
		/* Generic realization needs the same independent coverage payload
		 * as serialized BoTs. Terminal source adoption may never happen if
		 * the consumer cannot allocate storage for the leaves. */
		BObolCompactOccurrence overview =
		    bobol_database_source_coverage_overview(source,
			source->path.getValue().getString(), finalBounds,
			source->sourceRevision.getValue());
		if (overview.geometry) {
		    stream->setCoverageBounds(finalBounds);
		    stream->pushPriority(overview);
		    stream->setCoverageBoundsComplete(true);
		} else {
		    success = false;
		}
	    }
	}
	if (success &&
	    !item->warmManifest.load(std::memory_order_acquire) &&
	    !(stream && stream->hasWarmCensusComplete()) &&
	    item->storeManifest &&
	    !source_item_cancelled(job.get(), item))
	    item->manifestStored.store(
		item->storeManifest(item->database, source,
		    stream,
		    item->callbackContext.get()) ? true : false,
		std::memory_order_release);
    } else {
	success = false;
    }

    if (source_item_cancelled(job.get(), item)) {
	item->finish(BOBOL_SOURCE_REALIZATION_CANCELLED);
	return;
    }
    if (!success) {
	source_job_fail(job, itemIndex);
	return;
    }
    item->finish(BOBOL_SOURCE_REALIZATION_COMPLETE);
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

struct SourceRealizationWorker {
    std::thread thread;
    /* Accessed under the queue mutex so shutdown can cancel a job after its
     * work has left the queue and before its callback returns. */
    std::shared_ptr<BObolSourceRealizationJobPrivate> activeJob;
};

struct BObolSourceRealizationCoordinatorPrivate {
    BObolSourceRealizationCoordinatorPrivate(void) :
	stopped(false), stopping(false), active(0), activeBytes(0),
	queueRetirements(0),
	maxActiveBytes(bobol_lod_working_set_global_limit())
    {
    }

    std::mutex stopMutex;
    bool stopped;
    std::mutex mutex;
    std::condition_variable cv;
    bool stopping;
    size_t active;
    size_t activeBytes;
    /* Shutdown must also wait for queue owners already detached by callers,
     * whose source/context destructors may still use the cache registries. */
    size_t queueRetirements;
    size_t maxActiveBytes;
    std::list<SourceRealizationWork> queue;
    std::vector<SourceRealizationWorker> workers;
};

static void
source_job_retire_item(const std::shared_ptr<BObolSourceRealizationJobPrivate> &job)
{
    const size_t remaining = job ?
	job->remaining.fetch_sub(1, std::memory_order_acq_rel) : 0;
    if (job && remaining == 1) {
	const int terminal = job->failed.load(std::memory_order_acquire) ?
	    BOBOL_SOURCE_REALIZATION_FAILED :
	    (source_job_cancelled(job.get()) ?
		BOBOL_SOURCE_REALIZATION_CANCELLED :
		BOBOL_SOURCE_REALIZATION_COMPLETE);
	/* Callers retire queue ownership or active reservations first. Terminal
	 * publication certifies those effects and every item/source write. A
	 * failed item cancels siblings but remains an aggregate failure. */
	job->state.store(terminal, std::memory_order_release);
    }
}

void
bobol_source_realization_cancel_queued(
    const std::weak_ptr<BObolSourceRealizationCoordinatorPrivate> &coordinator)
{
    auto service = coordinator.lock();
    if (!service)
	return;
    std::list<SourceRealizationWork> retired;
    {
	std::lock_guard<std::mutex> lock(service->mutex);
	if (service->stopping)
	    return;
	for (auto entry = service->queue.begin(); entry != service->queue.end();) {
	    const auto &work = *entry;
	    if (source_item_cancelled(work.job.get(), work.job->items[work.itemIndex].get())) {
		/* Splicing transfers all selected ownership without allocating or
		 * destroying source/context storage under the queue mutex. */
		const auto cancelled = entry++;
		retired.splice(retired.end(), service->queue, cancelled);
	    } else {
		++entry;
	    }
	}
	if (retired.empty())
	    return;
	++service->queueRetirements;
    }
    /* Removing an aged root can admit healthy work before caller-owned
     * cleanup finishes. That cleanup must not withhold its peers' wakeup. */
    service->cv.notify_all();
    for (auto &work : retired) {
	work.job->items[work.itemIndex]->finish(BOBOL_SOURCE_REALIZATION_CANCELLED);
	source_job_retire_item(work.job);
    }
    /* Drop job ownership before certifying cleanup to process teardown. */
    retired.clear();
    {
	std::lock_guard<std::mutex> lock(service->mutex);
	--service->queueRetirements;
    }
    service->cv.notify_all();
}

static void
source_realization_worker(BObolSourceRealizationCoordinatorPrivate *service,
    size_t workerIndex)
{
    const auto fits = [service](const SourceRealizationWork &candidate) {
	const size_t estimate = source_work_estimated_bytes(candidate);
	return !service->active || service->maxActiveBytes == SIZE_MAX ||
	    (estimate <= service->maxActiveBytes &&
	     service->activeBytes <= service->maxActiveBytes - estimate);
    };
    for (;;) {
	SourceRealizationWork work;
	size_t admittedBytes = 0;
	{
	    std::unique_lock<std::mutex> lock(service->mutex);
	    service->cv.wait(lock, [service, &fits]() {
		if (service->stopping)
		    return true;
		for (const SourceRealizationWork &candidate : service->queue) {
		    if (fits(candidate))
			return true;
		    if (candidate.memoryBypassCount >=
			    SOURCE_REALIZATION_MAX_MEMORY_BYPASSES)
			return false;
		}
		return false;
	    });
	    if (service->stopping)
		return;
	    auto selected = service->queue.end();
	    for (auto candidate = service->queue.begin();
		    candidate != service->queue.end(); ++candidate) {
		if (fits(*candidate)) {
		    selected = candidate;
		    admittedBytes = source_work_estimated_bytes(*candidate);
		    break;
		}
		if (candidate->memoryBypassCount >=
			SOURCE_REALIZATION_MAX_MEMORY_BYPASSES)
		    break;
	    }
	    if (selected == service->queue.end()) {
		continue;
	    }
	    /* Only candidates actually skipped by this admission age.  A blocked
	     * candidate which later fits is consumed without carrying stale age,
	     * and unrelated work behind the selected item is unaffected. */
	    for (auto bypassed = service->queue.begin();
		    bypassed != selected; ++bypassed) {
		if (bypassed->memoryBypassCount != SIZE_MAX)
		    bypassed->memoryBypassCount++;
	    }
	    work = std::move(*selected);
	    service->queue.erase(selected);
	    service->workers[workerIndex].activeJob = work.job;
	    service->active++;
	    service->activeBytes =
		admittedBytes > SIZE_MAX - service->activeBytes ?
		SIZE_MAX : service->activeBytes + admittedBytes;
	}

	try {
	    /* Source discovery and mesh-LoD production may overlap during a cold
	     * draw.  Count both against the same process budget so each one's
	     * internal helpers cannot multiply the other's outer worker pool. */
	    BObolParallelBudgetLease cpuBudget;
	    cpuBudget.acquireOuter();
	    source_realize_item(work.job, work.itemIndex);
	} catch (const std::bad_alloc &) {
	    source_job_fail(work.job, work.itemIndex);
	    bu_log("BObol source realization ran out of memory\n");
	} catch (...) {
	    source_job_fail(work.job, work.itemIndex);
	    bu_log("BObol source realization failed with an unexpected exception\n");
	}

	{
	    std::lock_guard<std::mutex> lock(service->mutex);
	    service->active--;
	    service->activeBytes =
		admittedBytes >= service->activeBytes ?
		0 : service->activeBytes - admittedBytes;
	    service->workers[workerIndex].activeJob.reset();
	}

	source_job_retire_item(work.job);
	service->cv.notify_all();
    }
}

static void
source_realization_stop(BObolSourceRealizationCoordinatorPrivate *service)
{
    if (!service)
	return;
    std::unique_lock<std::mutex> stopLock(service->stopMutex);
    if (service->stopped)
	return;
    {
	std::lock_guard<std::mutex> lock(service->mutex);
	service->stopping = true;
    }
    service->cv.notify_all();
    for (SourceRealizationWorker &worker : service->workers) {
	std::shared_ptr<BObolSourceRealizationJobPrivate> activeJob;
	{
	    std::lock_guard<std::mutex> lock(service->mutex);
	    activeJob = worker.activeJob;
	}
	source_job_request_cancel(activeJob.get());
    }
    /* Stopping forbids further worker admission. Retire queued items here;
     * releasing their streams may invoke storage-owner cleanup, which must
     * run outside the coordinator lock just like worker callbacks. */
    for (;;) {
	SourceRealizationWork work;
	{
	    std::lock_guard<std::mutex> lock(service->mutex);
	    if (service->queue.empty())
		break;
	    work = std::move(service->queue.front());
	    service->queue.pop_front();
	}
	source_job_request_cancel(work.job.get());
	work.job->items[work.itemIndex]->finish(BOBOL_SOURCE_REALIZATION_CANCELLED);
	source_job_retire_item(work.job);
    }
    for (SourceRealizationWorker &worker : service->workers) {
	if (worker.thread.joinable())
	    worker.thread.join();
    }
    std::unique_lock<std::mutex> lock(service->mutex);
    service->cv.wait(lock, [service]() { return !service->queueRetirements; });
    service->stopped = true;
}

BObolSourceRealizationRequest::BObolSourceRealizationRequest(void) :
    source(NULL), snapshotSourceDatabase(NULL), clientToken(0),
    estimatedWorkingSetBytes(0), drawMode(0),
    allowWireFallback(FALSE), probeWarmManifest(NULL), storeManifest(NULL),
    callbackContext()
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
    /* Completed streams may have transferred their staging leases to a live
     * source. Dropping producer interest must preserve that consumer owner. */
    if (this->state() != BOBOL_SOURCE_REALIZATION_COMPLETE)
	this->cancel();
}

void
BObolSourceRealizationJob::cancel(void)
{
    source_job_request_cancel(this->p.get());
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
	current == BOBOL_SOURCE_REALIZATION_CANCELLED ||
	current == BOBOL_SOURCE_REALIZATION_CONSTRAINED ? TRUE : FALSE;
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
    const BObolSourceRealizationItemPrivate *item = this->p->items[index].get();
    result.stream = item->stream;
    result.clientToken = item->clientToken;
    result.state = item->state.load(std::memory_order_acquire);
    /* A running source is worker-exclusive. Only COMPLETE publishes an
     * immutable source pointer; other outcomes release it before publication. */
    if (result.state == BOBOL_SOURCE_REALIZATION_COMPLETE)
	result.source = item->source;
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
    p()
{
    bobol_draw_cache_runtime_prepare();
    bobol_mesh_lod_cache_runtime_prepare();
    std::shared_ptr<BObolSourceRealizationCoordinatorPrivate> candidate(
	new BObolSourceRealizationCoordinatorPrivate);
    const size_t count = source_realization_worker_count();
    candidate->workers.resize(count);
    const bool failWorkerStart = bobol_transaction_fault_requested(
	BObolTransactionFaultPoint::SOURCE_REALIZATION_WORKER_START);
    try {
	for (size_t i = 0; i < count; ++i) {
	    if (i > 0 && failWorkerStart)
		throw std::system_error(
		    std::make_error_code(std::errc::resource_unavailable_try_again));
	    candidate->workers[i].thread = std::thread(
		source_realization_worker, candidate.get(), i);
	}
    } catch (...) {
	/* Destruction of a joinable std::thread terminates the process.  Join
	 * every successfully started worker before discarding the candidate. */
	source_realization_stop(candidate.get());
	throw;
    }
    this->p = std::move(candidate);
}

BObolSourceRealizationCoordinator::~BObolSourceRealizationCoordinator(void)
{
    this->shutdown();
}

void
BObolSourceRealizationCoordinator::shutdown(void)
{
    source_realization_stop(this->p.get());
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

    try {
	std::shared_ptr<BObolSourceRealizationJobPrivate> job(
	    new BObolSourceRealizationJobPrivate);
	job->coordinator = this->p;
	if (bobol_transaction_fault_requested(
		BObolTransactionFaultPoint::SOURCE_REALIZATION_JOB_HANDLE))
	    throw std::bad_alloc();
	/* An uncommitted handle must not cancel the caller's streams when a
	 * later allocation fails.  Attach its interest only after transfer. */
	std::shared_ptr<BObolSourceRealizationJob> handle(
	    new BObolSourceRealizationJob(
		std::shared_ptr<BObolSourceRealizationJobPrivate>()));
	job->items.reserve(requests.size());
	/* Prepare every fallible allocation before crossing the ownership
	 * boundary.  A failed submission is specified to leave all request
	 * resources with the caller, including when process teardown has already
	 * closed the global queue. */
	for (const BObolSourceRealizationRequest &request : requests) {
	    std::unique_ptr<BObolSourceRealizationItemPrivate> item(
		new BObolSourceRealizationItemPrivate);
	    item->stream = request.stream;
	    item->clientToken = request.clientToken;
	    item->estimatedWorkingSetBytes =
		source_realization_request_working_set_bytes(request,
		    this->p->maxActiveBytes);
	    item->drawMode = request.drawMode;
	    item->allowWireFallback = request.allowWireFallback;
	    item->probeWarmManifest = request.probeWarmManifest;
	    item->storeManifest = request.storeManifest;
	    job->items.push_back(std::move(item));
	}

	/* A directory-derived estimate is a pre-import admission contract.  The
	 * historical "run one oversized root alone" escape still entered
	 * rt_db_get_internal, whose allocator reports an address-space failure
	 * with bu_bomb.  A display worker must instead leave the existing
	 * structural presentation intact and report an explicit constrained
	 * result.  Do this before any request can enter the worker queue. */
	bool admissionConstrained = false;
	if (this->p->maxActiveBytes != SIZE_MAX) {
	    for (const std::unique_ptr<BObolSourceRealizationItemPrivate> &item : job->items) {
		if (item && item->estimatedWorkingSetBytes >
		    this->p->maxActiveBytes) {
		    admissionConstrained = true;
		    break;
		}
	    }
	}

	bool cancelledBeforeBinding = false;
	{
	    std::lock_guard<std::mutex> lock(this->p->mutex);
	    if (this->p->stopping)
		return std::shared_ptr<BObolSourceRealizationJob>();
	    /* Prepare queue storage under the same lock as publication.  Roll back
	     * only this batch if queue growth fails, preserving earlier work and all
	     * caller resources.  Workers cannot observe the unfinished prefix. */
	    const size_t queuedBefore = this->p->queue.size();
	    try {
		if (!admissionConstrained) {
		    const bool failQueueAppend = bobol_transaction_fault_requested(
			BObolTransactionFaultPoint::SOURCE_REALIZATION_QUEUE_APPEND);
		    for (size_t i = 0; i < job->items.size(); ++i) {
			SourceRealizationWork work;
			work.job = job;
			work.itemIndex = i;
			/* Exercise rollback after a prefix already owns queue storage. */
			if (i > 0 && failQueueAppend)
			    throw std::bad_alloc();
			this->p->queue.push_back(std::move(work));
		    }
		}
	    } catch (...) {
		while (this->p->queue.size() > queuedBefore)
		    this->p->queue.pop_back();
		throw;
	    }
	    /* No allocating operation remains after this ownership boundary. */
	    for (size_t i = 0; i < requests.size(); ++i) {
		BObolSourceRealizationRequest &request = requests[i];
		BObolSourceRealizationItemPrivate *item = job->items[i].get();
		item->source = request.source;
		item->snapshotSourceDatabase = request.snapshotSourceDatabase;
		item->callbackContext = std::move(request.callbackContext);
		request.source = NULL;
		request.snapshotSourceDatabase = NULL;
		if (!admissionConstrained && item->stream &&
		    item->stream->bindRealizationCoordinator(this->p))
		    cancelledBeforeBinding = true;
	    }
	    handle->p = job;
	    if (!admissionConstrained) {
		job->remaining.store(job->items.size(), std::memory_order_release);
		job->state.store(BOBOL_SOURCE_REALIZATION_RUNNING,
		    std::memory_order_release);
	    }
	}
	if (admissionConstrained) {
	    /* Stream ownership release must not call out under the queue lock. */
	    for (const std::unique_ptr<BObolSourceRealizationItemPrivate> &item : job->items)
		item->finish(BOBOL_SOURCE_REALIZATION_CONSTRAINED);
	    /* A reused stream may also feed another queued job. Batch its queue
	     * notification after every constrained item's publication is closed. */
	    bobol_source_realization_cancel_queued(job->coordinator);
	    job->state.store(BOBOL_SOURCE_REALIZATION_CONSTRAINED,
		std::memory_order_release);
	    return handle;
	}
	if (cancelledBeforeBinding)
	    bobol_source_realization_cancel_queued(job->coordinator);
	this->p->cv.notify_all();
	return handle;
    } catch (const std::bad_alloc &) {
	bu_log("BObol source submission ran out of memory; preserving caller resources\n");
	return std::shared_ptr<BObolSourceRealizationJob>();
    }
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
