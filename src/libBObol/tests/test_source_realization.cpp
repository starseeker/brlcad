/*          T E S T _ S O U R C E _ R E A L I Z A T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BDrawCache.h"
#include "BObol/BInit.h"
#include "BObol/BSourceRealization.h"
#include "bu/app.h"
#include "bu/file.h"
#include "bu/str.h"
#include "rt/db_io.h"
#include "wdb.h"
#include "transaction_fault_test_private.h"

#include <Inventor/nodes/SoSeparator.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <memory>
#include <new>
#include <string>
#include <system_error>
#include <thread>
#include <vector>

namespace {

static const size_t source_admission_leaf_copy_count = 3;
static const size_t source_admission_leaf_fixed_bytes =
    8ULL * 1024ULL * 1024ULL;
static constexpr size_t pool_test_worker_count = 2;
/* Exceed the protected-root bypass allowance to establish a blocked frontier. */
static constexpr size_t fairness_small_request_count = 24;

struct ProbeGate {
    std::atomic<bool> entered{false};
    std::atomic<bool> release{false};
};

struct CountedProbeGate {
    std::atomic<size_t> entered{0};
    std::atomic<bool> release{false};
};

struct CompletionCounter {
    std::atomic<size_t> completed{0};
};

struct ShutdownCacheProbe {
    std::atomic<bool> entered{false};
};

static int
warm_complete_probe(SoBRLDatabaseSource *, struct db_i *, int, uint32_t,
    BObolCompactOccurrenceStream *, void *)
{
    return 2;
}

static int
blocking_warm_probe(SoBRLDatabaseSource *, struct db_i *, int, uint32_t,
    BObolCompactOccurrenceStream *, void *data)
{
    ProbeGate *gate = static_cast<ProbeGate *>(data);
    if (!gate)
	return 0;
    gate->entered.store(true, std::memory_order_release);
    while (!gate->release.load(std::memory_order_acquire))
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    return 2;
}

static int
counted_blocking_warm_probe(SoBRLDatabaseSource *, struct db_i *, int,
    uint32_t, BObolCompactOccurrenceStream *, void *data)
{
    CountedProbeGate *gate = static_cast<CountedProbeGate *>(data);
    if (!gate)
	return 0;
    gate->entered.fetch_add(1, std::memory_order_acq_rel);
    while (!gate->release.load(std::memory_order_acquire))
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    return 2;
}

static int
counting_complete_probe(SoBRLDatabaseSource *, struct db_i *, int, uint32_t,
    BObolCompactOccurrenceStream *, void *data)
{
    CompletionCounter *counter = static_cast<CompletionCounter *>(data);
    if (counter)
	counter->completed.fetch_add(1, std::memory_order_acq_rel);
    return 2;
}

static int
shutdown_cache_probe(SoBRLDatabaseSource *, struct db_i *database, int,
    uint32_t, BObolCompactOccurrenceStream *, void *data)
{
    ShutdownCacheProbe *probe = static_cast<ShutdownCacheProbe *>(data);
    if (!probe || !database)
	return 0;
    probe->entered.store(true, std::memory_order_release);
    /* Let main return and static teardown begin before touching the cache.
     * The coordinator destructor must be the barrier which keeps its lazy
     * registries alive until this callback has returned. */
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    BObolDrawLodAssetRecord record;
    (void)bobol_draw_lod_asset_cache_get(database,
	"__source_realization_shutdown_probe__", &record);
    return 2;
}

static bool
wait_until(const std::function<bool(void)> &predicate,
    std::chrono::milliseconds timeout)
{
    const std::chrono::steady_clock::time_point deadline =
	std::chrono::steady_clock::now() + timeout;
    while (!predicate()) {
	if (std::chrono::steady_clock::now() >= deadline)
	    return false;
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    return true;
}

static bool
make_request(BObolSourceRealizationRequest &request,
    struct db_i *database, BObolSourceWarmManifestProbe probe,
    const std::shared_ptr<void> &context, const char *sourcePath = NULL)
{
    if (!database)
	return false;
    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->ref();
    source->setDatabase(database);
    if (sourcePath && sourcePath[0])
	source->path = sourcePath;
    SoBRLDatabaseSource *detached =
	source->createDetachedRealizationTemplate();
    source->unref();
    if (!detached)
	return false;
    struct db_i *snapshot = db_clone_dbi(database, NULL);
    if (!snapshot) {
	detached->unref();
	return false;
    }
    request.source = detached;
    request.snapshotSourceDatabase = snapshot;
    request.stream = std::make_shared<BObolCompactOccurrenceStream>();
    request.probeWarmManifest = probe;
    request.callbackContext = context;
    return true;
}

static void
release_requests(std::vector<BObolSourceRealizationRequest> &requests)
{
    for (BObolSourceRealizationRequest &request : requests) {
	if (request.source)
	    request.source->unref();
	if (request.snapshotSourceDatabase)
	    db_close(request.snapshotSourceDatabase);
	request.source = NULL;
	request.snapshotSourceDatabase = NULL;
    }
}

static std::shared_ptr<BObolStagedSourceMesh>
make_retained_triangle(void)
{
    struct RetainedTriangle {
	point_t points[3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
	int faces[3] = {0, 1, 2};
    };
    auto triangle = std::make_shared<RetainedTriangle>();
    auto staged = std::make_shared<BObolStagedSourceMesh>();
    staged->owner = triangle;
    staged->points = triangle->points;
    staged->faces = triangle->faces;
    staged->pointCount = 3;
    staged->faceCount = 1;
    staged->assetName = "retained-import";
    staged->byteCount = sizeof(RetainedTriangle);
    return staged;
}

static int
test_completed_stream_transfer(BObolSourceRealizationCoordinator &coordinator,
    struct db_i *database)
{
    std::vector<BObolSourceRealizationRequest> requests(1);
    if (!make_request(requests[0], database, warm_complete_probe, {})) {
	std::fprintf(stderr, "FAIL: completed stream transfer setup\n");
	return 1;
    }
    auto stream = requests[0].stream;
    auto staged = make_retained_triangle();
    std::weak_ptr<const BObolStagedSourceMesh> retainedImport = staged;
    if (!stream->retainStagedSource(staged)) {
	release_requests(requests);
	std::fprintf(stderr, "FAIL: completed stream import setup\n");
	return 1;
    }
    staged.reset();
    auto job = coordinator.submit(requests);
    if (!job) {
	release_requests(requests);
	std::fprintf(stderr, "FAIL: completed stream transfer submission\n");
	return 1;
    }
    if (!wait_until([&]() { return job->isTerminal(); }, std::chrono::seconds(2)) ||
	job->state() != BOBOL_SOURCE_REALIZATION_COMPLETE) {
	std::fprintf(stderr, "FAIL: completed stream transfer did not finish\n");
	return 1;
    }
    job.reset();
    auto claimed = stream->claimStagedSource(retainedImport);
    if (stream->isCancelled() || !claimed || stream->stagedSourceByteCount()) {
	std::fprintf(stderr, "FAIL: completed job destruction revoked consumer staging\n");
	return 1;
    }
    stream->requestCancel();
    if (retainedImport.expired()) {
	std::fprintf(stderr, "FAIL: stream cancellation revoked an already claimed import\n");
	return 1;
    }
    claimed.reset();
    if (!retainedImport.expired()) {
	std::fprintf(stderr, "FAIL: claimed import retained an unexpected owner\n");
	return 1;
    }
    return 0;
}

struct ItemLifetimeEvidence {
    std::atomic<bool> sourceReleased{false};
    std::atomic<bool> contextReleased{false};
};

class ItemLifetimeNode : public SoSeparator {
public:
    explicit ItemLifetimeNode(const std::shared_ptr<ItemLifetimeEvidence> &value) :
	evidence(value) {}
    ~ItemLifetimeNode() override
    {
	evidence->sourceReleased.store(true, std::memory_order_release);
    }
private:
    std::shared_ptr<ItemLifetimeEvidence> evidence;
};

struct ItemLifetimeProbe {
    BObolSourceRealizationCoordinator *coordinator = nullptr;
    std::shared_ptr<ProbeGate> gate;
    std::shared_ptr<ItemLifetimeEvidence> evidence;
    int outcome = BOBOL_SOURCE_REALIZATION_COMPLETE;
    ~ItemLifetimeProbe()
    {
	/* User-owned callback storage may inspect accounting during cleanup. */
	(void)coordinator->activeWorkingSetBytesForDiagnostics();
	evidence->contextReleased.store(true, std::memory_order_release);
    }
};

static int
item_lifetime_probe(SoBRLDatabaseSource *source, struct db_i *database,
    int mode, uint32_t revision, BObolCompactOccurrenceStream *stream, void *data)
{
    auto *probe = static_cast<ItemLifetimeProbe *>(data);
    /* Install after cold realization has replaced the source's children. */
    source->addChild(new ItemLifetimeNode(probe->evidence));
    (void)blocking_warm_probe(source, database, mode, revision, stream,
	probe->gate.get());
    if (probe->outcome == BOBOL_SOURCE_REALIZATION_FAILED)
	throw std::bad_alloc();
    return 2;
}

static int
item_lifetime_store(struct db_i *database, const SoBRLDatabaseSource *source,
    BObolCompactOccurrenceStream *stream, void *data)
{
    return item_lifetime_probe(const_cast<SoBRLDatabaseSource *>(source),
	database, 0, 0, stream, data);
}

#ifdef __linux__
static size_t
database_open_file_count(const char *path)
{
    char **files = NULL;
    const size_t count = bu_file_list("/proc/self/fd", "[0-9]*", &files);
    size_t matches = 0;
    for (size_t i = 0; i < count; ++i) {
	const std::string descriptor = std::string("/proc/self/fd/") + files[i];
	char resolved[MAXPATHLEN] = {0};
	if (bu_file_realpath(descriptor.c_str(), resolved) &&
	    bu_strcmp(resolved, path) == 0)
	    ++matches;
    }
    bu_argv_free(count, files);
    return matches;
}
#endif

static int
test_item_resource_retirement(BObolSourceRealizationCoordinator &coordinator,
    struct db_i *siblingDatabase, int outcome, bool cold)
{
    char path[MAXPATHLEN] = {0};
    FILE *file = bu_temp_file(path, sizeof(path));
    if (!file)
	return 1;
    std::fclose(file);
    struct db_i *database = db_create(path, 5);
    if (!database) {
	(void)bu_file_delete(path);
	return 1;
    }
    auto evidence = std::make_shared<ItemLifetimeEvidence>();
    auto gate = std::make_shared<ProbeGate>();
    auto siblingGate = std::make_shared<ProbeGate>();
    auto probe = std::make_shared<ItemLifetimeProbe>();
    probe->coordinator = &coordinator;
    probe->gate = gate;
    probe->evidence = evidence;
    probe->outcome = outcome;
    std::vector<BObolSourceRealizationRequest> requests(2);
    struct rt_wdb *writer = wdb_dbopen(database, RT_WDB_TYPE_DB_DISK);
    point_t center = VINIT_ZERO;
    const bool setup = writer && mk_sph(writer, "retirement.s", center, 1.0) == 0 &&
	make_request(requests[0], database, cold ? NULL : item_lifetime_probe,
	    probe, "retirement.s") &&
	make_request(requests[1], siblingDatabase, blocking_warm_probe,
	    siblingGate, "admission.s");
    if (!setup) {
	release_requests(requests);
	db_close(database);
	(void)bu_file_delete(path);
	std::fprintf(stderr, "FAIL: item lifetime fixture setup\n");
	return 1;
    }
    if (cold)
	requests[0].storeManifest = item_lifetime_store;
    auto stream = requests[0].stream;
    auto job = coordinator.submit(requests);
    release_requests(requests);
    db_close(database);
    probe.reset();
    int failures = 0;
    const bool entered = job && wait_until([&]() {
	return gate->entered.load(std::memory_order_acquire) &&
	    siblingGate->entered.load(std::memory_order_acquire);
    }, std::chrono::seconds(2));
    BObolSourceRealizationItemResult result;
    if (!entered || !job->itemResult(0, result) || result.source ||
	evidence->sourceReleased.load(std::memory_order_acquire) ||
	evidence->contextReleased.load(std::memory_order_acquire)) {
	std::fprintf(stderr, "FAIL: running item exposed source or retired callback ownership\n");
	++failures;
    }
#ifdef __linux__
    if (entered && !database_open_file_count(path)) {
	std::fprintf(stderr, "FAIL: active item did not retain its database file\n");
	++failures;
    }
#endif
    if (outcome == BOBOL_SOURCE_REALIZATION_CANCELLED)
	stream->requestCancel();
    gate->release.store(true, std::memory_order_release);
    const bool retired = job && wait_until([&]() {
	return job->itemResult(0, result) && result.state == outcome &&
	    coordinator.activeItemCountForDiagnostics() == 1;
    }, std::chrono::seconds(2));
    const bool complete = outcome == BOBOL_SOURCE_REALIZATION_COMPLETE;
    if (!retired || job->isTerminal() ||
	!evidence->contextReleased.load(std::memory_order_acquire) ||
	evidence->sourceReleased.load(std::memory_order_acquire) == complete ||
	(result.source != NULL) != complete) {
	std::fprintf(stderr, "FAIL: item retained worker ownership behind its sibling "
	    "(cold=%d outcome=%d item=%d context-released=%d source-released=%d)\n",
	    cold, outcome, result.state,
	    evidence->contextReleased.load(std::memory_order_acquire),
	    evidence->sourceReleased.load(std::memory_order_acquire));
	++failures;
    }
#ifdef __linux__
    const bool needsDatabase = complete && cold;
    if (retired && (database_open_file_count(path) != 0) != needsDatabase) {
	std::fprintf(stderr, "FAIL: terminal item database ownership (cold=%d outcome=%d)\n",
	    cold, outcome);
	++failures;
    }
#endif
    SoBRLDatabaseSource *completedSource = result.source;
    if (complete && retired) {
	job->cancel();
	if (!job->itemResult(0, result) || !result.source || result.source != completedSource ||
	    evidence->sourceReleased.load(std::memory_order_acquire) ||
	    (cold && !db_lookup(result.source->getDatabase(), "retirement.s", LOOKUP_QUIET))) {
	    std::fprintf(stderr, "FAIL: cancellation revoked a borrowed completed result\n");
	    ++failures;
	}
    }
    siblingGate->release.store(true, std::memory_order_release);
    if (job && !wait_until([&]() { return job->isTerminal(); }, std::chrono::seconds(2))) {
	std::fprintf(stderr, "FAIL: item lifetime batch did not retire\n");
	++failures;
    }
    job.reset();
    if (!wait_until([&]() { return evidence->sourceReleased.load(std::memory_order_acquire); },
	std::chrono::seconds(2))) {
	std::fprintf(stderr, "FAIL: dropped item result retained its source\n");
	++failures;
    }
#ifdef __linux__
    if (database_open_file_count(path)) {
	std::fprintf(stderr, "FAIL: dropped item result retained its database file\n");
	++failures;
    }
#endif
    (void)bu_file_delete(path);
    return failures;
}

static int
test_individual_stream_cancellation(BObolSourceRealizationCoordinator &coordinator,
    struct db_i *database)
{
    auto cancelledGate = std::make_shared<ProbeGate>();
    auto healthyGate = std::make_shared<ProbeGate>();
    const auto coldProbe = [](SoBRLDatabaseSource *source, struct db_i *db,
	int mode, uint32_t revision, BObolCompactOccurrenceStream *stream, void *data) {
	(void)blocking_warm_probe(source, db, mode, revision, stream, data);
	return 0;
    };
    std::vector<BObolSourceRealizationRequest> requests(2);
    if (!make_request(requests[0], database, coldProbe, cancelledGate, "admission.bot") ||
	!make_request(requests[1], database, blocking_warm_probe, healthyGate, "admission.bot")) {
	release_requests(requests);
	std::fprintf(stderr, "FAIL: individual stream cancellation setup\n");
	return 1;
    }
    auto cancelledStream = requests[0].stream;
    auto healthyStream = requests[1].stream;
    auto job = coordinator.submit(requests);
    if (!job) {
	release_requests(requests);
	std::fprintf(stderr, "FAIL: individual stream cancellation submission\n");
	return 1;
    }
    int failures = 0;
    if (!wait_until([&]() { return cancelledGate->entered.load(std::memory_order_acquire); },
	std::chrono::seconds(2))) {
	std::fprintf(stderr, "FAIL: cancelled source did not enter its producer\n");
	failures++;
    }
    auto staged = make_retained_triangle();
    std::weak_ptr<const BObolStagedSourceMesh> retainedImport = staged;
    if (!cancelledStream->retainStagedSource(staged)) {
	std::fprintf(stderr, "FAIL: cancellation import lease setup\n");
	failures++;
    }
    staged.reset();
    BObolCompactOccurrence queued;
    cancelledStream->push(queued);
    cancelledStream->requestCancel();
    cancelledGate->release.store(true, std::memory_order_release);
    BObolSourceRealizationItemResult cancelledResult;
    if (!wait_until([&]() {
	return job->itemResult(0, cancelledResult) &&
	    cancelledResult.state != BOBOL_SOURCE_REALIZATION_PENDING &&
	    cancelledResult.state != BOBOL_SOURCE_REALIZATION_RUNNING &&
	    healthyGate->entered.load(std::memory_order_acquire) &&
	    coordinator.activeItemCountForDiagnostics() == 1 &&
	    coordinator.queuedItemCountForDiagnostics() == 0;
    }, std::chrono::seconds(2)) ||
	cancelledResult.state != BOBOL_SOURCE_REALIZATION_CANCELLED ||
	healthyStream->isCancelled() || job->isTerminal()) {
	std::fprintf(stderr, "FAIL: per-source cancellation escaped its item "
	    "(item=%d peer-cancelled=%d job=%d)\n", cancelledResult.state,
	    healthyStream->isCancelled(), job->state());
	failures++;
    }
    if (cancelledStream->size() || cancelledStream->stagedSourceByteCount() ||
	!retainedImport.expired()) {
	std::fprintf(stderr, "FAIL: cancelled item retained payload while its sibling runs "
	    "(queued=%zu staged=%zu import-live=%d active-bytes=%zu)\n",
	    cancelledStream->size(), cancelledStream->stagedSourceByteCount(),
	    !retainedImport.expired(), coordinator.activeWorkingSetBytesForDiagnostics());
	failures++;
    }
    if (cancelledStream->retainStagedSource(make_retained_triangle())) {
	std::fprintf(stderr, "FAIL: cancelled stream accepted a late valid import\n");
	failures++;
    }
    healthyGate->release.store(true, std::memory_order_release);
    BObolSourceRealizationItemResult healthyResult;
    if (!wait_until([&]() { return job->isTerminal(); }, std::chrono::seconds(2)) ||
	job->state() != BOBOL_SOURCE_REALIZATION_COMPLETE ||
	!job->itemResult(1, healthyResult) ||
	healthyResult.state != BOBOL_SOURCE_REALIZATION_COMPLETE ||
	!healthyGate->entered.load(std::memory_order_acquire) ||
	healthyStream->isCancelled() ||
	coordinator.activeItemCountForDiagnostics() ||
	coordinator.queuedItemCountForDiagnostics() ||
	coordinator.activeWorkingSetBytesForDiagnostics()) {
	std::fprintf(stderr, "FAIL: healthy source did not finish after peer cancellation\n");
	failures++;
    }
    return failures;
}

enum class QueuedCancellation {
    Stream,
    Job,
    Interest,
    Precancelled,
    StreamlessJob,
    ConcurrentStreams,
    SharedStreamJob,
    SharedStreamConstrained
};

static int
test_queued_cancellation(BObolSourceRealizationCoordinator &coordinator,
    struct db_i *database, bool memoryBlocked, QueuedCancellation cancellation)
{
    const size_t allowance = coordinator.workingSetLimitBytesForDiagnostics();
    if (memoryBlocked && (!allowance || allowance == SIZE_MAX))
	return 0;
    const bool constrained = cancellation == QueuedCancellation::SharedStreamConstrained;
    if (constrained && allowance == SIZE_MAX)
	return 0;
    const size_t blockerCount = memoryBlocked ? 1 :
	coordinator.workerCountForDiagnostics();
    auto blockerGate = std::make_shared<CountedProbeGate>();
    std::vector<BObolSourceRealizationRequest> blockers(blockerCount);
    for (auto &request : blockers) {
	if (!make_request(request, database, counted_blocking_warm_probe, blockerGate)) {
	    release_requests(blockers);
	    return 1;
	}
	request.estimatedWorkingSetBytes = memoryBlocked ? allowance : 1;
    }
    auto blockerJob = coordinator.submit(blockers);
    if (!blockerJob || !wait_until([&]() {
	return blockerGate->entered.load(std::memory_order_acquire) == blockerCount;
    }, std::chrono::seconds(2))) {
	blockerGate->release.store(true, std::memory_order_release);
	release_requests(blockers);
	std::fprintf(stderr, "FAIL: queued cancellation did not establish blocked admission\n");
	return 1;
    }

    constexpr size_t queuedItemCount = 3;
    constexpr size_t individualItem = 1;
    std::vector<BObolSourceRealizationRequest> requests(queuedItemCount);
    std::vector<std::shared_ptr<ItemLifetimeEvidence>> evidence;
    auto callbackGate = std::make_shared<ProbeGate>();
    callbackGate->release.store(true, std::memory_order_release);
    bool setup = true;
    for (auto &request : requests) {
	auto current = std::make_shared<ItemLifetimeEvidence>();
	auto probe = std::make_shared<ItemLifetimeProbe>();
	probe->coordinator = &coordinator;
	probe->gate = callbackGate;
	probe->evidence = current;
	if (!make_request(request, database, item_lifetime_probe, probe)) {
	    setup = false;
	    break;
	}
	request.source->addChild(new ItemLifetimeNode(current));
	request.estimatedWorkingSetBytes = constrained ? allowance + 1 : 1;
	if (cancellation == QueuedCancellation::StreamlessJob)
	    request.stream.reset();
	if (cancellation == QueuedCancellation::Precancelled)
	    request.stream->requestCancel();
	evidence.push_back(current);
    }
    std::shared_ptr<BObolSourceRealizationJob> sharedPeer;
    if (setup && (constrained || cancellation == QueuedCancellation::SharedStreamJob)) {
	std::vector<BObolSourceRealizationRequest> peer(1);
	setup = make_request(peer[0], database, warm_complete_probe, {});
	if (setup) {
	    peer[0].stream = requests[individualItem].stream;
	    peer[0].estimatedWorkingSetBytes = 1;
	    sharedPeer = coordinator.submit(peer);
	    setup = sharedPeer != nullptr;
	}
	release_requests(peer);
    }
    auto job = setup ? coordinator.submit(requests) :
	std::shared_ptr<BObolSourceRealizationJob>();
    release_requests(requests);
    int failures = 0;
    if (!job) {
	std::fprintf(stderr, "FAIL: queued cancellation submission\n");
	++failures;
    } else {
	switch (cancellation) {
	    case QueuedCancellation::Stream:
		requests[individualItem].stream->requestCancel();
		break;
	    case QueuedCancellation::Job:
	    case QueuedCancellation::StreamlessJob:
	    case QueuedCancellation::SharedStreamJob:
		job->cancel();
		break;
	    case QueuedCancellation::Interest:
		job.reset();
		break;
	    case QueuedCancellation::Precancelled:
	    case QueuedCancellation::SharedStreamConstrained:
		break;
	    case QueuedCancellation::ConcurrentStreams: {
		std::vector<std::thread> cancellers;
		try {
		    cancellers.reserve(requests.size());
		    for (auto &request : requests)
			cancellers.emplace_back([stream = request.stream]() { stream->requestCancel(); });
		} catch (const std::exception &error) {
		    std::fprintf(stderr, "FAIL: concurrent cancellation setup: %s\n", error.what());
		    ++failures;
		}
		for (auto &thread : cancellers)
		    thread.join();
		break;
	    }
	}
	const bool individual = cancellation == QueuedCancellation::Stream;
	const size_t remaining = individual ? queuedItemCount - 1 : 0;
	const bool retired = wait_until([&]() {
	    if (coordinator.queuedItemCountForDiagnostics() != remaining)
		return false;
	    for (size_t i = 0; i < evidence.size(); ++i) {
		if (individual && i != individualItem)
		    continue;
		if (!evidence[i]->sourceReleased.load(std::memory_order_acquire) ||
		    !evidence[i]->contextReleased.load(std::memory_order_acquire))
		    return false;
		if (job) {
		    BObolSourceRealizationItemResult result;
		    if (!job->itemResult(i, result) || result.source ||
			result.state != (constrained ? BOBOL_SOURCE_REALIZATION_CONSTRAINED :
			    BOBOL_SOURCE_REALIZATION_CANCELLED))
			return false;
		}
	    }
	    return true;
	}, std::chrono::seconds(2));
	if (!retired || callbackGate->entered.load(std::memory_order_acquire) ||
	    coordinator.activeItemCountForDiagnostics() != blockerCount ||
	    coordinator.activeWorkingSetBytesForDiagnostics() !=
		(memoryBlocked ? allowance : blockerCount) || blockerJob->isTerminal()) {
	    std::fprintf(stderr, "FAIL: cancellation waited for unrelated admission "
		"(memory-blocked=%d kind=%d queued=%zu expected=%zu)\n",
		memoryBlocked, static_cast<int>(cancellation),
		coordinator.queuedItemCountForDiagnostics(), remaining);
	    ++failures;
	}
	if (sharedPeer) {
	    BObolSourceRealizationItemResult peerResult;
	    if (!sharedPeer->isTerminal() || !sharedPeer->itemResult(0, peerResult) ||
		peerResult.state != BOBOL_SOURCE_REALIZATION_CANCELLED || peerResult.source) {
		std::fprintf(stderr, "FAIL: shared stream left a queued peer after cancellation\n");
		++failures;
	    }
	}
	if (individual && (job->isTerminal() ||
	    requests.front().stream->isCancelled() || requests.back().stream->isCancelled() ||
	    evidence.front()->contextReleased.load(std::memory_order_acquire) ||
	    evidence.back()->contextReleased.load(std::memory_order_acquire))) {
	    std::fprintf(stderr, "FAIL: queued stream cancellation retired healthy siblings\n");
	    ++failures;
	}
    }
    blockerGate->release.store(true, std::memory_order_release);
    if (!wait_until([&]() {
	return blockerJob->isTerminal() && (!job || job->isTerminal()) &&
	    (!sharedPeer || sharedPeer->isTerminal()) &&
	    !coordinator.activeItemCountForDiagnostics() &&
	    !coordinator.queuedItemCountForDiagnostics();
    }, std::chrono::seconds(2))) {
	std::fprintf(stderr, "FAIL: queued cancellation fixture did not drain\n");
	++failures;
    }
    if (job && cancellation == QueuedCancellation::Stream &&
	(job->state() != BOBOL_SOURCE_REALIZATION_COMPLETE ||
	 !callbackGate->entered.load(std::memory_order_acquire))) {
	std::fprintf(stderr, "FAIL: healthy queued siblings did not complete\n");
	++failures;
    }
    return failures;
}

static int
test_submission_allocation_failure(BObolSourceRealizationCoordinator &coordinator,
    struct db_i *database, BObolTransactionFaultPoint fault, bool constrained)
{
    const size_t limit = coordinator.workingSetLimitBytesForDiagnostics();
    if (constrained && limit == SIZE_MAX)
	return 0;
    const size_t blockerCount = limit == SIZE_MAX ?
	coordinator.workerCountForDiagnostics() : 1;
    auto gate = std::make_shared<CountedProbeGate>();
    std::vector<BObolSourceRealizationRequest> blockers(blockerCount);
    for (BObolSourceRealizationRequest &request : blockers) {
	if (!make_request(request, database, counted_blocking_warm_probe, gate)) {
	    release_requests(blockers);
	    std::fprintf(stderr, "FAIL: submission fault blocker setup\n");
	    return 1;
	}
	request.estimatedWorkingSetBytes = limit == SIZE_MAX ? 1 : limit;
    }
    std::shared_ptr<BObolSourceRealizationJob> blocker = coordinator.submit(blockers);
    if (!blocker || !wait_until([&]() {
	    return gate->entered.load(std::memory_order_acquire) == blockerCount;
	}, std::chrono::seconds(2))) {
	gate->release.store(true, std::memory_order_release);
	release_requests(blockers);
	std::fprintf(stderr, "FAIL: submission fault blocker did not start\n");
	return 1;
    }

    int failures = 0;
    auto sentinelCounter = std::make_shared<CompletionCounter>();
    auto rejectedCounter = std::make_shared<CompletionCounter>();
    std::vector<BObolSourceRealizationRequest> sentinelRequests(1);
    /* Two items force queue failure after one entry has been appended. */
    std::vector<BObolSourceRealizationRequest> requests(2);
    bool setup = make_request(sentinelRequests[0], database,
	counting_complete_probe, sentinelCounter);
    sentinelRequests[0].estimatedWorkingSetBytes = 1;
    for (BObolSourceRealizationRequest &request : requests) {
	if (!make_request(request, database, counting_complete_probe, rejectedCounter))
	    setup = false;
	request.estimatedWorkingSetBytes = constrained ? limit + 1 : 1;
    }
    std::shared_ptr<BObolSourceRealizationJob> sentinel;
    std::shared_ptr<BObolSourceRealizationJob> retry;
    if (!setup) {
	std::fprintf(stderr, "FAIL: submission fault request setup\n");
	failures++;
    } else {
	sentinel = coordinator.submit(sentinelRequests);
	const size_t queuedBefore = coordinator.queuedItemCountForDiagnostics();
	const size_t activeBytesBefore = coordinator.activeWorkingSetBytesForDiagnostics();
	const std::vector<BObolSourceRealizationRequest> original = requests;
	std::shared_ptr<BObolSourceRealizationJob> rejected;
	bool threw = false;
	{
	    ScopedTransactionFault injected(fault);
	    try {
		rejected = coordinator.submit(requests);
	    } catch (const std::bad_alloc &) {
		threw = true;
	    }
	}
	bool preserved = true;
	for (size_t i = 0; i < requests.size(); ++i) {
	    preserved = preserved && requests[i].source == original[i].source &&
		requests[i].snapshotSourceDatabase == original[i].snapshotSourceDatabase &&
		requests[i].stream == original[i].stream &&
		requests[i].callbackContext == original[i].callbackContext &&
		!requests[i].stream->isCancelled();
	}
	const size_t queuedAfter = coordinator.queuedItemCountForDiagnostics();
	if (!sentinel || queuedBefore != sentinelRequests.size() || rejected ||
	    threw || !preserved || queuedAfter != queuedBefore ||
	    coordinator.activeItemCountForDiagnostics() != blockerCount ||
	    coordinator.activeWorkingSetBytesForDiagnostics() != activeBytesBefore ||
	    rejectedCounter->completed.load(std::memory_order_acquire) != 0) {
	    std::fprintf(stderr,
		"FAIL: source submission fault %s (constrained=%d) escaped transaction "
		"(threw=%d preserved=%d queued=%zu/%zu)\n",
		bobol_transaction_fault_name(fault), constrained, threw, preserved,
		queuedAfter, queuedBefore);
	    failures++;
	} else {
	    retry = coordinator.submit(requests);
	    if (!retry) {
		std::fprintf(stderr, "FAIL: preserved source batch could not be retried\n");
		failures++;
	    }
	}
    }

    gate->release.store(true, std::memory_order_release);
    if (!wait_until([&]() {
	    return blocker->isTerminal() && (!sentinel || sentinel->isTerminal()) &&
		(!retry || retry->isTerminal()) &&
		coordinator.activeItemCountForDiagnostics() == 0 &&
		coordinator.queuedItemCountForDiagnostics() == 0;
	}, std::chrono::seconds(2))) {
	std::fprintf(stderr, "FAIL: source submission fault work did not retire\n");
	failures++;
    }
    const size_t expectedCallbacks = retry && !constrained ? requests.size() : 0;
    if ((sentinel && (sentinel->state() != BOBOL_SOURCE_REALIZATION_COMPLETE ||
	    sentinelCounter->completed.load(std::memory_order_acquire) != 1)) ||
	rejectedCounter->completed.load(std::memory_order_acquire) != expectedCallbacks ||
	coordinator.activeWorkingSetBytesForDiagnostics() != 0 ||
	(retry && retry->state() != (constrained ? BOBOL_SOURCE_REALIZATION_CONSTRAINED :
	    BOBOL_SOURCE_REALIZATION_COMPLETE))) {
	std::fprintf(stderr, "FAIL: source submission fault corrupted completion or reservations\n");
	failures++;
    }
    release_requests(requests);
    release_requests(sentinelRequests);
    return failures;
}

struct QueuedCleanupProgressProbe {
    std::shared_ptr<CompletionCounter> healthy;
    std::atomic<bool> *healthyCompleted = nullptr;
    ~QueuedCleanupProgressProbe()
    {
	healthyCompleted->store(wait_until([&]() {
	    return healthy->completed.load(std::memory_order_acquire) ==
		fairness_small_request_count;
	}, std::chrono::seconds(2)), std::memory_order_release);
    }
};

static int
test_cancelled_frontier_wakeup(BObolSourceRealizationCoordinator &coordinator,
    struct db_i *database)
{
    const size_t limit = coordinator.workingSetLimitBytesForDiagnostics();
    if (!limit || limit == SIZE_MAX || coordinator.workerCountForDiagnostics() < 2)
	return 0;
    auto blockerGate = std::make_shared<ProbeGate>();
    auto healthy = std::make_shared<CompletionCounter>();
    std::atomic<bool> healthyCompleted{false};
    auto cleanup = std::make_shared<QueuedCleanupProgressProbe>();
    cleanup->healthy = healthy;
    cleanup->healthyCompleted = &healthyCompleted;
    std::vector<BObolSourceRealizationRequest> blockers(1);
    std::vector<BObolSourceRealizationRequest> mixed(fairness_small_request_count + 1);
    bool setup = make_request(blockers[0], database, blocking_warm_probe, blockerGate) &&
	make_request(mixed[0], database, warm_complete_probe, cleanup);
    blockers[0].estimatedWorkingSetBytes = limit / 2 + 1;
    mixed[0].estimatedWorkingSetBytes = limit;
    for (size_t i = 1; setup && i < mixed.size(); ++i) {
	setup = make_request(mixed[i], database, counting_complete_probe, healthy);
	mixed[i].estimatedWorkingSetBytes = 1;
    }
    auto blocker = setup ? coordinator.submit(blockers) :
	std::shared_ptr<BObolSourceRealizationJob>();
    if (!blocker || !wait_until([&]() {
	return blockerGate->entered.load(std::memory_order_acquire);
    }, std::chrono::seconds(2))) {
	blockerGate->release.store(true, std::memory_order_release);
	release_requests(blockers);
	release_requests(mixed);
	std::fprintf(stderr, "FAIL: cancelled frontier blocker setup\n");
	return 1;
    }
    auto job = coordinator.submit(mixed);
    cleanup.reset();
    int failures = 0;
    if (!job) {
	release_requests(mixed);
	std::fprintf(stderr, "FAIL: cancelled frontier queue setup\n");
	++failures;
    } else {
	/* As in the fairness regression, let bounded bypasses exhaust before
	 * removing their blocking root. Its destructor needs a newly enabled peer. */
	std::this_thread::sleep_for(std::chrono::milliseconds(100));
	if (!healthy->completed.load(std::memory_order_acquire) ||
	    healthy->completed.load(std::memory_order_acquire) >= fairness_small_request_count ||
	    coordinator.activeItemCountForDiagnostics() != 1) {
	    std::fprintf(stderr, "FAIL: cancelled frontier was not blocked by fairness\n");
	    ++failures;
	}
	mixed[0].stream->requestCancel();
	if (!healthyCompleted.load(std::memory_order_acquire) || blocker->isTerminal()) {
	    std::fprintf(stderr, "FAIL: queue cleanup delayed newly admissible work\n");
	    ++failures;
	}
    }
    blockerGate->release.store(true, std::memory_order_release);
    if (!wait_until([&]() {
	return blocker->isTerminal() && (!job || job->isTerminal());
    }, std::chrono::seconds(2))) {
	std::fprintf(stderr, "FAIL: cancelled frontier did not drain\n");
	++failures;
    }
    return failures;
}

static size_t
process_thread_count(void)
{
    return bu_file_list("/proc/self/task", "[0-9]*", NULL);
}

static int
test_pool_construction_failure(void)
{
    const size_t baselineThreads = process_thread_count();
    if (!baselineThreads) {
	std::fprintf(stderr, "FAIL: construction test requires Linux thread enumeration\n");
	return 1;
    }
    constexpr size_t retryAttempts = 3;
    for (size_t attempt = 0; attempt < retryAttempts; ++attempt) {
	bool rejected = false;
	{
	    ScopedTransactionFault fault(
		BObolTransactionFaultPoint::SOURCE_REALIZATION_WORKER_START);
	    try {
		(void)BObolSourceRealizationCoordinator::global();
	    } catch (const std::system_error &error) {
		rejected = error.code() == std::errc::resource_unavailable_try_again;
	    }
	}
	if (!rejected || !wait_until([&]() {
		return process_thread_count() == baselineThreads;
	    }, std::chrono::seconds(2))) {
	    std::fprintf(stderr,
		"FAIL: partial source pool retained threads (rejected=%d threads=%zu/%zu)\n",
		rejected, process_thread_count(), baselineThreads);
	    return 1;
	}
    }
    BObolSourceRealizationCoordinator &coordinator =
	BObolSourceRealizationCoordinator::global();
    if (coordinator.workerCountForDiagnostics() != pool_test_worker_count ||
	!wait_until([&]() {
	    return process_thread_count() == baselineThreads + pool_test_worker_count;
	}, std::chrono::seconds(2))) {
	std::fprintf(stderr, "FAIL: source pool could not recover after failed construction\n");
	return 1;
    }
    return 0;
}

struct ShutdownCancellationProbe {
    std::atomic<bool> entered{false};
    std::atomic<bool> *observedCancellation = nullptr;
};

struct ShutdownQueuedProbe {
    std::atomic<size_t> *completed = nullptr;
};

struct ShutdownPayloadProbe {
    BObolSourceRealizationCoordinator *coordinator = nullptr;
    std::atomic<bool> *released = nullptr;
    std::shared_ptr<BObolStagedSourceMesh> source;
    ~ShutdownPayloadProbe(void)
    {
	/* A storage owner may inspect service accounting during release. */
	(void)coordinator->activeWorkingSetBytesForDiagnostics();
	released->store(true, std::memory_order_release);
    }
};

static int
shutdown_cancellation_probe(SoBRLDatabaseSource *, struct db_i *, int, uint32_t,
    BObolCompactOccurrenceStream *stream, void *data)
{
    auto *probe = static_cast<ShutdownCancellationProbe *>(data);
    if (!probe || !stream)
	return 0;
    probe->entered.store(true, std::memory_order_release);
    const bool cancelled = wait_until([&]() { return stream->isCancelled(); },
	std::chrono::seconds(2));
    probe->observedCancellation->store(cancelled, std::memory_order_release);
    return 2;
}

static int
shutdown_queued_probe(SoBRLDatabaseSource *, struct db_i *, int, uint32_t,
    BObolCompactOccurrenceStream *, void *data)
{
    auto *probe = static_cast<ShutdownQueuedProbe *>(data);
    probe->completed->fetch_add(1, std::memory_order_acq_rel);
    return 2;
}

struct ShutdownObserver {
    std::shared_ptr<BObolSourceRealizationJob> active;
    std::shared_ptr<BObolSourceRealizationJob> queued;
    std::shared_ptr<BObolSourceRealizationJob> retiring;
    std::shared_ptr<BObolCompactOccurrenceStream> completedStream;
    std::thread cancellationThread;
    std::weak_ptr<ShutdownCancellationProbe> activeProbe;
    std::weak_ptr<ShutdownQueuedProbe> queuedProbe;
    std::atomic<bool> observedCancellation{false};
    std::atomic<size_t> queuedCompletions{0};
    std::atomic<bool> payloadReleased{false};
    std::atomic<bool> cleanupStarted{false};
    std::atomic<bool> cleanupFinished{false};
    std::vector<std::shared_ptr<ItemLifetimeEvidence>> sourceLifetimes;

    ~ShutdownObserver()
    {
	if (!active || !queued)
	    return;
	const bool retired = active->state() == BOBOL_SOURCE_REALIZATION_CANCELLED &&
	    queued->state() == BOBOL_SOURCE_REALIZATION_CANCELLED &&
	    observedCancellation.load(std::memory_order_acquire) &&
	    queuedCompletions.load(std::memory_order_acquire) == 0 &&
	    activeProbe.expired() && queuedProbe.expired() &&
	    std::all_of(sourceLifetimes.begin(), sourceLifetimes.end(),
		[](const std::shared_ptr<ItemLifetimeEvidence> &evidence) {
		    return evidence->sourceReleased.load(std::memory_order_acquire);
		});
	/* Check before joining: coordinator shutdown must be the lifetime
	 * barrier for a caller which already detached queued resources. */
	const bool callerRetired = retiring && retiring->isTerminal() &&
	    retiring->state() == BOBOL_SOURCE_REALIZATION_CANCELLED &&
	    cleanupFinished.load(std::memory_order_acquire);
	if (cancellationThread.joinable())
	    cancellationThread.join();
	const bool streamSurvived = completedStream && !completedStream->isCancelled();
	if (completedStream)
	    completedStream->requestCancel();
	active.reset();
	queued.reset();
	if (!retired || !callerRetired || !streamSurvived ||
	    !payloadReleased.load(std::memory_order_acquire) ||
	    !activeProbe.expired() || !queuedProbe.expired()) {
	    std::fprintf(stderr, "FAIL: source pool shutdown lost cancellation or callback retirement\n");
	    std::_Exit(EXIT_FAILURE);
	}
    }
};

struct ShutdownQueueCleanupProbe {
    BObolSourceRealizationCoordinator *coordinator = nullptr;
    ShutdownObserver *observer = nullptr;
    ~ShutdownQueueCleanupProbe()
    {
	(void)coordinator->activeWorkingSetBytesForDiagnostics();
	observer->cleanupStarted.store(true, std::memory_order_release);
	const bool stopping = wait_until([&]() {
	    return observer->observedCancellation.load(std::memory_order_acquire);
	}, std::chrono::seconds(2));
	/* Keep caller-owned cleanup active after the pool's short callbacks
	 * finish, as in the existing cache-registry teardown probe. */
	std::this_thread::sleep_for(std::chrono::milliseconds(100));
	observer->cleanupFinished.store(stopping, std::memory_order_release);
    }
};

static int
test_pool_shutdown(struct db_i *database)
{
    /* Retain interest across coordinator destruction; dropping local handles
     * before main returns would test client cancellation instead of shutdown.
     * The retained import probes cleanup during actual pool shutdown. */
    static ShutdownObserver observer;
    BObolSourceRealizationCoordinator &coordinator =
	BObolSourceRealizationCoordinator::global();
    {
	std::vector<BObolSourceRealizationRequest> completed(1);
	if (!make_request(completed[0], database, warm_complete_probe, {}))
	    return 1;
	auto job = coordinator.submit(completed);
	if (!job || !wait_until([&]() { return job->isTerminal(); }, std::chrono::seconds(2))) {
	    release_requests(completed);
	    std::fprintf(stderr, "FAIL: shutdown completed stream setup\n");
	    return 1;
	}
	observer.completedStream = completed[0].stream;
    }
    auto activeProbe = std::make_shared<ShutdownCancellationProbe>();
    auto queuedProbe = std::make_shared<ShutdownQueuedProbe>();
    activeProbe->observedCancellation = &observer.observedCancellation;
    queuedProbe->completed = &observer.queuedCompletions;
    std::vector<BObolSourceRealizationRequest> activeRequests(1);
    std::vector<BObolSourceRealizationRequest> queuedRequests(2);
    bool setup = make_request(activeRequests[0], database,
	shutdown_cancellation_probe, activeProbe);
    const size_t allowance = coordinator.workingSetLimitBytesForDiagnostics();
    activeRequests[0].estimatedWorkingSetBytes = allowance;
    for (BObolSourceRealizationRequest &request : queuedRequests) {
	if (!make_request(request, database, shutdown_queued_probe, queuedProbe))
	    setup = false;
	request.estimatedWorkingSetBytes = 1;
    }
    if (!setup || !allowance || allowance == SIZE_MAX) {
	release_requests(activeRequests);
	release_requests(queuedRequests);
	std::fprintf(stderr, "FAIL: source shutdown pressure setup\n");
	return 1;
    }
    observer.activeProbe = activeProbe;
    observer.queuedProbe = queuedProbe;
    for (auto *requests : {&activeRequests, &queuedRequests}) {
	for (auto &request : *requests) {
	    auto evidence = std::make_shared<ItemLifetimeEvidence>();
	    request.source->addChild(new ItemLifetimeNode(evidence));
	    observer.sourceLifetimes.push_back(evidence);
	}
    }
    {
	auto releaseProbe = std::make_shared<ShutdownPayloadProbe>();
	releaseProbe->coordinator = &coordinator;
	releaseProbe->released = &observer.payloadReleased;
	releaseProbe->source = make_retained_triangle();
	auto staged = std::make_shared<BObolStagedSourceMesh>(*releaseProbe->source);
	staged->owner = releaseProbe;
	if (!activeRequests[0].stream->retainStagedSource(staged)) {
	    release_requests(activeRequests);
	    release_requests(queuedRequests);
	    std::fprintf(stderr, "FAIL: shutdown payload setup\n");
	    return 1;
	}
    }
    observer.active = coordinator.submit(activeRequests);
    if (!observer.active || !wait_until([&]() {
	    return activeProbe->entered.load(std::memory_order_acquire);
	}, std::chrono::seconds(2))) {
	release_requests(activeRequests);
	release_requests(queuedRequests);
	std::fprintf(stderr, "FAIL: shutdown source did not start\n");
	return 1;
    }
    observer.queued = coordinator.submit(queuedRequests);
    if (!observer.queued || coordinator.activeItemCountForDiagnostics() != 1 ||
	coordinator.activeWorkingSetBytesForDiagnostics() != allowance ||
	coordinator.queuedItemCountForDiagnostics() != queuedRequests.size()) {
	release_requests(queuedRequests);
	std::fprintf(stderr, "FAIL: shutdown test did not establish active and queued work\n");
	return 1;
    }
    {
	auto probe = std::make_shared<ShutdownQueueCleanupProbe>();
	probe->coordinator = &coordinator;
	probe->observer = &observer;
	std::vector<BObolSourceRealizationRequest> requests(1);
	if (!make_request(requests[0], database, warm_complete_probe, probe))
	    return 1;
	requests[0].estimatedWorkingSetBytes = 1;
	observer.retiring = coordinator.submit(requests);
	if (!observer.retiring) {
	    release_requests(requests);
	    return 1;
	}
    }
    observer.cancellationThread = std::thread([job = observer.retiring]() { job->cancel(); });
    if (!wait_until([&]() { return observer.cleanupStarted.load(std::memory_order_acquire); },
	std::chrono::seconds(2))) {
	std::fprintf(stderr, "FAIL: queued cancellation did not enter caller-owned cleanup\n");
	return 1;
    }
    return 0;
}

} // namespace

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    const bool poolConstruction = argc == 2 &&
	bu_strcmp(argv[1], "--pool-construction") == 0;
    const bool poolShutdown = argc == 2 &&
	bu_strcmp(argv[1], "--pool-shutdown") == 0;
    if (argc != 1 && !poolConstruction && !poolShutdown) {
	std::fprintf(stderr, "Usage: %s [--pool-construction|--pool-shutdown]\n", argv[0]);
	return 1;
    }
    bobol_init(NULL);
    if (poolConstruction || poolShutdown)
	bu_setenv("BOBOL_SOURCE_REALIZATION_WORKERS",
	    std::to_string(pool_test_worker_count).c_str(), 1);
    if (poolConstruction)
	return test_pool_construction_failure();

    char path[MAXPATHLEN] = {0};
    FILE *temporary = bu_temp_file(path, sizeof(path));
    if (!temporary) {
	std::fprintf(stderr, "FAIL: could not create temporary database path\n");
	return 1;
    }
    std::fclose(temporary);
    struct db_i *database = db_create(path, 5);
    if (!database) {
	std::fprintf(stderr, "FAIL: could not create temporary database\n");
	(void)bu_file_delete(path);
	return 1;
    }

    {
	struct rt_wdb *wdbp = wdb_dbopen(database, RT_WDB_TYPE_DB_DISK);
	point_t center = VINIT_ZERO;
	fastf_t vertices[] = {0, 0, 0, 1, 0, 0, 0, 1, 0};
	int faces[] = {0, 1, 2};
	if (!wdbp || mk_sph(wdbp, "admission.s", center, 1.0) != 0 ||
	    mk_bot(wdbp, "admission.bot", RT_BOT_SURFACE, RT_BOT_UNORIENTED,
		0, 3, 1, vertices, faces, NULL, NULL) != 0) {
	    std::fprintf(stderr,
		"FAIL: could not create source-admission fixture\n");
	    db_close(database);
	    (void)bu_file_delete(path);
	    return 1;
	}
    }

    if (poolShutdown) {
	const int status = test_pool_shutdown(database);
	db_close(database);
	(void)bu_file_delete(path);
	return status;
    }

    int failures = 0;
    BObolSourceRealizationCoordinator &coordinator =
	BObolSourceRealizationCoordinator::global();
    if (coordinator.workerCountForDiagnostics() < 1) {
	std::fprintf(stderr, "FAIL: realization coordinator has no workers\n");
	failures++;
    }

    /* Construct the cache registry during normal execution.  Without the
     * coordinator/cache lifetime ordering contract its destructor would run
     * before an in-flight realization worker at process exit. */
    BObolDrawLodAssetRecord primingRecord;
    (void)bobol_draw_lod_asset_cache_get(database,
	"__source_realization_cache_priming__", &primingRecord);

    failures += test_submission_allocation_failure(coordinator, database,
	BObolTransactionFaultPoint::SOURCE_REALIZATION_QUEUE_APPEND, false);
    failures += test_submission_allocation_failure(coordinator, database,
	BObolTransactionFaultPoint::SOURCE_REALIZATION_JOB_HANDLE, false);
    failures += test_submission_allocation_failure(coordinator, database,
	BObolTransactionFaultPoint::SOURCE_REALIZATION_JOB_HANDLE, true);

    {
	/* Batch validation is transactional.  A malformed later item must not
	 * consume the source, database, stream, or callback lifetime belonging to
	 * an earlier valid request. */
	std::shared_ptr<CompletionCounter> context =
	    std::make_shared<CompletionCounter>();
	std::weak_ptr<CompletionCounter> weakContext = context;
	std::vector<BObolSourceRealizationRequest> requests(2);
	if (!make_request(requests[0], database, counting_complete_probe,
		context)) {
	    std::fprintf(stderr, "FAIL: rejected batch request setup\n");
	    failures++;
	} else {
	    SoBRLDatabaseSource *source = requests[0].source;
	    struct db_i *snapshot = requests[0].snapshotSourceDatabase;
	    std::shared_ptr<BObolCompactOccurrenceStream> stream =
		requests[0].stream;
	    std::shared_ptr<BObolSourceRealizationJob> job =
		coordinator.submit(requests);
	    if (job || requests[0].source != source ||
		requests[0].snapshotSourceDatabase != snapshot ||
		requests[0].stream != stream || weakContext.expired()) {
		std::fprintf(stderr,
		    "FAIL: rejected realization batch consumed ownership\n");
		failures++;
	    }
	    requests[0].source->unref();
	    requests[0].source = NULL;
	    db_close(requests[0].snapshotSourceDatabase);
	    requests[0].snapshotSourceDatabase = NULL;
	}
    }

    {
	std::vector<BObolSourceRealizationRequest> requests(1);
	if (!make_request(requests[0], database, warm_complete_probe,
		std::shared_ptr<void>())) {
	    std::fprintf(stderr, "FAIL: normal request setup\n");
	    failures++;
	} else {
	    std::shared_ptr<BObolSourceRealizationJob> job =
		coordinator.submit(requests);
	    if (!job || !wait_until([&job]() { return job->isTerminal(); },
		    std::chrono::seconds(2))) {
		std::fprintf(stderr, "FAIL: warm request did not complete\n");
		failures++;
	    } else {
		BObolSourceRealizationItemResult result;
		if (job->state() != BOBOL_SOURCE_REALIZATION_COMPLETE ||
		    !job->itemResult(0, result) || !result.source ||
		    !result.warmManifest) {
		    std::fprintf(stderr,
			"FAIL: warm request result contract\n");
		    failures++;
		}
	    }
	}
    }

    {
	/* A zero caller estimate must be resolved from the immutable leaf
	 * directory before worker admission.  This is deliberately a warm probe:
	 * the test observes the outer reservation without paying unrelated mesh
	 * realization cost. */
	std::shared_ptr<ProbeGate> gate = std::make_shared<ProbeGate>();
	std::vector<BObolSourceRealizationRequest> requests(1);
	if (!make_request(requests[0], database, blocking_warm_probe, gate,
		"admission.s")) {
	    std::fprintf(stderr, "FAIL: automatic admission request setup\n");
	    failures++;
	} else {
	    struct directory *dp = db_lookup(database, "admission.s",
		LOOKUP_QUIET);
	    const size_t expectedMinimum = dp && dp->d_len <=
		(SIZE_MAX - source_admission_leaf_fixed_bytes) /
		source_admission_leaf_copy_count ?
		dp->d_len * source_admission_leaf_copy_count +
		source_admission_leaf_fixed_bytes : 0;
	    std::shared_ptr<BObolSourceRealizationJob> job =
		coordinator.submit(requests);
	    if (!job || !expectedMinimum || !wait_until([&gate]() {
		    return gate->entered.load(std::memory_order_acquire);
		}, std::chrono::seconds(2)) ||
		coordinator.activeWorkingSetBytesForDiagnostics() <
		expectedMinimum) {
		std::fprintf(stderr,
		    "FAIL: automatic leaf admission did not reserve source bytes\n");
		failures++;
	    }
	    gate->release.store(true, std::memory_order_release);
	    if (job && !wait_until([&job]() { return job->isTerminal(); },
		std::chrono::seconds(2))) {
		std::fprintf(stderr,
		    "FAIL: automatic admission request did not complete\n");
		failures++;
	    }
	}
    }

    for (bool streamedLod : {false, true}) {
	/* Admission must agree with execution: only the streamed serialized path
	 * may replace the primitive-import reservation with a bounded census. */
	std::shared_ptr<ProbeGate> gate = std::make_shared<ProbeGate>();
	std::vector<BObolSourceRealizationRequest> requests(1);
	if (!make_request(requests[0], database, blocking_warm_probe, gate,
		"admission.bot")) {
	    std::fprintf(stderr, "FAIL: streamed BoT admission setup\n");
	    failures++;
	    continue;
	}
	requests[0].source->drawMode = SoBRLDatabaseSource::SHADED;
	requests[0].source->lodBotThreshold = streamedLod ? 1 : 0;
	struct directory *dp = db_lookup(database, "admission.bot", LOOKUP_QUIET);
	const size_t expectedBytes = streamedLod ?
	    coordinator.workingSetLimitBytesForDiagnostics() :
	    (dp ? dp->d_len * source_admission_leaf_copy_count +
	     source_admission_leaf_fixed_bytes : 0);
	std::shared_ptr<BObolSourceRealizationJob> job = coordinator.submit(requests);
	if (!job || !wait_until([&gate]() {
		return gate->entered.load(std::memory_order_acquire);
	    }, std::chrono::seconds(2)) ||
	    !coordinator.activeWorkingSetBytesForDiagnostics() ||
	    (expectedBytes != SIZE_MAX &&
	     coordinator.activeWorkingSetBytesForDiagnostics() != expectedBytes)) {
	    std::fprintf(stderr, "FAIL: streamed BoT source reservation (LoD=%d)\n",
		streamedLod ? 1 : 0);
	    failures++;
	}
	gate->release.store(true, std::memory_order_release);
	if (job && !wait_until([&job]() { return job->isTerminal(); },
		std::chrono::seconds(2))) {
	    std::fprintf(stderr, "FAIL: streamed BoT admission did not drain\n");
	    failures++;
	}
    }

    {
	std::vector<BObolSourceRealizationRequest> requests(1);
	if (!make_request(requests[0], database, NULL, std::shared_ptr<void>(),
		"admission.bot")) {
	    std::fprintf(stderr, "FAIL: serialized BoT coverage setup\n");
	    failures++;
	} else {
	    requests[0].source->drawMode = SoBRLDatabaseSource::SHADED;
	    requests[0].source->lodBotThreshold = 1;
	    std::shared_ptr<BObolSourceRealizationJob> job = coordinator.submit(requests);
	    BObolSourceRealizationItemResult result;
	    if (!job || !wait_until([&job]() { return job->isTerminal(); },
		    std::chrono::seconds(2)) ||
		job->state() != BOBOL_SOURCE_REALIZATION_COMPLETE ||
		!job->itemResult(0, result) || !result.stream ||
		result.stream->getExpectedCount() != 1) {
		if (job)
		    job->itemResult(0, result);
		std::fprintf(stderr,
		    "FAIL: serialized BoT coverage did not complete (state=%d, "
		    "count=%zu, diagnostic=%s)\n",
		    job ? job->state() : -1,
		    result.stream ? result.stream->getExpectedCount() : 0,
		    result.source ? result.source->realizationDiagnostic.getValue().getString() : "no source");
		failures++;
	    } else {
		std::vector<BObolCompactOccurrence> occurrences;
		result.stream->drain(occurrences, result.stream->size());
		const bool meshContract = std::any_of(occurrences.begin(), occurrences.end(),
		    [](const BObolCompactOccurrence &occurrence) {
			return occurrence.sourceMeshRequestValid &&
			    occurrence.summary.path == SbString("/admission.bot") &&
			    occurrence.sourceMeshRequest.faceCount == 1 &&
			    occurrence.sourceMeshRequest.pointCount == 3 &&
			    !occurrence.sourceMeshRequest.bounds.isEmpty();
		    });
		if (!meshContract) {
		    std::fprintf(stderr, "FAIL: bare BoT lost its lazy mesh contract\n");
		    for (const BObolCompactOccurrence &occurrence : occurrences)
			std::fprintf(stderr, "  path=%s request=%d faces=%zu points=%zu empty=%d\n",
			    occurrence.summary.path.getString(), occurrence.sourceMeshRequestValid,
			    occurrence.sourceMeshRequest.faceCount,
			    occurrence.sourceMeshRequest.pointCount,
			    occurrence.sourceMeshRequest.bounds.isEmpty());
		    failures++;
		}
	    }
	}
    }

    failures += test_individual_stream_cancellation(coordinator, database);
    failures += test_completed_stream_transfer(coordinator, database);
    failures += test_cancelled_frontier_wakeup(coordinator, database);
    for (bool memoryBlocked : {false, true}) {
	for (auto cancellation : {QueuedCancellation::Stream, QueuedCancellation::Job,
		QueuedCancellation::Interest, QueuedCancellation::Precancelled,
		QueuedCancellation::StreamlessJob, QueuedCancellation::ConcurrentStreams,
		QueuedCancellation::SharedStreamJob, QueuedCancellation::SharedStreamConstrained})
	    failures += test_queued_cancellation(coordinator, database, memoryBlocked, cancellation);
    }
    for (bool cold : {false, true}) {
	for (int outcome : {BOBOL_SOURCE_REALIZATION_COMPLETE,
		BOBOL_SOURCE_REALIZATION_CANCELLED, BOBOL_SOURCE_REALIZATION_FAILED})
	    failures += test_item_resource_retirement(coordinator, database, outcome, cold);
    }

    {
	std::shared_ptr<ProbeGate> gate = std::make_shared<ProbeGate>();
	std::weak_ptr<ProbeGate> weakGate = gate;
	std::vector<BObolSourceRealizationRequest> requests(1);
	if (!make_request(requests[0], database, blocking_warm_probe, gate)) {
	    std::fprintf(stderr, "FAIL: cancellation request setup\n");
	    failures++;
	} else {
	    std::shared_ptr<BObolSourceRealizationJob> job =
		coordinator.submit(requests);
	    /* The request is now the only durable owner.  Client teardown may
	     * cancel immediately, but callback storage must survive until the
	     * already-running callback returns. */
	    gate.reset();
	    if (!job || !wait_until([&weakGate]() {
		    std::shared_ptr<ProbeGate> current = weakGate.lock();
		    return current &&
			current->entered.load(std::memory_order_acquire);
		}, std::chrono::seconds(2))) {
		std::fprintf(stderr, "FAIL: blocking request did not start\n");
		failures++;
	    } else {
		const std::chrono::steady_clock::time_point started =
		    std::chrono::steady_clock::now();
		job.reset();
		const std::chrono::milliseconds elapsed =
		    std::chrono::duration_cast<std::chrono::milliseconds>(
			std::chrono::steady_clock::now() - started);
		if (elapsed > std::chrono::milliseconds(50)) {
		    std::fprintf(stderr,
			"FAIL: job teardown blocked for %lld ms\n",
			static_cast<long long>(elapsed.count()));
			failures++;
		}
	    }
	    std::shared_ptr<ProbeGate> retained = weakGate.lock();
	    if (!retained) {
		std::fprintf(stderr,
		    "FAIL: cancelled worker callback context was destroyed\n");
		failures++;
	    } else {
		retained->release.store(true, std::memory_order_release);
		retained.reset();
	    }
	    if (!wait_until([&coordinator]() {
		    return coordinator.activeItemCountForDiagnostics() == 0;
		}, std::chrono::seconds(2))) {
		std::fprintf(stderr,
		    "FAIL: cancelled realization was not reaped\n");
		failures++;
	    }
	    if (!wait_until([&weakGate]() { return weakGate.expired(); },
		    std::chrono::seconds(2))) {
		std::fprintf(stderr,
		    "FAIL: completed worker retained callback context\n");
		failures++;
	    }
	}
    }

    {
	const size_t limit =
	    coordinator.workingSetLimitBytesForDiagnostics();
	const size_t workers = coordinator.workerCountForDiagnostics();
	if (limit != SIZE_MAX && limit > 0 && workers > 1) {
	    std::shared_ptr<CountedProbeGate> gate =
		std::make_shared<CountedProbeGate>();
	    const size_t requestCount = std::min<size_t>(workers, 3);
	    std::vector<BObolSourceRealizationRequest> requests(requestCount);
	    bool setup = true;
	    for (BObolSourceRealizationRequest &request : requests) {
		if (!make_request(request, database,
			counted_blocking_warm_probe, gate)) {
		    setup = false;
		    break;
		}
		request.estimatedWorkingSetBytes = limit;
	    }
	    if (!setup) {
		std::fprintf(stderr, "FAIL: memory admission request setup\n");
		failures++;
	    } else {
		std::shared_ptr<BObolSourceRealizationJob> job =
		    coordinator.submit(requests);
		if (!job || !wait_until([&gate]() {
			return gate->entered.load(std::memory_order_acquire) > 0;
		    }, std::chrono::seconds(2))) {
		    std::fprintf(stderr,
			"FAIL: memory admission request did not start\n");
		    failures++;
		} else {
		    std::this_thread::sleep_for(
			std::chrono::milliseconds(50));
		    const size_t active =
			coordinator.activeItemCountForDiagnostics();
		    const size_t activeBytes =
			coordinator.activeWorkingSetBytesForDiagnostics();
		    if (gate->entered.load(std::memory_order_acquire) != 1 ||
			active != 1 || activeBytes > limit) {
			std::fprintf(stderr,
			    "FAIL: memory governor admitted %zu probes, "
			    "%zu active items, %zu/%zu bytes\n",
			    gate->entered.load(std::memory_order_acquire),
			    active, activeBytes, limit);
			failures++;
		    }
		}
		gate->release.store(true, std::memory_order_release);
		if (job && !wait_until([&job]() {
			return job->isTerminal();
		    }, std::chrono::seconds(2))) {
		    std::fprintf(stderr,
			"FAIL: governed realization did not complete\n");
		    failures++;
		}
	    }
	}
    }

    {
	/* A large old root may be bypassed briefly while a half-budget root is
	 * active, but an unlimited first-fit queue would let every later small
	 * root run first.  Bounded aging must stop that stream, drain the active
	 * reservation, and admit the large root next. */
	const size_t limit =
	    coordinator.workingSetLimitBytesForDiagnostics();
	const size_t workers = coordinator.workerCountForDiagnostics();
	if (limit != SIZE_MAX && limit > 32 && workers > 1) {
	    std::shared_ptr<ProbeGate> blocker =
		std::make_shared<ProbeGate>();
	    std::vector<BObolSourceRealizationRequest> blockerRequests(1);
	    if (!make_request(blockerRequests[0], database,
		    blocking_warm_probe, blocker)) {
		std::fprintf(stderr, "FAIL: fairness blocker setup\n");
		failures++;
	    } else {
		blockerRequests[0].estimatedWorkingSetBytes = limit / 2 + 1;
		std::shared_ptr<BObolSourceRealizationJob> blockerJob =
		    coordinator.submit(blockerRequests);
		if (!blockerJob || !wait_until([&blocker]() {
			return blocker->entered.load(std::memory_order_acquire);
		    }, std::chrono::seconds(2))) {
		    std::fprintf(stderr, "FAIL: fairness blocker did not start\n");
		    failures++;
		} else {
		    std::shared_ptr<ProbeGate> large =
			std::make_shared<ProbeGate>();
		    std::shared_ptr<CompletionCounter> small =
			std::make_shared<CompletionCounter>();
		    static const size_t smallCount = fairness_small_request_count;
		    std::vector<BObolSourceRealizationRequest> mixed(
			smallCount + 1);
		    bool setup = make_request(mixed[0], database,
			blocking_warm_probe, large);
		    mixed[0].estimatedWorkingSetBytes = limit;
		    for (size_t i = 1; setup && i < mixed.size(); ++i) {
			setup = make_request(mixed[i], database,
			    counting_complete_probe, small);
			mixed[i].estimatedWorkingSetBytes = 1;
		    }
		    std::shared_ptr<BObolSourceRealizationJob> mixedJob =
			setup ? coordinator.submit(mixed) :
			std::shared_ptr<BObolSourceRealizationJob>();
		    if (!mixedJob) {
			std::fprintf(stderr, "FAIL: fairness queue setup\n");
			failures++;
		    } else {
			std::this_thread::sleep_for(
			    std::chrono::milliseconds(100));
			if (large->entered.load(std::memory_order_acquire) ||
			    small->completed.load(std::memory_order_acquire) >=
				smallCount) {
			    std::fprintf(stderr,
				"FAIL: large root was not protected from "
				"first-fit starvation (%zu/%zu small)\n",
				small->completed.load(std::memory_order_acquire),
				smallCount);
			    failures++;
			}
			blocker->release.store(true, std::memory_order_release);
			if (!wait_until([&large]() {
				return large->entered.load(
				    std::memory_order_acquire);
			    }, std::chrono::seconds(2))) {
			    std::fprintf(stderr,
				"FAIL: aged large root was not admitted\n");
			    failures++;
			}
			large->release.store(true, std::memory_order_release);
			if (!wait_until([&mixedJob]() {
				return mixedJob->isTerminal();
			    }, std::chrono::seconds(2))) {
			    std::fprintf(stderr,
				"FAIL: fairness queue did not drain\n");
			    failures++;
			}
		    }
		}
		blocker->release.store(true, std::memory_order_release);
		if (blockerJob && !wait_until([&blockerJob]() {
			return blockerJob->isTerminal();
		    }, std::chrono::seconds(2))) {
		    std::fprintf(stderr,
			"FAIL: fairness blocker did not drain\n");
		    failures++;
		}
	    }
	}
    }

    {
	/* Admission must refuse an oversized detached import before any worker
	 * reaches its callback.  This is distinct from ordinary serialization:
	 * librt import allocation can otherwise terminate a display worker. */
	const size_t limit =
	    coordinator.workingSetLimitBytesForDiagnostics();
	if (limit != SIZE_MAX) {
	    std::shared_ptr<CompletionCounter> counter =
		std::make_shared<CompletionCounter>();
	    std::vector<BObolSourceRealizationRequest> requests(1);
	    if (!make_request(requests[0], database, counting_complete_probe,
		counter)) {
		std::fprintf(stderr, "FAIL: constrained admission setup\n");
		failures++;
	    } else {
		requests[0].estimatedWorkingSetBytes = limit + 1;
		auto evidence = std::make_shared<ItemLifetimeEvidence>();
		requests[0].source->addChild(new ItemLifetimeNode(evidence));
		std::weak_ptr<CompletionCounter> weakCounter = counter;
		std::shared_ptr<BObolSourceRealizationJob> job =
		    coordinator.submit(requests);
		BObolSourceRealizationItemResult result;
		if (!job || !job->isTerminal() ||
		    job->state() != BOBOL_SOURCE_REALIZATION_CONSTRAINED ||
		    counter->completed.load(std::memory_order_acquire) != 0 ||
		    !job->itemResult(0, result) ||
		    result.state != BOBOL_SOURCE_REALIZATION_CONSTRAINED ||
		    requests[0].source || requests[0].snapshotSourceDatabase ||
		    coordinator.activeItemCountForDiagnostics() != 0) {
		    std::fprintf(stderr,
			"FAIL: oversized source admission entered realization "
			"(job=%d terminal=%d item=%d callbacks=%zu source=%p db=%p active=%zu)\n",
			job ? job->state() : -1,
			job && job->isTerminal() ? 1 : 0,
			job && job->itemResult(0, result) ? result.state : -1,
			counter->completed.load(std::memory_order_acquire),
			static_cast<void *>(requests[0].source),
			static_cast<void *>(requests[0].snapshotSourceDatabase),
			coordinator.activeItemCountForDiagnostics());
		    failures++;
		}
		counter.reset();
		if (!weakCounter.expired() || result.source ||
		    !evidence->sourceReleased.load(std::memory_order_acquire)) {
		    std::fprintf(stderr, "FAIL: constrained item retained source or callback ownership\n");
		    failures++;
		}
	    }
	}
    }

    /* Deliberately leave one callback active across return from main.  Its
     * request owns an independent database handle, and the coordinator's
     * process-lifetime destructor must join it before cache static teardown.
     * This makes the shutdown race an ordinary and sanitizer regression. */
    {
	std::shared_ptr<ShutdownCacheProbe> probe =
	    std::make_shared<ShutdownCacheProbe>();
	std::vector<BObolSourceRealizationRequest> requests(1);
	if (!make_request(requests[0], database, shutdown_cache_probe, probe)) {
	    std::fprintf(stderr, "FAIL: shutdown cache request setup\n");
	    failures++;
	} else {
	    std::shared_ptr<BObolSourceRealizationJob> job =
		coordinator.submit(requests);
	    if (!job || !wait_until([&probe]() {
		    return probe->entered.load(std::memory_order_acquire);
		}, std::chrono::seconds(2))) {
		std::fprintf(stderr,
		    "FAIL: shutdown cache request did not start\n");
		failures++;
	    }
	}
    }

    db_close(database);
    (void)bu_file_delete(path);
    if (failures)
	return 1;
    std::printf("source realization coordinator contract passed\n");
    return 0;
}
