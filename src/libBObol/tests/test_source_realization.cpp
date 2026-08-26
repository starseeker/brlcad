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
#include "rt/db_io.h"
#include "wdb.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <functional>
#include <memory>
#include <thread>
#include <vector>

namespace {

static const size_t source_admission_leaf_copy_count = 3;
static const size_t source_admission_leaf_fixed_bytes =
    8ULL * 1024ULL * 1024ULL;

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
    struct db_i *snapshot = db_clone_dbi(database, NULL);
    if (!snapshot) {
	source->unref();
	return false;
    }
    request.source = source;
    request.snapshotSourceDatabase = snapshot;
    request.stream = std::make_shared<BObolCompactOccurrenceStream>();
    request.probeWarmManifest = probe;
    request.callbackContext = context;
    return true;
}

} // namespace

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    if (argc != 1) {
	std::fprintf(stderr, "Usage: %s\n", argv[0]);
	return 1;
    }
    bobol_init(NULL);

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
	if (!wdbp || mk_sph(wdbp, "admission.s", center, 1.0) != 0) {
	    std::fprintf(stderr,
		"FAIL: could not create source-admission fixture\n");
	    db_close(database);
	    (void)bu_file_delete(path);
	    return 1;
	}
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
		    static const size_t smallCount = 24;
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
