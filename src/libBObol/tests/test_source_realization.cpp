/*          T E S T _ S O U R C E _ R E A L I Z A T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BInit.h"
#include "BObol/BSourceRealization.h"
#include "bu/app.h"
#include "bu/file.h"
#include "rt/db_io.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <functional>
#include <memory>
#include <thread>
#include <vector>

namespace {

struct ProbeGate {
    std::atomic<bool> entered{false};
    std::atomic<bool> release{false};
};

struct CountedProbeGate {
    std::atomic<size_t> entered{0};
    std::atomic<bool> release{false};
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
    struct db_i *database, BObolSourceWarmManifestProbe probe, void *data)
{
    if (!database)
	return false;
    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->ref();
    struct db_i *snapshot = db_clone_dbi(database, NULL);
    if (!snapshot) {
	source->unref();
	return false;
    }
    request.source = source;
    request.snapshotSourceDatabase = snapshot;
    request.stream = std::make_shared<BObolCompactOccurrenceStream>();
    request.probeWarmManifest = probe;
    request.callbackData = data;
    return true;
}

} // namespace

int
main(int argc, char **argv)
{
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

    int failures = 0;
    BObolSourceRealizationCoordinator &coordinator =
	BObolSourceRealizationCoordinator::global();
    if (coordinator.workerCountForDiagnostics() < 1) {
	std::fprintf(stderr, "FAIL: realization coordinator has no workers\n");
	failures++;
    }

    {
	std::vector<BObolSourceRealizationRequest> requests(1);
	if (!make_request(requests[0], database, warm_complete_probe, NULL)) {
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
	ProbeGate gate;
	std::vector<BObolSourceRealizationRequest> requests(1);
	if (!make_request(requests[0], database, blocking_warm_probe, &gate)) {
	    std::fprintf(stderr, "FAIL: cancellation request setup\n");
	    failures++;
	} else {
	    std::shared_ptr<BObolSourceRealizationJob> job =
		coordinator.submit(requests);
	    if (!job || !wait_until([&gate]() {
		    return gate.entered.load(std::memory_order_acquire);
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
	    gate.release.store(true, std::memory_order_release);
	    if (!wait_until([&coordinator]() {
		    return coordinator.activeItemCountForDiagnostics() == 0;
		}, std::chrono::seconds(2))) {
		std::fprintf(stderr,
		    "FAIL: cancelled realization was not reaped\n");
		failures++;
	    }
	}
    }

    {
	const size_t limit =
	    coordinator.workingSetLimitBytesForDiagnostics();
	const size_t workers = coordinator.workerCountForDiagnostics();
	if (limit != SIZE_MAX && limit > 0 && workers > 1) {
	    CountedProbeGate gate;
	    const size_t requestCount = std::min<size_t>(workers, 3);
	    std::vector<BObolSourceRealizationRequest> requests(requestCount);
	    bool setup = true;
	    for (BObolSourceRealizationRequest &request : requests) {
		if (!make_request(request, database,
			counted_blocking_warm_probe, &gate)) {
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
			return gate.entered.load(std::memory_order_acquire) > 0;
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
		    if (gate.entered.load(std::memory_order_acquire) != 1 ||
			active != 1 || activeBytes > limit) {
			std::fprintf(stderr,
			    "FAIL: memory governor admitted %zu probes, "
			    "%zu active items, %zu/%zu bytes\n",
			    gate.entered.load(std::memory_order_acquire),
			    active, activeBytes, limit);
			failures++;
		    }
		}
		gate.release.store(true, std::memory_order_release);
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

    db_close(database);
    (void)bu_file_delete(path);
    if (failures)
	return 1;
    std::printf("source realization coordinator contract passed\n");
    return 0;
}
