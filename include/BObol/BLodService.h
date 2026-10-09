/*                 B L O D S E R V I C E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BLodService.h */

#ifndef BOBOL_BLODSERVICE_H
#define BOBOL_BLODSERVICE_H

#include "BObol/BDefines.h"
#include "BObol/BLodRealization.h"
#include "BObol/BSourceMeshRequest.h"
#include "bv/view.h"

#include <Inventor/SbBasic.h>

#include <stddef.h>
#include <stdint.h>
#include <memory>
#include <vector>

class BObolLodService;
struct db_i;

/*
 * A background realization must not borrow the application's raw database
 * pointer.  This ref-counted lease holds one librt database use until every
 * copied provider and queued task has released it.
 */
class BOBOL_EXPORT BObolDatabaseLease {
public:
    static std::shared_ptr<BObolDatabaseLease> acquire(struct db_i *database);
    ~BObolDatabaseLease(void);

    struct db_i *get(void) const;

    BObolDatabaseLease(const BObolDatabaseLease &) = delete;
    BObolDatabaseLease &operator=(const BObolDatabaseLease &) = delete;

private:
    explicit BObolDatabaseLease(struct db_i *database);
    struct db_i *database;
};

/* A view-required refinement no larger than these aggregate-admitted deltas
 * is cheaper to load and publish directly than to manufacture intermediate
 * per-cut tasks.  Larger discontinuities (Lucy is the canonical example)
 * remain staged.  The submit action charges the complete selected delta to
 * its scene budget before the provider is allowed to make this jump. */
static constexpr uint64_t BOBOL_LOD_DIRECT_REFINEMENT_FACE_LIMIT = 16384;
static constexpr uint64_t BOBOL_LOD_DIRECT_REFINEMENT_POINT_LIMIT = 32768;

typedef BObolLodResult (*BObolLodTaskProc)(
	const BObolLodRequest &request, void *userData);
typedef void (*BObolLodTaskDataFreeProc)(void *userData);
typedef void (*BObolLodCacheWriteProc)(
	const BObolLodResult &result, void *userData);
typedef uint64_t BObolLodSubscriberId;
typedef void (*BObolLodResultReadyCB)(
	BObolLodService *service, void *userData);

struct BOBOL_EXPORT BObolMeshLodProvider {
    BObolLodService *service;
    uint64_t generation;
    std::shared_ptr<BObolDatabaseLease> databaseLease;
    std::shared_ptr<const BObolStagedSourceMesh> stagedSource;
    /* Cold BoTs which exceed practical whole-topology scratch limits use the
     * validated serialized V5 reader and bounded spatial pages.  This choice
     * is captured at task submission so an environment/configuration change
     * cannot alter a queued task's memory contract. */
    SbBool useSerializedSpatialSource;
    /* Exact persistent payload to reopen.  Request::sourceContentHash is a
     * general realization identity and must not be interpreted as an LMDB
     * key; authored sources and tests legitimately use non-cache hashes. */
    uint64_t meshAssetContentHash;
    /* Stable-view BREP representation refinement.  The standing source owns
     * the canonical band; a provider may materialize a finer immutable band
     * on a worker when it is absent from the persistent cache. */
    SbBool generateBrepVariant;
    double brepTessellationAbsTol;
    double brepTessellationRelTol;
    double brepTessellationNormTol;
    SbBool brepVariantMemoryLimited;
    SbBool refreshMissing;
    SbBool useForcedCut;
    SbBool shrinkAfterCopy;
    SbBool compactResident;
    /* Pace a cold/warm retained asset toward the view-selected target instead
     * of making the first visible replacement an unbounded prefix.  These are
     * delivery budgets, not terminal-quality caps: subsequent rendered frames
     * continue until requestedCut is resident. */
    SbBool progressiveDelivery;
    /* Render-cost allowance for the first cumulative prefix.  A single cost
     * currency avoids contradictory face/point ceilings and lets the worker
     * select the richest hierarchy cut which the scene allocator actually
     * admitted. */
    size_t initialRefinementCostBudget;
    double refinementGrowthFactor;
    SbBool useCurrentDrawCut;
    int currentDrawCut;
    /* Aggregate scene admission may allow less than the provider's local
     * per-asset growth rule.  Clamp this task's publication to the exact
     * cut whose complete delta was reserved by the submit action. */
    SbBool useDeliveryCutLimit;
    int deliveryCutLimit;
    /* The service's bounded transient allocation may be stricter than the
     * steady renderer/resident budget.  When submission lowers this task's
     * delivery cut to satisfy that bound, report a terminal constrained
     * presentation instead of repeatedly requesting an impossible prefix. */
    SbBool transientMemoryLimited;
    /* A zoom may need a richer resident prefix even when the measured frame
     * budget cannot present that complete prefix yet.  In that case the
     * delivery limit governs cache residency and this independent limit
     * governs the active draw cut published with the result.  It may lower an
     * already visible continuity cut only when presentationAdmissionCertified
     * proves that the scene allocator authorized the replacement. */
    SbBool usePresentationCutLimit;
    int presentationCutLimit;
    /* The owner-thread submit action has charged this task's visible cut to
     * one view/policy budget.  The shared mesh generation may invalidate the
     * allocation table before its result is drained; retain the admission
     * epoch separately so result adoption need not discard valid prepared
     * geometry and request the identical work again. */
    SbBool presentationAdmissionCertified;
    uint64_t presentationAdmissionViewRevision;
    uint64_t presentationAdmissionPolicyRevision;
    /* A finer adaptive representation must not replace a coherent standing
     * representation with its own minimum prefix.  Materialize the complete
     * requested prefix off the presentation path and publish it atomically.
     * This is deliberately distinct from cold first-content delivery: a cold
     * object has no useful mesh to preserve and should still expose its
     * minimum prefix as soon as possible. */
    SbBool atomicRepresentationHandoff;
    int forcedCut;
    SbBool resetExisting;

    BObolMeshLodProvider(void);
    void clear(void);
    SbBool setDatabase(struct db_i *database);
    struct db_i *getDatabase(void) const;
};

struct BOBOL_EXPORT BObolRtSourceFullDetailProvider {
    std::shared_ptr<BObolDatabaseLease> databaseLease;
    SbBool validateSourceMetrics;
    uint64_t maxFullDetailFaceCount;
    uint64_t maxFullDetailPointCount;

    BObolRtSourceFullDetailProvider(void);
    void clear(void);
    SbBool setDatabase(struct db_i *database);
    struct db_i *getDatabase(void) const;
};

struct BOBOL_EXPORT BObolRtProxyProvider {
    std::shared_ptr<BObolDatabaseLease> databaseLease;
    int proxyKind;
    SbBool useRequestBounds;

    BObolRtProxyProvider(void);
    void clear(void);
    SbBool setDatabase(struct db_i *database);
    struct db_i *getDatabase(void) const;
};

/* Host-visible latency intent for otherwise equivalent background work.
 * This is a dispatch hint only: it is not part of request/cache identity and
 * never bypasses task dependencies, transient-memory admission, or result
 * reservations.  Within one class the service retains its coarse-quality-
 * first, FIFO ordering. */
enum BObolLodTaskDispatchClass {
    BOBOL_LOD_TASK_DISPATCH_NORMAL = 0,
    BOBOL_LOD_TASK_DISPATCH_PREVIEW = 1
};

struct BOBOL_EXPORT BObolLodTask {
    uint64_t generation;
    BObolLodRequest request;
    std::vector<uint64_t> dependencies;
    BObolLodTaskProc realize;
    void *realizeData;
    BObolLodTaskDataFreeProc realizeDataFree;
    BObolLodCacheWriteProc cacheWrite;
    void *cacheWriteData;
    uint32_t debugDelayMilliseconds;
    /* Conservative peak bytes needed while realize() is executing.  Zero
     * asks the service to estimate from request.sourceCounts. */
    size_t estimatedWorkingSetBytes;
    /* Scheduling-only latency class.  A first-visible preview may pass queued
     * ordinary refinement, but cannot preempt executing work or weaken any
     * resource/dependency gate. */
    int dispatchClass;
    SbBool publishResult;
    SbBool writeCache;

    BObolLodTask(void);
    void clear(void);
    void addDependency(uint64_t taskId);
};

BOBOL_EXPORT BObolLodResult
bobol_mesh_lod_provider_task(const BObolLodRequest &request,
	void *userData);

BOBOL_EXPORT void
bobol_mesh_lod_provider_free(void *userData);

BOBOL_EXPORT BObolLodResult
bobol_mesh_lod_cache_provider_task(const BObolLodRequest &request,
	void *userData);

BOBOL_EXPORT BObolLodResult
bobol_rt_proxy_provider_task(const BObolLodRequest &request,
	void *userData);

BOBOL_EXPORT void
bobol_rt_proxy_provider_free(void *userData);

BOBOL_EXPORT BObolLodResult
bobol_rt_source_full_detail_provider_task(
	const BObolLodRequest &request, void *userData);

BOBOL_EXPORT void
bobol_rt_source_full_detail_provider_free(void *userData);

BOBOL_EXPORT SbBool
bobol_lod_rt_source_full_detail_request_from_source_mesh_request(
	BObolLodRequest &request,
	const BObolSourceMeshRequest &sourceRequest,
	const BObolLodRequest *templateRequest);

BOBOL_EXPORT uint64_t
bobol_lod_submit_rt_source_full_detail_request(
	BObolLodService *service,
	uint64_t generation,
	const BObolSourceMeshRequest &sourceRequest,
	struct db_i *dbip,
	const BObolLodRequest *templateRequest,
	uint64_t maxFullDetailFaceCount,
	uint64_t maxFullDetailPointCount);

/* Process-wide admission for transient LoD preparation memory.  Individual
 * services retain their own (possibly smaller) limits, while this shared
 * governor prevents several views or independent cold-preparation stages from
 * each consuming the full host-derived allowance at the same time.  FALSE
 * means the estimate can never fit the configured ceiling; callers must keep
 * their existing presentation and report a constrained result rather than
 * entering an allocation-heavy provider. */
BOBOL_EXPORT SbBool
bobol_lod_working_set_acquire(size_t estimatedBytes);

BOBOL_EXPORT void
bobol_lod_working_set_release(size_t estimatedBytes);

BOBOL_EXPORT size_t
bobol_lod_working_set_global_limit(void);

BOBOL_EXPORT size_t
bobol_lod_working_set_global_active_bytes(void);

BOBOL_EXPORT size_t
bobol_lod_working_set_global_peak_bytes(void);

BOBOL_EXPORT size_t
bobol_lod_working_set_global_active_tasks(void);

BOBOL_EXPORT size_t
bobol_lod_working_set_global_peak_tasks(void);

/* Select the bounded serialized spatial-page producer when a cold native PoP
 * build would exceed the service's transient working-set contract.  The
 * BOBOL_LOD_SPATIAL_LEAVES environment setting remains an explicit 0/1
 * diagnostic override; normal production selection is automatic. */
BOBOL_EXPORT SbBool
bobol_lod_spatial_source_enabled(const BObolLodRequest &request,
	size_t workingSetLimit);

BOBOL_EXPORT size_t
bobol_lod_spatial_task_working_set_bytes(void);

struct BObolLodServicePrivate;

/** One lock-consistent snapshot of transient service work.  Individual
 * counters remain available for phase diagnostics; callers deciding whether
 * a service is quiescent must use this value so new work classes cannot be
 * omitted from an independently reconstructed idle predicate. */
struct BOBOL_EXPORT BObolLodServiceWorkStatus {
    SbBool running = FALSE;
    SbBool stopping = FALSE;
    size_t pendingTasks = 0;
    size_t executingTasks = 0;
    /** Workers which have selected a task but are waiting for the shared
     * process-wide CPU budget.  These are included in executingTasks. */
    size_t cpuAdmissionWaitingTasks = 0;
    /** Workers which have acquired CPU admission but are waiting for the
     * process-wide transient-memory budget.  Ordinary tasks are included in
     * executingTasks; resident compactions are represented by the compaction
     * counters below. */
    size_t transientMemoryAdmissionWaitingTasks = 0;
    size_t inFlightTasks = 0;
    size_t resultReservations = 0;
    size_t cacheWriteReservations = 0;
    size_t activeRequests = 0;
    size_t queuedResults = 0;
    size_t queuedCacheWrites = 0;
    size_t delayedTasks = 0;
    size_t activeWorkingSetBytes = 0;
    size_t pendingResidentMeshCompactions = 0;
    size_t queuedResidentMeshCompactionResults = 0;
    size_t residentMeshCompactionResultReservations = 0;

    SbBool isIdle(void) const
    {
	return pendingTasks == 0 && executingTasks == 0 &&
	    cpuAdmissionWaitingTasks == 0 &&
	    transientMemoryAdmissionWaitingTasks == 0 &&
	    inFlightTasks == 0 && resultReservations == 0 &&
	    cacheWriteReservations == 0 && activeRequests == 0 &&
	    queuedResults == 0 && queuedCacheWrites == 0 &&
	    delayedTasks == 0 && activeWorkingSetBytes == 0 &&
	    pendingResidentMeshCompactions == 0 &&
	    queuedResidentMeshCompactionResults == 0 &&
	    residentMeshCompactionResultReservations == 0;
    }
};

/** One lock-consistent snapshot of work owned by a submission generation.
 * Shared-producer leases are explicit because a consumer generation can be
 * waiting for an immutable producer whose task belongs to another generation. */
struct BOBOL_EXPORT BObolLodGenerationWorkStatus {
    size_t activeTasks = 0;
    size_t pendingTasks = 0;
    size_t executingTasks = 0;
    /** Selected tasks waiting at the process-wide CPU gate.  This is a
     * subset of executingTasks, not an additional work obligation. */
    size_t cpuAdmissionWaitingTasks = 0;
    /** Selected tasks which passed the CPU gate and are waiting at the
     * process-wide transient-memory gate.  This is a subset of
     * executingTasks. */
    size_t transientMemoryAdmissionWaitingTasks = 0;
    size_t delayedTasks = 0;
    size_t queuedResults = 0;
    size_t queuedCacheWrites = 0;
    size_t sharedProducerLeases = 0;
    /* Exact stage of the earliest active mesh producer for this generation.
     * Multiple shared producers may occupy different stages; stageMask and
     * activeProducerCount expose that fact, while the task/unit fields
     * aggregate only producerStage.  A zero total is intentionally
     * indeterminate rather than an invented percentage. */
    uint32_t producerStageMask = 0;
    int producerStage = BOBOL_LOD_PRODUCER_STAGE_NONE;
    /* Number of relevant active producers in each exact stage.  Index zero
     * is the NONE sentinel and remains zero.  This preserves concurrent-stage
     * information which producerStage alone intentionally collapses to the
     * earliest stage. */
    size_t producerStageTaskCounts[BOBOL_LOD_PRODUCER_STAGE_COUNT] = {};
    size_t producerStageTaskCount = 0;
    size_t activeProducerCount = 0;
    uint64_t producerStageCompletedUnits = 0;
    uint64_t producerStageTotalUnits = 0;
    /** Age of the oldest task which has not yet acquired a worker.  This is
     * zero when no task owned by (or shared with) this generation is queued. */
    uint64_t oldestPendingTaskAgeMicroseconds = 0;
    /** Longest queue wait and execution lifetime among the active producers
     * represented by this snapshot. */
    uint64_t maximumProducerQueueWaitMicroseconds = 0;
    uint64_t maximumProducerElapsedMicroseconds = 0;
    /** Longest continuous residence in producerStage.  Together with the
     * unit counters this distinguishes slow work from a stalled queue. */
    uint64_t producerStageElapsedMicroseconds = 0;
    /** Source population represented by the active producer set.  Values are
     * saturating sums and are diagnostic only. */
    uint64_t activeProducerSourceFaceCount = 0;
    uint64_t activeProducerSourcePointCount = 0;
    uint64_t activeProducerSourceByteCount = 0;

    SbBool hasResultWork(void) const
    {
	return activeTasks != 0 || queuedResults != 0 ||
	    sharedProducerLeases != 0;
    }

    SbBool isIdle(void) const
    {
	return activeTasks == 0 && pendingTasks == 0 &&
	    executingTasks == 0 && delayedTasks == 0 &&
	    queuedResults == 0 && queuedCacheWrites == 0 &&
	    sharedProducerLeases == 0;
    }
};

/** Diagnostic classification of tasks which have not acquired a worker.
 * Counts are mutually exclusive and describe one generation at the instant
 * of the snapshot.  Submission-capacity flags are service-wide because a
 * full task or result reservation pool blocks every generation equally. */
struct BOBOL_EXPORT BObolLodQueueStatus {
    size_t runnableTasks = 0;
    size_t dependencyBlockedTasks = 0;
    size_t transientMemoryBlockedTasks = 0;
    SbBool taskCapacityBlocked = FALSE;
    SbBool resultCapacityBlocked = FALSE;

    size_t classifiedTaskCount(void) const
    {
	return runnableTasks + dependencyBlockedTasks +
	    transientMemoryBlockedTasks;
    }
};

/** Admission-safe resident-capacity observation.  Stable bytes are published
 * as one scalar before a completed growth releases its reservation.  A sample
 * may therefore conservatively count completed growth in both fields, but it
 * cannot expose that capacity to another producer twice. */
struct BOBOL_EXPORT BObolLodResidentCapacityStatus {
    size_t stableResidentBytes = 0;
    size_t reservedGrowthBytes = 0;
    size_t residentLimitBytes = 0;

    size_t occupiedBytes(void) const
    {
	return reservedGrowthBytes > SIZE_MAX - stableResidentBytes ?
	    SIZE_MAX : stableResidentBytes + reservedGrowthBytes;
    }

    SbBool exceedsLimit(void) const
    {
	return residentLimitBytes != SIZE_MAX &&
	    occupiedBytes() > residentLimitBytes;
    }

    SbBool fitsWithinLimitFraction(size_t divisor) const
    {
	if (divisor == 0)
	    return FALSE;
	if (residentLimitBytes == SIZE_MAX)
	    return TRUE;
	return occupiedBytes() <= residentLimitBytes / divisor;
    }
};

/**
 * Process-local bounded LoD execution and retained-residency service.
 *
 * Concurrency contract:
 *
 * - Public methods are thread safe unless their documentation explicitly
 *   names the presentation-owner thread.
 * - The service queue/generation/subscriber lock precedes an individual
 *   resident-asset lock whenever both are required.  Callers and providers
 *   must never invert that order.
 * - Provider realization, cache I/O, mesh preparation, result-ready
 *   callbacks, and subscriber callbacks execute with neither service nor
 *   resident-asset locks held.
 * - Task-local payloads are worker owned until immutable publication.
 *   Completed mesh arrays are shared immutable values; Coin nodes and fields
 *   remain presentation-owner-thread only.
 * - Result-ready callbacks may run on a worker.  They are edge notifications,
 *   not permission to mutate a view; clients must schedule their owner-thread
 *   result pump.
 * - stop() joins service workers and is the explicit process/service shutdown
 *   barrier.  View and endpoint teardown should instead unsubscribe or release
 *   its consumer interest, allowing unrelated shared work to continue.
 */
class BOBOL_EXPORT BObolLodService {
public:
    BObolLodService(void);
    ~BObolLodService(void);

    SbBool start(size_t workerCount = 1, SbBool startCacheWriter = TRUE);
    /* A shared service may acquire additional view clients over time.  Grow
     * the pool without cancelling resident assets or outstanding work. */
    SbBool ensureWorkerCount(size_t workerCount);
    void stop(void);
    SbBool isRunning(void) const;
    size_t workerCountForDiagnostics(void) const;

    uint64_t beginGeneration(void);
    uint64_t currentGeneration(void) const;
    void cancelGeneration(uint64_t generation);
    SbBool isGenerationCancelled(uint64_t generation) const;
    /* Provider-side cancellation for a potentially shared asset producer.
     * A cancelled submission generation does not stop immutable production
     * while another generation retains a lease on the same active key. */
    SbBool isProducerCancelled(uint64_t generation,
	const BObolLodRequest &request) const;

    void setQueueLimits(size_t maxActiveTasks,
	size_t maxQueuedResults, size_t maxQueuedCacheWrites);
    void getQueueLimits(size_t &maxActiveTasks,
	size_t &maxQueuedResults, size_t &maxQueuedCacheWrites) const;
    /* Bound concurrent realization scratch independently of worker count.
     * Passing zero restores the host-derived default; SIZE_MAX disables the
     * byte governor.  One task larger than the limit may run only when no
     * other realization is active, preventing deadlock while serializing
     * exceptional meshes. */
    void setWorkingSetLimit(size_t maxActiveBytes);
    size_t getWorkingSetLimit(void) const;
    /* Bound completed service-owned progressive mesh generations separately
     * from transient worker scratch.  Passing zero restores the host-derived
     * default; SIZE_MAX disables pressure eviction.  Stable compaction first
     * trims demanded assets to their view cuts, then releases wholly
     * undemanded assets when retained bytes exceed this limit. */
    void setResidentMeshLimit(size_t maxResidentBytes);
    size_t getResidentMeshLimit(void) const;
    /* Set the process-shared resident ceiling to a percentage of currently
     * available system memory.  The percentage is deliberately capped: the
     * retained CPU prefix is not the only copy of drawing data, and database
     * mappings, renderer records, GPU/OSMesa storage, and transient topology
     * still need headroom.  This is a snapshot-based explicit limit; passing
     * zero to setResidentMeshLimit() restores the automatic policy. */
    SbBool setResidentMeshAvailableMemoryPercent(double availableMemoryPercent);
    double getResidentMeshAvailableMemoryPercent(void) const;
    size_t getResidentMeshAvailableMemoryBasisBytes(void) const;
    static double getMaximumResidentMeshAvailableMemoryPercent(void);
    size_t activeWorkingSetBytesForDiagnostics(void) const;
    size_t executingTaskCountForDiagnostics(void) const;
    size_t peakWorkingSetBytesForDiagnostics(void) const;
    size_t peakExecutingTaskCountForDiagnostics(void) const;
    size_t availableResultTaskCapacity(void) const;

    uint64_t submit(const BObolLodTask &task);
    /* Suppress an equivalent request while either its producer is active or
     * its completed result is still awaiting publication. */
    uint64_t submitIfNotActive(const BObolLodTask &task);
    /* Offer an immutable intermediate result from the active provider for
     * generation.  This operation never waits for the service lock:
     * providers may call it while holding a resident-asset lock without
     * violating the documented lock order.  TRUE means the newest preview
     * was accepted for immediate or deferred publication; FALSE means it was
     * safely skipped.  Deferred values coalesce per producer, and the
     * ordinary task completion remains authoritative. */
    SbBool tryPublishIntermediateResult(
	uint64_t generation, BObolLodResult &&result);
    /*
     * Submit a producer wave under one service lock and wake workers once.
     * taskIds is resized to tasks.size(); rejected/duplicate entries are
     * reported as zero while accepted entries contain their task id.
     * Returns the number accepted.  The caller retains ownership of
     * realizeData for every zero-id entry.
     */
    size_t submitBatch(const std::vector<BObolLodTask> &tasks,
	std::vector<uint64_t> &taskIds,
	SbBool skipActiveDuplicates = FALSE);
    /** Record @p request as the newest demand for an equivalent active asset
     * producer.  Returns TRUE while either that producer is active or its
     * completed result is queued.  A cold producer uses the retained demand
     * to stamp immutable page previews with the current camera/policy epoch;
     * it never changes the producer's geometry identity or working-set
     * reservation. */
    SbBool updateActiveRequestDemand(const BObolLodRequest &request,
	uint64_t generation = 0);
    SbBool hasActiveRequest(const BObolLodRequest &request) const;
    /**
     * Drain a bounded presentation wave.  maxEstimatedBytes applies to mesh
     * payloads only and uses active point/index/normal populations; the first
     * result is always admitted so one exceptional mesh cannot starve.
     */
    size_t drainResults(std::vector<BObolLodResult> &results,
	size_t maxResults = 0, size_t maxEstimatedBytes = 0);
    /*
     * Drain only results owned by one submission generation.  A service may
     * be shared by several views; presentation consumers must use this entry
     * point so one view cannot consume another view's completed payloads.
     */
    size_t drainGenerationResults(std::vector<BObolLodResult> &results,
	uint64_t generation, size_t maxResults = 0,
	size_t maxEstimatedBytes = 0);
    size_t drainMatchingResults(std::vector<BObolLodResult> &results,
	const std::vector<BObolLodRequest> &requests,
	size_t maxResults = 0);
    BObolLodSubscriberId subscribeResultReady(
	BObolLodResultReadyCB callback, void *userData);
    void unsubscribeResultReady(BObolLodSubscriberId id);

    size_t inFlightCount(void) const;
    BObolLodServiceWorkStatus workStatus(void) const;
    size_t resultReservationCountForDiagnostics(void) const;
    size_t pendingTaskCountForDiagnostics(void) const;
    size_t queuedResultCountForDiagnostics(void) const;
    size_t queuedCacheWriteCountForDiagnostics(void) const;
    size_t delayedTaskCountForDiagnostics(void) const;
    /* O(1) per-generation diagnostics for shared-service consumers. */
    BObolLodGenerationWorkStatus generationWorkStatus(
	uint64_t generation) const;
    BObolLodQueueStatus generationQueueStatusForDiagnostics(
	uint64_t generation) const;
    size_t activeTaskCountForGeneration(uint64_t generation) const;
    size_t pendingTaskCountForGeneration(uint64_t generation) const;
    size_t executingTaskCountForGeneration(uint64_t generation) const;
    size_t queuedResultCountForGeneration(uint64_t generation) const;
    size_t queuedCacheWriteCountForGeneration(uint64_t generation) const;
    size_t delayedTaskCountForGeneration(uint64_t generation) const;
    uint64_t rejectedTaskCountForDiagnostics(void) const;
    uint64_t coalescedResultCountForDiagnostics(void) const;
    uint64_t coalescedCacheWriteCountForDiagnostics(void) const;
    uint64_t discardedStaleResultCountForDiagnostics(void) const;
    size_t activeRequestCountForDiagnostics(void) const;
    size_t sharedProducerCountForDiagnostics(void) const;
    size_t sharedProducerLeaseCountForDiagnostics(void) const;
    size_t sharedProducerLeaseCountForGeneration(
	uint64_t generation) const;
    size_t completedTaskCountForDiagnostics(void) const;
    size_t cancelledGenerationCountForDiagnostics(void) const;
    size_t residentMeshAssetCountForDiagnostics(void) const;
    size_t residentMeshBytesForDiagnostics(void) const;
    BObolLodResidentCapacityStatus residentCapacityStatus(void) const;
    /* Stable retained bytes exclude the reloadable cache-reader prefix which
     * is released after a quiet interval.  Optional suffix admission is
     * governed by this value; transient preparation has its own byte
     * governor. */
    size_t stableResidentMeshBytesForDiagnostics(void) const;
    size_t reservedResidentMeshGrowthBytesForDiagnostics(void) const;
    /* Advances only when stable resident capacity becomes newly available
     * (reclamation or a relaxed limit), not when another asset merely grows.
     * Memory-limited view bindings use it to suppress retry loops. */
    uint64_t residentMeshAdmissionRevision(void) const;
    uint64_t residentMeshCacheLoadCountForDiagnostics(void) const;
    uint64_t residentMeshHitCountForDiagnostics(void) const;
    uint64_t residentMeshCompactionCountForDiagnostics(void) const;
    uint64_t residentMeshEvictionCountForDiagnostics(void) const;
    size_t pendingResidentMeshCompactionCountForDiagnostics(void) const;
    /* Report the last complete demand-aware scan for one view.  TRUE means
     * the candidate count describes the current demand and resident
     * generation.  A zero count then certifies that retained bytes are
     * required by the visible working set (or its reloadable minimum), not
     * an unserviced reclamation request. */
    SbBool residentMeshCompactionPlanForDiagnostics(
	uint64_t consumerId, uint64_t *demandRevision,
	size_t *candidateCount) const;
    size_t queuedResidentMeshCompactionResultCountForDiagnostics(
	uint64_t consumerId) const;

    BObolLodResult realizeResidentMeshLod(
	const BObolLodRequest &request,
	const BObolMeshLodProvider &provider);
    /* Replace one view consumer's complete stable demand snapshot and queue
     * memory-bounded background trims against the aggregate demand.  This
     * call never reads or rewrites mesh arrays on its caller.  Assets absent
     * from current snapshots retain their richer prefixes while capacity is
     * available, making resident geometry a latency cache across view
     * changes.  Once the configured limit is exceeded, demanded assets may
     * be trimmed and undemanded assets evicted; richer cuts remain in the
     * on-disk cache.  Returns the number of newly queued assets.
     * planningComplete reports whether the bounded resident-asset scan
     * completed in this call; callers keep pumping quiet work when it is
     * FALSE. */
    size_t scheduleResidentMeshCompaction(
	uint64_t consumerId,
	uint64_t demandRevision,
	const std::vector<BObolLodResidentDemand> &demands,
	SbBool *planningComplete = NULL);
    /* Drain immutable, renderer-ready generations completed for one consumer.
     * The presentation owner adopts these handles directly; no vertex or
     * index arrays are copied on the GUI thread. */
    size_t drainResidentMeshCompactions(
	uint64_t consumerId,
	std::vector<BObolLodResidentCompaction> &results,
	size_t maxResults = 0);
    /* Retire one consumer's complete demand snapshot immediately.  This is
     * O(1) and invalidates every in-flight aggregate trim before its shared
     * mesh generation can commit.  The next quiet compaction pass installs a
     * replacement complete snapshot. */
    void invalidateResidentMeshConsumer(uint64_t consumerId);
    /* Invalidate an already planned quiet trim before promoting a drawable
     * prefix which is resident in memory.  Provider requests perform the same
     * use notification internally; this form covers zero-I/O retargets. */
    void noteResidentMeshUse(const BObolLodCacheKey &assetKey);
    void releaseResidentMeshConsumer(uint64_t consumerId);

private:
    BObolLodService(const BObolLodService &);
    BObolLodService &operator=(const BObolLodService &);

    BObolLodServicePrivate *p;
};

#endif /* BOBOL_BLODSERVICE_H */
