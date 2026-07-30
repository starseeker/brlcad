/*                 B L O D S E R V I C E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BLodService.h */

#ifndef BOBOL_BLODSERVICE_H
#define BOBOL_BLODSERVICE_H

#include "BObol/BLodRealization.h"
#include "BObol/BSourceMeshRequest.h"
#include "bv/view.h"

#include <Inventor/SbBasic.h>

#include <stddef.h>
#include <stdint.h>
#include <memory>
#include <vector>

class BObolLodService;

/* A view-required refinement no larger than these aggregate-admitted deltas
 * is cheaper to load and publish directly than to manufacture intermediate
 * per-level tasks.  Larger discontinuities (Lucy is the canonical example)
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
    struct db_i *dbip;
    std::shared_ptr<const BObolStagedSourceMesh> stagedSource;
    struct bv_view_info view;
    SbBool useView;
    SbBool refreshMissing;
    SbBool useForcedLevel;
    SbBool shrinkAfterCopy;
    SbBool compactResident;
    /* Pace a cold/warm retained asset toward the view-selected target instead
     * of making the first visible replacement an unbounded prefix.  These are
     * delivery budgets, not terminal-quality caps: subsequent rendered frames
     * continue until requestedLevel is resident. */
    SbBool progressiveDelivery;
    uint64_t initialRefinementFaceBudget;
    uint64_t initialRefinementPointBudget;
    double refinementGrowthFactor;
    SbBool useCurrentDrawLevel;
    int currentDrawLevel;
    /* Aggregate scene admission may allow less than the provider's local
     * per-asset growth rule.  Clamp this task's publication to the exact
     * level whose complete delta was reserved by the submit action. */
    SbBool useDeliveryLevelLimit;
    int deliveryLevelLimit;
    int forcedLevel;
    int reset;

    BObolMeshLodProvider(void);
    void clear(void);
};

struct BOBOL_EXPORT BObolRtSourceFullDetailProvider {
    struct db_i *dbip;
    SbBool validateSourceMetrics;
    uint64_t maxFullDetailFaceCount;
    uint64_t maxFullDetailPointCount;

    BObolRtSourceFullDetailProvider(void);
    void clear(void);
};

struct BOBOL_EXPORT BObolRtProxyProvider {
    struct db_i *dbip;
    int proxyKind;
    SbBool useRequestBounds;

    BObolRtProxyProvider(void);
    void clear(void);
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
 * each consuming the full host-derived allowance at the same time.  A request
 * larger than the limit is admitted only when it can run alone. */
BOBOL_EXPORT void
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

struct BObolLodServicePrivate;

class BOBOL_EXPORT BObolLodService {
public:
    BObolLodService(void);
    ~BObolLodService(void);

    SbBool start(size_t workerCount = 1, SbBool startCacheWriter = TRUE);
    void stop(void);
    SbBool isRunning(void) const;

    uint64_t beginGeneration(void);
    uint64_t currentGeneration(void) const;
    void cancelGeneration(uint64_t generation);
    SbBool isGenerationCancelled(uint64_t generation) const;

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
    size_t activeWorkingSetBytesForDiagnostics(void) const;
    size_t executingTaskCountForDiagnostics(void) const;
    size_t peakWorkingSetBytesForDiagnostics(void) const;
    size_t peakExecutingTaskCountForDiagnostics(void) const;
    size_t availableResultTaskCapacity(void) const;

    uint64_t submit(const BObolLodTask &task);
    uint64_t submitIfNotActive(const BObolLodTask &task);
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
    SbBool hasActiveRequest(const BObolLodRequest &request) const;
    /**
     * Drain a bounded presentation wave.  maxEstimatedBytes applies to mesh
     * payloads only and uses active point/index/normal populations; the first
     * result is always admitted so one exceptional mesh cannot starve.
     */
    size_t drainResults(std::vector<BObolLodResult> &results,
	size_t maxResults = 0, size_t maxEstimatedBytes = 0);
    size_t drainMatchingResults(std::vector<BObolLodResult> &results,
	const std::vector<BObolLodRequest> &requests,
	size_t maxResults = 0);
    BObolLodSubscriberId subscribeResultReady(
	BObolLodResultReadyCB callback, void *userData);
    void unsubscribeResultReady(BObolLodSubscriberId id);

    size_t inFlightCount(void) const;
    size_t pendingTaskCountForDiagnostics(void) const;
    size_t queuedResultCountForDiagnostics(void) const;
    size_t queuedCacheWriteCountForDiagnostics(void) const;
    size_t delayedTaskCountForDiagnostics(void) const;
    uint64_t rejectedTaskCountForDiagnostics(void) const;
    uint64_t coalescedResultCountForDiagnostics(void) const;
    uint64_t coalescedCacheWriteCountForDiagnostics(void) const;
    uint64_t discardedStaleResultCountForDiagnostics(void) const;
    size_t activeRequestCountForDiagnostics(void) const;
    size_t completedTaskCountForDiagnostics(void) const;
    size_t cancelledGenerationCountForDiagnostics(void) const;
    size_t residentMeshAssetCountForDiagnostics(void) const;
    size_t residentMeshBytesForDiagnostics(void) const;
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
    size_t queuedResidentMeshCompactionResultCountForDiagnostics(
	uint64_t consumerId) const;

    BObolLodResult realizeResidentMeshLod(
	const BObolLodRequest &request,
	const BObolMeshLodProvider &provider);
    /* Replace one view consumer's complete stable demand snapshot and queue
     * memory-bounded background trims against the aggregate demand.  This
     * call never reads or rewrites mesh arrays on its caller.  Assets absent
     * from all complete snapshots retain only their minimum useful prefix;
     * richer levels remain in the on-disk cache.  Returns the number of
     * newly queued assets.  planningComplete reports whether the bounded
     * resident-asset scan completed in this call; callers keep pumping quiet
     * work when it is FALSE. */
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
    void releaseResidentMeshConsumer(uint64_t consumerId);

private:
    BObolLodService(const BObolLodService &);
    BObolLodService &operator=(const BObolLodService &);

    BObolLodServicePrivate *p;
};

#endif /* BOBOL_BLODSERVICE_H */
