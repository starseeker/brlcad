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
#include <vector>

class BObolLodService;

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
    struct bv_view_info view;
    SbBool useView;
    SbBool refreshMissing;
    SbBool useForcedLevel;
    SbBool shrinkAfterCopy;
    SbBool compactResident;
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
    size_t availableResultTaskCapacity(void) const;

    uint64_t submit(const BObolLodTask &task);
    uint64_t submitIfNotActive(const BObolLodTask &task);
    SbBool hasActiveRequest(const BObolLodRequest &request) const;
    size_t drainResults(std::vector<BObolLodResult> &results,
	size_t maxResults = 0);
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
    uint64_t residentMeshCacheLoadCountForDiagnostics(void) const;
    uint64_t residentMeshHitCountForDiagnostics(void) const;
    uint64_t residentMeshCompactionCountForDiagnostics(void) const;

    BObolLodResult realizeResidentMeshLod(
	const BObolLodRequest &request,
	const BObolMeshLodProvider &provider);
    /* Replace one view consumer's complete stable demand snapshot, aggregate
     * it with every other consumer, and trim only assets whose resident prefix
     * exceeds the aggregate maximum plus headroom.  Returns the number of
     * assets compacted. */
    size_t compactResidentMeshes(
	uint64_t consumerId,
	uint64_t demandRevision,
	const std::vector<BObolLodResidentDemand> &demands,
	int headroomLevels = 1);
    void releaseResidentMeshConsumer(uint64_t consumerId);

private:
    BObolLodService(const BObolLodService &);
    BObolLodService &operator=(const BObolLodService &);

    BObolLodServicePrivate *p;
};

#endif /* BOBOL_BLODSERVICE_H */
