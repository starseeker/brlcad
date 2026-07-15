/*                  L O D _ S E R V I C E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/lod_service.h */

#ifndef BRLOBOL_LOD_SERVICE_H
#define BRLOBOL_LOD_SERVICE_H

#include "brlobol/lod_realization.h"
#include "brlobol/source_mesh_request.h"
#include "bv/view.h"

#include <Inventor/SbBasic.h>

#include <stddef.h>
#include <stdint.h>
#include <vector>

class BRLObolLodService;

typedef BRLObolLodResult (*BRLObolLodTaskProc)(
	const BRLObolLodRequest &request, void *userData);
typedef void (*BRLObolLodTaskDataFreeProc)(void *userData);
typedef void (*BRLObolLodCacheWriteProc)(
	const BRLObolLodResult &result, void *userData);
typedef uint64_t BRLObolLodSubscriberId;
typedef void (*BRLObolLodResultReadyCB)(
	BRLObolLodService *service, void *userData);

struct BRLOBOL_EXPORT BRLObolMeshLodProvider {
    struct db_i *dbip;
    struct bv_view_info view;
    SbBool useView;
    SbBool refreshMissing;
    SbBool useForcedLevel;
    SbBool shrinkAfterCopy;
    int forcedLevel;
    int reset;

    BRLObolMeshLodProvider(void);
    void clear(void);
};

struct BRLOBOL_EXPORT BRLObolRtSourceFullDetailProvider {
    struct db_i *dbip;
    SbBool validateSourceMetrics;
    uint64_t maxFullDetailFaceCount;
    uint64_t maxFullDetailPointCount;

    BRLObolRtSourceFullDetailProvider(void);
    void clear(void);
};

struct BRLOBOL_EXPORT BRLObolRtProxyProvider {
    struct db_i *dbip;
    int proxyKind;
    SbBool useRequestBounds;

    BRLObolRtProxyProvider(void);
    void clear(void);
};

struct BRLOBOL_EXPORT BRLObolLodTask {
    uint64_t generation;
    BRLObolLodRequest request;
    std::vector<uint64_t> dependencies;
    BRLObolLodTaskProc realize;
    void *realizeData;
    BRLObolLodTaskDataFreeProc realizeDataFree;
    BRLObolLodCacheWriteProc cacheWrite;
    void *cacheWriteData;
    uint32_t debugDelayMilliseconds;
    SbBool publishResult;
    SbBool writeCache;

    BRLObolLodTask(void);
    void clear(void);
    void addDependency(uint64_t taskId);
};

BRLOBOL_EXPORT BRLObolLodResult
brlobol_mesh_lod_provider_task(const BRLObolLodRequest &request,
	void *userData);

BRLOBOL_EXPORT void
brlobol_mesh_lod_provider_free(void *userData);

BRLOBOL_EXPORT BRLObolLodResult
brlobol_mesh_lod_cache_provider_task(const BRLObolLodRequest &request,
	void *userData);

BRLOBOL_EXPORT BRLObolLodResult
brlobol_rt_proxy_provider_task(const BRLObolLodRequest &request,
	void *userData);

BRLOBOL_EXPORT void
brlobol_rt_proxy_provider_free(void *userData);

BRLOBOL_EXPORT BRLObolLodResult
brlobol_rt_source_full_detail_provider_task(
	const BRLObolLodRequest &request, void *userData);

BRLOBOL_EXPORT void
brlobol_rt_source_full_detail_provider_free(void *userData);

BRLOBOL_EXPORT SbBool
brlobol_lod_rt_source_full_detail_request_from_source_mesh_request(
	BRLObolLodRequest &request,
	const BRLObolSourceMeshRequest &sourceRequest,
	const BRLObolLodRequest *templateRequest);

BRLOBOL_EXPORT uint64_t
brlobol_lod_submit_rt_source_full_detail_request(
	BRLObolLodService *service,
	uint64_t generation,
	const BRLObolSourceMeshRequest &sourceRequest,
	struct db_i *dbip,
	const BRLObolLodRequest *templateRequest,
	uint64_t maxFullDetailFaceCount,
	uint64_t maxFullDetailPointCount);

struct BRLObolLodServicePrivate;

class BRLOBOL_EXPORT BRLObolLodService {
public:
    BRLObolLodService(void);
    ~BRLObolLodService(void);

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

    uint64_t submit(const BRLObolLodTask &task);
    uint64_t submitIfNotActive(const BRLObolLodTask &task);
    SbBool hasActiveRequest(const BRLObolLodRequest &request) const;
    size_t drainResults(std::vector<BRLObolLodResult> &results,
	size_t maxResults = 0);
    size_t drainMatchingResults(std::vector<BRLObolLodResult> &results,
	const std::vector<BRLObolLodRequest> &requests,
	size_t maxResults = 0);
    BRLObolLodSubscriberId subscribeResultReady(
	BRLObolLodResultReadyCB callback, void *userData);
    void unsubscribeResultReady(BRLObolLodSubscriberId id);

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

private:
    BRLObolLodService(const BRLObolLodService &);
    BRLObolLodService &operator=(const BRLObolLodService &);

    BRLObolLodServicePrivate *p;
};

#endif /* BRLOBOL_LOD_SERVICE_H */
