/*          B S O U R C E R E A L I Z A T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BSourceRealization.h
 *
 * Process-wide, cancellable realization of detached CAD database sources.
 *
 * The coordinator owns worker execution and the resources transferred in a
 * request. A client-owned job is only an interest handle: destroying an
 * unfinished job requests cancellation and never joins a worker. Completed
 * streams can retain consumer-owned staging after the job is destroyed.
 * Detached source mutation
 * happens on workers; clients may only adopt a completed source on the scene
 * owner thread.
 */

#ifndef BOBOL_BSOURCEREALIZATION_H
#define BOBOL_BSOURCEREALIZATION_H

#include "BObol/BDefines.h"

#include <Inventor/SbBasic.h>

#include <stddef.h>
#include <stdint.h>
#include <memory>
#include <vector>

class SoBRLDatabaseSource;
class BObolCompactOccurrenceStream;
struct db_i;

typedef int (*BObolSourceWarmManifestProbe)(
	SoBRLDatabaseSource *source,
	struct db_i *database,
	int drawMode,
	uint32_t sourceRevision,
	BObolCompactOccurrenceStream *stream,
	void *userData);

typedef int (*BObolSourceManifestStore)(
	struct db_i *database,
	const SoBRLDatabaseSource *source,
	BObolCompactOccurrenceStream *stream,
	void *userData);

/**
 * One move-by-contract request.  submit() consumes source and
 * snapshotSourceDatabase on success and sets both caller fields to NULL.
 * source must have one caller-owned Coin reference.  The database must be an
 * independently owned handle suitable for worker reads.  callbackContext is
 * retained through the last worker callback. The coordinator releases its
 * reference before publishing the item's terminal state, including after the
 * client drops or cancels its job handle; callbacks must not borrow shorter-lived
 * client storage through another pointer. Long operations must observe cancellation
 * on the supplied stream so coordinator shutdown can join their workers.
 */
struct BOBOL_EXPORT BObolSourceRealizationRequest {
    BObolSourceRealizationRequest(void);

    SoBRLDatabaseSource *source;
    struct db_i *snapshotSourceDatabase;
    std::shared_ptr<BObolCompactOccurrenceStream> stream;
    uint64_t clientToken;
    /* Outer source-directory/traversal reservation.  Mesh-local PoP work has
     * its own finer-grained governor.  Zero asks the coordinator to derive a
     * bounded reservation from the immutable source path and snapshot
     * directory; unknown paths retain its conservative fallback. */
    size_t estimatedWorkingSetBytes;
    int drawMode;
    SbBool allowWireFallback;
    BObolSourceWarmManifestProbe probeWarmManifest;
    BObolSourceManifestStore storeManifest;
    std::shared_ptr<void> callbackContext;
};

enum BObolSourceRealizationState {
    BOBOL_SOURCE_REALIZATION_PENDING = 0,
    BOBOL_SOURCE_REALIZATION_RUNNING = 1,
    BOBOL_SOURCE_REALIZATION_COMPLETE = 2,
    BOBOL_SOURCE_REALIZATION_FAILED = 3,
    BOBOL_SOURCE_REALIZATION_CANCELLED = 4,
    /* Admission refused before a worker opens/imports the source.  Existing
     * structural coverage remains valid and a later capacity epoch may retry
     * realization. */
    BOBOL_SOURCE_REALIZATION_CONSTRAINED = 5
};

/** Item status and borrowed completed result. source is non-NULL only for
 * COMPLETE and remains valid while the job handle is retained, including
 * after cancellation. Other terminal outcomes release their detached source
 * and database resources before publication; pending/running sources remain
 * worker-exclusive. stream and scalar status remain available in every state. */
struct BOBOL_EXPORT BObolSourceRealizationItemResult {
    BObolSourceRealizationItemResult(void);

    SoBRLDatabaseSource *source;
    std::shared_ptr<BObolCompactOccurrenceStream> stream;
    uint64_t clientToken;
    int state;
    SbBool warmManifest;
    SbBool manifestStored;
};

struct BObolSourceRealizationJobPrivate;

class BOBOL_EXPORT BObolSourceRealizationJob {
public:
    ~BObolSourceRealizationJob(void);

    /** Cancel streams and retire unstarted items without waiting for worker
     * capacity. Active callbacks retain their resources until they return. */
    void cancel(void);
    /** Individual stream cancellation does not cancel siblings. A batch can
     * complete with both CANCELLED and COMPLETE item results. */
    int state(void) const;
    SbBool isTerminal(void) const;
    size_t itemCount(void) const;
    SbBool itemResult(size_t index,
	BObolSourceRealizationItemResult &result) const;

private:
    explicit BObolSourceRealizationJob(
	const std::shared_ptr<BObolSourceRealizationJobPrivate> &state);
    BObolSourceRealizationJob(const BObolSourceRealizationJob &) = delete;
    BObolSourceRealizationJob &operator=(
	const BObolSourceRealizationJob &) = delete;

    std::shared_ptr<BObolSourceRealizationJobPrivate> p;
    friend class BObolSourceRealizationCoordinator;
};

struct BObolSourceRealizationCoordinatorPrivate;

/**
 * Process-wide source realization service.  Tasks from independent views share
 * a bounded worker pool; cancellation and completion callbacks never execute
 * while the queue lock is held.
 */
class BOBOL_EXPORT BObolSourceRealizationCoordinator {
public:
    /** Initialize the shared pool. Construction errors propagate after every
     * partially started worker is joined; a later call may retry. */
    static BObolSourceRealizationCoordinator &global(void);

    /** Submit a complete batch atomically. Invalid input, a closed service,
     * or allocation failure returns an empty handle and preserves every
     * request resource and stream. Resource admission denial instead returns
     * a terminal CONSTRAINED job which owns the consumed requests. */
    std::shared_ptr<BObolSourceRealizationJob> submit(
	std::vector<BObolSourceRealizationRequest> &requests);

    /** Stop accepting work, cancel queued and active jobs, and wait for all
     * worker and caller-owned queue cleanup to finish.  Shutdown is
     * idempotent; a stopped process-wide coordinator cannot be restarted. */
    void shutdown(void);

    size_t workerCountForDiagnostics(void) const;
    size_t queuedItemCountForDiagnostics(void) const;
    size_t activeItemCountForDiagnostics(void) const;
    size_t activeWorkingSetBytesForDiagnostics(void) const;
    size_t workingSetLimitBytesForDiagnostics(void) const;

private:
    BObolSourceRealizationCoordinator(void);
    ~BObolSourceRealizationCoordinator(void);
    BObolSourceRealizationCoordinator(
	const BObolSourceRealizationCoordinator &) = delete;
    BObolSourceRealizationCoordinator &operator=(
	const BObolSourceRealizationCoordinator &) = delete;

    std::shared_ptr<BObolSourceRealizationCoordinatorPrivate> p;
};

#endif
