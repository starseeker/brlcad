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
 * request.  A client-owned job is only an interest handle: destroying it
 * requests cancellation and never joins a worker.  Detached source mutation
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
 * retained through the last worker callback, including after the client drops
 * or cancels its job handle; callbacks must not borrow shorter-lived client
 * storage through another pointer.
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

/** Borrowed result view.  It remains valid while the job handle is retained. */
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

    void cancel(void);
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
    static BObolSourceRealizationCoordinator &global(void);

    std::shared_ptr<BObolSourceRealizationJob> submit(
	std::vector<BObolSourceRealizationRequest> &requests);

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

    BObolSourceRealizationCoordinatorPrivate *p;
};

#endif
