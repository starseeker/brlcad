/*          T R A N S A C T I O N _ F A U L T _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_TRANSACTION_FAULT_PRIVATE_H
#define LIBBOBOL_TRANSACTION_FAULT_PRIVATE_H


#include "bu/str.h"
#include <cstdlib>
#include <cstring>

/* Deterministic failure points for executable transaction contracts.  These
 * are private test instrumentation: normal execution performs one unset
 * environment lookup and has no injected state or allocation. */
enum class BObolTransactionFaultPoint {
    RETAINED_SCENE_COMMIT = 0,
    PRESENTATION_COMMIT,
    DURABLE_CACHE_OPEN,
    DURABLE_CACHE_WRITE,
    DURABLE_CACHE_COMMIT,
    SOURCE_REALIZATION_QUEUE_APPEND,
    SOURCE_REALIZATION_JOB_HANDLE,
    SOURCE_REALIZATION_WORKER_START,
    SOURCE_STREAM_MERGE_AFTER_COMPLETION,
    SOURCE_STREAM_PARTIAL_MERGE_AFTER_COMPLETION,
    SOURCE_STREAM_RESERVE_AFTER_COMPLETION,
    SOURCE_TERMINAL_PREPARATION
};

static constexpr const char *BObolTransactionFaultEnvironment =
    "BOBOL_TEST_TRANSACTION_FAULT";

inline const char *
bobol_transaction_fault_name(BObolTransactionFaultPoint point)
{
    switch (point) {
	case BObolTransactionFaultPoint::RETAINED_SCENE_COMMIT:
	    return "retained-scene-commit";
	case BObolTransactionFaultPoint::PRESENTATION_COMMIT:
	    return "presentation-commit";
	case BObolTransactionFaultPoint::DURABLE_CACHE_OPEN:
	    return "durable-cache-open";
	case BObolTransactionFaultPoint::DURABLE_CACHE_WRITE:
	    return "durable-cache-write";
	case BObolTransactionFaultPoint::DURABLE_CACHE_COMMIT:
	    return "durable-cache-commit";
	case BObolTransactionFaultPoint::SOURCE_REALIZATION_QUEUE_APPEND:
	    return "source-realization-queue-append";
	case BObolTransactionFaultPoint::SOURCE_REALIZATION_JOB_HANDLE:
	    return "source-realization-job-handle";
	case BObolTransactionFaultPoint::SOURCE_REALIZATION_WORKER_START:
	    return "source-realization-worker-start";
	case BObolTransactionFaultPoint::SOURCE_STREAM_MERGE_AFTER_COMPLETION:
	    return "source-stream-merge-after-completion";
	case BObolTransactionFaultPoint::SOURCE_STREAM_PARTIAL_MERGE_AFTER_COMPLETION:
	    return "source-stream-partial-merge-after-completion";
	case BObolTransactionFaultPoint::SOURCE_STREAM_RESERVE_AFTER_COMPLETION:
	    return "source-stream-reserve-after-completion";
	case BObolTransactionFaultPoint::SOURCE_TERMINAL_PREPARATION:
	    return "source-terminal-preparation";
    }
    return "unknown";
}

inline bool
bobol_transaction_fault_requested(BObolTransactionFaultPoint point)
{
    const char *configured = std::getenv(BObolTransactionFaultEnvironment);
    return configured && configured[0] &&
	bu_strcmp(configured, bobol_transaction_fault_name(point)) == 0;
}

#endif /* LIBBOBOL_TRANSACTION_FAULT_PRIVATE_H */
