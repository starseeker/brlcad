/*          T R A N S A C T I O N _ F A U L T _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_TRANSACTION_FAULT_PRIVATE_H
#define LIBBOBOL_TRANSACTION_FAULT_PRIVATE_H

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
    DURABLE_CACHE_COMMIT
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
    }
    return "unknown";
}

inline bool
bobol_transaction_fault_requested(BObolTransactionFaultPoint point)
{
    const char *configured = std::getenv(BObolTransactionFaultEnvironment);
    return configured && configured[0] &&
	std::strcmp(configured, bobol_transaction_fault_name(point)) == 0;
}

#endif /* LIBBOBOL_TRANSACTION_FAULT_PRIVATE_H */
