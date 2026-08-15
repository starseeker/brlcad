/*              R E T A I N E D _ A L L O C A T I O N _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_RETAINED_ALLOCATION_PRIVATE_H
#define LIBBOBOL_RETAINED_ALLOCATION_PRIVATE_H

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

class BObolRetainedAllocationTransaction;
class BObolViewLodState;
class SoBRLDatabaseSource;

struct BObolRetainedAllocationInputs {
    const std::vector<SoBRLDatabaseSource *> *sources = NULL;
    BObolViewLodState *viewState = NULL;
    size_t sceneBudget = 0;
    bool allowProtectedFloor = false;
    size_t maximumProtectedBudget = 0;
    uint64_t viewRevision = 0;
    uint64_t policyRevision = 0;
    int progressiveCutCeiling = -1;
    float pointProxyPixelThreshold = 0.0f;
};

struct BObolRetainedAllocationResult {
    double normalizedError = std::numeric_limits<double>::infinity();
    double maximumProjectedErrorPixels =
	std::numeric_limits<double>::infinity();
    size_t protectedFloorBudget = 0;
    uint64_t protectedFloorSignature = 0;
};

enum BObolRetainedAllocationStatus {
    BOBOL_RETAINED_ALLOCATION_PENDING = 0,
    BOBOL_RETAINED_ALLOCATION_COMPLETE,
    BOBOL_RETAINED_ALLOCATION_STALE,
    BOBOL_RETAINED_ALLOCATION_FAILED
};

/**
 * Advance one owner-thread retained-allocation transaction.  A zero slice
 * executes synchronously; a nonzero slice preserves its cursor and returns
 * PENDING at a bounded safe point.  Only COMPLETE has atomically published
 * allocation metadata.
 */
BObolRetainedAllocationStatus bobol_retained_allocation_advance(
    std::shared_ptr<BObolRetainedAllocationTransaction> &transaction,
    const BObolRetainedAllocationInputs &inputs,
    uint64_t sliceMicroseconds,
    BObolRetainedAllocationResult &result);

/** Emit opt-in per-phase timing for a completed transaction. */
void bobol_retained_allocation_trace(
    const std::shared_ptr<BObolRetainedAllocationTransaction> &transaction);

#endif /* LIBBOBOL_RETAINED_ALLOCATION_PRIVATE_H */
