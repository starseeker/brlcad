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
    /* Render cost which belongs to the same scene frame but is not owned by
     * the CAD payloads enumerated by sources.  Scene budgets and deadline
     * samples cover the whole frame; omitting this fixed portion makes every
     * otherwise valid CAD allocation exceed its certificate by a constant
     * amount and produces repeated allocate/present/recover cycles. */
    size_t externalPresentationCost = 0;
    size_t sceneBudget = 0;
    /* A static view may spend a separately measured, hard-deadline-limited
     * remainder through occurrence-local marginal upgrades after an atomic
     * protected floor has been rejected.  This is deliberately distinct from
     * maximumProtectedBudget: it never authorizes re-testing the rejected
     * all-or-nothing floor. */
    size_t maximumMarginalBudget = 0;
    bool allowProtectedFloor = false;
    size_t maximumProtectedBudget = 0;
    uint64_t viewRevision = 0;
    uint64_t policyRevision = 0;
    float pointProxyPixelThreshold = 0.0f;
};

struct BObolRetainedAllocationResult {
    double normalizedError = std::numeric_limits<double>::infinity();
    double maximumProjectedErrorPixels =
	std::numeric_limits<double>::infinity();
    size_t protectedFloorBudget = 0;
    uint64_t protectedFloorSignature = 0;
    /*
     * Certificate for the complete occurrence-local presentation selected by
     * this transaction.  The temporary renderer-wide progressive ceiling is
     * deliberately absent from the allocator input: it protects the last
     * completed framebuffer while this plan is built, but is not part of the
     * population the plan must prove affordable.  A host may retire that
     * ceiling only while allocationPlanSerial remains the view state's active
     * plan and selectedPresentationCost fits certifiedPresentationBudget.
     */
    size_t selectedPresentationCost = 0;
    size_t certifiedPresentationBudget = 0;
    uint64_t allocationPlanSerial = 0;
    uint64_t cadRevision = 0;
    uint64_t residentDemandRevision = 0;
    uint64_t viewRevision = 0;
    uint64_t policyRevision = 0;
    float pointProxyPixelThreshold = 0.0f;
    size_t requestedSceneBudget = 0;
    size_t externalPresentationCost = 0;
    size_t fixedCadPresentationCost = 0;
    size_t maximumMarginalBudget = 0;
    size_t maximumProtectedBudget = 0;
    bool allowProtectedFloor = false;
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
 * PENDING at a bounded safe point.  STALE discards the invalid transaction so
 * the next call starts from current inputs.  Only COMPLETE has atomically
 * published allocation metadata.
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
