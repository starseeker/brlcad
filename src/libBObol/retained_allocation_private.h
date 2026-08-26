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

/* Canonical semantic identity of one retained-allocation request.  Fields
 * which are inactive under current policy are normalized before they enter
 * this value.  Transaction matching and completed-certificate reuse must use
 * this key rather than independently comparing raw controller inputs. */
struct BObolRetainedAllocationInputKey {
    size_t externalPresentationCost = 0;
    size_t sceneBudget = 0;
    size_t maximumMarginalBudget = 0;
    size_t maximumProtectedBudget = 0;
    uint64_t viewRevision = 0;
    uint64_t policyRevision = 0;
    float pointProxyPixelThreshold = 0.0f;
    bool allowProtectedFloor = false;

    bool operator==(const BObolRetainedAllocationInputKey &other) const;
    bool operator!=(const BObolRetainedAllocationInputKey &other) const
    {
	return !(*this == other);
    }
};

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

    /* Inactive policy inputs are not evidence.  Timing calibration may keep
     * updating the diagnostic protected-floor allowance after a deadline
     * handoff has disabled that trial.  Canonicalizing it prevents those
     * irrelevant samples from invalidating an otherwise identical retained
     * allocation transaction. */
    size_t effectiveMaximumProtectedBudget(void) const
    {
	return this->allowProtectedFloor ? this->maximumProtectedBudget : 0;
    }

    /* Reconciliation translates one completed renderer-wide safety cut into
     * occurrence-local cuts under that frame's exact certified allowance.
     * Throughput EMA updates and protected-floor estimates are not inputs to
     * this transaction.  Canonicalize them here so completed-plan reuse and
     * transaction matching cannot be invalidated by irrelevant frame jitter. */
    void setPresentationReconciliationBudget(size_t budget)
    {
	if (!budget || budget == SIZE_MAX)
	    return;
	this->sceneBudget = budget;
	this->maximumMarginalBudget = budget;
	this->allowProtectedFloor = false;
	this->maximumProtectedBudget = 0;
    }

    BObolRetainedAllocationInputKey inputKey(void) const;
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
    /* Complete view-dependent demand for the same candidate population.
     * This is the upper endpoint for renderer-capacity search; it is not the
     * mesh's full topology unless the current projection requires it. */
    size_t pixelDemandPresentationCost = 0;
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
    size_t pointProxyCandidateCount = 0;
    bool allowProtectedFloor = false;

    BObolRetainedAllocationInputKey inputKey(void) const;
};

enum BObolRetainedAllocationStatus {
    BOBOL_RETAINED_ALLOCATION_PENDING = 0,
    BOBOL_RETAINED_ALLOCATION_COMPLETE,
    BOBOL_RETAINED_ALLOCATION_STALE,
    BOBOL_RETAINED_ALLOCATION_FAILED
};

/* One coherent occurrence-local quality transition.  This lives in the
 * private interface so the exhaustive policy oracle can test the production
 * ordering directly rather than reproduce a look-alike comparator. */
struct BObolRetainedMarginalUpgrade {
    size_t candidateIndex = 0;
    int nextCut = -1;
    size_t nextCost = 0;
    size_t addedCost = 0;
    bool qualityFloorViolation = false;
    double weightedError = 0.0;
    double valuePerCost = 0.0;
};

/** True when a ranks below b in the retained allocator's marginal queue. */
bool bobol_retained_marginal_lower_priority(
    const BObolRetainedMarginalUpgrade &a,
    const BObolRetainedMarginalUpgrade &b);

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
