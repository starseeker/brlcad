/*              R E T A I N E D _ A L L O C A T I O N _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_RETAINED_ALLOCATION_PRIVATE_H
#define LIBBOBOL_RETAINED_ALLOCATION_PRIVATE_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

#include "lod_revision_private.h"

class BObolRetainedAllocationTransaction;
class BObolViewLodState;
class SoBRLDatabaseSource;

/* The service-wide resident epoch is allocation evidence only while this view
 * carries a memory-denial witness.  Ordinary cache growth and compaction in a
 * shared service must not invalidate an unrelated drawable scene plan. */
uint64_t bobol_retained_allocation_resident_admission_revision(
    const BObolViewLodState *viewState, uint64_t serviceRevision);

/* Canonical semantic identity of one retained-allocation request.  Fields
 * which are inactive under current policy are normalized before they enter
 * this value.  Transaction matching and completed-certificate reuse must use
 * this key rather than independently comparing raw controller inputs. */
struct BObolRetainedAllocationInputKey {
    size_t externalPresentationCost = 0;
    size_t sceneBudget = 0;
    size_t maximumMarginalBudget = 0;
    size_t maximumProtectedBudget = 0;
    BObolLodAdmissionRevisionStamp revisionStamp;
    /* Provider-memory denials are valid only for the admission epoch which
     * observed them.  Allocation identity must include that epoch so a
     * reclaimed-capacity edge can reopen richer resident candidates. */
    uint64_t residentAdmissionRevision = 0;
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
    BObolLodAdmissionRevisionStamp revisionStamp;
    uint64_t residentAdmissionRevision = 0;
    float pointProxyPixelThreshold = 0.0f;

    uint64_t viewRevision(void) const
    {
	return this->revisionStamp.view.value();
    }

    uint64_t policyRevision(void) const
    {
	return this->revisionStamp.policy.value();
    }

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

struct BObolRetainedProjectionRefreshPlan {
    SoBRLDatabaseSource *source = NULL;
    std::vector<uint32_t> compactEntryIndices;
    /* A payload without a current compact entry identity cannot be refreshed
     * selectively.  Its source must perform the ordinary dense pass. */
    bool denseRefreshRequired = false;
};

struct BObolRetainedAllocationResult {
    double maximumNormalizedError = std::numeric_limits<double>::infinity();
    double maximumProjectedErrorPixels =
	std::numeric_limits<double>::infinity();
    double visualImportanceDebt = 0.0;
    size_t protectedFloorBudget = 0;
    uint64_t protectedFloorSignature = 0;
    /* Ephemeral identity of the complete occurrence representation selected
     * by this transaction.  It is valid only inside the semantic revision
     * tuple below and lets capacity search recognize different numeric
     * budgets which map to the same discrete PoP population. */
    uint64_t selectedPopulationSignature = 0;
    /* The selected PoP population is discrete even though capacity search is
     * expressed in numeric render-cost units.  When known, this is the least
     * budget which can select a strictly richer population.  A zero value with
     * nextDistinctPresentationBudgetKnown set means the current population is
     * already the complete resident pixel-demand endpoint. */
    size_t nextDistinctPresentationBudget = 0;
    bool nextDistinctPresentationBudgetKnown = false;
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
    BObolLodAdmissionRevisionStamp revisionStamp;
    uint64_t residentAdmissionRevision = 0;
    float pointProxyPixelThreshold = 0.0f;
    size_t requestedSceneBudget = 0;
    size_t externalPresentationCost = 0;
    size_t fixedCadPresentationCost = 0;
    size_t maximumMarginalBudget = 0;
    size_t maximumProtectedBudget = 0;
    /* Candidates eligible at the allocation's current point threshold. */
    size_t pointProxyCandidateCount = 0;
    /* Candidates which could become eligible within the supported point
     * threshold range.  This is the structural proof that threshold
     * calibration has useful work to do; it must not be inferred from the
     * current-threshold population. */
    size_t reachablePointProxyCandidateCount = 0;
    /* Retained view-dependent payloads whose projection evidence does not
     * belong to this allocation epoch.  This makes the result a typed refresh
     * request rather than an allocation certificate: allocationPlanSerial and
     * all selected/certified population fields remain zero.  The source pass
     * must refresh these entries before one successor allocation may begin. */
    size_t unresolvedViewDependentPayloadCount = 0;
    std::vector<BObolRetainedProjectionRefreshPlan> projectionRefreshPlans;
    size_t selectedPointProxyCount = 0;
    /* The transaction changed at least one assembly's occurrence-level
     * point/mesh classification policy.  Unlike a PoP cut mutation this does
     * not stage CadPayload metadata, so it carries its own presentation
     * obligation. */
    bool pointProxyProtectionChanged = false;
    size_t prominentCandidateCount = 0;
    size_t prominentQualityFloorViolationCount = 0;
    bool allowProtectedFloor = false;

    uint64_t viewRevision(void) const
    {
	return this->revisionStamp.view.value();
    }

    uint64_t policyRevision(void) const
    {
	return this->revisionStamp.policy.value();
    }

    /** Return the plan only while it certifies the current control problem.
     * The renderer may safely keep displaying an older committed plan while
     * its replacement is being calculated; that fallback is not a current
     * planning certificate and must not be reported as one. */
    uint64_t currentPlanSerial(
	const BObolLodAdmissionRevisionStamp &currentRevision,
	uint64_t activePlanSerial) const
    {
	return this->allocationPlanSerial != 0 &&
	    this->allocationPlanSerial == activePlanSerial &&
	    this->revisionStamp.same(currentRevision) ?
	    this->allocationPlanSerial : 0;
    }

    BObolRetainedAllocationInputKey inputKey(void) const;

    /** Whether @p inputs can reuse this plan after it selected the complete
     * resident pixel-demand endpoint.  Enabling or disabling a protected-floor
     * trial cannot change a population which is already maximally rich; every
     * other allocation input remains semantic and must match exactly. */
    bool pixelDemandInputEquivalent(
	const BObolRetainedAllocationInputs &inputs) const;

    /* Every positive scene budget is a valid input to the deterministic
     * marginal allocator.  maximumMarginalBudget and maximumProtectedBudget
     * describe extra allowance enabled for this allocation pass; they are not
     * limits on what a later capacity candidate can request.  The richest
     * currently resident pixel-demand population is therefore the bounded
     * search endpoint.  Conflating the current pass allowance with that
     * endpoint made a conservative throughput estimate self-certifying and
     * prevented static views from ever testing visibly richer populations. */
    size_t maximumCapacitySearchBudget(void) const
    {
	return this->pixelDemandPresentationCost;
    }

    /** Whether this allocation selected every view-dependent representation
     * required by the current pixel-error demand.  Presentation and revision
     * validity are separate proofs owned by the controller. */
    bool selectsPixelDemand(void) const
    {
	return this->pixelDemandPresentationCost > 0 &&
	    this->selectedPresentationCost ==
		this->pixelDemandPresentationCost;
    }

    /** Whether the selected population is exactly the protected floor.
     * Marginal upgrades always have positive cost, so a selected cost above
     * protectedFloorBudget names a strictly richer population.  A deadline
     * miss from that richer population is not evidence against the floor. */
    bool selectsProtectedFloor(void) const
    {
	return this->protectedFloorBudget > 0 &&
	    this->selectedPresentationCost == this->protectedFloorBudget;
    }
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
    unsigned int visualEmphasis = 0;
    bool qualityFloorViolation = false;
    double normalizedError = 0.0;
    double visualBenefitPerCost = 0.0;
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

/** True only while the retained transaction has executable planning phases. */
bool bobol_retained_allocation_pending(
    const std::shared_ptr<BObolRetainedAllocationTransaction> &transaction);

/** Emit opt-in per-phase timing for a completed transaction. */
void bobol_retained_allocation_trace(
    const std::shared_ptr<BObolRetainedAllocationTransaction> &transaction);

#endif /* LIBBOBOL_RETAINED_ALLOCATION_PRIVATE_H */
