/*      L O D _ C A P A C I T Y _ P O L I C Y _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_CAPACITY_POLICY_PRIVATE_H
#define LIBBOBOL_LOD_CAPACITY_POLICY_PRIVATE_H

#include "common.h"

#include "lod_capacity_search_private.h"
#include "lod_revision_private.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>

class BObolLodAdmissionPlanner;

class BObolLodCapacityEvidence;

/*
 * A coverage or publication pass may need exactly one coherent framebuffer
 * before planning resumes.  That ordering edge is not a renderer-capacity
 * sample and must not share the bounded search's candidate counters.  A
 * quality-blocked ordinary pass additionally transfers its successor to the
 * complete importance allocator, which is the only path allowed to start a
 * numeric capacity search.
 */
class BObolLodCapacityFrameBarrier {
public:
    enum class Successor : uint8_t {
	NONE = 0,
	REPLAN,
	REALLOCATE
    };

    void requestReplan(void)
    {
	if (this->successorValue == Successor::NONE)
	    this->successorValue = Successor::REPLAN;
    }

    void requestReallocation(void)
    {
	this->successorValue = Successor::REALLOCATE;
    }

    Successor consume(void)
    {
	const Successor successor = this->successorValue;
	this->successorValue = Successor::NONE;
	return successor;
    }

    void reset(void)
    {
	this->successorValue = Successor::NONE;
    }

    bool pending(void) const
    {
	return this->successorValue != Successor::NONE;
    }

    Successor successor(void) const
    {
	return this->successorValue;
    }

private:
    Successor successorValue = Successor::NONE;
};

static_assert(std::is_trivially_copyable<
	BObolLodCapacityFrameBarrier>::value,
    "capacity frame barriers must remain allocation-free values");

/*
 * One retained occurrence-allocation request may be pending at a time.  The
 * former pending/preserve/reconciliation-budget field trio admitted invalid
 * combinations (for example, a reconciliation budget with no request) and
 * required callers to reset the companions together.  Keep that control edge
 * as a finite value; the numeric budget remains meaningful only for the
 * presentation-reconciliation alternative.
 */
class BObolLodRetainedAllocationRequest {
public:
    enum class Kind : uint8_t {
	NONE = 0,
	PRESERVE_BUDGET,
	RECOMPUTE_BUDGET,
	PRESENTATION_RECONCILIATION
    };

    void reset(void)
    {
	this->kindValue = Kind::NONE;
	this->reconciliationBudgetValue = SIZE_MAX;
    }

    void requestReallocation(bool preserveBudget)
    {
	/* A reconciliation request carries a completed-frame capacity proof.
	 * A weaker reallocation request cannot discard or reinterpret it. */
	if (this->kindValue == Kind::PRESENTATION_RECONCILIATION)
	    return;
	this->kindValue = preserveBudget ? Kind::PRESERVE_BUDGET :
	    Kind::RECOMPUTE_BUDGET;
	this->reconciliationBudgetValue = SIZE_MAX;
    }

    void requestPresentationReconciliation(size_t budget)
    {
	if (!budget || budget == SIZE_MAX)
	    return;
	if (this->kindValue != Kind::PRESENTATION_RECONCILIATION) {
	    this->kindValue = Kind::PRESENTATION_RECONCILIATION;
	    this->reconciliationBudgetValue = budget;
	    return;
	}
	this->reconciliationBudgetValue = std::min(
	    this->reconciliationBudgetValue, budget);
    }

    bool pending(void) const { return this->kindValue != Kind::NONE; }
    bool preservesBudget(void) const
    {
	return this->kindValue == Kind::PRESERVE_BUDGET;
    }
    bool reconcilesPresentation(void) const
    {
	return this->kindValue == Kind::PRESENTATION_RECONCILIATION;
    }
    Kind kind(void) const { return this->kindValue; }
    size_t reconciliationBudget(void) const
    {
	return this->reconcilesPresentation() ?
	    this->reconciliationBudgetValue : SIZE_MAX;
    }

private:
    Kind kindValue = Kind::NONE;
    size_t reconciliationBudgetValue = SIZE_MAX;
};

static_assert(std::is_trivially_copyable<
	BObolLodRetainedAllocationRequest>::value,
    "retained allocation requests must remain allocation-free values");

/*
 * Mechanical progress through one bounded admission plan.  This value owns no
 * renderer-capacity conclusion and survives no explicit reset.  Separating it
 * from completed-frame evidence lets a large scene resume across owner-thread
 * windows without making cursor movement look like new planning evidence.
 */
class BObolLodAdmissionCursor {
public:
    void reset(void)
    {
	this->initializedValue = false;
	this->singleOccurrenceBootstrapValue = false;
	this->activeCostValue = 0;
	this->minimumActiveCostValue = 0;
	this->refinementRemainingValue = SIZE_MAX;
	this->retainedAdmissionValue = false;
	this->retainedAdmissionRemainingValue = SIZE_MAX;
	this->revisionStampValue = BObolLodAdmissionRevisionStamp();
    }

    void consumeRefinement(size_t cost)
    {
	if (this->refinementRemainingValue == SIZE_MAX)
	    return;
	this->refinementRemainingValue =
	    cost >= this->refinementRemainingValue ?
		0 : this->refinementRemainingValue - cost;
    }

    void consumeRetainedAdmission(size_t cost)
    {
	if (this->retainedAdmissionRemainingValue == SIZE_MAX)
	    return;
	this->retainedAdmissionRemainingValue =
	    cost >= this->retainedAdmissionRemainingValue ?
		0 : this->retainedAdmissionRemainingValue - cost;
    }

    bool initialized(void) const { return this->initializedValue; }
    bool singleOccurrenceBootstrap(void) const
    {
	return this->singleOccurrenceBootstrapValue;
    }
    size_t activeCost(void) const { return this->activeCostValue; }
    size_t minimumActiveCost(void) const
    {
	return this->minimumActiveCostValue;
    }
    size_t refinementRemaining(void) const
    {
	return this->refinementRemainingValue;
    }
    bool retainedAdmission(void) const
    {
	return this->retainedAdmissionValue;
    }
    size_t retainedAdmissionRemaining(void) const
    {
	return this->retainedAdmissionRemainingValue;
    }
    bool matches(const BObolLodAdmissionRevisionStamp &stamp) const
    {
	return this->revisionStampValue.same(stamp);
    }
    const BObolLodAdmissionRevisionStamp &revisionStamp(void) const
    {
	return this->revisionStampValue;
    }

private:
    friend class BObolLodCapacityEvidence;
    friend class BObolLodAdmissionPlanner;

    void resetForRevision(const BObolLodAdmissionRevisionStamp &stamp)
    {
	this->reset();
	this->revisionStampValue = stamp;
    }

    bool initializedValue = false;
    bool singleOccurrenceBootstrapValue = false;
    size_t activeCostValue = 0;
    size_t minimumActiveCostValue = 0;
    size_t refinementRemainingValue = SIZE_MAX;
    bool retainedAdmissionValue = false;
    size_t retainedAdmissionRemainingValue = SIZE_MAX;
    BObolLodAdmissionRevisionStamp revisionStampValue;
};

static_assert(std::is_trivially_copyable<BObolLodAdmissionCursor>::value,
    "admission cursors must remain allocation-free values");

class BObolLodCapacityEvidence {
public:
    struct Inputs {
	size_t activeCost = 0;
	size_t minimumActiveCost = 0;
	float targetFps = 0.0f;
	long double calibratedCostPerSecond = 0.0L;
	uint64_t observedStableNanoseconds = 0;
	/* The preferred quiet cadence controls refinement pacing.  A completed
	 * frame below this separate hard limit is still a valid static
	 * presentation and must not be destructively coarsened merely for missing
	 * the preferred cadence. */
	uint64_t hardPresentationDeadlineNanoseconds = 0;
	uint64_t lastRenderNanoseconds = 0;
	uint64_t smoothedRenderNanoseconds = 0;
	bool interactive = false;
	bool scaleQualityProbe = false;
	/* A bounded static presentation may use the complete hard-frame allowance
	 * immediately.  This covers both a single-occurrence quality step and the
	 * minimum-mesh repair which replaces the last structural boxes in a
	 * many-part view.  The endpoint deadline can abort an optimistic trial, so
	 * imposing the ordinary 4x/1.25x refinement ladder here only creates
	 * visible intermediate stages or a no-progress repair loop. */
	bool hardDeadlinePresentation = false;
	/* Orthographic pose-only handoff starts from an already useful retained
	 * population.  Its reversible motion ceiling is not evidence that those
	 * occurrence cuts are unaffordable.  Preserve at least the currently
	 * presented population until a completed deadline miss establishes a real
	 * capacity ceiling. */
	bool preserveActivePopulation = false;
	/* Replacing a visible structural box with the minimum drawable mesh is a
	 * coverage obligation, not ordinary triangle-detail recovery.  A retained
	 * recovery ceiling may continue to constrain richer prefixes after this
	 * pass, but it must not make that minimum replacement impossible when the
	 * independently interruptible hard-frame allowance proves room for it. */
	bool structuralCoverageRepair = false;
	/* With exactly one visible occurrence and no retained presentation, a
	 * conservative many-object seed manufactures several blocky calibration
	 * frames.  Permit one larger, still deadline-protected bootstrap. */
	bool coldSingleOccurrence = false;
	bool forceTerminal = false;
	bool releaseCutFloor = false;
	bool stablePresentationHandoff = false;
	size_t stablePresentationCostFloor = 0;
    };

    struct Decision {
	bool initialized = false;
	bool overloadRecovery = false;
	size_t totalBudget = 0;
	size_t refinementBudget = 0;
	bool retainedAdmission = false;
	size_t retainedAdmissionBudget = SIZE_MAX;
    };

    struct CalibrationInputs {
	size_t activeCost = 0;
	uint64_t observedNanoseconds = 0;
	bool passAdmittedWork = false;
	/* The blocked cut came from a complete scene-wide minimax allocation.
	 * If calibration demonstrates a different allowance, its successor must
	 * recompute that allocation rather than letting ordinary per-occurrence
	 * refinement consume stale allocation stamps. */
	bool allocationCutsApplied = false;
	/* A complete retained allocation supplies the ordered candidate set and
	 * its pixel-demand endpoint.  Only this path may use the bounded search;
	 * structural/provider barriers retain their separate finite contracts. */
	bool boundedSearch = false;
	/* The occurrence-local cuts named by the allocation may be active while a
	 * temporary renderer-wide ceiling still hides their exact presentation.
	 * Such a population may advance ALLOCATING to MEASURING, but it may not
	 * consume a timing sample until the ceiling-free presentation is exact. */
	bool allocationPresentationRealized = false;
	BObolLodCapacitySearchKey searchKey;
	/* The allocator may raise a requested steady allowance to an atomic
	 * prominent-feature floor.  Name the actual certified candidate budget;
	 * it is not necessarily the scalar allowance which initiated the pass. */
	size_t candidateBudget = 0;
	size_t demandedBudget = 0;
	size_t knownSafeBudget = 0;
	uint64_t populationSignature = 0;
    };

    struct CompletedFrameInputs {
	enum class CandidateState : uint8_t {
	    CURRENT = 0,
	    PRESENTATION_PENDING,
	    REALLOCATION_REQUIRED
	};

	BObolLodCapacitySearchKey searchKey;
	size_t candidateBudget = 0;
	size_t presentedCost = 0;
	uint64_t populationSignature = 0;
	size_t knownSafeBudget = 0;
	uint64_t observedNanoseconds = 0;
	bool validSample = false;
	/* A current retained allocation may still be hidden by a reversible
	 * renderer ceiling.  That state waits for presentation.  A changed
	 * occurrence population or unapplied cut requires another allocation and
	 * must never be treated as an invalid timing sample. */
	CandidateState candidateState = CandidateState::CURRENT;
	/* Resident-result growth owns the successor allocation when it is pending.
	 * Retire the obsolete search without racing that coalesced producer. */
	bool reallocationProducerPending = false;
    };

    struct CalibrationDecision {
	bool candidateReallocation = false;
	bool searchActive = false;
	bool requestFrame = false;
	bool sampleFrame = false;
	bool restartSubmission = false;
	bool searchTerminal = false;
	BObolLodCapacitySearchCertificate::Result searchResult =
	    BObolLodCapacitySearchCertificate::Result::NONE;
    };

    struct CompletedFrameDecision {
	bool requestSampleFrame = false;
	bool restartSubmission = false;
	BObolLodCapacitySearchCertificate::Result searchResult =
	    BObolLodCapacitySearchCertificate::Result::NONE;
    };

private:
    friend class BObolLodAdmissionPlanner;

    Decision planPass(const Inputs &inputs, BObolLodAdmissionCursor &cursor)
    {
	Decision decision;
	decision.totalBudget = this->currentBudgetValue;
	decision.refinementBudget = cursor.refinementRemainingValue;
	decision.retainedAdmission = cursor.retainedAdmissionValue;
	decision.retainedAdmissionBudget =
	    cursor.retainedAdmissionRemainingValue;
	/* A capacity candidate owns the exact occurrence population throughout
	 * its finite measurement series.  No throughput estimate, ordinary quiet
	 * target, or retained-recovery request may rewrite that population between
	 * samples.  A semantic revision invalidates the certificate before this
	 * point; absent such an edge, initialize a no-op admission cursor and let
	 * the completed frame supply the next sample. */
	if (this->capacitySearchValue.phase() ==
		BObolLodCapacitySearchCertificate::Phase::PRESENTING ||
	    this->capacitySearchValue.phase() ==
		BObolLodCapacitySearchCertificate::Phase::MEASURING) {
	    const size_t candidateBudget =
		this->capacitySearchValue.candidateBudget();
	    this->currentBudgetValue = candidateBudget;
	    cursor.activeCostValue = inputs.activeCost;
	    cursor.minimumActiveCostValue = inputs.minimumActiveCost;
	    cursor.refinementRemainingValue = 0;
	    cursor.retainedAdmissionValue = false;
	    cursor.retainedAdmissionRemainingValue = SIZE_MAX;
	    cursor.initializedValue = true;
	    decision.initialized = true;
	    decision.totalBudget = candidateBudget;
	    decision.refinementBudget = 0;
	    decision.retainedAdmission = false;
	    decision.retainedAdmissionBudget = SIZE_MAX;
	    return decision;
	}
	if (cursor.initializedValue)
	    return decision;
	cursor.singleOccurrenceBootstrapValue =
	    inputs.coldSingleOccurrence && inputs.activeCost == 0 &&
	    inputs.calibratedCostPerSecond <= 0.0L;
	const bool requestedRetainedRecovery =
	    this->requestedRetainedRecoveryBudgetValue != SIZE_MAX &&
	    inputs.activeCost > this->requestedRetainedRecoveryBudgetValue;
	const bool requestedRetainedReallocation =
	    this->retainedAllocationRequestValue.pending() &&
	    !inputs.interactive && !inputs.forceTerminal &&
	    inputs.activeCost > 0;
	const bool requestedPresentationReconciliation =
	    requestedRetainedReallocation &&
	    this->retainedAllocationRequestValue.reconcilesPresentation();
	const size_t presentationReconciliationBudget =
	    this->retainedAllocationRequestValue.reconciliationBudget();
	const bool requestedCoverageCompletion =
	    inputs.structuralCoverageRepair &&
	    this->requestedCoverageCompletionAdditionalCostValue != SIZE_MAX;
	const size_t coverageCompletionBudget = requestedCoverageCompletion &&
	    this->requestedCoverageCompletionAdditionalCostValue >
		SIZE_MAX - inputs.activeCost ? SIZE_MAX :
	    inputs.activeCost +
		this->requestedCoverageCompletionAdditionalCostValue;
	const size_t retainedReallocationBudget = this->currentBudgetValue;

	const long double targetNanoseconds = inputs.targetFps > 0.0f ?
	    1000000000.0L / static_cast<long double>(inputs.targetFps) : 0.0L;
	const long double observedStableNanoseconds =
	    static_cast<long double>(inputs.observedStableNanoseconds);
	const bool severeStableOverload =
	    !inputs.interactive && targetNanoseconds > 0.0L &&
	    observedStableNanoseconds > targetNanoseconds * 2.0L;
	const bool observedWithinHardPresentationDeadline =
	    !inputs.interactive && inputs.activeCost > 0 &&
	    inputs.observedStableNanoseconds > 0 &&
	    inputs.hardPresentationDeadlineNanoseconds > 0 &&
	    inputs.observedStableNanoseconds <=
		inputs.hardPresentationDeadlineNanoseconds;
	const bool deadlineRecoveryRequested =
	    this->steadyDeadlineCapacityCeilingValue != SIZE_MAX &&
	    inputs.activeCost > this->steadyDeadlineCapacityCeilingValue;
	const bool overloadRecovery =
	    !inputs.interactive && !inputs.forceTerminal &&
	    !observedWithinHardPresentationDeadline &&
	    (!this->overloadRecoveryPerformedValue ||
	     this->overloadRecoveryActiveCostValue != inputs.activeCost) &&
	    inputs.activeCost > 0 &&
	    (deadlineRecoveryRequested || severeStableOverload) &&
	    targetNanoseconds > 0.0L &&
	    observedStableNanoseconds > targetNanoseconds * 1.20L;
	/* Three unchanged, exact presentation misses are stronger evidence than
	 * the triangle-throughput estimator.  The latter cannot represent the
	 * fixed per-occurrence/classification cost of a many-part scene: Hubble's
	 * measured triangle rate predicted headroom while three successive 75-85
	 * ms frames missed the 50 ms stable target.  Requiring both witnesses made
	 * the controller shave a few cost units, repeat three probes, and take
	 * tens of seconds to settle.  The probe count already rejects transient
	 * setup noise; use its direct duration with the safety-scaled recovery
	 * calculation below. */

	size_t costBudget = inputs.coldSingleOccurrence &&
	    inputs.activeCost == 0 &&
	    inputs.calibratedCostPerSecond <= 0.0L ?
		std::max(this->seedBudgetValue,
		    this->singleOccurrenceBootstrapBudget()) :
		this->seedBudgetValue;
	if (inputs.forceTerminal) {
	    costBudget = SIZE_MAX;
	} else if (inputs.targetFps > 0.0f &&
	    inputs.calibratedCostPerSecond > 0.0L) {
	    const long double affordable =
		inputs.calibratedCostPerSecond /
		static_cast<long double>(inputs.targetFps);
	    costBudget = affordable >= static_cast<long double>(SIZE_MAX) ?
		SIZE_MAX : std::max<size_t>(
		    1, static_cast<size_t>(affordable));
	    if (inputs.activeCost > 0 && costBudget > inputs.activeCost &&
		!inputs.hardDeadlinePresentation) {
		size_t growthNumerator = 4;
		size_t growthDenominator = 1;
		if (!inputs.interactive || inputs.scaleQualityProbe) {
		    if (targetNanoseconds > 0.0L &&
			observedStableNanoseconds >=
			    targetNanoseconds * 0.50L) {
			growthNumerator = 5;
			growthDenominator = 4;
		    } else if (targetNanoseconds > 0.0L &&
			observedStableNanoseconds >=
			    targetNanoseconds * 0.25L) {
			growthNumerator = 2;
		    }
		}
		const size_t growthBase =
		    !inputs.interactive && this->currentBudgetValue != SIZE_MAX ?
			std::max(inputs.activeCost, this->currentBudgetValue) :
			inputs.activeCost;
		const size_t growthIncrement =
		    growthBase > SIZE_MAX / growthNumerator ? SIZE_MAX :
		    growthBase * growthNumerator / growthDenominator;
		const size_t growthLimit =
		    std::max(this->seedBudgetValue, growthIncrement);
		costBudget = std::min(costBudget, growthLimit);
	    }
	    if (!inputs.interactive && !overloadRecovery &&
		!inputs.stablePresentationHandoff &&
		this->currentBudgetValue != SIZE_MAX)
		costBudget = std::max(costBudget, this->currentBudgetValue);
	    if (!inputs.interactive && !inputs.stablePresentationHandoff &&
		this->currentBudgetValue != SIZE_MAX &&
		costBudget > this->currentBudgetValue &&
		targetNanoseconds > 0.0L &&
		inputs.observedStableNanoseconds > 0 &&
		observedStableNanoseconds >= targetNanoseconds * 0.90L)
		costBudget = this->currentBudgetValue;
	    if (overloadRecovery && observedStableNanoseconds > 0.0L) {
		const long double recovered =
		    static_cast<long double>(inputs.activeCost) *
		    targetNanoseconds * 0.80L / observedStableNanoseconds;
		const size_t recoveredBudget =
		    recovered >= static_cast<long double>(SIZE_MAX) ? SIZE_MAX :
		    std::max<size_t>(1, static_cast<size_t>(recovered));
		costBudget = std::min(costBudget, recoveredBudget);
	    }
	}

	if (inputs.interactive && inputs.activeCost > 0 &&
	    targetNanoseconds > 0.0L) {
	    const uint64_t observedNanoseconds = std::max(
		inputs.lastRenderNanoseconds, inputs.smoothedRenderNanoseconds);
	    if (observedNanoseconds > 0 &&
		static_cast<long double>(observedNanoseconds) <=
		    targetNanoseconds * 1.05L) {
		const long double timingHeadroom =
		    inputs.scaleQualityProbe ? 1.0L : 0.80L;
		const long double affordable =
		    static_cast<long double>(inputs.activeCost) *
		    targetNanoseconds * timingHeadroom /
		    static_cast<long double>(observedNanoseconds);
		const size_t growthLimit = inputs.activeCost > SIZE_MAX / 4 ?
		    SIZE_MAX : inputs.activeCost * 4;
		size_t responsiveBudget =
		    affordable >= static_cast<long double>(growthLimit) ?
			growthLimit : static_cast<size_t>(
			    std::max<long double>(
				static_cast<long double>(inputs.activeCost),
				affordable));
		if (inputs.scaleQualityProbe &&
		    static_cast<long double>(observedNanoseconds) <
			targetNanoseconds * 0.50L)
		    responsiveBudget = growthLimit;
		costBudget = std::max(costBudget, responsiveBudget);
	    }
	}
	if (inputs.releaseCutFloor && this->currentBudgetValue != SIZE_MAX)
	    costBudget = std::max(costBudget, this->currentBudgetValue);
	/* observedStableNanoseconds is supplied only for an exact, reusable CAD
	 * presentation.  If that completed population met this pass's deadline,
	 * it is a direct proof that the population is affordable.  A triangle-rate
	 * estimate cannot contradict that proof: its linear currency omits fixed
	 * per-occurrence and command-dispatch costs and can otherwise place a
	 * many-part scene perpetually below its already displayed population. */
	const bool measuredActivePopulationDeadlineSafe =
	    !inputs.interactive && !inputs.stablePresentationHandoff &&
	    inputs.activeCost > 0 &&
	    targetNanoseconds > 0.0L && observedStableNanoseconds > 0.0L &&
	    observedStableNanoseconds <= targetNanoseconds;
	if (measuredActivePopulationDeadlineSafe && costBudget != SIZE_MAX)
	    costBudget = std::max(costBudget, inputs.activeCost);
	/* Static presentation is allowed to use the hard deadline.  Preserve an
	 * exact population which meets it even when the preferred quiet target is
	 * stricter; otherwise a pose-only rotation can replace a proven mesh view
	 * with boxes before structural repair has a chance to run. */
	if (observedWithinHardPresentationDeadline && costBudget != SIZE_MAX)
	    costBudget = std::max(costBudget, inputs.activeCost);
	/* A reusable presentation of this exact population is also the strongest
	 * available headroom estimate for the interruptible static-quality pass.
	 * The long-horizon throughput EMA can lag far behind after point
	 * aggregation or a large structural publication wave; using only that EMA
	 * left a 27 ms OSMesa frame with no allowance to replace its last few
	 * hundred boxes.  Extrapolate conservatively to eighty percent of the hard
	 * deadline and cap one step at four times the measured population.  The
	 * endpoint deadline remains the independent correctness guard if the cost
	 * model is non-linear. */
	if (inputs.hardDeadlinePresentation &&
	    measuredActivePopulationDeadlineSafe && costBudget != SIZE_MAX) {
	    const long double directAffordable =
		static_cast<long double>(inputs.activeCost) *
		targetNanoseconds * 0.80L / observedStableNanoseconds;
	    const size_t directLimit = inputs.activeCost > SIZE_MAX / 4 ?
		SIZE_MAX : inputs.activeCost * 4;
	    const size_t directBudget =
		!std::isfinite(directAffordable) || directAffordable <= 0.0L ?
		    inputs.activeCost :
		    (directAffordable >= static_cast<long double>(directLimit) ?
			directLimit : static_cast<size_t>(directAffordable));
	    costBudget = std::max(costBudget, directBudget);
	}
	/* Preserve only work which an exact completed presentation proved could
	 * meet the stable deadline.  This is deliberately a measured render-cost
	 * floor, not the previous pass allowance: the latter may include admission
	 * which was never submitted or presented. */
	if (inputs.stablePresentationCostFloor > 0 && costBudget != SIZE_MAX)
	    costBudget = std::max(
		costBudget, inputs.stablePresentationCostFloor);
	/* A view-importance reallocation deliberately changes only which
	 * occurrences consume the established allowance.  Re-running the
	 * throughput bootstrap here would replace that allowance with the cold
	 * seed before the first candidate is considered. */
	if (requestedRetainedReallocation &&
	    this->retainedAllocationRequestValue.preservesBudget())
	    costBudget = retainedReallocationBudget;
	/* Point aggregation is a last-resort population bound, not a substitute
	 * for PoP triangle allocation.  When the presentation controller is
	 * returning a temporary multi-pixel point cut to the one-pixel contract,
	 * it may explicitly ask this pass to compact reducible retained prefixes
	 * first.  The first request forces immediate retained admission; its
	 * measured ceiling remains in force for this view/policy epoch.  Without
	 * that persistent ceiling the cheaper recovery frame recalibrates a larger
	 * budget and immediately re-admits the exact discrete cut which just missed
	 * the frame deadline. */
	if (requestedRetainedRecovery)
	    costBudget = std::min(
		costBudget, this->requestedRetainedRecoveryBudgetValue);
	if (!inputs.forceTerminal && !inputs.structuralCoverageRepair &&
	    this->retainedRecoveryCeilingValue != SIZE_MAX)
	    costBudget = std::min(
		costBudget, this->retainedRecoveryCeilingValue);
	/* Recognizable prominent geometry is a harder quiet-view contract than
	 * the preferred 20 Hz cadence, provided its measured population remains
	 * inside the separate hard quality-frame deadline.  Persist that floor
	 * across overload passes; otherwise every 60-90 ms OSMesa wire frame
	 * alternates between the 20 Hz recovery budget and the same protected
	 * allocation forever. */
	if (!inputs.interactive && !inputs.forceTerminal &&
	    this->retainedQualityFloorBudgetValue > 0)
	    costBudget = std::max(
		costBudget, this->retainedQualityFloorBudgetValue);
	/* A completed hard-deadline abort is an unsafe upper-bound witness for
	 * this view/policy capacity epoch.  The constrained handoff below owns the
	 * immediate safe allocation, but its budget is intentionally one-shot.
	 * Keep a slightly lower strict ceiling after that handoff so ordinary
	 * throughput calibration cannot propose the identical failed population
	 * again.  A later static-quality pass may still try an intermediate budget
	 * below this bound, preserving useful fidelity without reopening a cycle. */
	const bool staticCapacitySearch =
	    this->capacitySearchValue.phase() !=
		BObolLodCapacitySearchCertificate::Phase::INACTIVE &&
	    this->capacitySearchValue.goal() ==
		BObolLodCapacitySearchCertificate::Goal::STATIC;
	const bool capacitySearchInactive =
	    this->capacitySearchValue.phase() ==
		BObolLodCapacitySearchCertificate::Phase::INACTIVE;
	const bool staticQualityBudget = staticCapacitySearch ||
	    (capacitySearchInactive &&
	     this->retainedQualityFloorBudgetValue > 0);
	const size_t deadlineCapacityCeiling = staticQualityBudget ?
	    this->staticDeadlineCapacityCeilingValue :
	    this->steadyDeadlineCapacityCeilingValue;
	if (!inputs.interactive && !inputs.forceTerminal &&
	    deadlineCapacityCeiling != SIZE_MAX)
	    costBudget = std::min(
		costBudget, deadlineCapacityCeiling);
	/* A pose-only handoff may redistribute the existing scene allowance after
	 * its current-camera visibility census, but it must not infer a smaller
	 * allowance from the temporary motion ceiling.  Explicit deadline and
	 * recovery ceilings are completed-frame evidence and therefore remain
	 * stronger than this continuity witness. */
	const bool posePopulationHasNoUnsafeWitness =
	    this->steadyDeadlineCapacityCeilingValue == SIZE_MAX &&
	    this->staticDeadlineCapacityCeilingValue == SIZE_MAX &&
	    this->retainedRecoveryCeilingValue == SIZE_MAX &&
	    this->requestedRetainedRecoveryBudgetValue == SIZE_MAX;
	if (inputs.preserveActivePopulation &&
	    posePopulationHasNoUnsafeWitness && costBudget != SIZE_MAX)
	    costBudget = std::max(costBudget, inputs.activeCost);
	/* A completed constrained framebuffer is a hard presentation-capacity
	 * witness, not another soft quality preference.  Reconcile the hidden
	 * retained prefixes at exactly that scene-wide allowance: a stale
	 * throughput seed may not undercut proven work, and a protected-quality
	 * floor may not exceed the frame deadline which created the handoff. */
	if (requestedPresentationReconciliation)
	    costBudget = presentationReconciliationBudget;
	/* Structural fallback above the maximum point-proxy threshold is an
	 * irreducible coverage population, not ordinary quality refinement.  An
	 * exact completed frame may transfer its separately bounded static-frame
	 * capacity into one coverage transaction.  Consume that certified budget
	 * verbatim after the heuristic/EMA clamps: allowing the preferred-cadence
	 * estimator to undercut it re-arms the same box frontier forever. */
	if (requestedCoverageCompletion)
	    costBudget = coverageCompletionBudget;

	this->currentBudgetValue = costBudget;
	/* Freeze the exact cost currencies which initialized this bounded pass.
	 * A 150k-occurrence pass may require hundreds of owner-thread windows;
	 * rediscovering the complete population before every window turns an O(N)
	 * allocation into O(N * windows) work and monopolizes the GUI thread.
	 * Occurrence changes are charged by the carried refinement remainders; a
	 * completed pass resets this snapshot before another census. */
	cursor.activeCostValue = inputs.activeCost;
	cursor.minimumActiveCostValue = inputs.minimumActiveCost;
	cursor.refinementRemainingValue = costBudget == SIZE_MAX ? SIZE_MAX :
	    (costBudget > inputs.activeCost ?
		costBudget - inputs.activeCost : 0);
	/* A presentation handoff removes a reversible renderer ceiling; it does
	 * not authorize rewriting retained occurrence cuts.  In particular, a
	 * zoom must begin at the last coherent cut and adjust from there instead
	 * of normalizing every occurrence to a cheap baseline.  Cold coverage and
	 * point-proxy recovery request retained admission explicitly, while a
	 * completed over-budget quiet frame uses overloadRecovery. */
	cursor.retainedAdmissionValue =
	    (requestedRetainedReallocation && costBudget != SIZE_MAX) ||
	    ((inputs.interactive || overloadRecovery ||
	      requestedRetainedRecovery) &&
	     costBudget != SIZE_MAX && inputs.activeCost > costBudget);
	if (overloadRecovery) {
	    this->overloadRecoveryPerformedValue = true;
	    this->overloadRecoveryActiveCostValue = inputs.activeCost;
	}
	/* The retained scene already owns one minimum drawable prefix per active
	 * occurrence.  Bounded recovery windows therefore consume only detail
	 * above that irreducible coverage floor.  Charging the complete budget in
	 * every window drove the whole scene to minimum before a second pass
	 * rebuilt it, producing visible cycling and destroying importance order. */
	cursor.retainedAdmissionRemainingValue =
	    cursor.retainedAdmissionValue ?
		(costBudget > inputs.minimumActiveCost ?
		    costBudget - inputs.minimumActiveCost : 0) : SIZE_MAX;
	cursor.initializedValue = true;
	this->requestedRetainedRecoveryBudgetValue = SIZE_MAX;
	this->retainedAllocationRequestValue.reset();
	this->requestedCoverageCompletionAdditionalCostValue = SIZE_MAX;

	decision.initialized = true;
	decision.overloadRecovery = overloadRecovery;
	decision.totalBudget = this->currentBudgetValue;
	decision.refinementBudget = cursor.refinementRemainingValue;
	decision.retainedAdmission = cursor.retainedAdmissionValue;
	decision.retainedAdmissionBudget =
	    cursor.retainedAdmissionRemainingValue;
	return decision;
    }

    CompletedFrameDecision applyCapacitySearchDecision(
	const BObolLodCapacitySearchCertificate::Decision &search,
	BObolLodAdmissionCursor *cursor)
    {
	CompletedFrameDecision decision;
	decision.searchResult = search.result;
	if (search.requestsFrame()) {
	    this->stableBudgetLimitedValue = false;
	    decision.requestSampleFrame = true;
	    return decision;
	}

	if (search.requestsReallocation()) {
	    this->currentBudgetValue = search.budget;
	    this->retainedAllocationRequestValue.requestReallocation(true);
	    /* The search result names the next candidate population exactly.  The
	     * retained allocator must apply that budget before any throughput EMA,
	     * deadline floor, or growth heuristic may propose another one.  Letting
	     * planPass() recompute it makes the following certificate observe a
	     * different candidate and correctly reject the transition as stale. */
	    this->stableBudgetLimitedValue = false;
	    if (cursor)
		cursor->reset();
	    decision.restartSubmission = true;
	    return decision;
	}

	if (!search.terminal())
	    return decision;
	if (search.result ==
		BObolLodCapacitySearchCertificate::Result::STALE_POPULATION) {
	    this->capacitySearchValue.reset();
	    this->stableBudgetLimitedValue = false;
	    if (cursor)
		cursor->reset();
	    decision.restartSubmission = true;
	    return decision;
	}

	this->stableBudgetLimitedValue = true;
	if (search.result ==
		BObolLodCapacitySearchCertificate::Result::CERTIFIED) {
	    const size_t certifiedBudget = std::max<size_t>(1, search.budget);
	    if (certifiedBudget != this->currentBudgetValue) {
		this->currentBudgetValue = certifiedBudget;
		this->retainedAllocationRequestValue.requestReallocation(true);
		if (cursor)
		    cursor->reset();
		decision.restartSubmission = true;
	    }
	}
	return decision;
    }

public:

    CalibrationDecision finishBlockedPass(
	const CalibrationInputs &inputs)
    {
	CalibrationDecision decision;
	if (inputs.boundedSearch && inputs.searchKey.valid() &&
	    inputs.candidateBudget > 0 &&
	    inputs.allocationCutsApplied && inputs.demandedBudget > 0) {
	    BObolLodCapacitySearchCertificate::Observation observation;
	    observation.key = inputs.searchKey;
	    observation.candidateBudget = inputs.candidateBudget;
	    observation.populationSignature = inputs.populationSignature;
	    /* A changed allocation pass freezes its pre-pass active cost in the
	     * mechanical cursor.  Bind the new candidate only when its completed
	     * frame supplies the post-commit population. */
	    const bool presentationPending = inputs.passAdmittedWork ||
		!inputs.allocationPresentationRealized;
	    observation.presentedCost = presentationPending ? 0 :
		inputs.activeCost;
	    observation.knownSafeBudget = inputs.knownSafeBudget;
	    observation.observedNanoseconds = inputs.observedNanoseconds;
	    observation.validSample = !presentationPending &&
		inputs.observedNanoseconds > 0;
	    const BObolLodCapacitySearchCertificate::Decision search =
		presentationPending ?
		    this->capacitySearchValue.prepare(observation) :
		    this->capacitySearchValue.observe(observation);
	    const CompletedFrameDecision transition =
		this->applyCapacitySearchDecision(search, NULL);
	    decision.searchActive = true;
	    decision.candidateReallocation = search.requestsReallocation();
	    decision.requestFrame = transition.requestSampleFrame;
	    decision.sampleFrame = transition.requestSampleFrame;
	    decision.restartSubmission = transition.restartSubmission;
	    decision.searchTerminal = search.terminal();
	    decision.searchResult = search.result;
	    return decision;
	}
	/* An ordinary per-occurrence pass does not own the complete ordered demand
	 * endpoint required by the capacity-search certificate.  If it changed the
	 * framebuffer, present that population exactly once.  Its successor is a
	 * complete retained importance allocation; if it changed nothing, begin
	 * that allocation immediately.  Replaying three untyped frames here was a
	 * second capacity controller and could alternate with the certified search. */
	decision.searchActive = false;
	decision.candidateReallocation = false;
	decision.sampleFrame = false;
	this->stableBudgetLimitedValue = false;
	if (inputs.passAdmittedWork) {
	    this->frameBarrierValue.requestReallocation();
	    decision.requestFrame = true;
	} else {
	    this->retainedAllocationRequestValue.requestReallocation(false);
	    decision.restartSubmission = true;
	}
	return decision;
    }

    CompletedFrameDecision completeCalibrationFrame(
	BObolLodAdmissionCursor &cursor)
    {
	CompletedFrameDecision decision;
	if (!this->frameBarrierValue.pending())
	    return decision;
	const BObolLodCapacityFrameBarrier::Successor successor =
	    this->frameBarrierValue.consume();
	if (successor ==
		BObolLodCapacityFrameBarrier::Successor::REALLOCATE) {
	    this->retainedAllocationRequestValue.requestReallocation(false);
	}
	cursor.reset();
	decision.restartSubmission = true;
	return decision;
    }

    CompletedFrameDecision completeCapacitySearchFrame(
	BObolLodAdmissionCursor &cursor,
	const CompletedFrameInputs &inputs)
    {
	if (!this->capacitySearchValue.awaitingSample())
	    return this->completeCalibrationFrame(cursor);
	if (inputs.candidateState ==
		CompletedFrameInputs::CandidateState::REALLOCATION_REQUIRED) {
	    CompletedFrameDecision decision;
	    decision.searchResult =
		BObolLodCapacitySearchCertificate::Result::STALE_POPULATION;
	    this->capacitySearchValue.reset();
	    this->stableBudgetLimitedValue = false;
	    cursor.reset();
	    if (!inputs.reallocationProducerPending) {
		this->retainedAllocationRequestValue.requestReallocation(true);
		decision.restartSubmission = true;
	    }
	    return decision;
	}
	if (inputs.candidateState ==
		CompletedFrameInputs::CandidateState::PRESENTATION_PENDING) {
	    /* A completed pass normally transfers this case to the presentation
	     * handoff before the search starts.  If asynchronous state ordering
	     * reaches the frame seam first, repainting cannot remove the ceiling.
	     * Retire the premature search and request one allocation-owner pass;
	     * its immutable owner selection performs the handoff. */
	    CompletedFrameDecision decision;
	    decision.searchResult =
		BObolLodCapacitySearchCertificate::Result::HANDOFF_REQUIRED;
	    this->capacitySearchValue.reset();
	    this->stableBudgetLimitedValue = false;
	    this->retainedAllocationRequestValue.requestReallocation(true);
	    cursor.reset();
	    decision.restartSubmission = true;
	    return decision;
	}
	BObolLodCapacitySearchCertificate::Observation observation;
	observation.key = inputs.searchKey;
	observation.candidateBudget = inputs.candidateBudget;
	observation.presentedCost = inputs.presentedCost;
	observation.populationSignature = inputs.populationSignature;
	observation.knownSafeBudget = inputs.knownSafeBudget;
	observation.observedNanoseconds = inputs.observedNanoseconds;
	observation.validSample = inputs.validSample &&
	    inputs.candidateState ==
		CompletedFrameInputs::CandidateState::CURRENT;
	return this->applyCapacitySearchDecision(
	    this->capacitySearchValue.observe(observation), &cursor);
    }

    /* A calibration barrier can outlive the population it was meant to
     * measure (most notably when a structural-proxy coverage wave temporarily
     * has no active managed mesh cost).  Such a frame is presentation-only:
     * replaying it cannot produce capacity evidence.  Retire the obsolete
     * probe series and return directly to admission instead of waiting on an
     * impossible sample forever. */
    CompletedFrameDecision retireUnmeasurableCalibrationFrame(
	BObolLodAdmissionCursor &cursor)
    {
	CompletedFrameDecision decision;
	if (!this->frameBarrierValue.pending())
	    return decision;
	decision = this->completeCalibrationFrame(cursor);
	return decision;
    }

    void requestRescanAfterFrame(void)
    {
	/* A coverage/publication edge changes the population being presented.
	 * Any numeric candidate awaiting a sample is therefore stale; the new
	 * population must obtain a fresh complete allocation after this barrier. */
	this->capacitySearchValue.reset();
	this->frameBarrierValue.requestReplan();
	this->stableBudgetLimitedValue = false;
    }

    void clearBudgetLimit(void)
    {
	this->stableBudgetLimitedValue = false;
    }

    void resetCalibration(void)
    {
	this->frameBarrierValue.reset();
	this->stableBudgetLimitedValue = false;
	this->capacitySearchValue.reset();
    }

    void resetOverloadRecovery(void)
    {
	this->overloadRecoveryPerformedValue = false;
	this->overloadRecoveryActiveCostValue = 0;
    }

    void reset(void)
    {
	this->currentBudgetValue = this->seedBudgetValue;
	this->requestedRetainedRecoveryBudgetValue = SIZE_MAX;
	this->retainedAllocationRequestValue.reset();
	this->requestedCoverageCompletionAdditionalCostValue = SIZE_MAX;
	this->retainedRecoveryCeilingValue = SIZE_MAX;
	this->steadyDeadlineCapacityCeilingValue = SIZE_MAX;
	this->staticDeadlineCapacityCeilingValue = SIZE_MAX;
	this->retainedQualityFloorBudgetValue = 0;
	this->retainedQualityFloorSignatureValue = 0;
	this->retainedQualityFloorRejectedValue = false;
	this->retainedQualityFloorMissCountValue = 0;
	this->resetOverloadRecovery();
	this->resetCalibration();
    }

    void reduceCurrentBudget(size_t budget)
    {
	this->currentBudgetValue = std::min(
	    this->currentBudgetValue, budget);
    }

    void requestRetainedRecovery(size_t budget)
    {
	if (!budget || budget == SIZE_MAX)
	    return;
	this->requestedRetainedRecoveryBudgetValue = std::min(
	    this->requestedRetainedRecoveryBudgetValue, budget);
	this->retainedRecoveryCeilingValue = std::min(
	    this->retainedRecoveryCeilingValue, budget);
	this->currentBudgetValue = std::min(
	    this->currentBudgetValue, budget);
    }

    /* Normalize a retained population for one admission pass without
     * recording a capacity ceiling.  Motion/cold handoff uses this to obtain
     * a coherent measured baseline; subsequent passes must be free to spend
     * the headroom demonstrated by that baseline. */
    void requestRetainedNormalization(size_t budget)
    {
	if (!budget || budget == SIZE_MAX)
	    return;
	this->requestedRetainedRecoveryBudgetValue = std::min(
	    this->requestedRetainedRecoveryBudgetValue, budget);
	this->currentBudgetValue = std::min(
	    this->currentBudgetValue, budget);
    }

    /* Reallocate an already resident population within the current measured
     * scene allowance.  Unlike overload recovery this is not evidence of a
     * smaller renderer capacity and therefore installs no persistent budget
     * ceiling.  The request is consumed by exactly one complete admission
     * pass; a later view census must explicitly request another one. */
    void requestRetainedReallocation(bool preserveCurrentBudget = true)
    {
	this->retainedAllocationRequestValue.requestReallocation(
	    preserveCurrentBudget);
    }

    /* Convert a completed renderer-wide constrained frame into one atomic
     * occurrence-local importance allocation.  Unlike normalization, this
     * must run even when point aggregation makes the exact presented cost
     * cheaper than the hidden retained prefixes. */
    void requestPresentationReconciliation(size_t budget)
    {
	if (!budget || budget == SIZE_MAX)
	    return;
	this->retainedAllocationRequestValue.
	    requestPresentationReconciliation(budget);
	this->currentBudgetValue = budget;
    }

    /* Admit an irreducible structural frontier once using capacity proven by
     * the exact framebuffer immediately preceding it.  The request does not
     * rewrite retained occurrence cuts and installs no long-lived throughput
     * floor; it is consumed by exactly one structural repair pass.  An older
     * deadline ceiling describes a different occurrence allocation and may
     * not serialize this coverage-first candidate.  The candidate budget is
     * derived from the immediately preceding exact frame and the endpoint
     * hard deadline independently accepts or rejects the resulting batch. */
    size_t requestCoverageCompletion(size_t activeCost,
	size_t certifiedBudget)
    {
	if (certifiedBudget == SIZE_MAX || certifiedBudget <= activeCost)
	    return 0;
	/* The exact frame proves marginal capacity above the population it drew.
	 * Retained presentation metadata may change before the admission plan
	 * freezes its
	 * accounting baseline, so carrying the old absolute total would either
	 * over-admit or strand the complete repair frontier. */
	this->requestedCoverageCompletionAdditionalCostValue =
	    certifiedBudget - activeCost;
	this->currentBudgetValue = certifiedBudget;
	this->stableBudgetLimitedValue = false;
	return certifiedBudget;
    }

    /* Record a strict upper bound after a quiet, capacity-relevant render
     * abort.  Five percent separates the successor from harmless timer jitter
     * and, more importantly, proves monotonic progress if another
     * intermediate trial also misses.  The immediate recovery may choose a
     * more conservative budget; this ceiling governs later calibration once
     * that one-shot handoff is complete. */
    void noteDeadlineCapacityMiss(size_t attemptedBudget,
	bool staticDeadline = false, BObolLodAdmissionCursor *cursor = NULL)
    {
	if (attemptedBudget <= 1 || attemptedBudget == SIZE_MAX)
	    return;
	const bool activeSearch = this->capacitySearchValue.phase() ==
	    BObolLodCapacitySearchCertificate::Phase::MEASURING;
	const bool activeStaticGoal = activeSearch &&
	    this->capacitySearchValue.goal() ==
		BObolLodCapacitySearchCertificate::Goal::STATIC;
	const bool effectiveStaticDeadline = staticDeadline || activeStaticGoal;
	/* Allocation and submitted-render costs may use different currencies.
	 * During a bounded search the allocation candidate is the only sound
	 * bracket coordinate; outside one, the renderer's estimate remains the
	 * best available recovery bound. */
	const size_t failedBudget = activeSearch ?
	    this->capacitySearchValue.candidateBudget() : attemptedBudget;
	const size_t reduction = std::max<size_t>(1, failedBudget / 20);
	const size_t provenSafeBudget = activeSearch ?
	    this->capacitySearchValue.safeBudget() : 0;
	const size_t strictCeiling = std::max(
	    provenSafeBudget, failedBudget - reduction);
	if (effectiveStaticDeadline) {
	    this->staticDeadlineCapacityCeilingValue = std::min(
		this->staticDeadlineCapacityCeilingValue, strictCeiling);
	    /* A population which misses the longer static deadline is also unsafe
	     * at the preferred cadence. */
	    this->steadyDeadlineCapacityCeilingValue = std::min(
		this->steadyDeadlineCapacityCeilingValue, strictCeiling);
	} else {
	    this->steadyDeadlineCapacityCeilingValue = std::min(
		this->steadyDeadlineCapacityCeilingValue, strictCeiling);
	}
	const size_t activeCeiling = effectiveStaticDeadline ?
	    this->staticDeadlineCapacityCeilingValue :
	    this->steadyDeadlineCapacityCeilingValue;
	if (activeSearch) {
	    BObolLodCapacitySearchKey updatedKey =
		this->capacitySearchValue.key();
	    updatedKey.preferredBudgetCeiling =
		this->steadyDeadlineCapacityCeilingValue;
	    updatedKey.maximumBudgetCeiling =
		this->staticDeadlineCapacityCeilingValue;
	    const BObolLodCapacitySearchCertificate::Decision search =
		this->capacitySearchValue.observeDeadlineMiss(updatedKey);
	    (void)this->applyCapacitySearchDecision(search, cursor);
	} else {
	    this->currentBudgetValue = std::min(
		this->currentBudgetValue, activeCeiling);
	}
	/* A protected floor above a hard failed-work bound cannot be defended by
	 * another soft allocation.  Retire it for this capacity epoch while
	 * keeping every immutable resident suffix available to a later view. */
	if (this->retainedQualityFloorBudgetValue >
	    this->staticDeadlineCapacityCeilingValue) {
	    this->retainedQualityFloorBudgetValue = 0;
	    this->retainedQualityFloorSignatureValue = 0;
	    this->retainedQualityFloorMissCountValue = 0;
	    this->retainedQualityFloorRejectedValue = true;
	}
    }

    /* A view-significance floor may deliberately trade some of the 20 Hz
     * stable target for recognizable prominent geometry, but only after the
     * controller proves that population fits its separately bounded hard
     * quality-frame allowance.  Update all three currencies of the already
     * initialized retained pass together; changing only currentBudgetValue
     * would leave its bounded actions consuming the former upgrade limit. */
    void setRetainedQualityFloorBudget(BObolLodAdmissionCursor &cursor,
	size_t budget,
	uint64_t populationSignature, size_t activeCost,
	size_t minimumActiveCost)
    {
	if (this->retainedQualityFloorRejectedValue)
	    return;
	const size_t nextBudget = budget == SIZE_MAX ? 0 : budget;
	if (nextBudget > 0 &&
	    this->staticDeadlineCapacityCeilingValue != SIZE_MAX &&
	    nextBudget > this->staticDeadlineCapacityCeilingValue) {
	    this->retainedQualityFloorBudgetValue = 0;
	    this->retainedQualityFloorSignatureValue = 0;
	    this->retainedQualityFloorMissCountValue = 0;
	    this->retainedQualityFloorRejectedValue = true;
	    return;
	}
	const uint64_t nextSignature = nextBudget ? populationSignature : 0;
	/* Deadline evidence is meaningful only for the exact protected
	 * occurrence/cut population which produced it.  During cold streaming or
	 * a retained reallocation, different floors can be attempted before an
	 * exact frame is presented.  Combining those misses can reject a later,
	 * affordable floor even though it has never missed once. */
	if (this->retainedQualityFloorSignatureValue != nextSignature)
	    this->retainedQualityFloorMissCountValue = 0;
	this->retainedQualityFloorBudgetValue = nextBudget;
	this->retainedQualityFloorSignatureValue = nextSignature;
	if (!cursor.initializedValue || !cursor.retainedAdmissionValue ||
	    !budget || budget == SIZE_MAX || budget <= this->currentBudgetValue)
	    return;
	this->currentBudgetValue = budget;
	cursor.refinementRemainingValue = budget > activeCost ?
	    budget - activeCost : 0;
	cursor.retainedAdmissionRemainingValue =
	    budget > minimumActiveCost ? budget - minimumActiveCost : 0;
	if (this->retainedRecoveryCeilingValue != SIZE_MAX &&
	    this->retainedRecoveryCeilingValue < budget)
	    this->retainedRecoveryCeilingValue = budget;
    }

    void clearRetainedQualityFloorBudget(void)
    {
	this->retainedQualityFloorBudgetValue = 0;
	this->retainedQualityFloorSignatureValue = 0;
	this->retainedQualityFloorRejectedValue = false;
	this->retainedQualityFloorMissCountValue = 0;
    }

    bool noteRetainedQualityFloorMiss(void)
    {
	if (!this->retainedQualityFloorBudgetValue ||
	    this->retainedQualityFloorRejectedValue)
	    return false;
	if (this->retainedQualityFloorMissCountValue < UINT_MAX)
	    this->retainedQualityFloorMissCountValue++;
	/* One cold upload, command-record transition, compositor stall, or
	 * software JIT event is not renderer-capacity evidence.  Require three
	 * independent ceiling-free attempts at the same view-local floor before
	 * allowing the soft allocator to sacrifice protected geometry. */
	if (this->retainedQualityFloorMissCountValue < 3)
	    return false;
	this->retainedQualityFloorBudgetValue = 0;
	this->retainedQualityFloorSignatureValue = 0;
	this->retainedQualityFloorRejectedValue = true;
	return true;
    }

    /* A protected floor is allowed to exceed the preferred quiet cadence
     * only while it remains inside the independently interruptible static
     * frame deadline.  A miss from the explicit static-quality trial is not
     * one of the transient setup samples counted above: that trial is the
     * endpoint contract itself.  Retire the floor for this policy epoch so a
     * recovery allocation cannot immediately restore the population which
     * just failed.  Resident PoP data and occurrence demand are untouched. */
    bool rejectRetainedQualityFloor(void)
    {
	const bool hadFloor = this->retainedQualityFloorBudgetValue > 0;
	this->retainedQualityFloorBudgetValue = 0;
	this->retainedQualityFloorSignatureValue = 0;
	this->retainedQualityFloorMissCountValue = 0;
	this->retainedQualityFloorRejectedValue = true;
	return hadFloor;
    }

    bool noteRetainedQualityFloorMet(bool exactProtectedPopulation,
	uint64_t populationSignature, size_t presentedCost)
    {
	if (!exactProtectedPopulation ||
	    !this->retainedQualityFloorBudgetValue ||
	    this->retainedQualityFloorRejectedValue ||
	    populationSignature != this->retainedQualityFloorSignatureValue ||
	    presentedCost < this->retainedQualityFloorBudgetValue)
	    return false;
	this->retainedQualityFloorMissCountValue = 0;
	return true;
    }

    bool retainedQualityFloorRejected(void) const
    {
	return this->retainedQualityFloorRejectedValue;
    }

    bool retainedQualityFloorActive(void) const
    {
	return this->retainedQualityFloorBudgetValue > 0;
    }

    size_t retainedQualityFloorBudget(void) const
    {
	return this->retainedQualityFloorBudgetValue;
    }

    uint64_t retainedQualityFloorSignature(void) const
    {
	return this->retainedQualityFloorSignatureValue;
    }

    unsigned int retainedQualityFloorMissCount(void) const
    {
	return this->retainedQualityFloorMissCountValue;
    }

    void clearRetainedRecoveryCeiling(void)
    {
	this->requestedRetainedRecoveryBudgetValue = SIZE_MAX;
	this->retainedRecoveryCeilingValue = SIZE_MAX;
    }

    void clearDeadlineCapacityCeiling(void)
    {
	this->steadyDeadlineCapacityCeilingValue = SIZE_MAX;
	this->staticDeadlineCapacityCeilingValue = SIZE_MAX;
    }

    /* A measured recovery ceiling protects the first coherent one-pixel
     * population from immediately re-admitting the cut which just missed its
     * deadline.  It must end once that population is actually ready for
     * presentation.  Keep this transition in the scalar policy so callers
     * cannot accidentally clear only the pass state, or retain the ceiling
     * forever when returning from a coarser point cut needs an extra frame. */
    bool confirmRetainedRecoveryPresentation(
	bool onePixelReady, BObolLodAdmissionCursor &cursor)
    {
	if (!onePixelReady ||
	    this->retainedRecoveryCeilingValue == SIZE_MAX)
	    return false;
	this->clearRetainedRecoveryCeiling();
	cursor.reset();
	return true;
    }

    void raiseCurrentBudget(size_t budget)
    {
	size_t raisedBudget = std::max(this->currentBudgetValue, budget);
	if (this->retainedRecoveryCeilingValue != SIZE_MAX)
	    raisedBudget = std::min(
		raisedBudget, this->retainedRecoveryCeilingValue);
	if (this->steadyDeadlineCapacityCeilingValue != SIZE_MAX)
	    raisedBudget = std::min(
		raisedBudget, this->steadyDeadlineCapacityCeilingValue);
	this->currentBudgetValue = raisedBudget;
    }

    size_t seedBudget(void) const { return this->seedBudgetValue; }
    static constexpr size_t singleOccurrenceBootstrapBudget(void)
    {
	/* A sole cold source has no competing occurrence whose first useful mesh
	 * could be starved.  Permit one globally classified, byte-capped PoP
	 * preview rather than rejecting it for a provisional 500k estimate and
	 * making the user wait for chunk persistence.  The endpoint's 100 ms hard
	 * deadline and completed-frame calibration remain the safety contract;
	 * this is only the bounded first-publication allowance. */
	return 1000000;
    }
    size_t currentBudget(void) const { return this->currentBudgetValue; }
    bool retainedRecoveryCeilingActive(void) const
    {
	return this->retainedRecoveryCeilingValue != SIZE_MAX;
    }
    size_t deadlineCapacityCeiling(void) const
    {
	return this->steadyDeadlineCapacityCeilingValue;
    }
    size_t staticDeadlineCapacityCeiling(void) const
    {
	return this->staticDeadlineCapacityCeilingValue;
    }
    bool overloadRecoveryPerformed(void) const
    {
	return this->overloadRecoveryPerformedValue;
    }
    size_t overloadRecoveryActiveCost(void) const
    {
	return this->overloadRecoveryActiveCostValue;
    }
    bool rescanAfterFrame(void) const
    {
	return this->frameBarrierValue.pending() ||
	    this->capacitySearchValue.awaitingSample();
    }
    bool stableBudgetLimited(void) const
    {
	return this->stableBudgetLimitedValue;
    }
    const BObolLodCapacitySearchCertificate &capacitySearch(void) const
    {
	return this->capacitySearchValue;
    }

private:
    size_t seedBudgetValue = 50000;
    size_t currentBudgetValue = 50000;
    size_t overloadRecoveryActiveCostValue = 0;
    size_t requestedRetainedRecoveryBudgetValue = SIZE_MAX;
    size_t retainedRecoveryCeilingValue = SIZE_MAX;
    size_t steadyDeadlineCapacityCeilingValue = SIZE_MAX;
    size_t staticDeadlineCapacityCeilingValue = SIZE_MAX;
    size_t retainedQualityFloorBudgetValue = 0;
    uint64_t retainedQualityFloorSignatureValue = 0;
    unsigned int retainedQualityFloorMissCountValue = 0;
    BObolLodRetainedAllocationRequest retainedAllocationRequestValue;
    size_t requestedCoverageCompletionAdditionalCostValue = SIZE_MAX;
    bool overloadRecoveryPerformedValue = false;
    bool stableBudgetLimitedValue = false;
    bool retainedQualityFloorRejectedValue = false;
    BObolLodCapacityFrameBarrier frameBarrierValue;
    BObolLodCapacitySearchCertificate capacitySearchValue;
};

static_assert(std::is_trivially_copyable<BObolLodCapacityEvidence>::value,
    "capacity evidence must remain an allocation-free value");

/*
 * Camera-local view-demand scheduling.  A zoom interaction may keep loading
 * the pixel-demanded immutable PoP suffix while exposing only one additional
 * coherent population per completed motion frame.  This allocation-free
 * policy owns that bounded-probe lifecycle and the proof that a coverage pass
 * is a scale-demand refresh rather than a cold scene restart.  Renderer
 * ceilings and scene measurements remain controller actions/inputs.
 */

#endif /* LIBBOBOL_LOD_CAPACITY_POLICY_PRIVATE_H */
