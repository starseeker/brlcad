/* L O D _ P R E S E N T A T I O N _ P O L I C Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "lod_presentation_policy_private.h"

namespace {

constexpr float MinimumPixelError = 0.25f;
constexpr float MaximumPixelError = 64.0f;
constexpr float DefaultPixelError = 1.0f;
constexpr float MinimumPointProxyThreshold = 1.0f;
constexpr float MaximumPointProxyThreshold = 64.0f;
constexpr float DefaultPointProxyThreshold = 1.0f;

} // namespace

bool
BObolLodPresentationPolicy::presentationLimitsReconciled(
    bool allocationCompleted, bool allocationCertified,
    size_t selectedPresentationCost, size_t certifiedPresentationBudget,
    size_t reconciliationBudget)
{
    if (!allocationCompleted || !allocationCertified ||
	    !certifiedPresentationBudget)
	return false;
    const size_t budget = reconciliationBudget > 0 ?
	std::min(certifiedPresentationBudget, reconciliationBudget) :
	certifiedPresentationBudget;
    return selectedPresentationCost <= budget;
}

bool
BObolLodPresentationPolicy::capacitySamplePending(
    bool searchAwaitingSample, bool allocationCurrent,
    bool allocationPresentationRealized)
{
    return searchAwaitingSample ||
	(allocationCurrent && !allocationPresentationRealized);
}

bool
BObolLodPresentationPolicy::nonterminalStructuralFrontier(
    bool occurrenceCoverageExact, size_t structuralFallbackCount,
    size_t terminalFailureCount)
{
    return occurrenceCoverageExact &&
	structuralFallbackCount > terminalFailureCount;
}

bool
BObolLodPresentationPolicy::presentationOnlyFrameAdvancesPlanning(
    bool lodPlanningRelevant, bool pointAdmissionFramePending,
    bool lodPresentationBarrierPending,
    bool nonterminalStructuralFrontierPending)
{
    return lodPlanningRelevant && (pointAdmissionFramePending ||
	lodPresentationBarrierPending || nonterminalStructuralFrontierPending);
}

bool
BObolLodPresentationPolicy::claimOverBudgetAllocation(
    bool allocationCurrent, bool allocationCutsApplied,
    size_t selectedPresentationCost, size_t certifiedPresentationBudget,
    size_t presentedRenderCost)
{
    if (!allocationCurrent || !allocationCutsApplied ||
	    !selectedPresentationCost || !certifiedPresentationBudget ||
	    selectedPresentationCost <= certifiedPresentationBudget)
	return false;
    if (!this->handoffActive())
	this->armHandoff(false, presentedRenderCost,
	    certifiedPresentationBudget);
    return true;
}

BObolLodPresentationPolicy::CompletedPassSelection
BObolLodPresentationPolicy::completedPassSelection(
    const CompletedPassInputs &inputs, bool capacityEligible,
    bool capacitySuccessorPending, bool capacityPresentationPending,
    int progressiveCeiling) const
{
    CompletedPassSelection selection;
    if (!inputs.completed || inputs.submissionPending)
	return selection;
    selection.consumePassAnnotations =
	this->capacitySampleRequiresCeilingFreeHandoff(inputs,
	    capacityPresentationPending, progressiveCeiling);
    selection.deferredCapacitySuccessor =
	selection.consumePassAnnotations && capacitySuccessorPending;

    /* Structural-only occurrences cannot participate in retained allocation
     * or capacity measurement.  Their exact projected frontier therefore
     * owns the completed-pass boundary before either a renderer-ceiling
     * handoff or a new capacity transaction.  The controller invalidates any
     * stale capacity sample and routes this owner to its finite classifier,
     * mesh-repair, or terminal-proxy successor. */
    if (inputs.structuralFrontierPending) {
	selection.consumePassAnnotations = false;
	selection.deferredCapacitySuccessor = false;
	selection.owner = CompletedPassOwner::STRUCTURAL_FRONTIER;
	return selection;
    }

    const bool cleanPass =
	(!inputs.capacityTransactionPending && !inputs.changedCut) ||
	selection.consumePassAnnotations;
    if (cleanPass && this->handoffActive() &&
	(!capacitySuccessorPending || selection.consumePassAnnotations)) {
	selection.owner = this->presentationHandoffPending() ?
	    CompletedPassOwner::PRESENTATION_HANDOFF :
	    CompletedPassOwner::ALLOCATION_HANDOFF;
	return selection;
    }
    if (capacityEligible) {
	selection.owner = CompletedPassOwner::CAPACITY;
	return selection;
    }
    /* Wanting a richer retained cut is an observation, not a producer.  Give
     * its bounded retirement an explicit owner so the completed-pass
     * transaction cannot become terminal while the planning ledger still
     * contains that observation. */
    if (cleanPass && !this->handoffActive() &&
	inputs.retainedRefinementPending &&
	!inputs.retainedRefinementBudgetBlocked)
	selection.owner = CompletedPassOwner::RETAINED_OBSERVATION;
    return selection;
}

bool
BObolLodPresentationPolicy::capacityProducerRequiredAfterHandoff(
    bool deferredCapacitySuccessor, bool handoffFinished,
    size_t progressivePayloadCount)
{
    return deferredCapacitySuccessor && handoffFinished &&
	progressivePayloadCount > 1;
}

bool
BObolLodPresentationPolicy::capacitySampleRequiresCeilingFreeHandoff(
    const CompletedPassInputs &inputs, bool capacityPresentationPending,
    int progressiveCeiling) const
{
    return capacityPresentationPending && progressiveCeiling >= 0 &&
	(inputs.capacityTransactionPending || inputs.changedCut) && inputs.completed &&
	!inputs.submissionPending && this->handoffActive() &&
	!this->presentationHandoffPending() &&
	inputs.retainedAllocationCompleted &&
	inputs.retainedAllocationCertified && inputs.populationQuiescent &&
	inputs.presentationLimitsReconciled;
}

void
BObolLodPresentationPolicy::capturePrior(
    float targetPixelError, int progressiveCeiling,
    float progressiveNextFraction, float pointProxyPixelThreshold,
    const Population &population, BObolLodViewEpoch viewEpoch)
{
    this->priorStableValue = makeSnapshot(
	targetPixelError, progressiveCeiling, progressiveNextFraction,
	pointProxyPixelThreshold, population, viewEpoch, 0);
    this->provenQualityValue = Snapshot();
    this->clearHandoff();
}

void
BObolLodPresentationPolicy::noteStableQuality(
    float targetPixelError, int progressiveCeiling,
    float progressiveNextFraction, float pointProxyPixelThreshold,
    const Population &population, BObolLodViewEpoch viewEpoch,
    size_t presentedRenderCost)
{
    if (!presentedRenderCost)
	return;
    const Snapshot candidate = makeSnapshot(
	targetPixelError, progressiveCeiling, progressiveNextFraction,
	pointProxyPixelThreshold, population, viewEpoch, presentedRenderCost);
    if (this->provenQualityValue.valid &&
	    this->provenQualityValue.viewEpoch == viewEpoch &&
	    populationMatches(this->provenQualityValue.population, population) &&
	    this->provenQualityValue.presentedRenderCost > presentedRenderCost)
	return;
    this->provenQualityValue = candidate;
}

BObolLodPresentationPolicy::RestoreDecision
BObolLodPresentationPolicy::beginQuiet(const QuietInputs &quietInputs)
{
    const RestoreInputs &inputs = quietInputs.presentation;
    if (quietInputs.unboundedTerminal) {
	RestoreDecision decision;
	decision.apply = true;
	decision.clearPresentationLimits = true;
	decision.targetPixelError = inputs.currentTargetPixelError;
	decision.progressiveCeiling = -1;
	decision.progressiveNextFraction = 0.0f;
	decision.pointProxyPixelThreshold = DefaultPointProxyThreshold;
	this->clearHandoff();
	return decision;
    }

    using Reducer = BObolLodQuietSuccessorReducer;
    Reducer::Inputs successorInputs;
    successorInputs.retainedMeshPayloads = inputs.retainedMeshPayloads;
    successorInputs.current.valid = true;
    successorInputs.current.targetPixelError = inputs.currentTargetPixelError;
    successorInputs.current.progressiveCeiling =
	inputs.currentProgressiveCeiling;
    successorInputs.current.progressiveNextFraction =
	inputs.currentProgressiveNextFraction;
    successorInputs.current.pointProxyPixelThreshold =
	inputs.currentPointProxyPixelThreshold;

    const bool priorControlsAvailable = inputs.orthographic &&
	this->priorStableValue.valid;
    const bool priorPopulationMatches = priorControlsAvailable &&
	populationMatches(this->priorStableValue.population, inputs.population);
    if (priorControlsAvailable) {
	successorInputs.priorStable = targetFromSnapshot(
	    this->priorStableValue);
	/* Resident publication may enlarge the drawable population during an
	 * interaction without changing the stable pixel/proxy targets captured
	 * at its start.  The old renderer ceiling is not transferable to that
	 * population, however.  Keep the current motion ceiling as a private
	 * presentation guard while the allocator recomputes demand from the
	 * restored semantic controls. */
	if (!priorPopulationMatches) {
	    successorInputs.priorStable.progressiveCeiling =
		inputs.currentProgressiveCeiling;
	    successorInputs.priorStable.progressiveNextFraction =
		inputs.currentProgressiveNextFraction;
	}
	successorInputs.priorStableCutsReusable =
	    priorPopulationMatches && !inputs.scaleChanged;
    }
    const bool provenMatches = inputs.scaleChanged &&
	this->provenQualityValue.valid &&
	this->provenQualityValue.viewEpoch == inputs.viewEpoch &&
	populationMatches(this->provenQualityValue.population,
	    inputs.population);
    if (provenMatches)
	successorInputs.provenScale = targetFromSnapshot(
	    this->provenQualityValue);
    if (quietInputs.exactView.valid) {
	successorInputs.exactView.valid = true;
	successorInputs.exactView.targetPixelError =
	    quietInputs.exactView.targetPixelError;
	successorInputs.exactView.progressiveCeiling =
	    quietInputs.exactView.progressiveCeiling;
	successorInputs.exactView.progressiveNextFraction =
	    quietInputs.exactView.progressiveNextFraction;
	successorInputs.exactView.pointProxyPixelThreshold =
	    quietInputs.exactView.pointProxyPixelThreshold;
	successorInputs.exactView.presentedRenderCost =
	    quietInputs.exactView.presentedRenderCost;
	successorInputs.exactView.provenRenderCostCapacity =
	    quietInputs.exactView.provenRenderCostCapacity;
    }

    const Reducer::Decision successor = Reducer::reduce(successorInputs);
    RestoreDecision decision;
    decision.apply = successor.apply;
    decision.restoredPriorStable =
	successor.source == Reducer::Source::PRIOR_STABLE;
    decision.restoredProvenQuality =
	successor.source == Reducer::Source::PROVEN_SCALE;
    decision.restoredExactView =
	successor.source == Reducer::Source::EXACT_VIEW;
    decision.handoff = successor.handoff;
    decision.clearPresentationLimits = successor.clearPresentationLimits;
    decision.targetPixelError = successor.target.targetPixelError;
    decision.progressiveCeiling = successor.target.progressiveCeiling;
    decision.progressiveNextFraction =
	successor.target.progressiveNextFraction;
    decision.pointProxyPixelThreshold =
	successor.target.pointProxyPixelThreshold;
    decision.provenRenderCostFloor = successor.target.presentedRenderCost;
    decision.provenRenderCostCapacity =
	successor.target.provenRenderCostCapacity;

    this->handoffStateValue =
	decision.handoff == Reducer::Handoff::PRESENTATION ?
	    HandoffState::PRESENTATION_REQUIRED :
	    (decision.handoff == Reducer::Handoff::ALLOCATION ?
		HandoffState::ALLOCATION_REQUIRED : HandoffState::INACTIVE);
    /* A proven scale target carries the exact cost of the completed frame it
     * restores.  Translate its renderer-wide (and potentially spatial-page)
     * ceiling into occurrence-local cuts at that exact allowance before
     * removing the ceiling.  Treating this witness only as a lower bound lets
     * the ordinary static allowance select a much richer global population;
     * a subsequent deadline miss can then discard the useful tiled frame and
     * fall all the way to the minimum cut.  A prior stable control vector is
     * only a starting point after zoom and carries no current-projection cost
     * proof.  Other handoffs establish their reconciliation budget through
     * completed-frame evidence. */
    const size_t provenScaleBudget =
        decision.handoff == Reducer::Handoff::ALLOCATION &&
        decision.restoredProvenQuality ? decision.provenRenderCostFloor : 0;
    this->handoffReconciliationBudgetValue = provenScaleBudget;
    this->handoffReconciliationBudgetLimitValue = provenScaleBudget;
    this->handoffCostFloorValue = provenScaleBudget;

    this->priorStableValue = Snapshot();
    this->provenQualityValue = Snapshot();
    return decision;
}

BObolLodPresentationPolicy::RestoreDecision
BObolLodPresentationPolicy::beginQuiet(const RestoreInputs &inputs)
{
    QuietInputs quietInputs;
    quietInputs.presentation = inputs;
    return this->beginQuiet(quietInputs);
}

BObolLodPresentationPolicy::CompletedPassDecision
BObolLodPresentationPolicy::completePass(const CompletedPassInputs &inputs)
{
    CompletedPassDecision decision;
    if (!inputs.completed || inputs.submissionPending)
	return decision;
    /* Structural repair precedes an allocation handoff.  The retained
     * allocator operates on published CAD occurrences, so asking it to
     * replace structural-only records produces an empty allocation and can
     * re-arm the same handoff indefinitely.  Its exact frontier also owns the
     * completed-pass boundary ahead of capacity and budget evidence, so the
     * ordinary pass's richer-cut observation must not block that successor. */
    if (inputs.structuralFrontierPending) {
	decision.retireRetainedObservation =
	    inputs.retainedRefinementPending;
	return decision;
    }
    if (inputs.capacityTransactionPending || inputs.changedCut ||
	    (inputs.retainedRefinementBudgetBlocked && !this->handoffActive()) ||
	    this->presentationHandoffPending())
	return decision;
    if (!this->handoffActive()) {
	/* Wanting a richer retained cut is an observation, not a liveness
	 * witness.  With no admitted cut, budget barrier, handoff, or rescan,
	 * the current calibrated cut is terminal until a new external edge
	 * starts another pass. */
	decision.retireRetainedObservation =
	    inputs.retainedRefinementPending;
	return decision;
    }
    if (!inputs.retainedAllocationCompleted ||
	    !inputs.retainedAllocationCertified) {
	decision.requestRetainedAllocation = inputs.populationQuiescent;
	return decision;
    }
    if (!inputs.populationQuiescent)
	return decision;
    if (!inputs.presentationLimitsReconciled) {
	decision.requestLocalPresentationReduction = true;
	return decision;
    }
    /* A hard-deadline ceiling may be removed only by an explicit complete-
     * population reconciliation at the constrained frame's proven budget.
     * Merely changing occurrence cuts does not prove that late resident work
     * or a larger retained budget still fits. */
    this->clearHandoff();
    decision.finishHandoff = true;
    decision.requestRetainedRescan = inputs.retainedRefinementPending &&
	!inputs.retainedRefinementBudgetBlocked;
    return decision;
}

void
BObolLodPresentationPolicy::armHandoff(
    bool presentationRequired, size_t provenRenderCost,
    size_t reconciliationBudgetLimit)
{
    this->handoffStateValue = presentationRequired ?
	HandoffState::PRESENTATION_REQUIRED : HandoffState::ALLOCATION_REQUIRED;
    /* A deadline recovery waits for its corrected constrained frame.  A
     * static predictor already has that frame and may expose its budget
     * immediately to the occurrence allocator. */
    this->handoffReconciliationBudgetValue = presentationRequired ? 0 :
	reconciliationBudgetLimit;
    this->handoffReconciliationBudgetLimitValue =
	reconciliationBudgetLimit;
    this->handoffCostFloorValue = provenRenderCost;
}

bool
BObolLodPresentationPolicy::noteFramePresented(
    size_t provenRenderCost, size_t reconciliationBudget)
{
    const bool released = this->handoffStateValue ==
	HandoffState::PRESENTATION_REQUIRED;
    if (!this->handoffActive())
	return released;

    if (released)
	this->handoffStateValue = HandoffState::ALLOCATION_REQUIRED;
    if (provenRenderCost)
	this->handoffCostFloorValue = std::max(
	    this->handoffCostFloorValue, provenRenderCost);
    if (!released)
	return released;

    /* A corrected frame can land just beyond the endpoint timer.  In that
     * case it supplies no new safe extrapolation, but it must retain the
     * recovery limit derived from the interrupted richer population. */
    size_t boundedBudget = reconciliationBudget;
    if (this->handoffReconciliationBudgetLimitValue > 0)
	boundedBudget = boundedBudget > 0 ?
	    std::min(boundedBudget,
		this->handoffReconciliationBudgetLimitValue) :
	    this->handoffReconciliationBudgetLimitValue;
    if (boundedBudget)
	this->handoffReconciliationBudgetValue = boundedBudget;
    return released;
}

void
BObolLodPresentationPolicy::acceptCapacityCertificate(void)
{
    this->clearHandoff();
}

void
BObolLodPresentationPolicy::cancelHandoff(void)
{
    this->clearHandoff();
}

void
BObolLodPresentationPolicy::viewInvalidated(void)
{
    this->clearHandoff();
    this->provenQualityValue = Snapshot();
}

void
BObolLodPresentationPolicy::reset(void)
{
    this->priorStableValue = Snapshot();
    this->provenQualityValue = Snapshot();
    this->clearHandoff();
}

void
BObolLodPresentationPolicy::clearHandoff(void)
{
    this->handoffStateValue = HandoffState::INACTIVE;
    this->handoffReconciliationBudgetValue = 0;
    this->handoffReconciliationBudgetLimitValue = 0;
    this->handoffCostFloorValue = 0;
}

float
BObolLodPresentationPolicy::sanitizePixelError(float pixelError)
{
    return std::isfinite(pixelError) ?
	std::max(MinimumPixelError, std::min(MaximumPixelError, pixelError)) :
	DefaultPixelError;
}

float
BObolLodPresentationPolicy::sanitizeThreshold(float threshold)
{
    return std::isfinite(threshold) ?
	std::max(MinimumPointProxyThreshold,
	    std::min(MaximumPointProxyThreshold, threshold)) :
	DefaultPointProxyThreshold;
}

float
BObolLodPresentationPolicy::sanitizeFraction(
    int progressiveCeiling, float fraction)
{
    if (progressiveCeiling < 0 || !std::isfinite(fraction))
	return 0.0f;
    return std::max(0.0f, std::min(1.0f, fraction));
}

BObolLodPresentationPolicy::Snapshot
BObolLodPresentationPolicy::makeSnapshot(
    float targetPixelError, int progressiveCeiling,
    float progressiveNextFraction, float pointProxyPixelThreshold,
    const Population &population, BObolLodViewEpoch viewEpoch,
    size_t presentedRenderCost)
{
    Snapshot snapshot;
    snapshot.valid = true;
    snapshot.viewEpoch = viewEpoch;
    snapshot.targetPixelError = sanitizePixelError(targetPixelError);
    snapshot.progressiveCeiling = progressiveCeiling < -1 ? -1 :
	std::min<int>(BOBOL_MESH_LOD_CUT_COUNT_MAX - 1,
	    progressiveCeiling);
    snapshot.progressiveNextFraction = sanitizeFraction(
	snapshot.progressiveCeiling, progressiveNextFraction);
    snapshot.pointProxyPixelThreshold =
	sanitizeThreshold(pointProxyPixelThreshold);
    snapshot.population = population;
    snapshot.presentedRenderCost = presentedRenderCost;
    return snapshot;
}

BObolLodQuietSuccessorReducer::Target
BObolLodPresentationPolicy::targetFromSnapshot(const Snapshot &snapshot)
{
    BObolLodQuietSuccessorReducer::Target target;
    target.valid = snapshot.valid;
    target.targetPixelError = snapshot.targetPixelError;
    target.progressiveCeiling = snapshot.progressiveCeiling;
    target.progressiveNextFraction = snapshot.progressiveNextFraction;
    target.pointProxyPixelThreshold = snapshot.pointProxyPixelThreshold;
    target.presentedRenderCost = snapshot.presentedRenderCost;
    return target;
}

bool
BObolLodPresentationPolicy::populationMatches(
    const Population &snapshot, const Population &current)
{
    return snapshot.available == current.available &&
	snapshot.sceneDomainRevision == current.sceneDomainRevision;
}
