/*          L O D _ A D M I S S I O N _ P O L I C Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "lod_admission_policy_private.h"

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::applyEvidenceAction(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, EvidenceAction action)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	switch (action) {
	    case EvidenceAction::RESET_CAPACITY:
		result.nextEvidence.capacityValue.reset();
		result.nextCursor.reset();
		break;
	    case EvidenceAction::CLEAR_CAPACITY_LIMIT:
		result.nextEvidence.capacityValue.clearBudgetLimit();
		break;
	    case EvidenceAction::INVALIDATE_CAPACITY_MEASUREMENT:
		result.nextEvidence.capacityValue.invalidateCalibration();
		break;
	    case EvidenceAction::RESET_CAPACITY_OVERLOAD:
		result.nextEvidence.capacityValue.resetOverloadRecovery();
		break;
	    case EvidenceAction::CLEAR_RETAINED_QUALITY_FLOOR:
		result.nextEvidence.capacityValue.clearRetainedQualityFloorBudget();
		break;
	    case EvidenceAction::CLEAR_RETAINED_RECOVERY_CEILING:
		result.nextEvidence.capacityValue.clearRetainedRecoveryCeiling();
		break;
	    case EvidenceAction::CLEAR_DEADLINE_CAPACITY_CEILING:
		result.nextEvidence.capacityValue.clearDeadlineCapacityCeiling();
		break;
	    case EvidenceAction::CANCEL_HEADROOM_RETRY:
		result.nextEvidence.headroomValue.cancelRetry();
		break;
	case EvidenceAction::RESET_POINT_PROXY:
		result.nextEvidence.responsivePointProxyValue.reset();
		result.nextEvidence.staticPointProxyValue.reset();
		break;
	    case EvidenceAction::RESET_STRUCTURAL_ADMISSION:
		result.nextEvidence.structuralValue.reset();
		break;
	    case EvidenceAction::RESET_RESOURCE_SERVICE:
		result.nextEvidence.resourcesValue.resetForServiceChange();
		break;
	    case EvidenceAction::MARK_RESOURCE_RECOVERY_HANDLED:
		result.nextEvidence.resourcesValue.markRecoveryHandled();
		break;
	}
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::completeAppliedAllocation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodCapacityEvidence::CompletedAllocationInputs &inputs)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.capacityValue.completeAppliedAllocation(inputs);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::requestCapacityRescan(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.capacityValue.requestPopulationBarrier();
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::acceptStaticPresentationConstraint(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.capacityValue.acceptStaticPresentationConstraint(budget);
	result.nextCursor.reset();
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::requestRetainedRecovery(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.capacityValue.requestRetainedRecovery(budget);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::requestRetainedNormalization(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.capacityValue.requestRetainedNormalization(budget);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::requestRetainedReallocation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	bool preserveCurrentBudget)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	const bool requestAccepted = result.nextEvidence.capacityValue.
	    requestRetainedReallocation(preserveCurrentBudget);
	/* A retained reallocation is a request for a new complete population
	 * traversal.  Its evidence latch is consumed only while initializing the
	 * admission cursor, so carrying an initialized cursor from the preceding
	 * pass makes the request permanently pending but non-executable.  Reset the
	 * mechanical cursor whenever this request is accepted.  A weaker request
	 * cannot displace an already pending presentation reconciliation and must
	 * leave that transaction untouched. */
	if (requestAccepted)
	    result.nextCursor.reset();
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::resumeCapacityCandidateAllocation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.transitionChanged = result.nextEvidence.capacityValue.
	    resumeCapacityCandidateAllocation();
	if (result.transitionChanged)
	    result.nextCursor.reset();
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::requestPresentationReconciliation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.transitionChanged = result.nextEvidence.capacityValue.
	    requestPresentationReconciliation(budget);
	if (result.transitionChanged)
	    result.nextCursor.reset();
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::requestCoverageCompletion(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t activeCost,
	size_t certifiedBudget)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.transitionValue =
	    result.nextEvidence.capacityValue.requestCoverageCompletion(
		activeCost, certifiedBudget);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::recordDeadlineCapacityMiss(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodCapacityEvidence::DeadlineMissInputs &inputs)
{
    BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
    result.completedFrameDecision = result.nextEvidence.capacityValue.
	noteDeadlineCapacityMiss(inputs, &result.nextCursor);
    return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::setRetainedQualityFloor(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget,
	uint64_t populationIdentity, size_t activeCost,
	size_t minimumActiveCost)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.capacityValue.setRetainedQualityFloorBudget(
	    result.nextCursor, budget, populationIdentity, activeCost,
	    minimumActiveCost);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::recordRetainedQualityFloorMiss(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.transitionChanged =
	    result.nextEvidence.capacityValue.noteRetainedQualityFloorMiss();
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::recordRetainedQualityFloorMet(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, bool exactProtectedPopulation,
	uint64_t populationIdentity, size_t presentedCost)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.transitionChanged =
	    result.nextEvidence.capacityValue.noteRetainedQualityFloorMet(
		exactProtectedPopulation, populationIdentity, presentedCost);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::confirmRetainedRecoveryPresentation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, bool onePixelReady)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.transitionChanged = result.nextEvidence.capacityValue.
	    confirmRetainedRecoveryPresentation(
		onePixelReady, result.nextCursor);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::raiseCapacityBudget(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.capacityValue.raiseCurrentBudget(budget);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::finishBlockedCapacityPass(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodCapacityEvidence::CalibrationInputs &inputs)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.calibrationDecision =
	    result.nextEvidence.capacityValue.finishBlockedPass(inputs);
	if (result.calibrationDecision.restartSubmission)
	    result.nextCursor.reset();
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::completeCapacityCalibrationFrame(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.completedFrameDecision = result.nextEvidence.capacityValue.
	    completeCalibrationFrame(result.nextCursor);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::completeCapacitySearchFrame(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodCapacityEvidence::CompletedFrameInputs &inputs)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.completedFrameDecision = result.nextEvidence.capacityValue.
	    completeCapacitySearchFrame(result.nextCursor, inputs);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::retireUnmeasurableCapacityFrame(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.completedFrameDecision = result.nextEvidence.capacityValue.
	    retireUnmeasurableCalibrationFrame(result.nextCursor);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::plan(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodAdmissionRevisionStamp &revisionStamp,
	const BObolLodCapacityEvidence::Inputs &inputs)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.revisionStamp = revisionStamp;
	if (!result.nextCursor.matches(revisionStamp))
	    result.nextCursor.resetForRevision(revisionStamp);
	result.capacityDecision =
	    result.nextEvidence.capacityValue.planPass(inputs, result.nextCursor);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::planResourceObservation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, bool cpuPressure,
	bool gpuPressure, bool recoveryEnabled)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.resourceDecision = result.nextEvidence.resourcesValue.observe(
	    cpuPressure, gpuPressure, recoveryEnabled);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::planHeadroomRetry(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
	size_t activePopulationCost, size_t currentBudget)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.headroomAccepted = result.nextEvidence.headroomValue.armRetry(
	    viewEpoch, policyEpoch, activePopulationCost, currentBudget);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::planHeadroomTransientDeferral(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
	size_t activePopulationCost)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.headroomAccepted =
	    result.nextEvidence.headroomValue.deferTransientReplay(
		viewEpoch, policyEpoch, activePopulationCost);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::planHeadroomConsumption(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
	size_t currentBudget, size_t activePopulationCost,
	long double demonstratedBudget, uint64_t elapsedNanoseconds,
	uint64_t targetNanoseconds, bool reusableSample)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.headroomAccepted = result.nextEvidence.headroomValue.consumeRetry(
	    viewEpoch, policyEpoch, currentBudget, activePopulationCost,
	    demonstratedBudget, elapsedNanoseconds, targetNanoseconds,
	    reusableSample);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::recordHeadroomRejection(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
	size_t activePopulationCost, size_t attemptedBudget)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.headroomValue.rejectRetry(viewEpoch, policyEpoch,
	    activePopulationCost, attemptedBudget);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::planPointStructuralDistribution(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, float currentThreshold,
	const std::array<size_t, 7> &cumulativeCount, size_t visibleCount,
	size_t maximumUncollapsedCount)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.pointProxyDecision =
	    result.nextEvidence.responsivePointProxyValue.
		seedFromStructuralDistribution(
		currentThreshold, cumulativeCount, visibleCount,
		maximumUncollapsedCount);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::planPointInterrupted(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodPointCalibrationGoal goal, float currentThreshold,
	uint64_t renderNanoseconds, float targetFps)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	BObolLodPointProxyEvidence &pointEvidence =
	    goal == BObolLodPointCalibrationGoal::STATIC ?
		result.nextEvidence.staticPointProxyValue :
		result.nextEvidence.responsivePointProxyValue;
	result.pointProxyDecision = pointEvidence.interrupted(
	    currentThreshold, renderNanoseconds, targetFps);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::planPointStructuralCoverageBlocked(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, float currentThreshold,
	size_t unresolvedStructuralCount)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.pointProxyDecision =
	    result.nextEvidence.responsivePointProxyValue.
		structuralCoverageBlocked(
		currentThreshold, unresolvedStructuralCount);
	return result;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::planPointCompleted(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodPointCalibrationGoal goal, float currentThreshold,
	uint64_t renderNanoseconds, float targetFps, bool reusableSample,
	bool allowCoarsening,
	size_t unresolvedStructuralCount)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	BObolLodPointProxyEvidence &pointEvidence =
	    goal == BObolLodPointCalibrationGoal::STATIC ?
		result.nextEvidence.staticPointProxyValue :
		result.nextEvidence.responsivePointProxyValue;
	result.pointProxyDecision = pointEvidence.completed(
	    currentThreshold, renderNanoseconds, targetFps, reusableSample,
	    allowCoarsening,
	    unresolvedStructuralCount);
	return result;
}

bool
BObolLodAdmissionPlanner::adaptivePointAggregationAllowed(
	bool meshFirstSceneSafe, bool staticQualityCapacityRejected)
{
	return !meshFirstSceneSafe || staticQualityCapacityRejected;
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::settlePointAtStructuralLimit(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodPointCalibrationGoal goal, float currentThreshold)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	BObolLodPointProxyEvidence &pointEvidence =
	    goal == BObolLodPointCalibrationGoal::STATIC ?
		result.nextEvidence.staticPointProxyValue :
		result.nextEvidence.responsivePointProxyValue;
	pointEvidence.settleAtStructuralLimit(currentThreshold);
	result.pointProxyDecision.threshold = currentThreshold;
	return result;
}

bool
BObolLodAdmissionPlanner::pointRequiresReusableConfirmation(float currentThreshold,
    size_t unresolvedStructuralCount)
{
	return BObolLodPointProxyEvidence::requiresReusableConfirmation(
	    currentThreshold, unresolvedStructuralCount);
}

float
BObolLodAdmissionPlanner::pointStructuralPreloadThreshold(
    float currentThreshold, float candidateThreshold)
{
    const float current = BObolLodPointProxyEvidence::sanitize(
	currentThreshold);
    const float candidate = BObolLodPointProxyEvidence::sanitize(
	candidateThreshold);
    if (candidate >= current)
	return candidate;

    float bucketLimit = 1.0f;
    while (bucketLimit * 2.0f <= candidate)
	bucketLimit *= 2.0f;
    return bucketLimit;
}

bool
BObolLodAdmissionPlanner::pointRelaxationFitsDeadline(
    size_t activeOccurrenceCount, size_t additionalOccurrenceCount,
    uint64_t measuredRenderNanoseconds,
    uint64_t presentationDeadlineNanoseconds)
{
    if (!additionalOccurrenceCount)
	return true;
    if (!activeOccurrenceCount || !measuredRenderNanoseconds)
	return false;
    if (!presentationDeadlineNanoseconds)
	return true;

    const long double projectedOccurrenceCount =
	static_cast<long double>(activeOccurrenceCount) +
	static_cast<long double>(additionalOccurrenceCount);
    const long double projectedNanoseconds =
	static_cast<long double>(measuredRenderNanoseconds) *
	projectedOccurrenceCount /
	static_cast<long double>(activeOccurrenceCount);
    return std::isfinite(projectedNanoseconds) &&
	projectedNanoseconds <=
	    static_cast<long double>(presentationDeadlineNanoseconds);
}

bool
BObolLodAdmissionPlanner::pointAggregationApplicable(size_t visibleOccurrenceCount)
{
	return BObolLodPointProxyEvidence::applicable(
	    visibleOccurrenceCount);
}

bool
BObolLodAdmissionPlanner::pointAggregationApplicableAcrossCameraInvalidation(
	bool currentCensusApplicable, float retainedThreshold)
{
	return BObolLodPointProxyEvidence::
	    applicableAcrossCameraInvalidation(
		currentCensusApplicable, retainedThreshold);
}

bool
BObolLodAdmissionPlanner::pointStreamingPopulationWorkPending(bool submissionPending,
	bool resultsPending, bool servicePending, size_t providerPending)
{
	return BObolLodPointProxyEvidence::streamingPopulationWorkPending(
	    submissionPending, resultsPending, servicePending,
	    providerPending);
}

bool
BObolLodAdmissionPlanner::pointCalibrationYieldsToStructuralRepair(
	bool calibrationPending, bool exactOccurrenceClassification,
	size_t unresolvedStructuralCount, bool producerOwnsFutureFrame)
{
	return BObolLodPointProxyEvidence::
	    calibrationYieldsToStructuralRepair(calibrationPending,
		exactOccurrenceClassification, unresolvedStructuralCount,
		producerOwnsFutureFrame);
}

bool
BObolLodAdmissionPlanner::onePixelPresentationReady(bool haveCadAssemblies,
	bool exactFrame, bool exactOccurrenceClassification,
	size_t presentedSubpixelOccurrences,
	size_t presentedStructuralBoxes, int progressiveCeiling,
	int maximumActiveProgressiveCut,
	float pointProxyPixelThreshold, size_t managedRenderCost)
{
	return BObolLodPointProxyEvidence::onePixelPresentationReady(
	    haveCadAssemblies, exactFrame, exactOccurrenceClassification,
	    presentedSubpixelOccurrences, presentedStructuralBoxes,
	    progressiveCeiling, maximumActiveProgressiveCut,
	    pointProxyPixelThreshold, managedRenderCost);
}

bool
BObolLodAdmissionPlanner::capacitySamplePopulationReady(bool haveCadAssemblies,
	bool exactFrame, bool exactOccurrenceClassification,
	size_t presentedStructuralBoxes, bool allocationPresentationRealized)
{
	return BObolLodPointProxyEvidence::capacitySamplePopulationReady(
	    haveCadAssemblies, exactFrame, exactOccurrenceClassification,
	    presentedStructuralBoxes, allocationPresentationRealized);
}

bool
BObolLodAdmissionPlanner::capacitySampleYieldsToStructuralFrontier(
	bool capacityPresentationPending, bool exactFrame,
	bool exactOccurrenceClassification, size_t presentedStructuralBoxes)
{
	return BObolLodPointProxyEvidence::
	    capacitySampleYieldsToStructuralFrontier(
		capacityPresentationPending, exactFrame,
		exactOccurrenceClassification, presentedStructuralBoxes);
}

bool
BObolLodAdmissionPlanner::pointDeadlineRequiresPopulationAggregation(
	size_t activeRenderCost, size_t minimumRenderCost,
	int presentedMaximum, int correctedCeiling,
	size_t correctedRenderCostBudget)
{
	return BObolLodPointProxyEvidence::
	    deadlineRequiresPopulationAggregation(activeRenderCost,
		minimumRenderCost, presentedMaximum, correctedCeiling,
		correctedRenderCostBudget);
}

bool
BObolLodAdmissionPlanner::pointProducerOwnsCalibrationFrame(
	const PointCalibrationProducerInputs &inputs)
{
	const bool submissionPaused = presentationPausesSubmission(
	    inputs.discoveryCalibrationPending,
	    inputs.stableCalibrationPending,
	    inputs.capacityAllocationPending,
	    inputs.capacitySamplePending,
	    inputs.stablePresentationAvailable);
	return BObolLodPointProxyEvidence::producerOwnsCalibrationFrame(
	    inputs.submissionPending, submissionPaused, inputs.providerPending,
	    inputs.servicePending, inputs.publicationAwaitingFrameRequest);
}

bool
BObolLodAdmissionPlanner::pointCalibrationOwnsCompletedFrame(
	bool pointCalibrationPending, bool capacitySamplePending)
{
	return pointCalibrationPending && !capacitySamplePending;
}

bool
BObolLodAdmissionPlanner::pointBlocksSourceAdmission(
	bool discoveryCalibrationPending, bool stableCalibrationPending)
{
	return BObolLodPointProxyEvidence::blocksSourceAdmission(
	    discoveryCalibrationPending, stableCalibrationPending);
}

bool
BObolLodAdmissionPlanner::structuralFallbackAdmissionDeferred(
	bool quietView, bool coverageComplete, bool exactFrame,
	size_t presentedStructuralBoxes, size_t terminalOccurrenceFailures,
	bool directPreviewAuthorized, bool structuralRepairActive)
{
    return quietView && coverageComplete && exactFrame &&
	presentedStructuralBoxes > terminalOccurrenceFailures &&
	!directPreviewAuthorized && !structuralRepairActive;
}

bool
BObolLodAdmissionPlanner::structuralCoverageOnly(
	bool terminalRefinement, bool coverageActive,
	bool demandRefreshActive, bool directPreviewAuthorized,
	bool fallbackAdmissionDeferred)
{
    return !terminalRefinement &&
	((coverageActive && !demandRefreshActive &&
	  !directPreviewAuthorized) || fallbackAdmissionDeferred);
}

bool
BObolLodAdmissionPlanner::presentationPausesSubmission(
	bool discoveryCalibrationPending, bool stableCalibrationPending,
	bool capacityAllocationPending, bool capacitySamplePending,
	bool stablePresentationAvailable)
{
    /* ALLOCATING is a producer phase, not a presentation phase.  A point
     * calibration may already be queued as the capacity transaction's
     * successor, but it cannot pause the retained-allocation cursor which
     * creates the candidate it is waiting behind. */
    return capacitySamplePending || discoveryCalibrationPending ||
	    (stableCalibrationPending && stablePresentationAvailable &&
	     !capacityAllocationPending);
}

bool
BObolLodAdmissionPlanner::maySeedPointStructuralDistribution(
	bool discoveryCalibrationPending, bool stableCalibrationPending,
	bool triangleRecoveryPending)
{
	return BObolLodPointProxyEvidence::maySeedStructuralDistribution(
	    discoveryCalibrationPending, stableCalibrationPending,
	    triangleRecoveryPending);
}

bool
BObolLodAdmissionPlanner::pointAtMaximumPixelThreshold(float threshold)
{
    return BObolLodPointProxyEvidence::atMaximumPixelThreshold(threshold);
}

float
BObolLodAdmissionPlanner::pointCalibrationTargetFps(
	BObolLodPointCalibrationGoal goal, float responsiveTargetFps,
	float staticTargetFps)
{
    return goal == BObolLodPointCalibrationGoal::STATIC ?
	staticTargetFps : responsiveTargetFps;
}

uint64_t
BObolLodAdmissionPlanner::pointCalibrationPresentationDeadline(
	BObolLodPointCalibrationGoal goal, uint64_t responsiveDeadline,
	uint64_t staticDeadline)
{
    return goal == BObolLodPointCalibrationGoal::STATIC ?
	staticDeadline : responsiveDeadline;
}

bool
BObolLodAdmissionPlanner::pointTerminalReplayRequired(
	const BObolLodAdmissionEvidence &evidence,
	BObolLodPointCalibrationGoal goal, float threshold)
{
	return evidence.pointProxy(goal).terminalReplayRequired(threshold);
}

bool
BObolLodAdmissionPlanner::shouldRecoverTriangleDetail(
	bool reducibleProgressiveDetail, bool stableSampleOverloaded,
	bool coarsePointCut, bool protectedQualityOwnsCuts)
{
	return BObolLodPointProxyEvidence::shouldRecoverTriangleDetail(
	    reducibleProgressiveDetail, stableSampleOverloaded, coarsePointCut,
	    protectedQualityOwnsCuts);
}

float
BObolLodAdmissionPlanner::triangleRecoveryPointThreshold(float currentThreshold,
	bool exactStructuralPopulation, size_t structuralFallbackCount)
{
	return BObolLodPointProxyEvidence::triangleRecoveryPointThreshold(
	    currentThreshold, exactStructuralPopulation,
	    structuralFallbackCount);
}

BObolLodAdmissionPlan
BObolLodAdmissionPlanner::planStructural(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodStructuralAdmissionEvidence::Inputs &inputs)
{
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.structuralDecision =
	    result.nextEvidence.structuralValue.planStructural(inputs);
	return result;
}

bool
BObolLodAdmissionPlanner::structuralCapacityFrameApplicable(bool exactFrame,
	bool renderCostObserved, uint64_t renderNanoseconds,
	uint64_t presentationDeadlineNanoseconds)
{
	/* A pure structural frame has a valid observed mesh cost of zero.  The
	 * scene allocator's standing seed budget supplies its first-mesh
	 * allowance; requiring positive mesh cost here strands that exact box
	 * frontier without a producer. */
	return exactFrame && renderCostObserved && renderNanoseconds > 0 &&
	    presentationDeadlineNanoseconds > 0;
}

uint64_t
BObolLodAdmissionPlanner::structuralRepairPresentationDeadline(
	bool directMeshScene, uint64_t stableDeadlineNanoseconds,
	uint64_t prominentDeadlineNanoseconds, uint64_t censusMultiplier)
{
	const uint64_t censusDeadline = stableDeadlineNanoseconds > 0 &&
	    censusMultiplier > UINT64_MAX / stableDeadlineNanoseconds ?
		UINT64_MAX : stableDeadlineNanoseconds * censusMultiplier;
	return directMeshScene ?
	    std::max(censusDeadline, prominentDeadlineNanoseconds) :
	    censusDeadline;
}

size_t
BObolLodAdmissionPlanner::unaggregatableStructuralCount(
	const std::array<size_t, 7> &cumulativeCount, size_t visibleCount)
{
	const size_t aggregatable = std::min(
	    visibleCount, cumulativeCount.back());
	return visibleCount - aggregatable;
}



    /* Divide only proven marginal capacity.  A zero result means the exact
     * frame did not establish even one scheduling unit per unresolved
     * occurrence, so the caller must not arm a partial repair transaction. */
size_t
BObolLodAdmissionPlanner::structuralPerOccurrenceReservation(size_t admittedBudget,
	size_t activeCost, size_t unresolvedCount)
{
	if (!unresolvedCount || admittedBudget <= activeCost)
	    return 0;
	const size_t available = admittedBudget - activeCost;
	return available >= unresolvedCount ?
	    available / unresolvedCount : 0;
}



    /* noteDeadlineCapacityMiss() records an already strict, view-local upper
     * bound: it is the immediate predecessor of the candidate which missed
     * the hard deadline.  Marginal recovery has its own throughput safety
     * factor before reaching this helper.  Applying another factor to this
     * exact bound would discard untested capacity and can skip a useful
     * discrete PoP population. */
size_t
BObolLodAdmissionPlanner::capBudgetAtDeadlineCeiling(size_t estimatedBudget,
	size_t deadlineCapacityCeiling)
{
	return deadlineCapacityCeiling == SIZE_MAX ? estimatedBudget :
	    std::min(estimatedBudget, deadlineCapacityCeiling);
}

bool
BObolLodAdmissionPlanner::rejectAfterInterruptedFrame(bool interactive,
	bool cadDrawAttempted, bool preparationOnlyRetry,
	bool staticTrialActive)
{
	return !interactive && cadDrawAttempted && !preparationOnlyRetry &&
	    staticTrialActive;
}



    /* The one-pixel point trial and its enclosing static-quality trial own one
     * candidate framebuffer.  If that frame is rejected, both owners must
     * release it together.  Retiring an unrelated point calibration would
     * instead skip a required exact-classification frame. */
bool
BObolLodAdmissionPlanner::rejectedOnePixelTrialReleasesCalibration(bool trialRejected,
	bool calibrationPending, float pointProxyPixelThreshold)
{
	return trialRejected && calibrationPending &&
	    std::isfinite(pointProxyPixelThreshold) &&
	    pointProxyPixelThreshold <= 1.01f;
}



    /* A deadline miss can close the static-quality phase only when it was an
     * explicit static-quality trial against a closed occurrence population.
     * An ordinary quiet-frame miss is capacity evidence for that exact
     * attempted population, but it must not retire a static phase which has
     * not started yet.  Doing so made a cold scene remember a transient
     * realization miss and permanently forgo the eventual event-driven 10 Hz
     * quality pass. */
bool
BObolLodAdmissionPlanner::terminalCapacityMiss(bool interactive,
	bool cadDrawAttempted, bool transientPresentationRetry,
	bool staticTrialActive, bool populationWorkPending)
{
	return !interactive && cadDrawAttempted &&
	    !transientPresentationRetry && staticTrialActive &&
	    !populationWorkPending;
}

bool
BObolLodAdmissionPlanner::ordinaryHeadroomAllowed(bool staticPhaseBlocked)
{
	return !staticPhaseBlocked;
}



    /* A bounded static phase may improve an unaffordable atomic quality floor
     * incrementally.  The same marginal path remains available after either
     * the complete floor or a richer presentation has been rejected, but it
     * must never leak into ordinary quiet convergence. */
bool
BObolLodAdmissionPlanner::marginalStaticCapacityAllowed(bool staticPhaseActive,
	bool staticPhaseRejected, bool protectedFloorRejected)
{
	return staticPhaseActive || staticPhaseRejected ||
	    protectedFloorRejected;
}

bool
BObolLodAdmissionPlanner::actionableProgressiveQualityDebt(
	size_t activePayloadCount, size_t satisfiedPayloadCount,
	size_t memoryLimitedPayloadCount, int activeMaximumCut,
	int presentationCeiling)
{
	const size_t unsatisfiedPayloadCount =
	    activePayloadCount > satisfiedPayloadCount ?
		activePayloadCount - satisfiedPayloadCount : 0;
	const bool residentDemandDebt =
	    unsatisfiedPayloadCount > memoryLimitedPayloadCount;
	const bool presentationDebt = progressivePresentationQualityDebt(
	    activeMaximumCut, presentationCeiling);
	return residentDemandDebt || presentationDebt;
}

bool
BObolLodAdmissionPlanner::progressivePresentationQualityDebt(int activeMaximumCut,
	int presentationCeiling)
{
	return presentationCeiling >= 0 &&
	    activeMaximumCut > presentationCeiling;
}



    /* A renderer-wide PoP ceiling is a transient capacity guard for a
	* multi-occurrence scene.  Once no live static staircase owns it (or that
	* staircase has rejected its next step), an exact quiet frame must be
	* translated into occurrence-local cuts before readiness.  For exactly one
	* progressive occurrence the ordinal is already occurrence-local, so the
	* proven ceiling is a valid terminal representation and avoids a redundant
	* allocator/repaint transaction. */
    /* Translate a renderer-wide CAD cut into the retained allocator's scene
     * currency.  Backend submitted-work records may expand indexed vertices,
     * normals, and triangles differently (OSMesa flat shading is the important
     * case), so they are not a valid occurrence-allocation budget.  Preserve
     * the scene's non-CAD floor and replace only its retained CAD component with
     * the canonical cost table for the exact presented ceiling. */
size_t
BObolLodAdmissionPlanner::canonicalSceneCostAtCadCeiling(size_t activeSceneCost,
	size_t activeCadCost, size_t cadCeilingCost)
{
	const size_t nonCadCost = activeSceneCost > activeCadCost ?
	    activeSceneCost - activeCadCost : 0;
	return cadCeilingCost > SIZE_MAX - nonCadCost ? SIZE_MAX :
	    nonCadCost + cadCeilingCost;
}



    /* A retained allocation certificate is the exact post-transaction
     * population.  The admission cursor deliberately freezes the pre-pass
     * scene cost for bounded traversal accounting, so using that older value
     * as the capacity candidate makes the completed frame look stale even
     * when no cut changed. */
size_t
BObolLodAdmissionPlanner::capacitySearchPresentedCost(bool currentAllocation,
	size_t frozenActiveCost, size_t selectedPresentationCost)
{
	return currentAllocation && selectedPresentationCost > 0 ?
	    selectedPresentationCost : frozenActiveCost;
}

bool
BObolLodAdmissionPlanner::onePixelTrialRequired(float pointProxyPixelThreshold)
{
	return std::isfinite(pointProxyPixelThreshold) &&
	    pointProxyPixelThreshold > 1.01f;
}



    /* A coarse point cut may hide either resident meshes or structural
     * fallback boxes.  Static headroom may expose the former under its hard
     * deadline, but exposing the latter is a visual regression, not quality
     * refinement.  Retain aggregation whenever structural records remain,
     * independent of whether a renderer-wide PoP ceiling is also installed. */
bool
BObolLodAdmissionPlanner::retainAggregatedPresentation(int progressiveCeiling,
	float pointProxyPixelThreshold, size_t activePayloadCount,
	bool structuralFallbackPopulation)
{
	return std::isfinite(pointProxyPixelThreshold) &&
	    pointProxyPixelThreshold > 1.01f && activePayloadCount > 1 &&
	    (progressiveCeiling >= 0 || structuralFallbackPopulation);
}

bool
BObolLodAdmissionPlanner::startOnePixelTrialFromSettledPointFrame(bool interactive,
	bool exactReusableFrame, bool producerWork, bool resourcePressure,
	bool structuralFallbackPopulation, bool staticPhaseBlocked,
	float pointProxyPixelThreshold,
	uint64_t renderNanoseconds, uint64_t preferredFrameNanoseconds,
	uint64_t staticDeadlineNanoseconds)
{
	return !interactive && exactReusableFrame && !producerWork &&
	    !resourcePressure && !structuralFallbackPopulation &&
	    !staticPhaseBlocked &&
	    onePixelTrialRequired(pointProxyPixelThreshold) &&
	    renderNanoseconds > 0 && preferredFrameNanoseconds > 0 &&
	    staticDeadlineNanoseconds > preferredFrameNanoseconds &&
	    renderNanoseconds <= staticDeadlineNanoseconds;
}

bool
BObolLodAdmissionPlanner::acceptSettledOnePixelFrame(bool interactive,
	bool exactReusableFrame, bool producerWork, bool resourcePressure,
	bool structuralFallbackPopulation, bool handoffPending,
	bool staticPhaseBlocked,
	int progressiveCeiling, float pointProxyPixelThreshold,
	bool reducibleProgressiveDetail, uint64_t renderNanoseconds,
	uint64_t preferredFrameNanoseconds,
	uint64_t staticDeadlineNanoseconds)
{
	return !interactive && exactReusableFrame && !producerWork &&
	    !resourcePressure && !structuralFallbackPopulation &&
	    !handoffPending && !staticPhaseBlocked &&
	    progressiveCeiling < 0 &&
	    std::isfinite(pointProxyPixelThreshold) &&
	    pointProxyPixelThreshold <= 1.01f &&
	    reducibleProgressiveDetail && renderNanoseconds > 0 &&
	    preferredFrameNanoseconds > 0 &&
	    staticDeadlineNanoseconds > preferredFrameNanoseconds &&
	    renderNanoseconds > preferredFrameNanoseconds &&
	    renderNanoseconds <= staticDeadlineNanoseconds;
}

int
BObolLodAdmissionPlanner::stagedProgressiveCeiling(int currentCeiling,
	    int activeMaximum, float pointProxyPixelThreshold,
	    size_t activePayloadCount)
{
	if (currentCeiling < 0 || activeMaximum <= currentCeiling ||
	    activePayloadCount != 1 ||
	    !std::isfinite(pointProxyPixelThreshold) ||
	    pointProxyPixelThreshold <= 0.0f)
	    return -1;
	return std::min(activeMaximum, currentCeiling + 1);
}
int
BObolLodAdmissionPlanner::staticProgressiveBaselineCeiling(int currentCeiling,
	    int activeMaximum, size_t activePayloadCount)
{
	return currentCeiling < 0 && activeMaximum >= 0 &&
	    activePayloadCount == 1 ? activeMaximum : -1;
}

bool
BObolLodAdmissionPlanner::staticOrdinalOverscanApplicable(size_t activePayloadCount)
{
	return activePayloadCount == 1;
}

bool
BObolLodAdmissionPlanner::retainedPointAggregationApplicable(
	bool logicalPopulationApplicable, size_t pointProxyCandidateCount)
{
    return logicalPopulationApplicable && pointProxyCandidateCount > 0;
}

bool
BObolLodAdmissionPlanner::deadlinePointAggregationApplicable(
	bool retainedPopulationApplicable, bool structuralAdmissionPending,
	size_t logicalOccurrenceCount)
{
    return retainedPopulationApplicable ||
	(structuralAdmissionPending && logicalOccurrenceCount > 1);
}

bool
BObolLodAdmissionPlanner::pointProxyThresholdInert(
	size_t retainedCandidateCount, bool exactStructuralProjection,
	size_t structuralOccurrenceCount, size_t terminalFailureCount)
{
    return retainedCandidateCount == 0 && exactStructuralProjection &&
	structuralOccurrenceCount <= terminalFailureCount;
}

bool
BObolLodAdmissionPlanner::structuralPointAggregationRequired(
	size_t logicalOccurrenceCount, size_t affordableOccurrenceCount)
{
	return affordableOccurrenceCount != SIZE_MAX &&
	    logicalOccurrenceCount > affordableOccurrenceCount;
}

bool
BObolLodAdmissionPlanner::pointPolicyExhaustedForStructuralFrontier(
	bool maximumThreshold, bool pointAggregationRequired,
	bool exactProjection, size_t unresolvedCount,
	size_t unaggregatableCount)
{
	if (!unresolvedCount)
	    return false;
	return maximumThreshold || !pointAggregationRequired ||
	    (exactProjection && unaggregatableCount >= unresolvedCount);
}

bool
BObolLodAdmissionPlanner::pointSuccessorRequiredForStructuralFrontier(
	bool exactProjection, size_t unresolvedCount,
	bool pointPolicyExhausted, bool anotherOwnerPending)
{
	/* A successful terminal presentation may not retain structural boxes.
	 * If an exact quiet frame still has an aggregatable box frontier, the
	 * point threshold owns the only successor which can reduce that frontier
	 * without reopening a whole-scene provider pass. */
	return exactProjection && unresolvedCount > 0 &&
	    !pointPolicyExhausted && !anotherOwnerPending;
}

size_t
BObolLodAdmissionPlanner::structuralFirstWaveOccurrenceLimit(size_t sceneBudget)
{
	static constexpr size_t occurrenceCostFloor = 64;
	static constexpr size_t minimumOccurrences = 512;
	static constexpr size_t maximumOccurrences = 8192;
	if (sceneBudget == SIZE_MAX)
	    return SIZE_MAX;
	return std::max(minimumOccurrences,
	    std::min(maximumOccurrences, sceneBudget / occurrenceCostFloor));
}
BObolLodAdmissionPlan
BObolLodAdmissionPlanner::beginPlan(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor)
{
	BObolLodAdmissionPlan result;
	result.nextEvidence = evidence;
	result.nextCursor = cursor;
	return result;
}
