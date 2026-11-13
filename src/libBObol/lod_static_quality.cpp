/*          L O D _ S T A T I C _ Q U A L I T Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"
#include "lod_view_policy_private.h"

BObolLodStaticQualityTrial::CompletedFrameDecision
BObolLodStaticQualityTrial::completeFrame(const CompletedFrameInputs &frame)
{
    using Action = CompletedFrameAction;
    if (!this->probing() || frame.interactive || frame.ceiling < 0 ||
	!frame.presentationReady)
	return {};

    /* Upload/command preparation can finish an exact frame without proving
     * its sustainable cost.  One unchanged traversal supplies the reusable
     * sample; a preparation-heavy duration cannot reject a retained cut. */
    if (frame.transientPresentation)
	return {Action::REPLAY_PREPARATION};

    if (!frame.exactCadFrame || !frame.measuredCadWork ||
	frame.renderNanoseconds > frame.deadlineNanoseconds ||
	frame.handoffPending)
	return {};

    if (frame.nextFraction > 0.0f) {
	/* Occurrence-wide reconciliation cannot encode this mixed page
	 * population.  Retain its exact renderer policy with its acceptance. */
	Acceptance acceptance;
	acceptance.revisionStamp = frame.revisionStamp;
	acceptance.ceiling = frame.ceiling;
	acceptance.nextFraction = frame.nextFraction;
	acceptance.presentedCost = frame.canonicalSceneCost;
	acceptance.allowedCost = std::max(acceptance.presentedCost,
	    BObolLodQualityPolicy::staticPresentationRenderCostLimit(
		frame.canonicalSceneCost, frame.renderNanoseconds,
		frame.deadlineNanoseconds));
	if (!this->acceptFractionalCeiling(acceptance))
	    return {};
	return {Action::ACCEPT_FRACTION, frame.ceiling, frame.nextFraction};
    }

    if (this->sampledCeiling() != frame.ceiling) {
	this->noteSampledCeiling(frame.ceiling);
	return {Action::REPLAY_STEADY};
    }

    /* Once all occurrence-local cuts fit under the guard, removing it is
     * exact.  Keep the static trial/deadline: retiring it here would restore
     * the older coarse cadence budget and reopen the same quality search. */
    if (frame.activeMaximum <= frame.ceiling)
	return {Action::RELEASE_REDUNDANT_CEILING};

    const CutPrediction &prediction = frame.prediction;
    if (prediction.ceiling > frame.ceiling) {
	/* Prediction selects the next candidate; renderer interruption remains
	 * the independent bound if that estimate is optimistic. */
	this->resetSample();
	return {Action::PRESENT_CUT, prediction.ceiling};
    }
    if (prediction.ceiling >= 0 && prediction.nextFraction > 0.0f)
	return {Action::PRESENT_FRACTION, frame.ceiling, prediction.nextFraction};

    if (!prediction.presentedCost || !prediction.allowedCost)
	return {};
    CompletedFrameDecision handoff;
    handoff.action = Action::HANDOFF_POPULATION;
    handoff.presentedCost = prediction.presentedCost;
    handoff.allowedCost = prediction.allowedCost;
    /* The single-occurrence query deliberately has no ordinal answer for a
     * multi-occurrence population.  Its -1 sentinel is not proof that the
     * renderer guard is redundant.  Transfer the measured budget to the
     * existing occurrence allocator, preserving the guard until handoff.
     * Absence of a prediction is not a capacity rejection either. */
    if (prediction.ceiling < 0)
	return handoff;

    /* The measured current cut is the safe side of this rejection.  The
     * occurrence handoff represents that proof before releasing the guard;
     * it must not restart a speculative rich frame. */
    Constraint constraint;
    constraint.reason = ConstraintReason::PREDICTED_NEXT_CUT;
    constraint.revisionStamp = frame.revisionStamp;
    constraint.committedCeiling = frame.ceiling;
    constraint.committedCost = prediction.presentedCost;
    constraint.candidateCeiling = frame.ceiling + 1;
    constraint.candidateCost = prediction.nextCutCost;
    constraint.allowedCost = prediction.allowedCost;
    if (!this->reject(constraint))
	return {};
    handoff.action = Action::RECONCILE_CUT;
    return handoff;
}
