/*          L O D _ V I E W _ P O L I C Y _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_VIEW_POLICY_PRIVATE_H
#define LIBBOBOL_LOD_VIEW_POLICY_PRIVATE_H

#include "common.h"

#include "BObol/BMeshLodCache.h"
#include "lod_revision_private.h"
#include "lod_view_snapshot_private.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

class BObolLodViewDemandPolicy {
public:
    static constexpr int64_t unbracketedQuietDebounceMicroseconds(void)
    {
	return 150000;
    }

    static constexpr uint64_t qualityFrameDurationNanoseconds(void)
    {
	return 100000000ULL;
    }

    /* A remembered quiet frame is a retained-image seed, not permission to
     * reuse the same expensive traversal during input.  Allow a small,
     * bounded scheduling tolerance when recording it: a software renderer
     * can cross the nominal static deadline by a few milliseconds because of
     * compositor or readback work even though the completed image remains a
     * useful immediate restore.  A later input epoch still uses the strict
     * interactive deadline and ordinary admission must revalidate the seed. */
    static constexpr uint64_t staticHistoryDeadlineNanoseconds(
	uint64_t deadline)
    {
	return !deadline ? 0 :
	    deadline > UINT64_MAX / 6 ? UINT64_MAX :
	    deadline + deadline / 5;
    }

    static constexpr float qualityTargetFramesPerSecond(void)
    {
	return 10.0f;
    }

    struct CameraChangeDecision {
	bool retainPriorQualityCeiling = false;
    };

    struct QualityProbeInputs {
	int activeMaximum = -1;
	int presentationCeiling = -1;
	uint64_t presentedWork = 0;
	bool refinementFramePending = false;
	bool publicationFramePending = false;
	bool motionFramePending = false;
    };

    struct QualityProbeDecision {
	bool consumed = false;
	bool begin = false;
	int progressiveCeiling = -1;
    };

    struct DemandPassInputs {
	bool submissionActive = false;
	bool rescanPending = false;
	bool selectivePass = false;
	bool structuralRepair = false;
	bool retainedAllocation = false;
	/* A stronger finite transaction may temporarily own the shared source
	 * cursor.  The demand obligation remains level-triggered until that owner
	 * retires, but it must not create a competing ordinary pass meanwhile. */
	bool strongerOwnerPending = false;
    };

    CameraChangeDecision observeCameraChange(bool scaleChanged,
	uint64_t lastRenderNanoseconds)
    {
	CameraChangeDecision decision;
	decision.retainPriorQualityCeiling =
	    scaleChanged && this->qualityBudgetActiveValue &&
	    this->qualityProbePresentedValue &&
	    lastRenderNanoseconds > 0 &&
	    lastRenderNanoseconds <= qualityFrameDurationNanoseconds();
	this->qualityProbeActiveValue = false;
	this->qualityProbePendingValue = false;
	this->qualityProbePresentedValue = false;
	this->viewScaleChangingValue = scaleChanged;
	return decision;
    }

    void seedQualityFloor(int progressiveCeiling)
    {
	this->qualityFloorValue = progressiveCeiling >= 0 ?
	    std::min(progressiveCeiling,
		BOBOL_MESH_LOD_CUT_COUNT_MAX - 1) : -1;
    }

    void beginGesture(bool newInteractionEpoch)
    {
	if (newInteractionEpoch) {
	    this->interactionScaleChangedValue = false;
	    this->qualityFloorValue = -1;
	    this->qualityCeilingLimitValue = -1;
	    this->qualityCeilingFailedWorkValue = 0;
	}
	this->viewScaleChangingValue = false;
	this->qualityProbeActiveValue = false;
	this->qualityProbePendingValue = false;
	this->qualityProbePresentedValue = false;
	this->qualityBudgetActiveValue = false;
    }

    void beginCameraInteraction(bool newInteractionEpoch,
	bool scaleChanged)
    {
	if (newInteractionEpoch) {
	    this->interactionScaleChangedValue = false;
	    this->qualityFloorValue = -1;
	    this->qualityCeilingLimitValue = -1;
	    this->qualityCeilingFailedWorkValue = 0;
	}
	if (scaleChanged)
	    this->interactionScaleChangedValue = true;
    }

    void endGesture(void)
    {
	this->qualityProbeActiveValue = false;
	this->qualityProbePendingValue = false;
	this->qualityProbePresentedValue = false;
    }

    void noteFramePresented(int progressiveCeiling, bool exact,
	uint64_t elapsedNanoseconds)
    {
	if (!this->qualityBudgetActiveValue)
	    return;
	const bool responsive = exact && elapsedNanoseconds > 0 &&
	    elapsedNanoseconds <= qualityFrameDurationNanoseconds();
	this->qualityProbePresentedValue = responsive;
	if (responsive && progressiveCeiling >= 0)
	    this->qualityFloorValue = std::max(this->qualityFloorValue,
		std::min(progressiveCeiling,
		    BOBOL_MESH_LOD_CUT_COUNT_MAX - 1));
    }

    void noteQualityMiss(int correctedCeiling, uint64_t failedWork = 0)
    {
	if (correctedCeiling < 0)
	    return;
	if (this->qualityFloorValue < 0)
	    this->qualityFloorValue = correctedCeiling;
	else
	    this->qualityFloorValue = std::min(
		this->qualityFloorValue, correctedCeiling);
	/* A hard deadline is direct evidence that the next richer population is
	 * not sustainable in this interaction epoch.  Retain the corrected
	 * renderer ceiling across subsequent camera samples.  Otherwise every
	 * wheel event rearms the same +1 quality probe, producing a visible
	 * rich/coarse cycle and repeatedly spending one long software frame on a
	 * result already known to miss its contract.  A still coarser correction
	 * tightens this bound; a new gesture/quiet epoch clears it. */
	if (this->qualityCeilingLimitValue < 0)
	    this->qualityCeilingLimitValue = correctedCeiling;
	else
	    this->qualityCeilingLimitValue = std::min(
		this->qualityCeilingLimitValue, correctedCeiling);
	if (failedWork)
	    this->qualityCeilingFailedWorkValue = failedWork;
	this->qualityProbePresentedValue = false;
    }

    bool rearmAfterQualityFrame(bool interactive, int activeMaximum,
	int presentedMaximum, bool exact, uint64_t elapsedNanoseconds)
    {
	/* A cache hit can make a richer occurrence cut active without publishing
	 * a worker result, so rearmAfterResidentGrowth() alone cannot guarantee
	 * that already resident detail is ever offered.  Permit one successor
	 * only after the current quality cut was exactly presented inside the hard
	 * deadline.  The completed-frame barrier serializes the steps, and a
	 * missed deadline takes the correction path instead of rearming. */
	if (!interactive || !this->viewScaleChangingValue ||
	    !this->interactionScaleChangedValue ||
	    !this->qualityBudgetActiveValue ||
	    !this->qualityProbeActiveValue || !exact ||
	    !elapsedNanoseconds ||
	    elapsedNanoseconds > qualityFrameDurationNanoseconds() ||
	    activeMaximum < 0 || presentedMaximum < 0 ||
	    activeMaximum <= presentedMaximum)
	    return false;
	this->qualityProbeActiveValue = false;
	this->qualityProbePendingValue = true;
	return true;
    }

    bool noteMotionFrameSettled(void)
    {
	if (!this->viewScaleChangingValue ||
	    !this->interactionScaleChangedValue ||
	    this->qualityProbeActiveValue)
	    return false;
	this->qualityProbePendingValue = true;
	return true;
    }

    QualityProbeDecision beginQualityProbe(
	const QualityProbeInputs &inputs)
    {
	QualityProbeDecision decision;
	if (!this->qualityProbePendingValue ||
	    !this->interactionScaleChangedValue ||
	    this->qualityProbeActiveValue ||
	    inputs.refinementFramePending ||
	    inputs.publicationFramePending ||
	    inputs.motionFramePending)
	    return decision;

	/* The pending edge belongs to this attempt even when no progressive
	 * population exists yet.  A later resident-growth result explicitly
	 * rearms it, avoiding an idle-loop retry against activeMaximum == -1. */
	this->qualityProbePendingValue = false;
	decision.consumed = true;
	if (inputs.activeMaximum < 0)
	    return decision;

	this->qualityProbeActiveValue = true;
	this->qualityProbePresentedValue = false;
	this->qualityBudgetActiveValue = true;
	/* Residency may advance several populations while a render-only ceiling
	 * protects interaction.  Probe from the cut actually being presented,
	 * not from that hidden resident/occurrence maximum; jumping back to the
	 * latter on every wheel event repeatedly recreates the same long frame
	 * and ratchets a previously safe ceiling downward. */
	const int presentedMaximum = inputs.presentationCeiling >= 0 ?
	    std::min(inputs.activeMaximum, inputs.presentationCeiling) :
	    inputs.activeMaximum;
	decision.progressiveCeiling =
	    (presentedMaximum >= BOBOL_MESH_LOD_CUT_COUNT_MAX - 1 ?
		BOBOL_MESH_LOD_CUT_COUNT_MAX - 1 : presentedMaximum + 1);
	/* A stable cut which was presented within the hard zoom-quality
	 * deadline is a stronger starting point than the deliberately coarser
	 * coarser motion cut.  Replaying it does not expose new immutable data and
	 * is still protected by the presentation deadline.  This avoids walking
	 * back through several PoP populations every time wheel input pauses. */
	if (this->qualityFloorValue >= 0)
	    decision.progressiveCeiling = std::max(
		decision.progressiveCeiling,
		std::min(inputs.activeMaximum, this->qualityFloorValue));
	/* Consume a rearmed edge at the proven motion ceiling without starting a
	 * presentation attempt.  Resident suffix growth changes availability, not
	 * the measured cost of the already-failed next prefix, and therefore is not
	 * by itself capacity evidence. */
	if (this->qualityCeilingLimitValue >= 0 &&
	    decision.progressiveCeiling > this->qualityCeilingLimitValue) {
	    /* The limit records the cost of a failed presentation, not a
	     * view-independent property of that global cut.  Spatially clustered
	     * leaves may expose only a small fraction of the same prefix after a
	     * zoom or pan.  Permit a new bounded probe once exact visible work has
	     * fallen materially below the failed population; the hard presentation
	     * deadline still protects the endpoint if the estimate is optimistic. */
	    const bool materiallyCheaperView =
		this->qualityCeilingFailedWorkValue > 0 &&
		inputs.presentedWork > 0 &&
		inputs.presentedWork <=
		    this->qualityCeilingFailedWorkValue / 2;
	    if (!materiallyCheaperView) {
		this->qualityProbeActiveValue = false;
		return decision;
	    }
	    /* This is a new, substantially cheaper presentation population.  The
	     * old ceiling is no longer evidence about it.  Retire that witness so
	     * a successful probe may continue through successively richer cuts;
	     * any new miss immediately installs a view-appropriate bound. */
	    this->qualityCeilingLimitValue = -1;
	    this->qualityCeilingFailedWorkValue = 0;
	}
	decision.begin = true;
	return decision;
    }

    bool rearmAfterResidentGrowth(bool interactive)
    {
	/* A suffix which arrives after an earlier wheel event remains useful,
	 * but it must not turn a later pose-only segment of the same debounce
	 * epoch back into a scale-quality probe.  The epoch-wide bit is the quiet
	 * convergence obligation; the camera-local bit says whether the current
	 * presentation is actually changing scale. */
	if (!interactive || !this->viewScaleChangingValue ||
	    !this->interactionScaleChangedValue)
	    return false;
	this->qualityProbeActiveValue = false;
	this->qualityProbePendingValue = true;
	this->qualityProbePresentedValue = false;
	this->qualityBudgetActiveValue = true;
	return true;
    }

    bool finishQuietInteraction(bool returnedToStartingScale = false)
    {
	const bool completedScaleInteraction =
	    this->interactionScaleChangedValue &&
	    !returnedToStartingScale;
	this->viewScaleChangingValue = false;
	this->interactionScaleChangedValue = false;
	this->qualityProbeActiveValue = false;
	this->qualityProbePendingValue = false;
	this->qualityProbePresentedValue = false;
	this->qualityBudgetActiveValue = false;
	this->qualityFloorValue = -1;
	this->qualityCeilingLimitValue = -1;
	this->qualityCeilingFailedWorkValue = 0;
	return completedScaleInteraction;
    }

    void refreshForViewRevision(bool interactive)
    {
	this->demandRefreshActiveValue =
	    this->viewScaleChangingValue ||
	    (interactive && this->interactionScaleChangedValue);
    }

    /* A policy epoch which changes physical mesh demand needs one complete
     * occurrence retarget pass.  This is broader than camera-scale demand:
     * quiet static-quality tiers and explicit forced-cut changes also alter
     * the requested prefix without invalidating the visibility census. */
    void refreshForPolicyRevision(bool demandRefreshRequired,
	bool interactive)
    {
	this->demandRefreshActiveValue =
	    demandRefreshRequired ||
	    (interactive && this->interactionScaleChangedValue);
    }

    void clearDemandRefresh(void)
    {
	this->demandRefreshActiveValue = false;
    }

    void requestDemandRefresh(void)
    {
	this->demandRefreshActiveValue = true;
    }

    bool completeDemandRefresh(bool completedPass, bool selectivePass,
	bool rescanPending, bool structuralRepair, bool retainedAllocation)
    {
	if (!this->demandRefreshActiveValue || !completedPass ||
	    selectivePass || rescanPending || structuralRepair ||
	    retainedAllocation)
	    return false;
	this->clearDemandRefresh();
	return true;
    }

    /* Demand refresh is level-triggered.  A policy transition normally
     * starts its dense pass through the revision-mismatch path, but a
     * stronger presentation owner may retire that cursor without consuming
     * the demand obligation.  Once no specialized pass owns the cursor, the
     * ordinary current-view census must be re-armed even if every revision
     * and source signature is unchanged. */
    bool demandPassRequired(const DemandPassInputs &inputs) const
    {
	return this->demandRefreshActiveValue && !inputs.submissionActive &&
	    !inputs.rescanPending && !inputs.selectivePass &&
	    !inputs.structuralRepair && !inputs.retainedAllocation &&
	    !inputs.strongerOwnerPending;
    }

    void reset(void)
    {
	this->demandRefreshActiveValue = false;
	this->viewScaleChangingValue = false;
	this->interactionScaleChangedValue = false;
	this->qualityProbeActiveValue = false;
	this->qualityProbePendingValue = false;
	this->qualityProbePresentedValue = false;
	this->qualityBudgetActiveValue = false;
	this->qualityFloorValue = -1;
	this->qualityCeilingLimitValue = -1;
	this->qualityCeilingFailedWorkValue = 0;
    }

    bool demandRefreshActive(void) const
    {
	return this->demandRefreshActiveValue;
    }
    bool viewScaleChanging(void) const
    {
	return this->viewScaleChangingValue;
    }
    bool interactionScaleChanged(void) const
    {
	return this->interactionScaleChangedValue;
    }
    bool qualityProbeActive(void) const
    {
	return this->qualityProbeActiveValue;
    }
    bool qualityProbePending(void) const
    {
	return this->qualityProbePendingValue;
    }
    bool qualityProbePresented(void) const
    {
	return this->qualityProbePresentedValue;
    }
    bool qualityBudgetActive(void) const
    {
	return this->qualityBudgetActiveValue;
    }
    int qualityFloor(void) const
    {
	return this->qualityFloorValue;
    }
    int qualityCeilingLimit(void) const
    {
	return this->qualityCeilingLimitValue;
    }
    uint64_t qualityCeilingFailedWork(void) const
    {
	return this->qualityCeilingFailedWorkValue;
    }
    bool scaleChangingInteraction(bool interactive) const
    {
	/* interactionScaleChangedValue intentionally survives a mixed wheel and
	 * mouse gesture until quiet convergence.  It cannot classify the current
	 * frame: after a wheel event, rotation/translation must keep the existing
	 * occurrence cuts and use only the reversible renderer ceiling for FPS
	 * control. */
	return interactive && this->viewScaleChangingValue;
    }

private:
    bool demandRefreshActiveValue = false;
    bool viewScaleChangingValue = false;
    bool interactionScaleChangedValue = false;
    bool qualityProbeActiveValue = false;
    bool qualityProbePendingValue = false;
    bool qualityProbePresentedValue = false;
    bool qualityBudgetActiveValue = false;
    int qualityFloorValue = -1;
    int qualityCeilingLimitValue = -1;
    uint64_t qualityCeilingFailedWorkValue = 0;
};

/*
 * One owner for the camera-interaction lifecycle and its motion-frame gate.
 * The former interactive/gesture flags admitted an impossible gesture-only
 * combination, while the separately writable render serial and timestamp
 * could outlive the interaction they described.  This value names the three
 * legal phases and owns every transition between them.
 */
class BObolLodInteractionSession {
public:
    enum class Phase : uint8_t {
	QUIET = 0,
	DEBOUNCING,
	GESTURE
    };

    bool active(void) const
    {
	return this->phaseValue != Phase::QUIET;
    }

    bool gestureActive(void) const
    {
	return this->phaseValue == Phase::GESTURE;
    }

    Phase phase(void) const
    {
	return this->phaseValue;
    }

    void beginGesture(int64_t nowMicroseconds)
    {
	this->phaseValue = Phase::GESTURE;
	this->lastViewChangeMicrosecondsValue = nowMicroseconds;
	this->settleAfterRenderSerialValue = 0;
	this->releaseCutFloorActiveValue = false;
    }

    bool endGesture(int64_t nowMicroseconds)
    {
	if (!this->gestureActive())
	    return false;
	this->phaseValue = Phase::DEBOUNCING;
	this->lastViewChangeMicrosecondsValue = nowMicroseconds;
	this->settleAfterRenderSerialValue = 0;
	this->releaseCutFloorActiveValue = true;
	return true;
    }

    void clearMotionFrameGate(void)
    {
	this->settleAfterRenderSerialValue = 0;
    }

    void observeCameraChange(int64_t nowMicroseconds,
	uint64_t renderCompletionSerial)
    {
	if (!this->active())
	    this->phaseValue = Phase::DEBOUNCING;
	this->lastViewChangeMicrosecondsValue = nowMicroseconds;
	this->releaseCutFloorActiveValue = false;
	this->settleAfterRenderSerialValue = renderCompletionSerial + 1;
	if (this->settleAfterRenderSerialValue == 0)
	    this->settleAfterRenderSerialValue = 1;
    }

    bool noteMotionFrameCompleted(uint64_t renderCompletionSerial,
	int64_t nowMicroseconds)
    {
	if (!this->active() || !this->settleAfterRenderSerialValue ||
	    renderCompletionSerial < this->settleAfterRenderSerialValue)
	    return false;
	this->lastViewChangeMicrosecondsValue = nowMicroseconds;
	this->settleAfterRenderSerialValue = 0;
	return true;
    }

    bool releaseExpiredMotionFrame(int64_t nowMicroseconds,
	int64_t quietDebounceMicroseconds)
    {
	if (!this->active() || this->gestureActive() ||
	    !this->settleAfterRenderSerialValue ||
	    this->lastViewChangeMicrosecondsValue <= 0 ||
	    elapsedMicroseconds(nowMicroseconds,
		this->lastViewChangeMicrosecondsValue) <
		quietDebounceMicroseconds)
	    return false;
	this->settleAfterRenderSerialValue = 0;
	return true;
    }

    bool quietReady(int64_t nowMicroseconds,
	int64_t quietDebounceMicroseconds) const
    {
	return this->active() && !this->gestureActive() &&
	    !this->settleAfterRenderSerialValue &&
	    this->lastViewChangeMicrosecondsValue > 0 &&
	    elapsedMicroseconds(nowMicroseconds,
		this->lastViewChangeMicrosecondsValue) >=
		quietDebounceMicroseconds;
    }

    void finishQuiet(void)
    {
	this->phaseValue = Phase::QUIET;
	this->settleAfterRenderSerialValue = 0;
	this->lastViewChangeMicrosecondsValue = 0;
	this->releaseCutFloorActiveValue = false;
    }

    void reset(void)
    {
	this->phaseValue = Phase::QUIET;
	this->settleAfterRenderSerialValue = 0;
	this->lastViewChangeMicrosecondsValue = 0;
	this->releaseCutFloorActiveValue = false;
	this->ceilingFeedbackRenderSerialValue = 0;
    }

    bool hasNewCeilingFeedback(uint64_t renderCompletionSerial) const
    {
	return renderCompletionSerial !=
	    this->ceilingFeedbackRenderSerialValue;
    }

    void noteCeilingFeedback(uint64_t renderCompletionSerial)
    {
	this->ceilingFeedbackRenderSerialValue = renderCompletionSerial;
    }

    void resetCeilingFeedback(void)
    {
	this->ceilingFeedbackRenderSerialValue = 0;
    }

    uint64_t settleAfterRenderSerial(void) const
    {
	return this->settleAfterRenderSerialValue;
    }

    bool releaseCutFloorActive(void) const
    {
	return this->releaseCutFloorActiveValue;
    }

private:
    static int64_t elapsedMicroseconds(int64_t nowMicroseconds,
	int64_t thenMicroseconds)
    {
	return nowMicroseconds >= thenMicroseconds ?
	    nowMicroseconds - thenMicroseconds : 0;
    }

    Phase phaseValue = Phase::QUIET;
    uint64_t settleAfterRenderSerialValue = 0;
    int64_t lastViewChangeMicrosecondsValue = 0;
    bool releaseCutFloorActiveValue = false;
    uint64_t ceilingFeedbackRenderSerialValue = 0;
};

static_assert(std::is_trivially_copyable<
    BObolLodInteractionSession>::value,
    "interaction sessions must remain allocation-free values");

/* Evidence captured once at the start of a continuous camera-interaction
 * epoch.  Scale identity, prior readiness, and the last stable capacity are
 * one proof: independently clearing their former validity flags admitted
 * impossible half-captured states during service and renderer resets. */
class BObolLodInteractionStartCertificate {
public:
    void capture(const BObolLodViewScaleSnapshot &scale,
	bool exactScale, bool startedFromReadyView, size_t stableBudget)
    {
	this->stateValue = exactScale ? State::EXACT_SCALE : State::NO_SCALE;
	this->scaleValue = scale;
	this->startedFromReadyViewValue = startedFromReadyView;
	this->stableBudgetValue = stableBudget;
    }

    void reset(void)
    {
	this->stateValue = State::EMPTY;
	this->scaleValue = BObolLodViewScaleSnapshot();
	this->startedFromReadyViewValue = false;
	this->stableBudgetValue = 0;
    }

    bool captured(void) const
    {
	return this->stateValue != State::EMPTY;
    }

    bool exactScale(void) const
    {
	return this->stateValue == State::EXACT_SCALE;
    }

    bool returnedToScale(const BObolLodViewScaleSnapshot &scale) const
    {
	return this->exactScale() && this->scaleValue.same(scale);
    }

    bool startedFromReadyView(void) const
    {
	return this->captured() && this->startedFromReadyViewValue;
    }

    size_t stableBudget(void) const
    {
	return this->captured() ? this->stableBudgetValue : 0;
    }

private:
    enum class State : uint8_t {
	EMPTY = 0,
	NO_SCALE,
	EXACT_SCALE
    };

    State stateValue = State::EMPTY;
    BObolLodViewScaleSnapshot scaleValue;
    bool startedFromReadyViewValue = false;
    size_t stableBudgetValue = 0;
};

static_assert(std::is_trivially_copyable<
    BObolLodInteractionStartCertificate>::value,
    "interaction start certificates must remain allocation-free values");

/*
 * Lifecycle of one event-driven static-quality trial.  The former active and
 * rejected booleans admitted four combinations with implicit meanings spread
 * across the controller.  Name those states here and keep the replay sample
 * with its sole owner so callers cannot manufacture a fifth combination.
 */
class BObolLodStaticQualityTrial {
public:
    enum class ConstraintReason : uint8_t {
	NONE = 0,
	PRESENTATION_DEADLINE,
	PREDICTED_NEXT_CUT,
	PROTECTED_MINIMUM,
	ALLOCATION_SATURATED
    };

    struct Constraint {
	ConstraintReason reason = ConstraintReason::NONE;
	BObolLodAdmissionRevisionStamp revisionStamp;
	int committedCeiling = -1;
	float committedNextFraction = 0.0f;
	int candidateCeiling = -1;
	float candidateNextFraction = 0.0f;
	size_t committedCost = 0;
	size_t candidateCost = 0;
	size_t allowedCost = 0;

	bool valid(void) const
	{
	    return this->reason != ConstraintReason::NONE &&
		!this->revisionStamp.empty() &&
		validCut(this->committedCeiling,
		    this->committedNextFraction) &&
		validCut(this->candidateCeiling,
		    this->candidateNextFraction) &&
		this->candidateCost > 0;
	}

    private:
	static bool validCut(int ceiling, float nextFraction)
	{
	    return ceiling >= -1 &&
		ceiling < BOBOL_MESH_LOD_CUT_COUNT_MAX &&
		std::isfinite(nextFraction) && nextFraction >= 0.0f &&
		nextFraction <= 1.0f &&
		(ceiling >= 0 || nextFraction <= 1.0e-6f);
	}
    };

    struct Acceptance {
	BObolLodAdmissionRevisionStamp revisionStamp;
	int ceiling = -1;
	float nextFraction = 0.0f;
	size_t presentedCost = 0;
	size_t allowedCost = 0;

	bool valid(void) const
	{
	    return !this->revisionStamp.empty() && this->ceiling >= 0 &&
		this->ceiling < BOBOL_MESH_LOD_CUT_COUNT_MAX &&
		std::isfinite(this->nextFraction) &&
		this->nextFraction > 0.0f && this->nextFraction <= 1.0f &&
		this->presentedCost > 0 && this->allowedCost > 0 &&
		this->presentedCost <= this->allowedCost;
	}
    };

    struct HandoffCompletion {
	bool completedRejectedConstraint = false;
	bool releaseRendererCeiling = false;
	bool measureCeilingFreeCandidate = false;
    };

    void reset(void)
    {
	this->stateValue = State::IDLE;
	this->sampledCeilingValue = -1;
	this->baselinePointProxyPixelThresholdValue = 1.0f;
	this->constraintValue = Constraint();
	this->acceptanceValue = Acceptance();
    }

    void begin(float pointProxyPixelThreshold = 1.0f)
    {
	if (this->blocksNewTrial())
	    return;
	this->stateValue = State::PROBING;
	this->sampledCeilingValue = -1;
	this->baselinePointProxyPixelThresholdValue =
	    sanitizePointProxyPixelThreshold(pointProxyPixelThreshold);
    }

    void restoreRendererCeiling(bool active)
    {
	this->reset();
	this->stateValue = active ? State::PROBING : State::IDLE;
    }

    void deactivate(void)
    {
	this->reset();
    }

    /* A rejected static candidate is not terminal until its renderer-wide
     * safety ceiling has been represented by a certified occurrence-local
     * allocation.  Keep that successful handoff as the sole transition from
     * RECONCILING to REJECTED.  Generic policy/view invalidation uses
     * deactivate(); it must not be able to masquerade as this certificate. */
    bool completeReconciliation(void)
    {
	if (this->stateValue != State::RECONCILING)
	    return false;
	this->stateValue = State::REJECTED;
	this->sampledCeilingValue = -1;
	return true;
    }

    /* A completed occurrence-local handoff is the sole commit edge for the
     * renderer guard which protected its allocation.  A probing trial keeps
     * its quality owner but must expose the newly allocated population in one
     * ceiling-free, deadline-bounded frame.  A rejected trial completes its
     * constraint reconciliation.  Only the narrow protected-minimum outcome
     * is allowed to retain the renderer ceiling as terminal presentation
     * policy.
     *
     * Keeping the old ceiling merely because the allocated population is
     * richer re-enters the same global-ordinal walk.  On large multi-object
     * scenes that can alternate allocation and reconciliation indefinitely
     * even though the occurrence-local handoff has already succeeded. */
    HandoffCompletion completeOccurrenceHandoff(
	const BObolLodAdmissionRevisionStamp &stamp, int rendererCeiling)
    {
	HandoffCompletion completion;
	completion.completedRejectedConstraint =
	    this->completeReconciliation();
	const bool retainRendererCeiling =
	    this->retainsRendererCeilingFor(stamp, rendererCeiling);
	completion.releaseRendererCeiling = rendererCeiling >= 0 &&
	    !retainRendererCeiling;
	completion.measureCeilingFreeCandidate = this->probing() &&
	    completion.releaseRendererCeiling;
	return completion;
    }

    bool reject(const Constraint &constraint)
    {
	if (!constraint.valid() || !this->inProgress())
	    return false;
	this->constraintValue = constraint;
	this->acceptanceValue = Acceptance();
	this->stateValue = State::RECONCILING;
	this->sampledCeilingValue = -1;
	return true;
    }

    /* A complete occurrence allocation may prove, without another rendered
     * trial, that its protected minimum exceeds the hard presentation
     * allowance.  Record that prediction directly as the terminal constraint.
     * Routing it through reject() silently does nothing from IDLE and lets the
     * identical over-budget allocation rearm its handoff forever. */
    bool constrain(const Constraint &constraint)
    {
	if (!constraint.valid() || this->stateValue != State::IDLE)
	    return false;
	this->constraintValue = constraint;
	this->acceptanceValue = Acceptance();
	this->stateValue = State::REJECTED;
	this->sampledCeilingValue = -1;
	return true;
    }

    /* The retained framebuffer may itself be the only presentation that
     * fits the hard deadline when an occurrence-local allocation cannot
     * encode its renderer-wide cut.  That case has no reconciliation effect
     * left to run: retain the exact completed presentation as the terminal
     * constraint instead of manufacturing an ownerless RECONCILING state. */
    bool constrainPresented(const Constraint &constraint)
    {
	if (!constraint.valid() || this->stateValue == State::ACCEPTED ||
	    this->stateValue == State::REJECTED)
	    return false;
	this->constraintValue = constraint;
	this->acceptanceValue = Acceptance();
	this->stateValue = State::REJECTED;
	this->sampledCeilingValue = -1;
	return true;
    }

    bool acceptFractionalCeiling(const Acceptance &acceptance)
    {
	if (!acceptance.valid() || !this->probing())
	    return false;
	this->acceptanceValue = acceptance;
	this->constraintValue = Constraint();
	this->stateValue = State::ACCEPTED;
	this->sampledCeilingValue = -1;
	return true;
    }
    bool probing(void) const
    {
	return this->stateValue == State::PROBING;
    }

    bool inProgress(void) const
    {
	return this->stateValue == State::PROBING ||
	    this->stateValue == State::RECONCILING;
	}

    bool accepted(void) const
    {
	return this->stateValue == State::ACCEPTED;
	}

    bool capacityRejected(void) const
    {
	return this->stateValue == State::RECONCILING ||
	    this->stateValue == State::REJECTED;
	}

    bool blocksNewTrial(void) const
    {
	return this->stateValue != State::IDLE;
	}

    bool usesStaticDeadline(void) const
    {
	return this->stateValue == State::PROBING ||
	    this->stateValue == State::RECONCILING ||
	    this->stateValue == State::ACCEPTED ||
	    this->stateValue == State::REJECTED;
    }

    bool rejectedFor(const BObolLodAdmissionRevisionStamp &stamp) const
    {
	if (this->stateValue != State::REJECTED ||
	    !this->constraintValue.valid())
	    return false;
	const BObolLodAdmissionRevisionStamp &rejected =
	    this->constraintValue.revisionStamp;
	/* The completed frame may itself advance the capacity ledger.  Geometry,
	 * availability, camera, or policy changes alter the represented
	 * population and must obtain a new constraint proof. */
	return rejected.inventory == stamp.inventory &&
	    rejected.availability == stamp.availability &&
	    rejected.view == stamp.view && rejected.policy == stamp.policy;
    }

    /* An occurrence-local allocation cannot always reproduce a renderer-wide
     * cut: the minimum cut of many independent meshes may itself exceed the
     * completed framebuffer's hard presentation allowance.  In that case the
     * exact framebuffer and its reversible renderer ceiling are the terminal
     * presentation for this semantic epoch.  Keep this exception narrow and
     * revision-bound; every other ceiling remains reconciliation debt. */
    bool retainsRendererCeilingFor(
	const BObolLodAdmissionRevisionStamp &stamp, int ceiling) const
    {
	return ceiling >= 0 && this->rejectedFor(stamp) &&
	    this->constraintValue.reason == ConstraintReason::PROTECTED_MINIMUM &&
	    this->constraintValue.committedCeiling == ceiling;
    }

    const Constraint &constraint(void) const
    {
	return this->constraintValue;
	}

    const Acceptance &acceptance(void) const
    {
	return this->acceptanceValue;
    }

    size_t acceptedPresentationCostFor(
	const BObolLodAdmissionRevisionStamp &stamp) const
    {
	if (!this->accepted() || !this->acceptanceValue.valid())
	    return 0;
	const BObolLodAdmissionRevisionStamp &accepted =
	    this->acceptanceValue.revisionStamp;
	/* Completed-frame capacity remains valid when the capacity ledger records
	 * that very frame.  Inventory, availability, view, or policy changes alter
	 * the submitted population and must obtain a new proof. */
	return accepted.inventory == stamp.inventory &&
	    accepted.availability == stamp.availability &&
	    accepted.view == stamp.view && accepted.policy == stamp.policy ?
		this->acceptanceValue.presentedCost : 0;
    }

    size_t constrainedPresentationBudgetFor(
	const BObolLodAdmissionRevisionStamp &stamp) const
    {
	/* A completed static rejection is a terminal presentation proof for the
	 * same immutable scene/view problem.  Its reconciled occurrence budget
	 * supersedes the older steady-cadence capacity certificate; otherwise the
	 * next ordinary plan can restore that certificate's cheaper predecessor
	 * and reopen the static handoff indefinitely. */
	return this->rejectedFor(stamp) ? this->constraintValue.allowedCost : 0;
    }

    int sampledCeiling(void) const
    {
	return this->sampledCeilingValue;
    }

    float baselinePointProxyPixelThreshold(void) const
    {
	return this->baselinePointProxyPixelThresholdValue;
    }

    void noteSampledCeiling(int ceiling)
    {
	if (this->stateValue == State::PROBING)
	    this->sampledCeilingValue = ceiling;
    }

    void resetSample(void)
    {
	this->sampledCeilingValue = -1;
    }

private:
    static float sanitizePointProxyPixelThreshold(float threshold)
    {
	if (!std::isfinite(threshold))
	    return 1.0f;
	return std::max(1.0f, std::min(64.0f, threshold));
    }

    enum class State : uint8_t {
	IDLE = 0,
	PROBING,
	ACCEPTED,
	RECONCILING,
	REJECTED
    };

    State stateValue = State::IDLE;
    int sampledCeilingValue = -1;
    float baselinePointProxyPixelThresholdValue = 1.0f;
    Constraint constraintValue;
    Acceptance acceptanceValue;
};

static_assert(std::is_trivially_copyable<
    BObolLodStaticQualityTrial>::value,
    "static-quality trial state must remain an allocation-free value");


/*
 * Renderer-independent presentation-quality calculations.  These functions
 * convert completed-frame evidence into a camera-local pixel error, small-
 * occurrence aggregation threshold, and reversible PoP prefix ceiling.  They
 * accept only scalar observations so the coordinator cannot acquire a scene
 * pointer or accidentally turn motion feedback into a hierarchy traversal.
 */
class BObolLodQualityPolicy {
public:
    static constexpr long double StaticDeadlineHeadroom = 0.95L;

    enum class DeadlineSuccessor : uint8_t {
	RETRY_TRANSACTION = 0,
	CONTINUE_POPULATION,
	RECOVER_PRESENTATION
    };

    struct StableCapacityEvidence {
	bool proven = false;
	long double renderCostPerSecond = 0.0L;
    };

    /* The first quiet traversal of a newly published retained population may
     * pay command construction, software-code warmup, or immutable-buffer
     * installation which precedes an exact finite-work certificate.  A hard
     * abort is still a responsiveness boundary, but it is not
     * yet proof that the retained cuts themselves are unsustainable.  Permit
     * exactly one unchanged retry while a typed population/frame obligation
     * or an explicit static-quality trial identifies that first presentation.
     * A second abort is ordinary capacity evidence.  Interactive input never
     * waits for this retry. */
    static bool retryTransientPresentation(bool interactive,
	unsigned int consecutiveInterruptedFrames,
	bool preparationAdvanced, bool publicationFramePending,
	bool refinementFramePending, bool pointCalibrationPending,
	bool staticQualityTrial)
    {
	if (preparationAdvanced)
	    return true;
	return !interactive && consecutiveInterruptedFrames == 1 &&
	    (publicationFramePending || refinementFramePending ||
	     pointCalibrationPending || staticQualityTrial);
    }

    /* A deadline is an observation, not an unconditional retry edge.  A
     * transaction which reduced its finite preparation obligation owns an
     * immediate replay.  Otherwise an existing population producer keeps its
     * cursor and supplies the next presentation edge when that work advances.
     * Capacity recovery is valid only when neither owner applies. */
    static DeadlineSuccessor deadlineSuccessor(
	bool transactionRetry, bool populationWorkPending)
    {
	if (transactionRetry)
	    return DeadlineSuccessor::RETRY_TRANSACTION;
	if (populationWorkPending)
	    return DeadlineSuccessor::CONTINUE_POPULATION;
	return DeadlineSuccessor::RECOVER_PRESENTATION;
    }

    /* A terminal static-quality certificate owns the allowance under which
     * its occurrence-local population was selected.  Transitional coverage
     * work may use a smaller extended quiet deadline only when no static
     * trial or terminal certificate is active.  Reversing this precedence
     * makes the smaller deadline reject its already-certified population and
     * reopens capacity recovery indefinitely. */
    static uint64_t quietPresentationDeadline(
	uint64_t staticQualityDeadline,
	uint64_t structuralRepairDeadline,
	bool structuralRepairPending,
	bool staticQualityOwnsDeadline)
    {
	if (staticQualityOwnsDeadline || !structuralRepairPending)
	    return staticQualityDeadline;
	return structuralRepairDeadline;
    }

    static bool quietPresentationIrreducible(bool interactive,
	bool presentationAvailable, int progressiveCeiling,
	size_t activeRenderCost, size_t minimumRenderCost,
	size_t progressivePayloadCount, bool pointTransitionPending)
    {
	return !interactive && presentationAvailable &&
	    (progressiveCeiling == 0 || progressivePayloadCount == 0 ||
	     (activeRenderCost > 0 && activeRenderCost <= minimumRenderCost)) &&
	    !pointTransitionPending;
    }

    static float staticFrameTargetFps(float preferredFps,
	uint64_t hardDeadlineNanoseconds)
    {
	if (!std::isfinite(preferredFps) || preferredFps <= 0.0f ||
	    !hardDeadlineNanoseconds)
	    return preferredFps;
	const long double deadlineFps = 1000000000.0L /
	    static_cast<long double>(hardDeadlineNanoseconds);
	if (!std::isfinite(deadlineFps) || deadlineFps <= 0.0L)
	    return preferredFps;
	return std::min(preferredFps, static_cast<float>(deadlineFps));
    }

    /* Convert one exact completed static presentation into a conservative
     * upper bound for the next discrete PoP population.  Retain five percent
     * scheduling headroom so harmless timing jitter does not turn a stable
     * framebuffer into a speculative deadline miss.  The current completed
     * work is always affordable by construction; callers use equality as the
     * proof that no richer cached population should be attempted. */
    static size_t staticPresentationRenderCostLimit(size_t presentedCost,
	uint64_t renderNanoseconds, uint64_t hardDeadlineNanoseconds)
    {
	if (!presentedCost || !renderNanoseconds ||
	    !hardDeadlineNanoseconds ||
	    renderNanoseconds > hardDeadlineNanoseconds)
	    return 0;
	const long double predicted =
	    static_cast<long double>(presentedCost) *
	    static_cast<long double>(hardDeadlineNanoseconds) *
	    StaticDeadlineHeadroom /
	    static_cast<long double>(renderNanoseconds);
	if (!std::isfinite(predicted) ||
	    predicted >= static_cast<long double>(SIZE_MAX))
	    return SIZE_MAX;
	return std::max(presentedCost, static_cast<size_t>(predicted));
    }

    /* A renderer-wide ceiling may fit the ordinary quiet cadence even though
     * the cheapest complete occurrence-local population does not.  A static
     * event-driven view can try that local minimum once under its longer hard
     * deadline.  Include conservatively extrapolated headroom, but never
     * exceed current pixel demand.  The endpoint deadline independently
     * rejects an optimistic estimate. */
    static size_t staticLocalMinimumRetryBudget(size_t presentedCost,
	size_t selectedLocalCost, size_t pixelDemandCost,
	uint64_t renderNanoseconds, uint64_t hardDeadlineNanoseconds)
    {
	if (!selectedLocalCost)
	    return 0;
	size_t budget = std::max(selectedLocalCost,
	    staticPresentationRenderCostLimit(presentedCost,
		renderNanoseconds, hardDeadlineNanoseconds));
	if (pixelDemandCost >= selectedLocalCost)
	    budget = std::min(budget, pixelDemandCost);
	return budget;
    }

    /* Correct an interrupted presentation in one cost-domain step.
     * The same five-percent scheduling margin used for successful headroom
     * probes is sufficient here: the elapsed/deadline ratio already supplies
     * the overload correction.  A second fixed twenty-percent penalty made a
     * harmless timer-edge miss discard exact-view capacity evidence and
     * visibly over-coarsen large meshes.  Interactive callers supply their
     * input-frame deadline; quiet callers supply their static deadline. */
    static size_t deadlineRecoveryCostLimit(size_t attemptedCost,
	uint64_t elapsedNanoseconds, uint64_t hardDeadlineNanoseconds)
    {
	if (!attemptedCost || !elapsedNanoseconds ||
		!hardDeadlineNanoseconds)
	    return 0;
	const long double predicted =
	    static_cast<long double>(attemptedCost) *
	    static_cast<long double>(hardDeadlineNanoseconds) *
	    StaticDeadlineHeadroom /
	    static_cast<long double>(elapsedNanoseconds);
	if (!std::isfinite(predicted) || predicted <= 0.0L)
	    return 0;
	if (predicted >= static_cast<long double>(SIZE_MAX))
	    return SIZE_MAX;
	return std::min(attemptedCost,
	    static_cast<size_t>(predicted));
    }

    /* Return the smallest scene-cost allowance which can represent the next
     * pixel-error tier under the same inverse-square work model used by
     * stablePixelError().  The value is a floor for one explicitly selected
     * tier, not permission to spend the complete capacity extrapolated from a
     * fast frame.  Keeping those two quantities distinct prevents a stale
     * streaming budget from rejecting proven static quality while avoiding an
     * unnecessarily large allocator jump. */
    static size_t pixelErrorRenderCostFloor(size_t activeCost,
	float currentPixelError, float nextPixelError)
    {
	if (!activeCost || !std::isfinite(currentPixelError) ||
	    !std::isfinite(nextPixelError) || currentPixelError <= 0.0f ||
	    nextPixelError <= 0.0f ||
	    nextPixelError + 1.0e-6f >= currentPixelError)
	    return activeCost;
	const long double ratio =
	    static_cast<long double>(currentPixelError) /
	    static_cast<long double>(nextPixelError);
	const long double predicted =
	    static_cast<long double>(activeCost) * ratio * ratio;
	if (!std::isfinite(predicted) ||
	    predicted >= static_cast<long double>(SIZE_MAX))
	    return SIZE_MAX;
	return std::max(activeCost,
	    static_cast<size_t>(std::ceil(predicted)));
    }

    /* Translate completed-frame headroom into the allocator's active-mesh
     * cost domain.  Presented work also contains the current frame's fixed
     * point-proxy, occurrence, and batch overhead; refining mesh prefixes does
     * not multiply that base.  Only the demonstrated excess above the exact
     * presented population is therefore added to the active-mesh allowance. */
    static size_t incrementalSceneCostBudget(size_t activeCost,
	size_t presentedCost, size_t presentationCostLimit)
    {
	if (!activeCost || !presentedCost ||
	    presentationCostLimit <= presentedCost)
	    return activeCost;
	const size_t additional = presentationCostLimit - presentedCost;
	return additional > SIZE_MAX - activeCost ?
	    SIZE_MAX : activeCost + additional;
    }

    /* An exact first-use framebuffer is not reusable timing evidence.  The
     * event-driven host must schedule one unchanged CAD replay before it may
     * classify a terminal static-quality cut.  Keeping this predicate in the
     * scalar policy makes the missing-sample successor explicit and testable;
     * the controller remains responsible only for arming its existing
     * presentation owner. */
    static bool staticQualityTimingReplayRequired(bool stableContext,
	bool exactCompletedFrame, uint64_t reusableCadNanoseconds)
    {
	return stableContext && exactCompletedFrame &&
	    reusableCadNanoseconds == 0;
    }

    static float interactivePixelError(uint64_t renderNanoseconds,
	float targetFps)
    {
	const double targetNanoseconds = targetFrameNanoseconds(targetFps);
	if (!renderNanoseconds)
	    return 4.0f;
	if (static_cast<double>(renderNanoseconds) <=
	    targetNanoseconds * 1.05)
	    return 1.0f;
	/* Surface work grows approximately with inverse squared screen error.
	 * Scaling the error by sqrt(observed / target) is therefore the direct
	 * correction supported by the completed frame.  The former additional
	 * factor of two made progressiveCeiling() discard one extra PoP
	 * population on every feedback sample.  During a held gesture that error
	 * compounded once per mouse event: Lucy reached a 6--10 ms cut while the
	 * configured target required only 16.7 ms.  Retain a minimum factor of two
	 * after a real miss so a discrete hierarchy always makes progress. */
	const double scale = std::sqrt(
	    static_cast<double>(renderNanoseconds) / targetNanoseconds);
	return static_cast<float>(std::max(2.0,
	    std::min(64.0, scale)));
    }

    /*
     * One pixel is the fast convergence contract, not necessarily the best
     * quiet terminal image.  A small scene can otherwise stop with visibly
     * faceted silhouettes while using a tiny fraction of both its measured
     * draw allowance and resident-memory allowance.  Admit 0.75-, 0.5-, and
     * 0.25-pixel quality tiers only from an exact, completed quiet-frame
     * witness.  Surface work scales approximately with inverse squared pixel
     * error, so a 1.0 -> 0.75 step requires about 1.78 times the work while a
     * 1.0 -> 0.5 step requires about four times the work.
     *
     * The quarter-pixel tier addresses a different failure mode than ordinary
     * faceting.  A PoP vertex displacement can be subpixel while a long,
     * skinny triangle collapses across its short axis, producing a conspicuous
     * one-pixel slit along a silhouette.  It is considered only after a
     * completed half-pixel frame.  That separate pass re-proves both measured
     * rendering capacity and the caller's one-quarter resident-memory
     * headroom, so its estimated fourfold growth fits the same bounded scene
     * and memory allowances.  We deliberately do not jump directly from one
     * pixel to one quarter pixel: the intermediate resident population is the
     * capacity witness.
     *
     * A first-use frame may include buffer preparation or upload; if that more
     * expensive frame still proves the required margin it is conservative
     * evidence, not a reason for cold and warm caches to converge differently.
     * The test is deliberately expressed in total scene-cost currency rather
     * than an object count or a model-specific cut ordinal.
     *
     * The caller owns the memory-headroom proof because resident limits and
     * reservations belong to the service.  Returning a scalar keeps this
     * policy allocation-free and directly testable.  Above one pixel, a
     * static view may use the largest safe refinement admitted by the same
     * inverse-square model.  Motion can legitimately leave a large object at
     * several pixels; treating that value as a terminal tier made a stopped
     * close zoom retain its interactive approximation even when the retained
     * framebuffer deadline had measured headroom.  Below one pixel the
     * discrete 0.75, 0.5, and 0.25 tiers keep the raster-stable terminal
     * contract.  Neither path is permission to load full geometry.
     */
    static float stablePixelError(float currentPixelError,
	size_t activeCost, size_t sceneBudget,
	uint64_t renderNanoseconds, float targetFps,
	bool exactCompletedFrame, bool residentMemoryHeadroom)
    {
	if (!std::isfinite(currentPixelError) ||
	    currentPixelError <= 0.2501f ||
	    !activeCost || !sceneBudget ||
	    !renderNanoseconds || !std::isfinite(targetFps) ||
	    targetFps <= 0.0f || !exactCompletedFrame ||
	    !residentMemoryHeadroom)
	    return currentPixelError;

	const long double targetNanoseconds =
	    1000000000.0L / static_cast<long double>(targetFps);
	const auto affordable = [&](float candidate) {
	    if (candidate + 1.0e-6f >= currentPixelError)
		return false;
	    const long double ratio =
		static_cast<long double>(currentPixelError) /
		static_cast<long double>(candidate);
	    const long double workMultiplier = ratio * ratio;
	    const bool timeHeadroom =
		static_cast<long double>(renderNanoseconds) * workMultiplier <=
		targetNanoseconds;
	    const bool costHeadroom = sceneBudget == SIZE_MAX ||
		static_cast<long double>(activeCost) * workMultiplier <=
		static_cast<long double>(sceneBudget);
	    return timeHeadroom && costHeadroom;
	};

	/* A quiet frame which is still coarser than one pixel is an interactive
	 * carry-over, not a canonical quality rung.  Find the finest pixel error
	 * justified by both independent measurements.  Rounding toward the
	 * coarser representable value preserves the proof when the float result is
	 * passed back through affordable(). */
	if (currentPixelError > 1.0f) {
	    const long double timeBound =
		static_cast<long double>(currentPixelError) * std::sqrt(
		    static_cast<long double>(renderNanoseconds) /
		    targetNanoseconds);
	    const long double costBound = sceneBudget == SIZE_MAX ? 0.0L :
		static_cast<long double>(currentPixelError) * std::sqrt(
		    static_cast<long double>(activeCost) /
		    static_cast<long double>(sceneBudget));
	    const long double minimumCandidate = std::max(1.0L,
		std::max(timeBound, costBound));
	    if (!std::isfinite(minimumCandidate) ||
		minimumCandidate >= static_cast<long double>(currentPixelError))
		return currentPixelError;
	    float candidate = static_cast<float>(minimumCandidate);
	    if (candidate + 1.0e-6f < currentPixelError && affordable(candidate))
		return candidate;
	    candidate = std::nextafter(candidate,
		std::numeric_limits<float>::infinity());
	    return candidate + 1.0e-6f < currentPixelError && affordable(candidate) ?
		candidate : currentPixelError;
	}

	/* Require a separate completed half-pixel witness before attempting the
	 * fourfold quarter-pixel population.  Besides making the draw estimate
	 * current, this lets the caller re-evaluate resident-memory headroom after
	 * the half-pixel suffix has actually arrived. */
	if (currentPixelError <= 0.5001f)
	    return affordable(0.25f) ? 0.25f : currentPixelError;
	if (affordable(0.5f))
	    return 0.5f;
	if (affordable(0.75f))
	    return 0.75f;
	return currentPixelError;
    }

    static float pointProxyThreshold(float currentThreshold,
	uint64_t renderNanoseconds, float targetFps)
    {
	const double targetNanoseconds = targetFrameNanoseconds(targetFps);
	const double observedThreshold =
	    static_cast<double>(currentThreshold);
	const double current = std::isfinite(observedThreshold) ?
	    std::max(1.0, std::min(64.0, observedThreshold)) : 1.0;
	if (!renderNanoseconds ||
	    static_cast<double>(renderNanoseconds) <=
		targetNanoseconds * 1.05)
	    return static_cast<float>(current);

	/* Apply the measured correction to the cut which produced the frame.
	 * Recomputing an absolute threshold reaches a false fixed point when
	 * per-object command work dominates. */
	const double correction = std::sqrt(
	    static_cast<double>(renderNanoseconds) / targetNanoseconds);
	const double absoluteSeed = interactivePixelError(
	    renderNanoseconds, targetFps);
	return static_cast<float>(std::min(64.0,
	    std::max(absoluteSeed, current * correction)));
    }

    static int progressiveCeiling(int activeMaximum, float pixelError,
	int currentCeiling = -1)
    {
	if (activeMaximum < 0)
	    return -1;
	const double observedError = static_cast<double>(pixelError);
	const double error = std::isfinite(observedError) ?
	    std::max(1.0, std::min(64.0, observedError)) :
	    (observedError > 0.0 ? 64.0 : 1.0);
	if (error <= 1.0 + std::numeric_limits<double>::epsilon())
	    return -1;

	/* PoP populations are discrete.  Round toward the known-safe coarser
	 * population, and apply repeated feedback relative to the ceiling which
	 * actually produced the measurement. */
	const int reduction = std::max(1,
	    static_cast<int>(std::ceil(std::log2(error))));
	const int measuredMaximum = currentCeiling >= 0 ?
	    std::min(activeMaximum, currentCeiling) : activeMaximum;
	return std::max(0, measuredMaximum - reduction);
    }

    static StableCapacityEvidence stableCapacityEvidence(
	size_t presentedRenderCost, uint64_t renderNanoseconds,
	float stableTargetFps, bool exactScaleQualityFrame)
    {
	StableCapacityEvidence evidence;
	if (!exactScaleQualityFrame || !presentedRenderCost ||
	    !renderNanoseconds || !std::isfinite(stableTargetFps) ||
	    stableTargetFps <= 0.0f)
	    return evidence;
	const long double deadline = 1000000000.0L /
	    static_cast<long double>(stableTargetFps);
	if (static_cast<long double>(renderNanoseconds) > deadline)
	    return evidence;
	const long double throughput =
	    static_cast<long double>(presentedRenderCost) * 1000000000.0L /
	    static_cast<long double>(renderNanoseconds);
	if (!std::isfinite(throughput) || throughput <= 0.0L)
	    return evidence;
	evidence.proven = true;
	evidence.renderCostPerSecond = throughput;
	return evidence;
    }

    /*
     * A globally ordered PoP prefix can have a discontinuous next population:
     * fixed per-frame overhead may make the current cut appear expensive when
     * extrapolated linearly even though the next complete prefix still fits.
     * Give one scene-global occurrence permission to measure that transition,
     * but bound the total experimental presentation by both four times the
     * active work and twice the predicted scene budget.  The returned value is
     * only the excess beyond the ordinary budget; the submit action consumes
     * its remaining ordinary allowance first and charges the full transition.
     */
    static size_t discreteTrialOverBudgetAllowance(
	size_t activeCost, size_t currentBudget)
    {
	if (!activeCost || !currentBudget || currentBudget == SIZE_MAX ||
	    activeCost > currentBudget)
	    return 0;
	const size_t activeLimit = activeCost > SIZE_MAX / 4 ?
	    SIZE_MAX : activeCost * 4;
	const size_t budgetLimit = currentBudget > SIZE_MAX / 2 ?
	    SIZE_MAX : currentBudget * 2;
	const size_t trialTotal = std::min(activeLimit, budgetLimit);
	return trialTotal > currentBudget ?
	    trialTotal - currentBudget : 0;
    }

private:
    static double targetFrameNanoseconds(float targetFps)
    {
	return std::isfinite(targetFps) && targetFps > 0.0f ?
	    1000000000.0 / static_cast<double>(targetFps) :
	    1000000000.0 / 60.0;
    }
};

/*
 * Geometry production is a pipeline, not just the worker pool.  A bounded
 * owner-thread submission cursor or database provider can be the only
 * evidence that another worker wave will exist.  Treating a momentarily empty
 * worker pool as end-of-stream publishes tiny result batches and serializes
 * large-scene realization behind whole-scene renders.
 *
 * Presentation-owned calibration may deliberately pause submission until a
 * changed cut has rendered.  Such a cursor cannot be counted as its own
 * future-frame witness: doing so closes a wait cycle in which publication
 * waits for submission and submission waits for publication.  Likewise, a
 * publication whose frame has already been requested is the presentation
 * obligation itself, not an independent producer of a later frame.  Only a
 * publication batch still waiting to request its frame is a future witness.
 */

#endif /* LIBBOBOL_LOD_VIEW_POLICY_PRIVATE_H */
