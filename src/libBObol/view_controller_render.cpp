/*       V I E W _ C O N T R O L L E R _ R E N D E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * @file view_controller_render.cpp
 *
 * Executes bounded render passes and owns frame interruption, timing,
 * static-frame reuse, image readback, and render-request bookkeeping.
 */

#include "common.h"

#include "bu/log.h"
#include "bu/str.h"
#include "bu/datetime.h"

#include "bv.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BEditPreview.h"
#include "BObol/BExportAction.h"
#include "BObol/BHUDLabelOverlay.h"
#include "BObol/BLineLayerOverlay.h"
#include "BObol/BLodService.h"
#include "BObol/BLodUpdateAction.h"
#include "BObol/BMeasureAction.h"
#include "BObol/BMeshLodSubmitAction.h"
#include "BObol/BMeshResidencyAction.h"
#include "BObol/BMeshShape.h"
#include "BObol/BPickDetail.h"
#include "BObol/BSceneController.h"
#include "BObol/BSnapAction.h"
#include "BObol/BViewAttachment.h"
#include "BObol/BViewController.h"
#include "BObol/BViewLod.h"
#include "BObol/BViewStore.h"
#include "cad_assembly_private.h"
#include "lod_coordinator_private.h"
#include "retained_allocation_private.h"
#include "view_controller_private.h"
#include "raytrace.h"
#include "rt/view.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <new>
#include <queue>
#include <string.h>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Inventor/SbName.h>
#include <Inventor/SbPlane.h>
#include <Inventor/SoDB.h>
#include <Inventor/SoRenderManager.h>
#include <Inventor/SoViewport.h>
#include <Inventor/SoOffscreenRenderer.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/elements/SoGLCacheContextElement.h>
#include <Inventor/gl.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoClipPlane.h>
#include <Inventor/nodes/SoDepthBuffer.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoEnvironment.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoLight.h>
#include <Inventor/nodes/SoLightModel.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoPointLight.h>
#include <Inventor/nodes/SoSpotLight.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
struct ControllerRenderDeadlineContext {
    BObolViewController *controller = NULL;
    uint64_t deadlineNanoseconds = 0;
    SoGLRenderAction::SoGLRenderAbortCB *previous = NULL;
    void *previousData = NULL;
};

static SoGLRenderAction::AbortCode
controller_render_deadline_cb(void *userData)
{
    ControllerRenderDeadlineContext *context =
	static_cast<ControllerRenderDeadlineContext *>(userData);
    if (!context)
	return SoGLRenderAction::CONTINUE;
    if (context->previous) {
	const SoGLRenderAction::AbortCode prior =
	    (*context->previous)(context->previousData);
	if (prior != SoGLRenderAction::CONTINUE)
	    return prior;
    }
    if (!context->controller || !context->deadlineNanoseconds)
	return SoGLRenderAction::CONTINUE;
    return context->controller->beginRenderTiming() >=
	context->deadlineNanoseconds ?
	SoGLRenderAction::ABORT : SoGLRenderAction::CONTINUE;
}

SbBool
BObolViewController::renderPending(SbBool clearWindow,
				     SbBool clearZBuffer,
				     SbString *reason)
{
    /* The caller owns and binds the render context.  In particular, an idle
     * poll must not advance work which can enqueue a frame and then enter Coin
     * on a context-free caller.  Progressive providers publish their own
     * frame request (and wake the host through the request callback), so an
     * already queued frame is the sole authority to begin preparation here. */
    if (!this->isRenderRequested())
	return FALSE;

    /* LoD off -> classic behavior: realize the whole scene before presenting.
     * LoD on -> stay on the progressive coarse-first path (no whole-tree
     * realize here; geometry streams in via advanceProgressiveWork). */
    if (this->isForceRealizeDisplay())
	(void)this->realizePending();

    /* A retained selection, highlight, or faceplate mutation requests a
     * presentation frame without requesting a geometry-capacity sample.  Do
     * not use that repaint as an opportunity to run an otherwise idle LoD
     * coordinator: a mouse-drag overlay could then reopen stable headroom or
     * quality policy and make an unchanged scene report that it was refining
     * again.  Real progressive work and capacity-relevant renders still get
     * the traditional pre-render pump, including the case where a
     * presentation request was merged with either one. */
    const BObolHostWorkSnapshot preRenderWork =
	this->getHostWorkSnapshot();
    if (preRenderWork.pumpPending() ||
	preRenderWork.capacitySampleRequested()) {
	(void)this->advanceProgressiveWork(NULL, NULL);
    }
    this->synchronizePresentation();


    if (!this->d->renderManager || !this->getRenderContextManager() ||
	!this->d->activeCamera || !this->getRenderRoot())
	return FALSE;

    SbString renderReasonCopy;
    SbBool lodCapacityRelevant = TRUE;
    if (!this->consumeRenderRequest(&renderReasonCopy,
	    &lodCapacityRelevant))
	return FALSE;

    if (reason)
	*reason = renderReasonCopy;

    const BObolViewLodState *presentationState = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;

    /* Full force-realize traversal is not an LoD capacity sample.  In
     * particular, switching LoD off must neither learn a scene budget from
     * the deliberately complete population nor arm an LoD deadline-recovery
     * barrier. */
    if (!this->isLodPresentationCapacityRelevant())
	lodCapacityRelevant = FALSE;

    const uint64_t started = this->beginRenderTiming();
    SoGLRenderAction *renderAction =
	this->d->renderManager->getGLRenderAction();
    ControllerRenderDeadlineContext deadlineContext;
    const uint64_t deadlineDuration = lodCapacityRelevant ?
	this->getCurrentPresentationFrameDeadline() : 0;
    if (renderAction && deadlineDuration) {
	deadlineContext.controller = this;
	deadlineContext.deadlineNanoseconds =
	    started > UINT64_MAX - deadlineDuration ? UINT64_MAX :
	    started + deadlineDuration;
	renderAction->getAbortCallback(
	    deadlineContext.previous, deadlineContext.previousData);
	renderAction->setAbortCallback(
	    controller_render_deadline_cb, &deadlineContext);
    }
    if (clearWindow)
	this->renderBackground();
    else if (clearZBuffer)
	glClear(GL_DEPTH_BUFFER_BIT);
    const uint64_t backgroundCompleted = this->beginRenderTiming();
    const uint64_t cadExecutionBefore = presentationState ?
	presentationState->cadPresentationExecutionSerial() : 0;
    const uint64_t cadPreparationBefore = presentationState ?
	presentationState->cadPresentationPreparationSerial() : 0;
    if (presentationState)
	presentationState->beginCadPresentationFrame();
    this->d->renderManager->render(static_cast<SbBool>(FALSE), static_cast<SbBool>(FALSE));
    const uint64_t sceneCompleted = this->beginRenderTiming();
    const uint64_t cadExecutionAfter = presentationState ?
	presentationState->cadPresentationExecutionSerial() : 0;
    const uint64_t cadPreparationAfter = presentationState ?
	presentationState->cadPresentationPreparationSerial() : 0;
    if (presentationState)
	presentationState->refreshCadPresentationFrameStatus();
    const SbBool cadFrameIncomplete = presentationState &&
	!presentationState->lastCadPresentationFrameExact();
    const SbBool interrupted =
	(renderAction && renderAction->hasTerminated()) || cadFrameIncomplete;
    if (renderAction && deadlineDuration)
	renderAction->setAbortCallback(
	    deadlineContext.previous, deadlineContext.previousData);
    this->d->lastBackgroundRenderTimeNanoseconds =
	backgroundCompleted >= started ?
	    backgroundCompleted - started : 0;
    this->d->lastSceneRenderTimeNanoseconds =
	sceneCompleted >= backgroundCompleted ?
	    sceneCompleted - backgroundCompleted : 0;
    if (interrupted) {
	this->notePresentationRenderInterrupted(
	    sceneCompleted > started ? sceneCompleted - started : 1,
	    cadExecutionAfter != cadExecutionBefore ? TRUE : FALSE,
	    cadPreparationAfter != cadPreparationBefore ? TRUE : FALSE,
	    lodCapacityRelevant);
	return FALSE;
    }
    this->completeRenderTiming(started, lodCapacityRelevant);
    return TRUE;
}

uint64_t
BObolViewController::beginRenderTiming(void) const
{
    return static_cast<uint64_t>(
	std::chrono::duration_cast<std::chrono::nanoseconds>(
	    std::chrono::steady_clock::now().time_since_epoch()).count());
}

void
BObolViewController::setPresentationFrameDeadlines(
    uint64_t interactiveNanoseconds, uint64_t stableNanoseconds)
{
    if (this->d->interactivePresentationFrameDeadlineNanoseconds ==
	    interactiveNanoseconds &&
	this->d->stablePresentationFrameDeadlineNanoseconds ==
	    stableNanoseconds)
	return;
    this->d->resetLodViewQualityHistory();
    this->d->interactivePresentationFrameDeadlineNanoseconds =
	interactiveNanoseconds;
    this->d->stablePresentationFrameDeadlineNanoseconds =
	stableNanoseconds;
    this->d->consecutiveInterruptedPresentationFrames = 0;
}

uint64_t
BObolViewController::getInteractivePresentationFrameDeadline(void) const
{
    return this->d->interactivePresentationFrameDeadlineNanoseconds;
}

uint64_t
BObolViewController::getStablePresentationFrameDeadline(void) const
{
    return this->d->stablePresentationFrameDeadlineNanoseconds;
}

uint64_t
BObolViewController::getCurrentPresentationFrameDeadline(void) const
{
    /* A force-realize view has no legal progressive presentation to fall
     * back to.  Applying the LoD quality deadline in that state can only
     * abort and retry the identical complete traversal forever (hidden-line
     * Hubble is a representative 103 ms frame against the ordinary 100 ms
     * quiet deadline).  LoD-off is an explicit request for the complete
     * representation, so let that traversal finish and exclude it from LoD
     * capacity calibration.  Interactive responsiveness for very large
     * models is supplied by the managed LoD policy; disabling that policy
     * deliberately opts out of its bounded-frame guarantee. */
    if (this->isForceRealizeDisplay())
	return 0;

    /* Deterministic/offline convergence is explicitly outside the interactive
     * presentation contract.  Applying the deadline here can replace a fully
     * resident terminal prefix with a responsiveness ceiling between the
     * final progressive pump and its capture. */
    if (this->d->forceTerminalLodRefinement)
	return 0;

    /* Button-up on a pose-only orthographic gesture has already restored an
     * exact stable presentation.  Judge that one release frame by the stable
     * hard deadline while the 150 ms semantic-census debounce remains active;
     * applying the motion deadline here immediately hides the snapshot we
     * just restored and recreates a visible refinement staircase. */
    const bool restoredPosePresentation =
	this->d->lodInteractive && !this->d->lodGestureActive &&
	this->d->lodPresentationPolicy.priorRestored() &&
	!this->d->lodViewDemandPolicy.interactionScaleChanged() &&
	this->getActiveLodMeshPayloadCount() > 0;
    if ((!this->d->lodInteractive && !this->d->lodGestureActive) ||
	restoredPosePresentation)
	return this->d->stablePresentationFrameDeadlineNanoseconds;

    uint64_t base =
	this->d->interactivePresentationFrameDeadlineNanoseconds;
    /* An ordinary motion frame targets the configured interactive cadence.
     * A deliberate zoom-quality probe has a separate 10 Hz hard floor;
     * aborting it at the ordinary 40 ms software deadline prevents the frame
     * which would calibrate and retain a useful richer cut from completing. */
    if (this->d->lodViewDemandPolicy.qualityBudgetActive())
	base = std::max<uint64_t>(base,
	    BObolLodViewDemandPolicy::qualityFrameDurationNanoseconds());
    if (!base || !this->d->consecutiveInterruptedPresentationFrames)
	return base;

    /* A hard deadline without forward-progress backoff can starve the first
     * useful software frame: every richer coherent retry may cost 45-60 ms,
     * so a 40 ms limit keeps presenting the old background forever.  Permit
     * bounded 50% steps after consecutive aborts, capped by the quiet-frame
     * deadline.  The first completed frame resets this allowance. */
    const uint64_t steps = std::min<uint64_t>(
	this->d->consecutiveInterruptedPresentationFrames, 4u);
    const uint64_t increment = base / 2u;
    const uint64_t candidate =
	increment && steps <= (UINT64_MAX - base) / increment ?
	base + steps * increment : UINT64_MAX;
    const uint64_t stableCap =
	this->d->stablePresentationFrameDeadlineNanoseconds;
    return stableCap ? std::min(candidate, std::max(base, stableCap)) :
	candidate;
}

SbBool
BObolViewController::isLodPresentationCapacityRelevant(void) const
{
    if (this->isForceRealizeDisplay() ||
	this->d->forceTerminalLodRefinement)
	return FALSE;

    const BObolViewLodState *state = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;

    /* A deadline is useful only when the retained scene contains an active
     * LoD-managed population which recovery can make cheaper.  Evaluated or
     * direct geometry may remain comparatively expensive after a policy
     * transition (mode-0 Hubble wire is a representative example), but
     * aborting that immutable traversal cannot alter its next frame.  A later
     * progressive publication requests another frame and supplies nonzero
     * managed cost, making subsequent capacity feedback relevant. */
    return state && state->activeRenderCost() > 0 ? TRUE : FALSE;
}

void
BObolViewController::notePresentationRenderInterrupted(
    uint64_t elapsedNanoseconds, SbBool cadDrawAttempted,
    SbBool cadPreparationChanged, SbBool lodCapacityRelevant)
{
    if (!elapsedNanoseconds)
	return;
    this->d->interruptedPresentationFrameCount++;
    this->d->lastInterruptedPresentationTimeNanoseconds =
	elapsedNanoseconds;
    if (cadDrawAttempted)
	this->d->lodInterruptedPresentationReplayPending = TRUE;
    if (!lodCapacityRelevant) {
	/* A presentation patch may still encounter SoCADAssembly's resumable
	 * command preparation.  Finish that patch, but do not let it cancel a
	 * geometry headroom witness, lower a PoP ceiling, or alter a scene budget. */
	this->requestPresentationRender("render-presentation-replay");
	return;
    }
    this->d->consecutiveInterruptedPresentationFrames++;
    const BObolViewLodState *state = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    const bool restoredPosePresentation =
	this->d->lodInteractive && !this->d->lodGestureActive &&
	this->d->lodPresentationPolicy.priorRestored() &&
	!this->d->lodViewDemandPolicy.interactionScaleChanged() &&
	this->getActiveLodMeshPayloadCount() > 0;
    const SbBool interactive =
	(this->d->lodInteractive || this->d->lodGestureActive) &&
	!restoredPosePresentation;
    const int interruptedActiveMaximum = state ?
	state->maximumActiveProgressiveCut() : -1;
    const int interruptedPresentedMaximum =
	interruptedActiveMaximum >= 0 &&
	this->d->lodInteractiveProgressiveCeiling >= 0 ?
	    std::min(interruptedActiveMaximum,
		this->d->lodInteractiveProgressiveCeiling) :
	    interruptedActiveMaximum;
    const size_t interruptedActiveCost = state ?
	state->activeRenderCost() : 0;
    const size_t interruptedMinimumCost = state ?
	state->minimumActiveRenderCost() : interruptedActiveCost;
    /* With no renderer-wide ceiling, a current retained-allocation
     * certificate is the exact population requested from the presentation
     * layer, including point-aggregated occurrences and immutable scene
     * cost.  activeRenderCost() deliberately describes richer retained mesh
     * cuts instead.  Using that retained value to correct an interrupted
     * presentation overstates the attempted work by the aggregated tail and
     * produces a long sequence of insufficient recovery allocations. */
    const BObolRetainedAllocationResult &allocationCertificate =
	this->d->lodRetainedAllocationCertificate;
    const bool certifiedUnconstrainedPresentation = state &&
	this->d->lodInteractiveProgressiveCeiling < 0 &&
	allocationCertificate.allocationPlanSerial != 0 &&
	allocationCertificate.allocationPlanSerial ==
	    state->activeCadAllocationPlan() &&
	allocationCertificate.viewRevision ==
	    this->d->lodViewRevision.value() &&
	allocationCertificate.policyRevision ==
	    this->d->lodPolicyRevision.value() &&
	std::isfinite(allocationCertificate.pointProxyPixelThreshold) &&
	std::fabs(allocationCertificate.pointProxyPixelThreshold -
	    this->d->lodPresentationPointProxyPixelThreshold) <= 1.0e-6f &&
	state->cadAllocationPlanCoversCurrentPopulation(
	    allocationCertificate.allocationPlanSerial,
	    allocationCertificate.viewRevision,
	    allocationCertificate.policyRevision,
	    allocationCertificate.fixedCadPresentationCost);
    int interruptedCorrectedCeiling = interruptedPresentedMaximum > 0 ?
	interruptedPresentedMaximum - 1 : interruptedPresentedMaximum;
    size_t interruptedEstimatedCost = 0;
    size_t interruptedTargetCost = 0;
    /*
     * A PoP cut ordinal is a precision boundary, not a cost ratio.  Once the
     * cache supports dozens of cuts, subtracting one ordinal per missed frame
     * can spend seconds walking past cuts which remove almost no work.  The
     * view state maintains the exact retained cost at every renderer-wide
     * ceiling.  For a quiet draw-capacity miss, scale the attempted cost by
     * the observed deadline ratio and select the richest prefix which meets
     * that corrected cost in one step.  Preparation-only retries remain
     * unchanged below, and the hard deadline still corrects an optimistic
     * prediction.
     *
     * This aggregate includes off-frustum retained occurrences and excludes
     * renderer point collapse, so it is deliberately conservative.  It does
     * not change occurrence cuts or resident data; later completed-frame
     * headroom policy may recover fidelity without cache I/O.
     */
    if (state && !interactive && !cadPreparationChanged &&
	interruptedPresentedMaximum > 0 && elapsedNanoseconds > 0 &&
	this->d->stablePresentationFrameDeadlineNanoseconds > 0) {
	if (certifiedUnconstrainedPresentation) {
	    interruptedEstimatedCost =
		allocationCertificate.selectedPresentationCost;
	} else {
	    const size_t activeCadCost = state->activeCadRenderCost();
	    const size_t nonCadCost = interruptedActiveCost > activeCadCost ?
		interruptedActiveCost - activeCadCost : 0;
	    const size_t presentedCadEstimate =
		state->cadRenderCostAtProgressiveCutCeiling(
		    interruptedPresentedMaximum);
	    interruptedEstimatedCost = nonCadCost >
		    SIZE_MAX - presentedCadEstimate ? SIZE_MAX :
		nonCadCost + presentedCadEstimate;
	}
	const long double deadlineRatio = std::min<long double>(
	    1.0L,
	    static_cast<long double>(
		this->d->stablePresentationFrameDeadlineNanoseconds) /
	    static_cast<long double>(elapsedNanoseconds));
	const long double corrected =
	    static_cast<long double>(interruptedEstimatedCost) *
	    deadlineRatio * 0.80L;
	interruptedTargetCost = corrected >=
		static_cast<long double>(SIZE_MAX) ? SIZE_MAX :
	    static_cast<size_t>(std::max<long double>(0.0L, corrected));
	const size_t sceneBudget = this->d->lodBudgetPolicy.currentBudget();
	if (sceneBudget != SIZE_MAX)
	    interruptedTargetCost = std::min(
		interruptedTargetCost, sceneBudget);
	const size_t activeCadCost = state->activeCadRenderCost();
	const size_t nonCadCost = interruptedActiveCost > activeCadCost ?
	    interruptedActiveCost - activeCadCost : 0;
	const size_t cadTarget = interruptedTargetCost > nonCadCost ?
	    interruptedTargetCost - nonCadCost : 0;
	const int predictedCeiling =
	    state->cadProgressiveCutWithinRenderCost(
		cadTarget, interruptedPresentedMaximum - 1);
	interruptedCorrectedCeiling = predictedCeiling >= 0 ?
	    predictedCeiling : 0;
	if (this->d->lodDeadlineSafeProgressiveCeiling >= 0 &&
	    this->d->lodDeadlineSafeViewRevision ==
		this->d->lodViewRevision.value() &&
	    this->d->lodDeadlineSafeQualityDomainRevision ==
		this->d->lodViewQualityDomainRevision)
	    /* The prior completed frame proves this cut is affordable for the
	     * unchanged camera and source population.  Preserve that visual floor
	     * when a richer speculative cut misses; the enclosing clamp still
	     * forces a retry below the cut which actually missed. */
	    interruptedCorrectedCeiling = std::max(
		interruptedCorrectedCeiling,
		this->d->lodDeadlineSafeProgressiveCeiling);
    }
    /* Any abort which advanced retained preparation is forward progress, not
	 * draw-capacity evidence.  Large command records are deliberately built in
	 * resumable deadline slices and may need more than one slice at 150k+
	 * occurrences.  Coarsening the second slice invalidates the work just
	 * prepared and creates a normalize/refine loop.  A genuinely stuck retry
	 * changes no preparation token and reaches the ordinary pressure path. */
    const SbBool transientPresentationRetry =
	BObolLodQualityPolicy::retryTransientPresentation(
	    interactive != FALSE,
	    this->d->consecutiveInterruptedPresentationFrames,
	    cadPreparationChanged != FALSE,
	this->d->lodPublicationPolicy.framePending(),
	this->d->lodFrameObligation.pending(),
	this->d->lodStablePointProxyCalibrationPending != FALSE,
	this->d->lodStaticOverscanActive != FALSE) ?
	TRUE : FALSE;
    /* A preparation/publication replay has not tested renderer capacity and
     * must retain its pending witness.  A steady quiet draw which crosses the
     * hard endpoint deadline has tested it: record a negative proof instead
     * of merely cancelling the pending retry.  Cancelling forgot the failed
     * population, so the cheap constrained successor immediately re-armed the
     * same rich frame and made large scenes visibly alternate coarse/fine.
     * Static quality is likewise one-shot for this camera/policy epoch after
     * such a miss; occurrence-local recovery remains free to redistribute
     * the demonstrated safe allowance. */
    const bool deadlinePopulationWorkPending =
	this->d->progressiveProviderPendingCount > 0 ||
	this->d->lodSubmissionPending ||
	this->d->lodResultsPending.load() != 0 ||
	this->d->lodPublicationPolicy.pending() ||
	this->d->lodResidentGrowthPolicy.pending() ||
	this->d->lodResidentGrowthResidencyDrainActive ||
	(this->d->lodService &&
	 (this->d->lodService->activeTaskCountForGeneration(
	      this->d->lodActiveGeneration) > 0 ||
	  this->d->lodService->queuedResultCountForGeneration(
	      this->d->lodActiveGeneration) > 0));
    if (!interactive && cadDrawAttempted && !transientPresentationRetry) {
	const size_t failedCapacityBudget = interruptedEstimatedCost > 1 ?
	    interruptedEstimatedCost :
	    this->d->lodBudgetPolicy.currentBudget();
	this->d->lodBudgetPolicy.noteDeadlineCapacityMiss(
	    failedCapacityBudget);
	this->d->lodHeadroomPolicy.rejectRetry(
	    this->d->lodViewRevision, this->d->lodPolicyRevision,
	    interruptedActiveCost, this->d->lodBudgetPolicy.currentBudget());
	if (BObolLodStaticQualityPolicy::terminalCapacityMiss(
		interactive != FALSE, cadDrawAttempted != FALSE,
		transientPresentationRetry != FALSE,
		this->d->lodStaticOverscanActive != FALSE,
		deadlinePopulationWorkPending))
	    this->d->lodStaticOverscanRejected = TRUE;
    } else if (interactive) {
	this->d->lodHeadroomPolicy.cancelRetry();
    }
    if (getenv("BOBOL_LOD_TRACE_DEADLINE"))
	bu_log("BObol LoD render deadline elapsed_ms=%.3f "
	       "draw=%d preparation=%d interactive=%d "
	       "active_max=%d presented_max=%d ceiling=%d "
	       "static_active=%d static_rejected=%d headroom=%d "
	       "floor_budget=%zu floor_misses=%u floor_rejected=%d "
	       "capacity_ceiling=%zu "
	       "certified=%d "
	       "active_cost=%zu minimum_cost=%zu handoff=%d "
	       "estimated_cost=%zu target_cost=%zu corrected=%d\n",
	       elapsedNanoseconds / 1000000.0,
	       cadDrawAttempted ? 1 : 0,
	       cadPreparationChanged ? 1 : 0,
	       interactive ? 1 : 0,
	       interruptedActiveMaximum,
	       interruptedPresentedMaximum,
	       this->d->lodInteractiveProgressiveCeiling,
	       this->d->lodStaticOverscanActive ? 1 : 0,
	       this->d->lodStaticOverscanRejected ? 1 : 0,
	       this->d->lodHeadroomPolicy.retryPending() ? 1 : 0,
	       this->d->lodBudgetPolicy.retainedQualityFloorBudget(),
	       this->d->lodBudgetPolicy.retainedQualityFloorMissCount(),
	       this->d->lodBudgetPolicy.retainedQualityFloorRejected() ? 1 : 0,
	       this->d->lodBudgetPolicy.deadlineCapacityCeiling(),
	       certifiedUnconstrainedPresentation ? 1 : 0,
	       state ? state->activeRenderCost() : 0,
	       state ? state->minimumActiveRenderCost() : 0,
	       this->d->lodPresentationPolicy.handoffPending() ? 1 : 0,
	       interruptedEstimatedCost, interruptedTargetCost,
	       interruptedCorrectedCeiling);
    if (!interactive && cadDrawAttempted && !transientPresentationRetry &&
	this->d->lodInteractiveProgressiveCeiling < 0 &&
	this->d->lodBudgetPolicy.retainedQualityFloorActive()) {
	const size_t floorBudget =
	    this->d->lodBudgetPolicy.retainedQualityFloorBudget();
	const uint64_t floorSignature =
	    this->d->lodBudgetPolicy.retainedQualityFloorSignature();
	const unsigned int missesBefore =
	    this->d->lodBudgetPolicy.retainedQualityFloorMissCount();
	const bool rejected =
	    this->d->lodBudgetPolicy.noteRetainedQualityFloorMiss();
	if (getenv("BOBOL_LOD_TRACE_BUDGET"))
	    bu_log("BObol LoD quality-floor miss budget=%zu "
		   "signature=%llu miss=%u elapsed_ms=%.3f "
		   "active_cost=%zu active_max=%d rejected=%d\n",
		   floorBudget, (unsigned long long)floorSignature,
		   missesBefore + 1, elapsedNanoseconds / 1000000.0,
		   state ? state->activeRenderCost() : 0,
		   interruptedActiveMaximum, rejected ? 1 : 0);
    }

    /* A Coin abort bounds total endpoint-thread traversal, not a completed
     * geometry presentation.  In particular, a newly rotated many-part CAD
     * scene may spend the whole deadline rebuilding command records before
     * the GPU cut can be measured.  activeRenderCost() is retained occurrence
     * demand (often richer than the renderer ceiling), so pairing it with
     * this elapsed time produced a fictitious low triangle rate and
     * destructively compacted otherwise proven resident PoP prefixes.
     *
     * Keep deadline recovery presentation-only below: lower the reversible
     * renderer ceiling and/or aggregate sub-pixel occurrences.  Only a
     * completed frame publishes exact presented work and may update the
     * persistent scene budget or throughput calibration. */

    /* A deadline abort is already sufficient evidence that the submitted
     * progressive cut is too expensive; waiting for a completed frame to
     * correct it creates a retry deadlock.  This applies both to a scale
     * interaction and to a quiet admission whose first coherent frame
     * exceeded the stable deadline.  Lower only the renderer-wide prefix
	 * ceiling.  Quiet frames use the retained per-cut cost aggregate above;
	 * interactive frames keep the conservative one-cut response.  The
	 * occurrence cut and resident suffix remain intact, so this is immediate,
	 * reversible, and cannot restart level walking.  A later bounded stable
	 * pass reconciles the retained cut before removing the ceiling.
     *
     * This applies to pose-only motion as well as zoom.  A prepared pose frame
     * normally lets completed-frame feedback choose the ceiling, but a hard
     * abort supplies no such sample; leaving pose motion excluded retries the
     * same over-budget cut until button-up.  Lowering only this global ceiling
     * keeps every occurrence cut and resident suffix available for immediate
     * quiet restoration. */
    if (state && cadDrawAttempted && !transientPresentationRetry) {
	BObolViewLodState *presentationState =
	    this->d->viewAttachment->getViewLodState();
	if (interruptedPresentedMaximum > 0) {
	    const int correctedCeiling = std::max(0,
		std::min(interruptedPresentedMaximum - 1,
		    interruptedCorrectedCeiling));
	    this->d->lodViewDemandPolicy.noteQualityMiss(
		correctedCeiling, state->activeRenderCost());
	    this->d->lodInteractiveProgressiveCeiling = correctedCeiling;
	    presentationState->setCadPresentationProgressiveCutCeiling(
		correctedCeiling);
	    if (!interactive) {
		/* Preserve the safe retained-population cost computed from the
		 * interrupted full frame.  Timing the newly installed coarse renderer
		 * ceiling is not evidence that the hidden full population fits; using
		 * that fallback frame to enlarge this budget recreates the same failed
		 * ceiling-release attempt forever. */
		this->d->lodPresentationPolicy.armHandoff(
		    true, 0, interruptedTargetCost);
		/* The camera and policy epochs did not change.  Start an explicit
		 * same-epoch pass; clearing the epoch witness while also setting a
		 * pending cursor makes the wrapper misclassify this as a view change
		 * during submission and append another complete rescan. */
		this->d->lodSubmissionSourceIndex = 0;
		this->d->lodSubmissionEntryOffset = 0;
		this->d->clearLodSubmissionPlan();
		this->d->lodSubmissionRescanPending = FALSE;
		this->d->lodSubmissionPending = TRUE;
	    }
	}
    }

    /* A static overscan miss is terminal evidence for this view/capacity
     * epoch.  Keep the richer occurrence prefix resident and the corrected
     * renderer ceiling above; do not let a later repaint reopen the same
     * failed quality staircase.  Keep the static allowance active only while
     * its handoff reconciles the retained occurrence cut.  Dropping to the
     * ordinary quiet allowance here discarded the last deadline-safe static
     * frame before that reconciliation could use its proof. */
    const bool rejectedStaticQualityTrial =
	BObolLodStaticQualityPolicy::terminalCapacityMiss(
	    interactive != FALSE, cadDrawAttempted != FALSE,
	    transientPresentationRetry != FALSE,
	    this->d->lodStaticOverscanActive != FALSE,
	    deadlinePopulationWorkPending);
    if (rejectedStaticQualityTrial) {
	this->d->lodStaticOverscanLeapAvailable = FALSE;
	this->d->lodStaticOverscanRejected = TRUE;
	/* The protected visual-significance floor is permitted to exceed the
	 * preferred quiet cadence only under this exact hard-deadline trial.  Once
	 * that trial fails, retaining the floor makes the recovery allocator
	 * restore the same failed population after every constrained frame.  Keep
	 * all resident PoP suffixes, but make the cheaper allocation a terminal
	 * fixed point for this view/capacity epoch. */
	const bool rejectedQualityFloor =
	    this->d->lodBudgetPolicy.rejectRetainedQualityFloor();
	if (rejectedQualityFloor && getenv("BOBOL_LOD_TRACE_BUDGET"))
	    bu_log("BObol LoD quality-floor rejected by static deadline\n");
    }

    /* A many-part view can hit a harder floor than PoP triangle detail: one
     * coherent minimum prefix per visible occurrence.  The normal small-part
     * aggregation controller learns only from completed frames, which is an
     * impossible prerequisite when every attempted minimum-prefix frame is
     * itself interrupted.  Treat a quiet deadline abort as conservative
     * feedback for the reversible presentation-only point cut as well.  This
     * keeps all immutable meshes and desired cuts resident, aggregates only
     * occurrences below the camera-local screen threshold, and lets the
     * ordinary completed-frame relaxation probe recover the finest sustainable
     * threshold afterward.
     *
     * Do not apply this to active input: pose/scale interaction already owns
     * a separately measured cut and repeated unapplied motion feedback could
     * otherwise over-coarsen the first post-input frame. */
    /* Unlike the renderer-wide PoP ceiling above, a point threshold changes
     * the occurrence population itself.  No preparation-heavy abort is valid
	 * evidence for that bracket: a lower PoP ceiling may make retained
	 * construction finish, while poisoning the
     * point bracket would preserve visibly coarse multi-pixel objects after
     * the one-time work is gone. */
    /* activeRenderCost() describes the retained occurrence cuts, while the
     * renderer-wide ceiling describes what this interrupted frame actually
     * attempted.  Once that reversible ceiling reaches the minimum PoP cut,
     * triangle population presented by every progressive occurrence is
     * already irreducible even if rich resident suffixes keep active cost
     * above minimumActiveRenderCost().  Requiring retained cost itself to be
     * minimal deadlocks a many-part software view at the minimum cut: no further
     * triangle correction is possible, yet the small-part aggregation
     * controller is never allowed to act. */
    const SbBool interruptedPopulationIrreducible =
	BObolLodPointProxyCalibrationPolicy::
	    deadlineRequiresPopulationAggregation(
		interruptedActiveCost, interruptedMinimumCost,
		interruptedPresentedMaximum, interruptedCorrectedCeiling,
		interruptedTargetCost) ? TRUE : FALSE;
    if (state && state->hasCadPresentationAssemblies() &&
	cadDrawAttempted && !cadPreparationChanged && !interactive &&
	interruptedPopulationIrreducible &&
	this->d->pointProxyAggregationApplicable() &&
	this->d->quietAllocationTargetFps() > 0.0f) {
	const BObolLodPointProxyCalibrationPolicy::Decision pressure =
	    this->d->lodPointProxyCalibrationPolicy.interrupted(
		this->d->lodPresentationPointProxyPixelThreshold,
		elapsedNanoseconds, this->d->quietAllocationTargetFps());
	if (pressure.changed) {
	    /* A coarser point cut changes only the presentation used to recover
	     * from the failed one-pixel trial.  It does not make that rejected
	     * one-pixel population cheaper.  Keep the rejection authoritative for
	     * this view/capacity epoch while the safe/unsafe point bracket settles;
	     * reopening the trial here produced an exact 1 -> coarse -> 1 cycle on
	     * 50k OSMesa wire views. */
	    this->d->lodPresentationPointProxyPixelThreshold =
		pressure.threshold;
	    BObolViewLodState *presentationState =
		this->d->viewAttachment->getViewLodState();
	    presentationState->setCadPresentationPointProxyPixelThreshold(
		pressure.threshold);
	    this->d->lodStablePointProxyCalibrationPending = TRUE;
	}
    }

    if (transientPresentationRetry) {
	/* Preserve every admission and calibration barrier until the first
	 * changed population has one unchanged replay.  A second non-preparation
	 * abort takes the ordinary bounded capacity-recovery path above. */
	this->requestRender(cadPreparationChanged ?
	    "render-preparation-replay" : "render-population-replay");
	return;
    }

    /* An interrupted traversal requires another bounded policy/presentation
     * attempt, but it is not a persistent capacity sample.  Keep the barriers
     * armed while the renderer-only corrections above make that retry
     * cheaper; retained-cut admission remains based on completed frames. */
    this->d->lodBudgetPolicy.resetOverloadRecovery();
    this->d->lodBudgetPolicy.setProbeCount(3);
    this->d->lodBudgetPolicy.resetPass();
    this->markProgressiveWorkPending();
    this->requestRender("render-deadline");
}

uint64_t
BObolViewController::getInterruptedPresentationFrameCount(void) const
{
    return this->d->interruptedPresentationFrameCount;
}

uint64_t
BObolViewController::getLastInterruptedPresentationTimeNanoseconds(void) const
{
    return this->d->lastInterruptedPresentationTimeNanoseconds;
}

void
BObolViewController::armStableLodHeadroomProbeIfReady(void)
{
    /* A threshold mutation is not authoritative until the renderer has
     * classified one exact frame at that threshold.  Provider pumping may
     * continue behind this barrier, but neither another threshold change nor
     * mesh admission may consume the preceding classifier result.  This
     * includes a one-pixel trial armed while completing a coarse structural
     * seed: the histogram below still describes the coarse frame and must not
     * overwrite the successor trial before it reaches the renderer. */
    if (!BObolLodPointProxyCalibrationPolicy::
	maySeedStructuralDistribution(
	    this->d->lodAdmissionPointProxyFramePending != FALSE,
	    this->d->lodStablePointProxyCalibrationPending != FALSE,
	    this->d->lodPointProxyTriangleRecoveryPending != FALSE))
	return;

    const BObolViewLodState *viewLodState = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    const size_t activePopulationCost = viewLodState ?
	viewLodState->activeRenderCost() : 0;
    const size_t activeTasks = this->d->lodService ?
	this->d->lodService->activeTaskCountForGeneration(
	    this->d->lodActiveGeneration) : 0;
    const size_t queuedResults = this->d->lodService ?
	this->d->lodService->queuedResultCountForGeneration(
	    this->d->lodActiveGeneration) : 0;
    size_t activePayloads = 0;
    size_t satisfiedPayloads = 0;
    size_t memoryLimitedPayloads = 0;
    if (viewLodState)
	viewLodState->convergencePayloadCounts(activePayloads,
	    satisfiedPayloads, memoryLimitedPayloads);

    size_t presentedSubpixelOccurrences = 0;
    size_t presentedStructuralBoxes = 0;
    const bool exactOccurrenceCoverage = viewLodState &&
	viewLodState->lastCadPresentationOccurrenceCoverage(
	    presentedSubpixelOccurrences, presentedStructuralBoxes);
    const size_t terminalOccurrenceFailures = viewLodState ?
	viewLodState->cadOccurrenceTerminalFailureCount() : 0;
    const bool externalProducersSettled =
	this->d->progressiveProviderPendingCount == 0;

    /* Initial structural coverage and terminal presentation coverage are
     * deliberately different proofs.  The former gets a complete model on
     * screen quickly; the latter permits only meshes, camera-valid subpixel
     * points, or an explicit provider failure.  A bounded source scan can
     * finish its box-first pass and later exhaust a refinement window without
     * visiting every remaining box.  If an exact quiet frame observes that
     * state, resume a normal (not structural-only) current-view pass from the
     * retained population.  No source data or PoP prefix is discarded.
     *
     * This is also the resize recovery edge: a resize may change which boxes
     * are physically subpixel while preserving the view/source epochs.  The
     * completed framebuffer is the authoritative classification, so it must
     * be able to re-arm admission without waiting for another mouse event. */
    const bool unresolvedStructuralPresentation =
	exactOccurrenceCoverage &&
	presentedStructuralBoxes > terminalOccurrenceFailures;
    /* The renderer has already projected every structural fallback in order
     * to draw this frame.  Reuse its exact cumulative size distribution to
     * bound first-mesh provider work before a cold/warm scene loads thousands
     * of tiny meshes which its eventual point presentation will not submit.
     * Obol owns view classification; libBObol owns this admission decision. */
    BObolViewLodState::CadStructuralProjectionHistogram structuralProjection;
    const bool exactStructuralProjection = viewLodState &&
	viewLodState->lastCadStructuralProjectionHistogram(
	    structuralProjection);
    const size_t sceneBudget = this->d->lodBudgetPolicy.currentBudget();
    const size_t firstWaveOccurrenceFloor = 64;
    const size_t maximumFirstWaveOccurrences = sceneBudget == SIZE_MAX ?
	SIZE_MAX : std::max<size_t>(512, std::min<size_t>(8192,
	    sceneBudget / firstWaveOccurrenceFloor));
    const bool haveStructuralProjectionPopulation =
	exactStructuralProjection &&
	structuralProjection.visibleCount > terminalOccurrenceFailures;

    /* A framebuffer histogram is exact for the occurrences installed when
     * that frame began, but it is not a whole-scene population proof while a
     * provider is still appending leaves.  Reseeding from each partial
     * population restarts the source cursor and turns parallel discovery into
     * a serial threshold/frame staircase.  The structural proxies already
     * provide immediate useful coverage during discovery; seed the first
     * mesh wave once, from the settled inventory which it is meant to bound. */
    if (haveStructuralProjectionPopulation &&
	externalProducersSettled &&
	!this->d->lodInteractive && !this->d->lodGestureActive &&
	BObolLodPointProxyCalibrationPolicy::applicable(
	    structuralProjection.visibleCount) &&
	maximumFirstWaveOccurrences != SIZE_MAX) {
	const BObolLodPointProxyCalibrationPolicy::Decision seed =
	    this->d->lodPointProxyCalibrationPolicy.
		seedFromStructuralDistribution(
		    this->d->lodPresentationPointProxyPixelThreshold,
		    structuralProjection.cumulativeCount,
		    structuralProjection.visibleCount,
		    maximumFirstWaveOccurrences);
	if (seed.changed) {
	    const float oldEffectivePointThreshold = std::max(
		this->d->lodPresentationPointProxyPixelThreshold,
		this->d->lodDiscoveryPointProxyPixelThreshold);
	    this->d->lodPresentationPointProxyPixelThreshold = seed.threshold;
	    this->d->lodDiscoveryPointProxyPixelThreshold = 1.0f;
	    this->d->lodDiscoveryPointProxyPolicy.reset();
	    const float newEffectivePointThreshold = seed.threshold;
	    const bool classifierChanged = std::fabs(
		newEffectivePointThreshold - oldEffectivePointThreshold) >
		1.0e-6f;
	    this->d->lodAdmissionPointProxyFramePending = classifierChanged ?
		TRUE : FALSE;
	    BObolViewLodState *presentationState =
		this->d->viewAttachment->getViewLodState();
	    if (presentationState) {
		presentationState->setCadPresentationPointProxyPixelThreshold(
		    seed.threshold);
		presentationState->
		    setCadPresentationDiscoveryPointProxyPixelThreshold(1.0f);
	    }
	    /* The old cursor was admitted against a different point/mesh split.
	     * Keep resident results, but restart planning from the new exact
	     * presentation after its bounded classifier frame completes. */
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionDeltaActive = FALSE;
	    this->d->lodSubmissionDeltaSources.clear();
	    this->d->lodSubmissionDeltaPlans.clear();
	    this->d->lodBudgetPolicy.resetPass();
	    /* This frame classifies which structural occurrences may enter the
	     * first mesh wave.  It is an admission barrier, not a stable timing
	     * calibration: the scene owns no mesh population to time yet.  Once
	     * completeRenderTiming() consumes the exact classifier frame,
	     * armStableLodHeadroomProbeIfReady() routes its remaining box frontier
	     * into the importance-ordered repair pass. */
	    this->d->lodStablePointProxyCalibrationPending = FALSE;
	    if (classifierChanged) {
		this->requestRender("lod-structural-distribution-seed");
		this->d->reconcilePhase(
		    BObolLodStateMachine::Event::WORK_SCHEDULED);
		return;
	    }
	    /* Discovery already presented this effective threshold.  Moving it
	     * into the stable-policy bracket is an ownership transfer, not a
	     * framebuffer mutation; fall through and admit the exact structural
	     * frontier without spending or waiting for a duplicate frame. */
	}
    }
    if (externalProducersSettled &&
	this->d->lodDiscoveryPointProxyPixelThreshold > 1.01f) {
	const float oldEffective = std::max(
	    this->d->lodPresentationPointProxyPixelThreshold,
	    this->d->lodDiscoveryPointProxyPixelThreshold);
	this->d->lodDiscoveryPointProxyPixelThreshold = 1.0f;
	this->d->lodDiscoveryPointProxyPolicy.reset();
	BObolViewLodState *presentationState =
	    this->d->viewAttachment->getViewLodState();
	if (presentationState)
	    presentationState->
		setCadPresentationDiscoveryPointProxyPixelThreshold(1.0f);
	const float newEffective =
	    this->d->lodPresentationPointProxyPixelThreshold;
	if (newEffective + 1.0e-6f < oldEffective) {
	    this->d->lodAdmissionPointProxyFramePending = TRUE;
	    this->requestRender("lod-discovery-point-release");
	    this->d->reconcilePhase(
		BObolLodStateMachine::Event::WORK_SCHEDULED);
	    return;
	}
    }
    const bool structuralRepairReady = unresolvedStructuralPresentation &&
	this->d->lodAutoSubmit && this->d->lodService &&
	!this->d->lodInteractive && !this->d->lodGestureActive &&
	this->d->lodCoveragePolicy.coverageComplete() &&
	!this->d->lodSubmissionPending &&
	!this->d->lodBudgetPolicy.rescanAfterFrame() &&
	!this->d->lodFrameObligation.pending() &&
	!this->d->lodRetainedRefinementPending &&
	!this->d->lodStablePointProxyCalibrationPending &&
	!this->d->lodPointProxyTriangleRecoveryPending &&
	externalProducersSettled &&
	this->d->lodResultsPending.load() == 0 && activeTasks == 0 &&
	queuedResults == 0;
    /* A presentation handoff cannot reconcile its renderer-wide ceiling while
     * visible structural boxes are absent from the retained occurrence
     * population.  Structural repair therefore precedes, rather than waits
     * for, handoff completion.  The existing ceiling remains installed while
     * the closed repair transaction runs, so continuity and its hard frame
     * deadline are preserved.  Making both phases wait on each other forced a
     * warm 150k view to admit the residual boxes through one ordinary provider
     * wave at a time before either phase could finish. */
    if (getenv("BOBOL_LOD_TRACE_STRUCTURAL_REPAIR") &&
	(exactOccurrenceCoverage || presentedStructuralBoxes))
	bu_log("BObol LoD structural repair unresolved=%d ready=%d "
	       "exact=%d boxes=%zu failures=%zu interactive=%d gesture=%d "
	       "coverage=%d submit=%d rescan=%d refinement_frame=%d "
	       "retained=%d handoff=%d point=%d recovery=%d results=%d "
	       "active_tasks=%zu queued_results=%zu\n",
	       unresolvedStructuralPresentation ? 1 : 0,
	       structuralRepairReady ? 1 : 0,
	       exactOccurrenceCoverage ? 1 : 0,
	       presentedStructuralBoxes, terminalOccurrenceFailures,
	       this->d->lodInteractive ? 1 : 0,
	       this->d->lodGestureActive ? 1 : 0,
	       this->d->lodCoveragePolicy.coverageComplete() ? 1 : 0,
	       this->d->lodSubmissionPending ? 1 : 0,
	       this->d->lodBudgetPolicy.rescanAfterFrame() ? 1 : 0,
	       this->d->lodFrameObligation.pending() ? 1 : 0,
	       this->d->lodRetainedRefinementPending ? 1 : 0,
	       this->d->lodPresentationPolicy.handoffPending() ? 1 : 0,
	       this->d->lodStablePointProxyCalibrationPending ? 1 : 0,
	       this->d->lodPointProxyTriangleRecoveryPending ? 1 : 0,
	       this->d->lodResultsPending.load(), activeTasks, queuedResults);

    /* Once the point classifier reaches its 64-pixel bound, every remaining
     * uncollapsed box is an irreducible coverage obligation.  Transfer the
     * immediately preceding exact frame's hard-deadline headroom into one
     * structural admission transaction.  The policy rejects an unchanged
     * population after that attempt, preventing the preferred 20 Hz budget
     * from replaying the same sparse frontier indefinitely. */
    const size_t unresolvedStructuralCount =
	presentedStructuralBoxes > terminalOccurrenceFailures ?
	    presentedStructuralBoxes - terminalOccurrenceFailures : 0;
    const size_t unaggregatableStructuralCount = exactStructuralProjection ?
	BObolLodStructuralAdmissionPolicy::unaggregatableCount(
	    structuralProjection.cumulativeCount,
	    structuralProjection.visibleCount) : 0;
    size_t structuralPresentedCost = 0;
    const bool exactStructuralCapacityFrame =
	exactOccurrenceCoverage && viewLodState &&
	viewLodState->lastCadPresentationFrameExact() &&
	viewLodState->lastCadPresentedRenderCost(structuralPresentedCost) &&
	structuralPresentedCost > 0 &&
	this->d->lastRenderTimeNanoseconds > 0 &&
	this->d->stablePresentationFrameDeadlineNanoseconds > 0;
    size_t structuralCertifiedBudget = 0;
    if (exactStructuralCapacityFrame) {
	const size_t presentationLimit =
	    BObolLodQualityPolicy::staticPresentationRenderCostLimit(
		structuralPresentedCost,
		this->d->lastRenderTimeNanoseconds,
		this->d->stablePresentationFrameDeadlineNanoseconds);
	const size_t additionalCapacity =
	    presentationLimit > structuralPresentedCost ?
		presentationLimit - structuralPresentedCost : 0;
	structuralCertifiedBudget =
	    additionalCapacity > SIZE_MAX - activePopulationCost ? SIZE_MAX :
		activePopulationCost + additionalCapacity;
	structuralCertifiedBudget = std::max(
	    structuralCertifiedBudget,
	    this->d->lodBudgetPolicy.currentBudget());
    }
    BObolLodStructuralAdmissionPolicy::Decision structuralAdmission;
    if (structuralRepairReady && exactStructuralCapacityFrame) {
	this->d->lodStructuralCoverageCostReservation = 0;
	BObolLodStructuralAdmissionPolicy::Inputs inputs;
	inputs.viewRevision = this->d->lodViewRevision.value();
	inputs.policyRevision = this->d->lodPolicyRevision.value();
	inputs.projectionRevision = structuralProjection.revision;
	inputs.unresolvedCount = unresolvedStructuralCount;
	inputs.unaggregatableCount = unaggregatableStructuralCount;
	inputs.activeCost = activePopulationCost;
	inputs.currentBudget = this->d->lodBudgetPolicy.currentBudget();
	inputs.certifiedBudget = structuralCertifiedBudget;
	inputs.exactProjection = exactStructuralProjection;
	inputs.maximumThreshold =
	    this->d->lodPresentationPointProxyPixelThreshold >= 63.99f;
	structuralAdmission =
	    this->d->lodStructuralAdmissionPolicy.evaluate(inputs);
	if (structuralAdmission.ownsFrontier) {
	    size_t admittedBudget = 0;
	    if (structuralAdmission.requestAdmission)
		admittedBudget = this->d->lodBudgetPolicy.
		    requestCoverageCompletion(activePopulationCost,
			structuralAdmission.budget);
	    const size_t provisionalReservation =
		BObolLodStructuralAdmissionPolicy::perOccurrenceReservation(
		    admittedBudget, activePopulationCost,
		    unresolvedStructuralCount);
	    if (getenv("BOBOL_LOD_TRACE_STRUCTURAL_BUDGET"))
		bu_log("BObol LoD structural capacity boxes=%zu "
		       "unaggregatable=%zu active=%zu current=%zu "
		       "certified=%zu admitted=%zu reservation=%zu "
		       "duplicate=%d limited=%d\n",
		       unresolvedStructuralCount,
		       unaggregatableStructuralCount, activePopulationCost,
		       inputs.currentBudget, structuralCertifiedBudget,
		       admittedBudget,
		       provisionalReservation,
		       structuralAdmission.duplicateRejected ? 1 : 0,
		       structuralAdmission.capacityLimited ? 1 : 0);
	    if (!structuralAdmission.requestAdmission || !admittedBudget ||
		!provisionalReservation) {
		if (this->d->lastLodDiagnostics.getLength() > 0)
		    this->d->lastLodDiagnostics += "\n";
		this->d->lastLodDiagnostics += "structural-admission: ";
		this->d->lastLodDiagnostics +=
		    structuralAdmission.duplicateRejected ?
			"scene LoD rejected an unchanged structural admission frontier" :
			"scene LoD structural coverage reached static-frame capacity";
		return;
	    }
	}
    }
    if (structuralRepairReady) {
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodSubmissionDeltaSources.clear();
	this->d->lodSubmissionDeltaPlans.clear();
	/* The completed framebuffer already performed the exact camera/frustum
	 * classification.  Route only its visible, uncollapsed structural
	 * occurrences back to the source loader.  A full compact scan here loads
	 * thousands of meshes whose subpixel point representations already satisfy
	 * the view and can turn a 293-box repair into 47k unnecessary payloads.
	 * Fall back to the complete scan if any retained identity cannot be mapped;
	 * incomplete selective repair would violate the no-box terminal contract. */
	size_t mappedStructuralEntries = 0;
	const std::vector<SoBRLDatabaseSource *> repairSources =
	    controller_render_database_source_roots(this);
	for (SoBRLDatabaseSource *source : repairSources) {
	    if (!source || !source->hasCompactInstanceIndex())
		continue;
	    SoBRLCadAssembly *assembly = dynamic_cast<SoBRLCadAssembly *>(
		viewLodState->findCadPresentation(source));
	    if (!assembly)
		continue;
	    std::vector<SbString> occurrenceKeys;
	    assembly->lastUncollapsedStructuralProxyOccurrenceKeys(
		occurrenceKeys);
	    std::vector<size_t> entryIndices;
	    entryIndices.reserve(occurrenceKeys.size());
	    for (const SbString &occurrenceKey : occurrenceKeys) {
		size_t entryIndex = 0;
		if (source->getCompactInstanceIndex(
			occurrenceKey.getString(), entryIndex))
		    entryIndices.push_back(entryIndex);
	    }
	    std::sort(entryIndices.begin(), entryIndices.end());
	    entryIndices.erase(
		std::unique(entryIndices.begin(), entryIndices.end()),
		entryIndices.end());
	    if (entryIndices.empty())
		continue;
	    mappedStructuralEntries += entryIndices.size();
	    this->d->lodSubmissionDeltaSources.push_back(source);
	    this->d->lodSubmissionDeltaPlans.emplace_back(
		source, std::move(entryIndices));
	}
	const bool exactStructuralFrontier =
	    mappedStructuralEntries == presentedStructuralBoxes &&
	    mappedStructuralEntries > 0;
	if (getenv("BOBOL_LOD_TRACE_STRUCTURAL_REPAIR"))
	    bu_log("BObol LoD structural repair frontier boxes=%zu "
		   "mapped=%zu exact=%d sources=%zu\n",
		   presentedStructuralBoxes, mappedStructuralEntries,
		   exactStructuralFrontier ? 1 : 0,
		   this->d->lodSubmissionDeltaPlans.size());
	this->d->lodSubmissionDeltaActive =
	    exactStructuralFrontier ? TRUE : FALSE;
	if (!exactStructuralFrontier) {
	    this->d->lodSubmissionDeltaSources.clear();
	    this->d->lodSubmissionDeltaPlans.clear();
	}
	this->d->lodSubmissionRefreshMissing = TRUE;
	this->d->lodSubmissionReset = 0;
	this->d->lodStructuralPresentationRepairPending = TRUE;
	this->d->lodStructuralRepairFrontierCount =
	    presentedStructuralBoxes - terminalOccurrenceFailures;
	this->d->lodStructuralCoverageCostReservation = 0;
	/* This accumulator belongs to the newly closed repair transaction.  A
	 * prior coverage/demand pass may have exhausted its budget while the
	 * renderer was still classifying boxes; carrying that debt into this exact
	 * frontier spuriously re-entered point aggregation after every box had
	 * already been replaced. */
	this->d->lodMissingMeshBudgetBlockedCount = 0;
	this->d->lodSubmissionPending = TRUE;
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodBudgetPolicy.clearBudgetLimit();
	this->d->lodBudgetPolicy.resetPass();
	this->markProgressiveWorkPending();
	this->requestRender("lod-structural-presentation-repair");
	return;
    }

    /* A complete one-pixel image is the first stable result.  When that exact
     * retained population uses less than one quarter of both the measured
     * scene allowance and resident-memory allowance, continue once to the
     * next subpixel tier.  This is intentionally distinct from the
     * budget-debt retry below: every current request is satisfied, so there is
     * no blocked occurrence capable of arming that mechanism.
     *
     * Honor an explicit LoD scale setting as the user's quality decision.
     * The automatic tier is only the default-scale policy; a scale above one
     * already requests subpixel data, while a scale below one deliberately
     * trades fidelity for capacity. */
    struct bv_lod_policy viewPolicy;
    bv_lod_policy_init(&viewPolicy);
    this->d->viewAttachment->getLodPolicy(&viewPolicy);
    const bool defaultQualityScale = std::isfinite(viewPolicy.scale) &&
	std::fabs(viewPolicy.scale - 1.0) <= 1.0e-6;
    size_t presentedRenderCost = 0;
    const bool exactCompletedFrame = viewLodState &&
	viewLodState->lastCadPresentationFrameExact() &&
	viewLodState->lastCadPresentedRenderCost(presentedRenderCost) &&
	presentedRenderCost > 0;
    bool residentMemoryHeadroom = false;
    if (this->d->lodService) {
	const size_t residentLimit =
	    this->d->lodService->getResidentMeshLimit();
	const size_t residentBytes = this->d->lodService->
	    stableResidentMeshBytesForDiagnostics();
	const size_t reservedGrowth = this->d->lodService->
	    reservedResidentMeshGrowthBytesForDiagnostics();
	residentMemoryHeadroom = residentLimit == SIZE_MAX ||
	    (residentBytes <= residentLimit / 4 &&
	     reservedGrowth <= residentLimit / 4 -
		std::min(residentBytes, residentLimit / 4));
    }
    const bool stableTerminalContext = this->d->lodAutoSubmit &&
	this->d->lodService && !this->d->lodInteractive &&
	!this->d->lodGestureActive && defaultQualityScale &&
	!this->d->lodUseForcedCut &&
	this->d->lodCoveragePolicy.coverageComplete() &&
	activePayloads > 0 && memoryLimitedPayloads == 0 &&
	!this->d->lodSubmissionPending &&
	!this->d->lodBudgetPolicy.rescanAfterFrame() &&
	!this->d->lodFrameObligation.pending() &&
	!this->d->lodRetainedRefinementPending &&
	!this->d->lodPresentationPolicy.handoffPending() &&
	!this->d->lodStablePointProxyCalibrationPending &&
	!this->d->lodPointProxyTriangleRecoveryPending &&
	externalProducersSettled &&
	this->d->lodResultsPending.load() == 0 && activeTasks == 0 &&
	queuedResults == 0;
    const bool stableSubpixelContext = stableTerminalContext &&
	activePayloads == satisfiedPayloads &&
	this->d->lodInteractiveProgressiveCeiling < 0;
    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
	    this->d->lodViewRevision.value()) &&
	this->d->lodTargetPixelError > 0.5001f &&
	this->d->lodTargetPixelError <= 1.01f)
	bu_log("BObol LoD subpixel arm eligible=%d exact=%d default=%d "
	       "coverage=%d submit=%d rescan=%d refinement_frame=%d "
	       "retained=%d handoff=%d point=%d recovery=%d ceiling=%d "
	       "results=%d active_tasks=%zu queued_results=%zu "
	       "active_cost=%zu budget=%zu active_payloads=%zu "
	       "satisfied=%zu memory_limited=%zu memory_headroom=%d "
	       "last_ms=%.3f\n",
	       stableSubpixelContext ? 1 : 0, exactCompletedFrame ? 1 : 0,
	       defaultQualityScale ? 1 : 0,
	       this->d->lodCoveragePolicy.coverageComplete() ? 1 : 0,
	       this->d->lodSubmissionPending ? 1 : 0,
	       this->d->lodBudgetPolicy.rescanAfterFrame() ? 1 : 0,
	       this->d->lodFrameObligation.pending() ? 1 : 0,
	       this->d->lodRetainedRefinementPending ? 1 : 0,
	       this->d->lodPresentationPolicy.handoffPending() ? 1 : 0,
	       this->d->lodStablePointProxyCalibrationPending ? 1 : 0,
	       this->d->lodPointProxyTriangleRecoveryPending ? 1 : 0,
	       this->d->lodInteractiveProgressiveCeiling,
	       this->d->lodResultsPending.load(), activeTasks, queuedResults,
	       activePopulationCost, this->d->lodBudgetPolicy.currentBudget(),
	       activePayloads, satisfiedPayloads, memoryLimitedPayloads,
	       residentMemoryHeadroom ? 1 : 0,
	       this->d->lastRenderTimeNanoseconds / 1000000.0);

    /* A renderer-wide prefix ceiling is an O(1) interruption mechanism, not
     * a stable scene policy.  A deadline/static path can occasionally retire
     * its handoff after exhausting point aggregation while the last coherent
     * ceiling remains installed.  Do not report that state ready or remember
     * it as exact-view quality.  Convert the exact presented population into
     * occurrence-local cuts under the work budget which the completed frame
     * itself proved, then remove the ceiling through the normal certified
     * handoff.  This is a single bounded minimax transaction; the static
     * rejected/headroom witnesses prevent the failed richer population from
     * being retried afterward. */
    if (BObolLodStaticQualityPolicy::
	    terminalGlobalCeilingRequiresReconciliation(
		stableTerminalContext, exactCompletedFrame,
		this->d->lodInteractiveProgressiveCeiling,
		this->d->lodStaticOverscanActive != FALSE,
		this->d->lodStaticOverscanRejected != FALSE)) {
	const std::vector<SoBRLDatabaseSource *> sources =
	    controller_render_database_source_roots(this);
	if (!sources.empty()) {
	    const size_t reconciliationBudget =
		std::max<size_t>(1, presentedRenderCost);
	    this->d->lodPresentationPolicy.armHandoff(
		false, presentedRenderCost);
	    this->d->lodBudgetPolicy.requestPresentationReconciliation(
		reconciliationBudget);
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionRescanPending = FALSE;
	    this->d->lodSubmissionPending = TRUE;
	    this->d->lodPassAdmittedWork = FALSE;
	    this->d->lodRetainedRefinementPending = FALSE;
	    this->d->lodRetainedResidencyPending = FALSE;
	    this->d->lodRetainedRefinementCutAdvanced = FALSE;
	    this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	    this->d->lodBudgetPolicy.resetPass();
	    this->d->lodStaticOverscanActive = FALSE;
	    this->markProgressiveWorkPending();
	    this->requestRender("lod-terminal-ceiling-reconciliation");
	    return;
	}
    }

    /* Preserve only an exact, terminal presentation which actually met the
     * hard static-frame contract.  This evidence is keyed by the complete
     * camera/viewport/LoD signature and the controller's source-policy
     * domain.  It may seed a later return to this exact view, but the live
     * allocator and deadline recovery still own admission and correction. */
    if (stableTerminalContext && exactCompletedFrame &&
	this->d->lodViewSignatureValid) {
	BObolLodViewQualityHistory::RememberInputs remembered;
	remembered.view = this->d->lodViewSignature;
	remembered.domainRevision =
	    this->d->lodViewQualityDomainRevision;
	remembered.sceneAvailable = viewLodState && activePayloads > 0;
	remembered.quality.targetPixelError =
	    this->d->lodTargetPixelError;
	remembered.quality.progressiveCeiling =
	    this->d->lodInteractiveProgressiveCeiling;
	remembered.quality.pointProxyPixelThreshold =
	    this->d->lodPresentationPointProxyPixelThreshold;
	remembered.quality.maximumProjectedErrorPixels =
	    this->d->lodRetainedAdmissionQualityViewRevision ==
		this->d->lodViewRevision.value() &&
	    this->d->lodRetainedAdmissionQualityPolicyRevision ==
		this->d->lodPolicyRevision.value() &&
	    std::isfinite(
		this->d->lodRetainedAdmissionQualityPointProxyPixelThreshold) &&
	    std::fabs(
		this->d->lodRetainedAdmissionQualityPointProxyPixelThreshold -
		this->d->lodPresentationPointProxyPixelThreshold) <= 1.0e-6f ?
		this->d->lodRetainedAdmissionMaximumProjectedErrorPixels :
		std::numeric_limits<double>::infinity();
	remembered.quality.presentedRenderCost = presentedRenderCost;
	remembered.exactCompletedFrame = true;
	remembered.terminalPresentationComplete =
	    exactOccurrenceCoverage && presentedStructuralBoxes == 0 &&
	    terminalOccurrenceFailures == 0;
	remembered.producersSettled = externalProducersSettled;
	remembered.presentationDeadlineMet =
	    !this->d->stablePresentationFrameDeadlineNanoseconds ||
	    this->d->lastRenderTimeNanoseconds <=
		this->d->stablePresentationFrameDeadlineNanoseconds;
	remembered.resourcePressure =
	    this->d->lodResourcePolicy.anyPressure();
	(void)this->d->lodViewQualityHistory.remember(remembered);
    }

    /* The ordinary allocator budget may have been calibrated against an
     * early streaming frame and can lag far behind the exact terminal
     * presentation.  That made provider batching change final fidelity: the
     * same warm 150k scene stopped at one pixel when it discovered leaves
     * quickly, but reached one quarter pixel when slower publication happened
     * to produce more calibration frames.  An exact presentation of the
     * current active population is direct, schedule-independent evidence of
     * how much the hard static deadline can draw.  Use that proof to choose
     * the next tier, then raise the allocator only to the cost floor of the
     * tier actually selected below. */
    size_t staticQualityBudget =
	this->d->lodBudgetPolicy.currentBudget();
    if (stableSubpixelContext && exactCompletedFrame) {
	const size_t demonstratedPresentationLimit =
	    BObolLodQualityPolicy::staticPresentationRenderCostLimit(
		presentedRenderCost, this->d->lastRenderTimeNanoseconds,
		this->d->stablePresentationFrameDeadlineNanoseconds);
	const size_t demonstratedStaticBudget =
	    BObolLodQualityPolicy::incrementalSceneCostBudget(
		activePopulationCost, presentedRenderCost,
		demonstratedPresentationLimit);
	staticQualityBudget = std::max(
	    staticQualityBudget, demonstratedStaticBudget);
    }
    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
	    this->d->lodViewRevision.value()) && stableSubpixelContext)
	bu_log("BObol LoD static budget proof exact=%d "
	       "presented_cost=%zu active_cost=%zu "
	       "current_budget=%zu static_budget=%zu deadline_ms=%.3f "
	       "last_ms=%.3f\n",
	       exactCompletedFrame ? 1 : 0, presentedRenderCost,
	       activePopulationCost,
	       this->d->lodBudgetPolicy.currentBudget(), staticQualityBudget,
	       this->d->stablePresentationFrameDeadlineNanoseconds /
		   1000000.0,
	       this->d->lastRenderTimeNanoseconds / 1000000.0);
    const float stablePixelError = stableSubpixelContext ?
	BObolLodQualityPolicy::stablePixelError(
	    this->d->lodTargetPixelError, activePopulationCost,
	    staticQualityBudget,
	    this->d->lastRenderTimeNanoseconds,
	    this->d->staticQualityTargetFps(),
	    exactCompletedFrame, residentMemoryHeadroom) :
	this->d->lodTargetPixelError;
    if (stablePixelError + 1.0e-6f < this->d->lodTargetPixelError) {
	const size_t nextTierBudget =
	    BObolLodQualityPolicy::pixelErrorRenderCostFloor(
		activePopulationCost, this->d->lodTargetPixelError,
		stablePixelError);
	if (nextTierBudget <= staticQualityBudget)
	    this->d->lodBudgetPolicy.raiseCurrentBudget(nextTierBudget);
	this->d->lodTargetPixelError = stablePixelError;
	this->advanceLodPolicyRevision();
	this->d->lodLastSubmittedPolicyRevision.reset();
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPending = TRUE;
	this->d->lodSubmissionRescanPending = FALSE;
	this->markProgressiveWorkPending();
	this->requestRender("lod-stable-subpixel");
	return;
    }
    const bool actionableQualityDebt =
	activePayloads > satisfiedPayloads &&
	activePayloads - satisfiedPayloads > memoryLimitedPayloads;

    /* The preferred quiet cadence is an initial convergence target, not a
     * permanent fidelity ceiling for an event-driven static image.  Once the
     * ordinary allocation is terminal, exact, fully covered, and free of
     * memory pressure, recompute the scene-wide importance allocation using
     * the separate hard static-frame deadline.  The completed framebuffer is
     * retained without redraw; any later input immediately leaves this mode
     * and restores the motion budget from the existing occurrence cuts.
     *
     * This is one phase transition, not an open-ended headroom probe.  The
     * ordinary deadline recovery bounds a discrete PoP jump which turns out
     * to be too expensive, while pixel demand and resident-memory limits
     * remain authoritative terminal conditions. */
    const uint64_t preferredQuietNanoseconds =
	this->d->lodStableTargetFps > 0.0f ?
	    static_cast<uint64_t>(1000000000.0L /
		static_cast<long double>(this->d->lodStableTargetFps)) : 0;
    const bool staticOverscanEligible = stableTerminalContext &&
	exactCompletedFrame && actionableQualityDebt &&
	this->d->lodBudgetPolicy.stableBudgetLimited() &&
	!this->d->lodStaticOverscanActive &&
	!this->d->lodStaticOverscanRejected &&
	!this->d->lodResourcePolicy.anyPressure() &&
	this->d->stablePresentationFrameDeadlineNanoseconds >
	    preferredQuietNanoseconds;
    if (staticOverscanEligible) {
	this->d->lodStaticOverscanActive = TRUE;
	this->d->lodStaticOverscanLeapAvailable = TRUE;
	BObolViewLodState *presentationState =
	    this->d->viewAttachment->getViewLodState();
	const int priorCeiling =
	    this->d->lodInteractiveProgressiveCeiling;
	const int stagedCeiling =
	    BObolLodStaticQualityPolicy::stagedProgressiveCeiling(
		priorCeiling,
		presentationState ? presentationState->
		    maximumActiveProgressiveCut() : -1,
		this->d->lodPresentationPointProxyPixelThreshold,
		activePayloads);
	const bool retainAggregatedPresentation =
	    BObolLodStaticQualityPolicy::retainAggregatedPresentation(
		priorCeiling,
		this->d->lodPresentationPointProxyPixelThreshold,
		activePayloads,
		exactStructuralProjection &&
		    structuralProjection.visibleCount > 0);
	if (stagedCeiling >= 0) {
	    /* The calibrated point population is part of the current exact frame,
	     * not a fidelity failure to undo before improving its visible meshes.
	     * Raise only the renderer prefix.  Successful frames continue this
	     * staircase in completeRenderTiming(); the first hard miss falls back
	     * to the preceding exact ceiling and reconciles that cost into the
	     * occurrence-local importance allocation. */
	    this->d->lodStaticOverscanLeapAvailable = FALSE;
	    this->d->lodInteractiveProgressiveCeiling = stagedCeiling;
	    presentationState->setCadPresentationProgressiveCutCeiling(
		stagedCeiling);
	    presentationState->setCadPresentationCameraMotionFrameReuse(FALSE);
	    this->d->lodHeadroomPolicy.cancelRetry();
	    this->d->lodDiscretePopulationTrialAvailable = TRUE;
	    this->markProgressiveWorkPending();
	    this->requestRender("lod-static-overscan-staged");
	    return;
	}
	/* A quiet handoff may have retained the responsive presentation ceiling
	 * learned during zoom.  That ceiling is valuable while input is active,
	 * but leaving it installed here prevents an already resident richer cut
	 * from ever participating in the bounded static-frame trial.  Remove only
	 * the renderer-side ceiling: occurrence demand and immutable resident PoP
	 * storage stay unchanged.  If the richer frame misses the hard static
	 * deadline, notePresentationRenderInterrupted() restores a one-cut-lower
	 * ceiling and the previous completed framebuffer remains available. */
	if (!retainAggregatedPresentation) {
	    this->d->lodInteractiveProgressiveCeiling = -1;
	    if (presentationState)
		presentationState->setCadPresentationProgressiveCutCeiling(-1);
	}
	if (presentationState)
	    presentationState->setCadPresentationCameraMotionFrameReuse(FALSE);
	/* The preferred-cadence point bracket is likewise a reversible
	 * presentation limit, not immutable geometry or a static fidelity policy.
	 * Test the one-pixel occurrence population directly under the interruptible
	 * hard static deadline.  If it is too expensive, the deadline handler keeps
	 * the preceding complete framebuffer and establishes the unsafe side of the
	 * point bracket in one observation.  Walking 64 -> 48 -> ... -> 1 after
	 * button-up is both slower and visibly distracting, and can strand a view at
	 * the 20 Hz cut even when its full one-pixel frame fits comfortably inside
	 * the 100 ms event-driven allowance. */
	const SbBool restoredOnePixelPopulation =
	    !retainAggregatedPresentation &&
	    BObolLodStaticQualityPolicy::onePixelTrialRequired(
		this->d->lodPresentationPointProxyPixelThreshold) ?
		TRUE : FALSE;
	if (restoredOnePixelPopulation) {
	    this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
	    this->d->lodPointProxyCalibrationPolicy.reset();
	    if (presentationState)
		presentationState->setCadPresentationPointProxyPixelThreshold(
		    1.0f);
	    this->d->lodStablePointProxyCalibrationPending = TRUE;
	}
	if (retainAggregatedPresentation)
	    this->d->lodPresentationPolicy.armHandoff(
		false, presentedRenderCost);
	else
	    this->d->lodPresentationPolicy.cancelHandoff();
	this->d->lodHeadroomPolicy.cancelRetry();
	this->d->lodBudgetPolicy.clearBudgetLimit();
	this->d->lodBudgetPolicy.resetProbeSeries();
	this->d->lodBudgetPolicy.resetOverloadRecovery();
	this->d->lodBudgetPolicy.resetPass();
	this->d->lodBudgetPolicy.requestRetainedReallocation(false);
	/* This is another pass in the same camera/policy/source epoch.  The
	 * explicit pending cursor is sufficient to bypass the completed-pass fast
	 * path.  Clearing the epoch witness as well makes the wrapper classify the
	 * already-pending cursor as a view change during submission and append an
	 * unnecessary full rescan after every allocation. */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPending = TRUE;
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodDiscretePopulationTrialAvailable = TRUE;
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
		this->d->lodViewRevision.value()))
	    bu_log("BObol LoD static overscan armed prior_ceiling=%d "
		   "active_max=%d active_cost=%zu budget=%zu "
		   "target_fps=%.3f last_ms=%.3f\n",
		   priorCeiling,
		   presentationState ? presentationState->
		       maximumActiveProgressiveCut() : -1,
		   activePopulationCost,
		   this->d->lodBudgetPolicy.currentBudget(),
		   this->d->quietAllocationTargetFps(),
		   this->d->lastRenderTimeNanoseconds / 1000000.0);
	this->markProgressiveWorkPending();
	this->requestRender("lod-static-overscan");
	return;
    }
    /* A terminal planning pass can finish while the motion-to-stable or
     * small-part presentation barrier still owns the next frame.  Arm one
     * exact, unchanged replay from either side of that barrier.  The
     * headroom policy makes the (view, policy, active population) witness
     * one-shot, so an actually capacity-limited population cannot repaint in
     * a loop. */
    const bool eligible = this->d->lodAutoSubmit && this->d->lodService &&
	!this->d->lodInteractive && !this->d->lodGestureActive &&
	BObolLodStaticQualityPolicy::ordinaryHeadroomAllowed(
	    this->d->lodStaticOverscanActive != FALSE,
	    this->d->lodStaticOverscanRejected != FALSE) &&
	this->d->lodBudgetPolicy.stableBudgetLimited() &&
	!this->d->lodSubmissionPending &&
	!this->d->lodBudgetPolicy.rescanAfterFrame() &&
	!this->d->lodFrameObligation.pending() &&
	!this->d->lodRetainedRefinementPending &&
	!this->d->lodPresentationPolicy.handoffPending() &&
	!this->d->lodStablePointProxyCalibrationPending &&
	!this->d->lodPointProxyTriangleRecoveryPending &&
	externalProducersSettled &&
	this->d->lodInteractiveProgressiveCeiling < 0 &&
	this->d->lodBudgetPolicy.currentBudget() != SIZE_MAX &&
	actionableQualityDebt &&
	this->d->lodResultsPending.load() == 0 && activeTasks == 0 &&
	queuedResults == 0;
    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
	    this->d->lodViewRevision.value()))
	bu_log("BObol LoD headroom arm eligible=%d interactive=%d gesture=%d "
	       "limited=%d submit=%d rescan=%d refinement_frame=%d "
	       "retained=%d handoff=%d point=%d recovery=%d ceiling=%d "
	       "results=%d active_tasks=%zu queued_results=%zu active_cost=%zu "
	       "budget=%zu quality_debt=%d active_payloads=%zu satisfied=%zu "
	       "memory_limited=%zu pending=%d\n",
	       eligible ? 1 : 0, this->d->lodInteractive ? 1 : 0,
	       this->d->lodGestureActive ? 1 : 0,
	       this->d->lodBudgetPolicy.stableBudgetLimited() ? 1 : 0,
	       this->d->lodSubmissionPending ? 1 : 0,
	       this->d->lodBudgetPolicy.rescanAfterFrame() ? 1 : 0,
	       this->d->lodFrameObligation.pending() ? 1 : 0,
	       this->d->lodRetainedRefinementPending ? 1 : 0,
	       this->d->lodPresentationPolicy.handoffPending() ? 1 : 0,
	       this->d->lodStablePointProxyCalibrationPending ? 1 : 0,
	       this->d->lodPointProxyTriangleRecoveryPending ? 1 : 0,
	       this->d->lodInteractiveProgressiveCeiling,
	       this->d->lodResultsPending.load(), activeTasks, queuedResults,
	       activePopulationCost, this->d->lodBudgetPolicy.currentBudget(),
	       actionableQualityDebt ? 1 : 0, activePayloads,
	       satisfiedPayloads, memoryLimitedPayloads,
	       this->d->lodHeadroomPolicy.retryPending() ? 1 : 0);
    if (!eligible)
	return;

    if (!activePopulationCost ||
	!this->d->lodHeadroomPolicy.armRetry(
	    this->d->lodViewRevision, this->d->lodPolicyRevision,
	    activePopulationCost,
	    this->d->lodBudgetPolicy.currentBudget()))
	return;

    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
	    this->d->lodViewRevision.value()))
	bu_log("BObol LoD headroom armed view=%llu policy=%llu "
	       "active_cost=%zu budget=%zu\n",
	       (unsigned long long)this->d->lodViewRevision.value(),
	       (unsigned long long)this->d->lodPolicyRevision.value(),
	       activePopulationCost, this->d->lodBudgetPolicy.currentBudget());

    this->markProgressiveWorkPending();
    this->requestRender("lod-calibrated-headroom-probe");
}

void
BObolViewController::resumeLodAfterRetainedRecovery(void)
{
    if (!this->d->lodBudgetPolicy.
	    confirmRetainedRecoveryPresentation(true))
	return;

    /* The recovery ceiling is a one-frame guard, not a fidelity policy.  Its
     * coherent successor is now visible at the current point/mesh split, so
     * let measured headroom grow from the retained population in bounded
     * screen-priority waves. */
    this->d->lodSubmissionSourceIndex = 0;
    this->d->lodSubmissionEntryOffset = 0;
    this->d->clearLodSubmissionPlan();
    this->d->lodSubmissionPending = TRUE;
    this->d->lodSubmissionRescanPending = FALSE;
    this->markProgressiveWorkPending();
}

void
BObolViewController::completeRenderTiming(uint64_t startedNanoseconds,
	SbBool lodCapacityRelevant)
{
    const uint64_t now = this->beginRenderTiming();
    if (!startedNanoseconds || now <= startedNanoseconds)
	return;

    const uint64_t elapsed = now - startedNanoseconds;
    if (elapsed >= 30000000000ULL)
	return;

    this->d->lastRenderTimeNanoseconds = elapsed;
    this->d->smoothedRenderTimeNanoseconds =
	this->d->smoothedRenderTimeNanoseconds ?
	(this->d->smoothedRenderTimeNanoseconds * 9 + elapsed) / 10 : elapsed;
    /* The completed traversal is the commit edge for any interrupted CAD
     * replay.  Release queued immutable provider/result publications before
     * capacity-specific early returns so presentation-only classifier frames
     * cannot strand the owner-thread mutation gate. */
    this->d->lodInterruptedPresentationReplayPending = FALSE;

    /* Presentation-only frames still contribute user-facing frame-time
     * telemetry and, when exact, are real presentations of the current CAD
     * population.  They therefore retire publication/refinement barriers,
     * but must not satisfy/cancel capacity probes or teach the retained
     * geometry allocator from a one-time style patch.
     *
     * There is one distinct no-sample transition: a budget calibration may
     * have been armed by a coverage wave whose managed mesh population has
     * since fallen to zero while structural proxies remain.  No future replay
     * of that population can become capacity relevant.  Retire the obsolete
     * probe and resume admission immediately; otherwise the generic
     * refinement-frame liveness guard repaints the same boxes forever. */
    if (!lodCapacityRelevant) {
	if (!this->isLodPresentationCapacityRelevant() &&
	    this->d->lodBudgetPolicy.rescanAfterFrame()) {
	    const BObolLodBudgetPolicy::CompletedFrameDecision calibration =
		this->d->lodBudgetPolicy.
		    retireUnmeasurableCalibrationFrame();
	    if (calibration.restartSubmission) {
		this->d->lodSubmissionSourceIndex = 0;
		this->d->lodSubmissionEntryOffset = 0;
		this->d->clearLodSubmissionPlan();
		this->d->lodSubmissionPending = TRUE;
		this->d->lodSubmissionRescanPending = FALSE;
		if (this->d->lodCoveragePolicy.required() &&
		    !this->d->lodCoveragePolicy.coverageComplete() &&
		    !this->d->lodCoveragePolicy.active())
		    this->d->lodCoveragePolicy.activate(false);
		if (this->d->lodCoveragePolicy.active())
		    this->d->lodCoveragePolicy.clearPassCounters();
		this->markProgressiveWorkPending();
	    }
	}
	/* The initial point/box classifier deliberately contains no managed mesh
	 * population, so it is a presentation transaction rather than a geometry
	 * capacity sample.  Its exact frame must nevertheless retire the admission
	 * barrier.  Leaving that acknowledgement in the capacity-only path below
	 * makes a large structural scene repaint forever with zero submitted mesh
	 * tasks: the very absence of meshes prevents the latch which admits them
	 * from being consumed. */
	const SbBool admissionFrameWasPending =
	    this->d->lodAdmissionPointProxyFramePending;
	if (admissionFrameWasPending) {
	    const BObolViewLodState *presentationState =
		this->d->viewAttachment ?
		this->d->viewAttachment->getViewLodState() : NULL;
	    size_t presentedSubpixelOccurrences = 0;
	    size_t presentedStructuralBoxes = 0;
	    const SbBool exactOccurrenceClassification =
		presentationState &&
		presentationState->lastCadPresentationOccurrenceCoverage(
		    presentedSubpixelOccurrences,
		    presentedStructuralBoxes);
	    const SbBool exactAdmissionClassification =
		presentationState &&
		presentationState->hasCadPresentationAssemblies() &&
		presentationState->lastCadPresentationFrameExact() &&
		exactOccurrenceClassification;
	    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_POINT",
		    this->d->lodViewRevision.value()))
		bu_log("BObol LoD presentation-only admission point frame "
		       "exact=%d occurrence_exact=%d subpixel=%zu boxes=%zu\n",
		       exactAdmissionClassification ? 1 : 0,
		       exactOccurrenceClassification ? 1 : 0,
		       presentedSubpixelOccurrences, presentedStructuralBoxes);
	    if (exactAdmissionClassification)
		this->d->lodAdmissionPointProxyFramePending = FALSE;
	    else
		this->requestRender("lod-admission-point-replay");
	}
	this->completePresentationBarrier(elapsed);
	if (admissionFrameWasPending &&
	    !this->d->lodAdmissionPointProxyFramePending)
	    this->armStableLodHeadroomProbeIfReady();
	return;
    }

    /* A completed presentation ends any deadline-backoff sequence.  Aborted
     * traversals never enter the persistent capacity estimator: they do not
     * publish an exact presented-work denominator. */
    this->d->consecutiveInterruptedPresentationFrames = 0;

    /* Calibrate geometry throughput from frames that actually traversed the
     * current aggregate cut.  This is deliberately scene-level: per-object
     * projected error says what would look good, but only observed aggregate
     * work says what this renderer, viewport, draw mode, and machine can
     * sustain.
     *
     * Motion and quiet views have separate calibrations.  A settled retained
     * frame can be exceptionally cheap (especially when one mesh is instanced
     * thousands of times), while a motion frame also pays for cut selection
     * and presentation updates.  Mixing those samples taught the interaction
     * policy a budget it could not sustain. */
    const BObolViewLodState *calibrationState = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    /* One complete-frame sample supplies triangle, GPU timer, replay, and
     * resource-pressure consumers.  Sampling here avoids four independent
     * source-presentation walks in the controller and HUD paths. */
    if (calibrationState)
	calibrationState->refreshCadPresentationFrameStatus();
    BObolViewLodState::CadGpuResourceStatus frameGpuResources;
    const SbBool haveFrameGpuResources = calibrationState &&
	calibrationState->cadGpuResourceStatus(frameGpuResources);
    const SbBool gpuPressure = haveFrameGpuResources &&
	frameGpuResources.memoryPressure ? TRUE : FALSE;
    SbBool cpuPressure = FALSE;
    if (this->d->lodService) {
	const size_t residentLimit =
	    this->d->lodService->getResidentMeshLimit();
	const size_t stableResident = this->d->lodService->
	    stableResidentMeshBytesForDiagnostics();
	const size_t reservedGrowth = this->d->lodService->
	    reservedResidentMeshGrowthBytesForDiagnostics();
	cpuPressure = residentLimit != SIZE_MAX &&
	    (stableResident > residentLimit ||
	     reservedGrowth > residentLimit -
		std::min(stableResident, residentLimit)) ? TRUE : FALSE;
    }
    const BObolLodResourcePolicy::Decision resourceDecision =
	this->d->lodResourcePolicy.observe(
	    cpuPressure != FALSE, gpuPressure != FALSE,
	    this->d->lodAutoSubmit && this->d->lodService);
    if (resourceDecision.queueRecovery) {
	/* Resource recovery obeys the same coverage and interaction safety gates
	 * as ordinary stable compaction.  The coordinator owns the pressure edge;
	 * this callback only schedules its decision. */
	const int64_t nowMicroseconds = bu_gettime();
	this->d->lodCompactionPolicy.requestImmediate(nowMicroseconds);
	this->markProgressiveWorkPending();
    }
    const size_t activeCalibrationFaces = calibrationState ?
	calibrationState->activeFaceCount() : 0;
    const size_t activeCalibrationCost = calibrationState ?
	calibrationState->activeRenderCost() : 0;
    size_t calibrationCost = activeCalibrationCost;
    const auto costForPresentedFaces = [activeCalibrationFaces,
	activeCalibrationCost](size_t presentedFaces) {
	if (!activeCalibrationFaces || !activeCalibrationCost)
	    return presentedFaces;
	const long double scaled =
	    static_cast<long double>(activeCalibrationCost) *
	    static_cast<long double>(presentedFaces) /
	    static_cast<long double>(activeCalibrationFaces);
	return scaled >= static_cast<long double>(SIZE_MAX) ? SIZE_MAX :
	    std::max(presentedFaces, static_cast<size_t>(scaled));
    };
    size_t presentedCadPrimitives = 0;
    const SbBool measuredCadPresentation = calibrationState &&
	calibrationState->lastCadPresentedPrimitiveCount(
	    presentedCadPrimitives);
    if (measuredCadPresentation)
	calibrationCost = costForPresentedFaces(presentedCadPrimitives);
    size_t presentedCadRenderCost = 0;
    const SbBool measuredCadRenderCost = calibrationState &&
	calibrationState->lastCadPresentedRenderCost(
	    presentedCadRenderCost);
    if (measuredCadRenderCost)
	calibrationCost = presentedCadRenderCost;
    size_t gpuCadFaces = 0;
    uint64_t gpuCadNanoseconds = 0;
    uint64_t gpuCadSampleSerial = 0;
    float gpuCadPointProxyPixelThreshold = 1.0f;
    const SbBool measuredGpuCadPresentation = calibrationState &&
	calibrationState->lastCadGpuMeasurement(
	    gpuCadFaces, gpuCadNanoseconds, gpuCadSampleSerial,
	    gpuCadPointProxyPixelThreshold);
    const SbBool newGpuCadMeasurement =
	measuredGpuCadPresentation &&
	gpuCadSampleSerial != this->d->lodLastCadGpuSampleSerial;
    if (newGpuCadMeasurement) {
	this->d->lodLastCadGpuSampleSerial = gpuCadSampleSerial;
	this->d->lodLastCadGpuTimeNanoseconds = gpuCadNanoseconds;
	calibrationCost = costForPresentedFaces(gpuCadFaces);
    }
    const SbBool haveCadPresentationAssemblies = calibrationState &&
	calibrationState->hasCadPresentationAssemblies();
    /* Admission budgets are denominated in BObolViewLodState's indexed-mesh
     * render-cost units.  SoCADAssembly's completed-work record intentionally
     * describes backend work instead: on the ordinary OSMesa route its
     * position count may be the expanded submitted stream.  Treating that
     * backend counter as allocator currency made a 75 ms Hubble frame appear
     * to contain nearly twice its retained cost; each recovery then removed
     * only a few entries and repeated for hundreds of frames.
     *
     * When no presentation-only PoP ceiling is hiding the retained cut, the
     * exact active population is the authoritative cost in allocator units.
     * A point threshold of one pixel is part of the normal presentation
     * contract, not a different retained allocation.  Keep the expanded work
     * record for diagnostics and GPU/backend analysis, but pair elapsed time
     * with activeCalibrationCost for admission calibration. */
    size_t presentedSubpixelOccurrences = 0;
    size_t presentedStructuralBoxes = 0;
    const SbBool exactOccurrenceClassification = calibrationState &&
	calibrationState->lastCadPresentationOccurrenceCoverage(
	    presentedSubpixelOccurrences, presentedStructuralBoxes);
    const SbBool admissionPointProxyFrameWasPending =
	this->d->lodAdmissionPointProxyFramePending;
    if (admissionPointProxyFrameWasPending) {
	const SbBool exactAdmissionClassification =
	    haveCadPresentationAssemblies &&
	    calibrationState->lastCadPresentationFrameExact() &&
	    exactOccurrenceClassification;
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_POINT",
		this->d->lodViewRevision.value())) {
	    static std::atomic<unsigned int> pointFrameTraceCount(0);
	    const unsigned int traceIndex = pointFrameTraceCount.fetch_add(1);
	    if (traceIndex < 256)
		bu_log("BObol LoD admission point frame exact=%d "
		       "assemblies=%d frame_exact=%d occurrence_exact=%d "
		       "subpixel=%zu boxes=%zu providers=%zu "
		       "discovery_threshold=%.3f stable_threshold=%.3f\n",
		       exactAdmissionClassification ? 1 : 0,
		       haveCadPresentationAssemblies ? 1 : 0,
		       calibrationState->lastCadPresentationFrameExact() ? 1 : 0,
		       exactOccurrenceClassification ? 1 : 0,
		       presentedSubpixelOccurrences, presentedStructuralBoxes,
		       this->d->progressiveProviderPendingCount,
		       this->d->lodDiscoveryPointProxyPixelThreshold,
		       this->d->lodPresentationPointProxyPixelThreshold);
	}
	if (exactAdmissionClassification) {
	    this->d->lodAdmissionPointProxyFramePending = FALSE;
	} else {
	    /* A publication or scene transaction may make the requested frame
	     * inexact.  Preserve an explicit successor-frame witness; otherwise
	     * mesh admission can remain paused after the provider timer becomes
	     * idle. */
	    this->requestRender("lod-discovery-point-replay");
	}
    }
    const size_t terminalOccurrenceFailures = calibrationState ?
	calibrationState->cadOccurrenceTerminalFailureCount() : 0;
    const size_t unresolvedStructuralOccurrences =
	exactOccurrenceClassification &&
	presentedStructuralBoxes > terminalOccurrenceFailures ?
	    presentedStructuralBoxes - terminalOccurrenceFailures : 0;
    BObolViewLodState::CadStructuralProjectionHistogram
	calibrationStructuralProjection;
    const SbBool exactStructuralProjection = calibrationState &&
	calibrationState->lastCadStructuralProjectionHistogram(
	    calibrationStructuralProjection);
    /*
     * An aggregate point backed only by a structural fallback is not a
     * hidden richer mesh.  Lowering the point threshold re-exposes its box,
     * so the complete projected structural population (not just the boxes
     * visible at the current threshold) is a quality-transition barrier.
     * Keep the last proven point cut while those records remain; mesh result
     * publication may replace them in place, after which ordinary relaxation
     * can continue without a box flash.
     */
    const size_t structuralFallbackOccurrences = exactStructuralProjection ?
	calibrationStructuralProjection.visibleCount :
	(presentedStructuralBoxes ? presentedStructuralBoxes : SIZE_MAX);
    const SbBool exactRetainedPopulation =
	BObolLodPointProxyCalibrationPolicy::onePixelPresentationReady(
	    haveCadPresentationAssemblies != FALSE,
	    calibrationState &&
		calibrationState->lastCadPresentationFrameExact(),
	    exactOccurrenceClassification != FALSE,
	    presentedSubpixelOccurrences, presentedStructuralBoxes,
	    this->d->lodInteractiveProgressiveCeiling,
	    this->d->lodPresentationPointProxyPixelThreshold,
	    activeCalibrationCost) ? TRUE : FALSE;
    if (exactRetainedPopulation) {
	calibrationCost = activeCalibrationCost;
	/* A box-free, one-pixel frame proves the aggregate presentation cut, but
	 * it does not by itself prove that the protected prominent-geometry floor
	 * was the population on screen.  Deadline recovery can temporarily
	 * present a cheaper occurrence plan with those same aggregate properties.
	 * Reset floor-miss evidence only when the exact current allocation
	 * certificate covers the complete occurrence population, the CAD draw
	 * supplied a completed-work record, and the certificate carries the
	 * protected-floor identity and cost which are being defended. */
	const BObolRetainedAllocationResult &floorAllocation =
	    this->d->lodRetainedAllocationCertificate;
	const bool exactProtectedFloor = calibrationState &&
	    measuredCadRenderCost &&
	    floorAllocation.allocationPlanSerial != 0 &&
	    floorAllocation.allocationPlanSerial ==
		calibrationState->activeCadAllocationPlan() &&
	    floorAllocation.viewRevision == this->d->lodViewRevision.value() &&
	    floorAllocation.policyRevision ==
		this->d->lodPolicyRevision.value() &&
	    std::isfinite(floorAllocation.pointProxyPixelThreshold) &&
	    std::fabs(floorAllocation.pointProxyPixelThreshold -
		this->d->lodPresentationPointProxyPixelThreshold) <= 1.0e-6f &&
	    calibrationState->cadAllocationPlanCoversCurrentPopulation(
		floorAllocation.allocationPlanSerial,
		floorAllocation.viewRevision, floorAllocation.policyRevision,
		floorAllocation.fixedCadPresentationCost);
	this->d->lodBudgetPolicy.noteRetainedQualityFloorMet(
	    exactProtectedFloor, floorAllocation.protectedFloorSignature,
	    activeCalibrationCost);
    }
    /* activeRenderCost is retained demand, not necessarily work submitted by
     * this completed traversal.  An empty/replaced CAD plan can legitimately
     * publish zero exact primitives for one frame while the prior occurrence
     * population is still active.  Dividing that retained population by the
     * tiny empty-frame time creates fictitious throughput and reopens a cut
     * which just missed its deadline. */
    const SbBool exactCalibrationWork =
	!haveCadPresentationAssemblies || measuredCadRenderCost ||
	newGpuCadMeasurement;
    this->d->lodLastRenderWasPreparedCadReplay =
	!measuredCadPresentation ||
	calibrationState->lastCadPresentationUsedPreparedReplay();
    const auto addUploadBytes = [](uint64_t left, uint64_t right) {
	return right > UINT64_MAX - left ? UINT64_MAX : left + right;
    };
    uint64_t geometryUploadBytes = 0;
    if (haveFrameGpuResources) {
	geometryUploadBytes = addUploadBytes(geometryUploadBytes,
	    frameGpuResources.ordinaryPartFullUploadBytes);
	geometryUploadBytes = addUploadBytes(geometryUploadBytes,
	    frameGpuResources.ordinaryPartSuffixUploadBytes);
	geometryUploadBytes = addUploadBytes(geometryUploadBytes,
	    frameGpuResources.triangleAtlasFullUploadBytes);
	geometryUploadBytes = addUploadBytes(geometryUploadBytes,
	    frameGpuResources.triangleAtlasSuffixUploadBytes);
    }
    const SbBool geometryUploadChanged = haveFrameGpuResources &&
	(!this->d->lodLastCadGpuGeometryUploadBytesValid ||
	 geometryUploadBytes != this->d->lodLastCadGpuGeometryUploadBytes) ?
	    TRUE : FALSE;
    if (haveFrameGpuResources) {
	this->d->lodLastCadGpuGeometryUploadBytes = geometryUploadBytes;
	this->d->lodLastCadGpuGeometryUploadBytesValid = TRUE;
    }
    /* Prepared multi-draw replay is one reusable presentation path, but it
     * is not the only one.  OSMesa deliberately uses retained ordinary VBOs;
     * once their cumulative upload counters stop changing, an unchanged
     * traversal is an equally valid sustainable-throughput sample. */
    this->d->lodLastRenderWasReusableCadPresentation =
	this->d->lodLastRenderWasPreparedCadReplay ||
	(haveFrameGpuResources && !geometryUploadChanged) ? TRUE : FALSE;
    uint64_t calibrationElapsed =
	newGpuCadMeasurement ? gpuCadNanoseconds : elapsed;
    const SbBool gpuPointProxySampleMatchesCurrent =
	newGpuCadMeasurement &&
	std::fabs(gpuCadPointProxyPixelThreshold -
	    this->d->lodPresentationPointProxyPixelThreshold) <= 0.01f ?
	    TRUE : FALSE;
    /* Timer queries complete on a later render.  Use them for the point-cut
     * bracket only when the query carries the same threshold as the cut being
     * evaluated; otherwise the just-completed CPU traversal is the current
     * evidence.  Scene throughput still uses the query's paired face count. */
    const uint64_t pointCalibrationElapsed =
	gpuPointProxySampleMatchesCurrent ? gpuCadNanoseconds : elapsed;
    /* Host presentation cadence is deliberately not renderer capacity.  A
     * slow mouse, scripted pause, compositor cap, or event-loop scheduling
     * gap cannot be improved by discarding geometry.  CPU traversal duration
     * and the asynchronous GL timer above are the actionable work samples. */
    /* Result publication may allocate/upload newly decoded geometry during
     * this traversal.  That one-time cost is relevant to the cooldown below,
     * but not to the sustainable FPS of the retained cut.  Quiet,
     * unchanged calibration probes measure the latter.  Interactive samples
     * keep the full cost because upload latency is part of the user's motion
     * experience. */
    const SbBool transientStablePublication =
	this->d->lodPublicationPolicy.framePending() &&
	!this->d->lodInteractive && !this->d->lodGestureActive;
    const SbBool transientStablePreparation =
	!this->d->lodLastRenderWasReusableCadPresentation &&
	!this->d->lodInteractive && !this->d->lodGestureActive;
    const SbBool reusableCadWorkSample =
	this->d->lodLastRenderWasReusableCadPresentation &&
	((haveCadPresentationAssemblies && measuredCadRenderCost) ||
	 (!haveCadPresentationAssemblies && activeCalibrationCost > 0)) ?
	    TRUE : FALSE;
    if (calibrationCost > 0 && calibrationElapsed > 0 &&
	exactCalibrationWork &&
	(newGpuCadMeasurement ||
	 (!transientStablePublication && !transientStablePreparation))) {
	const long double sample =
	    static_cast<long double>(calibrationCost) * 1000000000.0L /
	    static_cast<long double>(calibrationElapsed);
	if (std::isfinite(sample) && sample > 0.0L) {
	    const SbBool interactiveSample =
		this->d->lodInteractive || this->d->lodGestureActive;
	    long double &calibration = interactiveSample ?
		this->d->lodInteractiveCalibratedRenderCostPerSecond :
		this->d->lodStableCalibratedRenderCostPerSecond;
	    if (calibration <= 0.0L) {
		calibration = sample;
	    } else if (interactiveSample && sample < calibration) {
		/* Missing an interaction target is immediately visible to the
		 * user.  Adopt the observed lower ceiling immediately in that
		 * case.  A fast underloaded cut is not a lower throughput proof:
		 * fixed frame overhead divided into 44 or 118 triangles produces
		 * an arbitrarily small faces/second number and used to lock OSMesa
		 * permanently at that cut despite ample frame headroom. */
		const long double targetNanoseconds =
		    this->d->lodInteractiveTargetFps > 0.0f ?
		    1000000000.0L /
			static_cast<long double>(
			    this->d->lodInteractiveTargetFps) : 0.0L;
		const long double priorFrameBudget =
		    this->d->lodInteractiveTargetFps > 0.0f ?
		    calibration /
			static_cast<long double>(
			    this->d->lodInteractiveTargetFps) : 0.0L;
		const SbBool materiallyLoaded =
		    priorFrameBudget <= 0.0L ||
		    static_cast<long double>(calibrationCost) >=
			priorFrameBudget * 0.10L;
		if (targetNanoseconds > 0.0L &&
		    materiallyLoaded &&
		    static_cast<long double>(calibrationElapsed) >
			targetNanoseconds)
		    calibration = sample;
	    } else if (!interactiveSample && newGpuCadMeasurement &&
		sample < calibration) {
		/* A completed timer query measures the driver's actual queued CAD
		 * work rather than CPU command submission or an unrelated window
		 * composition delay.  A lower GPU ceiling is therefore actionable
		 * immediately; damping it recreates several known-bad long frames
		 * while the estimate slowly catches up. */
		calibration = sample;
	    } else if (!interactiveSample && sample < calibration) {
		/* A quiet event-driven view includes occasional screenshots,
		 * buffer readbacks, window composition, and cache activity.  One
		 * such outlier must not make admitted leaves disappear and then
		 * return on the next ordinary frame.  Reduce stable capacity only
		 * after a frame materially misses its target, and then do so with
		 * damping. */
		const long double targetNanoseconds =
		    this->d->quietAllocationTargetFps() > 0.0f ?
		    1000000000.0L /
			static_cast<long double>(
			    this->d->quietAllocationTargetFps()) : 0.0L;
		if (targetNanoseconds > 0.0L &&
		    static_cast<long double>(calibrationElapsed) >
			targetNanoseconds * 1.20L)
		    calibration = calibration * 0.80L + sample * 0.20L;
	    } else {
		const long double oldWeight = interactiveSample ? 0.98L : 0.90L;
		calibration =
		    calibration * oldWeight + sample * (1.0L - oldWeight);
	    }

	    /* A quality probe is a scale-changing presentation, but an
	     * completed frame which already meets the quieter stable deadline is
	     * a strict lower-bound proof for stable admission.  This remains true
	     * when that frame extended a resident VBO: total elapsed time includes
	     * the one-time upload, making the proof more conservative, and the
	     * completed prefix is resident afterward.  Transfer only exact
	     * submitted work: an occurrence may retain a much richer prefix behind
	     * the renderer ceiling, so activeRenderCost is not a safe substitute
	     * on compatibility renderers.  The scalar policy owns the deadline and
	     * throughput arithmetic. */
	    /* A triangle count is exact in face units but cannot safely stand in
	     * for the scheduler's weighted point/normal/occurrence units.  Direct
	     * renderers publish the complete logical work record for this
	     * cross-mode proof.  GPU-only tiers remain on their ordinary
	     * mode-local calibration until they publish the same record. */
	    const bool exactPresentedCost = measuredCadRenderCost;
	    const bool exactScaleQualityFrame =
		interactiveSample &&
		this->d->lodViewDemandPolicy.scaleChangingInteraction(
		    interactiveSample) &&
		exactPresentedCost;
	    int presentedMaximum = calibrationState ?
		calibrationState->maximumActiveProgressiveCut() : -1;
	    if (this->d->lodInteractiveProgressiveCeiling >= 0 &&
		presentedMaximum >= 0)
		presentedMaximum = std::min(presentedMaximum,
		    this->d->lodInteractiveProgressiveCeiling);
	    /* Stable throughput targets the configured quiet cadence, but the
	     * restoration witness has a different contract: it is the richest
	     * exact scale-quality frame which completed under the explicit 100 ms
	     * hard bound.  Conflating the two remembered an earlier 60 Hz motion
	     * cut and discarded a visibly richer, deadline-safe frame as soon as
	     * wheel input became quiet.  The richer cut is already resident and
	     * was actually presented; retaining it does not teach the 20 Hz
	     * allocator an inflated throughput. */
	    if (exactScaleQualityFrame && elapsed > 0 &&
		elapsed <= BObolLodViewDemandPolicy::
		    qualityFrameDurationNanoseconds())
		this->d->lodPresentationPolicy.noteStableQuality(
		    this->d->lodTargetPixelError,
		    presentedMaximum,
		    this->d->lodPresentationPointProxyPixelThreshold,
		    controller_lod_presentation_population(calibrationState,
			this->d->lodViewQualityDomainRevision),
		    this->d->lodViewRevision, calibrationCost);
	    const BObolLodQualityPolicy::StableCapacityEvidence
		stableEvidence = BObolLodQualityPolicy::stableCapacityEvidence(
		    calibrationCost, elapsed,
		    this->d->quietAllocationTargetFps(),
		    exactScaleQualityFrame);
		    if (stableEvidence.proven) {
			this->d->lodStableCalibratedRenderCostPerSecond = std::max(
			    this->d->lodStableCalibratedRenderCostPerSecond,
			    stableEvidence.renderCostPerSecond);
		    }
	}
    }

    /* A terminal budget pass explicitly requests an unchanged replay for its
     * late-headroom witness.  A retained VBO may still perform its one-time
     * upload or command preparation on that first replay; permit a strictly
     * bounded successor replay in that case, then consume only the first
     * reusable frame.  Opportunistic reuse of a later selection/HUD repaint
     * made a view report STABLE and then restart LoD work under unrelated user
     * interaction. */
    if (this->d->lodHeadroomPolicy.retryPending()) {
	const bool stableContext =
	    !this->d->lodInteractive && !this->d->lodGestureActive &&
	    this->d->lodBudgetPolicy.stableBudgetLimited() &&
	    !this->d->lodSubmissionPending &&
	    !this->d->lodBudgetPolicy.rescanAfterFrame() &&
	    !this->d->lodFrameObligation.pending() &&
	    this->d->quietAllocationTargetFps() > 0.0f &&
	    this->d->lodStableCalibratedRenderCostPerSecond > 0.0L &&
	    this->d->lodBudgetPolicy.currentBudget() != SIZE_MAX;
	const long double demonstratedBudget = stableContext ?
	    this->d->lodStableCalibratedRenderCostPerSecond /
		static_cast<long double>(
		    this->d->quietAllocationTargetFps()) : 0.0L;
	const uint64_t stableTargetNanoseconds = stableContext ?
	    static_cast<uint64_t>(1000000000.0L /
		static_cast<long double>(
		    this->d->quietAllocationTargetFps())) : 0;
	const bool matchingHeadroomWitness =
	    this->d->lodHeadroomPolicy.pendingMatches(
		this->d->lodViewRevision, this->d->lodPolicyRevision,
		activeCalibrationCost);
	const size_t discreteTrialExcess =
	    BObolLodQualityPolicy::discreteTrialOverBudgetAllowance(
		activeCalibrationCost,
		this->d->lodBudgetPolicy.currentBudget());
	const bool stableDiscreteTrial =
	    matchingHeadroomWitness && stableContext &&
	    measuredCadRenderCost && reusableCadWorkSample &&
	    !transientStablePublication && !transientStablePreparation &&
	    calibrationCost == activeCalibrationCost &&
	    stableTargetNanoseconds > 0 &&
	    elapsed <= stableTargetNanoseconds && discreteTrialExcess > 0;
	const bool transientHeadroomReplay =
	    matchingHeadroomWitness && stableContext &&
	    measuredCadRenderCost &&
	    (transientStablePublication || transientStablePreparation);
	const bool deferredHeadroomReplay = transientHeadroomReplay &&
	    this->d->lodHeadroomPolicy.deferTransientReplay(
		this->d->lodViewRevision, this->d->lodPolicyRevision,
		activeCalibrationCost);
	const bool reusableHeadroom = !deferredHeadroomReplay &&
	    this->d->lodHeadroomPolicy.consumeRetry(
		this->d->lodViewRevision, this->d->lodPolicyRevision,
		this->d->lodBudgetPolicy.currentBudget(), activeCalibrationCost,
		demonstratedBudget,
		calibrationElapsed, stableTargetNanoseconds,
		stableContext && reusableCadWorkSample &&
		    !transientStablePublication &&
		    !transientStablePreparation);
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
		this->d->lodViewRevision.value()))
	    bu_log("BObol LoD headroom complete stable=%d matching=%d "
		   "measured=%d reusable=%d transient_publication=%d "
		   "transient_preparation=%d deferred=%d admitted=%d "
		   "discrete=%d active_cost=%zu calibration_cost=%zu "
		   "budget=%zu demonstrated=%.3f elapsed_ms=%.3f "
		   "target_ms=%.3f view=%llu policy=%llu\n",
		   stableContext ? 1 : 0,
		   matchingHeadroomWitness ? 1 : 0,
		   measuredCadRenderCost ? 1 : 0,
		   reusableCadWorkSample ? 1 : 0,
		   transientStablePublication ? 1 : 0,
		   transientStablePreparation ? 1 : 0,
		   deferredHeadroomReplay ? 1 : 0,
		   reusableHeadroom ? 1 : 0,
		   stableDiscreteTrial ? 1 : 0,
		   activeCalibrationCost, calibrationCost,
		   this->d->lodBudgetPolicy.currentBudget(),
		   static_cast<double>(demonstratedBudget),
		   static_cast<double>(calibrationElapsed) / 1000000.0,
		   static_cast<double>(stableTargetNanoseconds) / 1000000.0,
		   (unsigned long long)this->d->lodViewRevision.value(),
		   (unsigned long long)this->d->lodPolicyRevision.value());
	if (deferredHeadroomReplay) {
	    this->markProgressiveWorkPending();
	    this->requestRender("lod-calibrated-headroom-replay");
	} else if (reusableHeadroom || stableDiscreteTrial) {
	    this->d->lodBudgetPolicy.clearBudgetLimit();
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionPending = TRUE;
	    this->d->lodSubmissionRescanPending = FALSE;
	    this->d->lodBudgetPolicy.resetPass();
	    /* The replay proved a different sustainable allowance.  Recompute one
	     * complete screen-importance allocation with that new capacity even
	     * when the old active population is already below it.  An ordinary
	     * first-come pass can strand the unused remainder on a discrete PoP
	     * population and makes cold and warm caches converge differently. */
	    this->d->lodBudgetPolicy.requestRetainedReallocation(false);
	    this->d->lodDiscretePopulationTrialAvailable =
		stableDiscreteTrial ? TRUE : FALSE;
	    this->markProgressiveWorkPending();
	    this->requestRender(stableDiscreteTrial ?
		"lod-discrete-headroom" : "lod-calibrated-headroom");
	}
    }

    /*
     * Consume a quiet small-part aggregation probe without invoking another
     * O(N) provider/admission scan.  If the measured reusable presentation
     * still misses the stable target, ratchet the camera-local threshold and
     * request exactly one successor frame.  The immutable per-part geometry
     * and desired PoP cuts are untouched throughout.
    */
    /* Point aggregation is implemented by SoCADAssembly's occurrence batch.
     * Ordinary retained VBO payloads have no point-cut population to
     * calibrate.  Waiting for one on that route creates an unwitnessable
     * progressive barrier even though all renderer/service work is idle.
     *
     * An assembly may also remain retained after every occurrence under its
     * drawn root has been hidden (for example, `erase all/r.stl`).  The latest
     * complete coverage census is the authority in that case: a zero visible
     * population has no point threshold to measure.  Retire the old bracket
     * and restore the one-pixel default while the scene is empty.  Otherwise
     * an interrupted pre-erase software frame can leave the controller asking
     * an exact zero-work frame for a CAD work record which cannot exist. */
    SbBool haveCurrentCadSourceTargets = FALSE;
    if (!haveCadPresentationAssemblies &&
	this->d->lodStablePointProxyCalibrationPending) {
	const std::vector<SoBRLDatabaseSource *> currentSources =
	    controller_render_database_source_roots(this);
	for (const SoBRLDatabaseSource *source : currentSources) {
	    if (source && source->getDisplayMeshLodRequestCount() > 0) {
		haveCurrentCadSourceTargets = TRUE;
		break;
	    }
	}
    }
    /* Candidate counts are invalidated at the start of a source/selection
     * transaction and are repopulated only by a complete projected census.
     * They may therefore be zero while the retained assembly still contains
     * thousands of visible instances.  Treating that transient telemetry gap
     * as an empty scene resets a measured point cut without restoring the
     * occurrence cuts hidden behind it; a selection/erase repaint can then
     * report ready at the preceding coarse rotation allocation.
     *
     * SoCADAssembly::instanceCount is the actual current presentation
     * population and hasCadPresentationAssemblies() samples it directly.  An
     * erase/expand transaction may briefly replace that instance array, so
     * the source-backed display request inventory is the second half of the
     * proof.  Only a zero on both sides may retire the bracket as an empty
     * scene. */
    const SbBool noVisibleCadPopulation =
	!haveCadPresentationAssemblies &&
	!haveCurrentCadSourceTargets &&
	!this->d->lodSubmissionPending ? TRUE : FALSE;
    if (this->d->lodStablePointProxyCalibrationPending &&
	noVisibleCadPopulation) {
	this->d->lodStablePointProxyCalibrationPending = FALSE;
	this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
	this->d->lodPointProxyCalibrationPolicy.reset();
	BObolViewLodState *presentationState = this->d->viewAttachment ?
	    this->d->viewAttachment->getViewLodState() : NULL;
	if (presentationState)
	    presentationState->setCadPresentationPointProxyPixelThreshold(1.0f);
	(void)this->d->lodBudgetPolicy.
	    confirmRetainedRecoveryPresentation(true);
    }
    if (this->d->lodStablePointProxyCalibrationPending &&
	!haveCadPresentationAssemblies &&
	!haveCurrentCadSourceTargets &&
	!this->d->lodSubmissionPending) {
	this->d->lodStablePointProxyCalibrationPending = FALSE;
	(void)this->d->lodBudgetPolicy.
	    confirmRetainedRecoveryPresentation(true);
    }
    if (this->d->lodStablePointProxyCalibrationPending &&
	!haveCadPresentationAssemblies && haveCurrentCadSourceTargets) {
	/* A path erase/expand transaction replaces the assembly occurrence table
	 * atomically.  Its source inventory already proves that a successor CAD
	 * presentation is expected, but this completed frame cannot calibrate the
	 * temporarily empty renderer record.  Preserve the obligation and wake the
	 * source producer which installs that presentation.  Requesting another
	 * empty render here while calibration also paused source admission formed a
	 * closed BALANCING/REFINING loop on 50k scenes. */
	if (!this->d->lodSubmissionPending) {
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionRescanPending = FALSE;
	    this->d->lodSubmissionPending = TRUE;
	    this->d->lodBudgetPolicy.resetPass();
	}
	this->markProgressiveWorkPending();
    }
    if (this->d->lodStablePointProxyCalibrationPending &&
	haveCadPresentationAssemblies) {
	this->d->lodStablePointProxyCalibrationPending = FALSE;
	const long double stableTargetNanoseconds =
	    this->d->quietAllocationTargetFps() > 0.0f ?
		1000000000.0L /
		    static_cast<long double>(
			this->d->quietAllocationTargetFps()) : 0.0L;
	const SbBool reusableStableSample =
	    !this->d->lodInteractive &&
	    !this->d->lodGestureActive &&
	    reusableCadWorkSample &&
	    !transientStablePublication &&
	    !transientStablePreparation &&
	    stableTargetNanoseconds > 0.0L &&
	    pointCalibrationElapsed > 0;
	/* A pending aggregation measurement is not a repaint timer.  While source
	 * submission or immutable results are still advancing, their publication
	 * edges will eventually provide a current frame and the progressive pump
	 * remains live independently.  Repainting the unchanged many-leaf scene on
	 * every pump starves that producer work on software renderers and makes
	 * post-rotation restoration look like a long refinement replay. */
	const bool pointCalibrationServiceWork =
	    this->hasPendingLodResults() ||
	    (this->d->lodService &&
	     (this->d->lodService->activeTaskCountForGeneration(
		 this->d->lodActiveGeneration) > 0 ||
	      this->d->lodService->queuedResultCountForGeneration(
		 this->d->lodActiveGeneration) > 0));
	const SbBool pointCalibrationProducerWork =
	    BObolLodPointProxyCalibrationPolicy::
		producerOwnsCalibrationFrame(
		    this->d->lodSubmissionPending != FALSE,
		    false,
		    this->d->progressiveProviderPendingCount > 0,
		    pointCalibrationServiceWork,
		    this->d->lodPublicationPolicy.pending()) ? TRUE : FALSE;
	const size_t minimumCalibrationCost = calibrationState ?
	    calibrationState->minimumActiveRenderCost() :
	    activeCalibrationCost;
	const uint64_t preferredStableNanoseconds =
	    this->d->lodStableTargetFps > 0.0f ?
	    static_cast<uint64_t>(1000000000.0L /
		static_cast<long double>(this->d->lodStableTargetFps)) : 0;
	const SbBool exactReusablePointFrame =
	    reusableStableSample && calibrationState &&
	    calibrationState->lastCadPresentationFrameExact() &&
	    exactOccurrenceClassification &&
	    unresolvedStructuralOccurrences == 0 ? TRUE : FALSE;
	const SbBool startStaticOnePixelTrial =
	    BObolLodStaticQualityPolicy::
		startOnePixelTrialFromSettledPointFrame(
		    this->d->lodInteractive || this->d->lodGestureActive,
		    exactReusablePointFrame != FALSE,
		    pointCalibrationProducerWork != FALSE,
		    this->d->lodResourcePolicy.anyPressure(),
		    structuralFallbackOccurrences != 0,
		    this->d->lodStaticOverscanRejected != FALSE,
		    this->d->lodPresentationPointProxyPixelThreshold,
		    pointCalibrationElapsed, preferredStableNanoseconds,
		    this->d->stablePresentationFrameDeadlineNanoseconds) ?
		TRUE : FALSE;
	const SbBool reducibleProgressiveDetail =
	    activeCalibrationCost > minimumCalibrationCost ? TRUE : FALSE;
	const SbBool stableSampleOverloaded =
	    reusableStableSample &&
	    static_cast<long double>(pointCalibrationElapsed) >
		stableTargetNanoseconds * 1.10L ? TRUE : FALSE;
	const SbBool coarsePointCut =
	    this->d->lodPresentationPointProxyPixelThreshold > 1.01f ?
		TRUE : FALSE;
	const SbBool acceptSettledOnePixelFrame =
	    BObolLodStaticQualityPolicy::acceptSettledOnePixelFrame(
		this->d->lodInteractive || this->d->lodGestureActive,
		exactReusablePointFrame != FALSE,
		pointCalibrationProducerWork != FALSE,
		this->d->lodResourcePolicy.anyPressure(),
		structuralFallbackOccurrences != 0,
		this->d->lodPresentationPolicy.handoffPending(),
		this->d->lodStaticOverscanActive != FALSE,
		this->d->lodStaticOverscanRejected != FALSE,
		this->d->lodInteractiveProgressiveCeiling,
		this->d->lodPresentationPointProxyPixelThreshold,
		reducibleProgressiveDetail != FALSE,
		pointCalibrationElapsed, preferredStableNanoseconds,
		this->d->stablePresentationFrameDeadlineNanoseconds) ?
		TRUE : FALSE;
	const BObolRetainedAllocationResult &pointAllocation =
	    this->d->lodRetainedAllocationCertificate;
	const bool triangleAllocationSaturated = calibrationState &&
	    this->d->lodPointProxyTriangleRecoverySaturatedPlanSerial != 0 &&
	    pointAllocation.allocationPlanSerial ==
		this->d->lodPointProxyTriangleRecoverySaturatedPlanSerial &&
	    pointAllocation.viewRevision == this->d->lodViewRevision.value() &&
	    pointAllocation.policyRevision ==
		this->d->lodPolicyRevision.value() &&
	    calibrationState->cadAllocationPlanCoversCurrentPopulation(
		pointAllocation.allocationPlanSerial,
		pointAllocation.viewRevision, pointAllocation.policyRevision,
		pointAllocation.fixedCadPresentationCost);
	SbBool scheduledTriangleRecovery = FALSE;
	SbBool restoredOnePixelCut =
	    startStaticOnePixelTrial || acceptSettledOnePixelFrame;

	if (acceptSettledOnePixelFrame) {
	    /* The current framebuffer is already the exact one-pixel static
	     * witness.  Preserve it while one occurrence-local minimax pass spends
	     * the capacity extrapolated to the hard deadline.  No point threshold,
	     * renderer ceiling, resident suffix, or visible cut is reset. */
	    this->d->lodStaticOverscanActive = TRUE;
	    this->d->lodStaticOverscanLeapAvailable = TRUE;
	    this->d->lodDiscretePopulationTrialAvailable = TRUE;
	    const size_t staticQualityBudget =
		BObolLodQualityPolicy::staticPresentationRenderCostLimit(
		    activeCalibrationCost, pointCalibrationElapsed,
		    this->d->stablePresentationFrameDeadlineNanoseconds);
	    this->d->lodHeadroomPolicy.cancelRetry();
	    this->d->lodBudgetPolicy.clearBudgetLimit();
	    this->d->lodBudgetPolicy.resetProbeSeries();
	    this->d->lodBudgetPolicy.resetOverloadRecovery();
	    this->d->lodBudgetPolicy.resetPass();
	    if (staticQualityBudget > 0)
		this->d->lodBudgetPolicy.raiseCurrentBudget(
		    staticQualityBudget);
	    this->d->lodBudgetPolicy.requestRetainedReallocation(false);
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionRescanPending = FALSE;
	    this->d->lodSubmissionPending = TRUE;
	    this->d->lodPassAdmittedWork = FALSE;
	    this->d->lodRetainedRefinementPending = FALSE;
	    this->d->lodRetainedResidencyPending = FALSE;
	    this->d->lodRetainedRefinementCutAdvanced = FALSE;
	    this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	    this->markProgressiveWorkPending();
	}

	if (startStaticOnePixelTrial) {
	    /* Preserve every retained occurrence cut and perform one bounded
	     * presentation-only population trial.  The endpoint deadline keeps the
	     * preceding complete framebuffer visible if this richer classification
	     * is too expensive; its miss establishes the unsafe side of the point
	     * bracket without a coarse triangle-allocation round trip. */
	    this->d->lodStaticOverscanActive = TRUE;
	    this->d->lodStaticOverscanLeapAvailable = TRUE;
	    this->d->lodDiscretePopulationTrialAvailable = TRUE;
	    this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
	    this->d->lodPointProxyCalibrationPolicy.reset();
	    BObolViewLodState *presentationState =
		this->d->viewAttachment->getViewLodState();
	    if (presentationState)
		presentationState->setCadPresentationPointProxyPixelThreshold(
		    1.0f);
	    this->d->lodStablePointProxyCalibrationPending = TRUE;
	    this->requestRender("lod-static-point-one-pixel");
	}

	/* Resolve the two quality dimensions from measured work.  A multi-pixel
	 * point cut represents only occurrences below that screen-space size; it
	 * is not permission to discard detail from large arrays, wheels, blades,
	 * tails, or other prominent parts.  The former coarsePointCut condition
	 * collapsed every reducible prefix even after a safe aggregated frame,
	 * restored the one-pixel population, and then repeated forever.  Recover
	 * triangle capacity only from an exact overloaded frame.  Otherwise let
	 * the point bracket converge while the importance allocator spends the
	 * retained triangle budget on visible multi-pixel geometry.
	 */
	if (!restoredOnePixelCut && reusableStableSample &&
	    !unresolvedStructuralOccurrences &&
	    BObolLodPointProxyCalibrationPolicy::shouldRecoverTriangleDetail(
		reducibleProgressiveDetail != FALSE,
		stableSampleOverloaded != FALSE,
		coarsePointCut != FALSE,
		this->d->lodBudgetPolicy.retainedQualityFloorActive() ||
		    this->d->lodBudgetPolicy.retainedQualityFloorRejected() ||
		    triangleAllocationSaturated)) {
	    /* recoveryBudget is consumed by the retained allocator, so scale its
	     * own active-cost population by the measured deadline ratio.  Backend
	     * submitted-work counts are not dimensionally interchangeable with
	     * this budget (notably for OSMesa's expanded position stream). */
	    long double affordable =
		static_cast<long double>(activeCalibrationCost) *
		stableTargetNanoseconds * 0.80L /
		static_cast<long double>(pointCalibrationElapsed);
	    if (!std::isfinite(affordable) || affordable <= 0.0L)
		affordable = static_cast<long double>(minimumCalibrationCost);
	    size_t recoveryBudget = affordable >=
		    static_cast<long double>(SIZE_MAX) ? SIZE_MAX :
		static_cast<size_t>(affordable);
	    recoveryBudget = std::max(
		minimumCalibrationCost, recoveryBudget);
	    if (recoveryBudget < activeCalibrationCost) {
		this->d->lodBudgetPolicy.requestRetainedRecovery(
		    recoveryBudget);
		this->d->lodBudgetPolicy.resetProbeSeries();
		this->d->lodBudgetPolicy.resetOverloadRecovery();
		this->d->lodBudgetPolicy.resetPass();
		this->d->lodHeadroomPolicy.cancelRetry();
		/* The handoff's cost floor was proven with the coarse point cut
		 * which hid this population.  It is not evidence for the one-pixel
		 * hierarchy and would immediately re-admit the detail being
		 * compacted here.  Retain the renderer ceiling itself until the
		 * coherent recovery frame below, but retire the stale proof. */
		this->d->lodPresentationPolicy.cancelHandoff();
		this->d->lodSubmissionSourceIndex = 0;
		this->d->lodSubmissionEntryOffset = 0;
		this->d->clearLodSubmissionPlan();
		this->d->lodSubmissionRescanPending = FALSE;
		this->d->lodSubmissionPending = TRUE;
		this->d->lodRetainedRefinementPending = FALSE;
		this->d->lodRetainedRefinementCutAdvanced = FALSE;
		this->d->lodRetainedRefinementBudgetBlocked = FALSE;
		this->d->lodPassAdmittedWork = FALSE;
		this->d->lodPointProxyTriangleRecoveryPending = TRUE;
		this->markProgressiveWorkPending();
		scheduledTriangleRecovery = TRUE;
	    } else if (coarsePointCut &&
		structuralFallbackOccurrences == 0) {
		/* The complete retained cut already fits the demonstrated capacity.
		 * Test the one-pixel population directly; no provider or admission
		 * traversal is necessary. */
		this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
		this->d->lodPointProxyCalibrationPolicy.reset();
		BObolViewLodState *presentationState =
		    this->d->viewAttachment->getViewLodState();
		if (presentationState)
		    presentationState->
			setCadPresentationPointProxyPixelThreshold(1.0f);
		this->d->lodStablePointProxyCalibrationPending = TRUE;
		this->requestRender("lod-stable-point-restore");
		restoredOnePixelCut = TRUE;
	    }
	}

	const BObolLodPointProxyCalibrationPolicy::Decision calibration =
	    (!scheduledTriangleRecovery && !restoredOnePixelCut) ?
	    this->d->lodPointProxyCalibrationPolicy.completed(
		this->d->lodPresentationPointProxyPixelThreshold,
		pointCalibrationElapsed,
		this->d->quietAllocationTargetFps(),
		reusableStableSample != FALSE,
		structuralFallbackOccurrences) :
	    BObolLodPointProxyCalibrationPolicy::Decision();
	if (!scheduledTriangleRecovery && !restoredOnePixelCut &&
	    calibration.changed) {
	    this->d->lodPresentationPointProxyPixelThreshold =
		calibration.threshold;
	    BObolViewLodState *presentationState =
		this->d->viewAttachment->getViewLodState();
	    if (presentationState)
		presentationState->
		    setCadPresentationPointProxyPixelThreshold(
			calibration.threshold);
	    this->d->lodStablePointProxyCalibrationPending = TRUE;
	    this->requestRender("lod-stable-point-calibration");
	} else if (!scheduledTriangleRecovery && !restoredOnePixelCut &&
	    !reusableStableSample &&
	    !this->d->lodInteractive &&
	    !this->d->lodGestureActive &&
	    haveCadPresentationAssemblies &&
	    stableTargetNanoseconds > 0.0L &&
	    this->d->lodPointProxyCalibrationPolicy.
		requiresReusableConfirmation(
		    this->d->lodPresentationPointProxyPixelThreshold,
		    structuralFallbackOccurrences)) {
	    /*
	     * Applying a different point threshold changes the retained
	     * occurrence classification.  Its first frame may rebuild/patch the
	     * prepared plan, publish worker results, or follow an interrupted
	     * traversal and therefore lack a reusable work record.  None of those
	     * cases is evidence that calibration converged.  Preserve the
	     * obligation until an unchanged replay measures the cut it actually
	     * drew.
	     *
	     * Result streaming may extend this handoff across several frames.
	     * Those frames already represent real pending work; retaining the
	     * idempotent frame request neither exposes a richer cut prematurely nor
	     * permits the view to report stable at an unmeasured emergency cut.
	     */
	    this->d->lodStablePointProxyCalibrationPending = TRUE;
	    if (pointCalibrationProducerWork) {
		this->markProgressiveWorkPending();
	    } else {
		this->requestRender("lod-stable-point-replay");
	    }
	}
	/* A triangle-prefix recovery which temporarily used a coarser aggregate
	 * point cut reaches this branch twice: once to request the one-pixel frame,
	 * and once after that frame completes.  The latter exact presentation is
	 * the missing transition which retires the measured recovery ceiling. */
	if (!scheduledTriangleRecovery && !restoredOnePixelCut &&
	    exactRetainedPopulation &&
	    this->d->lodPresentationPointProxyPixelThreshold <= 1.01f)
	    this->resumeLodAfterRetainedRecovery();
    }

	/* A retained-prefix recovery is complete only after its selected cut has
	 * been presented and no follow-up admission/handoff barrier remains.
	 * Restore one pixel only for a fully realized mesh population.  If
	 * structural records remain, retain their proven aggregate cut: lowering
	 * it would expose boxes rather than detail. */
    if (this->d->lodPointProxyTriangleRecoveryPending &&
	!this->d->lodInteractive && !this->d->lodGestureActive &&
	!this->d->lodSubmissionPending &&
	!this->d->lodSubmissionRescanPending &&
	!this->d->lodRetainedRefinementPending &&
	!this->d->lodBudgetPolicy.rescanAfterFrame() &&
	!this->d->lodPresentationPolicy.handoffPending()) {
	BObolViewLodState *presentationState =
	    this->d->viewAttachment->getViewLodState();
	const BObolRetainedAllocationResult &recoveryAllocation =
	    this->d->lodRetainedAllocationCertificate;
	if (presentationState && recoveryAllocation.allocationPlanSerial != 0 &&
	    recoveryAllocation.viewRevision ==
		this->d->lodViewRevision.value() &&
	    recoveryAllocation.policyRevision ==
		this->d->lodPolicyRevision.value() &&
	    presentationState->cadAllocationPlanCoversCurrentPopulation(
		recoveryAllocation.allocationPlanSerial,
		recoveryAllocation.viewRevision,
		recoveryAllocation.policyRevision,
		recoveryAllocation.fixedCadPresentationCost))
	    this->d->lodPointProxyTriangleRecoverySaturatedPlanSerial =
		recoveryAllocation.allocationPlanSerial;
	this->d->lodPointProxyTriangleRecoveryPending = FALSE;
	const float recoveredPointThreshold =
	    BObolLodPointProxyCalibrationPolicy::
		triangleRecoveryPointThreshold(
		    this->d->lodPresentationPointProxyPixelThreshold,
		    exactStructuralProjection != FALSE,
		    structuralFallbackOccurrences);
	const SbBool pointCutChanged = std::fabs(
	    recoveredPointThreshold -
	    this->d->lodPresentationPointProxyPixelThreshold) > 0.01f ?
		TRUE : FALSE;
	this->d->lodInteractiveProgressiveCeiling = -1;
	this->d->lodPresentationPointProxyPixelThreshold =
	    recoveredPointThreshold;
	this->d->lodPointProxyCalibrationPolicy.reset();
	if (presentationState)
	    presentationState->setCadPresentationProgressiveCutCeiling(-1);
	if (presentationState)
	    presentationState->setCadPresentationPointProxyPixelThreshold(
		recoveredPointThreshold);
	/* Point calibration measures a changed occurrence population.  A retained
	 * triangle-prefix recovery commonly starts and ends at the one-pixel cut;
	 * arming the point latch in that case asks an unchanged ordinary OSMesa
	 * VBO frame for CAD-batch evidence it cannot publish and can strand an
	 * otherwise idle view forever. */
	this->d->lodStablePointProxyCalibrationPending =
	    pointCutChanged && presentationState &&
	    presentationState->hasCadPresentationAssemblies() ? TRUE : FALSE;
	if (this->d->lodStablePointProxyCalibrationPending) {
	    this->requestRender("lod-stable-point-restore");
	} else {
	    /* The recovery ceiling prevents the cheaper preparation frame from
	     * immediately re-admitting the population which just missed its
	     * deadline.  It is a one-witness guard, not a permanent fidelity cap.
	     * Once the coherent recovered point/mesh cut has actually rendered, its
	     * measured throughput must be allowed to grow the budget again in
	     * bounded, screen-priority waves.  Keeping the ceiling indefinitely
	     * made a transient 84 ms Hubble frame settle forever at the 35k-face
	     * minimum even though subsequent 16 ms frames proved ample headroom. */
	    this->resumeLodAfterRetainedRecovery();
	}
    }

    /*
     * A cold or newly exposed many-leaf stream can spend seconds completing
     * its coverage pass.  Waiting until that pass ends before applying the
     * renderer's measured small-part cut leaves every intervening software
     * frame at 100+ ms and makes console/input handling feel stalled.
     *
     * Result publication time is deliberately excluded from sustainable
     * throughput calibration above, but it is still real user-visible frame
     * latency.  Under severe pressure, raise the presentation-only
     * aggregation cut for the next frame while coverage continues.  The
     * final quiet probe later calibrates it precisely and no retained mesh or
     * cache state is discarded.
     */
    const long double ongoingStableTargetNanoseconds =
	this->d->quietAllocationTargetFps() > 0.0f ?
	    1000000000.0L /
		static_cast<long double>(
		    this->d->quietAllocationTargetFps()) : 0.0L;
    const bool ongoingServiceWork = this->d->lodService &&
	(this->d->lodService->activeTaskCountForGeneration(
	     this->d->lodActiveGeneration) > 0 ||
	 this->d->lodService->queuedResultCountForGeneration(
	     this->d->lodActiveGeneration) > 0);
    const SbBool ongoingStableWork =
	BObolLodPointProxyCalibrationPolicy::
	    streamingPopulationWorkPending(
		this->d->lodSubmissionPending != FALSE,
		this->hasPendingLodResults(), ongoingServiceWork,
		this->d->progressiveProviderPendingCount) ? TRUE : FALSE;
    if (!admissionPointProxyFrameWasPending &&
	!this->d->lodAdmissionPointProxyFramePending &&
	!this->d->lodStablePointProxyCalibrationPending &&
	haveCadPresentationAssemblies &&
	this->d->pointProxyAggregationApplicable() &&
	!this->d->lodInteractive &&
	!this->d->lodGestureActive &&
	this->d->lodAutoSubmit &&
	this->d->lodDiscoveryPointProxyPixelThreshold <= 1.01f &&
	ongoingStableWork &&
	measuredCadPresentation &&
	ongoingStableTargetNanoseconds > 0.0L &&
	static_cast<long double>(elapsed) >
	    ongoingStableTargetNanoseconds * 1.50L) {
	const BObolLodPointProxyCalibrationPolicy::Decision pressure =
	    this->d->lodPointProxyCalibrationPolicy.interrupted(
		this->d->lodPresentationPointProxyPixelThreshold,
		elapsed, this->d->quietAllocationTargetFps());
	if (pressure.changed) {
	    this->d->lodPresentationPointProxyPixelThreshold =
		pressure.threshold;
	    BObolViewLodState *presentationState =
		this->d->viewAttachment->getViewLodState();
	    if (presentationState)
		presentationState->
		    setCadPresentationPointProxyPixelThreshold(
			pressure.threshold);
	    /* The changed SoCADAssembly occurrence population needs one exact
	     * successor frame before the bracket can be adjusted again. */
	    this->d->lodStablePointProxyCalibrationPending =
		this->d->lodPresentationPointProxyPixelThreshold > 1.01f ?
		    TRUE : FALSE;
	    /* This pressure step is deliberately conservative because its frame
	     * includes streaming/publication latency.  The unchanged successor
	     * frame will walk back toward the finest sustainable presentation when
	     * it demonstrates headroom. */
	    this->requestRender("lod-stream-point-pressure");
	}
    }

    /* A pressure step may have changed the aggregate threshold while source
     * or result streaming still owned the next frame.  If that producer goes
     * quiet with structural boxes in the exact completed population, the
     * calibration latch no longer has a frame it can legitimately measure.
     * Hand the frontier to armStableLodHeadroomProbeIfReady() below instead
     * of leaving both admission and repair permanently paused. */
    const bool pointCalibrationFutureProducer =
	ongoingStableWork || this->d->lodPublicationPolicy.pending();
    if (BObolLodPointProxyCalibrationPolicy::
	    calibrationYieldsToStructuralRepair(
		this->d->lodStablePointProxyCalibrationPending != FALSE,
		exactOccurrenceClassification != FALSE,
		unresolvedStructuralOccurrences,
		pointCalibrationFutureProducer)) {
	this->d->lodStablePointProxyCalibrationPending = FALSE;
	this->markProgressiveWorkPending();
    }

    /* Point aggregation and the stable-presentation handoff can consume the
     * frame which followed the terminal budget pass.  Revisit the explicit
     * late-headroom witness immediately after those barriers have either
     * converged or re-armed themselves.  This preserves an already resident
     * rich zoom cut and lets a quiet view continue from it instead of
     * compacting to the conservative seed budget. */
    this->armStableLodHeadroomProbeIfReady();

    /* A zoom-quality probe is allowed to spend more than an ordinary motion
     * frame, but never more than the 10 Hz hard responsiveness floor.  If a
     * newly exposed PoP cut misses that floor, correct the next frame with
     * the renderer-wide prefix ceiling.  This is an O(1), reversible
     * presentation change: the richer occurrence cut and its resident bytes
     * stay intact, so continued zoom can refine from them and quiet handoff
     * does not have to walk back up from the minimum mesh.
     *
     * Cut ordinals do not encode a population ratio.  Back off exactly one
     * producer-admissible cut after a proven miss and let the following
     * completed frame decide whether another step is necessary. */
    if ((this->d->lodInteractive || this->d->lodGestureActive) &&
	this->d->lodViewDemandPolicy.interactionScaleChanged() &&
	this->d->lodViewDemandPolicy.qualityProbeActive() &&
	elapsed >
	    BObolLodViewDemandPolicy::qualityFrameDurationNanoseconds()) {
	BObolViewLodState *presentationState =
	    this->d->viewAttachment->getViewLodState();
	const int activeMaximum = presentationState ?
	    presentationState->maximumActiveProgressiveCut() : -1;
	if (activeMaximum > 0) {
	    const int presentedMaximum =
		this->d->lodInteractiveProgressiveCeiling >= 0 ?
		    std::min(activeMaximum,
			this->d->lodInteractiveProgressiveCeiling) :
		    activeMaximum;
	    const int correctedCeiling =
		std::max(0, presentedMaximum - 1);
	    if (correctedCeiling < presentedMaximum) {
		this->d->lodViewDemandPolicy.noteQualityMiss(
		    correctedCeiling,
		    measuredCadRenderCost ? presentedCadRenderCost :
			(activeCalibrationCost ? activeCalibrationCost :
			 presentedCadPrimitives));
		this->d->lodInteractiveProgressiveCeiling =
		    correctedCeiling;
		presentationState->setCadPresentationProgressiveCutCeiling(
		    correctedCeiling);
		this->requestRender("lod-scale-quality-pressure");
	    }
	}
    }

    /* A budget-limited pass may have admitted a bounded first wave and then
     * scanned the remaining boxes without admitting them.  Re-plan from the
     * highest screen-value occurrence after that wave has actually rendered
     * and supplied a throughput sample.  Replanning before the frame would
     * repeatedly spend the same stale allowance. */
    if (this->d->lodBudgetPolicy.rescanAfterFrame()) {
	/*
	 * Calibration probes measure an unchanged retained cut.  They do not
	 * need an intervening O(N) occurrence-planning pass: that pass cannot
	 * admit anything until the samples have changed the aggregate budget.
	 * Present the bounded probe series back-to-back, then scan the sparse
	 * unsatisfied frontier once using the resulting calibration.
	 */
	const BObolLodBudgetPolicy::CompletedFrameDecision calibration =
	    this->d->lodBudgetPolicy.completeCalibrationFrame();
	if (calibration.requestCalibrationFrame) {
	    this->requestRender("lod-budget-calibration");
	} else if (calibration.restartSubmission) {
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionPending = TRUE;
	    this->d->lodSubmissionRescanPending = FALSE;
	    /* A missing bounded coverage census retires its counter pass while
	     * waiting for this measured-frame admission decision.  The calibrated
	     * retry is another complete current-view census, not merely a quality
	     * pass: reactivate its counters so a cold stream which has since
	     * received every minimum prefix can publish coverageComplete(). */
	    if (this->d->lodCoveragePolicy.required() &&
		!this->d->lodCoveragePolicy.coverageComplete() &&
		!this->d->lodCoveragePolicy.active())
		this->d->lodCoveragePolicy.activate(false);
	    if (this->d->lodCoveragePolicy.active())
		this->d->lodCoveragePolicy.clearPassCounters();
	    this->markProgressiveWorkPending();
	}
    }
    /* Once the scale-quality budget is active, every completed frame is a
	 * measurement of the current render ceiling even if a newer wheel event
	 * has already cleared the one-shot probe label.  Preserve that proof
	 * across the next camera epoch; otherwise the ordinary 60 Hz motion
	 * policy can immediately undo a 10 Hz quality cut which was just shown
	 * successfully. */
    int completedPresentedMaximum = calibrationState ?
	calibrationState->maximumActiveProgressiveCut() : -1;
    if (this->d->lodInteractiveProgressiveCeiling >= 0 &&
	completedPresentedMaximum >= 0)
	completedPresentedMaximum = std::min(completedPresentedMaximum,
	    this->d->lodInteractiveProgressiveCeiling);

    if (this->d->lodInteractiveProgressiveCeiling >= 0 &&
	completedPresentedMaximum >= 0 && calibrationState &&
	measuredCadRenderCost && presentedCadRenderCost > 0 &&
	(!this->d->stablePresentationFrameDeadlineNanoseconds ||
	 elapsed <= this->d->stablePresentationFrameDeadlineNanoseconds)) {
	/* The completed-work record belongs to the framebuffer just presented.
	 * A following submission may already have updated the mutable scene
	 * status, so it is not a reason to discard this direct capacity proof.
	 * Retain the richest completed cut for this unchanged camera/population;
	 * a later coarser recovery frame is not evidence that the richer frame
	 * stopped fitting. */
	const bool sameSafePresentation =
	    this->d->lodDeadlineSafeViewRevision ==
		this->d->lodViewRevision.value() &&
	    this->d->lodDeadlineSafeQualityDomainRevision ==
		this->d->lodViewQualityDomainRevision;
	this->d->lodDeadlineSafeProgressiveCeiling = sameSafePresentation ?
	    std::max(this->d->lodDeadlineSafeProgressiveCeiling,
		completedPresentedMaximum) : completedPresentedMaximum;
	this->d->lodDeadlineSafeViewRevision =
	    this->d->lodViewRevision.value();
	this->d->lodDeadlineSafeQualityDomainRevision =
	    this->d->lodViewQualityDomainRevision;
    }
    this->d->lodViewDemandPolicy.noteFramePresented(
	completedPresentedMaximum,
	measuredCadRenderCost != FALSE && presentedCadRenderCost > 0,
	elapsed);
    if (this->d->lodViewDemandPolicy.rearmAfterQualityFrame(
	    this->d->lodInteractive || this->d->lodGestureActive,
	    calibrationState ?
		calibrationState->maximumActiveProgressiveCut() : -1,
	    completedPresentedMaximum,
	    measuredCadRenderCost != FALSE && presentedCadRenderCost > 0,
	    elapsed))
	this->markProgressiveWorkPending();
    const size_t provenHandoffRenderCost =
	measuredCadRenderCost && presentedCadRenderCost > 0 &&
	calibrationState && calibrationState->lastCadPresentationFrameExact() &&
	(!this->d->stablePresentationFrameDeadlineNanoseconds ||
	 elapsed <= this->d->stablePresentationFrameDeadlineNanoseconds) ?
	    presentedCadRenderCost : 0;
    this->completePresentationBarrier(elapsed, provenHandoffRenderCost);

    /* A quiet single-occurrence handoff may have loaded several richer PoP
     * populations behind its last responsive renderer ceiling.  Expose them
     * only through completed-frame evidence.  The first hard-deadline miss
     * disables this sequence in notePresentationRenderInterrupted(), leaving
     * the last complete framebuffer and immutable resident suffix intact. */
    bool scheduledStaticPresentationStep = false;
    if (this->d->lodStaticOverscanActive &&
	!this->d->lodInteractive && !this->d->lodGestureActive &&
	this->d->lodInteractiveProgressiveCeiling >= 0 &&
	calibrationState && measuredCadRenderCost &&
	presentedCadRenderCost > 0 &&
	calibrationState->lastCadPresentationFrameExact() &&
	elapsed <= this->d->stablePresentationFrameDeadlineNanoseconds &&
	!this->d->lodPresentationPolicy.handoffPending() &&
	!this->d->lodSubmissionPending &&
	!this->d->lodBudgetPolicy.rescanAfterFrame() &&
	!this->d->lodFrameObligation.pending() &&
	!this->d->lodPublicationPolicy.pending()) {
	const int activeMaximum =
	    calibrationState->maximumActiveProgressiveCut();
	if (activeMaximum > this->d->lodInteractiveProgressiveCeiling) {
	    const size_t predictedCostLimit =
		BObolLodQualityPolicy::staticPresentationRenderCostLimit(
		    presentedCadRenderCost, elapsed,
		    this->d->stablePresentationFrameDeadlineNanoseconds);
	    const int boundedCeiling = calibrationState->
		singleCadProgressiveCutWithinRenderCost(predictedCostLimit);
	    if (boundedCeiling >= 0 && boundedCeiling <=
		    this->d->lodInteractiveProgressiveCeiling) {
		/* The exact cached next-cut cost does not fit the throughput
		 * demonstrated by this completed frame with a small jitter margin.
		 * The current framebuffer is therefore the terminal static-quality
		 * proof.  Reconcile hidden richer occurrence metadata to that cut
		 * under the same 10 Hz allowance before retiring the renderer-only
		 * ceiling; no speculative long frame or coarse restart is needed. */
		this->d->lodStaticOverscanLeapAvailable = FALSE;
		this->d->lodStaticOverscanRejected = TRUE;
		this->d->lodPresentationPolicy.armHandoff(
		    false, presentedCadRenderCost);
		this->d->lodSubmissionSourceIndex = 0;
		this->d->lodSubmissionEntryOffset = 0;
		this->d->clearLodSubmissionPlan();
		this->d->lodSubmissionPending = TRUE;
		this->d->lodSubmissionRescanPending = FALSE;
		this->d->lodBudgetPolicy.resetPass();
		this->markProgressiveWorkPending();
		scheduledStaticPresentationStep = true;
	    } else {
		int nextCeiling =
		    this->d->lodInteractiveProgressiveCeiling + 1;
		/* Avoid showing cheap intermediate prefixes after button-up, but
		 * never extrapolate beyond the exact draw-mode-aware cost which the
		 * completed frame predicts can retain five percent deadline
		 * headroom.  At most two ordinals are combined so a residual model
		 * mismatch still has one bounded midpoint fallback. */
		if (this->d->lodStaticOverscanLeapAvailable &&
		    boundedCeiling >=
			this->d->lodInteractiveProgressiveCeiling + 2)
		    nextCeiling =
			this->d->lodInteractiveProgressiveCeiling + 2;
		this->d->lodStaticOverscanLeapAvailable = FALSE;
		if (this->d->lodStaticOverscanActive) {
		    this->d->lodInteractiveProgressiveCeiling = nextCeiling;
		    calibrationState->setCadPresentationProgressiveCutCeiling(
			this->d->lodInteractiveProgressiveCeiling);
		    this->requestRender("lod-static-overscan-step");
		    scheduledStaticPresentationStep = true;
		}
	    }
	} else {
	    /* With a calibrated aggregate point population, catching the ceiling
	     * up to every currently active cut is not proof that all useful demand
	     * has been admitted.  Retain this exact safety ceiling while the
	     * headroom planner performs one occurrence-local allocation behind it;
	     * newly admitted cuts will then re-enter the staged sequence above.
	     * The allocation handoff removes the now-redundant ceiling when no
	     * richer local population can be selected. */
	    /* Once every occurrence-local cut is at or below the renderer limit,
	     * the limit is mathematically redundant.  Point aggregation is an
	     * independent, view-local batching policy and is not a reason to leave
	     * a scene-wide PoP ordinal installed in a stable state. */
	    this->d->lodInteractiveProgressiveCeiling = -1;
	    calibrationState->setCadPresentationProgressiveCutCeiling(-1);
	    /* Retiring a redundant renderer ceiling does not retire the accepted
	     * static-quality phase.  Its occurrence allocation was measured under
	     * the event-driven hard deadline and must remain the terminal policy
	     * until a camera, resource, or capacity edge invalidates it.  Clearing
	     * the phase here returned the next pass to the ordinary 20 Hz streaming
	     * budget; that coarsened the scene, re-armed the 10 Hz phase, and made
	     * large warm scenes alternate forever between the two allocations. */
	    this->d->lodStaticOverscanActive =
		BObolLodStaticQualityPolicy::retainAcceptedPhase(
		    this->d->lodStaticOverscanActive != FALSE,
		    this->d->lodStaticOverscanRejected != FALSE) ? TRUE : FALSE;
	    this->d->lodStaticOverscanLeapAvailable = FALSE;
	}
    }

    /* The frame-completion edge above is what retires the progressive-prefix
     * presentation barrier.  A terminal planning pass may already have tried
     * to arm headroom while that barrier was still active; if no other work
     * remains, there will be no later planning pass from which to try again.
     * Recheck after every completed-frame latch has been reconciled so a
     * cheap, newly presented population gets its one bounded richer probe
     * instead of being misreported as the stable final cut.  The headroom
     * policy keys that witness by view, policy, and active population, which
     * keeps an actually capacity-limited population from repainting in a
     * loop. */
    if (!scheduledStaticPresentationStep)
	this->armStableLodHeadroomProbeIfReady();
}

void
BObolViewController::completePresentationBarrier(uint64_t elapsedNanoseconds,
	size_t provenRenderCost)
{
    this->d->renderCompletionSerial++;
    if (this->d->renderCompletionSerial == 0)
	this->d->renderCompletionSerial++;

    const size_t reconciliationBudget =
	BObolLodQualityPolicy::staticPresentationRenderCostLimit(
	    provenRenderCost, elapsedNanoseconds,
	    this->d->stablePresentationFrameDeadlineNanoseconds);
    const bool deadlineHandoffReady =
	this->d->lodPresentationPolicy.noteFramePresented(
	    provenRenderCost, reconciliationBudget);
    if (deadlineHandoffReady) {
	/* The constrained retry frame is the proof which a deadline-created
	 * handoff was waiting for.  The admission pass normally ran before that
	 * frame and correctly left the ceiling armed; schedule the successor pass
	 * now so it can retire the handoff and restore the quiet presentation. */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPending = TRUE;
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodRetainedRefinementPending = FALSE;
	this->d->lodRetainedResidencyPending = FALSE;
	this->d->lodRetainedRefinementCutAdvanced = FALSE;
	this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodBudgetPolicy.resetPass();
	this->markProgressiveWorkPending();
    }

    /* Rendering time is not user-idle time.  Start the quiet clock only
     * after a frame for the newest camera has actually completed. */
    if (this->d->lodInteractive &&
	this->d->lodSettleAfterRenderSerial != 0 &&
	this->d->renderCompletionSerial >=
	    this->d->lodSettleAfterRenderSerial) {
	this->d->lodLastViewChangeMicroseconds = bu_gettime();
	this->d->lodSettleAfterRenderSerial = 0;
	if (this->d->lodViewDemandPolicy.noteMotionFrameSettled())
	    this->markProgressiveWorkPending();
    }

    /* A partial PoP result must be presented before the next, potentially
     * much larger prefix is requested.  This is presentation sequencing, not
     * capacity evidence: an exact selection/HUD frame which traversed the
     * current population satisfies the barrier just as an LoD-named frame
     * does. */
    const BObolLodFrameObligation::Completion frameCompletion =
	this->d->lodFrameObligation.complete(
	    this->d->renderCompletionSerial,
	    this->d->lodViewRevision.value(),
	    this->d->lodPolicyRevision.value());
    if (frameCompletion.retired) {
	const SbBool resumeResidentRefinement =
	    ((frameCompletion.reasons &
	      BObolLodFrameObligation::REASON_RESIDENT_REFINEMENT) != 0 ||
	     this->d->lodRetainedResidencyPending) ? TRUE : FALSE;
	this->d->lodBudgetPolicy.resetPass();
	const uint64_t responsiveFrame = 33333334ULL;
	int64_t cooldownMicroseconds = 0;
	if (elapsedNanoseconds > responsiveFrame) {
	    const uint64_t elapsedMicroseconds =
		elapsedNanoseconds / 1000ULL;
	    cooldownMicroseconds = static_cast<int64_t>(
		std::max<uint64_t>(50000ULL,
		    std::min<uint64_t>(2000000ULL, elapsedMicroseconds)));
	}
	const int64_t nowMicroseconds = bu_gettime();
	this->d->lodRefinementNotBeforeMicroseconds =
	    cooldownMicroseconds > 0 &&
	    nowMicroseconds <=
		std::numeric_limits<int64_t>::max() - cooldownMicroseconds ?
	    nowMicroseconds + cooldownMicroseconds : nowMicroseconds;
	/* Releasing a presentation barrier is another pass in the same camera,
	 * policy, and source epoch.  Its explicit pending cursor is the wakeup
	 * edge.  Invalidating the epoch witness while a calibration path had
	 * already armed that cursor made submitLodRequestsIfNeeded() classify the
	 * continuation as a view change during submission and append an O(N)
	 * source rescan.  Preserve an already pending pass selected by a stronger
	 * frame-completion policy. */
	if (resumeResidentRefinement)
	    this->d->lodRetainedResidencyPending = FALSE;
	if (!this->d->lodSubmissionPending) {
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionRescanPending = FALSE;
	    this->d->lodSubmissionPending = TRUE;
	}
	this->markProgressiveWorkPending();
    }

    this->d->lodPublicationPolicy.frameCompleted();
    this->scheduleResidentGrowthReallocationIfReady();

    /* A retained handoff pass can complete while its last result-publication
     * frame is still outstanding.  completePass() correctly refuses to use
     * that non-quiescent population as its final allocation proof, but after
     * this frame retires the publication barrier there may be no worker,
     * result, or submission cursor left to revisit the handoff.  The generic
     * refinement-frame liveness guard cannot discharge a planning
     * obligation: it merely redraws the same constrained population.
     *
     * Turn the publication-completion edge into the one missing retained
     * allocation transaction.  This changes only occurrence cuts within the
     * already resident population; it does not reload meshes or normalize
     * them to their minimum prefix.  The next complete pass either reconciles
     * the renderer ceiling or preserves it from measured deadline evidence. */
    const bool handoffServiceQuiescent = !this->d->lodService ||
	(this->d->lodService->activeTaskCountForGeneration(
	     this->d->lodActiveGeneration) == 0 &&
	 this->d->lodService->queuedResultCountForGeneration(
	     this->d->lodActiveGeneration) == 0 &&
	 this->d->lodResultsPending.load() == 0);
    const bool handoffPlanningReady =
	this->d->lodPresentationPolicy.handoffPending() &&
	!this->d->lodPresentationPolicy.handoffPresentationPending() &&
	!this->d->lodSubmissionPending &&
	!this->d->lodSubmissionRescanPending &&
	!this->d->lodFrameObligation.pending() &&
	!this->d->lodPublicationPolicy.pending() &&
	!this->d->lodResidentGrowthPolicy.pending() &&
	!this->d->lodBudgetPolicy.rescanAfterFrame() &&
	handoffServiceQuiescent;
    if (handoffPlanningReady) {
	const std::vector<SoBRLDatabaseSource *> sources =
	    controller_render_database_source_roots(this);
	if (!sources.empty()) {
	    const size_t delayedReconciliationBudget =
		this->d->lodPresentationPolicy.handoffReconciliationBudget();
	    if (delayedReconciliationBudget > 0)
		this->d->lodBudgetPolicy.requestPresentationReconciliation(
		    delayedReconciliationBudget);
	    else
		this->d->lodBudgetPolicy.requestRetainedReallocation();
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionRescanPending = FALSE;
	    this->d->lodSubmissionPending = TRUE;
	    this->d->lodPassAdmittedWork = FALSE;
	    this->d->lodRetainedRefinementPending = FALSE;
	    this->d->lodRetainedResidencyPending = FALSE;
	    this->d->lodRetainedRefinementCutAdvanced = FALSE;
	    this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	    this->d->lodBudgetPolicy.resetPass();
	    this->markProgressiveWorkPending();
	}
    }
    if (this->d->lodViewDemandPolicy.qualityProbePending()) {
	this->markProgressiveWorkPending();
    }
    /* synchronizePresentation() ran before this completed frame, so it
     * contains every immutable result accumulated by the owner thread. */
    this->d->reconcilePhase(BObolLodStateMachine::Event::FRAME_COMPLETED);
}

void
BObolViewController::scheduleLodRefinementFrame(const char *reason)
{
    /* This gate is useful only if a host presentation is guaranteed to
     * follow it.  Merely latching requestRender() is insufficient: a result
     * can be drained while progressiveWorkPending is already set, after
     * which the same advance clears that edge-triggered latch.  With no
     * unrelated Qt event, the controller would then wait forever for a frame
     * it never asked the host to schedule.
     *
     * Latch the render before invoking the callback so synchronous hosts also
     * observe a renderable request.  Do not move an existing gate forward;
     * every cut selected before the pending presentation belongs to the same
     * next frame. */
    const BObolLodFrameObligation::Reason obligationReason =
	this->d->lodRetainedResidencyPending ?
	    BObolLodFrameObligation::REASON_RESIDENT_REFINEMENT :
	    BObolLodFrameObligation::REASON_CUT_PRESENTATION;
    (void)this->d->lodFrameObligation.arm(obligationReason,
	this->d->renderCompletionSerial,
	this->d->lodViewRevision.value(),
	this->d->lodPolicyRevision.value());
    const char *frameReason = reason ? reason : "lod-refinement-frame";
    /* Preserve a more specific render reason already installed by result
     * application side effects such as residency eviction. */
    if (!this->isRenderRequested())
	this->requestRender(frameReason);
    this->d->reconcilePhase(BObolLodStateMachine::Event::WORK_SCHEDULED);
}

void
BObolViewController::scheduleResidentGrowthReallocationIfReady(void)
{
    if (!this->d->lodResidentGrowthPolicy.pending() ||
	!this->d->lodService)
	return;

    const bool streamIdle =
	this->d->lodService->activeTaskCountForGeneration(
	    this->d->lodActiveGeneration) == 0 &&
	this->d->lodService->queuedResultCountForGeneration(
	    this->d->lodActiveGeneration) == 0;
    /* Do not interrupt a structural/occurrence census or consume this edge
     * behind a presentation barrier.  The next progressive pump retries it
     * after the complete result wave and its coherent frame are visible. */
    const bool presentationReady =
	!this->d->lodSubmissionPending &&
	!this->d->lodFrameObligation.pending() &&
	!this->d->lodPublicationPolicy.pending() &&
	!this->d->lodBudgetPolicy.rescanAfterFrame();
    /* Pose-only motion owns a reversible renderer ceiling and deliberately
     * preserves retained occurrence cuts.  A zoom may spend newly available
     * residency during its bounded quality probes; otherwise retain this
     * obligation until the quiet transition. */
    const bool allocationAllowed =
	!this->d->lodInteractive ||
	this->d->lodViewDemandPolicy.interactionScaleChanged();
    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
	    this->d->lodViewRevision.value())) {
	static std::atomic<unsigned int> residentGrowthTraceCount(0);
	const unsigned int traceIndex = residentGrowthTraceCount.fetch_add(1);
	if (traceIndex < 256)
	    bu_log("BObol LoD resident growth pending stream_idle=%d "
		   "scene_ready=%d automatic=%d allocation_allowed=%d "
		   "submission=%d coverage=%d refinement_frame=%d "
		   "publication=%d budget_frame=%d\n",
		   streamIdle ? 1 : 0,
		   (presentationReady &&
		    this->d->lodCoveragePolicy.effectiveComplete()) ? 1 : 0,
		   this->d->lodAutoSubmit ? 1 : 0,
		   allocationAllowed ? 1 : 0,
		   this->d->lodSubmissionPending ? 1 : 0,
		   this->d->lodCoveragePolicy.effectiveComplete() ? 1 : 0,
		   this->d->lodFrameObligation.pending() ? 1 : 0,
		   this->d->lodPublicationPolicy.pending() ? 1 : 0,
		   this->d->lodBudgetPolicy.rescanAfterFrame() ? 1 : 0);
    }

    /* A newly resident suffix can arrive just after the pass which proved
     * minimum-mesh coverage.  If that pass ended with one or more missing
     * prefixes, coverage is deliberately false—but resident growth is now
     * the event which can satisfy it.  Waiting for coverage before scheduling
     * another census leaves a closed cycle: no task, result, frame, or cursor
     * remains to change either latch.  Resume (or start) the bounded coverage
     * pass first; only its completed proof may consume the scene-wide
     * reallocation below. */
    if (this->d->lodAutoSubmit && streamIdle && presentationReady &&
	allocationAllowed &&
	!this->d->lodCoveragePolicy.effectiveComplete()) {
	if (!this->d->lodCoveragePolicy.active()) {
	    this->d->lodCoveragePolicy.activate(false);
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	}
	this->d->lodSubmissionPending = TRUE;
	this->d->lodSubmissionRescanPending = FALSE;
	this->markProgressiveWorkPending();
	this->d->reconcilePhase(BObolLodStateMachine::Event::WORK_SCHEDULED);
	return;
    }
    const bool residentWorkReady = streamIdle && presentationReady &&
	this->d->lodCoveragePolicy.effectiveComplete();
    if (this->d->lodResidentGrowthPolicy.beginResidencyDrainIfReady(
	    this->d->lodAutoSubmit != FALSE, residentWorkReady,
	    allocationAllowed)) {
	/* Request every still-needed immutable suffix before recomputing the
	 * scene-wide cut distribution.  This pass uses the current unsatisfied
	 * frontier and preserves every visible cut; a bounded worker queue may
	 * require several such waves, but none performs an O(scene) allocation or
	 * invalidates the last complete framebuffer. */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	/* A residency drain is an ordinary suffix-request pass.  Carrying the
	 * preceding retained-admission mode through resetPass() immediately
	 * re-armed the O(scene) minimax allocator and turned every result batch
	 * into another BALANCING/REFINING cycle. */
	this->d->lodSubmissionRetainedAdmissionMode = FALSE;
	this->d->lodSubmissionRescanPending = FALSE;
	this->d->lodRetainedRefinementPending = FALSE;
	this->d->lodRetainedResidencyPending = FALSE;
	this->d->lodRetainedRefinementCutAdvanced = FALSE;
	this->d->lodRetainedRefinementBudgetBlocked = FALSE;
	this->d->lodPassAdmittedWork = FALSE;
	this->d->lodBudgetPolicy.resetPass();
	this->d->lodSubmissionPending = TRUE;
	this->d->lodResidentGrowthResidencyDrainActive = TRUE;
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
		this->d->lodViewRevision.value()))
	    bu_log("BObol LoD resident growth scheduled residency drain "
		   "view=%llu policy=%llu\n",
		   (unsigned long long)this->d->lodViewRevision.value(),
		   (unsigned long long)this->d->lodPolicyRevision.value());
	this->markProgressiveWorkPending();
	this->d->reconcilePhase(BObolLodStateMachine::Event::WORK_SCHEDULED);
	return;
    }
    if (!this->d->lodResidentGrowthPolicy.consumeIfReady(
	    this->d->lodAutoSubmit != FALSE,
	    residentWorkReady,
	    allocationAllowed))
	return;

    this->d->lodResidentGrowthResidencyDrainActive = FALSE;
    this->d->lodSubmissionSourceIndex = 0;
    this->d->lodSubmissionEntryOffset = 0;
    this->d->clearLodSubmissionPlan();
    this->d->lodSubmissionRescanPending = FALSE;
    this->d->lodRetainedRefinementPending = FALSE;
    this->d->lodRetainedResidencyPending = FALSE;
    this->d->lodRetainedRefinementCutAdvanced = FALSE;
    this->d->lodRetainedRefinementBudgetBlocked = FALSE;
    this->d->lodPassAdmittedWork = FALSE;
    this->d->lodBudgetPolicy.resetPass();
    this->d->lodBudgetPolicy.requestRetainedReallocation();
    /* This is an explicit pass in the current epoch.  Preserve the submitted
     * epoch witness so the wrapper cannot turn it into a view-change rescan. */
    this->d->lodSubmissionPending = TRUE;
    this->d->resetLodConvergenceFraction();
    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
	    this->d->lodViewRevision.value()))
	bu_log("BObol LoD resident growth scheduled scene-wide reallocation "
	       "view=%llu policy=%llu budget=%zu\n",
	       (unsigned long long)this->d->lodViewRevision.value(),
	       (unsigned long long)this->d->lodPolicyRevision.value(),
	       this->d->lodBudgetPolicy.currentBudget());
    this->markProgressiveWorkPending();
    this->d->reconcilePhase(BObolLodStateMachine::Event::WORK_SCHEDULED);
}

uint64_t
BObolViewController::getLastRenderTimeNanoseconds(void) const
{
    return this->d->lastRenderTimeNanoseconds;
}

uint64_t
BObolViewController::getSmoothedRenderTimeNanoseconds(void) const
{
    return this->d->smoothedRenderTimeNanoseconds;
}

uint64_t
BObolViewController::getLastBackgroundRenderTimeNanoseconds(void) const
{
    return this->d->lastBackgroundRenderTimeNanoseconds;
}

uint64_t
BObolViewController::getLastSceneRenderTimeNanoseconds(void) const
{
    return this->d->lastSceneRenderTimeNanoseconds;
}

uint64_t
BObolViewController::getLastProgressiveAdvanceTimeNanoseconds(void) const
{
    return this->d->lastProgressiveAdvanceTimeNanoseconds;
}

uint64_t
BObolViewController::getLastLodResultProcessingTimeNanoseconds(void) const
{
    return this->d->lastLodResultProcessingTimeNanoseconds;
}

uint64_t
BObolViewController::getLastProgressiveProviderTimeNanoseconds(void) const
{
    return this->d->lastProgressiveProviderTimeNanoseconds;
}

uint64_t
BObolViewController::getLastLodSubmissionTimeNanoseconds(void) const
{
    return this->d->lastLodSubmissionTimeNanoseconds;
}

uint64_t
BObolViewController::getLastPresentationSyncTimeNanoseconds(void) const
{
    return this->d->lastPresentationSyncTimeNanoseconds;
}

void
BObolViewController::noteFramePresented(void)
{
    const uint64_t now = this->beginRenderTiming();
    std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
    if (this->d->presentedFrameSerial != UINT64_MAX)
	this->d->presentedFrameSerial++;
    uint64_t generalDisplaySample = 0;
    const auto updateDisplayedCadence = [this](uint64_t interval) {
	uint64_t &displayed =
	    this->d->displayedPresentationIntervalNanoseconds;
	if (!displayed) {
	    displayed = interval;
	    return;
	}

	/* The control EMAs are deliberately quick: they let the LoD policy
	 * observe a change in renderer capacity within a few frames.  That same
	 * response makes a faceplate sampled four times per second look
	 * needlessly nervous.  Give displayed cadence a 750 ms elapsed-time
	 * constant.  Computing alpha from elapsed time, rather than using a
	 * fixed sample weight, gives the same visual response at 10, 60, or
	 * 240 FPS. */
	constexpr double timeConstantNanoseconds = 750000000.0;
	const double alpha = -std::expm1(
	    -static_cast<double>(interval) / timeConstantNanoseconds);
	const double next = static_cast<double>(displayed) +
	    alpha * (static_cast<double>(interval) -
		static_cast<double>(displayed));
	displayed = static_cast<uint64_t>(
	    std::max(1.0, std::floor(next + 0.5)));
    };

    if (this->d->lastPresentationTimestampNanoseconds &&
	now > this->d->lastPresentationTimestampNanoseconds) {
	const uint64_t interval =
	    now - this->d->lastPresentationTimestampNanoseconds;
	/* A retained view is event driven.  The first frame after an idle
	 * second is not a one-FPS rendering sample; it is a discontinuity
	 * between two bursts.  Excluding that gap keeps the faceplate useful
	 * while zooming/rotating without reintroducing an idle repaint timer. */
	if (interval > 1000 && interval <= 250000000ULL) {
	    this->d->smoothedPresentationIntervalNanoseconds =
		this->d->smoothedPresentationIntervalNanoseconds ?
		(this->d->smoothedPresentationIntervalNanoseconds * 4u +
		 interval) / 5u : interval;
	    if (!this->d->lodGestureActive)
		generalDisplaySample = interval;
	} else if (interval > 250000000ULL &&
		   !this->d->lodGestureActive) {
	    /* The next burst must establish its own display baseline.  Retaining
	     * the prior value here would make a newly exposed window or a gesture
	     * after a long idle period report stale historical cadence. */
	    this->d->displayedPresentationIntervalNanoseconds = 0;
	}
    }
    this->d->lastPresentationTimestampNanoseconds = now;

    /* Only an active input gesture supplies a continuous-demand FPS sample.
     * Event-driven stable views can legitimately go 100 ms between a timer,
     * screenshot, or overlay repaint; treating those idle gaps as renderer
     * capacity made the scene budget depend on test/poll cadence. */
    if (this->d->lodGestureActive) {
	if (this->d->lastInteractivePresentationTimestampNanoseconds &&
	    now > this->d->lastInteractivePresentationTimestampNanoseconds) {
	    const uint64_t interval =
		now - this->d->
		    lastInteractivePresentationTimestampNanoseconds;
	    if (interval > 1000 && interval <= 250000000ULL) {
		this->d->smoothedInteractivePresentationIntervalNanoseconds =
		    this->d->
			smoothedInteractivePresentationIntervalNanoseconds ?
		    (this->d->
			smoothedInteractivePresentationIntervalNanoseconds * 4u +
		     interval) / 5u : interval;
		/* Navigation brackets continuous demand explicitly.  Prefer that
		 * clean burst sample over the interval from an unrelated
		 * event-driven frame immediately before the gesture. */
		updateDisplayedCadence(interval);
	    }
	}
	this->d->lastInteractivePresentationTimestampNanoseconds = now;
    } else {
	this->d->lastInteractivePresentationTimestampNanoseconds = 0;
	if (generalDisplaySample)
	    updateDisplayedCadence(generalDisplaySample);
    }
}

uint64_t
BObolViewController::getPresentedFrameSerial(void) const
{
    std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
    return this->d->presentedFrameSerial;
}

uint64_t
BObolViewController::getRenderRequestSerial(void) const
{
    return this->renderRequestSerialGet();
}

uint64_t
BObolViewController::getRenderCompletionSerial(void) const
{
    return this->d->renderCompletionSerial;
}

uint64_t
BObolViewController::getLodSettleAfterRenderSerial(void) const
{
    return this->d->lodSettleAfterRenderSerial;
}

uint64_t
BObolViewController::getLodRefinementResumeAfterRenderSerial(void) const
{
    return this->d->lodFrameObligation.requiredRenderSerial();
}

uint64_t
BObolViewController::getSmoothedPresentationIntervalNanoseconds(void) const
{
    std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
    return this->d->smoothedPresentationIntervalNanoseconds;
}

uint64_t
BObolViewController::
getSmoothedInteractivePresentationIntervalNanoseconds(void) const
{
    std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
    return this->d->smoothedInteractivePresentationIntervalNanoseconds;
}

uint64_t
BObolViewController::getDisplayedPresentationIntervalNanoseconds(void) const
{
    std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
    return this->d->displayedPresentationIntervalNanoseconds;
}

int
BObolViewController::renderToImage(unsigned char **image,
				     int flip,
				     int alpha,
				     const SbColor *background,
				     SoDB::ContextManager *contextManager,
				     BObolProgressiveStatus *progressiveStatus)
{
    if (!image) {
	return BRLCAD_ERROR;
    }
    *image = NULL;

    if (!this->d->activeCamera || !this->getViewport() ||
	    !this->getRenderRoot()) {
	return BRLCAD_ERROR;
    }

    this->synchronizePresentation();
    /* LoD off -> force-realize the whole scene; LoD on -> coarse-first, letting
     * advanceProgressiveWork stream geometry in (matches the hardware and
     * headless render paths).  Default (no attachment / LoD off) is unchanged. */
    if (this->isForceRealizeDisplay())
	(void)this->realizePending();
    BObolProgressiveStatus localProgressiveStatus;
    /* Match renderPending() and the Qt endpoints: a presentation-only image
     * request must not opportunistically reopen an idle LoD coordinator.
     * In particular, selection and faceplate captures are observational once
     * their current view has settled.  Advancing unconditionally here could
     * start a fresh capacity/headroom probe after the caller had removed the
     * convergence HUD, causing the captured frame to contain a transient
     * responsiveness-limited indicator even though the geometry was already
     * terminal.  Real background work and capacity-relevant render requests
     * remain eligible for the normal bounded pre-render pump. */
    const BObolHostWorkSnapshot preRenderWork =
	this->getHostWorkSnapshot();
    if (preRenderWork.pumpPending() ||
	preRenderWork.capacitySampleRequested())
	(void)this->advanceProgressiveWork(NULL, &localProgressiveStatus);
    this->synchronizePresentation();
    if (progressiveStatus) {
	*progressiveStatus = localProgressiveStatus;
    }

    const SbViewportRegion &region = this->getViewportRegion();
    SbVec2s size = region.getViewportSizePixels();
    if (size[0] <= 0 || size[1] <= 0) {
	return BRLCAD_ERROR;
    }

    SoDB::ContextManager *resolvedManager = contextManager ?
	contextManager : this->getRenderContextManager();
    if (!resolvedManager)
	return BRLCAD_ERROR;

    const bool cacheRenderer =
	resolvedManager == this->getRenderContextManager();
    std::unique_ptr<SoOffscreenRenderer> overrideRenderer;
    SoOffscreenRenderer *renderer = NULL;
    bool configureRenderer = false;
    if (!cacheRenderer) {
	overrideRenderer.reset(new SoOffscreenRenderer(resolvedManager, region));
	renderer = overrideRenderer.get();
	configureRenderer = true;
    } else if (!this->d->imageRenderer ||
	this->d->imageRendererManager != resolvedManager) {
	delete this->d->imageRenderer;
	this->d->imageRenderer = new SoOffscreenRenderer(resolvedManager, region);
	this->d->imageRendererManager = resolvedManager;
	renderer = this->d->imageRenderer;
	configureRenderer = true;
    } else {
	this->d->imageRenderer->setViewportRegion(region);
	renderer = this->d->imageRenderer;
    }
    if (configureRenderer) {
	renderer->getGLRenderAction()->setTransparencyType(
	    this->d->transparencyEnabled ? SoGLRenderAction::BLEND :
	    SoGLRenderAction::NONE);
	renderer->getGLRenderAction()->setSmoothing(
	    this->d->antialiasingEnabled);
	renderer->getGLRenderAction()->setNumPasses(1);
    }
    const SbColor imageBottom = background ? *background :
	this->d->backgroundBottom;
    const SbColor imageTop =
	(background && *background != this->d->backgroundBottom) ?
	*background : this->d->backgroundTop;
    const SbBool gradient = imageBottom != imageTop;
    renderer->setComponents(SoOffscreenRenderer::RGB);
    renderer->setBackgroundColor(imageBottom);
    if (gradient)
	renderer->setBackgroundGradient(imageBottom, imageTop);
    else
	renderer->clearBackgroundGradient();

    const uint64_t requestSerial = this->renderRequestSerialGet();
    const uint64_t started = this->beginRenderTiming();
    const BObolViewLodState *presentationState = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    if (presentationState)
	presentationState->beginCadPresentationFrame();
    const SbBool rendered = renderer->render(this->getRenderRoot());
    const uint64_t completed = this->beginRenderTiming();
    if (presentationState)
	presentationState->refreshCadPresentationFrameStatus();
    const SbBool cadFrameIncomplete = presentationState &&
	!presentationState->lastCadPresentationFrameExact();
    if (!rendered) {
	return BRLCAD_ERROR;
    }
    if (cadFrameIncomplete) {
	/* Offscreen capture has no controller deadline and the renderer returned
	 * a complete pixel buffer.  A CAD assembly may nevertheless lack an exact
	 * submitted-work record (for example, a shared multi-view assembly which
	 * this view culled).  The image is presentable, but it cannot calibrate
	 * scene capacity.  Complete only the presentation barrier; classifying it
	 * as an interrupted frame would manufacture recovery work forever. */
	this->completeRenderTiming(started, FALSE);
    } else {
	this->d->lastBackgroundRenderTimeNanoseconds = 0;
	this->d->lastSceneRenderTimeNanoseconds =
	    completed >= started ? completed - started : 0;
	this->completeRenderTiming(started,
	    this->isForceRealizeDisplay() ? FALSE : TRUE);
    }

    const unsigned char *buffer = renderer->getBuffer();
    if (!buffer) {
	return BRLCAD_ERROR;
    }

    const size_t width = static_cast<size_t>(size[0]);
    const size_t height = static_cast<size_t>(size[1]);
    const size_t srcStride = width * 3;
    const size_t dstBpp = alpha ? 4 : 3;
    const size_t dstStride = width * dstBpp;
    unsigned char *out = static_cast<unsigned char *>(bu_calloc(
	    height * dstStride, sizeof(unsigned char),
	    "bobol viewport image"));

    for (size_t y = 0; y < height; y++) {
	const size_t srcY = flip ? (height - y - 1) : y;
	const unsigned char *src = buffer + srcY * srcStride;
	unsigned char *dst = out + y * dstStride;
	if (alpha) {
	    for (size_t x = 0; x < width; x++) {
		dst[x * 4 + 0] = src[x * 3 + 0];
		dst[x * 4 + 1] = src[x * 3 + 1];
		dst[x * 4 + 2] = src[x * 3 + 2];
		dst[x * 4 + 3] = 255;
	    }
	} else {
	    std::memcpy(dst, src, srcStride);
	}
    }

    if (!localProgressiveStatus.hasMore)
	this->clearRenderRequestIfUnchanged(requestSerial);
    *image = out;
    return BRLCAD_OK;
}

SbBool
BObolViewController::isRenderRequested(void) const
{
    if (!this->d->endpointGraphicalRenderingEnabled.load(
	    std::memory_order_acquire))
	return FALSE;

    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    return this->d->renderRequested;
}

SbString
BObolViewController::getRenderReason(void) const
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    return this->d->renderReason;
}
