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
#include "presentation_preparation_private.h"
#include "retained_allocation_private.h"
#include "view_controller_private.h"
#include "raytrace.h"
#include "rt/view.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <climits>
#include <cmath>
#include <condition_variable>
#include <cstdlib>
#include <cstdint>
#include <deque>
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
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoPointLight.h>
#include <Inventor/nodes/SoSpotLight.h>

static const uint64_t controller_structural_repair_deadline_multiplier = 2;
/* An exact point census must visit the complete occurrence population before
 * it can select a safe aggregation threshold.  On a software renderer that
 * finite, quiet-only traversal can exceed the ordinary presentation deadline
 * even at the minimum mesh cut.  Bound the one-shot census to 2.5 FPS rather
 * than restarting it forever at the normal 10 FPS hard limit. */
static const uint64_t controller_point_census_deadline_multiplier = 4;

struct ControllerStructuralFrontierIdentity {
    struct SourcePlan {
	uint64_t routingId = 0;
	uint64_t populationEpoch = 0;
	std::vector<size_t> entries;

	bool operator==(const SourcePlan &other) const
	{
	    return this->routingId == other.routingId &&
		this->populationEpoch == other.populationEpoch &&
		this->entries == other.entries;
	}
    };

    size_t visibleCount = 0;
    std::array<size_t,
	BObolViewLodState::CadStructuralProjectionHistogram::BucketCount>
	cumulativeCount = {};
    std::vector<SourcePlan> plans;

    bool operator==(const ControllerStructuralFrontierIdentity &other) const
    {
	return this->visibleCount == other.visibleCount &&
	    this->cumulativeCount == other.cumulativeCount &&
	    this->plans == other.plans;
    }
};

static uint64_t
controller_structural_frontier_identity(
    const BObolViewLodState::CadStructuralProjectionHistogram &projection,
    const std::vector<std::pair<SoBRLDatabaseSource *, std::vector<size_t>>>
	&plans)
{
    ControllerStructuralFrontierIdentity identity;
    identity.visibleCount = projection.visibleCount;
    identity.cumulativeCount = projection.cumulativeCount;
    identity.plans.reserve(plans.size());
    for (const auto &plan : plans) {
	ControllerStructuralFrontierIdentity::SourcePlan source;
	source.routingId = plan.first ?
	    plan.first->getCompactSourceRoutingId() : 0;
	source.populationEpoch = plan.first ?
	    plan.first->getCompactPopulationEpoch() : 0;
	source.entries = plan.second;
	identity.plans.push_back(std::move(source));
    }

    struct Record {
	uint64_t token = 0;
	ControllerStructuralFrontierIdentity identity;
    };
    static constexpr size_t retainedIdentityLimit = 64;
    static std::mutex registryMutex;
    static std::deque<Record> registry;
    static uint64_t nextToken = 1;
    std::lock_guard<std::mutex> lock(registryMutex);
    for (const Record &record : registry) {
	if (record.identity == identity)
	    return record.token;
    }
    Record record;
    record.identity = std::move(identity);
    record.token = bobol_nonzero_identity_take(nextToken);
    registry.push_back(std::move(record));
    if (registry.size() > retainedIdentityLimit)
	registry.pop_front();
    return registry.back().token;
}

enum class ControllerStructuralSelectionMode {
    UNCOLLAPSED,
    ABOVE_PIXEL_THRESHOLD
};

struct ControllerStructuralSelection {
    std::vector<std::pair<SoBRLDatabaseSource *, std::vector<size_t>>> plans;
    size_t totalOccurrenceCount = 0;
    size_t plannedOccurrenceCount = 0;
    bool exact = true;
};

static ControllerStructuralSelection
controller_select_structural_occurrences(
    BObolViewController *controller, const BObolViewLodState *viewLodState,
    ControllerStructuralSelectionMode mode, float pixelThreshold = 1.0f,
    size_t maximumPlannedOccurrences = SIZE_MAX)
{
    ControllerStructuralSelection selection;
    if (!controller || !viewLodState) {
	selection.exact = false;
	return selection;
    }

    const std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_source_roots(controller);
    for (SoBRLDatabaseSource *source : sources) {
	if (!source || !source->hasCompactInstanceIndex())
	    continue;
	SoBRLCadAssembly *assembly = dynamic_cast<SoBRLCadAssembly *>(
	    viewLodState->findCadPresentation(source));
	if (!assembly)
	    continue;

	std::vector<SbString> occurrenceKeys;
	if (mode == ControllerStructuralSelectionMode::UNCOLLAPSED) {
	    assembly->lastUncollapsedStructuralProxyOccurrenceKeys(
		occurrenceKeys);
	} else if (!assembly->lastStructuralProxyOccurrenceKeysAbovePixels(
		pixelThreshold, occurrenceKeys)) {
	    selection.exact = false;
	}

	std::vector<size_t> entryIndices;
	entryIndices.reserve(occurrenceKeys.size());
	for (const SbString &occurrenceKey : occurrenceKeys) {
	    size_t entryIndex = 0;
	    if (source->getCompactInstanceIndex(
		    occurrenceKey.getString(), entryIndex)) {
		entryIndices.push_back(entryIndex);
	    } else {
		selection.exact = false;
	    }
	}
	std::sort(entryIndices.begin(), entryIndices.end());
	entryIndices.erase(
	    std::unique(entryIndices.begin(), entryIndices.end()),
	    entryIndices.end());
	if (entryIndices.empty())
	    continue;
	selection.totalOccurrenceCount = entryIndices.size() >
	    SIZE_MAX - selection.totalOccurrenceCount ? SIZE_MAX :
	    selection.totalOccurrenceCount + entryIndices.size();
	const size_t remainingPlanCapacity =
	    selection.plannedOccurrenceCount >= maximumPlannedOccurrences ? 0 :
	    maximumPlannedOccurrences - selection.plannedOccurrenceCount;
	if (entryIndices.size() > remainingPlanCapacity)
	    entryIndices.resize(remainingPlanCapacity);
	if (entryIndices.empty())
	    continue;
	selection.plannedOccurrenceCount += entryIndices.size();
	selection.plans.emplace_back(source, std::move(entryIndices));
    }
    return selection;
}

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
    BObolLodControlTransitionScope controlTransition(this);
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
    SbBool lodPlanningRelevant = TRUE;
    if (!this->consumeRenderRequest(&renderReasonCopy,
	    &lodCapacityRelevant, &lodPlanningRelevant))
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
    if (presentationState)
	presentationState->beginCadPresentationFrame();
    this->d->renderManager->render(static_cast<SbBool>(FALSE), static_cast<SbBool>(FALSE));
    const uint64_t sceneCompleted = this->beginRenderTiming();
    if (presentationState)
	presentationState->refreshCadPresentationFrameStatus();
    const BObolCadPreparationProgress cadPreparation = presentationState ?
	presentationState->cadPresentationPreparationProgress() :
	BOBOL_CAD_PREPARATION_NONE;
    const SbBool cadFrameIncomplete = presentationState &&
	!presentationState->lastCadPresentationFrameExact();
    /* Coin checks the abort callback only between traversal steps.  Treat a
     * traversal which returns after its deadline as interrupted too: a single
     * long callback or native draw must not turn an over-budget presentation
     * into a successful calibration sample. */
    const SbBool deadlineExpired = renderAction && deadlineDuration &&
	sceneCompleted >= deadlineContext.deadlineNanoseconds;
    const SbBool interrupted = (renderAction && renderAction->hasTerminated()) ||
	deadlineExpired || cadFrameIncomplete;
    if (renderAction && deadlineDuration)
	renderAction->setAbortCallback(
	    deadlineContext.previous, deadlineContext.previousData);
    this->d->lastBackgroundRenderTimeNanoseconds =
	backgroundCompleted >= started ?
	    backgroundCompleted - started : 0;
    this->d->lastSceneRenderTimeNanoseconds =
	sceneCompleted >= backgroundCompleted ?
	    sceneCompleted - backgroundCompleted : 0;
    const BObolPresentationTimingContext timingContext(
	lodCapacityRelevant ?
	    BObolLodCapacityRelevance::RELEVANT :
	    BObolLodCapacityRelevance::EXCLUDED,
	lodPlanningRelevant ?
	    BObolLodPlanningRelevance::RELEVANT :
	    BObolLodPlanningRelevance::EXCLUDED,
	presentationState &&
	    presentationState->lastCadPresentationFrameExecuted() ?
	    BObolCadPresentationExecution::EXECUTED :
	    BObolCadPresentationExecution::NOT_EXECUTED,
	cadPreparation,
	cadFrameIncomplete ?
	    BObolCadPresentationCompleteness::INCOMPLETE :
	    BObolCadPresentationCompleteness::EXACT);
    if (interrupted) {
	this->notePresentationRenderInterrupted(
	    sceneCompleted > started ? sceneCompleted - started : 1,
	    timingContext);
	return FALSE;
    }
    this->completeRenderTiming(started, timingContext);
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
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    if (this->d->interactivePresentationFrameDeadlineNanoseconds ==
	    interactiveNanoseconds &&
	this->d->stablePresentationFrameDeadlineNanoseconds ==
	    stableNanoseconds &&
	this->d->prominentQualityFrameDeadlineNanoseconds ==
	    stableNanoseconds)
	return;
    this->d->resetLodViewQualityHistory();
    this->d->interactivePresentationFrameDeadlineNanoseconds =
	interactiveNanoseconds;
    this->d->stablePresentationFrameDeadlineNanoseconds =
	stableNanoseconds;
    /* This lower-level API is an embedding application's strict endpoint
     * contract.  Unlike the user-facing FPS policy, it also bounds the
     * prominent-quality exception. */
    this->d->prominentQualityFrameDeadlineNanoseconds = stableNanoseconds;
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

    /* A quiet view whose retained CAD population is already at the
     * irreducible presentation floor has no lower-cost successor for a
     * deadline abort to request.  Repeating a finite deadline in that state
     * only redraws and aborts the same image indefinitely (notably a
     * hidden-line software frame).  Let one static frame complete and retain
     * it; active input still uses the bounded interactive deadline below.
     *
     * Point aggregation is the one remaining population reduction after a
     * mesh cut reaches its floor, so preserve the deadline when that policy
     * is applicable. */
    const BObolViewLodState *presentationState = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    /* Mere eligibility is not a successor.  Large scenes are often eligible
     * for point aggregation long after the policy has settled; retaining the
     * frame deadline in that state aborts the same irreducible software frame
     * forever.  Only an armed classifier/recovery can still reduce this
     * population. */
    const bool pointAggregationCanReducePresentation =
	!BObolLodAdmissionPlanner::pointAtMaximumPixelThreshold(
	    this->d->lodPresentationPointProxyPixelThreshold);
    const bool pointAggregationTransitionPending =
	pointAggregationCanReducePresentation &&
	(this->d->lodPointAdmissionFrame.pending() ||
	 this->d->lodPointQualityPhase.pending() ||
	 this->d->lodDiscoveryPointProxyPixelThreshold > 1.0f);
    const bool quietIrreduciblePresentation =
	BObolLodQualityPolicy::quietPresentationIrreducible(
	    this->d->lodInteractionSession.active(),
	    presentationState != NULL,
	    this->d->lodInteractiveProgressiveCeiling,
	    presentationState ? presentationState->activeRenderCost() : 0,
	    presentationState ? presentationState->minimumActiveRenderCost() : 0,
	    presentationState ? presentationState->cadProgressivePayloadCount() : 0,
	    pointAggregationTransitionPending);
    if (quietIrreduciblePresentation)
	return 0;

    /* A handoff replaces the retained framebuffer only after a complete
     * successor frame.  It must obey the same static deadline as the
     * allocation which selected that successor; otherwise the HUD advertises
     * a bounded quality contract while the endpoint can still monopolize the
     * GUI thread indefinitely.  Deadline recovery preserves the preceding
     * complete framebuffer and reconciles the richest affordable cut. */

    /* Structural repair and point classification are bounded, one-shot
     * coverage transactions.  A large software frame can legitimately take
     * longer than the normal interactive-safe static deadline while visiting
     * its complete occurrence population.  Give those finite transactions
     * enough time to produce their exact certificate; ordinary static and
     * interactive frames retain their existing deadlines. */
    const bool completingStructuralCoverage =
	!this->d->lodInteractionSession.active() &&
	(this->d->lodStructuralRepair.active() ||
	 (this->d->lodPointQualityPhase.adaptiveCalibrationPending() &&
	  BObolLodAdmissionPlanner::pointAtMaximumPixelThreshold(
	      this->d->lodPresentationPointProxyPixelThreshold)));
    const bool completingPointCensus =
	!this->d->lodInteractionSession.active() &&
	(this->d->lodPointAdmissionFrame.pending() ||
	 this->d->lodPointQualityPhase.adaptiveCalibrationPending());
    if (!this->d->lodInteractionSession.active()) {
	const uint64_t staticQualityDeadline =
	    this->d->staticQualityPresentationDeadline();
	uint64_t finiteCensusDeadline = 0;
	if (completingPointCensus) {
	    finiteCensusDeadline = BObolLodAdmissionPlanner::
		structuralRepairPresentationDeadline(false,
		    this->d->stablePresentationFrameDeadlineNanoseconds, 0,
		    controller_point_census_deadline_multiplier);
	} else {
	    const bool directMeshScene = controller_lod_mesh_first_scene_safe(
		controller_render_database_source_roots(this));
	    finiteCensusDeadline = BObolLodAdmissionPlanner::
		structuralRepairPresentationDeadline(directMeshScene,
		    this->d->stablePresentationFrameDeadlineNanoseconds,
		    this->d->prominentQualityPresentationDeadline(),
		    controller_structural_repair_deadline_multiplier);
	}
	return BObolLodQualityPolicy::quietPresentationDeadline(
	    staticQualityDeadline, finiteCensusDeadline,
	    completingStructuralCoverage || completingPointCensus,
	    this->d->lodStaticQualityTrial.usesStaticDeadline() ||
		this->d->lodPointQualityPhase.staticCalibrationPending());
    }

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
    uint64_t elapsedNanoseconds,
    const BObolPresentationTimingContext &context)
{
    BObolLodControlTransitionScope controlTransition(this);
    if (context.cadPresentationCompleteness ==
	    BObolCadPresentationCompleteness::INCOMPLETE) {
	/* Population producers may supersede this target, but none may erase the
	 * requirement for a later exact current frame.  Keep that requirement
	 * separate from the exclusive interrupted-replay gate below: resumable
	 * command preparation freezes owner-thread mutation, whereas an assembly
	 * omitted from a partial traversal may safely be replaced. */
	this->d->requireExactPresentationFrame();
    }
    const SbBool cadDrawAttempted =
	context.cadPresentationExecution ==
	    BObolCadPresentationExecution::EXECUTED ? TRUE : FALSE;
    const BObolCadPreparationProgress cadPreparation =
	context.cadPreparation;
    const SbBool lodCapacityRelevant =
	context.lodCapacityRelevance == BObolLodCapacityRelevance::RELEVANT ?
	    TRUE : FALSE;
    const SbBool lodPlanningRelevant =
	context.lodPlanningRelevance == BObolLodPlanningRelevance::RELEVANT ?
	    TRUE : FALSE;

    if (!elapsedNanoseconds)
	return;
    BObolRenderClaimCompletionScope claimedRenderCompletion(this);
    bobol_saturating_counter_advance(
	this->d->interruptedPresentationFrameCount);
    this->d->lastInterruptedPresentationTimeNanoseconds =
	elapsedNanoseconds;
    if (!lodCapacityRelevant) {
	/* A presentation patch may still encounter SoCADAssembly's resumable
	 * command preparation.  Finish that patch, but do not let it cancel a
	 * geometry headroom witness, lower a PoP ceiling, or alter a scene budget. */
	if (lodPlanningRelevant)
	    this->requestLodPresentationRender("render-presentation-replay");
	else
	    this->requestPresentationRender("render-presentation-replay");
	return;
    }
    bobol_saturating_counter_advance(
	this->d->consecutiveInterruptedPresentationFrames);
    const BObolViewLodState *state = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    const SbBool interactive = this->d->lodInteractionSession.active();
    const SbBool cadPreparationAdvanced =
	BObolCadPreparationPolicy::permitsUnchangedRetry(cadPreparation) ?
	    TRUE : FALSE;
    const uint64_t staticQualityDeadline =
	this->d->staticQualityPresentationDeadline();
    const uint64_t recoveryDeadline = interactive ?
	this->d->interactivePresentationFrameDeadlineNanoseconds :
	staticQualityDeadline;
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
    const bool certifiedUnconstrainedCandidate = state &&
	this->d->lodInteractiveProgressiveCeiling < 0 &&
	this->d->retainedAllocationCertificateCurrent(state) &&
	this->d->retainedAllocationCutsApplied(state);
    const bool certifiedUnconstrainedPresentation =
	certifiedUnconstrainedCandidate &&
	this->d->retainedAllocationPresentationRealized(state);
    int interruptedCorrectedCeiling = interruptedPresentedMaximum > 0 ?
	interruptedPresentedMaximum - 1 : interruptedPresentedMaximum;
    size_t interruptedEstimatedCost = 0;
    size_t interruptedTargetCost = 0;
    /*
     * A PoP cut ordinal is a precision boundary, not a cost ratio.  Once the
     * cache supports dozens of cuts, subtracting one ordinal per missed frame
     * can spend seconds walking past cuts which remove almost no work.  The
     * view state maintains the exact retained cost at every renderer-wide
     * ceiling.  For a draw-capacity miss, scale the attempted cost by
     * the observed deadline ratio and select the richest prefix which meets
     * that corrected cost in one step.  This is equally important during
     * input: walking one ordinal at a time can consume several complete
     * deadlines without presenting a new camera pose.  Preparation-only
     * retries remain unchanged below, and the hard deadline still corrects
     * an optimistic prediction.
     *
     * This aggregate includes off-frustum retained occurrences and excludes
     * renderer point collapse, so it is deliberately conservative.  It does
     * not change occurrence cuts or resident data; later completed-frame
     * headroom policy may recover fidelity without cache I/O.
     */
    /* A retained allocation is an exact attempted-population coordinate even
     * when every occurrence is already at PoP cut zero.  The ordinal therefore
     * cannot gate deadline evidence: a 50k minimum-prefix wireframe may exceed
     * the endpoint deadline while maximumActiveProgressiveCut() is zero.  In
     * that case the allocation cost is precisely the evidence the bounded
     * search needs to reject the candidate.  Requiring a positive cut skipped
     * the search, requested the unchanged frame again, and starved the
     * occurrence-allocation producer indefinitely. */
    const bool deadlineCostAvailable = certifiedUnconstrainedCandidate ||
	interruptedPresentedMaximum > 0;
    if (state && !cadPreparationAdvanced && deadlineCostAvailable &&
	elapsedNanoseconds > 0 && recoveryDeadline > 0) {
	if (certifiedUnconstrainedCandidate) {
	    interruptedEstimatedCost =
		allocationCertificate.selectedPresentationCost;
	} else {
	    const size_t activeCadCost = state->activeCadRenderCost();
	    const size_t presentedCadEstimate =
		state->cadRenderCostAtProgressiveCutCeiling(
		    interruptedPresentedMaximum);
	    interruptedEstimatedCost = BObolLodAdmissionPlanner::
		canonicalSceneCostAtCadCeiling(interruptedActiveCost,
		    activeCadCost, presentedCadEstimate);
	}
	interruptedTargetCost =
	    BObolLodQualityPolicy::deadlineRecoveryCostLimit(
		interruptedEstimatedCost, elapsedNanoseconds,
		recoveryDeadline);
	const size_t sceneBudget = this->d->lodAdmissionEvidence.capacity().currentBudget();
	if (sceneBudget != SIZE_MAX)
	    interruptedTargetCost = std::min(
		interruptedTargetCost, sceneBudget);
	if (interruptedPresentedMaximum > 0) {
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
	    const int safeProgressiveCeiling =
		this->d->lodDeadlineSafePresentation.ceilingFor(
		    this->d->lodViewRevision.value(),
		    this->d->lodViewQualityDomainRevision);
	    if (!interactive && safeProgressiveCeiling >= 0)
		/* The prior completed frame proves this cut is affordable for the
		 * unchanged camera and source population.  Preserve that visual floor
		 * when a richer speculative cut misses; the enclosing clamp still
		 * forces a retry below the cut which actually missed. */
		interruptedCorrectedCeiling = std::max(
		    interruptedCorrectedCeiling,
		    safeProgressiveCeiling);
	}
    }
    /* Strict same-target retained progress is not draw-capacity evidence.
	 * Large command records are deliberately built in resumable deadline
	 * slices and may need more than one slice at 150k+ occurrences.  The typed
	 * certificate grants a retry only when finite remaining work decreases;
	 * activity counters and superseded targets have no policy authority. */
    const SbBool transientPresentationRetry =
	BObolLodQualityPolicy::retryTransientPresentation(
	    interactive != FALSE,
	    this->d->consecutiveInterruptedPresentationFrames,
	    cadPreparationAdvanced != FALSE,
	this->d->lodPresentationTransaction.publicationFramePending(),
	this->d->lodPresentationTransaction.barrierPending(),
	this->d->lodPointQualityPhase.presentationPending(),
	this->d->lodStaticQualityTrial.inProgress()) ?
	TRUE : FALSE;
    /* Keep a resumable static-quality traversal on the exact cut whose atlas
     * ranges it is preparing.  Resident workers may publish a richer suffix
     * meanwhile, but the append-only lineage preserves this cut and the
     * presentation ceiling prevents the next retry from chasing the newer
     * generation before the current one has produced a capacity sample.
     * A completed frame advances or retires this ceiling through the ordinary
     * static-overscan staircase below. */
    if (!interactive && cadDrawAttempted && cadPreparationAdvanced && state &&
	this->d->lodStaticQualityTrial.probing() &&
	interruptedPresentedMaximum >= 0) {
	const int currentCeiling =
	    this->d->lodInteractiveProgressiveCeiling;
	const int preparationCeiling = currentCeiling >= 0 ?
	    std::min(currentCeiling, interruptedPresentedMaximum) :
	    interruptedPresentedMaximum;
	this->d->publishCadProgressiveCeiling(preparationCeiling);
    }
    /* A preparation/publication replay has not tested renderer capacity and
     * must retain its pending witness.  A steady quiet draw which crosses the
     * hard endpoint deadline has tested it: record a negative proof instead
     * of merely cancelling the pending retry.  Cancelling forgot the failed
     * population, so the cheap constrained successor immediately re-armed the
     * same rich frame and made large scenes visibly alternate coarse/fine.
     * Static quality is likewise one-shot for this camera/policy epoch after
     * such a miss; occurrence-local recovery remains free to redistribute
     * the demonstrated safe allowance. */
    const bool deadlineServiceProducer =
	this->d->lodAvailabilityLedger.resultsPending() != 0 ||
	(this->d->lodService &&
	 (this->d->lodService->activeTaskCountForGeneration(
	      this->d->lodActiveGeneration) > 0 ||
	  this->d->lodService->queuedResultCountForGeneration(
	      this->d->lodActiveGeneration) > 0));
    const bool deadlineSubmissionPaused =
	BObolLodAdmissionPlanner::presentationPausesSubmission(
	    this->d->lodPointAdmissionFrame.pending(),
	    this->d->lodPointQualityPhase.presentationPending(),
	    this->d->lodAdmissionEvidence.capacity().capacityAllocationPending(),
	    this->d->lodAdmissionEvidence.capacity().presentationFramePending(),
	    state && state->hasCadPresentationAssemblies());
    /* Defer capacity recovery only to an owner which can publish a different
     * future population.  A requested publication frame, retained-growth
     * allocation, residency drain, or handoff is downstream of this exact
     * presentation barrier; counting one as a producer makes an over-deadline
     * minimum cut wait on itself forever. */
    const bool deadlinePopulationWorkPending =
	BObolLodProducerPolicy::ownsFutureFrame(
	    this->d->lodSubmissionPass.active(), deadlineSubmissionPaused,
	    this->d->lodAvailabilityLedger.providerPendingCount() > 0,
	    deadlineServiceProducer,
	    this->d->lodPresentationTransaction.
		publicationAwaitingFrameRequest());
    const BObolLodQualityPolicy::DeadlineSuccessor deadlineSuccessor =
	BObolLodQualityPolicy::deadlineSuccessor(
	    transientPresentationRetry != FALSE,
	    deadlinePopulationWorkPending);
    /* A replay is meaningful only when population producers are quiescent.
     * While source discovery, task submission, or result publication is
     * active, stopping the progressive owner thread to replay an obsolete
     * draw serializes all realization behind the renderer. */
    if (cadDrawAttempted && !deadlinePopulationWorkPending)
	this->d->lodInterruptedPresentationReplay.arm();
    /* An interrupted frame while discovery, source submission, or result
     * publication is still active does not measure a stable presentation.
     * In particular, a large hidden-line source can spend its deadline
     * compiling or publishing its next structural batch while the current
     * CAD population is nearly empty.  Feeding that interruption to the
     * capacity controller repeatedly shrinks its scene budget and eventually
     * makes the coverage pass unable to admit any useful mesh.  Capacity
     * evidence is valid only after all population producers are quiet. */
    const bool capacityRelevantDeadlineMiss =
	!interactive && cadDrawAttempted &&
	deadlineSuccessor == BObolLodQualityPolicy::DeadlineSuccessor::
	    RECOVER_PRESENTATION;
    BObolLodCapacityEvidence::CompletedFrameDecision
	deadlineCapacityDecision;
    if (capacityRelevantDeadlineMiss && interruptedEstimatedCost > 1) {
	BObolLodCapacityEvidence::DeadlineMissInputs miss;
	miss.attemptedBudget = interruptedEstimatedCost;
	const bool staticDeadlineMiss = staticQualityDeadline >
	    this->d->stablePresentationFrameDeadlineNanoseconds;
	miss.staticDeadline = staticDeadlineMiss;
		/* A current retained-allocation certificate is an exact coordinate in
		 * the bounded capacity search even though the interrupted traversal did
		 * not finish an exact framebuffer.  Enter the bounded certificate on the
		 * first such miss.  A scalar correction has no owner which can prove its
		 * endpoint, and near-deadline software frames can otherwise repeat that
		 * correction until they become ownerless at 99 percent.  The search still
		 * begins with the preferred target and may transition once to the longer
		 * static-image deadline. */
		if (certifiedUnconstrainedCandidate &&
		    allocationCertificate.certifiedPresentationBudget > 0 &&
	    allocationCertificate.pixelDemandPresentationCost > 0 &&
	    this->d->quietAllocationTargetFps() > 0.0f) {
	    const uint64_t preferredTargetNanoseconds =
		static_cast<uint64_t>(1000000000.0L /
		    static_cast<long double>(
			this->d->quietAllocationTargetFps()));
	    const size_t capacityMinimumBudget =
		this->d->lodAdmissionEvidence.capacity().
		    capacitySearchMinimumBudget(interruptedMinimumCost,
			allocationCertificate.protectedFloorBudget,
			allocationCertificate.protectedFloorIdentity);
	    miss.searchKey = BObolLodCapacitySearchCertificate::keyFor(
		this->d->admissionRevisionStamp(), preferredTargetNanoseconds,
		this->d->prominentQualityPresentationDeadline(),
		this->d->lodPresentationPointProxyPixelThreshold,
		allocationCertificate.maximumCapacitySearchBudget(),
		capacityMinimumBudget,
		this->d->lodRetainedViewContinuity.startCapacityAtStatic());
	    miss.searchKey.preferredBudgetCeiling =
		this->d->lodAdmissionEvidence.capacity().
		    deadlineCapacityCeiling();
	    miss.searchKey.maximumBudgetCeiling =
		this->d->lodAdmissionEvidence.capacity().
		    staticDeadlineCapacityCeiling();
	    miss.candidateBudget =
		allocationCertificate.certifiedPresentationBudget;
	}
	deadlineCapacityDecision = this->d->recordDeadlineCapacityMiss(miss);
	const BObolLodAdmissionPlan rejectionPlan =
	    BObolLodAdmissionPlanner::recordHeadroomRejection(
		this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor, this->d->lodViewRevision,
		this->d->lodPolicyRevision, interruptedActiveCost,
		this->d->lodAdmissionEvidence.capacity().currentBudget());
	this->d->commitAdmissionPlan(rejectionPlan);
    } else if (interactive) {
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::CANCEL_HEADROOM_RETRY);
    }
    if (getenv("BOBOL_LOD_TRACE_DEADLINE")) {
	const BObolLodCapacitySearchCertificate &deadlineSearch =
	    this->d->lodAdmissionEvidence.capacity().capacitySearch();
	bu_log("BObol LoD render deadline elapsed_ms=%.3f "
	       "deadline_ms=%.3f "
	       "draw=%d preparation=%d interactive=%d "
	       "active_max=%d presented_max=%d ceiling=%d "
	       "static_active=%d static_rejected=%d headroom=%d "
	       "floor_budget=%zu floor_misses=%u floor_rejected=%d "
	       "steady_capacity_ceiling=%zu static_capacity_ceiling=%zu "
	       "certified=%d "
	       "active_cost=%zu minimum_cost=%zu handoff=%d "
	       "estimated_cost=%zu target_cost=%zu corrected=%d "
	       "search_phase=%u goal=%u candidate=%zu safe=%zu unsafe=%zu "
	       "measured=%u samples_remaining=%u\n",
	       elapsedNanoseconds / 1000000.0,
	       staticQualityDeadline / 1000000.0,
	       cadDrawAttempted ? 1 : 0,
	       cadPreparationAdvanced ? 1 : 0,
	       interactive ? 1 : 0,
	       interruptedActiveMaximum,
	       interruptedPresentedMaximum,
	       this->d->lodInteractiveProgressiveCeiling,
	       this->d->lodStaticQualityTrial.inProgress() ? 1 : 0,
	       this->d->lodStaticQualityTrial.capacityRejected() ? 1 : 0,
	       this->d->lodAdmissionEvidence.headroom().retryPending() ? 1 : 0,
	       this->d->lodAdmissionEvidence.capacity().retainedQualityFloorBudget(),
	       this->d->lodAdmissionEvidence.capacity().retainedQualityFloorMissCount(),
	       this->d->lodAdmissionEvidence.capacity().retainedQualityFloorRejected() ? 1 : 0,
	       this->d->lodAdmissionEvidence.capacity().deadlineCapacityCeiling(),
	       this->d->lodAdmissionEvidence.capacity().staticDeadlineCapacityCeiling(),
	       certifiedUnconstrainedPresentation ? 1 : 0,
	       state ? state->activeRenderCost() : 0,
	       state ? state->minimumActiveRenderCost() : 0,
	       this->d->lodPresentationPolicy.handoffPending() ? 1 : 0,
	       interruptedEstimatedCost, interruptedTargetCost,
	       interruptedCorrectedCeiling,
	       static_cast<unsigned int>(deadlineSearch.phase()),
	       static_cast<unsigned int>(deadlineSearch.goal()),
	       deadlineSearch.candidateBudget(), deadlineSearch.safeBudget(),
	       deadlineSearch.unsafeBudget(),
	       deadlineSearch.measuredCandidateCount(),
	       deadlineSearch.samplesRemaining());
    }
    if (deadlineCapacityDecision.searchResult !=
	    BObolLodCapacitySearchCertificate::Result::NONE) {
	/* A bounded capacity search has consumed this hard miss and selected its
	 * sole successor.  Do not also install the generic renderer-ceiling
	 * recovery below: that second owner starts a handoff allocation, resets the
	 * search population, and can alternate the same two candidates forever.
	 * The last complete framebuffer remains visible while the selected
	 * occurrence allocation is planned. */
	this->d->lodInterruptedPresentationReplay.retire();
	this->d->lodPresentationTransaction.reset();
	this->d->lodPresentationPolicy.cancelHandoff();
	this->d->publishCadProgressiveCeiling(-1);
	this->d->lodSubmissionPass.retire();
	this->d->resetRetainedPassAnnotations();
	if (deadlineCapacityDecision.requestCeilingFreeFrame ||
	    deadlineCapacityDecision.requestSampleFrame) {
	    this->requestLodCapacityRender("lod-budget-sample");
	    return;
	}
	if (deadlineCapacityDecision.restartSubmission) {
	    this->restartLodCapacitySubmission(
		deadlineCapacityDecision.capacityCandidateChanged ? TRUE : FALSE);
	    return;
	}
	/* CERTIFIED and UNMEASURABLE are terminal for this unchanged capacity
	 * problem.  Retain the previously completed image until a real view,
	 * population, resource, or user-policy edge creates a new problem. */
	this->d->lodPointQualityPhase.completeCalibration();
	this->d->lodPointAdmissionFrame.retire();
	this->d->lastLodDiagnostics = "";
	/* Releasing the search owner must not erase an exact-frame obligation
	 * created by the ceiling change above.  Reconcile the complete work ledger
	 * so the retained image stays visible until its successor is presented. */
	this->synchronizeProgressiveWorkPending();
	return;
    }
    const bool exactProtectedFloorMiss =
	capacityRelevantDeadlineMiss && certifiedUnconstrainedCandidate &&
	allocationCertificate.selectsProtectedFloor() &&
	allocationCertificate.protectedFloorBudget ==
	    this->d->lodAdmissionEvidence.capacity().retainedQualityFloorBudget() &&
	allocationCertificate.protectedFloorIdentity ==
	    this->d->lodAdmissionEvidence.capacity().retainedQualityFloorIdentity();
    if (exactProtectedFloorMiss &&
	this->d->lodAdmissionEvidence.capacity().retainedQualityFloorActive()) {
	const size_t floorBudget =
	    this->d->lodAdmissionEvidence.capacity().retainedQualityFloorBudget();
	const uint64_t floorSignature =
	    this->d->lodAdmissionEvidence.capacity().retainedQualityFloorIdentity();
	const unsigned int missesBefore =
	    this->d->lodAdmissionEvidence.capacity().retainedQualityFloorMissCount();
	const bool rejected =
	    this->d->recordRetainedQualityFloorMiss();
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
	if (interruptedPresentedMaximum > 0) {
	    const int correctedCeiling = std::max(0,
		std::min(interruptedPresentedMaximum - 1,
		    interruptedCorrectedCeiling));
	    this->d->lodViewDemandPolicy.noteQualityMiss(
		correctedCeiling, interruptedEstimatedCost > 0 ?
		    interruptedEstimatedCost : state->activeRenderCost());
	    this->d->publishCadProgressiveCeiling(correctedCeiling);
	    if (!interactive &&
		deadlineSuccessor == BObolLodQualityPolicy::DeadlineSuccessor::
		    RECOVER_PRESENTATION) {
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
		this->d->rewindLodSubmissionCursor();
		this->d->lodSubmissionPass.beginFresh();
	    }
	}
    }

    /* Unbracketed wheel, trackpad, and resize interactions normally enter the
     * quiet state after their newest motion frame completes.  A renderer
     * deadline is equally decisive when that frame cannot complete: retaining
     * the completion latch otherwise leaves the controller permanently in
     * interactive mode and retries the same bounded frame forever.  Preserve
     * the debounce itself; this only lets its existing quiet-transition path
     * consume a hard capacity miss once input has been quiet long enough. */

    const int64_t deadlineObservationMicroseconds = bu_gettime();
    const int64_t quietDebounceMicroseconds =
	BObolLodViewDemandPolicy::unbracketedQuietDebounceMicroseconds();
    (void)this->d->lodInteractionSession.releaseExpiredMotionFrame(
	deadlineObservationMicroseconds, quietDebounceMicroseconds);
    if (this->d->lodInteractionSession.quietReady(
	    deadlineObservationMicroseconds, quietDebounceMicroseconds)) {
	/* The motion frame may have completed before the debounce expired, or
	 * this deadline miss may retire its still-standing frame gate.  In both
	 * states the quiet transition is now the sole successor.  Replaying the
	 * interactive traversal would preempt advanceLodProgressive() and keep a
	 * software view in deadline recovery forever.  Release only the replay
	 * latch; the request below wakes the established quiet transition and all
	 * queued source/result work remains intact. */
	this->d->lodInterruptedPresentationReplay.retire();
    }

    /* A static overscan miss is terminal evidence for this view/capacity
     * epoch.  Keep the richer occurrence prefix resident and the corrected
     * renderer ceiling above; do not let a later repaint reopen the same
     * failed quality staircase.  Keep the static allowance active only while
     * its handoff reconciles the retained occurrence cut.  Dropping to the
     * ordinary quiet allowance here discarded the last deadline-safe static
     * frame before that reconciliation could use its proof. */
    const bool rejectedStaticQualityTrial =
	BObolLodAdmissionPlanner::terminalCapacityMiss(
	    interactive != FALSE, cadDrawAttempted != FALSE,
	    transientPresentationRetry != FALSE,
	    this->d->lodStaticQualityTrial.inProgress(),
	    deadlinePopulationWorkPending);
    if (rejectedStaticQualityTrial) {
	BObolLodStaticQualityTrial::Constraint constraint;
	constraint.reason = BObolLodStaticQualityTrial::ConstraintReason::
	    PRESENTATION_DEADLINE;
	constraint.revisionStamp = this->d->admissionRevisionStamp();
	constraint.committedCeiling =
	    this->d->lodInteractiveProgressiveCeiling;
	constraint.candidateCeiling = interruptedPresentedMaximum;
	constraint.candidateCost = interruptedEstimatedCost > 0 ?
	    interruptedEstimatedCost : interruptedActiveCost;
	constraint.allowedCost = interruptedTargetCost;
	const bool trialRejected =
	    this->d->lodStaticQualityTrial.reject(constraint);
	/* A one-pixel point trial and the static-quality trial name the same
	 * candidate framebuffer.  Reject them atomically.  Leaving point
	 * calibration pending pauses the allocation cursor and requests the exact
	 * rejected frame forever, which a near-deadline 50k OSMesa population can
	 * expose nondeterministically.  Restore the point cut from the completed
	 * pre-trial framebuffer; the deadline handoff below may then reconcile
	 * triangle detail without another point-classification owner. */
	if (BObolLodAdmissionPlanner::
		rejectedOnePixelTrialReleasesCalibration(
		    trialRejected,
		    this->d->lodPointQualityPhase.presentationPending(),
		    this->d->lodPresentationPointProxyPixelThreshold)) {
	    const float retainedPointThreshold =
		this->d->lodStaticQualityTrial.
		    baselinePointProxyPixelThreshold();
	    this->d->publishCadPointProxyThreshold(retainedPointThreshold);
	    this->d->lodPointQualityPhase.completeCalibration();
	}
	/* Static overscan and the protected visual floor are distinct candidate
	 * populations.  This miss constrains the attempted overscan population;
	 * it is evidence against the floor only when the allocation certificate
	 * names that exact floor, which is handled by the population-specific miss
	 * accounting above.  Keeping the cheaper floor eligible prevents a
	 * transient rich-frame miss from collapsing prominent geometry to the
	 * ordinary cadence budget. */
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
	BObolLodAdmissionPlanner::pointDeadlineRequiresPopulationAggregation(
		interruptedActiveCost, interruptedMinimumCost,
		interruptedPresentedMaximum, interruptedCorrectedCeiling,
		interruptedTargetCost) ? TRUE : FALSE;
    const bool structuralPointAdmissionPending =
	this->d->lodPointAdmissionFrame.pending();
    const bool deadlinePointAggregationAllowed =
	controller_lod_adaptive_point_aggregation_allowed(this,
	    this->d->lodStaticQualityTrial.capacityRejected());
    if (state && state->hasCadPresentationAssemblies() &&
	cadDrawAttempted && !cadPreparationAdvanced && !interactive &&
	deadlineSuccessor == BObolLodQualityPolicy::DeadlineSuccessor::
	    RECOVER_PRESENTATION &&
	interruptedPopulationIrreducible &&
	deadlinePointAggregationAllowed &&
	!this->d->lodStructuralRepair.pointRelaxationPending() &&
	this->d->deadlinePointProxyAggregationApplicable() &&
	this->d->quietAllocationTargetFps() > 0.0f) {
	const BObolLodAdmissionPlan pressurePlan =
	    BObolLodAdmissionPlanner::planPointInterrupted(
		this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor,
		BObolLodPointCalibrationGoal::RESPONSIVE,
		this->d->lodPresentationPointProxyPixelThreshold,
		elapsedNanoseconds, this->d->pointCalibrationTargetFps(
		    BObolLodPointCalibrationGoal::RESPONSIVE));
	this->d->commitAdmissionPlan(pressurePlan);
	const BObolLodPointProxyEvidence::Decision &pressure =
	    pressurePlan.pointProxyDecision;
	if (pressure.changed) {
	    /* A coarser point cut changes only the presentation used to recover
	     * from the failed one-pixel trial.  It does not make that rejected
	     * one-pixel population cheaper.  Keep the rejection authoritative for
	     * this view/capacity epoch while the safe/unsafe point bracket settles;
	     * reopening the trial here produced an exact 1 -> coarse -> 1 cycle on
	     * 50k OSMesa wire views. */
	    this->d->publishCadPointProxyThreshold(pressure.threshold);
	    if (structuralPointAdmissionPending)
		this->requestStructuralPointAdmissionFrame(
		    "lod-deadline-structural-point-pressure");
	    else
		this->d->lodPointQualityPhase.requestCalibration();
	}
    }

    if (deadlineSuccessor == BObolLodQualityPolicy::DeadlineSuccessor::
	    RETRY_TRANSACTION) {
	/* Preserve every admission and calibration barrier until the first
	 * changed population has one unchanged replay.  A second non-preparation
	 * abort takes the ordinary bounded capacity-recovery path above. */
	this->requestLodCapacityRender(cadPreparationAdvanced ?
	    "render-preparation-replay" : "render-population-replay");
	return;
    }

    if (deadlineSuccessor == BObolLodQualityPolicy::DeadlineSuccessor::
	CONTINUE_POPULATION) {
	/* Preserve a genuine producer-owned continuation without restarting its
	 * finite 50k/150k cursor.  The pressure correction above may, however,
	 * have armed point calibration after deadlineSuccessor was classified.
	 * That calibration pauses the active cursor, so recompute ownership from
	 * the post-correction state.  Otherwise both transitions wait with only a
	 * pump level and no enabled producer or framebuffer edge. */
	this->d->lodInterruptedPresentationReplay.retire();
	this->markProgressiveWorkPending();
	BObolLodAdmissionPlanner::PointCalibrationProducerInputs
	    continuationInputs;
	continuationInputs.submissionPending =
	    this->d->lodSubmissionPass.active() != FALSE;
	continuationInputs.discoveryCalibrationPending =
	    this->d->lodPointAdmissionFrame.pending();
	continuationInputs.stableCalibrationPending =
	    this->d->lodPointQualityPhase.presentationPending();
	continuationInputs.capacityAllocationPending =
	    this->d->lodAdmissionEvidence.capacity().capacityAllocationPending();
	continuationInputs.capacitySamplePending =
	    this->d->lodAdmissionEvidence.capacity().presentationFramePending();
	continuationInputs.stablePresentationAvailable =
	    state && state->hasCadPresentationAssemblies();
	continuationInputs.providerPending =
	    this->d->lodAvailabilityLedger.providerPendingCount() > 0;
	continuationInputs.servicePending = deadlineServiceProducer;
	continuationInputs.publicationAwaitingFrameRequest =
	    this->d->lodPresentationTransaction.
		publicationAwaitingFrameRequest();
	const bool continuationProducer =
	    BObolLodAdmissionPlanner::pointProducerOwnsCalibrationFrame(
		continuationInputs);
	if (this->d->lodPointQualityPhase.presentationPending() &&
	    !continuationProducer)
	    this->requestLodCapacityRender("render-point-calibration");
	return;
    }

    /* An interrupted traversal requires another bounded policy/presentation
     * attempt, but it is not a persistent capacity sample.  Keep the barriers
     * armed while the renderer-only corrections above make that retry
     * cheaper; retained-cut admission remains based on completed frames. */
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_OVERLOAD);
    this->d->advanceAdmissionRevision(
        BObolLodAdmissionRevisionDomain::CAPACITY);
    this->markProgressiveWorkPending();
    this->requestLodCapacityRender("render-deadline");
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
BObolViewController::requestStructuralPointAdmissionFrame(const char *reason)
{
    BObolLodControlTransitionScope controlTransition(this);
    /* A structural classifier frame and a stable mesh-quality sample cannot
     * own the same presentation.  In particular, stable calibration may
     * pause an already installed CAD assembly while structural admission is
     * the producer which can replace its boxes.  Transfer ownership before
     * publishing the level-triggered render request, so there is always one
     * enabled successor and one frame consumer. */
    this->d->lodPointQualityPhase.completeCalibration();
    this->d->lodPointAdmissionFrame.request();
    this->requestLodCapacityRender(reason);
}

void
BObolViewController::armStableLodHeadroomProbeIfReady(void)
{
    BObolLodControlTransitionScope controlTransition(this);
    if (!this->automaticLodControlEnabled())
	return;

    /* Deterministic/offline convergence already owns an unbounded, pixel-
     * exact admission policy.  Static headroom probing is the interactive
     * policy's bounded attempt to exceed its preferred frame cadence; arming
     * that finite-policy transaction while terminal mode is active creates an
     * impossible successor.  The terminal allocator deliberately cannot
     * produce the finite retained-allocation certificate which a handoff
     * requires, so the cursor otherwise repeats an unchanged scene pass
     * forever.  Treat this as a mode-domain guard at the producer, rather than
     * repeatedly cancelling the invalid handoff in the consumer. */
    if (this->d->forceTerminalLodRefinement)
	return;

    /* A threshold mutation is not authoritative until the renderer has
     * classified one exact frame at that threshold.  Provider pumping may
     * continue behind this barrier, but neither another threshold change nor
     * mesh admission may consume the preceding classifier result.  This
     * includes a one-pixel trial armed while completing a coarse structural
     * seed: the histogram below still describes the coarse frame and must not
     * overwrite the successor trial before it reaches the renderer. */

    const bool saturatedPointCalibration =
	this->d->lodPointQualityPhase.presentationPending() &&
	BObolLodAdmissionPlanner::pointAtMaximumPixelThreshold(
	    this->d->lodPresentationPointProxyPixelThreshold);

    if (!BObolLodAdmissionPlanner::maySeedPointStructuralDistribution(
	    this->d->lodPointAdmissionFrame.pending(),
	    this->d->lodPointQualityPhase.presentationPending(),
	    this->d->lodPointQualityPhase.triangleRecoveryPending()) &&
	!saturatedPointCalibration)
	return;

    const BObolViewLodState *viewLodState = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    const size_t activePopulationCost = viewLodState ?
	viewLodState->activeRenderCost() : 0;
    size_t presentedRenderCost = 0;
    const bool presentedRenderCostObserved = viewLodState &&
	viewLodState->lastCadPresentedRenderCost(presentedRenderCost);
    const bool exactCompletedFrame = presentedRenderCostObserved &&
	presentedRenderCost > 0 &&
	viewLodState->lastCadPresentationFrameExact();
    /* Retained CAD work records outlive the host traversal which produced
     * them.  Point threshold alone is not a presentation identity: provider
     * publication or a renderer-route change can alter submitted work while
     * retaining that threshold.  Consume time only with the exact backend
     * cost reported by the same frame. */
    const uint64_t reusableCadSampleNanoseconds =
	presentedRenderCostObserved ?
	this->d->lodRendererPerformanceEvidence.cadPresentationNanosecondsAt(
	    this->d->lodPresentationPointProxyPixelThreshold,
	    presentedRenderCost) : 0;
    const uint64_t structuralSampleNanoseconds =
	presentedRenderCostObserved ?
	this->d->lodRendererPerformanceEvidence.
	    structuralPresentationNanosecondsAt(
		this->d->lodPresentationPointProxyPixelThreshold,
		presentedRenderCost) : 0;
    const uint64_t currentCadSampleNanoseconds =
	reusableCadSampleNanoseconds ? reusableCadSampleNanoseconds :
	structuralSampleNanoseconds;
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
    const size_t progressivePayloads = viewLodState ?
	viewLodState->cadProgressivePayloadCount() : 0;

    size_t presentedSubpixelOccurrences = 0;
    size_t presentedStructuralBoxes = 0;
    const bool exactOccurrenceCoverage = viewLodState &&
	viewLodState->lastCadPresentationOccurrenceCoverage(
	    presentedSubpixelOccurrences, presentedStructuralBoxes);
    const size_t terminalOccurrenceFailures = viewLodState ?
	viewLodState->cadOccurrenceTerminalFailureCount() : 0;

    const bool externalProducersSettled =
	this->d->lodAvailabilityLedger.providerPendingCount() == 0;
    const bool terminalAllocationPopulationSettled =
	BObolLodAvailabilityScheduler::allocationPopulationSettled(
	    externalProducersSettled,
	    activeTasks == 0 && queuedResults == 0,
	    this->d->lodAvailabilityLedger.resultsPending() == 0 &&
		!this->d->lodPresentationTransaction.publicationPending(),
	    this->d->lodAvailabilityLedger.residentGrowthPending(),
	    this->d->lodAvailabilityLedger.residencyDrainActive());

    /* Selective structural repair intentionally does not replay the complete
     * source census.  Once its exact, quiet framebuffer contains one useful
     * presentation (mesh, subpixel representative, or explicit terminal
     * failure) for every occurrence in the saved census, publish the stronger
     * coverage proof.  This enables resident suffix prefetch and the retained
     * allocator without manufacturing another source traversal. */
    const size_t usefulWithSubpixels = presentedSubpixelOccurrences >
	SIZE_MAX - activePayloads ? SIZE_MAX :
	activePayloads + presentedSubpixelOccurrences;
    const size_t usefulPresentationCount = terminalOccurrenceFailures >
	SIZE_MAX - usefulWithSubpixels ? SIZE_MAX :
	usefulWithSubpixels + terminalOccurrenceFailures;
    const bool rendererCoverageProofReady = exactOccurrenceCoverage &&
	presentedStructuralBoxes <= terminalOccurrenceFailures &&
	externalProducersSettled &&
	!this->d->lodPresentationTransaction.barrierPending() &&
	this->d->lodAvailabilityLedger.resultsPending() == 0 &&
	activeTasks == 0 && queuedResults == 0;
    if (rendererCoverageProofReady)
	(void)this->d->lodCoveragePolicy.confirmPresentedCoverage(
	    usefulPresentationCount);

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
    const size_t sceneBudget = this->d->lodAdmissionEvidence.capacity().currentBudget();
    const size_t maximumFirstWaveOccurrences =
	BObolLodAdmissionPlanner::structuralFirstWaveOccurrenceLimit(
	    sceneBudget);
    const std::vector<SoBRLDatabaseSource *> structuralSources =
	controller_render_database_source_roots(this);
    const bool meshFirstSceneSafe =
	controller_lod_mesh_first_scene_safe(structuralSources);
    const bool haveStructuralProjectionPopulation =
	exactStructuralProjection &&
	structuralProjection.visibleCount > terminalOccurrenceFailures;
    const bool pointAggregationOwnsStructuralFrontier =
	!meshFirstSceneSafe &&
	BObolLodAdmissionPlanner::pointAggregationApplicable(
	    exactStructuralProjection ? structuralProjection.visibleCount : 0) &&
	this->d->structuralPointAggregationRequired(
	    maximumFirstWaveOccurrences) &&
	maximumFirstWaveOccurrences != SIZE_MAX;

    /* A framebuffer histogram is exact for the occurrences installed when
     * that frame began, but it is not a whole-scene population proof while a
     * provider is still appending leaves.  Reseeding from each partial
     * population restarts the source cursor and turns parallel discovery into
     * a serial threshold/frame staircase.  The structural proxies already
     * provide immediate useful coverage during discovery; seed the first
     * mesh wave once, from the settled inventory which it is meant to bound. */
    if (haveStructuralProjectionPopulation &&
	externalProducersSettled &&
	!this->d->lodInteractionSession.active() &&
	pointAggregationOwnsStructuralFrontier) {
	const BObolLodAdmissionPlan seedPlan =
	    BObolLodAdmissionPlanner::planPointStructuralDistribution(
		this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor,
		this->d->lodPresentationPointProxyPixelThreshold,
		structuralProjection.cumulativeCount,
		structuralProjection.visibleCount,
		maximumFirstWaveOccurrences);
	this->d->commitAdmissionPlan(seedPlan);
	const BObolLodPointProxyEvidence::Decision &seed =
	    seedPlan.pointProxyDecision;
	if (seed.changed) {
	    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_POINT",
		    this->d->lodViewRevision.value()))
		bu_log("BObol LoD structural point seed stable=%.3f "
		       "discovery=%.3f next=%.3f visible=%zu limit=%zu\n",
		       this->d->lodPresentationPointProxyPixelThreshold,
		       this->d->lodDiscoveryPointProxyPixelThreshold,
		       seed.threshold, structuralProjection.visibleCount,
		       maximumFirstWaveOccurrences);
	    const float oldEffectivePointThreshold = std::max(
		this->d->lodPresentationPointProxyPixelThreshold,
		this->d->lodDiscoveryPointProxyPixelThreshold);
	    this->d->publishCadPointProxyThresholds(
		seed.threshold, BObolViewController::Impl::
		    POINT_PROXY_PIXEL_THRESHOLD_MINIMUM);
	    const float newEffectivePointThreshold = seed.threshold;
	    const bool classifierChanged = std::fabs(
		newEffectivePointThreshold - oldEffectivePointThreshold) >
		1.0e-6f;
	    /* The old cursor was admitted against a different point/mesh split.
	     * Keep resident results, but restart planning from the new exact
	     * presentation after its bounded classifier frame completes. */
	    this->d->rewindLodSubmissionCursor();
	    this->d->lodSubmissionDelta.reset();
	    this->d->advanceAdmissionRevision(
	        BObolLodAdmissionRevisionDomain::CAPACITY);
	    /* This frame classifies which structural occurrences may enter the
	     * first mesh wave.  It is an admission barrier, not a stable timing
	     * calibration: the scene owns no mesh population to time yet.  Once
	     * completeRenderTiming() consumes the exact classifier frame,
	     * armStableLodHeadroomProbeIfReady() routes its remaining box frontier
	     * into the importance-ordered repair pass. */
	    this->d->lodPointQualityPhase.completeCalibration();
	    if (classifierChanged) {
		this->requestStructuralPointAdmissionFrame(
		    "lod-structural-distribution-seed");
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
	this->d->publishCadDiscoveryPointProxyThreshold(
	    BObolViewController::Impl::POINT_PROXY_PIXEL_THRESHOLD_MINIMUM);
	const float newEffective =
	    this->d->lodPresentationPointProxyPixelThreshold;
	if (newEffective + 1.0e-6f < oldEffective) {
	    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_POINT",
		    this->d->lodViewRevision.value()))
		bu_log("BObol LoD discovery point release old=%.3f new=%.3f\n",
		       oldEffective, newEffective);
	    this->requestStructuralPointAdmissionFrame(
		"lod-discovery-point-release");
	    return;
	}
    }

    /* Structural repair owns an exact fallback frontier once point
     * aggregation has no useful successor.  This is true at the maximum
     * threshold, when every remaining proxy projects above that threshold, or
     * when the complete occurrence population already fits the scene
     * allowance and therefore should not be aggregated in the first place. */
    const bool maximumPointAggregation =
	BObolLodAdmissionPlanner::pointAtMaximumPixelThreshold(
	    this->d->lodPresentationPointProxyPixelThreshold);
    const size_t unresolvedStructuralCount =
	presentedStructuralBoxes > terminalOccurrenceFailures ?
	    presentedStructuralBoxes - terminalOccurrenceFailures : 0;
    const size_t unaggregatableStructuralCount = exactStructuralProjection ?
	BObolLodAdmissionPlanner::unaggregatableStructuralCount(
	    structuralProjection.cumulativeCount,
	    structuralProjection.visibleCount) : 0;
    bool pointPolicyExhausted = BObolLodAdmissionPlanner::
	pointPolicyExhaustedForStructuralFrontier(
	    maximumPointAggregation,
	    pointAggregationOwnsStructuralFrontier,
	    exactStructuralProjection, unresolvedStructuralCount,
	    unaggregatableStructuralCount);
    /* A retained-refinement observation is intentionally absent from these
     * readiness gates.  It is produced by planning against the incomplete
     * mesh population and is consumed when the exact structural owner starts
     * below.  Letting that subordinate observation block repair creates the
     * cycle structural -> allocation -> structural without changing a cut. */
    const bool structuralSuccessorBaseReady =
	unresolvedStructuralPresentation && exactStructuralProjection &&
	this->d->lodAutoSubmit && this->d->lodService &&
	!this->d->lodInteractionSession.active() &&
	this->d->lodCoveragePolicy.hasCompleteVisibleCount() &&
	!this->d->lodSubmissionPass.active() &&
	!this->d->lodPresentationTransaction.barrierPending() &&
	!this->d->lodPointQualityPhase.triangleRecoveryPending() &&
	BObolLodAvailabilityScheduler::structuralRepairMayOwn(
	    this->d->lodAvailabilityLedger.residentGrowthPending(),
	    this->d->lodAvailabilityLedger.residencyDrainActive()) &&
	externalProducersSettled &&
	this->d->lodAvailabilityLedger.resultsPending() == 0 &&
	activeTasks == 0 && queuedResults == 0;
    const bool pointCalibrationBlocksStructuralRepair =
	this->d->lodPointQualityPhase.presentationPending() &&
	!pointPolicyExhausted;
    const bool budgetCalibrationBlocksStructuralRepair =
	this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() &&
	!pointPolicyExhausted;

    /* The first-wave histogram bounds provider admission, but a completed
     * capacity search may prove that even that bounded tail does not fit.
     * With no producer, calibration, or frame owner left, advance the point
     * bracket from the exact unresolved frontier.  Previously this state had
     * no enabled transition: convergence remained below 100 percent with
     * thousands of boxes while every worker and control obligation was idle. */
    const bool pointSuccessorOwnerPending =
	pointCalibrationBlocksStructuralRepair ||
	budgetCalibrationBlocksStructuralRepair;
    if (structuralSuccessorBaseReady &&
	BObolLodAdmissionPlanner::
	    pointSuccessorRequiredForStructuralFrontier(
		exactStructuralProjection, unresolvedStructuralCount,
		pointPolicyExhausted, pointSuccessorOwnerPending)) {
	const BObolLodAdmissionPlan pointPlan =
	    BObolLodAdmissionPlanner::planPointStructuralCoverageBlocked(
		this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor,
		this->d->lodPresentationPointProxyPixelThreshold,
		unresolvedStructuralCount);
	this->d->commitAdmissionPlan(pointPlan);
	const BObolLodPointProxyEvidence::Decision &pointSuccessor =
	    pointPlan.pointProxyDecision;
	if (pointSuccessor.changed) {
	    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_POINT",
		    this->d->lodViewRevision.value()))
		bu_log("BObol LoD structural point successor current=%.3f "
		       "next=%.3f unresolved=%zu\n",
		       this->d->lodPresentationPointProxyPixelThreshold,
		       pointSuccessor.threshold, unresolvedStructuralCount);
	    this->d->publishCadPointProxyThreshold(pointSuccessor.threshold);
	    /* This threshold still classifies a structural-only population.  It
	     * must complete the admission-frame contract before source loading can
	     * use the new point/mesh split, but it cannot produce a stable mesh
	     * timing sample because no mesh has been admitted yet.  Assigning it to
	     * lodPointQualityPhase pauses the loader while waiting for evidence the
	     * paused loader alone can create. */
	    this->requestStructuralPointAdmissionFrame(
		"lod-structural-frontier-point-successor");
	    return;
	}
	/* A saturated point bracket has no further aggregation transition.  The
	 * same exact frontier now belongs to bounded structural repair. */
	pointPolicyExhausted = true;
    }
    const bool structuralRepairReady = structuralSuccessorBaseReady &&
	pointPolicyExhausted &&
	!budgetCalibrationBlocksStructuralRepair &&
	!pointCalibrationBlocksStructuralRepair;
    if (unresolvedStructuralPresentation && pointPolicyExhausted &&
	!structuralRepairReady) {
	SbString diagnostic;
	diagnostic.sprintf("structural repair waiting "
	    "(coverage=%d submit=%d budget=%d frame=%d retained=%d "
	    "producer=%d)",
	    this->d->lodCoveragePolicy.coverageComplete() ? 1 : 0,
	    this->d->lodSubmissionPass.active() ? 1 : 0,
	    this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() ? 1 : 0,
	    this->d->lodPresentationTransaction.barrierPending() ? 1 : 0,
	    this->d->lodRetainedPass.refinementPending() ? 1 : 0,
	    externalProducersSettled ? 0 : 1);
	if (this->d->lastLodDiagnostics.getLength() > 0)
	    this->d->lastLodDiagnostics += "\n";
	this->d->lastLodDiagnostics += diagnostic;
    }
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
	       "retained=%d handoff=%d point=%d point_exhausted=%d "
	       "recovery=%d results=%d "
	       "active_tasks=%zu queued_results=%zu\n",
	       unresolvedStructuralPresentation ? 1 : 0,
	       structuralRepairReady ? 1 : 0,
	       exactOccurrenceCoverage ? 1 : 0,
	       presentedStructuralBoxes, terminalOccurrenceFailures,
	       this->d->lodInteractionSession.active() ? 1 : 0,
	       this->d->lodInteractionSession.gestureActive() ? 1 : 0,
	       this->d->lodCoveragePolicy.coverageComplete() ? 1 : 0,
	       this->d->lodSubmissionPass.active() ? 1 : 0,
	       this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() ? 1 : 0,
	       this->d->lodPresentationTransaction.barrierPending() ? 1 : 0,
	       this->d->lodRetainedPass.refinementPending() ? 1 : 0,
	       this->d->lodPresentationPolicy.handoffPending() ? 1 : 0,
	       this->d->lodPointQualityPhase.presentationPending() ? 1 : 0,
	       pointPolicyExhausted ? 1 : 0,
	       this->d->lodPointQualityPhase.triangleRecoveryPending() ? 1 : 0,
	       this->d->lodAvailabilityLedger.resultsPending(), activeTasks, queuedResults);

    /* Once point aggregation has no useful successor, every remaining
     * uncollapsed box is a structural coverage obligation.  Transfer the
     * immediately preceding exact frame's hard-deadline headroom into one
     * structural admission transaction.  The policy rejects an unchanged
     * population after that attempt, preventing the preferred quiet budget
     * from replaying the same sparse frontier indefinitely. */
    size_t structuralPresentedCost = 0;
    const bool structuralPresentationFrameExact =
	exactOccurrenceCoverage && viewLodState &&
	viewLodState->lastCadPresentationFrameExact();
    const bool structuralRenderCostObserved = viewLodState &&
	viewLodState->lastCadPresentedRenderCost(structuralPresentedCost);
    const uint64_t structuralCadSampleNanoseconds =
	structuralRenderCostObserved ?
	    this->d->lodRendererPerformanceEvidence.
		cadPresentationNanosecondsAt(
		    this->d->lodPresentationPointProxyPixelThreshold,
		    structuralPresentedCost) : 0;
    const uint64_t structuralFallbackSampleNanoseconds =
	structuralRenderCostObserved ?
	    this->d->lodRendererPerformanceEvidence.
		structuralPresentationNanosecondsAt(
		    this->d->lodPresentationPointProxyPixelThreshold,
		    structuralPresentedCost) : 0;
    const uint64_t pairedStructuralSampleNanoseconds =
	structuralCadSampleNanoseconds ? structuralCadSampleNanoseconds :
	structuralFallbackSampleNanoseconds;
    const uint64_t structuralRepairDeadline =
	BObolLodAdmissionPlanner::structuralRepairPresentationDeadline(
	    meshFirstSceneSafe,
	    this->d->stablePresentationFrameDeadlineNanoseconds,
	    this->d->prominentQualityPresentationDeadline(),
	    controller_structural_repair_deadline_multiplier);
    const bool exactStructuralCapacityFrame =
	BObolLodAdmissionPlanner::structuralCapacityFrameApplicable(
	    structuralPresentationFrameExact,
	    structuralRenderCostObserved,
	    pairedStructuralSampleNanoseconds, structuralRepairDeadline);
    if (getenv("BOBOL_LOD_TRACE_STRUCTURAL_REPAIR") &&
	unresolvedStructuralPresentation)
	bu_log("BObol LoD structural frame exact=%d cost_observed=%d "
	       "cost=%zu render_ms=%.3f deadline_ms=%.3f applicable=%d\n",
	       structuralPresentationFrameExact ? 1 : 0,
	       structuralRenderCostObserved ? 1 : 0,
	       structuralPresentedCost,
	       pairedStructuralSampleNanoseconds / 1000000.0,
	       structuralRepairDeadline / 1000000.0,
	       exactStructuralCapacityFrame ? 1 : 0);
    size_t structuralCertifiedBudget = 0;
    if (exactStructuralCapacityFrame) {
	/* Structural placeholders are a time-to-first-image mechanism, never a
	 * stable-view representation.  A compact, fully profiled scene receives
	 * the existing one-shot prominent-quality allowance; an unbounded scene
	 * keeps the shorter finite-census allowance.  Admission and presentation
	 * consume the same deadline selected above. */
	const size_t presentationLimit =
	    BObolLodQualityPolicy::staticPresentationRenderCostLimit(
		structuralPresentedCost,
		pairedStructuralSampleNanoseconds,
		structuralRepairDeadline);
	const size_t additionalCapacity =
	    presentationLimit > structuralPresentedCost ?
		presentationLimit - structuralPresentedCost : 0;
	structuralCertifiedBudget =
	    additionalCapacity > SIZE_MAX - activePopulationCost ? SIZE_MAX :
		activePopulationCost + additionalCapacity;
	structuralCertifiedBudget = std::max(
	    structuralCertifiedBudget,
	    this->d->lodAdmissionEvidence.capacity().currentBudget());
    }
    /* The renderer's projection revision also changes when source metadata is
     * patched without changing the visible fallback population.  It is
     * therefore not a legal retry identity: keying admission to that revision
     * let a no-op source pass re-arm itself before the pending framebuffer was
     * presented.  Resolve and hash the actual occurrence frontier instead.
     * This same plan is transferred to submission below, so collection and
     * admission cannot disagree about which boxes the transaction owns. */
    ControllerStructuralSelection structuralRepairSelection;
    if (structuralRepairReady)
	structuralRepairSelection = controller_select_structural_occurrences(
	    this, viewLodState,
	    ControllerStructuralSelectionMode::UNCOLLAPSED);
    const bool exactStructuralFrontier =
	structuralRepairSelection.exact &&
	structuralRepairSelection.totalOccurrenceCount ==
	    presentedStructuralBoxes &&
	structuralRepairSelection.totalOccurrenceCount > 0;
    const uint64_t structuralFrontierIdentity =
	controller_structural_frontier_identity(
	    structuralProjection, structuralRepairSelection.plans);
    BObolLodStructuralAdmissionEvidence::Decision structuralAdmission;
    bool terminalProxyRepair = false;
    if (structuralRepairReady && exactStructuralCapacityFrame) {
	/* The exact structural repair transaction supersedes any point
	 * calibration or ordinary capacity sample which cannot change this
	 * frontier.  Consume those observations before installing the exclusive
	 * repair owner below. */
	if (this->d->lodPointQualityPhase.presentationPending())
	    this->d->lodPointQualityPhase.completeCalibration();
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::INVALIDATE_CAPACITY_MEASUREMENT);
	this->d->lodStructuralRepair.clearCoverageCostReservation();
	BObolLodStructuralAdmissionEvidence::Inputs inputs;
	inputs.viewRevision = this->d->lodViewRevision.value();
	inputs.policyRevision = this->d->lodPolicyRevision.value();
	inputs.frontierIdentity = structuralFrontierIdentity;
	inputs.unresolvedCount = unresolvedStructuralCount;
	inputs.activeCost = activePopulationCost;
	inputs.currentBudget = this->d->lodAdmissionEvidence.capacity().currentBudget();
	inputs.certifiedBudget = structuralCertifiedBudget;
	inputs.exactProjection = exactStructuralProjection;
	inputs.pointPolicyExhausted = pointPolicyExhausted;
	const BObolLodAdmissionPlan structuralPlan =
	    BObolLodAdmissionPlanner::planStructural(
		this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor, inputs);
	this->d->commitAdmissionPlan(structuralPlan);
	structuralAdmission = structuralPlan.structuralDecision;
	if (structuralAdmission.ownsFrontier) {
	    size_t admittedBudget = 0;
	    if (structuralAdmission.requestAdmission)
		admittedBudget = this->d->requestCoverageCompletion(
		    activePopulationCost, structuralAdmission.budget);
	    const size_t provisionalReservation =
		BObolLodAdmissionPlanner::structuralPerOccurrenceReservation(
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
		/* This exact, unchanged frontier has exhausted both point
		 * aggregation and mesh admission.  It still needs a legal terminal
		 * representation: convert the selected temporary AABBs to persistent
		 * OBBs in one bounded owner-thread pass.  An inexact selection cannot
		 * safely publish an occurrence-scoped constraint. */
		if (!exactStructuralFrontier)
		    return;
		terminalProxyRepair = true;
	    }
	}
    }
    if (structuralRepairReady) {
	/* The exact structural frontier is the higher-priority successor.  A
	 * richer-cut observation produced while discovering that frontier cannot
	 * survive into the closed repair transaction: it describes an incomplete
	 * occurrence population and otherwise re-arms ordinary allocation before
	 * repair can replace the boxes. */
	this->d->retireRetainedRefinementObservation();
	this->d->rewindLodSubmissionCursor();
	/* The completed framebuffer already performed the exact camera/frustum
	 * classification.  Route only its visible, uncollapsed structural
	 * occurrences back to the source loader.  A full compact scan here loads
	 * thousands of meshes whose subpixel point representations already satisfy
	 * the view and can turn a 293-box repair into 47k unnecessary payloads.
	 * Fall back to the complete scan if any retained identity cannot be mapped;
	 * incomplete selective repair would violate the no-box terminal contract. */
	if (exactStructuralFrontier)
	    this->d->lodSubmissionDelta.replaceSelectivePlans(
		std::move(structuralRepairSelection.plans));
	else
	    this->d->lodSubmissionDelta.reset();
	if (getenv("BOBOL_LOD_TRACE_STRUCTURAL_REPAIR"))
	bu_log("BObol LoD structural repair frontier boxes=%zu "
	       "mapped=%zu exact=%d sources=%zu\n",
	       presentedStructuralBoxes,
	       structuralRepairSelection.totalOccurrenceCount,
	       exactStructuralFrontier ? 1 : 0,
	       this->d->lodSubmissionDelta.planCount());
	this->d->lodSubmissionIntent.configure(true, false);
	if (terminalProxyRepair)
	    this->d->lodStructuralRepair.beginTerminalProxy(
		presentedStructuralBoxes - terminalOccurrenceFailures);
	else
	    this->d->lodStructuralRepair.begin(
		presentedStructuralBoxes - terminalOccurrenceFailures);
	/* This accumulator belongs to the newly closed repair transaction.  A
	 * prior coverage/demand pass may have exhausted its budget while the
	 * renderer was still classifying boxes; carrying that debt into this exact
	 * frontier spuriously re-entered point aggregation after every box had
	 * already been replaced. */
	this->d->lodRetainedPass.clearMissingMeshBudgetBlocked();
	this->d->lodSubmissionPass.beginFresh();
	this->d->lodRetainedPass.clearAdmittedWork();
	if (!terminalProxyRepair)
	    this->d->applyAdmissionEvidenceAction(
		BObolLodAdmissionPlanner::EvidenceAction::CLEAR_CAPACITY_LIMIT);
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::CAPACITY);
	this->markProgressiveWorkPending();
	this->requestLodCapacityRender("lod-structural-presentation-repair");
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
    const bool quietTerminalBoundaryContext = this->d->lodAutoSubmit &&
	this->d->lodService && !this->d->lodInteractionSession.active() &&
	!this->d->lodInteractionSession.gestureActive() && defaultQualityScale &&
	!this->d->lodForcedCut &&
	this->d->lodCoveragePolicy.coverageComplete() &&
	activePayloads > 0 &&
	!this->d->lodViewDemandPolicy.demandRefreshActive() &&
	!this->d->lodSubmissionPass.active() &&
	!this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() &&
	!this->d->lodPresentationTransaction.barrierPending() &&
	!this->d->lodRetainedPass.refinementPending() &&
	!this->d->lodPresentationPolicy.handoffPending() &&
	!this->d->lodPointQualityPhase.pending() &&
	terminalAllocationPopulationSettled;
    const bool stableTerminalContext = quietTerminalBoundaryContext &&
	memoryLimitedPayloads == 0;
    const bool stableSubpixelContext = stableTerminalContext &&
	activePayloads == satisfiedPayloads &&
	this->d->lodInteractiveProgressiveCeiling < 0;

    const auto schedulePointStructuralPreload =
	[this](ControllerStructuralSelection &&selection,
	    const char *renderReason) {
	    this->d->rewindLodSubmissionCursor();
	    this->d->lodSubmissionDelta.replaceSelectivePlans(
		std::move(selection.plans));
	    this->d->lodSubmissionIntent.configure(true, false);
	    this->d->lodRetainedPass.clearMissingMeshBudgetBlocked();
	    this->d->lodSubmissionPass.beginFresh();
	    this->d->lodRetainedPass.clearAdmittedWork();
	    this->markProgressiveWorkPending();
	    this->requestLodCapacityRender(renderReason);
	};

    /* Point relaxation is a two-phase presentation transaction.  The sparse
     * admission pass first replaces every structural occurrence which the
     * candidate threshold could reveal, while the renderer keeps drawing the
     * preceding coherent point image.  Only an exact post-publication frame
     * may commit the candidate threshold. */
    if (this->d->lodStructuralRepair.pointRelaxationPending() &&
	!this->d->lodStructuralRepair.active()) {
	const bool pointRelaxationFinalizeReady = this->d->lodAutoSubmit &&
	    this->d->lodService &&
	    !this->d->lodInteractionSession.active() &&
	    !this->d->lodInteractionSession.gestureActive() &&
	    this->d->lodCoveragePolicy.coverageComplete() &&
	    !this->d->lodSubmissionPass.active() &&
	    !this->d->lodPresentationTransaction.barrierPending() &&
	    !this->d->lodRetainedPass.refinementPending() &&
	    !this->d->lodPresentationPolicy.handoffPending() &&
	    externalProducersSettled &&
	    this->d->lodAvailabilityLedger.resultsPending() == 0 &&
	    activeTasks == 0 && queuedResults == 0;
	if (this->d->lodStructuralRepair.
		pointRelaxationPresentationPending() ||
	    !pointRelaxationFinalizeReady || !exactOccurrenceCoverage ||
	    !exactCompletedFrame)
	    return;

	const float relaxationTarget =
	    this->d->lodStructuralRepair.pointRelaxationTarget();
	ControllerStructuralSelection remaining =
	    controller_select_structural_occurrences(
		this, viewLodState,
		ControllerStructuralSelectionMode::ABOVE_PIXEL_THRESHOLD,
		relaxationTarget, maximumFirstWaveOccurrences);
	if (remaining.exact && remaining.totalOccurrenceCount == 0) {
	    this->d->lodStructuralRepair.reset();
	    this->d->publishCadPointProxyThreshold(relaxationTarget);
	    /* The candidate was private during preload.  Committing its point/mesh
	     * split is the semantic capacity edge: invalidate the old threshold's
	     * terminal certificate exactly here, after structural safety is proven,
	     * rather than while the preceding coherent image is still active. */
	    this->d->applyAdmissionEvidenceAction(
		BObolLodAdmissionPlanner::EvidenceAction::
		    INVALIDATE_CAPACITY_MEASUREMENT);
	    this->d->applyAdmissionEvidenceAction(
		BObolLodAdmissionPlanner::EvidenceAction::
		    CANCEL_HEADROOM_RETRY);
	    this->d->advanceAdmissionRevision(
		BObolLodAdmissionRevisionDomain::CAPACITY);
	    this->d->lodPointQualityPhase.requestStaticCalibration();
	    this->requestLodCapacityRender("lod-point-prefetch-handoff");
	    return;
	}
	const size_t precedingRemainingRank = this->d->lodStructuralRepair.
	    pointRelaxationRemainingRank();
	const bool retryAdmission = remaining.exact &&
	    this->d->lodStructuralRepair.retryPointRelaxationAdmission(
		remaining.plannedOccurrenceCount,
		remaining.totalOccurrenceCount);
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_POINT",
		this->d->lodViewRevision.value()))
	    bu_log("BObol LoD point preload finalize threshold=%.3f "
		   "preceding_rank=%zu remaining=%zu batch=%zu "
		   "exact=%d retry=%d\n",
		   relaxationTarget, precedingRemainingRank,
		   remaining.totalOccurrenceCount,
		   remaining.plannedOccurrenceCount,
		   remaining.exact ? 1 : 0,
		   retryAdmission ? 1 : 0);
	if (retryAdmission) {
	    schedulePointStructuralPreload(std::move(remaining),
		"lod-point-structural-prefetch-retry");
	    return;
	}

	/* Missing identity routing or a provider which could not replace its
	 * selected structural record makes the finer cut unsafe.  Preserve the
	 * current complete image and record a terminal capacity witness instead
	 * of exposing boxes or retrying an unchanged frontier indefinitely. */
	const BObolLodAdmissionPlan structuralLimitPlan =
	    BObolLodAdmissionPlanner::settlePointAtStructuralLimit(
		this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor,
		BObolLodPointCalibrationGoal::STATIC,
		this->d->lodPresentationPointProxyPixelThreshold);
	this->d->commitAdmissionPlan(structuralLimitPlan);
	this->d->lodStructuralRepair.reset();
	if (this->d->lastLodDiagnostics.getLength() > 0)
	    this->d->lastLodDiagnostics += "\n";
	this->d->lastLodDiagnostics += remaining.exact ?
	    "point relaxation retained its preceding cut because selected "
	    "structural geometry was unavailable" :
	    "point relaxation retained its preceding cut because structural "
	    "occurrence identity routing was incomplete";
	return;
    }

    /* Structural repair may consume the classifier frame which originally
     * carried a coarse point threshold, leaving a box-free presentation with
     * no quality-phase owner.  Whether cache hits or worker publication happen
     * to supply another frame must not decide final fidelity.  Require one
     * reusable sample of the exact coarse threshold; the point policy records
     * the settled witness and either selects a finer cut or proves this one is
     * terminal for the measured cadence. */
    const bool terminalPointWitnessRequired = exactCompletedFrame &&
	exactOccurrenceCoverage && presentedStructuralBoxes == 0 &&
	terminalOccurrenceFailures == 0 &&
	BObolLodAdmissionPlanner::pointTerminalReplayRequired(
	    this->d->lodAdmissionEvidence,
	    BObolLodPointCalibrationGoal::STATIC,
	    this->d->lodPresentationPointProxyPixelThreshold);

    BObolLodTerminalQualityReducer::Inputs terminalQualityInputs;
    terminalQualityInputs.stableTerminalContext = stableTerminalContext;
    terminalQualityInputs.exactCompletedFrame = exactCompletedFrame;
    terminalQualityInputs.progressiveCeiling =
	this->d->lodInteractiveProgressiveCeiling;
    terminalQualityInputs.progressivePayloadCount = progressivePayloads;
    terminalQualityInputs.staticPhaseActive =
	this->d->lodStaticQualityTrial.inProgress();
    terminalQualityInputs.staticPhaseRejected =
	this->d->lodStaticQualityTrial.capacityRejected();
    const BObolRetainedAllocationResult &terminalAllocation =
	this->d->lodRetainedAllocationCertificate;
    const bool terminalAllocationCurrent =
	this->d->retainedAllocationCertificateCurrent(viewLodState);
    const bool terminalAllocationPresentationRealized =
	this->d->retainedAllocationPresentationRealized(viewLodState);
    terminalQualityInputs.retainedCeilingConstraint =
	this->d->lodStaticQualityTrial.retainsRendererCeilingFor(
	    this->d->admissionRevisionStamp(),
	    this->d->lodInteractiveProgressiveCeiling);
    terminalQualityInputs.retainedAllocationReplacesCeiling =
	terminalAllocationCurrent && terminalAllocationPresentationRealized &&
	terminalAllocation.selectedPresentationCost > 0 &&
	terminalAllocation.selectedPresentationCost <=
	    terminalAllocation.certifiedPresentationBudget;
    terminalQualityInputs.structuralFrontierPending =
	unresolvedStructuralPresentation;
    terminalQualityInputs.pointPresentationRequired =
	terminalPointWitnessRequired;
    const bool terminalConstraintRecorded =
	BObolLodTerminalQualityReducer::allocationConstraintPresented(
	    this->d->lodAdmissionEvidence.capacity().stableBudgetLimited() ||
		this->d->lodStaticQualityTrial.blocksNewTrial(),
	    terminalAllocationCurrent,
	    terminalAllocationPresentationRealized);
    terminalQualityInputs.allocationPlanningRequired =
	BObolLodTerminalQualityReducer::allocationDemandNeedsPlanning(
	    progressivePayloads, terminalAllocationCurrent,
	    terminalAllocationPresentationRealized,
	    activePayloads > satisfiedPayloads,
	    terminalConstraintRecorded);
    terminalQualityInputs.allocationCapacityRequired =
	BObolLodTerminalQualityReducer::allocationDemandNeedsCapacity(
	    progressivePayloads, terminalAllocationCurrent,
	    terminalAllocationPresentationRealized,
	    terminalAllocation.selectedPresentationCost,
	    terminalAllocation.pixelDemandPresentationCost,
	    terminalAllocation.certifiedPresentationBudget,
	    terminalConstraintRecorded);
    const BObolLodTerminalQualityReducer::Decision terminalQualityDecision =
	BObolLodTerminalQualityReducer::reduce(terminalQualityInputs);

    if (terminalQualityDecision.owner ==
	    BObolLodTerminalQualityReducer::Owner::STATIC_QUALITY &&
	this->d->lodInteractiveProgressiveCeiling < 0) {
	/* A completed ceiling-free frame is the commit edge for the static
	 * trial.  Reaching this reducer means every producer and presentation
	 * barrier is already quiet.  The former code handled every other reducer
	 * owner but left this one in PROBING, which stranded an unsatisfied
	 * single-mesh allocation with no enabled successor.
	 *
	 * If the active population meets demand, no performance constraint is
	 * needed.  Otherwise the preceding static allocation and this exact frame
	 * jointly prove the terminal cut: keep its resident suffix, record the
	 * current-revision constraint, and let a future view/resource/policy edge
	 * reopen the problem. */
	if (activePayloads <= satisfiedPayloads) {
	    this->d->lodStaticQualityTrial.deactivate();
	    return;
	}

	BObolLodStaticQualityTrial::Constraint constraint;
	constraint.reason = BObolLodStaticQualityTrial::ConstraintReason::
	    ALLOCATION_SATURATED;
	constraint.revisionStamp = this->d->admissionRevisionStamp();
	constraint.committedCeiling = -1;
	constraint.committedCost = presentedRenderCost;
	const int activeMaximum = viewLodState ?
	    viewLodState->maximumActiveProgressiveCut() : -1;
	constraint.candidateCeiling = activeMaximum >= 0 ?
	    std::min(activeMaximum + 1,
		BOBOL_MESH_LOD_CUT_COUNT_MAX - 1) : -1;
	const size_t selectedCost =
	    terminalAllocation.selectedPresentationCost;
	const size_t minimumRicherCost = selectedCost == SIZE_MAX ? SIZE_MAX :
	    selectedCost + 1;
	constraint.candidateCost = std::max(
	    terminalAllocation.pixelDemandPresentationCost,
	    minimumRicherCost);
	constraint.allowedCost = std::max(
	    terminalAllocation.certifiedPresentationBudget,
	    this->d->lodAdmissionEvidence.capacity().currentBudget());
	if (!this->d->lodStaticQualityTrial.constrainPresented(constraint)) {
	    if (this->d->lastLodDiagnostics.getLength() > 0)
		this->d->lastLodDiagnostics += "\n";
	    this->d->lastLodDiagnostics +=
		"static-quality: could not record the exact terminal allocation constraint";
	    return;
	}
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
		this->d->lodViewRevision.value()))
	    bu_log("BObol LoD static allocation settled active=%zu "
		   "satisfied=%zu selected=%zu demand=%zu allowed=%zu\n",
		   activePayloads, satisfiedPayloads, selectedCost,
		   constraint.candidateCost, constraint.allowedCost);
	return;
    }

    if (quietTerminalBoundaryContext && memoryLimitedPayloads > 0 &&
	terminalPointWitnessRequired) {
	/* No finer point cut can be made coherent while its required mesh
	 * population is refused by the resident-memory boundary.  The current
	 * box-free frame is still a valid terminal presentation.  Record that
	 * resource witness once so READY remains reachable without pretending the
	 * missing population was affordable or replaying an unchanged frame. */
	const BObolLodAdmissionPlan structuralLimitPlan =
	    BObolLodAdmissionPlanner::settlePointAtStructuralLimit(
		this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor,
		BObolLodPointCalibrationGoal::STATIC,
		this->d->lodPresentationPointProxyPixelThreshold);
	this->d->commitAdmissionPlan(structuralLimitPlan);
	if (this->d->lastLodDiagnostics.getLength() > 0)
	    this->d->lastLodDiagnostics += "\n";
	this->d->lastLodDiagnostics +=
	    "point relaxation retained its current cut at the resident-memory boundary";
	return;
    }
    const bool terminalPointReplayRequired =
	terminalQualityDecision.owner ==
	    BObolLodTerminalQualityReducer::Owner::POINT_QUALITY;
    if (terminalPointReplayRequired) {
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_POINT",
		this->d->lodViewRevision.value())) {
	    static std::atomic<unsigned int> terminalPointTraceCount(0);
	    const unsigned int traceIndex =
		terminalPointTraceCount.fetch_add(1);
	    if (traceIndex < 64)
		bu_log("BObol LoD terminal point replay threshold=%.3f "
		       "active=%zu satisfied=%zu cost=%zu render_ms=%.3f\n",
		       this->d->lodPresentationPointProxyPixelThreshold,
		       activePayloads, satisfiedPayloads,
		       activePopulationCost,
		       currentCadSampleNanoseconds / 1000000.0);
	}
	const bool reusablePointSample =
	    this->d->lodRendererPerformanceEvidence.
		reusableCadPresentationAt(
		    this->d->lodPresentationPointProxyPixelThreshold) &&
	    currentCadSampleNanoseconds > 0;
	if (!reusablePointSample) {
	    this->d->lodPointQualityPhase.requestCalibration();
	    this->requestLodCapacityRender("lod-terminal-point-replay");
	    return;
	}

	/* Evaluate timing with a zero unresolved count deliberately: the
	 * resulting candidate is not exposed yet.  Its structural population is
	 * instead the sparse preload frontier installed below. */
	const BObolLodAdmissionPlan relaxationPlan =
	    BObolLodAdmissionPlanner::planPointCompleted(
		this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor,
		BObolLodPointCalibrationGoal::STATIC,
		this->d->lodPresentationPointProxyPixelThreshold,
		currentCadSampleNanoseconds,
		this->d->pointCalibrationTargetFps(
		    BObolLodPointCalibrationGoal::STATIC),
		true, true, 0);
	this->d->commitAdmissionPlan(relaxationPlan);
	const BObolLodPointProxyEvidence::Decision &relaxation =
	    relaxationPlan.pointProxyDecision;
	if (!relaxation.changed)
	    return;

	/* The projection inventory is bucketed at powers of two.  Every candidate
	 * within one bucket has the same conservative preload frontier, so commit
	 * its lower boundary directly rather than spending several unchanged
	 * calibration frames walking through that bucket. */
	const float preloadThreshold = structuralProjection.visibleCount > 0 ?
	    BObolLodAdmissionPlanner::pointStructuralPreloadThreshold(
		this->d->lodPresentationPointProxyPixelThreshold,
		relaxation.threshold) : relaxation.threshold;
	if (std::fabs(preloadThreshold -
		this->d->lodPresentationPointProxyPixelThreshold) <= 0.01f) {
	    /* A changed policy candidate must change the renderer control vector.
	     * Keep this runtime refinement guard even though the pure mapper is
	     * exhaustively tested: requesting the same threshold leaves the formal
	     * point-step rank unchanged and turns a finite calibration into an
	     * unbounded identical-frame loop. */
	    const BObolLodAdmissionPlan noProgressPlan =
		BObolLodAdmissionPlanner::settlePointAtStructuralLimit(
		    this->d->lodAdmissionEvidence,
		    this->d->lodAdmissionCursor,
		    BObolLodPointCalibrationGoal::STATIC,
		    this->d->lodPresentationPointProxyPixelThreshold);
	    this->d->commitAdmissionPlan(noProgressPlan);
	    if (this->d->lastLodDiagnostics.getLength() > 0)
		this->d->lastLodDiagnostics += "\n";
	    this->d->lastLodDiagnostics +=
		"point quality candidate retired because projection mapping "
		"did not change the presented threshold";
	    return;
	}
	ControllerStructuralSelection preload =
	    controller_select_structural_occurrences(
		this, viewLodState,
		ControllerStructuralSelectionMode::ABOVE_PIXEL_THRESHOLD,
		preloadThreshold, maximumFirstWaveOccurrences);
	if (!preload.exact) {
	    const BObolLodAdmissionPlan structuralLimitPlan =
		BObolLodAdmissionPlanner::settlePointAtStructuralLimit(
		    this->d->lodAdmissionEvidence,
		    this->d->lodAdmissionCursor,
		    BObolLodPointCalibrationGoal::STATIC,
		    this->d->lodPresentationPointProxyPixelThreshold);
	    this->d->commitAdmissionPlan(structuralLimitPlan);
	    if (this->d->lastLodDiagnostics.getLength() > 0)
		this->d->lastLodDiagnostics += "\n";
	    this->d->lastLodDiagnostics +=
		"point relaxation could not map its structural preload frontier";
	    return;
	}
	if (preload.totalOccurrenceCount == 0) {
	    this->d->publishCadPointProxyThreshold(preloadThreshold);
	    this->d->lodPointQualityPhase.requestStaticCalibration();
	    this->requestLodCapacityRender("lod-point-relaxation");
	    return;
	}
	const uint64_t staticPresentationDeadline =
	    this->d->pointCalibrationPresentationDeadline(
		BObolLodPointCalibrationGoal::STATIC);
	if (!BObolLodAdmissionPlanner::pointRelaxationFitsDeadline(
		activePayloads, preload.totalOccurrenceCount,
		currentCadSampleNanoseconds, staticPresentationDeadline)) {
	    const BObolLodAdmissionPlan structuralLimitPlan =
		BObolLodAdmissionPlanner::settlePointAtStructuralLimit(
		    this->d->lodAdmissionEvidence,
		    this->d->lodAdmissionCursor,
		    BObolLodPointCalibrationGoal::STATIC,
		    this->d->lodPresentationPointProxyPixelThreshold);
	    const bool terminalWitnessCommitted =
		this->d->commitAdmissionPlan(structuralLimitPlan);
	    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_POINT",
		    this->d->lodViewRevision.value())) {
		bu_log("BObol LoD point preload rejected threshold=%.3f "
		       "candidate=%.3f active=%zu additional=%zu "
		       "render_ms=%.3f deadline_ms=%.3f "
		       "witness_committed=%d replay_after=%d\n",
		       this->d->lodPresentationPointProxyPixelThreshold,
		       preloadThreshold, activePayloads,
		       preload.totalOccurrenceCount,
		       currentCadSampleNanoseconds / 1000000.0,
		       staticPresentationDeadline / 1000000.0,
		       terminalWitnessCommitted ? 1 : 0,
		       BObolLodAdmissionPlanner::pointTerminalReplayRequired(
			   this->d->lodAdmissionEvidence,
			   BObolLodPointCalibrationGoal::STATIC,
			   this->d->lodPresentationPointProxyPixelThreshold) ? 1 : 0);
	    }
	    return;
	}

	if (!this->d->lodStructuralRepair.beginPointRelaxation(
		preload.plannedOccurrenceCount,
		preload.totalOccurrenceCount, preloadThreshold)) {
	    const BObolLodAdmissionPlan invalidPlan =
		BObolLodAdmissionPlanner::settlePointAtStructuralLimit(
		    this->d->lodAdmissionEvidence,
		    this->d->lodAdmissionCursor,
		    BObolLodPointCalibrationGoal::STATIC,
		    this->d->lodPresentationPointProxyPixelThreshold);
	    this->d->commitAdmissionPlan(invalidPlan);
	    if (this->d->lastLodDiagnostics.getLength() > 0)
		this->d->lastLodDiagnostics += "\n";
	    this->d->lastLodDiagnostics +=
		"point relaxation produced an invalid bounded preload batch";
	    return;
	}
	/* Preload is deliberately invisible: the current threshold and its
	 * capacity certificate remain authoritative until the exact handoff above
	 * proves that lowering the threshold cannot reveal a box. */
	schedulePointStructuralPreload(std::move(preload),
	    "lod-point-structural-prefetch");
	return;
    }

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
	       this->d->lodSubmissionPass.active() ? 1 : 0,
	       this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() ? 1 : 0,
	       this->d->lodPresentationTransaction.barrierPending() ? 1 : 0,
	       this->d->lodRetainedPass.refinementPending() ? 1 : 0,
	       this->d->lodPresentationPolicy.handoffPending() ? 1 : 0,
	       this->d->lodPointQualityPhase.presentationPending() ? 1 : 0,
	       this->d->lodPointQualityPhase.triangleRecoveryPending() ? 1 : 0,
	       this->d->lodInteractiveProgressiveCeiling,
	       this->d->lodAvailabilityLedger.resultsPending(), activeTasks, queuedResults,
	       activePopulationCost, this->d->lodAdmissionEvidence.capacity().currentBudget(),
	       activePayloads, satisfiedPayloads, memoryLimitedPayloads,
	       residentMemoryHeadroom ? 1 : 0,
	       currentCadSampleNanoseconds / 1000000.0);

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
    if (terminalQualityDecision.owner ==
	    BObolLodTerminalQualityReducer::Owner::CEILING_RECONCILIATION) {
	const std::vector<SoBRLDatabaseSource *> sources =
	    controller_render_database_source_roots(this);
	if (!sources.empty()) {
	    const size_t cadCeilingCost = viewLodState->
		cadRenderCostAtProgressiveCutCeiling(
		    this->d->lodInteractiveProgressiveCeiling);
	    const size_t reconciliationBudget = std::max<size_t>(1,
		BObolLodAdmissionPlanner::canonicalSceneCostAtCadCeiling(
		    viewLodState->activeRenderCost(),
		    viewLodState->activeCadRenderCost(), cadCeilingCost));
	    this->d->lodPresentationPolicy.armHandoff(
		false, presentedRenderCost, reconciliationBudget);
	    this->d->requestPresentationReconciliation(
		reconciliationBudget);
	    this->d->rewindLodSubmissionCursor();
	    this->d->lodSubmissionPass.beginFresh();
	    this->d->resetRetainedPassAnnotations();
	    this->d->advanceAdmissionRevision(
	        BObolLodAdmissionRevisionDomain::CAPACITY);
	    /* Keep the terminal static-quality phase active while its occurrence
	     * allocation is translated from the renderer guard.  Releasing it here
	     * changes the next unchanged repaint back to the ordinary refinement
	     * deadline and coarsens a framebuffer which already met the explicit
	     * event-driven contract. */
	    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
		    this->d->lodViewRevision.value()))
		bu_log("BObol LoD terminal ceiling reconciliation "
		       "active=%d rejected=%d ceiling=%d budget=%zu\n",
		       this->d->lodStaticQualityTrial.inProgress() ? 1 : 0,
		       this->d->lodStaticQualityTrial.capacityRejected() ? 1 : 0,
		       this->d->lodInteractiveProgressiveCeiling,
		       reconciliationBudget);
	    this->markProgressiveWorkPending();
	    this->requestLodCapacityRender("lod-terminal-ceiling-reconciliation");
	    return;
	}
    }

    if (terminalQualityDecision.owner ==
	    BObolLodTerminalQualityReducer::Owner::ALLOCATION_PLANNING) {
	/* Point refinement invalidates the occurrence-allocation certificate by
	 * design.  Rebuild it once for the settled point population before asking
	 * the bounded capacity search to certify or constrain the remaining pixel
	 * demand.  Starting the complete scene cursor here gives this debt an
	 * executable owner instead of leaving an exact 99-percent view idle. */
	this->d->requestRetainedReallocation();
	this->beginSceneWideCapacitySubmission();
	this->markProgressiveWorkPending();
	this->requestLodCapacityRender("lod-terminal-allocation-planning");
	return;
    }

    if (terminalQualityDecision.owner ==
	    BObolLodTerminalQualityReducer::Owner::ALLOCATION_CAPACITY) {
	/* The exact occurrence allocation exhausted its certified allowance before
	 * reaching resident pixel demand.  Give that debt to the bounded capacity
	 * search.  Without this completed-frame edge, a zero-cut allocation pass
	 * can retire before its presentation becomes exact and leave no producer
	 * capable of either improving the cut or publishing a terminal constraint. */
	this->d->requestCapacityRescan();
	this->markProgressiveWorkPending();
	this->requestLodCapacityRender("lod-terminal-allocation-capacity");
	return;
    }

    /* Preserve only an exact, terminal presentation which actually met the
     * hard static-frame contract.  This evidence is keyed by the complete
     * camera/viewport/LoD signature and the controller's source-policy
     * domain.  It may seed a later return to this exact view, but the live
     * allocator and deadline recovery still own admission and correction. */
    if (stableTerminalContext && exactCompletedFrame &&
	this->d->lodViewSignature) {
	BObolLodViewQualityHistory::RememberInputs remembered;
	remembered.view = *this->d->lodViewSignature;
	remembered.domainRevision =
	    this->d->lodViewQualityDomainRevision;
	remembered.sceneAvailable = viewLodState && activePayloads > 0;
	remembered.quality.targetPixelError =
	    this->d->lodTargetPixelError;
	remembered.quality.progressiveCeiling =
	    this->d->lodInteractiveProgressiveCeiling;
	remembered.quality.progressiveNextFraction = viewLodState ?
	    viewLodState->cadPresentationProgressiveCutNextFraction() : 0.0f;
	remembered.quality.pointProxyPixelThreshold =
	    this->d->lodPresentationPointProxyPixelThreshold;
	remembered.quality.maximumProjectedErrorPixels =
	    this->d->lodRetainedAdmissionQualityEvidence.
		maximumProjectedErrorPixelsFor(
		    this->d->lodViewRevision.value(),
		    this->d->lodPolicyRevision.value(),
		    this->d->lodPresentationPointProxyPixelThreshold);
	remembered.quality.presentedRenderCost = presentedRenderCost;
	remembered.exactCompletedFrame = true;
	remembered.terminalPresentationComplete =
	    exactOccurrenceCoverage && presentedStructuralBoxes == 0 &&
	    terminalOccurrenceFailures == 0;
	remembered.producersSettled = externalProducersSettled;
	remembered.presentationDeadlineMet =
	    currentCadSampleNanoseconds > 0 &&
	    (!this->d->staticQualityPresentationDeadline() ||
	    currentCadSampleNanoseconds <=
		BObolLodViewDemandPolicy::staticHistoryDeadlineNanoseconds(
		    this->d->staticQualityPresentationDeadline()));
	remembered.resourcePressure =
	    this->d->lodAdmissionEvidence.resources().anyPressure();
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
	this->d->lodAdmissionEvidence.capacity().currentBudget();
    if (stableSubpixelContext && exactCompletedFrame &&
	currentCadSampleNanoseconds > 0) {
	const size_t demonstratedPresentationLimit =
	    BObolLodQualityPolicy::staticPresentationRenderCostLimit(
		presentedRenderCost, currentCadSampleNanoseconds,
		this->d->pointCalibrationPresentationDeadline(
		    BObolLodPointCalibrationGoal::STATIC));
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
	       this->d->lodAdmissionEvidence.capacity().currentBudget(), staticQualityBudget,
	       this->d->pointCalibrationPresentationDeadline(
		   BObolLodPointCalibrationGoal::STATIC) /
		   1000000.0,
	       currentCadSampleNanoseconds / 1000000.0);

    /* The first exact retained frame can still contain one-time VBO/atlas
     * preparation.  It establishes presentation correctness, but not the
     * reusable CAD timing sample needed to decide whether a finer pixel tier
     * fits.  In an event-driven viewport no later CAD traversal is implied by
     * that terminal-looking frame: HUD and screenshot repaints are explicitly
     * excluded from capacity evidence.  Give the existing point-quality owner
     * one unchanged replay so preparation/reuse accounting and asynchronous
     * GPU queries can advance.  The phase is consumed by completeRenderTiming
     * and cannot be rearmed once a threshold-matched reusable sample exists. */
    if (BObolLodQualityPolicy::staticQualityTimingReplayRequired(
	    stableSubpixelContext, exactCompletedFrame,
	    currentCadSampleNanoseconds)) {
	this->d->lodPointQualityPhase.requestCalibration();
	this->requestLodCapacityRender("lod-static-quality-timing-replay");
	return;
    }

    const float stablePixelError = stableSubpixelContext ?
	BObolLodQualityPolicy::stablePixelError(
	    this->d->lodTargetPixelError, activePopulationCost,
	    staticQualityBudget,
	    currentCadSampleNanoseconds,
	    this->d->pointCalibrationTargetFps(
		BObolLodPointCalibrationGoal::STATIC),
	    exactCompletedFrame, residentMemoryHeadroom) :
	this->d->lodTargetPixelError;
    if (stablePixelError + 1.0e-6f < this->d->lodTargetPixelError) {
	const size_t nextTierBudget =
	    BObolLodQualityPolicy::pixelErrorRenderCostFloor(
		activePopulationCost, this->d->lodTargetPixelError,
		stablePixelError);
	if (nextTierBudget <= staticQualityBudget)
	    this->d->raiseCapacityBudget(nextTierBudget);
	this->d->lodTargetPixelError = stablePixelError;
	this->advanceLodPolicyRevision(
	    LodPolicyTransition::CONTINUE_STATIC_QUALITY);
	this->d->lodLastSubmittedPolicyRevision.reset();
	this->d->rewindLodSubmissionCursor();
	this->d->lodSubmissionPass.beginFresh();
	this->markProgressiveWorkPending();
	this->requestLodCapacityRender("lod-stable-subpixel");
	return;
    }
    const int activeMaximumProgressiveCut = viewLodState ?
	viewLodState->maximumActiveProgressiveCut() : -1;
    const bool presentationQualityDebt = BObolLodAdmissionPlanner::
	progressivePresentationQualityDebt(activeMaximumProgressiveCut,
	    this->d->lodInteractiveProgressiveCeiling);
    const bool actionableQualityDebt = BObolLodAdmissionPlanner::
	actionableProgressiveQualityDebt(activePayloads, satisfiedPayloads,
	    memoryLimitedPayloads, activeMaximumProgressiveCut,
	    this->d->lodInteractiveProgressiveCeiling);

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
    /* The bounded capacity certificate already searches both the preferred
     * and event-driven static deadlines for a multi-occurrence allocation.
     * Re-entering the older renderer-ceiling trial afterward creates a second
     * owner for the same quality debt and can alternate a scene-wide coarse
     * prefix with the occurrence-local plan.  Keep the ordinal staircase only
     * for its valid domain: one progressive payload. */
    const bool staticOverscanEligible = BObolLodAdmissionPlanner::
	staticOrdinalOverscanApplicable(progressivePayloads) &&
	stableTerminalContext &&
	exactCompletedFrame && actionableQualityDebt &&
	!this->d->lodStaticQualityTrial.blocksNewTrial() &&
	!this->d->lodAdmissionEvidence.resources().anyPressure() &&
	this->d->pointCalibrationPresentationDeadline(
	    BObolLodPointCalibrationGoal::STATIC) >
	    preferredQuietNanoseconds;
    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
	    this->d->lodViewRevision.value()) &&
	(actionableQualityDebt || presentationQualityDebt ||
	 this->d->lodStaticQualityTrial.blocksNewTrial()))
	bu_log("BObol LoD static overscan eligible=%d terminal=%d exact=%d "
	       "actionable=%d presentation_debt=%d limited=%d active=%d "
	       "rejected=%d pressure=%d active_max=%d ceiling=%d "
	       "active_payloads=%zu satisfied=%zu memory_limited=%zu "
	       "deadline_ms=%.3f preferred_ms=%.3f\n",
	       staticOverscanEligible ? 1 : 0,
	       stableTerminalContext ? 1 : 0,
	       exactCompletedFrame ? 1 : 0,
	       actionableQualityDebt ? 1 : 0,
	       presentationQualityDebt ? 1 : 0,
	       this->d->lodAdmissionEvidence.capacity().stableBudgetLimited() ? 1 : 0,
	       this->d->lodStaticQualityTrial.inProgress() ? 1 : 0,
	       this->d->lodStaticQualityTrial.capacityRejected() ? 1 : 0,
	       this->d->lodAdmissionEvidence.resources().anyPressure() ? 1 : 0,
	       activeMaximumProgressiveCut,
	       this->d->lodInteractiveProgressiveCeiling,
	       activePayloads, satisfiedPayloads, memoryLimitedPayloads,
	       this->d->pointCalibrationPresentationDeadline(
		   BObolLodPointCalibrationGoal::STATIC) / 1000000.0,
	       preferredQuietNanoseconds / 1000000.0);
    if (staticOverscanEligible) {
	this->d->lodStaticQualityTrial.begin(
	    this->d->lodPresentationPointProxyPixelThreshold);
	BObolViewLodState *presentationState =
	    this->d->viewAttachment->getViewLodState();
	const int priorCeiling =
	    this->d->lodInteractiveProgressiveCeiling;
	const int activeMaximum = presentationState ? presentationState->
	    maximumActiveProgressiveCut() : -1;
	const int baselineCeiling =
	    BObolLodAdmissionPlanner::staticProgressiveBaselineCeiling(
		priorCeiling, activeMaximum, progressivePayloads);
	if (baselineCeiling >= 0) {
	    this->d->publishCadProgressiveCeiling(baselineCeiling);
	    this->d->publishCadCameraMotionFrameReuse(FALSE);
	}
	const int stagedCeiling =
	    BObolLodAdmissionPlanner::stagedProgressiveCeiling(
		priorCeiling,
		activeMaximum,
		this->d->lodPresentationPointProxyPixelThreshold,
		progressivePayloads);
	const bool retainAggregatedPresentation =
	    BObolLodAdmissionPlanner::retainAggregatedPresentation(
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
	    this->d->publishCadProgressiveCeiling(stagedCeiling);
	    this->d->publishCadCameraMotionFrameReuse(FALSE);
	    this->d->applyAdmissionEvidenceAction(
	        BObolLodAdmissionPlanner::EvidenceAction::CANCEL_HEADROOM_RETRY);
	    this->d->lodDiscretePopulationTrialPermit.grant();
	    this->markProgressiveWorkPending();
	    this->requestLodCapacityRender("lod-static-overscan-staged");
	    return;
	}
	/* A quiet handoff may have retained a responsive presentation ceiling.
	 * Multi-occurrence quality belongs to the allocator, so release that
	 * renderer-wide ordinal before its static pass.  A single occurrence keeps
	 * the exact baseline installed above while richer cuts are allocated behind
	 * it; completed frames will advance that ceiling through the bounded
	 * staircase. */
	if (!retainAggregatedPresentation && baselineCeiling < 0) {
	    this->d->publishCadProgressiveCeiling(-1);
	}
	if (presentationState)
	    this->d->publishCadCameraMotionFrameReuse(FALSE);
	/* The preferred-cadence point bracket is likewise a reversible
	 * presentation limit, not immutable geometry or a static fidelity policy.
	 * Test the one-pixel occurrence population directly under the interruptible
	 * hard static deadline.  If it is too expensive, the deadline handler keeps
	 * the preceding complete framebuffer and establishes the unsafe side of the
	 * point bracket in one observation.  Walking 64 -> 48 -> ... -> 1 after
	 * button-up is both slower and visibly distracting, and can strand a view at
	 * the ordinary quiet cut even when its full one-pixel frame fits inside
	 * the 100 ms event-driven allowance. */
	const SbBool restoredOnePixelPopulation =
	    !retainAggregatedPresentation &&
	    presentationState &&
	    presentationState->hasCadPresentationAssemblies() &&
	    BObolLodAdmissionPlanner::onePixelTrialRequired(
		this->d->lodPresentationPointProxyPixelThreshold) ?
		TRUE : FALSE;
	if (restoredOnePixelPopulation) {
	    this->d->publishCadPointProxyThreshold(
		BObolViewController::Impl::
		    POINT_PROXY_PIXEL_THRESHOLD_MINIMUM);
	    this->d->applyAdmissionEvidenceAction(
	        BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
	    this->d->lodPointQualityPhase.requestCalibration();
	}
	if (retainAggregatedPresentation)
	    this->d->lodPresentationPolicy.armHandoff(
		false, presentedRenderCost);
	else
	    this->d->lodPresentationPolicy.cancelHandoff();
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::CANCEL_HEADROOM_RETRY);
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::CLEAR_CAPACITY_LIMIT);
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::INVALIDATE_CAPACITY_MEASUREMENT);
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_OVERLOAD);
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::CAPACITY);
	this->d->requestRetainedReallocation(false);
	/* This is another pass in the same camera/policy/source epoch.  The
	 * explicit pending cursor is sufficient to bypass the completed-pass fast
	 * path.  Clearing the epoch witness as well makes the wrapper classify the
	 * already-pending cursor as a view change during submission and append an
	 * unnecessary full rescan after every allocation. */
	this->d->rewindLodSubmissionCursor();
	this->d->lodSubmissionPass.beginFresh();
	this->d->lodDiscretePopulationTrialPermit.grant();
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
		this->d->lodViewRevision.value()))
	    bu_log("BObol LoD static overscan armed prior_ceiling=%d "
		   "active_max=%d active_cost=%zu budget=%zu "
		   "target_fps=%.3f last_ms=%.3f\n",
		   priorCeiling,
		   activeMaximum,
		   activePopulationCost,
		   this->d->lodAdmissionEvidence.capacity().currentBudget(),
		   this->d->quietAllocationTargetFps(),
		   currentCadSampleNanoseconds / 1000000.0);
	this->markProgressiveWorkPending();
	this->requestLodCapacityRender("lod-static-overscan");
	return;
    }
    /* A terminal planning pass can finish while the motion-to-stable or
     * small-part presentation barrier still owns the next frame.  Arm one
     * exact, unchanged replay from either side of that barrier.  The
     * headroom policy makes the (view, policy, active population) witness
     * one-shot, so an actually capacity-limited population cannot repaint in
     * a loop. */
    const bool eligible = this->d->lodAutoSubmit && this->d->lodService &&
	!this->d->lodInteractionSession.active() &&
	BObolLodAdmissionPlanner::ordinaryHeadroomAllowed(
	    this->d->lodStaticQualityTrial.blocksNewTrial()) &&
	this->d->lodAdmissionEvidence.capacity().stableBudgetLimited() &&
	!this->d->lodSubmissionPass.active() &&
	!this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() &&
	!this->d->lodPresentationTransaction.barrierPending() &&
	!this->d->lodRetainedPass.refinementPending() &&
	!this->d->lodPresentationPolicy.handoffPending() &&
	!this->d->lodPointQualityPhase.pending() &&
	externalProducersSettled &&
	this->d->lodInteractiveProgressiveCeiling < 0 &&
	this->d->lodAdmissionEvidence.capacity().currentBudget() != SIZE_MAX &&
	actionableQualityDebt &&
	this->d->lodAvailabilityLedger.resultsPending() == 0 && activeTasks == 0 &&
	queuedResults == 0;
    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
	    this->d->lodViewRevision.value()))
	bu_log("BObol LoD headroom arm eligible=%d interactive=%d gesture=%d "
	       "limited=%d submit=%d rescan=%d refinement_frame=%d "
	       "retained=%d handoff=%d point=%d recovery=%d ceiling=%d "
	       "results=%d active_tasks=%zu queued_results=%zu active_cost=%zu "
	       "budget=%zu quality_debt=%d active_payloads=%zu satisfied=%zu "
	       "memory_limited=%zu pending=%d\n",
	       eligible ? 1 : 0, this->d->lodInteractionSession.active() ? 1 : 0,
	       this->d->lodInteractionSession.gestureActive() ? 1 : 0,
	       this->d->lodAdmissionEvidence.capacity().stableBudgetLimited() ? 1 : 0,
	       this->d->lodSubmissionPass.active() ? 1 : 0,
	       this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() ? 1 : 0,
	       this->d->lodPresentationTransaction.barrierPending() ? 1 : 0,
	       this->d->lodRetainedPass.refinementPending() ? 1 : 0,
	       this->d->lodPresentationPolicy.handoffPending() ? 1 : 0,
	       this->d->lodPointQualityPhase.presentationPending() ? 1 : 0,
	       this->d->lodPointQualityPhase.triangleRecoveryPending() ? 1 : 0,
	       this->d->lodInteractiveProgressiveCeiling,
	       this->d->lodAvailabilityLedger.resultsPending(), activeTasks, queuedResults,
	       activePopulationCost, this->d->lodAdmissionEvidence.capacity().currentBudget(),
	       actionableQualityDebt ? 1 : 0, activePayloads,
	       satisfiedPayloads, memoryLimitedPayloads,
	       this->d->lodAdmissionEvidence.headroom().retryPending() ? 1 : 0);
    if (!eligible)
	return;

    if (!activePopulationCost)
	return;
    const BObolLodAdmissionPlan headroomPlan =
	BObolLodAdmissionPlanner::planHeadroomRetry(
	    this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor, this->d->lodViewRevision,
	    this->d->lodPolicyRevision, activePopulationCost,
	    this->d->lodAdmissionEvidence.capacity().currentBudget());
    this->d->commitAdmissionPlan(headroomPlan);
    if (!headroomPlan.headroomAccepted)
	return;

    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
	    this->d->lodViewRevision.value()))
	bu_log("BObol LoD headroom armed view=%llu policy=%llu "
	       "active_cost=%zu budget=%zu\n",
	       (unsigned long long)this->d->lodViewRevision.value(),
	       (unsigned long long)this->d->lodPolicyRevision.value(),
	       activePopulationCost, this->d->lodAdmissionEvidence.capacity().currentBudget());

    this->markProgressiveWorkPending();
    this->requestLodCapacityRender("lod-calibrated-headroom-probe");
}

void
BObolViewController::finishLodRetainedRecovery(void)
{
    BObolLodControlTransitionScope controlTransition(this);
    if (!this->automaticLodControlEnabled())
	return;

    if (!this->d->confirmRetainedRecoveryPresentation(true))
	return;

    BObolViewLodState *presentationState = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    const BObolRetainedAllocationResult &allocation =
	this->d->lodRetainedAllocationCertificate;
    const bool allocationPresentationReady =
	this->d->retainedAllocationPresentationRealized(presentationState) &&
	allocation.certifiedPresentationBudget >=
	    allocation.selectedPresentationCost;
    if (allocationPresentationReady) {
	/* The recovery transaction already installed and presented the complete
	 * occurrence allocation which fits its certified allowance.  Remaining
	 * richer demand is quality debt beyond that allowance, not another
	 * recovery pass.  Re-admitting it here recomputes the scene budget from
	 * the just-completed frame, invalidates this certificate, and makes a
	 * large static scene alternate between the same two allocations. */
	this->d->retireRetainedRefinementObservation();
	return;
    }

    /* The recovery ceiling is a one-frame guard, not a fidelity policy.  A
     * changed classifier or population has made the prior allocation stale,
     * so rebuild it once from the coherent recovered presentation. */
    this->d->rewindLodSubmissionCursor();
    this->d->lodSubmissionPass.beginFresh();
    this->markProgressiveWorkPending();
}

SbBool
BObolViewController::completePointTriangleRecoveryIfReady(void)
{
    BObolLodControlTransitionScope controlTransition(this);
    if (!this->automaticLodControlEnabled() ||
	!this->d->lodPointQualityPhase.triangleRecoveryPending() ||
	this->d->lodInteractionSession.active() ||
	this->d->lodSubmissionPass.active() ||
	this->d->lodSubmissionPass.rescanPending() ||
	this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() ||
	this->d->lodPresentationPolicy.handoffPending() ||
	this->d->lodPresentationTransaction.barrierPending() ||
	this->d->lodPresentationTransaction.publicationPending())
	return FALSE;

    BObolViewLodState *presentationState = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    const BObolRetainedAllocationResult &recoveryAllocation =
	this->d->lodRetainedAllocationCertificate;
    const bool recoveryAllocationCurrent =
	this->d->retainedAllocationCertificateCurrent(presentationState);
    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_POINT",
	    this->d->lodViewRevision.value())) {
	static std::atomic<unsigned int> recoveryCompletionTraceCount(0);
	static constexpr unsigned int recoveryCompletionTraceLimit = 256;
	const unsigned int traceIndex =
	    recoveryCompletionTraceCount.fetch_add(1);
	if (traceIndex < recoveryCompletionTraceLimit)
	    bu_log("BObol LoD triangle recovery completion retained_pending=%d "
		   "allocation_current=%d plan=%llu active_plan=%llu\n",
		   this->d->lodRetainedPass.refinementPending() ? 1 : 0,
		   recoveryAllocationCurrent ? 1 : 0,
		   static_cast<unsigned long long>(
		       recoveryAllocation.allocationPlanSerial),
		   static_cast<unsigned long long>(presentationState ?
		       presentationState->activeCadAllocationPlan() : 0));
    }
    if (this->d->lodPointQualityPhase.recoveryAdmissionPending(
	    this->d->lodRetainedPass.refinementPending(),
	    recoveryAllocationCurrent))
	return FALSE;
    BObolViewLodState::CadStructuralProjectionHistogram structuralProjection;
    const SbBool exactStructuralProjection = presentationState &&
	presentationState->lastCadStructuralProjectionHistogram(
	    structuralProjection);
    const float recoveredPointThreshold =
	BObolLodAdmissionPlanner::triangleRecoveryPointThreshold(
	    this->d->lodPresentationPointProxyPixelThreshold,
	    exactStructuralProjection != FALSE,
	    exactStructuralProjection ? structuralProjection.visibleCount :
		SIZE_MAX);
    const SbBool pointCutChanged =
	!BObolLodRendererPerformanceEvidence::pointThresholdMatches(
	    recoveredPointThreshold,
	    this->d->lodPresentationPointProxyPixelThreshold) ? TRUE : FALSE;

    this->d->lodPointQualityPhase.completeTriangleRecovery(
	recoveryAllocationCurrent ? recoveryAllocation.allocationPlanSerial : 0);
    this->d->publishCadPresentationLimits(
	-1, 0.0f, recoveredPointThreshold);
    this->d->applyAdmissionEvidenceAction(
	BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
    if (pointCutChanged && presentationState &&
	presentationState->hasCadPresentationAssemblies())
	this->d->lodPointQualityPhase.requestCalibration();
    if (this->d->lodPointQualityPhase.presentationPending())
	this->requestLodCapacityRender("lod-stable-point-restore");
    else
	this->finishLodRetainedRecovery();
    return TRUE;
}

void
BObolViewController::completeRenderTiming(uint64_t startedNanoseconds,
	const BObolPresentationTimingContext &context)
{
    BObolLodControlTransitionScope controlTransition(this);
    const SbBool lodCapacityRelevant =
	context.lodCapacityRelevance == BObolLodCapacityRelevance::RELEVANT ?
	    TRUE : FALSE;
    const SbBool lodPlanningRelevant =
	context.lodPlanningRelevance == BObolLodPlanningRelevance::RELEVANT ?
	    TRUE : FALSE;
    const SbBool cadPresentationExecuted =
	context.cadPresentationExecution ==
	    BObolCadPresentationExecution::EXECUTED ? TRUE : FALSE;
    const bool automaticLod =
	this->automaticLodControlEnabled() != FALSE && this->d->lodService;

    const uint64_t now = this->beginRenderTiming();
    if (!startedNanoseconds || now <= startedNanoseconds)
	return;

    const uint64_t elapsed = now - startedNanoseconds;
    if (elapsed >= 30000000000ULL)
	return;
    BObolRenderClaimCompletionScope claimedRenderCompletion(this);

    if (context.cadPresentationCompleteness ==
	    BObolCadPresentationCompleteness::EXACT)
	(void)this->d->lodExactPresentationFrame.confirm(startedNanoseconds);

    this->d->lastRenderTimeNanoseconds = elapsed;

    this->d->smoothedRenderTimeNanoseconds =
	this->d->smoothedRenderTimeNanoseconds ?
	(this->d->smoothedRenderTimeNanoseconds * 9 + elapsed) / 10 : elapsed;
    /* The completed traversal is the commit edge for any interrupted CAD
     * replay.  Release queued immutable provider/result publications before
     * capacity-specific early returns so presentation-only classifier frames
     * cannot strand the owner-thread mutation gate. */
    this->d->lodInterruptedPresentationReplay.retire();

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
	    this->d->lodAdmissionEvidence.capacity().capacityTransactionPending()) {
	    const BObolLodCapacityEvidence::CompletedFrameDecision calibration =
		this->d->retireUnmeasurableCapacityFrame();
	    if (calibration.restartSubmission)
		this->restartLodCapacitySubmission(
		    calibration.capacityCandidateChanged ? TRUE : FALSE);
	}
	/* The initial point/box classifier deliberately contains no managed mesh
	 * population, so it is a presentation transaction rather than a geometry
	 * capacity sample.  Its exact frame must nevertheless retire the admission
	 * barrier.  Leaving that acknowledgement in the capacity-only path below
	 * makes a large structural scene repaint forever with zero submitted mesh
	 * tasks: the very absence of meshes prevents the latch which admits them
	 * from being consumed. */
	const BObolViewLodState *presentationState =
	    this->d->viewAttachment ?
	    this->d->viewAttachment->getViewLodState() : NULL;
	const bool admissionFrameWasPending =
	    this->d->lodPointAdmissionFrame.pending();
	const bool lodPresentationBarrierWasPending =
	    this->d->lodPresentationTransaction.barrierPending();
	size_t presentedSubpixelOccurrences = 0;
	size_t presentedStructuralBoxes = 0;
	const bool exactOccurrenceClassification = presentationState &&
	    presentationState->lastCadPresentationOccurrenceCoverage(
		presentedSubpixelOccurrences, presentedStructuralBoxes);
	const size_t terminalOccurrenceFailures = presentationState ?
	    presentationState->cadOccurrenceTerminalFailureCount() : 0;
	const bool nonterminalStructuralFrontier =
	    BObolLodPresentationPolicy::nonterminalStructuralFrontier(
		exactOccurrenceClassification, presentedStructuralBoxes,
		terminalOccurrenceFailures);
	if (admissionFrameWasPending) {
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
		this->d->lodPointAdmissionFrame.retire();
	    else
		this->requestLodCapacityRender("lod-admission-point-replay");
	}
	/* Presentation-only classifier frames deliberately bypass retained-mesh
	 * throughput calibration.  Preserve their exact duration separately so
	 * structural repair can derive a finite first-wave allowance instead of
	 * entering the generic zero-candidate reallocation loop. */
	if (lodPlanningRelevant && cadPresentationExecuted && presentationState &&
	    presentationState->lastCadPresentationFrameExact() &&
	    nonterminalStructuralFrontier) {
		size_t presentedRenderCost = 0;
		if (presentationState->lastCadPresentedRenderCost(
			presentedRenderCost))
		    this->d->lodRendererPerformanceEvidence.
			noteStructuralPresentation(elapsed,
			    this->d->lodPresentationPointProxyPixelThreshold,
			    presentedRenderCost);
	}
	this->completePresentationBarrier(elapsed,
	    context.cadPresentationCompleteness ==
		BObolCadPresentationCompleteness::EXACT ? TRUE : FALSE);
	/* An LoD-owned classifier/barrier frame may establish the last exact
	 * premise needed by the planning reducer, and cold discovery may expose a
	 * structural frontier without either latch.  Selection, highlighting, HUD,
	 * and other style-only frames are not implicit geometry successors: using
	 * their altered work identity to notice a missing CAD timing sample
	 * reopened a terminal OSMesa view during selection. */
	if (BObolLodPresentationPolicy::
		presentationOnlyFrameAdvancesPlanning(
		    lodPlanningRelevant != FALSE,
		    admissionFrameWasPending,
		    lodPresentationBarrierWasPending,
		    nonterminalStructuralFrontier))
	    this->armStableLodHeadroomProbeIfReady();
	/* Frame completion may either open a successor obligation or retire the
	 * last one.  Reconcile the complete level here instead of only setting it:
	 * otherwise a presentation-only selection repaint can leave the old pump
	 * bit standing after its exact frame commits.  Direct GL commonly receives
	 * an incidental later pump from Qt, but a synchronous software repaint has
	 * no such edge and then reports spurious background work indefinitely. */
	this->synchronizeProgressiveWorkPending();
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
    const BObolLodAdmissionPlan resourcePlan =
	BObolLodAdmissionPlanner::planResourceObservation(
	    this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor, cpuPressure != FALSE,
	    gpuPressure != FALSE,
	    automaticLod);
    this->d->commitAdmissionPlan(resourcePlan);
    const BObolLodResourceEvidence::Decision &resourceDecision =
	resourcePlan.resourceDecision;
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
	this->d->lodRendererPerformanceEvidence.acceptGpuMeasurement(
	    gpuCadSampleSerial, gpuCadNanoseconds);
    if (newGpuCadMeasurement) {
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
    if (this->d->lodStructuralRepair.pointRelaxationPending() &&
	this->d->lodStructuralRepair.
	    pointRelaxationPresentationPending() &&
	calibrationState && calibrationState->lastCadPresentationFrameExact() &&
	exactOccurrenceClassification)
	this->d->lodStructuralRepair.notePointRelaxationPresented();
    const SbBool admissionPointProxyFrameWasPending =
	this->d->lodPointAdmissionFrame.pending() ? TRUE : FALSE;
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
		       this->d->lodAvailabilityLedger.providerPendingCount(),
		       this->d->lodDiscoveryPointProxyPixelThreshold,
		       this->d->lodPresentationPointProxyPixelThreshold);
	}
	if (exactAdmissionClassification) {
	    this->d->lodPointAdmissionFrame.retire();
	} else {
	    /* A publication or scene transaction may make the requested frame
	     * inexact.  Preserve an explicit successor-frame witness; otherwise
	     * mesh admission can remain paused after the provider timer becomes
	     * idle. */
	    this->requestLodCapacityRender("lod-discovery-point-replay");
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
	BObolLodAdmissionPlanner::onePixelPresentationReady(
	    haveCadPresentationAssemblies != FALSE,
	    calibrationState &&
		calibrationState->lastCadPresentationFrameExact(),
	    exactOccurrenceClassification != FALSE,
	    presentedSubpixelOccurrences, presentedStructuralBoxes,
	    this->d->lodInteractiveProgressiveCeiling,
	    calibrationState ?
		calibrationState->maximumActiveProgressiveCut() : -1,
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
	    floorAllocation.revisionStamp.same(
		this->d->admissionRevisionStamp()) &&
	    std::isfinite(floorAllocation.pointProxyPixelThreshold) &&
	    std::fabs(floorAllocation.pointProxyPixelThreshold -
		this->d->lodPresentationPointProxyPixelThreshold) <= 1.0e-6f &&
	    calibrationState->cadAllocationPlanCoversCurrentPopulation(
		floorAllocation.allocationPlanSerial,
		floorAllocation.viewRevision(), floorAllocation.policyRevision(),
		floorAllocation.fixedCadPresentationCost);
	this->d->recordRetainedQualityFloorMet(
	    exactProtectedFloor, floorAllocation.protectedFloorIdentity,
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
    const bool preparedCadReplay = !measuredCadPresentation ||
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
    /* Prepared multi-draw replay is one reusable presentation path, but it
     * is not the only one.  OSMesa deliberately uses retained ordinary VBOs;
     * once their cumulative upload counters stop changing, an unchanged
     * traversal is an equally valid sustainable-throughput sample. */
    if (cadPresentationExecuted)
	this->d->lodRendererPerformanceEvidence.noteCadPresentation(
	    elapsed, this->d->lodPresentationPointProxyPixelThreshold,
	    measuredCadRenderCost ?
		std::optional<size_t>(presentedCadRenderCost) :
		std::optional<size_t>(),
	    preparedCadReplay,
	    haveFrameGpuResources ?
		std::optional<uint64_t>(geometryUploadBytes) :
		std::optional<uint64_t>());
    uint64_t calibrationElapsed =
	newGpuCadMeasurement ? gpuCadNanoseconds : elapsed;
    const SbBool gpuPointProxySampleMatchesCurrent =
	newGpuCadMeasurement &&
	BObolLodRendererPerformanceEvidence::pointThresholdMatches(
	    gpuCadPointProxyPixelThreshold,
	    this->d->lodPresentationPointProxyPixelThreshold) ?
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
	this->d->lodPresentationTransaction.publicationFramePending() &&
	!this->d->lodInteractionSession.active();
    const SbBool transientStablePreparation =
	(!cadPresentationExecuted ||
	 !this->d->lodRendererPerformanceEvidence.reusableCadPresentationAt(
	     this->d->lodPresentationPointProxyPixelThreshold)) &&
	!this->d->lodInteractionSession.active();
    const SbBool reusableCadWorkSample =
	cadPresentationExecuted &&
	this->d->lodRendererPerformanceEvidence.reusableCadPresentationAt(
	    this->d->lodPresentationPointProxyPixelThreshold) &&
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
		this->d->lodInteractionSession.active();
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
	     * hard bound.  Conflating the two remembered an earlier motion
	     * cut and discarded a visibly richer, deadline-safe frame as soon as
	     * wheel input became quiet.  The richer cut is already resident and
	     * was actually presented; retaining it does not teach the quiet
	     * allocator an inflated throughput. */
	    if (exactScaleQualityFrame && elapsed > 0 &&
		elapsed <= BObolLodViewDemandPolicy::
		    qualityFrameDurationNanoseconds())
		this->d->lodPresentationPolicy.noteStableQuality(
		    this->d->lodTargetPixelError,
		    presentedMaximum,
		    calibrationState ? calibrationState->
			cadPresentationProgressiveCutNextFraction() : 0.0f,
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
    if (this->d->lodAdmissionEvidence.headroom().retryPending()) {
	const bool stableContext =
	    !this->d->lodInteractionSession.active() &&
	    this->d->lodAdmissionEvidence.capacity().stableBudgetLimited() &&
	    !this->d->lodSubmissionPass.active() &&
	    !this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() &&
	    !this->d->lodPresentationTransaction.barrierPending() &&
	    this->d->quietAllocationTargetFps() > 0.0f &&
	    this->d->lodStableCalibratedRenderCostPerSecond > 0.0L &&
	    this->d->lodAdmissionEvidence.capacity().currentBudget() != SIZE_MAX;
	const long double demonstratedBudget = stableContext ?
	    this->d->lodStableCalibratedRenderCostPerSecond /
		static_cast<long double>(
		    this->d->quietAllocationTargetFps()) : 0.0L;
	const uint64_t stableTargetNanoseconds = stableContext ?
	    static_cast<uint64_t>(1000000000.0L /
		static_cast<long double>(
		    this->d->quietAllocationTargetFps())) : 0;
	const bool matchingHeadroomWitness =
	    this->d->lodAdmissionEvidence.headroom().pendingMatches(
		this->d->lodViewRevision, this->d->lodPolicyRevision,
		activeCalibrationCost);
	const size_t discreteTrialExcess =
	    BObolLodQualityPolicy::discreteTrialOverBudgetAllowance(
		activeCalibrationCost,
		this->d->lodAdmissionEvidence.capacity().currentBudget());
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
	bool deferredHeadroomReplay = false;
	if (transientHeadroomReplay) {
	    const BObolLodAdmissionPlan deferralPlan =
		BObolLodAdmissionPlanner::planHeadroomTransientDeferral(
		    this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor,
		    this->d->lodViewRevision, this->d->lodPolicyRevision,
		    activeCalibrationCost);
	    this->d->commitAdmissionPlan(deferralPlan);
	    deferredHeadroomReplay = deferralPlan.headroomAccepted;
	}
	bool reusableHeadroom = false;
	if (!deferredHeadroomReplay) {
	    const BObolLodAdmissionPlan consumptionPlan =
		BObolLodAdmissionPlanner::planHeadroomConsumption(
		    this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor,
		    this->d->lodViewRevision, this->d->lodPolicyRevision,
		    this->d->lodAdmissionEvidence.capacity().currentBudget(),
		    activeCalibrationCost, demonstratedBudget,
		    calibrationElapsed, stableTargetNanoseconds,
		    stableContext && reusableCadWorkSample &&
			!transientStablePublication &&
			!transientStablePreparation);
	    this->d->commitAdmissionPlan(consumptionPlan);
	    reusableHeadroom = consumptionPlan.headroomAccepted;
	}
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
		   this->d->lodAdmissionEvidence.capacity().currentBudget(),
		   static_cast<double>(demonstratedBudget),
		   static_cast<double>(calibrationElapsed) / 1000000.0,
		   static_cast<double>(stableTargetNanoseconds) / 1000000.0,
		   (unsigned long long)this->d->lodViewRevision.value(),
		   (unsigned long long)this->d->lodPolicyRevision.value());
	if (deferredHeadroomReplay) {
	    this->markProgressiveWorkPending();
	    this->requestLodCapacityRender("lod-calibrated-headroom-replay");
	} else if (reusableHeadroom || stableDiscreteTrial) {
	    this->d->applyAdmissionEvidenceAction(
	        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_CAPACITY_LIMIT);
	    this->d->rewindLodSubmissionCursor();
	    this->d->lodSubmissionPass.beginFresh();
	    this->d->advanceAdmissionRevision(
	        BObolLodAdmissionRevisionDomain::CAPACITY);
	    /* The replay proved a different sustainable allowance.  Recompute one
	     * complete screen-importance allocation with that new capacity even
	     * when the old active population is already below it.  An ordinary
	     * first-come pass can strand the unused remainder on a discrete PoP
	     * population and makes cold and warm caches converge differently. */
	    this->d->requestRetainedReallocation(false);
	    this->d->lodDiscretePopulationTrialPermit.setAvailable(
		stableDiscreteTrial);
	    this->markProgressiveWorkPending();
	    this->requestLodCapacityRender(stableDiscreteTrial ?
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
	this->d->lodPointQualityPhase.presentationPending()) {
	const std::vector<SoBRLDatabaseSource *> currentSources =
	    controller_render_database_source_roots(this);
	for (const SoBRLDatabaseSource *source : currentSources) {
	    if (source && source->getDisplayLodTargetCount() > 0) {
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
	!this->d->lodSubmissionPass.active() ? TRUE : FALSE;
    const bool pointCalibrationOwnsCompletedFrame =
	BObolLodAdmissionPlanner::pointCalibrationOwnsCompletedFrame(
	    this->d->lodPointQualityPhase.presentationPending(),
	    this->d->lodAdmissionEvidence.capacity().presentationFramePending());
    const bool pointCalibrationBlockedByGlobalCeiling =
	pointCalibrationOwnsCompletedFrame && haveCadPresentationAssemblies &&
	BObolLodTerminalQualityReducer::
	    globalCeilingBlocksPointPresentation(
		this->d->lodInteractiveProgressiveCeiling,
		calibrationState ?
		    calibrationState->cadProgressivePayloadCount() : 0);
    if (pointCalibrationOwnsCompletedFrame &&
	noVisibleCadPopulation) {
	this->d->lodPointQualityPhase.completeCalibration();
	this->d->publishCadPointProxyThreshold(
	    BObolViewController::Impl::POINT_PROXY_PIXEL_THRESHOLD_MINIMUM);
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
	(void)this->d->confirmRetainedRecoveryPresentation(true);
    }
    if (pointCalibrationOwnsCompletedFrame &&
	!haveCadPresentationAssemblies &&
	!haveCurrentCadSourceTargets &&
	!this->d->lodSubmissionPass.active()) {
	this->d->lodPointQualityPhase.completeCalibration();
	(void)this->d->confirmRetainedRecoveryPresentation(true);
    }
    if (pointCalibrationOwnsCompletedFrame &&
	!haveCadPresentationAssemblies && haveCurrentCadSourceTargets) {
	/* A path erase/expand transaction replaces the assembly occurrence table
	 * atomically.  Its source inventory already proves that a successor CAD
	 * presentation is expected, but this completed frame cannot calibrate the
	 * temporarily empty renderer record.  Preserve the obligation and wake the
	 * source producer which installs that presentation.  Requesting another
	 * empty render here while calibration also paused source admission formed a
	 * closed BALANCING/REFINING loop on 50k scenes. */
	if (!this->d->lodSubmissionPass.active()) {
	    this->d->rewindLodSubmissionCursor();
	    this->d->lodSubmissionPass.beginFresh();
	    this->d->advanceAdmissionRevision(
	        BObolLodAdmissionRevisionDomain::CAPACITY);
	}
	this->markProgressiveWorkPending();
    }
    if (pointCalibrationBlockedByGlobalCeiling) {
	/* This exact frame measured the renderer-wide ceiling, not the point
	 * population named by the pending threshold.  Retire the invalid point
	 * owner without changing its bracket.  The terminal-quality reducer below
	 * transfers the same completed-frame witness to ceiling reconciliation;
	 * rearming calibration here would keep both controls unchanged forever. */
	this->d->lodPointQualityPhase.completeCalibration();
	this->markProgressiveWorkPending();
    }
    if (pointCalibrationOwnsCompletedFrame &&
	haveCadPresentationAssemblies &&
	!pointCalibrationBlockedByGlobalCeiling) {
	const bool handoffConfirmation =
	    this->d->lodPointQualityPhase.handoffConfirmationPending();
	const bool staticPointCalibration =
	    this->d->lodPointQualityPhase.staticCalibrationPending();
	const BObolLodPointCalibrationGoal pointCalibrationGoal =
	    staticPointCalibration ? BObolLodPointCalibrationGoal::STATIC :
		BObolLodPointCalibrationGoal::RESPONSIVE;
	const float pointCalibrationTargetFps =
	    this->d->pointCalibrationTargetFps(pointCalibrationGoal);
	this->d->lodPointQualityPhase.completeCalibration();
	const long double stableTargetNanoseconds =
	    pointCalibrationTargetFps > 0.0f ?
		1000000000.0L /
		    static_cast<long double>(pointCalibrationTargetFps) : 0.0L;
	const auto requestPointCalibration = [this, staticPointCalibration]() {
	    if (staticPointCalibration)
		this->d->lodPointQualityPhase.requestStaticCalibration();
	    else
		this->d->lodPointQualityPhase.requestCalibration();
	};
	const SbBool reusableStableSample =
	    !this->d->lodInteractionSession.active() &&
	    !this->d->lodInteractionSession.gestureActive() &&
	    reusableCadWorkSample &&
	    !transientStablePublication &&
	    !transientStablePreparation &&
	    stableTargetNanoseconds > 0.0L &&
	    pointCalibrationElapsed > 0;
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_POINT",
		this->d->lodViewRevision.value())) {
	    static std::atomic<unsigned int> stablePointTraceCount(0);
	    const unsigned int traceIndex = stablePointTraceCount.fetch_add(1);
	    if (traceIndex < 64)
		bu_log("BObol LoD stable point frame threshold=%.3f "
		       "exact=%d occurrence=%d boxes=%zu structural=%zu "
		       "measured=%d render_cost=%d gpu=%d prepared=%d "
		       "reusable_path=%d transient_publication=%d "
		       "transient_preparation=%d reusable=%d elapsed_ms=%.3f\n",
		       this->d->lodPresentationPointProxyPixelThreshold,
		       calibrationState &&
			   calibrationState->lastCadPresentationFrameExact() ? 1 : 0,
		       exactOccurrenceClassification ? 1 : 0,
		       presentedStructuralBoxes, structuralFallbackOccurrences,
		       measuredCadPresentation ? 1 : 0,
		       measuredCadRenderCost ? 1 : 0,
		       newGpuCadMeasurement ? 1 : 0,
		       preparedCadReplay ? 1 : 0,
		   this->d->lodRendererPerformanceEvidence.
			   reusableCadPresentationAt(
			       this->d->lodPresentationPointProxyPixelThreshold) ?
			   1 : 0,
		       transientStablePublication ? 1 : 0,
		       transientStablePreparation ? 1 : 0,
		       reusableStableSample ? 1 : 0,
		       pointCalibrationElapsed / 1000000.0);
	}
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
	/* The completed frame does not make an active submission cursor an
	 * independent producer.  This callback has just released the current
	 * point-calibration phase, but any successor calibration requested below
	 * immediately freezes that cursor again.  Preserve the pause sampled at
	 * entry rather than asking the already-completed phase for its current
	 * state; otherwise calibration can wait for a cursor which is waiting for
	 * calibration, with neither side owning a host frame request. */
	BObolLodAdmissionPlanner::PointCalibrationProducerInputs producerInputs;
	producerInputs.submissionPending =
	    this->d->lodSubmissionPass.active() != FALSE;
	producerInputs.discoveryCalibrationPending =
	    this->d->lodPointAdmissionFrame.pending();
    producerInputs.stableCalibrationPending =
	    pointCalibrationOwnsCompletedFrame;
    producerInputs.capacityAllocationPending =
	this->d->lodAdmissionEvidence.capacity().capacityAllocationPending();
    producerInputs.capacitySamplePending =
	    this->d->lodAdmissionEvidence.capacity().presentationFramePending();
	producerInputs.stablePresentationAvailable =
	    haveCadPresentationAssemblies != FALSE;
	producerInputs.providerPending =
	    this->d->lodAvailabilityLedger.providerPendingCount() > 0;
	producerInputs.servicePending = pointCalibrationServiceWork;
	producerInputs.publicationAwaitingFrameRequest =
	    this->d->lodPresentationTransaction.
		publicationAwaitingFrameRequest();
	const SbBool pointCalibrationProducerWork =
	    BObolLodAdmissionPlanner::pointProducerOwnsCalibrationFrame(
		producerInputs) ? TRUE : FALSE;
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
	    !handoffConfirmation && BObolLodAdmissionPlanner::
		startOnePixelTrialFromSettledPointFrame(
		    this->d->lodInteractionSession.active(),
		    exactReusablePointFrame != FALSE,
		    pointCalibrationProducerWork != FALSE,
		    this->d->lodAdmissionEvidence.resources().anyPressure(),
		    structuralFallbackOccurrences != 0,
		    this->d->lodStaticQualityTrial.blocksNewTrial(),
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
	    !handoffConfirmation &&
	    BObolLodAdmissionPlanner::acceptSettledOnePixelFrame(
		this->d->lodInteractionSession.active(),
		exactReusablePointFrame != FALSE,
		pointCalibrationProducerWork != FALSE,
		this->d->lodAdmissionEvidence.resources().anyPressure(),
		structuralFallbackOccurrences != 0,
		this->d->lodPresentationPolicy.handoffPending(),
		this->d->lodStaticQualityTrial.blocksNewTrial(),
		this->d->lodInteractiveProgressiveCeiling,
		this->d->lodPresentationPointProxyPixelThreshold,
		reducibleProgressiveDetail != FALSE,
		pointCalibrationElapsed, preferredStableNanoseconds,
		this->d->stablePresentationFrameDeadlineNanoseconds) ?
		TRUE : FALSE;
	const BObolRetainedAllocationResult &pointAllocation =
	    this->d->lodRetainedAllocationCertificate;
	const bool triangleAllocationSaturated = calibrationState &&
	    this->d->lodPointQualityPhase.triangleRecoverySaturatedBy(
		pointAllocation.allocationPlanSerial) &&
	    pointAllocation.revisionStamp.same(
		this->d->admissionRevisionStamp()) &&
	    calibrationState->cadAllocationPlanCoversCurrentPopulation(
		pointAllocation.allocationPlanSerial,
		pointAllocation.viewRevision(), pointAllocation.policyRevision(),
		pointAllocation.fixedCadPresentationCost);
	SbBool scheduledTriangleRecovery = FALSE;
	SbBool restoredOnePixelCut =
	    startStaticOnePixelTrial || acceptSettledOnePixelFrame;

	if (acceptSettledOnePixelFrame) {
	    /* The current framebuffer is already the exact one-pixel static
	     * witness.  Preserve it while one occurrence-local minimax pass spends
	     * the capacity extrapolated to the hard deadline.  No point threshold,
	     * renderer ceiling, resident suffix, or visible cut is reset. */
	    this->d->lodStaticQualityTrial.begin(
		this->d->lodPresentationPointProxyPixelThreshold);
	    this->d->lodDiscretePopulationTrialPermit.grant();
	    const size_t staticQualityBudget =
		BObolLodQualityPolicy::staticPresentationRenderCostLimit(
		    activeCalibrationCost, pointCalibrationElapsed,
		    this->d->staticQualityPresentationDeadline());
	    this->d->applyAdmissionEvidenceAction(
	        BObolLodAdmissionPlanner::EvidenceAction::CANCEL_HEADROOM_RETRY);
	    this->d->applyAdmissionEvidenceAction(
	        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_CAPACITY_LIMIT);
	    this->d->applyAdmissionEvidenceAction(
	        BObolLodAdmissionPlanner::EvidenceAction::INVALIDATE_CAPACITY_MEASUREMENT);
	    this->d->applyAdmissionEvidenceAction(
	        BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_OVERLOAD);
	    this->d->advanceAdmissionRevision(
	        BObolLodAdmissionRevisionDomain::CAPACITY);
	    if (staticQualityBudget > 0)
		this->d->raiseCapacityBudget(
		    staticQualityBudget);
	    this->d->requestRetainedReallocation(false);
	    this->d->rewindLodSubmissionCursor();
	    this->d->lodSubmissionPass.beginFresh();
	    this->d->resetRetainedPassAnnotations();
	    this->markProgressiveWorkPending();
	}

	if (startStaticOnePixelTrial) {
	    /* Preserve every retained occurrence cut and perform one bounded
	     * presentation-only population trial.  The endpoint deadline keeps the
	     * preceding complete framebuffer visible if this richer classification
	     * is too expensive; its miss establishes the unsafe side of the point
	     * bracket without a coarse triangle-allocation round trip. */
	    this->d->lodStaticQualityTrial.begin(
		this->d->lodPresentationPointProxyPixelThreshold);
	    this->d->lodDiscretePopulationTrialPermit.grant();
	    this->d->publishCadPointProxyThreshold(
		BObolViewController::Impl::
		    POINT_PROXY_PIXEL_THRESHOLD_MINIMUM);
	    this->d->applyAdmissionEvidenceAction(
	        BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
	    requestPointCalibration();
	    this->requestLodCapacityRender("lod-static-point-one-pixel");
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
	if (!handoffConfirmation && !restoredOnePixelCut &&
	    reusableStableSample &&
	    !unresolvedStructuralOccurrences &&
	    BObolLodAdmissionPlanner::shouldRecoverTriangleDetail(
		reducibleProgressiveDetail != FALSE,
		stableSampleOverloaded != FALSE,
		coarsePointCut != FALSE,
		this->d->lodAdmissionEvidence.capacity().retainedQualityFloorActive() ||
		    this->d->lodAdmissionEvidence.capacity().retainedQualityFloorRejected() ||
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
		this->d->requestRetainedRecovery(
		    recoveryBudget);
		this->d->applyAdmissionEvidenceAction(
		    BObolLodAdmissionPlanner::EvidenceAction::INVALIDATE_CAPACITY_MEASUREMENT);
		this->d->applyAdmissionEvidenceAction(
		    BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_OVERLOAD);
		this->d->advanceAdmissionRevision(
		    BObolLodAdmissionRevisionDomain::CAPACITY);
		this->d->applyAdmissionEvidenceAction(
		    BObolLodAdmissionPlanner::EvidenceAction::CANCEL_HEADROOM_RETRY);
		/* The handoff's cost floor was proven with the coarse point cut
		 * which hid this population.  It is not evidence for the one-pixel
		 * hierarchy and would immediately re-admit the detail being
		 * compacted here.  Retain the renderer ceiling itself until the
		 * coherent recovery frame below, but retire the stale proof. */
		this->d->lodPresentationPolicy.cancelHandoff();
		this->d->rewindLodSubmissionCursor();
		this->d->lodSubmissionPass.beginFresh();
		this->d->retireRetainedRefinementObservation();
		this->d->lodRetainedPass.clearAdmittedWork();
		this->d->lodPointQualityPhase.beginTriangleRecovery();
		this->markProgressiveWorkPending();
		scheduledTriangleRecovery = TRUE;
	    } else if (coarsePointCut &&
		structuralFallbackOccurrences == 0) {
		/* The complete retained cut already fits the demonstrated capacity.
		 * Test the one-pixel population directly; no provider or admission
		 * traversal is necessary. */
		this->d->publishCadPointProxyThreshold(
		    BObolViewController::Impl::
			POINT_PROXY_PIXEL_THRESHOLD_MINIMUM);
		this->d->applyAdmissionEvidenceAction(
		    BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
		requestPointCalibration();
		this->requestLodCapacityRender("lod-stable-point-restore");
		restoredOnePixelCut = TRUE;
	    }
	}

	BObolLodPointProxyEvidence::Decision calibration;
	if (!handoffConfirmation && !scheduledTriangleRecovery &&
	    !restoredOnePixelCut) {
	    const bool pointCoarseningAllowed =
		controller_lod_adaptive_point_aggregation_allowed(this,
		    this->d->lodStaticQualityTrial.capacityRejected());
	    const BObolLodAdmissionPlan calibrationPlan =
		BObolLodAdmissionPlanner::planPointCompleted(
		    this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor,
		    pointCalibrationGoal,
		    this->d->lodPresentationPointProxyPixelThreshold,
		    pointCalibrationElapsed,
		    pointCalibrationTargetFps,
		    reusableStableSample != FALSE, pointCoarseningAllowed,
		    structuralFallbackOccurrences);
	    this->d->commitAdmissionPlan(calibrationPlan);
	    calibration = calibrationPlan.pointProxyDecision;
	}
	if (!scheduledTriangleRecovery && !restoredOnePixelCut &&
	    calibration.changed) {
	    this->d->publishCadPointProxyThreshold(calibration.threshold);
	    requestPointCalibration();
	    this->requestLodCapacityRender("lod-stable-point-calibration");
	} else if (!handoffConfirmation && !scheduledTriangleRecovery &&
	    !restoredOnePixelCut &&
	    !reusableStableSample &&
	    !this->d->lodInteractionSession.active() &&
	    !this->d->lodInteractionSession.gestureActive() &&
	    haveCadPresentationAssemblies &&
	    stableTargetNanoseconds > 0.0L &&
	    BObolLodAdmissionPlanner::pointRequiresReusableConfirmation(
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
	    requestPointCalibration();
	    if (pointCalibrationProducerWork) {
		this->markProgressiveWorkPending();
	    } else {
		this->requestLodCapacityRender("lod-stable-point-replay");
	    }
	}
	/* A triangle-prefix recovery which temporarily used a coarser aggregate
	 * point cut reaches this branch twice: once to request the one-pixel frame,
	 * and once after that frame completes.  The latter exact presentation is
	 * the missing transition which retires the measured recovery ceiling. */
	if (!scheduledTriangleRecovery && !restoredOnePixelCut &&
	    exactRetainedPopulation &&
	    this->d->lodPresentationPointProxyPixelThreshold <= 1.01f)
	    this->finishLodRetainedRecovery();
    }

	/* A retained-prefix recovery is complete only after its selected cut has
	 * been presented and no follow-up admission/handoff barrier remains.
	 * Restore one pixel only for a fully realized mesh population.  If
	 * structural records remain, retain their proven aggregate cut: lowering
	 * it would expose boxes rather than detail. */
    (void)this->completePointTriangleRecoveryIfReady();

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
	BObolLodAdmissionPlanner::pointStreamingPopulationWorkPending(
		this->d->lodSubmissionPass.active() != FALSE,
		this->hasPendingLodResults(), ongoingServiceWork,
		this->d->lodAvailabilityLedger.providerPendingCount()) ? TRUE : FALSE;
	if (!admissionPointProxyFrameWasPending &&
	    !this->d->lodPointAdmissionFrame.pending() &&
	!this->d->lodPointQualityPhase.pending() &&
	haveCadPresentationAssemblies &&
	this->d->pointProxyAggregationApplicable() &&
	this->d->structuralPointAggregationRequiredForBudget(
	    this->d->lodAdmissionEvidence.capacity().currentBudget()) &&
	!this->d->lodInteractionSession.active() &&
	!this->d->lodInteractionSession.gestureActive() &&
	automaticLod &&
	!this->d->lodCoveragePolicy.coverageComplete() &&
	!this->d->lodStructuralRepair.pointRelaxationPending() &&
	this->d->lodDiscoveryPointProxyPixelThreshold <= 1.01f &&
	controller_lod_adaptive_point_aggregation_allowed(this,
	    this->d->lodStaticQualityTrial.capacityRejected()) &&
	ongoingStableWork &&
	measuredCadPresentation &&
	ongoingStableTargetNanoseconds > 0.0L &&
	static_cast<long double>(elapsed) >
	    ongoingStableTargetNanoseconds * 1.50L) {
	const BObolLodAdmissionPlan pressurePlan =
	    BObolLodAdmissionPlanner::planPointInterrupted(
		this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor,
		BObolLodPointCalibrationGoal::RESPONSIVE,
		this->d->lodPresentationPointProxyPixelThreshold,
		elapsed, this->d->pointCalibrationTargetFps(
		    BObolLodPointCalibrationGoal::RESPONSIVE));
	this->d->commitAdmissionPlan(pressurePlan);
	const BObolLodPointProxyEvidence::Decision &pressure =
	    pressurePlan.pointProxyDecision;
	if (pressure.changed) {
	    this->d->publishCadPointProxyThreshold(pressure.threshold);
	    /* The changed SoCADAssembly occurrence population needs one exact
	     * successor frame before the bracket can be adjusted again.  This
	     * presentation-only calibration and retained triangle recovery are
	     * mutually exclusive phases: calibration pauses the retained producer,
	     * while recovery requires that producer to complete.  Arming both
	     * creates a no-owner liveness state. */
	    const bool structuralClassifierFrame =
		unresolvedStructuralOccurrences > 0;
	    bool classifierFrameRequested = false;
	    if (this->d->lodPresentationPointProxyPixelThreshold > 1.01f) {
		if (structuralClassifierFrame) {
		    this->requestStructuralPointAdmissionFrame(
			"lod-stream-point-pressure");
		    classifierFrameRequested = true;
		} else {
		    this->d->lodPointQualityPhase.requestCalibration();
		}
	    }
	    /* This pressure step is deliberately conservative because its frame
	     * includes streaming/publication latency.  The unchanged successor
	     * frame will walk back toward the finest sustainable presentation when
	     * it demonstrates headroom. */
	    if (!classifierFrameRequested)
		this->requestLodCapacityRender("lod-stream-point-pressure");
	}
    }

    /* A pressure step may have changed the aggregate threshold while source
     * or result streaming still owned the next frame.  If that producer goes
     * quiet with structural boxes in the exact completed population, the
     * calibration latch no longer has a frame it can legitimately measure.
     * Hand the frontier to armStableLodHeadroomProbeIfReady() below instead
     * of leaving both admission and repair permanently paused. */
    const bool pointCalibrationFutureProducer =
	ongoingStableWork ||
	this->d->lodPresentationTransaction.publicationPending();
    if (BObolLodAdmissionPlanner::pointCalibrationYieldsToStructuralRepair(
		this->d->lodPointQualityPhase.presentationPending(),
		exactOccurrenceClassification != FALSE,
		unresolvedStructuralOccurrences,
		pointCalibrationFutureProducer)) {
	this->d->lodPointQualityPhase.completeCalibration();
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
    if (this->d->lodInteractionSession.active() &&
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
		this->d->publishCadProgressiveCeiling(correctedCeiling);
		this->requestLodCapacityRender("lod-scale-quality-pressure");
	    }
	}
    }

    /* A budget-limited pass may have admitted a bounded first wave and then
     * scanned the remaining boxes without admitting them.  Re-plan from the
     * highest screen-value occurrence after that wave has actually rendered
     * and supplied a throughput sample.  Replanning before the frame would
     * repeatedly spend the same stale allowance. */
    const bool capacitySampleYieldsToStructuralFrontier =
	BObolLodAdmissionPlanner::capacitySampleYieldsToStructuralFrontier(
	    this->d->lodAdmissionEvidence.capacity().capacitySearch().
		awaitingPresentationFrame() &&
		!this->d->lodAdmissionEvidence.capacity().
		    calibrationFramePending(),
	    calibrationState &&
		calibrationState->lastCadPresentationFrameExact(),
	    exactOccurrenceClassification != FALSE,
	    presentedStructuralBoxes);
    if (capacitySampleYieldsToStructuralFrontier) {
	/* Structural placeholders are excluded from capacity timing.  This exact
	 * frame proves that an unchanged replay cannot satisfy the certificate;
	 * its successor is point classification or structural repair, not another
	 * capacity repaint.  Preserve the current allocation/budget and retire
	 * only the impossible measurement transaction. */
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::
		INVALIDATE_CAPACITY_MEASUREMENT);
	if (this->d->lodPointQualityPhase.presentationPending())
	    this->requestLodCapacityRender("lod-point-after-structural-capacity");
	else
	    this->markProgressiveWorkPending();
    } else if (this->d->lodAdmissionEvidence.capacity().
	    capacityTransactionPending()) {
	/*
	 * Calibration probes measure an unchanged retained cut.  They do not
	 * need an intervening O(N) occurrence-planning pass: that pass cannot
	 * admit anything until the samples have changed the aggregate budget.
	 * Present the bounded probe series back-to-back, then scan the sparse
	 * unsatisfied frontier once using the resulting calibration.
	 */
	const BObolLodCapacitySearchCertificate &capacitySearch =
	    this->d->lodAdmissionEvidence.capacity().capacitySearch();
	BObolLodCapacityEvidence::CompletedFrameDecision calibration;
	/* A population barrier predating an ALLOCATING capacity candidate owns
	 * this completed frame.  Consume it before inspecting the search phase;
	 * the transition preserves the exact candidate as the sole successor.
	 * Letting the broad awaitingSample() predicate mask this barrier pauses
	 * the candidate cursor behind the very frame which is already complete. */
	if (this->d->lodAdmissionEvidence.capacity().
		calibrationFramePending()) {
	    calibration = this->d->completeCapacityCalibrationFrame();
	} else if (capacitySearch.awaitingPresentationFrame()) {
	    const BObolLodCapacitySearchCertificate::Phase searchPhaseBefore =
		capacitySearch.phase();
	    const size_t searchCandidateBefore =
		capacitySearch.candidateBudget();
	    const unsigned int searchSamplesBefore =
		capacitySearch.samplesRemaining();
	    const unsigned int searchMeasuredBefore =
		capacitySearch.measuredCandidateCount();
	    BObolLodCapacityEvidence::CompletedFrameInputs inputs;
	    inputs.searchKey = capacitySearch.key();
	    inputs.candidateBudget =
		capacitySearch.candidateBudget();
	    const BObolRetainedAllocationResult &capacityAllocation =
		this->d->lodRetainedAllocationCertificate;
	    inputs.populationDigest =
		capacityAllocation.selectedPopulationDigest;
	    inputs.populationIdentity =
		capacityAllocation.selectedPopulationIdentity;
	    inputs.populationMinimumBudget =
		capacityAllocation.selectedPresentationCost;
	    inputs.nextDistinctPopulationBudget =
		capacityAllocation.nextDistinctPresentationBudget;
	    inputs.nextDistinctPopulationBudgetKnown =
		capacityAllocation.nextDistinctPresentationBudgetKnown;
	    const bool currentCapacityAllocation =
		this->d->retainedAllocationPresentationRealized(
		    calibrationState);
	    const bool capacityAllocationCertificateCurrent =
		this->d->retainedAllocationCertificateCurrent(calibrationState);
	    const bool capacityAllocationCutsApplied =
		this->d->retainedAllocationCutsApplied(calibrationState);
	    const BObolLodAdmissionRevisionStamp currentRevision =
		this->d->admissionRevisionStamp();
	    const bool capacitySearchRevisionCurrent =
		capacitySearch.key().inventory == currentRevision.inventory() &&
		capacitySearch.key().availability ==
		    currentRevision.availability() &&
		capacitySearch.key().view == currentRevision.view() &&
		capacitySearch.key().policy == currentRevision.policy();
	    const bool exactCapacityPopulation =
		BObolLodAdmissionPlanner::capacitySamplePopulationReady(
		    haveCadPresentationAssemblies != FALSE,
		    calibrationState &&
			calibrationState->lastCadPresentationFrameExact(),
		    exactOccurrenceClassification != FALSE,
		    presentedStructuralBoxes, currentCapacityAllocation);
	    /* The completed-work record may count a one-time upload, an expanded
	     * OSMesa stream, or a later prepared replay.  Those backend paths can
	     * report different work for the same immutable occurrence population.
	     * The certificate is keyed by that population, so every sample must use
	     * its retained allocator cost; elapsed time already captures backend
	     * performance differences. */
	    const size_t capacityPresentedCost = currentCapacityAllocation ?
		BObolLodAdmissionPlanner::capacitySearchPresentedCost(
		    true, calibrationCost,
		    capacityAllocation.selectedPresentationCost) : 0;
	    const bool validCapacitySample =
		!this->d->lodInteractionSession.active() &&
		!this->d->lodInteractionSession.gestureActive() &&
		exactCapacityPopulation && reusableCadWorkSample &&
		!transientStablePublication && !transientStablePreparation &&
		capacityPresentedCost > 0 &&
		calibrationElapsed > 0;
	    inputs.presentedCost = validCapacitySample ?
		capacityPresentedCost : 0;
	    inputs.observedNanoseconds = validCapacitySample ?
		calibrationElapsed : 0;
	    inputs.validSample = validCapacitySample;
	    using CandidateState = BObolLodCapacityEvidence::
		CompletedFrameInputs::CandidateState;
	    inputs.candidateState =
		!capacityAllocationCertificateCurrent ||
		!capacityAllocationCutsApplied ?
		    CandidateState::REALLOCATION_REQUIRED :
		    (currentCapacityAllocation ? CandidateState::CURRENT :
			CandidateState::PRESENTATION_PENDING);
	    inputs.reallocationProducerPending =
		this->d->lodAvailabilityLedger.residentGrowthPending();
	    /* The active certificate already owns its typed safe bracket.  Do not
	     * inject a presentation accepted under a different deadline. */
	    inputs.knownSafeBudget = capacitySearch.safeBudget();
	    calibration = this->d->completeCapacitySearchFrame(inputs);
	    if (calibration.searchResult ==
		    BObolLodCapacitySearchCertificate::Result::CERTIFIED) {
		/* The exact ceiling-free samples which published this terminal
		 * certificate supersede any older motion/deadline handoff.  Keeping
		 * that handoff would immediately replace the certified allocation
		 * with its stale, lower reconciliation budget. */
		this->d->lodPresentationPolicy.acceptCapacityCertificate();
	    }
	    static std::atomic<unsigned int> capacityTraceCount(0);
	    const bool capacityTraceEnabled =
		controller_lod_trace_enabled("BOBOL_LOD_TRACE_CAPACITY",
		    this->d->lodViewRevision.value()) ||
		controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
		    this->d->lodViewRevision.value());
	    const unsigned int capacityTraceIndex = capacityTraceEnabled ?
		capacityTraceCount.fetch_add(1) : 0;
	    if (capacityTraceEnabled && capacityTraceIndex < 256)
		bu_log("BObol LoD completed capacity frame "
		       "phase_before=%u phase_after=%u result=%u "
		       "candidate_before=%zu candidate_after=%zu "
		       "sample=%d restart=%d valid=%d presented=%zu "
		       "signature=%llu elapsed_ms=%.3f "
		       "remaining_before=%u remaining_after=%u "
		       "measured_before=%u measured_after=%u "
		       "allocation_current=%d "
		       "active_cost=%zu selected_cost=%zu one_pixel=%d "
		       "exact_capacity=%d reusable=%d "
		       "publication=%d preparation=%d frame_exact=%d "
		       "occurrence_exact=%d subpixel=%zu boxes=%zu "
		       "search_revision_current=%d search_epochs=%llu/%llu/%llu/%llu "
		       "current_epochs=%llu/%llu/%llu/%llu "
		       "allocation_plan=%llu active_plan=%llu "
		       "allocation_view=%llu current_view=%llu "
		       "allocation_policy=%llu current_policy=%llu "
		       "cuts_applied=%d progressive_ceiling=%d "
		       "active_progressive_max=%d payloads=%zu "
		       "candidate_state=%u producer_pending=%d\n",
		       static_cast<unsigned int>(searchPhaseBefore),
		       static_cast<unsigned int>(capacitySearch.phase()),
		       static_cast<unsigned int>(calibration.searchResult),
		       searchCandidateBefore,
		       capacitySearch.candidateBudget(),
		       calibration.requestSampleFrame ? 1 : 0,
		       calibration.restartSubmission ? 1 : 0,
		       validCapacitySample ? 1 : 0,
		       inputs.presentedCost,
		       static_cast<unsigned long long>(
			   inputs.populationDigest),
		       inputs.observedNanoseconds / 1000000.0,
		       searchSamplesBefore,
		       capacitySearch.samplesRemaining(),
		       searchMeasuredBefore,
		       capacitySearch.measuredCandidateCount(),
		       currentCapacityAllocation ? 1 : 0,
		       activeCalibrationCost,
		       capacityAllocation.selectedPresentationCost,
		       exactRetainedPopulation ? 1 : 0,
		       exactCapacityPopulation ? 1 : 0,
		       reusableCadWorkSample ? 1 : 0,
		       transientStablePublication ? 1 : 0,
		       transientStablePreparation ? 1 : 0,
		       calibrationState &&
			   calibrationState->lastCadPresentationFrameExact() ? 1 : 0,
		       exactOccurrenceClassification ? 1 : 0,
		       presentedSubpixelOccurrences,
		       presentedStructuralBoxes,
		       capacitySearchRevisionCurrent ? 1 : 0,
		       static_cast<unsigned long long>(
			   capacitySearch.key().inventory.value()),
		       static_cast<unsigned long long>(
			   capacitySearch.key().availability.value()),
		       static_cast<unsigned long long>(
			   capacitySearch.key().view.value()),
		       static_cast<unsigned long long>(
			   capacitySearch.key().policy.value()),
		       static_cast<unsigned long long>(currentRevision.inventory().value()),
		       static_cast<unsigned long long>(currentRevision.availability().value()),
		       static_cast<unsigned long long>(currentRevision.view().value()),
		       static_cast<unsigned long long>(currentRevision.policy().value()),
		       static_cast<unsigned long long>(
			   capacityAllocation.allocationPlanSerial),
		       static_cast<unsigned long long>(calibrationState ?
			   calibrationState->activeCadAllocationPlan() : 0),
		       static_cast<unsigned long long>(capacityAllocation.viewRevision()),
		       static_cast<unsigned long long>(this->d->lodViewRevision.value()),
		       static_cast<unsigned long long>(capacityAllocation.policyRevision()),
		       static_cast<unsigned long long>(this->d->lodPolicyRevision.value()),
		       capacityAllocationCutsApplied ? 1 : 0,
		       this->d->lodInteractiveProgressiveCeiling,
		       calibrationState ?
			   calibrationState->maximumActiveProgressiveCut() : -1,
		       calibrationState ?
			   calibrationState->cadProgressivePayloadCount() : 0,
		       static_cast<unsigned int>(inputs.candidateState),
		       inputs.reallocationProducerPending ? 1 : 0);
	}
	if (calibration.requestCeilingFreeFrame) {
	    /* The current candidate already fits its certified budget.  Its
	     * retained cuts must stay fixed while the renderer-only ceiling is
	     * removed for the first exact measurement frame. */
	    this->d->publishCadProgressiveCeiling(-1);
	    this->requestLodCapacityRender("lod-capacity-ceiling-free-sample");
	} else if (calibration.requestSampleFrame) {
	    this->requestLodCapacityRender("lod-budget-sample");
	} else if (calibration.restartSubmission) {
	    this->restartLodCapacitySubmission(
		calibration.capacityCandidateChanged ? TRUE : FALSE);
	} else if (calibration.capacityCandidateChanged) {
	    /* Resident growth may already own the replacement allocation.  It does
	     * not, however, publish the capacity-search invalidation which selected
	     * that replacement.  Record the semantic edge now so the producer's
	     * eventual allocation cannot share a revision with the stale plan. */
	    this->d->advanceAdmissionRevision(
		BObolLodAdmissionRevisionDomain::CAPACITY);
	}
    }
    /* Once the scale-quality budget is active, every completed frame is a
	 * measurement of the current render ceiling even if a newer wheel event
	 * has already cleared the one-shot probe label.  Preserve that proof
	 * across the next camera epoch; otherwise ordinary motion
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
	(!this->d->staticQualityPresentationDeadline() ||
	 elapsed <= this->d->staticQualityPresentationDeadline())) {
	/* The completed-work record belongs to the framebuffer just presented.
	 * A following submission may already have updated the mutable scene
	 * status, so it is not a reason to discard this direct capacity proof.
	 * Retain the richest completed cut for this unchanged camera/population;
	 * a later coarser recovery frame is not evidence that the richer frame
	 * stopped fitting. */
	this->d->lodDeadlineSafePresentation.remember(
	    completedPresentedMaximum, this->d->lodViewRevision.value(),
	    this->d->lodViewQualityDomainRevision);
    }
    this->d->lodViewDemandPolicy.noteFramePresented(
	completedPresentedMaximum,
	measuredCadRenderCost != FALSE && presentedCadRenderCost > 0,
	elapsed);
    if (this->d->lodViewDemandPolicy.rearmAfterQualityFrame(
	    this->d->lodInteractionSession.active(),
	    calibrationState ?
		calibrationState->maximumActiveProgressiveCut() : -1,
	    completedPresentedMaximum,
	    measuredCadRenderCost != FALSE && presentedCadRenderCost > 0,
	    elapsed))
	this->markProgressiveWorkPending();
    size_t provenHandoffRenderCost = 0;
    if (measuredCadRenderCost && presentedCadRenderCost > 0 &&
	calibrationState && calibrationState->lastCadPresentationFrameExact() &&
	(!this->d->staticQualityPresentationDeadline() ||
	 elapsed <= this->d->staticQualityPresentationDeadline())) {
	/* OSMesa's exact submitted-work record expands indexed positions and
	 * normals, while retained allocation budgets canonical indexed mesh work.
	 * A handoff is an allocator certificate, so translate the exact presented
	 * global prefix into that currency before extrapolating its deadline-safe
	 * allowance.  Feeding backend work directly into the allocator inflated a
	 * cut-zero multi-mesh proof several-fold and left rich hidden cuts behind a
	 * permanent renderer ceiling. */
	if (completedPresentedMaximum >= 0) {
	    const size_t presentedCadCost = calibrationState->
		cadRenderCostAtProgressiveCutCeiling(
		    completedPresentedMaximum);
	    provenHandoffRenderCost = BObolLodAdmissionPlanner::
		canonicalSceneCostAtCadCeiling(
		    calibrationState->activeRenderCost(),
		    calibrationState->activeCadRenderCost(),
		    presentedCadCost);
	} else if (this->d->retainedAllocationPresentationRealized(
		calibrationState)) {
	    provenHandoffRenderCost = this->d->
		lodRetainedAllocationCertificate.selectedPresentationCost;
	} else {
	    provenHandoffRenderCost = calibrationCost;
	}
    }
    this->completePresentationBarrier(elapsed,
	context.cadPresentationCompleteness ==
	    BObolCadPresentationCompleteness::EXACT ? TRUE : FALSE,
	provenHandoffRenderCost);

    /* A quiet single-occurrence handoff may have loaded several richer PoP
     * populations behind its last responsive renderer ceiling.  Expose them
     * only through completed-frame evidence.  The first hard-deadline miss
     * disables this sequence in notePresentationRenderInterrupted(), leaving
     * the last complete framebuffer and immutable resident suffix intact. */
    bool scheduledStaticPresentationStep = false;
    const bool transientStaticPresentation =
	transientStablePublication || transientStablePreparation;
    if (this->d->lodStaticQualityTrial.probing() &&
	!this->d->lodInteractionSession.active() &&
	this->d->lodInteractiveProgressiveCeiling >= 0 &&
	transientStaticPresentation &&
	!this->d->lodSubmissionPass.active() &&
	!this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() &&
	!this->d->lodPresentationTransaction.barrierPending() &&
	!this->d->lodPresentationTransaction.publicationPending()) {
	/* Completing an immutable atlas or VBO population can yield an exact
	 * framebuffer, but its elapsed time is not the sustainable cost of
	 * replaying that cut.  Keep the last completed image and request exactly
	 * one unchanged traversal; the cumulative upload/preparation counters make
	 * the successor a reusable sample.  Advancing or rejecting from this
	 * transient frame made large discrete meshes stop several cuts too coarse. */
	this->requestLodCapacityRender("lod-static-overscan-preparation-replay");
	scheduledStaticPresentationStep = true;
    } else if (this->d->lodStaticQualityTrial.probing() &&
	!this->d->lodInteractionSession.active() &&
	this->d->lodInteractiveProgressiveCeiling >= 0 &&
	calibrationState && measuredCadRenderCost &&
	presentedCadRenderCost > 0 &&
	calibrationState->lastCadPresentationFrameExact() &&
	elapsed <= this->d->staticQualityPresentationDeadline() &&
	!this->d->lodPresentationPolicy.handoffPending() &&
	!this->d->lodSubmissionPass.active() &&
	!this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() &&
	!this->d->lodPresentationTransaction.barrierPending() &&
	!this->d->lodPresentationTransaction.publicationPending()) {
	const int activeMaximum =
	    calibrationState->maximumActiveProgressiveCut();
	const int currentCeiling =
	    this->d->lodInteractiveProgressiveCeiling;
	const float currentNextFraction = calibrationState->
	    cadPresentationProgressiveCutNextFraction();
	if (currentNextFraction > 0.0f) {
	    /* This exact replay is the capacity proof for the deterministic
	     * page subset selected between two global PoP cuts.  Retain that
	     * framebuffer and render policy until a real view/policy/resource
	     * edge invalidates it; an occurrence-wide reconciliation cannot
	     * express this deliberately mixed page population. */
	    BObolLodStaticQualityTrial::Acceptance acceptance;
	    acceptance.revisionStamp = this->d->admissionRevisionStamp();
	    acceptance.ceiling = currentCeiling;
	    acceptance.nextFraction = currentNextFraction;
	    acceptance.presentedCost = calibrationCost;
	    acceptance.allowedCost = std::max(
		acceptance.presentedCost,
		BObolLodQualityPolicy::staticPresentationRenderCostLimit(
		    calibrationCost, elapsed,
		    this->d->staticQualityPresentationDeadline()));
	    (void)this->d->lodStaticQualityTrial.acceptFractionalCeiling(
		acceptance);
	    scheduledStaticPresentationStep = true;
	    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
		    this->d->lodViewRevision.value()))
		bu_log("BObol LoD static fractional ceiling accepted "
		       "ceiling=%d next_fraction=%.4f elapsed_ms=%.3f\n",
		       currentCeiling, currentNextFraction,
		       elapsed / 1000000.0);
	} else if (this->d->lodStaticQualityTrial.sampledCeiling() !=
		currentCeiling) {
	    this->d->lodStaticQualityTrial.noteSampledCeiling(currentCeiling);
	    this->requestLodCapacityRender("lod-static-overscan-steady-replay");
	    scheduledStaticPresentationStep = true;
	} else if (activeMaximum > currentCeiling) {
	    /* The renderer's completed-work counter describes backend-expanded
	     * flat-shading work on OSMesa, while the occurrence allocator uses
	     * canonical indexed-mesh cost units.  Query the exact guarded cut in
	     * allocator currency and pair that value with the measured elapsed
	     * time; mixing the two domains made a sustainable cut appear several
	     * times too expensive. */
	    const size_t presentedStaticCost =
		calibrationState->cadRenderCostAtProgressiveCutCeiling(
		    currentCeiling);
	    const size_t predictedCostLimit =
		BObolLodQualityPolicy::staticPresentationRenderCostLimit(
		    presentedStaticCost, elapsed,
		    this->d->staticQualityPresentationDeadline());
	    const int boundedCeiling = calibrationState->
		singleCadProgressiveCutWithinRenderCost(predictedCostLimit);
	    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
		    this->d->lodViewRevision.value())) {
		const size_t nextCeilingCost = currentCeiling < activeMaximum ?
		    calibrationState->cadRenderCostAtProgressiveCutCeiling(
			currentCeiling + 1) : presentedStaticCost;
		bu_log("BObol LoD static overscan frame ceiling=%d "
		       "active_max=%d bounded=%d presented_cost=%zu "
		       "predicted_cost=%zu current_cut_cost=%zu "
		       "next_cut_cost=%zu elapsed_ms=%.3f transient=%d\n",
		       this->d->lodInteractiveProgressiveCeiling,
		       activeMaximum, boundedCeiling, presentedStaticCost,
		       predictedCostLimit, presentedStaticCost,
		       nextCeilingCost, elapsed / 1000000.0,
		       transientStaticPresentation ? 1 : 0);
	    }
	    if (boundedCeiling >= 0 && boundedCeiling <=
		    this->d->lodInteractiveProgressiveCeiling) {
		const float nextFraction = calibrationState->
		    singleCadProgressiveNextFractionWithinRenderCost(
			predictedCostLimit, currentCeiling);
		if (nextFraction > 0.0f) {
		    /* Spatial pages are independent retained parts.  A fractional
		     * ceiling spends the measured residual capacity by promoting a
		     * deterministic subset one cut, avoiding an all-pages jump which
		     * would exceed the hard static deadline. */
		    this->d->publishCadProgressiveCeiling(
			currentCeiling, nextFraction);
		    this->requestLodCapacityRender("lod-static-overscan-fraction");
		    scheduledStaticPresentationStep = true;
		} else {
		/* The exact cached next-cut cost does not fit the throughput
		 * demonstrated by this completed frame with a small jitter margin.
		 * The current framebuffer is therefore the terminal static-quality
		 * proof.  Reconcile hidden richer occurrence metadata to that cut
		 * under the same 10 Hz allowance before retiring the renderer-only
		 * ceiling; no speculative long frame or coarse restart is needed. */
		BObolLodStaticQualityTrial::Constraint constraint;
		constraint.reason = BObolLodStaticQualityTrial::ConstraintReason::
		    PREDICTED_NEXT_CUT;
		constraint.revisionStamp = this->d->admissionRevisionStamp();
		constraint.committedCeiling = currentCeiling;
		constraint.committedCost = presentedStaticCost;
		constraint.candidateCeiling = currentCeiling + 1;
		constraint.candidateCost = calibrationState->
		    cadRenderCostAtProgressiveCutCeiling(currentCeiling + 1);
		constraint.allowedCost = predictedCostLimit;
		(void)this->d->lodStaticQualityTrial.reject(constraint);
		this->d->lodPresentationPolicy.armHandoff(
		    false, presentedStaticCost, predictedCostLimit);
		this->d->rewindLodSubmissionCursor();
		this->d->lodSubmissionPass.beginFresh();
		this->d->advanceAdmissionRevision(
		    BObolLodAdmissionRevisionDomain::CAPACITY);
		this->markProgressiveWorkPending();
		scheduledStaticPresentationStep = true;
		}
	    } else {
		/* The canonical per-cut table describes the exact layered population
		 * in allocator currency.  Advance directly to its richest predicted-
		 * safe cut; every renderer submission remains deadline-interruptible
		 * and preserves the preceding completed framebuffer if the timing
		 * model is optimistic.  Walking intermediate cuts needlessly rebuilds
		 * large software flat-shading ranges. */
		const int nextCeiling = boundedCeiling;
		if (this->d->lodStaticQualityTrial.probing()) {
		    this->d->publishCadProgressiveCeiling(nextCeiling);
		    this->d->lodStaticQualityTrial.resetSample();
		    this->requestLodCapacityRender("lod-static-overscan-step");
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
	    this->d->publishCadProgressiveCeiling(-1);
	    /* Retiring a redundant renderer ceiling does not retire the accepted
	     * static-quality phase.  Its occurrence allocation was measured under
	     * the event-driven hard deadline and must remain the terminal policy
	     * until a camera, resource, or capacity edge invalidates it.  Clearing
	     * the phase here returned the next pass to ordinary quiet streaming
	     * budget; that coarsened the scene, re-armed the 10 Hz phase, and made
	     * large warm scenes alternate forever between the two allocations. */
	    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
		    this->d->lodViewRevision.value()))
		bu_log("BObol LoD static overscan terminal active=%d "
		       "rejected=%d active_max=%d\n",
		       this->d->lodStaticQualityTrial.inProgress() ? 1 : 0,
		       this->d->lodStaticQualityTrial.capacityRejected() ? 1 : 0,
		       activeMaximum);
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
	SbBool exactFrame, size_t provenRenderCost)
{
    BObolLodControlTransitionScope controlTransition(this);
    const BObolViewLodState *initialPresentationState =
	this->d->viewAttachment ?
	    this->d->viewAttachment->getViewLodState() : NULL;
    const bool allocationCurrentBeforeBarrier =
	this->d->retainedAllocationPresentationRealized(
	    initialPresentationState);
    const uint64_t capacityRevisionBeforeBarrier =
	this->d->admissionRevisionStamp().capacity().value();
    bobol_identity_advance(this->d->renderCompletionSerial);

    const size_t reconciliationBudget =
	BObolLodQualityPolicy::staticPresentationRenderCostLimit(
	    provenRenderCost, elapsedNanoseconds,
	    this->d->staticQualityPresentationDeadline());
    const bool deadlineHandoffReady =
	this->d->lodPresentationPolicy.noteFramePresented(
	    provenRenderCost, reconciliationBudget);
    if (deadlineHandoffReady) {
	/* The constrained retry frame is the proof which a deadline-created
	 * handoff was waiting for.  The admission pass normally ran before that
	 * frame and correctly left the ceiling armed; schedule the successor pass
	 * now so it can retire the handoff and restore the quiet presentation. */
	this->d->rewindLodSubmissionCursor();
	this->d->lodSubmissionPass.beginFresh();
	this->d->resetRetainedPassAnnotations();
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::CAPACITY);
	this->markProgressiveWorkPending();
    }

    /* Rendering time is not user-idle time.  Start the quiet clock only
     * after a frame for the newest camera has actually completed. */
    if (this->d->lodInteractionSession.noteMotionFrameCompleted(
	    this->d->renderCompletionSerial, bu_gettime())) {
	if (this->d->lodViewDemandPolicy.noteMotionFrameSettled())
	    this->markProgressiveWorkPending();
    }

    /* A partial PoP result must be presented before the next, potentially
     * much larger prefix is requested.  This is presentation sequencing, not
     * capacity evidence: an exact selection/HUD frame which traversed the
     * current population satisfies the barrier just as an LoD-named frame
     * does. */
    const BObolLodPresentationTransaction::Completion frameCompletion =
	this->d->lodPresentationTransaction.complete(
	    this->d->renderCompletionSerial,
	    this->d->lodViewRevision.value(),
	    this->d->lodPolicyRevision.value(), exactFrame != FALSE);
    if (frameCompletion.retired) {
	/* A capacity candidate owns its complete sample sequence.  Result
	 * publication supplies its bounded frame cadence, so an additional
	 * render-duration cooldown would serialize independent asset loads behind
	 * those frames. */
	const bool capacityOwnsSuccessor =
	    this->d->lodAdmissionEvidence.capacity().
		capacityTransactionPending();
	const SbBool resumeResidentRefinement =
	    ((frameCompletion.reasons &
	      BObolLodPresentationTransaction::REASON_RESIDENT_REFINEMENT) != 0 ||
	     this->d->lodRetainedPass.residencyPending()) ? TRUE : FALSE;
	/* Retiring this latch acknowledges presentation ordering; it is not new
	 * capacity evidence.  Preserve an allocation which still certifies the
	 * exact presented occurrence population.  Ordinary result publication or
	 * cut mutation invalidates allocation population coverage independently,
	 * so only an already-stale plan needs a new capacity epoch here.  Advancing
	 * unconditionally after the final capacity sample made its freshly
	 * certified plan stale and stranded budget-limited meshes in REFINING with
	 * no enabled successor. */
	const BObolViewLodState *presentationState =
	    this->d->viewAttachment ?
		this->d->viewAttachment->getViewLodState() : NULL;
	if (!this->d->retainedAllocationPresentationRealized(
		presentationState))
	    this->d->advanceAdmissionRevision(
		BObolLodAdmissionRevisionDomain::CAPACITY);
	const int64_t cooldownMicroseconds =
	    BObolLodPresentationTransaction::refinementCooldownMicroseconds(
		elapsedNanoseconds, capacityOwnsSuccessor);
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
	    this->d->lodRetainedPass.clearResidencyPending();
	/* A capacity candidate owns its complete sample sequence.  The same
	 * framebuffer can also retire the generic cut-publication barrier, but
	 * that ordering edge must not reactivate occurrence planning between
	 * samples.  Doing both made a many-part scene alternate one unchanged
	 * measurement with one complete O(N) allocation pass.  The capacity
	 * transition above either requested its next frame or explicitly restarted
	 * submission; only the latter is entitled to make this cursor active. */
	if (!this->d->lodSubmissionPass.active() &&
	    !capacityOwnsSuccessor) {
	    this->d->rewindLodSubmissionCursor();
	    this->d->lodSubmissionPass.beginFresh();
	}
	if (this->d->lodSubmissionPass.active())
	    this->markProgressiveWorkPending();
    }

    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
	    this->d->lodViewRevision.value())) {
	static std::atomic<unsigned int> barrierTraceCount(0);
	const unsigned int traceIndex = barrierTraceCount.fetch_add(1);
	if (traceIndex < 256)
	    bu_log("BObol LoD presentation barrier allocation_before=%d "
		   "deadline_handoff=%d retired=%d reasons=%u "
		   "capacity_before=%llu capacity_after=%llu "
		   "allocation_after=%d\n",
		   allocationCurrentBeforeBarrier ? 1 : 0,
		   deadlineHandoffReady ? 1 : 0,
		   frameCompletion.retired ? 1 : 0,
		   static_cast<unsigned int>(frameCompletion.reasons),
		   static_cast<unsigned long long>(
		       capacityRevisionBeforeBarrier),
		   static_cast<unsigned long long>(
		       this->d->admissionRevisionStamp().capacity().value()),
		   this->d->retainedAllocationPresentationRealized(
		       initialPresentationState) ? 1 : 0);
    }

    this->scheduleResidentGrowthReallocationIfReady();
    (void)this->schedulePendingLodHandoffAllocationIfReady();
    if (this->d->lodViewDemandPolicy.qualityProbePending()) {
	this->markProgressiveWorkPending();
    }
    /* The completed-frame reducer may open exact presentation debt without a
     * provider callback edge.  Other successor policies publish their own
     * work while deciding them; this level-triggered exact latch is shared. */
    if (this->hasPendingLodRefinementFrame())
	this->markProgressiveWorkPending();
}

void
BObolViewController::beginSceneWideCapacitySubmission(SbBool active)
{
    BObolLodControlTransitionScope controlTransition(this);
    /* Capacity evidence and retained-allocation certificates describe the
     * complete current occurrence population.  A selective source delta is a
     * latency optimization for inventory or structural repair and cannot
     * apply such a certificate.  Retire that scope together with every
     * predecessor-pass annotation before exposing the new cursor. */
    this->d->lodSubmissionDelta.reset();
    this->d->rewindLodSubmissionCursor();
    this->d->resetRetainedPassAnnotations();
    this->d->lodSubmissionPass.beginFresh(active != FALSE);
}

void
BObolViewController::restartLodCapacitySubmission(
	SbBool capacityCandidateChanged)
{
    BObolLodControlTransitionScope controlTransition(this);
    if (capacityCandidateChanged) {
	(void)this->d->lodAvailabilityLedger.
	    transferResidentGrowthToCapacityCandidate();
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::CAPACITY);
    }
    this->beginSceneWideCapacitySubmission();
    /* A missing bounded coverage census retires its counter pass while
     * waiting for capacity evidence.  Its successor is another complete
     * current-view census, not merely a quality pass. */
    if (this->d->lodCoveragePolicy.required() &&
	!this->d->lodCoveragePolicy.coverageComplete() &&
	!this->d->lodCoveragePolicy.active())
	this->d->lodCoveragePolicy.activate(false);
    if (this->d->lodCoveragePolicy.active())
	this->d->lodCoveragePolicy.clearPassCounters();
    this->markProgressiveWorkPending();
}

SbBool
BObolViewController::scheduleCapacityCandidateAllocationIfReady(void)
{
    BObolLodControlTransitionScope controlTransition(this);
    if (!this->automaticLodControlEnabled() || !this->d->lodService ||
	!this->d->lodAdmissionEvidence.capacity().capacityAllocationPending() ||
	this->d->lodSubmissionPass.active() ||
	this->d->lodSubmissionPass.rescanPending() ||
	this->d->lodPresentationTransaction.barrierPending() ||
	this->d->lodPresentationTransaction.publicationPending())
	return FALSE;

    const std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_source_roots(this);
    if (sources.empty() || !this->d->resumeCapacityCandidateAllocation())
	return FALSE;

    (void)this->d->lodAvailabilityLedger.
	transferResidentGrowthToCapacityCandidate();
    this->beginSceneWideCapacitySubmission();
    if (!this->d->lodSubmissionPass.active())
	return FALSE;
    this->markProgressiveWorkPending();
    return TRUE;
}

SbBool
BObolViewController::schedulePendingLodHandoffAllocationIfReady(void)
{
    BObolLodControlTransitionScope controlTransition(this);
    if (!this->automaticLodControlEnabled())
	return FALSE;

    /* A handoff can become runnable at either of two asynchronous edges: its
     * last result-publication frame can complete, or background compaction can
     * release the service stream after that frame.  Keep the readiness test
     * and transition in one place so neither edge can strand allocation debt
     * behind an unchanged framebuffer. */
    const bool handoffServiceQuiescent = !this->d->lodService ||
	(this->d->lodService->activeTaskCountForGeneration(
	     this->d->lodActiveGeneration) == 0 &&
	 this->d->lodService->queuedResultCountForGeneration(
	     this->d->lodActiveGeneration) == 0 &&
	 this->d->lodAvailabilityLedger.resultsPending() == 0);
    const bool handoffPlanningReady =
	this->d->lodPresentationPolicy.handoffPending() &&
	!this->d->lodPresentationPolicy.handoffPresentationPending() &&
	!this->d->lodSubmissionPass.active() &&
	!this->d->lodSubmissionPass.rescanPending() &&
	!this->d->lodPresentationTransaction.barrierPending() &&
	!this->d->lodPresentationTransaction.publicationPending() &&
	!this->d->lodAvailabilityLedger.residentGrowthPending() &&
	!this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() &&
	handoffServiceQuiescent;
    if (!handoffPlanningReady)
	return FALSE;
    const std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_source_roots(this);
    if (sources.empty())
	return FALSE;
    const size_t delayedReconciliationBudget =
	this->d->lodPresentationPolicy.handoffReconciliationBudget();
    if (delayedReconciliationBudget > 0)
	this->d->requestPresentationReconciliation(
	    delayedReconciliationBudget);
    else
	this->d->requestRetainedReallocation();
    this->beginSceneWideCapacitySubmission();
    this->d->advanceAdmissionRevision(
	BObolLodAdmissionRevisionDomain::CAPACITY);
    this->markProgressiveWorkPending();
    return TRUE;
}

void
BObolViewController::scheduleLodRefinementFrame(const char *reason)
{
    BObolLodControlTransitionScope controlTransition(this);
    if (!this->automaticLodControlEnabled())
	return;

    /* This gate is useful only if a host presentation is guaranteed to
     * follow it.  Merely latching requestLodCapacityRender() is insufficient: a result
     * can be drained while progressiveWorkPending is already set, after
     * which the same advance clears that edge-triggered latch.  With no
     * unrelated Qt event, the controller would then wait forever for a frame
     * it never asked the host to schedule.
     *
     * Latch the render before invoking the callback so synchronous hosts also
     * observe a renderable request.  Do not move an existing gate forward;
     * every cut selected before the pending presentation belongs to the same
     * next frame. */
    const BObolLodPresentationTransaction::Reason obligationReason =
	this->d->lodRetainedPass.residencyPending() ?
	    BObolLodPresentationTransaction::REASON_RESIDENT_REFINEMENT :
	    BObolLodPresentationTransaction::REASON_CUT_PRESENTATION;
    (void)this->d->lodPresentationTransaction.arm(obligationReason,
	this->d->renderCompletionSerial,
	this->d->lodViewRevision.value(),
	this->d->lodPolicyRevision.value());
    const char *frameReason = reason ? reason : "lod-refinement-frame";
    /* Preserve a more specific render reason already installed by result
     * application side effects such as residency eviction. */
    if (!this->isRenderRequested())
	this->requestLodCapacityRender(frameReason);
}

void
BObolViewController::scheduleResidentGrowthReallocationIfReady(void)
{
    BObolLodControlTransitionScope controlTransition(this);
    if (!this->automaticLodControlEnabled() ||
	!this->d->lodService ||
	!this->d->lodAvailabilityLedger.residentGrowthPending())
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
	!this->d->lodSubmissionPass.active() &&
	!this->d->lodPresentationTransaction.barrierPending() &&
	!this->d->lodPresentationTransaction.publicationPending() &&
	!this->d->lodAdmissionEvidence.capacity().capacityTransactionPending();
    /* Pose-only motion owns a reversible renderer ceiling and deliberately
     * preserves retained occurrence cuts.  A zoom may spend newly available
     * residency during its bounded quality probes; otherwise retain this
     * obligation until the quiet transition. */
    const bool allocationAllowed =
	!this->d->lodInteractionSession.active() ||
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
		   this->d->lodSubmissionPass.active() ? 1 : 0,
		   this->d->lodCoveragePolicy.effectiveComplete() ? 1 : 0,
		   this->d->lodPresentationTransaction.barrierPending() ? 1 : 0,
		   this->d->lodPresentationTransaction.publicationPending() ? 1 : 0,
		   this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() ? 1 : 0);
    }

    /* A newly resident suffix can be the fact which makes coverage possible.
     * Coverage therefore cannot guard the drain itself: doing so makes an
     * incomplete coverage pass yield to resident growth while resident growth
     * restarts that same pass forever.  Drain the immutable suffix first;
     * after it reaches REALLOCATION_READY, consume the growth edge and hand
     * incomplete coverage to one ordinary successor below. */
    const bool residentWorkReady = streamIdle && presentationReady;
    if (this->d->lodAvailabilityLedger.beginResidencyDrainIfReady(
	    this->d->lodAutoSubmit != FALSE, residentWorkReady,
	    allocationAllowed)) {
	/* Request every still-needed immutable suffix before recomputing the
	 * scene-wide cut distribution.  This pass uses the current unsatisfied
	 * frontier and preserves every visible cut; a bounded worker queue may
	 * require several such waves, but none performs an O(scene) allocation or
	 * invalidates the last complete framebuffer. */
	this->d->rewindLodSubmissionCursor();
	/* A residency drain is an ordinary suffix-request pass.  Carrying the
	 * preceding retained-admission mode through resetPass() immediately
	 * re-armed the O(scene) minimax allocator and turned every result batch
	 * into another BALANCING/REFINING cycle. */
	this->d->lodSubmissionIntent.setRetainedAdmission(false);
	this->d->resetRetainedPassAnnotations();
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::CAPACITY);
	this->d->lodSubmissionPass.beginFresh();
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
		this->d->lodViewRevision.value()))
	    bu_log("BObol LoD resident growth scheduled residency drain "
		   "view=%llu policy=%llu\n",
		   (unsigned long long)this->d->lodViewRevision.value(),
		   (unsigned long long)this->d->lodPolicyRevision.value());
	this->markProgressiveWorkPending();
	return;
    }
    if (!this->d->lodAvailabilityLedger.consumeResidentGrowthIfReady(
	    this->d->lodAutoSubmit != FALSE,
	    residentWorkReady,
	    allocationAllowed))
	return;

    this->d->rewindLodSubmissionCursor();
    this->d->resetRetainedPassAnnotations();
    this->d->advanceAdmissionRevision(
        BObolLodAdmissionRevisionDomain::CAPACITY);
    const BObolLodAvailabilityScheduler::ResidentGrowthSuccessor successor =
	BObolLodAvailabilityScheduler::residentGrowthSuccessor(
	    this->d->lodCoveragePolicy.effectiveComplete());
    if (successor == BObolLodAvailabilityScheduler::
	    ResidentGrowthSuccessor::RETRY_COVERAGE) {
	this->d->lodCoveragePolicy.activate(false);
	this->d->lodSubmissionPass.beginFresh();
	this->markProgressiveWorkPending();
	return;
    }
    this->d->requestRetainedReallocation();
    /* This is an explicit pass in the current epoch.  Preserve the submitted
     * epoch witness so the wrapper cannot turn it into a view-change rescan. */
    this->d->lodSubmissionPass.beginFresh();
    this->d->resetLodConvergenceFraction();
    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
	    this->d->lodViewRevision.value()))
	bu_log("BObol LoD resident growth scheduled scene-wide reallocation "
	       "view=%llu policy=%llu budget=%zu\n",
	       (unsigned long long)this->d->lodViewRevision.value(),
	       (unsigned long long)this->d->lodPolicyRevision.value(),
	       this->d->lodAdmissionEvidence.capacity().currentBudget());
    this->markProgressiveWorkPending();
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
    BObolLodControlTransitionScope controlTransition(this);
    const uint64_t now = this->beginRenderTiming();
    std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
    bobol_saturating_counter_advance(this->d->presentedFrameSerial);
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
	    if (!this->d->lodInteractionSession.gestureActive())
		generalDisplaySample = interval;
	} else if (interval > 250000000ULL &&
		   !this->d->lodInteractionSession.gestureActive()) {
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
    if (this->d->lodInteractionSession.gestureActive()) {
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
    return this->d->lodInteractionSession.settleAfterRenderSerial();
}

uint64_t
BObolViewController::getLodRefinementResumeAfterRenderSerial(void) const
{
    return this->d->lodPresentationTransaction.requiredRenderSerial();
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

    /* Claim the pending transaction before traversal.  A retained mutation
     * published by a callback during traversal then creates a distinct
     * successor transaction. */
    SbBool lodCapacityRelevant = FALSE;
    SbBool lodPlanningRelevant = FALSE;
    (void)this->consumeRenderRequest(NULL, &lodCapacityRelevant,
	&lodPlanningRelevant);
    if (!this->isLodPresentationCapacityRelevant())
	lodCapacityRelevant = FALSE;
    const uint64_t started = this->beginRenderTiming();
    const BObolViewLodState *presentationState = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    if (presentationState)
	presentationState->beginCadPresentationFrame();
    const SbBool rendered = renderer->render(this->getRenderRoot());
    const uint64_t completed = this->beginRenderTiming();
    if (presentationState)
	presentationState->refreshCadPresentationFrameStatus();
    const BObolCadPreparationProgress cadPreparation = presentationState ?
	presentationState->cadPresentationPreparationProgress() :
	BOBOL_CAD_PREPARATION_NONE;
    const SbBool cadFrameIncomplete = presentationState &&
	!presentationState->lastCadPresentationFrameExact();
    const BObolCadPresentationExecution cadExecution = presentationState &&
	presentationState->lastCadPresentationFrameExecuted() ?
	    BObolCadPresentationExecution::EXECUTED :
	    BObolCadPresentationExecution::NOT_EXECUTED;
    if (!rendered) {
	return BRLCAD_ERROR;
    }
    if (cadFrameIncomplete) {
	/* SoOffscreenRenderer can return a valid partial image while a CAD
	 * assembly's bounded command preparation still has work remaining.  That
	 * image is a useful progressive preview, but it is not the completed frame
	 * required to retire publication, calibration, or convergence barriers.
	 * Treat it exactly like an interrupted on-screen traversal while excluding
	 * it from capacity feedback: notePresentationRenderInterrupted() retains
	 * the current population and level-triggers its presentation-only replay.
	 *
	 * A genuinely culled assembly must report an exact zero-work result.  If it
	 * does not, repeatedly accepting the stale report would strand visible
	 * occurrences just as surely as repeatedly requesting it would expose a
	 * renderer bug; keeping the obligation live makes that defect observable. */
	this->notePresentationRenderInterrupted(
	    completed > started ? completed - started : 1,
	    BObolPresentationTimingContext(
		BObolLodCapacityRelevance::EXCLUDED,
		lodPlanningRelevant ?
		    BObolLodPlanningRelevance::RELEVANT :
		    BObolLodPlanningRelevance::EXCLUDED,
		cadExecution,
		cadPreparation,
		BObolCadPresentationCompleteness::INCOMPLETE));
    } else {
	this->d->lastBackgroundRenderTimeNanoseconds = 0;
	this->d->lastSceneRenderTimeNanoseconds =
	    completed >= started ? completed - started : 0;
	this->completeRenderTiming(started,
	    BObolPresentationTimingContext(
		lodCapacityRelevant ?
		    BObolLodCapacityRelevance::RELEVANT :
		    BObolLodCapacityRelevance::EXCLUDED,
		lodPlanningRelevant ?
		    BObolLodPlanningRelevance::RELEVANT :
		    BObolLodPlanningRelevance::EXCLUDED,
		cadExecution,
		cadPreparation,
		BObolCadPresentationCompleteness::EXACT));
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
    return this->d->renderRequest.pending() ? TRUE : FALSE;
}

SbString
BObolViewController::getRenderReason(void) const
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    return this->d->renderRequest.reason();
}
