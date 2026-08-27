/*     V I E W _ C O N T R O L L E R _ P R O G R E S S I V E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file view_controller_progressive.cpp
 *
 * Progressive-provider scheduling and host-work state for a view controller.
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BLodService.h"
#include "BObol/BViewController.h"
#include "BObol/BViewLod.h"
#include "lod_coordinator_private.h"
#include "view_controller_private.h"

#include "bu/datetime.h"

#include <algorithm>
#include <limits>
#include <mutex>
#include <string>
#include <vector>

#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoOrthographicCamera.h>

static void
controller_accumulate_progressive_status(BObolProgressiveStatus &dst,
	const BObolProgressiveStatus &src)
{
    dst.providerCount += src.providerCount;
    dst.providerAdvanced += src.providerAdvanced;
    dst.lodResultsProcessed += src.lodResultsProcessed;
    dst.lodResultsApplied += src.lodResultsApplied;
    dst.submitted += src.submitted;
    dst.alreadyCached += src.alreadyCached;
    dst.expanded += src.expanded;
    dst.existing += src.existing;
    dst.remaining += src.remaining;
    dst.proxyPublished += src.proxyPublished;
    dst.metadataApplied += src.metadataApplied;
    dst.pendingTasks += src.pendingTasks;
    dst.inFlight += src.inFlight;
    dst.queuedResults += src.queuedResults;
    dst.queuedCacheWrites += src.queuedCacheWrites;
    if (src.changed)
	dst.changed = 1;
    if (src.hasMore)
	dst.hasMore = 1;
}

static bool
controller_has_cad_presentation(const BObolViewAttachment *attachment)
{
    const BObolViewLodState *state = attachment ?
	attachment->getViewLodState() : NULL;
    return state && state->hasCadPresentationAssemblies();
}

uint64_t
BObolViewController::registerProgressiveProvider(
    BObolProgressiveAdvanceCallback callback,
    void *userData,
    BObolProgressiveUserDataFreeCallback userDataFree)
{
    if (!callback)
	return 0;

    uint64_t token = this->d->progressiveProviderNextToken++;
    if (token == 0)
	token = this->d->progressiveProviderNextToken++;

    BObolProgressiveProviderRecord record;
    record.token = token;
    record.callback = callback;
    record.userData = userData;
    record.userDataFree = userDataFree;
    this->d->progressiveProviders.push_back(record);
    this->d->lodAvailabilityLedger.setProviderPendingCount(
	this->d->progressiveProviders.size());
    this->markProgressiveWorkPending();
    this->requestRender("progressive-provider");
    /* Registration composes into the standing host work level.  A conforming
     * host keeps its bounded service loop scheduled until that level drains. */
    return token;
}

void
BObolViewController::unregisterProgressiveProvider(uint64_t token)
{
    if (!token)
	return;

    for (std::vector<BObolProgressiveProviderRecord>::iterator it =
	 this->d->progressiveProviders.begin();
	 it != this->d->progressiveProviders.end(); ++it) {
	if (it->token != token)
	    continue;
	if (it->userDataFree)
	    (*it->userDataFree)(it->userData);
	this->d->progressiveProviders.erase(it);
	break;
    }
    this->d->lodAvailabilityLedger.setProviderPendingCount(std::min(
	this->d->lodAvailabilityLedger.providerPendingCount(),
	this->d->progressiveProviders.size()));
    if (this->d->progressiveProviders.empty()) {
	this->d->lodAvailabilityLedger.setProviderPendingCount(0);
	this->resetDiscoveryPointProxyFloor(TRUE);
	this->clearProgressiveWorkPending();
    }
}

void
BObolViewController::clearProgressiveProviders(void)
{
    for (const BObolProgressiveProviderRecord &record :
	 this->d->progressiveProviders)
	if (record.userDataFree)
	    (*record.userDataFree)(record.userData);
    this->d->progressiveProviders.clear();
    this->d->lodAvailabilityLedger.setProviderPendingCount(0);
    this->resetDiscoveryPointProxyFloor(TRUE);
    this->clearProgressiveWorkPending();
}

void *
BObolViewController::findProgressiveProviderData(
    BObolProgressiveAdvanceCallback callback) const
{
    if (!callback)
	return NULL;
    for (const BObolProgressiveProviderRecord &record :
	 this->d->progressiveProviders)
	if (record.callback == callback)
	    return record.userData;
    return NULL;
}

uint64_t
BObolViewController::findProgressiveProviderToken(
    BObolProgressiveAdvanceCallback callback) const
{
    if (!callback)
	return 0;
    for (const BObolProgressiveProviderRecord &record :
	 this->d->progressiveProviders)
	if (record.callback == callback)
	    return record.token;
    return 0;
}

SbBool
BObolViewController::hasProgressiveProviders(void) const
{
    return this->d->progressiveProviders.empty() ? FALSE : TRUE;
}

void
BObolViewController::setDefaultProgressiveOptions(
    const BObolProgressiveOptions *options)
{
    if (options)
	this->d->defaultProgressiveOptions = *options;
    else
	this->d->defaultProgressiveOptions = BObolProgressiveOptions();
    this->markProgressiveWorkPending();
    this->requestRender("progressive-options");
}

const BObolProgressiveOptions &
BObolViewController::getDefaultProgressiveOptions(void) const
{
    return this->d->defaultProgressiveOptions;
}

int
BObolViewController::advanceProgressiveWork(
    const BObolProgressiveOptions *options,
    BObolProgressiveStatus *status)
{
    const uint64_t advanceStarted = this->beginRenderTiming();

    const bool controllerOwnedDefault = options == NULL;
    if (!options)
	options = &this->d->defaultProgressiveOptions;
    BObolProgressiveOptions pacedOptions;
    SbBool terminalModeEntered = FALSE;
    if (controllerOwnedDefault) {
	pacedOptions = *options;
	/* A render initiated by a deterministic caller performs an internal
	 * controller-owned pump before traversing the scene.  Preserve the
	 * explicitly entered terminal mode across that pump; only another explicit
	 * options call may leave it. */
	if (this->d->forceTerminalLodRefinement)
	    pacedOptions.forceTerminalLodRefinement = TRUE;
	pacedOptions.maxProviderMicroseconds =
	    BObolLodAvailabilityScheduler::effectiveMicroseconds(
		options->maxProviderMicroseconds, true,
		this->d->lodInteractionSession.active() != FALSE);
	options = &pacedOptions;
    } else {
	terminalModeEntered =
	    !this->d->forceTerminalLodRefinement &&
	    options->forceTerminalLodRefinement ? TRUE : FALSE;
	const SbBool terminalModeChanged =
	    this->d->forceTerminalLodRefinement !=
	    options->forceTerminalLodRefinement ? TRUE : FALSE;
	this->d->forceTerminalLodRefinement =
	    options->forceTerminalLodRefinement;
	if (terminalModeChanged) {
	    /* The effective physical error target changes without changing the
	     * user's stored LoD policy.  Give that transient mode its own policy
	     * epoch so an already-settled view cannot take the unchanged-signature
	     * fast path and skip the terminal (or restored ordinary) demand pass.
	     * Do not arm a submission cursor directly: the policy may currently be
	     * disabled, in which case no legal consumer exists.  An enabled policy
	     * observes the new epoch and initializes its own cursor below. */
	    this->advanceLodPolicyRevision(
		LodPolicyTransition::PRESERVE_SCALE_DEMAND);
	    this->markProgressiveWorkPending();
	}
    }

    BObolProgressiveStatus localStatus;
    if (terminalModeEntered) {
	/* A deterministic/offline pump asks for the same view-dependent,
	 * pixel-exact terminal cut as an interactive host, but without waiting
	 * through responsiveness debounce tiers.  It also cannot inherit a
	 * temporary interaction/deadline presentation ceiling.  Such a handoff
	 * normally ends by certifying a finite occurrence-local allocation, while
	 * terminal mode deliberately installs an unbounded budget and therefore
	 * never enters retained allocation.  Leaving both policies active creates
	 * a closed no-work loop: the handoff repeatedly requests a certificate
	 * which the terminal budget can never produce.
	 *
	 * This is an entry action, not a mode invariant.  Repeating it on every
	 * host pump manufactures capacity revisions and discards the evidence
	 * needed to finish a finite allocation pass.
	 *
	 * Retire only renderer-local presentation constraints.  Resident PoP
	 * prefixes, occurrence cuts, source coverage, and provider work remain in
	 * place, so terminal refinement still starts from the useful retained
	 * scene and requests only missing suffixes. */
	this->d->lodRefinementNotBeforeMicroseconds = 0;
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::CAPACITY);
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_MEASUREMENT);
	this->d->lodPresentationPolicy.cancelHandoff();
	this->d->lodInteractiveProgressiveCeiling = -1;
	this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
	this->d->lodPointQualityPhase.reset();
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
	BObolViewLodState *terminalViewState =
	    this->d->viewAttachment->getViewLodState();
	if (terminalViewState) {
	    terminalViewState->setCadPresentationProgressiveCutCeiling(-1);
	    terminalViewState->setCadPresentationPointProxyPixelThreshold(1.0f);
	}
    }
    this->d->lastLodResultProcessingTimeNanoseconds = 0;
    this->d->lastProgressiveProviderTimeNanoseconds = 0;
    this->d->lastLodSubmissionTimeNanoseconds = 0;

    /* A deadline-interrupted CAD traversal is a closed presentation
     * transaction.  The renderer may have retained a resumable command-plan
     * cursor, so applying a provider batch, worker result, compaction, or cut
     * update here would invalidate work between replay slices.  Keep worker
     * queues bounded by their existing backpressure and let the already
     * requested successor frame commit before publishing more scene state. */
    if (this->d->lodInterruptedPresentationReplay.pending()) {
	localStatus.hasMore = 1;
	this->markProgressiveWorkPending();
	if (status)
	    *status = localStatus;
	const uint64_t advanceCompleted = this->beginRenderTiming();
	this->d->lastProgressiveAdvanceTimeNanoseconds =
	    advanceCompleted > advanceStarted ?
		advanceCompleted - advanceStarted : 0;
	return 1;
    }

    /* Camera motion keeps a responsive aggregate cut.  A completed scale
     * frame may nevertheless spend one bounded quality step immediately;
     * this keeps continuous zoom progressive rather than magnifying one
     * coarse image until input stops.  Ordinary stable convergence still
     * begins after the 150 ms quiet debounce. */
    if (this->d->lodAutoSubmit && this->d->lodInteractionSession.active()) {
	const int64_t now = bu_gettime();
	/* An unbracketed camera change normally retires its motion-frame latch in
	 * completePresentationBarrier().  A deadline-aborted frame has no such
	 * completion edge, and the renderer may then have no reason to run again.
	 * The progressive timer is the interaction session's remaining scheduler
	 * edge, so let it consume the same bounded timeout before testing
	 * quietReady().  Leaving this transition exclusively in render completion
	 * stranded wheel zoom in DEBOUNCING with an empty renderer queue. */
	(void)this->d->lodInteractionSession.releaseExpiredMotionFrame(
	    now, BObolLodViewDemandPolicy::
		unbracketedQuietDebounceMicroseconds());
	if (this->d->lodViewDemandPolicy.qualityProbePending()) {
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    const int activeMaximum = viewState ?
		viewState->maximumActiveProgressiveCut() : -1;
	    BObolLodViewDemandPolicy::QualityProbeInputs probeInputs;
	    probeInputs.activeMaximum = activeMaximum;
	    probeInputs.presentationCeiling =
		this->d->lodInteractiveProgressiveCeiling;
	    size_t presentedWork = 0;
	    if (viewState && viewState->lastCadPresentationFrameExact() &&
		viewState->lastCadPresentedRenderCost(presentedWork))
		probeInputs.presentedWork = presentedWork;
	    probeInputs.refinementFramePending =
		this->d->lodPresentationTransaction.barrierPending();
	    probeInputs.publicationFramePending =
		this->d->lodPresentationTransaction.publicationFramePending();
	    probeInputs.motionFramePending =
		this->d->lodInteractionSession.settleAfterRenderSerial() != 0;
	    const BObolLodViewDemandPolicy::QualityProbeDecision probe =
		this->d->lodViewDemandPolicy.beginQualityProbe(probeInputs);
	    if (probe.begin) {
		/* Offer one useful population per completed-frame witness.  Resident
		 * data may already be much richer, but opening the ceiling directly to
		 * that maximum makes a missed software frame back down through every
		 * ordinal and permits later repaint edges to repeat the staircase. */
		this->d->lodInteractiveProgressiveCeiling =
		    probe.progressiveCeiling;
		viewState->setCadPresentationProgressiveCutCeiling(
		    this->d->lodInteractiveProgressiveCeiling);
		this->advanceLodPolicyRevision(
		    LodPolicyTransition::PRESERVE_SCALE_DEMAND);
		this->d->lodDiscretePopulationTrialPermit.grant();
		this->markProgressiveWorkPending();
		this->requestRender("lod-scale-interaction-refine");
	    }
	}
	if (this->d->lodInteractionSession.quietReady(now,
		BObolLodViewDemandPolicy::
		    unbracketedQuietDebounceMicroseconds())) {
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    const bool returnedToStartingScale =
		this->d->lodViewDemandPolicy.interactionScaleChanged() &&
		this->d->lodInteractionStartCertificate.returnedToScale(
		    this->d->lodViewScaleSignature);
	    const bool scaleChanged =
		this->d->lodViewDemandPolicy.interactionScaleChanged() &&
		!returnedToStartingScale;
	    const bool haveRetainedMeshPayloads =
		this->getActiveLodMeshPayloadCount() > 0;
	    BObolLodPresentationPolicy::RestoreInputs restoreInputs;
	    restoreInputs.orthographic = this->d->activeCamera &&
		this->d->activeCamera->isOfType(
		    SoOrthographicCamera::getClassTypeId());
	    restoreInputs.scaleChanged = scaleChanged;
	    restoreInputs.retainedMeshPayloads = haveRetainedMeshPayloads;
	    restoreInputs.viewEpoch = this->d->lodViewRevision;
	    restoreInputs.population =
		controller_lod_presentation_population(viewState,
		    this->d->lodViewQualityDomainRevision);
	    restoreInputs.currentTargetPixelError =
		this->d->lodTargetPixelError;
	    restoreInputs.currentProgressiveCeiling =
		this->d->lodInteractiveProgressiveCeiling;
	    restoreInputs.currentProgressiveNextFraction = viewState ?
		viewState->cadPresentationProgressiveCutNextFraction() : 0.0f;
	    restoreInputs.currentPointProxyPixelThreshold =
		this->d->lodPresentationPointProxyPixelThreshold;
	    BObolLodViewQualityHistory::Quality recalledViewQuality;
	    BObolLodViewQualityHistory::RecallInputs recallInputs;
	    recallInputs.view = this->d->lodViewSignature.value_or(
		BObolLodViewSnapshot());
	    recallInputs.domainRevision =
		this->d->lodViewQualityDomainRevision;
	    /* History is also useful after resident-memory compaction has reduced
	     * a still-valid scene to proxies.  Source-domain identity proves which
	     * scene the entry belongs to; live provider and memory admission decide
	     * how much of the remembered target can be restored. */
	    recallInputs.sceneAvailable = viewState &&
		this->d->lodViewSignature &&
		(haveRetainedMeshPayloads ||
		 this->getActiveLodCadPayloadCount() > 0 ||
		 !this->d->lodLastSubmittedSources.empty());
	    recallInputs.resourcePressure =
		this->d->lodAdmissionEvidence.resources().anyPressure();
	    const bool recalledExactView =
		this->d->lodViewQualityHistory.recall(
		    recallInputs, recalledViewQuality);
	    /* A pose-only quiet pass still rechecks current visibility and
	     * resource pressure, but a fully covered orthographic population has
	     * no depth-dependent LoD demand.  Preserve its occurrence cuts unless
	     * measured FPS/resource admission explicitly requires coarsening. */
	    this->d->lodPoseContinuity.setRetainOccurrenceCuts(
		restoreInputs.orthographic && !scaleChanged &&
		this->d->lodInteractionStartCertificate.
		    startedFromReadyView() &&
		haveRetainedMeshPayloads);
	    /* Orthographic depth is irrelevant, but rotation still changes
	     * projected area and silhouette.  Perspective pose and scale changes
	     * can move importance as well.  Finish one exact demand census before
	     * authorizing a single redistribution of the retained scene budget. */
	    this->d->lodPlanningObligations.setImportanceCensus(
		this->d->lodInteractionStartCertificate.
		    startedFromReadyView() &&
		haveRetainedMeshPayloads);
	    if (this->d->forceTerminalLodRefinement) {
		this->d->lodTargetPixelError =
		    std::min(this->d->lodTargetPixelError, 0.25f);
		restoreInputs.currentTargetPixelError =
		    this->d->lodTargetPixelError;
		this->d->lodPlanningObligations.retireImportanceCensus();
	    }
	    BObolLodPresentationPolicy::QuietInputs quietInputs;
	    quietInputs.presentation = restoreInputs;
	    quietInputs.unboundedTerminal =
		this->d->forceTerminalLodRefinement != FALSE;
	    if (recalledExactView && !quietInputs.unboundedTerminal) {
		quietInputs.exactView.valid = true;
		quietInputs.exactView.targetPixelError =
		    recalledViewQuality.targetPixelError;
		quietInputs.exactView.progressiveCeiling =
		    recalledViewQuality.progressiveCeiling;
		quietInputs.exactView.progressiveNextFraction =
		    recalledViewQuality.progressiveNextFraction;
		quietInputs.exactView.pointProxyPixelThreshold =
		    recalledViewQuality.pointProxyPixelThreshold;
		quietInputs.exactView.presentedRenderCost =
		    recalledViewQuality.presentedRenderCost;
		quietInputs.exactView.provenRenderCostCapacity =
		    recalledViewQuality.provenRenderCostCapacity;
	    }
	    const BObolLodPresentationPolicy::RestoreDecision restore =
		this->d->lodPresentationPolicy.beginQuiet(quietInputs);
	    if (restore.apply) {
		this->d->lodTargetPixelError = restore.targetPixelError;
		this->d->lodInteractiveProgressiveCeiling =
		    restore.progressiveCeiling;
		this->d->lodPresentationPointProxyPixelThreshold =
		    restore.pointProxyPixelThreshold;
		if (viewState) {
		    viewState->setCadPresentationProgressiveCutCeiling(
			this->d->lodInteractiveProgressiveCeiling,
			restore.progressiveNextFraction);
		    viewState->setCadPresentationPointProxyPixelThreshold(
			this->d->lodPresentationPointProxyPixelThreshold);
		}
	    }
	    this->d->lodInteractionSession.finishQuiet();
	    const SbBool completedScaleInteraction =
		this->d->lodViewDemandPolicy.finishQuietInteraction(
		    returnedToStartingScale) ?
		    TRUE : FALSE;
	    /* Motion calibration targets 60 FPS and may intentionally lower the
	     * aggregate allowance behind a renderer-only ceiling.  Stable drawing
	     * targets 20 FPS; seed its handoff from the last proven stable scene
	     * budget, not from the transient motion floor.  Genuine quiet-frame
	     * overload evidence may still reduce this restored value. */
	    if (!this->d->forceTerminalLodRefinement &&
		this->d->lodInteractionStartCertificate.captured())
		this->d->commitAdmissionPlan(
		    BObolLodAdmissionPlanner::raiseCapacityBudget(
			this->d->lodAdmissionEvidence,
			this->d->lodAdmissionCursor,
			this->d->lodInteractionStartCertificate.
			    stableBudget()));
	    if (restore.restoredExactView)
		this->d->commitAdmissionPlan(
		    BObolLodAdmissionPlanner::raiseCapacityBudget(
			this->d->lodAdmissionEvidence,
			this->d->lodAdmissionCursor,
			restore.provenRenderCostCapacity));
	    this->d->lodInteractionStartCertificate.reset();
	    this->d->lodDiscretePopulationTrialPermit.revoke();
	    /* A motion ceiling also applies to native progressive wire stored in
	     * the standing CAD assembly.  It has no view-state mesh payload and
	     * therefore no occurrence pass capable of completing the mesh
	     * handoff below.  Restore its stable range directly; the policy
	     * revision already requests the one frame needed to present it. */
	    this->d->applyAdmissionEvidenceAction(
	        BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_MEASUREMENT);
	    this->d->applyAdmissionEvidenceAction(
	        BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_OVERLOAD);
	    this->d->lodPointQualityPhase.reset();
	    if (viewState) {
		viewState->setCadPresentationCameraMotionFrameReuse(FALSE);
		/*
		 * activeFaceCount() is the retained occurrence population, not
		 * necessarily the cut just presented under the motion ceiling.
		 * Promoting the stable budget to that value authorized the hidden
		 * rich population before the first quiet capacity check.  Keep the
		 * measured presentation and reconcile the retained cuts first.
		 */
	    }
	    /* Moving from the responsive pixel target to the stable target is a
	     * continuation of zoom demand, not a new cold-coverage epoch. */
	    this->advanceLodPolicyRevision(completedScaleInteraction ?
		LodPolicyTransition::PRESERVE_SCALE_DEMAND :
		LodPolicyTransition::ORDINARY);
	} else {
	    /* Keep the host frame pump alive through the idle debounce. */
	    localStatus.hasMore = 1;
	}
    }

    const uint64_t residentConsumerId =
	this->d->residentMeshConsumerId();
    const size_t queuedCompactionResults = this->d->lodService ?
	this->d->lodService->
	    queuedResidentMeshCompactionResultCountForDiagnostics(
		residentConsumerId) : 0;
    if (queuedCompactionResults) {
	if (this->d->lodInteractionSession.active()) {
	    /* Keep the last complete rich generation drawable while the user is
	     * moving.  The bounded completion queue back-pressures compaction
	     * workers until the quiet-view publication phase resumes. */
	    localStatus.hasMore = 1;
	} else {
	    const size_t applyLimit = options->maxLodResults ?
		std::min<size_t>(options->maxLodResults, 1024) : 1024;
	    std::vector<BObolLodResidentCompaction> compacted;
	    this->d->lodService->drainResidentMeshCompactions(
		residentConsumerId, compacted, applyLimit);
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    size_t published = 0;
	    if (viewState) {
		for (const BObolLodResidentCompaction &result : compacted)
		    published +=
			viewState->applyResidentMeshCompaction(result);
	    }
	    localStatus.lodResultsProcessed += compacted.size();
	    localStatus.lodResultsApplied += published;
	    if (published)
		localStatus.changed = 1;
	    if (this->d->lodService->
		    queuedResidentMeshCompactionResultCountForDiagnostics(
			residentConsumerId) > 0)
		localStatus.hasMore = 1;
	}
    }

    SbBool nonLodPresentationChanged = localStatus.changed ? TRUE : FALSE;
    const size_t queuedLodResults = this->d->lodService ?
	this->d->lodService->queuedResultCountForGeneration(
	    this->d->lodActiveGeneration) : 0;
    const SbBool havePendingLodResults =
	this->hasPendingLodResults() || queuedLodResults > 0;
    SbBool coalesceLodResults = FALSE;
    const SbBool holdRicherResultsDuringInteraction =
	havePendingLodResults && this->d->lodInteractionSession.active() &&
	!this->d->lodViewDemandPolicy.interactionScaleChanged() &&
	this->d->lodCoveragePolicy.coverageComplete() &&
	this->getActiveLodMeshPayloadCount() > 0;
    if (havePendingLodResults && this->d->lodService) {
	/* The first mesh is the time-to-useful-visual boundary and is never
	 * delayed.  Once a scene already has mesh data, coalesce a continuing
	 * many-asset stream until a full apply wave is ready or its latency
	 * bound expires.  A final/isolated result (no outstanding worker work)
	 * publishes immediately, preserving single-huge-mesh and sparse-scene
	 * behavior.
	 *
	 * Stable bulk population uses a 250 ms ceiling: four visibly progressive
	 * updates per second without paying hundreds of increasingly expensive
	 * whole-scene update/render traversals.  Motion uses 100 ms so newly
	 * visible occurrences catch up promptly while the interactive face
	 * budget protects FPS. */
	const size_t applyWave = options->maxLodResults ?
	    options->maxLodResults : 256;
	const SbBool submissionPausedByPresentation =
	    BObolLodAdmissionPlanner::presentationPausesSubmission(
		this->d->lodPointAdmissionFrame.pending(),
		this->d->lodPointQualityPhase.presentationPending(),
		this->d->lodAdmissionEvidence.capacity().rescanAfterFrame(),
		controller_has_cad_presentation(this->d->viewAttachment)) ?
		TRUE : FALSE;
	const SbBool continuingProducerStream =
	    BObolLodProducerPolicy::canProduceGeometry(
		this->d->lodSubmissionPass.active() != FALSE,
		submissionPausedByPresentation != FALSE,
		this->d->lodAvailabilityLedger.providerPendingCount() > 0,
		this->d->lodService->activeTaskCountForGeneration(
		    this->d->lodActiveGeneration) > 0) ? TRUE : FALSE;
	const SbBool sceneHasMesh =
	    this->getActiveLodMeshPayloadCount() > 0;
	int64_t firstReady =
	    this->d->lodAvailabilityLedger.firstResultReadyMicroseconds();
	const int64_t now = bu_gettime();
	if (firstReady <= 0) {
	    firstReady =
		this->d->lodAvailabilityLedger.ensureFirstResultReady(now);
	}
	const int64_t latencyLimit =
	    this->d->lodInteractionSession.active() ? 100000 : 250000;
	const int64_t resultAge =
	    now >= firstReady ? now - firstReady : latencyLimit;
	coalesceLodResults =
	    sceneHasMesh && continuingProducerStream &&
	    queuedLodResults < applyWave && resultAge < latencyLimit;
    }
    if (holdRicherResultsDuringInteraction) {
	/*
	 * Worker-produced arrays are immutable and the result queue is bounded
	 * and coalescing (currently 2048 slots).  Keeping richer generations
	 * there preserves the renderer's last prepared generation during
	 * pose-only camera input.  Scale-changing input is different: zoom must
	 * publish newly loaded suffixes as bounded result waves so detail improves
	 * under the pointer rather than only after the quiet debounce.  Coverage
	 * is reset on every view revision, so a newly exposed box-only leaf is
	 * never held behind this optimization.
	 */
	localStatus.hasMore = 1;
    } else if (coalesceLodResults) {
	localStatus.hasMore = 1;
    } else if (havePendingLodResults) {
	const uint64_t resultStarted = this->beginRenderTiming();
	(void)this->processPendingLodResults(options->maxLodResults,
	    options->maxLodApplyMicroseconds);
	const uint64_t resultCompleted = this->beginRenderTiming();
	this->d->lastLodResultProcessingTimeNanoseconds =
	    resultCompleted > resultStarted ? resultCompleted - resultStarted : 0;
	localStatus.lodResultsProcessed += this->d->lastLodResultCount;
	localStatus.lodResultsApplied += this->d->lastLodAppliedResultCount;
	if (this->d->lastLodAppliedResultCount > 0)
	    localStatus.changed = 1;
	if (!this->d->lodService ||
	    this->d->lodService->queuedResultCountForGeneration(
		this->d->lodActiveGeneration) == 0) {
	    this->d->lodAvailabilityLedger.clearFirstResultReady();
	    /* A result can race the empty observation without producing another
	     * empty->non-empty callback.  Recheck after clearing and install a
	     * fresh age origin if necessary. */
	    if (this->d->lodService &&
		this->d->lodService->queuedResultCountForGeneration(
		    this->d->lodActiveGeneration) > 0) {
		this->d->lodAvailabilityLedger.noteResultsReady(bu_gettime());
	    }
	}
    }

    const size_t priorProviderPendingCount =
	this->d->lodAvailabilityLedger.providerPendingCount();
    SbBool providerPresentationChanged = FALSE;
    size_t providerLimit = options->maxProviders;
    size_t providerIndex = 0;
    size_t providerPendingCount = 0;
    const uint64_t providerStarted = this->beginRenderTiming();
    for (const BObolProgressiveProviderRecord &record :
	 this->d->progressiveProviders) {
	if (!record.callback)
	    continue;
	if (providerLimit && providerIndex >= providerLimit) {
	    providerPendingCount +=
		this->d->progressiveProviders.size() - providerIndex;
	    localStatus.hasMore = 1;
	    break;
	}

	BObolProgressiveStatus providerStatus;
	providerStatus.providerCount = 1;
	int providerRet = (*record.callback)(this, record.userData, options,
			 &providerStatus);
	if (providerRet > 0) {
	    providerStatus.providerAdvanced++;
	} else if (providerRet < 0) {
	    providerStatus.hasMore = 0;
	}
	controller_accumulate_progressive_status(localStatus,
		providerStatus);
	if (providerStatus.changed) {
	    providerPresentationChanged = TRUE;
	    /*
	     * A provider mutation is already installed in the retained scene, but
	     * it does not require one presentation per bounded merge slice.  Feed
	     * it through the same owner-thread publication gate as immutable LoD
	     * results.  The host's independent progressive timer keeps draining
	     * providers while this gate coalesces several cheap mutations into one
	     * frame.  A standing render request (notably provider registration)
	     * still presents the first useful batch immediately, and the idle edge
	     * below publishes the last batch without waiting for the deadline.
	     */
	    this->d->lodPresentationTransaction.noteApplied(1, bu_gettime(),
		this->d->lodViewRevision.value(),
		this->d->lodPolicyRevision.value());
	}
	if (providerStatus.hasMore)
	    providerPendingCount++;
	providerIndex++;
    }
    this->d->lodAvailabilityLedger.setProviderPendingCount(
	providerPendingCount);

    /* A producer may publish its final geometry in one pass and retire its
     * staging lease in a later no-op pass.  That pending->settled edge can be
     * the last missing premise for terminal structural repair and exact-view
     * history.  Re-evaluate it here when the current pass did not mutate the
     * scene; a mutating pass first needs the ordinary requested presentation
     * to establish an exact completed-frame witness.  Without this edge the
     * host can correctly observe no more background work and go idle while
     * the last proven presentation still contains boxes. */
    if (priorProviderPendingCount > 0 && providerPendingCount == 0 &&
	!providerPresentationChanged)
	this->armStableLodHeadroomProbeIfReady();
    const uint64_t providerCompleted = this->beginRenderTiming();
    this->d->lastProgressiveProviderTimeNanoseconds =
	providerCompleted > providerStarted ?
	providerCompleted - providerStarted : 0;

    /* Providers publish the latest compact occurrences and their source-mesh
     * requests during this pump.  Submit view demand afterward so a provider
     * that completes in one pass cannot leave a newly published model at its
     * boxes forever simply because there is no reason to schedule a second
     * pass. */
    struct bv_lod_policy lodPolicy;
    bv_lod_policy_init(&lodPolicy);
    this->d->viewAttachment->getLodPolicy(&lodPolicy);
    this->d->lodCoveragePolicy.setRequired(
	this->d->lodAutoSubmit && this->d->lodService &&
	lodPolicy.policy != BV_LOD_OFF && lodPolicy.mesh_enabled);
    if (this->d->lodService) {
	const uint64_t admissionRevision =
	    this->d->lodService->residentMeshAdmissionRevision();
	if (this->d->lodResidentAdmissionRevision != 0 &&
	    admissionRevision !=
		this->d->lodResidentAdmissionRevision) {
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    if (viewState &&
		viewState->memoryLimitedPayloadCount() > 0) {
		/* Do not invalidate a cursor already proving coverage or current
		 * view demand.  Leave the observed revision unchanged; the first
		 * pump after that pass completes consumes the newest (coalesced)
		 * capacity edge.  Restarting the active pass at zero once per
		 * compacted asset made a four-item recovery rescan 150k entries
		 * dozens of times. */
		if (!this->d->lodSubmissionPass.active()) {
		    this->d->advanceAdmissionRevision(
			BObolLodAdmissionRevisionDomain::AVAILABILITY);
		    this->d->lodSubmissionSourceIndex = 0;
		    this->d->lodSubmissionEntryOffset = 0;
		    this->d->clearLodSubmissionPlan();
		    this->d->lodSubmissionPass.clearRescan();
		    this->d->lodPlanningObligations.
			setResidentAdmissionRetry(
			    this->d->lodCoveragePolicy.effectiveComplete());
		    this->d->lodSubmissionPass.activate();
		    this->d->lodResidentAdmissionRevision =
			admissionRevision;
		    this->markProgressiveWorkPending();
		}
	    } else {
		this->d->advanceAdmissionRevision(
		    BObolLodAdmissionRevisionDomain::AVAILABILITY);
		this->d->lodResidentAdmissionRevision =
		    admissionRevision;
	    }
	} else {
	    this->d->lodResidentAdmissionRevision =
		admissionRevision;
	}
    }
    this->scheduleResidentGrowthReallocationIfReady();
    /* Recovery can become enabled when its final capacity or handoff owner
     * retires without requesting another frame.  Treat the typed phase as a
     * level-triggered obligation so an idle renderer cannot strand it. */
    (void)this->completePointTriangleRecoveryIfReady();
    const int64_t refinementNow = bu_gettime();
    const bool refinementCooling =
	!options->forceTerminalLodRefinement &&
	this->d->lodRefinementNotBeforeMicroseconds > refinementNow &&
	!this->d->lodInteractionSession.active() && !localStatus.changed;
    if (refinementCooling)
	localStatus.hasMore = 1;
    else if (this->d->lodRefinementNotBeforeMicroseconds > 0 &&
	     refinementNow >=
		this->d->lodRefinementNotBeforeMicroseconds)
	this->d->lodRefinementNotBeforeMicroseconds = 0;

    /* Point aggregation controls a populated SoCADAssembly only.  A retained
     * structural placeholder assembly has no drawable CAD payload from which
     * to obtain its required timing witness; letting that barrier pause this
     * submission pass deadlocks discovery before the first mesh arrives. */
    if (this->d->lodPointQualityPhase.presentationPending() &&
	this->getActiveLodCadPayloadCount() == 0) {
	this->d->lodPointQualityPhase.completeCalibration();
	this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
	BObolViewLodState *presentationState = this->d->viewAttachment ?
	    this->d->viewAttachment->getViewLodState() : NULL;
	if (presentationState)
	    presentationState->setCadPresentationPointProxyPixelThreshold(1.0f);
	(void)this->d->confirmRetainedRecoveryPresentation(true);
    }

    if (this->d->lodAutoSubmit && lodPolicy.policy != BV_LOD_OFF &&
	lodPolicy.mesh_enabled && !refinementCooling) {
	const uint64_t submissionStarted = this->beginRenderTiming();
	(void)this->submitLodRequestsIfNeeded();
	/* A residency-drain cursor can retire in submitLodRequestsIfNeeded().
	 * Reevaluate its level-triggered handoff before deciding whether this pump
	 * has more work; otherwise no task, result, frame, or cursor may remain to
	 * wake the one retained reallocation which consumes the resident-growth
	 * debt. */
	this->scheduleResidentGrowthReallocationIfReady();
	const uint64_t submissionCompleted = this->beginRenderTiming();
	this->d->lastLodSubmissionTimeNanoseconds =
	    submissionCompleted > submissionStarted ?
	    submissionCompleted - submissionStarted : 0;
	/* A bounded compact-entry scan can have no worker tasks or results at
	 * the boundary between two chunks.  The submission cursor itself is
	 * pending work and must keep the frame pump alive.  Relying on
	 * submitLodRequests() calling markProgressiveWorkPending() is not
	 * sufficient: the status-based epilogue below otherwise clears that
	 * latch in the same advance.  The resulting scene refines only when an
	 * unrelated paint, checkpoint, or input event happens to arrive. */
	if (this->d->lodSubmissionPass.active())
	    localStatus.hasMore = 1;
	/* A capacity edge observed during an in-flight complete pass is
	 * deliberately coalesced rather than resetting that cursor.  If the pass
	 * completed in this same pump, retain one timer edge so the next pump can
	 * consume the newest admission revision and construct its sparse retry
	 * frontier. */
	if (this->d->lodService &&
	    this->d->lodService->residentMeshAdmissionRevision() !=
		this->d->lodResidentAdmissionRevision)
	    localStatus.hasMore = 1;
	/* A budget-limited quiet cut may intentionally request one unchanged
	 * frame to improve its throughput estimate.  Keep the host pump alive
	 * until completeRenderTiming() turns that presented probe into a new
	 * admission pass. */
	if (this->d->lodAdmissionEvidence.capacity().rescanAfterFrame())
	    localStatus.hasMore = 1;
    }

    /* Provider status describes database streaming.  Mesh refinement runs on
     * the controller-owned service and is independent of any one provider, so
     * sample it here.  Previously these fields remained zero and the common
     * host pump could declare a frame stable while a PoP task was still
     * queued or running. */
    if (this->d->lodService) {
	localStatus.pendingTasks +=
	    this->d->lodService->pendingTaskCountForGeneration(
		this->d->lodActiveGeneration);
	localStatus.pendingTasks +=
	    this->d->lodService->delayedTaskCountForGeneration(
		this->d->lodActiveGeneration);
	localStatus.inFlight +=
	    this->d->lodService->executingTaskCountForGeneration(
		this->d->lodActiveGeneration);
	localStatus.queuedResults +=
	    this->d->lodService->queuedResultCountForGeneration(
		this->d->lodActiveGeneration);
	localStatus.queuedCacheWrites +=
	    this->d->lodService->queuedCacheWriteCountForGeneration(
		this->d->lodActiveGeneration);
	/* Resident compaction planning is owner-thread work scoped to one view.
	 * lodCompactionPolicy below owns this controller's plan, and completed
	 * worker results wake only consumers which receive them.  The service's
	 * aggregate diagnostic also includes plans belonging to inactive sibling
	 * views; treating it as this controller's task kept an otherwise terminal
	 * view pumping forever while no local transition could retire that work. */
	localStatus.queuedResults +=
	    this->d->lodService->
		queuedResidentMeshCompactionResultCountForDiagnostics(
		    residentConsumerId);
    }

    const int pending_service_work =
	(localStatus.pendingTasks > 0 || localStatus.inFlight > 0 ||
	 localStatus.queuedResults > 0 || localStatus.queuedCacheWrites > 0) ?
	1 : 0;
    if (pending_service_work)
	localStatus.hasMore = 1;
    const int pending_lod_realization_work =
	this->d->lodService &&
	(this->d->lodService->activeTaskCountForGeneration(
	     this->d->lodActiveGeneration) > 0 ||
	 this->d->lodService->queuedResultCountForGeneration(
	     this->d->lodActiveGeneration) > 0 ||
	 this->d->lodService->queuedCacheWriteCountForGeneration(
	     this->d->lodActiveGeneration) > 0) ?
	1 : 0;

    /* Refinement and reclamation are separate phases.  A quiet view first
     * reaches its fast 1 px display target and may then enter the bounded
     * subpixel tiers when exact timing and memory witnesses prove sufficient
     * headroom.  Only after that terminal quality decision and a longer quiet
     * interval do we replace this consumer's complete demand snapshot and
     * trim shared CPU prefixes to the aggregate maximum. */
    const int64_t compactionNow = bu_gettime();
    BObolLodCompactionPolicy::Inputs compactionInputs;
    compactionInputs.automatic = this->d->lodAutoSubmit != FALSE;
    compactionInputs.interactive = this->d->lodInteractionSession.active() != FALSE;
    compactionInputs.coverageRequired =
	lodPolicy.policy != BV_LOD_OFF && lodPolicy.mesh_enabled;
    compactionInputs.coverageComplete =
	this->d->lodCoveragePolicy.coverageComplete();
    compactionInputs.coverageProgressPending = localStatus.hasMore != 0;
    compactionInputs.settlingPending =
	compactionInputs.coverageRequired &&
	(this->d->lodCoveragePolicy.active() ||
	 this->d->lodSubmissionPass.rescanPending() ||
	 this->d->lodRetainedPass.refinementPending() ||
	 this->d->lodRetainedPass.residencyPending() ||
	 this->d->lodAdmissionEvidence.capacity().rescanAfterFrame() ||
	 this->d->lodPresentationPolicy.handoffPending() ||
	 this->d->lodPresentationTransaction.barrierPending() ||
	 this->d->lodPresentationTransaction.publicationPending() ||
	 this->d->lodPointAdmissionFrame.pending() ||
	 this->d->lodPointQualityPhase.pending() ||
	 this->d->lodAvailabilityLedger.residentGrowthPending() ||
	 this->d->lodAdmissionEvidence.headroom().retryPending());
    compactionInputs.nowMicroseconds = compactionNow;
    compactionInputs.realizationPending =
	pending_lod_realization_work != 0;
    compactionInputs.submissionPending =
	this->d->lodSubmissionPass.active() != FALSE;
    compactionInputs.serviceAvailable = this->d->lodService != NULL;
    const BObolLodCompactionPolicy::Decision compactionDecision =
	this->d->lodCompactionPolicy.decide(compactionInputs);
    if (compactionDecision.retiredRequest &&
	this->d->lodAdmissionEvidence.resources().recoveryPending()) {
	/* The coverage producer disappeared (source erase, cancellation, LoD
	 * disable, or empty view).  Retiring compaction is the complete response
	 * to this pressure revision; leaving the resource edge unacknowledged
	 * otherwise advertises background work with no task, cursor, timer, or
	 * frame capable of discharging it. */
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::MARK_RESOURCE_RECOVERY_HANDLED);
    }
    if (compactionDecision.keepPumpAlive)
	localStatus.hasMore = 1;
    if (compactionDecision.admission ==
	BObolLodCompactionPolicy::Admission::PLAN) {
	    std::vector<BObolLodResidentDemand> demands;
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    const uint64_t demandRevision = viewState ?
		viewState->residentMeshDemandRevision() :
		this->d->lodViewRevision.value();
	    const SbBool continuePlanning =
		this->d->lodCompactionPolicy.continues(demandRevision) ?
		    TRUE : FALSE;
	    if (viewState && !continuePlanning)
		viewState->residentMeshDemands(demands);
	    SbBool planningComplete = FALSE;
	    const size_t queued =
		this->d->lodService->scheduleResidentMeshCompaction(
		    residentConsumerId,
		    demandRevision, demands,
		    &planningComplete);
	    this->d->lodCompactionPolicy.finishPlanning(
		planningComplete != FALSE, demandRevision, bu_gettime());
	    if (planningComplete &&
		this->d->lodAdmissionEvidence.resources().recoveryPending()) {
		/* One complete demand snapshot has now been admitted to the
		 * worker service for this pressure edge.  Persistent pressure can
		 * legitimately describe a visible working set larger than the
		 * configured limit; do not reschedule the identical compaction on
		 * every frame.  A falling/rising completed-frame edge starts a new
		 * recovery revision. */
		this->d->applyAdmissionEvidenceAction(
		    BObolLodAdmissionPlanner::EvidenceAction::MARK_RESOURCE_RECOVERY_HANDLED);
	    }
	    if (queued || !planningComplete)
		localStatus.hasMore = 1;
    }

    BObolLodPresentationTransaction::Inputs publicationInputs;
    publicationInputs.nowMicroseconds = bu_gettime();
    publicationInputs.observedRenderNanoseconds = std::max(
	this->d->lastSceneRenderTimeNanoseconds,
	this->d->lastRenderTimeNanoseconds);
    publicationInputs.interactive = this->d->lodInteractionSession.active() != FALSE;
    const bool publicationServiceProducer = this->d->lodService &&
	(this->d->lodService->activeTaskCountForGeneration(
	     this->d->lodActiveGeneration) > 0 ||
	 this->d->lodService->queuedResultCountForGeneration(
	     this->d->lodActiveGeneration) > 0);
    const bool publicationSubmissionPaused =
	BObolLodAdmissionPlanner::presentationPausesSubmission(
	    this->d->lodPointAdmissionFrame.pending(),
	    this->d->lodPointQualityPhase.presentationPending(),
	    this->d->lodAdmissionEvidence.capacity().rescanAfterFrame(),
	    controller_has_cad_presentation(this->d->viewAttachment));
    publicationInputs.streamIdle =
	!BObolLodProducerPolicy::canProduceGeometry(
	    this->d->lodSubmissionPass.active() != FALSE,
	    publicationSubmissionPaused,
	    providerPendingCount > 0, publicationServiceProducer);
    const BObolLodPresentationTransaction::Decision publicationDecision =
	this->d->lodPresentationTransaction.decide(publicationInputs);
    if (publicationDecision.keepPumpAlive)
	localStatus.hasMore = 1;
    if (publicationDecision.requestFrame) {
	this->requestRender("lod-result-batch");
	if (this->d->lodPresentationTransaction.barrierPending())
	    this->scheduleLodRefinementFrame("lod-result-batch");
    }

    /*
     * A refinement/calibration barrier is, by definition, waiting for a
     * completed presentation.  Result publication and Qt paint coalescing
     * can consume the request which originally installed the barrier before
     * the barrier's target render serial is reached.  Never allow that state
     * to advertise an idle pump with no frame request: reasserting the edge
     * is idempotent, and requestRender() wakes the host only on false->true.
     */
    if (this->hasPendingLodRefinementFrame()) {
	localStatus.hasMore = 1;
	/* A partial result installs its presentation barrier as soon as its
	 * immutable binding is applied, but sparse publication deliberately does
	 * not request that frame until the batch deadline.  Treating the barrier
	 * itself as a due frame collapses every batch to one result and restores
	 * the expensive whole-scene repaint stream this policy exists to avoid.
	 * The publication timer remains the liveness witness until decide()
	 * atomically replaces it with framePending(). */
	const SbBool publicationDeadlinePending =
	    this->d->lodPresentationTransaction.barrierPending() &&
	    this->d->lodPresentationTransaction.publicationAwaitingFrameRequest();
	/* Point aggregation calibration may remain logically pending while a
	 * source/result wave is still changing the population.  Those producers
	 * own their publication frames; forcing an additional unchanged render on
	 * each progressive timer tick serializes realization behind rasterization.
	 * When they become quiet, the still-pending latch reaches this branch on
	 * the next pump and requests its one explicit reusable replay. */
	const SbBool pointCalibrationWaitingForProducer =
	    this->d->lodPointQualityPhase.presentationPending() &&
	    BObolLodAdmissionPlanner::pointProducerOwnsCalibrationFrame(
		    this->d->lodSubmissionPass.active() != FALSE,
		    publicationSubmissionPaused,
		    providerPendingCount > 0,
		    pending_service_work != 0,
		    this->d->lodPresentationTransaction.
			publicationAwaitingFrameRequest());
	if (!publicationDeadlinePending &&
	    !pointCalibrationWaitingForProducer &&
	    !this->isRenderRequested()) {
	    this->requestRender("lod-refinement-pending");
	}
    }
    if (localStatus.hasMore)
	this->markProgressiveWorkPending();
    else
	this->clearProgressiveWorkPending();

    /* Pending background work needs a host timer, not a duplicate render of
     * unchanged pixels.  This distinction is especially important for
     * OSMesa, where merely repainting a multi-million-triangle stable cut can
     * consume seconds.  The Qt host pumps pending work independently and
     * actual result/cut changes install their own render requests. */
    if (localStatus.changed &&
	(nonLodPresentationChanged ||
	 this->d->lodPresentationTransaction.publicationFramePending()))
	this->requestRender("progressive-update");

    if (status)
	*status = localStatus;

    const uint64_t advanceCompleted = this->beginRenderTiming();
    this->d->lastProgressiveAdvanceTimeNanoseconds =
	advanceCompleted > advanceStarted ? advanceCompleted - advanceStarted : 0;
    return (localStatus.changed || localStatus.hasMore) ? 1 : 0;
}

void
BObolViewController::markProgressiveWorkPending(void)
{
    SbBool wakeEndpoint = FALSE;
    {
	std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
	if (!this->d->progressiveWorkPending) {
	    this->d->progressiveWorkPending = TRUE;
	    if (++this->d->hostWorkRevision == 0)
		++this->d->hostWorkRevision;
	    wakeEndpoint = TRUE;
	}
    }
    if (wakeEndpoint)
	this->notifyFrameRequest("progressive-work");
}

void
BObolViewController::clearProgressiveWorkPending(void)
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    if (!this->d->progressiveWorkPending)
	return;
    this->d->progressiveWorkPending = FALSE;
    if (++this->d->hostWorkRevision == 0)
	++this->d->hostWorkRevision;
}

SbBool
BObolViewController::hasProgressiveWorkPending(void) const
{
    std::lock_guard<std::mutex> lock(this->d->renderRequestMutex);
    return this->d->progressiveWorkPending;
}
