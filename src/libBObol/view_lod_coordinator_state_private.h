/*    V I E W _ L O D _ C O O R D I N A T O R _ S T A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_VIEW_LOD_COORDINATOR_STATE_PRIVATE_H
#define LIBBOBOL_VIEW_LOD_COORDINATOR_STATE_PRIVATE_H

#include "BObol/BPickDetail.h"
#include "BObol/BViewController.h"
#include "lod_coordinator_private.h"
#include "retained_allocation_private.h"

#include <limits>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <vector>

/**
 * Owner-thread progressive-display state.  This first extraction preserves
 * the exact field layout and algorithms formerly embedded in Impl; it adds
 * no allocation, virtual dispatch, or per-occurrence indirection.
 */
struct BObolLodCoordinator : BObolLodStateMachine {
    void resetRetainedAdmissionQualityProof(void)
    {
	lodRetainedAdmissionMaximumNormalizedError =
	    std::numeric_limits<double>::infinity();
	lodRetainedAdmissionMaximumProjectedErrorPixels =
	    std::numeric_limits<double>::infinity();
	lodRetainedAdmissionQualityViewRevision = 0;
	lodRetainedAdmissionQualityPolicyRevision = 0;
	lodRetainedAdmissionQualityPointProxyPixelThreshold =
	    std::numeric_limits<float>::quiet_NaN();
	lodRetainedAllocationCertificate = BObolRetainedAllocationResult();
    }

    void resetLodViewQualityHistory(void)
    {
	lodViewQualityHistory.reset();
	resetRetainedAdmissionQualityProof();
	lodViewQualityDomainRevision++;
	if (!lodViewQualityDomainRevision)
	    lodViewQualityDomainRevision = 1;
    }

    void resetRendererPerformanceEvidence(void)
    {
	/* A renderer epoch owns every timing-derived scalar, not just the
	 * exact-camera LRU.  Retained PoP buffers and occurrence cuts are scene
	 * data and deliberately survive this reset; they provide a coherent,
	 * conservative first frame from which the replacement renderer can be
	 * measured. */
	resetLodViewQualityHistory();
	lodInteractiveCalibratedRenderCostPerSecond = 0.0L;
	lodStableCalibratedRenderCostPerSecond = 0.0L;
	lastRenderTimeNanoseconds = 0;
	smoothedRenderTimeNanoseconds = 0;
	lodLastCadGpuSampleSerial = 0;
	lodLastCadGpuTimeNanoseconds = 0;
	lodLastCadGpuGeometryUploadBytes = 0;
	lodLastCadGpuGeometryUploadBytesValid = FALSE;
	lodLastRenderWasPreparedCadReplay = TRUE;
	lodLastRenderWasReusableCadPresentation = TRUE;
	lodStableBudgetBeforeInteraction = 0;
	lodStableBudgetBeforeInteractionValid = FALSE;
	lodPassAdmittedWork = FALSE;

	/* Capacity ceilings, overload witnesses, and handoff proofs were all
	 * established by the previous renderer.  Start a fresh bounded budget
	 * epoch, but keep the currently applied renderer-only ceiling and point
	 * cut for the first frame.  Their explicit recovery latches will relax
	 * those controls using measurements from the new renderer without
	 * exposing a potentially enormous retained suffix speculatively. */
	lodBudgetPolicy.reset();
	lodBudgetPolicy.requestRescanAfterFrame(true);
	lodHeadroomPolicy.cancelRetry();
	lodPresentationPolicy.reset();
	lodPointProxyCalibrationPolicy.reset();
	lodStablePointProxyCalibrationPending =
	    lodPresentationPointProxyPixelThreshold > 1.01f ? TRUE : FALSE;
	lodPointProxyTriangleRecoveryPending = FALSE;
	lodStaticOverscanActive =
	    lodInteractiveProgressiveCeiling >= 0 ? TRUE : FALSE;
	lodStaticOverscanLeapAvailable = FALSE;
	lodStaticOverscanRejected = FALSE;
	lodStaticOverscanRetryAfterPopulationChange = FALSE;
	lodDiscretePopulationTrialAvailable = FALSE;
	lodInteractiveCeilingFeedbackRenderSerial = 0;
	lodDeadlineSafeProgressiveCeiling = -1;
	lodDeadlineSafeViewRevision = 0;
	lodDeadlineSafePolicyRevision = 0;
	lodFrameObligation.reset();
	lodRefinementNotBeforeMicroseconds = 0;
	lodPublicationPolicy.reset();

	/* Delivery-rate telemetry is also renderer-specific.  It is protected
	 * independently because hosts may sample it while an endpoint property
	 * change is being applied on the owner thread. */
	{
	    std::lock_guard<std::mutex> lock(presentationTimingMutex);
	    lastPresentationTimestampNanoseconds = 0;
	    smoothedPresentationIntervalNanoseconds = 0;
	    displayedPresentationIntervalNanoseconds = 0;
	    lastInteractivePresentationTimestampNanoseconds = 0;
	    smoothedInteractivePresentationIntervalNanoseconds = 0;
	}
    }

    float quietAllocationTargetFps(void) const
    {
	/* The ordinary quiet cadence keeps cold streaming and multi-frame
	 * convergence responsive.  Once that phase has produced a terminal exact
	 * frame, an event-driven view may spend up to the separate hard static
	 * presentation deadline on a richer framebuffer which is then retained
	 * without redraw. */
	if (!lodStaticOverscanActive ||
	    !stablePresentationFrameDeadlineNanoseconds)
	    return lodStableTargetFps;
	const long double deadlineFps = 1000000000.0L /
	    static_cast<long double>(
		stablePresentationFrameDeadlineNanoseconds);
	if (!std::isfinite(deadlineFps) || deadlineFps <= 0.0L)
	    return lodStableTargetFps;
	const float bounded = static_cast<float>(deadlineFps);
	return lodStableTargetFps > 0.0f ?
	    std::min(lodStableTargetFps, bounded) : bounded;
    }

    float staticQualityTargetFps(void) const
    {
	/* Once the ordinary one-pixel frame is exact and terminal, a static
	 * event-driven view may use the separately configured hard frame deadline
	 * to test a finer pixel target.  The framebuffer is retained afterward,
	 * so requiring every such terminal quality step to meet the ordinary 20 Hz
	 * refinement cadence makes warm-cache convergence depend on one upload or
	 * command-preparation sample. */
	if (!stablePresentationFrameDeadlineNanoseconds)
	    return lodStableTargetFps;
	const long double deadlineFps = 1000000000.0L /
	    static_cast<long double>(
		stablePresentationFrameDeadlineNanoseconds);
	if (!std::isfinite(deadlineFps) || deadlineFps <= 0.0L)
	    return lodStableTargetFps;
	const float bounded = static_cast<float>(deadlineFps);
	return lodStableTargetFps > 0.0f ?
	    std::min(lodStableTargetFps, bounded) : bounded;
    }

    void reconcilePhase(Event event, size_t pendingCount = 0)
    {
	Inputs inputs;
	inputs.interactive = lodInteractive;
	inputs.gestureActive = lodGestureActive;
	/* Coverage is a prerequisite only while automatic mesh LoD is an
	 * active consumer.  A disabled/removed service has no mesh inventory to
	 * prove, and a planning cursor from the preceding source epoch must not
	 * be projected as concurrently compacting an invalidated inventory. */
	inputs.coverageComplete = lodCoveragePolicy.effectiveComplete();
	inputs.coverageActive = lodCoveragePolicy.effectiveActive();
	inputs.compacting = lodCompactionPolicy.planning() &&
	    lodCoveragePolicy.compactionAllowed();
	inputs.cpuMemoryPressure = lodResourcePolicy.cpuPressure();
	inputs.gpuMemoryPressure = lodResourcePolicy.gpuPressure();
	inputs.resourceRecoveryPending =
	    lodResourcePolicy.recoveryPending();
	inputs.generationActive = lodActiveGeneration != 0;
	inputs.settlingWork =
	    lodSubmissionPending || lodSubmissionRescanPending ||
	    lodRetainedRefinementPending ||
	    lodRetainedResidencyPending ||
	    lodBudgetPolicy.rescanAfterFrame() ||
	    lodPresentationPolicy.handoffPending() ||
	    lodFrameObligation.pending() ||
	    lodPublicationPolicy.pending() || pendingCount != 0 ||
	    lodBudgetPolicy.calibrationFramesRemaining() != 0 ||
	    lodDiscoveryPointProxyFramePending ||
	    lodStablePointProxyCalibrationPending ||
	    lodPointProxyTriangleRecoveryPending ||
	    lodHeadroomPolicy.retryPending();

	Witness witness;
	witness.viewEpoch = lodViewRevision.value();
	witness.policyEpoch = lodPolicyRevision.value();
	witness.renderSerial = renderCompletionSerial;
	witness.activeGeneration = lodActiveGeneration;
	witness.residentDemandRevision =
	    lodCompactionPolicy.demandRevision();
	witness.resourcePressureRevision = lodResourcePolicy.revision();
	witness.visibleCount = lodCoveragePolicy.visibleCount();
	witness.completedCount = lodCoveragePolicy.coveredCount();
	witness.pendingCount = pendingCount;
	(void)this->dispatch(event, inputs, witness);
    }

    uint64_t lastRenderTimeNanoseconds = 0;
    uint64_t smoothedRenderTimeNanoseconds = 0;
    uint64_t lastBackgroundRenderTimeNanoseconds = 0;
    uint64_t lastSceneRenderTimeNanoseconds = 0;
    /* Full CAD plan/atlas construction is a real latency problem, but not a
     * measurement of steady triangle throughput.  Only an unchanged prepared
     * replay may drive quiet-view capacity cuts. */
    SbBool lodLastRenderWasPreparedCadReplay = TRUE;
    SbBool lodLastRenderWasReusableCadPresentation = TRUE;
    uint64_t lodLastCadGpuSampleSerial = 0;
    uint64_t lodLastCadGpuTimeNanoseconds = 0;
    uint64_t lodLastCadGpuGeometryUploadBytes = 0;
    SbBool lodLastCadGpuGeometryUploadBytesValid = FALSE;
    uint64_t lastProgressiveAdvanceTimeNanoseconds = 0;
    uint64_t lastLodResultProcessingTimeNanoseconds = 0;
    uint64_t lastProgressiveProviderTimeNanoseconds = 0;
    uint64_t lastLodSubmissionTimeNanoseconds = 0;
    uint64_t lastPresentationSyncTimeNanoseconds = 0;
    uint64_t interactivePresentationFrameDeadlineNanoseconds = 40000000ULL;
    uint64_t stablePresentationFrameDeadlineNanoseconds = 100000000ULL;
    uint64_t interruptedPresentationFrameCount = 0;
    uint64_t consecutiveInterruptedPresentationFrames = 0;
    uint64_t lastInterruptedPresentationTimeNanoseconds = 0;
    uint64_t renderCompletionSerial = 0;
    mutable std::mutex presentationTimingMutex;
    uint64_t presentedFrameSerial = 0;
    uint64_t lastPresentationTimestampNanoseconds = 0;
    uint64_t smoothedPresentationIntervalNanoseconds = 0;
    uint64_t displayedPresentationIntervalNanoseconds = 0;
    uint64_t lastInteractivePresentationTimestampNanoseconds = 0;
    uint64_t smoothedInteractivePresentationIntervalNanoseconds = 0;
    std::vector<BObolProgressiveProviderRecord> progressiveProviders;
    uint64_t progressiveProviderNextToken = 1;
    /* Number of registered non-LoD geometry producers which reported work
     * remaining in the most recent complete provider pass.  Registration is
     * conservatively pending until the first pass proves otherwise. */
    size_t progressiveProviderPendingCount = 0;
    BObolProgressiveOptions defaultProgressiveOptions;
    BObolLodService *lodService = NULL;
    std::shared_ptr<BObolLodService> managedLodService;
    size_t managedLodWorkerCount = 0;
    uint64_t lodResultSubscriberId = 0;
    std::atomic<int> lodResultsPending {0};
    /* Timestamp the empty->non-empty worker-result edge.  Bulk warm-cache
     * streams use it to coalesce tiny timer-tick dribbles into bounded-latency
     * scene updates. */
    std::atomic<int64_t> lodResultsFirstReadyMicroseconds {0};
    BObolLodInventoryDeltaPolicy lodInventoryDeltaPolicy;
    SbBool lodAutoSubmit = FALSE;
    uint64_t lodActiveGeneration = 0;
    size_t lodSubmissionSourceIndex = 0;
    size_t lodSubmissionEntryOffset = 0;
    SoBRLDatabaseSource *lodSubmissionPlanSource = NULL;
    std::vector<size_t> lodSubmissionPlanEntries;
    SbBool lodSubmissionPlanValid = FALSE;
    SbBool lodSubmissionPlanRetainedAdmission = FALSE;
    /* Expensive stable-view minimax arithmetic is resumable.  The plan owns
     * no scene geometry and publishes nothing until its final epoch-checked
     * owner-thread commit. */
    std::shared_ptr<BObolRetainedAllocationTransaction>
	lodRetainedAllocationTransaction;
    SbBool lodSubmissionRetainedAdmissionMode = FALSE;
    double lodRetainedAdmissionMaximumNormalizedError =
	std::numeric_limits<double>::infinity();
    double lodRetainedAdmissionMaximumProjectedErrorPixels =
	std::numeric_limits<double>::infinity();
    uint64_t lodRetainedAdmissionQualityViewRevision = 0;
    uint64_t lodRetainedAdmissionQualityPolicyRevision = 0;
    float lodRetainedAdmissionQualityPointProxyPixelThreshold =
	std::numeric_limits<float>::quiet_NaN();
    BObolRetainedAllocationResult lodRetainedAllocationCertificate;
    SoBRLDatabaseSource *lodSubmissionVisibleCountSource = NULL;
    size_t lodSubmissionVisibleCount = 0;
    /* Large compact scenes are consumed in bounded GUI-thread windows.  A
     * coverage pass proves that every projected leaf has a useful structural
     * presentation before any of those leaves opens a PoP hierarchy.  The
     * following quality pass promotes view-significant leaves under the
     * aggregate budget, preventing early entries from starving later ones. */
    BObolLodCoveragePolicy lodCoveragePolicy;
    /*
     * A camera-scale refresh is not cold coverage.  Existing occurrence cuts
     * remain useful and a zoom must be allowed to retarget them, consume a
     * richer resident prefix, or request the missing cumulative suffix while
     * the bounded visibility scan continues.  New/source inventory coverage
     * keeps the stricter minimum-mesh-first rule below.
     */
    /*
     * The newest completed camera pass proved that every projected mesh
     * occurrence already has a retained mesh presentation.  While input is
     * active, richer worker results may then remain queued without exposing
     * boxes or missing newly visible geometry.
     */
    /* One byte per compact entry is deliberately used instead of a hash set:
     * entry indices are dense and source-stable, so a 150k-occurrence source
     * costs only 150 kB and an exact edit delta updates it in constant time. */
    BObolLodVisibilityCensus lodConvergenceCandidateCensus;
    /* Camera/geometric projection evidence is denser than the visibility
     * census but still bounded (roughly a few dozen bytes per occurrence).
     * It belongs to this view, never to the shared database source. */
    BObolLodProjectedDemandCache lodProjectedDemandCache;
    /*
     * Authoritative visible-mesh denominator from the newest completed
     * all-source coverage pass.  Per-source candidate counts are useful while
     * a bounded pass is still being assembled, but its final window used to
     * clear the accumulated coverage counters before convergence retained
     * their total.  A quality-only policy revision preserves this view proof;
     * camera and source-inventory revisions invalidate it.
    */
    mutable BObolLodConvergencePolicy lodConvergencePolicy;
    void clearLodSubmissionPlan(void)
    {
	lodSubmissionPlanSource = NULL;
	lodSubmissionPlanEntries.clear();
	lodSubmissionPlanValid = FALSE;
	lodSubmissionPlanRetainedAdmission = FALSE;
	lodRetainedAllocationTransaction.reset();
	lodSubmissionVisibleCountSource = NULL;
	lodSubmissionVisibleCount = 0;
    }
    void clearLodConvergenceCandidates(void)
    {
	lodConvergenceCandidateCensus.clear();
	lodCoveragePolicy.clearCompleteVisibleCount();
    }
    void resetLodConvergenceFraction(void)
    {
	lodConvergencePolicy.resetFraction();
    }
    static BObolLodVisibilityCensus::SourceKey lodConvergenceSourceKey(
	SoBRLDatabaseSource *source)
    {
	return source ? source->getCompactSourceRoutingId() : 0;
    }
    size_t lodConvergenceCandidateCensusTotal(void) const
    {
	return lodConvergenceCandidateCensus.total();
    }
    void publishExactLodConvergenceCandidateCount(void)
    {
	/* Preserve the coverage proof itself.  Exact deltas mutate only entries
	 * represented by the census; population-changing deltas clear it and run
	 * a complete pass instead. */
	if (lodCoveragePolicy.hasCompleteVisibleCount())
	    lodCoveragePolicy.setCompleteVisibleCount(
		lodConvergenceCandidateCensusTotal());
    }
    void setLodConvergenceCandidateCount(
	SoBRLDatabaseSource *source, size_t count)
    {
	if (!source)
	    return;
	lodConvergenceCandidateCensus.setCount(
	    lodConvergenceSourceKey(source), count);
	publishExactLodConvergenceCandidateCount();
    }
    void beginLodConvergenceCandidateCensus(
	SoBRLDatabaseSource *source, size_t entryCount)
    {
	if (!source)
	    return;
	lodConvergenceCandidateCensus.begin(
	    lodConvergenceSourceKey(source), entryCount);
    }
    void observeLodConvergenceCandidateVisibility(
	SoBRLDatabaseSource *source, size_t entryCount,
	const std::vector<std::pair<size_t, SbBool>> &observations)
    {
	if (!source)
	    return;
	for (const std::pair<size_t, SbBool> &observation : observations) {
	    lodConvergenceCandidateCensus.observe(
		lodConvergenceSourceKey(source), entryCount,
		observation.first, observation.second != FALSE);
	}
	if (lodConvergenceCandidateCensus.complete(
		lodConvergenceSourceKey(source)))
	    publishExactLodConvergenceCandidateCount();
    }
    void completeLodConvergenceCandidateCensus(
	SoBRLDatabaseSource *source, size_t observedVisibleCount)
    {
	if (!source)
	    return;
	/* The action's accumulated count is the coverage-pass authority.  The
	 * byte census records the same projection at entry granularity so later
	 * exact deltas can revise that total without a full traversal. */
	lodConvergenceCandidateCensus.complete(
	    lodConvergenceSourceKey(source), observedVisibleCount);
	publishExactLodConvergenceCandidateCount();
    }
    size_t lodConvergenceCandidateCount(void) const
    {
	if (lodCoveragePolicy.hasCompleteVisibleCount())
	    return lodCoveragePolicy.completeVisibleCount();
	return lodConvergenceCandidateCensusTotal();
    }
    bool hasCompleteLodConvergenceCandidateCensus(
	SoBRLDatabaseSource *source) const
    {
	return source && lodConvergenceCandidateCensus.complete(
	    lodConvergenceSourceKey(source));
    }
    bool isVisibleLodConvergenceCandidate(
	SoBRLDatabaseSource *source, size_t entryIndex) const
    {
	return source && lodConvergenceCandidateCensus.visible(
	    lodConvergenceSourceKey(source), entryIndex);
    }
    bool pointProxyAggregationApplicable(void) const
    {
	return BObolLodPointProxyCalibrationPolicy::applicable(
	    lodConvergenceCandidateCount());
    }
    SbBool lodSubmissionPending = FALSE;
    SbBool lodSubmissionRescanPending = FALSE;
    SbBool lodRetainedRefinementPending = FALSE;
    /* A retained minimax pass selected a level whose PoP suffix is not yet
     * resident.  This is provider work, not a performance-limited quality
     * observation, and must survive the coherent cut presentation barrier. */
    SbBool lodRetainedResidencyPending = FALSE;
    /* The last exact renderer frame observed structural boxes which the
     * predictive point classifier expected to collapse.  One bounded pass
     * must bypass that prediction and obtain their mesh presentations. */
    SbBool lodStructuralPresentationRepairPending = FALSE;
    /* Exact non-terminal box population which armed the current structural
     * repair.  Its completed budget witness controls how aggressively the
     * presentation-only point threshold aggregates the remaining tail. */
    size_t lodStructuralRepairTargetCount = 0;
    /* Accumulated across every bounded window of one scene pass.  Unlike the
     * general refinement-debt bit, this counts only visible box-to-first-mesh
     * admissions rejected by the finite scene allowance. */
    size_t lodMissingMeshBudgetBlockedCount = 0;
    SbBool lodRetainedRefinementCutAdvanced = FALSE;
    SbBool lodRetainedRefinementBudgetBlocked = FALSE;
    SbBool lodPassAdmittedWork = FALSE;
    SbBool forceTerminalLodRefinement = FALSE;
    SbBool lodSubmissionRefreshMissing = TRUE;
    int lodSubmissionReset = 0;
    BObolLodViewEpoch lodLastSubmittedViewRevision;
    BObolLodPolicyEpoch lodLastSubmittedPolicyRevision;
    struct LodSourceSnapshot {
	SoBRLDatabaseSource *source = NULL;
	struct db_i *database = NULL;
	BObolLodSourceRoutingId routingId;
	BObolLodInventoryEpoch inventoryRevision;
	SbString databaseId;
	SbString path;
	int drawMode = 0;
	int representationMode = 0;
	SbBool visible = FALSE;
	int lodBotThreshold = 0;
	uint32_t sourceRevision = 0;
	uint32_t inputsRevision = 0;

	bool sameIdentity(const LodSourceSnapshot &other) const
	{
	    return this->database == other.database &&
		this->routingId == other.routingId &&
		this->databaseId == other.databaseId &&
		this->path == other.path &&
		this->drawMode == other.drawMode &&
		this->representationMode == other.representationMode &&
		this->visible == other.visible &&
		this->lodBotThreshold == other.lodBotThreshold &&
		this->sourceRevision == other.sourceRevision &&
		this->inputsRevision == other.inputsRevision;
	}
    };
    std::vector<LodSourceSnapshot> lodLastSubmittedSources;
    SbBool lodSubmissionDeltaActive = FALSE;
    std::vector<SoBRLDatabaseSource *> lodSubmissionDeltaSources;
    std::vector<std::pair<SoBRLDatabaseSource *, std::vector<size_t>>>
	lodSubmissionDeltaPlans;
    BObolLodViewSnapshot lodViewSignature;
    BObolLodViewScaleSnapshot lodViewScaleSignature;
    SbBool lodViewSignatureValid = FALSE;

    /* Net scale, not the mere presence of a wheel event, decides whether a
     * quiet orthographic view must retarget every occurrence.  Keep the
     * signature and ready-view proof from the start of the continuous input
     * epoch.  A bracketed pose gesture may begin before an unbracketed wheel
     * epoch's debounce expires; in that case these values deliberately keep
     * describing the original stable view. */
    BObolLodViewScaleSnapshot lodInteractionStartScaleSignature;
    SbBool lodInteractionStartScaleSignatureValid = FALSE;
    SbBool lodInteractionStartedFromReadyView = FALSE;
    size_t lodStableBudgetBeforeInteraction = 0;
    SbBool lodStableBudgetBeforeInteractionValid = FALSE;
    /* A completed, fully covered orthographic pose may verify visibility,
     * selection, and resource pressure, but may not rewrite retained PoP
     * cuts merely because the interaction debounce ended. */
    SbBool lodRetainPoseOccurrenceCuts = FALSE;
    /* A settled retained camera epoch first records exact projected demand
     * for every visible occurrence, then performs one scene-budgeted
     * importance reallocation.  Keeping this as a census-completion latch
     * prevents stale previous-view metrics and makes the redistribution
     * explicitly one-shot. */
    SbBool lodRetainedImportanceCensusPending = FALSE;
    /* Worker completions extend immutable residency, but active occurrence
     * cuts remain an aggregate scene-budget decision.  Coalesce a whole
     * completion wave before invoking that allocator. */
    BObolLodResidentGrowthPolicy lodResidentGrowthPolicy;
    /* While the resident-growth policy drains further cache/provider work,
     * source actions may request immutable suffixes but must not rewrite the
     * visible occurrence allocation. */
    SbBool lodResidentGrowthResidencyDrainActive = FALSE;
    /* A resident-capacity revision reopens only occurrences whose provider
     * was denied by an older admission epoch.  Keep this distinct from the
     * ordinary unsatisfied-quality frontier: reclaimed bytes must not restart
     * a 150k-entry view census for a handful of denied assets. */
    SbBool lodResidentAdmissionRetryActive = FALSE;
    BObolLodViewEpoch lodViewRevision {1};
    BObolLodPolicyEpoch lodPolicyRevision {1};
    BObolLodViewDemandPolicy lodViewDemandPolicy;
    /* One scale-quality or exact stable-headroom probe may measure one
     * otherwise-unaffordable populated PoP transition across the entire
     * scene.  Individual submit actions are time-sliced and source-local, so
     * uniqueness must live here. */
    SbBool lodDiscretePopulationTrialAvailable = FALSE;
    int64_t lodLastViewChangeMicroseconds = 0;
    SbBool lodInteractive = FALSE;
    SbBool lodGestureActive = FALSE;
    /* A scale gesture whose camera has paused is still interactive, but it
     * should not remain a magnified copy of the motion cut until button-up.
     * One stable-budget quality probe may expose already resident detail;
     * the next camera event immediately returns to motion calibration. */
    /* A scale-changing frame may be followed by one bounded quality frame.
     * This is distinct from quiet/stable refinement: continuous wheel or
     * trackpad input must be able to expose a newly resident PoP suffix
     * without first ending the interaction. */
    /* Keep the complete held-button presentation through button-release
     * debounce.  This preserves the measured aggregate presentation, not the
     * richer retained per-occurrence cuts which may be hidden by the
     * renderer-wide motion ceiling. */
    SbBool lodReleaseCutFloorActive = FALSE;
    /* Renderer-only presentation continuity and motion-to-stable handoff.
     * The policy owns its snapshot/latch proof; this coordinator only applies
     * returned ceiling and point-proxy values to the retained view state. */
    BObolLodPresentationPolicy lodPresentationPolicy;
    /* Exact-view history survives camera epochs, but never source inventory,
     * service, or user quality-contract changes.  The separate domain makes
     * that lifetime explicit and prevents a coincidentally identical camera
     * from accepting evidence belonging to a different scene. */
    BObolLodViewQualityHistory lodViewQualityHistory;
    uint64_t lodViewQualityDomainRevision = 1;
    uint64_t lodSettleAfterRenderSerial = 0;
    BObolLodCompactionPolicy lodCompactionPolicy;
    /* Complete-frame resource pressure is retained independently of HUD
     * queries.  A pressure edge requests one safe compaction pass; if the
     * current visible working set itself exceeds a renderer limit, that pass
     * may finish in a stable memory-limited presentation rather than loop. */
    BObolLodResourcePolicy lodResourcePolicy;
    uint64_t lodResidentAdmissionRevision = 0;
    BObolLodFrameObligation lodFrameObligation;
    /* A result-publication frame may include one-time CPU/GPU buffer
     * construction.  It controls refinement pacing, but must not teach the
     * steady-state face-capacity estimator that every later frame costs the
     * same amount. */
    BObolLodPublicationPolicy lodPublicationPolicy;
    int64_t lodRefinementNotBeforeMicroseconds = 0;
    float lodTargetPixelError = 1.0f;
    float lodInteractiveTargetFps = 60.0f;
    float lodStableTargetFps = 20.0f;
    /* A terminal ordinary quiet frame may prove capacity for one or more
     * bounded static-image overscan allocations.  This never changes motion
     * calibration and is reset by every camera, policy, service, or input
     * epoch. */
    SbBool lodStaticOverscanActive = FALSE;
    /* A single-occurrence static handoff may combine its first two modest PoP
     * populations into one presentation.  Every later probe is one cut. */
    SbBool lodStaticOverscanLeapAvailable = FALSE;
    /* A hard static-frame miss is a capacity witness for the current view.
     * Preserve it across internal presentation/policy bookkeeping so an
     * ordinary repaint cannot reopen the same rejected quality staircase.
     * Genuine camera, user-policy, service, and cadence epochs clear it. */
    SbBool lodStaticOverscanRejected = FALSE;
    /* A hard-trial miss is valid only for the occurrence population which
     * produced it.  If point aggregation subsequently removes independent
     * draws, remember one bounded retry to arm after the current handoff has
     * finished.  Keeping this separate from the active/rejected bits retains
     * the 10 Hz reconciliation allowance without making the stale miss
     * permanent. */
    SbBool lodStaticOverscanRetryAfterPopulationChange = FALSE;
    BObolLodBudgetPolicy lodBudgetPolicy;
    long double lodInteractiveCalibratedRenderCostPerSecond = 0.0L;
    long double lodStableCalibratedRenderCostPerSecond = 0.0L;
    void seedInteractiveCalibrationFromStable(void)
    {
	if (!(lodStableCalibratedRenderCostPerSecond > 0.0L))
	    return;
	/* Stable retained rendering is a measured geometry-throughput witness,
	 * but motion also pays camera/cut/update overhead.  Carry half of that
	 * known capacity into a new gesture.  This is deliberately a floor, not
	 * a reset: an established faster interaction estimate is retained, while
	 * a 118-triangle underloaded frame can no longer become a permanent
	 * self-fulfilling capacity ceiling. */
	const long double seed =
	    lodStableCalibratedRenderCostPerSecond * 0.5L;
	lodInteractiveCalibratedRenderCostPerSecond = std::max(
	    lodInteractiveCalibratedRenderCostPerSecond, seed);
    }
    /* A late reusable frame may demonstrate headroom after the ordinary
     * bounded probe series has ended.  Remember the exact camera/policy and
     * budget already retried so each newly enlarged allowance gets one
     * admission pass without permitting an unchanged discrete-cut loop. */
    BObolLodHeadroomPolicy lodHeadroomPolicy;
    int lodInteractiveProgressiveCeiling = -1;
    /* Richest renderer-only ceiling which completed inside the hard deadline
     * for the current view/policy.  If a later allocator handoff probes the
     * unconstrained retained population and misses, return directly to this
     * proven cut instead of replaying every intervening PoP ordinal. */
    int lodDeadlineSafeProgressiveCeiling = -1;
    uint64_t lodDeadlineSafeViewRevision = 0;
    uint64_t lodDeadlineSafePolicyRevision = 0;
    /*
     * Current renderer-side small-occurrence aggregation cut.  Interaction
     * raises it immediately when measured frames miss their target.  A quiet
     * performance-limited view may also raise it when every retained mesh is
     * already at its minimum prefix and that irreducible population still
     * exceeds the calibrated scene budget.
     */
    float lodPresentationPointProxyPixelThreshold = 1.0f;
    BObolLodPointProxyCalibrationPolicy lodPointProxyCalibrationPolicy;
    /* Streaming discovery may temporarily aggregate tiny structural leaf
     * proxies more aggressively than the measured stable-view policy.  It is
     * renderer-only and is retired as soon as the producer inventory settles;
     * keeping a separate value prevents discovery pacing from contaminating
     * the terminal point-calibration bracket. */
    float lodDiscoveryPointProxyPixelThreshold = 1.0f;
    BObolLodPointProxyCalibrationPolicy lodDiscoveryPointProxyPolicy;
    SbBool lodDiscoveryPointProxyFramePending = FALSE;
    /*
     * A changed quiet aggregation threshold needs one measured presentation
     * before convergence is authoritative.  Keep this distinct from PoP
     * admission: no provider scan or retained-array mutation is required.
     */
    SbBool lodStablePointProxyCalibrationPending = FALSE;
    /* A temporary multi-pixel point cut may be relaxed only after reducible
     * PoP triangle detail has been compacted to the measured scene capacity.
     * This latch keeps convergence behind that retained-prefix pass. */
    SbBool lodPointProxyTriangleRecoveryPending = FALSE;
    /* A completed retained-recovery pass owns one authoritative triangle
     * allocation plan.  Point-threshold calibration may change only the
     * renderer-side small-occurrence classification afterward; it must not
     * submit the same no-op triangle recovery again.  The plan serial plus
     * cadAllocationPlanCoversCurrentPopulation() make this witness expire
     * automatically on a view, policy, occurrence, or resident-mesh change. */
    uint64_t lodPointProxyTriangleRecoverySaturatedPlanSerial = 0;
    uint64_t lodInteractiveCeilingFeedbackRenderSerial = 0;
    SbBool lodUseForcedCut = FALSE;
    int lodForcedCut = 0;
    uint64_t maxExactFullDetailFaceCount = 0;
    uint64_t maxExactFullDetailPointCount = 0;
    std::vector<BObolRtPickCache *> rtPickCaches;
    std::vector<SbString> rtPickCachePaths;
    std::vector<struct db_i *> rtPickCacheDatabases;
    std::vector<uint32_t> rtPickCacheSourceRevisions;
    SbBool meshResidencyBudgetEnabled = FALSE;
    size_t maxResidentMeshBytes = 0;
    SbBool meshResidencyEvictDisplayPayloads = TRUE;
    size_t lastMeshBudgetInitialResidentBytes = 0;
    size_t lastMeshBudgetFinalResidentBytes = 0;
    size_t lastMeshBudgetFreedFullDetailBytes = 0;
    size_t lastMeshBudgetFreedDisplayBytes = 0;
    unsigned int lastMeshBudgetVisitedMeshCount = 0;
    unsigned int lastMeshBudgetEvictedFullDetailMeshCount = 0;
    unsigned int lastMeshBudgetEvictedDisplayMeshCount = 0;
    unsigned int lastLodVisitedMeshCount = 0;
    unsigned int lastLodSubmittedTaskCount = 0;
    unsigned int lastLodUpdatedCutCount = 0;
    unsigned int lastLodSkippedMeshCount = 0;
    size_t lastLodResultCount = 0;
    unsigned int lastLodMatchedResultCount = 0;
    unsigned int lastLodAppliedResultCount = 0;
    unsigned int lastLodRejectedResultCount = 0;
    unsigned int lastLodUnmatchedResultCount = 0;
    SbString lastLodDiagnostics = SbString("");
};


#endif /* LIBBOBOL_VIEW_LOD_COORDINATOR_STATE_PRIVATE_H */
