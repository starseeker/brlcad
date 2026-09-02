/*                L O D _ C O N T R O L _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_CONTROL_PRIVATE_H
#define LIBBOBOL_LOD_CONTROL_PRIVATE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>

/*
 * One submission cursor may be active while a complete successor rescan is
 * owed.  Discovery can also pause the cursor while preserving that rescan
 * until a later inventory revision supplies more entries.  These four states
 * were formerly represented by two freely written booleans in three
 * translation units.
 */
class BObolLodSubmissionPass {
public:
    enum class State : uint8_t {
	IDLE = 0,
	ACTIVE,
	IDLE_RESCAN,
	ACTIVE_RESCAN
    };

    void activate(void)
    {
	this->stateValue = this->rescanPending() ?
	    State::ACTIVE_RESCAN : State::ACTIVE;
    }

    void deactivate(void)
    {
	this->stateValue = this->rescanPending() ?
	    State::IDLE_RESCAN : State::IDLE;
    }

    void setActive(bool active)
    {
	if (active)
	    this->activate();
	else
	    this->deactivate();
    }

    void requestRescan(void)
    {
	this->stateValue = this->active() ?
	    State::ACTIVE_RESCAN : State::IDLE_RESCAN;
    }

    void clearRescan(void)
    {
	this->stateValue = this->active() ? State::ACTIVE : State::IDLE;
    }

    void setRescanPending(bool pending)
    {
	if (pending)
	    this->requestRescan();
	else
	    this->clearRescan();
    }

    /* A new bounded pass never inherits rescan debt from its predecessor.
     * Keep this as one transition: separate clearRescan()/activate() writes
     * expose an invalid intermediate state to later maintenance edits. */
    void beginFresh(bool active = true)
    {
	this->stateValue = active ? State::ACTIVE : State::IDLE;
    }

    /* Discovery may finish a predecessor pass while compact inventory is
     * still open.  Once the producer closes, consume the preserved successor
     * obligation as the new active pass.  activate() is intentionally not
     * used here: ACTIVE_RESCAN means an active predecessor still owes a later
     * full pass, whereas this transition starts that full pass itself. */
    bool beginPendingRescan(void)
    {
	if (this->stateValue != State::IDLE_RESCAN)
	    return false;
	this->stateValue = State::ACTIVE;
	return true;
    }

    /* Retiring a completed transaction consumes both its activity and rescan
     * obligation.  deactivate() remains available only for deliberate pauses
     * which must preserve a pending rescan. */
    void retire(void)
    {
	this->stateValue = State::IDLE;
    }

    bool active(void) const
    {
	return this->stateValue == State::ACTIVE ||
	    this->stateValue == State::ACTIVE_RESCAN;
    }

    bool rescanPending(void) const
    {
	return this->stateValue == State::IDLE_RESCAN ||
	    this->stateValue == State::ACTIVE_RESCAN;
    }

    State state(void) const
    {
	return this->stateValue;
    }

private:
    State stateValue = State::IDLE;
};

/* A completed renderer presentation is useful capacity evidence only for the
 * exact view and source-quality domain which produced it.  Keep the key and
 * value inseparable so a partial controller update cannot authorize a stale
 * progressive cut. */
class BObolLodDeadlineSafePresentation {
public:
    void reset(void)
    {
	this->progressiveCeiling = -1;
	this->viewRevision = 0;
	this->qualityDomainRevision = 0;
    }

    void remember(int candidateCeiling, uint64_t candidateViewRevision,
	uint64_t candidateQualityDomainRevision)
    {
	if (candidateCeiling < 0)
	    return;
	if (this->validFor(candidateViewRevision,
		candidateQualityDomainRevision)) {
	    if (candidateCeiling > this->progressiveCeiling)
		this->progressiveCeiling = candidateCeiling;
	    return;
	}
	this->progressiveCeiling = candidateCeiling;
	this->viewRevision = candidateViewRevision;
	this->qualityDomainRevision = candidateQualityDomainRevision;
    }

    int ceilingFor(uint64_t candidateViewRevision,
	uint64_t candidateQualityDomainRevision) const
    {
	return this->validFor(candidateViewRevision,
		candidateQualityDomainRevision) ? this->progressiveCeiling : -1;
    }

private:
    bool validFor(uint64_t candidateViewRevision,
	uint64_t candidateQualityDomainRevision) const
    {
	return this->progressiveCeiling >= 0 &&
	    this->viewRevision == candidateViewRevision &&
	    this->qualityDomainRevision == candidateQualityDomainRevision;
    }

    int progressiveCeiling = -1;
    uint64_t viewRevision = 0;
    uint64_t qualityDomainRevision = 0;
};

/* The retained allocator's error bounds are meaningful only with the view,
 * policy, and point-proxy threshold which produced them.  The normalized
 * error also drives the next retained update, while the projected error is
 * exported to exact-view history only when the complete key still matches. */
class BObolLodRetainedQualityEvidence {
public:
    void reset(void)
    {
	this->maximumNormalizedErrorValue =
	    std::numeric_limits<double>::infinity();
	this->maximumProjectedErrorPixelsValue =
	    std::numeric_limits<double>::infinity();
	this->viewRevision = 0;
	this->policyRevision = 0;
	this->pointProxyPixelThreshold =
	    std::numeric_limits<float>::quiet_NaN();
    }

    void remember(double maximumNormalizedError,
	double maximumProjectedErrorPixels, uint64_t currentViewRevision,
	uint64_t currentPolicyRevision, float currentPointProxyPixelThreshold)
    {
	this->maximumNormalizedErrorValue = maximumNormalizedError;
	this->maximumProjectedErrorPixelsValue = maximumProjectedErrorPixels;
	this->viewRevision = currentViewRevision;
	this->policyRevision = currentPolicyRevision;
	this->pointProxyPixelThreshold = currentPointProxyPixelThreshold;
    }

    double maximumNormalizedError(void) const
    {
	return this->maximumNormalizedErrorValue;
    }

    double maximumProjectedErrorPixelsFor(uint64_t currentViewRevision,
	uint64_t currentPolicyRevision,
	float currentPointProxyPixelThreshold) const
    {
	if (this->viewRevision != currentViewRevision ||
	    this->policyRevision != currentPolicyRevision ||
	    !std::isfinite(this->pointProxyPixelThreshold) ||
	    !std::isfinite(currentPointProxyPixelThreshold) ||
	    std::fabs(this->pointProxyPixelThreshold -
		currentPointProxyPixelThreshold) > POINT_PROXY_THRESHOLD_EPSILON)
	    return std::numeric_limits<double>::infinity();
	return this->maximumProjectedErrorPixelsValue;
    }

private:
    static constexpr float POINT_PROXY_THRESHOLD_EPSILON = 1.0e-6f;

    double maximumNormalizedErrorValue =
	std::numeric_limits<double>::infinity();
    double maximumProjectedErrorPixelsValue =
	std::numeric_limits<double>::infinity();
    uint64_t viewRevision = 0;
    uint64_t policyRevision = 0;
    float pointProxyPixelThreshold =
	std::numeric_limits<float>::quiet_NaN();
};

/* Renderer timing, point threshold, and retained-upload observations form one
 * epoch-local sample.  Time and upload/reuse bookkeeping must advance
 * together; otherwise an old duration can be paired with a new threshold or
 * a first upload can be mistaken for sustainable retained replay. */
class BObolLodRendererPerformanceEvidence {
public:
    static bool pointThresholdMatches(float first, float second)
    {
	return std::isfinite(first) && std::isfinite(second) &&
	    std::fabs(first - second) <= POINT_THRESHOLD_SAMPLE_EPSILON;
    }

    void reset(void)
    {
	this->lastGpuSampleSerial = 0;
	this->lastGpuTimeNanosecondsValue = 0;
	this->lastGeometryUploadBytes.reset();
	this->lastPresentationReusable = false;
	this->lastCadTimeNanosecondsValue = 0;
	this->lastCadRenderCostValue.reset();
	this->lastCadPointThreshold =
	    std::numeric_limits<float>::quiet_NaN();
	this->lastStructuralTimeNanosecondsValue = 0;
	this->lastStructuralRenderCostValue.reset();
	this->lastStructuralPointThreshold =
	    std::numeric_limits<float>::quiet_NaN();
    }

    bool acceptGpuMeasurement(uint64_t sampleSerial,
	uint64_t elapsedNanoseconds)
    {
	if (sampleSerial == this->lastGpuSampleSerial)
	    return false;
	this->lastGpuSampleSerial = sampleSerial;
	this->lastGpuTimeNanosecondsValue = elapsedNanoseconds;
	return true;
    }

    void noteCadPresentation(uint64_t elapsedNanoseconds,
	float pointProxyPixelThreshold,
	const std::optional<size_t> &presentedRenderCost, bool preparedReplay,
	const std::optional<uint64_t> &geometryUploadBytes)
    {
	const bool uploadChanged = geometryUploadBytes &&
	    (!this->lastGeometryUploadBytes ||
	     *geometryUploadBytes != *this->lastGeometryUploadBytes);
	if (geometryUploadBytes)
	    this->lastGeometryUploadBytes = geometryUploadBytes;
	this->lastPresentationReusable = preparedReplay ||
	    (geometryUploadBytes && !uploadChanged);
	this->lastCadTimeNanosecondsValue = elapsedNanoseconds;
	this->lastCadRenderCostValue = presentedRenderCost;
	this->lastCadPointThreshold = pointProxyPixelThreshold;
    }

    uint64_t lastGpuTimeNanoseconds(void) const
    {
	return this->lastGpuTimeNanosecondsValue;
    }

    bool reusableCadPresentationAt(float pointProxyPixelThreshold) const
    {
	return this->lastPresentationReusable &&
	    this->lastCadTimeNanosecondsValue > 0 &&
	    pointThresholdMatches(this->lastCadPointThreshold,
		pointProxyPixelThreshold);
    }

    uint64_t cadPresentationNanosecondsAt(
	float pointProxyPixelThreshold) const
    {
	return this->reusableCadPresentationAt(pointProxyPixelThreshold) ?
	    this->lastCadTimeNanosecondsValue : 0;
    }

    uint64_t cadPresentationNanosecondsAt(
	float pointProxyPixelThreshold, size_t presentedRenderCost) const
    {
	return this->reusableCadPresentationAt(pointProxyPixelThreshold) &&
	    this->lastCadRenderCostValue &&
	    *this->lastCadRenderCostValue == presentedRenderCost ?
		this->lastCadTimeNanosecondsValue : 0;
    }

    /* A point/box classifier frame contains no managed mesh population, so
     * it must not calibrate sustainable mesh throughput.  It is nevertheless
     * the only measured cost available when structural repair decides how
     * many boxes may be replaced in its first mesh wave.  Keep that evidence
     * in a separate channel so the repair planner can consume it without
     * teaching the ordinary retained allocator from a transient frame. */
    void noteStructuralPresentation(uint64_t elapsedNanoseconds,
	float pointProxyPixelThreshold, size_t presentedRenderCost)
    {
	this->lastStructuralTimeNanosecondsValue = elapsedNanoseconds;
	this->lastStructuralRenderCostValue = presentedRenderCost;
	this->lastStructuralPointThreshold = pointProxyPixelThreshold;
    }

    uint64_t structuralPresentationNanosecondsAt(
	float pointProxyPixelThreshold) const
    {
	return this->lastStructuralTimeNanosecondsValue > 0 &&
	    pointThresholdMatches(this->lastStructuralPointThreshold,
		pointProxyPixelThreshold) ?
	    this->lastStructuralTimeNanosecondsValue : 0;
    }

    uint64_t structuralPresentationNanosecondsAt(
	float pointProxyPixelThreshold, size_t presentedRenderCost) const
    {
	return this->lastStructuralTimeNanosecondsValue > 0 &&
	    this->lastStructuralRenderCostValue &&
	    *this->lastStructuralRenderCostValue == presentedRenderCost &&
	    pointThresholdMatches(this->lastStructuralPointThreshold,
		pointProxyPixelThreshold) ?
	    this->lastStructuralTimeNanosecondsValue : 0;
    }

private:
    static constexpr float POINT_THRESHOLD_SAMPLE_EPSILON = 0.01f;

    uint64_t lastGpuSampleSerial = 0;
    uint64_t lastGpuTimeNanosecondsValue = 0;
    std::optional<uint64_t> lastGeometryUploadBytes;
    bool lastPresentationReusable = false;
    uint64_t lastCadTimeNanosecondsValue = 0;
    std::optional<size_t> lastCadRenderCostValue;
    float lastCadPointThreshold =
	std::numeric_limits<float>::quiet_NaN();
    uint64_t lastStructuralTimeNanosecondsValue = 0;
    std::optional<size_t> lastStructuralRenderCostValue;
    float lastStructuralPointThreshold =
	std::numeric_limits<float>::quiet_NaN();
};

/* Immutable-in-meaning configuration retained while a bounded submission
 * cursor is active.  Mode, refresh behavior, and reset intent used to be
 * independent fields, allowing a resumed cursor to inherit only part of the
 * request which created it. */
class BObolLodSubmissionIntent {
public:
    enum class Mode : uint8_t {
	ORDINARY = 0,
	RETAINED_ADMISSION
    };

    void configure(bool refreshMissing, bool resetExisting)
    {
	this->refreshMissingValue = refreshMissing;
	this->resetExistingValue = resetExisting;
    }

    void setRetainedAdmission(bool enabled)
    {
	this->modeValue = enabled ?
	    Mode::RETAINED_ADMISSION : Mode::ORDINARY;
    }

    void reset(void)
    {
	this->modeValue = Mode::ORDINARY;
	this->refreshMissingValue = true;
	this->resetExistingValue = false;
    }

    bool retainedAdmission(void) const
    {
	return this->modeValue == Mode::RETAINED_ADMISSION;
    }

    bool refreshMissing(void) const
    {
	return this->refreshMissingValue;
    }

    bool resetExisting(void) const
    {
	return this->resetExistingValue;
    }

    Mode mode(void) const
    {
	return this->modeValue;
    }

private:
    Mode modeValue = Mode::ORDINARY;
    bool refreshMissingValue = true;
    bool resetExistingValue = false;
};

/* Retained-presentation proof for one interaction-to-quiet handoff.
 *
 * A ready view supplies two independent facts.  An orthographic pose-only
 * change may reuse its exact occurrence cuts because projected depth does not
 * change demand.  Any ready view with retained meshes, including a zoom, may
 * start capacity calibration at the event-driven static deadline: its
 * coherent presentation is already useful and must not first be discarded to
 * rediscover the preferred redraw cadence.  The current-view visibility
 * census remains a separate readiness obligation. */
class BObolLodRetainedViewContinuity {
public:
    void beginQuiet(bool orthographic, bool scaleChanged,
	bool startedFromReadyView, bool retainedMeshPayloads)
    {
	const bool retainedReadyPresentation =
	    startedFromReadyView && retainedMeshPayloads;
	this->retainOccurrenceCutsValue = retainedReadyPresentation &&
	    orthographic && !scaleChanged;
	this->startCapacityAtStaticValue = retainedReadyPresentation;
    }

    void beginStableMutation(bool retainedMeshPayloads)
    {
	/* A visibility-only hierarchy edit keeps the camera and immutable mesh
	 * population unchanged.  Its current occurrence cuts are therefore the
	 * coherent baseline, and any replacement capacity question belongs to the
	 * event-driven static deadline rather than the interactive redraw target. */
	this->retainOccurrenceCutsValue = retainedMeshPayloads;
	this->startCapacityAtStaticValue = retainedMeshPayloads;
    }

    void clearHandoff(void)
    {
	this->retainOccurrenceCutsValue = false;
	this->startCapacityAtStaticValue = false;
    }

    bool retainOccurrenceCuts(void) const
    {
	return this->retainOccurrenceCutsValue;
    }

    bool startCapacityAtStatic(void) const
    {
	return this->startCapacityAtStaticValue;
    }

    void deferVisibilityCensus(void)
    {
	this->visibilityCensusDeferredValue = true;
    }

    void completeVisibilityCensus(void)
    {
	this->visibilityCensusDeferredValue = false;
    }

    bool visibilityCensusDeferred(void) const
    {
	return this->visibilityCensusDeferredValue;
    }

    void reset(void)
    {
	this->clearHandoff();
	this->visibilityCensusDeferredValue = false;
    }

private:
    bool retainOccurrenceCutsValue = false;
    bool startCapacityAtStaticValue = false;
    bool visibilityCensusDeferredValue = false;
};

/* Complete annotations produced by one bounded retained-admission pass.
 * These facts may coexist, but they share one pass lifetime.  Keeping their
 * mutation private prevents a successor pass from inheriting a convenient
 * subset of the preceding pass's cut, residency, or budget evidence. */
class BObolLodRetainedPassAnnotations {
public:
    void reset(void)
    {
	this->refinementPendingValue = false;
	this->residencyPendingValue = false;
	this->cutAdvancedValue = false;
	this->budgetBlockedValue = false;
	this->admittedWorkValue = false;
	this->missingMeshBudgetBlockedCountValue = 0;
    }

    void retireRefinement(void)
    {
	this->refinementPendingValue = false;
	this->cutAdvancedValue = false;
	this->budgetBlockedValue = false;
    }

    void noteRefinementPending(void) { this->refinementPendingValue = true; }
    void clearRefinementPending(void)
    {
	this->refinementPendingValue = false;
    }
    bool refinementPending(void) const { return this->refinementPendingValue; }

    void noteResidencyPending(void) { this->residencyPendingValue = true; }
    void clearResidencyPending(void) { this->residencyPendingValue = false; }
    bool residencyPending(void) const { return this->residencyPendingValue; }

    void noteCutAdvanced(void) { this->cutAdvancedValue = true; }
    bool cutAdvanced(void) const { return this->cutAdvancedValue; }

    void noteBudgetBlocked(void) { this->budgetBlockedValue = true; }
    bool budgetBlocked(void) const { return this->budgetBlockedValue; }

    void noteAdmittedWork(void) { this->admittedWorkValue = true; }
    void clearAdmittedWork(void) { this->admittedWorkValue = false; }
    bool admittedWork(void) const { return this->admittedWorkValue; }

    void addMissingMeshBudgetBlocked(size_t count)
    {
	const size_t maximum = (std::numeric_limits<size_t>::max)();
	this->missingMeshBudgetBlockedCountValue =
	    count > maximum - this->missingMeshBudgetBlockedCountValue ?
		maximum : this->missingMeshBudgetBlockedCountValue + count;
    }
    void clearMissingMeshBudgetBlocked(void)
    {
	this->missingMeshBudgetBlockedCountValue = 0;
    }
    size_t missingMeshBudgetBlockedCount(void) const
    {
	return this->missingMeshBudgetBlockedCountValue;
    }

private:
    bool refinementPendingValue = false;
    bool residencyPendingValue = false;
    bool cutAdvancedValue = false;
    bool budgetBlockedValue = false;
    bool admittedWorkValue = false;
    size_t missingMeshBudgetBlockedCountValue = 0;
};

/* Exclusive owner-thread gate for one interrupted retained traversal.  The
 * replay is either absent or waiting for its exact successor frame; source
 * publication cannot infer or clear this state indirectly. */
class BObolLodInterruptedPresentationReplay {
public:
    enum class State : uint8_t {
	IDLE = 0,
	AWAITING_FRAME
    };

    void arm(void) { this->stateValue = State::AWAITING_FRAME; }
    void retire(void) { this->stateValue = State::IDLE; }
    void reset(void) { this->retire(); }
    bool pending(void) const
    {
	return this->stateValue == State::AWAITING_FRAME;
    }
    State state(void) const { return this->stateValue; }

private:
    State stateValue = State::IDLE;
};

/* Level-triggered proof obligation for the current CAD presentation and its
 * renderer controls.  Unlike BObolLodInterruptedPresentationReplay this latch
 * does not freeze scene publication: a producer may supersede the incomplete
 * target, but some exact successor frame must still commit before the view can
 * become terminal. */
class BObolLodExactPresentationFrame {
public:
    enum class State : uint8_t {
	CURRENT = 0,
	REQUEST_REQUIRED,
	AWAITING_FRAME
    };

    /* A semantic mutation supersedes every frame which began at or before
     * the mutation boundary.  The timestamp is captured from the same steady
     * clock as presentation timing, so an older in-flight frame cannot retire
     * a newer style or geometry obligation merely because it completed last. */
    void require(uint64_t mutationBoundaryNanoseconds)
    {
	this->stateValue = State::REQUEST_REQUIRED;
	this->mutationNanoseconds = std::max(
	    this->mutationNanoseconds, mutationBoundaryNanoseconds);
    }
    void noteFrameRequested(void)
    {
	if (this->stateValue == State::REQUEST_REQUIRED)
	    this->stateValue = State::AWAITING_FRAME;
    }
    void noteRequestRetired(void)
    {
	if (this->stateValue == State::AWAITING_FRAME)
	    this->stateValue = State::REQUEST_REQUIRED;
    }
    bool confirm(uint64_t frameStartedNanoseconds)
    {
	if (this->stateValue != State::CURRENT &&
	    frameStartedNanoseconds <= this->mutationNanoseconds)
	    return false;
	this->stateValue = State::CURRENT;
	this->mutationNanoseconds = 0;
	return true;
    }
    void reset(void)
    {
	this->stateValue = State::CURRENT;
	this->mutationNanoseconds = 0;
    }
    bool pending(void) const
    {
	return this->stateValue != State::CURRENT;
    }
    bool requestPending(void) const
    {
	return this->stateValue == State::REQUEST_REQUIRED;
    }
    bool framePending(void) const
    {
	return this->stateValue == State::AWAITING_FRAME;
    }
    static bool recoveryRequired(bool automatic, bool hasPresentation,
	bool exactPresentationCurrent, bool foregroundWorkPending)
    {
	return automatic && hasPresentation && !exactPresentationCurrent &&
	    !foregroundWorkPending;
    }
    State state(void) const { return this->stateValue; }

private:
    State stateValue = State::CURRENT;
    uint64_t mutationNanoseconds = 0;
};

/* One classifier-changing point threshold must be presented before source
 * admission may consume it.  This is a frame obligation, not a free boolean
 * or a stable point-quality calibration phase. */
class BObolLodPointAdmissionFrame {
public:
    enum class State : uint8_t {
	IDLE = 0,
	AWAITING_FRAME
    };

    void request(void) { this->stateValue = State::AWAITING_FRAME; }
    void setPending(bool pending)
    {
	this->stateValue = pending ? State::AWAITING_FRAME : State::IDLE;
    }
    void retire(void) { this->stateValue = State::IDLE; }
    void reset(void) { this->retire(); }
    bool pending(void) const
    {
	return this->stateValue == State::AWAITING_FRAME;
    }
    State state(void) const { return this->stateValue; }

private:
    State stateValue = State::IDLE;
};

/* One exact non-terminal structural-box frontier and the per-occurrence cost
 * reservation derived for its bounded repair pass.  Point relaxation is one
 * finite transaction composed of bounded admission batches: each batch
 * replaces part of the occurrence population the candidate threshold would
 * expose, an exact frame proves the remaining population is strictly smaller,
 * and only a zero remainder permits publication.  Encoding those stages and
 * that decreasing rank here prevents a partial Boolean update from skipping
 * or repeating a phase. */
class BObolLodStructuralRepair {
public:
    enum class Representation : uint8_t {
	MESH = 0,
	TERMINAL_PROXY
    };

    enum class PointRelaxationState : uint8_t {
	INACTIVE = 0,
	ADMISSION_PENDING,
	PRESENTATION_PENDING,
	FINALIZATION_PENDING
    };

    void begin(size_t frontierCount)
    {
	this->frontierCountValue = frontierCount;
	this->representationValue = Representation::MESH;
	this->coverageCostReservationValue = 0;
	this->pointRelaxationTargetValue = 0.0f;
	this->pointRelaxationRemainingRankValue = 0;
	this->pointRelaxationStateValue = PointRelaxationState::INACTIVE;
    }

    void beginTerminalProxy(size_t frontierCount)
    {
	this->frontierCountValue = frontierCount;
	this->representationValue = Representation::TERMINAL_PROXY;
	this->coverageCostReservationValue = 0;
	this->pointRelaxationTargetValue = 0.0f;
	this->pointRelaxationRemainingRankValue = 0;
	this->pointRelaxationStateValue = PointRelaxationState::INACTIVE;
    }

    bool beginPointRelaxation(size_t batchCount, size_t remainingRank,
	float targetThreshold)
    {
	if (!batchCount || batchCount > remainingRank)
	    return false;
	this->frontierCountValue = batchCount;
	this->pointRelaxationRemainingRankValue = remainingRank;
	this->representationValue = Representation::MESH;
	this->coverageCostReservationValue = 0;
	this->pointRelaxationTargetValue = targetThreshold;
	this->pointRelaxationStateValue =
	    PointRelaxationState::ADMISSION_PENDING;
	return true;
    }

    void reset(void)
    {
	this->frontierCountValue = 0;
	this->representationValue = Representation::MESH;
	this->coverageCostReservationValue = 0;
	this->pointRelaxationTargetValue = 0.0f;
	this->pointRelaxationRemainingRankValue = 0;
	this->pointRelaxationStateValue = PointRelaxationState::INACTIVE;
    }

    void completePointRelaxationAdmission(void)
    {
	if (this->pointRelaxationStateValue !=
		PointRelaxationState::ADMISSION_PENDING)
	    return;
	this->frontierCountValue = 0;
	this->coverageCostReservationValue = 0;
	this->pointRelaxationStateValue =
	    PointRelaxationState::PRESENTATION_PENDING;
    }

    bool retryPointRelaxationAdmission(size_t batchCount,
	size_t remainingRank)
    {
	if (this->pointRelaxationStateValue !=
		PointRelaxationState::FINALIZATION_PENDING ||
	    !batchCount || batchCount > remainingRank ||
	    remainingRank >= this->pointRelaxationRemainingRankValue)
	    return false;
	this->frontierCountValue = batchCount;
	this->pointRelaxationRemainingRankValue = remainingRank;
	this->coverageCostReservationValue = 0;
	this->pointRelaxationStateValue =
	    PointRelaxationState::ADMISSION_PENDING;
	return true;
    }

    static bool pointRelaxationDomainChanged(bool viewOrPolicyChanged,
	bool visibilityChanged, bool sourceCoverageInvalidated)
    {
	/* Resident-result publication changes inventory revision while fulfilling
	 * this transaction.  It is expected progress, not a new candidate domain.
	 * Only changes which alter the selected occurrence population or its
	 * projection invalidate the private threshold. */
	return viewOrPolicyChanged || visibilityChanged ||
	    sourceCoverageInvalidated;
    }

    void notePointRelaxationPresented(void)
    {
	if (this->pointRelaxationStateValue ==
		PointRelaxationState::PRESENTATION_PENDING)
	    this->pointRelaxationStateValue =
		PointRelaxationState::FINALIZATION_PENDING;
    }

    void cancelPointRelaxation(void)
    {
	if (this->pointRelaxationStateValue ==
		PointRelaxationState::INACTIVE)
	    return;
	this->pointRelaxationTargetValue = 0.0f;
	this->pointRelaxationRemainingRankValue = 0;
	this->pointRelaxationStateValue = PointRelaxationState::INACTIVE;
	/* A submitted selective frontier still owns finite useful mesh work.
	 * Let it finish as ordinary structural admission, but an admitted
	 * candidate awaiting presentation has no work left to preserve. */
	if (!this->active()) {
	    this->frontierCountValue = 0;
	    this->coverageCostReservationValue = 0;
	}
    }

    bool active(void) const
    {
	return this->frontierCountValue > 0;
    }

    bool terminalProxy(void) const
    {
	return this->active() &&
	    this->representationValue == Representation::TERMINAL_PROXY;
    }

    size_t frontierCount(void) const
    {
	return this->frontierCountValue;
    }

    bool pointRelaxationPending(void) const
    {
	return this->pointRelaxationStateValue !=
	    PointRelaxationState::INACTIVE;
    }

    float pointRelaxationTarget(void) const
    {
	return this->pointRelaxationTargetValue;
    }

    size_t pointRelaxationRemainingRank(void) const
    {
	return this->pointRelaxationRemainingRankValue;
    }

    bool pointRelaxationPresentationPending(void) const
    {
	return this->pointRelaxationStateValue ==
	    PointRelaxationState::PRESENTATION_PENDING;
    }

    PointRelaxationState pointRelaxationState(void) const
    {
	return this->pointRelaxationStateValue;
    }

    size_t coverageCostReservation(void) const
    {
	return this->coverageCostReservationValue;
    }

    void clearCoverageCostReservation(void)
    {
	this->coverageCostReservationValue = 0;
    }

    void reserveCoverageCost(size_t reservation)
    {
	if (this->active() && !this->coverageCostReservationValue)
	    this->coverageCostReservationValue = reservation;
    }

private:
    size_t frontierCountValue = 0;
    Representation representationValue = Representation::MESH;
    size_t coverageCostReservationValue = 0;
    float pointRelaxationTargetValue = 0.0f;
    size_t pointRelaxationRemainingRankValue = 0;
    PointRelaxationState pointRelaxationStateValue =
	PointRelaxationState::INACTIVE;
};

/* Concurrent one-shot planning obligations owned by the current evidence
 * epoch.  These are work requests, not outcome flags: each must be explicitly
 * requested and retired by its sole effect owner. */
class BObolLodPlanningObligations {
private:
    enum class Work : uint8_t {
	RETAINED_IMPORTANCE_CENSUS = 1u << 0,
	RESIDENT_ADMISSION_RETRY = 1u << 1,
	EXACT_VISIBILITY_REALLOCATION = 1u << 2
    };

    static constexpr uint8_t bit(Work work)
    {
	return static_cast<uint8_t>(work);
    }

    void request(Work work)
    {
	this->workValue |= bit(work);
    }

    void retire(Work work)
    {
	this->workValue &= static_cast<uint8_t>(~bit(work));
    }

    void setPending(Work work, bool pending)
    {
	if (pending)
	    this->request(work);
	else
	    this->retire(work);
    }

    void clearAll(void)
    {
	this->workValue = 0;
    }

    bool pending(Work work) const
    {
	return (this->workValue & bit(work)) != 0;
    }

public:
    void reset(void)
    {
	this->clearAll();
    }

    void requestImportanceCensus(void)
    {
	this->request(Work::RETAINED_IMPORTANCE_CENSUS);
    }

    void retireImportanceCensus(void)
    {
	this->retire(Work::RETAINED_IMPORTANCE_CENSUS);
    }

    void setImportanceCensus(bool pending)
    {
	this->setPending(Work::RETAINED_IMPORTANCE_CENSUS, pending);
    }

    bool importanceCensusPending(void) const
    {
	return this->pending(Work::RETAINED_IMPORTANCE_CENSUS);
    }

    void setResidentAdmissionRetry(bool pending)
    {
	this->setPending(Work::RESIDENT_ADMISSION_RETRY, pending);
    }

    void retireResidentAdmissionRetry(void)
    {
	this->retire(Work::RESIDENT_ADMISSION_RETRY);
    }

    bool residentAdmissionRetryPending(void) const
    {
	return this->pending(Work::RESIDENT_ADMISSION_RETRY);
    }

    void requestExactVisibilityReallocation(void)
    {
	this->request(Work::EXACT_VISIBILITY_REALLOCATION);
    }

    void retireExactVisibilityReallocation(void)
    {
	this->retire(Work::EXACT_VISIBILITY_REALLOCATION);
    }

    bool exactVisibilityReallocationPending(void) const
    {
	return this->pending(Work::EXACT_VISIBILITY_REALLOCATION);
    }

    /* The exact visibility census and the exact framebuffer classification
     * are separate prerequisites.  Reallocation must also yield to structural
     * repair and existing capacity/presentation owners. */
    bool exactVisibilityReallocationReady(bool submissionPending,
	bool exactPresentationPending, bool structuralRepairPending,
	bool capacityTransactionPending, bool presentationBarrierPending) const
    {
	return this->exactVisibilityReallocationPending() &&
	    !submissionPending && !exactPresentationPending &&
	    !structuralRepairPending && !capacityTransactionPending &&
	    !presentationBarrierPending;
    }

private:
    uint8_t workValue = 0;
};

/* Scene-wide permission to try one otherwise-unaffordable discrete
 * population.  Granting and consuming the permit are explicit so multiple
 * time-sliced source actions cannot each infer their own allowance. */
class BObolLodDiscreteTrialPermit {
public:
    enum class State : uint8_t {
	UNAVAILABLE = 0,
	AVAILABLE
    };

    void grant(void) { this->stateValue = State::AVAILABLE; }
    void revoke(void) { this->stateValue = State::UNAVAILABLE; }
    void consume(void) { this->revoke(); }
    void reset(void) { this->revoke(); }
    void setAvailable(bool available)
    {
	this->stateValue = available ? State::AVAILABLE : State::UNAVAILABLE;
    }
    bool available(void) const
    {
	return this->stateValue == State::AVAILABLE;
    }
    State state(void) const { return this->stateValue; }

private:
    State stateValue = State::UNAVAILABLE;
};

/*
 * Allocation-free refinement map from the remaining concrete controller
 * latches to the finite work ledger in ObolProgressivePipeline.tla.  This is
 * deliberately a projection, not another scheduler: effect ownership moves
 * behind the reducer in deletion-producing slices, while readiness and trace
 * diagnostics use this single mapping immediately.
 */
class BObolLodControlRefinement {
public:
    /* Concrete production facts are kept distinct even when they refine to
     * the same abstract work class.  The stable diagnostic mask is the
     * executable bridge to ObolControlRefinement.tla and makes a stuck alias
     * identifiable without reconstructing private controller state. */
    enum class Fact : uint32_t {
	INTERACTION = 1u << 0,
	INVENTORY = 1u << 1,
	AVAILABILITY = 1u << 2,
	RESULT = 1u << 3,
	PUBLICATION = 1u << 4,
	SUBMISSION = 1u << 5,
	SUBMISSION_RESCAN = 1u << 6,
	SUBMISSION_DELTA = 1u << 7,
	QUALITY_PROBE = 1u << 8,
	RETAINED_ALLOCATION = 1u << 9,
	RETAINED_ALLOCATION_TRANSACTION = 1u << 10,
	IMPORTANCE_CENSUS = 1u << 11,
	RESIDENT_ADMISSION_RETRY = 1u << 12,
	CAPACITY_ALLOCATION = 1u << 13,
	RESIDENT_GROWTH = 1u << 14,
	POINT_TRIANGLE_RECOVERY = 1u << 15,
	STRUCTURAL_FRONTIER = 1u << 16,
	PRESENTATION_REPLAY = 1u << 17,
	PRESENTATION_BARRIER = 1u << 18,
	CAPACITY_FRAME = 1u << 19,
	POINT_ADMISSION_FRAME = 1u << 20,
	POINT_CALIBRATION = 1u << 21,
	CAPACITY_CALIBRATION = 1u << 22,
	HEADROOM_PROBE = 1u << 23,
	HANDOFF = 1u << 24,
	COMPACTION = 1u << 25,
	CACHE_WRITE = 1u << 26,
	DEMAND_REFRESH = 1u << 27,
	EXACT_PRESENTATION = 1u << 28
    };

    enum class Work : uint32_t {
	INTERACTION = 1u << 0,
	INVENTORY = 1u << 1,
	AVAILABILITY = 1u << 2,
	PUBLICATION = 1u << 3,
	PLANNING = 1u << 4,
	PRESENTATION = 1u << 5,
	HANDOFF = 1u << 6,
	COMPACTION = 1u << 7,
	CACHE_WRITE = 1u << 8
    };

    enum class Owner : uint8_t {
	NONE = 0,
	INTERACTION,
	INVENTORY,
	AVAILABILITY,
	PUBLICATION,
	PLANNING,
	PRESENTATION,
	HANDOFF,
	COMPACTION,
	CACHE_WRITE
    };

    enum class Violation : uint32_t {
	OWNERLESS_WORK = 1u << 0,
	TERMINAL_WITH_WORK = 1u << 1,
	INVALID_READINESS = 1u << 2,
	INVALID_OWNER = 1u << 3,
	UNWITNESSED_PRESENTATION = 1u << 4,
	UNWITNESSED_CONSTRAINT = 1u << 5,
	UNWITNESSED_PLANNING = 1u << 6,
	NONTERMINAL_WITHOUT_PROGRESS = 1u << 7
    };

    enum class PresentationWitness : uint32_t {
	RENDER = 1u << 0,
	CONTROLLER_PUMP = 1u << 1,
	TIMER = 1u << 2,
	INDEPENDENT_PRODUCER = 1u << 3
    };

    struct PresentationProgress {
	bool renderPending = false;
	bool controllerPumpPending = false;
	bool finiteTimerPending = false;
	bool independentProducerPending = false;

	uint32_t witnessMask(void) const
	{
	    uint32_t mask = 0;
	    if (this->renderPending)
		mask |= bit(PresentationWitness::RENDER);
	    if (this->controllerPumpPending)
		mask |= bit(PresentationWitness::CONTROLLER_PUMP);
	    if (this->finiteTimerPending)
		mask |= bit(PresentationWitness::TIMER);
	    if (this->independentProducerPending)
		mask |= bit(PresentationWitness::INDEPENDENT_PRODUCER);
	    return mask;
	}

	bool witnessed(void) const
	{
	    return this->witnessMask() != 0;
	}
    };

    struct Inputs {
	bool interaction = false;
	bool inventory = false;
	bool availability = false;
	bool result = false;
	bool publication = false;
	bool submission = false;
	bool demandRefresh = false;
	bool submissionRescan = false;
	bool submissionDelta = false;
	bool qualityProbe = false;
	bool retainedAllocation = false;
	bool retainedAllocationTransaction = false;
	bool importanceCensus = false;
	bool residentAdmissionRetry = false;
	/* A bounded capacity candidate is planning work until its occurrence
	 * allocation has been committed.  It must not masquerade as a frame
	 * obligation merely because the same certificate will later require a
	 * presentation sample. */
	bool capacityAllocation = false;
	bool residentGrowth = false;
	bool pointTriangleRecovery = false;
	bool structuralFrontier = false;
	/* An interrupted replay is an exclusive transaction.  Exact-frame debt
	 * is not: a newer capacity allocation may supersede its target before the
	 * eventual exact presentation.  Keep the facts distinct so owner
	 * precedence cannot turn downstream frame debt into an allocation block. */
	bool presentationReplay = false;
	bool exactPresentation = false;
	bool presentationBarrier = false;
	bool capacityFrame = false;
	bool pointAdmissionFrame = false;
	bool pointCalibration = false;
	bool capacityCalibration = false;
	bool headroomProbe = false;
	bool handoff = false;
	bool compaction = false;
	bool cacheWrite = false;
    };

    struct Snapshot {
	uint32_t obligations = 0;
	Owner owner = Owner::NONE;

	bool has(Work obligation) const
	{
	    return (this->obligations & bit(obligation)) != 0;
	}

	bool calibrationPending(void) const
	{
	    return this->has(Work::PRESENTATION);
	}

	bool controlPending(void) const
	{
	    return this->has(Work::PLANNING) || this->has(Work::HANDOFF);
	}

	bool foregroundPending(void) const
	{
	    const uint32_t background =
		bit(Work::COMPACTION) | bit(Work::CACHE_WRITE);
	    return (this->obligations & ~background) != 0;
	}

	bool refinementFramePending(void) const
	{
	    return this->has(Work::PRESENTATION);
	}
    };

    static constexpr uint32_t bit(Work work)
    {
	return static_cast<uint32_t>(work);
    }

    static constexpr uint32_t bit(Fact fact)
    {
	return static_cast<uint32_t>(fact);
    }

    static uint32_t factMask(const Inputs &inputs)
    {
	uint32_t mask = 0;
	const auto add = [&mask](bool active, Fact fact) {
	    if (active)
		mask |= bit(fact);
	};
	add(inputs.interaction, Fact::INTERACTION);
	add(inputs.inventory, Fact::INVENTORY);
	add(inputs.availability, Fact::AVAILABILITY);
	add(inputs.result, Fact::RESULT);
	add(inputs.publication, Fact::PUBLICATION);
	add(inputs.submission, Fact::SUBMISSION);
	add(inputs.demandRefresh, Fact::DEMAND_REFRESH);
	add(inputs.submissionRescan, Fact::SUBMISSION_RESCAN);
	add(inputs.submissionDelta, Fact::SUBMISSION_DELTA);
	add(inputs.qualityProbe, Fact::QUALITY_PROBE);
	add(inputs.retainedAllocation, Fact::RETAINED_ALLOCATION);
	add(inputs.retainedAllocationTransaction,
	    Fact::RETAINED_ALLOCATION_TRANSACTION);
	add(inputs.importanceCensus, Fact::IMPORTANCE_CENSUS);
	add(inputs.residentAdmissionRetry, Fact::RESIDENT_ADMISSION_RETRY);
	add(inputs.capacityAllocation, Fact::CAPACITY_ALLOCATION);
	add(inputs.residentGrowth, Fact::RESIDENT_GROWTH);
	add(inputs.pointTriangleRecovery, Fact::POINT_TRIANGLE_RECOVERY);
	add(inputs.structuralFrontier, Fact::STRUCTURAL_FRONTIER);
	add(inputs.presentationReplay, Fact::PRESENTATION_REPLAY);
	add(inputs.exactPresentation, Fact::EXACT_PRESENTATION);
	add(inputs.presentationBarrier, Fact::PRESENTATION_BARRIER);
	add(inputs.capacityFrame, Fact::CAPACITY_FRAME);
	add(inputs.pointAdmissionFrame, Fact::POINT_ADMISSION_FRAME);
	add(inputs.pointCalibration, Fact::POINT_CALIBRATION);
	add(inputs.capacityCalibration, Fact::CAPACITY_CALIBRATION);
	add(inputs.headroomProbe, Fact::HEADROOM_PROBE);
	add(inputs.handoff, Fact::HANDOFF);
	add(inputs.compaction, Fact::COMPACTION);
	add(inputs.cacheWrite, Fact::CACHE_WRITE);
	return mask;
    }

    static constexpr uint32_t bit(Violation violation)
    {
	return static_cast<uint32_t>(violation);
    }

    static constexpr uint32_t bit(PresentationWitness witness)
    {
	return static_cast<uint32_t>(witness);
    }

    static constexpr uint32_t ownerObligation(Owner owner)
    {
	return owner == Owner::INTERACTION ? bit(Work::INTERACTION) :
	    owner == Owner::INVENTORY ? bit(Work::INVENTORY) :
	    owner == Owner::AVAILABILITY ? bit(Work::AVAILABILITY) :
	    owner == Owner::PUBLICATION ? bit(Work::PUBLICATION) :
	    owner == Owner::PLANNING ? bit(Work::PLANNING) :
	    owner == Owner::PRESENTATION ? bit(Work::PRESENTATION) :
	    owner == Owner::HANDOFF ? bit(Work::HANDOFF) :
	    owner == Owner::COMPACTION ? bit(Work::COMPACTION) :
	    owner == Owner::CACHE_WRITE ? bit(Work::CACHE_WRITE) : 0;
    }

    static Snapshot evaluate(const Inputs &inputs)
    {
	Snapshot snapshot;
	if (inputs.interaction)
	    snapshot.obligations |= bit(Work::INTERACTION);
	if (inputs.inventory)
	    snapshot.obligations |= bit(Work::INVENTORY);
	if (inputs.availability)
	    snapshot.obligations |= bit(Work::AVAILABILITY);
	if (inputs.result || inputs.publication)
	    snapshot.obligations |= bit(Work::PUBLICATION);
	if (inputs.submission || inputs.demandRefresh ||
	    inputs.submissionRescan ||
	    inputs.submissionDelta || inputs.qualityProbe ||
	    inputs.retainedAllocation ||
	    inputs.retainedAllocationTransaction ||
	    inputs.importanceCensus || inputs.residentAdmissionRetry ||
	    inputs.capacityAllocation || inputs.residentGrowth ||
	    inputs.pointTriangleRecovery || inputs.structuralFrontier)
	    snapshot.obligations |= bit(Work::PLANNING);
	if (inputs.presentationReplay || inputs.exactPresentation ||
	    inputs.presentationBarrier ||
	    inputs.capacityFrame || inputs.pointAdmissionFrame ||
	    inputs.pointCalibration || inputs.capacityCalibration ||
	    inputs.headroomProbe)
	    snapshot.obligations |= bit(Work::PRESENTATION);
	if (inputs.handoff)
	    snapshot.obligations |= bit(Work::HANDOFF);
	if (inputs.compaction)
	    snapshot.obligations |= bit(Work::COMPACTION);
	if (inputs.cacheWrite)
	    snapshot.obligations |= bit(Work::CACHE_WRITE);

	/* A closed interrupted presentation transaction blocks every owner-thread
	 * scene mutation until its exact target commits or returns typed capacity
	 * evidence.  The remaining order mirrors the current bounded pump and is
	 * the explicit seam along which those effects will be migrated. */
	if (inputs.presentationReplay) {
	    snapshot.owner = Owner::PRESENTATION;
	} else if (inputs.interaction) {
	    snapshot.owner = Owner::INTERACTION;
	} else if (inputs.result || inputs.publication) {
	    snapshot.owner = Owner::PUBLICATION;
	} else if (inputs.inventory) {
	    snapshot.owner = Owner::INVENTORY;
	} else if (inputs.availability) {
	    snapshot.owner = Owner::AVAILABILITY;
	} else if (inputs.capacityAllocation) {
	    /* The candidate's exact frame is a successor of its occurrence
	     * allocation.  An older exact-frame or presentation-barrier latch may
	     * remain visible while the candidate is installed, but it cannot be
	     * the enabled owner until this planning edge completes. */
	    snapshot.owner = Owner::PLANNING;
	} else if (snapshot.has(Work::PRESENTATION)) {
	    snapshot.owner = Owner::PRESENTATION;
	} else if (snapshot.has(Work::PLANNING)) {
	    snapshot.owner = Owner::PLANNING;
	} else if (inputs.handoff) {
	    snapshot.owner = Owner::HANDOFF;
	} else if (inputs.compaction) {
	    snapshot.owner = Owner::COMPACTION;
	} else if (inputs.cacheWrite) {
	    snapshot.owner = Owner::CACHE_WRITE;
	}
	return snapshot;
    }

    static uint32_t validateProducers(const Inputs &inputs)
    {
	/* A source-delta fact is a mode of the bounded submission cursor, not
	 * durable semantic debt.  Once the cursor retires the selective scope must
	 * retire atomically; otherwise it can claim planning ownership while
	 * blocking the broader demand pass which should follow it. */
	uint32_t violations = inputs.submissionDelta && !inputs.submission ?
	    bit(Violation::UNWITNESSED_PLANNING) : 0;
	/* The retained importance census is consumed by either a complete
	 * current-view demand pass or the bounded inventory/coverage pass which
	 * proves that every current occurrence was considered.  Submission by
	 * itself is not a witness: it may be a selective source-delta pass which
	 * cannot retire scene-wide census debt. */
	const bool importanceProducer = inputs.demandRefresh ||
	    (inputs.submission && inputs.inventory);
	if (inputs.importanceCensus && !importanceProducer)
	    violations |= bit(Violation::UNWITNESSED_PLANNING);
	return violations;
    }

    static uint32_t validate(const Snapshot &snapshot, bool terminal,
	bool viewReady, bool terminalError,
	bool presentationProgressWitness, bool constrainedOutcome,
	bool constraintEvidence)
    {
	uint32_t violations = 0;
	if (snapshot.obligations != 0 && snapshot.owner == Owner::NONE)
	    violations |= bit(Violation::OWNERLESS_WORK);
	const uint32_t ownedObligation = ownerObligation(snapshot.owner);
	if (snapshot.owner != Owner::NONE &&
	    !(snapshot.obligations & ownedObligation))
	    violations |= bit(Violation::INVALID_OWNER);
	if (terminal && snapshot.foregroundPending())
	    violations |= bit(Violation::TERMINAL_WITH_WORK);
	if (viewReady && (!terminal || terminalError))
	    violations |= bit(Violation::INVALID_READINESS);
	/* PRESENTATION is an active transition, not a descriptive phase label.
	 * It must own a requested frame, an independent producer which can publish
	 * that frame, a finite publication timer, or the controller-scoped pump
	 * which converts its standing obligation into a render request.  A generic
	 * shared host pump is not sufficient: an unrelated provider wakeup must not
	 * conceal a stalled presentation reducer. */
	if (snapshot.owner == Owner::PRESENTATION &&
	    !presentationProgressWitness)
	    violations |= bit(Violation::UNWITNESSED_PRESENTATION);
	/* CONSTRAINED is a terminal proof, not a synonym for "the image is
	 * coarse."  It must name completed capacity, presentation, or resource
	 * evidence which prevents further foreground work in this epoch. */
	if (constrainedOutcome && !constraintEvidence)
	    violations |= bit(Violation::UNWITNESSED_CONSTRAINT);
	return violations;
    }

    /* The abstract convergence machine gives every nonterminal state a
     * finite successor.  Service workers are independent producers and may
     * legitimately own that successor outside this controller ledger; absent
     * either kind of witness, a nonterminal LoD state is a liveness defect. */
    static uint32_t validateLiveness(const Snapshot &snapshot, bool terminal,
	bool hasLodState, bool externalProgressWitness)
    {
	return hasLodState && !terminal && !snapshot.foregroundPending() &&
	    !externalProgressWitness ?
	    bit(Violation::NONTERMINAL_WITHOUT_PROGRESS) : 0;
    }
};

#endif /* LIBBOBOL_LOD_CONTROL_PRIVATE_H */
