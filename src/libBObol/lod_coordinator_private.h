/*          L O D _ C O O R D I N A T O R _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_COORDINATOR_PRIVATE_H
#define LIBBOBOL_LOD_COORDINATOR_PRIVATE_H

#include "common.h"
#include "BObol/BLodIdentifiers.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

/*
 * Use the same strong epoch domains as requests and results.  Keeping a
 * second coordinator-only implementation previously allowed the state
 * machine and service boundary to drift while appearing type safe.
 */
using BObolLodViewEpoch = BObolViewEpoch;
using BObolLodPolicyEpoch = BObolPolicyEpoch;
using BObolLodInventoryEpoch = BObolInventoryEpoch;
using BObolLodSourceRoutingId = BObolSourceRoutingId;

/*
 * Completed-frame resource-pressure edge policy.  The controller supplies
 * measurements and executes the resulting compaction request, but it may not
 * independently reinterpret the revision/one-shot contract.  Keeping this
 * allocation-free value beside the phase machine makes persistent pressure a
 * provable memory-limited terminal state instead of a collection of latches.
 */
class BObolLodResourcePolicy {
public:
    struct Decision {
	bool changed = false;
	bool queueRecovery = false;
	bool pressureCleared = false;
	uint64_t revision = 0;
    };

    Decision observe(bool cpuPressure, bool gpuPressure,
	bool recoveryEnabled)
    {
	Decision decision;
	decision.revision = this->pressureRevision;
	if (cpuPressure == this->cpuPressureValue &&
	    gpuPressure == this->gpuPressureValue) {
	    /* The pressure sample may predate service/automatic submission.
	     * Enabling recovery later must consume that still-current edge once,
	     * without requiring pressure to fall and rise again. */
	    if ((cpuPressure || gpuPressure) && recoveryEnabled &&
		this->handledRevision < this->pressureRevision &&
		!this->recoveryPendingValue) {
		this->recoveryPendingValue = true;
		decision.queueRecovery = true;
	    }
	    return decision;
	}

	this->cpuPressureValue = cpuPressure;
	this->gpuPressureValue = gpuPressure;
	this->pressureRevision++;
	if (!this->pressureRevision)
	    this->pressureRevision++;
	decision.changed = true;
	decision.revision = this->pressureRevision;
	if (cpuPressure || gpuPressure) {
	    this->recoveryPendingValue = recoveryEnabled;
	    decision.queueRecovery = this->recoveryPendingValue;
	} else {
	    this->recoveryPendingValue = false;
	    this->handledRevision = this->pressureRevision;
	    decision.pressureCleared = true;
	}
	return decision;
    }

    void markRecoveryHandled(void)
    {
	this->handledRevision = this->pressureRevision;
	this->recoveryPendingValue = false;
    }

    void resetForServiceChange(void)
    {
	this->cpuPressureValue = false;
	this->gpuPressureValue = false;
	this->recoveryPendingValue = false;
	this->handledRevision = this->pressureRevision;
    }

    bool cpuPressure(void) const { return this->cpuPressureValue; }
    bool gpuPressure(void) const { return this->gpuPressureValue; }
    bool anyPressure(void) const
    {
	return this->cpuPressureValue || this->gpuPressureValue;
    }
    bool recoveryPending(void) const
    {
	return this->recoveryPendingValue &&
	    this->handledRevision < this->pressureRevision;
    }
    uint64_t revision(void) const { return this->pressureRevision; }
    uint64_t recoveryHandledRevision(void) const
    {
	return this->handledRevision;
    }

private:
    bool cpuPressureValue = false;
    bool gpuPressureValue = false;
    bool recoveryPendingValue = false;
    uint64_t pressureRevision = 0;
    uint64_t handledRevision = 0;
};

/*
 * One-shot proof for a late calibrated-headroom admission retry.  Normal
 * probe counts remain a bounded short-horizon mechanism; this witness covers
 * the ordering where reusable timing becomes available only after those
 * probes ended.  A retry is unique for (view, policy, current budget), so an
 * impossible discrete next population cannot repaint forever.
 */
class BObolLodHeadroomPolicy {
public:
    bool armRetry(BObolLodViewEpoch viewEpoch,
	BObolLodPolicyEpoch policyEpoch, size_t currentBudget)
    {
	if ((this->pendingValue &&
		this->pendingViewRevision == viewEpoch.value() &&
		this->pendingPolicyRevision == policyEpoch.value() &&
		this->pendingBudget == currentBudget) ||
	    (this->viewRevision == viewEpoch.value() &&
		this->policyRevision == policyEpoch.value() &&
		this->budget == currentBudget))
	    return false;
	this->pendingValue = true;
	this->pendingViewRevision = viewEpoch.value();
	this->pendingPolicyRevision = policyEpoch.value();
	this->pendingBudget = currentBudget;
	return true;
    }

    bool consumeRetry(BObolLodViewEpoch viewEpoch,
	BObolLodPolicyEpoch policyEpoch, size_t currentBudget,
	long double demonstratedBudget, uint64_t elapsedNanoseconds,
	uint64_t targetNanoseconds, bool reusableSample)
    {
	if (!this->pendingValue)
	    return false;
	const bool matchingWitness =
	    this->pendingViewRevision == viewEpoch.value() &&
	    this->pendingPolicyRevision == policyEpoch.value() &&
	    this->pendingBudget == currentBudget;
	this->pendingValue = false;
	this->pendingViewRevision = 0;
	this->pendingPolicyRevision = 0;
	this->pendingBudget = SIZE_MAX;
	if (!matchingWitness)
	    return false;
	/* Consuming the explicitly requested replay retires this exact witness
	 * even if the evidence is insufficient.  Otherwise an unchanged discrete
	 * population can request calibration frames forever. */
	this->viewRevision = viewEpoch.value();
	this->policyRevision = policyEpoch.value();
	this->budget = currentBudget;
	if (!reusableSample || !targetNanoseconds ||
	    !(demonstratedBudget * 20.0L >
		static_cast<long double>(currentBudget) * 21.0L) ||
	    static_cast<long double>(elapsedNanoseconds) * 5.0L >=
		static_cast<long double>(targetNanoseconds) * 6.0L)
	    return false;
	return true;
    }

    void cancelRetry(void)
    {
	this->pendingValue = false;
	this->pendingViewRevision = 0;
	this->pendingPolicyRevision = 0;
	this->pendingBudget = SIZE_MAX;
    }

    bool retryPending(void) const { return this->pendingValue; }
    bool pendingMatches(BObolLodViewEpoch viewEpoch,
	BObolLodPolicyEpoch policyEpoch, size_t currentBudget) const
    {
	return this->pendingValue &&
	    this->pendingViewRevision == viewEpoch.value() &&
	    this->pendingPolicyRevision == policyEpoch.value() &&
	    this->pendingBudget == currentBudget;
    }

private:
    uint64_t viewRevision = 0;
    uint64_t policyRevision = 0;
    size_t budget = SIZE_MAX;
    bool pendingValue = false;
    uint64_t pendingViewRevision = 0;
    uint64_t pendingPolicyRevision = 0;
    size_t pendingBudget = SIZE_MAX;
};

/*
 * Current-view minimum-mesh coverage proof.  Large compact sources are
 * visited in bounded owner-thread windows, but their counters describe one
 * logical pass and must therefore survive between windows.  Completing a
 * pass atomically publishes its visible denominator and retires the counters;
 * camera/inventory invalidation explicitly clears that proof, whereas a
 * quality-only policy pass may preserve it.
 */
class BObolLodCoveragePolicy {
public:
    struct Completion {
	bool completed = false;
	bool bounded = false;
	bool missing = false;
	size_t visibleCount = 0;
	size_t coveredCount = 0;
    };

    void activate(bool invalidateCoverage)
    {
	this->activeValue = true;
	if (invalidateCoverage)
	    this->coverageCompleteValue = false;
    }

    /* The semantic view policy owns whether minimum-mesh coverage is a
     * prerequisite at all.  Keep that fact beside the coverage proof instead
     * of adding a second controller latch: LoD-disabled and detached-service
     * views are vacuously covered, while a stale compaction cursor may never
     * overlap an invalidated required inventory. */
    void setRequired(bool required)
    {
	this->requiredValue = required;
    }

    bool required(void) const { return this->requiredValue; }
    bool effectiveActive(void) const
    {
	return this->requiredValue && this->activeValue;
    }
    bool effectiveComplete(void) const
    {
	return !this->requiredValue || this->coverageCompleteValue;
    }
    bool compactionAllowed(void) const
    {
	return this->effectiveComplete();
    }

    void deactivate(void)
    {
	this->activeValue = false;
    }

    void markBoundedSource(void)
    {
	if (this->activeValue)
	    this->sawBoundedSourceValue = true;
    }

    void observe(size_t visibleCount, size_t coveredCount)
    {
	if (!this->activeValue)
	    return;
	this->visibleCountValue = saturatingAdd(
	    this->visibleCountValue, visibleCount);
	this->coveredCountValue = saturatingAdd(
	    this->coveredCountValue, coveredCount);
    }

    void clearPassCounters(void)
    {
	this->sawBoundedSourceValue = false;
	this->visibleCountValue = 0;
	this->coveredCountValue = 0;
    }

    Completion completeIfReady(bool completedPass, bool rescanPending)
    {
	Completion completion;
	if (!completedPass || !this->activeValue || rescanPending)
	    return completion;
	completion.completed = true;
	completion.bounded = this->sawBoundedSourceValue;
	completion.visibleCount = this->visibleCountValue;
	completion.coveredCount = this->coveredCountValue;
	completion.missing = completion.coveredCount < completion.visibleCount;
	this->coverageCompleteValue = !completion.missing;
	this->completeVisibleCountValue = completion.visibleCount;
	this->completeVisibleCountValidValue = true;
	this->activeValue = false;
	this->clearPassCounters();
	return completion;
    }

    void clearCompleteVisibleCount(void)
    {
	this->completeVisibleCountValue = 0;
	this->completeVisibleCountValidValue = false;
    }

    void reset(void)
    {
	this->activeValue = false;
	this->coverageCompleteValue = false;
	this->clearPassCounters();
	this->clearCompleteVisibleCount();
    }

    bool active(void) const { return this->activeValue; }
    bool sawBoundedSource(void) const
    {
	return this->sawBoundedSourceValue;
    }
    bool coverageComplete(void) const
    {
	return this->coverageCompleteValue;
    }
    size_t visibleCount(void) const { return this->visibleCountValue; }
    size_t coveredCount(void) const { return this->coveredCountValue; }
    bool hasCompleteVisibleCount(void) const
    {
	return this->completeVisibleCountValidValue;
    }
    size_t completeVisibleCount(void) const
    {
	return this->completeVisibleCountValue;
    }

private:
    static size_t saturatingAdd(size_t left, size_t right)
    {
	return right > SIZE_MAX - left ? SIZE_MAX : left + right;
    }

    bool requiredValue = false;
    bool activeValue = true;
    bool sawBoundedSourceValue = false;
    bool coverageCompleteValue = false;
    bool completeVisibleCountValidValue = false;
    size_t visibleCountValue = 0;
    size_t coveredCountValue = 0;
    size_t completeVisibleCountValue = 0;
};

/*
 * User-visible convergence proof.  The controller samples detailed source,
 * worker, renderer, and memory observations; this policy alone decides
 * readiness, terminal limitations, HUD phase, and the monotonic progress
 * fraction for a fixed view/policy epoch.  It contains no scene pointers and
 * performs no allocation, so diagnostics cannot accidentally become another
 * scheduler or hierarchy traversal.
 */
class BObolLodConvergencePolicy {
public:
    enum class Phase : uint8_t {
	IDLE = 0,
	DISCOVERING,
	INTERACTIVE,
	REFINING,
	CALIBRATING,
	BACKGROUND,
	ERROR
    };

    struct Inputs {
	BObolLodViewEpoch viewEpoch;
	BObolLodPolicyEpoch policyEpoch;
	size_t expectedLeafCount = 0;
	size_t availableLeafCount = 0;
	size_t visibleTargetCount = 0;
	size_t activePayloadCount = 0;
	size_t satisfiedPayloadCount = 0;
	size_t memoryLimitedPayloadCount = 0;
	size_t pendingTasks = 0;
	size_t inFlight = 0;
	size_t queuedCacheWrites = 0;
	size_t stableResidentMeshBytes = 0;
	size_t residentMeshLimitBytes = SIZE_MAX;
	size_t gpuTrackedBufferBytes = 0;
	unsigned int failedSourceCount = 0;
	bool structuralDiscovery = false;
	bool submissionPending = false;
	bool resultPending = false;
	bool publicationPending = false;
	bool calibrationPending = false;
	bool interactive = false;
	bool compactionPending = false;
	bool progressiveWorkPending = false;
	bool gpuMemoryPressure = false;
	bool stableBudgetLimited = false;
    };

    struct Decision {
	Phase phase = Phase::IDLE;
	float fraction = 1.0f;
	bool viewReady = false;
	bool backgroundPending = false;
	bool performanceLimited = false;
	bool memoryLimited = false;
	bool visualPending = false;
	bool hasLodState = false;
    };

    Decision evaluate(const Inputs &inputs)
    {
	Decision decision;
	const bool structuralPending =
	    inputs.expectedLeafCount > inputs.availableLeafCount;
	decision.visualPending =
	    structuralPending || inputs.submissionPending ||
	    inputs.resultPending || inputs.publicationPending ||
	    inputs.pendingTasks > 0 || inputs.inFlight > 0 ||
	    inputs.calibrationPending;
	decision.viewReady =
	    !decision.visualPending && !inputs.interactive;
	decision.backgroundPending =
	    inputs.queuedCacheWrites > 0 || inputs.compactionPending ||
	    (decision.viewReady && inputs.progressiveWorkPending);
	decision.memoryLimited =
	    inputs.memoryLimitedPayloadCount > 0 ||
	    inputs.gpuMemoryPressure ||
	    (inputs.residentMeshLimitBytes != SIZE_MAX &&
	     inputs.stableResidentMeshBytes > inputs.residentMeshLimitBytes);
	decision.performanceLimited = decision.viewReady &&
	    (inputs.stableBudgetLimited || decision.memoryLimited ||
	     (inputs.visibleTargetCount > 0 &&
	      (inputs.activePayloadCount < inputs.visibleTargetCount ||
	       inputs.satisfiedPayloadCount < inputs.visibleTargetCount)));
	decision.hasLodState =
	    inputs.expectedLeafCount > 0 || inputs.availableLeafCount > 0 ||
	    inputs.visibleTargetCount > 0 || inputs.activePayloadCount > 0 ||
	    inputs.gpuTrackedBufferBytes > 0 || decision.visualPending ||
	    decision.backgroundPending;
	if (!decision.hasLodState)
	    return decision;

	if (structuralPending || inputs.structuralDiscovery) {
	    decision.phase = Phase::DISCOVERING;
	    const long double coverage =
		structuralPending && inputs.expectedLeafCount > 0 ?
		    static_cast<long double>(inputs.availableLeafCount) /
			static_cast<long double>(inputs.expectedLeafCount) : 0.0L;
	    decision.fraction = static_cast<float>(
		std::min<long double>(0.40L, 0.40L * coverage));
	} else if (inputs.interactive) {
	    decision.phase = Phase::INTERACTIVE;
	    if (inputs.visibleTargetCount == 0) {
		decision.fraction = 0.40f;
	    } else {
		const long double quality =
		    static_cast<long double>(std::min(
			inputs.satisfiedPayloadCount,
			inputs.visibleTargetCount)) /
		    static_cast<long double>(inputs.visibleTargetCount);
		decision.fraction = static_cast<float>(
		    0.40L + std::min<long double>(0.45L, 0.45L * quality));
	    }
	} else if (inputs.calibrationPending) {
	    decision.phase = Phase::CALIBRATING;
	    decision.fraction = 0.95f;
	} else if (decision.visualPending) {
	    decision.phase = Phase::REFINING;
	    size_t target = std::max(
		inputs.visibleTargetCount, inputs.activePayloadCount);
	    if (!target) {
		target = saturatingAdd(inputs.pendingTasks, inputs.inFlight);
		target = std::max<size_t>(1, target);
	    }
	    const long double quality =
		static_cast<long double>(std::min(
		    inputs.satisfiedPayloadCount, target)) /
		static_cast<long double>(target);
	    decision.fraction = static_cast<float>(
		0.40L + std::min<long double>(0.55L, 0.55L * quality));
	} else if (decision.backgroundPending) {
	    decision.phase = Phase::BACKGROUND;
	    decision.fraction = 1.0f;
	} else {
	    decision.phase = Phase::IDLE;
	    decision.fraction = 1.0f;
	}

	if (inputs.failedSourceCount > 0)
	    decision.phase = Phase::ERROR;
	decision.fraction = std::max(0.0f,
	    std::min(1.0f, decision.fraction));
	if (this->fractionViewEpoch != inputs.viewEpoch ||
	    this->fractionPolicyEpoch != inputs.policyEpoch) {
	    this->fractionViewEpoch = inputs.viewEpoch;
	    this->fractionPolicyEpoch = inputs.policyEpoch;
	    this->fractionFloorValue = 0.0f;
	}
	if (decision.phase != Phase::ERROR) {
	    decision.fraction = std::max(
		decision.fraction, this->fractionFloorValue);
	    this->fractionFloorValue = decision.fraction;
	}
	return decision;
    }

    void resetFraction(void)
    {
	this->fractionViewEpoch.reset();
	this->fractionPolicyEpoch.reset();
	this->fractionFloorValue = 0.0f;
    }

    float fractionFloor(void) const { return this->fractionFloorValue; }

private:
    static size_t saturatingAdd(size_t left, size_t right)
    {
	return right > SIZE_MAX - left ? SIZE_MAX : left + right;
    }

    BObolLodViewEpoch fractionViewEpoch;
    BObolLodPolicyEpoch fractionPolicyEpoch;
    float fractionFloorValue = 0.0f;
};

/*
 * Aggregate scene render-cost admission.  Per-occurrence projected error is
 * demand; this value owns how much of that demand one measured renderer may
 * expose in a frame.  It retains the current allowance, bounded pass
 * remainders, retained-cut recovery admission, and the one-shot overload
 * witness.  No occurrence data or renderer objects enter this policy.
 */
class BObolLodBudgetPolicy {
public:
    struct Inputs {
	size_t activeCost = 0;
	float targetFps = 0.0f;
	long double calibratedCostPerSecond = 0.0L;
	uint64_t observedStableNanoseconds = 0;
	uint64_t lastRenderNanoseconds = 0;
	uint64_t smoothedRenderNanoseconds = 0;
	bool interactive = false;
	bool scaleQualityProbe = false;
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
	uint64_t targetNanoseconds = 0;
	uint64_t observedNanoseconds = 0;
	long double calibratedBudget = 0.0L;
	bool interactive = false;
	bool stablePresentationHandoff = false;
	bool passAdmittedWork = false;
    };

    struct CalibrationDecision {
	bool probeCandidate = false;
	bool probeEligible = false;
	bool requestFrame = false;
	bool calibrationFrame = false;
    };

    struct CompletedFrameDecision {
	bool requestCalibrationFrame = false;
	bool restartSubmission = false;
    };

    Decision beginPass(const Inputs &inputs)
    {
	Decision decision;
	decision.totalBudget = this->currentBudgetValue;
	decision.refinementBudget = this->refinementRemainingValue;
	decision.retainedAdmission = this->retainedAdmissionValue;
	decision.retainedAdmissionBudget =
	    this->retainedAdmissionRemainingValue;
	if (this->passInitializedValue)
	    return decision;

	const long double targetNanoseconds = inputs.targetFps > 0.0f ?
	    1000000000.0L / static_cast<long double>(inputs.targetFps) : 0.0L;
	const long double calibratedBudget =
	    inputs.targetFps > 0.0f &&
	    inputs.calibratedCostPerSecond > 0.0L ?
		inputs.calibratedCostPerSecond /
		    static_cast<long double>(inputs.targetFps) : 0.0L;
	const long double observedStableNanoseconds =
	    static_cast<long double>(inputs.observedStableNanoseconds);
	const bool severeStableOverload =
	    !inputs.interactive && targetNanoseconds > 0.0L &&
	    observedStableNanoseconds > targetNanoseconds * 2.0L;
	const bool overloadRecovery =
	    !inputs.interactive && !inputs.forceTerminal &&
	    (!this->overloadRecoveryPerformedValue ||
	     this->overloadRecoveryActiveCostValue != inputs.activeCost) &&
	    inputs.activeCost > 0 &&
	    (this->probeCountValue >= 3 ||
	     severeStableOverload) &&
	    targetNanoseconds > 0.0L &&
	    observedStableNanoseconds > targetNanoseconds * 1.20L &&
	    (severeStableOverload ||
	     (calibratedBudget > 0.0L && calibratedBudget <
		static_cast<long double>(inputs.activeCost) * 0.95L));

	size_t costBudget = this->seedBudgetValue;
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
	    if (inputs.activeCost > 0 && costBudget > inputs.activeCost) {
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
		this->currentBudgetValue != SIZE_MAX)
		costBudget = std::max(costBudget, this->currentBudgetValue);
	    if (!inputs.interactive && this->currentBudgetValue != SIZE_MAX &&
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
	/* Preserve only work which an exact completed presentation proved could
	 * meet the stable deadline.  This is deliberately a measured render-cost
	 * floor, not the previous pass allowance: the latter may include admission
	 * which was never submitted or presented. */
	if (inputs.stablePresentationCostFloor > 0 && costBudget != SIZE_MAX)
	    costBudget = std::max(
		costBudget, inputs.stablePresentationCostFloor);

	this->currentBudgetValue = costBudget;
	this->refinementRemainingValue = costBudget == SIZE_MAX ? SIZE_MAX :
	    (costBudget > inputs.activeCost ?
		costBudget - inputs.activeCost : 0);
	this->retainedAdmissionValue =
	    (inputs.interactive || overloadRecovery ||
	     inputs.stablePresentationHandoff) &&
	    costBudget != SIZE_MAX && inputs.activeCost > costBudget;
	if (overloadRecovery) {
	    this->overloadRecoveryPerformedValue = true;
	    this->overloadRecoveryActiveCostValue = inputs.activeCost;
	}
	this->retainedAdmissionRemainingValue =
	    this->retainedAdmissionValue ? costBudget : SIZE_MAX;
	this->passInitializedValue = true;

	decision.initialized = true;
	decision.overloadRecovery = overloadRecovery;
	decision.totalBudget = this->currentBudgetValue;
	decision.refinementBudget = this->refinementRemainingValue;
	decision.retainedAdmission = this->retainedAdmissionValue;
	decision.retainedAdmissionBudget =
	    this->retainedAdmissionRemainingValue;
	return decision;
    }

    CalibrationDecision finishBlockedPass(
	const CalibrationInputs &inputs)
    {
	CalibrationDecision decision;
	const bool calibratedHeadroom =
	    inputs.calibratedBudget >
		static_cast<long double>(this->currentBudgetValue) * 1.05L &&
	    inputs.observedNanoseconds < static_cast<uint64_t>(
		static_cast<long double>(inputs.targetNanoseconds) * 1.20L);
	decision.probeCandidate =
	    !inputs.interactive && !inputs.stablePresentationHandoff &&
	    !inputs.passAdmittedWork && inputs.activeCost > 0 &&
	    this->currentBudgetValue != SIZE_MAX &&
	    inputs.targetNanoseconds > 0 &&
	    (inputs.observedNanoseconds < static_cast<uint64_t>(
		static_cast<long double>(inputs.targetNanoseconds) * 0.90L) ||
	     calibratedHeadroom);
	if (inputs.passAdmittedWork) {
	    this->resetProbeSeries();
	} else if (this->probeActiveCostValue != inputs.activeCost) {
	    this->probeActiveCostValue = inputs.activeCost;
	    this->probeCountValue = 0;
	    this->calibrationFramesRemainingValue = 0;
	}
	static constexpr unsigned int minimumProbes = 3;
	static constexpr unsigned int maximumProbes = 12;
	decision.probeEligible =
	    !inputs.interactive && !inputs.stablePresentationHandoff &&
	    !inputs.passAdmittedWork && inputs.activeCost > 0 &&
	    this->currentBudgetValue != SIZE_MAX &&
	    inputs.targetNanoseconds > 0;
	decision.calibrationFrame = decision.probeEligible &&
	    (this->probeCountValue < minimumProbes ||
	     decision.probeCandidate) &&
	    this->probeCountValue < maximumProbes;
	if (decision.calibrationFrame) {
	    this->probeCountValue++;
	    const unsigned int targetProbeCount =
		decision.probeCandidate ? maximumProbes : minimumProbes;
	    this->calibrationFramesRemainingValue =
		targetProbeCount > this->probeCountValue ?
		    targetProbeCount - this->probeCountValue : 0;
	} else if (!decision.probeCandidate) {
	    this->probeCountValue = 0;
	    this->calibrationFramesRemainingValue = 0;
	}
	this->rescanAfterFrameValue =
	    inputs.passAdmittedWork || decision.calibrationFrame;
	this->stableBudgetLimitedValue = !this->rescanAfterFrameValue;
	decision.requestFrame = this->rescanAfterFrameValue;
	return decision;
    }

    CompletedFrameDecision completeCalibrationFrame(void)
    {
	CompletedFrameDecision decision;
	if (!this->rescanAfterFrameValue)
	    return decision;
	if (this->calibrationFramesRemainingValue > 0) {
	    this->calibrationFramesRemainingValue--;
	    this->probeCountValue++;
	    decision.requestCalibrationFrame = true;
	    return decision;
	}
	this->rescanAfterFrameValue = false;
	this->resetPass();
	decision.restartSubmission = true;
	return decision;
    }

    void requestRescanAfterFrame(bool clearCalibrationFrames = false)
    {
	this->rescanAfterFrameValue = true;
	this->stableBudgetLimitedValue = false;
	if (clearCalibrationFrames)
	    this->calibrationFramesRemainingValue = 0;
    }

    void clearBudgetLimit(void)
    {
	this->stableBudgetLimitedValue = false;
    }

    void resetProbeSeries(void)
    {
	this->probeActiveCostValue = 0;
	this->probeCountValue = 0;
	this->calibrationFramesRemainingValue = 0;
    }

    void setProbeCount(unsigned int count)
    {
	this->probeCountValue = count;
    }

    void resetCalibration(void)
    {
	this->resetProbeSeries();
	this->rescanAfterFrameValue = false;
	this->stableBudgetLimitedValue = false;
    }

    void resetPass(void)
    {
	this->passInitializedValue = false;
	this->refinementRemainingValue = SIZE_MAX;
	this->retainedAdmissionValue = false;
	this->retainedAdmissionRemainingValue = SIZE_MAX;
    }

    void resetOverloadRecovery(void)
    {
	this->overloadRecoveryPerformedValue = false;
	this->overloadRecoveryActiveCostValue = 0;
    }

    void reset(void)
    {
	this->currentBudgetValue = this->seedBudgetValue;
	this->resetPass();
	this->resetOverloadRecovery();
	this->resetCalibration();
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

    void reduceCurrentBudget(size_t budget)
    {
	this->currentBudgetValue = std::min(
	    this->currentBudgetValue, budget);
    }

    void raiseCurrentBudget(size_t budget)
    {
	this->currentBudgetValue = std::max(
	    this->currentBudgetValue, budget);
    }

    size_t seedBudget(void) const { return this->seedBudgetValue; }
    size_t currentBudget(void) const { return this->currentBudgetValue; }
    bool passInitialized(void) const { return this->passInitializedValue; }
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
    bool overloadRecoveryPerformed(void) const
    {
	return this->overloadRecoveryPerformedValue;
    }
    size_t overloadRecoveryActiveCost(void) const
    {
	return this->overloadRecoveryActiveCostValue;
    }
    size_t probeActiveCost(void) const { return this->probeActiveCostValue; }
    unsigned int probeCount(void) const { return this->probeCountValue; }
    unsigned int calibrationFramesRemaining(void) const
    {
	return this->calibrationFramesRemainingValue;
    }
    bool rescanAfterFrame(void) const
    {
	return this->rescanAfterFrameValue;
    }
    bool stableBudgetLimited(void) const
    {
	return this->stableBudgetLimitedValue;
    }

private:
    size_t seedBudgetValue = 50000;
    size_t currentBudgetValue = 50000;
    size_t refinementRemainingValue = SIZE_MAX;
    size_t retainedAdmissionRemainingValue = SIZE_MAX;
    size_t overloadRecoveryActiveCostValue = 0;
    size_t probeActiveCostValue = 0;
    unsigned int probeCountValue = 0;
    unsigned int calibrationFramesRemainingValue = 0;
    bool passInitializedValue = false;
    bool retainedAdmissionValue = false;
    bool overloadRecoveryPerformedValue = false;
    bool rescanAfterFrameValue = false;
    bool stableBudgetLimitedValue = false;
};

/*
 * Camera-local view-demand scheduling.  A zoom interaction may keep loading
 * the pixel-demanded immutable PoP suffix while exposing only one additional
 * coherent population per completed motion frame.  This allocation-free
 * policy owns that bounded-probe lifecycle and the proof that a coverage pass
 * is a scale-demand refresh rather than a cold scene restart.  Renderer
 * ceilings and scene measurements remain controller actions/inputs.
 */
class BObolLodViewDemandPolicy {
public:
    static constexpr uint64_t qualityFrameDurationNanoseconds(void)
    {
	return 100000000ULL;
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
	bool refinementFramePending = false;
	bool publicationFramePending = false;
	bool motionFramePending = false;
    };

    struct QualityProbeDecision {
	bool consumed = false;
	bool begin = false;
	int progressiveCeiling = -1;
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

    void beginGesture(bool newInteractionEpoch)
    {
	if (newInteractionEpoch)
	    this->interactionScaleChangedValue = false;
	this->viewScaleChangingValue = false;
	this->qualityProbeActiveValue = false;
	this->qualityProbePendingValue = false;
	this->qualityProbePresentedValue = false;
	this->qualityBudgetActiveValue = false;
    }

    void beginCameraInteraction(bool newInteractionEpoch,
	bool scaleChanged)
    {
	if (newInteractionEpoch)
	    this->interactionScaleChangedValue = false;
	if (scaleChanged)
	    this->interactionScaleChangedValue = true;
    }

    void endGesture(void)
    {
	this->qualityProbeActiveValue = false;
	this->qualityProbePendingValue = false;
	this->qualityProbePresentedValue = false;
    }

    void noteFramePresented(void)
    {
	if (this->qualityBudgetActiveValue)
	    this->qualityProbePresentedValue = true;
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
	decision.begin = true;
	/* Residency may advance several populations while a render-only ceiling
	 * protects interaction.  Probe from the cut actually being presented,
	 * not from that hidden resident/occurrence maximum; jumping back to the
	 * latter on every wheel event repeatedly recreates the same long frame
	 * and ratchets a previously safe ceiling downward. */
	const int presentedMaximum = inputs.presentationCeiling >= 0 ?
	    std::min(inputs.activeMaximum, inputs.presentationCeiling) :
	    inputs.activeMaximum;
	decision.progressiveCeiling = presentedMaximum >= 15 ?
	    15 : presentedMaximum + 1;
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
	return completedScaleInteraction;
    }

    void refreshForViewRevision(bool interactive)
    {
	this->scaleDemandRefreshActiveValue =
	    this->viewScaleChangingValue ||
	    (interactive && this->interactionScaleChangedValue);
    }

    void refreshForPolicyRevision(bool preserveScaleDemandRefresh,
	bool interactive)
    {
	this->scaleDemandRefreshActiveValue =
	    preserveScaleDemandRefresh ||
	    (interactive && this->interactionScaleChangedValue);
    }

    void clearDemandRefresh(void)
    {
	this->scaleDemandRefreshActiveValue = false;
    }

    void reset(void)
    {
	this->scaleDemandRefreshActiveValue = false;
	this->viewScaleChangingValue = false;
	this->interactionScaleChangedValue = false;
	this->qualityProbeActiveValue = false;
	this->qualityProbePendingValue = false;
	this->qualityProbePresentedValue = false;
	this->qualityBudgetActiveValue = false;
    }

    bool scaleDemandRefreshActive(void) const
    {
	return this->scaleDemandRefreshActiveValue;
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
    bool scaleDemandRefreshActiveValue = false;
    bool viewScaleChangingValue = false;
    bool interactionScaleChangedValue = false;
    bool qualityProbeActiveValue = false;
    bool qualityProbePendingValue = false;
    bool qualityProbePresentedValue = false;
    bool qualityBudgetActiveValue = false;
};

/*
 * Renderer-only presentation continuity across camera interaction.  The
 * retained occurrence cuts and resident PoP arrays remain authoritative; this
 * allocation-free value owns only the reversible global prefix/point-proxy
 * presentation and its motion-to-stable handoff proof.
 *
 * Two snapshots are intentionally distinct:
 *
 *  - the pre-interaction stable snapshot may be restored after pose-only
 *    orthographic motion when the retained occurrence population is unchanged;
 *  - an exact scale-quality frame which met the stable deadline is a proven
 *    safe lower bound for the newest view.  Quiet handoff starts from that
 *    presentation instead of discarding it and calibrating upward again.
 *
 * No scene pointers, containers, clocks, or callbacks enter this policy.
 */
class BObolLodPresentationPolicy {
public:
    struct Population {
	bool available = false;
	size_t activeFaces = 0;
	uint64_t cadRevision = 0;
    };

    struct Snapshot {
	bool valid = false;
	BObolLodViewEpoch viewEpoch;
	int progressiveCeiling = -1;
	float pointProxyPixelThreshold = 1.0f;
	Population population;
	size_t presentedRenderCost = 0;
    };

    struct RestoreInputs {
	bool orthographic = false;
	bool scaleChanged = false;
	bool retainedMeshPayloads = false;
	BObolLodViewEpoch viewEpoch;
	Population population;
	int currentProgressiveCeiling = -1;
	float currentPointProxyPixelThreshold = 1.0f;
    };

    struct RestoreDecision {
	bool apply = false;
	bool restoredPriorStable = false;
	bool restoredProvenQuality = false;
	bool startHandoff = false;
	bool clearPresentationLimits = false;
	int progressiveCeiling = -1;
	float pointProxyPixelThreshold = 1.0f;
	size_t provenRenderCostFloor = 0;
    };

    struct CompletedPassInputs {
	bool completed = false;
	bool submissionPending = false;
	bool rescanAfterFrame = false;
	bool changedCut = false;
	bool retainedRefinementPending = false;
	bool retainedRefinementBudgetBlocked = false;
    };

    struct CompletedPassDecision {
	bool finishHandoff = false;
	bool requestRetainedRescan = false;
	bool retireRetainedObservation = false;
    };

    void capturePrior(int progressiveCeiling,
	float pointProxyPixelThreshold, const Population &population,
	BObolLodViewEpoch viewEpoch)
    {
	this->priorStableValue = makeSnapshot(
	    progressiveCeiling, pointProxyPixelThreshold, population,
	    viewEpoch, 0);
	this->priorRestoredValue = false;
	this->provenQualityValue = Snapshot();
	this->handoffActiveValue = false;
	this->handoffPresentationRequiredValue = false;
	this->handoffCostFloorValue = 0;
    }

    void noteStableQuality(int progressiveCeiling,
	float pointProxyPixelThreshold, const Population &population,
	BObolLodViewEpoch viewEpoch, size_t presentedRenderCost)
    {
	if (!presentedRenderCost)
	    return;
	const Snapshot candidate = makeSnapshot(
	    progressiveCeiling, pointProxyPixelThreshold, population,
	    viewEpoch, presentedRenderCost);
	if (this->provenQualityValue.valid &&
	    this->provenQualityValue.viewEpoch == viewEpoch &&
	    populationMatches(this->provenQualityValue.population, population) &&
	    this->provenQualityValue.presentedRenderCost > presentedRenderCost)
	    return;
	this->provenQualityValue = candidate;
    }

    RestoreDecision restorePrior(const RestoreInputs &inputs)
    {
	RestoreDecision decision;
	decision.progressiveCeiling = inputs.currentProgressiveCeiling;
	decision.pointProxyPixelThreshold = sanitizeThreshold(
	    inputs.currentPointProxyPixelThreshold);
	if (this->priorRestoredValue) {
	    if (inputs.orthographic && !inputs.scaleChanged &&
		this->priorStableValue.valid &&
		populationMatches(
		    this->priorStableValue.population, inputs.population)) {
		decision.restoredPriorStable = true;
		/* Button-up may restore this snapshot before the 150 ms quiet
		 * debounce.  A frame which then misses the stricter interactive
		 * deadline can install another temporary renderer ceiling.  Reapply
		 * the proven stable limits here; otherwise the prior-restored flag
		 * suppresses both restoration and handoff, leaving that motion
		 * ceiling in an apparently ready quiet view. */
		if (inputs.currentProgressiveCeiling !=
			this->priorStableValue.progressiveCeiling ||
		    std::fabs(sanitizeThreshold(
			inputs.currentPointProxyPixelThreshold) -
			this->priorStableValue.pointProxyPixelThreshold) > 1.0e-6f) {
		    decision.apply = true;
		    decision.progressiveCeiling =
			this->priorStableValue.progressiveCeiling;
		    decision.pointProxyPixelThreshold =
			this->priorStableValue.pointProxyPixelThreshold;
		}
	    } else
		this->priorRestoredValue = false;
	    return decision;
	}
	if (!inputs.orthographic || inputs.scaleChanged ||
	    !this->priorStableValue.valid ||
	    !populationMatches(
		this->priorStableValue.population, inputs.population))
	    return decision;
	decision.apply = true;
	decision.restoredPriorStable = true;
	decision.progressiveCeiling =
	    this->priorStableValue.progressiveCeiling;
	decision.pointProxyPixelThreshold =
	    this->priorStableValue.pointProxyPixelThreshold;
	this->priorRestoredValue = true;
	return decision;
    }

    RestoreDecision beginQuiet(const RestoreInputs &inputs)
    {
	RestoreDecision decision = this->restorePrior(inputs);
	if (!decision.restoredPriorStable && inputs.scaleChanged &&
	    this->provenQualityValue.valid &&
	    this->provenQualityValue.viewEpoch == inputs.viewEpoch &&
	    populationMatches(
		this->provenQualityValue.population, inputs.population)) {
	    decision.apply = true;
	    decision.restoredProvenQuality = true;
	    decision.progressiveCeiling =
		this->provenQualityValue.progressiveCeiling;
	    decision.pointProxyPixelThreshold =
		this->provenQualityValue.pointProxyPixelThreshold;
	    decision.provenRenderCostFloor =
		this->provenQualityValue.presentedRenderCost;
	}

	if (!inputs.retainedMeshPayloads) {
	    decision.clearPresentationLimits = true;
	    decision.apply = true;
	    decision.progressiveCeiling = -1;
	    decision.pointProxyPixelThreshold = 1.0f;
	    this->handoffActiveValue = false;
	    this->handoffPresentationRequiredValue = false;
	    this->handoffCostFloorValue = 0;
	} else {
	    const bool constrained = decision.progressiveCeiling >= 0 ||
		decision.pointProxyPixelThreshold > 1.01f;
	    this->handoffActiveValue =
		!decision.restoredPriorStable && constrained;
	    /* The constrained interaction cut was already presented before the
	     * quiet debounce called beginQuiet().  A deadline recovery uses
	     * armHandoff(true) instead and must prove one successful constrained
	     * frame before its ceiling can be removed. */
	    this->handoffPresentationRequiredValue = false;
	    this->handoffCostFloorValue =
		this->handoffActiveValue && decision.restoredProvenQuality ?
		    decision.provenRenderCostFloor : 0;
	    decision.startHandoff = this->handoffActiveValue;
	}

	this->priorStableValue = Snapshot();
	this->provenQualityValue = Snapshot();
	this->priorRestoredValue = false;
	return decision;
    }

    CompletedPassDecision completePass(const CompletedPassInputs &inputs)
    {
	CompletedPassDecision decision;
	if (!inputs.completed ||
	    inputs.submissionPending || inputs.rescanAfterFrame ||
	    inputs.changedCut || inputs.retainedRefinementBudgetBlocked ||
	    this->handoffPresentationRequiredValue)
	    return decision;
	if (!this->handoffActiveValue) {
	    /* Wanting a richer retained cut is an observation, not a liveness
	     * witness.  With no admitted cut, budget barrier, handoff, or rescan,
	     * the current calibrated cut is terminal until a new external edge
	     * starts another pass. */
	    decision.retireRetainedObservation =
		inputs.retainedRefinementPending;
	    return decision;
	}
	this->handoffActiveValue = false;
	this->handoffPresentationRequiredValue = false;
	this->handoffCostFloorValue = 0;
	decision.finishHandoff = true;
	decision.requestRetainedRescan =
	    inputs.retainedRefinementPending;
	return decision;
    }

    void armHandoff(bool presentationRequired)
    {
	this->handoffActiveValue = true;
	this->handoffPresentationRequiredValue = presentationRequired;
	this->handoffCostFloorValue = 0;
    }
    bool noteFramePresented(void)
    {
	const bool released = this->handoffActiveValue &&
	    this->handoffPresentationRequiredValue;
	if (this->handoffActiveValue)
	    this->handoffPresentationRequiredValue = false;
	return released;
    }
    void cancelHandoff(void)
    {
	this->handoffActiveValue = false;
	this->handoffPresentationRequiredValue = false;
	this->handoffCostFloorValue = 0;
    }

    void viewInvalidated(void)
    {
	this->handoffActiveValue = false;
	this->handoffPresentationRequiredValue = false;
	this->handoffCostFloorValue = 0;
	this->provenQualityValue = Snapshot();
    }

    void reset(void)
    {
	this->priorStableValue = Snapshot();
	this->provenQualityValue = Snapshot();
	this->priorRestoredValue = false;
	this->handoffActiveValue = false;
	this->handoffPresentationRequiredValue = false;
	this->handoffCostFloorValue = 0;
    }

    bool handoffPending(void) const { return this->handoffActiveValue; }
    bool handoffPresentationPending(void) const
    {
	return this->handoffActiveValue &&
	    this->handoffPresentationRequiredValue;
    }
    size_t handoffCostFloor(void) const
    {
	return this->handoffCostFloorValue;
    }
    bool priorSnapshotValid(void) const
    {
	return this->priorStableValue.valid;
    }
    bool priorRestored(void) const { return this->priorRestoredValue; }
    bool provenQualityValid(void) const
    {
	return this->provenQualityValue.valid;
    }
    int provenQualityCeiling(void) const
    {
	return this->provenQualityValue.progressiveCeiling;
    }

private:
    static float sanitizeThreshold(float threshold)
    {
	return std::isfinite(threshold) ?
	    std::max(1.0f, std::min(64.0f, threshold)) : 1.0f;
    }

    static Snapshot makeSnapshot(int progressiveCeiling,
	float pointProxyPixelThreshold, const Population &population,
	BObolLodViewEpoch viewEpoch, size_t presentedRenderCost)
    {
	Snapshot snapshot;
	snapshot.valid = true;
	snapshot.viewEpoch = viewEpoch;
	snapshot.progressiveCeiling =
	    progressiveCeiling < -1 ? -1 :
	    std::min(15, progressiveCeiling);
	snapshot.pointProxyPixelThreshold =
	    sanitizeThreshold(pointProxyPixelThreshold);
	snapshot.population = population;
	snapshot.presentedRenderCost = presentedRenderCost;
	return snapshot;
    }

    static bool populationMatches(const Population &snapshot,
	const Population &current)
    {
	if (current.available)
	    return snapshot.activeFaces == current.activeFaces &&
		snapshot.cadRevision == current.cadRevision;
	return snapshot.activeFaces == 0;
    }

    Snapshot priorStableValue;
    Snapshot provenQualityValue;
    bool priorRestoredValue = false;
    bool handoffActiveValue = false;
    bool handoffPresentationRequiredValue = false;
    size_t handoffCostFloorValue = 0;
};

static_assert(std::is_trivially_copyable<
    BObolLodPresentationPolicy>::value,
    "presentation handoff policy must remain an allocation-free value");

/*
 * Renderer-independent presentation-quality calculations.  These functions
 * convert completed-frame evidence into a camera-local pixel error, small-
 * occurrence aggregation threshold, and reversible PoP prefix ceiling.  They
 * accept only scalar observations so the coordinator cannot acquire a scene
 * pointer or accidentally turn motion feedback into a hierarchy traversal.
 */
class BObolLodQualityPolicy {
public:
    struct StableCapacityEvidence {
	bool proven = false;
	long double renderCostPerSecond = 0.0L;
    };

    static float interactivePixelError(uint64_t renderNanoseconds,
	float targetFps)
    {
	const double targetNanoseconds = targetFrameNanoseconds(targetFps);
	if (!renderNanoseconds)
	    return 4.0f;
	if (static_cast<double>(renderNanoseconds) <=
	    targetNanoseconds * 1.05)
	    return 1.0f;
	const double scale = std::sqrt(
	    static_cast<double>(renderNanoseconds) / targetNanoseconds);
	return static_cast<float>(std::max(2.0,
	    std::min(64.0, 2.0 * scale)));
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
 * Owner-thread result-publication batching.  Applied worker results already
 * own immutable CPU bindings, so several drain quanta may be presented by one
 * frame; nevertheless every unpresented batch must retain either its bounded
 * timer or an explicitly requested frame as a liveness witness.
 */
class BObolLodPublicationPolicy {
public:
    struct Inputs {
	int64_t nowMicroseconds = 0;
	uint64_t observedRenderNanoseconds = 0;
	bool interactive = false;
	bool firstUseful = false;
	bool streamIdle = false;
    };

    struct Decision {
	bool keepPumpAlive = false;
	bool requestFrame = false;
    };

    void noteApplied(size_t count, int64_t nowMicroseconds)
    {
	if (!count)
	    return;
	this->unpresentedCountValue =
	    count > SIZE_MAX - this->unpresentedCountValue ?
		SIZE_MAX : this->unpresentedCountValue + count;
	if (this->firstUnpresentedMicrosecondsValue <= 0)
	    this->firstUnpresentedMicrosecondsValue =
		nowMicroseconds > 0 ? nowMicroseconds : 1;
    }

    Decision decide(const Inputs &inputs)
    {
	Decision decision;
	if (!this->unpresentedCountValue)
	    return decision;
	decision.keepPumpAlive = true;
	if (this->framePendingValue)
	    return decision;
	if (!inputs.firstUseful && !inputs.streamIdle && !due(inputs))
	    return decision;
	this->framePendingValue = true;
	decision.requestFrame = true;
	return decision;
    }

    void frameCompleted(void)
    {
	this->framePendingValue = false;
	this->unpresentedCountValue = 0;
	this->firstUnpresentedMicrosecondsValue = 0;
    }

    void reset(void)
    {
	this->frameCompleted();
    }

    bool pending(void) const
    {
	return this->unpresentedCountValue != 0;
    }

    bool framePending(void) const
    {
	return this->framePendingValue;
    }

    bool awaitingDeadline(void) const
    {
	return this->unpresentedCountValue != 0 &&
	    !this->framePendingValue;
    }

    size_t unpresentedCount(void) const
    {
	return this->unpresentedCountValue;
    }

    int64_t firstUnpresentedMicroseconds(void) const
    {
	return this->firstUnpresentedMicrosecondsValue;
    }

    static int64_t presentationIntervalMicroseconds(
	uint64_t observedRenderNanoseconds, bool interactive)
    {
	const int64_t minimum = interactive ? 33333 : 50000;
	const int64_t maximum = interactive ? 100000 : 250000;
	if (!observedRenderNanoseconds)
	    return minimum;
	const uint64_t observedMicroseconds =
	    observedRenderNanoseconds / 1000ULL;
	const uint64_t scaled = observedMicroseconds >
		static_cast<uint64_t>(maximum) / 2ULL ?
	    static_cast<uint64_t>(maximum) : observedMicroseconds * 2ULL;
	return static_cast<int64_t>(std::max<uint64_t>(
	    static_cast<uint64_t>(minimum),
	    std::min<uint64_t>(static_cast<uint64_t>(maximum), scaled)));
    }

private:
    bool due(const Inputs &inputs) const
    {
	if (!this->unpresentedCountValue ||
	    this->firstUnpresentedMicrosecondsValue <= 0)
	    return false;
	const int64_t interval = presentationIntervalMicroseconds(
	    inputs.observedRenderNanoseconds, inputs.interactive);
	const int64_t elapsed =
	    inputs.nowMicroseconds >= this->firstUnpresentedMicrosecondsValue ?
		inputs.nowMicroseconds -
		    this->firstUnpresentedMicrosecondsValue : interval;
	return elapsed >= interval;
    }

    bool framePendingValue = false;
    size_t unpresentedCountValue = 0;
    int64_t firstUnpresentedMicrosecondsValue = 0;
};

/*
 * Stable resident-prefix compaction admission and continuation state.  The
 * controller owns the demand snapshot and invokes the bounded service cursor;
 * this value owns when that cursor may run and whether its next slice is a
 * continuation of the same demand revision.  Keeping the four related
 * latches together prevents a deadline reset, camera coverage pass, or partial
 * service result from creating an invalid COMPACTING transition.
 */
class BObolLodCompactionPolicy {
public:
    enum class Admission : uint8_t {
	INACTIVE = 0,
	DEFER,
	PLAN
    };

    struct Decision {
	Admission admission = Admission::INACTIVE;
	bool keepPumpAlive = false;
	bool retiredRequest = false;
    };

    struct Inputs {
	bool automatic = false;
	bool interactive = false;
	bool coverageRequired = false;
	bool coverageComplete = false;
	bool coverageProgressPending = false;
	bool settlingPending = false;
	bool realizationPending = false;
	bool submissionPending = false;
	bool serviceAvailable = false;
	int64_t nowMicroseconds = 0;
    };

    void requestAfter(int64_t nowMicroseconds, int64_t delayMicroseconds)
    {
	this->pendingValue = true;
	this->deadlineValue = deadlineAfter(
	    nowMicroseconds, delayMicroseconds);
    }

    void requestImmediate(int64_t nowMicroseconds)
    {
	this->pendingValue = true;
	if (this->deadlineValue <= 0 ||
	    this->deadlineValue > nowMicroseconds)
	    this->deadlineValue = nowMicroseconds;
    }

    void resetForServiceChange(bool serviceAvailable,
	int64_t nowMicroseconds, int64_t delayMicroseconds)
    {
	this->pendingValue = serviceAvailable;
	this->planningValue = false;
	this->demandRevisionValue = 0;
	this->deadlineValue = serviceAvailable ?
	    deadlineAfter(nowMicroseconds, delayMicroseconds) : 0;
    }

    Decision decide(const Inputs &inputs)
    {
	Decision decision;
	if (!inputs.automatic || !this->pendingValue || inputs.interactive ||
	    this->deadlineValue <= 0)
	    return decision;
	/* Incomplete coverage is a prerequisite, not compaction-owned work.  If
	 * its producer has an explicit progress witness, keep the pump alive; if
	 * it has none (LoD disabled, empty view, or definitive fallback), leave
	 * the request dormant until a later view/source event supplies proof.
	 * Otherwise compaction itself becomes an unfulfillable background latch. */
	if (inputs.coverageRequired && !inputs.coverageComplete) {
	    if (inputs.coverageProgressPending ||
		inputs.realizationPending || inputs.submissionPending) {
		decision.admission = Admission::DEFER;
		decision.keepPumpAlive = true;
	    } else {
		this->pendingValue = false;
		this->planningValue = false;
		this->demandRevisionValue = 0;
		this->deadlineValue = 0;
		decision.retiredRequest = true;
	    }
	    return decision;
	}
	if (inputs.settlingPending ||
	    inputs.nowMicroseconds < this->deadlineValue ||
	    inputs.realizationPending || inputs.submissionPending) {
	    decision.admission = Admission::DEFER;
	    decision.keepPumpAlive = true;
	    return decision;
	}
	if (inputs.serviceAvailable)
	    decision.admission = Admission::PLAN;
	return decision;
    }

    bool continues(uint64_t demandRevision) const
    {
	return this->planningValue &&
	    this->demandRevisionValue == demandRevision;
    }

    void finishPlanning(bool complete, uint64_t demandRevision,
	int64_t nowMicroseconds)
    {
	this->pendingValue = !complete;
	this->planningValue = !complete;
	this->demandRevisionValue = complete ? 0 : demandRevision;
	this->deadlineValue = complete ? 0 : nowMicroseconds;
    }

    bool pending(void) const { return this->pendingValue; }
    bool planning(void) const { return this->planningValue; }
    uint64_t demandRevision(void) const
    {
	return this->demandRevisionValue;
    }
    int64_t deadlineMicroseconds(void) const
    {
	return this->deadlineValue;
    }

private:
    static int64_t deadlineAfter(int64_t nowMicroseconds,
	int64_t delayMicroseconds)
    {
	if (delayMicroseconds <= 0)
	    return nowMicroseconds;
	if (nowMicroseconds > INT64_MAX - delayMicroseconds)
	    return INT64_MAX;
	return nowMicroseconds + delayMicroseconds;
    }

    bool pendingValue = false;
    bool planningValue = false;
    uint64_t demandRevisionValue = 0;
    int64_t deadlineValue = 0;
};

/*
 * Owner-thread coordinator state machine.  It stores no strings, allocates no
 * memory, and is never consulted by per-occurrence planning.  Phase mutation
 * has one entry point: dispatch().  Every transition therefore has an explicit
 * cause, a retained input observation, and a progress witness.  This is
 * deliberately stronger than the former phase tracker, whose callers could
 * enter arbitrary phases after independently interpreting a collection of
 * Boolean flags.
 *
 * The view controller owns detailed bounded cursors and renderer policy.  It
 * projects observations into Inputs at named owner-thread boundaries.  The
 * state machine retains both that raw observation and a canonical state:
 * lifecycle events own their edge transitions, pending witnesses cannot be
 * reported as stable, and contradictory compaction/coverage input degrades to
 * coverage or fallback.  Invalid caller combinations are therefore visible in
 * diagnostics without becoming the externally reported state.
 */
class BObolLodStateMachine {
public:
    enum class Phase : uint8_t {
	FALLBACK = 0,
	COVERAGE,
	INTERACTIVE,
	SETTLING,
	STABLE,
	COMPACTING,
	COUNT
    };

    enum class Event : uint8_t {
	INITIALIZE = 0,
	FRAME_COMPLETED,
	WORK_SCHEDULED,
	WORK_PUMPED,
	RESULT_PUBLISHED,
	SERVICE_CHANGED,
	GENERATION_CANCELLED,
	AUTO_SUBMIT_CHANGED,
	VIEW_INVALIDATED,
	POLICY_CHANGED,
	INTERACTION_STARTED,
	INTERACTION_ENDED,
	VIEW_OBSERVED,
	COUNT
    };

    enum Invariant : uint32_t {
	INVARIANT_NONE = 0,
	INVARIANT_COMPLETED_EXCEEDS_VISIBLE = 1u << 0,
	INVARIANT_COMPACTION_WITHOUT_COVERAGE = 1u << 1,
	INVARIANT_STABLE_WITHOUT_COVERAGE = 1u << 2,
	INVARIANT_STABLE_WITH_PENDING_WORK = 1u << 3,
	INVARIANT_INTERACTION_EVENT_WITHOUT_INTERACTION = 1u << 4,
	INVARIANT_CANCEL_EVENT_WITH_GENERATION = 1u << 5,
	INVARIANT_INTERACTION_END_WITH_GESTURE = 1u << 6,
	INVARIANT_RESOURCE_RECOVERY_WITHOUT_PRESSURE = 1u << 7
    };

    /*
     * Allocation-free phase decision inputs.  The view coordinator projects
     * its detailed counters and latches into this compact contract; tests can
     * then exercise the state graph without constructing a GL view or relying
     * on wall-clock timing.
     */
    struct Inputs {
	bool interactive = false;
	bool compacting = false;
	bool coverageActive = false;
	bool coverageComplete = false;
	bool generationActive = false;
	bool settlingWork = false;
	/* A pointer gesture is narrower than the interaction phase: release ends
	 * the gesture while the quiet debounce intentionally remains interactive. */
	bool gestureActive = false;
	/* Completed-frame CPU/GPU observations are coordinator inputs, not HUD-
	 * only diagnostics.  Recovery is a one-shot request for the current
	 * pressure revision; persistent unavoidable pressure may subsequently be
	 * reported as a stable, memory-limited terminal state. */
	bool cpuMemoryPressure = false;
	bool gpuMemoryPressure = false;
	bool resourceRecoveryPending = false;
    };

    struct Witness {
	uint64_t sequence = 0;
	uint64_t viewEpoch = 0;
	uint64_t policyEpoch = 0;
	uint64_t renderSerial = 0;
	uint64_t activeGeneration = 0;
	uint64_t residentDemandRevision = 0;
	uint64_t resourcePressureRevision = 0;
	size_t visibleCount = 0;
	size_t completedCount = 0;
	size_t pendingCount = 0;
    };

    struct Transition {
	Phase previous = Phase::FALLBACK;
	Phase current = Phase::FALLBACK;
	Event event = Event::INITIALIZE;
	uint32_t invariantMask = INVARIANT_NONE;
	bool changed = false;
	bool progressed = false;
	bool canonicalized = false;
    };

    BObolLodStateMachine(void) = default;

    Phase currentPhase(void) const
    {
	return this->phase;
    }

    uint64_t phaseTransitionSerial(void) const
    {
	return this->transitionSerial;
    }

    uint64_t dispatchSerial(void) const
    {
	return this->dispatchCount;
    }

    Event lastDispatchedEvent(void) const
    {
	return this->lastEvent;
    }

    const Inputs &lastInputs(void) const
    {
	return this->stateInputsValue;
    }

    const Inputs &lastObservedInputs(void) const
    {
	return this->observedInputsValue;
    }

    const Witness &lastObservedWitness(void) const
    {
	return this->observedWitnessValue;
    }

    uint32_t lastInvariantMask(void) const
    {
	return this->invariantMask;
    }

    uint32_t invariantHistoryMask(void) const
    {
	return this->violationHistoryMask;
    }

    uint64_t invariantViolationCount(void) const
    {
	return this->violationCount;
    }

    uint64_t stagnantDispatchCount(void) const
    {
	return this->stagnantCount;
    }

    const Witness &phaseWitness(Phase requestedPhase) const
    {
	const size_t index = static_cast<size_t>(requestedPhase);
	return this->witnesses[index < this->witnesses.size() ?
	    index : static_cast<size_t>(Phase::FALLBACK)];
    }

    Transition dispatch(Event event, const Inputs &nextInputs,
	const Witness &nextWitness)
    {
	Transition result;
	result.previous = this->phase;
	result.event = event;
	result.invariantMask = validate(event, nextInputs, nextWitness);
	Inputs authoritativeInputs = canonicalInputs(event, nextInputs,
	    nextWitness);
	Witness authoritativeWitness = canonicalWitness(nextWitness);
	result.canonicalized = !sameInputs(nextInputs, authoritativeInputs) ||
	    !sameProgress(nextWitness, authoritativeWitness);
	this->dispatchCount++;
	if (this->dispatchCount == 0)
	    this->dispatchCount++;
	this->lastEvent = event;
	this->observedInputsValue = nextInputs;
	this->observedWitnessValue = nextWitness;
	this->stateInputsValue = authoritativeInputs;
	this->invariantMask = result.invariantMask;
	if (result.invariantMask != INVARIANT_NONE) {
	    this->violationHistoryMask |= result.invariantMask;
	    this->violationCount++;
	    if (this->violationCount == 0)
		this->violationCount++;
	}

	const Phase nextPhase = phaseFor(
	    authoritativeInputs, authoritativeWitness);
	result.changed = this->phase != nextPhase;
	result.progressed = this->transitionTo(
	    nextPhase, authoritativeWitness);
	result.current = this->phase;
	if (result.progressed)
	    this->stagnantCount = 0;
	else if (this->stagnantCount != UINT64_MAX)
	    this->stagnantCount++;
	return result;
    }

private:
    static Phase phaseFor(const Inputs &inputs, const Witness &witness)
    {
	/* Interaction is the only phase which may interrupt every other phase.
	 * Minimum coverage precedes compaction and richer refinement. */
	if (inputs.interactive || inputs.gestureActive)
	    return Phase::INTERACTIVE;
	if (inputs.coverageActive || !inputs.coverageComplete)
	    return inputs.generationActive ? Phase::COVERAGE :
		Phase::FALLBACK;
	if (inputs.compacting || inputs.resourceRecoveryPending)
	    return Phase::COMPACTING;
	if (inputs.settlingWork || witness.pendingCount != 0)
	    return Phase::SETTLING;
	return Phase::STABLE;
    }

    static Inputs canonicalInputs(Event event, const Inputs &observed,
	const Witness &witness)
    {
	Inputs state = observed;
	/* These named lifecycle edges are commands, not advisory labels. */
	if (event == Event::INTERACTION_STARTED) {
	    state.gestureActive = true;
	    state.interactive = true;
	} else if (event == Event::INTERACTION_ENDED) {
	    state.gestureActive = false;
	}
	else if (event == Event::GENERATION_CANCELLED)
	    state.generationActive = false;

	/* Coverage is a safety prerequisite for memory compaction. */
	if (!state.coverageComplete)
	    state.compacting = false;
	if (!state.cpuMemoryPressure && !state.gpuMemoryPressure)
	    state.resourceRecoveryPending = false;
	/* A retained progress witness is sufficient proof that the view is not
	 * stable, even if a caller omitted its redundant settling latch. */
	if (witness.pendingCount != 0)
	    state.settlingWork = true;
	return state;
    }

    static Witness canonicalWitness(const Witness &observed)
    {
	Witness state = observed;
	if (state.completedCount > state.visibleCount)
	    state.completedCount = state.visibleCount;
	return state;
    }

    static uint32_t validate(Event event, const Inputs &inputs,
	const Witness &witness)
    {
	uint32_t violations = INVARIANT_NONE;
	if (witness.completedCount > witness.visibleCount)
	    violations |= INVARIANT_COMPLETED_EXCEEDS_VISIBLE;

	if (inputs.compacting && !inputs.coverageComplete)
	    violations |= INVARIANT_COMPACTION_WITHOUT_COVERAGE;
	const Phase target = phaseFor(inputs, witness);
	if (target == Phase::STABLE && !inputs.coverageComplete)
	    violations |= INVARIANT_STABLE_WITHOUT_COVERAGE;
	if (!inputs.interactive && !inputs.gestureActive && !inputs.compacting &&
	    inputs.coverageComplete && !inputs.coverageActive &&
	    !inputs.settlingWork && witness.pendingCount != 0)
	    violations |= INVARIANT_STABLE_WITH_PENDING_WORK;
	if (event == Event::INTERACTION_STARTED &&
	    (!inputs.interactive || !inputs.gestureActive))
	    violations |= INVARIANT_INTERACTION_EVENT_WITHOUT_INTERACTION;
	if (event == Event::INTERACTION_ENDED && inputs.gestureActive)
	    violations |= INVARIANT_INTERACTION_END_WITH_GESTURE;
	if (event == Event::GENERATION_CANCELLED && inputs.generationActive)
	    violations |= INVARIANT_CANCEL_EVENT_WITH_GENERATION;
	if (inputs.resourceRecoveryPending &&
	    !inputs.cpuMemoryPressure && !inputs.gpuMemoryPressure)
	    violations |= INVARIANT_RESOURCE_RECOVERY_WITHOUT_PRESSURE;
	return violations;
    }

    static bool sameProgress(const Witness &lhs, const Witness &rhs)
    {
	return lhs.viewEpoch == rhs.viewEpoch &&
	    lhs.policyEpoch == rhs.policyEpoch &&
	    lhs.renderSerial == rhs.renderSerial &&
	    lhs.activeGeneration == rhs.activeGeneration &&
	    lhs.residentDemandRevision == rhs.residentDemandRevision &&
	    lhs.resourcePressureRevision == rhs.resourcePressureRevision &&
	    lhs.visibleCount == rhs.visibleCount &&
	    lhs.completedCount == rhs.completedCount &&
	    lhs.pendingCount == rhs.pendingCount;
    }

    static bool sameInputs(const Inputs &lhs, const Inputs &rhs)
    {
	return lhs.interactive == rhs.interactive &&
	    lhs.compacting == rhs.compacting &&
	    lhs.coverageActive == rhs.coverageActive &&
	    lhs.coverageComplete == rhs.coverageComplete &&
	    lhs.generationActive == rhs.generationActive &&
	    lhs.settlingWork == rhs.settlingWork &&
	    lhs.gestureActive == rhs.gestureActive &&
	    lhs.cpuMemoryPressure == rhs.cpuMemoryPressure &&
	    lhs.gpuMemoryPressure == rhs.gpuMemoryPressure &&
	    lhs.resourceRecoveryPending == rhs.resourceRecoveryPending;
    }

    bool transitionTo(Phase nextPhase, const Witness &nextWitness)
    {
	const size_t index = static_cast<size_t>(nextPhase);
	if (index >= this->witnesses.size())
	    return false;
	bool progressed = false;

	if (this->phase != nextPhase) {
	    this->phase = nextPhase;
	    this->transitionSerial++;
	    if (this->transitionSerial == 0)
		this->transitionSerial++;
	    progressed = true;
	}

	Witness &recorded = this->witnesses[index];
	if (!recorded.sequence ||
	    !BObolLodStateMachine::sameProgress(recorded, nextWitness)) {
	    recorded = nextWitness;
	    this->progressSerial++;
	    if (this->progressSerial == 0)
		this->progressSerial++;
	    recorded.sequence = this->progressSerial;
	    progressed = true;
	}
	return progressed;
    }

    Phase phase = Phase::FALLBACK;
    Event lastEvent = Event::INITIALIZE;
    Inputs observedInputsValue {};
    Inputs stateInputsValue {};
    Witness observedWitnessValue {};
    uint64_t transitionSerial = 0;
    uint64_t progressSerial = 0;
    uint64_t dispatchCount = 0;
    uint64_t violationCount = 0;
    uint64_t stagnantCount = 0;
    uint32_t invariantMask = INVARIANT_NONE;
    uint32_t violationHistoryMask = INVARIANT_NONE;
    std::array<Witness, static_cast<size_t>(Phase::COUNT)> witnesses {};
};

#endif /* LIBBOBOL_LOD_COORDINATOR_PRIVATE_H */
