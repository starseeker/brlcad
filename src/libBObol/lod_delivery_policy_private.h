/*        L O D _ D E L I V E R Y _ P O L I C Y _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_DELIVERY_POLICY_PRIVATE_H
#define LIBBOBOL_LOD_DELIVERY_POLICY_PRIVATE_H

#include "common.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>

class BObolLodDeliveryPolicy {
public:
    static int visibleFirstCut(int currentPresentationCut,
	int admittedPresentationCut, int desiredResidentCut,
	bool expensiveSpeculativeDelivery)
    {
	if (expensiveSpeculativeDelivery &&
	    admittedPresentationCut > currentPresentationCut &&
	    desiredResidentCut > admittedPresentationCut)
	    return admittedPresentationCut;
	return desiredResidentCut;
    }
};

/*
 * Database discovery is producer/consumer work rather than a presentation
 * frame.  A fixed 4 ms owner-thread slice leaves more than half of a fast GL
 * frame unused and makes the 16 ms host timer the throughput ceiling at
 * 50k--150k occurrences.  Conversely, spending that larger slice while the
 * user is moving would directly worsen input latency.  Expand only the
 * controller-owned quiet default; an 8 ms cooperative slice remains below the
 * host's 16 ms timer even when a software presentation is independently slow.
 */
class BObolLodAvailabilityScheduler {
public:
    enum class CompletedPassSuccessor : uint8_t {
	NONE = 0,
	COMPLETE_RESIDENCY_DRAIN,
	YIELD_TO_RESIDENT_GROWTH,
	PRESENT_POINT_CALIBRATION,
	CALIBRATE_CAPACITY
    };

    static uint64_t effectiveMicroseconds(uint64_t configuredMicroseconds,
	bool controllerOwnedDefault, bool interactive)
    {
	if (!controllerOwnedDefault ||
	    configuredMicroseconds <= interactiveSliceMicroseconds())
	    return configuredMicroseconds;
	if (interactive)
	    return std::min(configuredMicroseconds,
		interactiveSliceMicroseconds());
	return configuredMicroseconds;
    }

    static bool canRetainPresentation(bool residencyDrainActive,
	bool retainedPayload, bool sameAsset, bool activeCutPreserved,
	bool richerResidentPrefix)
    {
	return residencyDrainActive && retainedPayload && sameAsset &&
	    activeCutPreserved && richerResidentPrefix;
    }

    static bool allocationPopulationSettled(bool providerInventorySettled,
	bool serviceStreamIdle, bool resultDeliveryIdle,
	bool growthTransactionPending, bool residencyDrainActive)
    {
	return providerInventorySettled && serviceStreamIdle &&
	    resultDeliveryIdle && !growthTransactionPending &&
	    !residencyDrainActive;
    }

    /* A resident drain owns immutable-prefix availability, not the cut drawn
     * from an already available prefix.  A completed occurrence-allocation
     * plan must therefore remain able to reconcile presentation cuts while a
     * drain is active; otherwise the global renderer ceiling and the
     * occurrence plan can wait on one another indefinitely.  Ordinary demand
     * changes still defer downgrades until the drain has retired. */
    static bool presentationCutDowngradeAllowed(bool interactive,
	bool gestureActive, bool residencyDrainActive,
	bool scaleDemandChanged, bool retainedAllocation)
    {
	if (interactive || gestureActive)
	    return false;
	if (retainedAllocation)
	    return true;
	return scaleDemandChanged && !residencyDrainActive;
    }

    /* Structural repair and resident growth both consume the shared source
     * cursor.  A newly loaded immutable suffix must first receive its one
     * occurrence-allocation pass; otherwise structural repair repeatedly
     * inspects the pre-growth framebuffer while preventing the transaction
     * capable of making that frontier drawable. */
    static bool structuralRepairMayOwn(bool residentGrowthPending,
	bool residencyDrainActive)
    {
	return !residentGrowthPending && !residencyDrainActive;
    }

    /* Resident growth owns availability population before ordinary capacity
     * calibration may restart the shared submission cursor.  Triangle
     * recovery then owns its finite retained-allocation budget until its
     * coherent pass is presented or proved unchanged.  Treating residual
     * recovery quality debt as an ordinary capacity block restarts the same
     * all-scene pass indefinitely. */
    static CompletedPassSuccessor completedPassSuccessor(bool completed,
	bool residencyDrainActive, bool residentGrowthPending,
	bool pointTriangleRecovery, bool pointCalibrationPending,
	bool capacityBlocked)
    {
	if (!completed)
	    return CompletedPassSuccessor::NONE;
	if (residencyDrainActive)
	    return CompletedPassSuccessor::COMPLETE_RESIDENCY_DRAIN;
	if (residentGrowthPending)
	    return CompletedPassSuccessor::YIELD_TO_RESIDENT_GROWTH;
	if (pointTriangleRecovery)
	    return CompletedPassSuccessor::NONE;
	if (pointCalibrationPending)
	    return CompletedPassSuccessor::PRESENT_POINT_CALIBRATION;
	if (capacityBlocked)
	    return CompletedPassSuccessor::CALIBRATE_CAPACITY;
	return CompletedPassSuccessor::NONE;
    }

    static bool residentGrowthOwnsSuccessor(
	CompletedPassSuccessor successor)
    {
	return successor == CompletedPassSuccessor::COMPLETE_RESIDENCY_DRAIN ||
	    successor == CompletedPassSuccessor::YIELD_TO_RESIDENT_GROWTH;
    }

private:
    static constexpr uint64_t interactiveSliceMicroseconds(void)
    {
	return 4000;
    }
};

/* Own every mutable fact associated with presenting an applied LoD change.
 * Applied worker results already own immutable CPU bindings, so several drain
 * quanta may share one frame.  Result batching and frame barriers nevertheless
 * name the same revision, retain one timer or requested-frame witness, and
 * retire on the same completed CAD frame. */
class BObolLodPresentationTransaction {
public:
    enum Reason : uint32_t {
	REASON_NONE = 0,
	REASON_CUT_PRESENTATION = 1u << 0,
	REASON_RESULT_PUBLICATION = 1u << 1,
	REASON_RESIDENT_REFINEMENT = 1u << 2,
	REASON_CALIBRATION = 1u << 3,
	REASON_HANDOFF = 1u << 4
    };

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

    struct Completion {
	bool retired = false;
	bool stale = false;
	uint32_t reasons = REASON_NONE;
    };

    void noteApplied(size_t count, int64_t nowMicroseconds,
	uint64_t viewEpoch, uint64_t policyEpoch)
    {
	if (!count)
	    return;
	this->begin(viewEpoch, policyEpoch);
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
	if (this->publicationFramePendingValue)
	    return decision;
	if (!inputs.firstUseful && !inputs.streamIdle && !due(inputs))
	    return decision;
	this->publicationFramePendingValue = true;
	decision.requestFrame = true;
	return decision;
    }

    bool arm(Reason reason, uint64_t completedRenderSerial,
	uint64_t viewEpoch, uint64_t policyEpoch)
    {
	if (reason == REASON_NONE)
	    return false;
	this->begin(viewEpoch, policyEpoch);
	this->reasonMask |= static_cast<uint32_t>(reason);
	if (this->barrierPendingValue)
	    return false;
	this->barrierPendingValue = true;
	this->requiredRenderSerialValue = completedRenderSerial == UINT64_MAX ?
	    UINT64_MAX : completedRenderSerial + 1;
	return true;
    }

    Completion complete(uint64_t completedRenderSerial,
	uint64_t viewEpoch, uint64_t policyEpoch)
    {
	Completion result;
	if (!this->active())
	    return result;
	if (this->viewEpochValue != viewEpoch ||
	    this->policyEpochValue != policyEpoch) {
	    result.stale = true;
	    this->clear();
	    return result;
	}
	if (this->barrierPendingValue &&
	    completedRenderSerial < this->requiredRenderSerialValue)
	    return result;
	result.retired = this->barrierPendingValue;
	result.reasons = this->reasonMask;
	this->clear();
	return result;
    }

    void reset(void)
    {
	this->clear();
    }

    bool barrierPending(void) const
    {
	return this->barrierPendingValue;
    }

    bool publicationPending(void) const
    {
	return this->unpresentedCountValue != 0;
    }

    bool publicationFramePending(void) const
    {
	return this->publicationFramePendingValue;
    }

    bool publicationAwaitingFrameRequest(void) const
    {
	/* Once publicationFramePendingValue is true, the render request is the
	 * only remaining liveness witness.  It must not be classified as a future
	 * producer which permits that same request to disappear. */
	return this->unpresentedCountValue != 0 &&
	    !this->publicationFramePendingValue;
    }

    size_t unpresentedCount(void) const
    {
	return this->unpresentedCountValue;
    }

    int64_t firstUnpresentedMicroseconds(void) const
    {
	return this->firstUnpresentedMicrosecondsValue;
    }

    uint32_t reasons(void) const
    {
	return this->reasonMask;
    }

    uint64_t requiredRenderSerial(void) const
    {
	return this->barrierPendingValue ?
	    this->requiredRenderSerialValue : 0;
    }

    uint64_t viewEpoch(void) const
    {
	return this->viewEpochValue;
    }

    uint64_t policyEpoch(void) const
    {
	return this->policyEpochValue;
    }

    uint64_t sequence(void) const
    {
	return this->sequenceValue;
    }

    static int64_t presentationIntervalMicroseconds(
	uint64_t observedRenderNanoseconds, bool interactive)
    {
	const int64_t minimum = interactive ?
	    interactiveMinimumIntervalMicroseconds() :
	    quietMinimumIntervalMicroseconds();
	const int64_t maximum = interactive ?
	    interactiveMaximumIntervalMicroseconds() :
	    quietMaximumIntervalMicroseconds();
	if (!observedRenderNanoseconds)
	    return minimum;
	const uint64_t observedMicroseconds =
	    observedRenderNanoseconds / nanosecondsPerMicrosecond();
	const uint64_t scaled = observedMicroseconds >
		static_cast<uint64_t>(maximum) / intervalScale() ?
	    static_cast<uint64_t>(maximum) :
	    observedMicroseconds * intervalScale();
	return static_cast<int64_t>(std::max<uint64_t>(
	    static_cast<uint64_t>(minimum),
	    std::min<uint64_t>(static_cast<uint64_t>(maximum), scaled)));
    }

private:
    static constexpr int64_t interactiveMinimumIntervalMicroseconds()
    {
	return 33333;
    }

    static constexpr int64_t quietMinimumIntervalMicroseconds()
    {
	return 50000;
    }

    static constexpr int64_t interactiveMaximumIntervalMicroseconds()
    {
	return 100000;
    }

    static constexpr int64_t quietMaximumIntervalMicroseconds()
    {
	return 250000;
    }

    static constexpr uint64_t intervalScale()
    {
	return 2;
    }

    static constexpr uint64_t nanosecondsPerMicrosecond()
    {
	return 1000;
    }

    bool active(void) const
    {
	return this->barrierPendingValue || this->unpresentedCountValue != 0;
    }

    void begin(uint64_t viewEpoch, uint64_t policyEpoch)
    {
	if (this->active() &&
	    (this->viewEpochValue != viewEpoch ||
	     this->policyEpochValue != policyEpoch))
	    this->clear();
	if (this->active())
	    return;
	this->viewEpochValue = viewEpoch;
	this->policyEpochValue = policyEpoch;
	if (++this->sequenceValue == 0)
	    ++this->sequenceValue;
    }

    void clear(void)
    {
	this->barrierPendingValue = false;
	this->publicationFramePendingValue = false;
	this->reasonMask = REASON_NONE;
	this->requiredRenderSerialValue = 0;
	this->viewEpochValue = 0;
	this->policyEpochValue = 0;
	this->unpresentedCountValue = 0;
	this->firstUnpresentedMicrosecondsValue = 0;
    }

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
    bool barrierPendingValue = false;
    bool publicationFramePendingValue = false;
    uint32_t reasonMask = REASON_NONE;
    uint64_t requiredRenderSerialValue = 0;
    uint64_t viewEpochValue = 0;
    uint64_t policyEpochValue = 0;
    uint64_t sequenceValue = 0;
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
	/* This is an admission edge, not a debounce timer.  Camera settling,
	 * coverage, realization, and submission are independent gates in
	 * decide().  Sliding the deadline on every equivalent retained-demand
	 * pass can postpone stable reclamation forever when an idempotent
	 * retarget continues to report an updated cut. */
	if (!this->pendingValue)
	    this->deadlineValue = deadlineAfter(
		nowMicroseconds, delayMicroseconds);
	this->pendingValue = true;
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

#endif /* LIBBOBOL_LOD_DELIVERY_POLICY_PRIVATE_H */
