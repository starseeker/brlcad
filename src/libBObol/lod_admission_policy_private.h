/*      L O D _ A D M I S S I O N _ P O L I C Y _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_ADMISSION_POLICY_PRIVATE_H
#define LIBBOBOL_LOD_ADMISSION_POLICY_PRIVATE_H

#include "common.h"

#include "lod_capacity_policy_private.h"
#include "lod_revision_private.h"
#include "lod_scene_evidence_private.h"
#include "lod_view_policy_private.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>

class BObolLodProducerPolicy {
public:
    static bool canProduceGeometry(bool submissionPending,
	bool submissionPausedByPresentation, bool providerPending,
	bool servicePending)
    {
	return providerPending || servicePending ||
	    (submissionPending && !submissionPausedByPresentation);
    }

    static bool ownsFutureFrame(bool submissionPending,
	bool submissionPausedByPresentation, bool providerPending,
	bool servicePending, bool publicationAwaitingFrameRequest)
    {
	return publicationAwaitingFrameRequest || canProduceGeometry(
	    submissionPending, submissionPausedByPresentation,
	    providerPending, servicePending);
    }

};

/*
 * Single owner for availability edges arriving from source providers and the
 * LoD service.  The atomics are the worker-to-owner result notification; all
 * remaining values are presentation-thread state.  Keeping result age,
 * provider terminality, inventory coalescing, and resident-growth obligation
 * together makes the scheduler's progress witness observable in one place.
 *
 * This ledger delays only a pure inventory continuation.  A camera/policy or
 * source-identity change, an active submission cursor, and the producer's
 * final pending->idle edge all bypass it.  The bounded age means a producer
 * which never fills a large batch still makes progress.
 */
class BObolLodAvailabilityLedger {
public:
    void noteResultsReady(int64_t nowMicroseconds)
    {
	this->resultsPendingValue.store(true);
	(void)this->ensureFirstResultReady(nowMicroseconds);
    }

    void setResultsPending(bool pending)
    {
	this->resultsPendingValue.store(pending);
    }

    bool resultsPending(void) const
    {
	return this->resultsPendingValue.load();
    }

    int64_t ensureFirstResultReady(int64_t nowMicroseconds)
    {
	const int64_t now = nowMicroseconds > 0 ? nowMicroseconds :
	    minimumTimestamp();
	int64_t expected = 0;
	(void)this->firstResultReadyMicrosecondsValue.
	    compare_exchange_strong(expected, now);
	return this->firstResultReadyMicrosecondsValue.load();
    }

    int64_t firstResultReadyMicroseconds(void) const
    {
	return this->firstResultReadyMicrosecondsValue.load();
    }

    void clearFirstResultReady(void)
    {
	this->firstResultReadyMicrosecondsValue.store(0);
    }

    void resetResultQueue(void)
    {
	this->resultsPendingValue.store(false);
	this->clearFirstResultReady();
    }

    void setProviderPendingCount(size_t count)
    {
	this->providerPendingCountValue = count;
    }

    size_t providerPendingCount(void) const
    {
	return this->providerPendingCountValue;
    }

    bool deferInventoryDelta(bool inventoryChanged, bool providerPending,
	bool submissionPending, bool initialSubmissionComplete,
	bool interactive, int64_t nowMicroseconds)
    {
	if (!inventoryChanged || !providerPending || submissionPending ||
	    !initialSubmissionComplete) {
	    this->inventoryFirstPendingMicrosecondsValue = 0;
	    return false;
	}
	const int64_t now = nowMicroseconds > 0 ? nowMicroseconds :
	    minimumTimestamp();
	if (this->inventoryFirstPendingMicrosecondsValue <= 0)
	    this->inventoryFirstPendingMicrosecondsValue = now;
	const int64_t limit = inventoryCoalescingLimit(interactive);
	const int64_t age = now >= this->inventoryFirstPendingMicrosecondsValue ?
	    now - this->inventoryFirstPendingMicrosecondsValue : limit;
	return age < limit;
    }

    void commitInventoryDelta(void)
    {
	this->inventoryFirstPendingMicrosecondsValue = 0;
    }

    int64_t inventoryFirstPendingMicroseconds(void) const
    {
	return this->inventoryFirstPendingMicrosecondsValue;
    }

    void noteRicherPrefixAvailable(void)
    {
	switch (this->residentGrowthPhaseValue) {
	    case ResidentGrowthPhase::IDLE:
	    case ResidentGrowthPhase::REALLOCATION_READY:
		this->residentGrowthPhaseValue =
		    ResidentGrowthPhase::DRAIN_REQUIRED;
		break;
	    case ResidentGrowthPhase::DRAIN_ACTIVE:
		this->residentGrowthPhaseValue =
		    ResidentGrowthPhase::DRAIN_ACTIVE_DIRTY;
		break;
	    case ResidentGrowthPhase::DRAIN_REQUIRED:
	    case ResidentGrowthPhase::DRAIN_ACTIVE_DIRTY:
		break;
	}
    }

    bool beginResidencyDrainIfReady(bool automatic, bool streamIdle,
	bool workAllowed)
    {
	if (this->residentGrowthPhaseValue !=
		ResidentGrowthPhase::DRAIN_REQUIRED ||
	    !automatic || !streamIdle ||
	    !workAllowed)
	    return false;
	this->residentGrowthPhaseValue = ResidentGrowthPhase::DRAIN_ACTIVE;
	return true;
    }

    bool completeResidencyDrain(void)
    {
	if (this->residentGrowthPhaseValue ==
		ResidentGrowthPhase::DRAIN_ACTIVE) {
	    this->residentGrowthPhaseValue =
		ResidentGrowthPhase::REALLOCATION_READY;
	    return true;
	}
	if (this->residentGrowthPhaseValue ==
		ResidentGrowthPhase::DRAIN_ACTIVE_DIRTY) {
	    this->residentGrowthPhaseValue =
		ResidentGrowthPhase::DRAIN_REQUIRED;
	    return true;
	}
	return false;
    }

    void interruptResidencyDrain(void)
    {
	if (this->residencyDrainActive())
	    this->residentGrowthPhaseValue =
		ResidentGrowthPhase::DRAIN_REQUIRED;
    }

    bool consumeResidentGrowthIfReady(bool automatic, bool streamIdle,
	bool allocationAllowed)
    {
	if (this->residentGrowthPhaseValue !=
		ResidentGrowthPhase::REALLOCATION_READY ||
	    !automatic || !streamIdle || !allocationAllowed)
	    return false;
	this->residentGrowthPhaseValue = ResidentGrowthPhase::IDLE;
	return true;
    }

    bool residentGrowthPending(void) const
    {
	return this->residentGrowthPhaseValue != ResidentGrowthPhase::IDLE;
    }

    bool transferResidentGrowthToCapacityCandidate(void)
    {
	if (!this->residentGrowthPending())
	    return false;
	/* The capacity candidate now owns both the already published prefix and
	 * any suffix task still in flight.  Its population-settlement gate waits
	 * for the service stream, so retaining a second growth phase adds no safety
	 * and can only create a circular scheduler dependency. */
	this->residentGrowthPhaseValue = ResidentGrowthPhase::IDLE;
	return true;
    }

    bool residencyDrainRequired(void) const
    {
	return this->residentGrowthPhaseValue ==
		ResidentGrowthPhase::DRAIN_REQUIRED ||
	    this->residentGrowthPhaseValue ==
		ResidentGrowthPhase::DRAIN_ACTIVE_DIRTY;
    }

    bool residencyDrainActive(void) const
    {
	return this->residentGrowthPhaseValue ==
		ResidentGrowthPhase::DRAIN_ACTIVE ||
	    this->residentGrowthPhaseValue ==
		ResidentGrowthPhase::DRAIN_ACTIVE_DIRTY;
    }

    void resetResidentGrowth(void)
    {
	this->residentGrowthPhaseValue = ResidentGrowthPhase::IDLE;
    }

private:
    enum class ResidentGrowthPhase : uint8_t {
	IDLE,
	DRAIN_REQUIRED,
	DRAIN_ACTIVE,
	DRAIN_ACTIVE_DIRTY,
	REALLOCATION_READY
    };

    static constexpr int64_t minimumTimestamp(void) { return 1; }
    static constexpr int64_t inventoryCoalescingLimit(bool interactive)
    {
	return interactive ? 100000 : 250000;
    }

    std::atomic<bool> resultsPendingValue {false};
    std::atomic<int64_t> firstResultReadyMicrosecondsValue {0};
    size_t providerPendingCountValue = 0;
    int64_t inventoryFirstPendingMicrosecondsValue = 0;
    ResidentGrowthPhase residentGrowthPhaseValue = ResidentGrowthPhase::IDLE;
};

/*
 * Bracket the camera-local small-occurrence aggregation cut using actual
 * completed/interrupted draw evidence.  A larger threshold is cheaper.  The
 * old controller alternated an aggressive pressure multiplication with a
 * fixed 0.75 relaxation step; near the renderer deadline those two rules
 * formed a permanent coarse/fine oscillation.  Remembering the finest proven
 * safe cut and the coarsest proven unsafe cut turns the same feedback into a
 * convergent geometric search.
 */
/**
 * Sole owner of the quiet point-presentation quality transition.
 *
 * Adaptive calibration measures a candidate point cut.  Handoff confirmation
 * presents an already selected point cut without adapting it.  Triangle
 * recovery needs retained admission to run.  Representing those purposes with
 * independent booleans allowed more than one to be armed, leaving a pending
 * controller with no enabled producer or letting confirmation undo the cut it
 * was meant to prove.  Recovery has priority because its coherent retained
 * allocation must finish before a new presentation-only calibration can be
 * meaningful.
 */
class BObolLodPointQualityPhase {
public:
    enum class Phase {
	IDLE,
	ADAPTIVE_CALIBRATION,
	STATIC_CALIBRATION,
	HANDOFF_CONFIRMATION,
	TRIANGLE_RECOVERY
    };

    bool presentationPending(void) const
    {
	return this->phase == Phase::ADAPTIVE_CALIBRATION ||
	    this->phase == Phase::STATIC_CALIBRATION ||
	    this->phase == Phase::HANDOFF_CONFIRMATION;
    }

    bool adaptiveCalibrationPending(void) const
    {
	return this->phase == Phase::ADAPTIVE_CALIBRATION;
    }

    bool staticCalibrationPending(void) const
    {
	return this->phase == Phase::STATIC_CALIBRATION;
    }

    bool handoffConfirmationPending(void) const
    {
	return this->phase == Phase::HANDOFF_CONFIRMATION;
    }

    bool triangleRecoveryPending(void) const
    {
	return this->phase == Phase::TRIANGLE_RECOVERY;
    }

    bool pending(void) const
    {
	return this->phase != Phase::IDLE;
    }

    bool requiresRecoveryPresentation(bool completedPass, bool changedCut,
	bool submissionPending, bool presentationPending) const
    {
	return this->phase == Phase::TRIANGLE_RECOVERY && completedPass &&
	    changedCut && !submissionPending && !presentationPending;
    }

    bool recoveryAdmissionPending(bool retainedRefinementPending,
	bool allocationCoversCurrentPopulation) const
    {
	/* A constrained recovery plan may deliberately leave richer pixel-target
	 * demand.  Once that exact plan covers the current occurrence population,
	 * the remaining annotation is quality debt beyond the recovery budget—not
	 * unfinished recovery work.  Requiring it to clear starts the same no-op
	 * scene scan forever. */
	return this->phase == Phase::TRIANGLE_RECOVERY &&
	    retainedRefinementPending && !allocationCoversCurrentPopulation;
    }

    void requestCalibration(void)
    {
	if (this->phase != Phase::TRIANGLE_RECOVERY &&
	    this->phase != Phase::STATIC_CALIBRATION)
	    this->phase = Phase::ADAPTIVE_CALIBRATION;
    }

    /* Event-driven static refinement is allowed to retain a framebuffer
     * whose one-time render exceeds the ordinary streaming cadence.  Keep
     * that deadline class in the sum-typed point owner until its exact
     * candidate frame either settles or schedules its next bounded step. */
    void requestStaticCalibration(void)
    {
	if (this->phase != Phase::TRIANGLE_RECOVERY)
	    this->phase = Phase::STATIC_CALIBRATION;
    }

    void requestHandoffConfirmation(void)
    {
	if (this->phase != Phase::TRIANGLE_RECOVERY)
	    this->phase = Phase::HANDOFF_CONFIRMATION;
    }

    void completeCalibration(void)
    {
	if (this->presentationPending())
	    this->phase = Phase::IDLE;
    }

    void beginTriangleRecovery(void)
    {
	this->phase = Phase::TRIANGLE_RECOVERY;
    }

    void completeTriangleRecovery(uint64_t saturatedPlanSerial = 0)
    {
	if (this->phase == Phase::TRIANGLE_RECOVERY) {
	    this->phase = Phase::IDLE;
	    this->saturatedPlanSerialValue = saturatedPlanSerial;
	}
    }

    bool triangleRecoverySaturatedBy(uint64_t planSerial) const
    {
	return planSerial != 0 &&
	    this->saturatedPlanSerialValue == planSerial;
    }

    void reset(void)
    {
	this->phase = Phase::IDLE;
	this->saturatedPlanSerialValue = 0;
    }

private:
    Phase phase = Phase::IDLE;
    /* A completed constrained allocation must not be resubmitted merely
     * because presentation-only point calibration follows it.  The active
     * allocation's own revision stamp supplies freshness; keeping the serial
     * with the point-quality owner prevents a detached controller latch from
     * outliving this lifecycle. */
    uint64_t saturatedPlanSerialValue = 0;
};

class BObolLodPointProxyEvidence {
private:
    friend class BObolLodAdmissionPlanner;

    static constexpr float maximumPixelThreshold(void)
    {
	return 64.0f;
    }

    static constexpr float narrowBracketRatio(void)
    {
	return 1.08f;
    }

    static bool atMaximumPixelThreshold(float threshold)
    {
	return threshold >= maximumPixelThreshold() - 0.01f;
    }

    bool terminalReplayRequired(float threshold) const
    {
	const float current = sanitize(threshold);
	return current > 1.01f &&
	    std::fabs(this->settledThreshold - current) > 0.01f;
    }

    void settleAtStructuralLimit(float threshold)
    {
	this->settledThreshold = sanitize(threshold);
    }

public:
    struct Decision {
	float threshold = 1.0f;
	bool changed = false;
	bool continueRelaxation = false;
    };

private:
    static bool applicable(size_t visibleOccurrenceCount)
    {
	/* This controller exists to bound the per-occurrence floor of a
	 * many-part scene.  One prominent mesh has no dispatch population to
	 * aggregate; increasing its point threshold cannot help until it destroys
	 * the very silhouette LoD is required to preserve. */
	return visibleOccurrenceCount > 1;
    }

    /* A camera revision invalidates the projected visibility census, not the
     * renderer capacity proof represented by an already active aggregate
     * point cut.  The replacement census is intentionally empty for a few
     * bounded source windows; treating that transient zero as a one-object
     * scene exposes every structural fallback box and makes pose/zoom changes
     * repeat cold-start calibration.  Retain the cut as the starting
     * presentation and let exact frames in the new view relax or tighten it.
     * Source/population domain changes reset the threshold independently. */
    static bool applicableAcrossCameraInvalidation(
	bool currentCensusApplicable, float retainedThreshold)
    {
	return currentCensusApplicable || sanitize(retainedThreshold) > 1.01f;
    }

    /* Source providers append structural records asynchronously.  That work
     * is precisely when completed frames need presentation-only aggregation:
     * excluding it leaves a large scene showing an ever-growing box forest
     * until discovery ends.  Keep this producer test separate from calibration
     * frame ownership, since a provider may both own future frames and remain
     * eligible for an immediate pressure correction. */
    static bool streamingPopulationWorkPending(bool submissionPending,
	bool resultsPending, bool servicePending, size_t providerPending)
    {
	return submissionPending || resultsPending || servicePending ||
	    providerPending > 0;
    }

    /* Point calibration cannot consume a frame which still contains an
     * unresolved structural fallback: its cheap timing describes boxes, not
     * the mesh/point population governed by the calibration bracket.  While
     * a provider, submission, or result publication is active, that producer
     * owns the successor frame at the changed threshold.  Once all such work
     * is quiet, however, retaining the calibration latch forms a closed wait:
     * stable calibration waits for a realized population while structural
     * repair is forbidden to realize it.  Yield the latch to the exact
     * structural frontier; the repair result supplies the next valid timing
     * sample without changing the proven aggregate threshold. */
    static bool calibrationYieldsToStructuralRepair(
	bool calibrationPending, bool exactOccurrenceClassification,
	size_t unresolvedStructuralCount, bool producerOwnsFutureFrame)
    {
	return calibrationPending && exactOccurrenceClassification &&
	    unresolvedStructuralCount > 0 && !producerOwnsFutureFrame;
    }

    static bool progressivePresentationUnconstrained(
	int progressiveCeiling, int maximumActiveProgressiveCut)
    {
	return progressiveCeiling < 0 || maximumActiveProgressiveCut < 0 ||
	    progressiveCeiling >= maximumActiveProgressiveCut;
    }

    /* The one-pixel presentation contract deliberately permits genuinely
     * subpixel occurrences to remain in the assembly's point batch.  They
     * are the bounded replacement for thousands of independent minimum-mesh
     * draws, not evidence that a temporary coarse population cut remains.
     * Recovery is complete once the exact classifier frame contains no
     * structural boxes and no renderer-wide PoP constraint.  A ceiling at or
     * above every active cut is inert and therefore does not make the frame
     * inexact. */
    static bool onePixelPresentationReady(bool haveCadAssemblies,
	bool exactFrame, bool exactOccurrenceClassification,
	size_t presentedSubpixelOccurrences,
	size_t presentedStructuralBoxes, int progressiveCeiling,
	int maximumActiveProgressiveCut,
	float pointProxyPixelThreshold, size_t managedRenderCost)
    {
	(void)presentedSubpixelOccurrences;
	return haveCadAssemblies && exactFrame &&
	    exactOccurrenceClassification &&
	    presentedStructuralBoxes == 0 &&
	    progressivePresentationUnconstrained(
		progressiveCeiling, maximumActiveProgressiveCut) &&
	    std::isfinite(pointProxyPixelThreshold) &&
	    pointProxyPixelThreshold <= 1.01f && managedRenderCost > 0;
    }

    /* A capacity sample measures the exact presentation named by its
     * occurrence-allocation certificate.  Unlike the terminal one-pixel
     * contract, that allocation may intentionally aggregate small
     * occurrences.  The presentation proof accounts for any renderer ceiling;
     * structural boxes may not stand in for the mesh work being timed. */
    static bool capacitySamplePopulationReady(bool haveCadAssemblies,
	bool exactFrame, bool exactOccurrenceClassification,
	size_t presentedStructuralBoxes, bool allocationPresentationRealized)
    {
	return haveCadAssemblies && exactFrame &&
	    exactOccurrenceClassification && presentedStructuralBoxes == 0 &&
	    allocationPresentationRealized;
    }

    static bool capacitySampleYieldsToStructuralFrontier(
	bool capacityPresentationPending, bool exactFrame,
	bool exactOccurrenceClassification, size_t presentedStructuralBoxes)
    {
	/* Capacity timing deliberately excludes structural placeholders.  An
	 * exact unchanged frame containing one can therefore never satisfy the
	 * pending sample.  Release that sample so point classification or
	 * structural repair can own the frontier instead of repainting it. */
	return capacityPresentationPending && exactFrame &&
	    exactOccurrenceClassification && presentedStructuralBoxes > 0;
    }

    /* Decide when a quiet deadline miss has reached the per-occurrence floor
     * rather than merely an expensive triangle prefix.  The renderer may
     * have attempted a rich prefix, but its one-step cost predictor can prove
     * that even cut zero exceeds the safe budget.  In that case another PoP
     * handoff cannot succeed; the reversible small-occurrence population cut
     * is the remaining presentation control. */
    static bool deadlineRequiresPopulationAggregation(
	    size_t activeRenderCost, size_t minimumRenderCost,
	    int presentedMaximum, int correctedCeiling,
	    size_t correctedRenderCostBudget)
    {
	return activeRenderCost <= minimumRenderCost ||
	    presentedMaximum <= 0 ||
	    (correctedCeiling == 0 && correctedRenderCostBudget > 0 &&
	     minimumRenderCost > correctedRenderCostBudget);
    }

    /* Decide whether an independent producer already owns the next frame
     * needed to measure a changed point/mesh classification.  A submission
     * cursor normally counts as such a producer, except while point
     * calibration itself blocks submitLodRequestsIfNeeded(): counting that
     * paused cursor creates a closed wait in which admission waits for the
     * frame and the frame waits for admission. */
    static bool producerOwnsCalibrationFrame(bool submissionPending,
	bool submissionPausedByCalibration, bool providerPending,
	bool servicePending, bool publicationAwaitingFrameRequest)
    {
	return BObolLodProducerPolicy::ownsFutureFrame(
	    submissionPending, submissionPausedByCalibration,
	    providerPending, servicePending,
	    publicationAwaitingFrameRequest);
    }

    /* Discovery needs an exact structural-size census before it can bound
     * first-wave provider work.  Stable calibration measures already
     * discovered occurrences and cannot block the source pass which creates
     * the CAD presentation it needs to observe. */
    static bool blocksSourceAdmission(bool discoveryCalibrationPending,
	bool stableCalibrationPending)
    {
	(void)stableCalibrationPending;
	return discoveryCalibrationPending;
    }

    /* A point-population mutation is a presentation transaction: its next
     * exact frame must classify the threshold which was actually requested.
     * In particular, a completed coarse structural seed may arm a one-pixel
     * static-quality trial.  Reseeding from the just-completed coarse frame
     * before that trial renders replaces its threshold and creates an
     * unbounded coarse/one-pixel request loop with no producer work. */
    static bool maySeedStructuralDistribution(
	bool discoveryCalibrationPending,
	bool stableCalibrationPending,
	bool triangleRecoveryPending)
    {
	return !discoveryCalibrationPending &&
	    !stableCalibrationPending && !triangleRecoveryPending;
    }

public:
    void reset(void)
    {
	this->unsafeThreshold = 0.0f;
	this->safeThreshold = 0.0f;
	this->settledThreshold = 0.0f;
    }

private:

    /* Use the renderer's exact projected-size census to bound the number of
     * structural occurrences which may start independent mesh-provider work.
     * This is a lower-bound cost proof, not a timing guess: every uncollapsed
     * occurrence has irreducible instance/dispatch work even at its minimum
     * PoP prefix.  Pick the finest power-of-two threshold which brings that
     * floor within the current first-wave allowance.  A completed frame then
     * supplies the ordinary safe-side timing witness and the same bracket can
     * relax toward pixel exactness without reloading discarded geometry. */
    Decision seedFromStructuralDistribution(float currentThreshold,
	const std::array<size_t, 7> &cumulativeCount,
	size_t visibleCount, size_t maximumUncollapsedCount)
    {
	Decision decision;
	const float current = sanitize(currentThreshold);
	decision.threshold = current;
	if (!visibleCount || !maximumUncollapsedCount)
	    return decision;

	float limit = 1.0f;
	float next = current;
	for (size_t bucket = 0; bucket < cumulativeCount.size();
		++bucket, limit *= 2.0f) {
	    const size_t collapsed = std::min(
		visibleCount, cumulativeCount[bucket]);
	    const size_t uncollapsed = visibleCount - collapsed;
	    if (uncollapsed <= maximumUncollapsedCount) {
		next = std::max(current, limit);
		break;
	    }
	    if (bucket + 1 == cumulativeCount.size())
		next = std::max(current, limit);
	}

	next = sanitize(next);
	if (next <= current + 0.01f)
	    return decision;
	/* Replaying an unchanged, exact structural census is an idempotent
	 * observation.  In particular, it must not erase the terminal witness
	 * established when the current cut could not afford its finer preload.
	 * Only a threshold transition changes the point/mesh population named by
	 * that witness. */
	this->settledThreshold = 0.0f;
	/* Every threshold below the selected census boundary leaves too many
	 * independently submitted occurrences for the first-wave allowance. */
	this->unsafeThreshold = std::max(this->unsafeThreshold, current);
	if (this->safeThreshold > 0.0f &&
		this->safeThreshold <= this->unsafeThreshold + 0.0001f)
	    this->safeThreshold = 0.0f;
	decision.threshold = next;
	decision.changed = true;
	decision.continueRelaxation = next > 1.01f;
	return decision;
    }

    Decision interrupted(float currentThreshold,
	uint64_t renderNanoseconds, float targetFps)
    {
	Decision decision;
	const float current = sanitize(currentThreshold);
	decision.threshold = current;
	this->unsafeThreshold = std::max(this->unsafeThreshold, current);
	/* A formerly safe threshold which now misses describes a larger/newer
	 * scene population.  Its old upper bracket is no longer evidence. */
	if (this->safeThreshold > 0.0f &&
	    this->safeThreshold <= this->unsafeThreshold + 0.0001f)
	    this->safeThreshold = 0.0f;

	float next = BObolLodQualityPolicy::pointProxyThreshold(
	    current, renderNanoseconds, targetFps);
	if (this->safeThreshold > this->unsafeThreshold) {
	    /* The timing correction and the safe/unsafe bracket are independent
	     * evidence about the same coarsening direction.  Honor the stronger
	     * correction while keeping the candidate inside the proven bracket.
	     * Falling back to geometric bisection alone can spend many expensive
	     * frames on thresholds which classify the same large occurrence
	     * population. */
	    if (this->safeThreshold / this->unsafeThreshold <=
		    narrowBracketRatio()) {
		next = this->safeThreshold;
	    } else {
		const float midpoint = static_cast<float>(std::sqrt(
		    static_cast<double>(this->safeThreshold) *
		    static_cast<double>(this->unsafeThreshold)));
		next = std::min(this->safeThreshold,
		    std::max(midpoint, next));
	    }
	}
	else if (next <= current + 0.01f && current < 64.0f) {
	    /* pointProxyThreshold() deliberately ignores completed-frame timing
	     * within five percent of its preferred cadence.  This entry point is
	     * different: the renderer has already crossed its hard endpoint
	     * deadline and aborted an incomplete frame after PoP detail reached its
	     * irreducible prefix.  Applying the completed-frame jitter deadband here
	     * leaves the point cut unchanged and requests the identical frame
	     * forever (a first-warm 50k OSMesa wire view reliably exposed it at
	     * roughly 101 ms against a 100 ms deadline).
	     *
	     * Advance the unsafe side by one bounded step.  The successor exact
	     * frame establishes the safe side and resumes the ordinary geometric
	     * bracket, so this does not turn small timing noise into an unbounded
	     * quality collapse. */
	    next = std::min(64.0f, current * 1.25f);
	}
	next = sanitize(next);
	decision.threshold = next;
	decision.changed = std::fabs(next - current) > 0.01f;
	decision.continueRelaxation = decision.changed && next > 1.01f;
	this->settledThreshold = decision.changed ? 0.0f : current;
	return decision;
    }

    /*
     * An exact structural-repair pass is ordered by screen-space importance.
     * If its finite scene allowance cannot replace every remaining box, the
     * current point cut is not a terminal presentation even when each frame
     * is individually fast.  Admit the important mesh wave, then aggregate
     * the unresolved small-object tail before attempting another pass.
     *
     * This is a second kind of unsafe witness: it proves cumulative coverage
     * cannot fit the scene budget, rather than proving one rendered frame
     * missed its deadline.  It therefore advances deterministically from the
     * unresolved population and uses the same safe/unsafe bracket as timing
     * calibration.  Once the bracket is narrow, return directly to its proven
     * safe side instead of flickering around the boundary.
     */
    Decision structuralCoverageBlocked(float currentThreshold,
	size_t unresolvedStructuralCount)
    {
	Decision decision;
	const float current = sanitize(currentThreshold);
	decision.threshold = current;
	if (!unresolvedStructuralCount)
	    return decision;
	this->settledThreshold = 0.0f;

	this->unsafeThreshold = std::max(this->unsafeThreshold, current);
	if (this->safeThreshold > 0.0f &&
	    this->safeThreshold <= this->unsafeThreshold + 0.0001f)
	    this->safeThreshold = 0.0f;

	float next = current;
	if (this->safeThreshold > this->unsafeThreshold) {
	    if (this->safeThreshold / this->unsafeThreshold <=
		    narrowBracketRatio())
		next = this->safeThreshold;
	    else
		next = static_cast<float>(std::sqrt(
		    static_cast<double>(this->safeThreshold) *
		    static_cast<double>(this->unsafeThreshold)));
	} else {
	    const float factor = unresolvedStructuralCount >= 4096 ? 2.0f :
		(unresolvedStructuralCount >= 512 ?
		 static_cast<float>(std::sqrt(2.0)) : 1.25f);
	    next = current * factor;
	}

	next = sanitize(next);
	decision.threshold = next;
	decision.changed = std::fabs(next - current) > 0.01f;
	decision.continueRelaxation = decision.changed && next > 1.01f;
	return decision;
    }

    Decision completed(float currentThreshold,
	uint64_t renderNanoseconds, float targetFps,
	bool reusableSample, bool allowCoarsening,
	size_t unresolvedStructuralCount = 0)
    {
	Decision decision;
	const float current = sanitize(currentThreshold);
	decision.threshold = current;
	if (!reusableSample || !renderNanoseconds || targetFps <= 0.0f)
	    return decision;

	const long double target = 1000000000.0L /
	    static_cast<long double>(targetFps);
	if (static_cast<long double>(renderNanoseconds) > target * 1.10L) {
	    if (!allowCoarsening) {
		/* A compact scene is governed by the independently bounded static
		 * deadline, not by the preferred quiet redraw cadence.  Do not
		 * contaminate its point bracket with an ordinary-frame miss before
		 * that hard trial has supplied negative capacity evidence. */
		this->settledThreshold = current;
		return decision;
	    }
	    return this->interrupted(
		current, renderNanoseconds, targetFps);
	}
	/* A cheap frame containing unresolved structural boxes is not evidence
	 * that this point cut is safe.  Let the importance-ordered repair pass
	 * either satisfy coverage or produce the explicit budget-blocked witness
	 * above.  In particular, do not lower the threshold and re-expose boxes
	 * before that proof exists. */
	if (unresolvedStructuralCount) {
	    this->settledThreshold = 0.0f;
	    return decision;
	}

	if (this->safeThreshold <= 0.0f || current < this->safeThreshold)
	    this->safeThreshold = current;
	/* A coarse presentation cut is never terminal merely because the
	 * controller's original pressure-probe latch was consumed by an
	 * intervening publication or plan-preparation frame.  Any unchanged,
	 * reusable frame with material headroom is fresh proof that a finer cut
	 * should be tried.  The remembered safe/unsafe bracket still makes the
	 * search bounded and prevents threshold chatter near the deadline. */
	if (current <= 1.01f ||
	    static_cast<long double>(renderNanoseconds) >= target * 0.80L) {
	    this->settledThreshold = current;
	    return decision;
	}

	float next = current;
	if (this->unsafeThreshold > 0.0f &&
	    this->safeThreshold > this->unsafeThreshold) {
	    /* Once the bracket is within eight percent, retaining the proven
	     * safe side avoids visible threshold chatter for negligible quality. */
	    if (this->safeThreshold / this->unsafeThreshold <=
		    narrowBracketRatio()) {
		this->settledThreshold = current;
		return decision;
	    }
	    next = static_cast<float>(std::sqrt(
		static_cast<double>(this->safeThreshold) *
		static_cast<double>(this->unsafeThreshold)));
	} else {
	    const long double ratio = std::sqrt(
		static_cast<long double>(renderNanoseconds) /
		(target * 0.80L));
	    const long double factor = std::max<long double>(0.75L,
		std::min<long double>(1.0L, ratio));
	    next = static_cast<float>(std::max<long double>(1.0L,
		static_cast<long double>(current) * factor));
	}

	next = sanitize(next);
	decision.threshold = next;
	decision.changed = std::fabs(next - current) > 0.01f;
	decision.continueRelaxation = decision.changed && next > 1.01f;
	this->settledThreshold = decision.changed ? 0.0f : current;
	return decision;
    }

    /* A transient exact frame at the one-pixel fidelity floor is already a
     * terminal quality witness.  Requiring an unchanged reusable replay can
     * strand an otherwise idle software view when faceplate publication or
     * retained-plan preparation keeps that replay transient.  Coarser point
     * cuts still require reusable evidence before they are accepted. */
    static bool requiresReusableConfirmation(float currentThreshold,
	size_t unresolvedStructuralCount = 0)
    {
	/* An unchanged box frame measures fallback presentation, not the realized
	 * CAD population governed by this bracket.  Let source repair consume the
	 * structural frontier; its resulting population supplies the next valid
	 * calibration frame. */
	return unresolvedStructuralCount == 0 &&
	    sanitize(currentThreshold) > 1.01f;
    }

    /*
     * Point aggregation and triangle-prefix admission are independent
     * quality dimensions.  A multi-pixel point cut is not, by itself,
     * evidence that the retained mesh population is too expensive: it may
     * be the finest sustainable way to represent thousands of genuinely
     * sub-pixel occurrences while leaving the triangle budget available for
     * large, visually prominent parts.  Compact PoP prefixes only when the
     * exact completed presentation is itself overloaded and those prefixes
     * can actually reduce its cost.
     *
     * Keeping this decision in the policy (rather than open-coding it in the
     * controller) makes the no-oscillation contract directly testable.
     */
    static bool shouldRecoverTriangleDetail(bool reducibleProgressiveDetail,
	bool stableSampleOverloaded, bool coarsePointCut,
	bool protectedQualityPolicyOwnsCuts)
    {
	(void)coarsePointCut;
	/* The protected-quality policy owns triangle redistribution both while a
	 * floor is active and after repeated hard misses reject that floor.  In
	 * the former state, compacting makes the next allocation restore the same
	 * cuts; in the latter, a separate recovery would bypass the soft minimax
	 * allocator and collapse the whole scene to minimum.  Either produces a
	 * sparse/cycling wire view.  Point aggregation may still calibrate, while
	 * the ordinary scene allocator alone chooses the retained PoP population. */
	return reducibleProgressiveDetail && stableSampleOverloaded &&
	    !protectedQualityPolicyOwnsCuts;
    }

    /* Triangle-prefix recovery and point/box aggregation are independent
     * controls.  Once triangle cost has been redistributed, return to one
     * pixel only when the exact assembly contains no structural fallback
     * population.  Otherwise a lower threshold exposes boxes rather than
     * richer geometry; retain the proven aggregate cut until mesh
     * publication replaces those records. */
    static float triangleRecoveryPointThreshold(float currentThreshold,
	bool exactStructuralPopulation, size_t structuralFallbackCount)
    {
	const float current = sanitize(currentThreshold);
	return exactStructuralPopulation && !structuralFallbackCount ?
	    1.0f : current;
    }

private:
    static float sanitize(float threshold)
    {
	return std::isfinite(threshold) ?
	    std::max(1.0f, std::min(64.0f, threshold)) : 1.0f;
    }

    float unsafeThreshold = 0.0f;
    float safeThreshold = 0.0f;
    /* A coarse point cut is terminal only after a reusable frame evaluated
     * that exact threshold.  Keeping this witness in the same reset domain as
     * the safe/unsafe bracket prevents source batching or cache warmth from
     * deciding whether the final relaxation sample happens to be requested. */
    float settledThreshold = 0.0f;
};

/* Point aggregation serves two different presentation contracts.  The
 * responsive bracket targets the ordinary quiet redraw cadence, while the
 * static bracket targets the longer interruptible endpoint deadline.  A
 * threshold which misses the former may be safely inside the latter, so
 * their safe/unsafe observations must never share one evidence domain. */
enum class BObolLodPointCalibrationGoal : uint8_t {
    RESPONSIVE,
    STATIC
};

/* Structural proxies with no useful point-aggregation successor are a
 * coverage floor: they can disappear only when their first mesh prefix is
 * admitted.  The point policy is exhausted either because the remaining
 * proxies exceed its maximum screen-space threshold or because the complete
 * occurrence population already fits the scene allowance.  Keep that
 * transition separate from ordinary quality calibration so a preferred-FPS
 * budget cannot repeatedly reopen the same exact frontier.  The controller
 * supplies a budget certified by the immediately preceding exact frame and
 * this policy grants at most one attempt for an unchanged
 * view/policy/projection population. */
class BObolLodStructuralAdmissionEvidence {
public:
    struct Inputs {
	uint64_t viewRevision = 0;
	uint64_t policyRevision = 0;
	uint64_t frontierDigest = 0;
	size_t unresolvedCount = 0;
	size_t activeCost = 0;
	size_t currentBudget = 0;
	size_t certifiedBudget = 0;
	bool exactProjection = false;
	bool pointPolicyExhausted = false;
    };

    struct Decision {
	bool ownsFrontier = false;
	bool requestAdmission = false;
	bool duplicateRejected = false;
	bool capacityLimited = false;
	size_t budget = 0;
    };

private:
    friend class BObolLodAdmissionPlanner;

    Decision planStructural(const Inputs &inputs)
    {
	Decision decision;
	decision.ownsFrontier = inputs.exactProjection &&
	    inputs.pointPolicyExhausted && inputs.unresolvedCount > 0;
	if (!decision.ownsFrontier)
	    return decision;

	const bool duplicate = this->validValue &&
	    this->viewRevisionValue == inputs.viewRevision &&
	    this->policyRevisionValue == inputs.policyRevision &&
	    this->frontierDigestValue == inputs.frontierDigest &&
	    this->unresolvedCountValue == inputs.unresolvedCount &&
	    this->activeCostValue == inputs.activeCost;
	if (duplicate) {
	    decision.duplicateRejected = true;
	    decision.capacityLimited = true;
	    decision.budget = this->attemptedBudgetValue;
	    return decision;
	}

	this->validValue = true;
	this->viewRevisionValue = inputs.viewRevision;
	this->policyRevisionValue = inputs.policyRevision;
	this->frontierDigestValue = inputs.frontierDigest;
	this->unresolvedCountValue = inputs.unresolvedCount;
	this->activeCostValue = inputs.activeCost;
	this->attemptedBudgetValue = inputs.certifiedBudget;
	decision.budget = inputs.certifiedBudget;
	/* currentBudget is a soft policy allowance and may describe a different
	 * retained allocation.  This transaction needs only positive capacity
	 * above the exact population which certified the structural frontier. */
	decision.requestAdmission = inputs.certifiedBudget > inputs.activeCost;
	decision.capacityLimited = !decision.requestAdmission;
	return decision;
    }

public:

    void reset(void)
    {
	this->validValue = false;
	this->viewRevisionValue = 0;
	this->policyRevisionValue = 0;
	this->frontierDigestValue = 0;
	this->unresolvedCountValue = 0;
	this->activeCostValue = 0;
	this->attemptedBudgetValue = 0;
    }

private:
    bool validValue = false;
    uint64_t viewRevisionValue = 0;
    uint64_t policyRevisionValue = 0;
    uint64_t frontierDigestValue = 0;
    size_t unresolvedCountValue = 0;
    size_t activeCostValue = 0;
    size_t attemptedBudgetValue = 0;
};

/*
 * All renderer-capacity facts consumed by admission have one owner.  The
 * component values remain deliberately small and allocation-free; grouping
 * them prevents a caller from committing one decision while accidentally
 * retaining stale evidence from another capacity dimension.
 */
class BObolLodAdmissionEvidence {
public:
    BObolLodAdmissionEvidence(void) = default;
    explicit BObolLodAdmissionEvidence(
	const BObolLodCapacityEvidence &value) : capacityValue(value) {}
    explicit BObolLodAdmissionEvidence(
	const BObolLodResourceEvidence &value) : resourcesValue(value) {}
    explicit BObolLodAdmissionEvidence(
	const BObolLodHeadroomEvidence &value) : headroomValue(value) {}
    explicit BObolLodAdmissionEvidence(
	const BObolLodPointProxyEvidence &value) :
	responsivePointProxyValue(value), staticPointProxyValue(value) {}
    explicit BObolLodAdmissionEvidence(
	const BObolLodStructuralAdmissionEvidence &value) :
	structuralValue(value) {}

    const BObolLodCapacityEvidence &capacity(void) const
    {
	return this->capacityValue;
    }
    const BObolLodResourceEvidence &resources(void) const
    {
	return this->resourcesValue;
    }
    const BObolLodHeadroomEvidence &headroom(void) const
    {
	return this->headroomValue;
    }
    const BObolLodPointProxyEvidence &pointProxy(void) const
    {
	return this->responsivePointProxyValue;
    }
    const BObolLodPointProxyEvidence &pointProxy(
	BObolLodPointCalibrationGoal goal) const
    {
	return goal == BObolLodPointCalibrationGoal::STATIC ?
	    this->staticPointProxyValue : this->responsivePointProxyValue;
    }
    const BObolLodStructuralAdmissionEvidence &structural(void) const
    {
	return this->structuralValue;
    }

private:
    friend class BObolLodAdmissionPlanner;

    BObolLodCapacityEvidence capacityValue;
    BObolLodResourceEvidence resourcesValue;
    BObolLodHeadroomEvidence headroomValue;
    BObolLodPointProxyEvidence responsivePointProxyValue;
    BObolLodPointProxyEvidence staticPointProxyValue;
    BObolLodStructuralAdmissionEvidence structuralValue;
};

/*
 * One admission calculation consumes an immutable evidence snapshot and
 * returns both its typed decision and complete successor evidence.  Planning
 * therefore has no hidden mutation: callers may compare, discard, or commit
 * the result atomically with the revision tuple which supplied its inputs.
 */
struct BObolLodAdmissionPlan {
    BObolLodAdmissionRevisionStamp revisionStamp;
    BObolLodCapacityEvidence::Decision capacityDecision;
    BObolLodResourceEvidence::Decision resourceDecision;
    BObolLodPointProxyEvidence::Decision pointProxyDecision;
    BObolLodStructuralAdmissionEvidence::Decision structuralDecision;
    BObolLodCapacityEvidence::CalibrationDecision calibrationDecision;
    BObolLodCapacityEvidence::CompletedFrameDecision completedFrameDecision;
    bool headroomAccepted = false;
    bool transitionChanged = false;
    size_t transitionValue = 0;
    BObolLodAdmissionEvidence nextEvidence;
    BObolLodAdmissionCursor nextCursor;

    bool certifiedFor(const BObolLodAdmissionRevisionStamp &stamp) const
    {
	return !this->revisionStamp.empty() &&
	    this->revisionStamp.same(stamp);
    }
};

static_assert(std::is_trivially_copyable<BObolLodAdmissionEvidence>::value,
    "admission evidence must remain allocation-free");
static_assert(std::is_trivially_copyable<BObolLodAdmissionPlan>::value,
    "admission plans must remain allocation-free");

/* Event-driven static quality is a separate presentation phase, not an
 * annotation which any ordinary quiet-frame failure may retire.  Keep its
 * two release-critical predicates independent of controller latch ordering so
 * unit tests (and the formal coordinator model) can enforce that boundary. */
class BObolLodAdmissionPlanner {
public:
    struct PointCalibrationProducerInputs {
	bool submissionPending = false;
	bool discoveryCalibrationPending = false;
	bool stableCalibrationPending = false;
	bool capacityAllocationPending = false;
	bool capacitySamplePending = false;
	bool stablePresentationAvailable = false;
	bool providerPending = false;
	bool servicePending = false;
	bool publicationAwaitingFrameRequest = false;
    };

    enum class EvidenceAction : uint8_t {
	RESET_CAPACITY = 0,
	CLEAR_CAPACITY_LIMIT,
	INVALIDATE_CAPACITY_MEASUREMENT,
	RESET_CAPACITY_OVERLOAD,
	CLEAR_RETAINED_QUALITY_FLOOR,
	CLEAR_RETAINED_RECOVERY_CEILING,
	CLEAR_DEADLINE_CAPACITY_CEILING,
	CANCEL_HEADROOM_RETRY,
	RESET_POINT_PROXY,
	RESET_STRUCTURAL_ADMISSION,
	RESET_RESOURCE_SERVICE,
	MARK_RESOURCE_RECOVERY_HANDLED
    };

    static BObolLodAdmissionPlan applyEvidenceAction(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, EvidenceAction action);

    static BObolLodAdmissionPlan requestCapacityRescan(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor);

    static BObolLodAdmissionPlan acceptStaticPresentationConstraint(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget);

    static BObolLodAdmissionPlan completeAppliedAllocation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodCapacityEvidence::CompletedAllocationInputs &inputs);

    static BObolLodAdmissionPlan requestRetainedRecovery(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget);

    static BObolLodAdmissionPlan requestRetainedNormalization(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget);

    static BObolLodAdmissionPlan requestRetainedReallocation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	bool preserveCurrentBudget = true);

    static BObolLodAdmissionPlan resumeCapacityCandidateAllocation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor);

    static BObolLodAdmissionPlan requestPresentationReconciliation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget);

    static BObolLodAdmissionPlan requestCoverageCompletion(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t activeCost,
	size_t certifiedBudget);

    static BObolLodAdmissionPlan recordDeadlineCapacityMiss(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodCapacityEvidence::DeadlineMissInputs &inputs);

    static BObolLodAdmissionPlan setRetainedQualityFloor(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget,
	uint64_t populationSignature, size_t activeCost,
	size_t minimumActiveCost);

    static BObolLodAdmissionPlan recordRetainedQualityFloorMiss(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor);

    static BObolLodAdmissionPlan recordRetainedQualityFloorMet(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, bool exactProtectedPopulation,
	uint64_t populationSignature, size_t presentedCost);

    static BObolLodAdmissionPlan confirmRetainedRecoveryPresentation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, bool onePixelReady);

    static BObolLodAdmissionPlan raiseCapacityBudget(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget);

    static BObolLodAdmissionPlan finishBlockedCapacityPass(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodCapacityEvidence::CalibrationInputs &inputs);

    static BObolLodAdmissionPlan completeCapacityCalibrationFrame(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor);

    static BObolLodAdmissionPlan completeCapacitySearchFrame(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodCapacityEvidence::CompletedFrameInputs &inputs);

    static BObolLodAdmissionPlan retireUnmeasurableCapacityFrame(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor);

    static BObolLodAdmissionPlan plan(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodAdmissionRevisionStamp &revisionStamp,
	const BObolLodCapacityEvidence::Inputs &inputs);

    static BObolLodAdmissionPlan planResourceObservation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, bool cpuPressure,
	bool gpuPressure, bool recoveryEnabled);

    static BObolLodAdmissionPlan planHeadroomRetry(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
	size_t activePopulationCost, size_t currentBudget);

    static BObolLodAdmissionPlan planHeadroomTransientDeferral(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
	size_t activePopulationCost);

    static BObolLodAdmissionPlan planHeadroomConsumption(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
	size_t currentBudget, size_t activePopulationCost,
	long double demonstratedBudget, uint64_t elapsedNanoseconds,
	uint64_t targetNanoseconds, bool reusableSample);

    static BObolLodAdmissionPlan recordHeadroomRejection(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
	size_t activePopulationCost, size_t attemptedBudget);

    static BObolLodAdmissionPlan planPointStructuralDistribution(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, float currentThreshold,
	const std::array<size_t, 7> &cumulativeCount, size_t visibleCount,
	size_t maximumUncollapsedCount);

    static BObolLodAdmissionPlan planPointInterrupted(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodPointCalibrationGoal goal, float currentThreshold,
	uint64_t renderNanoseconds, float targetFps);

    static BObolLodAdmissionPlan planPointStructuralCoverageBlocked(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, float currentThreshold,
	size_t unresolvedStructuralCount);

    static BObolLodAdmissionPlan planPointCompleted(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodPointCalibrationGoal goal, float currentThreshold,
	uint64_t renderNanoseconds, float targetFps, bool reusableSample,
	bool allowCoarsening,
	size_t unresolvedStructuralCount = 0);

    /* A compact, fully classified source profile earns one direct static
     * quality attempt before presentation-only small-part aggregation may
     * hide multi-pixel occurrences.  Large or incomplete profiles retain the
     * adaptive path immediately; a compact profile enters it only after the
     * hard static deadline supplies a negative capacity witness. */
    static bool adaptivePointAggregationAllowed(bool meshFirstSceneSafe,
	bool staticQualityCapacityRejected);

    static BObolLodAdmissionPlan settlePointAtStructuralLimit(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodPointCalibrationGoal goal, float currentThreshold);

    static bool pointRequiresReusableConfirmation(float currentThreshold,
	size_t unresolvedStructuralCount = 0);

    /* Coarsening hides additional occurrences and therefore cannot expose a
     * structural box; preserve the exact timing-selected threshold.  A finer
     * candidate can reveal meshes, so map it to the conservative lower
     * projection bucket whose complete frontier must be preloaded first. */
    static float pointStructuralPreloadThreshold(float currentThreshold,
	float candidateThreshold);

    /* A finer aggregate cut may require first meshes for every newly exposed
     * structural occurrence.  Reject a candidate whose complete projected
     * population exceeds the static presentation deadline.  Admission itself
     * is independently split into bounded first-wave batches. */
    static bool pointRelaxationFitsDeadline(
	size_t activeOccurrenceCount, size_t additionalOccurrenceCount,
	uint64_t measuredRenderNanoseconds,
	uint64_t presentationDeadlineNanoseconds);

    static bool pointAggregationApplicable(size_t visibleOccurrenceCount);

    static bool pointAggregationApplicableAcrossCameraInvalidation(
	bool currentCensusApplicable, float retainedThreshold);

    static bool pointStreamingPopulationWorkPending(bool submissionPending,
	bool resultsPending, bool servicePending, size_t providerPending);

    static bool pointCalibrationYieldsToStructuralRepair(
	bool calibrationPending, bool exactOccurrenceClassification,
	size_t unresolvedStructuralCount, bool producerOwnsFutureFrame);

    static bool onePixelPresentationReady(bool haveCadAssemblies,
	bool exactFrame, bool exactOccurrenceClassification,
	size_t presentedSubpixelOccurrences,
	size_t presentedStructuralBoxes, int progressiveCeiling,
	int maximumActiveProgressiveCut,
	float pointProxyPixelThreshold, size_t managedRenderCost);

    static bool capacitySamplePopulationReady(bool haveCadAssemblies,
	bool exactFrame, bool exactOccurrenceClassification,
	size_t presentedStructuralBoxes, bool allocationPresentationRealized);

    static bool capacitySampleYieldsToStructuralFrontier(
	bool capacityPresentationPending, bool exactFrame,
	bool exactOccurrenceClassification, size_t presentedStructuralBoxes);

    static bool pointDeadlineRequiresPopulationAggregation(
	size_t activeRenderCost, size_t minimumRenderCost,
	int presentedMaximum, int correctedCeiling,
	size_t correctedRenderCostBudget);

    static bool pointProducerOwnsCalibrationFrame(
	const PointCalibrationProducerInputs &inputs);

    /* One completed framebuffer may advance only one measurement owner.  A
     * bounded capacity search has already frozen the point/mesh population it
     * is measuring; consuming point calibration from the same frame would
     * reset that certificate before it can collect its finite sample set.
     * The still-pending point phase supplies the successor frame after the
     * capacity search reaches a terminal decision. */
    static bool pointCalibrationOwnsCompletedFrame(
	bool pointCalibrationPending, bool capacitySamplePending);

    static bool pointBlocksSourceAdmission(
	bool discoveryCalibrationPending, bool stableCalibrationPending);

    /* Missing meshes in an exact, quiet structural frontier are admitted only
     * by the structural-repair transaction.  Ordinary source passes may still
     * retarget resident meshes, but must not race point classification by
     * replacing whatever boxes happen to be exposed by an intermediate
     * threshold. */
    static bool structuralFallbackAdmissionDeferred(
	bool quietView, bool coverageComplete, bool exactFrame,
	size_t presentedStructuralBoxes, size_t terminalOccurrenceFailures,
	bool directPreviewAuthorized, bool structuralRepairActive);

    /* Terminal/offline refinement has no interactive structural-repair
     * successor.  Keep the complete structural-only decision in one tested
     * mode-domain predicate so it cannot create ownerless mesh debt. */
    static bool structuralCoverageOnly(
	bool terminalRefinement, bool coverageActive,
	bool demandRefreshActive, bool directPreviewAuthorized,
	bool fallbackAdmissionDeferred);

    /* A presentation-owned measurement freezes the occurrence population it
     * is measuring.  An unchanged submission cursor must not run, or count as
     * a future geometry producer, until that frame completes.  Source and
     * inventory invalidation are handled before callers apply this gate and
     * explicitly retire stale capacity evidence. */
    static bool presentationPausesSubmission(
	bool discoveryCalibrationPending, bool stableCalibrationPending,
	bool capacityAllocationPending, bool capacitySamplePending,
	bool stablePresentationAvailable);

    static bool maySeedPointStructuralDistribution(
	bool discoveryCalibrationPending, bool stableCalibrationPending,
	bool triangleRecoveryPending);

    static bool pointAtMaximumPixelThreshold(float threshold);

    /* Keep timing evidence and its deadline class in the same mapping.  A
     * transient controller phase may retire before terminal reduction, but
     * that must not change which target governs the evidence domain. */
    static float pointCalibrationTargetFps(
	BObolLodPointCalibrationGoal goal, float responsiveTargetFps,
	float staticTargetFps);

    static uint64_t pointCalibrationPresentationDeadline(
	BObolLodPointCalibrationGoal goal, uint64_t responsiveDeadline,
	uint64_t staticDeadline);

    static float maximumPointProxyPixelThreshold(void)
    {
	return BObolLodPointProxyEvidence::maximumPixelThreshold();
    }

    static bool pointTerminalReplayRequired(
	const BObolLodAdmissionEvidence &evidence,
	BObolLodPointCalibrationGoal goal, float threshold);

    static bool shouldRecoverTriangleDetail(
	bool reducibleProgressiveDetail, bool stableSampleOverloaded,
	bool coarsePointCut, bool protectedQualityOwnsCuts);

    static float triangleRecoveryPointThreshold(float currentThreshold,
	bool exactStructuralPopulation, size_t structuralFallbackCount);

    static BObolLodAdmissionPlan planStructural(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodStructuralAdmissionEvidence::Inputs &inputs);

    /* An exact structural-only frame is valid admission evidence even though
     * its CAD mesh cost is zero.  The current scene allowance supplies the
     * bounded first-mesh budget; requiring a positive mesh measurement here
     * leaves a cold all-box population with no transition capable of creating
     * that first mesh. */
    static bool structuralCapacityFrameApplicable(bool exactFrame,
	bool renderCostObserved, uint64_t renderNanoseconds,
	uint64_t presentationDeadlineNanoseconds);

    /* Compact scenes which passed the complete source-profile gate may use
     * the one-shot prominent-quality deadline for structural replacement.
     * Unbounded scenes retain the tighter finite-census deadline. */
    static uint64_t structuralRepairPresentationDeadline(
	bool directMeshScene, uint64_t stableDeadlineNanoseconds,
	uint64_t prominentDeadlineNanoseconds, uint64_t censusMultiplier);

    static size_t unaggregatableStructuralCount(
	const std::array<size_t, 7> &cumulativeCount, size_t visibleCount);

    /* Divide only proven marginal capacity.  A zero result means the exact
     * frame did not establish even one scheduling unit per unresolved
     * occurrence, so the caller must not arm a partial repair transaction. */
    static size_t structuralPerOccurrenceReservation(size_t admittedBudget,
	size_t activeCost, size_t unresolvedCount);

    /* noteDeadlineCapacityMiss() records an already strict, view-local upper
     * bound: it is the immediate predecessor of the candidate which missed
     * the hard deadline.  Marginal recovery has its own throughput safety
     * factor before reaching this helper.  Applying another factor to this
     * exact bound would discard untested capacity and can skip a useful
     * discrete PoP population. */
    static size_t capBudgetAtDeadlineCeiling(size_t estimatedBudget,
	size_t deadlineCapacityCeiling);

    static bool rejectAfterInterruptedFrame(bool interactive,
	bool cadDrawAttempted, bool preparationOnlyRetry,
	bool staticTrialActive);

    /* The one-pixel point trial and its enclosing static-quality trial own one
     * candidate framebuffer.  If that frame is rejected, both owners must
     * release it together.  Retiring an unrelated point calibration would
     * instead skip a required exact-classification frame. */
    static bool rejectedOnePixelTrialReleasesCalibration(bool trialRejected,
	bool calibrationPending, float pointProxyPixelThreshold);

    /* A deadline miss can close the static-quality phase only when it was an
     * explicit static-quality trial against a closed occurrence population.
     * An ordinary quiet-frame miss is capacity evidence for that exact
     * attempted population, but it must not retire a static phase which has
     * not started yet.  Doing so made a cold scene remember a transient
     * realization miss and permanently forgo the eventual event-driven 10 Hz
     * quality pass. */
    static bool terminalCapacityMiss(bool interactive,
	bool cadDrawAttempted, bool transientPresentationRetry,
	bool staticTrialActive, bool populationWorkPending);

    static bool ordinaryHeadroomAllowed(bool staticPhaseBlocked);

    /* A bounded static phase may improve an unaffordable atomic quality floor
     * incrementally.  The same marginal path remains available after either
     * the complete floor or a richer presentation has been rejected, but it
     * must never leak into ordinary quiet convergence. */
    static bool marginalStaticCapacityAllowed(bool staticPhaseActive,
	bool staticPhaseRejected, bool protectedFloorRejected);

    /* Retained occurrence cuts and the renderer's submitted cut are separate
     * quality dimensions.  A single progressive object can have every
     * requested suffix resident while a motion/quiet-cadence ceiling still
     * hides most of it.  Treat that reversible presentation gap as actionable
     * static quality debt; otherwise provider satisfaction falsely terminates
     * the event-driven quality pass after a retained-view interaction. */
    static bool actionableProgressiveQualityDebt(
	size_t activePayloadCount, size_t satisfiedPayloadCount,
	size_t memoryLimitedPayloadCount, int activeMaximumCut,
	int presentationCeiling);

    static bool progressivePresentationQualityDebt(int activeMaximumCut,
	int presentationCeiling);

    /* Translate a renderer-wide CAD cut into the retained allocator's scene
     * currency.  Backend submitted-work records may expand indexed vertices,
     * normals, and triangles differently (OSMesa flat shading is the important
     * case), so they are not a valid occurrence-allocation budget.  Preserve
     * the scene's non-CAD floor and replace only its retained CAD component with
     * the canonical cost table for the exact presented ceiling. */
    static size_t canonicalSceneCostAtCadCeiling(size_t activeSceneCost,
	size_t activeCadCost, size_t cadCeilingCost);

    /* A retained allocation certificate is the exact post-transaction
     * population.  The admission cursor deliberately freezes the pre-pass
     * scene cost for bounded traversal accounting, so using that older value
     * as the capacity candidate makes the completed frame look stale even
     * when no cut changed. */
    static size_t capacitySearchPresentedCost(bool currentAllocation,
	size_t frozenActiveCost, size_t selectedPresentationCost);

    static bool onePixelTrialRequired(float pointProxyPixelThreshold);

    /* A coarse point cut may hide either resident meshes or structural
     * fallback boxes.  Static headroom may expose the former under its hard
     * deadline, but exposing the latter is a visual regression, not quality
     * refinement.  Retain aggregation whenever structural records remain,
     * independent of whether a renderer-wide PoP ceiling is also installed. */
    static bool retainAggregatedPresentation(int progressiveCeiling,
	float pointProxyPixelThreshold, size_t activePayloadCount,
	bool structuralFallbackPopulation);

    /* Point aggregation is a convergence aid for a changing or streaming
	 * many-part scene, not the final fidelity policy of a fully realized
	 * event-driven framebuffer.  Once an exact, reusable point population has
	 * neither a producer nor structural fallback data behind it, test the
	 * normal one-pixel occurrence contract directly under the independently
	 * interruptible static deadline.  Structural points retain their proven
	 * cut until mesh publication replaces them; revealing their boxes is not a
	 * static-quality improvement.
     *
     * Keep this decision scalar and explicit: the controller owns state
     * transitions, while this predicate documents the proof required to skip
     * the destructive intermediate phase.  A rejected trial is terminal for
     * the current view/capacity epoch and must not be repeated. */
    static bool startOnePixelTrialFromSettledPointFrame(bool interactive,
	bool exactReusableFrame, bool producerWork, bool resourcePressure,
	bool structuralFallbackPopulation, bool staticPhaseBlocked,
	float pointProxyPixelThreshold,
	uint64_t renderNanoseconds, uint64_t preferredFrameNanoseconds,
	uint64_t staticDeadlineNanoseconds);

    /* A quiet exact one-pixel frame which misses the preferred streaming
     * cadence but fits the independently interruptible static deadline has
     * already proved the transition which a coarse point frame would need to
     * trial.  Enter the event-driven quality phase directly.  Sending this
     * frame through ordinary point/triangle overload recovery first coarsens
     * it to the preferred quiet budget; the later static pass then restores
     * the same population at its longer deadline, producing a visible cycle
     * balance/refine cycle on software renderers.
     *
     * A renderer-wide ceiling, live handoff, producer, structural fallback,
     * or resource edge makes the sample incomplete capacity evidence and must
     * retain the conservative path. */
    static bool acceptSettledOnePixelFrame(bool interactive,
	bool exactReusableFrame, bool producerWork, bool resourcePressure,
	bool structuralFallbackPopulation, bool handoffPending,
	bool staticPhaseBlocked,
	int progressiveCeiling, float pointProxyPixelThreshold,
	bool reducibleProgressiveDetail, uint64_t renderNanoseconds,
	uint64_t preferredFrameNanoseconds,
	uint64_t staticDeadlineNanoseconds);

    /* A global PoP ordinal is occurrence-local only when the scene contains
	* exactly one progressive payload.  In that case a one-cut completed-frame
	* probe bounds a discrete Lucy-like jump without distorting relative scene
	* importance.  A many-part scene must return -1 and spend static headroom
	* through the occurrence allocator; walking a scene-wide ordinal makes all
	* parts equally coarse and starves prominent shapes. */
    static int stagedProgressiveCeiling(int currentCeiling,
	    int activeMaximum, float pointProxyPixelThreshold,
	    size_t activePayloadCount);

    /* A single progressive occurrence can preserve its last exact cut while
     * a richer occurrence allocation is prepared behind it.  Without this
     * baseline ceiling, publication exposes the richest resident cut before
     * the completed-frame staircase can measure any neighboring cut.  A
     * global ordinal is not a valid relative-quality policy for multiple
     * occurrences, so those scenes continue to rely on the allocator. */
    static int staticProgressiveBaselineCeiling(int currentCeiling,
	    int activeMaximum, size_t activePayloadCount);

    static bool staticOrdinalOverscanApplicable(size_t activePayloadCount);

    static bool retainedPointAggregationApplicable(
	bool logicalPopulationApplicable, size_t pointProxyCandidateCount);

    /* A structural admission frame is itself a certificate that a preceding
     * exact census found a multi-occurrence population which the point cut can
     * change.  Preserve that certificate while the replacement frame is
     * incomplete: the live visibility census is intentionally inexact during
     * an interrupted traversal, and consulting it alone can disable the only
     * recovery control which makes the traversal cheaper. */
    static bool deadlinePointAggregationApplicable(
	bool retainedPopulationApplicable, bool structuralAdmissionPending,
	size_t logicalOccurrenceCount);

    /* A retained allocation counts only mesh-backed point candidates.  The
     * point threshold is inert only when an exact renderer census also proves
     * that no non-terminal structural fallback can use it. */
    static bool pointProxyThresholdInert(size_t retainedCandidateCount,
	bool exactStructuralProjection, size_t structuralOccurrenceCount,
	size_t terminalFailureCount);

    static bool structuralPointAggregationRequired(
	size_t logicalOccurrenceCount, size_t affordableOccurrenceCount);

    static bool pointPolicyExhaustedForStructuralFrontier(
	bool maximumThreshold, bool pointAggregationRequired,
	bool exactProjection, size_t unresolvedCount,
	size_t unaggregatableCount);

    static bool pointSuccessorRequiredForStructuralFrontier(
	bool exactProjection, size_t unresolvedCount,
	bool pointPolicyExhausted, bool anotherOwnerPending);

    static size_t structuralFirstWaveOccurrenceLimit(size_t sceneBudget);

private:
    static BObolLodAdmissionPlan beginPlan(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor);
};

/*
 * A provider result is published atomically.  If one task prepares both an
 * admitted presentation cut and a richer speculative resident suffix, the
 * useful visible improvement cannot be installed until all speculative work
 * has completed.  Large spatial meshes and wire preparation can turn that
 * coupling into several seconds of apparent post-motion level walking.
 *
 * Keep presentation and residency independent at the scheduling boundary:
 * an admitted visible advance is always delivered first.  A later pass may
 * continue prefetching.  When presentation cannot advance, resident-only
 * prefetch remains legal so memory and draw budgets cannot strand zoom data.
 */

#endif /* LIBBOBOL_LOD_ADMISSION_POLICY_PRIVATE_H */
