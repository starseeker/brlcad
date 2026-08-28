/* L O D _ P R E S E N T A T I O N _ P O L I C Y _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_PRESENTATION_POLICY_PRIVATE_H
#define LIBBOBOL_LOD_PRESENTATION_POLICY_PRIVATE_H

#include "common.h"
#include "BObol/BMeshLodCache.h"
#include "lod_quiet_successor_private.h"
#include "lod_revision_private.h"
#include "lod_view_snapshot_private.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

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
	/* Identity of the logical source/occurrence domain.  Resident prefix
	 * publication and view-local cut changes must not change this value: they
	 * are precisely the activity across which a pose-only presentation is
	 * restored. */
	uint64_t sceneDomainRevision = 0;
    };

    struct Snapshot {
	bool valid = false;
	BObolLodViewEpoch viewEpoch;
	float targetPixelError = 1.0f;
	int progressiveCeiling = -1;
	float progressiveNextFraction = 0.0f;
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
	float currentTargetPixelError = 1.0f;
	int currentProgressiveCeiling = -1;
	float currentProgressiveNextFraction = 0.0f;
	float currentPointProxyPixelThreshold = 1.0f;
    };

    /* Exact-view history is validated by its own revision-bound store before
     * it reaches the quiet successor.  Keeping the complete control vector in
     * one value prevents a later caller from restoring only an integer cut or
     * dropping the capacity proof associated with the completed frame. */
    struct QuietCertificate {
	bool valid = false;
	float targetPixelError = 1.0f;
	int progressiveCeiling = -1;
	float progressiveNextFraction = 0.0f;
	float pointProxyPixelThreshold = 1.0f;
	size_t presentedRenderCost = 0;
	size_t provenRenderCostCapacity = 0;
    };

    struct QuietInputs {
	RestoreInputs presentation;
	QuietCertificate exactView;
	/* Deterministic/offline convergence has no frame-capacity handoff.
	 * Interactions may finish after that mode is entered, so the quiet
	 * transition itself must reject restoration of a finite presentation
	 * ceiling rather than relying on the caller to cancel it afterward. */
	bool unboundedTerminal = false;
    };

    struct RestoreDecision {
	using Handoff = BObolLodQuietSuccessorReducer::Handoff;

	bool apply = false;
	bool restoredPriorStable = false;
	bool restoredProvenQuality = false;
	bool restoredExactView = false;
	Handoff handoff = Handoff::NONE;
	bool clearPresentationLimits = false;
	float targetPixelError = 1.0f;
	int progressiveCeiling = -1;
	float progressiveNextFraction = 0.0f;
	float pointProxyPixelThreshold = 1.0f;
	size_t provenRenderCostFloor = 0;
	size_t provenRenderCostCapacity = 0;

	bool needsHandoff(void) const
	{
	    return this->handoff != Handoff::NONE;
	}
    };

    struct CompletedPassInputs {
	bool completed = false;
	bool submissionPending = false;
	bool rescanAfterFrame = false;
	bool changedCut = false;
	/* A bounded source/delta scan is not a scene allocation proof.  Handoff
	 * may consume only a completed retained-allocation transaction covering
	 * the current occurrence population. */
	bool retainedAllocationCompleted = false;
	/* The completed transaction still names the active allocation plan and
	 * exact view/policy/point-classification epoch. */
	bool retainedAllocationCertified = false;
	/* The allocation transaction is atomic with respect to the owner-thread
	 * CAD revision, but a worker result may still enlarge the population
	 * immediately afterward.  Do not release a renderer-wide safety limit
	 * until no task/result/publication or coalesced resident-growth edge can
	 * invalidate the proof. */
	bool populationQuiescent = false;
	/* A complete scene-wide allocation represented the constrained
	 * framebuffer in occurrence-local cuts under its measured hard budget.
	 * This is stronger than a no-change observation: even when every cut was
	 * already suitable, the renderer-wide ceiling is now redundant. */
	bool presentationLimitsReconciled = false;
	bool retainedRefinementPending = false;
	bool retainedRefinementBudgetBlocked = false;
    };

    struct CompletedPassDecision {
	bool finishHandoff = false;
	/* The handoff observed a clean bounded/delta pass, but still needs the
	 * authoritative scene-wide retained allocation.  This is an immediate
	 * planning edge, not a presentation-frame request. */
	bool requestRetainedAllocation = false;
	/* The complete occurrence-local minimum still exceeds its certified
	 * presentation budget.  The controller must reduce independent small-part
	 * submissions and allocate again; a global PoP ceiling is not a terminal
	 * substitute for that operation. */
	bool requestLocalPresentationReduction = false;
	bool requestRetainedRescan = false;
	bool retireRetainedObservation = false;
    };

    enum class CompletedPassOwner : uint8_t {
	NONE = 0,
	CAPACITY,
	PRESENTATION_HANDOFF,
	ALLOCATION_HANDOFF
    };

    struct CompletedPassSelection {
	CompletedPassOwner owner = CompletedPassOwner::NONE;
	/* A current allocation certificate may prove that the next useful
	 * capacity sample must be taken without the temporary renderer ceiling.
	 * In that case the handoff owns the completed-pass snapshot and consumes
	 * its local cut/rescan annotations.  The still-pending capacity frame is
	 * evidence for the successor presentation, not a second owner. */
	bool consumePassAnnotations = false;

	bool capacityOwns(void) const
	{
	    return this->owner == CompletedPassOwner::CAPACITY;
	}

	bool handoffOwns(void) const
	{
	    return this->owner == CompletedPassOwner::PRESENTATION_HANDOFF ||
		this->owner == CompletedPassOwner::ALLOCATION_HANDOFF;
	}
    };

    /* A retained allocation is an authoritative replacement for a global
     * renderer ceiling only when the population it committed actually fits
     * the constrained frame's proven budget.  Merely completing the pass is
     * not sufficient: the minimum occurrence population may itself exceed
     * that budget, or a late resident suffix may leave the active cost above
     * it.  Keep this predicate here so the controller and policy tests share
     * one explicit handoff contract. */
    static bool presentationLimitsReconciled(bool allocationCompleted,
	    bool allocationCertified, size_t selectedPresentationCost,
	    size_t certifiedPresentationBudget,
	    size_t reconciliationBudget = 0);

    /* Capacity ownership begins with the first exact presentation, not with
     * construction of the search object.  Predict that pending sample when a
     * complete allocation is still hidden by a renderer ceiling so completed
     * pass ownership can remove the ceiling before calibration starts. */
    static bool capacitySamplePending(bool searchAwaitingSample,
	    bool allocationCurrent, bool allocationPresentationRealized);

    /* A complete occurrence allocation which exceeds its own certified
     * budget is not another capacity-search candidate.  It is the concrete
     * lower bound which the presentation-reconciliation policy must either
     * reduce (by aggregating independent small occurrences) or accept as a
     * bounded static constraint.  Claim that successor here so an inactive
     * handoff cannot leave the generic capacity fallback resubmitting the
     * identical allocation forever. */
    bool claimOverBudgetAllocation(bool allocationCurrent,
	    bool allocationCutsApplied, size_t selectedPresentationCost,
	    size_t certifiedPresentationBudget, size_t presentedRenderCost);

    /* Normalize the immutable completed-pass snapshot and choose its sole
     * effect-producing owner before capacity calibration or handoff mutates
     * evidence.  Callers execute only the selected owner. */
    CompletedPassSelection completedPassSelection(
	    const CompletedPassInputs &inputs, bool capacityEligible,
	    bool capacitySamplePending, int progressiveCeiling) const;

    /* A capacity-search sample must describe the exact occurrence-local
     * allocation named by its certificate.  A renderer-wide handoff ceiling
     * makes that sample inexact, while treating the pending sample as a
     * reason to retain the ceiling creates a circular wait.  Once the
     * allocation has supplied every other handoff proof, let ceiling removal
     * precede the sample.  The caller deliberately retains the capacity
     * search's frame latch, so the next ceiling-free presentation consumes
     * the same candidate rather than starting another allocation. */
    bool capacitySampleRequiresCeilingFreeHandoff(
	    const CompletedPassInputs &inputs, bool capacitySamplePending,
	    int progressiveCeiling) const;

    void capturePrior(float targetPixelError, int progressiveCeiling,
	float progressiveNextFraction, float pointProxyPixelThreshold,
	const Population &population,
	BObolLodViewEpoch viewEpoch);

    void noteStableQuality(float targetPixelError, int progressiveCeiling,
	float progressiveNextFraction, float pointProxyPixelThreshold,
	const Population &population,
	BObolLodViewEpoch viewEpoch, size_t presentedRenderCost);

    RestoreDecision beginQuiet(const QuietInputs &quietInputs);

    RestoreDecision beginQuiet(const RestoreInputs &inputs);

    CompletedPassDecision completePass(const CompletedPassInputs &inputs);

    void armHandoff(bool presentationRequired, size_t provenRenderCost = 0,
	size_t reconciliationBudgetLimit = 0);
    bool noteFramePresented(size_t provenRenderCost = 0,
	size_t reconciliationBudget = 0);
    void cancelHandoff(void);
    void viewInvalidated(void);
    void reset(void);

    bool handoffPending(void) const { return this->handoffActive(); }
    bool handoffPresentationPending(void) const
    {
	return this->presentationHandoffPending();
    }
    size_t handoffCostFloor(void) const
    {
	return this->handoffCostFloorValue;
    }
    size_t handoffReconciliationBudget(void) const
    {
	return this->handoffReconciliationBudgetValue;
    }
    bool priorSnapshotValid(void) const
    {
	return this->priorStableValue.valid;
    }
    bool provenQualityValid(void) const
    {
	return this->provenQualityValue.valid;
    }
    int provenQualityCeiling(void) const
    {
	return this->provenQualityValue.progressiveCeiling;
    }

private:
    enum class HandoffState : uint8_t {
	INACTIVE = 0,
	ALLOCATION_REQUIRED,
	PRESENTATION_REQUIRED
    };

    bool handoffActive(void) const
    {
	return this->handoffStateValue != HandoffState::INACTIVE;
    }

    bool presentationHandoffPending(void) const
    {
	return this->handoffStateValue ==
	    HandoffState::PRESENTATION_REQUIRED;
    }

    void clearHandoff(void);

    static float sanitizePixelError(float pixelError);

    static float sanitizeThreshold(float threshold);

    static float sanitizeFraction(int progressiveCeiling, float fraction);

    static Snapshot makeSnapshot(float targetPixelError,
	int progressiveCeiling, float progressiveNextFraction,
	float pointProxyPixelThreshold, const Population &population,
	BObolLodViewEpoch viewEpoch, size_t presentedRenderCost);

    static BObolLodQuietSuccessorReducer::Target targetFromSnapshot(
	const Snapshot &snapshot);

    static bool populationMatches(const Population &snapshot,
	const Population &current);

    Snapshot priorStableValue;
    Snapshot provenQualityValue;
    HandoffState handoffStateValue = HandoffState::INACTIVE;
    size_t handoffReconciliationBudgetValue = 0;
    size_t handoffReconciliationBudgetLimitValue = 0;
    size_t handoffCostFloorValue = 0;
};

static_assert(std::is_trivially_copyable<
    BObolLodPresentationPolicy>::value,
    "presentation handoff policy must remain an allocation-free value");

/*
 * Bounded history of exact camera states which have actually produced a
 * complete, deadline-safe CAD presentation.  This is deliberately separate
 * from BObolLodPresentationPolicy: that policy owns continuity within one
 * interaction epoch, while this value recognizes a later return to a proven
 * view.
 *
 * The history is a seed, never an admission authority.  Source inventory and
 * semantic policy changes advance/reset the controller-owned domain, resource
 * pressure vetoes recall, and the ordinary visibility, resident-memory, and
 * frame-deadline policies still validate the resulting frame.  A fixed LRU
 * avoids allocator work in the GUI/render path and exact snapshots avoid hash
 * collisions or approximate-pose surprises.
 */
class BObolLodViewQualityHistory {
public:
    /* A wheel/trackpad session can settle at more than eight useful scales
     * before returning to its starting view.  Thirty-two exact snapshots are
     * still only a few KiB, keep lookup bounded, and require no allocator in
     * the owner/render path. */
    static constexpr size_t capacityValue = 32;
    static constexpr size_t capacity(void) { return capacityValue; }

    struct Quality {
	float targetPixelError = 1.0f;
	int progressiveCeiling = -1;
	float progressiveNextFraction = 0.0f;
	float pointProxyPixelThreshold = 1.0f;
	/* Conservative scene-wide projected-pixel error bound produced by the
	 * retained importance allocator.  Unlike its allocation objective, this
	 * witness is not divided by the current requested target and can therefore
	 * compare completed frames.  Infinity means the frame did not carry that
	 * optional proof; it is not an error. */
	double maximumProjectedErrorPixels =
	    std::numeric_limits<double>::infinity();
	/* Work actually presented by the frame which supplied the fidelity
	 * controls above. */
	size_t presentedRenderCost = 0;
	/* Largest exact terminal workload proved deadline-safe at this same
	 * view.  This is capacity evidence, not a visual-quality ordering. */
	size_t provenRenderCostCapacity = 0;
    };

    struct RememberInputs {
	BObolLodViewSnapshot view;
	uint64_t domainRevision = 0;
	bool sceneAvailable = false;
	Quality quality;
	bool exactCompletedFrame = false;
	bool terminalPresentationComplete = false;
	/* Database/CSG and other registered geometry producers must have
	 * completed the exact source-policy epoch represented by this frame.
	 * Mesh-service idleness alone cannot prove that condition. */
	bool producersSettled = false;
	bool presentationDeadlineMet = false;
	bool resourcePressure = false;
    };

    struct RecallInputs {
	BObolLodViewSnapshot view;
	uint64_t domainRevision = 0;
	bool sceneAvailable = false;
	bool resourcePressure = false;
    };

    bool remember(const RememberInputs &inputs);
    bool recall(const RecallInputs &inputs, Quality &quality);
    void reset(void);

    size_t size(void) const;
    size_t rememberCount(void) const;
    size_t recallCount(void) const;

private:
    struct Entry {
	bool valid = false;
	BObolLodViewSnapshot view;
	uint64_t domainRevision = 0;
	Quality quality;
	size_t provenRenderCostCapacity = 0;
    };

    static bool validQuality(const Quality &quality);
    static Quality sanitizeQuality(const Quality &quality);
    static double progressiveRank(int ceiling, float nextFraction);
    static bool controlsDominate(const Quality &candidate,
	const Quality &current);
    static bool controlsEquivalent(const Quality &candidate,
	const Quality &current);
    static bool preferFidelityCandidate(const Quality &candidate,
	const Quality &current);
    void promote(size_t index);

    std::array<Entry, capacityValue> entries = {};
    size_t entryCount = 0;
    size_t rememberCountValue = 0;
    size_t recallCountValue = 0;
};

static_assert(std::is_trivially_copyable<
    BObolLodViewQualityHistory>::value,
    "exact-view quality history must remain an allocation-free value");

#endif /* LIBBOBOL_LOD_PRESENTATION_POLICY_PRIVATE_H */
