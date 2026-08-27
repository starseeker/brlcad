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
#include "BObol/BMeshLodCache.h"
#include "lod_capacity_search_private.h"
#include "lod_quiet_successor_private.h"
#include "lod_revision_private.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>
#include <unordered_map>
#include <vector>

/* Private provider parameters for the budgeted terminal-display path.  These
 * are request behavior, not persistent asset identity or public API. */
static constexpr const char *BOBOL_LOD_PREPARED_CAD_ONLY_PARAM =
    "display.prepared_cad_only";

static inline bool
bobol_lod_exact_scalar(double left, double right)
{
    if (!std::isfinite(left) || !std::isfinite(right))
	return false;
    uint64_t leftBits = 0;
    uint64_t rightBits = 0;
    std::memcpy(&leftBits, &left, sizeof(leftBits));
    std::memcpy(&rightBits, &right, sizeof(rightBits));
    if ((leftBits << 1) == 0)
	leftBits = 0;
    if ((rightBits << 1) == 0)
	rightBits = 0;
    return leftBits == rightBits;
}

static inline bool
bobol_lod_exact_scalar(float left, float right)
{
    if (!std::isfinite(left) || !std::isfinite(right))
	return false;
    uint32_t leftBits = 0;
    uint32_t rightBits = 0;
    std::memcpy(&leftBits, &left, sizeof(leftBits));
    std::memcpy(&rightBits, &right, sizeof(rightBits));
    if ((leftBits << 1) == 0)
	leftBits = 0;
    if ((rightBits << 1) == 0)
	rightBits = 0;
    return leftBits == rightBits;
}

/*
 * Exact camera and camera-scale identities used by the LoD coordinator.
 * Keep these allocation-free values beside the policies which consume them:
 * an exact-view quality proof must compare the same semantic inputs as a
 * camera epoch, rather than a lossy hash or an approximate pose.
 */
struct BObolLodViewSnapshot {
    bool same(const BObolLodViewSnapshot &other) const
    {
	/* Exact numeric equality intentionally treats +0 and -0 as the same
	 * camera while rejecting NaNs.  A byte comparison did the opposite and
	 * caused semantically identical restored views to miss history. */
	if (this->haveCamera != other.haveCamera ||
	    this->width != other.width || this->height != other.height ||
	    !bobol_lod_exact_scalar(this->size, other.size) ||
	    !bobol_lod_exact_scalar(this->lodScale, other.lodScale) ||
	    !bobol_lod_exact_scalar(this->curveScale, other.curveScale) ||
	    !bobol_lod_exact_scalar(this->pointScale, other.pointScale) ||
	    this->botThreshold != other.botThreshold)
	    return false;
	for (size_t i = 0; i < 16; ++i) {
	    if (!bobol_lod_exact_scalar(this->viewVolumeMatrix[i],
		    other.viewVolumeMatrix[i]))
		return false;
	}
	return true;
    }

    uint8_t haveCamera = 0;
    int32_t width = 0;
    int32_t height = 0;
    double size = 0.0;
    double lodScale = 0.0;
    double curveScale = 0.0;
    double pointScale = 0.0;
    uint32_t botThreshold = 0;
    float viewVolumeMatrix[16] = {};
};

struct BObolLodViewScaleSnapshot {
    bool same(const BObolLodViewScaleSnapshot &other) const
    {
	return this->haveCamera == other.haveCamera &&
	    this->width == other.width &&
	    this->height == other.height &&
	    bobol_lod_exact_scalar(this->size, other.size) &&
	    bobol_lod_exact_scalar(this->lodScale, other.lodScale) &&
	    bobol_lod_exact_scalar(this->curveScale, other.curveScale) &&
	    bobol_lod_exact_scalar(this->pointScale, other.pointScale) &&
	    this->botThreshold == other.botThreshold &&
	    this->viewportWidth == other.viewportWidth &&
	    this->viewportHeight == other.viewportHeight &&
	    this->cameraTypeKey == other.cameraTypeKey &&
	    bobol_lod_exact_scalar(this->aspectRatio, other.aspectRatio) &&
	    bobol_lod_exact_scalar(this->focalDistance,
		other.focalDistance) &&
	    bobol_lod_exact_scalar(this->projectionScale,
		other.projectionScale);
    }

    uint8_t haveCamera = 0;
    int32_t width = 0;
    int32_t height = 0;
    double size = 0.0;
    double lodScale = 0.0;
    double curveScale = 0.0;
    double pointScale = 0.0;
    uint32_t botThreshold = 0;
    int16_t viewportWidth = 0;
    int16_t viewportHeight = 0;
    uint64_t cameraTypeKey = 0;
    float aspectRatio = 0.0f;
    float focalDistance = 0.0f;
    float projectionScale = 0.0f;
};

static_assert(std::is_trivially_copyable<BObolLodViewSnapshot>::value,
    "view signatures must remain allocation-free values");
static_assert(std::is_trivially_copyable<BObolLodViewScaleSnapshot>::value,
    "view scale signatures must remain allocation-free values");

class BObolLodAdmissionPlanner;

/*
 * Completed-frame resource-pressure edge policy.  The controller supplies
 * measurements and executes the resulting compaction request, but it may not
 * independently reinterpret the revision/one-shot contract.  Keeping this
 * allocation-free value beside the phase machine makes persistent pressure a
 * provable memory-limited terminal state instead of a collection of latches.
 */
class BObolLodResourceEvidence {
public:
    struct Decision {
	bool changed = false;
	bool queueRecovery = false;
	bool pressureCleared = false;
	uint64_t revision = 0;
    };

private:
    friend class BObolLodAdmissionPlanner;

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

public:
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
 * Bounded proof for a late calibrated-headroom admission retry.  Normal probe
 * counts remain a short-horizon mechanism; this witness covers the ordering
 * where reusable timing becomes available only after those probes ended.
 *
 * The first retry is unique for a view/policy capacity epoch.  The active
 * render cost is retained only to prove that the explicitly requested replay
 * presented the same population; it is not the persistent epoch identity.
 * A minimax redistribution naturally changes that scalar even when the
 * occurrence population and renderer capacity are unchanged.  Treating it as
 * an identity lets a failed rich allocation recover to a slightly different
 * cost and immediately re-arm the same failed allowance forever.
 *
 * A later retry in the same capacity epoch is allowed only when calibration
 * has raised the finite render budget by more than 25 percent above the
 * largest budget already witnessed.  Source/camera/policy invalidation starts
 * a new epoch through the controller reset/revision contracts.  Consequently
 * repeated retries are logarithmically bounded without confusing allocation
 * output with a source-population revision.
 */
class BObolLodHeadroomEvidence {
private:
    friend class BObolLodAdmissionPlanner;

    bool armRetry(BObolLodViewEpoch viewEpoch,
	BObolLodPolicyEpoch policyEpoch, size_t activePopulationCost,
	size_t currentBudget)
    {
	if (currentBudget == SIZE_MAX)
	    return false;
	if ((this->pendingValue &&
		this->pendingViewRevision == viewEpoch.value() &&
		this->pendingPolicyRevision == policyEpoch.value() &&
		this->pendingPopulationCost == activePopulationCost))
	    return false;
	const bool witnessedCapacityEpoch =
	    this->viewRevision == viewEpoch.value() &&
	    this->policyRevision == policyEpoch.value();
	if (witnessedCapacityEpoch &&
	    !materialBudgetGrowth(currentBudget, this->budgetHighWater))
	    return false;
	this->pendingValue = true;
	this->pendingViewRevision = viewEpoch.value();
	this->pendingPolicyRevision = policyEpoch.value();
	this->pendingPopulationCost = activePopulationCost;
	this->pendingBudget = currentBudget;
	this->pendingTransientReplayCount = 0;
	return true;
    }

    bool deferTransientReplay(BObolLodViewEpoch viewEpoch,
	BObolLodPolicyEpoch policyEpoch, size_t activePopulationCost)
    {
	if (!this->pendingMatches(viewEpoch, policyEpoch,
		activePopulationCost) ||
	    this->pendingTransientReplayCount >=
		maxTransientReplayDeferrals())
	    return false;
	this->pendingTransientReplayCount++;
	return true;
    }

    static constexpr unsigned int maxTransientReplayDeferrals(void)
    {
	return 2;
    }

    bool consumeRetry(BObolLodViewEpoch viewEpoch,
	BObolLodPolicyEpoch policyEpoch, size_t currentBudget,
	size_t activePopulationCost,
	long double demonstratedBudget, uint64_t elapsedNanoseconds,
	uint64_t targetNanoseconds, bool reusableSample)
    {
	if (!this->pendingValue)
	    return false;
	const bool matchingWitness =
	    this->pendingViewRevision == viewEpoch.value() &&
	    this->pendingPolicyRevision == policyEpoch.value() &&
	    this->pendingPopulationCost == activePopulationCost;
	const size_t witnessedBudget = this->pendingBudget;
	this->pendingValue = false;
	this->pendingViewRevision = 0;
	this->pendingPolicyRevision = 0;
	this->pendingPopulationCost = SIZE_MAX;
	this->pendingBudget = 0;
	this->pendingTransientReplayCount = 0;
	if (!matchingWitness)
	    return false;
	/* Consuming the explicitly requested replay retires this exact witness
	 * even if the evidence is insufficient.  Otherwise an unchanged discrete
	 * population can request calibration frames forever. */
	const bool sameCapacityEpoch =
	    this->viewRevision == viewEpoch.value() &&
	    this->policyRevision == policyEpoch.value();
	this->viewRevision = viewEpoch.value();
	this->policyRevision = policyEpoch.value();
	this->populationCost = activePopulationCost;
	const size_t observedBudget = currentBudget > witnessedBudget ?
	    currentBudget : witnessedBudget;
	this->budgetHighWater = sameCapacityEpoch &&
	    this->budgetHighWater > observedBudget ?
	    this->budgetHighWater : observedBudget;
	/* There are two legitimate reasons to recompute the scene allocation:
	 * newly measured capacity exceeds its allowance, or the authoritative
	 * allowance is already much larger than the population stamped into the
	 * occurrences.  The latter happens when cold result timing changes which
	 * wave first calibrated throughput.  Requiring demonstratedBudget to beat
	 * the already-large allowance then strands stale minimum-cut allocation
	 * stamps even though the unchanged frame proves substantial headroom over
	 * what is actually drawn. */
	const long double population =
	    static_cast<long double>(activePopulationCost);
	const long double allowance =
	    static_cast<long double>(currentBudget);
	const bool measuredAllowanceGrowth =
	    demonstratedBudget * 20.0L > allowance * 21.0L;
	const bool staleUnderallocation =
	    population * 5.0L < allowance * 4.0L &&
	    demonstratedBudget * 5.0L > population * 6.0L;
	if (!reusableSample || !targetNanoseconds ||
	    (!measuredAllowanceGrowth && !staleUnderallocation) ||
	    static_cast<long double>(elapsedNanoseconds) * 5.0L >=
		static_cast<long double>(targetNanoseconds) * 6.0L)
	    return false;
	return true;
    }

public:
    void cancelRetry(void)
    {
	this->pendingValue = false;
	this->pendingViewRevision = 0;
	this->pendingPolicyRevision = 0;
	this->pendingPopulationCost = SIZE_MAX;
	this->pendingBudget = 0;
	this->pendingTransientReplayCount = 0;
    }

    /* A completed hard-deadline miss is a negative headroom proof, not a
     * cancellation.  Remember the attempted allowance before retiring its
     * pending replay so the same camera/policy capacity epoch cannot arm the
     * same rich frame after its constrained recovery frame.  The active cost
     * recorded below is diagnostic/replay state, not an epoch key: recovery's
     * occurrence-local redistribution is expected to change it.  Materially
     * larger measured capacity or a new view/policy epoch remains an explicit
     * progress witness accepted by armRetry(). */
private:
    void rejectRetry(BObolLodViewEpoch viewEpoch,
	BObolLodPolicyEpoch policyEpoch, size_t activePopulationCost,
	size_t attemptedBudget)
    {
	const bool sameCapacityEpoch =
	    this->viewRevision == viewEpoch.value() &&
	    this->policyRevision == policyEpoch.value();
	const size_t pendingAttempt = this->pendingValue &&
	    this->pendingViewRevision == viewEpoch.value() &&
	    this->pendingPolicyRevision == policyEpoch.value() ?
		this->pendingBudget : 0;
	this->viewRevision = viewEpoch.value();
	this->policyRevision = policyEpoch.value();
	this->populationCost = activePopulationCost;
	const size_t rejectedBudget = std::max(attemptedBudget,
	    pendingAttempt);
	this->budgetHighWater = sameCapacityEpoch ?
	    std::max(this->budgetHighWater, rejectedBudget) : rejectedBudget;
	this->cancelRetry();
    }

public:
    bool retryPending(void) const { return this->pendingValue; }
    bool pendingMatches(BObolLodViewEpoch viewEpoch,
	BObolLodPolicyEpoch policyEpoch, size_t activePopulationCost) const
    {
	return this->pendingValue &&
	    this->pendingViewRevision == viewEpoch.value() &&
	    this->pendingPolicyRevision == policyEpoch.value() &&
	    this->pendingPopulationCost == activePopulationCost;
    }

private:
    static bool materialBudgetGrowth(size_t currentBudget,
	size_t previousBudget)
    {
	if (currentBudget <= previousBudget)
	    return false;
	/* Subtract before comparing so this remains valid at SIZE_MAX. */
	return currentBudget - previousBudget > previousBudget / 4;
    }

    uint64_t viewRevision = 0;
    uint64_t policyRevision = 0;
    size_t populationCost = SIZE_MAX;
    size_t budgetHighWater = 0;
    bool pendingValue = false;
    uint64_t pendingViewRevision = 0;
    uint64_t pendingPolicyRevision = 0;
    size_t pendingPopulationCost = SIZE_MAX;
    size_t pendingBudget = 0;
    unsigned int pendingTransientReplayCount = 0;
};

/*
 * Current-view useful-presentation coverage proof.  Large compact sources are
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

    /* Coverage may defer quality until every visible occurrence has a useful
     * presentation.  Once that proof completes, only an unowned ordinary
     * coverage pass needs the generic deferred-quality successor.  A retained
     * allocation, point recovery, or presentation handoff already owns the
     * next transition and must reach its terminal reducer.  Starting another
     * submission pass in front of that owner prevents the owner from
     * consuming its complete-pass proof and creates a zero-cut rescan loop. */
    static bool needsDeferredQualitySuccessor(bool coverageCompleted,
	bool coverageMissing, bool retainedRefinementPending,
	bool qualitySuccessorOwned)
    {
	return coverageCompleted && !coverageMissing &&
	    retainedRefinementPending && !qualitySuccessorOwned;
    }

    void activate(bool invalidateCoverage)
    {
	this->activeValue = true;
	if (invalidateCoverage) {
	    this->coverageCompleteValue = false;
	    /* The counters belong to one exact view/inventory pass.  A camera or
	     * population epoch may begin while the preceding bounded pass is only
	     * partly consumed; carrying those observations into the replacement
	     * pass double-counts visible occurrences and can leave an otherwise
	     * complete scene permanently responsiveness-limited.  Quality-only
	     * resumptions use activate(false) and deliberately retain them. */
	    this->clearPassCounters();
	}
    }

    /* The semantic view policy owns whether useful source coverage is a
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

    /* An exact source delta may update the completed projected-visibility
     * census entry by entry while preserving the unchanged population's
     * coverage proof.  Publish that revised denominator without pretending a
     * new all-source coverage pass completed. */
    void setCompleteVisibleCount(size_t visibleCount)
    {
	this->completeVisibleCountValue = visibleCount;
	this->completeVisibleCountValidValue = true;
    }

    /* A closed structural-repair transaction is confirmed by the renderer,
     * not by another source scan.  Its exact useful-presentation population
     * supersedes a partially consumed retry census when it equals the
     * denominator published by the last complete visibility census. */
    bool confirmPresentedCoverage(size_t usefulPresentationCount)
    {
	if (!this->completeVisibleCountValidValue ||
	    usefulPresentationCount != this->completeVisibleCountValue)
	    return false;
	this->coverageCompleteValue = true;
	this->activeValue = false;
	this->clearPassCounters();
	return true;
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

/* Conservative source-profile gate for presenting meshes during initial
 * coverage instead of first completing a box-only pass.  This does not admit
 * bytes or drawing work: the service working-set governor, resident-memory
 * limit, scene allocator, and frame deadline remain authoritative.  It only
 * says the complete known scene is small enough that bounded mesh admission
 * may compete immediately with its structural proxies. */
class BObolLodSourceProfilePolicy {
public:
    static constexpr uint64_t meshFirstOccurrenceLimit = 4096;
    /* The profile counts every unique leaf asset, not only terminal-mesh
     * candidates.  Keep its cardinality bound aligned with the complete
     * occurrence bound; projection, mesh type, scene cost, service working
     * set, and resident memory independently gate each direct request. */
    static constexpr uint64_t meshFirstUniqueAssetLimit =
	meshFirstOccurrenceLimit;
    static constexpr uint64_t meshFirstEncodedByteLimit =
	64ULL * 1024ULL * 1024ULL;
    static constexpr uint64_t meshFirstLargestAssetByteLimit =
	16ULL * 1024ULL * 1024ULL;

    static bool safeMeshFirstPreview(bool profileComplete,
	uint64_t occurrenceCount, uint64_t uniqueAssetCount,
	uint64_t encodedSourceBytes, uint64_t largestAssetBytes,
	uint64_t meshRequestCount)
    {
	return profileComplete && occurrenceCount > 0 &&
	    occurrenceCount <= meshFirstOccurrenceLimit &&
	    uniqueAssetCount > 0 &&
	    uniqueAssetCount <= meshFirstUniqueAssetLimit &&
	    encodedSourceBytes > 0 &&
	    encodedSourceBytes <= meshFirstEncodedByteLimit &&
	    largestAssetBytes > 0 &&
	    largestAssetBytes <= meshFirstLargestAssetByteLimit &&
	    meshRequestCount > 0 && meshRequestCount <= occurrenceCount;
    }

    /* A safe source profile only enables this decision.  Projection and the
     * live aggregate draw allowance remain the authoritative per-occurrence
     * admission checks. */
    static bool admitTerminalMesh(bool safeScene, bool visible,
	bool subpixelPresentation, bool botSource, bool supportedDrawMode,
	size_t terminalCost, size_t remainingCost)
    {
	return safeScene && visible && !subpixelPresentation && botSource &&
	    supportedDrawMode && terminalCost > 0 &&
	    terminalCost <= remainingCost;
    }
};

/* Dense, view-owned projection evidence for compact CAD occurrences.
 *
 * A camera epoch projects every candidate box exactly once.  Subsequent
 * result, allocation, calibration, and handoff passes often change only the
 * LoD policy; repeating eight corner transforms plus hull construction for
 * every one of 50k--150k entries made those otherwise incremental passes the
 * dominant GUI-thread cost.  This cache stores only camera/geometric facts.
 * Target pixel error, cut selection, residency, budget, selection emphasis,
 * and visibility authoring remain live policy inputs on every visit.
 *
 * The controller owns this object and all access is on its presentation
 * thread.  Source routing identity plus compact-population epoch prevents
 * pointer/index recycling, while per-entry geometry and placement revisions
 * make exact edits miss without clearing unaffected evidence. */
class BObolLodProjectedDemandCache {
public:
    struct Evidence {
	bool visible = false;
	float pixelDiameter = 0.0f;
	float pixelArea = 0.0f;
	float pixelPerimeter = 0.0f;
	bool boundsContained = false;
	bool presentationValid = false;
	bool presentationVisible = false;
	bool presentationContained = false;
	float presentationPixelWidth = 0.0f;
	float presentationPixelHeight = 0.0f;
    };

    struct Entry {
	uint64_t generation = 0;
	uint64_t geometryRevision = 0;
	uint64_t placementRevision = 0;
	Evidence evidence;
    };

    struct Source {
	uint64_t routingId = 0;
	uint64_t populationEpoch = 0;
	uint64_t viewRevision = 0;
	uint64_t generation = 1;
	std::vector<Entry> entries;
    };

    void clear(void)
    {
	this->sources.clear();
    }

    Source *bind(uint64_t routingId, uint64_t populationEpoch,
	uint64_t viewRevision, size_t entryCount)
    {
	if (!routingId || !populationEpoch || !viewRevision)
	    return NULL;
	Source *bound = NULL;
	for (Source &source : this->sources) {
	    if (source.routingId == routingId) {
		bound = &source;
		break;
	    }
	}
	if (!bound) {
	    this->sources.push_back(Source());
	    bound = &this->sources.back();
	    bound->routingId = routingId;
	}
	if (bound->populationEpoch != populationEpoch ||
	    bound->viewRevision != viewRevision) {
	    bound->populationEpoch = populationEpoch;
	    bound->viewRevision = viewRevision;
	    if (++bound->generation == 0) {
		bound->generation = 1;
		bound->entries.clear();
	    }
	}
	if (bound->entries.size() < entryCount)
	    bound->entries.resize(entryCount);
	return bound;
    }

    static bool lookup(const Source *source, size_t entryIndex,
	uint64_t geometryRevision, uint64_t placementRevision,
	Evidence &evidence)
    {
	if (!source || entryIndex >= source->entries.size())
	    return false;
	const Entry &entry = source->entries[entryIndex];
	if (entry.generation != source->generation ||
	    entry.geometryRevision != geometryRevision ||
	    entry.placementRevision != placementRevision)
	    return false;
	evidence = entry.evidence;
	return true;
    }

    static void store(Source *source, size_t entryIndex,
	uint64_t geometryRevision, uint64_t placementRevision,
	const Evidence &evidence)
    {
	if (!source || entryIndex >= source->entries.size())
	    return;
	Entry &entry = source->entries[entryIndex];
	entry.generation = source->generation;
	entry.geometryRevision = geometryRevision;
	entry.placementRevision = placementRevision;
	entry.evidence = evidence;
    }

private:
    std::vector<Source> sources;
};

/* Projected mesh-target census for exact source deltas.
 *
 * A complete camera pass establishes one dense visibility bit per compact
 * source entry.  Visibility, selection, and in-place edit transactions then
 * revise only their changed indices.  This is the bridge between a durable
 * all-scene convergence denominator and O(delta) interactive scene updates:
 * no source or scene object is owned or dereferenced here. */
class BObolLodVisibilityCensus {
public:
    /* Fixed-width, object-lifetime routing identity.  A raw node address can
     * be recycled without an observable empty scene, causing a replacement
     * source to inherit the retired source's visibility bits. */
    using SourceKey = uint64_t;

    /* Reconcile the census identity domain with the controller's current LoD
     * source contracts.  Scene synchronization may replace one draw-mode
     * source with another without exposing an intermediate empty scene or a
     * separately observable source-change callback.  Marking and sweeping at
     * the already O(source-count) signature boundary prevents retired source
     * totals from contributing to the current-view denominator. */
    void beginSourceSetUpdate(void)
    {
	for (auto &entry : this->sources)
	    entry.second.retained = false;
    }

    bool retainSource(SourceKey source)
    {
	Source *entry = this->find(source);
	if (entry) {
	    entry->retained = true;
	    return false;
	}
	Source &inserted = this->sources[source];
	inserted.retained = true;
	return true;
    }

    bool endSourceSetUpdate(void)
    {
	const size_t priorSize = this->sources.size();
	for (auto entry = this->sources.begin(); entry != this->sources.end();) {
	    if (!entry->second.retained)
		entry = this->sources.erase(entry);
	    else
		++entry;
	}
	return this->sources.size() != priorSize;
    }

    void clear(void)
    {
	this->sources.clear();
    }

    void begin(SourceKey source, size_t entryCount)
    {
	Source &entry = this->source(source);
	entry.visibleCount = 0;
	entry.entryVisible.assign(entryCount, 0);
	entry.complete = false;
    }

    void observe(SourceKey source, size_t entryCount, size_t entryIndex,
	bool visible)
    {
	Source &entry = this->source(source);
	if (entry.entryVisible.size() < entryCount)
	    entry.entryVisible.resize(entryCount, 0);
	if (entryIndex >= entry.entryVisible.size())
	    entry.entryVisible.resize(entryIndex + 1, 0);
	const unsigned char next = visible ? 1 : 0;
	unsigned char &prior = entry.entryVisible[entryIndex];
	if (prior == next)
	    return;
	if (next) {
	    if (entry.visibleCount < SIZE_MAX)
		entry.visibleCount++;
	} else if (entry.visibleCount > 0) {
	    entry.visibleCount--;
	}
	prior = next;
    }

    void finish(SourceKey source)
    {
	Source &entry = this->source(source);
	entry.complete = true;
    }

    void setCount(SourceKey source, size_t visibleCount)
    {
	this->source(source).visibleCount = visibleCount;
    }

    bool complete(SourceKey source) const
    {
	const Source *entry = this->find(source);
	return entry && entry->complete;
    }

    bool visible(SourceKey source, size_t entryIndex) const
    {
	const Source *entry = this->find(source);
	return entry && entry->complete &&
	    entryIndex < entry->entryVisible.size() &&
	    entry->entryVisible[entryIndex] != 0;
    }

    size_t sourceCount(SourceKey source) const
    {
	const Source *entry = this->find(source);
	return entry ? entry->visibleCount : 0;
    }

    size_t total(void) const
    {
	size_t total = 0;
	for (const auto &entry : this->sources)
	    total = entry.second.visibleCount > SIZE_MAX - total ?
		SIZE_MAX : total + entry.second.visibleCount;
	return total;
    }

    size_t sourceEntryCount(void) const
    {
	return this->sources.size();
    }

private:
    struct Source {
	size_t visibleCount = 0;
	std::vector<unsigned char> entryVisible;
	bool complete = false;
	bool retained = true;
    };

    Source *find(SourceKey key)
    {
	auto entry = this->sources.find(key);
	return entry != this->sources.end() ? &entry->second : NULL;
    }

    const Source *find(SourceKey key) const
    {
	auto entry = this->sources.find(key);
	return entry != this->sources.end() ? &entry->second : NULL;
    }

    Source &source(SourceKey key)
    {
	return this->sources[key];
    }

    std::unordered_map<SourceKey, Source> sources;
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
	CONVERROR,
	PREPARING
    };

    enum class Outcome : uint8_t {
	ACTIVE = 0,
	READY,
	CONSTRAINED,
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
	size_t presentedSubpixelOccurrenceCount = 0;
	size_t presentedStructuralBoxCount = 0;
	size_t terminalOccurrenceFailureCount = 0;
	size_t memoryLimitedPayloadCount = 0;
	size_t pendingTasks = 0;
	size_t inFlight = 0;
	size_t queuedCacheWrites = 0;
	size_t stableResidentMeshBytes = 0;
	size_t residentMeshLimitBytes = SIZE_MAX;
	size_t gpuTrackedBufferBytes = 0;
	unsigned int failedSourceCount = 0;
	bool structuralDiscovery = false;
	/* TRUE only after the bounded projection scan has visited the complete
	 * current source population.  Before then visibleTargetCount is a prefix,
	 * not a safe progress denominator. */
	bool visibilityCensusComplete = false;
	/* An append-only source or PoP producer has not closed its immutable
	 * population.  Published prefixes remain useful, but calibration against
	 * the partial population is not terminal work. */
	bool sourcePreparationPending = false;
	bool submissionPending = false;
	bool resultPending = false;
	bool publicationPending = false;
	bool calibrationPending = false;
	bool controlPending = false;
	bool interactive = false;
	bool compactionPending = false;
	bool progressiveWorkPending = false;
	bool gpuMemoryPressure = false;
	bool stableBudgetLimited = false;
	/* The retained occurrence population may be richer than the completed
	 * framebuffer when a renderer-only PoP ceiling or point-proxy threshold
	 * protects the presentation deadline.  This is a valid terminal view, but
	 * it is still explicitly performance limited and must be reported as such
	 * to progress/HUD consumers. */
	bool presentationLimited = false;
    };

    struct Decision {
	Phase phase = Phase::IDLE;
	Outcome outcome = Outcome::READY;
	float fraction = 1.0f;
	bool terminal = true;
	bool terminalError = false;
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
	/* Once a source reports a terminal failure, its missing inventory has no
	 * enabled discovery edge.  Other sources remain active through their
	 * explicit discovery, preparation, and task witnesses below; treating the
	 * failed source's expected count as pending work would manufacture an
	 * ERROR state which can never become terminal. */
	const bool structuralPending = inputs.failedSourceCount == 0 &&
	    inputs.expectedLeafCount > inputs.availableLeafCount;
	const size_t unresolvedStructuralBoxes =
	    inputs.presentedStructuralBoxCount >
		inputs.terminalOccurrenceFailureCount ?
	    inputs.presentedStructuralBoxCount -
		inputs.terminalOccurrenceFailureCount : 0;
	const bool foregroundPending =
	    structuralPending || inputs.structuralDiscovery ||
	    inputs.sourcePreparationPending ||
	    inputs.submissionPending ||
	    inputs.resultPending || inputs.publicationPending ||
	    inputs.calibrationPending || inputs.controlPending ||
	    unresolvedStructuralBoxes > 0;
	const size_t representedPayloads = saturatingAdd(
	    inputs.activePayloadCount, inputs.presentedSubpixelOccurrenceCount);
	const size_t satisfiedPayloads = saturatingAdd(
	    inputs.satisfiedPayloadCount,
	    inputs.presentedSubpixelOccurrenceCount);
	/* A serialized PoP builder can keep one worker occupied persisting its
	 * reusable hierarchy after the current view is already completely
	 * represented.  That work must not make the HUD claim the visible scene is
	 * still refining.  Conversely, retain task work as foreground until every
	 * visible target has both a presentation and a current-quality witness. */
	/* A completed visibility census may legitimately find an empty view after
	 * a subpath erase or frustum change.  Its retained framebuffer is already
	 * complete; an asset worker left alive for cache reuse is background work,
	 * just as it is when all nonempty visible targets are satisfied.  Requiring
	 * a positive target count made that worker hold an empty view in REFINING. */
	const bool emptyViewResolved = inputs.visibilityCensusComplete &&
	    inputs.visibleTargetCount == 0;
	const bool populatedViewResolved = inputs.visibleTargetCount > 0 &&
	    representedPayloads >= inputs.visibleTargetCount &&
	    satisfiedPayloads >= inputs.visibleTargetCount;
	const bool viewPresentationResolved =
	    !foregroundPending && !inputs.interactive &&
	    (emptyViewResolved || populatedViewResolved);
	const bool backgroundTaskWork = viewPresentationResolved &&
	    (inputs.pendingTasks > 0 || inputs.inFlight > 0);
	decision.visualPending = foregroundPending ||
	    (!backgroundTaskWork &&
	     (inputs.pendingTasks > 0 || inputs.inFlight > 0));
	decision.terminal =
	    !decision.visualPending && !inputs.interactive;
	const bool hasTerminalError = inputs.failedSourceCount > 0 ||
	    inputs.terminalOccurrenceFailureCount > 0;
	decision.terminalError = decision.terminal && hasTerminalError;
	decision.viewReady = decision.terminal && !hasTerminalError;
	decision.backgroundPending =
	    backgroundTaskWork || inputs.queuedCacheWrites > 0 ||
	    inputs.compactionPending ||
	    inputs.structuralDiscovery ||
	    (decision.viewReady && inputs.progressiveWorkPending);
	decision.memoryLimited =
	    inputs.memoryLimitedPayloadCount > 0 ||
	    inputs.gpuMemoryPressure ||
	    (inputs.residentMeshLimitBytes != SIZE_MAX &&
	     inputs.stableResidentMeshBytes > inputs.residentMeshLimitBytes);
	decision.performanceLimited = decision.viewReady &&
	    (inputs.stableBudgetLimited || inputs.presentationLimited ||
	     decision.memoryLimited ||
	     (inputs.visibleTargetCount > 0 &&
	      (saturatingAdd(inputs.activePayloadCount,
		   inputs.presentedSubpixelOccurrenceCount) <
		   inputs.visibleTargetCount ||
	       saturatingAdd(inputs.satisfiedPayloadCount,
		   inputs.presentedSubpixelOccurrenceCount) <
		   inputs.visibleTargetCount)));
	decision.outcome = hasTerminalError ? Outcome::ERROR :
	    !decision.terminal ? Outcome::ACTIVE :
	    decision.performanceLimited ? Outcome::CONSTRAINED : Outcome::READY;
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
	} else if (inputs.sourcePreparationPending) {
	    decision.phase = Phase::PREPARING;
	    /* The visibility census may itself still be streaming while source
	     * preparation is active.  Using that partial census as the denominator
	     * made the fraction jump to 75 percent as soon as the first batch was
	     * represented, regardless of how much of a large model remained.  Use
	     * the known structural population until the producer closes; the exact
	     * visible denominator takes over in the refining phase. */
	    const size_t preparationTarget = std::max(
		inputs.visibleTargetCount, inputs.availableLeafCount);
	    if (preparationTarget == 0) {
		decision.fraction = 0.40f;
	    } else {
		const size_t prepared = std::min(
		    saturatingAdd(inputs.activePayloadCount,
			inputs.presentedSubpixelOccurrenceCount),
		    preparationTarget);
		const long double coverage =
		    static_cast<long double>(prepared) /
		    static_cast<long double>(preparationTarget);
		/* Reserve the final quarter for coherent presentation, allocation,
		 * and calibration after the producer closes. */
		decision.fraction = static_cast<float>(
		    0.40L + std::min<long double>(0.35L, 0.35L * coverage));
	    }
	} else if (inputs.interactive) {
	    decision.phase = Phase::INTERACTIVE;
	    if (inputs.visibleTargetCount == 0) {
		decision.fraction = 0.40f;
	    } else {
		const long double quality =
		    static_cast<long double>(std::min(
			saturatingAdd(inputs.satisfiedPayloadCount,
			    inputs.presentedSubpixelOccurrenceCount),
			inputs.visibleTargetCount)) /
		    static_cast<long double>(inputs.visibleTargetCount);
		decision.fraction = static_cast<float>(
		    0.40L + std::min<long double>(0.45L, 0.45L * quality));
	    }
	} else if (inputs.calibrationPending) {
	    decision.phase = Phase::CALIBRATING;
	    decision.fraction = visualWorkFraction(inputs);
	} else if (decision.visualPending) {
	    decision.phase = Phase::REFINING;
	    decision.fraction = visualWorkFraction(inputs);
	} else if (decision.backgroundPending) {
	    decision.phase = Phase::BACKGROUND;
	    decision.fraction = 1.0f;
	} else {
	    decision.phase = Phase::IDLE;
	    decision.fraction = 1.0f;
	}

	if (hasTerminalError)
	    decision.phase = Phase::CONVERROR;
	decision.fraction = std::max(0.0f,
	    std::min(1.0f, decision.fraction));
	if (this->fractionViewEpoch != inputs.viewEpoch ||
	    this->fractionPolicyEpoch != inputs.policyEpoch) {
	    this->fractionViewEpoch = inputs.viewEpoch;
	    this->fractionPolicyEpoch = inputs.policyEpoch;
	    this->fractionFloorValue = 0.0f;
	}
	if (decision.phase != Phase::CONVERROR) {
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
    static float visualWorkFraction(const Inputs &inputs)
    {
	size_t target = std::max(
	    inputs.visibleTargetCount, inputs.activePayloadCount);
	if (!inputs.visibilityCensusComplete)
	    target = std::max(target, inputs.availableLeafCount);
	if (!target) {
	    target = saturatingAdd(inputs.pendingTasks, inputs.inFlight);
	    target = std::max<size_t>(1, target);
	}
	const size_t represented = std::min(
	    saturatingAdd(inputs.activePayloadCount,
		inputs.presentedSubpixelOccurrenceCount), target);
	const long double representationCoverage =
	    static_cast<long double>(represented) /
	    static_cast<long double>(target);
	const size_t unresolvedStructuralBoxes =
	    inputs.presentedStructuralBoxCount >
		inputs.terminalOccurrenceFailureCount ?
	    inputs.presentedStructuralBoxCount -
		inputs.terminalOccurrenceFailureCount : 0;
	const size_t richTarget = saturatingAdd(
	    inputs.activePayloadCount, unresolvedStructuralBoxes);
	long double richCoverage = 0.0L;
	if (richTarget > 0) {
	    richCoverage =
		static_cast<long double>(std::min(
		    inputs.satisfiedPayloadCount, richTarget)) /
		static_cast<long double>(richTarget);
	} else if (inputs.visibilityCensusComplete && represented == target) {
	    /* An exact all-subpixel view has no rich tail to resolve. */
	    richCoverage = 1.0L;
	}
	/* Structural discovery owns the first 40 percent and initial useful
	 * representations the next 35.  The final visible-mesh tail is weighted
	 * separately: classifying one terminal subpixel proxy is much cheaper than
	 * loading, publishing, and validating one mesh, so treating both as equal
	 * objects made a 150k-part scene jump directly to 99 percent with hundreds
	 * of expensive mesh replacements still outstanding.  Gate rich progress
	 * by representation coverage so a partial visibility census cannot claim a
	 * nearly complete tail.  The last percent remains the coherent-frame
	 * handoff obligation; only a proven ready view reports 100 percent. */
	return static_cast<float>(
	    0.40L +
	    std::min<long double>(0.35L,
		0.35L * representationCoverage) +
	    std::min<long double>(0.24L,
		0.24L * representationCoverage * richCoverage));
    }

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
class BObolLodCapacityEvidence;

/*
 * A coverage or publication pass may need exactly one coherent framebuffer
 * before planning resumes.  That ordering edge is not a renderer-capacity
 * sample and must not share the bounded search's candidate counters.  A
 * quality-blocked ordinary pass additionally transfers its successor to the
 * complete importance allocator, which is the only path allowed to start a
 * numeric capacity search.
 */
class BObolLodCapacityFrameBarrier {
public:
    enum class Successor : uint8_t {
	NONE = 0,
	REPLAN,
	REALLOCATE
    };

    void requestReplan(void)
    {
	if (this->successorValue == Successor::NONE)
	    this->successorValue = Successor::REPLAN;
    }

    void requestReallocation(void)
    {
	this->successorValue = Successor::REALLOCATE;
    }

    Successor consume(void)
    {
	const Successor successor = this->successorValue;
	this->successorValue = Successor::NONE;
	return successor;
    }

    void reset(void)
    {
	this->successorValue = Successor::NONE;
    }

    bool pending(void) const
    {
	return this->successorValue != Successor::NONE;
    }

    Successor successor(void) const
    {
	return this->successorValue;
    }

private:
    Successor successorValue = Successor::NONE;
};

static_assert(std::is_trivially_copyable<
	BObolLodCapacityFrameBarrier>::value,
    "capacity frame barriers must remain allocation-free values");

/*
 * One retained occurrence-allocation request may be pending at a time.  The
 * former pending/preserve/reconciliation-budget field trio admitted invalid
 * combinations (for example, a reconciliation budget with no request) and
 * required callers to reset the companions together.  Keep that control edge
 * as a finite value; the numeric budget remains meaningful only for the
 * presentation-reconciliation alternative.
 */
class BObolLodRetainedAllocationRequest {
public:
    enum class Kind : uint8_t {
	NONE = 0,
	PRESERVE_BUDGET,
	RECOMPUTE_BUDGET,
	PRESENTATION_RECONCILIATION
    };

    void reset(void)
    {
	this->kindValue = Kind::NONE;
	this->reconciliationBudgetValue = SIZE_MAX;
    }

    void requestReallocation(bool preserveBudget)
    {
	/* A reconciliation request carries a completed-frame capacity proof.
	 * A weaker reallocation request cannot discard or reinterpret it. */
	if (this->kindValue == Kind::PRESENTATION_RECONCILIATION)
	    return;
	this->kindValue = preserveBudget ? Kind::PRESERVE_BUDGET :
	    Kind::RECOMPUTE_BUDGET;
	this->reconciliationBudgetValue = SIZE_MAX;
    }

    void requestPresentationReconciliation(size_t budget)
    {
	if (!budget || budget == SIZE_MAX)
	    return;
	if (this->kindValue != Kind::PRESENTATION_RECONCILIATION) {
	    this->kindValue = Kind::PRESENTATION_RECONCILIATION;
	    this->reconciliationBudgetValue = budget;
	    return;
	}
	this->reconciliationBudgetValue = std::min(
	    this->reconciliationBudgetValue, budget);
    }

    bool pending(void) const { return this->kindValue != Kind::NONE; }
    bool preservesBudget(void) const
    {
	return this->kindValue == Kind::PRESERVE_BUDGET;
    }
    bool reconcilesPresentation(void) const
    {
	return this->kindValue == Kind::PRESENTATION_RECONCILIATION;
    }
    Kind kind(void) const { return this->kindValue; }
    size_t reconciliationBudget(void) const
    {
	return this->reconcilesPresentation() ?
	    this->reconciliationBudgetValue : SIZE_MAX;
    }

private:
    Kind kindValue = Kind::NONE;
    size_t reconciliationBudgetValue = SIZE_MAX;
};

static_assert(std::is_trivially_copyable<
	BObolLodRetainedAllocationRequest>::value,
    "retained allocation requests must remain allocation-free values");

/*
 * Mechanical progress through one bounded admission plan.  This value owns no
 * renderer-capacity conclusion and survives no explicit reset.  Separating it
 * from completed-frame evidence lets a large scene resume across owner-thread
 * windows without making cursor movement look like new planning evidence.
 */
class BObolLodAdmissionCursor {
public:
    void reset(void)
    {
	this->initializedValue = false;
	this->singleOccurrenceBootstrapValue = false;
	this->activeCostValue = 0;
	this->minimumActiveCostValue = 0;
	this->refinementRemainingValue = SIZE_MAX;
	this->retainedAdmissionValue = false;
	this->retainedAdmissionRemainingValue = SIZE_MAX;
	this->revisionStampValue = BObolLodAdmissionRevisionStamp();
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

    bool initialized(void) const { return this->initializedValue; }
    bool singleOccurrenceBootstrap(void) const
    {
	return this->singleOccurrenceBootstrapValue;
    }
    size_t activeCost(void) const { return this->activeCostValue; }
    size_t minimumActiveCost(void) const
    {
	return this->minimumActiveCostValue;
    }
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
    bool matches(const BObolLodAdmissionRevisionStamp &stamp) const
    {
	return this->revisionStampValue.same(stamp);
    }
    const BObolLodAdmissionRevisionStamp &revisionStamp(void) const
    {
	return this->revisionStampValue;
    }

private:
    friend class BObolLodCapacityEvidence;
    friend class BObolLodAdmissionPlanner;

    void resetForRevision(const BObolLodAdmissionRevisionStamp &stamp)
    {
	this->reset();
	this->revisionStampValue = stamp;
    }

    bool initializedValue = false;
    bool singleOccurrenceBootstrapValue = false;
    size_t activeCostValue = 0;
    size_t minimumActiveCostValue = 0;
    size_t refinementRemainingValue = SIZE_MAX;
    bool retainedAdmissionValue = false;
    size_t retainedAdmissionRemainingValue = SIZE_MAX;
    BObolLodAdmissionRevisionStamp revisionStampValue;
};

static_assert(std::is_trivially_copyable<BObolLodAdmissionCursor>::value,
    "admission cursors must remain allocation-free values");

class BObolLodCapacityEvidence {
public:
    struct Inputs {
	size_t activeCost = 0;
	size_t minimumActiveCost = 0;
	float targetFps = 0.0f;
	long double calibratedCostPerSecond = 0.0L;
	uint64_t observedStableNanoseconds = 0;
	/* The preferred quiet cadence controls refinement pacing.  A completed
	 * frame below this separate hard limit is still a valid static
	 * presentation and must not be destructively coarsened merely for missing
	 * the preferred cadence. */
	uint64_t hardPresentationDeadlineNanoseconds = 0;
	uint64_t lastRenderNanoseconds = 0;
	uint64_t smoothedRenderNanoseconds = 0;
	bool interactive = false;
	bool scaleQualityProbe = false;
	/* A bounded static presentation may use the complete hard-frame allowance
	 * immediately.  This covers both a single-occurrence quality step and the
	 * minimum-mesh repair which replaces the last structural boxes in a
	 * many-part view.  The endpoint deadline can abort an optimistic trial, so
	 * imposing the ordinary 4x/1.25x refinement ladder here only creates
	 * visible intermediate stages or a no-progress repair loop. */
	bool hardDeadlinePresentation = false;
	/* Replacing a visible structural box with the minimum drawable mesh is a
	 * coverage obligation, not ordinary triangle-detail recovery.  A retained
	 * recovery ceiling may continue to constrain richer prefixes after this
	 * pass, but it must not make that minimum replacement impossible when the
	 * independently interruptible hard-frame allowance proves room for it. */
	bool structuralCoverageRepair = false;
	/* With exactly one visible occurrence and no retained presentation, a
	 * conservative many-object seed manufactures several blocky calibration
	 * frames.  Permit one larger, still deadline-protected bootstrap. */
	bool coldSingleOccurrence = false;
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
	uint64_t observedNanoseconds = 0;
	bool passAdmittedWork = false;
	/* The blocked cut came from a complete scene-wide minimax allocation.
	 * If calibration demonstrates a different allowance, its successor must
	 * recompute that allocation rather than letting ordinary per-occurrence
	 * refinement consume stale allocation stamps. */
	bool allocationCutsApplied = false;
	/* A complete retained allocation supplies the ordered candidate set and
	 * its pixel-demand endpoint.  Only this path may use the bounded search;
	 * structural/provider barriers retain their separate finite contracts. */
	bool boundedSearch = false;
	/* The occurrence-local cuts named by the allocation may be active while a
	 * temporary renderer-wide ceiling still hides their exact presentation.
	 * Such a population may advance ALLOCATING to MEASURING, but it may not
	 * consume a timing sample until the ceiling-free presentation is exact. */
	bool allocationPresentationRealized = false;
	BObolLodCapacitySearchKey searchKey;
	/* The allocator may raise a requested steady allowance to an atomic
	 * prominent-feature floor.  Name the actual certified candidate budget;
	 * it is not necessarily the scalar allowance which initiated the pass. */
	size_t candidateBudget = 0;
	size_t demandedBudget = 0;
	size_t knownSafeBudget = 0;
    };

    struct CompletedFrameInputs {
	BObolLodCapacitySearchKey searchKey;
	size_t candidateBudget = 0;
	size_t presentedCost = 0;
	size_t knownSafeBudget = 0;
	uint64_t observedNanoseconds = 0;
	bool validSample = false;
    };

    struct CalibrationDecision {
	bool candidateReallocation = false;
	bool searchActive = false;
	bool requestFrame = false;
	bool sampleFrame = false;
	bool restartSubmission = false;
	bool searchTerminal = false;
	BObolLodCapacitySearchCertificate::Result searchResult =
	    BObolLodCapacitySearchCertificate::Result::NONE;
    };

    struct CompletedFrameDecision {
	bool requestSampleFrame = false;
	bool restartSubmission = false;
    };

private:
    friend class BObolLodAdmissionPlanner;

    Decision planPass(const Inputs &inputs, BObolLodAdmissionCursor &cursor)
    {
	Decision decision;
	decision.totalBudget = this->currentBudgetValue;
	decision.refinementBudget = cursor.refinementRemainingValue;
	decision.retainedAdmission = cursor.retainedAdmissionValue;
	decision.retainedAdmissionBudget =
	    cursor.retainedAdmissionRemainingValue;
	if (cursor.initializedValue)
	    return decision;
	cursor.singleOccurrenceBootstrapValue =
	    inputs.coldSingleOccurrence && inputs.activeCost == 0 &&
	    inputs.calibratedCostPerSecond <= 0.0L;
	const bool requestedRetainedRecovery =
	    this->requestedRetainedRecoveryBudgetValue != SIZE_MAX &&
	    inputs.activeCost > this->requestedRetainedRecoveryBudgetValue;
	const bool requestedRetainedReallocation =
	    this->retainedAllocationRequestValue.pending() &&
	    !inputs.interactive && !inputs.forceTerminal &&
	    inputs.activeCost > 0;
	const bool requestedPresentationReconciliation =
	    requestedRetainedReallocation &&
	    this->retainedAllocationRequestValue.reconcilesPresentation();
	const size_t presentationReconciliationBudget =
	    this->retainedAllocationRequestValue.reconciliationBudget();
	const bool requestedCoverageCompletion =
	    inputs.structuralCoverageRepair &&
	    this->requestedCoverageCompletionAdditionalCostValue != SIZE_MAX;
	const size_t coverageCompletionBudget = requestedCoverageCompletion &&
	    this->requestedCoverageCompletionAdditionalCostValue >
		SIZE_MAX - inputs.activeCost ? SIZE_MAX :
	    inputs.activeCost +
		this->requestedCoverageCompletionAdditionalCostValue;
	const size_t retainedReallocationBudget = this->currentBudgetValue;

	const long double targetNanoseconds = inputs.targetFps > 0.0f ?
	    1000000000.0L / static_cast<long double>(inputs.targetFps) : 0.0L;
	const long double observedStableNanoseconds =
	    static_cast<long double>(inputs.observedStableNanoseconds);
	const bool severeStableOverload =
	    !inputs.interactive && targetNanoseconds > 0.0L &&
	    observedStableNanoseconds > targetNanoseconds * 2.0L;
	const bool observedWithinHardPresentationDeadline =
	    !inputs.interactive && inputs.activeCost > 0 &&
	    inputs.observedStableNanoseconds > 0 &&
	    inputs.hardPresentationDeadlineNanoseconds > 0 &&
	    inputs.observedStableNanoseconds <=
		inputs.hardPresentationDeadlineNanoseconds;
	const bool deadlineRecoveryRequested =
	    this->steadyDeadlineCapacityCeilingValue != SIZE_MAX &&
	    inputs.activeCost > this->steadyDeadlineCapacityCeilingValue;
	const bool overloadRecovery =
	    !inputs.interactive && !inputs.forceTerminal &&
	    !observedWithinHardPresentationDeadline &&
	    (!this->overloadRecoveryPerformedValue ||
	     this->overloadRecoveryActiveCostValue != inputs.activeCost) &&
	    inputs.activeCost > 0 &&
	    (deadlineRecoveryRequested || severeStableOverload) &&
	    targetNanoseconds > 0.0L &&
	    observedStableNanoseconds > targetNanoseconds * 1.20L;
	/* Three unchanged, exact presentation misses are stronger evidence than
	 * the triangle-throughput estimator.  The latter cannot represent the
	 * fixed per-occurrence/classification cost of a many-part scene: Hubble's
	 * measured triangle rate predicted headroom while three successive 75-85
	 * ms frames missed the 50 ms stable target.  Requiring both witnesses made
	 * the controller shave a few cost units, repeat three probes, and take
	 * tens of seconds to settle.  The probe count already rejects transient
	 * setup noise; use its direct duration with the safety-scaled recovery
	 * calculation below. */

	size_t costBudget = inputs.coldSingleOccurrence &&
	    inputs.activeCost == 0 &&
	    inputs.calibratedCostPerSecond <= 0.0L ?
		std::max(this->seedBudgetValue,
		    this->singleOccurrenceBootstrapBudget()) :
		this->seedBudgetValue;
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
	    if (inputs.activeCost > 0 && costBudget > inputs.activeCost &&
		!inputs.hardDeadlinePresentation) {
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
		!inputs.stablePresentationHandoff &&
		this->currentBudgetValue != SIZE_MAX)
		costBudget = std::max(costBudget, this->currentBudgetValue);
	    if (!inputs.interactive && !inputs.stablePresentationHandoff &&
		this->currentBudgetValue != SIZE_MAX &&
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
	/* observedStableNanoseconds is supplied only for an exact, reusable CAD
	 * presentation.  If that completed population met this pass's deadline,
	 * it is a direct proof that the population is affordable.  A triangle-rate
	 * estimate cannot contradict that proof: its linear currency omits fixed
	 * per-occurrence and command-dispatch costs and can otherwise place a
	 * many-part scene perpetually below its already displayed population. */
	const bool measuredActivePopulationDeadlineSafe =
	    !inputs.interactive && !inputs.stablePresentationHandoff &&
	    inputs.activeCost > 0 &&
	    targetNanoseconds > 0.0L && observedStableNanoseconds > 0.0L &&
	    observedStableNanoseconds <= targetNanoseconds;
	if (measuredActivePopulationDeadlineSafe && costBudget != SIZE_MAX)
	    costBudget = std::max(costBudget, inputs.activeCost);
	/* Static presentation is allowed to use the hard deadline.  Preserve an
	 * exact population which meets it even when the preferred quiet target is
	 * stricter; otherwise a pose-only rotation can replace a proven mesh view
	 * with boxes before structural repair has a chance to run. */
	if (observedWithinHardPresentationDeadline && costBudget != SIZE_MAX)
	    costBudget = std::max(costBudget, inputs.activeCost);
	/* A reusable presentation of this exact population is also the strongest
	 * available headroom estimate for the interruptible static-quality pass.
	 * The long-horizon throughput EMA can lag far behind after point
	 * aggregation or a large structural publication wave; using only that EMA
	 * left a 27 ms OSMesa frame with no allowance to replace its last few
	 * hundred boxes.  Extrapolate conservatively to eighty percent of the hard
	 * deadline and cap one step at four times the measured population.  The
	 * endpoint deadline remains the independent correctness guard if the cost
	 * model is non-linear. */
	if (inputs.hardDeadlinePresentation &&
	    measuredActivePopulationDeadlineSafe && costBudget != SIZE_MAX) {
	    const long double directAffordable =
		static_cast<long double>(inputs.activeCost) *
		targetNanoseconds * 0.80L / observedStableNanoseconds;
	    const size_t directLimit = inputs.activeCost > SIZE_MAX / 4 ?
		SIZE_MAX : inputs.activeCost * 4;
	    const size_t directBudget =
		!std::isfinite(directAffordable) || directAffordable <= 0.0L ?
		    inputs.activeCost :
		    (directAffordable >= static_cast<long double>(directLimit) ?
			directLimit : static_cast<size_t>(directAffordable));
	    costBudget = std::max(costBudget, directBudget);
	}
	/* Preserve only work which an exact completed presentation proved could
	 * meet the stable deadline.  This is deliberately a measured render-cost
	 * floor, not the previous pass allowance: the latter may include admission
	 * which was never submitted or presented. */
	if (inputs.stablePresentationCostFloor > 0 && costBudget != SIZE_MAX)
	    costBudget = std::max(
		costBudget, inputs.stablePresentationCostFloor);
	/* A view-importance reallocation deliberately changes only which
	 * occurrences consume the established allowance.  Re-running the
	 * throughput bootstrap here would replace that allowance with the cold
	 * seed before the first candidate is considered. */
	if (requestedRetainedReallocation &&
	    this->retainedAllocationRequestValue.preservesBudget())
	    costBudget = retainedReallocationBudget;
	/* Point aggregation is a last-resort population bound, not a substitute
	 * for PoP triangle allocation.  When the presentation controller is
	 * returning a temporary multi-pixel point cut to the one-pixel contract,
	 * it may explicitly ask this pass to compact reducible retained prefixes
	 * first.  The first request forces immediate retained admission; its
	 * measured ceiling remains in force for this view/policy epoch.  Without
	 * that persistent ceiling the cheaper recovery frame recalibrates a larger
	 * budget and immediately re-admits the exact discrete cut which just missed
	 * the frame deadline. */
	if (requestedRetainedRecovery)
	    costBudget = std::min(
		costBudget, this->requestedRetainedRecoveryBudgetValue);
	if (!inputs.forceTerminal && !inputs.structuralCoverageRepair &&
	    this->retainedRecoveryCeilingValue != SIZE_MAX)
	    costBudget = std::min(
		costBudget, this->retainedRecoveryCeilingValue);
	/* Recognizable prominent geometry is a harder quiet-view contract than
	 * the preferred 20 Hz cadence, provided its measured population remains
	 * inside the separate hard quality-frame deadline.  Persist that floor
	 * across overload passes; otherwise every 60-90 ms OSMesa wire frame
	 * alternates between the 20 Hz recovery budget and the same protected
	 * allocation forever. */
	if (!inputs.interactive && !inputs.forceTerminal &&
	    this->retainedQualityFloorBudgetValue > 0)
	    costBudget = std::max(
		costBudget, this->retainedQualityFloorBudgetValue);
	/* A completed hard-deadline abort is an unsafe upper-bound witness for
	 * this view/policy capacity epoch.  The constrained handoff below owns the
	 * immediate safe allocation, but its budget is intentionally one-shot.
	 * Keep a slightly lower strict ceiling after that handoff so ordinary
	 * throughput calibration cannot propose the identical failed population
	 * again.  A later static-quality pass may still try an intermediate budget
	 * below this bound, preserving useful fidelity without reopening a cycle. */
	const bool staticCapacitySearch =
	    this->capacitySearchValue.phase() !=
		BObolLodCapacitySearchCertificate::Phase::INACTIVE &&
	    this->capacitySearchValue.goal() ==
		BObolLodCapacitySearchCertificate::Goal::STATIC;
	const bool capacitySearchInactive =
	    this->capacitySearchValue.phase() ==
		BObolLodCapacitySearchCertificate::Phase::INACTIVE;
	const bool staticQualityBudget = staticCapacitySearch ||
	    (capacitySearchInactive &&
	     this->retainedQualityFloorBudgetValue > 0);
	const size_t deadlineCapacityCeiling = staticQualityBudget ?
	    this->staticDeadlineCapacityCeilingValue :
	    this->steadyDeadlineCapacityCeilingValue;
	if (!inputs.interactive && !inputs.forceTerminal &&
	    deadlineCapacityCeiling != SIZE_MAX)
	    costBudget = std::min(
		costBudget, deadlineCapacityCeiling);
	/* A completed constrained framebuffer is a hard presentation-capacity
	 * witness, not another soft quality preference.  Reconcile the hidden
	 * retained prefixes at exactly that scene-wide allowance: a stale
	 * throughput seed may not undercut proven work, and a protected-quality
	 * floor may not exceed the frame deadline which created the handoff. */
	if (requestedPresentationReconciliation)
	    costBudget = presentationReconciliationBudget;
	/* Structural fallback above the maximum point-proxy threshold is an
	 * irreducible coverage population, not ordinary quality refinement.  An
	 * exact completed frame may transfer its separately bounded static-frame
	 * capacity into one coverage transaction.  Consume that certified budget
	 * verbatim after the heuristic/EMA clamps: allowing the preferred-cadence
	 * estimator to undercut it re-arms the same box frontier forever. */
	if (requestedCoverageCompletion)
	    costBudget = coverageCompletionBudget;

	this->currentBudgetValue = costBudget;
	/* Freeze the exact cost currencies which initialized this bounded pass.
	 * A 150k-occurrence pass may require hundreds of owner-thread windows;
	 * rediscovering the complete population before every window turns an O(N)
	 * allocation into O(N * windows) work and monopolizes the GUI thread.
	 * Occurrence changes are charged by the carried refinement remainders; a
	 * completed pass resets this snapshot before another census. */
	cursor.activeCostValue = inputs.activeCost;
	cursor.minimumActiveCostValue = inputs.minimumActiveCost;
	cursor.refinementRemainingValue = costBudget == SIZE_MAX ? SIZE_MAX :
	    (costBudget > inputs.activeCost ?
		costBudget - inputs.activeCost : 0);
	/* A presentation handoff removes a reversible renderer ceiling; it does
	 * not authorize rewriting retained occurrence cuts.  In particular, a
	 * zoom must begin at the last coherent cut and adjust from there instead
	 * of normalizing every occurrence to a cheap baseline.  Cold coverage and
	 * point-proxy recovery request retained admission explicitly, while a
	 * completed over-budget quiet frame uses overloadRecovery. */
	cursor.retainedAdmissionValue =
	    (requestedRetainedReallocation && costBudget != SIZE_MAX) ||
	    ((inputs.interactive || overloadRecovery ||
	      requestedRetainedRecovery) &&
	     costBudget != SIZE_MAX && inputs.activeCost > costBudget);
	if (overloadRecovery) {
	    this->overloadRecoveryPerformedValue = true;
	    this->overloadRecoveryActiveCostValue = inputs.activeCost;
	}
	/* The retained scene already owns one minimum drawable prefix per active
	 * occurrence.  Bounded recovery windows therefore consume only detail
	 * above that irreducible coverage floor.  Charging the complete budget in
	 * every window drove the whole scene to minimum before a second pass
	 * rebuilt it, producing visible cycling and destroying importance order. */
	cursor.retainedAdmissionRemainingValue =
	    cursor.retainedAdmissionValue ?
		(costBudget > inputs.minimumActiveCost ?
		    costBudget - inputs.minimumActiveCost : 0) : SIZE_MAX;
	cursor.initializedValue = true;
	this->requestedRetainedRecoveryBudgetValue = SIZE_MAX;
	this->retainedAllocationRequestValue.reset();
	this->requestedCoverageCompletionAdditionalCostValue = SIZE_MAX;

	decision.initialized = true;
	decision.overloadRecovery = overloadRecovery;
	decision.totalBudget = this->currentBudgetValue;
	decision.refinementBudget = cursor.refinementRemainingValue;
	decision.retainedAdmission = cursor.retainedAdmissionValue;
	decision.retainedAdmissionBudget =
	    cursor.retainedAdmissionRemainingValue;
	return decision;
    }

    CompletedFrameDecision applyCapacitySearchDecision(
	const BObolLodCapacitySearchCertificate::Decision &search,
	BObolLodAdmissionCursor *cursor)
    {
	CompletedFrameDecision decision;
	if (search.requestsFrame()) {
	    this->stableBudgetLimitedValue = false;
	    decision.requestSampleFrame = true;
	    return decision;
	}

	if (search.requestsReallocation()) {
	    this->currentBudgetValue = search.budget;
	    this->retainedAllocationRequestValue.requestReallocation(true);
	    /* The search result names the next candidate population exactly.  The
	     * retained allocator must apply that budget before any throughput EMA,
	     * deadline floor, or growth heuristic may propose another one.  Letting
	     * planPass() recompute it makes the following certificate observe a
	     * different candidate and correctly reject the transition as stale. */
	    this->stableBudgetLimitedValue = false;
	    if (cursor)
		cursor->reset();
	    decision.restartSubmission = true;
	    return decision;
	}

	if (!search.terminal())
	    return decision;
	if (search.result ==
		BObolLodCapacitySearchCertificate::Result::STALE_POPULATION) {
	    this->capacitySearchValue.reset();
	    this->stableBudgetLimitedValue = false;
	    if (cursor)
		cursor->reset();
	    decision.restartSubmission = true;
	    return decision;
	}

	this->stableBudgetLimitedValue = true;
	if (search.result ==
		BObolLodCapacitySearchCertificate::Result::CERTIFIED) {
	    const size_t certifiedBudget = std::max<size_t>(1, search.budget);
	    if (certifiedBudget != this->currentBudgetValue) {
		this->currentBudgetValue = certifiedBudget;
		this->retainedAllocationRequestValue.requestReallocation(true);
		if (cursor)
		    cursor->reset();
		decision.restartSubmission = true;
	    }
	}
	return decision;
    }

public:

    CalibrationDecision finishBlockedPass(
	const CalibrationInputs &inputs)
    {
	CalibrationDecision decision;
	if (inputs.boundedSearch && inputs.searchKey.valid() &&
	    inputs.candidateBudget > 0 &&
	    inputs.allocationCutsApplied && inputs.demandedBudget > 0) {
	    BObolLodCapacitySearchCertificate::Observation observation;
	    observation.key = inputs.searchKey;
	    observation.candidateBudget = inputs.candidateBudget;
	    /* A changed allocation pass freezes its pre-pass active cost in the
	     * mechanical cursor.  Bind the new candidate only when its completed
	     * frame supplies the post-commit population. */
	    const bool presentationPending = inputs.passAdmittedWork ||
		!inputs.allocationPresentationRealized;
	    observation.presentedCost = presentationPending ? 0 :
		inputs.activeCost;
	    observation.knownSafeBudget = inputs.knownSafeBudget;
	    observation.observedNanoseconds = inputs.observedNanoseconds;
	    observation.validSample = !presentationPending &&
		inputs.observedNanoseconds > 0;
	    const BObolLodCapacitySearchCertificate::Decision search =
		presentationPending ?
		    this->capacitySearchValue.prepare(observation) :
		    this->capacitySearchValue.observe(observation);
	    const CompletedFrameDecision transition =
		this->applyCapacitySearchDecision(search, NULL);
	    decision.searchActive = true;
	    decision.candidateReallocation = search.requestsReallocation();
	    decision.requestFrame = transition.requestSampleFrame;
	    decision.sampleFrame = transition.requestSampleFrame;
	    decision.restartSubmission = transition.restartSubmission;
	    decision.searchTerminal = search.terminal();
	    decision.searchResult = search.result;
	    return decision;
	}
	/* An ordinary per-occurrence pass does not own the complete ordered demand
	 * endpoint required by the capacity-search certificate.  If it changed the
	 * framebuffer, present that population exactly once.  Its successor is a
	 * complete retained importance allocation; if it changed nothing, begin
	 * that allocation immediately.  Replaying three untyped frames here was a
	 * second capacity controller and could alternate with the certified search. */
	decision.searchActive = false;
	decision.candidateReallocation = false;
	decision.sampleFrame = false;
	this->stableBudgetLimitedValue = false;
	if (inputs.passAdmittedWork) {
	    this->frameBarrierValue.requestReallocation();
	    decision.requestFrame = true;
	} else {
	    this->retainedAllocationRequestValue.requestReallocation(false);
	    decision.restartSubmission = true;
	}
	return decision;
    }

    CompletedFrameDecision completeCalibrationFrame(
	BObolLodAdmissionCursor &cursor)
    {
	CompletedFrameDecision decision;
	if (!this->frameBarrierValue.pending())
	    return decision;
	const BObolLodCapacityFrameBarrier::Successor successor =
	    this->frameBarrierValue.consume();
	if (successor ==
		BObolLodCapacityFrameBarrier::Successor::REALLOCATE) {
	    this->retainedAllocationRequestValue.requestReallocation(false);
	}
	cursor.reset();
	decision.restartSubmission = true;
	return decision;
    }

    CompletedFrameDecision completeCapacitySearchFrame(
	BObolLodAdmissionCursor &cursor,
	const CompletedFrameInputs &inputs)
    {
	if (this->capacitySearchValue.phase() !=
		BObolLodCapacitySearchCertificate::Phase::MEASURING)
	    return this->completeCalibrationFrame(cursor);
	BObolLodCapacitySearchCertificate::Observation observation;
	observation.key = inputs.searchKey;
	observation.candidateBudget = inputs.candidateBudget;
	observation.presentedCost = inputs.presentedCost;
	observation.knownSafeBudget = inputs.knownSafeBudget;
	observation.observedNanoseconds = inputs.observedNanoseconds;
	observation.validSample = inputs.validSample;
	return this->applyCapacitySearchDecision(
	    this->capacitySearchValue.observe(observation), &cursor);
    }

    /* A calibration barrier can outlive the population it was meant to
     * measure (most notably when a structural-proxy coverage wave temporarily
     * has no active managed mesh cost).  Such a frame is presentation-only:
     * replaying it cannot produce capacity evidence.  Retire the obsolete
     * probe series and return directly to admission instead of waiting on an
     * impossible sample forever. */
    CompletedFrameDecision retireUnmeasurableCalibrationFrame(
	BObolLodAdmissionCursor &cursor)
    {
	CompletedFrameDecision decision;
	if (!this->frameBarrierValue.pending())
	    return decision;
	decision = this->completeCalibrationFrame(cursor);
	return decision;
    }

    void requestRescanAfterFrame(void)
    {
	/* A coverage/publication edge changes the population being presented.
	 * Any numeric candidate awaiting a sample is therefore stale; the new
	 * population must obtain a fresh complete allocation after this barrier. */
	this->capacitySearchValue.reset();
	this->frameBarrierValue.requestReplan();
	this->stableBudgetLimitedValue = false;
    }

    void clearBudgetLimit(void)
    {
	this->stableBudgetLimitedValue = false;
    }

    void resetCalibration(void)
    {
	this->frameBarrierValue.reset();
	this->stableBudgetLimitedValue = false;
	this->capacitySearchValue.reset();
    }

    void resetOverloadRecovery(void)
    {
	this->overloadRecoveryPerformedValue = false;
	this->overloadRecoveryActiveCostValue = 0;
    }

    void reset(void)
    {
	this->currentBudgetValue = this->seedBudgetValue;
	this->requestedRetainedRecoveryBudgetValue = SIZE_MAX;
	this->retainedAllocationRequestValue.reset();
	this->requestedCoverageCompletionAdditionalCostValue = SIZE_MAX;
	this->retainedRecoveryCeilingValue = SIZE_MAX;
	this->steadyDeadlineCapacityCeilingValue = SIZE_MAX;
	this->staticDeadlineCapacityCeilingValue = SIZE_MAX;
	this->retainedQualityFloorBudgetValue = 0;
	this->retainedQualityFloorSignatureValue = 0;
	this->retainedQualityFloorRejectedValue = false;
	this->retainedQualityFloorMissCountValue = 0;
	this->resetOverloadRecovery();
	this->resetCalibration();
    }

    void reduceCurrentBudget(size_t budget)
    {
	this->currentBudgetValue = std::min(
	    this->currentBudgetValue, budget);
    }

    void requestRetainedRecovery(size_t budget)
    {
	if (!budget || budget == SIZE_MAX)
	    return;
	this->requestedRetainedRecoveryBudgetValue = std::min(
	    this->requestedRetainedRecoveryBudgetValue, budget);
	this->retainedRecoveryCeilingValue = std::min(
	    this->retainedRecoveryCeilingValue, budget);
	this->currentBudgetValue = std::min(
	    this->currentBudgetValue, budget);
    }

    /* Normalize a retained population for one admission pass without
     * recording a capacity ceiling.  Motion/cold handoff uses this to obtain
     * a coherent measured baseline; subsequent passes must be free to spend
     * the headroom demonstrated by that baseline. */
    void requestRetainedNormalization(size_t budget)
    {
	if (!budget || budget == SIZE_MAX)
	    return;
	this->requestedRetainedRecoveryBudgetValue = std::min(
	    this->requestedRetainedRecoveryBudgetValue, budget);
	this->currentBudgetValue = std::min(
	    this->currentBudgetValue, budget);
    }

    /* Reallocate an already resident population within the current measured
     * scene allowance.  Unlike overload recovery this is not evidence of a
     * smaller renderer capacity and therefore installs no persistent budget
     * ceiling.  The request is consumed by exactly one complete admission
     * pass; a later view census must explicitly request another one. */
    void requestRetainedReallocation(bool preserveCurrentBudget = true)
    {
	this->retainedAllocationRequestValue.requestReallocation(
	    preserveCurrentBudget);
    }

    /* Convert a completed renderer-wide constrained frame into one atomic
     * occurrence-local importance allocation.  Unlike normalization, this
     * must run even when point aggregation makes the exact presented cost
     * cheaper than the hidden retained prefixes. */
    void requestPresentationReconciliation(size_t budget)
    {
	if (!budget || budget == SIZE_MAX)
	    return;
	this->retainedAllocationRequestValue.
	    requestPresentationReconciliation(budget);
	this->currentBudgetValue = budget;
    }

    /* Admit an irreducible structural frontier once using capacity proven by
     * the exact framebuffer immediately preceding it.  The request does not
     * rewrite retained occurrence cuts and installs no long-lived throughput
     * floor; it is consumed by exactly one structural repair pass.  An older
     * deadline ceiling describes a different occurrence allocation and may
     * not serialize this coverage-first candidate.  The candidate budget is
     * derived from the immediately preceding exact frame and the endpoint
     * hard deadline independently accepts or rejects the resulting batch. */
    size_t requestCoverageCompletion(size_t activeCost,
	size_t certifiedBudget)
    {
	if (certifiedBudget == SIZE_MAX || certifiedBudget <= activeCost)
	    return 0;
	/* The exact frame proves marginal capacity above the population it drew.
	 * Retained presentation metadata may change before the admission plan
	 * freezes its
	 * accounting baseline, so carrying the old absolute total would either
	 * over-admit or strand the complete repair frontier. */
	this->requestedCoverageCompletionAdditionalCostValue =
	    certifiedBudget - activeCost;
	this->currentBudgetValue = certifiedBudget;
	this->stableBudgetLimitedValue = false;
	return certifiedBudget;
    }

    /* Record a strict upper bound after a quiet, capacity-relevant render
     * abort.  Five percent separates the successor from harmless timer jitter
     * and, more importantly, proves monotonic progress if another
     * intermediate trial also misses.  The immediate recovery may choose a
     * more conservative budget; this ceiling governs later calibration once
     * that one-shot handoff is complete. */
    void noteDeadlineCapacityMiss(size_t attemptedBudget,
	bool staticDeadline = false)
    {
	if (attemptedBudget <= 1 || attemptedBudget == SIZE_MAX)
	    return;
	/* The interrupted frame is a typed unsafe witness owned by deadline
	 * recovery.  Retire any unchanged-frame candidate so its pending samples
	 * cannot later overwrite the stricter ceiling with stale safe evidence. */
	this->capacitySearchValue.reset();
	const size_t reduction = std::max<size_t>(1, attemptedBudget / 20);
	const size_t strictCeiling = attemptedBudget - reduction;
	if (staticDeadline) {
	    this->staticDeadlineCapacityCeilingValue = std::min(
		this->staticDeadlineCapacityCeilingValue, strictCeiling);
	    /* A population which misses the longer static deadline is also unsafe
	     * at the preferred cadence. */
	    this->steadyDeadlineCapacityCeilingValue = std::min(
		this->steadyDeadlineCapacityCeilingValue, strictCeiling);
	} else {
	    this->steadyDeadlineCapacityCeilingValue = std::min(
		this->steadyDeadlineCapacityCeilingValue, strictCeiling);
	}
	const size_t activeCeiling = staticDeadline ?
	    this->staticDeadlineCapacityCeilingValue :
	    this->steadyDeadlineCapacityCeilingValue;
	this->currentBudgetValue = std::min(
	    this->currentBudgetValue, activeCeiling);
	/* A protected floor above a hard failed-work bound cannot be defended by
	 * another soft allocation.  Retire it for this capacity epoch while
	 * keeping every immutable resident suffix available to a later view. */
	if (this->retainedQualityFloorBudgetValue >
	    this->staticDeadlineCapacityCeilingValue) {
	    this->retainedQualityFloorBudgetValue = 0;
	    this->retainedQualityFloorSignatureValue = 0;
	    this->retainedQualityFloorMissCountValue = 0;
	    this->retainedQualityFloorRejectedValue = true;
	}
    }

    /* A view-significance floor may deliberately trade some of the 20 Hz
     * stable target for recognizable prominent geometry, but only after the
     * controller proves that population fits its separately bounded hard
     * quality-frame allowance.  Update all three currencies of the already
     * initialized retained pass together; changing only currentBudgetValue
     * would leave its bounded actions consuming the former upgrade limit. */
    void setRetainedQualityFloorBudget(BObolLodAdmissionCursor &cursor,
	size_t budget,
	uint64_t populationSignature, size_t activeCost,
	size_t minimumActiveCost)
    {
	if (this->retainedQualityFloorRejectedValue)
	    return;
	const size_t nextBudget = budget == SIZE_MAX ? 0 : budget;
	if (nextBudget > 0 &&
	    this->staticDeadlineCapacityCeilingValue != SIZE_MAX &&
	    nextBudget > this->staticDeadlineCapacityCeilingValue) {
	    this->retainedQualityFloorBudgetValue = 0;
	    this->retainedQualityFloorSignatureValue = 0;
	    this->retainedQualityFloorMissCountValue = 0;
	    this->retainedQualityFloorRejectedValue = true;
	    return;
	}
	const uint64_t nextSignature = nextBudget ? populationSignature : 0;
	/* Deadline evidence is meaningful only for the exact protected
	 * occurrence/cut population which produced it.  During cold streaming or
	 * a retained reallocation, different floors can be attempted before an
	 * exact frame is presented.  Combining those misses can reject a later,
	 * affordable floor even though it has never missed once. */
	if (this->retainedQualityFloorSignatureValue != nextSignature)
	    this->retainedQualityFloorMissCountValue = 0;
	this->retainedQualityFloorBudgetValue = nextBudget;
	this->retainedQualityFloorSignatureValue = nextSignature;
	if (!cursor.initializedValue || !cursor.retainedAdmissionValue ||
	    !budget || budget == SIZE_MAX || budget <= this->currentBudgetValue)
	    return;
	this->currentBudgetValue = budget;
	cursor.refinementRemainingValue = budget > activeCost ?
	    budget - activeCost : 0;
	cursor.retainedAdmissionRemainingValue =
	    budget > minimumActiveCost ? budget - minimumActiveCost : 0;
	if (this->retainedRecoveryCeilingValue != SIZE_MAX &&
	    this->retainedRecoveryCeilingValue < budget)
	    this->retainedRecoveryCeilingValue = budget;
    }

    void clearRetainedQualityFloorBudget(void)
    {
	this->retainedQualityFloorBudgetValue = 0;
	this->retainedQualityFloorSignatureValue = 0;
	this->retainedQualityFloorRejectedValue = false;
	this->retainedQualityFloorMissCountValue = 0;
    }

    bool noteRetainedQualityFloorMiss(void)
    {
	if (!this->retainedQualityFloorBudgetValue ||
	    this->retainedQualityFloorRejectedValue)
	    return false;
	if (this->retainedQualityFloorMissCountValue < UINT_MAX)
	    this->retainedQualityFloorMissCountValue++;
	/* One cold upload, command-record transition, compositor stall, or
	 * software JIT event is not renderer-capacity evidence.  Require three
	 * independent ceiling-free attempts at the same view-local floor before
	 * allowing the soft allocator to sacrifice protected geometry. */
	if (this->retainedQualityFloorMissCountValue < 3)
	    return false;
	this->retainedQualityFloorBudgetValue = 0;
	this->retainedQualityFloorSignatureValue = 0;
	this->retainedQualityFloorRejectedValue = true;
	return true;
    }

    /* A protected floor is allowed to exceed the preferred quiet cadence
     * only while it remains inside the independently interruptible static
     * frame deadline.  A miss from the explicit static-quality trial is not
     * one of the transient setup samples counted above: that trial is the
     * endpoint contract itself.  Retire the floor for this policy epoch so a
     * recovery allocation cannot immediately restore the population which
     * just failed.  Resident PoP data and occurrence demand are untouched. */
    bool rejectRetainedQualityFloor(void)
    {
	const bool hadFloor = this->retainedQualityFloorBudgetValue > 0;
	this->retainedQualityFloorBudgetValue = 0;
	this->retainedQualityFloorSignatureValue = 0;
	this->retainedQualityFloorMissCountValue = 0;
	this->retainedQualityFloorRejectedValue = true;
	return hadFloor;
    }

    bool noteRetainedQualityFloorMet(bool exactProtectedPopulation,
	uint64_t populationSignature, size_t presentedCost)
    {
	if (!exactProtectedPopulation ||
	    !this->retainedQualityFloorBudgetValue ||
	    this->retainedQualityFloorRejectedValue ||
	    populationSignature != this->retainedQualityFloorSignatureValue ||
	    presentedCost < this->retainedQualityFloorBudgetValue)
	    return false;
	this->retainedQualityFloorMissCountValue = 0;
	return true;
    }

    bool retainedQualityFloorRejected(void) const
    {
	return this->retainedQualityFloorRejectedValue;
    }

    bool retainedQualityFloorActive(void) const
    {
	return this->retainedQualityFloorBudgetValue > 0;
    }

    size_t retainedQualityFloorBudget(void) const
    {
	return this->retainedQualityFloorBudgetValue;
    }

    uint64_t retainedQualityFloorSignature(void) const
    {
	return this->retainedQualityFloorSignatureValue;
    }

    unsigned int retainedQualityFloorMissCount(void) const
    {
	return this->retainedQualityFloorMissCountValue;
    }

    void clearRetainedRecoveryCeiling(void)
    {
	this->requestedRetainedRecoveryBudgetValue = SIZE_MAX;
	this->retainedRecoveryCeilingValue = SIZE_MAX;
    }

    void clearDeadlineCapacityCeiling(void)
    {
	this->steadyDeadlineCapacityCeilingValue = SIZE_MAX;
	this->staticDeadlineCapacityCeilingValue = SIZE_MAX;
    }

    /* A measured recovery ceiling protects the first coherent one-pixel
     * population from immediately re-admitting the cut which just missed its
     * deadline.  It must end once that population is actually ready for
     * presentation.  Keep this transition in the scalar policy so callers
     * cannot accidentally clear only the pass state, or retain the ceiling
     * forever when returning from a coarser point cut needs an extra frame. */
    bool confirmRetainedRecoveryPresentation(
	bool onePixelReady, BObolLodAdmissionCursor &cursor)
    {
	if (!onePixelReady ||
	    this->retainedRecoveryCeilingValue == SIZE_MAX)
	    return false;
	this->clearRetainedRecoveryCeiling();
	cursor.reset();
	return true;
    }

    void raiseCurrentBudget(size_t budget)
    {
	size_t raisedBudget = std::max(this->currentBudgetValue, budget);
	if (this->retainedRecoveryCeilingValue != SIZE_MAX)
	    raisedBudget = std::min(
		raisedBudget, this->retainedRecoveryCeilingValue);
	if (this->steadyDeadlineCapacityCeilingValue != SIZE_MAX)
	    raisedBudget = std::min(
		raisedBudget, this->steadyDeadlineCapacityCeilingValue);
	this->currentBudgetValue = raisedBudget;
    }

    size_t seedBudget(void) const { return this->seedBudgetValue; }
    static constexpr size_t singleOccurrenceBootstrapBudget(void)
    {
	/* A sole cold source has no competing occurrence whose first useful mesh
	 * could be starved.  Permit one globally classified, byte-capped PoP
	 * preview rather than rejecting it for a provisional 500k estimate and
	 * making the user wait for chunk persistence.  The endpoint's 100 ms hard
	 * deadline and completed-frame calibration remain the safety contract;
	 * this is only the bounded first-publication allowance. */
	return 1000000;
    }
    size_t currentBudget(void) const { return this->currentBudgetValue; }
    bool retainedRecoveryCeilingActive(void) const
    {
	return this->retainedRecoveryCeilingValue != SIZE_MAX;
    }
    size_t deadlineCapacityCeiling(void) const
    {
	return this->steadyDeadlineCapacityCeilingValue;
    }
    size_t staticDeadlineCapacityCeiling(void) const
    {
	return this->staticDeadlineCapacityCeilingValue;
    }
    bool overloadRecoveryPerformed(void) const
    {
	return this->overloadRecoveryPerformedValue;
    }
    size_t overloadRecoveryActiveCost(void) const
    {
	return this->overloadRecoveryActiveCostValue;
    }
    bool rescanAfterFrame(void) const
    {
	return this->frameBarrierValue.pending() ||
	    this->capacitySearchValue.awaitingSample();
    }
    bool stableBudgetLimited(void) const
    {
	return this->stableBudgetLimitedValue;
    }
    const BObolLodCapacitySearchCertificate &capacitySearch(void) const
    {
	return this->capacitySearchValue;
    }

private:
    size_t seedBudgetValue = 50000;
    size_t currentBudgetValue = 50000;
    size_t overloadRecoveryActiveCostValue = 0;
    size_t requestedRetainedRecoveryBudgetValue = SIZE_MAX;
    size_t retainedRecoveryCeilingValue = SIZE_MAX;
    size_t steadyDeadlineCapacityCeilingValue = SIZE_MAX;
    size_t staticDeadlineCapacityCeilingValue = SIZE_MAX;
    size_t retainedQualityFloorBudgetValue = 0;
    uint64_t retainedQualityFloorSignatureValue = 0;
    unsigned int retainedQualityFloorMissCountValue = 0;
    BObolLodRetainedAllocationRequest retainedAllocationRequestValue;
    size_t requestedCoverageCompletionAdditionalCostValue = SIZE_MAX;
    bool overloadRecoveryPerformedValue = false;
    bool stableBudgetLimitedValue = false;
    bool retainedQualityFloorRejectedValue = false;
    BObolLodCapacityFrameBarrier frameBarrierValue;
    BObolLodCapacitySearchCertificate capacitySearchValue;
};

static_assert(std::is_trivially_copyable<BObolLodCapacityEvidence>::value,
    "capacity evidence must remain an allocation-free value");

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
	 * 60 Hz motion cut.  Replaying it does not expose new immutable data and
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
	this->qualityFloorValue = -1;
	this->qualityCeilingLimitValue = -1;
	this->qualityCeilingFailedWorkValue = 0;
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
    bool scaleDemandRefreshActiveValue = false;
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
	PROTECTED_MINIMUM
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
	    size_t reconciliationBudget = 0)
    {
	if (!allocationCompleted || !allocationCertified ||
	    !certifiedPresentationBudget)
	    return false;
	const size_t budget = reconciliationBudget > 0 ?
	    std::min(certifiedPresentationBudget, reconciliationBudget) :
	    certifiedPresentationBudget;
	return selectedPresentationCost <= budget;
    }

    /* A complete occurrence allocation which exceeds its own certified
     * budget is not another capacity-search candidate.  It is the concrete
     * lower bound which the presentation-reconciliation policy must either
     * reduce (by aggregating independent small occurrences) or accept as a
     * bounded static constraint.  Claim that successor here so an inactive
     * handoff cannot leave the generic capacity fallback resubmitting the
     * identical allocation forever. */
    bool claimOverBudgetAllocation(bool allocationCurrent,
	    bool allocationCutsApplied, size_t selectedPresentationCost,
	    size_t certifiedPresentationBudget, size_t presentedRenderCost)
    {
	if (!allocationCurrent || !allocationCutsApplied ||
	    !selectedPresentationCost || !certifiedPresentationBudget ||
	    selectedPresentationCost <= certifiedPresentationBudget)
	    return false;
	if (!this->handoffActive())
	    this->armHandoff(false, presentedRenderCost,
		certifiedPresentationBudget);
	return true;
    }

    /* Capacity calibration may speculatively schedule another allocation
     * before the handoff reducer sees the pass which just completed.  Once
     * that pass supplies a current, quiescent occurrence-local certificate
     * for the constrained presentation, the speculative successor has no new
     * evidence to seek.  Let the stronger handoff proof consume it; concrete
     * cut, rescan, publication, or population work remains authoritative. */
    bool currentHandoffAllocationSupersedesCapacityRestart(
	    const CompletedPassInputs &inputs,
	    bool capacityRestartScheduled) const
    {
	/* Once a current, quiescent occurrence allocation exists, the handoff
	 * reducer owns the next decision.  A fitting allocation completes the
	 * handoff; an over-budget minimum requests local presentation reduction.
	 * Letting the generic capacity fallback restart first leaves
	 * submissionPending set, which prevents either handoff transition and
	 * creates an unbounded no-op allocation loop. */
	return capacityRestartScheduled && inputs.completed &&
	    inputs.submissionPending && !inputs.rescanAfterFrame &&
	    !inputs.changedCut && this->handoffActive() &&
	    !this->presentationHandoffPending() &&
	    inputs.retainedAllocationCompleted &&
	    inputs.retainedAllocationCertified &&
	    inputs.populationQuiescent;
    }

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
	    int progressiveCeiling) const
    {
	return capacitySamplePending && progressiveCeiling >= 0 &&
	    inputs.rescanAfterFrame && inputs.completed &&
	    !inputs.submissionPending &&
	    this->handoffActive() && !this->presentationHandoffPending() &&
	    inputs.retainedAllocationCompleted &&
	    inputs.retainedAllocationCertified &&
	    inputs.populationQuiescent &&
	    inputs.presentationLimitsReconciled;
    }

    void capturePrior(float targetPixelError, int progressiveCeiling,
	float progressiveNextFraction, float pointProxyPixelThreshold,
	const Population &population,
	BObolLodViewEpoch viewEpoch)
    {
	this->priorStableValue = makeSnapshot(
	    targetPixelError, progressiveCeiling, progressiveNextFraction,
	    pointProxyPixelThreshold, population,
	    viewEpoch, 0);
	this->provenQualityValue = Snapshot();
	this->clearHandoff();
    }

    void noteStableQuality(float targetPixelError, int progressiveCeiling,
	float progressiveNextFraction, float pointProxyPixelThreshold,
	const Population &population,
	BObolLodViewEpoch viewEpoch, size_t presentedRenderCost)
    {
	if (!presentedRenderCost)
	    return;
	const Snapshot candidate = makeSnapshot(
	    targetPixelError, progressiveCeiling, progressiveNextFraction,
	    pointProxyPixelThreshold, population,
	    viewEpoch, presentedRenderCost);
	if (this->provenQualityValue.valid &&
	    this->provenQualityValue.viewEpoch == viewEpoch &&
	    populationMatches(this->provenQualityValue.population, population) &&
	    this->provenQualityValue.presentedRenderCost > presentedRenderCost)
	    return;
	this->provenQualityValue = candidate;
    }

    RestoreDecision beginQuiet(const QuietInputs &quietInputs)
    {
	const RestoreInputs &inputs = quietInputs.presentation;
	using Reducer = BObolLodQuietSuccessorReducer;
	Reducer::Inputs successorInputs;
	successorInputs.retainedMeshPayloads = inputs.retainedMeshPayloads;
	successorInputs.current.valid = true;
	successorInputs.current.targetPixelError =
	    inputs.currentTargetPixelError;
	successorInputs.current.progressiveCeiling =
	    inputs.currentProgressiveCeiling;
	successorInputs.current.progressiveNextFraction =
	    inputs.currentProgressiveNextFraction;
	successorInputs.current.pointProxyPixelThreshold =
	    inputs.currentPointProxyPixelThreshold;

	const bool priorMatches = inputs.orthographic && !inputs.scaleChanged &&
	    this->priorStableValue.valid && populationMatches(
		this->priorStableValue.population, inputs.population);
	if (priorMatches)
	    successorInputs.priorStable = targetFromSnapshot(
		this->priorStableValue);
	const bool provenMatches = !priorMatches && inputs.scaleChanged &&
	    this->provenQualityValue.valid &&
	    this->provenQualityValue.viewEpoch == inputs.viewEpoch &&
	    populationMatches(
		this->provenQualityValue.population, inputs.population);
	if (provenMatches)
	    successorInputs.provenScale = targetFromSnapshot(
		this->provenQualityValue);
	if (quietInputs.exactView.valid) {
	    successorInputs.exactView.valid = true;
	    successorInputs.exactView.targetPixelError =
		quietInputs.exactView.targetPixelError;
	    successorInputs.exactView.progressiveCeiling =
		quietInputs.exactView.progressiveCeiling;
	    successorInputs.exactView.progressiveNextFraction =
		quietInputs.exactView.progressiveNextFraction;
	    successorInputs.exactView.pointProxyPixelThreshold =
		quietInputs.exactView.pointProxyPixelThreshold;
	    successorInputs.exactView.presentedRenderCost =
		quietInputs.exactView.presentedRenderCost;
	    successorInputs.exactView.provenRenderCostCapacity =
		quietInputs.exactView.provenRenderCostCapacity;
	}

	const Reducer::Decision successor = Reducer::reduce(successorInputs);
	RestoreDecision decision;
	decision.apply = successor.apply;
	decision.restoredPriorStable =
	    successor.source == Reducer::Source::PRIOR_STABLE;
	decision.restoredProvenQuality =
	    successor.source == Reducer::Source::PROVEN_SCALE;
	decision.restoredExactView =
	    successor.source == Reducer::Source::EXACT_VIEW;
	decision.handoff = successor.handoff;
	decision.clearPresentationLimits =
	    successor.clearPresentationLimits;
	decision.targetPixelError = successor.target.targetPixelError;
	decision.progressiveCeiling = successor.target.progressiveCeiling;
	decision.progressiveNextFraction =
	    successor.target.progressiveNextFraction;
	decision.pointProxyPixelThreshold =
	    successor.target.pointProxyPixelThreshold;
	decision.provenRenderCostFloor =
	    successor.target.presentedRenderCost;
	decision.provenRenderCostCapacity =
	    successor.target.provenRenderCostCapacity;

	this->handoffStateValue = decision.handoff == Reducer::Handoff::PRESENTATION ?
	    HandoffState::PRESENTATION_REQUIRED :
	    (decision.handoff == Reducer::Handoff::ALLOCATION ?
		HandoffState::ALLOCATION_REQUIRED : HandoffState::INACTIVE);
	this->handoffReconciliationBudgetValue = 0;
	this->handoffReconciliationBudgetLimitValue = 0;
	/* A proven scale target carries the exact cost of the completed frame it
	 * restores.  Other handoffs must establish their cost through normal
	 * completed-frame evidence. */
	this->handoffCostFloorValue =
	    decision.handoff == Reducer::Handoff::ALLOCATION &&
	    decision.restoredProvenQuality ?
		decision.provenRenderCostFloor : 0;

	this->priorStableValue = Snapshot();
	this->provenQualityValue = Snapshot();
	return decision;
    }

    RestoreDecision beginQuiet(const RestoreInputs &inputs)
    {
	QuietInputs quietInputs;
	quietInputs.presentation = inputs;
	return this->beginQuiet(quietInputs);
    }

    CompletedPassDecision completePass(const CompletedPassInputs &inputs)
    {
	CompletedPassDecision decision;
	if (!inputs.completed ||
	    inputs.submissionPending || inputs.rescanAfterFrame ||
	    inputs.changedCut ||
	    (inputs.retainedRefinementBudgetBlocked &&
	     !this->handoffActive()) ||
	    this->presentationHandoffPending())
	    return decision;
	if (!this->handoffActive()) {
	    /* Wanting a richer retained cut is an observation, not a liveness
	     * witness.  With no admitted cut, budget barrier, handoff, or rescan,
	     * the current calibrated cut is terminal until a new external edge
	     * starts another pass. */
	    decision.retireRetainedObservation =
		inputs.retainedRefinementPending;
	    return decision;
	}
	if (!inputs.retainedAllocationCompleted ||
	    !inputs.retainedAllocationCertified) {
	    decision.requestRetainedAllocation =
		inputs.populationQuiescent;
	    return decision;
	}
	if (!inputs.populationQuiescent)
	    return decision;
	if (!inputs.presentationLimitsReconciled) {
	    decision.requestLocalPresentationReduction = true;
	    return decision;
	}
	/* A hard-deadline ceiling may be removed only by an explicit complete-
	 * population reconciliation at the constrained frame's proven budget.
	 * Merely changing one or more occurrence cuts is not that proof: late
	 * resident results or a larger retained budget can leave the hidden
	 * ceiling-free population more expensive than the frame which completed.
	 * Clearing on "changed" created an exact cycle of constrained completion,
	 * ceiling release, deadline abort, and another full-scene allocation.
	 *
	 * An ordinary quiet handoff has no failed-frame safety witness to retain;
	 * there, a changed coherent allocation is sufficient and the historical
	 * budget-blocked/no-change rule remains appropriate. */
	this->clearHandoff();
	decision.finishHandoff = true;
	decision.requestRetainedRescan =
	    inputs.retainedRefinementPending &&
	    !inputs.retainedRefinementBudgetBlocked;
	return decision;
    }

    void armHandoff(bool presentationRequired, size_t provenRenderCost = 0,
	size_t reconciliationBudgetLimit = 0)
    {
	this->handoffStateValue = presentationRequired ?
	    HandoffState::PRESENTATION_REQUIRED :
	    HandoffState::ALLOCATION_REQUIRED;
	/* A presentation-required handoff follows a hard frame deadline.  Its
	 * first completed constrained frame supplies the cost proof below; do not
	 * immediately remove that cut and retry the frame which just exceeded the
	 * deadline.  A static predictor may instead arm an already-presented
	 * handoff with its exact cost floor and no additional presentation gate. */
	/* A deadline recovery must first present its corrected renderer cut, so
	 * noteFramePresented() derives the final reconciliation allowance while
	 * respecting this upper bound.  A static staircase has already presented
	 * the cut it is retaining; its completed-frame prediction is immediately
	 * usable by the occurrence allocator.  Keeping that value here prevents
	 * the allocator from falling back to a broader scene budget and selecting
	 * a hidden population which the renderer ceiling cannot safely release. */
	this->handoffReconciliationBudgetValue = presentationRequired ? 0 :
	    reconciliationBudgetLimit;
	this->handoffReconciliationBudgetLimitValue =
	    reconciliationBudgetLimit;
	this->handoffCostFloorValue = provenRenderCost;
    }
    bool noteFramePresented(size_t provenRenderCost = 0,
	size_t reconciliationBudget = 0)
    {
	const bool released = this->handoffStateValue ==
	    HandoffState::PRESENTATION_REQUIRED;
	if (this->handoffActive()) {
	    if (released)
		this->handoffStateValue = HandoffState::ALLOCATION_REQUIRED;
	    if (provenRenderCost)
		this->handoffCostFloorValue = std::max(
		    this->handoffCostFloorValue, provenRenderCost);
	    /* A hard-deadline handoff may retire its global renderer limit only
	     * after the retained allocator has encoded an equivalent bounded
	     * presentation in occurrence-local cuts.  Keep this exact-frame
	     * capacity witness distinct from the affordable-work floor above. */
	    if (released) {
		/* The constrained frame measures only the renderer-limited cut.  A
		 * preceding failed full frame may already have computed a stricter
		 * budget in the retained-population cost domain.  Never let fast
		 * timing of the coarse fallback erase that safety proof.
		 *
		 * A corrected frame can itself land just beyond the endpoint timer.
		 * In that case it supplies no new safe extrapolation, but completing
		 * the frame must not discard the recovery limit derived from the
		 * interrupted richer population.  Losing that limit lets the handoff
		 * accept the known-failing population and repeat the same deadline
		 * correction indefinitely. */
		size_t boundedBudget = reconciliationBudget;
		if (this->handoffReconciliationBudgetLimitValue > 0)
		    boundedBudget = boundedBudget > 0 ?
			std::min(boundedBudget,
			    this->handoffReconciliationBudgetLimitValue) :
			this->handoffReconciliationBudgetLimitValue;
		if (boundedBudget)
		    this->handoffReconciliationBudgetValue = boundedBudget;
	    }
	}
	return released;
    }
    void cancelHandoff(void)
    {
	this->clearHandoff();
    }

    void viewInvalidated(void)
    {
	this->clearHandoff();
	this->provenQualityValue = Snapshot();
    }

    void reset(void)
    {
	this->priorStableValue = Snapshot();
	this->provenQualityValue = Snapshot();
	this->clearHandoff();
    }

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

    void clearHandoff(void)
    {
	this->handoffStateValue = HandoffState::INACTIVE;
	this->handoffReconciliationBudgetValue = 0;
	this->handoffReconciliationBudgetLimitValue = 0;
	this->handoffCostFloorValue = 0;
    }

    static float sanitizePixelError(float pixelError)
    {
	return std::isfinite(pixelError) ?
	    std::max(0.25f, std::min(64.0f, pixelError)) : 1.0f;
    }

    static float sanitizeThreshold(float threshold)
    {
	return std::isfinite(threshold) ?
	    std::max(1.0f, std::min(64.0f, threshold)) : 1.0f;
    }

    static float sanitizeFraction(int progressiveCeiling, float fraction)
    {
	if (progressiveCeiling < 0 || !std::isfinite(fraction))
	    return 0.0f;
	return std::max(0.0f, std::min(1.0f, fraction));
    }

    static Snapshot makeSnapshot(float targetPixelError,
	int progressiveCeiling, float progressiveNextFraction,
	float pointProxyPixelThreshold, const Population &population,
	BObolLodViewEpoch viewEpoch, size_t presentedRenderCost)
    {
	Snapshot snapshot;
	snapshot.valid = true;
	snapshot.viewEpoch = viewEpoch;
	snapshot.targetPixelError = sanitizePixelError(targetPixelError);
	snapshot.progressiveCeiling =
	    progressiveCeiling < -1 ? -1 :
	    std::min<int>(BOBOL_MESH_LOD_CUT_COUNT_MAX - 1,
		progressiveCeiling);
	snapshot.progressiveNextFraction = sanitizeFraction(
	    snapshot.progressiveCeiling, progressiveNextFraction);
	snapshot.pointProxyPixelThreshold =
	    sanitizeThreshold(pointProxyPixelThreshold);
	snapshot.population = population;
	snapshot.presentedRenderCost = presentedRenderCost;
	return snapshot;
    }

    static BObolLodQuietSuccessorReducer::Target targetFromSnapshot(
	const Snapshot &snapshot)
    {
	BObolLodQuietSuccessorReducer::Target target;
	target.valid = snapshot.valid;
	target.targetPixelError = snapshot.targetPixelError;
	target.progressiveCeiling = snapshot.progressiveCeiling;
	target.progressiveNextFraction = snapshot.progressiveNextFraction;
	target.pointProxyPixelThreshold = snapshot.pointProxyPixelThreshold;
	target.presentedRenderCost = snapshot.presentedRenderCost;
	return target;
    }

    static bool populationMatches(const Population &snapshot,
	const Population &current)
    {
	return snapshot.available == current.available &&
	    snapshot.sceneDomainRevision == current.sceneDomainRevision;
    }

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

    bool remember(const RememberInputs &inputs)
    {
	if (!inputs.view.haveCamera || !inputs.domainRevision ||
	    !inputs.sceneAvailable || !inputs.exactCompletedFrame ||
	    !inputs.terminalPresentationComplete ||
	    !inputs.producersSettled ||
	    !inputs.presentationDeadlineMet || inputs.resourcePressure ||
	    !validQuality(inputs.quality))
	    return false;

	size_t index = this->entryCount;
	for (size_t i = 0; i < this->entryCount; ++i) {
	    if (this->entries[i].domainRevision == inputs.domainRevision &&
		this->entries[i].view.same(inputs.view)) {
		index = i;
		break;
	    }
	}

	/* Renderer throughput and visual fidelity are independent evidence.  A
	 * high-cost frame containing many minor triangles may still leave a wheel
	 * or blade visibly coarser than a lower-cost importance allocation.  Keep
	 * the largest completed workload as a capacity proof, but replace the
	 * complete fidelity snapshot only when screen-error evidence (or the
	 * conservative control-vector fallback) says it is better. */
	if (index < this->entryCount) {
	    Entry &current = this->entries[index];
	    current.provenRenderCostCapacity = std::max(
		current.provenRenderCostCapacity,
		inputs.quality.presentedRenderCost);
	    if (!preferFidelityCandidate(inputs.quality, current.quality)) {
		this->promote(index);
		return false;
	    }
	}

	Entry candidate;
	candidate.valid = true;
	candidate.view = inputs.view;
	candidate.domainRevision = inputs.domainRevision;
	candidate.quality = sanitizeQuality(inputs.quality);
	candidate.provenRenderCostCapacity =
	    index < this->entryCount ?
		this->entries[index].provenRenderCostCapacity :
		candidate.quality.presentedRenderCost;

	if (index == this->entryCount) {
	    if (this->entryCount < capacity())
		this->entryCount++;
	    index = this->entryCount - 1;
	}
	for (size_t i = index; i > 0; --i)
	    this->entries[i] = this->entries[i - 1];
	this->entries[0] = candidate;
	if (this->rememberCountValue != SIZE_MAX)
	    this->rememberCountValue++;
	return true;
    }

    bool recall(const RecallInputs &inputs, Quality &quality)
    {
	if (!inputs.view.haveCamera || !inputs.domainRevision ||
	    !inputs.sceneAvailable || inputs.resourcePressure)
	    return false;
	for (size_t i = 0; i < this->entryCount; ++i) {
	    if (!this->entries[i].valid ||
		this->entries[i].domainRevision != inputs.domainRevision ||
		!this->entries[i].view.same(inputs.view))
		continue;
	    quality = this->entries[i].quality;
	    quality.provenRenderCostCapacity =
		this->entries[i].provenRenderCostCapacity;
	    this->promote(i);
	    if (this->recallCountValue != SIZE_MAX)
		this->recallCountValue++;
	    return true;
	}
	return false;
    }

    void reset(void)
    {
	this->entries = {};
	this->entryCount = 0;
	this->rememberCountValue = 0;
	this->recallCountValue = 0;
    }

    size_t size(void) const { return this->entryCount; }
    size_t rememberCount(void) const { return this->rememberCountValue; }
    size_t recallCount(void) const { return this->recallCountValue; }

private:
    struct Entry {
	bool valid = false;
	BObolLodViewSnapshot view;
	uint64_t domainRevision = 0;
	Quality quality;
	size_t provenRenderCostCapacity = 0;
    };

    static bool validQuality(const Quality &quality)
    {
	return std::isfinite(quality.targetPixelError) &&
	    quality.targetPixelError >= 0.249f &&
	    quality.targetPixelError <= 1.01f &&
	    quality.progressiveCeiling >= -1 &&
	    quality.progressiveCeiling < BOBOL_MESH_LOD_CUT_COUNT_MAX &&
	    std::isfinite(quality.progressiveNextFraction) &&
	    quality.progressiveNextFraction >= 0.0f &&
	    quality.progressiveNextFraction <= 1.0001f &&
	    (quality.progressiveCeiling >= 0 ||
	     quality.progressiveNextFraction <= 1.0e-6f) &&
	    std::isfinite(quality.pointProxyPixelThreshold) &&
	    quality.pointProxyPixelThreshold >= 1.0f &&
	    quality.pointProxyPixelThreshold <= 64.01f &&
	    (std::isfinite(quality.maximumProjectedErrorPixels) ?
		quality.maximumProjectedErrorPixels >= 0.0 :
		std::isinf(quality.maximumProjectedErrorPixels) &&
		quality.maximumProjectedErrorPixels > 0.0) &&
	    quality.presentedRenderCost > 0;
    }

    static Quality sanitizeQuality(const Quality &quality)
    {
	Quality result = quality;
	result.targetPixelError = std::max(0.25f,
	    std::min(1.0f, result.targetPixelError));
	result.progressiveCeiling = result.progressiveCeiling < -1 ? -1 :
	    std::min<int>(BOBOL_MESH_LOD_CUT_COUNT_MAX - 1,
		result.progressiveCeiling);
	result.progressiveNextFraction = result.progressiveCeiling < 0 ? 0.0f :
	    std::max(0.0f, std::min(1.0f,
		result.progressiveNextFraction));
	result.pointProxyPixelThreshold = std::max(1.0f,
	    std::min(64.0f, result.pointProxyPixelThreshold));
	if (!std::isfinite(result.maximumProjectedErrorPixels) ||
	    result.maximumProjectedErrorPixels < 0.0)
	    result.maximumProjectedErrorPixels =
		std::numeric_limits<double>::infinity();
	result.provenRenderCostCapacity = 0;
	return result;
    }

    static double progressiveRank(int ceiling, float nextFraction)
    {
	return ceiling < 0 ?
	    static_cast<double>(BOBOL_MESH_LOD_CUT_COUNT_MAX) :
	    static_cast<double>(ceiling) +
		static_cast<double>(std::max(0.0f,
		    std::min(1.0f, nextFraction)));
    }

    static bool controlsDominate(const Quality &candidate,
	const Quality &current)
    {
	const bool pixelNoWorse = candidate.targetPixelError <=
	    current.targetPixelError + 1.0e-6f;
	const bool proxyNoWorse = candidate.pointProxyPixelThreshold <=
	    current.pointProxyPixelThreshold + 1.0e-6f;
	const bool ceilingNoWorse = progressiveRank(
	    candidate.progressiveCeiling,
	    candidate.progressiveNextFraction) + 1.0e-6 >= progressiveRank(
		current.progressiveCeiling,
		current.progressiveNextFraction);
	const bool strictlyBetter = candidate.targetPixelError <
		current.targetPixelError - 1.0e-6f ||
	    candidate.pointProxyPixelThreshold <
		current.pointProxyPixelThreshold - 1.0e-6f ||
	    progressiveRank(candidate.progressiveCeiling,
		candidate.progressiveNextFraction) > progressiveRank(
		current.progressiveCeiling,
		current.progressiveNextFraction) + 1.0e-6;
	return pixelNoWorse && proxyNoWorse && ceilingNoWorse &&
	    strictlyBetter;
    }

    static bool controlsEquivalent(const Quality &candidate,
	const Quality &current)
    {
	return std::fabs(candidate.targetPixelError -
		current.targetPixelError) <= 1.0e-6f &&
	    std::fabs(candidate.pointProxyPixelThreshold -
		current.pointProxyPixelThreshold) <= 1.0e-6f &&
	    std::fabs(progressiveRank(candidate.progressiveCeiling,
		candidate.progressiveNextFraction) - progressiveRank(
		current.progressiveCeiling,
		current.progressiveNextFraction)) <= 1.0e-6;
    }

    static bool preferFidelityCandidate(const Quality &candidate,
	const Quality &current)
    {
	const bool candidateHasBound =
	    std::isfinite(candidate.maximumProjectedErrorPixels);
	const bool currentHasBound =
	    std::isfinite(current.maximumProjectedErrorPixels);
	/* A measured candidate is more informative, but not intrinsically more
	 * detailed than an unmeasured one.  Accept it across that evidence
	 * boundary only when its controls are identical or conservatively dominate
	 * the current controls. */
	if (candidateHasBound != currentHasBound)
	    return controlsEquivalent(candidate, current) ||
		controlsDominate(candidate, current);
	if (candidateHasBound && std::fabs(
		candidate.maximumProjectedErrorPixels -
		current.maximumProjectedErrorPixels) > 1.0e-9)
	    return candidate.maximumProjectedErrorPixels <
		current.maximumProjectedErrorPixels;
	return controlsDominate(candidate, current);
    }

    void promote(size_t index)
    {
	if (index == 0 || index >= this->entryCount)
	    return;
	const Entry entry = this->entries[index];
	for (size_t i = index; i > 0; --i)
	    this->entries[i] = this->entries[i - 1];
	this->entries[0] = entry;
    }

    std::array<Entry, capacityValue> entries = {};
    size_t entryCount = 0;
    size_t rememberCountValue = 0;
    size_t recallCountValue = 0;
};

static_assert(std::is_trivially_copyable<
    BObolLodViewQualityHistory>::value,
    "exact-view quality history must remain an allocation-free value");

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
	HANDOFF_CONFIRMATION,
	TRIANGLE_RECOVERY
    };

    bool presentationPending(void) const
    {
	return this->phase == Phase::ADAPTIVE_CALIBRATION ||
	    this->phase == Phase::HANDOFF_CONFIRMATION;
    }

    bool adaptiveCalibrationPending(void) const
    {
	return this->phase == Phase::ADAPTIVE_CALIBRATION;
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
	if (this->phase != Phase::TRIANGLE_RECOVERY)
	    this->phase = Phase::ADAPTIVE_CALIBRATION;
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

    void completeTriangleRecovery(void)
    {
	if (this->phase == Phase::TRIANGLE_RECOVERY)
	    this->phase = Phase::IDLE;
    }

    void reset(void)
    {
	this->phase = Phase::IDLE;
    }

private:
    Phase phase = Phase::IDLE;
};

class BObolLodPointProxyEvidence {
private:
    friend class BObolLodAdmissionPlanner;

    static constexpr float maximumPixelThreshold(void)
    {
	return 64.0f;
    }

    static bool atMaximumPixelThreshold(float threshold)
    {
	return threshold >= maximumPixelThreshold() - 0.01f;
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
	    presentedMaximum == 0 ||
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
	if (this->safeThreshold > this->unsafeThreshold)
	    next = static_cast<float>(std::sqrt(
		static_cast<double>(this->safeThreshold) *
		static_cast<double>(this->unsafeThreshold)));
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

	this->unsafeThreshold = std::max(this->unsafeThreshold, current);
	if (this->safeThreshold > 0.0f &&
	    this->safeThreshold <= this->unsafeThreshold + 0.0001f)
	    this->safeThreshold = 0.0f;

	float next = current;
	if (this->safeThreshold > this->unsafeThreshold) {
	    if (this->safeThreshold / this->unsafeThreshold <= 1.08f)
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
	bool reusableSample, size_t unresolvedStructuralCount = 0)
    {
	Decision decision;
	const float current = sanitize(currentThreshold);
	decision.threshold = current;
	if (!reusableSample || !renderNanoseconds || targetFps <= 0.0f)
	    return decision;

	const long double target = 1000000000.0L /
	    static_cast<long double>(targetFps);
	if (static_cast<long double>(renderNanoseconds) > target * 1.10L) {
	    return this->interrupted(
		current, renderNanoseconds, targetFps);
	}
	/* A cheap frame containing unresolved structural boxes is not evidence
	 * that this point cut is safe.  Let the importance-ordered repair pass
	 * either satisfy coverage or produce the explicit budget-blocked witness
	 * above.  In particular, do not lower the threshold and re-expose boxes
	 * before that proof exists. */
	if (unresolvedStructuralCount)
	    return decision;

	if (this->safeThreshold <= 0.0f || current < this->safeThreshold)
	    this->safeThreshold = current;
	/* A coarse presentation cut is never terminal merely because the
	 * controller's original pressure-probe latch was consumed by an
	 * intervening publication or plan-preparation frame.  Any unchanged,
	 * reusable frame with material headroom is fresh proof that a finer cut
	 * should be tried.  The remembered safe/unsafe bracket still makes the
	 * search bounded and prevents threshold chatter near the deadline. */
	if (current <= 1.01f ||
	    static_cast<long double>(renderNanoseconds) >= target * 0.80L)
	    return decision;

	float next = current;
	if (this->unsafeThreshold > 0.0f &&
	    this->safeThreshold > this->unsafeThreshold) {
	    /* Once the bracket is within eight percent, retaining the proven
	     * safe side avoids visible threshold chatter for negligible quality. */
	    if (this->safeThreshold / this->unsafeThreshold <= 1.08f)
		return decision;
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
	const BObolLodPointProxyEvidence &value) : pointProxyValue(value) {}
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
	return this->pointProxyValue;
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
    BObolLodPointProxyEvidence pointProxyValue;
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
    enum class EvidenceAction : uint8_t {
	RESET_CAPACITY = 0,
	CLEAR_CAPACITY_LIMIT,
	RESET_CAPACITY_MEASUREMENT,
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
	const BObolLodAdmissionCursor &cursor, EvidenceAction action)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	switch (action) {
	    case EvidenceAction::RESET_CAPACITY:
		result.nextEvidence.capacityValue.reset();
		result.nextCursor.reset();
		break;
	    case EvidenceAction::CLEAR_CAPACITY_LIMIT:
		result.nextEvidence.capacityValue.clearBudgetLimit();
		break;
	    case EvidenceAction::RESET_CAPACITY_MEASUREMENT:
		result.nextEvidence.capacityValue.resetCalibration();
		break;
	    case EvidenceAction::RESET_CAPACITY_OVERLOAD:
		result.nextEvidence.capacityValue.resetOverloadRecovery();
		break;
	    case EvidenceAction::CLEAR_RETAINED_QUALITY_FLOOR:
		result.nextEvidence.capacityValue.clearRetainedQualityFloorBudget();
		break;
	    case EvidenceAction::CLEAR_RETAINED_RECOVERY_CEILING:
		result.nextEvidence.capacityValue.clearRetainedRecoveryCeiling();
		break;
	    case EvidenceAction::CLEAR_DEADLINE_CAPACITY_CEILING:
		result.nextEvidence.capacityValue.clearDeadlineCapacityCeiling();
		break;
	    case EvidenceAction::CANCEL_HEADROOM_RETRY:
		result.nextEvidence.headroomValue.cancelRetry();
		break;
	    case EvidenceAction::RESET_POINT_PROXY:
		result.nextEvidence.pointProxyValue.reset();
		break;
	    case EvidenceAction::RESET_STRUCTURAL_ADMISSION:
		result.nextEvidence.structuralValue.reset();
		break;
	    case EvidenceAction::RESET_RESOURCE_SERVICE:
		result.nextEvidence.resourcesValue.resetForServiceChange();
		break;
	    case EvidenceAction::MARK_RESOURCE_RECOVERY_HANDLED:
		result.nextEvidence.resourcesValue.markRecoveryHandled();
		break;
	}
	return result;
    }

    static BObolLodAdmissionPlan requestCapacityRescan(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.capacityValue.requestRescanAfterFrame();
	return result;
    }

    static BObolLodAdmissionPlan requestRetainedRecovery(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.capacityValue.requestRetainedRecovery(budget);
	return result;
    }

    static BObolLodAdmissionPlan requestRetainedNormalization(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.capacityValue.requestRetainedNormalization(budget);
	return result;
    }

    static BObolLodAdmissionPlan requestRetainedReallocation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	bool preserveCurrentBudget = true)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.capacityValue.requestRetainedReallocation(
	    preserveCurrentBudget);
	return result;
    }

    static BObolLodAdmissionPlan requestPresentationReconciliation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.capacityValue.requestPresentationReconciliation(budget);
	return result;
    }

    static BObolLodAdmissionPlan requestCoverageCompletion(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t activeCost,
	size_t certifiedBudget)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.transitionValue =
	    result.nextEvidence.capacityValue.requestCoverageCompletion(
		activeCost, certifiedBudget);
	return result;
    }

    static BObolLodAdmissionPlan recordDeadlineCapacityMiss(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t attemptedBudget,
	bool staticDeadline = false)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.capacityValue.noteDeadlineCapacityMiss(
	    attemptedBudget, staticDeadline);
	return result;
    }

    static BObolLodAdmissionPlan setRetainedQualityFloor(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget,
	uint64_t populationSignature, size_t activeCost,
	size_t minimumActiveCost)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.capacityValue.setRetainedQualityFloorBudget(
	    result.nextCursor, budget, populationSignature, activeCost,
	    minimumActiveCost);
	return result;
    }

    static BObolLodAdmissionPlan recordRetainedQualityFloorMiss(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.transitionChanged =
	    result.nextEvidence.capacityValue.noteRetainedQualityFloorMiss();
	return result;
    }

    static BObolLodAdmissionPlan rejectRetainedQualityFloor(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.transitionChanged =
	    result.nextEvidence.capacityValue.rejectRetainedQualityFloor();
	return result;
    }

    static BObolLodAdmissionPlan recordRetainedQualityFloorMet(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, bool exactProtectedPopulation,
	uint64_t populationSignature, size_t presentedCost)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.transitionChanged =
	    result.nextEvidence.capacityValue.noteRetainedQualityFloorMet(
		exactProtectedPopulation, populationSignature, presentedCost);
	return result;
    }

    static BObolLodAdmissionPlan confirmRetainedRecoveryPresentation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, bool onePixelReady)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.transitionChanged = result.nextEvidence.capacityValue.
	    confirmRetainedRecoveryPresentation(
		onePixelReady, result.nextCursor);
	return result;
    }

    static BObolLodAdmissionPlan raiseCapacityBudget(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, size_t budget)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.capacityValue.raiseCurrentBudget(budget);
	return result;
    }

    static BObolLodAdmissionPlan finishBlockedCapacityPass(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodCapacityEvidence::CalibrationInputs &inputs)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.calibrationDecision =
	    result.nextEvidence.capacityValue.finishBlockedPass(inputs);
	if (result.calibrationDecision.restartSubmission)
	    result.nextCursor.reset();
	return result;
    }

    static BObolLodAdmissionPlan completeCapacityCalibrationFrame(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.completedFrameDecision = result.nextEvidence.capacityValue.
	    completeCalibrationFrame(result.nextCursor);
	return result;
    }

    static BObolLodAdmissionPlan completeCapacitySearchFrame(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodCapacityEvidence::CompletedFrameInputs &inputs)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.completedFrameDecision = result.nextEvidence.capacityValue.
	    completeCapacitySearchFrame(result.nextCursor, inputs);
	return result;
    }

    static BObolLodAdmissionPlan retireUnmeasurableCapacityFrame(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.completedFrameDecision = result.nextEvidence.capacityValue.
	    retireUnmeasurableCalibrationFrame(result.nextCursor);
	return result;
    }

    static BObolLodAdmissionPlan plan(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodAdmissionRevisionStamp &revisionStamp,
	const BObolLodCapacityEvidence::Inputs &inputs)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.revisionStamp = revisionStamp;
	if (!result.nextCursor.matches(revisionStamp))
	    result.nextCursor.resetForRevision(revisionStamp);
	result.capacityDecision =
	    result.nextEvidence.capacityValue.planPass(inputs, result.nextCursor);
	return result;
    }

    static BObolLodAdmissionPlan planResourceObservation(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, bool cpuPressure,
	bool gpuPressure, bool recoveryEnabled)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.resourceDecision = result.nextEvidence.resourcesValue.observe(
	    cpuPressure, gpuPressure, recoveryEnabled);
	return result;
    }

    static BObolLodAdmissionPlan planHeadroomRetry(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
	size_t activePopulationCost, size_t currentBudget)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.headroomAccepted = result.nextEvidence.headroomValue.armRetry(
	    viewEpoch, policyEpoch, activePopulationCost, currentBudget);
	return result;
    }

    static BObolLodAdmissionPlan planHeadroomTransientDeferral(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
	size_t activePopulationCost)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.headroomAccepted =
	    result.nextEvidence.headroomValue.deferTransientReplay(
		viewEpoch, policyEpoch, activePopulationCost);
	return result;
    }

    static BObolLodAdmissionPlan planHeadroomConsumption(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
	size_t currentBudget, size_t activePopulationCost,
	long double demonstratedBudget, uint64_t elapsedNanoseconds,
	uint64_t targetNanoseconds, bool reusableSample)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.headroomAccepted = result.nextEvidence.headroomValue.consumeRetry(
	    viewEpoch, policyEpoch, currentBudget, activePopulationCost,
	    demonstratedBudget, elapsedNanoseconds, targetNanoseconds,
	    reusableSample);
	return result;
    }

    static BObolLodAdmissionPlan recordHeadroomRejection(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	BObolLodViewEpoch viewEpoch, BObolLodPolicyEpoch policyEpoch,
	size_t activePopulationCost, size_t attemptedBudget)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.nextEvidence.headroomValue.rejectRetry(viewEpoch, policyEpoch,
	    activePopulationCost, attemptedBudget);
	return result;
    }

    static BObolLodAdmissionPlan planPointStructuralDistribution(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, float currentThreshold,
	const std::array<size_t, 7> &cumulativeCount, size_t visibleCount,
	size_t maximumUncollapsedCount)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.pointProxyDecision =
	    result.nextEvidence.pointProxyValue.seedFromStructuralDistribution(
		currentThreshold, cumulativeCount, visibleCount,
		maximumUncollapsedCount);
	return result;
    }

    static BObolLodAdmissionPlan planPointInterrupted(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, float currentThreshold,
	uint64_t renderNanoseconds, float targetFps)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.pointProxyDecision = result.nextEvidence.pointProxyValue.interrupted(
	    currentThreshold, renderNanoseconds, targetFps);
	return result;
    }

    static BObolLodAdmissionPlan planPointStructuralCoverageBlocked(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, float currentThreshold,
	size_t unresolvedStructuralCount)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.pointProxyDecision =
	    result.nextEvidence.pointProxyValue.structuralCoverageBlocked(
		currentThreshold, unresolvedStructuralCount);
	return result;
    }

    static BObolLodAdmissionPlan planPointCompleted(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor, float currentThreshold,
	uint64_t renderNanoseconds, float targetFps, bool reusableSample,
	size_t unresolvedStructuralCount = 0)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.pointProxyDecision = result.nextEvidence.pointProxyValue.completed(
	    currentThreshold, renderNanoseconds, targetFps, reusableSample,
	    unresolvedStructuralCount);
	return result;
    }

    static bool pointRequiresReusableConfirmation(float currentThreshold,
	size_t unresolvedStructuralCount = 0)
    {
	return BObolLodPointProxyEvidence::requiresReusableConfirmation(
	    currentThreshold, unresolvedStructuralCount);
    }

    static bool pointAggregationApplicable(size_t visibleOccurrenceCount)
    {
	return BObolLodPointProxyEvidence::applicable(
	    visibleOccurrenceCount);
    }

    static bool pointAggregationApplicableAcrossCameraInvalidation(
	bool currentCensusApplicable, float retainedThreshold)
    {
	return BObolLodPointProxyEvidence::
	    applicableAcrossCameraInvalidation(
		currentCensusApplicable, retainedThreshold);
    }

    static bool pointStreamingPopulationWorkPending(bool submissionPending,
	bool resultsPending, bool servicePending, size_t providerPending)
    {
	return BObolLodPointProxyEvidence::streamingPopulationWorkPending(
	    submissionPending, resultsPending, servicePending,
	    providerPending);
    }

    static bool pointCalibrationYieldsToStructuralRepair(
	bool calibrationPending, bool exactOccurrenceClassification,
	size_t unresolvedStructuralCount, bool producerOwnsFutureFrame)
    {
	return BObolLodPointProxyEvidence::
	    calibrationYieldsToStructuralRepair(calibrationPending,
		exactOccurrenceClassification, unresolvedStructuralCount,
		producerOwnsFutureFrame);
    }

    static bool onePixelPresentationReady(bool haveCadAssemblies,
	bool exactFrame, bool exactOccurrenceClassification,
	size_t presentedSubpixelOccurrences,
	size_t presentedStructuralBoxes, int progressiveCeiling,
	int maximumActiveProgressiveCut,
	float pointProxyPixelThreshold, size_t managedRenderCost)
    {
	return BObolLodPointProxyEvidence::onePixelPresentationReady(
	    haveCadAssemblies, exactFrame, exactOccurrenceClassification,
	    presentedSubpixelOccurrences, presentedStructuralBoxes,
	    progressiveCeiling, maximumActiveProgressiveCut,
	    pointProxyPixelThreshold, managedRenderCost);
    }

    static bool capacitySamplePopulationReady(bool haveCadAssemblies,
	bool exactFrame, bool exactOccurrenceClassification,
	size_t presentedStructuralBoxes, bool allocationPresentationRealized)
    {
	return BObolLodPointProxyEvidence::capacitySamplePopulationReady(
	    haveCadAssemblies, exactFrame, exactOccurrenceClassification,
	    presentedStructuralBoxes, allocationPresentationRealized);
    }

    static bool pointDeadlineRequiresPopulationAggregation(
	size_t activeRenderCost, size_t minimumRenderCost,
	int presentedMaximum, int correctedCeiling,
	size_t correctedRenderCostBudget)
    {
	return BObolLodPointProxyEvidence::
	    deadlineRequiresPopulationAggregation(activeRenderCost,
		minimumRenderCost, presentedMaximum, correctedCeiling,
		correctedRenderCostBudget);
    }

    static bool pointProducerOwnsCalibrationFrame(bool submissionPending,
	bool submissionPausedByCalibration, bool providerPending,
	bool servicePending, bool publicationAwaitingFrameRequest)
    {
	return BObolLodPointProxyEvidence::producerOwnsCalibrationFrame(
	    submissionPending, submissionPausedByCalibration, providerPending,
	    servicePending, publicationAwaitingFrameRequest);
    }

    /* One completed framebuffer may advance only one measurement owner.  A
     * bounded capacity search has already frozen the point/mesh population it
     * is measuring; consuming point calibration from the same frame would
     * reset that certificate before it can collect its finite sample set.
     * The still-pending point phase supplies the successor frame after the
     * capacity search reaches a terminal decision. */
    static bool pointCalibrationOwnsCompletedFrame(
	bool pointCalibrationPending, bool capacitySamplePending)
    {
	return pointCalibrationPending && !capacitySamplePending;
    }

    static bool pointBlocksSourceAdmission(
	bool discoveryCalibrationPending, bool stableCalibrationPending)
    {
	return BObolLodPointProxyEvidence::blocksSourceAdmission(
	    discoveryCalibrationPending, stableCalibrationPending);
    }

    /* A presentation-owned measurement freezes the occurrence population it
     * is measuring.  An unchanged submission cursor must not run, or count as
     * a future geometry producer, until that frame completes.  Source and
     * inventory invalidation are handled before callers apply this gate and
     * explicitly retire stale capacity evidence. */
    static bool presentationPausesSubmission(
	bool discoveryCalibrationPending, bool stableCalibrationPending,
	bool capacitySamplePending, bool stablePresentationAvailable)
    {
	return capacitySamplePending || discoveryCalibrationPending ||
	    (stableCalibrationPending && stablePresentationAvailable);
    }

    static bool maySeedPointStructuralDistribution(
	bool discoveryCalibrationPending, bool stableCalibrationPending,
	bool triangleRecoveryPending)
    {
	return BObolLodPointProxyEvidence::maySeedStructuralDistribution(
	    discoveryCalibrationPending, stableCalibrationPending,
	    triangleRecoveryPending);
    }

    static bool pointAtMaximumPixelThreshold(float threshold)
    {
	return BObolLodPointProxyEvidence::atMaximumPixelThreshold(threshold);
    }

    static bool shouldRecoverTriangleDetail(
	bool reducibleProgressiveDetail, bool stableSampleOverloaded,
	bool coarsePointCut, bool protectedQualityOwnsCuts)
    {
	return BObolLodPointProxyEvidence::shouldRecoverTriangleDetail(
	    reducibleProgressiveDetail, stableSampleOverloaded, coarsePointCut,
	    protectedQualityOwnsCuts);
    }

    static float triangleRecoveryPointThreshold(float currentThreshold,
	bool exactStructuralPopulation, size_t structuralFallbackCount)
    {
	return BObolLodPointProxyEvidence::triangleRecoveryPointThreshold(
	    currentThreshold, exactStructuralPopulation,
	    structuralFallbackCount);
    }

    static BObolLodAdmissionPlan planStructural(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor,
	const BObolLodStructuralAdmissionEvidence::Inputs &inputs)
    {
	BObolLodAdmissionPlan result = beginPlan(evidence, cursor);
	result.structuralDecision =
	    result.nextEvidence.structuralValue.planStructural(inputs);
	return result;
    }

    static size_t unaggregatableStructuralCount(
	const std::array<size_t, 7> &cumulativeCount, size_t visibleCount)
    {
	const size_t aggregatable = std::min(
	    visibleCount, cumulativeCount.back());
	return visibleCount - aggregatable;
    }

    /* Divide only proven marginal capacity.  A zero result means the exact
     * frame did not establish even one scheduling unit per unresolved
     * occurrence, so the caller must not arm a partial repair transaction. */
    static size_t structuralPerOccurrenceReservation(size_t admittedBudget,
	size_t activeCost, size_t unresolvedCount)
    {
	if (!unresolvedCount || admittedBudget <= activeCost)
	    return 0;
	const size_t available = admittedBudget - activeCost;
	return available >= unresolvedCount ?
	    available / unresolvedCount : 0;
    }

    /* noteDeadlineCapacityMiss() records an already strict, view-local upper
     * bound: it is five percent below the population which missed the hard
     * deadline.  Marginal recovery has its own throughput safety factor
     * before reaching this helper.  Applying that factor to this bound again
     * needlessly discards a second fifth of capacity and leaves prominent
     * geometry coarse even though a smaller allocation is still admissible. */
    static size_t capBudgetAtDeadlineCeiling(size_t estimatedBudget,
	size_t deadlineCapacityCeiling)
    {
	return deadlineCapacityCeiling == SIZE_MAX ? estimatedBudget :
	    std::min(estimatedBudget, deadlineCapacityCeiling);
    }

    static bool rejectAfterInterruptedFrame(bool interactive,
	bool cadDrawAttempted, bool preparationOnlyRetry,
	bool staticTrialActive)
    {
	return !interactive && cadDrawAttempted && !preparationOnlyRetry &&
	    staticTrialActive;
    }

    /* The one-pixel point trial and its enclosing static-quality trial own one
     * candidate framebuffer.  If that frame is rejected, both owners must
     * release it together.  Retiring an unrelated point calibration would
     * instead skip a required exact-classification frame. */
    static bool rejectedOnePixelTrialReleasesCalibration(bool trialRejected,
	bool calibrationPending, float pointProxyPixelThreshold)
    {
	return trialRejected && calibrationPending &&
	    std::isfinite(pointProxyPixelThreshold) &&
	    pointProxyPixelThreshold <= 1.01f;
    }

    /* A deadline miss can close the static-quality phase only when it was an
     * explicit static-quality trial against a closed occurrence population.
     * An ordinary quiet-frame miss is capacity evidence for that exact
     * attempted population, but it must not retire a static phase which has
     * not started yet.  Doing so made a cold scene remember a transient
     * realization miss and permanently forgo the eventual event-driven 10 Hz
     * quality pass. */
    static bool terminalCapacityMiss(bool interactive,
	bool cadDrawAttempted, bool transientPresentationRetry,
	bool staticTrialActive, bool populationWorkPending)
    {
	return !interactive && cadDrawAttempted &&
	    !transientPresentationRetry && staticTrialActive &&
	    !populationWorkPending;
    }

    static bool ordinaryHeadroomAllowed(bool staticPhaseBlocked)
    {
	return !staticPhaseBlocked;
    }

    /* A bounded static phase may improve an unaffordable atomic quality floor
     * incrementally.  The same marginal path remains available after either
     * the complete floor or a richer presentation has been rejected, but it
     * must never leak into ordinary quiet convergence. */
    static bool marginalStaticCapacityAllowed(bool staticPhaseActive,
	bool staticPhaseRejected, bool protectedFloorRejected)
    {
	return staticPhaseActive || staticPhaseRejected ||
	    protectedFloorRejected;
    }

    /* Retained occurrence cuts and the renderer's submitted cut are separate
     * quality dimensions.  A single progressive object can have every
     * requested suffix resident while a motion/quiet-cadence ceiling still
     * hides most of it.  Treat that reversible presentation gap as actionable
     * static quality debt; otherwise provider satisfaction falsely terminates
     * the event-driven quality pass after a pose-only interaction. */
    static bool actionableProgressiveQualityDebt(
	size_t activePayloadCount, size_t satisfiedPayloadCount,
	size_t memoryLimitedPayloadCount, int activeMaximumCut,
	int presentationCeiling)
    {
	const size_t unsatisfiedPayloadCount =
	    activePayloadCount > satisfiedPayloadCount ?
		activePayloadCount - satisfiedPayloadCount : 0;
	const bool residentDemandDebt =
	    unsatisfiedPayloadCount > memoryLimitedPayloadCount;
	const bool presentationDebt = progressivePresentationQualityDebt(
	    activeMaximumCut, presentationCeiling);
	return residentDemandDebt || presentationDebt;
    }

    static bool progressivePresentationQualityDebt(int activeMaximumCut,
	int presentationCeiling)
    {
	return presentationCeiling >= 0 &&
	    activeMaximumCut > presentationCeiling;
    }

    /* A renderer-wide PoP ceiling is a transient capacity guard for a
	* multi-occurrence scene.  Once no live static staircase owns it (or that
	* staircase has rejected its next step), an exact quiet frame must be
	* translated into occurrence-local cuts before readiness.  For exactly one
	* progressive occurrence the ordinal is already occurrence-local, so the
	* proven ceiling is a valid terminal representation and avoids a redundant
	* allocator/repaint transaction. */
    static bool terminalGlobalCeilingRequiresReconciliation(
	bool stableTerminalContext, bool exactCompletedFrame,
	int progressiveCeiling, bool staticPhaseActive,
	bool staticPhaseRejected, size_t progressivePayloadCount)
    {
	return stableTerminalContext && exactCompletedFrame &&
	    progressiveCeiling >= 0 &&
	    progressivePayloadCount != 1 &&
	    (!staticPhaseActive || staticPhaseRejected);
    }

    /* Translate a renderer-wide CAD cut into the retained allocator's scene
     * currency.  Backend submitted-work records may expand indexed vertices,
     * normals, and triangles differently (OSMesa flat shading is the important
     * case), so they are not a valid occurrence-allocation budget.  Preserve
     * the scene's non-CAD floor and replace only its retained CAD component with
     * the canonical cost table for the exact presented ceiling. */
    static size_t canonicalSceneCostAtCadCeiling(size_t activeSceneCost,
	size_t activeCadCost, size_t cadCeilingCost)
    {
	const size_t nonCadCost = activeSceneCost > activeCadCost ?
	    activeSceneCost - activeCadCost : 0;
	return cadCeilingCost > SIZE_MAX - nonCadCost ? SIZE_MAX :
	    nonCadCost + cadCeilingCost;
    }

    /* A retained allocation certificate is the exact post-transaction
     * population.  The admission cursor deliberately freezes the pre-pass
     * scene cost for bounded traversal accounting, so using that older value
     * as the capacity candidate makes the completed frame look stale even
     * when no cut changed. */
    static size_t capacitySearchPresentedCost(bool currentAllocation,
	size_t frozenActiveCost, size_t selectedPresentationCost)
    {
	return currentAllocation && selectedPresentationCost > 0 ?
	    selectedPresentationCost : frozenActiveCost;
    }

    static bool onePixelTrialRequired(float pointProxyPixelThreshold)
    {
	return std::isfinite(pointProxyPixelThreshold) &&
	    pointProxyPixelThreshold > 1.01f;
    }

    /* A coarse point cut may hide either resident meshes or structural
     * fallback boxes.  Static headroom may expose the former under its hard
     * deadline, but exposing the latter is a visual regression, not quality
     * refinement.  Retain aggregation whenever structural records remain,
     * independent of whether a renderer-wide PoP ceiling is also installed. */
    static bool retainAggregatedPresentation(int progressiveCeiling,
	float pointProxyPixelThreshold, size_t activePayloadCount,
	bool structuralFallbackPopulation)
    {
	return std::isfinite(pointProxyPixelThreshold) &&
	    pointProxyPixelThreshold > 1.01f && activePayloadCount > 1 &&
	    (progressiveCeiling >= 0 || structuralFallbackPopulation);
    }

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
	uint64_t staticDeadlineNanoseconds)
    {
	return !interactive && exactReusableFrame && !producerWork &&
	    !resourcePressure && !structuralFallbackPopulation &&
	    !staticPhaseBlocked &&
	    onePixelTrialRequired(pointProxyPixelThreshold) &&
	    renderNanoseconds > 0 && preferredFrameNanoseconds > 0 &&
	    staticDeadlineNanoseconds > preferredFrameNanoseconds &&
	    renderNanoseconds <= staticDeadlineNanoseconds;
    }

    /* A quiet exact one-pixel frame which misses the preferred streaming
     * cadence but fits the independently interruptible static deadline has
     * already proved the transition which a coarse point frame would need to
     * trial.  Enter the event-driven quality phase directly.  Sending this
     * frame through ordinary point/triangle overload recovery first coarsens
     * it to the 20 Hz budget; the later static pass then restores the same
     * population at 10 Hz, producing a visible and potentially unbounded
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
	uint64_t staticDeadlineNanoseconds)
    {
	return !interactive && exactReusableFrame && !producerWork &&
	    !resourcePressure && !structuralFallbackPopulation &&
	    !handoffPending && !staticPhaseBlocked &&
	    progressiveCeiling < 0 &&
	    std::isfinite(pointProxyPixelThreshold) &&
	    pointProxyPixelThreshold <= 1.01f &&
	    reducibleProgressiveDetail && renderNanoseconds > 0 &&
	    preferredFrameNanoseconds > 0 &&
	    staticDeadlineNanoseconds > preferredFrameNanoseconds &&
	    renderNanoseconds > preferredFrameNanoseconds &&
	    renderNanoseconds <= staticDeadlineNanoseconds;
    }

    /* A global PoP ordinal is occurrence-local only when the scene contains
	* exactly one progressive payload.  In that case a one-cut completed-frame
	* probe bounds a discrete Lucy-like jump without distorting relative scene
	* importance.  A many-part scene must return -1 and spend static headroom
	* through the occurrence allocator; walking a scene-wide ordinal makes all
	* parts equally coarse and starves prominent shapes. */
    static int stagedProgressiveCeiling(int currentCeiling,
	    int activeMaximum, float pointProxyPixelThreshold,
	    size_t activePayloadCount)
    {
	if (currentCeiling < 0 || activeMaximum <= currentCeiling ||
	    activePayloadCount != 1 ||
	    !std::isfinite(pointProxyPixelThreshold) ||
	    pointProxyPixelThreshold <= 0.0f)
	    return -1;
	return std::min(activeMaximum, currentCeiling + 1);
    }

    /* A single progressive occurrence can preserve its last exact cut while
     * a richer occurrence allocation is prepared behind it.  Without this
     * baseline ceiling, publication exposes the richest resident cut before
     * the completed-frame staircase can measure any neighboring cut.  A
     * global ordinal is not a valid relative-quality policy for multiple
     * occurrences, so those scenes continue to rely on the allocator. */
    static int staticProgressiveBaselineCeiling(int currentCeiling,
	    int activeMaximum, size_t activePayloadCount)
    {
	return currentCeiling < 0 && activeMaximum >= 0 &&
	    activePayloadCount == 1 ? activeMaximum : -1;
    }

    static bool staticOrdinalOverscanApplicable(size_t activePayloadCount)
    {
	return activePayloadCount == 1;
    }

    static bool retainedPointAggregationApplicable(
	bool logicalPopulationApplicable, size_t pointProxyCandidateCount)
    {
	return logicalPopulationApplicable && pointProxyCandidateCount > 0;
    }

    static bool structuralPointAggregationRequired(
	size_t logicalOccurrenceCount, size_t affordableOccurrenceCount)
    {
	return affordableOccurrenceCount != SIZE_MAX &&
	    logicalOccurrenceCount > affordableOccurrenceCount;
    }

    static bool pointPolicyExhaustedForStructuralFrontier(
	bool maximumThreshold, bool pointAggregationRequired,
	bool exactProjection, size_t unresolvedCount,
	size_t unaggregatableCount)
    {
	if (!unresolvedCount)
	    return false;
	return maximumThreshold || !pointAggregationRequired ||
	    (exactProjection && unaggregatableCount >= unresolvedCount);
    }

    static bool pointSuccessorRequiredForStructuralFrontier(
	bool exactProjection, size_t unresolvedCount,
	bool pointPolicyExhausted, bool anotherOwnerPending)
    {
	/* A successful terminal presentation may not retain structural boxes.
	 * If an exact quiet frame still has an aggregatable box frontier, the
	 * point threshold owns the only successor which can reduce that frontier
	 * without reopening a whole-scene provider pass. */
	return exactProjection && unresolvedCount > 0 &&
	    !pointPolicyExhausted && !anotherOwnerPending;
    }

    static size_t structuralFirstWaveOccurrenceLimit(size_t sceneBudget)
    {
	static constexpr size_t occurrenceCostFloor = 64;
	static constexpr size_t minimumOccurrences = 512;
	static constexpr size_t maximumOccurrences = 8192;
	if (sceneBudget == SIZE_MAX)
	    return SIZE_MAX;
	return std::max(minimumOccurrences,
	    std::min(maximumOccurrences, sceneBudget / occurrenceCostFloor));
    }

private:
    static BObolLodAdmissionPlan beginPlan(
	const BObolLodAdmissionEvidence &evidence,
	const BObolLodAdmissionCursor &cursor)
    {
	BObolLodAdmissionPlan result;
	result.nextEvidence = evidence;
	result.nextCursor = cursor;
	return result;
    }
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

#endif /* LIBBOBOL_LOD_COORDINATOR_PRIVATE_H */
