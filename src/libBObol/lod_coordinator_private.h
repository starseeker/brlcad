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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>
#include <unordered_map>
#include <vector>

/*
 * Use the same strong epoch domains as requests and results.  Keeping a
 * second coordinator-only implementation previously allowed the state
 * machine and service boundary to drift while appearing type safe.
 */
using BObolLodViewEpoch = BObolViewEpoch;
using BObolLodPolicyEpoch = BObolPolicyEpoch;
using BObolLodInventoryEpoch = BObolInventoryEpoch;
using BObolLodSourceRoutingId = BObolSourceRoutingId;

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
 * Bounded proof for a late calibrated-headroom admission retry.  Normal probe
 * counts remain a short-horizon mechanism; this witness covers the ordering
 * where reusable timing becomes available only after those probes ended.
 *
 * The first retry is unique for (view, policy, active population).  A later
 * retry for that same immutable population is allowed only when calibration
 * has raised the finite render budget by more than 25 percent above the
 * largest budget already witnessed.  The monotonic high-water mark permits a
 * newly learned capacity result to refine a stranded cut without making the
 * floating budget itself an oscillating identity.  Consequently repeated
 * retries for an unchanged population are logarithmically bounded.
 */
class BObolLodHeadroomPolicy {
public:
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
	const bool witnessedPopulation =
	    this->viewRevision == viewEpoch.value() &&
		this->policyRevision == policyEpoch.value() &&
		this->populationCost == activePopulationCost;
	if (witnessedPopulation &&
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
	const bool samePopulation =
	    this->viewRevision == viewEpoch.value() &&
	    this->policyRevision == policyEpoch.value() &&
	    this->populationCost == activePopulationCost;
	this->viewRevision = viewEpoch.value();
	this->policyRevision = policyEpoch.value();
	this->populationCost = activePopulationCost;
	const size_t observedBudget = currentBudget > witnessedBudget ?
	    currentBudget : witnessedBudget;
	this->budgetHighWater = samePopulation &&
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

    void cancelRetry(void)
    {
	this->pendingValue = false;
	this->pendingViewRevision = 0;
	this->pendingPolicyRevision = 0;
	this->pendingPopulationCost = SIZE_MAX;
	this->pendingBudget = 0;
	this->pendingTransientReplayCount = 0;
    }

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
     * presentation.  Once that proof completes, only an ordinary coverage
     * pass needs the deferred-quality successor: a retained-allocation pass
     * is already the authoritative scene-wide successor and must be allowed
     * to reach budget/handoff completion. */
    static bool needsDeferredQualitySuccessor(bool coverageCompleted,
	bool coverageMissing, bool retainedRefinementPending,
	bool retainedAllocationCompleted)
    {
	return coverageCompleted && !coverageMissing &&
	    retainedRefinementPending && !retainedAllocationCompleted;
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

    void complete(SourceKey source, size_t visibleCount)
    {
	Source &entry = this->source(source);
	entry.visibleCount = visibleCount;
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
	const size_t unresolvedStructuralBoxes =
	    inputs.presentedStructuralBoxCount >
		inputs.terminalOccurrenceFailureCount ?
	    inputs.presentedStructuralBoxCount -
		inputs.terminalOccurrenceFailureCount : 0;
	decision.visualPending =
	    structuralPending || inputs.structuralDiscovery ||
	    inputs.sourcePreparationPending ||
	    inputs.submissionPending ||
	    inputs.resultPending || inputs.publicationPending ||
	    inputs.pendingTasks > 0 || inputs.inFlight > 0 ||
	    inputs.calibrationPending || unresolvedStructuralBoxes > 0;
	decision.viewReady =
	    !decision.visualPending && !inputs.interactive;
	decision.backgroundPending =
	    inputs.queuedCacheWrites > 0 || inputs.compactionPending ||
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

	if (inputs.failedSourceCount > 0 ||
	    inputs.terminalOccurrenceFailureCount > 0)
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
class BObolLodBudgetPolicy {
public:
    struct Inputs {
	size_t activeCost = 0;
	size_t minimumActiveCost = 0;
	float targetFps = 0.0f;
	long double calibratedCostPerSecond = 0.0L;
	uint64_t observedStableNanoseconds = 0;
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
	uint64_t targetNanoseconds = 0;
	uint64_t observedNanoseconds = 0;
	long double calibratedBudget = 0.0L;
	bool interactive = false;
	bool stablePresentationHandoff = false;
	bool passAdmittedWork = false;
	/* The blocked cut came from a complete scene-wide minimax allocation.
	 * If calibration demonstrates a different allowance, its successor must
	 * recompute that allocation rather than letting ordinary per-occurrence
	 * refinement consume stale allocation stamps. */
	bool retainedAllocation = false;
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
	this->singleOccurrenceBootstrapValue =
	    inputs.coldSingleOccurrence && inputs.activeCost == 0 &&
	    inputs.calibratedCostPerSecond <= 0.0L;
	const bool requestedRetainedRecovery =
	    this->requestedRetainedRecoveryBudgetValue != SIZE_MAX &&
	    inputs.activeCost > this->requestedRetainedRecoveryBudgetValue;
	const bool requestedRetainedReallocation =
	    this->requestedRetainedReallocationValue &&
	    !inputs.interactive && !inputs.forceTerminal &&
	    inputs.activeCost > 0;
	const bool requestedPresentationReconciliation =
	    requestedRetainedReallocation &&
	    this->requestedPresentationReconciliationBudgetValue != SIZE_MAX;
	const size_t presentationReconciliationBudget =
	    this->requestedPresentationReconciliationBudgetValue;
	const size_t retainedReallocationBudget = this->currentBudgetValue;

	const long double targetNanoseconds = inputs.targetFps > 0.0f ?
	    1000000000.0L / static_cast<long double>(inputs.targetFps) : 0.0L;
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
	    this->requestedRetainedReallocationPreserveBudgetValue)
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
	/* A completed constrained framebuffer is a hard presentation-capacity
	 * witness, not another soft quality preference.  Reconcile the hidden
	 * retained prefixes at exactly that scene-wide allowance: a stale
	 * throughput seed may not undercut proven work, and a protected-quality
	 * floor may not exceed the frame deadline which created the handoff. */
	if (requestedPresentationReconciliation)
	    costBudget = presentationReconciliationBudget;

	this->currentBudgetValue = costBudget;
	/* Freeze the exact cost currencies which initialized this bounded pass.
	 * A 150k-occurrence pass may require hundreds of owner-thread windows;
	 * rediscovering the complete population before every window turns an O(N)
	 * allocation into O(N * windows) work and monopolizes the GUI thread.
	 * Occurrence changes are charged by the carried refinement remainders; a
	 * completed pass resets this snapshot before another census. */
	this->passActiveCostValue = inputs.activeCost;
	this->passMinimumActiveCostValue = inputs.minimumActiveCost;
	this->refinementRemainingValue = costBudget == SIZE_MAX ? SIZE_MAX :
	    (costBudget > inputs.activeCost ?
		costBudget - inputs.activeCost : 0);
	/* A presentation handoff removes a reversible renderer ceiling; it does
	 * not authorize rewriting retained occurrence cuts.  In particular, a
	 * zoom must begin at the last coherent cut and adjust from there instead
	 * of normalizing every occurrence to a cheap baseline.  Cold coverage and
	 * point-proxy recovery request retained admission explicitly, while a
	 * completed over-budget quiet frame uses overloadRecovery. */
	this->retainedAdmissionValue =
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
	this->retainedAdmissionRemainingValue =
	    this->retainedAdmissionValue ?
		(costBudget > inputs.minimumActiveCost ?
		    costBudget - inputs.minimumActiveCost : 0) : SIZE_MAX;
	this->passInitializedValue = true;
	this->requestedRetainedRecoveryBudgetValue = SIZE_MAX;
	this->requestedRetainedReallocationValue = false;
	this->requestedRetainedReallocationPreserveBudgetValue = true;
	this->requestedPresentationReconciliationBudgetValue = SIZE_MAX;

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
	/* Three unchanged presentations are enough to reject a one-frame setup
	 * or compositor outlier.  Replaying the same cut twelve times does not
	 * add new capacity information; on very large scenes it needlessly holds
	 * convergence behind repeated full presentation work.  The resulting
	 * calibrated budget is retained and the next admission pass supplies the
	 * next distinct sample. */
	static constexpr unsigned int minimumProbes = 3;
	static constexpr unsigned int maximumProbes = minimumProbes;
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
	this->reallocateAfterCalibrationValue =
	    this->rescanAfterFrameValue && inputs.retainedAllocation;
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
	if (this->reallocateAfterCalibrationValue) {
	    this->requestedRetainedReallocationValue = true;
	    this->requestedRetainedReallocationPreserveBudgetValue = false;
	}
	this->reallocateAfterCalibrationValue = false;
	this->resetPass();
	decision.restartSubmission = true;
	return decision;
    }

    /* A calibration barrier can outlive the population it was meant to
     * measure (most notably when a structural-proxy coverage wave temporarily
     * has no active managed mesh cost).  Such a frame is presentation-only:
     * replaying it cannot produce capacity evidence.  Retire the obsolete
     * probe series and return directly to admission instead of waiting on an
     * impossible sample forever. */
    CompletedFrameDecision retireUnmeasurableCalibrationFrame(void)
    {
	CompletedFrameDecision decision;
	if (!this->rescanAfterFrameValue)
	    return decision;
	this->calibrationFramesRemainingValue = 0;
	decision = this->completeCalibrationFrame();
	this->resetProbeSeries();
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
	this->singleOccurrenceBootstrapValue = false;
	this->passActiveCostValue = 0;
	this->passMinimumActiveCostValue = 0;
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
	this->requestedRetainedRecoveryBudgetValue = SIZE_MAX;
	this->requestedRetainedReallocationValue = false;
	this->requestedRetainedReallocationPreserveBudgetValue = true;
	this->requestedPresentationReconciliationBudgetValue = SIZE_MAX;
	this->reallocateAfterCalibrationValue = false;
	this->retainedRecoveryCeilingValue = SIZE_MAX;
	this->retainedQualityFloorBudgetValue = 0;
	this->retainedQualityFloorSignatureValue = 0;
	this->retainedQualityFloorRejectedValue = false;
	this->retainedQualityFloorMissCountValue = 0;
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
	this->requestedRetainedReallocationValue = true;
	this->requestedRetainedReallocationPreserveBudgetValue =
	    preserveCurrentBudget;
    }

    /* Convert a completed renderer-wide constrained frame into one atomic
     * occurrence-local importance allocation.  Unlike normalization, this
     * must run even when point aggregation makes the exact presented cost
     * cheaper than the hidden retained prefixes. */
    void requestPresentationReconciliation(size_t budget)
    {
	if (!budget || budget == SIZE_MAX)
	    return;
	this->requestedRetainedReallocationValue = true;
	this->requestedRetainedReallocationPreserveBudgetValue = false;
	this->requestedPresentationReconciliationBudgetValue = std::min(
	    this->requestedPresentationReconciliationBudgetValue, budget);
	this->currentBudgetValue = budget;
    }

    /* A view-significance floor may deliberately trade some of the 20 Hz
     * stable target for recognizable prominent geometry, but only after the
     * controller proves that population fits its separately bounded hard
     * quality-frame allowance.  Update all three currencies of the already
     * initialized retained pass together; changing only currentBudgetValue
     * would leave its bounded actions consuming the former upgrade limit. */
    void setRetainedQualityFloorBudget(size_t budget,
	uint64_t populationSignature, size_t activeCost,
	size_t minimumActiveCost)
    {
	if (this->retainedQualityFloorRejectedValue)
	    return;
	const size_t nextBudget = budget == SIZE_MAX ? 0 : budget;
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
	if (!this->passInitializedValue || !this->retainedAdmissionValue ||
	    !budget || budget == SIZE_MAX || budget <= this->currentBudgetValue)
	    return;
	this->currentBudgetValue = budget;
	this->refinementRemainingValue = budget > activeCost ?
	    budget - activeCost : 0;
	this->retainedAdmissionRemainingValue =
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

    void noteRetainedQualityFloorMet(void)
    {
	if (this->retainedQualityFloorBudgetValue > 0 &&
	    !this->retainedQualityFloorRejectedValue)
	    this->retainedQualityFloorMissCountValue = 0;
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

    /* A measured recovery ceiling protects the first coherent one-pixel
     * population from immediately re-admitting the cut which just missed its
     * deadline.  It must end once that population is actually ready for
     * presentation.  Keep this transition in the scalar policy so callers
     * cannot accidentally clear only the pass state, or retain the ceiling
     * forever when returning from a coarser point cut needs an extra frame. */
    bool confirmRetainedRecoveryPresentation(bool onePixelReady)
    {
	if (!onePixelReady ||
	    this->retainedRecoveryCeilingValue == SIZE_MAX)
	    return false;
	this->clearRetainedRecoveryCeiling();
	this->resetPass();
	return true;
    }

    void raiseCurrentBudget(size_t budget)
    {
	size_t raisedBudget = std::max(this->currentBudgetValue, budget);
	if (this->retainedRecoveryCeilingValue != SIZE_MAX)
	    raisedBudget = std::min(
		raisedBudget, this->retainedRecoveryCeilingValue);
	this->currentBudgetValue = raisedBudget;
    }

    size_t seedBudget(void) const { return this->seedBudgetValue; }
    static constexpr size_t singleOccurrenceBootstrapBudget(void)
    {
	/* About one medium software-rendered mesh frame.  The endpoint's 100 ms
	 * hard deadline and normal completed-frame calibration remain the actual
	 * safety contract; this value only avoids 50k-face level walking before
	 * the first useful timing sample exists. */
	return 500000;
    }
    size_t currentBudget(void) const { return this->currentBudgetValue; }
    size_t passActiveCost(void) const { return this->passActiveCostValue; }
    size_t passMinimumActiveCost(void) const
    {
	return this->passMinimumActiveCostValue;
    }
    bool retainedRecoveryCeilingActive(void) const
    {
	return this->retainedRecoveryCeilingValue != SIZE_MAX;
    }
    bool passInitialized(void) const { return this->passInitializedValue; }
    bool singleOccurrenceBootstrap(void) const
    {
	return this->singleOccurrenceBootstrapValue;
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
    size_t passActiveCostValue = 0;
    size_t passMinimumActiveCostValue = 0;
    size_t refinementRemainingValue = SIZE_MAX;
    size_t retainedAdmissionRemainingValue = SIZE_MAX;
    size_t overloadRecoveryActiveCostValue = 0;
    size_t requestedRetainedRecoveryBudgetValue = SIZE_MAX;
    size_t retainedRecoveryCeilingValue = SIZE_MAX;
    size_t retainedQualityFloorBudgetValue = 0;
    uint64_t retainedQualityFloorSignatureValue = 0;
    unsigned int retainedQualityFloorMissCountValue = 0;
    size_t probeActiveCostValue = 0;
    unsigned int probeCountValue = 0;
    unsigned int calibrationFramesRemainingValue = 0;
    bool passInitializedValue = false;
    bool singleOccurrenceBootstrapValue = false;
    bool requestedRetainedReallocationValue = false;
    bool requestedRetainedReallocationPreserveBudgetValue = true;
    size_t requestedPresentationReconciliationBudgetValue = SIZE_MAX;
    bool reallocateAfterCalibrationValue = false;
    bool retainedAdmissionValue = false;
    bool overloadRecoveryPerformedValue = false;
    bool rescanAfterFrameValue = false;
    bool stableBudgetLimitedValue = false;
    bool retainedQualityFloorRejectedValue = false;
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
	float currentPointProxyPixelThreshold = 1.0f;
    };

    struct RestoreDecision {
	bool apply = false;
	bool restoredPriorStable = false;
	bool restoredProvenQuality = false;
	bool startHandoff = false;
	bool clearPresentationLimits = false;
	float targetPixelError = 1.0f;
	int progressiveCeiling = -1;
	float pointProxyPixelThreshold = 1.0f;
	size_t provenRenderCostFloor = 0;
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
	bool preservePresentationLimits = false;
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

    void capturePrior(float targetPixelError, int progressiveCeiling,
	float pointProxyPixelThreshold, const Population &population,
	BObolLodViewEpoch viewEpoch)
    {
	this->priorStableValue = makeSnapshot(
	    targetPixelError, progressiveCeiling,
	    pointProxyPixelThreshold, population,
	    viewEpoch, 0);
	this->priorRestoredValue = false;
	this->provenQualityValue = Snapshot();
	this->handoffActiveValue = false;
	this->handoffPresentationRequiredValue = false;
	this->handoffPreservePresentationLimitsValue = false;
	this->handoffChangedCutValue = false;
	this->handoffReconciledValue = false;
	this->handoffReconciliationBudgetValue = 0;
	this->handoffReconciliationBudgetLimitValue = 0;
	this->handoffCostFloorValue = 0;
    }

    void noteStableQuality(float targetPixelError, int progressiveCeiling,
	float pointProxyPixelThreshold, const Population &population,
	BObolLodViewEpoch viewEpoch, size_t presentedRenderCost)
    {
	if (!presentedRenderCost)
	    return;
	const Snapshot candidate = makeSnapshot(
	    targetPixelError, progressiveCeiling,
	    pointProxyPixelThreshold, population,
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
	decision.targetPixelError = sanitizePixelError(
	    inputs.currentTargetPixelError);
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
		if (std::fabs(sanitizePixelError(
			inputs.currentTargetPixelError) -
			this->priorStableValue.targetPixelError) > 1.0e-6f ||
		    inputs.currentProgressiveCeiling !=
			this->priorStableValue.progressiveCeiling ||
		    std::fabs(sanitizeThreshold(
			inputs.currentPointProxyPixelThreshold) -
			this->priorStableValue.pointProxyPixelThreshold) > 1.0e-6f) {
		    decision.apply = true;
		    decision.targetPixelError =
			this->priorStableValue.targetPixelError;
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
	decision.targetPixelError = this->priorStableValue.targetPixelError;
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
	    decision.targetPixelError =
		this->provenQualityValue.targetPixelError;
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
	    decision.targetPixelError = 1.0f;
	    decision.progressiveCeiling = -1;
	    decision.pointProxyPixelThreshold = 1.0f;
	    this->handoffActiveValue = false;
	    this->handoffPresentationRequiredValue = false;
	    this->handoffPreservePresentationLimitsValue = false;
	    this->handoffChangedCutValue = false;
	    this->handoffReconciledValue = false;
	    this->handoffReconciliationBudgetValue = 0;
	    this->handoffReconciliationBudgetLimitValue = 0;
	    this->handoffCostFloorValue = 0;
	} else {
	    const bool constrained = decision.progressiveCeiling >= 0 ||
		decision.pointProxyPixelThreshold > 1.01f;
	    this->handoffActiveValue =
		!decision.restoredPriorStable && constrained;
	    this->handoffPreservePresentationLimitsValue = false;
	    this->handoffChangedCutValue = false;
	    this->handoffReconciledValue = false;
	    this->handoffReconciliationBudgetValue = 0;
	    this->handoffReconciliationBudgetLimitValue = 0;
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
	if (this->handoffActiveValue && inputs.completed && inputs.changedCut)
	    this->handoffChangedCutValue = true;
	if (this->handoffActiveValue && inputs.completed &&
	    inputs.retainedAllocationCompleted &&
	    inputs.populationQuiescent &&
	    inputs.presentationLimitsReconciled)
	    this->handoffReconciledValue = true;
	if (!inputs.completed ||
	    inputs.submissionPending || inputs.rescanAfterFrame ||
	    inputs.changedCut ||
	    (inputs.retainedRefinementBudgetBlocked &&
	     !this->handoffActiveValue) ||
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
	const bool preservePresentationLimits = false;
	this->handoffActiveValue = false;
	this->handoffPresentationRequiredValue = false;
	this->handoffPreservePresentationLimitsValue = false;
	this->handoffChangedCutValue = false;
	this->handoffReconciledValue = false;
	this->handoffReconciliationBudgetValue = 0;
	this->handoffReconciliationBudgetLimitValue = 0;
	this->handoffCostFloorValue = 0;
	decision.finishHandoff = true;
	decision.preservePresentationLimits = preservePresentationLimits;
	decision.requestRetainedRescan =
	    inputs.retainedRefinementPending &&
	    !inputs.retainedRefinementBudgetBlocked &&
	    !preservePresentationLimits;
	return decision;
    }

    void armHandoff(bool presentationRequired, size_t provenRenderCost = 0,
	size_t reconciliationBudgetLimit = 0)
    {
	this->handoffActiveValue = true;
	this->handoffPresentationRequiredValue = presentationRequired;
	/* A presentation-required handoff follows a hard frame deadline.  Its
	 * first completed constrained frame supplies the cost proof below; do not
	 * immediately remove that cut and retry the frame which just exceeded the
	 * deadline.  A static predictor may instead arm an already-presented
	 * handoff with its exact cost floor and no additional presentation gate. */
	this->handoffPreservePresentationLimitsValue = presentationRequired;
	this->handoffChangedCutValue = false;
	this->handoffReconciledValue = false;
	this->handoffReconciliationBudgetValue = 0;
	this->handoffReconciliationBudgetLimitValue =
	    presentationRequired ? reconciliationBudgetLimit : 0;
	this->handoffCostFloorValue = provenRenderCost;
    }
    bool noteFramePresented(size_t provenRenderCost = 0,
	size_t reconciliationBudget = 0)
    {
	const bool released = this->handoffActiveValue &&
	    this->handoffPresentationRequiredValue;
	if (this->handoffActiveValue) {
	    this->handoffPresentationRequiredValue = false;
	    if (provenRenderCost)
		this->handoffCostFloorValue = std::max(
		    this->handoffCostFloorValue, provenRenderCost);
	    /* A hard-deadline handoff may retire its global renderer limit only
	     * after the retained allocator has encoded an equivalent bounded
	     * presentation in occurrence-local cuts.  Keep this exact-frame
	     * capacity witness distinct from the affordable-work floor above. */
	    if (released && this->handoffPreservePresentationLimitsValue &&
		reconciliationBudget) {
		/* The constrained frame measures only the renderer-limited cut.  A
		 * preceding failed full frame may already have computed a stricter
		 * budget in the retained-population cost domain.  Never let fast
		 * timing of the coarse fallback erase that safety proof. */
		const size_t boundedBudget =
		    this->handoffReconciliationBudgetLimitValue > 0 ?
			std::min(reconciliationBudget,
			    this->handoffReconciliationBudgetLimitValue) :
			reconciliationBudget;
		this->handoffReconciliationBudgetValue = boundedBudget;
	    }
	}
	return released;
    }
    void cancelHandoff(void)
    {
	this->handoffActiveValue = false;
	this->handoffPresentationRequiredValue = false;
	this->handoffPreservePresentationLimitsValue = false;
	this->handoffChangedCutValue = false;
	this->handoffReconciledValue = false;
	this->handoffReconciliationBudgetValue = 0;
	this->handoffReconciliationBudgetLimitValue = 0;
	this->handoffCostFloorValue = 0;
    }

    void viewInvalidated(void)
    {
	this->handoffActiveValue = false;
	this->handoffPresentationRequiredValue = false;
	this->handoffPreservePresentationLimitsValue = false;
	this->handoffChangedCutValue = false;
	this->handoffReconciledValue = false;
	this->handoffReconciliationBudgetValue = 0;
	this->handoffReconciliationBudgetLimitValue = 0;
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
	this->handoffPreservePresentationLimitsValue = false;
	this->handoffChangedCutValue = false;
	this->handoffReconciledValue = false;
	this->handoffReconciliationBudgetValue = 0;
	this->handoffReconciliationBudgetLimitValue = 0;
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
    size_t handoffReconciliationBudget(void) const
    {
	return this->handoffReconciliationBudgetValue;
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

    static Snapshot makeSnapshot(float targetPixelError,
	int progressiveCeiling,
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
	snapshot.pointProxyPixelThreshold =
	    sanitizeThreshold(pointProxyPixelThreshold);
	snapshot.population = population;
	snapshot.presentedRenderCost = presentedRenderCost;
	return snapshot;
    }

    static bool populationMatches(const Population &snapshot,
	const Population &current)
    {
	return snapshot.available == current.available &&
	    snapshot.sceneDomainRevision == current.sceneDomainRevision;
    }

    Snapshot priorStableValue;
    Snapshot provenQualityValue;
    bool priorRestoredValue = false;
    bool handoffActiveValue = false;
    bool handoffPresentationRequiredValue = false;
    bool handoffPreservePresentationLimitsValue = false;
    bool handoffChangedCutValue = false;
    bool handoffReconciledValue = false;
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
	result.pointProxyPixelThreshold = std::max(1.0f,
	    std::min(64.0f, result.pointProxyPixelThreshold));
	if (!std::isfinite(result.maximumProjectedErrorPixels) ||
	    result.maximumProjectedErrorPixels < 0.0)
	    result.maximumProjectedErrorPixels =
		std::numeric_limits<double>::infinity();
	result.provenRenderCostCapacity = 0;
	return result;
    }

    static int progressiveRank(int ceiling)
    {
	return ceiling < 0 ? BOBOL_MESH_LOD_CUT_COUNT_MAX : ceiling;
    }

    static bool controlsDominate(const Quality &candidate,
	const Quality &current)
    {
	const bool pixelNoWorse = candidate.targetPixelError <=
	    current.targetPixelError + 1.0e-6f;
	const bool proxyNoWorse = candidate.pointProxyPixelThreshold <=
	    current.pointProxyPixelThreshold + 1.0e-6f;
	const bool ceilingNoWorse = progressiveRank(
	    candidate.progressiveCeiling) >= progressiveRank(
		current.progressiveCeiling);
	const bool strictlyBetter = candidate.targetPixelError <
		current.targetPixelError - 1.0e-6f ||
	    candidate.pointProxyPixelThreshold <
		current.pointProxyPixelThreshold - 1.0e-6f ||
	    progressiveRank(candidate.progressiveCeiling) >
		progressiveRank(current.progressiveCeiling);
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
	    progressiveRank(candidate.progressiveCeiling) ==
		progressiveRank(current.progressiveCeiling);
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
    struct StableCapacityEvidence {
	bool proven = false;
	long double renderCostPerSecond = 0.0L;
    };

    /* The first quiet traversal of a newly published retained population may
     * pay command construction, software-code warmup, or immutable-buffer
     * installation which is not reflected by the assembly preparation
     * serial.  A hard abort is still a responsiveness boundary, but it is not
     * yet proof that the retained cuts themselves are unsustainable.  Permit
     * exactly one unchanged retry while a typed population/frame obligation
     * identifies that first presentation.  A second abort is ordinary
     * capacity evidence.  Interactive input never waits for this retry. */
    static bool retryTransientPresentation(bool interactive,
	unsigned int consecutiveInterruptedFrames,
	bool preparationChanged, bool publicationFramePending,
	bool refinementFramePending, bool pointCalibrationPending)
    {
	if (preparationChanged)
	    return true;
	return !interactive && consecutiveInterruptedFrames == 1 &&
	    (publicationFramePending || refinementFramePending ||
	     pointCalibrationPending);
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
	    static_cast<long double>(hardDeadlineNanoseconds) * 0.95L /
	    static_cast<long double>(renderNanoseconds);
	if (!std::isfinite(predicted) ||
	    predicted >= static_cast<long double>(SIZE_MAX))
	    return SIZE_MAX;
	return std::max(presentedCost, static_cast<size_t>(predicted));
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
     * policy allocation-free and directly testable.  No further tiers are
     * generated: 0.25 pixels is a raster-stable terminal target, not a route
     * to unconditionally loading full geometry.
     */
    static float stablePixelError(float currentPixelError,
	size_t activeCost, size_t sceneBudget,
	uint64_t renderNanoseconds, float targetFps,
	bool exactCompletedFrame, bool residentMemoryHeadroom)
    {
	if (!std::isfinite(currentPixelError) ||
	    currentPixelError <= 0.2501f || currentPixelError > 1.01f ||
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
 * waits for submission and submission waits for publication.
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
	bool servicePending, bool publicationPending)
    {
	return publicationPending || canProduceGeometry(
	    submissionPending, submissionPausedByPresentation,
	    providerPending, servicePending);
    }
};

/*
 * Coalesce append-only source-inventory deltas before starting another LoD
 * planning transaction.  The source provider has already installed useful
 * boxes in the retained scene; submitting its first published mesh contract
 * immediately establishes time-to-first-mesh, but running one completed LoD
 * pass for every later 16 ms provider drain serializes a parallel 50k source
 * behind hundreds of owner-thread passes.
 *
 * This policy delays only a pure inventory continuation.  A camera/policy or
 * source-identity change, an active submission cursor, and the producer's
 * final pending->idle edge all bypass it.  The bounded age means a producer
 * which never fills a large batch still makes progress.
 */
class BObolLodInventoryDeltaPolicy {
public:
    bool defer(bool inventoryChanged, bool providerPending,
	bool submissionPending, bool initialSubmissionComplete,
	bool interactive, int64_t nowMicroseconds)
    {
	if (!inventoryChanged || !providerPending || submissionPending ||
	    !initialSubmissionComplete) {
	    this->firstPendingMicrosecondsValue = 0;
	    return false;
	}
	const int64_t now = nowMicroseconds > 0 ? nowMicroseconds : 1;
	if (this->firstPendingMicrosecondsValue <= 0)
	    this->firstPendingMicrosecondsValue = now;
	const int64_t limit = interactive ? 100000 : 250000;
	const int64_t age = now >= this->firstPendingMicrosecondsValue ?
	    now - this->firstPendingMicrosecondsValue : limit;
	return age < limit;
    }

    void committed(void)
    {
	this->firstPendingMicrosecondsValue = 0;
    }

    void reset(void)
    {
	this->committed();
    }

    int64_t firstPendingMicroseconds(void) const
    {
	return this->firstPendingMicrosecondsValue;
    }

private:
    int64_t firstPendingMicrosecondsValue = 0;
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
class BObolLodPointProxyCalibrationPolicy {
public:
    struct Decision {
	float threshold = 1.0f;
	bool changed = false;
	bool continueRelaxation = false;
    };

    static bool applicable(size_t visibleOccurrenceCount)
    {
	/* This controller exists to bound the per-occurrence floor of a
	 * many-part scene.  One prominent mesh has no dispatch population to
	 * aggregate; increasing its point threshold cannot help until it destroys
	 * the very silhouette LoD is required to preserve. */
	return visibleOccurrenceCount > 1;
    }

    /* The one-pixel presentation contract deliberately permits genuinely
     * subpixel occurrences to remain in the assembly's point batch.  They
     * are the bounded replacement for thousands of independent minimum-mesh
     * draws, not evidence that a temporary coarse population cut remains.
     * Recovery is complete once the exact classifier frame contains no
     * structural boxes and no renderer-wide PoP ceiling. */
    static bool onePixelPresentationReady(bool haveCadAssemblies,
	bool exactFrame, bool exactOccurrenceClassification,
	size_t presentedSubpixelOccurrences,
	size_t presentedStructuralBoxes, int progressiveCeiling,
	float pointProxyPixelThreshold, size_t managedRenderCost)
    {
	(void)presentedSubpixelOccurrences;
	return haveCadAssemblies && exactFrame &&
	    exactOccurrenceClassification &&
	    presentedStructuralBoxes == 0 && progressiveCeiling < 0 &&
	    std::isfinite(pointProxyPixelThreshold) &&
	    pointProxyPixelThreshold <= 1.01f && managedRenderCost > 0;
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
	bool servicePending, bool publicationPending)
    {
	return BObolLodProducerPolicy::ownsFutureFrame(
	    submissionPending, submissionPausedByCalibration,
	    providerPending, servicePending, publicationPending);
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

    void reset(void)
    {
	this->unsafeThreshold = 0.0f;
	this->safeThreshold = 0.0f;
    }

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
    bool requiresReusableConfirmation(float currentThreshold) const
    {
	return sanitize(currentThreshold) > 1.01f;
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

private:
    static float sanitize(float threshold)
    {
	return std::isfinite(threshold) ?
	    std::max(1.0f, std::min(64.0f, threshold)) : 1.0f;
    }

    float unsafeThreshold = 0.0f;
    float safeThreshold = 0.0f;
};

/* Event-driven static quality is a separate presentation phase, not an
 * annotation which any ordinary quiet-frame failure may retire.  Keep its
 * two release-critical predicates independent of controller latch ordering so
 * unit tests (and the formal coordinator model) can enforce that boundary. */
class BObolLodStaticQualityPolicy {
public:
    static bool rejectAfterInterruptedFrame(bool interactive,
	bool cadDrawAttempted, bool preparationOnlyRetry,
	bool staticTrialActive)
    {
	return !interactive && cadDrawAttempted && !preparationOnlyRetry &&
	    staticTrialActive;
    }

    static bool onePixelTrialRequired(float pointProxyPixelThreshold)
    {
	return std::isfinite(pointProxyPixelThreshold) &&
	    pointProxyPixelThreshold > 1.01f;
    }

    /* Point aggregation is a convergence aid for a changing or streaming
     * many-part scene, not the final fidelity policy of an event-driven
     * framebuffer.  Once an exact, reusable point population has no producer
     * behind it, test the normal one-pixel occurrence contract directly under
     * the independently interruptible static deadline.  Otherwise the point
     * controller first compacts retained triangle prefixes to satisfy the
     * preferred multi-frame cadence, only for the static-quality pass to
     * restore those same prefixes immediately afterward.
     *
     * Keep this decision scalar and explicit: the controller owns state
     * transitions, while this predicate documents the proof required to skip
     * the destructive intermediate phase.  A rejected trial is terminal for
     * the current view/capacity epoch and must not be repeated. */
    static bool startOnePixelTrialFromSettledPointFrame(bool interactive,
	bool exactReusableFrame, bool producerWork, bool resourcePressure,
	bool staticTrialRejected, float pointProxyPixelThreshold,
	uint64_t renderNanoseconds, uint64_t preferredFrameNanoseconds,
	uint64_t staticDeadlineNanoseconds)
    {
	return !interactive && exactReusableFrame && !producerWork &&
	    !resourcePressure && !staticTrialRejected &&
	    onePixelTrialRequired(pointProxyPixelThreshold) &&
	    renderNanoseconds > 0 && preferredFrameNanoseconds > 0 &&
	    staticDeadlineNanoseconds > preferredFrameNanoseconds &&
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
class BObolLodProviderPacingPolicy {
public:
    static uint64_t effectiveMicroseconds(uint64_t configuredMicroseconds,
	bool controllerOwnedDefault, bool interactive)
    {
	if (!controllerOwnedDefault || configuredMicroseconds <= 4000)
	    return configuredMicroseconds;
	if (interactive)
	    return std::min<uint64_t>(configuredMicroseconds, 4000);
	return configuredMicroseconds;
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

/* A provider completion can append a richer immutable PoP suffix without
 * changing the active draw cut.  That is deliberately not a per-result
 * admission decision: only the scene-wide allocator knows which occurrences
 * should spend the calibrated render budget.  Coalesce all such completions
 * into one allocation edge after the worker/result wave becomes idle.
 *
 * Keeping this obligation separate from result publication is important.
 * Publication controls when newly bound arrays become visible; resident
 * growth controls when the complete occurrence population is rebalanced.
 * Conflating the two either makes result order the quality policy or runs an
 * O(scene) allocator once per result batch. */
class BObolLodResidentGrowthPolicy {
public:
    /* A quiet residency-drain result may replace the immutable backing
     * generation without changing the pixels selected by the owner-thread
     * allocation.  Such a result is not a publication edge: keep the last
     * coherent framebuffer while the rest of the resident wave drains, then
     * let the scene-wide allocator publish one complete population.
     *
     * Keep every premise explicit.  In particular, ordinary interaction,
     * first-useful coverage, asset replacement, and a changed active cut must
     * remain immediately presentable rather than being hidden behind a
     * background residency transaction. */
    static bool canRetainPresentation(bool residencyDrainActive,
	bool retainedPayload, bool sameAsset, bool activeCutPreserved,
	bool richerResidentPrefix)
    {
	return residencyDrainActive && retainedPayload && sameAsset &&
	    activeCutPreserved && richerResidentPrefix;
    }

    void noteRicherPrefixAvailable(void)
    {
	this->pendingValue = true;
	this->residencyDrainRequiredValue = true;
    }

    /* A completed resident batch is not necessarily the terminal resident
     * population.  Before spending O(scene) work on importance allocation,
     * let one ordinary source pass refill the bounded worker queue without
     * changing presentation cuts.  Results arriving during that pass re-arm
     * this edge, so an arbitrary number of cache batches collapse into one
     * terminal allocation. */
    bool beginResidencyDrainIfReady(bool automatic, bool streamIdle,
	bool workAllowed)
    {
	if (!this->pendingValue || !this->residencyDrainRequiredValue ||
	    !automatic || !streamIdle || !workAllowed)
	    return false;
	this->residencyDrainRequiredValue = false;
	return true;
    }

    bool consumeIfReady(bool automatic, bool streamIdle,
	bool allocationAllowed)
    {
	if (!this->pendingValue || !automatic || !streamIdle ||
	    !allocationAllowed || this->residencyDrainRequiredValue)
	    return false;
	this->pendingValue = false;
	return true;
    }

    bool pending(void) const
    {
	return this->pendingValue;
    }

    bool residencyDrainRequired(void) const
    {
	return this->residencyDrainRequiredValue;
    }

    void reset(void)
    {
	this->pendingValue = false;
	this->residencyDrainRequiredValue = false;
    }

private:
    bool pendingValue = false;
    bool residencyDrainRequiredValue = false;
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

/*
 * Typed presentation acknowledgement owned by the LoD coordinator.
 *
 * A frame barrier is a transaction, not a pair of loosely related Boolean and
 * serial fields.  Every obligation records the view/policy epoch for which it
 * was created, the first render completion which may satisfy it, and all
 * reasons sharing that frame.  Adding another reason never moves an existing
 * barrier forward; an unrelated HUD/selection frame may satisfy it only when
 * it is a complete CAD presentation in the same epochs.
 */
class BObolLodFrameObligation {
public:
    enum Reason : uint32_t {
	REASON_NONE = 0,
	REASON_CUT_PRESENTATION = 1u << 0,
	REASON_RESULT_PUBLICATION = 1u << 1,
	REASON_RESIDENT_REFINEMENT = 1u << 2,
	REASON_CALIBRATION = 1u << 3,
	REASON_HANDOFF = 1u << 4
    };

    struct Completion {
	bool retired = false;
	bool stale = false;
	uint32_t reasons = REASON_NONE;
    };

    bool arm(Reason reason, uint64_t completedRenderSerial,
	uint64_t viewEpoch, uint64_t policyEpoch)
    {
	if (reason == REASON_NONE)
	    return false;
	const uint64_t required = completedRenderSerial == UINT64_MAX ?
	    UINT64_MAX : completedRenderSerial + 1;
	const bool sameTransaction = this->pendingValue &&
	    this->viewEpochValue == viewEpoch &&
	    this->policyEpochValue == policyEpoch;
	if (sameTransaction) {
	    this->reasonMask |= static_cast<uint32_t>(reason);
	    return false;
	}
	this->pendingValue = true;
	this->reasonMask = static_cast<uint32_t>(reason);
	this->requiredRenderSerialValue = required;
	this->viewEpochValue = viewEpoch;
	this->policyEpochValue = policyEpoch;
	if (++this->sequenceValue == 0)
	    ++this->sequenceValue;
	return true;
    }

    Completion complete(uint64_t completedRenderSerial,
	uint64_t viewEpoch, uint64_t policyEpoch)
    {
	Completion result;
	if (!this->pendingValue)
	    return result;
	if (this->viewEpochValue != viewEpoch ||
	    this->policyEpochValue != policyEpoch) {
	    result.stale = true;
	    this->reset();
	    return result;
	}
	if (completedRenderSerial < this->requiredRenderSerialValue)
	    return result;
	result.retired = true;
	result.reasons = this->reasonMask;
	this->reset();
	return result;
    }

    void reset(void)
    {
	this->pendingValue = false;
	this->reasonMask = REASON_NONE;
	this->requiredRenderSerialValue = 0;
	this->viewEpochValue = 0;
	this->policyEpochValue = 0;
    }

    bool pending(void) const
    {
	return this->pendingValue;
    }

    uint32_t reasons(void) const
    {
	return this->reasonMask;
    }

    uint64_t requiredRenderSerial(void) const
    {
	return this->pendingValue ? this->requiredRenderSerialValue : 0;
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

private:
    bool pendingValue = false;
    uint32_t reasonMask = REASON_NONE;
    uint64_t requiredRenderSerialValue = 0;
    uint64_t viewEpochValue = 0;
    uint64_t policyEpochValue = 0;
    uint64_t sequenceValue = 0;
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
