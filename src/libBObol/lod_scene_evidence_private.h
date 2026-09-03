/*       L O D _ S C E N E _ E V I D E N C E _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_SCENE_EVIDENCE_PRIVATE_H
#define LIBBOBOL_LOD_SCENE_EVIDENCE_PRIVATE_H

#include "common.h"

#include "identity_counter_private.h"
#include "lod_revision_private.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <unordered_map>
#include <vector>

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
	bobol_identity_advance(this->pressureRevision);
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
	bool demandDeferred = false;
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
	    this->demandCensusRequiredValue = false;
	    this->demandDeferredValue = false;
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

    bool shouldResumeSubmission(bool interactive, bool submissionActive,
	bool providerPending) const
    {
	return this->effectiveActive() && !interactive &&
	    !submissionActive && !providerPending;
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

    /* A structural-only coverage pass proves where visible occurrences are,
     * but deliberately does not create their mesh demand.  Keep that fact
     * separate from the visibility proof: a sparse unsatisfied frontier is
     * authoritative only after one later dense demand census has visited the
     * complete visible population. */
    void noteDemandDeferred(void)
    {
	if (this->activeValue || this->demandCensusRequiredValue)
	    this->demandDeferredValue = true;
    }

    bool demandCensusRequired(void) const
    {
	return this->demandCensusRequiredValue;
    }

    /* Return true only when this completed pass actually established demand.
     * A point-classification owner may defer mesh admission again; retain the
     * requirement in that case so a later ordinary pass cannot use a sparse
     * frontier whose missing-payload entries were never recorded. */
    bool completeDemandCensus(void)
    {
	if (!this->demandCensusRequiredValue)
	    return false;
	if (this->demandDeferredValue) {
	    this->demandDeferredValue = false;
	    return false;
	}
	this->demandCensusRequiredValue = false;
	return true;
    }

    void clearPassCounters(void)
    {
	this->sawBoundedSourceValue = false;
	this->demandDeferredValue = false;
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
	completion.demandDeferred = this->demandDeferredValue;
	completion.visibleCount = this->visibleCountValue;
	completion.coveredCount = this->coveredCountValue;
	completion.missing = completion.coveredCount < completion.visibleCount;
	this->coverageCompleteValue = !completion.missing;
	this->demandCensusRequiredValue =
	    !completion.missing && completion.demandDeferred;
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
	this->demandCensusRequiredValue = false;
	this->demandDeferredValue = false;
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

    /* Coverage and demand are distinct proofs.  The first flag survives the
     * completed structural pass; the second records whether an attempted
     * successor was itself deferred by a stronger presentation owner. */
    bool demandCensusRequiredValue = false;
    bool demandDeferredValue = false;
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
	size_t terminalCost, size_t currentCost, size_t remainingCost)
    {
	if (!safeScene || !visible || subpixelPresentation || !botSource ||
	    !supportedDrawMode || terminalCost == 0)
	    return false;
	const size_t additionalCost = terminalCost > currentCost ?
	    terminalCost - currentCost : 0;
	return additionalCost <= remainingCost;
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
	    bobol_identity_advance(bound->generation);
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
	bool enabled = true;
	size_t expectedLeafCount = 0;
	size_t availableLeafCount = 0;
	size_t visibleTargetCount = 0;
	size_t activePayloadCount = 0;
	size_t satisfiedPayloadCount = 0;
	size_t presentedSubpixelOccurrenceCount = 0;
	size_t presentedStructuralBoxCount = 0;
	size_t terminalProxyOccurrenceCount = 0;
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
	uint64_t sourcePreparationCompletedUnits = 0;
	uint64_t sourcePreparationTotalUnits = 0;
	bool submissionPending = false;
	bool resultPending = false;
	bool publicationPending = false;
	bool calibrationPending = false;
	uint64_t capacitySearchCompletedUnits = 0;
	uint64_t capacitySearchTotalUnits = 0;
	bool controlPending = false;
	bool interactive = false;
	bool compactionPending = false;
	bool progressiveWorkPending = false;
	bool gpuMemoryPressure = false;
	bool stableBudgetLimited = false;
	/* A current allocation which selected the complete pixel-demand
	 * population has been presented exactly.  This is stronger than the
	 * occurrence-local satisfied counter, which cannot express fractional PoP
	 * cuts between two integer hierarchy ordinals. */
	bool pixelDemandPresentationProven = false;
	/* A completed capacity allocation or static hard-deadline trial may
	 * deliberately retain cuts below pixel demand.  This is terminal quality
	 * only when exact current-revision evidence proves the represented
	 * population was presented. */
	bool constrainedPresentationProven = false;
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
	/* Full-detail rendering does not publish into the LoD occurrence ledger.
	 * Disabling automatic LoD therefore normally terminates this coordinator
	 * rather than asking absent LoD counters to prove the ordinary draw path.
	 * A caller may still submit an explicit progressive generation, however,
	 * and its immutable result is not presented until the publication barrier
	 * observes a completed frame.  Keep that concrete transaction and any
	 * downstream reducer obligation visible even though automatic scheduling
	 * remains disabled. */
	if (!inputs.enabled) {
	    if (inputs.publicationPending || inputs.controlPending) {
		decision.phase = Phase::REFINING;
		decision.outcome = Outcome::ACTIVE;
		decision.fraction = 0.0f;
		decision.terminal = false;
		decision.viewReady = false;
		decision.visualPending = true;
		decision.hasLodState = true;
		return decision;
	    }
	    decision.viewReady = true;
	    return decision;
	}
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
	const bool nonResultForegroundPending =
	    structuralPending || inputs.structuralDiscovery ||
	    inputs.sourcePreparationPending ||
	    inputs.submissionPending ||
	    inputs.publicationPending ||
	    inputs.calibrationPending || inputs.controlPending ||
	    unresolvedStructuralBoxes > 0;
	const size_t representedPayloads = saturatingAdd(
	    saturatingAdd(inputs.activePayloadCount,
		inputs.presentedSubpixelOccurrenceCount),
	    saturatingAdd(inputs.terminalProxyOccurrenceCount,
		inputs.terminalOccurrenceFailureCount));
	const size_t satisfiedPayloads = saturatingAdd(
	    saturatingAdd(inputs.satisfiedPayloadCount,
		inputs.presentedSubpixelOccurrenceCount),
	    saturatingAdd(inputs.memoryLimitedPayloadCount,
		saturatingAdd(inputs.terminalProxyOccurrenceCount,
		    inputs.terminalOccurrenceFailureCount)));
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
	const bool allocationQualityProven =
	    inputs.pixelDemandPresentationProven ||
	    inputs.constrainedPresentationProven;
	const size_t terminalSatisfiedPayloads = allocationQualityProven ?
	    representedPayloads : satisfiedPayloads;
	const bool populatedViewResolved = inputs.visibleTargetCount > 0 &&
	    representedPayloads >= inputs.visibleTargetCount &&
	    terminalSatisfiedPayloads >= inputs.visibleTargetCount;
	const bool noVisiblePopulation = inputs.expectedLeafCount == 0 &&
	    inputs.availableLeafCount == 0 &&
	    inputs.visibleTargetCount == 0 &&
	    inputs.presentedStructuralBoxCount == 0;
	const bool failedSourceResolved = inputs.failedSourceCount > 0;
	const bool viewPresentationResolved =
	    !nonResultForegroundPending && !inputs.interactive &&
	    (noVisiblePopulation || emptyViewResolved || populatedViewResolved ||
	     failedSourceResolved);
	/* A queued result is foreground only while the current view still needs
	 * it for coverage or quality.  Once the current presentation is complete,
	 * an in-flight reusable-cache producer and its eventual result are one
	 * background transaction.  Reclassifying the result as foreground merely
	 * because it crossed the worker/owner-thread boundary reopens a terminal
	 * view without a semantic revision.  Applying a result which changes the
	 * represented scene advances the appropriate revision and starts the next
	 * foreground epoch normally. */
	const bool backgroundTaskWork = viewPresentationResolved &&
	    (inputs.pendingTasks > 0 || inputs.inFlight > 0 ||
	     inputs.resultPending);
	decision.visualPending = nonResultForegroundPending ||
	    !viewPresentationResolved ||
	    (!backgroundTaskWork &&
	     (inputs.pendingTasks > 0 || inputs.inFlight > 0 ||
	      inputs.resultPending));
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
	     decision.memoryLimited);
	decision.outcome = hasTerminalError ? Outcome::ERROR :
	    !decision.terminal ? Outcome::ACTIVE :
	    decision.performanceLimited ? Outcome::CONSTRAINED : Outcome::READY;
	decision.hasLodState =
	    inputs.expectedLeafCount > 0 || inputs.availableLeafCount > 0 ||
	    inputs.visibleTargetCount > 0 || inputs.activePayloadCount > 0 ||
	    inputs.terminalProxyOccurrenceCount > 0 ||
	    inputs.gpuTrackedBufferBytes > 0 || decision.visualPending ||
	    decision.backgroundPending;
	if (!decision.hasLodState)
	    return decision;

	/* An unfinished leaf census owns discovery.  Once that census is complete,
	 * the source producer is the more useful user-facing owner even if an
	 * autoview revision has not yet rebuilt its visibility census.  Reporting
	 * the derived structuralDiscovery alias first made a cold draw appear to
	 * restart at zero while its exact source-preparation rank remained intact. */
	if (structuralPending ||
	    (inputs.structuralDiscovery && !inputs.sourcePreparationPending)) {
	    decision.phase = Phase::DISCOVERING;
	    const long double coverage =
		structuralPending && inputs.expectedLeafCount > 0 ?
		    static_cast<long double>(inputs.availableLeafCount) /
			static_cast<long double>(inputs.expectedLeafCount) : 0.0L;
	    decision.fraction = static_cast<float>(
		std::min<long double>(0.40L, 0.40L * coverage));
	} else if (inputs.sourcePreparationPending) {
	    decision.phase = Phase::PREPARING;
	    if (inputs.sourcePreparationTotalUnits > 0) {
		const uint64_t completed = std::min(
		    inputs.sourcePreparationCompletedUnits,
		    inputs.sourcePreparationTotalUnits);
		const long double coverage =
		    static_cast<long double>(completed) /
		    static_cast<long double>(
			inputs.sourcePreparationTotalUnits);
		/* Structural coverage owns the first 40 percent.  Exact finite
		 * source preparation owns the next 35 percent; coherent view
		 * allocation and presentation retain the final quarter. */
		decision.fraction = static_cast<float>(
		    0.40L + std::min<long double>(0.35L, 0.35L * coverage));
	    } else {
		/* The visibility census may itself still be streaming while source
		 * preparation is active.  Using that partial census as the
		 * denominator made the fraction jump to 75 percent as soon as the
		 * first batch was represented, regardless of how much of a large
		 * model remained.  Fall back to the known structural population
		 * only until the producer publishes an exact work rank. */
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
		    decision.fraction = static_cast<float>(
			0.40L + std::min<long double>(
			    0.35L, 0.35L * coverage));
		}
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
	static constexpr long double discoveryWeight = 0.40L;
	static constexpr long double representationWeight = 0.35L;
	static constexpr long double richQualityWeight = 0.19L;
	static constexpr long double capacitySearchWeight = 0.05L;
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
	    const size_t terminalSatisfiedPayloads =
		inputs.constrainedPresentationProven ?
		    inputs.activePayloadCount : inputs.satisfiedPayloadCount;
	    richCoverage =
		static_cast<long double>(std::min(
		    terminalSatisfiedPayloads, richTarget)) /
		static_cast<long double>(richTarget);
	} else if (inputs.visibilityCensusComplete && represented == target) {
	    /* An exact all-subpixel view has no rich tail to resolve. */
	    richCoverage = 1.0L;
	}
	/* Structural discovery owns the first 40 percent and initial useful
	 * representations the next 35.  The visible-mesh tail is weighted
	 * separately: classifying one terminal subpixel proxy is much cheaper than
	 * loading, publishing, and validating one mesh, so treating both as equal
	 * objects made a 150k-part scene jump directly to 99 percent with hundreds
	 * of expensive mesh replacements still outstanding.  Gate rich progress
	 * by representation coverage so a partial visibility census cannot claim a
	 * nearly complete tail.  A bounded capacity proof owns the next five
	 * percent, using its exact allocation/presentation/sample rank instead of
	 * leaving the bar at 99 percent throughout the search.  The last percent
	 * remains the coherent-frame handoff obligation; only a proven ready view
	 * reports 100 percent. */
	long double capacitySearchCoverage = 0.0L;
	if (inputs.capacitySearchTotalUnits > 0) {
	    capacitySearchCoverage =
		static_cast<long double>(std::min(
		    inputs.capacitySearchCompletedUnits,
		    inputs.capacitySearchTotalUnits)) /
		static_cast<long double>(inputs.capacitySearchTotalUnits);
	}
	return static_cast<float>(
	    discoveryWeight +
	    std::min(representationWeight,
		representationWeight * representationCoverage) +
	    std::min(richQualityWeight,
		richQualityWeight * representationCoverage * richCoverage) +
	    std::min(capacitySearchWeight,
		capacitySearchWeight * capacitySearchCoverage));
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

#endif /* LIBBOBOL_LOD_SCENE_EVIDENCE_PRIVATE_H */
