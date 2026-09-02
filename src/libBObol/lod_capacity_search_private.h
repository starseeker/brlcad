/*       L O D _ C A P A C I T Y _ S E A R C H _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_CAPACITY_SEARCH_PRIVATE_H
#define LIBBOBOL_LOD_CAPACITY_SEARCH_PRIVATE_H

#include "common.h"
#include "lod_revision_private.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>

/*
 * Capacity evidence is the output of a search, so the search cannot use the
 * pipeline's changing capacity-output epoch as one of its own inputs.  The
 * Inventory, availability, view, and policy, preferred and maximum target
 * cadences, and the allocator's maximum realizable budget identify the frozen
 * problem.  Exact visibility edits deliberately do not: they require a new
 * occurrence allocation, but a terminal renderer-cost budget remains a useful
 * conservative certificate for that same immutable scene and endpoint.
 * Publishing a successor or terminal result advances the ordinary capacity
 * epoch outside this value.
 */
struct BObolLodCapacitySearchKey {
    BObolLodInventoryEpoch inventory;
    BObolLodAvailabilityEpoch availability;
    BObolLodViewEpoch view;
    BObolLodPolicyEpoch policy;
    uint64_t preferredTargetNanoseconds = 0;
    uint64_t maximumTargetNanoseconds = 0;
    /* A ready view change retains a previously completed occurrence
     * population.  After pose or zoom, its first capacity question is whether
     * that useful presentation fits the event-driven static deadline, not
     * whether it meets the preferred redraw cadence.  Zoom still recomputes
     * pixel demand; this flag governs only the capacity-search lower bound.
     * Keep the choice explicit: equal target values hide STATIC ownership from
     * the renderer and let a shorter endpoint abort corrupt the bracket. */
    bool startAtStatic = false;
    /* Point aggregation changes both occurrence membership and renderer
     * cost.  A certificate measured at another threshold belongs to a
     * different search problem even when the scene epochs are unchanged. */
    float pointProxyPixelThreshold = 1.0f;
    /* Largest complete retained allocation which the current allocator
     * policy can realize.  The view's pixel demand may be larger and is not
     * part of this finite search domain. */
    size_t maximumCandidateBudget = 0;
    /* One complete occurrence-local minimum is the lowest candidate the
     * allocator can represent.  A search below this floor cannot produce a
     * different framebuffer and must advance to its next deadline goal. */
    size_t minimumBudget = 0;
    /* A miss at the preferred deadline narrows only the steady search.  It
     * must not erase the richer candidate domain available to the longer,
     * event-driven static deadline.  A miss at that longer deadline supplies
     * the independent maximum-goal ceiling. */
    size_t preferredBudgetCeiling = SIZE_MAX;
    size_t maximumBudgetCeiling = SIZE_MAX;

    bool sameProblem(const BObolLodCapacitySearchKey &other) const
    {
	return this->inventory == other.inventory &&
	    this->availability == other.availability &&
	    this->view == other.view && this->policy == other.policy &&
	    this->preferredTargetNanoseconds ==
		other.preferredTargetNanoseconds &&
	    this->maximumTargetNanoseconds ==
		other.maximumTargetNanoseconds &&
	    this->startAtStatic == other.startAtStatic &&
	    std::memcmp(&this->pointProxyPixelThreshold,
		&other.pointProxyPixelThreshold, sizeof(float)) == 0 &&
	    this->maximumCandidateBudget == other.maximumCandidateBudget &&
	    this->minimumBudget == other.minimumBudget;
    }

    bool same(const BObolLodCapacitySearchKey &other) const
    {
	return this->sameProblem(other) &&
	    this->preferredBudgetCeiling ==
		other.preferredBudgetCeiling &&
	    this->maximumBudgetCeiling == other.maximumBudgetCeiling;
    }

    bool coversRevisions(const BObolLodAdmissionRevisionStamp &stamp) const
    {
	/* Capacity is an output of this search and exact visibility is an
	 * allocation input.  Both are deliberately absent from this frozen
	 * renderer-capacity tuple.  Their revisions may advance without changing
	 * the measured budget this candidate is meant to bound. */
	return this->inventory == stamp.inventory &&
	    this->availability == stamp.availability &&
	    this->view == stamp.view && this->policy == stamp.policy;
    }

    bool valid(void) const
    {
	return this->preferredTargetNanoseconds > 0 &&
	    this->maximumTargetNanoseconds >=
		this->preferredTargetNanoseconds &&
	    std::isfinite(this->pointProxyPixelThreshold) &&
	    this->pointProxyPixelThreshold >= 1.0f &&
	    this->maximumCandidateBudget > 0 &&
	    this->minimumBudget <= this->maximumCandidateBudget;
    }
};

/*
 * Bounded search over complete, visual-importance-ordered scene allocations.
 * A candidate is an allocation budget, not a triangle or object ordinal.  It
 * consumes exactly sampleLimit valid unchanged frames, then strictly moves a
 * known-safe lower bound or a known-unsafe upper bound.  At most
 * candidateLimit distinct budgets may be measured for each of the key's two
 * ordered goals.  Once the final goal has both a safe population and a
 * rejected richer population, it may measure one midpoint bridge before
 * settling.  That bounded bridge recovers large amounts of otherwise unused
 * static capacity without exposing an open-ended binary search or a visible
 * quality staircase.  The optional static goal inherits the steady-safe
 * lower bound and cannot transition back.
 *
 * Throughput estimates are deliberately confined to proposing the next
 * candidate.  Only completed-frame classifications change the bracket.  This
 * is the executable counterpart of ObolCapacitySearch.tla.
 */
class BObolLodCapacitySearchCertificate {
private:
    static constexpr unsigned int SampleLimit = 3;
    static constexpr unsigned int CandidateLimit = 4;
    static constexpr unsigned int InvalidSampleLimit = 2;

public:
    enum class Goal : uint8_t {
	STEADY = 0,
	STATIC
    };

    enum class Phase : uint8_t {
	INACTIVE = 0,
	ALLOCATING,
	PRESENTING,
	MEASURING,
	TERMINAL
    };

    enum class Result : uint8_t {
	NONE = 0,
	REQUEST_SAMPLE,
	REALLOCATE,
	CERTIFIED,
	UNMEASURABLE,
	STALE_POPULATION,
	PRESENTATION_REQUIRED
    };

    struct Observation {
	BObolLodCapacitySearchKey key;
	size_t candidateBudget = 0;
	size_t presentedCost = 0;
	/* Exact retained-allocation identity for this semantic epoch.  Different
	 * numeric budgets may select the same discrete set of PoP cuts. */
	uint64_t populationSignature = 0;
	/* Exact numeric interval represented by this discrete allocator
	 * population.  The minimum is its selected render cost.  When the
	 * successor flag is set, zero means this is the richest resident
	 * population and any positive value is the first budget which can change
	 * the selected cuts. */
	size_t populationMinimumBudget = 0;
	size_t nextDistinctPopulationBudget = 0;
	bool nextDistinctPopulationBudgetKnown = false;
	size_t knownSafeBudget = 0;
	uint64_t observedNanoseconds = 0;
	bool validSample = false;
    };

    struct Decision {
	Result result = Result::NONE;
	size_t budget = 0;
	size_t safeBudget = 0;
	size_t unsafeBudget = SIZE_MAX;
	unsigned int samplesRemaining = 0;
	unsigned int measuredCandidates = 0;
	unsigned int totalMeasuredCandidates = 0;

	bool requestsFrame(void) const
	{
	    return this->result == Result::REQUEST_SAMPLE;
	}

	bool requestsReallocation(void) const
	{
	    return this->result == Result::REALLOCATE;
	}

	bool terminal(void) const
	{
	    return this->result == Result::CERTIFIED ||
		this->result == Result::UNMEASURABLE ||
		this->result == Result::STALE_POPULATION;
	}
    };

    static constexpr unsigned int sampleLimit(void) { return SampleLimit; }
    static constexpr unsigned int candidateLimit(void)
    {
	return CandidateLimit;
    }
    static constexpr unsigned int totalCandidateLimit(void)
    {
	return 2 * CandidateLimit;
    }
    static constexpr unsigned int invalidSampleLimit(void)
    {
	return InvalidSampleLimit;
    }

    static BObolLodCapacitySearchKey keyFor(
	const BObolLodAdmissionRevisionStamp &stamp,
	uint64_t preferredTargetNanoseconds,
	uint64_t maximumTargetNanoseconds, float pointProxyPixelThreshold,
	size_t maximumCandidateBudget,
	size_t minimumBudget = 0, bool startAtStatic = false)
    {
	BObolLodCapacitySearchKey key;
	key.inventory = stamp.inventory;
	key.availability = stamp.availability;
	key.view = stamp.view;
	key.policy = stamp.policy;
	key.preferredTargetNanoseconds = preferredTargetNanoseconds;
	key.maximumTargetNanoseconds = maximumTargetNanoseconds;
	key.startAtStatic = startAtStatic;
	key.pointProxyPixelThreshold = pointProxyPixelThreshold;
	key.maximumCandidateBudget = maximumCandidateBudget;
	key.minimumBudget = minimumBudget;
	return key;
    }

    void reset(void)
    {
	*this = BObolLodCapacitySearchCertificate();
    }

    Decision prepare(const Observation &observation)
    {
	if (!observation.key.valid() || !observation.candidateBudget)
	    return Decision();

	bool began = false;
	if (this->phaseValue == Phase::INACTIVE ||
	    !this->keyValue.same(observation.key)) {
	    this->begin(observation);
	    began = true;
	}
	if (this->phaseValue == Phase::TERMINAL) {
	    return this->decision(this->terminalResultValue,
		this->safeBudgetValue);
	}
	if (!this->bindCandidatePopulationInterval(observation))
	    return this->finish(Result::STALE_POPULATION);
	if (began && this->phaseValue == Phase::ALLOCATING)
	    return this->decision(Result::REALLOCATE,
		this->candidateBudgetValue);
	if (began && this->phaseValue == Phase::PRESENTING)
	    return this->decision(Result::REQUEST_SAMPLE,
		this->candidateBudgetValue);

	const size_t candidate = std::max(this->goalMinimumBudget(),
	    std::min(observation.candidateBudget,
		this->goalMaximumBudget()));
	if (candidate != this->candidateBudgetValue) {
	    if (this->measured(candidate))
		return this->finish(Result::STALE_POPULATION);
	    this->startCandidate(candidate, observation.presentedCost,
		observation.validSample ? Phase::MEASURING :
		    Phase::PRESENTING);
	    if (!this->bindCandidatePopulationInterval(observation))
		return this->finish(Result::STALE_POPULATION);
	} else if (this->phaseValue == Phase::ALLOCATING) {
	    /* A reallocation decision names the next complete population, but it
	     * does not own a framebuffer sample.  Admission describes the selected
	     * retained-allocation cost, while the renderer reports the canonical
	     * submitted-work cost; those are deliberately different currencies on
	     * paths such as OSMesa's expanded position stream.  Bind the population
	     * identity from the first completed frame below, not from this
	     * allocation barrier.  Otherwise a discrete candidate which changes no
	     * drawable cut is rejected as stale and reopens the same search forever. */
	    /* Applying the occurrence cuts is not sufficient when a temporary
	     * renderer ceiling still hides them.  Keep allocation ownership until
	     * an exact completed presentation names the candidate.  Otherwise the
	     * handoff frames consume the invalid-sample allowance and terminate a
	     * perfectly measurable search before its population reaches screen. */
	    this->phaseValue = observation.validSample ? Phase::MEASURING :
		Phase::PRESENTING;
	    this->candidatePresentedCostValue = 0;
	    this->candidatePopulationSignatureValue = 0;
	    if (this->phaseValue == Phase::PRESENTING)
		return this->decision(Result::REQUEST_SAMPLE,
		    this->candidateBudgetValue);
	} else if (this->phaseValue == Phase::PRESENTING) {
	    if (!observation.validSample)
		return this->decision(Result::REQUEST_SAMPLE,
		    this->candidateBudgetValue);
	    this->phaseValue = Phase::MEASURING;
	}
	if (this->candidatePresentedCostValue && observation.presentedCost &&
	    this->candidatePresentedCostValue != observation.presentedCost)
	    return this->finish(Result::STALE_POPULATION);
	if (!this->candidatePresentedCostValue)
	    this->candidatePresentedCostValue = observation.presentedCost;
	if (this->candidatePopulationSignatureValue &&
	    observation.populationSignature &&
	    this->candidatePopulationSignatureValue !=
		observation.populationSignature)
	    return this->finish(Result::STALE_POPULATION);
	if (!this->candidatePopulationSignatureValue)
	    this->candidatePopulationSignatureValue =
		observation.populationSignature;
	const int equivalent = this->measuredPopulation(
	    this->candidatePopulationSignatureValue,
	    this->candidatePresentedCostValue);
	if (equivalent >= 0)
	    return this->classifyEquivalentPopulation(
		static_cast<unsigned int>(equivalent));
	return this->decision(Result::REQUEST_SAMPLE,
	    this->candidateBudgetValue);
    }

    Decision observe(const Observation &observation)
    {
	const Decision prepared = this->prepare(observation);
	if (prepared.result != Result::REQUEST_SAMPLE)
	    return prepared;
	/* ALLOCATING requests the first exact presentation of the selected
	 * population.  It is a producer barrier, not a failed timing sample. */
	if (this->phaseValue == Phase::ALLOCATING)
	    return prepared;
	/* prepare() is also used by planning passes which cannot present a frame,
	 * so it deliberately leaves PRESENTING retries uncounted.  observe(), on
	 * the other hand, names a completed renderer attempt.  Bound completed
	 * attempts which still cannot supply exact CAD work: interrupted software
	 * frames can otherwise leave an unchanged view repainting forever. */
	if (this->phaseValue == Phase::PRESENTING) {
	    if (++this->invalidSampleCountValue <= invalidSampleLimit())
		return prepared;
	    return this->finish(Result::UNMEASURABLE);
	}

	if (!observation.validSample || !observation.observedNanoseconds ||
	    !observation.presentedCost) {
	    if (++this->invalidSampleCountValue <= invalidSampleLimit())
		return this->decision(Result::REQUEST_SAMPLE,
		    this->candidateBudgetValue);
	    return this->finish(Result::UNMEASURABLE);
	}

	this->invalidSampleCountValue = 0;
	const bool metTarget = observation.observedNanoseconds <=
	    this->activeTargetNanoseconds();
	if (metTarget)
	    this->safeSampleCountValue++;
	if (observation.observedNanoseconds <=
		observation.key.maximumTargetNanoseconds)
	    this->maximumSafeSampleCountValue++;
	this->proposalValues[this->sampleCountValue] = proposedBudget(
	    observation.presentedCost, this->activeTargetNanoseconds(),
	    observation.observedNanoseconds,
	    observation.key.maximumCandidateBudget, this->candidateBudgetValue,
	    metTarget);
	this->sampleCountValue++;
	if (this->sampleCountValue < sampleLimit())
	    return this->decision(Result::REQUEST_SAMPLE,
		this->candidateBudgetValue);

	return this->classifyCandidate();
    }

    /* A hard renderer deadline is an unsafe operational observation of the
     * active candidate transition, even when its exact framebuffer did not
     * finish.  Fold that conservative bound into the current bracket.
     * Resetting the certificate here would turn every abort into a new search
     * and permit an unchanged scene to alternate between the same allocations
     * forever. */
    Decision observeDeadlineMiss(const BObolLodCapacitySearchKey &updatedKey)
    {
	if (!this->awaitingSample() ||
	    !this->candidateBudgetValue)
	    return Decision();
	if (!this->keyValue.sameProblem(updatedKey))
	    return this->finish(Result::STALE_POPULATION);

	/* Deadline ceilings are monotone evidence within one semantic problem.
	 * Keep the certificate key synchronized with the externally retained
	 * evidence so the successor allocation cannot be mistaken for a rekey. */
	this->keyValue.preferredBudgetCeiling = std::min(
	    this->keyValue.preferredBudgetCeiling,
	    updatedKey.preferredBudgetCeiling);
	this->keyValue.maximumBudgetCeiling = std::min(
	    this->keyValue.maximumBudgetCeiling,
	    updatedKey.maximumBudgetCeiling);
	const size_t activeCeiling = this->goalValue == Goal::STATIC ?
	    this->keyValue.maximumBudgetCeiling :
	    this->keyValue.preferredBudgetCeiling;
	this->unsafeBudgetValue = std::min(this->unsafeBudgetValue,
	    strictUpperBound(this->keyValue.maximumCandidateBudget,
		activeCeiling));

	/* No timing extrapolation is available from an aborted frame.  Bisect the
	 * remaining proven bracket; candidateLimit still bounds all subsequent
	 * observations. */
	const size_t maximum = this->goalMaximumBudget();
	const size_t proposal = maximum > this->safeBudgetValue ?
	    this->safeBudgetValue +
		(maximum - this->safeBudgetValue) / 2 :
	    this->safeBudgetValue;
	return this->classifyCandidate(false, proposal);
    }

    Phase phase(void) const { return this->phaseValue; }
    Goal goal(void) const { return this->goalValue; }
    Result terminalResult(void) const { return this->terminalResultValue; }
    uint64_t activeTargetNanoseconds(void) const
    {
	return this->goalValue == Goal::STATIC ?
	    this->keyValue.maximumTargetNanoseconds :
	    this->keyValue.preferredTargetNanoseconds;
    }
    const BObolLodCapacitySearchKey &key(void) const
    {
	return this->keyValue;
    }
    size_t safeBudget(void) const { return this->safeBudgetValue; }
    size_t unsafeBudget(void) const { return this->unsafeBudgetValue; }
    size_t candidateBudget(void) const
    {
	return this->candidateBudgetValue;
    }
    unsigned int measuredCandidateCount(void) const
    {
	return this->measuredCandidateCountValue;
    }
    unsigned int totalMeasuredCandidateCount(void) const
    {
	return this->totalMeasuredCandidateCountValue;
    }
    unsigned int samplesRemaining(void) const
    {
	return this->phaseValue == Phase::MEASURING &&
	    this->sampleCountValue < sampleLimit() ?
	    sampleLimit() - this->sampleCountValue : 0;
    }

    unsigned int maximumCandidateCount(void) const
    {
	if (this->phaseValue == Phase::INACTIVE)
	    return 0;
	const bool hasSteadyAndStaticGoals = !this->keyValue.startAtStatic &&
	    this->keyValue.maximumTargetNanoseconds >
		this->keyValue.preferredTargetNanoseconds;
	return hasSteadyAndStaticGoals ? totalCandidateLimit() :
	    candidateLimit();
    }

    uint64_t progressTotalUnits(void) const
    {
	/* Each candidate has one allocation edge, one exact-presentation edge,
	 * and a fixed number of timing samples.  This is an upper bound: a proof
	 * may terminate before consuming every candidate. */
	static constexpr uint64_t candidateSetupUnits = 2;
	return static_cast<uint64_t>(this->maximumCandidateCount()) *
	    (candidateSetupUnits + sampleLimit());
    }

    uint64_t progressCompletedUnits(void) const
    {
	const uint64_t total = this->progressTotalUnits();
	if (!total || this->phaseValue == Phase::INACTIVE)
	    return 0;
	if (this->phaseValue == Phase::TERMINAL)
	    return total;
	static constexpr uint64_t candidateSetupUnits = 2;
	const uint64_t candidateUnits = candidateSetupUnits + sampleLimit();
	uint64_t completed =
	    static_cast<uint64_t>(this->totalMeasuredCandidateCountValue) *
	    candidateUnits;
	if (this->phaseValue == Phase::PRESENTING) {
	    completed += 1;
	} else if (this->phaseValue == Phase::MEASURING) {
	    completed += candidateSetupUnits +
		std::min<unsigned int>(this->sampleCountValue, sampleLimit());
	}
	return std::min(completed, total);
    }

    bool awaitingSample(void) const
    {
	/* ALLOCATING waits for occurrence cuts, PRESENTING waits for their first
	 * exact framebuffer, and MEASURING waits for unchanged timing frames.  All
	 * three retain the capacity successor and its presentation edge. */
	return this->phaseValue == Phase::ALLOCATING ||
	    this->phaseValue == Phase::PRESENTING ||
	    this->phaseValue == Phase::MEASURING;
    }

    bool awaitingPresentationFrame(void) const
    {
	return this->phaseValue == Phase::PRESENTING ||
	    this->phaseValue == Phase::MEASURING;
    }

private:
    static size_t saturatingAdd(size_t left, size_t right)
    {
	return right > SIZE_MAX - left ? SIZE_MAX : left + right;
    }

    static size_t strictUpperBound(size_t demand, size_t ceiling)
    {
	const size_t maximum = std::min(demand, ceiling);
	return maximum == SIZE_MAX ? SIZE_MAX : maximum + 1;
    }

    size_t goalMaximumBudget(void) const
    {
	return this->unsafeBudgetValue == SIZE_MAX ?
	    this->keyValue.maximumCandidateBudget :
	    this->unsafeBudgetValue - 1;
    }

    size_t goalMinimumBudget(void) const
    {
	return std::min(this->keyValue.minimumBudget,
	    this->keyValue.maximumCandidateBudget);
    }

    static size_t proposedBudget(size_t presentedCost,
	uint64_t targetNanoseconds, uint64_t observedNanoseconds,
	size_t demand, size_t candidate, bool metTarget)
    {
	/* Leave twenty percent timing headroom.  The endpoint deadline remains
	 * the independent guard for non-linear renderer costs. */
	static constexpr uint64_t headroomNumerator = 4;
	static constexpr uint64_t headroomDenominator = 5;
	const long double proposal =
	    static_cast<long double>(presentedCost) *
	    static_cast<long double>(targetNanoseconds) *
	    static_cast<long double>(headroomNumerator) /
	    (static_cast<long double>(observedNanoseconds) *
		static_cast<long double>(headroomDenominator));
	const size_t growthLimit = candidate > SIZE_MAX / 4 ?
	    SIZE_MAX : candidate * 4;
	size_t result = proposal >= static_cast<long double>(SIZE_MAX) ?
	    SIZE_MAX : static_cast<size_t>(std::max<long double>(1.0L,
		proposal));
	result = std::min(result, growthLimit);
	result = std::min(result, demand);
	/* A safe candidate is itself stronger evidence than a conservative
	 * extrapolation from its fixed per-frame overhead. */
	return metTarget ? std::max(candidate, result) : result;
    }

    static size_t median(std::array<size_t, SampleLimit> values)
    {
	std::sort(values.begin(), values.end());
	return values[values.size() / 2];
    }

    void begin(const Observation &observation)
    {
	this->reset();
	this->keyValue = observation.key;
	this->goalValue = observation.key.startAtStatic ?
	    Goal::STATIC : Goal::STEADY;
	this->unsafeBudgetValue = strictUpperBound(
	    observation.key.maximumCandidateBudget,
	    this->goalValue == Goal::STATIC ?
		observation.key.maximumBudgetCeiling :
		observation.key.preferredBudgetCeiling);
	const size_t maximum = this->goalMaximumBudget();
	const size_t minimum = this->goalMinimumBudget();
	if (maximum < minimum) {
	    if (this->goalValue == Goal::STEADY &&
		observation.key.maximumTargetNanoseconds >
		    observation.key.preferredTargetNanoseconds)
		(void)this->advanceToStaticGoal();
	    else
		(void)this->finish(Result::CERTIFIED);
	    return;
	}
	const size_t candidate = std::max(minimum,
	    std::min(observation.candidateBudget, maximum));
	const size_t knownSafe = std::min(
	    observation.knownSafeBudget, candidate);
	this->safeBudgetValue = knownSafe >= minimum ? knownSafe : 0;
	const bool allocationPending =
	    candidate != observation.candidateBudget;
	const Phase initialPhase = allocationPending ? Phase::ALLOCATING :
	    (observation.validSample ? Phase::MEASURING : Phase::PRESENTING);
	this->startCandidate(candidate,
	    observation.validSample ? observation.presentedCost : 0,
	    initialPhase);
    }

    void startCandidate(size_t budget, size_t presentedCost,
	Phase phase = Phase::MEASURING)
    {
	this->phaseValue = phase;
	this->candidateBudgetValue = budget;
	this->candidatePresentedCostValue = presentedCost;
	this->candidatePopulationSignatureValue = 0;
	this->candidatePopulationMinimumBudgetValue = 0;
	this->candidateNextDistinctPopulationBudgetValue = 0;
	this->candidateNextDistinctPopulationBudgetKnownValue = false;
	this->sampleCountValue = 0;
	this->safeSampleCountValue = 0;
	this->maximumSafeSampleCountValue = 0;
	this->invalidSampleCountValue = 0;
	this->proposalValues.fill(0);
    }

    bool bindCandidatePopulationInterval(const Observation &observation)
    {
	if (observation.candidateBudget != this->candidateBudgetValue)
	    return true;
	if (observation.populationMinimumBudget) {
	    if (observation.populationMinimumBudget >
		    this->candidateBudgetValue ||
		(this->candidatePopulationMinimumBudgetValue &&
		 this->candidatePopulationMinimumBudgetValue !=
		    observation.populationMinimumBudget))
		return false;
	    this->candidatePopulationMinimumBudgetValue =
		observation.populationMinimumBudget;
	}
	if (!observation.nextDistinctPopulationBudgetKnown)
	    return true;
	if (observation.nextDistinctPopulationBudget &&
	    observation.nextDistinctPopulationBudget <=
		this->candidateBudgetValue)
	    return false;
	if (this->candidateNextDistinctPopulationBudgetKnownValue &&
	    this->candidateNextDistinctPopulationBudgetValue !=
		observation.nextDistinctPopulationBudget)
	    return false;
	this->candidateNextDistinctPopulationBudgetValue =
	    observation.nextDistinctPopulationBudget;
	this->candidateNextDistinctPopulationBudgetKnownValue = true;
	return true;
    }

    bool measured(size_t budget) const
    {
	for (unsigned int i = 0; i < this->measuredCandidateCountValue; ++i) {
	    if (this->measuredCandidatesValue[i] == budget)
		return true;
	}
	return false;
    }

    int measuredPopulation(uint64_t signature, size_t presentedCost) const
    {
	if (!signature || !presentedCost)
	    return -1;
	for (unsigned int i = 0; i < this->measuredCandidateCountValue; ++i)
	    if (this->measuredPopulationSignaturesValue[i] == signature &&
		this->measuredPopulationCostsValue[i] == presentedCost)
		return static_cast<int>(i);
	return -1;
    }

    size_t unmeasuredCandidate(size_t preferred) const
    {
	const size_t upper = this->goalMaximumBudget();
	const size_t lower = std::max(this->goalMinimumBudget(),
	    saturatingAdd(this->safeBudgetValue, 1));
	if (upper < lower)
	    return 0;
	preferred = std::max(lower, std::min(preferred, upper));
	if (!this->measured(preferred))
	    return preferred;

	const size_t midpoint = this->safeBudgetValue +
	    (upper - this->safeBudgetValue) / 2;
	if (midpoint > this->safeBudgetValue && !this->measured(midpoint))
	    return midpoint;
	if (!this->measured(upper))
	    return upper;
	return !this->measured(lower) ? lower : 0;
    }

    Decision classifyEquivalentPopulation(unsigned int measuredIndex)
    {
	return this->classifyCandidate(
	    this->measuredPopulationSafeValue[measuredIndex],
	    this->measuredPopulationProposalsValue[measuredIndex], true,
	    this->measuredPopulationMaximumSafeValue[measuredIndex]);
    }

    Decision classifyCandidate(void)
    {
	const bool safe = this->safeSampleCountValue * 2 >= sampleLimit();
	const bool maximumSafe =
	    this->maximumSafeSampleCountValue * 2 >= sampleLimit();
	return this->classifyCandidate(safe, median(this->proposalValues),
	    false, maximumSafe);
    }

    Decision classifyCandidate(bool safe, size_t proposal,
	bool equivalentPopulation = false, bool maximumSafe = false)
    {
	const size_t candidate = this->candidateBudgetValue;
	const size_t populationMinimum =
	    this->candidatePopulationMinimumBudgetValue ?
		this->candidatePopulationMinimumBudgetValue : candidate;
	const size_t populationMaximum =
	    this->candidateNextDistinctPopulationBudgetKnownValue ?
		(this->candidateNextDistinctPopulationBudgetValue ?
		    this->candidateNextDistinctPopulationBudgetValue - 1 :
		    this->keyValue.maximumCandidateBudget) : candidate;
	const size_t classifiedMinimum = std::max(this->goalMinimumBudget(),
	    std::min(populationMinimum, this->goalMaximumBudget()));
	const size_t classifiedMaximum = std::max(classifiedMinimum,
	    std::min(populationMaximum, this->goalMaximumBudget()));
	const unsigned int measuredIndex = this->measuredCandidateCountValue++;
	this->totalMeasuredCandidateCountValue++;
	this->measuredCandidatesValue[measuredIndex] = candidate;
	this->measuredPopulationSignaturesValue[measuredIndex] =
	    this->candidatePopulationSignatureValue;
	this->measuredPopulationCostsValue[measuredIndex] =
	    this->candidatePresentedCostValue;
	this->measuredPopulationSafeValue[measuredIndex] = safe;
	this->measuredPopulationMaximumSafeValue[measuredIndex] = maximumSafe;
	this->measuredPopulationProposalsValue[measuredIndex] = proposal;
	/* The same completed samples classify both ordered deadlines.  When a
	 * candidate misses the preferred cadence but already meets the static hard
	 * deadline, it is the static goal's proven safe visual.  Carry it forward
	 * directly instead of first coarsening to the steady lower bound and then
	 * rebuilding toward the population we just measured. */
	if (!safe && maximumSafe && this->goalValue == Goal::STEADY &&
	    this->keyValue.maximumTargetNanoseconds >
		this->keyValue.preferredTargetNanoseconds) {
	    this->safeBudgetValue = std::max(
		this->safeBudgetValue, classifiedMaximum);
	    return this->advanceToStaticGoal();
	}
	if (safe)
	    this->safeBudgetValue = std::max(
		this->safeBudgetValue, classifiedMaximum);
	else
	    this->unsafeBudgetValue = std::min(
		this->unsafeBudgetValue, classifiedMinimum);

	const size_t goalMaximum = this->goalMaximumBudget();
	const size_t goalMinimum = this->goalMinimumBudget();
	const bool exhaustedAtMinimum =
	    this->unsafeBudgetValue != SIZE_MAX &&
	    this->unsafeBudgetValue <= goalMinimum &&
	    this->safeBudgetValue < goalMinimum;
    /* Capacity is an operational guard, not a quality objective.  Do not
     * search for an adjacent integer budget after establishing a bracket.
     * The final goal may spend one midpoint bridge below; all other brackets
     * settle immediately and let the visual-importance allocator spend the
     * proven budget.  A new semantic epoch may calibrate again from this
     * evidence. */
	const bool usefulBracketEstablished =
	    this->safeBudgetValue > 0 &&
	    this->unsafeBudgetValue <=
		this->keyValue.maximumCandidateBudget;
	const bool finalGoal = this->goalValue == Goal::STATIC ||
	    this->keyValue.maximumTargetNanoseconds ==
		this->keyValue.preferredTargetNanoseconds;
	const bool bracketBridgeAvailable = usefulBracketEstablished &&
	    finalGoal && !this->finalBracketBridgeAttemptedValue &&
	    this->measuredCandidateCountValue < candidateLimit() &&
	    this->unsafeBudgetValue >
		saturatingAdd(this->safeBudgetValue, 1);
	if (this->safeBudgetValue >= goalMaximum ||
	    exhaustedAtMinimum ||
	    (usefulBracketEstablished && !bracketBridgeAvailable) ||
	    this->measuredCandidateCountValue >= candidateLimit() ||
	    (this->unsafeBudgetValue != SIZE_MAX &&
	     this->unsafeBudgetValue <=
		saturatingAdd(this->safeBudgetValue, 1))) {
	    if (this->goalValue == Goal::STEADY &&
		this->keyValue.maximumTargetNanoseconds >
		    this->keyValue.preferredTargetNanoseconds &&
		this->safeBudgetValue <
		    this->keyValue.maximumCandidateBudget)
		return this->advanceToStaticGoal();
	    return this->finish(Result::CERTIFIED);
	}

	/* Every completed sample contributes a conservative throughput proposal.
	 * Use it in both search directions.  Falling back to a midpoint after a
	 * miss discarded that evidence and forced a large discrete mesh through a
	 * visible logarithmic staircase before the independent static goal could
	 * restore useful quality.  unmeasuredCandidate() still clamps the proposal
	 * to the proven bracket and supplies a midpoint when rounding names an
	 * already measured endpoint, so progress and bracket soundness are
	 * unchanged. */
	/* A numeric budget which selects an already measured population adds no
	 * new timing evidence.  Reusing that population's proposal can name the
	 * immediately adjacent integer budget, which the discrete allocator may
	 * map to the same cuts again.  Move to the middle of the remaining proven
	 * bracket instead.  This preserves safety while ensuring every reused
	 * population consumes a meaningful part of the finite search interval. */
	size_t next = bracketBridgeAvailable ?
	    this->safeBudgetValue +
		(this->unsafeBudgetValue - this->safeBudgetValue) / 2 :
	    equivalentPopulation ?
	    this->safeBudgetValue +
		(this->goalMaximumBudget() - this->safeBudgetValue) / 2 :
	    proposal;
	next = this->unmeasuredCandidate(next);
	if (!next)
	    return this->finish(Result::CERTIFIED);
	if (bracketBridgeAvailable)
	    this->finalBracketBridgeAttemptedValue = true;

	this->startCandidate(next, 0, Phase::ALLOCATING);
	return this->decision(Result::REALLOCATE, next);
    }

    Decision advanceToStaticGoal(void)
    {
	/* Preserve the population proven at the preferred cadence and reset only
	 * the per-goal unsafe bracket and candidate history.  This is the sole
	 * legal goal transition and directly implements AdvanceToStaticGoal in
	 * ObolCapacitySearch.tla. */
	const size_t preferredProposal = median(this->proposalValues);
	this->goalValue = Goal::STATIC;
	this->unsafeBudgetValue = strictUpperBound(
	    this->keyValue.maximumCandidateBudget,
	    this->keyValue.maximumBudgetCeiling);
	this->measuredCandidateCountValue = 0;
	this->measuredCandidatesValue.fill(0);
	this->finalBracketBridgeAttemptedValue = false;
	const size_t upper = this->goalMaximumBudget();
	const size_t lower = std::max(this->goalMinimumBudget(),
	    saturatingAdd(this->safeBudgetValue, 1));
	if (upper < lower)
	    return this->finish(Result::CERTIFIED);
	/* The steady samples already measured fixed per-frame and per-occurrence
	 * cost.  Scale their conservative proposal by the deadline ratio instead
	 * of bisecting the entire pixel-demand range.  A midpoint can be orders of
	 * magnitude larger for a few detailed meshes and wastes the first static
	 * attempt on a predictable hard-deadline abort. */
	size_t proposed = lower;
	if (preferredProposal > 0 &&
	    this->keyValue.preferredTargetNanoseconds > 0) {
	    const long double scaled =
		static_cast<long double>(preferredProposal) *
		static_cast<long double>(
		    this->keyValue.maximumTargetNanoseconds) /
		static_cast<long double>(
		    this->keyValue.preferredTargetNanoseconds);
	    proposed = scaled >= static_cast<long double>(SIZE_MAX) ?
		SIZE_MAX : static_cast<size_t>(scaled);
	}
	const size_t candidate = std::max(lower, std::min(proposed, upper));
	this->startCandidate(candidate, 0, Phase::ALLOCATING);
	return this->decision(Result::REALLOCATE, candidate);
    }

    Decision finish(Result result)
    {
	this->phaseValue = Phase::TERMINAL;
	this->terminalResultValue = result;
	this->candidateBudgetValue = 0;
	this->candidatePresentedCostValue = 0;
	this->candidatePopulationSignatureValue = 0;
	this->candidatePopulationMinimumBudgetValue = 0;
	this->candidateNextDistinctPopulationBudgetValue = 0;
	this->candidateNextDistinctPopulationBudgetKnownValue = false;
	this->sampleCountValue = 0;
	this->safeSampleCountValue = 0;
	this->maximumSafeSampleCountValue = 0;
	this->invalidSampleCountValue = 0;
	return this->decision(result, this->safeBudgetValue);
    }

    Decision decision(Result result, size_t budget) const
    {
	Decision value;
	value.result = result;
	value.budget = budget;
	value.safeBudget = this->safeBudgetValue;
	value.unsafeBudget = this->unsafeBudgetValue;
	value.samplesRemaining = this->samplesRemaining();
	value.measuredCandidates = this->measuredCandidateCountValue;
	value.totalMeasuredCandidates = this->totalMeasuredCandidateCountValue;
	return value;
    }

    BObolLodCapacitySearchKey keyValue;
    Goal goalValue = Goal::STEADY;
    Phase phaseValue = Phase::INACTIVE;
    Result terminalResultValue = Result::NONE;
    size_t safeBudgetValue = 0;
    size_t unsafeBudgetValue = SIZE_MAX;
    size_t candidateBudgetValue = 0;
    size_t candidatePresentedCostValue = 0;
    uint64_t candidatePopulationSignatureValue = 0;
    size_t candidatePopulationMinimumBudgetValue = 0;
    size_t candidateNextDistinctPopulationBudgetValue = 0;
    bool candidateNextDistinctPopulationBudgetKnownValue = false;
    unsigned int sampleCountValue = 0;
    unsigned int safeSampleCountValue = 0;
    unsigned int maximumSafeSampleCountValue = 0;
    unsigned int invalidSampleCountValue = 0;
    unsigned int measuredCandidateCountValue = 0;
    unsigned int totalMeasuredCandidateCountValue = 0;
    bool finalBracketBridgeAttemptedValue = false;
    std::array<size_t, SampleLimit> proposalValues {{0, 0, 0}};
    std::array<size_t, CandidateLimit> measuredCandidatesValue {{}};
    std::array<uint64_t, CandidateLimit>
	measuredPopulationSignaturesValue {{}};
    std::array<size_t, CandidateLimit> measuredPopulationCostsValue {{}};
    std::array<bool, CandidateLimit> measuredPopulationSafeValue {{}};
    std::array<bool, CandidateLimit>
	measuredPopulationMaximumSafeValue {{}};
    std::array<size_t, CandidateLimit> measuredPopulationProposalsValue {{}};
};

static_assert(std::is_trivially_copyable<
	BObolLodCapacitySearchCertificate>::value,
    "capacity-search certificates must remain allocation-free values");

#endif /* LIBBOBOL_LOD_CAPACITY_SEARCH_PRIVATE_H */
