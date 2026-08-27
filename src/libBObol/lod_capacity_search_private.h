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
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

/*
 * Capacity evidence is the output of a search, so the search cannot use the
 * pipeline's changing capacity-output epoch as one of its own inputs.  The
 * other four semantic epochs, preferred and maximum target cadences, and
 * complete pixel-demand cost identify the frozen problem.  Publishing a
 * successor or terminal result advances the ordinary capacity epoch outside
 * this value.
 */
struct BObolLodCapacitySearchKey {
    BObolLodInventoryEpoch inventory;
    BObolLodAvailabilityEpoch availability;
    BObolLodViewEpoch view;
    BObolLodPolicyEpoch policy;
    uint64_t preferredTargetNanoseconds = 0;
    uint64_t maximumTargetNanoseconds = 0;
    size_t demandedBudget = 0;
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

    bool same(const BObolLodCapacitySearchKey &other) const
    {
	return this->inventory == other.inventory &&
	    this->availability == other.availability &&
	    this->view == other.view && this->policy == other.policy &&
	    this->preferredTargetNanoseconds ==
		other.preferredTargetNanoseconds &&
	    this->maximumTargetNanoseconds ==
		other.maximumTargetNanoseconds &&
	    this->demandedBudget == other.demandedBudget &&
	    this->minimumBudget == other.minimumBudget &&
	    this->preferredBudgetCeiling ==
		other.preferredBudgetCeiling &&
	    this->maximumBudgetCeiling == other.maximumBudgetCeiling;
    }

    bool valid(void) const
    {
	return this->preferredTargetNanoseconds > 0 &&
	    this->maximumTargetNanoseconds >=
		this->preferredTargetNanoseconds &&
	    this->demandedBudget > 0 &&
	    this->minimumBudget <= this->demandedBudget;
    }
};

/*
 * Bounded search over complete, visual-importance-ordered scene allocations.
 * A candidate is an allocation budget, not a triangle or object ordinal.  It
 * consumes exactly sampleLimit valid unchanged frames, then strictly moves a
 * known-safe lower bound or a known-unsafe upper bound.  At most
 * candidateLimit distinct budgets may be measured for each of the key's two
 * ordered goals.  The optional static goal inherits the steady-safe lower
 * bound and cannot transition back.
 *
 * Throughput estimates are deliberately confined to proposing the next
 * candidate.  Only completed-frame classifications change the bracket.  This
 * is the executable counterpart of ObolCapacitySearch.tla.
 */
class BObolLodCapacitySearchCertificate {
private:
    static constexpr unsigned int SampleLimit = 3;
    static constexpr unsigned int CandidateLimit = 8;
    static constexpr unsigned int InvalidSampleLimit = 2;

public:
    enum class Goal : uint8_t {
	STEADY = 0,
	STATIC
    };

    enum class Phase : uint8_t {
	INACTIVE = 0,
	ALLOCATING,
	MEASURING,
	TERMINAL
    };

    enum class Result : uint8_t {
	NONE = 0,
	REQUEST_SAMPLE,
	REALLOCATE,
	CERTIFIED,
	UNMEASURABLE,
	STALE_POPULATION
    };

    struct Observation {
	BObolLodCapacitySearchKey key;
	size_t candidateBudget = 0;
	size_t presentedCost = 0;
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
    static constexpr unsigned int invalidSampleLimit(void)
    {
	return InvalidSampleLimit;
    }

    static BObolLodCapacitySearchKey keyFor(
	const BObolLodAdmissionRevisionStamp &stamp,
	uint64_t preferredTargetNanoseconds,
	uint64_t maximumTargetNanoseconds, size_t demandedBudget,
	size_t minimumBudget = 0)
    {
	BObolLodCapacitySearchKey key;
	key.inventory = stamp.inventory;
	key.availability = stamp.availability;
	key.view = stamp.view;
	key.policy = stamp.policy;
	key.preferredTargetNanoseconds = preferredTargetNanoseconds;
	key.maximumTargetNanoseconds = maximumTargetNanoseconds;
	key.demandedBudget = demandedBudget;
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
	if (began && this->phaseValue == Phase::ALLOCATING)
	    return this->decision(Result::REALLOCATE,
		this->candidateBudgetValue);

	const size_t candidate = std::max(this->goalMinimumBudget(),
	    std::min(observation.candidateBudget,
		this->goalMaximumBudget()));
	if (candidate != this->candidateBudgetValue) {
	    if (this->measured(candidate))
		return this->finish(Result::STALE_POPULATION);
	    this->startCandidate(candidate, observation.presentedCost);
	} else if (this->phaseValue == Phase::ALLOCATING) {
	    /* A reallocation decision names the next complete population, but it
	     * does not own a framebuffer sample.  Only the successor observation
	     * proves that the population was applied and transfers the certificate
	     * into its finite measurement phase. */
	    this->phaseValue = Phase::MEASURING;
	    this->candidatePresentedCostValue = observation.presentedCost;
	} else if (this->candidatePresentedCostValue &&
	    observation.presentedCost &&
	    this->candidatePresentedCostValue != observation.presentedCost) {
	    return this->finish(Result::STALE_POPULATION);
	} else if (!this->candidatePresentedCostValue) {
	    this->candidatePresentedCostValue = observation.presentedCost;
	}
	return this->decision(Result::REQUEST_SAMPLE,
	    this->candidateBudgetValue);
    }

    Decision observe(const Observation &observation)
    {
	const Decision prepared = this->prepare(observation);
	if (prepared.result != Result::REQUEST_SAMPLE)
	    return prepared;

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
	this->proposalValues[this->sampleCountValue] = proposedBudget(
	    observation.presentedCost, this->activeTargetNanoseconds(),
	    observation.observedNanoseconds,
	    observation.key.demandedBudget, this->candidateBudgetValue,
	    metTarget);
	this->sampleCountValue++;
	if (this->sampleCountValue < sampleLimit())
	    return this->decision(Result::REQUEST_SAMPLE,
		this->candidateBudgetValue);

	return this->classifyCandidate();
    }

    Phase phase(void) const { return this->phaseValue; }
    Goal goal(void) const { return this->goalValue; }
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
    unsigned int samplesRemaining(void) const
    {
	return this->phaseValue == Phase::MEASURING &&
	    this->sampleCountValue < sampleLimit() ?
	    sampleLimit() - this->sampleCountValue : 0;
    }

    bool awaitingSample(void) const
    {
	return this->phaseValue == Phase::MEASURING;
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
	    this->keyValue.demandedBudget : this->unsafeBudgetValue - 1;
    }

    size_t goalMinimumBudget(void) const
    {
	return std::min(this->keyValue.minimumBudget,
	    this->keyValue.demandedBudget);
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
	this->goalValue = Goal::STEADY;
	this->unsafeBudgetValue = strictUpperBound(
	    observation.key.demandedBudget,
	    observation.key.preferredBudgetCeiling);
	const size_t maximum = this->goalMaximumBudget();
	const size_t minimum = this->goalMinimumBudget();
	if (maximum < minimum) {
	    if (observation.key.maximumTargetNanoseconds >
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
	const bool allocationPending = candidate != observation.candidateBudget;
	this->startCandidate(candidate,
	    allocationPending ? 0 : observation.presentedCost,
	    allocationPending);
    }

    void startCandidate(size_t budget, size_t presentedCost,
	bool allocationPending = false)
    {
	this->phaseValue = allocationPending ?
	    Phase::ALLOCATING : Phase::MEASURING;
	this->candidateBudgetValue = budget;
	this->candidatePresentedCostValue = presentedCost;
	this->sampleCountValue = 0;
	this->safeSampleCountValue = 0;
	this->invalidSampleCountValue = 0;
	this->proposalValues.fill(0);
    }

    bool measured(size_t budget) const
    {
	for (unsigned int i = 0; i < this->measuredCandidateCountValue; ++i) {
	    if (this->measuredCandidatesValue[i] == budget)
		return true;
	}
	return false;
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

    Decision classifyCandidate(void)
    {
	const size_t candidate = this->candidateBudgetValue;
	this->measuredCandidatesValue[this->measuredCandidateCountValue++] =
	    candidate;
	const bool safe = this->safeSampleCountValue * 2 >= sampleLimit();
	if (safe)
	    this->safeBudgetValue = std::max(this->safeBudgetValue, candidate);
	else
	    this->unsafeBudgetValue = std::min(
		this->unsafeBudgetValue, candidate);

	const size_t goalMaximum = this->goalMaximumBudget();
	const size_t goalMinimum = this->goalMinimumBudget();
	const bool exhaustedAtMinimum =
	    this->unsafeBudgetValue != SIZE_MAX &&
	    this->unsafeBudgetValue <= goalMinimum &&
	    this->safeBudgetValue < goalMinimum;
	if (this->safeBudgetValue >= goalMaximum ||
	    exhaustedAtMinimum ||
	    this->measuredCandidateCountValue >= candidateLimit() ||
	    (this->unsafeBudgetValue != SIZE_MAX &&
	     this->unsafeBudgetValue <=
		saturatingAdd(this->safeBudgetValue, 1))) {
	    if (this->goalValue == Goal::STEADY &&
		this->keyValue.maximumTargetNanoseconds >
		    this->keyValue.preferredTargetNanoseconds &&
		this->safeBudgetValue < this->keyValue.demandedBudget)
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
	const size_t proposal = median(this->proposalValues);
	size_t next = proposal;
	next = this->unmeasuredCandidate(next);
	if (!next)
	    return this->finish(Result::CERTIFIED);

	this->startCandidate(next, 0, true);
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
	    this->keyValue.demandedBudget,
	    this->keyValue.maximumBudgetCeiling);
	this->measuredCandidateCountValue = 0;
	this->measuredCandidatesValue.fill(0);
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
	this->startCandidate(candidate, 0, true);
	return this->decision(Result::REALLOCATE, candidate);
    }

    Decision finish(Result result)
    {
	this->phaseValue = Phase::TERMINAL;
	this->terminalResultValue = result;
	this->candidateBudgetValue = 0;
	this->candidatePresentedCostValue = 0;
	this->sampleCountValue = 0;
	this->safeSampleCountValue = 0;
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
    unsigned int sampleCountValue = 0;
    unsigned int safeSampleCountValue = 0;
    unsigned int invalidSampleCountValue = 0;
    unsigned int measuredCandidateCountValue = 0;
    std::array<size_t, SampleLimit> proposalValues {{0, 0, 0}};
    std::array<size_t, CandidateLimit> measuredCandidatesValue {{}};
};

static_assert(std::is_trivially_copyable<
	BObolLodCapacitySearchCertificate>::value,
    "capacity-search certificates must remain allocation-free values");

#endif /* LIBBOBOL_LOD_CAPACITY_SEARCH_PRIVATE_H */
