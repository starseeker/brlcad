/*       L O D _ P R O G R E S S _ E S T I M A T O R _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_PROGRESS_ESTIMATOR_PRIVATE_H
#define LIBBOBOL_LOD_PROGRESS_ESTIMATOR_PRIVATE_H

#include "common.h"

#include "lod_revision_private.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>

/*
 * Observed-time projection of the finite convergence ranks.  This is a HUD
 * estimator, not a scheduler: no production decision may consume its output.
 * Ordinary ranks use observed unit rates.  Visible refinement instead uses
 * conservative completed-frame cycles so fast result publication cannot imply
 * a completion time which ignores presentation and replanning.  Unknown work
 * or an unqualified refinement frontier makes the result indeterminate.
 */
class BObolLodProgressEstimator {
public:
    enum class Rank : size_t {
	DISCOVERY = 0,
	SOURCE_PREPARATION,
	RENDERER_PREPARATION,
	VISIBLE_RESOLUTION,
	CAPACITY_SEARCH,
	COUNT
    };

    struct WorkRank {
	bool present = false;
	uint64_t completed = 0;
	uint64_t total = 0;
    };

    struct Inputs {
	uint64_t episodeRevision = 0;
	/* Admission policy epochs identify internally selected refinement tiers.
	 * They reset only the cycle forecast; the user-visible episode is owned by
	 * episodeRevision. */
	BObolLodPolicyEpoch refinementTierEpoch;
	uint64_t renderCompletionSerial = 0;
	int64_t episodeStartMicroseconds = 0;
	int64_t observationMicroseconds = 0;
	bool terminal = false;
	bool unknownForegroundWork = false;
	uint64_t finalPresentationMicroseconds = 0;
	std::array<WorkRank, static_cast<size_t>(Rank::COUNT)> ranks;

	WorkRank &rank(Rank value)
	{
	    return this->ranks[static_cast<size_t>(value)];
	}

	const WorkRank &rank(Rank value) const
	{
	    return this->ranks[static_cast<size_t>(value)];
	}
    };

    struct Estimate {
	bool available = false;
	bool refinementCycleBased = false;
	float fraction = 0.0f;
	uint64_t remainingMilliseconds = 0;
	uint64_t remainingRefinementCycles = 0;
    };

    Estimate evaluate(const Inputs &inputs)
    {
	const bool newEpisode =
	    this->episodeRevision != inputs.episodeRevision;
	if (newEpisode)
	    this->beginEpisode(inputs);
	const bool refinementTierChanged = !newEpisode &&
	    this->refinementTierEpoch != inputs.refinementTierEpoch;
	if (refinementTierChanged) {
	    this->refinementTierEpoch = inputs.refinementTierEpoch;
	    this->fractionFloor = 0.0f;
	    this->cycleForecast.reset();
	}
	const bool terminalChanged = inputs.terminal != this->lastTerminal;
	this->lastTerminal = inputs.terminal;

	if (inputs.terminal) {
	    this->lastEstimate.available = true;
	    this->lastEstimate.refinementCycleBased = false;
	    this->lastEstimate.fraction = 1.0f;
	    this->lastEstimate.remainingMilliseconds = 0;
	    this->lastEstimate.remainingRefinementCycles = 0;
	    return this->lastEstimate;
	}

	bool rankChanged = false;
	for (size_t i = 0; i < this->rates.size(); ++i)
	    rankChanged = this->rates[i].observe(inputs.ranks[i],
		inputs.observationMicroseconds, this->episodeStartMicroseconds) ||
		rankChanged;
	const CycleForecast::Estimate cycleEstimate =
	    this->cycleForecast.observe(inputs.refinementTierEpoch,
		inputs.renderCompletionSerial,
		inputs.rank(Rank::VISIBLE_RESOLUTION),
		inputs.observationMicroseconds);
	if (cycleEstimate.confidenceReset)
	    this->fractionFloor = 0.0f;
	const bool estimateInputsChanged = newEpisode || terminalChanged ||
	    refinementTierChanged || rankChanged || cycleEstimate.changed ||
	    inputs.unknownForegroundWork != this->unknownForegroundWork ||
	    inputs.finalPresentationMicroseconds !=
		this->finalPresentationMicroseconds;
	this->unknownForegroundWork = inputs.unknownForegroundWork;
	this->finalPresentationMicroseconds =
	    inputs.finalPresentationMicroseconds;
	if (!estimateInputsChanged)
	    return this->lastEstimate;

	long double remainingMicroseconds = 0.0L;
	bool hasIncompleteRank = false;
	bool allIncompleteRanksMeasured = !inputs.unknownForegroundWork;
	bool refinementCycleBased = false;
	uint64_t remainingRefinementCycles = 0;
	for (size_t i = 0; i < this->rates.size(); ++i) {
	    const WorkRank &rank = inputs.ranks[i];
	    if (!rank.present || rank.total == 0)
		continue;
	    const uint64_t completed = std::min(rank.completed, rank.total);
	    if (completed == rank.total)
		continue;
	    hasIncompleteRank = true;
	    if (i == static_cast<size_t>(Rank::VISIBLE_RESOLUTION)) {
		if (!cycleEstimate.available) {
		    allIncompleteRanksMeasured = false;
		    continue;
		}
		remainingMicroseconds += cycleEstimate.remainingMicroseconds;
		refinementCycleBased = true;
		remainingRefinementCycles = cycleEstimate.remainingCycles;
		continue;
	    }
	    long double rate = this->rates[i].microsecondsPerUnit();
	    /* Capacity-search units include a bounded presentation and timing
	     * sample.  The current completed-frame duration is a useful seed before
	     * the first candidate advances; unlike a general source unit, its
	     * physical meaning is known. */
	    if (rate <= 0.0L && i == static_cast<size_t>(Rank::CAPACITY_SEARCH))
		rate = static_cast<long double>(
		    inputs.finalPresentationMicroseconds);
	    if (rate <= 0.0L) {
		allIncompleteRanksMeasured = false;
		continue;
	    }
	    remainingMicroseconds += static_cast<long double>(
		rank.total - completed) * rate;
	}

	/* A completed-frame cycle already includes its presentation.  Other ranks
	 * still need one final exact frame after their last unit is published. */
	if (inputs.finalPresentationMicroseconds > 0 && !refinementCycleBased)
	    remainingMicroseconds += static_cast<long double>(
		inputs.finalPresentationMicroseconds);
	const bool hasEstimatedWork = hasIncompleteRank ||
	    inputs.finalPresentationMicroseconds > 0;
	if (!allIncompleteRanksMeasured || !hasEstimatedWork) {
	    this->lastEstimate = Estimate();
	    return this->lastEstimate;
	}

	const int64_t elapsed = std::max<int64_t>(0,
	    inputs.observationMicroseconds - this->episodeStartMicroseconds);
	const long double estimatedTotal =
	    static_cast<long double>(elapsed) + remainingMicroseconds;
	float fraction = estimatedTotal > 0.0L ? static_cast<float>(
	    static_cast<long double>(elapsed) / estimatedTotal) : 0.0f;
	/* One coherent presentation is still outstanding, so only the terminal
	 * convergence decision may publish 100 percent. */
	static constexpr float maximumActiveFraction = 0.99f;
	fraction = std::max(this->fractionFloor,
	    std::min(maximumActiveFraction, fraction));
	this->fractionFloor = fraction;
	this->lastEstimate.available = true;
	this->lastEstimate.refinementCycleBased = refinementCycleBased;
	this->lastEstimate.fraction = fraction;
	this->lastEstimate.remainingRefinementCycles =
	    remainingRefinementCycles;
	this->lastEstimate.remainingMilliseconds =
	    remainingMicroseconds >= static_cast<long double>(UINT64_MAX) *
		microsecondsPerMillisecond ? UINT64_MAX :
	    static_cast<uint64_t>((remainingMicroseconds +
		microsecondsPerMillisecond - 1.0L) /
		microsecondsPerMillisecond);
	return this->lastEstimate;
    }

    void resetEpisode(void)
    {
	this->episodeRevision = 0;
	this->refinementTierEpoch.reset();
	this->episodeStartMicroseconds = 0;
	this->unknownForegroundWork = false;
	this->finalPresentationMicroseconds = 0;
	this->lastTerminal = false;
	this->fractionFloor = 0.0f;
	this->lastEstimate = Estimate();
	this->cycleForecast.reset();
	for (Rate &rate : this->rates)
	    rate.resetObservation();
    }

private:
    class CycleForecast {
    public:
	struct Estimate {
	    bool changed = false;
	    bool confidenceReset = false;
	    bool available = false;
	    uint64_t remainingCycles = 0;
	    long double remainingMicroseconds = 0.0L;
	};

	Estimate observe(BObolLodPolicyEpoch tier, uint64_t renderSerial,
	    const WorkRank &rank, int64_t nowMicroseconds)
	{
	    Estimate result;
	    if (!rank.present || rank.total == 0 ||
		std::min(rank.completed, rank.total) == rank.total) {
		result.changed = this->active;
		this->reset();
		return result;
	    }

	    const uint64_t completed = std::min(rank.completed, rank.total);
	    const uint64_t unresolved = rank.total - completed;
	    if (!this->active || this->tier != tier ||
		this->totalUnits != rank.total) {
		result.confidenceReset = true;
		this->begin(tier, renderSerial, rank.total, unresolved,
		    nowMicroseconds);
		result.changed = true;
		return result;
	    }

	    /* Newly discovered debt invalidates every projection immediately, even
	     * when it appears between completed frames.  The next completed frame
	     * establishes a fresh comparable baseline. */
	    if (unresolved > this->lastObservedUnresolved) {
		this->resetSamples();
		this->setBoundary(renderSerial, unresolved, nowMicroseconds);
		this->lastObservedUnresolved = unresolved;
		result.changed = true;
		result.confidenceReset = true;
		return result;
	    }
	    this->lastObservedUnresolved = unresolved;

	    if (!renderSerial)
		return this->estimate(unresolved, result);
	    if (!this->boundarySerial) {
		this->setBoundary(renderSerial, unresolved, nowMicroseconds);
		result.changed = true;
		return result;
	    }
	    if (renderSerial != this->boundarySerial) {
		const bool consecutive = this->boundarySerial != UINT64_MAX &&
		    renderSerial == this->boundarySerial + 1;
		const int64_t duration = nowMicroseconds -
		    this->boundaryMicroseconds;
		if (!consecutive || duration <= 0 ||
		    unresolved >= this->boundaryUnresolved) {
		    this->resetSamples();
		    result.confidenceReset = true;
		} else {
		    this->appendSample(this->boundaryUnresolved - unresolved,
			duration);
		}
		this->setBoundary(renderSerial, unresolved, nowMicroseconds);
		result.changed = true;
	    }
	    return this->estimate(unresolved, result);
	}

	void reset(void)
	{
	    this->active = false;
	    this->tier.reset();
	    this->totalUnits = 0;
	    this->lastObservedUnresolved = 0;
	    this->boundarySerial = 0;
	    this->boundaryUnresolved = 0;
	    this->boundaryMicroseconds = 0;
	    this->resetSamples();
	}

    private:
	void begin(BObolLodPolicyEpoch newTier, uint64_t renderSerial,
	    uint64_t total, uint64_t unresolved, int64_t nowMicroseconds)
	{
	    this->reset();
	    this->active = true;
	    this->tier = newTier;
	    this->totalUnits = total;
	    this->lastObservedUnresolved = unresolved;
	    this->setBoundary(renderSerial, unresolved, nowMicroseconds);
	}

	void setBoundary(uint64_t renderSerial, uint64_t unresolved,
	    int64_t nowMicroseconds)
	{
	    this->boundarySerial = renderSerial;
	    this->boundaryUnresolved = unresolved;
	    this->boundaryMicroseconds = nowMicroseconds;
	}

	void resetSamples(void)
	{
	    this->sampleCount = 0;
	    this->reductions.fill(0);
	    this->durations.fill(0);
	}

	void appendSample(uint64_t reduction, int64_t duration)
	{
	    if (!reduction || duration <= 0)
		return;
	    if (this->sampleCount < sampleCapacity) {
		this->reductions[this->sampleCount] = reduction;
		this->durations[this->sampleCount] =
		    static_cast<uint64_t>(duration);
		this->sampleCount++;
		return;
	    }
	    for (size_t i = 1; i < sampleCapacity; ++i) {
		this->reductions[i - 1] = this->reductions[i];
		this->durations[i - 1] = this->durations[i];
	    }
	    this->reductions[sampleCapacity - 1] = reduction;
	    this->durations[sampleCapacity - 1] =
		static_cast<uint64_t>(duration);
	}

	Estimate estimate(uint64_t unresolved, Estimate result) const
	{
	    if (this->sampleCount < minimumQualifiedSamples)
		return result;
	    long double totalReduction = 0.0L;
	    long double totalDuration = 0.0L;
	    for (size_t i = 0; i < this->sampleCount; ++i) {
		totalReduction += static_cast<long double>(
		    this->reductions[i]);
		totalDuration += static_cast<long double>(
		    this->durations[i]);
	    }
	    if (totalReduction <= 0.0L || totalDuration <= 0.0L)
		return result;

	    /* Preserve the measured relationship between work and elapsed time.
	     * Combining the smallest reduction with the longest duration, when they
	     * came from different frames, manufactured a rate which had never been
	     * observed.  On a recorded 37k-item frontier that amplified a roughly
	     * 20-second tail into a 15-minute forecast.  Aggregate paired samples
	     * instead, retaining a modest margin for the usual convergence taper. */
	    static constexpr long double safetyNumerator = 5.0L;
	    static constexpr long double safetyDenominator = 4.0L;
	    const long double safety = safetyNumerator / safetyDenominator;
	    const long double projectedCycles =
		static_cast<long double>(unresolved) *
		static_cast<long double>(this->sampleCount) /
		totalReduction * safety;
	    if (projectedCycles >= static_cast<long double>(UINT64_MAX)) {
		result.remainingCycles = UINT64_MAX;
	    } else {
		result.remainingCycles = static_cast<uint64_t>(projectedCycles);
		if (static_cast<long double>(result.remainingCycles) <
		    projectedCycles)
		    result.remainingCycles++;
	    }
	    result.remainingMicroseconds =
		static_cast<long double>(unresolved) * totalDuration /
		totalReduction * safety;
	    result.available = true;
	    return result;
	}

	static constexpr size_t sampleCapacity = 8;
	static constexpr size_t minimumQualifiedSamples = 3;
	bool active = false;
	BObolLodPolicyEpoch tier;
	uint64_t totalUnits = 0;
	uint64_t lastObservedUnresolved = 0;
	uint64_t boundarySerial = 0;
	uint64_t boundaryUnresolved = 0;
	int64_t boundaryMicroseconds = 0;
	size_t sampleCount = 0;
	std::array<uint64_t, sampleCapacity> reductions = {};
	std::array<uint64_t, sampleCapacity> durations = {};
    };

    class Rate {
    public:
	bool observe(const WorkRank &rank, int64_t nowMicroseconds,
	    int64_t episodeStartMicroseconds)
	{
	    if (!rank.present || rank.total == 0) {
		const bool changed = this->observing;
		this->resetObservation();
		return changed;
	    }

	    const uint64_t completed = std::min(rank.completed, rank.total);
	    if (!this->observing || rank.total != this->totalUnits ||
		completed < this->completedUnits) {
		this->observing = true;
		this->totalUnits = rank.total;
		this->completedUnits = completed;
		this->sampleCompleted = completed;
		this->sampleMicroseconds = nowMicroseconds;
		if (!this->rateAvailable && completed > 0 &&
		    episodeStartMicroseconds > 0 &&
		    nowMicroseconds - episodeStartMicroseconds >=
			minimumRateIntervalMicroseconds) {
		    this->updateRate(static_cast<long double>(
			nowMicroseconds - episodeStartMicroseconds) /
			static_cast<long double>(completed));
		}
		return true;
	    }

	    if (completed == this->completedUnits)
		return false;
	    this->completedUnits = completed;
	    const int64_t duration = nowMicroseconds - this->sampleMicroseconds;
	    const uint64_t units = completed - this->sampleCompleted;
	    if (duration >= minimumRateIntervalMicroseconds && units > 0) {
		this->updateRate(static_cast<long double>(duration) /
		    static_cast<long double>(units));
		this->sampleCompleted = completed;
		this->sampleMicroseconds = nowMicroseconds;
	    }
	    return true;
	}

	long double microsecondsPerUnit(void) const
	{
	    return this->rateAvailable ? this->rateMicrosecondsPerUnit : 0.0L;
	}

	void resetObservation(void)
	{
	    this->observing = false;
	    this->totalUnits = 0;
	    this->completedUnits = 0;
	    this->sampleCompleted = 0;
	    this->sampleMicroseconds = 0;
	}

    private:
	void updateRate(long double rate)
	{
	    if (!(rate > 0.0L))
		return;
	    static constexpr long double recentSampleWeight = 0.25L;
	    this->rateMicrosecondsPerUnit = this->rateAvailable ?
		(1.0L - recentSampleWeight) * this->rateMicrosecondsPerUnit +
		    recentSampleWeight * rate : rate;
	    this->rateAvailable = true;
	}

	static constexpr int64_t minimumRateIntervalMicroseconds = 1000;
	bool observing = false;
	bool rateAvailable = false;
	uint64_t totalUnits = 0;
	uint64_t completedUnits = 0;
	uint64_t sampleCompleted = 0;
	int64_t sampleMicroseconds = 0;
	long double rateMicrosecondsPerUnit = 0.0L;
    };

    void beginEpisode(const Inputs &inputs)
    {
	this->episodeRevision = inputs.episodeRevision;
	this->refinementTierEpoch = inputs.refinementTierEpoch;
	this->episodeStartMicroseconds =
	    inputs.episodeStartMicroseconds > 0 &&
	    inputs.episodeStartMicroseconds <= inputs.observationMicroseconds ?
	    inputs.episodeStartMicroseconds : inputs.observationMicroseconds;
	this->unknownForegroundWork = inputs.unknownForegroundWork;
	this->finalPresentationMicroseconds =
	    inputs.finalPresentationMicroseconds;
	this->lastTerminal = false;
	this->fractionFloor = 0.0f;
	this->lastEstimate = Estimate();
	this->cycleForecast.reset();
	for (Rate &rate : this->rates)
	    rate.resetObservation();
    }

    uint64_t episodeRevision = 0;
    BObolLodPolicyEpoch refinementTierEpoch;
    int64_t episodeStartMicroseconds = 0;
    bool unknownForegroundWork = false;
    bool lastTerminal = false;
    uint64_t finalPresentationMicroseconds = 0;
    float fractionFloor = 0.0f;
    Estimate lastEstimate;
    CycleForecast cycleForecast;
    std::array<Rate, static_cast<size_t>(Rank::COUNT)> rates;
    static constexpr long double microsecondsPerMillisecond = 1000.0L;
};

#endif /* LIBBOBOL_LOD_PROGRESS_ESTIMATOR_PRIVATE_H */
