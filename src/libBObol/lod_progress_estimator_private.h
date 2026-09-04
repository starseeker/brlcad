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
 * A rank without a denominator or an observed rate makes the result
 * indeterminate.  Keeping that distinction explicit avoids a confident 99%
 * display while one expensive source or renderer-preparation item remains.
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
	BObolLodViewEpoch viewEpoch;
	BObolLodPolicyEpoch policyEpoch;
	int64_t observationMicroseconds = 0;
	bool terminal = false;
	bool unknownForegroundWork = false;
	uint64_t finalPresentationMicroseconds = 0;
	std::array<WorkRank, static_cast<size_t>(Rank::COUNT)> ranks;

	WorkRank &rank(Rank value)
	{
	    return this->ranks[static_cast<size_t>(value)];
	}
    };

    struct Estimate {
	bool available = false;
	float fraction = 0.0f;
	uint64_t remainingMilliseconds = 0;
    };

    Estimate evaluate(const Inputs &inputs)
    {
	const bool newEpisode = this->viewEpoch != inputs.viewEpoch ||
	    this->policyEpoch != inputs.policyEpoch;
	if (newEpisode)
	    this->beginEpisode(inputs);
	const bool terminalChanged = inputs.terminal != this->lastTerminal;
	this->lastTerminal = inputs.terminal;

	if (inputs.terminal) {
	    this->lastEstimate.available = true;
	    this->lastEstimate.fraction = 1.0f;
	    this->lastEstimate.remainingMilliseconds = 0;
	    return this->lastEstimate;
	}

	bool rankChanged = false;
	for (size_t i = 0; i < this->rates.size(); ++i)
	    rankChanged = this->rates[i].observe(inputs.ranks[i],
		inputs.observationMicroseconds, this->episodeStartMicroseconds) ||
		rankChanged;

	const bool estimateInputsChanged = newEpisode || terminalChanged ||
	    rankChanged ||
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
	for (size_t i = 0; i < this->rates.size(); ++i) {
	    const WorkRank &rank = inputs.ranks[i];
	    if (!rank.present || rank.total == 0)
		continue;
	    const uint64_t completed = std::min(rank.completed, rank.total);
	    if (completed == rank.total)
		continue;
	    hasIncompleteRank = true;
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

	if (inputs.finalPresentationMicroseconds > 0)
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
	this->lastEstimate.fraction = fraction;
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
	this->viewEpoch.reset();
	this->policyEpoch.reset();
	this->episodeStartMicroseconds = 0;
	this->unknownForegroundWork = false;
	this->finalPresentationMicroseconds = 0;
	this->lastTerminal = false;
	this->fractionFloor = 0.0f;
	this->lastEstimate = Estimate();
	for (Rate &rate : this->rates)
	    rate.resetObservation();
    }

private:
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
	this->viewEpoch = inputs.viewEpoch;
	this->policyEpoch = inputs.policyEpoch;
	this->episodeStartMicroseconds = inputs.observationMicroseconds;
	this->unknownForegroundWork = inputs.unknownForegroundWork;
	this->finalPresentationMicroseconds =
	    inputs.finalPresentationMicroseconds;
	this->lastTerminal = false;
	this->fractionFloor = 0.0f;
	this->lastEstimate = Estimate();
	for (Rate &rate : this->rates)
	    rate.resetObservation();
    }

    BObolLodViewEpoch viewEpoch;
    BObolLodPolicyEpoch policyEpoch;
    int64_t episodeStartMicroseconds = 0;
    bool unknownForegroundWork = false;
    bool lastTerminal = false;
    uint64_t finalPresentationMicroseconds = 0;
    float fractionFloor = 0.0f;
    Estimate lastEstimate;
    std::array<Rate, static_cast<size_t>(Rank::COUNT)> rates;
    static constexpr long double microsecondsPerMillisecond = 1000.0L;
};

#endif /* LIBBOBOL_LOD_PROGRESS_ESTIMATOR_PRIVATE_H */
