/*       T E S T _ L O D _ P R O G R E S S _ E S T I M A T O R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "lod_progress_estimator_private.h"
#include "bu/app.h"

#include <cmath>
#include <cstdio>
#include <type_traits>

static_assert(std::is_trivially_copyable<BObolLodProgressEstimator>::value,
    "progress estimator must remain an allocation-free value");

int
main(int argc, char **argv)
{
    (void)argc;
    bu_setprogname(argv[0]);
    using Estimator = BObolLodProgressEstimator;

    Estimator finalizationEstimator;
    Estimator::Inputs finalizationInput;
    finalizationInput.episodeRevision = 1;
    finalizationInput.refinementTierEpoch.set(1);
    finalizationInput.observationMicroseconds = 1000000;
    finalizationInput.finalPresentationMicroseconds = 20000;
    const Estimator::Estimate finalizationEstimate =
	finalizationEstimator.evaluate(finalizationInput);
    if (!finalizationEstimate.available ||
	finalizationEstimate.fraction > 1.0e-6f ||
	finalizationEstimate.remainingMilliseconds != 20) {
	std::fprintf(stderr,
	    "first finite finalization estimate was unavailable\n");
	return 1;
    }

    Estimator estimator;
    Estimator::Inputs input;
    input.episodeRevision = 1;
    input.refinementTierEpoch.set(1);
    input.observationMicroseconds = 1000000;
    input.episodeStartMicroseconds = input.observationMicroseconds;
    input.finalPresentationMicroseconds = 100000;
    Estimator::WorkRank &discovery = input.rank(Estimator::Rank::DISCOVERY);
    discovery.present = true;
    discovery.total = 100;

    Estimator::Estimate estimate = estimator.evaluate(input);
    if (estimate.available) {
	std::fprintf(stderr,
	    "progress estimate invented an unobserved discovery rate\n");
	return 1;
    }

    input.observationMicroseconds = 2000000;
    discovery.completed = 50;
    estimate = estimator.evaluate(input);
    if (!estimate.available || estimate.remainingMilliseconds != 1100 ||
	estimate.fraction < 0.47f || estimate.fraction > 0.48f) {
	std::fprintf(stderr,
	    "progress estimate did not use the observed finite rank\n");
	return 1;
    }
    const Estimator::Estimate unchanged = estimator.evaluate(input);
    if (!unchanged.available ||
	std::fabs(unchanged.fraction - estimate.fraction) > 1.0e-6f ||
	unchanged.remainingMilliseconds != estimate.remainingMilliseconds) {
	std::fprintf(stderr,
	    "observer time advanced progress without completed work\n");
	return 1;
    }

    input.observationMicroseconds = 3000000;
    discovery.completed = 75;
    const Estimator::Estimate advanced = estimator.evaluate(input);
    if (!advanced.available || advanced.fraction <= estimate.fraction ||
	advanced.remainingMilliseconds >= estimate.remainingMilliseconds) {
	std::fprintf(stderr,
	    "finite-rank progress did not improve the estimate\n");
	return 1;
    }

    input.unknownForegroundWork = true;
    if (estimator.evaluate(input).available) {
	std::fprintf(stderr,
	    "unknown foreground work produced a determinate estimate\n");
	return 1;
    }
    input.unknownForegroundWork = false;
    input.terminal = true;
    estimate = estimator.evaluate(input);
    if (!estimate.available || std::fabs(estimate.fraction - 1.0f) >
	1.0e-6f || estimate.remainingMilliseconds != 0) {
	std::fprintf(stderr, "terminal progress estimate was invalid\n");
	return 1;
    }

    input.terminal = false;
    estimate = estimator.evaluate(input);
    if (!estimate.available || estimate.fraction >= 1.0f ||
	estimate.remainingMilliseconds == 0) {
	std::fprintf(stderr,
	    "reopened progress estimate retained terminal readiness\n");
	return 1;
    }

    /* A compact stream may publish a newer inventory revision while the
     * source/view/policy episode remains unchanged.  A later transition time
     * must not reset either the measured rate or the monotonic fraction. */
    input.episodeStartMicroseconds = input.observationMicroseconds;
    const Estimator::Estimate appended = estimator.evaluate(input);
    if (!appended.available ||
	std::fabs(appended.fraction - advanced.fraction) > 1.0e-6f ||
	appended.remainingMilliseconds != advanced.remainingMilliseconds) {
	std::fprintf(stderr,
	    "append-only inventory restarted the progress episode\n");
	return 1;
    }

    /* Loading a different source in the same camera/policy tuple is a new
     * user-visible episode.  Append-only inventory publications deliberately
     * retain this identity and therefore cannot restart the progress clock. */
    input.episodeRevision = 2;
    input.observationMicroseconds = 3500000;
    input.episodeStartMicroseconds = input.observationMicroseconds;
    discovery.completed = 0;
    estimate = estimator.evaluate(input);
    if (!estimate.available || estimate.fraction > 1.0e-6f ||
	estimate.remainingMilliseconds == 0) {
	std::fprintf(stderr,
	    "source inventory did not reset the progress episode\n");
	return 1;
    }

    /* Visible-detail ETA is based on completed refinement cycles, not the
     * transient rate at which worker results enter the publication queue. */
    Estimator cycleEstimator;
    Estimator::Inputs cycleInput;
    cycleInput.episodeRevision = 10;
    cycleInput.refinementTierEpoch.set(40);
    cycleInput.episodeStartMicroseconds = 10000000;
    cycleInput.observationMicroseconds = 10000000;
    cycleInput.renderCompletionSerial = 100;
    cycleInput.finalPresentationMicroseconds = 50000;
    Estimator::WorkRank &visible = cycleInput.rank(
	Estimator::Rank::VISIBLE_RESOLUTION);
    visible.present = true;
    visible.total = 1000;
    visible.completed = 100;
    if (cycleEstimator.evaluate(cycleInput).available) {
	std::fprintf(stderr,
	    "refinement forecast qualified without completed cycles\n");
	return 1;
    }

    cycleInput.renderCompletionSerial = 101;
    cycleInput.observationMicroseconds = 10100000;
    visible.completed = 200;
    if (cycleEstimator.evaluate(cycleInput).available) {
	std::fprintf(stderr,
	    "refinement forecast qualified after only one cycle\n");
	return 1;
    }
    cycleInput.renderCompletionSerial = 102;
    cycleInput.observationMicroseconds = 10220000;
    visible.completed = 280;
    if (cycleEstimator.evaluate(cycleInput).available) {
	std::fprintf(stderr,
	    "refinement forecast qualified after only two cycles\n");
	return 1;
    }
    cycleInput.renderCompletionSerial = 103;
    cycleInput.observationMicroseconds = 10370000;
    visible.completed = 340;
    estimate = cycleEstimator.evaluate(cycleInput);
    if (!estimate.available || !estimate.refinementCycleBased ||
	estimate.remainingRefinementCycles != 11 ||
	estimate.remainingMilliseconds != 1272) {
	std::fprintf(stderr,
	    "refinement forecast did not preserve aggregate cycle throughput\n");
	return 1;
    }

    /* Discovering more unresolved work, including between frame boundaries,
     * immediately retracts the forecast. */
    cycleInput.observationMicroseconds = 10380000;
    visible.completed = 300;
    if (cycleEstimator.evaluate(cycleInput).available) {
	std::fprintf(stderr,
	    "rising refinement frontier retained a stale forecast\n");
	return 1;
    }

    cycleInput.renderCompletionSerial = 104;
    cycleInput.observationMicroseconds = 10530000;
    visible.completed = 360;
    (void)cycleEstimator.evaluate(cycleInput);
    cycleInput.renderCompletionSerial = 105;
    cycleInput.observationMicroseconds = 10680000;
    visible.completed = 420;
    (void)cycleEstimator.evaluate(cycleInput);
    cycleInput.renderCompletionSerial = 106;
    cycleInput.observationMicroseconds = 10830000;
    visible.completed = 470;
    if (!cycleEstimator.evaluate(cycleInput).available) {
	std::fprintf(stderr,
	    "improving refinement cycles did not restore the forecast\n");
	return 1;
    }
    cycleInput.renderCompletionSerial = 107;
    cycleInput.observationMicroseconds = 10980000;
    if (cycleEstimator.evaluate(cycleInput).available) {
	std::fprintf(stderr,
	    "flat refinement cycle retained a stale forecast\n");
	return 1;
    }

    /* An internal quality tier keeps the episode clock, but its newly opened
     * frontier must establish independent cycle evidence. */
    cycleInput.refinementTierEpoch.advance();
    cycleInput.renderCompletionSerial = 108;
    cycleInput.observationMicroseconds = 11080000;
    visible.completed = 520;
    if (cycleEstimator.evaluate(cycleInput).available) {
	std::fprintf(stderr,
	    "new refinement tier reused predecessor cycle evidence\n");
	return 1;
    }
    cycleInput.renderCompletionSerial = 109;
    cycleInput.observationMicroseconds = 11180000;
    visible.completed = 570;
    (void)cycleEstimator.evaluate(cycleInput);
    cycleInput.renderCompletionSerial = 110;
    cycleInput.observationMicroseconds = 11280000;
    visible.completed = 620;
    (void)cycleEstimator.evaluate(cycleInput);
    cycleInput.renderCompletionSerial = 111;
    cycleInput.observationMicroseconds = 11380000;
    visible.completed = 670;
    estimate = cycleEstimator.evaluate(cycleInput);
    if (!estimate.available || !estimate.refinementCycleBased ||
	estimate.remainingRefinementCycles != 9 ||
	estimate.remainingMilliseconds != 825 || estimate.fraction < 0.6f) {
	std::fprintf(stderr,
	    "new refinement tier did not preserve the enclosing episode clock\n");
	return 1;
    }

    visible.total = 1100;
    if (cycleEstimator.evaluate(cycleInput).available) {
	std::fprintf(stderr,
	    "changing refinement denominator retained a stale forecast\n");
	return 1;
    }
    return 0;
}
