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
    finalizationInput.viewEpoch.set(1);
    finalizationInput.policyEpoch.set(1);
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
    input.viewEpoch.set(1);
    input.policyEpoch.set(1);
    input.observationMicroseconds = 1000000;
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

    input.viewEpoch.set(2);
    input.observationMicroseconds = 4000000;
    discovery.completed = 0;
    estimate = estimator.evaluate(input);
    if (!estimate.available || estimate.fraction > 1.0e-6f ||
	estimate.remainingMilliseconds == 0) {
	std::fprintf(stderr,
	    "learned rate was not reusable in a new view epoch\n");
	return 1;
    }
    return 0;
}
