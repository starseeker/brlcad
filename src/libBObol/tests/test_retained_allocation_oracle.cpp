/*     T E S T _ R E T A I N E D _ A L L O C A T I O N _ O R A C L E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BViewLod.h"
#include "retained_allocation_private.h"

#include <cmath>
#include <cstdio>
#include <vector>

static bool
oracle_lower_priority(const BObolRetainedMarginalUpgrade &a,
    const BObolRetainedMarginalUpgrade &b)
{
    if (a.qualityFloorViolation != b.qualityFloorViolation)
        return !a.qualityFloorViolation;
    if (a.qualityFloorViolation &&
        (a.weightedError < b.weightedError ||
         a.weightedError > b.weightedError))
        return a.weightedError < b.weightedError;
    if (a.valuePerCost < b.valuePerCost || a.valuePerCost > b.valuePerCost)
        return a.valuePerCost < b.valuePerCost;
    if (a.weightedError < b.weightedError ||
        a.weightedError > b.weightedError)
        return a.weightedError < b.weightedError;
    return a.candidateIndex > b.candidateIndex;
}

static size_t
oracle_winner(const std::vector<BObolRetainedMarginalUpgrade> &upgrades)
{
    size_t winner = 0;
    for (size_t i = 1; i < upgrades.size(); ++i) {
	if (oracle_lower_priority(upgrades[winner], upgrades[i]))
	    winner = i;
    }
    return winner;
}

int
main()
{
    struct ErrorContractCase {
	double target;
	double projectedLimit;
    };
    const ErrorContractCase errorContractCases[] = {
	{0.25, 0.75},
	{0.5, 1.5},
	{1.0, 3.0},
	{2.0, 6.0}
    };
    for (const ErrorContractCase &testCase : errorContractCases) {
	if (std::fabs(BObolViewLodState::prominentMaximumProjectedError(
		testCase.target) - testCase.projectedLimit) > 1.0e-12) {
	    fprintf(stderr,
		"fractional-pixel prominent error contract mismatch\n");
	    return 1;
	}
    }

    const double errors[] = {0.25, 1.0, 4.0};
    const double values[] = {0.25, 1.0, 4.0};
    std::vector<BObolRetainedMarginalUpgrade> domain;
    for (int floor = 0; floor < 2; ++floor) {
	for (double error : errors) {
	    for (double value : values) {
		for (size_t index = 0; index < 3; ++index) {
		    BObolRetainedMarginalUpgrade upgrade;
		    upgrade.qualityFloorViolation = floor ? true : false;
		    upgrade.weightedError = error;
		    upgrade.valuePerCost = value;
		    upgrade.candidateIndex = index;
		    domain.push_back(upgrade);
		}
	    }
	}
    }

    for (size_t a = 0; a < domain.size(); ++a) {
	for (size_t b = 0; b < domain.size(); ++b) {
	    if (bobol_retained_marginal_lower_priority(domain[a], domain[b]) !=
		oracle_lower_priority(domain[a], domain[b])) {
		fprintf(stderr, "pairwise marginal ordering mismatch\n");
		return 1;
	    }
	}
    }

    /* Exhaustively check every bounded three-option frontier.  This covers
     * competing floors, equal-ratio minimax ties, and the deterministic index
     * tie break used to prevent a stable view from oscillating. */
    for (size_t a = 0; a < domain.size(); ++a) {
	for (size_t b = 0; b < domain.size(); ++b) {
	    for (size_t c = 0; c < domain.size(); ++c) {
		const std::vector<BObolRetainedMarginalUpgrade> frontier = {
		    domain[a], domain[b], domain[c]
		};
		const size_t expected = oracle_winner(frontier);
		size_t actual = 0;
		for (size_t i = 1; i < frontier.size(); ++i) {
		    if (bobol_retained_marginal_lower_priority(
			frontier[actual], frontier[i]))
		actual = i;
		}
		if (actual != expected) {
		    fprintf(stderr, "frontier marginal ordering mismatch\n");
		    return 1;
		}
	    }
	}
    }

    return 0;
}
