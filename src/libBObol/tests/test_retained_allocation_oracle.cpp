/*     T E S T _ R E T A I N E D _ A L L O C A T I O N _ O R A C L E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BViewLod.h"
#include "lod_visual_importance_private.h"
#include "retained_allocation_private.h"

#include <cmath>
#include <cstdio>
#include <vector>

static bool
oracle_lower_priority(const BObolRetainedMarginalUpgrade &a,
    const BObolRetainedMarginalUpgrade &b)
{
    if (a.visualEmphasis != b.visualEmphasis)
	return a.visualEmphasis < b.visualEmphasis;
    if (a.qualityFloorViolation != b.qualityFloorViolation)
        return !a.qualityFloorViolation;
    if (a.qualityFloorViolation &&
        (a.normalizedError < b.normalizedError ||
         a.normalizedError > b.normalizedError))
        return a.normalizedError < b.normalizedError;
    if (a.visualBenefitPerCost < b.visualBenefitPerCost ||
	a.visualBenefitPerCost > b.visualBenefitPerCost)
        return a.visualBenefitPerCost < b.visualBenefitPerCost;
    if (a.normalizedError < b.normalizedError ||
        a.normalizedError > b.normalizedError)
        return a.normalizedError < b.normalizedError;
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

    const double thinFootprint = bobol_lod_visual_footprint(
	0.0, 64.0, 4.0);
    if (std::fabs(thinFootprint - 16.0) > 1.0e-12 ||
	!bobol_lod_visual_prominent(thinFootprint) ||
	bobol_lod_visual_prominent(thinFootprint - 0.01)) {
	fprintf(stderr, "projected visual-footprint contract mismatch\n");
	return 1;
    }
    const double ordinaryProtected = bobol_lod_protected_visual_error(
	BOBOL_LOD_VISUAL_ORDINARY, thinFootprint, 0.25);
    const double smallProtected = bobol_lod_protected_visual_error(
	BOBOL_LOD_VISUAL_ORDINARY, thinFootprint - 0.01, 0.25);
    if (std::fabs(ordinaryProtected - 0.75) > 1.0e-12 ||
	!std::isinf(smallProtected)) {
	fprintf(stderr, "recognizable-feature floor contract mismatch\n");
	return 1;
    }
    const double smallBenefit = bobol_lod_marginal_visual_benefit(
	16.0, 2.0, 1.0, 0.25);
    const double largeBenefit = bobol_lod_marginal_visual_benefit(
	64.0, 2.0, 1.0, 0.25);
    if (std::fabs(smallBenefit - 64.0) > 1.0e-12 ||
	std::fabs(largeBenefit - 256.0) > 1.0e-12) {
	fprintf(stderr, "marginal visual-benefit contract mismatch\n");
	return 1;
    }

    const double errors[] = {0.25, 4.0};
    const double values[] = {0.25, 4.0};
    std::vector<BObolRetainedMarginalUpgrade> domain;
    for (unsigned int emphasis = BOBOL_LOD_VISUAL_ORDINARY;
	emphasis <= BOBOL_LOD_VISUAL_SELECTED; ++emphasis) {
	for (int floor = 0; floor < 2; ++floor) {
	    for (double error : errors) {
		for (double value : values) {
		    for (size_t index = 0; index < 2; ++index) {
		    BObolRetainedMarginalUpgrade upgrade;
		    upgrade.visualEmphasis = emphasis;
		    upgrade.qualityFloorViolation = floor ? true : false;
		    upgrade.normalizedError = error;
		    upgrade.visualBenefitPerCost = value;
		    upgrade.candidateIndex = index;
		    domain.push_back(upgrade);
		    }
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
