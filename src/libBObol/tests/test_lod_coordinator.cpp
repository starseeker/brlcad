/*              T E S T _ L O D _ C O O R D I N A T O R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "lod_coordinator_private.h"

#include <cstdio>
#include <array>
#include <type_traits>

static_assert(sizeof(BObolLodViewEpoch) == sizeof(uint64_t),
    "typed LoD epochs must remain zero-overhead");
static_assert(sizeof(BObolLodPolicyEpoch) == sizeof(uint64_t),
    "typed LoD epochs must remain zero-overhead");
static_assert(std::is_trivially_copyable<BObolLodViewEpoch>::value,
    "view epochs must remain trivially copyable");
static_assert(std::is_trivially_copyable<BObolLodPolicyEpoch>::value,
    "policy epochs must remain trivially copyable");
static_assert(!std::is_same<BObolLodViewEpoch, BObolLodPolicyEpoch>::value,
    "view and policy epochs must be distinct types");
static_assert(!std::is_convertible<BObolLodViewEpoch, uint64_t>::value,
    "epoch serialization must be explicit");

static int
test_epochs(void)
{
    BObolLodViewEpoch view(1);
    BObolLodPolicyEpoch policy(9);

    view.advance();
    if (view.value() != 2 || policy.value() != 9) {
	std::fprintf(stderr, "FAIL: independent epoch advance\n");
	return 1;
    }

    view.reset();
    policy.set(17);
    if (view.value() != 0 || policy.value() != 17) {
	std::fprintf(stderr, "FAIL: epoch reset/set\n");
	return 1;
    }

    return 0;
}

static int
test_phase_witnesses(void)
{
    BObolLodPhaseTracker tracker;
    using Phase = BObolLodPhaseTracker::Phase;
    using Witness = BObolLodPhaseTracker::Witness;

    if (tracker.currentPhase() != Phase::FALLBACK ||
	tracker.phaseTransitionSerial() != 0) {
	std::fprintf(stderr, "FAIL: initial phase\n");
	return 1;
    }

    Witness fallback;
    fallback.viewEpoch = 1;
    fallback.policyEpoch = 1;
    tracker.enterFallback(fallback);
    const uint64_t fallbackSequence =
	tracker.phaseWitness(Phase::FALLBACK).sequence;
    if (!fallbackSequence || tracker.phaseTransitionSerial() != 0) {
	std::fprintf(stderr, "FAIL: initial fallback witness\n");
	return 1;
    }

    /* Reobserving identical state is not false progress. */
    tracker.enterFallback(fallback);
    if (tracker.phaseWitness(Phase::FALLBACK).sequence !=
	    fallbackSequence) {
	std::fprintf(stderr, "FAIL: duplicate witness advanced\n");
	return 1;
    }

    Witness coverage = fallback;
    coverage.activeGeneration = 41;
    coverage.visibleCount = 100;
    coverage.completedCount = 12;
    coverage.pendingCount = 88;
    tracker.enterCoverage(coverage);
    if (tracker.currentPhase() != Phase::COVERAGE ||
	tracker.phaseTransitionSerial() != 1 ||
	!tracker.phaseWitness(Phase::COVERAGE).sequence) {
	std::fprintf(stderr, "FAIL: coverage transition\n");
	return 1;
    }

    const uint64_t coverageSequence =
	tracker.phaseWitness(Phase::COVERAGE).sequence;
    coverage.completedCount = 37;
    coverage.pendingCount = 63;
    tracker.enterCoverage(coverage);
    if (tracker.phaseTransitionSerial() != 1 ||
	tracker.phaseWitness(Phase::COVERAGE).sequence <=
	    coverageSequence) {
	std::fprintf(stderr, "FAIL: coverage progress witness\n");
	return 1;
    }

    Witness interactive = coverage;
    interactive.viewEpoch = 2;
    interactive.renderSerial = 7;
    tracker.enterInteractive(interactive);
    Witness settling = interactive;
    settling.policyEpoch = 2;
    settling.pendingCount = 11;
    tracker.enterSettling(settling);
    Witness stable = settling;
    stable.pendingCount = 0;
    stable.completedCount = stable.visibleCount;
    tracker.enterStable(stable);
    Witness compacting = stable;
    compacting.residentDemandRevision = 83;
    compacting.pendingCount = 4;
    tracker.enterCompacting(compacting);

    if (tracker.currentPhase() != Phase::COMPACTING ||
	tracker.phaseTransitionSerial() != 5 ||
	tracker.phaseWitness(Phase::INTERACTIVE).viewEpoch != 2 ||
	tracker.phaseWitness(Phase::SETTLING).pendingCount != 11 ||
	tracker.phaseWitness(Phase::STABLE).pendingCount != 0 ||
	tracker.phaseWitness(Phase::COMPACTING).
	    residentDemandRevision != 83) {
	std::fprintf(stderr, "FAIL: retained per-phase witnesses\n");
	return 1;
    }

    tracker.enterFallback(fallback);
    if (tracker.currentPhase() != Phase::FALLBACK ||
	tracker.phaseTransitionSerial() != 6) {
	std::fprintf(stderr, "FAIL: reset to fallback\n");
	return 1;
    }

    return 0;
}

static int
test_phase_decision_contract(void)
{
    using Tracker = BObolLodPhaseTracker;
    using Phase = Tracker::Phase;

    struct Scenario {
	const char *name;
	Tracker::Inputs inputs;
	Phase expected;
    };

    const std::array<Scenario, 10> scenarios = {{
	{"cold-no-generation", {false, false, true, false, false, false},
	    Phase::FALLBACK},
	{"cold-coverage", {false, false, true, false, true, false},
	    Phase::COVERAGE},
	{"camera-gesture", {true, false, true, false, true, true},
	    Phase::INTERACTIVE},
	{"quiet-refinement", {false, false, false, true, true, true},
	    Phase::SETTLING},
	{"stable", {false, false, false, true, true, false},
	    Phase::STABLE},
	{"resident-compaction", {false, true, false, true, true, false},
	    Phase::COMPACTING},
	/* Interaction always wins: memory reclamation cannot block input. */
	{"interaction-during-compaction",
	    {true, true, false, true, true, true}, Phase::INTERACTIVE},
	/* Coverage always precedes richer refinement. */
	{"coverage-before-settling",
	    {false, false, true, false, true, true}, Phase::COVERAGE},
	/* An invalidated inventory without a generation returns to fallback. */
	{"erased-generation", {false, false, false, false, false, true},
	    Phase::FALLBACK},
	/* Selection/style-only mutations do not disturb a stable LoD proof. */
	{"selection-only", {false, false, false, true, true, false},
	    Phase::STABLE}
    }};

    for (const Scenario &scenario : scenarios) {
	if (Tracker::phaseFor(scenario.inputs) == scenario.expected)
	    continue;
	std::fprintf(stderr, "FAIL: phase decision %s\n", scenario.name);
	return 1;
    }
    return 0;
}

/*
 * Model the admission edge independently of worker timing.  Only results from
 * the active generation and exact view/policy epochs may advance a witness;
 * reordered stale completion, cancellation, and a source erase are no-ops
 * until a new generation is explicitly begun.
 */
static int
test_adversarial_result_order(void)
{
    using Phase = BObolLodPhaseTracker::Phase;
    BObolLodPhaseTracker tracker;
    uint64_t generation = 7;
    uint64_t view = 3;
    uint64_t policy = 5;
    size_t completed = 0;

    auto publish = [&](uint64_t resultGeneration, uint64_t resultView,
		       uint64_t resultPolicy) {
	if (resultGeneration != generation || resultView != view ||
	    resultPolicy != policy)
	    return false;
	BObolLodPhaseTracker::Witness witness;
	witness.activeGeneration = generation;
	witness.viewEpoch = view;
	witness.policyEpoch = policy;
	witness.visibleCount = 3;
	witness.completedCount = ++completed;
	witness.pendingCount = 3 - completed;
	if (completed == 3)
	    tracker.enterSettling(witness);
	else
	    tracker.enterCoverage(witness);
	return true;
    };

    if (publish(6, view, policy) ||
	publish(generation, view - 1, policy) ||
	publish(generation, view, policy - 1) ||
	tracker.phaseTransitionSerial() != 0) {
	std::fprintf(stderr, "FAIL: stale result advanced phase witness\n");
	return 1;
    }
    if (!publish(generation, view, policy) ||
	tracker.currentPhase() != Phase::COVERAGE) {
	std::fprintf(stderr, "FAIL: active result did not enter coverage\n");
	return 1;
    }

    /* A camera change invalidates completions already in flight. */
    ++view;
    if (publish(generation, view - 1, policy) || completed != 1) {
	std::fprintf(stderr, "FAIL: reordered camera result was admitted\n");
	return 1;
    }
    if (!publish(generation, view, policy) ||
	!publish(generation, view, policy) ||
	tracker.currentPhase() != Phase::SETTLING) {
	std::fprintf(stderr, "FAIL: current results did not settle\n");
	return 1;
    }

    /* Cancellation/erase retires the generation before any late callback. */
    const uint64_t cancelledGeneration = generation;
    generation = 0;
    BObolLodPhaseTracker::Witness fallback;
    fallback.viewEpoch = view;
    fallback.policyEpoch = policy;
    tracker.enterFallback(fallback);
    const uint64_t fallbackSequence =
	tracker.phaseWitness(Phase::FALLBACK).sequence;
    if (publish(cancelledGeneration, view, policy) ||
	tracker.phaseWitness(Phase::FALLBACK).sequence != fallbackSequence) {
	std::fprintf(stderr, "FAIL: cancelled result changed fallback proof\n");
	return 1;
    }

    return 0;
}

int
main(void)
{
    if (test_epochs())
	return 1;
    if (test_phase_witnesses())
	return 1;
    if (test_phase_decision_contract())
	return 1;
    if (test_adversarial_result_order())
	return 1;
    return 0;
}
