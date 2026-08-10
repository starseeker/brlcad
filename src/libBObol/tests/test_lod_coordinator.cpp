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
#include <cstdint>
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
static_assert(std::is_trivially_copyable<BObolLodResourcePolicy>::value,
    "resource policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodHeadroomPolicy>::value,
    "headroom policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodCoveragePolicy>::value,
    "coverage policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodConvergencePolicy>::value,
    "convergence policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodBudgetPolicy>::value,
    "budget policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodViewDemandPolicy>::value,
    "view-demand policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodPresentationPolicy>::value,
    "presentation policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodQualityPolicy>::value,
    "quality policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodPublicationPolicy>::value,
    "publication policy must remain an allocation-free value");
static_assert(std::is_trivially_copyable<BObolLodCompactionPolicy>::value,
    "compaction policy must remain an allocation-free value");

using LodMachine = BObolLodStateMachine;

static LodMachine::Inputs
inputs_for_phase(LodMachine::Phase phase)
{
    LodMachine::Inputs inputs;
    switch (phase) {
	case LodMachine::Phase::FALLBACK:
	    inputs.coverageActive = true;
	    break;
	case LodMachine::Phase::COVERAGE:
	    inputs.coverageActive = true;
	    inputs.generationActive = true;
	    break;
	case LodMachine::Phase::INTERACTIVE:
	    inputs.interactive = true;
	    inputs.gestureActive = true;
	    break;
	case LodMachine::Phase::SETTLING:
	    inputs.coverageComplete = true;
	    inputs.generationActive = true;
	    inputs.settlingWork = true;
	    break;
	case LodMachine::Phase::STABLE:
	    inputs.coverageComplete = true;
	    inputs.generationActive = true;
	    break;
	case LodMachine::Phase::COMPACTING:
	    inputs.coverageComplete = true;
	    inputs.compacting = true;
	    break;
	case LodMachine::Phase::COUNT:
	    break;
    }
    return inputs;
}

static LodMachine::Transition
dispatch_phase(LodMachine &machine, LodMachine::Event event,
    LodMachine::Phase phase, const LodMachine::Witness &witness)
{
    return machine.dispatch(event, inputs_for_phase(phase), witness);
}

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
    LodMachine tracker;
    using Phase = LodMachine::Phase;
    using Witness = LodMachine::Witness;

    if (tracker.currentPhase() != Phase::FALLBACK ||
	tracker.phaseTransitionSerial() != 0) {
	std::fprintf(stderr, "FAIL: initial phase\n");
	return 1;
    }

    Witness fallback;
    fallback.viewEpoch = 1;
    fallback.policyEpoch = 1;
    dispatch_phase(tracker, LodMachine::Event::INITIALIZE,
	Phase::FALLBACK, fallback);
    const uint64_t fallbackSequence =
	tracker.phaseWitness(Phase::FALLBACK).sequence;
    if (!fallbackSequence || tracker.phaseTransitionSerial() != 0) {
	std::fprintf(stderr, "FAIL: initial fallback witness\n");
	return 1;
    }

    /* Reobserving identical state is not false progress. */
    dispatch_phase(tracker, LodMachine::Event::WORK_PUMPED,
	Phase::FALLBACK, fallback);
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
    dispatch_phase(tracker, LodMachine::Event::WORK_SCHEDULED,
	Phase::COVERAGE, coverage);
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
    dispatch_phase(tracker, LodMachine::Event::RESULT_PUBLISHED,
	Phase::COVERAGE, coverage);
    if (tracker.phaseTransitionSerial() != 1 ||
	tracker.phaseWitness(Phase::COVERAGE).sequence <=
	    coverageSequence) {
	std::fprintf(stderr, "FAIL: coverage progress witness\n");
	return 1;
    }

    Witness interactive = coverage;
    interactive.viewEpoch = 2;
    interactive.renderSerial = 7;
    dispatch_phase(tracker, LodMachine::Event::INTERACTION_STARTED,
	Phase::INTERACTIVE, interactive);
    Witness settling = interactive;
    settling.policyEpoch = 2;
    settling.pendingCount = 11;
    dispatch_phase(tracker, LodMachine::Event::INTERACTION_ENDED,
	Phase::SETTLING, settling);
    Witness stable = settling;
    stable.pendingCount = 0;
    stable.completedCount = stable.visibleCount;
    dispatch_phase(tracker, LodMachine::Event::FRAME_COMPLETED,
	Phase::STABLE, stable);
    Witness compacting = stable;
    compacting.residentDemandRevision = 83;
    compacting.pendingCount = 4;
    dispatch_phase(tracker, LodMachine::Event::WORK_PUMPED,
	Phase::COMPACTING, compacting);

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

    dispatch_phase(tracker, LodMachine::Event::VIEW_INVALIDATED,
	Phase::FALLBACK, fallback);
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
    using Tracker = BObolLodStateMachine;
    using Phase = Tracker::Phase;

    struct Scenario {
	const char *name;
	Tracker::Inputs inputs;
	Phase expected;
    };

    const std::array<Scenario, 14> scenarios = {{
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
	/* A stale compaction latch cannot bypass minimum coverage. */
	{"coverage-preempts-compaction",
	    {false, true, true, false, true, false}, Phase::COVERAGE},
	/* Coverage always precedes richer refinement. */
	{"coverage-before-settling",
	    {false, false, true, false, true, true}, Phase::COVERAGE},
	/* An invalidated inventory without a generation returns to fallback. */
	{"erased-generation", {false, false, false, false, false, true},
	    Phase::FALLBACK},
	/* Selection/style-only mutations do not disturb a stable LoD proof. */
	{"selection-only", {false, false, false, true, true, false},
	    Phase::STABLE},
	/* Persistent pressure after one recovery pass is a valid, explicitly
	 * memory-limited stable terminal state. */
	{"cpu-pressure-observed",
	    {false, false, false, true, true, false, false, true, false,
		false}, Phase::STABLE},
	{"cpu-pressure-recovery",
	    {false, false, false, true, true, false, false, true, false,
		true}, Phase::COMPACTING},
	/* Input remains preemptive while a safe pressure recovery is pending. */
	{"gpu-pressure-interaction",
	    {true, false, false, true, true, false, true, false, true,
		true}, Phase::INTERACTIVE}
    }};

    for (const Scenario &scenario : scenarios) {
	Tracker tracker;
	Tracker::Witness witness;
	witness.visibleCount = 1;
	witness.completedCount = scenario.inputs.coverageComplete ? 1 : 0;
	tracker.dispatch(Tracker::Event::WORK_PUMPED, scenario.inputs, witness);
	if (tracker.currentPhase() == scenario.expected)
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
    using Phase = LodMachine::Phase;
    LodMachine tracker;
    uint64_t generation = 7;
    uint64_t view = 3;
    uint64_t policy = 5;
    size_t completed = 0;

    auto publish = [&](uint64_t resultGeneration, uint64_t resultView,
		       uint64_t resultPolicy) {
	if (resultGeneration != generation || resultView != view ||
	    resultPolicy != policy)
	    return false;
	LodMachine::Witness witness;
	witness.activeGeneration = generation;
	witness.viewEpoch = view;
	witness.policyEpoch = policy;
	witness.visibleCount = 3;
	witness.completedCount = ++completed;
	witness.pendingCount = 3 - completed;
	if (completed == 3)
	    dispatch_phase(tracker, LodMachine::Event::RESULT_PUBLISHED,
		Phase::SETTLING, witness);
	else
	    dispatch_phase(tracker, LodMachine::Event::RESULT_PUBLISHED,
		Phase::COVERAGE, witness);
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
    LodMachine::Witness fallback;
    fallback.viewEpoch = view;
    fallback.policyEpoch = policy;
    dispatch_phase(tracker, LodMachine::Event::GENERATION_CANCELLED,
	Phase::FALLBACK, fallback);
    const uint64_t fallbackSequence =
	tracker.phaseWitness(Phase::FALLBACK).sequence;
    if (publish(cancelledGeneration, view, policy) ||
	tracker.phaseWitness(Phase::FALLBACK).sequence != fallbackSequence) {
	std::fprintf(stderr, "FAIL: cancelled result changed fallback proof\n");
	return 1;
    }

    return 0;
}

static bool
same_machine_inputs(const LodMachine::Inputs &lhs,
    const LodMachine::Inputs &rhs)
{
    return lhs.interactive == rhs.interactive &&
	lhs.compacting == rhs.compacting &&
	lhs.coverageActive == rhs.coverageActive &&
	lhs.coverageComplete == rhs.coverageComplete &&
	lhs.generationActive == rhs.generationActive &&
	lhs.settlingWork == rhs.settlingWork &&
	lhs.gestureActive == rhs.gestureActive &&
	lhs.cpuMemoryPressure == rhs.cpuMemoryPressure &&
	lhs.gpuMemoryPressure == rhs.gpuMemoryPressure &&
	lhs.resourceRecoveryPending == rhs.resourceRecoveryPending;
}

static LodMachine::Phase
phase_for_canonical_inputs(const LodMachine::Inputs &inputs,
    const LodMachine::Witness &witness)
{
    if (inputs.interactive || inputs.gestureActive)
	return LodMachine::Phase::INTERACTIVE;
    if (inputs.coverageActive || !inputs.coverageComplete)
	return inputs.generationActive ? LodMachine::Phase::COVERAGE :
	    LodMachine::Phase::FALLBACK;
    if (inputs.compacting || inputs.resourceRecoveryPending)
	return LodMachine::Phase::COMPACTING;
    if (inputs.settlingWork || witness.pendingCount != 0)
	return LodMachine::Phase::SETTLING;
    return LodMachine::Phase::STABLE;
}

/*
 * Exhaust the complete scalar input domain for every lifecycle event.  This
 * is deliberately deterministic rather than a timing fuzzer: a failing mask
 * is an exact replay token, all contradictory caller observations are
 * covered, and the follow-up valid observation proves that canonicalization
 * never installs a sticky recovery latch.
 */
static int
test_exhaustive_state_contract(void)
{
    using Event = LodMachine::Event;
    using Invariant = LodMachine::Invariant;
    using Phase = LodMachine::Phase;
    static const size_t inputDomain = 1u << 10;
    static const size_t witnessDomain = 3;

    for (size_t inputMask = 0; inputMask < inputDomain; ++inputMask) {
	LodMachine::Inputs observed;
	observed.interactive = (inputMask & (1u << 0)) != 0;
	observed.compacting = (inputMask & (1u << 1)) != 0;
	observed.coverageActive = (inputMask & (1u << 2)) != 0;
	observed.coverageComplete = (inputMask & (1u << 3)) != 0;
	observed.generationActive = (inputMask & (1u << 4)) != 0;
	observed.settlingWork = (inputMask & (1u << 5)) != 0;
	observed.gestureActive = (inputMask & (1u << 6)) != 0;
	observed.cpuMemoryPressure = (inputMask & (1u << 7)) != 0;
	observed.gpuMemoryPressure = (inputMask & (1u << 8)) != 0;
	observed.resourceRecoveryPending =
	    (inputMask & (1u << 9)) != 0;

	for (size_t witnessCase = 0; witnessCase < witnessDomain;
		++witnessCase) {
	    LodMachine::Witness observedWitness;
	    observedWitness.viewEpoch = 11;
	    observedWitness.policyEpoch = 17;
	    observedWitness.activeGeneration = 23;
	    observedWitness.visibleCount = witnessCase ? 5 : 0;
	    observedWitness.completedCount = witnessCase == 1 ? 7 :
		(witnessCase == 2 ? 3 : 0);
	    observedWitness.pendingCount = witnessCase == 2 ? 2 : 0;

	    for (size_t eventIndex = 0;
		    eventIndex < static_cast<size_t>(Event::COUNT);
		    ++eventIndex) {
		const Event event = static_cast<Event>(eventIndex);
		LodMachine::Inputs canonical = observed;
		if (event == Event::INTERACTION_STARTED) {
		    canonical.interactive = true;
		    canonical.gestureActive = true;
		} else if (event == Event::INTERACTION_ENDED) {
		    canonical.gestureActive = false;
		} else if (event == Event::GENERATION_CANCELLED) {
		    canonical.generationActive = false;
		}
		if (!canonical.coverageComplete)
		    canonical.compacting = false;
		if (!canonical.cpuMemoryPressure &&
			!canonical.gpuMemoryPressure)
		    canonical.resourceRecoveryPending = false;
		if (observedWitness.pendingCount != 0)
		    canonical.settlingWork = true;

		LodMachine::Witness canonicalWitness = observedWitness;
		if (canonicalWitness.completedCount >
			canonicalWitness.visibleCount)
		    canonicalWitness.completedCount =
			canonicalWitness.visibleCount;

		uint32_t expectedInvariant = Invariant::INVARIANT_NONE;
		if (observedWitness.completedCount >
			observedWitness.visibleCount)
		    expectedInvariant |=
			Invariant::INVARIANT_COMPLETED_EXCEEDS_VISIBLE;
		if (observed.compacting && !observed.coverageComplete)
		    expectedInvariant |=
			Invariant::INVARIANT_COMPACTION_WITHOUT_COVERAGE;
		if (!observed.interactive && !observed.gestureActive &&
			!observed.compacting && observed.coverageComplete &&
			!observed.coverageActive && !observed.settlingWork &&
			observedWitness.pendingCount != 0)
		    expectedInvariant |=
			Invariant::INVARIANT_STABLE_WITH_PENDING_WORK;
		if (event == Event::INTERACTION_STARTED &&
			(!observed.interactive || !observed.gestureActive))
		    expectedInvariant |=
			Invariant::INVARIANT_INTERACTION_EVENT_WITHOUT_INTERACTION;
		if (event == Event::INTERACTION_ENDED &&
			observed.gestureActive)
		    expectedInvariant |=
			Invariant::INVARIANT_INTERACTION_END_WITH_GESTURE;
		if (event == Event::GENERATION_CANCELLED &&
			observed.generationActive)
		    expectedInvariant |=
			Invariant::INVARIANT_CANCEL_EVENT_WITH_GENERATION;
		if (observed.resourceRecoveryPending &&
			!observed.cpuMemoryPressure &&
			!observed.gpuMemoryPressure)
		    expectedInvariant |=
			Invariant::INVARIANT_RESOURCE_RECOVERY_WITHOUT_PRESSURE;

		const Phase expectedPhase = phase_for_canonical_inputs(
		    canonical, canonicalWitness);
		const bool expectedCanonicalization =
		    !same_machine_inputs(observed, canonical) ||
		    observedWitness.completedCount !=
			canonicalWitness.completedCount;

		LodMachine machine;
		const LodMachine::Transition transition = machine.dispatch(
		    event, observed, observedWitness);
		if (transition.current != expectedPhase ||
		    transition.invariantMask != expectedInvariant ||
		    transition.canonicalized != expectedCanonicalization ||
		    !same_machine_inputs(machine.lastInputs(), canonical) ||
		    !same_machine_inputs(
			machine.lastObservedInputs(), observed) ||
		    machine.phaseWitness(expectedPhase).completedCount !=
			canonicalWitness.completedCount ||
		    machine.invariantViolationCount() !=
			(expectedInvariant ? 1u : 0u)) {
		    std::fprintf(stderr,
			"FAIL: exhaustive state contract mask=%zu witness=%zu "
			"event=%zu phase=%u/%u invariant=0x%x/0x%x\n",
			inputMask, witnessCase, eventIndex,
			static_cast<unsigned int>(transition.current),
			static_cast<unsigned int>(expectedPhase),
			transition.invariantMask, expectedInvariant);
		    return 1;
		}

		/* No contradictory observation may poison the next valid quiet
		 * proof.  This is the bounded discharge guarantee used after
		 * cancellation, provider failure, and pressure recovery. */
		LodMachine::Inputs stable;
		stable.coverageComplete = true;
		stable.generationActive = true;
		LodMachine::Witness stableWitness;
		stableWitness.viewEpoch = 12;
		stableWitness.policyEpoch = 18;
		stableWitness.visibleCount = 5;
		stableWitness.completedCount = 5;
		const LodMachine::Transition recovered = machine.dispatch(
		    Event::VIEW_OBSERVED, stable, stableWitness);
		if (recovered.current != Phase::STABLE ||
		    recovered.invariantMask != Invariant::INVARIANT_NONE ||
		    machine.lastInvariantMask() != Invariant::INVARIANT_NONE ||
		    !machine.lastInputs().coverageComplete) {
		    std::fprintf(stderr,
			"FAIL: exhaustive state recovery mask=%zu witness=%zu "
			"event=%zu\n", inputMask, witnessCase, eventIndex);
		    return 1;
		}
	    }
	}
    }
    return 0;
}

static int
test_event_and_invariant_contract(void)
{
    using Event = LodMachine::Event;
    using Invariant = LodMachine::Invariant;
    using Phase = LodMachine::Phase;
    LodMachine machine;
    LodMachine::Witness witness;
    witness.viewEpoch = 4;
    witness.policyEpoch = 9;

    LodMachine::Transition initial = dispatch_phase(machine,
	Event::INITIALIZE, Phase::FALLBACK, witness);
    if (initial.changed || !initial.progressed ||
	machine.dispatchSerial() != 1 ||
	machine.lastDispatchedEvent() != Event::INITIALIZE ||
	machine.lastInvariantMask() != Invariant::INVARIANT_NONE) {
	std::fprintf(stderr, "FAIL: initial dispatch diagnostics\n");
	return 1;
    }

    LodMachine::Transition duplicate = dispatch_phase(machine,
	Event::WORK_PUMPED, Phase::FALLBACK, witness);
    if (duplicate.progressed || machine.stagnantDispatchCount() != 1 ||
	machine.dispatchSerial() != 2) {
	std::fprintf(stderr, "FAIL: stagnant dispatch accounting\n");
	return 1;
    }

    LodMachine::Inputs invalidInteraction = inputs_for_phase(Phase::FALLBACK);
    LodMachine::Transition invalid = machine.dispatch(
	Event::INTERACTION_STARTED, invalidInteraction, witness);
    if (!(invalid.invariantMask &
	    Invariant::INVARIANT_INTERACTION_EVENT_WITHOUT_INTERACTION) ||
	!invalid.canonicalized || invalid.current != Phase::INTERACTIVE ||
	!machine.lastInputs().interactive ||
	!machine.lastInputs().gestureActive ||
	machine.lastObservedInputs().interactive ||
	machine.lastObservedInputs().gestureActive ||
	machine.invariantViolationCount() != 1) {
	std::fprintf(stderr, "FAIL: interaction invariant\n");
	return 1;
    }

    LodMachine::Inputs invalidCompaction =
	inputs_for_phase(Phase::COMPACTING);
    invalidCompaction.coverageComplete = false;
    witness.visibleCount = 2;
    witness.completedCount = 3;
    invalid = machine.dispatch(Event::WORK_PUMPED, invalidCompaction,
	witness);
    if (!(invalid.invariantMask &
	    Invariant::INVARIANT_COMPACTION_WITHOUT_COVERAGE) ||
	!(invalid.invariantMask &
	    Invariant::INVARIANT_COMPLETED_EXCEEDS_VISIBLE) ||
	!invalid.canonicalized || invalid.current != Phase::FALLBACK ||
	machine.lastInputs().compacting ||
	machine.lastObservedWitness().completedCount != 3 ||
	machine.phaseWitness(Phase::FALLBACK).completedCount != 2 ||
	machine.invariantViolationCount() != 2 ||
	!(machine.invariantHistoryMask() &
	    Invariant::INVARIANT_INTERACTION_EVENT_WITHOUT_INTERACTION) ||
	!(machine.invariantHistoryMask() &
	    Invariant::INVARIANT_COMPACTION_WITHOUT_COVERAGE)) {
	std::fprintf(stderr, "FAIL: compaction/count invariants\n");
	return 1;
    }

    witness.completedCount = witness.visibleCount;
    LodMachine::Transition recovered = dispatch_phase(machine,
	Event::INTERACTION_STARTED, Phase::INTERACTIVE, witness);
    if (recovered.current != Phase::INTERACTIVE ||
	machine.lastInvariantMask() != Invariant::INVARIANT_NONE ||
	machine.stagnantDispatchCount() != 0 ||
	machine.invariantViolationCount() != 2 ||
	machine.invariantHistoryMask() == Invariant::INVARIANT_NONE) {
	std::fprintf(stderr, "FAIL: valid transition did not recover diagnostics\n");
	return 1;
    }

    LodMachine::Inputs cancelled = inputs_for_phase(Phase::COVERAGE);
    invalid = machine.dispatch(Event::GENERATION_CANCELLED, cancelled,
	witness);
    if (!(invalid.invariantMask &
	    Invariant::INVARIANT_CANCEL_EVENT_WITH_GENERATION) ||
	!invalid.canonicalized || invalid.current != Phase::FALLBACK ||
	machine.lastInputs().generationActive ||
	machine.invariantViolationCount() != 3) {
	std::fprintf(stderr, "FAIL: cancellation invariant\n");
	return 1;
    }

    LodMachine pendingMachine;
    LodMachine::Inputs omittedSettling = inputs_for_phase(Phase::STABLE);
    LodMachine::Witness pendingWitness;
    pendingWitness.visibleCount = 4;
    pendingWitness.completedCount = 4;
    pendingWitness.pendingCount = 1;
    LodMachine::Transition pending = pendingMachine.dispatch(
	Event::WORK_SCHEDULED, omittedSettling, pendingWitness);
    if (pending.current != Phase::SETTLING || !pending.canonicalized ||
	!pendingMachine.lastInputs().settlingWork ||
	!(pending.invariantMask &
	    Invariant::INVARIANT_STABLE_WITH_PENDING_WORK)) {
	std::fprintf(stderr, "FAIL: pending witness remained stable\n");
	return 1;
    }

    LodMachine::Inputs invalidEnd = inputs_for_phase(Phase::INTERACTIVE);
    invalid = machine.dispatch(Event::INTERACTION_ENDED, invalidEnd, witness);
    if (!(invalid.invariantMask &
	    Invariant::INVARIANT_INTERACTION_END_WITH_GESTURE) ||
	!invalid.canonicalized || invalid.current != Phase::INTERACTIVE ||
	!machine.lastInputs().interactive ||
	machine.lastInputs().gestureActive) {
	std::fprintf(stderr, "FAIL: interaction-end invariant\n");
	return 1;
    }

    LodMachine pressureMachine;
    LodMachine::Inputs invalidRecovery = inputs_for_phase(Phase::STABLE);
    invalidRecovery.resourceRecoveryPending = true;
    LodMachine::Transition recovery = pressureMachine.dispatch(
	Event::WORK_SCHEDULED, invalidRecovery, witness);
    if (!(recovery.invariantMask &
	    Invariant::INVARIANT_RESOURCE_RECOVERY_WITHOUT_PRESSURE) ||
	!recovery.canonicalized || recovery.current != Phase::STABLE ||
	pressureMachine.lastInputs().resourceRecoveryPending) {
	std::fprintf(stderr, "FAIL: resource recovery invariant\n");
	return 1;
    }

    return 0;
}

static int
test_resource_policy(void)
{
    BObolLodResourcePolicy policy;
    BObolLodResourcePolicy::Decision decision =
	policy.observe(false, false, true);
    if (decision.changed || decision.queueRecovery || policy.anyPressure() ||
	policy.recoveryPending() || policy.revision() != 0) {
	std::fprintf(stderr, "FAIL: initial resource observation\n");
	return 1;
    }

    decision = policy.observe(true, false, true);
    if (!decision.changed || !decision.queueRecovery ||
	!policy.cpuPressure() || policy.gpuPressure() ||
	!policy.recoveryPending() || policy.revision() != 1 ||
	policy.recoveryHandledRevision() != 0) {
	std::fprintf(stderr, "FAIL: CPU pressure edge\n");
	return 1;
    }

    decision = policy.observe(true, false, true);
    if (decision.changed || decision.queueRecovery || policy.revision() != 1) {
	std::fprintf(stderr, "FAIL: duplicate pressure edge\n");
	return 1;
    }

    policy.markRecoveryHandled();
    if (policy.recoveryPending() || !policy.anyPressure() ||
	policy.recoveryHandledRevision() != policy.revision()) {
	std::fprintf(stderr, "FAIL: handled persistent pressure\n");
	return 1;
    }

    decision = policy.observe(true, true, true);
    if (!decision.changed || !decision.queueRecovery ||
	!policy.recoveryPending() || policy.revision() != 2) {
	std::fprintf(stderr, "FAIL: GPU pressure edge\n");
	return 1;
    }

    decision = policy.observe(false, false, true);
    if (!decision.changed || !decision.pressureCleared ||
	decision.queueRecovery || policy.anyPressure() ||
	policy.recoveryPending() || policy.revision() != 3 ||
	policy.recoveryHandledRevision() != 3) {
	std::fprintf(stderr, "FAIL: resource pressure clear\n");
	return 1;
    }

    (void)policy.observe(true, false, true);
    policy.resetForServiceChange();
    if (policy.anyPressure() || policy.recoveryPending() ||
	policy.recoveryHandledRevision() != policy.revision()) {
	std::fprintf(stderr, "FAIL: resource service reset\n");
	return 1;
    }

    decision = policy.observe(true, false, false);
    if (!decision.changed || decision.queueRecovery ||
	policy.recoveryPending()) {
	std::fprintf(stderr, "FAIL: disabled recovery pressure edge\n");
	return 1;
    }
    decision = policy.observe(true, false, true);
    if (decision.changed || !decision.queueRecovery ||
	!policy.recoveryPending()) {
	std::fprintf(stderr, "FAIL: late resource recovery enable\n");
	return 1;
    }
    policy.markRecoveryHandled();
    decision = policy.observe(true, false, true);
    if (decision.changed || decision.queueRecovery ||
	policy.recoveryPending()) {
	std::fprintf(stderr, "FAIL: handled resource recovery repeated\n");
	return 1;
    }
    return 0;
}

static int
test_headroom_policy(void)
{
    BObolLodHeadroomPolicy policy;
    const BObolLodViewEpoch view(4);
    const BObolLodPolicyEpoch lodPolicy(9);
    if (policy.retryPending() ||
	!policy.armRetry(view, lodPolicy, 100) ||
	!policy.retryPending() ||
	!policy.pendingMatches(view, lodPolicy, 100) ||
	policy.pendingMatches(BObolLodViewEpoch(5), lodPolicy, 100) ||
	policy.pendingMatches(view, lodPolicy, 101) ||
	policy.armRetry(view, lodPolicy, 100) ||
	!policy.consumeRetry(view, lodPolicy, 100, 106.0L, 90, 100,
	    true) || policy.retryPending() ||
	policy.armRetry(view, lodPolicy, 100)) {
	std::fprintf(stderr, "FAIL: explicit headroom witness lifecycle\n");
	return 1;
    }
    if (!policy.armRetry(view, lodPolicy, 120) ||
	!policy.consumeRetry(view, lodPolicy, 120, 200.0L, 10, 100,
	    true) ||
	!policy.armRetry(BObolLodViewEpoch(5), lodPolicy, 120) ||
	!policy.consumeRetry(BObolLodViewEpoch(5), lodPolicy, 120,
	    200.0L, 10, 100, true) ||
	!policy.armRetry(BObolLodViewEpoch(5), BObolLodPolicyEpoch(10),
	    120) ||
	!policy.consumeRetry(BObolLodViewEpoch(5),
	    BObolLodPolicyEpoch(10), 120, 200.0L, 10, 100, true)) {
	std::fprintf(stderr, "FAIL: headroom progress witnesses\n");
	return 1;
    }
    if (!policy.armRetry(BObolLodViewEpoch(6), lodPolicy, 200) ||
	policy.consumeRetry(BObolLodViewEpoch(6), lodPolicy, 200,
	    210.0L, 10, 100, true) ||
	policy.armRetry(BObolLodViewEpoch(6), lodPolicy, 200)) {
	std::fprintf(stderr, "FAIL: exact threshold headroom admitted\n");
	return 1;
    }
    if (!policy.armRetry(BObolLodViewEpoch(7), lodPolicy, 200) ||
	policy.consumeRetry(BObolLodViewEpoch(7), lodPolicy, 200,
	    300.0L, 120, 100, true)) {
	std::fprintf(stderr, "FAIL: slow headroom sample admitted\n");
	return 1;
    }
    if (!policy.armRetry(BObolLodViewEpoch(8), lodPolicy, 200) ||
	policy.consumeRetry(BObolLodViewEpoch(8), lodPolicy, 200,
	    300.0L, 10, 100, false)) {
	std::fprintf(stderr, "FAIL: non-reusable headroom sample admitted\n");
	return 1;
    }
    if (!policy.armRetry(BObolLodViewEpoch(9), lodPolicy, 200) ||
	policy.consumeRetry(BObolLodViewEpoch(10), lodPolicy, 200,
	    300.0L, 10, 100, true) || policy.retryPending() ||
	!policy.armRetry(BObolLodViewEpoch(9), lodPolicy, 200)) {
	std::fprintf(stderr, "FAIL: stale headroom witness handling\n");
	return 1;
    }
    policy.cancelRetry();
    if (policy.retryPending() || policy.pendingMatches(
	    BObolLodViewEpoch(9), lodPolicy, 200)) {
	std::fprintf(stderr, "FAIL: headroom witness cancellation\n");
	return 1;
    }
    return 0;
}

static int
test_coverage_policy(void)
{
    using Policy = BObolLodCoveragePolicy;
    Policy policy;
    if (!policy.active() || policy.sawBoundedSource() ||
	policy.coverageComplete() || policy.visibleCount() != 0 ||
	policy.coveredCount() != 0 || policy.hasCompleteVisibleCount() ||
	policy.required() || policy.effectiveActive() ||
	!policy.effectiveComplete() || !policy.compactionAllowed()) {
	std::fprintf(stderr, "FAIL: initial coverage policy\n");
	return 1;
    }
    policy.setRequired(true);
    if (!policy.required() || !policy.effectiveActive() ||
	policy.effectiveComplete() || policy.compactionAllowed()) {
	std::fprintf(stderr, "FAIL: required coverage gate\n");
	return 1;
    }

    policy.markBoundedSource();
    policy.observe(4, 3);
    policy.observe(SIZE_MAX, SIZE_MAX);
    if (!policy.sawBoundedSource() || policy.visibleCount() != SIZE_MAX ||
	policy.coveredCount() != SIZE_MAX) {
	std::fprintf(stderr, "FAIL: bounded coverage accumulation\n");
	return 1;
    }
    Policy::Completion completion = policy.completeIfReady(false, false);
    if (completion.completed || !policy.active()) {
	std::fprintf(stderr, "FAIL: incomplete coverage pass retired\n");
	return 1;
    }
    completion = policy.completeIfReady(true, true);
    if (completion.completed || !policy.active()) {
	std::fprintf(stderr, "FAIL: rescan coverage pass retired\n");
	return 1;
    }

    policy.clearPassCounters();
    policy.observe(10, 7);
    policy.markBoundedSource();
    completion = policy.completeIfReady(true, false);
    if (!completion.completed || !completion.bounded ||
	!completion.missing || completion.visibleCount != 10 ||
	completion.coveredCount != 7 || policy.active() ||
	policy.coverageComplete() || !policy.hasCompleteVisibleCount() ||
	policy.completeVisibleCount() != 10 || policy.sawBoundedSource() ||
	policy.visibleCount() != 0 || policy.coveredCount() != 0) {
	std::fprintf(stderr, "FAIL: missing bounded coverage proof\n");
	return 1;
    }

    /* A quality-only pass preserves the latest view denominator and its
     * minimum-mesh proof while collecting a fresh set of counters. */
    policy.activate(false);
    policy.observe(5, 5);
    completion = policy.completeIfReady(true, false);
    if (!completion.completed || completion.bounded || completion.missing ||
	!policy.coverageComplete() || policy.completeVisibleCount() != 5) {
	std::fprintf(stderr, "FAIL: complete unbounded coverage proof\n");
	return 1;
    }

    policy.activate(true);
    if (!policy.active() || policy.coverageComplete() ||
	!policy.hasCompleteVisibleCount()) {
	std::fprintf(stderr, "FAIL: view invalidation semantics\n");
	return 1;
    }
    policy.clearCompleteVisibleCount();
    if (policy.hasCompleteVisibleCount() ||
	policy.completeVisibleCount() != 0) {
	std::fprintf(stderr, "FAIL: visible denominator invalidation\n");
	return 1;
    }
    policy.reset();
    if (policy.active() || policy.coverageComplete() ||
	policy.sawBoundedSource() || policy.hasCompleteVisibleCount() ||
	!policy.required() || policy.effectiveComplete() ||
	policy.compactionAllowed()) {
	std::fprintf(stderr, "FAIL: coverage reset\n");
	return 1;
    }
    policy.setRequired(false);
    if (policy.effectiveActive() || !policy.effectiveComplete() ||
	!policy.compactionAllowed()) {
	std::fprintf(stderr, "FAIL: optional coverage gate\n");
	return 1;
    }
    return 0;
}

static int
test_convergence_policy(void)
{
    using Policy = BObolLodConvergencePolicy;
    Policy policy;
    const auto baseInput = []() {
	Policy::Inputs input;
	input.viewEpoch.set(1);
	input.policyEpoch.set(1);
	return input;
    };
    Policy::Inputs input = baseInput();
    Policy::Decision decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::IDLE ||
	decision.fraction < 0.999f || decision.fraction > 1.001f ||
	!decision.viewReady ||
	decision.backgroundPending || decision.visualPending ||
	decision.hasLodState) {
	std::fprintf(stderr, "FAIL: empty convergence state\n");
	return 1;
    }

    input = baseInput();
    input.expectedLeafCount = 10;
    input.availableLeafCount = 4;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::DISCOVERING ||
	decision.viewReady || !decision.visualPending ||
	decision.fraction < 0.159f || decision.fraction > 0.161f) {
	std::fprintf(stderr, "FAIL: structural convergence phase\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(2);
    input.structuralDiscovery = true;
    input.progressiveWorkPending = true;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::DISCOVERING ||
	!decision.backgroundPending || decision.fraction > 0.001f) {
	std::fprintf(stderr, "FAIL: provider discovery phase\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(3);
    input.interactive = true;
    input.visibleTargetCount = 4;
    input.activePayloadCount = 4;
    input.satisfiedPayloadCount = 2;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::INTERACTIVE ||
	decision.viewReady || decision.fraction < 0.624f ||
	decision.fraction > 0.626f) {
	std::fprintf(stderr, "FAIL: interactive convergence phase\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(4);
    input.calibrationPending = true;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::CALIBRATING ||
	decision.viewReady || decision.fraction < 0.949f ||
	decision.fraction > 0.951f) {
	std::fprintf(stderr, "FAIL: calibration convergence phase\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(5);
    input.resultPending = true;
    input.visibleTargetCount = 4;
    input.activePayloadCount = 4;
    input.satisfiedPayloadCount = 2;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::REFINING ||
	decision.viewReady || !decision.visualPending ||
	decision.fraction < 0.674f || decision.fraction > 0.676f) {
	std::fprintf(stderr, "FAIL: refining convergence phase\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(6);
    input.activePayloadCount = 1;
    input.queuedCacheWrites = 1;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::BACKGROUND ||
	!decision.viewReady || !decision.backgroundPending ||
	decision.fraction < 0.999f || decision.fraction > 1.001f) {
	std::fprintf(stderr, "FAIL: background convergence phase\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(7);
    input.visibleTargetCount = 2;
    input.activePayloadCount = 1;
    input.satisfiedPayloadCount = 1;
    input.memoryLimitedPayloadCount = 1;
    decision = policy.evaluate(input);
    if (!decision.viewReady || !decision.memoryLimited ||
	!decision.performanceLimited ||
	decision.phase != Policy::Phase::IDLE) {
	std::fprintf(stderr, "FAIL: limited terminal convergence\n");
	return 1;
    }

    input.failedSourceCount = 1;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::ERROR) {
	std::fprintf(stderr, "FAIL: convergence error phase\n");
	return 1;
    }

    input = baseInput();
    input.viewEpoch.set(8);
    input.calibrationPending = true;
    decision = policy.evaluate(input);
    if (decision.fraction < 0.949f || decision.fraction > 0.951f) {
	std::fprintf(stderr, "FAIL: convergence fraction seed\n");
	return 1;
    }
    input.calibrationPending = false;
    input.expectedLeafCount = 10;
    input.availableLeafCount = 1;
    decision = policy.evaluate(input);
    if (decision.phase != Policy::Phase::DISCOVERING ||
	decision.fraction < 0.949f || decision.fraction > 0.951f ||
	policy.fractionFloor() < 0.949f ||
	policy.fractionFloor() > 0.951f) {
	std::fprintf(stderr, "FAIL: monotonic convergence fraction\n");
	return 1;
    }
    input.viewEpoch.set(9);
    decision = policy.evaluate(input);
    if (decision.fraction < 0.039f || decision.fraction > 0.041f) {
	std::fprintf(stderr, "FAIL: convergence epoch reset\n");
	return 1;
    }
    policy.resetFraction();
    if (policy.fractionFloor() < -0.001f ||
	policy.fractionFloor() > 0.001f) {
	std::fprintf(stderr, "FAIL: explicit convergence reset\n");
	return 1;
    }
    return 0;
}

static int
test_budget_policy(void)
{
    using Policy = BObolLodBudgetPolicy;
    Policy policy;
    Policy::Inputs input;
    Policy::Decision decision = policy.beginPass(input);
    if (!decision.initialized || decision.totalBudget != 50000 ||
	decision.refinementBudget != 50000 || decision.retainedAdmission ||
	!policy.passInitialized() || policy.currentBudget() != 50000) {
	std::fprintf(stderr, "FAIL: seed scene budget\n");
	return 1;
    }
    decision = policy.beginPass(input);
    if (decision.initialized || decision.totalBudget != 50000) {
	std::fprintf(stderr, "FAIL: duplicate scene budget initialization\n");
	return 1;
    }
    policy.consumeRefinement(12000);
    if (policy.refinementRemaining() != 38000) {
	std::fprintf(stderr, "FAIL: refinement budget consumption\n");
	return 1;
    }

    policy.resetPass();
    input.activeCost = 100000;
    input.targetFps = 20.0f;
    input.calibratedCostPerSecond = 10000000.0L;
    input.observedStableNanoseconds = 10000000ULL;
    decision = policy.beginPass(input);
    if (!decision.initialized || decision.totalBudget != 400000 ||
	decision.refinementBudget != 300000 ||
	decision.retainedAdmission) {
	std::fprintf(stderr, "FAIL: bounded fast-frame budget growth\n");
	return 1;
    }

    policy.resetPass();
    input.activeCost = 400000;
    input.calibratedCostPerSecond = 4000000.0L;
    input.observedStableNanoseconds = 70000000ULL;
    policy.setProbeCount(3);
    decision = policy.beginPass(input);
    if (!decision.overloadRecovery || decision.totalBudget != 200000 ||
	!decision.retainedAdmission ||
	decision.retainedAdmissionBudget != 200000 ||
	!policy.overloadRecoveryPerformed() ||
	policy.overloadRecoveryActiveCost() != 400000) {
	std::fprintf(stderr, "FAIL: stable overload recovery\n");
	return 1;
    }
    policy.consumeRetainedAdmission(75000);
    if (policy.retainedAdmissionRemaining() != 125000) {
	std::fprintf(stderr, "FAIL: retained admission consumption\n");
	return 1;
    }

    policy.reset();
    input = Policy::Inputs();
    input.interactive = true;
    input.activeCost = 100000;
    input.targetFps = 60.0f;
    input.lastRenderNanoseconds = 4000000ULL;
    decision = policy.beginPass(input);
    if (decision.totalBudget < 332000 || decision.totalBudget > 334000 ||
	decision.refinementBudget < 232000 ||
	decision.refinementBudget > 234000) {
	std::fprintf(stderr, "FAIL: interactive retained-cut bootstrap\n");
	return 1;
    }

    policy.reset();
    input = Policy::Inputs();
    input.interactive = true;
    input.scaleQualityProbe = true;
    input.activeCost = 200000;
    input.targetFps = 10.0f;
    input.lastRenderNanoseconds = 40000000ULL;
    decision = policy.beginPass(input);
    if (decision.totalBudget != 800000 ||
	decision.refinementBudget != 600000) {
	std::fprintf(stderr, "FAIL: scale-quality discrete probe\n");
	return 1;
    }

    policy.reset();
    input = Policy::Inputs();
    input.activeCost = 100000;
    input.targetFps = 20.0f;
    input.calibratedCostPerSecond = 1000000.0L;
    input.stablePresentationHandoff = true;
    input.stablePresentationCostFloor = 90000;
    decision = policy.beginPass(input);
    if (decision.totalBudget != 90000 ||
	!decision.retainedAdmission ||
	decision.retainedAdmissionBudget != 90000) {
	std::fprintf(stderr, "FAIL: proven presentation cost floor\n");
	return 1;
    }

    policy.resetPass();
    input = Policy::Inputs();
    input.forceTerminal = true;
    input.activeCost = 123;
    decision = policy.beginPass(input);
    if (decision.totalBudget != SIZE_MAX ||
	decision.refinementBudget != SIZE_MAX ||
	decision.retainedAdmission) {
	std::fprintf(stderr, "FAIL: unbounded terminal policy\n");
	return 1;
    }

    policy.reduceCurrentBudget(90000);
    if (policy.currentBudget() != 90000) {
	std::fprintf(stderr, "FAIL: deadline budget backoff\n");
	return 1;
    }
    policy.raiseCurrentBudget(120000);
    policy.raiseCurrentBudget(100000);
    if (policy.currentBudget() != 120000) {
	std::fprintf(stderr, "FAIL: stable budget restoration\n");
	return 1;
    }
    policy.resetOverloadRecovery();
    if (policy.overloadRecoveryPerformed() ||
	policy.overloadRecoveryActiveCost() != 0) {
	std::fprintf(stderr, "FAIL: overload witness reset\n");
	return 1;
    }

    policy.reset();
    Policy::CalibrationInputs calibrationInput;
    calibrationInput.activeCost = 100000;
    calibrationInput.targetNanoseconds = 50000000ULL;
    calibrationInput.observedNanoseconds = 10000000ULL;
    calibrationInput.calibratedBudget = 100000.0L;
    Policy::CalibrationDecision calibration =
	policy.finishBlockedPass(calibrationInput);
    if (!calibration.probeCandidate || !calibration.probeEligible ||
	!calibration.requestFrame || !calibration.calibrationFrame ||
	policy.probeActiveCost() != 100000 || policy.probeCount() != 1 ||
	policy.calibrationFramesRemaining() != 11 ||
	!policy.rescanAfterFrame() || policy.stableBudgetLimited()) {
	std::fprintf(stderr, "FAIL: stable headroom probe admission\n");
	return 1;
    }
    for (unsigned int i = 0; i < 11; i++) {
	Policy::CompletedFrameDecision completed =
	    policy.completeCalibrationFrame();
	if (!completed.requestCalibrationFrame ||
	    completed.restartSubmission) {
	    std::fprintf(stderr, "FAIL: bounded calibration series\n");
	    return 1;
	}
    }
    Policy::CompletedFrameDecision completed =
	policy.completeCalibrationFrame();
    if (completed.requestCalibrationFrame ||
	!completed.restartSubmission || policy.rescanAfterFrame() ||
	policy.probeCount() != 12 || policy.passInitialized()) {
	std::fprintf(stderr, "FAIL: calibration series completion\n");
	return 1;
    }

    policy.resetCalibration();
    calibrationInput.observedNanoseconds = 46000000ULL;
    calibrationInput.calibratedBudget = 50000.0L;
    calibration = policy.finishBlockedPass(calibrationInput);
    if (calibration.probeCandidate || !calibration.calibrationFrame ||
	policy.probeCount() != 1 ||
	policy.calibrationFramesRemaining() != 2) {
	std::fprintf(stderr, "FAIL: minimum calibration sample series\n");
	return 1;
    }
    for (unsigned int i = 0; i < 2; i++) {
	completed = policy.completeCalibrationFrame();
	if (!completed.requestCalibrationFrame ||
	    completed.restartSubmission) {
	    std::fprintf(stderr, "FAIL: minimum calibration frame\n");
	    return 1;
	}
    }
    completed = policy.completeCalibrationFrame();
    if (!completed.restartSubmission || policy.rescanAfterFrame()) {
	std::fprintf(stderr, "FAIL: minimum calibration rescan\n");
	return 1;
    }
    calibration = policy.finishBlockedPass(calibrationInput);
    if (calibration.requestFrame || calibration.calibrationFrame ||
	!policy.stableBudgetLimited() || policy.probeCount() != 0) {
	std::fprintf(stderr, "FAIL: stable calibrated terminal cut\n");
	return 1;
    }

    policy.resetCalibration();
    calibrationInput.passAdmittedWork = true;
    calibration = policy.finishBlockedPass(calibrationInput);
    if (!calibration.requestFrame || calibration.calibrationFrame ||
	!policy.rescanAfterFrame() || policy.stableBudgetLimited()) {
	std::fprintf(stderr, "FAIL: admitted-work frame barrier\n");
	return 1;
    }
    completed = policy.completeCalibrationFrame();
    if (!completed.restartSubmission || policy.rescanAfterFrame()) {
	std::fprintf(stderr, "FAIL: admitted-work rescan\n");
	return 1;
    }

    policy.requestRescanAfterFrame(true);
    if (!policy.rescanAfterFrame() ||
	policy.calibrationFramesRemaining() != 0 ||
	policy.stableBudgetLimited()) {
	std::fprintf(stderr, "FAIL: explicit calibration frame barrier\n");
	return 1;
    }
    policy.resetCalibration();
    if (policy.rescanAfterFrame() || policy.stableBudgetLimited() ||
	policy.probeActiveCost() != 0 || policy.probeCount() != 0) {
	std::fprintf(stderr, "FAIL: full calibration reset\n");
	return 1;
    }
    return 0;
}

static int
test_quality_policy(void)
{
    using Policy = BObolLodQualityPolicy;
    const auto near = [](float left, float right) {
	return std::fabs(left - right) <= 0.0001f;
    };
    if (!near(Policy::interactivePixelError(0, 60.0f), 4.0f) ||
	!near(Policy::interactivePixelError(10000000ULL, 60.0f), 1.0f) ||
	!near(Policy::interactivePixelError(10000000ULL, 0.0f), 1.0f) ||
	!near(Policy::interactivePixelError(UINT64_MAX, 60.0f), 64.0f)) {
	std::fprintf(stderr, "FAIL: interactive pixel-error policy\n");
	return 1;
    }
    const float error =
	Policy::interactivePixelError(66666668ULL, 60.0f);
    if (error < 3.99f || error > 4.01f) {
	std::fprintf(stderr, "FAIL: measured pixel-error scaling\n");
	return 1;
    }


    if (!near(Policy::pointProxyThreshold(
	    0.0f, 10000000ULL, 60.0f), 1.0f) ||
	!near(Policy::pointProxyThreshold(
	    70.0f, 10000000ULL, 60.0f), 64.0f) ||
	!near(Policy::pointProxyThreshold(
	    std::numeric_limits<float>::quiet_NaN(),
	    10000000ULL, 60.0f), 1.0f)) {
	std::fprintf(stderr, "FAIL: point-proxy threshold bounds\n");
	return 1;
    }
    const float proxy =
	Policy::pointProxyThreshold(4.0f, 66666668ULL, 60.0f);
    if (proxy < 7.99f || proxy > 8.01f) {
	std::fprintf(stderr, "FAIL: cumulative point-proxy correction\n");
	return 1;
    }

    if (Policy::progressiveCeiling(-1, 4.0f) != -1 ||
	Policy::progressiveCeiling(10, 1.0f) != -1 ||
	Policy::progressiveCeiling(10, 3.99f) != 8 ||
	Policy::progressiveCeiling(10, 2.0f, 7) != 6 ||
	Policy::progressiveCeiling(3, 64.0f) != 0 ||
	Policy::progressiveCeiling(
	    3, std::numeric_limits<float>::infinity()) != 0 ||
	Policy::progressiveCeiling(
	    3, std::numeric_limits<float>::quiet_NaN()) != -1) {
	std::fprintf(stderr, "FAIL: progressive ceiling policy\n");
	return 1;
    }

    const Policy::StableCapacityEvidence exactStable =
	Policy::stableCapacityEvidence(1000, 50000000ULL, 20.0f, true);
    const Policy::StableCapacityEvidence fasterStable =
	Policy::stableCapacityEvidence(1000, 40000000ULL, 20.0f, true);
    if (!exactStable.proven ||
	std::fabs(exactStable.renderCostPerSecond - 20000.0L) > 0.001L ||
	!fasterStable.proven ||
	std::fabs(fasterStable.renderCostPerSecond - 25000.0L) > 0.001L ||
	Policy::stableCapacityEvidence(
	    1000, 50000001ULL, 20.0f, true).proven ||
	Policy::stableCapacityEvidence(
	    1000, 40000000ULL, 20.0f, false).proven ||
	Policy::stableCapacityEvidence(
	    0, 40000000ULL, 20.0f, true).proven ||
	Policy::stableCapacityEvidence(
	    1000, 0, 20.0f, true).proven ||
	Policy::stableCapacityEvidence(
	    1000, 40000000ULL,
	    std::numeric_limits<float>::infinity(), true).proven ||
	Policy::stableCapacityEvidence(
	    1000, 40000000ULL,
	    std::numeric_limits<float>::quiet_NaN(), true).proven) {
	std::fprintf(stderr, "FAIL: reusable quality-to-stable evidence\n");
	return 1;
    }

    const size_t lucyTrial =
	Policy::discreteTrialOverBudgetAllowance(222780, 431301);
    const size_t overflowTrial =
	Policy::discreteTrialOverBudgetAllowance(
	    SIZE_MAX / 2, SIZE_MAX / 2);
    if (lucyTrial != 431301 ||
	Policy::discreteTrialOverBudgetAllowance(100, 150) != 150 ||
	Policy::discreteTrialOverBudgetAllowance(100, 500) != 0 ||
	Policy::discreteTrialOverBudgetAllowance(0, 100) != 0 ||
	Policy::discreteTrialOverBudgetAllowance(100, 0) != 0 ||
	Policy::discreteTrialOverBudgetAllowance(101, 100) != 0 ||
	Policy::discreteTrialOverBudgetAllowance(100, SIZE_MAX) != 0 ||
	overflowTrial == 0 || overflowTrial > SIZE_MAX / 2) {
	std::fprintf(stderr, "FAIL: bounded discrete population trial policy "
	    "(lucy=%zu overflow=%zu)\n", lucyTrial, overflowTrial);
	return 1;
    }
    return 0;
}

static int
test_view_demand_policy(void)
{
    using Policy = BObolLodViewDemandPolicy;
    Policy policy;
    if (policy.scaleDemandRefreshActive() || policy.viewScaleChanging() ||
	policy.interactionScaleChanged() || policy.qualityProbeActive() ||
	policy.qualityProbePending() || policy.qualityProbePresented() ||
	policy.qualityBudgetActive() ||
	policy.scaleChangingInteraction(true)) {
	std::fprintf(stderr, "FAIL: initial view-demand state\n");
	return 1;
    }

    Policy::CameraChangeDecision camera =
	policy.observeCameraChange(true, 10000000ULL);
    policy.refreshForViewRevision(true);
    policy.beginCameraInteraction(true, true);
    if (camera.retainPriorQualityCeiling ||
	!policy.viewScaleChanging() ||
	!policy.interactionScaleChanged() ||
	!policy.scaleDemandRefreshActive() ||
	!policy.scaleChangingInteraction(true) ||
	!policy.noteMotionFrameSettled() ||
	!policy.qualityProbePending()) {
	std::fprintf(stderr, "FAIL: scale-demand camera edge\n");
	return 1;
    }

    Policy::QualityProbeInputs probeInput;
    probeInput.activeMaximum = 6;
    probeInput.motionFramePending = true;
    Policy::QualityProbeDecision probe =
	policy.beginQualityProbe(probeInput);
    if (probe.consumed || probe.begin ||
	!policy.qualityProbePending()) {
	std::fprintf(stderr, "FAIL: quality probe bypassed frame barrier\n");
	return 1;
    }
    probeInput.motionFramePending = false;
    probe = policy.beginQualityProbe(probeInput);
    if (!probe.consumed || !probe.begin ||
	probe.progressiveCeiling != 7 ||
	!policy.qualityProbeActive() ||
	policy.qualityProbePending() ||
	policy.qualityProbePresented() ||
	!policy.qualityBudgetActive()) {
	std::fprintf(stderr, "FAIL: bounded zoom-quality probe\n");
	return 1;
    }

    policy.noteFramePresented();
    camera = policy.observeCameraChange(true, 100000000ULL);
    if (!camera.retainPriorQualityCeiling ||
	policy.qualityProbeActive() || policy.qualityProbePending() ||
	policy.qualityProbePresented() ||
	!policy.qualityBudgetActive() || !policy.viewScaleChanging()) {
	std::fprintf(stderr, "FAIL: responsive quality-floor retention\n");
	return 1;
    }
    policy.noteFramePresented();
    camera = policy.observeCameraChange(true, 100000001ULL);
    if (camera.retainPriorQualityCeiling) {
	std::fprintf(stderr, "FAIL: slow quality floor retained\n");
	return 1;
    }

    policy.beginCameraInteraction(false, true);
    if (!policy.rearmAfterResidentGrowth(true) ||
	policy.qualityProbeActive() || !policy.qualityProbePending() ||
	policy.qualityProbePresented() || !policy.qualityBudgetActive()) {
	std::fprintf(stderr, "FAIL: resident-growth quality rearm\n");
	return 1;
    }
    probeInput.activeMaximum = 15;
    probeInput.presentationCeiling = 6;
    probe = policy.beginQualityProbe(probeInput);
    if (!probe.begin || probe.progressiveCeiling != 7) {
	std::fprintf(stderr, "FAIL: quality probe bypassed presented cut\n");
	return 1;
    }
    policy.endGesture();
    if (policy.qualityProbeActive() || policy.qualityProbePending() ||
	policy.qualityProbePresented() || !policy.qualityBudgetActive()) {
	std::fprintf(stderr, "FAIL: gesture-end view-demand state\n");
	return 1;
    }

    /* A bracketed rotation may begin before an unbracketed wheel epoch's
     * quiet debounce.  The earlier scale change remains a stable-convergence
     * obligation, but it must not classify the rotation frame or a late
     * resident suffix as another scale-quality probe. */
    policy.beginGesture(false);
    camera = policy.observeCameraChange(false, 10000000ULL);
    policy.beginCameraInteraction(false, false);
    if (camera.retainPriorQualityCeiling ||
	!policy.interactionScaleChanged() ||
	policy.scaleChangingInteraction(true) ||
	policy.rearmAfterResidentGrowth(true)) {
	std::fprintf(stderr, "FAIL: mixed zoom/pose demand isolation\n");
	return 1;
    }
    if (!policy.finishQuietInteraction() ||
	policy.viewScaleChanging() || policy.interactionScaleChanged() ||
	policy.qualityProbeActive() || policy.qualityProbePending() ||
	policy.qualityProbePresented() || policy.qualityBudgetActive()) {
	std::fprintf(stderr, "FAIL: quiet scale handoff\n");
	return 1;
    }

    policy.reset();
    (void)policy.observeCameraChange(true, 10000000ULL);
    policy.beginCameraInteraction(true, true);
    if (policy.finishQuietInteraction(true)) {
	std::fprintf(stderr, "FAIL: scale round trip required retarget\n");
	return 1;
    }

    policy.refreshForPolicyRevision(true, false);
    if (!policy.scaleDemandRefreshActive()) {
	std::fprintf(stderr, "FAIL: scale-demand policy continuation\n");
	return 1;
    }
    policy.clearDemandRefresh();
    policy.beginGesture(true);
    camera = policy.observeCameraChange(false, 10000000ULL);
    policy.refreshForViewRevision(true);
    policy.beginCameraInteraction(true, false);
    if (camera.retainPriorQualityCeiling ||
	policy.scaleDemandRefreshActive() || policy.viewScaleChanging() ||
	policy.interactionScaleChanged() ||
	policy.scaleChangingInteraction(true) ||
	policy.rearmAfterResidentGrowth(true)) {
	std::fprintf(stderr, "FAIL: pose-only demand isolation\n");
	return 1;
    }

    policy.reset();
    (void)policy.observeCameraChange(true, 0);
    policy.beginCameraInteraction(true, true);
    (void)policy.noteMotionFrameSettled();
    probeInput = Policy::QualityProbeInputs();
    probe = policy.beginQualityProbe(probeInput);
    if (!probe.consumed || probe.begin || policy.qualityProbePending()) {
	std::fprintf(stderr, "FAIL: unavailable population probe retirement\n");
	return 1;
    }
    policy.reset();
    if (policy.scaleDemandRefreshActive() || policy.viewScaleChanging() ||
	policy.interactionScaleChanged() || policy.qualityBudgetActive()) {
	std::fprintf(stderr, "FAIL: view-demand reset\n");
	return 1;
    }
    return 0;
}

static int
test_presentation_policy(void)
{
    using Policy = BObolLodPresentationPolicy;
    Policy policy;
    Policy::Population population;
    population.available = true;
    population.activeFaces = 100;
    population.cadRevision = 3;

    Policy::RestoreInputs input;
    input.orthographic = true;
    input.retainedMeshPayloads = true;
    input.viewEpoch.set(11);
    input.population = population;
    input.currentProgressiveCeiling = 4;
    input.currentPointProxyPixelThreshold = 3.0f;

    policy.capturePrior(8, 2.0f, population, BObolLodViewEpoch(10));
    Policy::RestoreDecision restore = policy.restorePrior(input);
    if (!restore.apply || !restore.restoredPriorStable ||
	restore.progressiveCeiling != 8 ||
	std::fabs(restore.pointProxyPixelThreshold - 2.0f) > 0.0001f ||
	!policy.priorRestored()) {
	std::fprintf(stderr, "FAIL: pose-only stable presentation restore\n");
	return 1;
    }
    restore = policy.beginQuiet(input);
	if (!restore.apply || !restore.restoredPriorStable ||
	restore.progressiveCeiling != 8 ||
	std::fabs(restore.pointProxyPixelThreshold - 2.0f) > 0.0001f ||
	restore.startHandoff ||
	policy.handoffPending() || policy.handoffCostFloor() != 0) {
	std::fprintf(stderr,
	    "FAIL: quiet debounce did not reapply restored presentation\n");
	return 1;
    }

    policy.reset();
    policy.capturePrior(8, 2.0f, population, BObolLodViewEpoch(10));
    input.population.activeFaces++;
    restore = policy.beginQuiet(input);
    if (restore.apply || restore.restoredPriorStable ||
	!restore.startHandoff || !policy.handoffPending() ||
	policy.handoffCostFloor() != 0) {
	std::fprintf(stderr, "FAIL: changed population presentation handoff\n");
	return 1;
    }
    Policy::CompletedPassInputs completed;
    completed.completed = true;
    completed.changedCut = true;
    if (policy.completePass(completed).finishHandoff) {
	std::fprintf(stderr, "FAIL: changed cut completed handoff\n");
	return 1;
    }
    completed.changedCut = false;
    completed.retainedRefinementPending = true;
    Policy::CompletedPassDecision completion = policy.completePass(completed);
    if (!completion.finishHandoff || !completion.requestRetainedRescan ||
	completion.retireRetainedObservation || policy.handoffPending()) {
	std::fprintf(stderr, "FAIL: unchanged pass did not complete handoff\n");
	return 1;
    }

    /* A quiet render-deadline correction has not presented its constrained
     * cut yet.  An unchanged planning pass cannot remove that ceiling; only a
     * completed presentation supplies the missing proof. */
    policy.reset();
    policy.armHandoff(true);
    completed = Policy::CompletedPassInputs();
    completed.completed = true;
    completion = policy.completePass(completed);
    if (completion.finishHandoff || !policy.handoffPending() ||
	!policy.handoffPresentationPending()) {
	std::fprintf(stderr, "FAIL: unpresented deadline handoff completed\n");
	return 1;
    }
    if (!policy.noteFramePresented()) {
	std::fprintf(stderr, "FAIL: deadline frame did not release handoff barrier\n");
	return 1;
    }
    completion = policy.completePass(completed);
    if (!completion.finishHandoff || policy.handoffPending() ||
	policy.handoffPresentationPending()) {
	std::fprintf(stderr, "FAIL: presented deadline handoff not completed\n");
	return 1;
    }

    /* A richer retained target with no handoff, cut, budget barrier, or
     * rescan is terminal quality state, not self-sustaining background work. */
    completed = Policy::CompletedPassInputs();
    completed.completed = true;
    completed.retainedRefinementPending = true;
    completion = policy.completePass(completed);
    if (!completion.retireRetainedObservation ||
	completion.finishHandoff || completion.requestRetainedRescan) {
	std::fprintf(stderr, "FAIL: unwitnessed retained observation not retired\n");
	return 1;
    }
    completed.retainedRefinementBudgetBlocked = true;
    completion = policy.completePass(completed);
    if (completion.retireRetainedObservation) {
	std::fprintf(stderr, "FAIL: budget-owned retained observation retired\n");
	return 1;
    }
    completed.retainedRefinementBudgetBlocked = false;
    completed.changedCut = true;
    completion = policy.completePass(completed);
    if (completion.retireRetainedObservation) {
	std::fprintf(stderr, "FAIL: frame-owned retained observation retired\n");
	return 1;
    }

    policy.reset();
    input.population = population;
    input.scaleChanged = true;
    input.viewEpoch.set(11);
    input.currentProgressiveCeiling = 6;
    input.currentPointProxyPixelThreshold = 2.0f;
    policy.noteStableQuality(7, 1.5f, population,
	BObolLodViewEpoch(11), 7000);
    restore = policy.beginQuiet(input);
    if (!restore.apply || !restore.restoredProvenQuality ||
	!restore.startHandoff || restore.progressiveCeiling != 7 ||
	std::fabs(restore.pointProxyPixelThreshold - 1.5f) > 0.0001f ||
	restore.provenRenderCostFloor != 7000 ||
	policy.handoffCostFloor() != 7000) {
	std::fprintf(stderr, "FAIL: stable zoom-quality proof restore\n");
	return 1;
    }
    completed = Policy::CompletedPassInputs();
    completed.completed = true;
    completed.submissionPending = true;
    if (policy.completePass(completed).finishHandoff ||
	policy.handoffCostFloor() != 7000) {
	std::fprintf(stderr, "FAIL: pending submission released quality proof\n");
	return 1;
    }
    completed.submissionPending = false;
    completion = policy.completePass(completed);
    if (!completion.finishHandoff || policy.handoffCostFloor() != 0) {
	std::fprintf(stderr, "FAIL: quality proof handoff completion\n");
	return 1;
    }

    policy.reset();
    policy.noteStableQuality(7, 1.5f, population,
	BObolLodViewEpoch(11), 7000);
    input.viewEpoch.set(12);
    restore = policy.beginQuiet(input);
    if (restore.apply || restore.restoredProvenQuality ||
	policy.handoffCostFloor() != 0) {
	std::fprintf(stderr, "FAIL: stale view reused quality proof\n");
	return 1;
    }

    policy.reset();
    input.retainedMeshPayloads = false;
    restore = policy.beginQuiet(input);
    if (!restore.apply || !restore.clearPresentationLimits ||
	restore.progressiveCeiling != -1 ||
	std::fabs(restore.pointProxyPixelThreshold - 1.0f) > 0.0001f ||
	policy.handoffPending()) {
	std::fprintf(stderr, "FAIL: payload-free presentation reset\n");
	return 1;
    }
    return 0;
}

static int
test_publication_policy(void)
{
    using Policy = BObolLodPublicationPolicy;
    Policy policy;
    Policy::Inputs input;
    input.nowMicroseconds = 100;
    Policy::Decision decision = policy.decide(input);
    if (decision.keepPumpAlive || decision.requestFrame || policy.pending() ||
	policy.framePending() || policy.awaitingDeadline() ||
	policy.unpresentedCount() != 0) {
	std::fprintf(stderr, "FAIL: initial publication policy\n");
	return 1;
    }

    policy.noteApplied(2, 100);
    policy.noteApplied(3, 110);
    input.nowMicroseconds = 50099;
    decision = policy.decide(input);
    if (!decision.keepPumpAlive || decision.requestFrame ||
	!policy.pending() || policy.framePending() ||
	!policy.awaitingDeadline() ||
	policy.unpresentedCount() != 5 ||
	policy.firstUnpresentedMicroseconds() != 100) {
	std::fprintf(stderr, "FAIL: bounded publication batch\n");
	return 1;
    }
    input.nowMicroseconds = 50100;
    decision = policy.decide(input);
    if (!decision.keepPumpAlive || !decision.requestFrame ||
	!policy.framePending() || policy.awaitingDeadline()) {
	std::fprintf(stderr, "FAIL: stable publication deadline\n");
	return 1;
    }
    decision = policy.decide(input);
    if (!decision.keepPumpAlive || decision.requestFrame ||
	!policy.framePending()) {
	std::fprintf(stderr, "FAIL: duplicate publication frame\n");
	return 1;
    }
    policy.frameCompleted();
    if (policy.pending() || policy.framePending() ||
	policy.awaitingDeadline() ||
	policy.unpresentedCount() != 0 ||
	policy.firstUnpresentedMicroseconds() != 0) {
	std::fprintf(stderr, "FAIL: publication completion\n");
	return 1;
    }

    policy.noteApplied(1, 200);
    input = Policy::Inputs();
    input.nowMicroseconds = 200;
    input.firstUseful = true;
    decision = policy.decide(input);
    if (!decision.requestFrame || !policy.framePending()) {
	std::fprintf(stderr, "FAIL: first-useful publication\n");
	return 1;
    }
    policy.reset();
    policy.noteApplied(1, 300);
    input = Policy::Inputs();
    input.nowMicroseconds = 300;
    input.streamIdle = true;
    decision = policy.decide(input);
    if (!decision.requestFrame || !policy.framePending()) {
	std::fprintf(stderr, "FAIL: idle-stream publication\n");
	return 1;
    }

    policy.reset();
    policy.noteApplied(1, 1000);
    input = Policy::Inputs();
    input.interactive = true;
    input.nowMicroseconds = 34332;
    decision = policy.decide(input);
    if (decision.requestFrame) {
	std::fprintf(stderr, "FAIL: early interactive publication\n");
	return 1;
    }
    input.nowMicroseconds = 34333;
    decision = policy.decide(input);
    if (!decision.requestFrame) {
	std::fprintf(stderr, "FAIL: interactive publication deadline\n");
	return 1;
    }

    policy.reset();
    policy.noteApplied(1, 1000);
    input = Policy::Inputs();
    input.observedRenderNanoseconds = 100000000ULL;
    input.nowMicroseconds = 200999;
    decision = policy.decide(input);
    if (decision.requestFrame) {
	std::fprintf(stderr, "FAIL: adaptive publication fired early\n");
	return 1;
    }
    input.nowMicroseconds = 201000;
    decision = policy.decide(input);
    if (!decision.requestFrame) {
	std::fprintf(stderr, "FAIL: adaptive publication deadline\n");
	return 1;
    }

    policy.reset();
    policy.noteApplied(SIZE_MAX - 1, 500);
    policy.noteApplied(10, 600);
    input = Policy::Inputs();
    input.nowMicroseconds = 499;
    decision = policy.decide(input);
    if (!decision.requestFrame || policy.unpresentedCount() != SIZE_MAX) {
	std::fprintf(stderr, "FAIL: publication saturation/backward clock\n");
	return 1;
    }
    policy.reset();
    if (policy.pending() || policy.framePending()) {
	std::fprintf(stderr, "FAIL: publication reset\n");
	return 1;
    }
    return 0;
}

static int
test_compaction_policy(void)
{
    using Policy = BObolLodCompactionPolicy;
    Policy policy;
    const auto readyInput = []() {
	Policy::Inputs input;
	input.automatic = true;
	input.coverageRequired = true;
	input.coverageComplete = true;
	input.serviceAvailable = true;
	input.nowMicroseconds = 150;
	return input;
    };
    Policy::Decision decision = policy.decide(readyInput());
    if (decision.admission != Policy::Admission::INACTIVE ||
	decision.keepPumpAlive || decision.retiredRequest ||
	policy.pending() || policy.planning()) {
	std::fprintf(stderr, "FAIL: initial compaction policy\n");
	return 1;
    }

    policy.requestAfter(100, 50);
    if (!policy.pending() || policy.planning() ||
	policy.deadlineMicroseconds() != 150) {
	std::fprintf(stderr, "FAIL: delayed compaction request\n");
	return 1;
    }
    Policy::Inputs input = readyInput();
    input.coverageComplete = false;
    input.coverageProgressPending = true;
    input.nowMicroseconds = 200;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::DEFER ||
	!decision.keepPumpAlive || decision.retiredRequest ||
	!policy.pending()) {
	std::fprintf(stderr, "FAIL: active coverage not deferred\n");
	return 1;
    }
    input.coverageProgressPending = false;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::INACTIVE ||
	decision.keepPumpAlive || !decision.retiredRequest ||
	policy.pending() || policy.deadlineMicroseconds() != 0) {
	std::fprintf(stderr, "FAIL: unwitnessed coverage request retained\n");
	return 1;
    }

    /* With mesh LoD disabled there is no mesh-coverage prerequisite.  An
	 * empty demand pass must be allowed to retire service-initialization
	 * state instead of becoming an impossible background latch. */
    policy.requestAfter(100, 50);
    input = readyInput();
    input.coverageRequired = false;
    input.coverageComplete = false;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::PLAN ||
	decision.keepPumpAlive || decision.retiredRequest) {
	std::fprintf(stderr, "FAIL: optional coverage blocked compaction\n");
	return 1;
    }
    input = readyInput();
    input.nowMicroseconds = 149;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::DEFER ||
	!decision.keepPumpAlive) {
	std::fprintf(stderr, "FAIL: compaction bypassed quiet deadline\n");
	return 1;
    }
    input = readyInput();
    input.realizationPending = true;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::DEFER ||
	!decision.keepPumpAlive) {
	std::fprintf(stderr, "FAIL: compaction bypassed realization\n");
	return 1;
    }
    input = readyInput();
    input.submissionPending = true;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::DEFER ||
	!decision.keepPumpAlive) {
	std::fprintf(stderr, "FAIL: compaction bypassed submission\n");
	return 1;
    }
    input = readyInput();
    input.settlingPending = true;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::DEFER ||
	!decision.keepPumpAlive) {
	std::fprintf(stderr, "FAIL: compaction bypassed visual settling\n");
	return 1;
    }
    input = readyInput();
    input.automatic = false;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::INACTIVE ||
	decision.keepPumpAlive) {
	std::fprintf(stderr, "FAIL: disabled compaction kept pump alive\n");
	return 1;
    }
    input = readyInput();
    input.interactive = true;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::INACTIVE ||
	decision.keepPumpAlive) {
	std::fprintf(stderr, "FAIL: interactive compaction admitted\n");
	return 1;
    }
    input = readyInput();
    input.serviceAvailable = false;
    decision = policy.decide(input);
    if (decision.admission != Policy::Admission::INACTIVE ||
	decision.keepPumpAlive) {
	std::fprintf(stderr, "FAIL: service-free compaction admitted\n");
	return 1;
    }
    decision = policy.decide(readyInput());
    if (decision.admission != Policy::Admission::PLAN ||
	decision.keepPumpAlive || policy.continues(17)) {
	std::fprintf(stderr, "FAIL: ready compaction not admitted\n");
	return 1;
    }

    policy.finishPlanning(false, 17, 151);
    if (!policy.pending() || !policy.planning() ||
	!policy.continues(17) || policy.continues(18) ||
	policy.demandRevision() != 17 ||
	policy.deadlineMicroseconds() != 151) {
	std::fprintf(stderr, "FAIL: partial compaction continuation\n");
	return 1;
    }
    policy.finishPlanning(true, 17, 152);
    if (policy.pending() || policy.planning() ||
	policy.demandRevision() != 0 ||
	policy.deadlineMicroseconds() != 0) {
	std::fprintf(stderr, "FAIL: completed compaction retained state\n");
	return 1;
    }

    policy.requestAfter(200, 100);
    policy.requestImmediate(250);
    if (policy.deadlineMicroseconds() != 250) {
	std::fprintf(stderr, "FAIL: immediate compaction deadline\n");
	return 1;
    }
    policy.requestImmediate(275);
    if (policy.deadlineMicroseconds() != 250) {
	std::fprintf(stderr, "FAIL: immediate compaction delayed request\n");
	return 1;
    }
    policy.resetForServiceChange(false, 300, 50);
    if (policy.pending() || policy.planning() ||
	policy.deadlineMicroseconds() != 0) {
	std::fprintf(stderr, "FAIL: compaction service removal\n");
	return 1;
    }
    policy.resetForServiceChange(true, INT64_MAX - 10, 100);
    if (!policy.pending() || policy.planning() ||
	policy.deadlineMicroseconds() != INT64_MAX) {
	std::fprintf(stderr, "FAIL: compaction deadline saturation\n");
	return 1;
    }
    return 0;
}

/*
 * A deterministic composed-lifecycle model for the coordinator policies.
 *
 * The direct tests above prove each allocation-free value at its boundaries;
 * this runner proves that their witnesses compose.  It deliberately models
 * the owner-thread contracts rather than sleeping or launching workers: every
 * seed is an exact replay token and FakeClock is the only source of time.
 * Random perturbations include repeated partial coverage checkpoints,
 * provider failure during interaction, cancellation under memory pressure,
 * delayed publication, view invalidation, and bounded compaction.
 */
class BObolLodSeededSequence {
public:
    enum Witness : uint32_t {
	WITNESS_NONE = 0,
	WITNESS_SERVICE_TASK = 1u << 0,
	WITNESS_COVERAGE_CURSOR = 1u << 1,
	WITNESS_PUBLICATION_DEADLINE = 1u << 2,
	WITNESS_PRESENTATION_FRAME = 1u << 3,
	WITNESS_INTERACTION_INPUT = 1u << 4,
	WITNESS_QUIET_TIMER = 1u << 5,
	WITNESS_COMPACTION_DEADLINE = 1u << 6,
	WITNESS_COMPACTION_CURSOR = 1u << 7
    };

    explicit BObolLodSeededSequence(uint64_t seedValue) :
	seed(seedValue ? seedValue : 1), randomState(seedValue ? seedValue : 1)
    {
	coverage.reset();
	this->reconcile(LodMachine::Event::INITIALIZE);
    }

    uint32_t random(void)
    {
	/* xorshift64* has no platform library state and therefore preserves the
	 * seed as a stable failure reproducer. */
	uint64_t value = randomState;
	value ^= value >> 12;
	value ^= value << 25;
	value ^= value >> 27;
	randomState = value;
	return static_cast<uint32_t>(
	    (value * UINT64_C(2685821657736338717)) >> 32);
    }

    void perturb(unsigned int action)
    {
	switch (action % 13) {
	    case 0: this->startSource(); break;
	    case 1: this->coverageCheckpoint(); break;
	    case 2: this->completeProvider(); break;
	    case 3: this->startInteraction(); break;
	    case 4: this->endGesture(); break;
	    case 5: this->togglePressure(); break;
	    case 6: this->cancelGeneration(); break;
	    case 7: this->failProvider(); break;
	    case 8: this->applyResult(); break;
	    case 9: this->advanceClock(); break;
	    case 10: this->requestCompaction(); break;
	    case 11: this->invalidateView(); break;
	    case 12: this->changePolicy(); break;
	}
    }

    bool progressOne(void)
    {
	if (gestureActive) {
	    this->endGesture();
	    return true;
	}
	if (interactive || quietTimer) {
	    nowMicroseconds += 200000;
	    interactive = false;
	    quietTimer = false;
	    policyEpoch.advance();
	    this->reconcile(LodMachine::Event::POLICY_CHANGED);
	    return true;
	}
	if (servicePending) {
	    this->completeProvider();
	    return true;
	}
	if (coverageCursor) {
	    this->coverageCheckpoint();
	    return true;
	}
	if (publication.framePending()) {
	    publication.frameCompleted();
	    renderSerial++;
	    this->reconcile(LodMachine::Event::FRAME_COMPLETED);
	    return true;
	}
	if (publication.awaitingDeadline()) {
	    nowMicroseconds += 300000;
	    BObolLodPublicationPolicy::Inputs inputs;
	    inputs.nowMicroseconds = nowMicroseconds;
	    inputs.interactive = interactive;
	    inputs.streamIdle = !servicePending;
	    (void)publication.decide(inputs);
	    this->reconcile(LodMachine::Event::WORK_PUMPED);
	    return true;
	}
	if (compactionCursor) {
	    compaction.finishPlanning(true, residentDemandRevision,
		nowMicroseconds);
	    compactionCursor = false;
	    if (resources.recoveryPending())
		resources.markRecoveryHandled();
	    this->reconcile(LodMachine::Event::WORK_PUMPED);
	    return true;
	}
	if (compaction.pending() || resources.recoveryPending()) {
	    BObolLodCompactionPolicy::Inputs inputs =
		this->compactionInputs();
	    if (nowMicroseconds < compaction.deadlineMicroseconds())
		nowMicroseconds = compaction.deadlineMicroseconds();
	    inputs.nowMicroseconds = nowMicroseconds;
	    const BObolLodCompactionPolicy::Decision decision =
		compaction.decide(inputs);
	    if (decision.admission ==
		    BObolLodCompactionPolicy::Admission::PLAN) {
		residentDemandRevision++;
		if (!residentDemandRevision)
		    residentDemandRevision++;
		compaction.finishPlanning(false, residentDemandRevision,
		    nowMicroseconds);
		compactionCursor = true;
	    } else if (decision.retiredRequest &&
		resources.recoveryPending()) {
		resources.markRecoveryHandled();
	    } else if (decision.admission ==
		BObolLodCompactionPolicy::Admission::DEFER) {
		/* Its prerequisite must be one of the earlier named witnesses. */
		return false;
	    }
	    this->reconcile(LodMachine::Event::WORK_PUMPED);
	    return true;
	}
	return false;
    }

    uint32_t witnesses(void) const
    {
	uint32_t result = WITNESS_NONE;
	if (servicePending)
	    result |= WITNESS_SERVICE_TASK;
	if (coverageCursor)
	    result |= WITNESS_COVERAGE_CURSOR;
	if (publication.awaitingDeadline())
	    result |= WITNESS_PUBLICATION_DEADLINE;
	if (publication.framePending())
	    result |= WITNESS_PRESENTATION_FRAME;
	if (gestureActive)
	    result |= WITNESS_INTERACTION_INPUT;
	else if (interactive || quietTimer)
	    result |= WITNESS_QUIET_TIMER;
	if (compactionCursor)
	    result |= WITNESS_COMPACTION_CURSOR;
	else if (compaction.pending() || resources.recoveryPending())
	    result |= WITNESS_COMPACTION_DEADLINE;
	return result;
    }

    bool terminal(void)
    {
	const BObolLodConvergencePolicy::Decision decision =
	    this->convergenceDecision();
	if (decision.phase == BObolLodConvergencePolicy::Phase::ERROR)
	    return true;
	if (this->witnesses() != WITNESS_NONE || decision.visualPending)
	    return false;
	return machine.currentPhase() == LodMachine::Phase::STABLE ||
	    machine.currentPhase() == LodMachine::Phase::FALLBACK;
    }

    bool valid(void)
    {
	if (machine.lastInvariantMask() != LodMachine::INVARIANT_NONE)
	    return false;
	const BObolLodConvergencePolicy::Decision decision =
	    this->convergenceDecision();
	if (!this->terminal() && this->witnesses() == WITNESS_NONE)
	    return false;
	if (decision.viewReady && decision.visualPending)
	    return false;
	const uint32_t foregroundWitnesses = this->witnesses() &
	    ~(WITNESS_COMPACTION_DEADLINE | WITNESS_COMPACTION_CURSOR);
	/* Resident maintenance is explicitly background work: the visible
	 * coordinator may be STABLE while convergence reports BACKGROUND, but
	 * every foreground witness must keep it out of the terminal phase. */
	if (machine.currentPhase() == LodMachine::Phase::STABLE &&
	    foregroundWitnesses != WITNESS_NONE)
	    return false;
	return true;
    }

    uint64_t replaySeed(void) const { return seed; }
    LodMachine::Phase phase(void) const { return machine.currentPhase(); }
    uint32_t invariantMask(void) const
    {
	return machine.lastInvariantMask();
    }

private:
    void startSource(void)
    {
	if (sourceActive || servicePending)
	    return;
	providerFailed = false;
	providerCancelled = false;
	sourceActive = true;
	coverage.setRequired(true);
	servicePending = true;
	generationActive = true;
	visibleCount = 7 + random() % 57;
	coveredCount = 0;
	coverage.reset();
	coverage.activate(true);
	coverageCursor = true;
	viewEpoch.advance();
	this->reconcile(LodMachine::Event::WORK_SCHEDULED);
    }

    void coverageCheckpoint(void)
    {
	if (!coverageCursor)
	    return;
	const size_t remaining = visibleCount > coveredCount ?
	    visibleCount - coveredCount : 0;
	const size_t quantum = std::min<size_t>(
	    remaining, 1 + random() % 11);
	coverage.observe(quantum, quantum);
	coveredCount += quantum;
	if (quantum)
	    publication.noteApplied(quantum, nowMicroseconds);
	if (coveredCount >= visibleCount) {
	    (void)coverage.completeIfReady(true, false);
	    coverageCursor = false;
	}
	this->publicationDecision(!coveredCount || coveredCount == quantum,
	    !servicePending && !coverageCursor);
	this->reconcile(LodMachine::Event::RESULT_PUBLISHED);
    }

    void completeProvider(void)
    {
	if (!servicePending)
	    return;
	servicePending = false;
	this->publicationDecision(false, true);
	this->reconcile(LodMachine::Event::RESULT_PUBLISHED);
    }

    void startInteraction(void)
    {
	if (gestureActive)
	    return;
	interactive = true;
	gestureActive = true;
	quietTimer = false;
	this->reconcile(LodMachine::Event::INTERACTION_STARTED);
    }

    void endGesture(void)
    {
	if (!gestureActive)
	    return;
	gestureActive = false;
	interactive = true;
	quietTimer = true;
	this->reconcile(LodMachine::Event::INTERACTION_ENDED);
    }

    void togglePressure(void)
    {
	const bool cpu = !resources.cpuPressure();
	const bool gpu = cpu && ((random() & 1u) != 0);
	const BObolLodResourcePolicy::Decision decision =
	    resources.observe(cpu, gpu, true);
	if (decision.queueRecovery)
	    compaction.requestImmediate(nowMicroseconds);
	this->reconcile(LodMachine::Event::FRAME_COMPLETED);
    }

    void cancelGeneration(void)
    {
	if (!sourceActive && !servicePending)
	    return;
	servicePending = false;
	sourceActive = false;
	coverage.setRequired(false);
	generationActive = false;
	providerCancelled = true;
	providerFailed = false;
	coverageCursor = false;
	coverage.reset();
	publication.reset();
	this->reconcile(LodMachine::Event::GENERATION_CANCELLED);
    }

    void failProvider(void)
    {
	if (!servicePending)
	    return;
	servicePending = false;
	sourceActive = false;
	coverage.setRequired(false);
	generationActive = false;
	providerFailed = true;
	providerCancelled = false;
	coverageCursor = false;
	coverage.reset();
	publication.reset();
	this->reconcile(LodMachine::Event::RESULT_PUBLISHED);
    }

    void applyResult(void)
    {
	if (!sourceActive || providerFailed || providerCancelled)
	    return;
	publication.noteApplied(1 + random() % 8, nowMicroseconds);
	this->publicationDecision(false, !servicePending);
	this->reconcile(LodMachine::Event::RESULT_PUBLISHED);
    }

    void advanceClock(void)
    {
	nowMicroseconds += 1 + random() % 400000;
	if (quietTimer && nowMicroseconds > 150000) {
	    interactive = false;
	    quietTimer = false;
	    policyEpoch.advance();
	}
	this->publicationDecision(false, !servicePending);
	this->reconcile(LodMachine::Event::WORK_PUMPED);
    }

    void requestCompaction(void)
    {
	compaction.requestAfter(nowMicroseconds, 50000);
	this->reconcile(LodMachine::Event::WORK_SCHEDULED);
    }

    void invalidateView(void)
    {
	viewEpoch.advance();
	if (sourceActive) {
	    coverage.reset();
	    coverage.activate(true);
	    coverageCursor = true;
	    coveredCount = 0;
	}
	this->reconcile(LodMachine::Event::VIEW_INVALIDATED);
    }

    void changePolicy(void)
    {
	policyEpoch.advance();
	this->reconcile(LodMachine::Event::POLICY_CHANGED);
    }

    void publicationDecision(bool firstUseful, bool streamIdle)
    {
	BObolLodPublicationPolicy::Inputs inputs;
	inputs.nowMicroseconds = nowMicroseconds;
	inputs.interactive = interactive;
	inputs.firstUseful = firstUseful;
	inputs.streamIdle = streamIdle;
	(void)publication.decide(inputs);
    }

    BObolLodCompactionPolicy::Inputs compactionInputs(void) const
    {
	BObolLodCompactionPolicy::Inputs inputs;
	inputs.automatic = true;
	inputs.interactive = interactive;
	inputs.coverageRequired = sourceActive;
	inputs.coverageComplete =
	    !sourceActive || coverage.coverageComplete();
	inputs.coverageProgressPending = coverageCursor;
	inputs.settlingPending = publication.pending() || quietTimer;
	inputs.realizationPending = servicePending;
	inputs.submissionPending = coverageCursor;
	inputs.serviceAvailable = true;
	inputs.nowMicroseconds = nowMicroseconds;
	return inputs;
    }

    BObolLodConvergencePolicy::Decision convergenceDecision(void)
    {
	BObolLodConvergencePolicy::Inputs inputs;
	inputs.viewEpoch = viewEpoch;
	inputs.policyEpoch = policyEpoch;
	inputs.expectedLeafCount = sourceActive || providerFailed ?
	    visibleCount : 0;
	inputs.availableLeafCount = sourceActive ? coveredCount : 0;
	inputs.visibleTargetCount = sourceActive ? visibleCount : 0;
	inputs.activePayloadCount = sourceActive ? coveredCount : 0;
	inputs.satisfiedPayloadCount = inputs.activePayloadCount;
	inputs.pendingTasks = servicePending ? 1 : 0;
	inputs.submissionPending = coverageCursor;
	inputs.publicationPending = publication.pending();
	inputs.interactive = interactive;
	inputs.compactionPending =
	    compaction.pending() || compactionCursor ||
	    resources.recoveryPending();
	inputs.structuralDiscovery = servicePending;
	inputs.failedSourceCount = providerFailed ? 1 : 0;
	inputs.gpuMemoryPressure = resources.gpuPressure();
	return convergence.evaluate(inputs);
    }

    void reconcile(LodMachine::Event event)
    {
	LodMachine::Inputs inputs;
	inputs.interactive = interactive;
	inputs.gestureActive = gestureActive;
	inputs.coverageActive = coverage.effectiveActive();
	inputs.coverageComplete = coverage.effectiveComplete();
	inputs.compacting = (compaction.planning() || compactionCursor) &&
	    coverage.compactionAllowed();
	inputs.generationActive = generationActive;
	inputs.cpuMemoryPressure = resources.cpuPressure();
	inputs.gpuMemoryPressure = resources.gpuPressure();
	inputs.resourceRecoveryPending = resources.recoveryPending();
	inputs.settlingWork = servicePending || quietTimer ||
	    publication.pending();

	LodMachine::Witness witness;
	witness.viewEpoch = viewEpoch.value();
	witness.policyEpoch = policyEpoch.value();
	witness.renderSerial = renderSerial;
	witness.activeGeneration = generationActive ? 1 : 0;
	witness.residentDemandRevision = residentDemandRevision;
	witness.resourcePressureRevision = resources.revision();
	witness.visibleCount = visibleCount;
	witness.completedCount = coveredCount;
	witness.pendingCount =
	    static_cast<size_t>(servicePending) +
	    static_cast<size_t>(coverageCursor) +
	    static_cast<size_t>(publication.pending()) +
	    static_cast<size_t>(quietTimer);
	(void)machine.dispatch(event, inputs, witness);
    }

    const uint64_t seed;
    uint64_t randomState;
    int64_t nowMicroseconds = 1;
    BObolLodViewEpoch viewEpoch {1};
    BObolLodPolicyEpoch policyEpoch {1};
    uint64_t renderSerial = 0;
    uint64_t residentDemandRevision = 0;
    size_t visibleCount = 0;
    size_t coveredCount = 0;
    bool sourceActive = false;
    bool servicePending = false;
    bool generationActive = false;
    bool coverageCursor = false;
    bool providerFailed = false;
    bool providerCancelled = false;
    bool interactive = false;
    bool gestureActive = false;
    bool quietTimer = false;
    bool compactionCursor = false;
    LodMachine machine;
    BObolLodResourcePolicy resources;
    BObolLodCoveragePolicy coverage;
    BObolLodPublicationPolicy publication;
    BObolLodCompactionPolicy compaction;
    BObolLodConvergencePolicy convergence;
};

static int
test_seeded_composed_lifecycle(void)
{
    static const size_t seedCount = 512;
    static const size_t perturbationsPerSeed = 96;
    static const size_t dischargeBound = 128;
    const auto discharge = [](BObolLodSeededSequence &sequence,
	    const char *scenario) {
	for (size_t step = 0;
		step < dischargeBound && !sequence.terminal(); ++step) {
	    if (sequence.progressOne() && sequence.valid())
		continue;
	    std::fprintf(stderr,
		"FAIL: composed lifecycle %s seed=%llu step=%zu "
		"phase=%u witnesses=0x%x invariant=0x%x\n",
		scenario,
		static_cast<unsigned long long>(sequence.replaySeed()),
		step, static_cast<unsigned int>(sequence.phase()),
		sequence.witnesses(), sequence.invariantMask());
	    return 1;
	}
	if (sequence.terminal())
	    return 0;
	std::fprintf(stderr,
	    "FAIL: composed lifecycle %s exceeded discharge bound seed=%llu "
	    "phase=%u witnesses=0x%x\n", scenario,
	    static_cast<unsigned long long>(sequence.replaySeed()),
	    static_cast<unsigned int>(sequence.phase()),
	    sequence.witnesses());
	return 1;
    };

    /* Keep the release-critical compositions explicit even though the seeded
     * corpus also reaches them.  The action numbers are part of this fake
     * service's stable replay vocabulary. */
    {
	BObolLodSeededSequence repeatedCheckpoints(UINT64_C(0x43504f494e54));
	repeatedCheckpoints.perturb(0);  /* source/service start */
	repeatedCheckpoints.perturb(1);
	repeatedCheckpoints.perturb(1);
	repeatedCheckpoints.perturb(1);
	repeatedCheckpoints.perturb(2);  /* provider stream idle */
	if (!repeatedCheckpoints.valid() ||
	    discharge(repeatedCheckpoints, "repeated-checkpoints"))
	    return 1;
    }
    {
	BObolLodSeededSequence failureDuringInput(UINT64_C(0x4641494c5552));
	failureDuringInput.perturb(0); /* source/service start */
	failureDuringInput.perturb(3); /* interaction start */
	failureDuringInput.perturb(7); /* terminal provider failure */
	if (!failureDuringInput.valid() || !failureDuringInput.terminal()) {
	    std::fprintf(stderr,
		"FAIL: provider failure during interaction did not terminate\n");
	    return 1;
	}
    }
    {
	BObolLodSeededSequence pressureCancellation(UINT64_C(0x505245535355));
	pressureCancellation.perturb(0); /* source/service start */
	pressureCancellation.perturb(5); /* memory pressure edge */
	pressureCancellation.perturb(6); /* generation cancellation */
	if (!pressureCancellation.valid() ||
	    discharge(pressureCancellation, "pressure-during-cancellation"))
	    return 1;
    }

    for (uint64_t seed = 1; seed <= seedCount; ++seed) {
	BObolLodSeededSequence sequence(seed);
	for (size_t step = 0; step < perturbationsPerSeed; ++step) {
	    sequence.perturb(sequence.random());
	    const size_t opportunistic = sequence.random() % 3;
	    for (size_t i = 0; i < opportunistic; ++i)
		(void)sequence.progressOne();
	    if (!sequence.valid()) {
		std::fprintf(stderr,
		    "FAIL: composed lifecycle seed=%llu step=%zu "
		    "phase=%u witnesses=0x%x invariant=0x%x\n",
		    static_cast<unsigned long long>(sequence.replaySeed()),
		    step, static_cast<unsigned int>(sequence.phase()),
		    sequence.witnesses(), sequence.invariantMask());
		return 1;
	    }
	}
	if (discharge(sequence, "seeded"))
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
    if (test_exhaustive_state_contract())
	return 1;
    if (test_event_and_invariant_contract())
	return 1;
    if (test_resource_policy())
	return 1;
    if (test_headroom_policy())
	return 1;
    if (test_coverage_policy())
	return 1;
    if (test_convergence_policy())
	return 1;
    if (test_budget_policy())
	return 1;
    if (test_quality_policy())
	return 1;
    if (test_view_demand_policy())
	return 1;
    if (test_presentation_policy())
	return 1;
    if (test_publication_policy())
	return 1;
    if (test_compaction_policy())
	return 1;
    if (test_seeded_composed_lifecycle())
	return 1;
    return 0;
}
