---------------------------- MODULE LodCoordinator ----------------------------
EXTENDS LodCommon

CONSTANT InjectPrematureStable

ASSUME InjectPrematureStable \in BOOLEAN

VARIABLES
    step,
    phase,
    lastEvent,
    observedInputs,
    inputs,
    observedWitness,
    witness,
    invariantMask,
    invariantHistory,
    violationCount,
    dispatchCount,
    transitionSerial,
    progressSerial,
    stagnantCount,
    phaseWitness

vars == <<
    step, phase, lastEvent, observedInputs, inputs, observedWitness, witness,
    invariantMask, invariantHistory, violationCount, dispatchCount,
    transitionSerial, progressSerial, stagnantCount, phaseWitness
>>

NoInputs == [
    interactive |-> FALSE,
    compacting |-> FALSE,
    coverageActive |-> FALSE,
    coverageComplete |-> FALSE,
    generationActive |-> FALSE,
    settlingWork |-> FALSE,
    gestureActive |-> FALSE,
    cpuMemoryPressure |-> FALSE,
    gpuMemoryPressure |-> FALSE,
    resourceRecoveryPending |-> FALSE
]

NoWitness == [
    sequence |-> 0,
    viewEpoch |-> 0,
    policyEpoch |-> 0,
    renderSerial |-> 0,
    activeGeneration |-> 0,
    residentDemandRevision |-> 0,
    resourcePressureRevision |-> 0,
    visibleCount |-> 0,
    completedCount |-> 0,
    pendingCount |-> 0
]

InputDomain == {
    [ interactive |-> i,
      compacting |-> c,
      coverageActive |-> ca,
      coverageComplete |-> cc,
      generationActive |-> g,
      settlingWork |-> s,
      gestureActive |-> ga,
      cpuMemoryPressure |-> cp,
      gpuMemoryPressure |-> gp,
      resourceRecoveryPending |-> rp ] :
        i \in BOOLEAN, c \in BOOLEAN, ca \in BOOLEAN, cc \in BOOLEAN,
        g \in BOOLEAN, s \in BOOLEAN, ga \in BOOLEAN, cp \in BOOLEAN,
        gp \in BOOLEAN, rp \in BOOLEAN
}

\* These are the three boundary witnesses used by the exhaustive C++ test:
\* empty, completed greater than visible, and partially complete with work.
WitnessCases == {
    [NoWitness EXCEPT
        !.viewEpoch = 11, !.policyEpoch = 17,
        !.activeGeneration = 23],
    [NoWitness EXCEPT
        !.viewEpoch = 11, !.policyEpoch = 17,
        !.activeGeneration = 23,
        !.visibleCount = 5, !.completedCount = 7],
    [NoWitness EXCEPT
        !.viewEpoch = 11, !.policyEpoch = 17,
        !.activeGeneration = 23,
        !.visibleCount = 5, !.completedCount = 3, !.pendingCount = 2]
}

StableInputs == [NoInputs EXCEPT
    !.coverageComplete = TRUE,
    !.generationActive = TRUE
]

StableWitness == [NoWitness EXCEPT
    !.viewEpoch = 12,
    !.policyEpoch = 18,
    !.activeGeneration = 23,
    !.visibleCount = 5,
    !.completedCount = 5
]

InputType(i) ==
    i.interactive \in BOOLEAN
    /\ i.compacting \in BOOLEAN
    /\ i.coverageActive \in BOOLEAN
    /\ i.coverageComplete \in BOOLEAN
    /\ i.generationActive \in BOOLEAN
    /\ i.settlingWork \in BOOLEAN
    /\ i.gestureActive \in BOOLEAN
    /\ i.cpuMemoryPressure \in BOOLEAN
    /\ i.gpuMemoryPressure \in BOOLEAN
    /\ i.resourceRecoveryPending \in BOOLEAN

WitnessType(w) ==
    /\ w.sequence \in Nat
    /\ w.viewEpoch \in Nat
    /\ w.policyEpoch \in Nat
    /\ w.renderSerial \in Nat
    /\ w.activeGeneration \in Nat
    /\ w.residentDemandRevision \in Nat
    /\ w.resourcePressureRevision \in Nat
    /\ w.visibleCount \in Nat
    /\ w.completedCount \in Nat
    /\ w.pendingCount \in Nat

PhaseFor(i, w) ==
    IF i.interactive \/ i.gestureActive THEN Interactive
    ELSE IF i.coverageActive \/ ~i.coverageComplete
        THEN IF i.generationActive THEN Coverage ELSE Fallback
    ELSE IF i.compacting \/ i.resourceRecoveryPending THEN Compacting
    ELSE IF i.settlingWork \/ w.pendingCount # 0 THEN Settling
    ELSE Stable

CanonicalInputs(e, observed, w) ==
    LET edge ==
        CASE e = InteractionStarted ->
                [observed EXCEPT
                    !.interactive = TRUE,
                    !.gestureActive = TRUE]
          [] e = InteractionEnded ->
                [observed EXCEPT !.gestureActive = FALSE]
          [] e = GenerationCancelled ->
                [observed EXCEPT !.generationActive = FALSE]
          [] OTHER -> observed
        coverageSafe ==
            IF ~edge.coverageComplete
            THEN [edge EXCEPT !.compacting = FALSE]
            ELSE edge
        pressureSafe ==
            IF ~coverageSafe.cpuMemoryPressure /\
               ~coverageSafe.gpuMemoryPressure
            THEN [coverageSafe EXCEPT !.resourceRecoveryPending = FALSE]
            ELSE coverageSafe
    IN IF w.pendingCount # 0
       THEN [pressureSafe EXCEPT !.settlingWork = TRUE]
       ELSE pressureSafe

CanonicalWitness(observed) ==
    [observed EXCEPT
        !.completedCount = MinNat(@, observed.visibleCount)]

InvariantSet(e, i, w) ==
    (IF w.completedCount > w.visibleCount
        THEN {"COMPLETED_EXCEEDS_VISIBLE"} ELSE {})
    \cup (IF i.compacting /\ ~i.coverageComplete
        THEN {"COMPACTION_WITHOUT_COVERAGE"} ELSE {})
    \cup (IF PhaseFor(i, w) = Stable /\ ~i.coverageComplete
        THEN {"STABLE_WITHOUT_COVERAGE"} ELSE {})
    \cup (IF ~i.interactive /\ ~i.gestureActive /\ ~i.compacting
              /\ i.coverageComplete /\ ~i.coverageActive
              /\ ~i.settlingWork /\ w.pendingCount # 0
        THEN {"STABLE_WITH_PENDING_WORK"} ELSE {})
    \cup (IF e = InteractionStarted
              /\ (~i.interactive \/ ~i.gestureActive)
        THEN {"INTERACTION_EVENT_WITHOUT_INTERACTION"} ELSE {})
    \cup (IF e = InteractionEnded /\ i.gestureActive
        THEN {"INTERACTION_END_WITH_GESTURE"} ELSE {})
    \cup (IF e = GenerationCancelled /\ i.generationActive
        THEN {"CANCEL_EVENT_WITH_GENERATION"} ELSE {})
    \cup (IF i.resourceRecoveryPending
              /\ ~i.cpuMemoryPressure /\ ~i.gpuMemoryPressure
        THEN {"RESOURCE_RECOVERY_WITHOUT_PRESSURE"} ELSE {})

SameProgress(a, b) ==
    /\ a.viewEpoch = b.viewEpoch
    /\ a.policyEpoch = b.policyEpoch
    /\ a.renderSerial = b.renderSerial
    /\ a.activeGeneration = b.activeGeneration
    /\ a.residentDemandRevision = b.residentDemandRevision
    /\ a.resourcePressureRevision = b.resourcePressureRevision
    /\ a.visibleCount = b.visibleCount
    /\ a.completedCount = b.completedCount
    /\ a.pendingCount = b.pendingCount

\* Used only by the expected-failure mutation test.  Normal configurations
\* select PhaseFor exactly.
SelectedPhase(i, w) ==
    IF InjectPrematureStable
       /\ i.coverageComplete
       /\ ~i.interactive
       /\ ~i.gestureActive
       /\ ~i.compacting
       /\ ~i.resourceRecoveryPending
       /\ w.pendingCount # 0
    THEN Stable
    ELSE PhaseFor(i, w)

ApplyDispatch(e, observed, observedW, nextStep) ==
    LET canonical == CanonicalInputs(e, observed, observedW)
        canonicalW == CanonicalWitness(observedW)
        nextPhase == SelectedPhase(canonical, canonicalW)
        violations == InvariantSet(e, observed, observedW)
        phaseChanged == phase # nextPhase
        recordProgress ==
            phaseWitness[nextPhase].sequence = 0
            \/ ~SameProgress(phaseWitness[nextPhase], canonicalW)
        nextProgressSerial ==
            IF recordProgress THEN progressSerial + 1 ELSE progressSerial
        nextWitness ==
            [canonicalW EXCEPT !.sequence = nextProgressSerial]
    IN
        /\ step' = nextStep
        /\ phase' = nextPhase
        /\ lastEvent' = e
        /\ observedInputs' = observed
        /\ inputs' = canonical
        /\ observedWitness' = observedW
        /\ witness' = canonicalW
        /\ invariantMask' = violations
        /\ invariantHistory' = invariantHistory \cup violations
        /\ violationCount' = violationCount +
                (IF violations = {} THEN 0 ELSE 1)
        /\ dispatchCount' = dispatchCount + 1
        /\ transitionSerial' = transitionSerial +
                (IF phaseChanged THEN 1 ELSE 0)
        /\ progressSerial' = nextProgressSerial
        /\ stagnantCount' =
                IF phaseChanged \/ recordProgress THEN 0
                ELSE stagnantCount + 1
        /\ phaseWitness' =
                IF recordProgress
                THEN [phaseWitness EXCEPT ![nextPhase] = nextWitness]
                ELSE phaseWitness

Init ==
    /\ step = 0
    /\ phase = Fallback
    /\ lastEvent = Initialize
    /\ observedInputs = NoInputs
    /\ inputs = NoInputs
    /\ observedWitness = NoWitness
    /\ witness = NoWitness
    /\ invariantMask = {}
    /\ invariantHistory = {}
    /\ violationCount = 0
    /\ dispatchCount = 0
    /\ transitionSerial = 0
    /\ progressSerial = 0
    /\ stagnantCount = 0
    /\ phaseWitness = [p \in Phases |-> NoWitness]

DispatchObserved ==
    /\ step = 0
    /\ \E e \in Events, i \in InputDomain, w \in WitnessCases:
        ApplyDispatch(e, i, w, 1)

RecoverQuietly ==
    /\ step = 1
    /\ ApplyDispatch(ViewObserved, StableInputs, StableWitness, 2)

Done ==
    /\ step = 2
    /\ UNCHANGED vars

Next == DispatchObserved \/ RecoverQuietly \/ Done

Spec == Init /\ [][Next]_vars

TypeOK ==
    /\ step \in 0..2
    /\ phase \in Phases
    /\ lastEvent \in Events
    /\ InputType(observedInputs)
    /\ InputType(inputs)
    /\ WitnessType(observedWitness)
    /\ WitnessType(witness)
    /\ invariantMask \subseteq InvariantSet(lastEvent,
            observedInputs, observedWitness)
    /\ invariantHistory \subseteq {
        "COMPLETED_EXCEEDS_VISIBLE", "COMPACTION_WITHOUT_COVERAGE",
        "STABLE_WITHOUT_COVERAGE", "STABLE_WITH_PENDING_WORK",
        "INTERACTION_EVENT_WITHOUT_INTERACTION",
        "INTERACTION_END_WITH_GESTURE", "CANCEL_EVENT_WITH_GENERATION",
        "RESOURCE_RECOVERY_WITHOUT_PRESSURE"
       }
    /\ violationCount \in Nat
    /\ dispatchCount \in Nat
    /\ transitionSerial \in Nat
    /\ progressSerial \in Nat
    /\ stagnantCount \in Nat
    /\ \A p \in Phases: WitnessType(phaseWitness[p])

CanonicalState ==
    /\ inputs = CanonicalInputs(lastEvent, observedInputs, observedWitness)
    /\ witness = CanonicalWitness(observedWitness)
    /\ witness.completedCount <= witness.visibleCount
    /\ (~inputs.coverageComplete => ~inputs.compacting)
    /\ (~inputs.cpuMemoryPressure /\ ~inputs.gpuMemoryPressure =>
            ~inputs.resourceRecoveryPending)
    /\ (witness.pendingCount # 0 => inputs.settlingWork)

ExactPhaseDecision == phase = SelectedPhase(inputs, witness)

InvariantDiagnosticsExact ==
    invariantMask = InvariantSet(lastEvent, observedInputs, observedWitness)

CurrentPhaseWitnessTracksProgress ==
    step = 0 \/
        (phaseWitness[phase].sequence > 0
         /\ SameProgress(phaseWitness[phase], witness))

StableHasNoForeground ==
    phase = Stable =>
        inputs.coverageComplete
        /\ ~inputs.interactive
        /\ ~inputs.gestureActive
        /\ ~inputs.compacting
        /\ ~inputs.resourceRecoveryPending
        /\ ~inputs.settlingWork
        /\ witness.pendingCount = 0

ContradictoryInputDischarges ==
    step = 2 =>
        phase = Stable
        /\ invariantMask = {}
        /\ inputs.coverageComplete
        /\ witness.completedCount = witness.visibleCount

=============================================================================
