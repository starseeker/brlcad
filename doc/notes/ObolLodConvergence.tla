------------------------- MODULE ObolLodConvergence -------------------------
\* Control-plane model for the two extreme progressive-display workloads.
\*
\* "lucy" represents one logical occurrence with several cumulative resident
\* suffix and quality transitions.  "scene" represents many occurrences with
\* a long coverage census and a comparatively short terminal quality choice.
\* Small counters are sufficient: TLC is checking the protocol, not geometry
\* cardinality, allocation scoring, or wall-clock performance.
\* `physicalDebtRemaining` is known pixel demand; `residencyRemaining` is only
\* the suffix admitted by the current allocation.  A constrained terminal
\* state may retain the former, but never the latter.
\*
\* The contract under test is that camera input may interrupt any bounded
\* pass, but a finite final view epoch never restarts coverage or quality work
\* without a new input/source epoch.  Every submission level has an enabled
\* owner-thread action, and a finite workload eventually reaches either its
\* requested quality or an explicit performance-limited stable state.

EXTENDS Naturals, TLC

CONSTANT MaxInputEpoch

Profiles == {"lucy", "scene"}
WorkerStates == {"idle", "queued", "inflight", "result"}
ConstraintWitnesses == {"none", "deadline", "memory", "presentation"}

ProfileCoverage(p) == IF p = "lucy" THEN 1 ELSE 3
ProfileResidency(p) == IF p = "lucy" THEN 2 ELSE 1
ProfileQuality(p) == IF p = "lucy" THEN 3 ELSE 1
ProfilePhysicalDebt(p) == IF p = "lucy" THEN 3 ELSE 2

VARIABLES profile,
          inputEpoch,
          inputOpen,
          coverageRemaining,
          residencyRemaining,
          physicalDebtRemaining,
          qualityRemaining,
          scanActive,
          workerState,
          allocationPending,
          framePending,
          handoffPending,
          submissionPending,
          performanceLimited,
          constraintWitness

vars == <<profile, inputEpoch, inputOpen, coverageRemaining,
          residencyRemaining, physicalDebtRemaining, qualityRemaining,
          scanActive, workerState,
          allocationPending, framePending, handoffPending,
          submissionPending, performanceLimited, constraintWitness>>

TypeOK ==
    /\ profile \in Profiles
    /\ inputEpoch \in 0..MaxInputEpoch
    /\ inputOpen \in BOOLEAN
    /\ coverageRemaining \in 0..3
    /\ residencyRemaining \in 0..2
    /\ physicalDebtRemaining \in 0..3
    /\ qualityRemaining \in 0..3
    /\ scanActive \in BOOLEAN
    /\ workerState \in WorkerStates
    /\ allocationPending \in BOOLEAN
    /\ framePending \in BOOLEAN
    /\ handoffPending \in BOOLEAN
    /\ submissionPending \in BOOLEAN
    /\ performanceLimited \in BOOLEAN
    /\ constraintWitness \in ConstraintWitnesses

Init ==
    /\ profile \in Profiles
    /\ inputEpoch = 0
    /\ inputOpen = FALSE
    /\ coverageRemaining = ProfileCoverage(profile)
    /\ residencyRemaining = ProfileResidency(profile)
    /\ physicalDebtRemaining = ProfilePhysicalDebt(profile)
    /\ qualityRemaining = ProfileQuality(profile)
    /\ scanActive = FALSE
    /\ workerState = "idle"
    /\ allocationPending = FALSE
    /\ framePending = FALSE
    /\ handoffPending = FALSE
    /\ submissionPending = TRUE
    /\ performanceLimited = FALSE
    /\ constraintWitness = "none"

\* A camera event starts a new bounded view epoch.  Resident mesh suffixes are
\* cumulative and survive, while the visibility census and quality demand are
\* recomputed.  An obsolete scan, allocation, or frame is invalidated; an
\* already-running cumulative loader remains useful and owns the next wakeup.
BeginInput ==
    /\ ~inputOpen
    /\ inputEpoch < MaxInputEpoch
    /\ inputEpoch' = inputEpoch + 1
    /\ inputOpen' = TRUE
    /\ coverageRemaining' = ProfileCoverage(profile)
    /\ UNCHANGED <<residencyRemaining, physicalDebtRemaining>>
    /\ qualityRemaining' = ProfileQuality(profile)
    /\ scanActive' = FALSE
    /\ UNCHANGED workerState
    /\ allocationPending' = FALSE
    /\ framePending' = FALSE
    /\ handoffPending' = FALSE
    /\ submissionPending' = (workerState = "idle")
    /\ performanceLimited' = FALSE
    /\ constraintWitness' = "none"
    /\ UNCHANGED profile

EndInput ==
    /\ inputOpen
    /\ inputOpen' = FALSE
    /\ submissionPending' =
           (workerState = "idle") /\
           ~scanActive /\
           ~allocationPending /\
           ~framePending
    /\ UNCHANGED <<profile, inputEpoch, coverageRemaining,
                    residencyRemaining, physicalDebtRemaining,
                    qualityRemaining, scanActive,
                    workerState, allocationPending, framePending,
                    handoffPending, performanceLimited,
                    constraintWitness>>

StartCoverageScan ==
    /\ submissionPending
    /\ coverageRemaining > 0
    /\ workerState = "idle"
    /\ ~allocationPending
    /\ ~framePending
    /\ scanActive' = TRUE
    /\ submissionPending' = FALSE
    /\ UNCHANGED <<profile, inputEpoch, inputOpen, coverageRemaining,
                    residencyRemaining, physicalDebtRemaining,
                    qualityRemaining, workerState,
                    allocationPending, framePending, handoffPending,
                    performanceLimited, constraintWitness>>

AdvanceCoverageScan ==
    /\ scanActive
    /\ coverageRemaining > 0
    /\ coverageRemaining' = coverageRemaining - 1
    /\ scanActive' = (coverageRemaining > 1)
    /\ framePending' = (coverageRemaining = 1)
    /\ handoffPending' = (handoffPending \/ (coverageRemaining = 1))
    /\ UNCHANGED <<profile, inputEpoch, inputOpen, residencyRemaining,
                    physicalDebtRemaining, qualityRemaining, workerState,
                    allocationPending,
                    submissionPending, performanceLimited,
                    constraintWitness>>

QueueResidentLoad ==
    /\ submissionPending
    /\ coverageRemaining = 0
    /\ residencyRemaining > 0
    /\ ~scanActive
    /\ workerState = "idle"
    /\ ~allocationPending
    /\ ~framePending
    /\ workerState' = "queued"
    /\ submissionPending' = FALSE
    /\ UNCHANGED <<profile, inputEpoch, inputOpen, coverageRemaining,
                    residencyRemaining, physicalDebtRemaining,
                    qualityRemaining, scanActive,
                    allocationPending, framePending, handoffPending,
                    performanceLimited, constraintWitness>>

StartResidentLoad ==
    /\ workerState = "queued"
    /\ workerState' = "inflight"
    /\ UNCHANGED <<profile, inputEpoch, inputOpen, coverageRemaining,
                    residencyRemaining, physicalDebtRemaining,
                    qualityRemaining, scanActive,
                    allocationPending, framePending, handoffPending,
                    submissionPending, performanceLimited,
                    constraintWitness>>

FinishResidentLoad ==
    /\ workerState = "inflight"
    /\ workerState' = "result"
    /\ UNCHANGED <<profile, inputEpoch, inputOpen, coverageRemaining,
                    residencyRemaining, physicalDebtRemaining,
                    qualityRemaining, scanActive,
                    allocationPending, framePending, handoffPending,
                    submissionPending, performanceLimited,
                    constraintWitness>>

ApplyResidentResult ==
    /\ workerState = "result"
    /\ residencyRemaining > 0
    /\ residencyRemaining' = residencyRemaining - 1
    /\ physicalDebtRemaining' =
           (IF physicalDebtRemaining > 0
            THEN physicalDebtRemaining - 1
            ELSE 0)
    /\ workerState' = "idle"
    /\ framePending' = TRUE
    /\ handoffPending' = TRUE
    /\ submissionPending' = FALSE
    /\ UNCHANGED <<profile, inputEpoch, inputOpen, coverageRemaining,
                    qualityRemaining, scanActive, allocationPending,
                    performanceLimited, constraintWitness>>

StartAllocation ==
    /\ submissionPending
    /\ ~inputOpen
    /\ coverageRemaining = 0
    /\ residencyRemaining = 0
    /\ qualityRemaining > 0
    /\ ~scanActive
    /\ workerState = "idle"
    /\ ~allocationPending
    /\ ~framePending
    /\ allocationPending' = TRUE
    /\ submissionPending' = FALSE
    /\ UNCHANGED <<profile, inputEpoch, inputOpen, coverageRemaining,
                    residencyRemaining, physicalDebtRemaining,
                    qualityRemaining, scanActive,
                    workerState, framePending, handoffPending,
                    performanceLimited, constraintWitness>>

\* The allocator either admits one richer coherent population or proves that
\* the current one is the richest sustainable population.  Both outcomes are
\* terminal progress; repeating the same allocation is not an outcome.
FinishAllocation ==
    /\ allocationPending
    /\ qualityRemaining > 0
    /\ allocationPending' = FALSE
    /\ framePending' = TRUE
    /\ handoffPending' = TRUE
    /\ submissionPending' = FALSE
    /\ \/ /\ qualityRemaining' = qualityRemaining - 1
           /\ physicalDebtRemaining' =
                  (IF physicalDebtRemaining > 0
                   THEN physicalDebtRemaining - 1
                   ELSE 0)
           /\ UNCHANGED <<performanceLimited, constraintWitness>>
       \/ /\ qualityRemaining' = 0
           /\ UNCHANGED physicalDebtRemaining
           /\ performanceLimited' = TRUE
           /\ constraintWitness' \in
                  (ConstraintWitnesses \ {"none"})
    /\ UNCHANGED <<profile, inputEpoch, inputOpen, coverageRemaining,
                    residencyRemaining, scanActive, workerState>>

MoreOwnerWork ==
    (coverageRemaining > 0) \/
    (residencyRemaining > 0) \/
    (~inputOpen /\ qualityRemaining > 0)

\* The implementation exposes several HUD names, but they must reduce to one
\* controller state.  Keep the mode derived from the work ownership state:
\* a second mutable mode flag would permit the historical "balancing" /
\* "refining" report to diverge from the work that can actually run.
\*
\* DISCOVERING has priority over settling because source coverage is the
\* user-visible first obligation.  Interactive input deliberately has highest
\* priority: retained prefixes remain useful, but no quiet-view quality claim
\* may be made while a camera epoch is changing.
ViewMode ==
    IF inputOpen THEN "interactive"
    ELSE IF coverageRemaining > 0 THEN "discovering"
    ELSE IF residencyRemaining > 0 \/ qualityRemaining > 0 \/
            scanActive \/ workerState # "idle" \/ allocationPending \/
            framePending \/ handoffPending \/ submissionPending
         THEN "settling"
    ELSE IF performanceLimited THEN "constrained"
    ELSE "stable"

CompleteFrame ==
    /\ framePending
    /\ framePending' = FALSE
    /\ submissionPending' = MoreOwnerWork
    /\ handoffPending' = IF MoreOwnerWork THEN TRUE ELSE FALSE
    /\ UNCHANGED <<profile, inputEpoch, inputOpen, coverageRemaining,
                    residencyRemaining, physicalDebtRemaining,
                    qualityRemaining, scanActive,
                    workerState, allocationPending, performanceLimited,
                    constraintWitness>>

\* A cursor can reach the end after the last bounded window or after a
\* performance-limited allocation.  Retiring that empty level is a real
\* owner-thread action; leaving it set is the 150k-style no-producer loop.
RetireEmptySubmission ==
    /\ submissionPending
    /\ ~inputOpen
    /\ coverageRemaining = 0
    /\ residencyRemaining = 0
    /\ qualityRemaining = 0
    /\ ~scanActive
    /\ workerState = "idle"
    /\ ~allocationPending
    /\ ~framePending
    /\ submissionPending' = FALSE
    /\ handoffPending' = FALSE
    /\ UNCHANGED <<profile, inputEpoch, inputOpen, coverageRemaining,
                    residencyRemaining, physicalDebtRemaining,
                    qualityRemaining, scanActive,
                    workerState, allocationPending, framePending,
                    performanceLimited, constraintWitness>>

Next ==
    \/ BeginInput
    \/ EndInput
    \/ StartCoverageScan
    \/ AdvanceCoverageScan
    \/ QueueResidentLoad
    \/ StartResidentLoad
    \/ FinishResidentLoad
    \/ ApplyResidentResult
    \/ StartAllocation
    \/ FinishAllocation
    \/ CompleteFrame
    \/ RetireEmptySubmission

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(BeginInput)
    /\ WF_vars(EndInput)
    /\ WF_vars(StartCoverageScan)
    /\ WF_vars(AdvanceCoverageScan)
    /\ WF_vars(QueueResidentLoad)
    /\ WF_vars(StartResidentLoad)
    /\ WF_vars(FinishResidentLoad)
    /\ WF_vars(ApplyResidentResult)
    /\ WF_vars(StartAllocation)
    /\ WF_vars(FinishAllocation)
    /\ WF_vars(CompleteFrame)
    /\ WF_vars(RetireEmptySubmission)

SubmissionCanProgress ==
    ENABLED StartCoverageScan \/
    ENABLED QueueResidentLoad \/
    ENABLED StartAllocation \/
    ENABLED RetireEmptySubmission

\* A submission level is a bounded cursor, not an abstract desire.  It must
\* always have a legal next owner-thread action in the current state.
PendingSubmissionHasWitness ==
    submissionPending => SubmissionCanProgress

\* Progress ownership is level-triggered but exclusive.  A transition may
\* leave more work desired; it must not assert a new submission cursor while
\* an existing scan, loader/result, allocation, or frame owns the next edge.
SingleProgressOwner ==
    /\ submissionPending =>
        /\ ~scanActive
        /\ workerState = "idle"
        /\ ~allocationPending
        /\ ~framePending
    /\ scanActive =>
        /\ ~submissionPending
        /\ workerState = "idle"
        /\ ~allocationPending
        /\ ~framePending
    /\ workerState # "idle" =>
        /\ ~submissionPending
        /\ ~scanActive
        /\ ~allocationPending
        /\ ~framePending
    /\ allocationPending =>
        /\ ~submissionPending
        /\ ~scanActive
        /\ workerState = "idle"
        /\ ~framePending
    /\ framePending =>
        /\ ~submissionPending
        /\ ~scanActive
        /\ workerState = "idle"
        /\ ~allocationPending

\* A constrained allocation is a real terminal quality result for its epoch.
\* Physical demand remains evidence for a future capacity/input edge, but it
\* cannot silently recreate residency or allocation work in the quiet epoch.
PerformanceLimitIsTerminal ==
    performanceLimited =>
        /\ residencyRemaining = 0
        /\ qualityRemaining = 0

\* Quality debt alone cannot justify a constrained terminal result.  The
\* controller must retain completed deadline, memory, or presentation
\* evidence for the unchanged view epoch.
PerformanceLimitHasWitness ==
    performanceLimited <=> constraintWitness # "none"

\* A constrained scene may retain known physical pixel debt.  That debt is an
\* input to a later capacity/view revision, not an enabled quiet-epoch loader.
\* Requiring physicalDebtRemaining = 0 here recreates the 150k residency loop.
ConstrainedDebtHasNoOwner ==
    (performanceLimited /\ physicalDebtRemaining > 0) =>
        /\ residencyRemaining = 0
        /\ qualityRemaining = 0

Stable ==
    /\ inputEpoch = MaxInputEpoch
    /\ ~inputOpen
    /\ coverageRemaining = 0
    /\ residencyRemaining = 0
    /\ qualityRemaining = 0
    /\ (physicalDebtRemaining = 0 \/ performanceLimited)
    /\ ~scanActive
    /\ workerState = "idle"
    /\ ~allocationPending
    /\ ~framePending
    /\ ~handoffPending
    /\ ~submissionPending

StableIsQuiescent ==
    Stable =>
        /\ ~submissionPending
        /\ workerState = "idle"
        /\ ~allocationPending
        /\ ~framePending
        /\ ~handoffPending

\* The named behavior modes are a total, disjoint presentation of the
\* controller state.  They give TLC an explicit check that a terminal HUD
\* cannot describe background work as visual refinement, nor describe an
\* unfinished coverage census as stable.
ModeIsWellDefined ==
    ViewMode \in {"interactive", "discovering", "settling", "constrained", "stable"}

StableModeIsQuiescent ==
    ViewMode = "stable" =>
        /\ ~inputOpen
        /\ coverageRemaining = 0
        /\ residencyRemaining = 0
        /\ qualityRemaining = 0
        /\ ~scanActive
        /\ workerState = "idle"
        /\ ~allocationPending
        /\ ~framePending
        /\ ~handoffPending
        /\ ~submissionPending

ConstrainedModeIsTerminal ==
    ViewMode = "constrained" =>
        /\ performanceLimited
        /\ ~inputOpen
        /\ coverageRemaining = 0
        /\ residencyRemaining = 0
        /\ qualityRemaining = 0

DiscoveringModeHasCoverageWork ==
    ViewMode = "discovering" => coverageRemaining > 0

\* MaxInputEpoch bounds user input for model checking.  Once that final input
\* ends, weak fairness of each named progress witness must drain the protocol.
EventuallyStableAfterFinalInput ==
    []((inputEpoch = MaxInputEpoch /\ ~inputOpen) => <>Stable)

EventuallyStable == <>Stable

=============================================================================
