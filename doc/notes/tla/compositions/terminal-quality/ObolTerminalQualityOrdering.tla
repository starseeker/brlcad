-------------------- MODULE ObolTerminalQualityOrdering --------------------
\* Composition contract for the structural frontier, terminal
\* renderer-ceiling, static-quality, and point-quality owners.
\*
\* A renderer-wide progressive ceiling is a transient presentation guard.
\* Point refinement changes the occurrence population hidden behind that
\* guard.  The point owner therefore cannot commit or measure a finer point
\* cut until the static owner has retired and the ceiling has been translated
\* into occurrence-local cuts.  When the occurrence-local minimum cannot
\* encode the completed safe framebuffer, that framebuffer and its ceiling
\* remain as an explicit constrained outcome.  Invalid completed
\* presentations consume a finite allowance and end at the same typed state.

EXTENDS Naturals, TLC

CONSTANTS PointStepLimit, PresentationAttemptLimit, StructuralRecordLimit

ASSUME PointStepLimit > 0
ASSUME PresentationAttemptLimit > 0
ASSUME StructuralRecordLimit > 0

Phases == {"select", "structural", "static", "reconcile", "point", "allocate",
           "capacity", "terminal"}
Outcomes == {"none", "ready", "constrained"}
ConstraintKinds == {"none", "static", "reconcile", "point", "capacity"}

VARIABLES phase,
          structuralFallbackCount,
          terminalFailureCount,
          terminalAggregateCount,
          ceilingEffective,
          staticOwner,
          pointStepsRemaining,
          initialPointSteps,
          capacityDemand,
          presentedAboveCertificate,
          allocationCurrent,
          initialCapacityDemand,
          allocationCanReplaceCeiling,
          attemptsRemaining,
          constraintKind,
          outcome

vars == <<phase, structuralFallbackCount, terminalFailureCount,
          terminalAggregateCount, ceilingEffective, staticOwner, pointStepsRemaining,
          initialPointSteps, capacityDemand, presentedAboveCertificate,
          allocationCurrent, initialCapacityDemand,
          allocationCanReplaceCeiling,
          attemptsRemaining, constraintKind, outcome>>

\* This refinement starts from an exact completed framebuffer.  The renderer's
\* broader projection histogram may still contain source records represented
\* by terminal subpixel aggregates; only uncollapsed fallbacks above the
\* explicit failure count constitute repair debt.
structuralFrontier == structuralFallbackCount > terminalFailureCount

TypeOK ==
    /\ phase \in Phases
    /\ structuralFallbackCount \in 0..StructuralRecordLimit
    /\ terminalFailureCount \in 0..structuralFallbackCount
    /\ terminalAggregateCount \in 0..StructuralRecordLimit
    /\ ceilingEffective \in BOOLEAN
    /\ staticOwner \in BOOLEAN
    /\ pointStepsRemaining \in 0..PointStepLimit
    /\ initialPointSteps \in 0..PointStepLimit
    /\ capacityDemand \in BOOLEAN
    /\ presentedAboveCertificate \in BOOLEAN
    /\ allocationCurrent \in BOOLEAN
    /\ initialCapacityDemand \in BOOLEAN
    /\ allocationCanReplaceCeiling \in BOOLEAN
    /\ attemptsRemaining \in 0..PresentationAttemptLimit
    /\ constraintKind \in ConstraintKinds
    /\ outcome \in Outcomes

Init ==
    /\ phase = "select"
    /\ structuralFallbackCount \in 0..StructuralRecordLimit
    /\ terminalFailureCount \in 0..structuralFallbackCount
    /\ terminalAggregateCount \in 0..StructuralRecordLimit
    /\ ceilingEffective \in BOOLEAN
    /\ staticOwner \in BOOLEAN
    /\ pointStepsRemaining \in 0..PointStepLimit
    /\ initialPointSteps = pointStepsRemaining
    /\ capacityDemand \in BOOLEAN
    /\ presentedAboveCertificate \in BOOLEAN
    /\ allocationCurrent \in BOOLEAN
    /\ initialCapacityDemand = capacityDemand
    /\ allocationCanReplaceCeiling \in BOOLEAN
    /\ attemptsRemaining = PresentationAttemptLimit
    /\ constraintKind = "none"
    /\ outcome = "none"

\* An exact structural fallback population is not an allocation-quality
\* problem.  It owns selection first and invalidates the occurrence allocation
\* when repair changes that population.  A stale retained-refinement
\* observation is deliberately absent from the guard: it is subordinate
\* evidence produced while the frontier was incomplete.
SelectStructural ==
    /\ phase = "select"
    /\ structuralFrontier
    /\ phase' = "structural"
    /\ UNCHANGED <<structuralFallbackCount, terminalFailureCount,
                    terminalAggregateCount, ceilingEffective, staticOwner,
                    pointStepsRemaining, initialPointSteps, capacityDemand,
                    presentedAboveCertificate, allocationCurrent,
                    initialCapacityDemand, allocationCanReplaceCeiling,
                    attemptsRemaining, constraintKind, outcome>>

CompleteStructural ==
    /\ phase = "structural"
    /\ phase' = "select"
    /\ structuralFallbackCount' = terminalFailureCount
    /\ allocationCurrent' = FALSE
    /\ UNCHANGED <<terminalFailureCount, terminalAggregateCount,
                    ceilingEffective, staticOwner, pointStepsRemaining,
                    initialPointSteps, capacityDemand,
                    presentedAboveCertificate, initialCapacityDemand,
                    allocationCanReplaceCeiling, attemptsRemaining,
                    constraintKind, outcome>>

\* Static quality may own a bounded candidate after the renderer-wide ceiling
\* has been translated.  The modeled ceiling is a multi-occurrence global
\* control, not the occurrence-local ordinal used by a one-object staircase.
SelectStatic ==
    /\ phase = "select"
    /\ ~structuralFrontier
    /\ ~ceilingEffective
    /\ staticOwner
    /\ pointStepsRemaining = 0
    /\ ~capacityDemand
    /\ ~presentedAboveCertificate
    /\ phase' = "static"
    /\ UNCHANGED <<ceilingEffective, staticOwner, pointStepsRemaining,
                    initialPointSteps, capacityDemand,
                    presentedAboveCertificate,
                    allocationCurrent,
                    initialCapacityDemand, allocationCanReplaceCeiling,
                    attemptsRemaining,
                    constraintKind, outcome>>

CompleteStatic ==
    /\ phase = "static"
    /\ phase' = "select"
    /\ staticOwner' = FALSE
    /\ UNCHANGED <<ceilingEffective, pointStepsRemaining,
                    initialPointSteps, capacityDemand,
                    presentedAboveCertificate,
                    allocationCurrent,
                    initialCapacityDemand, allocationCanReplaceCeiling,
                    attemptsRemaining,
                    constraintKind, outcome>>

\* Static quality may prove that no richer complete framebuffer fits its hard
\* deadline.  That is a typed terminal capacity result, not permission to
\* silently clear the owner and publish ready with unresolved quality debt.
ConstrainStatic ==
    /\ phase = "static"
    /\ phase' = "terminal"
    /\ staticOwner' = FALSE
    /\ attemptsRemaining' = 0
    /\ constraintKind' = "static"
    /\ outcome' = "constrained"
    /\ UNCHANGED <<ceilingEffective, pointStepsRemaining,
                    initialPointSteps, capacityDemand,
                    presentedAboveCertificate, allocationCurrent,
                    initialCapacityDemand,
                    allocationCanReplaceCeiling>>

\* An effective multi-occurrence renderer ceiling owns the sole successor.
\* Point and static demand remain recorded but cannot own a framebuffer.
SelectReconciliation ==
    /\ phase = "select"
    /\ ~structuralFrontier
    /\ ceilingEffective
    /\ phase' = "reconcile"
    /\ attemptsRemaining' = PresentationAttemptLimit
    /\ UNCHANGED <<ceilingEffective, staticOwner, pointStepsRemaining,
                    initialPointSteps, capacityDemand,
                    presentedAboveCertificate,
                    allocationCurrent,
                    initialCapacityDemand, allocationCanReplaceCeiling,
                    constraintKind, outcome>>

CompleteReconciliation ==
    /\ phase = "reconcile"
    /\ ceilingEffective
    /\ phase' = "select"
    /\ ceilingEffective' = FALSE
    /\ attemptsRemaining' = PresentationAttemptLimit
    /\ UNCHANGED <<staticOwner, pointStepsRemaining,
                    initialPointSteps, capacityDemand,
                    presentedAboveCertificate,
                    allocationCurrent,
                    initialCapacityDemand, allocationCanReplaceCeiling,
                    constraintKind, outcome>>

RetryReconciliation ==
    /\ phase = "reconcile"
    /\ attemptsRemaining > 1
    /\ attemptsRemaining' = attemptsRemaining - 1
    /\ UNCHANGED <<phase, ceilingEffective, staticOwner,
                    pointStepsRemaining, initialPointSteps, capacityDemand,
                    presentedAboveCertificate,
                    allocationCurrent,
                    initialCapacityDemand, allocationCanReplaceCeiling,
                    constraintKind, outcome>>

ConstrainReconciliation ==
    /\ phase = "reconcile"
    /\ attemptsRemaining = 1
    /\ ~allocationCanReplaceCeiling
    /\ phase' = "terminal"
    /\ staticOwner' = FALSE
    /\ attemptsRemaining' = 0
    /\ constraintKind' = "reconcile"
    /\ outcome' = "constrained"
    /\ UNCHANGED <<ceilingEffective, pointStepsRemaining,
                    initialPointSteps, capacityDemand,
                    presentedAboveCertificate,
                    allocationCurrent,
                    initialCapacityDemand, allocationCanReplaceCeiling>>

\* Point refinement is enabled only after the complete renderer control vector
\* is ceiling-free.  Each exact presentation consumes one finite bracket step.
SelectPoint ==
    /\ phase = "select"
    /\ ~structuralFrontier
    /\ ~ceilingEffective
    /\ pointStepsRemaining > 0
    /\ phase' = "point"
    /\ attemptsRemaining' = PresentationAttemptLimit
    /\ UNCHANGED <<ceilingEffective, staticOwner, pointStepsRemaining,
                    initialPointSteps, capacityDemand,
                    presentedAboveCertificate,
                    allocationCurrent,
                    initialCapacityDemand, allocationCanReplaceCeiling,
                    constraintKind, outcome>>

CompletePoint ==
    /\ phase = "point"
    /\ ~ceilingEffective
    /\ pointStepsRemaining > 0
    /\ phase' = "select"
    /\ pointStepsRemaining' = pointStepsRemaining - 1
    /\ allocationCurrent' = FALSE
    /\ attemptsRemaining' = PresentationAttemptLimit
    /\ UNCHANGED <<ceilingEffective, staticOwner, initialPointSteps,
                    capacityDemand, presentedAboveCertificate,
                    initialCapacityDemand, allocationCanReplaceCeiling,
                    constraintKind, outcome>>

RetryPoint ==
    /\ phase = "point"
    /\ attemptsRemaining > 1
    /\ attemptsRemaining' = attemptsRemaining - 1
    /\ UNCHANGED <<phase, ceilingEffective, staticOwner,
                    pointStepsRemaining, initialPointSteps, capacityDemand,
                    presentedAboveCertificate,
                    allocationCurrent,
                    initialCapacityDemand, allocationCanReplaceCeiling,
                    constraintKind, outcome>>

ConstrainPoint ==
    /\ phase = "point"
    /\ attemptsRemaining = 1
    /\ phase' = "terminal"
    /\ staticOwner' = FALSE
    /\ attemptsRemaining' = 0
    /\ constraintKind' = "point"
    /\ outcome' = "constrained"
    /\ UNCHANGED <<ceilingEffective, pointStepsRemaining,
                    initialPointSteps, capacityDemand,
                    presentedAboveCertificate,
                    allocationCurrent,
                    initialCapacityDemand, allocationCanReplaceCeiling>>

\* A point-population transition invalidates its predecessor allocation.
\* Rebuild that occurrence-local plan before capacity may measure or constrain
\* the new population.  Planning is deterministic bounded work and owns no
\* presentation retry by itself.
SelectAllocation ==
    /\ phase = "select"
    /\ ~structuralFrontier
    /\ ~ceilingEffective
    /\ pointStepsRemaining = 0
    /\ ~allocationCurrent
    /\ phase' = "allocate"
    /\ UNCHANGED <<ceilingEffective, staticOwner, pointStepsRemaining,
                    initialPointSteps, capacityDemand,
                    presentedAboveCertificate, allocationCurrent,
                    initialCapacityDemand, allocationCanReplaceCeiling,
                    attemptsRemaining, constraintKind, outcome>>

CompleteAllocation ==
    /\ phase = "allocate"
    /\ phase' = "select"
    /\ allocationCurrent' = TRUE
    /\ UNCHANGED <<ceilingEffective, staticOwner, pointStepsRemaining,
                    initialPointSteps, capacityDemand,
                    presentedAboveCertificate,
                    initialCapacityDemand, allocationCanReplaceCeiling,
                    attemptsRemaining, constraintKind, outcome>>

\* A current exact multi-occurrence allocation owns the bounded capacity
\* search when either pixel demand exceeds its selection or structural
\* publication made the presented selection richer than its prior
\* certificate.  The latter is recertification work, not permission to drop
\* the successor.  Neither state may settle as ready or defer to the
\* single-occurrence static ordinal owner.
SelectCapacity ==
    /\ phase = "select"
    /\ ~structuralFrontier
    /\ ~ceilingEffective
    /\ pointStepsRemaining = 0
    /\ allocationCurrent
    /\ (capacityDemand \/ presentedAboveCertificate)
    /\ phase' = "capacity"
    /\ attemptsRemaining' = PresentationAttemptLimit
    /\ UNCHANGED <<ceilingEffective, staticOwner, pointStepsRemaining,
                    initialPointSteps, capacityDemand,
                    presentedAboveCertificate,
                    allocationCurrent,
                    initialCapacityDemand, allocationCanReplaceCeiling,
                    constraintKind, outcome>>

CompleteCapacity ==
    /\ phase = "capacity"
    /\ (capacityDemand \/ presentedAboveCertificate)
    /\ phase' = "select"
    /\ capacityDemand' = FALSE
    /\ presentedAboveCertificate' = FALSE
    /\ attemptsRemaining' = PresentationAttemptLimit
    /\ UNCHANGED <<ceilingEffective, staticOwner, pointStepsRemaining,
                    initialPointSteps, allocationCurrent,
                    initialCapacityDemand,
                    allocationCanReplaceCeiling, constraintKind, outcome>>

RetryCapacity ==
    /\ phase = "capacity"
    /\ attemptsRemaining > 1
    /\ attemptsRemaining' = attemptsRemaining - 1
    /\ UNCHANGED <<phase, ceilingEffective, staticOwner,
                    pointStepsRemaining, initialPointSteps, capacityDemand,
                    presentedAboveCertificate,
                    allocationCurrent,
                    initialCapacityDemand, allocationCanReplaceCeiling,
                    constraintKind, outcome>>

ConstrainCapacity ==
    /\ phase = "capacity"
    /\ attemptsRemaining = 1
    /\ phase' = "terminal"
    /\ staticOwner' = FALSE
    /\ attemptsRemaining' = 0
    /\ constraintKind' = "capacity"
    /\ outcome' = "constrained"
    /\ UNCHANGED <<ceilingEffective, pointStepsRemaining,
                    initialPointSteps, capacityDemand,
                    presentedAboveCertificate,
                    allocationCurrent,
                    initialCapacityDemand, allocationCanReplaceCeiling>>

Settle ==
    /\ phase = "select"
    /\ ~structuralFrontier
    /\ ~ceilingEffective
    /\ pointStepsRemaining = 0
    /\ ~capacityDemand
    /\ ~presentedAboveCertificate
    /\ allocationCurrent
    /\ ~staticOwner
    /\ phase' = "terminal"
    /\ staticOwner' = FALSE
    /\ outcome' = "ready"
    /\ UNCHANGED <<ceilingEffective, pointStepsRemaining,
                    initialPointSteps, capacityDemand,
                    presentedAboveCertificate,
                    allocationCurrent,
                    initialCapacityDemand, allocationCanReplaceCeiling,
                    attemptsRemaining,
                    constraintKind>>

QualityNext ==
    \/ SelectStatic
    \/ CompleteStatic
    \/ ConstrainStatic
    \/ SelectReconciliation
    \/ CompleteReconciliation
    \/ RetryReconciliation
    \/ ConstrainReconciliation
    \/ SelectPoint
    \/ CompletePoint
    \/ RetryPoint
    \/ ConstrainPoint
    \/ SelectAllocation
    \/ CompleteAllocation
    \/ SelectCapacity
    \/ CompleteCapacity
    \/ RetryCapacity
    \/ ConstrainCapacity
    \/ Settle

\* Quality actions deliberately do not repeat the three structural counters.
\* This wrapper completes each action before it is used by Next or fairness;
\* applying WF_vars to the partial action would leave primed counters
\* unconstrained and is not the production transition being promised.
WithStableStructuralCounters(Action) ==
    /\ Action
    /\ UNCHANGED <<structuralFallbackCount, terminalFailureCount,
                    terminalAggregateCount>>

Next ==
    \/ SelectStructural
    \/ CompleteStructural
    \/ WithStableStructuralCounters(QualityNext)

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(SelectStructural)
    /\ WF_vars(CompleteStructural)
    /\ WF_vars(WithStableStructuralCounters(SelectStatic))
    /\ WF_vars(WithStableStructuralCounters(CompleteStatic))
    /\ WF_vars(WithStableStructuralCounters(ConstrainStatic))
    /\ WF_vars(WithStableStructuralCounters(SelectReconciliation))
    /\ WF_vars(WithStableStructuralCounters(CompleteReconciliation))
    /\ WF_vars(WithStableStructuralCounters(RetryReconciliation))
    /\ WF_vars(WithStableStructuralCounters(ConstrainReconciliation))
    /\ WF_vars(WithStableStructuralCounters(SelectPoint))
    /\ WF_vars(WithStableStructuralCounters(CompletePoint))
    /\ WF_vars(WithStableStructuralCounters(RetryPoint))
    /\ WF_vars(WithStableStructuralCounters(ConstrainPoint))
    /\ WF_vars(WithStableStructuralCounters(SelectAllocation))
    /\ WF_vars(WithStableStructuralCounters(CompleteAllocation))
    /\ WF_vars(WithStableStructuralCounters(SelectCapacity))
    /\ WF_vars(WithStableStructuralCounters(CompleteCapacity))
    /\ WF_vars(WithStableStructuralCounters(RetryCapacity))
    /\ WF_vars(WithStableStructuralCounters(ConstrainCapacity))
    /\ WF_vars(WithStableStructuralCounters(Settle))

PointRequiresCeilingFree == phase = "point" => ~ceilingEffective

StructuralDemandPrecedesQuality ==
    phase = "select" /\ structuralFrontier =>
        /\ ENABLED SelectStructural
        /\ ~ENABLED SelectReconciliation
        /\ ~ENABLED SelectPoint
        /\ ~ENABLED SelectAllocation
        /\ ~ENABLED SelectCapacity
        /\ ~ENABLED SelectStatic
        /\ ~ENABLED Settle

StructuralOwnerHasExecutableSuccessor ==
    phase = "structural" => ENABLED CompleteStructural

TerminalAggregatesDoNotCreateStructuralDebt ==
    terminalAggregateCount > 0 /\
    structuralFallbackCount <= terminalFailureCount =>
        /\ ~structuralFrontier
        /\ ~ENABLED SelectStructural

\* Pending point-population work owns selection before either allocation or
\* capacity may begin.  The C++ refinement permits an already active capacity
\* sample to finish atomically, but completion must not derive a new capacity
\* search while this owner remains pending.
PointDemandPrecedesNewCapacity ==
    phase = "select" /\ ~structuralFrontier /\ ~ceilingEffective /\
    pointStepsRemaining > 0 =>
        /\ ENABLED SelectPoint
        /\ ~ENABLED SelectAllocation
        /\ ~ENABLED SelectCapacity

PointOwnerHasExecutableSuccessor ==
    phase = "point" => ENABLED CompletePoint

\* No point-bracket progress is possible while a global ceiling remains.
CeilingPreservesPointDemand ==
    ceilingEffective => pointStepsRemaining = initialPointSteps

CeilingPreservesCapacityDemand ==
    ceilingEffective => capacityDemand = initialCapacityDemand

UncertifiedPresentationHasCapacitySuccessor ==
    phase = "select" /\ ~structuralFrontier /\ ~ceilingEffective /\
    pointStepsRemaining = 0 /\ presentedAboveCertificate =>
        ENABLED SelectCapacity \/ ENABLED SelectAllocation

StaleAllocationHasPlanningSuccessor ==
    phase = "select" /\ ~structuralFrontier /\ ~ceilingEffective /\
    pointStepsRemaining = 0 /\ ~allocationCurrent =>
        ENABLED SelectAllocation

\* Selecting the allocation owner is not merely a status transition.  It
\* must expose an executable complete-population producer.  In the C++
\* refinement this means the retained-allocation request also resets the
\* mechanical admission cursor; otherwise an initialized predecessor cursor
\* can leave the owner selected while every pass remains ordinary work.
AllocationOwnerHasExecutableSuccessor ==
    phase = "allocate" => ENABLED CompleteAllocation

TerminalIsQuiescent ==
    phase = "terminal" =>
        /\ outcome \in {"ready", "constrained"}
        /\ ~structuralFrontier
        /\ (~ceilingEffective \/
             (outcome = "constrained" /\ constraintKind = "reconcile"))
        /\ ~staticOwner

ReadyHasNoQualityDebt ==
    outcome = "ready" =>
        /\ phase = "terminal"
        /\ ~ceilingEffective
        /\ ~structuralFrontier
        /\ ~staticOwner
        /\ pointStepsRemaining = 0
        /\ allocationCurrent
        /\ ~capacityDemand
        /\ ~presentedAboveCertificate
        /\ constraintKind = "none"

ConstraintIsWitnessed ==
    outcome = "constrained" =>
        /\ phase = "terminal"
        /\ attemptsRemaining = 0
        /\ constraintKind \in ConstraintKinds \ {"none"}

ConstraintKindMatchesOwner ==
    /\ constraintKind = "static" => outcome = "constrained"
    /\ constraintKind = "reconcile" => outcome = "constrained"
    /\ constraintKind = "point" => outcome = "constrained"
    /\ constraintKind = "capacity" => outcome = "constrained"

CapacityConstraintHasCurrentAllocation ==
    constraintKind = "capacity" => allocationCurrent

TerminalCeilingIsWitnessed ==
    phase = "terminal" /\ ceilingEffective =>
        /\ outcome = "constrained"
        /\ constraintKind = "reconcile"
        /\ ~allocationCanReplaceCeiling

Ready == phase = "terminal" /\ outcome = "ready"
Constrained == phase = "terminal" /\ outcome = "constrained"
Failed == FALSE
Terminal == Ready \/ Constrained \/ Failed

\* Every action in Next changes at least one variable.  Testing Next directly
\* is therefore equivalent to testing its non-stuttering form and avoids a
\* TLC 2.19 evaluator failure on ENABLED <<Next>>_vars for this wide tuple.
DeadlockOnlyAtTerminal == ~ENABLED Next => Terminal

EventuallyTerminal == <> (phase = "terminal")

=============================================================================
