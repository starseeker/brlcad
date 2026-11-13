---------------- MODULE ObolTerminalConvergenceComposition ----------------
\* Composition contract for the closed-input tail of one LoD view epoch.
\*
\* The focused capacity, static-quality, presentation-handoff, and compaction
\* models prove that their own transactions terminate.  That is not enough:
\* production owners may create successor work for one another.  In
\* particular, a static headroom trial may reopen capacity allocation, an
\* over-budget allocation may enter local representation reconciliation, and
\* quiet compaction follows foreground convergence.  This model makes those
\* cross-owner edges explicit and requires every re-entry to consume a member
\* of a finite semantic set.
\*
\* MaxCapacityCandidates is the bounded candidate set for one capacity epoch.
\* MaxStaticTrials bounds distinct static headroom populations.  A static
\* trial may start a new capacity epoch, but consumes itself when doing so.
\* MaxProjectionSteps bounds the immutable occurrence projection census which
\* may precede an allocation.  That census is a prerequisite request, not an
\* allocation plan, and therefore cannot consume the one-plan-per-revision
\* authority.  MaxPreparationSteps bounds the immutable command-preparation units for one
\* exact allocation/presentation target.  MaxReductionSteps bounds strictly
\* coarser local representations for an over-budget allocation.
\* MaxCompactionSteps bounds the quiet resident cursor.  None of these
\* counters is a wall-clock timeout.
\*
\* Geometry quality and frame duration remain outside this model.  The
\* production refinement is responsible for mapping each decrement to a real
\* candidate, population, threshold, or cursor advance rather than an
\* arbitrary retry count.
\*
\* allocationRevision is exclusively the six-domain admission revision.
\* Applying a plan may advance concrete CAD-binding or resident-demand
\* bookkeeping revisions, but those are transaction-staleness witnesses, not
\* authority to plan again.  They therefore do not appear in this composition:
\* CompleteAllocation commits exactly once until a real successor action
\* advances allocationRevision.

EXTENDS Naturals, TLC

CONSTANTS MaxVisibilitySteps, MaxCapacityCandidates, MaxStaticTrials,
          MaxProjectionSteps,
          MaxPreparationSteps, MaxReductionSteps, MaxCompactionSteps

ASSUME MaxVisibilitySteps > 0
ASSUME MaxCapacityCandidates > 0
ASSUME MaxStaticTrials > 0
ASSUME MaxProjectionSteps > 0
ASSUME MaxPreparationSteps > 0
ASSUME MaxReductionSteps > 0
ASSUME MaxCompactionSteps > 0

Phases == {"visibility", "visibilityPresentation", "allocation",
           "projection", "presentation", "capacity", "static", "reconcile",
           "compaction", "terminal"}
Outcomes == {"active", "ready", "constrained"}

VARIABLES phase,
          visibilityStepsRemaining,
          visibilityCurrent,
          capacityCandidatesRemaining,
          staticTrialsRemaining,
          projectionStepsRemaining,
          preparationStepsRemaining,
          reductionStepsRemaining,
          compactionStepsRemaining,
          allocationOverBudget,
          ceilingEffective,
          ceilingConstraint,
          staticCandidateActive,
          staticConstraintCommitted,
          exactPresentation,
          projectionCurrent,
          allocationRevision,
          committedAllocationRevision,
          presentedRevision,
          outcome

vars == <<phase, visibilityStepsRemaining, visibilityCurrent,
          capacityCandidatesRemaining, staticTrialsRemaining,
          projectionStepsRemaining, preparationStepsRemaining,
          reductionStepsRemaining,
          compactionStepsRemaining,
          allocationOverBudget, ceilingEffective, ceilingConstraint,
          staticCandidateActive, staticConstraintCommitted,
          exactPresentation, projectionCurrent, allocationRevision,
          committedAllocationRevision, presentedRevision, outcome>>

Ready == phase = "terminal" /\ outcome = "ready"
Constrained == phase = "terminal" /\ outcome = "constrained"
Failed == FALSE
Terminal == Ready \/ Constrained \/ Failed

ForegroundPhase ==
    phase \in {"visibility", "visibilityPresentation", "allocation",
               "projection", "presentation", "capacity", "static", "reconcile"}

TypeOK ==
    /\ phase \in Phases
    /\ visibilityStepsRemaining \in 0..MaxVisibilitySteps
    /\ visibilityCurrent \in BOOLEAN
    /\ capacityCandidatesRemaining \in 0..MaxCapacityCandidates
    /\ staticTrialsRemaining \in 0..MaxStaticTrials
    /\ projectionStepsRemaining \in 0..MaxProjectionSteps
    /\ preparationStepsRemaining \in 0..MaxPreparationSteps
    /\ reductionStepsRemaining \in 0..MaxReductionSteps
    /\ compactionStepsRemaining \in 0..MaxCompactionSteps
    /\ allocationOverBudget \in BOOLEAN
    /\ ceilingEffective \in BOOLEAN
    /\ ceilingConstraint \in BOOLEAN
    /\ ceilingConstraint => ceilingEffective
    /\ staticCandidateActive \in BOOLEAN
    /\ staticConstraintCommitted \in BOOLEAN
    /\ staticConstraintCommitted => ~staticCandidateActive
    /\ exactPresentation \in BOOLEAN
    /\ projectionCurrent \in BOOLEAN
    /\ allocationRevision \in Nat
    /\ committedAllocationRevision \in Nat
    /\ committedAllocationRevision <= allocationRevision
    /\ presentedRevision \in Nat
    /\ presentedRevision <= allocationRevision
    /\ outcome \in Outcomes

Init ==
    /\ visibilityCurrent \in BOOLEAN
    /\ phase = IF visibilityCurrent THEN "allocation" ELSE "visibility"
    /\ visibilityStepsRemaining =
        IF visibilityCurrent THEN 0 ELSE MaxVisibilitySteps
    /\ capacityCandidatesRemaining = MaxCapacityCandidates
    /\ staticTrialsRemaining = MaxStaticTrials
    /\ projectionStepsRemaining = 0
    /\ preparationStepsRemaining = 0
    /\ reductionStepsRemaining = MaxReductionSteps
    /\ compactionStepsRemaining = MaxCompactionSteps
    /\ allocationOverBudget \in BOOLEAN
    /\ ceilingEffective \in BOOLEAN
    /\ ceilingConstraint = FALSE
    /\ staticCandidateActive = FALSE
    /\ staticConstraintCommitted = FALSE
    /\ exactPresentation = FALSE
    /\ projectionCurrent \in BOOLEAN
    /\ allocationRevision = 1
    /\ committedAllocationRevision = 0
    /\ presentedRevision = 0
    /\ outcome = "active"

\* An exact visibility delta is a finite metadata scan followed by an exact
\* framebuffer classification.  Allocation is downstream of both.  The
\* production failure represented here restored the database census and then
\* allocated immediately from the predecessor frame; meshes retired while the
\* subpath was hidden consequently had no representation and no producer.
AdvanceVisibilityDelta ==
    /\ phase = "visibility"
    /\ visibilityStepsRemaining > 0
    /\ visibilityStepsRemaining' = visibilityStepsRemaining - 1
    /\ UNCHANGED <<phase, visibilityCurrent,
                    capacityCandidatesRemaining, staticTrialsRemaining,
                    projectionStepsRemaining, preparationStepsRemaining,
                    reductionStepsRemaining, compactionStepsRemaining,
                    allocationOverBudget, ceilingEffective, ceilingConstraint,
                    staticCandidateActive, staticConstraintCommitted,
                    exactPresentation, projectionCurrent, allocationRevision,
                    committedAllocationRevision, presentedRevision, outcome>>

CompleteVisibilityDelta ==
    /\ phase = "visibility"
    /\ visibilityStepsRemaining = 0
    /\ phase' = "visibilityPresentation"
    /\ preparationStepsRemaining' = MaxPreparationSteps
    /\ exactPresentation' = FALSE
    /\ UNCHANGED <<visibilityStepsRemaining, visibilityCurrent,
                    capacityCandidatesRemaining, staticTrialsRemaining,
                    projectionStepsRemaining, reductionStepsRemaining,
                    compactionStepsRemaining, allocationOverBudget,
                    ceilingEffective, ceilingConstraint,
                    staticCandidateActive, staticConstraintCommitted,
                    projectionCurrent, allocationRevision,
                    committedAllocationRevision, presentedRevision, outcome>>

AdvanceVisibilityPresentation ==
    /\ phase = "visibilityPresentation"
    /\ preparationStepsRemaining > 0
    /\ preparationStepsRemaining' = preparationStepsRemaining - 1
    /\ UNCHANGED <<phase, visibilityStepsRemaining, visibilityCurrent,
                    capacityCandidatesRemaining, staticTrialsRemaining,
                    projectionStepsRemaining, reductionStepsRemaining,
                    compactionStepsRemaining, allocationOverBudget,
                    ceilingEffective, ceilingConstraint,
                    staticCandidateActive, staticConstraintCommitted,
                    exactPresentation, projectionCurrent, allocationRevision,
                    committedAllocationRevision, presentedRevision, outcome>>

CompleteVisibilityPresentation ==
    /\ phase = "visibilityPresentation"
    /\ preparationStepsRemaining = 0
    /\ phase' = "allocation"
    /\ visibilityCurrent' = TRUE
    /\ exactPresentation' = TRUE
    /\ UNCHANGED <<visibilityStepsRemaining, capacityCandidatesRemaining,
                    staticTrialsRemaining, projectionStepsRemaining,
                    preparationStepsRemaining, reductionStepsRemaining,
                    compactionStepsRemaining, allocationOverBudget,
                    ceilingEffective, ceilingConstraint,
                    staticCandidateActive, staticConstraintCommitted,
                    projectionCurrent, allocationRevision,
                    committedAllocationRevision, presentedRevision, outcome>>

CompleteAllocation ==
    /\ phase = "allocation"
    /\ visibilityCurrent
    /\ projectionCurrent
    /\ committedAllocationRevision < allocationRevision
    /\ phase' = "presentation"
    /\ preparationStepsRemaining' = MaxPreparationSteps
    /\ exactPresentation' = FALSE
    /\ committedAllocationRevision' = allocationRevision
    /\ UNCHANGED <<visibilityStepsRemaining, visibilityCurrent,
                    capacityCandidatesRemaining, staticTrialsRemaining,
                    projectionStepsRemaining, reductionStepsRemaining,
                    compactionStepsRemaining,
                    allocationOverBudget, ceilingEffective,
                    ceilingConstraint, staticCandidateActive,
                    staticConstraintCommitted, projectionCurrent,
                    allocationRevision,
                    presentedRevision, outcome>>

\* A retained payload with stale view-dependent evidence yields a typed
\* projection-refresh frontier.  It does not begin, stage, or commit an
\* allocation plan.  Once the finite census completes, exactly one successor
\* allocation owns the unchanged semantic revision.
RequestProjectionRefresh ==
    /\ phase = "allocation"
    /\ ~projectionCurrent
    /\ committedAllocationRevision < allocationRevision
    /\ phase' = "projection"
    /\ projectionStepsRemaining' = MaxProjectionSteps
    /\ UNCHANGED <<visibilityStepsRemaining, visibilityCurrent,
                    capacityCandidatesRemaining, staticTrialsRemaining,
                    preparationStepsRemaining, reductionStepsRemaining,
                    compactionStepsRemaining, allocationOverBudget,
                    ceilingEffective, ceilingConstraint,
                    staticCandidateActive, staticConstraintCommitted,
                    exactPresentation,
                    projectionCurrent, allocationRevision,
                    committedAllocationRevision, presentedRevision, outcome>>

AdvanceProjectionRefresh ==
    /\ phase = "projection"
    /\ projectionStepsRemaining > 0
    /\ projectionStepsRemaining' = projectionStepsRemaining - 1
    /\ UNCHANGED <<phase, visibilityStepsRemaining, visibilityCurrent,
                    capacityCandidatesRemaining,
                    staticTrialsRemaining, preparationStepsRemaining,
                    reductionStepsRemaining, compactionStepsRemaining,
                    allocationOverBudget, ceilingEffective,
                    ceilingConstraint, staticCandidateActive,
                    staticConstraintCommitted, exactPresentation,
                    projectionCurrent,
                    allocationRevision, committedAllocationRevision,
                    presentedRevision, outcome>>

CompleteProjectionRefresh ==
    /\ phase = "projection"
    /\ projectionStepsRemaining = 0
    /\ phase' = "allocation"
    /\ projectionCurrent' = TRUE
    /\ UNCHANGED <<visibilityStepsRemaining, visibilityCurrent,
                    capacityCandidatesRemaining, staticTrialsRemaining,
                    projectionStepsRemaining, preparationStepsRemaining,
                    reductionStepsRemaining, compactionStepsRemaining,
                    allocationOverBudget, ceilingEffective,
                    ceilingConstraint, staticCandidateActive,
                    staticConstraintCommitted, exactPresentation,
                    allocationRevision, committedAllocationRevision,
                    presentedRevision, outcome>>

\* SoCADAssembly command preparation is resumable across bounded endpoint
\* traversals.  An unchanged retry is enabled only when its exact-target
\* certificate consumes one remaining unit; time or frame serial growth alone
\* cannot justify another retry.
AdvancePresentationPreparation ==
    /\ phase = "presentation"
    /\ preparationStepsRemaining > 0
    /\ preparationStepsRemaining' = preparationStepsRemaining - 1
    /\ UNCHANGED <<phase, visibilityStepsRemaining, visibilityCurrent,
                    capacityCandidatesRemaining,
                    staticTrialsRemaining, projectionStepsRemaining,
                    reductionStepsRemaining,
                    compactionStepsRemaining, allocationOverBudget,
                    ceilingEffective, ceilingConstraint,
                    staticCandidateActive, staticConstraintCommitted,
                    exactPresentation,
                    projectionCurrent, allocationRevision,
                    committedAllocationRevision, presentedRevision, outcome>>

CompletePresentation ==
    /\ phase = "presentation"
    /\ preparationStepsRemaining = 0
    /\ exactPresentation' = TRUE
    /\ presentedRevision' = allocationRevision
    /\ phase' =
        IF ceilingEffective /\ ~ceilingConstraint THEN "reconcile"
        ELSE IF capacityCandidatesRemaining > 0 THEN "capacity"
        ELSE IF staticTrialsRemaining > 0 THEN "static"
        ELSE IF compactionStepsRemaining > 0 THEN "compaction"
        ELSE "terminal"
    /\ staticCandidateActive' =
        IF phase' \in {"compaction", "terminal"} THEN FALSE
        ELSE staticCandidateActive
    /\ outcome' =
        IF phase' = "terminal"
        THEN IF ceilingConstraint THEN "constrained" ELSE "ready"
        ELSE "active"
    /\ UNCHANGED <<visibilityStepsRemaining, visibilityCurrent,
                    capacityCandidatesRemaining, staticTrialsRemaining,
                    projectionStepsRemaining, preparationStepsRemaining,
                    reductionStepsRemaining, compactionStepsRemaining,
                    allocationOverBudget, ceilingEffective,
                    ceilingConstraint, staticConstraintCommitted,
                    projectionCurrent, allocationRevision,
                    committedAllocationRevision>>

\* A capacity successor names a distinct member of the bounded search set.
\* It may select an over-budget population, but cannot retry the same member.
TryCapacityCandidate ==
    /\ phase = "capacity"
    /\ capacityCandidatesRemaining > 0
    /\ capacityCandidatesRemaining' = capacityCandidatesRemaining - 1
    /\ preparationStepsRemaining' = 0
    /\ allocationRevision' = allocationRevision + 1
    /\ allocationOverBudget' \in BOOLEAN
    /\ projectionCurrent' \in BOOLEAN
    /\ phase' = "allocation"
    /\ exactPresentation' = FALSE
    /\ outcome' = "active"
    /\ UNCHANGED <<visibilityStepsRemaining, visibilityCurrent,
                    staticTrialsRemaining, projectionStepsRemaining,
                    reductionStepsRemaining,
                    compactionStepsRemaining, ceilingEffective,
                    ceilingConstraint, staticCandidateActive,
                    staticConstraintCommitted, committedAllocationRevision,
                    presentedRevision>>

\* One distinct static population may reopen a fresh bounded capacity search.
\* Consuming the static token before the re-entry is the cross-owner rank.
TryStaticPopulation ==
    /\ phase = "static"
    /\ staticTrialsRemaining > 0
    /\ staticTrialsRemaining' = staticTrialsRemaining - 1
    /\ capacityCandidatesRemaining' = MaxCapacityCandidates
    /\ preparationStepsRemaining' = 0
    /\ allocationRevision' = allocationRevision + 1
    /\ allocationOverBudget' \in BOOLEAN
    /\ projectionCurrent' \in BOOLEAN
    /\ staticCandidateActive' = TRUE
    /\ staticConstraintCommitted' = FALSE
    /\ phase' = "allocation"
    /\ exactPresentation' = FALSE
    /\ outcome' = "active"
    /\ UNCHANGED <<visibilityStepsRemaining, visibilityCurrent,
                    projectionStepsRemaining, reductionStepsRemaining,
                    compactionStepsRemaining,
                    ceilingEffective, ceilingConstraint,
                    committedAllocationRevision, presentedRevision>>

\* Reconciliation either reaches a representation which fits, consumes one
\* strictly coarser local representation, or publishes the retained-ceiling
\* constraint after the finite representation set is exhausted.
FinishFittingReconciliation ==
    /\ phase = "reconcile"
    /\ ~allocationOverBudget
    /\ ceilingEffective' = FALSE
    /\ ceilingConstraint' = FALSE
    /\ capacityCandidatesRemaining' =
        IF staticCandidateActive THEN 0 ELSE capacityCandidatesRemaining
    /\ staticTrialsRemaining' =
        IF staticCandidateActive THEN 0 ELSE staticTrialsRemaining
    /\ staticConstraintCommitted' = staticCandidateActive
    /\ staticCandidateActive' = FALSE
    /\ phase' = "presentation"
    /\ preparationStepsRemaining' = MaxPreparationSteps
    /\ exactPresentation' = FALSE
    /\ outcome' = "active"
    /\ UNCHANGED <<visibilityStepsRemaining, visibilityCurrent,
                    projectionStepsRemaining, reductionStepsRemaining,
                    compactionStepsRemaining, allocationOverBudget,
                    projectionCurrent, allocationRevision,
                    committedAllocationRevision,
                    presentedRevision>>

ReduceLocalRepresentation ==
    /\ phase = "reconcile"
    /\ allocationOverBudget
    /\ reductionStepsRemaining > 0
    /\ reductionStepsRemaining' = reductionStepsRemaining - 1
    /\ preparationStepsRemaining' = 0
    /\ allocationRevision' = allocationRevision + 1
    /\ allocationOverBudget' \in BOOLEAN
    /\ projectionCurrent' \in BOOLEAN
    /\ phase' = "allocation"
    /\ exactPresentation' = FALSE
    /\ outcome' = "active"
    /\ UNCHANGED <<visibilityStepsRemaining, visibilityCurrent,
                    capacityCandidatesRemaining, staticTrialsRemaining,
                    projectionStepsRemaining, compactionStepsRemaining,
                    ceilingEffective, ceilingConstraint,
                    staticCandidateActive, staticConstraintCommitted,
                    committedAllocationRevision, presentedRevision>>

PublishRetainedCeilingConstraint ==
    /\ phase = "reconcile"
    /\ allocationOverBudget
    /\ reductionStepsRemaining = 0
    /\ ceilingEffective
    /\ ceilingConstraint' = TRUE
    /\ capacityCandidatesRemaining' =
        IF staticCandidateActive THEN 0 ELSE capacityCandidatesRemaining
    /\ staticTrialsRemaining' =
        IF staticCandidateActive THEN 0 ELSE staticTrialsRemaining
    /\ staticConstraintCommitted' = staticCandidateActive
    /\ staticCandidateActive' = FALSE
    /\ phase' = "presentation"
    /\ preparationStepsRemaining' = MaxPreparationSteps
    /\ exactPresentation' = FALSE
    /\ outcome' = "active"
    /\ UNCHANGED <<visibilityStepsRemaining, visibilityCurrent,
                    projectionStepsRemaining, reductionStepsRemaining,
                    compactionStepsRemaining,
                    allocationOverBudget, ceilingEffective,
                    projectionCurrent, allocationRevision,
                    committedAllocationRevision, presentedRevision>>

\* Compaction is background work over the terminal foreground population.  It
\* cannot change allocation, presentation, capacity, or static-quality state.
AdvanceCompaction ==
    /\ phase = "compaction"
    /\ compactionStepsRemaining > 0
    /\ compactionStepsRemaining' = compactionStepsRemaining - 1
    /\ phase' = IF compactionStepsRemaining = 1
                 THEN "terminal" ELSE "compaction"
    /\ outcome' =
        IF phase' = "terminal"
        THEN IF ceilingConstraint THEN "constrained" ELSE "ready"
        ELSE "active"
    /\ UNCHANGED <<visibilityStepsRemaining, visibilityCurrent,
                    capacityCandidatesRemaining, staticTrialsRemaining,
                    projectionStepsRemaining, preparationStepsRemaining,
                    reductionStepsRemaining,
                    allocationOverBudget,
                    ceilingEffective, ceilingConstraint,
                    staticCandidateActive, staticConstraintCommitted,
                    exactPresentation,
                    projectionCurrent, allocationRevision,
                    committedAllocationRevision, presentedRevision>>

Next ==
    \/ AdvanceVisibilityDelta
    \/ CompleteVisibilityDelta
    \/ AdvanceVisibilityPresentation
    \/ CompleteVisibilityPresentation
    \/ CompleteAllocation
    \/ RequestProjectionRefresh
    \/ AdvanceProjectionRefresh
    \/ CompleteProjectionRefresh
    \/ AdvancePresentationPreparation
    \/ CompletePresentation
    \/ TryCapacityCandidate
    \/ TryStaticPopulation
    \/ FinishFittingReconciliation
    \/ ReduceLocalRepresentation
    \/ PublishRetainedCeilingConstraint
    \/ AdvanceCompaction

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(AdvanceVisibilityDelta)
    /\ WF_vars(CompleteVisibilityDelta)
    /\ WF_vars(AdvanceVisibilityPresentation)
    /\ WF_vars(CompleteVisibilityPresentation)
    /\ WF_vars(CompleteAllocation)
    /\ WF_vars(RequestProjectionRefresh)
    /\ WF_vars(AdvanceProjectionRefresh)
    /\ WF_vars(CompleteProjectionRefresh)
    /\ WF_vars(AdvancePresentationPreparation)
    /\ WF_vars(CompletePresentation)
    /\ WF_vars(TryCapacityCandidate)
    /\ WF_vars(TryStaticPopulation)
    /\ WF_vars(FinishFittingReconciliation)
    /\ WF_vars(ReduceLocalRepresentation)
    /\ WF_vars(PublishRetainedCeilingConstraint)
    /\ WF_vars(AdvanceCompaction)

ExactTerminalPresentation ==
    outcome # "active" =>
        /\ phase = "terminal"
        /\ visibilityCurrent
        /\ exactPresentation
        /\ projectionCurrent
        /\ committedAllocationRevision = allocationRevision
        /\ presentedRevision = allocationRevision

ProjectionRefreshDoesNotCommit ==
    phase = "projection" =>
        committedAllocationRevision < allocationRevision

VisibilityPresentationPrecedesAllocation ==
    /\ ~visibilityCurrent =>
        phase \in {"visibility", "visibilityPresentation"}
    /\ committedAllocationRevision > 0 => visibilityCurrent

PresentationHasOneCurrentPlan ==
    phase \in {"presentation", "capacity", "static", "reconcile",
               "compaction", "terminal"} =>
        committedAllocationRevision = allocationRevision

RetainedCeilingHasConstraint ==
    phase = "terminal" /\ ceilingEffective => ceilingConstraint

StaticConstraintSupersedesOlderSearch ==
    staticConstraintCommitted =>
        /\ capacityCandidatesRemaining = 0
        /\ staticTrialsRemaining = 0

TerminalHasNoActiveStaticCandidate ==
    phase = "terminal" => ~staticCandidateActive

CompactionDoesNotReopenForeground ==
    [](phase = "compaction" => [](~ForegroundPhase))

TerminalDoesNotReopen ==
    [](phase = "terminal" => [](phase = "terminal"))

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Terminal

EventuallyTerminal == <>(phase = "terminal")

=============================================================================
