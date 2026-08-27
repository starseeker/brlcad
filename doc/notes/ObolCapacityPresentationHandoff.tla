------------------- MODULE ObolCapacityPresentationHandoff -------------------
\* Ordering contract between one frozen occurrence allocation and the
\* renderer-capacity sample which measures it.
\*
\* The production allocator may change many occurrence-local cuts in one
\* bounded pass while an older renderer-wide ceiling still protects the last
\* coherent framebuffer.  The changed pass is valid application evidence: it
\* must be allowed to hand off to a ceiling-free successor frame.  A sample is
\* consumable only after every assigned cut is active and the ceiling is
\* absent or provably inert.  Frames observed before that point are not samples
\* and cannot spend the bounded sample allowance.
\*
\* A complete allocation may itself exceed its certified budget.  That state
\* belongs to the same handoff: a finite local-representation reduction either
\* makes it fit or publishes an explicit bounded constraint.  It must never
\* fall back into another identical generic allocation.
\*
\* This model proves ownership, safety, and liveness for a frozen allocation.
\* It deliberately does not model frame-time classification, aggregate cost,
\* or visual quality; those remain pure C++ and graphical tests.

EXTENDS Naturals, TLC

CONSTANTS SampleLimit, ReductionLimit

ASSUME SampleLimit > 0
ASSUME ReductionLimit > 0

Phases == {"apply", "handoff", "measure", "terminal"}
CeilingStates == {"effective", "inert", "none"}
BudgetStates == {"fits", "over"}

VARIABLES phase,
          cutsApplied,
          assignedCutsDrawable,
          ceilingState,
          allocationBudgetState,
          reductionAttempts,
          samplesRemaining,
          samplesConsumed,
          invalidFramesObserved,
          certificatePublished,
          constraintPublished

vars == <<phase, cutsApplied, assignedCutsDrawable, ceilingState,
          allocationBudgetState,
          reductionAttempts, samplesRemaining, samplesConsumed,
          invalidFramesObserved, certificatePublished,
          constraintPublished>>

PresentationExact == cutsApplied /\ ceilingState # "effective"

TypeOK ==
    /\ phase \in Phases
    /\ cutsApplied \in BOOLEAN
    /\ assignedCutsDrawable \in BOOLEAN
    /\ ceilingState \in CeilingStates
    /\ allocationBudgetState \in BudgetStates
    /\ reductionAttempts \in 0..ReductionLimit
    /\ samplesRemaining \in 0..SampleLimit
    /\ samplesConsumed \in 0..SampleLimit
    /\ invalidFramesObserved \in 0..2
    /\ certificatePublished \in BOOLEAN
    /\ constraintPublished \in BOOLEAN

SampleAccounting ==
    samplesRemaining + samplesConsumed = SampleLimit

AppliedCutsAreDrawable ==
    cutsApplied => assignedCutsDrawable

SampleRequiresExactPresentation ==
    samplesConsumed > 0 =>
        /\ PresentationExact
        /\ allocationBudgetState = "fits"

TerminalCertificateIsExact ==
    certificatePublished =>
        /\ phase = "terminal"
        /\ PresentationExact
        /\ allocationBudgetState = "fits"
        /\ samplesRemaining = 0
        /\ samplesConsumed = SampleLimit

TerminalConstraintIsBounded ==
    constraintPublished =>
        /\ phase = "terminal"
        /\ cutsApplied
        /\ allocationBudgetState = "over"
        /\ reductionAttempts = ReductionLimit

NoPrematureTerminal ==
    phase = "terminal" => certificatePublished \/ constraintPublished

Init ==
    /\ assignedCutsDrawable \in BOOLEAN
    /\ cutsApplied \in BOOLEAN
    /\ cutsApplied => assignedCutsDrawable
    /\ ceilingState \in {"effective", "inert"}
    /\ allocationBudgetState \in BudgetStates
    /\ reductionAttempts = 0
    /\ phase = IF cutsApplied THEN "handoff" ELSE "apply"
    /\ samplesRemaining = SampleLimit
    /\ samplesConsumed = 0
    /\ invalidFramesObserved = 0
    /\ certificatePublished = FALSE
    /\ constraintPublished = FALSE

ConstrainAllocationToAvailability ==
    /\ phase = "apply"
    /\ ~cutsApplied
    /\ ~assignedCutsDrawable
    /\ assignedCutsDrawable' = TRUE
    /\ UNCHANGED <<phase, cutsApplied, ceilingState,
                    allocationBudgetState, reductionAttempts,
                    samplesRemaining, samplesConsumed,
                    invalidFramesObserved, certificatePublished,
                    constraintPublished>>

ApplyAssignedCuts ==
    /\ phase = "apply"
    /\ ~cutsApplied
    /\ assignedCutsDrawable
    /\ cutsApplied' = TRUE
    /\ phase' = "handoff"
    /\ UNCHANGED <<assignedCutsDrawable, ceilingState,
                    allocationBudgetState, reductionAttempts,
                    samplesRemaining, samplesConsumed,
                    invalidFramesObserved, certificatePublished,
                    constraintPublished>>

ReduceOverBudgetAllocation ==
    /\ phase = "handoff"
    /\ cutsApplied
    /\ allocationBudgetState = "over"
    /\ reductionAttempts < ReductionLimit
    /\ reductionAttempts' = reductionAttempts + 1
    /\ allocationBudgetState' \in BudgetStates
    /\ UNCHANGED <<phase, cutsApplied, assignedCutsDrawable, ceilingState,
                    samplesRemaining,
                    samplesConsumed, invalidFramesObserved,
                    certificatePublished, constraintPublished>>

ReconcileCeiling ==
    /\ phase = "handoff"
    /\ cutsApplied
    /\ allocationBudgetState = "fits"
    /\ ceilingState' = IF ceilingState = "effective"
                           THEN "none"
                           ELSE ceilingState
    /\ phase' = "measure"
    /\ UNCHANGED <<cutsApplied, assignedCutsDrawable,
                    allocationBudgetState, reductionAttempts,
                    samplesRemaining, samplesConsumed,
                    invalidFramesObserved, certificatePublished,
                    constraintPublished>>

PublishBoundedConstraint ==
    /\ phase = "handoff"
    /\ cutsApplied
    /\ allocationBudgetState = "over"
    /\ reductionAttempts = ReductionLimit
    /\ phase' = "terminal"
    /\ constraintPublished' = TRUE
    /\ UNCHANGED <<cutsApplied, assignedCutsDrawable, ceilingState,
                    allocationBudgetState,
                    reductionAttempts, samplesRemaining, samplesConsumed,
                    invalidFramesObserved, certificatePublished>>

ObserveInvalidFrame ==
    /\ phase \in {"apply", "handoff"}
    /\ invalidFramesObserved < 2
    /\ invalidFramesObserved' = invalidFramesObserved + 1
    /\ UNCHANGED <<phase, cutsApplied, assignedCutsDrawable, ceilingState,
                    allocationBudgetState, reductionAttempts,
                    samplesRemaining, samplesConsumed,
                    certificatePublished, constraintPublished>>

ConsumeSample ==
    /\ phase = "measure"
    /\ PresentationExact
    /\ allocationBudgetState = "fits"
    /\ samplesRemaining > 0
    /\ samplesRemaining' = samplesRemaining - 1
    /\ samplesConsumed' = samplesConsumed + 1
    /\ IF samplesRemaining = 1
          THEN /\ phase' = "terminal"
               /\ certificatePublished' = TRUE
          ELSE /\ phase' = phase
               /\ certificatePublished' = FALSE
    /\ UNCHANGED <<cutsApplied, assignedCutsDrawable, ceilingState,
                    allocationBudgetState,
                    reductionAttempts, invalidFramesObserved,
                    constraintPublished>>

Next ==
    \/ ConstrainAllocationToAvailability
    \/ ApplyAssignedCuts
    \/ ReduceOverBudgetAllocation
    \/ ReconcileCeiling
    \/ PublishBoundedConstraint
    \/ ObserveInvalidFrame
    \/ ConsumeSample

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(ConstrainAllocationToAvailability)
    /\ WF_vars(ApplyAssignedCuts)
    /\ WF_vars(ReduceOverBudgetAllocation)
    /\ WF_vars(ReconcileCeiling)
    /\ WF_vars(PublishBoundedConstraint)
    /\ WF_vars(ConsumeSample)

EventuallyReady == <>(certificatePublished \/ constraintPublished)

=============================================================================
