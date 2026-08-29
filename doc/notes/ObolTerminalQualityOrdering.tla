-------------------- MODULE ObolTerminalQualityOrdering --------------------
\* Composition contract for the terminal renderer-ceiling, static-quality,
\* and point-quality owners.
\*
\* A renderer-wide progressive ceiling is a transient presentation guard.
\* Point refinement changes the occurrence population hidden behind that
\* guard.  The point owner therefore cannot commit or measure a finer point
\* cut until the static owner has retired and the ceiling has been translated
\* into occurrence-local cuts.  Invalid completed presentations consume a
\* finite allowance and end at an explicit constrained outcome.

EXTENDS Naturals, TLC

CONSTANTS PointStepLimit, PresentationAttemptLimit

ASSUME PointStepLimit > 0
ASSUME PresentationAttemptLimit > 0

Phases == {"select", "static", "reconcile", "point", "terminal"}
Outcomes == {"none", "ready", "constrained"}

VARIABLES phase,
          ceilingEffective,
          staticOwner,
          pointStepsRemaining,
          initialPointSteps,
          attemptsRemaining,
          outcome

vars == <<phase, ceilingEffective, staticOwner, pointStepsRemaining,
          initialPointSteps, attemptsRemaining, outcome>>

TypeOK ==
    /\ phase \in Phases
    /\ ceilingEffective \in BOOLEAN
    /\ staticOwner \in BOOLEAN
    /\ pointStepsRemaining \in 0..PointStepLimit
    /\ initialPointSteps \in 0..PointStepLimit
    /\ attemptsRemaining \in 0..PresentationAttemptLimit
    /\ outcome \in Outcomes

Init ==
    /\ phase = "select"
    /\ ceilingEffective \in BOOLEAN
    /\ staticOwner \in BOOLEAN
    /\ pointStepsRemaining \in 0..PointStepLimit
    /\ initialPointSteps = pointStepsRemaining
    /\ attemptsRemaining = PresentationAttemptLimit
    /\ outcome = "none"

\* Static quality may own a bounded candidate after the renderer-wide ceiling
\* has been translated.  The modeled ceiling is a multi-occurrence global
\* control, not the occurrence-local ordinal used by a one-object staircase.
SelectStatic ==
    /\ phase = "select"
    /\ ~ceilingEffective
    /\ staticOwner
    /\ pointStepsRemaining = 0
    /\ phase' = "static"
    /\ UNCHANGED <<ceilingEffective, staticOwner, pointStepsRemaining,
                    initialPointSteps, attemptsRemaining, outcome>>

CompleteStatic ==
    /\ phase = "static"
    /\ phase' = "select"
    /\ staticOwner' = FALSE
    /\ UNCHANGED <<ceilingEffective, pointStepsRemaining,
                    initialPointSteps, attemptsRemaining, outcome>>

\* An effective multi-occurrence renderer ceiling owns the sole successor.
\* Point and static demand remain recorded but cannot own a framebuffer.
SelectReconciliation ==
    /\ phase = "select"
    /\ ceilingEffective
    /\ phase' = "reconcile"
    /\ attemptsRemaining' = PresentationAttemptLimit
    /\ UNCHANGED <<ceilingEffective, staticOwner, pointStepsRemaining,
                    initialPointSteps, outcome>>

CompleteReconciliation ==
    /\ phase = "reconcile"
    /\ ceilingEffective
    /\ phase' = "select"
    /\ ceilingEffective' = FALSE
    /\ attemptsRemaining' = PresentationAttemptLimit
    /\ UNCHANGED <<staticOwner, pointStepsRemaining,
                    initialPointSteps, outcome>>

RetryReconciliation ==
    /\ phase = "reconcile"
    /\ attemptsRemaining > 1
    /\ attemptsRemaining' = attemptsRemaining - 1
    /\ UNCHANGED <<phase, ceilingEffective, staticOwner,
                    pointStepsRemaining, initialPointSteps, outcome>>

ConstrainReconciliation ==
    /\ phase = "reconcile"
    /\ attemptsRemaining = 1
    /\ phase' = "terminal"
    /\ attemptsRemaining' = 0
    /\ outcome' = "constrained"
    /\ UNCHANGED <<ceilingEffective, staticOwner, pointStepsRemaining,
                    initialPointSteps>>

\* Point refinement is enabled only after the complete renderer control vector
\* is ceiling-free.  Each exact presentation consumes one finite bracket step.
SelectPoint ==
    /\ phase = "select"
    /\ ~ceilingEffective
    /\ pointStepsRemaining > 0
    /\ phase' = "point"
    /\ attemptsRemaining' = PresentationAttemptLimit
    /\ UNCHANGED <<ceilingEffective, staticOwner, pointStepsRemaining,
                    initialPointSteps, outcome>>

CompletePoint ==
    /\ phase = "point"
    /\ ~ceilingEffective
    /\ pointStepsRemaining > 0
    /\ phase' = "select"
    /\ pointStepsRemaining' = pointStepsRemaining - 1
    /\ attemptsRemaining' = PresentationAttemptLimit
    /\ UNCHANGED <<ceilingEffective, staticOwner, initialPointSteps, outcome>>

RetryPoint ==
    /\ phase = "point"
    /\ attemptsRemaining > 1
    /\ attemptsRemaining' = attemptsRemaining - 1
    /\ UNCHANGED <<phase, ceilingEffective, staticOwner,
                    pointStepsRemaining, initialPointSteps, outcome>>

ConstrainPoint ==
    /\ phase = "point"
    /\ attemptsRemaining = 1
    /\ phase' = "terminal"
    /\ attemptsRemaining' = 0
    /\ outcome' = "constrained"
    /\ UNCHANGED <<ceilingEffective, staticOwner, pointStepsRemaining,
                    initialPointSteps>>

Settle ==
    /\ phase = "select"
    /\ ~ceilingEffective
    /\ pointStepsRemaining = 0
    /\ ~staticOwner
    /\ phase' = "terminal"
    /\ staticOwner' = FALSE
    /\ outcome' = "ready"
    /\ UNCHANGED <<ceilingEffective, pointStepsRemaining,
                    initialPointSteps, attemptsRemaining>>

Next ==
    \/ SelectStatic
    \/ CompleteStatic
    \/ SelectReconciliation
    \/ CompleteReconciliation
    \/ RetryReconciliation
    \/ ConstrainReconciliation
    \/ SelectPoint
    \/ CompletePoint
    \/ RetryPoint
    \/ ConstrainPoint
    \/ Settle

Spec == Init /\ [][Next]_vars /\ WF_vars(Next)

PointRequiresCeilingFree == phase = "point" => ~ceilingEffective

\* No point-bracket progress is possible while a global ceiling remains.
CeilingPreservesPointDemand ==
    ceilingEffective => pointStepsRemaining = initialPointSteps

TerminalIsQuiescent ==
    phase = "terminal" => outcome \in {"ready", "constrained"}

ReadyHasNoQualityDebt ==
    outcome = "ready" =>
        /\ phase = "terminal"
        /\ ~ceilingEffective
        /\ ~staticOwner
        /\ pointStepsRemaining = 0

ConstraintIsWitnessed ==
    outcome = "constrained" =>
        /\ phase = "terminal"
        /\ attemptsRemaining = 0

EventuallyTerminal == <> (phase = "terminal")

=============================================================================
