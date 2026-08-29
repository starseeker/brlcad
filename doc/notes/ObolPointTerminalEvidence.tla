-------------------- MODULE ObolPointTerminalEvidence --------------------
\* Focused contract for a constrained point-quality endpoint.  Rejecting a
\* finer structural preload consumes the one decision available for the
\* current semantic population and records a terminal witness.  Re-observing
\* the same exact structural census is a no-op: it cannot erase that witness
\* and leave no owner capable of producing a successor.

EXTENDS Naturals

CONSTANTS SemanticChangeLimit, NoopObservationLimit

ASSUME SemanticChangeLimit >= 0
ASSUME NoopObservationLimit >= 0

Phases == {"active", "terminal"}
Outcomes == {"none", "constrained"}

VARIABLES phase,
          semanticChangesRemaining,
          noopObservationsRemaining,
          rejectionAvailable,
          terminalWitness,
          outcome

vars == <<phase, semanticChangesRemaining, noopObservationsRemaining,
          rejectionAvailable, terminalWitness, outcome>>

TypeOK ==
    /\ phase \in Phases
    /\ semanticChangesRemaining \in 0..SemanticChangeLimit
    /\ noopObservationsRemaining \in 0..NoopObservationLimit
    /\ rejectionAvailable \in BOOLEAN
    /\ terminalWitness \in BOOLEAN
    /\ outcome \in Outcomes

Init ==
    /\ phase = "active"
    /\ semanticChangesRemaining = SemanticChangeLimit
    /\ noopObservationsRemaining = NoopObservationLimit
    /\ rejectionAvailable = TRUE
    /\ terminalWitness = FALSE
    /\ outcome = "none"

\* The current population cannot afford its finer preload.  This transition
\* consumes the decision for that semantic epoch and records its constrained
\* endpoint atomically.
RejectFinerCandidate ==
    /\ phase = "active"
    /\ rejectionAvailable
    /\ ~terminalWitness
    /\ rejectionAvailable' = FALSE
    /\ terminalWitness' = TRUE
    /\ UNCHANGED <<phase, semanticChangesRemaining,
                    noopObservationsRemaining, outcome>>

\* A repeated exact census carries no new semantic information.  It may be
\* observed any bounded number of times, but all point-decision state is
\* unchanged.  Clearing terminalWitness here reproduces the production
\* liveness defect: rejectionAvailable is already consumed, so no transition
\* can restore the witness.
ObserveUnchangedCensus ==
    /\ phase = "active"
    /\ noopObservationsRemaining > 0
    /\ noopObservationsRemaining' = noopObservationsRemaining - 1
    /\ UNCHANGED <<phase, semanticChangesRemaining, rejectionAvailable,
                    terminalWitness, outcome>>

\* A changed view/policy/population epoch invalidates the old witness and
\* supplies one fresh decision.  The finite environment allowance keeps the
\* liveness claim honest without assuming that an endlessly changing view
\* must converge.
ChangeSemanticPopulation ==
    /\ phase = "active"
    /\ semanticChangesRemaining > 0
    /\ semanticChangesRemaining' = semanticChangesRemaining - 1
    /\ noopObservationsRemaining' = NoopObservationLimit
    /\ rejectionAvailable' = TRUE
    /\ terminalWitness' = FALSE
    /\ UNCHANGED <<phase, outcome>>

Finish ==
    /\ phase = "active"
    /\ terminalWitness
    /\ phase' = "terminal"
    /\ outcome' = "constrained"
    /\ UNCHANGED <<semanticChangesRemaining, noopObservationsRemaining,
                    rejectionAvailable, terminalWitness>>

Next ==
    \/ RejectFinerCandidate
    \/ ObserveUnchangedCensus
    \/ ChangeSemanticPopulation
    \/ Finish

Spec == Init /\ [][Next]_vars /\ WF_vars(Next)

\* Once the epoch's rejection has been consumed, its terminal witness cannot
\* disappear until a semantic population change supplies a fresh decision.
ConsumedRejectionHasWitness ==
    phase = "active" /\ ~rejectionAvailable => terminalWitness

WitnessOwnsConsumedRejection ==
    terminalWitness => ~rejectionAvailable

TerminalIsConstrainedAndWitnessed ==
    phase = "terminal" =>
        /\ outcome = "constrained"
        /\ terminalWitness
        /\ ~rejectionAvailable

EventuallyTerminal == <> (phase = "terminal")

=============================================================================
