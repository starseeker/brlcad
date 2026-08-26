------------------- MODULE ObolPresentationPreparation -------------------
\* Finite-work contract for renderer-side presentation preparation.
\*
\* A changed diagnostic serial is not a progress witness.  For one exact
\* revision-bound target, progress means that a finite remaining-work measure
\* strictly decreases while its admitted transient reservation stays bounded.
\* Superseding demand cancels the old obligation.  An interrupted attempt
\* which neither completes nor decreases remaining work must return a typed
\* capacity constraint; it cannot request the same frame indefinitely.

EXTENDS Naturals, TLC

CONSTANT MaxRevision, MaxUnits, MaxReservation

PreparationStates == {"idle", "preparing", "complete",
                      "constrained", "failed"}
TerminalStates == {"complete", "constrained", "failed"}

VARIABLES demandRevision,
          inputClosed,
          state,
          targetRevision,
          totalUnits,
          remainingUnits,
          reservation,
          attemptStartRemaining,
          attemptInFlight,
          committedRevision

vars == <<demandRevision, inputClosed, state, targetRevision, totalUnits,
          remainingUnits, reservation, attemptStartRemaining,
          attemptInFlight, committedRevision>>

TypeOK ==
    /\ demandRevision \in 1..MaxRevision
    /\ inputClosed \in BOOLEAN
    /\ state \in PreparationStates
    /\ targetRevision \in 0..MaxRevision
    /\ totalUnits \in 0..MaxUnits
    /\ remainingUnits \in 0..MaxUnits
    /\ remainingUnits <= totalUnits
    /\ reservation \in 0..MaxReservation
    /\ attemptStartRemaining \in 0..MaxUnits
    /\ attemptInFlight \in BOOLEAN
    /\ committedRevision \in 0..MaxRevision

Init ==
    /\ demandRevision = 1
    /\ inputClosed = FALSE
    /\ state = "idle"
    /\ targetRevision = 0
    /\ totalUnits = 0
    /\ remainingUnits = 0
    /\ reservation = 0
    /\ attemptStartRemaining = 0
    /\ attemptInFlight = FALSE
    /\ committedRevision = 0

PublishDemand ==
    /\ ~inputClosed
    /\ demandRevision < MaxRevision
    /\ demandRevision' = demandRevision + 1
    /\ state' = "idle"
    /\ targetRevision' = 0
    /\ totalUnits' = 0
    /\ remainingUnits' = 0
    /\ reservation' = 0
    /\ attemptStartRemaining' = 0
    /\ attemptInFlight' = FALSE
    /\ UNCHANGED <<inputClosed, committedRevision>>

CloseInput ==
    /\ ~inputClosed
    /\ inputClosed' = TRUE
    /\ UNCHANGED <<demandRevision, state, targetRevision, totalUnits,
                    remainingUnits, reservation, attemptStartRemaining,
                    attemptInFlight, committedRevision>>

\* Admission chooses the finite work and reservation for the exact target.
\* Zero work is already complete; a nonzero target must own memory.
AdmitPreparation(units, bytes) ==
    /\ state = "idle"
    /\ units \in 0..MaxUnits
    /\ bytes \in 0..MaxReservation
    /\ (units = 0) \/ (bytes > 0)
    /\ targetRevision' = demandRevision
    /\ totalUnits' = units
    /\ remainingUnits' = units
    /\ reservation' = bytes
    /\ state' = IF units = 0 THEN "complete" ELSE "preparing"
    /\ attemptStartRemaining' = 0
    /\ attemptInFlight' = FALSE
    /\ UNCHANGED <<demandRevision, inputClosed, committedRevision>>

RejectAdmission ==
    /\ state = "idle"
    /\ targetRevision' = demandRevision
    /\ state' = "constrained"
    /\ totalUnits' = 0
    /\ remainingUnits' = 0
    /\ reservation' = 0
    /\ attemptStartRemaining' = 0
    /\ attemptInFlight' = FALSE
    /\ UNCHANGED <<demandRevision, inputClosed, committedRevision>>

StartPreparation ==
    \E units \in 0..MaxUnits, bytes \in 0..MaxReservation:
        AdmitPreparation(units, bytes)

BeginAttempt ==
    /\ state = "preparing"
    /\ targetRevision = demandRevision
    /\ ~attemptInFlight
    /\ remainingUnits > 0
    /\ attemptInFlight' = TRUE
    /\ attemptStartRemaining' = remainingUnits
    /\ UNCHANGED <<demandRevision, inputClosed, state, targetRevision,
                    totalUnits, remainingUnits, reservation,
                    committedRevision>>

\* A unit is immutable and complete before it decrements the ledger.
PrepareUnit ==
    /\ state = "preparing"
    /\ targetRevision = demandRevision
    /\ attemptInFlight
    /\ remainingUnits > 0
    /\ remainingUnits' = remainingUnits - 1
    /\ UNCHANGED <<demandRevision, inputClosed, state, targetRevision,
                    totalUnits, reservation, attemptStartRemaining,
                    attemptInFlight, committedRevision>>

CompleteAttempt ==
    /\ state = "preparing"
    /\ attemptInFlight
    /\ remainingUnits = 0
    /\ state' = "complete"
    /\ attemptInFlight' = FALSE
    /\ attemptStartRemaining' = 0
    /\ UNCHANGED <<demandRevision, inputClosed, targetRevision, totalUnits,
                    remainingUnits, reservation, committedRevision>>

\* An interrupted attempt may retry only when its exact finite measure
\* decreased.  This is the contract a renderer preparation certificate must
\* establish; incrementing an activity serial is insufficient.
InterruptWithProgress ==
    /\ state = "preparing"
    /\ attemptInFlight
    /\ remainingUnits > 0
    /\ remainingUnits < attemptStartRemaining
    /\ attemptInFlight' = FALSE
    /\ attemptStartRemaining' = 0
    /\ UNCHANGED <<demandRevision, inputClosed, state, targetRevision,
                    totalUnits, remainingUnits, reservation,
                    committedRevision>>

\* A no-progress interruption terminally constrains this target.  A later
\* demand/capacity revision may supersede it through PublishDemand.
InterruptWithoutProgress ==
    /\ state = "preparing"
    /\ attemptInFlight
    /\ remainingUnits = attemptStartRemaining
    /\ state' = "constrained"
    /\ attemptInFlight' = FALSE
    /\ attemptStartRemaining' = 0
    /\ UNCHANGED <<demandRevision, inputClosed, targetRevision, totalUnits,
                    remainingUnits, reservation, committedRevision>>

FailPreparation ==
    /\ state = "preparing"
    /\ state' = "failed"
    /\ attemptInFlight' = FALSE
    /\ attemptStartRemaining' = 0
    /\ UNCHANGED <<demandRevision, inputClosed, targetRevision, totalUnits,
                    remainingUnits, reservation, committedRevision>>

Commit ==
    /\ state = "complete"
    /\ targetRevision = demandRevision
    /\ committedRevision' = targetRevision
    /\ reservation' = 0
    /\ UNCHANGED <<demandRevision, inputClosed, state, targetRevision,
                    totalUnits, remainingUnits, attemptStartRemaining,
                    attemptInFlight>>

Next ==
    \/ PublishDemand
    \/ CloseInput
    \/ RejectAdmission
    \/ StartPreparation
    \/ BeginAttempt
    \/ PrepareUnit
    \/ CompleteAttempt
    \/ InterruptWithProgress
    \/ InterruptWithoutProgress
    \/ FailPreparation
    \/ Commit

Spec == Init /\ [][Next]_vars
        /\ WF_vars(StartPreparation)
        /\ WF_vars(RejectAdmission)
        /\ WF_vars(BeginAttempt)
        /\ WF_vars(PrepareUnit)
        /\ WF_vars(CompleteAttempt)
        /\ WF_vars(InterruptWithProgress)
        /\ WF_vars(InterruptWithoutProgress)
        /\ WF_vars(Commit)

TargetIsRevisionBound ==
    state # "idle" => targetRevision > 0

ActiveWorkHasReservation ==
    state = "preparing" => reservation > 0

CompletedWorkIsFinite ==
    state = "complete" => remainingUnits = 0

NoStaleCommit ==
    committedRevision > 0 => committedRevision <= demandRevision

EventuallyTerminalAfterFinalDemand ==
    [](inputClosed => <>(state \in TerminalStates))

=============================================================================
