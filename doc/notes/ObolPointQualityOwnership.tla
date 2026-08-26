--------------------- MODULE ObolPointQualityOwnership ---------------------
\* Focused ownership model for quiet small-occurrence calibration and retained
\* triangle recovery.  Calibration pauses retained admission while it waits
\* for an exact frame.  Recovery needs that admission producer.  They are one
\* sum-typed phase, never independent pending flags.

EXTENDS Naturals

Phases == {"idle", "calibration", "triangleRecovery"}
Producers == {"idle", "submission"}

VARIABLES phase, producer, framePending, recoveryPresented, pumpPending,
          terminal

vars == <<phase, producer, framePending, recoveryPresented, pumpPending,
          terminal>>

TypeOK ==
    /\ phase \in Phases
    /\ producer \in Producers
    /\ framePending \in BOOLEAN
    /\ recoveryPresented \in BOOLEAN
    /\ pumpPending \in BOOLEAN
    /\ terminal \in BOOLEAN

Init ==
    /\ phase = "idle"
    /\ producer = "idle"
    /\ framePending = FALSE
    /\ recoveryPresented = FALSE
    /\ pumpPending = FALSE
    /\ terminal = FALSE

StartCalibration ==
    /\ ~terminal
    /\ phase = "idle"
    /\ producer = "idle"
    /\ ~framePending
    /\ phase' = "calibration"
    /\ framePending' = TRUE
    /\ UNCHANGED <<producer, recoveryPresented, pumpPending, terminal>>

\* The exact calibration frame either settles the presentation or determines
\* that a coherent retained-prefix recovery is required.  The latter transfer
\* atomically retires calibration and enables its producer.
CompleteCalibration ==
    /\ phase = "calibration"
    /\ framePending
    /\ framePending' = FALSE
    /\ \/ /\ phase' = "idle"
          /\ terminal' = TRUE
          /\ UNCHANGED <<producer, recoveryPresented, pumpPending>>
       \/ /\ phase' = "triangleRecovery"
          /\ producer' = "submission"
          /\ recoveryPresented' = FALSE
          /\ pumpPending' = FALSE
          /\ UNCHANGED terminal

RunRecovery ==
    /\ phase = "triangleRecovery"
    /\ producer = "submission"
    /\ ~framePending
    /\ producer' = "idle"
    /\ framePending' = TRUE
    /\ UNCHANGED <<phase, recoveryPresented, pumpPending, terminal>>

CompleteRecoveryFrame ==
    /\ phase = "triangleRecovery"
    /\ producer = "idle"
    /\ framePending
    /\ phase' = "idle"
    /\ framePending' = FALSE
    /\ recoveryPresented' = FALSE
    /\ pumpPending' = FALSE
    /\ terminal' = TRUE
    /\ UNCHANGED producer

\* The recovery frame may also retire a stronger capacity or handoff barrier.
\* Its pixels remain valid evidence.  Preserve a level-triggered pump witness
\* rather than requiring a second frame merely to revisit the recovery owner.
DeferRecoveryCompletion ==
    /\ phase = "triangleRecovery"
    /\ producer = "idle"
    /\ framePending
    /\ framePending' = FALSE
    /\ recoveryPresented' = TRUE
    /\ pumpPending' = TRUE
    /\ UNCHANGED <<phase, producer, terminal>>

CompleteRecoveryPump ==
    /\ phase = "triangleRecovery"
    /\ producer = "idle"
    /\ ~framePending
    /\ recoveryPresented
    /\ pumpPending
    /\ phase' = "idle"
    /\ recoveryPresented' = FALSE
    /\ pumpPending' = FALSE
    /\ terminal' = TRUE
    /\ UNCHANGED <<producer, framePending>>

Next ==
    \/ StartCalibration
    \/ CompleteCalibration
    \/ RunRecovery
    \/ CompleteRecoveryFrame
    \/ DeferRecoveryCompletion
    \/ CompleteRecoveryPump

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(StartCalibration)
    /\ WF_vars(CompleteCalibration)
    /\ WF_vars(RunRecovery)
    /\ WF_vars(CompleteRecoveryFrame)
    /\ WF_vars(DeferRecoveryCompletion)
    /\ WF_vars(CompleteRecoveryPump)

\* Every active phase has exactly the kind of progress witness it can consume.
\* In particular, calibration can never coexist with the recovery producer.
PhaseHasWitness ==
    /\ (phase = "calibration" =>
           /\ producer = "idle"
           /\ framePending)
    /\ (phase = "triangleRecovery" =>
           \/ /\ producer = "submission"
                 /\ ~framePending
           \/ /\ producer = "idle"
                 /\ framePending
              \/ /\ producer = "idle"
                 /\ ~framePending
                 /\ recoveryPresented
                 /\ pumpPending)

TerminalIsQuiescent ==
    terminal =>
        /\ phase = "idle"
        /\ producer = "idle"
        /\ ~framePending
        /\ ~recoveryPresented
        /\ ~pumpPending

EventuallyTerminal == <>terminal

=============================================================================
