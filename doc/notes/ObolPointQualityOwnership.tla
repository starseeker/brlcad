--------------------- MODULE ObolPointQualityOwnership ---------------------
\* Focused ownership model for quiet small-occurrence calibration and retained
\* triangle recovery.  Calibration pauses retained admission while it waits
\* for an exact frame.  Recovery needs that admission producer.  They are one
\* sum-typed phase, never independent pending flags.

EXTENDS Naturals

Phases == {"idle", "adaptiveCalibration", "handoffConfirmation",
           "triangleRecovery"}
Producers == {"idle", "submission"}
PointCuts == {"retained", "onePixel"}

VARIABLES phase, producer, framePending, recoveryPresented, pumpPending,
          allocationCurrent, residualQualityDebt, capacityCalibration,
          staticTrial, pointCut, rejectionOccurred, terminal

vars == <<phase, producer, framePending, recoveryPresented, pumpPending,
          allocationCurrent, residualQualityDebt, capacityCalibration,
          staticTrial, pointCut, rejectionOccurred, terminal>>

TypeOK ==
    /\ phase \in Phases
    /\ producer \in Producers
    /\ framePending \in BOOLEAN
    /\ recoveryPresented \in BOOLEAN
    /\ pumpPending \in BOOLEAN
    /\ allocationCurrent \in BOOLEAN
    /\ residualQualityDebt \in BOOLEAN
    /\ capacityCalibration \in BOOLEAN
    /\ staticTrial \in BOOLEAN
    /\ pointCut \in PointCuts
    /\ rejectionOccurred \in BOOLEAN
    /\ terminal \in BOOLEAN

Init ==
    /\ phase = "idle"
    /\ producer = "idle"
    /\ framePending = FALSE
    /\ recoveryPresented = FALSE
    /\ pumpPending = FALSE
    /\ allocationCurrent = FALSE
    /\ residualQualityDebt = FALSE
    /\ capacityCalibration = FALSE
    /\ staticTrial = FALSE
    /\ pointCut = "retained"
    /\ rejectionOccurred = FALSE
    /\ terminal = FALSE

StartCalibration ==
    /\ ~terminal
    /\ phase = "idle"
    /\ producer = "idle"
    /\ ~framePending
    /\ ~capacityCalibration
    /\ phase' = "adaptiveCalibration"
    /\ framePending' = TRUE
    /\ UNCHANGED <<producer, recoveryPresented, pumpPending,
                    allocationCurrent, residualQualityDebt,
                    capacityCalibration, staticTrial, pointCut,
                    rejectionOccurred, terminal>>

\* A static point-quality trial replaces the completed retained point cut with
\* an exact one-pixel candidate.  It shares the calibration frame owner; a
\* deadline rejection must therefore release both owners atomically.
StartStaticCalibration ==
    /\ ~terminal
    /\ phase = "idle"
    /\ producer = "idle"
    /\ ~framePending
    /\ ~capacityCalibration
    /\ pointCut = "retained"
    /\ phase' = "adaptiveCalibration"
    /\ framePending' = TRUE
    /\ staticTrial' = TRUE
    /\ pointCut' = "onePixel"
    /\ UNCHANGED <<producer, recoveryPresented, pumpPending,
                    allocationCurrent, residualQualityDebt,
                    capacityCalibration, rejectionOccurred, terminal>>

StartCapacityCalibration ==
    /\ ~terminal
    /\ phase = "idle"
    /\ producer = "idle"
    /\ ~framePending
    /\ ~capacityCalibration
    /\ framePending' = TRUE
    /\ capacityCalibration' = TRUE
    /\ UNCHANGED <<phase, producer, recoveryPresented, pumpPending,
                    allocationCurrent, residualQualityDebt, staticTrial,
                    pointCut, rejectionOccurred, terminal>>

\* Point quality may change while a finite capacity sample is already in
\* flight.  The capacity owner consumes that frame, then leaves the same
\* level-triggered frame witness for point calibration.  No transition may
\* start a new capacity sample while either calibration phase owns it.
StartCalibrationDuringCapacity ==
    /\ ~terminal
    /\ phase = "idle"
    /\ producer = "idle"
    /\ framePending
    /\ capacityCalibration
    /\ phase' = "adaptiveCalibration"
    /\ UNCHANGED <<producer, framePending, recoveryPresented, pumpPending,
                    allocationCurrent, residualQualityDebt,
                    capacityCalibration, staticTrial, pointCut,
                    rejectionOccurred, terminal>>

\* A presentation handoff may supersede an adaptive classification frame.  It
\* retains that already scheduled framebuffer but changes its semantic owner:
\* the completed frame confirms the chosen cut and must not run adaptation.
StartHandoffConfirmation ==
    /\ ~terminal
    /\ phase \in {"idle", "adaptiveCalibration"}
    /\ producer = "idle"
    /\ ~capacityCalibration
    /\ ~staticTrial
    /\ (phase = "idle" => ~framePending)
    /\ phase' = "handoffConfirmation"
    /\ framePending' = TRUE
    /\ UNCHANGED <<producer, recoveryPresented, pumpPending,
                    allocationCurrent, residualQualityDebt,
                    capacityCalibration, staticTrial, pointCut,
                    rejectionOccurred, terminal>>

CompleteCapacityCalibration ==
    /\ phase = "idle"
    /\ producer = "idle"
    /\ framePending
    /\ capacityCalibration
    /\ framePending' = FALSE
    /\ capacityCalibration' = FALSE
    /\ terminal' = TRUE
    /\ UNCHANGED <<phase, producer, recoveryPresented, pumpPending,
                    allocationCurrent, residualQualityDebt, staticTrial,
                    pointCut, rejectionOccurred>>

CompleteCapacityBeforePoint ==
    /\ phase = "adaptiveCalibration"
    /\ producer = "idle"
    /\ framePending
    /\ capacityCalibration
    /\ capacityCalibration' = FALSE
    /\ UNCHANGED <<phase, producer, framePending, recoveryPresented,
                    pumpPending, allocationCurrent, residualQualityDebt,
                    staticTrial, pointCut, rejectionOccurred, terminal>>

\* The exact calibration frame either settles the presentation or determines
\* that a coherent retained-prefix recovery is required.  The latter transfer
\* atomically retires calibration and enables its producer.
CompleteCalibration ==
    /\ phase \in {"adaptiveCalibration", "handoffConfirmation"}
    /\ framePending
    /\ ~capacityCalibration
    /\ framePending' = FALSE
    /\ staticTrial' = FALSE
    /\ \/ /\ phase' = "idle"
          /\ terminal' = TRUE
          /\ allocationCurrent' = FALSE
          /\ residualQualityDebt' = FALSE
          /\ UNCHANGED <<producer, recoveryPresented, pumpPending,
                          capacityCalibration, pointCut, rejectionOccurred>>
       \/ /\ phase = "adaptiveCalibration"
          /\ phase' = "triangleRecovery"
          /\ producer' = "submission"
          /\ recoveryPresented' = FALSE
          /\ pumpPending' = FALSE
          /\ allocationCurrent' = FALSE
          /\ residualQualityDebt' = FALSE
          /\ UNCHANGED <<capacityCalibration, pointCut,
                          rejectionOccurred, terminal>>

RejectStaticCalibration ==
    /\ phase = "adaptiveCalibration"
    /\ producer = "idle"
    /\ framePending
    /\ staticTrial
    /\ pointCut = "onePixel"
    /\ phase' = "idle"
    /\ framePending' = FALSE
    /\ staticTrial' = FALSE
    /\ pointCut' = "retained"
    /\ rejectionOccurred' = TRUE
    /\ terminal' = TRUE
    /\ UNCHANGED <<producer, recoveryPresented, pumpPending,
                    allocationCurrent, residualQualityDebt,
                    capacityCalibration>>

RunRecovery ==
    /\ phase = "triangleRecovery"
    /\ producer = "submission"
    /\ ~framePending
    /\ producer' = "idle"
    /\ framePending' = TRUE
    /\ allocationCurrent' = TRUE
    /\ residualQualityDebt' = FALSE
    /\ UNCHANGED <<phase, recoveryPresented, pumpPending,
                    capacityCalibration, staticTrial, pointCut,
                    rejectionOccurred, terminal>>

\* A finite recovery allocation may cover the complete occurrence population
\* while deliberately leaving richer pixel-target demand.  That annotation is
\* residual quality debt, not another recovery producer.  Preserve a pump edge
\* so the owner can retire the phase without reopening the same scene scan.
CompleteNoOpRecoveryPass ==
    /\ phase = "triangleRecovery"
    /\ producer = "submission"
    /\ ~framePending
    /\ producer' = "idle"
    /\ recoveryPresented' = TRUE
    /\ pumpPending' = TRUE
    /\ allocationCurrent' = TRUE
    /\ residualQualityDebt' = TRUE
    /\ UNCHANGED <<phase, framePending, capacityCalibration, staticTrial,
                    pointCut, rejectionOccurred, terminal>>

CompleteRecoveryFrame ==
    /\ phase = "triangleRecovery"
    /\ producer = "idle"
    /\ framePending
    /\ phase' = "idle"
    /\ framePending' = FALSE
    /\ recoveryPresented' = FALSE
    /\ pumpPending' = FALSE
    /\ allocationCurrent' = FALSE
    /\ residualQualityDebt' = FALSE
    /\ terminal' = TRUE
    /\ UNCHANGED <<producer, capacityCalibration, staticTrial, pointCut,
                    rejectionOccurred>>

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
    /\ UNCHANGED <<phase, producer, allocationCurrent,
                    residualQualityDebt, capacityCalibration, staticTrial,
                    pointCut, rejectionOccurred, terminal>>

CompleteRecoveryPump ==
    /\ phase = "triangleRecovery"
    /\ producer = "idle"
    /\ ~framePending
    /\ recoveryPresented
    /\ pumpPending
    /\ (~residualQualityDebt \/ allocationCurrent)
    /\ phase' = "idle"
    /\ recoveryPresented' = FALSE
    /\ pumpPending' = FALSE
    /\ allocationCurrent' = FALSE
    /\ residualQualityDebt' = FALSE
    /\ terminal' = TRUE
    /\ UNCHANGED <<producer, framePending, capacityCalibration,
                    staticTrial, pointCut, rejectionOccurred>>

Next ==
    \/ StartCalibration
    \/ StartStaticCalibration
    \/ StartCapacityCalibration
    \/ StartCalibrationDuringCapacity
    \/ StartHandoffConfirmation
    \/ CompleteCapacityCalibration
    \/ CompleteCapacityBeforePoint
    \/ CompleteCalibration
    \/ RejectStaticCalibration
    \/ RunRecovery
    \/ CompleteNoOpRecoveryPass
    \/ CompleteRecoveryFrame
    \/ DeferRecoveryCompletion
    \/ CompleteRecoveryPump

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(StartCalibration)
    /\ WF_vars(StartStaticCalibration)
    /\ WF_vars(StartCapacityCalibration)
    /\ WF_vars(StartCalibrationDuringCapacity)
    /\ WF_vars(StartHandoffConfirmation)
    /\ WF_vars(CompleteCapacityCalibration)
    /\ WF_vars(CompleteCapacityBeforePoint)
    /\ WF_vars(CompleteCalibration)
    /\ WF_vars(RejectStaticCalibration)
    /\ WF_vars(RunRecovery)
    /\ WF_vars(CompleteNoOpRecoveryPass)
    /\ WF_vars(CompleteRecoveryFrame)
    /\ WF_vars(DeferRecoveryCompletion)
    /\ WF_vars(CompleteRecoveryPump)

\* Every active phase has exactly the kind of progress witness it can consume.
\* In particular, calibration can never coexist with the recovery producer.
PhaseHasWitness ==
    /\ (phase \in {"adaptiveCalibration", "handoffConfirmation"} =>
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

ResidualDebtHasCurrentAllocation ==
    (phase = "triangleRecovery" /\ producer = "idle" /\
     residualQualityDebt) => allocationCurrent

PointWaitsForExistingCapacity ==
    (phase \in {"adaptiveCalibration", "handoffConfirmation"} /\
     capacityCalibration) =>
        /\ producer = "idle"
        /\ framePending

HandoffConfirmationIsNotStatic ==
    (phase = "handoffConfirmation") => ~staticTrial

StaticTrialOwnsOnePixelFrame ==
    staticTrial =>
        /\ phase = "adaptiveCalibration"
        /\ framePending
        /\ pointCut = "onePixel"

RejectedTrialRestoresBaseline ==
    rejectionOccurred =>
        /\ terminal
        /\ pointCut = "retained"

TerminalIsQuiescent ==
    terminal =>
        /\ phase = "idle"
        /\ producer = "idle"
        /\ ~framePending
        /\ ~recoveryPresented
        /\ ~pumpPending
        /\ ~allocationCurrent
        /\ ~residualQualityDebt
        /\ ~capacityCalibration
        /\ ~staticTrial

EventuallyTerminal == <>terminal

=============================================================================
