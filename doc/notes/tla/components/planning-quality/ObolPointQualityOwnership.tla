--------------------- MODULE ObolPointQualityOwnership ---------------------
\* Focused ownership model for quiet small-occurrence calibration and retained
\* triangle recovery.  Calibration pauses retained admission while it waits
\* for an exact frame.  Recovery needs that admission producer.  They are one
\* sum-typed phase, never independent pending flags.

EXTENDS Naturals

Phases == {"idle", "adaptiveCalibration", "staticCalibration",
           "handoffConfirmation", "triangleRecovery"}
Producers == {"idle", "submission"}
PointCuts == {"retained", "onePixel"}

VARIABLES phase, producer, framePending, recoveryPresented, pumpPending,
          allocationCurrent, residualQualityDebt, capacityAllocation,
          capacityCalibration,
          capacityMeasurable, staticTrial, pointCut, rejectionOccurred,
          pointCandidate, terminal

vars == <<phase, producer, framePending, recoveryPresented, pumpPending,
          allocationCurrent, residualQualityDebt, capacityAllocation,
          capacityCalibration,
          capacityMeasurable, staticTrial, pointCut, rejectionOccurred,
          pointCandidate, terminal>>

TypeOK ==
    /\ phase \in Phases
    /\ producer \in Producers
    /\ framePending \in BOOLEAN
    /\ recoveryPresented \in BOOLEAN
    /\ pumpPending \in BOOLEAN
    /\ allocationCurrent \in BOOLEAN
    /\ residualQualityDebt \in BOOLEAN
    /\ capacityAllocation \in BOOLEAN
    /\ capacityCalibration \in BOOLEAN
    /\ capacityMeasurable \in BOOLEAN
    /\ staticTrial \in BOOLEAN
    /\ pointCut \in PointCuts
    /\ rejectionOccurred \in BOOLEAN
    /\ pointCandidate \in BOOLEAN
    /\ terminal \in BOOLEAN

Init ==
    /\ phase = "idle"
    /\ producer = "idle"
    /\ framePending = FALSE
    /\ recoveryPresented = FALSE
    /\ pumpPending = FALSE
    /\ allocationCurrent = FALSE
    /\ residualQualityDebt = FALSE
    /\ capacityAllocation = FALSE
    /\ capacityCalibration = FALSE
    /\ capacityMeasurable = TRUE
    /\ staticTrial = FALSE
    /\ pointCut = "retained"
    /\ rejectionOccurred = FALSE
    /\ pointCandidate = FALSE
    /\ terminal = FALSE

\* Projected demand may identify a terminal point-quality opportunity before
\* the stronger capacity transaction has retired.  Discovery records planning
\* input only: it neither owns a framebuffer nor pauses the capacity producer.
DiscoverPointCandidate ==
    /\ ~terminal
    /\ phase = "idle"
    /\ ~pointCandidate
    /\ capacityCalibration
    /\ pointCandidate' = TRUE
    /\ UNCHANGED <<phase, producer, framePending, recoveryPresented,
                    pumpPending, allocationCurrent, residualQualityDebt,
                    capacityAllocation, capacityCalibration,
                    capacityMeasurable, staticTrial, pointCut,
                    rejectionOccurred, terminal>>

\* Only this reducer transition turns a derived candidate into a typed point
\* presentation obligation.  Stronger capacity work must have retired first.
CommitPointCandidate ==
    /\ ~terminal
    /\ phase = "idle"
    /\ pointCandidate
    /\ producer = "idle"
    /\ ~framePending
    /\ ~capacityAllocation
    /\ ~capacityCalibration
    /\ phase' = "adaptiveCalibration"
    /\ framePending' = TRUE
    /\ pointCandidate' = FALSE
    /\ UNCHANGED <<producer, recoveryPresented, pumpPending,
                    allocationCurrent, residualQualityDebt,
                    capacityAllocation, capacityCalibration,
                    capacityMeasurable, staticTrial, pointCut,
                    rejectionOccurred, terminal>>

StartCalibration ==
    /\ ~terminal
    /\ phase = "idle"
    /\ ~pointCandidate
    /\ producer = "idle"
    /\ ~framePending
    /\ ~capacityCalibration
    /\ phase' = "adaptiveCalibration"
    /\ framePending' = TRUE
    /\ UNCHANGED <<producer, recoveryPresented, pumpPending,
                    allocationCurrent, residualQualityDebt,
                    capacityAllocation, capacityCalibration,
                    capacityMeasurable, staticTrial, pointCut,
                    rejectionOccurred, pointCandidate, terminal>>

\* A candidate selected from a retained event-driven framebuffer keeps the
\* static deadline through every exact point-bracket successor.  Ordinary
\* calibration requests cannot downgrade this owner to the streaming cadence.
StartStaticPointCalibration ==
    /\ ~terminal
    /\ phase = "idle"
    /\ ~pointCandidate
    /\ producer = "idle"
    /\ ~framePending
    /\ ~capacityCalibration
    /\ phase' = "staticCalibration"
    /\ framePending' = TRUE
    /\ UNCHANGED <<producer, recoveryPresented, pumpPending,
                    allocationCurrent, residualQualityDebt,
                    capacityAllocation, capacityCalibration,
                    capacityMeasurable, staticTrial, pointCut,
                    rejectionOccurred, pointCandidate, terminal>>

\* A static point-quality trial replaces the completed retained point cut with
\* an exact one-pixel candidate.  It shares the calibration frame owner; a
\* deadline rejection must therefore release both owners atomically.
StartStaticCalibration ==
    /\ ~terminal
    /\ phase = "idle"
    /\ ~pointCandidate
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
                    capacityAllocation, capacityCalibration,
                    capacityMeasurable,
                    rejectionOccurred, pointCandidate, terminal>>

StartCapacityCalibration ==
    /\ ~terminal
    /\ phase = "idle"
    /\ ~pointCandidate
    /\ producer = "idle"
    /\ ~framePending
    /\ ~capacityCalibration
    /\ producer' = "submission"
    /\ capacityAllocation' = TRUE
    /\ capacityCalibration' = TRUE
    /\ capacityMeasurable' = TRUE
    /\ UNCHANGED <<phase, framePending, recoveryPresented, pumpPending,
                    allocationCurrent, residualQualityDebt, staticTrial,
                    pointCut, rejectionOccurred, pointCandidate, terminal>>

\* The capacity candidate must be allocated before either capacity or point
\* calibration can own its exact presentation.  A point phase queued during
\* this transition is a successor, not permission to pause the sole
\* submission producer.
CompleteCapacityAllocation ==
    /\ ~terminal
    /\ capacityCalibration
    /\ capacityAllocation
    /\ producer = "submission"
    /\ ~framePending
    /\ producer' = "idle"
    /\ capacityAllocation' = FALSE
    /\ framePending' = TRUE
    /\ allocationCurrent' =
        IF phase = "triangleRecovery" THEN TRUE ELSE allocationCurrent
    /\ UNCHANGED <<phase, recoveryPresented, pumpPending,
                    residualQualityDebt,
                    capacityCalibration, capacityMeasurable, staticTrial,
                    pointCut, rejectionOccurred, pointCandidate, terminal>>

\* A recovery pass may expose a bounded capacity candidate while it is
\* allocating the coherent retained population.  That candidate is the
\* recovery transaction's stronger planning owner.  It uses the existing
\* submission producer and must complete before recovery consumes a frame.
StartRecoveryCapacityPlanning ==
    /\ ~terminal
    /\ phase = "triangleRecovery"
    /\ producer = "submission"
    /\ ~framePending
    /\ ~capacityAllocation
    /\ ~capacityCalibration
    /\ capacityAllocation' = TRUE
    /\ capacityCalibration' = TRUE
    /\ capacityMeasurable' = TRUE
    /\ UNCHANGED <<phase, producer, framePending, recoveryPresented,
                    pumpPending, allocationCurrent, residualQualityDebt,
                    staticTrial, pointCut, rejectionOccurred,
                    pointCandidate, terminal>>

\* Point quality may change while a finite capacity sample is already in
\* flight.  The capacity owner consumes that frame, then leaves the same
\* level-triggered frame witness for point calibration.  No transition may
\* start a new capacity sample while either calibration phase owns it.
StartCalibrationDuringCapacity ==
    /\ ~terminal
    /\ phase = "idle"
    /\ capacityCalibration
    /\ \/ /\ capacityAllocation
           /\ producer = "submission"
           /\ ~framePending
       \/ /\ ~capacityAllocation
           /\ producer = "idle"
           /\ framePending
    /\ phase' = "adaptiveCalibration"
    /\ pointCandidate' = FALSE
    /\ capacityMeasurable' \in BOOLEAN
    /\ UNCHANGED <<producer, framePending, recoveryPresented, pumpPending,
                    allocationCurrent, residualQualityDebt,
                    capacityAllocation, capacityCalibration, staticTrial,
                    pointCut,
                    rejectionOccurred, terminal>>

\* A presentation handoff may supersede an adaptive classification frame.  It
\* retains that already scheduled framebuffer but changes its semantic owner:
\* the completed frame confirms the chosen cut and must not run adaptation.
StartHandoffConfirmation ==
    /\ ~terminal
    /\ phase \in {"idle", "adaptiveCalibration", "staticCalibration"}
    /\ producer = "idle"
    /\ ~capacityCalibration
    /\ ~staticTrial
    /\ ~pointCandidate
    /\ (phase = "idle" => ~framePending)
    /\ phase' = "handoffConfirmation"
    /\ framePending' = TRUE
    /\ UNCHANGED <<producer, recoveryPresented, pumpPending,
                    allocationCurrent, residualQualityDebt,
                    capacityAllocation, capacityCalibration,
                    capacityMeasurable, staticTrial, pointCut,
                    rejectionOccurred, pointCandidate, terminal>>

CompleteCapacityCalibration ==
    /\ phase = "idle"
    /\ producer = "idle"
    /\ framePending
    /\ ~capacityAllocation
    /\ capacityCalibration
    /\ capacityMeasurable
    /\ framePending' = FALSE
    /\ capacityCalibration' = FALSE
    /\ terminal' = ~pointCandidate
    /\ UNCHANGED <<phase, producer, recoveryPresented, pumpPending,
                    allocationCurrent, residualQualityDebt,
                    capacityAllocation, staticTrial,
                    pointCut, rejectionOccurred, capacityMeasurable,
                    pointCandidate>>

CompleteCapacityBeforePoint ==
    /\ phase = "adaptiveCalibration"
    /\ producer = "idle"
    /\ framePending
    /\ capacityCalibration
    /\ ~capacityAllocation
    /\ capacityMeasurable
    /\ capacityCalibration' = FALSE
    /\ UNCHANGED <<phase, producer, framePending, recoveryPresented,
                    pumpPending, allocationCurrent, residualQualityDebt,
                    capacityAllocation,
                    capacityMeasurable, staticTrial, pointCut,
                    rejectionOccurred, pointCandidate, terminal>>

\* The exact capacity frame belongs to the stronger transaction even while
\* triangle recovery is pending.  Once capacity consumes its evidence, the
\* same completed framebuffer is available to retire recovery.
CompleteCapacityBeforeRecovery ==
    /\ phase = "triangleRecovery"
    /\ producer = "idle"
    /\ framePending
    /\ capacityCalibration
    /\ ~capacityAllocation
    /\ capacityMeasurable
    /\ capacityCalibration' = FALSE
    /\ UNCHANGED <<phase, producer, framePending, recoveryPresented,
                    pumpPending, allocationCurrent, residualQualityDebt,
                    capacityAllocation, capacityMeasurable, staticTrial,
                    pointCut, rejectionOccurred, pointCandidate, terminal>>

\* A complete structural-placeholder frame cannot become capacity evidence by
\* repainting the same population: capacity timing intentionally excludes the
\* placeholders.  Retire that measurement owner and leave the frame witness
\* for point classification (or the subsequent structural-repair owner).
YieldUnmeasurableCapacityToPoint ==
    /\ phase = "adaptiveCalibration"
    /\ producer = "idle"
    /\ framePending
    /\ capacityCalibration
    /\ ~capacityAllocation
    /\ ~capacityMeasurable
    /\ capacityCalibration' = FALSE
    /\ UNCHANGED <<phase, producer, framePending, recoveryPresented,
                    pumpPending, allocationCurrent, residualQualityDebt,
                    capacityAllocation,
                    capacityMeasurable, staticTrial, pointCut,
                    rejectionOccurred, pointCandidate, terminal>>

\* The exact calibration frame either settles the presentation or determines
\* that a coherent retained-prefix recovery is required.  The latter transfer
\* atomically retires calibration and enables its producer.
CompleteCalibration ==
    /\ phase \in {"adaptiveCalibration", "staticCalibration",
                    "handoffConfirmation"}
    /\ framePending
    /\ ~capacityCalibration
    /\ framePending' = FALSE
    /\ staticTrial' = FALSE
    /\ \/ /\ phase' = "idle"
          /\ terminal' = TRUE
          /\ allocationCurrent' = FALSE
          /\ residualQualityDebt' = FALSE
          /\ UNCHANGED <<producer, recoveryPresented, pumpPending,
                          capacityAllocation, capacityCalibration,
                          capacityMeasurable, pointCut,
                          rejectionOccurred, pointCandidate>>
       \/ /\ phase = "adaptiveCalibration"
          /\ phase' = "triangleRecovery"
          /\ producer' = "submission"
          /\ recoveryPresented' = FALSE
          /\ pumpPending' = FALSE
          /\ allocationCurrent' = FALSE
          /\ residualQualityDebt' = FALSE
          /\ UNCHANGED <<capacityAllocation, capacityCalibration,
                          capacityMeasurable, pointCut,
                          rejectionOccurred, pointCandidate, terminal>>

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
                    capacityAllocation, capacityCalibration,
                    capacityMeasurable, pointCandidate>>

RunRecovery ==
    /\ phase = "triangleRecovery"
    /\ producer = "submission"
    /\ ~framePending
    /\ ~capacityAllocation
    /\ ~capacityCalibration
    /\ producer' = "idle"
    /\ framePending' = TRUE
    /\ allocationCurrent' = TRUE
    /\ residualQualityDebt' = FALSE
    /\ UNCHANGED <<phase, recoveryPresented, pumpPending,
                    capacityAllocation, capacityCalibration,
                    capacityMeasurable, staticTrial, pointCut,
                    rejectionOccurred, pointCandidate, terminal>>

\* A finite recovery allocation may cover the complete occurrence population
\* while deliberately leaving richer pixel-target demand.  That annotation is
\* residual quality debt, not another recovery producer.  Preserve a pump edge
\* so the owner can retire the phase without reopening the same scene scan.
CompleteNoOpRecoveryPass ==
    /\ phase = "triangleRecovery"
    /\ producer = "submission"
    /\ ~framePending
    /\ ~capacityAllocation
    /\ ~capacityCalibration
    /\ producer' = "idle"
    /\ recoveryPresented' = TRUE
    /\ pumpPending' = TRUE
    /\ allocationCurrent' = TRUE
    /\ residualQualityDebt' = TRUE
    /\ UNCHANGED <<phase, framePending, capacityAllocation,
                    capacityCalibration,
                    capacityMeasurable, staticTrial, pointCut,
                    rejectionOccurred, pointCandidate, terminal>>

CompleteRecoveryFrame ==
    /\ phase = "triangleRecovery"
    /\ producer = "idle"
    /\ framePending
    /\ ~capacityAllocation
    /\ ~capacityCalibration
    /\ phase' = "idle"
    /\ framePending' = FALSE
    /\ recoveryPresented' = FALSE
    /\ pumpPending' = FALSE
    /\ allocationCurrent' = FALSE
    /\ residualQualityDebt' = FALSE
    /\ terminal' = TRUE
    /\ UNCHANGED <<producer, capacityAllocation, capacityCalibration,
                    capacityMeasurable,
                    staticTrial, pointCut, rejectionOccurred,
                    pointCandidate>>

\* The recovery frame may also retire a stronger capacity or handoff barrier.
\* Its pixels remain valid evidence.  Preserve a level-triggered pump witness
\* rather than requiring a second frame merely to revisit the recovery owner.
DeferRecoveryCompletion ==
    /\ phase = "triangleRecovery"
    /\ producer = "idle"
    /\ framePending
    /\ ~capacityAllocation
    /\ ~capacityCalibration
    /\ framePending' = FALSE
    /\ recoveryPresented' = TRUE
    /\ pumpPending' = TRUE
    /\ UNCHANGED <<phase, producer, allocationCurrent,
                    residualQualityDebt, capacityAllocation,
                    capacityCalibration,
                    capacityMeasurable, staticTrial, pointCut,
                    rejectionOccurred, pointCandidate, terminal>>

CompleteRecoveryPump ==
    /\ phase = "triangleRecovery"
    /\ producer = "idle"
    /\ ~framePending
    /\ recoveryPresented
    /\ pumpPending
    /\ ~capacityAllocation
    /\ ~capacityCalibration
    /\ (~residualQualityDebt \/ allocationCurrent)
    /\ phase' = "idle"
    /\ recoveryPresented' = FALSE
    /\ pumpPending' = FALSE
    /\ allocationCurrent' = FALSE
    /\ residualQualityDebt' = FALSE
    /\ terminal' = TRUE
    /\ UNCHANGED <<producer, framePending, capacityAllocation,
                    capacityCalibration,
                    capacityMeasurable, staticTrial, pointCut,
                    rejectionOccurred, pointCandidate>>

Next ==
    \/ DiscoverPointCandidate
    \/ CommitPointCandidate
    \/ StartCalibration
    \/ StartStaticPointCalibration
    \/ StartStaticCalibration
    \/ StartCapacityCalibration
    \/ CompleteCapacityAllocation
    \/ StartRecoveryCapacityPlanning
    \/ StartCalibrationDuringCapacity
    \/ StartHandoffConfirmation
    \/ CompleteCapacityCalibration
    \/ CompleteCapacityBeforePoint
    \/ CompleteCapacityBeforeRecovery
    \/ YieldUnmeasurableCapacityToPoint
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
    /\ WF_vars(DiscoverPointCandidate)
    /\ WF_vars(CommitPointCandidate)
    /\ WF_vars(StartCalibration)
    /\ WF_vars(StartStaticPointCalibration)
    /\ WF_vars(StartStaticCalibration)
    /\ WF_vars(StartCapacityCalibration)
    /\ WF_vars(CompleteCapacityAllocation)
    /\ WF_vars(StartRecoveryCapacityPlanning)
    /\ WF_vars(StartCalibrationDuringCapacity)
    /\ WF_vars(StartHandoffConfirmation)
    /\ WF_vars(CompleteCapacityCalibration)
    /\ WF_vars(CompleteCapacityBeforePoint)
    /\ WF_vars(CompleteCapacityBeforeRecovery)
    /\ WF_vars(YieldUnmeasurableCapacityToPoint)
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
    /\ (phase \in {"adaptiveCalibration", "staticCalibration",
                     "handoffConfirmation"} =>
           \/ /\ capacityAllocation
                 /\ producer = "submission"
                 /\ ~framePending
              \/ /\ ~capacityAllocation
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
    (phase \in {"adaptiveCalibration", "staticCalibration",
                 "handoffConfirmation"} /\
     capacityCalibration) =>
        \/ /\ capacityAllocation
              /\ producer = "submission"
              /\ ~framePending
           \/ /\ ~capacityAllocation
              /\ producer = "idle"
              /\ framePending

RecoveryCapacityHasPriority ==
    (phase = "triangleRecovery" /\ capacityCalibration) =>
        \/ /\ capacityAllocation
              /\ producer = "submission"
              /\ ~framePending
           \/ /\ ~capacityAllocation
              /\ producer = "idle"
              /\ framePending

\* A discovered opportunity is not a presentation owner.  If some frame is
\* pending before the reducer commits the candidate, that frame is evidence
\* for the stronger capacity transaction.
CandidateFrameBelongsToCapacity ==
    (pointCandidate /\ framePending) =>
        /\ phase = "idle"
        /\ capacityCalibration

HandoffConfirmationIsNotStatic ==
    (phase = "handoffConfirmation") => ~staticTrial

StaticPointOwnsFrame ==
    (phase = "staticCalibration") =>
        /\ framePending
        /\ producer = "idle"
        /\ ~capacityAllocation
        /\ ~capacityCalibration

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
        /\ ~capacityAllocation
        /\ ~capacityCalibration
        /\ ~staticTrial
        /\ ~pointCandidate

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => terminal

EventuallyTerminal == <>terminal

=============================================================================
