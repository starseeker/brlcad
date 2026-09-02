--------------------------- MODULE ObolStaticQuality -----------------------
\* Control-plane model for quiet-view quality admission.
\*
\* A rich static overscan population and the protected visual floor are
\* distinct requests.  Rejecting overscan is not evidence that the cheaper
\* floor is unaffordable.  The floor itself is atomic and may still be too
\* expensive even when one smaller occurrence-local improvement fits.  The
\* production controller must then keep the last complete framebuffer, enter
\* its marginal path, and either present bounded improvements or retain an
\* explicit capacity witness.  It must not restart ordinary convergence in the
\* same quiet epoch.
\*
\* Quality units stand for ordered renderer cuts, not equal triangle counts or
\* equal elapsed time.  TLC is checking ownership, preservation, termination,
\* and liveness; projection and performance heuristics remain executable tests.

EXTENDS Naturals, TLC

CONSTANTS MaxInputEpoch, MaxQuality

Phases == {"input", "ordinary", "overscan", "floor", "marginal", "present", "reconcile",
           "ready", "accepted", "constrained"}
Origins == {"none", "overscan", "floor", "marginal"}

VARIABLES inputEpoch,
          inputOpen,
          capacity,
          phase,
          committedQuality,
          presentedQuality,
          candidateQuality,
          candidateOrigin,
          overscanRejected,
          protectedFloorRejected,
          nextCutRejected,
          globalCeiling,
          terminalCeiling,
          staticStarts

vars == <<inputEpoch, inputOpen, capacity, phase, committedQuality,
          presentedQuality, candidateQuality, candidateOrigin,
          overscanRejected, protectedFloorRejected, nextCutRejected, globalCeiling,
          terminalCeiling, staticStarts>>

Terminal == phase \in {"ready", "accepted", "constrained"}

FloorTarget ==
    IF committedQuality + 2 > MaxQuality
    THEN MaxQuality
    ELSE committedQuality + 2

TypeOK ==
    /\ inputEpoch \in 0..MaxInputEpoch
    /\ inputOpen \in BOOLEAN
    /\ capacity \in 0..MaxQuality
    /\ phase \in Phases
    /\ committedQuality \in 0..MaxQuality
    /\ presentedQuality \in 0..MaxQuality
    /\ candidateQuality \in 0..MaxQuality
    /\ candidateOrigin \in Origins
    /\ overscanRejected \in BOOLEAN
    /\ protectedFloorRejected \in BOOLEAN
    /\ nextCutRejected \in BOOLEAN
    /\ globalCeiling \in BOOLEAN
    /\ terminalCeiling \in BOOLEAN
    /\ staticStarts \in 0..1

Init ==
    /\ inputEpoch = 0
    /\ inputOpen = FALSE
    /\ capacity \in 0..MaxQuality
    /\ phase = "ordinary"
    /\ committedQuality = 0
    /\ presentedQuality = 0
    /\ candidateQuality = 0
    /\ candidateOrigin = "none"
    /\ overscanRejected = FALSE
    /\ protectedFloorRejected = FALSE
    /\ nextCutRejected = FALSE
    /\ globalCeiling = FALSE
    /\ terminalCeiling = FALSE
    /\ staticStarts = 0

\* Input invalidates view-local trials, but cumulative resident/presented
\* quality remains a legal starting point for the next view epoch.
BeginInput ==
    /\ ~inputOpen
    /\ inputEpoch < MaxInputEpoch
    /\ inputEpoch' = inputEpoch + 1
    /\ inputOpen' = TRUE
    /\ phase' = "input"
    /\ candidateQuality' = 0
    /\ candidateOrigin' = "none"
    /\ overscanRejected' = FALSE
    /\ protectedFloorRejected' = FALSE
    /\ nextCutRejected' = FALSE
    /\ globalCeiling' = FALSE
    /\ terminalCeiling' = FALSE
    /\ staticStarts' = 0
    /\ UNCHANGED <<capacity, committedQuality, presentedQuality>>

EndInput ==
    /\ inputOpen
    /\ inputOpen' = FALSE
    /\ phase' = "ordinary"
    /\ UNCHANGED <<inputEpoch, capacity, committedQuality,
                    presentedQuality, candidateQuality, candidateOrigin,
                    overscanRejected, protectedFloorRejected, nextCutRejected, globalCeiling,
                    terminalCeiling, staticStarts>>

StartStatic ==
    /\ ~inputOpen
    /\ phase = "ordinary"
    /\ committedQuality < MaxQuality
    /\ phase' = "overscan"
    /\ staticStarts' = staticStarts + 1
    /\ UNCHANGED <<inputEpoch, capacity, committedQuality,
                    presentedQuality, candidateQuality, candidateOrigin,
                    overscanRejected, protectedFloorRejected, nextCutRejected, globalCeiling,
                    terminalCeiling, inputOpen>>

FinishAlreadyRich ==
    /\ ~inputOpen
    /\ phase = "ordinary"
    /\ committedQuality = MaxQuality
    /\ phase' = "ready"
    /\ UNCHANGED <<inputEpoch, capacity, committedQuality,
                    presentedQuality, candidateQuality, candidateOrigin,
                    overscanRejected, protectedFloorRejected, nextCutRejected, staticStarts,
                    globalCeiling, terminalCeiling, inputOpen>>

\* The richest static population is a separate opportunistic candidate.  Its
\* rejection advances to the cheaper protected floor without changing the
\* floor's evidence state.
AdmitStaticOverscan ==
    /\ phase = "overscan"
    /\ capacity >= MaxQuality
    /\ candidateQuality' = MaxQuality
    /\ candidateOrigin' = "overscan"
    /\ phase' = "present"
    /\ UNCHANGED <<inputEpoch, inputOpen, capacity, committedQuality,
                    presentedQuality, overscanRejected,
                    protectedFloorRejected, nextCutRejected, globalCeiling,
                    terminalCeiling, staticStarts>>

RejectStaticOverscan ==
    /\ phase = "overscan"
    /\ capacity < MaxQuality
    /\ overscanRejected' = TRUE
    /\ phase' = "floor"
    /\ UNCHANGED <<inputEpoch, inputOpen, capacity, committedQuality,
                    presentedQuality, candidateQuality, candidateOrigin,
                    protectedFloorRejected, nextCutRejected, globalCeiling,
                    terminalCeiling, staticStarts>>

AdmitProtectedFloor ==
    /\ phase = "floor"
    /\ capacity >= FloorTarget
    /\ candidateQuality' = FloorTarget
    /\ candidateOrigin' = "floor"
    /\ phase' = "present"
    /\ UNCHANGED <<inputEpoch, inputOpen, capacity, committedQuality,
                    presentedQuality, overscanRejected, protectedFloorRejected,
                    nextCutRejected, globalCeiling, staticStarts>>
    /\ UNCHANGED terminalCeiling

RejectProtectedFloor ==
    /\ phase = "floor"
    /\ capacity < FloorTarget
    /\ protectedFloorRejected' = TRUE
    /\ phase' = "marginal"
    /\ UNCHANGED <<inputEpoch, inputOpen, capacity, committedQuality,
                    presentedQuality, candidateQuality, candidateOrigin,
                    overscanRejected, nextCutRejected, globalCeiling, staticStarts>>
    /\ UNCHANGED terminalCeiling

\* This is the critical fallback.  It must remain enabled while the static
\* phase is active even though the complete protected floor did not fit.
AdmitMarginal ==
    /\ phase = "marginal"
    /\ committedQuality < MaxQuality
    /\ capacity >= committedQuality + 1
    /\ candidateQuality' = committedQuality + 1
    /\ candidateOrigin' = "marginal"
    /\ phase' = "present"
    /\ UNCHANGED <<inputEpoch, inputOpen, capacity, committedQuality,
                    presentedQuality, overscanRejected, protectedFloorRejected,
                    nextCutRejected, globalCeiling, staticStarts>>
    /\ UNCHANGED terminalCeiling

RejectMarginal ==
    /\ phase = "marginal"
    /\ committedQuality < MaxQuality
    /\ capacity < committedQuality + 1
    /\ nextCutRejected' = TRUE
    /\ globalCeiling' = TRUE
    /\ phase' = "reconcile"
    /\ UNCHANGED <<inputEpoch, inputOpen, capacity, committedQuality,
                    presentedQuality, candidateQuality, candidateOrigin,
                    overscanRejected, protectedFloorRejected, staticStarts>>
    /\ UNCHANGED terminalCeiling

\* A rejected richer renderer cut leaves a temporary global safety ceiling.
\* One certified occurrence-local allocation consumes that guard and makes
\* the capacity result terminal.  No unchanged frame may recreate this action.
CompleteRejectedReconciliation ==
    /\ phase = "reconcile"
    /\ nextCutRejected
    /\ globalCeiling
    /\ phase' = "constrained"
    /\ globalCeiling' = FALSE
    /\ UNCHANGED <<inputEpoch, inputOpen, capacity, committedQuality,
                    presentedQuality, candidateQuality, candidateOrigin,
                    overscanRejected, protectedFloorRejected, nextCutRejected, staticStarts>>
    /\ UNCHANGED terminalCeiling

\* Some renderer-wide cuts cannot be represented by a different occurrence-
\* local allocation.  The last completed framebuffer is then the exact hard-
\* deadline proof.  Retain its ceiling as terminal policy rather than leaving
\* the reconciliation phase without an effect owner.
RetainPresentedMinimum ==
    /\ phase = "reconcile"
    /\ nextCutRejected
    /\ globalCeiling
    /\ phase' = "constrained"
    /\ globalCeiling' = FALSE
    /\ terminalCeiling' = TRUE
    /\ UNCHANGED <<inputEpoch, inputOpen, capacity, committedQuality,
                    presentedQuality, candidateQuality, candidateOrigin,
                    overscanRejected, protectedFloorRejected, nextCutRejected, staticStarts>>

\* Only a completed frame commits a richer cut.  Until this edge, the previous
\* complete framebuffer remains the presentation authority.
CompleteFrame ==
    /\ phase = "present"
    /\ candidateOrigin \in {"overscan", "floor", "marginal"}
    /\ candidateQuality > committedQuality
    /\ committedQuality' = candidateQuality
    /\ presentedQuality' = candidateQuality
    /\ candidateQuality' = 0
    /\ candidateOrigin' = "none"
    /\ phase' = IF candidateQuality = MaxQuality
                 THEN "ready"
                 ELSE "marginal"
    /\ UNCHANGED <<inputEpoch, inputOpen, capacity,
                    overscanRejected, protectedFloorRejected, nextCutRejected, staticStarts>>
    /\ UNCHANGED <<globalCeiling, terminalCeiling>>

\* A completed fractional population may be the richest presentation proved
\* affordable by the hard static deadline even though the next complete cut
\* does not fit.  This is a successful terminal certificate, not a rejected
\* trial.  Keeping it distinct prevents ordinary convergence from disabling
\* its quality floor or reopening the same candidate.
AcceptCompletedBoundedFrame ==
    /\ phase = "present"
    /\ candidateOrigin \in {"floor", "marginal"}
    /\ candidateQuality > committedQuality
    /\ candidateQuality < MaxQuality
    /\ committedQuality' = candidateQuality
    /\ presentedQuality' = candidateQuality
    /\ candidateQuality' = 0
    /\ candidateOrigin' = "none"
    /\ phase' = "accepted"
    /\ UNCHANGED <<inputEpoch, inputOpen, capacity,
                    overscanRejected, protectedFloorRejected, nextCutRejected, staticStarts>>
    /\ UNCHANGED <<globalCeiling, terminalCeiling>>

Next ==
    \/ BeginInput
    \/ EndInput
    \/ StartStatic
    \/ FinishAlreadyRich
    \/ AdmitStaticOverscan
    \/ RejectStaticOverscan
    \/ AdmitProtectedFloor
    \/ RejectProtectedFloor
    \/ AdmitMarginal
    \/ RejectMarginal
    \/ CompleteRejectedReconciliation
    \/ RetainPresentedMinimum
    \/ CompleteFrame
    \/ AcceptCompletedBoundedFrame

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(BeginInput)
    /\ WF_vars(EndInput)
    /\ WF_vars(StartStatic)
    /\ WF_vars(FinishAlreadyRich)
    /\ WF_vars(AdmitStaticOverscan)
    /\ WF_vars(RejectStaticOverscan)
    /\ WF_vars(AdmitProtectedFloor)
    /\ WF_vars(RejectProtectedFloor)
    /\ WF_vars(AdmitMarginal)
    /\ WF_vars(RejectMarginal)
    /\ WF_vars(CompleteRejectedReconciliation)
    /\ WF_vars(RetainPresentedMinimum)
    /\ WF_vars(CompleteFrame)
    /\ WF_vars(AcceptCompletedBoundedFrame)

\* Renderer state is committed only by a completed frame.  A failed trial or
\* input interruption therefore cannot make the visible image regress.
CompleteFrameIsCommitAuthority == committedQuality = presentedQuality

CandidateDoesNotReplacePresentation ==
    phase = "present" => presentedQuality < candidateQuality

OverscanRejectionPreservesAffordableFloor ==
    overscanRejected /\ capacity >= FloorTarget =>
        ~protectedFloorRejected

ProtectedFloorRejectionHasExactEvidence ==
    protectedFloorRejected => capacity < FloorTarget

ProtectedFloorFallbackHasOwner ==
    protectedFloorRejected /\ ~nextCutRejected /\ ~Terminal =>
        phase \in {"marginal", "present"}

MarginalCapacityRemainsEnabled ==
    /\ protectedFloorRejected
    /\ phase = "marginal"
    /\ committedQuality < MaxQuality
    /\ capacity >= committedQuality + 1
    => ENABLED AdmitMarginal

ConstrainedHasCapacityWitness ==
    phase = "constrained" =>
        /\ nextCutRejected
        /\ committedQuality < MaxQuality
        /\ capacity < committedQuality + 1
        /\ ~globalCeiling

RejectedCeilingHasOneOwner ==
    globalCeiling => phase = "reconcile" /\ nextCutRejected

TerminalCeilingIsPolicy ==
    terminalCeiling => phase = "constrained" /\ nextCutRejected

ReadyIsRichest == phase = "ready" => committedQuality = MaxQuality

AcceptedHasCompletedQuality ==
    phase = "accepted" =>
        /\ committedQuality = presentedQuality
        /\ committedQuality > 0
        /\ committedQuality < MaxQuality
        /\ ~nextCutRejected

\* A quiet epoch may start the static protocol once.  Only BeginInput resets
\* this counter, preventing reject/reconcile/rearm cycles.
AtMostOneStaticStartPerEpoch == staticStarts <= 1

TerminalIsQuiescent ==
    Terminal =>
        /\ candidateOrigin = "none"
        /\ candidateQuality = 0
        /\ ~globalCeiling

EventuallyTerminalAfterFinalInput ==
    []((inputEpoch = MaxInputEpoch /\ ~inputOpen) => <>Terminal)

RejectedFloorEventuallyTerminates ==
    []((inputEpoch = MaxInputEpoch /\ ~inputOpen /\
        protectedFloorRejected) => <>Terminal)

RejectedCeilingEventuallyConsumed ==
    []((inputEpoch = MaxInputEpoch /\ ~inputOpen /\ globalCeiling) =>
        <>~globalCeiling)

=============================================================================
