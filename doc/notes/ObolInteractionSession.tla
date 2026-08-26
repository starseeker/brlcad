----------------------- MODULE ObolInteractionSession -----------------------
\* Camera-interaction lifecycle below the canonical progressive pipeline.
\*
\* Time is represented only by explicit deadline/debounce events.  A Qt timer
\* may deliver either event, but firing a timer is not itself control evidence.
\* The model proves that the three legal phases cannot encode the former
\* gesture-only combination and that a completed or deadline-retired motion
\* frame cannot strand a closed input session outside quiet state.

EXTENDS Naturals, TLC

CONSTANT MaxRenderSerial

Phases == {"quiet", "debouncing", "gesture"}

VARIABLES phase,
          inputClosed,
          renderSerial,
          requiredRenderSerial,
          debounceExpired

vars == <<phase, inputClosed, renderSerial, requiredRenderSerial,
          debounceExpired>>

TypeOK ==
    /\ phase \in Phases
    /\ inputClosed \in BOOLEAN
    /\ renderSerial \in 0..MaxRenderSerial
    /\ requiredRenderSerial \in 0..MaxRenderSerial
    /\ debounceExpired \in BOOLEAN

Init ==
    /\ phase = "quiet"
    /\ inputClosed = FALSE
    /\ renderSerial = 0
    /\ requiredRenderSerial = 0
    /\ debounceExpired = FALSE

BeginGesture ==
    /\ ~inputClosed
    /\ phase # "gesture"
    /\ phase' = "gesture"
    /\ requiredRenderSerial' = 0
    /\ debounceExpired' = FALSE
    /\ UNCHANGED <<inputClosed, renderSerial>>

CameraChanged ==
    /\ ~inputClosed
    /\ renderSerial < MaxRenderSerial
    /\ phase' = IF phase = "gesture" THEN "gesture" ELSE "debouncing"
    /\ requiredRenderSerial' = renderSerial + 1
    /\ debounceExpired' = FALSE
    /\ UNCHANGED <<inputClosed, renderSerial>>

EndGesture ==
    /\ phase = "gesture"
    /\ phase' = "debouncing"
    /\ requiredRenderSerial' = 0
    /\ debounceExpired' = FALSE
    /\ UNCHANGED <<inputClosed, renderSerial>>

CompleteMotionFrame ==
    /\ phase = "debouncing"
    /\ requiredRenderSerial > renderSerial
    /\ renderSerial' = requiredRenderSerial
    /\ requiredRenderSerial' = 0
    /\ debounceExpired' = FALSE
    /\ UNCHANGED <<phase, inputClosed>>

RetireExpiredMotionFrame ==
    /\ phase = "debouncing"
    /\ requiredRenderSerial > 0
    /\ requiredRenderSerial' = 0
    /\ debounceExpired' = TRUE
    /\ UNCHANGED <<phase, inputClosed, renderSerial>>

ExpireDebounce ==
    /\ phase = "debouncing"
    /\ requiredRenderSerial = 0
    /\ ~debounceExpired
    /\ debounceExpired' = TRUE
    /\ UNCHANGED <<phase, inputClosed, renderSerial,
                    requiredRenderSerial>>

FinishQuiet ==
    /\ phase = "debouncing"
    /\ requiredRenderSerial = 0
    /\ debounceExpired
    /\ phase' = "quiet"
    /\ debounceExpired' = FALSE
    /\ UNCHANGED <<inputClosed, renderSerial, requiredRenderSerial>>

CloseInput ==
    /\ ~inputClosed
    /\ inputClosed' = TRUE
    /\ UNCHANGED <<phase, renderSerial, requiredRenderSerial,
                    debounceExpired>>

Next ==
    \/ BeginGesture
    \/ CameraChanged
    \/ EndGesture
    \/ CompleteMotionFrame
    \/ RetireExpiredMotionFrame
    \/ ExpireDebounce
    \/ FinishQuiet
    \/ CloseInput

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(EndGesture)
    /\ WF_vars(CompleteMotionFrame)
    /\ WF_vars(RetireExpiredMotionFrame)
    /\ WF_vars(ExpireDebounce)
    /\ WF_vars(FinishQuiet)

QuietHasNoMotionGate ==
    phase = "quiet" =>
        requiredRenderSerial = 0 /\ ~debounceExpired

GestureCannotBeQuiet ==
    phase = "gesture" => phase # "quiet"

MotionGateNamesAFutureFrame ==
    requiredRenderSerial > 0 =>
        requiredRenderSerial > renderSerial

DebouncedSessionCanFinish ==
    phase = "debouncing" /\ requiredRenderSerial = 0 /\
    debounceExpired => ENABLED FinishQuiet

ActiveSessionHasProgressWitness ==
    phase # "quiet" =>
        \/ ENABLED EndGesture
        \/ ENABLED CompleteMotionFrame
        \/ ENABLED RetireExpiredMotionFrame
        \/ ENABLED ExpireDebounce
        \/ ENABLED FinishQuiet

EventuallyQuietAfterFinalInput ==
    [](inputClosed => <> (phase = "quiet"))

=============================================================================
