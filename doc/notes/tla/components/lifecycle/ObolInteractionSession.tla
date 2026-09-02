----------------------- MODULE ObolInteractionSession -----------------------
\* Camera-interaction lifecycle below the canonical progressive pipeline.
\*
\* Time is represented only by explicit deadline/debounce events.  A Qt timer
\* may deliver either event, but firing a timer is not itself control evidence.
\* The model proves that the three legal phases cannot encode the former
\* gesture-only combination and that a completed or deadline-retired motion
\* frame cannot strand a closed input session outside quiet state.  A
\* deadline-interrupted replay is explicit: expiry retires it so replay
\* priority cannot preempt the enabled quiet transition.  The model also
\* covers the mode edge which caused a real ownerless loop: if deterministic
\* unbounded terminal convergence is entered before a late quiet transition,
\* that transition cannot restore a finite presentation handoff.

EXTENDS Naturals, TLC

CONSTANT MaxRenderSerial

Phases == {"quiet", "debouncing", "gesture"}

VARIABLES phase,
          inputClosed,
          renderSerial,
          requiredRenderSerial,
          debounceExpired,
          interruptedReplay,
          terminalMode,
          finiteHandoff

vars == <<phase, inputClosed, renderSerial, requiredRenderSerial,
          debounceExpired, interruptedReplay, terminalMode, finiteHandoff>>

TypeOK ==
    /\ phase \in Phases
    /\ inputClosed \in BOOLEAN
    /\ renderSerial \in 0..MaxRenderSerial
    /\ requiredRenderSerial \in 0..MaxRenderSerial
    /\ debounceExpired \in BOOLEAN
    /\ interruptedReplay \in BOOLEAN
    /\ terminalMode \in BOOLEAN
    /\ finiteHandoff \in BOOLEAN

Init ==
    /\ phase = "quiet"
    /\ inputClosed = FALSE
    /\ renderSerial = 0
    /\ requiredRenderSerial = 0
    /\ debounceExpired = FALSE
    /\ interruptedReplay = FALSE
    /\ terminalMode = FALSE
    /\ finiteHandoff = FALSE

BeginGesture ==
    /\ ~inputClosed
    /\ phase # "gesture"
    /\ phase' = "gesture"
    /\ requiredRenderSerial' = 0
    /\ debounceExpired' = FALSE
    /\ interruptedReplay' = FALSE
    /\ UNCHANGED <<inputClosed, renderSerial, terminalMode,
                    finiteHandoff>>

CameraChanged ==
    /\ ~inputClosed
    /\ renderSerial < MaxRenderSerial
    /\ phase' = IF phase = "gesture" THEN "gesture" ELSE "debouncing"
    /\ requiredRenderSerial' = renderSerial + 1
    /\ debounceExpired' = FALSE
    /\ interruptedReplay' = FALSE
    /\ UNCHANGED <<inputClosed, renderSerial, terminalMode,
                    finiteHandoff>>

EndGesture ==
    /\ phase = "gesture"
    /\ phase' = "debouncing"
    /\ requiredRenderSerial' = 0
    /\ debounceExpired' = FALSE
    /\ UNCHANGED <<inputClosed, renderSerial, interruptedReplay, terminalMode,
                    finiteHandoff>>

CompleteMotionFrame ==
    /\ phase = "debouncing"
    /\ requiredRenderSerial > renderSerial
    /\ renderSerial' = requiredRenderSerial
    /\ requiredRenderSerial' = 0
    /\ debounceExpired' = FALSE
    /\ interruptedReplay' = FALSE
    /\ UNCHANGED <<phase, inputClosed, terminalMode, finiteHandoff>>

RetireExpiredMotionFrame ==
    /\ phase = "debouncing"
    /\ requiredRenderSerial > 0
    /\ requiredRenderSerial' = 0
    /\ debounceExpired' = TRUE
    /\ interruptedReplay' = FALSE
    /\ UNCHANGED <<phase, inputClosed, renderSerial, terminalMode,
                    finiteHandoff>>

ExpireDebounce ==
    /\ phase = "debouncing"
    /\ requiredRenderSerial = 0
    /\ ~debounceExpired
    /\ debounceExpired' = TRUE
    /\ interruptedReplay' = FALSE
    /\ UNCHANGED <<phase, inputClosed, renderSerial,
                    requiredRenderSerial, terminalMode, finiteHandoff>>

DeadlineInterrupted ==
    /\ phase # "quiet"
    /\ ~interruptedReplay
    /\ ~(phase = "debouncing" /\ requiredRenderSerial = 0 /\
         debounceExpired)
    /\ interruptedReplay' = TRUE
    /\ UNCHANGED <<phase, inputClosed, renderSerial,
                    requiredRenderSerial, debounceExpired,
                    terminalMode, finiteHandoff>>

CompleteInterruptedReplay ==
    /\ interruptedReplay
    /\ interruptedReplay' = FALSE
    /\ UNCHANGED <<phase, inputClosed, renderSerial,
                    requiredRenderSerial, debounceExpired,
                    terminalMode, finiteHandoff>>

FinishQuiet ==
    /\ phase = "debouncing"
    /\ requiredRenderSerial = 0
    /\ debounceExpired
    /\ ~interruptedReplay
    /\ phase' = "quiet"
    /\ debounceExpired' = FALSE
    /\ finiteHandoff' = IF terminalMode THEN FALSE ELSE TRUE
    /\ UNCHANGED <<inputClosed, renderSerial, requiredRenderSerial,
                    interruptedReplay, terminalMode>>

EnterTerminal ==
    /\ ~terminalMode
    /\ terminalMode' = TRUE
    /\ finiteHandoff' = FALSE
    /\ interruptedReplay' = FALSE
    /\ UNCHANGED <<phase, inputClosed, renderSerial,
                    requiredRenderSerial, debounceExpired>>

LeaveTerminal ==
    /\ terminalMode
    /\ terminalMode' = FALSE
    /\ UNCHANGED <<phase, inputClosed, renderSerial,
                    requiredRenderSerial, debounceExpired,
                    interruptedReplay, finiteHandoff>>

CompleteFiniteHandoff ==
    /\ finiteHandoff
    /\ finiteHandoff' = FALSE
    /\ UNCHANGED <<phase, inputClosed, renderSerial,
                    requiredRenderSerial, debounceExpired,
                    interruptedReplay, terminalMode>>

CloseInput ==
    /\ ~inputClosed
    /\ inputClosed' = TRUE
    \* Source/session teardown may arrive without a final gesture-up event.
    \* Retire that environment-owned edge atomically so the remaining timer,
    \* frame, and handoff owners are all schedulable production actions.
    /\ phase' = IF phase = "gesture" THEN "debouncing" ELSE phase
    /\ requiredRenderSerial' =
         IF phase = "gesture" THEN 0 ELSE requiredRenderSerial
    /\ debounceExpired' =
         IF phase = "gesture" THEN FALSE ELSE debounceExpired
    /\ interruptedReplay' =
         IF phase = "gesture" THEN FALSE ELSE interruptedReplay
    /\ UNCHANGED <<renderSerial, terminalMode, finiteHandoff>>

Next ==
    \/ BeginGesture
    \/ CameraChanged
    \/ EndGesture
    \/ CompleteMotionFrame
    \/ RetireExpiredMotionFrame
    \/ ExpireDebounce
    \/ DeadlineInterrupted
    \/ CompleteInterruptedReplay
    \/ FinishQuiet
    \/ EnterTerminal
    \/ LeaveTerminal
    \/ CompleteFiniteHandoff
    \/ CloseInput

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(CompleteMotionFrame)
    /\ WF_vars(RetireExpiredMotionFrame)
    /\ WF_vars(ExpireDebounce)
    /\ WF_vars(CompleteInterruptedReplay)
    /\ WF_vars(FinishQuiet)
    /\ WF_vars(CompleteFiniteHandoff)

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
    debounceExpired => ~interruptedReplay /\ ENABLED FinishQuiet

QuietReadyHasNoInterruptedReplay ==
    phase = "debouncing" /\ requiredRenderSerial = 0 /\
    debounceExpired => ~interruptedReplay

ActiveSessionHasInternalOrExternalWitness ==
    phase # "quiet" =>
        \/ /\ phase = "gesture"
           /\ ~inputClosed
        \/ ENABLED CompleteMotionFrame
        \/ ENABLED RetireExpiredMotionFrame
        \/ ENABLED ExpireDebounce
        \/ ENABLED CompleteInterruptedReplay
        \/ ENABLED FinishQuiet

TerminalModeHasNoFiniteHandoff ==
    terminalMode => ~finiteHandoff

FiniteHandoffHasProgressWitness ==
    finiteHandoff => ENABLED CompleteFiniteHandoff

Ready == phase = "quiet" /\ ~finiteHandoff
Constrained == FALSE
Failed == FALSE
Terminal == Ready \/ Constrained \/ Failed

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Terminal

EventuallyQuietAfterFinalInput ==
    [](inputClosed => <> (phase = "quiet"))

EventuallyDrainedAfterFinalInput ==
    [](inputClosed => <> (phase = "quiet" /\ ~finiteHandoff))

=============================================================================
