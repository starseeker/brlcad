----------------------- MODULE ObolDeferredAutoview ------------------------
\* Bounded ownership model for deferred progressive autoview.
\*
\* This models control ownership only.  It deliberately does not model Qt,
\* Coin, geometry bounds, or floating point camera matrices.  The production
\* contract is implemented by ged_draw_obol_progressive_autoview_* in
\* src/libged/draw_obol.cpp: a deferred fit owns the requested view's framing
\* (center/scale), never its orientation, and has exactly one terminal result.

EXTENDS Naturals, TLC

CONSTANT MaxRequests

Views == {"primary", "secondary"}
None == "none"

VARIABLES live,
          pendingOwner,
          boundsReady,
          producerPending,
          hostWakePending,
          presentationPending,
          requestSerial,
          terminalSerial,
          framing,
          orientation,
          appliedOwner

vars == <<live, pendingOwner, boundsReady, producerPending,
          hostWakePending, presentationPending, requestSerial, terminalSerial,
          framing, orientation, appliedOwner>>

TypeOK ==
    /\ live \subseteq Views
    /\ pendingOwner \in Views \cup {None}
    /\ boundsReady \in BOOLEAN
    /\ producerPending \in BOOLEAN
    /\ hostWakePending \in BOOLEAN
    /\ presentationPending \in BOOLEAN
    /\ requestSerial \in 0..MaxRequests
    /\ terminalSerial \in 0..requestSerial
    /\ framing \in [Views -> 0..1]
    /\ orientation \in [Views -> 0..1]
    /\ appliedOwner \in Views \cup {None}

Init ==
    /\ live = Views
    /\ pendingOwner = None
    /\ boundsReady = FALSE
    /\ producerPending = FALSE
    /\ hostWakePending = FALSE
    /\ presentationPending = FALSE
    /\ requestSerial = 0
    /\ terminalSerial = 0
    /\ framing = [view \in Views |-> 0]
    /\ orientation = [view \in Views |-> 0]
    /\ appliedOwner = None

\* A new request replaces any older pending request.  This is the only action
\* that creates pending ownership; provider ticks cannot rearm a completed
\* request on their own.
Request(view) ==
    /\ view \in live
    /\ requestSerial < MaxRequests
    /\ pendingOwner' = view
    /\ boundsReady' = FALSE
    /\ producerPending' = TRUE
    /\ hostWakePending' = FALSE
    /\ presentationPending' = FALSE
    /\ requestSerial' = requestSerial + 1
    /\ appliedOwner' = None
    /\ UNCHANGED <<live, terminalSerial, framing, orientation>>

BoundsBecomeReady ==
    /\ pendingOwner # None
    /\ producerPending
    /\ boundsReady' = TRUE
    /\ producerPending' = FALSE
    /\ hostWakePending' = TRUE
    /\ UNCHANGED <<live, pendingOwner, presentationPending, requestSerial,
                   terminalSerial, framing, orientation, appliedOwner>>

\* The producer-to-host callback may wake the owner thread without needing a
\* direct paint.  The host pump turns that wake into a presentation obligation;
\* it must not publish idle while a ready deferred camera operation remains.
HostPump ==
    /\ pendingOwner # None
    /\ hostWakePending
    /\ hostWakePending' = FALSE
    /\ presentationPending' = boundsReady
    /\ UNCHANGED <<live, pendingOwner, boundsReady, producerPending,
                   requestSerial, terminalSerial, framing, orientation,
                   appliedOwner>>

\* An A/E/orbit orientation change is compatible with a pending autoview.
\* Applying the fit must retain this later orientation.
OrientationInput(view) ==
    /\ view \in live
    /\ orientation' = [orientation EXCEPT ![view] = 1 - @]
    /\ UNCHANGED <<live, pendingOwner, boundsReady, producerPending,
                   hostWakePending, presentationPending, requestSerial,
                   terminalSerial, framing, appliedOwner>>

\* A pan, zoom, or explicit size/center owns framing and cancels only the
\* request for the same view.  Input in another view cannot cancel or redirect
\* the original request.
FramingInput(view) ==
    /\ view \in live
    /\ framing' = [framing EXCEPT ![view] = 1 - @]
    /\ pendingOwner' = IF pendingOwner = view THEN None ELSE pendingOwner
    /\ boundsReady' = IF pendingOwner = view THEN FALSE ELSE boundsReady
    /\ producerPending' = IF pendingOwner = view THEN FALSE ELSE producerPending
    /\ hostWakePending' = IF pendingOwner = view THEN FALSE ELSE hostWakePending
    /\ presentationPending' =
        IF pendingOwner = view THEN FALSE ELSE presentationPending
    /\ terminalSerial' = IF pendingOwner = view THEN requestSerial
                         ELSE terminalSerial
    /\ appliedOwner' = None
    /\ UNCHANGED <<live, requestSerial, orientation>>

\* A completed fit writes only the requested live view's framing.  It neither
\* changes orientation nor touches a different view.  Completion consumes the
\* request atomically, so an idle provider cannot apply it twice.
PresentAndApply ==
    /\ pendingOwner # None
    /\ pendingOwner \in live
    /\ boundsReady
    /\ presentationPending
    /\ framing' = [framing EXCEPT ![pendingOwner] = 1 - @]
    /\ terminalSerial' = requestSerial
    /\ appliedOwner' = pendingOwner
    /\ pendingOwner' = None
    /\ boundsReady' = FALSE
    /\ producerPending' = FALSE
    /\ hostWakePending' = FALSE
    /\ presentationPending' = FALSE
    /\ UNCHANGED <<live, requestSerial, orientation>>

\* View disposal is terminal for a request owned by that view.  A completion
\* must never write through a dead or replacement endpoint.
Teardown(view) ==
    /\ view \in live
    /\ live' = live \ {view}
    /\ pendingOwner' = IF pendingOwner = view THEN None ELSE pendingOwner
    /\ boundsReady' = IF pendingOwner = view THEN FALSE ELSE boundsReady
    /\ producerPending' = IF pendingOwner = view THEN FALSE ELSE producerPending
    /\ hostWakePending' = IF pendingOwner = view THEN FALSE ELSE hostWakePending
    /\ presentationPending' =
        IF pendingOwner = view THEN FALSE ELSE presentationPending
    /\ terminalSerial' = IF pendingOwner = view THEN requestSerial
                         ELSE terminalSerial
    /\ appliedOwner' = None
    /\ UNCHANGED <<requestSerial, framing, orientation>>

Next ==
    \/ \E view \in Views: Request(view)
    \/ BoundsBecomeReady
    \/ HostPump
    \/ \E view \in Views: OrientationInput(view)
    \/ \E view \in Views: FramingInput(view)
    \/ PresentAndApply
    \/ \E view \in Views: Teardown(view)

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(BoundsBecomeReady)
    /\ WF_vars(HostPump)
    /\ WF_vars(PresentAndApply)

\* A provider completion cannot be redirected to a different endpoint.
NoCrossViewApply ==
    appliedOwner # None =>
        /\ appliedOwner \in live
        /\ terminalSerial = requestSerial

\* No fit may be applied before a request has reached a terminal result.
NoPhantomCompletion == terminalSerial <= requestSerial

\* Once a request is quiescent, it has exactly one terminal outcome: apply,
\* cancellation by framing input, replacement, or teardown.  An orientation
\* change alone is not a cancellation path.
NoPendingAfterTerminal ==
    pendingOwner = None =>
        /\ ~boundsReady
        /\ ~producerPending
        /\ ~hostWakePending
        /\ ~presentationPending

\* A live deferred fit always has one liveness witness: the source producer,
\* a coalesced host wake, or a presentation obligation.  This is the contract
\* the fast-warm OSMesa regression violated when it dropped its final frame
\* request before the camera had been synchronized.
PendingFitHasLivenessWitness ==
    pendingOwner # None =>
        producerPending \/ hostWakePending \/ presentationPending

\* In a stable live view, fairness over ready data, host pumping, and the
\* presentation step guarantees the deferred request does not remain pending
\* indefinitely.
DeferredFitEventuallyTerminates ==
    [](pendingOwner # None /\ pendingOwner \in live =>
       <>(pendingOwner = None \/ terminalSerial = requestSerial))

=============================================================================
