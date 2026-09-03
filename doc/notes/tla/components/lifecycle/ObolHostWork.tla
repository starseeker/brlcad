----------------------------- MODULE ObolHostWork -----------------------------
\* Host/controller work-level protocol for an event-driven Obol endpoint.
\*
\* This deliberately models only the concurrency seam.  Mesh selection,
\* residency, and rendering algorithms are data-plane concerns.  The control
\* contract here is that controller work is level-triggered, every standing
\* level has a progress witness, duplicate pending requests are idempotent, and
\* claiming a request before traversal makes work published during that
\* traversal a distinct successor transaction.

EXTENDS Naturals, TLC

CONSTANT MaxRevision

VARIABLES pumpPending,
          controllerPumpPending,
          renderPending,
          renderCapacity,
          renderRevision,
          notificationPending,
          hostLoopScheduled,
          renderInFlight,
          capturedRenderRevision,
          capturedRenderCapacity,
          acknowledgedRenderRevision,
          providerRegistered,
          providerPending,
          sourceInputsPending,
          exactFrameRequired,
          exactPumpPending,
          capturedExactFrameRequired,
          publicationState,
          capturedPublicationFrame,
          producersClosed

vars == <<pumpPending, controllerPumpPending, renderPending, renderCapacity, renderRevision,
          notificationPending, hostLoopScheduled, renderInFlight,
          capturedRenderRevision, capturedRenderCapacity,
          acknowledgedRenderRevision, providerRegistered, providerPending,
          sourceInputsPending, exactFrameRequired, exactPumpPending,
          capturedExactFrameRequired, publicationState,
          capturedPublicationFrame, producersClosed>>

WorkPending == pumpPending \/ renderPending

TypeOK == /\ pumpPending \in BOOLEAN
          /\ controllerPumpPending \in BOOLEAN
          /\ renderPending \in BOOLEAN
          /\ renderCapacity \in BOOLEAN
          /\ (~renderPending => ~renderCapacity)
          /\ renderRevision \in 0..MaxRevision
          /\ notificationPending \in BOOLEAN
          /\ hostLoopScheduled \in BOOLEAN
          /\ renderInFlight \in BOOLEAN
          /\ capturedRenderRevision \in 0..MaxRevision
          /\ capturedRenderCapacity \in BOOLEAN
          /\ acknowledgedRenderRevision \in 0..MaxRevision
          /\ acknowledgedRenderRevision <= renderRevision
          /\ providerRegistered \in BOOLEAN
          /\ providerPending \in BOOLEAN
          /\ sourceInputsPending \in BOOLEAN
          /\ exactFrameRequired \in BOOLEAN
          /\ exactPumpPending \in BOOLEAN
          /\ capturedExactFrameRequired \in BOOLEAN
          /\ publicationState \in {"idle", "timer", "frame"}
          /\ capturedPublicationFrame \in BOOLEAN
          /\ producersClosed \in BOOLEAN

Init == /\ pumpPending = FALSE
        /\ controllerPumpPending = FALSE
        /\ renderPending = FALSE
        /\ renderCapacity = FALSE
        /\ renderRevision = 0
        /\ notificationPending = FALSE
        /\ hostLoopScheduled = FALSE
        /\ renderInFlight = FALSE
        /\ capturedRenderRevision = 0
        /\ capturedRenderCapacity = FALSE
        /\ acknowledgedRenderRevision = 0
        /\ providerRegistered = FALSE
        /\ providerPending = FALSE
        /\ sourceInputsPending = FALSE
        /\ exactFrameRequired = FALSE
        /\ exactPumpPending = FALSE
        /\ capturedExactFrameRequired = FALSE
        /\ publicationState = "idle"
        /\ capturedPublicationFrame = FALSE
        /\ producersClosed = FALSE

\* Publishing either work level also supplies the edge which starts a host
\* loop.  Repeated publications coalesce into the standing level.
PublishPump == /\ ~producersClosed
               /\ ~controllerPumpPending
               /\ controllerPumpPending' = TRUE
               /\ pumpPending' = TRUE
               /\ notificationPending' =
                   IF pumpPending THEN notificationPending ELSE TRUE
               /\ UNCHANGED <<publicationState, capturedPublicationFrame>>
               /\ UNCHANGED <<renderPending, renderCapacity, renderRevision,
                               hostLoopScheduled, renderInFlight,
                               capturedRenderRevision,
                               capturedRenderCapacity,
                               acknowledgedRenderRevision,
                               providerRegistered, providerPending,
                               sourceInputsPending,
                               exactFrameRequired, exactPumpPending,
                               capturedExactFrameRequired, producersClosed>>

PublishPresentation ==
    /\ ~producersClosed
    /\ ~renderPending
    /\ renderRevision < MaxRevision
    /\ renderPending' = TRUE
    /\ renderCapacity' = FALSE
    /\ renderRevision' = renderRevision + 1
    /\ publicationState' =
        IF publicationState = "timer" THEN "frame" ELSE publicationState
    /\ pumpPending' =
        IF publicationState = "timer"
        THEN (controllerPumpPending \/ providerPending \/
              sourceInputsPending \/ exactPumpPending)
        ELSE pumpPending
    /\ notificationPending' = TRUE
    /\ UNCHANGED capturedPublicationFrame
    /\ UNCHANGED <<controllerPumpPending, hostLoopScheduled, renderInFlight,
                    capturedRenderRevision, capturedRenderCapacity,
                    acknowledgedRenderRevision, providerRegistered,
                    providerPending, sourceInputsPending,
                    exactFrameRequired, exactPumpPending,
                    capturedExactFrameRequired,
                    producersClosed>>

\* A capacity request either starts a transaction or upgrades one pending
\* presentation transaction exactly once.  Same-strength duplicate requests
\* are stuttering steps and do not manufacture evidence revisions.
PublishCapacity ==
    /\ ~producersClosed
    /\ (~renderPending \/ ~renderCapacity)
    /\ renderRevision < MaxRevision
    /\ renderPending' = TRUE
    /\ renderCapacity' = TRUE
    /\ renderRevision' = renderRevision + 1
    /\ publicationState' =
        IF publicationState = "timer" THEN "frame" ELSE publicationState
    /\ pumpPending' =
        IF publicationState = "timer"
        THEN (controllerPumpPending \/ providerPending \/
              sourceInputsPending \/ exactPumpPending)
        ELSE pumpPending
    /\ notificationPending' =
        IF renderPending THEN notificationPending ELSE TRUE
    /\ UNCHANGED capturedPublicationFrame
    /\ UNCHANGED <<controllerPumpPending, hostLoopScheduled, renderInFlight,
                    capturedRenderRevision, capturedRenderCapacity,
                    acknowledgedRenderRevision, providerRegistered,
                    providerPending, sourceInputsPending,
                    exactFrameRequired, exactPumpPending,
                    capturedExactFrameRequired,
                    producersClosed>>

\* Registration initially raises a conservative provider pump.  Whether the
\* callback reports more work is a data-plane decision outside this model.
RegisterProvider == /\ ~producersClosed
                    /\ ~providerRegistered
                    /\ providerRegistered' = TRUE
                    /\ providerPending' = TRUE
                    /\ pumpPending' = TRUE
                    /\ notificationPending' =
                        IF pumpPending THEN notificationPending ELSE TRUE
                    /\ UNCHANGED <<publicationState,
                                    capturedPublicationFrame>>
                    /\ UNCHANGED <<controllerPumpPending,
                                    renderPending, renderCapacity,
                                    renderRevision, hostLoopScheduled,
                                    renderInFlight, capturedRenderRevision,
                                    capturedRenderCapacity,
                                    acknowledgedRenderRevision,
                                    sourceInputsPending,
                                    exactFrameRequired, exactPumpPending,
                                    capturedExactFrameRequired,
                                    producersClosed>>

\* A registered provider may be idle indefinitely and later publish another
\* level.  Registration is capability; providerPending is the work fact.
PublishProviderWork ==
    /\ ~producersClosed
    /\ providerRegistered
    /\ ~providerPending
    /\ providerPending' = TRUE
    /\ pumpPending' = TRUE
    /\ notificationPending' =
        IF pumpPending THEN notificationPending ELSE TRUE
    /\ UNCHANGED <<publicationState, capturedPublicationFrame>>
    /\ UNCHANGED <<controllerPumpPending, renderPending, renderCapacity,
                    renderRevision, hostLoopScheduled, renderInFlight,
                    capturedRenderRevision, capturedRenderCapacity,
                    acknowledgedRenderRevision, providerRegistered,
                    sourceInputsPending,
                    exactFrameRequired, exactPumpPending,
                    capturedExactFrameRequired,
                    producersClosed>>

\* An interrupted or superseded CAD traversal leaves level-triggered exact
\* presentation debt.  The next bounded pump converts it into a render.
RequireExactFrame == /\ ~producersClosed
                     /\ ~renderInFlight
                     /\ ~exactFrameRequired
                     /\ renderRevision < MaxRevision
                     /\ exactFrameRequired' = TRUE
                     /\ exactPumpPending' = TRUE
                     /\ pumpPending' = TRUE
                     /\ notificationPending' =
                         IF pumpPending THEN notificationPending ELSE TRUE
                     /\ UNCHANGED <<publicationState,
                                     capturedPublicationFrame>>
                     /\ UNCHANGED <<controllerPumpPending,
                                     renderPending, renderCapacity,
                                     renderRevision, hostLoopScheduled,
                                     renderInFlight, capturedRenderRevision,
                                     capturedRenderCapacity,
                                     acknowledgedRenderRevision,
                                     providerRegistered, providerPending,
                                     sourceInputsPending,
                                     capturedExactFrameRequired,
                                     producersClosed>>

\* A compact hierarchy edit can change effective visibility without changing
\* immutable mesh inventory.  Presentation may already be correct, but the
\* planner still owes this exact source delta before readiness is current.
PublishSourceInputs ==
    /\ ~producersClosed
    /\ ~sourceInputsPending
    /\ sourceInputsPending' = TRUE
    /\ pumpPending' = TRUE
    /\ notificationPending' =
        IF pumpPending THEN notificationPending ELSE TRUE
    /\ UNCHANGED <<publicationState, capturedPublicationFrame>>
    /\ UNCHANGED <<controllerPumpPending, renderPending, renderCapacity,
                    renderRevision, hostLoopScheduled, renderInFlight,
                    capturedRenderRevision, capturedRenderCapacity,
                    acknowledgedRenderRevision, providerRegistered,
                    providerPending, exactFrameRequired, exactPumpPending,
                    capturedExactFrameRequired, producersClosed>>

\* Provider teardown owns only provider debt.  Synchronization may clear its
\* conservative pump only when no independent exact frame is outstanding.
UnregisterProvider == /\ providerRegistered
                      /\ providerRegistered' = FALSE
                      /\ providerPending' = FALSE
                      /\ pumpPending' =
                          (controllerPumpPending \/ sourceInputsPending \/
                           exactPumpPending \/
                           publicationState = "timer")
                      /\ UNCHANGED <<publicationState,
                                      capturedPublicationFrame>>
                      /\ UNCHANGED <<controllerPumpPending,
                                      renderPending, renderCapacity,
                                      renderRevision, notificationPending,
                                      hostLoopScheduled, renderInFlight,
                                      capturedRenderRevision,
                                      capturedRenderCapacity,
                                      acknowledgedRenderRevision,
                                      sourceInputsPending,
                                      exactFrameRequired, exactPumpPending,
                                      capturedExactFrameRequired,
                                      producersClosed>>

CloseProducers == /\ ~producersClosed
                  /\ producersClosed' = TRUE
                  /\ providerRegistered' = FALSE
                  /\ providerPending' = FALSE
                  /\ sourceInputsPending' = FALSE
                  /\ pumpPending' =
                      (controllerPumpPending \/ exactPumpPending \/
                       publicationState = "timer")
                  /\ UNCHANGED <<publicationState,
                                  capturedPublicationFrame>>
                  /\ UNCHANGED <<controllerPumpPending,
                                  renderPending, renderCapacity,
                                  renderRevision, notificationPending,
                                  hostLoopScheduled, renderInFlight,
                                  capturedRenderRevision,
                                  capturedRenderCapacity,
                                  acknowledgedRenderRevision,
                                  exactFrameRequired, exactPumpPending,
                                  capturedExactFrameRequired>>

\* Losing the endpoint is stronger than an orderly producer close.  There is
\* no host left to drain a notification, timer, pending request, or claimed
\* frame, so the close transaction cancels every such owner atomically.
CloseEndpoint ==
    /\ ~producersClosed
    /\ producersClosed' = TRUE
    /\ pumpPending' = FALSE
    /\ controllerPumpPending' = FALSE
    /\ renderPending' = FALSE
    /\ renderCapacity' = FALSE
    /\ notificationPending' = FALSE
    /\ hostLoopScheduled' = FALSE
    /\ renderInFlight' = FALSE
    /\ capturedRenderCapacity' = FALSE
    /\ providerRegistered' = FALSE
    /\ providerPending' = FALSE
    /\ sourceInputsPending' = FALSE
    /\ exactFrameRequired' = FALSE
    /\ exactPumpPending' = FALSE
    /\ capturedExactFrameRequired' = FALSE
    /\ publicationState' = "idle"
    /\ capturedPublicationFrame' = FALSE
    /\ UNCHANGED <<renderRevision, capturedRenderRevision,
                    acknowledgedRenderRevision>>

\* The toolkit callback never performs work recursively.  It coalesces onto
\* one owner-thread loop and releases the callback edge.
DispatchNotification == /\ notificationPending
                        /\ notificationPending' = FALSE
                        /\ hostLoopScheduled' = TRUE
                        /\ UNCHANGED <<publicationState,
                                        capturedPublicationFrame>>
                        /\ UNCHANGED <<pumpPending, controllerPumpPending,
                                      renderPending,
                                      renderCapacity, renderRevision,
                                      renderInFlight,
                                      capturedRenderRevision,
                                      capturedRenderCapacity,
                                      acknowledgedRenderRevision,
                                      providerRegistered, providerPending,
                                      sourceInputsPending,
                                      exactFrameRequired, exactPumpPending,
                                      capturedExactFrameRequired,
                                      producersClosed>>

\* Applied immutable results keep one finite batching timer.  Once it expires,
\* the same bounded host pump transfers that obligation to a coalesced render
\* and is no longer itself a runnable reason.
PublishPublicationBatch ==
    /\ ~producersClosed
    /\ ~renderInFlight
    /\ publicationState = "idle"
    /\ (renderRevision < MaxRevision \/
        (renderPending /\ renderCapacity))
    /\ publicationState' = "timer"
    /\ pumpPending' = TRUE
    /\ notificationPending' =
        IF pumpPending THEN notificationPending ELSE TRUE
    /\ UNCHANGED <<controllerPumpPending, renderPending, renderCapacity,
                    renderRevision, hostLoopScheduled, renderInFlight,
                    capturedRenderRevision, capturedRenderCapacity,
                    acknowledgedRenderRevision, providerRegistered,
                    providerPending, sourceInputsPending,
                    exactFrameRequired, exactPumpPending,
                    capturedExactFrameRequired, capturedPublicationFrame,
                    producersClosed>>

QueuePublicationFrame ==
    /\ hostLoopScheduled
    /\ pumpPending
    /\ publicationState = "timer"
    /\ (renderRevision < MaxRevision \/
        (renderPending /\ renderCapacity))
    /\ publicationState' = "frame"
    /\ renderPending' = TRUE
    /\ renderCapacity' = TRUE
    /\ renderRevision' =
        IF ~renderPending \/ ~renderCapacity
        THEN renderRevision + 1 ELSE renderRevision
    /\ notificationPending' =
        IF renderPending THEN notificationPending ELSE TRUE
    /\ pumpPending' =
        (controllerPumpPending \/ providerPending \/ sourceInputsPending \/
         exactPumpPending)
    /\ hostLoopScheduled' = TRUE
    /\ UNCHANGED <<renderInFlight, capturedRenderRevision,
                    capturedRenderCapacity, acknowledgedRenderRevision,
                    providerRegistered, providerPending,
                    sourceInputsPending,
                    exactFrameRequired, exactPumpPending,
                    capturedExactFrameRequired, controllerPumpPending,
                    capturedPublicationFrame, producersClosed>>

\* A bounded controller pump discharges its current obligation.  If a render
\* is already standing, the same host loop proceeds to it without requiring a
\* second edge.
Pump == /\ hostLoopScheduled
        /\ pumpPending
        /\ ~sourceInputsPending
        /\ ~exactPumpPending
        /\ publicationState # "timer"
        /\ pumpPending' = FALSE
        /\ controllerPumpPending' = FALSE
        /\ providerPending' = FALSE
        /\ hostLoopScheduled' = renderPending
        /\ UNCHANGED <<publicationState, capturedPublicationFrame>>
        /\ UNCHANGED <<renderPending, renderCapacity, renderRevision,
                       notificationPending, renderInFlight,
                       capturedRenderRevision, capturedRenderCapacity,
                       acknowledgedRenderRevision, providerRegistered,
                       sourceInputsPending, exactFrameRequired,
                       exactPumpPending,
                       capturedExactFrameRequired,
                       producersClosed>>

\* Source-input debt has its own reducer action.  A render cooldown may delay
\* this transition, but the generic pump is not allowed to clear the shared
\* wakeup level while the source revision remains unconsumed.
ConsumeSourceInputs ==
    /\ hostLoopScheduled
    /\ pumpPending
    /\ sourceInputsPending
    /\ sourceInputsPending' = FALSE
    /\ pumpPending' =
        (controllerPumpPending \/ providerPending \/ exactPumpPending \/
         publicationState = "timer")
    /\ hostLoopScheduled' = (pumpPending' \/ renderPending)
    /\ UNCHANGED <<publicationState, capturedPublicationFrame>>
    /\ UNCHANGED <<controllerPumpPending, renderPending, renderCapacity,
                    renderRevision, notificationPending, renderInFlight,
                    capturedRenderRevision, capturedRenderCapacity,
                    acknowledgedRenderRevision, providerRegistered,
                    providerPending, exactFrameRequired, exactPumpPending,
                    capturedExactFrameRequired, producersClosed>>

PumpExactIntoPendingRender ==
    /\ hostLoopScheduled
    /\ pumpPending
    /\ exactPumpPending
    /\ renderPending
    /\ pumpPending' =
        (sourceInputsPending \/ publicationState = "timer")
    /\ exactPumpPending' = FALSE
    /\ controllerPumpPending' = FALSE
    /\ providerPending' = FALSE
    /\ hostLoopScheduled' = TRUE
    /\ UNCHANGED <<publicationState, capturedPublicationFrame>>
    /\ UNCHANGED <<renderPending, renderCapacity, renderRevision,
                    notificationPending, renderInFlight,
                    capturedRenderRevision, capturedRenderCapacity,
                    acknowledgedRenderRevision, providerRegistered,
                    sourceInputsPending,
                    exactFrameRequired, capturedExactFrameRequired,
                    producersClosed>>

PumpExactRequestRender ==
    /\ hostLoopScheduled
    /\ pumpPending
    /\ exactPumpPending
    /\ ~renderPending
    /\ publicationState # "timer"
    /\ renderRevision < MaxRevision
    /\ pumpPending' =
        (sourceInputsPending \/ publicationState = "timer")
    /\ exactPumpPending' = FALSE
    /\ controllerPumpPending' = FALSE
    /\ providerPending' = FALSE
    /\ renderPending' = TRUE
    /\ renderCapacity' = FALSE
    /\ renderRevision' = renderRevision + 1
    /\ hostLoopScheduled' = TRUE
    /\ UNCHANGED <<publicationState, capturedPublicationFrame>>
    /\ UNCHANGED <<notificationPending, renderInFlight,
                    capturedRenderRevision, capturedRenderCapacity,
                    acknowledgedRenderRevision, providerRegistered,
                    sourceInputsPending,
                    exactFrameRequired, capturedExactFrameRequired,
                    producersClosed>>

\* Pumping precedes rendering.  Claiming the transaction clears its pending
\* level before traversal; a later publication therefore owns a new revision.
BeginRender == /\ hostLoopScheduled
               /\ ~pumpPending
               /\ renderPending
               /\ ~renderInFlight
               /\ renderInFlight' = TRUE
               /\ capturedRenderRevision' = renderRevision
               /\ capturedRenderCapacity' = renderCapacity
               /\ capturedExactFrameRequired' = exactFrameRequired
               /\ capturedPublicationFrame' =
                   (publicationState = "frame")
               /\ renderPending' = FALSE
               /\ renderCapacity' = FALSE
               /\ hostLoopScheduled' = FALSE
               /\ UNCHANGED <<pumpPending, controllerPumpPending,
                               renderRevision,
                               notificationPending,
                               acknowledgedRenderRevision,
                               providerRegistered, providerPending,
                               sourceInputsPending,
                               exactFrameRequired, exactPumpPending,
                               publicationState,
                               producersClosed>>

CompleteRender == /\ renderInFlight
                  /\ renderInFlight' = FALSE
                  /\ acknowledgedRenderRevision' = capturedRenderRevision
                  /\ exactFrameRequired' =
                      IF capturedExactFrameRequired THEN FALSE
                      ELSE exactFrameRequired
                  /\ capturedExactFrameRequired' = FALSE
                  /\ publicationState' =
                      IF capturedPublicationFrame THEN "idle"
                      ELSE publicationState
                  /\ capturedPublicationFrame' = FALSE
                  /\ hostLoopScheduled' = (pumpPending \/ renderPending)
                  /\ UNCHANGED <<pumpPending, controllerPumpPending,
                                  renderPending, renderCapacity,
                                  renderRevision, notificationPending,
                                  capturedRenderRevision,
                                  capturedRenderCapacity,
                                  providerRegistered, providerPending,
                                  sourceInputsPending,
                                  exactPumpPending,
                                  producersClosed>>

\* Consuming a completed frame may itself open controller work when the timing
\* sample changes a capacity bracket, presentation handoff, or demand scan.
\* This transition runs on the host thread, so it arms the next host loop
\* directly instead of relying on the false-to-true producer callback edge.
CompleteRenderAndOpenPump ==
    /\ renderInFlight
    /\ ~producersClosed
    /\ ~pumpPending
    /\ renderInFlight' = FALSE
    /\ acknowledgedRenderRevision' = capturedRenderRevision
    /\ exactFrameRequired' =
        IF capturedExactFrameRequired THEN FALSE ELSE exactFrameRequired
    /\ capturedExactFrameRequired' = FALSE
    /\ publicationState' =
        IF capturedPublicationFrame THEN "idle" ELSE publicationState
    /\ capturedPublicationFrame' = FALSE
    /\ pumpPending' = TRUE
    /\ controllerPumpPending' = TRUE
    /\ hostLoopScheduled' = TRUE
    /\ UNCHANGED <<renderPending, renderCapacity, renderRevision,
                    notificationPending, capturedRenderRevision,
                    capturedRenderCapacity, providerRegistered,
                    providerPending, sourceInputsPending,
                    exactPumpPending,
                    producersClosed>>

\* A completed-frame reducer may change the effective CAD presentation
\* controls after consuming the frame it just measured.  That side effect
\* creates a new exact successor and must publish its pump in the same owner-
\* thread transition; no provider callback is guaranteed to follow it.
CompleteRenderAndRequireExact ==
    /\ renderInFlight
    /\ ~producersClosed
    /\ renderRevision < MaxRevision
    /\ renderInFlight' = FALSE
    /\ acknowledgedRenderRevision' = capturedRenderRevision
    /\ exactFrameRequired' = TRUE
    /\ exactPumpPending' = TRUE
    /\ capturedExactFrameRequired' = FALSE
    /\ publicationState' =
        IF capturedPublicationFrame THEN "idle" ELSE publicationState
    /\ capturedPublicationFrame' = FALSE
    /\ pumpPending' = TRUE
    /\ hostLoopScheduled' = TRUE
    /\ UNCHANGED <<controllerPumpPending,
                    renderPending, renderCapacity, renderRevision,
                    notificationPending, capturedRenderRevision,
                    capturedRenderCapacity, providerRegistered,
                    providerPending, sourceInputsPending,
                    producersClosed>>

RetireEmptyLoop == /\ hostLoopScheduled
                   /\ ~WorkPending
                   /\ ~renderInFlight
                   /\ hostLoopScheduled' = FALSE
                   /\ UNCHANGED <<publicationState,
                                   capturedPublicationFrame>>
                   /\ UNCHANGED <<pumpPending, controllerPumpPending,
                                   renderPending, renderCapacity,
                                   renderRevision, notificationPending,
                                   renderInFlight, capturedRenderRevision,
                                   capturedRenderCapacity,
                                   acknowledgedRenderRevision,
                                   providerRegistered, providerPending,
                                   sourceInputsPending,
                                   exactFrameRequired, exactPumpPending,
                                   capturedExactFrameRequired,
                                   producersClosed>>

Next == \/ PublishPump
        \/ PublishPresentation
        \/ PublishCapacity
        \/ RegisterProvider
        \/ PublishProviderWork
        \/ PublishSourceInputs
        \/ PublishPublicationBatch
        \/ RequireExactFrame
        \/ UnregisterProvider
        \/ CloseProducers
        \/ CloseEndpoint
        \/ DispatchNotification
        \/ QueuePublicationFrame
        \/ Pump
        \/ ConsumeSourceInputs
        \/ PumpExactIntoPendingRender
        \/ PumpExactRequestRender
        \/ BeginRender
        \/ CompleteRender
        \/ CompleteRenderAndOpenPump
        \/ CompleteRenderAndRequireExact
        \/ RetireEmptyLoop

Spec == Init /\ [][Next]_vars
        /\ WF_vars(DispatchNotification)
        /\ WF_vars(QueuePublicationFrame)
        /\ WF_vars(Pump)
        /\ WF_vars(ConsumeSourceInputs)
        /\ WF_vars(PumpExactIntoPendingRender)
        /\ WF_vars(PumpExactRequestRender)
        /\ WF_vars(BeginRender)
        /\ WF_vars(CompleteRender)
        /\ WF_vars(RetireEmptyLoop)

\* A standing work level is never an orphaned Boolean.  It always has a
\* callback edge, an owner-thread loop, or a frame which is consuming it.
PendingWorkHasWitness ==
    (WorkPending \/ renderInFlight) =>
        (notificationPending \/ hostLoopScheduled \/ renderInFlight)

\* A frame never acknowledges a request which was published after it began.
NoFutureAcknowledgement ==
    acknowledgedRenderRevision <= renderRevision

\* Provider retirement may remove its own conservative pump but cannot orphan
\* the exact CAD frame owed by a mutation it just published.
ExactPresentationHasWitness ==
    exactFrameRequired =>
        (exactPumpPending \/ renderPending \/ renderInFlight)

PublicationHasWitness ==
    /\ publicationState = "timer" => pumpPending
    /\ publicationState = "frame" => (renderPending \/ renderInFlight)

\* The shared pump bit is an implementation wakeup level.  It must be exactly
\* the projection of current controller/provider work or exact-frame debt not
\* already represented by a pending or in-flight render.
PumpLevelMatchesReasons ==
    pumpPending =
        (controllerPumpPending \/ providerPending \/ sourceInputsPending \/
         exactPumpPending \/
         publicationState = "timer")

RegisteredProviderMayBeIdle ==
    providerRegistered /\ ~providerPending /\ ~controllerPumpPending /\
    ~sourceInputsPending /\ ~exactPumpPending /\ publicationState # "timer"
        => ~pumpPending

Terminal ==
    /\ producersClosed
    /\ ~WorkPending
    /\ ~notificationPending
    /\ ~hostLoopScheduled
    /\ ~renderInFlight
    /\ ~exactFrameRequired
    /\ ~exactPumpPending
    /\ publicationState = "idle"
    /\ ~capturedPublicationFrame
    /\ ~providerRegistered
    /\ ~providerPending
    /\ ~sourceInputsPending

\* Action-level postcondition exported to the lifecycle composition contract.
EndpointClosureRetiresAllWork == CloseEndpoint => Terminal'

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Terminal

\* Once producers close, fair host service drains all levels and witnesses.
EventuallyQuiescent ==
    [](producersClosed =>
       <>Terminal)

=============================================================================
