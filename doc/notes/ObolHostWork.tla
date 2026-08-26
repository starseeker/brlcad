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
          renderPending,
          renderCapacity,
          renderRevision,
          notificationPending,
          hostLoopScheduled,
          renderInFlight,
          capturedRenderRevision,
          capturedRenderCapacity,
          acknowledgedRenderRevision,
          producersClosed

vars == <<pumpPending, renderPending, renderCapacity, renderRevision,
          notificationPending, hostLoopScheduled, renderInFlight,
          capturedRenderRevision, capturedRenderCapacity,
          acknowledgedRenderRevision, producersClosed>>

WorkPending == pumpPending \/ renderPending

TypeOK == /\ pumpPending \in BOOLEAN
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
          /\ producersClosed \in BOOLEAN

Init == /\ pumpPending = FALSE
        /\ renderPending = FALSE
        /\ renderCapacity = FALSE
        /\ renderRevision = 0
        /\ notificationPending = FALSE
        /\ hostLoopScheduled = FALSE
        /\ renderInFlight = FALSE
        /\ capturedRenderRevision = 0
        /\ capturedRenderCapacity = FALSE
        /\ acknowledgedRenderRevision = 0
        /\ producersClosed = FALSE

\* Publishing either work level also supplies the edge which starts a host
\* loop.  Repeated publications coalesce into the standing level.
PublishPump == /\ ~producersClosed
               /\ ~pumpPending
               /\ pumpPending' = TRUE
               /\ notificationPending' = TRUE
               /\ UNCHANGED <<renderPending, renderCapacity, renderRevision,
                               hostLoopScheduled, renderInFlight,
                               capturedRenderRevision,
                               capturedRenderCapacity,
                               acknowledgedRenderRevision, producersClosed>>

PublishPresentation ==
    /\ ~producersClosed
    /\ ~renderPending
    /\ renderRevision < MaxRevision
    /\ renderPending' = TRUE
    /\ renderCapacity' = FALSE
    /\ renderRevision' = renderRevision + 1
    /\ notificationPending' = TRUE
    /\ UNCHANGED <<pumpPending, hostLoopScheduled, renderInFlight,
                    capturedRenderRevision, capturedRenderCapacity,
                    acknowledgedRenderRevision, producersClosed>>

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
    /\ notificationPending' =
        IF renderPending THEN notificationPending ELSE TRUE
    /\ UNCHANGED <<pumpPending, hostLoopScheduled, renderInFlight,
                    capturedRenderRevision, capturedRenderCapacity,
                    acknowledgedRenderRevision, producersClosed>>

CloseProducers == /\ ~producersClosed
                  /\ producersClosed' = TRUE
                  /\ UNCHANGED <<pumpPending, renderPending, renderCapacity,
                                  renderRevision, notificationPending,
                                  hostLoopScheduled, renderInFlight,
                                  capturedRenderRevision,
                                  capturedRenderCapacity,
                                  acknowledgedRenderRevision>>

\* The toolkit callback never performs work recursively.  It coalesces onto
\* one owner-thread loop and releases the callback edge.
DispatchNotification == /\ notificationPending
                        /\ notificationPending' = FALSE
                        /\ hostLoopScheduled' = TRUE
                        /\ UNCHANGED <<pumpPending, renderPending,
                                      renderCapacity, renderRevision,
                                      renderInFlight,
                                      capturedRenderRevision,
                                      capturedRenderCapacity,
                                      acknowledgedRenderRevision,
                                      producersClosed>>

\* A bounded controller pump discharges its current obligation.  If a render
\* is already standing, the same host loop proceeds to it without requiring a
\* second edge.
Pump == /\ hostLoopScheduled
        /\ pumpPending
        /\ pumpPending' = FALSE
        /\ hostLoopScheduled' = renderPending
        /\ UNCHANGED <<renderPending, renderCapacity, renderRevision,
                       notificationPending, renderInFlight,
                       capturedRenderRevision, capturedRenderCapacity,
                       acknowledgedRenderRevision, producersClosed>>

\* Pumping precedes rendering.  Claiming the transaction clears its pending
\* level before traversal; a later publication therefore owns a new revision.
BeginRender == /\ hostLoopScheduled
               /\ ~pumpPending
               /\ renderPending
               /\ ~renderInFlight
               /\ renderInFlight' = TRUE
               /\ capturedRenderRevision' = renderRevision
               /\ capturedRenderCapacity' = renderCapacity
               /\ renderPending' = FALSE
               /\ renderCapacity' = FALSE
               /\ hostLoopScheduled' = FALSE
               /\ UNCHANGED <<pumpPending, renderRevision,
                               notificationPending,
                               acknowledgedRenderRevision, producersClosed>>

CompleteRender == /\ renderInFlight
                  /\ renderInFlight' = FALSE
                  /\ acknowledgedRenderRevision' = capturedRenderRevision
                  /\ hostLoopScheduled' = (pumpPending \/ renderPending)
                  /\ UNCHANGED <<pumpPending, renderPending, renderCapacity,
                                  renderRevision, notificationPending,
                                  capturedRenderRevision,
                                  capturedRenderCapacity, producersClosed>>

RetireEmptyLoop == /\ hostLoopScheduled
                   /\ ~WorkPending
                   /\ ~renderInFlight
                   /\ hostLoopScheduled' = FALSE
                   /\ UNCHANGED <<pumpPending, renderPending, renderCapacity,
                                   renderRevision, notificationPending,
                                   renderInFlight, capturedRenderRevision,
                                   capturedRenderCapacity,
                                   acknowledgedRenderRevision,
                                   producersClosed>>

Next == \/ PublishPump
        \/ PublishPresentation
        \/ PublishCapacity
        \/ CloseProducers
        \/ DispatchNotification
        \/ Pump
        \/ BeginRender
        \/ CompleteRender
        \/ RetireEmptyLoop

Spec == Init /\ [][Next]_vars
        /\ WF_vars(DispatchNotification)
        /\ WF_vars(Pump)
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

\* Once producers close, fair host service drains all levels and witnesses.
EventuallyQuiescent ==
    [](producersClosed =>
       <>(~WorkPending /\ ~notificationPending /\
          ~hostLoopScheduled /\ ~renderInFlight))

=============================================================================
