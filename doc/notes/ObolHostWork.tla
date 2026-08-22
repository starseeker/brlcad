----------------------------- MODULE ObolHostWork -----------------------------
\* Host/controller work-level protocol for an event-driven Obol endpoint.
\*
\* This deliberately models only the concurrency seam.  Mesh selection,
\* residency, and rendering algorithms are data-plane concerns.  The control
\* contract here is that controller work is level-triggered, every standing
\* level has a progress witness, and a frame completion may retire only the
\* render revision it actually captured.

EXTENDS Naturals, TLC

CONSTANT MaxRevision

VARIABLES pumpPending,
          renderPending,
          renderRevision,
          notificationPending,
          hostLoopScheduled,
          renderInFlight,
          capturedRenderRevision,
          acknowledgedRenderRevision,
          producersClosed

vars == <<pumpPending, renderPending, renderRevision,
          notificationPending, hostLoopScheduled, renderInFlight,
          capturedRenderRevision, acknowledgedRenderRevision,
          producersClosed>>

WorkPending == pumpPending \/ renderPending

TypeOK == /\ pumpPending \in BOOLEAN
          /\ renderPending \in BOOLEAN
          /\ renderRevision \in 0..MaxRevision
          /\ notificationPending \in BOOLEAN
          /\ hostLoopScheduled \in BOOLEAN
          /\ renderInFlight \in BOOLEAN
          /\ capturedRenderRevision \in 0..MaxRevision
          /\ acknowledgedRenderRevision \in 0..MaxRevision
          /\ acknowledgedRenderRevision <= renderRevision
          /\ producersClosed \in BOOLEAN

Init == /\ pumpPending = FALSE
        /\ renderPending = FALSE
        /\ renderRevision = 0
        /\ notificationPending = FALSE
        /\ hostLoopScheduled = FALSE
        /\ renderInFlight = FALSE
        /\ capturedRenderRevision = 0
        /\ acknowledgedRenderRevision = 0
        /\ producersClosed = FALSE

\* Publishing either work level also supplies the edge which starts a host
\* loop.  Repeated publications may coalesce into the standing notification.
PublishPump == /\ producersClosed = FALSE
               /\ pumpPending' = TRUE
               /\ notificationPending' = TRUE
               /\ UNCHANGED <<renderPending, renderRevision,
                               hostLoopScheduled, renderInFlight,
                               capturedRenderRevision,
                               acknowledgedRenderRevision, producersClosed>>

PublishRender == /\ producersClosed = FALSE
                 /\ renderRevision < MaxRevision
                 /\ renderPending' = TRUE
                 /\ renderRevision' = renderRevision + 1
                 /\ notificationPending' = TRUE
                 /\ UNCHANGED <<pumpPending, hostLoopScheduled,
                                 renderInFlight, capturedRenderRevision,
                                 acknowledgedRenderRevision, producersClosed>>

CloseProducers == /\ producersClosed = FALSE
                  /\ producersClosed' = TRUE
                  /\ UNCHANGED <<pumpPending, renderPending, renderRevision,
                                  notificationPending, hostLoopScheduled,
                                  renderInFlight, capturedRenderRevision,
                                  acknowledgedRenderRevision>>

\* The toolkit callback never performs work recursively.  It coalesces onto
\* one owner-thread loop and releases the callback edge.
DispatchNotification == /\ notificationPending
                        /\ notificationPending' = FALSE
                        /\ hostLoopScheduled' = TRUE
                        /\ UNCHANGED <<pumpPending, renderPending,
                                      renderRevision, renderInFlight,
                                      capturedRenderRevision,
                                      acknowledgedRenderRevision,
                                      producersClosed>>

\* A bounded controller pump discharges its current obligation.  If a render
\* is already standing, the same host loop proceeds to it without requiring a
\* second edge.
Pump == /\ hostLoopScheduled
        /\ pumpPending
        /\ pumpPending' = FALSE
        /\ hostLoopScheduled' = renderPending
        /\ UNCHANGED <<renderPending, renderRevision,
                       notificationPending, renderInFlight,
                       capturedRenderRevision,
                       acknowledgedRenderRevision, producersClosed>>

\* Pumping precedes rendering.  The captured revision is the transaction
\* token which protects a newer request published during the traversal.
BeginRender == /\ hostLoopScheduled
               /\ pumpPending = FALSE
               /\ renderPending
               /\ renderInFlight = FALSE
               /\ renderInFlight' = TRUE
               /\ capturedRenderRevision' = renderRevision
               /\ hostLoopScheduled' = FALSE
               /\ UNCHANGED <<pumpPending, renderPending, renderRevision,
                               notificationPending,
                               acknowledgedRenderRevision, producersClosed>>

CompleteRender == /\ renderInFlight
                  /\ renderInFlight' = FALSE
                  /\ acknowledgedRenderRevision' =
                        capturedRenderRevision
                  /\ renderPending' =
                        (renderRevision # capturedRenderRevision)
                  /\ hostLoopScheduled' =
                        (pumpPending \/
                         (renderRevision # capturedRenderRevision))
                  /\ UNCHANGED <<pumpPending, renderRevision,
                                  notificationPending,
                                  capturedRenderRevision, producersClosed>>

RetireEmptyLoop == /\ hostLoopScheduled
                   /\ WorkPending = FALSE
                   /\ hostLoopScheduled' = FALSE
                   /\ UNCHANGED <<pumpPending, renderPending,
                                   renderRevision, notificationPending,
                                   renderInFlight, capturedRenderRevision,
                                   acknowledgedRenderRevision,
                                   producersClosed>>

Next == \/ PublishPump
        \/ PublishRender
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
\* callback edge, an owner-thread loop, or the frame which is consuming it.
PendingWorkHasWitness ==
    WorkPending =>
        (notificationPending \/ hostLoopScheduled \/ renderInFlight)

\* A stale frame cannot acknowledge or clear a later render revision.
NoFutureAcknowledgement ==
    acknowledgedRenderRevision <= renderRevision

\* Once producers close, fair host service drains all levels and witnesses.
EventuallyQuiescent ==
    [](producersClosed =>
       <>(WorkPending = FALSE /\ notificationPending = FALSE /\
          hostLoopScheduled = FALSE /\ renderInFlight = FALSE))

=============================================================================
