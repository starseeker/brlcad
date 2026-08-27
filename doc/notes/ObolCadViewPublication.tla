---------------------- MODULE ObolCadViewPublication ----------------------
\* View-qualified preparation and completed-frame publication for one retained
\* CAD presentation node.  Geometry may be shared between presentation nodes,
\* but renderer preparation and completed work belong to the exact view which
\* traversed one node.  A changed view or demand revision cancels, rather than
\* completes, the old preparation.

EXTENDS Naturals, TLC

CONSTANT Views, MaxRevision, MaxUnits

PreparationStates == {"idle", "preparing", "complete", "constrained"}
TerminalStates == {"complete", "constrained"}

VARIABLES activeView,
          inputClosed,
          demandRevision,
          preparationState,
          preparationView,
          preparationRevision,
          remainingUnits,
          publishedView,
          publishedRevision,
          consumedView,
          consumedRevision

vars == <<activeView, inputClosed, demandRevision, preparationState,
          preparationView, preparationRevision, remainingUnits,
          publishedView, publishedRevision, consumedView, consumedRevision>>

TypeOK ==
    /\ activeView \in Views
    /\ inputClosed \in BOOLEAN
    /\ demandRevision \in [Views -> 1..MaxRevision]
    /\ preparationState \in PreparationStates
    /\ preparationView \in Views \union {"none"}
    /\ preparationRevision \in 0..MaxRevision
    /\ remainingUnits \in 0..MaxUnits
    /\ publishedView \in Views \union {"none"}
    /\ publishedRevision \in 0..MaxRevision
    /\ consumedView \in Views \union {"none"}
    /\ consumedRevision \in 0..MaxRevision

Init ==
    /\ activeView \in Views
    /\ inputClosed = FALSE
    /\ demandRevision = [view \in Views |-> 1]
    /\ preparationState = "idle"
    /\ preparationView = "none"
    /\ preparationRevision = 0
    /\ remainingUnits = 0
    /\ publishedView = "none"
    /\ publishedRevision = 0
    /\ consumedView = "none"
    /\ consumedRevision = 0

ResetPreparation ==
    /\ preparationState' = "idle"
    /\ preparationView' = "none"
    /\ preparationRevision' = 0
    /\ remainingUnits' = 0

ChangeView ==
    /\ ~inputClosed
    /\ \E view \in Views: view # activeView /\ activeView' = view
    /\ ResetPreparation
    /\ consumedView' = "none"
    /\ consumedRevision' = 0
    /\ UNCHANGED <<inputClosed, demandRevision,
                    publishedView, publishedRevision>>

SupersedeDemand ==
    /\ ~inputClosed
    /\ demandRevision[activeView] < MaxRevision
    /\ demandRevision' =
        [demandRevision EXCEPT ![activeView] = @ + 1]
    /\ ResetPreparation
    /\ consumedView' = "none"
    /\ consumedRevision' = 0
    /\ UNCHANGED <<activeView, inputClosed,
                    publishedView, publishedRevision>>

CloseInput ==
    /\ ~inputClosed
    /\ inputClosed' = TRUE
    /\ UNCHANGED <<activeView, demandRevision, preparationState,
                    preparationView, preparationRevision, remainingUnits,
                    publishedView, publishedRevision,
                    consumedView, consumedRevision>>

Admit(units) ==
    /\ preparationState = "idle"
    /\ units \in 0..MaxUnits
    /\ preparationView' = activeView
    /\ preparationRevision' = demandRevision[activeView]
    /\ remainingUnits' = units
    /\ preparationState' = IF units = 0 THEN "complete" ELSE "preparing"
    /\ UNCHANGED <<activeView, inputClosed, demandRevision,
                    publishedView, publishedRevision,
                    consumedView, consumedRevision>>

Start == \E units \in 0..MaxUnits: Admit(units)

PrepareUnit ==
    /\ preparationState = "preparing"
    /\ preparationView = activeView
    /\ preparationRevision = demandRevision[activeView]
    /\ remainingUnits > 0
    /\ remainingUnits' = remainingUnits - 1
    /\ preparationState' =
        IF remainingUnits = 1 THEN "complete" ELSE "preparing"
    /\ UNCHANGED <<activeView, inputClosed, demandRevision,
                    preparationView, preparationRevision,
                    publishedView, publishedRevision,
                    consumedView, consumedRevision>>

Constrain ==
    /\ preparationState = "preparing"
    /\ preparationView = activeView
    /\ preparationRevision = demandRevision[activeView]
    /\ preparationState' = "constrained"
    /\ UNCHANGED <<activeView, inputClosed, demandRevision,
                    preparationView, preparationRevision, remainingUnits,
                    publishedView, publishedRevision,
                    consumedView, consumedRevision>>

Publish ==
    /\ preparationState = "complete"
    /\ preparationView = activeView
    /\ preparationRevision = demandRevision[activeView]
    /\ remainingUnits = 0
    /\ publishedView' = preparationView
    /\ publishedRevision' = preparationRevision
    /\ UNCHANGED <<activeView, inputClosed, demandRevision,
                    preparationState, preparationView,
                    preparationRevision, remainingUnits,
                    consumedView, consumedRevision>>

\* A completed report remains an immutable historical sample when the view
\* changes.  Consumers may accept it only when its complete stamp matches the
\* current target; otherwise the stale report is harmless and detectable.
AcceptPublished ==
    /\ publishedView = activeView
    /\ publishedRevision = demandRevision[activeView]
    /\ consumedView' = publishedView
    /\ consumedRevision' = publishedRevision
    /\ UNCHANGED <<activeView, inputClosed, demandRevision,
                    preparationState, preparationView,
                    preparationRevision, remainingUnits,
                    publishedView, publishedRevision>>

Next ==
    \/ ChangeView
    \/ SupersedeDemand
    \/ CloseInput
    \/ Start
    \/ PrepareUnit
    \/ Constrain
    \/ Publish
    \/ AcceptPublished

Spec == Init /\ [][Next]_vars
        /\ WF_vars(Start)
        /\ WF_vars(PrepareUnit)
        /\ WF_vars(Constrain)
        /\ WF_vars(Publish)

PreparationIsViewQualified ==
    preparationState # "idle" =>
        /\ preparationView = activeView
        /\ preparationRevision = demandRevision[activeView]

CompleteHasNoRemainingWork ==
    preparationState = "complete" => remainingUnits = 0

PublishedWorkIsStamped ==
    publishedView # "none" =>
        publishedRevision <= demandRevision[publishedView]

ConsumedWorkMatchesCurrentView ==
    consumedView # "none" =>
        /\ consumedView = activeView
        /\ consumedRevision = demandRevision[activeView]

EventuallyTerminalAfterStableInput ==
    [](inputClosed => <>(preparationState \in TerminalStates))

=============================================================================
