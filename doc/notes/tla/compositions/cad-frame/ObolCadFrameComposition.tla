------------------- MODULE ObolCadFrameComposition -------------------
\* Composition seam from atomic retained-scene mutation through exact-target
\* renderer preparation, host frame ownership, and accepted publication.
\*
\* Scene mutation and rendering execute on one owner thread in production,
\* but a mutation may supersede an already prepared or requested target
\* between frames.  A completed historical preparation/report is harmless;
\* only a report whose complete scene/view-control stamp is still current may
\* replace the committed framebuffer.  viewRevision abstracts both camera
\* changes and renderer-local presentation controls (progressive ceiling,
\* fractional cut, point threshold, and motion reuse); production may prove
\* the latter by direct value equality rather than a separate numeric stamp.

EXTENDS Naturals, TLC

CONSTANT MaxSceneRevision, MaxViewRevision, MaxUnits, MaxBatches

PreparationStates == {"idle", "preparing", "complete", "constrained"}

VARIABLES inputOpen,
          sceneRevision,
          viewRevision,
          mutationStaged,
          updateOpen,
          mutationRemaining,
          preparationState,
          targetSceneRevision,
          targetViewRevision,
          preparationRemaining,
          renderPending,
          renderInFlight,
          exactFrameRequired,
          reportSceneRevision,
          reportViewRevision,
          committedSceneRevision,
          committedViewRevision

vars == <<inputOpen, sceneRevision, viewRevision, mutationStaged,
          updateOpen, mutationRemaining, preparationState,
          targetSceneRevision, targetViewRevision, preparationRemaining,
          renderPending, renderInFlight, exactFrameRequired,
          reportSceneRevision,
          reportViewRevision, committedSceneRevision,
          committedViewRevision>>

TypeOK ==
    /\ inputOpen \in BOOLEAN
    /\ sceneRevision \in 1..MaxSceneRevision
    /\ viewRevision \in 1..MaxViewRevision
    /\ mutationStaged \in BOOLEAN
    /\ updateOpen \in BOOLEAN
    /\ mutationRemaining \in 0..MaxBatches
    /\ preparationState \in PreparationStates
    /\ targetSceneRevision \in 0..MaxSceneRevision
    /\ targetViewRevision \in 0..MaxViewRevision
    /\ preparationRemaining \in 0..MaxUnits
    /\ renderPending \in BOOLEAN
    /\ renderInFlight \in BOOLEAN
    /\ exactFrameRequired \in BOOLEAN
    /\ reportSceneRevision \in 0..MaxSceneRevision
    /\ reportViewRevision \in 0..MaxViewRevision
    /\ committedSceneRevision \in 0..MaxSceneRevision
    /\ committedViewRevision \in 0..MaxViewRevision

Init ==
    /\ inputOpen = TRUE
    /\ sceneRevision = 1
    /\ viewRevision = 1
    /\ mutationStaged = FALSE
    /\ updateOpen = FALSE
    /\ mutationRemaining = 0
    /\ preparationState = "idle"
    /\ targetSceneRevision = 0
    /\ targetViewRevision = 0
    /\ preparationRemaining = 0
    /\ renderPending = FALSE
    /\ renderInFlight = FALSE
    /\ exactFrameRequired = FALSE
    /\ reportSceneRevision = 0
    /\ reportViewRevision = 0
    /\ committedSceneRevision = 0
    /\ committedViewRevision = 0

ResetPreparation ==
    /\ preparationState' = "idle"
    /\ targetSceneRevision' = 0
    /\ targetViewRevision' = 0
    /\ preparationRemaining' = 0

StageMutation ==
    /\ inputOpen
    /\ sceneRevision < MaxSceneRevision
    /\ ~mutationStaged
    /\ ~updateOpen
    /\ mutationStaged' = TRUE
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, updateOpen,
                    mutationRemaining, preparationState,
                    targetSceneRevision, targetViewRevision,
                    preparationRemaining, renderPending, renderInFlight,
                    exactFrameRequired,
                    reportSceneRevision, reportViewRevision,
                    committedSceneRevision, committedViewRevision>>

BeginMutation ==
    /\ mutationStaged
    /\ ~updateOpen
    /\ updateOpen' = TRUE
    /\ mutationRemaining' \in 1..MaxBatches
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, mutationStaged,
                    preparationState, targetSceneRevision,
                    targetViewRevision, preparationRemaining,
                    renderPending, renderInFlight, reportSceneRevision,
                    exactFrameRequired,
                    reportViewRevision, committedSceneRevision,
                    committedViewRevision>>

ApplyMutationBatch ==
    /\ updateOpen
    /\ mutationRemaining > 0
    /\ mutationRemaining' = mutationRemaining - 1
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, mutationStaged,
                    updateOpen, preparationState, targetSceneRevision,
                    targetViewRevision, preparationRemaining,
                    renderPending, renderInFlight, reportSceneRevision,
                    exactFrameRequired,
                    reportViewRevision, committedSceneRevision,
                    committedViewRevision>>

CommitMutation ==
    /\ updateOpen
    /\ mutationRemaining = 0
    /\ sceneRevision < MaxSceneRevision
    /\ sceneRevision' = sceneRevision + 1
    /\ mutationStaged' = FALSE
    /\ updateOpen' = FALSE
    /\ mutationRemaining' = 0
    /\ ResetPreparation
    /\ renderPending' = TRUE
    /\ exactFrameRequired' = TRUE
    /\ UNCHANGED <<inputOpen, viewRevision, renderInFlight,
                    reportSceneRevision, reportViewRevision,
                    committedSceneRevision, committedViewRevision>>

ChangeViewOrPresentationControls ==
    /\ inputOpen
    /\ viewRevision < MaxViewRevision
    /\ viewRevision' = viewRevision + 1
    /\ ResetPreparation
    /\ renderPending' = TRUE
    /\ exactFrameRequired' = TRUE
    /\ UNCHANGED <<inputOpen, sceneRevision, mutationStaged, updateOpen,
                    mutationRemaining, renderInFlight, reportSceneRevision,
                    reportViewRevision, committedSceneRevision,
                    committedViewRevision>>

StartPreparation(units) ==
    /\ ~updateOpen
    /\ preparationState = "idle"
    /\ units \in 0..MaxUnits
    /\ targetSceneRevision' = sceneRevision
    /\ targetViewRevision' = viewRevision
    /\ preparationRemaining' = units
    /\ preparationState' = IF units = 0 THEN "complete" ELSE "preparing"
    /\ renderPending' = TRUE
    /\ exactFrameRequired' = TRUE
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, mutationStaged,
                    updateOpen, mutationRemaining, renderInFlight,
                    reportSceneRevision, reportViewRevision,
                    committedSceneRevision, committedViewRevision>>

Start == \E units \in 0..MaxUnits: StartPreparation(units)


\* SoCADAssembly command preparation is render-sliced: an inexact traversal
\* may advance one bounded unit, but it cannot publish a report or consume the
\* exact-presentation obligation.  The next frame edge is level-triggered.
IncompleteFrame ==
    /\ renderInFlight
    /\ preparationState = "preparing"
    /\ targetSceneRevision = sceneRevision
    /\ targetViewRevision = viewRevision
    /\ preparationRemaining > 0
    /\ preparationRemaining' = preparationRemaining - 1
    /\ preparationState' =
        IF preparationRemaining = 1 THEN "complete" ELSE "preparing"
    /\ renderInFlight' = FALSE
    /\ renderPending' = TRUE
    /\ exactFrameRequired' = TRUE
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, mutationStaged,
                    updateOpen, mutationRemaining, targetSceneRevision,
                    targetViewRevision,
                    reportSceneRevision, reportViewRevision,
                    committedSceneRevision, committedViewRevision>>

\* A scene/view supersession may invalidate the preparation target while an
\* older host frame is in flight.  Its completion is not a report; it merely
\* transfers ownership back to the current level-triggered frame request.
AbandonStaleFrame ==
    /\ renderInFlight
    /\ \/ preparationState = "idle"
       \/ targetSceneRevision # sceneRevision
       \/ targetViewRevision # viewRevision
    /\ renderInFlight' = FALSE
    /\ renderPending' = TRUE
    /\ exactFrameRequired' = TRUE
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, mutationStaged,
                    updateOpen, mutationRemaining, preparationState,
                    targetSceneRevision, targetViewRevision,
                    preparationRemaining, reportSceneRevision,
                    reportViewRevision, committedSceneRevision,
                    committedViewRevision>>

ConstrainPreparation ==
    /\ preparationState = "preparing"
    /\ ~renderInFlight
    /\ targetSceneRevision = sceneRevision
    /\ targetViewRevision = viewRevision
    /\ preparationState' = "constrained"
    /\ renderPending' = TRUE
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, mutationStaged,
                    updateOpen, mutationRemaining, targetSceneRevision,
                    targetViewRevision, preparationRemaining, renderInFlight,
                    exactFrameRequired,
                    reportSceneRevision, reportViewRevision,
                    committedSceneRevision, committedViewRevision>>

RequestFrame ==
    /\ exactFrameRequired
    /\ preparationState \in {"complete", "constrained"}
    /\ targetSceneRevision = sceneRevision
    /\ targetViewRevision = viewRevision
    /\ ~renderInFlight
    /\ ~renderPending
    /\ renderPending' = TRUE
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, mutationStaged,
                    updateOpen, mutationRemaining, preparationState,
                    targetSceneRevision, targetViewRevision,
                    preparationRemaining, renderInFlight,
                    exactFrameRequired,
                    reportSceneRevision, reportViewRevision,
                    committedSceneRevision, committedViewRevision>>

BeginFrame ==
    /\ renderPending
    /\ ~renderInFlight
    /\ renderPending' = FALSE
    /\ renderInFlight' = TRUE
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, mutationStaged,
                    updateOpen, mutationRemaining, preparationState,
                    targetSceneRevision, targetViewRevision,
                    preparationRemaining, reportSceneRevision,
                    exactFrameRequired,
                    reportViewRevision, committedSceneRevision,
                    committedViewRevision>>

\* A traversal records the target it actually used.  Acceptance is a separate
\* exact-stamp edge, so a scene/view change during an in-flight host frame can
\* never acknowledge future or stale work.
CompleteFrame ==
    /\ renderInFlight
    /\ preparationState \in {"complete", "constrained"}
    /\ targetSceneRevision = sceneRevision
    /\ targetViewRevision = viewRevision
    /\ reportSceneRevision' = targetSceneRevision
    /\ reportViewRevision' = targetViewRevision
    /\ renderInFlight' = FALSE
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, mutationStaged,
                    updateOpen, mutationRemaining, preparationState,
                    targetSceneRevision, targetViewRevision,
                    preparationRemaining, renderPending,
                    exactFrameRequired,
                    committedSceneRevision, committedViewRevision>>

AcceptReport ==
    /\ reportSceneRevision = sceneRevision
    /\ reportViewRevision = viewRevision
    /\ preparationState \in {"complete", "constrained"}
    /\ targetSceneRevision = reportSceneRevision
    /\ targetViewRevision = reportViewRevision
    /\ committedSceneRevision' = reportSceneRevision
    /\ committedViewRevision' = reportViewRevision
    /\ exactFrameRequired' = FALSE
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, mutationStaged,
                    updateOpen, mutationRemaining, preparationState,
                    targetSceneRevision, targetViewRevision,
                    preparationRemaining, renderPending, renderInFlight,
                    reportSceneRevision, reportViewRevision>>

CloseInput ==
    /\ inputOpen
    /\ inputOpen' = FALSE
    /\ UNCHANGED <<sceneRevision, viewRevision, mutationStaged, updateOpen,
                    mutationRemaining, preparationState,
                    targetSceneRevision, targetViewRevision,
                    preparationRemaining, renderPending, renderInFlight,
                    exactFrameRequired,
                    reportSceneRevision, reportViewRevision,
                    committedSceneRevision, committedViewRevision>>

Quiescent ==
    /\ ~inputOpen
    /\ ~mutationStaged
    /\ ~updateOpen
    /\ ~renderPending
    /\ ~renderInFlight
    /\ committedSceneRevision = sceneRevision
    /\ committedViewRevision = viewRevision
    /\ UNCHANGED vars

Next ==
    \/ StageMutation
    \/ BeginMutation
    \/ ApplyMutationBatch
    \/ CommitMutation
    \/ ChangeViewOrPresentationControls
    \/ Start
    \/ ConstrainPreparation
    \/ RequestFrame
    \/ BeginFrame
    \/ IncompleteFrame
    \/ AbandonStaleFrame
    \/ CompleteFrame
    \/ AcceptReport
    \/ CloseInput
    \/ Quiescent

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(BeginMutation)
    /\ WF_vars(ApplyMutationBatch)
    /\ WF_vars(CommitMutation)
    /\ WF_vars(Start)
    /\ WF_vars(ConstrainPreparation)
    /\ WF_vars(RequestFrame)
    /\ WF_vars(BeginFrame)
    /\ WF_vars(IncompleteFrame)
    /\ WF_vars(AbandonStaleFrame)
    /\ WF_vars(CompleteFrame)
    /\ WF_vars(AcceptReport)

OpenMutationOwnsFiniteWork ==
    updateOpen => mutationStaged /\ mutationRemaining \in 0..MaxBatches

PreparationIsExactTarget ==
    preparationState # "idle" =>
        /\ targetSceneRevision = sceneRevision
        /\ targetViewRevision = viewRevision

CompletePreparationHasNoWork ==
    preparationState = "complete" => preparationRemaining = 0

CommittedFrameIsCurrent ==
    committedSceneRevision # 0 =>
        /\ committedSceneRevision <= sceneRevision
        /\ committedViewRevision <= viewRevision

NoUnpreparedCurrentCommit ==
    committedSceneRevision = sceneRevision /\
    committedViewRevision = viewRevision =>
        /\ preparationState \in {"complete", "constrained"}
        /\ targetSceneRevision = sceneRevision
        /\ targetViewRevision = viewRevision

ExactPresentationHasSuccessor ==
    exactFrameRequired =>
        \/ renderPending
        \/ renderInFlight
        \/ /\ reportSceneRevision = sceneRevision
           /\ reportViewRevision = viewRevision
           /\ preparationState \in {"complete", "constrained"}
           /\ targetSceneRevision = reportSceneRevision
           /\ targetViewRevision = reportViewRevision

Terminal ==
    /\ ~inputOpen
    /\ ~mutationStaged
    /\ ~updateOpen
    /\ ~renderPending
    /\ ~renderInFlight
    /\ ~exactFrameRequired
    /\ committedSceneRevision = sceneRevision
    /\ committedViewRevision = viewRevision

Ready == Terminal /\ preparationState = "complete"
Constrained == Terminal /\ preparationState = "constrained"
Failed == FALSE

EventuallyTerminalAfterInputCloses ==
    [](~inputOpen => <>Terminal)

=============================================================================
