------------------- MODULE ObolCadFrameComposition -------------------
\* Composition seam from atomic retained-scene mutation through exact-target
\* renderer preparation, host frame ownership, and accepted publication.
\*
\* Scene mutation and rendering execute on one owner thread in production,
\* but a mutation may supersede an already prepared or requested target
\* between frames.  A completed historical preparation/report is harmless;
\* only a report whose complete scene/view stamp is still current may replace
\* the committed framebuffer.

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
          reportSceneRevision,
          reportViewRevision,
          committedSceneRevision,
          committedViewRevision

vars == <<inputOpen, sceneRevision, viewRevision, mutationStaged,
          updateOpen, mutationRemaining, preparationState,
          targetSceneRevision, targetViewRevision, preparationRemaining,
          renderPending, renderInFlight, reportSceneRevision,
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
    /\ UNCHANGED <<inputOpen, viewRevision, renderInFlight,
                    reportSceneRevision, reportViewRevision,
                    committedSceneRevision, committedViewRevision>>

ChangeView ==
    /\ inputOpen
    /\ viewRevision < MaxViewRevision
    /\ viewRevision' = viewRevision + 1
    /\ ResetPreparation
    /\ renderPending' = TRUE
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
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, mutationStaged,
                    updateOpen, mutationRemaining, renderPending,
                    renderInFlight, reportSceneRevision, reportViewRevision,
                    committedSceneRevision, committedViewRevision>>

Start == \E units \in 0..MaxUnits: StartPreparation(units)

PrepareUnit ==
    /\ preparationState = "preparing"
    /\ targetSceneRevision = sceneRevision
    /\ targetViewRevision = viewRevision
    /\ preparationRemaining > 0
    /\ preparationRemaining' = preparationRemaining - 1
    /\ preparationState' =
        IF preparationRemaining = 1 THEN "complete" ELSE "preparing"
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, mutationStaged,
                    updateOpen, mutationRemaining, targetSceneRevision,
                    targetViewRevision, renderPending, renderInFlight,
                    reportSceneRevision, reportViewRevision,
                    committedSceneRevision, committedViewRevision>>

ConstrainPreparation ==
    /\ preparationState = "preparing"
    /\ targetSceneRevision = sceneRevision
    /\ targetViewRevision = viewRevision
    /\ preparationState' = "constrained"
    /\ renderPending' = TRUE
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, mutationStaged,
                    updateOpen, mutationRemaining, targetSceneRevision,
                    targetViewRevision, preparationRemaining, renderInFlight,
                    reportSceneRevision, reportViewRevision,
                    committedSceneRevision, committedViewRevision>>

RequestFrame ==
    /\ preparationState \in {"complete", "constrained"}
    /\ targetSceneRevision = sceneRevision
    /\ targetViewRevision = viewRevision
    /\ ~renderInFlight
    /\ renderPending' = TRUE
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, mutationStaged,
                    updateOpen, mutationRemaining, preparationState,
                    targetSceneRevision, targetViewRevision,
                    preparationRemaining, renderInFlight,
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
                    reportViewRevision, committedSceneRevision,
                    committedViewRevision>>

\* A traversal records the target it actually used.  Acceptance is a separate
\* exact-stamp edge, so a scene/view change during an in-flight host frame can
\* never acknowledge future or stale work.
CompleteFrame ==
    /\ renderInFlight
    /\ reportSceneRevision' = targetSceneRevision
    /\ reportViewRevision' = targetViewRevision
    /\ renderInFlight' = FALSE
    /\ UNCHANGED <<inputOpen, sceneRevision, viewRevision, mutationStaged,
                    updateOpen, mutationRemaining, preparationState,
                    targetSceneRevision, targetViewRevision,
                    preparationRemaining, renderPending,
                    committedSceneRevision, committedViewRevision>>

AcceptReport ==
    /\ reportSceneRevision = sceneRevision
    /\ reportViewRevision = viewRevision
    /\ preparationState \in {"complete", "constrained"}
    /\ targetSceneRevision = reportSceneRevision
    /\ targetViewRevision = reportViewRevision
    /\ committedSceneRevision' = reportSceneRevision
    /\ committedViewRevision' = reportViewRevision
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
    \/ ChangeView
    \/ Start
    \/ PrepareUnit
    \/ ConstrainPreparation
    \/ RequestFrame
    \/ BeginFrame
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
    /\ WF_vars(PrepareUnit)
    /\ WF_vars(ConstrainPreparation)
    /\ WF_vars(RequestFrame)
    /\ WF_vars(BeginFrame)
    /\ WF_vars(CompleteFrame)
    /\ WF_vars(AcceptReport)
    /\ WF_vars(CloseInput)

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

Terminal ==
    /\ ~inputOpen
    /\ ~mutationStaged
    /\ ~updateOpen
    /\ ~renderPending
    /\ ~renderInFlight
    /\ committedSceneRevision = sceneRevision
    /\ committedViewRevision = viewRevision

EventuallyTerminal == <>Terminal

=============================================================================
